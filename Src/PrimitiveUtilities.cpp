#include "Node.h"
#include "Solver.h"

#include <tetgen.h>

#include <cstdlib>

namespace Pies {
namespace {
float randf() { return static_cast<float>(double(std::rand()) / RAND_MAX); }

glm::vec3 randColor() { return glm::vec3(randf(), randf(), randf()); }

struct GridId {
  uint32_t x;
  uint32_t y;
  uint32_t z;
};

struct Grid {
  uint32_t width;
  uint32_t height;
  uint32_t depth;

  GridId nodeIdToGridId(uint32_t nodeId) {
    GridId gridId;

    gridId.z = nodeId % width;
    gridId.y = (nodeId / depth) % height;
    gridId.x = (nodeId / depth) / height;

    return gridId;
  }

  uint32_t gridIdToNodeId(size_t nodeIdOffset, const GridId& gridId) {
    return gridId.z + depth * (gridId.y + height * gridId.x) +
           static_cast<uint32_t>(nodeIdOffset);
  }
};
} // namespace

void Solver::addNodes(const std::vector<glm::vec3>& vertices) {

  // TODO: parameterize more of these
  float mass = 1.0f;
  float radius = 0.5f;
  glm::vec3 initialVelocity(0.0f);

  glm::vec3 color = randColor();
  float roughness = randf();
  float metallic = static_cast<float>(std::rand() % 2);

  size_t currentNodeCount = this->_nodes.size();
  this->_nodes.reserve(currentNodeCount + vertices.size());
  for (uint32_t i = 0; i < static_cast<uint32_t>(vertices.size()); ++i) {
    Node& node = this->_nodes.emplace_back();
    node.id = currentNodeCount + i;
    node.position = vertices[i];
    node.prevPosition = node.position;
    node.velocity = initialVelocity;
    node.radius = radius;
    node.invMass = 1.0f / mass;
  }

  this->_vertices.resize(this->_nodes.size());
  for (size_t i = currentNodeCount; i < this->_nodes.size(); ++i) {
    this->_vertices[i].position = this->_nodes[i].position;
    this->_vertices[i].radius = this->_nodes[i].radius;
    this->_vertices[i].baseColor = color;
    this->_vertices[i].roughness = roughness;
    this->_vertices[i].metallic = metallic;
  }

  this->renderStateDirty = true;
}

void Solver::addFixedRegions(
    const std::vector<glm::mat4>& regionMatrices,
    float w) {
  // Adds fixed position constraints to all nodes inside any of the bounding
  // boxes specified by the region matrices

  // This is not particularly efficient, but it should only need to be done once
  // during setup
  // TODO: Could speed up with spatial hash if this is super slow
  size_t currentGoalConstraintCount = this->_goalConstraints.size();
  this->_goalConstraints.reserve(
      currentGoalConstraintCount + regionMatrices.size());

  size_t currentFixedRegionsCount = this->_fixedRegions.size();
  this->_fixedRegions.reserve(currentFixedRegionsCount + regionMatrices.size());
  for (const glm::mat4& regionToWorld : regionMatrices) {
    FixedRegion region;
    region.initialTransform = regionToWorld;
    region.invInitialTransform = glm::inverse(regionToWorld);
    region.goalMatchingConstraint = this->_goalConstraints.size();

    std::vector<uint32_t> constrainedNodes;

    for (const Node& node : this->_nodes) {
      glm::vec4 pos = glm::vec4(node.position, 1.0f);
      glm::vec4 local = region.invInitialTransform * pos;
      if (-1.0f <= local.x && local.x <= 1.0f && -1.0f <= local.y &&
          local.y <= 1.0f && -1.0f <= local.z && local.z <= 1.0f) {
        constrainedNodes.push_back(node.id);
      }
    }

    this->_goalConstraints.emplace_back(this->_nodes, constrainedNodes, w);
    this->_fixedRegions.push_back(region);
  }
}

void Solver::updateFixedRegions(const std::vector<glm::mat4>& regionMatrices) {
  if (regionMatrices.size() != this->_fixedRegions.size()) {
    assert(false);
    return;
  }

  for (uint32_t i = 0; i < this->_fixedRegions.size(); ++i) {
    const glm::mat4& currentRegionTransform = regionMatrices[i];
    const FixedRegion& region = this->_fixedRegions[i];
    // transform = regionToCurrentWorld * initialWorldToRegion
    glm::mat4 transform = currentRegionTransform * region.invInitialTransform;
    this->_goalConstraints[region.goalMatchingConstraint].setTransform(
        transform);
  }
}

void Solver::addLinkedRegions(
    const std::vector<glm::mat4>& regionMatrices,
    float w) {
  // Adds shape matching constraints to all nodes inside each of the
  // bounding boxes specified by the region matrices

  // This is not particularly efficient, but it should only need to be done once
  // during setup
  // TODO: Could speed up with spatial hash if this is super slow
  std::vector<glm::vec3> materialCoords;
  std::vector<uint32_t> nodeIndices;
  for (const glm::mat4& region : regionMatrices) {
    glm::mat4 worldToRegion = glm::inverse(region);

    for (const Node& node : this->_nodes) {
      glm::vec4 pos = glm::vec4(node.position, 1.0f);
      glm::vec4 local = worldToRegion * pos;
      if (-1.0f <= local.x && local.x <= 1.0f && -1.0f <= local.y &&
          local.y <= 1.0f && -1.0f <= local.z && local.z <= 1.0f) {
        materialCoords.push_back(node.position);
        nodeIndices.push_back(node.id);
      }
    }

    if (materialCoords.size() >= 3) {
      this->_shapeConstraints
          .emplace_back(this->_nodes, nodeIndices, materialCoords, w);
    }

    materialCoords.clear();
    nodeIndices.clear();
  }
}

void Solver::addTriMeshVolume(
    const std::vector<glm::vec3>& vertices,
    const std::vector<uint32_t>& indices,
    const glm::vec3& initialVelocity,
    float density,
    float strainStiffness,
    float minStrain,
    float maxStrain,
    float volumeStiffness,
    float compression,
    float stretching) {

  float mass = density;
  float radius = 0.5f;

  glm::vec3 color = randColor();
  float roughness = randf();
  float metallic = static_cast<float>(std::rand() % 2);

  tetgenio tetgenInput{};
  tetgenInput.numberofpoints = static_cast<int>(vertices.size());
  tetgenInput.pointlist = new double[vertices.size() * 3];
  for (size_t i = 0; i < vertices.size(); ++i) {
    const glm::vec3& vertex = vertices[i];
    tetgenInput.pointlist[3 * i + 0] = vertex.x;
    tetgenInput.pointlist[3 * i + 1] = vertex.y;
    tetgenInput.pointlist[3 * i + 2] = vertex.z;
  }

  tetgenInput.numberoffacets = static_cast<int>(indices.size() / 3);
  tetgenInput.facetlist = new tetgenio::facet[tetgenInput.numberoffacets];
  for (int i = 0; i < tetgenInput.numberoffacets; ++i) {
    tetgenio::facet& face = tetgenInput.facetlist[i];
    face.numberofpolygons = 1;
    face.polygonlist = new tetgenio::polygon[1];

    face.polygonlist[0].numberofvertices = 3;
    face.polygonlist[0].vertexlist = new int[3];
    face.polygonlist[0].vertexlist[0] = indices[3 * i + 0];
    face.polygonlist[0].vertexlist[1] = indices[3 * i + 1];
    face.polygonlist[0].vertexlist[2] = indices[3 * i + 2];

    face.numberofholes = 0;
    face.holelist = nullptr;
  }

  tetgenio tetgenInterm{};
  tetgenio tetgenOutput{};
  tetgenbehavior behavior{};
  behavior.plc = 1;
  behavior.facesout = 1;
  behavior.neighout = 2;
  behavior.zeroindex = 1;
  behavior.quality = 1;
  behavior.minratio = 1.5;

  behavior.regionattrib = 1;
  // 1414;

  // tetgenbehavior behavior{};
  // // behavior.addinfilename = "pq1.414a0.1aA";
  // behavior.cdt = 1;
  // behavior.plc = 1;
  // behavior.refine = 1;
  // behavior.quality = 1;
  // behavior.minratio = 1.1;
  // //414;

  // behavior.varvolume = 1;
  // behavior.maxvolume = 0.1;
  // //behavior.nomergefacet = 1;
  // behavior.facesout = 1;
  // behavior.zeroindex = 1;
  // //behavior.regionattrib = 1;
  // //behavior.plc = 1;
  // //behavior.nofacewritten = 1;

  tetrahedralize(&behavior, &tetgenInput, &tetgenOutput);

  size_t existingNodesCount = this->_nodes.size();
  size_t existingTetsCount = this->_tetConstraints.size();
  size_t existingTrisCount = this->_triangles.size();

  this->_triangles.reserve(existingTrisCount + tetgenOutput.numberoftrifaces);
  for (size_t i = 0; i < tetgenOutput.numberoftrifaces; ++i) {
    int v0 = tetgenOutput.trifacelist[3 * i];
    int v1 = tetgenOutput.trifacelist[3 * i + 1];
    int v2 = tetgenOutput.trifacelist[3 * i + 2];

    int t0 = tetgenOutput.face2tetlist[2 * i];
    int t1 = tetgenOutput.face2tetlist[2 * i + 1];

    if (t0 >= 0 && t1 >= 0) {
      continue;
    }

    Triangle& tri = this->_triangles.emplace_back();

    // Switch winding so normals point outward
    tri.nodeIds[0] = static_cast<uint32_t>(existingNodesCount + v0);
    tri.nodeIds[1] = static_cast<uint32_t>(existingNodesCount + v2);
    tri.nodeIds[2] = static_cast<uint32_t>(existingNodesCount + v1);
  }

  this->_nodes.reserve(existingNodesCount + tetgenOutput.numberofpoints);
  for (int i = 0; i < tetgenOutput.numberofpoints; i++) {
    Node& node = this->_nodes.emplace_back();
    node.id = static_cast<uint32_t>(this->_nodes.size() - 1);
    node.position = glm::vec3(
        static_cast<float>(tetgenOutput.pointlist[3 * i + 0]),
        static_cast<float>(tetgenOutput.pointlist[3 * i + 1]),
        static_cast<float>(tetgenOutput.pointlist[3 * i + 2]));
    node.prevPosition = node.position;
    node.velocity = initialVelocity;
    node.radius = radius;
    node.invMass = 1.0f / mass;
  }

  this->_tetConstraints.reserve(
      existingTetsCount + tetgenOutput.numberoftetrahedra);
  for (int i = 0; i < tetgenOutput.numberoftetrahedra; ++i) {
    int v1 = tetgenOutput.tetrahedronlist[4 * i + 0];
    int v2 = tetgenOutput.tetrahedronlist[4 * i + 1];
    int v3 = tetgenOutput.tetrahedronlist[4 * i + 2];
    int v4 = tetgenOutput.tetrahedronlist[4 * i + 3];

    if (strainStiffness != 0.0f) {
      this->_tetConstraints.emplace_back(
          this->_nodes[existingNodesCount + v1],
          this->_nodes[existingNodesCount + v2],
          this->_nodes[existingNodesCount + v3],
          this->_nodes[existingNodesCount + v4],
          strainStiffness,
          minStrain,
          maxStrain,
          compression,
          stretching);
    }
  }

  this->_vertices.resize(this->_nodes.size());
  for (size_t i = existingNodesCount; i < this->_nodes.size(); ++i) {
    this->_vertices[i].position = this->_nodes[i].position;
    this->_vertices[i].radius = this->_nodes[i].radius;
    this->_vertices[i].baseColor = color;
    this->_vertices[i].roughness = roughness;
    this->_vertices[i].metallic = metallic;
  }

  this->renderStateDirty = true;
}
} // Pies
