#include "Node.h"
#include "Solver.h"

#include <tetgen.h>

#include <cstdlib>
#include <unordered_map>

namespace {
float randf() { return static_cast<float>(double(std::rand()) / RAND_MAX); }

glm::vec3 randColor() { return glm::vec3(randf(), randf(), randf()); }

struct Edge {
  uint32_t vertexId1;
  uint32_t vertexId2; // where vertexId1 < vertexId2

  bool operator==(const Edge& other) const {
    return vertexId1 == other.vertexId1 && vertexId2 == other.vertexId2;
  }
};

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

namespace std {
  template <> struct hash<Edge> { 
    std::size_t operator()(const Edge& k) const {

      std::hash<uint32_t> h;
      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:

      return (h(k.vertexId1) << 1) ^ h(k.vertexId2);     
    }
  };
}


namespace Pies {

void Solver::addClothMesh(
    const std::vector<glm::vec3>& vertices,
    const std::vector<uint32_t>& indices,
    float w) {

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

  struct Adjacency {
    uint32_t triId; // where 3 * triId is the offset to triangle indicies
    std::optional<uint32_t> triId2;
  };

  std::unordered_map<Edge, Adjacency> adjacencyMap;

  for (uint32_t i = 0; i < indices.size() - 1; ++i) {
    Edge e{};
    e.vertexId1 = indices[i];
    e.vertexId2 = indices[i + 1];

    if (e.vertexId1 > e.vertexId2) {
      std::swap(e.vertexId1, e.vertexId2);
    }
    
    auto eIt = adjacencyMap.find(e);
    if (eIt == adjacencyMap.end()) {
      adjacencyMap.emplace(e, Adjacency{i/3, std::nullopt});
    } else {
      eIt->second.triId2 = i/3;
    }
  }

  for (uint32_t i = 0; i < indices.size(); i += 3) {
    Triangle t;
    t.nodeIds[0] = indices[i] + currentNodeCount;
    t.nodeIds[1] = indices[i + 1] + currentNodeCount;
    t.nodeIds[2] = indices[i + 2] + currentNodeCount;

    this->_triangles.push_back(t);
  }
    
  //for each over adjacencyMap , for auto adjIt in adjMap

  for (auto adjIt : adjacencyMap) {
    uint32_t vId1 = adjIt.first.vertexId1;
    uint32_t vId2 = adjIt.first.vertexId2;
    
    this->_distanceConstraints.push_back(createDistanceConstraint(
        this->_constraintId++,
        this->_nodes[currentNodeCount + vId1],
        this->_nodes[currentNodeCount + vId2],
		w));
    
    //TODO: bend constraint
    if (adjIt.second.triId2.has_value()) {
      
      
      uint32_t triId1 = adjIt.second.triId;
      uint32_t triId2 = adjIt.second.triId2.value();
      
      if (triId1 == vId1 || triId1 == vId2) {
        continue;
      }

      this->_bendConstraints.push_back(createBendConstraint(
          this->_constraintId++,
          w,
          this->_nodes[triId1],
          this->_nodes[vId1],
          this->_nodes[vId2],
          this->_nodes[triId2]));
    }

  }
    
  //add triangles
 
  //list of pairs, <nodeid, nodeid>

  this->renderStateDirty = true;

  return;
}

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
  std::vector<glm::mat4> worldToRegionMatrices;
  worldToRegionMatrices.reserve(regionMatrices.size());
  for (const glm::mat4& region : regionMatrices) {
    worldToRegionMatrices.push_back(glm::inverse(region));
  }

  for (const Node& node : this->_nodes) {
    glm::vec4 pos = glm::vec4(node.position, 1.0f);

    for (const glm::mat4& worldToRegion : worldToRegionMatrices) {
      glm::vec4 local = worldToRegion * pos;
      if (-1.0f <= local.x && local.x <= 1.0f && -1.0f <= local.y &&
          local.y <= 1.0f && -1.0f <= local.z && local.z <= 1.0f) {
        this->_positionConstraints.push_back(
            createPositionConstraint(this->_constraintId++, node, w));
        break;
      }
    }
  }
}

void Solver::addTriMeshVolume(
    const std::vector<glm::vec3>& vertices,
    const std::vector<uint32_t>& indices,
    float w) {

  // TODO: parameterize more of these
  float mass = 1.0f;
  float radius = 0.5f;
  glm::vec3 initialVelocity(0.0f);

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
  size_t existingVolumesCount = this->_volumeConstraints.size();
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
  this->_volumeConstraints.reserve(
      existingVolumesCount + tetgenOutput.numberoftetrahedra);
  for (int i = 0; i < tetgenOutput.numberoftetrahedra; ++i) {
    int v1 = tetgenOutput.tetrahedronlist[4 * i + 0];
    int v2 = tetgenOutput.tetrahedronlist[4 * i + 1];
    int v3 = tetgenOutput.tetrahedronlist[4 * i + 2];
    int v4 = tetgenOutput.tetrahedronlist[4 * i + 3];

    this->_tetConstraints.push_back(createTetrahedralConstraint(
        this->_constraintId++,
        w,
        this->_nodes[existingNodesCount + v1],
        this->_nodes[existingNodesCount + v2],
        this->_nodes[existingNodesCount + v3],
        this->_nodes[existingNodesCount + v4]));

    this->_volumeConstraints.push_back(createVolumeConstraint(
        this->_constraintId++,
        w,
        this->_nodes[existingNodesCount + v1],
        this->_nodes[existingNodesCount + v2],
        this->_nodes[existingNodesCount + v3],
        this->_nodes[existingNodesCount + v4]));
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

void Solver::createTetBox(
    const glm::vec3& translation,
    float scale,
    const glm::vec3& initialVelocity,
    float w,
    float mass,
    bool hinged) {
  Grid grid{3, 3, 3};

  if (hinged) {
    grid = Grid{10, 2, 10};
  }

  size_t currentNodeCount = this->_nodes.size();
  size_t currentTetConstraintsCount = this->_tetConstraints.size();
  size_t currentLinesCount = this->_lines.size();
  size_t currentTriCount = this->_triangles.size();

  glm::vec3 boxColor = randColor();
  float boxRoughness = randf();
  float boxMetallic = static_cast<float>(std::rand() % 2);

  // Add nodes in a grid
  this->_nodes.reserve(
      currentNodeCount + grid.width * grid.height * grid.depth);
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      for (uint32_t k = 0; k < grid.depth; ++k) {
        uint32_t nodeId = grid.gridIdToNodeId(currentNodeCount, {i, j, k});

        Node& node = this->_nodes.emplace_back();
        node.id = nodeId;
        node.position = scale * glm::vec3(i, j, k) + translation;
        node.prevPosition = node.position;
        node.velocity = initialVelocity;
        node.radius = 0.95f * 0.5f * scale;
        node.invMass = 1.0f / mass;

        // if (hinged && i == 0) { //} && j == 0) {
        //   this->_positionConstraints.push_back(
        //       createPositionConstraint(this->_constraintId++, node,
        //       stiffness));
        // }
      }
    }
  }

  // Add tetrahedral constraints in each grid cell
  this->_tetConstraints.reserve(
      currentTetConstraintsCount +
      6 * (grid.width - 1) * (grid.height - 1) * (grid.depth - 1));
  this->_tets.reserve(
      this->_tets.size() +
      6 * (grid.width - 1) * (grid.height - 1) * (grid.depth - 1));
  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t j = 0; j < grid.height - 1; ++j) {
      for (uint32_t k = 0; k < grid.depth - 1; ++k) {
        uint32_t node000 = grid.gridIdToNodeId(currentNodeCount, {i, j, k});
        uint32_t node001 = grid.gridIdToNodeId(currentNodeCount, {i, j, k + 1});
        uint32_t node010 = grid.gridIdToNodeId(currentNodeCount, {i, j + 1, k});
        uint32_t node011 =
            grid.gridIdToNodeId(currentNodeCount, {i, j + 1, k + 1});
        uint32_t node100 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j, k});
        uint32_t node101 =
            grid.gridIdToNodeId(currentNodeCount, {i + 1, j, k + 1});
        uint32_t node110 =
            grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, k});
        uint32_t node111 =
            grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, k + 1});

        // TODO: consolidate tets and tet constraints
        this->_tetConstraints.push_back(createTetrahedralConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node001],
            this->_nodes[node011],
            this->_nodes[node111]));
        this->_volumeConstraints.push_back(createVolumeConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node001],
            this->_nodes[node011],
            this->_nodes[node111]));
        this->_tets.push_back(
            {{this->_nodes[node000].id,
              this->_nodes[node001].id,
              this->_nodes[node011].id,
              this->_nodes[node111].id}});
        this->_tetConstraints.push_back(createTetrahedralConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node010],
            this->_nodes[node011],
            this->_nodes[node111]));
        this->_volumeConstraints.push_back(createVolumeConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node010],
            this->_nodes[node011],
            this->_nodes[node111]));
        this->_tets.push_back(
            {{this->_nodes[node000].id,
              this->_nodes[node010].id,
              this->_nodes[node011].id,
              this->_nodes[node111].id}});
        this->_tetConstraints.push_back(createTetrahedralConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node001],
            this->_nodes[node101],
            this->_nodes[node111]));
        this->_volumeConstraints.push_back(createVolumeConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node001],
            this->_nodes[node101],
            this->_nodes[node111]));
        this->_tets.push_back(
            {{this->_nodes[node000].id,
              this->_nodes[node001].id,
              this->_nodes[node101].id,
              this->_nodes[node111].id}});
        this->_tetConstraints.push_back(createTetrahedralConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node100],
            this->_nodes[node101],
            this->_nodes[node111]));
        this->_volumeConstraints.push_back(createVolumeConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node100],
            this->_nodes[node101],
            this->_nodes[node111]));
        this->_tets.push_back(
            {{this->_nodes[node000].id,
              this->_nodes[node100].id,
              this->_nodes[node101].id,
              this->_nodes[node111].id}});
        this->_tetConstraints.push_back(createTetrahedralConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node010],
            this->_nodes[node110],
            this->_nodes[node111]));
        this->_volumeConstraints.push_back(createVolumeConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node010],
            this->_nodes[node110],
            this->_nodes[node111]));
        this->_tets.push_back(
            {{this->_nodes[node000].id,
              this->_nodes[node010].id,
              this->_nodes[node110].id,
              this->_nodes[node111].id}});
        this->_tetConstraints.push_back(createTetrahedralConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node100],
            this->_nodes[node110],
            this->_nodes[node111]));
        this->_volumeConstraints.push_back(createVolumeConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node000],
            this->_nodes[node100],
            this->_nodes[node110],
            this->_nodes[node111]));
        this->_tets.push_back(
            {{this->_nodes[node000].id,
              this->_nodes[node100].id,
              this->_nodes[node110].id,
              this->_nodes[node111].id}});
      }
    }
  }

  this->_triangles.reserve(
      currentTriCount +
      3 * 2 *
          (2 * grid.width * grid.height + 2 * grid.width * grid.depth +
           2 * grid.height * grid.depth));
  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t j = 0; j < grid.height - 1; ++j) {
      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, grid.depth - 1}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j, grid.depth - 1}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {i + 1, j + 1, grid.depth - 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, grid.depth - 1}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {i + 1, j + 1, grid.depth - 1}),
           grid.gridIdToNodeId(currentNodeCount, {i, j + 1, grid.depth - 1})});
    }
  }

  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t k = 0; k < grid.depth - 1; ++k) {
      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, 0, k}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k + 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, 0, k}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {i, 0, k + 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {i + 1, grid.height - 1, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, grid.height - 1, k})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k}),
           grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k + 1}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {i + 1, grid.height - 1, k + 1})});
    }
  }

  for (uint32_t j = 0; j < grid.height - 1; ++j) {
    for (uint32_t k = 0; k < grid.depth - 1; ++k) {
      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {0, j, k}),
           grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {0, j, k}),
           grid.gridIdToNodeId(currentNodeCount, {0, j, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k + 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k}),
           grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j + 1, k}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {grid.width - 1, j + 1, k + 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {grid.width - 1, j + 1, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k + 1})});
    }
  }

  this->_vertices.resize(this->_nodes.size());
  for (size_t i = currentNodeCount; i < this->_nodes.size(); ++i) {
    this->_vertices[i].position = this->_nodes[i].position;
    this->_vertices[i].radius = this->_nodes[i].radius;
    this->_vertices[i].baseColor = boxColor;
    this->_vertices[i].roughness = boxRoughness;
    this->_vertices[i].metallic = boxMetallic;
  }

  this->renderStateDirty = true;
}

void Solver::createBox(const glm::vec3& translation, float scale, float w) {
  Grid grid{5, 5, 5};

  size_t currentNodeCount = this->_nodes.size();
  size_t currentDistConstraintsCount = this->_distanceConstraints.size();
  size_t currentLinesCount = this->_lines.size();
  size_t currentTriCount = this->_triangles.size();

  glm::vec3 boxColor = randColor();
  float boxRoughness = randf();
  float boxMetallic = static_cast<float>(std::rand() % 2);

  // Add nodes in a grid
  this->_nodes.reserve(
      currentNodeCount + grid.width * grid.height * grid.depth);
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      for (uint32_t k = 0; k < grid.depth; ++k) {
        uint32_t nodeId = grid.gridIdToNodeId(currentNodeCount, {i, j, k});

        Node& node = this->_nodes.emplace_back();
        node.id = nodeId;
        node.position = scale * glm::vec3(i, j, k) + translation;
        node.prevPosition = node.position;
        node.velocity = glm::vec3(0.0f);
        node.radius = 0.5f * scale;
        node.invMass = 1.0f;

        // if (i == 0 && j == 0) {
        //   this->_positionConstraints.push_back(
        //       createPositionConstraint(this->_constraintId++, &node));
        // }
      }
    }
  }

  // Add distance constraints in each grid cell
  this->_distanceConstraints.reserve(
      currentDistConstraintsCount +
      8 * (grid.width - 1) * (grid.height - 1) * (grid.depth - 1));
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      for (uint32_t k = 0; k < grid.depth; ++k) {
        uint32_t node000 = grid.gridIdToNodeId(currentNodeCount, {i, j, k});
        uint32_t node001 = grid.gridIdToNodeId(currentNodeCount, {i, j, k + 1});
        uint32_t node010 = grid.gridIdToNodeId(currentNodeCount, {i, j + 1, k});
        uint32_t node011 =
            grid.gridIdToNodeId(currentNodeCount, {i, j + 1, k + 1});
        uint32_t node100 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j, k});
        uint32_t node101 =
            grid.gridIdToNodeId(currentNodeCount, {i + 1, j, k + 1});
        uint32_t node110 =
            grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, k});
        uint32_t node111 =
            grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, k + 1});

        // Grid-aligned constraints
        if (i < (grid.width - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              this->_nodes[node000],
              this->_nodes[node100],
              w));
        }

        if (j < (grid.height - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              this->_nodes[node000],
              this->_nodes[node010],
              w));
        }

        if (k < (grid.depth - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              this->_nodes[node000],
              this->_nodes[node001],
              w));
        }

        // Long diagonal constraints
        if (i < (grid.width - 1) && j < (grid.height - 1) &&
            k < (grid.depth - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              this->_nodes[node000],
              this->_nodes[node111],
              w));
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              this->_nodes[node100],
              this->_nodes[node011],
              w));
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              this->_nodes[node010],
              this->_nodes[node101],
              w));
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              this->_nodes[node001],
              this->_nodes[node110],
              w));
        }
      }
    }
  }

  this->_triangles.reserve(
      currentTriCount +
      3 * 2 *
          (2 * grid.width * grid.height + 2 * grid.width * grid.depth +
           2 * grid.height * grid.depth));
  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t j = 0; j < grid.height - 1; ++j) {
      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, grid.depth - 1}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j, grid.depth - 1}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {i + 1, j + 1, grid.depth - 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, grid.depth - 1}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {i + 1, j + 1, grid.depth - 1}),
           grid.gridIdToNodeId(currentNodeCount, {i, j + 1, grid.depth - 1})});
    }
  }

  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t k = 0; k < grid.depth - 1; ++k) {
      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, 0, k}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k + 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, 0, k}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {i, 0, k + 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {i + 1, grid.height - 1, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, grid.height - 1, k})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k}),
           grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k + 1}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {i + 1, grid.height - 1, k + 1})});
    }
  }

  for (uint32_t j = 0; j < grid.height - 1; ++j) {
    for (uint32_t k = 0; k < grid.depth - 1; ++k) {
      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {0, j, k}),
           grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {0, j, k}),
           grid.gridIdToNodeId(currentNodeCount, {0, j, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k + 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k}),
           grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j + 1, k}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {grid.width - 1, j + 1, k + 1})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k}),
           grid.gridIdToNodeId(
               currentNodeCount,
               {grid.width - 1, j + 1, k + 1}),
           grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k + 1})});
    }
  }

  // TODO: This is only for PBD
  // for (size_t i = currentDistConstraintsCount;
  //      i < this->_distanceConstraints.size();
  //      ++i) {
  //   this->_distanceConstraints[i].setWeight(
  //       1.0f - powf(1.0f - w, 1.0f / this->_options.iterations));
  // }

  this->_lines.reserve(
      currentLinesCount +
      2 * (this->_distanceConstraints.size() - currentDistConstraintsCount));
  for (size_t i = currentDistConstraintsCount;
       i < this->_distanceConstraints.size();
       ++i) {
    const DistanceConstraint& constraint = this->_distanceConstraints[i];
    this->_lines.push_back(constraint.getNodeId(0));
    this->_lines.push_back(constraint.getNodeId(1));
  }

  this->_vertices.resize(this->_nodes.size());
  for (size_t i = currentNodeCount; i < this->_nodes.size(); ++i) {
    this->_vertices[i].position = this->_nodes[i].position;
    this->_vertices[i].radius = this->_nodes[i].radius;
    this->_vertices[i].baseColor = boxColor;
    this->_vertices[i].roughness = boxRoughness;
    this->_vertices[i].metallic = boxMetallic;
  }

  this->renderStateDirty = true;
}

void Solver::createSheet(
    const glm::vec3& translation,
    float scale,
    float mass,
    float w) {
  Grid grid{10, 10, 1};

  size_t currentNodeCount = this->_nodes.size();
  size_t currentDistConstraintsCount = this->_distanceConstraints.size();
  size_t currentLinesCount = this->_lines.size();
  size_t currentTriCount = this->_triangles.size();

  glm::vec3 boxColor = randColor();
  float boxRoughness = randf();
  float boxMetallic = static_cast<float>(std::rand() % 2);

  // Add nodes in a grid
  this->_nodes.reserve(currentNodeCount + grid.width * grid.height);
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      uint32_t nodeId = grid.gridIdToNodeId(currentNodeCount, {i, j, 0});

      Node& node = this->_nodes.emplace_back();
      node.id = nodeId;
      node.position = scale * glm::vec3(i, 0, j) + translation;
      node.prevPosition = node.position;
      node.velocity = glm::vec3(0.0f);
      node.radius = 0.5f * scale;
      node.invMass = 1.0f / mass;

      if (i == 0 || i == (grid.width - 1) || j == 0 || j == (grid.height - 1)) {
        this->_positionConstraints.push_back(
            createPositionConstraint(this->_constraintId++, node, w));
      }
    }
  }

  // Add distance constraints in each grid cell
  this->_distanceConstraints.reserve(
      currentDistConstraintsCount + 4 * (grid.width - 1) * (grid.height - 1));
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      uint32_t node00 = grid.gridIdToNodeId(currentNodeCount, {i, j, 0});
      uint32_t node01 = grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0});
      uint32_t node10 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0});
      uint32_t node11 =
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0});

      // Grid-aligned constraints
      if (i < (grid.width - 1)) {
        this->_distanceConstraints.push_back(createDistanceConstraint(
            this->_constraintId++,
            this->_nodes[node00],
            this->_nodes[node10],
            w));
      }

      if (j < (grid.height - 1)) {
        this->_distanceConstraints.push_back(createDistanceConstraint(
            this->_constraintId++,
            this->_nodes[node00],
            this->_nodes[node01],
            w));
      }

      // Long diagonal constraints
      if (i < (grid.width - 1) && j < (grid.height - 1)) {
        this->_distanceConstraints.push_back(createDistanceConstraint(
            this->_constraintId++,
            this->_nodes[node00],
            this->_nodes[node11],
            w));
        this->_distanceConstraints.push_back(createDistanceConstraint(
            this->_constraintId++,
            this->_nodes[node10],
            this->_nodes[node01],
            w));
      }

      // TODO: Add constraints for other diagonal
    }
  }

  this->_triangles.reserve(currentTriCount + (2 * grid.width * grid.height));
  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t j = 0; j < grid.height - 1; ++j) {
      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0})});
    }
  }

  // TODO: This should only be for PBD
  // for (size_t i = currentDistConstraintsCount;
  //      i < this->_distanceConstraints.size();
  //      ++i) {
  //   this->_distanceConstraints[i].setWeight(
  //       1.0f - powf(1.0f - k, 1.0f / this->_options.iterations));
  // }

  this->_lines.reserve(
      currentLinesCount +
      2 * (this->_distanceConstraints.size() - currentDistConstraintsCount));
  for (size_t i = currentDistConstraintsCount;
       i < this->_distanceConstraints.size();
       ++i) {
    const DistanceConstraint& constraint = this->_distanceConstraints[i];
    this->_lines.push_back(constraint.getNodeId(0));
    this->_lines.push_back(constraint.getNodeId(1));
  }

  this->_vertices.resize(this->_nodes.size());
  for (size_t i = currentNodeCount; i < this->_nodes.size(); ++i) {
    this->_vertices[i].position = this->_nodes[i].position;
    this->_vertices[i].radius = this->_nodes[i].radius;
    this->_vertices[i].baseColor = boxColor;
    this->_vertices[i].roughness = boxRoughness;
    this->_vertices[i].metallic = boxMetallic;
  }

  this->renderStateDirty = true;
}

namespace {
struct ShapeMatchingPatch {
  std::vector<glm::vec3> materialCoords;
  std::vector<uint32_t> indices;
};
} // namespace

void Solver::createShapeMatchingBox(
    const glm::vec3& translation,
    uint32_t countX,
    uint32_t countY,
    uint32_t countZ,
    float scale,
    const glm::vec3& initialVelocity,
    float w) {
  Grid grid{countX, countY, countZ};

  scale = 0.5f;

  size_t currentNodeCount = this->_nodes.size();

  glm::vec3 boxColor = randColor();
  float boxRoughness = randf();
  float boxMetallic = static_cast<float>(std::rand() % 2);

  // Add nodes in a grid
  this->_nodes.reserve(
      currentNodeCount + grid.width * grid.height * grid.depth);
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      for (uint32_t k = 0; k < grid.depth; ++k) {
        uint32_t nodeId = grid.gridIdToNodeId(currentNodeCount, {i, j, k});

        Node& node = this->_nodes.emplace_back();
        node.id = nodeId;
        node.position = scale * glm::vec3(i, j, k) + translation;
        node.prevPosition = node.position;
        node.velocity = glm::vec3(0.0f);
        node.radius = 0.5f * scale;
        node.invMass = 1.0f / 10.0f;

        // if (i == 0 && j == 0) {
        //   this->_positionConstraints.push_back(
        //       createPositionConstraint(this->_constraintId++, &node));
        // }
      }
    }
  }

  std::vector<uint32_t> nodeIndices(grid.width * grid.height * grid.depth);
  std::vector<glm::vec3> materialCoords(nodeIndices.size());
  for (uint32_t i = 0; i < nodeIndices.size(); ++i) {
    uint32_t nodeId = static_cast<uint32_t>(currentNodeCount) + i;
    nodeIndices[i] = nodeId;
    materialCoords[i] = this->_nodes[nodeId].position;
  }

  this->_shapeConstraints.emplace_back(nodeIndices, materialCoords, w);

  this->_vertices.resize(this->_nodes.size());
  for (size_t i = currentNodeCount; i < this->_nodes.size(); ++i) {
    this->_vertices[i].position = this->_nodes[i].position;
    this->_vertices[i].radius = this->_nodes[i].radius;
    this->_vertices[i].baseColor = boxColor;
    this->_vertices[i].roughness = boxRoughness;
    this->_vertices[i].metallic = boxMetallic;
  }

  this->renderStateDirty = true;
}

void Solver::createShapeMatchingSheet(
    const glm::vec3& translation,
    float scale,
    const glm::vec3& initialVelocity,
    float w) {
  Grid grid{50, 50, 1};

  size_t currentNodeCount = this->_nodes.size();

  glm::vec3 boxColor = randColor();
  float boxRoughness = randf();
  float boxMetallic = static_cast<float>(std::rand() % 2);

  constexpr uint32_t patchWidth = 3;
  constexpr uint32_t patchHeight = 3;
  std::vector<ShapeMatchingPatch> patches(
      (grid.width / patchWidth) * (grid.height / patchHeight));

  for (ShapeMatchingPatch& patch : patches) {
    patch.indices.reserve(patchWidth * patchHeight);
    patch.materialCoords.reserve(patchWidth * patchHeight);
  }

  // Add nodes in a grid
  this->_nodes.reserve(
      currentNodeCount + grid.width * grid.height * grid.depth);
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      for (uint32_t k = 0; k < grid.depth; ++k) {
        uint32_t nodeId = grid.gridIdToNodeId(currentNodeCount, {i, j, k});

        Node& node = this->_nodes.emplace_back();
        node.id = nodeId;
        node.position = scale * glm::vec3(i, j, k) + translation;
        node.prevPosition = node.position;
        node.velocity = glm::vec3(0.0f);
        node.radius = 0.5f * scale;
        node.invMass = 1.0f;

        uint32_t patchId = i / patchWidth * patchHeight + j / patchHeight;
        patches[patchId].materialCoords.push_back(node.position);
        patches[patchId].indices.push_back(nodeId);

        if ((i % patchWidth) == (patchWidth - 1) && i < (grid.width - 1)) {
          patchId = (1 + i / patchWidth) * patchHeight + j / patchHeight;
          patches[patchId].materialCoords.push_back(node.position);
          patches[patchId].indices.push_back(nodeId);
        }

        if ((j % patchHeight) == (patchHeight - 1) && j < (grid.height - 1)) {
          patchId = i / patchWidth * patchHeight + j / patchHeight + 1;
          patches[patchId].materialCoords.push_back(node.position);
          patches[patchId].indices.push_back(nodeId);
        }
      }
    }
  }

  this->_shapeConstraints.reserve(
      this->_shapeConstraints.size() + patchWidth * patchHeight);
  for (const ShapeMatchingPatch& patch : patches) {
    this->_shapeConstraints.emplace_back(
        patch.indices,
        patch.materialCoords,
        w);
  }

  this->_vertices.resize(this->_nodes.size());
  for (size_t i = currentNodeCount; i < this->_nodes.size(); ++i) {
    this->_vertices[i].position = this->_nodes[i].position;
    this->_vertices[i].radius = this->_nodes[i].radius;
    this->_vertices[i].baseColor = boxColor;
    this->_vertices[i].roughness = boxRoughness;
    this->_vertices[i].metallic = boxMetallic;
  }

  this->renderStateDirty = true;
}

void Solver::createBendSheet(
    const glm::vec3& translation,
    float scale,
    float w) {
  Grid grid{10, 10, 1};

  glm::vec3 boxColor = randColor();
  float boxRoughness = randf();
  float boxMetallic = static_cast<float>(std::rand() % 2);

  size_t currentNodeCount = this->_nodes.size();
  size_t currentDistConstraintsCount = this->_distanceConstraints.size();
  size_t currentBendConstraintsCount = this->_bendConstraints.size();
  size_t currentLinesCount = this->_lines.size();
  size_t currentTriCount = this->_triangles.size();
  this->_nodes.reserve(currentNodeCount + grid.width * grid.height);
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      uint32_t nodeId = grid.gridIdToNodeId(currentNodeCount, {i, j, 0});

      Node& node = this->_nodes.emplace_back();
      node.id = nodeId;
      node.position = scale * glm::vec3(i, 0, j) + translation;
      node.prevPosition = node.position;
      node.velocity = glm::vec3(0.0f);
      node.radius = 0.5f * scale;
      node.invMass = 1.0f;

      if (i < 3) {
        this->_positionConstraints.push_back(
            createPositionConstraint(this->_constraintId++, node, w));
      }
    }
  }

  // Add distance constraints in each grid cell
  this->_distanceConstraints.reserve(
      currentDistConstraintsCount + 8 * (grid.width - 1) * (grid.height - 1));
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      uint32_t node00 = grid.gridIdToNodeId(currentNodeCount, {i, j, 0});
      uint32_t node01 = grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0});
      uint32_t node10 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0});
      uint32_t node11 =
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0});

      // Grid-aligned constraints
      if (i < (grid.width - 1)) {
        this->_distanceConstraints.push_back(createDistanceConstraint(
            this->_constraintId++,
            this->_nodes[node00],
            this->_nodes[node10],
            w));
      }

      if (j < (grid.height - 1)) {
        this->_distanceConstraints.push_back(createDistanceConstraint(
            this->_constraintId++,
            this->_nodes[node00],
            this->_nodes[node01],
            w));
      }

      // Long diagonal constraints
      if (i < (grid.width - 1) && j < (grid.height - 1)) {
        this->_distanceConstraints.push_back(createDistanceConstraint(
            this->_constraintId++,
            this->_nodes[node00],
            this->_nodes[node11],
            w));
      }
    }
  }

  // TODO: HARD-CODED FOR GRID SIZE
  // Add distance constraints in each grid cell
  this->_bendConstraints.reserve(
      currentBendConstraintsCount + 2 * (grid.width - 1) * (grid.height - 1));
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
      uint32_t node00 = grid.gridIdToNodeId(currentNodeCount, {i, j, 0});
      uint32_t node01 = grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0});
      uint32_t node10 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0});
      uint32_t node11 =
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0});

      // Diagonal constraints
      if (i < (grid.width - 1) && j < (grid.height - 1)) {
        this->_bendConstraints.push_back(createBendConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node00],
            this->_nodes[node11],
            this->_nodes[node10],
            this->_nodes[node01]));
      }

      if (i < (grid.width - 2) && j < (grid.height - 2)) {
        uint32_t node02 = grid.gridIdToNodeId(currentNodeCount, {i, j + 2, 0});
        uint32_t node12 =
            grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 2, 0});
        uint32_t node20 = grid.gridIdToNodeId(currentNodeCount, {i + 2, j, 0});
        uint32_t node21 =
            grid.gridIdToNodeId(currentNodeCount, {i + 2, j + 1, 0});

        // shared edge to the right of current grid square
        this->_bendConstraints.push_back(createBendConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node10],
            this->_nodes[node11],
            this->_nodes[node00],
            this->_nodes[node21]));

        // shared edge to the right of current grid square
        this->_bendConstraints.push_back(createBendConstraint(
            this->_constraintId++,
            w,
            this->_nodes[node01],
            this->_nodes[node11],
            this->_nodes[node00],
            this->_nodes[node12]));
      }
    }
  }

  this->_triangles.reserve(currentTriCount + (2 * grid.width * grid.height));
  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t j = 0; j < grid.height - 1; ++j) {
      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0})});

      this->_triangles.push_back(
          {grid.gridIdToNodeId(currentNodeCount, {i, j, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0}),
           grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0})});
    }
  }

  this->_lines.reserve(
      currentLinesCount +
      2 * (this->_distanceConstraints.size() - currentDistConstraintsCount));
  for (size_t i = currentDistConstraintsCount;
       i < this->_distanceConstraints.size();
       ++i) {
    const DistanceConstraint& constraint = this->_distanceConstraints[i];
    this->_lines.push_back(constraint.getNodeId(0));
    this->_lines.push_back(constraint.getNodeId(1));
  }

  this->_vertices.resize(this->_nodes.size());
  for (size_t i = currentNodeCount; i < this->_nodes.size(); ++i) {
    this->_vertices[i].position = this->_nodes[i].position;
    this->_vertices[i].radius = this->_nodes[i].radius;
    this->_vertices[i].baseColor = boxColor;
    this->_vertices[i].roughness = boxRoughness;
    this->_vertices[i].metallic = boxMetallic;
  }

  this->renderStateDirty = true;
}
} // namespace Pies