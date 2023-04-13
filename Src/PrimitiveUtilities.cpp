#include "Node.h"
#include "Solver.h"

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

void Solver::addNodes(
    const std::vector<glm::vec3>& vertices) {
  
  // TODO: parameterize more of these
  float mass = 1.0f;
  float radius = 0.5f;
  glm::vec3 initialVelocity(0.0f);

  glm::vec3 color = randColor();
  float roughness = randf();
  float metallic = static_cast<float>(std::rand() % 2);

  size_t currentNodeCount = this->_nodes.size();
  this->_nodes.reserve(
      currentNodeCount + vertices.size());
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

void Solver::createTetBox(
    const glm::vec3& translation,
    float scale,
    const glm::vec3& initialVelocity,
    float w,
    float mass,
    bool hinged) {
  Grid grid{3, 3, 3};

  if (hinged) {
    grid = Grid{15, 2, 2};
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
        node.radius = 0.5f * scale;
        node.invMass = 1.0f / mass;

        if (hinged && i == 0) { //} && j == 0) {
          this->_positionConstraints.push_back(
              createPositionConstraint(this->_constraintId++, node, w));
        }
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
          grid.gridIdToNodeId(currentNodeCount, {i, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j, grid.depth - 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j, grid.depth - 1}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {i + 1, j + 1, grid.depth - 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j, grid.depth - 1}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {i + 1, j + 1, grid.depth - 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j + 1, grid.depth - 1}));
    }
  }

  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t k = 0; k < grid.depth - 1; ++k) {
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, 0, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, 0, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k + 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, 0, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, grid.height - 1, k}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {i + 1, grid.height - 1, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {i + 1, grid.height - 1, k + 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k + 1}));
    }
  }

  for (uint32_t j = 0; j < grid.height - 1; ++j) {
    for (uint32_t k = 0; k < grid.depth - 1; ++k) {
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k + 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j + 1, k}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {grid.width - 1, j + 1, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {grid.width - 1, j + 1, k + 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k + 1}));
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

void Solver::createBox(
    const glm::vec3& translation,
    float scale,
    float w) {
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
          grid.gridIdToNodeId(currentNodeCount, {i, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j, grid.depth - 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j, grid.depth - 1}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {i + 1, j + 1, grid.depth - 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j, grid.depth - 1}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {i + 1, j + 1, grid.depth - 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j + 1, grid.depth - 1}));
    }
  }

  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t k = 0; k < grid.depth - 1; ++k) {
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, 0, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, 0, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, 0, k + 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, 0, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, grid.height - 1, k}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {i + 1, grid.height - 1, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {i + 1, grid.height - 1, k + 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, grid.height - 1, k + 1}));
    }
  }

  for (uint32_t j = 0; j < grid.height - 1; ++j) {
    for (uint32_t k = 0; k < grid.depth - 1; ++k) {
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j + 1, k + 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {0, j, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j + 1, k}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {grid.width - 1, j + 1, k + 1}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k}));
      this->_triangles.push_back(grid.gridIdToNodeId(
          currentNodeCount,
          {grid.width - 1, j + 1, k + 1}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {grid.width - 1, j, k + 1}));
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
  Grid grid{40, 40, 1};

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
          grid.gridIdToNodeId(currentNodeCount, {i, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0}));
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

void Solver::createBendSheet(const glm::vec3& translation, float scale, float w) {
  Grid grid{10, 10, 1};

  size_t currentNodeCount = this->_nodes.size();
  size_t currentDistConstraintsCount = this->_distanceConstraints.size();
  size_t currentBendConstraintsCount = this->_bendConstraints.size();
  size_t currentLinesCount = this->_lines.size();
  size_t currentTriCount = this->_triangles.size();

  glm::vec3 boxColor = randColor();
  float boxRoughness = randf();
  float boxMetallic = static_cast<float>(std::rand() % 2);

  // Add nodes in a grid
  this->_nodes.reserve(
      currentNodeCount + grid.width * grid.height);
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
      currentDistConstraintsCount +
      8 * (grid.width - 1) * (grid.height - 1));
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
        uint32_t node00 = grid.gridIdToNodeId(currentNodeCount, {i, j, 0});
        uint32_t node01 = grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0});
        uint32_t node10 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0});
        uint32_t node11 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0});

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
      currentBendConstraintsCount +
      2 * (grid.width - 1) * (grid.height - 1));
  for (uint32_t i = 0; i < grid.width; ++i) {
    for (uint32_t j = 0; j < grid.height; ++j) {
        uint32_t node00 = grid.gridIdToNodeId(currentNodeCount, {i, j, 0});
        uint32_t node01 = grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0});
        uint32_t node10 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0});
        uint32_t node11 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0});

        // Diagonal constraints
        if (i < (grid.width - 1) && j < (grid.height - 1)) {
            this->_bendConstraints.push_back(createBendConstraint(
                this->_constraintId++,
                w,
                this->_nodes[node00],
                this->_nodes[node11],
                this->_nodes[node10],
                this->_nodes[node01]
                ));
        }

        if (i < (grid.width - 2) && j < (grid.height - 2)) {
            uint32_t node02 = grid.gridIdToNodeId(currentNodeCount, {i, j + 2, 0});
            uint32_t node12 = grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 2, 0});
            uint32_t node20 = grid.gridIdToNodeId(currentNodeCount, {i + 2, j, 0});
            uint32_t node21 = grid.gridIdToNodeId(currentNodeCount, {i + 2, j + 1, 0});

            //shared edge to the right of current grid square
            this->_bendConstraints.push_back(createBendConstraint(
                this->_constraintId++,
                w,
                this->_nodes[node10],
                this->_nodes[node11],
                this->_nodes[node00],
                this->_nodes[node21]
            ));

            //shared edge to the right of current grid square
            this->_bendConstraints.push_back(createBendConstraint(
                this->_constraintId++,
                w,
                this->_nodes[node01],
                this->_nodes[node11],
                this->_nodes[node00],
                this->_nodes[node12]
            ));
        }
    }
  }

  this->_triangles.reserve(
      currentTriCount + (2 * grid.width * grid.height));
  for (uint32_t i = 0; i < grid.width - 1; ++i) {
    for (uint32_t j = 0; j < grid.height - 1; ++j) {
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}));

      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i + 1, j + 1, 0}));
      this->_triangles.push_back(
          grid.gridIdToNodeId(currentNodeCount, {i, j + 1, 0}));
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