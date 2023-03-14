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

// TODO: ... this is the worst hack yet...
// Need to rethink constraints keeping Node pointers, given that the Node list
// can grow and potentially reallocate. This fixes up the constraint node
// pointers, by computing their old and new indices.
template <typename TConstraint>
void fixupNodePointers(
    std::vector<TConstraint>& constraints,
    size_t previousConstraintCount,
    Node* pPrevNodeStart,
    std::vector<Node>& currentNodes) {
  for (size_t i = 0; i < previousConstraintCount; ++i) {
    TConstraint& constraint = constraints[i];
    for (Node*& pNode : constraint._nodes) {
      size_t index = static_cast<size_t>(pNode - pPrevNodeStart);
      pNode = &currentNodes[index];
    }
  }
}
} // namespace

void Solver::createBox(const glm::vec3& translation, float scale, float k) {
  Node* pPrevNodeStart = this->_nodes.data();

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
        node.mass = 1.0f;

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
              &this->_nodes[node000],
              &this->_nodes[node100]));
        }

        if (j < (grid.height - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              &this->_nodes[node000],
              &this->_nodes[node010]));
        }

        if (k < (grid.depth - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              &this->_nodes[node000],
              &this->_nodes[node001]));
        }

        // Long diagonal constraints
        if (i < (grid.width - 1) && j < (grid.height - 1) &&
            k < (grid.depth - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              &this->_nodes[node000],
              &this->_nodes[node111]));
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              &this->_nodes[node100],
              &this->_nodes[node011]));
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              &this->_nodes[node010],
              &this->_nodes[node101]));
          this->_distanceConstraints.push_back(createDistanceConstraint(
              this->_constraintId++,
              &this->_nodes[node001],
              &this->_nodes[node110]));
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

  for (size_t i = currentDistConstraintsCount;
       i < this->_distanceConstraints.size();
       ++i) {
    this->_distanceConstraints[i].setWeight(
        1.0f - powf(1.0f - k, 1.0f / this->_options.iterations));
  }

  this->_lines.reserve(
      currentLinesCount +
      2 * (this->_distanceConstraints.size() - currentDistConstraintsCount));
  for (size_t i = currentDistConstraintsCount;
       i < this->_distanceConstraints.size();
       ++i) {
    const DistanceConstraint& constraint = this->_distanceConstraints[i];
    this->_lines.push_back(constraint.getNode(0).id);
    this->_lines.push_back(constraint.getNode(1).id);
  }

  this->_vertices.resize(this->_nodes.size());
  for (uint32_t i = currentNodeCount; i < this->_nodes.size(); ++i) {
    this->_vertices[i].position = this->_nodes[i].position;
    this->_vertices[i].baseColor = boxColor;
    this->_vertices[i].roughness = boxRoughness;
    this->_vertices[i].metallic = boxMetallic;
  }

  // TODO: remove hack
  fixupNodePointers(
      this->_distanceConstraints,
      currentDistConstraintsCount,
      pPrevNodeStart,
      this->_nodes);

  this->renderStateDirty = true;
}
} // namespace Pies