#include "Solver.h"

#define GRID_WIDTH 10
#define GRID_HEIGHT 10
#define GRID_DEPTH 10

#define GRAVITY 10.0f

#define DAMPING 0.001f
#define FRICTION 0.1f

#define SOLVER_ITERS 4
#define SOLVER_SUBSTEPS 10

namespace Pies {
namespace {
struct GridId {
  uint32_t x;
  uint32_t y;
  uint32_t z;
};

GridId nodeIdToGridId(uint32_t nodeId) {
  GridId gridId;

  gridId.z = nodeId % GRID_DEPTH;
  gridId.y = (nodeId / GRID_DEPTH) % GRID_HEIGHT;
  gridId.x = (nodeId / GRID_DEPTH) / GRID_HEIGHT;

  return gridId;
}

uint32_t gridIdToNodeId(const GridId& gridId) {
  return gridId.z + GRID_DEPTH * (gridId.y + GRID_HEIGHT * gridId.x);
}
} // namespace

Solver::Solver() {
  glm::vec3 posOffs = glm::vec3(-10.0f, 5.0f, 0.0f);
  float scale = 0.5f;

  uint32_t constraintId = 0;

  // Add nodes in a grid
  this->_nodes.reserve(GRID_WIDTH * GRID_HEIGHT * GRID_DEPTH);
  for (uint32_t i = 0; i < GRID_WIDTH; ++i) {
    for (uint32_t j = 0; j < GRID_HEIGHT; ++j) {
      for (uint32_t k = 0; k < GRID_DEPTH; ++k) {
        uint32_t nodeId = gridIdToNodeId({i, j, k});

        Node& node = this->_nodes.emplace_back();
        node.id = nodeId;
        node.position = scale * glm::vec3(i, j, k) + posOffs;
        node.velocity = glm::vec3(0.0f);
        node.mass = 1.0f;

        if (i == 0 && j == 0) {
          this->_positionConstraints.push_back(
              createPositionConstraint(constraintId++, &node));
        }
      }
    }
  }

  // Add distance constraints in each grid cell
  this->_distanceConstraints.reserve(
      8 * (GRID_WIDTH - 1) * (GRID_HEIGHT - 1) * (GRID_DEPTH - 1));
  for (uint32_t i = 0; i < GRID_WIDTH; ++i) {
    for (uint32_t j = 0; j < GRID_HEIGHT; ++j) {
      for (uint32_t k = 0; k < GRID_DEPTH; ++k) {
        uint32_t node000 = gridIdToNodeId({i, j, k});
        uint32_t node001 = gridIdToNodeId({i, j, k + 1});
        uint32_t node010 = gridIdToNodeId({i, j + 1, k});
        uint32_t node011 = gridIdToNodeId({i, j + 1, k + 1});
        uint32_t node100 = gridIdToNodeId({i + 1, j, k});
        uint32_t node101 = gridIdToNodeId({i + 1, j, k + 1});
        uint32_t node110 = gridIdToNodeId({i + 1, j + 1, k});
        uint32_t node111 = gridIdToNodeId({i + 1, j + 1, k + 1});

        // Grid-aligned constraints
        if (i < (GRID_WIDTH - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              constraintId++,
              &this->_nodes[node000],
              &this->_nodes[node100]));
        }

        if (j < (GRID_HEIGHT - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              constraintId++,
              &this->_nodes[node000],
              &this->_nodes[node010]));
        }

        if (k < (GRID_DEPTH - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              constraintId++,
              &this->_nodes[node000],
              &this->_nodes[node001]));
        }

        // Long diagonal constraints
        if (i < (GRID_WIDTH - 1) && j < (GRID_HEIGHT - 1) &&
            k < (GRID_DEPTH - 1)) {
          this->_distanceConstraints.push_back(createDistanceConstraint(
              constraintId++,
              &this->_nodes[node000],
              &this->_nodes[node111]));
          this->_distanceConstraints.push_back(createDistanceConstraint(
              constraintId++,
              &this->_nodes[node100],
              &this->_nodes[node011]));
          this->_distanceConstraints.push_back(createDistanceConstraint(
              constraintId++,
              &this->_nodes[node010],
              &this->_nodes[node101]));
          this->_distanceConstraints.push_back(createDistanceConstraint(
              constraintId++,
              &this->_nodes[node001],
              &this->_nodes[node110]));
        }
      }
    }
  }

  this->_triangles.reserve(
      3 * 2 *
      (2 * GRID_WIDTH * GRID_HEIGHT + 2 * GRID_WIDTH * GRID_DEPTH +
       2 * GRID_HEIGHT * GRID_DEPTH));
  for (uint32_t i = 0; i < GRID_WIDTH - 1; ++i) {
    for (uint32_t j = 0; j < GRID_HEIGHT - 1; ++j) {
      this->_triangles.push_back(gridIdToNodeId({i, j, 0}));
      this->_triangles.push_back(gridIdToNodeId({i + 1, j, 0}));
      this->_triangles.push_back(gridIdToNodeId({i + 1, j + 1, 0}));

      this->_triangles.push_back(gridIdToNodeId({i, j, 0}));
      this->_triangles.push_back(gridIdToNodeId({i + 1, j + 1, 0}));
      this->_triangles.push_back(gridIdToNodeId({i, j + 1, 0}));

      this->_triangles.push_back(gridIdToNodeId({i, j, GRID_DEPTH - 1}));
      this->_triangles.push_back(gridIdToNodeId({i + 1, j, GRID_DEPTH - 1}));
      this->_triangles.push_back(
          gridIdToNodeId({i + 1, j + 1, GRID_DEPTH - 1}));

      this->_triangles.push_back(gridIdToNodeId({i, j, GRID_DEPTH - 1}));
      this->_triangles.push_back(
          gridIdToNodeId({i + 1, j + 1, GRID_DEPTH - 1}));
      this->_triangles.push_back(gridIdToNodeId({i, j + 1, GRID_DEPTH - 1}));
    }
  }

  for (uint32_t i = 0; i < GRID_WIDTH - 1; ++i) {
    for (uint32_t k = 0; k < GRID_DEPTH - 1; ++k) {
      this->_triangles.push_back(gridIdToNodeId({i, 0, k}));
      this->_triangles.push_back(gridIdToNodeId({i + 1, 0, k}));
      this->_triangles.push_back(gridIdToNodeId({i + 1, 0, k + 1}));

      this->_triangles.push_back(gridIdToNodeId({i, 0, k}));
      this->_triangles.push_back(gridIdToNodeId({i + 1, 0, k + 1}));
      this->_triangles.push_back(gridIdToNodeId({i, 0, k + 1}));

      this->_triangles.push_back(gridIdToNodeId({i, GRID_HEIGHT - 1, k}));
      this->_triangles.push_back(gridIdToNodeId({i + 1, GRID_HEIGHT - 1, k}));
      this->_triangles.push_back(
          gridIdToNodeId({i + 1, GRID_HEIGHT - 1, k + 1}));

      this->_triangles.push_back(gridIdToNodeId({i, GRID_HEIGHT - 1, k}));
      this->_triangles.push_back(
          gridIdToNodeId({i + 1, GRID_HEIGHT - 1, k + 1}));
      this->_triangles.push_back(gridIdToNodeId({i, GRID_HEIGHT - 1, k + 1}));
    }
  }

  for (uint32_t j = 0; j < GRID_HEIGHT - 1; ++j) {
    for (uint32_t k = 0; k < GRID_DEPTH - 1; ++k) {
      this->_triangles.push_back(gridIdToNodeId({0, j, k}));
      this->_triangles.push_back(gridIdToNodeId({0, j + 1, k}));
      this->_triangles.push_back(gridIdToNodeId({0, j + 1, k + 1}));

      this->_triangles.push_back(gridIdToNodeId({0, j, k}));
      this->_triangles.push_back(gridIdToNodeId({0, j + 1, k + 1}));
      this->_triangles.push_back(gridIdToNodeId({0, j, k + 1}));

      this->_triangles.push_back(gridIdToNodeId({GRID_WIDTH - 1, j, k}));
      this->_triangles.push_back(gridIdToNodeId({GRID_WIDTH - 1, j + 1, k}));
      this->_triangles.push_back(
          gridIdToNodeId({GRID_WIDTH - 1, j + 1, k + 1}));

      this->_triangles.push_back(gridIdToNodeId({GRID_WIDTH - 1, j, k}));
      this->_triangles.push_back(
          gridIdToNodeId({GRID_WIDTH - 1, j + 1, k + 1}));
      this->_triangles.push_back(gridIdToNodeId({GRID_WIDTH - 1, j, k + 1}));
    }
  }

  for (DistanceConstraint& constraint : this->_distanceConstraints) {
    float k = 0.85f;
    constraint.setWeight(1.0f - powf(1.0f - k, 1.0f / SOLVER_ITERS));
  }

  this->_distanceConstraintLines.resize(2 * this->_distanceConstraints.size());
  for (size_t i = 0; i < this->_distanceConstraints.size(); ++i) {
    const DistanceConstraint& constraint = this->_distanceConstraints[i];
    this->_distanceConstraintLines[2 * i] = constraint.getNode(0).id;
    this->_distanceConstraintLines[2 * i + 1] = constraint.getNode(1).id;
  }

  this->_vertices.resize(this->_nodes.size());
  for (uint32_t i = 0; i < this->_nodes.size(); ++i) {
    this->_vertices[i] = this->_nodes[i].position;
  }
}

void Solver::tick(float timestep) {
  float deltaTime = timestep / SOLVER_SUBSTEPS;

  // Time substeps
  for (int substep = 0; substep < SOLVER_SUBSTEPS; ++substep) {
    // Apply external forces and advect nodes
    for (Node& node : this->_nodes) {
      node.velocity += glm::vec3(0.0f, -GRAVITY, 0.0f) * deltaTime;
      node.position += node.velocity * deltaTime;
    }

    for (uint32_t i = 0; i < SOLVER_ITERS; ++i) {
      if (!releaseHinge) {
        for (PositionConstraint& constraint : this->_positionConstraints) {
          constraint.projectNodePositions();
        }
      }

      for (DistanceConstraint& constraint : this->_distanceConstraints) {
        constraint.projectNodePositions();
      }

      // Floor constraint
      for (Node& node : this->_nodes) {
        if (node.position.y < -8.0f) {
          node.position.y = -8.0f;
        }
      }
    }

    // Compute new velocity and construct new vertex positions
    for (uint32_t i = 0; i < this->_nodes.size(); ++i) {
      this->_nodes[i].velocity =
          (1.0f - DAMPING) * (this->_nodes[i].position - this->_vertices[i]) /
          deltaTime;

      if (this->_nodes[i].position.y <= -8.0f) {
        this->_nodes[i].velocity.x *= 1.0f - FRICTION;
        this->_nodes[i].velocity.z *= 1.0f - FRICTION;
      }

      this->_vertices[i] = this->_nodes[i].position;
    }
  }
}
} // namespace Pies