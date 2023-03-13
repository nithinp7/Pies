#include "Solver.h"

#define GRID_WIDTH 10
#define GRID_HEIGHT 10
#define GRID_DEPTH 10

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

Solver::Solver(const SolverOptions& options) : _options(options) {
  createBox(glm::vec3(-10.0f, 5.0f, 0.0f), 0.5f, 0.85f);
}

void Solver::tick(float /*timestep*/) {
  float deltaTime =
      this->_options.fixedTimestepSize / this->_options.timeSubsteps;

  // Time substeps
  for (int substep = 0; substep < this->_options.timeSubsteps; ++substep) {
    // Apply external forces and advect nodes
    for (Node& node : this->_nodes) {
      node.position += glm::vec3(0.0f, -this->_options.gravity, 0.0f) *
                       deltaTime * deltaTime;
    }

    for (uint32_t i = 0; i < this->_options.iterations; ++i) {
      if (!releaseHinge) {
        for (PositionConstraint& constraint : this->_positionConstraints) {
          constraint.projectNodePositions();
        }
      }

      for (DistanceConstraint& constraint : this->_distanceConstraints) {
        constraint.projectNodePositions();
      }

      // TODO: Collision solver
      // TODO: Apply projections to collision constraints

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
          (1.0f - this->_options.damping) *
          (this->_nodes[i].position - this->_vertices[i]) / deltaTime;

      // TODO: friction for dynamically generated constraints
      if (this->_nodes[i].position.y <= -8.0f) {
        this->_nodes[i].velocity.x *= 1.0f - this->_options.friction;
        this->_nodes[i].velocity.z *= 1.0f - this->_options.friction;
      }

      this->_vertices[i] = this->_nodes[i].position;
    }
  }
}
} // namespace Pies