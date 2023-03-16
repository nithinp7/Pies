#include "Solver.h"

#include <glm/gtc/matrix_transform.hpp>

namespace Pies {

Solver::Solver(const SolverOptions& options) : 
    _options(options),
    _spatialHashNodes(glm::scale(glm::mat4(1.0f), glm::vec3(0.01f))) {}

void Solver::tick(float /*timestep*/) {
  float deltaTime =
      this->_options.fixedTimestepSize / this->_options.timeSubsteps;

  // Time substeps
  for (uint32_t substep = 0; substep < this->_options.timeSubsteps; ++substep) {
    // Apply external forces and advect nodes
    for (Node& node : this->_nodes) {
      node.prevPosition = node.position;
      node.position += node.velocity * deltaTime +
                       glm::vec3(0.0f, -this->_options.gravity, 0.0f) *
                           deltaTime * deltaTime;
    }

    // TODO: Collision solver
    // this->_spatialHashNodes.clear();
    // this->_spatialHashNodes.parallelBulkInsert(this->_nodes);

    for (uint32_t i = 0; i < this->_options.iterations; ++i) {
      if (!releaseHinge) {
        for (PositionConstraint& constraint : this->_positionConstraints) {
          constraint.projectNodePositions(this->_nodes);
        }
      }

      for (DistanceConstraint& constraint : this->_distanceConstraints) {
        constraint.projectNodePositions(this->_nodes);
      }

      for (TetrahedralConstraint& constraint : this->_tetConstraints) {
        constraint.projectNodePositions(this->_nodes);
      }

      // TODO: Collision solver
      this->_spatialHashNodes.clear();
      this->_spatialHashNodes.parallelBulkInsert(this->_nodes);
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
          (this->_nodes[i].position - this->_nodes[i].prevPosition) / deltaTime;

      // TODO: friction for dynamically generated constraints
      if (this->_nodes[i].position.y <= -8.0f) {
        this->_nodes[i].velocity.x *= 1.0f - this->_options.friction;
        this->_nodes[i].velocity.z *= 1.0f - this->_options.friction;
      }

      this->_vertices[i].position = this->_nodes[i].position;
    }
  }
}

SpatialHashGridCellRange
Solver::NodeCompRange::operator()(const Node& node) const {
  glm::vec3 gridLocalPos =
      glm::vec3(this->grid.worldToGrid * glm::vec4(node.position, 1.0f));

  SpatialHashGridCellRange range{};
  range.minX = static_cast<int64_t>(gridLocalPos.x);
  range.minY = static_cast<int64_t>(gridLocalPos.y);
  range.minZ = static_cast<int64_t>(gridLocalPos.z);

  // Assume infinitely small
  range.lengthX = 1;
  range.lengthY = 1;
  range.lengthZ = 1;

  return range;
}
} // namespace Pies