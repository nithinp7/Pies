#include "Solver.h"

namespace Pies {

Solver::Solver(const SolverOptions& options) : _options(options) {}

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

    for (uint32_t i = 0; i < this->_options.iterations; ++i) {
      if (!releaseHinge) {
        for (PositionConstraint& constraint : this->_positionConstraints) {
          constraint.projectNodePositions(this->_nodes);
        }
      }

      for (DistanceConstraint& constraint : this->_distanceConstraints) {
        constraint.projectNodePositions(this->_nodes);
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
} // namespace Pies