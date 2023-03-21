#include "Solver.h"

#include <glm/gtc/matrix_transform.hpp>

namespace Pies {

Solver::Solver(const SolverOptions& options)
    : _options(options),
      _spatialHashNodes(glm::vec3(0.0f), 0.5f),
      _spatialHashTets(glm::vec3(0.0f), 0.5f) {}

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
    // this->_spatialHashNodes.parallelBulkInsert(this->_nodes, {});

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

      // // TODO: Collision solver
      // this->_spatialHashTets.clear();
      // this->_spatialHashTets.parallelBulkInsert(this->_tets, {this->_nodes});

      this->_spatialHashNodes.clear();
      this->_spatialHashNodes.parallelBulkInsert(this->_nodes, {});

      // Detect and resolve collisions
      std::vector<SpatialHashGridCellBucket<Node>*> scratchBuckets;
      for (Node& node : this->_nodes) {
        this->_spatialHashNodes.findCollisions(node, {}, scratchBuckets);
        for (SpatialHashGridCellBucket<Node>* pBucket : scratchBuckets) {
          // Check each node within each bucket
          for (Node* pOtherNode : pBucket->values) {
            // Check intersection
            glm::vec3 diff = pOtherNode->position - node.position;
            float dist = glm::length(diff);

            float disp = node.radius + pOtherNode->radius - dist;

            if (disp <= 0.0) {
              continue;
            }

            glm::vec3 dir(1.0f, 0.0f, 0.0f);
            if (dist > 0.00001f) {
              dir = diff / dist;
            }

            float massSum = node.mass + pOtherNode->mass;

            node.position += 0.85f * -disp * dir * node.mass / massSum;
            pOtherNode->position += 0.85f * disp * dir * pOtherNode->mass / massSum;
          }
        }

        scratchBuckets.clear();
      }

      for (Node& node : this->_nodes) {
        if (node.position.y - node.radius < this->_options.floorHeight) {
          node.position.y = this->_options.floorHeight + node.radius;
        }
      }
    }

    // Compute new velocity and construct new vertex positions
    for (uint32_t i = 0; i < this->_nodes.size(); ++i) {
      Node& node = this->_nodes[i];

      node.velocity =
          (1.0f - this->_options.damping) *
          (node.position - node.prevPosition) / deltaTime;

      // TODO: friction for dynamically generated constraints
      if (node.position.y - node.radius <= this->_options.floorHeight) {
        node.velocity.x *= 1.0f - this->_options.friction;
        node.velocity.z *= 1.0f - this->_options.friction;
      }

      this->_vertices[i].position = this->_nodes[i].position;
    }
  }
}

SpatialHashGridCellRange Solver::NodeCompRange::operator()(
    const Node& node,
    const SpatialHashGrid& grid) const {
  float radiusPadding = 0.0f;
  float gridLocalRadius = (node.radius + radiusPadding) / grid.scale;
  glm::vec3 gridLocalPos = (node.position - grid.translation) / grid.scale;
  glm::vec3 gridLocalMin = gridLocalPos - glm::vec3(gridLocalRadius);

  SpatialHashGridCellRange range{};
  range.minX = static_cast<int64_t>(glm::floor(gridLocalMin.x));
  range.minY = static_cast<int64_t>(glm::floor(gridLocalMin.y));
  range.minZ = static_cast<int64_t>(glm::floor(gridLocalMin.z));

  glm::vec3 adjustedDiameter =
      glm::round(glm::fract(gridLocalMin) + glm::vec3(2 * gridLocalRadius));
  range.lengthX = static_cast<uint32_t>(adjustedDiameter.r);
  range.lengthY = static_cast<uint32_t>(adjustedDiameter.g);
  range.lengthZ = static_cast<uint32_t>(adjustedDiameter.b);

  return range;
}

SpatialHashGridCellRange Solver::TetCompRange::operator()(
    const Tetrahedron& tet,
    const SpatialHashGrid& grid) const {
  glm::vec3 x1 =
      (nodes[tet.nodeIds[0]].position - grid.translation) / grid.scale;
  glm::vec3 x2 =
      (nodes[tet.nodeIds[1]].position - grid.translation) / grid.scale;
  glm::vec3 x3 =
      (nodes[tet.nodeIds[2]].position - grid.translation) / grid.scale;
  glm::vec3 x4 =
      (nodes[tet.nodeIds[3]].position - grid.translation) / grid.scale;

  glm::vec3 min(std::numeric_limits<float>::max());
  glm::vec3 max(std::numeric_limits<float>::lowest());

  SpatialHashGridCellRange range{};
  range.minX = static_cast<int64_t>(
      glm::floor(glm::min(glm::min(glm::min(x1.x, x2.x), x3.x), x4.x)));
  range.minY = static_cast<int64_t>(
      glm::floor(glm::min(glm::min(glm::min(x1.y, x2.y), x3.y), x4.y)));
  range.minZ = static_cast<int64_t>(
      glm::floor(glm::min(glm::min(glm::min(x1.z, x2.z), x3.z), x4.z)));

  range.lengthX = glm::max(
      static_cast<uint32_t>(
          glm::max(glm::max(glm::max(x1.x, x2.x), x3.x), x4.x) - range.minX),
      1u);
  range.lengthY = glm::max(
      static_cast<uint32_t>(
          glm::max(glm::max(glm::max(x1.y, x2.y), x3.y), x4.y) - range.minY),
      1u);
  range.lengthZ = glm::max(
      static_cast<uint32_t>(
          glm::max(glm::max(glm::max(x1.z, x2.z), x3.z), x4.z) - range.minZ),
      1u);

  return range;
}
} // namespace Pies