#pragma once

#include "Node.h"

#include <array>
#include <cstdint>
#include <stdexcept>

namespace Pies {

/**
 * @param NodeCount The number of nodes involved in this constraint.
 * @param TProjection Projects the current node configuration onto the
 * constraint manifold.
 */
template <unsigned int NodeCount, typename TProjection> class Constraint {
protected:
  // This ID represents the row index of this constraint within the global
  // stiffness matrix.
  uint32_t _id;

  // A constant weight for this constraint.
  float _w = 1.0f;

  // The constant matrix A.
  // Matrix<NodeCount, NodeCount> _A;

  // The constant matrix B.
  // Matrix<NodeCount, NodeCount> _B;

  // The auxiliary variable containing projected node configurations.
  std::array<glm::vec3, NodeCount> _projectedConfig;

  TProjection _projection;

public:
  // TODO: rethink constraints owning Node pointers... global node list may
  // reallocate between frames. _nodes is currently public to fixup the node
  // list on reallocation - it is a dirty hack and should be removed soon.

  // A list of nodes involved in this constraint.
  std::array<Node*, NodeCount> _nodes;

  Constraint(
      uint32_t id,
      float w,
      const TProjection& projection,
      const std::array<Node*, NodeCount>& nodes)
      : _id(id), _w(w), _projection(std::move(projection)), _nodes(nodes) {}

  /**
   * @brief Initialize the constraint within the global stiffness matrix. Only
   * needs to be called during initialization or if node connectivity changes.
   *
   * @param stiffnessMatrix The global matrix describing the relationship
   * between nodal displacements and nodal forces.
   */
  // void setupGlobalStiffnessMatrix(Matrix& systemMatrix) const {

  // }

  /**
   * @brief Add nodal forces due to this constraint to the global force vector.
   * This is called during every solver iteration.
   *
   * @param forceVector The global vector of forces on each axis of each node.
   */
  // void setupGlobalForceVector(Vector& forceVector) const {

  // }

  /**
   * @brief Projects the current node configuration onto the constraint
   * manifold. Stores the projected configuration in the auxiliary variable.
   * Uses the TProjection template parameter to do the projection.
   */
  void projectToAuxiliaryVariable() {
    this->_projection(this->_nodes, this->_projectedConfig);
  }

  /**
   * @brief Projects the current node configuration onto the constraint
   * manifold and directly updates the node positions. This is only used in the
   * PBD solver. Uses the TProjection template parameter to do the projection.
   */
  void projectNodePositions() {
    std::array<glm::vec3, NodeCount> fixedPositions;
    this->_projection(this->_nodes, fixedPositions);

    for (uint32_t i = 0; i < NodeCount; ++i) {
      this->_nodes[i]->position +=
          this->_w * (fixedPositions[i] - this->_nodes[i]->position);
    }
  }

  const Node& getNode(uint32_t nodeIndex) const {
    if (nodeIndex >= NodeCount) {
      throw std::runtime_error(
          "Invalid nodeIndex given to Constraint::getNode.");
    }

    return *this->_nodes[nodeIndex];
  }

  void setWeight(float w) { this->_w = w; }
};

struct DistanceConstraintProjection {
  float targetDistance;

  void operator()(
      const std::array<Node*, 2>& nodes,
      std::array<glm::vec3, 2>& projected) const;
};
typedef Constraint<2, DistanceConstraintProjection> DistanceConstraint;
DistanceConstraint createDistanceConstraint(uint32_t id, Node* a, Node* b);

struct PositionConstraintProjection {
  glm::vec3 fixedPosition;

  void operator()(
      const std::array<Node*, 1>& nodes,
      std::array<glm::vec3, 1>& projected) const;
};
typedef Constraint<1, PositionConstraintProjection> PositionConstraint;
PositionConstraint createPositionConstraint(uint32_t id, Node* node);

// struct CollisionConstraintProjection {
//   glm::vec3 intersection;
//   glm::vec3 normal;

//   void operator()(
//       const std::array<Node*, 2>& nodes,
//       std::array<glm::vec3, 2>& projected) const;
// };
// typedef Constraint<2, CollisionConstraintProjection> CollisionConstraint;
// CollisionConstraint createCollisionConstraint(uint32_t id, Node* a, Node* b);
} // namespace Pies