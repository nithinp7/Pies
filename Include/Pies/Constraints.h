#pragma once

#include "Node.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

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

  // TODO: This should be defined by the templated constraint implementation,
  // individual constraints do not need to differ wrt matrix A and B.

  // TODO: This should be defined by the templated constraint implementation,
  // individual constraints do not need to differ wrt matrix A and B.

  // Precomputed values defining the energy potential for the constraint
  // See the Projective Dynamics paper for more information about A and B.
  Eigen::Matrix<float, NodeCount, NodeCount> _AtA;
  Eigen::Matrix<float, NodeCount, NodeCount> _AtB;
  Eigen::Matrix<float, NodeCount, 3> _wAtBp;

  // A list of nodes involved in this constraint.
  std::array<uint32_t, NodeCount> _nodeIds;

  // The auxiliary variable containing projected node configurations.
  std::array<glm::vec3, NodeCount> _projectedConfig;

  TProjection _projection;

public:
  Constraint(
      uint32_t id,
      float w,
      const Eigen::Matrix<float, NodeCount, NodeCount>& A,
      const Eigen::Matrix<float, NodeCount, NodeCount>& B,
      const TProjection& projection,
      const std::array<uint32_t, NodeCount>& nodeIds)
      : _id(id),
        _w(w),
        _AtA(A.transpose() * A),
        _AtB(A.transpose() * B),
        _projection(std::move(projection)),
        _nodeIds(nodeIds) {}

  /**
   * @brief Initialize the constraint within the global stiffness matrix. Only
   * needs to be called during initialization or if node connectivity changes.
   *
   * @param stiffnessMatrix The global matrix describing the relationship
   * between nodal displacements and nodal forces.
   */
  void
  setupGlobalStiffnessMatrix(Eigen::SparseMatrix<float>& systemMatrix) const {
    // Note: We only fill in the lower triangular part of the matrix.
    for (uint32_t i = 0; i < NodeCount; ++i) {
      uint32_t nodeId_i = this->_nodeIds[i];
      for (uint32_t j = 0; j < NodeCount; ++j) {
        uint32_t nodeId_j = this->_nodeIds[j];
        systemMatrix.coeffRef(nodeId_i, nodeId_j) +=
            this->_w * this->_AtA.coeff(i, j);
      }
    }
  }

  /**
   * @brief Add nodal forces due to this constraint to the global force vector.
   * This is called during every solver iteration.
   *
   * @param forceVector The global vector of forces on each axis of each node.
   */
  void setupGlobalForceVector(Eigen::MatrixXf& forceVector) const {
    for (uint32_t i = 0; i < NodeCount; ++i) {
      uint32_t nodeId_i = this->_nodeIds[i];
      forceVector.coeffRef(nodeId_i, 0) += this->_wAtBp.coeff(i, 0);
      forceVector.coeffRef(nodeId_i, 1) += this->_wAtBp.coeff(i, 1);
      forceVector.coeffRef(nodeId_i, 2) += this->_wAtBp.coeff(i, 2);
    }
  }

  /**
   * @brief Projects the current node configuration onto the constraint
   * manifold. Stores the projected configuration in the auxiliary variable.
   * Uses the TProjection template parameter to do the projection.
   */
  void projectToAuxiliaryVariable(const std::vector<Node>& nodes) {
    this->_projection(nodes, this->_nodeIds, this->_projectedConfig);    
    
    // Set up projected nodes as eigen matrix
    Eigen::Matrix<float, NodeCount, 3> p;
    for (uint32_t i = 0; i < NodeCount; ++i) {
      p.coeffRef(i, 0) = this->_projectedConfig[i].x;
      p.coeffRef(i, 1) = this->_projectedConfig[i].y;
      p.coeffRef(i, 2) = this->_projectedConfig[i].z;
    }

    this->_wAtBp = this->_w * this->_AtB * p;
  }

  /**
   * @brief Projects the current node configuration onto the constraint
   * manifold and directly updates the node positions. This is only used in the
   * PBD solver. Uses the TProjection template parameter to do the projection.
   */
  void projectNodePositions(std::vector<Node>& nodes) {
    std::array<glm::vec3, NodeCount> fixedPositions;
    this->_projection(nodes, this->_nodeIds, fixedPositions);

    for (uint32_t i = 0; i < NodeCount; ++i) {
      Node& node = nodes[this->_nodeIds[i]];
      node.position += this->_w * (fixedPositions[i] - node.position);
    }
  }

  uint32_t getNodeId(uint32_t nodeIndex) const {
    if (nodeIndex >= NodeCount) {
      throw std::runtime_error(
          "Invalid nodeIndex given to Constraint::getNode.");
    }

    return this->_nodeIds[nodeIndex];
  }

  void setWeight(float w) { this->_w = w; }

  TProjection& getProjection() {
    return this->_projection;
  }
};

struct DistanceConstraintProjection {
  float targetDistance;

  void operator()(
      const std::vector<Node>& nodes,
      const std::array<uint32_t, 2>& nodeIds,
      std::array<glm::vec3, 2>& projected) const;
};
typedef Constraint<2, DistanceConstraintProjection> DistanceConstraint;
DistanceConstraint
createDistanceConstraint(uint32_t id, const Node& a, const Node& b, float w);

struct PositionConstraintProjection {
  glm::vec3 fixedPosition;
  glm::vec3 offset;

  void operator()(
      const std::vector<Node>& nodes,
      const std::array<uint32_t, 1>& nodeIds,
      std::array<glm::vec3, 1>& projected) const;
};
typedef Constraint<1, PositionConstraintProjection> PositionConstraint;
PositionConstraint createPositionConstraint(uint32_t id, const Node& node, float w);

struct TetrahedralConstraintProjection {
  glm::mat3 Q;
  glm::mat3 Qinv;

  float minStrain;
  float maxStrain;

  void operator()(
      const std::vector<Node>& nodes,
      const std::array<uint32_t, 4>& nodeIds,
      std::array<glm::vec3, 4>& projected) const;
};
typedef Constraint<4, TetrahedralConstraintProjection> TetrahedralConstraint;
TetrahedralConstraint createTetrahedralConstraint(
    uint32_t id,
    float w,
    const Node& a,
    const Node& b,
    const Node& c,
    const Node& d,
    float minStrain = 0.8f,
    float maxStrain = 1.0f);

struct VolumeConstraintProjection {
  glm::mat3 Qinv;
  float minOmega;
  float maxOmega;

  void operator()(
      const std::vector<Node>& nodes,
      const std::array<uint32_t, 4>& nodeIds,
      std::array<glm::vec3, 4>& projected) const;
};
typedef Constraint<4, VolumeConstraintProjection> VolumeConstraint;
VolumeConstraint createVolumeConstraint(
    uint32_t id,
    float w,
    const Node& a,
    const Node& b,
    const Node& c,
    const Node& d,
    float compression = 1.0f,
    float stretching = 1.0f);

struct BendConstraintProjection {
  float initialAngle;

  void operator()(
      const std::vector<Node>& nodes,
      const std::array<uint32_t, 4>& nodeIds,
      std::array<glm::vec3, 4>& projected) const;
};
typedef Constraint<4, BendConstraintProjection> BendConstraint;
BendConstraint createBendConstraint(
    uint32_t id,
    float w,
    const Node& a,
    const Node& b,
    const Node& c,
    const Node& d);

} // namespace Pies
