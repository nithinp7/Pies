#pragma once

#include "Node.h"
#include "eig3.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
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

  TProjection& getProjection() { return this->_projection; }
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
PositionConstraint
createPositionConstraint(uint32_t id, const Node& node, float w);

struct TetrahedralConstraintProjection {
  glm::mat3 Q;
  glm::mat3 Qinv;

  float minStrain;
  float maxStrain;

  float minOmega;
  float maxOmega;

  inline void operator()(
      const std::vector<Node>& nodes,
      const std::array<uint32_t, 4>& nodeIds,
      std::array<glm::vec3, 4>& projected) const {
    const Node& x1 = nodes[nodeIds[0]];
    const Node& x2 = nodes[nodeIds[1]];
    const Node& x3 = nodes[nodeIds[2]];
    const Node& x4 = nodes[nodeIds[3]];

    glm::mat3 P(
        x2.position - x1.position,
        x3.position - x1.position,
        x4.position - x1.position);

    // Deformation gradient
    glm::mat3 F = P * this->Qinv;

    Eigen::Matrix3d F_;
    F_ << F[0][0], F[0][1], F[0][2], F[1][0], F[1][1], F[1][2], F[2][0],
        F[2][1], F[2][2];

    double A[3][3];
    A[0][0] = glm::dot(F[0], F[0]);
    A[1][1] = glm::dot(F[1], F[1]);
    A[2][2] = glm::dot(F[2], F[2]);

    A[0][1] = A[1][0] = glm::dot(F[0], F[1]);
    A[0][2] = A[2][0] = glm::dot(F[0], F[2]);
    A[1][2] = A[2][1] = glm::dot(F[1], F[2]);

    double d_[3];
    double V_[3][3];
    eigen_decomposition(A, V_, d_);

    glm::mat3 V(
        static_cast<float>(V_[0][0]),
        static_cast<float>(V_[1][0]),
        static_cast<float>(V_[2][0]),
        static_cast<float>(V_[0][1]),
        static_cast<float>(V_[1][1]),
        static_cast<float>(V_[2][1]),
        static_cast<float>(V_[0][2]),
        static_cast<float>(V_[1][2]),
        static_cast<float>(V_[2][2]));
    
    glm::vec3 sigma(
        static_cast<float>(sqrt(d_[0])),
        static_cast<float>(sqrt(d_[1])),
        static_cast<float>(sqrt(d_[2])));

    const uint32_t COMP_D_ITERS = 10;
    glm::vec3 D(0.0f);
    for (uint32_t i = 0; i < COMP_D_ITERS; ++i) {
      glm::vec3 sigmaPlusD = sigma + D;
      float product = sigmaPlusD.x * sigmaPlusD.y * sigmaPlusD.z;
      float omega = glm::clamp(product, this->minOmega, this->maxOmega);
      float C = product - omega;
      glm::vec3 gradC(
          sigmaPlusD.y * sigmaPlusD.z,
          sigmaPlusD.x * sigmaPlusD.z,
          sigmaPlusD.x * sigmaPlusD.y);
      D = (glm::dot(gradC, D) - C) * gradC / glm::dot(gradC, gradC);
    }

    sigma += D;
    sigma = glm::clamp(sigma, this->minStrain, this->maxStrain);

    if (glm::determinant(F) < 0.0f) {
      sigma[2] *= -1.0f;
    }

    glm::vec3 d = sigma * sigma;

    // TODO: This is not right, U needs to be computed and all sorts of corner cases need
    // to be handled to use eigen decomp in this context - see:
    // Invertible Finite Elements For Robust Simulation of Large Deformation
    // Irving et. al, 2004

    // The "fixed" deformation gradient
    glm::mat3 Fhat = glm::mat3(V[0] * d[0], V[1] * d[1], V[2] * d[2]) * glm::transpose(V);
    glm::mat3 P1 = glm::transpose(Fhat);

    projected[0] = glm::vec3(0.0f);
    projected[1] = P1[0];
    projected[2] = P1[1];
    projected[3] = P1[2];
  }
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
    float maxStrain = 1.0f,
    float compression = 1.0f,
    float stretching = 1.0f);

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
