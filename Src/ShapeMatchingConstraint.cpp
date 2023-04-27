#include "ShapeMatchingConstraint.h"

#include <Eigen/Geometry>

namespace Pies {
ShapeMatchingConstraint::ShapeMatchingConstraint(
    const std::vector<Node>& nodes,
    const std::vector<uint32_t>& indices,
    const std::vector<glm::vec3>& materialCoordinates,
    float w)
    : _nodeIndices(indices),
      _materialCoords(3, materialCoordinates.size()),
      _projectedPositions(3, materialCoordinates.size()),
      _currentRotation(Eigen::Matrix3d::Identity()),
      _w(w) {
  size_t vertexCount = materialCoordinates.size();

  glm::vec3 centerOfMass(0.0f);
  float weight = 1.0f / static_cast<float>(vertexCount);
  for (const glm::vec3& coord : materialCoordinates) {
    centerOfMass += weight * coord;
  }

  Eigen::Matrix3d Q = Eigen::Matrix3d::Zero();
  for (size_t i = 0; i < vertexCount; ++i) {
    glm::vec3 matCoord = materialCoordinates[i] - centerOfMass;
    this->_materialCoords.coeffRef(0, i) = matCoord.x;
    this->_materialCoords.coeffRef(1, i) = matCoord.y;
    this->_materialCoords.coeffRef(2, i) = matCoord.z;

    glm::mat3 dQ =
        glm::outerProduct(matCoord, matCoord) / nodes[indices[i]].invMass;

    Q.coeffRef(0, 0) += dQ[0][0];
    Q.coeffRef(0, 1) += dQ[1][0];
    Q.coeffRef(0, 2) += dQ[2][0];

    Q.coeffRef(1, 0) += dQ[0][1];
    Q.coeffRef(1, 1) += dQ[1][1];
    Q.coeffRef(1, 2) += dQ[2][1];

    Q.coeffRef(2, 0) += dQ[0][2];
    Q.coeffRef(2, 1) += dQ[1][2];
    Q.coeffRef(2, 2) += dQ[2][2];
  }

  this->_Qinv = Q.inverse();
}

void ShapeMatchingConstraint::setupGlobalStiffnessMatrix(
    Eigen::SparseMatrix<float>& systemMatrix) const {
  // Assumes A and B are Identity.
  for (uint32_t nodeId : this->_nodeIndices) {
    systemMatrix.coeffRef(nodeId, nodeId) += this->_w;
  }
}

void ShapeMatchingConstraint::setupGlobalForceVector(
    Eigen::MatrixXf& forceVector) const {
  // Assumes A and B are Identity.
  for (uint32_t i = 0; i < this->_nodeIndices.size(); ++i) {
    uint32_t nodeId = this->_nodeIndices[i];

    forceVector.coeffRef(nodeId, 0) +=
        this->_w * this->_projectedPositions.coeff(0, i);
    forceVector.coeffRef(nodeId, 1) +=
        this->_w * this->_projectedPositions.coeff(1, i);
    forceVector.coeffRef(nodeId, 2) +=
        this->_w * this->_projectedPositions.coeff(2, i);
  }
}

// This routine is from the paper:
// A Robust Method to Extract the Rotational Part of Deformations, Muller et al.
void extractRotation(
    const Eigen::Matrix3d& A,
    Eigen::Quaterniond& q,
    const unsigned int maxIter) {
  for (unsigned int iter = 0; iter < maxIter; iter++) {
    Eigen::Matrix3d R = q.matrix();
    Eigen::Vector3d omega =
        (R.col(0).cross(A.col(0)) + R.col(1).cross(A.col(1)) +
         R.col(2).cross(A.col(2))) *
        (1.0 / fabs(
                   R.col(0).dot(A.col(0)) + R.col(1).dot(A.col(1)) +
                   R.col(2).dot(A.col(2))) +
         1.0e-9);
    double w = omega.norm();
    if (w < 1.0e-9)
      break;
    q = Eigen::Quaterniond(Eigen::AngleAxisd(w, (1.0 / w) * omega)) * q;
    q.normalize();
  }
}

void ShapeMatchingConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {
  glm::vec3 centerOfMass(0.0f);
  float weight = 1.0f / static_cast<float>(this->_nodeIndices.size());
  for (uint32_t nodeId : this->_nodeIndices) {
    centerOfMass += weight * nodes[nodeId].position;
  }

  Eigen::Matrix3d P = Eigen::Matrix3d::Zero();
  for (size_t i = 0; i < this->_nodeIndices.size(); ++i) {
    const Node& node = nodes[this->_nodeIndices[i]];
    glm::vec3 localCoord = node.position - centerOfMass;
    P += Eigen::Vector3d(localCoord.x, localCoord.y, localCoord.z) *
         Eigen::Vector3d(this->_materialCoords.col(i)).transpose() /
         node.invMass;
  }

  // Deformation gradient
  Eigen::Matrix3d F = P * this->_Qinv;

  // Want to extract rotation from
  extractRotation(F, this->_currentRotation, 100);
  Eigen::Matrix3d R = this->_currentRotation.toRotationMatrix();
  Eigen::Vector3d T(centerOfMass.x, centerOfMass.y, centerOfMass.z);

  this->_projectedPositions.noalias() =
      (R * this->_materialCoords).colwise() + T;
}
} // namespace Pies