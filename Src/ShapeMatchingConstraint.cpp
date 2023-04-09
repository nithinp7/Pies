#include "ShapeMatchingConstraint.h"

#include <Eigen/Geometry>

namespace Pies {
ShapeMatchingConstraint::ShapeMatchingConstraint(
    const std::vector<uint32_t>& indices,
    const std::vector<glm::vec3>& materialCoordinates,
    float w)
    : _nodeIndices(indices),
      _currentPositions(3, materialCoordinates.size()),
      _materialCoordinates(3, materialCoordinates.size()),
      _projectedPositions(3, materialCoordinates.size()),
      _w(w) {
  size_t vertexCount = materialCoordinates.size();

  glm::vec3 centerOfMass(0.0f);
  float weight = 1.0f / static_cast<float>(vertexCount);
  for (const glm::vec3& coord : materialCoordinates) {
    centerOfMass += weight * coord;
  }

  for (size_t i = 0; i < vertexCount; ++i) {
    const glm::vec3& coordinate = materialCoordinates[i];
    this->_materialCoordinates.coeffRef(0, i) = coordinate.x - centerOfMass.x;
    this->_materialCoordinates.coeffRef(1, i) = coordinate.y - centerOfMass.y;
    this->_materialCoordinates.coeffRef(2, i) = coordinate.z - centerOfMass.z;
  }
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

void ShapeMatchingConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {
  glm::vec3 centerOfMass(0.0f);
  float weight = 1.0f / static_cast<float>(this->_nodeIndices.size());
  for (uint32_t nodeId : this->_nodeIndices) {
    centerOfMass += weight * nodes[nodeId].position;
  }

  for (uint32_t i = 0; i < this->_nodeIndices.size(); ++i) {
    uint32_t nodeId = this->_nodeIndices[i];
    const glm::vec3& currentPosition = nodes[nodeId].position;
    this->_currentPositions.coeffRef(0, i) = currentPosition.x - centerOfMass.x;
    this->_currentPositions.coeffRef(1, i) = currentPosition.y - centerOfMass.y;
    this->_currentPositions.coeffRef(2, i) = currentPosition.z - centerOfMass.z;
  }

  Eigen::Matrix4f T = Eigen::umeyama(
      this->_materialCoordinates,
      this->_currentPositions,
      false);
  Eigen::Matrix3f R = T.block(0, 0, 3, 3);
  Eigen::Vector3f t(centerOfMass.x, centerOfMass.y, centerOfMass.z);//T.block(0, 3, 3, 1);

  this->_projectedPositions.noalias() =
      (R * this->_materialCoordinates).colwise() + t;
}
} // namespace Pies