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

  for (size_t i = 0; i < this->_materialCoordinates.size(); ++i) {
    const glm::vec3& coordinate = materialCoordinates[i];
    this->_materialCoordinates.coeffRef(0, i) = coordinate.x;
    this->_materialCoordinates.coeffRef(1, i) = coordinate.y;
    this->_materialCoordinates.coeffRef(2, i) = coordinate.z;
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
        this->_w * this->_projectedPositions.coeff(i, 0);
    forceVector.coeffRef(nodeId, 1) +=
        this->_w * this->_projectedPositions.coeff(i, 1);
    forceVector.coeffRef(nodeId, 2) +=
        this->_w * this->_projectedPositions.coeff(i, 2);
  }
}

void ShapeMatchingConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {
  for (uint32_t i = 0; i < this->_nodeIndices.size(); ++i) {
    uint32_t nodeId = this->_nodeIndices[i];
    const glm::vec3& currentPosition = nodes[nodeId].position;
    this->_currentPositions.coeffRef(0, i) = currentPosition.x;
    this->_currentPositions.coeffRef(1, i) = currentPosition.y;
    this->_currentPositions.coeffRef(2, i) = currentPosition.z;
  }

  Eigen::MatrixXf T = Eigen::umeyama(
      this->_materialCoordinates,
      this->_currentPositions,
      false);
  Eigen::Matrix3f R = T.block(0, 0, 3, 3);
  Eigen::Vector3f t = T.block(0, 3, 3, 1);

  // TODO: Check that this is doing column-wise addition of translation
  // as intended.
  this->_projectedPositions = R * this->_materialCoordinates + t;
}
} // namespace Pies