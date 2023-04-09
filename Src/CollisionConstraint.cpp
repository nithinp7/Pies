#include "CollisionConstraint.h"

#include <algorithm>

namespace Pies {
CollisionConstraint::CollisionConstraint(
    const Node& a,
    const Node& b,
    const glm::vec3& disp)
    : nodeIds{a.id, b.id}, 
      n(0.0f) {
        
  float wSum = a.invMass + b.invMass;

  this->projectedPositions[0] = a.position - disp * a.invMass / wSum;
  this->projectedPositions[1] = b.position + disp * b.invMass / wSum;

  float dispLength = glm::length(disp);
  if (dispLength > 0.00001f) {
    this->n = disp / dispLength;
  }
}

void CollisionConstraint::setupCollisionMatrix(
    Eigen::SparseMatrix<float>& systemMatrix) const {
  systemMatrix.coeffRef(nodeIds[0], nodeIds[0]) += this->w;
  systemMatrix.coeffRef(nodeIds[1], nodeIds[1]) += this->w;  
}

void CollisionConstraint::setupGlobalForceVector(
    Eigen::MatrixXf& forceVector,
    uint32_t threadId,
    uint32_t threadCount) const {
  // TODO: Do we need to declar no alias for projectedPositions or nodeIds??

  if (threadId == nodeIds[0] % threadCount) {
    forceVector.coeffRef(nodeIds[0], 0) += this->w * projectedPositions[0].x;
    forceVector.coeffRef(nodeIds[0], 1) += this->w * projectedPositions[0].y;
    forceVector.coeffRef(nodeIds[0], 2) += this->w * projectedPositions[0].z;
  }

  if (threadId == nodeIds[1] % threadCount) {
    forceVector.coeffRef(nodeIds[1], 0) += this->w * projectedPositions[1].x;
    forceVector.coeffRef(nodeIds[1], 1) += this->w * projectedPositions[1].y;
    forceVector.coeffRef(nodeIds[1], 2) += this->w * projectedPositions[1].z;
  }
}

StaticCollisionConstraint::StaticCollisionConstraint(
    const Node& node,
    const glm::vec3& projectedPosition_)
    : nodeId(node.id), 
      projectedPosition(projectedPosition_),
      n(0.0f) {
  glm::vec3 diff = projectedPosition - node.position;
  float dist = glm::length(diff);
  if (dist > 0.00001f) {
    n = diff / dist;
  }
}

void StaticCollisionConstraint::setupCollisionMatrix(
    Eigen::SparseMatrix<float>& systemMatrix) const {
  systemMatrix.coeffRef(nodeId, nodeId) += this->w;
}

void StaticCollisionConstraint::setupGlobalForceVector(
    Eigen::MatrixXf& forceVector) const {
  // TODO: Do we need to declar no alias for projectedPositions or nodeIds??
  forceVector.coeffRef(nodeId, 0) += this->w * projectedPosition.x;
  forceVector.coeffRef(nodeId, 1) += this->w * projectedPosition.y;
  forceVector.coeffRef(nodeId, 2) += this->w * projectedPosition.z;
}
} // namespace Pies