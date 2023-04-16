#include "CollisionConstraint.h"

#include <algorithm>
#include <optional>

namespace Pies {
CollisionConstraint::CollisionConstraint(const Node& a, const Node& b)
    : nodeIds{a.id, b.id}, n(0.0f) {}

void CollisionConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {
  const Node& nodeA = nodes[this->nodeIds[0]];
  const Node& nodeB = nodes[this->nodeIds[1]];

  this->projectedPositions[0] = nodeA.position;
  this->projectedPositions[1] = nodeB.position;

  glm::vec3 diff = nodeB.position - nodeA.position;
  float distSq = glm::dot(diff, diff);
  float r = nodeA.radius + nodeB.radius;
  float rSq = r * r;
  if (distSq >= rSq) {
    return;
  }

  float dist = sqrt(distSq);
  float dispLength = r - dist;

  glm::vec3 disp;
  if (dist > 0.00001f) {
    disp = dispLength * diff / dist;
  } else {
    disp = glm::vec3(dispLength, 0.0f, 0.0f);
  }

  float wSum = nodeA.invMass + nodeB.invMass;

  this->projectedPositions[0] -= disp * nodeA.invMass / wSum;
  this->projectedPositions[1] += disp * nodeB.invMass / wSum;
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

PointTriangleCollisionConstraint::PointTriangleCollisionConstraint(
    const Node& a,
    const Node& b,
    const Node& c,
    const Node& d)
    : nodeIds{a.id, b.id, c.id, d.id}, n(0.0f) {}

void PointTriangleCollisionConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {

  const Node& nodeA = nodes[nodeIds[0]];
  const Node& nodeB = nodes[nodeIds[1]];
  const Node& nodeC = nodes[nodeIds[2]];
  const Node& nodeD = nodes[nodeIds[3]];

  this->projectedPositions[0] = nodeA.position;
  this->projectedPositions[1] = nodeB.position;
  this->projectedPositions[2] = nodeC.position;
  this->projectedPositions[3] = nodeD.position;

  glm::vec3 c = nodeC.position - nodeB.position;
  glm::vec3 d = nodeD.position - nodeB.position;
  glm::vec3 p = nodeA.position - nodeB.position;

  glm::vec3 n = glm::normalize(glm::cross(
      nodeC.position - nodeB.position,
      nodeD.position - nodeB.position));

  float nDotP = glm::dot(n, p);
  if (nDotP < thickness) {
    glm::vec3 disp = (thickness - nDotP) * n;

    float wTriSum = nodeB.invMass + nodeC.invMass + nodeD.invMass;
    float wSum = nodeA.invMass + wTriSum;

    this->projectedPositions[0] += disp * nodeA.invMass / wSum;
    this->projectedPositions[1] -= disp * wTriSum / wSum;
    this->projectedPositions[2] -= disp * wTriSum / wSum;
    this->projectedPositions[3] -= disp * wTriSum / wSum;
  }
}

void PointTriangleCollisionConstraint::stabilizeCollisions(
    std::vector<Node>& nodes) {
  Node& nodeA = nodes[nodeIds[0]];
  Node& nodeB = nodes[nodeIds[1]];
  Node& nodeC = nodes[nodeIds[2]];
  Node& nodeD = nodes[nodeIds[3]];

  // If the point is behind the triangle but within a desired thickness
  // just push it out.
  glm::vec3 c = nodeC.position - nodeB.position;
  glm::vec3 d = nodeD.position - nodeB.position;
  glm::vec3 p = nodeA.position - nodeB.position;

  glm::vec3 n = glm::normalize(glm::cross(
      nodeC.position - nodeB.position,
      nodeD.position - nodeB.position));

  float nDotP = glm::dot(n, p);
  if (nDotP < thickness) {
    glm::vec3 disp = (thickness - nDotP) * n;

    float wTriSum = nodeB.invMass + nodeC.invMass + nodeD.invMass;
    float wSum = nodeA.invMass + wTriSum;

    nodeA.position += disp * nodeA.invMass / wSum;
    nodeB.position -= disp * wTriSum / wSum;
    nodeC.position -= disp * wTriSum / wSum;
    nodeD.position -= disp * wTriSum / wSum;
  }
}

void PointTriangleCollisionConstraint::setupCollisionMatrix(
    Eigen::SparseMatrix<float>& systemMatrix) const {
  systemMatrix.coeffRef(nodeIds[0], nodeIds[0]) += this->w;
  systemMatrix.coeffRef(nodeIds[1], nodeIds[1]) += this->w;
  systemMatrix.coeffRef(nodeIds[2], nodeIds[2]) += this->w;
  systemMatrix.coeffRef(nodeIds[3], nodeIds[3]) += this->w;
}

void PointTriangleCollisionConstraint::setupGlobalForceVector(
    Eigen::MatrixXf& forceVector) const {
  // TODO: How do we make this "unilateral" - is that covered in the projection
  // step?

  forceVector.coeffRef(nodeIds[0], 0) += this->w * projectedPositions[0].x;
  forceVector.coeffRef(nodeIds[0], 1) += this->w * projectedPositions[0].y;
  forceVector.coeffRef(nodeIds[0], 2) += this->w * projectedPositions[0].z;

  forceVector.coeffRef(nodeIds[1], 0) += this->w * projectedPositions[1].x;
  forceVector.coeffRef(nodeIds[1], 1) += this->w * projectedPositions[1].y;
  forceVector.coeffRef(nodeIds[1], 2) += this->w * projectedPositions[1].z;

  forceVector.coeffRef(nodeIds[2], 0) += this->w * projectedPositions[2].x;
  forceVector.coeffRef(nodeIds[2], 1) += this->w * projectedPositions[2].y;
  forceVector.coeffRef(nodeIds[2], 2) += this->w * projectedPositions[2].z;

  forceVector.coeffRef(nodeIds[3], 0) += this->w * projectedPositions[3].x;
  forceVector.coeffRef(nodeIds[3], 1) += this->w * projectedPositions[3].y;
  forceVector.coeffRef(nodeIds[3], 2) += this->w * projectedPositions[3].z;
}

StaticCollisionConstraint::StaticCollisionConstraint(
    const Node& node,
    const glm::vec3& projectedPosition_)
    : nodeId(node.id), projectedPosition(projectedPosition_), n(0.0f) {
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