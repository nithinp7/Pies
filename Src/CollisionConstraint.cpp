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
    const Node& d,
    float thickness_)
    : nodeIds{a.id, b.id, c.id, d.id}, n(0.0f), thickness(thickness_) {
  Eigen::Matrix4f A = Eigen::Matrix4f::Zero();
  A.coeffRef(1, 0) = -1.0f;
  A.coeffRef(2, 0) = -1.0f;
  A.coeffRef(3, 0) = -1.0f;

  A.coeffRef(1, 1) = 1.0f;
  A.coeffRef(2, 2) = 1.0f;
  A.coeffRef(3, 3) = 1.0f;

  this->AtA = A.transpose() * A;
}

void PointTriangleCollisionConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {
  colliding = false;

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

    // TODO: Use relative coords

    this->projectedPositions[0] += disp;
    // this->projectedPositions[1] = glm::vec3(0.0f);
    // this->projectedPositions[2] = glm::vec3(0.0f);
    // this->projectedPositions[3] = glm::vec3(0.0f);

    colliding = true;
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

    // TODO: Use relative coords
    nodeA.position += disp * nodeA.invMass / wSum;
    nodeB.position -= disp * wTriSum / wSum;
    nodeC.position -= disp * wTriSum / wSum;
    nodeD.position -= disp * wTriSum / wSum;

    // This prevents spuriously adding velocity to the system
    nodeA.prevPosition += disp * nodeA.invMass / wSum;
    nodeB.prevPosition -= disp * wTriSum / wSum;
    nodeC.prevPosition -= disp * wTriSum / wSum;
    nodeD.prevPosition -= disp * wTriSum / wSum;
  }
}

void PointTriangleCollisionConstraint::setupCollisionMatrix(
    Eigen::SparseMatrix<float>& systemMatrix) const {
  for (uint32_t i = 0; i < 4; ++i) {
    uint32_t nodeId_i = this->nodeIds[i];
    for (uint32_t j = 0; j < 4; ++j) {
      uint32_t nodeId_j = this->nodeIds[j];
      systemMatrix.coeffRef(nodeId_i, nodeId_j) +=
          this->w * this->AtA.coeff(i, j);
    }
  }
}

void PointTriangleCollisionConstraint::setupGlobalForceVector(
    Eigen::MatrixXf& forceVector) const {
  // Set up projected nodes as eigen matrix
  Eigen::Matrix<float, 4, 3> p;
  for (uint32_t i = 0; i < 4; ++i) {
    p.coeffRef(i, 0) = this->projectedPositions[i].x;
    p.coeffRef(i, 1) = this->projectedPositions[i].y;
    p.coeffRef(i, 2) = this->projectedPositions[i].z;
  }

  // Reminder: A == B in this case
  Eigen::Matrix<float, 4, 3> AtBp = this->AtA * p;
  for (uint32_t i = 0; i < 4; ++i) {
    uint32_t nodeId_i = this->nodeIds[i];
    forceVector.coeffRef(nodeId_i, 0) += this->w * AtBp.coeff(i, 0);
    forceVector.coeffRef(nodeId_i, 1) += this->w * AtBp.coeff(i, 1);
    forceVector.coeffRef(nodeId_i, 2) += this->w * AtBp.coeff(i, 2);
  }
}

EdgeCollisionConstraint::EdgeCollisionConstraint(
    const Node& a,
    const Node& b,
    const Node& c,
    const Node& d)
    : nodeIds{a.id, b.id, c.id, d.id} {
  Eigen::Matrix4f A = Eigen::Matrix4f::Zero();
  A.coeffRef(1, 0) = -1.0f;
  A.coeffRef(2, 0) = -1.0f;
  A.coeffRef(3, 0) = -1.0f;

  A.coeffRef(1, 1) = 1.0f;
  A.coeffRef(2, 2) = 1.0f;
  A.coeffRef(3, 3) = 1.0f;

  this->AtA = A.transpose() * A;

  glm::vec3 ab = b.position - a.position;
  glm::vec3 ac = c.position - a.position;
  glm::vec3 ad = d.position - a.position;

  glm::vec3 n0 = glm::normalize(glm::cross(ab, ad - ac));

  float n0DotL1 = glm::max(glm::dot(n0, ab), 0.0f);
  float n0DotL2 = glm::min(glm::dot(n0, ac), glm::dot(n0, ad));

  orientation = n0DotL1 < n0DotL2 ? 1.0f : -1.0f;
}

void EdgeCollisionConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {
  const Node& nodeA = nodes[nodeIds[0]];
  const Node& nodeB = nodes[nodeIds[1]];
  const Node& nodeC = nodes[nodeIds[2]];
  const Node& nodeD = nodes[nodeIds[3]];

  this->projectedPositions[0] = nodeA.position;
  this->projectedPositions[1] = nodeB.position;
  this->projectedPositions[2] = nodeC.position;
  this->projectedPositions[3] = nodeD.position;

  glm::vec3 ab = nodeB.position - nodeA.position;
  glm::vec3 ac = nodeC.position - nodeA.position;
  glm::vec3 ad = nodeD.position - nodeA.position;

  glm::vec3 cd = nodeD.position - nodeC.position;

  float abMagSq = glm::dot(ab, ab);
  float cdMagSq = glm::dot(cd, cd);
  float abDotCd = glm::dot(ab, cd);

  float acDotAb = glm::dot(ac, ab);
  float acDotCd = glm::dot(ac, cd);

  float det = abMagSq * -cdMagSq + abDotCd * abDotCd;
  float u = 0.0f;
  float v = 0.0f;
  if (det != 0.0f) {
    det = 1.0f / det;
    float u = (acDotAb * -cdMagSq + abDotCd * acDotCd) * det;
    float v = (abMagSq * acDotCd - acDotAb * abDotCd) * det;
  } else {
    float u0 = glm::dot(nodeA.position, ab);
    float u1 = glm::dot(nodeB.position, ab);
    float v0 = glm::dot(nodeC.position, ab);
    float v1 = glm::dot(nodeD.position, ab);

    bool flip0 = false;
    bool flip1 = false;

    if (u0 > u1) {
      std::swap(u0, u1);
      flip0 = true;
    }

    if (v0 > v1) {
      std::swap(v0, v1);
      flip1 = true;
    }

    if (u0 >= v1) {
      u = flip0 ? 1.0f : 0.0f;
      v = flip1 ? 0.0f : 1.0f;
    } else if (v0 >= u1) {
      u = flip0 ? 0.0f : 1.0f;
      v = flip1 ? 1.0f : 0.0f;
    } else {
      float mid = (u0 > v0) ? (u0 + v1) * 0.5f : (v0 + u1) * 0.5f;
      u = (u0 == u1) ? 0.5f : (mid - u0) / (u1 - u0);
      v = (v0 == v1) ? 0.5f : (mid - v0) / (v1 - v0);
    }
  }

  u = glm::clamp(u, 0.0f, 1.0f);
  v = glm::clamp(v, 0.0f, 1.0f);

  glm::vec3 q0 = glm::mix(glm::vec3(0.0f), ab, u);
  glm::vec3 q1 = glm::mix(ac, ad, v);

  glm::vec3 n = q0 - q1;
  float dist = glm::length(n);
  n /= dist;

  if (dist < thickness) {
    glm::vec3 disp = -(thickness - dist) * n;

    float s = nodeA.invMass * (1.0f - u) * (1.0f - u) + nodeB.invMass * u * u +
              nodeC.invMass * (1.0f - v) * (1.0f - v) + nodeD.invMass * v * v;

    if (s == 0.0f) {
      return;
    }

    this->projectedPositions[0] += disp * nodeA.invMass * (1.0f - u) / s;
    this->projectedPositions[1] += disp * nodeB.invMass * u / s;
    this->projectedPositions[2] -= disp * nodeC.invMass * (1.0f - v) / s;
    this->projectedPositions[3] -= disp * nodeD.invMass * v / s;
  }
}

void EdgeCollisionConstraint::stabilizeCollisions(std::vector<Node>& nodes) {
  Node& nodeA = nodes[nodeIds[0]];
  Node& nodeB = nodes[nodeIds[1]];
  Node& nodeC = nodes[nodeIds[2]];
  Node& nodeD = nodes[nodeIds[3]];


  glm::vec3 ab = nodeB.position - nodeA.position;
  glm::vec3 ac = nodeC.position - nodeA.position;
  glm::vec3 ad = nodeD.position - nodeA.position;

  glm::vec3 cd = nodeD.position - nodeC.position;

  float abMagSq = glm::dot(ab, ab);
  float cdMagSq = glm::dot(cd, cd);
  float abDotCd = glm::dot(ab, cd);

  float acDotAb = glm::dot(ac, ab);
  float acDotCd = glm::dot(ac, cd);

  float det = abMagSq * -cdMagSq + abDotCd * abDotCd;
  float u = 0.0f;
  float v = 0.0f;
  if (det != 0.0f) {
    det = 1.0f / det;
    float u = (acDotAb * -cdMagSq + abDotCd * acDotCd) * det;
    float v = (abMagSq * acDotCd - acDotAb * abDotCd) * det;
  } else {
    float u0 = glm::dot(nodeA.position, ab);
    float u1 = glm::dot(nodeB.position, ab);
    float v0 = glm::dot(nodeC.position, ab);
    float v1 = glm::dot(nodeD.position, ab);

    bool flip0 = false;
    bool flip1 = false;

    if (u0 > u1) {
      std::swap(u0, u1);
      flip0 = true;
    }

    if (v0 > v1) {
      std::swap(v0, v1);
      flip1 = true;
    }

    if (u0 >= v1) {
      u = flip0 ? 1.0f : 0.0f;
      v = flip1 ? 0.0f : 1.0f;
    } else if (v0 >= u1) {
      u = flip0 ? 0.0f : 1.0f;
      v = flip1 ? 1.0f : 0.0f;
    } else {
      float mid = (u0 > v0) ? (u0 + v1) * 0.5f : (v0 + u1) * 0.5f;
      u = (u0 == u1) ? 0.5f : (mid - u0) / (u1 - u0);
      v = (v0 == v1) ? 0.5f : (mid - v0) / (v1 - v0);
    }
  }

  u = glm::clamp(u, 0.0f, 1.0f);
  v = glm::clamp(v, 0.0f, 1.0f);

  glm::vec3 q0 = glm::mix(glm::vec3(0.0f), ab, u);
  glm::vec3 q1 = glm::mix(ac, ad, v);

  glm::vec3 n = q0 - q1;
  float dist = glm::length(n);
  n /= dist;

  if (dist < thickness) {
    glm::vec3 disp = (thickness - dist) * n;

    float s = nodeA.invMass * (1.0f - u) * (1.0f - u) + nodeB.invMass * u * u +
              nodeC.invMass * (1.0f - v) * (1.0f - v) + nodeD.invMass * v * v;

    if (s == 0.0f) {
      return;
    }

    nodeA.position += disp * nodeA.invMass * (1.0f - u) / s;
    nodeB.position += disp * nodeB.invMass * u / s;
    nodeC.position -= disp * nodeC.invMass * (1.0f - v) / s;
    nodeD.position -= disp * nodeD.invMass * v / s;

    // This prevents spuriously adding velocity to the system
    nodeA.prevPosition += disp * nodeA.invMass * (1.0f - u) / s;
    nodeB.prevPosition += disp * nodeB.invMass * u / s;
    nodeC.prevPosition -= disp * nodeC.invMass * (1.0f - v) / s;
    nodeD.prevPosition -= disp * nodeD.invMass * v / s;
  }
}

void EdgeCollisionConstraint::setupCollisionMatrix(
    Eigen::SparseMatrix<float>& systemMatrix) const {
  for (uint32_t i = 0; i < 4; ++i) {
    uint32_t nodeId_i = this->nodeIds[i];
    for (uint32_t j = 0; j < 4; ++j) {
      uint32_t nodeId_j = this->nodeIds[j];
      systemMatrix.coeffRef(nodeId_i, nodeId_j) +=
          this->w * this->AtA.coeff(i, j);
    }
  }
}

void EdgeCollisionConstraint::setupGlobalForceVector(
    Eigen::MatrixXf& forceVector) const {
  // Set up projected nodes as eigen matrix
  Eigen::Matrix<float, 4, 3> p;
  for (uint32_t i = 0; i < 4; ++i) {
    p.coeffRef(i, 0) = this->projectedPositions[i].x;
    p.coeffRef(i, 1) = this->projectedPositions[i].y;
    p.coeffRef(i, 2) = this->projectedPositions[i].z;
  }

  // Reminder: A == B in this case
  Eigen::Matrix<float, 4, 3> AtBp = this->AtA * p;
  for (uint32_t i = 0; i < 4; ++i) {
    uint32_t nodeId_i = this->nodeIds[i];
    forceVector.coeffRef(nodeId_i, 0) += this->w * AtBp.coeff(i, 0);
    forceVector.coeffRef(nodeId_i, 1) += this->w * AtBp.coeff(i, 1);
    forceVector.coeffRef(nodeId_i, 2) += this->w * AtBp.coeff(i, 2);
  }
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