#include "CollisionConstraint.h"

#include <glm/common.hpp>

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
    float thickness_,
    const glm::vec3& barycentric,
    float toi_)
    : toi(toi_),
      nodeIds{a.id, b.id, c.id, d.id},
      n(0.0f),
      thickness(thickness_) {
  Eigen::Matrix4f A = Eigen::Matrix4f::Zero();
  A.coeffRef(1, 0) = -1.0f;
  A.coeffRef(2, 0) = -1.0f;
  A.coeffRef(3, 0) = -1.0f;

  A.coeffRef(1, 1) = 1.0f;
  A.coeffRef(2, 2) = 1.0f;
  A.coeffRef(3, 3) = 1.0f;

  this->AtA = A.transpose() * A;

  // Determine which side of the triangle the point starts on
  glm::vec3 n = glm::cross(
      c.prevPosition - b.prevPosition,
      d.prevPosition - b.prevPosition);

  this->side =
      (glm::dot(n, a.prevPosition - b.prevPosition) >= 0.0f) ? 1.0f : -1.0f;
  n *= this->side;

  glm::vec3 triVel = barycentric.x * (b.position - b.prevPosition) +
                     barycentric.y * (c.position - c.prevPosition) +
                     barycentric.z * (d.position - d.prevPosition);
  glm::vec3 relativeVelocity = a.position - a.prevPosition - triVel;

  disp = toi * relativeVelocity +
         (1.0f - toi) * (glm::mat3(1.0f) - 1.5f * glm::outerProduct(n, n)) *
             relativeVelocity;

  this->projectedPositions[0] = a.prevPosition - 0.5f * disp;
  this->projectedPositions[1] = b.prevPosition + 0.5f * disp * barycentric.x;
  this->projectedPositions[2] = c.prevPosition + 0.5f * disp * barycentric.y;
  this->projectedPositions[3] = d.prevPosition + 0.5f * disp * barycentric.z;
}

void PointTriangleCollisionConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {
  colliding = false;

  const Node& nodeA = nodes[nodeIds[0]];
  const Node& nodeB = nodes[nodeIds[1]];
  const Node& nodeC = nodes[nodeIds[2]];
  const Node& nodeD = nodes[nodeIds[3]];

  // If the toi is 1.0f, there is going to be no rebound, so just
  // push the nodes apart to the threshold distance.

  if (toi == 1.0f) {
    this->projectedPositions[0] = nodeA.position;
    this->projectedPositions[1] = nodeB.position;
    this->projectedPositions[2] = nodeC.position;
    this->projectedPositions[3] = nodeD.position;

    glm::vec3 c = nodeC.position - nodeB.position;
    glm::vec3 d = nodeD.position - nodeB.position;
    glm::vec3 p = nodeA.position - nodeB.position;

    glm::vec3 n = this->side * glm::normalize(glm::cross(
                                   nodeC.position - nodeB.position,
                                   nodeD.position - nodeB.position));

    float nDotP = glm::dot(n, p);
    if (nDotP < thickness) {
      glm::vec3 disp = (thickness - nDotP) * n;

      this->projectedPositions[0] += disp;
      // this->projectedPositions[1] = glm::vec3(0.0f);
      // this->projectedPositions[2] = glm::vec3(0.0f);
      // this->projectedPositions[3] = glm::vec3(0.0f);

      colliding = true;
    }
  }
}

void PointTriangleCollisionConstraint::stabilizeCollisions(
    std::vector<Node>& nodes) {
  Node& nodeA = nodes[nodeIds[0]];
  Node& nodeB = nodes[nodeIds[1]];
  Node& nodeC = nodes[nodeIds[2]];
  Node& nodeD = nodes[nodeIds[3]];
}

void PointTriangleCollisionConstraint::setupTriplets(
    std::vector<Eigen::Triplet<float>>& triplets) const {
  for (uint32_t i = 0; i < 4; ++i) {
    uint32_t nodeId_i = this->nodeIds[i];
    for (uint32_t j = 0; j < 4; ++j) {
      uint32_t nodeId_j = this->nodeIds[j];
      triplets.emplace_back(
          nodeId_i,
          nodeId_j,
          this->w * this->AtA.coeff(i, j));
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
    const Node& d,
    float thickness_,
    float toi_,
    const glm::vec2& uv_,
    const glm::vec3& n_)
    : toi(toi_), nodeIds{a.id, b.id, c.id, d.id}, uv(uv_), n(n_), thickness(thickness_) {
  Eigen::Matrix4f A = Eigen::Matrix4f::Zero();
  A.coeffRef(1, 0) = -1.0f;
  A.coeffRef(2, 0) = -1.0f;
  A.coeffRef(3, 0) = -1.0f;

  A.coeffRef(1, 1) = 1.0f;
  A.coeffRef(2, 2) = 1.0f;
  A.coeffRef(3, 3) = 1.0f;

  this->AtA = A.transpose() * A;

  glm::vec3 int0Vel =
      glm::mix(a.position - a.prevPosition, b.position - b.prevPosition, uv_.x);
  glm::vec3 int1Vel =
      glm::mix(c.position - c.prevPosition, d.position - d.prevPosition, uv_.y);
  glm::vec3 relativeVelocity = int1Vel - int0Vel;

  disp = toi * relativeVelocity +
         (1.0f - toi) * (glm::mat3(1.0f) - 1.5f * glm::outerProduct(n, n)) *
             relativeVelocity;

  this->projectedPositions[0] = a.prevPosition - 0.5f * disp * (1.0f - uv_.x);
  this->projectedPositions[1] = b.prevPosition - 0.5f * disp * uv_.x;
  this->projectedPositions[2] = c.prevPosition + 0.5f * disp * (1.0f - uv_.y);
  this->projectedPositions[3] = d.prevPosition + 0.5f * disp * uv_.y;
}

void EdgeCollisionConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {
  const Node& nodeA = nodes[nodeIds[0]];
  const Node& nodeB = nodes[nodeIds[1]];
  const Node& nodeC = nodes[nodeIds[2]];
  const Node& nodeD = nodes[nodeIds[3]];

  if (toi == 1.0f) {
    this->projectedPositions[0] = nodeA.position;
    this->projectedPositions[1] = nodeB.position;
    this->projectedPositions[2] = nodeC.position;
    this->projectedPositions[3] = nodeD.position;

    glm::vec3 ab = nodeB.position - nodeA.position;
    glm::vec3 ac = nodeC.position - nodeA.position;
    glm::vec3 ad = nodeD.position - nodeA.position;

    float l0 = glm::dot(ab * uv.x, n);//glm::max(glm::dot(ab, n), 0.0f);
    float l1 = glm::dot(glm::mix(ac, ad, uv.y), n);////glm::min(glm::dot(ac, n), glm::dot(ad, n));

    if (l0 < l1 + thickness) {
      glm::vec3 disp = (l1 + thickness - l0) * n;
      this->projectedPositions[0] += 0.5f * disp * (1.0f - uv.x);
      this->projectedPositions[1] += 0.5f * disp * uv.x;
      this->projectedPositions[2] -= 0.5f * disp * (1.0f - uv.y);
      this->projectedPositions[3] -= 0.5f * disp * uv.y;
    }
  }
}

void EdgeCollisionConstraint::stabilizeCollisions(std::vector<Node>& nodes) {}

void EdgeCollisionConstraint::setupTriplets(
    std::vector<Eigen::Triplet<float>>& triplets) const {
  for (uint32_t i = 0; i < 4; ++i) {
    uint32_t nodeId_i = this->nodeIds[i];
    for (uint32_t j = 0; j < 4; ++j) {
      uint32_t nodeId_j = this->nodeIds[j];
      triplets.emplace_back(
          nodeId_i,
          nodeId_j,
          this->w * this->AtA.coeff(i, j));
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
    float thickness_,
    float toi_)
    : thickness(thickness_), toi(toi_), nodeId(node.id) {

  glm::vec3 vt = -(node.position - node.prevPosition);
  glm::vec3 n(0.0f, 1.0f, 0.0f);
  glm::vec3 disp =
      toi * vt +
      (1.0f - toi) * (glm::mat3(1.0f) - 1.5f * glm::outerProduct(n, n)) * vt;
  projectedPosition = node.prevPosition - disp;
  // // projectedPosition.y = glm::max(projectedPosition.y, thickness);
}

void StaticCollisionConstraint::setupTriplets(
    std::vector<Eigen::Triplet<float>>& triplets) const {
  triplets.emplace_back(nodeId, nodeId, this->w);
}

void StaticCollisionConstraint::projectToAuxiliaryVariable(
    const std::vector<Node>& nodes) {
  // const Node& node = nodes[nodeId];
  // projectedPosition = node.position;
  if (toi == 1.0f && nodes[nodeId].position.y < thickness) {
    projectedPosition.y = thickness;
  }
}

void StaticCollisionConstraint::setupGlobalForceVector(
    Eigen::MatrixXf& forceVector) const {
  // TODO: Do we need to declar no alias for projectedPositions or nodeIds??
  forceVector.coeffRef(nodeId, 0) += this->w * projectedPosition.x;
  forceVector.coeffRef(nodeId, 1) += this->w * projectedPosition.y;
  forceVector.coeffRef(nodeId, 2) += this->w * projectedPosition.z;
}
} // namespace Pies