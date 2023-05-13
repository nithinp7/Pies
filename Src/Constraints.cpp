#include "Constraints.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace Pies {
void DistanceConstraintProjection::operator()(
    const std::vector<Node>& nodes,
    const std::array<uint32_t, 2>& nodeIds,
    std::array<glm::vec3, 2>& projected) const {
  const Node& a = nodes[nodeIds[0]];
  const Node& b = nodes[nodeIds[1]];

  projected[0] = a.position;
  projected[1] = b.position;

  glm::vec3 diff = b.position - a.position;
  float dist = glm::length(diff);

  glm::vec3 dir(1.0f, 0.0f, 0.0f);
  if (dist > 0.00001f) {
    dir = diff / dist;
  }

  float disp = this->targetDistance -
               dist; // glm::min(this->targetDistance - dist, 0.0f);

  float wSum = a.invMass + b.invMass;

  projected[0] += -disp * dir; // a.position - disp * dir * a.invMass / wSum;
  // projected[1] += glm::vec3(0.0f); // b.position + disp * dir * b.invMass /
  // wSum;
}

DistanceConstraint
createDistanceConstraint(uint32_t id, const Node& a, const Node& b, float w) {
  // A == B
  Eigen::Matrix2f A;
  //  = Eigen::Matrix2f::Zero();
  A.coeffRef(0, 0) = 0.5f;
  A.coeffRef(0, 1) = -0.5f;
  A.coeffRef(1, 0) = -0.5f;
  A.coeffRef(1, 1) = 0.5f;

  return DistanceConstraint(
      id,
      w,
      A,
      A,
      {glm::length(b.position - a.position)},
      {a.id, b.id});
}

void PositionConstraintProjection::operator()(
    const std::vector<Node>& nodes,
    const std::array<uint32_t, 1>& nodeIds,
    std::array<glm::vec3, 1>& projected) const {
  projected[0] = this->fixedPosition;
}

PositionConstraint
createPositionConstraint(uint32_t id, const Node& node, float w) {
  return PositionConstraint(
      id,
      w,
      Eigen::Matrix<float, 1, 1>::Identity(),
      Eigen::Matrix<float, 1, 1>::Identity(),
      {node.position, glm::vec3(0.0f)},
      {node.id});
}

static glm::vec3
computeD(const glm::vec3& sigma, float omegaMin, float omegaMax) {
  const uint32_t COMP_D_ITERS = 10;
  glm::vec3 D(0.0f);
  for (uint32_t i = 0; i < COMP_D_ITERS; ++i) {
    glm::vec3 sigmaPlusD = sigma + D;
    float product = sigmaPlusD.x * sigmaPlusD.y * sigmaPlusD.z;
    float omega = glm::clamp(product, omegaMin, omegaMax);
    float C = product - omega;
    glm::vec3 gradC(
        sigmaPlusD.y * sigmaPlusD.z,
        sigmaPlusD.x * sigmaPlusD.z,
        sigmaPlusD.x * sigmaPlusD.y);
    D = (glm::dot(gradC, D) - C) * gradC / glm::dot(gradC, gradC);
  }

  return D;
}

/*void TetrahedralConstraintProjection::operator()(
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

  Eigen::Matrix3f F_;
  F_ << F[0][0], F[0][1], F[0][2], F[1][0], F[1][1], F[1][2], F[2][0], F[2][1],
      F[2][2];

  Eigen::JacobiSVD<Eigen::Matrix3f> svdF(
      F_,
      Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Vector3f singularValues = svdF.singularValues();
  glm::vec3 D = computeD(
      glm::vec3(singularValues[0], singularValues[1], singularValues[2]),
      this->minOmega,
      this->maxOmega);
  singularValues[0] += D.x;
  singularValues[1] += D.y;
  singularValues[2] += D.z;
  
  for (uint32_t i = 0; i < 3; ++i) {
    singularValues[i] =
        glm::clamp(singularValues[i], this->minStrain, this->maxStrain);
  }
  
  if (glm::determinant(F) < 0.0f) {
    singularValues[2] *= -1.0f;
  }

  // The "fixed" deformation gradient
  Eigen::Matrix3f Fhat =
      svdF.matrixU() * singularValues.asDiagonal() * svdF.matrixV().transpose();
  glm::mat3 P1 = glm::transpose(glm::mat3(
      Fhat.coeff(0, 0),
      Fhat.coeff(1, 0),
      Fhat.coeff(2, 0),
      Fhat.coeff(0, 1),
      Fhat.coeff(1, 1),
      Fhat.coeff(2, 1),
      Fhat.coeff(0, 2),
      Fhat.coeff(1, 2),
      Fhat.coeff(2, 2)));

  projected[0] = glm::vec3(0.0f);
  projected[1] = P1[0];
  projected[2] = P1[1];
  projected[3] = P1[2];
}*/

TetrahedralConstraint createTetrahedralConstraint(
    uint32_t id,
    float w,
    const Node& x1,
    const Node& x2,
    const Node& x3,
    const Node& x4,
    float minStrain,
    float maxStrain,
    float compression,
    float stretching) {

  // Converts world positions to differential coords
  Eigen::Matrix<float, 3, 4> worldToDiff = Eigen::Matrix<float, 3, 4>::Zero();
  worldToDiff.coeffRef(0, 0) = -1.0f;
  worldToDiff.coeffRef(1, 0) = -1.0f;
  worldToDiff.coeffRef(2, 0) = -1.0f;

  worldToDiff.coeffRef(0, 1) = 1.0f;
  worldToDiff.coeffRef(1, 2) = 1.0f;
  worldToDiff.coeffRef(2, 3) = 1.0f;

  // Converts barycentric coords to world differential cords
  glm::mat3 baryToDiff(
      x2.position - x1.position,
      x3.position - x1.position,
      x4.position - x1.position);
  glm::mat3 diffToBary = glm::inverse(baryToDiff);

  Eigen::Matrix3f diffToBary_;
  diffToBary_ << diffToBary[0][0], diffToBary[0][1], diffToBary[0][2],
      diffToBary[1][0], diffToBary[1][1], diffToBary[1][2], diffToBary[2][0],
      diffToBary[2][1], diffToBary[2][2];

  Eigen::Matrix<float, 3, 4> A_ = diffToBary_ * worldToDiff;
  Eigen::Matrix4f A = Eigen::Matrix4f::Zero();
  // A.coeffRef(1, 0) = -1.0f;
  // A.coeffRef(2, 0) = -1.0f;
  // A.coeffRef(3, 0) = -1.0f;

  // A.coeffRef(1, 1) = 1.0f;
  // A.coeffRef(2, 2) = 1.0f;
  // A.coeffRef(3, 3) = 1.0f;

  A.row(0) << 0.0f, 0.0f, 0.0f, 0.0f;
  A.row(1) = A_.row(0);
  A.row(2) = A_.row(1);
  A.row(3) = A_.row(2);

  return TetrahedralConstraint(
      id,
      w,
      A,
      Eigen::Matrix4f::Identity(),
      {baryToDiff, diffToBary, minStrain, maxStrain, compression, stretching},
      {x1.id, x2.id, x3.id, x4.id});
}

void VolumeConstraintProjection::operator()(
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
  Eigen::Matrix3f F_;
  F_ << F[0][0], F[0][1], F[0][2], F[1][0], F[1][1], F[1][2], F[2][0], F[2][1],
      F[2][2];

  Eigen::JacobiSVD<Eigen::Matrix3f> svdF(
      F_,
      Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Vector3f singularValues = svdF.singularValues();
  glm::vec3 D = computeD(
      glm::vec3(singularValues[0], singularValues[1], singularValues[2]),
      this->minOmega,
      this->maxOmega);
  singularValues[0] += D.x;
  singularValues[1] += D.y;
  singularValues[2] += D.z;

  // The "fixed" deformation gradient
  Eigen::Matrix3f Fhat =
      svdF.matrixU() * singularValues.asDiagonal() * svdF.matrixV().transpose();
  glm::mat3 P1 = glm::transpose(glm::mat3(
      Fhat.coeff(0, 0),
      Fhat.coeff(1, 0),
      Fhat.coeff(2, 0),
      Fhat.coeff(0, 1),
      Fhat.coeff(1, 1),
      Fhat.coeff(2, 1),
      Fhat.coeff(0, 2),
      Fhat.coeff(1, 2),
      Fhat.coeff(2, 2)));

  projected[0] = glm::vec3(0.0f);
  projected[1] = P1[0];
  projected[2] = P1[1];
  projected[3] = P1[2];
}

VolumeConstraint createVolumeConstraint(
    uint32_t id,
    float w,
    const Node& x1,
    const Node& x2,
    const Node& x3,
    const Node& x4,
    float compression,
    float stretching) {
  // Converts world positions to differential coords
  Eigen::Matrix<float, 3, 4> worldToDiff = Eigen::Matrix<float, 3, 4>::Zero();
  worldToDiff.coeffRef(0, 0) = -1.0f;
  worldToDiff.coeffRef(1, 0) = -1.0f;
  worldToDiff.coeffRef(2, 0) = -1.0f;

  worldToDiff.coeffRef(0, 1) = 1.0f;
  worldToDiff.coeffRef(1, 2) = 1.0f;
  worldToDiff.coeffRef(2, 3) = 1.0f;

  // Converts barycentric coords to world differential cords
  glm::mat3 baryToDiff(
      x2.position - x1.position,
      x3.position - x1.position,
      x4.position - x1.position);
  glm::mat3 diffToBary = glm::inverse(baryToDiff);

  Eigen::Matrix3f diffToBary_;
  diffToBary_ << diffToBary[0][0], diffToBary[0][1], diffToBary[0][2],
      diffToBary[1][0], diffToBary[1][1], diffToBary[1][2], diffToBary[2][0],
      diffToBary[2][1], diffToBary[2][2];

  Eigen::Matrix<float, 3, 4> A_ = diffToBary_ * worldToDiff;
  Eigen::Matrix4f A = Eigen::Matrix4f::Zero();
  // A.coeffRef(1, 0) = -1.0f;
  // A.coeffRef(2, 0) = -1.0f;
  // A.coeffRef(3, 0) = -1.0f;

  // A.coeffRef(1, 1) = 1.0f;
  // A.coeffRef(2, 2) = 1.0f;
  // A.coeffRef(3, 3) = 1.0f;

  A.row(0) << 0.0f, 0.0f, 0.0f, 0.0f;
  A.row(1) = A_.row(0);
  A.row(2) = A_.row(1);
  A.row(3) = A_.row(2);

  return VolumeConstraint(
      id,
      w,
      A,
      Eigen::Matrix4f::Identity(),
      {diffToBary, compression, stretching},
      {x1.id, x2.id, x3.id, x4.id});
}

void BendConstraintProjection::operator()(
    const std::vector<Node>& nodes,
    const std::array<uint32_t, 4>& nodeIds,
    std::array<glm::vec3, 4>& projected) const {

  const Node& x1 = nodes[nodeIds[0]];
  const Node& x2 = nodes[nodeIds[1]];
  const Node& x3 = nodes[nodeIds[2]];
  const Node& x4 = nodes[nodeIds[3]];

  glm::vec3 p2 = x2.position - x1.position;
  glm::vec3 p3 = x3.position - x1.position;
  glm::vec3 p4 = x4.position - x1.position;

  glm::vec3 p2Xp3 = glm::cross(p2, p3);
  glm::vec3 p2Xp4 = glm::cross(p2, p4);

  float p2Xp3_len = glm::length(p2Xp3);
  float p2Xp4_len = glm::length(p2Xp4);

  // TODO: Divide by zero check for degenerate triangles
  glm::vec3 n1 = p2Xp3 / p2Xp3_len;
  glm::vec3 n2 = p2Xp4 / p2Xp4_len;

  float d = glm::dot(n1, n2);
  float d2 = d * d;

  float C = acos(d) - this->initialAngle;
  projected[0] = x1.position;
  projected[1] = x2.position;
  projected[2] = x3.position;
  projected[3] = x4.position;

  glm::vec3 q3 = (glm::cross(p2, n2) + (glm::cross(n1, p2) * d)) / p2Xp3_len;
  glm::vec3 q4 = (glm::cross(p2, n1) + (glm::cross(n2, p2) * d)) / p2Xp4_len;
  glm::vec3 q2 =
      -((glm::cross(p3, n2) + (glm::cross(n1, p3) * d)) / p2Xp3_len) -
      ((glm::cross(p4, n1) + (glm::cross(n2, p4) * d)) / p2Xp4_len);
  glm::vec3 q1 = -q2 - q3 - q4;

  float wSum = x1.invMass + x2.invMass + x3.invMass + x4.invMass;
  float qSquaredSum =
      glm::dot(q1, q1) + glm::dot(q2, q2) + glm::dot(q3, q3) + glm::dot(q4, q4);
  float projNumerator = sqrt(glm::max(1.0f - d2, 0.0f)) * C;

  if (qSquaredSum < 0.00001f) {
    return;
  }

  // Based on Bending Constraint Projection in Appendix A of PBD 2007 Paper
  projected[0] += -q1 * (4 * x1.invMass / wSum) * (projNumerator) / qSquaredSum;
  projected[1] += -q2 * (4 * x2.invMass / wSum) * (projNumerator) / qSquaredSum;
  projected[2] += -q3 * (4 * x3.invMass / wSum) * (projNumerator) / qSquaredSum;
  projected[3] += -q4 * (4 * x4.invMass / wSum) * (projNumerator) / qSquaredSum;
}

BendConstraint createBendConstraint(
    uint32_t id,
    float w,
    const Node& x1,
    const Node& x2,
    const Node& x3,
    const Node& x4) {

  // Node x2 and x3 comprise the shared edge of the adjacent triangles

  glm::vec3 p2 = x2.position - x1.position;
  glm::vec3 p3 = x3.position - x1.position;
  glm::vec3 p4 = x4.position - x1.position;

  glm::vec3 n1 = glm::normalize(glm::cross(p2, p3));
  glm::vec3 n2 = glm::normalize(glm::cross(p2, p4));

  float targetAngle = acos(glm::dot(n1, n2));

  return BendConstraint(
      id,
      w,
      Eigen::Matrix4f::Identity(),
      Eigen::Matrix4f::Identity(),
      {targetAngle},
      {x1.id, x2.id, x3.id, x4.id});
}

} // namespace Pies