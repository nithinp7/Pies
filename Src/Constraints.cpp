#include "Constraints.h"

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

  glm::vec3 diff = b.position - a.position;
  float dist = glm::length(diff);

  glm::vec3 dir(1.0f, 0.0f, 0.0f);
  if (dist > 0.00001f) {
    dir = diff / dist;
  }

  float disp = this->targetDistance - dist;
  float massSum = a.mass + b.mass;

  projected[0] = a.position - disp * dir * a.mass / massSum;
  projected[1] = b.position + disp * dir * b.mass / massSum;
}

DistanceConstraint
createDistanceConstraint(uint32_t id, const Node& a, const Node& b, float w) {
  return DistanceConstraint(
      id,
      w,
      Eigen::Matrix2f::Identity(),
      Eigen::Matrix2f::Identity(),
      {glm::length(b.position - a.position)},
      {a.id, b.id});
}

void PositionConstraintProjection::operator()(
    const std::vector<Node>& nodes,
    const std::array<uint32_t, 1>& nodeIds,
    std::array<glm::vec3, 1>& projected) const {
  projected[0] = this->fixedPosition;
}

PositionConstraint createPositionConstraint(uint32_t id, const Node& node, float w) {
  return PositionConstraint(
      id, 
      w, 
      Eigen::Matrix<float, 1, 1>::Identity(),
      Eigen::Matrix<float, 1, 1>::Identity(),
      {node.position}, 
      {node.id});
}

void TetrahedralConstraintProjection::operator()(
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
  glm::mat3 Ft = glm::transpose(F);

  float stretchResistance = 1.0f;
  glm::mat3 S = Ft * F;

  // Tensor of stretch and strain constraints
  // Green - St Vernant strain formulation
  // TODO: Use sqrt version
  glm::mat3 C = S - glm::mat3(stretchResistance * stretchResistance);

  projected[0] = x1.position;
  projected[1] = x2.position;
  projected[2] = x3.position;
  projected[3] = x4.position;

  // Project the constraints one after another
  for (uint32_t i = 0; i < 3; ++i) {
    for (uint32_t j = 0; j <= i; ++j) {
      glm::mat3 dCijdX =
          glm::outerProduct(glm::column(F, j), glm::column(this->Qinv, i)) +
          glm::outerProduct(glm::column(F, i), glm::column(this->Qinv, j));
      glm::vec3 dCijdx1 = -glm::column(dCijdX, 0) - glm::column(dCijdX, 1) -
                          glm::column(dCijdX, 2);
      float denom =
          x1.mass * glm::dot(dCijdx1, dCijdx1) +
          x2.mass * glm::dot(glm::column(dCijdX, 0), glm::column(dCijdX, 0)) +
          x3.mass * glm::dot(glm::column(dCijdX, 1), glm::column(dCijdX, 1)) +
          x4.mass * glm::dot(glm::column(dCijdX, 2), glm::column(dCijdX, 2));
      float lambda = C[i][j] / denom;

      // Recompute F from projected positions?
      projected[0] += -lambda * x1.mass * dCijdx1;
      projected[1] += -lambda * x2.mass * glm::column(dCijdX, 0);
      projected[2] += -lambda * x3.mass * glm::column(dCijdX, 1);
      projected[3] += -lambda * x4.mass * glm::column(dCijdX, 2);
    }
  }
}

TetrahedralConstraint createTetrahedralConstraint(
    uint32_t id,
    float k,
    const Node& x1,
    const Node& x2,
    const Node& x3,
    const Node& x4) {

  glm::mat3 Q(
      x2.position - x1.position,
      x3.position - x1.position,
      x4.position - x1.position);

  return TetrahedralConstraint(
      id,
      k,
      Eigen::Matrix4f::Identity(),
      Eigen::Matrix4f::Identity(),
      {glm::inverse(Q)},
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

  glm::vec3 x21 = x2.position - x1.position;
  glm::vec3 x31 = x3.position - x1.position;
  glm::vec3 x41 = x4.position - x1.position;

  glm::vec3 dCdx2 = glm::cross(x31, x41) / 6.0f;
  glm::vec3 dCdx3 = glm::cross(x41, x21) / 6.0f;
  glm::vec3 dCdx4 = glm::cross(x21, x31) / 6.0f;
  glm::vec3 dCdx1 = -dCdx2 - dCdx3 - dCdx4;

  float C = glm::dot(dCdx4, x41) - this->targetVolume;
  float dCdX_magsq = glm::dot(dCdx1, dCdx1) + glm::dot(dCdx2, dCdx2) +
                     glm::dot(dCdx3, dCdx3) + glm::dot(dCdx4, dCdx4);

  float lambda = -C / dCdX_magsq;
  float massSum = x1.mass + x2.mass + x3.mass + x4.mass;

  // TODO: add in mass / massSum
  projected[0] = x1.position + lambda * dCdx1; // * x1.mass / massSum;
  projected[1] = x2.position + lambda * dCdx2; // * x2.mass / massSum;
  projected[2] = x3.position + lambda * dCdx3; // * x3.mass / massSum;
  projected[3] = x4.position + lambda * dCdx4; // * x4.mass / massSum;
}

VolumeConstraint createVolumeConstraint(
    uint32_t id,
    float k,
    const Node& x1,
    const Node& x2,
    const Node& x3,
    const Node& x4) {
  glm::vec3 x21 = x2.position - x1.position;
  glm::vec3 x31 = x3.position - x1.position;
  glm::vec3 x41 = x4.position - x1.position;

  float targetVolume = glm::dot(glm::cross(x21, x31), x41) / 6.0f;
  return VolumeConstraint(
      id, 
      k, 
      Eigen::Matrix4f::Identity(),
      Eigen::Matrix4f::Identity(),
      {targetVolume}, 
      {x1.id, x2.id, x3.id, x4.id});
}
} // namespace Pies