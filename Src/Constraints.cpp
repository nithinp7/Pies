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