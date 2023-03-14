#include "Constraints.h"

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
createDistanceConstraint(uint32_t id, const Node& a, const Node& b) {
  return DistanceConstraint(
      id,
      1.0f,
      {glm::length(b.position - a.position)},
      {a.id, b.id});
}

void PositionConstraintProjection::operator()(
    const std::vector<Node>& nodes,
    const std::array<uint32_t, 1>& nodeIds,
    std::array<glm::vec3, 1>& projected) const {
  projected[0] = this->fixedPosition;
}

PositionConstraint createPositionConstraint(uint32_t id, const Node& node) {
  return PositionConstraint(id, 1.0f, {node.position}, {node.id});
}
} // namespace Pies