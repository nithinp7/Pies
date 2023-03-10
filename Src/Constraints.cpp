#include "Constraints.h"

namespace PiesForAlthea {
void DistanceConstraintProjection::operator()(
    const std::array<Node*, 2>& nodes,
    std::array<glm::vec3, 2>& projected) const {
  glm::vec3 diff = nodes[1]->position - nodes[0]->position;
  float dist = glm::length(diff);

  glm::vec3 dir(1.0f, 0.0f, 0.0f);
  if (dist > 0.00001f) {
    dir = diff / dist;
  }

  float disp = this->targetDistance - dist;
  float massSum = nodes[0]->mass + nodes[1]->mass;

  projected[0] = nodes[0]->position - disp * dir * nodes[0]->mass / massSum;
  projected[1] = nodes[1]->position + disp * dir * nodes[1]->mass / massSum;
}

DistanceConstraint createDistanceConstraint(uint32_t id, Node* a, Node* b) {
  return DistanceConstraint(
      id,
      1.0f,
      {glm::length(b->position - a->position)},
      {a, b});
}

void PositionConstraintProjection::operator()(
    const std::array<Node*, 1>& nodes,
    std::array<glm::vec3, 1>& projected) const {
  projected[0] = this->fixedPosition;
}

PositionConstraint createPositionConstraint(uint32_t id, Node* node) {
  return PositionConstraint(id, 1.0f, {node->position}, {node});
}
} // namespace PiesForAlthea