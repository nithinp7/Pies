#pragma once

#include <glm/glm.hpp>

#include <cstdint>

namespace Pies {
struct Node {
  // This ID represents the index of this node within
  // the global position and force vectors.
  // E.g., indexX = 3*id, indexY = 3*id+1, indexZ = 3*id+2
  uint32_t id = 0;

  glm::vec3 position{};
  glm::vec3 prevPosition{};
  glm::vec3 velocity{};
  glm::vec3 force{};
  float radius = 0.1f;
  float invMass = 1.0f;
};
} // namespace Pies