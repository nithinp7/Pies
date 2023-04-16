#pragma once

#include <glm/glm.hpp>

#include <optional>

namespace Pies {
namespace CollisionDetection {

std::optional<float> linearCCD(
    const glm::vec3& ap0,
    const glm::vec3& ab0,
    const glm::vec3& ac0,
    const glm::vec3& ap1,
    const glm::vec3& ab1,
    const glm::vec3& ac1);
} // namespace CollisionDetection
} // namespace Pies