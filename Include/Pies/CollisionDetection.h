#pragma once

#include <glm/glm.hpp>

#include <optional>

namespace Pies {
namespace CollisionDetection {
std::optional<float> pointTriangleCCD(
    const glm::vec3& ap0,
    const glm::vec3& ab0,
    const glm::vec3& ac0,
    const glm::vec3& ap1,
    const glm::vec3& ab1,
    const glm::vec3& ac1,
    float thresholdDistance);
    
std::optional<float> edgeEdgeCCD(
    const glm::vec3& ab0,
    const glm::vec3& ac0,
    const glm::vec3& ad0,
    const glm::vec3& ab1,
    const glm::vec3& ac1,
    const glm::vec3& ad1);
} // namespace CollisionDetection
} // namespace Pies