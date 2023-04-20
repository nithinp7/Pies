#include "CollisionDetection.h"

#include <unsupported/Eigen/Polynomials>

#include <algorithm>
#include <optional>

namespace {
struct CubicExpression {
  float cubicCoeff;
  float quadCoeff;
  float linearCoeff;
  float constCoeff;

  float eval(float t) const {
    float t2 = t * t;
    return cubicCoeff * t2 * t + quadCoeff * t2 + linearCoeff * t + constCoeff;
  }

  std::optional<float> findRootInInterval() const {
    // TODO: Epsilon 0 check??
    if (cubicCoeff == 0.0f) {
      if (quadCoeff == 0.0f) {
        if (linearCoeff == 0.0f) {
          if (constCoeff == 0.0f) {
            return 0.0f;
          }

          return std::nullopt;
        }

        float t = -constCoeff / linearCoeff;
        if (t >= 0.0f && t <= 1.0f) {
          return t;
        }

        return std::nullopt;
      } else {
        // Quadratic expression
        float b2_4ac =
            linearCoeff * linearCoeff - 4.0f * quadCoeff * constCoeff;
        if (b2_4ac < 0.0f) {
          // Imaginary roots
          return std::nullopt;
        }

        float sqrtb2_4ac = sqrt(b2_4ac);

        float t = (-linearCoeff - sqrtb2_4ac) / (2.0f * quadCoeff);
        if (t > 1.0f) {
          return std::nullopt;
        } else if (t < 0.0f) {
          // Consider the second root if the first root falls before the
          // interval.
          t = (-linearCoeff + sqrtb2_4ac) / (2.0f * quadCoeff);
        }

        if (t >= 0.0f && t <= 1.0f) {
          return t;
        }

        return std::nullopt;
      }
    }

    Eigen::Vector4f coeffs(constCoeff, linearCoeff, quadCoeff, cubicCoeff);
    Eigen::PolynomialSolver<float, 3> solver(coeffs);
    const Eigen::PolynomialSolver<float, 3>::RootsType& roots = solver.roots();

    std::optional<float> validRoot = std::nullopt;
    for (uint32_t i = 0; i < 3; ++i) {
      // Is it better to just check if the imaginary part is eq to 0?
      if (std::abs(roots[i].imag()) < 0.0000001f && roots[i].real() >= 0.0f &&
          roots[i].real() <= 1.0f) {
        if (!validRoot || roots[i].real() < *validRoot) {
          validRoot = roots[i].real();
        }
      }
    }

    return validRoot;
  }
};

// TODO: Better name??
inline void expandTerm(
    float a0,
    float b0,
    float c0,
    float ad,
    float bd,
    float cd,
    CubicExpression& expression) {
  expression.cubicCoeff += ad * bd * cd;
  expression.quadCoeff += ad * bd * c0 + a0 * bd * cd + ad * b0 * cd;
  expression.linearCoeff += ad * b0 * c0 + a0 * bd * c0 + a0 * b0 * cd;
  expression.constCoeff += a0 * b0 * c0;
}
} // namespace

namespace Pies {
namespace CollisionDetection {

std::optional<float> linearCCD(
    const glm::vec3& ap0,
    const glm::vec3& ab0,
    const glm::vec3& ac0,
    const glm::vec3& ap1,
    const glm::vec3& ab1,
    const glm::vec3& ac1) {

  // Early check to see if the point ever crosses the triangle plane
  // Assumes normal doesn't rotate significantly over the short interval.
  glm::vec3 n0 = glm::normalize(glm::cross(ab0, ac0));
  glm::vec3 n1 = glm::normalize(glm::cross(ab1, ac1));
  float nDotP0 = glm::dot(n0, ap0);
  float nDotP1 = glm::dot(n1, ap1);

  glm::vec3 apd = ap1 - ap0;
  glm::vec3 abd = ab1 - ab0;
  glm::vec3 acd = ac1 - ac0;

  std::optional<float> t = std::nullopt;
  if (nDotP0 * nDotP1 < 0.0f) { 
    // Look for a ray-plane intersection throughout the interval
    CubicExpression expression{};
    expandTerm(ap0.x, ab0.y, ac0.z, apd.x, abd.y, acd.z, expression);
    expandTerm(-ap0.x, ac0.y, ab0.z, -apd.x, acd.y, abd.z, expression);
    expandTerm(-ab0.x, ap0.y, ac0.z, -abd.x, apd.y, acd.z, expression);
    expandTerm(ab0.x, ac0.y, ap0.z, abd.x, acd.y, apd.z, expression);
    expandTerm(ac0.x, ap0.y, ab0.z, acd.x, apd.y, abd.z, expression);
    expandTerm(-ac0.x, ab0.y, ap0.z, -acd.x, abd.y, apd.z, expression);

    t = expression.findRootInInterval();
  }

  if (t) {
    // Now that we have a t-value for the ray-plane intersection, check if
    // the ray hits the triangle particularly
    glm::vec3 apt = ap0 + *t * apd;
    glm::vec3 abt = ab0 + *t * abd;
    glm::vec3 act = ac0 + *t * acd;

    glm::vec3 n = glm::normalize(glm::cross(abt, act));

    glm::vec3 barycentricCoords = glm::inverse(glm::mat3(abt, act, n)) * apt;

    if ((0.0 > barycentricCoords.x) || (barycentricCoords.x > 1.0) ||
        (0.0 > barycentricCoords.y) || (barycentricCoords.y > 1.0) ||
        (barycentricCoords.x + barycentricCoords.y > 1.0)) {
      return std::nullopt;
    }

    return t;
  } else {
    // Linear CCD failed, test if the point is closer than a desired thickness
    // perpendicular to the triangle normal _at the end of the interval_.

    // TODO: Should we consider points _behind_ the triangle??
    float thickness = 0.4f;
    if (nDotP1 >= 0.0f && nDotP1 < thickness) {
      glm::vec3 barycentricCoords = glm::inverse(glm::mat3(ab1, ac1, n1)) * ap1;

      if ((0.0 > barycentricCoords.x) || (barycentricCoords.x > 1.0) ||
          (0.0 > barycentricCoords.y) || (barycentricCoords.y > 1.0) ||
          (barycentricCoords.x + barycentricCoords.y > 1.0)) {
        return std::nullopt;
      }

      return 0.0f;
    }

    return std::nullopt;
  }
}
} // namespace CollisionDetection
} // namespace Pies