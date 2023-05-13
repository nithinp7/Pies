#include "CollisionDetection.h"

#include <unsupported/Eigen/Polynomials>

#include <algorithm>
#include <optional>

namespace {
struct Interval {
  float start = 0.0f;
  float end = 1.0f;
};

struct CubicExpression {
  float cubicCoeff;
  float quadCoeff;
  float linearCoeff;
  float constCoeff;

  float eval(float t) const {
    float t2 = t * t;
    return cubicCoeff * t2 * t + quadCoeff * t2 + linearCoeff * t + constCoeff;
  }

  std::optional<Interval> findEarliestIntervalOfRoot() const {
    float f0 = this->eval(0.0f);
    float f1 = this->eval(1.0f);

    // Compute critical points (roots of derivative) in [0, 1]
    // Use quadratic formula on the derivative
    float a = 3 * cubicCoeff;
    float b = 2 * quadCoeff;
    float c = linearCoeff;

    // Quadratic expression
    float b2_4ac = b * b - 4.0f * a * c;
    if (a < 0.00001f || b2_4ac < 0.0f) {
      // No real critical points, the function is monotonically
      // increasing or decreasing.

      // Intermediate value theorem
      return (f0 * f1 <= 0.0f) ? std::optional<Interval>({0.0f, 1.0f})
                               : std::nullopt;
    }

    float sqrtb2_4ac = sqrt(b2_4ac);

    float t0 = (-b - sqrtb2_4ac) / (2.0f * a);
    float t1 = (-b + sqrtb2_4ac) / (2.0f * a);

    if (t0 >= 1.0f) {
      // The earliest critical point is after the interval
      // So the interval is monotonic.

      // Intermediate value theorem
      return (f0 * f1 <= 0.0f) ? std::optional<Interval>({0.0f, 1.0f})
                               : std::nullopt;

    } else if (t0 <= 0.0f) {
      // The first critical point is before the interval, so consider
      // the second critical point.
      if (t1 >= 1.0f) {
        // Second critical point is after the interval, so the interval
        // itself is monotonic.

        return (f0 * f1 <= 0.0f) ? std::optional<Interval>({0.0f, 1.0f})
                                 : std::nullopt;
      } else if (t1 <= 0.0f) {
        // All critical points are before the interval, so we have no roots
        return std::nullopt;
      } else  {
        // We have two intervals to check here, [0,t1] and [t1,1]
        float ft1 = this->eval(t1);

        if (f0 * ft1 <= 0.0f) {
          return Interval{0.0f, t1};
        } else {
          return (ft1 * f1 <= 0.0f) ? std::optional<Interval>({t1, 1.0f})
                                    : std::nullopt;
        }
      }
    } else {
      // The first critical point is in the interval
      // First check [0,t0]
      float ft0 = this->eval(t0);
      if (f0 * ft0 <= 0.0f) {
        return Interval{0, t0};
      }

      if (t1 >= 1.0f) {
        // Second critical point is after the interval, so check [t0,1]
        return (ft0 * f1 <= 0.0f) ? std::optional<Interval>({t0, 1.0f})
                                  : std::nullopt;
      } else {
        // Second critical point is also within the interval, need to check
        // [t0,t1] and [t1,1]
        float ft1 = this->eval(t1);

        if (ft0 * ft1 <= 0.0f) {
          return Interval{t0, t1};
        } else if (ft1 * f1 <= 0.0f) {
          return Interval{t1, 1.0f};
        } else {
          return std::nullopt;
        }
      }
    }
  }

  // Based on:
  // High-Performance Polynomial Root Finding For Graphics, Cem Yuksel
  std::optional<float> fastFindRootInInterval() const {
    std::optional<Interval> interval = findEarliestIntervalOfRoot();
    if (!interval) {
      return std::nullopt;
    }

    float derivQuadCoeff = 3 * this->cubicCoeff;
    float derivLinearCoeff = 2 * this->quadCoeff;
    float derivConst = this->linearCoeff;

    const uint32_t NEWTON_ITERS = 20;
    const float EPS = 0.001f;

    // Initial guess half-way in the interval, hopefully sufficiently far away from
    // critical points.
    float t = interval->start;
    for (uint32_t i = 0; i < NEWTON_ITERS; ++i) {
      // Derivative evaluated at t
      float fpt =
          derivQuadCoeff * t * t + derivLinearCoeff * t + derivConst;
      float ft = this->eval(t);

      // TODO: Handle fpt close to 0??
      float tn = glm::clamp(t - ft / fpt, interval->start, interval->end);
      if (glm::abs(tn - t) < EPS) {
        return tn;
      }

      t = tn;
    }

    return t;
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

std::optional<float> pointTriangleCCD(
    const glm::vec3& ap0,
    const glm::vec3& ab0,
    const glm::vec3& ac0,
    const glm::vec3& ap1,
    const glm::vec3& ab1,
    const glm::vec3& ac1,
    float thresholdDistance) {

  // Early check to see if the point ever crosses the triangle plane
  // Assumes normal doesn't rotate significantly over the short interval.
  glm::vec3 n0 = glm::cross(ab0, ac0);
  glm::vec3 n1 = glm::cross(ab1, ac1);
  float nDotP0 = glm::dot(n0, ap0);
  float nDotP1 = glm::dot(n1, ap1);

  if (nDotP0 * nDotP1 >= 0.0f) {
    // The point never crosses the plane of the triangle. We still consider it
    // to be colliding if it is within the thickness.

    float n1mag = glm::length(n1);
    nDotP1 /= n1mag;
    n1 /= n1mag;

    // Considers points behind the triangle as well...
    if (glm::abs(nDotP1) < thresholdDistance) {
      glm::vec3 barycentricCoords = glm::inverse(glm::mat3(ab1, ac1, n1)) * ap1;

      if ((0.0 > barycentricCoords.x) || (barycentricCoords.x > 1.0) ||
          (0.0 > barycentricCoords.y) || (barycentricCoords.y > 1.0) ||
          (barycentricCoords.x + barycentricCoords.y > 1.0)) {
        return std::nullopt;
      }

      return 1.0f;
    }

    return std::nullopt;
  }

  glm::vec3 apd = ap1 - ap0;
  glm::vec3 abd = ab1 - ab0;
  glm::vec3 acd = ac1 - ac0;

  // Look for a ray-plane intersection throughout the interval
  CubicExpression expression{};
  expandTerm(ap0.x, ab0.y, ac0.z, apd.x, abd.y, acd.z, expression);
  expandTerm(-ap0.x, ac0.y, ab0.z, -apd.x, acd.y, abd.z, expression);
  expandTerm(-ab0.x, ap0.y, ac0.z, -abd.x, apd.y, acd.z, expression);
  expandTerm(ab0.x, ac0.y, ap0.z, abd.x, acd.y, apd.z, expression);
  expandTerm(ac0.x, ap0.y, ab0.z, acd.x, apd.y, abd.z, expression);
  expandTerm(-ac0.x, ab0.y, ap0.z, -acd.x, abd.y, apd.z, expression);

  std::optional<float> t = expression.fastFindRootInInterval();
  //std::optional<float> t = expression.findRootInInterval();

  if (!t) {
    // CCD failed to find a point-plane intersection. Static collision has
    // already previously failed.
    return std::nullopt;
  }

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
}

std::optional<float> edgeEdgeCCD(
    const glm::vec3& ab0,
    const glm::vec3& ac0,
    const glm::vec3& ad0,
    const glm::vec3& ab1,
    const glm::vec3& ac1,
    const glm::vec3& ad1) {
  // TODO: Conservative test for early exit??
  // Find the shortest distance between two lines
  {
    glm::vec3 cd1 = ad1 - ac1;

    float abMagSq = glm::dot(ab1, ab1);
    float cdMagSq = glm::dot(cd1, cd1);
    float abDotCd = glm::dot(ab1, cd1);

    float acDotAb = glm::dot(ac1, ab1);
    float acDotCd = glm::dot(ac1, cd1);

    float det = abMagSq * -cdMagSq + abDotCd * abDotCd;
    float u = 0.0f;
    float v = 0.0f;
    if (det != 0.0f) {
      det = 1.0f / det;
      float u = (acDotAb * -cdMagSq + abDotCd * acDotCd) * det;
      float v = (abMagSq * acDotCd - acDotAb * abDotCd) * det;
    } else {
      float u0 = 0.0f;
      float u1 = 1.0f;
      float v0 = glm::dot(ac1, ab1);
      float v1 = glm::dot(ad1, ab1);

      bool flip0 = false;
      bool flip1 = false;

      if (u0 > u1) {
        std::swap(u0, u1);
        flip0 = true;
      }

      if (v0 > v1) {
        std::swap(v0, v1);
        flip1 = true;
      }

      if (u0 >= v1) {
        u = flip0 ? 1.0f : 0.0f;
        v = flip1 ? 0.0f : 1.0f;
      } else if (v0 >= u1) {
        u = flip0 ? 0.0f : 1.0f;
        v = flip1 ? 1.0f : 0.0f;
      } else {
        float mid = (u0 > v0) ? (u0 + v1) * 0.5f : (v0 + u1) * 0.5f;
        u = (u0 == u1) ? 0.5f : (mid - u0) / (u1 - u0);
        v = (v0 == v1) ? 0.5f : (mid - v0) / (v1 - v0);
      }
    }

    u = glm::clamp(u, 0.0f, 1.0f);
    v = glm::clamp(v, 0.0f, 1.0f);

    glm::vec3 q0 = glm::mix(glm::vec3(0.0f), ab1, u);
    glm::vec3 q1 = glm::mix(ac1, ad1, v);

    glm::vec3 n = q0 - q1;
    float dist = glm::length(n);
    n /= dist;

    float thickness = 0.5f;
    if (dist < thickness) {
      return 1.0f;
    }
  }

  glm::vec3 abd = ab1 - ab0;
  glm::vec3 acd = ac1 - ac0;
  glm::vec3 add = ad1 - ad0;

  // Look for a ray-plane intersection throughout the interval
  CubicExpression expression{};
  expandTerm(ab0.x, ac0.y, ad0.z, abd.x, acd.y, add.z, expression);
  expandTerm(-ab0.x, ad0.y, ac0.z, -abd.x, add.y, acd.z, expression);
  expandTerm(-ac0.x, ab0.y, ad0.z, -acd.x, abd.y, add.z, expression);
  expandTerm(ac0.x, ad0.y, ab0.z, acd.x, add.y, abd.z, expression);
  expandTerm(ad0.x, ab0.y, ac0.z, add.x, abd.y, acd.z, expression);
  expandTerm(-ad0.x, ac0.y, ab0.z, -add.x, acd.y, abd.z, expression);

  std::optional<float> t = expression.findRootInInterval();

  // Check
  if (!t) {
    return std::nullopt;
  }

  glm::vec3 abt = ab0 + abd * *t;
  glm::vec3 act = ac0 + acd * *t;
  glm::vec3 adt = ad0 + add * *t;
  glm::vec3 cdt = adt - act;

  glm::vec3 nt = glm::normalize(glm::cross(abt, cdt));

  // Check for 2D line intersection at time t

  // abt * u = act + cdt * v
  // abt * u - cdt * v = act
  // [abt, -cdt] [u,v]t = act
  glm::mat3 M(abt, -cdt, nt);
  glm::vec2 uv(glm::inverse(M) * act);

  if (uv.x < 0.0f || uv.x > 1.0f || uv.y < 0.0f || uv.y > 1.0f) {
    return std::nullopt;
  }

  return t;
}
} // namespace CollisionDetection
} // namespace Pies