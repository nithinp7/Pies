#pragma once

#include "Node.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <glm/glm.hpp>

#include <cstdint>

namespace Pies {
// Generalize to different types of collision
struct CollisionConstraint {
  float w = 100000.0f;
  uint32_t nodeIds[2];
  glm::vec3 projectedPositions[2];
  glm::vec3 n;

  CollisionConstraint(const Node& a, const Node& b);

  void projectToAuxiliaryVariable(const std::vector<Node>& nodes);

  void setupCollisionMatrix(Eigen::SparseMatrix<float>& systemMatrix) const;

  void setupGlobalForceVector(
      Eigen::MatrixXf& forceVector,
      uint32_t threadId,
      uint32_t threadCount) const;
};

struct PointTriangleCollisionConstraint {
  float w = 10000.0f;
  uint32_t nodeIds[4];
  glm::vec3 projectedPositions[4];
  glm::vec3 n;
  bool colliding = false;

  PointTriangleCollisionConstraint(
      const Node& a,
      const Node& b,
      const Node& c,
      const Node& d);

  void projectToAuxiliaryVariable(const std::vector<Node>& nodes);
  void setupCollisionMatrix(Eigen::SparseMatrix<float>& systemMatrix) const;
  void setupGlobalForceVector(Eigen::MatrixXf& forceVector) const;
};

struct StaticCollisionConstraint {
  float w = 1000000.0f;
  uint32_t nodeId;
  glm::vec3 projectedPosition;

  glm::vec3 n;

  StaticCollisionConstraint(
      const Node& node,
      const glm::vec3& projectedPosition);

  void setupCollisionMatrix(Eigen::SparseMatrix<float>& systemMatrix) const;

  void setupGlobalForceVector(Eigen::MatrixXf& forceVector) const;
};
} // namespace Pies