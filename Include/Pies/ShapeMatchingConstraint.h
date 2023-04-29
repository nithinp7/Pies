#pragma once

#include "Node.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <glm/glm.hpp>

#include <cstdint>
#include <vector>

namespace Pies {
class ShapeMatchingConstraint {
public:
  ShapeMatchingConstraint(
      const std::vector<Node>& nodes,
      const std::vector<uint32_t>& indices,
      const std::vector<glm::vec3>& materialCoordinates,
      float w);

  void
  setupGlobalStiffnessMatrix(Eigen::SparseMatrix<float>& systemMatrix) const;
  void setupGlobalForceVector(Eigen::MatrixXf& forceVector) const;
  void projectToAuxiliaryVariable(const std::vector<Node>& nodes);

private:
  std::vector<uint32_t> _nodeIndices;
  Eigen::MatrixXd _materialCoords;

  Eigen::Matrix3d _Qinv;

  Eigen::MatrixXd _projectedPositions;
  Eigen::Quaterniond _currentRotation;
  
  float _w;
};

class GoalMatchingConstraint {
public:
  GoalMatchingConstraint(
      const std::vector<Node>& nodes,
      const std::vector<uint32_t>& indices,
      float w);

  void
  setupGlobalStiffnessMatrix(Eigen::SparseMatrix<float>& systemMatrix) const;
  void setupGlobalForceVector(Eigen::MatrixXf& forceVector) const;
  void projectToAuxiliaryVariable(const std::vector<Node>& nodes);
  void setTransform(const glm::mat4& transform);

private:
  std::vector<uint32_t> _nodeIndices;
  std::vector<glm::vec3> _materialCoords;
  Eigen::MatrixXd _projectedPositions;
  glm::mat4 _currentTransform;
  
  float _w;
};
} // namespace Pies