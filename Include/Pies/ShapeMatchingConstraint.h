#pragma once

#include "Node.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <glm/glm.hpp>

#include <cstdint>
#include <vector>

namespace Pies {
class ShapeMatchingConstraint {
public:
  ShapeMatchingConstraint(
      const std::vector<uint32_t>& indices,
      const std::vector<glm::vec3>& materialCoordinates,
      float w);

  void
  setupGlobalStiffnessMatrix(Eigen::SparseMatrix<float>& systemMatrix) const;
  void setupGlobalForceVector(Eigen::MatrixXf& forceVector) const;
  void projectToAuxiliaryVariable(const std::vector<Node>& nodes);

private:
  std::vector<uint32_t> _nodeIndices;

  Eigen::MatrixXf _currentPositions;
  Eigen::MatrixXf _materialCoordinates;

  Eigen::MatrixXf _projectedPositions;
  
  float _w;
};
} // namespace Pies