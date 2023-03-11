#pragma once

#include "Constraints.h"
#include "Node.h"

#include <cstdint>
#include <vector>

namespace Pies {
class Solver {
public:
  Solver();
  void tick(float deltaTime);

  const std::vector<glm::vec3>& getVertices() const { return this->_vertices; }

  const std::vector<uint32_t>& getDistanceConstraintLines() const {
    return this->_distanceConstraintLines;
  }

private:
  std::vector<Node> _nodes;
  std::vector<PositionConstraint> _positionConstraints;
  std::vector<DistanceConstraint> _distanceConstraints;

  // Scratch vertices vector to avoid reallocation
  std::vector<glm::vec3> _vertices;

  // Indices for lines representing the distance constraints.
  std::vector<uint32_t> _distanceConstraintLines;
};
} // namespace Pies