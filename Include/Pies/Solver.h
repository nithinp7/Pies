#pragma once

#include "Constraints.h"
#include "Node.h"

#include <cstdint>
#include <vector>

namespace Pies {
struct SolverOptions {
  uint32_t iterations = 4;
  uint32_t timeSubsteps = 10;
  float fixedTimestepSize = 0.05f;
  float gravity = 10.0f;
  float damping = 0.001f;
  float friction = 0.1f;
};

class Solver {
public:
  struct Vertex {
    glm::vec3 position{};
    glm::vec3 baseColor{};
    float roughness{};
    float metallic{};
  };

  bool renderStateDirty = true;
  bool releaseHinge = false;

  Solver() = default;
  Solver(const SolverOptions& options);

  void tick(float deltaTime);

  const std::vector<Vertex>& getVertices() const { return this->_vertices; }

  const std::vector<uint32_t>& getLines() const { return this->_lines; }

  const std::vector<uint32_t>& getTriangles() const { return this->_triangles; }

  const SolverOptions& getOptions() const { return this->_options; }

  // Utilities for spawning primitives
  void createBox(const glm::vec3& translation, float scale, float k);

private:
  SolverOptions _options;
  uint32_t _constraintId = 0;

  std::vector<Node> _nodes;
  std::vector<PositionConstraint> _positionConstraints;
  std::vector<DistanceConstraint> _distanceConstraints;

  // TODO: Seperate into individual buffers for each object??
  std::vector<uint32_t> _triangles;
  std::vector<uint32_t> _lines;
  std::vector<Vertex> _vertices;
};
} // namespace Pies