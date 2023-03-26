#pragma once

#include "Constraints.h"
#include "Node.h"
#include "SpatialHash.h"
#include "Tetrahedron.h"

#include <Eigen/Core>

#include <cstdint>
#include <vector>

namespace Pies {
enum class SolverName {
  PBD,
  PD
};

struct SolverOptions {
  uint32_t iterations = 4;
  uint32_t timeSubsteps = 10;
  float fixedTimestepSize = 0.05f;
  float gravity = 10.0f;
  float damping = 0.0005f;
  float friction = 0.25f;
  float floorHeight = 0.0f;
  float gridSpacing = 1.0f;
  SolverName solver = SolverName::PD;
};

class Solver {
public:
  struct Vertex {
    glm::vec3 position{};
    float radius{};

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

  void clear();
  
  // Utilities for spawning primitives
  void createBox(const glm::vec3& translation, float scale, float stiffness);
  void createTetBox(
      const glm::vec3& translation,
      float scale,
      const glm::vec3& initialVelocity,
      float stiffness);
  void createSheet(const glm::vec3& translation, float scale, float k);

private:
  struct NodeCompRange {
    SpatialHashGridCellRange
    operator()(const Node& node, const SpatialHashGrid& grid) const;
  };

  struct TetCompRange {
    const std::vector<Node>& nodes;

    SpatialHashGridCellRange
    operator()(const Tetrahedron& node, const SpatialHashGrid& grid) const;
  };

  SolverOptions _options;
  uint32_t _constraintId = 0;

  SpatialHash<Node, NodeCompRange> _spatialHashNodes;
  SpatialHash<Tetrahedron, TetCompRange> _spatialHashTets;

  std::vector<Node> _nodes;
  std::vector<Tetrahedron> _tets;
  std::vector<PositionConstraint> _positionConstraints;
  std::vector<DistanceConstraint> _distanceConstraints;
  std::vector<TetrahedralConstraint> _tetConstraints;

  Eigen::MatrixXf _stateVector;
  Eigen::MatrixXf _stiffnessMatrix;

  // TODO: Seperate into individual buffers for each object??
  std::vector<uint32_t> _triangles;
  std::vector<uint32_t> _lines;
  std::vector<Vertex> _vertices;
};
} // namespace Pies