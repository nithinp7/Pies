#pragma once

#include "CollisionConstraint.h"
#include "Constraints.h"
#include "Node.h"
#include "SpatialHash.h"
#include "Tetrahedron.h"

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <cstdint>
#include <memory>
#include <vector>

namespace Pies {
enum class SolverName { PBD, PD };

struct SolverOptions {
  uint32_t iterations = 4;
  uint32_t timeSubsteps = 3;
  float fixedTimestepSize = 0.012f;
  float gravity = 10.0f;
  float damping = 0.0005f;
  float friction = 0.01f;
  float staticFrictionThreshold = 0.0f;
  float floorHeight = 0.0f;
  float gridSpacing = 1.0f;
  uint32_t threadCount = 12;
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
  void tickPBD(float deltaTime);
  void tickPD(float deltaTime);

  const std::vector<Vertex>& getVertices() const { return this->_vertices; }

  const std::vector<uint32_t>& getLines() const { return this->_lines; }

  const std::vector<uint32_t>& getTriangles() const { return this->_triangles; }

  const SolverOptions& getOptions() const { return this->_options; }

  void clear();

  // Utilities for importing meshes
  void addNodes(const std::vector<glm::vec3>& vertices);
  
  // Utilities for spawning primitives
  void createBox(const glm::vec3& translation, float scale, float w);
  void createTetBox(
      const glm::vec3& translation,
      float scale,
      const glm::vec3& initialVelocity,
      float w,
      float mass,
      bool hinged);
  void createSheet(const glm::vec3& translation, float scale, float mass, float w);
  void createBendSheet(const glm::vec3& translation, float scale, float w);

private:
  void _computeCollisions();
  void _parallelComputeCollisions();

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
  std::vector<BendConstraint> _bendConstraints;

  uint32_t _previousNodeCount = 0;
  Eigen::MatrixXf _stateVector;
  Eigen::MatrixXf _forceVector;
  Eigen::MatrixXf _Msn_h2;
  Eigen::SparseMatrix<float> _stiffnessMatrix;
  Eigen::SparseMatrix<float> _collisionMatrix;
  Eigen::SparseMatrix<float> _stiffnessAndCollisionMatrix;
  // This isn't movable, so we keep it on the heap
  std::unique_ptr<Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>> _pLltDecomp;

  struct ThreadData {
    std::vector<CollisionConstraint> collisions;
    std::vector<StaticCollisionConstraint> staticCollisions;
  };

  std::vector<ThreadData> _threadData;

  std::vector<CollisionConstraint> _collisions;
  std::vector<StaticCollisionConstraint> _staticCollisions;

  // TODO: Seperate into individual buffers for each object??
  std::vector<uint32_t> _triangles;
  std::vector<uint32_t> _lines;
  std::vector<Vertex> _vertices;
};
} // namespace Pies