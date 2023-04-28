#pragma once

#include "CollisionConstraint.h"
#include "Constraints.h"
#include "Node.h"
#include "ShapeMatchingConstraint.h"
#include "SpatialHash.h"
#include "Tetrahedron.h"
#include "Triangle.h"

#include <Eigen/Core>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <cstdint>
#include <memory>
#include <thread>
#include <vector>

namespace Pies {
enum class SolverName { PBD, PD };

struct SolverOptions {
  float fixedTimestepSize = 0.012f;
  uint32_t timeSubsteps = 1;
  uint32_t iterations = 4;
  uint32_t collisionStabilizationIterations = 4;
  float collisionThresholdDistance = 0.1f;
  float collisionThickness = 0.05f;
  float gravity = 10.0f;
  float damping = 0.006f;
  float friction = 0.01f;              // 1f;
  float staticFrictionThreshold = 0.f; // 1.0f;
  float floorHeight = 0.0f;
  float gridSpacing = 2.0f;
  uint32_t threadCount = 8;
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
  Solver(Solver&& rhs) = default;
  Solver& operator=(Solver&& rhs) = default;

  ~Solver();

  void tick(float deltaTime);
  void tickPBD(float deltaTime);
  void tickPD(float deltaTime);

  const std::vector<Vertex>& getVertices() const { return this->_vertices; }

  const std::vector<uint32_t>& getLines() const { return this->_lines; }

  const std::vector<Triangle>& getTriangles() const { return this->_triangles; }

  const SolverOptions& getOptions() const { return this->_options; }

  void clear();

  // Utilities for importing meshes
  void addNodes(const std::vector<glm::vec3>& vertices);
  void addTriMeshVolume(
      const std::vector<glm::vec3>& vertices,
      const std::vector<uint32_t>& triIndices,
      const glm::vec3& initialVelocity,
      float density,
      float strainStiffness,
      float minStrain,
      float maxStrain,
      float volumeStiffness,
      float compression,
      float stretching);
  void addFixedRegions(const std::vector<glm::mat4>& regionMatrices, float w);
  void addLinkedRegions(const std::vector<glm::mat4>& regionsMatrices, float w);

  // Utilities for spawning primitives
  void createBox(const glm::vec3& translation, float scale, float w);
  void createTetBox(
      const glm::vec3& translation,
      float scale,
      const glm::vec3& initialVelocity,
      float w,
      float mass,
      bool hinged);
  void
  createSheet(const glm::vec3& translation, float scale, float mass, float k);
  void createShapeMatchingBox(
      const glm::vec3& translation,
      uint32_t countX,
      uint32_t countY,
      uint32_t countZ,
      float scale,
      const glm::vec3& initialVelocity,
      float w);
  void createShapeMatchingSheet(
      const glm::vec3& translation,
      float scale,
      const glm::vec3& initialVelocity,
      float w);
  void createBendSheet(const glm::vec3& translation, float scale, float w);

private:
  void _parallelComputeCollisions();
  void _parallelPointTriangleCollisions();

  struct NodeCompRange {
    SpatialHashGridCellRange
    operator()(const Node& node, const SpatialHashGrid& grid) const;
  };

  struct TetCompRange {
    const std::vector<Node>& nodes;

    SpatialHashGridCellRange
    operator()(const Tetrahedron& node, const SpatialHashGrid& grid) const;
  };

  struct TriCompRange {
    const std::vector<Node>& nodes;

    SpatialHashGridCellRange
    operator()(const Triangle& triangle, const SpatialHashGrid& grid) const;
  };

  SolverOptions _options;
  uint32_t _constraintId = 0;

  SpatialHash<Node, NodeCompRange> _spatialHashNodes;
  SpatialHash<Triangle, TriCompRange> _spatialHashTris;
  SpatialHash<Tetrahedron, TetCompRange> _spatialHashTets;

  std::vector<Node> _nodes;
  std::vector<Tetrahedron> _tets;
  std::vector<PositionConstraint> _positionConstraints;
  std::vector<DistanceConstraint> _distanceConstraints;
  std::vector<TetrahedralConstraint> _tetConstraints;
  std::vector<VolumeConstraint> _volumeConstraints;
  std::vector<ShapeMatchingConstraint> _shapeConstraints;
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
    std::vector<PointTriangleCollisionConstraint> triCollisions;
    std::vector<EdgeCollisionConstraint> edgeCollisions;
    std::vector<StaticCollisionConstraint> staticCollisions;
  };

  std::vector<ThreadData> _threadData;

  std::vector<CollisionConstraint> _collisions;
  std::vector<PointTriangleCollisionConstraint> _triCollisions;
  std::vector<EdgeCollisionConstraint> _edgeCollisions;
  std::vector<StaticCollisionConstraint> _staticCollisions;

  std::thread _clearSpatialHashThread;

  // TODO: Seperate into individual buffers for each object??
  std::vector<Triangle> _triangles;
  std::vector<uint32_t> _lines;
  std::vector<Vertex> _vertices;
};
} // namespace Pies