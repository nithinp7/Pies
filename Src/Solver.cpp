#include "Solver.h"

#include "CollisionDetection.h"

#include <Eigen/IterativeLinearSolvers>
#include <glm/gtc/matrix_transform.hpp>
#include <omp.h>

#include <cstdint>

namespace Pies {

Solver::Solver(const SolverOptions& options)
    : _options(options),
      _spatialHashNodes(options.gridSpacing),
      _spatialHashTris(options.gridSpacing) {
  // TODO: This is temporary, fix
  assert(options.threadCount == 16);

  Eigen::initParallel();
  omp_set_num_threads(options.threadCount);
  Eigen::setNbThreads(options.threadCount);

  this->_threadData.resize(options.threadCount);
}

Solver::~Solver() {
  if (this->_clearSpatialHashThread.joinable()) {
    this->_clearSpatialHashThread.join();
  }
}

void Solver::tick(float timestep) {
  if (this->_simFailed) {
    return;
  }

  switch (this->_options.solver) {
  case SolverName::PD:
    this->tickPD(timestep);
    break;
  };
}

void Solver::tickPD(float /*timestep*/) {
  uint32_t nodeCount = static_cast<uint32_t>(this->_nodes.size());

  float h = this->_options.fixedTimestepSize / this->_options.timeSubsteps;
  float h2 = h * h;

  if (this->_previousNodeCount != nodeCount) {
    this->_previousNodeCount = nodeCount;

    // TODO: Make smart factorization updates instead of re-running the
    // precomputation
    // Construct stiffness matrix
    this->_stiffnessMatrix = Eigen::SparseMatrix<float>(nodeCount, nodeCount);
    this->_collisionMatrix = Eigen::SparseMatrix<float>(nodeCount, nodeCount);
    this->_stiffnessAndCollisionMatrix =
        Eigen::SparseMatrix<float>(nodeCount, nodeCount);

    for (uint32_t i = 0; i < nodeCount; ++i) {
      this->_stiffnessMatrix.coeffRef(i, i) =
          1.0f / (this->_nodes[i].invMass * h2);
    }

    for (PositionConstraint& constraint : this->_positionConstraints) {
      constraint.setupGlobalStiffnessMatrix(this->_stiffnessMatrix);
    }

    for (DistanceConstraint& constraint : this->_distanceConstraints) {
      constraint.setupGlobalStiffnessMatrix(this->_stiffnessMatrix);
    }

    for (TetrahedralConstraint& constraint : this->_tetConstraints) {
      constraint.setupGlobalStiffnessMatrix(this->_stiffnessMatrix);
    }

    for (VolumeConstraint& constraint : this->_volumeConstraints) {
      constraint.setupGlobalStiffnessMatrix(this->_stiffnessMatrix);
    }

    for (ShapeMatchingConstraint& constraint : this->_shapeConstraints) {
      constraint.setupGlobalStiffnessMatrix(this->_stiffnessMatrix);
    }

    for (GoalMatchingConstraint& constraint : this->_goalConstraints) {
      constraint.setupGlobalStiffnessMatrix(this->_stiffnessMatrix);
    }

    for (BendConstraint& constraint : this->_bendConstraints) {
      constraint.setupGlobalStiffnessMatrix(this->_stiffnessMatrix);
    }

    // Perform Sparse Cholesky LLT factorization
    this->_pLltDecomp =
        std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>>(
            this->_stiffnessMatrix);

    // Set the correct allocations for the force and state vectors
    this->_stateVectorX = Eigen::VectorXf(nodeCount);
    this->_stateVectorY = Eigen::VectorXf(nodeCount);
    this->_stateVectorZ = Eigen::VectorXf(nodeCount);

    this->_forceVector = Eigen::MatrixXf(nodeCount, 3);
    this->_Msn_h2 = Eigen::MatrixXf(nodeCount, 3);
  }

  // Apply global forces
  for (Node& node : this->_nodes) {
    node.force = glm::vec3(0.0f, -this->_options.gravity, 0.0f) / node.invMass;
  }

  for (uint32_t substep = 0; substep < this->_options.timeSubsteps; ++substep) {
    for (uint32_t i = 0; i < nodeCount; ++i) {
      Node& node = this->_nodes[i];
      // Construct momentum estimate for qn+1
      node.position += h * node.velocity;

      glm::vec3 Msn_h2 = node.position / node.invMass / h2;
      this->_Msn_h2.coeffRef(i, 0) = Msn_h2.x;
      this->_Msn_h2.coeffRef(i, 1) = Msn_h2.y;
      this->_Msn_h2.coeffRef(i, 2) = Msn_h2.z;
    }

    this->_parallelPointTriangleCollisions();

    // this->_collisionMatrix.setZero();
    this->_stiffnessAndCollisionMatrix.setZero();

    // TODO: Change this to be a member variable
    static std::vector<Eigen::Triplet<float>> collisionTriplets;
    collisionTriplets.clear();

    for (const PointTriangleCollisionConstraint& collision :
         this->_triCollisions) {
      collision.setupTriplets(collisionTriplets);
      // collision.setupCollisionMatrix(this->_collisionMatrix);
    }

    // for (const EdgeCollisionConstraint& collision : this->_edgeCollisions) {
    //   collision.setupCollisionMatrix(this->_collisionMatrix);
    // }

    for (const StaticCollisionConstraint& collision : this->_staticCollisions) {
      collision.setupTriplets(collisionTriplets);
      // collision.setupCollisionMatrix(this->_collisionMatrix);
    }

    this->_collisionMatrix.setFromTriplets(
        collisionTriplets.begin(),
        collisionTriplets.end());

    this->_stiffnessAndCollisionMatrix =
        this->_stiffnessMatrix + this->_collisionMatrix;

    for (uint32_t iter = 0; iter < this->_options.iterations; ++iter) {
      // Construct global force vector and set initial node position estimate
      this->_forceVector = this->_Msn_h2;

      // Project all constraints onto auxiliary variable (local step)

#pragma omp parallel
      {
#pragma omp single nowait
        {
          if (iter == 0) {
            // TODO: Use CHOLMOD rank-updates??

            // We hide the latency of the below computation behind the first PD
            // iteration since the decomposition result is not needed until the
            // end of the first iteration.
            this->_pLltDecomp = std::make_unique<
                Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>>(
                this->_stiffnessAndCollisionMatrix);
          }
        }

#pragma omp for
        for (int i = 0; i < this->_volumeConstraints.size(); ++i) {
          this->_volumeConstraints[i].projectToAuxiliaryVariable(this->_nodes);
        }

#pragma omp for
        for (int i = 0; i < this->_tetConstraints.size(); ++i) {
          this->_tetConstraints[i].projectToAuxiliaryVariable(this->_nodes);
        }

// Parallelize the rest of the constraints??
#pragma omp single
        {
          for (PositionConstraint& constraint : this->_positionConstraints) {
            constraint.projectToAuxiliaryVariable(this->_nodes);
          }

          for (DistanceConstraint& constraint : this->_distanceConstraints) {
            constraint.projectToAuxiliaryVariable(this->_nodes);
          }

          for (BendConstraint& constraint : this->_bendConstraints) {
            constraint.projectToAuxiliaryVariable(this->_nodes);
          }

          for (ShapeMatchingConstraint& constraint : this->_shapeConstraints) {
            constraint.projectToAuxiliaryVariable(this->_nodes);
          }

          for (GoalMatchingConstraint& constraint : this->_goalConstraints) {
            constraint.projectToAuxiliaryVariable(this->_nodes);
          }

          for (PointTriangleCollisionConstraint& collision :
               this->_triCollisions) {
            collision.projectToAuxiliaryVariable(this->_nodes);
          }

          for (EdgeCollisionConstraint& collision : this->_edgeCollisions) {
            collision.projectToAuxiliaryVariable(this->_nodes);
          }

          for (StaticCollisionConstraint& collision : this->_staticCollisions) {
            collision.projectToAuxiliaryVariable(this->_nodes);
          }
        }

#pragma omp single
        {
          for (PositionConstraint& constraint : this->_positionConstraints) {
            constraint.setupGlobalForceVector(this->_forceVector);
          }

          for (DistanceConstraint& constraint : this->_distanceConstraints) {
            constraint.setupGlobalForceVector(this->_forceVector);
          }

          for (TetrahedralConstraint& constraint : this->_tetConstraints) {
            constraint.setupGlobalForceVector(this->_forceVector);
          }

          for (VolumeConstraint& constraint : this->_volumeConstraints) {
            constraint.setupGlobalForceVector(this->_forceVector);
          }
          for (BendConstraint& constraint : this->_bendConstraints) {
            constraint.setupGlobalForceVector(this->_forceVector);
          }

          for (ShapeMatchingConstraint& constraint : this->_shapeConstraints) {
            constraint.setupGlobalForceVector(this->_forceVector);
          }

          for (GoalMatchingConstraint& constraint : this->_goalConstraints) {
            constraint.setupGlobalForceVector(this->_forceVector);
          }

          for (const PointTriangleCollisionConstraint& collision :
               this->_triCollisions) {
            collision.setupGlobalForceVector(this->_forceVector);
          }

          for (const EdgeCollisionConstraint& collision :
               this->_edgeCollisions) {
            collision.setupGlobalForceVector(this->_forceVector);
          }

          for (const StaticCollisionConstraint& collision :
               this->_staticCollisions) {
            collision.setupGlobalForceVector(this->_forceVector);
          }
        }

#pragma omp barrier

        // Solve x,y,z coordinates in parallel
#pragma omp sections
        {
#pragma omp section
          {
            this->_stateVectorX =
                this->_pLltDecomp->solve(this->_forceVector.col(0));
          }
#pragma omp section
          {
            this->_stateVectorY =
                this->_pLltDecomp->solve(this->_forceVector.col(1));
          }
#pragma omp section
          {
            this->_stateVectorZ =
                this->_pLltDecomp->solve(this->_forceVector.col(2));
          }
        }

        // Update node positions from the global solve results
#pragma omp for
        for (int i = 0; i < nodeCount; ++i) {
          Node& node = this->_nodes[i];
          node.position.x = this->_stateVectorX.coeff(i);
          node.position.y = this->_stateVectorY.coeff(i);
          node.position.z = this->_stateVectorZ.coeff(i);
        }
      }
    }

    for (uint32_t collisionIter = 0;
         collisionIter < this->_options.collisionStabilizationIterations;
         ++collisionIter) {

      for (PointTriangleCollisionConstraint& collision : this->_triCollisions) {
        collision.stabilizeCollisions(this->_nodes);
      }

      for (EdgeCollisionConstraint& collision : this->_edgeCollisions) {
        collision.stabilizeCollisions(this->_nodes);
      }

      for (StaticCollisionConstraint& collision : this->_staticCollisions) {
        // TODO: stiffness
        this->_nodes[collision.nodeId].position = collision.projectedPosition;
      }
    }

    // Update node velocities
    for (uint32_t i = 0; i < nodeCount; ++i) {
      Node& node = this->_nodes[i];
      node.velocity = (1.0f - this->_options.damping) *
                          (node.position - node.prevPosition) / h +
                      h * node.force * node.invMass;

      node.prevPosition = node.position;
      this->_vertices[i].position = node.position;
      // this->_vertices[i].baseColor = glm::vec3(0.0f, 1.0f, 0.0f);
    }

    // Update friction
    for (const PointTriangleCollisionConstraint& collision :
         this->_triCollisions) {
      // if (!collision.colliding) {
      //   continue;
      // }

      // for (uint32_t i = 0; i < 4; ++i) {
      //   this->_vertices[collision.nodeIds[i]].baseColor =
      //       glm::vec3(1.0f, 0.0f, 0.0f);
      // }

      Node& a = this->_nodes[collision.nodeIds[0]];
      Node& b = this->_nodes[collision.nodeIds[1]];
      Node& c = this->_nodes[collision.nodeIds[2]];
      Node& d = this->_nodes[collision.nodeIds[3]];

      glm::vec3 avgTriVelocity = (b.velocity + c.velocity + d.velocity) / 3.0f;

      glm::vec3 n = glm::normalize(
          glm::cross(c.position - b.position, d.position - b.position));

      glm::vec3 relativeVelocity = a.velocity - avgTriVelocity;
      float vDotN = glm::dot(relativeVelocity, n);
      glm::vec3 normVel = vDotN * n;
      glm::vec3 perpVel = relativeVelocity - normVel;

      float friction = this->_options.friction;
      if (glm::length(perpVel) < this->_options.staticFrictionThreshold) {
        friction = 1.0f;
      }

      float triWSum = b.invMass + c.invMass + d.invMass;
      float wSum = a.invMass + triWSum;

      glm::vec3 dv = -friction * perpVel - 1.1f * glm::min(vDotN, 0.0f) * n;

      a.velocity += dv * a.invMass / wSum;
      b.velocity += -dv * triWSum / wSum;
      c.velocity += -dv * triWSum / wSum;
      d.velocity += -dv * triWSum / wSum;
    }

    for (const StaticCollisionConstraint& collision : this->_staticCollisions) {
      Node& node = this->_nodes[collision.nodeId];

      glm::vec3 perpVel = glm::vec3(node.velocity.x, 0.0f, node.velocity.z);

      float friction = this->_options.friction;
      if (glm::length(perpVel) < this->_options.staticFrictionThreshold) {
        friction = 1.0f;
      }

      node.velocity += -friction * perpVel;
    }
  }
} // namespace Pies

void Solver::clear() {
  this->_lines.clear();
  this->_triangles.clear();
  this->_nodes.clear();
  this->_distanceConstraints.clear();
  this->_tetConstraints.clear();
  this->_volumeConstraints.clear();
  this->_shapeConstraints.clear();
  this->_goalConstraints.clear();
  this->_bendConstraints.clear();
  this->_tets.clear();
  this->_positionConstraints.clear();
  this->_vertices.clear();
  this->_constraintId = 0;

  this->_stateVectorX = {};
  this->_stateVectorY = {};
  this->_stateVectorZ = {};
  this->_stiffnessMatrix = {};

  this->renderStateDirty = true;
}

namespace {
SpatialHashGridCellRange sweptTriRange(
    const Triangle& triangle,
    const std::vector<Node>& nodes,
    const SpatialHashGrid& grid,
    float threshold,
    bool& simFailed) {
  const Node& n1 = nodes[triangle.nodeIds[0]];
  const Node& n2 = nodes[triangle.nodeIds[1]];
  const Node& n3 = nodes[triangle.nodeIds[2]];

  glm::vec3 max = n1.position;
  glm::vec3 min = n1.position;

  for (uint32_t i = 0; i < 3; ++i) {
    const Node& node = nodes[triangle.nodeIds[i]];
    max = glm::max(node.position, max);
    max = glm::max(node.prevPosition, max);

    min = glm::min(node.position, min);
    min = glm::min(node.prevPosition, min);
  }

  glm::vec3 padding(0.5f * threshold);
  min -= padding;
  max += padding;

  min /= grid.scale;
  max /= grid.scale;

  SpatialHashGridCellRange range{};
  range.minX = static_cast<int64_t>(glm::floor(min.x));
  range.minY = static_cast<int64_t>(glm::floor(min.y));
  range.minZ = static_cast<int64_t>(glm::floor(min.z));

  range.lengthX = static_cast<uint32_t>(glm::ceil(max.x) - range.minX);
  range.lengthY = static_cast<uint32_t>(glm::ceil(max.y) - range.minY);
  range.lengthZ = static_cast<uint32_t>(glm::ceil(max.z) - range.minZ);

  if (range.lengthX > 60 || range.lengthY > 60 || range.lengthZ > 60) {
    simFailed = true;
    return {};
  }

  return range;
}

bool trianglesFacingIn(
    const Triangle& tri1,
    const Triangle& tri2,
    const std::vector<Node>& nodes) {
  const glm::vec3& a = nodes[tri1.nodeIds[0]].position;
  const glm::vec3& b = nodes[tri1.nodeIds[1]].position;
  const glm::vec3& c = nodes[tri1.nodeIds[2]].position;

  const glm::vec3& d = nodes[tri2.nodeIds[0]].position;
  const glm::vec3& e = nodes[tri2.nodeIds[1]].position;
  const glm::vec3& f = nodes[tri2.nodeIds[2]].position;

  glm::vec3 n1 = glm::cross(b - a, c - a);
  glm::vec3 n2 = glm::cross(e - d, f - d);

  // Have further limited cone?
  return glm::dot(n1, n2) < 0.0f;
}
} // namespace

void Solver::_parallelPointTriangleCollisions() {
  // Update collisions
  this->_triCollisions.clear();
  this->_edgeCollisions.clear();
  this->_staticCollisions.clear();
  this->_collisionMatrix.setZero();
  this->_stiffnessAndCollisionMatrix.setZero();

  if (this->_clearSpatialHashThread.joinable()) {
    this->_clearSpatialHashThread.join();
  }

  this->_spatialHashTris.parallelBulkInsert(
      this->_triangles,
      {this->_nodes, this->_options.collisionThresholdDistance});

  // Detect and resolve collisions
  auto fnComputeCollisions =
      [&nodes = this->_nodes,
       &triangles = this->_triangles,
       threadCount = this->_options.threadCount,
       &spatialHash = this->_spatialHashTris,
       &threadData = this->_threadData,
       threshold = this->_options.collisionThresholdDistance,
       thickness = this->_options.collisionThickness,
       floorHeight = this->_options.floorHeight](size_t threadId) {
        ThreadData& data = threadData[threadId];

        std::vector<const SpatialHashGridCellBucket<Triangle>*> buckets;

        data.triCollisions.clear();
        data.edgeCollisions.clear();
        data.staticCollisions.clear();

        for (size_t triId = threadId; triId < triangles.size();
             triId += threadCount) {
          const Triangle& tri = triangles[triId];

          SpatialHashGridCellRange range = sweptTriRange(
              tri,
              nodes,
              spatialHash.getGrid(),
              threshold,
              data.failed);
          if (data.failed) {
            return;
          }

          //     ccdRange(nodeA.prevPosition, nodeA.position,
          //     spatialHash.getGrid());

          for (uint32_t dx = 0; dx < range.lengthX; ++dx) {
            for (uint32_t dy = 0; dy < range.lengthY; ++dy) {
              for (uint32_t dz = 0; dz < range.lengthZ; ++dz) {
                // Compute grid cell id
                SpatialHashGridCellId id{
                    range.minX + dx,
                    range.minY + dy,
                    range.minZ + dz};

                const SpatialHashGridCellBucket<Triangle>* pBucket =
                    spatialHash.findCollisions(id);
                if (pBucket) {
                  buckets.push_back(pBucket);
                }
              }
            }
          }

          if (buckets.size() > 200) {
            // Safety check to avoid simulation hangs
            data.failed = true;
            return;
          }

          for (const SpatialHashGridCellBucket<Triangle>* pBucket : buckets) {
            // Check all triangles in the bucket
            for (const Triangle* pOtherTri : pBucket->values) {

              if (pBucket->values.size() > 500) {
                // Safety check to avoid simulation hangs
                data.failed = true;
                return;
              }

              bool containsCommonNode = false;
              // TODO: Is there any case where we want to collide between two
              // connected triangles??
              for (uint32_t i = 0; i < 3; ++i) {
                for (uint32_t j = 0; j < 3; ++j) {
                  if (tri.nodeIds[i] == pOtherTri->nodeIds[j]) {
                    containsCommonNode = true;
                  }
                }
              }

              if (containsCommonNode) {
                continue;
              }

              if (!trianglesFacingIn(tri, *pOtherTri, nodes)) {
                continue;
              }

              const Node& nodeB = nodes[pOtherTri->nodeIds[0]];
              const Node& nodeC = nodes[pOtherTri->nodeIds[1]];
              const Node& nodeD = nodes[pOtherTri->nodeIds[2]];

              // Point-triangle collisions
              for (uint32_t i = 0; i < 3; ++i) {
                const Node& nodeA = nodes[tri.nodeIds[i]];

                std::optional<float> optT =
                    CollisionDetection::pointTriangleCCD(
                        nodeA.prevPosition - nodeB.prevPosition,
                        nodeC.prevPosition - nodeB.prevPosition,
                        nodeD.prevPosition - nodeB.prevPosition,
                        nodeA.position - nodeB.position,
                        nodeC.position - nodeB.position,
                        nodeD.position - nodeB.position,
                        threshold);

                if (!optT) {
                  // CCD did not find intersection
                  continue;
                }

                data.triCollisions
                    .emplace_back(nodeA, nodeB, nodeC, nodeD, thickness, *optT);
              }

              // for (uint32_t i = 0; i < 3; ++i) {
              //   for (uint32_t j = 0; j < 3; ++j) {
              //     const Node& a = nodes[tri.nodeIds[i]];
              //     const Node& b = nodes[tri.nodeIds[(i + 1) % 3]];

              //     const Node& c = nodes[pOtherTri->nodeIds[j]];
              //     const Node& d = nodes[pOtherTri->nodeIds[(j + 1) % 3]];

              //     std::optional<float> optT =
              //     CollisionDetection::edgeEdgeCCD(
              //         b.prevPosition - a.prevPosition,
              //         c.prevPosition - a.prevPosition,
              //         d.prevPosition - a.prevPosition,
              //         b.position - a.position,
              //         c.position - a.position,
              //         d.position - a.position);

              //     if (!optT) {
              //       // CCD did not find intersection
              //       continue;
              //     }

              //     data.edgeCollisions.emplace_back(a, b, c, d);
              //   }
              // }
            }
          }

          buckets.clear();

          for (uint32_t i = 0; i < 3; ++i) {
            const Node& node = nodes[tri.nodeIds[i]];
            if (node.position.y < floorHeight + thickness) {
              data.staticCollisions.emplace_back(node);
            }
          }
        }
      };

// TODO: No need for lambda function fnComputeCollisions
#pragma omp parallel
  {
    size_t thread_idx = static_cast<size_t>(omp_get_thread_num());
    fnComputeCollisions(thread_idx);
  }

  // TODO: Use OpenMP instead of std thread here as well??
  this->_clearSpatialHashThread =
      std::thread([this]() { this->_spatialHashTris.clear(); });

  // Aggregate across per-thread results
  for (const ThreadData& data : this->_threadData) {
    if (data.failed) {
      this->_simFailed = true;
      return;
    }

    for (const PointTriangleCollisionConstraint& collision :
         data.triCollisions) {
      // collision.setupGlobalForceVector(this->_forceVector);
      // collision.setupCollisionMatrix(this->_collisionMatrix);
      this->_triCollisions.push_back(collision);
    }

    for (const EdgeCollisionConstraint& collision : data.edgeCollisions) {
      this->_edgeCollisions.push_back(collision);
    }

    for (const StaticCollisionConstraint& collision : data.staticCollisions) {
      // collision.setupCollisionMatrix(this->_collisionMatrix);
      // collision.setupGlobalForceVector(this->_forceVector);
      this->_staticCollisions.push_back(collision);
    }
  }
}

SpatialHashGridCellRange Solver::NodeCompRange::operator()(
    const Node& node,
    const SpatialHashGrid& grid) const {
  float radiusPadding = 0.5f;
  float gridLocalRadius = (node.radius + radiusPadding) / grid.scale;
  glm::vec3 gridLocalPos = node.position / grid.scale;
  glm::vec3 gridLocalMin = gridLocalPos - glm::vec3(gridLocalRadius);

  SpatialHashGridCellRange range{};
  range.minX = static_cast<int64_t>(glm::floor(gridLocalMin.x));
  range.minY = static_cast<int64_t>(glm::floor(gridLocalMin.y));
  range.minZ = static_cast<int64_t>(glm::floor(gridLocalMin.z));

  glm::vec3 adjustedDiameter =
      glm::ceil(glm::fract(gridLocalMin) + glm::vec3(2 * gridLocalRadius));
  range.lengthX = static_cast<uint32_t>(adjustedDiameter.x);
  range.lengthY = static_cast<uint32_t>(adjustedDiameter.y);
  range.lengthZ = static_cast<uint32_t>(adjustedDiameter.z);

  if (range.lengthX > 50 || range.lengthY > 50 || range.lengthZ > 50) {
    return {};
  }

  return range;
}

SpatialHashGridCellRange Solver::TriCompRange::operator()(
    const Triangle& triangle,
    const SpatialHashGrid& grid) const {
  const Node& n1 = nodes[triangle.nodeIds[0]];
  const Node& n2 = nodes[triangle.nodeIds[1]];
  const Node& n3 = nodes[triangle.nodeIds[2]];

  glm::vec3 max = n1.position;
  glm::vec3 min = n1.position;

  for (uint32_t i = 0; i < 3; ++i) {
    const Node& node = nodes[triangle.nodeIds[i]];
    max = glm::max(node.position, max);
    max = glm::max(node.prevPosition, max);

    min = glm::min(node.position, min);
    min = glm::min(node.prevPosition, min);
  }

  glm::vec3 padding(0.5f * this->threshold);

  min -= padding;
  max += padding;

  min /= grid.scale;
  max /= grid.scale;

  SpatialHashGridCellRange range{};
  range.minX = static_cast<int64_t>(glm::floor(min.x));
  range.minY = static_cast<int64_t>(glm::floor(min.y));
  range.minZ = static_cast<int64_t>(glm::floor(min.z));

  range.lengthX = static_cast<uint32_t>(glm::ceil(max.x) - range.minX);
  range.lengthY = static_cast<uint32_t>(glm::ceil(max.y) - range.minY);
  range.lengthZ = static_cast<uint32_t>(glm::ceil(max.z) - range.minZ);

  if (range.lengthX > 50 || range.lengthY > 50 || range.lengthZ > 50) {
    return {};
  }

  return range;
}
} // namespace Pies