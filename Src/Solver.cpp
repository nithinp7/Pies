#include "Solver.h"

#include "CollisionDetection.h"

#include <glm/gtc/matrix_transform.hpp>

#include <cstdint>

namespace Pies {

Solver::Solver(const SolverOptions& options)
    : _options(options),
      _spatialHashNodes(options.gridSpacing),
      _spatialHashTets(options.gridSpacing),
      _spatialHashTris(options.gridSpacing) {
  this->_threadData.resize(this->_options.threadCount);
}

Solver::~Solver() {
  if (this->_clearSpatialHashThread.joinable()) {
    this->_clearSpatialHashThread.join();
  }
}

void Solver::tick(float timestep) {
  switch (this->_options.solver) {
  case SolverName::PBD:
    this->tickPBD(timestep);
    break;
  case SolverName::PD:
    this->tickPD(timestep);
    break;
  };
}

void Solver::tickPBD(float /*timestep*/) {
  float deltaTime =
      this->_options.fixedTimestepSize / this->_options.timeSubsteps;

  // Time substeps
  for (uint32_t substep = 0; substep < this->_options.timeSubsteps; ++substep) {
    // Apply external forces and advect nodes
    for (Node& node : this->_nodes) {
      node.prevPosition = node.position;
      node.position += node.velocity * deltaTime +
                       glm::vec3(0.0f, -this->_options.gravity, 0.0f) *
                           deltaTime * deltaTime;
    }

    // TODO: Collision solver
    // this->_spatialHashNodes.clear();
    // this->_spatialHashNodes.parallelBulkInsert(this->_nodes, {});

    for (uint32_t i = 0; i < this->_options.iterations; ++i) {
      if (!releaseHinge) {
        for (PositionConstraint& constraint : this->_positionConstraints) {
          constraint.projectNodePositions(this->_nodes);
        }
      }

      for (DistanceConstraint& constraint : this->_distanceConstraints) {
        constraint.projectNodePositions(this->_nodes);
      }

      for (TetrahedralConstraint& constraint : this->_tetConstraints) {
        constraint.projectNodePositions(this->_nodes);
      }

      // // TODO: Collision solver
      // this->_spatialHashTets.clear();
      // this->_spatialHashTets.parallelBulkInsert(this->_tets, {this->_nodes});

      this->_spatialHashNodes.clear();
      this->_spatialHashNodes.parallelBulkInsert(this->_nodes, {});

      // Detect and resolve collisions
      std::vector<const SpatialHashGridCellBucket<Node>*> scratchBuckets;
      for (Node& node : this->_nodes) {
        this->_spatialHashNodes.findCollisions(node, {}, scratchBuckets);
        for (const SpatialHashGridCellBucket<Node>* pBucket : scratchBuckets) {
          // Check each node within each bucket
          for (Node* pOtherNode : pBucket->values) {
            // Check intersection
            glm::vec3 diff = pOtherNode->position - node.position;
            float dist = glm::length(diff);

            float disp = node.radius + pOtherNode->radius - dist;

            if (disp <= 0.0) {
              continue;
            }

            glm::vec3 dir(1.0f, 0.0f, 0.0f);
            if (dist > 0.00001f) {
              dir = diff / dist;
            }

            float wSum = node.invMass + pOtherNode->invMass;

            node.position += 0.85f * -disp * dir * node.invMass / wSum;
            pOtherNode->position +=
                0.85f * disp * dir * pOtherNode->invMass / wSum;

            // TODO: Add friction between dynamic objects
            glm::vec3 relativeVelocity = pOtherNode->velocity - node.velocity;
            glm::vec3 perpVel =
                relativeVelocity - glm::dot(relativeVelocity, dir) * dir;

            float friction = this->_options.friction;
            if (glm::length(perpVel) < this->_options.staticFrictionThreshold) {
              friction = 1.0f;
            }

            // TODO: Decouple friction from solver iteration count
            node.velocity += -friction * perpVel * node.invMass / wSum;
            pOtherNode->velocity +=
                friction * perpVel * pOtherNode->invMass / wSum;
          }
        }

        scratchBuckets.clear();
      }

      for (Node& node : this->_nodes) {
        if (node.position.y - node.radius < this->_options.floorHeight) {
          node.position.y = this->_options.floorHeight + node.radius;
        }
      }
    }

    // Compute new velocity and construct new vertex positions
    for (uint32_t i = 0; i < this->_nodes.size(); ++i) {
      Node& node = this->_nodes[i];

      node.velocity = (1.0f - this->_options.damping) *
                      (node.position - node.prevPosition) / deltaTime;

      // TODO: friction for dynamically generated constraints
      if (node.position.y - node.radius <= this->_options.floorHeight) {
        if (glm::length(glm::vec2(node.velocity.x, node.velocity.z)) < 5.0f) {
          node.velocity.x = 0.0f;
          node.velocity.z = 0.0f;
        } else {
          node.velocity.x *= 1.0f - this->_options.friction;
          node.velocity.z *= 1.0f - this->_options.friction;
        }
      }

      this->_vertices[i].position = this->_nodes[i].position;
    }
  }
}

void Solver::tickPD(float /*timestep*/) {
  uint32_t nodeCount = static_cast<uint32_t>(this->_nodes.size());

  // TODO: Remove time substeps for PD solver
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

    // Perform Sparse Cholesky LLT factorization
    this->_pLltDecomp =
        std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>>(
            this->_stiffnessMatrix);

    // Set the correct allocations for the force and state vectors
    this->_stateVector = Eigen::MatrixXf(nodeCount, 3);
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
    }

    this->_parallelPointTriangleCollisions();
#if 1
    for (uint32_t collisionIter = 0;
         collisionIter < this->_options.collisionIterations;
         ++collisionIter) {
      // this->_parallelComputeCollisions();
      // this->_parallelPointTriangleCollisions();
      for (CollisionConstraint& collision : this->_collisions) {
        collision.projectToAuxiliaryVariable(this->_nodes);
        Node& nodeA = this->_nodes[collision.nodeIds[0]];
        Node& nodeB = this->_nodes[collision.nodeIds[1]];
        nodeA.position += this->_options.collionStiffness *
                          (collision.projectedPositions[0] - nodeA.position);
        nodeB.position += this->_options.collionStiffness *
                          (collision.projectedPositions[1] - nodeB.position);
      }

      for (PointTriangleCollisionConstraint& collision : this->_triCollisions) {
        collision.stabilizeCollisions(this->_nodes);
      }

      for (StaticCollisionConstraint& collision : this->_staticCollisions) {
        // TODO: stiffness
        this->_nodes[collision.nodeId].position = collision.projectedPosition;
      }
    }
#endif

    for (uint32_t i = 0; i < nodeCount; ++i) {
      const Node& node = this->_nodes[i];
      glm::vec3 Msn_h2 = node.position / node.invMass / h2;
      this->_Msn_h2.coeffRef(i, 0) = Msn_h2.x;
      this->_Msn_h2.coeffRef(i, 1) = Msn_h2.y;
      this->_Msn_h2.coeffRef(i, 2) = Msn_h2.z;
    }

  
    this->_collisionMatrix.setZero();
    this->_stiffnessAndCollisionMatrix.setZero();

    for (const PointTriangleCollisionConstraint& collision :
          this->_triCollisions) {
      collision.setupCollisionMatrix(this->_collisionMatrix);
    }

    for (const StaticCollisionConstraint& collision :
          this->_staticCollisions) {
      collision.setupCollisionMatrix(this->_collisionMatrix);
    }

    this->_stiffnessAndCollisionMatrix =
        this->_stiffnessMatrix + this->_collisionMatrix;
    this->_pLltDecomp =
        std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>>(
            this->_stiffnessAndCollisionMatrix);


    for (uint32_t iter = 0; iter < this->_options.iterations; ++iter) {
      // Construct global force vector and set initial node position estimate
      this->_forceVector = this->_Msn_h2;

      // Project all constraints onto auxiliary variable (local step)
      // TODO: Parallelize this step
      for (PositionConstraint& constraint : this->_positionConstraints) {
        constraint.projectToAuxiliaryVariable(this->_nodes);
      }

      for (DistanceConstraint& constraint : this->_distanceConstraints) {
        constraint.projectToAuxiliaryVariable(this->_nodes);
      }

      for (TetrahedralConstraint& constraint : this->_tetConstraints) {
        constraint.projectToAuxiliaryVariable(this->_nodes);
      }

      for (VolumeConstraint& constraint : this->_volumeConstraints) {
        constraint.projectToAuxiliaryVariable(this->_nodes);
      }

      for (ShapeMatchingConstraint& constraint : this->_shapeConstraints) {
        constraint.projectToAuxiliaryVariable(this->_nodes);
      }

      for (PointTriangleCollisionConstraint& collision : this->_triCollisions) {
        collision.projectToAuxiliaryVariable(this->_nodes);
      }

      for (StaticCollisionConstraint& collision : this->_staticCollisions) {
        const Node& node = this->_nodes[collision.nodeId];
        collision.projectedPosition = node.position;
        collision.projectedPosition.y = this->_options.floorHeight;
      }

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

      for (ShapeMatchingConstraint& constraint : this->_shapeConstraints) {
        constraint.setupGlobalForceVector(this->_forceVector);
      }

      for (const PointTriangleCollisionConstraint& collision :
           this->_triCollisions) {
        collision.setupGlobalForceVector(this->_forceVector);
      }

      for (const StaticCollisionConstraint& collision :
           this->_staticCollisions) {
        collision.setupGlobalForceVector(this->_forceVector);
      }

      // Parallelize solving x,y,z coordinates
      // Note: Make sure to offset the starting node index of each thread
      // to avoid false-sharing when writing to the state vector

      // Solve the global system using the precomputed factorization
      this->_stateVector = this->_pLltDecomp->solve(this->_forceVector);

      // Update node positions from the global solve results
      for (uint32_t i = 0; i < nodeCount; ++i) {
        Node& node = this->_nodes[i];
        node.position.x = this->_stateVector.coeff(i, 0);
        node.position.y = this->_stateVector.coeff(i, 1);
        node.position.z = this->_stateVector.coeff(i, 2);
      }
    }

    for (uint32_t collisionIter = 0;
         collisionIter < this->_options.collisionStabilizationIterations;
         ++collisionIter) {

      for (PointTriangleCollisionConstraint& collision : this->_triCollisions) {
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
    for (const CollisionConstraint& collision : this->_collisions) {
      Node& a = this->_nodes[collision.nodeIds[0]];
      Node& b = this->_nodes[collision.nodeIds[1]];

      glm::vec3 diff = b.position - a.position;
      float dist = glm::length(diff);

      if (dist > a.radius + b.radius) {
        continue;
      }

      glm::vec3 n = diff / dist;

      glm::vec3 relativeVelocity = b.velocity - a.velocity;
      glm::vec3 perpVel = relativeVelocity - glm::dot(relativeVelocity, n) * n;

      float friction = -this->_options.friction;
      if (glm::length(perpVel) < this->_options.staticFrictionThreshold) {
        friction = 1.0f;
      }

      float wSum = a.invMass + b.invMass;

      a.velocity += -friction * perpVel * a.invMass / wSum;
      b.velocity += friction * perpVel * b.invMass / wSum;

      // this->_vertices[collision.nodeIds[0]].baseColor =
      //     glm::vec3(1.0f, 0.0f, 0.0f);
      // this->_vertices[collision.nodeIds[1]].baseColor =
      //     glm::vec3(1.0f, 0.0f, 0.0f);
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

      glm::vec3 dv = -friction * perpVel - glm::min(vDotN, 0.0f) * n;

      a.velocity += dv * a.invMass / wSum;
      b.velocity += -dv * triWSum / wSum;
      c.velocity += -dv * triWSum / wSum;
      d.velocity += -dv * triWSum / wSum;
    }

    for (const StaticCollisionConstraint& collision : this->_staticCollisions) {
      Node& node = this->_nodes[collision.nodeId];

      glm::vec3 perpVel =
          node.velocity - glm::dot(node.velocity, collision.n) * collision.n;

      float friction = this->_options.friction;
      if (glm::length(perpVel) < this->_options.staticFrictionThreshold) {
        friction = 1.0f;
      }

      node.velocity += -friction * perpVel;
    }
  }
}

void Solver::clear() {
  this->_lines.clear();
  this->_triangles.clear();
  this->_nodes.clear();
  this->_distanceConstraints.clear();
  this->_tetConstraints.clear();
  this->_volumeConstraints.clear();
  this->_shapeConstraints.clear();
  this->_tets.clear();
  this->_positionConstraints.clear();
  this->_vertices.clear();
  this->_constraintId = 0;

  this->_stateVector = {};
  this->_stiffnessMatrix = {};

  this->renderStateDirty = true;
}

void Solver::_computeCollisions() {
  // Update collisions
  this->_collisions.clear();
  this->_triCollisions.clear();
  this->_staticCollisions.clear();
  this->_collisionMatrix.setZero();
  this->_stiffnessAndCollisionMatrix.setZero();

  this->_spatialHashNodes.clear();
  this->_spatialHashNodes.parallelBulkInsert(this->_nodes, {});

  // Detect and resolve collisions
  // TODO: Parallelize this!!
  std::vector<const SpatialHashGridCellBucket<Node>*> scratchBuckets;
  for (Node& node : this->_nodes) {
    this->_spatialHashNodes.findCollisions(node, {}, scratchBuckets);
    for (const SpatialHashGridCellBucket<Node>* pBucket : scratchBuckets) {
      // Check each node within each bucket
      for (Node* pOtherNode : pBucket->values) {
        if (node.id <= pOtherNode->id) {
          // Only fill lower-triangular matrix
          continue;
        }

        this->_collisions.emplace_back(node, *pOtherNode);
      }
    }

    if (node.position.y - node.radius < this->_options.floorHeight) {
      this->_staticCollisions.emplace_back(
          node,
          glm::vec3(
              node.position.x,
              this->_options.floorHeight + node.radius,
              node.position.z));
    }
  }
}

void Solver::_parallelComputeCollisions() {
  // Update collisions
  this->_collisions.clear();
  this->_staticCollisions.clear();
  this->_collisionMatrix.setZero();
  this->_stiffnessAndCollisionMatrix.setZero();

  if (this->_clearSpatialHashThread.joinable()) {
    this->_clearSpatialHashThread.join();
  }

  this->_spatialHashNodes.parallelBulkInsert(this->_nodes, {});

  // Detect and resolve collisions
  auto fnComputeCollisions = [&nodes = this->_nodes,
                              threadCount = this->_options.threadCount,
                              &spatialHash = this->_spatialHashNodes,
                              &threadData = this->_threadData,
                              floorHeight =
                                  this->_options.floorHeight](size_t threadId) {
    std::vector<const SpatialHashGridCellBucket<Node>*> scratchBuckets;
    ThreadData& data = threadData[threadId];

    data.collisions.clear();
    data.staticCollisions.clear();

    for (size_t nodeId = threadId; nodeId < nodes.size();
         nodeId += threadCount) {
      const Node& node = nodes[nodeId];

      spatialHash.findCollisions(node, {}, scratchBuckets);
      for (const SpatialHashGridCellBucket<Node>* pBucket : scratchBuckets) {
        // Check each node within each bucket
        for (Node* pOtherNode : pBucket->values) {
          if (node.id <= pOtherNode->id) {
            // Only fill lower-triangular matrix
            continue;
          }

          data.collisions.emplace_back(node, *pOtherNode);
        }
      }

      if (node.position.y - node.radius < floorHeight) {
        data.staticCollisions.emplace_back(
            node,
            glm::vec3(
                node.position.x,
                floorHeight + node.radius,
                node.position.z));
      }
    }
  };

  std::vector<std::thread> threads;
  for (uint32_t threadId = 0; threadId < this->_options.threadCount;
       ++threadId) {
    threads.emplace_back(fnComputeCollisions, threadId);
  }

  for (std::thread& thread : threads) {
    thread.join();
  }

  this->_clearSpatialHashThread =
      std::thread([this]() { this->_spatialHashNodes.clear(); });

  // Parallelize the setup of dynamic collision matrix and global force vector.
  // We assign each thread to all nodes with the same nodeId % threadCount.
  /*auto fnSetupCollision = [this](uint32_t threadId) {
    for (const ThreadData& data : this->_threadData) {
      for (const CollisionConstraint& collision : data.collisions) {
        collision.setupGlobalForceVector(
            this->_forceVector,
            threadId,
            this->_options.threadCount);
      }
    }
  };

  for (uint32_t threadId = 0; threadId < this->_options.threadCount;
       ++threadId) {
    threads[threadId] = std::thread(fnSetupCollision, threadId);
  }

  for (std::thread& thread : threads) {
    thread.join();
  }*/

  // Aggregate across per-thread results
  for (const ThreadData& data : this->_threadData) {
    for (const CollisionConstraint& collision : data.collisions) {
      // collision.setupGlobalForceVector(this->_forceVector, 0, 1);
      // collision.setupCollisionMatrix(this->_collisionMatrix);
      this->_collisions.push_back(collision);
    }

    for (const StaticCollisionConstraint& collision : data.staticCollisions) {
      // collision.setupCollisionMatrix(this->_collisionMatrix);
      // collision.setupGlobalForceVector(this->_forceVector);
      this->_staticCollisions.push_back(collision);
    }
  }
}

namespace {
SpatialHashGridCellRange ccdRange(
    const glm::vec3& prevPos,
    const glm::vec3& pos,
    const SpatialHashGrid& grid) {
  float radius = 0.5f * glm::length(pos - prevPos);
  glm::vec3 center = 0.5f * pos + 0.5f * prevPos;
  float radiusPadding = 1.0f;
  float gridLocalRadius = (radius + radiusPadding) / grid.scale;
  glm::vec3 gridLocalPos = center / grid.scale;
  glm::vec3 gridLocalMin = gridLocalPos - glm::vec3(gridLocalRadius);

  SpatialHashGridCellRange range{};
  range.minX = static_cast<int64_t>(glm::floor(gridLocalMin.x));
  range.minY = static_cast<int64_t>(glm::floor(gridLocalMin.y));
  range.minZ = static_cast<int64_t>(glm::floor(gridLocalMin.z));

  glm::vec3 adjustedDiameter =
      glm::ceil(glm::fract(gridLocalMin) + glm::vec3(2 * gridLocalRadius));
  range.lengthX = static_cast<uint32_t>(adjustedDiameter.r);
  range.lengthY = static_cast<uint32_t>(adjustedDiameter.g);
  range.lengthZ = static_cast<uint32_t>(adjustedDiameter.b);

  if (range.lengthX > 100 || range.lengthY > 100 || range.lengthZ > 100) {
    return {};
  }

  return range;
}
} // namespace

void Solver::_parallelPointTriangleCollisions() {
  // Update collisions
  this->_collisions.clear();
  this->_triCollisions.clear();
  this->_staticCollisions.clear();
  this->_collisionMatrix.setZero();
  this->_stiffnessAndCollisionMatrix.setZero();

  if (this->_clearSpatialHashThread.joinable()) {
    this->_clearSpatialHashThread.join();
  }

  this->_spatialHashTris.parallelBulkInsert(this->_triangles, {this->_nodes});

  // Detect and resolve collisions
  auto fnComputeCollisions = [&nodes = this->_nodes,
                              threadCount = this->_options.threadCount,
                              &spatialHash = this->_spatialHashTris,
                              &threadData = this->_threadData,
                              floorHeight =
                                  this->_options.floorHeight](size_t threadId) {
    ThreadData& data = threadData[threadId];

    std::vector<const SpatialHashGridCellBucket<Triangle>*> buckets;

    data.collisions.clear();
    data.triCollisions.clear();
    data.staticCollisions.clear();

    for (size_t nodeId = threadId; nodeId < nodes.size();
         nodeId += threadCount) {
      const Node& nodeA = nodes[nodeId];

      SpatialHashGridCellRange range =
          ccdRange(nodeA.prevPosition, nodeA.position, spatialHash.getGrid());

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

      for (const SpatialHashGridCellBucket<Triangle>* pBucket : buckets) {
        // Check all triangles in the bucket
        for (const Triangle* pTriangle : pBucket->values) {
          if (nodeId == pTriangle->nodeIds[0] ||
              nodeId == pTriangle->nodeIds[1] ||
              nodeId == pTriangle->nodeIds[2]) {
            continue;
          }

          const Node& nodeB = nodes[pTriangle->nodeIds[0]];
          const Node& nodeC = nodes[pTriangle->nodeIds[1]];
          const Node& nodeD = nodes[pTriangle->nodeIds[2]];

          std::optional<float> optT = CollisionDetection::linearCCD(
              nodeA.prevPosition - nodeB.prevPosition,
              nodeC.prevPosition - nodeB.prevPosition,
              nodeD.prevPosition - nodeB.prevPosition,
              nodeA.position - nodeB.position,
              nodeC.position - nodeB.position,
              nodeD.position - nodeB.position);

          if (!optT) {
            // CCD did not find intersection
            // TODO: Still should resolve static collisions?
            continue;
          }

          data.triCollisions.emplace_back(nodeA, nodeB, nodeC, nodeD);
        }
      }

      buckets.clear();

      if (nodeA.position.y < floorHeight) {
        data.staticCollisions.emplace_back(
            nodeA,
            glm::vec3(nodeA.position.x, floorHeight, nodeA.position.z));
      }
    }
  };

  std::vector<std::thread> threads;
  for (uint32_t threadId = 0; threadId < this->_options.threadCount;
       ++threadId) {
    threads.emplace_back(fnComputeCollisions, threadId);
  }

  for (std::thread& thread : threads) {
    thread.join();
  }

  this->_clearSpatialHashThread =
      std::thread([this]() { this->_spatialHashTris.clear(); });

  // Aggregate across per-thread results
  for (const ThreadData& data : this->_threadData) {
    for (const PointTriangleCollisionConstraint& collision :
         data.triCollisions) {
      // collision.setupGlobalForceVector(this->_forceVector);
      // collision.setupCollisionMatrix(this->_collisionMatrix);
      this->_triCollisions.push_back(collision);
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
  range.lengthX = static_cast<uint32_t>(adjustedDiameter.r);
  range.lengthY = static_cast<uint32_t>(adjustedDiameter.g);
  range.lengthZ = static_cast<uint32_t>(adjustedDiameter.b);

  return range;
}

SpatialHashGridCellRange Solver::TetCompRange::operator()(
    const Tetrahedron& tet,
    const SpatialHashGrid& grid) const {
  glm::vec3 x1 = nodes[tet.nodeIds[0]].position / grid.scale;
  glm::vec3 x2 = nodes[tet.nodeIds[1]].position / grid.scale;
  glm::vec3 x3 = nodes[tet.nodeIds[2]].position / grid.scale;
  glm::vec3 x4 = nodes[tet.nodeIds[3]].position / grid.scale;

  glm::vec3 min(std::numeric_limits<float>::max());
  glm::vec3 max(std::numeric_limits<float>::lowest());

  SpatialHashGridCellRange range{};
  range.minX = static_cast<int64_t>(
      glm::floor(glm::min(glm::min(glm::min(x1.x, x2.x), x3.x), x4.x)));
  range.minY = static_cast<int64_t>(
      glm::floor(glm::min(glm::min(glm::min(x1.y, x2.y), x3.y), x4.y)));
  range.minZ = static_cast<int64_t>(
      glm::floor(glm::min(glm::min(glm::min(x1.z, x2.z), x3.z), x4.z)));

  range.lengthX = glm::max(
      static_cast<uint32_t>(
          glm::max(glm::max(glm::max(x1.x, x2.x), x3.x), x4.x) - range.minX),
      1u);
  range.lengthY = glm::max(
      static_cast<uint32_t>(
          glm::max(glm::max(glm::max(x1.y, x2.y), x3.y), x4.y) - range.minY),
      1u);
  range.lengthZ = glm::max(
      static_cast<uint32_t>(
          glm::max(glm::max(glm::max(x1.z, x2.z), x3.z), x4.z) - range.minZ),
      1u);

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

  glm::vec3 x1 = nodes[triangle.nodeIds[0]].position / grid.scale;
  glm::vec3 x2 = nodes[triangle.nodeIds[1]].position / grid.scale;
  glm::vec3 x3 = nodes[triangle.nodeIds[2]].position / grid.scale;

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