#include "Solver.h"

#include <glm/gtc/matrix_transform.hpp>

#include <cstdint>
#include <thread>

namespace Pies {

Solver::Solver(const SolverOptions& options)
    : _options(options),
      _spatialHashNodes(glm::vec3(0.0f), options.gridSpacing),
      _spatialHashTets(glm::vec3(0.0f), options.gridSpacing) {
  this->_threadData.resize(this->_options.threadCount);
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
      node.position += h * node.velocity + h2 * node.force * node.invMass;

      glm::vec3 Msn_h2 = node.position / node.invMass / h2;
      this->_Msn_h2.coeffRef(i, 0) = Msn_h2.x;
      this->_Msn_h2.coeffRef(i, 1) = Msn_h2.y;
      this->_Msn_h2.coeffRef(i, 2) = Msn_h2.z;
    }

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

      for (ShapeMatchingConstraint& constraint : this->_shapeConstraints) {
        constraint.projectToAuxiliaryVariable(this->_nodes);
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

      for (ShapeMatchingConstraint& constraint : this->_shapeConstraints) {
        constraint.setupGlobalForceVector(this->_forceVector);
      }

      // this->_computeCollisions();
      this->_parallelComputeCollisions();

      this->_stiffnessAndCollisionMatrix =
          this->_stiffnessMatrix + this->_collisionMatrix;
      // TODO: Just for testing, don't solve collision this way.
      this->_pLltDecomp =
          std::make_unique<Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>>(
              this->_stiffnessAndCollisionMatrix);

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

      // float p; //???
      // omega = 4.0 / (4.0 - p * p * omega);
    }

    // Update node velocities
    for (uint32_t i = 0; i < nodeCount; ++i) {
      Node& node = this->_nodes[i];
      node.velocity = (1.0f - this->_options.damping) *
                      (node.position - node.prevPosition) / h;

      node.prevPosition = node.position;
      this->_vertices[i].position = node.position;
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

        glm::vec3 diff = pOtherNode->position - node.position;
        float dist = glm::length(diff);
        float dispLength = node.radius + pOtherNode->radius - dist;

        if (dispLength <= 0.0f) {
          continue;
        }

        glm::vec3 disp;
        if (dist > 0.00001f) {
          disp = dispLength * diff / dist;
        } else {
          disp = glm::vec3(dispLength, 0.0f, 0.0f);
        }

        CollisionConstraint& collision =
            this->_collisions.emplace_back(node, *pOtherNode, disp);
        // collision.setupCollisionMatrix(this->_collisionMatrix);
        // collision.setupGlobalForceVector(this->_forceVector);
      }
    }

    if (node.position.y - node.radius < this->_options.floorHeight) {
      StaticCollisionConstraint& collision =
          this->_staticCollisions.emplace_back(
              node,
              glm::vec3(
                  node.position.x,
                  this->_options.floorHeight + node.radius,
                  node.position.z));
      // collision.setupCollisionMatrix(this->_collisionMatrix);
      // collision.setupGlobalForceVector(this->_forceVector);
    }
  }

  for (CollisionConstraint& collision : this->_collisions) {
    collision.setupCollisionMatrix(this->_collisionMatrix);
    collision.setupGlobalForceVector(this->_forceVector);
  }

  for (StaticCollisionConstraint& collision : this->_staticCollisions) {
    collision.setupCollisionMatrix(this->_collisionMatrix);
    collision.setupGlobalForceVector(this->_forceVector);
  }
}

void Solver::_parallelComputeCollisions() {
  // Update collisions
  this->_collisions.clear();
  this->_staticCollisions.clear();
  this->_collisionMatrix.setZero();
  this->_stiffnessAndCollisionMatrix.setZero();

  this->_spatialHashNodes.clear();
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

          glm::vec3 diff = pOtherNode->position - node.position;
          float dist = glm::length(diff);
          float dispLength = node.radius + pOtherNode->radius - dist;

          if (dispLength <= 0.0f) {
            continue;
          }

          glm::vec3 disp;
          if (dist > 0.00001f) {
            disp = dispLength * diff / dist;
          } else {
            disp = glm::vec3(dispLength, 0.0f, 0.0f);
          }

          data.collisions.emplace_back(node, *pOtherNode, disp);
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

  // Aggregate across per-thread results
  for (const ThreadData& data : this->_threadData) {
    for (const CollisionConstraint& collision : data.collisions) {
      collision.setupCollisionMatrix(this->_collisionMatrix);
      collision.setupGlobalForceVector(this->_forceVector);
      this->_collisions.push_back(collision);
    }

    for (const StaticCollisionConstraint& collision : data.staticCollisions) {
      collision.setupCollisionMatrix(this->_collisionMatrix);
      collision.setupGlobalForceVector(this->_forceVector);
      this->_staticCollisions.push_back(collision);
    }
  }
}

SpatialHashGridCellRange Solver::NodeCompRange::operator()(
    const Node& node,
    const SpatialHashGrid& grid) const {
  float radiusPadding = 0.0f;
  float gridLocalRadius = (node.radius + radiusPadding) / grid.scale;
  glm::vec3 gridLocalPos = (node.position - grid.translation) / grid.scale;
  glm::vec3 gridLocalMin = gridLocalPos - glm::vec3(gridLocalRadius);

  SpatialHashGridCellRange range{};
  range.minX = static_cast<int64_t>(glm::floor(gridLocalMin.x));
  range.minY = static_cast<int64_t>(glm::floor(gridLocalMin.y));
  range.minZ = static_cast<int64_t>(glm::floor(gridLocalMin.z));

  glm::vec3 adjustedDiameter =
      glm::round(glm::fract(gridLocalMin) + glm::vec3(2 * gridLocalRadius));
  range.lengthX = static_cast<uint32_t>(adjustedDiameter.r);
  range.lengthY = static_cast<uint32_t>(adjustedDiameter.g);
  range.lengthZ = static_cast<uint32_t>(adjustedDiameter.b);

  return range;
}

SpatialHashGridCellRange Solver::TetCompRange::operator()(
    const Tetrahedron& tet,
    const SpatialHashGrid& grid) const {
  glm::vec3 x1 =
      (nodes[tet.nodeIds[0]].position - grid.translation) / grid.scale;
  glm::vec3 x2 =
      (nodes[tet.nodeIds[1]].position - grid.translation) / grid.scale;
  glm::vec3 x3 =
      (nodes[tet.nodeIds[2]].position - grid.translation) / grid.scale;
  glm::vec3 x4 =
      (nodes[tet.nodeIds[3]].position - grid.translation) / grid.scale;

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
} // namespace Pies