
#include "TetrahedralConstraintCollection.h"

#include <svd3x3/svd3_cuda.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define GLM_FORCE_CUDA
#include <glm/glm.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_transform.hpp>

__global__ void projectTets(
    Pies::TetrahedralConstraint* devTets,
    glm::vec3* devNodePositions,
    Eigen::Matrix<float, 4, 3>* wAtBp,
    int tetCount) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= tetCount) {
    return;
  }

  // TODO: Avoid the awkward switching back and forth between
  // Eigen and glm here

  const Pies::TetrahedralConstraint& tet = devTets[i];

  const glm::vec3& x1 = devNodePositions[tet.nodeIds[0]];
  const glm::vec3& x2 = devNodePositions[tet.nodeIds[1]];
  const glm::vec3& x3 = devNodePositions[tet.nodeIds[2]];
  const glm::vec3& x4 = devNodePositions[tet.nodeIds[3]];

  glm::mat3 P(x2 - x1, x3 - x1, x4 - x1);

  // Deformation gradient
  glm::mat3 F = P * tet.Qinv;

  Eigen::Matrix3f F_;
  F_ << F[0][0], F[0][1], F[0][2], F[1][0], F[1][1], F[1][2], F[2][0], F[2][1],
      F[2][2];

  Eigen::Matrix3f U;
  Eigen::Vector3f singularValues;
  Eigen::Matrix3f V;

  svd(
      // Remember glm is stored column major
      F[0][0],
      F[1][0],
      F[2][0],
      F[0][1],
      F[1][1],
      F[2][1],
      F[0][2],
      F[1][2],
      F[2][2],
      U.coeffRef(0, 0),
      U.coeffRef(0, 1),
      U.coeffRef(0, 2),
      U.coeffRef(1, 0),
      U.coeffRef(1, 1),
      U.coeffRef(1, 2),
      U.coeffRef(2, 0),
      U.coeffRef(2, 1),
      U.coeffRef(2, 2),
      singularValues[0],
      singularValues[1],
      singularValues[2],
      V.coeffRef(0, 0),
      V.coeffRef(0, 1),
      V.coeffRef(0, 2),
      V.coeffRef(1, 0),
      V.coeffRef(1, 1),
      V.coeffRef(1, 2),
      V.coeffRef(2, 0),
      V.coeffRef(2, 1),
      V.coeffRef(2, 2));

  glm::vec3 sigma(singularValues[0], singularValues[1], singularValues[2]);

  const uint32_t COMP_D_ITERS = 10;
  glm::vec3 D(0.0f);
  for (uint32_t i = 0; i < COMP_D_ITERS; ++i) {
    glm::vec3 sigmaPlusD = sigma + D;
    float product = sigmaPlusD.x * sigmaPlusD.y * sigmaPlusD.z;
    float omega = glm::clamp(product, tet.minOmega, tet.maxOmega);
    float C = product - omega;
    glm::vec3 gradC(
        sigmaPlusD.y * sigmaPlusD.z,
        sigmaPlusD.x * sigmaPlusD.z,
        sigmaPlusD.x * sigmaPlusD.y);
    D = (glm::dot(gradC, D) - C) * gradC / glm::dot(gradC, gradC);
  }

  singularValues[0] += D.x;
  singularValues[1] += D.y;
  singularValues[2] += D.z;

  for (uint32_t i = 0; i < 3; ++i) {
    singularValues[i] =
        glm::clamp(singularValues[i], tet.minStrain, tet.maxStrain);
  }

  if (glm::determinant(F) < 0.0f) {
    singularValues[2] *= -1.0f;
  }

  // The "fixed" deformation gradient
  Eigen::Matrix3f Fhat =
      U * singularValues.asDiagonal() * V.transpose();

  Eigen::Matrix<float, 4, 3> P1 = Eigen::Matrix<float, 4, 3>::Zero();
  // TODO: Fhat.transpose()???
  P1.block(1, 0, 3, 3) = Fhat;

  wAtBp[i] = tet.w * tet.AtB * P1;
}

namespace Pies {
TetrahedralConstraint::TetrahedralConstraint(
    const Node& a,
    const Node& b,
    const Node& c,
    const Node& d,
    float w_,
    float minStrain_,
    float maxStrain_,
    float minOmega_,
    float maxOmega_)
    : nodeIds{static_cast<int>(a.id), static_cast<int>(b.id), static_cast<int>(c.id), static_cast<int>(d.id)},
      w(w_),
      minStrain(minStrain_),
      maxStrain(maxStrain_),
      minOmega(minOmega_),
      maxOmega(maxOmega_) {

  // Converts world positions to differential coords
  Eigen::Matrix<float, 3, 4> worldToDiff = Eigen::Matrix<float, 3, 4>::Zero();
  worldToDiff.coeffRef(0, 0) = -1.0f;
  worldToDiff.coeffRef(1, 0) = -1.0f;
  worldToDiff.coeffRef(2, 0) = -1.0f;

  worldToDiff.coeffRef(0, 1) = 1.0f;
  worldToDiff.coeffRef(1, 2) = 1.0f;
  worldToDiff.coeffRef(2, 3) = 1.0f;

  // Converts barycentric coords to world differential cords
  glm::mat3 baryToDiff(
      b.position - a.position,
      c.position - a.position,
      d.position - a.position);
  glm::mat3 diffToBary = glm::inverse(baryToDiff);

  Eigen::Matrix3f diffToBary_;
  diffToBary_ << diffToBary[0][0], diffToBary[0][1], diffToBary[0][2],
      diffToBary[1][0], diffToBary[1][1], diffToBary[1][2], diffToBary[2][0],
      diffToBary[2][1], diffToBary[2][2];

  Eigen::Matrix<float, 3, 4> A_ = diffToBary_ * worldToDiff;
  Eigen::Matrix4f A = Eigen::Matrix4f::Zero();
  // A.coeffRef(1, 0) = -1.0f;
  // A.coeffRef(2, 0) = -1.0f;
  // A.coeffRef(3, 0) = -1.0f;

  // A.coeffRef(1, 1) = 1.0f;
  // A.coeffRef(2, 2) = 1.0f;
  // A.coeffRef(3, 3) = 1.0f;

  A.row(0) << 0.0f, 0.0f, 0.0f, 0.0f;
  A.row(1) = A_.row(0);
  A.row(2) = A_.row(1);
  A.row(3) = A_.row(2);

  // B is identity

  this->Qinv = diffToBary;
  this->AtB = A;
}

TetrahedralConstraintCollection::TetrahedralConstraintCollection(
    std::vector<TetrahedralConstraint>&& tetConstraints,
    Eigen::SparseMatrix<float>& systemMatrix)
    : _tets(std::move(tetConstraints)) {

  this->_wAtBp.resize(this->_tets.size());

  for (const TetrahedralConstraint& tet : this->_tets) {
    // B is identity here
    Eigen::Matrix<float, 4, 4> AtA = tet.AtB.transpose() * tet.AtB;
    for (int i = 0; i < 4; ++i) {
      int nodeId_i = tet.nodeIds[i];
      for (int j = 0; j < 4; ++j) {
        int nodeId_j = tet.nodeIds[j];
        systemMatrix.coeffRef(nodeId_i, nodeId_j) += tet.w * AtA.coeff(i, j);
      }
    }
  }

  // Create device memory for the tet constraints
  cudaMalloc(
      &this->_dev_tets,
      sizeof(TetrahedralConstraint) * this->_tets.size());
  cudaMemcpy(
      this->_dev_tets,
      this->_tets.data(),
      sizeof(TetrahedralConstraint) * this->_tets.size(),
      cudaMemcpyHostToDevice);

  // Create device memory to hold the projection output.
  cudaMalloc(
      &this->_dev_wAtBp,
      sizeof(Eigen::Matrix<float, 4, 3>) * this->_tets.size());
}

TetrahedralConstraintCollection::~TetrahedralConstraintCollection() {
  cudaFree(this->_dev_tets);
  cudaFree(this->_dev_wAtBp);
}

void TetrahedralConstraintCollection::project(glm::vec3* devNodePositions) {
  int tetCount = static_cast<int>(this->_tets.size());

  int threadCount = 504;
  int pblks = int((tetCount + threadCount - 1) / threadCount);
  projectTets<<<pblks, threadCount>>>(
      this->_dev_tets,
      devNodePositions,
      this->_dev_wAtBp,
      tetCount);

  cudaMemcpy(
      this->_wAtBp.data(),
      this->_dev_wAtBp,
      sizeof(Eigen::Matrix<float, 4, 3>) * this->_tets.size(),
      cudaMemcpyDeviceToHost);
}

void TetrahedralConstraintCollection::setupGlobalForceVector(
    Eigen::MatrixXf& forceVector) const {
  for (size_t tetId = 0; tetId < this->_tets.size(); ++tetId) {
    const TetrahedralConstraint& tet = this->_tets[tetId];
    const Eigen::Matrix<float, 4, 3>& wAtBp = this->_wAtBp[tetId];
    for (uint32_t i = 0; i < 4; ++i) {
      uint32_t nodeId_i = tet.nodeIds[i];
      forceVector.coeffRef(nodeId_i, 0) += wAtBp.coeff(i, 0);
      forceVector.coeffRef(nodeId_i, 1) += wAtBp.coeff(i, 1);
      forceVector.coeffRef(nodeId_i, 2) += wAtBp.coeff(i, 2);
    }
  }
}
} // namespace Pies