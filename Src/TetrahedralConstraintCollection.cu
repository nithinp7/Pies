
#include "Node.h"
#include "TetrahedralConstraintCollection.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// #include <svd3x3/svd3_cuda.h>
#include <SVD.h>

#define GLM_FORCE_CUDA
#include <glm/glm.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_transform.hpp>

__global__ void projectTets(
    Pies::TetrahedralConstraint* devTets,
    glm::vec3* devNodePositions,
    glm::mat3x4* wAtBp,
    int tetCount) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index >= tetCount) {
    return;
  }

  const Pies::TetrahedralConstraint& tet = devTets[index];

  const glm::vec3& x1 = devNodePositions[tet.nodeIds[0]];
  const glm::vec3& x2 = devNodePositions[tet.nodeIds[1]];
  const glm::vec3& x3 = devNodePositions[tet.nodeIds[2]];
  const glm::vec3& x4 = devNodePositions[tet.nodeIds[3]];

  glm::mat3 Xg(x2 - x1, x3 - x1, x4 - x1);

  // Deformation gradient
  glm::mat3 F = Xg * tet.Xf_inv;

  SVD::Mat3x3 F_(
      F[0][0],
      F[1][0],
      F[2][0],
      F[0][1],
      F[1][1],
      F[2][1],
      F[0][2],
      F[1][2],
      F[2][2]);
  SVD::SVDSet svdRes = SVD::svd(F_);

  glm::mat3 U(
      svdRes.U.m_00,
      svdRes.U.m_10,
      svdRes.U.m_20,
      svdRes.U.m_01,
      svdRes.U.m_11,
      svdRes.U.m_21,
      svdRes.U.m_02,
      svdRes.U.m_12,
      svdRes.U.m_22);
  glm::vec3 sigma(svdRes.S.m_00, svdRes.S.m_11, svdRes.V.m_22);
  glm::mat3 V(
      svdRes.V.m_00,
      svdRes.V.m_10,
      svdRes.V.m_20,
      svdRes.V.m_01,
      svdRes.V.m_11,
      svdRes.V.m_21,
      svdRes.V.m_02,
      svdRes.V.m_12,
      svdRes.V.m_22);

  // glm::mat3 U(0.0f);
  // glm::vec3 sigma(0.0f);
  // glm::mat3 V(0.0f);

  // svd(
  //     // Remember glm is stored column major
  //     F[0][0],
  //     F[1][0],
  //     F[2][0],
  //     F[0][1],
  //     F[1][1],
  //     F[2][1],
  //     F[0][2],
  //     F[1][2],
  //     F[2][2],
  //     U[0][0],
  //     U[1][0],
  //     U[2][0],
  //     U[0][1],
  //     U[1][1],
  //     U[2][1],
  //     U[0][2],
  //     U[1][2],
  //     U[2][2],
  //     sigma[0],
  //     sigma[1],
  //     sigma[2],
  //     V[0][0],
  //     V[1][0],
  //     V[2][0],
  //     V[0][1],
  //     V[1][1],
  //     V[2][1],
  //     V[0][2],
  //     V[1][2],
  //     V[2][2]);

  if (glm::determinant(F) < 0.0f) {
    sigma[2] *= -1.0f;
  }

  sigma = glm::clamp(sigma, tet.minStrain, tet.maxStrain);

  // const uint32_t COMP_D_ITERS = 10;
  // glm::vec3 D(0.0f);
  // for (uint32_t i = 0; i < COMP_D_ITERS; ++i) {
  //   glm::vec3 sigmaPlusD = sigma + D;
  //   float product = sigmaPlusD.x * sigmaPlusD.y * sigmaPlusD.z;
  //   float omega = glm::clamp(product, tet.minOmega, tet.maxOmega);
  //   float C = product - omega;
  //   glm::vec3 gradC(
  //       sigmaPlusD.y * sigmaPlusD.z,
  //       sigmaPlusD.x * sigmaPlusD.z,
  //       sigmaPlusD.x * sigmaPlusD.y);
  //   D = (glm::dot(gradC, D) - C) * gradC / glm::dot(gradC, gradC);
  // }

  // sigma += D;

  // The "fixed" deformation gradient
  glm::mat3
      sigma_(sigma.x, 0.0f, 0.0f, 0.0f, sigma.y, 0.0f, 0.0f, 0.0f, sigma.z);
  glm::mat3 T = U * sigma_ * glm::transpose(V);

  // Note that Tt = Bp in our construction

  wAtBp[index] = tet.w * tet.At * glm::transpose(T);
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
  glm::mat4x3 M(0.0f);
  M[0][0] = -1.0f;
  M[0][1] = -1.0f;
  M[0][2] = -1.0f;

  M[1][0] = 1.0f;
  M[2][1] = 1.0f;
  M[3][2] = 1.0f;

  // Converts barycentric coords to world differential cords
  glm::mat3 Xf(
      b.position - a.position,
      c.position - a.position,
      d.position - a.position);
  this->Xf_inv = glm::inverse(Xf);

  this->At = glm::transpose(glm::transpose(this->Xf_inv) * M);
}

TetrahedralConstraintCollection::TetrahedralConstraintCollection(
    std::vector<TetrahedralConstraint>&& tetConstraints,
    Eigen::SparseMatrix<float>& systemMatrix)
    : _tets(std::move(tetConstraints)) {

  this->_wAtBp.resize(this->_tets.size());

  for (const TetrahedralConstraint& tet : this->_tets) {
    // B is identity here
    glm::mat4 AtA = tet.At * glm::transpose(tet.At);
    for (int i = 0; i < 4; ++i) {
      int nodeId_i = tet.nodeIds[i];
      for (int j = 0; j < 4; ++j) {
        int nodeId_j = tet.nodeIds[j];
        systemMatrix.coeffRef(nodeId_i, nodeId_j) += tet.w * AtA[j][i];
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
  cudaMalloc(&this->_dev_wAtBp, sizeof(glm::mat3x4) * this->_tets.size());

  cudaDeviceSynchronize();
}

TetrahedralConstraintCollection::TetrahedralConstraintCollection(
    TetrahedralConstraintCollection&& rhs)
    : _tets(std::move(rhs._tets)), _wAtBp(std::move(rhs._wAtBp)) {
  if (this->_dev_tets) {
    cudaFree(this->_dev_tets);
  }

  if (this->_dev_wAtBp) {
    cudaFree(this->_dev_wAtBp);
  }

  cudaDeviceSynchronize();

  this->_dev_tets = rhs._dev_tets;
  this->_dev_wAtBp = rhs._dev_wAtBp;

  rhs._dev_tets = nullptr;
  rhs._dev_wAtBp = nullptr;
}

TetrahedralConstraintCollection& TetrahedralConstraintCollection::operator=(
    TetrahedralConstraintCollection&& rhs) {
  if (this->_dev_tets) {
    cudaFree(this->_dev_tets);
  }

  if (this->_dev_wAtBp) {
    cudaFree(this->_dev_wAtBp);
  }

  cudaDeviceSynchronize();

  this->_tets = std::move(rhs._tets);
  this->_dev_tets = rhs._dev_tets;
  this->_dev_wAtBp = rhs._dev_wAtBp;
  this->_wAtBp = std::move(rhs._wAtBp);

  rhs._dev_tets = nullptr;
  rhs._dev_wAtBp = nullptr;

  return *this;
}

TetrahedralConstraintCollection::~TetrahedralConstraintCollection() {
  if (this->_dev_tets) {
    cudaFree(this->_dev_tets);
  }

  if (this->_dev_wAtBp) {
    cudaFree(this->_dev_wAtBp);
  }

  cudaDeviceSynchronize();
}

void TetrahedralConstraintCollection::project(glm::vec3* devNodePositions) {
  int tetCount = static_cast<int>(this->_tets.size());

  if (tetCount == 0) {
    return;
  }

  int threadCount = 256;
  int pblks = int((tetCount + threadCount - 1) / threadCount);
  projectTets<<<pblks, threadCount>>>(
      this->_dev_tets,
      devNodePositions,
      this->_dev_wAtBp,
      tetCount);

  cudaMemcpy(
      this->_wAtBp.data(),
      this->_dev_wAtBp,
      sizeof(glm::mat3x4) * this->_tets.size(),
      cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
}

void TetrahedralConstraintCollection::setupGlobalForceVector(
    Eigen::MatrixXf& forceVector) const {
  for (size_t tetId = 0; tetId < this->_tets.size(); ++tetId) {
    const TetrahedralConstraint& tet = this->_tets[tetId];
    const glm::mat3x4& wAtBp = this->_wAtBp[tetId];
    for (uint32_t i = 0; i < 4; ++i) {
      uint32_t nodeId_i = tet.nodeIds[i];
      forceVector.coeffRef(nodeId_i, 0) += wAtBp[0][i];
      forceVector.coeffRef(nodeId_i, 1) += wAtBp[1][i];
      forceVector.coeffRef(nodeId_i, 2) += wAtBp[2][i];
    }
  }
}
} // namespace Pies