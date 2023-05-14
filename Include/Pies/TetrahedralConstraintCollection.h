#pragma once

#include "Node.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

namespace Pies {
struct alignas(16) TetrahedralConstraint {
  // TODO: Can we avoid making A specific to each tet's rest pose?
  Eigen::Matrix<float, 4, 4> AtB;
  glm::mat3 Qinv;

  int nodeIds[4];

  float w;
  float minStrain;
  float maxStrain;
  float minOmega;
  float maxOmega;

  TetrahedralConstraint(
      const Node& a,
      const Node& b,
      const Node& c,
      const Node& d,
      float w,
      float minStrain,
      float maxStrain,
      float minOmega,
      float maxOmega);
};

class TetrahedralConstraintCollection {
public:
  TetrahedralConstraintCollection() = default;
  TetrahedralConstraintCollection(
      std::vector<TetrahedralConstraint>&& tetConstraints,
      Eigen::SparseMatrix<float>& systemMatrix);
  ~TetrahedralConstraintCollection();

  TetrahedralConstraintCollection(TetrahedralConstraintCollection&& rhs);
  TetrahedralConstraintCollection&
  operator=(TetrahedralConstraintCollection&& rhs);
  TetrahedralConstraintCollection(const TetrahedralConstraintCollection& rhs) =
      delete;
  TetrahedralConstraintCollection&
  operator=(const TetrahedralConstraintCollection& rhs) = delete;

  void project(glm::vec3* devNodePositions);

  void setupGlobalForceVector(Eigen::MatrixXf& forceVector) const;

  size_t getTetCount() const { return this->_tets.size(); }

private:
  std::vector<TetrahedralConstraint> _tets;

  // SHARED CUDA RESOURCES:
  // Stacked buffer of tetrahedral constraint structs, set once at the beginning
  TetrahedralConstraint* _dev_tets = nullptr;

  // TODO: Remove
  // Stacked buffer of individual SVD results, each tet constraint will
  // correspond to 21 floats representing 3x3 matrix U (9 floats), vec3 S (3
  // floats), and 3x3 matrix V (9 floats). The matrices are in row-major order.
  // float* svdOut = nullp;

  // Stacked buffer of the result of the projection
  // This is a little obfuscated to understand from the code, but this is a
  // quantity needed for the serial global force vector setup in PD. Everything
  // up until this quantity can be done in parallel.
  Eigen::Matrix<float, 4, 3>* _dev_wAtBp = nullptr;
  std::vector<Eigen::Matrix<float, 4, 3>> _wAtBp;
};
} // namespace Pies