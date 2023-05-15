#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

namespace Pies {
struct Node;

struct alignas(16) TetrahedralConstraint {
  glm::mat3x4 At;
  glm::mat3 Xf_inv;

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

  // Stacked buffer of the result of the projection
  // This is a little obfuscated to understand from the code, but this is a
  // quantity needed for the serial global force vector setup in PD. Everything
  // up until this quantity can be done in parallel.
  // Note: These are 4(rows)x3(cols) matrices.
  glm::mat3x4* _dev_wAtBp = nullptr;
  std::vector<glm::mat3x4> _wAtBp;
};
} // namespace Pies