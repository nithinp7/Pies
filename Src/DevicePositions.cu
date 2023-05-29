#include "DevicePositions.h"
#include "Node.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cassert>

namespace Pies {
DevicePositions::DevicePositions(size_t count) 
  : _count(count) {
  this->_scratch.resize(count);
  cudaError_t err = cudaMalloc(&this->_devPositions, sizeof(glm::vec3) * count);
  
  cudaDeviceSynchronize();
}

DevicePositions::~DevicePositions() {
  if (this->_devPositions) {
    cudaFree(this->_devPositions);
  }

  cudaDeviceSynchronize();
}

DevicePositions::DevicePositions(DevicePositions&& rhs) 
  : _count(rhs._count),
    _devPositions(rhs._devPositions),
    _scratch(std::move(rhs._scratch)) {
  rhs._devPositions = nullptr;
  rhs._count = 0;
}

DevicePositions& DevicePositions::operator=(DevicePositions&& rhs) {
  this->_count = rhs._count;
  this->_devPositions = rhs._devPositions;
  this->_scratch = std::move(rhs._scratch);

  rhs._devPositions = nullptr;
  rhs._count = 0;

  return *this;
}

void DevicePositions::upload(const std::vector<Node>& nodes) {
  assert(nodes.size() == this->_count);

  if (nodes.size() == 0) {
    return;
  }
  
  // Scalarize / parallelize this??
  for (size_t i = 0; i < this->_count; ++i) {
    this->_scratch[i] = nodes[i].position;
  }

  cudaMemcpy(this->_devPositions, this->_scratch.data(), sizeof(glm::vec3) * this->_count, cudaMemcpyHostToDevice);
  
  cudaDeviceSynchronize();
}
} // namespace Pies