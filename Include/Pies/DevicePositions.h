#pragma once

#include <glm/glm.hpp>

#include <vector>
#include <cstdint>

namespace Pies {
struct Node;

class DevicePositions {
public:
  DevicePositions() = default;
  DevicePositions(size_t count);
  ~DevicePositions();

  DevicePositions(DevicePositions&& rhs);
  DevicePositions& operator=(DevicePositions&& rhs);
  
  DevicePositions(const DevicePositions& rhs) = delete;
  DevicePositions& operator=(const DevicePositions& rhs) = delete;

  void upload(const std::vector<Node>& nodes);

  glm::vec3* getDeviceBuffer() const {
    return this->_devPositions;
  }

private:
  size_t _count = 0;
  glm::vec3* _devPositions = nullptr;

  std::vector<glm::vec3> _scratch;
};
} // namespace Pies