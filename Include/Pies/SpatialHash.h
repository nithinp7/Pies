#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <math.h>
#include <parallel_hashmap/phmap.h>

#include <cstdint>
#include <functional>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

namespace Pies {
struct SpatialHashGridCellId {
  int64_t x = 0;
  int64_t y = 0;
  int64_t z = 0;

  bool operator==(const SpatialHashGridCellId& rhs) const {
    return this->x == rhs.x && this->y == rhs.y && this->z == rhs.z;
  }
};
} // namespace Pies

namespace std {
template <> struct hash<Pies::SpatialHashGridCellId> {
  size_t operator()(const Pies::SpatialHashGridCellId& id) const noexcept {
    int64_t hash = (id.x * 92837111) ^ (id.y * 689287499) ^ (id.z * 283923481);
    return static_cast<size_t>(std::abs(hash));
  }
};
} // namespace std

namespace Pies {

template <typename TValue> struct SpatialHashGridCellBucket {
  std::vector<TValue*> values;
};

struct SpatialHashGridCellRange {
  int64_t minX;
  int64_t minY;
  int64_t minZ;

  uint32_t lengthX;
  uint32_t lengthY;
  uint32_t lengthZ;
};

/**
 * @brief Defined to be a uniform grid such that, in local space, the cell
 * lengths are 1 and the cells are aligned to the XYZ axes. Transforms between
 * grid space and world space are provided.
 */
struct SpatialHashGrid {
  glm::mat4 gridToWorld;
  glm::mat4 worldToGrid;
};

template <typename TValue, typename TCompRange> class SpatialHash {
public:
  SpatialHash(const glm::mat4& transform = glm::mat4(1.0f))
      : _grid({transform, glm::affineInverse(transform)}) {}

  const glm::mat4& localToWorld() const { return this->_gridTransform; }

  const glm::mat4& worldToLocal() const { return this->_invGridTransform; }

  void parallelBulkInsert(
      std::vector<TValue>& values,
      const TCompRange& compRangeFn) {
    // NOTE: Based on example given in:
    // https://greg7mdp.github.io/parallel-hashmap/
    constexpr size_t numThreads = 8; // has to be a power of two
    std::unique_ptr<std::thread> threads[numThreads];
    auto threadFn = [&grid = this->_grid,
                     &hash = this->_hashMap,
                     &values,
                     &compRangeFn,
                     numThreads](size_t thread_idx) {
      size_t modulo =
          hash.subcnt() / numThreads; // subcnt() returns the number of submaps

      for (size_t i = 0; i < values.size(); ++i) {
        TValue& value = values[i];
        // const GridCellRange& range = gridCellRanges[i];
        SpatialHashGridCellRange range = compRangeFn(value, grid);

        for (uint32_t dx = 0; dx < range.lengthX; ++dx) {
          for (uint32_t dy = 0; dy < range.lengthY; ++dy) {
            for (uint32_t dz = 0; dz < range.lengthZ; ++dz) {
              // Compute grid cell id
              SpatialHashGridCellId id{
                  range.minX + dx,
                  range.minY + dy,
                  range.minZ + dz};
              // Compute its hash
              size_t hashVal = hash.hash(id);
              // Compute the submap index for this hash
              size_t idx = hash.subidx(hashVal);
              if (idx / modulo == thread_idx) {
                // If the submap is suitable for this thread, add the value to
                // the grid cell bucket.
                auto it = hash.find(id, hashVal);
                if (it != hash.end()) {
                  it->second.values.push_back(&value);
                } else {
                  hash.emplace_with_hash(
                      hashVal,
                      std::make_pair(
                          id,
                          SpatialHashGridCellBucket<TValue>{{&value}}));
                }
              }
            }
          }
        }
      }
    };

    for (size_t i = 0; i < numThreads; ++i) {
      threads[i].reset(new std::thread(threadFn, i));
    }

    // wait for the threads to finish their work and exit
    for (size_t i = 0; i < numThreads; ++i) {
      threads[i]->join();
    }
  }

  void clear() { this->_hashMap.clear(); }

private:
  SpatialHashGrid _grid;
  phmap::parallel_flat_hash_map<
      SpatialHashGridCellId,
      SpatialHashGridCellBucket<TValue>>
      _hashMap;
};
} // namespace Pies