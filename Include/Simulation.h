#pragma once

#include "Node.h"
#include "Constraints.h"

#include <Althea/Application.h>
#include <Althea/DrawContext.h>
#include <Althea/DynamicVertexBuffer.h>
#include <Althea/GraphicsPipeline.h>
#include <Althea/IndexBuffer.h>
#include <Althea/SingleTimeCommandBuffer.h>
#include <Althea/InputManager.h>

#include <vector>

using namespace AltheaEngine;

namespace PiesForAlthea {
class Simulation {
public:
  static void initInputBindings(InputManager& inputManager);
  static void buildPipeline(GraphicsPipelineBuilder& builder);

  Simulation(
      const Application& app,
      SingleTimeCommandBuffer& commandBuffer);

  void tick(
      const Application& app,
      float deltaTime);
  void draw(const DrawContext& context) const;

private:
  std::vector<Node> _nodes;
  std::vector<PositionConstraint> _positionConstraints;
  std::vector<DistanceConstraint> _distanceConstraints;
  DynamicVertexBuffer<glm::vec3> _vertexBuffer;
  IndexBuffer _indexBuffer;
  
  // Scratch vertices vector to avoid reallocation
  std::vector<glm::vec3> _vertices;
};
} // namespace PiesForAlthea