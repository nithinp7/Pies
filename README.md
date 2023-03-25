# Pies

## Current Features
- PBD Solver
- Constraints
  - Fixed-Position
  - Distance
  - Tetrahedral Volume Preservation
  - Tetrahedral Strain Limiting
- Dynamic Collisions
  - Per-iteration collision detection using spatial hash (implemented on top of the parallel hashmap library)
  - Node-Node collision detection / resolution

<img src="https://github.com/nithinp7/Pies/blob/main/Media/PBDCollisions.gif" width=125%>
