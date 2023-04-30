# Pies

Pies is a constraint- and particle-based, soft-body physics engine based on the paper _"Projective Dynamics: Fusing Constraint Projections for Fast Simulation", Bouaziz et. al 2014_. 

<img src="https://github.com/nithinp7/Pies/blob/main/Media/PiesMaya4.gif">

<img src="https://github.com/nithinp7/Pies/blob/main/Media/PiesMaya1.gif">

## Current Features
- Projective dynamics solver
- Position-based dynamics solver
- Constraints
  - Fixed-Position
  - Distance
  - Tetrahedral Volume Preservation
  - Tetrahedral Strain Limiting
  - Shape matching
  - Goal matching
  - Collision constraints
  - Friction
- Dynamic Collisions
  - Per-iteration collision detection using spatial hash (implemented on top of the parallel hashmap library)
    - Parallelized spatial hash construction
    - Parallelized collision detection
  - Triangle-Triangle continuous collision detection.
  - Node-Node collision detection / resolution

## Integrations
Pies has been integrated successfully in two places so far, check out the linked repositories for more information:
- Pies was originally prototyped and tested with the [Pies for Althea](https://github.com/nithinp7/PiesForAlthea) test application built on top of the [Althea](https://github.com/nithinp7/Althea) rendering engine.
- Pies has since been integrated into Autodesk Maya with the [Pies for Maya](https://github.com/nithinp7/PiesForMaya) plugin.

## Pies for Maya

- [Pies for Maya Repository](https://github.com/nithinp7/PiesForMaya)

<img src="https://github.com/nithinp7/Pies/blob/main/Media/PiesMaya2.gif" width=800>

<img src="https://github.com/nithinp7/Pies/blob/main/Media/PiesMaya3.gif" width=800>

## Pies for Althea

- [Pies for Althea Repository](https://github.com/nithinp7/PiesForAlthea)

<img src="https://github.com/nithinp7/Pies/blob/main/Media/PiesAlthea1.gif" width=800>

<img src="https://github.com/nithinp7/Pies/blob/main/Media/PiesAlthea2.png" width=800>

## Original Prototype (Old)
### Projective Dynamics
<img src="https://github.com/nithinp7/Pies/blob/main/Media/PD.gif" width=800>

### Position-Based Dynamics
<img src="https://github.com/nithinp7/Pies/blob/main/Media/PBDCollisions.gif" width=800>
