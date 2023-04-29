# Pies

Pies is a constraint- and particle-based, soft-body physics engine based on the paper _"Projective Dynamics: Fusing Constraint Projections for Fast Simulation", Bouaziz et. al 2014_. 

<img src="https://github.com/nithinp7/Pies/blob/main/Media/PiesMaya1.gif">

## Current Features
- Projective dynamics solver [Bouaziz et. al 2014]
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

<img src="https://github.com/nithinp7/Pies/blob/main/Media/PiesMaya2.gif">

## Pies for Althea


## Original Prototype
### Projective Dynamics
YouTube links:
- [Projective Dynamics video 1](https://youtu.be/3i7e_A1Btxo)
- [Projective Dynamics video 2](https://youtu.be/SwnAaLc8OYM)
<img src="https://github.com/nithinp7/Pies/blob/main/Media/PD.gif" width=800>

### Position-Based Dynamics
YouTube link:
- [Position Based Dynamics video](https://youtu.be/IaQ6OIHqq7Y)
<img src="https://github.com/nithinp7/Pies/blob/main/Media/PBDCollisions.gif" width=800>

## Screenshots
<img src="https://github.com/nithinp7/Pies/blob/main/Media/PBDScene.png" width=600>
<img src="https://github.com/nithinp7/Pies/blob/main/Media/Reflections.png" width=600>
