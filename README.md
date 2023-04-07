# Pies

## Current Features
- Projective dynamics solver [Bouaziz et. al 2014]
- Position-based dynamics solver
- Constraints
  - Fixed-Position
  - Distance
  - Tetrahedral Volume Preservation
  - Tetrahedral Strain Limiting
- Dynamic Collisions
  - Per-iteration collision detection using spatial hash (implemented on top of the parallel hashmap library)
    - Parallelized spatial hash construction
    - Parallelized collision detection
  - Node-Node collision detection / resolution

## Plans
- Currently Pies is being tested with the [Pies for Althea](https://github.com/nithinp7/PiesForAlthea) test application built on top of the [Althea](https://github.com/nithinp7/Althea) engine.
- We plan to integrate Pies into Maya soon as well.
- More details will be added once core features are complete.

## Projective Dynamics
YouTube links:
- [Projective Dynamics video 1](https://youtu.be/3i7e_A1Btxo)
- [Projective Dynamics video 2](https://youtu.be/SwnAaLc8OYM)
<img src="https://github.com/nithinp7/Pies/blob/main/Media/PD.gif" width=800>

## Position-Based Dynamics
YouTube link:
- [Position Based Dynamics video](https://youtu.be/IaQ6OIHqq7Y)
<img src="https://github.com/nithinp7/Pies/blob/main/Media/PBDCollisions.gif" width=800>

## Screenshots
<img src="https://github.com/nithinp7/Pies/blob/main/Media/PBDScene.png" width=600>
<img src="https://github.com/nithinp7/Pies/blob/main/Media/Reflections.png" width=600>
