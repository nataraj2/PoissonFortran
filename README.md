# Parallel Poisson solver using PETSc (C++ and Fortran)

<img src="Images/PetscSolution.png?raw=true&v=50" alt="your_alternative_text" width="50%" height="50%"><img src="Images/ExactSolution.png?raw=true&v=50" alt="your_alternative_text" width="50%" height="50%">
<img src="Images/Comparison.png?raw=true&v=50" alt="your_alternative_text" width="50%" height="50%">

This repository containes the code for solving a Poisson equation in parallel in two dimensions in PetSc. It solves the inhomogeneous 
Poisson equation $\nabla^2 u = f$ with Neumann boundary conditions on all boundaries.

## Installation and compilation

```1. git clone https://github.com/nataraj2/petsc.git```
```2. ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90```
3. Then set PETSC_DIR = <path>. The output of the above on the screen will tell the exact command. i
Just copy paste it. Add it to the `bashrc` file.
4. sh run_poisson2d_neumann.sh
5. python PoissonSolution.py - plots the contours of the exact and petsc solutions
