# Parallel Poisson solver using PETSc (C++ and Fortran)

# 
<img src="Images/PetscSolution.png?raw=true&v=50" alt="your_alternative_text" width="50%" height="50%"><img src="Images/ExactSolution.png?raw=true&v=50" alt="your_alternative_text" width="50%" height="50%">
<img src="Images/Comparison.png?raw=true&v=50" alt="your_alternative_text" width="50%" height="50%">

This repository containes the code for solving a Poisson equation in parallel in two dimensions in PetSc. It solves the inhomogeneous 
Poisson equation with Neumann boundary conditions on all boundaries.

## Governing equation, domain and boundary conditions
The Poisson equation in two-dimensions is given by  
$\cfrac{\partial^2 u}{\partial x^2} + \cfrac{\partial^2 u}{\partial y^2} = f$.  
In this example $f = 8\pi^2\cos(2\pi x)cos(2\pi y)$ and the domain is given by 
$(0,1)\times(0,1)$. Neumann boundary conditions are imposed on all boundaries 
$\cfrac{\partial u}{\partial n} = 0$. The exact solution is given by 
$u(x,y) = \cos(2\pi x)cos(2\pi y)$. 


## Installation and compilation

```1. git clone https://github.com/nataraj2/petsc.git```
```2. ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90```
3. Then set PETSC_DIR = <path>. The output of the above on the screen will tell the exact command. i
Just copy paste it. Add it to the `bashrc` file.
4. sh run_poisson2d_neumann.sh
5. python PoissonSolution.py - plots the contours of the exact and petsc solutions
