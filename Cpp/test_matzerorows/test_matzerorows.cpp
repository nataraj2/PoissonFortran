#include <petscdmda.h>
#include <petscmat.h>
#include <petscvec.h>
#include <vector>
#include <iostream>
using namespace std;

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    DM da;
    Mat A;
    Vec b;
    PetscInt mx = 5, my = 5;  // Global grid size (6x6)
    PetscInt rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Create 2D DMDA with a 5-point stencil
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                 DMDA_STENCIL_BOX, mx, my, PETSC_DECIDE, PETSC_DECIDE, 
                 1, 1, NULL, NULL, &da);
    DMSetFromOptions(da);
    DMSetUp(da);

    // Create parallel matrix and vector
    DMCreateMatrix(da, &A);
    DMCreateGlobalVector(da, &b);

    // Get local ownership range
    PetscInt xs, ys, xm, ym;  // Local grid corners and sizes
    DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);

    // Initialize matrix A and vector b
    for (PetscInt j = ys; j < ys + ym; j++) {
        for (PetscInt i = xs; i < xs + xm; i++) {
            PetscInt row = j * mx + i;  // Application ordering
            PetscScalar v = 4.0;
            MatSetValue(A, row, row, v, INSERT_VALUES);
            VecSetValue(b, row, 1.0, INSERT_VALUES);
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);


    // Get AO mapping (Application Order to PETSc order)
    AO ao;
    DMDAGetAO(da, &ao);

    // Identify boundary rows in application order
    std::vector<PetscInt> zeroedRows;
    for (PetscInt j = 0; j < my; j++) {
        for (PetscInt i = 0; i < mx; i++) {
            if (i == 0 || i == mx - 1 || j == 0 || j == my - 1) { // Boundary points
                PetscInt row = j * mx + i;
                zeroedRows.push_back(row);
            }
        }
    }

    // Convert to PETSc ordering
    AOApplicationToPetsc(ao, zeroedRows.size(), zeroedRows.data());

    // Apply MatZeroRows
    MatZeroRows(A, zeroedRows.size(), zeroedRows.data(), 1.0, NULL, NULL);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    // View the modified system
    PetscPrintf(PETSC_COMM_WORLD, "Modified Matrix A:\n");
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Modified Vector b:\n");
    VecView(b, PETSC_VIEWER_STDOUT_WORLD);

    // Cleanup
    MatDestroy(&A);
    VecDestroy(&b);
    DMDestroy(&da);

    PetscFinalize();
    return 0;
}
