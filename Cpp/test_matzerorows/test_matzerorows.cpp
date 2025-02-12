#include <petscdmda.h>
#include <petscmat.h>
#include <petscvec.h>
using namespace std;

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    DM da;
    Mat A;
    Vec b;
    PetscInt mx = 8, my = 8;  // Global grid size (6x6)
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

    // Initialize A and b
    PetscInt i, j, row;
    PetscScalar v;
    for (j = ys; j < ys + ym; j++) {
        for (i = xs; i < xs + xm; i++) {
            row = j * mx + i;  // Convert (i, j) to global row index
            v = 4.0;
            MatSetValue(A, row, row, v, INSERT_VALUES);
            VecSetValue(b, row, 1.0, INSERT_VALUES); // RHS initialization
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    // Identify boundary rows to zero
    PetscInt zeroedRows[mx * 2 + my * 2 - 4];  // Perimeter points
    //std::vector<int> ibkRows;
    PetscInt count = 0;

    for (j = 0; j < my; j++) {
        for (i = 0; i < mx; i++) {
            if (i == 0 || i == mx - 1 || j == 0 || j == my - 1) { // Boundary points
                row = j * mx + i;
                if (row >= xs * mx + ys && row < (xs + xm) * mx + (ys + ym)) { // Local check
                    zeroedRows[count++] = row;
                }
            }
        }
    }

    // Apply MatZeroRows in parallel
    MatZeroRows(A, count, zeroedRows, 1.0, NULL, b);

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
