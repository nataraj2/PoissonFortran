#include <iostream>
#include <vector>

#include <petscksp.h>

int my_rank, nprocs;
int my_nx, my_ny;
int offset_i = 0, offset_j = 0;
double xmin = 0.0, xmax = 1.0, ymin = 0.0, ymax = 1.0;
PetscInt idx_start, idx_end;
std::vector<int> my_nx_vec, my_ny_vec, offset_i_vec, offset_j_vec;
std::vector<int> istart_vec, iend_vec, jstart_vec, jend_vec, idx_start_vec;
double dx, dy;

void output_solution_to_file(Vec &vec1, Vec &vec2);
int get_idx_glo(const int i_glo, const int j_glo);
PetscErrorCode MyKSPMonitorShort(KSP ksp, PetscInt its, PetscReal rnorm, void *ctx);

int main(int argc, char **args) {
  Vec x, b, b_exact; /* approx solution, RHS, exact solution */
  Mat A;             /* linear system matrix */
  KSP ksp;           /* linear solver context */
  PetscErrorCode ierr;

  ierr = PetscInitialize(NULL, NULL, NULL, NULL);CHKERRQ(ierr);

  // SETUP BEGINS

  // Setting up a problem to demonstrate for just 4 processors
  int nx = 200, ny = 200;

  MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);

  // Make the layout available for all processors. Hence all processors
  // will know the imin, imax, jmin, jmax, offset_i, offset_j etc.

  my_nx_vec.resize(nprocs);
  my_ny_vec.resize(nprocs);
  offset_i_vec.resize(nprocs);
  offset_j_vec.resize(nprocs);
  istart_vec.resize(nprocs);
  iend_vec.resize(nprocs);
  jstart_vec.resize(nprocs);
  jend_vec.resize(nprocs);
  idx_start_vec.resize(nprocs, 0.0);

  my_nx_vec[0] = 100;
  my_ny_vec[0] = 100;
  offset_i_vec[0] = 0;
  offset_j_vec[0] = 0;

  my_nx_vec[1] = 100;
  my_ny_vec[1] = 100;
  offset_i_vec[1] = 100;
  offset_j_vec[1] = 0;

  my_nx_vec[2] = 100;
  my_ny_vec[2] = 100;
  offset_i_vec[2] = 0;
  offset_j_vec[2] = 100;

  my_nx_vec[3] = 100;
  my_ny_vec[3] = 100;
  offset_i_vec[3] = 100;
  offset_j_vec[3] = 100;


  for (int rank = 0; rank <= nprocs - 1; rank++) {
    istart_vec[rank] = offset_i_vec[rank];
    iend_vec[rank] = offset_i_vec[rank] + my_nx_vec[rank] - 1;
    jstart_vec[rank] = offset_j_vec[rank];
    jend_vec[rank] = offset_j_vec[rank] + my_ny_vec[rank] - 1;
  }

  for (int rank = 0; rank <= nprocs - 1; rank++) {
    for (int p = 0; p < rank; p++) {
      idx_start_vec[rank] = idx_start_vec[rank] + my_nx_vec[p] * my_ny_vec[p];
    }
  }

  /*for (int rank = 0; rank <= nprocs - 1; rank++) {
	if(my_rank == rank){
		std::cout << "Rank is " << my_rank << "\n";
		std::cout << "istart, iend, jstart, jend " << "\n";
		std::cout << istart_vec[rank] << " " << iend_vec[rank]	<< " " << 
					 jstart_vec[rank] << " " << jend_vec[rank] << "\n";
		std::cout << "idx_start" << "\n";
		std::cout << idx_start_vec[rank] << "\n";
	}
	MPI_Barrier(PETSC_COMM_WORLD);
  }*/

  dx = (xmax - xmin) / nx;
  dy = (ymax - ymin) / ny;

  // SETUP ENDS

  ierr = MatCreate(PETSC_COMM_WORLD, &A);CHKERRQ(ierr);

  PetscInt matrix_size = nx * ny;

  my_nx = my_nx_vec[my_rank];
  my_ny = my_ny_vec[my_rank];

  ierr = MatSetSizes(A, my_nx * my_ny, my_nx * my_ny, PETSC_DECIDE, PETSC_DECIDE);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
  */
  ierr = MatGetOwnershipRange(A, &idx_start, &idx_end);CHKERRQ(ierr);

  /*for (int rank = 0; rank < nprocs; rank++) {
    if (my_rank == rank) {
      std::cout << "rank, Istart, Iend = " << my_rank << " " << idx_start << " "
                << idx_end << "\n";
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  MPI_Barrier(PETSC_COMM_WORLD);*/

  /*
     Set matrix elements for the 2-D, five-point stencil in parallel.
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly).
      - Always specify global rows and columns of matrix entries.
   */

  PetscInt row;
  PetscInt col[5];
  PetscScalar val[5];

  int i_glo, j_glo;

  offset_i = offset_i_vec[my_rank];
  offset_j = offset_j_vec[my_rank];

  for (int i = 0; i < my_nx; i++) {
    for (int j = 0; j < my_ny; j++) {

      i_glo = offset_i + i;
      j_glo = offset_j + j;

      row = get_idx_glo(i_glo, j_glo);

      if (i_glo == 0) {
		col[0] = row;
        col[1] = get_idx_glo(i_glo + 1, j_glo);
		val[0] = -1.0;
		val[1] = 1.0;
		MatSetValues(A, 1, &row, 2, col, val, INSERT_VALUES);
      }
      else if (i_glo == nx - 1) {
		col[0] = row;
        col[1] = get_idx_glo(i_glo - 1, j_glo);
		val[0] = -1.0;
		val[1] = 1.0;
		MatSetValues(A, 1, &row, 2, col, val, INSERT_VALUES);
      }
      else if (j_glo == 0) {
		col[0] = row;
        col[1] = get_idx_glo(i_glo, j_glo+1);
		val[0] = -1.0;
		val[1] = 1.0;
		MatSetValues(A, 1, &row, 2, col, val, INSERT_VALUES);
      }
      else if (j_glo == ny - 1) {
		col[0] = row;
        col[1] = get_idx_glo(i_glo, j_glo-1);
		val[0] = -1.0;
		val[1] = 1.0;
		MatSetValues(A, 1, &row, 2, col, val, INSERT_VALUES);
      }
	  else{
		col[0] = row;
		col[1] = get_idx_glo(i_glo-1, j_glo);
		col[2] = get_idx_glo(i_glo+1, j_glo);
		col[3] = get_idx_glo(i_glo, j_glo-1);
		col[4] = get_idx_glo(i_glo, j_glo+1);
		
		val[0] = -2.0*(1.0/(dx*dx)+ 1.0/(dy*dy));
		val[1] = 1.0/(dx*dx);
		val[2] = 1.0/(dx*dx);
		val[3] = 1.0/(dy*dy);
		val[4] = 1.0/(dy*dy);

      	MatSetValues(A, 1, &row, 5, col, val, INSERT_VALUES);
	  }
    }
  }

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /*MatNullSpace nullspace;

  ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
  CHKERRQ(ierr);
  ierr = MatSetNullSpace(A, nullspace);
  CHKERRQ(ierr);
  ierr = MatNullSpaceDestroy(&nullspace);
  CHKERRQ(ierr); */

  /*PetscViewer viewer;
  PetscViewerASCIIOpen(
      PETSC_COMM_WORLD, "matrix_data.m",
      &viewer); // Create a viewer to write to a file "matrix_data.m"

  PetscViewerSetFormat(viewer,
                       PETSC_VIEWER_DEFAULT); // Set the viewer format to MATLAB
  // PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); // Set the viewer
  // format to MATLAB
  MatView(A, viewer); // View the matrix in the MATLAB format

  PetscViewerDestroy(&viewer); // Destroy the viewer*/

  ierr = VecCreate(PETSC_COMM_WORLD, &b);CHKERRQ(ierr);
  ierr = VecSetSizes(b, my_nx * my_ny, PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(b);CHKERRQ(ierr);
  ierr = VecDuplicate(b, &x);CHKERRQ(ierr);
  ierr = VecDuplicate(b, &b_exact);CHKERRQ(ierr);

  /*
     Set the right-hand-side vector correspoding to the exact solution
  */

  PetscScalar value[1];
  for (int i = 0; i < my_nx; i++) {
    for (int j = 0; j < my_ny; j++) {

      i_glo = offset_i + i;
      j_glo = offset_j + j;
      row = get_idx_glo(i_glo, j_glo);

      double x = xmin + (i_glo + 0.5) * dx;
      double y = ymin + (j_glo + 0.5) * dy;

	  if(i_glo == 0 or i_glo==nx-1 or j_glo==0 or j_glo == ny-1){
		value[0] = 0.0;
	  }
	  else{
      value[0] = -4.0 * 3.1416 * 3.1416 * cos(2.0 * 3.1416 * x);//*
                 //cos(2.0 * 3.1416 * y);
	  }
      VecSetValues(b, 1, &row, value, INSERT_VALUES);
      //value[0] = cos(2.0 * 3.1416 * x) * cos(2.0 * 3.1416 * y);
      //VecSetValues(b_exact, 1, &row, value, INSERT_VALUES);
    }
  }

  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  VecAssemblyBegin(b_exact);
  VecAssemblyEnd(b_exact);

  /*ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
  CHKERRQ(ierr);
  ierr = MatNullSpaceRemove(nullspace, b);
  CHKERRQ(ierr);
  ierr = MatNullSpaceDestroy(&nullspace);
  CHKERRQ(ierr);*/

  /*PetscViewerASCIIOpen(
      PETSC_COMM_WORLD, "vector_data.m",
      &viewer); // Create a viewer to write to a file "matrix_data.m"

  // PetscViewerSetFormat(viewer, PETSC_VIEWER_DEFAULT); // Set the viewer
  // format to MATLAB
  PetscViewerSetFormat(
      viewer, PETSC_VIEWER_ASCII_MATLAB); // Set the viewer format to MATLAB
  VecView(b, viewer); // View the matrix in the MATLAB format

  PetscViewerDestroy(&viewer); // Destroy the viewer*/

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);

  ierr = KSPMonitorSet(ksp, MyKSPMonitorShort, NULL, NULL); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1000); CHKERRQ(ierr);


  ierr = KSPSolve(ksp, b, x);CHKERRQ(ierr);

  output_solution_to_file(x, b_exact);

  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_view).
  */
  ierr = PetscFinalize();
  return ierr;
}

int get_idx_glo(const int i_glo, const int j_glo) {
  // Find to which rank the i and j belong to
  int rank_curr=-1;

  for (int rank = 0; rank <= nprocs - 1; rank++) {
    if ((i_glo - istart_vec[rank]) * (i_glo - iend_vec[rank]) <= 0 and
        (j_glo - jstart_vec[rank]) * (j_glo - jend_vec[rank]) <= 0) {
      rank_curr = rank;
      break;
    }
  }

  if(rank_curr == -1){
		std::cout << "Could not find " << i_glo << " , " << j_glo << " in any rank " << "\n";
  }

  int i_loc = i_glo - offset_i_vec[rank_curr];
  int j_loc = j_glo - offset_j_vec[rank_curr];

  int idx_loc = j_loc * my_nx_vec[rank_curr] + i_loc;
  int idx_glo = idx_loc + idx_start_vec[rank_curr];

  return idx_glo;
}

void output_solution_to_file(Vec &vec1, Vec &vec2) {
  FILE *sol_file;

  std::string filename;
  PetscInt my_rank;

  MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);

  filename = "sol_file" + std::to_string(my_rank);
  filename = filename + ".txt";

  PetscScalar values[1], value_exact[1];
  PetscInt row[1];

  sol_file = fopen(filename.c_str(), "w");
  for (int i = 0; i < my_nx; i++) {
    for (int j = 0; j < my_ny; j++) {
	  int i_glo = offset_i + i;
	  int j_glo = offset_j + j;	
      row[0] = get_idx_glo(i_glo, j_glo);
      double x = xmin + (i_glo + 0.5) * dx;
      double y = ymin + (j_glo + 0.5) * dy;
      VecGetValues(vec1, 1, row, values);
      fprintf(sol_file, "%d %d %0.15g %0.15g %0.15g\n", i_glo, j_glo, x, y,
              values[0]);
    }
  }

  fclose(sol_file);
}

PetscErrorCode MyKSPMonitorShort(KSP ksp, PetscInt its, PetscReal rnorm, void *ctx) {
    PetscPrintf(PETSC_COMM_WORLD, "Iteration %D: Residual Norm %e\n", its, (double)rnorm);
    return 0;
}
