program PoissonFortran2D
#include <petsc/finclude/petscksp.h>
	use petscksp	
	use petscdmda
	use ModuleVariables
	implicit none
	PetscErrorCode ierr
	DM dm
	PetscInt  xm,ym,xs,ys, mx, my
	PetscInt  nx, ny, i, j
	double precision xval, yval, Hx, Hy
	PetscInt  rank, ierror
	character(len=10) :: file_id
	character(len=50) :: file_name

	! Global dimensions of the matrix
	nx = 400
	ny = 200

	! Petsc stuff to initialize the parallel blocks

	call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
	call DMDACreate2D(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC,DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR,nx,ny,&
	 &			 PETSC_DECIDE,PETSC_DECIDE,1,1, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, dm, ierr)
	call DMSetUp(dm,ierr)
	! Get the block limits on this processor
	call DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)

	call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierror)
	print*, "Rank ", rank, "xlo, xhi " , xs, xs+xm-1, " ylo, yhi = ", ys, ys+ym-1

	! Solving on a domain (0,1) x (0,1)
	Hx = 1.0/nx
	Hy = 1.0/ny
	
	! Allocate the rhs and pressure on this processor with appropriate index limits
	allocate(poisson_rhs(xs:xs+xm-1,ys:ys+ym-1))
	allocate(pressure(xs:xs+xm-1,ys:ys+ym-1))

	do j=ys,ys+ym-1
		do i=xs,xs+xm-1
			xval = (i+0.5)*Hx
			yval = (j+0.5)*Hy
			poisson_rhs(i,j) = 8.0*3.1415*3.1415*cos(2.0*3.1415*xval)*cos(2.0*3.1415*yval)*Hx*Hy;
	   	end do
	end do

	! Call Petsc to solve - writes the solution to "pressure"
	call solve_poisson(dm)

	! Write out solution
	write(file_id, '(i0)') rank
	file_name = 'sol_file_' // trim(adjustl(file_id)) // '.txt'
	open(unit=10,file = trim(file_name), status="new")

	do j=ys,ys+ym-1
	   do i=xs,xs+xm-1
			write(10,*)i, j, (i+0.5)*Hx, (j+0.5)*Hy, pressure(i,j)
	   end do
	end do

	close(10)

	deallocate(poisson_rhs)
	deallocate(pressure)
	
	call DMDestroy(dm,ierr)
	call PetscFinalize(ierr)	

end program PoissonFortran2D
