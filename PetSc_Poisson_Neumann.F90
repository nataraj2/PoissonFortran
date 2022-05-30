subroutine solve_poisson(dm)
#include <petsc/finclude/petscksp.h>
	use ModuleVariables
    use petscdmda
    use petscksp
    implicit none

    PetscInt is,js,iw,jw
    PetscInt one,three
    PetscErrorCode ierr
    KSP ksp
    DM dm
    Vec x
    external ComputeRHS,ComputeMatrix,ComputeInitialGuess, OutputSolutionToFile

    one = 1
    three = 3

    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
    call DMSetFromOptions(dm,ierr)
    call KSPSetDM(ksp,dm,ierr)
    call KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,0,ierr)
    call KSPSetComputeRHS(ksp,ComputeRHS,0,ierr)
    call KSPSetComputeOperators(ksp,ComputeMatrix,0,ierr)
    call DMDAGetCorners(dm,is,js,PETSC_NULL_INTEGER,iw,jw,PETSC_NULL_INTEGER,ierr)
    call KSPSetFromOptions(ksp,ierr)
    call KSPSetUp(ksp,ierr)
    call KSPSolve(ksp,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr)

    call KSPGetSolution(ksp,x,ierr)
    call OutputSolutionToFile(ksp,x);

    call KSPDestroy(ksp,ierr)
    end subroutine solve_poisson

    subroutine OutputSolutionToFile(ksp,x)
	use ModuleVariables
    use petscksp
    use petscdmda
    implicit none
    PetscErrorCode ierr
    KSP ksp
    Vec x
    DM           dm
    PetscScalar  Hx, Hy
    PetscInt     i,j,mx,my,xm,ym,xs,ys
    PetscScalar, pointer :: xx(:,:)

	integer :: m
	integer :: fu
	character(len=10) :: file_id
	character(len=50) :: file_name
    PetscInt :: rank, ierror

	call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierror)

	write(file_id, '(i0)') rank

	file_name = 'sol_file_' // trim(adjustl(file_id)) // '.txt'

    call KSPGetDM(ksp,dm,ierr)
    call DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
    call DMDAVecGetArrayF90(dm,x,xx,ierr);CHKERRQ(ierr)

    print*, "xs, ys = " , xs, ys
    print*, "xe, ye = " , xs+xm-1, ys+ym-1 
 
    do j=ys,ys+ym-1
       do i=xs,xs+xm-1
			pressure(i,j) = xx(i,j)
       end do
    end do

        
    end subroutine OutputSolutionToFile


    subroutine ComputeInitialGuess(ksp,b,ctx,ierr)
    use petscksp
    implicit none
    PetscErrorCode  ierr
    KSP ksp
    PetscInt ctx(*)
    Vec b
    PetscScalar  h

    h=0.0
    call VecSet(b,h,ierr)
    end subroutine

    subroutine ComputeRHS(ksp,b,dummy,ierr)
    use petscksp
    use petscdmda
    implicit none

    PetscErrorCode  ierr
    KSP ksp
    Vec b
    integer dummy(*)
    PetscScalar  h,Hx,Hy
    PetscInt  mx,my, xs, ys, xm, ym, i, j
    DM dm
	PetscScalar xval, yval;
	PetscScalar, pointer :: xx(:,:)
	MatNullSpace nullspace;
	Vec dummyVecs(1)

    call KSPGetDM(ksp,dm,ierr)
    call DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,          &
     &            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,          &
     &            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)

    Hx = 1.0 / real(mx)
    Hy = 1.0 / real(my)
    h = Hx*Hy
	!call ComputePoissonRHS(check_rhs)

    call DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
    call DMDAVecGetArrayF90(dm,b,xx,ierr);CHKERRQ(ierr)

	do j=ys,ys+ym-1
       do i=xs,xs+xm-1
			xval = (i+0.5)*Hx
			yval = (j+0.5)*Hy	
			!if(i.eq.0 .and. j.eq. 0)then
			!	xx(i,j) = cos(2.0*3.1415*xval)*cos(2.0*3.1415*yval)
			!else
				xx(i,j) = 8.0*3.1415*3.1415*cos(2.0*3.1415*xval)*cos(2.0*3.1415*yval)*Hx*Hy;
			!end if
			
       end do
    end do
	
	call DMDAVecRestoreArrayF90(dm,b,xx,ierr);CHKERRQ(ierr)
	call VecAssemblyBegin(b, ierr);CHKERRQ(ierr);
    call VecAssemblyEnd(b, ierr);CHKERRQ(ierr);


    call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,dummyVecs,nullspace, ierr);
    call MatNullSpaceRemove(nullspace,b, ierr);
    call MatNullSpaceDestroy(nullspace, ierr);

    end subroutine

      subroutine ComputeMatrix(ksp,A,B,dummy,ierr)
      use petscksp
    implicit none
    PetscErrorCode  ierr
    KSP ksp
    Mat A,B
    integer dummy(*)
    DM dm

      PetscInt    i,j,mx,my,xm
      PetscInt    ym,xs,ys,i1,i5, num, numi, numj
      PetscScalar  v(5),Hx,Hy
      PetscScalar  HxdHy,HydHx
      MatStencil   row(4),col(4,5)
	  MatNullSpace nullspace;
	  Vec dummyVecs(1)

      i1 = 1
      i5 = 5
      call KSPGetDM(ksp,dm,ierr)
      call DMDAGetInfo(dm,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,          &
     &           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,          &
     &           PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)

      Hx = 1.0 / real(mx)
      Hy = 1.0 / real(my)
      HxdHy = Hx/Hy
      HydHx = Hy/Hx
      call DMDAGetCorners(dm,xs,ys,PETSC_NULL_INTEGER,xm,ym,PETSC_NULL_INTEGER,ierr)
   do 10,j=ys,ys+ym-1
     do 20,i=xs,xs+xm-1
       row(MatStencil_i) = i
       row(MatStencil_j) = j
	
	   !if(i.eq.0 .and. j.eq.0)then
		!v(1) = 1.0;
	    !col(MatStencil_i,1) = i
        !col(MatStencil_j,1) = j
		!call MatSetValuesStencil(B,1,row,1,col,v,INSERT_VALUES, ierr)
	   	

       !else if (i.eq.0 .or. j.eq.0 .or. i.eq.mx-1 .or. j.eq.my-1 .and. (i+j .ne. 0)) then
       if (i.eq.0 .or. j.eq.0 .or. i.eq.mx-1 .or. j.eq.my-1) then

		  num = 1; numi=0; numj=0;
          if (j.ne.0) then
            v(num) = -HxdHy
            col(MatStencil_i,num) = i
            col(MatStencil_j,num) = j-1
            num = num+1
			numj = numj+1
		  end if

          if (i.ne.0) then
            v(num)   = -HydHx
            col(MatStencil_i,num) = i-1
            col(MatStencil_j,num) = j
            num = num+1
			numi = numi+1
		  end if
          
          if (i.ne.mx-1) then
            v(num)   = -HydHx
            col(MatStencil_i,num) = i+1
            col(MatStencil_j,num) = j
            num = num+1
			numi = numi+1
		  end if
          
          if (j.ne.my-1) then
            v(num)   = -HxdHy
            col(MatStencil_i,num) = i
            col(MatStencil_j,num) = j+1
            num = num + 1;
			numj = numj + 1;
		  end if
          
          v(num) = numj*HxdHy + numi*HydHx
		  col(MatStencil_i,num) = i
		  col(MatStencil_j,num) = j
          !num = num+1
          call MatSetValuesStencil(B,1,row,num,col,v,INSERT_VALUES, ierr)

       else
         v(1) = -HxdHy
         col(MatStencil_i,1) = i
         col(MatStencil_j,1) = j-1
         v(2) = -HydHx
         col(MatStencil_i,2) = i-1
         col(MatStencil_j,2) = j
         v(3) = 2.0*(HxdHy + HydHx)
         col(MatStencil_i,3) = i
         col(MatStencil_j,3) = j
         v(4) = -HydHx
         col(MatStencil_i,4) = i+1
         col(MatStencil_j,4) = j
         v(5) = -HxdHy
         col(MatStencil_i,5) = i
         col(MatStencil_j,5) = j+1
         call MatSetValuesStencil(B,i1,row,i5,col,v,INSERT_VALUES,ierr)
         endif
 20      continue
 10   continue
    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)


    call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,dummyVecs,nullspace,ierr)
    call MatSetNullSpace(A,nullspace,ierr);
    call MatNullSpaceDestroy(nullspace,ierr);

    if (A .ne. B) then
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    endif
    end subroutine

!/*TEST
!
!   test:
!      nsize: 4
!      args: -ksp_monitor_short -da_refine 5 -pc_type mg -pc_mg_levels 5 -mg_levels_ksp_type chebyshev -mg_levels_ksp_max_it 2 -mg_levels_pc_type jacobi -ksp_pc_side right
!
!TEST*/
