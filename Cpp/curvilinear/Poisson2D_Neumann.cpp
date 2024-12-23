/*T
   Concepts: KSP^solving a system of linear equations
   Concepts: KSP^Laplacian, 2d
   Processors: n
T*/

/*
Laplacian in 2D. Modeled by the partial differential equation

   div  grad u = -f,  0 < x,y < 1,

with forcing function

   f = 8 * pi**2 * cos(2*pi*x) * cos(2*pi*y) 

with pure Neumann boundary conditions

The functions are cell-centered

This uses multigrid to solve the linear system

	   Contributed by Andrei Draganescu <aidraga@sandia.gov>

Note the nice multigrid convergence despite the fact it is only using
peicewise constant interpolation/restriction. This is because cell-centered multigrid
does not need the same rule:

	polynomial degree(interpolation) + polynomial degree(restriction) + 2 > degree of PDE

that vertex based multigrid needs.
*/

static char help[] = "Solves 2D inhomogeneous Laplacian using multigrid.\n\n";
#include "header.h"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <iostream>
#include <vector>
using namespace std;

vector<vector<double>> rhs;
extern PetscErrorCode ComputeMatrix(KSP,Mat,Mat,void*);
extern PetscErrorCode ComputeRHS(KSP,Vec,void*);
void apply_face_dpdn(int k,int nxc,int nyc,MatStencil &row,MatStencil col[],PetscScalar v[],double** kk_f,double** xi_c,double** eta_c,double** inv_Jac,int face);
void output_solution_to_file(KSP &ksp, Vec &x);
void output_p3d_file(KSP &ksp, Vec &x,double** Xc, double** Yc);

typedef enum {DIRICHLET, NEUMANN} BCType;

typedef struct {
  PetscScalar nu;
  BCType	  bcType;
} UserContext;

double** xi_xc;
double** xi_yc;
double** eta_xc;
double** eta_yc;
double** inv_Jac;
double** xi_xf;
double** xi_yf;
double** eta_xf;
double** eta_yf;

int btm = 1; int lft = 2; int tp = 3; int rght = 4;
int Ncols = 6;
int main(int argc,char **argv)
{
  KSP			 ksp;
  DM			 da;	// NS: data structure for managing data on a structured grid
  UserContext	 user;
  const char	 *bcTypes[2] = {"dirichlet","neumann"};
  PetscErrorCode ierr;
  PetscInt	     bc;
  Vec 			 x;
  PetscInt       xm,ym,xs,ys;

  PetscInt nx = 51;
  PetscInt ny = 51;
  PetscInt lx[] = {40,60};
  PetscInt ly[] = {40,60};
  
  /////////////////// GENERATE GRID ////////////////////////////////
  double** X = Create2DMatrix(nx,ny);
  double** Y = Create2DMatrix(nx,ny);
   
  get_curv_grid(X,Y,nx,ny);
  
  // Cell-center coordinates
  double** Xc = Create2DMatrix(nx-1,ny-1);
  double** Yc = Create2DMatrix(nx-1,ny-1);
  for(int i=0;i<nx-1;i++){
	for(int j=0;j<ny-1;j++)
	{
		Xc[i][j] = 0.25 * (X[i][j]+X[i+1][j]+X[i][j+1]+X[i+1][j+1]);
		Yc[i][j] = 0.25 * (Y[i][j]+Y[i+1][j]+Y[i][j+1]+Y[i+1][j+1]);
	}}
  
  // Face/Edge-center coordinates	  
  double** Xfx = Create2DMatrix(nx,ny-1);
  double** Yfx = Create2DMatrix(nx,ny-1);
  for(int i=0;i<nx;i++){
	for(int j=0;j<ny-1;j++)
	{
		Xfx[i][j] = 0.5 * (X[i][j]+X[i][j+1]);
		Yfx[i][j] = 0.5 * (Y[i][j]+Y[i][j+1]);
	}}
	
  double** Xfy = Create2DMatrix(nx-1,ny);
  double** Yfy = Create2DMatrix(nx-1,ny);
  for(int i=0;i<nx-1;i++){
    for(int j=0;j<ny;j++)
    {
   	    Xfy[i][j] = 0.5 * (X[i][j]+X[i+1][j]);
  		Yfy[i][j] = 0.5 * (Y[i][j]+Y[i+1][j]);
  	}}
	
  //////////////////////  METRIC TERMS //////////////////////////////
  double** xc_xi = differentiate_fn(Xc,nx-1,ny-1,1,0);
  double** xc_eta = differentiate_fn(Xc,nx-1,ny-1,2,0);
  double** yc_xi = differentiate_fn(Yc,nx-1,ny-1,1,0);
  double** yc_eta = differentiate_fn(Yc,nx-1,ny-1,2,0);
  
  inv_Jac = Create2DMatrix(nx-1,ny-1);
  double** Jac = Create2DMatrix(nx-1,ny-1);
  xi_xc = Create2DMatrix(nx-1,ny-1);
  xi_yc = Create2DMatrix(nx-1,ny-1);
  eta_xc = Create2DMatrix(nx-1,ny-1);
  eta_yc = Create2DMatrix(nx-1,ny-1);
  for(int i=0;i<nx-1;i++){	
	for(int j=0;j<ny-1;j++)
	{
		inv_Jac[i][j] = xc_xi[i][j]*yc_eta[i][j] - xc_eta[i][j]*yc_xi[i][j];
		Jac[i][j] = 1.0/inv_Jac[i][j];
		xi_xc[i][j] = Jac[i][j]*yc_eta[i][j];
		xi_yc[i][j] = -Jac[i][j]*xc_eta[i][j];
		eta_xc[i][j] = -Jac[i][j]*yc_xi[i][j];
		eta_yc[i][j] = Jac[i][j]*xc_xi[i][j];
	}}

  xi_xf = Create2DMatrix(nx,ny-1);
  xi_yf = Create2DMatrix(nx,ny-1);
  eta_xf = Create2DMatrix(nx-1,ny);
  eta_yf = Create2DMatrix(nx-1,ny);
  get_face_metric_terms(xi_xf,xi_yf,eta_xf,eta_yf,Xfx,Yfx,Xfy,Yfy,nx,ny);
  
  ///////////////////////////////////////////////////////////////////

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr); // Create a context used to to solve a linear system
  
  // NS: Create the grid
  ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,DMDA_STENCIL_BOX,nx-1,ny-1,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,&da);CHKERRQ(ierr);
  //ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,nx,ny,2,2,1,1,lx,ly,&da);CHKERRQ(ierr);
  
  ierr = DMSetFromOptions(da);CHKERRQ(ierr); // NS: Set da using the command-line arguments
  ierr = DMSetUp(da);CHKERRQ(ierr);			// NS: Setup the da object
  ierr = DMDASetInterpolationType(da, DMDA_Q0);CHKERRQ(ierr);

  ierr = KSPSetDM(ksp,da);CHKERRQ(ierr);

  PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for the inhomogeneous Poisson equation", "DM");
  bc		  = (PetscInt)NEUMANN;
  user.bcType = (BCType)bc;
  PetscOptionsEnd();

  //PetscScalar Hx   = 1.0 / (PetscReal)(nx);
  //PetscScalar Hy   = 1.0 / (PetscReal)(ny);
  ierr = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  //ierr = DMDAGetCorners(da,&ys,&xs,0,&ym,&xm,0);CHKERRQ(ierr);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  cout << "Rank is " << rank <<  " xs = " << xs << " ys = " << ys << " xm = " << xm << " ym = " << ym << endl;
  // cout << "div1 = " << 4/5 <<  " div2 = " << 9/5 << endl;

  rhs.resize(xm,vector<double>(ym));

  /*for (int i=0; i<xm; i++) {
    for (int j=0; j<ym; j++) {
        PetscReal x = ((PetscReal)(xs+i)+0.5)*Hx;
        PetscReal y = ((PetscReal)(ys+j)+0.5)*Hy;
        rhs[i][j] = 8.0*3.1416*3.1416*cos(2.0*3.1416*x)*cos(2.0*3.1416*y)*Hx*Hy;
    }
  }*/
  
  for (int k=0; k<(nx-1)*(ny-1); k++) {
	  int j = k/(nx-1);
	  int i = k - j*(nx-1);
	  rhs[i][j] = -8.0*pi*pi*cos(2.0*pi*Xc[i][j])*cos(2.0*pi*Yc[i][j]) * inv_Jac[i][j];
  }
  rhs[0][0] = cos(2.0*pi*Xc[0][0])*cos(2.0*pi*Yc[0][0]);	// To fix a point value for avoiding singularity with periodic BCs
	  
  ierr = KSPSetComputeRHS(ksp,ComputeRHS,&user);CHKERRQ(ierr);
  ierr = KSPSetComputeOperators(ksp,ComputeMatrix,&user);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);	// NS: Set the algorithm, preconditioner, and the associated parameters, using the command-line arguments
  ierr = KSPSolve(ksp,NULL,NULL);CHKERRQ(ierr);	// NS: Solve the system of linear equations
  ierr = KSPGetSolution(ksp,&x);CHKERRQ(ierr);

  //output_solution_to_file(ksp,x);
  output_p3d_file(ksp,x,Xc,Yc);
  
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);	// NS: Free the KSP context and all storage associated with it
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = PetscFinalize();
  
  //////////////////////  DEALLOCATE MEMORY //////////////////////////////
  Delete2DMatrix(X,nx,ny);
  Delete2DMatrix(Y,nx,ny);
  Delete2DMatrix(Xc,nx-1,ny-1);
  Delete2DMatrix(Yc,nx-1,ny-1);
  Delete2DMatrix(Xfx,nx,ny-1);
  Delete2DMatrix(Yfx,nx,ny-1);
  Delete2DMatrix(Xfy,nx-1,ny);
  Delete2DMatrix(Yfy,nx-1,ny);
  
  Delete2DMatrix(xc_xi,nx-1,ny-1);
  Delete2DMatrix(xc_eta,nx-1,ny-1);
  Delete2DMatrix(yc_xi,nx-1,ny-1);
  Delete2DMatrix(yc_eta,nx-1,ny-1);
  Delete2DMatrix(inv_Jac,nx-1,ny-1);
  Delete2DMatrix(Jac,nx-1,ny-1);
  Delete2DMatrix(xi_xc,nx-1,ny-1);
  Delete2DMatrix(xi_yc,nx-1,ny-1);
  Delete2DMatrix(eta_xc,nx-1,ny-1);
  Delete2DMatrix(eta_yc,nx-1,ny-1);
  Delete2DMatrix(xi_xf,nx,ny-1);
  Delete2DMatrix(xi_yf,nx,ny-1);
  Delete2DMatrix(eta_xf,nx-1,ny);
  Delete2DMatrix(eta_yf,nx-1,ny);
  
  ///////////////////////////////////////////////////////////////////

  return ierr;
}

PetscErrorCode ComputeRHS(KSP ksp,Vec b,void *ctx)
{
  UserContext   *user = (UserContext*)ctx;
  PetscErrorCode ierr;
  PetscInt     i,j,mx,my,xm,ym,xs,ys,nxc,nyc;
  PetscScalar   Hx,Hy;
  PetscScalar   **array;
  DM             da;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&da);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da, 0, &mx, &my, 0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  
  cout << "In ComputeRHS: mx = " << mx <<  " my = " << my << endl;
  nxc = mx;
  nyc = my;
  
  Hx   = 1.0 / (PetscReal)(mx);
  Hy   = 1.0 / (PetscReal)(my);
  ierr = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, b, &array);CHKERRQ(ierr);
  /*for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
         array[j][i] = rhs[i-xs][j-ys];
    }
  }*/
  for(int i=0;i<nxc;i++){	
	for(int j=0;j<nyc;j++)
	{
		array[j][i] = rhs[i][j];
	}}
  ierr = DMDAVecRestoreArray(da, b, &array);CHKERRQ(ierr);  
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  /* force right hand side to be consistent for singular matrix */
  /* note this is really a hack, normally the model would provide you with a consistent right handside */
  /*if (user->bcType == NEUMANN) {
    MatNullSpace nullspace;

    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
    ierr = MatNullSpaceRemove(nullspace,b);CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
  }*/
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeMatrix(KSP ksp, Mat J,Mat jac, void *ctx)
{
  UserContext	*user = (UserContext*)ctx;
  PetscErrorCode ierr;
  PetscInt	   i,j,mx,my,xm,ym,xs,ys,num, numi, numj,nxc,nyc;
  PetscScalar	v[Ncols];
  MatStencil	 row, col[Ncols];
  DM			 da;
  
  PetscFunctionBeginUser;
  ierr  = KSPGetDM(ksp,&da);CHKERRQ(ierr);
  ierr  = DMDAGetInfo(da,0,&mx,&my,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr  = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
	
  cout << "In ComputeMatrix: mx = " << mx <<  " my = " << my << endl;
  nxc = mx;
  nyc = my;
  
  for (int k=0; k<nxc*nyc; k++) 
  {
	  apply_face_dpdn(k,nxc,nyc,row,col,v,xi_xf,xi_xc,eta_xc,inv_Jac,rght);
  	  //cout << "row.i = " << row.i << "   row.j = " << row.j << "   col[0].i = " << col[0].i << "   col[0].j = " << col[0].j << endl;
	  ierr = MatSetValuesStencil(jac,1,&row,Ncols,col,v,ADD_VALUES);CHKERRQ(ierr);
	    
	  apply_face_dpdn(k,nxc,nyc,row,col,v,xi_xf,xi_xc,eta_xc,inv_Jac,lft);
	  ierr = MatSetValuesStencil(jac,1,&row,Ncols,col,v,ADD_VALUES);CHKERRQ(ierr);
	  
	  apply_face_dpdn(k,nxc,nyc,row,col,v,xi_yf,xi_yc,eta_yc,inv_Jac,rght); 
	  ierr = MatSetValuesStencil(jac,1,&row,Ncols,col,v,ADD_VALUES);CHKERRQ(ierr);
	  
	  apply_face_dpdn(k,nxc,nyc,row,col,v,xi_yf,xi_yc,eta_yc,inv_Jac,lft);
	  ierr = MatSetValuesStencil(jac,1,&row,Ncols,col,v,ADD_VALUES);CHKERRQ(ierr);
	  
	  apply_face_dpdn(k,nxc,nyc,row,col,v,eta_xf,eta_xc,xi_xc,inv_Jac,tp);
	  ierr = MatSetValuesStencil(jac,1,&row,Ncols,col,v,ADD_VALUES);CHKERRQ(ierr);
	  
	  apply_face_dpdn(k,nxc,nyc,row,col,v,eta_xf,eta_xc,xi_xc,inv_Jac,btm);
	  ierr = MatSetValuesStencil(jac,1,&row,Ncols,col,v,ADD_VALUES);CHKERRQ(ierr);
	  
	  apply_face_dpdn(k,nxc,nyc,row,col,v,eta_yf,eta_yc,xi_yc,inv_Jac,tp);
	  ierr = MatSetValuesStencil(jac,1,&row,Ncols,col,v,ADD_VALUES);CHKERRQ(ierr);
	  
	  apply_face_dpdn(k,nxc,nyc,row,col,v,eta_yf,eta_yc,xi_yc,inv_Jac,btm);
	  ierr = MatSetValuesStencil(jac,1,&row,Ncols,col,v,ADD_VALUES);CHKERRQ(ierr);
  }
    
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  PetscInt fixedIdx = 0;
  ierr = MatZeroRows(jac, 1, &fixedIdx, 1.0, NULL, NULL);CHKERRQ(ierr);	  // To fix a point value for avoiding singularity with periodic BCs
  
  /*if (user->bcType == NEUMANN) {
	MatNullSpace nullspace;

	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
	ierr = MatSetNullSpace(J,nullspace);CHKERRQ(ierr);
	ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
  }*/
  PetscFunctionReturn(0);
}

void apply_face_dpdn(int k,int nxc,int nyc,MatStencil &row,MatStencil col[],PetscScalar v[],double** kk_f,double** xi_c,double** eta_c,double** inv_Jac,int face)
{
  	int j = floor(k/nxc);
  	int i = k - j*nxc;
  
	// Set indices for face normal derivative calculation
	int iface, jface, nml;
	if (face==rght){
		iface = i+1;    jface = j;
		nml = 1;
	}
	else if (face==lft){
		iface = i;    jface = j;
	    nml = -1;
	}
	else if (face==tp){
		iface = i;    jface = j+1;
	    nml = 1;
	}
	else if (face==btm){
		iface = i;    jface = j;
		nml = -1;
	}
	
	// Set indices for face tangential derivative calculation
	int iout, iin, jout, jin, it[4], jt[4];
	if (face==rght || face==lft){
		iout = i+nml;     jout = j;
	    iin = i;     	jin = j;		
		/// MODIFY INDEX FOR PERIODIC BOUNDARY ///
		if (iout<0){
		    iout = nxc + iout;
		}
		else if (iout>(nxc-1)){
			iout = iout - nxc;
		}
		//////////////////////////////////////////
		
		it[0] = iout;     jt[0] = jout+1;
	    it[1] = iin;      jt[1] = jin+1;
	    it[2] = iout;     jt[2] = jout-1;
	    it[3] = iin;      jt[3] = jin-1;
	}
	else if (face==tp || face==btm){
		iout = i;     jout = j+nml;
	    iin = i;     jin = j;		
		/// MODIFY INDEX FOR PERIODIC BOUNDARY ///
		if (jout<0){
			jout = nyc + jout;
		}
		else if (jout>(nyc-1)){
		    jout = jout - nyc;
		}
		//////////////////////////////////////////
		
		it[0] = iout+1;     jt[0] = jout;
	    it[1] = iin+1;      jt[1] = jin;
	    it[2] = iout-1;     jt[2] = jout;
	    it[3] = iin-1;      jt[3] = jin;
	}
	
	/// MODIFY INDEX FOR PERIODIC BOUNDARY ///
    for (int idx=0; idx<4; idx++){
		if (it[idx]<0){
			it[idx] = nxc + it[idx];
		}
		else if (it[idx]>(nxc-1)){
			it[idx] = it[idx] - nxc;
		}
		
		if (jt[idx]<0){
			jt[idx] = nyc + jt[idx];
		}
		else if (jt[idx]>(nyc-1)){
			jt[idx] = jt[idx] - nyc;
		}
	}
	//////////////////////////////////////////
	
	/*cout << "In apply_face_dpdn: k = " << k << "   row.i = " << i << "   row.j = " << j << "   col[0].i = " << iout << "   col[0].j = " << jout <<  
		"   col[1].i = " << iin << "   col[1].j = " << jin << endl;
	cout << "   col[2].i = " << it[0] << "   col[2].j = " << jt[0] << "   col[3].i = " << it[1] << "   col[3].j = " << jt[1] 
		 << "   col[4].i = " << it[2] << "   col[4].j = " << jt[2] << "   col[5].i = " << it[3] << "   col[5].j = " << jt[3] << endl;*/
	
	row.i = i; row.j = j;
	// normal derivative on the face
	col[0].i = iout;   col[0].j = jout;		v[0] = kk_f[iface][jface] * xi_c[iout][jout] * inv_Jac[iout][jout];
 	col[1].i = iin;   col[1].j = jin;		v[1] = - kk_f[iface][jface] * xi_c[iin][jin] * inv_Jac[iin][jin];
	
    // tangential derivative on the face
 	col[2].i = it[0];   col[2].j = jt[0];		v[2] = nml*0.25*kk_f[iface][jface] * eta_c[it[0]][jt[0]] * inv_Jac[it[0]][jt[0]];
 	col[3].i = it[1];   col[3].j = jt[1];		v[3] = nml*0.25*kk_f[iface][jface] * eta_c[it[1]][jt[1]] * inv_Jac[it[1]][jt[1]];
 	col[4].i = it[2];   col[4].j = jt[2];		v[4] = - nml*0.25*kk_f[iface][jface] * eta_c[it[2]][jt[2]] * inv_Jac[it[2]][jt[2]];
 	col[5].i = it[3];   col[5].j = jt[3];		v[5] = - nml*0.25*kk_f[iface][jface] * eta_c[it[3]][jt[3]] * inv_Jac[it[3]][jt[3]];
}

PetscErrorCode ComputeMatrix_old(KSP ksp, Mat J,Mat jac, void *ctx)
{
  UserContext	*user = (UserContext*)ctx;
  PetscErrorCode ierr;
  PetscInt	   i,j,mx,my,xm,ym,xs,ys,num, numi, numj;
  PetscScalar	v[5],Hx,Hy,HydHx,HxdHy;
  MatStencil	 row, col[5];
  DM			 da;

  PetscFunctionBeginUser;
  ierr  = KSPGetDM(ksp,&da);CHKERRQ(ierr);
  ierr  = DMDAGetInfo(da,0,&mx,&my,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  Hx	= 1.0 / (PetscReal)(mx);
  Hy	= 1.0 / (PetscReal)(my);
  HxdHy = Hx/Hy;
  HydHx = Hy/Hx;
  ierr  = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
	for (i=xs; i<xs+xm; i++) {
	  row.i = i; row.j = j;
	  if ((i==0 || j==0 || i==mx-1 || j==my-1)) {
		if (user->bcType == DIRICHLET) {
		  v[0] = 2.0*(HxdHy + HydHx);
		  ierr = MatSetValuesStencil(jac,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
		  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Dirichlet boundary conditions not supported !\n");
		} else if (user->bcType == NEUMANN) {
		  num = 0; numi=0; numj=0;
		  if (j!=0) {
			v[num] = -HxdHy;
			col[num].i = i;
			col[num].j = j-1;
			num++; numj++;
		  }
		  if (i!=0) {
			v[num]	 = -HydHx;
			col[num].i = i-1;
			col[num].j = j;
			num++; numi++;
		  }
		  if (i!=mx-1) {
			v[num]	 = -HydHx;
			col[num].i = i+1;
			col[num].j = j;
			num++; numi++;
		  }
		  if (j!=my-1) {
			v[num]	 = -HxdHy;
			col[num].i = i;
			col[num].j = j+1;
			num++; numj++;
		  }
		  v[num] = (PetscReal)(numj)*HxdHy + (PetscReal)(numi)*HydHx; col[num].i = i;   col[num].j = j;
		  num++;
		  ierr = MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);CHKERRQ(ierr);
		}
	  } else {
		v[0] = -HxdHy;			  col[0].i = i;   col[0].j = j-1;
		v[1] = -HydHx;			  col[1].i = i-1; col[1].j = j;
		v[2] = 2.0*(HxdHy + HydHx); col[2].i = i;   col[2].j = j;
		v[3] = -HydHx;			  col[3].i = i+1; col[3].j = j;
		v[4] = -HxdHy;			  col[4].i = i;   col[4].j = j+1;
		ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
	  }
	}
  }
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (user->bcType == NEUMANN) {
	MatNullSpace nullspace;

	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
	ierr = MatSetNullSpace(J,nullspace);CHKERRQ(ierr);
	ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

void output_solution_to_file(KSP &ksp, Vec &x)
{
	PetscErrorCode ierr;
	DM			 dm;
	PetscScalar	**barray;
	PetscScalar	Hx, Hy;
	PetscInt	   i,j,mx,my,xm,ym,xs,ys;
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	PetscFunctionBeginUser;
	ierr  = KSPGetDM(ksp,&dm);
	ierr  = DMDAGetInfo(dm, 0, &mx, &my, 0,0,0,0,0,0,0,0,0,0);
	Hx	= 1.0 / (PetscReal)(mx);
	Hy	= 1.0 / (PetscReal)(my);
	ierr  = DMDAGetCorners(dm,&xs,&ys,0,&xm,&ym,0);
	ierr  = DMDAVecGetArray(dm,x,&barray);

	std::cout << "Rank is " << rank << " " << xs << " " << xs+xm << " " << ys << " " << ys+ym << "\n";	

	FILE *sol_file;
	
	std::string filename;
	
	filename = "sol_file" + std::to_string(rank);
	filename = filename + ".txt";

	sol_file = fopen(filename.c_str(),"w");

	for (j=ys; j<ys+ym; j++) {
	   for (i=xs; i<xs+xm; i++) {
		   fprintf(sol_file, "%d %d %0.15g %0.15g %0.15g\n", i, j, ((PetscReal)i+0.5)*Hx, ((PetscReal)j+0.5)*Hy, barray[j][i]);
	   }
	}

	fclose(sol_file);

	DMDAVecRestoreArray(dm,x,&barray);


}

void output_p3d_file(KSP &ksp, Vec &x,double** Xc, double** Yc)
{
	PetscErrorCode ierr;
	DM			 dm;
	PetscScalar	**barray;
	PetscScalar	Hx, Hy;
	PetscInt	   i,j,mx,my,xm,ym,xs,ys;
	int rank;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	PetscFunctionBeginUser;
	ierr  = KSPGetDM(ksp,&dm);
	ierr  = DMDAGetInfo(dm, 0, &mx, &my, 0,0,0,0,0,0,0,0,0,0);
	ierr  = DMDAVecGetArray(dm,x,&barray);
	
	double exactSol[mx][my];
	
    char filename[20];
	sprintf(filename, "t=%i.dat",0);
	FILE *fid = fopen (filename, "w");
	fprintf(fid, "title = \"sample mesh\"\n");
    fprintf(fid, "variables = \"x\", \"y\", \"sol\", \"exactSol\", \"error\"\n");
    fprintf(fid, "zone i=%d, j=%d, f=point\n",mx,my);
	for(int j=0;j<my;j++){
       for(int i=0;i<mx;i++){
		   exactSol[i][j] = cos(2*pi*Xc[i][j])*cos(2*pi*Yc[i][j]);
		   double error = exactSol[i][j] - barray[j][i];
           fprintf(fid, "%e   %e   %e   %e   %e\n",Xc[i][j],Yc[i][j],barray[j][i],exactSol[i][j],error);
       }
    }
	fclose(fid);
	
	DMDAVecRestoreArray(dm,x,&barray);
}

/*TEST

   test:
	  args: -pc_type mg -pc_mg_type full -ksp_type fgmres -ksp_monitor_short -pc_mg_levels 3 -mg_coarse_pc_factor_shift_type nonzero

TEST*/
