#include <iostream>
#include <cmath>
#include "header.h"
using namespace std;

void get_curv_grid(double** &x,double** &y,int Nx,int Ny)
{
	double xmin = 0.0;
	double xmax = 1.0;
	double ymin = 0.0;
	double ymax = 1.0;
	
	double X_cart[Nx][Ny], Y_cart[Nx][Ny];
	double dx = (xmax-xmin)/(Nx-1);
	double dy = (ymax-ymin)/(Ny-1);
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++)
		{
			X_cart[i][j] = xmin+i*dx;	
			Y_cart[i][j] = ymin+j*dy;
		}}
		
	double a = 0.1;
	double b = 2.0;
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++)
		{
			x[i][j] = X_cart[i][j] + a * sin(b*pi*Y_cart[i][j]);
			y[i][j] = Y_cart[i][j] + a * sin(b*pi*X_cart[i][j]);
		}}
		
    char filename[20];
	sprintf(filename, "grid.dat");
	// sprintf(filename, "t=%i.dat",0);
	FILE *fid = fopen (filename, "w");
	fprintf(fid, "title = \"sample mesh\"\n");
    fprintf(fid, "variables = \"x\", \"y\"\n");
    fprintf(fid, "zone i=%d, j=%d, f=point\n",Nx,Ny);
	for(int j=0;j<Ny;j++){
       for(int i=0;i<Nx;i++){
           fprintf(fid, "%e   %e\n",x[i][j],y[i][j]);
       }
    }
	fclose(fid);
}