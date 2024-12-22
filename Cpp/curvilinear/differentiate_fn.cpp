#include <iostream>
#include <cmath>
#include "header.h"
using namespace std;

double** differentiate_fn(double** f,int Nx,int Ny,int dir,int periodicity)
{
	double** dfdxi = Create2DMatrix(Nx,Ny);
	int i,j;
	double dxi = 1.0;

	//////////////////////////////////////////////////////////////////////////
	///////////////////// SPACE ORDER = 2 ////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////

	if (space_order==2)
	{
		if (dir==1)
		{
			// INTERIOR DERIVATIVE APPROXIMATION
			for(i=1;i<Nx-1;i++){	
				for(j=0;j<Ny;j++)
				{
					dfdxi[i][j] = (f[i+1][j]-f[i-1][j])/(2.0*dxi);
				}}

			// LEFT BOUNDARY
			i = 0;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = (f[i+1][j]-f[i][j])/(dxi);
				}
				else if (periodicity==1){
					dfdxi[i][j] = (f[i+1][j]-f[Nx-2][j])/(2.0*dxi);
				}
			}

			// RIGHT BOUNDARY
			i = Nx-1;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = (f[i][j]-f[i-1][j])/(dxi);
				}
				else if (periodicity==1){
					dfdxi[i][j] = (f[1][j]-f[i-1][j])/(2.0*dxi);
				}
			}
		}

		else if (dir==2)
		{
			// INTERIOR DERIVATIVE APPROXIMATION
			for(i=0;i<Nx;i++){	
				for(j=1;j<Ny-1;j++)
				{
					dfdxi[i][j] = (f[i][j+1]-f[i][j-1])/(2*dxi);
				}}

			// BOTTOM BOUNDARY
			j = 0;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = (f[i][j+1]-f[i][j])/(dxi);
			}

			// TOP BOUNDARY
			j = Ny-1;
			for(int i=0;i<Nx;i++)
			{
				dfdxi[i][j] = (f[i][j]-f[i][j-1])/(dxi);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	///////////////////// SPACE ORDER = 3 ////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////

	else if (space_order==3)
	{
		if (dir==1)
		{
			// INTERIOR DERIVATIVE APPROXIMATION
			for(i=4;i<Nx-4;i++){	
				for(j=0;j<Ny;j++)
				{
					dfdxi[i][j] = -1.0/12*f[i+2][j] + 2.0/3*f[i+1][j] - 2.0/3*f[i-1][j] + 1.0/12*f[i-2][j];
				}}

			// LEFT BOUNDARY
			i = 0;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -24.0/17*f[i][j] + 59.0/34*f[i+1][j] - 4.0/17*f[i+2][j] -3.0/34*f[i+3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = -1.0/12*f[i+2][j] + 2.0/3*f[i+1][j] - 2.0/3*f[Nx-2][j] + 1.0/12*f[Nx-3][j];
				}
			}

			i = 1;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -0.5*f[i-1][j] + 0.5*f[i+1][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = -1.0/12*f[i+2][j] + 2.0/3*f[i+1][j] - 2.0/3*f[i-1][j] + 1.0/12*f[Nx-2][j];
				}
			}

			i = 2;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 4.0/43*f[i-2][j] - 59.0/86*f[i-1][j] + 59.0/86*f[i+1][j] - 4.0/43*f[i+2][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = -1.0/12*f[i+2][j] + 2.0/3*f[i+1][j] - 2.0/3*f[i-1][j] + 1.0/12*f[i-2][j];
				}
			}

			i = 3;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 3.0/98*f[i-3][j] - 59.0/98*f[i-1][j] + 32.0/49*f[i+1][j] - 4.0/49*f[i+2][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = -1.0/12*f[i+2][j] + 2.0/3*f[i+1][j] - 2.0/3*f[i-1][j] + 1.0/12*f[i-2][j];
				}
			}

			// RIGHT BOUNDARY
			i = Nx-1;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 24.0/17*f[i][j] - 59.0/34*f[i-1][j] + 4.0/17*f[i-2][j] + 3.0/34*f[i-3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = -1.0/12*f[2][j] + 2.0/3*f[1][j] - 2.0/3*f[Nx-2][j] + 1.0/12*f[Nx-3][j];
				}
			}

			i = Nx-2;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -0.5*f[i-1][j] + 0.5*f[i+1][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = -1.0/12*f[1][j] + 2.0/3*f[i+1][j] - 2.0/3*f[i-1][j] + 1.0/12*f[i-2][j];
				}
			}

			i = Nx-3;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 4.0/43*f[i-2][j] - 59.0/86*f[i-1][j] + 59.0/86*f[i+1][j] - 4.0/43*f[i+2][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = -1.0/12*f[i+2][j] + 2.0/3*f[i+1][j] - 2.0/3*f[i-1][j] + 1.0/12*f[i-2][j];
				}
			}

			i = Nx-4;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -3.0/98*f[i+3][j] + 59.0/98*f[i+1][j] - 32.0/49*f[i-1][j] + 4.0/49*f[i-2][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = -1.0/12*f[i+2][j] + 2.0/3*f[i+1][j] - 2.0/3*f[i-1][j] + 1.0/12*f[i-2][j];
				}
			}
		}
	
		else if (dir==2)
		{
			// INTERIOR DERIVATIVE APPROXIMATION
			for(i=0;i<Nx;i++){	
				for(j=4;j<Ny-4;j++)
				{
					dfdxi[i][j] = -1.0/12*f[i][j+2] + 2.0/3*f[i][j+1] - 2.0/3*f[i][j-1] + 1.0/12*f[i][j-2];
				}}

			// BOTTOM BOUNDARY
			j = 0;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -24.0/17*f[i][j] + 59.0/34*f[i][j+1] - 4.0/17*f[i][j+2] -3.0/34*f[i][j+3];
			}

			j = 1;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -0.5*f[i][j-1] + 0.5*f[i][j+1];			
			}

			j = 2;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = 4.0/43*f[i][j-2] - 59.0/86*f[i][j-1] + 59.0/86*f[i][j+1] - 4.0/43*f[i][j+2];
			}

			j = 3;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = 3.0/98*f[i][j-3] - 59.0/98*f[i][j-1] + 32.0/49*f[i][j+1] - 4.0/49*f[i][j+2];			
			}

			// TOP BOUNDARY
			j = Ny-1;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = 24.0/17*f[i][j] - 59.0/34*f[i][j-1] + 4.0/17*f[i][j-2] + 3.0/34*f[i][j-3];
			}

			j = Ny-2;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -0.5*f[i][j-1] + 0.5*f[i][j+1];
			}

			j = Ny-3;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = 4.0/43*f[i][j-2] - 59.0/86*f[i][j-1] + 59.0/86*f[i][j+1] - 4.0/43*f[i][j+2];
			}

			j = Ny-4;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -3.0/98*f[i][j+3] + 59.0/98*f[i][j+1] - 32.0/49*f[i][j-1] + 4.0/49*f[i][j-2];
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	///////////////////// SPACE ORDER = 4 ////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////

	else if (space_order==4)
	{
		double x1 = 0.660731095679012;
		double D1[6][9];
		D1[0][0] = -21600.0/13649;
		D1[0][1] = 8.0*(16200*x1-953)/40947;
		D1[0][2] = (-1036800.0*x1+715489)/81894;
		D1[0][3] = 3.0*(86400*x1-62639)/13649;
		D1[0][4] = 5.0*(-207360*x1+147127)/81894; 
		D1[0][5] = (129600.0*x1-89387)/40947;
	
		D1[1][0] = 8.0*(-16200*x1+953)/180195;
		D1[1][2] = (86400.0*x1-57139)/12013;
		D1[1][3] = (-1036800.0*x1+745733)/72078;
		D1[1][4] = 5.0*(25920*x1-18343)/12013;
		D1[1][5] = (-345600.0*x1+240569)/120130;
	
		D1[2][0] = (1036800.0*x1-715489)/162660;
		D1[2][1] = (-86400.0*x1+57139)/5422;
		D1[2][3] = (259200.0*x1-176839)/8133;
		D1[2][4] = (-345600.0*x1+242111)/10844;
		D1[2][5] = (259200.0*x1-182261)/27110;
	
		D1[3][0] = 3.0*(-86400.0*x1+62639)/53590;
		D1[3][1] = (1036800.0*x1-745733)/64308;
		D1[3][2] = (-259200.0*x1+176839)/16077;
		D1[3][4] = (259200.0*x1-165041)/32154;
		D1[3][5] = (-1036800.0*x1+710473)/321540;
		D1[3][6] = 72.0/5359;
	
		D1[4][0] = (207360.0*x1-147127)/47262;
		D1[4][1] = 5.0*(-25920.0*x1+18343)/7877;
		D1[4][2] = (345600.0*x1-242111)/15754;
		D1[4][3] = (-259200.0*x1+165041)/23631;
		D1[4][5] = 8640.0*x1/7877;
		D1[4][6] = -1296.0/7877;
		D1[4][7] = 144.0/7877;
	
		D1[5][0] = (-129600.0*x1+89387)/131403;
		D1[5][1] = (345600.0*x1-240569)/87602;
		D1[5][2] = (-259200.0*x1+182261)/43801;
		D1[5][3] = (1036800.0*x1-710473)/262806;
		D1[5][4] = -43200.0*x1/43801;
		D1[5][6] = 32400.0/43801;
		D1[5][7] = -6480.0/43801;
		D1[5][8] = 720.0/43801;

		if (dir==1)
		{
			// INTERIOR DERIVATIVE APPROXIMATION
			for(i=4;i<Nx-4;i++){	
				for(j=0;j<Ny;j++)
				{
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[i-3][j];
				}}

			// LEFT BOUNDARY
			i = 0;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -21600.0/13649*f[i][j] + 104009.0/54596*f[i+1][j] + 30443.0/81894*f[i+2][j] - 33311.0/27298*f[i+3][j] + 16863.0/27298*f[i+4][j] - 15025.0/163788*f[i+5][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[Nx-2][j] + 3.0/20*f[Nx-3][j] - 1.0/60*f[Nx-4][j];
				}
			}

			i = 1;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -104009.0/240260*f[i-1][j] - 311.0/72078*f[i+1][j] + 20229.0/24026*f[i+2][j] - 24337.0/48052*f[i+3][j] + 36661.0/360390*f[i+4][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[Nx-2][j] - 1.0/60*f[Nx-3][j];
				}
			}

			i = 2;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -30443.0/162660*f[i-2][j] + 311.0/32532*f[i-1][j] - 11155.0/16266*f[i+1][j] + 41287.0/32532*f[i+2][j] - 21999.0/54220*f[i+3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[Nx-2][j];
				}
			}

			i = 3;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 33311.0/107180*f[i-3][j] - 20229.0/21436*f[i-2][j] + 485.0/1398*f[i-1][j] + 4147.0/21436*f[i+1][j] + 25427.0/321540*f[i+2][j] + 72.0/5359*f[i+3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[i-3][j];
				}
			}

			i = 4;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -16863.0/78770*f[i-4][j] + 24337.0/31508*f[i-3][j] - 41287.0/47262*f[i-2][j] - 4147.0/15754*f[i-1][j] + 342523.0/472620*f[i+1][j] - 1296.0/7877*f[i+2][j] + 144.0/7877*f[i+3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[i-3][j];
				}
			}

			i = 5;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 15025.0/525612*f[i-5][j] - 36661.0/262806*f[i-4][j] + 21999.0/87602*f[i-3][j] - 25427.0/262806*f[i-2][j] - 342523.0/525612*f[i-1][j] + 32400.0/43801*f[i+1][j] - 6480.0/43801*f[i+2][j] + 720.0/43801*f[i+3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[i-3][j];
				}
			}

			// RIGHT BOUNDARY
			i = Nx-1;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 21600.0/13649*f[i][j] - 104009.0/54596*f[i-1][j] - 30443.0/81894*f[i-2][j] + 33311.0/27298*f[i-3][j] - 16863.0/27298*f[i-4][j] + 15025.0/163788*f[i-5][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[3][j] - 3.0/20*f[2][j] + 3.0/4*f[1][j] - 3.0/4*f[Nx-2][j] + 3.0/20*f[Nx-3][j] - 1.0/60*f[Nx-4][j];
				}
			}

			i = Nx-2;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 104009.0/240260*f[i+1][j] + 311.0/72078*f[i-1][j] - 20229.0/24026*f[i-2][j] + 24337.0/48052*f[i-3][j] - 36661.0/360390*f[i-4][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[2][j] - 3.0/20*f[1][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[i-3][j];
				}
			}

			i = Nx-3;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 30443.0/162660*f[i+2][j] - 311.0/32532*f[i+1][j] + 11155.0/16266*f[i-1][j] - 41287.0/32532*f[i-2][j] + 21999.0/54220*f[i-3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[1][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[i-3][j];
				}
			}

			i = Nx-4;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -33311.0/107180*f[i+3][j] + 20229.0/21436*f[i+2][j] - 485.0/1398*f[i+1][j] - 4147.0/21436*f[i-1][j] - 25427.0/321540*f[i-2][j] - 72.0/5359*f[i-3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[i-3][j];
				}
			}

			i = Nx-5;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = 16863.0/78770*f[i+4][j] - 24337.0/31508*f[i+3][j] + 41287.0/47262*f[i+2][j] + 4147.0/15754*f[i+1][j] - 342523.0/472620*f[i-1][j] + 1296.0/7877*f[i-2][j] - 144.0/7877*f[i-3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[i-3][j];
				}
			}

			i = Nx-6;
			for(int j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = -15025.0/525612*f[i+5][j] + 36661.0/262806*f[i+4][j] - 21999.0/87602*f[i+3][j] + 25427.0/262806*f[i+2][j] + 342523.0/525612*f[i+1][j] - 32400.0/43801*f[i-1][j] + 6480.0/43801*f[i-2][j] - 720.0/43801*f[i-3][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[i-2][j] - 1.0/60*f[i-3][j];
				}
			}
		}
	
		else if (dir==2)
		{
			// INTERIOR DERIVATIVE APPROXIMATION
			for(i=0;i<Nx;i++){	
				for(j=4;j<Ny-4;j++)
				{
					dfdxi[i][j] = 1.0/60*f[i][j+3] - 3.0/20*f[i][j+2] + 3.0/4*f[i][j+1] - 3.0/4*f[i][j-1] + 3.0/20*f[i][j-2] - 1.0/60*f[i][j-3];
				}}

			// BOTTOM BOUNDARY
			j = 0;
			for(i=0;i<Nx;i++)
			{
					dfdxi[i][j] = -21600.0/13649*f[i][j] + 104009.0/54596*f[i][j+1] + 30443.0/81894*f[i][j+2] - 33311.0/27298*f[i][j+3] + 16863.0/27298*f[i][j+4] - 15025.0/163788*f[i][j+5];
			}

			j = 1;
			for(i=0;i<Nx;i++)
			{
					dfdxi[i][j] = -104009.0/240260*f[i][j-1] - 311.0/72078*f[i][j+1] + 20229.0/24026*f[i][j+2] - 24337.0/48052*f[i][j+3] + 36661.0/360390*f[i][j+4];
			}

			j = 2;
			for(i=0;i<Nx;i++)
			{
					dfdxi[i][j] = -30443.0/162660*f[i][j-2] + 311.0/32532*f[i][j-1] - 11155.0/16266*f[i][j+1] + 41287.0/32532*f[i][j+2] - 21999.0/54220*f[i][j+3];
			}

			j = 3;
			for(i=0;i<Nx;i++)
			{
					dfdxi[i][j] = 33311.0/107180*f[i][j-3] - 20229.0/21436*f[i][j-2] + 485.0/1398*f[i][j-1] + 4147.0/21436*f[i][j+1] + 25427.0/321540*f[i][j+2] + 72.0/5359*f[i][j+3];
			}

			j = 4;
			for(i=0;i<Nx;i++)
			{
					dfdxi[i][j] = -16863.0/78770*f[i][j-4] + 24337.0/31508*f[i][j-3] - 41287.0/47262*f[i][j-2] - 4147.0/15754*f[i][j-1] + 342523.0/472620*f[i][j+1] - 1296.0/7877*f[i][j+2] + 144.0/7877*f[i][j+3];
			}

			j = 5;
			for(i=0;i<Nx;i++)
			{
					dfdxi[i][j] = 15025.0/525612*f[i][j-5] - 36661.0/262806*f[i][j-4] + 21999.0/87602*f[i][j-3] - 25427.0/262806*f[i][j-2] - 342523.0/525612*f[i][j-1] + 32400.0/43801*f[i][j+1] - 6480.0/43801*f[i][j+2] + 720.0/43801*f[i][j+3];
			}


			// TOP BOUNDARY
			j = Ny-1;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = 21600.0/13649*f[i][j] - 104009.0/54596*f[i][j-1] - 30443.0/81894*f[i][j-2] + 33311.0/27298*f[i][j-3] - 16863.0/27298*f[i][j-4] + 15025.0/163788*f[i][j-5];
			}

			j = Ny-2;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = 104009.0/240260*f[i][j+1] + 311.0/72078*f[i][j-1] - 20229.0/24026*f[i][j-2] + 24337.0/48052*f[i][j-3] - 36661.0/360390*f[i][j-4];
			}

			j = Ny-3;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = 30443.0/162660*f[i][j+2] - 311.0/32532*f[i][j+1] + 11155.0/16266*f[i][j-1] - 41287.0/32532*f[i][j-2] + 21999.0/54220*f[i][j-3];
			}

			j = Ny-4;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -33311.0/107180*f[i][j+3] + 20229.0/21436*f[i][j+2] - 485.0/1398*f[i][j+1] - 4147.0/21436*f[i][j-1] - 25427.0/321540*f[i][j-2] - 72.0/5359*f[i][j-3];
			}

			j = Ny-5;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = 16863.0/78770*f[i][j+4] - 24337.0/31508*f[i][j+3] + 41287.0/47262*f[i][j+2] + 4147.0/15754*f[i][j+1] - 342523.0/472620*f[i][j-1] + 1296.0/7877*f[i][j-2] - 144.0/7877*f[i][j-3];
			}

			j = Ny-6;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -15025.0/525612*f[i][j+5] + 36661.0/262806*f[i][j+4] - 21999.0/87602*f[i][j+3] + 25427.0/262806*f[i][j+2] + 342523.0/525612*f[i][j+1] - 32400.0/43801*f[i][j-1] + 6480.0/43801*f[i][j-2] - 720.0/43801*f[i][j-3];
			}

		}
/*			// LEFT BOUNDARY
			i = 0;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = D1[0][i]*f[i][j] + D1[0][i+1]*f[i+1][j] + D1[0][i+2]*f[i+2][j] + D1[0][i+3]*f[i+3][j] + D1[0][i+4]*f[i+4][j] + D1[0][i+5]*f[i+5][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[Nx-2][j] + 3.0/20*f[Nx-3][j] - 1.0/60*f[Nx-4][j];
				}
			}

			i = 1;
			for(j=0;j<Ny;j++)
			{
				if (periodicity==0){
					dfdxi[i][j] = D1[1][i-1]*f[i-1][j] + D1[1][i+1]*f[i+1][j] + D1[1][i+2]*f[i+2][j] + D1[1][i+3]*f[i+3][j] + D1[1][i+4]*f[i+4][j];
				}
				else if (periodicity==1){
					dfdxi[i][j] = 1.0/60*f[i+3][j] - 3.0/20*f[i+2][j] + 3.0/4*f[i+1][j] - 3.0/4*f[i-1][j] + 3.0/20*f[Nx-2][j] - 1.0/60*f[Nx-3][j];
				}
			}

			i = 2;
			for(j=0;j<Ny;j++)
			{
				dfdxi[i][j] = D1[2][i-2]*f[i-2][j] + D1[2][i-1]*f[i-1][j] + D1[2][i+1]*f[i+1][j] + D1[2][i+2]*f[i+2][j] + D1[2][i+3]*f[i+3][j];
			}

			i = 3;
			for(j=0;j<Ny;j++)
			{
				dfdxi[i][j] = D1[3][i-3]*f[i-3][j] + D1[3][i-2]*f[i-2][j] + D1[3][i-1]*f[i-1][j] + D1[3][i+1]*f[i+1][j] + D1[3][i+2]*f[i+2][j] + D1[3][i+3]*f[i+3][j];
			}

			i = 4;
			for(j=0;j<Ny;j++)
			{
				dfdxi[i][j] = D1[4][i-4]*f[i-4][j] + D1[4][i-3]*f[i-3][j] + D1[4][i-2]*f[i-2][j] + D1[4][i-1]*f[i-1][j] + D1[4][i+1]*f[i+1][j] + D1[4][i+2]*f[i+2][j] + D1[4][i+3]*f[i+3][j];
			}

			i = 5;
			for(j=0;j<Ny;j++)
			{
				dfdxi[i][j] = D1[5][i-5]*f[i-5][j] + D1[5][i-4]*f[i-4][j] + D1[5][i-3]*f[i-3][j] + D1[5][i-2]*f[i-2][j] + D1[5][i-1]*f[i-1][j] + D1[5][i+1]*f[i+1][j] + D1[5][i+2]*f[i+2][j] + D1[5][i+3]*f[i+3][j];
			}

			// RIGHT BOUNDARY
			i = Nx-1;
			for(int j=0;j<Ny;j++)
			{
				dfdxi[i][j] = -(D1[0][0]*f[i][j] + D1[0][1]*f[i-1][j] + D1[0][2]*f[i-2][j] + D1[0][3]*f[i-3][j] + D1[0][4]*f[i-4][j] + D1[0][5]*f[i-5][j]);
			}

			i = Nx-2;
			for(int j=0;j<Ny;j++)
			{
				dfdxi[i][j] = -(D1[1][0]*f[i+1][j] + D1[1][2]*f[i-1][j] + D1[1][3]*f[i-2][j] + D1[1][4]*f[i-3][j] + D1[1][5]*f[i-4][j]);
			}

			i = Nx-3;
			for(int j=0;j<Ny;j++)
			{
				dfdxi[i][j] = -(D1[2][0]*f[i+2][j] + D1[2][1]*f[i+1][j] + D1[2][3]*f[i-1][j] + D1[2][4]*f[i-2][j] + D1[2][5]*f[i-3][j]);
			}

			i = Nx-4;
			for(int j=0;j<Ny;j++)
			{
				dfdxi[i][j] = -(D1[3][0]*f[i+3][j] + D1[3][1]*f[i+2][j] + D1[3][2]*f[i+1][j] + D1[3][4]*f[i-1][j] + D1[3][5]*f[i-2][j] + D1[3][6]*f[i-3][j]);
			}

			i = Nx-5;
			for(int j=0;j<Ny;j++)
			{
				dfdxi[i][j] = -(D1[4][0]*f[i+4][j] + D1[4][1]*f[i+3][j] + D1[4][2]*f[i+2][j] + D1[4][3]*f[i+1][j] + D1[4][5]*f[i-1][j] + D1[4][6]*f[i-2][j] + D1[4][7]*f[i-3][j]);
			}

			i = Nx-6;
			for(int j=0;j<Ny;j++)
			{
				dfdxi[i][j] = -(D1[5][0]*f[i+5][j] + D1[5][1]*f[i+4][j] + D1[5][2]*f[i+3][j] + D1[5][3]*f[i+2][j] + D1[5][4]*f[i+1][j] + D1[5][6]*f[i-1][j] + D1[5][7]*f[i-2][j] + D1[5][8]*f[i-3][j]);
			}
		}
	
		else if (dir==2)
		{
			// INTERIOR DERIVATIVE APPROXIMATION
			for(i=0;i<Nx;i++){	
				for(j=4;j<Ny-4;j++)
				{
					dfdxi[i][j] = 1.0/60*f[i][j+3] - 3.0/20*f[i][j+2] + 3.0/4*f[i][j+1] - 3.0/4*f[i][j-1] + 3.0/20*f[i][j-2] - 1.0/60*f[i][j-3];
				}}

			// BOTTOM BOUNDARY
			j = 0;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = D1[0][j]*f[i][j] + D1[0][j+1]*f[i][j+1] + D1[0][j+2]*f[i][j+2] + D1[0][j+3]*f[i][j+3] + D1[0][j+4]*f[i][j+4] + D1[0][j+5]*f[i][j+5];
			}

			j = 1;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = D1[1][j-1]*f[i][j-1] + D1[1][j+1]*f[i][j+1] + D1[1][j+2]*f[i][j+2] + D1[1][j+3]*f[i][j+3] + D1[1][j+4]*f[i][j+4];
			}

			j = 2;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = D1[2][j-2]*f[i][j-2] + D1[2][j-1]*f[i][j-1] + D1[2][j+1]*f[i][j+1] + D1[2][j+2]*f[i][j+2] + D1[2][j+3]*f[i][j+3];
			}

			j = 3;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = D1[3][j-3]*f[i][j-3] + D1[3][j-2]*f[i][j-2] + D1[3][j-1]*f[i][j-1] + D1[3][j+1]*f[i][j+1] + D1[3][j+2]*f[i][j+2] + D1[3][j+3]*f[i][j+3];
			}

			j = 4;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = D1[4][j-4]*f[i][j-4] + D1[4][j-3]*f[i][j-3] + D1[4][j-2]*f[i][j-2] + D1[4][j-1]*f[i][j-1] + D1[4][j+1]*f[i][j+1] + D1[4][j+2]*f[i][j+2] + D1[4][j+3]*f[i][j+3];
			}

			j = 5;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = D1[5][j-5]*f[i][j-5] + D1[5][j-4]*f[i][j-4] + D1[5][j-3]*f[i][j-3] + D1[5][j-2]*f[i][j-2] + D1[5][j-1]*f[i][j-1] + D1[5][j+1]*f[i][j+1] + D1[5][j+2]*f[i][j+2] + D1[5][j+3]*f[i][j+3];
			}


			// TOP BOUNDARY
			j = Ny-1;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -(D1[0][0]*f[i][j] + D1[0][1]*f[i][j-1] + D1[0][2]*f[i][j-2] + D1[0][3]*f[i][j-3] + D1[0][4]*f[i][j-4] + D1[0][5]*f[i][j-5]);
			}

			j = Ny-2;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -(D1[1][0]*f[i][j+1] + D1[1][2]*f[i][j-1] + D1[1][3]*f[i][j-2] + D1[1][4]*f[i][j-3] + D1[1][5]*f[i][j-4]);
			}

			j = Ny-3;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -(D1[2][0]*f[i][j+2] + D1[2][1]*f[i][j+1] + D1[2][3]*f[i][j-1] + D1[2][4]*f[i][j-2] + D1[2][5]*f[i][j-3]);
			}

			j = Ny-4;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -(D1[3][0]*f[i][j+3] + D1[3][1]*f[i][j+2] + D1[3][2]*f[i][j+1] + D1[3][4]*f[i][j-1] + D1[3][5]*f[i][j-2] + D1[3][6]*f[i][j-3]);
			}

			j = Ny-5;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -(D1[4][0]*f[i][j+4] + D1[4][1]*f[i][j+3] + D1[4][2]*f[i][j+2] + D1[4][3]*f[i][j+1] + D1[4][5]*f[i][j-1] + D1[4][6]*f[i][j-2] + D1[4][7]*f[i][j-3]);
			}

			j = Ny-6;
			for(i=0;i<Nx;i++)
			{
				dfdxi[i][j] = -(D1[5][0]*f[i][j+5] + D1[5][1]*f[i][j+4] + D1[5][2]*f[i][j+3] + D1[5][3]*f[i][j+2] + D1[5][4]*f[i][j+1] + D1[5][6]*f[i][j-1] + D1[5][7]*f[i][j-2] + D1[5][8]*f[i][j-3]);
			}

		}*/

	}

	return dfdxi;
}

void get_face_metric_terms(double** &xi_xf,double** &xi_yf,double** &eta_xf,double** &eta_yf,double** Xfx,double** Yfx,double** Xfy,double** Yfy,
							int nx,int ny)
{
	// Calculate xi_xf & xi_yf at the X-FACE-POINTS
    double** x_xi = differentiate_fn(Xfx,nx,ny-1,1,0);
    double** x_eta = differentiate_fn(Xfx,nx,ny-1,2,0);
    double** y_xi = differentiate_fn(Yfx,nx,ny-1,1,0);
    double** y_eta = differentiate_fn(Yfx,nx,ny-1,2,0);
	
    for(int i=0;i<nx;i++){	
  		for(int j=0;j<ny-1;j++)
  		{
			double inv_Jac = x_xi[i][j]*y_eta[i][j] - x_eta[i][j]*y_xi[i][j];
			double Jac = 1.0/inv_Jac;
			xi_xf[i][j] = Jac * y_eta[i][j];
			xi_yf[i][j] = -Jac * x_eta[i][j];
		}}
	
	Delete2DMatrix(x_xi,nx,ny-1);
	Delete2DMatrix(x_eta,nx,ny-1);
	Delete2DMatrix(y_xi,nx,ny-1);
	Delete2DMatrix(y_eta,nx,ny-1);
	
	// Calculate eta_xf & eta_yf at the Y-FACE-POINTS
    x_xi = differentiate_fn(Xfy,nx-1,ny,1,0);
    x_eta = differentiate_fn(Xfy,nx-1,ny,2,0);
    y_xi = differentiate_fn(Yfy,nx-1,ny,1,0);
    y_eta = differentiate_fn(Yfy,nx-1,ny,2,0);
	
    for(int i=0;i<nx-1;i++){	
  		for(int j=0;j<ny;j++)
  		{
			double inv_Jac = x_xi[i][j]*y_eta[i][j] - x_eta[i][j]*y_xi[i][j];
			double Jac = 1.0/inv_Jac;
			eta_xf[i][j] = -Jac * y_xi[i][j];
			eta_yf[i][j] = Jac * x_xi[i][j];
		}}
		
	Delete2DMatrix(x_xi,nx-1,ny);
	Delete2DMatrix(x_eta,nx-1,ny);
	Delete2DMatrix(y_xi,nx-1,ny);
	Delete2DMatrix(y_eta,nx-1,ny);	
}



