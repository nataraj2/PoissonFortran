#include <iostream>
#include <cmath>
#include "header.h"
using namespace std;

double** Absolute_of_Matrix(double** theMatrix,int dim)  // ONLY DIAGONAL ELEMENTS
{
	double** AbsMatrix=Create2DMatrix(dim,dim);
	for (int i=0;i<dim;i++)
	{
		AbsMatrix[i][i] = fabs(theMatrix[i][i]);
	}
	return AbsMatrix;
}
		

double** AddMatrix( double** MatrixOne, double** MatrixTwo,int dim)
{
	double** sumMatrix=Create2DMatrix(dim,dim);
	for (int i=0;i<dim;i++)
	{
		for (int j=0;j<dim;j++)
		{
			 sumMatrix[i][j]=MatrixOne[i][j]+MatrixTwo[i][j];
		}
	}
	return sumMatrix;
}

double** SubtractMatrix(double** MatrixOne, double** MatrixTwo,int dim)   // MatrixOne - MatrixTwo
{
	double** subMatrix=Create2DMatrix(dim,dim);
	for (int i=0;i<dim;i++)
	{
		for (int j=0;j<dim;j++)
		{
			 subMatrix[i][j]=MatrixOne[i][j]-MatrixTwo[i][j];
		}
	}
	return subMatrix;
}

double** MultiplyMatrix(double** MatrixOne, double** MatrixTwo,int Nrow_mat1,int Ncol_mat1,int Ncol_mat2)
 {
     double** ProductMatrix=Create2DMatrix(Nrow_mat1,Ncol_mat2);
     for (int i=0;i<Nrow_mat1;i++)
     {
         for (int j=0;j<Ncol_mat2;j++)
         {
			for (int count=0;count<Ncol_mat1;count++)
			{			
				ProductMatrix[i][j] = ProductMatrix[i][j] + MatrixOne[i][count]*MatrixTwo[count][j];
			}
         }
     }
     return ProductMatrix;
 }


double*** Create3DMatrix(int Nk, int Ni, int Nj)
{
    double*** the_array = new double**[Nk];
    for(int k(0); k < Nk; ++k)
    {
        the_array[k] = new double*[Ni];

        for(int i(0); i < Ni; ++i)
        {
            the_array[k][i] = new double[Nj];

            for(int j(0); j < Nj; ++j)
            {
                the_array[k][i][j]= 0.;
            }
        }
    }
    return the_array;
}

double** Create2DMatrix(int Ni, int Nj)
{
    double** the_array = new double* [Ni];
    for(int i(0); i < Ni; ++i)
    {
        the_array[i] = new double[Nj];

        for(int j(0); j < Nj; ++j)
        {
            the_array[i][j] = 0;            
        }
    }
    return the_array;
}


void Delete3DMatrix(double***& the_array, int Nk, int Ni, int Nj)
{
    for (int k = 0; k < Nk; ++k) 
    {
        for (int i = 0; i < Ni; ++i)
        {
            delete [] the_array[k][i];
        }
        delete [] the_array[k];
    }
    delete [] the_array;
}

void Delete2DMatrix(double**& the_array, int Ni, int Nj)
{
    for (int i = 0; i < Ni; ++i)
    {
        delete [] the_array[i];
    }        
	delete [] the_array;
}


