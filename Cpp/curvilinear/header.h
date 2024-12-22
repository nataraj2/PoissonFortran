#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;

const double pi = M_PI;
const double e = M_E;

const int space_order = 3;

void get_curv_grid(double**&,double**&,int,int);
	
double** differentiate_fn(double**,int,int,int,int);
void get_face_metric_terms(double**&,double**&,double**&,double**&,double**,double**,double**,double**,int,int);
	
void inner_boundary_derivatives(double***&,double***&,double***,double***);

double** Absolute_of_Matrix(double**,int);
double** AddMatrix(double**,double**,int);
double** SubtractMatrix(double**,double**,int);
double** MultiplyMatrix(double**,double**,int,int,int);

double*** Create3DMatrix(int, int, int);
double** Create2DMatrix(int,int);
void Delete3DMatrix(double***&,int,int,int);
void Delete2DMatrix(double**&,int,int);