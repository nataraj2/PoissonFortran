import numpy as np
from numpy import *
from matplotlib import *
from pylab import *
import math
from scipy.sparse.linalg import gmres


# Define the size of the matrix (n x n)
m = 200

xmin = 0.0; xmax = 1.0;

dx = (xmax-xmin)/m;

matrix_size = m;

# Create an empty matrix A of size m x m
A = np.zeros((matrix_size, matrix_size))

# Set values in the matrix A within a for loop
for i in arange(0,m,1):
	if(i==0):
		A[0,0] = -1.0/dx**2;
		A[0,1] =  1.0/dx**2 
	elif(i==m-1):
		A[m-1,m-1] = -1.0/dx**2
		A[m-1,m-2] =  1.0/dx**2
	else:
		A[i, i] = -2.0/dx**2;
		A[i, i-1] = 1.0/dx**2
		A[i, i+1] = 1.0/dx**2

#print(A)

# Right-hand side vector b
b = np.zeros((matrix_size))
pi = 3.1416

X = zeros(m);
for i in arange(0,m,1):
	X[i] = xmin + (i+0.5)*dx


for i in range(m):
	b[i] = -4.0*pi**2*cos(2.0*pi*X[i])

			
# Solve the system of linear equations Ax = b
sol, info = gmres(A, b)

exact = zeros(m)
for i in arange(0,m,1):
	exact[i] = cos(2*pi*X[i])

for idx in range(m):
	res = sol[idx] - exact[idx]	
	print(res)

pressure = zeros(m);
for i in arange(0,m,1):
	pressure[i] = sol[i] 

plt.figure(1)
plt.plot(X,pressure,'b',label="Numerical")
plt.plot(X,exact,'or',label="Exact")
plt.legend()
plt.show()

