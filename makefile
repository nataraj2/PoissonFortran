IDIR =-I$(PETSC_DIR)/include
IDIR+=-I$(PETSC_DIR)/arch-darwin-c-debug/include 

FC=mpif90
FFLAGS=-O2

ODIR=.
LDIR=-L$(PETSC_DIR)/arch-darwin-c-debug/lib

LIBS=-lpetsc -llapack -lblas  -lc++ -ldl -lmpifort -lmpi -lpmpi -lgfortran -lquadmath -lm -lc++ -ldl

_DEPS = 

OBJ1 = ModuleVariables.o PetSc_Poisson_Neumann.o PoissonFortran2D.o
OBJ2 = ModuleVariables.o PetSc_Poisson_Periodic.o PoissonFortran2D.o

%.o: %.F90 $(DEPS)
	$(FC) -c -o $@ $< $(IDIR) $(FFLAGS)

poissonfortran2D_Neumann: $(OBJ1)
	$(FC) -o $@ $^ $(IDIR) $(FFLAGS) $(LDIR) $(LIBS)

poissonfortran2D_Periodic: $(OBJ2)
	$(FC) -o $@ $^ $(IDIR) $(FFLAGS) $(LDIR) $(LIBS)


clean:
	rm -f *.o *.mod poissonfortran2D_Neumann poissonfortran2D_Periodic
