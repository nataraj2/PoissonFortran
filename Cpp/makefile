IDIR =-I$(PETSC_DIR)/include
IDIR+=-I$(PETSC_DIR)/arch-darwin-c-debug/include 

CXX=mpicxx
CXXFLAGS=-std=c++14

ODIR=.
LDIR=-L$(PETSC_DIR)/arch-darwin-c-debug/lib

LIBS=-lpetsc -llapack -lblas  -lc++ -ldl -lmpifort -lmpi -lpmpi -lgfortran -lquadmath -lm -lc++ -ldl

_DEPS = 

OBJ1 = Poisson2D_Neumann.o

%.o: %.cpp $(DEPS)
	mpicxx -c -o $@ $< $(IDIR) $(CXXFLAGS)

poisson2d_neumann: $(OBJ1)
	mpicxx -o $@ $^ $(IDIR) $(CXXFLAGS) $(LDIR) $(LIBS)

clean:
	rm -f *.o  poisson2d_neumann
