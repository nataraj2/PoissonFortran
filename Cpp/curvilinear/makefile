IDIR =-I$(PETSC_DIR)/include
IDIR+=-I$(PETSC_DIR)/arch-darwin-c-debug/include 

CXX=mpicxx
CFLAGS = -mmacosx-version-min=14.5 
LDFLAGS = -mmacosx-version-min=14.5
CXXFLAGS=-std=c++14

ODIR=.
LDIR=-L$(PETSC_DIR)/arch-darwin-c-debug/lib

LIBS=-lpetsc -llapack -lblas  -lc++ -ldl -lmpifort -lmpi -lpmpi -lgfortran -lquadmath -lm -lc++ -ldl

_DEPS = header.h

OBJ1 = generate_grid.o differentiate_fn.o matrix_ops.o Poisson2D_Neumann.o 

%.o: %.cpp $(DEPS)
	mpicxx -c -o $@ $< $(IDIR) $(CXXFLAGS) $(CFLAGS) $(LDFLAGS)

output.out: $(OBJ1)
	mpicxx -o $@ $^ $(IDIR) $(CXXFLAGS) $(CFLAGS) $(LDFLAGS) $(LDIR) $(LIBS)

clean:
	rm -f *.o output.out
