# ----- Make Macros -----

CXX = mpicxx
CXXFLAGS = -std=c++11 -qreport -qlistfmt=html -qsmp=omp
OPTFLAGS = -O3 -qarch=pwr9 -qtune=pwr9 -qsimd=noauto -qstrict 

NVCC = nvcc
NVCCFLAGS = -O3 -std=c++11 -gencode arch=compute_70,code=sm_70 -ccbin=mpicxx -Xcompiler -qsmp=omp

LD_FLAGS = -ccbin=mpicxx -Xcompiler -qsmp=omp

TARGETS = hpcg
OBJECTS = main.o spmv.o jacobi.o gaussseidel.o gradientdescent.o conjugategradient.o setpositions.o

# ----- Make Rules -----

all:	$(TARGETS)

%.o:	%.cpp
	${CXX} ${CXXFLAGS} ${OPTFLAGS} $^ -c -o $@

%.o:    %.cu
	${NVCC} ${NVCCFLAGS} $^ -c -o $@

hpcg:	$(OBJECTS)
	$(NVCC) -o $@ $(OBJECTS) $(LD_FLAGS)

clean:
	rm -f $(TARGETS) *.o *.o.* *.txt *.bin core
