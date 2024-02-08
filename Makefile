GCC=gcc
MPICC=mpicc

CFILES=src/SparseMatrix.c
CFILESLINEAR=src/linear/main.c src/linear/solve_funcs_linear.c
CFILESMPI=src/mpi/main.c src/mpi/solve_funcs_mpi.c
CFILESPTH=src/pthreads/main.c src/pthreads/solve_funcs_pthread.c
CFILESOMP=src/omp/main.c src/omp/solve_funcs_omp.c
HDIR="src/include"
DEBUGADD="_dbg"

OUTNAME_BASE="prog"

OUTNAME_LIN=$(OUTNAME_BASE)"_lin"

# mpich required
# sudo apt-get install mpich
MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)
OUTNAME_MPI=$(OUTNAME_BASE)"_mpi"

PTH_COMPILE_FLAGS = -pthread
PTH_LINK_FLAGS = -lpthread
OUTNAME_PTH=$(OUTNAME_BASE)"_pth"

OMP_COMPILE_FLAGS =
OMP_LINK_FLAGS = -fopenmp
OUTNAME_OMP=$(OUTNAME_BASE)"_omp"

# there is no python compile target so...
# python ./src/main.py 		#for linear

# python -m pip install mpi4py 		#to install mpi4py
# mpiexec -n PROCESS_COUNT python ./src/main_mpi.py 	#for python MPI


all: linear mpi pthread omp

linear:
	$(GCC) -I $(HDIR) $(CFILES) $(CFILESLINEAR) -o $(OUTNAME_LIN)

linear_debug:
	$(GCC) -I $(HDIR) -g $(CFILES) $(CFILESLINEAR) -o $(OUTNAME_LIN)$(DEBUGADD)

#launch with "mpiexec -n count_of_threads ./OUTNAME"
mpi:
	$(GCC) -I $(HDIR) $(MPI_COMPILE_FLAGS) $(CFILES) $(CFILESMPI) $(MPI_LINK_FLAGS) -o $(OUTNAME_MPI)

mpi_debug:
	$(GCC) -I $(HDIR) -g -D DEBUG_MPI $(MPI_COMPILE_FLAGS) $(CFILES) $(CFILESMPI) $(MPI_LINK_FLAGS) -o $(OUTNAME_MPI)$(DEBUGADD)

pthread:
	$(GCC) -I $(HDIR) $(PTH_COMPILE_FLAGS) $(CFILES) $(CFILESPTH) $(PTH_LINK_FLAGS) -o $(OUTNAME_PTH)

pthread_debug:
	$(GCC) -I $(HDIR) -g -D DEBUG_PTH $(PTH_COMPILE_FLAGS) $(CFILES) $(CFILESPTH) $(PTH_LINK_FLAGS) -o $(OUTNAME_PTH)$(DEBUGADD)

omp:
	$(GCC) -I $(HDIR) $(OMP_COMPILE_FLAGS) $(CFILES) $(CFILESOMP) $(OMP_LINK_FLAGS) -o $(OUTNAME_OMP)

omp_debug:
	$(GCC) -I $(HDIR) -g -D DEBUG_OMP $(OMP_COMPILE_FLAGS) $(CFILES) $(CFILESOMP) $(OMP_LINK_FLAGS) -o $(OUTNAME_OMP)$(DEBUGADD)

clean:
	rm ./$(OUTNAME_BASE)*