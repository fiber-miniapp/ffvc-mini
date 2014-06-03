#CC = gcc
CC = mpicc
CFLAGS = -Wall -g 

#FC = gfortran
FC = mpif90
FFLAGS = -Wall -g

#CXX = g++
CXX = mpic++
CXXFLAGS = -Wall -g

C_OBJS = c_test.o
F_OBJS = f_test.o
CXX_OBJS = cxx_test.o

all: c_test

export

c_lib:
	$(MAKE) -C ../src c MAPROF_C="CC CFLAGS" OUTPUT=c.yaml

c_mpi_lib:
	$(MAKE) -C ../src c_mpi MAPROF_C="CC CFLAGS" OUTPUT=c.yaml

f_lib:
	$(MAKE) -C ../src f MAPROF_F="FC FFLAGS" OUTPUT=f.yaml

f_mpi_lib:
	$(MAKE) -C ../src f_mpi MAPROF_F="FC FFLAGS" OUTPUT=f.yaml

cxx_lib:
	$(MAKE) -C ../src cxx MAPROF_CXX="CXX CXXFLAGS" OUTPUT=cxx.yaml

cxx_mpi_lib:
	$(MAKE) -C ../src cxx_mpi MAPROF_CXX="CXX CXXFLAGS" OUTPUT=cxx.yaml


CFLAGS += -I../src
FFLAGS += -I../src
CXXFLAGS += -I../src

LDFLAGS += -L../src

c_test c_mpi_test cxx_test cxx_mpi_test: LIBS += -lmaprof
f_test f_mpi_test: LIBS += -lmaprof_f

c_mpi_test: CFLAGS += -DUSE_MPI
f_mpi_test: FFLAGS += -DUSE_MPI
cxx_mpi_test: CXXFLAGS += -DUSE_MPI

c_test c_mpi_test cxx_test cxx_mpi_test: LIBS += -lm

c_test: c_lib $(C_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(C_OBJS) -o $@ $(LIBS)

c_mpi_test: c_mpi_lib $(C_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(C_OBJS) -o $@ $(LIBS)

f_test: f_lib $(F_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) $(F_OBJS) -o $@ $(LIBS)

f_mpi_test: f_mpi_lib $(F_OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) $(F_OBJS) -o $@ $(LIBS)

cxx_test: cxx_lib $(CXX_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CXX_OBJS) -o $@ $(LIBS)

cxx_mpi_test: cxx_mpi_lib $(CXX_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CXX_OBJS) -o $@ $(LIBS)

.SUFFIXES:
.SUFFIXES: .F90 .c .cpp .o .mod

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

.F90.o:
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod *.i90
	$(MAKE) -C ../src clean

distclean: clean
	rm -f c_test c_mpi_test f_test f_mpi_test cxx_test cxx_mpi_test \
	      c.yaml f.yaml cxx.yaml


