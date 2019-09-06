FCOMP = ${MPI_PATH}/bin/mpif90
CCOMP = ${MPI_PATH}/bin/mpicc
CPPCOMP = ${MPI_PATH}/bin/mpiCC
LINK = ${MPI_PATH}/bin/mpif90 FFLAGS_OPT = -c -r8 -i4 -O3 -real_size 64# -unroll -align -prefetch -pad -ip
FFLAGS_DEBUG = -c -g -r8 -i4 -check bounds -check format -check output_conversion -warn all -real_size 64
FFLAGS_TEST = -c -r8 -i4 -O2 -real_size 64
CFLAGS_OPT = -c -O3 -D_LARGEFILE64_SOURCE
CFLAGS_DEBUG = -c -g -debug extended -D_LARGEFILE64_SOURCE
CFLAGS_TEST = -c -O2 -D_LARGEFILE64_SOURCE
CFLAGS_HDF5 = -I $(HDF5)/include -std=c++11 
CFLAGS_NCMPI = -I $(NCMPI_PATH)/include
CFLAGS_MPI = -I$(MPI_PATH)/include
LFLAGS_OPT = -r8 -i4 -Vaxlib -lsvml -Ur -o
LFLAGS_DEBUG = -r8 -i4 -Vaxlib -g -o
LFLAGS_TEST = -r8 -i4 -Vaxlib -o
LIB_HDF5 = -L $(HDF5)/lib -lhdf5 -lz
LIB_MPI = -L$(MPI_PATH)/lib -lfmpich -lmpich
LIB_NCMPI = -L$(NCMPI_PATH)/lib -lpnetcdf 

MPICXX=mpicxx
CC=xlc++
CXX=xlc++
CFLAGS=-g -I$(HDF5)/include -Wall

parallel: readFileHDF51

HDF5CLIBS=-L$(HDF5)/lib -lhdf5 

readFileHDF5.o: readFileHDF5.cpp
		$(MPICXX) -g $(CFLAGS_OPT) readFileHDF5.cpp $(CFLAGS_HDF5)	
readFileHDF51: readFileHDF5.o
		$(MPICXX) -g readFileHDF5.o -o readFileHDF5 $(LIB_HDF5)

clean: 
		rm -f *.o readFileHDF5
