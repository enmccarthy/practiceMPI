
#include <iostream> 
#include "mpi.h" 
 
int main(int argc, char *argv[]) 
{ 
    int bufsize, *buf;
	int rank, nprocs;
	MPI_File fh;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	int FILESIZE;
	bufsize = FILESIZE/nprocs;
	int nints = bufsize/sizeof(int);

	MPI_File_open(MPI_COMM_WORLD, "File", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	//two things I need to do
	//NPY files
	//1. ignore the header and try splitting the file
	//2. assume I have the dimensions of the array
	//	2.a skip the header
	//	3.b get the shape from the header and modify the split of the file
	//NPZ files
	//1. if .npz then get each .npy portion of this or look into how it is implemented
	MPI_File_seek(fh, rank*bufsize, MPI_SEEK_SET);
	MPI_File_read(fh, buf, nints, MPI_INT, &status);
	MPI_File_close(&fh);
}
