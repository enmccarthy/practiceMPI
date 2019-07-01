#include <iostream> 
#include "mpi.h" 
#include <fstream> 
#include <sstream>
int main(int argc, char *argv[]) 
{
	int bufsize;
	int rank, nprocs;
	MPI_Init(0,0);
	MPI_File fh;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Offset FILESIZE;
	MPI_File_open(MPI_COMM_WORLD, "test.npy", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_get_size(fh, &FILESIZE);
	std::cout<<FILESIZE;
	std::cout<<" " << sizeof(long long int);
	 
	bufsize = (FILESIZE-128)/nprocs;
	int nints = bufsize/sizeof(long long int);		
	long long int buf[nints];
	//two things I need to do
	//NPY files
	//1. ignore the header and try splitting the file
	//2. assume I have the dimensions of the array
	//	2.a skip the header
	//	3.b get the shape from the header and modify the split of the file
	//      NPZ files
	//      1. if .npz then get each .npy portion of this or look into how it is implemented
	std::cout<<" " << nints << " ";
	std::cout<<rank<< " rank ";
	MPI_File_seek(fh, (rank*bufsize)+128, MPI_SEEK_SET);
	MPI_File_read(fh, buf, nints, MPI_LONG_LONG, &status);
	MPI_File_close(&fh);
	
	std::ofstream output;
	std::string outname = "output";
	std::string out;
	std::stringstream ss;
	ss << rank;
	out = ss.str();
	outname.append(out);
	outname.append(".txt");
	output.open(outname.c_str());

	for(int iter =0; iter < nints; iter++) {
		output << buf[iter] << " ";
		//std::cout<<"here";
	} 
	output.close();
	MPI_Finalize();
	return 0;
		
}
