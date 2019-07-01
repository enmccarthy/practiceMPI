#include <iostream> 
#include "mpi.h" 
#include <fstream> 
#include <sstream>
#include <vector>
#include <cstring>
#include <string>
//read file using MPI
// filename of input file
// vector of ints which represent the shape of the datastructure (ex: (1,2,3), vector<int> vect{1,2,3})
// 	did this to accomodate different shapes
// dtype should be an int representing number of bytes of possibly MPI_data_type (**come back to this**)
// debug, should it write out the values to check for acc or whatever debug is implemented
int main(int argc, char *argv[]) 
{
	// from the command line pass in a path to the file
	//program name, file name, dtype, shape 
	MPI_Init(&argc,&argv);
	bool debug = true;
	std::string filename = "";
	std::string dtype = "";
	std::vector<int> shape;
	std::cout<<"argc "<< argc ;
	for (int i =0; i <argc; i++) { 
		std::cout<<"in for loop "<<argv[i]<<std::endl;
		std::cout<<(strcmp(argv[i], "-s")==0)<<" string compare " << argv[i]<<std::endl;
		if (strcmp(argv[i], "-s")==0) {
			std::cout<<i<< " before ";
			i++;
			std::cout<<i<<" after ";
			std::cout<<"here ";
			std::cout<<argv[i]<< " what ";
			while((i < argc) && !(strcmp(argv[i], "-f")==0) && !(strcmp(argv[i], "-d")==0)) {
				std::cout<<"here"<<" "<< argv[i];
				std::string str = argv[i];
				shape.push_back(std::stoi(str));
				i++;
			} 
		} else if (strcmp(argv[i], "-f")==0) {
			std::cout<<i<< "file flag i "<<std::endl;
			i++;
			filename = argv[i];
			std::cout<<filename<<" got the filename ";
			//TODO: make sure it ends in .npy
		} else if (strcmp(argv[i],"-d")==0) {
			i++;
			dtype = argv[i]; 
		}
	}
	if ((dtype == "") && (shape.size() == 0)) {
		//TODO: parse header
	}
	// TODO: get the right datatype
	int bufsize;
	int rank, nprocs;
	MPI_File fh;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Offset FILESIZE;
	MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_get_size(fh, &FILESIZE);
	//TODO: think about if this should be shape*shape*shape*dtype to validate the size of the header
	//TODO: fix this
	bufsize = (FILESIZE-128)/nprocs;
	int nints = bufsize/sizeof(long long int);
	//TODO: can this be changed to the MPI type ????
	long long int buf[nints];
	

	// reading in TODO: change the variables to not be hardcoded
	MPI_File_seek(fh, (rank*bufsize)+128, MPI_SEEK_SET);
	MPI_File_read(fh, buf, nints, MPI_LONG_LONG, &status);
	MPI_File_close(&fh);
	
	if (debug) {
		std::ofstream output;
		std::string outname = "output";
		std::string out;
		std::stringstream ss;
		ss << rank;
		out = ss.str();
		outname.append(out);
		outname.append(".txt");
		output.open(outname.c_str());
		for (int iter = 0; iter < nints; iter++) {
			output << buf[iter] << " ";
		}
		output.close();
	}
	MPI_Finalize();
	return 0;
}
