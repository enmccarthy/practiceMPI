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
// dtype should be an int representing number of bytes of possibly MPI_data_type (*come back to this*)
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
	for (int i =0; i <argc; i++) { 
		if (strcmp(argv[i], "-s")==0) {
			i++;
			while((i < argc) 
					&& !(strcmp(argv[i], "-f")==0) 
					&& !(strcmp(argv[i], "-d")==0)) {
				std::string str = argv[i];
				shape.push_back(std::stoi(str));
				i++;
			} 
		} else if (strcmp(argv[i], "-f")==0) {
			filename = argv[++i];
			//TODO: make sure it ends in .npy
		} else if (strcmp(argv[i],"-d")==0) {
			dtype = argv[++i]; 
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
	//TODO: can this be changed to the MPI type ????, probably fix w templating
	long long int buf[nints];
	//TODO: maybe this should be changed to have a shape size of 1
	// this should probably be changed I dont like how the positions change with channel
	// but it should probably match the other implementation and is trivial here
	int channel = 1;
	int numSamples = 0;
	int x = shape[0];
	int y = shape[1];
	if (shape.size()==3) {
		numSamples = shape[2];
	} else if (shape.size() ==4) {
		channels = shape[2];
		numSamples = shape[3];
	}
	int headerlen = 128;
	// reading in TODO: change the variables to not be hardcoded
	// this is going to need to be in the for loop
	// should x and y be in bytes
	int seekvalue = 0;
	for (int iterS = 0; iterS < numSamples; iterS++) {
		for (int iterC = 0; iterC < channels; iterC++) {
			for (int iterY = 0; iterY < y; iterY++) {
				// this needs to be the full x or y 
				seekvalue = (numsamples*x*y*channels) + (iterC*x*y) + (iterY*x)
				MPI_File_seek(fh, seekvalue+headerlen, MPI_SEEK_SET);
				// I will need to think about buf, possibly pass a pointer
				// to the position of buf 
				// this should be the split x and y just advance the x pointer 
				MPI_File_read(fh, buf[x*y], x, MPI_LONG_LONG, &status);
			}
		}
	}
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
// Things to think about
// Templating the type depending on what is in the file 
// how to get the specific shape divide throughout the sample
// just split off of number of processes/nodes? this will probably be dictated elsewhere 
// 2D trivial, just split based off # of processes
// 3D 
// Odd case 11x7x30
//  ___________
// |     |     |   topL is 6x4
// |     |     |   bottomL is 6x3 	
// |     |     |   tobR is 5x4
// |_____|_____|   bottomR is 5x3
// |     |     |   int division,num of proc and subtraction should help get this
// |     |     |   reading bytes in it would be 6xnumbytes to read in then next line
// |_____|_____|   which is 5xnumbytes away, and do this 4 times
// 				   Next you need to read in the 3rd dimension which means skip the whole sec
// 				   ond half which would be 11x3xnumbytes for the number in the 3rd dim
// 4D --- each dimension adds a for loop will need to think more about it
// what would 4D look like and what is the important part to split
// need to think about how to do this for n dimensions without knowing the number of 
// dimensions before hand
// need to benchmark at somepoint
// x = 11
// y = 7 
// nodes = 2
// the options to give x and y, assume this has been optimized else where or
// just divide for now
//
// portionX = x/nodes
// startX = portionX*rank
// how to deal with the odd
// idea: structure as 4d so we have an x,y,z and layers where z is the number of datapoints
