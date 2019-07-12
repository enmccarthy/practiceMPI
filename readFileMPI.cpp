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
	//program name, file name (-f), dtype (-d), shape (-s) , word size (-w) 
	MPI_Init(&argc,&argv);
	bool debug = true;
	std::string filename = "";
	std::string dtype = "";
	std::vector<int> shape;
	int wordsize;
	for (int i =0; i <argc; i++) { 
		if (strcmp(argv[i], "-s")==0) {
			i++;
			while((i < argc) 
					&& !(strcmp(argv[i], "-f")==0) 
					&& !(strcmp(argv[i], "-d")==0)
					&& !(strcmp(argv[i], "-w")==0)) {
				std::string str = argv[i];
				shape.push_back(std::stoi(str));
				i++;
			} 
		} 
		if (strcmp(argv[i], "-f")==0) {
			filename = argv[++i];
			//TODO: make sure it ends in .npy
		} else if (strcmp(argv[i],"-d")==0) {
			dtype = argv[++i]; 
		} else if (strcmp(argv[i], "-w")==0) {
			wordsize = std::stoi(argv[++i]); 
		}
	}
	//if ((dtype == "") && (shape.size() == 0)) {
		//TODO: parse header
	//}
	int bufsize;
	int rank, nprocs;
	MPI_File fh;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Offset FILESIZE;
	MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_get_size(fh, &FILESIZE);
	// get total buf size
	std::cout<<FILESIZE<< " FILESIZE \n";
	int tempbufsize = 1;
	for (int shapeI = 0; shapeI < shape.size(); shapeI++) {
		bufsize = tempbufsize*shape[shapeI];
		tempbufsize = bufsize;
	} 
	bufsize = (tempbufsize)/nprocs;
	std::cout<<bufsize<< " bufsize \n";
	int nints = bufsize;
	std::cout<< nints << "nints \n";
	long long int  buf[nints];
	long long int *bufP = buf;
	
	//TODO: maybe this should be changed to have a shape size of 1
	// this should probably be changed I dont like how the positions change with channel
	// but it should probably match the other implementation and is trivial here
	int channels = 1;
	int numSamples = 1;
	std::cout<<shape.size()<< " shape size \n";
	// this needs to be relative to the number of processes 
	int x = shape[0];
	int y = shape[1];
	// something hardcoded
	int xPerNode = 1000;
	int yPerNode = 500;
	if (shape.size()==3) {
		numSamples = shape[2];
	} else if (shape.size() ==4) {
		channels = shape[2];
		numSamples = shape[3];
	}
	int ylines = 1;
	int xlines = 0;
	//TODO: should this be set from the header?
	int headerlen = 128;
	int seekvalue = 0;
	for (int iterS = 0; iterS < numSamples; iterS++) {
		for (int iterC = 0; iterC < channels; iterC++) {
			for (int iterY = (yPerNode*rank); iterY < ((yPerNode*rank)+yPerNode); iterY++) {
				// seek value for the file
				// try seek curr
				if (seekvalue == 0) {
					seekvalue = (yPerNode*rank)+headerlen;
				} else {
					seekvalue = (xPerNode*wordsize);
				}
				std::cout<<(seekvalue)<<" seek value \n";
				MPI_File_seek(fh,(seekvalue*rank), MPI_SEEK_CUR);
				if (xPerNode < x) {
					bufP = bufP + (xPerNode*rank);
				} else { 
					MPI_File_read(fh, bufP, xPerNode, MPI_LONG_LONG, &status);
				}
				bufP = bufP + xPerNode;
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
