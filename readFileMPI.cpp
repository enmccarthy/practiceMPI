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
	double start = MPI_Wtime();
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
	int tempbufsize = 1;
	for (int shapeI = 0; shapeI < shape.size(); shapeI++) {
		bufsize = tempbufsize*shape[shapeI];
		tempbufsize = bufsize;
	} 
	bufsize = (tempbufsize)/nprocs;
	int nints = bufsize;
	
	//TODO: as input
	int channels = 1;
	int numSamples = 1; 
	int x = shape[0];
	int y = shape[1];
	int z = 1;
	int s = 1; // this might be channels
	// in this case the samples are in different files
	if (shape.size()==3) {
		numSamples = shape[2];
	} else if (shape.size() ==4) {
		channels = shape[2];
		numSamples = shape[3];
	}
    // default split only along the x but have ways to split across other dimensions
    // first x%xPerNode partitions will have one extra
    // the rest will have xPerNode amount

	// these would be input of the ways you wanted it split
	// look up conditional setting to make this prettier
	int ylines = 2;
	int xlines = 2;
    int zlines = 1;
    int samplelines = 1;
    int xPerNode = x/xlines;
	int yPerNode = y/ylines;
	int zPerNode = z/zlines;
	int sPerNode = s/samplelines;
	int iterS = 0; 
	int iterZ = 0;
	int iterX = 0;
	int iterY = 0;
	// can also check odd here and maybe add one to xPerNode??? 
	if (xlines > 1) {
		// if an  odd number then add one to the ranks below
		// the remainder
		if(rank < (x%xPerNode)) {
			iterX = ((xPerNode+1)*rank);
			xPerNode++; 
		} else {
			iterX = ((xPerNode+1)*(x%xPerNode))+(xPerNode*(rank-(x%xPerNode))*(rank%2));
		}
		//std::cout<<"iterX after \n" << iterX<< "\n";
	}
	// if things are in chunks ????????????????
	if (ylines > 1) {
		if(rank < (y%yPerNode)) {
			iterY = ((yPerNode+1)*rank);
			yPerNode++;
		} else {
			iterY = ((yPerNode+1)*(y%yPerNode)) + (yPerNode*(rank-(y%yPerNode)));
		}
	}
	if (zlines > 1) {
		if(rank < (z%zPerNode)) {
			iterZ = (zPerNode+1)*rank;
			zPerNode++;
		} else {
			iterZ = ((iterZ+1)*(z%zPerNode)) + (zPerNode*(rank-(z%zPerNode)));
		}
	}
	if (samplelines > 1) {
		if (rank < (s%sPerNode)) {
			iterS = ((sPerNode+1)*rank);
			sPerNode++;
		} else {
			iterS = ((iterS+1)*(s%sPerNode)) + (sPerNode*(rank-(s%sPerNode)));
		}
	}
	std::cout<<iterX<<" iterX \n";
	std::cout<<iterY<<"iterY \n";
	std::cout<<xPerNode<<" xPerNode \n";
	std::cout<<yPerNode<<" yPerNode \n";

	//std::cout<<xPerNode <<"xPerNode \n";
	//std::cout<<"buf size " << (xPerNode*yPerNode*zPerNode*sPerNode)<<"\n";
	double  buf[(xPerNode*yPerNode*zPerNode*sPerNode)];
	double *bufP = buf;
	//TODO: should be set from header value to confirm 128 (as is most of the time)
	int headerlen = 128;
	int seekvalue = -1;
	int gotoS = (iterS+sPerNode);
	int gotoZ = (iterZ+zPerNode);
	int gotoY = (iterY+yPerNode);
	//std::cout<<gotoY<<" gotoY \n";
	// what if y is split and it is all consecutive for x
	for (; iterS < gotoS; iterS++) {
		for (; iterZ < gotoZ; iterZ++) {
			for (; iterY < gotoY; iterY++) {
				// change seek value to match
				// is this correct for higher dimensions, what if I split S 
				if (seekvalue == -1) {
					seekvalue = (((iterS*x*y*z)+(iterZ*x*y)+(iterY*x)+iterX)*wordsize)+headerlen;
					std::cout<<(seekvalue)/8<<" seek value \n";
				} else if (seekvalue == -2) {
					// when it is resetting, is this right?????
					seekvalue = ((iterS*x*y*z)+(iterZ*x*y)+(iterY*x)+iterX)*wordsize;
				//	std::cout<< "should be resetting " << seekvalue/8<<"\n";
				} else {
					if(xPerNode < x) {
						seekvalue = xPerNode*wordsize;
					} else {
						seekvalue = 0;
					}
				}
				// does reading move the pointer I dont think so
				//std::cout<<(seekvalue)<<" seek value \n";
				MPI_Offset offset;
				MPI_File_get_position(fh, &offset);
				//std::cout<<"offset before" << offset <<"\n";
				MPI_File_seek(fh,seekvalue, MPI_SEEK_CUR);
				MPI_File_get_position(fh, &offset);
				//std::cout<<"offset after" << offset << "\n";
				MPI_File_read(fh, bufP, xPerNode, MPI_DOUBLE, &status);
				bufP = bufP + xPerNode;
			}
			seekvalue = -2;
		}
	}
	MPI_File_close(&fh);
	double end = MPI_Wtime();
	
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
	std::cout<< "The process took " << end - start << " seconds to run. \n";
	return 0;
}
