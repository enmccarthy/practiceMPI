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
	
	//TODO: as input
	int channels = 1;
	int numSamples = 1;
	std::cout<<shape.size()<< " shape size \n"; 
	int x = shape[0];
	int y = shape[1];
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
	int ylines = 1;
	int xlines = nprocs;
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
		if(rank < (x%xPerNode)) {
			iterX = ((xPerNode+1)*rank);
			xPerNode++; 
		} else {
			iterX = ((xPerNode+1)*(x%xPerNode))+(xPerNode*(rank-(x%xPerNode)));
		}
	}
	if (ylines > 1) {
		if(rank < (y%yPerNode)) {
			iterY = ((yPerNode+1)*rank);
			yPerNode++;
		} else {
			iterY = ((yPerNode+1)*(y%yPerNode)) + (yPerNode*(rank-(y*yPerNode)));
		}
	}
	if (zlines > 1) {
		if(rank < (z%zPerNodel)) {
			iterZ = (zPerNode+1)*rank;
			zPerNode++;
		} else {
			iterZ = ((iterZ+1)*(z%zPerNode)) + (zPerNode*(rank-(z*zPerNode)));
		}
	}
	if (slines > 1) {
		if (rank < (s%sPerNode)) {
			iterS = ((sPerNode+1)*rank);
			sPerNode++;
		} else {
			iterS = ((iterS+1)*(s%sPerNode)) + (sPerNode*(rank-(s*sPerNode)));
		}
	}
	  
	long long int  buf[(xPerNode*yPerNode*zPerNode*sPerNode)];
	long long int *bufP = buf;
	//TODO: should be set from header value to confirm 128 (as is most of the time)
	int headerlen = 128;
	int seekvalue = 0;
	for (; iterS < (iterS+sPerNode); iterS++) {
		for (; iterZ < (iterZ+zPerNode); iterC++) {
			for (; iterY < (iterY+yPerNode); iterY++) {
				// change seek value to match
				// is this correct for higher dimensions, what if I split S
				// 
				if (seekvalue == -1) {
					seekvalue = (((iterS*x*y*z)+(iterZ*x*y)+(iterY*x))*wordsize)+headerlen;
				} else if (seekvalue == 0) {
					// when it is resetting
					seekvalue = ((iterS*x*y*z)+(iterZ*x*y)+(iterY*x))*wordsize;
				} else {
					seekvalue = ((iterY*x)+iterX)*wordsize;
				}
				// does reading move the pointer I dont think so
				std::cout<<(seekvalue)<<" seek value \n";
				MPI_File_seek(fh,seekvalue, MPI_SEEK_CUR);
				MPI_File_read(fh, bufP, xPerNode, MPI_LONG_LONG, &status);
				bufP = bufP + xPerNode;
			}
			seekvalue = 0;
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
