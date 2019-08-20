#include <iostream> 
#include "mpi.h" 
#include <fstream> 
#include <sstream>
#include <vector>
#include <cstring>
#include <string>
#include <dirent.h>
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
    std::vector<std::string> files;
    std::string path = "/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/dim512_cube_nt4_npz/train/";
    struct dirent *entry;
    DIR *dir = opendir(path.c_str());
    if (dir == NULL) {
        return 0;
    }

    while ((entry = readdir(dir)) != NULL) {
        std::string temppath = path;
        files.push_back(temppath.append(entry->d_name));
    }

    closedir(dir);
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
	int bufsize;
	int rank, nprocs;
    files.erase(files.begin(), files.begin()+2); 
    int numsamples = files.size();
	//int numsamples = 1;
    MPI_File fh;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int nux = 0;
   for(int nux = 0; nux<(numsamples/(nprocs/4)); nux++) {
        
        filename = files[((rank/4)+nux)%numsamples];
        //filename = files[10];
        //std::cout<<filename<<" numsamples\n";
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	    //MPI_File_get_size(fh, &FILESIZE);
	    // get total buf size
	    int tempbufsize = 1;
	    for (int shapeI = 0; shapeI < shape.size(); shapeI++) {
		    bufsize = tempbufsize*shape[shapeI];
		    tempbufsize = bufsize;
	    } 
        //TODO change this hardcoded 2
	    bufsize = (tempbufsize)/4;
	    int nints = bufsize;
	
	    //TODO: as input
        int dims[4];
        dims[0] = shape[3];
        dims[1] = shape[2];
        dims[2] = shape[1];
        dims[3] = shape[0];
        // default split only along the x but have ways to split across other dimensions
        	   
	    short int  buf[bufsize];
        // std::cout<<(xPerNode*yPerNode*zPerNode*sPerNode)<<"\n";
	    //TODO: should be set from header value to confirm 128 (as is most of the time)
    	//int headerlen = 128;
	    int headerlen = 186; 
        int seekvalue = -1;
        int dist [4];
        MPI_Datatype filetype;
        dist[0] = MPI_DISTRIBUTE_BLOCK;
        dist[1] = MPI_DISTRIBUTE_BLOCK;
        dist[2] = MPI_DISTRIBUTE_BLOCK;
        dist[3] = MPI_DISTRIBUTE_BLOCK;

        int dargs[4];
        dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
        dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
        dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
        dargs[3] = MPI_DISTRIBUTE_DFLT_DARG;

        int psize[4];
       
       //These are line corresponding
        psize[0] = 1;
        psize[1] = 1;
        psize[2] = 2;
        psize[3] = 2;

        MPI_Type_create_darray(4,(rank%4), 4,dims,dist, dargs, psize, MPI_ORDER_C, MPI_SHORT, &filetype);   
        MPI_Type_commit(&filetype);
        MPI_Offset disp = headerlen;
        MPI_Datatype etype = MPI_SHORT;
        MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
        MPI_File_read(fh,&buf, nints, MPI_SHORT, MPI_STATUS_IGNORE);
					
	    MPI_File_close(&fh);
        MPI_Type_free(&filetype);
   }
	double end = MPI_Wtime();
    //std::cout<<*(bufP-xPerNode)<<"final buf \n";
    //std::cout<<nints<<"nints \n";
	if (false) {
		std::ofstream output;
		std::string outname = "output";
		std::string out;
		std::stringstream ss;
		ss << rank;
		out = ss.str();
		outname.append(out);
		outname.append(".txt");
		output.open(outname.c_str());
        //std::cout<<(xPerNode*yPerNode*zPerNode*sPerNode) << "\n";
	//	for (int iter = 0; iter < (xPerNode*yPerNode*zPerNode*sPerNode); iter++) {
	//	  output << buf[iter] << " ";
	//	}
		output.close();
	}
	MPI_Finalize();
	std::cout<< "The process took " << end - start << " seconds to run. \n";
	return 0;
}
