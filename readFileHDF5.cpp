#include <iostream> 
#include "mpi.h" 
#include <fstream> 
#include <sstream>
#include <vector>
#include <cstring>
#include <string>
#include "hdf5.h"
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
    //TODO change these to hdf5
    std::string files[] = {"/p/gpfs1/emccarth/test0.npy","/p/gpfs1/emccarth/test1.npy","/p/gpfs1/emccarth/test2.npy","/p/gpfs1/emccarth/test3.npy","/p/gpfs1/emccarth/test4.npy","/p/gpfs1/emccarth/test5.npy","/p/gpfs1/emccarth/test6.npy","/p/gpfs1/emccarth/test7.npy","/p/gpfs1/emccarth/test8.npy","/p/gpfs1/emccarth/test9.npy", "/p/gpfs1/emccarth/test10.npy", "/p/gpfs1/emccarth/test11.npy"};
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
		} else if (strcmp(argv[i],"-d")==0) {
			dtype = argv[++i]; 
		} else if (strcmp(argv[i], "-w")==0) {
			wordsize = std::stoi(argv[++i]); 
		}
	}

	int bufsize;
	int rank, nprocs;
    int numsamples = 1000;

 
    MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   
    // int nux = 0;
   for(int nux = 0; nux<(numsamples/(nprocs/2)); nux++) {
        // idk where this should be
	    hid_t file;
        hid_t dataset;
        hid_t filespace;
        hid_t memspace;
        hid_t cparms;
        //TODO change these to be the correct dimensions
        hsize_t dims[2];
        hsize_t chunk_dims[2];
        hsize_t col_dims[1];
        hsize_t count[2];
        hsize_t offset[2];

        herr_t status, status_n; 
        // TODO define dims
        // Should I keep all the dims like this 
	    int data_out[dim1][dim2];
        int chunk_out[2][5];
        int column[10];
        int rank_data, rank_chunk;
        hsize_t i, j;
        filename = files[((rank/2)+nux)%12];
        // open the file and the dataset
        // TODO how to get the dataset name?
        file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        dataset = H5Dopen(file, DATASETNAME);

        // Get dataset ran and dimensions
        // TODO ^ this will probably change the rest of the code
        filespace = H5Dget_space(dataset);
        rank_data = H5Sget_simple_extent_ndims(filespace);
        status_n = H5Sget_simple_extent_dims(filespace, dims, NULL); 
        printf("dataset rank %d, dimensions %lu x %lu\n", 
                rank, (unsigned long) (dims[0]), (unsigned long) (dims[1]));

        // properties list?????????
        cparms = H5Dget_create_plist(dataset);

        // check if dataset is chunked???? 
        if(H5D_chunked == H5Pget_layout(cparms)) {
            rank_chunk = H5Pget_chunk(cparms,2, chunk_dims);
            printf("chunk rank %d, dimensions %lu x %lu\n", rank_chunk, 
                    (unsigned long) (chunk_dims[0]), (unsigned long) (chunk_dims[1]));
        }

        // define memory space to read dataset???
        memspace = H5Screate_simple(RANK, dims, NULL);

        // read dataset back and "display"
        status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, data_out);

	    //MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	    //MPI_File_get_size(fh, &FILESIZE);
	    // get total buf size
	    int tempbufsize = 1;
	    for (int shapeI = 0; shapeI < shape.size(); shapeI++) {
		    bufsize = tempbufsize*shape[shapeI];
		    tempbufsize = bufsize;
	    } 
        //TODO change this hardcoded 2
	    bufsize = (tempbufsize)/2;
	    int nints = bufsize;
	
	    //TODO: as input
	    int x = shape[0];
	    int y = shape[1];
	    int z = shape[2];
	    int s = shape[3]; // this might be channels
        // default split only along the x but have ways to split across other dimensions

	    // these would be input of the ways you wanted it split
	    // look up conditional setting to make this prettier
	    int ylines = 2;
	    int xlines = 1;
        int zlines = 1;
        int samplelines = 1;
        // per node per iter of the for loop
        int xPerNode = x/xlines;
	    int yPerNode = y/ylines;
	    int zPerNode = z/zlines;
	    int sPerNode = s/samplelines;
	    int iterS = 0; 
	    int iterZ = 0;
	    int iterX = 0;
	    int iterY = 0;
        // this is setting where the iter should start for the for loop
        // as well as take care of the odd case theoretically
	    if (xlines > 1) {
		    // if an  odd number then add one to the ranks below
		    // the remainder
		    if(rank < (x%xPerNode)) {
			    iterX = ((xPerNode+1)*(rank%2));
			    xPerNode++; 
		    } else {
                // in the even case x%xPerNode is 0
			    iterX = ((xPerNode+1)*(x%xPerNode))+(xPerNode*(((rank%2)-(x%xPerNode))%xlines));
		    }
	    }
	    // if things are in chunks ????????????????
	    if (ylines > 1) {
		    if((rank%2) < (y%yPerNode)) {
			    iterY = ((yPerNode+1)*(rank%2));
			    yPerNode++;
		    } else {
			    iterY = ((yPerNode+1)*(y%yPerNode)) + (yPerNode*(((rank%2)-(y%yPerNode))/xlines));
	    	}   
	    }
	    if (zlines > 1) {
		    if(rank < (z%zPerNode)) {
			    iterZ = (zPerNode+1)*(rank%2);
			    zPerNode++;
		    } else {
			    iterZ = ((zPerNode+1)*(z%zPerNode)) + (zPerNode*(((rank%2)-(z%zPerNode))%zlines));
		    }
	    }
	    if (samplelines > 1) {
		    if (rank < (s%sPerNode)) {
			    iterS = ((sPerNode+1)*(rank%2));
		    	sPerNode++;
		    } else {
			    iterS = ((sPerNode+1)*(s%sPerNode)) + 
                    (sPerNode*(((rank%2)-(s%sPerNode))%samplelines));
		    }
	    }
        
        //TODO: I am going to need to do this and define the offset
        offset[0] = 0;
        offset[1] = 2;
        // I think this will be the pernode
        count[0] = 10;
        count[1] = 1;

	    short int  buf[(xPerNode*yPerNode*zPerNode*sPerNode)];
	    short int *bufP = buf;
	    int seekvalue = -1;
	    int gotoS = (iterS+sPerNode);
	    int gotoZ = (iterZ+zPerNode);
	    int gotoY = (iterY+yPerNode);
	    // what if y is split and it is all consecutive for x
        // TODO: rethink seek value
        // TODO: change the for loop
	    for (; iterS < gotoS; iterS++) {
            for (int niterZ = iterZ; niterZ < gotoZ; niterZ++) {
                if(xlines == 1) {
                    //TODO I dont need this if statement 
                    status =  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, 
                                        count, NULL);
                    status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace, 
                            H5P_DEFAULT, column);
                // this should just read in the chunk y * x
                    if (seekvalue == -1) {
                        seekvalue = ((iterS*x*y*z) + (niterZ*x*y) + (iterY*x) + iterX)*wordsize;
                    } else {
                        seekvalue = ((x*y)-(xPerNode*yPerNode))*wordsize;
                    }
                    MPI_File_seek(fh, seekvalue, MPI_SEEK_CUR);
                    MPI_File_read(fh, bufP, (xPerNode*yPerNode), MPI_SHORT, &status);
                    bufP = bufP + (xPerNode*yPerNode);
                } else {
                    for (; iterY < gotoY; iterY++) {
				        if (seekvalue == -1) {
					        seekvalue = (((iterS*x*y*z)+(iterZ*x*y)+(iterY*x)+iterX)*wordsize);
				         } else {
					         if(xPerNode < x) {
					            seekvalue = xPerNode*wordsize;
					         } else {
					             seekvalue = 0;
					        }
				        }
				        MPI_Offset offset;
				        MPI_File_get_position(fh, &offset);
				        MPI_File_seek(fh,seekvalue, MPI_SEEK_CUR);
				        MPI_File_get_position(fh, &offset);
				        MPI_File_read(fh, bufP, xPerNode, MPI_SHORT_INT, &status);
				        bufP = bufP + xPerNode;
			        }
		         }
	         }
        }
        int lp = *(bufP-xPerNode);
       // for(int i =0; i<10;i++) {
         //   std::cout<<*(bufP-i)<<"final buf \n";
       // }
	    MPI_File_close(&fh);
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
		//for (int iter = 0; iter < (xPerNode*yPerNode*zPerNode*sPerNode); iter++) {
		  //  output << buf[iter] << " ";
		//}
		output.close();
	}
	MPI_Finalize();
	std::cout<< "The process took " << end - start << " seconds to run. \n";
	return 0;
}
