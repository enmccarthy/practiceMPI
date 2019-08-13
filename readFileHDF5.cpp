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
	int RANK = 4;
    std::cout<< "hello, its me \n";
    double start = MPI_Wtime();
	MPI_Init(&argc,&argv);
	bool debug = true;
	std::string filename = "";
	std::string dtype = "";
	std::vector<int> shape;
	//TODO change these to hdf5
	std::string files[] = {"/p/gpfs1/emccarth/test0.h5","/p/gpfs1/emccarth/test1.h5","/p/gpfs1/emccarth/test2.h5",
        "/p/gpfs1/emccarth/test3.h5","/p/gpfs1/emccarth/test4.h5",
        "/p/gpfs1/emccarth/test5.h5","/p/gpfs1/emccarth/test6.h5",
        "/p/gpfs1/emccarth/test7.h5","/p/gpfs1/emccarth/test8.h5",
        "/p/gpfs1/emccarth/test9.h5","/p/gpfs1/emccarth/test10.h5",
        "/p/gpfs1/emccarth/test11.h5","/p/gpfs1/emccarth/test12.h5",
        "/p/gpfs1/emccarth/test13.h5","/p/gpfs1/emccarth/test14.h5",
        "/p/gpfs1/emccarth/test15.h5","/p/gpfs1/emccarth/test16.h5",
        "/p/gpfs1/emccarth/test17.h5","/p/gpfs1/emccarth/test18.h5",
        "/p/gpfs1/emccarth/test19.h5","/p/gpfs1/emccarth/test20.h5",
        "/p/gpfs1/emccarth/test21.h5","/p/gpfs1/emccarth/test22.h5",
        "/p/gpfs1/emccarth/test23.h5","/p/gpfs1/emccarth/test24.h5",
        "/p/gpfs1/emccarth/test25.h5","/p/gpfs1/emccarth/test26.h5",
        "/p/gpfs1/emccarth/test27.h5","/p/gpfs1/emccarth/test28.h5",
        "/p/gpfs1/emccarth/test29.h5","/p/gpfs1/emccarth/test30.h5",
        "/p/gpfs1/emccarth/test31.h5","/p/gpfs1/emccarth/test32.h5",
        "/p/gpfs1/emccarth/test33.h5","/p/gpfs1/emccarth/test34.h5",
        "/p/gpfs1/emccarth/test35.h5","/p/gpfs1/emccarth/test36.h5",
        "/p/gpfs1/emccarth/test37.h5","/p/gpfs1/emccarth/test38.h5",
        "/p/gpfs1/emccarth/test39.h5","/p/gpfs1/emccarth/test40.h5",
        "/p/gpfs1/emccarth/test41.h5","/p/gpfs1/emccarth/test42.h5",
        "/p/gpfs1/emccarth/test43.h5","/p/gpfs1/emccarth/test44.h5",
        "/p/gpfs1/emccarth/test45.h5","/p/gpfs1/emccarth/test46.h5",
        "/p/gpfs1/emccarth/test47.h5","/p/gpfs1/emccarth/test48.h5",
        "/p/gpfs1/emccarth/test49.h5"};
	int wordsize;
	for (int i =0; i <argc; i++) { 
		if (strcmp(argv[i], "-s")==0) {
			i++;
			while((i < argc) && !(strcmp(argv[i], "-f")==0) 
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
	int numsamples = 2;
 
	MPI_Status status1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 
    MPI_Info mpi_info;
    MPI_Info_create(&mpi_info);

	//int nux = 0;
    std::cout<< "nprocs" << nprocs <<"\n";
    std::cout<< " " << (numsamples/(nprocs/2)) << "\n";
    for(int nux = 0; nux<(numsamples/(nprocs/2)); nux++) {
    //for(int nux = 0; nux < 1; nux++) {	
    // idk where this should be
		hid_t file;
        hid_t dataset;
        hid_t filespace;
        hid_t memspace;
        hid_t cparms;
        //TODO changese to be the correct dimensions
        //TODO HEREHERHERHEREHEREHERE 
        //does it want MPI_COMM_WORLD, MPI_COMM_SELF???? 
        //idk what H5Pcreate does, read about it
        
        hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD,mpi_info)
        herr_t status, status_n;  
		filename = files[((rank/2)+nux)%2];
        // open the file and the dataset
        // TODO how to get the dataset name?
        std::cout<<"file name " << filename.c_str() <<"\n";
        file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		hsize_t num_obj;
		char name[100];
		H5Gget_num_objs(file, &num_obj);
        if(num_obj ==1) {
            //std::cout<< "here \n";
		    H5Gget_objname_by_idx(file, 0,name, 10000);
		}
	    hid_t dapl_id;	
        dataset = H5Dopen(file, name);
        int rank_data;
        // Get dataset ran and dimensions
        // TODO ^ this will probably change the rest of the code
        
        filespace = H5Dget_space(dataset);
        rank_data = H5Sget_simple_extent_ndims(filespace);
		
        // properties list????????
        //std::cout<<"properties list \n";
        //std::cout<<"rank data " << rank_data << "\n";
        cparms = H5Dget_create_plist(dataset);

		// get total buf size
		hsize_t dims[shape.size()];
        hsize_t count[shape.size()];
        hsize_t offset[shape.size()];
		// unsure if I need this
        int rank_chunk;
        hsize_t i, j;
		
		int x = shape[0];
		int y = shape[1];
		int z = shape[2];
		int s = shape[3]; 
    
		int ylines = 2;
		int xlines = 1;
        int zlines = 1;
        int samplelines = 1;
        // per node per iter of the for loop
        
		int xPerNode = x/xlines;
		int yPerNode = y/ylines;
		int zPerNode = z/zlines;
		int sPerNode = s/samplelines;
        
        double data_out[yPerNode*xPerNode*zPerNode*sPerNode];
		status_n = H5Sget_simple_extent_dims(filespace, dims, NULL); 	
        
        std::cout<<" dims " << dims[0] << " " << dims[1] << " " << dims[2] << " "  <<dims[3] << "\n";
        // define memory space to read dataset???
        // idk what RANK should be
        dims[0] = dims[0]/2;
		memspace = H5Screate_simple(RANK, dims, NULL);
        //std::cout<<" memspace " << memspace << "\n";
		
        // read dataset back and "display"
        // status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, data_out);
        // this is setting where the iter should start for the for loop
        // as well as take care of the odd case theoretically
		if (xlines > 1) {
			// if an  odd number then add one to the ranks below
			// the remainder
			if(rank < (x%xPerNode)) {
				offset[0] = ((xPerNode+1)*(rank%2));
				xPerNode++; 
			} else {
                // in the even case x%xPerNode is 0
				offset[0] = ((xPerNode+1)*(x%xPerNode))+(xPerNode*(((rank%2)-(x%xPerNode))%xlines));
			}
		} else {
            offset[1] = 0;
        }
		// if things are in chunks ??????
		if (ylines > 1) {
			if((rank%2) < (y%yPerNode)) {
				offset[0] = ((yPerNode+1)*(rank%2));
				yPerNode++;
			} else {
				offset[0] = ((yPerNode+1)*(y%yPerNode)) + (yPerNode*(((rank%2)-(y%yPerNode))/xlines));
			}   
		} else {
            offset[0] = 0;
        }
		if (zlines > 1) {
			if(rank < (z%zPerNode)) {
				offset[2] = (zPerNode+1)*(rank%2);
				zPerNode++;
			} else {
				offset[2] = ((zPerNode+1)*(z%zPerNode)) + (zPerNode*(((rank%2)-(z%zPerNode))%zlines));
			}
		} else {
            offset[2] = 0;
        }
		if (samplelines > 1) {
			if (rank < (s%sPerNode)) {
				offset[3] = ((sPerNode+1)*(rank%2));
				sPerNode++;
			} else {
				offset[3] = ((sPerNode+1)*(s%sPerNode)) + 
                    (sPerNode*(((rank%2)-(s%sPerNode))%samplelines));
			}
		} else {
            offset[3] = 0; 
        }
        
        // I dont think count is correct
        hsize_t block[4];
        block[0] = yPerNode;
        block[1] = xPerNode;
		block[2] = zPerNode;
		block[3] = sPerNode;
        
	    count[0] = 1;
        count[1] = 1;
        count[2] = 1;
        count[3] = 1;     
        //maybe change this data type
		double  buf[(xPerNode*yPerNode*zPerNode*sPerNode)];
		// change this pointer
		double *bufP = buf;
		int seekvalue = -1;
		// start -> a starting location for the hyperslab
        // stride -> the number of elements to separate each element or block to be selected
        // count -> the number of elemenets or blocks to select along each dimension
        // block -> the size of the block selected from the dataspace 

        // what if y is split and it is all consecutive for x
        int stride[4];
        //I think stride can be null
        std::cout<< "offset " << offset[0] <<" "<< offset[1] << " "<< offset [2] << " " << offset[3] <<"\n";
        std::cout<< "block " << block[0] <<" "<< block[1] << " "<< block[2] << " " << block[3] <<"\n";
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, 
									  count, block);
        std::cout<<status<< " status\n";
        status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, 
						 H5P_DEFAULT, data_out);
        //std::cout<<"buf "<< buf[0]<<"\n";
        H5Dclose(dataset);
        H5Fclose(file);
	}
	double end = MPI_Wtime();
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
    //    std::cout<<(xPerNode*yPerNode*zPerNode*sPerNode) << "\n";
		//for (int iter = 0; iter < (xPerNode*yPerNode*zPerNode*sPerNode); iter++) {
		 // output << data_out[iter] << " ";
		//}
		output.close();
	}
	MPI_Finalize();
	std::cout<< "The process took " << end - start << " seconds to run. \n";
	return 0;}

