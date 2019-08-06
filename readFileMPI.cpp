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
			//TODO: make sure it ends in .npy
		} else if (strcmp(argv[i],"-d")==0) {
			dtype = argv[++i]; 
		} else if (strcmp(argv[i], "-w")==0) {
			wordsize = std::stoi(argv[++i]); 
		}
	}
	int bufsize;
	int rank, nprocs;
    int numsamples = 1000;
	MPI_File fh;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   // int nux = 0;
   for(int nux = 0; nux<(numsamples/(nprocs/2)); nux++) {
        filename = files[((rank/2)+nux)%12];
	    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
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
	    //std::cout<<iterX<<" iterX \n";
	    //std::cout<<iterY<<"iterY \n";
	    //std::cout<<xPerNode<<" xPerNode \n";
	    //std::cout<<yPerNode<<" yPerNode \n";

	    //std::cout<<xPerNode <<"xPerNode \n";
	    //std::cout<<"buf size " << (xPerNode*yPerNode*zPerNode*sPerNode)<<"\n";
	    short int  buf[(xPerNode*yPerNode*zPerNode*sPerNode)];
	    short int *bufP = buf;
        // std::cout<<(xPerNode*yPerNode*zPerNode*sPerNode)<<"\n";
	    //TODO: should be set from header value to confirm 128 (as is most of the time)
    	int headerlen = 128;
	    int seekvalue = -1;
	    int gotoS = (iterS+sPerNode);
	    int gotoZ = (iterZ+zPerNode);
	    int gotoY = (iterY+yPerNode);
        //std::cout<<gotoS<<"gotoS \n";
	    //std::cout<<gotoY<<" gotoY \n";
        //std::cout<<gotoZ<<"gotoZ \n";
	    // what if y is split and it is all consecutive for x
        // TODO: rethink seek value
	    for (; iterS < gotoS; iterS++) {
           // std::cout<<iterS<<"\n";
           // std::cout<<gotoS<<"\n";
            for (int niterZ = iterZ; niterZ < gotoZ; niterZ++) {
                if(xlines == 1) {
                // this should just read in the chunk y * x
                    if (seekvalue == -1) {
                        seekvalue = ((iterS*x*y*z) + (niterZ*x*y) + (iterY*x) + iterX)*wordsize+headerlen;
                    } else {
                        seekvalue = ((x*y)-(xPerNode*yPerNode))*wordsize;
                    }
                    //std::cout<<iterS<<"here\n";
                    //std::cout<<iterZ<<"here1\n";
                    //std::cout<<seekvalue<<"seekvalue div by 2 \n";
                    //std::cout<<(xPerNode*yPerNode) <<"\n";
                    MPI_File_seek(fh, seekvalue, MPI_SEEK_CUR);
                   // std::cout<<*bufP<<"\n";
                    MPI_File_read(fh, bufP, (xPerNode*yPerNode), MPI_SHORT, &status);
                   // std::cout<<*bufP<<"\n";
                    bufP = bufP + (xPerNode*yPerNode);
                } else {
                    //FIX THIS FOR LOOP 
                    for (; iterY < gotoY; iterY++) {
				        // change seek value to match
				        // is this correct for higher dimensions, what if I split S 
				        if (seekvalue == -1) {
					        seekvalue = (((iterS*x*y*z)+(iterZ*x*y)+(iterY*x)+iterX)*wordsize)+headerlen;
					    // std::cout<<(seekvalue)/8<<" seek value \n";
				         } else {
					        //seekvalue = ((iterS*x*y*z)+(iterZ*x*y)+(iterY*x)+iterX)*wordsize;
					         if(xPerNode < x) {
					            seekvalue = xPerNode*wordsize;
					         } else {
					             seekvalue = 0;
					        }
				        }
				        //std::cout<<(seekvalue)<<" seek value \n";
				        MPI_Offset offset;
				        MPI_File_get_position(fh, &offset);
				        //std::cout<<"offset before" << offset <<"\n";
				        MPI_File_seek(fh,seekvalue, MPI_SEEK_CUR);
				        MPI_File_get_position(fh, &offset);
				        // std::cout<<"offset after" << offset << "\n";
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
