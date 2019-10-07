#ifndef READHDF5_CPP
#define READHDF5_CPP
#include <iostream> 
#include "mpi.h" 
#include <fstream> 
#include <sstream>
#include <vector>
#include <cstring>
#include <string>
//read file using MPI
#include "hdf5.h"
#include "dirent.h"
#ifdef H5_HAVE_PARALLEL
#warning yes
#else
#warning no
#endif
#ifdef H5F_DEBUG
#warning yes
#else
#warning no
#endif
int main(int argc, char *argv[]) 
{ 
  //TODO: replace this I might be able to get it from the metadata
  int RANK = 4;
  MPI_Init(&argc,&argv);
  std::string filename = "";
  std::vector<std::string> files;
  std::string path = "/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/hdf5/";
  struct dirent *entry;
    //std::string directory[] = {"21688988/","21922619/","21997469/","22059249/","22098324/","22309462/",
    //                            "21812950/","21929749/","22021490/","22074825/","22118427/"};
  std::string directory[] = {"21688988/"};
    //ski[ the first two bc they are . and ..
  for(int dirInd = 0; dirInd < 1; dirInd++) {
    std::string dirPath = path;
    dirPath.append(directory[dirInd]);
    int skip = 0;
    DIR *dir = opendir(dirPath.c_str());
    while ((entry = readdir(dir)) != NULL) {
      std::string temppath = dirPath;
      if(skip > 1) {
        files.push_back(temppath.append(entry->d_name));
        
       }
       skip++;
    }
    closedir(dir);
  }
  std::cout<<files.size()<< "files \n";
  int rank, nprocs;
  int numsamples = files.size(); 
  //double start = MPI_Wtime();
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 
  MPI_Info mpi_info = MPI_INFO_NULL;
  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  int n = nprocs/4;
  int ranks[n];
  int count = 0;
  for(int i = rank/4; i< (rank/4+4); i++){
    ranks[count] = i;
    count++;
  }
    //Group test
    //MPI_Group file_group;
    //MPI_Group_incl(world_group, 4, ranks, &file_group);
    //MPI_Comm file_comm;
    //MPI_Comm_create_group(MPI_COMM_WORLD, file_group, 0, &file_comm);
    
    //int nux = 0;
   MPI_Comm m_comm;
   int color = rank/4;
   MPI_Comm_split(MPI_COMM_WORLD, color, rank, &m_comm); 
   hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(fapl_id, m_comm, mpi_info);
   short int data_out[512*512*(512/4)*4];
   double start = MPI_Wtime();
  for(int nux = 0; nux<numsamples; nux+= nprocs/4) {
    hid_t file, dataset, filespace, memspace;
    herr_t status, status_n;  
    filename = files[nux];
    file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hsize_t num_obj;
    char name[100];
    
    H5Gget_num_objs(file, &num_obj);
    H5Gget_objname_by_idx(file, 0, name, 100);
          
    dataset = H5Dopen(file, name, H5P_DEFAULT);
    int rank_data;
    filespace = H5Dget_space(dataset);
    // can I use this to get the rank???/
    rank_data = H5Sget_simple_extent_ndims(filespace);

    // get total buf size
    hsize_t dims[RANK];
    hsize_t count[RANK];
    hsize_t offset[RANK];
    // unsure if I need this
    int rank_chunk;
    hsize_t i, j;
    
    status_n = H5Sget_simple_extent_dims(filespace, dims, NULL);
    int x = dims[0];
    int y = dims[1];
    int z = dims[2];
    int s = dims[3]; 
    
    int ylines = 1;
    int xlines = 1;
    int zlines = 4;
    int samplelines = 1;
        
    int xPerNode = x/xlines;
    int yPerNode = y/ylines;
    int zPerNode = z/zlines;
    int sPerNode = s/samplelines;
        
    //short int data_out[yPerNode*xPerNode*zPerNode*sPerNode];  
    hsize_t dims_local[4];
    dims_local[3] = sPerNode;
    dims_local[2] = zPerNode;
    dims_local[1] = yPerNode;
    dims_local[0] = xPerNode;
    memspace = H5Screate_simple(RANK, dims_local, NULL);
    int num_io_parts = 4;
    int odd_offset;
    if (xlines > 1) {
      // this is theoretically the odd case but it has not been tested yet
      // so everything with odd should be ignored
      if((rank%num_io_parts) < int(dims[0]%xPerNode)) {
        offset[0] = ((xPerNode+1)*(rank%num_io_parts));
      } else {
        // offset of this x dim for this rank;
        // in most cases odd_offset will be 0
        odd_offset = (xPerNode+1)*(dims[0]%xPerNode);
        offset[0] = odd_offset +(xPerNode*((rank%num_io_parts)%xlines));
      }
    } else {
      offset[0] = 0;
    }

    if (ylines > 1) {
      if((rank%num_io_parts) < int(dims[1]%yPerNode)) {
        offset[1] = ((yPerNode+1)*(rank%num_io_parts));
      } else {
        // offset of the y dim for this rank
        odd_offset = (yPerNode+1)*(dims[1]%yPerNode);
        offset[1] = odd_offset + (yPerNode*((rank%num_io_parts)/ylines));
      }
    } else {
      offset[1] = 0;
    }

    if (zlines > 1) {
      offset[2] = zPerNode*(rank%num_io_parts);
    } else {
      offset[2] = 0;
    }
    // I have channel dim splits 
    // all combos arent really tested
    //add the rest later
    offset[3] = 0;
 
    count[0] = 1;
    count[1] = 1;
    count[2] = 1;
    count[3] = 1;     
    // start -> a starting location for the hyperslab
    // stride -> the number of elements to separate each element or block to be selected
    // count -> the number of elemenets or blocks to select along each dimension
    // block -> the size of the block selected from the dataspace 
    //I think stride can be null
    //std::cout<< "offset " << offset[0] <<" "<< offset[1] << " "<< offset [2] << " " << offset[3] <<"\n";
    //std::cout<< "block " << block[0] <<" "<< block[1] << " "<< block[2] << " " << block[3] <<"\n";
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, 
                    count, dims_local);
    status = H5Dread(dataset, H5T_NATIVE_SHORT, memspace, filespace, 
                       H5P_DEFAULT, data_out);
      //std::cout<<data_out[10]<<"\n";
      //MPI_Comm_free(&file_com);
       
    H5Dclose(dataset);
    H5Fclose(file);
  }
  double end = MPI_Wtime();
    //MPI_Group_free(&world_group);
    //  MPI_Group_free(&file_group);
    //  MPI_Comm_free(&file_comm);
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
      //   output << data_out[iter] << " ";
    //}
    output.close();
  }
  MPI_Finalize();
  std::cout<< "Rank " << rank << " \n";
  std::cout<< "The process took " << end - start << " seconds to run. \n";
  return 0;}
#endif
