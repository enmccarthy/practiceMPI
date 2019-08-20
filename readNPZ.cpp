#include<cstring>
#include<string>
#include<vector>
#include<iostream>
#include<cstdio>
#include <fstream> 
#include <sstream>
int main() {
    std::string fname = "test1.npz";
    FILE* fp = fopen("/p/gpfs1/brainusr/datasets/cosmoflow/cosmoUniverse_2019_05_4parE/dim512_cube_nt4_npz/train/train_a1430483_int16.npz","rb");
    std::vector<char> local_header(30);
    size_t header_res = fread(&local_header[0],sizeof(char),30,fp);
    
    uint16_t name_len = *(uint16_t*) &local_header[26];
    std::string varname(name_len,' ');
    size_t vname_res = fread(&varname[0],sizeof(char),name_len,fp);
    std::cout<<varname<<"\n";
    uint16_t extra_field_len = *(uint16_t*) &local_header[28];
    if(extra_field_len > 0) {
        std::vector<char> buff(extra_field_len);
        size_t efield_res = fread(&buff[0],sizeof(char),extra_field_len,fp);
        if(efield_res != extra_field_len)
          throw std::runtime_error("npz_load: failed fread");      
    }
    char buffer[256];
    size_t res = fread(buffer, sizeof(char), 11, fp);
    if(res != 11)
       std::cout<<"hmmm\n";
    std::string head = fgets(buffer, 256, fp);
    std::cout<< ftell(fp) << "where I am \n";
    short int first;
    fread(&first, sizeof(short int), 1, fp);
    std::cout<<first<<" first int?? \n"; 
    return 0;
}
