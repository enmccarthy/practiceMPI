#include <iostream>  
#include <fstream>
#include <vector> 
#include <sstream> 
int main(int argc, char *argv[]) 
{
	// read in byte by byte
	// first 6 bytes are (not needed)
	// next byte is version number of file format (needed)	
	// next two byes are the header length (needed, unsigned int, if 128 nothing unusual)
	// next header len bytes are describing arrays format 
		// end in new line
	// do I just need the word size and then can convert 
	// for now just get the shape of the data, I guess type is needed for this to understand how many
	// bytes to read in
	std::ifstream infile;
	infile.open("test.npy", std::ios::binary);
	char version;
	infile.seekg(6, std::ios::beg);
	infile.read(&version, 1);
	std::cout<<(int)version<< " version before casting"; 
	unsigned int headerLen;
	// there is possible a third version but it is not clear if the only difference is ascii and
	// if it defaults the rest to one or two
	if ((int)version == 1){
		unsigned short int len;
		unsigned short int holder = 10;
		// seek then read
		// or unsigned short int
		infile.seekg(1, std::ios::cur);
		std::cout<< "\n here \n"; 
		// this is either two of four depending on the version number
		infile.read((char *) (&len), 2);
		std::cout<< (int) len<< " len ";
		//std::cout<<(unsigned int) len << " casted len ";
		headerLen = holder  + len; 
		std::cout<< headerLen << " headerlen ";
	} else {
		// fix this
		unsigned int holder = 12;
		unsigned int len;
		infile.read(reinterpret_cast<char *> (&len), 4);
		headerLen = holder + len;
	} 
	//std::cout<<len<<" len ";  
	char header[headerLen];
	infile.read(header, headerLen);
	std::string headStr(header);
	// get this as a string instead for simplicity 
	if ((headerLen/64)==2) {
		std::cout<< headStr;
		size_t loc1, loc2;
		loc1 = headStr.find("descr");
		std::cout<< " " << headStr[loc1+11];
		int word_size = atoi(&headStr[loc1+11]);
		std::cout<< " " << word_size <<" word_size ";
		std::vector<int> input_shape;
		loc1 = headStr.find("(");
		loc2 = headStr.find(")");
		std::string numberStr;
		std::cout<< loc1 << " " << loc2 << " ";
		std::string prac = headStr.substr(loc1+1, (loc2-loc1));
		std::cout<<"\n" << prac << " \n";
		std::stringstream substring(prac);
		getline(substring, numberStr, ',');
		do {
			int num = atoi(numberStr.c_str());
			input_shape.push_back(num);
			getline(substring, numberStr, ',');
		} while (substring.good());
		for(int i=0; i < input_shape.size(); i++){
			std::cout<< " " << input_shape[i];
		}		
		
	}   
}
