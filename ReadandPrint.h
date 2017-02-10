#ifndef ReadandPrint_H
#define ReadandPrint_H
#include <fstream>
#include <iostream>
#include <vector>

void outputWriter(std::vector<std::vector<double> > pointers){//This function writes any number of vectors to a file, to use it create a vector of the vectors you want to print and put it in the function, it will only print for the size of the first vector
	
	std::ofstream newfile;
	
	std::string filename;
	std::cout << "Enter Output file name..." << std::endl;
	std::cin >> filename;
	
	newfile.open(filename.c_str());
	for(std::vector<std::vector<double> >::iterator i=pointers.begin();i!=pointers.end();++i){
		
		for(size_t j=0;j<pointers[1].size();++j){
				newfile<<(*i)[j]<<",";
		}
		newfile<<std::endl;
	}
	newfile.close();
}

/*void inputReader(){
	
	std::string filename;
	std::cout << "Enter Output file name..." << std::endl;
	std::cin >> filename;
	
	std::ifstream in_file;
	in_file.open(filename);
	
}
*/
#endif // ReadandPrint_H