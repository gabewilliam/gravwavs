/*
 * program to find the convolution of two functions f(t) and g(t)
 * h(t) = [f*g](t).
 * g must be at the start of the array and the data must come after it.
 *
 *
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

//computes the convolution of the data with the template
void convolution(std::vector<double> *signal, std::vector<double> *filter){ 


	size_t N = signal->size();
	//if the data sets are different sizes then doesnt compute the result
	if (N != filter->size()){
		std::cerr << "unmatched data size"<< std::endl;	
	}
	else{
		//array for the convolution
		double h[N];
		
		
		for (size_t t = 0; t < N; t++){
	
			h[t] = 0;
			
			//multiplies the correct elements of the arrays together and adds the result ot the correct element of the convolution array
			//move the template to the left by equivalently moving the signal to the right
			//the only elements of f that are accesible are the elements that correspond to g(t=0) 
			//the number of accesible elements decreases as this is looped through
			for(size_t tau = 0; tau <= N-t; tau++){
				h[t] += (*signal)[t+tau]*(*filter)[tau];
			}
		}
		
		//saves the result to a file
		FILE *fout;
		fout = fopen("convolution.csv","w");
		//FILE *fout2;
		//fout2 = fopen("f.dat","w");
		//FILE *fout3;
		//fout3 = fopen("g.dat","w");
	
		for (size_t t = 0; t < N; t++){
			fprintf(fout, "%f,%lu\n",h[t],",",t);
			//fprintf(fout2, "%f,%lu\n",(*f)[t],t);
			//fprintf(fout3, "%f,%lu\n",(*g)[t],t);
		}
	
	}
	
	

}

int main(){	//main is used to load the data


	//User input for filename
	std::string filename;
	std::cout << "Enter signal file name..." << std::endl;
	std::cin >> filename;
	
	//Setting up input stream
	std::ifstream in_file;
	in_file.open(filename.c_str());
	
	//Checking if stream is setup succesfully
	if(in_file.fail())
	{
		std::cout << "Input file not found" << std::endl;
		return false;
	}
	
	std::vector<double> signalTime;
	std::vector<double> signal;
	
	double data;
	

	
	while (in_file.eof()){	
		//if(){
			in_file >> data;
			signalTime.push_back(data);
			in_file >> data;
			signal.push_back(data);
		//}else{
		//	loop = false;
		//}			
	}
	
	//User input for filename
	std::string filename2;
	std::cout << "Enter Filter file name..." << std::endl;
	std::cin >> filename2;
	
	//Setting up input stream
	std::ifstream in_file2;
	in_file2.open(filename2.c_str());
	
	//Checking if stream is setup succesfully
	if(in_file.fail())
	{
		std::cout << "Input file not found" << std::endl;
		return false;
	}
	
	std::vector<double> filterTime;
	std::vector<double> filter;
	
	

	while (in_file.eof()){	
		//if(){
			in_file2 >> data;
			filterTime.push_back(data);
			in_file2 >> data;
			filter.push_back(data);
		//}else{
		//	loop = false;
		//}			
	}

	convolution(&signal,&filter);	
	
	return 0;
}