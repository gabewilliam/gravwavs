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



void outputWriter(std::vector<double>Frequency,std::vector<double> fourierAmplitude){//This function writes out the two vectors to a file which can be plotted in MATLAB
	
	std::ofstream newfile;
	newfile.open("Outputforgraphing.csv");
	const int m=fourierAmplitude.size();

for(int j=0;j<m;j++){
	
	newfile<<Frequency[j]<<","<<fourierAmplitude[j]<<std::endl;
	//std::cout<<fTime[i]<<"Ampl  "<<fAmplitude[i]<<std::endl;
	//i++;
}
}
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
		std::vector<double> Convolution;
		std::vector<double> Time;
		for (size_t t = 0; t < N; t++){
	
			h[t] = 0;
			
			//multiplies the correct elements of the arrays together and adds the result ot the correct element of the convolution array
			//move the template to the left by equivalently moving the signal to the right
			//the only elements of f that are accesible are the elements that correspond to g(t=0) 
			//the number of accesible elements decreases as this is looped through
			for(size_t tau = 0; tau <= N-t; tau++){
				h[t] += (*signal)[t+tau]*(*filter)[tau];
				
			}
			Convolution.push_back(h[t]);
			Time.push_back(t);
		}
		

		outputWriter(Time,Convolution);
	
	}
	
	

}

int main(){	//main is used to load the data


	//User input for filename

	
	//Setting up input stream
	std::ifstream in_file;
	in_file.open("hidden.txt");
	
	//Checking if stream is setup succesfully
	if(in_file.fail())
	{
		std::cout << "Input file not found" << std::endl;
		//return false;
	}
	
	std::vector<double> signalTime;
	std::vector<double> signal;
	
	double a,b;
	
	
	
	while (!in_file.eof()){	
		//if(){
			in_file>>a>>b;
			signalTime.push_back(a);
			signal.push_back(b);
		//}else{
		//	loop = false;
		//}			
	}
	
	//User input for filename

	
	//Setting up input stream
	std::ifstream in_file2;
	in_file2.open("filter.txt");
	
	//Checking if stream is setup succesfully
	if(in_file2.fail())
	{
		std::cout << "Input file not found" << std::endl;
		//return false;
	}
	
	std::vector<double> filterTime;
	std::vector<double> filter;
	
	double c,d;

	while (!in_file2.eof()){	
		//if(){
			in_file2>>c>>d;
			filterTime.push_back(c);
			filter.push_back(d);
		//}else{
		//	loop = false;
		//}			
	}
	std::cout<<filter.size();
	convolution(&signal,&filter);	
	
	return 0;
}