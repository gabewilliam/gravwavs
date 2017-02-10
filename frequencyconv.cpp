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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "ReadandPrint.h"//print function



//computes the convolution of the data with the template
void convolution(std::vector<double> *signal, std::vector<double> *filter,std::vector<double> Time){ 

	
	size_t N = signal->size();
	//if the data sets are different sizes then doesnt compute the result
	
		//array for the convolution
		double h[N];
		std::vector<double> Convolution;

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
			
		}
		std::vector<std::vector<double> > printer;//create 
		printer.push_back(Time);
		printer.push_back(Convolution);
	
		
		outputWriter(printer);
	
	
	
	

}
void frequencyDomainConvolution (std::vector<double> *signal, std::vector<double> *filter,std::vector<double> signalTime){
	size_t N = (*signal).size();
	double stride = signalTime[N-1]/N;
	
	//computes the real FFT for the signal and filter
	gsl_fft_real_wavetable * real;
	gsl_fft_real_workspace * work;
	work = gsl_fft_real_workspace_alloc (N);
	real = gsl_fft_real_wavetable_alloc (N);
	
	gsl_fft_real_transform(&(*signal)[0], stride, N, real, work);
	gsl_fft_real_transform(&(*filter)[0], stride, N, real, work); 
	gsl_fft_real_wavetable_free (real);
	
	//initialises the convolution vector as the complex conjugate 
	std::vector <double> convolution;
	
	for(size_t i = 0; i < N/2; i++){ 
		convolution.push_back((*signal)[i]*(*filter)[i]);
	}
	for(size_t i = N/2; i < N; i++){
		convolution.push_back(-(*signal)[i]*(*filter)[i]);
	}
	
	//the inverse transform is taken to get back to the real convolution
	gsl_fft_halfcomplex_wavetable *hc;
	hc = gsl_fft_halfcomplex_wavetable_alloc (N);
	gsl_fft_halfcomplex_inverse (&convolution[0], stride, N, hc, work);
	gsl_fft_halfcomplex_wavetable_free (hc);
	gsl_fft_real_workspace_free (work);
	
	std::vector<std::vector<double> > printer;
	printer.push_back(signalTime);
	printer.push_back(convolution);
	
	outputWriter(printer);
	
}


int main(){	//main is used to load the data


	//User input for filename

	
	//Setting up input stream
	std::ifstream in_file;
	in_file.open("signalA.dat");
	
	//Checking if stream is setup succesfully
	if(in_file.fail())
	{
		std::cout << "Input file not found" << std::endl;
		return false;
	}
	
	std::vector<double> signalTime;
	std::vector<double> signal;
	
	double a,b;

	while (!in_file.eof()){	

			in_file>>a>>b;
			signalTime.push_back(a);
			signal.push_back(b);
		
	}
	
	
	//Setting up input stream
	std::ifstream in_file2;
	in_file2.open("filterB.dat");
	
	//Checking if stream is setup succesfully
	if(in_file2.fail())
	{
		std::cout << "Input file not found" << std::endl;
		return false;
	}
	
	std::vector<double> filterTime;
	std::vector<double> filter;
	
	double c,d;
	
	while (!in_file2.eof()){	

		in_file2>>c>>d;
		filterTime.push_back(c);
		filter.push_back(d);
		
	}
	
	while (filter.size()<signal.size()){
		filter.push_back(0.0);
	}
	int pushcounter;
	while (filter.size()>signal.size()){
		signal.insert(signal.begin() , 0.0);
		++pushcounter;
	}
	convolution(&signal,&filter,signalTime);	
	frequencyDomainConvolution (&signal,&filter,signalTime);
	return 0;
}