#include <iostream>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_errno.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "ReadandPrint.h"//print function

#define REAL(z, i) (z[2*(i)])
#define IMAG(z, i) (z[2*(i)+1])

int main(){

	
	
	//Setting up input stream
	std::ifstream in_file;
	in_file.open("a.dat");
	
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
	in_file2.open("a.dat");
	
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
		in_file2 >> c >> d;	
		filterTime.push_back(c);
		filter.push_back(d);
		
	}
	
	//________________________________________________________________

	//set up the needed factors
	int i = 0;
	int N = filter.size();
	//new arrays are 2*N as they contain imaginary data 
	double signalData[2*N];
	double filterData[2*N];
	double convolutionData[2*N];
	double stride = 1;
	
	gsl_fft_complex_wavetable *wavetable;
	gsl_fft_complex_workspace *workspace;
	
	//initialise to 0
	for (i = 0; i < N; i++){
		REAL(signalData, i) = 0.0;
		IMAG(signalData, i) = 0.0;
		REAL(filterData, i) = 0.0;
		IMAG(filterData, i) = 0.0;
		REAL(convolutionData, i) = 0.0;
		IMAG(convolutionData, i) = 0.0;	
	}
	//set the real parts of the array to be the filter and signal
	for (i = 0; i < N; i++){
		REAL(signalData, i) = signal[i];
		REAL(filterData, i) = filter[i];
	}

	wavetable = gsl_fft_complex_wavetable_alloc(N);
	workspace = gsl_fft_complex_workspace_alloc(N);
	
	
	gsl_fft_complex_forward(signalData, stride, N, wavetable, workspace);
	gsl_fft_complex_forward(filterData, stride, N, wavetable, workspace);

	for (i = 0; i < N; i++){
		REAL(convolutionData, i) = REAL(filterData, i)*REAL(signalData, i);
		IMAG(convolutionData, i) = -IMAG(filterData, i)*IMAG(signalData, i);
	}
	
	gsl_fft_complex_backward (convolutionData, stride, N, wavetable, workspace);

	gsl_fft_complex_wavetable_free(wavetable);
	gsl_fft_complex_workspace_free(workspace);
	
	
	std::vector<double> convolution;
	
	for (i = 0; i < N; i++){
		convolution.push_back(REAL(convolutionData, i)*REAL(convolutionData, i)+IMAG(convolutionData,i)*IMAG(convolutionData,i));
	}
	
	
	
	std::vector<std::vector<double> > printer;//create 
	printer.push_back(signalTime);
	printer.push_back(convolution);
	
	std::string filename="Output2.csv";
	outputWriter(printer,filename);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


	return 0;
}
