#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <stdio.h>

#include <vector>
#include <iostream>
#include <fstream>

#define REAL(z, i) (z[2*(i)])
#define IMAG(z, i) (z[2*(i)+1])

void InverseComplex(){

	//Setting up input stream
	std::ifstream in_file;
	in_file.open("noise.csv");
	
	//Checking if stream is setup succesfully
	if(in_file.fail())
	{
		std::cout << "Input file not found" << std::endl;
		//return false;
	}

	std::vector<double> freq;
	std::vector<double> ampRe;
	std::vector<double> ampIm;

	double a,b,c;

	while (!in_file.eof()){	

		in_file>>a>>b>>c;
		freq.push_back(a);
		ampRe.push_back(b);
		ampIm.push_back(c);
		
	}

	//______________________________

	//set up the needed factors
	int i = 0;
	int N = freq.size();
	//new arrays are 2*N as they contain imaginary data 
	double amp[2*N];
	double stride = 1;

	//initialise to 0
	for (i = 0; i < N; i++){
		REAL(amp, i) = 0.0;
		IMAG(amp, i) = 0.0;	
	}

	//set the array to contain real and imaginary
	for (i = 0; i < N; i++){
		REAL(amp, i) = ampRe[i];
		IMAG(amp, i) = ampIm[i];
	}


	//Allocateing a work space and look up tables for the gsl fft function
	gsl_fft_complex_workspace* compWS = gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N);

	gsl_fft_complex_inverse(amp, stride, N, compWT, compWS);
	
	gsl_fft_complex_workspace_free(compWS);
	gsl_fft_complex_wavetable_free(compWT);


	//convert frequency into time 
	double df = freq[1];
	double dt = 1/(N*df);
	
	double time[N];
	
	for (i=0; i<N; i++){
		time[i] = i*dt;
	}
	
	
	//output data to file
	std::ofstream newfile;
	newfile.open("TimeDomainComplexInversion.csv");

	for (i=0; i<N; i++){
		newfile<<time[i]<<","<<amp[i]<<std::endl;
	}
	
	return;

}