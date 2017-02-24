#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include <string>
#include <sstream>



#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i)+1])

using namespace std;
 
int main(){
	
	//Setting up input stream
	std::ifstream in_file;
	in_file.open("noise.csv");	
	
	//Checking if stream is setup succesfully
	if(in_file.fail())
	{
		cout << "Input file not found" << endl;
		return false;
	}
	
	vector<double> freq;
	vector<double> ampRe;
	vector<double> ampIm;
	
	string line;
	
	string a,b,c;
	
	double fr,re,im;	
	
	while(getline(in_file, line)){
		
		stringstream linestream(line);
		
		//loading parametres from first line of template
		getline(linestream,a, ',');
		getline(linestream,b, ',');
		getline(linestream,c, ',');
		
		//convert string to double
		stringstream converta(a);
		converta >> fr;
		stringstream convertb(b);
		convertb >> re;
		stringstream convertc(c);
		convertc >> im;
		
		//add to vector
		freq.push_back(fr);
		ampRe.push_back(re);
		ampIm.push_back(im);
		
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
		REAL(amp,i) = 0.0;
		IMAG(amp,i) = 0.0;	
	}

	//set the array to contain real and imaginary
	for (i = 0; i < N; i++){
		REAL(amp,i) = ampRe[i];
		IMAG(amp,i) = ampIm[i];
	}
	

	//Allocateing a work space and look up tables for the gsl fft function
	gsl_fft_complex_workspace* compWS = gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N);

	gsl_fft_complex_inverse(amp, stride, N, compWT, compWS);

	gsl_fft_complex_workspace_free(compWS);
	gsl_fft_complex_wavetable_free(compWT);


	//convert frequency into time 
	double df = freq[1];
	double dt = fabs(1/(N*df));
	
	double time[N];
	
	for (i=0; i<N; i++){
		time[i] = fabs(i*dt);
	}
		
	ofstream oFile;
	oFile.open("InverseNoise.txt", std::ios_base::app);
	
	for (i=0; i<N; i++){
		oFile << time[i] << "," << amp[i] << "\n" ;
	}
	
	cout << "Output finished and written to 'InverseNoise.txt' \n" << endl;
	
	oFile.close();
	
}
