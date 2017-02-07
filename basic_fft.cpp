#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

//Function to load the data from file "filename" to the data[] pointer
bool data_read(int* n, double* data[])
{
	//User input for filename
	std::string filename;
	std::cout << "Enter input file name..." << std::endl;
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

	//Counting number of rows before actual read. This is so the data arrays may be initialized.
	//Could use std::vector but gsl requires input of double[].
	int N=0;

	//Filename used as dummy as it is no longer needed
	while(std::getline(in_file, filename))
		N++;

	*n = N;
	data[0] = new double[N];
	data[1] = new double[N];

	//Returns to the start of the file for actual data read
	in_file.clear();
	in_file.seekg(0, std::ios::beg);

	//Actually Reading the data
	for(int i=0; i<N; i++)
	{
		in_file >> data[0][i];
		in_file >> data[1][i];
	}

	in_file.close(); 
	
	return true;
}

void data_write(int N, double* data[])
{
	//User input for output file name
	std::string filename;
	std::cout << "Enter output file name..." << std::endl;
	std::cin >> filename;

	//setting up out stream
	std::ofstream out_file;
	out_file.open(filename.c_str());

	//Writing results
	for(int i=0; i<N; i++)
	{
		out_file << std::setprecision(17) << data[0][i] << "\t" 
				 <<  std::setprecision(17) << data[1][i] << std::endl;
	}
}

int main ()
{
	//Setting up for data_read
	double* data[2];
	int N;

	//If data_read returns false the program exits
	if(!data_read(&N, data))
	{
		std::cout << "Failed to load data" << std::endl;
		return 0;
	}									

	//Allocateing a work space and look up tables for the gsl fft function
	gsl_fft_real_workspace* realWS = gsl_fft_real_workspace_alloc(N);
	gsl_fft_real_wavetable* realWT = gsl_fft_real_wavetable_alloc(N);

	gsl_fft_real_transform(data[1], 1, N, realWT, realWS);

	//spectral width is 1/dt where dt is the time spacing between data points
	double sw = N/(2*data[0][N-1]);	

	//new data set for outputing frequency data
	double freq[N/2];
	double real[N/2];
	double imag[N/2];

	for(int i=0; i<N/2; i++)
	{
		//Unpacking data
		freq[i] = i*sw/N;
		real[i] = data[1][i];
		imag[i] = data[1][N/2 + i];
	}

	//Packing data for output
	double* out_data[2];

	out_data[0] = freq;
	out_data[1] = real;

	//Input data will always be real so second half of output
	//(negative frequencies) may be ignored
	for(int i=0; i<N/2; i++)
	{
		out_data[0][i] = i*sw/N;
		out_data[1][i] = data[1][i];
	}

	data_write(N/2, out_data);

	return 0;
}
