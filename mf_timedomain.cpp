#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>

#define pi (atan(1)*4)

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

	//If data_read return false the program exits
	if(!data_read(&N, data))
	{
		std::cout << "Failed to load data" << std::endl;
		return 0;
	}							

	//Number of trial templates
	int I =1000;

	double h = 1e-3;
	double result[I];
	double f[I];
	double ip;
	
	//Looping through all templates
	for(int i=0; i < I; i++)
	{
		ip = 0;
		f[i] = i*500/I;	

		for(int n=0; n<N; n++)
		{
			ip += (h/2)*(sin(2*pi*f[i]*data[0][n])*data[1][n] + sin(2*pi*f[i]*data[0][n+1])*data[1][n+1]);
		}
		
		result[i] = 4*ip;
	}

	//Packaging data for witing
	double* data_out[2];
	data_out[0] = f;
	data_out[1] = snr;	

	data_write(I,data_out);

	return 0;
}
