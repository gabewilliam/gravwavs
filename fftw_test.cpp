#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <fftw3.h>
#include "ReadandPrint.h"//print function
#include "templates.h" //for Template struct

#define REAL(z, i) (z[i][0])
#define IMAG(z, i) (z[i][1])

int main(){

	//Setting up input stream for signal
	std::ifstream in_file;
	in_file.open("signalB.dat");
	
	//Checking if stream is setup succesfully
	if(in_file.fail())
	{
		std::cout << "Input file not found" << std::endl;
		return 0;
	}
	
	//transfer signal to vector
	std::vector<double> signalTime;
	std::vector<double> inputSignal;

	double a,b;

	while (!in_file.eof()){	

		in_file>>a>>b;
		signalTime.push_back(a);
		inputSignal.push_back(b);
	}

	//should loop through all templates and set up to use this
	std::vector<std::vector<double> > printer;
	printer.push_back(signalTime);
	int N=inputSignal.size();

	//import filters
	std::vector	<double> filterTime;
	std::vector<Template>filters;

	//load templates
	if (!load_templates("templates2.csv",&filters,csv)){
		std::cout<<"Template file not found"<<std::endl;
		return 0;
	}

	for(size_t j=0;j<filters.size();++j){
		while (filters[j].waveform.size()<inputSignal.size()){
			filters[j].waveform.push_back(0.0);
		}
	}

	//create output - currently messy
	std::vector<double> output(N);
	std::vector<double> outputB(2*N-1);

	//single filter - set up to loop over all
	std::vector<double> inputFilter(N);
	for (int i=0; i<N; i++) {
		inputFilter[i]=filters[0].waveform[i];
	}

	//create complex arrays for use in transforms
	fftw_complex * signal, * filter, * signalOutput, * filterOutput, * product, * convolution;

	//allocate space for arrays
	//do not initialise until after creation of plan
	signal=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));
	filter=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));
	signalOutput=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));
	filterOutput=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));	
	product=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));
	convolution=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));

	//create plans for each transform
	//plans can be re-used but I need to figure out how
	fftw_plan planSig=fftw_plan_dft_1d(2*N-1,signal,signalOutput,FFTW_FORWARD,FFTW_MEASURE);
	fftw_plan planFil=fftw_plan_dft_1d(2*N-1,filter,filterOutput,FFTW_FORWARD,FFTW_MEASURE);
	fftw_plan planCon=fftw_plan_dft_1d(2*N-1,product,convolution,FFTW_BACKWARD,FFTW_MEASURE);

	//set real and imaginary parts of all elements of signal and filter
	for (int i=0; i<2*N-1; i++){
		if (i<N) {
			signal[i][0]= inputSignal[i];
			filter[i][0] = inputFilter[i];
		}
		else {
			signal[i][0]=0.0;
			filter[i][0]=0.0;
		}
		signal[i][1]=0.0;
		filter[i][1]=0.0;
	}

	//carry out forward transforms
	fftw_execute(planSig);
	fftw_execute(planFil);
	
	//calculate product of real and imaginary parts
	for (int i=0; i < 2*N-1; i++) {
		product[i][0]=signalOutput[i][0]*filterOutput[i][0];
		product[i][1]=-signalOutput[i][1]*filterOutput[i][1];
	}
	
	//carry out reverse transform
	fftw_execute(planCon);
	
	//set up output - needs fixing so length is back to N if possible - will look into
	for (int i=0; i<2*N-1; i++) {
		outputB[i]=convolution[(i+N)%(2*N-1)][0];
	}

	//print - should use printer but it segfaulted once and I'm lazy
	FILE * fileout;
	fileout=fopen("testout2.csv","w");
	for (int i=0; i<N; i++){
		//real part used for now because it looks nice
		fprintf(fileout,"%g, ",convolution[i][0]/*sqrt(pow(convolution[i][0],2)+pow(convolution[i][1],2))*/);
	}

	//free memory
	fftw_destroy_plan(planSig);
	fftw_destroy_plan(planFil);
	fftw_destroy_plan(planCon);
	fftw_free(signal);
	fftw_free(filter);
	fftw_free(signalOutput);
	fftw_free(filterOutput);
	fftw_free(product);
	fftw_free(convolution);

	return 0;
}
