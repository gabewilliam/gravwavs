#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <fftw3.h>
#include "gwReadWrite.h"//print function
#include "gwDataTypes.h" //for Template struct

//I don't use these but you could adapt it so it does if you so desire
#define REAL(z, i) (z[i][0])
#define IMAG(z, i) (z[i][1])

int main(){

	//import signals
	std::vector<Signal> signals;

	if (!loadSignals("signal.csv",&signals,csv)){
		std::cout<<"Signal file not found"<<std::endl;
		return 0;
	}

	std::cout<<"Successfully imported signal."<<std::endl;
	
	std::vector< Signal > printer;
	
	//could extend to multiple signals later
	int N=signals[0].waveform[1].size();

	//import filters
	std::vector< Template > filters;

	//load templates
	if (!loadTemplates("templates.csv",&filters,csv)){
		std::cout<<"Template file not found"<<std::endl;
		return 0;
	}

	std::cout<<"Successfully imported filters."<<std::endl;

	//create output
	vec_d output(N);	
	vec_d outputB(2*N-1);

	//single filter
	vec_d inputFilter(N);

	//create complex arrays for use in transforms
	fftw_complex * signal, * filter, * signalOutput, * filterOutput, * product, * convolution;

	//allocate space for arrays
	//do not initialise them until after creation of plan
	signal=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));
	filter=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));
	signalOutput=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));
	filterOutput=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));	
	product=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));
	convolution=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(2*N-1));

	//create plans for each transform
	fftw_plan planSig=fftw_plan_dft_1d(2*N-1,signal,signalOutput,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan planFil=fftw_plan_dft_1d(2*N-1,filter,filterOutput,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan planCon=fftw_plan_dft_1d(2*N-1,product,convolution,FFTW_BACKWARD,FFTW_ESTIMATE);

	//set signal
	for (int i=0; i<2*N-1; i++){
		if (i<N){
			signal[i][0]=signals[0].waveform[1][i];
		}
		else {
			signal[i][0]=0.0;
		}
		signal[i][1]=0.0;
	}
	
	//forward transform of signal
	fftw_execute(planSig);

	//loop over all templates
	for(int j=0; j<filters.size(); ++j){
		//pad filter to be same length as signal
		while (filters[j].waveform[1].size()<N){
			filters[j].waveform[1].push_back(0.0);
		}
		//set real and imaginary parts of complex filter
		for (int i=0; i<2*N-1; i++){
			if (i<N) {
				filter[i][0] = filters[j].waveform[1][i];
			}
			else {
				filter[i][0]=0.0;
			}
			filter[i][1]=0.0;

		}

		//forward transform of filter
		fftw_execute(planFil);
	
		//calculate product of real and imaginary parts
		for (int i=0; i < 2*N-1; i++) {
			product[i][0]=signalOutput[i][0]*filterOutput[i][0];
			product[i][1]=-signalOutput[i][1]*filterOutput[i][1];
		}
	
		//carry out reverse transform
		fftw_execute(planCon);
	
		//set up output - needs fixing so length is back to N if possible - will look into
		//can play around with real parts/absolute values of convolution elements here
		//convolution[i][0]   or  sqrt(pow(convolution[i][0],2)+pow(convolution[i][1],2))

		for (int i=0; i<N; i++) {
			output[i]=convolution[i][0]+convolution[2*N-1-i][0];
		}
		
		//output
		Signal convolOutput;
		for (int i=0; i<N; i++){
			convolOutput.waveform[0].push_back(signals[0].waveform[0][i]);
			convolOutput.waveform[1].push_back(output[i]);
		}
		printer.push_back(convolOutput);			
	
		//reset filter, product etc to zero
		for (int i=0; i<2*N-1; i++) {
			filter[i][0]=0.0;
			filter[i][1]=0.0;
			filterOutput[i][0]=0.0;
			filterOutput[i][1]=0.0;
			product[i][0]=0.0;
			product[i][1]=0.0;
			convolution[i][0]=0.0;
			convolution[i][1]=0.0;
		}

		//progress output
		std::cout<<j<<std::endl;
	}
	
	//save to file
	if(!saveSignals("testout.csv",&printer,csv)){
		std::cout<<"Save failed."<<std::endl;
		return 0;
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
