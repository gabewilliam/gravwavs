#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "pe_gwLikelihood.h"
#include "gwSigGen.h"
#include "pe_gwSaveToFile.h"

long double likelihood( double m1, double m2, double d, std::string signalFile, Signal dataSignal, vec_d noiseAmplitude){

	//Finds total time elapsed
	double T = 1/(dataSignal.waveform[0][1]-dataSignal.waveform[0][0]);

	//Creates model function for given masses
	vec_d modelAmplitude = ParameterFunction( m1, m2, d, dataSignal.waveform[0], T );

	//Finds number of elements in freq domain signal
	int n = dataSignal.waveform[0].size();

	long double pdh = 0.0;
	long double SNR = 0.0;

	//Extracts data amplitudes in a single vector
	vec_d dataAmplitude = dataSignal.waveform[1];

	//Calculates the sum for the likelihood function and the signal to noise ratio of the data
	for( int i = 0; i < n; i++ ){

		pdh += pow( std::abs( dataAmplitude[i] - modelAmplitude[i] ), 2 )/noiseAmplitude[i];	
		SNR += pow( std::abs( dataAmplitude[i] ), 2 )/noiseAmplitude[i];
	}

	SNR = sqrt(SNR/T);	
	pdh = (-2/T) * pdh;

	return pdh;
}

//Function that computes the model signal
vec_d ParameterFunction( double m1, double m2, double d, vec_d f, double T ){ 

	//Finds number of frequency values in the signal
	int n = f.size();

	vec_d modelAmplitude;	

	//Sets up the necessary variables to create a signal
	parameters PARAMS;
	parameters *P = &PARAMS;

	setFundamentalParameters(0.0,T,0.0,10.0,m1,m2,d,0.0,0.0,0.0,P);
	
	complex<double> GRAVWAV = 0.0;
	complex<double> *G = &GRAVWAV;

	//Computes amplitudes for a signal of given parameter values
	for(int i=0; i<n; i++){
		modelAmplitude.push_back(updatedAmplitude(P,f[i]/C_CONST,G));
	}
	
	return modelAmplitude;
}

