#include "pe_gwLikelihood.h"
#include "gwReadWrite.h"
#include "gwSigGenBasic.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>


//AIM: reads in d(t) from signal, generates h(t) from known formula -> fourier transforms both to get d(f) and h(f), s(f) is hopefully already known and in frequency domain.
//Uses d(f), h(f) and s(f) using equation A20 in "Veich, Vecchio (2010)" to calculate likelihood using two parameters m1 and m2.
//Returns a NON-NORMALISED likelihood.

long double likelihood( double m1, double m2, double d, std::string signalFile ){

	Signal dataSignal;//Data signal (frequency domain)
	std::vector< Signal > idataSignal; //loadSignals requires std::vector

	loadSignals(signalFile, &idataSignal, csv); //Loads signal into dt

	dataSignal = idataSignal[0];//converts from vector of signals to signal

	vec_d modelAmplitude = ParameterFunction( m1, m2, d, dataSignal.waveform[0] );//Creates model function for given masses.

	vec_d noiseAmplitude = NoiseFunction( dataSignal.waveform[0] ); //creates noise probability function
	int n = dataSignal.waveform[0].size();//finds number of elements in freq domain signal

	long double sum = 0.0;
	long double SNR = 0.0;
	vec_d dataAmplitude = dataSignal.waveform[1];//splits signal into std::vectors for use in loop

	int nan=0;

	for( int i = 0; i < n; i++ ){ //evaluates the sum part in equation A20
		if(!(modelAmplitude[i]==modelAmplitude[i])) {
			modelAmplitude[i]=dataAmplitude[i];
			nan++;
		}

		sum += pow( std::abs( dataAmplitude[i] - modelAmplitude[i] ), 2 )/noiseAmplitude[i];
		SNR += pow( std::abs( dataAmplitude[i] ), 2 )/noiseAmplitude[i];
	}

	double T = 100.0;//total time elapsed
	SNR = sqrt(SNR/T);	
	
	long double pdh = (-2/T) * sum;//calculates p(d|h)

	return pdh;

}

vec_d ParameterFunction( double m1, double m2, double d, vec_d t ){ 

	Parameters PARAMS;
	Parameters *P = &PARAMS;

	NoisySignal NOISY;
	NoisySignal *N = &NOISY;

	Signal sig;
	Signal *C = &sig;

	Signal amp;
	Signal* A = &amp;

	vector<Signal> sigVect;
	vector<Signal> *S = &sigVect;

	gwSimulateDetection(m1,m2,d,0.0,0.0,10.0,100.0,P,A,C,N,S);

	return amp.waveform[1];
}

vec_d NoiseFunction( vec_d f ){
	
	int n = f.size();	
	vec_d sf;

	for( int i = 0; i < n; i++ ){
		sf.push_back( 1e-39 ); //PLACEHOLDER: INSERT NOISE FUNCTION IN PARENTHESIS
	}

	return sf;
}

