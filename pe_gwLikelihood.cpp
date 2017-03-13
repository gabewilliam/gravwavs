#include "pe_gwLikelihood.h"
#include "gwReadWrite.h"
//#include "gwFFT.h"
//#include "generate.h"
#include "gwSigGen.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>


//AIM: reads in d(t) from signal, generates h(t) from known formula -> fourier transforms both to get d(f) and h(f), s(f) is hopefully already known and in frequency domain.
//Uses d(f), h(f) and s(f) using equation A20 in "Veich, Vecchio (2010)" to calculate likelihood using two parameters m1 and m2.
//Returns a NON-NORMALISED likelihood.

long double likelihood( double m1, double m2, double d, std::string signalFile ){

	Signal dt;//Data signal (time domain)
	std::vector< Signal > idt; //loadSignals requires std::vector
	Signal df;//Data signal (freq domain)
	Signal ht;//Model signal (time domain)
	Signal hf;//Model signal (freq domain)

	loadSignals("signal.csv", &idt, csv); //Loads signal into dt

	dt = idt[0];//converts from vector of signals to signal

	vec_d ht2 = ParameterFunction( m1, m2, d, dt.waveform[0] );//Creates model function for given masses.

	ht.waveform[0] = dt.waveform[0]; //ht uses same time scale as dt
	ht.waveform[1] = ht2;	//sets ht equal to vector produced by function
	
	//Extractor Fourier = Extractor(); //performs fourier transforms on each of the time domain signals 
	//Fourier.setSignal(&dt);
	//Fourier.fft(&df);
	//Fourier.setSignal(&ht);
	//Fourier.fft(&hf);
	

	vec_d sf = NoiseFunction( dt.waveform[0] ); //creates noise probability function
	int n = dt.waveform[0].size();//finds number of elements in freq domain signal

	long double sum = 0.0;
	long double SNR = 0.0;
	vec_d vdf = dt.waveform[1]; //splits signal into std::vectors for use in loop
	vec_d vhf = ht.waveform[1];

	int nan=0;

	for( int i = 0; i < n; i++ ){ //evaluates the sum part in equation A20
		//if (i==500) std::cout << vdf[i] << "\t" << vhf[i] << std::endl;	
		if(!(vhf[i]==vhf[i])) {
			vhf[i]=vdf[i];
			nan++;
		}

		sum += pow( std::abs( vdf[i] - vhf[i] ), 2 )/sf[i];
		SNR += pow( std::abs( vdf[i] ), 2 )/sf[i];
	}
	
	//std::cout<<"SNR="<<SNR<<std::endl;

	vec_d t = dt.waveform[0];//takes the time vector
	int n2 = t.size();//finds the number of elements	
	double T = t[n2-1];//takes the last element for T (total time elapsed)
	SNR = sqrt(SNR/T);	
	
	long double pdh = (-2/T) * sum;//calculates p(d|h)
	return pdh;

}

vec_d ParameterFunction( double m1, double m2, double d, vec_d t ){ 
	
	int n = t.size();
	vec_d ht;

	parameters PARAMS;
	parameters *P = &PARAMS;

	// Set the characteristic parameters
	setFundamentalParameters(10.0, 
							 20.0,
							 0.0, 
							 10.0, 
							 m1, 
							 m2, 
							 d,
							 0.0,
							 0.0,
							 0.0, 
							 P);

	complex<double> gravWav = 0.0;
	complex<double> * gw = &gravWav;

	double df = P->df;
	double f;

	for( int i = 0; i < n; i++ ){

		f = double(i)*df;

		ht.push_back(updatedAmplitude(P,f,gw));
	}

	return ht;
}

vec_d NoiseFunction( vec_d f ){
	
	int n = f.size();	
	vec_d sf;

	for( int i = 0; i < n; i++ ){
		sf.push_back( 1e-39 ); //PLACEHOLDER: INSERT NOISE FUNCTION IN PARENTHESIS
	}

	return sf;
}

