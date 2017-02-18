#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "gwDataTypes.h"
#include "gwReadWrite.h"
#include "gwFFT.h"

//AIM: reads in d(t) from signal, generates h(t) from known formula -> fourier transforms both to get d(f) and h(f), s(f) is hopefully already known and in frequency domain.
//Uses d(f), h(f) and s(f) using equation A20 in "Veich, Vecchio (2010)" to calculate likelihood using two parameters m1 and m2.

vec_d ParameterFunction( double, double, vec_d );//This function is used to create model function ht to compare with the data
vec_d NoiseFunction( vec_d ); //This is used to create the noise probability as a function of frequency

double PdhFunction( double m1, double m2, std::string signalFile ){
	
	Signal *dt;//Data signal (time domain)
	std::vector< Signal >* idt; //loadSignals requires std::vector
	Signal *df;//Data signal (freq domain)
	Signal *ht;//Model signal (time domain)
	Signal *hf;//Model signal (freq domain)

	bool loadSignals(signalFile, idt, csv); //Loads signal into dt
	idt[0] = dt;//converts from vector of signals to signal

	vec_d ht2 = ParameterFunction( m1, m2, dt->waveform[0] );//Creates model function for given masses.

	ht->waveform[0] = dt->waveform[0]; //ht uses same time scale as dt
	ht->waveform[1] = ht2;	//sets ht equal to vector produced by function
	
	Extractor Fourier = Extractor(); //performs fourier transforms on each of the time domain signals 
	Fourier.setSignal(dt);
	Fourier.fft(df);
	Fourier.setSignal(ht);
	Fourier.fft(hf);

	vec_d sf = NoiseFunction( df->waveform[0] ); //creates noise probability function

	int n = df->waveform[0].size();//finds number of elements in freq domain signal

	double sum = 0.0;
	vec_d vdf = df->waveform[1]; //splits signal into std::vectors for use in loop
	vec_d vhf = hf->waveform[1];
	
	for( int i = 0; i < n; i++ ){ //evaluates the sum part in equation A20	
		double x = pow( std::abs( vdf[i] - vhf[i] ), 2 )/sf[i]; 		
		sum = sum + x;	
	}
	
	vec_d t = dt->waveform[0];//takes the time vector
	int n2 = t.size();//finds the number of elements	
	double T = t[n2];//takes the last element for T (total time elapsed)
	
	double pdh = exp( (-2/T) * sum );//calculates p(d|h)
	
	return pdh;
}

vec_d ParameterFunction( double m1, double m2, vec_d t ){ 
	
	int n = t.size();
	vec_d ht;

	for( int i = 0; i < n; i++ ){
		ht.push_back( m1*m2*t[i] ); //PLACEHOLDER: INSERT PARAMETER FUNCTION IN PARENTHESIS
	}

	return ht;
}

vec_d NoiseFunction( vec_d f ){
	
	int n = f.size();	
	vec_d sf;

	for( int i = 0; i < n; i++ ){
		sf.push_back( f[i] ); //PLACEHOLDER: INSERT NOISE FUNCTION IN PARENTHESIS
	}

	return sf;
}

