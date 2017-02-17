#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <iostream>

#include "gwFFT.h"

Extractor::Extractor(){}

void Extractor::setSignal(Signal* sig)
{
	mSignalT = sig;

	return;
}

void Extractor::fft(Signal* sigFFT)
{
	vec_d time = mSignalT->waveform[0];
	vec_d amp = mSignalT->waveform[1];

	int N = time.size();

	//Allocateing a work space and look up tables for the gsl fft function
	gsl_fft_real_workspace* realWS = gsl_fft_real_workspace_alloc(N);
	gsl_fft_real_wavetable* realWT = gsl_fft_real_wavetable_alloc(N);

	gsl_fft_real_transform(&amp[0], 1, N, realWT, realWS);

	gsl_fft_real_workspace_free(realWS);
	gsl_fft_real_wavetable_free(realWT);

	//spectral width is 1/dt where dt is the time spacing between data points
	double sw = N/(2*time[N-1]);	

	vec_d freq;

	for(int j=0; j<N; j++)
	{
		if(j<N/2)
			freq.push_back(sw*j);
		else
			freq.push_back(-sw*(j-N/2));
	} 

	sigFFT->waveform[0] = freq;
	sigFFT->waveform[1] = amp;

	return;
}
