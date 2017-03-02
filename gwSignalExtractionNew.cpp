#include <iostream>
#include <gsl/gsl_fft_complex.h>
#include "gwSignalExtraction.h"

#define REAL(z, i)  (z[i][0])
#define IMAG(z, i) (z[i][1])

Extractor::Extractor(){}

void Extractor::setSignalT(Signal* sig)
{
	mSignalT = sig;
	return;
}

void Extractor::setSignalF(Signal* sig)
{
	mSignalF = sig;
	return;
}

void Extractor::setTemplates(std::vector<Template>* temps)
{
	mTemplates = temps;
	return;
}

void Extractor::fft(std::vector<Template>* output)
{
	//Left as exercise for the reader
	return;
}

void Extractor::Convolution(std::vector<Template>* output)
{
	int I, J, pn;
	double result;
	
	I = mTemplates->size();
	J =mSignalF->waveform[0].size();

	pn = 1;
	
	for(int i=0; i<I; i++)
	{	
		Template* temp = &mTemplates[0][i];
		
		Template op;

		op.param[0] = temp->param[0];
		op.param[1] = temp->param[1];
		
		for(int j=0; j<J; j++)
		{
			op.waveform[0].push_back(mSignalF->waveform[0][j]);

			result = pn * mSignalF->waveform[1][j] * temp->waveform[1][j];

			op.waveform[1].push_back(result);

			pn = -pn;
		}

		output->push_back(op);
	}

	mConResults = output;
	
	return;
}

void Extractor::fftInverse(std::vector<Template>* output)
{
	int I, N;
	vec_d* amp; 
	vec_d* freq;
	double sw;

	I = mConResults->size();	
	N = (*mConResults)[0].waveform[0].size();	

	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N/2);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N/2);
	
	for(int i=0; i<I; i++)
	{		
		freq = &(*mConResults)[i].waveform[0];

		sw = N/(2*(*freq)[N-2]);
		
		amp = &(*mConResults)[i].waveform[1];

		gsl_fft_complex_inverse(&(*amp)[0], 1, N/2, compWT, complexWS);		

		Template temp;

		temp.param[0] = (*mConResults)[i].param[0];
		temp.param[1] = (*mConResults)[i].param[1];
 
		for(int n=0; n<N/2; n++)
		{
			temp.waveform[0].push_back(n/sw);
			temp.waveform[1].push_back((*amp)[2*n]);
		}

		output->push_back(temp);
	}

	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	

	return;
}
