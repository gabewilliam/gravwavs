#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "gwSignalExtraction.h"

Extractor::Extractor(){}

void Extractor::setSignal(Signal* sig)
{
	mSignalT = sig;

	return;
}

void Extractor::setTemplates(std::vector<Template>* temps)
{
	mTemplatesT = temps;

	return;
}

void Extractor::fft(std::vector<Template>* tempsFFT)
{
	int I = mTemplatesT->size();

	for(int i=0; i<I; i++)
	{
		Template* temp = &mTemplatesT[0][i];
		Template tempFFT;

		vec_d time = temp->waveform[0];
		vec_d amp = temp->waveform[1];

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

		tempFFT.param[0] = temp->param[0];
		tempFFT.param[1] = temp->param[1];
		tempFFT.waveform[0] = freq;
		tempFFT.waveform[1] = amp;

		tempsFFT->push_back(tempFFT);
	}

	mTemplatesF = tempsFFT;
	
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

	mSignalF = sigFFT;

	return;
}

void Extractor::fftInverse(std::vector<Signal>* sigsFFTI)
{
	int I = sigsFFTI->size();

	for(int i=0; i<I; i++)
	{
		Signal* sig = &sigsFFTI[0][i];

		vec_d freq = sig->waveform[0];
		vec_d amp = sig->waveform[1];

		int N = freq.size();

		//Allocateing a work space and look up tables for the gsl fft function
		gsl_fft_real_workspace* realWS = gsl_fft_real_workspace_alloc(N);
		gsl_fft_halfcomplex_wavetable* halfcompWT = gsl_fft_halfcomplex_wavetable_alloc(N);

		gsl_fft_halfcomplex_inverse(&amp[0], 1, N, halfcompWT, realWS);

		gsl_fft_real_workspace_free(realWS);
		gsl_fft_halfcomplex_wavetable_free(halfcompWT);	

		vec_d time;

		///////////////////////////wrong
		for(int j=0; j<N; j++)
		{
			if(j<N/2)
				time.push_back(j);
			else
				time.push_back(j);
		} 
		///////////////////////////wrong

		sig->waveform[0] = time;
		sig->waveform[1] = amp;
	}

	return;
}


void Extractor::tConvolution(std::vector<Signal>* output)
{
	int I, J, K;
	double result;
	vec_d op;

	I = mTemplatesT->size();

	//Looping over all templates
	for(int i=0; i<I; i++)
	{
		Template* temp = &mTemplatesT[0][i];

		J = mSignalT->waveform[0].size();
		K = J;
		
		Signal op;

		for(int j=0; j<J; j++)
		{
			for(int k=0; k<K; k++)
			{
				result += mSignalT->waveform[1][j-k] * temp->waveform[1][k];
			}

			op.waveform[0].push_back(temp->waveform[0][j]);
			op.waveform[1].push_back(result);
		}

		

		output->push_back(op);
	}
	
	return;
}

void Extractor::fConvolution(std::vector<Signal>* output)
{
	int I, J, K;
	double result;
	vec_d op;

	I = mTemplatesF->size();

	

	//Looping over all templates
	for(int i=0; i<I; i++)
	{
		Template* temp = &mTemplatesF[0][i];

		J = mSignalF->waveform[0].size();
		K = J;
		
		Signal op;

		for(int k=0; k<K; k++)
		{
			if(k<K/2)
				result = mSignalF->waveform[1][k] * temp->waveform[1][k];
			else
				result = -(mSignalF->waveform[1][k] * temp->waveform[1][k]);

			op.waveform[0].push_back(temp->waveform[0][k]);
			op.waveform[1].push_back(result);
			
		}

		output->push_back(op);
	}
	
	return;
}
