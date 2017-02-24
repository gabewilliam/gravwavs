#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <iostream>

#include "gwSignalExtraction.h"
#define REAL(z, i)  ((z[2*(i)]))
#define IMAG(z, i) ((z[2*(i)+1]))//redundant

Extractor::Extractor(){}

void Extractor::setSignal(Signal* sig){
	mSignalT = sig;
	mOriginalTime=sig->waveform[0];
	return;
}
void Extractor::setTemplates(std::vector<Template>* temps){
	mTemplatesT = temps;

	return;
}


void Extractor::fft(std::vector<Template>* tempsFFT){
	int I = mTemplatesT->size();

	for(int k=0; k<I; k++){
		Template* temp = &mTemplatesT[0][k];
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
void Extractor::fft(Signal* sigFFT){
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
void Extractor::fftComplex(std::vector<Template>* tempsFFT){
	//Help I am been forced to comment against my will
	
	//Calculate number of Templates
	int I = tempsFFT->size();
	
	//Used for converting to frequency domain 
	Template* tempa = &tempsFFT[0][0];
	vec_d time = tempa->waveform[0];
	
	//calculate the length of the signal
	size_t L = (mSignalT->waveform[1]).size();
	//std::cout<<L<<std::endl;
	
	//calculate the next highest power of two and zero pad to double 
	//increasing to the next power of two increases run speed of fft algorithm
	//increasing to the one beyond that prevents reflections from interfering with each other	
	int M =(int)(log2(L)+2);
	size_t N= pow(2,M);
	
	
	//this isn't working or essential at the moment; can come back and fix if necessary
	vec_d freq;
	for(size_t j=0; j<N/2; j++){
		freq.push_back(j);
		freq.push_back(N-j);
	} 
	////////////////////////////////

	//Allocateing a work space and look up tables for the gsl fft function
	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* complexWT = gsl_fft_complex_wavetable_alloc(N);
	
	for(int k=0; k<I; k++){
		
		//cycle through the templates	
		//std::cout<<"FFt1  "<<k<<std::endl;
		Template* temp = &tempsFFT[0][k];
		Template* tempFFT=new Template;

		
		std::cout<<"FFt  "<<k<<std::endl;
		
		
		//load in the template		
		vec_d amp = temp->waveform[1];

		size_t K = amp.size();
		
		//Create vector for transform on the heap.	
		vec_d *filterData=new vec_d;

		

		//input data
		for (size_t i = 0; i < N; i++){
			//fill "real part" of the vector with template. if template shorter than signal, fills with zeroes so they are the same length
			if(i<K)
				filterData->push_back(amp[i]);
			else
				filterData->push_back(0);
			
			//put an 0 in for the "imaginary" part
			filterData->push_back(0);
		}
		


		gsl_fft_complex_forward(&(*filterData)[0], 1, N, complexWT, complexWS);
		


		//fill up template with data	tempFFT->param[0] = temp->param[0];
		tempFFT->param[1] = temp->param[1];
		tempFFT->waveform[0] = freq;
		tempFFT->waveform[1]=*filterData;

		//add template to the vector		
		tempsFFT[0][k]=*tempFFT;

		
		
	}
	//free work space
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(complexWT);
	
	//save as the data member for use in other functions	
	mTemplatesF = tempsFFT;

	/**/	return;
	
	
}
void Extractor::fftComplex(Signal* sigFFT){
	
	
	//get the signal from the data member
	vec_d time = mSignalT->waveform[0];
	vec_d amp = mSignalT->waveform[1];

	//
	size_t L=amp.size();
	int M =(int)(log2(L)+2);
	size_t N= pow(2,M);
	size_t P = time.size();
		//calculate length of signal
	
	vec_d signalData;
	//decalre another holder
	
	//initialise to 0
		for (size_t i = 0; i < L; i++){//fill up vector
		
			signalData.push_back(amp[i]);
			//for real parts put in the signal
			signalData.push_back(0);
			//for imaginary parts put in 0
		}
		
		for (size_t j=L;j<N;j++){
			signalData.push_back(0);
			signalData.push_back(0);
		}
		
	//Allocateing a work space and look up tables for the gsl fft function
	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* complexWT = gsl_fft_complex_wavetable_alloc(N);

	gsl_fft_complex_forward (&signalData[0], 1, N, complexWT, complexWS);
	//perform transform
	
	
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(complexWT);
	//free workspace
	
	//this isn't working or essential at the moment; can come back and fix if necessary
	double sw = N/(2*time[N-1]);	

	vec_d freq;


	for(size_t j=0; j<N/2; j++){
		freq.push_back(j);
		freq.push_back(N-j);	
	}
	
	//////////////////////////////
	sigFFT->waveform[0] = freq;
	sigFFT->waveform[1]=signalData;
	//load tranform into data member
	mSignalF = sigFFT;


	return;

}

void Extractor::fftInverse(std::vector<Signal>* sigsFFTI){
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

void Extractor::fftInverseComplex(std::vector<Signal>* convsFFTI){
	int I = convsFFTI->size();
	//calculate number of filters	

	vec_d time,timedouubled;
	//set the time to that of the original signal
	time = mSignalT->waveform[0];	
	size_t P=mOriginalTime.size();
	
	
	size_t N=convsFFTI[0][0].waveform[1].size();
	
	
	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N/2);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N/2);
	//Allocateing a work space and look up tables for the gsl fft function
	
	for(int k=0; k<I; k++){
		// run transform for each signal
		//set the current signal
		Signal* conv = &convsFFTI[0][k];
		//std::cout<<"FFtInv  "<<k<<std::endl;
		
		//set the vector to be current convolution
		vec_d convData= conv->waveform[1];
		std::cout<<"FFtInv  "<<k<<std::endl;
		size_t N=convData.size();

		gsl_fft_complex_inverse(&convData[0], 1, N/2, compWT, complexWS);		
		
		vec_d convTime;
		//declare a vector to hold the magntiude of the convolution in time domain 
		for(int l=0;l<P;l++){
			convTime.push_back(convData[2*l]);//+convData[l+1]);
			//l++;
		}//i think this is right it might need correcting
		
		
		(*convsFFTI)[k].waveform[0] = mOriginalTime;
		(*convsFFTI)[k].waveform[1]=convTime;//set data to be back into time domain
		std::cout<<P<<std::endl;
	}
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	
	//free the work space
	return;
}

void Extractor::tConvolution(std::vector<Signal>* output)
{
	int I, J, K;
	double result;

	I = mTemplatesT->size();

	//Looping over all templates
	for(int i=0; i<I; i++)
	{
		Template* temp = &mTemplatesT[0][i];

		J = mSignalT->waveform[0].size();
		K = temp->waveform[0].size();
		
		Signal op;

		for(int j=0; j<J; j++)
		{
			result = 0;
			
			for(int k=0; k < K && j+k < J; k++)
			{
				result += mSignalT->waveform[1][j+k] * temp->waveform[1][k];
			}

			op.waveform[0].push_back(mSignalT->waveform[0][j]);
			op.waveform[1].push_back(result);
		}

		

		output->push_back(op);
	}
	
	return;
}

void Extractor::fConvolution(std::vector<Signal>* output){
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
			if(k<K/2){
				result = mSignalF->waveform[1][k] * temp->waveform[1][k];
			}
			
			else{
				result = -(mSignalF->waveform[1][k] * temp->waveform[1][k]);
			}
			
			op.waveform[0].push_back(temp->waveform[0][k]);
			op.waveform[1].push_back(result);
			
		}

		output->push_back(op);
	}
	
	return;
}
void Extractor::fConvolutionComplex(std::vector<Signal>* output){
	int I, K;//declare size holders
	double result;//declare something to hold tempory results
	vec_d op;//decalre a vector to hold convolution for each run
	
	I = mTemplatesF->size();//calculate number of templates 
	
	K =mSignalF->waveform[0].size();//declare number of data points
	//std::cout<<"Conv  "<<I<<std::endl;

	//Looping over all templates
	for(int i=0; i<I; i++){
		
		Template* temp = &mTemplatesF[0][i];
		//decalre a template to hold the current filters data
		std::cout<<"conv  "<<i<<std::endl;
		Signal op;
		//redeclare the signal vector
		for(int k=0; k<K; k++){
			//calculate and push back the data with a waveform  
			result = mSignalF->waveform[1][2*k] * temp->waveform[1][2*k];	
			op.waveform[0].push_back(2*k);
			op.waveform[1].push_back(result);
			//and its conjugate if it is imaginary
			result = -(mSignalF->waveform[1][2*k+1] * temp->waveform[1][2*k+1]);
			op.waveform[0].push_back(2*k);
			op.waveform[1].push_back(result);
		
		}
		//put this into the output vector
		output->push_back(op);
	}
	
	return;
}
