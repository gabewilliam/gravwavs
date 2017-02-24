#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <iostream>

#include "gwSignalExtraction.h"

#define REAL(z, i)  (z[i][0])
#define IMAG(z, i) (z[i][1])//not redundant

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
void Extractor::fftComplex(Signal* sigFFT)
{
	//get the signal from the data member
	vec_d time = mSignalT->waveform[0];
	vec_d amp = mSignalT->waveform[1];

	size_t L=amp.size();
	int M =(int)(log2(L)+2);
	size_t N= pow(2,M);
	size_t P = time.size();
	
	vec_d signalData;

	for (size_t i=0; i<L; i++)
	{
		signalData.push_back(amp[i]);
		signalData.push_back(0);
	}
		
	for (size_t j=L; j<N; j++)
	{
		signalData.push_back(0);
		signalData.push_back(0);
	}

	//Fourier Transform
	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* complexWT = gsl_fft_complex_wavetable_alloc(N);

	gsl_fft_complex_forward (&signalData[0], 1, N, complexWT, complexWS);
	
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(complexWT);

	
	//Frequency Scaling
	double sw = N/(2*time[N-1]);	
	vec_d freq;

	for(size_t j=0; j<N/2; j++)
	{
		freq.push_back(j*sw);
		freq.push_back(-j*sw));	
	}

	sigFFT->waveform[0] = freq;
	sigFFT->waveform[1]=signalData;
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

void Extractor::gateFilter(Signal* data, double freq)
{
	for(int i=0; data->waveform[0][2*i] < freq; i++)
		data->waveform[1][2*i] = 0;
	
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

void Extractor::fConvolutionComplex(std::vector<Signal>* output)
{
	int I, K;//declare size holders
	double result, freq;//declare something to hold tempory results
	vec_d op;//decalre a vector to hold convolution for each run
	
	I = mTemplatesF->size();//calculate number of templates 
	
	K =mSignalF->waveform[0].size();//declare number of data points
	
	//Looping over all templates
	for(int i=0; i<I; i++)
	{	
		Template* temp = &mTemplatesF[0][i];
		
		Signal op;
		
		for(int k=0; k<K; k++)
		{
			//calculate and push back the data with a waveform  
			result = mSignalF->waveform[1][2*k] * temp->waveform[1][2*k];	
			op.waveform[0].push_back(2*k);
			op.waveform[1].push_back(result);
			
			//and its conjugate
			result = -(mSignalF->waveform[1][2*k+1] * temp->waveform[1][2*k+1]);
			op.waveform[0].push_back(2*k);
			op.waveform[1].push_back(result);
		}
		
		gateFilter(op, 20);
		
		//put this into the output vector
		output->push_back(op);
	}
	
	return;
}

//FFTW FUNCTIONS__________________________________

void Extractor::fftw(std::vector<Template>* templatesFFT){
	//number of templates
	int T = mTemplatesT->size();
	//size of templates
	int N;
	
	int size = mSignalT->waveform[0].size();
	int allocSize = 2*size - 1;
	int j = 0;
	double sw;
	

	//used to input into mTemplatesF
	vec_d freq;
	vec_d result;
	
	//templates for ease
	Template *tempT;
	Template tempF;

	fftw_complex *filterT;
	fftw_complex *filterF;
	
	filterT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	filterF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	
	fftw_plan forwardPlan = fftw_plan_dft_1d(allocSize, filterT, filterF, FFTW_FORWARD, FFTW_ESTIMATE);
	
	//loop through templates
	for (int i = 0; i < T; i++){
	
		//sets the current template	
		tempT = &mTemplatesT[0][i];
		//gets the size of the template
		N = tempT->waveform[1].size();
		
		//loops up the the size of the template and then fills the reast with 0 zeros up to 2 * signal size
		for (j = 0; j < N; j++){
			REAL(filterT, j) = tempT->waveform[1][j];
			IMAG(filterT, j) = 0.0;		
		}
		for (j = N; j < allocSize; j++){
			REAL(filterT, j) = 0.0;
			IMAG(filterT, j) = 0.0;			
		}
		
		fftw_execute(forwardPlan);
		
		//the frequency stuff is mostly wrong
		sw = N/(2.0*(tempT->waveform[0])[N-1]);
		
		for(j = 0; j < size; j++){
			if (j < size/2){
				freq.push_back(sw*j);
			}else{
				freq.push_back(-sw*(j - size/2));
			}		
		}
		//the size of this might have to be adjusted
		//should save the result as a vec_d packed as in the gsl complex routines		
		for (j = 0; j < allocSize; j++){
			result.push_back(REAL(filterF, j));
			result.push_back(IMAG(filterF, j));
		}
		
		//allocate each 
		tempF.param[0] = tempT->param[0];
		tempF.param[1] = tempT->param[1];
		tempF.waveform[0] = freq;
		tempF.waveform[1] = result;
		
		templatesFFT->push_back(tempF);
			
	}
	mTemplatesF = templatesFFT;
	
	fftw_free(filterT);
	fftw_free(filterF);
	fftw_destroy_plan(forwardPlan);

} 

void Extractor::fftw(Signal* signalFFT){

	vec_d freq;
	vec_d result;
	
	//the FFTW stuff
	fftw_complex *signalT;
	fftw_complex *signalF;
	
	int N = mSignalT->waveform[0].size();
	int allocSize = 2*N - 1;
	
	signalT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	signalF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	
	fftw_plan forwardPlan = fftw_plan_dft_1d(allocSize, signalT, signalF, FFTW_FORWARD, FFTW_ESTIMATE);
	
	//loops to the 2 * signal size -1
	//zero pads for the later half
	for (int i = 0; i < N; i++){
		REAL(signalT, i) = mSignalT->waveform[0][i];
		IMAG(signalT, i) = 0.0;
	}
	for (int i = N; i < allocSize; i++){
		REAL(signalT, i) = 0.0;
		IMAG(signalT, i) = 0.0;
	}
	
	fftw_execute(forwardPlan);
	
	//pack the vectors as was done for the gsl complex routine
	for (int i = 0; i < allocSize; i++){
		result.push_back(REAL(signalF, i));
		result.push_back(IMAG(signalF, i));
	}
	
	signalFFT->waveform[1] = result;

	fftw_free(signalT);
	fftw_free(signalF);
	fftw_destroy_plan(forwardPlan);

}


//THIS FUNCTION IS NOT SCALED SO THE Y AXIS WILL NOT BE CORRECT
void Extractor::fftwInverse(std::vector<Signal>* signalsFFTI){
	
	//number of templates
	int N = signalsFFTI->size();	
	
	//gets the time of the initial signal
	vec_d time;
	time = mSignalT->waveform[0];
	vec_d result;
	
	int j = 0;
	
	//gets the size of each of the signals
	int allocSize = (*signalsFFTI)[0].waveform[1].size()/2;
	
	//FFTW stuff
	fftw_complex *convoT;
	fftw_complex *convoF;	
	
	convoT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	convoF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	
	fftw_plan backwardPlan = fftw_plan_dft_1d(allocSize, convoF, convoT, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	for (int i = 0; i < N; i++){
		//unpacks the complex vectors from the initial routines
		for (j = 0; j < allocSize; j++){
			REAL(convoF, j) = (*signalsFFTI)[i].waveform[1][2*j];
			IMAG(convoF, j) = (*signalsFFTI)[i].waveform[1][2*j + 1];
		}
		fftw_execute(backwardPlan);
		
		for (j = 0; j < allocSize; j++){
			result.push_back(REAL(convoT, j));
		}
		//change the input to the time domain
		(*signalsFFTI)[i].waveform[0] = time;
		(*signalsFFTI)[i].waveform[1] = result;

	}
	

	fftw_free(convoT);
	fftw_free(convoF);
	fftw_destroy_plan(backwardPlan);

}

