#include <gsl/gsl_fft_complex.h>
#include "gwSignalExtraction.h"

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

void Extractor::fft(std::vector<Template>* output){
	size_t I=output->size();
	int J = (mSignalT->waveform[0]).size();
	int M =(int)(log2(J)+2);
	size_t N= pow(2,M);
	
	double sampleFreq=J/((mSignalT->waveform[0])[J-1]-(mSignalT->waveform[0])[0]);
	
	vec_d freq;
	for(int j=0; j<N/4; j++){
		freq.push_back(j*sampleFreq/N);
		freq.push_back(j*sampleFreq/N);
	} 
	for(int j=N/4; j<N/2; j++){
		freq.push_back((N/4-j)*sampleFreq/N);
		freq.push_back((N/4-j)*sampleFreq/N);
	} 
	
	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* complexWT = gsl_fft_complex_wavetable_alloc(N);
	;
	for(size_t i=0;i<I;i++){
		
		vec_d* Amp=new vec_d;
	
		for(size_t j=0;j<J;j++){
			Amp->push_back((output[0][i].waveform[1])[j]);
			Amp->push_back(0);
			
			//(output[0][i].waveform[1]).push_back(0);
		}	
		(output[0][i].waveform[1])=*Amp;
		for(size_t j=J;j<N;j++){
			(output[0][i].waveform[1]).push_back(0);
			(output[0][i].waveform[1]).push_back(0);
		}
		
		gsl_fft_complex_forward(&(output[0][i].waveform[1])[0], 1, N, complexWT, complexWS);
		output[0][i].waveform[0]=freq;/**/
	}
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(complexWT);
	
	mTemplates=output;
	
	return;
	
}
void Extractor::fft(Signal* output){

	size_t J = (output->waveform[0]).size();
	int M =(int)(log2(J)+2);
	size_t N= pow(2,M);
	
	double sampleFreq=J/((output->waveform[0])[J-1]-(output->waveform[0])[0]);
	
	vec_d freq;
	for(int j=0; j<N/4; j++){
		freq.push_back(j*sampleFreq/N);
		freq.push_back(j*sampleFreq/N);
	} 
	for(int j=N/4; j<N/2; j++){
		freq.push_back((N/4-j)*sampleFreq/N);
		freq.push_back((N/4-j)*sampleFreq/N);
	} 
	
		vec_d* Amp=new vec_d;
	
		for(size_t j=0;j<J;j++){
			Amp->push_back((output->waveform[1])[j]);
			Amp->push_back(0);
			
			//(output[0][i].waveform[1]).push_back(0);
		}	
		(output->waveform[1])=*Amp;
		for(size_t j=J;j<N;j++){
			(output->waveform[1]).push_back(0);
			(output->waveform[1]).push_back(0);
		}
	

	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* complexWT = gsl_fft_complex_wavetable_alloc(N);

	gsl_fft_complex_forward (&(output->waveform[1])[0], 1, N, complexWT, complexWS);
	
	output->waveform[0]=freq;
	
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(complexWT);
	
	mSignalF=output;
	return;
}

void Extractor::Convolution(std::vector<Template>* output)
{
	int I, J, pn;
	double result;
	double noise;
	
	I = mTemplates->size();
	J = mSignalF->waveform[1].size();

	pn = 1;
	
	for(int i=0; i<I; i++)
	{	
		Template* temp = &mTemplates[0][i];
		
		Template op;

		op.param[0] = temp->param[0];
		op.param[1] = temp->param[1];
		
		for(int j=0; j<J; j++)
		{
			//noise = getASD(mSignalF->waveform[1][j]);

			op.waveform[0].push_back(mSignalF->waveform[0][j]);

			//std::cout << noise << std::endl;
		
			result = pn * mSignalF->waveform[1][j] * temp->waveform[1][j];// / (noise * noise);

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
	vec_d amp; 
	vec_d freq;
	double dt;

	I = (*mConResults).size();	
	N = (*mConResults)[0].waveform[0].size();	

	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N/2);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N/2);
	
	for(int i=0; i<I; i++)
	{		
		freq = (*mConResults)[i].waveform[0];

		dt = 1/(freq[N/2-1]);
	
		amp = (*mConResults)[i].waveform[1];

		gsl_fft_complex_inverse(&amp[0], 1, N/2, compWT, complexWS);		

		Template temp;

		temp.param[0] = (*mConResults)[i].param[0];
		temp.param[1] = (*mConResults)[i].param[1];
 
		for(int n=0; n<N/2; n++)
		{
			temp.waveform[0].push_back(n*dt);
			temp.waveform[1].push_back(amp[2*n]);
		}

		output->push_back(temp);
	}

	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	

	return;
}

bool Extractor::loadCurve(std::string filename)
{
	std::ifstream inFile;
	
	inFile.open(filename.c_str());

	double d;
	
	std::string line, element;

	NoiseCurve curve;

	while(getline(inFile,line))
	{
		std::istringstream iss1(line);
		
		while(getline(iss1, element, ','))
		{	
			std::istringstream is(element);
			is >> d;
			curve.freq.push_back(d);	
		}
		
		getline(inFile, line);
		
		std::istringstream iss2(line);
		
		while(getline(iss2, element, ','))
		{
			std::istringstream is(element);
			is >> d;
			curve.asd.push_back(d);		
		}
	}

	curve.fMin = curve.freq.front();
	curve.fMax = curve.freq.back();

	fNoiseCurve = curve;
	
	return true;
}

double Extractor::getASD(double f)
{
	
	if((f < fNoiseCurve.fMin) || (f > fNoiseCurve.fMax))
		return 0;
	
	double fTest=0.0, fPrev=0.0, asdTest=0.0, asdPrev=0.0, grad=0.0, asd=0.0;
	
	//Find best asd to use for given frequency
	for(int i = 0; i <= fNoiseCurve.freq.size(); i++)
	{
		fPrev = fTest;
		fTest = fNoiseCurve.freq[i];
		
		asdPrev = asdTest;
		asdTest = fNoiseCurve.asd[i];
		
		if(f==fTest)
		{
			asd = fNoiseCurve.asd[i];
			
			i = fNoiseCurve.freq.size() + 1;	
		}

		if(f<fTest)
		{			
			grad = (asdTest-asdPrev) / (fTest-fPrev);
			
			asd = asdPrev + ((f-fPrev) * grad);
			
			i = fNoiseCurve.freq.size() + 1;	
		}
	}

	return asd;
}

double Extractor::fAutoCorrComplex(Template* temp)
{
	int N, pn; 
	double result;
	vec_d op;

	N = temp->waveform[0].size();
	pn = 1;

	for(int n=0; n<N; n++)
	{  
		result = pn * temp->waveform[1][n] * temp->waveform[1][n];	
		op.push_back(result);

		pn = -pn;
	}

	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N);

	gsl_fft_complex_inverse(&op[0], 1, N, compWT, complexWS);

	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	
	
	return op[0];
}

