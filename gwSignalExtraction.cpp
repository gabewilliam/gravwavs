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

void Extractor::fft(std::vector<Template>* output)
{
	//Left as exercise for the reader
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

