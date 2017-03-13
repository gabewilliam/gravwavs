#include <gsl/gsl_fft_complex.h>
#include "gwSignalExtraction.h"

#include <cmath>
#include <iostream>

Extractor::Extractor(){}

void Extractor::setSignalT(Signal* sig){
	mSignalT = sig;
	return;
}

void Extractor::setSignalF(Signal* sig){
	mSignalF = sig;
	return;
}

void Extractor::setTemplates(std::vector<Template>* temps){
	mTemplates = temps;
	return;
}

void Extractor::splitSignal(){

	int N = mSignalT->waveform[0].size();
	double base = log2(N);
	std::vector<Signal*> *dummy = new std::vector<Signal*>;
	
	int size = ceil(base);
	
	
	size = pow(2,size);
	
	double dt = mSignalT->waveform[0][2] - mSignalT->waveform[0][0];
	std::cout<<"dt = "<<dt<<std::endl;
	
	for (int i = N; i < size; i++){
		mSignalT->waveform[1].push_back(0.0);
	}
	mSignalT->waveform[0].clear();
	
	for(int i = 0; i < size/2; i++){
		mSignalT->waveform[0].push_back(dt*i);
		mSignalT->waveform[0].push_back(dt*i);
	}
	
	std::cout<<mSignalT->waveform[0].size()<<"  ,  "<<mSignalT->waveform[1].size()<<std::endl;
	
	int chunkSize = (*mTemplates)[0].waveform[0].size();
	
	int nChunks = 2*size/chunkSize;
	
	int start;
	
	std::cout<<nChunks<<"  ,  "<<log2(chunkSize)<<std::endl;		
		
	for (int i = 0; i < nChunks; i++){
		
		Signal *sig = new Signal;
		start = i*(chunkSize)/2;
		if (i < nChunks){
			for(int j = start; j < start + chunkSize; j++){			
				
				sig->waveform[0].push_back(mSignalT->waveform[0][j]);
				sig->waveform[1].push_back(mSignalT->waveform[1][j]);			
					
			}	
		}else{
			for(int j = start; j < start + chunkSize/2; j++){
				sig->waveform[0].push_back(mSignalT->waveform[0][j]);
				sig->waveform[1].push_back(mSignalT->waveform[1][j]);
			}
			for(int j = 0; j < chunkSize/2; j++){
				sig->waveform[0].push_back(mSignalT->waveform[0][j]);
				sig->waveform[1].push_back(mSignalT->waveform[1][j]);
			}
		}
		dummy->push_back(sig);
	}	
	
	mSignalsT = dummy;

	
}

void Extractor::fft(){
	
	
	std::cout<<mSignalsT->size()<<std::endl;
	
	int N = (*mSignalsT)[0]->waveform[0].size()/2;
	int I = (*mSignalsT).size();


	
	std::vector<Signal*>* dummy = new std::vector<Signal*>;
	
	double dt = ((*mSignalsT)[0]->waveform[0][2] - (*mSignalsT)[0]->waveform[0][0]);
	
	double sampleFreq = 1/(2.0*N*dt);
	
	std::cout<<"dt = "<<dt<<std::endl;
	std::cout<<"df = "<<sampleFreq<<std::endl;
	
	std::cout<<"max freq = "<<sampleFreq*N<<std::endl;
	
	gsl_fft_complex_workspace* complexWS =  gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* complexWT = gsl_fft_complex_wavetable_alloc(N);
	
	vec_d freq;
	for(int j=0; j<N/2; j++){
		freq.push_back(j*sampleFreq);
		freq.push_back(j*sampleFreq);
	} 
	for(int j=N/2; j<N; j++){
		freq.push_back((N/2-j)*sampleFreq);
		freq.push_back((N/2-j)*sampleFreq);
	} 
	
	vec_d amp;
		
	for (int i = 0; i < I; i++){
				
		Signal *sig = new Signal;
		
		amp = (*mSignalsT)[i]->waveform[1];
	
		gsl_fft_complex_forward (&amp[0], 1, N, complexWT, complexWS);
		
		
		sig->waveform[0] = freq;
		
		sig->waveform[1] = amp;
		dummy->push_back(sig);
	
	}
		
	mSignalsF = dummy;
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(complexWT);
	
	
	return;
}
void Extractor::Convolution/*crossCorrelation*/()
{
	int I, J;
	double imagResult;
	double realResult;
	double noise;
	double noisePower;
	double realTemp = 0;
	double imagTemp = 0;
	double df = fabs((*mSignalsF)[0]->waveform[0][2] - (*mSignalsF)[0]->waveform[0][0]);
	double norm;
	double freq;
	
	
	
	I = mTemplates->size();
	J = (*mSignalsF)[0]->waveform[1].size();
	
	
	
	std::vector<std::vector<Template*> *>* result = new std::vector<std::vector<Template*>* >;
	
	int signals = mSignalsF->size();
	
	for (int i = 0; i < I; i++){
		Template temp = (*mTemplates)[i];
		std::cout<<"noise number: "<<i<<std::endl;
		for(int j = 0; j < J/2; j++){
			
			freq = temp.waveform[0][2*j];
			if(freq < 9){
				noise = 1.0;
			}else{
				noise = getASD(freq);
			}
		
			noisePower = noise*noise;
			//std::cout<<noisePower<<"  ,  "<<freq<<std::endl;
			
			realTemp += temp.waveform[1][2*j] * temp.waveform[1][2*j] / noisePower;
			imagTemp += temp.waveform[1][2*j] * temp.waveform[1][2*j] / noisePower;
		
		}
		
		
		realTemp = realTemp * 2* df;
		imagTemp = imagTemp * 2* df;
		
		norm = realTemp * realTemp + imagTemp * imagTemp;
		
		norm = sqrt(norm);
		norm = sqrt(norm);
		
		mNorms.push_back(1/norm);
		
		
	}	
	
	std::cout<<"noise here guys"<<std::endl;
	for (int k = 0; k < signals; k++){
		
		std::vector<Template*>* dummy = new std::vector<Template*>;
		
		std::cout<<"signal: "<<k<<std::endl;
				
		for(int i=0; i<I; i++)
		{	
			Template temp = mTemplates[0][i];
		
			Template *op = new Template;

			op->param[0] = temp.param[0];
			op->param[1] = temp.param[1];
		
			realTemp = 0;
			imagTemp = 0;
		
			//not sure if this is looping through all of the elements correctly;
			//the alternative for the more sensible frequency packing is commented out
			for(int j = 0 ; j < J/2; j++){
			
				/*
				
				noise = getASD(mSignalF->waveform[0][j]);
				op.waveform[0].push_back(mSignalF->waveform[0][j]);
			
				noisePower = noise*noise;
				*/	
				freq = fabs((*mSignalsF)[k]->waveform[0][2*j]);
				
				
				
				
				noise = getASD((*mSignalsF)[k]->waveform[0][2*j]);
				
				if (noise == 0){
					noise = 1;
				}
				

				op->waveform[0].push_back((*mSignalsF)[k]->waveform[0][2*j]);
				op->waveform[0].push_back((*mSignalsF)[k]->waveform[0][2*j+1]);
				
				
				noisePower = 1/(noise*noise);
				realResult = (*mSignalsF)[k]->waveform[1][2*j] * temp.waveform[1][2*j] * noisePower;
				imagResult = -(*mSignalsF)[k]->waveform[1][2*j+1] * temp.waveform[1][2*j+1] * noisePower;
			
				op->waveform[1].push_back(realResult); //this should maybe be changed to SNR
				op->waveform[1].push_back(imagResult);
			}	
		
		

			dummy->push_back(op);
		}
		result->push_back(dummy);
	}
	mConResults = result;

	
	
	
	
	return;
}



void Extractor::fftInverse(std::vector<Template>* out1, std::vector<Template>* out2)
{
	int I, N;
	vec_d amp; 
	//vec_d freq;
	double dt;
	double norm;
	int index;
	int fullSignal = mSignalT->waveform[1].size();
	std::cout<<fullSignal<<"  ,  "<<mSignalT->waveform[0].size()<<std::endl;
	int signals = (*mSignalsF).size();
	vec_d *time1 = &((*mSignalT).waveform[0]);
	vec_d time2;
	
	//set times
	int limit = (*time1).size()/signals;
	for(int i = limit; i < (*time1).size(); i++){
		time2.push_back((*time1)[i]);
	} 
	for(int i = 0; i < limit; i++){
		time2.push_back((*time1)[i]);
	}
	

	I = (*mConResults)[0]->size();	
	N = (*(*mConResults)[0])[0]->waveform[0].size();	

	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N/2);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N/2);
	
	
	std::vector<std::vector<Template*>* >* final = new std::vector<std::vector<Template*>* >;
	

	
	std::cout<<mNorms[0]<<std::endl;
	
	
	
	std::cout<<(*mSignalsT)[0]->waveform[0].size()<<"  ,  "<<N<<std::endl;
	
	for(int k = 0; k < signals; k++){
	
		std::vector<Template*> *result = new std::vector<Template*>;
		for(int i=0; i<I; i++)
		{		
			Template *temp = new Template;
			//freq = (*mConResults)[i].waveform[0];
			//dt = 2.0/(freq[N/2-1]);
			amp = (*(*mConResults)[k])[i]->waveform[1];
			gsl_fft_complex_inverse(&amp[0], 1, N/2, compWT, complexWS);		
			norm = mNorms[i];
			std::cout<<norm<<std::endl;
			temp->param[0] = (*(*mConResults)[k])[i]->param[0];
			temp->param[1] = (*(*mConResults)[k])[i]->param[1];
			
			for(int n=0; n<N/2; n++)
			{
				if (2*n > (*mSignalsT)[k]->waveform[0].size()){
					std::cout<<"tooo biiiiiig"<<std::endl;
				}
				temp->waveform[0].push_back((*mSignalsT)[k]->waveform[0][2*n]);
				temp->waveform[1].push_back(fabs(amp[2*n] / norm));
			}
			
			
			result->push_back(temp);
			
		}
		
		final->push_back(result);		
		
	}

	clearSignalMemory();
	clearConvolutionMemory();
	
	mToBeDeleted = final;
	
	
	std::cout<<"hello"<<std::endl;
	
	for (int i = 0; i < I; i++){
		std::cout<<i<<std::endl;
		Template temp1;
		Template temp2;
		temp1.param[0] = (*mTemplates)[i].param[0];
		temp1.param[1] = (*mTemplates)[i].param[1];
		for (int j = 0; j < signals/2; j++){
		std::cout<<j<<std::endl;
			for(int k = 0; k < N/2; k++){
				
				
				temp1.waveform[0].push_back((*((*final)[2*j]))[i]->waveform[0][k]);
				temp1.waveform[1].push_back((*((*final)[2*j]))[i]->waveform[1][k]);
				temp2.waveform[0].push_back((*((*final)[2*j]))[i]->waveform[0][k]);
				temp2.waveform[1].push_back((*((*final)[2*j]))[i]->waveform[1][k]);
			}
		}
		
		
		out1->push_back(temp1);
		out2->push_back(temp2);

	}
	


	//out2 = outs2;

	std::cout<<"asdfvrfusbecd"<<std::endl;
	std::cout<<(*out1)[0].waveform[0][0]<<std::endl;
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	
	
	
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
double Extractor::getASD(double f){
	
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

double Extractor::fAutoCorrComplex(Template* temp){
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

	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N/2);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N/2);

	gsl_fft_complex_inverse(&op[0], 1, N/2, compWT, complexWS);

	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	

	return op[0];
}


void Extractor::clearSignalMemory(){

	int signals = mSignalsT->size();
	for(int i = signals-1; i >= 0; i--){
		delete ((*mSignalsT)[i]);
		delete ((*mSignalsF)[i]);
		mSignalsT->erase(mSignalsT->begin() + i);
		mSignalsF->erase(mSignalsF->begin() + i);
	}
	delete (mSignalsT);
	delete (mSignalsF);


}
void Extractor::clearConvolutionMemory(){

	int signals = mConResults->size();
	int templates = (*mConResults)[0]->size();
	
	for(int i = signals-1; i >= 0; i--){
		for(int j = templates-1; j >= 0; j--){
			delete ((*(*mConResults)[i])[j]);
		}
		delete ((*mConResults)[i]);
	}
	delete (mConResults);
}

void Extractor::clearOtherMemory(){
	int signals=mToBeDeleted->size();
	int templates=(*mToBeDeleted)[0]->size();

	for(int i = signals-1; i >= 0; i--){
		for(int j = templates-1; j >= 0; j--){
			delete ((*(*mToBeDeleted)[i])[j]);
		}
		delete ((*mToBeDeleted)[i]);
	}
	delete(mToBeDeleted);
}

void Extractor::clear(){
	clearSignalMemory();
	clearConvolutionMemory();
	clearOtherMemory();
}

