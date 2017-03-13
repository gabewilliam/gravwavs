#include <gsl/gsl_fft_complex.h>
#include "gwSignalExtraction.h"
#include <cassert>

#include <cmath>
#include <iostream>

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
<<<<<<< HEAD

void Extractor::splitSignal(){

	int N = mSignalT->waveform[0].size();
	double base = log2(N);

	std::vector<Signal*> *dummy = new std::vector<Signal*>;
	int size = ceil(base);
	
	size = pow(2,size);
	
	double dt = mSignalT->waveform[0][2] - mSignalT->waveform[0][0];
	
	for (int i = N; i < size; i++){
		mSignalT->waveform[1].push_back(0.0);
	}
	mSignalT->waveform[0].clear();
	
	for(int i = 0; i < size/2; i++){
		mSignalT->waveform[0].push_back(dt*i);
		mSignalT->waveform[0].push_back(dt*i);
	}

	int chunkSize = (*mTemplates)[0].waveform[0].size();
	
	int nChunks = size/chunkSize;
	
	int start;
	
	std::cout<<chunkSize<<"  ,  "<<size<<"  ,  "<<mSignalT->waveform[0].size()<<"  ,  "<<mSignalT->waveform[1].size()<<std::endl;
		
	for (int i = 0; i < nChunks; i++){
		
		Signal *sig = new Signal;
		start = i*(chunkSize);
		if (i < nChunks){
			for(int j = start; j < start + chunkSize; j++){			
				
				sig->waveform[0].push_back(mSignalT->waveform[0][j]);
				sig->waveform[1].push_back(mSignalT->waveform[1][j]);			
					
			}	
		}
		dummy->push_back(sig);
	}	
	
	mSignalsT = dummy;

||||||| merged common ancestors
void Extractor::setTemplatesT(std::vector<Template>* temps)
{
	mTemplatesT = temps;
	return;
=======

void Extractor::splitSignal(){

	int N = mSignalT->waveform[0].size();
	double base = log2(N);
	int size = base;
	std::vector<Signal*> *dummy = new std::vector<Signal*>;
	
	size = ceil(base);
	
	
	size = pow(2,size);
	
	double dt = mSignalT->waveform[1][1] - mSignalT->waveform[1][0];
	
	for (int i = N; i < size; i++){
		mSignalT->waveform[1].push_back(0.0);
	}
	mSignalT->waveform[0].clear();
	
	for(int i = 0; i < size; i++){
		mSignalT->waveform[0].push_back(dt*i);
		mSignalT->waveform[0].push_back(dt*i);
	}
	
	int chunkSize = (*mTemplates)[0].waveform[0].size();
	
	int nChunks = size/chunkSize;
	
	int start;
	
	std::cout<<chunkSize<<"  ,  "<<size<<"  ,  "<<mSignalT->waveform[0].size()<<"  ,  "<<mSignalT->waveform[1].size()<<std::endl;
		
	for (int i = 0; i < nChunks; i++){
		
		Signal *sig = new Signal;
		start = i*(chunkSize);
		if (i < nChunks){
			for(int j = start; j < start + chunkSize; j++){			
				
				sig->waveform[0].push_back(mSignalT->waveform[0][j]);
				sig->waveform[1].push_back(mSignalT->waveform[1][j]);			
					
			}	
		}
		dummy->push_back(sig);
	}	
	
	mSignalsT = dummy;

>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
}
/*
void Extractor::fft(std::vector<Template>* output){
	size_t I=mTemplatesT->size();
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

	
	for(size_t i=0;i<I;i++){
		
		vec_d* Amp=new vec_d;
		Template* Temp=new Template;
		J=mTemplatesT[0][i].waveform[1].size();
		
		for(size_t j=0;j<J;j++){
			Amp->push_back((mTemplatesT[0][i].waveform[1])[j]);
			Amp->push_back(0);
		}	
		//(output[0][i].waveform[1])=*Amp;
		for(size_t j=J;j<N;j++){
			Amp->push_back(0);
			Amp->push_back(0);
		}
		
		gsl_fft_complex_forward(&(*Amp)[0], 1, N, complexWT, complexWS);
		
		Temp->param[0]=mTemplatesT[0][i].param[0];
			Temp->param[1]=mTemplatesT[0][i].param[1];
		Temp->waveform[0]=freq;/*
		Temp->waveform[1]=*Amp;
		output->push_back(*Temp);
		
	}
	setTemplates(output);
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(complexWT);
	
	return;
}*/
void Extractor::fft(){
	
	
	std::cout<<mSignalsT->size()<<std::endl;
	
	int N = (*mSignalsT)[0]->waveform[0].size()/2;
	int I = (*mSignalsT).size();


	
<<<<<<< HEAD
	std::vector<Signal*>* dummy = new std::vector<Signal*>;
	
	double sampleFreq = 1/(2.0*N*((*mSignalsT)[0]->waveform[0][2] - (*mSignalsT)[0]->waveform[0][0]));

	
	
	
	gsl_fft_complex_workspace* complexWS =  gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* complexWT = gsl_fft_complex_wavetable_alloc(N);
	
||||||| merged common ancestors
	double sampleFreq=J/((mSignalT->waveform[0])[J-1]-(mSignalT->waveform[0])[0]);
=======
	std::vector<Signal*>* dummy = new std::vector<Signal*>;
	
	double sampleFreq = 1/(2.0*I*((*mSignalsT)[0]->waveform[0][2] - (*mSignalsT)[0]->waveform[0][0]));

	
	
	
	gsl_fft_complex_workspace* complexWS =  gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* complexWT = gsl_fft_complex_wavetable_alloc(N);
	
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
	vec_d freq;
	for(int j=0; j<N/2; j++){
		freq.push_back(j*sampleFreq/N);
		freq.push_back(j*sampleFreq/N);
	} 
	for(int j=N/2; j<N; j++){
		freq.push_back((N/2-j)*sampleFreq/N);
		freq.push_back((N/2-j)*sampleFreq/N);
	} 
	
	vec_d amp;
	
		
<<<<<<< HEAD
	for (int i = 0; i < I; i++){
		
		Signal *sig = new Signal;
		amp = (*mSignalsT)[i]->waveform[1];
||||||| merged common ancestors
		for(size_t j=0;j<J;j++){
			Amp.push_back((mSignalT->waveform[1])[j]);
			Amp.push_back(0);
			
			//(output[0][i].waveform[1]).push_back(0);
		}	
		//(output->waveform[1])=*Amp;
		for(size_t j=J;j<N;j++){
			Amp.push_back(0);
			Amp.push_back(0);
		}
=======
	for (int i = 0; i < I; i++){
		
		Signal *sig = new Signal;
		amp = (*mSignalsT)[0]->waveform[1];
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
	
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
	
	for (int k = 0; k < signals; k++){
		
<<<<<<< HEAD
		std::vector<Template*>* dummy = new std::vector<Template*>;
		vec_d norms;
				
		for(int i=0; i<I; i++)
		{	
			Template temp = mTemplates[0][i];
		
			Template *op = new Template;
||||||| merged common ancestors
		Template* op=new Template;

		SNR = 0;

=======
		std::vector<Template*>* dummy = new std::vector<Template*>;
		vec_d norms;
				
		for(int i=0; i<I; i++)
		{	
			Template* temp = &mTemplates[0][i];
		
			Template *op = new Template;
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091

<<<<<<< HEAD
			op->param[0] = temp.param[0];
			op->param[1] = temp.param[1];
		
			realTemp = 0;
			imagTemp = 0;
		
			//not sure if this is looping through all of the elements correctly;
			//the alternative for the more sensible frequency packing is commented out
			for(int j = 0 ; j < J/2; j++){
||||||| merged common ancestors
		op->param[0] = temp->param[0];
		op->param[1] = temp->param[1];
		for(int j = 0; j <J/2; j++){
=======
			op->param[0] = temp->param[0];
			op->param[1] = temp->param[1];
		
			realTemp = 0;
			imagTemp = 0;
		
			//not sure if this is looping through all of the elements correctly;
			//the alternative for the more sensible frequency packing is commented out
			for(int j = 0 ; j < J/2; j++){
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
			
				/*
				
				noise = getASD(mSignalF->waveform[0][j]);
				op.waveform[0].push_back(mSignalF->waveform[0][j]);
			
				noisePower = noise*noise;
				*/	
				freq = (*mSignalsF)[k]->waveform[0][2*j];
				if(freq < 9){
					noise = 1.0;
				}else{
					noise = getASD((*mSignalsF)[k]->waveform[0][2*j]);
				}
				
				

				op->waveform[0].push_back((*mSignalsF)[k]->waveform[0][2*j]);
				op->waveform[0].push_back((*mSignalsF)[k]->waveform[0][2*j+1]);
			
				noisePower = noise*noise;
		
<<<<<<< HEAD
				realResult = (*mSignalsF)[k]->waveform[1][2*j] * temp.waveform[1][2*j] / noisePower;
				imagResult = -(*mSignalsF)[k]->waveform[1][2*j+1] * temp.waveform[1][2*j+1] / noisePower;
||||||| merged common ancestors
			realResult = mSignalF->waveform[1][2*j] * temp->waveform[1][2*j];
			imagResult = -mSignalF->waveform[1][2*j+1] * temp->waveform[1][2*j+1];
=======
				realResult = (*mSignalsF)[k]->waveform[1][2*j] * temp->waveform[1][2*j] / noisePower;
				imagResult = -(*mSignalsF)[k]->waveform[1][2*j+1] * temp->waveform[1][2*j+1] / noisePower;
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
			
<<<<<<< HEAD
				realTemp += temp.waveform[1][2*j] * temp.waveform[1][2*j] / noisePower;
				imagTemp += -temp.waveform[1][2*j] * temp.waveform[1][2*j+1] / noisePower;
||||||| merged common ancestors
			power = realResult*realResult + imagResult*imagResult;
=======
				realTemp += temp->waveform[1][2*j] * temp->waveform[1][2*j] / noisePower;
				imagTemp += -temp->waveform[1][2*j] * temp->waveform[1][2*j+1] / noisePower;
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
			
			
				op->waveform[1].push_back(realResult); //this should maybe be changed to SNR
				op->waveform[1].push_back(imagResult);
			}	
		
			realTemp = realTemp * 2 * df;
			imagTemp = imagTemp * 2 * df;
		
<<<<<<< HEAD
			norm = realTemp * realTemp + imagTemp * imagTemp;
			norm = sqrt(norm);
			norm=sqrt(norm);
			norms.push_back(norm);
||||||| merged common ancestors
		mSNR.push_back(SNR);
=======
			norm = realTemp * realTemp + imagTemp * imagTemp;
			norm = sqrt(norm);
		
			norms.push_back(norm);
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091

			dummy->push_back(op);
		}
		mNorms.push_back(norms);
		result->push_back(dummy);
	}
	mConResults = result;
	
	return;
}

void Extractor::fftInverse(std::vector<Template>* out1)
{
	int I, N;
	vec_d amp; 
	//vec_d freq;
	double dt;
	double norm;
<<<<<<< HEAD
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
||||||| merged common ancestors
	int K=mOriginalTime.size();
	I = (*mConResults).size();	
	N = (*mConResults)[0].waveform[0].size();	
=======
	int index;
	int fullSignal = mSignalT->waveform[1].size();
	std::cout<<fullSignal<<"  ,  "<<mSignalT->waveform[0].size()<<std::endl;
	int signals = (*mSignalsF).size();
	vec_d *time1 = &((*mSignalT).waveform[0]);
	vec_d time2;
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
	
<<<<<<< HEAD
	I = (*mConResults)[0]->size();	
	N = (*(*mConResults)[0])[0]->waveform[0].size();	

||||||| merged common ancestors
=======
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

>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N/2);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N/2);
	
	std::vector<std::vector<Template*>* > * final=new std::vector<std::vector<Template*>* >;
	
<<<<<<< HEAD
	
	std::cout<<mNorms[0][0]<<std::endl;
	
	for(int k = 0; k < signals; k++){
	
		std::vector<Template*> *result = new std::vector<Template*>;
		for(int i=0; i<I; i++)
		{		
			Template *temp = new Template;
			//freq = (*mConResults)[i].waveform[0];
			//dt = 2.0/(freq[N/2-1]);
			amp = (*(*mConResults)[k])[i]->waveform[1];
			gsl_fft_complex_inverse(&amp[0], 1, N/2, compWT, complexWS);		
			norm = mNorms[k][i];
			temp->param[0] = (*(*mConResults)[k])[i]->param[0];
			temp->param[1] = (*(*mConResults)[k])[i]->param[1];
			
			for(int n=0; n<N/2; n++)
			{
				temp->waveform[0].push_back((*mSignalsT)[k]->waveform[0][2*n]);
				temp->waveform[1].push_back((amp[2*n] * amp[2*n]) / norm);
			}
			
			
			result->push_back(temp);
			
		}
||||||| merged common ancestors
		amp = (*mConResults)[i].waveform[1];

		gsl_fft_complex_inverse(&amp[0], 1, N/2, compWT, complexWS);		

		Template *temp=new Template;

		norm = fAutoCorrComplex(&(*mTemplates)[i]);

		//std::cout<<norm<<std::endl;
=======
	clearSignalMemory();
	
	std::cout<<mNorms[0][0]<<std::endl;
	
	for(int k = 0; k < signals; k++){
	
		std::vector<Template*> *result = new std::vector<Template*>;
		for(int i=0; i<I; i++)
		{		
			Template *temp = new Template;
			//freq = (*mConResults)[i].waveform[0];
			//dt = 2.0/(freq[N/2-1]);
			amp = (*(*mConResults)[k])[i]->waveform[1];
			gsl_fft_complex_inverse(&amp[0], 1, N/2, compWT, complexWS);		
			norm = mNorms[k][i];
			temp->param[0] = (*(*mConResults)[k])[i]->param[0];
			temp->param[1] = (*(*mConResults)[k])[i]->param[1];
			
			for(int n=0; n<N/2; n++)
			{
				temp->waveform[0].push_back((*mSignalsT)[k]->waveform[0][2*n]);
				temp->waveform[1].push_back((amp[2*n] * amp[2*n]) / norm);
			}
			
			
			result->push_back(temp);
			
		}
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
		
		final->push_back(result);		
		
<<<<<<< HEAD
	}

	clearSignalMemory();
	clearConvolutionMemory();
	
	mToBeDeleted = final;
	
	for (int i = 0; i < I; i++){
		std::cout<<i<<std::endl;
		Template temp1;
		//Template *temp2 = new Template;
		temp1.param[0] = (*mTemplates)[i].param[0];
		temp1.param[1] = (*mTemplates)[i].param[1];
		for (int j = 0; j < signals; j++){
		std::cout<<j<<std::endl;
			for(int k = 0; k < N/2; k++){
				temp1.waveform[0].push_back((*((*final)[j]))[i]->waveform[0][k]);
				temp1.waveform[1].push_back((*((*final)[j]))[i]->waveform[1][k]);
				//temp2->waveform[0].push_back((*(final[2*j+1]))[i]->waveform[0][k]);
				//temp2->waveform[1].push_back((*(final[2*j+1]))[i]->waveform[1][k]);
			}
||||||| merged common ancestors
		for(int n=0; n<K; n++)
		{
			temp->waveform[1].push_back(amp[2*n]/norm);
=======
	}

	clearConvolutionMemory();
	
	mToBeDeleted = final;
	
	for (int i = 0; i < I; i++){
		std::cout<<i<<std::endl;
		Template *temp1 = new Template;
		//Template *temp2 = new Template;
		temp1->param[0] = (*mTemplates)[i].param[0];
		temp1->param[1] = (*mTemplates)[i].param[1];
		for (int j = 0; j < signals; j++){
		std::cout<<j<<std::endl;
			for(int k = 0; k < N/2; k++){
				temp1->waveform[0].push_back((*((*final)[j]))[i]->waveform[0][k]);
				temp1->waveform[1].push_back((*((*final)[j]))[i]->waveform[1][k]);
				//temp2->waveform[0].push_back((*(final[2*j+1]))[i]->waveform[0][k]);
				//temp2->waveform[1].push_back((*(final[2*j+1]))[i]->waveform[1][k]);
			}
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
		}
		
<<<<<<< HEAD
		out1->push_back(temp1);
		//outs2->push_back(*temp2);
||||||| merged common ancestors
		temp->waveform[0]=mOriginalTime;
		output->push_back(*temp);
=======
		out1->push_back(*temp1);
		delete temp1;
		//outs2->push_back(*temp2);
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
	}

<<<<<<< HEAD
	//mToBeDeleted2=out1;
	//out2 = outs2;
||||||| merged common ancestors
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	
=======
	//mToBeDeleted2=out1;
	//out2 = outs2;

	std::cout<<"asdfvrfusbecd"<<std::endl;
	std::cout<<(*out1)[0].waveform[0][349]<<std::endl;
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091

<<<<<<< HEAD
	std::cout<<"asdfvrfusbecd"<<std::endl;
	std::cout<<(*out1)[0].waveform[0][0]<<std::endl;
	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);		
	
||||||| merged common ancestors
=======
	
	
>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
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

	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N/2);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N/2);

	gsl_fft_complex_inverse(&op[0], 1, N/2, compWT, complexWS);

	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	

	return op[0];
}
<<<<<<< HEAD

void Extractor::clearSignalMemory(){
	//if (mSignalsT && mSignalsF){
		int signals = mSignalsT->size();
		for(int i = signals-1; i >= 0; i--){
			delete ((*mSignalsT)[i]);
			delete ((*mSignalsF)[i]);
			mSignalsT->erase(mSignalsT->begin() + i);
			mSignalsF->erase(mSignalsF->begin() + i);
		}
		delete (mSignalsT);
		delete (mSignalsF);
	//}
	//else if (mSignalsT){
	//	
	//}
	//else if (mSignalsF){
	//}
	return;
}

void Extractor::clearConvolutionMemory(){
	//if(mConResults)
	//{
		int signals = mConResults->size();
		int templates = (*mConResults)[0]->size();
	
		for(int i = signals-1; i >= 0; i--){
			for(int j = templates-1; j >= 0; j--){
				delete ((*(*mConResults)[i])[j]);
			}
			delete ((*mConResults)[i]);
		}
		delete(mConResults);
	//}
	return;
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

||||||| merged common ancestors
=======

void Extractor::clearSignalMemory(){
	if (mSignalsT && mSignalsF){
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
	else if (mSignalsT){
		
	}
	else if (mSignalsF){
	}
	return;
}

void Extractor::clearConvolutionMemory(){
	if(mConResults)
	{
		int signals = mConResults->size();
		int templates = (*mConResults)[0]->size();
	
		for(int i = signals-1; i >= 0; i--){
			for(int j = templates-1; j >= 0; j--){
				delete ((*(*mConResults)[i])[j]);
			}
			delete ((*mConResults)[i]);
		}
		delete(mConResults);
	}
	return;
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

>>>>>>> 79fd5d780623b369c13fa0895728bab5aab56091
