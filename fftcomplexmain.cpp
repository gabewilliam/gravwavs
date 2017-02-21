
#include "gwSignalExtraction.h"
#include "gwReadWrite.h"

int main(){
	
	delimiter f=tab;
	//declare the type of save file
	std::string file1="filterA.dat";//template file name
	std::vector<Template> Templates;
	loadTemplates(file1,&Templates,f);//load the templates
	
	std::string file2="SignalA.dat";//signal file name 
	std::vector<Signal> signal;	
	loadSignals(file2,&signal,f);//load signals



	
	Extractor Bob;//declare the extractor
	
	Bob.setSignal(&signal[0]);//load the extractor with signals and templates
	Bob.setTemplates(&Templates);
	
	Bob.fftComplex(&Templates);//transform the vectors into frequency domain
	Bob.fftComplex(&signal[0]);
	
	std::vector<Signal> Convolution;//decalre a signal vector to hold the convolutions
	
	Bob.fConvolutionComplex(&Convolution);//convolute the signal and templates
	Bob.fftInverseComplex(&Convolution);//transform back to time domain
	
	
	std::string file3="Out1.dat";//output vector
	saveSignals(file3,&Convolution,f);
	/*
	std::string file3="Out1.dat";
	saveSignals(file3,&Convolution,f);
	/*
	
	delimiter f=tab;
	std::string file3="Out1.dat";
	saveSignals(file3,&signal,f);
	
	delimiter f=tab;
	std::string file4="Out2.dat";
	saveTemplates(file4,&Templates,f);
*/	
	return 0;
}