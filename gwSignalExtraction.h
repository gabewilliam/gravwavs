#ifndef GWSIGNALEXTRACTION_H
#define GWSIGNALEXTRACTION_H

#include "gwDataTypes.h"

#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <zlib.h>

struct NoiseCurve 
{
	double fMin, fMax;
	std::vector<double> freq;
	std::vector<double> asd;
};

class Extractor
{
	public:
		Extractor();

		void setSignalT(Signal* input);
		void setSignalF(Signal* input);
		void splitSignal();
		void setTemplates(std::vector<Template>* templates);	
		

		//void fft(std::vector<Template>* output);
		void fft();

		void Convolution();

		void fftInverse(std::vector<Template>* out1, std::vector<Template>* out2);

		double getASD(double f);

		double fAutoCorrComplex(Template* temp);	

		bool loadCurve(std::string filename);
		
		void clearSignalMemory();
		void clearConvolutionMemory();
		void clearOtherMemory();
		void clear();
		
	private:
		Signal *mSignalT;
		Signal *mSignalF;
		std::vector<Signal*>* mSignalsT;
		std::vector<Signal*>* mSignalsF;
		std::vector<Template>* mTemplates;
		std::vector<std::vector<Template*>* >* mConResults;
		std::vector<Signal*>* mNoiseConvolution;
		std::vector<std::vector<Template*>* >* mToBeDeleted;
		std::vector<double> mNorms;
		NoiseCurve fNoiseCurve;
};

#endif //GWSIGNALEXTRACTION_H
