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
		void setTemplates(std::vector<Template>* templates);	

		void fft(std::vector<Template>* output);

		void Convolution(std::vector<Template>* output);

		void fftInverse(std::vector<Template>* output);

		double getASD(double f);

		double fAutoCorrComplex(Template* temp);	

		bool loadCurve(std::string filename);
		
	private:
		Signal* mSignalT;
		Signal* mSignalF;
		std::vector<Template>* mTemplates;
		std::vector<Template>* mConResults;
		NoiseCurve fNoiseCurve;
};

#endif //GWSIGNALEXTRACTION_H
