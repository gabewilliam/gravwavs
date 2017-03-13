#include "gwDataTypes.h"

struct NoiseCurve{//from gabe's noise curves
	double fMin, fMax;
	std::vector<double> freq;
	std::vector<double> asd;
};

class Extractor
{
	public:
		Extractor();

		//Setup memeber functions
		void setSignal(Signal* input);
		void setTemplates(std::vector<Template>* templates);
		std::vector<Template>* getTemplates() const{return mTemplatesF;};
		Signal* getSignal() const{return mSignalF;};
		
		
		//Overloaded function to perform fft on a 
		//set of tenplates or signals
		void fft(std::vector<Template>* templatesFFT);
		void fft(Signal* signalFFT);
		void fftComplex(std::vector<Template>* templatesFFT);
		void fftComplex(Signal* signalFFT);
		
		void fftInverse(std::vector<Signal>* signalsFFTI);
		void fftInverseComplex(std::vector<Signal>* signalsFFTI);

		void gateFilter(Signal*, double);
	
		void fftw(std::vector<Template>* templatesFFT);//redundant
		void fftw(Signal* signalFFT);
		void fftwInverse(std::vector<Signal>* signalFFT);
		
		//Performs time/frequency domain convolutions of the respective signal 
		//against all templates in the respective template bank
		void tConvolution(std::vector<Signal>* output);
		void fConvolution(std::vector<Signal>* output);
		void fConvolutionComplex(std::vector<Signal>* output);
	
		//from gabe's noise curves
		bool loadCurve(std::string);
		double getASD(double);
	
	private:
		Signal* mSignalT;
		Signal* mSignalF;
		
		vec_d mOriginalTime;
		std::vector<Template>* mTemplatesT;
		std::vector<Template>* mTemplatesF;
	
		NoiseCurve fNoiseCurve;//you can probably guess where this is from
	
		
};
