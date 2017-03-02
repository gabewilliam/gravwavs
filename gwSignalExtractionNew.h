#include "gwDataTypes.h"

class Extractor
{
	public:
		Extractor();

		//Setup functions
		void setSignalT(Signal* input);
		void setSignalF(Signal* input);
		void setTemplates(std::vector<Template>* templates);	

		void fft(std::vector<Template>* output);

		//Frequency convolution function
		void Convolution(std::vector<Template>* output);

		//FFT functions
		void fftInverse(std::vector<Template>* output);
		
	private:
		Signal* mSignalT;
		Signal* mSignalF;
		std::vector<Template>* mTemplates;
		std::vector<Template>* mConResults;
};
