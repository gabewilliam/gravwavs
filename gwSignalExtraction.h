#include "gwDataTypes.h"

class Extractor
{
	public:
		Extractor();

		//Setup memeber functions
		void setSignal(Signal* input);
		void setTemplates(std::vector<Template>* templates);

		//Performs a time domain convolution of the signal 
		//against all templates in the template bank
		void tConvolution(std::vector<Signal>* output);

		//Performs a frequency domain convolution of the signal 
		//against all templates in the template bank
		void fConvolution(std::vector<Signal>* output);

		//Overloaded function to perform fft on a 
		//set of tenplates or signals
		void fft(std::vector<Template>* templatesFFT);
		void fft(Signal* signalFFT);

	private:
		Signal* mSignalT;
		Signal* mSignalF;

		std::vector<Template>* mTemplatesT;
		std::vector<Template>* mTemplatesF;
};
