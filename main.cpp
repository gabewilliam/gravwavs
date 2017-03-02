#include "gwDataTypes.h"
#include "gwReadWrite.h"

#include "gwSignalExtraction.h"

int main()
{
	std::cout << "Loading signal and template data..." << std::endl;
	
	std::vector<Signal> sigs;
	std::vector<Template> temps;
	
	if(!loadSignals("../data/sig_30_26.csv", &sigs, csv))
		return 0;

	if(!loadTemplates("../data/temp_30_26.csv", &temps, csv))
		return 0;

	int N, n;
	Signal sig;

	N = sigs.size();

	if(N > 1)
	{
		std::cout << N << " signals found." << std::endl;
		std::cout << "Which signal would you like to investigate (" << 1 << "-" << N << ")..." << std::endl;
		std::cin >> n;

		sig = sigs[n];
	}
	else
		sig = sigs[0];

	Extractor bill;
	
	bill.setSignalF(&sig);
	bill.setTemplates(&temps);

	std::cout << "Performing frequency domain convolution against all templates..." << std::endl;

	std::vector<Template> conOut;

	bill.Convolution(&conOut);

	saveTemplates("../data/conout.dat", &conOut, csv);

	std::cout << "Performing inverse fourier transform" << std::endl;
	
	std::vector<Template> output;

	bill.fftInverse(&output);
	
	std::cout << "Saving output to output.dat..." << std::endl;
	
	saveTemplates("../data/output.dat", &output, csv);

	return 1;
}

