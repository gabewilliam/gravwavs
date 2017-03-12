#include "gwSignalExtraction.h"
#include "gwReadWrite.h"

int main()
{
	std::cout << "Loading signal and template data..." << std::endl;
	
	std::vector<Signal> sigs;
	std::vector<Template> temps;
	
	if(!loadSignals("../data/SpecialNoiseSignalForEverybody.csv", &sigs, csv))
		return 0;

	if(!loadTemplates("../data/temp_50_30.csv", &temps, csv))
		return 0;

	int N, n;
	Signal sig;

	N = sigs.size();

	if(N > 1)
	{
		std::cout << "Multiple signals found." << std::endl;
		std::cout << "Which signal would you like to investigate (" << 1 << "-" << N << ")..." << std::endl;
		std::cin >> n;

		sig = sigs[n];
	}
	else
		sig = sigs[0];

	Extractor bill;
	
	if(bill.loadCurve("noisecurve.csv")){
		std::cout<<"noise successfully loaded"<<std::endl;
	}else{
		std::cout<<"no noise I'm afraid"<<std::endl;
		return 0;
	}
	
	bill.setSignalT(&sig);
	bill.setTemplates(&temps);
	
	std::cout << "splitting" << std::endl;
	
	bill.splitSignal();
	
	std::cout << "FFT" << std::endl;
	
	bill.fft();

	std::cout << "Performing testuency domain convolution against all templates..." << std::endl;

	std::vector<Template> conOut;

	bill.Convolution();

	saveTemplates("../data/conout.dat", &conOut, csv);

	std::cout << "Performing inverse fourier transform" << std::endl;
	
	
	std::vector<Template> output1;
	std::vector<Template> output2;

	bill.fftInverse(&output1);
	std::cout<<"doe breah"<<std::endl;
	std::cout<<output1[0].waveform[0][345]<<std::endl;
	
	std::cout << "Saving output to output.dat..." << std::endl;
	
	
	
	
	
	saveTemplates("../data/output1.dat", &output1, csv);
	//saveTemplates("../data/output2.dat", output2, csv);
	std::cout<<"here"<<std::endl;
	
	bill.clearOtherMemory();

	return 1;
}
