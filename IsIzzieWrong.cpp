#include <iostream>
#include "gwReadWrite.h"
#include "gwSignalExtraction.h"

int main()
{
	std::vector<Signal> sigs;

	loadSignals("NonAdaptiveSignal.dat", &sigs, tab);

	Signal sig = sigs[0];

	int tl = sig.waveform[0].size();
	int al = sig.waveform[1].size();

	if(tl == al)
	{
		std::cout << "Izzie is right\t" << tl << "\t" << al << std::endl;
	}
	if(!(tl == al))
	{
		std::cout << "Izzie is wrong\t" << tl << "\t" << al << std::endl;
	}

	return 0;
}

