#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <cmath>
#include <zlib.h>
#include "gwReadWrite.h"
#include "gwDataTypes.h"
#include "gwSignalExtraction.h"

int main() 
{
	std::vector<Signal> printer;
	Signal signal,output;
	double Pi=4*atan(1);
	double N=10000;
	double f=10;
	double x0=300;
	double b=100;
	for (int i=1; i<N; i++)
	{
		signal.waveform[0].push_back(i);
		double a=sin(2*Pi*f*i)*exp(-pow((i-x0)/b,2)); //sinusoidal packet centred on 300
		signal.waveform[1].push_back(a);
	}


	Extractor * extractor=new Extractor();
	extractor->setSignal(&signal);

	extractor->fftComplex(&output);

	printer.push_back(output);

	saveSignals("testout.csv",&printer,csv);
	
	delete extractor;
	return 0;
}
