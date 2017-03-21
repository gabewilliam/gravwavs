/*
		GRAVITATIONAL WAVES GROUP STUDIES
			DATA GENERATION GROUP
		    NOISY SIGNAL GENERATOR
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "NoiseGen/gwNoiseGen.h"
#include "WaveGen/gwSigGenBasic.h"

struct IniParams{
	double fMin, m1, m2, lumDist, tArrival, fMax, phaseOffset;
};

IniParams inputParams();

int main(){
	
	IniParams inputs = inputParams();
	
	Parameters oParams;
	Parameters *P = &oParams;
	
	Signal sAmp;
	Signal *A = &sAmp;
	
	Signal sComp;
	Signal *C = &sComp;
	
	NoisySignal sNoise;
	NoisySignal *N = &sNoise;
	
	std::vector<Signal> sigs;
	std::vector<Signal> *S = &sigs;
	
	gwSimulateDetection(inputs.m1,
						 inputs.m2,
						 inputs.lumDist,
						 inputs.tArrival,
						 inputs.phaseOffset,
						 inputs.fMin,
						 inputs.fMax,
						 P,
						 A,
						 C,
						 N,
						 S);

	srand(time(NULL));

	std::cout<<"\r\n\r\n\tNOISE SETTINGS [* * * *]\r\n\r\nSelect type of noise:\r\n\t"
	<< "1) aLIGO Schutz\r\n\t"
	<< "2) aLIGO Zero Det High P\r\n\t"
	<< "3) aLIGO Zero Det Low P\r\n\t"
	<< "4) aLIGO NSNS Opt\r\n\t"
	<< "5) aLIGO No SRM\r\n\t"
	<< "6) aLIGO High Freq\r\n\t"
	<< "7) aLIGO BHBH 20 Deg\r\n\t"
	<< "8) White\r\n";
	
	int choice;
	std::cin >> choice;
	
	NoiseGenerator *nGen;
	
		 if(choice == 1){	nGen = new AligoSchutz();		}
	else if(choice == 2){	nGen = new AligoZeroDetHighP();	}
	else if(choice == 3){	nGen = new AligoZeroDetLowP();	}
	else if(choice == 4){	nGen = new AligoNsnsOpt();		}
	else if(choice == 5){	nGen = new AligoNoSrm();		}
	else if(choice == 6){	nGen = new AligoHighFreq();		}
	else if(choice == 7){	nGen = new AligoBhbh20Deg();	}
	else if(choice == 8){	nGen = new WhiteNoise();		}
	else				{	nGen = new WhiteNoise();		}	
	
	std::vector<double> *freq = new std::vector<double>;
	std::vector<Complex> *noise = new std::vector<Complex>;
	
	double fMax, fStep;
	
	fStep = oParams.fDF * C_CONST;
	
	nGen->setMinFreq(20.0);
	
	nGen->genSpectrum(freq, noise, 5000.0, fStep);
	std::ofstream oFile;
	oFile.open("NoisySignal.csv");
	
	double ff,rr,ii;
	
	for(int i=0; i < freq->size(); i++){

		ff=freq->at(i);
		
		if(i<sNoise.frequency.size()){
			rr=sNoise.compWave.at(i).real() += (noise->at(i)).real;
			ii=sNoise.compWave.at(i).imag() += (noise->at(i)).imag;
		}
		else{
			rr= (noise->at(i)).real;
			ii= (noise->at(i)).imag;
		}
		
		oFile << ff << "," << rr << "," << ii << "\r\n";

	}
	
	oFile.close();
	
	return 0;
}

IniParams inputParams(){ //Get parameter input
	
	double fMin, m1, m2, lumDist, tArrival, fMax, phaseOffset;
	
	//		SOURCE
	
	// Set masses of components (Solar Masses)
	std::cout << "\r\n\tSOURCE SETTINGS [* - - -]\r\n\r\n(1/3)\tEnter first mass (Solar Masses):\r\n";
	std::cin >> m1;
	std::cout << "(2/3)\tEnter second mass (Solar Masses):\r\n";
	std::cin >> m2;
	// Set luminosity distance (rads)
	std::cout << "(3/3)\tEnter luminosity distance (MPc):\r\n";
	std::cin >> lumDist;
	
	//		SIGNAL SETTINGS
	
	// Set total time of signal (s)
	std::cout << "\r\n\r\n\tSIGNAL SETTINGS [* * - -]\r\n\r\n(1/3)\tEnter max freq (Hz):\r\n";
	std::cin >> fMax;
	// Set time of arrival of signal (s)
	std::cout << "(2/3)\tEnter time of arrival of signal (s):\r\n";
	std::cin >> tArrival;
	// Set phase of wave at arrival (rads)
	std::cout << "(3/3)\tEnter phase offset of wave at arrival (rads):\r\n";
	std::cin >> phaseOffset;
	
	// 		FREQUENCIES
	
	// Set lower frequency sensitivity (Hz)
	std::cout << "\r\n\r\n\tFREQUENCY SETTINGS [* * * -]\r\n\r\n(1/1)\tEnter minimum detector frequency (Hz):\r\n";
	std::cin >> fMin;

	
	//Assign parameter inputs
	IniParams params;

	params.fMin=fMin;
	params.m1=m1;
	params.m2=m2;
	params.lumDist=lumDist;
	params.tArrival=tArrival;
	params.fMax=fMax;
	params.phaseOffset=phaseOffset;
	
	return params;
	
}

