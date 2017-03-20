#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "NoiseGen/gwNoiseGen.h"

using namespace std;

int main(){

		srand(time(NULL));

	cout<<"\r\n\r\n\tNOISE SETTINGS [* * * *]\r\n\r\nSelect type of noise:\r\n\t"
	<< "1) aLIGO Schutz\r\n\t"
	<< "2) aLIGO Zero Det High P\r\n\t"
	<< "3) aLIGO Zero Det Low P\r\n\t"
	<< "4) aLIGO NSNS Opt\r\n\t"
	<< "5) aLIGO No SRM\r\n\t"
	<< "6) aLIGO High Freq\r\n\t"
	<< "7) aLIGO BHBH 20 Deg\r\n\t"
	<< "8) White\r\n";
	
	int choice;
	cin >> choice;
	
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
	
	vector<double> *freq = new vector<double>;
	vector<Complex> *noise = new vector<Complex>;
	
	double fMax, fStep, fMin;
	
	cout<<"Enter time of signal: ";
	cin >> fStep;
	
	fStep = 1/fStep;
	
	cout<<"Enter max frequency: ";
	cin >> fMax;
	
	cout<<"Enter min frequency: ";
	cin >> fMin;
	nGen->setMinFreq(fMin);
	nGen->genSpectrum(freq, noise, fMax, fStep);

	ofstream oFile;
	oFile.open("Noise.csv");
	
	double ff,rr,ii;
	
	for(int i=0; i < freq->size(); i++){

		ff=freq->at(i);
		rr= (noise->at(i)).real;
		ii= (noise->at(i)).imag;
		
		oFile << ff << "," << rr << "," << ii << "\r\n";

	}
	
	oFile.close();
	
	return 0;
	
}