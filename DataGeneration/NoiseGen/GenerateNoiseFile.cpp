#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "NoiseGeneration.h"

using namespace std;

int main(){

	srand(time(NULL));

	cout<<"Select type of noise:\r\n\t"
	<< "1) aLIGO Schutz\r\n\t"
	<< "2) aLIGO Zero Det High P\r\n\t"
	<< "3) aLIGO Zero Det Low P\r\n\t"
	<< "4) aLIGO NSNS Opt\r\n\t"
	<< "5) aLIGO No SRM\r\n\t"
	<< "6) aLIGO High Freq\r\n\t"
	<< "7) aLIGO BHBH 20 Deg\r\n";
	
	int choice;
	cin >> choice;
	
	NoiseGenerator *nGen;
	
		 if(choice == 1){	nGen =  new AligoSchutz();		}
	else if(choice == 2){	nGen = new AligoZeroDetHighP();	}
	else if(choice == 3){	nGen = new AligoZeroDetLowP();	}
	else if(choice == 4){	nGen = new AligoNsnsOpt();		}
	else if(choice == 5){	nGen = new AligoNoSrm();		}
	else if(choice == 6){	nGen = new AligoHighFreq();		}
	else if(choice == 7){	nGen = new AligoBhbh20Deg();	}
	else				{	nGen = new AligoSchutz();		}	
	
	ofstream oFile;
	oFile.open("noise.csv");
	
	ofstream oFile2;
	oFile2.open("noise2.csv");
	
	double fMax, fInc;
	
	vector<double> *freq = new vector<double>;
	vector<Complex> *noise = new vector<Complex>;

	cout << "Enter maximum frequency: \r\n";
	cin >> fMax;
	cout << "Enter frequency increment: \r\n";
	cin >> fInc;
	
	nGen->genSpectrum(freq, noise, fMax, fInc);
	
	for(int k=freq->size()/2; k < freq->size(); k++){
		
		oFile << freq->at(k) << "," << (noise->at(k)).real << "\r\n" << freq->at(k) << "," << (noise->at(k)).imag << "\r\n";

	}
	for(int k=0; k <= freq->size()/2; k++){
		
		oFile << freq->at(k) << "," << (noise->at(k)).real << "\r\n" << freq->at(k) << "," << (noise->at(k)).imag << "\r\n";

	}
	for(int k=0; k < freq->size(); k++){
		
		oFile2 << freq->at(k) << "," << (noise->at(k)).real << "," << (noise->at(k)).imag << "\r\n";

	}
	
	cout << "Output finished and written to 'noise.csv' \r\n";
	
	oFile.close();
	
}