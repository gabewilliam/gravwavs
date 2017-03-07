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

	
	AligoNoSrm nGen;

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
	
	nGen.genSpectrum(freq, noise, fMax, fInc);
	
	

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