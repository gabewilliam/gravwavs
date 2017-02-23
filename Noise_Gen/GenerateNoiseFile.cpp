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
	
	NoiseGenerator nGen;

	ofstream oFile;
	oFile.open("noise.csv");

	double fsamp, fInc;

	vector<double> *freq = new vector<double>;
	vector<Complex> *noise = new vector<Complex>;
	
	nGen.genSpectrum(freq, noise, 10000, 0.01);

	for(int k=0; k < freq->size(); k++){
		
		oFile << freq->at(k) << "," << (noise->at(k)).real << "," << (noise->at(k)).imag << "\r\n";
		
	}
	
	oFile.close();
	
}
