#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

//#include "gwReadWrite.h"
#include "NoiseGeneration.h"

using namespace std;

int main(){
	
	srand(time(NULL));

	Complex sample;
	
	ofstream outFile;
	
	NoiseGenerator ngen;
	
	outFile.open("noise.csv");
	
	for(int i=0; i<200; i++){
		
		outFile << 0 << "," << 0 << "\r\n";
		
	}
	for(int j=200; j<300000; j++){
		
		sample=ngen.getSample(0.01*j);
		
		outFile << sample.real << "," << sample.imag << "\r\n";
		
	}
	outFile << 0;
	outFile.close();
	
}
