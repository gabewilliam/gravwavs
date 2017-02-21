#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

//#include "gwReadWrite.h"
#include "NoiseGeneration.h"

using namespace std;

int main(){

	srand(time(NULL));

	Complex sample;
	
	NoiseGenerator ngen;
	
	ofstream oFile;
	oFile.open("noise.csv");
	
	double fs, df;
	
	fs = 10000;
	df = 0.01;
	
	vector<double> freq;
	vector<double> im;
	vector<double> re;
	
	
	
	int ind;
	for(int i=0; i<(20*df); i++){
		
		re[i] = 0.0;
		im[i] = 0.0;
		
		ind = (fs/df)-(20*df)+i;
		
		re[ind]=0.0;
		im[ind]=0.0;
		
	}
	for(int j=(20*df); j < (fs/df); j++){
		
		sample=ngen.getSample(0.01*j);
		
		re[j]=sample.real;
		re[fs/df-j]=sample.real;
		
		im[j]=sample.imag;
		im[fs/df-j]=-sample.imag;
		
	}
	
	for(int k=0; k < (fs/df); k++){
		
		oFile << freq[k] << "," << re[k] << im[k] << "\r\n";
		
	}
	
	oFile.close();
	
}
