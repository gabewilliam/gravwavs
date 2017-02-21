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

	double fsamp, df;
	
	fsamp = 10000;
	df = 0.01;
	
	int cutoff = 20 / df;
	int N = (int) (fsamp / df);

	vector<double> freq;
	vector<double> im;
	vector<double> re;

	int ind;
	for(int i=0; i < cutoff; i++){
		cout << i << " Nsoasdf \r\n";
		re[i] = 0.0;
		im[i] = 0.0;			

		ind = ( N-cutoff+i );
		
		re[ind]=0.0;
		im[ind]=0.0;
		
	}
	for(int j=cutoff; j < N; j++){
		
		sample=ngen.getSample(0.01*j);
		
		re[j]=sample.real;
		re[N-j]=sample.real;
		
		im[j]=sample.imag;
		im[N-j]=-sample.imag;
		
	}
	for(int k=0; k < N; k++){
		
		oFile << freq[k] << "," << re[k] << im[k] << "\r\n";
		
	}
	
	oFile.close();
	
}
