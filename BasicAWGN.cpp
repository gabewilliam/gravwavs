#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>

#include "gwReadWrite.h"

#define C_PI 3.1415926536

using namespace std;

double AWGN(double mu, double sigma, double relativeAmplitude);

int main(){
	
	srand(time(NULL));
	
	vector<Signal> *sigs = new vector<Signal>;
	
	string filename;
	
	cout << "Enter filename: \n"; 
	
	cin >> filename;
	
	cout << loadSignals(filename, sigs, csv) << "\n";
	
	vector<Signal> signals = *sigs;
	Signal sig;
	double sample;
	
	for(int i=0; i<signals.size(); i++){
	
		sig = signals[i];
	
		for(int j=0; j<sig.waveform[0].size(); j++){
		
			sig.waveform[1][j] += AWGN(0, 1, 1);
			cout << AWGN(0, 1, 1) <<"\n";
		
		}
		
		signals[i] = sig;
	}

	sigs = &signals;
	
	saveSignals("NoisySig.csv",sigs,csv);
	
}

double AWGN(double mu, double sigma, double relativeAmplitude){
 
	double x;
	double y;
	double result;
	int p;

	y = 0;

	while( y == 0 ){y = ( rand() / ( (double)RAND_MAX ) );}

	x = cos( ( 2.0 * (double)C_PI ) * rand() / ( (double)RAND_MAX ) );
	result = sqrt( -2.0 * log( y ) ) * x;

	return (mu + relativeAmplitude * sigma * result);

}