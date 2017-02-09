//                         GRAVITATIONAL WAVE ASTRONOMY                            //
//                                Data Generation                                  //

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <time.h>

using namespace std;

// We should have these as globals somewhere at the root of the final code structure
int C_PRECISION  = 17;
const int C_DATASIZE = 1000;

typedef double dataVect[C_DATASIZE];

double gaussianNoise(double mean, 
					 double stdDev, 
					 double value){

		return (1.0/(stdDev*sqrt(2.0*M_PI)))*exp(-(value - mean)/(2.0*pow(stdDev, 2.0)));			   		
}

int main(){
	
	// Set output precision
	cout.precision(C_PRECISION);
    // Initialise output file
	ofstream Noise;
	Noise.open("Noise.txt", std::ios_base::app);

	// Create or import bare signal
	dataVect signalArray;
	dataVect timeArray;
	// Just create a simple siusoid for now
	double amp = 0.5;
	double freq = 20.0;
	for(int i=0; i<C_DATASIZE; i++){
		timeArray[i] = i*0.001;
		signalArray[i] = amp*sin(2.0*M_PI*freq*timeArray[i]);
	}
	// Create vector to store noise to add to signal
	dataVect addedNoise;
	// Create vector to store final noisy signal
	dataVect noisyData;
	
	// Initialise mean (parameter vector) and standard deviation (sigma) for Gaussian distribution
	double paramVect = 1.0;
	double sigma = 1.0;
	
	// Initialise added noise array of all zeroes
	for (int i=0; i<C_DATASIZE; i++){
		addedNoise[i] = 0.0;
	}
	
	bool addGaussian = true;
	bool addRandom = false;
	
	// Create random seed
	srand (time(NULL));
	double largestNoise = 5.0;
	double scalingFactor = double(RAND_MAX);
	
	// Gaussian noise
	if (addGaussian){
		for(int i=0; i<C_DATASIZE; i++){
			double randomVal = ((double(rand())/scalingFactor) - 0.5);
			addedNoise[i] += gaussianNoise(paramVect, sigma, randomVal);
		}
	}
	
	// Random noise
	if (addRandom){
		for(int i=0; i<C_DATASIZE; i++){
			addedNoise[i] += ((double(rand())/scalingFactor) - 0.5) * largestNoise;
		}
	}
	
	// Add complete noise vector to signal
	for(int i=0; i<C_DATASIZE; i++){
		noisyData[i] = signalArray[i] + addedNoise[i];
		// Send final results to data file
		Noise << timeArray[i] << "\t" 
			  << signalArray[i] << "\t" 
			  << addedNoise[i] << "\t" 
			  << noisyData[i] << endl;
	}
	
}
