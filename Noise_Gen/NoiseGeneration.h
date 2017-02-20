#ifndef NOISEGENERATION_H
#define NOISEGENERATION_H

#define C_PI 3.1415926536

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

struct Complex{
	double real, imag;
};

class NoiseGenerator{
	
	public:
		
		NoiseGenerator();
		//NoiseGenerator(double (*func)(double), double);
		
		//void setNoiseCurve(double (*func)(double)); //Set noise curve function
		
		double getASD(double); //Return amplitude spectral density for given frequency
		double getMag(double); //Generate spectral magnitude for given frequency
		double getPhase();
		
		Complex getSample(double); //Get both phase and magnitude of frequency bin
		
		double noiseCurveALIGO(double);		
		
	private:
	
		//double (*fNoiseCurve)(double); //Function to use as the noise curve
		double fRelativeAmplitude; //Relative amplitude to scale the noise
	
};

//Returns normally distributed value using basic Box-Muller routine
double gaussianSample(double, double);

#endif //NOISEGENERATION_H