/*
		GRAVITATIONAL WAVES GROUP STUDIES
			DATA GENERATION GROUP
				NOISE GENERATOR
				
		Generates noise samples in the frequency domain with a given
		amplitude spectral density curve.
*/
#ifndef NOISEGENERATION_H
#define NOISEGENERATION_H

#define C_PI 3.1415926536

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

//Sruct for complex number storage
struct Complex{
	double real, imag;
};

class NoiseGenerator{
	
	public:
		
		NoiseGenerator();
		//NoiseGenerator(double (*func)(double), double);
		
		//void setNoiseCurve(double (*func)(double)); //Set noise curve function
		
		double getASD(double); //Return amplitude spectral density for given frequency
		double genMag(double); //Generate spectral magnitude for given frequency
		double genPhase();
		
		Complex genSample(double); //Get both phase and magnitude of frequency bin
		
		bool genSpectrum(std::vector<double>*, std::vector<Complex>* ,double, double); //Generate entire frequency spectrum
		
		double noiseCurveALIGO(double);		
		
	private:
	
		//double (*fNoiseCurve)(double); //Function to use as the noise curve
		double fRelativeAmplitude; //Relative amplitude to scale the noise
	
};

//Returns normally distributed value using basic Box-Muller routine
double gaussianSample(double, double);

#endif //NOISEGENERATION_H