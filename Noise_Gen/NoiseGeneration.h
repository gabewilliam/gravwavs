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
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
#include <istream>

//Sruct for complex number storage
struct Complex{
	double real, imag;
};
struct NoiseCurve{
	double fMin, fMax;
	std::vector<double> freq;
	std::vector<double> asd;
};

class NoiseGenerator{
	
	public:
		
		NoiseGenerator();
		
		virtual double getASD(double); //Return amplitude spectral density for given frequency
		double genMag(double); //Generate spectral magnitude for given frequency
		double genPhase();
		
		Complex genSample(double); //Get both phase and magnitude of frequency bin
		
		bool genSpectrum(std::vector<double>*, std::vector<Complex>* ,double, double); //Generate entire frequency spectrum
		
		bool loadCurve(std::string);
		
	protected:
	
		double fRelativeAmplitude; //Relative amplitude to scale the noise
		NoiseCurve fNoiseCurve;
		
};

class ALIGOSchutz: public NoiseGenerator{
	
	public:
		
		double getASD(double);
	
};

class ALIGOZeroDetHighP: public NoiseGenerator{
	
	public:
	
		ALIGOZeroDetHighP();
	
		double getASD(double);
	
};
//Returns normally distributed value using basic Box-Muller routine
double gaussianSample(double, double);

#endif //NOISEGENERATION_H