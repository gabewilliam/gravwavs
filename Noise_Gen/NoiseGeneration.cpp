#include "NoiseGeneration.h"

NoiseGenerator::NoiseGenerator(){
		fRelativeAmplitude=1;
}

double NoiseGenerator::getASD(double freq){
	
	//Get amplitude spectral density (standard deviation) from noise curve
	double sigma = this->noiseCurveALIGO(freq); //Standard deviation
	return sigma;
	
}

double NoiseGenerator::getMag(double freq){
	
	double sigma = this->getASD(freq);
	double noise = gaussianSample(freq, sigma);
	
	return noise * fRelativeAmplitude;
	
}

double NoiseGenerator::getPhase(){ //Generate phase information
	
	double phase = 2.0 * (double)C_PI * ( rand() / (double)RAND_MAX );
	return phase;
	
}

Complex NoiseGenerator::getSample(double freq){
	
	Complex sample;
	
	double mod = this->getMag(freq);
	double arg = this->getPhase();
	
	sample.real = mod * cos( arg );
	sample.imag = mod * sin( arg );
	
	return sample;
	
}

double NoiseGenerator::noiseCurveALIGO(double f){
	
	double fs = 20;
	double f0 = 215;
	double x = f/f0;
	
	double psd0 = 1E-49;
	double asd0 = sqrt(psd0);
	
	double x1,x2,x3,x3_1,x3_2,x3_3,asd,psd;
	
	x1 = pow(x,-4.14);
	x2 = -5.0 * pow(x, -2);
	
	x3_1 = -(x * x);
	x3_2 = 0.5 * pow(x,4);
	x3_3 = 0.5 * x * x;
	
	x3 = 111*(1+x3_1 +x3_2) / (1+x3_3);
	
	psd = x1 + x2 + x3;
	asd = sqrt(psd);
	
	return asd*asd0;
	
}

double gaussianSample(double f, double sigma){
 
	double x;
	double y;
	
	double result;
	
	y = 0; //If y=0, the log function explodes
	while( y == 0 ){y = ( rand() / ( (double)RAND_MAX ) );}

	x = cos( ( 2.0 * (double)C_PI ) * rand() / ( (double)RAND_MAX ) );
	
	//Result will be normally distributed random value
	result = sqrt( -2.0 * log( y ) ) * x;

	return (sigma * result);

}