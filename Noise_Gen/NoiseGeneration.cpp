#include "NoiseGeneration.h"

NoiseGenerator::NoiseGenerator(){
		fRelativeAmplitude=1;
}

double NoiseGenerator::getASD(double freq){
	
	return 0;
}

double NoiseGenerator::genMag(double freq){ //Get a sample magnitude
	
	double sigma = this->getASD(freq); //Use ASD to get frequency dependent S.D.
	double noise = gaussianSample(freq, sigma); //Generate single gaussian value
	
	return noise * fRelativeAmplitude;
	
}

double NoiseGenerator::genPhase(){ //Generate phase information
	
	//Phase is uniformly random about 360 degs
	double phase = 2.0 * (double)C_PI * ( rand() / (double)RAND_MAX );
	return phase;
	
}

Complex NoiseGenerator::genSample(double freq){ //Get real and imaginary parts of sample
	
	Complex sample;
	
	//Get phase and amplitude
	double mod = this->genMag(freq);
	double arg = this->genPhase();
	
	//Calculate real and imaginary parts from phase & amplitude
	sample.real = mod * cos( arg );
	sample.imag = mod * sin( arg );
	
	return sample;
	
}

bool NoiseGenerator::genSpectrum(std::vector<double>* freqs, std::vector<Complex>* noise, double fMax, double fInc){
	
	//Calculate maximum number of elements for given sampling frequency fMax
	int N = floor( (fMax / fInc) );
	std::cout<<N<<"\r\n";
	Complex sample;

	//Populate all elements with zeroes and fill frequency;
	sample.real = 0.0;
	sample.imag = 0.0;
	
	double freq = - fMax;
	
	for(int i=0; i < 2*N; i++){
		
		noise->push_back(sample);
		
		freqs->push_back(freq);
		
		freq += fInc;

	}
	std::cout<<freqs->at(N)<<"\r\n";
	for(int j=N; j < 2*N; j++){
		
		

		//Get a random complex sample
		sample=this->genSample(freqs->at(j));
		//Assign it to the spectrum
		noise->at(j)=sample;
		//Calculate complex conjugate of sample
		sample.imag *= -1;
		//Assign complex conjugate to 'positive' frequency axis
		noise->at(2*N-j)=sample;
		
	}
	
	return true;
	
}

bool NoiseGenerator::loadCurve(std::string filename){
	
	std::ifstream inFile;
	
	inFile.open(filename.c_str());
	
	double d;
	
	std::string line, element;

	NoiseCurve curve;

	while(getline(inFile,line)){
		
		std::istringstream iss1(line);
		
		while(getline(iss1, element, ',')){
			
			std::istringstream is(element);
			is >> d;
			curve.freq.push_back(d);
			
		}
		
		getline(inFile, line);
		
		std::istringstream iss2(line);
		
		while(getline(iss2, element, ',')){
			
			std::istringstream is(element);
			is >> d;
			curve.asd.push_back(d);	
			
		}
	}
	
	curve.fMin = curve.freq.front();
	curve.fMax = curve.freq.back();
	
	fNoiseCurve = curve;
	
	return true;
	
}

//Analytic fit for ALIGO noise curve
//	from Sathyaprakash & Schutz ~  https://arxiv.org/abs/0903.0338v1
double ALIGOSchutz::getASD(double f){
	
	double fs = 20.0; //Lowest frequency before fit blows up
	
	if(f > fs){
	
		double f0 = 215.0;
		double x = f/f0;
		
		double psd0 = 1E-49;
		double asd0 = sqrt(psd0);
		
		double x1,x2,x3,x3_1,x3_2,x3_3,asd,psd;
		
		x1 = pow(x,-4.14);
		x2 = -5.0 * pow(x, -2);
		
		x3_1 = -(x * x);
		x3_2 = 0.5 * pow(x,4);
		x3_3 = 0.5 * x * x;
		
		x3 = 111.0*(1.0+x3_1 +x3_2) / (1+x3_3);
		
		psd = x1 + x2 + x3;
		asd = sqrt(psd);
		
		return asd*asd0;
		
	}
	else{ return 0.0; }
	
}

ALIGOZeroDetHighP::ALIGOZeroDetHighP(){
	
	this->loadCurve("ZERO_DET_high_P.csv");
	
}

double ALIGOZeroDetHighP::getASD(double f){
	
	if((f < fNoiseCurve.fMin) || (f > fNoiseCurve.fMax)){
		return 0;
	}
	
	double fTest=0.0, fPrev=0.0, asdTest=0.0, asdPrev=0.0, grad=0.0, asd=0.0;
	
	//Find best asd to use for given frequency
	for(int i = 0; i <= fNoiseCurve.freq.size(); i++){
		
		fPrev = fTest;
		fTest = fNoiseCurve.freq[i];
		
		asdPrev = asdTest;
		asdTest = fNoiseCurve.asd[i];
		
		if(f=fTest){
			
			asd = fNoiseCurve.asd[i];
			
			i = fNoiseCurve.freq.size() + 1;
			
		}
		if(f<fTest){
			
			grad = (asdTest - asdPrev) / (fTest - fPrev);
			
			asd = asdPrev + ( (f-fPrev) * grad );
			
			i = fNoiseCurve.freq.size() + 1;
			
		}
		
	}
	
	return asd;
	
}

//Basic Box-Muller method for normally distributed random values
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