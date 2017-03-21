//------------------- GRAVITATIONAL WAVE ASTRONOMY -------------------//
//------------------------ DATA GENERATION ---------------------------//

//----------- FREQUENCY DOMAIN GRAVITATIONAL WAVE GENERATOR ----------//
//----(Simplest model, no redshift or incident angle dependencies)----//
//------- Model simplified from predecessor to IMRPhenomB model, -----//
//-------- as described by P. Ajith et al. in paper available at -----// 
//---------------- [https://arxiv.org/pdf/0710.2335.pdf] -------------//

#ifndef GWSIGGEN_H
#define GWSIGGEN_H

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <complex>

//-- Include files
#include "gwReadWrite.h"

//-- Definitions for calculations
#define C_CONST 2.99792458E8
#define G_CONST 6.6740831E-11
#define SM_CONST 1.989E30
#define MPc_CONST 3.086E22
#define MHz_CONST 1E6

//-- Definitions for error warnings
#define SUCCESS 0
#define FAILURE 1


enum OutputType{
	 SIG,
	 TEMP
};

//----------------------------- STORAGE ------------------------------//

struct Parameters{
	
	//-- Mass parameters
	double fM1, fM2, fTotalMass, fSymmetricMass;
	
	//-- Distance parameter
	double fLumDistance;
	
	//-- IMRPhenomB model vector components:	
	//-- Transition frequencies and sigma
	double fMerge, fRing, fSigma, fCutoff;
	//-- Phase expansions
	double fPhaseExpansions[8];
	
	//-- Initial values
	double fArrivalTime, fPhaseOffset, fMinFreq;
	
	//-- Final values
	double fTotalTime, fMaxFreq;
	
	//-- Frequency increment
	double fDF;
	
	//-- Scaling amplitude of the wave	
	double fScalingAmplitude;
	
};	

//-- Struct for storing signal in format convenient for noise addition
struct NoisySignal{
	
	std::vector<double> frequency;
	std::vector< std::complex<double> > compWave;
	
};
//--------------------------------------------------------------------//
//--------------------- CONSTANT CALCULATIONS ------------------------//

//-- Calculate phase components
double computePhaseComponents(double x, 
							  double y, 
							  double z, 
							  double k,
							  Parameters *P){
	 
	double symmetricMass = P->fSymmetricMass;
	double totalMass = P->fTotalMass;
								  
	return (x*pow(symmetricMass, 2.0) + y*symmetricMass + z)
					  /(symmetricMass*pow(M_PI*totalMass, (5.0-k)/3.0));							  
								  
}

//-- Calculate transition frequencies and sigma
double computeTransitionComponents(double a, 
								   double b, 
								   double c, 
								   Parameters *P){
								  
	double symmetricMass = P->fSymmetricMass;
	double totalMass = P->fTotalMass;
	
	return (a*pow(symmetricMass, 2.0) + b*symmetricMass + c)
													  /(M_PI*totalMass);
}

//-- Calculate constant scaling amplitude
double scalingAmplitude(Parameters *P){
	
	double symmetricMass = P->fSymmetricMass;
	double totalMass = P->fTotalMass;
	double mergerFreq = P->fMerge;
	double lumDistance = P->fLumDistance;
	
	return pow(mergerFreq, -7.0/6.0)*pow(totalMass, 5.0/6.0)
								*pow(5.0*symmetricMass/24.0, 1.0/2.0)
									   /(lumDistance*pow(M_PI, 2.0/3.0));

}

//--------------------------------------------------------------------//
//------------------- WAVE LIFETIME CALCULATIONS ---------------------//

//-- Consistency factor
double wR(Parameters *P){
	
	double Sigma = P->fSigma;
	double mergeFreq = P->fMerge;
	double ringFreq = P->fRing;
	
	return (M_PI/2.0)*Sigma*pow(ringFreq/mergeFreq, -2.0/3.0);
	
}

//-- Lorentzian function for ringdown section
double lorentzianFunction(Parameters *P, double currentFreq){
	
	double Sigma = P->fSigma;
	double ringFreq = P->fRing;
	
	return Sigma/(2.0*M_PI*(pow((currentFreq - ringFreq), 2.0) 
											+ pow((Sigma/2.0), 2.0)));
	
}

//--------------------- Amplitude corrections:

//-- Inspiral
double inspiralAmp(Parameters *P, double currentFreq){
	double fDash = currentFreq/P->fMerge;	
	return pow(fDash, -7.0/6.0);	
}

//-- Merger
double mergerAmp(Parameters *P, double currentFreq){
	double fDash = currentFreq/P->fMerge;
	return pow(fDash, -2.0/3.0);
}

//-- Ringdown
double ringdownAmp(Parameters *P, double currentFreq){
	return wR(P)*lorentzianFunction(P, currentFreq);
}

//--------------------- Phase corrections:

double phaseTerms(Parameters *P, double currentFreq){
	
	double sumPhase = 2.0*M_PI*P->fArrivalTime*currentFreq
													  + P->fPhaseOffset;
	
	for(int i = 0; i<8; i++){
		
		sumPhase += P->fPhaseExpansions[i]
								 *pow(currentFreq, (double(i)-5.0)/3.0);
	}
	
	return sumPhase;
}

//--------------------------------------------------------------------//
//---------------------- AMPLITUDE CALCULATION -----------------------//

//-- Calculate amplitude for any given input frequency
//-- (Useful for Parameter Extraction subgroup)
double updatedAmplitude(Parameters *P, 
						double currentFreq,
						std::complex<double> *gWave){
	
	//-- Find the amplitude correction for the wave, windowing the 
	//-- signal between the lower frequency sensitivity and the cutoff
	double correctedAmp;
	
	if (currentFreq >= P->fMinFreq){
		
		if (currentFreq < P->fMerge){
			correctedAmp = inspiralAmp(P, currentFreq);
		}
		else if (currentFreq >= P->fMerge 
				 && currentFreq < P->fRing){
			correctedAmp = mergerAmp(P, currentFreq);
		}
		else if (currentFreq >= P->fMerge
				 && currentFreq < P->fCutoff){
			correctedAmp = ringdownAmp(P, currentFreq);
		}				
		else{
			correctedAmp = 0.0;
		}
		
	}
	else{
		correctedAmp = 0.0;
	}
	
	//-- Set up imaginary root of -1						
	const std::complex<double> J(0.0, 1.0);
	
	double wavePhase = phaseTerms(P, currentFreq);
	
	if (currentFreq == 0.0){			
		*gWave = 0.0;
	}else{		
		*gWave = P->fScalingAmplitude*correctedAmp*exp(-J*wavePhase);
	}
	
	double gwAmp = std::abs(*gWave);
	
	// Return the absolute magnitude in SI units (Hz^-1)
	return gwAmp/C_CONST;							
}
//--------------------------------------------------------------------//
//--------------------- WHOLE WAVE CALCULATION -----------------------//

//-- Generate entire signal 
void simulateGravitationalWave(Parameters *P, 
							   Signal *signalAmplitude,
							   Signal *signalComplex,
							   NoisySignal *signalForNoise,
							   std::vector<Signal> *signalVector){
	
	//-- Number of positive frequencies
	double nPosFreq = P->fMaxFreq/P->fDF + 1.0;
	//-- Number of data points
	double nDataPoints = 4.0*nPosFreq + 1.0;
	
	//-- Populate waveform with zeros to allow direct index reference
	for (int p=0; p < int(nDataPoints); p++){
		
		signalComplex->waveform[0].push_back(0.0);
		signalComplex->waveform[1].push_back(0.0);
		
		signalForNoise->frequency.push_back(0.0);
		signalForNoise->compWave.push_back(0.0);
		
	}
	
	//-- Set up complex number for the wave itself and its complex 
	//-- amplitude in Hz-1
	std::complex<double> mWave, hzWave;
	//-- Pointer to wave
	std::complex<double> * gW = &mWave;
	
	//-- Frequency in Hz and amplitude in Hz^-1
	double hzFreq, hzAmplitude;
	
	//-- Initialise frequency in m^-1
	double mFreq = 0.0;
	
	for (int i = 0; i <= int(nPosFreq); i++){
		
		mFreq = double(i)*P->fDF;
		
		//-- Update the wave amplitude
        hzAmplitude = updatedAmplitude(P, mFreq, gW);
				
		//-- Convert output back to SI units
		hzFreq = mFreq*C_CONST;	
		hzWave = mWave/C_CONST;
		
		//-- Store pure amplitude for Parameter Extraction
		signalAmplitude->waveform[0].push_back(hzFreq);
		signalAmplitude->waveform[1].push_back(hzAmplitude);
		
		//-- Store signal in format convenient for noise addition(1)
		signalForNoise->frequency[i+nPosFreq] = hzFreq;
		signalForNoise->compWave[i+nPosFreq] = hzWave;
		
		//-- Store signal in format convenient for Signal Extraction(1)
		signalComplex->waveform[0][2*i] = hzFreq;
		signalComplex->waveform[0][2*i+1] = hzFreq;
		
		signalComplex->waveform[1][2*i] = hzWave.real();
		signalComplex->waveform[1][2*i+1] = hzWave.imag();
		
		//-- Repeat the procedure for negative frequencies
		int k = nPosFreq - i;
		
		//-- Store signal in format convenient for noise addition(2)
		signalForNoise->frequency[k] = -hzFreq;
		signalForNoise->compWave[k] = conj(hzWave);
		
		//-- Store signal in format convenient for Signal Extraction(2)
		signalComplex->waveform[0][k] = -hzFreq;
		signalComplex->waveform[0][k+1] = -hzFreq;
		
		signalComplex->waveform[1][k] = hzWave.real();
		signalComplex->waveform[1][k+1] = -hzWave.imag();		
			
	}
	
	signalVector->push_back(*signalComplex);
	
}
//--------------------------------------------------------------------//
//--------------------------- CONSTRUCT ------------------------------//

void setParameters(double M1,
				   double M2,
				   double lumDistance,
				   double arrivalTime,
				   double phaseOffset,
				   double minFreq,
				   double maxFreq,
				   Parameters *P){

	//-- Assign stored variables in natural geometric units
	
	double massConversionFactor = SM_CONST*G_CONST/pow(C_CONST, 2.0);
	
	
	//-- Number of data points
	double nPoints = pow(2.0, 22.0);
	
	P->fM1 = M1*massConversionFactor;
	P->fM2 = M2*massConversionFactor;
	P->fLumDistance = lumDistance*MPc_CONST;
	P->fArrivalTime = arrivalTime*C_CONST;
	P->fPhaseOffset = phaseOffset;
	P->fMinFreq = minFreq/C_CONST;
	P->fMaxFreq = maxFreq/C_CONST;
	
	//-- Set frequency increment
	P->fDF = (4.0*P->fMaxFreq)/(nPoints - 5.0);
	
	P->fTotalTime = 1.0/P->fDF;
	
	//-- Set total mass
	P->fTotalMass = P->fM1 + P->fM2;
	
	//-- Set symmetric mass ratio
	P->fSymmetricMass = (P->fM1*P->fM2)/pow(P->fTotalMass, 2.0);
	
	//-- Set coefficients
	//-- Values taken from P. Ajith et al., accessed at 
	//-- [http://link.aps.org/doi/10.1103/PhysRevD.77.104017]
	
	//-- Phase expansions
	P->fPhaseExpansions[0] = 
	computePhaseComponents(1.7516E-1, 7.9483E-2, -7.2390E-2, 0., P);
	P->fPhaseExpansions[1] = 0.0;
	P->fPhaseExpansions[2] = 
	computePhaseComponents(-5.1571E1, -1.7595E1, 1.3253E1, 2., P);
	P->fPhaseExpansions[3] = 
	computePhaseComponents(6.5866E2, 1.7803E2, -1.5972E2, 3., P);
	P->fPhaseExpansions[4] = 
	computePhaseComponents(-3.9031E3, -7.7493E2, 8.8195E2, 4., P);
	P->fPhaseExpansions[5] = 0.0;
	P->fPhaseExpansions[6] = 
	computePhaseComponents(-2.4874E4, -1.4892E3, 4.4588E3, 6., P);
	P->fPhaseExpansions[7] = 
	computePhaseComponents(2.5196E4, 3.3970E2, -3.9573E3, 7., P);
	
	//-- Transitional frequencies and sigma
	P->fMerge = 
	computeTransitionComponents(2.9740E-1, 4.4810E-2, 9.5560E-2, P);
	P->fRing = 
	computeTransitionComponents(5.9411E-1, 8.9794E-2, 1.9111E-1, P);
	P->fSigma = 
	computeTransitionComponents(5.0801E-1, 7.7515E-2, 2.2369E-2, P);
	P->fCutoff = 
	computeTransitionComponents(8.4845E-1, 1.2848E-1, 2.7299E-1, P);
	
	P->fScalingAmplitude = scalingAmplitude(P);
	
}

//--------------------------------------------------------------------//
//------------------------- ENTIRE PROCESS ---------------------------//

//-- Simulate the detection, passing pointers to storage for the wave 
//-- amplitude and complex wave
void gwSimulateDetection(double M1,
						 double M2,
						 double lumDistance,
						 double arrivalTime,
						 double phaseOffset,
						 double minFreq,
						 double maxFreq,
						 Parameters *P,
						 Signal *A,
						 Signal *C,
						 NoisySignal *N,
						 std::vector<Signal> *S){
						  
	//-- Construct the system	
	setParameters(M1, M2, 
				  lumDistance, 
				  arrivalTime, 
				  phaseOffset, 
				  minFreq, 
				  maxFreq, P);
	
	//-- Compute the graviational wave signature
	simulateGravitationalWave(P, A, C, N, S);
	
}

//-- Simulate the detection and output the data to a file of the 
//-- designated type
int gwGenerateAndPack(double M1,
					  double M2,
					  double lumDistance,
					  double arrivalTime,
					  double phaseOffset,
					  double minFreq,
					  double maxFreq,
					  Parameters *P,
					  OutputType OUT){
	
	//-- Assign a vector of signals for passing to gwReadWrite functions
	std::vector<Signal> vectorSignal;
	//-- Pointer to vector
	std::vector<Signal> *S = &vectorSignal;
	
	//-- Assign a vector of templates for passing to gwReadWrite functions
	//-- (if OUT == TEMP)
	std::vector<Template> vectorTemplate;
	std::vector<Template> *T = &vectorTemplate;
	
	//-- Signals for passing through the simulation	
	Signal sigAmplitude, sigComplex;
	Signal *A = &sigAmplitude;
	Signal *C = &sigComplex;
	
	//-- Signal intended to be embedded in noise
	NoisySignal sigNoisy;
	NoisySignal *N = &sigNoisy;
	
    gwSimulateDetection(M1, M2,
						lumDistance,
						arrivalTime,
						phaseOffset,
						minFreq,
						maxFreq,
						P, A, C, N, S);
	
	if (OUT == TEMP){
		
		//-- Set up a template to pass the signal into
		Template outputTemplate;
		outputTemplate.param[0] = M1;
		outputTemplate.param[1] = M2;
		
		//-- Harvest the last generated signal
		Signal lastSignal = vectorSignal.back();
		
		for (unsigned int indx = 0; 
			 indx < lastSignal.waveform[0].size(); 
											indx++){
			
			outputTemplate.waveform[0].push_back
											(lastSignal.waveform[0][indx]);
			outputTemplate.waveform[1].push_back
											(lastSignal.waveform[1][indx]);
			
		}
		
		vectorTemplate.push_back(outputTemplate);
	}
		
	//-- Set up the name of the data file
	std::ostringstream s1, s2, s3, s4;
	
	if (OUT == SIG){
		s1 << "SIG";
	}
	else if (OUT == TEMP){
		s1 << "TEMP";
	}
	
	s2 << M1 + M2;
	s3 << P->fMinFreq*C_CONST;
	std::string del = "csv";
	s4 << del;
	std::string fileIdentity = s1.str() + "_"
							   + s2.str() + "_" 
							   + s3.str() + "." 
							   + s4.str();
	
	//-- Save the signal/template using gwReadWrite functions
	if (OUT == SIG){
		
		if (saveSignals(fileIdentity, S, csv)){
			
			std::cout << "Signal stored successfully in the file " 
				      << fileIdentity << "." << std::endl;
				 
		return SUCCESS;
		
		}
	}
	else if (OUT == TEMP){
		
		if (saveTemplates(fileIdentity, T, csv)){
			
			std::cout << "Template stored successfully in the file " 
					  << fileIdentity << "." << std::endl;
				 
		return SUCCESS;
		
		}
	}

	std::cout << "There was an error" << std::endl;
	
	return FAILURE;	
						  
}

//--------------------------------------------------------------------//


#endif // GWSIGGEN_H