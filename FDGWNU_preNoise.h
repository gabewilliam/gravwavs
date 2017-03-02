//                             GRAVITATIONAL WAVE ASTRONOMY                            //
//                                    Data Generation                                  //

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <complex>

#include "gwReadWrite.h"

#define C_CONST 2.99792458E8
#define G_CONST 6.6740831E-11
#define SM_CONST 1.989E30
#define MPc_CONST 3.086E22

#define SUCCESS 0
#define FAILURE 1

using namespace std;

//-------------------------------------------------------------------------------------//

// Total mass
double comboMass(double m1, double m2){
	return m1 + m2;
}

// Symmetric mass ratio
double nu(double m1, double m2){
	return (m1*m2)/pow(comboMass(m1, m2), 2.0);
}

// Velocity 
double v(double M, double f){
	return pow(M_PI*M*f, 1.0/3.0);
}


//-------------------------------------------------------------------------------------//

// Struct for storage of mass and distance parameters 
// and intial parameters
typedef double theta_vec [3];
// Vector for storage of phase parameters
typedef double psi_vec [8];
// Vector for storage of transition parameters
typedef double mu_vec [4];

enum mode{
	
	Naked,
	PreNoise
};

struct parameters{
	
	// Vector of initial values:
	// init[0] = time of signal arrival
	// init[1] = phase of signal at arrival time
	// init[2] = minimum frequency of generated signal
	theta_vec init;
	
	// Vector of mass and distance parameters:
	// theta[0] = distance, D 
	// theta[1] = combined mass, M 
	// theta[2] = symmetrical mass ratio, n 
	theta_vec theta;
	
	// Vector of phase components
	psi_vec psi;
	
	// Vector of transition frequency components:
	// mu[0] = inspiral -> merger transition frequency, fMerge 
	// mu[1] = merger -> ringdown transition frequency, fRing
	// mu[2] = width of Lorentzian frequency distribution, sigma
	// mu[3] = maximum cutoff frequency, fCut
	mu_vec mu;
	
};

// Calculate phase components
double computePhaseComps(double x, double y, double z, double n, double M, double k){
	return (x*pow(n, 2.0) + y*n + z)/(n*pow(M_PI*M, (5.0-k)/3.0));
}

// Calculate transition frequencies and sigma
double computeTransitionComps(double a, double b, double c, double n, double M){
	return (a*pow(n, 2.0) + b*n + c)/(M_PI*M);
}

// Set the fundamental parameters of the system
void setFundamentalParameters(
	double init_time, 
	double init_phase, 
	double min_freq,
	double m_1, 
	double m_2, 
	double dist_MPc, 
	parameters *params
	){
		
	// Conversion factor for masses
	double convFact = G_CONST/pow(C_CONST, 2.0);
	// Convert masses into units of distance
	double m1 = m_1*convFact*SM_CONST;
	double m2 = m_2*convFact*SM_CONST;
	double M = comboMass(m1, m2);
	double n = nu(m1, m2);
	
	// Convert distance into meters
	double dist = dist_MPc*MPc_CONST;
	
	// Set initial time and phase parameters
	// (time in meters, frequency in meters^-1)
	params->init[0] = init_time*C_CONST;
	params->init[1] = init_phase;
	params->init[2] = min_freq/C_CONST;
	
	// Set mass and distance parameters
	params->theta[0] = dist;
	params->theta[1] = M;
	params->theta[2] = n;
	
	// Set phase parameters
	params->psi[0] = computePhaseComps(1.7516E-1, 7.9483E-2, -7.2390E-2, n, M, 0.);
	params->psi[1] = 0.0;
	params->psi[2] = computePhaseComps(-5.1571E1, -1.7595E1, 1.3253E1, n, M, 2.);
	params->psi[3] = computePhaseComps(6.5866E2, 1.7803E2, -1.5972E2, n, M, 3.);
	params->psi[4] = computePhaseComps(-3.9031E3, -7.7493E2, 8.8195E2, n, M, 4.);
	params->psi[5] = 0.0;
	params->psi[6] = computePhaseComps(-2.4874E4, -1.4892E3, 4.4588E3, n, M, 6.);
	params->psi[7] = computePhaseComps(2.5196E4, 3.3970E2, -3.9573E3, n, M, 7.);
	
	// Set transition parameters
	// fMerge
	params->mu[0] = computeTransitionComps(2.9740E-1, 4.4810E-2, 9.5560E-2, n, M);
	// fRing
	params->mu[1] = computeTransitionComps(5.9411E-1, 8.9794E-2, 1.9111E-1, n, M);
	// Sigma
	params->mu[2] = computeTransitionComps(5.0801E-1, 7.7515E-2, 2.2369E-2, n, M);
	// fCut
	params->mu[3] = computeTransitionComps(8.4845E-1, 1.2848E-1, 2.7299E-1, n, M);
	
}
//-------------------------------------------------------------------------------------//
// Lorentz equation
double Lorentzian(parameters *params, double f){
	return params->mu[2]/(2.0*M_PI*(pow((f - params->mu[1]), 2.0) + pow((params->mu[2]/2.0), 2.0)));
}
// Merger -> Ringdown
double omega_r(parameters *params){
	return pow(params->mu[1]/params->mu[0], -2.0/3.0)/Lorentzian(params, params->mu[1]);
}

//-------------------------------------------------------------------------------------//

// Preceding constant
double effAmp(parameters *params){
	return pow(params->mu[0], -7.0/6.0)*pow(params->theta[1], 5.0/6.0)*pow(5.0*params->theta[2]/24.0, 1.0/2.0)
			/(params->theta[0]*pow(M_PI, 2.0/3.0));
}

// Amplitude:
// Inspiral
double inspiralPhenomBAmp(parameters *params, double f){
	double fDash = f/params->mu[0];	
	return pow(fDash, -7.0/6.0);	
}
// Merger
double mergerPhenomBAmp(parameters *params, double f){
	double fDash = f/params->mu[0];
	return pow(fDash, -2.0/3.0);
}
// Ringdown
double ringdownPhenomBAmp(parameters *params, double f){
	double w_r = omega_r(params);
	return w_r*Lorentzian(params, f);
}

// Phase
double PhenomBPhase(parameters *params, double f){
	
	double sumPhase = 2.0*M_PI*params->init[0]*f;
	sumPhase += params->init[1];
	
	for(int i = 0; i<8; i++){
		sumPhase += params->psi[i]*pow(f, (double(i)-5.0)/3.0);
	}
	
	return sumPhase;
}

// Whole wave - signal
int gravitationalWave(parameters *params, vector<Signal> *sigs, int nData, mode outMode){
		
	// Set up storage of both real and imaginary components of the wave
	Signal totSig, reSig, imSig;
		
	vec_d addNoise[3];
	
	int nPoints = 2*nData;
	
	// Populate waveforms of arrays with all zeros to allow future direct index reference
	for (int k = 0; k < nPoints; k++){
		
		totSig.waveform[0].push_back(0.0);
		totSig.waveform[1].push_back(0.0);
		
		reSig.waveform[0].push_back(0.0);
		reSig.waveform[1].push_back(0.0);
		
		imSig.waveform[0].push_back(0.0);
		imSig.waveform[1].push_back(0.0);
		
		addNoise[0].push_back(0.0);
		addNoise[1].push_back(0.0);
		addNoise[2].push_back(0.0);
		
	}
	
	ofstream forNoise;
		
	// Set up the names of the file
	ostringstream s1, s2, s3;
	s1 << params->theta[0]/MPc_CONST;
	s2 << params->theta[1]*pow(C_CONST, 2.0)/(SM_CONST*G_CONST);
	string del = "csv";
	s3 << del;
	std::string sigFileIdentity = "dist_" + s1.str() + "_mass_" + s2.str() + "." + s3.str();		
	forNoise.open(sigFileIdentity.c_str(), std::ios_base::app);
	
	double df = 2.0*(params->mu[3] - params->init[2])/double(nData);
	
	double amp0 = effAmp(params);
	double freq = 0.0;
	
	double finalAmp, wavePhase;
	
	// Set up imaginary number J
	const complex<double> J(0.0, 1.0);
	// Set up complex number for the wave itself
	complex<double> gravWav;
	
	// Storage for frequency in Hz
	double freq_Hz;
			
	for (int i = 0; i < nData; i++){
		
		freq = i*df + params->init[2];
				
		if (freq < params->mu[0]){
			finalAmp = inspiralPhenomBAmp(params, freq);
		}
		else if (freq >= params->mu[0] && freq < params->mu[1]){
			finalAmp = mergerPhenomBAmp(params, freq);
		}
		else if (freq >= params->mu[1] && freq < params->mu[3]){
			finalAmp = ringdownPhenomBAmp(params, freq);
		}				
		else {
			finalAmp = 0.0;
		}
		
		wavePhase = PhenomBPhase(params, freq);
		
		gravWav = amp0*finalAmp*(cos(wavePhase)-J*sin(wavePhase));
		
		// Convert frequency back to Hz
		freq_Hz = freq*C_CONST;
		
		if (outMode == Naked){
			
			totSig.waveform[0][i] = freq_Hz;
			totSig.waveform[0][i] = freq_Hz;
			
			totSig.waveform[1][i] = gravWav.real();
			totSig.waveform[1][i] = gravWav.imag();
			
			int k = nData + i;
			
			totSig.waveform[0][k] = -freq_Hz;
			totSig.waveform[0][k] = -freq_Hz;
			
			totSig.waveform[1][k] = gravWav.real();
			totSig.waveform[1][k] = -gravWav.imag();			
			
		}else{
			
			reSig.waveform[0][i] = freq_Hz;		
			imSig.waveform[0][i] = freq_Hz;		
			
			reSig.waveform[1][i] = gravWav.real();		
			imSig.waveform[1][i] = gravWav.imag();
			
			addNoise[0][i] = freq_Hz;
			addNoise[1][i] = gravWav.real();
			addNoise[2][i] = gravWav.imag();
			
			int k = nData + i;
			
			reSig.waveform[0][k] = -freq_Hz;		
			imSig.waveform[0][k] = -freq_Hz;		
			
			reSig.waveform[1][k] = gravWav.real();		
			imSig.waveform[1][k] = -gravWav.imag();
			
			addNoise[0][k] = -freq_Hz;
			addNoise[1][k] = gravWav.real();
			addNoise[2][k] = -gravWav.imag();
		
		}
	}
		
	if (outMode == Naked){
		sigs->push_back(totSig);
	}else{
		sigs->push_back(reSig);
		sigs->push_back(imSig);
	}
	
	for (int l = 0; l < nPoints; l++){
		
		forNoise << addNoise[0][l] << ", " << addNoise[1][l] << ", " << addNoise[2][l] << endl;
		
	}
	
	return SUCCESS;
}

// Whole wave - template
int gravitationalWaveTemps(parameters *params, vector<Template> *temps, int nData, double m_1, double m_2){
		
	// Set up storage of both real and imaginary components of the wave
	Template totSig;
	
	// Set parameters 
	totSig.param[0] = m_1;
	totSig.param[1] = m_2;
	
	double df = (params->mu[3] - params->init[2])/double(nData);
	
	double amp0 = effAmp(params);
	double freq = 0.0;
	
	double finalAmp, wavePhase;
	
	// Set up imaginary number J
	const complex<double> J(0.0, 1.0);
	// Set up complex number for the wave itself
	complex<double> gravWav;
	
	// Storage for frequency in Hz
	double freq_Hz;
	
	for (int i = 0; i < nData; i++){
		
		freq = i*df + params->init[2];
		
		if (freq < params->mu[0]){
			finalAmp = inspiralPhenomBAmp(params, freq);
		}
		else if (freq >= params->mu[0] && freq < params->mu[1]){
			finalAmp = mergerPhenomBAmp(params, freq);
		}
		else if (freq >= params->mu[1] && freq < params->mu[3]){
			finalAmp = ringdownPhenomBAmp(params, freq);
		}				
		else {
			finalAmp = 0.0;
		}
		
		wavePhase = PhenomBPhase(params, freq);
		
		gravWav = amp0*finalAmp*(cos(wavePhase)-J*sin(wavePhase));
		
		// Convert frequency back to Hz
		freq_Hz = freq*C_CONST;
		
		totSig.waveform[0].push_back(freq_Hz);	
		totSig.waveform[0].push_back(-freq_Hz);
		
		totSig.waveform[1].push_back(gravWav.real());		
		totSig.waveform[1].push_back(gravWav.imag());
		
	}
	
	temps->push_back(totSig);
	
	return SUCCESS;
}
	

	

