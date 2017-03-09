//                             GRAVITATIONAL WAVE ASTRONOMY                            //
//                                    Data Generation                                  //

// This file loops over frequency, writing the subsequent Fourier-domain gravitational //
//                                   waveform to a file.                               //

#include "gwSigGen.h"

int main(){
	
	// Set masses of components (Solar Masses)
	double M1 = 30.0;
	double M2 = 26.0;
	// Set distance of system (MPc)
	double Dist = 500.0;
	// Set total time of signal (s)
	double totalTime = 20.0;
	// Set time of arrival of signal (s)
	double initTime = 10.0;
	// Set phase of wave at arrival (rads)
	double initPhase = 0.0;
	// Set minimum detector frequency (Hz)
	double fMin = 10.0;
	
	// Set up struct to store characteristic information of signal
	parameters PARAMS;
	parameters *P = &PARAMS;
	
	// Set the characteristic parameters
	setFundamentalParameters(initTime, 
							 totalTime,
							 initPhase, 
							 fMin, 
							 M1, 
							 M2, 
							 Dist, 
							 P);
							 
	// Set up struct to store vector of final data
	vector<Signal> signalVector;
	vector<Signal> *S = &signalVector;
	
	// Set up signal struct to store single signal
	Signal singleSignal;
	
	// Set up complex number representing the wave itself
	complex<double> complexWave = 0.0;
	complex<double> * gWave = &complexWave;
	
	// Initialise frequency and amplitude
	double freq = 0.0;
	double amp_Hz = 0.0;
	
	// Set variable for storing frequency in Hz
	double freq_Hz = 0.0;
	
	// Set maximum frequency
	double maxFreq = 4095.0/C_CONST;
	
	// Set frequency increment
	double df = 1.0/(totalTime*C_CONST);
	
	// Set the number of frequency bins to loop over
	double nFreq = maxFreq/df + 1.0;
	
	// Calculate the scaling amplitude of the wave
	double amp0 = effAmp(P);
	
	for (int i = 0; i < nFreq; i++){
		
		freq = double(i)*df;
		
		amp_Hz = updatedAmplitude(P, freq, gWave, amp0);
		
		// Convert frequency back to Hz
		freq_Hz = freq*C_CONST;		
			
		singleSignal.waveform[0].push_back(freq_Hz);	
		singleSignal.waveform[1].push_back(amp_Hz);
					
			
	}
	
	signalVector.push_back(singleSignal);
	
	// Set up the name of the data file
	ostringstream s1, s2, s3;
	s1 << M1;
	s2 << M2;
	string del = "csv";
	s3 << del;
	std::string fileIdentity = s1.str() + "_" + s2.str() + "." + s3.str();
	
	if(saveSignals(fileIdentity, S, csv)){
		cout << "Signal stored successfully in the file " << fileIdentity << "." << endl;
	return SUCCESS;
	}
	else{
		cout << "There was an error" << endl;
		return FAILURE;
	} 	
	
}