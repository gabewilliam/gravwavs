//                             GRAVITATIONAL WAVE ASTRONOMY                            //
//                                    Data Generation                                  //

#include "gwSigGen.h"

int main(){
	
	// Set masses of components
	double M1 = 30.0;
	double M2 = 26.0;
	// Set distance of system
	double Dist = 500.0;
	// Set time of arrival of signal
	double initTime = 0.5;
	// Set phase of wave at arrival
	double initPhase = 0.0;
	// Set minimum detector frequency (Hz)
	double fMin = 20;
	
	// Set up struct to store characteristic information of signal
	parameters PARAMS;
	parameters *P = &PARAMS;
	
	// Set the characteristic parameters
	setFundamentalParameters(initTime, initPhase, fMin, M1, M2, Dist, P);
	
	cout << P->mu[0]*C_CONST << " " << P->mu[1]*C_CONST << " " << P->mu[3]*C_CONST << endl;
	cout << P->theta[2] << endl;
    // Set up struct to store vector of final data
	vector<Signal> signalVector;
	vector<Signal> *S = &signalVector;
	
	// Compute the gravitational wave
	int waveDone = gravitationalWave(P, S);
	
	
	if (waveDone == SUCCESS){
		
		// Set up the name of the data file
		ostringstream s1, s2, s3;
		s1 << M1;
		s2 << M2;
		string del = "csv";
		s3 << del;
		std::string fileIdentity = s1.str() + "_" + s2.str() + "." + s3.str();
		
		if(saveSignals(fileIdentity, S, csv)){
			cout << "Signal stored successfully" << endl;
		return SUCCESS;
		}
		else{
			cout << "There was an error" << endl;
			return FAILURE;
		} 	
	}else{
		return FAILURE;
	}
	
	
	return SUCCESS;
}