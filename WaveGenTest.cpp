//                             GRAVITATIONAL WAVE ASTRONOMY                            //
//                                    Data Generation                                  //

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
	double fMin = 30.0;
	// Angle of inclination (rads)
	double incline = M_PI/3.0;
	
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
							 incline, 
							 P);
	
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
			cout << "Signal stored successfully in the file " << fileIdentity << "." << endl;
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