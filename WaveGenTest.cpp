//                             GRAVITATIONAL WAVE ASTRONOMY                            //
//                                    Data Generation                                  //

#include "gwSigGen.h"
#include <iostream>

int main(){
	
	//		MASSES
	// Set masses of components (Solar Masses)
	double M1;
	
	std::cout << "\tMASSES\r\n\r\n(1/2)\tEnter first mass (Solar Masses):\r\n";
	std::cin >> M1;
	
	double M2;
	
	std::cout << "(2/2)\tEnter second mass (Solar Masses):\r\n";
	std::cin >> M2;
	
	//		FURTHER SOURCE SETTINGS
	// Set distance of system (MPc)
	double Dist;
	
	std::cout << "\r\n\r\n\tFURTHER SOURCE SETTINGS\r\n\r\n(1/4)\tEnter distance of system (MPc):\r\n";
	std::cin >> Dist;
	
	// Angle of inclination of binary (rads)
	double psi;
	
	std::cout << "(2/4)\tEnter inclination angle of binary (radians):\r\n";
	cin >> psi;
	
	// Angle of incident wave (1)
	double theta;
	
	std::cout << "(3/4)\tEnter sky angle of incident wave (radians):\r\n";
	cin >> theta;
	
	// Angle of incident wave (2)
	double phi;
	
	std::cout << "(4/4)\tEnter plane angle of incident wave (radians):\r\n";
	cin >> phi;
	
	//		SIGNAL SETTINGS
	// Set total time of signal (s)
	double totalTime;
	
	std::cout << "\r\n\r\n\tSIGNAL SETTINGS\r\n\r\n(1/4)\tEnter total time of signal (s):\r\n";
	std::cin >> totalTime;
	
	// Set time of arrival of signal (s)
	double initTime;
	
	std::cout << "(2/4)\tEnter time of arrival of signal (s):\r\n";
	cin >> initTime;
	
	// Set phase of wave at arrival (rads)
	double initPhase;
	
	std::cout << "(3/4)\tPhase of wave at arrival (rads):\r\n";
	cin >> initPhase;
	
	// Set minimum detector frequency (Hz)
	double fMin;
	
	std::cout << "(4/4)\tEnter minimum detector frequency (Hz):\r\n";
	cin >> fMin;
	
	std::cout << "\r\nGenerating wave...\r\n";
	
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
							 theta,
							 psi,
							 phi,
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