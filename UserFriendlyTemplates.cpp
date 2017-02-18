//                         GRAVITATIONAL WAVE ASTRONOMY                            //
//                                Data Generation                                  //

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "templates.h"

using namespace std;

typedef double twoVect[2];

#define C_PRECISION 17
#define G_CONST 6.67E-11
#define c_CONST 3.0E8

double newtonianFrequency(double mA, double mB, double t){
	double A = 5.0/256.0;
	double B = pow((mA + mB), (1.0/3.0));
	double C = pow(c_CONST, 5.0);
	double D = pow(G_CONST, (5.0/3.0));
	double E = (A*B*C)/(mA*mB*t*D);
	
	double f = (1.0/M_PI)*pow(E, (3.0/8.0));
	
	return f;	
	
}

double initialSeparation(double mA, double mB, double lowLimit){	
	return pow((G_CONST*(mA + mB)/pow(2.0*M_PI*lowLimit, 2.0)), (1.0/3.0));
}

double waveAmplitude(double chirpM, double nextFreq, double lumD){
	double A = G_CONST/pow(c_CONST, 2.0);
    double B = (chirpM/lumD);
	double C = G_CONST/pow(c_CONST, 3.0);
	
	return 4.0*A*B*pow(C*M_PI*nextFreq*chirpM, (2.0/3.0));
}

double chirpMagnitude(double chirpM, double nextFreq){
	double A = (96.0/5.0);
	double B = pow(c_CONST, 3.0)/G_CONST;
	double C = nextFreq/chirpM;
	
	return A*B*C*pow(C*M_PI*nextFreq*chirpM, (8.0/3.0));
}

double timeTilMerge(double mA, double mB, double sep){
	double A = 5.0/256.0;
	double B = pow(sep, 4.0);
	double C = pow(c_CONST, 5.0);
	double D = pow(G_CONST, 3.0);
	double E = mA*mB*(mA + mB);
	
    return (A*B*C)/(D*E);
}

void binaryInspiral(double mA_z, 
				    double mB_z,
					double sFreq,
				    double lumD,
				    double sense,
					std::vector<Template> * tmps){
	
	// Initialise storage for template
	Template temp;
	
	temp.param[0] = mA_z;
	temp.param[1] = mB_z;
	
	double MPc = 3.086E24;
	double dLum = lumD*MPc;
	
	// Convert to real mass 
	double solarMass = 1.989E30;
	double mA = mA_z*solarMass;
	double mB = mB_z*solarMass;

	// Calculate the Chirp Mass
	double chirpMass = pow(mA*mB, (3.0/5.0))/pow((mA + mB), (1.0/5.0));
	
	// Find separation of binary for frequency within LIGO band
	double a = initialSeparation(mA, mB, sense/2.0);
	
	// Find time that binary signal will spend within LIGO band before merging
	double mergingTime = timeTilMerge(mA, mB, a);
	double tMerge = mergingTime;
	
	// Set the initial time of the signal
	double tLast = 0.0;
	double tCurrent = 0.0;
	
		// Loop over increasing time
	while (tMerge > 0.0){	
	
		// Update the time until merge
		tMerge = mergingTime - tCurrent;
		
		// Generate oscillating function with increasing frequency and amplitude
		double wFreq = newtonianFrequency(mA, mB, tMerge);
	
		// Define the variable time increment
		double dt = 1.0/sFreq;
	
		// Find the current time for the output
		tCurrent = tLast + dt;
		// Update the variable storing the previous time
		tLast = tCurrent;

		double waveForm = waveAmplitude(chirpMass, wFreq, dLum)*cos(2.0*M_PI*wFreq*tCurrent);
		
		temp.waveform[0].push_back(tCurrent);
		temp.waveform[1].push_back(waveForm);				 		
		
	}
	
	// Add to the vector of templates
	tmps->push_back(temp);
	
}

int main(){
		
	double lowM, highM, samplingFreq;
	cout << "Please enter the lower mass limit of the binary in solar masses: ";
	cin >> lowM;
	cout << endl << "Please enter the higher mass limit of the binary in solar masses: ";
	cin >> highM;
	cout << endl << "Please specify the desired sampling frequency:";
	cin >> samplingFreq;
	
	if (lowM > highM || lowM == 0.0 || samplingFreq == 0.0){
		cout << "Incorrect input parameters" << endl;
		return 1;
	}
	else{
		int Max = 5;
		int Min = 1;
		double M1, M2;	
		// Set luminosity distance in MPc
		double lumDist = 500.0;	
		// Set lower sensitivity of LIGO band
		double lowBound = 20.0;
		
		// Initialise vector of templates to save all generated data to
		std::vector<Template> templ;
		std::vector<Template> * t = &templ;
		
		// Loop over all mass pairs
		for (int Mp = Min; Mp <= Max; Mp++){
			// Primary mass
			M1 = double(Mp)*(highM-lowM)/double(Max) + lowM;
			for (int Ms = Min; Ms <= Max-Mp; Ms++){
				// Secondary mass
				M2 = double(Ms)*(highM-lowM)/double(Max) + lowM;
				binaryInspiral(M1, M2, samplingFreq, lumDist, lowBound, t);
				cout << "Mass 1: " << M1 << ", Mass 2: " << M2 << endl;
			}
		}
		
		// Set up name for output file
		ostringstream s1, s2, s3;
		s1 << lowM;
		s2 << highM;
		string del = "csv";
		s3 << del;
		std::string fileIdentity = s1.str() + "_" + s2.str() + "." + s3.str();
		
		// Save the templates, check the success of the program and exit
		if(saveTemplates(fileIdentity.c_str(), t, tab)){
			cout << "Templates stored successfully in file " << fileIdentity.c_str() << "." << endl;
			return 0;
		}
		else{
			cout << "There was an error" << endl;
			return 2;
		}
	}
	
}