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

int C_PRECISION  = 17;

double newtonianFrequency(double mA, double mB, double t){
	double G = 6.67E-11;
	double c = 3.0E8;
	double A = 5.0/256.0;
	double B = pow((mA + mB), (1.0/3.0));
	double C = pow(c, 5.0);
	double D = pow(G, (5.0/3.0));
	double E = (A*B*C)/(mA*mB*t*D);
	
	double f = (1.0/M_PI)*pow(E, (3.0/8.0));
	
	return f;	
	
}

double initialSeparation(double mA, double mB, double lowLimit){
	double G = 6.67E-11;	
	return pow((G*(mA + mB)/pow(2.0*M_PI*lowLimit, 2.0)), (1.0/3.0));
}

double waveAmplitude(double chirpM, double nextFreq, double lumD){
	double G = 6.67E-11;
	double c = 3.0E8;
	double A = G/pow(c, 2.0);
    double B = (chirpM/lumD);
	double C = G/pow(c, 3.0);
	
	return 4.0*A*B*pow(C*M_PI*nextFreq*chirpM, (2.0/3.0));
}

double chirpMagnitude(double chirpM, double nextFreq){
	double G = 6.67E-11;
	double c = 3.0E8;
	double A = (96.0/5.0);
	double B = pow(c, 3.0)/G;
	double C = nextFreq/chirpM;
	
	return A*B*C*pow(C*M_PI*nextFreq*chirpM, (8.0/3.0));
}

double timeTilMerge(double mA, double mB, double sep){
	double G = 6.67E-11;
	double c = 3.0E8;
	double A = 5.0/256.0;
	double B = pow(sep, 4.0);
	double C = pow(c, 5.0);
	double D = pow(G, 3.0);
	double E = mA*mB*(mA + mB);
	
    return (A*B*C)/(D*E);
}

void binaryInspiral(double mA_z, 
				    double mB_z,
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
	
	// Declare the granularity of the signal
	double granularity = 10.0;
	
		// Loop over increasing time
	while (tMerge > 0.0){	
	
		// Update the time until merge
		tMerge = mergingTime - tCurrent;
		
		// Generate oscillating function with increasing frequency and amplitude
		double freq = newtonianFrequency(mA, mB, tMerge);
	
		// Define the variable time increment
		double dt = 1.0/(freq*granularity);
	
		// Find the current time for the output
		tCurrent = tLast + dt;
		// Update the variable storing the previous time
		tLast = tCurrent;

		double waveForm = waveAmplitude(chirpMass, freq, dLum)*cos(2.0*M_PI*freq*tCurrent);
		
		temp.waveform[0].push_back(tCurrent);
		temp.waveform[1].push_back(waveForm);				 		
		
	}
	
	// Add to the vector of templates
	tmps->push_back(temp);
	
}

int main(){
		
	// Initialise masses
	double maxM = 200.0;
	double minM = 1.0;
	double M1, M2;	
	// Set luminosity distance in MPc
	double lumDist = 500.0;	
	// Set lower sensitivity of LIGO band
	double lowBound = 20.0;
	
	// Initialise vector of templates to save all generated data to
	std::vector<Template> templ;
	std::vector<Template> * t = &templ;
	
	// Loop over all mass pairs
	for (int Mp = int(minM); Mp <= int(maxM); Mp++){
		// Primary mass
		M1 = double(Mp);
		for (int Ms = int(minM); Ms <= int(maxM - M1); Ms++){
			// Secondary mass
			M2 = double(Ms);
			binaryInspiral(M1, M2, lumDist, lowBound, t);
			cout << "Mass 1: " << M1 << ", Mass 2: " << M2 << endl;
		}
	}
	
	if(saveTemplates("Templates.dat", t, tab)){
		cout << "Templates stored successfully" << endl;
		return 0;
	}
	else{
		cout << "There was an error" << endl;
		return 1;
	}
	
}