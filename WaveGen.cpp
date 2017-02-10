//                         GRAVITATIONAL WAVE ASTRONOMY                            //
//                                Data Generation                                  //

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

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

double arbitrarilyIncreasingAmplitude(double lastAmp){
	return 2.0*lastAmp;
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

int main(){
	
	// Initialise output file
	ofstream WaveData;
	WaveData.open("WaveData.csv", std::ios_base::app);
	// Set output precision
	cout.precision(C_PRECISION);
	WaveData.precision(C_PRECISION);
	
	// Initialise masses
	double solarMass = 1.989E30;
	double M1 = 30.0*solarMass;
	double M2 = 30.0*solarMass;
	
	// Calculate the Chirp Mass
	double chirpMass = pow(M1*M2, (3.0/5.0))/pow((M1 + M2), (1.0/5.0));
	
	// Set the luminosity distance
	double MPc = 3.086E22;
	double lumDist = 500.0*MPc;
	
	// Set lower sensitivity of LIGO band
	double lowBound = 20;
	
	// Find separation of binary for frequency within LIGO band
	double a = initialSeparation(M1, M2, lowBound/2.0);
	
	cout << "Initial separation is " << a << endl;
	
	// Find time that binary signal will spend within LIGO band before merging
	double mergingTime = timeTilMerge(M1, M2, a);
	double tMerge = mergingTime;
	
	cout << "Merging time is " << mergingTime << endl;
	
	// Set the initial time of the signal
	double tLast = 0.0;
	double tCurrent = 0.0;
	
	// Declare the granularity of the signal
	double granularity = 100000.0;
	
	// Loop over increasing time
	while (tMerge >= 0.0){	
	
		// Update the time until merge
		tMerge = mergingTime - tCurrent;
		
		// Generate oscillating function with increasing frequency and amplitude
		double freq = newtonianFrequency(M1, M2, tMerge);
	
		// Define the variable time increment
		double dt = 1.0/(freq*granularity);
	
		// Find the current time for the output
		tCurrent = tLast + dt;
		// Update the variable storing the previous time
		tLast = tCurrent;

		double waveForm = waveAmplitude(chirpMass, freq, lumDist)*cos(2.0*M_PI*freq*tCurrent);
		
		// Send output to data file
		WaveData << tCurrent << "," 
			     << waveForm << endl;
		
		
	}
	
	return 0;
	
}