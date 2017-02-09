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

double initialSeparation(double mA, double mB){
	double G = 6.67E-11;
	return pow((G*(mA + mB)/pow(20.0*M_PI, 2.0)), (1.0/3.0));
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
	
	// Find separation of binary for frequency within LIGO band
	double a = initialSeparation(M1, M2);
	
	cout << "Initial separation is " << a << endl;
	
	// Find time that binary signal will spend within LIGO band before merging
	double mergingTime = timeTilMerge(M1, M2, a);
	
	cout << "Merging time is " << mergingTime << endl;
	
	// Set number of data points to generate
	int nData = 10000;
	
	// Initialise the starting amplitude of the wave
	double A = 1.0;
	
	// Loop over increasing time
	for (int i=0; i<=nData; i++){	
	
		// Find the current time for the output
		double tCurrent = i*(mergingTime/nData);

		// Update the time until merge
		double tMerge = mergingTime - tCurrent;
		
		// Generate oscillating function with increasing frequency and amplitude
		double freq = newtonianFrequency(M1, M2, tMerge);
		double waveForm = A*sin(2.0*M_PI*tCurrent*freq);
		
		// Send output to data file
		WaveData << tCurrent << "," 
			     << waveForm << endl;
		
		
	}
	
	return 0;
	
}