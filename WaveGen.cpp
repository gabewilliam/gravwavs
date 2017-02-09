//                         GRAVITATIONAL WAVE ASTRONOMY                            //
//                                Data Generation                                  //

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

include namespace std;

int C_PRECISION  = 17;

double newtonianFrequency(double mA, double mB, double t){
	double G = 6.67E-11.0;
	double c = 3.0E8.0;
	double A = 5.0/256.0;
	double B = pow((mA + mB), (1.0/3.0));
	double C = pow(c, 5.0));
	double D = pow(G, (5.0/3.0));
	double E = (A*B*C)/(C*mA*mB*t);
	
	double f = (1.0/M_PI)*pow(E, (3.0/8.0));
	
	return f;	
	
}

int main(){
	
	// Set output precision
	cout.precision(C_PRECISION);
	
	
	
}