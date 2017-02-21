#include "generate.h"

#include<cmath>
#include<stdio.h>

int main() {
	
	FILE * output;
	output = fopen("signal.txt", "w");
	
	
	double //define parameters in natural units
	ma = 40*1.989e30, //mass of star a in kg
	mb = 30*1.989e30, //mass of star a in kg
	r = 140e6*3.0857e16, //distance from source in m
	theta = 0, //angle of emission measured from axis of orbit (i.e. orientation of orbit)
	phi = 0, //angle of emission around axis of orbit (this is just an arbitrary phase shift)
	dt = 1e-6, //sample rate to be used in seconds
	Tm = 1; //Time at which black hole merge signal reaches earth
	
	//fprintf(output ,"%20.15s\t%20.15s\n", "time", "amplitude");
	
	for(int i = 0; i*dt<=2*Tm; i++){
		//printf("%g\n",r);
		fprintf(output ,"%20.15g\t%20.15g\n",i*dt, generate(ma, mb, i*dt, r, theta, phi, Tm));
	}

	fclose("signal.txt")
	
}
