#include "generate.h"
#include "gwReadWrite.h"

#include<cmath>
#include<stdio.h>

int main() {
	
	/*FILE * output;
	output = fopen("signal.txt", "w");*/

	std::vector<Signal> signalV;
	Signal signal;
	
	double //define parameters in natural units
	ma = 40*1.989e30, //mass of star a in kg
	mb = 30*1.989e30, //mass of star a in kg
	r = 140e6*3.0857e16, //distance from source in m
	theta = 0, //angle of emission measured from axis of orbit (i.e. orientation of orbit)
	phi = 0, //angle of emission around axis of orbit (this is just an arbitrary phase shift
	dt = 1e-2, //sample rate to be used in seconds
	Tm = 1; //Time at which black hole merge signal reaches earth
	
	//fprintf(output ,"%20.15s\t%20.15s\n", "time", "amplitude");
	
	for(int i = 0; i*dt<=Tm; i++){
		//printf("%g\n",r);
		//fprintf(output ,"%20.15g\t%20.15g\n",i*dt, generate(ma, mb, i*dt, r, theta, phi, Tm));
		signal.waveform[0].push_back(i*dt);
		signal.waveform[1].push_back(generate(ma, mb, i*dt, r, theta, phi, Tm));
	}

	std::cout << signal.waveform[0].size() << "\t" << signal.waveform[1].size() << std::endl; 

	signalV.push_back(signal);

	saveSignals("signal.txt", &signalV, csv);

	//fclose(output);
	
}
