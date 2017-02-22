/*	Simulates a single gravitational wave
 *	with given parameters over a given period
 */

#include "generateDG.h"
#include "gwReadWrite.h"

#include<cmath>
//#include<stdio.h>
#include<iostream>

int main() {

	std::vector<Template> templateV;
	Template temp;
	
	double //define parameters in SI units
	ma = 0, //mass of star a in kg
	mb = 0, //mass of star a in kg
	r = 140e6*3.0857e16*1e0, //distance from source in m
	theta = 0, //angle of emission measured from axis of orbit (i.e. orientation of orbit)
	phi = 0, //angle of emission around axis of orbit (this is just an arbitrary phase shift
	power = 0, //power of 2 for sampleing frequency
	dt = 0, //sample rate to be used in seconds
	Tm = 0, //Time at which black hole merge signal reaches earth
	dm = 0, //Mass incriments in kg
	m = 0, //mass lower limit
	M = 0, //mass upper limit
	solarMass = 1.989e30;
	
	bool twat = false;
	
	while(true){
		std::cout << "Please enter the lower mass limit of the binary in solar masses: ";
		std::cin >> m;
		std::cout << std::endl << "Please enter the higher mass limit of the binary in solar masses: ";
		std::cin >> M;
		std::cout << std::endl << "Please specify the desired sampling frequency as a power of 2: ";
		std::cin >> power;
		std::cout << std::endl << "Please enter the required mass increment in solar masses: ";
		std::cin >> dm;
		std::cout << std::endl << "Please enter the duration of each template in seconds: ";
		std::cin >> Tm;
		std::cout << std::endl;
	
		if(m>M || m<=0 || dm <= 0 || m+dm > M){
			if(!twat){
				std::cout<<"Nice try buddy, but that aint gonna fly here."<<std::endl;
				twat = true;
			}
			else std::cout<<"Look, I can do this all day."<<std::endl;
		}
		else break;
	}
	
	m *= solarMass;
	M *= solarMass;
	dt = pow(2,-power);
	dm *= solarMass;
	
	for(int k = 0; k*dm<=M-m; k ++){
		
		ma = m+k*dm;
	
		temp.param[0] = ma;
		
		for(int j = k; j*dm<=M-m; j ++){
			
			mb = m+j*dm;
			
			temp.param[1] = mb;
		
			for(int i = 0; i*dt<=Tm; i++){
				temp.waveform[0].push_back(i*dt);
				temp.waveform[1].push_back(generate(ma, mb, i*dt, r, theta, phi, Tm));
			}

			templateV.push_back(temp);
			
			temp.waveform[0].clear();
			temp.waveform[1].clear();
			
		}
		
	}
	
	saveTemplates("templates.csv", &templateV, csv);

	
}
