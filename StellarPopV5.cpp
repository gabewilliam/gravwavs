#include "Binary.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


int main() {

	
	//Defines a few key constants in SI units
	double mSolar = 1.989e30; //kg
	double rSolar = 6.975e8; //m
	double AU = 1.496e11; //m
	double c = 3.0e8; //ms^-2
	double G = 6.67408e-11;//m^3 kg^-1 s^-2
	double yr = 365.25*24*60*60;//s


	/*Sets up a random number generator to seed the other two generators, and
	/ seeds this based on the system time.*/
	gsl_rng * seedGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(seedGen,(long unsigned int) time(NULL));
			

	/*Initialises random number generators to be used in sampling the masses
	/ and separartions the binary memebers.*/
	gsl_rng * m1Gen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(m1Gen, gsl_rng_get(seedGen));
	gsl_rng * qGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(qGen, gsl_rng_get(seedGen));
	gsl_rng * aGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(aGen, gsl_rng_get(seedGen));


	/*Defines the exponent of the m1 distribution.*/
	double alpha = 1.35;


	//Defines the upper and lower limits on mass ratio
	double qMax = 1;
	double qMin = 0.1;

	//Defines the upper limit on separation (and its logarithms)
	double aMax = 100*AU; //m
	double laMax = std::log10(aMax);

	
	/*Declares variables for the masses, mass ratio, radii, separation and 
	/ minimum separation.*/
	double m1, m2, q, r1, r2, la, a, aMin, laMin;


	FILE * popFile;
	popFile = fopen("pop.csv","w");
	std::vector <Binary> binaries;
	
		
	
	//Asks the user how many systems they want to generate	
	int N;
	std::cout << "How many binaries would you like to generate?" <<std::endl;
	std::cin >> N;
	
	int n = 0;
	
	std::cout << "Generating binaries..." << std::endl;

	int test = 0;

	/*Generates N binaries with component masses governed by the expressions
	/ cited in the Mandel and deMink paper. First m1 is sampled. Then the
	/ mass ratio is sampled and used to find m2.  Then, using m1 and m2, the 
	/ stellar radii are calculated from this, a lower limit on initial 
	/ separation is asigned, which is used to assign a separation value based
	/ on its logarithmic distribution. The process continues until N binaries
	/ are created which fit in the limit 20 < m << 300.*/
	while(n < N) {

		//Finds m1 in solar masses
		m1 = gsl_ran_pareto(m1Gen,alpha,1);
		
		//Finds the mass ratio, and thus m2
		q = gsl_ran_flat(qGen,qMin,qMax);
		m2 = m1*q;

		//Finds the stellar radii
		r1 = pow(m1,0.6);
		r2 = pow(m2,0.6);
		
		//Finds the lower limits on separation and its logarithm
		aMin = (r1*rSolar) + (r2*rSolar);
		laMin = std::log10(aMin);

		//Randomly samples a separation
		la = gsl_ran_flat(aGen,laMin,laMax);
		a = pow(10,la);	
		a = a/AU;
				
		binaries.push_back(Binary(m1, m2, a));
		n++;

		std::cout << "\r"<< std::setw(9) << std::right
					  << n << " binaries generated." << std::flush;
				
		
	}
	
	std::cout << std::endl;
	std::cout << "Sorting binaries..." << std::endl;
	int Nr = binaries.size();
	int check;
	int whatWentWrong[7] = {0, 0, 0, 0, 0, 0, 0};
	double tMerge;

	for(int i = N-1; i >= 0; i--) {

		check = binaries[i].checkCandidate();
		whatWentWrong[check]++;

		if(check > 0) {
			binaries.erase(binaries.begin() + i);
			Nr=Nr-1;
		}
		else {
			tMerge = binaries[i].mergeTime()/yr;
			fprintf(popFile,"%.15g,%.15g,%.15g\n",
					binaries[i].getMass(1),binaries[i].getMass(2),tMerge);	
		}

		std::cout << "\r" << std::setw(9) << std::right
				  << Nr << " binaries remaining." << std::flush;
		
	}

	
	double nCandidates = binaries.size();

	if(nCandidates==0) {
		std::cout<< "Sorry there are no binaries left." <<std::endl;
	}
	

	double tm, tmMin;
	int nextMerge;

	for(int i=0; i<nCandidates; i++) {
		
		tm = binaries[i].mergeTime();

		if((tm<tmMin) || (i==0)){
			nextMerge=i;
			tmMin=tm;
		}
	}

	std::cout << std::endl << std::endl;	
	std::cout << "Next binary merger will be: "<< std::endl;
	binaries[nextMerge].printGets();


	std::cout << std::endl;
	std::cout << "Binaries kept: " << whatWentWrong[0] << std::endl;
	std::cout << "Binaries remaining after each sorting phase:" << std::endl;

	int nRemaining = N;

	for(int i=1; i<7; i++) {

		nRemaining -= whatWentWrong[i];
		std::cout << i << ": " << nRemaining << std::endl;

	}

	std::cout << "Binaries lost at each sorting phase: " << std::endl;

	for(int i=1; i<7; i++) {
		std::cout << i << ": " << whatWentWrong[i] << std::endl;
	}
			

	fclose(popFile);
			

	//Frees the memory associated with the RNGs
	gsl_rng_free(seedGen);
	gsl_rng_free(m1Gen);
	gsl_rng_free(qGen);
	gsl_rng_free(aGen);


}