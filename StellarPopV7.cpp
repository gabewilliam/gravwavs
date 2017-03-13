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

	
	/*Declares variables for the masses, mass ratio, radii, separation, 
	/ minimum separation and total mass generated.*/
	double m1, m2, q, r1, r2, la, a, aMin, laMin, mTotal;
	mTotal = 0;


	//FILE * popFile0;
	//popFile0 = fopen("pop0.csv","w");
	FILE * popFile1;
	popFile1 = fopen("pop1.csv","w");
	FILE * popFile2;
	popFile2 = fopen("pop2.csv","w");
	FILE * popFile3;
	popFile3 = fopen("pop3.csv","w");
	FILE * popFile3a;
	popFile3a = fopen("pop3a.csv","w");
	FILE * popFile4;
	popFile4 = fopen("pop4.csv","w");
	FILE * popFile4a;
	popFile4a = fopen("pop4a.csv","w");
	FILE * popFile5;
	popFile5 = fopen("pop5.csv","w");
	FILE * popFile5a;
	popFile5a = fopen("pop5a.csv","w");
	FILE * popFile6;
	popFile6 = fopen("pop6.csv","w");
	
	std::vector <Binary> binaries;
	
		
	
	//Asks the user how many systems they want to generate	
	int N;
	std::cout << "How many binaries would you like to generate?" <<std::endl;
	std::cin >> N;
	
	int n = 0;
	
	std::cout << "Generating binaries..." << std::endl;

	
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

		mTotal += (m1+m2);
				
		
	}
	
	std::cout << std::endl;
	std::cout << "Sorting binaries..." << std::endl;
	int Nr = binaries.size();
	int check;
	int whatWentWrong[7] = {0, 0, 0, 0, 0, 0, 0};
	//double a;
	double mCandidates = 0;

	for(int i = N-1; i >= 0; i--) {
		
		
		a = binaries[i].getSeparation();
		m1 = binaries[i].getMass(1);
		m2 = binaries[i].getMass(2);
		//fprintf(popFile0,"%.15g,%.15g,%.15g\n", m1, m2, a);

		if(binaries[i].checkMassRange() == false) {
			binaries.erase(binaries.begin() + i);
			Nr=Nr-1;
		}
		else if(true) {
			a = binaries[i].getSeparation();
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile1,"%.15g,%.15g,%.15g\n", m1, m2, a);
		}
		else if(binaries[i].checkRocheLobe() == false) {
			binaries.erase(binaries.begin() + i);
			Nr=Nr-1;
		}
		else if(true) {
			a = binaries[i].getSeparation();
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile2,"%.15g,%.15g,%.15g\n", m1, m2, a);
		}
		else if(binaries[i].checkHomogeneousMixing() == false) {
			binaries.erase(binaries.begin() + i);
			Nr=Nr-1;
		}	
		else if(true) {
			a = binaries[i].getSeparation();
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile3,"%.15g,%.15g,%.15g\n", m1, m2, a);
			binaries[i].evolveMainSequence();
			a = binaries[i].getSeparation();
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile3a,"%.15g,%.15g,%.15g\n", m1, m2, a);
		}
		else if(binaries[i].checkHomogeneousMixing() == false) {
			binaries.erase(binaries.begin() + i);
			Nr=Nr-1;
		}
		else if(true) {
			a = binaries[i].getSeparation();
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile4,"%.15g,%.15g,%.15g\n", m1, m2, a);
			binaries[i].evolveWolfRayet();
			a = binaries[i].getSeparation();
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile4a,"%.15g,%.15g,%.15g\n", m1, m2, a);
		}
		else if(binaries[i].checkPairInstability() == false) {
			binaries.erase(binaries.begin() + i);
			Nr=Nr-1;
		}
		else if(true) {
			a = binaries[i].getSeparation();
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile5,"%.15g,%.15g,%.15g\n", m1, m2, a);
			binaries[i].evolveSupernova();
			a = binaries[i].getSeparation();
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile5a,"%.15g,%.15g,%.15g\n", m1, m2, a);
		}
		else if(binaries[i].checkMergeTime() == false){
			binaries.erase(binaries.begin() + i);
			Nr=Nr-1;
		}	
		else {
			a = binaries[i].getSeparation();
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile6,"%.15g,%.15g,%.15g\n", m1, m2, a);
			mCandidates += (m1+m2);
		}

		/*
		check = binaries[i].checkCandidate();
		whatWentWrong[check]++;

		if(check > 0) {
			binaries.erase(binaries.begin() + i);
			Nr=Nr-1;
		}
		else {
			tm = binaries[i].mergeTime()/yr;
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile,"%.15g,%.15g,%.15g\n", m1, m2, tm);
			mCandidates += (m1+m2);
		}
		*/

		std::cout << "\r" << std::setw(9) << std::right
				  << Nr << " binaries remaining." << std::flush;
		
	}

	
	double nCandidates = binaries.size();

	if(nCandidates==0) {
		std::cout<< "Sorry there are no binaries left." <<std::endl;
	}
	
	double tm;
	double tmMin;
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

	std::cout << std::endl;
	std::cout << "Total Mass of Binaries = " << mTotal << std::endl;
	std::cout << "Total Mass of Candidates = " << mCandidates << std::endl;
			

	//fclose(popFile0);
	fclose(popFile1);
	fclose(popFile2);
	fclose(popFile3);
	fclose(popFile3a);
	fclose(popFile4);
	fclose(popFile4a);
	fclose(popFile5);
	fclose(popFile5a);
	fclose(popFile6);
			

	//Frees the memory associated with the RNGs
	gsl_rng_free(seedGen);
	gsl_rng_free(m1Gen);
	gsl_rng_free(qGen);
	gsl_rng_free(aGen);


}
