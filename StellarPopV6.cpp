#include "Binary.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>


double CDF(double Z, double z);
double dMSFRdtdV(double z);
double Ez(double z);
double tIntegrand(double z, void * params);
double tLookback(double z);
double birthRate(double Z, double z, double M);


int main() {

	
	//Defines a few key constants in SI units
	const double mSolar = 1.989e30; //kg
	const double rSolar = 6.975e8; //m
	const double AU = 1.496e11; //m
	const double c = 3.0e8; //ms^-2
	const double G = 6.67408e-11;//m^3 kg^-1 s^-2
	const double yr = 365.25*24*60*60;//s
	


	//GENERATING A POPULATION OF BINARIES

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


	FILE * popFile;
	popFile = fopen("pop.csv","w");
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
	double tm;
	double mCandidates = 0;

	for(int i = N-1; i >= 0; i--) {

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

		std::cout << "\r" << std::setw(9) << std::right
				  << Nr << " binaries remaining." << std::flush;
		
	}

	
	double nCandidates = binaries.size();

	if(nCandidates==0) {
		std::cout<< "Sorry there are no binaries left." <<std::endl;
	}
	

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
			

	fclose(popFile);
			

	//Frees the memory associated with the RNGs
	gsl_rng_free(seedGen);
	gsl_rng_free(m1Gen);
	gsl_rng_free(qGen);
	gsl_rng_free(aGen);



	//FINDING THE MERGER RATE

	gsl_set_error_handler_off();//UNLIMITED POWER
	
	//Defines the maximum redshift which will be looked at
	double zMax = 25;

	//Defines the number of bins the universe is to be split into
	size_t nBins = 25;

	gsl_histogram * zHist = gsl_histogram_alloc(nBins);
	gsl_histogram_set_ranges_uniform(zHist,0,zMax);
	gsl_histogram * tHist = gsl_histogram_alloc(nBins);
	
	double tRanges[nBins+1];
	tRanges[0] = 0;

	double zLower, zUpper, tUpper;
	
	for(int j = 0; j < nBins; j++) {

		gsl_histogram_get_range(zHist, j, &zLower, &zUpper);
		tUpper = tLookback(zUpper);
		tRanges[j+1] = tUpper;

	}

	std::cout << tRanges[nBins] <<std::endl;

	gsl_histogram_set_ranges(tHist, tRanges, nBins+1);

	double tjMax, tjMin, tj, tiMax, tiMin, ti,
			zjMax, zjMin, zj, ziMax, ziMin, zi;
	double dtj, dti;
	size_t mergeIndexMax;
	size_t mergeIndexMin;
	double dNMergedtdV;

	double Zmax = 0.004;
	
	FILE * rateFile;
	rateFile = fopen("rates.csv","w");

	//Merge time bins
	for(int j = 0; j < nBins; j++) {

		dNMergedtdV = 0;

		gsl_histogram_get_range(zHist, j, &zjMin, &zjMax);
		zj = (zjMax +zjMin)/2;
		gsl_histogram_get_range(tHist, j, &tjMin, &tjMax);
		dtj = tjMax - tjMin;
	
		//Binaries
		for(int k = 0; k < nCandidates; k++) {
		
			tm = binaries[k].mergeTime()/yr;

			//Birth time bins
			for(int i = 0; i < nBins; i++) {

				//std::cout<<j<<" "<<k<<" "<<i<<std::endl;

				gsl_histogram_get_range(zHist, i, &ziMin, &ziMax);
				zi = (ziMax +ziMin)/2;
				ti = tLookback(zi);

				gsl_histogram_get_range(tHist, i, &tiMin, &tiMax);
				dti = tiMax - tiMin;
				tiMin += tm;
				tiMax += tm;
		
				//std::cout<<tiMax<<std::endl;	
				//gsl_histogram_find(tHist, tiMin, &mergeIndexMin);
				//gsl_histogram_find(tHist, tiMax, &mergeIndexMax);

				if(((tiMin>=tjMin)&&(tiMin<=tjMax))||((tiMax>=tjMin)&&(tiMax<=tjMax))) {
					dNMergedtdV += birthRate(Zmax, zi, mTotal)*(dti/dtj);
				}
											
			}
			
		}

		fprintf(rateFile, "%.15g,%.15g\n", zj, dNMergedtdV);

	}

	fclose(rateFile);

	gsl_histogram_free(zHist);
	gsl_histogram_free(tHist);

}


double CDF(double Z, double z) {

	const double ZSolar = 0.0134;
	
	double alpha = -1.16;
	double beta = 2;
	
	double a = alpha + 2;
	double Zbeta = pow((Z/ZSolar),beta);
	double tenbetaz = pow(10,(0.15*beta*z));
	double x = Zbeta*tenbetaz;

	return gsl_sf_gamma_inc_P(a,x);

}

double dMSFRdtdV(double z) {

	double num = pow((1+z),2.7);
	double denom = pow((1+(1+z)/2.9),5.6);
	
	return 0.015*num/denom;
	
}

double Ez(double z) {

	const double omegaLambda = 0.718;
	const double omegaM = 1 - omegaLambda;

	double Esquared = omegaM*pow((1+z),3) +omegaLambda;
	return sqrt(Esquared);

}

double tIntegrand(double z, void * params) {

	return 1/(Ez(z)*(1+z));

}

double tLookback(double z) {

	double tH =  4.55e17;//s
	const double yr = 365.25*24*60*60;//s
	tH = tH/yr; //Converts the Hubble time to years

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	
	double tL, tError;
	double param = 1; //I think this has to be done :(
	double epsRel = 1e-5;

	gsl_function tInt;
	tInt.function = &tIntegrand;
	tInt.params = &param;

	gsl_integration_qag(&tInt, 0, z, 0, epsRel, 1000, 6, w, &tL, &tError);

	gsl_integration_workspace_free(w);

	return tH*tL;

}

double birthRate(double Z, double z, double M) {

	return (CDF(Z,z)*dMSFRdtdV(z))/M;

}
