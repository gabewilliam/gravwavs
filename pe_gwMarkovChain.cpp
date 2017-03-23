#include <iostream>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include "gwNoiseGen.h"
#include "pe_gwLikelihood.h"
#include "Interpolator.h"
#include "pe_gwSaveToFile.h"

double gaussian(double, double, double, double, double);
double prior(double, double, double, double, double, double, double);
double distancePrior(double, double, double);

int main() {

	//INSERT DATA FILE NAME HERE
	std::string fileName = "28_32_over_age.csv";

	//Sets up a random number generator to seed the other two generators, and
	//seeds this based on the system time
	gsl_rng * seedGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(seedGen,(long unsigned int) time(NULL));
			
	//Initialises random number generators to be used in the Markov-Chain
	//routine, and seeds them using values from the previous RNG
	gsl_rng * normGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(normGen, gsl_rng_get(seedGen));
	gsl_rng * rGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rGen, gsl_rng_get(seedGen));
	gsl_rng * startGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(startGen, gsl_rng_get(seedGen));

	//Takes mass limits for the prior function
	double mLower,mUpper,dLower,dUpper;

	//Sets prior limits for mass and distance
	mLower = 10;
	mUpper = 63;
	dLower = 450;
	dUpper = 550;
	
	//Takes an input for the number of samples used in the Monte Carlo routine
	std::cout<< "Enter the number of Monte-Carlo samples:" << std::endl;
	int N;
	std::cin >> N;

	//Declares the variables used throughout the routine
	double ma, mb, maProposal, mbProposal;
	double mChirp, mRatio, mChirpProposal, mRatioProposal;
	double distance, distanceProposal;
	long double p, pProposal, alpha;
	double nZeroMChirp, nZeroMRatio, nZeroDistance, rZero;
	double* mChirpArray = new double[N];
	double* mRatioArray = new double[N];
	double* distanceArray = new double[N];
	double acceptance=0;

	//Initialises the noise curve
	AligoZeroDetHighP noise;
	
	//Data signal (frequency domain)
	Signal dataSignal;

	//loadSignals requires std::vector
	std::vector< Signal > idataSignal; 

	//Loads signal into dt
	loadSignals(fileName, &idataSignal, csv);

	//Converts from vector of signals to signal
	dataSignal = idataSignal[0];

	int nFreq = dataSignal.waveform[0].size();	
	vec_d sf;

	//Fills a std::vector with the noise curve values
	for( int i = 0; i < nFreq; i++ ){
		sf.push_back( pow(noise.getASD(dataSignal.waveform[0][i]),2) );
	}

	//Sets the starting ma, mb and distance values for the routine as random integers
	/*
	ma = gsl_rng_uniform(startGen)*(mUpper-mLower)+mLower;
	mb = gsl_rng_uniform(startGen)*(mUpper-mLower)+mLower;
	distance = gsl_rng_uniform(startGen)*(dUpper-dLower)+dLower;
	*/
	
	ma = 32.0;
	mb = 28.0;
	distance = 500.0;
	
	//Constraint to ensure that ma is greater than mb
	if (mb > ma) {
		double mDummy = ma;
		ma = mb;
		mb = mDummy;
	}
	
	//Calculates initial chirp mass and mass ratio
	mChirp = pow((ma*mb),3./5)/pow((ma+mb),1./5);
	mRatio = mb/ma;

	//Initialises the interpolator for the chirp mass/mass ratio prior from CHE
	Interpolator chirpRatPrior = Interpolator("ChirpRatKernelGrid.csv","ChirpRatKernelProb.txt");

	//Start of MCMC Routine
	for(int i = 1; i <= N; i++) {

		//Calculates the posterior for the current point
		p = likelihood(ma,mb,500.0,fileName,dataSignal,sf)+log(distancePrior(distanceProposal,dUpper,dLower)*chirpRatPrior.estimateProb(mChirp,mRatio));	

		//Uses RNGs to samples random noise from proposal distributions for each parameter
		nZeroMChirp = gsl_ran_gaussian(normGen, 0.5);
		nZeroMRatio = gsl_ran_gaussian(normGen, 0.05);
		nZeroDistance =	gsl_ran_gaussian(normGen, 5);

		//Adds the random noise to the current parameter values
		mChirpProposal = mChirp + nZeroMChirp;
		mRatioProposal = mRatio + nZeroMRatio;
		distanceProposal = distance + nZeroDistance;

		//Constraint to ensure that mass ratio is always less than 1
		if (mRatioProposal > 1) {
			mRatioProposal = 1./mRatioProposal;
		}

		//Calculates ma and mb for the proposed point
		maProposal = mChirpProposal*pow((1+mRatioProposal),1./5)*pow(mRatioProposal,-3./5);
		mbProposal = mChirpProposal*pow((1+mRatioProposal),1./5)*pow(mRatioProposal,2./5);

		//Calculates the posterior for the proposed point
		pProposal = likelihood(maProposal,mbProposal,500.0,fileName,dataSignal,sf)+log(distancePrior(distanceProposal,dUpper,dLower)*chirpRatPrior.estimateProb(mChirpProposal,mRatioProposal));

		//Finds the acceptance ratio (need to convert from logs which have been carried over to avoid rounding errors)
		alpha = exp(pProposal-p);

		//Random number between 0 and 1
		rZero = gsl_rng_uniform_pos(rGen);
		
		//Accepts new point if the acceptance ratio is greater than the random number
		if(alpha > rZero) {
		
			ma = maProposal;
			mb = mbProposal;
			mChirp = mChirpProposal;
			mRatio = mRatioProposal;
			//distance = distanceProposal;

			acceptance++;
		}
	
		//Prints a progress percentage to the screen
		if (i%(N/100)==0) {
			std::cout << (i*100)/N << "% complete." << std::endl;
		}

		//Records current parameter values in arrays
		mChirpArray[i-1] = mChirp;
		mRatioArray[i-1] = mRatio;
		distanceArray[i-1] = distance;
		
	}	

	//Prints the number of accepted proposals
	std::cout<<"Accepted: "<<acceptance<<std::endl;
	
	//Saves the arrays as columns to a file
	saveArraysToFile(mRatioArray,mChirpArray,distanceArray,1,N,"RawData.txt");

	//Frees the memory associated with the random
	//number generators and deallocates memory
	delete [] mChirpArray;
	delete [] mRatioArray;
	gsl_rng_free(seedGen);
	gsl_rng_free(normGen);
	gsl_rng_free(rGen);

}

//Defines a flat prior function which is a step of height 1, centred on the origin, with a width of w.
//Since obtaining the new prior from CHE, we no longer use this
double prior(double ma, double mb, double mUpper, double mLower, double d, double dUpper, double dLower) {

	if((ma < mUpper && ma > mLower) 
	&& (mb < mUpper && mb > mLower)
	&& (d < dUpper && d > dLower)) return 1;

	else return 0;

}

//The flat prior for the distance parameter only
double distancePrior(double d, double dUpper, double dLower) {
	if (d < dUpper && d > dLower) return 1;
	else return 0;
}






