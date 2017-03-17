#include <iostream>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <string>
#include <iomanip>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include "pe_gwLikelihood.h"

double gaussian(double, double, double, double, double);
double prior(double, double, double, double, double, double, double);
void saveToFile(double [],double [],double [],int,int,std::string);

int main() {

	/*Sets up a random number generator to seed the other two generators, and
	/ seeds this based on the system time.*/
	gsl_rng * seedGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(seedGen,(long unsigned int) time(NULL));
			
	/*Initialises random number generators to be used in the Markov-Chain
	/ routine, and seeds them using values from the previous RNG.*/
	gsl_rng * normGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(normGen, gsl_rng_get(seedGen));
	gsl_rng * rGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rGen, gsl_rng_get(seedGen));
	gsl_rng * startGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(startGen, gsl_rng_get(seedGen));

	//Takes mass limits for the prior function
	double mLower,mUpper,dLower,dUpper;
	/*
	std::cout<< "Enter the mass upper limit in solar masses:" << std::endl;
	std::cin >> mUpper;
	std::cout<< "Enter the mass lower limit in solar masses:" << std::endl;
	std::cin >> mLower;	
	std::cout<< "Enter the distance upper limit in MPc:" << std::endl;
	std::cin >> dUpper;	
	std::cout<< "Enter the distance lower limit in MPc:" << std::endl;
	std::cin >> dLower;	
	*/
	mLower = 10;
	mUpper = 63;
	dLower = 45;
	dUpper = 55;
	
	//Takes an input for the number of samples used in the Monte Carlo routine
	std::cout<< "Enter the number of Monte-Carlo samples:" << std::endl;
	int N;
	std::cin >> N;

	//Declares the variables used throughout the routine
	std::string fileName = "26_30.csv";
	double ma, mb, maProposal, mbProposal;
	double mChirp, mRatio, mChirpProposal, mRatioProposal;
	double distance, distanceProposal;
	long double p, pProposal, alpha;
	double nZeroMChirp, nZeroMRatio, nZeroDistance, rZero;
	double* mChirpArray = new double[N];
	double* mRatioArray = new double[N];
	double* distanceArray = new double[N];
	double acceptance=0;

	//Sets the starting ma and mb values for the routine as random integers.
	/*
	ma = gsl_rng_uniform(startGen)*(mUpper-mLower)+mLower;
	mb = gsl_rng_uniform(startGen)*(mUpper-mLower)+mLower;
	distance = gsl_rng_uniform(startGen)*(dUpper-dLower)+dLower;
	*/

	ma = 30.0;
	mb = 26.0;
	distance = 50.0;

	if (mb > ma) {
		double mDummy = ma;
		ma = mb;
		mb = mDummy;
	}

	mChirp = pow((ma*mb),3./5)/pow((ma+mb),1./5);
	mRatio = mb/ma;

	/*Loops over the number of iterations specified by the input. In each run,
	/ the likelihood function is evaluated at (ma,mb). Then, one of the RNGs is
	/ used to draw values from a normal distribution of width nSigma. The
	/ values obtained in this way are then used to find a new trial (ma,mb) pair
	/ [(ma',mb') in Veitch's notation]. The likelihood at this new point is
	/ then found. This is used to evaluate the acceptance value, alpha, which
	/ is compared to a random number drawn from a uniform distribution, to
	/ decide if (ma',mb') is accepted as the new (ma,mb). Either way, the (ma,mb)
	/ is output to the file.*/
	for(int i = 1; i <= N; i++) {

		p = likelihood(ma,mb,distance,fileName,true)+log((pow(ma,2)/mChirp)*prior(ma,mb,mUpper,mLower,distance,dUpper,dLower));

		nZeroMChirp = gsl_ran_gaussian(normGen, 0.5);
		nZeroMRatio = gsl_ran_gaussian(normGen, 0.05);
		nZeroDistance =	gsl_ran_gaussian(normGen, 1);

		mChirpProposal = mChirp + nZeroMChirp;
		mRatioProposal = mRatio + nZeroMRatio;
		distanceProposal = distance + nZeroDistance;

		maProposal = mChirpProposal*pow((1+mRatioProposal),1./5)*pow(mRatioProposal,-3./5);
		mbProposal = mChirpProposal*pow((1+mRatioProposal),1./5)*pow(mRatioProposal,2./5);

		if (mRatioProposal > 1) {
			mRatioProposal = 1./mRatioProposal;
		}
		
		pProposal = likelihood(maProposal,mbProposal,distanceProposal,fileName,false)+log((pow(maProposal,2)/mChirpProposal)*prior(maProposal,mbProposal,mUpper,mLower,distanceProposal,dUpper,dLower));

		alpha = exp(pProposal-p);

		rZero = gsl_rng_uniform_pos(rGen);
		
		if(alpha > rZero) {
		
			ma = maProposal;
			mb = mbProposal;

			mChirp = mChirpProposal;
			mRatio = mRatioProposal;

			distance = distanceProposal;

			acceptance++;
		}

		if (i%(N/100)==0) {
			std::cout << (i*100)/N << "% complete." << std::endl;
		}

		mChirpArray[i-1] = mChirp;
		mRatioArray[i-1] = mRatio;
		distanceArray[i-1] = distance;
		
	}	

	std::cout<<"Accepted: "<<acceptance<<std::endl;
	
	saveToFile(mRatioArray,mChirpArray,distanceArray,1,N,"RawData.txt");
	//saveToFile(mRatioArray,mChirpArray,distanceArray,1000,N,"ThinnedData.txt");

	/*Frees the memory associated with the random
	/ number generators and deallocates memory.*/
	delete [] mChirpArray;
	delete [] mRatioArray;
	gsl_rng_free(seedGen);
	gsl_rng_free(normGen);
	gsl_rng_free(rGen);

}

/*Defines a prior function which is a step of height 1, centred on the origin,
/ with a width of w.*/
double prior(double ma, double mb, double mUpper, double mLower, double d, double dUpper, double dLower) {

	if((ma < mUpper && ma > mLower) 
	&& (mb < mUpper && mb > mLower)
	&& (d < dUpper && d > dLower)) return 1;

	else return 0;

}

void saveToFile(double parameterA[], double parameterB[], double parameterC[], int lag, int size, std::string fileName) {
	
	//Opens the output text file
	FILE * outFile;
	outFile = fopen(fileName.c_str(),"w");

	for(int i = 0; i < size; i++){
		if (i%lag==0){
			fprintf(outFile,"%.15g,%.15g,%.15g\n",parameterA[i],parameterB[i],parameterC[i]);
		}
	}
	
	//Closes the output file
	fclose(outFile);
}






