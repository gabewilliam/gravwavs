#include <iostream>
#include <cmath>
#include <cstdio>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "pe_gwLikelihood.h"
//#include "gwReadWrite.h"


//Prototype for the likelihood calculation function
double likelihood(double, double, std::string, double (*function)(double,double,std::string));
double gaussian(double, double, double, double, double);
double prior(double, double, double, double);


int main() {


	/*Sets up a random number generator to seed the other two generators, and
	/ seeds this based on the system time.*/
	gsl_rng * seedGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(seedGen,(long unsigned int) time);
	

	/*Initialises random number generators to be used in the Markov-Chain
	/ routine, and seeds them using values from the previous RNG.*/
	gsl_rng * normGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(normGen, gsl_rng_get(seedGen));
	gsl_rng * rGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rGen, gsl_rng_get(seedGen));
	gsl_rng * startGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(startGen, gsl_rng_get(seedGen));

	/*Takes an input for the standard deviation, sigma of the likelihood
	/ distribution, and then calculates a value to use in the distribution
	/ which is sampled to find the next x and y value [N(0,nSigma)]. The
	/ distribution from which the random number is drawn is 1/10 the width of
	/ the likelihood distribution. This is somewhat arbitrary, but means it
	/ won't be too wide. A slightly wider distribution would maybe give faster
	/ convergence.*/
	std::cout << "Enter the standard deviation, sigma, of the likelihood:" 
			  << std::endl;
	double sigma;
	std::cin >> sigma;
	
	double nSigma = sigma*1.989e30;

	/*
	//Takes inputs for the centre coordinates of the likelihood
	std::cout << "Enter the values (a,b) the likelihood is centred on:" 
			  << std::endl;
	double a,b;
	std::cout << "a:"<< std::endl;
	std::cin >> a;
	std::cout << "b:"<< std::endl;
	std::cin >> b;
	*/

	//Takes mass limits for the prior function
	double mLower,mUpper;
	
	std::cout<< "Enter the mass upper limit in solar masses:" << std::endl;
	std::cin >> mUpper;
	std::cout<< "Enter the mass lower limit in solar masses:" << std::endl;
	std::cin >> mLower;	
	
	mLower*=1.989e30;
	mUpper*=1.989e30;
	
	//Takes an input for the number of samples used in the Monte Carlo routine
	std::cout<< "Enter the number of Monte-Carlo samples:" << std::endl;
	int N;
	std::cin >> N;


	//Declares the variables used throughout the routine
	std::string fileName = "signal.txt";
	double ma, maPrime, mb, mbPrime, p, pPrime, alpha;
	double nZeroMA, nZeroMB, rZero;


	//Opens the output text file
	FILE * outFile;
	outFile = fopen("2DMonte.txt","w");
	fprintf(outFile,"%.15g\n",sigma);


	//Sets the starting ma and mb values for the routine as random integers in the range [-50,50]
	ma = gsl_rng_uniform(startGen)*(mUpper-mLower)+mLower;
	mb = gsl_rng_uniform(startGen)*(mUpper-mLower)+mLower;


	/*Loops over the number of iterations specified by the input. In each run,
	/ the likelihood function is evaluated at (x,y). Then, one of the RNGs is
	/ used to draw values from a normal distribution of width nSigma. The
	/ values obtained in this way are then used to find a new trial (x,y) pair
	/ [(x',y') in Veitch's notation]. The likelihood at this new point is
	/ then found. This is used to evaluate the acceptance value, alpha, which
	/ is compared to a random number drawn from a uniform distribution, to
	/ decide if (x',y') is accepted as the new (x,y). Either way, the (x,y)
	/ is output to the file.*/
	for(int i = 1; i <= N; i++) {

		//std::cout<<x<<std::endl;
		p = likelihood (ma,mb,fileName,PdhFunction)*prior(ma,mb,mUpper,mLower);

		if(i==1) {
			std::cout<<ma<<std::endl;
			std::cout<<likelihood(ma,mb,fileName,PdhFunction)<<std::endl;
		}

		nZeroMA = gsl_ran_gaussian(normGen, nSigma);
		nZeroMB = gsl_ran_gaussian(normGen, nSigma);

		maPrime = ma + nZeroMA;
		mbPrime = mb + nZeroMB;
		pPrime = likelihood(maPrime,mbPrime,fileName,PdhFunction)*prior(maPrime,mbPrime,mUpper,mLower);

		alpha = pPrime/p;

		//std::cout<<alpha<<std::endl;
		//std::cout<<pPrime<<std::endl;

		rZero = gsl_rng_uniform_pos(rGen);
		//if(i==1) std::cout<<alpha<<"\t"<<rZero<<std::endl;
		if(alpha > rZero) {
		
			ma = maPrime;
			mb = mbPrime;

		}
		
		if (i%100==0){
			fprintf(outFile,"%.15g,%.15g\n",ma,mb);
		}

	}


	/*Closes the output file, then frees the memory associated with the random
	/ number generators.*/
	fclose(outFile);

	gsl_rng_free(seedGen);
	gsl_rng_free(normGen);
	gsl_rng_free(rGen);
	

}

/*Defines a general likelihood function, which can take any function containing the appropriate arguments as 
/ input.*/
double likelihood(double ma, double mb, std::string file, double (*function)(double,double,std::string)) {

	double g;
    g = (*function)(ma,mb,file);
    return (g);

}


/*Defines a prior function which is a step of height 1, centred on the origin,
/ with a width of w.*/
double prior(double ma, double mb, double mUpper, double mLower) {

	if((ma < mUpper && ma > mLower) 
	&& (mb < mUpper && mb > mLower)) { return 1/((mUpper-mLower)*(mUpper-mLower)); }

	else { return 0; }

}
