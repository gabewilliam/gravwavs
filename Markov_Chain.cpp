#include <iostream>
#include <cmath>
#include <cstdio>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "gwReadWrite.h"  //Commented out two functions in header file that we don't use but cause errors


//Prototype for the likelihood calculation function
double likelihood(double, double, double, double, double, double (*function)(double,double,double,double,double));
double gaussian(double, double, double, double, double);
double prior(double, double, double, double, double, double);


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

	double nSigma = sigma/10;


	//Takes inputs for the centre coordinates of the likelihood
	std::cout << "Enter the values (a,b) the likelihood is centred on:" 
			  << std::endl;
	double a,b;
	std::cout << "a:"<< std::endl;
	std::cin >> a;
	std::cout << "b:"<< std::endl;
	std::cin >> b;


	//Takes x limits for the prior function
	double xLower,xUpper;
	std::cout<< "Enter the x upper limit:" << std::endl;
	std::cin >> xUpper;
	std::cout<< "Enter the x lower limit:" << std::endl;
	std::cin >> xLower;

	//Takes y limits for the prior function
	double yLower,yUpper;
	std::cout<< "Enter the y upper limit:" << std::endl;
	std::cin >> yUpper;
	std::cout<< "Enter the y lower limit:" << std::endl;
	std::cin >> yLower;

	
	//Takes an input for the number of samples used in the Monte Carlo routine
	std::cout<< "Enter the number of Monte-Carlo samples:" << std::endl;
	int N;
	std::cin >> N;


	//Declares the variables used throughout the routine
	double x, xPrime, y, yPrime, p, pPrime, alpha;
	double nZeroX, nZeroY, rZero;


	//Opens the output text file
	FILE * outFile;
	outFile = fopen("2DMonte.txt","w");
	fprintf(outFile,"%.15g,%.15g,%.15g\n",sigma,a,b);


	//std::vector<Signal>* data;

	//loadSignals(" ", data, tab);


	//Sets the starting x and y values for the routine as random integers in the range [-50,50]
	x = gsl_rng_uniform(startGen)*(xUpper-xLower)+xLower;
	y = gsl_rng_uniform(startGen)*(yUpper-yLower)+yLower;


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
		

		p = likelihood (x,y,a,b,sigma,gaussian)*prior(x,y,xUpper,xLower,yUpper,yLower);

		nZeroX = gsl_ran_gaussian(normGen, nSigma);
		nZeroY = gsl_ran_gaussian(normGen, nSigma);

		xPrime = x + nZeroX;
		yPrime = y + nZeroY;
		pPrime = likelihood(xPrime,yPrime,a,b,sigma,gaussian)*prior(xPrime,yPrime,xUpper,xLower,yUpper,yLower);

		alpha = pPrime/p;

		rZero = gsl_rng_uniform_pos(rGen);

		if(alpha > rZero) {
		
			x = xPrime;
			y = yPrime;

		}
		
		if (i%100==0){
			fprintf(outFile,"%.15g,%.15g\n",x,y);
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
double likelihood(double x, double y, double a, double b, double s, double (*function)(double,double,double,double,double)) {

	double g;
    g = (*function)(x,y,a,b,s);
    return (g);

}

/*Defines an example of a normalised likelihood function, which is a Gaussian of (x,y) centred on (a,b)
/ and with a standard deviation of s.*/
double gaussian(double x, double y, double a, double b, double s) {

	return 1/(2*M_PI*s*s)*exp(-(pow((x-a),2) + pow((y-b),2))/(2*s*s));

}


/*Defines a prior function which is a step of height 1, centred on the origin,
/ with a width of w.*/
double prior(double x, double y, double xUpper, double xLower, double yUpper, double yLower) {

	if((x < xUpper && x > xLower) 
	&& (y < yUpper && y > yLower)) { return 1/((xUpper-xLower)*(yUpper-yLower)); }

	else { return 0; }

}
