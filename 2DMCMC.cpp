#include <iostream>
#include <cmath>
#include <cstdio>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double likelihood(double x, double y, double sigma);

int main() {

	gsl_rng * normGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng * rGen = gsl_rng_alloc(gsl_rng_taus);
	

	std::cout<< "Enter the value of sigma:" << std::endl;
	double sigma;
	std::cin >> sigma;

	double nSigma = sigma/10;

	std::cout<< "Enter the number of Monte-Carlo samples:" << std::endl;
	double N;
	std::cin >> N;


	double x, xPrime, y, yPrime, p, pPrime, alpha;
	double nZeroX, nZeroY, rZero;


	FILE * outFile;
	outFile = fopen("2DMonte.txt","w");


	x = 0;
	y = 0;

	for(int i = 1; i <= N; i++) {

		p = likelihood(x,y,sigma);


		nZeroX = gsl_ran_gaussian(normGen, nSigma);
		nZeroY = gsl_ran_gaussian(normGen, nSigma);

		xPrime = x + nZeroX;
		yPrime = y + nZeroY;
		pPrime = likelihood(xPrime,yPrime,sigma);

		alpha = pPrime/p;

		rZero = gsl_rng_uniform_pos(rGen);

		if(alpha > rZero) {
		
			x = xPrime;
			y = yPrime;

		}

		fprintf(outFile,"%25.15g%25.15g\n",x,y);

	}

	fclose(outFile);

	gsl_rng_free(normGen);
	gsl_rng_free(rGen);
	

}

double likelihood(double x, double y, double sigma) {

	return exp(-(x*x + y*y)/(2*sigma*sigma));

}
