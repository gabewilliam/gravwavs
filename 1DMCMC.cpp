#include <iostream>
#include <cmath>
#include <cstdio>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double likelihood(double x, double sigma);

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


	double x, xPrime, p, pPrime, alpha;
	double nZero, rZero;


	FILE * outFile;
	outFile = fopen("1DMonte.txt","w");


	x = 0;

	for(int i = 1; i <= N; i++) {

		p = likelihood(x,sigma);


		nZero = gsl_ran_gaussian(normGen, nSigma);

		xPrime = x +nZero;
		pPrime = likelihood(xPrime,sigma);

		alpha = pPrime/p;

		rZero = gsl_rng_uniform_pos(rGen);

		if(alpha > rZero) {
		
			x = xPrime;

		}

		fprintf(outFile,"%20.15g\n",x);

	}

	fclose(outFile);

	gsl_rng_free(normGen);
	gsl_rng_free(rGen);
	

}

double likelihood(double x, double sigma) {

	return exp(-(x*x)/(2*sigma*sigma));

}
