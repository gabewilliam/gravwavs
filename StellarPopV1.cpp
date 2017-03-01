#include <iostream>
#include <cstdio>
#include <cmath>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


int main() {


	/*Sets up a random number generator to seed the other two generators, and
	/ seeds this based on the system time.*/
	gsl_rng * seedGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(seedGen,(long unsigned int) time(NULL));
			

	/*Initialises random number generators to be used in sampling the masses
	/ of the binary memebers.*/
	gsl_rng * m1Gen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(m1Gen, gsl_rng_get(seedGen));
	gsl_rng * qGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(qGen, gsl_rng_get(seedGen));
	

	/*Defines the exponent of the m1 distribution, and a constant used in
	/ sampling from it.*/
	double a = 1.35;
	double b = 1;


	//Defines the upper and lower limits on mass ratio
	double qMax = 1;
	double qMin = 0.1;


	//Defines the upper and lower limits on stellar mass
	double mMax = 300;
	double mMin = 20;

	
	//Declares variables for the masses and mass ratio.
	double m1, m2, q;


	FILE * popFile1;
	popFile1 = fopen("fullpop.csv","w");
	FILE * popFile2;
	popFile2 = fopen("masslimitpop.csv","w");

	
	//Asks the user how many systems they want to generate	
	int N;
	std::cout << "How many binaries would you like to generate?" <<std::endl;
	std::cin >> N;
	

	/*Generates N binaries with component masses governed by the expressions
	/ cited in the Mandel and deMink paper. First m1 is sampled. Then the
	/ mass ratio is sampled and used to find m2. All of the generated systems
	/ are output to one file. Those satisfying the 20 < m < 300 component
	/ mass limits are output to a second file.*/
	for(int n = 0; n < N; n++) {


		m1 = gsl_ran_pareto(m1Gen,a,b);
		
		q = gsl_ran_flat(qGen,qMin,qMax);
		m2 = m1*q;


		fprintf(popFile1,"%.15g,%.15g\n",m1,m2);

		if((m1 >= mMin) && (m2 >= mMin)
			&& (m1 <= mMax) && (m2 <= mMax)) {

			fprintf(popFile2,"%.15g,%.15g\n",m1,m2);

		}


	}	

	
	fclose(popFile1);
	fclose(popFile2);
	

	//Frees the memory associated with the RNGs
	gsl_rng_free(seedGen);
	gsl_rng_free(m1Gen);
	gsl_rng_free(qGen);


}
