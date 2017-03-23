
#include "Bin.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


double Ez(double z);
double tIntegrand(double z, void * params);
double dLIntegrand(double z, void * params);
double VcIntegrand(double z, void * params);


//Defines a few key constants in SI units
const double mSolar = 1.989e30; //kg
const double rSolar = 6.975e8; //m
const double AU = 1.496e11; //m
const double c = 3.0e8; //ms^-2
const double G = 6.67408e-11;//m^3 kg^-1 s^-2
const double yr = 365.25*24*60*60;//s
const double tH = 4.34e17;//s
const double H0 = 69.7;//kms^-1Mpc^-1
const double pi = 4*atan(1); 


int main() {

	//Probabilities of the three IMF regions
	size_t K = 3;
	const double IMF[3] = {0.3714, 0.4780, 0.1506};
	

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
	gsl_rng * PGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(PGen, gsl_rng_get(seedGen));
	gsl_ran_discrete_t * mAGen = gsl_ran_discrete_preproc(K, IMF);

	
	//Defines the upper and lower limits on mass ratio
	double qMax = 1;
	double qMin = 0.1;

	//Defines the upper limit on separation (and its logarithms)
	double aMax = 100*AU; //m
	double laMax = std::log10(aMax);

	
	/*Declares variables for the masses, mass ratio, radii, separation, 
	/ minimum separation, total mass generated and period.*/
	double m1, m2, q, r1, r2, la, a, aMin, laMin, mTotal, lP, P;
	mTotal = 0;


	FILE * popFile;
	popFile = fopen("pop.csv","w");
	std::vector <Binary> binaries;
			
	
	//Asks the user how many systems they want to generate	
	int N;
	std::cout << "How many binaries would you like to generate?" <<std::endl;
	std::cin >> N;

	int Nall=N;
	int n = 0;
	

	std::cout << "Generating binaries..." << std::endl;

	//Declares two variables to be used in drawing from the IMF
	size_t b;
	double alpha;

	/*Generates N binaries with component masses governed by the expressions
	/ cited in the Mandel and deMink paper. Also assigns period (and thus
	/ separation based on the same paper.*/
	while(n < N) {//1

		//Picks one of the three IMF regions
		b = gsl_ran_discrete(m1Gen, mAGen);
				
		//Asings m1 appropriately based on which region it fals into
		if(b==0) {//2
			alpha = -0.7;
			m1 = 0.01549;//mSolar
		}
		else if(b==1) {//3
			alpha = 0.3;
			m1 = 0.1010;//mSolar
		}
		else if (b==2) {//4
			alpha = 1.3;
			m1 = -1;
			//Draws a mass from the Kroupa.
			while(m1<0.5 || m1>150){

				m1 = gsl_ran_pareto(m1Gen,alpha,0.5);//mSolar

			}
			
		} 
		
		//Finds the mass ratio, and thus m2
		q = gsl_ran_flat(qGen,qMin,qMax);
		m2 = m1*q;//mSolar

		//Finds the stellar radii (in solar units)
		r1 = pow(m1,0.6);//rSolar
		r2 = pow(m2,0.6);//rSolar
		
		//Finds log(P), and uses this to find P and the separation
		lP = -1;
		while(lP < 0.075 || lP >3.5) {
			lP = gsl_ran_pareto(PGen,-0.5,3.5);
		}
		P = pow(10,lP)*24*60*60;//s
		a = pow((G*(m1*mSolar+m2*mSolar)*P*P)/(4*pi*pi),(1.0/3.0))/AU;//AU
		
		/*
		//Alternative sampler, which samples a instead
		//Finds the lower limits on separation and its logarithm
		aMin = (r1*rSolar) + (r2*rSolar);
		laMin = std::log10(aMin);

		//Randomly samples a separation
		la = gsl_ran_flat(aGen,laMin,laMax);
		a = pow(10,la);	
		a = a/AU;//AU
		*/
		//Pushes the generated binary to the vector		
		binaries.push_back(Binary(m1, m2, a));
		n++;

		std::cout << "\r"<< std::setw(9) << std::right
					  << n << " binaries generated." << std::flush;

		//Accumulates the total mass
		mTotal += (m1+m2);//mSolar
				
		
	}
	

	//Prepares for the sorting process
	std::cout << std::endl;
	std::cout << "Sorting binaries..." << std::endl;
	int Nr = binaries.size();
	int check;
	int whatWentWrong[7] = {0, 0, 0, 0, 0, 0, 0};
	double tm;
	double mCandidates = 0;


	//Cycles through the binaries and discards them if they aren't candidates
	for(int i = N-1; i >= 0; i--) {//5

		//Checks each candidate and records whether it was succesful.
		check = binaries[i].checkCandidate();
		whatWentWrong[check]++;
		
		//Gets rid of it, if it failed
		if(check > 0) {//6
			binaries.erase(binaries.begin() + i);
			Nr=Nr-1;
		}
		//Keeps it if it passed (and writes it to file)
		else {//7
			a = binaries[i].getSeparation();
			tm = binaries[i].mergeTime()/yr;
			m1 = binaries[i].getMass(1);
			m2 = binaries[i].getMass(2);
			fprintf(popFile,"%.15g,%.15g,%.15g,%.15g\n", m1, m2, a, tm);
			mCandidates += (m1+m2);
		}

		std::cout << "\r" << std::setw(9) << std::right
				  << Nr << " binaries remaining." << std::flush;
		
	}


	//Gets the number which survived
	double nCandidates = binaries.size();

	if(nCandidates==0) {
		std::cout<< "Sorry there are no binaries left." <<std::endl;
	}
	

	//Finds the "next" one which would merge
	double tmMin;
	int nextMerge;

	for(int i=0; i<nCandidates; i++) {
		
		tm = binaries[i].mergeTime();

		if((tm<tmMin) || (i==0)){
			nextMerge=i;
			tmMin=tm;
		}
	}


	//Prints the next merger
	std::cout << std::endl << std::endl;	
	std::cout << "Next binary merger will be: "<< std::endl;
	binaries[nextMerge].printGets();


	//Prints details about the progress of the routine
	std::cout << std::endl;
	std::cout << "Binaries kept: " << whatWentWrong[0] << std::endl;
	std::cout << "Binaries remaining after each sorting phase:" << std::endl;

	int nRemaining = N;

	for(int i=1; i<7; i++) {

		nRemaining -= whatWentWrong[i];
		std::cout << i << ": " << nRemaining << std::endl;

	}

	std::cout << "Binaries ''lost'' at each sorting phase: " << std::endl;

	for(int i=1; i<7; i++) {
		std::cout << i << ": " << whatWentWrong[i] << std::endl;
	}


	//Prints the average and total masses
	std::cout << std::endl;
	std::cout << "Average Mass of Binaries = " << mTotal/Nall << std::endl;
	std::cout << "Total Mass of Candidates = " << mCandidates << std::endl;
			

	fclose(popFile);
 	FILE * rateFile;
	rateFile = fopen("rate.csv","w");
			

	//Frees the memory associated with the RNGs
	gsl_rng_free(seedGen);
	gsl_rng_free(m1Gen);
	gsl_rng_free(qGen);
	gsl_rng_free(PGen);
	gsl_rng_free(aGen);


	/*Used from above: a vector of binaries which can evolve via CHE for
	/ which Tmerge has been calculated.*/ 


	//Declares the number of redshift bins, and their range
	double binTot = 200.0;
	double minRedshift = 0;
	double maxRedshift = 2;

	
	/*Declares the metalicity limit, and the fraction of stars which are 
	/ in binaries.*/
	double metallicity = 0.004;
	double binaryPercent = 0.7; //Fraction of stars which are in binaries.	


	/*Finds the redshift bin width and declares variables for the birth and
	/ merge time bin widths.*/
	double binWidth = (maxRedshift-minRedshift)/binTot;
	double birthWidth, mergeWidth;


	//Some kind of thing
	bool massiveError;

	
	/*Creating a vector of redshift bins, each of width 1 and at redshifts 
	/ ranging from 0-maxZ.*/
	std::vector <Bin> redshiftBins; 
	double z;


	//Fills the vector of redshift bins
	for (int i=0; i<=binTot; i++){//8

		z = i*binWidth;

		Bin ithRedShiftBin(z,&tIntegrand);

		redshiftBins.push_back(ithRedShiftBin);

	}


	//Declares variables used in the next loop
	double birthRate, desiredMergeLookbackTime, dtiDivdtj;


	/*Performing the calculation. Requires second for-loop accross the same 
	/ thing as all redshift bins must have been created before can do the 
	/ rest.*/
	//Outside, cycles over birth time
	for (int i=0; i<binTot; i++){//9

		//Gets the current birth bin's redshift
		z = redshiftBins[i].getRedshift();
	
		//Step one: calculate the birth rate for binaries in bin i.
		birthRate=(redshiftBins[i].CDF(z,metallicity))*(redshiftBins[i].dMSFRdtdV(z))*binaryPercent/mTotal; //Gpc^_3yr^-1 
		
		//Finds the width of the birth time bin
		birthWidth = redshiftBins[i+1].gettLookback() -
						 redshiftBins[i].gettLookback();//yr

		
		/*Step 2: for each binary, j, calculate the lookback time at the 
		/ point of merging.*/ 
		for (int j=0; j<Nr; j++){//10

			desiredMergeLookbackTime=(redshiftBins[i].gettLookback()) -
											binaries[j].gettMerge()/yr;//yr 
			
			/*Step 3: if the binary has merged, loop through all the 
			/ redshift bins with tLookback less than that of bin i, to find 
			/ the first one where the lookback time is greater than the 
			/ desired look back time. The birth rate is then added to the 
			/ sum of the merge rates in the bin before that with the greater
			/ lookback time than desired.*/
			//Checks the binary merges in the past or present
			if (desiredMergeLookbackTime>=0){//11
								
				//This is fine
				massiveError = true;

				//Loops over all merge time bins
				for (int k=0; k<=i; k++){//12
	
					//Picks out the right binary
					if (desiredMergeLookbackTime<redshiftBins[k].gettLookback() && massiveError == true){//13

						//Finds the width of the merge time bin
						mergeWidth = redshiftBins[k].gettLookback() - 
										redshiftBins[k-1].gettLookback();//yr
						//Takes a ratio of the widths
						dtiDivdtj = birthWidth/mergeWidth;

						//Adds to the sum and tracks which bin we're on
						redshiftBins[k-1].addToMergeRateSum(birthRate*dtiDivdtj);
						redshiftBins[k-1].addToCounter();

						//Makes sure the binary isn't counted again
						massiveError = false;

					}

				}

			}

		}

	}


	//Declares variables for the quantities of interest
	double mergeRate, dL, Vc;

	//Loops back over all of the bins to find quantities of interest
	for (int i=0; i<binTot; i++){//14

		//Finds the total merge rate
		mergeRate=(redshiftBins[i].getMergeRateSum());//Gpc^_3yr^-1
		//Finds the redshift
		z = redshiftBins[i].getRedshift();
		//Finds the luminosity distance in Mpc
		dL=(c*pow(10,-3)*(1+z))*redshiftBins[i].Integrator(z,&dLIntegrand)*(1/H0);
		//Finds the comoving volume in Mpc^3
		Vc=redshiftBins[i].Integrator(z,&VcIntegrand);
		//Otputs the rates to the console			
		std::cout<<"For redshift bin "<< i 
				 <<" the merge rate is "<< mergeRate << std::endl;
		//Outputs everything to file
		fprintf(rateFile,"%.15g,%.15g,%.15g,%.15g,%.15g\n",
				z,dL,Vc,redshiftBins[i].gettLookback(),mergeRate);

	}

	//Declares some variables used to find the cumulative merger rate
	double zMax = 2;
	int binMax = 1+zMax/binWidth;
	double lastVc, dVc, dz;
	void * fudge;
	double cumRate = 0;//rofl
	Vc=0;

	

	//Loops up to the required redshift and accumulates the merger rate
	for(int i=0; i<binMax; i++) {
		lastVc=Vc;
		z = redshiftBins[i+1].getRedshift();

		if(z<zMax){
			mergeRate=(redshiftBins[i+1].getMergeRateSum());//Gpc^{-3}yr^{-1}
			Vc = redshiftBins[i+1].Integrator(z, &VcIntegrand); //Mpc^3
			dVc = Vc-lastVc;
			
			cumRate += mergeRate*(1e-9)*(dVc)/(1+z);//yr^-1
		}
		
	}

	//Outputs the result
	std::cout<<std::endl;
	std::cout << "Merge rate within z = " <<zMax<< " : " 
			  << cumRate << "yr^-1" << std::endl;

		
	//Closes the file
	
	fclose(rateFile);
	std::cout << "Massive Death" <<std::endl;
	
}


//FUNCTION DEFINITIONS: 

double Ez(double z){
	//Defines the density parameters
	const double omegaLambda = 0.718;
	const double omegaM = 1-omegaLambda;
	
	//Finds E^2(z)S
	double Esquared = omegaM*pow((1+z),3) + omegaLambda;

	//Returns E(z)
	return sqrt(Esquared);

}


double tIntegrand(double z, void * params) {
	
	//Returns the integrand of lookback time
	return (tH/yr)/(Ez(z)*(1+z));//yr
}

 
double dLIntegrand(double z, void * params){
	
	//Returns the integrand of luminsoity distance
	return (1/Ez(z));//Mpc
}


double VcIntegrand(double z, void * params){
	/*Creates a dummy bin to fudge the GSL, with its lookback time function
	/ overwritten by the luminosity distance function*/
	Bin dummybin(z,&dLIntegrand);

	//Finds the luminosity distance of the redshift given to the dummy bin
	double dLsquared=pow(((c*pow(10,-3)*(1+z))*dummybin.gettLookback()*(1/H0)),2);	 
	
	//Returns the integrand of comoving volume
	return (4*pi*c*pow(10,-3)*dLsquared)/(H0*(1+z)*Ez(z));//Mpc^3
}




