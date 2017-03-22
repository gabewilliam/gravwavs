#include "Bin.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>



//CONSTRUCTORS

Bin::Bin() {}

Bin::Bin(double redshift,double (*func)(double,void*)):fRedshift(redshift){

	//Resets the merge rate sum and counter
	fMergeRateSum=0;
	fCounter=0;

	//Finds the lookback time
	ftLookback=Integrator(fRedshift,func);

}



//DESTRUCTOR

Bin::~Bin() {}



//GET METHODS

double Bin::getRedshift(){

	return fRedshift;

}

double Bin::gettLookback(){

	return ftLookback;

}

double Bin::getMergeRateSum(){

	return fMergeRateSum;

}

double Bin::getCounter(){

	return fCounter;

}



//CALCULATIONS

double Bin::CDF(double Z, double z){

	//Defines constants
	const double ZSolar = 0.0134;
	double alpha = -1.16;
	double beta = 2;

	//Finds the variables which are input into the gamma
	double a = alpha + 2;
	double Zbeta = pow((Z/ZSolar),beta);
	double tenbetaz=pow(10,(0.15*beta*z));
	double x=Zbeta*tenbetaz; 

	//Finds the CDF
	return gsl_sf_gamma_inc_P(a,x);

}

double Bin::dMSFRdtdV(double z){

	//Finds the numerator and denominator
	double num=pow((1+z), 2.7);
	double denom=pow((1+(1+z)/2.9),5.6);
	
	return 0.015*num/denom*1e9; //Returns in Gpc^_3yr^-1

}

double Bin::Integrator(double z,double (*func)(double,void*)){
	
	//Creates an integration workspace
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

	//Creates variables for the return values
	double ans, tError;

	/*Defines a dummy paramter and the relative error which is to be
	/calculated to.*/
	double param = 1;
	double epsRel=1e-5;

	//Takes the input function and turns it into a GSL function
	gsl_function fInt;
	fInt.function = func;
	fInt.params = &param;

	//Integrates from 0 to z
	gsl_integration_qag(&fInt, 0, z, 0, epsRel, 1000, 6, w, &ans, &tError);

	//Frees the workspace memory
	gsl_integration_workspace_free(w);

	return ans;

}

void Bin::addToMergeRateSum(double birthRate){
	
	//Adds an input value to the merge rate sum
	fMergeRateSum=fMergeRateSum + birthRate;

}

void Bin::addToCounter(){

	fCounter=fCounter + 1;

}

