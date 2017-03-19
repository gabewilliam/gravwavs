#include "Bin.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

//CONSTRUCTORS
Bin::Bin() {}

Bin::Bin(double redshift):fredshift(redshift){
fmergeRateSum=0;
fcounter=0;
ftLookback=tLookback(fredshift);
}

//DESTRUCTOR
Bin::~Bin() {}

double Bin::getRedshift(){

	return fredshift;

}

double Bin::gettLookback(){

	return ftLookback;

}

double Bin::getMergeRateSum(){

	return fmergeRateSum;

}

double Bin::getCounter(){

	return fcounter;

}


double Bin::CDF(double Z, double z){
	const double ZSolar = 0.0134;
	double alpha = -1.16;
	double beta = 2;

	double a = alpha + 2;
	double Zbeta = pow((Z/ZSolar),beta);
	double tenbetaz=pow(10,(0.15*beta*z));
	double x=Zbeta*tenbetaz; 

	return gsl_sf_gamma_inc_P(a,x);
}

double Bin::dMSFRdtdV(double z){

	double num=pow((1+z), 2.7);
	double denom=pow((1+(1+z)/2.9),5.6);
	
	return 0.015*num/denom;

}

double Bin::Ez(double z){
	const double omegaLambda = 0.718;
	const double omegaM = 1-omegaLambda;
	
	double Esquared = omegaM*pow((1+z),3) + omegaLambda;

	return sqrt(Esquared);

}

double Bin::tIntegrand(double z, void * params) {

	return 1/(Ez(z)*(1+z));
}

double Bin::tLookback(double z){

	double tH= 4.55e17; //s
	const double yr=365.25*24*60*60; //s
	tH=tH/yr; //converts the Hubble time to years

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

	double tL, tError;
	double param = 1;
	double epsRel=1e-5;

	gsl_function tInt;
	tInt.function = &tIntegrand;
	tInt.params = &param;

	gsl_integration_qag(&tInt, 0, z, 0, epsRel, 1000, 6, w, &tL, &tError);

	gsl_integration_workspace_free(w);

	return tH*tL;
}

void Bin::addToMergeRateSum(double birthrate){

	fmergeRateSum=fmergeRateSum + birthrate;

}

void Bin::addToCounter(){

	fcounter=fcounter + 1;

}

