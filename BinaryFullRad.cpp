#include "Binary.h"

#include <cmath>
#include <iostream>



//CONSTANTS
const double c = 3.0e8; //ms^-2 
const double G = 6.6748e-11; //m^3 kg^-1 s^-2
const double Gs = 1.2e19; //m^3 Ms^-1 s^-2
const double mSolar = 1.989e30; //kg
const double rSolar = 6.975e8; //m
const double AU = 1.496e11; //m
const double pi = 3.14;



//CONSTRUCTORS

Binary::Binary() {}

Binary::Binary(double mass1, double mass2,  double separation)
	:fm1(mass1), fm2(mass2), fa(separation){

	//Sets the radii based on the masses
	this->updateRadii();
		
	if (fm2>fm1){
		std::cout<<"are you sure about those masses?"<<std::endl;
	}

}



//DESTRUCTOR
Binary::~Binary() {}



//GET METHODS

double Binary::getMass(int n) {
    
    if (n==1){
        return fm1;//In solar masses
    }
    
    else if (n==2){
        return fm2;//In solar masses
    }
    
    else {
        std::cout<<"you fucked it (getMass)"<<std::endl;
        return 0;
    }
    
}

double Binary::getRadius(int n) {    
   
  this->updateRadii();//this is my favourite line of code <3
  if (n==1){
        return fr1;//In solar radii
    }
    
    else if (n==2){
        return fr2;//In solar radii
    }
    
    else {
        std::cout<<"you fucked it (getRadius)"<<std::endl;
        return 0;
    }
    
}

double Binary::getSeparation(){
    
    return fa;//In AU//kinda reads like there's no output
    
}

double Binary::getRatio(int n) {
    
    if (n==1){
        double q=fm1/fm2;
        return q;
    }
    
    else if (n==2){
        double q=fm2/fm1;
        return q;
    }
    
    else {
        std::cout<<"you fucked it (getRatio)"<<std::endl;
        return 0;
    }
    
}

double Binary::gettMerge(){

	fTMerge = this->mergeTime();//s
	return fTMerge;//s

}

void Binary::printGets(){
	
	//Prints the data memebers of the binary's variables

	for(int n=1; n<3; n++) {

		std::cout<<"Star "<<n<<": ";
		std::cout<<"mass = "<<this->getMass(n);
		std::cout<<", radius = "<<this->getRadius(n)<<std::endl;

	}
	std::cout<<"Binary separation = "<<this->getSeparation()<<std::endl;
	std::cout<<"Binary mass ratio (m2/m1) = "<<this->getRatio(2)<<std::endl;
	std::cout<<"Time to merger = "<<this->mergeTime()<<std::endl;
	return;

}



//SET METHODS

void Binary::setMass(int n, double m){
    
    if(n==1) {
        fm1=m;//mSolar
    }
    
    else if(n==2) {
        fm2=m;//mSolar
    }
    
    else {
        std::cout<<"you fucked it (setMass)"<<std::endl;
    }
    
	//Updates the radii based on the new masses
	this->updateRadii();
    return;    

}

void Binary::setSeparation(double a){
    
    fa=a;//AU
    return;    

}



//COMPUTING BINARY PROPERTIES

void Binary::updateRadii() {
    
	//Finds the radii based on the memebrs masses
	fr1=this->radius(1);
    fr2=this->radius(2);
    
}

double Binary::radius(int n) {
	
	//Defines the constants for the equation
	double theta = 1.1368;
	double iota = 4.0392;
	double kappa = 6.5741;
	double lambda = 0.7869;
	double mu = 0.0826;
	double nu = 0.01077422;
	double xi = 2.0678;
	double o = 10.6115;
	double pie = 0.0029;
	
	double m;

	if(n==1) {
        m=fm1;//mSolar
    }
    
    else if(n==2) {
        m=fm2;//mSolar
    }
    
    else {
        std::cout<<"you fucked it (radius)"<<std::endl;
    }
	
	//Calculates the numerator
	double num = theta*pow(m,2.5) + iota*pow(m,6.5) + kappa*pow(m,11) + 
					lambda*pow(m,19) + mu*pow(m,19.5);
	//Calculates the denominator
	double den = nu + xi*pow(m,2) + o*pow(m,8.5) + pow(m,18.5)
					+pie*pow(m,19.5);

	//Returns the answer
	return num/den;//rSolar
	

	return pow(m,0.6);
}

double Binary::rocheLobe(int n) {
    
	//Finds the appropriate mass ratio and raises it to some powers
    double q = this->getRatio(n);
    double qOneThird = pow(q,1.0/3.0);
    double qTwoThird = pow(q,2.0/3.0);

	//Finds a in rSolar
	double a = this->getSeparation()*(AU/rSolar); 

	//Calculates the roche lobe radius
    double rl = 0.49*qTwoThird*a/((0.6*qTwoThird)+log(1+qOneThird));
    return rl;//rSolar
    
}

double Binary::mixingRatio(int n) {

	double mr = 0.0; 
	double m;

	if(n==1) {
		m = fm1;//In solar masses
	}
	else if(n==2) {
		m = fm2;//In solar masses
	}
	else {
		std::cout<<"You fucked it (mixingRatio)"<<std::endl;
		return 0;
	}

	//Finds the mixing ratio for the star of interest
	if(m<50.0) {
		mr=0.2+((2.7e-4)*(m-50)*(m-50));
	}
	else if(m>=50.0) {
		mr=0.2;
	}

	//Returns the minimum required mixing ratio
	return mr;

}

double Binary::omegaC (int n){
	
	//Checks the radii are right
	this->updateRadii();

	//Gets the masses and radii of the stars
	double m1 = this->getMass(1)*mSolar;
	double m2 = this->getMass(2)*mSolar;//kg
	double r1 = this->getRadius(1)*rSolar;
	double r2 = this->getRadius(2)*rSolar;//m
	double a = this->getSeparation()*AU;//m
	double M = m1 + m2;//kg
	
	double omegaC,omegaCsquared;

	if (n==1){
		omegaCsquared = (M*pow(r1,3))/(m1*pow(a,3));
		omegaC = pow(omegaCsquared,0.5);
	}
	if (n==2){
		omegaCsquared = (M*pow(r2,3))/(m2*pow(a,3));
		omegaC = pow(omegaCsquared,0.5);
	}

	//Returns the actual mixing ratio
	return omegaC;

}

double Binary::mergeTime() {

	double cFive = pow(c,5);//(ms^-1)^5
	double GThree = pow(G,3);//(m^3 kg^-1 s^-2)^3
	double a = fa*AU;//m
	double m1 = fm1*mSolar;//kg
	double m2 = fm2*mSolar;//kg
	double b = pow(a,4);//m^4

	//Finds the merge tme and sets the data member
	fTMerge = (5.0/256.0)*(cFive/(GThree*m1*m2*(m1+m2)))*b;//s 
	return fTMerge;//In seconds

}



//EVOLUTION OF STARS

void Binary::evolveMainSequence () {
	
	//Evolves the masses and separation
	double dm = fm1*0.1+fm2*0.1;
	fa = fa*(fm1+fm2)/(fm1+fm2-dm);
	fm1 = 0.9*fm1;
	fm2 = 0.9*fm2;
	return;

}

void Binary::evolveWolfRayet() {

	//Evolves the masses and separation
	double dm = fm1*0.25+fm2*0.25;
	fa = fa*(fm1+fm2)/(fm1+fm2-dm);
	fm1 = 0.75*fm1;
	fm2 = 0.75*fm2;
	return;	

}

void Binary::evolveSupernova() {

	//Evolves the masses and separation (?)
	double dm = fm1*0.1+fm2*0.1;
	fa = (fa*(fm1+fm2-dm))/(2*(fm1+fm2-dm)-(fm1+fm2));
	fm1 = 0.9*fm1;
	fm2 = 0.9*fm2;
	return;	

}



//CHECK FUNCTIONS (true - still a candidate for CHE,
//				   false - will not undergo CHE, can be deleted)

bool Binary::checkMassRange() {
	
	//Defines the upper and lower limits on stellar mass
	double mMax = 200;//mSolar
	double mMin = 20;//mSolar

	double m1 = this->getMass(1);//mSolar
	double m2 = this->getMass(2);//mSolar

	if((m1 <= mMax) && (m2 >= mMin)) {
		return true;
	}
	else {
		return false;
	}

}

bool Binary::checkRocheLobe() {
	
	this->updateRadii();
	double r1 = this->getRadius(1);//rSolar
	double r2 = this->getRadius(2);//rSolar
	double rl1 = this->rocheLobe(1);//rSolar
	double rl2 = this->rocheLobe(2);//rSolar
	if(r1<=rl1 && r2<=rl2) {
		return true;
	}
	else {
		return false;
	}
	
}

bool Binary::checkHomogeneousMixing() {

	this->updateRadii();

	double mr1 = this->mixingRatio(1);
	double mr2 = this->mixingRatio(2);
	double omegaC1 = this->omegaC(1);
	double omegaC2 = this->omegaC(2);
	if(omegaC1 >= mr1 && omegaC2 >= mr2 && this->checkBreakup()) {
		return true;
	}
	else {
		return false;
	}
	
}

bool Binary::checkBreakup(){
	
	double omegaC1 = this->omegaC(1);
	double omegaC2 = this->omegaC(2);
	if(omegaC1>1 || omegaC2>1) {
		return false;
	}
	return true;
}

bool Binary::checkPairInstability() {

	if(fm1 < 63.0 && fm2 < 63.0) {
		return true;
	}
	else {
		return false;
	}

}

bool Binary::checkMergeTime() {

	double tm = this->mergeTime();//s
	double tH =  4.34e17;// s, check using worksheet (this is the correct value for hubble time)
	if (tm==0){
		std::cout<<"broken"<<std::endl;
	}
	if(tm < tH) {
		return true;
	}
	else {
		return false;
	}

	
}

int Binary::checkCandidate() {

	if(this->checkMassRange() == false) {
		return 1;
	}
	if(this->checkRocheLobe() == false) {
		return 2;
	}
	if(this->checkHomogeneousMixing() == false) {
		return 3;
	}	
	else {
		this->evolveMainSequence();
	}
	if(this->checkHomogeneousMixing() == false) {
		return 4;
	}
	else {
		this->evolveWolfRayet();
	}
	if(this->checkPairInstability() == false) {
		return 5;
	}
	else {
		this->evolveSupernova();
	}
	if(this->checkMergeTime() == false) {
		return 6;
	}	
	else {
		return 0;
	}
}
	
	

