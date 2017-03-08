//Write functions for Binary.h in here
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
   
  this->updateRadii();//this is my favourite line of code
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

void Binary::printGets(){//easy to check information on a binary (could probably change to access datamembers directly)
	
	//std::cout<<"printing all gets"<<std::endl;

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
        fm1=m;
    }
    
    else if(n==2) {
        fm2=m;
    }
    
    else {
        std::cout<<"you fucked it (setMass)"<<std::endl;
    }
    
	this->updateRadii();
    return;    

}

void Binary::setSeparation(double a){
    
    fa=a;
    return;    

}



//COMPUTING BINARY PROPERTIES

void Binary::updateRadii() {
    
    fr1=pow(fm1,0.6);
    fr2=pow(fm2,0.6);
    
}

double Binary::rocheLobe(int n) {
    
    double q = this->getRatio(n);
    double qOneThird = pow(q,1.0/3.0);
    double qTwoThird = pow(q,2.0/3.0);
	double a = this->getSeparation()*(AU/rSolar); //Finds a in rSolar

    double rl = 0.49*qTwoThird*a/((0.6*qTwoThird)+log(1+qOneThird));
    return rl;//In units of solar radii
    
}

double Binary::keplerFrequency() {

	double a = this->getSeparation()*AU;//Gets a in m
	double m1 = this->getMass(1)*mSolar;
	double m2 = this->getMass(2)*mSolar;//Gets both masses in kg
	double b = pow(a,3);
	double freqSquared = ((m1+m2))/(b);//removed factor of G
	double freq = pow(freqSquared,0.5);
	//std::cout<<std::endl<<freq;
	return freq;

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
		std::cout<<"You fucked it (mixingFrequency)"<<std::endl;
		return 0;
	}

	if(m<50.0) {
		mr=0.2+((2.7e-4)*(m-50)*(m-50));
	}
	else if(m>=50.0) {
		mr=0.2;
	}

	//std::cout<<"  ,  "<<mr<<std::endl;

	return mr;

}

double Binary::omegaC (int n){

	double omegaC = this->keplerFrequency();
	double r;
	
	if (n==1){
		r=pow(fr1,1.5)/pow(fm1,0.5);
	}
	if (n==2){
		r=pow(fr2,1.5)/pow(fm2,0.5);
	}

	omegaC = omegaC*r;
	return omegaC;

}

double Binary::mergeTime() {

	double cFive = pow(c,5);
	double GThree = pow(G,3);
	double a = fa*AU;//Gets the separation in m
	double m1 = fm1*mSolar;
	double m2 = fm2*mSolar;
	double b = pow(a,4);
	double Tmerge = (5.0/256.0)*(cFive/(GThree*m1*m2*(m1+m2)))*b;//In s 
	return Tmerge;//In seconds

}



//EVOLUTION OF STARS

void Binary::evolveMainSequence () {

	double dm = fm1*0.1+fm2*0.1;
	fa = fa*(fm1+fm2)/(fm1+fm2-dm);
	fm1 = 0.9*fm1;
	fm2 = 0.9*fm2;
	return;

}

void Binary::evolveWolfRayet() {

	double dm = fm1*0.25+fm2*0.25;
	fa = fa*(fm1+fm2)/(fm1+fm2-dm);
	fm1 = 0.75*fm1;
	fm2 = 0.75*fm2;
	return;	

}

void Binary::evolveSupernova() {

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
	double mMax = 300;
	double mMin = 20;

	double m1 = this->getMass(1);
	double m2 = this->getMass(2);

	if((m1 >= mMin) && (m1 <= mMax) && (m2 >= mMin) && (m2 <=mMax)) {
		return true;
	}
	else {
		return false;
	}

}

bool Binary::checkRocheLobe() {
	
	double r1 = this->getRadius(1);
	double r2 = this->getRadius(2);
	double rl1 = this->rocheLobe(1);
	double rl2 = this->rocheLobe(2);
	if(r1<=rl1 && r2<=rl2) {
		return true;
	}
	else {
		return false;
	}
	
}

bool Binary::checkHomogeneousMixing() {

	double mr1 = this->mixingRatio(1);
	double mr2 = this->mixingRatio(2);
	double omegaC1 = this->omegaC(1);
	double omegaC2 = this->omegaC(2);
	if(omegaC1 >= mr1 && omegaC2 >= mr2) {
		return true;
	}
	else {
		return false;
	}

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

	double tm = this->mergeTime();
	double tH =  4.55e17;// s, check using worksheet (this is the correct value for hubble time)
	if (tm==0){
		std::cout<<"broken"<<std::endl;
	}
	if(tm < tH) {
		return true;
	}
	else {
		return false;
		std::cout<<"slow merger"<<std::endl;
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
	if(this->checkMergeTime() == false){
		return 6;
	}	
	else {
		return 0;
	}
}
	
	

bool Binary::steveIsATotalLoserAndDeservesToDie(){
	return true;
}

