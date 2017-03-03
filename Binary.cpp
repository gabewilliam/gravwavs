//Write functions for Binary.h in here
#include "Binary.h"
#include <cmath>
#include <iostream>

//const double, can be taken out if declared in main
const double c = 3.0e8; //ms^-2 
const double G = 6.67408e-11; //m^3 kg^-1 s^-2

//constructors

Binary::Binary() {}

Binary::Binary(double mass1, double mass2,  double separation)
	:fm1(mass1), fm2(mass2), fa(separation){

	this->updateRadii();

}

Binary::~Binary(){}


//GET METHODS

double Binary::getMass (int n){
    
    if (n==1){
        return fm1;
    }
    
    else if (n==2){
        return fm2;
    }
    
    else {
        std::cout<<"you fucked it"<<std::endl;
        return 0;
    }
    
}

double Binary::getRadius (int n){    
   
  this->updateRadii();//this is my favourite line of code
  if (n==1){
        return fr1;
    }
    
    else if (n==2){
        return fr2;
    }
    
    else {
        std::cout<<"you fucked it (getRadius)"<<std::endl;
        return 0;
    }
    
}

double Binary::getSeparation (){
    
    return fa;//kinda reads like there's no output
    
}

double Binary::getRatio (int n){
    
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


void Binary::printGets(){//easy to check information on a binary
	
	std::cout<<"printing all gets"<<std::endl;

	for (int n=1; n<3; n++){

		std::cout<<"Star "<<n<<": ";
		std::cout<<"mass = "<<this->getMass(n);
		std::cout<<", radius = "<<this->getRadius(n)<<std::endl;

	}
	std::cout<<"Binary separation = "<<this->getSeparation()<<std::endl;
	std::cout<<"Binary mass ratio (m2/m1) = "<<this->getRatio(2)<<std::endl;
	return;

}

//SET METHODS

void Binary::setMass (int n, double m){
    
    if (n==1){
        fm1=m;
    }
    
    else if (n==2){
        fm2=m;
    }
    
    else {
        std::cout<<"you fucked it (setMass)"<<std::endl;
    }
    
    return;    

}

void Binary::setSeparation (double a){
    
    fa=a;
    return;    

}


//COMPUTING BINARY VARIABLES

void Binary::updateRadii (){
    
    fr1=pow(fm1,0.6);
    fr2=pow(fm2,0.6);
    
}

double Binary::rocheLobe(int n){
    
    double q = this->getRatio(n);
    double qonethird = pow(q,1/3);
    double qtwothird = pow(q,2/3);
    double rl = 0.49*qtwothird*fa/((0.6*qtwothird)+log(1+qonethird));
    return rl;
    
}


//DOUBLES

double Binary::keplerFrequency(){

	double b=pow(fa,3);
	double freqsquared=(G*(fm1+fm2))/b;
	double freq=pow(freqsquared,1/2);
	return freq;

}

double Binary::mixingFrequency(int n){

	double mSolar = 1.989e30;
	double omegaC = 0.0; 
	double m;	
	if 	(n==1){
		m=fm1; 
	}
	else if (n==2){
		m=fm2;
	}
	else {
		std::cout<<"You fucked it (mixingFrequency)"<<std::endl;
		return 0;
	}

	if (m<50*mSolar){
		omegaC=0.2+(2.7e-4*((m/mSolar)-50)*((m/mSolar)-50));
	}
	else if (m>=50*mSolar){
		omegaC=0.2;
	}
	return omegaC;

}

double Binary::mergeTime(){

	double cfive=pow(c,5);
	double Gthree=pow(G,3);
	double b=pow(fa,4);
	double Tmerge=(5/256)*(cfive/(Gthree*fm1*fm2*(fm1+fm2)))*b;
	return Tmerge;

}


//EVOLUTION OF STARS

void Binary::evolveMainSequence (){

	double dm=fm1*0.1+fm2*0.1;
	fa=fa*(fm1+fm2)/(fm1+fm2-dm);
	fm1=0.9*fm1;
	fm2=0.9*fm2;
	return;

}

void Binary::evolveWolfRayet (){

	double dm=fm1*0.25+fm2*0.25;
	fa=fa*(fm1+fm2)/(fm1+fm2-dm);
	fm1=0.75*fm1;
	fm2=0.75*fm2;
	return;	

}

void Binary::evolveSupernova (){

	double dm=fm1*0.1+fm2*0.1;
	fa=(fa*(fm1+fm2-dm))/(2*(fm1+fm2-dm)-(fm1+fm2));
	fm1=0.9*fm1;
	fm2=0.9*fm2;
	return;	

}


//CHECK FUNCTIONS

bool Binary:: checkRocheLobe (){
	
	this->updateRadii();
	double r1 = this->getRadius(1);
	double r2 = this->getRadius(2);
	double rl1 = this->rocheLobe(1);
	double rl2 = this->rocheLobe(2);
	if (r1<=rl1 && r2<=rl2){
		return true;
	}
	else {
		return false;
	}
	
}

bool Binary:: checkHomogeneousMixing (){

	double omegaK = this-> keplerFrequency();
	double omegaM1 = this->mixingFrequency(1);
	double omegaM2 = this->mixingFrequency(2);
	if (omegaK >= omegaM1 && omegaK >= omegaM2) {
		return true;
	}
	else {
		return false;
	}

}

bool Binary:: checkPairInstability (){

	if (fm1 < 63 && fm2 < 63){
		return true;
	}
	else {
		return false;
	}

}

bool Binary::checkMergeTime (){

	double tm = this->checkMergeTime();
	double tH =  4.55e17;// s, check using worksheet
	if (tm<tH){
		return true;
	}
	else {
		return false;
	}
}


