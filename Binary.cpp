//Write functions for Binary.h in here
#include "binary.h";
#include <math.h>;


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

    this->updateRadii();     
    if (n==1){
        return fr1;
    }
    
    else if (n==2){
        return fr2;
    }
    
    else {
        std::cout<<"you fucked it"<<std::endl;
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
        std::cout<<"you fucked it"<<std::endl;
        return 0;
    }
    
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
        std::cout<<"you fucked it"<<std::endl;
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
    double a = this->getSeparation();
    double rl = 0.49*qtwothird*a/((0.6*qtwothird)+log(1+qonethird));
    return rl;
    
}

//EVOLUTION OF STARS

void Binary::evolveMainSequence (){
	double dm=fm1*0.1+fm2*0.1;
	a=a*(fm1+fm2)/(fm1+fm2-dm);
	fm1=0.9*fm1;
	fm2=0.9*fm2;
	return;
}

void Binary::evolveWolfRayet (){
	double dm=fm1*0.25+fm2*0.25;
	a=a*(fm1+fm2)/(fm1+fm2-dm);
	fm1=0.75*fm1;
	fm2=0.75*fm2;
	return;	
}

void Binary::evolveSupernova (){
	double dm=fm1*0.1+fm2*0.1;
	a=(a*(fm1+fm2-dm))/(2*(fm1+fm2-dm)-(fm1+fm2))
	fm1=0.9*fm1;
	fm2=0.9*fm2;
	return;	
}

//CHECK FUNCTIONS

bool Binary:: checkRocheLobe (int n){
	
	this->updateRadii();
	double r = this->getRadius(n);
	double rl = this->rocheLobe(n);
	if (r<=rl){
		return true;
	}
	else if (r>rl){
		return false;
	}
	
}

bool Binary:: checkHomogeneousMixing (int n){

	double omegaK = this-> keplerFrequency();
	double omegaM = this->mixingFrequency(int n);
	if (omegaK >= omegaM) {
		return true;
	}
	else {
		return false
	}

}

bool Binary:: checkPairInstability (){

	double m1 = this->getMass(1);
	double m2 = this->getMass(2);
	if (m1<63 && m2<63){
		return true;
	}
	else {
		return false;
	}

}

bool binary::checkMergeTime (){

	double tm = this->checkMergeTime();
	tH =  4.55e17;// s, check using worksheet
	if (tm<tH){
		return true
	}
	else {
		return false;
	}
}

	
 
