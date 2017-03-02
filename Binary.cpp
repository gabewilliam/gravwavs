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
    double rl = 0.49*qtwothird/((0.6*qtwothird)+log(1+qonethird));
    return rl;
    
}

double Binary::keplerFrequency(){
	double b=pow(fa,3);
	double freqsquared=(G(fm1+fm2))/b;
	double freq=pow(freqsquared,1/2);
	return freq;
}

double Binary::mixingFrequency(int n){
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
		omegaC=0.2+(2.7e-4*((m/mSolar)-50)*((m/mSolar)-50))
	}
	else if (m>=50*mSolar){
		omegaC=0.2;
	}
	return omegaC;
}

double Binary::mergeTime(){
	cfive=pow(c,5);
	Gthree=pow(G,3);
	b=pow(fa,4)
	return Tmerge=(5/256)*(cfive/(Gthree*fm1*fm2*(fm1+fm2)))*b;

}
