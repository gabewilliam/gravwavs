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
  double rl = 0.49*qtwothird/((0.6*qtwothird)+log(1+qonethird));
  return rl;
  
}
  
