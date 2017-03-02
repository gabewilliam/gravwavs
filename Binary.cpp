//Write functions for Binary.h in here
#include "binary.h";


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

double Binary::getRatio (){
  
  double q=fm2/fm1;//assume you want the ratio m2:m1 so that it's <1 but idk
  return q;
  
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

void Binary::setSeparation (double ){
  
  fa=a;
  return;  

}
