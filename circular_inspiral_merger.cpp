#include <cmath>
#include <iostream>

using namespace std;

const double c_G = 6.674E-11;
const double c_c = 2.998E8;
const double c4 = pow(c_c, 4);
const double G2 = pow(c_G, 2);
const double c5 = pow(c_c, 5);
const double G3 = pow(c_G, 3);

struct twoVect{
	double x, y;
};

int main(){
	
	return 0;
}

double lambda(m1, m2){
	
	double M = m1*m1*m2*m2*(m1+m2);

	return(G3*M)/(c5);
	
}

double dadt(double a, double m1, double m2){
	
	double a3 = pow(a,3);
	
	return -(64/5)*lambda(m1,m2)*(1/a3);
	
}

double timeToMerge(double a, double m1, double m2){

	double a4 = pow(a,4);
	
	return (5/256)*(a4/lambda(m1,m2));
	
}

double sepAtTime(double a0, double m1, double m2, double t0, double dt){
	
	double a04 = pow(a0,4);
	double lamb = lamda(m1,m2);
	
	double a4 = (a04 - (256/5)*lamb*dt);
	
	return pow(a4, 0.25);
	
}

double hPlus(double R, double th, double om, double t, double M, double r){
	
	return -(G2*2*M*(1+cos(th)*cos(th))*cos(2*om*(t-(R/c_c))))/(R*c4*r);
	
}

double hCross(double R, double th, double om, double t, double M, double r){
	
	return -(G2*4*M*(cos(th))*sin(2*om*(t-(R/c_c))))/(R*c4*r);
	
}

twoVect hTot(double R, double th, double om, double t, double m1, double m2, double r){
	
	double M = m1*m2
	double x = hPlus(R, th, om, t, M, r);
	double y = hCross(R, th, om, t, M, r);
	
	twoVect plrn;
	plrn.x = x;
	plrn.y = y;
	
	return plrn;
	
}

double circInspiral(double m1, double m2, double a1, double a2){
	
}