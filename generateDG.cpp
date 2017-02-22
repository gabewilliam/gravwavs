#include "generateDG.h"

#include<cmath>

double generate(double ma, double mb, double t, double r=140e6*3.0857e16, double theta=0, double phi=0, double Tm=1){

	static const double G = 6.67408e-11;
	static const double c = 2.99792458e8;
	ma = ma * G/pow(c,2);
	mb = mb * G/pow(c,2);
	Tm = c*Tm;
	t = c*t;
		//convert all passed SI values to natural geometric units
	double M = ma + mb;
	double mu = (ma*mb)/M;
		//calculate total and reduced mass
	double z = r/(0.7*9.26e25);
		//calculate redshift (small distance approx)
	double omega = pow((256./5.)*(mu*pow(M,2./3)*(Tm-t)),-3./8);
		//calculate the angular velocity of the orbit (to convert to freq in Hz, omega*c/(4*atan))
	double hp = -2*(1+pow(cos(theta),2))*(mu*pow(M*omega,2./3.)/r)*cos(2*(omega*t/(1+z) - phi));
		//calculate the component of the wave polarized with the axis
	double hc = -4*cos(theta)*(mu*pow(M*omega,2./3.)/r)*cos(2*(omega*t/(1+z) - phi));
		//calculate the component of the wave polarized diagonal to the axis
	
	double out = hp+hc;
		//compute total polarization of the wave
	
	if (out != out)	return 0;
		//currently code doesnt compute post merger wave, calculating NaN, so return 0 in such a case (for plotting)
	else return out;
	
}
