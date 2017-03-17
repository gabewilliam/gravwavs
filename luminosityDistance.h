#include <cmath>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

//-- Definitions for error warnings
#define SUCCESS 0
#define FAILURE 1

double Dl(double z){
	
	double
	dz = 0.001, //bin width
	I = 0, //the integral
	fn1 = 0, //f(x_{n})
	fn0 = 0, //f(x_{n-1})
		//variables for calculating function value at both ends of trapezium
	H0 = 67.8e3/3.0857e22, //hubble constant in units of per second
	Om = 0.31, //current fractional energy density of matter
	Oa = 1 - Om, //current fractional energy density of dark matter
	c = 3e8; //speed on light in m/s
	
	
	
	for(unsigned int i = 1; i*dz<=z; i++){
		fn1 = pow((Om*pow(1+(i*dz),3))+Oa,-0.5);
		I += (fn1+fn0)*dz/2;
		fn0 = fn1;
	}
	
	return c*(1+z)*I/H0;
}

double rs(double lumDis){
	
	double
	rs0 = 0., //redshift lower limit
	rsN = 40., //redshift upper limit (size of observal universe)
	drs = 0.01; //initial variation in redshift
	
	const size_t N = (int)((rsN-rs0)/drs) + 1;
	
	double * z = new double[N];
	double * lD = new double[N];
	
	std::ifstream inFile;
	inFile.open("lumDis.csv");

	char de = ',';

	if(inFile.fail()){
		std::cout << "lumDis file not found." << std::endl;
/* 		FILE * output;
		output = fopen("lumDis.csv", "w");
	
		for(unsigned int i = 0; i<N; i++){
			z[i] = rs0 + i*drs;
			lD[i] = Dl(z[i]);
		}
		
		for(unsigned int i = 0; i<N; i++){
			if(i!=0)//fprintf(output, ", ", z[i]);
			fprintf(output, "%g", z[i]);
		}
			fprintf(output, "\n");
		
		for(unsigned int i = 0; i<N; i++){
			if(i!=0)//fprintf(output, ", ", lD[i]);
			fprintf(output, "%g", lD[i]);
		}
		
		fclose(output);
		
		std::cout<<"file created"<<std::endl;*/
		
		return FAILURE;

	}
	
	else{
		double d;
		std::string line, element;
		
		int i = 0;

		while(getline(inFile, line))
		{
			int j = 0;
			
			std::istringstream iss(line);

			//looping over the next two lines to extract the time and waveform data
			while(getline(iss, element, de))
			{
				std::istringstream is(element);
				is >> d;
				if(i==0) z[j] = d;
				if(i==1) lD[j] = d;
				j++;
			}
			i++;
		}
	}
		

	const gsl_interp_type *t = gsl_interp_linear;
	gsl_interp *lin=gsl_interp_alloc(t, N);
	gsl_interp_accel* acc=gsl_interp_accel_alloc();
	gsl_interp_init(lin, lD,z,N );
	double yi=gsl_interp_eval(lin, &lD[0], &z[0], lumDis, acc);
   
	delete [] z;
	delete [] lD;
	
	return yi;
	
}