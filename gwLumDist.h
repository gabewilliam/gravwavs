#include <cmath>
#include <iomanip>

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
	
	
	
	for(int i = 1; i*dz<=z; i++){
		fn1 = pow((Om*pow(1+(i*dz),3))+Oa,-0.5);
		I += (fn1+fn0)*dz/2;
		fn0 = fn1;
	}
	
	return c*(1+z)*I/H0;
}

double rs(double lumDis){
	
	double
	rs0 = 0, //redshift lower limit
	rsN = 40, //redshift upper limit (size of observal universe)
	drs = 0.001; //initial variation in redshift
	
	const size_t N = (int)((rsN-rs0)/drs) + 1;
	
	double * z = new double[N];
	double * lD = new double[N];
	
	for(int i = 0; i<N; i++){
		z[i] = rs0 + i*drs;
		lD[i] = Dl(z[i]);
	}

	const gsl_interp_type *t = gsl_interp_linear;
	gsl_interp *lin=gsl_interp_alloc(t, N);
	gsl_interp_accel* acc=gsl_interp_accel_alloc();
	gsl_interp_init(lin, lD,z,N );
	double xi=4.5;
	double yi=gsl_interp_eval(lin, &lD[0], &z[0], lumDis, acc);
	   
	delete [] z;
	delete [] lD;
	
	return yi;
	
}