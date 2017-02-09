#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <iomanip>

#define C_PI 3.1415926536

using namespace std;

void sigDriver( double (*generator)(), double sigLength, int samples);
double noiseGenGaussian();

int main(){
	
	sigDriver(noiseGenGaussian, 100, 1000);
	
	return 0;
	
}

void sigDriver( double (*generator)(), double sigLength, int samples){

	double t=0;
	double dt = sigLength/samples;
	double time[samples];
	double amp[samples];
	
	ofstream sigFile;
	sigFile.open("signal.csv");
	
	for(int i = 0; i<samples; i++){
		
		amp[i] += (*generator)();
		
		time[i] = t;
		
		t+=dt;
		
		sigFile<<time[i]<<","<<amp[i]<<"\n";
		
	}

}

/*
		LITERALLY COPIED FROM SOMEBODY ONLINE
		
		MUST BE CHANGED SO AS NOT TO VIOLATE UNIVERSITY/HUMAN RIGHTS LAW
		
		DON'T LET IT BE SEEN BY ANY MEMBERS OF STAFF
*/
double noiseGenGaussian(){/* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */
 
  double temp1;
  double temp2;
  double result;
  int p;

  p = 1;

  while( p > 0 )
  {
	temp2 = ( rand() / ( (double)RAND_MAX ) ); /*  rand() function generates an
                                                       integer between 0 and  RAND_MAX,
                                                       which is defined in stdlib.h.
                                                   */

    if ( temp2 == 0 )
    {// temp2 is >= (RAND_MAX / 2)
      p = 1;
    }// end if
    else
    {// temp2 is < (RAND_MAX / 2)
       p = -1;
    }// end else

  }// end while()

  temp1 = cos( ( 2.0 * (double)C_PI ) * rand() / ( (double)RAND_MAX ) );
  result = sqrt( -2.0 * log( temp2 ) ) * temp1;

  return result;	// return the generated random sample to the caller

}// end AWGN_generator()
/*
		CONSIDER YOURSELF WARNED
*/
