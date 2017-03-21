#include "pe_gwLikelihood.h"
#include "gwSigGen.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

void saveToFile(vec_d,vec_d,vec_d,int,int,std::string);


//AIM: reads in d(t) from signal, generates h(t) from known formula -> fourier transforms both to get d(f) and h(f), s(f) is hopefully already known and in frequency domain.
//Uses d(f), h(f) and s(f) using equation A20 in "Veich, Vecchio (2010)" to calculate likelihood using two parameters m1 and m2.
//Returns a NON-NORMALISED likelihood.

long double likelihood( double m1, double m2, double d, std::string signalFile, Signal dataSignal, vec_d noiseAmplitude){

	vec_d modelAmplitude = ParameterFunction( m1, m2, d, dataSignal.waveform[0] );//Creates model function for given masses.

	int n = dataSignal.waveform[0].size();//finds number of elements in freq domain signal
	//std::cout<<n<<"\t"<<modelAmplitude.size()<<std::endl;
	long double sum = 0.0;
	long double SNR = 0.0;
	vec_d dataAmplitude = dataSignal.waveform[1];//splits signal into std::vectors for use in loop	

	for( int i = 0; i < n; i++ ){ //evaluates the sum part in equation A20
		//std::cout<<noiseAmplitude[i]<<std::endl;
		//std::cout<<dataAmplitude[50]<<"\t"<<noiseAmplitude[i]<<std::endl;
		sum += pow( std::abs( dataAmplitude[i] - modelAmplitude[i] ), 2 )/noiseAmplitude[i];
		//std::cout<<modelAmplitude[i]<<std::endl;
		SNR += pow( std::abs( dataAmplitude[i] ), 2 )/noiseAmplitude[i];
		//std::cout<<SNR<<std::endl;
	}

	double T = 1/(dataSignal.waveform[0][1]-dataSignal.waveform[0][0]);//total time elapsed
	//std::cout<<T<<std::endl;
	SNR = sqrt(SNR/T);	
	//std::cout<<SNR<<std::endl;
	
	long double pdh = (-2/T) * sum;//calculates p(d|h)

	return pdh;

}

vec_d ParameterFunction( double m1, double m2, double d, vec_d f ){ 

	int n = f.size();
	vec_d modelAmplitude;	

	parameters PARAMS;
	parameters *P = &PARAMS;

	setFundamentalParameters(0.0,3.0,0.0,10.0,m1,m2,d,0.0,0.0,0.0,P);
	
	complex<double> GRAVWAV = 0.0;
	complex<double> *G = &GRAVWAV;

	for(int i=0; i<n; i++){
		modelAmplitude.push_back(updatedAmplitude(P,f[i]/C_CONST,G));
		//std::cout<<updatedAmplitude(P,f[i]/C_CONST,G)<<std::endl;
	}
	
	return modelAmplitude;
}


void saveToFile(vec_d parameterA, vec_d parameterB, vec_d parameterC, int lag, int size, std::string fileName) {
	
	//Opens the output text file
	FILE * outFile;
	outFile = fopen(fileName.c_str(),"w");

	for(int i = 0; i < size; i++){
		if (i%lag==0 && i!=0){
			fprintf(outFile,"%.15g,%.15g,%.15g\n",parameterA[i],parameterB[i],parameterC[i]);
		}
	}
	
	//Closes the output file
	fclose(outFile);
}

