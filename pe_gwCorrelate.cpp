#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <gsl/gsl_statistics_double.h>

double autoCorrelation(double [], int);
double checkForNaN(double);
void saveToFile(double [],double [],double [],int,int,std::string);

int main(){
	
	std::string fileName = "RawData.txt";
	
	std::ifstream lineFile;
	lineFile.open( fileName.c_str() );

	int N = 0;
	std::string line;
	std::vector< std::string > input;
	
	while( std::getline( lineFile, line )){
		input.push_back( line );
		++N;
	}
	std::cout << N << std::endl;

	FILE * inFile;
	inFile=fopen( fileName.c_str() ,"r" );
	
	double * chirpMass = new double[N];
	double * massRatio = new double[N];
	double * distance = new double[N];

	for( int i = 0; i < N; i++ ){
		fscanf( inFile, "%lf %*c %lf %*c %lf\n", &chirpMass[i], &massRatio[i], &distance[i] );
	}

	double corrFactorChirp = checkForNaN( autoCorrelation( chirpMass, N ) );
	double corrFactorRatio = checkForNaN( autoCorrelation( massRatio, N ) );
	double corrFactorDistance = checkForNaN( autoCorrelation( distance, N ) );
	
	std::cout << "Chirp Mass: " << corrFactorChirp << std::endl;
	std::cout << "Mass Ratio: " << corrFactorRatio << std::endl;
	std::cout << "Distance: " << corrFactorDistance << std::endl;

	
	double acls [] = {corrFactorChirp, corrFactorRatio, corrFactorDistance};
	int aclMax = *(std::max_element(acls,acls+3));

	saveToFile( chirpMass, massRatio, distance, aclMax, N, "UncorrelatedData.txt" );

	delete [] chirpMass;
	delete [] massRatio;	
	delete [] distance;

	return 0;
}


double autoCorrelation( double inArray[], int size ){

	int lag = 1;
	double acl = 0;
	double autoCorr=1;

	while(autoCorr==std::abs(autoCorr)) {

		double * lagArray = new double[size-lag];
		double * leadArray = new double[size-lag];
		std::copy(inArray+lag,inArray+size,leadArray);
		std::copy(inArray,inArray+(size-lag),lagArray);
		autoCorr=gsl_stats_correlation(lagArray, 1, leadArray, 1, (size-lag));
		if(!(autoCorr==autoCorr)){
			for(int i=0;i<size;i++){
				std::cout<<inArray[i]<<std::endl;
			}
		}
		if(lag%100==0) std::cout<<"Lag = "<<lag<<std::endl;

		lag+=1;
		acl+=2*autoCorr;

		delete [] lagArray;
		delete [] leadArray;
	}

	return acl;
}

void saveToFile(double parameterA[], double parameterB[], double parameterC[], int lag, int size, std::string fileName) {
	
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

double checkForNaN(double ACL) {
	if(!(ACL==ACL)) return 1;
	else return ACL;
}


