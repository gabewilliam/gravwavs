#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <gsl/gsl_statistics_double.h>

#include "pe_gwSaveToFile.h"

double autoCorrelation(double [], int);

int main(){
	
	//INSERT PARAMETER FILE NAME HERE
	std::string fileName = "RawData.txt";
	
	//Finds number of rows in file
	std::ifstream lineFile;
	lineFile.open( fileName.c_str() );

	int N = 0;
	std::string line;
	std::vector< std::string > input;
	
	while( std::getline( lineFile, line )){
		input.push_back( line );
		++N;
	}

	//Opens file to be read to arrays
	FILE * inFile;
	inFile=fopen( fileName.c_str() ,"r" );
	
	//Initialises parameter arrays
	double * chirpMass = new double[N];
	double * massRatio = new double[N];
	double * distance = new double[N];

	//Reads raw data file into parameter arrays
	for( int i = 0; i < N; i++ ){
		fscanf( inFile, "%lf %*c %lf %*c %lf\n", &chirpMass[i], &massRatio[i], &distance[i] );
	}

	//Calculates autocorrelation for each parameter
	double corrFactorChirp = autoCorrelation( chirpMass, N );
	double corrFactorRatio = autoCorrelation( massRatio, N );
	double corrFactorDistance = autoCorrelation( distance, N );
	
	//Prints ACL values
	std::cout << "Chirp Mass: " << corrFactorChirp << std::endl;
	std::cout << "Mass Ratio: " << corrFactorRatio << std::endl;
	std::cout << "Distance: " << corrFactorDistance << std::endl;

	//Picks the largest ACL value
	double acls [] = {corrFactorChirp, corrFactorRatio, corrFactorDistance};
	int aclMax = *(std::max_element(acls,acls+3));
	
	//Saves parameters to a file using the maximum ACL value as a thinning factor
	saveArraysToFile( chirpMass, massRatio, distance, aclMax, N, "ThinnedData.txt" );

	//Deallocates memory
	delete [] chirpMass;
	delete [] massRatio;	
	delete [] distance;

	return 0;
}


double autoCorrelation( double inArray[], int size ){

	//Initialises variables
	int lag = 1;
	double acl = 0;
	double autoCorr=1;

	//Keeps looping until the autocorrelation is no longer positive
	while(autoCorr==std::abs(autoCorr)) {

		//Redefines arrays (the arrays decrease in size each iteration so have to be redeclared)
		double * lagArray = new double[size-lag];
		double * leadArray = new double[size-lag];

		//Removes elements from the end of one array and from the start of the other then cross-correlates them
		std::copy(inArray+lag,inArray+size,leadArray);
		std::copy(inArray,inArray+(size-lag),lagArray);
		autoCorr=gsl_stats_correlation(lagArray, 1, leadArray, 1, (size-lag));

		//Increments the lag and calculates the ACL
		lag+=1;
		acl+=2*autoCorr;

		//Deallocates memory
		delete [] lagArray;
		delete [] leadArray;
	}

	return acl;
}
