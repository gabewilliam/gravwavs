#include <vector>
#include <cstdio>
#include <iostream>
#include <cmath>

#include "Interpolator.h"


//CONSTRUCTORS
//============

Interpolator::Interpolator() {}

Interpolator::Interpolator(char * gridName, char * probName) {

	this->setGrid(gridName);
	this->setProb(probName);

}


//DESTRUCTOR
//==========

Interpolator::~Interpolator() {}


//SET METHODS
//===========

void Interpolator::setGrid(char * gridName) {

	FILE * gridFile;
	gridFile = fopen(gridName,"r");
	if(gridFile == NULL) {
		std::cout << "File not found: " << gridName <<std::endl;
		return;
	}

	double x, y;
	std::vector <double> xTemp;
	std::vector <double> yTemp;
	
	while(feof(gridFile) == 0) {
		fscanf(gridFile,"%lf,%lf\n",&x,&y);
		xTemp.push_back(x);
		yTemp.push_back(y);
	}
	
	fX = xTemp;
	fY = yTemp;

	fclose(gridFile);

}

void Interpolator::setProb(char * probName) {
	
	FILE * probFile;
	probFile = fopen(probName,"r");
	if(probFile == NULL) {
		std::cout << "File not found: " << probName <<std::endl;
		return;
	}

	double z;
	std::vector <double> zCol; //Column vector to store z values in
	
	while(feof(probFile) == 0) {
		fscanf(probFile,"%lf,%lf\n",&z);
		zCol.push_back(z);
	}

	int length = zCol.size();
	
	fclose(probFile);
	length = sqrt(length);
	
	std::vector <std::vector <double> > zTemp(length); //z array

	int k = 0;
	
	for(int i = 0; i < length; i++) {

		for(int j = 0; j < length; j++) {
			zTemp[i].push_back(zCol[k]);
			k++;
					
		}

	}

	fZ = zTemp;

}


//FUNCTIONS
//=========

int Interpolator::rAbove(int var, double r) {

	int nPts = fX.size();
		
	if(var==1){
		
		for(int i = 0; i < nPts; i++) {
			
			if(((r < fX[i]) && (i==0)) || ((r > fX[i]) && (i==(nPts-1)))) {
				return -1;
			}
			else if(r < fX[i]) {
				return i;
			}
	
		}
	}
	else if(var==2) {
				
		for(int i = 0; i < nPts; i++) {

			if(((r < fY[i]) && (i==0)) || ((r > fY[i]) && (i==(nPts-1)))) {
				return -1;
			}
			else if(r < fY[i]) {
				return i;
			}
	
		}
	
	}

}

int Interpolator::quadrangilate(double x, double y, double * xBelow, 
								double * xAbove, double * yBelow,
								double * yAbove,
								int * xIndex, int * yIndex) {

	int i,j;
	i = this->rAbove(1,x);
	j = this->rAbove(2,y);
	//std::cout << i << " " << j << std::endl;
	if((i<0)||(j<0)) { return -1; }
	else {
		
		* xBelow = fX[i - 1];
		* xAbove = fX[i];
		* yBelow = fY[j - 1];
		* yAbove = fY[j];
		* xIndex = i;
		* yIndex = j;
		return 1;

	}

}

double Interpolator::interpolate(double x, double y, 
								double xBelow, double xAbove,
								double yBelow, double yAbove,
								int i, int j) {

	//z11 = z(xBelow,yBelow)
	//z12 = z(xBelow,yAbove)
	//z21 = z(xAbove,yBelow)
	//z22 = z(xAbove,yAbove)
	double z11, z12, z21, z22;
	z11 = fZ[i-1][j-1];
	z12 = fZ[i-1][j];
	z21 = fZ[i][j-1];
	z22 = fZ[i][j];

	double xGap = xAbove - xBelow;
	double yGap = yAbove - yBelow;

	double zx1, zx2, zxy;

	zx1 = ((xAbove - x)/xGap)*z11 + ((x - xBelow)/xGap)*z21;
	zx2 = ((xAbove - x)/xGap)*z12 + ((x - xBelow)/xGap)*z22;

	zxy = ((yAbove - y)/yGap)*zx1 + ((y - yBelow)/yGap)*zx2;

	return zxy;	

}

double Interpolator::estimateProb(double x, double y) {

	double xBelow, xAbove, yBelow, yAbove, pInterp;
	int quadSuccess, i, j;
	quadSuccess = this->quadrangilate(x,y,&xBelow,&xAbove,&yBelow,&yAbove,&i,&j);
	
	if(quadSuccess < 0) { pInterp = 0; }
	else {
		pInterp = this->interpolate(x,y,xBelow,xAbove,yBelow,yAbove,i,j);
	}

	return pInterp;

}






