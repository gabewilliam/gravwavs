#include <vector>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <string>

#include "Interpolator.h"


//CONSTRUCTORS
//============

Interpolator::Interpolator() {}

Interpolator::Interpolator(std::string gridName, std::string probName) {

	//These read the two input files. The first should be a two column csv,
	//the second should be a one column txt
	this->setGrid(gridName);
	this->setProb(probName);

}


//DESTRUCTOR
//==========

Interpolator::~Interpolator() {}


//SET METHODS
//===========

void Interpolator::setGrid(std::string gridName) {

	//Opens the input files and checks it is there. It will return a file
	//not found message and segfault if it is not found.
	FILE * gridFile;
	gridFile = fopen(gridName.c_str(),"r");
	if(gridFile == NULL) {
		std::cout << "File not found: " << gridName <<std::endl;
		return;
	}

	//Declares some temporary variables used whilst reading in
	double x, y;
	std::vector <double> xTemp;
	std::vector <double> yTemp;
	
	//Reads the file line by line, extracting the x and y values
	while(feof(gridFile) == 0) {
		fscanf(gridFile,"%lf,%lf\n",&x,&y);
		xTemp.push_back(x);
		yTemp.push_back(y);
	}
	
	//Assigns the values to the object's data members
	fX = xTemp;
	fY = yTemp;

	fclose(gridFile);

}

void Interpolator::setProb(std::string probName) {
	
	//Opens the file and checks it is there
	FILE * probFile;
	probFile = fopen(probName.c_str(),"r");
	if(probFile == NULL) {
		std::cout << "File not found: " << probName <<std::endl;
		return;
	}

	//Declares some temporary storage variables
	double z;
	std::vector <double> zCol; //Column vector to store z values in
	
	//Reads the file into a single column vector
	while(feof(probFile) == 0) {
		fscanf(probFile,"%lf,%lf\n",&z);
		zCol.push_back(z);
	}

	//Declares a variable for the number of data points extracted. Usually
	//this will be 900.
	int length = zCol.size();
	
	fclose(probFile);

	//SQRTs the number of data points, to find the number of x and y
	//rows or columns. This value should equal the number of rows in the
	//other input file.
	length = sqrt(length);
	
	//Declares a temporary storage "matrix"
	std::vector <std::vector <double> > zTemp(length); //z array

	int k = 0;
	
	//Fills the matrix in the appropriate sequence
	for(int i = 0; i < length; i++) {

		for(int j = 0; j < length; j++) {
			zTemp[i].push_back(zCol[k]);
			k++;
					
		}

	}

	//Passes the extracted values to the data member
	fZ = zTemp;

}


//FUNCTIONS
//=========

int Interpolator::rAbove(int var, double r) {

	//Finds the number of points in a column (and assumes this is equals the
	//number in a row)
	int nPts = fX.size();
	
	//Finds the grid point higher than the input variable in the x direction
	if(var==1){
		
		//Cycles through x coordinates unti a result is found
		for(int i = 0; i < nPts; i++) {
			
			//Returns -1 if the input point lies outside the grid
			if(((r < fX[i]) && (i==0)) || ((r > fX[i]) && (i==(nPts-1)))) {
				return -1;
			}
			//Returns the column number of the x coordinate above the input
			//x coordinate 
			else if(r < fX[i]) {
				return i;
			}
	
		}
	}

	//Finds the grid point above the input variable in the y direction
	else if(var==2) {
				
		//Cycles through y coordinates until a result is found
		for(int i = 0; i < nPts; i++) {

			//Returns -1 if the input point lies outside the grid.
			if(((r < fY[i]) && (i==0)) || ((r > fY[i]) && (i==(nPts-1)))) {
				return -1;
			}
			//Returns the row number of the y coordinate above the input
			//y value
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
	
	//Finds the column and row indices of the upper bounds of the grid
	//square which contains the input point:
	/*If the point is represented by the dot, then these will be:
			
			  	i-1		i
			j	_________ 
				|		|
				|	.	|
			j-1	|_______| 
		
	*/
	int i,j;
	i = this->rAbove(1,x);
	j = this->rAbove(2,y);

	//Returns -1 if the input point doesn't fit into any grid square
	if((i<0)||(j<0)) { return -1; }
	//Sets the upper and lower x and y bounds of the grid square, along
	//with the column and row indices of the upper corner. Returns 1 to
	//represent success.
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

	//Fins the z values at the corner of the grid square defined by the
	//inputs
	//z11 = z(xBelow,yBelow)
	//z12 = z(xBelow,yAbove)
	//z21 = z(xAbove,yBelow)
	//z22 = z(xAbove,yAbove)
	double z11, z12, z21, z22;
	z11 = fZ[i-1][j-1];
	z12 = fZ[i-1][j];
	z21 = fZ[i][j-1];
	z22 = fZ[i][j];

	//Defines the width of the square in the x and y directions
	double xGap = xAbove - xBelow;
	double yGap = yAbove - yBelow;

	//Defines a few variables to use in claculational steps
	double zx1, zx2, zxy;

	//Interpolates in the x direction, along both edges
	zx1 = ((xAbove - x)/xGap)*z11 + ((x - xBelow)/xGap)*z21;
	zx2 = ((xAbove - x)/xGap)*z12 + ((x - xBelow)/xGap)*z22;

	//Interpolates in the y direction, between the two interpolated x points
	zxy = ((yAbove - y)/yGap)*zx1 + ((y - yBelow)/yGap)*zx2;

	//Returns the interpolated z value
	return zxy;	

}

double Interpolator::estimateProb(double x, double y) {

	//Defines some variables
	double xBelow, xAbove, yBelow, yAbove, pInterp;
	int quadSuccess, i, j;

	//Finds the grid square in which the input point lies
	quadSuccess = this->quadrangilate(x,y,&xBelow,&xAbove,
										&yBelow,&yAbove,&i,&j);
	
	//If the point lays outside the grid, sets the result to zero
	if(quadSuccess < 0) { pInterp = 0; }
	//Otherwise, finds the interpolated value of the distributon at the
	//point
	else {
		pInterp = this->interpolate(x,y,xBelow,xAbove,yBelow,yAbove,i,j);
	}

	//Returns the result
	return pInterp;

}






