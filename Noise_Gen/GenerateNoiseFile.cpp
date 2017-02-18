#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "gwReadWrite.h"
#include "NoiseGeneration.h"

using namespace std;

int main(){
	
	srand(time(NULL));

	double sample;
	
	ofstream outFile;
	
	NoiseGenerator ngen;
	
	outFile.open("noise.csv");
	double s;
	for(int i=0; i<100; i++){
		
		outFile << 0 << ",";
		
	}
	for(int j=100; j<300000; j++){
		
		s=ngen.getSample(0.01*j);
		
		outFile << s << ",";
		
	}
	outFile << 0;
	outFile.close();
	
}
