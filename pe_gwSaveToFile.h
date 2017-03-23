#ifndef PE_GWSAVETOFILE_H
#define PE_GWSAVETOFILE_H

#include <fstream>
#include <vector>

//Useful function used to save the parameter arrays to a file as columns
inline void saveArraysToFile(double parameterA[], double parameterB[], double parameterC[], int thinningFactor, int size, std::string fileName) {
	
	//Opens the output text file
	FILE * outFile;
	outFile = fopen(fileName.c_str(),"w");

	for(int i = 0; i < size; i++){
		if (i%thinningFactor==0){
			fprintf(outFile,"%.15g,%.15g,%.15g\n",parameterA[i],parameterB[i],parameterC[i]);
		}
	}
	
	//Closes the output file
	fclose(outFile);
}

//Useful function used to save parameter vectors (rather than arrays) to a file as columns
inline void saveVectorsToFile(vec_d parameterA, vec_d parameterB, vec_d parameterC, int thinningFactor, int size, std::string fileName) {
	
	//Opens the output text file
	FILE * outFile;
	outFile = fopen(fileName.c_str(),"w");

	for(int i = 0; i < size; i++){
		if (i%thinningFactor==0 && i!=0){
			fprintf(outFile,"%.15g,%.15g,%.15g\n",parameterA[i],parameterB[i],parameterC[i]);
		}
	}
	
	//Closes the output file
	fclose(outFile);
}

#endif //PE_GWSAVETOFILE_H
