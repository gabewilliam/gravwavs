//                             GRAVITATIONAL WAVE ASTRONOMY                            //
//                                    Data Generation                                  //

#include "gwSigGen.h"

int main(){
	
	// Set masses of components
	double M1 = 30.0;
	double M2 = 26.0;
	// Set distance of system
	double Dist = 500.0;
	// Set time of arrival of signal
	double initTime = 0.5;
	// Set phase of wave at arrival
	double initPhase = 0.0;
	// Set minimum detector frequency (MHz)
	double fMin = 20;
	
	// Set up struct to store characteristic information of signal
	parameters PARAMS;
	parameters *P = &PARAMS;
	
	// Set the characteristic parameters
	setFundamentalParameters(initTime, initPhase, fMin, M1, M2, Dist, P);
	
	
    // Set up struct to store vector of final data
	vector<Signal> signalVector;
	vector<Signal> *S = &signalVector;
	
	// Compute the gravitational wave
	int waveDone = gravitationalWave(P, S);
	
	// Set a vector storing the 'template' version of the data
	Template temp;	
	
	// Assign template vector to signal vector
	temp.param[0] = M1;
	temp.param[1] = M2;

	for(unsigned int i = 0; i < signalVector[0].waveform[0].size(); i++){
		temp.waveform[0].push_back(signalVector[0].waveform[0][i]);
		temp.waveform[1].push_back(signalVector[0].waveform[1][i]);
	}
	// Push into template storage vector
	vector<Template> templateVector;
	templateVector.push_back(temp);
	
	// Assign pointer to template vector
	vector<Template> *T = &templateVector;
	
	int OUTPUT;
	
	if (waveDone == SUCCESS){
		
		// Set up the names of the file
		ostringstream s1, s2, s3;
		s1 << M1;
		s2 << M2;
		string del = "csv";
		s3 << del;
		std::string sigFileIdentity = "sig_" + s1.str() + "_" + s2.str() + "." + s3.str();
		std::string tempFileIdentity = "temp_" + s1.str() + "_" + s2.str() + "." + s3.str();
		
		if(saveSignals(sigFileIdentity, S, csv)){
			cout << "Signal stored successfully" << endl;
			OUTPUT = SUCCESS;
		}
		else{
			cout << "There was an error" << endl;
			OUTPUT = FAILURE;
		} 
		
		if(saveTemplates(tempFileIdentity, T, csv)){
			cout << "Template stored successfully" << endl;
			OUTPUT = SUCCESS;
		}
		else{
			cout << "There was an error" << endl;
			OUTPUT = FAILURE;
		} 
		
	}else{
		OUTPUT = FAILURE;
	}
	
	
	return OUTPUT;
}