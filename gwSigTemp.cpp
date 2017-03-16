//                             GRAVITATIONAL WAVE ASTRONOMY                            //
//                                    Data Generation                                  //

#include "gwSigGen.h"

int main(){
	
	
	// Generate all data from user inputs:
	double M1, M2, Dist, totalTime, initTime, initPhase, fMin, incline;
	// Set masses of components
	cout << "Please enter the first mass (in Solar Masses): ";
	cin >> M1;
	cout << endl << "Please enter the second mass (in Solar Masses): ";
	cin >> M2;
	// Set distance of system
	cout << endl << "Please enter the distance to the system (in Megaparsecs): ";
	cin >> Dist;
	// Set time of arrival of signal
	cout << endl << "Please enter the total time of the sample (in Seconds): ";
	cin >> totalTime;
	// Ensure signal is in the centre of the data set
	initTime = totalTime/2.0;	
	// Set phase of wave at arrival
	cout << endl << "Please enter the phase of the wave upon arrival (in Radians): ";
	cin >> initPhase;
	// Set minimum detector frequency (MHz)
	cout << endl << "Please enter the minimum frequency of the detector sensitivity (in Hz): ";
	cin >> fMin;
	// Angle of inclination
	cout << endl << "Please enter the angle of inclination of the binary system (in Radians): ";
	cin >> incline;
	cout << endl;
	
	// ADD SAFETY CHECKS
	
	// Set up struct to store characteristic information of signal
	parameters PARAMS;
	parameters *P = &PARAMS;
	
	// Set the characteristic parameters
	setFundamentalParameters(initTime, 
							 totalTime,
							 initPhase, 
							 fMin, 
							 M1, 
							 M2, 
							 Dist, 
							 incline, 
							 P);
		
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
			cout << "Signal stored successfully in the file " << sigFileIdentity << "." << endl;
			OUTPUT = SUCCESS;
		}
		else{
			cout << "There was an error storing the signal." << endl;
			OUTPUT = FAILURE;
		} 
		
		if(saveTemplates(tempFileIdentity, T, csv)){
			cout << "Template stored successfully in the file " << tempFileIdentity << "." << endl;
			OUTPUT = SUCCESS;
		}
		else{
			cout << "There was an error storing the template." << endl;
			OUTPUT = FAILURE;
		} 
		
	}else{
		OUTPUT = FAILURE;
	}
	
	
	return OUTPUT;
}