//                             GRAVITATIONAL WAVE ASTRONOMY                            //
//                                    Data Generation                                  //

//   This program allows the user to generate the Fourier-domain gravitational wave    //
//     signature of a coalescing binary system with specific masses at a specific      //
//                                     distance.                                       //

#include "gwSigGen.h"

int main(){
	
	// Obtain parameters from user
	double M_min, M_max, Dist, nSteps;
	// Mass 1
	cout << "Please enter minumum mass of one body (in Solar Masses): ";
	cin >> M_min;
	// Mass 2
	cout << endl << "Please enter the maximum mass of one body (in Solar Masses) : ";
	cin >> M_max;
	// Distance
	cout << endl << "Please enter the distance to the binary systems (in MPc) : ";
	cin >> Dist;
	cout << endl << "Please enter the number of steps between minimum and maximum mass: ";
	cin >> nSteps;
	cout << endl;
	
	// Safety checking
	if (M_min == 0.0 || M_max == 0.0 || Dist == 0.0 || M_min >= M_max){
		cout << "Cannot compute" << endl;
		return FAILURE;
	}
	
	// Set time of arrival of signal
	double initTime = 10.0;
	// Set phase of wave at arrival
	double initPhase = 0.0;
	// Set minimum detector frequency (MHz)
	double fMin = 20;	
	// Angle of inclination
	double incline = M_PI/3.0;
	
	double dM = (M_max - M_min)/nSteps;
	
	// Storage vector for all signals
	vector<Signal> signalStorage;
	vector<Signal> * S = &signalStorage;
	
	// Storage vector for all templates
	vector<Template> templateStorage;
	vector<Template> * T = &templateStorage;	
	
	for (double M1 = M_min; M1 <= M_max; M1+=dM){
		
		for (double M2 = M1; M2 <= M_max; M2+=dM){
		
			// Set up struct to store characteristic information of signal
			parameters PARAMS;
			parameters *P = &PARAMS;
			
			// Set the characteristic parameters
			setFundamentalParameters(initTime, initPhase, fMin, M1, M2, Dist, incline, P);
			
			// Compute the gravitational wave
			if (gravitationalWave(P, S) == SUCCESS){
				
				Template oneTemp;
				oneTemp.param[0] = M1;
				oneTemp.param[1] = M2;
				
				Signal oneSig = signalStorage.back();
				
				for (unsigned int indx = 0; indx < oneSig.waveform[0].size(); indx++){
					
					oneTemp.waveform[0].push_back(oneSig.waveform[0][indx]);
					oneTemp.waveform[1].push_back(oneSig.waveform[1][indx]);
					
				}					
				templateStorage.push_back(oneTemp);
				
				cout << "Mass 1: " << M1 << ", Mass 2: " << M2 << endl;
			}
			else{
				cout << "Error with masses " << M1 << " and " << M2 << endl;			
			}
		}
	}
		
	// Set up the name of the signal file
	ostringstream s1, s2, s3;
	s1 << M_max;
	s2 << M_min;
	string del = "csv";
	s3 << del;
	std::string tempFileIdentity = "temp_" + s1.str() + "_" + s2.str() + "." + s3.str();
	
	if(saveTemplates(tempFileIdentity, T, csv)){
		cout << "Templates stored successfully in the file " << tempFileIdentity << "." << endl;
	return SUCCESS;
	}
	else{
		cout << "There was an error with storage." << endl;
		return FAILURE;
	} 		
}