//                             GRAVITATIONAL WAVE ASTRONOMY                            //
//                                    Data Generation                                  //

//   This program allows the user to generate the Fourier-domain gravitational wave    //
//     signature of a coalescing binary system with specific masses at a specific      //
//                                     distance.                                       //

#include "FDGWNU_2.h"

int main(){
	
	// Obtain parameters from user
	double M_min, M_max, Dist;
	// Mass 1
	cout << "Please enter minumum mass of one body (in Solar Masses): ";
	cin >> M_min;
	// Mass 2
	cout << endl << "Please enter the maximum mass of one body (in Solar Masses) : ";
	cin >> M_max;
	// Distance
	cout << endl << "Please enter the distance to the binary systems (in MPc) : ";
	cin >> Dist;
	cout << endl;
	
	// Safety checking
	if (M_min == 0.0 || M_max == 0.0 || Dist == 0.0 || M_min >= M_max){
		cout << "Cannot compute" << endl;
		return FAILURE;
	}
	
	// Set time of arrival of signal
	double initTime = 0.0;
	// Set phase of wave at arrival
	double initPhase = 0.0;
	// Set minimum detector frequency (MHz)
	double fMin = 2E-5;	
	
	double dM = (M_max - M_min)/10.0;
	
	// Declare number of data points to generate for each template
	int N = 1E5;
	
	// Storage vector for all templates
	vector<Template> templateStorage;
	vector<Template> * T = &templateStorage;	
	
	M_max += dM;
	
	for (double M1 = M_min; M1 <= M_max; M1+=dM){
		
		for (double M2 = M_min; M2 <= (M_max - M1 + M_min); M2+=dM){
		
		// Set up struct to store characteristic information of signal
		parameters PARAMS;
		parameters *P = &PARAMS;
		
		// Set the characteristic parameters
		setFundamentalParameters(initTime, initPhase, fMin, M1, M2, Dist, P);
		
		// Compute the gravitational wave
	    if (gravitationalWaveTemps(P, T, N, M1, M2) == SUCCESS){
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
		cout << "Templates stored successfully" << endl;
	return SUCCESS;
	}
	else{
		cout << "There was an error" << endl;
		return FAILURE;
	} 		
}