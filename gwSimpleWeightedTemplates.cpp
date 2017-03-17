//------------------- GRAVITATIONAL WAVE ASTRONOMY -------------------//
//------------------------ DATA GENERATION ---------------------------//

//----------- FREQUENCY DOMAIN TEMPLATE DATA BANK GENERATOR ----------//
//----- Model purely weighted based on mass pairs likely to form -----//
//-------------- through Chemically Homogeneous Evolution ------------// 

//--- Uses a density threshold and centrally located mass to reduce --//
//-------------------- number of templates computed ------------------//

#include "gwSigGenRedshift.h"

int main(){
	
	//-- Set up output file to store visual representation of template 
	//-- distribution
	ofstream visualRep;
	visualRep.open("visualRep.txt", std::ios_base::app);
	
	//-- Set up parameters constant to all templates
	double lumDistance = 500.0;
	double arrivalTime = 5.0;
	double phaseOffset = 0.0;
	double minFreq = 10.0;
	double totalTime = 20.0;
	
	//-- Define things needed for later calculation
	Parameters Params;
	Parameters *P = &Params;
	
	Signal Amplitude;
	Signal *A = &Amplitude;
	
	Signal Complex;
	Signal *C = &Complex;
	
	NoisySignal Noisy;
	NoisySignal *N = &Noisy;
	
	vector<Signal> signalVector;
	vector<Signal> *S = &signalVector;
	
	std::ifstream MassDistribution("filledbins.txt");
	
	if(MassDistribution.fail()){
		std::cerr << "Input file not found" << std::endl;
		return FAILURE;
	}
	
	vector<double> lowMA, highMA, lowMB, highMB, filledBins;
	
	while (!MassDistribution.eof()){
			
		double mA1, mA2, mB1, mB2, nBins;
		
		MassDistribution >> mA1 >>
							mA2 >>
							mB1 >> 
							mB2 >>
							nBins;			
							
		lowMA.push_back(mA1);
		highMA.push_back(mA2);
		lowMB.push_back(mB1);
		highMB.push_back(mB2);
		filledBins.push_back(nBins);
		
	}
		 
	cout << "Masses successfully extracted" << endl;
	//-- Declare vector of templates
	vector<Template> weightedTemplateVector;
	vector<Template> *W = &weightedTemplateVector;
	
	double totalTemplates = 0.0;
	
	//-- Decide density threshold for investigation of region 
	//-- of parameter space
	
	int densityThreshold = 12;
	double dense = 10.0;
	
	for (unsigned int l = 0; l < filledBins.size(); l++){
		
		if (filledBins[l] > densityThreshold){
			totalTemplates += int(filledBins[l]/dense);
		}
		
	}
	
	cout << "Templates to compute: " << totalTemplates << endl;
	
	double thisTemplate = 0.0;
	
	for (unsigned int i = 0; i < lowMA.size(); i++){
		
		if (filledBins[i] > densityThreshold){
			
			int nTemplates = int(filledBins[i]/dense);
			
			for (int t = 0; t < nTemplates; t++){
				
				thisTemplate += 1.0;
				
				double M1 = lowMA[i] + t*(highMA[i] + lowMA[i])/nTemplates;
				
				double M2 = lowMB[i] + t*(highMB[i] + lowMB[i])/nTemplates;
						
				gwSimulateDetection(M1, M2,
									lumDistance,
									arrivalTime,
									phaseOffset,
									minFreq,
									totalTime,
									P, A, C, N, S);
								
				visualRep << M1 << "\t" << M2 << "\t" << nTemplates << endl;
				
				Template oneTemplate;
				oneTemplate.param[0] = M1;
				oneTemplate.param[1] = M2;
						
				Signal lastSignal = signalVector.back();
				
				for (unsigned int indx = 0;
						indx < lastSignal.waveform[0].size();
														indx++){
				oneTemplate.waveform[0].push_back
										(lastSignal.waveform[0][indx]);
				oneTemplate.waveform[1].push_back
										(lastSignal.waveform[1][indx]);
				}
									
				W->push_back(oneTemplate);

				cout << "Progress: " 
					 << thisTemplate*100.0/totalTemplates 
					 << "% of templates computed..." << endl;
			}					
		}
	}

	visualRep.close();
	
	cout << "Signal templates generated." << endl;
	
	//-- Set up the name of the data file
	ostringstream s1, s2, s3;
	
	s1 << "WEIGHTEDTEMP";
	s2 << lowMA[0];
	s3 << highMA.back();
	
	std::string fileIdentity = s1.str() + "_"
							   + s2.str() + "_" 
							   + s3.str() + ".csv";
	
	cout << "Setting up storage file..." << endl;	

	if (saveTemplates(fileIdentity, W, csv)){
		
		cout << "Weighted templates stored successfully in the file " 
			 << fileIdentity << "." << endl;
			 
		return SUCCESS;
	
	}
	else {
		cout << "Failed to store data." << endl;
		
		return FAILURE;
	}
}
