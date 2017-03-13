//File to combine noise and wave

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "gwReadWrite.h"
#include "gwDataTypes.h"
#include "gwSigGen.h"
#include "NoiseGeneration.h"

using namespace std;
 
int main(){
	//iterator
	int i = 0;
	
	//Generate wave
	
	cout << "Enter first mass of component (Solar Masses) 26/30: \r\n";
	double M1, M2;
	cin >> M1;
	cout << "Enter other mass of component (Solar Masses) 30/26: \r\n";
	cin >> M2;
	
	cout << "Enter distance of system (MPc) 500: \r\n";
	double Dist;
	cin >> Dist;
	
	cout << "Enter total time of signal (s) 20/100: \r\n";
	double totalTime;
	cin >> totalTime;
	
	cout << "Enter time of arrival of signal (s) 10: \r\n";
	double initTime;
	cin >> initTime;
	
	cout << "Enter phase of wave at arrival (rads) 0: \r\n";
	double initPhase;
	cin >> initPhase;
	
	cout << "Enter minimum detector frequency (Hz) 10: \r\n";
	double fMin;
	cin >> fMin;
	
	cout << "Enter angle of inclination of binary (rads) 0/10: \r\n";
	double psi;
	cin >> psi;
	
	cout << "Enter angle of incident wave theta 0/3: \r\n";
	double theta, phi;
	cin >> theta;
	cout << "Enter angle of incident wave phi 0/5: \r\n";
	cin >> phi;
	
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
							 theta,
							 psi,
							 phi,
							 P);
		
	// Set up signal struct to store single signal
	Signal singleSignal;
	
	// Set up complex number representing the wave itself
	complex<double> complexWave = 0.0;
	complex<double> * gWave = &complexWave;
	
	// Initialise frequency and amplitude
	double freqw = 0.0;
	double amp_Hz = 0.0;
	
	// Set variable for storing frequency in Hz
	double freq_Hz = 0.0;
	
	cout << "Enter maximum frequency 4095/?: \r\n";
	double maxFreq;
	cin >> maxFreq;
	maxFreq = maxFreq/C_CONST;
	
	// Set the number of frequency bins to loop over
	double nFreq = maxFreq/P->df + 1.0;
	
	for (i = 0; i < nFreq; i++){
		
		freqw = double(i)*P->df;
		
		amp_Hz = updatedAmplitude(P, freqw, gWave);
		
		// Convert frequency back to Hz
		freq_Hz = freqw*C_CONST;		
		
		singleSignal.waveform[0].push_back(freq_Hz);	
		singleSignal.waveform[1].push_back(amp_Hz);
		
	}
	
	vector<double> Wfreq;
	vector<double> Wre;
	vector<double> Wim;
	
	for(i=0; i < singleSignal.waveform[0].size(); i+=2){
		double f = singleSignal.waveform[0][i];
		double r = singleSignal.waveform[1][i];
		double m = singleSignal.waveform[1][i+1];
		
		Wfreq.push_back(f);
		Wre.push_back(r);
		Wim.push_back(m);
	}
	
	ofstream oFileW;
	oFileW.open("AVWavey.csv");
	for (int i=0; i < Wfreq.size(); i++){
		oFileW << Wfreq[i] << "," << Wre[i] << "," << Wim[i] << "\r\n" ;
	}
	cout << "wave finished and written to 'AVWavey.csv' \n" << endl;
	
	//______________________________
	//Generate Noise
	
	srand(time(NULL));

	cout<<"Select type of noise:\r\n\t"
	<< "1) ALIGO Schutz\r\n\t"
	<< "2) ALIGO Zero Det High P\r\n\t"
	<< "3) ALIGO Zero Det Low P\r\n\t"
	<< "4) ALIGO NSNS Opt\r\n\t"
	<< "5) ALIGO No SRM\r\n\t"
	<< "6) ALIGO High Freq\r\n\t"
	<< "7) ALIGO BHBH 20 Deg\r\n";
	
	int choice;
	cin >> choice;
	
	NoiseGenerator *nGen;
	
	if(choice == 1)		{	nGen =  new AligoSchutz();		}
	else if(choice == 2){	nGen = new AligoZeroDetHighP();	}
	else if(choice == 3){	nGen = new AligoZeroDetLowP();	}
	else if(choice == 4){	nGen = new AligoNsnsOpt();		}
	else if(choice == 5){	nGen = new AligoNoSrm();		}
	else if(choice == 6){	nGen = new AligoHighFreq();		}
	else if(choice == 7){	nGen = new AligoBhbh20Deg();	}
	else				{	nGen = new AligoSchutz();		}
	
	vector<double> *freq = new vector<double>;
	vector<Complex> *noise = new vector<Complex>;
	
	double fMax, fInc;
	fMax = maxFreq;
	fInc = P->df;
	
	nGen->genSpectrum(freq, noise, fMax, fInc);
	
	vector<double> Nre;
	vector<double> Nim;
	vector<double> Nfreq;
	cout << "a"<<endl;
	for(i=0; i < freq->size(); i++){
		Nfreq.push_back(freq->at(i));
	}
	cout << "b"<<endl;
	//Change units of noise into same units as wave
	for (i=0; i < noise->size(); i++){
		Nre.push_back(((noise->at(i)).real)/sqrt(fabs(Nfreq[i])));
		Nim.push_back(((noise->at(i)).imag)/sqrt(fabs(Nfreq[i])));
	}
	
	cout << "c"<<endl;
	
	ofstream oFileN;
	oFileN.open("AVNoisy.csv");
	for (int i=0; i < Nfreq.size(); i++){
		oFileN << Nfreq[i] << "," << Nre[i] << "," << Nim[i] << "\r\n" ;
	}
	cout << "noise finished and written to 'AVNoisy.csv' \n" << endl;
	
	//______________________________
	//Addition of values
	
	vector<double> frq;
	vector<double> re;
	vector<double> im;
	cout << "d"<<endl;
	//Goes through vectors until both empty
	while( !Wfreq.empty() || !Nfreq.empty() ){
		
		//if Wfreq and Nfreq are same, add values
		if( !(fabs(Wfreq[0]-Nfreq[0]) > 0 )){
			frq.push_back(Wfreq[0]);
			re.push_back(Wre[0]+Nre[0]);
			im.push_back(Wim[0]+Nim[0]);
			
			Wfreq.erase(Wfreq.begin());
			Wre.erase(Wre.begin());
			Wim.erase(Wim.begin());
		
			Nfreq.erase(Nfreq.begin());
			Nre.erase(Nre.begin());
			Nim.erase(Nim.begin());
		}
	
		//if empty vectors adds other vector
		else if(Nfreq.empty()){
			frq.push_back(Wfreq[0]);
			re.push_back(Wre[0]);
			im.push_back(Wim[0]);
		
			Wfreq.erase(Wfreq.begin());
			Wre.erase(Wre.begin());
			Wim.erase(Wim.begin());
		}
	
		else if(Wfreq.empty()){
			frq.push_back(Nfreq[0]);
			re.push_back(Nre[0]);
			im.push_back(Nim[0]);
		
			Nfreq.erase(Nfreq.begin());
			Nre.erase(Nre.begin());
			Nim.erase(Nim.begin());
		}
	
		//Adds lowest time value to end of vector for increasing time
		else if(Wfreq[0]<Nfreq[0]){
			frq.push_back(Wfreq[0]);
			re.push_back(Wre[0]);
			im.push_back(Wim[0]);
		
			Wfreq.erase(Wfreq.begin());
			Wre.erase(Wre.begin());
			Wim.erase(Wim.begin());
		}
	
		else if(Nfreq[0]<Wfreq[0]){
			frq.push_back(Nfreq[0]);
			re.push_back(Nre[0]);
			im.push_back(Nim[0]);
		
			Nfreq.erase(Nfreq.begin());
			Nre.erase(Nre.begin());
			Nim.erase(Nim.begin());
		}
	}
	cout << "e"<<endl;
	//______________________________
	//vector<Signal> Noisy;
	//vector<Signal>* noisy;
	
//	vector<Signal> Wave;
//	vector<Signal> *wavey = &Wave;
//	Signal sig;
//	int k=0;
//	for(i=0; i < freq.size(); i+=2){
//		sig.waveform[0].push_back(freq[k]);
//		sig.waveform[1].push_back(re[k]);
//		noisy->push_back(sig);
		
//		k = k+1;
//	}
	
	//saveSignals("NoisyWave.csv", noisy, csv);
	//above file seemingly too big 
	
	//Write to output file
	ofstream oFile;
	oFile.open("AVNoisyWave.csv");
	cout << "f"<<endl;
	for (int i=0; i < frq.size(); i++){
		oFile << frq[i] << "," << re[i] << "," << im[i] << "\r\n" ;
	}
	
	cout << "Output finished and written to 'AVNoisyWave.csv' \n" << endl;
	
	oFile.close();
	
}
