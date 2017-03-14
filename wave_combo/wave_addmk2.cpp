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
		
	// Set up struct to store vector of final data
	vector<Signal> signalVector;
	vector<Signal> *S = &signalVector;
	
	// Compute the gravitational wave
	int waveDone = gravitationalWave(P, S);
	
	vector<double> Wfreq;
	vector<double> Wre;
	vector<double> Wim;
	
	for(i=0; i < signalVector[0].waveform[0].size(); i+=2){
		double f = signalVector[0].waveform[0][i];
		double r = signalVector[0].waveform[1][i];
		double m = signalVector[0].waveform[1][i+1];
		
		Wfreq.push_back(f);
		Wre.push_back(r);
		Wim.push_back(m);
	}
	
	//wave intial data with columns down page, freq, real, imag
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
	
	
//	ofstream oFile2;
//	oFile2.open("avnoisey2.csv");
	
	vector<double> *freq = new vector<double>;
	vector<Complex> *noise = new vector<Complex>;
	
	//fMax and fInc defined within gravitationalWave method
	double fMax = 4096;
	double fInc = (P->df)*C_CONST;
	
	nGen->genSpectrum(freq, noise, fMax, fInc);
	
//	for(int k=0; k < freq->size(); k++){
//		oFile2 << freq->at(k) << "," << (noise->at(k)).real << "," << (noise->at(k)).imag << "\r\n";
//	}
	
	vector<double> Nre;
	vector<double> Nim;
	vector<double> Nfreq;
	
	for(i=0; i < freq->size(); i++){
		Nfreq.push_back(freq->at(i));
	}
	
	//Change units of noise into same units as wave
	for (i=0; i < noise->size(); i++){
		Nre.push_back(((noise->at(i)).real)/sqrt(fabs(Nfreq[i])));
		Nim.push_back(((noise->at(i)).imag)/sqrt(fabs(Nfreq[i])));
	}
	
	//noise intial data with columns down page, freq, real, imag
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
	
	//______________________________
	//Write to output file
	
	ofstream oFile;
	oFile.open("AVNoisyWave.csv");
	
	//file with columns down page, freq, real, imag
//	for (int i=0; i < frq.size(); i++){
//		oFile << frq[i] << "," << re[i] << "," << im[i] << "\r\n" ;
//	}
	
	//file with rows across, freq first row, real, imag
//	for (int i=0; i < frq.size(); i++){
//		oFile << frq[i] << ",";
//	}
//	oFile << "\r\n";
	
//	for (int i=0; i < re.size(); i++){
//		oFile << re[i] << ",";
//	}
//	oFile << "\r\n";
	
//	for (int i=0; i < im.size(); i++){
//		oFile << im[i] << ",";
//	}
//	oFile << "\r\n";
	
	
	//file with freq1, freq1, next row real, imag (izzie style pack)
	for (int i=0; i < frq.size(); i++){
		oFile << frq[i] << "," << frq[i] << ",";
	}
	oFile << "\r\n";
	
	for (int i=0; i < frq.size(); i++){
		oFile << re[i] << "," << im[i] << ",";
	}
	oFile << "\r\n";
	
	cout << "Output finished and written to 'AVNoisyWave.csv' \n" << endl;
	
	oFile.close();
	
}
