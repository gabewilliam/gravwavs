#include <stdio.h>
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

using namespace std;
 
int main(){
	
	//Setting up input stream for Noise file
	std::ifstream in_fileN;
	in_fileN.open("InverseNoise.txt");	
	
	//Checking if stream is setup succesfully
	if(in_fileN.fail())
	{
		cout << "Noise file not found" << endl;
		return false;
	}
	
	
	
	
	
	vector<double> NTime;
	vector<double> NAmp;
	
	string lineN;
	
	string a,b;
	
	double Nt,Na;	
	
	while(getline(in_fileN, lineN)){
		
		stringstream linestream(lineN);
		
		//loading parametres from first line of template
		getline(linestream,a, ',');
		getline(linestream,b, ',');
		
		//convert string to double
		stringstream converta(a);
		converta >> Nt;
		stringstream convertb(b);
		convertb >> Na;
		
		//add to vector
		NTime.push_back(Nt);
		NAmp.push_back(Na);
		
	}
	
	//______________________________
	
	vector<Signal> Wave;
	vector<Signal> *wave = &Wave;
	
	loadSignals("signal.csv", wave, csv);
	
	vector<double> WTime;
	vector<double> WAmp;
	
	Signal wsig;
	
	for(int j=0; j < Wave[0].waveform[0].size(); j++){
		double gh = Wave[0].waveform[0][j];
		double dh = Wave[0].waveform[1][j];
		
		WTime.push_back(gh);
		WAmp.push_back(dh);
	}
	
	//______________________________
	
	vector<double> Time;
	vector<double> Amp;
	
	//Goes through vectors until both empty
	while( !WTime.empty() || !NTime.empty() ){
		
		if( !(fabs(WTime[0]-NTime[0]) > 0 )){
			Time.push_back(WTime[0]);
			Amp.push_back(WAmp[0]+NAmp[0]);
		
			WTime.erase(WTime.begin());
			WAmp.erase(WAmp.begin());
		
			NTime.erase(NTime.begin());
			NAmp.erase(NAmp.begin());
		}
	
		//if empty vectors adds other vector
		else if(NTime.empty()){
			Time.push_back(WTime[0]);
			Amp.push_back(WAmp[0]);
		
			WTime.erase(WTime.begin());
			WAmp.erase(WAmp.begin());
		}
	
		else if(WTime.empty()){
			Time.push_back(NTime[0]);
			Amp.push_back(NAmp[0]);
		
			NTime.erase(NTime.begin());
			NAmp.erase(NAmp.begin());
		}
	
		//Adds lowest time value to end of vector for increasing time
		else if(WTime[0]<NTime[0]){
			Time.push_back(WTime[0]);
			Amp.push_back(WAmp[0]);
		
			WTime.erase(WTime.begin());
			WAmp.erase(WAmp.begin());
		}
	
		else if(NTime[0]<WTime[0]){
			Time.push_back(NTime[0]);
			Amp.push_back(NAmp[0]);
		
			NTime.erase(NTime.begin());
			NAmp.erase(NAmp.begin());
		}
	}
	
	//______________________________
	//vector<Signal> Noisy;
	//vector<Signal>* noisy;
	
	vector<Signal> Noise;
	vector<Signal> *noisy = &Noise;
	Signal sig;
	
	for(int k=0; k < Time.size(); k++){
		sig.waveform[0].push_back(Time[k]);
		sig.waveform[1].push_back(Amp[k]);
		noisy->push_back(sig);
		
	}
	
	//saveSignals("NoisyWave.csv", noisy, csv);
	//above file seemingly too big 
	
	//Write to output file
	ofstream oFile;
	oFile.open("noisewave2.txt", std::ios_base::app);
	
	for (int i=0; i < Time.size(); i++){
		oFile << Time[i] << "," << Amp[i] << "\n" ;
	}
	
	cout << "Output finished and written to 'NoisyWave.txt' \n" << endl;
	
	oFile.close();
	
}
