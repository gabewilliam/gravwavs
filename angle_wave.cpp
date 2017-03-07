#include <cmath>

#include "gwReadWrite.h"
#include "gwDataTypes.h"

using namespace std;
 
int main(){
	
	
	Signal wave;
	
	vector<double> freq;
	vector<double> re;
	vector<double> im;
	int i=0;
	
	for(i=0; i < wave.waveform[0].size(); i+=2){
		double f = wave.waveform[0][i];
		double r = wave.waveform[1][i];
		double m = wave.waveform[1][i+1];
		
		freq.push_back(f);
		re.push_back(r);
		im.push_back(m);
	}
	
	
	double fplus,fcross;
	double theta,phi,psi;
	
	cout << "Check page 32 for angle diagram ' https://arxiv.org/pdf/0903.0338.pdf '. \r\n";
	cout << "Enter theta, angle between z axis and wave incoming: \r\n";
	cin >> theta;
	cout << "Enter phi, angle between x axis and wave incoming: \r\n";
	cin >> phi;
	cout << "Enter psi, angle on wave polarization thingy: \r\n";
	cin >> psi;
	
	fplus = (0.5*(1+cos(theta)*cos(theta))*cos(2*phi)*cos(2*psi)) - (cos(theta)*sin(2*phi)*sin(2*psi));
	fcross = (0.5*(1+cos(theta)*cos(theta))*cos(2*phi)*cos(2*psi)) + (cos(theta)*sin(2*phi)*cos(2*psi));
	
	
	vector<double> amp;
	double h;

	for(i=0; i < re.size(); i++){
		//h+ is imag, hx is real
		h = fplus*im[i] + fcross*re[i];
		
		//h+ is real, hx is img
		//h = fplus*re[i] + fcross*im[i];
		
		amp.push_back(h);
	}
	
	
	Signal waveout;
	
	for(i=0; i < freq.size(); i++){
		waveout.waveform[0].push_back(freq[i]);
		waveout.waveform[1].push_back(amp[i]);
	
	}
}
