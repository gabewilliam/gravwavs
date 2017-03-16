#ifndef GWDATATYPES_H
#define GWDATATYPES_H

#include <vector>

//Typedef to keep code readable
typedef std::vector<double> vec_d;

//Template Struct to handle template waveforms & parameters
struct Template 
{
	vec_d waveform[2]; 
	double param[2];
};

//Signal struct for just waveforms without parameters
struct Signal
{
	vec_d waveform[2]; 
};

#endif //GWDATATYPES_H