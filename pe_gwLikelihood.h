#ifndef PE_GWLIKELIHOOD_H
#define PE_GWLIKELIHOOD_H

#include "gwDataTypes.h"
#include "gwNoiseGen.h"
#include <string>

long double likelihood( double, double, double, std::string, AligoZeroDetHighP );//This calculates the likelihood function for 2 parameters m1 and m2
vec_d ParameterFunction( double, double, double, vec_d );//This function is used to create model function ht to compare with the data
vec_d NoiseFunction( vec_d, AligoZeroDetHighP ); //This is used to create the noise probability as a function of frequency

#endif //PE_GWLIKELIHOOD_H
