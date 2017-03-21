#ifndef PE_GWLIKELIHOOD_H
#define PE_GWLIKELIHOOD_H

#include <string>

#include "gwReadWrite.h"

long double likelihood( double, double, double, std::string, Signal, vec_d);//This calculates the likelihood function for 2 parameters m1 and m2
vec_d ParameterFunction( double, double, double, vec_d );//This function is used to create model function ht to compare with the data

#endif //PE_GWLIKELIHOOD_H
