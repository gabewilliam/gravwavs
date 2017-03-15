#include "gwSigGenBasic.h"

int main(){
	
	Parameters params;
	Parameters * p = &params;
	
	return gwGenerateAndPack(30.0,
							 31.0,
							 400.0,
							 10.0,
							 0.0,
							 10.0,
							 100.0,
							 p,
							 SIG);
}