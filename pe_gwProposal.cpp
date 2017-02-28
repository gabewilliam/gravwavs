//as per "Veitch et al (2015)"
#include "pe_gwProposal.h"

#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

double Proposal ( double current, int index, int gaRatio, int deRatio, int ejRatio, int agRatio, int gwRatio, double gaSigma, double dePenultimate, double deAntepenultimate, int deN ){
	
	std::vector( int ) proposalOrder = ProposalOrder( deRatio, ejRatio, agRatio, gwRatio );
	int size = proposalOrder.size();
	
	int i = index;
	while( i > size ) i = i - size;
	int proposalNum = proposalOrder[ i - 1 ];
	
	double proposed;
	switch( proposalNum ){
		case 1 : 
			proposed = Gaussian( current, gaSigma );
			break;

		case 2 : 
			proposed = EigenvectorJump( current );
			break;

		case 3 : 
			proposed = AdaptiveGaussian( current );
			break;

		case 4 : 
			proposed = GravWavesSpecific( current );
			break;
		
		case 5 : 
			proposed = DifferentialEvolution( current, dePenultimate, deAntepenultimate );
			break;	
	}
	return proposed;
}

double Gaussian( double current, double sigma ){
	double gaussian = gsl_ran_gaussian(normGen, sigma);
	double proposed = current + gaussian;

	return proposed;
}

double EigenvectorJump( double );
double AdaptiveGaussian( double );
double GravWavesSpecific( double );

double DifferentialEvolution( double current, double penultimate, double antepenultimate, int N ){
	srand( time( NULL ) );	
	double r = ( ( double ) rand() / ( RAND_MAX ) );
	double gamma;
	if( r > 0.5 ) gamma = 1;
	else gamma = gsl_ran_gaussian( normGen, (2.38 / sqrt( 2 * N )) );

	double proposed = current + gamma * ( penultimate - antepenultimate )
}

std::vector( int ) ProposalOrder( int gaRatio, int deRatio, int ejRatio, int agRatio, int gwRatio ){
	
	std::vector( int ) findMax;
	findMax.push_back( gaRatio );	
	findMax.push_back( deRatio );
	findMax.push_back( ejRatio );
	findMax.push_back( agRatio );
	findMax.push_back( gwRatio );
	int max = *max_element( findMax.begin(), findMax.end() );
	

	std::vector( int ) proposalOrder;
	for( int i = 0; i < max; i++ ){
		if( gaRatio >= i ) proposalOrder.push_back( 1 );
		if( ejRatio >= i ) proposalOrder.push_back( 2 );
		if( agRatio >= i ) proposalOrder.push_back( 3 );
		if( gwRatio >= i ) proposalOrder.push_back( 4 );
		if( deRatio >= i ) proposalOrder.push_back( 5 );
	}
return proposalOrder;
}
