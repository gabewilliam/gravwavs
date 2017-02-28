#ifndef PE_GWPROPOSAL_H
#define PE_GWPROPOSAL_H

double Proposal( double, int, int, int, int, int, int, double, double, double, int );
double Gaussian( double, double );
double EigenvectorJump( double );
double AdaptiveGaussian( double );
double GravWavesSpecific( double );
double DifferentialEvolution( double, double, double, int );
double ProposalOrder( int, int, int, int, int );

#endif //PE_GWPROPOSAL_H
