#include "Binary.h"
#include <vector>

class Bin {

	public:

		Bin();
		
		Bin(double redshift,double (*func)(double,void*));

		virtual ~Bin();

		double getRedshift();

		double gettLookback();
		
		double getMergeRateSum();

		double getCounter();

		double CDF(double Z, double z);

		double dMSFRdtdV(double z);	

		double Integrator(double z,double (*func)(double,void*)); 
		
		void addToMergeRateSum(double birthrate);

		void addToCounter();

	private:

		std::vector<Binary>fbinaries;

		double fRedshift;

		double fMergeRateSum;

		double ftLookback;

		double fCounter;

};
