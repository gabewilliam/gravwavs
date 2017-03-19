#include "Binary.h"
#include <vector>

class Bin {

	public:

		Bin();
		
		Bin(double redshift);

		virtual ~Bin();

		double getRedshift();

		double gettLookback();
		
		double getMergeRateSum();

		double getCounter();

		double CDF(double Z, double z);

		double dMSFRdtdV(double z);	
		
		double Ez(double z);

		double tIntegrand(double z, void * params);

		double tLookback(double z); 
		
		void addToMergeRateSum(double birthrate);

		void addToCounter();

	private:

		std::vector<Binary>fbinaries;

		double fredshift;

		double fmergeRateSum;

		double ftLookback;

		double fcounter;

};
