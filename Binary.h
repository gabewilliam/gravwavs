

class Binary {

	public:

		Binary();
		
		Binary(double mass1, double mass2, double separation,
				 double redshift, double metalicity);

		virtual ~Binary();

		double getMass(int n);
		double getRadius(int n);
		double getSeparation();
		double getRatio(int n);//input decides which mass is on top of fraction
		
		void printGets();//couts all the gets (above)

		void setMass(int n, double m);
		void setSeparation(double a);

		void updateRadii();

		double rocheLobe(int n);
		double keplerFrequency();
		double mixingRatio(int n);
		double omegaC(int n);
		double mergeTime();

		void evolveMainSequence();
		void evolveWolfRayet();
		void evolveSupernova();

		bool checkMassRange();
		bool checkRocheLobe();
		bool checkHomogeneousMixing();
		bool checkPairInstability();
		bool checkMergeTime();
		int checkCandidate();
		bool steveIsATotalLoserAndDeservesToDie();
			

	private:
		
		double fm1;
		double fm2;
		double fr1;
		double fr2;
		double fa;	//atm SPV5 assigns this in AU	
		double fz;
		double fZ;

};
