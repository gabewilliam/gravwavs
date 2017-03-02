

class Binary {

	public:

		Binary();
		
		Binary(double mass1, double mass2, double separation);

		virtual ~Binary();

		double getMass(int n);
		double getRadius(int n);
		double getSeparation();
		double getRatio(int n);//input decides which mass is on top of fraction
		
		void setMass(int n, double m);
		void setSeparation(double a);

		void updateRadii();

		double rocheLobe(int n);
		double keplerFrequency();
		double mixingFrequency(int n);
		double mergeTime();

		void evolveMainSequence();
		void evolveWolfRayet();
		void evolveSupernova();

		bool checkRocheLobe(int n);
		bool checkHomogeneousMixing(int n);
		bool checkPairInstability();
		bool checkMergeTime();
			

	private:
		
		double fm1;
		double fm2;
		double fr1;
		double fr2;
		double fa;		

}
