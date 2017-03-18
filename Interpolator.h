
class Interpolator {

	public:

		Interpolator();
		Interpolator(char * gridName, char * probName);

		virtual ~Interpolator();

		void setGrid(char * gridName);
		void setProb(char * probName);

		int rAbove(int var, double r);
		int quadrangilate(double x, double y, double * xBelow, double * xAbove, 
							double * yBelow, double * yAbove,
							int * xIndex, int * yIndex);
		double interpolate(double x, double y, double xBelow, double xAbove,
							double yBelow, double yAbove, int i, int j);
		double estimateProb(double x, double y);

	private:

		std::vector <double> fX; //Stores the X values which form the input grid
		std::vector <double> fY; //Stores the Y values which form the input grid
		std::vector <std::vector <double> > fZ; //Stores the probabilities for the grid points

		/*
		The storage format is:
		
		The values X and Y are stored in the two column vectors.

		The probabilities of (X,Y) pairs are stored in the array Z. The
		format they are stored in is:
		   |	
		  i|       X
		_j_|_______________
		   |
		   |
		Y  |
		   |
		   |
		
		So Z[i][j] would return p(X[i],Y[j]).
		*/


};
