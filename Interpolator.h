
class Interpolator {

	public:

		//Constructors, the latter of which takes the input filenames as
		//arguments
		Interpolator();
		Interpolator(std::string gridName, std::string probName);

		//Destructor
		virtual ~Interpolator();

		//Functions to set the input files
		void setGrid(std::string gridName);
		void setProb(std::string probName);

		//Function which finds the index of the coordinate (in one direction)
		//which lies above the input point
		int rAbove(int var, double r);
		//Function which finds the grid square in which the input point lies
		int quadrangilate(double x, double y, double * xBelow, double * xAbove, 
							double * yBelow, double * yAbove,
							int * xIndex, int * yIndex);
		//Function whcih finds the interpolated value of the distribution at
		//the input point, given the grid square it lies in
		double interpolate(double x, double y, double xBelow, double xAbove,
							double yBelow, double yAbove, int i, int j);
		//Top level function which takes the input point and finds the
		//interpolated value of the function there.
		double estimateProb(double x, double y);

	private:

		//Stores the X values which form the input grid
		std::vector <double> fX;
		//Stores the Y values which form the input grid
		std::vector <double> fY;
		//Stores the probabilities for the grid points
		std::vector <std::vector <double> > fZ; 

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

		The number of X and Y elements should equal the square root of the
		number of elements in Z. There is no error checking for this,
		however. It is taken on trust that the correct input file format will
		be supplied.
		*/


};
