#include "gwDataTypes.h"

/*NOTE FOR PARAMETER EXTRACTION PEOPLE:
  This basically uses signal extraction's code to do a FFT on a signal (of
  data type signal). To use it, you'll have to create an Extractor. You can
  then pass the signal you want to transform to it, using the setSignal
  function. You will then need to use the fft function to transform your
  signal. The input to the fft function is the signal to which the output will
  be passed (which could be your original signal, if you want to overwrite it,
  or could be a new object of type signal, if you want to preserve the
  original).
*/

class Extractor
{
	public:
		Extractor();

		//Setup memeber functions
		void setSignal(Signal* input);
		
		//Overloaded function to perform fft on a signal
		void fft(Signal* signalFFT);

	private:
		Signal* mSignalT;
		Signal* mSignalF;

};
