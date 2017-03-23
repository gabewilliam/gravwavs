PARAMETER EXTRACTION

-----------------------------------------------------------------------------------

pe_gwSaveToFile.h:

- Contains functions which save three sets of parameters stored in either arrays or std::vectors as columns to a text file
- The thinning factor determines how many elements to skip before saving the next value

- Reads in and out a text file of three columns with equal numbers of rows

-----------------------------------------------------------------------------------

pe_gwMarkovChain.cpp

- Chooses random starting parameter values using gsl RNGs within the prior limits
- Sets up the AligoZeroDetHighP noise curve
- Implements the chirp mass - mass ratio prior from CHE
- Carries out the MCMC routine for chirp mass, mass ratio and luminosity distance
- Saves the raw parameter data using the pe_gwSaveToFile.h function

- Reads in a data file of two rows: frequency in row 1 and amplitude in row 2
	- We trim the data up to a frequency value of 9Hz as Gabe's noise curves are zero padded up to this value
	- This zero padding causes NaNs in our code when we later compute the likelihood
	- We also only take in the positive frequencies and real parts of the waveform as these also cause errors in the code
- Reads out a text file of three columns with equal numbers of rows
- Takes in the number of samples as input (minimum of 100 but ideally a value of 1 million is used time permitting)

-----------------------------------------------------------------------------------

pe_gwLikelihood.cpp / pe_gwLikelihood.h


- Computes the likelihood function for the posterior calculation in the MCMC
- Also calculates the signal to noise ratio
- Produces the model data signal using Izzie's gwSigGen.h functions
- The inputs and outputs for this program are implicitly carried out within the pe_gwMarkovChain.cpp file

-----------------------------------------------------------------------------------

pe_gwCorrelate.cpp

- Stand-alone program which takes in the raw data from the pe_gwMarkovChain.cpp file and calculates a thinning factor
- This thinning factor is obtained through cross-correlation of the parameter arrays with themselves

- Reads in and out a text file of three columns with equal numbers of rows using the pe_gwSaveToFile.h file

-----------------------------------------------------------------------------------