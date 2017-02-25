-------------------------------------------------------------------------------------------------------
						NOISE GENERATION
-------------------------------------------------------------------------------------------------------

UNIX compile:

	'g++ -o GenerateNoiseFile GenerateNoiseFile.cpp NoiseGeneration.cpp'
	
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
	
Noise generators can be used to generate a spectrum of noise samples with the method .genSpectrum,
Parameters:

	*freqs
		- Pointer to vector<double> of frequencies
	*noise
		- Pointer to vector<Complex> of real and imaginary data
	fMax
		- Maximum frequency to go up to
	fInc
		- Frequency increment

-------------------------------------------------------------------------------------------------------
	
Noise Generator Classes:

	ALIGOSchutz
		- Analytical fit of ALIGO noise curve,
		from Sathyaprakash & Schutz ~  https://arxiv.org/abs/0903.0338v1
	
	ALIGOZeroDetHighP
		- Numerical noise asd for zero-detuned recycling mirror and high laser power, 
		from https://dcc.ligo.org/LIGO-T0900288/public
		
-------------------------------------------------------------------------------------------------------
