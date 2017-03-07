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

	AligoSchutz
		- Analytical fit of ALIGO noise curve,
		from Sathyaprakash & Schutz ~  https://arxiv.org/abs/0903.0338v1
	
	AligoZeroDetHighP
		- Numerical noise asd for zero-detuned recycling mirror and high laser power, 
		from https://dcc.ligo.org/LIGO-T0900288/public
		
	AligoZeroDetLowP
		- Numerical noise asd for zero-detuned recycling mirror and low laser power, 
		from https://dcc.ligo.org/LIGO-T0900288/public
		
	AligoHighFreq
		- Numerical noise asd for high frequency, 
		from https://dcc.ligo.org/LIGO-T0900288/public
		
	AligoBhbh20Deg
		- Numerical noise asd for black hole inspiral, 
		from https://dcc.ligo.org/LIGO-T0900288/public
		
	AligoNoSrm
		- Numerical noise asd with no recycling mirror, 
		from https://dcc.ligo.org/LIGO-T0900288/public
		
	AligoNsnsOpt
		- Numerical noise asd for neutron star inspiral, 
		from https://dcc.ligo.org/LIGO-T0900288/public
		
-------------------------------------------------------------------------------------------------------
