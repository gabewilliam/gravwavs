void Extractor::Convolution(std::vector<Template>* output)
{
	int I, J, pn;
	double imagResult;
	double realResult;
	double noise;
	double noisePower;
	double power;
	double SNR;
	
	I = mTemplates->size();
	J = mSignalF->waveform[1].size();

	pn = 1;
	
	for(int i=0; i<I; i++)
	{	
		Template* temp = &mTemplates[0][i];
		
		Template op;

		op.param[0] = temp->param[0];
		op.param[1] = temp->param[1];
		
		SNR = 0;
		
		//not sure if this is looping through all of the elements correctly;
		//the alternative for the more sensible frequency packing is commented out
		for(int j = 0 ; j < J/2; j++){
			
			/*
				
			noise = getASD(mSignalF->waveform[0][j]);

			op.waveform[0].push_back(mSignalF->waveform[0][j]);
			
			noisePower = noise*noise;

			*/	
			noise = getASD(mSignalF->waveform[0][2*j]);

			op.waveform[0].push_back(mSignalF->waveform[0][2*j]);
			op.waveform[0].push_back(mSignalF->waveform[0][2*j+1]);
			
			noisePower = noise*noise;
		
			realResult = mSignalF->waveform[1][2*j] * temp->waveform[1][2*j];
			imagResult = -mSignalF->waveform[1][2*j+1] * temp->waveform[1][2*j+1];
			
			power = realResult*realResult + imagResult*imagResult;
			
			SNR += (power/noisePower);
			
			

			op.waveform[1].push_back(realResult); //this should maybe be changed to SNR
			op.waveform[1].push_back(imagResult);
				
			
		}	
		
		
		mSNR.push_back(SNR);

		output->push_back(op);
	}

	mConResults = output;
	
	return;
}
void Extractor::fftInverse(std::vector<Template>* output)
{
	int I, N;
	vec_d amp; 
	vec_d freq;
	double dt;
	double norm;

	I = (*mConResults).size();	
	N = (*mConResults)[0].waveform[0].size();	

	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N/2);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N/2);
	
	for(int i=0; i<I; i++)
	{		
		freq = (*mConResults)[i].waveform[0];

		dt = 1/(freq[N/2-1]);
	
		amp = (*mConResults)[i].waveform[1];

		gsl_fft_complex_inverse(&amp[0], 1, N/2, compWT, complexWS);		

		Template temp;
		
		norm = fAutoCorrComplex(&mTemplates[0][i]);

		temp.param[0] = (*mConResults)[i].param[0];
		temp.param[1] = (*mConResults)[i].param[1];
 
		for(int n=0; n<N/2; n++)
		{
			temp.waveform[0].push_back(n*dt);
			temp.waveform[1].push_back(amp[2*n]/norm);
		}

		output->push_back(temp);
	}

	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	

	return;
}	
