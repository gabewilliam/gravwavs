//for putting in gwSignal...h
double fAutoCorrComplex(Template* temp);





double Extractor::fAutoCorrComplex(Template* temp){
	int K =temp->waveform[0].size();//declare number of data points
	Signal op;//redeclare the signal vector
	double result;

		for(int k=0; k<K; k++){
			//calculate and push back the data with a waveform (multiply template by itself in f domain)  
			result = temp->waveform[1][2*k] * temp->waveform[1][2*k];	
			op.waveform[0].push_back(k);
			op.waveform[1].push_back(result);
			//and its conjugate if it is imaginary
			result = -(temp->waveform[1][2*k+1] * temp->waveform[1][2*k+1]);
			op.waveform[0].push_back(k);
			op.waveform[1].push_back(result);
		
		}
	size_t N=K;

	//perform inverse fft to find correlation between template and itself at 0 displacement
	gsl_fft_complex_workspace* complexWS = gsl_fft_complex_workspace_alloc(N);
	gsl_fft_complex_wavetable* compWT = gsl_fft_complex_wavetable_alloc(N);

	gsl_fft_complex_inverse(&(op.waveform[1][0]), 1, N, compWT, complexWS);

	gsl_fft_complex_workspace_free(complexWS);
	gsl_fft_complex_wavetable_free(compWT);	
	
	//output the first element of the autocorrelation
	return op.waveform[1][0];
}