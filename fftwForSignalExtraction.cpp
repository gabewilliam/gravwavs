
#define REAL(z, i)  (z[i][0])
#define IMAG(z, i) (z[i][1])//not redundant

//FFTW FUNCTIONS__________________________________

void Extractor::fftw(std::vector<Template>* templatesFFT){
	//number of templates
	int T = mTemplatesT->size();
	//size of templates
	int N;
	
	int size = mSignalT->waveform[0].size();
	int allocSize = 2*size - 1;
	int j = 0;
	double sw;
	

	//used to input into mTemplatesF
	vec_d freq;
	vec_d result;
	
	//templates for ease
	Template *tempT;
	Template tempF;

	fftw_complex *filterT;
	fftw_complex *filterF;
	
	filterT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	filterF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	
	fftw_plan forwardPlan = fftw_plan_dft_1d(allocSize, filterT, filterF, FFTW_FORWARD, FFTW_ESTIMATE);
	
	//loop through templates
	for (int i = 0; i < T; i++){
	
		//sets the current template	
		tempT = &mTemplatesT[0][i];
		//gets the size of the template
		N = tempT->waveform[1].size();
		
		//loops up the the size of the template and then fills the reast with 0 zeros up to 2 * signal size
		for (j = 0; j < N; j++){
			REAL(filterT, j) = tempT->waveform[1][j];
			IMAG(filterT, j) = 0.0;		
		}
		for (j = N; j < allocSize; j++){
			REAL(filterT, j) = 0.0;
			IMAG(filterT, j) = 0.0;			
		}
		
		fftw_execute(forwardPlan);
		
		//the frequency stuff is mostly wrong
		sw = N/(2.0*(tempT->waveform[0])[N-1]);
		
		for(j = 0; j < size; j++){
			if (j < size/2){
				freq.push_back(sw*j);
			}else{
				freq.push_back(-sw*(j - size/2));
			}		
		}
		//the size of this might have to be adjusted
		//should save the result as a vec_d packed as in the gsl complex routines		
		for (j = 0; j < allocSize; j++){
			result.push_back(REAL(filterF, j));
			result.push_back(IMAG(filterF, j));
		}
		
		//allocate each 
		tempF.param[0] = tempT->param[0];
		tempF.param[1] = tempT->param[1];
		tempF.waveform[0] = freq;
		tempF.waveform[1] = result;
		
		templatesFFT->push_back(tempF);
			
	}
	mTemplatesF = templatesFFT;
	
	fftw_free(filterT);
	fftw_free(filterF);
	fftw_destroy_plan(forwardPlan);

} 

void Extractor::fftw(Signal* signalFFT){

	vec_d freq;
	vec_d result;
	
	//the FFTW stuff
	fftw_complex *signalT;
	fftw_complex *signalF;
	
	int N = mSignalT->waveform[0].size();
	int allocSize = 2*N - 1;
	
	signalT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	signalF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	
	fftw_plan forwardPlan = fftw_plan_dft_1d(allocSize, signalT, signalF, FFTW_FORWARD, FFTW_ESTIMATE);
	
	//loops to the 2 * signal size -1
	//zero pads for the later half
	for (int i = 0; i < N; i++){
		REAL(signalT, i) = mSignalT->waveform[0][i];
		IMAG(signalT, i) = 0.0;
	}
	for (int i = N; i < allocSize; i++){
		REAL(signalT, i) = 0.0;
		IMAG(signalT, i) = 0.0;
	}
	
	fftw_execute(forwardPlan);
	
	//pack the vectors as was done for the gsl complex routine
	for (int i = 0; i < allocSize; i++){
		result.push_back(REAL(signalF, i));
		result.push_back(IMAG(signalF, i));
	}
	
	signalFFT->waveform[1] = result;

	fftw_free(signalT);
	fftw_free(signalF);
	fftw_destroy_plan(forwardPlan);

}


//THIS FUNCTION IS NOT SCALED SO THE Y AXIS WILL NOT BE CORRECT
void Extractor::fftwInverse(std::vector<Signal>* signalsFFTI){
	
	//number of templates
	int N = signalsFFTI->size();	
	
	//gets the time of the initial signal
	vec_d time;
	time = mSignalT->waveform[0];
	vec_d result;
	
	int j = 0;
	
	//gets the size of each of the signals
	int allocSize = (*signalsFFTI)[0].waveform[1].size()/2;
	
	//FFTW stuff
	fftw_complex *convoT;
	fftw_complex *convoF;	
	
	convoT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	convoF = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*allocSize);
	
	fftw_plan backwardPlan = fftw_plan_dft_1d(allocSize, convoF, convoT, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	for (int i = 0; i < N; i++){
		//unpacks the complex vectors from the initial routines
		for (j = 0; j < allocSize; j++){
			REAL(convoF, j) = (*signalsFFTI)[i].waveform[1][2*j];
			IMAG(convoF, j) = (*signalsFFTI)[i].waveform[1][2*j + 1];
		}
		fftw_execute(backwardPlan);
		
		for (j = 0; j < allocSize; j++){
			result.push_back(REAL(convoT, j));
		}
		//change the input to the time domain
		(*signalsFFTI)[i].waveform[0] = time;
		(*signalsFFTI)[i].waveform[1] = result;

	}
	

	fftw_free(convoT);
	fftw_free(convoF);
	fftw_destroy_plan(backwardPlan);

}

