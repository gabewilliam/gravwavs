/*
 * program to find the convolution of two functions f(t) and g(t)
 * h(t) = [f*g](t).
 * g must be at the start of the array and the data must come after it.
 *
 *
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "ReadandPrint.h"//print function



//computes the convolution of the data with the template
void convolution(std::vector<double> *signal, std::vector<std::vector<double> > filters,std::vector<double> Time){ 
	
	std::vector<std::vector<double> > printer;//create 
	size_t NoFilters=filters.size();

	size_t N = signal->size();
	//if the data sets are different sizes then doesnt compute the result
	printer.push_back(Time);
		//array for the convolution
		double h[N];
		std::vector<std::vector<double> >Convolution;
		
	for(size_t i=0;i<NoFilters;++i){
			
		std::vector<double> *conv= new std::vector<double>;
		Convolution.push_back(*conv);
	
		for (size_t t = 0; t < N; t++){

			h[t] = 0;
				
			//multiplies the correct elements of the arrays together and adds the result ot the correct element of the convolution array
			//move the template to the left by equivalently moving the signal to the right
			//the only elements of f that are accesible are the elements that correspond to g(t=0)
			//the number of accesible elements decreases as this is looped through
			for(size_t tau = 0; tau <= N-t; tau++){
				h[t] += (*signal)[t+tau]*(filters[i])[tau];
					
			}
			Convolution[i].push_back(h[t]);
				
		}
		printer.push_back(Convolution[i]);
	}
	/*std::string filename;
	std::cout << "Enter Output file name..." << std::endl;
	std::cin >> filename;*/
	std::string filename="Output1.csv";

	outputWriter(printer,filename);
}

void frequencyDomainConvolution (std::vector<double> *signal, std::vector<std::vector<double> >filters,std::vector<double> signalTime){
	
	size_t N = (*signal).size();
	double stride = 1;
	size_t NoFilters=filters.size();
	std::vector<std::vector<double> > printer;
	
	gsl_fft_real_wavetable * real;
	gsl_fft_real_workspace * work;

	work = gsl_fft_real_workspace_alloc (N);
    real = gsl_fft_real_wavetable_alloc (N);

	gsl_fft_real_transform(&(*signal)[0], stride, N, real, work);
	  
	printer.push_back(signalTime);
	std::vector<double> filter;	
	
	gsl_fft_real_wavetable_free (real);
	gsl_fft_real_workspace_free (work);
	
	
	for(size_t r=0;r < NoFilters;r++){
		
		
		work = gsl_fft_real_workspace_alloc (N);
		real = gsl_fft_real_wavetable_alloc (N);
	    filter=filters[r];
	  
		gsl_fft_real_transform(&filter[0], stride, N, real, work); 
		gsl_fft_real_wavetable_free (real);
	
	//initialises the convolution vector as the complex conjugate
	  std::vector <double> *convolution=new std::vector<double>;
	  
		for(size_t i = 0; i < N/2; i++){ 
			convolution->push_back((*signal)[i]*filter[i]);
		}
		for(size_t i = N/2; i < N; i++){
			convolution->push_back(-(*signal)[i]*filter[i]);
		}
	
	//the inverse transform is taken to get back to the real convolution
		gsl_fft_halfcomplex_wavetable *hc;
		hc = gsl_fft_halfcomplex_wavetable_alloc (N);
		gsl_fft_halfcomplex_inverse (&(*convolution)[0], stride, N, hc, work);
		
		gsl_fft_halfcomplex_wavetable_free (hc);
		gsl_fft_real_workspace_free (work);
	
		
	printer.push_back(*convolution);
	
	
	}
	/*std::string filename;
	std::cout << "Enter Output file name..." << std::endl;
	std::cin >> filename;*/
	std::string filename="Output2.csv";
	outputWriter(printer, filename);
	
}


int main(){	//main is used to load the data


	//User input for filename


	//Setting up input stream
	std::ifstream in_file;
	in_file.open("signalA.dat");
	
	//Checking if stream is setup succesfully
	if(in_file.fail())
	{
		std::cout << "Input file not found" << std::endl;
		return false;
	}
	
	
	std::vector<double> signalTime;
	std::vector<double> signal;
	
	double a,b;

	while (!in_file.eof()){	

			in_file>>a>>b;
			signalTime.push_back(a);
			signal.push_back(b);
		
	}
	

	
	//Setting up input stream
	std::ifstream in_file2;
	in_file2.open("o.txt");
	
	//Checking if stream is setup succesfully
	if(in_file2.fail())
	{
		std::cout << "Input file not found" << std::endl;
		return false;
	}
	
	std::vector	<double> filterTime;
	std::vector<std::vector<double> >filters;
	
	int c=0;
	const size_t noofFilters=3;
	for(size_t l=0;l<noofFilters;l++){
			std::vector<double> *filter= new std::vector<double>;
			filters.push_back(*filter);
	}
	
	double d;
	std::cout<<"hello?";
	while (!in_file2.eof()){	

		for(size_t c=0;c<noofFilters;++c){
			in_file2>>d;		
			filters[c].push_back(d);
		}
	}
		
	for(size_t j=0;j<noofFilters;++j){
		while (filters[j].size()<signal.size()){
			filters[j].push_back(0.0);
		}
	}
	
	int pushcounter;
	/*
	while (filter.size()>signal.size()){
		signal.insert(signal.begin() , 0.0);
		++pushcounter;
	}
	*/
	//std::string filename="Output2.csv";
	//outputWriter(filters,filename);
	convolution(&signal,filters,signalTime);	
	frequencyDomainConvolution(&signal,filters,signalTime);
	
	return 0;

}
