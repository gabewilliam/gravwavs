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

//computes the convolution of the data with the template
void convolution(std::vector<double> *f, std::vector<double> *g){ 


	size_t N = f->size();
	//if the data sets are different sizes then doesnt compute the result
	if (N != g->size()){
		std::cerr << "unmatched data size"<< std::endl;
	}else{
		//array for the convolution
		double h[N];
		
		
		for (size_t t = 0; t < N; t++){
	
			h[t] = 0;
			
			//multiplies the correct elements of the arrays together and adds the result ot the correct element of the convolution array
			//move the template to the left by equivalently moving the signal to the right
			//the only elements of f that are accesible are the elements that correspond to g(t=0) 
			//the number of accesible elements decreases as this is looped through
			for(size_t tau = 0; tau <= N-t; tau++){
				h[t] += (*f)[t+tau]*(*g)[tau];
			}
		}
		
		//saves the result to a file
		FILE *fout;
		fout = fopen("convolution.dat","w");
		FILE *fout2;
		fout2 = fopen("f.dat","w");
		FILE *fout3;
		fout3 = fopen("g.dat","w");
	
		for (size_t t = 0; t < N; t++){
			fprintf(fout, "%f,%lu\n",h[t],t);
			fprintf(fout2, "%f,%lu\n",(*f)[t],t);
			fprintf(fout3, "%f,%lu\n",(*g)[t],t);
		}
	
	}
	
	

}

int main(){


	//main is used to load the data
	std::fstream filter("hidden.txt", std::ios_base::in);
	
	std::vector<double> filterTime;
	std::vector<double> f;
	
	double data = 0;
	
	bool loop = true;
	
	
	
	while (loop){	
		if(filter >> data){
			filterTime.push_back(data);
			filter >> data;
			f.push_back(data);
		}else{
			loop = false;
		}			
	}

	std::fstream hidden("filter.txt", std::ios_base::in);
	
	std::vector<double> hiddenTime;
	std::vector<double> g;
	
	
	loop = true;

	while (loop){	
		if(hidden >> data){
			hiddenTime.push_back(data);
			hidden >> data;
			g.push_back(data);
		}else{
			loop = false;
		}			
	}

	convolution(&f,&g);	
	
	return 0;
}





