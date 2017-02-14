#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

typedef std::vector<double> vec_d;

enum delimiter {tab, csv};

//Template Struct to handle template waveforms & parameters
struct Template 
{
	vec_d waveform[2]; 
	double param[2];
};

bool load_templates(std::string filename, std::vector<Template>* temps, delimiter delim)
{
	std::ifstream in_file;
	in_file.open(filename.c_str());

	char de;

	switch(delim)
	{
		case tab : de = '\t';
				   break;
	
		case csv : de = ',';
				   break;
	}

	if(in_file.fail())
	{
		std::cerr << "Input file not found" << std::endl;
		return false;
	}

	double d;
	std::string line, element;

	while(getline(in_file, line))
	{
		Template temp;
		std::istringstream iss(line);

		//loading parametres from first line of template
		for(int i=0; i<2; i++)
		{
			getline(iss, element, de);

			std::istringstream is(element);
			is >> d;
			temp.param[i] = d;
		}

		//looping over the next two lines to extract the time and waveform data
		for(int i=0; i<2; i++)
		{
			getline(in_file, line);

			std::istringstream iss(line);
			
			//looping until the end of the line whilst extracting elements
			while(getline(iss, element, de))
			{
		   	  	std::istringstream is(element);
			   	is >> d;
			  	temp.waveform[i].push_back(d);	
			}
	  	}

		temps->push_back(temp);
	}

	return true;
}

bool save_templates(std::string filename, std::vector<Template>* temps, delimiter delim)
{
	std::ofstream out_file;
	out_file.open(filename.c_str());	

	if(out_file.fail())
	{
		std::cerr << "Output file could not be opened" << std::endl;
		return false;
	}

	char de;

	switch(delim)
	{
		case tab : de = '\t';
				   break;
	
		case csv : de = ',';
				   break;
	}

	std::vector<Template> templates = *temps;

	int I, K;
	Template temp;

	I = templates.size();

	for(int i=0; i<I; i++)
	{
		temp = templates[i];

		out_file << temp.param[0] << de << temp.param[1] << "\n";


		K = temp.waveform[0].size();

		for(int j=0; j<2; j++)
		{
			for(int k=0; k<K; k++)
			{
				out_file << temp.waveform[j][k] << de;		
			}
			
			out_file << "\n";
		}
	}

	return true;
}


