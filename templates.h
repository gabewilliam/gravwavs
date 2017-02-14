#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

typedef std::vector<double> vec_d;

enum delim {tab, csv};

//Template Struct to handle template waveforms & parameters
struct Template 
{
	vec_d waveform; 
	double param[2];
};

bool load_templates(std::string filename, std::vector<Template>* temps, delim de)
{
	std::ifstream in_file;
	in_file.open(filename.c_str());	

	if(in_file.fail())
	{
		std::cerr << "Input file not found" << std::endl;
		return false;
	}

	double d;
	std::string line, element;
	std::string::size_type st;

	while(getline(in_file, line))
	{
		Template temp;

		std::istringstream iss(line);

		for(int i=0; i<2; i++)
		{
			switch(de)
			{
				case tab : getline(iss, element, '\t');
						   break;
			
				case csv : getline(iss, element, ',');
						   break;
			}

			std::istringstream is(element);
			is >> d;
			temp.param[i] = d;
		}

		switch(de)
		{
			case tab : while(getline(iss, element, '\t'))
					   {
						   std::istringstream is(element);
						   is >> d;
						   temp.waveform.push_back(d);	
					   };
					   break;
		
			case csv : while(getline(iss, element, ','))
					   {
						   std::istringstream is(element);
						   is >> d;
						   temp.waveform.push_back(d);	
					   };
					   break;
		}

		temps->push_back(temp);
	}

	return true;
}

bool save_templates(std::string filename, std::vector<Template>* temps, delim de)
{
	std::ofstream out_file;
	out_file.open(filename.c_str());	

	if(out_file.fail())
	{
		std::cerr << "Output file could not be opened" << std::endl;
		return false;
	}

	int I, J;
	Template temp;

	std::vector<Template> templates = *temps;

	I = templates.size();

	for(int i=0; i<I; i++)
	{
		temp = templates[i];

		switch(de)
		{
			case tab : out_file << temp.param[0] << "\t" << temp.param[1];
					   break;
			
			case csv : out_file << temp.param[0] << "," << temp.param[1];
					   break;
		}

		J = temp.waveform.size();

		for(int j=0; j<J; j++)
		{
			switch(de)
			{
				case tab : out_file << "\t" << temp.waveform[j];
						   break;
			
				case csv : out_file << "," << temp.waveform[j];
						   break;
			}
			
		}

		out_file << std::endl;
	}

	return true;
}


