#ifndef GWREADWRITE_H
#define GWREADWRITE_H

#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>

#include "gwDataTypes.h"

enum delimiter {tab, csv};

bool loadTemplates(std::string filename, std::vector<Template>* temps, delimiter delim)
{
	std::ifstream inFile;
	inFile.open(filename.c_str());

	char de;

	switch(delim)
	{
		case tab : de = '\t';
				   break;
	
		case csv : de = ',';
				   break;
	}

	if(inFile.fail())
	{
		std::cerr << "Input file not found" << std::endl;
		return false;
	}

	double d;
	std::string line, element;
	int n = 0;

	while(getline(inFile, line))
	{
		std::cout<<"Template "<<n<<std::endl;
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
			getline(inFile, line);

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
		n++;
	}

	inFile.close();

	return true;
}

bool saveTemplates(std::string filename, std::vector<Template>* temps, delimiter delim)
{
	std::ofstream outFile;
	outFile.open(filename.c_str());	

	if(outFile.fail())
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
		std::cout<<"Template "<<i<<std::endl;
		temp = templates[i];

		outFile << temp.param[0] << de << temp.param[1] << "\n";

		K = temp.waveform[0].size();

		for(int j=0; j<2; j++)
		{
			for(int k=0; k<(K-1); k++)
			{
				outFile << temp.waveform[j][k] << de;		
			}
			outFile << temp.waveform[j][K-1];
			outFile << "\n";
		}
	}

	outFile.close();

	return true;
}

bool loadSignals(std::string filename, std::vector<Signal>* sigs, delimiter delim)
{
	std::ifstream inFile;
	inFile.open(filename.c_str());

	char de;

	switch(delim)
	{
		case tab : de = '\t';
				   break;
	
		case csv : de = ',';
				   break;
	}

	if(inFile.fail())
	{
		std::cerr << "Input file not found" << std::endl;
		return false;
	}

	double d;
	std::string line, element;

	while(getline(inFile, line))
	{
		Signal sig;
		std::istringstream iss1(line);

		while(getline(iss1, element, de))
		{
	   	  	std::istringstream is(element);
		   	is >> d;
		  	sig.waveform[0].push_back(d);	
		}

		getline(inFile, line);

		std::istringstream iss2(line);

		while(getline(iss2, element, de))
		{
	   	  	std::istringstream is(element);
		   	is >> d;
		  	sig.waveform[1].push_back(d);	
		}
 
		sigs->push_back(sig);
	}

	inFile.close();

	return true;
}

bool saveSignals(std::string filename, std::vector<Signal>* sigs, delimiter delim)
{
	std::ofstream outFile;
	outFile.open(filename.c_str());	

	if(outFile.fail())
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

	std::vector<Signal> signals = *sigs;

	int I, K;
	Signal sig;

	I = signals.size();

	for(int i=0; i<I; i++)
	{
		sig = signals[i];

		K = sig.waveform[0].size();

		for(int j=0; j<2; j++)
		{
			for(int k=0; k<(K-1); k++)
			{
				outFile << sig.waveform[j][k] << de;		
			}
			outFile << sig.waveform[j][K-1];
			outFile << "\n";
		}
	}

	outFile.close();

	return true;
}

bool saveSNRs(std::string filename, std::vector<Template>* temps, delimiter delim)
{
	std::ofstream outFile;
	outFile.open(filename.c_str());	

	if(outFile.fail())
	{
		std::cerr << "Output file could not be opened" << std::endl;
		return false;
	}

	std::ofstream sampleFile;
	sampleFile.open("../data/sample.dat");
	if(sampleFile.fail())
	{
		std::cerr<<"Sample output file could not be opened"<<std::endl;
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

	//print 10 template outputs as a sample
	//as all templates give similar shapes, the first 10 outputs should be satisfactory to confirm expected shape
	for(int i=0; i<10; i++)
	{
		std::cout<<"Sample output "<<i<<std::endl;
		temp = templates[i];

		sampleFile << temp.param[0] << de << temp.param[1] << "\n";

		K = temp.waveform[0].size();

		for(int j=0; j<2; j++)
		{
			for(int k=0; k<(K-1); k++)
			{
				sampleFile << temp.waveform[j][k] << de;		
			}
			sampleFile << temp.waveform[j][K-1];
			sampleFile << "\n";
		}
	}

	std::cout<<"Sample outputs saved to sample.dat"<<std::endl;
	sampleFile.close();	

	for(int i=0; i<I; i++)
	{
		std::cout<<"SNR "<<i<<std::endl;
		temp = templates[i];

		std::vector<double>::iterator it=std::max_element(temp.waveform[1].begin(),temp.waveform[1].end());
		int maxIndex=std::distance(temp.waveform[1].begin(),it);

		outFile << temp.param[0] << de << temp.param[1] << de << temp.waveform[0][maxIndex] << de << temp.waveform[1][maxIndex] << "\n";

	}

	outFile.close();

	return true;

}

#endif //GWREADWRITE_H
