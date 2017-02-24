#ifndef GWREADWRITE_H
#define GWREADWRITE_H

#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
//#include <boost/iostreams/filtering_stream.hpp>    
//#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/copy.hpp>
//#include <boost/iostreams/filter/gzip.hpp>
#include <zlib.h>

#include "gwDataTypes.h"

//namespace bo=boost::iostreams;

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

	while(getline(inFile, line))
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
	}

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
			outFile << "\r\n";
		}
	}

	return true;
}

bool saveSignalsCompressed(std::string filename, std::vector<Signal>* sigs, delimiter delim)
{
	//filename must end in .gz for this to work
    gzFile outFile=gzopen(filename.c_str(),"wb");    

	/*if(outFile.fail())
	{
		std::cerr << "Output file could not be opened" << std::endl;
		return false;
	}*/

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
				gzprintf(outFile,"%g%c",sig.waveform[j][k],de);				
				//outFile << sig.waveform[j][k] << de;		
			}
			gzprintf(outFile,"%g\r\n",sig.waveform[j][K-1]);
			//outFile << sig.waveform[j][K-1];
			//outFile << "\r\n";
		}
	}

	gzclose(outFile);

	return true;
}

#endif //GWREADWRITE_H
