#ifndef GWREADWRITE_H
#define GWREADWRITE_H

#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <zlib.h>

#include "gwDataTypes.h"

#define CHUNK 16384
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
	FILE * outFile=fopen(filename.c_str(),"wb");

	std::string inStr;
	char* element;
	std::vector<Signal> signals=*sigs;

	char de;

	switch(delim)
	{
		case tab : de = '\t';
				   break;
	
		case csv : de = ',';
				   break;
	}

	//put signals into a string for output
	for (int i=0; i<signals.size(); i++)
	{
		int J=signals[i].waveform[0].size();
		for (int j=0; j<J-1; j++)
		{
			sprintf(element,"%f,",signals[i].waveform[0][j]);
			inStr.append(element);
		}

		sprintf(element,"%f\r\n",signals[i].waveform[0][J-1]);
		inStr.append(element);
		for (int j=0; j<J-1; j++)
		{
			sprintf(element,"%f,",signals[i].waveform[1][j]);
			inStr.append(element);
		}
		sprintf(element,"%f\r\n",signals[i].waveform[1][J-1]);
		inStr.append(element);
		
		std::cout<<"test "<<i<<std::endl;
	}

	//lifting from zlib docs begins here
	int ret, flush;
	unsigned have;
	z_stream strm;
	unsigned char in[CHUNK];
	unsigned char out[CHUNK];

	//allocate deflate state
	strm.zalloc=Z_NULL;
	strm.zfree=Z_NULL;
	strm.opaque=Z_NULL;
	ret=deflateInit(&strm,Z_DEFAULT_COMPRESSION);
	if (ret!=Z_OK)
	{
		std::cout<<"Error in ret initialisation"<<std::endl;
		return false;
	}

	std::stringstream ss(inStr);
	int k=0;
	//compress until end of input - create an input!!
	do
	{
		char chunk[CHUNK];
		ss.read(chunk,CHUNK);
		if(ss.eof())
		{
			std::cout<<"eof"<<std::endl;		
		}
		else if(ss.fail())
		{
			std::cout<<"fail"<<std::endl;		
		}
		else if(ss.bad())
		{
			std::cout<<"bad"<<std::endl;		
		}
		strm.avail_in=ss.gcount();
		flush=ss.eof() ? Z_FINISH : Z_NO_FLUSH;
		strm.next_in=(unsigned char*)chunk;
		//run deflate on input until output buffer not full, finish compression if all of source has been read in
		do
		{
			strm.avail_out=CHUNK;
			std::cout<<out<<" ";
			strm.next_out=out;
			ret=deflate(&strm,flush); //no bad return value
			std::cout<<strm.next_out<<std::endl;
			if (ret==Z_STREAM_ERROR) //state not clobbered
			{
				std::cout<<"Error: Z_STREAM_ERROR"<<std::endl;
				return false;
			}
			have=CHUNK-strm.avail_out;
			if(fwrite(out,1,have,outFile)!=have || ferror(outFile))
			{
				(void)deflateEnd(&strm);
				std::cout<<"Error: "<<Z_ERRNO<<std::endl;
				return false;
			}
		}
		while (strm.avail_out==0);
		if(strm.avail_in!=0)
		{
			std::cout<<"Error: strm.avail_in not zero"<<std::endl;
			return false;
		}
		std::cout<<k<<std::endl;
		k++;
	}
	while (flush!=Z_FINISH);

	
	if (ret!=Z_STREAM_END)
	{
		std::cout<<"Error: ret!=Z_STREAM_END"<<std::endl;
		return false;
	}

	(void)deflateEnd(&strm);

	return true;

	/*
	//filename must end in .gz for this to work
    gzFile outFile=gzopen(filename.c_str(),"wb");    

	//if(outFile.fail())
	//{
	//	std::cerr << "Output file could not be opened" << std::endl;
	//	return false;
	//}

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
	*/
}

bool saveTemplatesCompressed(std::string filename, std::vector<Template>* temps, delimiter delim)
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

	std::vector<Template> templates = *temps;

	int I, K;
	Template temp;

	I = templates.size();

	for(int i=0; i<I; i++)
	{
		temp = templates[i];

		K = temp.waveform[0].size();

		gzprintf(outFile,"%g%c %g\r\n",temp.param[0],de,temp.param[1]);

		for(int j=0; j<2; j++)
		{
			for(int k=0; k<(K-1); k++)
			{
				gzprintf(outFile,"%g%c",temp.waveform[j][k],de);				
				//outFile << sig.waveform[j][k] << de;		
			}
			gzprintf(outFile,"%g\r\n",temp.waveform[j][K-1]);
			//outFile << sig.waveform[j][K-1];
			//outFile << "\r\n";
		}
	}

	gzclose(outFile);

	return true;
}

//source for the below: zlib documentation - Usage example
bool loadSignalsCompressed(std::string filename, std::vector<Signal>* sigs, delimiter delim)
{
	//check that file in correct format
	/*bool gzipped=endsWith(filename, ".gz");
	if(!gzipped)
	{
		std::cout<<"Incorrect input file format"<<std::endl;
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

	FILE * inFile=fopen(filename.c_str(),"r");    
	std::string inStr;

	//THIS SECTION TAKEN FROM ZLIB DOCUMENTATION
	//local variables
	int ret;
	unsigned have;
	z_stream strm;
	unsigned char in[CHUNK];
	unsigned char out[CHUNK];

	//allocate inflate state
	strm.zalloc=Z_NULL;
	strm.zfree=Z_NULL;
	strm.opaque=Z_NULL;
	strm.avail_in=0;
	strm.next_in=Z_NULL;
	ret=inflateInit2(&strm,16+MAX_WBITS);
	if (ret!=Z_OK)
	{
		std::cout<<"Loading failed: ret"<<std::endl;
		return false;
	}
	int i=0;
	//decompress until deflate stream ends or end of file
	do
	{
		int j=0;
		strm.avail_in=fread(in,1,CHUNK,inFile);
		if (ferror(inFile))
		{
			(void)inflateEnd(&strm);
			std::cout<<"Loading failed: error in input file"<<std::endl;
			return false;
		}
		if (strm.avail_in==0)
			break;
		strm.next_in=in;
		//run inflate on input until output buffer not full
		do
		{
			std::cout<<j<<std::endl;
			strm.avail_out=CHUNK;
			strm.next_out=out;
			ret=inflate(&strm, Z_NO_FLUSH);
			if (ret==Z_STREAM_ERROR)
			{
				std::cout<<"ret: Z_STREAM_ERROR"<<std::endl;
				return false;
			}
			switch(ret)
			{
				case Z_NEED_DICT:
					ret=Z_DATA_ERROR;
					(void)inflateEnd(&strm);
					std::cout<<"Z_NEED_DICT"<<std::endl;
					return false;
				case Z_DATA_ERROR:
					(void)inflateEnd(&strm);
					std::cout<<"Z_DATA_ERROR"<<std::endl;
					return false;
				case Z_MEM_ERROR:
					(void)inflateEnd(&strm);
					std::cout<<"ret: Z_NEED_DICT, Z_DATA_ERROR or Z_MEM_ERROR"<<std::endl;
					return false;
			}
			have=CHUNK-strm.avail_out;

			//writing to string
			inStr.append(reinterpret_cast<const char*>(out));
			std::cout<<inStr;
			j++;

			/*
			if(fwrite(out,1,have,dest)!=have || ferror(dest))
			{
				(void)inflateEnd(&strm);
				std::cout<<"Error in loading to variable"<<std::endl;
				return false;
			}
			*/
		}
		while(strm.avail_out==0);
		std::cout<<"i "<<i<<std::endl;
		i++;
	}
	while (ret!=Z_STREAM_END);

	//clean up
	(void)inflateEnd(&strm);
	return ret==Z_STREAM_END ? true : false;

	//END ZLIB DOCUMENTATION LIFTING
	std::cout<<inStr[12]<<std::endl;
	double d;
	std::string line, element;
	std::istringstream ss(inStr);

	while(getline(ss,line,'\n'))
	{
		Signal sig;
		std::istringstream iss1(line);

		while(getline(iss1, element, de))
		{
	   	  	std::istringstream is(element);
		   	is >> d;
		  	sig.waveform[0].push_back(d);	
		}

		getline(ss,line,'\n');

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


#endif //GWREADWRITE_H
