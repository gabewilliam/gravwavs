#include "templates.h"

int main()
{
	std::vector<Template> temps;

	if(!load_templates("templates.txt", &temps, csv))
	{
		std::cout << "Program exiting..." << std::endl;
		return 0;
	}

	std::cout << temps.size() << std::endl;

	if(!save_templates("templates.dat", &temps, tab))
	{
		std::cout << "Program exiting..." << std::endl;
		return 0;
	}

	return 1;
}
