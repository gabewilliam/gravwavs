#include "Binary.cpp"
#include <iostream>


int main() {

	Binary TomAndJerry (75.0, 50.0, 17);

	TomAndJerry.printGets();	

	std::cout<<std::endl<<"change m1, m2 and a:"<< std::endl;
	double m[2], a;
	std::cin>>m[0]>>m[1]>>a;

	for(int n=1; n<3; n++){

		TomAndJerry.setMass(n, m[n-1]);

	}

	TomAndJerry.setSeparation(a);

	TomAndJerry.printGets();

	return 0;

}
