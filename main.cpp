#include "Binary.cpp"
#include <iostream>
//so far this tests a new printGets function (itself a test of all the gets) and that the sets work

int main() {

	Binary TomAndJerry (75.0, 50.0, 17);//create a sample Binary

	TomAndJerry.printGets();//use printGets() to cout all the Binary::get...() methods

	std::cout<<std::endl<<"change m1, m2 and a:"<< std::endl;//can change any of the fundemental properties of binary
	double m[2], a;
	std::cin>>m[0]>>m[1]>>a;

	for(int n=1; n<3; n++){

		TomAndJerry.setMass(n, m[n-1]);

	}

	TomAndJerry.setSeparation(a);

	TomAndJerry.printGets();//checking that Binary::set...() methods actually worked
	
	//next time... test evolutions & CHE checks, maybe test building a std::vector of Binaries

	return 0;

}
