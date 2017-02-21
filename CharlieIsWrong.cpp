#include <iostream>
#include <cstdlib>

using namespace std;

bool CharlieIsWrong(){
	return true;
}

int main(){
	
	if (CharlieIsWrong){
		cout << "Obviously Charlie is wrong" << endl;
	}
	else{
		cout << "Program error" << endl;		
	}
		
	return 0;
	
}