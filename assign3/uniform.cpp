#include <iostream>
#include <fstream>
#include <cmath>
#include "rand.h"
using namespace std;



int main(){
	Ran myran(42);
	int lim(pow(10,5));
	double random[lim] = {};
	for(int n(0); n<lim ;n++) {
		double X = myran.doub();
		random[n] = X;
	}
	ofstream outfile1("uniform.dat");
	if ( ! outfile1.is_open() ) {
    	    	cout << "Error opening file." << endl;
        	return 1;
    	}	
	for(int n(0); n<lim; n++){
		outfile1 << random[n] << endl;
	}
	outfile1.close();
	return 0;
}

