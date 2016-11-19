// Computational Methods Assignment 3 - Random Numbers
// Martik Aghajanian
//
// This is a program that calls the random number generator structure
// called Ran, and produces N=10^5 samples of a uniform deviate, from 
// this generator. This is then outputted into a data file to be 
// plotted as a histogram with 100 bins.
#include <iostream>
#include <fstream>
#include <cmath>
#include "rand.h"
using namespace std;

int main(){
	// Here the C++ struct Ran assigned to variable "myran" is initialized
	// with seed 42, and the number of samples required, "lim" is specified.
	Ran myran(42);
	int lim(pow(10,5));
	double random[lim] = {};

	// For each recurrence of the loop, the member function doub() of Ran is
	// called which both modifies the three integers of "myran" that produce
	// the random numbers, and outputs a new double-precision floating point
	// between 0.0 and 1.0. This is stored in a double-type array "random".
	for(int n(0); n<lim ;n++) {
		double X = myran.doub();
		random[n] = X;
	}

	// This section of the program outputs the array of random numbers to 
	// an external data file.
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

