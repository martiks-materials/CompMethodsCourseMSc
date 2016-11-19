// Computational Methods Assignment 3 - Random Numbers - Program 3
// Martik Aghajanian
//
// This program samples a uniform deviate x from a random number generator
// structure Ran and uses the Rejection method to obtain from it a random
// variable y which comes from the probability distribution (2/pi)(sin(y))^2
// with a comparision function of (2/pi)sin(y). This is timed for comparison
// with the transformation method.
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include "rand.h"
#include "compare.h"
using namespace std;

int main() {
	// Here the C++ struct Ran assigned to variable "myran" is initialized
        // with seed 42, and the number of samples required, "lim" is specified.
        int lim(pow(10,5));
  	double pi = 4*atan(1.0);
        double random[lim] = {};
	double A = 4.0/pi;	  
  
        // The timer is started here as this is where the random number generator
        // is still seeded and no unnecessary loops or conditionals are present		
	clock_t start = clock();
	Ran myran(42);
	
	// For each recurrence of the loop, the member function doub() of Ran is
        // called which both modifies the three integers of "myran" that produce
        // the random numbers, and outputs a new double-precision floating point
        // between 0.0 and 1.0. This is multiplied by the size of the range over
	// which the cummulative distribution function of the comparison function
	// can take value of. This value is inverted to get Y, and the doub()
	// function produces another random number to decided whether or not the
	// value Y should be rejected (using the conditional statement below).
	// if successful, the non-rejected Y value is stored in array "random".
        for(int n(0); n<lim ;n++) {
		bool reject(true);
                while(reject){
			double X = myran.doub()*A;
			double Y = invcdf(X);
			double Z = myran.doub()*comp(Y);
			if(Z<=pdf(Y)){
				reject=false;
				random[n] = Y;
			}
     		}
        }

	// The timer is stopped here, as all random number generation in the program,
	// and the rejection conditionals  have ceased. The time difference has to be 
	// divided by the "clock time" to  obtain it in seconds. The data from the array 
	// "random" is then outputted to  an external data file for plotting.
	clock_t end = clock();
	cout << "Rejection Method Runtine: " << (end-start)/(double)CLOCKS_PER_SEC << "s for " << lim << " samples." << endl;
        ofstream outfile1("reject.dat");
        if ( !outfile1.is_open() ){
                cout << "Error opening file." << endl;
                return 1;
        }
        for(int n(0); n<lim; n++){
                outfile1 << random[n] <<  endl;
        }
        outfile1.close();
        return 0;
}	
