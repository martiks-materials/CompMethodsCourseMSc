// Computational Methods Assignment 3 - Random Numbers - Program 2
// Martik Aghajanian
// 
// This program uses the random number generator structure Ran to generate
// samples of a uniform deviate which are transformed into a random variable
// y which follows a probability distribution function 0.5*sin(y). This is
// outputted to a file to plot on a frequency distribution histogram and is
// also timed to compare to the rejection method.


#include <iostream>
#include <fstream> 
#include <cmath>
#include <ctime>
#include "rand.h"
using namespace std;

int main() {
	// Here the C++ struct Ran assigned to variable "myran" is initialized
        // with seed 42, and the number of samples required, "lim" is specified.
        int lim(pow(10,5));
        double random[lim] = {};

	// The timer is started here as this is where the random number generator
        // is still seeded and no unnecessary loops or conditionals are present.
	time_t start = clock();
	Ran myran(42);

	// For each recurrence of the loop, the member function doub() of Ran is
        // called which both modifies the three integers of "myran" that produce
        // the random numbers, and outputs a new double-precision floating point
        // between 0.0 and 1.0. This is then fed into the appropriate function
	// that transforms it into random variable y from distribution 0.5*sin(y).
	// This is then stored in the double-type array "random".
        for(int n(0); n<lim ;n++) {
                double X = myran.doub();
		random[n] = acos(1-(2*X));
        }
	
	// The timer is stopped here, as all random number generation in the program
	// has ceased. The time difference has to be divided by the "clock time" to
	// obtain it in seconds. The data from the array "random" is then outputted to
	// an external data file for plotting.
	time_t end = clock();
	cout << "Transformation Method Runtime: " << (end-start)/(double)CLOCKS_PER_SEC << "s for " << lim << " samples." << endl;
        ofstream outfile1("transform.dat");
        if ( !outfile1.is_open() ){
     	        cout << "Error opening file." << endl;
                return 1;
        }
        for(int n(0); n<lim; n++){
                outfile1 << random[n] << endl;
        }
        outfile1.close();
	return 0;
}

	
