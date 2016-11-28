// Computational Methods Assignment 4 - Numerical Integration
// Martik Aghajanian, Cohort 8
//
// Program which estimates the value of the error function, which
// is the integral from 0 to argument x of the error function 
// integrand erfi(x), using both the trapezoid rule and Simpson's
// rule to compare the number of necessary function evaluations.
#include <iostream>
#include <cmath>
#include "integrate.h"
using namespace std;

double erfi(double x) {
	// Given double-preicision float "x", returns the value of the
	// error function integrand (erfi), using 4*arctan(1) as pi.
	// This simple function was not declared separately as it is not 
	// detrimental to the readability of the programs general purpose.
        return (2/sqrt(atan(1.0)*4))*exp(-x*x);
}

int main() {
	// First this program instantiates two "Integrate" objects which are identical
	// in setup, but progress using separate methods, one using the trapezoid rule
	// and one using Simpson's rule, each with function "erfi", a precision 1e-6
	// and limits from 0 to 2.
        Integrate Result1(&erfi, 0, 2, 1e-6);
	Integrate Result2(&erfi, 0, 2, 1e-6);

	// Both values of the integral using the separate methods are displayed, each 
	// with their respective number of function evaluations for comparison.
	cout << "Trapezium rule: " << Result1.trap();
	cout <<  " using " << Result1.get_evals() << " evaluations." << endl;
	cout << "Simpsons's rule: " << Result2.simp();
	cout << " using " << Result2.get_evals() << " evaluations." << endl;
	return 0;
}
