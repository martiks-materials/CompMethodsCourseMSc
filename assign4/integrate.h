// Computational Methods Assignment 4 - Numerical Integration
// Martik Aghajanian Cohort 8
//
// Class declaration for integration class which takes a user-specified function
// as well as limits and a user-specified precision to which the method of 
// integration is expected to converge. This can be used to integrate functions
// defined on an interval,  within a specified accuracy.
#ifndef INTEGRATE_H
#define INTEGRATE_H
#include <cmath>
using namespace std;

class Integrate{
// Class for integrating function "fnc" within limits "a", "b" within given precision "eps"
private:
	// For a function "fnc" which takes argument "x", with limits "a" and "b"
	// iteration the number of points "n" used in the calculation 
	// in addition to the boundary points. The convergence of the integral is
	// assessed using double precision float parameter "eps"
	double (*fnc)(double x);
	double a;
	double b;  	
	int n;
	double eps;

	// This class continually updates the variables "n" and "evals", the number
	// function evaluations made in total, upon calling the member function 
	// "iterate". The variables "I1" and "I2" are updatable values of the integral
	// values which are compared to assess convergence. The samples taken to 
	// evaluate the integral are spaced with spacing "h" which is halved upon
	// each iteration of doubling the number of points used.
	int evals;
	double h;
	double I2;
	double I1;
public:
	// See associated integrate.cpp file containing member functions for descriptions.
	Integrate(double (*func)(double x), double low, double up, double epsilon);

	void iterate();
	
	double trap();

	double simp();
		
	// This member function is a simple inline function and thus does not need to
	// be placed in an external .cpp file. This outputs the number of function evalulations 
	// required to determine the integral of "fnc" over interval between "a" and
	// "b" using either the trapeziod rule or Simpson's rule.
	int get_evals() { return evals; }
};

#endif
