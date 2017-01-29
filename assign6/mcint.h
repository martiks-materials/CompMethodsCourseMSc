// Computational Methods Assignment 6 - Monte Carlo Integrations
// Martik Aghajanian
//
// Contains the class declaration of the Monte Carlo integrator
// that is used in this assignment. 
#include <iostream>
#include <vector>
#include <cmath>
#include "rand.h"
using namespace std;


// These typedefs allow for less messy declaration in this program.
typedef vector<double> vec;
typedef double (*func)(vector<double>);
typedef double (*unifun)(double);

class Mcint {
// This class takes in an arbitraty N-variable function and a vector of limits
// to integrate using Monte Carlo Methods to a certain accuracy. This can be done 
// using an arbitrary probability function to implement importance sampling. The 
// means of achieving the arbitrarily distributed random variables can be specified
// by the user. 
private:	
	func f, p; 
	unifun cy;
	int Nd;
	vec x_up, x_low;;
	double eps;
	Ran myran;  
public:
	// To avoid overflow, the number of steps M must be a long unsigned integer
	// and the initial number of steps M_init must also be since the variables 
	// are compared.
	vec steps, errors, vals, X;
	unsigned long int M, M_init, M_max;
	double R, V, fsum_old, fsum_new, fsq_old, fsq_new; 	
	bool importance, reject;
	Mcint(func fnc, vec a, vec b, int N, int seed, double tol, int Mstart, 
	      bool imp, func pdf, unifun cdf, unsigned long int stepmax, bool rej=false);
	
	void sample();

	double error1();

	double error2();

	void integrate();

	void generator();
};
			
					
			
		
		
		



	
