// Computational Methods Assignment 6 - Monte Carlo Integrations
// Martik Aghajanian
//
// Contains the class declaration of the Monte Carlo integrator
// that is used in assignment 6. 
#include <iostream>
#include <vector>
#include <cmath>
#include "rand.h"
using namespace std;


// These typedefs allow for less messy declaration in this program.
typedef vector<double> vec;
typedef double (*func)(vector<double>);
typedef double (*unifun)(double);
typedef vector< vector<double> > vecvec;

class Mcint {
// This class takes in an arbitraty N-variable function and a vector of limits
// to integrate using Monte Carlo Methods to a certain accuracy. This can be done 
// using an arbitrary probability function to implement importance sampling. The 
// means of achieving the arbitrarily distributed random variables can be specified
// by the user. 
private:	
	// The function to be integrated 'f' and the probability distribution that
	// the random variables are picked from (if importance sampling is used), are
	// specified by the user. For the case of a single independent variable in the 
	// integrand, an inverse cummulative probability distribution 'cy' corresponding
	// to 'p' can be specified so that the transformation method can be used. The
	// number of independent variables in the function is gives the dimension of the
	// integration 'Nd'. The user must also specify bounds in the form of vectors
	// 'x_up' and 'x_low' which give upper and lower bounds respectively. The tolerance
	// or accuracy desired by the user is 'eps'. The random numbers for this can be used
	// from the uniform random generator 'myran' and transformed using 'cy', or can be 
	// supplied by the user using a vector of vectors 'disvec' which will contain a number
	// of Nd-component vectors equal to the maximum number of samples in this integration.
	// This means for Gaussian distributions and non-invertible many dimension integrals, 
	// supplying the random numbers will be the users task, though the corresponding 
	// probability distribution must also be supplied.

	func f, p; 
	unifun cy;
	int Nd;
	vec x_up, x_low;;
	double eps;
	Ran myran;
	vecvec disvec;  
public:
	// The vectors 'steps', 'errors' and 'vals' are the data structures where the number
	// of iterations of the method, the errors at each iteration, and the value of the
	// integral at each iteration respectively. The updatable vector 'X' is used to store
	// the random numbers in to pass to the function. The number of samples in the 
	// integral is 'M', whilst 'M_init' gives the number of samples used initially in
	// the first iteration and 'M_max' puts a constraint on how many samples can be added
	// before the MC integrator will cease. To avoid overflow, the number of steps M must 
	// be a long unsigned integer as must be M_init and M_max must also be since there are
	// functions which compare these variables. The bool 'importance' acts as the switch 
	// to implement importance sampling, whilst 'transform' (defaulted to true) is the 
	// switch for the transformation method, which if false will instead rely on the user-
	// supplied vector of vectors containing arbitrariy distributed random numbers.

	vec steps, errors, vals, X;
	unsigned long int M, M_init, M_max;
	double R, V, fsum, fsq; 	
	bool importance, transform;
	Mcint(func fnc, vec a, vec b, int N, int seed, double tol, int Mstart, 
	      bool imp, func pdf, unifun cdf, unsigned long int stepmax, bool transform=true, vecvec distri={});
	
	void sample();

	double error();

	void integrate();

	void generator();
};
			
					
			
		
		
		



	
