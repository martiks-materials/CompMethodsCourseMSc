#include <iostream>
#include <vector>
#include <cmath>
#include "rand.h"
#include "mcint.h"
using namespace std;

// These typedefs allow for less messy declaration in this program.
typedef vector<double> vec;
typedef double (*func)(vector<double>);
typedef double (*unifun)(double);


Mcint::Mcint(func fnc, vec a, vec b, int N, int seed, double tol, int Mstart, 
      bool imp, func pdf, unifun cdf, unsigned long int stepmax){
	f = fnc;
	p = pdf;
	cy = cdf;
	x_low = a;
	x_up = b;
	Nd = N;
	// For uniform sampling, the volume is needed, so is constructed using
	// the bounds of the 'x_low' and 'x_up' for each dimension.
	V=1;
	for(int i(0); i<Nd; i++){
		V *= x_up[i] - x_low[i];
		X.push_back(0);
	}

	myran.seed(seed);
	eps = tol;
	M = 0;
	M_init = Mstart;
	M_max = stepmax;
	importance = imp;
	fsum = 0;
	fsq = 0;
	R = 1;
	// Intialises the first integral value by collecting 'M_init' samples
	// and estimating the sum of value and the sum of squared values.
	while(M < M_init){
		generator();
		double quant = (*f)(X)/((importance)?((*p)(X)):(1./V));
		fsum += quant;
		fsq += quant*quant;	
		M++;
	}
	
	// To get averages, they must be divided by the current number of samples
	fsum /= M;
	fsq /= M;
	steps.push_back(M);
	errors.push_back(R);
	vals.push_back(fsum); 
}

void Mcint::sample(){
	// Takes a sample of values of the function, with the amount of samples being
	// taken equal to the number of samples currently used in 'fsum' and 'fsquare',
	// such that the total sample size doubles. This creates dummies of the current
	// variables, modifies them and updates them so that none of the member variables
	// are accidentally changed during the samples (e.g. M_old needs to be fixed)
	// The temporary sum and sum of squares must me multiplied by current steps.
	double fsum_temp = M*fsum, fsq_temp = M*fsq ;
	int M_old = M;
	for(int i(0); i<M_old; i++){
		// The generator is called which updates X and depending on whether or not
		// importance sampling is being implemented, the quantities 'quant' will 
		// give the function value for this random X or the function divided by the
		// corresponding probability distribution. For each sample, the sum, sum of
		// squares, and the number of samples 'M' are updated accordingly.
		generator();
                double quant = (*f)(X)/((importance)?((*p)(X)):(1./V));
		fsum_temp += quant;
		fsq_temp += quant*quant;
		M++;
	}
	// The updated sum and sum of squares are now divided by the updated M.
	fsum = fsum_temp/M;
	fsq = fsq_temp/M;
}	

double Mcint::error(){
	// Returns 2.58 times the standard error on the mean to give the bounds 
	// around the mean (the integral's value) of the 99% confidence interval
	// of the integrals value. With this value, it can be stated with 99% 
	// confidence that the true integral value lies in that range.
	return 2.58*sqrt((fsq-pow(fsum, 2))/M);
}

void Mcint::integrate(){
	// Performs one iteration (sample double) of the MC integration, ceasing
	// the calculation if the maximum number of samples are exceeded, or if
	// the desired accuracy is reached.
	while(R>eps){
		sample();
		R = error();
		// Output the number of doublings to the terminal, i.e. the number
		// of iteration calls, and the value of the integral at this point.
		int doublings = (int)log2(M/M_init);
		cout << doublings << " steps: " << fsum << endl;
		errors.push_back(R);
		steps.push_back(M);
		vals.push_back(fsum);
		if(M>=M_max){
			// Condition for maximum sample limit.
			cout << "Threshold of steps reached. Integration terminated." << endl;
			break;
		}
	}
}

void Mcint::generator(){
	// Uses a random number generator to update the vector of random variables
	// to pass to the function each time a sample is needed.
	for(int i(0); i<Nd; i++){
		// For each variable, calculates the random deviate (uniform) and then
		// applies the transformation method if cummulative distribution function
		// is supplied.
		double rx = myran.doub();
		X[i] = x_low[i] + (importance)?((*cy)(rx)):((x_up[i]-x_low[i])*rx);	
	}
}
	



