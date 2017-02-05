// Computational Methods Assignment 7 - Markov Chain Monte Carlo
// Martik Aghajanian, Cohort 8
// 
// Class declaration for the Markov Chain Monte Carlo function exploration
// and maximisation class which takes in an initialised random state and specified
// initial variables, mapping out and optimising a specific user-defined function.
// Options are available for the proposal function to be fixed after burn-in period
// and also to have a multivariate Gaussian as the proposal function too.

#include <iostream>
#include <cmath>
#include <vector>
#include "rand.h"
using namespace std;

// Useful typedefs for functions, vectors and vectors of vectors, since
// such are initialised several times in this program. 
typedef vector<double> vec;
typedef vector< vector<double> > vecvec;
typedef double (*func)(vector<double>);

class MarkovChain {
private:
	// Variable Definitions:
	// 'f': Function of interest, 'xi': Current x vector, 'sig': Current x vector std. dev.
	// 'eps': Tolerance for function standard deviation convergence
	// 'stored': 2nd normal deviate obtained via Box-Muller Transform
	// 'fsum': Updatable sum of function evaluations
	// 'fsq': Updatable sum of squared evaluations
	// 'extra': Indicator of stored variable for Box-Muller Transform
	// 'burnin': Indicator of system state being in burn-in period
	// 'fixprop': Indicator of fixing proposal function after burn-in period
	// 'seed': of the uniform deviate generator 'myran'

	func f;
	vec xi, sig;
	double  eps, stored, fsum, fsq;
	bool extra, burnin, fixprop;
	int seed;
	Ran myran;
public:
	// 'maxf': Current maximum function evaluation, 'varold': storage variable for variance
	// 'converged': Indicator of algorithm termination due to convergence of variance
	// 'multivar': Indicator of multivariate Gaussian proposal function
	// 'fvals': Container for function evaluations after burn-in period exceeded
	// 'fburn': Container for function evaluations during burn-in period
	// 'xsum': Container for sum of x vector values over a specified 'sig_period'
	// 'xsq': Container for sum of squared x vector components over specified 'sig_period'
	// 'xsum_b': Container for sum of x vector values over a specified 'sig_period' during burn-in
	// 'xsq_b': Container for sum of squared x vector components over specified 'sig_period' during burn-in
	// 'xmax': Current x vector for which function is maximal, corresponding to 'maxf'
	// 'xvals': Container for all post-burn-in period x vector values
	// 'xburn': Container for x vector values during burn-in
	// 'sigvals': Container for standard deviations of the proposal function for xi
	// 'cov_mat': Updatable covariance matrix
	// 'Nd': Dimension of 'xi', 'burntime': Length of burn-in period, 'counter': of steps
	// 'maxstep': Maximum number of successful Metropolis steps before non-convergence termination
	// 'minstep': Minimum number of successful Metropolis steps before convergence-based termination
	// 'fvalcount': Size of 'fvals' vector, 'sig_period': Steps over which x_variance() is assessed
	// 'check_period': Steps before auto-correlation function evaluated.
	// 'func_evals': Number of function evaluations, 'accepts': Number of accepted Metropolis steps

	double maxf, varold;
	bool converged, multivar;
	vec fvals, fburn, xsum, xsq, xsum_b, xsq_b, xmax;
	vecvec xvals, xburn, sigvals, cov_mat;
	int Nd, burntime, counter, maxstep, minstep, fvalcount, sig_period, check_period, func_evals, accepts;
	
	MarkovChain(func fn, int dim, vec init, int b, int seed, vec sig_init, double epsy, int maxs, 
		    int mins, bool fix, int sigper, int checkper, bool multi);

	double normran();

	vec multivarnorm(); 

	bool iterate(); 

	int optimise(); 
	
	double x_variance(int i); 

	double variance();
	
	double autocorr(int k, int comp);
	
	double covariance(int i, int j);

};
