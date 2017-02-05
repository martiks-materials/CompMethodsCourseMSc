#include <iostream>
#include <vector>
#include <cmath>
#include "rand.h"
#include "mcmc.h"
using namespace std;

MarkovChain::MarkovChain(func fn, int dim, vec init, int b, int seed, vec sig_init, double epsy, int maxs, 
		    int mins, bool fix, int sigper, int checkper, bool multi) {
	// Initialize member function and variables for optimisation, seeding generator 
	f = fn;
	Nd = dim;
	myran.seed(seed);
	xi = init;	
	sig = sig_init;

	// Variables for Gaussian Distribution via Box-Muller or for multivariate
	// normal distribution.
	stored = 0.;
	extra = false;
	multivar = multi;
		
	// Variables for counting steps and identifying regions of chain
	counter = 0;
	fvalcount = 0;
	maxstep = maxs;
	minstep = mins;
	func_evals = 0;
	accepts = 0;
	burntime = b;
	burnin = true;
	// Variables for determining the variance of the function and x points
	fsum = 0;
	fsq = 0;
	maxf = 0;
	xmax = xi;
	eps = (multivar)?(epsy*10):epsy;
	sig_period = sigper;
	check_period=checkper;
	fixprop = fix;

	// For multivariate Gaussian, initialises covariance matrix with initial 
	// standard deviation (squared) which lasts until 'sig_period' has been exceeded.
	if(multivar){
		for(int i(0); i< Nd; i++){
			vec cov_column;
			for(int j(0); j< Nd; j++){
				double value = (i==j)?(sig_init[i]*sig_init[i]):0;
				cov_column.push_back(value);
			}
			cov_mat.push_back(cov_column);
		}
	}
	// Makes room for Nd vectors, each tracks the value of a particular variable.
	for(int i(0); i< Nd; i++){
		xvals.push_back({});
		sigvals.push_back({});
		xsum.push_back(0);
		xsum_b.push_back(0);
		xsq.push_back(0);
		xsq_b.push_back(0);
	}
}



double MarkovChain::normran(){
	// Implements the Box-Muller transform to generate 2 normally ditributed numbers
	// from 2 uniform deviates. Each first call of the function stores the second 
	// generated number in a member variable 'stored', indicated by the bool 'extra'
	// whilst every second call will return the 'stored' number.
	if(!extra){
		// Jacobian of this transformation is a product of two independent
		// normal distribution.
		double x1 = myran.doub();
		double x2 = myran.doub();
		double y1 = sqrt(-2*log(x1))*cos(2*M_PI*x2);
		double y2 = sqrt(-2*log(x1))*sin(2*M_PI*x2);
		stored = y2;
		extra = true;
		return y1;
	}
	else {
		extra = false;
		return stored;
	}
}
	
vec MarkovChain::multivarnorm() {
	// Cholesky Decomposition to then obtain a multivariate Gaussian distribution 
	// from a uniform deviate generator. L is a lower triangular matrix for which
	// L multiplied by it's transpose will produce the covariance matrix. It 
	// essentially acts as the square root of the covariance matrix.

	vecvec L;

	// Multivariate normal distribution and Box-Muller generated x vectors
	vec x_mv, x_bm;
	for(int i(0); i< Nd; i++){
		x_bm.push_back(normran());
	}
	// Sets empty L matrix and sets up covariance matrix
	for(int i(0); i< Nd; i++) {
		// For a vector of vectors, the [i][j] element is the jth element of the ith
		// vector, so each vector is a column.
		vec L_column;
		for(int j(0); j<Nd; j++) {
			L_column.push_back(0);
		}
		L.push_back(L_column);
	}
	// Now implements the algorithm for Cholesky decomposition of the Covariance matrix.
	for(int i(0); i< Nd; i++) {
		for(int j(0); j < i+1; j++) {
			double sum = 0;
			for(int k(0); k < Nd; k++) {
				sum += L[k][i]*L[k][j];
			}
			L[j][i] = (i==j)?sqrt(cov_mat[i][i]-sum):((cov_mat[j][i]-sum)/L[j][j]);
		}
	}
	// Matrix multiplication for new multivariate normal distributed vector. This is done in
	// such a way that the mean of this distribution is the vector of current x values and the
	// covariance matrix of the current 'sig_period' represents the spread of this distribution.
	// This is the multidimensional equivalent of rescaling the normal distributed random variable
	// to have a specific mean and standard deviation.
	for(int i(0); i< Nd; i++){
		double x_comp = xi[i];
		for(int j(0); j <Nd; j++){
			x_comp += L[j][i]*x_bm[i];
		}
		x_mv.push_back(x_comp);
	}
	return x_mv;
}

bool MarkovChain::iterate() {
	// Performs one step of the Metropolis algorithm to find the next point in
	// the Markov Chain, which will return a boolean to indicate the success
	// of acceptance. The newly generated 'x_temp' is tested.
	vec x_temp;
	bool accept;
	if(!multivar){
		// For the case of a proposal function which is just a product of
		// independent gaussian, the Box-Muller Transform is used.
		for(int i(0); i < Nd; i++) {
			// The Box-Muller transform produces a Gaussian for which the mean is
			// zero and the standard deviation is 1, so it must be re-scaled here.
			x_temp.push_back(xi[i] + sig[i]*normran());
		}
	}
	else {
		// For the case of a multivariate proposal function, the trial x vector is produced 
		// entirely by a separate function here.
		x_temp = multivarnorm();	
	}

	// Metropolis Algorithm Step: Calculated the relative value of the function which works 
	// as a probability if less than unity, causing automatic acceptance
	double quant = (*f)(x_temp);
	double alpha = quant/(*f)(xi);
	func_evals += 2;
	if(alpha>1.){
		xi = x_temp;
		accept = true;
	}
	else {
		// Produces uniform deviate to determine acceptance if 'alpha' < 1
		double beta = myran.doub();
		xi = (beta<alpha)?x_temp:xi;
		accept = (beta<alpha)?true:false;
	}

	// Now the function assesses convergence and updates member variables 
	if(accept){
		counter++;
		accepts++;
		// When evaluating the variance for x values, the position is
		// used to remove old terms in the 'xsum(_b)' and 'xsq(_b)' values
		int position = (burnin)?(counter):(counter-burntime);	

		if(burnin){
			// I.e. if the system state is still in the burn-in period.
			for(int i(0); i<Nd; i++){
				// Re-evaluates the sum of each x components over the chain and the sum of 
				// the squares of each x component over the chain so that their respective
				// variances are only taken over a finite 'sig_period'.
				xsum_b[i] += xi[i];
				xsq_b[i] += xi[i]*xi[i];
				xsum_b[i] -= (position > sig_period)? xvals[i][position-1-sig_period]:0; 
				xsq_b[i] -= (position > sig_period)? pow(xvals[i][position-1-sig_period], 2):0; 
			}
			for(int i(0); i<Nd; i++) {
				// Update the values of x recorded by the class and update the proposal function
				// if necessary. The proposal function is not changed until the number of steps taken
				// is sufficient to calculate a suitable variance in each x component to change it.
				xvals[i].push_back(xi[i]);
				if(!multivar){
					// If the normal distribution is a simple product of  Gaussian's, then the 
					// individual standard deviations of the Gaussiansin each direction is updated.
					sig[i]  = (counter>=sig_period)?sqrt(x_variance(i)):sig[i];
				}
				else {
					// If the proposal function is a multivariate Gaussian with coupled variables, then
					// the covariance matrix used is updated instead of the standard deviation of x vectors
					for(int i(0); i< Nd; i++) {	
						for(int j(0); j<Nd; j++) {
							if(i==j){
								cov_mat[i][i] = (counter>=sig_period)?x_variance(i):cov_mat[i][i];;
							}
							else {
								cov_mat[j][i] = (counter>=sig_period)?covariance(i, j):cov_mat[j][i];
							}
					
						}
					}
				}
			}
			fsum += quant;
			fsq += quant*quant;
			// Here, if the counter reaches a whole period of 'check_period', then the autocorrelation
			// is checked for all values to see if any of them are below the threshold. If ALL components
			// of the auto-correlation function are below the 0.05 threshold then the burn-in period ceases.
			if(counter%check_period == 0){ 
				for(int k(1); k< counter; k++) {
					vec rk;
					for(int i(0); i< Nd; i++) {
						double rki = abs(autocorr(k, i));
						rk.push_back(rki);
					}
					bool corr = false;
					for(int i(0); i<Nd ; i++){
						if(rk[i]>0.05){
							corr=true;
						}
					}
					if(!corr){
						// If the autocorrelation function produces a values for ALL components
						// that are below the 0.05 threshold, then the burn-in period is ended.

						cout << "Autocorrelation func. r(" << k << ") = ";
						for(int i(0); i< Nd; i++){	
							cout << rk[i] << "  ";  
						}
						cout << endl;
						burnin=false;
						burntime = counter;
						varold = variance();
	
						// Here it resets the recorded values of function and points visited
						// with the burn-in points stored in an additional vector. This is so 
						// the burn-in values and the non-burn-in values can be outputted to
						// different files and to easily distinguish between them.
						xburn = xvals;
						for(int i(0); i<Nd; i++){
							xvals[i].clear();
						}
						fsum = 0;
						fsq = 0;
						break;
					}	
				}
			}						
			if(burnin&&(counter == burntime)){
				// This condition describes the case if the auto-correlation function does not reduce
				// sufficiently and the maximum burn-in period is reached.

				burnin=false;
				varold = variance();

				// Here it resets the recorded values of function and points visited
				// with the burn-in points stored in an additional vector. This is so 
				// the burn-in values and the non-burn-in values can be outputted to
				// different files and to easily distinguish between them.
				xburn = xvals;
				for(int i(0); i<Nd; i++){
					xvals[i].clear();
				}
				fsum = 0;
				fsq = 0;
				cout <<  "Final sig after burnin  = " << sig[0] << ", " << sig[1] << endl;
			}
		}
		else {
			// If the system state is not in the burn-in period, then the variances and values are 
			// still updated and the convergence of the standard deviation of the function evaluations
			// is assessed to determine stopping the MCMC optimisation.
			for(int i(0); i<Nd; i++){
				xsum[i] += xi[i];
				xsq[i] += xi[i]*xi[i];
				xsum[i] -= (position > sig_period)? xvals[i][position-1-sig_period]:0; 
				xsq[i] -= (position > sig_period)? pow(xvals[i][position-1-sig_period], 2):0; 
			}
			for(int i(0); i < Nd; i++){
				// Update the xvals travelled to, and the sigmas with them.
				xvals[i].push_back(xi[i]);
				sigvals[i].push_back(sig[i]);
				if(!fixprop) {
					// If the proposal function is not fixed after the burn-in period then either the
					// standard deviation or covariance matrix is updated depending on whether or not
					// the proposal function is multivariate or not.	
					if(!multivar){
						// If the proposal function is not fixed after burn-in and is a product
						// of independent Gaussians.
						sig[i] = (counter-burntime>=sig_period)?sqrt(x_variance(i)):sig[i];
					}
					else {
						// If the proposal function is not fixed after burn-in and is a multivariate 
						// Gaussian with coupled variables, then the covariance matrix is updated 
						for(int i(0); i< Nd; i++) {	
							for(int j(0); j<Nd; j++) {
								if(i==j){
									cov_mat[i][i] = (counter-burntime>=sig_period)?x_variance(i):cov_mat[i][i];
								}
								else {
									cov_mat[j][i] = (counter-burntime>=sig_period)?covariance(i, j):cov_mat[j][i];
								}
							}
						}
					}
				}	
			}
			fvals.push_back(quant);
			fvalcount++;
			// Remove the value from fsum so we keep the window of convergence assessment moving
			double removed = (counter > (2*burntime))?fvals[counter-1-2*burntime]:0;
			fsum -= (counter > (2*burntime))?removed:0;
			fsq -= (counter > (2*burntime))?(removed*removed):0;
			fsum += quant;
			fsq += quant*quant;
			if(quant > maxf){
				// If the current function evaluation is higher than the current maximum.
				maxf = quant;
				for(int i(0); i< Nd; i++){
					xmax[i] = xi[i];

				}
			}
			/*
			int div = int(0.01*maxstep);
		
			if(counter%div==0){		
				cout << "Step " << counter-burntime << " after burn-in. X_variance = ( ";
				for(int i(0); i< Nd; i++){
					cout << sqrt(x_variance(i)) << " ";
				}
				cout << "). F_Variance: " << variance(); 
				cout << " with f(x) =  "<< (*f)(xi) << endl;
			}
			*/
		}		
	}
	return accept;
}


int MarkovChain::optimise() {
	// Performs the iteration of MCMC until convergence of the standard deviation of the function evaluations
	// is obtained. If this does not occur then the MCMC method is stopped after a 'maxstep' point.
	while(counter<maxstep){
		bool accept = iterate();
		if(accept){
			if(counter>=(2*burntime)){
				// The width of the window is the same as the burn-in-period for the moving average over 
				// which the variance is taken, is equal to the burn-in period. This means that the 
				// variance can only be assessed once  a full period of values is recorded and the 
				// variance can be compared properly.

				double v = variance();
				// If the standard deviation over the interval converges then algorithm is stopped.
				if(abs(sqrt(varold)-sqrt(v))/sqrt(abs(v))<eps){
					cout << "Converged with variance: " << v << " in " << counter << " steps. " <<  endl;
					converged = true;
					return 0;
				}
				// Stores the old variance for comparison in the next step.
				varold = v;
			}		
		}
	}
	// This is outputted if no convergence has been reached and the maximum step number is reached.
	converged = false;
	cout << "Maximum steps reached: " << counter <<  endl;
	return 0;	
}



double MarkovChain::x_variance(int i) {
	// Returns the variance of the x vector component 'i', and selects the data set depending on whether or not
	// the system is in a burn-in period or not. This is always over 'sig_period' number of values.
	int numvals = sig_period;
	double mean = ((burnin)?xsum_b[i]:xsum[i])/double(numvals);
	double meansq = ((burnin)?xsq_b[i]:xsq[i])/double(numvals);
	return meansq - (mean*mean);
}

double MarkovChain::variance() {
	// Returns the variance of the function evaluations over the period 'burntime'.
	double M = double(burntime);
	return (fsq/M)- pow((fsum/M), 2);
}


double MarkovChain::autocorr(int k, int comp){
	// Returns the 'k'th auto-correlation function to determine when to stop the burn-in period.
	// This is for the 'comp'th component of the x vector.
	double mean, denom, numer;
	for(int i(0); i< counter; i++){
		mean += xvals[comp][i]/double(counter);
	}
	for(int i(0); i< counter-k ; i++){
		denom += pow((xvals[comp][i]-mean), 2);
		numer += (xvals[comp][i]-mean)*(xvals[comp][i+k]-mean);
	}
	return numer/denom;
}

double MarkovChain::covariance(int i, int j) {
	// Produces the covariance matrix of a windwo of particular x values taken over
	// 'sig_period' steps.
	double M = double(sig_period);
	// Number of vals;
	int position = (burnin)?(counter):(counter-burntime);	
	double mean_i = ((burnin)?xsum_b[i]:xsum[i])/M;
	double mean_j = ((burnin)?xsum_b[j]:xsum[j])/M;
	double covar(0);
	for(int k = (position-sig_period);k<position; k++){
		covar += ((xvals[i][k]-mean_i)*(xvals[j][k]-mean_j))/M;
	}
	return covar;
}


