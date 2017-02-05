#include <iostream>
#include <cmath>
#include <vector>
#include "rand.h"
using namespace std;

typedef vector<double> vec;
typedef vector< vector<double> > vecvec;
typedef double (*func)(vector<double>);

class MarkovChain {
private:
	// Function f, xvector xi, stdev vector sig
	func f;
	vec xi, sig;
	double  eps, stored, fsum, fsq;
	bool extra, burnin, fixprop;
	int seed;
	Ran myran;
public:
	// max of function maxf, storage varaible varold, 
	// vectors fvals stores function values, xsig stored the sum of 
	double maxf, varold;
	bool converged;
	vec fvals, fburn, xsum, xsq, xsum_b, xsq_b, xmax;
	vecvec xvals, xburn, sigvals;
	int Nd, burntime, counter, maxstep, minstep, fvalcount, sig_period, check_period;
	MarkovChain(func fn, int dim, vec init, int b, int seed, vec sig_init, double epsy, int maxs, int mins, bool fix, int sigper, int checkper) {
		f = fn;
		Nd = dim;
		myran.seed(seed);
		burntime = b;
		burnin = true;
		fixprop = fix;
		xi = init;
		sig = sig_init;
		sig_period = sigper;
		check_period=checkper;
		// Variables for Gaussian Distribution via Box-Muller
		stored = 0.;
		extra = false;
		
		// Variables for counting steps and identifying regions of chain
		counter = 0;
		fvalcount = 0;
		maxstep = maxs;
		minstep = mins;

		// Variables for determining the variance of the function and x points
		fsum = 0;
		fsq = 0;
		maxf = 0;
		xmax = xi;
		eps = epsy;

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

	double normran(){
		// Implements the Box-Muller transform to generate 2 normally ditributed numbers
		// from 2 uniform deviates. Each first call of the function stores the second 
		// generated number in a member variable, whilst every second call will return
		// the stored number.
		if(!extra){
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
	

	bool iterate() {
		// ACCEPTING STEP: WE DONT COUNT THE ONES THAT ARE REJECTED
		vec x_temp;
		bool accept;
		//_________________________METROPOLIS ALGORITHM____________________
		for(int i(0); i < Nd; i++) {
			x_temp.push_back(xi[i] + sig[i]*normran());
		}
		double alpha = (*f)(x_temp)/(*f)(xi);
		if(alpha>1.){
			xi = x_temp;
			accept = true;
		}
		else {
			double beta = myran.doub();
			xi = (beta<alpha)?x_temp:xi;
			accept = (beta<alpha)?true:false;
		}
		//______________________________________________________________
		// CONVERGENCE AND RECORDING STEPS
		if(accept){
			counter++;
			int position = (burnin)?(counter):(counter-burntime);	
			double quant = (*f)(xi);
			if(burnin){	
				for(int i(0); i<Nd; i++){
					// Re-evaluates the sum of each x components over the chain and the sum of 
					// the squares of each x component over the chain so that their respective
					// variances are only taken over a finite 'sig_period'
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
					sig[i]  = (counter>=sig_period)?sqrt(x_variance(i)):sig[i];
				}
				fsum += quant;
				fsq += quant*quant;
				// Here, if the counter reaches a whole period of 'burn-in check', then the autocorrelation
				// is checked for all values to see if any of them are below the threshold.
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
							for(int i(0); i< Nd; i++){	
								cout << "r(" << k << ") = " << rk[i] << " for component " << i+1<< endl; 
							}
							burnin=false;
							burntime = counter;
							//varold = fmean();  //
							varold = variance();
							cout << "Burnin Variance: " << varold << endl;
		
							// Here it resets the recorded values of function and points visited
							// with the burn-in points stored in an additional vector.
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
					burnin=false;
					//varold = fmean();  //
					varold = variance();
					cout << "Burnin Variance: " << varold << endl;

					// Here it resets the recorded values of function and points visited
					// with the burn-in points stored in an additional vector.
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
				for(int i(0); i<Nd; i++){
					xsum[i] += xi[i];
					xsq[i] += xi[i]*xi[i];
					xsum[i] -= (position > sig_period)? xvals[i][position-1-sig_period]:0; 
					xsq[i] -= (position > sig_period)? pow(xvals[i][position-1-sig_period], 2):0; 
				}
				for(int i(0); i < Nd; i++){
					// Update the xvals we travel to and the sigmas with them too
					// xvals is reset when we exit burnin 
					xvals[i].push_back(xi[i]);
					sigvals[i].push_back(sig[i]);
					if(!fixprop) {	
						sig[i] = (counter-burntime>=sig_period)?sqrt(x_variance(i)):sig[i];
					}	
				}
				fvals.push_back(quant);
				fvalcount++;
				// Remove the value from fsum so we keep the window moving
				double removed = (counter > (2*burntime))?fvals[counter-1-2*burntime]:0;
				fsum -= (counter > (2*burntime))?removed:0;
				fsq -= (counter > (2*burntime))?(removed*removed):0;
				fsum += quant;
				fsq += quant*quant;
				if(quant > maxf){
					maxf = quant;
					for(int i(0); i< Nd; i++){
						xmax[i] = xi[i];
					}
				}
				int div = int(0.01*maxstep);
				/*
				if(counter%div==0){		
					cout << "Step " << counter-burntime << " after burn-in. X_variance = ( ";
					for(int i(0); i< Nd; i++){
						cout << sqrt(x_variance(i)) << " ";
					}
					cout << "). F_Variance: " << variance(); 
					cout << " Mean: " << fmean() << " with f(x) =  "<< (*f)(xi) << endl;
				}
				*/
			}		
		}
		return accept;
	}


	int optimise() {
		while(counter<maxstep){
			bool accept = iterate();
			if(accept){
				if(counter==burntime){
				}
				else if(counter>=(2*burntime)){
					// The width of the window for the moving average over which the variance is taken,
					// is equal to the burn-in period. This means that the variance can only be assessed once
					// a full period of values is recorded and the variance can be compared properly.
					//double v = fmean(); // 

					double v = variance();
					//if(sqrt(v)<10){
					if(abs(sqrt(varold)-sqrt(v))/sqrt(abs(v))<eps){
						cout << "Converged with variance: " << v << " in " << counter << " steps. " <<  endl;
						converged = true;
						return 0;
					}
					varold = v;
				}		
			}
		}
		converged = false;
		cout << "Maximum steps reached: " << counter <<  endl;
		return 0;	
	}
	
	double x_variance(int i) {
		int numvals = sig_period;
		double mean = ((burnin)?xsum_b[i]:xsum[i])/double(numvals);
		double meansq = ((burnin)?xsq_b[i]:xsq[i])/double(numvals);
		return meansq - (mean*mean);
	}

	double variance() {
		double M = double(burntime);
		return (fsq/M)- pow((fsum/M), 2);
	}
	
	double fmean(){
		double M = double(burntime);
		return fsum/M;
	}

	double autocorr(int k, int comp){
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
};














