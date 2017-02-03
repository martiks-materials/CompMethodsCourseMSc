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
	func f;
	vec xi, sig;
	double  eps, stored, fsum, fsq, fsum_b, fsq_b;
	bool extra, burnin;
	int seed;
	Ran myran;
public:
	double maxf, varold;
	vec fvals, xsig, xsig2, xmax;
	vecvec xvals, sigvals;
	int Nd, bperiod, counter, maxstep, minstep, fvalcount;
	MarkovChain(func fn, int dim, vec init, int b, int seed, vec sig_init, double epsy, int maxs, int mins) {
		f = fn;
		Nd = dim;
		myran.seed(seed);
		bperiod = b;
		burnin = true;
		xi = init;
		sig = sig_init;
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
			xsig.push_back(0);
			xsig2.push_back(0);
		}
	}

	double normran(){
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
		
		// CONVERGENCE AND RECORDING STEPS
		if(accept){
			counter++;
				
			for(int i(0); i<Nd; i++){
				xsig[i] += xi[i];
				xsig2[i] += xi[i]*xi[i];
			}
			double quant = (*f)(xi);
			if(burnin){	
				for(int i(0); i<Nd; i++) {
					double sigsig = (counter>5)?x_variance(i):sig[i];
					sig[i] = sqrt(sigsig);
				}
				fsum += quant;
				fsq += quant*quant;
				// Here there will be 'counter' terms in each fsum and fsq
				cout <<  "sig = " << sig[0] << ", " << sig[1] << endl;
			}
			else {
				for(int i(0); i < Nd; i++){
					// Update the xvals we travel to and the sigmas with them too
					xvals[i].push_back(xi[i]);
					sigvals[i].push_back(sig[i]);		
				}
				fvals.push_back(quant);
				fvalcount++;
				// Remove the value from fsum so we keep the window moving
				double removed =(counter > (2*bperiod))?fvals[counter-2*bperiod]:0;
				fsum -= (counter > (2*bperiod))?removed:0;
				fsq -= (counter > (2*bperiod))?(removed*removed):0;
				fsum += quant;
				fsq += quant*quant;
				cout << "Fsum = " << fsum << ", Fsq = "  << fsq << endl;
				if(quant > maxf){
					maxf = quant;
					for(int i(0); i< Nd; i++){
						xmax[i] = xi[i];
					}
				}		
				cout << "Variance: " << variance()<< ", sig = " << sig[0] << ", " << sig[1] << endl;
			}		
		}
		return accept;
	}


	int optimise() {
		while(counter<maxstep){
			bool accept = iterate();
			if(accept){
				if(counter==bperiod){
					burnin=false;
					varold = variance();
					cout << "Burnin Variance: " << varold << endl;
					fsum = 0;
					fsq = 0;
				}
				else if(counter>=(2*bperiod)){
					// FIX THIS, MAKE OLD VARIANCE VARIABLE TO UPDATE
					// YOU SHOULD HAVE SOME WINDOW OF ASSESSMENT
					double v = variance();
					if(abs(sqrt(varold)-sqrt(v))/sqrt(abs(v))<eps){
						cout << "Converged with variance: " << v << " in " << counter << " steps. " <<  endl;
						return 0;
					}
					varold = v;
				}		
			}
		}
		cout << "Maximum steps reached: " << counter <<  endl;
		return 0;	
	}
	
	double x_variance(int i) {
		return (xsig2[i]/double(counter))-pow((xsig[i]/double(counter)), 2);
	}

	double variance() {
		double M = double(bperiod);
		return (fsq/M)- pow((fsum/M), 2);
	}

			
};














