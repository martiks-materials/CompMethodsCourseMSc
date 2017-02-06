// Computational Methods - Assignment 7: Markov Chain Monte Carlo
// Martik Aghajanian
//
// Question 1:
// Program for mapping and optimising the Rosenbrock function using 
// Markov Chain Monte Carlo. This corresponds to answering question
// 1 of the assignment and outputs files for the mapped function 
// distribution and also the burn-in points. 
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "mcmc.h"
#include "rand.h"

double rosenbrock(vec x) {
	// Returns the Rosenbrock function if the output is positive and
	// returns 0 if negative, which aids the Metropolis Algorithm.
	double result =  1000 - (1-x[0])*(1-x[0]) -100*pow(((x[0]*x[0]) - x[1]), 2);
	return (result>0)?result:0;
}

int main() {
	// This is the master 'distributor' which is used to seed the generators of all 
	// Markov chains and also decide the initial points of each chain
	Ran distributor;
	distributor.seed(424242); 

	// Vector 'funcmaxes' for storing the maximum points of each of the Markov chains
	vec funcmaxes, startsig = {5, 5};
	vecvec xmaxes;
	int numchains(12),  sig_period(1E2), check_period(1E2), burn(1E4), maxi(1E6), Ndim(2);

	// Counters for the total steps, total function evaluations and number of converged MC's
	int totstep(0), totvals(0), numcon(0);
	
	// Convergence tolerance and initial point range.
	double epsig(1E-8), range(1.0);
	for(int i(0); i< Ndim; i++){
		xmaxes.push_back({});
	}
	ofstream outfile1("markovdata.dat");
	ofstream outfile2("burndata.dat");	
	for(int i(0); i<numchains; i++){
		// Produces 'numchains' distinct Markov chains at various starting points with various seeds.
		int see = distributor.int64();
		vec xinit;
		for(int j(0); j< Ndim; j++){
			double starter = -range + (2*range*distributor.doub());
			xinit.push_back(starter);
		}
		cout << "__________________________________________________________" << endl;
		cout << "Markov Chain " << i+1<< " starting at ("<< xinit[0]<< ", " << xinit[1]<< ") " << endl;
		
		// For each randomly initialised state (starting point and seed), the Markov Chain is set off
		// to optimise with the proposal function fixed after the burn-in period (fixprop=true)
		// Also the multivariate bool at the end is set to false to have a proposal function which is
		// a product of independent Gaussians.

		MarkovChain Marko = {rosenbrock, Ndim, xinit, burn, see, startsig, epsig, maxi, 10, true, sig_period, check_period, false};
		Marko.optimise();
		for(int i(0); i < Marko.fvalcount; i++) {
			for(int j(0); j < Marko.Nd; j++) {
				outfile1 << Marko.xvals[j][i] << " ";
			}
			outfile1 << Marko.fvals[i] << " " << endl;
		}
		for(int i(0); i< Marko.burntime; i++){
			for(int j(0); j< Marko.Nd; j++){
				outfile2 << Marko.xburn[j][i] << " ";
			}
			outfile2 << endl;
		}
		cout << "Max(f(x)) =  " << Marko.maxf << " at (" << Marko.xmax[0] << ", " << Marko.xmax[1] << " )" <<  " with variance " << Marko.variance() << endl;
		cout << Marko.accepts << " steps with  " << Marko.func_evals << " function evaluations" << endl;
		numcon += (Marko.converged)?1:0;
		totvals += Marko.func_evals;
		totstep += Marko.counter;
		funcmaxes.push_back(Marko.maxf);
		for(int i(0) ; i<Ndim; i++){	
			xmaxes[i].push_back(Marko.xmax[i]);
		}
	}
	// Of all the chains, the maximum one is outputted to the terminal along with the value which maximises it.
	auto result = max_element(begin(funcmaxes), end(funcmaxes));
	cout << "Largest Function maximum: " << double(*result) << " at ";
	for(int i(0);i<numchains; i++){
		if(abs(funcmaxes[i]-double(*result))<1E-10){
			cout << "(" <<  xmaxes[0][i] << ", " << xmaxes[1][i] << ") ";
		}
	}
	cout << endl << numcon << "/" << numchains << " converged." << endl;
	cout << totstep << " total steps needed with " << totvals << " function evaluations." <<  endl;
	outfile1.close();
	outfile2.close();
	return 0;
}
