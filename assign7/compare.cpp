// Computational Methods - Assignment 7: Markov Chain Monte Carlo
// Martik Aghajanian
//
// Question 2: Comparing fixed proposal functions
// This computes the number of function evaluations of the Markov
// Chain Monte Carlo method for mapping and optimising a function, for
// both a fixed proposal function and a proposal function which varies
// after the burn-in period. This can be applied to a proposal function
// which is a product of independent Gaussians or a multivariate normal
// distribution.
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
	double result = 1000 - (1-x[0])*(1-x[0]) -100*pow(((x[0]*x[0]) - x[1]), 2);
	return (result>0)?result:0;
}

int main() {
	// This is the master 'distributor' which is used to seed the generators of all 
	// Markov chains and also decide the initial points of each chain
	Ran distributor;
	distributor.seed(42); 
	vec funcmaxes1, funcmaxes2, startsig = {5, 5};
	vecvec xmaxes1, xmaxes2;
	int numchains(12), sig_period(1E2), check_period(1E2),  burn(1E4), maxi(1E6), Ndim(2);
	int  numcon1(0), numcon2(0), numval1(0), numval2(0), numstep1(0), numstep2(0);
	double epsig(1E-8), range(1.0);
	for(int i(0); i< Ndim; i++){
		xmaxes1.push_back({});
		xmaxes2.push_back({});
	}
	for(int i(0); i<numchains; i++){
		// Produces 'numchains' distinct Markov chains at various starting points with various seeds.
		int see = distributor.int64();
		vec xinit1;
		vec xinit2;
		for(int j(0); j< Ndim; j++){
			double starter1 = -range + (2*range*distributor.doub());
			double starter2 = -range + (2*range*distributor.doub());
			xinit1.push_back(starter1);
			xinit2.push_back(starter2);
		}
		cout << "__________________________________________________________________" << endl;	
		cout << "Markov Chain " << i+1<< " (Fixed Prop. Function) starting at ("<< xinit1[0]<< ", " << xinit1[1]<< ")" << endl;
		MarkovChain Marko1 = {rosenbrock, Ndim, xinit1, burn, see, startsig, epsig, maxi, 10, true, sig_period, check_period, false};
		Marko1.optimise();
		cout << "Max(f(x)) =  " << Marko1.maxf << " at (" << Marko1.xmax[0] << ", " << Marko1.xmax[1] << " )" <<  " with variance " << Marko1.variance() << endl;
		cout << Marko1.accepts << " steps with " << Marko1.func_evals << " function evaluations" << endl;
		numval1 += Marko1.func_evals;	
		numstep1 += Marko1.counter;
	
		
		cout << "_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ __" << endl;	
		cout << "Markov Chain " << i+1<< " (Variable Prop. Function) starting at ("<< xinit2[0]<< ", " << xinit2[1]<< ")" << endl;
		MarkovChain Marko2 = {rosenbrock, Ndim, xinit2, burn, see, startsig, epsig, maxi, 10, false, sig_period, check_period, false};
		Marko2.optimise();
		cout << "Max(f(x)) =  " << Marko2.maxf << " at (" << Marko2.xmax[0] << ", " << Marko2.xmax[1] << " )" <<  " with variance " << Marko2.variance() << endl;
		cout << Marko2.accepts << " steps with " << Marko2.func_evals << "  function evaluations" << endl;
		numval2 += Marko2.func_evals;
		numstep2 += Marko2.counter;	
	
		numcon1 += (Marko1.converged)?1:0;
		numcon2 += (Marko2.converged)?1:0;
		funcmaxes1.push_back(Marko1.maxf);
		for(int i(0) ; i<Ndim; i++){	
			xmaxes1[i].push_back(Marko1.xmax[i]);
		}
		funcmaxes2.push_back(Marko2.maxf);
		for(int i(0) ; i<Ndim; i++){	
			xmaxes2[i].push_back(Marko2.xmax[i]);
		}
	}
	// Of all the chains, the maximum one is outputted to the terminal along with the value which maximises it.
	cout << "__________________________________________________________________" << endl;	
	auto result1 = max_element(begin(funcmaxes1), end(funcmaxes1));
	cout << "Largest Function maximum (Fixed Proposal Function): " << double(*result1) << " at ";
	for(int i(0);i<numchains; i++){
		if(abs(funcmaxes1[i]-double(*result1))<1E-10){
			cout << "(" <<  xmaxes1[0][i] << ", " << xmaxes1[1][i] << ") ";
		}
	}
	cout << endl << numval1 << " function evaluations." << endl;
	cout << numstep1 << " steps taken." << endl;
	cout << numcon1 << "/" << numchains << " converged." << endl;
	
	cout << "_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ __" << endl;	
	auto result2 = max_element(begin(funcmaxes2), end(funcmaxes2));
	cout << "Largest Function maximum (Variable Proposal Function): " << double(*result2) << " at ";
	for(int i(0);i<numchains; i++){
		if(abs(funcmaxes2[i]-double(*result2))<1E-10){
			cout << "(" <<  xmaxes2[0][i] << ", " << xmaxes2[1][i] << ") ";
		}
	}
	cout << endl << numval2 << " function evaluations." << endl;
	cout << numstep2 << " steps taken." << endl;
	cout << numcon2 << "/" << numchains << " converged." << endl;
	return 0;
}
