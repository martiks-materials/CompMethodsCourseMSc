#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "mcmc.h"

typedef vector<double> vec;

double rosenbrock(vec x) {
	return 1000 - (1-x[0])*(1-x[0]) -100*pow(((x[0]*x[0]) - x[1]), 2);
}

int main() {
	//MarkovChain(func fn, int dim, vec init, int b, int seed, vec sig_init, double thresh, int maxs) {
	vec xinit = {0, 0}, startsig = {50, 50};
	int burn(1E5), see(4242), maxi(1E6);
	double epsig(1E-7);	
	MarkovChain Marko = {rosenbrock, 2, xinit, burn, see, startsig, epsig, maxi, 10};
	Marko.optimise();
	ofstream outfile("markovdata.dat");	
	
	for(int i(0); i<Marko.fvalcount; i++) {
		for(int j(0); j < Marko.Nd; j++){
			outfile << Marko.xvals[j][i] << " ";
		}
		outfile << Marko.fvals[i] << " " << endl;
	}	
	cout << "Max(f(x)) =  " << Marko.maxf << " at (" << Marko.xmax[0] << ", " << Marko.xmax[1] << " )" <<  " with variance " << Marko.variance() << endl;
	outfile.close();
	return 0;
}
