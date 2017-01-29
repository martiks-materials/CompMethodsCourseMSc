// Computational Methods Assignment 6 - Monte Carlo Integration
// Martik Aghajanian, Cohort 8
//
// Program to integrate the error function evaluated at 2 using Monte Carlo
// integration, using both a uniform sampling and importance sampling for
// two separate probability distribution functions. Outputs the data for each 
// case into respective data sets. 
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "mcint.h"
#include "multifunc.h"
using namespace std;


int main() {
	// 'Ndim' gives the number of dimensions of the problem, i.e. the number of independent
	// variables the integrand depends on. The first iterations to estimate the first guess
	// of the integration is done using 'mstart' samples. This is done with some random
	// for which 'see' gives the seed. The integration will stop once a maximum of 'maxsteps' 
	// iterations (each step doubles the number of samples in the integration). The variables 
	// 'mstart' and 'maxsteps' can be combined to produce the variables 'mmax' which represents
	// how many samples this requires to integrate the function. The desired tolerance
	// that the integration aspires to achieve, and for which the integration will cease if the 
	// maximum steps is not reached, is 'eps'. 

	int Ndim(1), mstart(500), see(142), maxsteps(22);
	unsigned long int mmax = int(mstart*pow(2, maxsteps));
	double eps(1E-6);

	// Three separate intgration objects are instantiated between the same vector of lower and 
	// upper limits, represented by 'a' and 'b' respectively. The first has the boolean variable
	// 'importance' set to false whilst the other two have it set to true, so importance sampling
	// can be implemented.
        vec a = {0}, b = {2};
        Mcint Monte1 = {erfi, a, b, Ndim, see, eps, mstart, false, pdf, invcdf, mmax};
        Mcint Monte2 = {erfi, a, b, Ndim, see, eps, mstart, true, pdf, invcdf, mmax};
        Mcint Monte3 = {erfi, a, b, Ndim, see, eps, mstart, true, pdf2, invcdf2, mmax};
        
	// All three are integrated until either the maximum number of sample doublings has been 
	// reached or that the tolerance has been achieved. The data of convergence of both value
	// and error (to 0) is recorded into .dat files to be plotted.
	Monte1.integrate();
	Monte2.integrate();
	Monte3.integrate();
	ofstream outfile1a("ertrack1.dat");
	ofstream outfile1b("numtrack1.dat");
	ofstream outfile2a("ertrack2.dat");
	ofstream outfile2b("numtrack2.dat");
	ofstream outfile3a("ertrack3.dat");
	ofstream outfile3b("numtrack3.dat");
	for(int i(1); i <= maxsteps; i++){
		outfile1a << log2(Monte1.steps[i]/mstart) << "   " << log(Monte1.errors[i]) << endl;
		outfile1b << log2(Monte1.steps[i]/mstart) << "   " << Monte1.vals[i] << endl;
		outfile2a << log2(Monte2.steps[i]/mstart) << "   " << log(Monte2.errors[i]) << endl;
		outfile2b << log2(Monte2.steps[i]/mstart) << "   " << Monte2.vals[i] << endl;
		outfile3a << log2(Monte3.steps[i]/mstart) << "   " << log(Monte3.errors[i]) << endl;
		outfile3b << log2(Monte3.steps[i]/mstart) << "   " << Monte3.vals[i] << endl;
	}
	outfile1a.close();
	outfile1b.close();
	outfile2a.close();
	outfile2b.close();
	outfile3a.close();
	outfile3b.close();
        cout << "Uniform Sampling: " << Monte1.fsum_new << " +/- " << Monte1.error2() <<  endl;
	cout << "Importance Sampling (PDF1): " << Monte2.fsum_new << " +/- " << Monte2.error2() <<  endl;
	cout << "Importance Sampling (PDF2): " << Monte3.fsum_new << " +/- " << Monte3.error2() <<  endl;
        return 0;
}
