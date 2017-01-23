// Computational Methods Assignment 5 - ODE Integration
// Martik Aghajanian Cohort 8
// 
// Program for determining the error function evaluated at x=2 using 
// both the Runge-Kutta 45 with adaptive step size, and the collapsed
// 3-point adaptive Simpson's rule. For both methods, the values, number
// of function evaluations, number of steps, and number of repeated
// steps are outputted. The error function integrand is included at the
// top of this program since it does not uptake significant space. 
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <valarray>
#include "rk4.h"
using namespace std;


valarray<double> erfint(valarray<double> y, double x){
	// Integrand for the error function. Takes a valarray of doubles, 
	// but does not use it, as the Runge-Kutta 45 class requires this 
	// form. To avoid 'unused variable' warnings, 'y' is made to be void.
	(void)y;
	valarray<double> result = {(2./sqrt(M_PI))*exp(-1*x*x)};
	return result;
}

int main(){
	// Initiate a valarray (vector with element wise addition and scalar
	// multiplication) 'y', having one element, hence setting the number
	// of dimensions 'n_d' to 1. The double precision values
	// 'xmax' and 'xmin' give the limits of integration whilst 'dx' is
	// the inital step size and 'eps' is the desired tolerance used in
	// the adaptive step size. To attain a step size that is slightly
	// more likely to be more optimal, a 'safety' factor is added too. 

	valarray<double> y = {0};
	int n_d(1);
	double xmax(2.0), xmin(0.0), dx(1E-2), eps(1E-6), safety(0.98);
	
	// The difference between these two instances of the Rungekutta4 objects,
	// is that 'rk_simps' is set as a 3-point adaptive Simpson's rule via the
	// boolean 'collapse' (set to true), whilst the original Runge-Kutta 
	// method without collapsing is implemented through 'rk_full' with 
	// 'collapse' set to false.

	Rungekutta4 rk_simps = {y, &erfint, dx, eps, xmax, xmin, n_d, true, safety};
	Rungekutta4 rk_full = {y, &erfint, dx, eps, xmax, xmin, n_d, false, safety};
	
	// Iterates the methods until the current independent variable reaches
	// upper integral limit.

	while(rk_simps.x_now<xmax){
		rk_simps.iterate();
	}
	while(rk_full.x_now<xmax){
		rk_full.iterate();
	}
	cout.precision(11);

	cout << "RK45 Method Collapsed to 3-point Simpson's Rule:" << endl;
	cout << "erf(x=" << rk_simps.xlist.back() << ") = " <<  rk_simps.yn[0].back() << endl;
	cout << "No. Steps = " << rk_simps.numvals-1 << endl;
	cout << "Fnc. Evaluations = " << rk_simps.evals << endl;
	cout << "Repeated steps = " << rk_simps.repeats << endl << endl;

	cout << "RK45 Original Method (No collapse):" << endl;
        cout << "erf(x=" << rk_full.xlist.back() << ") = " <<  rk_full.yn[0].back() << endl;
        cout << "No. Steps = " << rk_full.numvals-1 << endl;
	cout << "Fnc. Evaluations = " << rk_full.evals << endl;
	cout << "Repeated steps = " << rk_full.repeats << endl;
	return 0;
}
		
