#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <valarray>
#include "rk4.h"
using namespace std;

valarray<double> erfint(valarray<double> y, double x){
	(void)y;
	valarray<double> result = {(2./sqrt(M_PI))*exp(-1*x*x)};
	return result;
}

int main(){
	valarray<double> y = {0};
	double xrange(2.0), xmin(0.0), dx(1E-2), eps(1E-6);
	Rungekutta4 rk_simps = {y, &erfint, dx, eps, xrange, xmin, 1, true, 0.98};
	Rungekutta4 rk_full = {y, &erfint, dx, eps, xrange, xmin, 1, false, 0.98};
	while(rk_simps.x_now<xrange){
		rk_simps.iterate();
	}
	while(rk_full.x_now<xrange){
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
		
