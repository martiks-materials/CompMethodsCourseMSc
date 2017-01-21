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
	double xrange(2.0), dx(1E-2), eps(1E-6);
	Rungekutta4 Rk4 = {y, &erfint, dx, eps, xrange, 1};
	while(Rk4.x_now<xrange){
		Rk4.iterate();
	}
	cout << "x = " << Rk4.xlist.back() << endl;
	cout << "y = " << Rk4.yn[0].back() << endl;
	return 0;
}
		
