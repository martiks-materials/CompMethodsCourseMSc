// Contains the functions necessary to answer question 1 in assignment 9.
#include <valarray>
#include <cmath>
#include "rk4.h"

std::valarray<double> f(std::valarray<double> y, double x){
	// Vector function used for assignment 5 question 1 which
	// depends only on the dependent variable vector 'y'.
	// While this function demands an independent variable due
	// to the requirement of the function in the Runge-Kutta
	// class, it is not needed in this particular function so
	// to avoid 'unused variable' messages, it is voided.
	(void)x;
	std::valarray<double> yn = {0, 0};
	yn[0] = -y[0]*y[0] - y[1];
	yn[1] = 5*y[0] -y[1];
	return yn;
} 

double solution(double guess){
	double x_start = guess;
	valdoub y_start = {1.5, x_start};
	//y0[0] = 1.5;
	//y0[1] = guess;	
	int n_d(2);
	double dx(0.1), xmax(10), xmin(0), eps(1E-6), safety(0.98);
	// The 'collapse' argument for this instance of Rungekutta4 is set
	// to false since the vector function depends on vector of dependent
	// variables 'y'. This RK45 method is iterated until the desired 
	// upper limit is reached.
	Rungekutta4 Rk4(y_start, &f, dx, eps, xmax, xmin, n_d, false, safety);	
	while(Rk4.x_now<xmax){
		Rk4.iterate();
	}
	//cout << Rk4.y[0]<< endl;
	return Rk4.y[0];
}
