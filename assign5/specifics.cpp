#include <valarray>
using namespace std;

valarray<double> f(valarray<double> y, double x){
	// Vector function used for assignment 5 question 1 which
	// depends only on the dependent variable vector 'y'.
	// While this function demands an independent variable due
	// to the requirement of the function in the Runge-Kutta
	// class, it is not needed in this particular function so
	// to avoid 'unused variable' messages, it is voided.
	(void)x;
	valarray<double> yn = {0, 0};
	yn[0] = -y[0]*y[0] - y[1];
	yn[1] = 5*y[0] -y[1];
	return yn;
} 
