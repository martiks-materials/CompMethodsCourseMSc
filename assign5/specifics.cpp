#include <valarray>
using namespace std;

valarray<double> f(valarray<double> y, double x){
	(void)x;
	valarray<double> yn = {0, 0};
	yn[0] = -y[0]*y[0] - y[1];
	yn[1] = 5*y[0] -y[1];
	return yn;
} 
