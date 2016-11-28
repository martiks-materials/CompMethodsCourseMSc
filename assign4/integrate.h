// Computational Methods Assignment 4 - Numerical Integration
// Martik Aghajanian Cohort 8
//
// Class definition for integration class which takes a user-specified function
// as well as limits and a user-specified precision to which the method of 
// integration is expected to converge. This can be used to integrate functions
// defined on an interval, to within a specified accuracy.

#include <iostream>
#include <cmath>
using namespace std;

class Integrate{
// Class for integrating function "fnc" within limits "a", "b" within given precision "eps"
private:
	// For a function "fnc" which takes argument "x", with limits "a" and "b"
	// iteration the number of points "n" used in the calculation 
	// in addition to the boundary points. The convergence of the integral is
	// assessed using double precision float parameter "eps"
	double (*fnc)(double x);
	double a;
	double b;  	
	int n;
	double eps;

	// This class continually updates the variables "n" and "evals", the number
	// function evaluations made in total, upon calling the member function 
	// "iterate". The variables "I1" and "I2" are updatable values of the integral
	// values which are compared to assess convergence. The samples taken to 
	// evaluate the integral are spaced with spacing "h" which is halved upon
	// each iteration of doubling the number of points used.
	int evals;
	double h;
	double I2;
	double I1;
public:
	Integrate(double (*func)(double x), double low, double up, double epsilon){
		// Upon instantiation, three function evaluations are carried out to
		// provide the estimation for two integrals "I1" and "I2" which form
		// the first two stages in either the trapezoid or Simpson's rule.
		// The values of "n" and "h" are updated accordingly.
		a = low;
		b = up;
		h = 0.5*(b-a);
		n = 1;
		eps = epsilon;
		fnc = func;
		I1 = h*((*fnc)(a)+(*fnc)(b));
		I2 = ((0.5*I1)+((*fnc)(0.5*(a+b))*h));
		evals = 3;
	}

	void iterate() {
		// Member function which iterates the procedure of doubling the number of
		// points, and calculates the new value of the integral "I2" and compares
		// it to the previous value "I1" to assess convergence using relative error
		// "eps". By doubling the number of points "n", the spacing between points
		// "h" halves.
		I1 = I2;
		I2 = 0.5*I1;
		n *= 2;

		// Instead of recalculating the values at all points whose spacing is half
		// that of the previous iteration, the previous value is just halved, and 
		// the extra sampling points intermediate between those previously selected
		// are calculated and added to the new value of the integral estimation. This
		// is to minimise the number of function evaluations needed. The number of 
		// function evaluations "evals" is recorded and updated each time this iterate()
		// member function is called.
		for(int i(0); i<n; i++){
			double x_new = a + (i+0.5)*h;
			I2 += 0.5*h*(*fnc)(x_new);
			evals++;
		}
		h *= 0.5;
	}

	double trap(){
		// Using the iterate() member function, implements the trapezoid rule,
		// by increasing the number of points sampled in the integration until
		// the relative change in the integral evaluation upon iteration is
		// below the required precision.
		while(abs(I2-I1)>eps*abs(I1)){
			Integrate::iterate();
		}
		return I2;
	}
	
	double simp(){
		// Using the iterate() member function, implements the Simpson's rule 
		// by removing 1/3 times the previous evaluation of the trapezoid rule
		// from 4/3 times the first iteration of the trapezoid rule to get a value
		// for the integral evaluated using Simpson's rule "s2". This value is then
		// compared to the previous value using the "eps" parameter.
		double s1 = (4*I2/3.)-(I1/3.);
		Integrate::iterate();
		double s2 = (4*I2/3.)-(I1/3.);
		while(abs(s2-s1)>eps*abs(s1)){
			s1 = s2;
			Integrate::iterate();
			s2 = (4*I2/3.)-(I1/3.);
		}
		return s2;
	}
		
	// This member function outputs the number of function evalulations "evals"  
	// required to determine the integral of "fnc" over interval between "a" and
	// "b" using either the trapeziod rule or Simpson's rule.
	int get_evals() { return evals; }

};