// Computational Methods Assignment 9 - Shooting Method
// Martik Aghajanian, Cohort 8
//
// Class declarations for Brent's Method class which takes
// two initial root-bracketing points and iterates until given 
// convergence criteria is met.
#pragma once
#include <iostream>
#include <cmath>
#include <vector>

typedef double (*func)(double x);
typedef std::vector<double> vec;


class BrentMethod {
private:
	// Class which implements Brent's Method, by taking a function 'f'
	// to find the root of. The boolean 'converged' indicates to the class
	// that the method has reached the target tolerance and 'bisect' is a 
	// marker that the last step was a bisection to determine whether or not
	// an interpolation step is accepted.
	func f;
	bool converged, bisect;
public:
	// The class initially takes points 'a' and 'b' (best guess) which are
	// expected to bracket the root, whilst to begin with 'c' is set the same as
	// 'a' here. The variables 'fa', 'fb' and 'fc' are the corresponding function
	// evaluations. The convergence is assessed using 'tol' whilst the small 'delta'
	// is used to determine whether or not the interpolation step is accepted, which is
	// also the purpose of the updatable variables 'pq1' and 'pq2' which are the 
	// addtions to the best guess from the last and second to last step. The number
	// of steps is counted with 'count', terminating the algorithm if it reaches
	// a maximum 'maxsteps'.
	double a, b, c, fa, fb, fc, tol, pq1, pq2, delta;
	int count, maxsteps;
	BrentMethod(func fun, vec init, double eps, int maxstep, double del) {
		f = fun;
		a = init[1];
		b = init[0];
		c = init[1];
		fa = (*f)(a);
		fb = (*f)(b);
		fc = (*f)(c);
		if((fb*fc>0)&&(fa*fb<0)){
			// Here if the contra point c and the 'best guess' b do not bracket
			// the root but the a point does, then switch a and c, so that the next
			// steps of the algorithm can work. This is here for the first step
			// of the algorithm.
			double xtemp = a, ftemp = fa;
			fa = fc;
			a = c;
			c = xtemp;
			fc = ftemp;
		}
		tol = eps;
		delta = del;
		maxsteps=maxstep;
		converged = false;
		count = 0;
		bisect = false;
	}
	
	double operate();
};
