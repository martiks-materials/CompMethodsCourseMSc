// Computational Methods Assignment 8 - Root Finding
// Martik Aghajanian, Cohort 8
//
// Class declarations for the bisection method and Brent's Method
// classes which each take two initial root-bracketing points and
// iterate until given convergence criteria is met.
#pragma once
#include <iostream>
#include <cmath>
#include <vector>

typedef double (*func)(double x);
typedef std::vector<double> vec;

class Bisector{
private:
	// Bisection method class which takes the function for
	// which the root is to be found, 'f' and two initial
	// points 'x1' and 'x2' which bracket the root. The 
	// three updatable function values 'f1', f2', and 'f3'
	// are used as the evaluations of the two bracketing
	// points and the new bisected point respectively. The 
	// algorithm is iterated until convergence determined
	// by 'tol' or until a maximum number of steps 'maxstep'
	// is reached.
	func f;
	double x1, x2, f1, f2, f3, tol;
	int maxstep;
public:
	// The number of steps is counted using 'count'. 
	int count;
	Bisector(func fun, vec xinit, double eps, int maxsteps){
		f = fun;
		x1 = xinit[0];
		x2 = xinit[1];
		f1 = (*f)(x1);
		f2 = (*f)(x2);
		f3 = 0;
		tol = eps;
		count = 0;
		maxstep = maxsteps;
	}
	double operate();
};


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
