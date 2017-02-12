#pragma once
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

typedef double (*func)(double);
typedef vector<double> vec;

class BrentMethod {
private:
	func f;
	bool converged, bisect;
public:
	double a, b, c, fa, fb, fc, tol, pq1, pq2, delta;
	vec fvals;
	vec bestvals;
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
	
	double operate() {
		while(!converged){
			if((fa*fb>0)&&(fb*fc>0)){
				// If it is found that the root is not bracketed either side by the
				// current set of x values, then the algorithm can not function and
				// here it is terminated. This condition can be found using the sign
				// of the product of fb with the function evaluated at the neighbouring
				// points a and c.
				cout << "fa = " << fa << ", fb = " << fb << ", fc = " << fc << endl;
				cout << "Failed to bracket root." << endl;
				break;
			}
			
			if(abs(fc)< abs(fb)){
				// Here if the function at the contrapoint c is closer to the root than
				// it is at b, b and c are switched here, so that b is now the best guess
				// and if the secant method is triggered, the result will have less of the
				// function (y) axis to traverse in the finite difference approximation of
				// the gradient in the secant method.	
				double xtemp = b, ftemp = fb;
				fb = fc;
				b = c;
				c = xtemp;
				fc = ftemp;
			}

			double x, pq0;
			if(abs(fa-fc)>1E-16){
				// If it is not the case that the double precision numbers fa and fc
				// are not within a certain precision of each other (equal), then use
				// inverse quadratic interpolation. This fits an inverse parabola and 
				// sets the dependent variable to zero. Since this interpolates, the
				// root found between the contrapoint and the best guess should be 
				// close to the true root.
				double R = fb/fc, S = fb/fa, T = fa/fc;
				double P = S*(T*(R-T)*(c-b) - (1-R)*(b-a));
				double Q = (T-1)*(R-1)*(S-1);
				pq0 = P/Q;
				x = b + P/Q;
			}
			else {
				// If the values of a and c are the same then the inverse quadratic
				// interpolation fails and the secant method (finite difference Newton
				// Raphson Method) must be used to estimate the root linearly.
				cout << "secant" << endl;
				x = b -fb*(b-c)/(fb-fc);
				pq0 = x-b;
			}
			// Here a bool for allowing the interpolation is determined so that Brent's
			// conditions are satisfied. If the previous step was a bisection, then the 
			// addition to the best guess to get the new best guess (P/Q) must be smaller
			// than the addition to the best guess in the previous step, 'pq1', provided that said
			// previous addition was not ridiculuously small (smaller than delta), which is to
			// avoid slow convergence but to avoid overshooting the root.
			//
			// If the previous step was an interpolation, then the same condition holds but 
			// instead of the previous addition, the method looks at the 2nd to last addition 
			// which is 'pq2'.
			bool allowed=true;
			if(bisect) {
				allowed = ((abs(x-b)<0.5*abs(pq1))&&(abs(pq1)>delta));
			}
			else {
				allowed = ((abs(x-b)<0.5*abs(pq2))&&(abs(pq2)>delta));
			}	

			if ((abs(x-b) <= 0.75*abs(b-c))&&allowed){
				cout << "brent" << endl;
				// Provided that the addition to the current best guess does not exceed 3/4 the
				// way between the best guess and the contrapoint, and that the previous conditions
				// allow it, the inverse quadratic interpolation step is accepted and the values 
				// are updated. The previous best guess becomes the 'a' value for the next step. and
				// a marker 'bisect' is set up for the next steps conditions on allowing the acceptance
				// of that steps inverse quadatic interpolated guess.
				a = b;
				fa = fb;
				b = x;
				fb = (*f)(b);
				bisect = false;
			}
			else {	
				// If the conditions for accepting the inverse quadratic interpolated step are not met,
				// then the next best guess 'b' is obtained by a bisection step, by halfing the bracketed
				// area. A marker 'bisect' is set so that the next iteration is aware of this and can
				// decide whether or not to accept that next steps interpolated guess.
				a = b;
				fa = fb;
				b = 0.5*(b+c);
				fb = (*f)(b);
				bisect = true;
			}

			// These are the stored values of the additions to the best guess over the last two iterations,
			// used to decide the fate of the next interpolation. Here they are updated and the step is 
			// counted towards the total number of iterations.
			pq1 = pq0;
			pq2 = pq1;
			count++;

			if((fb*fc>0)&&(fa*fb<0)){
				// Here if the contra point c and the 'best guess' b do not bracket
				// the root but the a point does, then switch a and c, so that the next
				// steps of the algorithm can work.
				double xtemp = a, ftemp = fa;
				fa = fc;
				a = c;
				c = xtemp;
				fc = ftemp;
			}
			// The check for convergence here is assessed and found to be apparent if it is the case that 
			// either the function is within machine accuracy of zero, or that the region of bracketing between
			// the current best guess and 'contrapoint' is below a certain tolerance.
			if((0.5*abs(b-c)<tol)||abs(fb)<1E-16){
				cout << "Converged in " << count << " steps" << endl;
				converged = true;
			}
			else if(count > maxsteps) {
				cout << "Maximum iterations reached" << endl;
				break;
			}
		}
		// The bracketing value which has the lowest absolute function value is outputted.
		return (abs(fb)<abs(fc))?b:c;
	}


};
