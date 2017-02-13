#include <iostream>
#include <vector>
#include "rootmethods.h"

double Bisector::operate() {
	// Iteration of the bisection method, where the region of interest in half
	// whilst still bracketing the root, to converge to a given tolerance
	bool chose1;
	while(count < maxstep){
		if(f1*f2>=0) {
			// If the values of the function at either sides of the
			// bracketing region have opposite sign, then the bisection
			// method can not work as it can only move to somewhere 
			// inside the initially given region.
			std::cout << "Failed to bracket." << std::endl;
			break;
		}
		// Calculated the midpoint as the next best guess
		double x3 = 0.5*(x1+x2);
		f3 = (*f)(x3);
		// 'dx' will be the separation between the bracketing points which
		// will be used to determine convergence
		double dx;
		if(f3*f2<0) {
			// If 'x3' and 'x2' bracket the root then change the old
			// 'x1' to be the new calculated 'x3'.
			x1 = x3;
			f1 = f3;
			dx = 0.5*std::abs(x1-x2);
			chose1 = false;
		}
		else if(f1*f3<0){
			// If 'x3' and 'x1' bracket the root then change the old 
			// 'x1' to be the new calculated 'x3'.
			x2 = x3;
			f2 = f3;
			dx = 0.5*std::abs(x1-x2);
			chose1=true;
		}
		count++;
		// If it is the case that the function evaluated at the current 
		// best guess is less than a small number corresponding to machine
		// tolerance (i.e. it is zero), or the separation between bracketed points
		// is below the specified tolerance of the class, then the algorithm 
		// terminates.
		if((std::abs(f3)<1E-16)||dx<tol){
			return x3;
			break;
		}
	}
	double result = chose1? x2:x1;
	return result;
}


double BrentMethod::operate() {
	while(!converged){
		if((fa*fb>0)&&(fb*fc>0)){
			// If it is found that the root is not bracketed either side by the
			// current set of x values, then the algorithm can not function and
			// here it is terminated. This condition can be found using the sign
			// of the product of fb with the function evaluated at the neighbouring
			// points a and c.
			std::cout << "fa = " << fa << ", fb = " << fb << ", fc = " << fc << std::endl;
			std::cout << "Failed to bracket root." << std::endl;
			break;
		}
		
		if(std::abs(fc)< std::abs(fb)){
			// Here if the function at the contrapoint c is closer to the root than
			// it is at b, then b and c are switched here, so that c now becomes the 
			// new best guess b, meaning that if the secant method is triggered, the 
			// result will have less of the function (y) axis to traverse in the finite 
			// difference approximation of the gradient in the secant method.	
			double xtemp = b, ftemp = fb;
			fb = fc;
			b = c;
			c = xtemp;
			fc = ftemp;
		}
		// 'x' will be defined as the calculated value for the potential new best guess, whilst
		// 'pq0' will be the current addition to the best guess to get the new one, whose size 
		// will determine if the interpolating step is allowed.
		double x, pq0;
		if(std::abs(fa-fc)>1E-16){
			// If it is not the case that the double precision numbers fa and fc
			// are not within a certain precision of each other (i.e. equal), then use
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
			x = b -fb*(b-c)/(fb-fc);
			pq0 = x-b;
		}
		// Here a bool for allowing the interpolation is determined, so that Brent's
		// conditions are satisfied. If the previous step was a bisection, then the 
		// addition to the best guess to get the new best guess (P/Q), must be smaller
		// than the previous addition to the best guess ('pq1'). as well as the conditions 
		// that 'pq1' was not ridiculuously small (smaller than delta), which is to
		// avoid slow convergence but to avoid overshooting the root at the same time.
		//
		// If the previous step was an interpolation, then the same condition applies but 
		// instead of the previous addition, the method looks at the 2nd to last addition 
		// which is 'pq2'.
		bool allowed=true;
		if(bisect) {
			allowed = ((std::abs(x-b)<0.5*std::abs(pq1))&&(std::abs(pq1)>delta));
		}
		else {
			allowed = ((std::abs(x-b)<0.5*std::abs(pq2))&&(std::abs(pq2)>delta));
		}	

		if ((std::abs(x-b) <= 0.75*std::abs(b-c))&&allowed){
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
		// The check for convergence here is assessed and found to be apparent if it
		// is the case that  either the function is within machine accuracy of zero,
		// or that the region of bracketing between the current best guess and 'contrapoint'
		// is below a certain tolerance.
		if((0.5*std::abs(b-c)<tol)||std::abs(fb)<1E-16){
			converged = true;
		}
		else if(count > maxsteps) {
			std::cout << "Maximum iterations reached" << std::endl;
			break;
		}
	}
	// The bracketing value which has the lowest absolute function value is outputted.
	return (std::abs(fb)<std::abs(fc))?b:c;
}

