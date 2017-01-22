// Computational Methods Assignment 5 - ODE Integration
// Martik Aghajanian, Cohort 8
//
// This file contains the class member functions of the "RungeKutta4"
// class used for assignment 5."
#include <vector>
#include <valarray>
#include <cmath>
#include "rk4.h"
using namespace std;


Rungekutta4::Rungekutta4(valdoub yin, function func, double dx, double epsilon, 
		double xmax, double xmin, int ns, bool col, double safe){
	y = yin;
	fnc = func;
	N = ns;
	h = dx;
	limit = xmax;
        eps = epsilon;
	x_now = xmin;
	safety = safe;
	collapse = col;
	evals = 0;														               numvals = 1;
	repeats = 0;
	xlist.push_back(x_now);
	for(int i(0); i<N ; i++){
		vecdoub temp = {yin[i]};
		yn.push_back(temp);
	}
}

void Rungekutta4::step(valarray<double> &ys , double xs, double d, bool reuse){
        // Implements the 4th order Runge-Kutta Method for a vector function
        // of both independent variable 'x' and vector of dependent variables
        // 'y', over a step size of 'd'. Updates the current member variable
        // in the 'y' vector, opposed to returning a result. Since the function
        // evaluation in the 'k1' can be re-used, then the 'reuse' bool
        // allows the old value to be placed instead of re-evaluating 'fnc'.
        valdoub k1 = ((reuse)?k_val:(*fnc)(ys,xs))*d;
        valdoub k2 = (*fnc)(ys+0.5*k1, xs+0.5*d)*d;
        valdoub k3 = (collapse)? k2 : ((*fnc)(ys+0.5*k2, xs+0.5*d)*d);
        valdoub k4 = (*fnc)(ys+k3, xs+d)*d;
        evals += ((collapse)? 3 : 4) - ((reuse)? 1 : 0);
        for(int i(0); i<N ; i++){
        	ys[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
        }       
}

void Rungekutta4::iterate(){
        // Performs an iteration of the Runge-Kutta45 method, in which step
        // doubling is used.
        // If the step h takes the independent variable over the maximum limit
        // then the step size is changed to take it exactly to the limit.

        h = (x_now+h>limit)? limit-x_now: h;

        // The step-doubling is performed using copies of the updatable valarray.
        // The k_val is updated prior to this so that it can be reused in the
        // first two calls to the iterate() function for efficiency. After both
        // full-step and half-steps are performed, the error 'delta1' is estimated.

        valdoub y_c = y, y_half = y;
        k_val = (*fnc)(y, x_now);
        Rungekutta4::step(y_c, x_now, h, true);
        Rungekutta4::step(y_half, x_now, 0.5*h, true);
        Rungekutta4::step(y_half, x_now+0.5*h, 0.5*h, false);
        valdoub delta1 = y_half-y_c;

        // The current step size is then altered to that which gives the desired
        // accuracy, additionally using a'safety factor.

        double h_new = h*safety*pow(abs(eps*y_c.max()/delta1.max()), 0.2);

	if(h>h_new){
                // If the calculated optimal step size is less than the step size
                // then the step must be redone with the step-doubling to determine
                // new 'delta1' error and make the calculation 5th order accurate.

                h = h_new;
                valdoub y_new = y,
                y_half = y;
                Rungekutta4::step(y_new, x_now, h, true);
                Rungekutta4::step(y_half, x_now, 0.5*h, true);
                Rungekutta4::step(y_half, x_now+0.5*h, 0.5*h, false);
                delta1 = y_half-y_new;
                y = y_half + (1./15)*delta1;
                x_now += h;
		repeats++;
	}
	else {
                // If the step taken originally was smaller than the step size needed,
                // then the step is not re-done, but the optimal step size is used as
                // the initial step size used for the next step.

                x_now += h;
                y = y_half + (1./15)*delta1;
		h = h_new;

        }

        // Once the step has been determined to be either optimal or to have smaller
        // error than necessary, the member data structures of this classes are
	// appropriately updated.

        xlist.push_back(x_now);
	for(int i(0); i<N ; i++){
                yn[i].push_back(y[i]);
        }
	numvals++;
}

