// Computational Methods Assignment 5 - ODE Integration
// Martik Aghajanian, Cohort 8
//
// This file contains the class member functions of the "RungeKutta4"
// class used for assignment 5.
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
	numvals = 1;
	repeats = 0;

	// If function has no dependence on the dependent variable vector then the
	// method collapses to 3-point Simpson's rule. In this case then 'k_rep' is
	// used to reduce function evaluations. If regular RK45 is needed then 
	// 'k_rep' is not used and is initially set to any valarray.
	
	k_rep = (collapse)? ((*fnc)(yin, x_now)) : yin;
	evals = (collapse)? 1:0;
	xlist.push_back(x_now);
	for(int i(0); i<N ; i++){
		vecdoub temp = {yin[i]};
		yn.push_back(temp);
	}
}

void Rungekutta4::step(valarray<double> &ys , double xs, double d, bool reuse, bool single, bool endhalf){
        // Implements the 4th order Runge-Kutta Method for a vector function
        // of both independent variable 'x' and vector of dependent variables
        // 'y', over a step size of 'd'. Updates the current member variable
        // in the 'y' vector, opposed to returning a result. Since the function
        // evaluation in the 'k1' can be re-used, then the 'reuse' bool
        // allows the old value to be placed instead of re-evaluating 'fnc'.
	
        valdoub k1 = ((reuse)?k_val:(*fnc)(ys,xs))*d;
        valdoub k2 = (*fnc)(ys+0.5*k1, xs+0.5*d)*d;
        valdoub k3 = (collapse)? k2 : ((*fnc)(ys+0.5*k2, xs+0.5*d)*d);
	
        if(single&&collapse) {
		// If collapsing to a 3-point rule, then the k_rep is updated so
		// it can be used as the starting k_val for the next step. It is
		// only desirable to update this if it is the single step so that
		// this function evaluation can be used at the end of the second
		// half step. Hence the 'single' bool is used.
		k_rep = (*fnc)(ys+k3, xs+d);
		evals++;
	}
	
	// For the collapsed 3-point rule, the updated 'k_rep' can be used for the 
	// 'k4' of the first step ('single'=true') and also for the latter half step
	// ('endhalf'=true). For the first half step this value cannot be used. 

	valdoub k4 = (((collapse)&&(single||endhalf))? k_rep: ((*fnc)(ys+k3, xs+d)))*d;
	for(int i(0); i<N ; i++){
        	ys[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
        }
	
	// The number of function evaluations in each call of iterate(), is reduced 
	// from 4 to 3 if k2 = k3 in the collapsed procedure, and is reduced by 1
	// further if the first value k_val is reused. Finally, reusing the 'k_rep'
	// value reduces by a further 1.
       
	evals += ((collapse)? 3 : 4) + ((reuse)? -1 : 0) + (((collapse)&&(single||endhalf))?-1:0);
        
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

	// If this procedure is to collapse to 3-point rule, then the new k_val
	// is the old k4 from the last successful step. This reduces the function
	// evaluations used.

        k_val = (collapse)? k_rep: ((*fnc)(y, x_now));
	evals += (collapse)? 0 : 1;
	
	Rungekutta4::step(y_c, x_now, h, true, true, false);
        Rungekutta4::step(y_half, x_now, 0.5*h, true, false, false);
        Rungekutta4::step(y_half, x_now+0.5*h, 0.5*h, false, false, true);
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
                Rungekutta4::step(y_new, x_now, h, true, true, false);
                Rungekutta4::step(y_half, x_now, 0.5*h, true, false, false);
                Rungekutta4::step(y_half, x_now+0.5*h, 0.5*h, false, false, true);
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

