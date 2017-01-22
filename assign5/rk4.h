// Computational Methods - Assignment 5
// This contains the class definition and the member function prototypes
// for the Rungekutta45 class which carries out generally coupled ODE solving,
// used in the assignment.
#pragma once
using namespace std;

typedef vector<double> vecdoub;
typedef valarray<double> valdoub;
typedef valarray<double> (*function) (valarray<double>, double); 

class Rungekutta4{
// This class takes an initial N element solution vector, and a vector function
// that defines the ODE problem, and iterates in dependent variable steps to a 
// solution. This was used to gradually update the solution and provide simplicity 
// when attemping to store and obtain the independent variable and vector of 
// dependent variables at each step.
private:
	// The variable 'y' is the updatable current value of the solution vector
	// at each time step, and the function 'fnc' is that which provides the
	// ODE problem (d(y)/dx = fnc). Also a desired tolerance 'eps' is used
	// for selecting the desired relative accuracy of the solution at each
	// step of the dependent variable. This class is also supplied a 'limit'
	// which prevents the independent variable from exceeding a particular
	// value when implementing an adaptive step size.

	valdoub y, k_val;
        function fnc;
	double eps;
	double limit;
public:	
	// For N coupled ODE's, 'yn' is the vector of vectors, for which each
	// nested vector stores the evolution of one of the dependent variables
	// in the initial vector for all steps. The evolution of the independent
	// variable (x) is stored in 'xlist'.

	vector< vector<double> > yn;
	vecdoub xlist;

	// The current value of the independent variable is 'x_now', whilst the 
	// current step size is 'h'. To aid the adapive step size calculation,
	// a safety factor 'S' can be used. The integers 'counter' and 'N' 
	// represent the number of steps taken and number of dependent variables
	// in the solution vector respectively.
	double x_now, h, safety;
	int  evals, numvals, repeats, N;
	bool collapse;
	Rungekutta4(valdoub yin, function func, double dx, double epsilon, double xmax, double xmin=0, int ns=2, bool col=false, double safe=1.){
		y = yin;
		N = ns;
		fnc = func;
		h = dx;
		for(int i(0); i<N ; i++){
			vecdoub temp = {yin[i]};
			yn.push_back(temp);
		}
		limit = xmax;
		eps = epsilon;
		x_now = xmin;
		safety = safe;
		collapse = col;
		evals = 0;
		numvals = 1;
		repeats = 0;
		xlist.push_back(x_now);
	}

	void step(valarray<double> &ys , double xs, double d, bool reuse){
		// Implements the 4th order Runge-Kutta Method for a vector function
		// of both independent variable 'x' and vector of dependent variables
		// 'y', over a step size of 'd'. Updates the current member variable
		// in the 'y' vector, opposed to returning a result. Since the full
		// is always carried out first in this class, then the bool 'half' 
		// allows for one less function evaluation. If 'half' is false then
		// class member variable 'k1' is updated, whilst 'half' being true
		// means the current 'k1' is used.

		valdoub k1 = ((reuse)?k_val:(*fnc)(ys,xs))*d;
                valdoub k2 = (*fnc)(ys+0.5*k1, xs+0.5*d)*d;
                valdoub k3 = (collapse)? k2 : ((*fnc)(ys+0.5*k2, xs+0.5*d)*d);
		valdoub k4 = (*fnc)(ys+k3, xs+d)*d;
		evals += ((collapse)? 3 : 4) - ((reuse)? 1 : 0);
		for(int i(0); i<N ; i++){
                        ys[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
                }
	}
	
	void iterate(){
		// Performs an iteration of the Runge-Kutta45 method, in which step
		// doubling is used.
		h = (x_now+h>limit)? limit-x_now: h;
		valdoub y_c = y, y_half = y;
		k_val = (*fnc)(y, x_now);
		step(y_c, x_now, h, true);
		step(y_half, x_now, 0.5*h, true);
		step(y_half, x_now+0.5*h, 0.5*h, false);
		valdoub delta1 = y_half-y_c;
		// The current step size is then altered to that which gives the desired
		// accuracy, additionally using a'safety factor' S.
		double h_new = h*safety*pow(abs(eps*y_c.max()/delta1.max()), 0.2);
		if(h>h_new){
			// If the calculated optimal step size is less than the step size
			// then the step must be redone with the step-doubling to make the
			// calculation 5th order accurate.
			h = h_new;
			valdoub y_new = y,
	                y_half = y;
	                step(y_new, x_now, h, true);
			step(y_half, x_now, 0.5*h, true);
	                step(y_half, x_now+0.5*h, 0.5*h, false);
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
	
		xlist.push_back(x_now);
		for(int i(0); i<N ; i++){
                        yn[i].push_back(y[i]);
		}
		numvals++;
	}
};

