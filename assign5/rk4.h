// Computational Methods Assignment 5 - ODE Integration
// Martik Aghajanian Cohort 8
//
// Class declaration for Runge-Kutta '45' class which takes a user-specified function
// as well as initial values for the independent and dependent variables. This class
// can be used to iterate the adaptive step size method in order to obtain a 5th order
// accurate result from a 4th order method.  

#pragma once
using namespace std;

// Several valarrays, vectors and functions are initialised in this class so it is 
// sensible to define some typedefs 

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
	// at each time step and the function 'fnc' is that which provides the
	// ODE problem (d(y)/dx = fnc). To reduce the number of function evaluations
	// in the step-doubling procedure, 'k_val' is an updatable value of the 
	// function at the current step. Also a desired tolerance 'eps' is used
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
	// in the initial vector, for all steps. The evolution of the independent
	// variable (x) is stored in 'xlist'.

	vector< vector<double> > yn;
	vecdoub xlist;

	// The current value of the independent variable is 'x_now', whilst the 
	// current step size is 'h'. To aid the adapive step size calculation,
	// a 'safety' factor can be used. The integer 'N' is the number of dependent
	// variables in the vector. Respectively, the 'evals', 'numvals' and 'repeats' 
	// variables count the number of function evaluations, number of values in the 
	// 'xlist', and number of steps which had to be repeated. Finally, the bool
	// 'collapse' is used if the function does not depende on any of the dependent
	// variables and as a result the method collapsed to a 3-point Simpson's rule.

	double x_now, h, safety;
	int  evals, numvals, repeats, N;
	bool collapse;

	Rungekutta4(valdoub yin, function func, double dx, double epsilon, double xmax, 
				double xmin=0, int ns=2, bool col=false, double safe=1.);

	void step(valarray<double> &ys , double xs, double d, bool reuse);
		
	void iterate();
};
