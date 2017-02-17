// Computational Methods Assignment 9 - Shooting Method
// Martik Aghajanian Cohort 8
//
// Class declaration for Runge-Kutta '45' class which takes a user-specified function
// as well as initial values for the independent and dependent variables. This class
// can be used to iterate the adaptive step size method in order to obtain a 5th order
// accurate result from a 4th order method.  

#pragma once
// Several valarrays, vectors and functions are initialised in this class so it is 
// sensible to define some typedefs 

typedef std::vector<double> vecdoub;
typedef std::valarray<double> valdoub;
typedef std::valarray<double> (*function) (std::valarray<double>, double); 

class Rungekutta4{
// This class takes an initial N element solution vector, and a vector function
// that defines the ODE problem, and iterates the dependent variable in steps to a 
// solution. This was used to gradually update the solution and provide simplicity 
// when attemping to store and obtain the independent variable and vector of 
// dependent variables at each step.
private:
	// The variable 'y' is the updatable current value of the solution vector
	// at each time step and the function 'fnc' is that which provides the
	// ODE problem (d(y)/dx = fnc). To reduce the number of function evaluations
	// in the step-doubling procedure, 'k_val' is an updatable value of the 
	// function at the current step. In the case of the function not depending on
	// the vector of dependent variables, the updatable valarray 'k_rep' is used to 
	// further reduce the function evaluation count by reusing values from a previous 
	// step for a 3-point Simpson's rule. Also a desired tolerance 'eps' is used
	// for selecting the desired relative accuracy of the solution at each
	// step of the dependent variable. This class is also supplied a 'limit'
	// which prevents the independent variable from exceeding a particular
	// value when implementing an adaptive step size.

	valdoub  k_val, k_rep;
        function fnc;
	double eps;
	double limit;

public:	
	// For N coupled ODE's, 'yn' is the vector of vectors, for which each
	// nested vector stores the evolution of one of the dependent variables
	// in the initial vector, for all steps. The evolution of the independent
	// variable (x) is stored in 'xlist'.
	valdoub y;
	std::vector< std::vector<double> > yn;
	vecdoub xlist;

	// The current value of the independent variable is 'x_now', whilst the 
	// current step size is 'h'. To aid the adapive step size calculation,
	// a 'safety' factor can be used. The integer 'N' is the number of dependent
	// variables in the vector. Respectively, the 'evals', 'numvals' and 'repeats' 
	// variables count the number of function evaluations, number of values in the 
	// 'xlist', and number of steps which had to be repeated. Finally, the bool
	// 'collapse' is used if the function does not depend on any of the dependent
	// variables and as a result the method collapsed to a 3-point Simpson's rule.

	double x_now, h, safety;
	int  evals, numvals, repeats, N;
	bool collapse;
	
	// See rk4.cpp for descriptions of the member functions below.

	Rungekutta4(valdoub yin, function func, double dx, double epsilon, double xmax, 
				double xmin=0, int ns=2, bool col=false, double safe=1.);

	void step(std::valarray<double> &ys , double xs, double d, bool reuse, bool single, bool endhald);
		
	void iterate();
};
