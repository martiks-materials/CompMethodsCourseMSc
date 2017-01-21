// Computational Methods - Assignment 5
// This contains the class definition and the member function prototypes
// for the Rungekutta45 class which carries out generally coupled ODE solving,
// used in the assignment.
#pragma once
using namespace std;

//typedef valarray<double> vec;
//typedef valarray<double> (*mypointertype) (valarray<double>, double); 

class Rungekutta4{
// This class takes an initial N element solution vector, and a vector function
// that defines the ODE problem, and iterates in dependent variable steps to a 
// solution. This was used to gradually update the solution and provide simplicity 
// when attemping to store and obtain the independent variable and vector of 
//dependent variables at each step.
private:
	// The variable 'y' is the updatable current value of the solution vector
	// at each time step, and the function 'fnc' is that which provides the
	// ODE problem (d(y)/dx = fnc). Also a desired tolerance 'eps' is used
	// for selecting the desired relative accuracy of the solution at each
	// step of the dependent variable. This class is also supplied a 'limit'
	// which prevents the independent variable from exceeding a particular
	// value when implementing an adaptive step size.

	valarray<double> y;
        // mypointertype fnc;
	valarray<double> (*fnc)(valarray<double>, double );
	double eps;
	double limit;
public:	
	// For N coupled ODE's, 'yn' is the vector of vectors, for which each
	// nested vector stores the evolution of one of the dependent variables
	// in the initial vector for all steps. The evolution of the independent
	// variable (x) is stored in 'xlist'.

	vector< vector<double> > yn;
	vector<double> xlist;

	// The current value of the independent variable is 'x_now', whilst the 
	// current step size is 'h'. To aid the adapive step size calculation,
	// a safety factor 'S' can be used. The integers 'counter' and 'N' 
	// represent the number of steps taken and number of dependent variables
	// in the solution vector respectively.
	double x_now, h, S;
	int counter, N;
	Rungekutta4(valarray<double> yin, valarray<double> (*func)(valarray<double>, double), double dx, double epsilon, double xmax, int Np=2){
		y = yin;
		N = Np;
		fnc = func;
		h = dx;
		for(int i(0); i<N ; i++){
			vector<double> temp = {yin[i]};
			yn.push_back(temp);
		}
		/*
		for(int j(0); j<N; j++){	
			yn[j].push_back(yin[j]);
		} 
		*/
		limit = xmax;
		eps = epsilon;
		x_now = 0;
		counter = 1;
		S = 0.95; 
		xlist.push_back(x_now);
	}

	void step(valarray<double> &ys , double xs, double d){
		// Implements the 4th order Runge-Kutta Method for a vector function
		// of both independent variable 'x' and vector of dependent variables
		// 'y', over a step size of 'd'. Updates the current member variable
		// in the 'y' vector, opposed to returning a result.
		valarray<double> k1 = (*fnc)(ys,xs)*d;
                valarray<double> k2 = (*fnc)(ys+0.5*k1, xs+0.5*d)*d;
                valarray<double> k3 = (*fnc)(ys+0.5*k2, xs+0.5*d)*d;
		valarray<double> k4 = (*fnc)(ys+k3, xs+d)*d;
		for(int i(0); i<N ; i++){
                        ys[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
                }
	}
	
	void iterate(){
		// Performs an iteration of the Runge-Kutta45 method, in which step
		// doubling is used.
		valarray<double> y_half = y;
		if(x_now+h>limit){
			h = limit-x_now;
		}
		step(y, x_now, h);
		step(y_half, x_now, 0.5*h);
		step(y_half, x_now+0.5*h, 0.5*h);
		valarray<double> delta1 = y_half-y;
		x_now += h;
		xlist.push_back(x_now);
		h *= S*pow(abs(eps*y.max())/abs(delta1.max()), 0.2);
		y = y_half + (1./15)*delta1;
		for(int i(0); i<N ; i++){
                        yn[i].push_back(y[i]);
		}
		counter++;
	}
};

