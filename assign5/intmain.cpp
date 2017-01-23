// Computational Methods Assignment 5 - ODE Integration
// Martik Aghajanian Cohort 8
// 
// Program for determining the solution of question 1 in assignment 5
// using the adaptive step-size Runge-kutta 45 method. This outputs 
// the values of both dependent variables and the independent variable
// to separate files which can be plotted in gnuplot. 
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <valarray>
#include "specifics.h"
#include "rk4.h"
using namespace std;

int main(){
	// Initiate a valarray (vector with element wise addition and scalar
        // multiplication) 'y', number of elements/dimensions 'n_d' set to 2.
	//  The double precision values
        // 'xmax' and 'xmin' give the limits of integration whilst 'dx' is
        // the inital step size and 'eps' is the desired tolerance used in
        // the adaptive step size. To attain a step size that is slightly
        // more likely to be more optimal, a 'safety' factor is added too. 
	
	valarray<double> y = {1.5, 1.5};
	int n_d(2);
	double dx(0.1), xmax(10), xmin(0), eps(1E-6), safety(0.98);
	
	// The 'collapse' argument for this instance of Rungekutta4 is set
	// to false since the vector function depends on vector of dependent
	// variables 'y'. This RK45 method is iterated until the desired 
	// upper limit is reached.

	Rungekutta4 Rk4(y, &f, dx, eps, xmax, xmin, n_d, false, safety);	
	while(Rk4.x_now<xmax){
		Rk4.iterate();
	}
	
	
	// Data for (x, y1), (x, y2), and (y1, y2) is stored for plotting.

	ofstream outfile1("runge1.dat");
	ofstream outfile2("runge2.dat");
	ofstream outfile3("runge3.dat");
	if(((!outfile1.is_open())||(!(outfile2.is_open()))||(!outfile3.is_open()))){
		cout << "Error opening file..." << endl;
		return 1;
	}
	

	cout << "No. Steps = " << Rk4.numvals-1 << endl;
        cout << "Fnc. Evaluations = " << Rk4.evals << endl;
        cout << "Repeated steps = " << Rk4.repeats << endl;

	for(int i(0); i<Rk4.numvals ; i++){
		outfile1 << Rk4.xlist.at(i) << " " << Rk4.yn[0].at(i) << endl;
		outfile2 << Rk4.xlist.at(i) << " " << Rk4.yn[1].at(i) << endl;
		outfile3 << Rk4.yn[0].at(i) << " " << Rk4.yn[1].at(i) << endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	return 0;
}
			
	
