// Computational Methods Assignment 8 - Root Finding
// Martik Aghajanian Cohort 8
//
// Program which finds the roots of the interpolation function
// of the data in assigment 2, via both the bisection method
// and Brents Method
#include <iostream>
#include <cmath>
#include <vector>
#include "splining.h"
#include "bisect.h"
#include <iomanip>
#include "brent.h"

int main(){
	// The range of the x-data from the spline is from 'xmin' to 'xmax'
	// 'epsilon' is the convergence tolerance convergence for either
	// method, whilst 'maxi' is the maximum amount of steps taken before the
	// algorithms are terminated. The variable 'delta' is used for the bisection 
	// method when searching for the second root.
	double epsilon(1E-8), delta(1E-8);
	int maxi(10000);
	// The bisection method is implemented here, with the starting points for the
	// 
	vec xinit1 = {2.0, 2.5};
	Bisector Splisect1 = {oldspline, xinit1, epsilon, maxi};
	double result1 = Splisect1.operate();
	vec xinit2 = {2.75, 3.0};
	Bisector Splisect2 = {oldspline, xinit2, epsilon, maxi};
	double result2 = Splisect2.operate();
	std::cout << "Bisection Root 1: " << std::setprecision(9)<< result1 << endl;
	std::cout << "Bisection Root 2: " << std::setprecision(9) << result2 << endl;
	// Brent's Method
	vec xinit3 = {2.0, 2.5};
	BrentMethod Sprent1 = {oldspline, xinit3, epsilon, maxi, delta};
	double result3 = Sprent1.operate();
	vec xinit4 = {2.75, 3.0}; 
	BrentMethod Sprent2 = {oldspline, xinit4, epsilon, maxi, delta};
	double result4 = Sprent2.operate();
	std::cout << "Brent Method Root 1: " << result3 << endl;
	std::cout << "Brent Method Root 2: " << result4 << endl;
	return 0;
}
