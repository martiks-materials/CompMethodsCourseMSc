// Computational Methods - Assigment 2: Interpolation
// Martik Aghajanian, Cohort 8.
// This program uses data from the sheet to interpolate values of x within
// the data range using linear interpolation and cubic spline interpolation
// with either natural boundary conditions or fixed first derivative conditions.
#include <iostream>
#include <fstream>
#include <cmath>
#include "interp_tools.h"
using namespace std;

int main(){
	// The "xdat" and "ydat" are arrays containing the known tabulated data in Q1.
	double xdat[] = {-2.1, -1.45, -1.3, -0.2, 0.1, 0.15, 0.8, 1.1, 1.5, 2.8, 3.8};
	double ydat[] = {0.012155, 0.122151, 0.184520, 0.960789, 0.990050, 0.977751, 
			  0.527292, 0.298197, 0.105399, 3.936690e-4, 5.355348e-7};
	int n = 11;
	double delta_1(0), delta_n(0);
	double x[4] = {0.4, -0.128, -2, 3.2};
	// This adds the raw data to a file to be added to a plot for comparison
	// to the interpolated values.
	ofstream outfile1("raw.dat");
	if(!outfile1.is_open()){
		cout << "Error opening raw.dat" << endl;
		return 1;
	}
	for(int i(0); i<n; i++){
		outfile1 << xdat[i] << " " << ydat[i] << endl;
	}
	
	// Uses the linear interpolation function (linear.cpp) on the specific values of x:
	cout << "Linear Interpolation Results:" << endl;
	for(int i(0); i<4; i++){
		cout << "x = " <<  x[i] << ", y = " << linear(x[i], xdat, ydat, n) << endl;
	}	
	

	// This defines two arrays for the second derivatives evaluated at each known point,
	// name "y_nat" and "y_zfd", corresponding to calculations performed with natural 
	// boundary conditions and zero-first-derivative boundary conditions respectively.
	// These are then obtained using triag_solve function which derives the second
	// derivatives from the tabulated known data and imposing continuity of the first
	// derivative.
	double y_nat[n], y_zfd[n];
	triag_solve(xdat, ydat, y_nat, n, true, delta_1, delta_n);
	triag_solve(xdat, ydat, y_zfd, n, false);
	for(int i(0); i<4; i++){
		cout << "Natural: y(x = " << x[i] << ") = " << spline(x[i], xdat, ydat, y_nat, n) << endl; 	
		cout << "1st Deriv. = 0: y(x = " << x[i] << ") = " << spline(x[i], xdat, ydat, y_zfd, n) << endl;
	}


	// For N points equally spaced between the maximum and minimum x values in the known
	// data, both linear and cubic spline (with natural boundary conditions) interpolation
	// is used to compare their effectiveness at representing the function. The variable
	// "dx" denotes the spacing between successive points in the array of x points.
	int N = 1000;
	double dx = (xdat[n-1]-xdat[0])/N;
	double xr[N], y_lin[N], y_cub[N];
	for(int i(0); i<N; i++){
		xr[i] = xdat[0] + i*dx;
		y_cub[i] = spline(xr[i], xdat, ydat, y_nat, n);
		y_lin[i] = linear(xr[i], xdat, ydat, n);
	}	

	// Loads the calculated interpolated values for the N points into data files for plotting.
	ofstream outfile2("linear.dat");
	ofstream outfile3("cubicspline.dat");
	if((!outfile2.is_open())||(!outfile3.is_open())){
		cout << "Error opening..." << endl;
		return 1;
	}
	
	for(int i(0); i<N; i++){
		outfile2 << xr[i] << "\t" << y_lin[i]  <<  endl;
		outfile3 << xr[i] << "\t" << y_cub[i] << endl;
	}

	return 0;

}
	
	
