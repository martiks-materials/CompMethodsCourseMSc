// Computational Methods - Assigment 2: Interpolation
// Martik Aghajanian, Cohort 8.
#include <iostream>
#include <fstream>
#include <cmath>
#include "interp_tools.h"
using namespace std;

int main(){
	// The "xdat" and "ydat" are arrays containing the tabulated data in Q1.
	double xdat[] = {-2.1, -1.45, -1.3, -0.2, 0.1, 0.15, 0.8, 1.1, 1.5, 2.8, 3.8};
	double ydat[] = {0.012155, 0.122151, 0.184520, 0.960789, 0.990050, 0.977751, 
			  0.527292, 0.298197, 0.105399, 3.936690e-4, 5.355348e-7};
	ofstream outfile1("raw.dat");
	if(!outfile1.is_open()){
		cout << "Error opening raw.dat" << endl;
		return 1;
	}
	for(int i(0); i<11; i++){
		outfile1 << xdat[i] << " " << ydat[i] << endl;
	}

	cout << "Linear Interpolation Results:" << endl;	
	cout << "x = 0.4: y = " << linear(0.4, xdat, ydat, 11) << endl;
	cout << "x = -0.128: y = " << linear(-0.128, xdat, ydat, 11) << endl;
	cout << "x = -2.0: y = " << linear(-2.0, xdat, ydat, 11) << endl;
	cout << "x = 3.2: y = " << linear(3.2, xdat, ydat, 11) << endl;
	
	int n = 11;
	double delta(0);
	double y_nat[n], y_zfd[n];
	triag_solve(xdat, ydat, y_nat, n, true, delta);
	triag_solve(xdat, ydat, y_zfd, n, false, delta);
	
	// For natural splines
	double y2_nat[n] = {0};
	for(int i(1); i<n-1; i++){
		y2_nat[i] = y_nat[i-1];
	}
	double x[4] = {0.4, -0.128, -2, 3.2};
	for(int i(0); i<4; i++){
		cout << "Natural: y(x = " << x[i] << ") = " << spline(x[i], xdat, ydat, y2_nat, n) << endl; 	
		cout << "1st Deriv. = 0: y(x = " << x[i] << ") = " << spline(x[i], xdat, ydat, y_zfd, n) << endl;
	}


	// Now to get data for plotting
	int N = 1000;
	double dx = (xdat[n-1]-xdat[0])/N;
	double xr[N], y_lin[N], y_cub1[N], y_cub2[N];
	for(int i(0); i<N; i++){
		xr[i] = xdat[0] + i*dx;
		y_cub1[i] = spline(xr[i], xdat, ydat, y2_nat, n);
		y_cub2[i] = spline(xr[i], xdat, ydat, y_zfd, n);
		y_lin[i] = linear(xr[i], xdat, ydat, 11);
	}	


	ofstream outfile2("linear.dat");
	ofstream outfile3("nat_spline.dat");
	ofstream outfile4("zero_spline.dat");
	if((!outfile2.is_open())||(!outfile3.is_open())||(!outfile4.is_open())){
		cout << "Error opening..." << endl;
		return 1;
	}
	
	for(int i(0); i<N; i++){
		outfile2 << xr[i] << "\t" << y_lin[i]  <<  endl;
		outfile3 << xr[i] << "\t" << y_cub1[i] << endl;
		outfile4 << xr[i] << "\t" << y_cub2[i] << endl;
	}

	return 0;

}
	
	
