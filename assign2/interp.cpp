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
	
	int n = 11, start = 0, shift = 1;
	bool natural = true;
	double a[n], b[n], c[n], F[n], y[n];
	double delta_1(0), delta_n(0);
	// Cases for Boundary Conditions.
	if(natural==true){  	
	} else{
		start += 1;
		shift -= 1;
		a[0] = 0;        
		b[0] = (xdat[1]-xdat[0])/3;	 
		c[0] = b[0]/2;		   
		F[0] = (ydat[1]-ydat[0])/(xdat[1]-xdat[0])-delta_1;    
		c[n-1] = 0;
		b[n-1] = (xdat[n-1]-xdat[n-2])/3;
		a[n-1] = b[n-1]/2;
		F[n-1] = delta_n - (ydat[n-1]-ydat[n-2])/(xdat[n-1]-xdat[n-2]);
	}
	
	// Cases for i = 1,2,...,n-2. 
	for(int j(start); j<n-1; j++){
		int i = j + shift;
		a[j] = (xdat[i] - xdat[i-1])/6;
		b[j] = (xdat[i+1] - xdat[i-1])/3;
		c[j] = (xdat[i+1] - xdat[i])/6;
		F[j] = (ydat[i+1]-ydat[i])/(xdat[i+1]-xdat[i])-(ydat[i]-ydat[i-1])/(xdat[i]-xdat[i-1]);
	};
	triag(a, b, c, F, y, n);
	
	//triag_solve(xdat, ydat, y, n, natural);
	for(int j(0); j<n; j++){
		cout << y[j] << ", ";
	}
	cout << endl << endl;
	
	// For natural splines
	double y_dd[n];
	if(natural==true){
		y_dd[0] = 0;
		for(int i(1); i<n-1; i++){
			y_dd[i] = y[i-1];
		}
	} else {
		for(int i(0); i<n; i++){
			y_dd[i] = y[i];
		}
	}
	double x[4] = {0.4, -0.128, -2, 3.2};
	for(int i(0); i<4; i++){
		cout << spline(x[i], xdat, ydat, y_dd, n) << endl; 	
	}

	// Now to get data for plotting
	int N = 500;
	double dx = (xdat[n-1]-xdat[0])/N;
	double xr[N], y_lin[N], y_cub1[N];
	for(int i(0); i<N; i++){
		xr[i] = xdat[0] + i*dx;
		y_cub1[i] = spline(xr[i], xdat, ydat, y_dd, n);
		y_lin[i] = linear(xr[i], xdat, ydat, 11);
	}


	ofstream outfile2("linear.dat");
	ofstream outfile3("cubspline.dat");
	if((!outfile2.is_open())||(!outfile3.is_open())) {
		cout << "Error opening..." << endl;
		return 1;
	}
	
	for(int i(0); i<N; i++){
		outfile2 << xr[i] << "\t" << y_lin[i]  <<  endl;
		outfile3 << xr[i] << "\t" << y_cub1[i] << endl;
	}

	return 0;

}
	
	
