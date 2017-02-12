#include <iostream>
#include <cmath>
#include "splining.h"
using namespace std;
typedef double (*func)(double x);
	
double oldspline(double x) { //, double *xdat, double *ydat, double *yderiv, int length){
	// Given known (x, y) data points and a calculated set of second derivatives,
	// performs cubic spline interpolation, returning the estimate of the function.
	// The boundary conditions for the spline  are already contained within the 
	// array of second derivatives and this function just produces the result to 
	// avoid recalculation of the coefficients and solving of the matrix used
	// to obtain the second derivatives.
	// Function parameters:
	// 	- x:	  (double)   Value to compute using the cubic spline.
	// 	- xdat:   (double[]) Known x data points.
	// 	- ydat:   (double[]) Known y data points.
	// 	- yderiv: (double[]) Calculated second derivatives array.
	// 	- length: (int)	     Number of (x, y) known data points.
	
	double xdat[] = {-2.1, -1.45, -1.3, -0.2, 0.1, 0.15, 0.8, 1.1, 1.5, 2.8, 3.8};
	double ydat[] = {0.012155, 0.122151, 0.184520, 0.960789, 0.990050, 0.977751, 
			  0.527292, 0.298197, 0.105399, 3.936690e-4, 5.355348e-7};
	int n = 11;
	double yderiv[n];
	triag_solve(xdat, ydat, yderiv, n, false, 0, 0);
	// This finds which interval between known x data points that the give
        // x value resides in. The index we require is "j" and this is the index
        // of the highest x data point for which the given x is more than. 
	int j(0);
        for(int i(0); i<n; i++) {
                if(x<xdat[i]) {
			// If x is less than xdat[i], then 1 must be removed to
                        // give the index j of the closest data point less than x. 
                        j--;
                        break;
                }
                j++;

        }
	
        // Using the correct index j, calculates the cubic spline interpolation coefficients. 	
	double A = (xdat[j+1]-x)/(xdat[j+1]-xdat[j]);
	double B = (x-xdat[j])/(xdat[j+1]-xdat[j]);
	double C = (pow(A, 3)-A)*pow((xdat[j+1]-xdat[j]), 2)/6;
	double D = (pow(B, 3)-B)*pow((xdat[j+1]-xdat[j]), 2)/6;

	return (A*ydat[j] + B*ydat[j+1] + C*yderiv[j] + D*yderiv[j+1]);
}




void triag(double *a, double *b, double *c, double *F, double *u, int n){
        // Adaptation of the tridag routine from Numerical Recipes 3rd Edition,
	// which takes a tridiagonal matrix M defined by a[n], b[n] and c[n] arrays,
	// and an inhomogenous vector to solve a set of linear equations for 
	// solution vector u[n].
	// Function parameters:
	// 	- a, b, c: (double[]) Coefficients for tridiagonal matrix M.
	//	- u:       (double[]) Solution array which is modified not returned.
	// 	- F:	   (double[]) Inhomogenous vector array such that Mu = F.
	// 	- n: 	   (int)      Length of the arrays used, dimensions of matrix M.

        double beta=b[0];
        double gamma[n];

	// Defines an epsilon to compare double precision floating points to 0.
        double eps(1e-10);
        if(abs(b[0])<eps){
                cout << "Error 1.\n";
        }

        u[0] = F[0]/beta;
        for(int j(1); j<n; j++){
                gamma[j] = c[j-1]/beta;
                beta = b[j]-a[j]*gamma[j];
                if(abs(beta)<eps){
                        cout << "Error 2.\n";
                }
                u[j] = (F[j]-a[j]*u[j-1])/beta;
        }
        for(int j=n-2; j>=0; j--){
                u[j] -= gamma[j+1]*u[j+1];
        }
}


void triag_solve(double *xdat, double *ydat, double *u, int n, bool natural, double delta_1, double delta_n){
	// From given known x and y data points, derives a tridiagonal matrix arising from the
	// master equation for cubic spline interpolation, derived from imposing continuity of
	// the first derivative at each data point. Boundary conditions can either be natural,
	// where the second derivative is zero at the boundaries, or by specifying the first
	// derivative at the boundaries via the "delta_1/n" variables. This matrix is then solved
	// and modifies an (empty) solution vector u[n], a set of second derivatives at each point.
	// Function parameters (not specified in the triag function):
	// 	- xdat:    (double[]) Known x data points.
	// 	- ydat:    (double[]) Known y data points.
	// 	- natural: (bool)     Turns on natural boundary conditions if true.
	//	- delta_1: (double)   First derivative at the lower boundary.
	//	- delta_n: (double)   First derivative at the upper boundary.
	

	// Constructs the matrix elements for the bulk of the interval (i = 1, 2,... , n-2).
	double a[n], b[n], c[n], F[n];
        for(int j(1); j<n-1; j++){
                a[j] = (xdat[j] - xdat[j-1])/6;
                b[j] = (xdat[j+1] - xdat[j-1])/3;
                c[j] = (xdat[j+1] - xdat[j])/6;
                F[j] = (ydat[j+1]-ydat[j])/(xdat[j+1]-xdat[j])-(ydat[j]-ydat[j-1])/(xdat[j]-xdat[j-1]);
        }
	// Boundary cases specified by the "natural" bool variable.
	if(natural){
		a[0]=0;
		a[1]=0;
		b[0]=1;
		c[0]=0;
		F[0]=0;
		a[n-1]=0;
		b[n-1]=1;
		c[n-1]=0;
		c[n-2]=0;
		F[n-1]=0;
	} else{
		a[0] = 0;
	        b[0] = (xdat[1]-xdat[0])/3;
        	c[0] = b[0]/2;
        	F[0] = (ydat[1]-ydat[0])/(xdat[1]-xdat[0])-delta_1;
        	c[n-1] = 0;
        	b[n-1] = (xdat[n-1]-xdat[n-2])/3;
        	a[n-1] = b[n-1]/2;
        	F[n-1] = delta_n - (ydat[n-1]-ydat[n-2])/(xdat[n-1]-xdat[n-2]);
	}
        // Inputs this into the triag function to correctly modify the solution vector u.
	triag(a, b, c, F, u, n);
}


