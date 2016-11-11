#include <iostream>
#include <cmath>
using namespace std;

double spline(double x, double *xdat, double *ydat, double *yderiv, int length){
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



	// This finds which interval between known x data points that the give
        // x value resides in. The index we require is "j" and this is the index
        // of the highest x data point for which the given x is more than. 
	int j(0);
        for(int i(0); i<length; i++) {
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
