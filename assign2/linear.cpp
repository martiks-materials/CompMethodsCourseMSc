#include <iostream>
using namespace std;

double linear(double x, double x_dat[], double y_dat[], int length){
        // This function takes known (x, y) data point arrays and for a given
	// value of x, will use linear interpolation to estimate the y value.
	// Function Parameters:
	// 	- x:      (double) Value of x to interpolate
	//  	- xdat:   (double[]) Known x data
	// 	- ydat:   (double[]) Known y data
	// 	- length: (int) Number of known (x, y) data points
	
	int j(0);
	// This finds which interval between known x data points that the give
	// x value resides in. The index we require is "j" and this is the index
	// of the highest x data point for which the given x is more than. 
        for(int i(0); i<length; i++) {
                if(x<x_dat[i]) {
			// If x is less than xdat[i], then 1 must be removed to
			// give the index j of the closest data point less than x. 
                        j--;
                        break;
                }
                j++;

        }

	// Using the correct index j, calculates the linear interpolation coefficients. 
        double A = (x_dat[j+1] - x)/(x_dat[j+1] - x_dat[j]);
        double B = (x - x_dat[j])/(x_dat[j+1] - x_dat[j]);

        return (A*y_dat[j]) + (B*y_dat[j+1]);
}
