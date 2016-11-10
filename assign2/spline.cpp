#include <iostream>
#include <cmath>
using namespace std;

double spline(double x, double *xdat, double *ydat, double *yderiv, int length){
	int j(0);
        for(int i(0); i<length; i++) {
                if(x<xdat[i]) {
                        j--;
                        break;
                }
                j++;

        }
	
	double A = (xdat[j+1]-x)/(xdat[j+1]-xdat[j]);
	double B = (x-xdat[j])/(xdat[j+1]-xdat[j]);
	double C = (pow(A, 3)-A)*pow((xdat[j+1]-xdat[j]), 2)/6;
	double D = (pow(B, 3)-B)*pow((xdat[j+1]-xdat[j]), 2)/6;

	return (A*ydat[j] + B*ydat[j+1] + C*yderiv[j] + D*yderiv[j+1]);
}
