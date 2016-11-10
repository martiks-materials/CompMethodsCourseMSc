#include <iostream>
#include <cmath>

using namespace std;

void triag(double *a, double *b, double *c, double *F, double *u, int n){
        // Adaptation of the tridag routine from numerical recipes
        // where a, b and c are the sets of coefficients, F is the
        // vector of functions and u is the vector being solved for.

        double beta=b[0];
        double gamma[n];
        double eps = 1e-10;
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


void triag_solve(double *xdat, double *ydat, double *u, int n, bool natural, double delta){
	int shift = (natural)? 1 : 0;
	int ny = (natural)? n-2 : n;
	double a[ny], b[ny], c[ny], F[ny];
	double delta_1(delta), delta_n(delta);
	if(!natural){
		a[0] = 0;
	        b[0] = (xdat[1]-xdat[0])/3;
        	c[0] = b[0]/2;
        	F[0] = (ydat[1]-ydat[0])/(xdat[1]-xdat[0])-delta_1;
        	c[n-1] = 0;
        	b[n-1] = (xdat[n-1]-xdat[n-2])/3;
        	a[n-1] = b[n-1]/2;
        	F[n-1] = delta_n - (ydat[n-1]-ydat[n-2])/(xdat[n-1]-xdat[n-2]);
	}
        // Cases for Boundary Conditions.
       	// Cases for i = 1,2,...,n-2. 
        for(int j(1); j<n-1; j++){
		a[j-shift] = (xdat[j] - xdat[j-1])/6;
         	b[j-shift] = (xdat[j+1] - xdat[j-1])/3;
                c[j-shift] = (xdat[j+1] - xdat[j])/6;
                F[j-shift] = (ydat[j+1]-ydat[j])/(xdat[j+1]-xdat[j])-(ydat[j]-ydat[j-1])/(xdat[j]-xdat[j-1]);
        }
	triag(a, b, c, F, u, ny);
}
	
