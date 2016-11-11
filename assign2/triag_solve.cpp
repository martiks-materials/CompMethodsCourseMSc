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


void triag_solve(double *xdat, double *ydat, double *u, int n, bool natural, double delta_1=0, double delta_n=0){
	double a[n], b[n], c[n], F[n];
	// Cases for i = 1,2,...,n-2. 
        for(int j(1); j<n-1; j++){
                a[j] = (xdat[j] - xdat[j-1])/6;
                b[j] = (xdat[j+1] - xdat[j-1])/3;
                c[j] = (xdat[j+1] - xdat[j])/6;
                F[j] = (ydat[j+1]-ydat[j])/(xdat[j+1]-xdat[j])-(ydat[j]-ydat[j-1])/(xdat[j]-xdat[j-1]);
        }
	// Boundary Conditions
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
        // Cases for Boundary Conditions.
	triag(a, b, c, F, u, n);
}
	
