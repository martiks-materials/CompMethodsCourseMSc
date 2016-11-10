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

	
