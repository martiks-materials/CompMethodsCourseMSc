#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "mcint.h"
using namespace std;


double erfi(vector<double> x){
	return (2./sqrt(M_PI))*exp(-x[0]*x[0]);
}

double pdf(vector<double> y){
        return 0.98-(0.48*y[0]);
}

double invcdf(double x){
	double X = (0.98-sqrt((0.98*0.98)-(0.96*x)))/0.48;
	return X;
}

double pdf2(vector<double> y){
	return 2*exp(-2*y[0])/(1-pow(M_E, -4));
}

double invcdf2(double x){
	return -0.5*log(1-((1-pow(M_E,-4))*x));
}

int main() {
	bool importan;
	cout << "Enter 0 for Uniform, 1 for Importance Sampling: ";
	cin >> importan;
	cout << endl;
	int Ndim(1), mstart(500), see(142), maxsteps(22);
	unsigned long int mmax = int(mstart*pow(2, maxsteps));
	double eps(1E-6);
        vec a = {0}, b = {2};
        Mcint Monte = {erfi, a, b, Ndim, see, eps, mstart, importan, pdf, invcdf, mmax};
        Monte.integrate();
	if(importan){
		ofstream outfile1("ertrack2.dat");
		ofstream outfile2("numtrack2.dat");
		for(int i(0); i <= maxsteps; i++){
			outfile1 << log2(Monte.steps[i]/mstart) << "   " << Monte.errors[i] << endl;
			outfile2 << log2(Monte.steps[i]/mstart) << "   " << Monte.vals[i] << endl;
		}
		outfile1.close();
		outfile2.close();
        	cout << "Importance Sampling: " << Monte.fsum_new << " +/- " << Monte.error2() <<  endl;
	}
	else {
		ofstream outfile1("ertrack1.dat");
		ofstream outfile2("numtrack1.dat");
		for(int i(0); i <= maxsteps; i++){
			outfile1 << log2(Monte.steps[i]/mstart) << "   " << Monte.errors[i] << endl;
			outfile2 << log2(Monte.steps[i]/mstart) << "   " << Monte.vals[i] << endl;
		}	
		outfile1.close();
		outfile2.close();
        	cout << "Uniform Sampling: " << Monte.fsum_new << " +/- " << Monte.error2() <<  endl;
	}
        return 0;
}
