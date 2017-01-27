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
	return 1.5*exp(-1.5*y[0])/(1-pow(M_E, -3));
}

double invcdf2(double x){
	return -(2./3)*log(1-((1-pow(M_E,-3))*x));
}

int main() {
	int Ndim(1), mstart(500), see(142), maxsteps(22);
	unsigned long int mmax = int(mstart*pow(2, maxsteps));
	double eps(1E-6);
        vec a = {0}, b = {2};
        Mcint Monte1 = {erfi, a, b, Ndim, see, eps, mstart, false, pdf, invcdf, mmax};
        Mcint Monte2 = {erfi, a, b, Ndim, see, eps, mstart, true, pdf, invcdf, mmax};
        Mcint Monte3 = {erfi, a, b, Ndim, see, eps, mstart, true, pdf2, invcdf2, mmax};
        Monte1.integrate();
	Monte2.integrate();
	Monte3.integrate();
	ofstream outfile1a("ertrack1.dat");
	ofstream outfile1b("numtrack1.dat");
	ofstream outfile2a("ertrack2.dat");
	ofstream outfile2b("numtrack2.dat");
	ofstream outfile3a("ertrack3.dat");
	ofstream outfile3b("numtrack3.dat");
	for(int i(1); i <= maxsteps; i++){
		outfile1a << log2(Monte1.steps[i]/mstart) << "   " << log(Monte1.errors[i]) << endl;
		outfile1b << log2(Monte1.steps[i]/mstart) << "   " << Monte1.vals[i] << endl;
		outfile2a << log2(Monte2.steps[i]/mstart) << "   " << log(Monte2.errors[i]) << endl;
		outfile2b << log2(Monte2.steps[i]/mstart) << "   " << Monte2.vals[i] << endl;
		outfile3a << log2(Monte3.steps[i]/mstart) << "   " << log(Monte3.errors[i]) << endl;
		outfile3b << log2(Monte3.steps[i]/mstart) << "   " << Monte3.vals[i] << endl;
	}
	outfile1a.close();
	outfile1b.close();
	outfile2a.close();
	outfile2b.close();
	outfile3a.close();
	outfile3b.close();
        cout << "Uniform Sampling: " << Monte1.fsum_new << " +/- " << Monte1.error2() <<  endl;
	cout << "Importance Sampling (PDF1): " << Monte2.fsum_new << " +/- " << Monte2.error2() <<  endl;
	cout << "Importance Sampling (PDF2): " << Monte3.fsum_new << " +/- " << Monte3.error2() <<  endl;
        return 0;
}
