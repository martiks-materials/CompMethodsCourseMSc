#include <iostream>
#include <vector>
#include <cmath>
#include "rand.h"
#include "mcint.h"
using namespace std;

// These typedefs allow for less messy declaration in this program.
typedef vector<double> vec;
typedef double (*func)(vector<double>);
typedef double (*unifun)(double);


Mcint::Mcint(func fnc, vec a, vec b, int N, int seed, double tol, int Mstart, 
      bool imp, func pdf, unifun cdf, unsigned long int stepmax, bool rej){
	f = fnc;
	p = pdf;
	cy = cdf;
	x_low = a;
	x_up = b;
	Nd = N;
	V=1;
	for(int i(0); i<Nd; i++){
		V *= x_up[i] - x_low[i];
		X.push_back(0);
	}
	myran.seed(seed);
	eps = tol;
	M = 0;
	M_init = Mstart;
	M_max = stepmax;
	importance = imp;
	fsum_old = 0;
	fsum_new = 0;
	fsq_old = 0;
	fsq_new = 0;
	R = 1;
	reject = rej;
	while(M < M_init){
		generator();
		double quant = (*f)(X)/((importance)?((*p)(X)):(1./V));
		fsum_old += quant;
		fsq_old += quant*quant;	
		M++;
	}
	fsum_old /= M;
	fsq_old /= M;
	steps.push_back(M);
	errors.push_back(R);
	vals.push_back(fsum_old); 
}

void Mcint::sample(){
	fsum_new = M*fsum_old;
	fsq_new = M*fsq_old;
	int M_old = M;
	for(int i(0); i<M_old; i++){
		generator();
                double quant = (*f)(X)/((importance)?((*p)(X)):(1./V));
		fsum_new += quant;
		fsq_new += quant*quant;
		M++;
	}
	fsum_new /= M;
	fsq_new /= M;
}

double Mcint::error1(){
	return abs(fsum_new-fsum_old)/abs(fsum_new);
}
	
double Mcint::error2(){
	 return sqrt((fsq_new-pow(fsum_new, 2))/M);
}

void Mcint::integrate(){
	while(R>eps){
		sample();
		R = error2();
		//cout << R << endl;
		fsum_old = fsum_new;
		fsq_old = fsq_new;
		int doublings = (int)log2(M/M_init);
		cout << doublings << " steps: " << fsum_new << endl;
		errors.push_back(R);
		steps.push_back(M);
		vals.push_back(fsum_new);
		if(M>=M_max){
			cout << "Threshold of steps reached. Integration terminated." << endl;
			break;
		}
	}
}

void Mcint::generator(){
	for(int i(0); i<Nd; i++){
		double rx = myran.doub();
		X[i] = x_low[i] + (importance)?((*cy)(rx)):((x_up[i]-x_low[i])*rx);	
	}
}
	



