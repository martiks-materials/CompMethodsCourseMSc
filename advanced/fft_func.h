#pragma once
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;
typedef complex<double> dcomp;
typedef double (*func)(double);
typedef vector<double> vec;


class FFT {
private:
	func f;
	int nmax;
	double L;
public:
	int N;
	double h, dk;
	vec f_x, F_kr, F_ki, posvec, wvec; 
	FFT(func fun, int nnn, double length) {
		nmax = nnn;
		N = int(pow(2, nnn));
		L = length;
		h = L/double(N);
		dk =  1./(h*double(N));
		f = fun;
		for(int j=-N/2; j<N/2+1; j++){
			posvec.push_back(j*h);
			f_x.push_back((*f)(double(j)*h));
		}	
	}
	
	
	void transform() {
		int i2, j1, k1;
		i2 = N >> 1;
		j1 = 0;
		for(int i=0;i<N-1;i++) {
			if(i<j1) {
				double tx = f_x[i];
				f_x[i] = f_x[j1];
				f_x[j1] = tx;
			}
			k1 = i2;
			while(k1 <= j1) {
				j1 -= k1;
				k1 >>= 1;
			}
			j1 += k1;
		}
		// For each k inbetween the positive and negative Nyquist Frequencies, do transform
		dcomp I(0.0, 1.0);
		for(int k=-N/2+1; k<N/2+1;k++){
			wvec.push_back(k*dk);
			// Make vector of rearranged function times the appropriate exponential
			vector<dcomp> Fouri;
			for(int j(0);j<N;j++){
				dcomp fx(f_x[j], 0);
				Fouri.push_back(fx);
			}
			int N_new = N;
			for(int n=1;n<=nmax; n++){
				vector<dcomp> Fouri_temp;
				N_new = int(pow(2,nmax-n)) ;
				dcomp Wk, theta(2*M_PI*k*N_new/double(N), 0);
				Wk = exp(I*theta);
				
				for(int j=0; j<N_new; j++) {
					Fouri_temp.push_back( Fouri[2*j] +Wk*Fouri[2*j + 1] );
				}
				Fouri.resize(N_new);
				for(int j(0); j<N_new; j++){
					Fouri[j] = Fouri_temp[j];
				};
			}
			if(k%2==0) {
				F_kr.push_back(h*Fouri[0].real());
				F_ki.push_back(h*Fouri[0].imag());
			}	
			else {
				F_kr.push_back(-h*Fouri[0].real());
				F_ki.push_back(-h*Fouri[0].imag());
			}
		}
	}
};










