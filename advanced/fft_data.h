#pragma once
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;
typedef complex<double> dcomp;
typedef vector< complex<double> > vecomp;
typedef double (*func)(double);
typedef vector<double> vec;


class FFT {
private:
	int nmax;
	double L;
public:
	int N;
	dcomp h, dk, Nfac;
	vecomp f_x, F_k;
	vec wvec, posvec;
	bool inverse;
	FFT(vecomp rawdat, int nnn, double length, bool inv) {
		nmax = nnn;
		N = int(pow(2, nnn));
		Nfac.real(1./N);
		Nfac.imag(0.0);
		L = length;
		h.real(L/double(N));
		h.imag(0.0);
		dk.real(1./(h.real()*double(N)));
		dk.imag(0.0);
		f_x = rawdat;
		inverse = inv;
	}
	
	
	void transform() {
		int i2, j1, k1;
		i2 = N >> 1;
		j1 = 0;
		for(int i=0;i<N-1;i++) {
			if(i<j1) {
				dcomp tx = f_x[i];
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
		for(int k=-N/2; k<N/2;k++) {
			wvec.push_back(k*dk.real());
			posvec.push_back(k*h.real());
			// Make vector of rearranged function times the appropriate exponential
			vecomp Fouri;
			for(int j(0);j<N;j++){
                                Fouri.push_back(f_x[j]);
                        }
			int N_new = N;
			for(int n=1;n<=nmax; n++) {
				vecomp Fouri_temp;
				N_new = int(pow(2,nmax-n)) ;
				dcomp Wk, theta(2*M_PI*k*N_new/double(N), 0);
				Wk = (inverse)?exp(-I*theta):exp(I*theta);
				
				for(int j=0; j<N_new; j++) {
					Fouri_temp.push_back( Fouri[2*j] +Wk*Fouri[2*j + 1] );
				}
				Fouri.resize(N_new);
				for(int j(0); j<N_new; j++){
					Fouri[j] = Fouri_temp[j];
				}
			}	
			if(inverse){	
				if(k%2==0) {
					F_k.push_back(Nfac*Fouri[0]);
				}	
				else {
					F_k.push_back(-Nfac*Fouri[0]);
				}
			}
			else {
				if(k%2==0) {
					F_k.push_back(h*Fouri[0]);
				}	
				else {
					F_k.push_back(-h*Fouri[0]);
				}
			}
		}
	}

	void invform() {
		f_x = F_k;
		for(int j(0); j<N; j++){
			f_x[j] /= h;
		}
		F_k.resize(0);
		inverse = true;
		transform();
		vecomp f_temp = f_x;
		f_x = F_k;
		F_k = f_temp;
		for(int j(0); j<N; j++) {
			F_k[j] /= Nfac;
		}
		inverse = false;
	}
};










