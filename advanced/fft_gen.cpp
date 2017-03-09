// Computational Methods - Advanced Topics: Fast Fourier Transforms
// 
// Class member functions for the FFT algorithm to compute for either 
// a complex data set sampling, or to receive a complex function
// which the algorithm can sample for the user.
#include "fft_gen.h"
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;
typedef complex<double> dcomp;
typedef vector< complex<double> > vecomp;
typedef double (*func)(double);
typedef vector<double> vec;


FFT::FFT(int nnn, double length, bool inv) {
	nmax = nnn;
	N = int(pow(2, nnn));
	Nfac.real(1./N);
	Nfac.imag(0.0);
	L = length;
	h.real(L/double(N));
	h.imag(0.0);
	dk.real(1./(h.real()*double(N)));
	dk.imag(0.0);
	inverse = inv;
}
	
void FFT::init_data(vecomp rawdat) {
	// Initialise the data to be Fourier transformed using a given dataset 'rawdat.'
	for(int j=-N/2; j<N/2; j++){
		posvec.push_back(j*h.real());
	}
	f_x = rawdat;
}
	
void FFT::init_func(func fun) {
	// Initialise the data to be Fourier transformed using a given function 'fun'.
	for(int j=-N/2; j<N/2+1; j++){
		posvec.push_back(j*h.real());
		f_x.push_back((*fun)(double(j)*h.real()));
	}		
}	
	
void FFT::transform() {
	// Perform the FWD FFT.
	int i2, j1, k1;
	i2 = N >> 1;
	j1 = 0;
	for(int i=0;i<N-1;i++) {
		if(i<j1) {
			// This performs the switch if the bit-reversed label is smaller.
			dcomp tx = f_x[i];
			f_x[i] = f_x[j1];
			f_x[j1] = tx;
		}
		k1 = i2;
		// Bit reverse the labels of the data.
		while(k1 <= j1) {
			j1 -= k1;
			k1 >>= 1;
		}
		j1 += k1;
	}
	// For each k inbetween the positive and negative Nyquist frequencies, transform.
	dcomp I(0.0, 1.0);
	for(int k=-N/2; k<N/2;k++) {
		// Record the wavevector element being transformed.
		wvec.push_back(k*dk.real());
		// 'Fouri' is a temporary array which successively reduces in size, and is
		// no initialised as a copy of the bit-reversed labels data.
		vecomp Fouri;
		for(int j(0);j<N;j++){
			Fouri.push_back(f_x[j]);
		}
		int N_new = N;
		// Halve the dataset and perform complex multiplication on the odd part of the data.
		for(int n=1;n<=nmax; n++) {
			vecomp Fouri_temp;
			N_new = int(pow(2,nmax-n));
			// Exponential factor we must apply that depends on k. Sign depends on 
			// bool 'inverse'
			dcomp Wk, theta(2*M_PI*k*N_new/double(N), 0);
			Wk = (inverse)?exp(-I*theta):exp(I*theta);
			
			for(int j=0; j<N_new; j++) {
				Fouri_temp.push_back( Fouri[2*j] +Wk*Fouri[2*j + 1] );
			}
			// Resize the array to half and update according to 'Fouri_temp'.
			Fouri.resize(N_new);
			for(int j(0); j<N_new; j++){
				Fouri[j] = Fouri_temp[j];
			}
		}
		// Once the array has been reduced to the whole FFT (unit size), then update
		// the final list of FTs for each k. The factor of this will depend on whether
		// or not we will be computing the inverse.
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
	
void FFT::invform() {
	// Modifies the input and output of the FFT algorithm to compute the inverse instead.
	// The input for this will be the previously outputted 'F_k'
	f_x = F_k;
	// Remove the factors of step size from the input before calculating.
	for(int j(0); j<N; j++){
		f_x[j] /= h;
	}
	F_k.resize(0);
	// Activate the inverse bool and go ahead with the transform. 
	inverse = true;
	transform();
	vecomp f_temp = f_x;
	f_x = F_k;
	F_k = f_temp;
	// Apply the factor 1/N to each element when computing the inverse.
	for(int j(0); j<N; j++) {
		F_k[j] /= Nfac;
	}
	// Reset the inverse bool to false if the FT wanted to be computed again.
	inverse = false;
}

