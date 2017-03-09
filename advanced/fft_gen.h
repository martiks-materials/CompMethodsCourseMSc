// Computational Methods - Advanced Topics: Fast Fourier Transforms
// 
// Class declaration for the FFT algorithm to compute for either 
// a complex data set sampling, or to receive a complex function
// which the algorithm can sample for the user.

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
	// nmax - Power of 2 for which the dataset contains 2^{nmax} points.
	// L    - Physical length of the spatial domain. 
	int nmax;
	double L;
public:
	// N      - Number of points in the dataset
	// h      - Step between adjacent points (Complex to simplify multiplication)
	// dk     - Frequency domain step (Complex to -\\-)
	// Nfac   - Factor used in converting sum to inverse FFT (Complex to -\\-)
	// wvec   - Array of wavevectors for plotting
	// posvec - Array of positions for plotting
	// inverse- Boolean to signify the inverse FFT
	int N;
	dcomp h, dk, Nfac;
	vecomp f_x, F_k;
	vec wvec, posvec;
	bool inverse;

	FFT(int nnn, double length, bool inv);
	
	void init_data(vecomp rawdat);

	void init_func(func fun);

	
	void transform();

	void invform();

};
