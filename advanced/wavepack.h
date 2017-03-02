#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

#include "fft_data.h"

class WavePack {
private:
public:
	int N, nmax;
	double h, h_bar, mass, L, tau, t,  V_0, V_width, V_center; 	
	double E, x_0, psi_norm, sigma_0, k_0, velocity, pi;
	double T, frames_per_second;
	vector<double> x;	
	vector< complex<double> > psi, T_exp_factor, V_exp_factor;   // precomputed phase rotations
	WavePack() {
		T = 5;
		frames_per_second = 50;
		t = 0;
		h_bar = 1;
		mass = 1;
		nmax = 7;
		N = int(pow(2, nmax));
		L = 100;
		h = L/double(N);
		tau = 0.1;
		V_0 = 0.5;
		V_width = 10;
		V_center = 0.75*L;
		pi = 4*atan(1.0);
		x_0 = L/ 4; 
		E = 1;
		sigma_0 = L/10.;
		// reset vectors
		// reset the lattice
		h = L/double(N);
		for (int j = 0; j < N; j++){
			x.push_back(j*h);
		}
		// inititalize the packet
		k_0 = sqrt(2*mass*E - h_bar*h_bar/(2*sigma_0*sigma_0))/h_bar;
		velocity = k_0/mass;
		psi_norm = 1./sqrt(sigma_0*sqrt(pi));
		for (int j = 0; j < N; j++) {
		    double exp_factor=exp(-(x[j] - x_0)*(x[j] - x_0)/(2*sigma_0*sigma_0));
		    psi.push_back( complex<double>(psi_norm*cos(k_0*x[j])*exp_factor, psi_norm*sin(k_0*x[j])*exp_factor));
		}
		// initialize the phase rotation factors
		for (int j = 0; j < N; j++) {
			// kinetic factor exp[-iT/h_bar tau]
			double p = ( j < N/2)? j : (j - N);
			p *= h_bar*2*pi/L;
			double theta = -p*p/(2*mass)/h_bar*tau;
			complex <double> kin(cos(theta), -sin(theta));
			T_exp_factor.push_back(kin);
			// potential factor exp[-iV(x)/(2h_bar) tau]
			theta = -V(x[j])/(2*h_bar*tau);
			complex<double> pot(cos(theta), -sin(theta));
			V_exp_factor.push_back(pot);
		}
		t = 0; 
	}


	double V(double x) {
		double half_width = abs(0.5 * V_width);
		if (abs(x - V_center) <= half_width){
		     return V_0;
		}
		else {
			return 0; 
		}
	}


	void take_step() {
		// first half potential phase rotation
		for (int j(0); j < N; j++){
			psi[j] *= V_exp_factor[j];
		}	// FFT to momentum space
		FFT fft = {psi, nmax, L, false };
		fft.transform();
		// kinetic phase rotation
		for (int j(0); j < N; j++){
			fft.F_k[j] *= T_exp_factor[j];
		}
		// FFT back to position space
		fft.invform();
		psi = fft.f_x;
		// second half potential phase rotation
		for (int j(0); j < N; j++){
			psi[j] *= V_exp_factor[j];
		}
		t += tau; 
	}

/*
	void time_step() {
		static clock_t clock_start;
		static bool done;
		if (!done) {
			double t0 = t;
			do {
				take_step();
			} while (abs(velocity * (t - t0)) < L / T / frames_per_second);
			done = true;
		}
		clock_t clock_now = clock();
		double seconds = (clock_now - clock_start) / double(CLOCKS_PER_SEC);
		if ( seconds < 1 / frames_per_second ) {
			return;
		} 
		else {
			clock_start = clock_now;
		done = false;
		}
		glutPostRedisplay();
		glFlush(); 
	}
*/

};

