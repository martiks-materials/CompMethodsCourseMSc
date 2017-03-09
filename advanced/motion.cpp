// Computational Methods - Advanced Topics: Fast Fourier Transforms
//
// Create data of the wavefunction for each timestep in a simulation
// of a Gaussian wavepacket propagating towards a potential barrier.
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "fft_gen.h"
using namespace std;

double V(double x, double V_0, double V_center, double V_width) {
	// Function regarind the potential used for this wavepacket simulation.
	// Basically a barrier.
	double half_width = abs(0.5 * V_width);
	if (abs(x - V_center) <= half_width){
	     return V_0;
	}
	else {
		return 0; 
	}
}

void save_psi(int plot_number, vector<double> xdat, vector< complex<double> > psidat, int Num) {
	// Plot the data into a file for a given plot number for labelling, using an array of spatial
	// data and a complex vector of wavefunction elements with a certain 'Num' of elements.
	ostringstream os;
	// Label data based on plot number (Note directories need to exist first, so if they do not 
	// just create them)
	os << "Data/psi_" << plot_number << ".data"; string file_name(os.str());
	ofstream file(file_name.c_str());
	if (file){
		cout << " writing " << file_name << endl;
	}
	else {
		cerr << " cannot open " << file_name << endl;
	}
	for (int i = 0; i < Num; i++){
		file << xdat[i] << "\t" << norm(psidat[i]) << "\t" << psidat[i].real() << "\t" << psidat[i].imag() << "\n";
	}
	file.close();
}

int main() {	
	// N 			- Number of points in spatial domain
	// nmax			- Power of two such that N=2^{nmax}
	// h			- Spatial increment
	// h_bar		- Fundamental constant set to 1
	// mass 		- Particle mass set to 1
	// L			- Physical length of system
	// tau			- Time increment
	// t			- Updatable time
	// V_0			- Potential amplitude
	// V_width		- Potential width
	// V_center		- Potential centre
	// E			- Energy of the wavepacket, set to 1
	// x_0			- Mean position of wavepacket
	// psi_norm		- Norm of the wavefunction vector
	// sigma_0		- Width of wavepacket
	// k_0			- Central wavevector of gaussian (initial momentum)
	// pi			- Fundamental constant
	// x			- Data storage of the position for plotting
	// psi 			- Our vector approximation to the wavefunction
	// T_exp_factor		- Kinetic energy exponential
	// V_exp_factor		- Potential energy exponential
	int nmax = 9;
	int N = int(pow(2, nmax));
	double L = 100;
	double h = L/double(N), h_bar = 1, mass = 1, tau = 0.01, t = 0;
	double V_0 = 0.3, V_width = 10, V_center = 0.75*L, E = 1, x_0 = L/4;
	double psi_norm, k_0, sigma_0 = L/10., pi = 4*atan(1.0);
	vector<double> x;	
	vector< complex<double> > psi, T_exp_factor, V_exp_factor;
	// Reset all vectors and the lattice
	for (int j = 0; j < N; j++){
		x.push_back(j*h);
	}
	// Inititalize the packet the wavefunction packet
	k_0 = sqrt(2*mass*E - h_bar*h_bar/(2*sigma_0*sigma_0))/h_bar;
	psi_norm = 1./sqrt(sigma_0*sqrt(pi));
	for (int j = 0; j < N; j++) {
	    double exp_factor=exp(-(x[j] - x_0)*(x[j] - x_0)/(2*sigma_0*sigma_0));
	    psi.push_back( complex<double>(psi_norm*cos(k_0*x[j])*exp_factor, psi_norm*sin(k_0*x[j])*exp_factor));
	}
	// Initialize the phase rotation factors
	for (int j = 0; j < N; j++) {
		double p = ( j < N/2)? j : (j - N);
		p *= h_bar*2*pi/L;
		double theta = -p*p/(2*mass)/h_bar*tau;
		// Initialise the kinetic energy terms.
		complex <double> kin(cos(theta), -sin(theta));
		T_exp_factor.push_back(kin);
		theta = -V(x[j], V_0, V_center, V_width)/(2*h_bar*tau);
		// Initialise the potential energy terms.
		complex<double> pot(cos(theta), -sin(theta));
		V_exp_factor.push_back(pot);
	}
	cout << "Quantum Wavepacket Motion" << endl;
	ofstream file("potential.data");
	for (int i = 0; i < N; i++) {
		file << x[i] << "\t" << V(x[i], V_0, V_center, V_width) << "\n";
	}
	file.close();
	cout << " saved V(x) in file potential.data" << endl;
	save_psi(0, x, psi, N);
	int plots = 2499;
	// For each plot, take a time step of the wavepacket motion and update the time
	// saving the data in a separate file to the storage directory.
	for (int plot = 1; plot <= plots; plot++) {
		double delta_t = 0;
		while (delta_t < tau) {
			for (int j(0); j < N; j++){
				psi[j] *= V_exp_factor[j];
			}
			// Fourier transform into momentum space.
			FFT fft = {nmax, L, false };
			fft.init_data(psi);
			fft.transform();
			// Now we are in momentum space, multiply by the quadratic KE operator exponential term.
			for (int j(0); j < N; j++){
				fft.F_k[j] *= T_exp_factor[j];
			}
			// Inverse FFT to get back to the wavefunction so we can apply the second square root
			// of the potential energy part.
			fft.invform();
			psi = fft.f_x;
			for (int j(0); j < N; j++){
				psi[j] *= V_exp_factor[j];
			}
			// Update the time variable.
			t += tau; 
			delta_t += tau;
		}
		save_psi(plot, x, psi, N);
	}
	return 0;
}
