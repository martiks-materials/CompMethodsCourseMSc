
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

#include "wavepack.h"


void save_psi(int plot_number, vector<double> xdat, vector< complex<double> > psidat, int Num) {
	ostringstream os;
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
	WavePack WP;
	cout << " Quantum Wavepacket Motion" << endl;
	ofstream file("potential.data");
	for (int i = 0; i < WP.N; i++) {
		file << WP.x[i] << "\t" << WP.V(WP.x[i]) << "\n";
	}
	file.close();
	cout << " saved V(x) in file potential.data" << endl;
	save_psi(0, WP.x, WP.psi, WP.N);
	int plots = 999;
	for (int plot = 1; plot <= plots; plot++) {
		double delta_t = 0;
		while (delta_t < WP.tau) {//WP.L/(plots*WP.velocity)) {
			WP.take_step();
			delta_t += WP.tau;
		}
		save_psi(plot, WP.x, WP.psi, WP.N);
	}
	return 0;
}
