// Computational Methods - Advanced Topics: Fast Fourier Transforms.
//
// Using data obtained from an audiofile, find the FFT and outputs it
// to a separate directory for plotting
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "fft_gen.h"
using namespace std;

void save_psi(int plot_number, bool space, bool first,  vector<double> xdat, vector< complex<double> > psidat, int Num) {
	ostringstream os;
	double thresh = 25;
	if(space){
		if(first) {
			os << "seepdat/space/seepspace1/space_" << plot_number << ".data";
		}
		else {
			os << "seepdat/space/seepspace2/space_" << plot_number << ".data";
		}
	}
	else {
		if(first) {
			os << "seepdat/mom/seepmom1/mom_" << plot_number << ".data";
		}
		else {
			os << "seepdat/mom/seepmom2/mom_" << plot_number << ".data";
		}
	}
	string file_name(os.str());
	ofstream file(file_name.c_str());
	if (file){
		cout << " writing " << file_name << endl;
	}
	else {
		cerr << " cannot open " << file_name << endl;
	}
	for (int i = 0; i < Num; i++){
		if((space==false)&&(norm(psidat[i]) < thresh)) {
			file << xdat[i] << "\t" << 0 << "\t" << psidat[i].real() << "\t" << psidat[i].imag() << "\n";
		}
		else {	
			file << xdat[i] << "\t" << norm(psidat[i]) << "\t" << psidat[i].real() << "\t" << psidat[i].imag() << "\n";
		}
	}
	file.close();
}


int main() {
	// Initial Parameters for sampling.
	int plots = 19523;
	int powmax = 9;
	int N = int(pow(2, powmax));
	double L = double(N);
	cout << "L = " << L << endl;
	double h = L/double(N);
	int count = 0;
	vector<double> fulldat;	
	vector<double> fulldat2;
	string input_line;
	double samp1, samp2;
	ifstream infile;
	infile.open("dragon/allmusic.dat");
	while (infile >> samp1 >> samp2 ) {
        	 fulldat.push_back(samp1);
        	 fulldat2.push_back(samp2);
    	} 
	for (int time = 1; time <= plots; time++) {
		vector< complex<double> > rawdata;
		vector< complex<double> > rawdata2;
		int shift = 64*(time-1);
		for(int j=0; j <N; j++) {
			complex<double> sample(fulldat[j+shift], 0.0);
			complex<double> sample2(fulldat2[j+shift], 0.0);
			rawdata.push_back(sample);
			rawdata2.push_back(sample2);
		}
		FFT fft1 = {powmax, L, false}; 
		fft1.init_data(rawdata);
		FFT fft2 = {powmax, L, false}; 
		fft2.init_data(rawdata2);
		save_psi(time, true, true, fft1.posvec, fft1.f_x, fft1.N);
		fft1.transform();	
		save_psi(time, false, true, fft1.wvec, fft1.F_k, fft1.N);
		save_psi(time, true, false, fft2.posvec, fft2.f_x, fft2.N);
		fft2.transform();	
		save_psi(time, false, false, fft2.wvec, fft2.F_k, fft2.N);
	}
	return 0;
}
