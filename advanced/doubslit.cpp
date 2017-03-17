
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

#include "fft_gen.h"

typedef complex<double> dcomp;
typedef vector< complex<double> > vecomp;
typedef double (*func)(double);
typedef vector<double> vec;


double slit(double x, double w) {
        // double slit
	double a = 0.75;
        if(((x>-a-w)&&(x<-w+a))||((x>w-a)&&(x<w+a))){
                return 2;
        }
        else {
                return 0;
        }
}

void save_psi(int plot_number, bool space,  vector<double> xdat, vector< complex<double> > psidat, int Num) {
	ostringstream os;
	if(space){
		os << "gratspace/space_" << plot_number << ".data";
	}
	else {
		os << "gratmom/mom_" << plot_number << ".data";
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
		file << xdat[i] << "\t" << norm(psidat[i]) << "\t" << psidat[i].real() << "\t" << psidat[i].imag() << "\n";
	}
	file.close();
}

int main() {
	int plots = 999;
	int powmax = 10;
	int N = int(pow(2, powmax));
	double L = 80;
	double h = L/double(N);
	double width = 0.75;
	int count = 0;
	for (int plot = 1; plot <= plots; plot++) {
		vector< complex<double> > rawdata;
		width += plot*0.0001;
		count++;
		for(int j=-N/2; j < N/2; j++) {
			rawdata.push_back(slit(j*h, width));
		}
		FFT fft = {powmax, L, false}; 
		fft.init_data(rawdata);
		save_psi(plot, true, fft.posvec, fft.f_x, fft.N);
		fft.transform();	
		save_psi(plot, false, fft.wvec, fft.F_k, fft.N);
	}
	return 0;
}
