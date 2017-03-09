#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "fft_gen.h"
#include "func_zoo.h"

using namespace std;
typedef complex<double> dcomp;
typedef double (*func)(double);
typedef vector<double> vec;


int main() {
	// Initiliase functions
	int nm = 10;
	double L = 40.;
	FFT gaussian = {nm, L, false};
	FFT triangle = {nm, L, false};
	FFT square = {nm, L, false};
	FFT doubsquare = {nm, L, false};
	gaussian.init_func(func1);
	triangle.init_func(func2);
	square.init_func(func3a);
	doubsquare.init_func(func3b);
	ofstream outfile1("space_gauss.dat");
	ofstream outfile2("space_tri.dat");
	ofstream outfile3a("space_squ.dat");
	ofstream outfile3b("space_doub.dat");
	for (int index = 0; index < gaussian.N; index++){
		outfile1 << setw(20) << left << gaussian.posvec[index] << setw(20) << left <<  gaussian.f_x[index].real()  << setw(20) << gaussian.f_x[index].imag() << endl;
		outfile2 << setw(20) << left << triangle.posvec[index] << setw(20) << left <<  triangle.f_x[index].real()  << setw(20) << triangle.f_x[index].imag() << endl;
		outfile3a << setw(20) << left << square.posvec[index] << setw(20) << left <<  square.f_x[index].real()  << setw(20) << square.f_x[index].imag() << endl;
		outfile3b << setw(20) << left << doubsquare.posvec[index] << setw(20) << left <<  doubsquare.f_x[index].real()  << setw(20) << doubsquare.f_x[index].imag() << endl;
	}	
	outfile1.close();
	outfile2.close();
	outfile3a.close();
	outfile3b.close();

	gaussian.transform();
	triangle.transform();
	square.transform();
	doubsquare.transform();
	
	ofstream outfilef1("freq_gauss.dat");
	ofstream outfilef2("freq_tri.dat");
	ofstream outfilef3a("freq_squ.dat");
	ofstream outfilef3b("freq_doub.dat");
	for (int index = 0; index < gaussian.N; index++){
		outfilef1 << setw(20) << left << gaussian.wvec[index] << setw(20) << left <<  gaussian.F_k[index].real()  << setw(20) << gaussian.F_k[index].imag() << endl;
		outfilef2 << setw(20) << left << triangle.wvec[index] << setw(20) << left <<  triangle.F_k[index].real()  << setw(20) << triangle.F_k[index].imag() << endl;
		outfilef3a << setw(20) << left << square.wvec[index] << setw(20) << left <<  square.F_k[index].real()  << setw(20) << square.F_k[index].imag() << endl;
		outfilef3b << setw(20) << left << doubsquare.wvec[index] << setw(20) << left <<  doubsquare.F_k[index].real()  << setw(20) << doubsquare.F_k[index].imag() << endl;
	}	
	outfilef1.close();
	outfilef2.close();
	outfilef3a.close();
	outfilef3b.close();
	return 0;
}

