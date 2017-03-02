#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
using namespace std;
typedef complex<double> dcomp;

double func1(double x, double a){
	// Gaussian
	return exp(-x*x/a);
}

double func2(double x, double a) {
	// Triangle
	if((x>-a)&&(x<a)){
		return a-abs(x);;
	}
	else {
		return 0;
	}
}
double func3(double x, double a) {
	// square
	if((x>-a)&&(x<a)){
		return 1/a;
	}
	else {
		return 0;
	}
}

double func4(double x, double a) {
	// abs exponential
	return exp(-a*abs(x));
}

int binaryflip( int inp, int maxpow){
	int result = 0;
	for(int n=maxpow-1; n>=0; n--){
		int n2 = int(pow(2, n));
		if(inp-n2>=0){
			inp -= n2;
			result += int(pow(10, maxpow-1-n));
		}
	}
	return result;
}


// Trial is out now and then turn it into a class.
int main(){
	// Number of real-space samples
	dcomp I(0.0, 1.0);
	int nmax = 12;
	int N = int(pow(2, nmax));
	cout << N << endl;
	// Spatial step
	double L = 70.0;
	double A = 0.5;
	double h = L/double(N);
	vector<double> f_x; 
	vector<long unsigned int> binum;
	// Sample the function
	ofstream outfile2("space.dat");
	for(int j=-N/2; j<N/2+1; j++){
		double ffs = func1(double(j)*h, A);
		outfile2 << h*j << right << setw(20) << ffs << endl;
		f_x.push_back(ffs);
	}
	outfile2.close();
	for(int n(0); n<N; n++){
                binum.push_back(binaryflip(n, nmax));
        }
	cout << endl;
	// write the sample indices in terms of reversed bits of each number	
	int i2, j1, k1;
	i2 = N >> 1;
	j1 = 0;
	for(int i=0;i<N-1;i++) {
		if(i<j1) {
			double tx = f_x[i];
			int bt = binum[i]; 
			f_x[i] = f_x[j1];
			f_x[j1] = tx;
			int nt = binum[i];
			binum[i] = binum[j1];
			binum[j1] = bt;
		}
		k1 = i2;
		while(k1 <= j1) {
			j1 -= k1;
			k1 >>= 1;
		}
		j1 += k1;
	}
	/*
	for(int n(0); n<N; n++){
		binum.push_back(binaryflip(n, nmax));
	}

	// Order the sample f_x with respect to the bit flipped numbers
	for(int start(0); start<N; start++){
		int smallest = start;
		for(int current = start+1; current<N; current++){
			if(binum[current]<binum[smallest]){
				smallest = current;
			}
		}
		// Swap binary numbers but also swap function values
		swap(binum[start], binum[smallest]);
		swap(f_x[start], f_x[smallest]);
	}

	for (int index = 0; index < N; index++){
		cout << f_x[index] << ' ';
	}	
	cout << endl;
	for (int index = 0; index < N; index++){
		cout << binum[index] << ' ';
	}
	cout << endl;
	*/



	// Create vector of fourier transforms.
	vector<double> F_kr, F_ki;
	double dk = 1./(h*double(N));
	double PI = 4*atan(1);
	// For each k inbetween the positive and negative Nyquist Frequencies, do transform
	for(int k=-N/2+1; k<N/2+1;k++){
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
			dcomp Wk, theta(2*PI*k*N_new/double(N), 0);
			Wk = exp(I*theta);
			
			for(int j=0; j<N_new; j++) {
				Fouri_temp.push_back( Fouri[2*j] +Wk*Fouri[2*j + 1] );
			}
			Fouri.resize(N_new);
			for(int j(0); j<N_new; j++){
				Fouri[j] = Fouri_temp[j];
			};
		}
		//cout << Fouri.size() << endl;
		if(k%2==0) {
			F_kr.push_back(h*Fouri[0].real());
			F_ki.push_back(h*Fouri[0].imag());
		}
	
		else {
			F_kr.push_back(-h*Fouri[0].real());
                        F_ki.push_back(-h*Fouri[0].imag());
		}
	}
	ofstream outfile("freq.dat");
	for (int index = 0; index < N; index++){
		outfile << setw(20) << left << (index+1-double(N)/2)*dk << setw(20) << left <<  F_kr[index]  << setw(20) << F_ki[index] << endl;
	}	
	outfile.close();
	return 0;
	
	
}










