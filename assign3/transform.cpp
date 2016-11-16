#include <iostream>
#include <fstream> 
#include <cmath>
#include <ctime>
#include "rand.h"
using namespace std;

int main() {
	time_t start = clock();
	Ran myran(42);
        int lim(pow(10,5));
        double inc(0.01);
	double pi(4*atan(1.0));
        double random[lim] = {};
        for(int n(0); n<lim ;n++) {
                double X = myran.doub();
		double Y = acos(1-(2*X));
		random[n] = Y;
        }
	time_t end = clock();
	cout << "Time take to run: " << (end-start)/(double)CLOCKS_PER_SEC << "s" << endl;
        ofstream outfile1("transform.dat");
        if ( !outfile1.is_open() ){
     	        cout << "Error opening file." << endl;
                return 1;
        }
        for(int n(0); n<lim; n++){
                outfile1 << random[n] << endl;
        }
        outfile1.close();
	return 0;
}

	
