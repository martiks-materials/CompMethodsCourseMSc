#include <iostream>
#include <fstream> 
#include <cmath>
#include "rand.h"
using namespace std;

int main() {
	Ran myran(34);
        int lim(pow(10,5));
        double inc(0.01);
	double pi(3.14159265359);
	int bins(pow(10,2));
        int hist[bins] = {};
        double ranges[bins] = {};
        for(int m(0); m<bins; m++){
                ranges[m] = m*inc*pi;
        }

        for(int n(0); n<lim ;n++) {
                double X = myran.doub();
		double Y = acos(1-(2*X));
                for(int l(0); l<bins; l++){
                        if(ranges[l]>Y) {
                                hist[l] += 1;
                                break;
                        }
                }
        }
        ofstream outfile1("transform.dat");
        if ( !outfile1.is_open() ){
     	        cout << "Error opening file." << endl;
                return 1;
        }
        for(int n(0); n<bins; n++){
                outfile1 << ranges[n] << " " << hist[n] << endl;
        }
        outfile1.close();
	return 0;
}

	
