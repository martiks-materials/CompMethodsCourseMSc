#include <iostream>
#include <fstream>
#include <cmath>
#include "rand.h"
#include "compare.h"
using namespace std;

int main() {
        Ran myran(34);
        int lim(pow(10,5));
        double inc(0.01);
  	double pi = 4*atan(1.0);
        int bins(pow(10,2));
        int hist[bins] = {};
        double ranges[bins] = {};
	double A = 4.0/pi;
        for(int m(0); m<bins; m++){
                ranges[m] = m*inc*pi;
        }

        for(int n(0); n<lim ;n++) {
		bool reject(true);
                while(reject){
			double X = myran.doub()*A;
                	double Y = invcdf(X);
			double Z = myran.doub()*comp(Y);
			if(Z<=pdf(Y)){
				reject=false;
		                for(int l(0); l<bins; l++){
                	        	if(ranges[l]>Y) {
                        	        	hist[l] += 1;
                                		break;
                       			 }
				}
			}
     		 }
        }
        ofstream outfile1("reject.dat");
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
