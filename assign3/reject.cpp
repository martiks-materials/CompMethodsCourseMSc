#include <iostream>
#include <fstream>
#include <cmath>
#include "rand.h"
#include "compare.h"
using namespace std;

int main() {
        Ran myran(42);
        int lim(pow(10,5));
  	double pi = 4*atan(1.0);
        double random[lim] = {};
	double A = 4.0/pi;
        for(int n(0); n<lim ;n++) {
		bool reject(true);
                while(reject){
			double X = myran.doub()*A;
                	double Y = invcdf(X);
			double Z = myran.doub()*comp(Y);
			if(Z<=pdf(Y)){
				reject=false;
				random[n] = Y;
			}
			
     		 }
        }
        ofstream outfile1("reject.dat");
        if ( !outfile1.is_open() ){
                cout << "Error opening file." << endl;
                return 1;
        }
        for(int n(0); n<lim; n++){
                outfile1 << random[n] <<  endl;
        }
        outfile1.close();
        return 0;
}	
