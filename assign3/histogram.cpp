#include <iostream>
#include <fstream>
#include <cmath>
#include "rand.h"
using namespace std;



int main(){
	Ran myran(34);
	int lim(pow(10,5));
	double inc(0.01);
	int bins(pow(10,2));
	int hist[bins]={};
	double ranges[bins]={};
	for(int j(0);j<bins;j++) {	
		ranges[j] = j*inc;
	}
	for(int n(0); n<lim ;n++) {
		double X = myran.doub();
		for(int j(0); j<bins; j++) {
			if(ranges[j]>X){
				hist[j]+= 1;
				break;
			}
		}
	}
	ofstream outfile1("hist.dat");
	if ( ! outfile1.is_open() ) {
    	    cout << "Error opening file." << endl;
        	return 1;
    	}	
	for(int n(0); n<bins; n++){
		outfile1 << ranges[n] << " " << hist[n] << endl;
	}
	outfile1.close();
	return 0;
}

