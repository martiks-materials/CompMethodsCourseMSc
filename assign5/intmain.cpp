#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <valarray>
#include "specifics.h"
#include "rk4.h"
using namespace std;

int main(){
	valarray<double> y = {1.5, 1.5};
	double dx(0.1), xmax(10), eps(1E-6);
	Rungekutta4 Rk4(y, &f, dx, eps, xmax, 2);	
	while(Rk4.x_now<xmax){
		Rk4.iterate();
	}
	
	ofstream outfile1("runge1.dat");
	ofstream outfile2("runge2.dat");
	ofstream outfile3("runge3.dat");
	if(((!outfile1.is_open())||(!(outfile2.is_open()))||(!outfile3.is_open()))){
		cout << "Error opening file..." << endl;
		return 1;
	}
	
	for(int i(0); i<Rk4.counter ; i++){
		outfile1 << Rk4.xlist.at(i) << " " << Rk4.yn[0].at(i) << endl;
		outfile2 << Rk4.xlist.at(i) << " " << Rk4.yn[1].at(i) << endl;
		outfile3 << Rk4.yn[0].at(i) << " " << Rk4.yn[1].at(i) << endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	return 0;
}
			
	
