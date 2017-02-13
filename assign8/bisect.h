#pragma once

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

typedef double (*func)(double x);
typedef vector<double> vec;

class Bisector{
private:
	func f;
	double x1, x2, f1, f2, f3, tol;
	int maxstep, numspace;
public:
	int count;
	vec fvals, xvals1, xvals2;
	Bisector(func fun, vec xinit, double eps, int maxsteps, int N=15){
		f = fun;
		x1 = xinit[0];
		x2 = xinit[1];
		f1 = (*f)(x1);
		f2 = (*f)(x2);
		f3 = 0;
		tol = eps;
		count = 0;
		maxstep = maxsteps;
		numspace = N;
	}

	double operate() {
		bool chose1;
		while(count < maxstep){
			if(sign(f1*f2)>=0) {
				cout << "Failed to bracket." << endl;
				break;
			}
			double x3 = 0.5*(x1+x2);
			f3 = (*f)(x3);
			double dx;
			if(f3*f2<0) {
				x1 = x3;
				f1 = f3;
				dx = 0.5*abs(x1-x2);
				chose1 = false;
			}
			else if(f1*f3<0){
				// x3 and x1 bracket the root
				x2 = x3;
				f2 = f3;
				dx = 0.5*abs(x1-x2);
				chose1=true;
			}
			count++;
			if((abs(f3)<1E-16)||dx<tol){
				return x3;
				break;
			}
		}
		double result = chose1? x2:x1;
		return result;
		
	}
		
	int sign(double x) {
		return ((x>0)-(x<0));
	}

};
