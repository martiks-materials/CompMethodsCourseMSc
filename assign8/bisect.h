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
	int step, maxstep, numspace;
public:
	vec fvals, xvals1, xvals2;
	Bisector(func fun, vec xinit, double eps, int maxsteps, int N=15){
		f = fun;
		x1 = xinit[0];
		x2 = xinit[1];
		f1 = (*f)(x1);
		f2 = (*f)(x2);
		f3 = 0;
		tol = eps;
		step = 0;
		maxstep = maxsteps;
		numspace = N;
	}

	double operate() {
		bool chose1;
		while(step < maxstep){
			if(sign(f1*f2)>=0) {
				//bool found=false;
				//double h = 1/double(numspace);
				//for(int i(0); i<numspace; i++){
				//	double trialx = x1+(x2-x1)*h*i;
				//	double trialf = (*f)(trialx);
				//	if(sign(trialf*f1)<0){
				//		f2 = trialf;
				//		x2 = trialx;
				//		found = true;
				//		break;
				//	}
				//}
				// If same sign, then not guaranteed to work, so a better
				// initial guess for the root bracketing must be given.
				//if(found){
				//	continue;
				//}
				//else {
					cout << "Failed to bracket." << endl;
					break;
				//}
			}
			// cout << "x1 = " << x1 << ", x2 = " << x2 << endl;
			double x3 = 0.5*(x1+x2);
			f3 = (*f)(x3);
			double dx;
			// cout << "f = " << f3 << endl;
			if(sign(f3*f2)<0) {
				// x3 ad x2 bracket the root
				x1 = x3;
				f1 = f3;
				dx = 0.5*abs(x1-x2);
				chose1 = false;
			}
			else if(sign(f1*f3)<0){
				// x3 and x1 bracket the root
				x2 = x3;
				f2 = f3;
				dx = 0.5*abs(x1-x2);
				chose1=true;
			}
			step++;
			if((abs(f3)<1E-16)||dx<tol){
				cout << "Converged to: " << x3 << ", f = " << f3 << " in " << step << " iterations." << endl;
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
