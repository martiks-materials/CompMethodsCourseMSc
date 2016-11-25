#include <iostream>
#include <cmath>
using namespace std;

class Integrate{
private:
	double (*fnc)(double x);
	double a;
	double b;   //should these all be on the same line???
	double t;	
	int n;
	int evals;
	double h;
	double eps;
	double I2;
	double I1;
public:
	Integrate(double (*func)(double x), double low, double up, double epsilon){
		a = low;
		b = up;
		h = 0.5*(b-a); // times by two since we perform 2 iterations in initialization
		n = 1;
		eps = epsilon;
		fnc = func;	 //Correct way to pass on functions????
		I1 = h*((*fnc)(a)+(*fnc)(b));
		I2 = ((0.5*I1)+((*fnc)(0.5*(a+b))*h));
		evals = 3;
	}

	void iterate() {
		I1 = I2;
		I2 = 0.5*I1;
		n *= 2;
		for(int i(0); i<n; i++){
			double x_new = a + (i+0.5)*h;
			I2 += 0.5*h*(*fnc)(x_new);
			evals++;
		}
		h *= 0.5;
	}

	double trap(){
		while(abs(I2-I1)>eps*abs(I1)){
			Integrate::iterate();
		}
		return I2;
	}
	
	double simp(){
		double s1 = (4*I2/3.)-(I1/3.);
		Integrate::iterate();
		double s2 = (4*I2/3.)-(I1/3.);
		while(abs(s2-s1)>eps*abs(s1)){
			s1 = s2;
			Integrate::iterate();
			s2 = (4*I2/3.)-(I1/3.);
		}
		return s2;
	}
	
	int get_evals() { return evals; }
};
