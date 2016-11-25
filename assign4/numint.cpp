#include <iostream>
#include <cmath>
#include "integrate.h" // do i need a header for each function?
// do i need a file for each function?
using namespace std;

double func1(double x) {
        return (2/sqrt(atan(1.0)*4))*exp(-x*x);
}

int main() {
        Integrate Result1(&func1, 0, 2, 1e-6);
	Integrate Result2(&func1, 0, 2, 1e-6);
	cout << "Trapezium rule: " << Result1.trap();
	cout <<  " using " << Result1.get_evals() << " evaluations." << endl;
	cout << "Simpsons's rule: " << Result2.simp();
	cout << " using " << Result2.get_evals() << " evaluations." << endl;
	return 0;
}
