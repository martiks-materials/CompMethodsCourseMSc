// Program for Assignment 1 of Computational Methods:
// This program, when run, will calculate the machine accuracy
// of this system, for both single,  double, and extended 
// precision floating point variables.
#include <iostream>
using namespace std;


int main(){
	// The variables "subject_s/d/e" are reduced until they are
	// equal to the machine accuracy of the system for each
	// floating point variable precision. This is done by
	// adding it to 1 (to give "combo_s/d/e") and compared to
	// unity. The "_s", "_d" and "_e" extensions to variables 
	// refer to "single", "double" and "extended" precision respectively. 
	float unity_s(1.0), subject_s(1.0), combo_s(2.0);
	double unity_d(1.0), subject_d(1.0), combo_d(2.0);
	long double unity_e(1.0), subject_e(1.0), combo_e(2.0);

	// This while loop will continue to divide the "subject" by
	// two and add it to unity, so long as the resulting "combo"
	// is larger that unity. This ensures the loop ends when the 
	// "subject" has been reduced to a value below the  machine 
	// accuracy. This is why, after the loop, the "subject_s/d/e"
	// is multiplied by 2 so that the outputted value is the
	// smallest associated floating point which can be added to
	// unity such that the result is larger than unity.

	// Single precision floating points:
	while(combo_s>unity_s){
		subject_s *= 0.5;
		combo_s = unity_s+subject_s;
	}
	subject_s *= 2;
	
	// Double precision floating points:
	while(combo_d>unity_d){
		subject_d *= 0.5;
		combo_d = unity_d+subject_d;
	}
	subject_d *= 2;

	// Extended precision floating points
	while(combo_e>unity_e){
		subject_e *= 0.5;
		combo_e = unity_e+subject_e;
	}
	subject_e *= 2;
	
	cout << "Floating Point Machine Accuracies." << endl;
	cout << "Single precision (float): " << subject_s << endl;
	cout << "Double precision (double): " << subject_d << endl;
	cout << "Extended precision (long double): " << subject_e << endl;

	return 0;
}
	
		
		
