#include <iostream>
#include "func_zoo.h"
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>

double func1(double x){
        // Gaussian function.
	double a = 1;
        return exp(-x*x/a);
}

double func2(double x) {
        // Triangle function.
	double a = 1;
        if((x>-a)&&(x<a)){
                return a-abs(x);
        }
        else {
                return 0;
        }
}

double func3a(double x) {
        // Square function.
	double a = 1;
        if((x>-a)&&(x<a)){
                return 1/a;
        }
        else {
                return 0;
        }
}


double func3b(double x) {
        // Double square function.
	double a = 0.5;
	double w = 1;
        if(((x>-w-a)&&(x<-w+a))||((x>w-a)&&(x<w+a))){
                return 1/a;
        }
        else {
                return 0;
        }
}

