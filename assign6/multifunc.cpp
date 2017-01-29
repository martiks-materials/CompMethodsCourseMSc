// Collection of functions required for the Monte Carlo
// integration of the error function at x=2. Includes two 
// two different sets of probability distribution and 
// corresponding inverse cummulative distribution functions
// for the importance sampling.
#include <cmath>
#include <vector>
using namespace std;

double erfi(vector<double> x){
	// Integrand of the error function.
        return (2./sqrt(M_PI))*exp(-x[0]*x[0]);
}

double pdf(vector<double> y){
	// Probability Density function no. 1 for the 
	// implementation of the importance sampling.
	// for MC integrating the error function.
        return 0.98-(0.48*y[0]);
}

double invcdf(double x){
       	// Inverse cummulative distribution function no. 1
	// for the implementation of the importance sampling.
	// for MC integrating the error function.
        double X = (0.98-sqrt((0.98*0.98)-(0.96*x)))/0.48;
        return X;
}

double pdf2(vector<double> y){
       	// Probability Density function no. 2 for the 
	// implementation of the importance sampling.
	// for MC integrating the error function.
         return 1.5*exp(-1.5*y[0])/(1-pow(M_E, -3));
}

double invcdf2(double x){
       	// Inverse cummulative distribution function no. 2
	// for the implementation of the importance sampling.
	// for MC integrating the error function.
        return -(2./3)*log(1-((1-pow(M_E,-3))*x));
}

