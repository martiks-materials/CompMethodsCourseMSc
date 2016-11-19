#include <cmath>

double comp(double y){
	// Returns a function evaluation of the comparison function needed
	// in the rejection method, where pi is obtained through 4*arctan(1.0).

	double pi=4*atan(1.0);
	return (2.0/pi)*sin(y);
}

double pdf(double y){
	// Returns a function evaluation of the desired pdf for the rejection
	// method exercise. Pi is obtained through 4*arctan(1.0).
	
	double pi=4*atan(1.0);
	return (2.0/pi)*sin(y)*sin(y);
}

double invcdf(double y){
	// Returns a function evaluation of the inverse of the cummulative
	// distribution function of the comparison function. Pi is obtained 
	// through 4*arctan(1.0).

	double pi=4*atan(1.0);
	return acos(1-(pi*y/2));
}

