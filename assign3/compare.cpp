#include <cmath>

double comp(double y){
	double pi=4*atan(1.0);
	return (2.0/pi)*sin(y);
}

double pdf(double y){
	double pi=4*atan(1.0);
	return (2.0/pi)*sin(y)*sin(y);
}

double invcdf(double y){
	double pi=4*atan(1.0);
	return acos(1-(pi*y/2));
}

