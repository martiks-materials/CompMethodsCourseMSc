#include <fstream>
#include <iostream>
#include <vector>
#include <valarray>
#include <cmath>
#include "brent.h"
#include "rk4.h"
#include "specifics.h"


int main() {
	double epsilon(1E-8), delta(1E-8);
        int maxi(10000);
	vec yinit = {-10, 0};
	BrentMethod Hitman = {solution, yinit, epsilon, maxi, delta};
	double result = Hitman.operate();
	std::cout << "Best guess for y2(0) = " << result << " in " << Hitman.count << " iterations" << std::endl;	
	int n_d(2);
	double dx(0.1), xmax(10), xmin(0), eps(1E-6), safety(0.98);
	std::valarray<double> y = {1.5, result};
	Rungekutta4 Rk4(y, &f, dx, eps, xmax, xmin, n_d, false, safety);	
	while(Rk4.x_now<xmax){
		Rk4.iterate();
	}
	
	
	// Data for (x, y1), (x, y2), and (y1, y2) is stored for plotting.

	std::fstream outfile1("runge1.dat");
	std::ofstream outfile2("runge2.dat");
	std::ofstream outfile3("runge3.dat");
	if(((!outfile1.is_open())||(!(outfile2.is_open()))||(!outfile3.is_open()))){
		std::cout << "Error opening file..." << std::endl;
		return 1;
	}
	

	std::cout << "No. Steps = " << Rk4.numvals-1 << std::endl;
        std::cout << "Fnc. Evaluations = " << Rk4.evals << std::endl;
        std::cout << "Repeated steps = " << Rk4.repeats << std::endl;

	for(int i(0); i<Rk4.numvals ; i++){
		outfile1 << Rk4.xlist.at(i) << " " << Rk4.yn[0].at(i) << std::endl;
		outfile2 << Rk4.xlist.at(i) << " " << Rk4.yn[1].at(i) << std::endl;
		outfile3 << Rk4.yn[0].at(i) << " " << Rk4.yn[1].at(i) << std::endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	return 0;
}
		
	













