#include <iostream>
using namespace std;

double linear(double x, double x_dat[], double y_dat[], int length){
        int j(0);

        for(int i(0); i<length; i++) {
                if(x<x_dat[i]) {
                        j--;
                        break;
                }
                j++;

        }

        double A = (x_dat[j+1] - x)/(x_dat[j+1] - x_dat[j]);
        double B = (x - x_dat[j])/(x_dat[j+1] - x_dat[j]);

        return (A*y_dat[j]) + (B*y_dat[j+1]);
}
