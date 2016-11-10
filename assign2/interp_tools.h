#ifndef INTERP_TOOLS_H
#define INTERP_TOOLS_H
double linear(double x, double x_dat[], double y_dat[], int length);

void triag(double *a, double *b, double *c, double *F, double *u, int n);

void triag_solve(double *xdat, double *ydat, double *u, double n, double natural);

double spline(double x, double *xdat, double *ydat, double *yderiv, int length);

#endif 
