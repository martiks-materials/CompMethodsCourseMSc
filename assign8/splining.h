#pragma once
void triag(double *a, double *b, double *c, double *F, double *u, int n);

void triag_solve(double *xdat, double *ydat, double *u, int n, bool natural, double delta_1=0, double delta_n=0);

double oldspline(double x);

