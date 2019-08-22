// poly34.h : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com

#ifndef poly34
#define poly34

int   SolveP3(double *x,double a,double b,double c);			// solve cubic equation x^3 + a*x^2 + b*x + c = 0
int   SolveP4(double *x,double a,double b,double c,double d);	// solve equation x^4 + a*x^3 + b*x^2 + c*x + d = 0 by Dekart-Euler method

int   SolveP4Bi(double *x, double b, double d);				// solve equation x^4 + b*x^2 + d = 0
int   SolveP4De(double *x, double b, double c, double d);	// solve equation x^4 + b*x^2 + c*x + d = 0
void  CSqrt( double x, double y, double &a, double &b);		// returns as a+i*s,  sqrt(x+i*y)
int NewtonStepP4(double *x, double a,double b,double c,double d);// one Newton step for x^4 + a*x^3 + b*x^2 + c*x + d

double SolveP5_1(double a,double b,double c,double d,double e);	// return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0

int SolveP2(double *x, double a, double b);	// solve equation x^2 + a*x + b 

#endif