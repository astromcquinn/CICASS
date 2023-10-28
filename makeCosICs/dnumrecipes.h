/* numrecipes.h
 * Header file containing prototypes for functions used from 
 * Numerical Recipes
 */

#ifndef NUM_RECIPES_H
#define NUM_RECIPES_H

// Utility functions
#define NR_END 1 // used by vector
#define FREE_ARG char*

double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void spline(double x[], double y[], int n, double yp1, double ypn,
	    double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x,
	    double *y);
void splie2d(double x1a[], double x2a[], double **ya, int m, int n,
	     double **y2a);
void splin2d(double x1a[], double x2a[], double **ya, double **y2a, int m,
	     int n, double x1, double x2, double *y);
float ran1(long *idum);
extern double gasdev1(long *idum);
double qromb(double (*func)(double), double min, double max, double acc);
double fmin(double A, double B);
double qsimp(double (*func)(double), double a, double b,
	     double accuracy);
double trapzd(double (*func)(double), double lowerEndpoint, 
	      double upperEndpoint, int n, double *s);

void chsone(double *bins, double *ebins,int nbins, int knstrn,
	    double *df, double *chsq, double *prob);
void chstwo(double *bins1, double *bins2, int nbins, int knstrn,
	    double *df, double *chsq, double *prob);
double gammq(double a, double x);
void gser(double *gamser,double a,double x, double *gln);
void gcf(double *gammcf,double a, double x, double *gln);
double gammln(double xx);
double bessj0(double x);
double bessj1(double x);
#endif
