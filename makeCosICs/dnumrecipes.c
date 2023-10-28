#define JMAX 30   // maximum number of steps for standard
#define JMAXM 14  // maximum number of steps for open
#define MAXJ 5     // order of interpolation
#define EPS 1.0e-4

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "dnumrecipes.h"


double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i;
  long nrow(nrh-nrl+1);
  long ncol=(nch-ncl+1);
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m)
    {  fprintf(stderr, "allocation failure 2 in dmatrix()"); exit(-9);}
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl])
    {  fprintf(stderr, "allocation failure 2 in dmatrix()"); exit(-9);}
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i=nrl+1;i<=nrh;i++) 
    m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch) 
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

double *dvector(long nl, long nh)
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

/* Constructs spline table for 1-D interpolation.  x is dependent variable
 * vector, y is independent variable vector, and table is returned
 * as y2.  n is the array size.  yp1, yp2 are set to 1.0e30 for natural
 * splines. */
void spline(double x[], double y[], int n, double yp1, double ypn,
	    double y2[])
{
  int i,k;
  double p,qn,sig,un,*u;

  u = dvector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1] = u[1] = 0.0;
  else {
    y2[1] = -0.5;
    u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k] = y2[k]*y2[k+1]+u[k];
  free_dvector(u,1,n-1);
  return;
}

/* Spline evaluation function.  xa is the dependent variable vector
 * and ya is the independent variable vector.  y2a is the spline table
 * from a previous call to spline.  n is the array size.  x is the point
 * at which the interpolation is to be performed.  Result is returned
 * in y. NOTE xa must be in increasing order! */
void splint(double xa[], double ya[], double y2a[], int n, double x,
	    double *y)
{
  int klo,khi,k;
  double h,b,a;

  klo = 1;
  khi = n;
  while (khi-klo > 1) {
    k = (khi + klo) >> 1;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }
  h = xa[khi] - xa[klo];
  if (h == 0.0)
    {
      fprintf(stderr, "Bad xa input into routine splint");
      exit(-1);
    }
  a = (xa[khi] - x)/h;
  b = (x - xa[klo])/h;
  *y = a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] + 
				(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return;
}

/*********************************************************************
 ******************* Random Number Generators ************************
 ********************************************************************/

float ran1(long *idum)
{
  const int IA = 16807,IM =2147483647,IQ=127773,IR=2836;
  const float AM=1.0/IM,REPS=1.2e-7,RNMX=1.0-REPS;
  const int NTAB=32,NDIV= (1+(IM-1)/NTAB);

  int j;
  long k;
  static long iy= (0);
  static long iv[32];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) 
      *idum=1;
    else 
      *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k = (*idum)/IQ;
      *idum = IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0)
	*idum += IM;
      if (j < NTAB)
	iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0)
    *idum += IM;
  j = iy/NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX)
    return RNMX;
  else
    return temp;
}

double gasdev1(long *idum)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if (*idum < 0)
    iset = 0;
  if (iset == 0) {
    do {
      v1 = 2.0*ran1(idum) - 1.0;
      v2 = 2.0*ran1(idum) - 1.0;
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}

/*********************************************************************
 ********************** Integration Routines *************************
 ********************************************************************/
// Integral from min to max using Romberg integration.  Converges provided
// f is continuous on the interval.  func is double precision.  
// accuracy is desired relative accuracy in the integral.
// Taken ultimately from CMBFAST (Seljak & Zaldarriaga 1996).
double qromb(double (*func)(double), double min, double max, double acc)
{
  double h,gmax,error,g0 = 0.0,g1,fourj;
  double g[MAXJ+1];
  int jmaxx,nint,i,j, k;
  
  h = 0.5*(max - min);
  gmax = h*((*func)(min) + (*func)(max));
  g[1] = gmax;
  nint = 1;
  error = 1.0e20;
  for (i=1;i <= JMAX;i++) {
    if (i > 5 && fabs(error) < acc)
      break;
    
    // Calculate next trapezoidal rule approx for integral
    g0 = 0.0;
    for (k=1;k<=nint;k++) {
      g0 = g0 + (*func)(min+ (double) (k+k-1)*h);
    }

    g0 = 0.5*g[1]+h*g0;
    h = 0.5*h;
    nint *= 2;
    jmaxx = fmin(i,MAXJ);
    fourj = 1.0;

    // Richardson extrapolation
    for (j=1; j<=jmaxx; j++) {
      fourj *= 4.0;
      g1 = g0 + (g0-g[j])/(fourj-1.0);
      g[j]=g0;
      g0=g1;
    }
    
    if (fabs(g0) > acc) {
      error = 1.0 - gmax/g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmaxx+1]=g0;
  }
  if (i > JMAX && fabs(error) > acc)
  {
      fprintf(stderr,"rombint failed to converge");
      exit(-3);
  }
  return g0;
}

double qsimp(double (*func)(double), double a, double b,
	     double accuracy)
{
  int j;
  double s,st,ost,os;
  double sCurrent = 0.0;

  ost = os = -1.0e30;
  for (j=1;j<=JMAX;j++) {
    st = trapzd(func,a,b,j,&sCurrent);
    s = (4.0*st - ost)/3.0;
    if (j > 5) {
      if (fabs(s-os) < accuracy*fabs(os) || (s==0.0 && os == 0.0))
	return s;
    }
    os = s;
    ost = st;
  }
  fprintf(stderr, "Too many steps in routine qsimp");
  return 0.0;
}

// Performs nth iteration of trapezoid rule; must be called in sequential 
// order of n's.  (*func) is the function to be integrated.
double trapzd(double (*func)(double), double lowerEndpoint, 
	     double upperEndpoint, int n, double *s)
{
  double x,tnm, sum, spacing;
  int numberOfEvals, j;

  if (n==1) {
    return (*s = 0.5*(upperEndpoint - lowerEndpoint)*
	    (((*func)(lowerEndpoint)) + ((*func)(upperEndpoint))));
  } else {
    for (numberOfEvals=1,j=1;j<(n-1);j++) 
      numberOfEvals <<= 1;
    tnm = numberOfEvals;
    spacing = (upperEndpoint - lowerEndpoint)/tnm;
    x = lowerEndpoint + 0.5*spacing;
    for (sum=0.0,j=1;j<=numberOfEvals;j++,x+=spacing)
      sum += (*func)(x);
    *s = 0.5*((*s) + (upperEndpoint - lowerEndpoint)*
			sum/tnm);
    return *s;
  }
}

double fmin(double A, double B)
{
  return ((A < B) ? A : B);

}

void chsone(double *bins, double *ebins,int nbins, int knstrn,
	    double *df, double *chsq, double *prob)
{
	int j;
	double temp;
	double gammq(double a, double x);
	void nrerror();

	*df=nbins-1-knstrn;
	*chsq=0.0;
	for (j=1;j<=nbins;j++) {
		if (ebins[j] <= 0.0){
		    fprintf(stderr, "Bad expected number in CHSONE");
		    exit(-1);
		}
		temp=bins[j]-ebins[j];
		*chsq += temp*temp/ebins[j];
	}
	*prob=gammq(0.5*(*df),0.5*(*chsq));
}


void chstwo(double *bins1, double *bins2, int nbins, int knstrn,
	    double *df, double *chsq, double *prob)
{
	int j;
        double temp;
	//double gammq();

	*df=nbins-1-knstrn;
	*chsq=0.0;
	for (j=1;j<=nbins;j++)
		if (bins1[j] == 0.0 && bins2[j] == 0.0)
			*df -= 1.0;
		else {
			temp=bins1[j]-bins2[j];
			*chsq += temp*temp/(bins1[j]+bins2[j]);
		}
	*prob=gammq(0.5*(*df),0.5*(*chsq));
}

#define ITMAX 100

double gammq(double a, double x)
{
	double gamser,gammcf,gln;
	//void gcf(),gser();

	if (x < 0.0 || a <= 0.0) {
	    fprintf(stderr, "Invalid arguments in routine GAMMQ");
	    exit(-1);
	}if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

void gser(double *gamser,double a,double x, double *gln)
{
	int n;
	double sum,del,ap;
	//double gammln();

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) {
		    fprintf(stderr, "x less than 0 in routine GSER");
		    exit(-1);
		}
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		fprintf(stderr,"a too large, ITMAX too small in routine GSER");
		exit(-1);
	}
}


void gcf(double *gammcf,double a, double x, double *gln)
{
	int n;
	double gold=0.0,g,fac=1.0,b1=1.0;
	double b0=0.0,anf,ana,an,a1,a0=1.0;
	//double gammln();

	*gln=gammln(a);
	a1=x;
	for (n=1;n<=ITMAX;n++) {
		an=(float) n;
		ana=an-a;
		a0=(a1+a0*ana)*fac;
		b0=(b1+b0*ana)*fac;
		anf=an*fac;
		a1=x*a0+anf*a1;
		b1=x*b0+anf*b1;
		if (a1) {
			fac=1.0/a1;
			g=b1*fac;
			if (fabs((g-gold)/g) < EPS) {
				*gammcf=exp(-x+a*log(x)-(*gln))*g;
				return;
			}
			gold=g;
		}
	}
	fprintf(stderr, "a too large, ITMAX too small in routine GCF");
	exit(-1);
}

#undef ITMAX

double gammln(double xx)
{
	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}


double bessj0(double x)
{
  double ax,z;
  double xx,y,ans,ans1,ans2;

  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=(57568490574.0
	  +y*(-13362590354.0
	      +y*(651619640.7
		  +y*(-11214424.18+y*(77392.33017+y*(-184.9052456))))));
    ans2=(57568490411.0
	  +y*(1029532985.0
	      +y*(9494680.718
		  +y*(59272.64853+y*(267.8532712+y*1.0)))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-0.785398164;
    ans1=1.0+y*(-0.1098628627e-2
		+y*(0.2734510407e-4
		    +y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = (-0.1562499995e-1
	    +y*(0.1430488765e-3
		+y*(-0.6911147651e-5
		    +y*(0.7621095161e-6-y*0.934935152e-7))));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}


double bessj1(double x)
{
  double ax,z;
  double xx,y,ans,ans1,ans2;

  if ((ax=fabs(x)) < 8.0) {
    y = x*x;
    ans1 = x*(72362614232.0 + 
	      y*(-7895059235.0 + 
		 y*(242396853.1 + 
		    y*(-2972611.439 + 
		       y*(15704.48260 + 
			  y*(-30.16036606))))));
    ans2 = (144725228442.0 + 
	    y*(2300535178.0 + 
	       y*(18583304.74 + 
		  y*(99447.43394 + 
		     y*(376.9991397 + y*1.0)))));
    ans = ans1/ans2;
  } else {
    z = 8.0/ax;
    y = z*z;
    xx = ax - 2.356194491;
    ans1 = (1.0 + 
	    y*(0.183105e-2 + 
	       y*(-0.3516396496e-4 + 
		  y*(0.2457520174e-5 +
		     y*(-0.240337019e-6)))));
    ans2 = (0.04687499995 + 
	    y*(-0.2002690873e-3 + 
	       y*(0.8449199096e-5 +
		  y*(-0.88228987e-6 + 
		     y*0.105787412e-6))));
    ans = sqrt(0.636619772/ax)*(cos(xx)*ans1 - z*sin(xx)*ans2);
    if (x < 0.0)
      ans = -ans;
  }
  return ans;
}


/* Function to create 2D spline table.  x1a,x2a are the dependent variable
 * vectors.  ya is the matrix of function evaluations, of size m x n.
 * y2a is where the spline table is put. */
void splie2d(double x1a[], double x2a[], double **ya, int m, int n,
	    double **y2a)
{
  int j;

  for (j=1;j<=m;j++)
    spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
  return;
}

/* Function to perform spline interpolation in 2D.  x1a, x2a are 
 * dependent variable vectors of length m and n, respectively.  ya is the
 * table of function evaluations.  y2a is the spline table, from a 
 * previous call to splie2.  x1, x2 is the location of the interpolation
 * point, and the result is returned in y. */
void splin2d(double x1a[], double x2a[], double **ya, double **y2a, int m,
	    int n, double x1, double x2, double *y)
{
  int j;
  double *ytmp,*yytmp;

  ytmp = dvector(1,m);
  yytmp = dvector(1,m);
  for (j=1;j<=m;j++)
    splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
  splint(x1a,yytmp,ytmp,m,x1,y);
  free_dvector(yytmp,1,m);
  free_dvector(ytmp,1,m);
}


