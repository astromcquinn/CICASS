
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern double BOXSIZE, VSTREAM, ZSTREAM, ZINIT, AINIT, TAVG, NEAVG, hconst;
extern int SIZE;
extern double DIM_ARRAY, NUMK_PAR_ARRAY;

#ifdef GSL_INTERP

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/**********************************************************
Program to return growth factors
flag == 0 initializes
********************************************************/
void interpPk(double k, double costh, double G[], gsl_interp_accel **acc[], int flag)
{
  int i, j, n, s;
  static int numk_par, numk_perp, dim;
  static double *kparvec, *kperpvec, **Garr[20]; //the 20s here are just bad programming style as (# of fields used) << 20
  static gsl_spline **Garr2[20];
  void splie2dgsl(double x1a[], double x2a[], double **ya, gsl_spline **spline, int m, int n);
  void splin2dgsl(double x1a[], double x2a[], gsl_spline **spline, gsl_interp_accel **acc, int m,
		int n, double x1, double x2, double *y);

  if(flag == 0)//read in theta grid
    {    
      FILE *inf;
      char filename[200];

      s = SIZE;
#ifdef USE_OTHER_GRID_SIZE
      s= USE_OTHER_GRID_SIZE;
#endif

      sprintf(filename, "initSimCartZI%.1lf_Vbc%.1lf_%d_%.1lf.dat", ZINIT, VSTREAM, s, BOXSIZE);
      fprintf(stdout, "Using %s!!!\n", filename);

      if((inf = fopen(filename, "r")) == NULL)
	{
	  fprintf(stderr, "Could not open power spectrum file %s\n", filename);
	  exit(-5);
	}
      fprintf(stdout, "#Opened power spectrum file %s\n", filename);

      fread(&dim, sizeof(int), 1, inf);  
      fread(&numk_par, sizeof(int), 1, inf);
      fread(&numk_perp, sizeof(int), 1, inf);
	    

      DIM_ARRAY = dim;
      NUMK_PAR_ARRAY = numk_par;
      //fprintf(stderr, "#dim %d numkpar %d numkperp %d\n", dim, numk_par, numk_perp);

      kparvec = (double *) malloc(sizeof(double)*numk_par); //dvector(1, numk_par);
      kperpvec =  (double *) malloc(sizeof(double)*numk_perp);

      fread(kparvec, sizeof(double), numk_par, inf);

      fread(kperpvec, sizeof(double), numk_perp, inf);



      fread(&TAVG, sizeof(double), 1, inf);
      fread(&NEAVG, sizeof(double), 1, inf);
	    
      fprintf(stdout, "#dim %d numkpar %d numkperp %d TAVG %g NEAVG %g\n", dim, numk_par, 
	      numk_perp, TAVG, NEAVG);
      


      for(n=0; n< dim; n++)    
	{
	  Garr[n] = (double **) malloc(sizeof(double *)*numk_par);
	  Garr2[n] = (gsl_spline **) malloc(sizeof(gsl_spline *)*numk_par);
	  // acc[n] = (gsl_interp_accel **) malloc(sizeof(gsl_interp_accel *)*numk_par);
	  for(i=0; i< numk_par;i++)
	    {
	      Garr[n][i] = (double *) malloc(sizeof(double)*numk_perp);
	      Garr2[n][i] = gsl_spline_alloc(gsl_interp_cspline, numk_perp);
	      //     acc[n][i]  = gsl_interp_accel_alloc ();
	    }
	}


      for(i=0; i< numk_par; i++)
	for(j=0; j< numk_perp; j++)
	  for(n=0; n< dim; n++)
	    {
	      fread(&Garr[n][i][j], sizeof(double), 1, inf);

	      if(n%2==0)
		{
		  if(i !=0 || j != 0) //a strange convention I ended up using for the zero mode
		    Garr[n][i][j]*= pow(kparvec[i]*kparvec[i]+kperpvec[j]*kperpvec[j], .75);
		}
	    } 



      for(n=0; n< dim; n++)    
	splie2dgsl(kparvec, kperpvec, Garr[n], Garr2[n], numk_par, numk_perp);


      fclose(inf);

      fprintf(stdout, "#Read in growth grid.\n");

      return;
    }else if(flag == 2)
    {
      fprintf(stdout, "Freeing growthfac grid.\n");
      free(kperpvec);
      free(kparvec); 

      for(n=0; n< dim; n++)    
	for(i=0; i< numk_par; i++)
	{
	  free(Garr[n][i]);
	  gsl_spline_free (Garr2[n][i]);
	}
      return;
    }

  double kpar, kperp;

  if(numk_par >1){
    kpar = k*costh; 
    kperp = k*sqrt(1.-costh*costh); //note different usage
  }else{
    kpar = 0; 
    kperp = k; 
  }


    if((numk_par!=1 && kpar < kparvec[0]) || kpar > kparvec[numk_par-1] || kperp < kperpvec[0] || kperp > kperpvec[numk_perp-1])
    {
      fprintf(stderr, "kpar, kperp is not in bounds %le %le (%le %le %le %le).\n", kpar, kperp, 
	      kparvec[0], kparvec[numk_par-1],  kperpvec[0], kperpvec[numk_perp-1]);
      exit(-5);
    }

  for(n=0; n< dim; n++) 
    {
      splin2dgsl(kparvec, kperpvec, Garr2[n], acc[n], numk_par, numk_perp, kpar,
	      kperp, G+n+1);

      if(n%2==0)
      	G[n+1] /= pow(k, 1.5);
    }

   if(numk_par== 1){G[11] =0; G[12]=0;} //this case does not handle temperature
}

/* Function to create 2D spline table.  x1a,x2a are the dependent variable
 * vectors.  ya is the matrix of function evaluations, of size m x n.
 * y2a is where the spline table is put. */
void splie2dgsl(double x1a[], double x2a[], double **ya, gsl_spline **spline, int m, int n)
{
  int j;

  for (j=0;j<m;j++)
    {
      gsl_spline_init(spline[j], x2a, ya[j], n);
    }
    //    spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
  return;
}

/* Function to perform spline interpolation in 2D.  x1a, x2a are 
 * dependent variable vectors of length m and n, respectively.  ya is the
 * table of function evaluations.  y2a is the spline table, from a 
 * previous call to splie2.  x1, x2 is the location of the interpolation
 * point, and the result is returned in y. */
void splin2dgsl(double x1a[], double x2a[], gsl_spline **spline, gsl_interp_accel **acc, int m,
	    int n, double x1, double x2, double *y)
{
  int j;
  double *yytmp;

  if(m==1) //if format has kpar=1
    {
      *y= gsl_spline_eval(spline[0], sqrt(x1*x1+x2*x2), acc[0]);
      return;
    }

  gsl_spline *spl = gsl_spline_alloc (gsl_interp_cspline, m);
  gsl_interp_accel *acc1 = gsl_interp_accel_alloc(); //I'm sure this is slow

  yytmp = (double *) malloc(sizeof(double)*(m));
  for (j=0;j<m;j++)
    yytmp[j] = gsl_spline_eval(spline[j], x2, acc[j]);
  //splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);

  gsl_spline_init(spl, x1a, yytmp, m);
    //  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
  *y = gsl_spline_eval(spl, x1, acc1);
//splint(x1a,yytmp,ytmp,m,x1,y);
  free(yytmp);
  gsl_spline_free (spl);
  gsl_interp_accel_free(acc1);
}


#else

/*****************************************************************************
Below is my code that works with numerical recipes rather than GSL libraries.  Please ignore it.
 *****************************************************************************/

/**********************************************************
#include "dnumrecipes.h"
WILL DELETE THIS ONCE I'VE TESTED THINGS A BIT MORE!!! 

Program to return growth factors
flag == 0 initializes
********************************************************/
void interpPk(double k, double costh, double G[], void *unused, int flag)
{
  int i, j, n, s;
  static int numk_par, numk_perp, dim;
  static double *kparvec, *kperpvec, **Garr[20], **Garr2[20];

  if(flag == 0)//read in theta grid
    {    
      FILE *inf;
      char filename[200];


      s= (128 > SIZE) ? 128 : SIZE;
      sprintf(filename, "transfer_func_vbc/initSimCartZI%.1lf_Vbc%.1lf_%d_%.1lf.dat", ZINIT, VSTREAM, s, BOXSIZE);
      fprintf(stdout, "Using %s!!!\n", filename);

      if((inf = fopen(filename, "r")) == NULL)
	{
	  fprintf(stderr, "Could not open %s\n", filename);
	  exit(-5);
	}
      fprintf(stdout, "#Opened %s\n", filename);

      fread(&dim, sizeof(int), 1, inf);  
      fread(&numk_par, sizeof(int), 1, inf);
      fread(&numk_perp, sizeof(int), 1, inf);
	    
      //  fprintf(stderr, "#dim %d numkpar %d numkperp %d\n", dim, numk_par, numk_perp);

      kparvec = dvector(1, numk_par);
      kperpvec = dvector(1, numk_perp);

      fread(kparvec+1, sizeof(double), numk_par, inf);

      fread(kperpvec+1, sizeof(double), numk_perp, inf);

      fread(&TAVG, sizeof(double), 1, inf);
      fread(&NEAVG, sizeof(double), 1, inf);
	    
      fprintf(stdout, "#dim %d numk_par %d numk_perp %d TAVG %g NEAVG %g\n", dim, numk_par, 
	      numk_perp, TAVG, NEAVG);
      
      for(n=0; n< dim; n++)    
	{
	  Garr[n] =  dmatrix(1, numk_par, 1, numk_perp);
	  Garr2[n] =  dmatrix(1, numk_par, 1, numk_perp);
	}

      for(i=1; i<= numk_par; i++)
	for(j=1; j<= numk_perp; j++)
	  for(n=0; n< dim; n++)
	    {
	      fread(&Garr[n][i][j], sizeof(double), 1, inf);

	      //if(n==0 && j <=3 && i <= 3)
	      //	printf("%g %g %le\n", kparvec[i], kperpvec[j], Garr[n][i][j]);

	      if(n%2==0)
		{
		  if(i !=1 || j != 1) //I set k = 1 in this case
		    Garr[n][i][j]*= pow(kparvec[i]*kparvec[i]+kperpvec[j]*kperpvec[j], .75);
		}
	    } 

    
      for(n=0; n< dim; n++)    
	splie2d(kparvec, kperpvec, Garr[n], numk_par, numk_perp, Garr2[n]);

      fclose(inf);

      fprintf(stdout, "#Read in growth grid.\n");

      return;
    }else if(flag == 2)
    {
      fprintf(stdout, "Freeing growthfac grid.\n");
      free_dvector(kperpvec, 1, numk_perp);
      free_dvector(kparvec, 1, numk_par); 

      for(n=0; n< dim; n++)    
	{
	  free_dmatrix(Garr[n], 1, numk_par, 1, numk_perp);
	  free_dmatrix(Garr2[n], 1, numk_par, 1, numk_perp);
	}
      return;
    }

  double kpar = k*costh; 
  double kperp = k*sqrt(1.-costh*costh); //note different usage

 

  if(kpar < kparvec[1] || kpar > kparvec[numk_par] || kperp < kperpvec[1] || kperp > kperpvec[numk_perp])
    {
      fprintf(stderr, "kpar, kperp is not in bounds %le %le (%le %le %le %le).\n", kpar, kperp, 
	      kparvec[1], kparvec[numk_par],  kperpvec[1], kperpvec[numk_perp]);
      exit(-5);
    }

  for(n=0; n< dim; n++) 
    {
      splin2d(kparvec, kperpvec, Garr[n], Garr2[n], numk_par, numk_perp, kpar,
	      kperp, G+n+1);

      if(n%2==0)
      	G[n+1] /= pow(k, 1.5);
      //if(k< 100)
      //printf("vals %g %g %g %g %g\n", logk, k, costh, G[0], G[1]);
    }
}

#endif
