#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "init_cond.h"

extern double BOXSIZE, ZINIT, Mpart[2];
extern int SIZE;
extern char *BASEOUT;


/**************************gridAndPrintPowerOfDisplacements*****************************
flm == 0 is dark matter
flm == 1 is baryons
****************************************************************************************/
void gridAndPrintPowerOfDisplacements(struct dispGrid *grid, int N, int Nf, int flm)
{
  int i;
#ifndef NGP
  int ii, jj, kk, ll;
#endif
  double rho = Mpart[flm], vm;//really is mass  
  double *dens = (double *) malloc(sizeof(double)*N);
  double *vel = (double *) malloc(sizeof(double)*N);

  for(i=0; i< N; i++)
    {
      dens[i] = 0.; vel[i]=0.;
    }


  //put particles on grid with NGP/CIC
  for(i=0; i< N; i++)
    {      

#ifndef NGP
      vm = grid->v[flm][0][i]*rho;
      gridCIC(grid->dis[flm][0][i], grid->dis[flm][1][i], grid->dis[flm][2][i], &rho, &dens, 1);
      gridCIC(grid->dis[flm][0][i], grid->dis[flm][1][i], grid->dis[flm][2][i], &(vm), &vel, 1);
#else
      ii = ((int) (grid->dis[flm][0][i]/BOXSIZE*SIZE + SIZE))%SIZE;
      jj = ((int) (grid->dis[flm][1][i]/BOXSIZE*SIZE + SIZE))%SIZE;
      kk = ((int) (grid->dis[flm][2][i]/BOXSIZE*SIZE + SIZE))%SIZE;

      ll = ii*SIZE*SIZE + jj*SIZE + kk;
      dens[ll] += rho;  //in Msun/Mpc^3 h^2
      vel[ll] += grid->v[flm][0][i]*rho; 
#endif
 
    }

  for(i=0; i< N; i++)
    if(dens[i] > 1e-20)
      vel[i]/= dens[i];

  free(dens); free(vel);
}



/*************************************************************************************************
Routines below print the power spectrum of whatever field one wants.  An example
of how to call these routines is shown at the end of the function main().  These
are primarily useful for debugging
 ************************************************************************************************/

void makePk(double *K, double *Pk, int *Count, fftw_complex *Rf,  
	    fftw_complex *Rf2)
{
  int i, j, k, l, a, n;
  double m;
  double dt = 2.*M_PI/BOXSIZE;
  double volume;

  for(n = 0; n< SIZE*sqrt(3)/2.; n++)
    {
      K[n] = 0;
      Pk[n] = 0;
      Count[n] = 0;
    }
  volume = BOXSIZE*BOXSIZE*BOXSIZE;
  for(i =  0; i< SIZE; i++)
    for(k =  0; k< SIZE; k++)
      for(j =  0; j< SIZE/2 + 1; j++)
	{
	  a = 2;
	  if(j == 0 || j == SIZE/2 + 1)
	    a = 1;
	  l = (i*SIZE + k)*(SIZE/2 +1) + j;
	  m = sqrt(knum(i)*knum(i)+ knum(k)*knum(k) + j*j);
	  n = (int) (m + .5);
	  K[n] += a*dt*m;
	  Pk[n] += a*(Rf[l][0]*Rf2[l][0] + Rf[l][1]*Rf2[l][1]);
	  Count[n]+= a;
	}
 
  for(n = 0; n< SIZE*sqrt(3)/2.; n++)
    {
      K[n] /= Count[n];
      Pk[n] /= (Count[n]*volume);
  }
}


void printStatistics(double *K, double *Pk, char *out) 
{
  int i;
  double coef = 1.;
  FILE *outfile;
  char filename[200];

  coef = 1./(2.*M_PI*M_PI);

  sprintf(filename, "%s.dat", out);
  printf("OUTFILE = %s\n", filename);
  if((outfile = fopen(filename, "w")) == NULL)
  {
      fprintf(stderr, "Could not open output file.\n");
      exit(1);
  }
  for(i =  1; i< SIZE*sqrt(3.)/2; i++)
      fprintf(outfile, "%e %e\n", K[i], coef*pow(K[i],3.)*Pk[i]);
    
  fclose(outfile);
}

void initialize(double **K,double **Pk, int **Count, int flag)
{
  int i;

  if(flag == 0)
  {
    *K = (double *) malloc(sizeof(double)*(1e4)); //to be safe
    *Pk = (double *) malloc(sizeof(double)*(1e4));
    *Count = (int *) malloc(sizeof(int)*(1e4));

    for(i =  0; i< 1e4; i++)
      {
	  Pk[0][i] = 0.;
	  Count[0][i] = 0;
	  K[0][i] = 0.;
      }
  }else {
      free(*K);
      free(*Pk);
      free(*Count);
  }
}
