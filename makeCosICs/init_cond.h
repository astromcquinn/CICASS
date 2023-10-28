#ifndef INIT_COND_H
#define INIT_COND_H

#define MPCTOKM 3.08568e19
#define TCMB 2.73
#define CRITDENSITYMPC 2.7755e11  // solar masses per Mpc
#define CONVERTTOMASSUNITS_MSUNPH 1e-10    //converts MSUN to our mass unit (same as default GADGET -- 1e10 Msun)

#include <gsl/gsl_spline.h>



struct dispGrid{
  double *delta[2], **dis[2], **v[2], *T;
};

//function prototypes
void  printInputPower();
void setParams(int argc, char *argv[]);
void readInArray(double *R, char *file, double B);
void initialize(double **K,double **Pk, int **Count, int flag);
void printStatistics(double *k, double *Pk, char *out);

void doFFT(double *R, fftw_complex *Rf, long Nf);
int knum(int i);
void doInverseFFT(fftw_complex *Rf, double *R, long N);
void addGaussianField(fftw_complex *phif, fftw_complex *phidotf, fftw_complex *delta, fftw_complex *temp, int flag);
void getNoise(double *L, double *noise,  double time);
void calculateVariance(double *L, double *noise);
void buildFluctuation(double *array, long num);
void outputICs(struct dispGrid *grid);
void makePk(double *K, double *Pk, int *Count, fftw_complex *Rf,  
	    fftw_complex *Rf2);
void setdphi(fftw_complex **dphif, fftw_complex *phif);
double returnPk(double k, double theta, int species, int flag);
void calculateDispField(struct dispGrid *grid, double *phi, double **dphi, double *phi2);
double getHz(double zg);
#ifdef GSL_INTERP
void interpPk(double k, double costh, double G[],  gsl_interp_accel **acc[], int flag);
#else
void interpPk(double k, double costh, double G[], void *unused, int flag);
#endif
void generateDisplacements(int flag);
void calculateDispField(struct dispGrid *grid, double **dphi, double **dphidot, int species);
void generateDisplacements(struct dispGrid *grid, int flag);
void gridAndPrintPowerOfDisplacements(struct dispGrid *grid, int N, int Nf, int flm);
void gridCIC(double u, double v, double w, double val[], double **g, int Naxes);
double periodic_wrap(double x);
void gridWithGlass(struct dispGrid *grid, double **dphi, double **dphidot, int s, double fac, int Np);
void gridDeltabAndTemperatureCIC(dispGrid *grid, int flag);

#ifdef OUTPUT_ENZO
int SaveEnzoSnapshot(double **disdm, double **veldm, double *deltab, double **velb, double *temp);
#endif

#endif
