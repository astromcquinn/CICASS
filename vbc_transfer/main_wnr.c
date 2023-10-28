/*********************************************************************
Generates transfer function used by initial conditions code
or reads in CAMB file
*********************************************************************/

#define MPCTOKM 3.086e19
#define MPCTOCM 3.086e24
#define DELTAC 1
#define DELTAB 5
#define REAL 0
#define IMAG 1
#define VELREAL 2
#define VELIMAG 3

#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include "dnumrecipes.h"

double setInterpTF(double k, double z, double *TFdir, int flag);
double interpolateRecFast(double z, double *xe, int flag);
void initialConditionsCartesian(double kmin);
double getPowerZ(double k, double z, double G[], int flag);
double TF_Exact(double k);
double alphaB_recomb(double T, int species);
void returnGrowth(double G[], double z, double k, double costh);
void setICs(double *y, double *dy, double k, double costh, double vstream);
void setCalcGrowthInt(double a, double *y, double *deriv, double k1, 
		      double vstream1, double costh1,  int flag);
void calcGrowthInt(double a, double *y, double *dy);
double hubbleZ(double z);
void resetPowerSpectrum();


const double TCMB(2.726);                   // CMB temperature at z=0
const double MPROTON(1.6726e-24);           // in grams
const double MELECTRON(9.11e-28);
const double MEL_EV(5.11e5);
const double LIGHTSPEED(3.0e10);

const double CRITDENSITY(1.8791e-29);   // Current critical density, 
                                        // in g/cm^3, not including h^2!
const double SIGMAT(0.665e-24);   //Thomson cross section in cm^2
const double YHE(0.25);
const double BOLTZK(1.3806e-16);            // in erg/K
const double MUB(1.22);           // Mean molecular weight, for primordial
const double aSB(7.56e-15); //u = aSB T^4 

/*******These all will be externs***************/
double BOXSIZE = 0.2;
double SIZE =32;
double VSTREAM = 0.;
double ZSTART = 1000.;     //The redshift VSTREAM is valid at
double ZINIT = 50;

double hconst = 0.71; 
double OmegaM = 0.27;        //assumes flat universe
double OmegaB = 0.046;
double Sigma8 = 0.8;
double Tilt = 0.95;
double OmegaR =  4.15e-5/(hconst*hconst); //assumes relativistic neutrinos, 
                                // from Dodelson, eqn 2.87


#ifndef ISOTHERMAL
int  NUMEQNS = 12;
#else
int  NUMEQNS = 8; //not doing temperature yet
#endif 

double sNorm; //sigma8


using namespace std;

int main()
{
  interpolateRecFast(0., NULL, 0);
  setInterpTF(0, 0, NULL, 0);
  resetPowerSpectrum();

  double G[100], k=1, costh=1;

  // cout << sNorm << endl;

  for(k =.1; k< 10000; k*=1.2)
    {
      returnGrowth(G, ZINIT, k, costh);
      cout << k << " " << k*k*k/19.73*getPowerZ(k, 50., G, 1) << " " << k*k*k/19.73*getPowerZ(k, 50., G, 2)  << endl;
    }
  //initialConditionsCartesian(2.*PI/BOXSIZE);

  return 1;
}



void initialConditionsCartesian(double kmin)
{
  int NUMK_VALS = SIZE;
  int NUMK_PAR = (NUMK_VALS + 1);
  int NUMK_PERP = ((int) ((NUMK_VALS + 1)*sqrt(1.25)));  //THIS IS ODD!!!
  double   boxsize = BOXSIZE;

  int i, j, n, DIM = NUMEQNS;
  double costh, Tk, xe, k, g[20], gold[20], growth[20];
  double *kvec_par, *kvec_perp;    
  char filename[200];
  double delta0, deltarmsvec[100];
  FILE *outfile;


  kvec_par = (double *) malloc(sizeof(double)*NUMK_PAR);
  kvec_perp = (double *) malloc(sizeof(double)*NUMK_PERP);

  sprintf(filename, "initSimCartZI%.1lf_Vbc%.1lf_%d_%.1lf.dat", ZINIT, VSTREAM, NUMK_VALS, boxsize); 
  if((outfile = fopen(filename, "w"))==NULL)
    {
      fprintf(stderr, "makeGFGrid Cartesian:  Could not open %s for writing.\n",
	      filename);
      exit(-9);
    }
  fprintf(stdout, "#makeGFGrid  Cartesian: opened %s\n", filename);

  for(i=0; i < NUMK_PAR; i++)
    kvec_par[i] = 2.*M_PI/boxsize*(i);//let this one go to zero 
  for(i=0; i < NUMK_PERP; i++)
    kvec_perp[i] = 2.*M_PI/boxsize*(i);
    
  fwrite(&DIM, sizeof(int), 1, outfile);
  fwrite(&(i = NUMK_PAR), sizeof(int), 1, outfile);
  fwrite(&(i = NUMK_PERP), sizeof(int), 1, outfile);
  fwrite(kvec_par, sizeof(double), NUMK_PAR, outfile);
  fwrite(kvec_perp, sizeof(double), NUMK_PERP, outfile);

  Tk = interpolateRecFast(ZINIT, &xe, 1);
  cout << "# " << Tk << " " << xe <<  endl;
  fwrite(&Tk, sizeof(double), 1, outfile);
  fwrite(&xe, sizeof(double), 1, outfile);
  
  for(i=0; i < NUMK_PAR; i++)
    {  
      for(j=0; j< NUMK_PERP;j++)
	{
	  k= sqrt(kvec_par[i]*kvec_par[i]+kvec_perp[j]*kvec_perp[j]);

	  if(i == 0 && j==0) k = kmin; 

	  costh = kvec_par[i]/k;

	  delta0 = sqrt(getPowerZ(k*hconst, ZSTART, NULL, 1)); //dark matter power

	  returnGrowth(growth, ZINIT, k*hconst, costh);

	  /* cout << i << " " << j << " k = " << k << " " << costh << " " << growth[5] << " " << growth[6] 
	     << " " << growth[7] << " " << growth[8] <<endl;*/
	  //cout << kvec_par[i] << " " << kvec_p[i] << endl;

	  //convert to amplitude and phase here

	  for(n=0;n<DIM;n+=2)//interpolate using amplitude and phase
	    {
	      g[n] = sqrt(growth[n+1]*growth[n+1]+ growth[n+2]*growth[n+2]);
	      g[n+1] = acos(growth[n+1]/g[n]);
	      //  if(n == 0)
	      //cout << " " << g[n] << " " << g[n+1] << endl;
	      if(growth[n+2] < 0)
		g[n+1] = -g[n+1];

	      if(j!= 0)  //to make interpolation continous
		{
		  if(g[n+1] - gold[n+1] > .6)
		    {
		      g[n+1] -= 2.*PI;
		    }else if(g[n+1] - gold[n+1] < -.6)
		    {
		      g[n+1] += 2.*PI;
		    }

		  if(fabs(g[n+1] - gold[n+1]) > .6)
		    {
		      fprintf(stderr, "Too big of a jump in phase of %le at %d %d %d.  Try more k-bins!\n", 
			      fabs(g[n+1] - gold[n+1]), n, i, j);
		      cout <<  k << " " << g[n+1] << " " << gold[n+1] << endl;
		      // exit(-5);
		    }
		}
	      gold[n]=g[n]; gold[n+1]=g[n+1];
	    }

	  for(n=0; n< DIM; n+=2)//only mutiply those that aren't phases
	    {	
	      deltarmsvec[n] = g[n]*delta0*pow(hconst, 1.5);
	      deltarmsvec[n+1] = g[n+1];
	      //    cout << delta0 << " " << deltarmsvec[n] << endl; exit(-5);
	    }

	  fwrite(deltarmsvec, sizeof(double), DIM, outfile);
	}
    }

  fclose(outfile);
   
  cout << "#Done outputing " << filename << " for IC generator." << endl; 

  free( kvec_par );  free( kvec_perp ); 
  return;
}




/***************************************************************
calculates growth factor as a function of z, k and costh
 **************************************************************/
void returnGrowth(double G[], double z, double k, double costh)
{
  double *y = dvector(1, NUMEQNS), *dy = dvector(1, NUMEQNS);
  int i, goodSteps,badSteps;
  
  setICs(y, dy, k, costh, VSTREAM); //set initial conditions

  setCalcGrowthInt(0., NULL, NULL, k, VSTREAM, costh, 0);

  odeint(y, NUMEQNS, 1./(1.+ZSTART), 1./(1.+ z), 1e-3, 1e-6, 0.,
	 &goodSteps,&badSteps, calcGrowthInt, bsstep); //ODEs


  //fprintf(stderr, "G1 = %le\n", y[DELTAC]);

  for(i=1; i<= NUMEQNS; i++)
    G[i] = y[i];

  free_dvector(y, 1, NUMEQNS); free_dvector(dy, 1, NUMEQNS);
}


/************************************************************************
sets initial conditions for calculation
 **********************************************************************/
void setICs(double *y, double *dy, double k, double costh, double vstream)
{
  const double DELTA0 = 1.; //sets amplitude to unity

  double a = 1./(1+ZSTART);
  double H = 100.*hconst*hubbleZ(ZSTART)/MPCTOKM;
  double deltab, deltacdot, deltabdot;  
  double TFc, TFdirc, TFb, TFdirb; 
  TFc = setInterpTF(k, 0., &TFdirc, 6); //dz of dark matter at z=1000
  TFb = setInterpTF(k, 0., &TFdirb, 7); //dz of dark matter at z=1000
  deltab = TFb/TFc*DELTA0;
  deltacdot = TFdirc/TFc*H*(1+ZSTART)*DELTA0;
  deltabdot = TFdirb/TFb*deltab*H*(1+ZSTART);

  /********************Initial Conditions*************************/
  y[1] = DELTA0; //real part of dark matter density
  y[2] = 0.;  //imaginary part
  y[3] = -deltacdot; 
  y[4] = 0.; //after much thought, phase falls out so this is correct
  y[5] = deltab;  //baryonic delta, set to zero 
  y[6] = 0.;
  y[7] = -deltabdot; //baryonic velocity
  y[8] = 0.;

#ifndef ISOTHERMAL
  y[9]  = TCMB*(1+ZSTART); //temperature -- unused by other eqns unless CALCULATE_MEAN_TEMP
  y[10] = 1.;  //electron fraction after recombination.  
               /***********Even though this is a total kludge because I don't include radition, I find that this does pretty well--electron       
                   fraction is enough of an attractor solution.***********************************************************/
  y[11] = 0.; // initially assumed isothermal (which should be a good assumption)
  y[12] = 0.;
#endif
}

/***********************************************************
Dummy function for differential equation solver
 **********************************************************/
void calcGrowthInt(double a, double *y, double *dy)
{
  setCalcGrowthInt(a, y, dy, 0., 0., 0., 1);
}




/****************************************************************
solves for dy/dt for delta_c, delta_b, etc 
 ***************************************************************/
void setCalcGrowthInt(double a, double *y, double *deriv, double k1, 
		      double vstream1, double costh1,  int flag)
{
  const double fHe = .25*YHE/((1. - YHE) + .25*YHE);
  const double a1 = 0.008403, a2 = 0.008696; //fit to recfast in TH11--eventually remove
  static double k, costh, h, OC0, OB0, OR0, nH, ucmb0, tgamma0, vstream;
  double H, vbck, OMEGAC, OMEGAB, cs, T, xe, eta1;
  int i;
  double xe_z(double z);//poor man's interpolation--eventually remove

  if(flag == 0)
    {
      k = k1; costh = costh1; vstream = vstream1;
      OC0 = OmegaM-OmegaB;
      OB0 = OmegaB;
      h = hconst;
      OR0 = OmegaR;
 

      nH = CRITDENSITY*OB0*h*h*(1. - YHE)/(MPROTON);  //Hydrogen density
                              //hydrogen density, assuming helium is neutral
      ucmb0 = aSB*TCMB*TCMB*TCMB*TCMB;
      tgamma0 = 3.*MELECTRON*LIGHTSPEED/(8.*SIGMAT*ucmb0);
      //    fprintf(stdout, "Initalized with %le %le %le %le %le\n",
      //	      k, costh, OC0, OB0, h);
      return;
    }

  // H = 100.*h*sqrt((OC0+OB0)/(a*a*a))/MPCTOKM;
  H = 100.*hconst*hubbleZ(1./a - 1.)/MPCTOKM;
  vbck = vstream/MPCTOKM/(a*(1.+ZSTART)); //VSTREAM is in km/s
  OMEGAC = OC0/(OC0 +OB0 + OR0/a);
  OMEGAB = OB0/(OC0 +OB0 + OR0/a);

#ifdef TURNONVBK
  //if(a < 1./101) vbck = 0.;
  vbck*=exp(-pow(101.*a, -500.)); //goes to zero abvoe z = 200 (currently)
#endif

#ifdef CALCULATE_MEAN_TEMP
  T = y[9];
  xe = y[10];
#elif RECFAST //fit in TH
  T = interpolateRecFast(1./a-1., &xe, 1);
#elif REIONIZATION 
  if(a > 0.1)
    {
      xe = 1;
      T = y[9] = 1e4;  
    }
  else
    T = interpolateRecFast(1./a-1., &xe, 1);
#else  //use Hirata interpolation formula
  T= 2.726/a/(1+ a/a1/(1+ a2/a*sqrt(a2/a)));
  // interpolateRecFast(1./a-1., &xe, 1);//try to switch to fit

  xe = y[10];
  //cout << 1./a-1 << " " << T << " " << y[9] <<  " " << interpolateRecFast(1./a-1., NULL, 1) << endl; 
#endif

  cs = sqrt(BOLTZK*T/MUB/MPROTON)/MPCTOCM; //isothermal sound speed (What appears in eqn)

  // cout << vbck/cs << endl;
  // cout << H << " " << vbck << " " << OMEGAC << " " << sqrt(BOLTZK*T/MUB/MPROTON) << endl;

  if(NUMEQNS == 4) { y[5] = 0.; y[6] = 0.;}
  if(NUMEQNS < 8){y[11] = 0.; y[12] = 0.;}

  /*************************DARK MATTER*********************************/
  deriv[DELTAC+REAL] =  -y[DELTAC+VELREAL];

  //cout <<a << " " << k << " " << vbck << " " << costh << " " <<  vbck*k*costh*y[DELTAC+REAL]/a << endl;

  // fprintf(stdout, "hi %d %e %e\n",  DELTAC+VELREAL, deriv[DELTAC+REAL],  y[DELTAC+VELREAL]);
  deriv[DELTAC+IMAG] =  -y[DELTAC+VELIMAG];
  deriv[DELTAC+VELREAL] =  
    - 1.5*H*H*(OMEGAC*y[DELTAC+REAL]+OMEGAB*y[DELTAB+REAL])
    - 2.*H*y[DELTAC+VELREAL];
  deriv[DELTAC+VELIMAG] =  
    - 1.5*H*H*(OMEGAC*y[DELTAC+IMAG]+OMEGAB*y[DELTAB+IMAG])
    - 2.*H*y[DELTAC+VELIMAG];

//baryons
  if(NUMEQNS != 4)
    {
      deriv[DELTAB+REAL] = vbck*k*costh*y[DELTAB+IMAG]/a - y[DELTAB+VELREAL];
      deriv[DELTAB+IMAG] = -vbck*k*costh*y[DELTAB+REAL]/a - y[DELTAB+VELIMAG];
      deriv[DELTAB+VELREAL] =  vbck*k*costh*y[DELTAB+VELIMAG]/a 
	- 1.5*H*H*(OMEGAC*y[DELTAC+REAL]+OMEGAB*y[DELTAB+REAL])
	- 2.*H*y[DELTAB+VELREAL] + cs*cs*k*k/(a*a)*(y[DELTAB+REAL] + y[11]);
      deriv[DELTAB+VELIMAG] =  -vbck*k*costh*y[DELTAB+VELREAL]/a  
	- 1.5*H*H*(OMEGAC*y[DELTAC+IMAG]+OMEGAB*y[DELTAB+IMAG])
	- 2.*H*y[DELTAB+VELIMAG] + cs*cs*k*k/(a*a)*(y[DELTAB+IMAG] + y[12]); 
    

#ifndef ISOTHERMAL //evolve temperature and electrons
      eta1 = (1.+ fHe + y[10]);

      //EQNS 9 and 10 are not currently used
      deriv[9] =-2.*H*y[9] + y[10]/(eta1*tgamma0*(a*a*a*a))
		   *(TCMB/a - y[9]);//mean temperature (adiabatic and compton cooling)
 
      deriv[10] = -alphaB_recomb(T, 0)*y[10]*y[10]*nH/(a*a*a); 
      deriv[11] =  vbck*k*costh*y[12]/a + 2./3.*(deriv[DELTAB+REAL]- vbck*k*costh*y[DELTAB+IMAG]/a) - xe/(eta1*tgamma0*(a*a*a*a))
	*(TCMB/a/T)*y[11]; //temperature fluctuations --uses either recfast temp or y[9]
      deriv[12] = -vbck*k*costh*y[11]/a + 2./3.*(deriv[DELTAB+IMAG] + vbck*k*costh*y[DELTAB+REAL]/a) - xe/(eta1*tgamma0*(a*a*a*a))
	*(TCMB/a/T)*y[12]; //temperature fluctuations (imaginary)

      /******Currently do not track electron fraction fluctuations***********/
#endif
    }
  
  for(i=1; i<=NUMEQNS; i++)
    {
      deriv[i] /= (H*a); //to switch from time to scale factor derivative
      //cout << i << " " << y[i] << " " <<deriv[i] <<  endl;
      //   fprintf(stdout, "%d %le %le %le\n", i, y[i], deriv[i], a);
    }
}



/***********************************************************
Reads in CAMB TFs for z= 0 and 1000
currently z does nothing -- should solve for growth factors
Tdir gives the derivative at z= 1000
 ***********************************************************/
double setInterpTF(double k, double z, double *TFdir, int flag)
{
  int i, j;
  static int init = 0, numk, numk2;
  static double *lkarr, *lkarr2, *T[12], *T2[12], *Tdir[12];
  const double natural(1.0e30);
  double TF;
  

  if(init == 0)
    {      
      cout << "#Initializing TFs" << endl;
      init++;
      int numk3, numk4;
      double *lkarr3;
      double t[10];

      //temporary
      char *BASETRANSFER = "initSB_transfer_out";
      int z = (int) 1000, dz=3;

      char filename0[200], filename1[200], filename1p[200], filename1m[200];
      sprintf(filename0, "TFs/%s_z%03d.dat", BASETRANSFER, 0);
      sprintf(filename1, "TFs/%s_z%03d.dat", BASETRANSFER, z);
      sprintf(filename1m, "TFs/%s_z%03d.dat", BASETRANSFER, z-dz);
      sprintf(filename1p, "TFs/%s_z%03d.dat", BASETRANSFER, z+dz);

      // char *infz0_name = (char *)"TFs/initSB_transfer_out_z000.dat";
      FILE *infz0 = fopen(filename0, "r");
      //char *infz1000_name = (char *)"TFs/initSB_transfer_out_z1000.dat";
      FILE *infzinit = fopen(filename1, "r");
      //char *infz1003_name = (char *)"TFs/initSB_transfer_out_z1003.dat";
      FILE *infzinitp = fopen(filename1p, "r");
      // char *infz997_name = (char *)"TFs/initSB_transfer_out_z997.dat";
      FILE *infzinitm = fopen(filename1m, "r");

      if(infz0 == NULL || infzinit == NULL || infzinitp == NULL || infzinitm == NULL)
	{
	  fprintf(stderr, "Could not open one of transfer function files:  %s %s %s %s\n", 
		  filename0, filename1, filename1m, filename1p);
	  exit(-6);
	}


     fscanf(infz0, "%d\n", &numk);
      fscanf(infzinit, "%d\n", &numk2);
      fscanf(infzinitp, "%d\n", &numk3);
      fscanf(infzinitm, "%d\n", &numk4);

      //cout << numk << " " << numk2 << " " << numk3 << endl;

      if(numk !=numk2 || numk != numk3 || numk3 != numk4){
	cerr << "Arrays not equal!!!\n" << endl; 
	exit(-6);
      }

      lkarr = dvector(1, numk);
      lkarr2 = dvector(1, numk);
      lkarr3 = dvector(1, numk);      
      for(i=0;i<12;i++)
	{
	  T[i] = dvector(1, numk);
	  T2[i] = dvector(1, numk);
	  Tdir[i] = dvector(1, numk);
	}

      /***********************************Read in TRANSFERFUCNTIONS**********************/
      for(i=1;i<=numk;i++)
	{
	  fscanf(infz0, "   %lg   %lg   %lg   %lg   %lg   %lg   %lg\n", lkarr+i,
	       T[0]+i, T[1]+i, T[2]+i,T[3]+i,T[4]+i,T[5]+i);
	  fscanf(infzinit, "%le %le %le %le %le %le %le\n", lkarr2+i,
	       T2[0]+i, T2[1]+i, T2[2]+i,T2[3]+i,T2[4]+i,T2[5]+i);
	  fscanf(infzinitp, "%le %le %le %le %le %le %le\n", lkarr3+i,
	       Tdir[0]+i, Tdir[1]+i, Tdir[2]+i,Tdir[3]+i,Tdir[4]+i,Tdir[5]+i);
	  fscanf(infzinitm, "%le %le %le %le %le %le %le\n", &t[9],
	       &t[0], &t[1], &t[2], &t[3], &t[4], &t[5]);
	  for(j=0;j<6;j++)
	    Tdir[j][i] = (t[j] - Tdir[j][i])/(2.*dz);  //Not terribly general

	  lkarr[i] = log10(lkarr[i]);
	  lkarr2[i] = log10(lkarr2[i]);
	}


      
      for(i=0;i<6;i++)
	{
	  spline(lkarr, T[i], numk, natural, natural, T[6+i]);  //change to GSL
	  spline(lkarr2, T2[i], numk, natural, natural, T2[6+i]);
	  spline(lkarr2, Tdir[i], numk, natural, natural, Tdir[6+i]);
	}
      cout << "#done!" << endl;

      fclose(infz0); fclose(infzinit); fclose(infzinitm); fclose(infzinitp);
      free_dvector(lkarr3, 1,numk3);
      return -1.;
    }else if(flag == 2){
    free_dvector(lkarr, 1,numk2);
    free_dvector(lkarr2, 1,numk2);

    for(i=0;i<12;i++)
      {
	free_dvector(T[i], 1, numk);
	free_dvector(T2[i], 1, numk);
	free_dvector(Tdir[i], 1, numk);
      }
    return 1;
  }
  
  double lgk = log10(k/hconst); //to put in h/Mpc

  if(lkarr[numk] < lgk || lkarr[1] > lgk)
    {
      fprintf(stderr, "k is out of range\n");
	 return 0.;
    }  

  if(flag >=6)//z=1000 TF and dirivative
    {
      splint(lkarr2, T2[flag-6], T2[flag], numk2, lgk, &TF);  
      if(TFdir != NULL)
	{
	  splint(lkarr2, Tdir[flag-6], Tdir[flag], numk2, lgk, TFdir);  
	}
    }
  else //otherwise
    splint(lkarr, T[flag], T[6+flag], numk, lgk, &TF);
 
  return TF;
}


/***********************************************
Returns density power spectrum
flag = 0 total
flag = 1 dark matter
flag = 2 baryon
flag = 3 temperature
*************************************************/
double getPowerZ(double k, double z, double G[], int flag)
{
  double temp, g1, g2, g3;

  /***********Initializes directly from z=1000 transfer function*********/
  if((z -1000.)*(z -1000) < 1e-5)//will use this part to init sims
    {

      if(flag == 1)  //dark matter (z=1000)
	temp = 2.0*pow(PI*sNorm*setInterpTF(k, 0., NULL, 6),2.0);
      else if(flag == 2) //baryons
	temp = 2.0*pow(PI*sNorm*setInterpTF(k, 0., NULL, 7),2.0);
      else if(flag == 0) //total
	temp = 2.0*pow(PI*sNorm*setInterpTF(k, 0., NULL, 11),2.0);	
      else
	{
	  cout << "#powerZ: This option does not work here" << endl;
	  exit(-5);
	}
    }else
    {	 
      g1= sqrt(G[2]*G[2]+G[1]*G[1]); //dm
      g2 =  sqrt(G[5]*G[5]+G[6]*G[6]);//baryons
      g3 =  sqrt(G[11]*G[11]+G[12]*G[12]);//temperature

      if(flag == 0) //tot
	{
	  double fb = OmegaB/OmegaM;
	  temp = 2.0*pow(PI*sNorm*(g1*(1.-fb)*setInterpTF(k, 0., NULL, 6)
				   +g2*fb*setInterpTF(k, 0., NULL, 6)),2.0);
	}
      else if (flag == 1)
	temp = 2.0*pow(PI*sNorm*g1*setInterpTF(k, 0., NULL, 6), 2.);
      else if (flag == 2)
	temp = 2.0*pow(PI*sNorm*g2*setInterpTF(k, 0., NULL, 6), 2.);
      else if (flag == 3)
	temp = 2.0*pow(PI*sNorm*g3*setInterpTF(k, 0., NULL, 6), 2.);
      else
	{
	  cout << "powerZ: This option does not work here" << endl;
	  exit(-6);
	}
    }

  temp *= pow(k,Tilt);
  return temp;
}


//returns total matter z= 0 transfer function using CAMB
double TF_Exact(double k)
{
  // double setInterpTF(double k, double z, double *Tdir, int flag);
 
  return setInterpTF(k, 0., NULL, 5); 
}


/********************************************************
Interpolates a recfast file to return temperature
(and a pointer to the electron fraction) as a function of z
*********************************************************/
double interpolateRecFast(double z, double *xe, int flag)
{
  static int numz;
  static double *zarr, *Tarr, *xearr, *Tarr2, *xearr2;
  int i;
  const double natural(1.0e30);

  if(flag == 0)
    {
      char *fname = "recfast/xeTrecfast.out";
      FILE *infile = fopen(fname, "r");
      fscanf(infile, "%d\n", &numz);

      zarr = dvector(1, numz);
      Tarr = dvector(1, numz);
      xearr = dvector(1, numz);
      Tarr2 = dvector(1, numz);
      xearr2 = dvector(1, numz);

      for(i=numz; i>=1; i--) //need to read backwards so z is increasing
	fscanf(infile, "%le %le %le\n", &zarr[i], &xearr[i], &Tarr[i]);

      spline(zarr, Tarr, numz, natural, natural, Tarr2);
      spline(zarr, xearr, numz, natural, natural, xearr2);
      fclose(infile);
      return -1.;
    }
  else if(flag ==2)
    {
      free_dvector(Tarr, 1, numz);
      free_dvector(Tarr2, 1, numz);
      free_dvector(xearr, 1, numz);
      free_dvector(xearr2, 1, numz);
      free_dvector(zarr, 1, numz);
      return -1.;
    }

  double T;

  if(xe != NULL)
    splint(zarr, xearr, xearr2, numz, z, xe);
  splint(zarr, Tarr, Tarr2, numz, z, &T);

  return T;
}


/************************************************************
Returns value of Hubble constant at zCurrent in units of H0. 
************************************************************/
double hubbleZ(double zCurrent)
{
  double temp;

  zCurrent += 1.0;
  temp = ((1.0 - OmegaM-OmegaR) 
	  + OmegaM*pow(zCurrent,3.0) + OmegaR*pow(zCurrent,4.0));
  return sqrt(temp);
}


/********************************************************************
Below are routines to correctly normalize power spectrum
 *******************************************************************/

void resetPowerSpectrum()
{
  double acc,sig8;
  double setSigmatop(double kl, int flag);
  double sigmatop2(double k); double sigmatop(double kl);

  // Normalize sigma8 at present time
  double scale = 8.0/hconst;
  acc = 1.0e-6;
  setSigmatop(0.0,1);
  //Matt: switched to qrombnr as the other was not working properly
  //when I upgraded to lion gcc
  sig8 = qromb(sigmatop2,1.0e-5,0.001/scale,acc);
  sig8 += qromb(sigmatop,log(0.001/scale),log(0.1/scale),acc);
  sig8 += qromb(sigmatop,log(0.1/scale),log(1.0/scale),acc);
  sig8 += qromb(sigmatop,log(1.0/scale),log(10.0/scale),acc); 
  sig8 += qromb(sigmatop,log(10.0/scale),log(100.0/scale),acc);
  sig8 = sqrt(sig8);
  sNorm = Sigma8/sig8;

  cout << "# " << sNorm << endl;
}


/* Calculates sigmatophat, assuming k is entered as ln(k) */
double setSigmatop(double kl, int flag)
{
  double k,x,sigtop;
  double dummy;
  double scale = 8.0/hconst; //for sigma8

  if (flag == 1) {
    return 0.0;
  }

  k = exp(kl);
  x = scale*k;
 
   sigtop = pow(k,3.0+Tilt)*pow(TF_Exact(k),2.0)*
           pow(3.0*(x*cos(x) - sin(x))/pow(x,3.0),2.0);

  return sigtop;
}

double sigmatop(double kl)
{
  return setSigmatop(kl,0);
}

/* Calculates sigmatophat, assuming k is entered as linear. */
double sigmatop2(double k) 
{
  return sigmatop(log(k))/k;
}


//case A and B rec coefficients
//updated H to be valid at ultra low temperatures (using recfast rate)
//currently only H is used
double alphaB_recomb(double T, int species)
{
   double lambda;

   //Hydrogen
   if(species == 0)
     {
     //eqn 70 in RECFAST paper from Hummer (1994) --valid at low temperatures
     const double a = 4.309, b=-0.6166, c=0.6703, d=0.5300;
     return 1.e-13*a*pow(T/1.e4,b)/(1.+c*pow(T/1e4,d));
     }
   
   cerr << "Species not supported!\n" << endl;
   exit(-5);
   return 0.;
}
