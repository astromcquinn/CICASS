/**********************************************************************************
Reads gadget glass file and (also) outputs Gadget format

**********************************************************************************/



#define  PI          3.14159265358979323846 
#define  GRAVITY     6.672e-8
#define  HUBBLE      3.2407789e-18   /* in h/sec */
#define  BOLTZMANN   1.3806e-16
#define  PROTONMASS 1.6726e-24

//internal units I use for Gadget
#define UnitLength_in_cm          3.085678e21   /* defines length unit of output (in cm/h) */
#define UnitMass_in_g             1.989e43      /* defines mass unit of output (in g/cm)*/
#define UnitVelocity_in_cm_per_s  1e5           /* defines velocity unit of output (in cm/sec)*/


typedef int int4byte;
typedef unsigned int uint4byte;

struct part_data 
{
  float Pos[3];
  float Vel[3];
  long long ID;
} *P;



struct io_header_1
{
  uint4byte npart[6];      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
  double mass[6];          /*!< mass[1] gives the particle mass */
  double time;             /*!< time (=cosmological scale factor) of snapshot */
  double redshift;         /*!< redshift of snapshot */
  int4byte flag_sfr;       /*!< flags whether star formation is used (not available in L-Gadget2) */
  int4byte flag_feedback;  /*!< flags whether feedback from star formation is included */
  uint4byte npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
                                the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
  int4byte flag_cooling;   /*!< flags whether radiative cooling is included */
  int4byte num_files;      /*!< determines the number of files that are used for a snapshot */
  double BoxSize;          /*!< Simulation box size (in code units) */
  double Omega0;           /*!< matter density */
  double OmegaLambda;      /*!< vacuum energy density */
  double HubbleParam;      /*!< little 'h' */
  int4byte flag_stellarage;     /*!< flags whether the age of newly formed stars is recorded and saved */
  int4byte flag_metals;         /*!< flags whether metal enrichment is included */
  int4byte hashtabsize;         /*!< gives the size of the hashtable belonging to this snapshot file */
  char fill[84];		/*!< fills to 256 Bytes */
}
header, header1;



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef TWO_ARRAY_SIZES
extern int DMSIZE, BARSIZE;
#endif
extern int SIZE;
extern double BOXSIZE;
extern char GlassFile[];
extern int GlassTileFac;
extern double TAVG, NEAVG;
extern double AINIT, Mpart[20];
extern double hconst, OmegaM, OmegaB;
extern char BASEOUT[], OutputDir[];

int find_files(char *fname);


void read_glass(double ***dis, int species)
{
  int i, j, k, n, m, count, type;
  unsigned int dummy, dummy2;
  float *pos = 0;
  float x, y, z;
  FILE *fd = NULL;
  size_t bytes;
  int num, numfiles, skip, nlocal, Nglass;
  char buf[500], fname[200];

  int NumPart = SIZE*SIZE*SIZE;
  double periodic_wrap(double x);


#define SKIP {fread(&dummy, sizeof(int), 1, fd);}
#define SKIP2 {fread(&dummy2, sizeof(int), 1, fd);}

  sprintf(fname, "%s", GlassFile);

  if(1)
    {
      printf("\nreading Lagrangian glass file...\n");
      fflush(stdout);

      numfiles = find_files(fname);

      for(num = 0, skip = 0; num < numfiles; num++)
	{
	  if(numfiles > 1)
	    sprintf(buf, "%s.%d", fname, num);
	  else
	    sprintf(buf, "%s", fname);

	  if(!(fd = fopen(buf, "r")))
	    {
	      fprintf(stderr, "can't open file `%s' for reading glass file.\n", buf);
	      exit(-1);
	    }

	  SKIP;
	  fread(&header1, sizeof(header1), 1, fd);
	  SKIP2;

	  if(dummy != sizeof(header1) || dummy2 != sizeof(header1))
	    {
	      fprintf(stderr, "incorrect header size!\n");
	      exit(-2);
	    }

	  nlocal = 0;

	  for(k = 0; k < 6; k++)
	    nlocal += header1.npart[k];


	  printf("reading '%s' with %d particles\n", fname, nlocal);

	  if(num == 0)
	    {
	      Nglass = 0;

	      for(k = 0; k < 6; k++)
		Nglass += header1.npartTotal[k];

	      printf("\nNglass= %d\n\n", Nglass);
	      pos = (float *) malloc(sizeof(float) * Nglass * 3);

	      if(!(pos))
		{
		  fprintf(stderr, "failed to allocate %g Mbyte for glass file\n",
			 sizeof(float) * Nglass * 3.0 / (1024.0 * 1024.0));
		  exit(-112);
		}
	    }

	  SKIP;
	  fread(&pos[3 * skip], sizeof(float), 3 * nlocal, fd);
	  SKIP2;

	  // for(n=0; n<Nglass;n++)
	  //  printf("%e\n", pos[3*n]);

	  if(dummy != sizeof(float) * 3 * nlocal || dummy2 != sizeof(float) * 3 * nlocal)
	    {
	      printf("incorrect block structure in positions block!\n");
	      exit(-3);
	    }
	  skip += nlocal;

	  fclose(fd);
	}
    }

  count = 0;
  
  for(i = 0; i < GlassTileFac; i++)
    for(j = 0; j < GlassTileFac; j++)
      for(k = 0; k < GlassTileFac; k++)
	{
	  for(type = 0, n = 0; type < 6; type++)
	    {
	      for(m = 0; m < header1.npartTotal[type]; m++, n++)
		{

		  x = pos[3 * n] / header1.BoxSize * (BOXSIZE / GlassTileFac) + i * (BOXSIZE / GlassTileFac);

		  y = pos[3 * n + 1] / header1.BoxSize * (BOXSIZE / GlassTileFac) + j * (BOXSIZE / GlassTileFac);
		  z = pos[3 * n + 2] / header1.BoxSize * (BOXSIZE / GlassTileFac) + k * (BOXSIZE / GlassTileFac);

		  if(species == 1) //shift baryonic glass
		    {
#ifdef SMALL_SHIFT  //shifts the position of the baryons by half a mean interparticle spacing in each direction
		      x += .5*BOXSIZE/SIZE; //move by a cell size
		      y += .5*BOXSIZE/SIZE;
		      z += .5*BOXSIZE/SIZE;
#else
		      x += .5*BOXSIZE/GlassTileFac;  //move by half a tiled glass box in each direction
		      y += .5*BOXSIZE/GlassTileFac;  //this choice is arbitrary 
		      z += .5*BOXSIZE/GlassTileFac;	
#endif
		      x = periodic_wrap(x); y = periodic_wrap(y); z = periodic_wrap(z); //enforce periodicity
		    }

		  dis[species][0][count] = (double) x; 
		  dis[species][1][count] = (double) y;
		  dis[species][2][count] = (double) z;

		  count++;
		}
	    }
	}

  if(count != NumPart)
    {
      printf("fatal mismatch (%d %d)\n", count, NumPart);
      exit(-1);
    }

  free(pos); 
}




#ifdef OUTPUT_GADGET

const double MPCHTOGADGETUNITS = 1000;

void save_gadget_snap(double ***dis, double ***vel, double *delta, 
		      double *temp, int Ndm, int Ngas, int count)
{
#define BUFFER 10
  size_t bytes;
  float *block;
  int *blockid;
  int blockmaxlen,  maxlongidlen;
  int4byte dummy;
  FILE *fd;
  char buf[300];
  int i, k, pc, numf;
  //int NumPart = SIZE*SIZE*SIZE;
  int TotNumPartGAS, TotNumPartDM;
#ifdef TWO_ARRAY_SIZES
  TotNumPartGAS = BARSIZE*BARSIZE*BARSIZE;
  TotNumPartDM = DMSIZE*DMSIZE*DMSIZE;
#else
  TotNumPartGAS = TotNumPartDM = SIZE*SIZE*SIZE;
#endif

  double rhoavg, sqainit = sqrt(AINIT);
  double Box = BOXSIZE*MPCHTOGADGETUNITS;

  if(Ngas+Ndm == 0)
    return;

  if(Ngas == Ndm)
    {
      numf = 1;
      sprintf(buf, "%s/%s_IC.dat", OutputDir, BASEOUT);
    }else
    {
      numf = 2;
      sprintf(buf, "%s/%s_IC.dat.%d", OutputDir, BASEOUT, count);
    }
  if(!(fd = fopen(buf, "w")))
    {
      printf("Error. Can't write in file '%s'\n", buf);
      exit(-10);
    }
  fprintf(stderr, "Opened '%s' to output IC.\n", buf);


  for(i = 0; i < 6; i++)
    {
      header.npart[i] = 0;
      header.npartTotal[i] = 0;
      header.mass[i] = 0;
    }


  double UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  double UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
  double G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  double Hubble = HUBBLE * UnitTime_in_s;

  header.npart[1] = Ndm;
  header.npartTotal[1] = TotNumPartDM;
  header.npartTotal[2] = (TotNumPartDM >> 32); //I have no idea what this does

  header.npart[0] = Ngas;
  header.npartTotal[0] = TotNumPartGAS;
  header.mass[0] = (OmegaB) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / TotNumPartGAS;
  header.mass[1] = (OmegaM - OmegaB) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / TotNumPartDM;

  //compare with what I calculated elsewhere
  //fprintf(stdout, "%le %le %le %le\n", header.mass[0], header.mass[1], Mpart[1], Mpart[0]); exit(-10);


  header.time = AINIT;
  header.redshift = 1.0 / AINIT - 1;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 1;
  header.flag_stellarage = 0;
  header.flag_metals = 0;

  header.num_files = numf;

  header.BoxSize = Box;
  header.Omega0 = OmegaM;
  header.OmegaLambda = 1-OmegaM;
  header.HubbleParam = hconst;

  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.hashtabsize = 0;
  memset(header.fill, 0, 84);

  dummy = sizeof(header);
  fwrite(&dummy, sizeof(dummy), 1, fd);
  fwrite(&header, sizeof(header), 1, fd);
  fwrite(&dummy, sizeof(dummy), 1, fd);


  // meanspacing = Box / pow(TotNumPart, 1.0 / 3);


  if(!(block = (float *) malloc(bytes = BUFFER * 1024 * 1024)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double)bytes);
      exit(-24);
    }

  blockmaxlen = bytes / (3 * sizeof(float));

  blockid = (int *) block;
  maxlongidlen = bytes / (sizeof(long long));

  /* write coordinates */
  dummy = sizeof(float) * 3 * (Ndm +Ngas);

  //gas positions
  fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < Ngas; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  block[3 * pc + k] = (float) dis[1][k][i]*MPCHTOGADGETUNITS;
	}

      pc++;

      if(pc == blockmaxlen)
	{
	  fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    fwrite(block, sizeof(float), 3 * pc, fd);


  //dark matter positions
  for(i = 0, pc = 0; i < Ndm; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  block[3 * pc + k] = (float) dis[0][k][i]*MPCHTOGADGETUNITS;
	}

      pc++;

      if(pc == blockmaxlen)
	{
	  fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    fwrite(block, sizeof(float), 3 * pc, fd);

  fwrite(&dummy, sizeof(dummy), 1, fd);


  /* write velocities */
  dummy = sizeof(float) * 3 * (Ndm+Ngas);

  //  dummy *= 2;//with gas

  fwrite(&dummy, sizeof(dummy), 1, fd);

  for(i = 0, pc = 0; i < Ngas; i++)
    {
      for(k = 0; k < 3; k++)
	block[3 * pc + k] = (float) vel[1][k][i]/sqainit;

      pc++;

      if(pc == blockmaxlen)
	{
	  fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    fwrite(block, sizeof(float), 3 * pc, fd);


  //write gas
  for(i = 0, pc = 0; i < Ndm; i++)
    {
      for(k = 0; k < 3; k++)
	block[3 * pc + k] = (float) vel[0][k][i]/sqainit;

      pc++;

      if(pc == blockmaxlen)
	{
	  fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    fwrite(block, sizeof(float), 3 * pc, fd);

  fwrite(&dummy, sizeof(dummy), 1, fd);


  /* write particle ID */
  dummy = sizeof(int) * (Ndm+Ngas);

  // dummy *= 2;//there is gas

  fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < Ngas; i++)
    {
      blockid[pc] = i+1;

      pc++;

      if(pc == maxlongidlen)
	{
	  fwrite(blockid, sizeof(int), pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    {
      fwrite(blockid, sizeof(int), pc, fd);
    }

  for(i = 0, pc = 0; i < Ndm; i++)
    {
      blockid[pc] = i + TotNumPartGAS +1;

      pc++;

      if(pc == maxlongidlen)
	{
	  fwrite(blockid, sizeof(int), pc, fd);

	  pc = 0;
	}
    }
  if(pc > 0)
    {
      fwrite(blockid, sizeof(int), pc, fd);
    }

  fwrite(&dummy, sizeof(dummy), 1, fd);


  /* write temperatures if needed */
  dummy = sizeof(float) * Ngas;
  fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < Ngas; i++)
    {
#ifndef NO_TEMPERATURE
      block[pc] = TAVG;
      //need to CIC temperature to put this in here
      block[pc] *= (1+temp[i]);
      // printf("Havent' checked temp--gadget units?\n");
      // exit(-5);
#else
      block[pc] = 0.;
#endif
      //convert to internal energy
      // MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * PROTONMASS;	 
	
      block[pc] /= (1.22*PROTONMASS/BOLTZMANN * (5./3.-1)); //gamma=5/3 and primordial composition
      block[pc] /=  (UnitEnergy_in_cgs/ UnitMass_in_g);

      pc++;

      if(pc == blockmaxlen)
	{
	  fwrite(block, sizeof(float), pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    fwrite(block, sizeof(float), pc, fd);
  fwrite(&dummy, sizeof(dummy), 1, fd);

#ifndef NO_DELTA //don't ned to calculate these if set temperature in gadget
  /* write densities */
  dummy = sizeof(float) * Ngas;
  fwrite(&dummy, sizeof(dummy), 1, fd);
  rhoavg = (OmegaB) * 3 * Hubble * Hubble / (8 * PI * G);
  for(i = 0, pc = 0; i < Ngas; i++)
    {
      block[pc] = rhoavg*(1+delta[i]); 

      pc++;

      if(pc == blockmaxlen)
	{
	  fwrite(block, sizeof(float), pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    fwrite(block, sizeof(float), pc, fd);
  fwrite(&dummy, sizeof(dummy), 1, fd);

  /* write Ne --needed to go from U to temperature*/
  dummy = sizeof(float) * Ngas;
  fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < Ngas; i++)
    {
      block[pc] = NEAVG;

      pc++;

      if(pc == blockmaxlen)
	{
	  fwrite(block, sizeof(float), pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    fwrite(block, sizeof(float), pc, fd);
  fwrite(&dummy, sizeof(dummy), 1, fd);
#endif


  free(block);

  fclose(fd);
}

#endif


/*****************************************************************************************
routine used to find file names of gadget snapshots....ripped from GADGET2 by Volker Springel
*****************************************************************************************/
int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;

  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if((fd = fopen(buf, "r")))
    {
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header, sizeof(header), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      fclose(fd);

      return header.num_files;
    }

  if((fd = fopen(buf1, "r")))
    {
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header, sizeof(header), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      fclose(fd);
      header.num_files = 1;

      return header.num_files;
    }

  exit(-121);
  return 0;
}
