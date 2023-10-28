/***********************************************************************
/
/  STRUCTURE FOR PARAMETERS
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

struct parmstruct {

  /* Size parameters */

  EEINT Rank;
  EEINT GridDims[3];
  EEINT ParticleDims[3];
  EEINT MaxDims[3];
  EEINT NewCenter[3];
  EEINT StartIndex[3];
  EEINT GridRefinement;
  EEINT ParticleRefinement;

  /* Temporary parameters (used to set other parameters, in
     ReadParameterFile and then not used after). */

  FLOAT NewCenterFloat[3];
  EEINT StartIndexInNewCenterTopGridSystem[3];
  EEINT EndIndexInNewCenterTopGridSystem[3];
  EEINT TopGridStart[3];
  EEINT TopGridEnd[3];
  EEINT RootGridDims[3];

  /* Parameters used to automatically generate hierarchy. */

  FLOAT RefineRegionLeftEdge[3];
  FLOAT RefineRegionRightEdge[3];
  EEINT RefineBy;
  EEINT MaximumInitialRefinementLevel;

  /* Boolean flags. */

  EEINT InitializeParticles;
  EEINT InitializeGrids;
  EEINT RandomNumberGenerator;

  /* Names. */

  char ParticlePositionName[100];
  char ParticleVelocityName[100];
  char ParticleMassName[100];
  char ParticleTypeName[100];
  char GridDensityName[100];
  char GridVelocityName[100];
  char GridGasEnergyName[100];

  /* Power spectrum. */

  EEINT WaveNumberCutoff;

};

/* HDF5 definitions */


#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I8 H5T_STD_I64BE
#define HDF5_FILE_R4 H5T_IEEE_F32BE
#define HDF5_FILE_R8 H5T_IEEE_F64BE
#define HDF5_FILE_B8 H5T_STD_B8BE


/* #ifdef HDF5_LE */
/* #define HDF5_FILE_I4 H5T_STD_I32LE */
/* #define HDF5_FILE_I8 H5T_STD_I64LE */
/* #define HDF5_FILE_R4 H5T_IEEE_F32LE */
/* #define HDF5_FILE_R8 H5T_IEEE_F64LE */
/* #define HDF5_FILE_B8 H5T_STD_B8LE */
/* #endif */

#define HDF5_I4 H5T_NATIVE_INT
#define HDF5_I8 H5T_NATIVE_LLONG
#define HDF5_R4 H5T_NATIVE_FLOAT
#define HDF5_R8 H5T_NATIVE_DOUBLE
#define HDF5_R16 H5T_NATIVE_LDOUBLE

/* Precision-dependent definitions */

#define Eint int
#define ISYM "d"
#define HDF5_INT HDF5_I4
#define HDF5_FILE_INT HDF5_FILE_I4
#define nint(A) ( (int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define ABS(A) abs((int) (A))


#define FALSE 0
#define TRUE 1

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define POW(X,Y) pow((double) (X), (double) (Y))

#ifndef OLD_HDF5
#define hssize_t hsize_t
#endif
