

#ifdef OUTPUT_ENZO 

#define FLOAT double //float
#define EEINT long int //int
#define FLOAT_UNDEFINED -99999.
#define INT_UNDEFINED -99999
#define debug 1
#define SUCCESS 1
#define LINUX 1

#define CRITDENSITY 1.8791e-29    //cgs
#define GNEWTON 6.67259e-8       //cgs  -- old value
#define MPCTOKM 3.08568e19
#define MPC 3.08568e24
#define M_H 1.67353258e-24// 1.67//262177e-24 //cgs
#define K_B 1.38065e-16  //cgs
#define MU 1.214 //0.6 //from enzo - should be 1.22
#define GAMMA 1.6667

#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "enzofiles/enzo_param.h"
#include "enzofiles/extern_hdf5.h"
 
int WriteField(EEINT Rank, EEINT Dims[3], FLOAT *Field, char *Name, 
	       EEINT Part, EEINT Npart,
               EEINT GridRank, EEINT Starts[3], EEINT Ends[3], EEINT RootGridDims[3]);
int SetParameterDefaults(parmstruct *Parameters);
void fcol(FLOAT *x, int n, int m, FILE *log_fptr);


double MPART_ENZO;
double DensityUnits, LengthUnits, TimeUnits, VelocityUnits,TempUnits;
 
int         io_log = 1;
int         dump_ok = 0;


extern double TAVG;
extern double BOXSIZE;
extern double OmegaM, OmegaB, hconst;
extern int SIZE;
extern double Mpart[];
extern double ZINIT;

void setEnzoUnits()
{
  DensityUnits = CRITDENSITY*OmegaM*hconst*hconst*pow(1.+ZINIT, 3.);
  LengthUnits = BOXSIZE*MPC/(hconst)/(1.+ZINIT);  //something is wrong with this, need to divide by boxsize in kpc?  should be between 0. and 1.
  TimeUnits = 1./sqrt(4.*M_PI*GNEWTON*DensityUnits);
  VelocityUnits = (LengthUnits/TimeUnits); //1+z cancel out


  MPART_ENZO = DensityUnits*pow(LengthUnits, 3.)/(Mpart[0]*1e10/hconst);//dark matter

  fprintf(stdout, "Enzo Units: density %le Length %le Time %le Velocity %le\n",
	 DensityUnits, LengthUnits, TimeUnits, VelocityUnits);
  fprintf(stdout, "Enzo unit change the length by 100.*1_z and time by 1+z\n");
  TempUnits = M_H*VelocityUnits*VelocityUnits/K_B;  //missing MU just like in ENZO
  
  fprintf(stdout, "New Enzo Units: density %le Length %le Time %le Velocity %le Temp %le\n",
	  DensityUnits, LengthUnits, TimeUnits, VelocityUnits, TempUnits);

  double OD,OL,OT,OV,OTEMP;
  OD = 1.8788e-29*OmegaM*POW(hconst,2.)*POW(1. + ZINIT,3.);
  OL = 3.085678e24*BOXSIZE/hconst/
                      (1. + ZINIT);
  OTEMP =   1.88e6*POW(BOXSIZE,2.)*OmegaM*
                      (1. +  ZINIT);

  OT = 2.519445e17/sqrt(OmegaM)/hconst/
                      POW(1. +  ZINIT,FLOAT(1.5));
  
  OV =  1.22475e7*BOXSIZE*sqrt(OmegaM)*
                      sqrt(1. +  ZINIT);


  fprintf(stdout, "using OT instead\n");
  TempUnits=OTEMP;
  fprintf(stdout,"EnzoDeft Units density %le Length %le Time %le Velocity %le Temp %le\n", OD,OL,OT,OV,OTEMP );
}


 

int SaveEnzoSnapshot(double **disdm, double **veldm, double *deltab, double **velb, double *temp)
{
  parmstruct *Parameters;
  FLOAT *ParticleField, *GridField;
  // int *ParticleTypeField;
 
  //FLOAT GrowthFunction, aye = 1.0, ayed, Temp;
  int i, dim, size=1;
  int j,k,f1,f2;
  int *FI; //fortran order
  EEINT NumberOfParticles=1;
  FILE *dumpfile;
  double x;
  //FI is used to flip the arrays into fortran format from c.
  FI = (int *) malloc(sizeof(int)*SIZE*SIZE*SIZE);
  for(k=0; k<SIZE; k++){
    for(j=0; j<SIZE; j++){
      for (i = 0; i <SIZE; i++){
	f1 = i+j*SIZE+k*SIZE*SIZE;
	f2 = k+j*SIZE+i*SIZE*SIZE;
	FI[f1]=f2;
      }
    }
  }
  


  Parameters = (parmstruct *) malloc(sizeof(struct parmstruct));

  SetParameterDefaults(Parameters);

  setEnzoUnits();

  
  for (dim = 0; dim < Parameters->Rank; dim++) {
      size *= (Parameters->ParticleDims[dim] + ((dim == 0) ? 2 : 0));
      NumberOfParticles *= Parameters->ParticleDims[dim];
        fprintf(stderr, "PARTICLEDIMS %"ISYM" are %"ISYM"\n", dim,
		Parameters->ParticleDims[dim]);
  }

  //if((ParticleField = (float *) malloc(sizeof(FLOAT)*NumberOfParticles)) == NULL){
  if ((ParticleField = new FLOAT[size]) == NULL) {
    fprintf(stderr, "GenerateRealization: malloc failure (%d).\n", size);
    exit(EXIT_FAILURE);
  }
  

  if (debug) printf("NumberOfParticles = %ld\n", NumberOfParticles);
 

  // Loop over dimensions
  //need to use fortran ordering for grid.  fine for particles!
  for (dim = 0; dim < 3; dim++) 
    {

      double x = (1.+ZINIT)*LengthUnits/(MPC/hconst); 
      for(i=0; i< SIZE*SIZE*SIZE; i++)	
	ParticleField[i] = disdm[dim][i]/x; //convert to code units
      

      
      printf("%ld %ld %ld %ld\n", NumberOfParticles, Parameters->Rank, Parameters->TopGridStart[0], Parameters->RootGridDims[0]);
      //   NumberOfParticles=10;
      fprintf(stderr,"#writefield1\n");
      WriteField(1, &NumberOfParticles,
	       ParticleField, Parameters->ParticlePositionName, dim, 3,
	       Parameters->Rank,
	       Parameters->TopGridStart,
	       Parameters->TopGridEnd,
	       Parameters->RootGridDims);

    x = (VelocityUnits/1e5);
    fprintf(stderr,"#pf1 %g %g\n", ParticleField[10],ParticleField[100]);

    for(i=0; i< SIZE*SIZE*SIZE; i++)  //was in units of km/s
      ParticleField[i] = veldm[dim][i]/x;
    
    fprintf(stderr,"#writefield2\n");
  fprintf(stderr,"#pf1 %g %g\n", ParticleField[10],ParticleField[100]);
    
    WriteField(1, &NumberOfParticles,
	       ParticleField, Parameters->ParticleVelocityName, dim, 3,
	       Parameters->Rank,
	       Parameters->TopGridStart,
	       Parameters->TopGridEnd,
	       Parameters->RootGridDims);
 
  } // end: loop over dims
 
    // Generate and write out particle mass field (dark matter particles)

#ifdef OUTPUT_DMMASS  
  for(i=0; i< SIZE*SIZE*SIZE; i++) 
    ParticleField[i] = MPART_ENZO; //dark matter particles
    
    WriteField(1, &NumberOfParticles,
                 ParticleField,
		 Parameters->ParticleMassName, 0, 1,
                 Parameters->Rank,
                 Parameters->TopGridStart,
                 Parameters->TopGridEnd,
                 Parameters->RootGridDims);
#endif
 
 
    delete ParticleField;
 

 
  // end: if (InitializeParticles)
 
  /* ------------------------------------------------------------------- */
 
  // Set grids
 
    // Compute required size and allocate field
    size = 1;
    for (dim = 0; dim < Parameters->Rank; dim++)
      size *= (Parameters->GridDims[dim] + ((dim == 0) ? 2 : 0));
 
 
    if ((GridField = new FLOAT[size]) == NULL) {
      fprintf(stderr, "GenerateRealization: malloc failure (%d).\n", size);
      exit(EXIT_FAILURE);
    }
 
    /* 1) density (add one and multiply by mean density). */
 
    for(i=0; i< SIZE*SIZE*SIZE; i++)
      GridField[FI[i]] = OmegaB/OmegaM*(1+deltab[i]); //enzo density units
    
    fprintf(stderr,"#db %g %g\n", deltab[10],deltab[100]);

    WriteField(Parameters->Rank, Parameters->GridDims,
	       GridField, Parameters->GridDensityName, 0, 1,
               Parameters->Rank, Parameters->TopGridStart,
               Parameters->TopGridEnd,
               Parameters->RootGridDims);
 
    if (dump_ok) fprintf(dumpfile, "Density, size=%"ISYM"\n", size);
    if (dump_ok) fcol(GridField, size, 8, dumpfile);
    fprintf(stderr,"#griddensity\n");
    /* 2) velocities. */
 
    for (dim = 0; dim < Parameters->Rank; dim++) 
      {
 
	/* 1) velocities
	   -generate the displacement field (f delta_k -i vec(k)/k^2).
	   -multiply adot to get velocity */
	
	x=(VelocityUnits/1e5);
	for(i=0; i< SIZE*SIZE*SIZE; i++)
	  GridField[FI[i]] = velb[dim][i]/x;
	
	fprintf(stderr,"#grid velcoity %ld\n",dim);
	fprintf(stderr,"#gf1 %g %g\n", GridField[10],GridField[100]);

	WriteField(Parameters->Rank, Parameters->GridDims,
		   GridField, Parameters->GridVelocityName, dim,3,
		   Parameters->Rank, Parameters->TopGridStart,
		   Parameters->TopGridEnd,
		   Parameters->RootGridDims);
      }

    //#ifndef NO_TEMPERATURE 
    for(i=0; i< SIZE*SIZE*SIZE; i++)
      GridField[FI[i]] = ((temp[i]+1.)*TAVG)/(GAMMA-1.)/(MU)/TempUnits;
    fprintf(stdout,"#Temperature Average %le : in units %le %le \n",TAVG,(TAVG/(GAMMA-1.)/(MU)/TempUnits),TempUnits);
        WriteField(Parameters->Rank, Parameters->GridDims,
		   GridField, Parameters->GridGasEnergyName, 0, 1,
               Parameters->Rank, Parameters->TopGridStart,
               Parameters->TopGridEnd,
               Parameters->RootGridDims);
	//#endif
    
	fprintf(stderr,"#finished enzo_out.c\n");

  return 0.;
}
 
int SetParameterDefaults(parmstruct *Parameters)
{ 

  int dim;
 
  char ppos_name[] = "ParticlePositions";
  char pvel_name[] = "ParticleVelocities";
  char pmass_name[] = "ParticleMasses";
  char gden_name[] = "GridDensity";
  char gvel_name[] = "GridVelocities";
  char gene_name[] = "GridEnergy";
 
  // Set Defaults
 
  Parameters->Rank                = 3;
  Parameters->GridRefinement      = 1;
  Parameters->ParticleRefinement  = 1;
  Parameters->WaveNumberCutoff    = INT_UNDEFINED;
  Parameters->InitializeParticles = TRUE;
  Parameters->InitializeGrids     = TRUE;
  //Parameters->RandomNumberGenerator = 0;
 
  strcpy(Parameters->ParticlePositionName, ppos_name); 
  strcpy(Parameters->ParticleVelocityName, pvel_name);
  strcpy(Parameters->ParticleMassName, pmass_name);
  //Parameters->ParticleTypeName     = NULL;
  strcpy(Parameters->GridDensityName, gden_name);
  strcpy(Parameters->GridVelocityName, gvel_name);
  strcpy(Parameters->GridGasEnergyName, gene_name);
 
  for (dim = 0; dim < Parameters->Rank; dim++) {
    Parameters->GridDims[dim]     = SIZE;
    Parameters->ParticleDims[dim] = SIZE;
    Parameters->MaxDims[dim]      = INT_UNDEFINED;
    Parameters->NewCenter[dim]    = INT_UNDEFINED;
    Parameters->NewCenterFloat[dim] = FLOAT_UNDEFINED;
    Parameters->StartIndex[dim]   = 0;
    Parameters->StartIndexInNewCenterTopGridSystem[dim] = INT_UNDEFINED;
    Parameters->EndIndexInNewCenterTopGridSystem[dim] = INT_UNDEFINED;
    Parameters->TopGridStart[dim] = INT_UNDEFINED;
    Parameters->TopGridEnd[dim] = INT_UNDEFINED;
    Parameters->RootGridDims[dim] = INT_UNDEFINED;
    Parameters->RefineRegionLeftEdge[dim] = FLOAT_UNDEFINED;
    Parameters->RefineRegionRightEdge[dim] = FLOAT_UNDEFINED;
  }

  //Parameters->RefineBy = 2;
  //Parameters->MaximumInitialRefinementLevel = INT_UNDEFINED;
 
  return SUCCESS;
}




int WriteField(EEINT Rank, EEINT Dims[3], FLOAT *Field, char *Name, EEINT Part, EEINT Npart,
               EEINT GridRank, EEINT Starts[3], EEINT Ends[3], EEINT Tops[3])
{
  

  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, mem_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id;
  hid_t       int_file_type_id, int_mem_type_id;

  hsize_t     out_dims[3];
  hsize_t     dimm;
  hsize_t     slab_dims[4];
  hsize_t     slab_rank;
  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[4], file_count[4], file_block[4];
  hsize_t     attr_count;

  
  hssize_t    mem_offset;
  hssize_t    file_offset[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  int         dim;

  EEINT         component_rank_attr;
  EEINT         component_size_attr;
  EEINT         field_rank_attr;
  EEINT         field_dims_attr[3];

  FILE        *dumpfile;
  FILE        *log;

 fprintf(stdout, "  Name %s\n", Name);

  if (dump_ok) dumpfile = fopen("DumpWF","a");
  if (io_log) log = fopen("IO_Log","a");

  if (io_log) fprintf(log, "On entry to WriteField\n");
  if (io_log) fprintf(log, "  Rank %"ISYM"\n", Rank);
  if (io_log) fprintf(log, "  Dims %"ISYM"  %"ISYM"  %"ISYM"\n", Dims[0], Dims[1], Dims[2]);
  if (io_log) fprintf(log, "  Name %s\n", Name);
  if (io_log) fprintf(log, "  Part %"ISYM" of %"ISYM"\n", Part, Npart);

//  GB: Reverse dim ordering since we are using fortran array ordering
//      (actually, we don't have to do this here).
//
//  out_dims[Rank-dim-1] = Dims[dim];

  for ( dim = 0; dim < Rank; dim++ )
  {
    out_dims[dim] = Dims[dim];
  }

//  reverse these also?

  for ( dim =0; dim < GridRank; dim++ )
  {
    if ( Starts[dim] == INT_UNDEFINED )
      Starts[dim] = 0;
    if ( Ends[dim] == INT_UNDEFINED )
      Ends[dim] = Tops[dim] - 1;
  }

  slab_rank = Rank+1;

  slab_dims[0] = Npart;

  for ( dim = 1; dim < slab_rank; dim++ )
  {
    slab_dims[dim] = out_dims[dim-1];
  }

  if (io_log) fprintf(log, "  Extended Rank %"ISYM"\n", (int) slab_rank);

  for ( dim = 0; dim < slab_rank; dim++ )
  {
    if (io_log) fprintf(log, "    %"ISYM":  %"ISYM"\n", dim, (int) slab_dims[dim]);
  }

  dimm = 1;

  for ( dim = 0; dim < Rank; dim++ )
  {
    dimm = dimm * Dims[dim];
  }

  if (io_log) fprintf(log, "  Grid Elements %"ISYM"\n", (int) dimm);

  if (dump_ok) fcol(Field, (int) dimm, 8, dumpfile);

  component_rank_attr = Npart;
  component_size_attr = dimm;

  field_rank_attr = Rank;

  for ( dim = 0; dim < Rank; dim++ )
  {
    field_dims_attr[dim] = Dims[dim];
  }

  /* HDF5
   
     The HDF4 DFSD interface does not have a direct analogue in HDF5.
     In particular, there is no "append" on open and there are no
     specific features to store rank, dimensions, scales etc.
     These are replaced by dataset attributes, which can be named.
     To eliminate redundancy and to keep the data structure close
     to the original, each ENZO file will contain a single dataset
     of the same name (as opposed to separate datasets for each 
     component).

     The dataspace uses the number of components for each of the
     the dimensions of the dataset as an additional "dimension".
     For example, a 3D scalar field of 4x4x4 points will have a 
     dataspace of {1,4,4,4}, while a 3D vector field of, say,
     {Vx,Vy,Vz} and 4x4x4 points will have a dataspace of {3,4,4,4}.
     
     Slab I/O is used in preparation for parallel I/O.

     It is ASSUMED that this routine is called with
     Part = {0,1,2,...,Npart-1} to write Npart components of equal length
     (the product of the field dimensions {Dims[i] for i=0,Rank-1}.
     Each component is offset in the file by {0,1,2,...,Npart-1} * the field size.

     The rank and dimensions of the field and the number of field
     components are HDF5 attributes.

  */

  int ll = sizeof(EEINT);
  printf("ll  %ld  changed to? 8\n",ll);
  printf("ll:  int %ld long %d long long %d   hint %d\n",sizeof(int),sizeof(long),sizeof(long long),sizeof(dimm));
  //ll = 8;
  switch(ll)
  {
 
    case 4:
      int_mem_type_id = HDF5_I4;
      int_file_type_id = HDF5_FILE_I4;
      break;
    case 8:
      int_mem_type_id = HDF5_I8;
      int_file_type_id = HDF5_FILE_I8;
      break;
    default:
      int_mem_type_id = HDF5_I8;
      int_file_type_id = HDF5_FILE_I8;
  }
 
  int ii = sizeof(FLOAT);
  
  switch(ii)
  {

    case 4:
      mem_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
      break;

    case 8:
      mem_type_id = HDF5_R8;
      file_type_id = HDF5_FILE_R8;
      break;

    default:
      mem_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;

  }

// Data in memory is considered 1D, stride 1, with zero offset

  mem_stride = 1;      // contiguous elements
  mem_count = dimm;    // number of elements in field
  mem_offset = 0;      // zero offset in buffer
  mem_block = 1;       // single element blocks

// 1D memory model

  mem_dsp_id = H5Screate_simple(1, &dimm, NULL);
    if (io_log) fprintf(log, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );

  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect mem slab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );


//  Data in the file is (1+Rank)D with Npart components per grid point.
//  Offset[0] is the component Part of Npart components.  Data for each
//  Part are contiguous in the file, so stride = 1.

  file_stride[0] = 1;      // contiguous elements
  file_count[0] = 1;       // one component per call
  file_offset[0] = Part;   // component Part of Npart
  file_block[0] = 1;       // single element blocks

  for ( dim = 1; dim < slab_rank; dim++ )
  {
    file_stride[dim] = 1;                   // contiguous elements
    file_count[dim] = out_dims[dim-1];       // field dimensions
    file_offset[dim] = 0;                   // complete field, no offset
    file_block[dim] = 1;                    // single element blocks
  }

  file_dsp_id = H5Screate_simple(slab_rank, slab_dims, NULL);
    if (io_log) fprintf(log, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
    assert( file_dsp_id != h5_error );

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
    if (io_log) fprintf(log, "H5Sselect file slab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

//  If Part is zero, create an HDF5 file with a single dataset of 
//  the same name and attach the dataset attributes, otherwise
//  open an existing file and dataset. */

  if ( Part == 0 )
  {
    if (io_log) fprintf(log, "Calling H5Fcreate with Name = %s\n", Name);

    file_id = H5Fcreate(Name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fcreate id: %"ISYM"\n", file_id);
      assert( file_id != h5_error );

    if (io_log) fprintf(log, "Calling H5Dcreate with Name = %s\n", Name);
    if (io_log) fprintf(log, "%"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",file_id,file_type_id,file_dsp_id,H5P_DEFAULT);

    //   fprintf(stdout, "%d %s %d %d %d\n", 
    dset_id =  H5Dcreate(file_id, Name, file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Dcreate id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );


    attr_count = 1;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Component_Rank",  int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  int_mem_type_id, &component_rank_attr);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );


    attr_count = 1;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Component_Size",  int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  int_mem_type_id, &component_size_attr);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );


    attr_count = 1;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Rank", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  int_mem_type_id, &field_rank_attr);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );


    attr_count = Rank;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "Dimensions", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  int_mem_type_id, field_dims_attr);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );


    attr_count = GridRank;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "TopGridStart", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  int_mem_type_id, Starts);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    attr_count = GridRank;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "TopGridEnd", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  int_mem_type_id, Ends);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    attr_count = GridRank;

    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );

    attr_id = H5Acreate(dset_id, "TopGridDims", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );

    h5_status = H5Awrite(attr_id,  int_mem_type_id, Tops);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );

  }

  else

  {
    if (io_log) fprintf(log, "Calling H5Fopen with Name = %s\n", Name);

    file_id = H5Fopen(Name, H5F_ACC_RDWR, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fopen id: %"ISYM"\n", file_id);
      assert( file_id != h5_error );

    if (io_log) fprintf(log, "Calling H5Dopen with Name = %s\n", Name);

    dset_id =  H5Dopen(file_id, Name);
      if (io_log) fprintf(log, "H5Dopen id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );

  }


  h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, Field);
    if (io_log) fprintf(log, "H5Dwrite: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  if (io_log) fprintf(log, "Exit WriteField\n");

  if (dump_ok) fclose(dumpfile);
  if (io_log) fclose(log);

  return SUCCESS;
}

void fcol(FLOAT *x, int n, int m, FILE *log_fptr)
{
 
  int nrow,mrow;
  int i,j;
 
  nrow = n/m;
  mrow = n - nrow * m;
 
  if( mrow > 0 )
  {
    nrow = nrow+1;
  }
 
  fprintf(log_fptr, "\n");
 
  for(j=0;j<n;j=j+m)
  {
    for(i=j;i<min(j+m,n);i++)
    {
      fprintf(log_fptr, "%12.4e", x[i]);
    }
    fprintf(log_fptr, "\n");
  }
 
}

#endif
