#----------------------------------------------------------------------
# From the list below, please activate/deactivate the options that     
# apply to your run. If you modify any of these options, make sure     
# that you recompile the whole code by typing "make clean; make".      
#
# The format of this Makefile copies the excellent Makefile for GADGET2, by Volker Springel                                                                      
#----------------------------------------------------------------------


#-----------------------------Output operations
OPT += -DOUTPUT_GADGET      	    #OUTPUT a gadget snapshot file
#OPT += -DTWO_GADGET_OUTFILES       #do dark matter and gas calculation separately and output 2 files (to conserve memory) 
#OPT += -DOUTPUT_ENZO               #OUTPUT enzo grid file (cannot output gadget and enzo for same run)


#-----------------------------Glass file options (program uses a single glass file that is tiled appropriately accross the box)
#---------------------------------------------if not set, uses a grid.  It would be unusual to initialize DM with glass but not gas
OPT += -DGLASS_DM                    #Generates DM positions from a glass file
OPT += -DGLASS_GAS                   #Generate gas positions from a glass file (the default is so that they are staggered at .5 mean interparticle spacings from dark matter particles) 
OPT += -DLARGE_SPACING               #Offset gas glass positions by dimension_glass/2 in each direction (This is the preferred method in Yoshida et al (2003)
                                     # but we do not recommend it as we find it leads to larger particle coupling.

#------------------------------Miscellaneous
#OPT += -DNO_DELTA                   #do not calculate/output baryonic density  CANNOT BE SET FOR ENZO OUTPUTS. Saves a minor amount of memory.
#OPT += -DNO_TEMPERATURE             #does not output temperature FOR GADGET ONLY Saves a minor amount of memory.
#OPT += -DNGP                        # The default is CIC gridding for interpolation onto glass locations, but one can do NGP
#OPT += -DBOOST_EVERYTHING_BYVKMSEC=3  #velocity boost in the i direction at start of simulation (tests that results do not depend on frame)
#OPT += -DUSE_OTHER_GRID_SIZE=256    #use an initial conditions grid that is different than the output grid
#OPT += -DTWO_ARRAY_SIZES            #currently this is unfinished

#------------------------------Performance
#OPT += -DOPENMP -fopenmp             #parallelizes over shared memory part of caclulation that is slowest, which is the 2D spline

                   
OPT += -DGSL_INTERP                  #THIS ALWAYS NEEDS TO BE SET SINCE OTHERWISE CODE USES NUMBERICAL RECIPES FUNCTIONS THAT ARE NOT PROVIDED

CC=g++
OPTIMIZE =  -Wall -g -O3 


SYSTYPE="quasar"

ifeq ($(SYSTYPE),"henyey")
GSL_INCL =  -I/global/software/centos-5.x86_64/modules/gsl/1.9-intel/gsl-1.9/include
GSL_LIBS =  -L/global/software/centos-5.x86_64/modules/gsl/1.9-intel/gsl-1.9/lib
FFTW_INCL=  -I/global/software/centos-5.x86_64/modules/fftw/3.2.2-gcc/include
FFTW_LIBS=  -L/global/software/centos-5.x86_64/modules/fftw/3.2.2-gcc/lib
HDF5INCL =  -I/global/home/users/mqk/local/hdf5/1.8.4-intel/include           #only need for ENZO outs
HDF5LIB  =  -L/global/home/users/mqk/local/hdf5/1.8.4-intel/lib -lhdf5 -lz    #only need for ENZO outs, but comment out -l if you don't need hdf5
endif

ifeq ($(SYSTYPE),"quasar")
GSL_INCL =  -I/opt/local/include/gsl
GSL_LIBS =  -L/opt/local/lib/
FFTW_INCL=  -I/opt/local/include/           
FFTW_LIBS=  #-L/opt/local/lib         
HDF5INCL =     			      #hdf5 only needed for ENZO outs
HDF5LIB  =  #-lhdf5 -lz               #hdf5 only needed for ENZO outs
endif

CFLAGS = $(OPTIMIZE) $(OPT) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL)

LIBS =   $(HDF5LIB) $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIB) -lfftw3 -lm -DH5_USE_16_API -DLINUX -DOPTIMIZED_CTP


RM=/bin/rm -f

EXE = genICs.x
OBJ = main.o io.o read_glass.o getPk.o enzo_out.o  #dnumrecipes.o

$(EXE): $(OBJ) Makefile
	$(CC) $(CFLAGS) -o $(EXE) $(OBJ) $(LIBS)

.c.o: MAKEFILE
	$(CC) $(CFLAGS) -o $*.o -c $*.c

clean:
	$(RM) $(EXE) $(OBJ)
