CC=g++
CFLAGS=-Wall -g -O3 -I/opt/local/include -I/opt/local/include/gsl  -DRECFAST -DPRINT_PK -DGSL #DONT_PRINT_ICS
LIB=-lm -L/opt/local/lib -lgsl -lgslcblas


RM=/bin/rm -f

EXE = transfer.x
OBJ = main.o #dnumrecipes.o dcomplex.o 

$(EXE): $(OBJ) Makefile
	$(CC) $(CFLAGS) -o $(EXE) $(OBJ) $(LIB)

.cc.o: MAKEFILE
	$(CC) $(CFLAGS) -o $*.o -c $*.cc

clean:
	$(RM) $(EXE) $(OBJ)

