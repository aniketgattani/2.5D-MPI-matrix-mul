EXEC=reduce_avg

OBJ =  $(EXEC) $(EXEC)-debug $(EXEC)-trace

VIEWER=jumpshot

# flags
OPT=-O2 -g
DEBUG=-O0 -g

all: $(OBJ)

# build the debug parallel version of the program
$(EXEC)-debug: $(EXEC).cpp
	mpicxx $(DEBUG) $(OMP) -o $(EXEC)-debug $(EXEC).cpp -lrt


# build the version of the program that records MPE traces for jumpshot
$(EXEC)-trace: $(EXEC).cpp
	mpecxx -mpilog $(OPT) -o $(EXEC)-trace $(EXEC).cpp -lrt 

# build the optimized parallel version of the program
$(EXEC): $(EXEC).cpp
	mpicxx $(OPT) $(OMP) -o $(EXEC) $(EXEC).cpp -lrt

runp:
	echo try running like this when on an interactive compute node:
	echo "srun -n <number of processors> yourprogram <arguments>"

#view a trace with jumpshot 
view:
	$(VIEWER) Unknown.clog2

clean:
	/bin/rm -rf $(OBJ) 
