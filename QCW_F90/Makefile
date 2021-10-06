#####################################################################
#
# Typing the command below => results in the following actions
#    make             =>      Executable a.out
#    make clean       =>      Remove object files during compilation
#    make realclean   =>      Remove object files, executables
#
######################################################################

FF = mpiifort
PAL = -D PARALLEL
LIB = -mkl

FILES = mt19937.o FKQCW.o 

all: FKQCW
FKQCW: $(FILES)
	$(FF) $(LIB) $(FILES) -o $@

clean:
	rm -f *.o *.mod

.SUFFIXES:

mt19937.o: mt19937.f90
	$(FF) -c mt19937.f90

FKQCW.o: FKQCW.F90
	$(FF) $(PAL) -c FKQCW.F90
