# Makefile for SHDOMPP in distribution
BIN = bin
SRC = src
OBJ = shdompp  \
      disortsh \
		  ppmieprp \
			put

FFLAGS= -O2
FC= gfortran
GCC=gcc

# Main
main: $(OBJ)

shdompp: $(SRC)/shdompp.f90 $(SRC)/fftpack.o
	$(FC) $(FFLAGS) $(SRC)/shdompp.f90 $(SRC)/fftpack.o -o $(BIN)/shdompp.e

disortsh: $(SRC)/disortsh.f90 $(SRC)/disort2.o
	$(FC) $(FFLAGS) $(SRC)/disortsh.f90 $(SRC)/disort2.o  -o $(BIN)/disortsh.e

ppmieprp: $(SRC)/ppmieprp.f90
	$(FC) $(FFLAGS) $(SRC)/ppmieprp.f90 -o $(BIN)/ppmieprp.e

.f.o :; $(FC) $(FFLAGS) -o $@ -c $<

put: $(SRC)/put.c
	$(GCC) $(FFLAGS) $(SRC)/put.c -o $(BIN)/put.e

.PHONY: test
test:
	cd tests ; ./run_shdompp_examples.sh

.PHONY: clean
clean:
	$(RM) $(SRC)/*.o tests/*.pp tests/*.out tests/*.prp tests/*.sfc
