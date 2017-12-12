# Makefile for SHDOMPP in distribution

#FFLAGS= -64 -ansi -O2
FFLAGS= -fast
FC= pgf90


all: shdompp disortsh ppmieprp

shdompp: shdompp.f90 fftpack.o
	$(FC) $(FFLAGS) shdompp.f90 fftpack.o -o shdompp

disortsh: disortsh.f90 disort2.o
	$(FC) $(FFLAGS) disortsh.f90 disort2.o  -o disortsh

ppmieprp: ppmieprp.f90 
	$(FC) $(FFLAGS) ppmieprp.f90 -o ppmieprp


.f.o :; $(FC) -c $(FFLAGS) $*.f

