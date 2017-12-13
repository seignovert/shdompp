#!/bin/csh
# Runs solar radiative transfer in a single layer cloud with a
# Henyey-Greenstein phase function

echo "==========================================="
echo "Runs solar radiative transfer in a single"
echo " layer cloud with a Henyey-Greenstein"
echo "phase function"
echo "==========================================="

# First make a single layer SHDOMPP property file with a H-G phase function
set tau=10.0
set omega=0.99
set g=0.85
set Nleg=128
set temp1=260
set temp2=280
set prpfile=test1hg.pp
awk -v tau=$tau -v omega=$omega -v g=$g -v nl=$Nleg -v T1=$temp1 -v T2=$temp2 \
   'BEGIN {print 1; printf "%5.1f %5.1f\n",1.0,T1; \
           printf "%5.1f %5.1f %6.2f %6.4f %d",0.0,T2,tau,omega,nl; \
           for (l=1; l<=nl; l++) printf " %8.5f",(2*l+1)*g^l; printf "\n";}' \
   >! $prpfile


# Set some parameters for the radiative transfer
set solarflux=1.0      # solar flux on a horizontal surface
set solarmu=0.5        # cosine solar zenith angles
set sfcalb=0.2        # Lambertian surface albedo
set skyrad=0.0         # isotropic radiation incident on top boundary

# Set SHDOMPP parameters and run it
set Nmu=8  ;   set Nphi=16
set wavelen=0.5
set splitacc=0.0003 ;  set solacc=1.0E-4
set accel=T  ;   set maxiter=100
set outbase=test1shpp

put.e $prpfile NONE "$Nmu $Nphi" S "$solarflux $solarmu" $skyrad \
    L $sfcalb $wavelen $splitacc $solacc "$accel $maxiter" \
    ${outbase}f.out 1  ${outbase}r.out 1 1.0 \
    6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
   | shdompp.e

# Run DISORTSH
set outbase=test1dis
put.e $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
    ${outbase}f.out  ${outbase}r.out 1 1.0 \
    6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
   | disortsh.e
