#!/bin/csh
# Runs 10.7 micron thermal radiative transfer in an atmosphere with
# water vapor and a three layer thin stratus cloud.

echo "==========================================="
echo "Runs 10.7 micron thermal radiative transfer"
echo " in an atmosphere with water vapor and a"
echo " three layer thin stratus cloud."
echo "==========================================="

# Make the SHDOMPP property file using ppmieprp.e
set prpfile=test3.pp
set Ncomp=1
set Nlay=7
set heights=(6.0 4.0 2.0 0.65  0.60  0.55  0.50  0.0)
set temps = (249 262 275 284.0 284.3 284.6 284.9 288.1)
set molabs = (0.0026 0.014 0.030 0.0015 0.0015 0.0015 0.020)
set raylcoef=0

put.e $Ncomp data/mie/water_w10.7_mie.dat $Nlay "$heights" "$temps" \
     0 0  0 0  0 0  13.90 9.77  8.34 8.24   2.78 5.72  0 0 \
   "$molabs"  $raylcoef $prpfile | ppmieprp.e > /dev/null


# Set up for SHDOMPP and DISORTSH runs
set waveno1=935.7 ; set waveno2=935.8
set Tsfc=283.0  ;  set Tsky=0.0
set sfcalb=0.02

# Run SHDOMPP
set Nmu=8 ;   set Nphi=16
set splitacc=0.001 ;  set solacc=1.0E-4
set accel=T  ;   set maxiter=30
set outbase=test3shpp
put.e $prpfile NONE "$Nmu $Nphi" T $Tsfc $Tsky  \
   L $sfcalb -1 "$waveno1 $waveno2" $splitacc $solacc "$accel $maxiter" \
   ${outbase}f.out 1   ${outbase}r.out 1 $heights[1] \
     6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000" 1 0 \
 | shdompp.e


# Run DISORTSH
set outbase=test3dis
put.e $prpfile NONE $Nmu T $Tsfc $Tsky $sfcalb  "$waveno1 $waveno2" \
   ${outbase}f.out  ${outbase}r.out 1 $heights[1] \
   6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000" 1 0  | disortsh.e
