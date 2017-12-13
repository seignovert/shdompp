#!/bin/csh
# Run solar radiative transfer in a dust aerosol atmosphere

echo "==========================================="
echo "Run solar radiative transfer in a dust"
echo "aerosol atmosphere"
echo "==========================================="

# FIXME: Missing dependancies (cloudprp). Input `dust_w0.55_mie.dat` is provided
# if (0) then
#    # Make Mie table with cloudprp (from SHDOM distribution)
#   set miefile="data/mie/dust_w0.55_mie.dat"
#   set wavelen=0.55
#   set index="(1.50,-0.002)"
#   set alpha=0.70
#   set Nretab=20 ;  set Sretab=0.1 ;  set Eretab=2.0
#   set maxNleg=500
#   put.e O $miefile A "$index" "$wavelen $wavelen" L $alpha \
#       "$Nretab $Sretab $Eretab" $maxNleg | cloudprp
# endif


# Make the SHDOMPP property file using ppmieprp.e
set prpfile=test2.pp
set Ncomp=1
set Nlay=3
set heights=(11.0 4.0 2.0 0.0)
set temps = (217  262 275 288)
set molabs = (0 0 0)
set raylcoef=0.00332

put.e $Ncomp data/mie/dust_w0.55_mie.dat $Nlay "$heights" "$temps" \
     0 0  0.10 1.0  0.20 1.5  "$molabs"  $raylcoef $prpfile | ppmieprp.e


# Set some parameters for the radiative transfer
set solarflux=1.0      # solar flux on a horizontal surface
set solarmu=0.866      # cosine solar zenith angles
set sfcalb=0.10        # Lambertian surface albedo
set skyrad=0.0         # isotropic radiation incident on top boundary

# Run SHDOMPP
set Nmu=16  ;   set Nphi=32
set wavelen=0.5
set splitacc=0.001 ;  set solacc=1.0E-5
set accel=T  ;   set maxiter=30
set outbase=test2shpp
put.e $prpfile NONE "$Nmu $Nphi" S "$solarflux $solarmu" $skyrad \
    L $sfcalb $wavelen $splitacc $solacc "$accel $maxiter" \
    ${outbase}f.out 1  ${outbase}r.out 1 $heights[1] \
    6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
   | shdompp.e


# # Run DISORTSH
# # NOTE: This test fail: DISORT--input and/or dimension errors
# set outbase=test2dis
# put.e $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
#     ${outbase}f.out  ${outbase}r.out 1 $heights[1]  \
#     6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
#    | disortsh.e
