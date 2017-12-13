#!/bin/csh
# Compares SHDOMPP and SHDOM RPV surface BRDF with an aerosol atmosphere.

echo "==========================================="
echo "Compares SHDOMPP and SHDOM RPV surface BRDF"
echo " with an aerosol atmosphere."
echo "========================================"

#   The agreement is limited by a potential SHDOM bug for integrating
#   the downwelling radiance at the surface for the BRDF calculation.

# Make a single layer SHDOMPP property file using ppmieprp.e
set prpfile=test4.pp
set Ncomp=1  ;   set Nlay=1
set heights=(2.0 0.0)
set temps = (275 288)
set molabs = (0)
set raylcoef=0.00332
put.e $Ncomp data/mie/dust_w0.55_mie.dat $Nlay "$heights" "$temps" \
    0.1 1.0  "$molabs"  $raylcoef $prpfile | ppmieprp.e


# Set some parameters for the radiative transfer
set solarflux=1.0      # solar flux on a horizontal surface
set solarmu=0.866      # cosine solar zenith angles
set skyrad=0.0         # isotropic radiation incident on top boundary
set rho0=0.2 ;   set k=1.0 ;  set Theta=-0.24  # RPV parameters

# Run SHDOMPP
#  set Nmu=32  ;   set Nphi=64
set Nmu=16  ;   set Nphi=32
set wavelen=0.5
set splitacc=0.001 ;  set solacc=1.0E-5
set accel=T  ;   set maxiter=30
set outbase=test4shpp
put.e $prpfile NONE "$Nmu $Nphi" S "$solarflux $solarmu" $skyrad \
    R "$rho0 $k $Theta" $wavelen $splitacc $solacc "$accel $maxiter" \
    /dev/null 1  ${outbase}r.out 1 $heights[1] \
    6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
   | shdompp.e


# FIXME: Missing dependancies (shdom)

# # Translate SHDOMPP property file to SHDOM property file (for single layer)
# set shdomprpfile=test4.prp
# cat $prpfile | awk '{if (NR==1) {if ($1 != 1) exit;} \
#         if (NR==2) {Ztop=$1; Ttop=$2;} \
#         if (NR==3) {Zbot=$1; Tbot=$2; tau=$3; alb=$4; nleg=$5; \
#              for (l=1; l<=nleg; l++) leg[l]=$(l+5);} } \
#      END {ext=tau/(Ztop-Zbot); print "1 1 2"; print "1.0 1.0 ",Zbot,Ztop;\
#           printf "%d %d %d %5.1f %7.3f %7.5f %4d",1,1,2,Ttop,ext,alb,nleg;\
#             for (l=1; l<=nleg; l++) printf " %8.5f",leg[l]; printf "\n"; \
#           printf "%d %d %d %5.1f %7.3f %7.5f %4d",1,1,1,Tbot,ext,alb,nleg;\
#             for (l=1; l<=nleg; l++) printf " %8.5f",leg[l]; printf "\n"; }' \
#      >! $shdomprpfile
#
#
# # Make the SHDOM surface file:
# set sfcfile=test4.sfc
# put.e R "1 1 2.0 2.0" "1 1 290 $rho0 $k $Theta" >! $sfcfile
#
# set outbase=test4shdom
# set sfcalb=0.1
# set splitacc=0.0001
# set deltaM=T
# put.e $shdomprpfile $sfcfile NONE NONE NONE "1 1 3" "$Nmu $Nphi" 0 0 $deltaM E \
#     S "$solarflux $solarmu 0" $skyrad $sfcalb \
#     "$splitacc 0.0" "$accel $solacc $maxiter" \
#     1 R "1 0 0 0 0 30 \
#      0.2588 0 0.2588 45 0.2588 90 0.2588 135 0.2588 180 \
#      0.5000 0 0.5000 45 0.5000 90 0.5000 135 0.5000 180 \
#      0.7071 0 0.7071 45 0.7071 90 0.7071 135 0.7071 180 \
#      0.8660 0 0.8660 45 0.8660 90 0.8660 135 0.8660 180 \
#      0.9659 0 0.9659 45 0.9659 90 0.9659 135 0.9659 180 \
#      1.0000 0 1.0000 45 1.0000 90 1.0000 135 1.0000 180" \
#     ${outbase}r.out  | ~/shdom/shdom
#
# # Reformat radiance outputs and compare SHDOMPP with SHDOM
# cat test4shdomr.out | awk '{if (NR>19) {if ($1=="!") {mu=$2; phi=$3;} \
#    else printf "%7.4f %5.1f %12.5E\n",mu,phi,$3; }}' >! test4shdomr.t
# cat test4shppr.out | awk '{if ($1!="!") \
#    {printf "%7.4f %5.1f %12.5E\n",$2,$3,$4;}}' >! test4shppr.t
# paste test4shdomr.t test4shppr.t | awk '{sum+=(($6-$3)/$3)^2; n++;} \
#    END {print "rms fractional radiance difference=",sqrt(sum/n);}'
# rm -f test4shdomr.t test4shppr.t
