#!/bin/csh
# Script with examples of running SHDOMPP.
#
# To run all sections requires ppmieprp, shdompp, disortsh, and put in
# the current path and Mie table files dust_w0.55.mie and water_w10.7.mie
# in the current directory.  The BRDF test also requires shdom.  The
# k-distribution test requires swrrtm_mls.ckd and swrrtm_liq_8.mie.


set SolarHGcloudTest=0
set SolarMieAerosolTest=0
set ThermalCloudTest=0
set BRDFtest=0
set RunKdistTest=1


if ($SolarHGcloudTest) then
 # Runs solar radiative transfer in a single layer cloud with a
 #   Henyey-Greenstein phase function

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

  put $prpfile NONE "$Nmu $Nphi" S "$solarflux $solarmu" $skyrad \
      L $sfcalb $wavelen $splitacc $solacc "$accel $maxiter" \
      ${outbase}f.out 1  ${outbase}r.out 1 1.0 \
      6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
     | shdompp

  # Run DISORTSH
  set outbase=test1dis
  put $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
      ${outbase}f.out  ${outbase}r.out 1 1.0 \
      6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
     | disortsh
endif




if ($SolarMieAerosolTest) then
  # Run solar radiative transfer in a dust aerosol atmosphere

  if (0) then
     # Make Mie table with cloudprp (from SHDOM distribution)
    set miefile=dust_w0.55.mie
    set wavelen=0.55
    set index="(1.50,-0.002)"
    set alpha=0.70
    set Nretab=20 ;  set Sretab=0.1 ;  set Eretab=2.0
    set maxNleg=500
    put O $miefile A "$index" "$wavelen $wavelen" L $alpha \
        "$Nretab $Sretab $Eretab" $maxNleg | cloudprp
  endif


  # Make the SHDOMPP property file using ppmieprp
  set prpfile=test2.pp
  set Ncomp=1
  set Nlay=3
  set heights=(11.0 4.0 2.0 0.0)
  set temps = (217  262 275 288)
  set molabs = (0 0 0)
  set raylcoef=0.00332

  put $Ncomp dust_w0.55.mie $Nlay "$heights" "$temps" \
       0 0  0.10 1.0  0.20 1.5  "$molabs"  $raylcoef $prpfile | ppmieprp


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
  put $prpfile NONE "$Nmu $Nphi" S "$solarflux $solarmu" $skyrad \
      L $sfcalb $wavelen $splitacc $solacc "$accel $maxiter" \
      ${outbase}f.out 1  ${outbase}r.out 1 $heights[1] \
      6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
     | shdompp


  # Run DISORTSH
  set outbase=test2dis
  put $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
      ${outbase}f.out  ${outbase}r.out 1 $heights[1]  \
      6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
     | disortsh
endif





if ($ThermalCloudTest) then
 # Runs 10.7 micron thermal radiative transfer in an atmosphere with
 # water vapor and a three layer thin stratus cloud.

  # Make the SHDOMPP property file using ppmieprp
  set prpfile=test3.pp
  set Ncomp=1
  set Nlay=7
  set heights=(6.0 4.0 2.0 0.65  0.60  0.55  0.50  0.0)
  set temps = (249 262 275 284.0 284.3 284.6 284.9 288.1)
  set molabs = (0.0026 0.014 0.030 0.0015 0.0015 0.0015 0.020)
  set raylcoef=0

  put $Ncomp water_w10.7.mie $Nlay "$heights" "$temps" \
        0 0  0 0  0 0  13.90 9.77  8.34 8.24   2.78 5.72  0 0 \
      "$molabs"  $raylcoef $prpfile | ppmieprp > /dev/null


  # Set up for SHDOMPP and DISORTSH runs
  set waveno1=935.7 ; set waveno2=935.8
  set Tsfc=283.0  ;  set Tsky=0.0
  set sfcalb=0.02

  # Run SHDOMPP
  set Nmu=8 ;   set Nphi=16
  set splitacc=0.001 ;  set solacc=1.0E-4
  set accel=T  ;   set maxiter=30
  set outbase=test3shpp
  put $prpfile NONE "$Nmu $Nphi" T $Tsfc $Tsky  \
      L $sfcalb -1 "$waveno1 $waveno2" $splitacc $solacc "$accel $maxiter" \
      ${outbase}f.out 1   ${outbase}r.out 1 $heights[1] \
        6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000" 1 0 \
    | shdompp


  # Run DISORTSH
  set outbase=test3dis
  put $prpfile NONE $Nmu T $Tsfc $Tsky $sfcalb  "$waveno1 $waveno2" \
      ${outbase}f.out  ${outbase}r.out 1 $heights[1] \
      6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000" 1 0  | disortsh
endif





if ($BRDFtest) then
  # Compares SHDOMPP and SHDOM RPV surface BRDF with an aerosol atmosphere.
  #   The agreement is limited by a potential SHDOM bug for integrating
  #   the downwelling radiance at the surface for the BRDF calculation.

  # Make a single layer SHDOMPP property file using ppmieprp
  set prpfile=test4.pp
  set Ncomp=1  ;   set Nlay=1
  set heights=(2.0 0.0)
  set temps = (275 288)
  set molabs = (0)
  set raylcoef=0.00332
  put $Ncomp dust_w0.55.mie $Nlay "$heights" "$temps" \
      0.1 1.0  "$molabs"  $raylcoef $prpfile | ppmieprp


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
  put $prpfile NONE "$Nmu $Nphi" S "$solarflux $solarmu" $skyrad \
      R "$rho0 $k $Theta" $wavelen $splitacc $solacc "$accel $maxiter" \
      /dev/null 1  ${outbase}r.out 1 $heights[1] \
      6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
     | shdompp


  # Translate SHDOMPP property file to SHDOM property file (for single layer)
  set shdomprpfile=test4.prp
  cat $prpfile | awk '{if (NR==1) {if ($1 != 1) exit;} \
          if (NR==2) {Ztop=$1; Ttop=$2;} \
          if (NR==3) {Zbot=$1; Tbot=$2; tau=$3; alb=$4; nleg=$5; \
               for (l=1; l<=nleg; l++) leg[l]=$(l+5);} } \
       END {ext=tau/(Ztop-Zbot); print "1 1 2"; print "1.0 1.0 ",Zbot,Ztop;\
            printf "%d %d %d %5.1f %7.3f %7.5f %4d",1,1,2,Ttop,ext,alb,nleg;\
              for (l=1; l<=nleg; l++) printf " %8.5f",leg[l]; printf "\n"; \
            printf "%d %d %d %5.1f %7.3f %7.5f %4d",1,1,1,Tbot,ext,alb,nleg;\
              for (l=1; l<=nleg; l++) printf " %8.5f",leg[l]; printf "\n"; }' \
       >! $shdomprpfile


  # Make the SHDOM surface file:
  set sfcfile=test4.sfc
  put R "1 1 2.0 2.0" "1 1 290 $rho0 $k $Theta" >! $sfcfile

  set outbase=test4shdom
  set sfcalb=0.1
  set splitacc=0.0001
  set deltaM=T
  put $shdomprpfile $sfcfile NONE NONE NONE "1 1 3" "$Nmu $Nphi" 0 0 $deltaM E \
      S "$solarflux $solarmu 0" $skyrad $sfcalb \
      "$splitacc 0.0" "$accel $solacc $maxiter" \
      1 R "1 0 0 0 0 30 \
       0.2588 0 0.2588 45 0.2588 90 0.2588 135 0.2588 180 \
       0.5000 0 0.5000 45 0.5000 90 0.5000 135 0.5000 180 \
       0.7071 0 0.7071 45 0.7071 90 0.7071 135 0.7071 180 \
       0.8660 0 0.8660 45 0.8660 90 0.8660 135 0.8660 180 \
       0.9659 0 0.9659 45 0.9659 90 0.9659 135 0.9659 180 \
       1.0000 0 1.0000 45 1.0000 90 1.0000 135 1.0000 180" \
      ${outbase}r.out  | ~/shdom/shdom

  # Reformat radiance outputs and compare SHDOMPP with SHDOM
  cat test4shdomr.out | awk '{if (NR>19) {if ($1=="!") {mu=$2; phi=$3;} \
     else printf "%7.4f %5.1f %12.5E\n",mu,phi,$3; }}' >! test4shdomr.t
  cat test4shppr.out | awk '{if ($1!="!") \
     {printf "%7.4f %5.1f %12.5E\n",$2,$3,$4;}}' >! test4shppr.t
  paste test4shdomr.t test4shppr.t | awk '{sum+=(($6-$3)/$3)^2; n++;} \
     END {print "rms fractional radiance difference=",sqrt(sum/n);}'
  rm -f test4shdomr.t test4shppr.t
endif





if ($RunKdistTest) then
  # Runs SHDOMPP in k-distribution mode using the shortwave RRTM k-distribution.

  #    Set up the atmospheric profile and the CKD file names
  set atmfile=mls.atm
  set ckdfile=swrrtm_mls.ckd
  set miebase="swrrtm_liq_"

  #  Cloudprp stuff:
  #   Set the starting and ending bands of the shortwave RRTM k-distribution
  set Sband=8
  set Eband=8
  #    Set the cloud microphysical properties that aren't in the LWC file
  set dropconc=100              # droplet concentration (/cm^3)
  set clddisttype=G             # G for gamma distribution
  set alpha = 7                 # gamma dist shape param
  set cldphase="W"              # water phase: I for ice, W for water
  #   Set the scattering property table effective radius parameters
  set Nretab=50;  set Sretab=0.5;  set Eretab=25.0;
  set maxNleg=2000
  #   Set the Mie spectral averaging parameters
  set specavgflag="A";  set plancktemp=5800
  set delwave=(.03 .04 .04 .05 .05 .05 .05 .05 .05 .01 .05 .05 .05 .1)

  #  Set the Rayleigh molecular scattering coefficients for swrrtm bands
  set molecoef=(0.1463370 3.5541952E-02 1.3056834E-02 4.2549581E-03 1.3035705E-03 3.8487033E-04 1.1156387E-04 6.8169866E-05 2.9846953E-05 1.6589693E-05 1.0069725E-05 4.9696632E-06 2.1081225E-06 6.7372980e-7)

  #  Set the solar flux in each swrrtm band
  set solarflux=(3.0799 48.3722 129.4950 347.1923 218.1870 345.7425 24.2936 102.9314 55.6266 22.4277 23.7297 20.3651 12.1096 12.8894)
  set solarconst=1366.44

  #    The wavenumber and wavelength ranges for swrrtm's bands
  set waveno1=(38000 29000 22650 16000 12850 8050 7700 6150 5150 4650 4000 3250 2600 820)
  set waveno2=(50000 38000 29000 22650 16000 12850 8050 7700 6150 5150 4650 4000 3250 2600)
  set wavelen1=(.2   .263 .344 .44  .625 .78 1.24 1.30 1.63 1.94 2.15 2.50  3.077 3.846)
  set wavelen2=(.263 .344 .440 .625 .78 1.24 1.30 1.63 1.94 2.15 2.50 3.077 3.846 12.195)

  #    If needed make the Mie tables
  if (0) then
    set ib=$Eband
    while ($ib >= $Sband)
      put O ${miebase}${ib}.mie $specavgflag $delwave[$ib] $plancktemp \
          $cldphase  "$wavelen1[$ib] $wavelen2[$ib]" $clddisttype \
          $alpha "$Nretab $Sretab $Eretab" $maxNleg | cloudprp
      @ ib--
    end
    exit
  endif


  #    Make the correlated k-distribution data file
  if (0) then
    # 360, 1.7 and 0.3 are the ppmv of CO2, CH4, and N2O
    put $atmfile $ckdfile 1.0 "360 1.7 0.3" | ~/kdist/rrtm/ckdswrrtm
  endif


  # Loop over the RRTM bands
  set ib=$Sband
  while ($ib <= $Eband)

    # Make the property file
    set prpfile=stcu_b${ib}.pp
    set Nlay=9
    set heights=(13.0 6.0  0.80  0.75  0.70  0.65  0.60  0.55  0.50 0.0)
    set temps = (216  261  280.0 278.4 278.7 279.0 279.4 279.7 280  285)
    set LWPreff = (0 0  0 0  30.57 12.70  25.01 11.88  19.45 10.93 \
                   13.90 9.77  8.34 8.24  2.78 5.72   0 0)
    set molabs = (0 0 0 0 0 0 0 0 0)

    put 1 ${miebase}${ib}.mie  $Nlay "$heights" "$temps" $LWPreff \
        "$molabs"  $molecoef[$ib] $prpfile | ppmieprp > /dev/null


    # Run SHDOMPP
    set solarmu=0.707
    set sfcalb=0.02         # Lambertian surface albedo
    set solarconstfac=1.0   # solar flux factor (correction to CKD file flux)
    set skyrad=0.0          # isotropic radiation incident on top boundary
    set Nmu=8  ;   set Nphi=16
    set splitacc=0.001 ;  set solacc=1.0E-5
    set accel=T  ;   set maxiter=100

    set outbase=kdist_shpp
    put $prpfile $ckdfile "$Nmu $Nphi" S "$solarconstfac $solarmu" $skyrad \
        L $sfcalb -1 "$waveno1[$ib] $waveno2[$ib]" $splitacc $solacc "$accel $maxiter" \
        ${outbase}${ib}f.out 1   /dev/null 0 $heights[1]  1 1.0  1 0  \
      | /usr/bin/time shdompp

    set outbase=kdist_dis
    put $prpfile $ckdfile $Nmu S "$solarconstfac $solarmu" $skyrad \
        $sfcalb "$waveno1[$ib] $waveno2[$ib]" \
        ${outbase}${ib}f.out  /dev/null 0 $heights[1]  1 1.0  1 0  \
      | /usr/bin/time disortsh

    @ ib++
  end

endif
