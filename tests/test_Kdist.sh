#!/bin/csh
# Runs SHDOMPP in k-distribution mode using the shortwave RRTM k-distribution.

echo "==========================================="
echo "Runs SHDOMPP in k-distribution mode using"
echo " the shortwave RRTM k-distribution."
echo "==========================================="

#    Set up the atmospheric profile and the CKD file names
set atmfile=mls.atm
set ckdfile="data/ckd/swrrtm_mls_ckd.dat"
set miebase="data/mie/swrrtm_liq_"

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
    put.e O ${miebase}${ib}_mie.dat $specavgflag $delwave[$ib] $plancktemp \
        $cldphase  "$wavelen1[$ib] $wavelen2[$ib]" $clddisttype \
        $alpha "$Nretab $Sretab $Eretab" $maxNleg | cloudprp
    @ ib--
  end
  exit
endif


#    Make the correlated k-distribution data file
if (0) then
  # 360, 1.7 and 0.3 are the ppmv of CO2, CH4, and N2O
  put.e $atmfile $ckdfile 1.0 "360 1.7 0.3" | ~/kdist/rrtm/ckdswrrtm
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

  put.e 1 ${miebase}${ib}_mie.dat  $Nlay "$heights" "$temps" $LWPreff \
      "$molabs"  $molecoef[$ib] $prpfile | ppmieprp.e > /dev/null


  # Run SHDOMPP
  set solarmu=0.707
  set sfcalb=0.02         # Lambertian surface albedo
  set solarconstfac=1.0   # solar flux factor (correction to CKD file flux)
  set skyrad=0.0          # isotropic radiation incident on top boundary
  set Nmu=8  ;   set Nphi=16
  set splitacc=0.001 ;  set solacc=1.0E-5
  set accel=T  ;   set maxiter=100

  set outbase=kdist_shpp
  put.e $prpfile $ckdfile "$Nmu $Nphi" S "$solarconstfac $solarmu" $skyrad \
      L $sfcalb -1 "$waveno1[$ib] $waveno2[$ib]" $splitacc $solacc "$accel $maxiter" \
      ${outbase}${ib}f.out 1   /dev/null 0 $heights[1]  1 1.0  1 0  \
    | /usr/bin/time shdompp.e

  set outbase=kdist_dis
  put.e $prpfile $ckdfile $Nmu S "$solarconstfac $solarmu" $skyrad \
      $sfcalb "$waveno1[$ib] $waveno2[$ib]" \
      ${outbase}${ib}f.out  /dev/null 0 $heights[1]  1 1.0  1 0  \
    | /usr/bin/time disortsh.e

  @ ib++
end
