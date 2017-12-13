#!/bin/csh
# Script to compare SHDOMPP to DISORT for a number of cases.
# The executable shdomppt should be compiled with the timing
# loop uncommented.  The executable disort${Nmu} (e.g. disort16)
# has the timing loop enabled and the MXCMU parameter in DISORT.f
# set to Nmu.


set MakeMieTables=0
set SolarAerosolCase=0
set SolarStCuCloudCase=0
set SolarThickCloudCase=0
set ThermalStCuCloudCase=1



if ($MakeMieTables) then
  # Runs cloudprp to make Mie table files for ppmieprp.e.
     # Water cloud droplets at 0.86 um
    set miefile="../data/mie/water_w0.86_mie.dat"
    set wavelen=0.86
    set alpha=7
    set Nretab=15 ;  set Sretab=1 ;  set Eretab=15
    set maxNleg=2000
    put.e O $miefile W "$wavelen $wavelen" G $alpha \
        "$Nretab $Sretab $Eretab" $maxNleg | cloudprp

     # Water cloud droplets at 2.13 um
    set miefile="../data/mie/water_w2.13_mie.dat"
    set wavelen=2.13
    set alpha=7
    set Nretab=15 ;  set Sretab=1 ;  set Eretab=15
    set maxNleg=2000
    put.e O $miefile W "$wavelen $wavelen" G $alpha \
        "$Nretab $Sretab $Eretab" $maxNleg | cloudprp

     # Water cloud droplets at 10.7 um
    set miefile="../data/mie/water_w10.7_mie.dat"
    set wavelen=10.7
    set alpha=7
    set Nretab=15 ;  set Sretab=1 ;  set Eretab=15
    set maxNleg=2000
    put.e O $miefile W "$wavelen $wavelen" G $alpha \
        "$Nretab $Sretab $Eretab" $maxNleg | cloudprp

     # Sea salt aerosols at 0.86 um
    set miefile="../data/mie/seasalt_w0.86_mie.dat"
    set wavelen=0.86
    set index="(1.28,-3.0e-6)"
    set alpha=0.70
    set Nretab=20 ;  set Sretab=0.05 ;  set Eretab=1.0
    set maxNleg=1000
    put.e O $miefile A "$index" "$wavelen $wavelen" L $alpha \
        "$Nretab $Sretab $Eretab" $maxNleg | cloudprp

     # Sea salt aerosols at 2.13 um
    set miefile="../data/mie/seasalt_w2.13_mie.dat"
    set wavelen=2.13
    set index="(1.445,-0.0015)"
    set alpha=0.70
    set Nretab=20 ;  set Sretab=0.05 ;  set Eretab=1.0
    set maxNleg=1000
    put.e O $miefile A "$index" "$wavelen $wavelen" L $alpha \
        "$Nretab $Sretab $Eretab" $maxNleg | cloudprp

     # Dust aerosols
    set miefile="../data/mie/dust_w0.55_mie.dat"
    set wavelen=0.55
    set index="(1.50,-0.002)"
    set alpha=0.70
    set Nretab=20 ;  set Sretab=0.1 ;  set Eretab=2.0
    set maxNleg=1000
    put.e O $miefile A "$index" "$wavelen $wavelen" L $alpha \
        "$Nretab $Sretab $Eretab" $maxNleg | cloudprp
endif





if ($SolarAerosolCase) then
  # Make the property file at 0.55 um
  set prpfile=dust_w0.55.pp
  set Ncomp=1
  set Nlay=7
  set heights=(40.0  20.0  11.0  4.0   3.0   2.0   1.0   0.0)
  set temps = (250.4 216.6 216.8 262.2 268.7 275.1 281.6 288.1)
  set molabs = (0.020 0.009 0 0 0 0 0 0)
  set raylcoef=0.00332

  put.e $Ncomp ../data/mie/dust_w0.55_mie.dat $Nlay "$heights" "$temps" \
       0 0  0 0  0 0   0.04 0.6  0.08 0.8  0.15 1.0  0.30 1.5 \
      "$molabs"  $raylcoef $prpfile | ppmieprp.e  > /dev/null



  # Run the radiative transfer comparison
  set solarmu=0.5
  set wavelen=0.55
  echo " "
  # Set up for SHDOMPP and DISORT runs
  set prpfile=dust_w${wavelen}.pp
  set solarflux=1.0  ;  set skyrad=0.0
  set sfcalb=0.15

  # Run the DISORT reference case (high accuracy)
  set Nmu=128
  set fluxref=dust${wavelen}disref_flux.out
  set radref=dust${wavelen}disref_rad.out
  if (1) then
    put.e $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
        $fluxref $radref 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180"  \
      1 | /usr/bin/time disortsh${Nmu} >&! disreflog${wavelen}.t
  endif
  set CPUref=`tail -2 disreflog${wavelen}.t | awk '{if (NR==1) print $1+0;}'`
  echo "DISORT Nmu=$Nmu reference CPU time: $CPUref" ; echo " "


  # Run the DISORT lower resolution cases
  echo "DISORT accuracy for solar reflectance from dust atmos  wavelength=$wavelen  mu0=$solarmu"
  echo " Nmu   CPU rmsFluxDif rmsRadFdif"
  foreach Nmu (4 8 16 32)
    @ Ntimes = 4096 / ( $Nmu * $Nmu)

    set fluxfile=dust${wavelen}dis_flux.out
    set radfile=dust${wavelen}dis_rad.out
    put.e $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
      $fluxfile  $radfile 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180"  \
      $Ntimes | /usr/bin/time disortsh${Nmu} >&! disortlog.t
    set CPU=`tail -2 disortlog.t | awk '{if (NR==1) print $1+0;}'`

    # Compare hemispheric flux and radiance files: rms difference
    paste $fluxref $fluxfile | awk -v Nmu=$Nmu -v cpu=$CPU  -v nt=$Ntimes \
       '{if ($1 != "!") {sum+=($6-$2)^2+($7-$3)^2; n+=2;}} \
       END {printf "%4d  %5.3f   %6.4f\n",Nmu,cpu/nt,sqrt(sum/n);}' >! fd.t
    paste $radref $radfile | awk -v Nmu=$Nmu '{if ($1 != "!") \
       {rmsdif+=(($8-$4)/$4)^2; n++;}}  END {rmsdif=sqrt(rmsdif/n); \
        printf "  %6.4f\n",rmsdif;}' >! rd.t
    paste -d " " fd.t rd.t
  end  # DISORT Nmu loop
  echo " "


  # Run the SHDOMPP cases
  echo "SHDOMPP accuracy for solar reflectance from dust atmos  wavelength=$wavelen  mu0=$solarmu"
  cp /dev/null shdompplog.t
  foreach solacc (1.0E-5 1.0E-4 1.0E-3)
    foreach splitacc (0.001 0)
      echo " Nmu   CPU rmsFluxDif rmsRadFdif iters Npts   splitacc=$splitacc  solacc=$solacc"
      foreach Nmu (4 8 16 32)
        @ Nphi = 2 * $Nmu
        set accel=T ;  set maxiter=300
        @ Ntimes = 16384 / ( $Nmu * $Nmu )

        set fluxfile=dust${wavelen}shpp_flux.out
        set radfile=dust${wavelen}shpp_rad.out
        put.e $prpfile NONE "$Nmu $Nphi"  S "$solarflux $solarmu" $skyrad \
           L $sfcalb $wavelen $splitacc $solacc "$accel $maxiter" \
           $fluxfile 1 $radfile 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
           $Ntimes | /usr/bin/time shdomppt >>& shdompplog.t
        set CPU=`tail -2 shdompplog.t | awk '{if (NR==1) print $1+0;}'`
        set iters=`cat $fluxfile | awk '{if (NR==14) {i=index($0,"NUMBER_ITERATIONS="); print substr($0,i+18,4);}}'`
        set npts=`cat $fluxfile | awk '{if (NR==3) {i=index($0,"NPTS="); print substr($0,i+5,4);}}'`

        # Compare hemispheric flux and radiance files: rms difference
        cat $fluxref | awk '{if ($1 != "!") print $0;}' >! ref.t
        cat $fluxfile | awk '{if ($1 != "!") print $0;}' >! exp.t
        paste ref.t exp.t | awk -v Nmu=$Nmu -v cpu=$CPU -v nt=$Ntimes \
          '{sum+=($6-$2)^2+($7-$3)^2; n+=2;} \
           END {printf "%4d  %5.3f   %6.4f\n",Nmu,cpu/nt,sqrt(sum/n);}' >! fd.t
        cat $radref | awk '{if ($1 != "!") print $0;}' >! ref.t
        cat $radfile | awk '{if ($1 != "!") print $0;}' >! exp.t
        paste ref.t exp.t | awk -v iters=$iters -v npts=$npts \
           '{rmsdif+=(($8-$4)/$4)^2; n++;}  END {rmsdif=sqrt(rmsdif/n); \
             printf "  %6.4f  %4d %4d\n",rmsdif,iters,npts;}' >! rd.t
        paste -d " " fd.t rd.t
      end  # SHDOMPP Nmu loop
    end  # SHDOMPP splitacc loop
  end  # SHDOMPP solacc loop
  echo " "; echo " "; echo " "
  rm -f  fd.t rd.t ref.t exp.t disortlog.t shdompplog.t
endif





if ($SolarStCuCloudCase) then

  # Make StCu cloud LWP and Reff profile: Linear LWC increase with height,
  #    constant number concentration (N)
  if (0) then
    awk 'BEGIN {nl=6; dz=0.050; lwcmax=0.667; N=100; a=lwcmax/(nl*dz);
      for (i=1; i<=nl; i++) {z=i*dz; W=a*z; P=0.5*a*(z^2-(z-dz)^2);
        re=100*(P/(dz*N*3))^0.333; tau=1500*P/re; printf "%5.2f %5.3f %6.2f %5.2f %5.2f\n",
           z,W,1000*P,re,tau; lwp+=1000*P; taut+=tau;}  print "LWP=",lwp,"   tau=",taut; }'
  endif

  # Make the property file at 0.86 um
  set prpfile=stcu_w0.86.pp
  set Ncomp=2
  set Nlay=8
  set heights=(10.0  0.80  0.75  0.70  0.65  0.60  0.55  0.50 0.0)
  set temps = (220   280.0 278.4 278.7 279.0 279.4 279.7 280  285)
  set molabs = (0 0 0 0 0 0 0 0)
  set raylcoef=0.00054

  put.e $Ncomp ../data/mie/water_w0.86_mie.dat seasalt_w0.86_mie.dat \
      $Nlay "$heights" "$temps" \
      "0 0"    "0 0" \
      "30.57 0" "12.70 0" \
      "25.01 0" "11.88 0" \
      "19.45 0" "10.93 0" \
      "13.90 0"  "9.77 0" \
       "8.34 0"  "8.24 0" \
       "2.78 0"  "5.72 0" \
       "0 0.02"  "0 0.5" \
      "$molabs"  $raylcoef $prpfile | ppmieprp.e > /dev/null

  # Make the property file at 2.13 um
  #  Molecular absorption from MODTRAN3 for US std atmos
  set prpfile=stcu_w2.13.pp
  set Ncomp=2
  set Nlay=8
  set heights=(10.0  0.80  0.75  0.70  0.65  0.60  0.55  0.50 0.0)
  set temps = (220   280.0 278.4 278.7 279.0 279.4 279.7 280  285)
  set molabs = (0.0173 0.0004 0.0004 0.0004 0.0004 0.0004 0.0004 0.0044)
  set raylcoef=0.0

  put.e $Ncomp ../data/mie/water_w2.13_mie.dat ../data/mie/seasalt_w2.13_mie.dat \
      $Nlay "$heights" "$temps" \
      "0 0"    "0 0" \
      "30.57 0" "12.70 0" \
      "25.01 0" "11.88 0" \
      "19.45 0" "10.93 0" \
      "13.90 0"  "9.77 0" \
       "8.34 0"  "8.24 0" \
       "2.78 0"  "5.72 0" \
       "0 0.02"  "0 0.5" \
      "$molabs"  $raylcoef $prpfile | ppmieprp.e > /dev/null


  # Run the radiative transfer comparison
  set solarmu=0.5

  foreach wavelen (0.86 2.13)
    echo " "
    # Set up for SHDOMPP and DISORT runs
    set prpfile=stcu_w${wavelen}.pp
    set solarflux=1.0  ;  set skyrad=0.0
    set sfcalb=0.0

    # Run the DISORT reference case (high accuracy)
    set Nmu=128
    set fluxref=stcu${wavelen}disref_flux.out
    set radref=stcu${wavelen}disref_rad.out
    if (1) then
     put.e $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
        $fluxref $radref 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180"  \
      1 | /usr/bin/time disortsh${Nmu} >&! disreflog${wavelen}.t
    endif
    set CPUref=`tail -2 disreflog${wavelen}.t | awk '{if (NR==1) print $1+0;}'`
    echo "DISORT Nmu=$Nmu reference CPU time: $CPUref" ; echo " "


    # Run the DISORT lower resolution cases
    echo "DISORT accuracy for solar reflectance from stcu cloud  wavelength=$wavelen  mu0=$solarmu"
    echo " Nmu   CPU rmsFluxDif rmsRadFdif"
    foreach Nmu (4 8 16 32)
      @ Ntimes = 16384 / ( $Nmu * $Nmu)

      set fluxfile=stcu${wavelen}dis_flux.out
      set radfile=stcu${wavelen}dis_rad.out
      put.e $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
        $fluxfile  $radfile 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180"  \
        $Ntimes | /usr/bin/time disortsh${Nmu} >&! disortlog.t
      set CPU=`tail -2 disortlog.t | awk '{if (NR==1) print $1+0;}'`

      # Compare hemispheric flux and radiance files: rms difference
      paste $fluxref $fluxfile | awk -v Nmu=$Nmu -v cpu=$CPU  -v nt=$Ntimes \
         '{if ($1 != "!") {sum+=($6-$2)^2+($7-$3)^2; n+=2;}} \
         END {printf "%4d  %5.3f   %6.4f\n",Nmu,cpu/nt,sqrt(sum/n);}' >! fd.t
      paste $radref $radfile | awk -v Nmu=$Nmu '{if ($1 != "!") \
         {rmsdif+=(($8-$4)/$4)^2; n++;}}  END {rmsdif=sqrt(rmsdif/n); \
          printf "  %6.4f\n",rmsdif;}' >! rd.t
      paste -d " " fd.t rd.t
    end  # DISORT Nmu loop
    echo " "


    # Run the SHDOMPP cases
    echo "SHDOMPP accuracy for solar reflectance from stcu cloud  wavelength=$wavelen  mu0=$solarmu"
    cp /dev/null shdompplog.t
    foreach solacc (1.0E-5 1.0E-4)
     foreach splitacc (0.003 0.01 0.03)
      echo " Nmu   CPU rmsFluxDif rmsRadFdif iters Npts   splitacc=$splitacc  solacc=$solacc"
      foreach Nmu (4 8 16 32)
        @ Nphi = 2 * $Nmu
        set accel=T ;  set maxiter=300
        @ Ntimes = 16384 / ( $Nmu * $Nmu )

        set fluxfile=stcu${wavelen}shpp_flux.out
        set radfile=stcu${wavelen}shpp_rad.out
        put.e $prpfile NONE "$Nmu $Nphi"  S "$solarflux $solarmu" $skyrad \
           L $sfcalb $wavelen $splitacc $solacc "$accel $maxiter" \
           $fluxfile 1 $radfile 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
           $Ntimes | /usr/bin/time shdomppt  >>& shdompplog.t
        set CPU=`tail -2 shdompplog.t | awk '{if (NR==1) print $1+0;}'`
        set iters=`cat $fluxfile | awk '{if (NR==14) {i=index($0,"NUMBER_ITERATIONS="); print substr($0,i+18,4);}}'`
        set npts=`cat $fluxfile | awk '{if (NR==3) {i=index($0,"NPTS="); print substr($0,i+5,4);}}'`

        # Compare hemispheric flux and radiance files: rms difference
        cat $fluxref | awk '{if ($1 != "!") print $0;}' >! ref.t
        cat $fluxfile | awk '{if ($1 != "!") print $0;}' >! exp.t
        paste ref.t exp.t | awk -v Nmu=$Nmu -v cpu=$CPU -v nt=$Ntimes \
          '{sum+=($6-$2)^2+($7-$3)^2; n+=2;} \
           END {printf "%4d  %5.3f   %6.4f\n",Nmu,cpu/nt,sqrt(sum/n);}' >! fd.t
        cat $radref | awk '{if ($1 != "!") print $0;}' >! ref.t
        cat $radfile | awk '{if ($1 != "!") print $0;}' >! exp.t
        paste ref.t exp.t | awk -v iters=$iters -v npts=$npts \
           '{rmsdif+=(($8-$4)/$4)^2; n++;}  END {rmsdif=sqrt(rmsdif/n); \
             printf "  %6.4f  %4d %4d\n",rmsdif,iters,npts;}' >! rd.t
        paste -d " " fd.t rd.t
      end  # SHDOMPP Nmu loop
     end  # SHDOMPP splitacc loop
    end  # SHDOMPP solacc loop

  end  # wavelength loop

  echo " "; echo " "; echo " "
  rm -f  fd.t rd.t ref.t exp.t disortlog.t  shdompplog.t
endif





if ($SolarThickCloudCase) then
  # Make the property file at 0.86 um
  set prpfile=thick_w0.86.pp
  set Ncomp=1
  set Nlay=1
  set heights=(1.0 0.0)
  set temps = (285 290)
  set molabs = (0)
  set raylcoef=0.0

  put.e $Ncomp ../data/mie/water_w0.86_mie.dat \
      $Nlay "$heights" "$temps"  600 10 \
      "$molabs"  $raylcoef $prpfile | ppmieprp.e > /dev/null

  # Run the radiative transfer comparison
  set solarmu=0.5
  set wavelen=0.86
  echo " "
  # Set up for SHDOMPP and DISORT runs
  set prpfile=thick_w${wavelen}.pp
  set solarflux=1.0  ;  set skyrad=0.0
  set sfcalb=0.0

  # Run the DISORT reference case (high accuracy)
  set Nmu=128
  set fluxref=thick${wavelen}disref_flux.out
  set radref=thick${wavelen}disref_rad.out
  if (1) then
     put.e $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
        $fluxref $radref 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180"  \
      1 | /usr/bin/time disortsh${Nmu} >&! disreflog${wavelen}t.t
  endif
  set CPUref=`tail -2 disreflog${wavelen}t.t | awk '{if (NR==1) print $1+0;}'`
  echo "DISORT Nmu=$Nmu reference CPU time: $CPUref" ; echo " "


  # Run the DISORT lower resolution cases
  echo "DISORT accuracy for solar reflectance from thick cloud  wavelength=$wavelen  mu0=$solarmu"
  echo " Nmu   CPU rmsFluxDif rmsRadFdif"
  foreach Nmu (4 8 16)
    @ Ntimes = 32768 / ( $Nmu * $Nmu)

    set fluxfile=thick${wavelen}dis_flux.out
    set radfile=thick${wavelen}dis_rad.out
    put.e $prpfile NONE $Nmu S "$solarflux $solarmu" $skyrad $sfcalb  \
        $fluxfile  $radfile 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180"  \
        $Ntimes | /usr/bin/time disortsh${Nmu} >&! disortlog.t
    set CPU=`tail -2 disortlog.t | awk '{if (NR==1) print $1+0;}'`

    # Compare hemispheric flux and radiance files: rms difference
    paste $fluxref $fluxfile | awk -v Nmu=$Nmu -v cpu=$CPU  -v nt=$Ntimes \
       '{if ($1 != "!") {sum+=($6-$2)^2+($7-$3)^2; n+=2;}} \
       END {printf "%4d %6.4f   %6.4f\n",Nmu,cpu/nt,sqrt(sum/n);}' >! fd.t
    paste $radref $radfile | awk -v Nmu=$Nmu '{if ($1 != "!") \
       {rmsdif+=(($8-$4)/$4)^2; n++;}}  END {rmsdif=sqrt(rmsdif/n); \
        printf "  %6.4f\n",rmsdif;}' >! rd.t
    paste -d " " fd.t rd.t
  end  # DISORT Nmu loop
  echo " "


  # Run the SHDOMPP cases
  echo "SHDOMPP accuracy for solar reflectance from thick cloud wavelength=$wavelen  mu0=$solarmu"
  cp /dev/null shdompplog.t
  foreach solacc (1.0E-5 1.0E-4)
    foreach splitacc (0.003 0.01 0.03)
      echo " Nmu   CPU rmsFluxDif rmsRadFdif iters Npts   splitacc=$splitacc  solacc=$solacc"
      foreach Nmu (4 8 16)
        @ Nphi = 2 * $Nmu
        set accel=T ;  set maxiter=300
        @ Ntimes = 4096 / ( $Nmu * $Nmu )

        set fluxfile=thick${wavelen}shpp_flux.out
        set radfile=thick${wavelen}shpp_rad.out
        put.e $prpfile NONE "$Nmu $Nphi"  S "$solarflux $solarmu" $skyrad \
           L $sfcalb $wavelen $splitacc $solacc "$accel $maxiter" \
           $fluxfile 1 $radfile 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
           $Ntimes | /usr/bin/time shdomppt  >>& shdompplog.t
        set CPU=`tail -2 shdompplog.t | awk '{if (NR==1) print $1+0;}'`
        set iters=`cat $fluxfile | awk '{if (NR==14) {i=index($0,"NUMBER_ITERATIONS="); print substr($0,i+18,4);}}'`
        set npts=`cat $fluxfile | awk '{if (NR==3) {i=index($0,"NPTS="); print substr($0,i+5,4);}}'`

        # Compare hemispheric flux and radiance files: rms difference
        cat $fluxref | awk '{if ($1 != "!") print $0;}' >! ref.t
        cat $fluxfile | awk '{if ($1 != "!") print $0;}' >! exp.t
        paste ref.t exp.t | awk -v Nmu=$Nmu -v cpu=$CPU -v nt=$Ntimes \
          '{sum+=($6-$2)^2+($7-$3)^2; n+=2;} \
           END {printf "%4d  %5.3f   %6.4f\n",Nmu,cpu/nt,sqrt(sum/n);}' >! fd.t
        cat $radref | awk '{if ($1 != "!") print $0;}' >! ref.t
        cat $radfile | awk '{if ($1 != "!") print $0;}' >! exp.t
        paste ref.t exp.t | awk -v iters=$iters -v npts=$npts \
           '{rmsdif+=(($8-$4)/$4)^2; n++;}  END {rmsdif=sqrt(rmsdif/n); \
             printf "  %6.4f  %4d %4d\n",rmsdif,iters,npts;}' >! rd.t
        paste -d " " fd.t rd.t
      end  # SHDOMPP Nmu loop
    end  # SHDOMPP splitacc loop
  end  # SHDOMPP solacc loop

  echo " "; echo " "; echo " "
  rm -f  fd.t rd.t ref.t exp.t disortlog.t shdompplog.t
endif






if ($ThermalStCuCloudCase) then

  # Make the property file at 10.7 um
  set prpfile=stcu_w10.7.pp
  set Ncomp=1
  set Nlay=10
  set heights=(6.0 4.0 2.0 0.80  0.75  0.70  0.65  0.60  0.55  0.50 0.0)
  set temps = (249 262 275 282.9 283.3 283.6 284.0 284.3 284.6 284.9 288.1)
  set molabs = (0.0026 0.014 0.025 0.0015 0.0015 0.0015 0.0015 0.0015 0.0015 0.020)
  set raylcoef=0

  put.e $Ncomp ../data/mie/water_w10.7_mie.dat $Nlay "$heights" "$temps" \
        0 0  0 0  0 0   30.57 12.70   25.01 11.88  19.45 10.93 \
       13.90  9.77  8.34  8.24   2.78  5.72  0 0 \
      "$molabs"  $raylcoef $prpfile | ppmieprp.e > /dev/null
  echo " "


  # Set up for SHDOMPP and DISORT runs
  set wavelen=10.7
  set prpfile=stcu_w${wavelen}.pp
  set waveno1=935.7 ; set waveno2=935.8
  set Tsfc=283.0  ;  set Tsky=0.0
  set sfcalb=0.02

  # Run the DISORT reference case (high accuracy)
  set Nmu=64
  set fluxref=stcu${wavelen}disref_flux.out
  set radref=stcu${wavelen}disref_rad.out
  put.e $prpfile NONE $Nmu T $Tsfc $Tsky $sfcalb  "$waveno1 $waveno2" \
      $fluxref  $radref 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180"  \
     10 | /usr/bin/time disortsh${Nmu} >&! disreflog${wavelen}.t
  set CPUref=`tail -2 disreflog${wavelen}.t | awk '{if (NR==1) print $1/10;}'`
  echo "Reference CPU time: $CPUref"

  # Run the DISORT lower resolution cases
  echo "DISORT accuracy for thermal emission from stcu cloud  wavelength=$wavelen"
  echo " Nmu   CPU rmsFluxFdif rmsRadFdif"
  foreach Nmu (4 8 16)
    @ Ntimes = 16384 / $Nmu

    set fluxfile=stcu${wavelen}dis_flux.out
    set radfile=stcu${wavelen}dis_rad.out
    put.e $prpfile NONE $Nmu T $Tsfc $Tsky $sfcalb "$waveno1 $waveno2"  \
        $fluxfile  $radfile 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180"  \
        $Ntimes | /usr/bin/time disortsh${Nmu} >&! disortlog.t
    set CPU=`tail -2 disortlog.t | awk '{if (NR==1) print $1+0;}'`

    # Compare hemispheric flux and radiance files: rms difference
    paste $fluxref $fluxfile | awk -v Nmu=$Nmu -v cpu=$CPU -v nt=$Ntimes \
       '{if ($1 != "!") {sum+=(($6-$2)/$2)^2; n++; if ($3>0) {sum+=(($7-$3)/$3)^2; n++;}}} \
       END {printf "%4d %6.4f   %6.4f\n",Nmu,cpu/nt,sqrt(sum/n);}' >! fd.t
    paste $radref $radfile | awk -v Nmu=$Nmu '{if ($1 != "!") \
       {rmsdif+=(($8-$4)/$4)^2; n++;}}  END {rmsdif=sqrt(rmsdif/n); \
        printf "  %6.4f\n",rmsdif;}' >! rd.t
    paste -d " " fd.t rd.t
  end  # DISORT Nmu loop

  echo " "

  # Run the SHDOMPP cases
  echo "SHDOMPP accuracy for thermal emission from stcu cloud  wavelength=$wavelen"
  cp /dev/null shdompplog.t
  foreach solacc (1.0E-5 1.0E-4)
   foreach splitacc (0.001 0.003 0.01)
    echo " Nmu   CPU rmsFluxFdif rmsRadFdif iters Npts   splitacc=$splitacc  solacc=$solacc"
    foreach Nmu (4 8 16)
      set Nphi = 1
      set accel=T ;  set maxiter=30
      @ Ntimes = 4096 / $Nmu

      set fluxfile=stcu${wavelen}shpp_flux.out
      set radfile=stcu${wavelen}shpp_rad.out
      put.e $prpfile NONE "$Nmu $Nphi" T $Tsfc $Tsky  \
          L $sfcalb -1 "$waveno1 $waveno2" $splitacc $solacc "$accel $maxiter" \
          $fluxfile 1 $radfile 1 $heights[1]  6 "0.2588 0.5000 0.7071 0.8660 0.9659 1.000"  5 "0 45 90 135 180" \
         $Ntimes | /usr/bin/time shdomppt >>& shdompplog.t
      set CPU=`tail -2 shdompplog.t | awk '{if (NR==1) print $1+0;}'`
      set iters=`cat $fluxfile | awk '{if (NR==13) {i=index($0,"NUMBER_ITERATIONS="); print substr($0,i+18,4);}}'`
      set npts=`cat $fluxfile | awk '{if (NR==3) {i=index($0,"NPTS="); print substr($0,i+5,4);}}'`

      # Compare hemispheric flux and radiance files: rms difference
      cat $fluxref | awk '{if ($1 != "!") print $0;}' >! ref.t
      cat $fluxfile | awk '{if ($1 != "!") print $0;}' >! exp.t
      paste ref.t exp.t | awk -v Nmu=$Nmu -v cpu=$CPU -v nt=$Ntimes \
         '{sum+=(($6-$2)/$2)^2; n++; if ($3>0) {sum+=(($7-$3)/$3)^2; n++;}} \
          END {printf "%4d %6.4f   %6.4f\n",Nmu,cpu/nt,sqrt(sum/n);}' >! fd.t
      cat $radref | awk '{if ($1 != "!") print $0;}' >! ref.t
      cat $radfile | awk '{if ($1 != "!") print $0;}' >! exp.t
      paste ref.t exp.t | awk -v iters=$iters -v npts=$npts \
         '{rmsdif+=(($8-$4)/$4)^2; n++;}  END {rmsdif=sqrt(rmsdif/n); \
           printf "  %6.4f  %4d %4d\n",rmsdif,iters,npts;}' >! rd.t
      paste -d " " fd.t rd.t
    end  # SHDOMPP Nmu loop
    end  # SHDOMPP splitacc loop
   end  # SHDOMPP solacc loop

  echo " "; echo " "; echo " "
  rm -f  fd.t rd.t ref.t exp.t disortlog.t  shdompplog.t
endif
