! SHDOMPP: Spherical harmonic discrete ordinate radiative transfer method
!          for plane-parallel atmospheres.
!     See shdompp.doc for documentation.

PROGRAM SHDOMPP
IMPLICIT NONE

  !  MAXDELTAU  maximum scattering optical depth allowed in initial sublayers
  !  SPLITTING_FACTOR  max number of internal grid points is this factor
  !                    times the number of input layers + optical depth
  REAL, PARAMETER :: MAXDELTAU = 1.0
  REAL, PARAMETER :: SPLITTING_FACTOR = 10.0
  LOGICAL, PARAMETER :: PRNT=.TRUE.   ! True for interation printing
!  LOGICAL, PARAMETER :: PRNT=.FALSE.

  INTEGER, PARAMETER   :: MAXPAR=100  ! Max array size for user input arrays
  REAL                 :: ZOUT(MAXPAR), MUOUT(MAXPAR), PHIOUT(MAXPAR)
  REAL                 :: SFCPARMS(MAXPAR)

  INTEGER :: NLAY, MAXLEG, NPTS, MAXIG
  INTEGER :: NMU, NPHI, MAXITER, ITER, NUMITER
  INTEGER :: NZOUT, NMUOUT, NPHIOUT, FLUXPARM
  INTEGER :: NCASES, ICASE, NG, NZCKD, IG
  LOGICAL :: KDIST, ACCELFLAG
  REAL    :: SOLARFLUX, SOLFLUX, SOLARMU
  REAL    :: SFCTEMP, SKYRAD, WAVELEN, WAVENO(2)
  REAL    :: SOLACC, SOLCRIT, SPLITACC
  CHARACTER(LEN=1) :: SRCTYPE, SFCTYPE
  CHARACTER(LEN=80) :: PROPFILE, CKDFILE, FLUXFILE, RADFILE

  INTEGER, ALLOCATABLE :: NLEGP(:), IXP(:,:)
  REAL,    ALLOCATABLE :: TAUP(:), ALBEDOP(:), LEGENP(:,:)
  REAL,    ALLOCATABLE :: TAUPM(:), ALBEDOPM(:)
  REAL,    ALLOCATABLE :: HEIGHTP(:), TEMPP(:)
  REAL,    ALLOCATABLE :: DELG(:), ZCKD(:), KABS(:,:)
  REAL,    ALLOCATABLE :: HEIGHTG(:), TAUG(:)
  REAL,    ALLOCATABLE :: FLUXUP(:), FLUXDN(:), MEANRAD(:), RADOUT(:,:,:)
  REAL,    ALLOCATABLE :: FLUXUPM(:), FLUXDNM(:), MEANRADM(:), RADOUTM(:,:,:)


  CALL USER_INPUT (PROPFILE, CKDFILE, KDIST, NMU, NPHI, &
                   SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                   SFCTEMP, SFCTYPE, SFCPARMS, WAVELEN, WAVENO, &
                   ACCELFLAG, SOLACC, MAXITER, SPLITACC, &
                   FLUXFILE, FLUXPARM,  RADFILE, MAXPAR, NZOUT, ZOUT, &
                   NMUOUT, MUOUT, NPHIOUT, PHIOUT)

  ALLOCATE (RADOUT(NZOUT,NMUOUT,NPHIOUT), RADOUTM(NZOUT,NMUOUT,NPHIOUT))

  ! Find out the size of the property grid arrays for allocating arrays
  CALL READ_PROPERTY_SIZE (PROPFILE, NLAY, MAXLEG)
  ALLOCATE (HEIGHTP(NLAY+1), TEMPP(NLAY+1))
  ALLOCATE (TAUP(NLAY), ALBEDOP(NLAY), TAUPM(NLAY), ALBEDOPM(NLAY))
  ALLOCATE (NLEGP(NLAY), LEGENP(0:MAXLEG,NLAY),  IXP(2,NLAY))

  ! Read in the properties of the medium
  CALL READ_PROPERTIES (PROPFILE, NLAY, MAXLEG, &
                        HEIGHTP, TEMPP, TAUP, ALBEDOP, NLEGP, LEGENP)

  ! Guess the maximum size of the adaptive grid and allocate arrays
  MAXIG = NINT(SPLITTING_FACTOR*(NLAY+SUM(TAUP(:)*ALBEDOP(:))/MAXDELTAU))
  ALLOCATE (HEIGHTG(MAXIG), TAUG(MAXIG))
  ALLOCATE (FLUXUP(MAXIG), FLUXDN(MAXIG), MEANRAD(MAXIG))
  ALLOCATE (FLUXUPM(MAXIG), FLUXDNM(MAXIG), MEANRADM(MAXIG))


  ! If doing a k-distribution then get the band info from the CKD file
  IF (KDIST) THEN
    CALL READ_CKD_SIZE (CKDFILE, WAVENO, NG, NZCKD)
    ALLOCATE (DELG(NG), ZCKD(NZCKD), KABS(NZCKD,NG))
    CALL READ_CKD (CKDFILE, WAVENO, NG, NZCKD, SOLFLUX, DELG, ZCKD, KABS)
    SOLARFLUX = SOLFLUX*SOLARFLUX*ABS(SOLARMU) 
  ELSE
    NG = 1
    NZCKD = 0
    ALLOCATE (DELG(NG), ZCKD(NZCKD), KABS(NZCKD,NG))
    DELG(1) = 1.0
  ENDIF
      

!  WRITE (*,*) 'Number of times to run SHDOMPP routine for timing'
!  READ (*,*) NCASES
!    WRITE (*,*) NCASES
!  DO ICASE = 1, NCASES

  ! Loop over the k-distribution g's from low to high absorption
  NUMITER = 0
  DO IG = 1, NG
    ! Add in the layer gaseous absorption if doing a k-distribution
    CALL ADD_IN_KDIST (NLAY, HEIGHTP, TAUP, ALBEDOP, &
                       NZCKD, ZCKD, KABS(:,IG),  TAUPM, ALBEDOPM)

    CALL SOLVE_SHDOMPP (NLAY, MAXLEG, HEIGHTP, TEMPP, TAUPM, ALBEDOPM, &
                        NLEGP, LEGENP, SRCTYPE, SOLARFLUX, SOLARMU, &
                        WAVELEN, WAVENO, SKYRAD, SFCTYPE, SFCPARMS, SFCTEMP, &
                        NMU, NPHI, MAXIG, MAXDELTAU, SPLITACC, SOLACC, &
                        MAXITER, ACCELFLAG, PRNT,  SOLCRIT, ITER, &
                        NPTS, HEIGHTG, TAUG, IXP, FLUXUPM, FLUXDNM, MEANRADM, &
                        NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                        RADOUTM)

    NUMITER = NUMITER + ITER
    IF (KDIST) THEN
      CALL SUM_KDIST_OUTPUT (IG, DELG(IG), NLAY, TAUP, IXP, NPTS, TAUG, &
                             HEIGHTP, HEIGHTG, FLUXUPM, FLUXUP, &
                             FLUXDNM, FLUXDN, MEANRADM, MEANRAD, &
                             NZOUT, NMUOUT, NPHIOUT, RADOUTM, RADOUT)
    ELSE
      FLUXUP(1:NPTS) = FLUXUPM(1:NPTS)
      FLUXDN(1:NPTS) = FLUXDNM(1:NPTS)
      MEANRAD(1:NPTS) = MEANRADM(1:NPTS)
      RADOUT(:,:,:) = RADOUTM(:,:,:)
    ENDIF
  ENDDO

!  ENDDO

  CALL OUTPUT_RESULTS (NLAY, NPTS, NMU, NPHI,  NG, &
                     PROPFILE, CKDFILE, SRCTYPE, SOLARFLUX, SOLARMU, &
                     SKYRAD, SFCTYPE, SFCTEMP, SFCPARMS, WAVELEN, WAVENO, & 
                     SPLITACC, SOLACC, MAXITER, NUMITER, &
                     IXP, HEIGHTG, FLUXUP, FLUXDN, MEANRAD, FLUXPARM, &
                     NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
                     'F', FLUXFILE)

  IF (NZOUT*NMUOUT*NPHIOUT > 0) THEN
    CALL OUTPUT_RESULTS (NLAY, NPTS, NMU, NPHI,  NG, &
                     PROPFILE, CKDFILE, SRCTYPE, SOLARFLUX, SOLARMU, &
                     SKYRAD, SFCTYPE, SFCTEMP, SFCPARMS, WAVELEN, WAVENO, & 
                     SPLITACC, SOLACC, MAXITER, NUMITER, &
                     IXP, HEIGHTG, FLUXUP, FLUXDN, MEANRAD, FLUXPARM, &
                     NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
                     'R', RADFILE)
  ENDIF

  DEALLOCATE (HEIGHTP, TEMPP, TAUP, ALBEDOP, TAUPM, ALBEDOPM, NLEGP, LEGENP)
  DEALLOCATE (IXP, HEIGHTG, TAUG, FLUXUP, FLUXDN, MEANRAD, RADOUT)
  DEALLOCATE (FLUXUPM, FLUXDNM, MEANRADM, RADOUTM)
  DEALLOCATE (DELG, ZCKD, KABS)
END




SUBROUTINE SOLVE_SHDOMPP (NLAY, MAXLEG, HEIGHTP, TEMPP, TAUP, ALBEDOP, &
                          NLEGP, LEGENP, SRCTYPE, SOLARFLUX, SOLARMU, &
                          WAVELEN, WAVENO, SKYRAD, SFCTYPE, SFCPARMS, SFCTEMP, &
                          NMU, NPHI, MAXIG, MAXDELTAU, SPLITACC, SOLACC, &
                          MAXITER, ACCELFLAG, PRNT,   SOLCRIT, ITER, &
                          NPTS, HEIGHTG, TAUG, IXP, FLUXUP, FLUXDN, MEANRAD, &
                          NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                          RADOUT)
 ! Performs the SHDOMPP solution procedure to calculate plane-parallel
 ! monochromatic, unpolarized radiative transfer.  The optical properties,
 ! uniform in each layer, are input.  The output are the upwelling and 
 ! downwelling fluxes and mean radiance at the adaptive grid optical depth
 ! levels (TAUG), and the radiances in the NMUOUT*NPHIOUT specified
 ! directions at each of NZOUT radiance output levels.
 !
 ! Inputs:
 !  NLAY      number of input property layers (TAUP, ALBEDOP, etc.)
 !  MAXLEG    maximum Legendre series order for phase functions
 !  HEIGHTP   heights (e.g. km) of input layer interfaces (NLAY+1 values)
 !  TEMPP     temperatures (K) of input layer interfaces (NLAY+1 values)
 !              temperature is linearly interpolated in tau within each layer
 !  TAUP      optical depths of the input layers (each layer from top down)
 !  ALBEDOP   single scattering albedo of the inputs layers
 !  NLEGP     order of each layer Legendre series phase function
 !  LEGENP    phase function Legendre series coefficients for each layer
 !  SRCTYPE   'S' for solar source, 'T' for thermal source, 'B' for both
 !  SOLARFLUX solar flux on a *horizontal* surface
 !  SOLARMU   negative of cosine of solar zenith angle (< 0)
 !              direct beam is assumed to be going in phi=0 azimuth
 !  WAVELEN   wavelength (microns) for thermal emission
 !  WAVENO(2) starting and ending wavenumber (cm^-1) for thermal (if WAVELEN=0)
 !  SKYRAD    top of atmosphere downward sky radiance (or temp for SRCTYPE='T')
 !  SFCTYPE   surface type (currently 'L' for Lambertian, 'R' for RPV
 !  SFCPARMS  array of surface parameters (first is albedo for Lambertian)
 !  SFCTEMP   surface temperature (K)
 !  NMU       number of discrete ordinate zenith angles in both hemispheres
 !  NPHI      number of discrete ordinate azimuth angles in 2\pi
 !  MAXIG     maximum possible size of adaptive grid arrays
 !  MAXDELTAU max sublayer scattering optical depth in adaptive grid
 !  SPLITACC  adaptive layer splitting accuracy
 !  SOLACC    desired solution accuracy (try 1E-4 or less)
 !  MAXITER   maximum number of iterations allowed
 !  ACCELFLAG logical flag: true to perform iteration acceleration
 !  PRNT      logical flat: true to printer iteration information
 !  NZOUT     number of radiance output height levels
 !  ZOUT      height levels of radiance output (set to closest adaptive level)
 !  NMUOUT    number of cosine zenith angles for radiance output
 !  MUOUT     cosine zenith angles for radiance output
 !  NPHIOUT   number of azimuth angles for radiance output
 !  PHIOUT    azimuth angles (degrees) for radiance output
 !
 ! Outputs: 
 !  SOLCRIT   actual solution accuracy achieved
 !  ITER      actual number of iterations
 !  NPTS      number of adaptive grid points
 !  HEIGHTG   heights of the adaptive grid levels (same units as HEIGHTP)
 !  TAUG      sorted grid point optical depth array for the adaptive grid
 !  IXP       index of starting and ending point in TAUG for each input layer
 !  FLUXUP    upward hemispheric flux at each adaptive grid level
 !  FLUXDN    downward hemispheric fluxes (includes direct solar) at each level
 !  MEANRAD   mean radiance (includes direct solar) at each level
 !  RADOUT    array of output radiances
 ! 
 !  TAUG(IXP(1,layer)) is the optical depth at the top of the input layer
 !  TAUG(IXP(2,layer)) is the optical depth at the bottom of the input layer
 !  FLUXUP(IXP(1,layer)) is the upwelling flux at the top of the layer
 !  FLUXDN(IXP(2,layer)) is the downwelling flux at the bottom of the layer
 !
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NLAY, MAXLEG, NLEGP(NLAY), MAXIG
  INTEGER, INTENT(IN)  :: NMU, NPHI, MAXITER
  INTEGER, INTENT(IN)  :: NZOUT, NMUOUT, NPHIOUT
  LOGICAL, INTENT(IN)  :: ACCELFLAG, PRNT
  REAL,    INTENT(IN)  :: TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG,NLAY)
  REAL,    INTENT(IN)  :: HEIGHTP(NLAY+1), TEMPP(NLAY+1)
  REAL,    INTENT(IN)  :: SOLARFLUX, SOLARMU, WAVELEN, WAVENO(2), SKYRAD
  REAL,    INTENT(IN)  :: SFCPARMS(*), SFCTEMP
  REAL,    INTENT(IN)  :: MAXDELTAU, SPLITACC, SOLACC
  REAL,    INTENT(IN)  :: ZOUT(NZOUT), MUOUT(NMUOUT), PHIOUT(NPHIOUT)
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE, SFCTYPE
  INTEGER, INTENT(OUT) :: ITER, NPTS, IXP(2,NLAY)
  REAL,    INTENT(OUT) :: SOLCRIT
  REAL,    INTENT(OUT) :: HEIGHTG(MAXIG), TAUG(MAXIG)
  REAL,    INTENT(OUT) :: FLUXUP(MAXIG), FLUXDN(MAXIG), MEANRAD(MAXIG)
  REAL,    INTENT(OUT) :: RADOUT(NZOUT,NMUOUT,NPHIOUT)

  INTEGER :: ML, MM, NLM, NPHI0MAX
  INTEGER :: ORDINATESET, OLDNPTS, PREVNPTS, LASTNPTS, NANG, I
  LOGICAL, PARAMETER :: DELTAM=.TRUE.
  LOGICAL :: FFTFLAG(NMU)
  REAL,    PARAMETER :: MINDELTAUSPLIT=1.0E-3
  REAL    :: DELJDOT, DELJOLD, DELJNEW
  REAL    :: SOLARAZ, PI, GNDALB, DIRFLUX
  REAL    :: TRUESOLCRIT, CURSPLITACC
  INTEGER, ALLOCATABLE :: IXG(:), NPHI0(:)
  REAL, ALLOCATABLE :: TAUPSC(:), ALBEDOSC(:), LEGENSC(:,:)
  REAL, ALLOCATABLE :: TAUGSC(:), RADIANCE(:,:,:), SOURCE(:,:)
  REAL, ALLOCATABLE :: MU(:), WTDO(:,:), PHI(:,:), WTMU(:), YLMSUN(:)
  REAL, ALLOCATABLE :: CMU1(:,:), CMU2(:,:), CPHI1(:,:,:), CPHI2(:,:,:)
  REAL, ALLOCATABLE :: WPHISAVE(:,:), SRCINT(:,:,:)
  REAL, ALLOCATABLE :: CONSTSOURCE(:,:), SHRADIANCE(:,:), DOSOURCE(:,:,:)
  REAL, ALLOCATABLE :: DELSOURCE(:,:)

  ! Make NPHI0MAX, ML, MM, and NLM from NMU and NPHI
  IF (NMU /= MAX(2,2*INT((NMU+1)/2)) ) STOP 'SOLVE_SHDOMPP: Bad NMU'
  IF (NPHI < 1) STOP 'SOLVE_SHDOMPP: NPHI must be greater than 0'
  NPHI0MAX = INT((NPHI+2)/2)
  ML = NMU-1
  MM = MAX(0,INT(NPHI/2)-1)
  NLM = (MM+1)*(ML+1) - (MM*(MM+1))/2

  ALLOCATE (TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML,NLAY))
  ALLOCATE (TAUGSC(MAXIG), IXG(MAXIG))
  ALLOCATE (MU(NMU), PHI(NPHI0MAX,NMU), NPHI0(NMU), WTDO(NMU,NPHI0MAX))
  ALLOCATE (YLMSUN(NLM), WTMU(NMU), CMU1(NLM,NMU), CMU2(NMU,NLM))
  ALLOCATE (CPHI1(0:16,32,NMU), CPHI2(32,0:16,NMU), WPHISAVE(3*NPHI0MAX+15,NMU))
  ALLOCATE (SOURCE(NLM,MAXIG), RADIANCE(NPHI0MAX,NMU,MAXIG))
  ALLOCATE (CONSTSOURCE(NLM,MAXIG), DELSOURCE(NLM,MAXIG))
  ALLOCATE (SRCINT(0:4,NMU,MAXIG))


  ! Set up some things before solution loop
 
  ! If Delta-M then scale the optical depth, albedo, and Legendre terms.
  IF (DELTAM) THEN
    CALL DELTA_SCALE (NLAY, ML, TAUP, ALBEDOP, MAXLEG, NLEGP, LEGENP, &
                      TAUPSC, ALBEDOSC, LEGENSC)
  ELSE
    TAUPSC(:) = TAUP(:)
    ALBEDOSC(:) = ALBEDOP(:)
    LEGENSC(0:ML,:) = LEGENP(0:ML,:)
  ENDIF

  CALL MAKE_TAU_GRID (NLAY, TAUPSC, ALBEDOSC, MAXDELTAU, MINDELTAUSPLIT,&
                      NPTS, TAUGSC, IXP, IXG)

  ! Precompute Ylm's for solar direction
  IF (SRCTYPE /= 'T') THEN
    SOLARAZ = 0.0
    CALL YLMALL (SOLARMU, SOLARAZ, ML, MM, 1, YLMSUN)
  ENDIF
  ! Prepare the constant part of the level source function, which are
  !   the Planck blackbody source and the solar attenuated single scattered
  !   radiance, as a function of level.
  CONSTSOURCE(1,1:NPTS) = -1.0
  CALL PREPARE_SOURCE (NLAY, TAUPSC, ALBEDOSC, TEMPP, LEGENSC, &
                       SRCTYPE, ML, MM, NLM, NPTS, TAUGSC, IXP, IXG, &
                       SOLARFLUX, SOLARMU, YLMSUN, WAVELEN, WAVENO, &
                       CONSTSOURCE)

  ! Make the discrete ordinates (angles) 
  !   (set 2 is reduced gaussian, 3 is reduced double gauss)
  ORDINATESET = 3
  CALL MAKE_ANGLE_SET (NMU, NPHI, 1, NPHI0MAX, ORDINATESET, &
                       NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

  ! Make the Ylm transform coefficients
  CALL MAKE_SH_DO_COEF (ML, MM, NLM, NMU, NPHI0, &
                        NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG, &
                        CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)

  ! Initialize the radiance on the base grid using Eddington 
  !   two-stream plane-parallel
  IF (SFCTYPE == 'L') THEN
    GNDALB = SFCPARMS(1)
  ELSE
    GNDALB = 0.2
  ENDIF
  ALLOCATE (SHRADIANCE(NLM,NPTS))
  CALL INIT_RADIANCE (NLAY, ML, TAUPSC, ALBEDOSC, LEGENSC, TEMPP, &
                      SRCTYPE, SOLARFLUX, SOLARMU, GNDALB, SFCTEMP, & 
                      SKYRAD, WAVELEN,  WAVENO, NPTS, TAUGSC, IXP, NLM, &
                      SHRADIANCE)

  ! Initialize the SH source function from the SH radiance field
  OLDNPTS = 0
  SOURCE(:,:) = 0.0
  CALL COMPUTE_SOURCE (ML,MM, NLM, NLAY, NPTS, IXP, IXG, &
                       ALBEDOSC, LEGENSC, &
                       SHRADIANCE, CONSTSOURCE, SOURCE, &
                       OLDNPTS, DELSOURCE, &
                       DELJDOT, DELJOLD, DELJNEW, SOLCRIT)
  DEALLOCATE (SHRADIANCE)


  SOLCRIT = 1.0
  TRUESOLCRIT = 1.0
  CURSPLITACC = 1.0
  LASTNPTS = -1
  OLDNPTS = NPTS
  RADIANCE(:,:,:) = 0.0
  DELSOURCE(:,:) = 0.0
  ITER = 0
  IF (MAXITER <= ITER) RETURN

  IF (PRNT) WRITE (*,*) '! Iter Log(Sol)  Log(True)  Npoints'
 
  ! Main solution loop
    ! Iterate until the iterations exceed the limit, or the desired
    ! solution criterion is reached 
  DO WHILE (ITER < MAXITER .AND. (SOLCRIT > SOLACC .OR. NPTS /= OLDNPTS))
    ITER = ITER + 1

    ! Transform the source from spherical harmonics to discrete ordinates
    ALLOCATE (DOSOURCE(NPHI0MAX,NMU,MAXIG))
    CALL SH_TO_DO (NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                   FFTFLAG, CMU1, CPHI1, WPHISAVE, SOURCE, DOSOURCE)
      ! Enforce the discrete ordinate source function to be non-negative
    DOSOURCE(:,:,:) = MAX(0.0,DOSOURCE(:,:,:))


    ! Keep on doing integrations as long as new levels are added
    OLDNPTS = NPTS
    PREVNPTS = 0
    DO WHILE (NPTS > PREVNPTS)
      ! Integrate the source function along discrete ordinates to
      ! compute the radiance field.  
      CALL PATH_INTEGRATION (NLAY, NPTS, IXP, IXG, TAUGSC, &
                             NMU, NPHI0MAX, NPHI0, &
                             SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                             SFCTYPE, SFCPARMS, SFCTEMP, WAVELEN, WAVENO, &
                             MU, PHI, WTDO, SRCINT, LASTNPTS, &
                             FLUXUP, FLUXDN, MEANRAD, DOSOURCE, RADIANCE)
      PREVNPTS = NPTS
      ! Do the adaptive grid cell splitting stuff if desired
      CALL SPLIT_CELLS (NLAY, NPTS, MAXIG, IXP, IXG, TAUGSC, &
                        NMU, NPHI0MAX, NPHI0, MU, DOSOURCE, RADIANCE, &
                        SPLITACC, TRUESOLCRIT, CURSPLITACC)
      IF (NPTS > PREVNPTS) THEN
        CONSTSOURCE(1,PREVNPTS+1:NPTS) = -1.0
        CALL PREPARE_SOURCE (NLAY, TAUPSC, ALBEDOSC, TEMPP, LEGENSC, &
                             SRCTYPE, ML, MM, NLM, NPTS, TAUGSC, IXP, IXG, &
                             SOLARFLUX, SOLARMU, YLMSUN, WAVELEN, WAVENO, &
                             CONSTSOURCE)
      ENDIF
    ENDDO
    DEALLOCATE (DOSOURCE)

    ! Transform the radiance from discrete ordinates to spherical harmonics
    ALLOCATE (SHRADIANCE(NLM,NPTS))
    CALL DO_TO_SH (NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                   FFTFLAG, CMU2, CPHI2, WPHISAVE, RADIANCE, SHRADIANCE)

    ! Calculate the source function from the radiance in spherical harmonics.
    !   Also computes the solution criterion and the series acceleration
    !   stuff (DELJDOT, DELJOLD, DELJNEW, DELSOURCE).
    CALL COMPUTE_SOURCE (ML,MM, NLM, NLAY, NPTS, IXP, IXG, &
                         ALBEDOSC, LEGENSC, &
                         SHRADIANCE, CONSTSOURCE, SOURCE, &
                         OLDNPTS, DELSOURCE, &
                         DELJDOT, DELJOLD, DELJNEW, SOLCRIT)
    DEALLOCATE (SHRADIANCE)

    ! Accelerate the convergence of the source function vector
    IF (ACCELFLAG) THEN
      CALL ACCELERATE_SOLUTION (NLM, NPTS, OLDNPTS, SOURCE, DELSOURCE, &
                                DELJDOT, DELJOLD, DELJNEW)
    ENDIF
    IF (DELJNEW < DELJOLD) &
        TRUESOLCRIT = SOLCRIT/(1.0-SQRT(DELJNEW/DELJOLD))

    ! Print out the iteration, log solution criterion, and number of points
    IF (PRNT) WRITE (*,'(2X,I4,2F8.3,1X,I6)')  ITER, &
        LOG10(MAX(SOLCRIT,1.0E-20)), LOG10(MAX(TRUESOLCRIT,1.0E-20)), NPTS
  ENDDO

  IF (PRNT) WRITE (*,'(1X,A,I6,A,F9.6)') &
              '! Iterations: ', ITER, '     Final Criterion: ', SOLCRIT


  ! Add the direct solar component to the downwelling flux and mean radiance
  IF (SRCTYPE /= 'T') THEN
    PI = ACOS(-1.0)
    DO I = 1, NPTS
      DIRFLUX = SOLARFLUX*EXP(-TAUGSC(I)/ABS(SOLARMU))
      FLUXDN(I) = FLUXDN(I) + DIRFLUX
      MEANRAD(I) = MEANRAD(I) + DIRFLUX/(4*PI*ABS(SOLARMU))
    ENDDO
  ENDIF

  ! Make an unscaled optical depth grid for the output file
  IF (DELTAM) THEN
    CALL UNDELTA_SCALE (NLAY, TAUPSC, TAUP, IXP, NPTS, TAUGSC, TAUG)
  ELSE
    TAUG(:) = TAUGSC(:)
  ENDIF
  CALL OUTPUT_HEIGHT_GRID (NLAY, TAUP, IXP, NPTS, TAUG, HEIGHTP, HEIGHTG)

  IF (NZOUT*NMUOUT*NPHIOUT > 0) THEN
    ! Compute the radiance output and put in array
    CALL COMPUTE_RADIANCE (NLAY, ALBEDOSC, LEGENSC, &
                         MAXLEG, NLEGP, LEGENP,  DELTAM, &
                         NPTS, HEIGHTG, TAUG, TAUGSC, IXP, IXG, &
                         ML, MM, NLM, NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                         SRCTYPE, SOLARFLUX, SOLARMU, &
                         SKYRAD, WAVELEN, WAVENO, SFCTEMP, &
                         SFCTYPE, SFCPARMS, SOURCE, RADIANCE, &
                         NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                         RADOUT)
  ENDIF

  DEALLOCATE (TAUPSC, ALBEDOSC, LEGENSC, TAUGSC)
  DEALLOCATE (MU, PHI, NPHI0, WTDO, YLMSUN, WTMU)
  DEALLOCATE (CPHI1, CPHI2, WPHISAVE, CMU1, CMU2)
  DEALLOCATE (SOURCE, RADIANCE, CONSTSOURCE, DELSOURCE, SRCINT)
END SUBROUTINE SOLVE_SHDOMPP

 


SUBROUTINE DELTA_SCALE (NLAY, ML, TAUP, ALBEDOP, MAXLEG, NLEGP, LEGENP, &
                        TAUPSC, ALBEDOSC, LEGENSC)
 ! Delta scales the layer optical depth, albedo, and Legendre terms;
 ! only the 0 to ML LEGEN terms are scaled.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, ML, MAXLEG, NLEGP(NLAY)
  REAL,    INTENT(IN) ::  TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG,NLAY)
  REAL,    INTENT(OUT) ::  TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML,NLAY)
  INTEGER :: I, L
  REAL    :: F

  DO I = 1, NLAY
    IF (MAXLEG < ML+1) THEN
      F = 0.0
      LEGENSC(:,I) = 0.0
    ELSE
      F = LEGENP(ML+1,I)/(2*(ML+1)+1)
    ENDIF
    DO L = 0, MIN(ML,MAXLEG)
      LEGENSC(L,I) = (2*L+1)*(LEGENP(L,I)/(2*L+1) - F)/(1-F)
    ENDDO
    TAUPSC(I) = (1.0-ALBEDOP(I)*F)*TAUP(I)
    ALBEDOSC(I) = (1.0-F)*ALBEDOP(I)/(1.0-ALBEDOP(I)*F)
  ENDDO
END SUBROUTINE DELTA_SCALE


SUBROUTINE UNDELTA_SCALE (NLAY, TAUPSC, TAUP, IXP, NPTS, TAUGSC, TAUG)
 ! Undoes the delta scaling to get the optical depth grid in the 
 ! original scale.  Uses the scaled and unscaled input layer optical depths.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, IXP(2,NLAY), NPTS
  REAL,    INTENT(IN) :: TAUPSC(NLAY), TAUP(NLAY), TAUGSC(NPTS)
  REAL,    INTENT(OUT) :: TAUG(NPTS)
  INTEGER :: LAY, I
  REAL    :: R

  TAUG(1) = 0.0
  DO LAY = 1, NLAY
    IF (TAUPSC(LAY) > 0.0) THEN
      R = TAUP(LAY)/TAUPSC(LAY)
    ELSE
      R = 1.0
    ENDIF
    DO I = IXP(1,LAY), IXP(2,LAY)-1
      TAUG(I+1) = TAUG(I) + R*(TAUGSC(I+1)-TAUGSC(I))
    ENDDO
    IF (LAY < NLAY) TAUG(IXP(1,LAY+1)) = TAUG(IXP(2,LAY))
  ENDDO
END SUBROUTINE UNDELTA_SCALE



SUBROUTINE OUTPUT_HEIGHT_GRID (NLAY, TAUP, IXP, NPTS, TAUG, HEIGHTP, HEIGHTG)
 ! Makes the output height grid (HEIGHTG) by interpolating between the 
 ! levels in the input height grid (HEIGHTP) according to the optical depth
 ! grid (TAUG) assuming constant extinction in each input layer.
 ! The levels in HEIGHTG correspond to the optical depths in TAUG.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, IXP(2,NLAY), NPTS
  REAL,    INTENT(IN) :: TAUP(NLAY), TAUG(NPTS), HEIGHTP(NLAY+1)
  REAL,    INTENT(OUT) :: HEIGHTG(NPTS)
  INTEGER :: LAY, I
  REAL    :: E

  DO LAY = 1, NLAY
    HEIGHTG(IXP(1,LAY)) = HEIGHTP(LAY)
    HEIGHTG(IXP(2,LAY)) = HEIGHTP(LAY+1)
    IF (TAUP(LAY) > 0.0) THEN
      E = (HEIGHTP(LAY)-HEIGHTP(LAY+1))/TAUP(LAY)
      DO I = IXP(1,LAY)+1, IXP(2,LAY)-1
        HEIGHTG(I) = HEIGHTP(LAY) - E*(TAUG(I)-TAUG(IXP(1,LAY)))
      ENDDO
    ELSE
      DO I = IXP(1,LAY)+1, IXP(2,LAY)-1
        HEIGHTG(I) = HEIGHTP(LAY) + (HEIGHTP(LAY+1)-HEIGHTP(LAY))* &
                    (I-IXP(1,LAY))/FLOAT(IXP(2,LAY)-IXP(1,LAY))
      ENDDO
    ENDIF
  ENDDO
END SUBROUTINE OUTPUT_HEIGHT_GRID



SUBROUTINE MAKE_TAU_GRID (NLAY, TAUPSC, ALBEDOSC, MAXDELTAU, MINDELTAUSPLIT, &
                          NPTS, TAUGSC, IXP, IXG)
 ! Makes the initial optical depth grid, TAUGSC, of increasing values
 ! (representing decreasing altitudes in the atmosphere). TAUGSC(1) is zero.
 ! The "base" tau grid is the input property grid boundaries and even 
 ! subdivision of each of these layers so that the sublayers are less 
 ! than MAXDELTAU scattering optical depth.   All layers with scaled
 ! optical depth greater than MINDELTAUSPLIT are split into two sublayers.
 ! The IXP index array points to the top and bottom TAUGSC levels for each 
 ! property layer. The IXG index array points to the source and 
 ! radiance arrays for each TAUGSC level.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NLAY
  REAL,    INTENT(IN)  :: TAUPSC(NLAY), ALBEDOSC(NLAY)
  REAL,    INTENT(IN)  :: MAXDELTAU, MINDELTAUSPLIT
  REAL,    INTENT(OUT) :: TAUGSC(*)
  INTEGER, INTENT(OUT) :: NPTS, IXP(2,NLAY), IXG(*)
  INTEGER :: I, L, N
  REAL    :: TAU

  NPTS = 0
  TAU = 0.0
  DO L = 1, NLAY
    NPTS = NPTS + 1
    TAUGSC(NPTS) = TAU
    IXP(1,L) = NPTS
    IF (TAUPSC(L) > MINDELTAUSPLIT) THEN
      N = MAX(1,INT(ALBEDOSC(L)*TAUPSC(L)/MAXDELTAU))
    ELSE
      N = 0
    ENDIF
    DO I = 1, N
      NPTS = NPTS + 1
      TAUGSC(NPTS) = TAUGSC(IXP(1,L)) + TAUPSC(L)*FLOAT(I)/(N+1)
    ENDDO
    NPTS = NPTS + 1
    IXP(2,L) = NPTS
    TAUGSC(NPTS) = TAUGSC(IXP(1,L)) + TAUPSC(L)
    TAU = TAUGSC(NPTS)
  ENDDO
  IXG(1:NPTS) = (/ (I, I=1,NPTS) /)
END SUBROUTINE MAKE_TAU_GRID




SUBROUTINE PREPARE_SOURCE (NLAY, TAUPSC, ALBEDOSC, TEMPP, LEGENSC, &
                           SRCTYPE, ML, MM, NLM, NPTS, TAUGSC, IXP, IXG, &
                           SOLARFLUX, SOLARMU, YLMSUN, WAVELEN, WAVENO, &
                           CONSTSOURCE)
 ! Prepares the constant part of the source function for each level.
 ! This is either or both thermal emission or the single scattering
 ! solar pseudo-source.  Only those levels (points) with CONSTSOURCE(1,:)
 ! less than zero are changed, which allows this to be run for each time
 ! adaptive grid points are added.
 ! The thermal emission source is \sqrt{4\pi}(1-\omega)B(T), where the
 ! temperature T is assumed to vary linearly in optical depth.
 ! The solar source is 
 !   S_0/\mu_0 e^{-\tau/\mu_0} Y_{lm}(\mu_0,\phi_0) {\omega \chi_l \over 2l+1}
 ! where S_0 is the solar flux on a horizontal plane, \chi_l are the Legendre
 ! coefficients, and \tau is the optical depth at the level.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, ML, MM, NLM, NPTS
  INTEGER, INTENT(IN) :: IXP(2,NLAY), IXG(NPTS)
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, YLMSUN(NLM), WAVELEN, WAVENO(2)
  REAL,    INTENT(IN) :: TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML,NLAY)
  REAL,    INTENT(IN) :: TEMPP(NLAY+1), TAUGSC(NPTS)
  REAL,    INTENT(INOUT) :: CONSTSOURCE(NLM,NPTS)
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE
  INTEGER :: I, J, K, L, LAY, ME, M, LOFJ(NLM)
  REAL    :: S0, C, TAU1, U, TEMP, BB, D, PLANCKFUNC

  S0 = SOLARFLUX/ABS(SOLARMU)
  C = SQRT(4.0*ACOS(-1.0))
  ! Make the l index as a function of SH term (J)
  J = 0
  DO L = 0, ML
    ME = MIN(L,MM)
    DO M = 0, ME
      J = J + 1
      LOFJ(J) = L
    ENDDO
  ENDDO

  DO LAY = 1, NLAY
    TAU1 = TAUGSC(IXP(1,LAY))
    DO I = IXP(1,LAY), IXP(2,LAY)
      K = IXG(I)
      IF (CONSTSOURCE(1,K) < 0) THEN
        CONSTSOURCE(:,K) = 0.0
        IF (SRCTYPE /= 'S') THEN
          IF (TAUPSC(LAY) > 0.0) THEN
            U = (TAUGSC(I)-TAU1)/TAUPSC(LAY)
          ELSE
            U = (I-IXP(1,LAY))/FLOAT(IXP(2,LAY)-IXP(1,LAY))
          ENDIF
          TEMP = TEMPP(LAY) + U*(TEMPP(LAY+1)-TEMPP(LAY))
          CONSTSOURCE(1,K) = C*(1.0-ALBEDOSC(LAY)) &
                            *PLANCKFUNC(WAVELEN,WAVENO,TEMP)
        ENDIF
        IF (SRCTYPE /= 'T') THEN
          D = S0*EXP(-TAUGSC(I)/ABS(SOLARMU))
          DO J = 1, NLM
            L = LOFJ(J)
            CONSTSOURCE(J,K) = CONSTSOURCE(J,K) &
                           + D*YLMSUN(J)*ALBEDOSC(LAY)*LEGENSC(L,LAY)/(2*L+1)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE PREPARE_SOURCE


FUNCTION PLANCKFUNC (WAVELEN, WAVENO, TEMP)
 ! Calculates the Planck blackbody radiance. If WAVELEN<=0 then doing a 
 ! spectral band integration and the Planck blackbody radiance in 
 ! [Watts /(meter^2 ster)] over a wavenumber range [cm^-1] (WAVENO) is 
 ! returned.  Otherwise, the Planck blackbody radiance in 
 ! [Watts /(meter^2 ster micron)] for a temperature in [Kelvins] at 
 ! a wavelength (WAVELEN) in [microns] is returned.
  IMPLICIT NONE
  REAL, INTENT(IN) :: WAVELEN, WAVENO(2), TEMP
  REAL :: PLANCKFUNC
  DOUBLE PRECISION :: X1, X2, F

  IF (TEMP > 0.0) THEN
    IF (WAVELEN>0) THEN
      PLANCKFUNC = 1.1911E8 / WAVELEN**5 / (EXP(1.4388E4/(WAVELEN*TEMP)) - 1)
    ELSE
      X1 = 1.4388D0*WAVENO(1)/TEMP
      X2 = 1.4388D0*WAVENO(2)/TEMP
      CALL INTEGRATE_PLANCK (X1, X2, F)
      PLANCKFUNC = 1.1911D-8*(TEMP/1.4388D0)**4 *F
    ENDIF
  ELSE
    PLANCKFUNC = 0.0
  ENDIF
END FUNCTION PLANCKFUNC



SUBROUTINE MAKE_ANGLE_SET (NMU, NPHI, NCS, NPHI0MAX, ITYPE, &
                           NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)
 ! Make the set of angles for the discrete space representation.
 ! The number of mu's (cosine zenith angles) and maximum number of
 ! phi's (azimuth angles) is input.  The actual number of azimuthal 
 ! angles for each mu is output (NPHI0).  There are three types of
 ! discrete ordinate sets: ITYPE=1 is a gaussian grid, 2 is a reduced
 ! gaussian grid, and 3 is a reduced double gaussian set.
 ! If NCS=1 then only about half the azimuthal angles (from 0 to pi) 
 ! are used because the radiance is even in phi (cosine terms).  
 ! The output is the NMU mu values, the NPHI0 phi values for each mu,
 ! and the integration weight for each ordinate. The first NMU/2 mu 
 ! angles are the downwelling (mu<0) angles.  Also output are the
 ! the flags for doing an azimuthal FFT for each mu and the total
 ! number of angles (NANG).
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NMU, NPHI, NCS, NPHI0MAX, ITYPE
  INTEGER, INTENT(OUT) :: NPHI0(NMU), NANG
  LOGICAL, INTENT(OUT) :: FFTFLAG(NMU)
  REAL,    INTENT(OUT) :: MU(NMU), WTMU(NMU)
  REAL,    INTENT(OUT) :: PHI(NPHI0MAX,NMU), WTDO(NMU,NPHI0MAX)
  INTEGER :: J, K, MM
  REAL    :: DELPHI
  INTEGER, PARAMETER :: MAXNPHI=1024
  INTEGER :: GOODNFFT(MAXNPHI)

  CALL GOOD_FFT_LENGTHS (MAXNPHI, GOODNFFT)
  GOODNFFT(1:14) = (/ (K, K=1,14) /)

  MM = MAX(0,INT(NPHI/2)-1)
  NANG = 0
  IF (ITYPE .GE. 1 .OR. ITYPE .LE. 3) THEN
    IF (ITYPE .LE. 2) THEN
      CALL GAUSQUADS (NMU, MU, WTMU)
    ELSE IF (ITYPE .EQ. 3) THEN
      CALL DGAUSQUADS (NMU, MU, WTMU)
    ENDIF
    DO J = 1, NMU
      IF (NCS .EQ. 1) THEN
        IF (NPHI0MAX .NE. INT((NPHI+2)/2)) STOP 'MAKE_ANGLE_SET: bad NPHI0MAX'
        IF (ITYPE .EQ. 1) THEN
          NPHI0(J) = NPHI0MAX
        ELSE
        ! For the reduced gaussian set, make the smaller NPHI0 values
        ! still be good for the FFT (but don't let NPHI0MAX be exceeded).
          IF (NPHI0MAX .GT. MAXNPHI) STOP 'MAKE_ANGLE_SET: exceeded GOODNFFT'
          NPHI0(J) = INT(0.9+1+(NPHI0MAX-1)*SQRT(1-MU(J)**2))
          IF (NPHI0(J) .GT. 1) NPHI0(J) = MIN(NPHI0MAX,GOODNFFT(NPHI0(J)-1)+1)
        ENDIF
        ! Compute the azimuth angles and weights
        DELPHI = ACOS(-1.0)/MAX(1,NPHI0(J)-1)
        DO K = 1, NPHI0(J)
          PHI(K,J) = (K-1)*DELPHI
          IF ((K.EQ.1 .OR. K.EQ.NPHI0(J)) .AND. NPHI0(J).NE.1) THEN
            WTDO(J,K) = DELPHI*WTMU(J)
          ELSE
            WTDO(J,K) = 2.0*DELPHI*WTMU(J)
          ENDIF
        ENDDO
      ELSE
        IF (NPHI0MAX .NE. NPHI) STOP 'MAKE_ANGLE_SET: bad NPHI0MAX'
        IF (ITYPE .EQ. 1) THEN
          NPHI0(J) = NPHI0MAX
        ELSE
          IF (NPHI0MAX .GT. MAXNPHI) STOP 'MAKE_ANGLE_SET: exceeded GOODNFFT'
          NPHI0(J) = INT(0.9+NPHI0MAX*SQRT(1-MU(J)**2))
          NPHI0(J) = MIN(NPHI0MAX,GOODNFFT(NPHI0(J)))
        ENDIF
        DELPHI = 2.0*ACOS(-1.0)/NPHI0(J)
        DO K = 1, NPHI0(J)
          PHI(K,J) = (K-1)*DELPHI
          WTDO(J,K) = DELPHI*WTMU(J)
        ENDDO
      ENDIF
      NANG = NANG + NPHI0(J)
      FFTFLAG(J) = (NPHI0(J) .GT. 14) .OR. (MM .GT. 16)
    ENDDO
        
  ELSE
    STOP 'MAKE_ANGLE_SET: invalid discrete ordinate type'
  ENDIF
END SUBROUTINE MAKE_ANGLE_SET



SUBROUTINE GOOD_FFT_LENGTHS (NMAX, GOODNFFT)
 ! Returns a list in GOODNFFT(I) that rounds up to the nearest I that
 ! has only divisors of 2 or 3.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NMAX
  INTEGER, INTENT(OUT) :: GOODNFFT(NMAX)
  INTEGER :: I, J, ILAST
  INTEGER, PARAMETER :: PRIMES(310) = &
          (/  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, &
             31, 37, 41, 43, 47, 53, 59, 61, 67, 71, &
             73, 79, 83, 89, 97,101,103,107,109,113, &
            127,131,137,139,149,151,157,163,167,173, &
            179,181,191,193,197,199,211,223,227,229, &
            233,239,241,251,257,263,269,271,277,281, &
            283,293,307,311,313,317,331,337,347,349, &
            353,359,367,373,379,383,389,397,401,409, &
            419,421,431,433,439,443,449,457,461,463, &
            467,479,487,491,499,503,509,521,523,541, &
            547,557,563,569,571,577,587,593,599,601, &
            607,613,617,619,631,641,643,647,653,659, &
            661,673,677,683,691,701,709,719,727,733, &
            739,743,751,757,761,769,773,787,797,809, &
            811,821,823,827,829,839,853,857,859,863, &
            877,881,883,887,907,911,919,929,937,941, &
            947,953,967,971,977,983,991,997,1009,1013, &
            1019,1021,1031,1033,1039,1049,1051,1061,1063,1069, &
            1087,1091,1093,1097,1103,1109,1117,1123,1129,1151, &
            1153,1163,1171,1181,1187,1193,1201,1213,1217,1223, &
            1229,1231,1237,1249,1259,1277,1279,1283,1289,1291, &
            1297,1301,1303,1307,1319,1321,1327,1361,1367,1373, &
            1381,1399,1409,1423,1427,1429,1433,1439,1447,1451, &
            1453,1459,1471,1481,1483,1487,1489,1493,1499,1511, &
            1523,1531,1543,1549,1553,1559,1567,1571,1579,1583, &
            1597,1601,1607,1609,1613,1619,1621,1627,1637,1657, &
            1663,1667,1669,1693,1697,1699,1709,1721,1723,1733, &
            1741,1747,1753,1759,1777,1783,1787,1789,1801,1811, &
            1823,1831,1847,1861,1867,1871,1873,1877,1879,1889, &
            1901,1907,1913,1931,1933,1949,1951,1973,1979,1987, &
            1993,1997,1999,2003,2011,2017,2027,2029,2039,2053 /)

  GOODNFFT(:) = (/ (I, I=1,NMAX) /)
  DO I = 3, SIZE(PRIMES)
    DO J = PRIMES(I), NMAX, PRIMES(I)
      GOODNFFT(J) = 0
    ENDDO
  ENDDO
  ILAST = NMAX
  IF (GOODNFFT(ILAST) == 0)  STOP 'GOOD_FFT_LENGTHS: NMAX not good'
  DO I = NMAX, 1, -1
    IF (GOODNFFT(I) == 0) THEN
      GOODNFFT(I) = ILAST 
    ELSE
      ILAST = I
    ENDIF
  ENDDO
END SUBROUTINE GOOD_FFT_LENGTHS




SUBROUTINE MAKE_SH_DO_COEF (ML, MM, NLM, NMU, NPHI0, &
                            NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG, &
                            CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)
 ! Makes the transformation coefficients for the spherical harmonic
 ! transform.  The main coefficients are output in four arrays: 1 is for
 ! the SH_TO_DO forward transform, 2 is for the DO_TO_SH back transform,
 ! which contains the discrete angle integration weights.
 ! The coefficients are divided up into the mu dependent set CMUn 
 ! (function of l, m, mu_j), and the phi dependent set CPHIn 
 ! (function of m, phi_k) for each mu_j.
 ! The FFTPACK phase coefficients for the FFT in azimuth are also 
 ! output in WPHISAVE.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ML, MM, NLM, NMU, NPHI0(NMU), NPHI0MAX
  LOGICAL, INTENT(IN) :: FFTFLAG(NMU)
  REAL,    INTENT(IN) :: MU(NMU), WTMU(NMU)
  REAL,    INTENT(IN) :: PHI(NPHI0MAX,NMU), WTDO(NMU,NPHI0MAX)
  REAL,    INTENT(OUT) :: CMU1(NLM,NMU), CMU2(NMU,NLM)
  REAL,    INTENT(OUT) :: CPHI1(0:16,32,NMU), CPHI2(32,0:16,NMU)
  REAL,    INTENT(OUT) :: WPHISAVE(3*NPHI0MAX+15,NMU)
  INTEGER :: I, K, M
  REAL    :: W

  ! Make the to and from associate Legendre coefficients       
  DO I = 1, NMU
    CALL YLMALL (MU(I), 0.0, ML, MM, -1, CMU1(1,I))
  ENDDO
  DO I = 1, NMU
    CMU2(I,1:NLM) = CMU1(1:NLM,I)*WTMU(I)
  ENDDO

  ! Make the to and from Fourier coefficients for each mu
  DO I = 1, NMU
    IF (.NOT. FFTFLAG(I)) THEN
      ! If not doing an FFT for this mu then make the DFT coefficient
      W = 1.0/WTMU(I)
      DO K = 1, NPHI0(I)
        IF (NPHI0(I) .GT. 32 .OR. MM .GT. 16)  &
          STOP 'MAKE_SH_DO_COEF: Fourier coefficients array exceeded'
        DO M = 0, MM
          CPHI1(M,K,I)  = COS(M*PHI(K,I))
          CPHI2(K,M,I)  = CPHI1(M,K,I)*WTDO(I,K)*W
        ENDDO
      ENDDO
    ELSE
      ! Precompute the phase factors for the FFTs
      CALL COSTI (NPHI0(I),WPHISAVE(1,I))
    ENDIF
  ENDDO
END SUBROUTINE MAKE_SH_DO_COEF
 




SUBROUTINE SH_TO_DO (NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                     FFTFLAG, CMU1, CPHI1, WSAVE, INDATA, OUTDATA)
 ! Transforms the input data from spherical harmonic space to discrete 
 ! ordinate space for all points in the array.  Successive transforms are 
 ! done in zenith and then azimuth angle.  The Fourier transform in 
 ! azimuth is done using a matrix multiply DFT for smaller number of 
 ! angles (NPHI0*(MM+1)<160) and with an FFT for larger number of 
 ! angles.  The FFT is used for any number of angles above the limit, 
 ! but is much more efficient for NPHI0 whose factors are small primes
 ! (powers of 2 are best).  The number of floating point operations for 
 ! the DFT is 2*Nmu*Nphi0*Nm. The number of floating point operations for 
 ! the zenith angle transform is 2*Nmu*Nlm.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0(NMU)
  LOGICAL, INTENT(IN) :: FFTFLAG(NMU)
  REAL,    INTENT(IN) :: INDATA(NLM,NPTS)
  REAL,    INTENT(IN) :: CMU1(NLM,NMU), CPHI1(0:16,32,NMU)
  REAL,    INTENT(IN) :: WSAVE(3*NPHI0MAX+15,NMU)
  REAL,    INTENT(OUT) :: OUTDATA(NPHI0MAX,NMU,NPTS)
  INTEGER :: I, J, L, M, ME, IMU, IPHI
  INTEGER :: MOFJ(NLM)
  REAL    :: S(0:MM)

  ! Make the M for each J array
  J = 0
  DO L = 0, ML
    DO M = 0, MIN(L,MM)
      J = J + 1
      MOFJ(J) = M
    ENDDO
  ENDDO

  ! Do the transform for each grid point
  DO I = 1, NPTS
    DO IMU = 1, NMU
      ! Figure the max Fourier azimuthal mode we can do for this Nphi0
      ME = MAX(0,MIN(NPHI0(IMU)-2,MOFJ(NLM)))
      ! First do Legendre transform by summing over l for each m.
      S(:) = 0.0
      DO J = 1, NLM
        M = MOFJ(J)
        S(M) = S(M) + CMU1(J,IMU)*INDATA(J,I)
      ENDDO
      ! Then do the Fourier transform from m to phi for each mu
      IF (FFTFLAG(IMU)) THEN
        ! For cosine modes only use a cosine FFT
        OUTDATA(1,IMU,I) = S(0)
        OUTDATA(2:ME+1,IMU,I) = 0.5*S(1:ME)
        OUTDATA(ME+2:NPHI0(IMU),IMU,I) = 0.0
        CALL COST (NPHI0(IMU),OUTDATA(:,IMU,I),WSAVE(:,IMU))
      ELSE
        ! Else do a slow DFT
        DO IPHI = 1, NPHI0(IMU)
          OUTDATA(IPHI,IMU,I) = SUM( CPHI1(0:ME,IPHI,IMU)*S(0:ME) )
        ENDDO
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE SH_TO_DO



SUBROUTINE DO_TO_SH (NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                     FFTFLAG, CMU2, CPHI2, WSAVE, INDATA, OUTDATA)
 ! Transforms the input field from discrete ordinate space to spherical 
 ! harmonic space for all the points.  (See SH_TO_DO for more).
 ! The SELECTPNTS logical array indicates whether a point should be 
 ! transformed or to use the previous points (for the input layer interfaces).
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0(NMU)
  LOGICAL, INTENT(IN) :: FFTFLAG(NMU)
  REAL,    INTENT(IN) :: INDATA(NPHI0MAX,NMU,NPTS)
  REAL,    INTENT(IN) :: CMU2(NMU,NLM), CPHI2(32,0:16,NMU)
  REAL,    INTENT(IN) :: WSAVE(3*NPHI0MAX+15,NMU)
  REAL,    INTENT(OUT) :: OUTDATA(NLM,NPTS)
  INTEGER :: I, IMU, IPHI, J, L, M, ME, N
  REAL    :: PI, DELPHI
  INTEGER :: MOFJ(NLM)
  REAL    :: S(0:NPHI0MAX-1)

  PI = ACOS(-1.0)
  ! Make the M for each J array
  J = 0
  DO L = 0, ML
    DO M = 0, MIN(L,MM)
      J = J + 1
      MOFJ(J) = M
    ENDDO
  ENDDO

  ! Do the transform for each grid point
  DO I = 1, NPTS
    OUTDATA(:,I) = 0.0
    DO IMU = 1, NMU
      N = NPHI0(IMU)
      ! Figure the max Fourier azimuthal mode we can do for this Nphi0
      ME = MAX(0,MIN(N-2,MOFJ(NLM)))
      ! First do Fourier transform from phi to m for each mu
      IF (FFTFLAG(IMU)) THEN
        ! For cosine modes only use a cosine FFT
        S(0:N-1) = INDATA(1:N,IMU,I)
        CALL COST (N, S(0:N-1), WSAVE(1,IMU))
        DELPHI = PI/MAX(1,N-1)
        S(0:ME) = S(0:ME)*DELPHI
      ELSE
        ! Else do a slow DFT
        S(:) = 0.0
        DO M = 0, ME
          S(M) = SUM( CPHI2(1:N,M,IMU)*INDATA(1:N,IMU,I) )
        ENDDO
      ENDIF
      ! Then do Legendre transform by adding to output for each l, m
      DO J = 1, NLM
        M = MOFJ(J)
        OUTDATA(J,I) = OUTDATA(J,I) + CMU2(IMU,J)*S(M)
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE DO_TO_SH



SUBROUTINE COMPUTE_SOURCE (ML,MM, NLM, NLAY, NPTS, IXP, IXG, &
                           ALBEDOSC, LEGENSC, &
                           SHRADIANCE, CONSTSOURCE, SOURCE, &
                           OLDNPTS, DELSOURCE, &
                           DELJDOT, DELJOLD, DELJNEW, SOLCRIT)
 ! Computes the source function (SOURCE) in spherical harmonic space 
 ! for all the grid levels in all the input layers.  
 ! The thermal source and/or solar pseudo-sources in CONSTSOURCE is
 ! added to the scattering source (computed from LEGENP and the spherical 
 ! harmonic expansion in SHRADIANCE).
 ! To save memory this routine also computes some things for the sequence
 ! acceleration: the dot products of the difference in source function 
 ! between successive iterations and the new difference source function 
 ! vector (DELSOURCE).  
 ! Computes the new solution criterion (SOLCRIT, the RMS difference
 ! in succesive source function fields normalized by the RMS of the field
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ML, MM, NLM, NLAY, NPTS, OLDNPTS
  INTEGER, INTENT(IN) :: IXP(2,NLAY), IXG(NPTS)
  REAL,    INTENT(IN) :: ALBEDOSC(NLAY), LEGENSC(0:ML,NLAY)
  REAL,    INTENT(IN) :: SHRADIANCE(NLM,NPTS), CONSTSOURCE(NLM,NPTS)
  REAL,    INTENT(INOUT) :: SOURCE(NLM,NPTS), DELSOURCE(NLM,OLDNPTS)
  REAL,    INTENT(OUT) :: DELJDOT, DELJOLD, DELJNEW, SOLCRIT
  INTEGER :: LOFJ(NLM), I, J, K, L, LAY, M, ME
  REAL    :: SOURCET, NORM

  ! Make the l index as a function of SH term (J)
  J = 0
  DO L = 0, ML
    ME = MIN(L,MM)
    DO M = 0, ME
      J = J + 1
      LOFJ(J) = L
    ENDDO
  ENDDO

  
  DELJDOT = 0.0
  DELJOLD = 0.0
  DELJNEW = 0.0
  NORM = 0.0
  DO LAY = 1, NLAY
    DO I = IXP(1,LAY), IXP(2,LAY)
      K = IXG(I)
      DO J = 1, NLM
        L = LOFJ(J)
        SOURCET = CONSTSOURCE(J,K) &
                 + SHRADIANCE(J,K)*ALBEDOSC(LAY)*LEGENSC(L,LAY)/(2*L+1)
        IF (K <= OLDNPTS) THEN
          DELJDOT = DELJDOT + (SOURCET-SOURCE(J,K))*DELSOURCE(J,K)
          DELJOLD = DELJOLD + DELSOURCE(J,K)**2
          DELSOURCE(J,K) = SOURCET - SOURCE(J,K)
          DELJNEW = DELJNEW + (SOURCET-SOURCE(J,K))**2
          NORM = NORM + SOURCET**2
        ENDIF
        SOURCE(J,K) = SOURCET
      ENDDO
    ENDDO
  ENDDO
  IF (NORM > 0.0) THEN
    SOLCRIT = SQRT(DELJNEW/NORM)
  ELSE
    SOLCRIT = 1.0
  ENDIF
END SUBROUTINE COMPUTE_SOURCE
 




SUBROUTINE ACCELERATE_SOLUTION (NLM, NPTS, OLDNPTS, SOURCE, DELSOURCE, &
                               DELJDOT, DELJOLD, DELJNEW)
 ! Accelerates the successive order of scattering series.  Uses
 ! information about the ratio in lengths and angle between successive
 ! source function differences to decide how far along the current
 ! source function difference to extrapolate.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLM, NPTS, OLDNPTS
  REAL,    INTENT(IN) :: DELSOURCE(NLM,OLDNPTS)
  REAL,    INTENT(INOUT) :: SOURCE(NLM,NPTS)
  REAL,    INTENT(IN) :: DELJDOT, DELJOLD, DELJNEW
  INTEGER :: I, J, IS, ISD, NS, NSD, NSC
  REAL, SAVE :: A, R, THETA

  ! Accelerate if the grids are the same, didn't last time, and 
  !   things are converging. 
  IF (NPTS == OLDNPTS .AND. A == 0.0 .AND. DELJNEW < DELJOLD) THEN
    ! Compute the acceleration extrapolation factor and apply it.
    R = SQRT(DELJNEW/DELJOLD)
    THETA = ACOS(DELJDOT/SQRT(DELJOLD*DELJNEW))
    A = (1 - R*COS(THETA) + R**(1+0.5*3.14159/THETA)) &
         /(1 + R**2  - 2*R*COS(THETA))  - 1.0
    A = MAX(0.0,A)
    ! WRITE (*,'(1X,A,3(1X,F7.3))') '! Acceleration: ', A,R,THETA
    SOURCE(:,:) = SOURCE(:,:) + A*DELSOURCE(:,:)
  ELSE
    A = 0.0
  ENDIF
END SUBROUTINE ACCELERATE_SOLUTION




SUBROUTINE PATH_INTEGRATION (NLAY, NPTS, IXP, IXG, TAUGSC, &
                             NMU, NPHI0MAX, NPHI0, &
                             SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                             SFCTYPE, SFCPARMS, SFCTEMP, WAVELEN, WAVENO, &
                             MU, PHI, WTDO, SRCINT, LASTNPTS, &
                             FLUXUP, FLUXDN, MEANRAD, DOSOURCE, RADIANCE)
 ! Performs the path integrations through the medium specified by
 ! the optical depth (TAUGSC) and source function (DOSOURCE) in
 ! discrete ordinates.  The integrations are done from the top boundary
 ! (with no reflection) down to the bottom boundary for the downwelling 
 ! angles, and then up from the bottom to the top after computing the 
 ! surface reflection/emission.   The hemispheric flux (FLUXES) and mean
 ! radiance (MEANRAD) are also computed.
 ! The discrete ordinates in MU, PHI must be a complete set with the
 ! downwelling (mu<0) angles first and downwelling and upwelling matching
 ! [MU(J) = -MU(NMU/2+J)].
 ! The surface reflection is handled differently for Lambertian surfaces 
 ! and more general surfaces specified by bidirectional reflection 
 ! distribution functions (BRDF).  For the general BRDF, the bottom 
 ! boundary radiance must be computed for each upwelling ordinate, while 
 ! for the Lambertian surface the boundary radiances can be computed just 
 ! once, since they are isotropic.  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, NPTS, IXP(2,NLAY), IXG(NPTS)
  INTEGER, INTENT(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
  INTEGER, INTENT(INOUT) :: LASTNPTS
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, SKYRAD
  REAL,    INTENT(IN) :: SFCPARMS(*), SFCTEMP, WAVELEN, WAVENO(2)
  REAL,    INTENT(IN) :: MU(NMU), PHI(NPHI0MAX,NMU), WTDO(NMU,NPHI0MAX)
  REAL,    INTENT(INOUT) :: SRCINT(0:4,NMU,NPTS)
  REAL,    INTENT(IN) :: TAUGSC(NPTS), DOSOURCE(NPHI0MAX,NMU,NPTS)
  REAL,    INTENT(OUT) :: FLUXUP(NPTS), FLUXDN(NPTS), MEANRAD(NPTS)
  REAL,    INTENT(OUT) :: RADIANCE(NPHI0MAX,NMU,NPTS)
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE, SFCTYPE
 
  INTEGER :: LAY, I, I1, I2, J, K, IPHI, IMU, M, N
  REAL    :: PI, TEMP, SFCPLANCK, PLANCKFUNC, ALB
  REAL    :: DIRFLUXSFC, TOPRAD, BOTRAD
  DOUBLE PRECISION, ALLOCATABLE :: A(:,:), B(:,:), C(:,:), D(:,:), E(:)
  DOUBLE PRECISION, ALLOCATABLE :: T(:), T1(:), T3(:), T4(:), TR(:)

  PI = ACOS(-1.0)

  M=NMU/2
  
  IF (LASTNPTS /= NPTS) THEN
    ! Calculate and store the source function integration coefficients
    ALLOCATE (A(4,M),B(4,M),C(4,M),D(4,M),E(M), T(M),T1(M),T3(M),T4(M), TR(M))
    DO LAY = 1, NLAY
      I1 = IXP(1,LAY)
      I2 = IXP(2,LAY)-1
      IF (I1 == I2) THEN   ! special case for layer with just two points
        T(:) = (TAUGSC(I1+1) - TAUGSC(I1))/ABS(MU(1:M))
        TR(:) = EXP(-T(:))        ! sublayer transmission along ordinate
        IF (T(1) < 1.0E-4) THEN
          D(1,:) = 0.5*T(:) - (1.0D0/3)*T(:)**2
          D(2,:) = 0.5*T(:) - (1.0D0/6)*T(:)**2
        ELSE
          D(1,:) = (1-TR(:)*(1+T(:)))/T(:)
          D(2,:) = (T(:)-(1-TR(:)))/T(:)
        ENDIF
        SRCINT(0,1:M,I1) = TR(:)
        SRCINT(0,M+1:NMU,I1) = TR(:)
        SRCINT(1,1:M,I1) = D(2,:)
        SRCINT(2,1:M,I1) = D(1,:)
        SRCINT(1,M+1:NMU,I1) = D(1,:)
        SRCINT(2,M+1:NMU,I1) = D(2,:)
      ELSE
        DO I = I1, I2
          A(:,:)=0.0 ; B(:,:)=0.0 ; C(:,:)=0.0 ; D(:,:)=0.0
          T3(:) = (TAUGSC(I+1)-TAUGSC(I))/ABS(MU(1:M))
          IF (I == I1) THEN
            T4(:) = (TAUGSC(I+2)-TAUGSC(I))/ABS(MU(1:M))
            E=T3*T4      ; B(2,:)=1/E ; C(2,:)=-(T3+T4)/E ; D(2,:)=1.0
            E=T3*(T3-T4) ; B(3,:)=1/E ; C(3,:)=-T4/E      ; D(3,:)=0.0
            E=T4*(T4-T3) ; B(4,:)=1/E ; C(4,:)=-T3/E      ; D(4,:)=0.0
          ELSE IF (I < I2) THEN
            T1(:) = (TAUGSC(I-1)-TAUGSC(I))/ABS(MU(1:M))
            T4(:) = (TAUGSC(I+2)-TAUGSC(I))/ABS(MU(1:M))
            E=T1*(T1-T3)*(T1-T4) ; A(1,:)=1/E
            B(1,:)=-(T3+T4)/E ; C(1,:)=T3*T4/E ; D(1,:)=0.0
            E=-T1*T3*T4 ; A(2,:)=1/E 
            B(2,:)=-(T1+T3+T4)/E ; C(2,:)=(T1*T3+T1*T4+T3*T4)/E ; D(2,:)=1.0
            E=T3*(T3-T1)*(T3-T4) ; A(3,:)=1/E
            B(3,:)=-(T1+T4)/E ; C(3,:)=T1*T4/E ; D(3,:)=0.0
            E=T4*(T4-T1)*(T4-T3) ; A(4,:)=1/E
            B(4,:)=-(T1+T3)/E ; C(4,:)=T1*T3/E ; D(4,:)=0.0
          ELSE
            T1(:) = (TAUGSC(I-1)-TAUGSC(I))/ABS(MU(1:M))
            E=T1*(T1-T3) ; B(1,:)=1/E ; C(1,:)=-T3/E      ; D(1,:)=0.0
            E=T1*T3      ; B(2,:)=1/E ; C(2,:)=-(T1+T3)/E ; D(2,:)=1.0
            E=T3*(T3-T1) ; B(3,:)=1/E ; C(3,:)=-T1/E      ; D(3,:)=0.0
          ENDIF
          T(:) = T3(:)
          TR(:) = EXP(-T(:))
          SRCINT(0,1:M,I) = TR(:)
          SRCINT(0,M+1:NMU,I) = TR(:)
          DO J = 1, 4
            SRCINT(J,1:M,I) = A(J,:)*(6*TR+T**3-3*T**2+6*T-6) &
                            + B(J,:)*(-2*TR+T**2-2*T+2) &
                            + C(J,:)*(TR+T-1) + D(J,:)*(1-TR)
            SRCINT(J,M+1:NMU,I) = A(J,:)*(6-TR*(T**3+3*T**2+6*T+6)) &
                                + B(J,:)*(2-TR*(T**2+2*T+2)) &
                                + C(J,:)*(1-TR*(T+1)) + D(J,:)*(1-TR)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    DEALLOCATE (A, B, C, D, E, T, T1, T3, T4, TR)
  ENDIF

  ! Zero the fluxes and mean radiance arrays
  FLUXUP(:) = 0.0  ;  FLUXDN(:) = 0.0  ;  MEANRAD(:) = 0.0


  ! Make the isotropic radiances for the top boundary
  IF (SRCTYPE .EQ. 'T') THEN
    TEMP = SKYRAD
    TOPRAD = PLANCKFUNC(WAVELEN,WAVENO,TEMP)
  ELSE
    TOPRAD = SKYRAD
  ENDIF
  RADIANCE(:,1:NMU/2,IXG(IXP(1,1))) = TOPRAD

  ! First integrate the downwelling zenith angles (these angles must be first)
  ! For downwelling ordinates integrate source function from the top down
  !   Integrations are performed over the vector of azimuths
  DO IMU = 1, NMU/2
    N = NPHI0(IMU)
    DO LAY = 1, NLAY
      I1 = IXP(1,LAY)
      I2 = IXP(2,LAY)-1
      IF (I1 == I2) THEN
        RADIANCE(1:N,IMU,IXG(I1+1)) = SRCINT(0,IMU,I1)*RADIANCE(1:N,IMU,IXG(I1)) &
                             + SRCINT(1,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1)) &
                             + SRCINT(2,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1+1))
      ELSE
        RADIANCE(1:N,IMU,IXG(I1+1)) = SRCINT(0,IMU,I1)*RADIANCE(1:N,IMU,IXG(I1)) &
                             + SRCINT(2,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1)) &
                             + SRCINT(3,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1+1)) &
                             + SRCINT(4,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1+2))
        DO I = I1+1, I2-1
          RADIANCE(1:N,IMU,IXG(I+1)) = SRCINT(0,IMU,I)*RADIANCE(1:N,IMU,IXG(I)) &
                             + SRCINT(1,IMU,I)*DOSOURCE(1:N,IMU,IXG(I-1)) &
                             + SRCINT(2,IMU,I)*DOSOURCE(1:N,IMU,IXG(I)) &
                             + SRCINT(3,IMU,I)*DOSOURCE(1:N,IMU,IXG(I+1)) &
                             + SRCINT(4,IMU,I)*DOSOURCE(1:N,IMU,IXG(I+2))
        ENDDO
        RADIANCE(1:N,IMU,IXG(I2+1)) = SRCINT(0,IMU,I2)*RADIANCE(1:N,IMU,IXG(I2)) &
                             + SRCINT(1,IMU,I2)*DOSOURCE(1:N,IMU,IXG(I2-1)) &
                             + SRCINT(2,IMU,I2)*DOSOURCE(1:N,IMU,IXG(I2)) &
                             + SRCINT(3,IMU,I2)*DOSOURCE(1:N,IMU,IXG(I2+1))
      ENDIF
      IF (LAY<NLAY) RADIANCE(:,IMU,IXG(IXP(1,LAY+1))) = RADIANCE(:,IMU,IXG(IXP(2,LAY)))
    ENDDO

    ! Integrate over the DO radiances to get the downwelling fluxes
    DO I = 1, NPTS
      DO K = 1, NPHI0(IMU)
        FLUXDN(I) = FLUXDN(I) + ABS(MU(IMU))*WTDO(IMU,K)*RADIANCE(K,IMU,IXG(I))
        MEANRAD(I) =  MEANRAD(I) + WTDO(IMU,K)*RADIANCE(K,IMU,IXG(I))/(4*PI)
      ENDDO
    ENDDO
  ENDDO

 
    ! Now do all the surface reflection stuff
  IF (SRCTYPE /= 'T') THEN
    DIRFLUXSFC = SOLARFLUX*EXP(-TAUGSC(IXP(2,NLAY))/ABS(SOLARMU))
  ENDIF
  SFCPLANCK = 0.0
  IF (SRCTYPE /= 'S') THEN
    SFCPLANCK = PLANCKFUNC(WAVELEN,WAVENO,SFCTEMP)
  ENDIF

  ! For a Lambertian surface calculate the isotropic surface radiance
  !   from the thermal emission or direct reflection and the reflected 
  !   downwelling flux at the surface.
  IF (SFCTYPE == 'L') THEN
    ALB = SFCPARMS(1)
    IF (SRCTYPE == 'T') THEN
      BOTRAD = (1.0-ALB)*SFCPLANCK
    ELSE IF (SRCTYPE == 'S') THEN
      BOTRAD = ALB*DIRFLUXSFC/PI
    ELSE IF (SRCTYPE .EQ. 'B') THEN
      BOTRAD = (1.0-ALB)*SFCPLANCK + ALB*DIRFLUXSFC/PI
    ENDIF
    BOTRAD = BOTRAD + FLUXDN(IXP(2,NLAY))*ALB/PI
    RADIANCE(:,NMU/2+1:NMU,IXG(IXP(2,NLAY))) = BOTRAD
  ELSE
    ! If not a Lambertian surface, compute the radiance for each of the
    !   upwelling ordinates by integrating the BRDF over the stored 
    !   downwelling radiances.
    DO IMU = NMU/2+1, NMU
      DO IPHI = 1, NPHI0(IMU)
        CALL VARIABLE_BRDF_SURFACE (NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                                    MU(IMU), PHI(IPHI,IMU), &
                                    SRCTYPE, SOLARMU, DIRFLUXSFC, &
                                    SFCTYPE, SFCPARMS, SFCPLANCK, &
                                    RADIANCE(:,:,IXG(IXP(2,NLAY))), BOTRAD)
        RADIANCE(IPHI,IMU,IXG(IXP(2,NLAY))) = BOTRAD
      ENDDO
    ENDDO
  ENDIF


  ! Then integrate the upwelling zenith angles
      ! For upwelling ordinates integrate source function from the bottom up
  DO IMU = NMU/2+1, NMU
    N = NPHI0(IMU)
    DO LAY = NLAY, 1, -1
      I1 = IXP(1,LAY)
      I2 = IXP(2,LAY)-1
      IF (I1 == I2) THEN
        RADIANCE(1:N,IMU,IXG(I1)) = SRCINT(0,IMU,I1)*RADIANCE(1:N,IMU,IXG(I1+1)) &
                             + SRCINT(1,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1)) &
                             + SRCINT(2,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1+1))
      ELSE
        RADIANCE(1:N,IMU,IXG(I2)) = SRCINT(0,IMU,I2)*RADIANCE(1:N,IMU,IXG(I2+1)) &
                             + SRCINT(1,IMU,I2)*DOSOURCE(1:N,IMU,IXG(I2-1)) &
                             + SRCINT(2,IMU,I2)*DOSOURCE(1:N,IMU,IXG(I2)) &
                             + SRCINT(3,IMU,I2)*DOSOURCE(1:N,IMU,IXG(I2+1))
        DO I = I2-1, I1+1, -1
          RADIANCE(1:N,IMU,IXG(I)) = SRCINT(0,IMU,I)*RADIANCE(1:N,IMU,IXG(I+1)) &
                             + SRCINT(1,IMU,I)*DOSOURCE(1:N,IMU,IXG(I-1)) &
                             + SRCINT(2,IMU,I)*DOSOURCE(1:N,IMU,IXG(I)) &
                             + SRCINT(3,IMU,I)*DOSOURCE(1:N,IMU,IXG(I+1)) &
                             + SRCINT(4,IMU,I)*DOSOURCE(1:N,IMU,IXG(I+2))
        ENDDO
        RADIANCE(1:N,IMU,IXG(I1)) = SRCINT(0,IMU,I1)*RADIANCE(1:N,IMU,IXG(I1+1)) &
                             + SRCINT(2,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1)) &
                             + SRCINT(3,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1+1)) &
                             + SRCINT(4,IMU,I1)*DOSOURCE(1:N,IMU,IXG(I1+2))
      ENDIF
      IF (LAY>1) RADIANCE(:,IMU,IXG(IXP(2,LAY-1))) = RADIANCE(:,IMU,IXG(IXP(1,LAY)))
    ENDDO

    ! Integrate over the DO radiances to get the upwelling fluxes
    DO I = 1, NPTS
      DO K = 1, NPHI0(IMU)
        FLUXUP(I) = FLUXUP(I) + ABS(MU(IMU))*WTDO(IMU,K)*RADIANCE(K,IMU,IXG(I))
        MEANRAD(I) =  MEANRAD(I) + WTDO(IMU,K)*RADIANCE(K,IMU,IXG(I))/(4*PI)
      ENDDO
    ENDDO
  ENDDO

  LASTNPTS=NPTS
END SUBROUTINE PATH_INTEGRATION
 



SUBROUTINE SPLIT_CELLS (NLAY, NPTS, MAXIG, IXP, IXG, TAUGSC, &
                        NMU, NPHI0MAX, NPHI0, MU, DOSOURCE, RADIANCE, &
                        SPLITACC, TRUESOLCRIT, CURSPLITACC)
 ! Determines the current cell splitting accuracy (CURSPLITACC) from the
 ! solution accuracy (SOLCRIT) and the desired final cell splitting accuracy
 ! (SPLITACC).  Calculates the cell splitting criterion for each sublayer, 
 ! The cell splitting criterion for a sublayer is the rms difference over 
 ! outgoing ordinates in radiance between the original cubic source function
 ! interpolation integration (in RADIANCE) and a linear interpolation 
 ! integration, normalized by the rms of RADIANCE over all ordinates and 
 ! levels.  Splits cells with criterion greater than the current splitting
 ! accuracy.  Interpolates the discrete ordinate source function to the new
 ! levels.  The number of levels (NPTS), the adaptive index arrays (IXP, IXG),
 ! and the optical depth grid array (TAUGSC) are also updated.
 ! The boundary condition information is not needed for this routine because
 ! it has already been computed by PATH_INTEGRATION and stored in RADIANCE.
 ! SPLITCRITRATIO is the factor needed to divide the cell splitting criterion
 ! to convert it to a typical radiance accuracy of the split sublayer with
 ! the cubic interpolation integration.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, MAXIG
  INTEGER, INTENT(INOUT) :: NPTS, IXP(2,NLAY), IXG(*)
  INTEGER, INTENT(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
  REAL,    INTENT(IN) :: MU(NMU)
  REAL,    INTENT(IN) :: SPLITACC, TRUESOLCRIT
  REAL,    INTENT(IN) :: RADIANCE(NPHI0MAX,NMU,NPTS)
  REAL,    INTENT(INOUT) :: CURSPLITACC
  REAL,    INTENT(INOUT) :: TAUGSC(*), DOSOURCE(NPHI0MAX,NMU,*)
  INTEGER :: LAY, I, IMU, IMU2, N
  INTEGER :: IXPN(2,NLAY), NEWNPTS, I1, I2, I3, I4
  REAL    :: SPLITCRIT, NORM, T1, T2, T3, T4
  DOUBLE PRECISION  :: TAU(NMU/2), TR(NMU/2)
  DOUBLE PRECISION, ALLOCATABLE :: SRCINT1(:,:,:), DIFFRAD(:,:)

  NORM = SQRT( SUM(RADIANCE(:,:,1:NPTS)**2)/(NPTS*SUM(NPHI0(:))) )
   ! Leave immediately if splitting is not desired
  IF (SPLITACC <= 0.0 .OR. NORM == 0.0) RETURN

  ALLOCATE (SRCINT1(0:2,NMU/2,NPTS), DIFFRAD(NPHI0MAX,NMU))
  DIFFRAD(:,:) = 0.0

  CURSPLITACC = MAX(SPLITACC,MIN(CURSPLITACC,TRUESOLCRIT))

  ! Calculate and store the linear source function integration coefficients
  DO LAY = 1, NLAY
    DO I = IXP(1,LAY), IXP(2,LAY)-1
      TAU(:) = (TAUGSC(I+1) - TAUGSC(I))/ABS(MU(1:NMU/2))
      TR(:) = EXP(-TAU(:))        ! sublayer transmission along ordinate
      SRCINT1(0,:,I) = TR(:)
      IF (TAU(1) < 1.0E-4) THEN
        SRCINT1(1,:,I) = 0.5*TAU(:) - (1.0D0/3)*TAU(:)**2
        SRCINT1(2,:,I) = 0.5*TAU(:) - (1.0D0/6)*TAU(:)**2
      ELSE
        SRCINT1(1,:,I) = (1-TR(:)*(1+TAU(:)))/TAU(:)
        SRCINT1(2,:,I) = (TAU(:)-(1-TR(:)))/TAU(:)
      ENDIF
    ENDDO
  ENDDO

  ! For each sublayer integrate the downwelling zenith angles from top to
  ! bottom and the upwelling angles from bottom to top using linear 
  ! source function interpolation and difference from the previous
  ! cubic/quadratic interpolation results.
  NEWNPTS = NPTS
  IXPN(:,:) = IXP(:,:)
  DO LAY = 1, NLAY
   IF (IXP(1,LAY) < IXP(2,LAY)-1) THEN
    DO I = IXP(1,LAY), IXP(2,LAY)-1
      DO IMU = 1, NMU/2
        IMU2 = IMU + NMU/2
        N = NPHI0(IMU)
        DIFFRAD(1:N,IMU) = -RADIANCE(1:N,IMU,IXG(I+1)) &
                           + SRCINT1(0,IMU,I)*RADIANCE(1:N,IMU,IXG(I)) &
                           + SRCINT1(1,IMU,I)*DOSOURCE(1:N,IMU,IXG(I)) &
                           + SRCINT1(2,IMU,I)*DOSOURCE(1:N,IMU,IXG(I+1))
        DIFFRAD(1:N,IMU2) = -RADIANCE(1:N,IMU2,IXG(I)) &
                           + SRCINT1(0,IMU,I)*RADIANCE(1:N,IMU2,IXG(I+1)) &
                           + SRCINT1(1,IMU,I)*DOSOURCE(1:N,IMU2,IXG(I+1)) &
                           + SRCINT1(2,IMU,I)*DOSOURCE(1:N,IMU2,IXG(I))
      ENDDO
      SPLITCRIT = SQRT( SUM(DIFFRAD(:,:)**2)/SUM(NPHI0(:)) ) /NORM
        ! If the splitting criterion is large enough then add a level in
        !   the center of the current sublayer.  Put the new level at the
        !   end of the radiance array and update the pointers.
      IF (SPLITCRIT>CURSPLITACC) THEN
        IF (NEWNPTS>=MAXIG) THEN
          WRITE (*,*) 
          WRITE (*,*) 'Adaptive grid exceeded MAXIG=',MAXIG
          WRITE (*,*) 'Increase SPLITTING_FACTOR'
          STOP
        ENDIF
        NEWNPTS = NEWNPTS + 1
        TAUGSC(NEWNPTS) = 0.5*(TAUGSC(I)+TAUGSC(I+1))
        ! print *, newnpts,taug(newnpts), splitcrit
        IXG(NEWNPTS) = NEWNPTS
        IXPN(2,LAY:NLAY) = IXPN(2,LAY:NLAY) + 1
        IXPN(1,LAY+1:NLAY) = IXPN(1,LAY+1:NLAY) + 1
        I1 = MAX(I-1,IXP(1,LAY))  ;   I2 = I
        I3 = I+1  ; I4 = MIN(I+2,IXP(2,LAY))
        CALL INTERP_SOURCE (NPHI0MAX, NMU, TAUGSC(I1), DOSOURCE(:,:,IXG(I1)), &
                            TAUGSC(I2), DOSOURCE(:,:,IXG(I2)), &
                            TAUGSC(I3), DOSOURCE(:,:,IXG(I3)), &
                            TAUGSC(I4), DOSOURCE(:,:,IXG(I4)), &
                            TAUGSC(NEWNPTS), DOSOURCE(:,:,IXG(NEWNPTS)) )
      ENDIF
    ENDDO
   ENDIF
  ENDDO

  IF (NEWNPTS > NPTS) THEN
    CALL SSORT (TAUGSC, IXG, NEWNPTS, 2)
    NPTS = NEWNPTS
    IXP(:,:) = IXPN(:,:)
     ! Make sure that the IXG pointers are in the correct order
    DO I = 1, NPTS-1
      IF (TAUGSC(I) == TAUGSC(I+1) .AND. IXG(I) > IXG(I+1)) THEN
        I1 = IXG(I)
        IXG(I) = IXG(I+1)
        IXG(I+1) = I1
      ENDIF
    ENDDO
  ENDIF
  DEALLOCATE (SRCINT1, DIFFRAD)
END SUBROUTINE SPLIT_CELLS
 


SUBROUTINE INTERP_SOURCE (NPH, NMU, T1,S1, T2,S2, T3,S3, T4,S4, T,S)
 ! Does either quadratic or cubic polynomial interpolation of the discrete
 ! ordinate source function as a function of optical depth.  The source 
 ! function (S1,S2,S3,S4) is specified at three or four optical depths 
 ! (T1,T2,T3,T4), interpolated to T, and output in S.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NPH, NMU
  REAL,    INTENT(IN) :: T1, T2, T3, T4, T
  REAL,    INTENT(IN) :: S1(NPH,NMU), S2(NPH,NMU), S3(NPH,NMU), S4(NPH,NMU)
  REAL,    INTENT(OUT) :: S(NPH,NMU)

  IF (ABS(T1-T2) < 4*SPACING(T1)) THEN
    S = (T-T3)*(T-T4)*S2/((T2-T3)*(T2-T4)) &
      + (T-T2)*(T-T4)*S3/((T3-T2)*(T3-T4)) &
      + (T-T2)*(T-T3)*S4/((T4-T2)*(T4-T3))
  ELSE IF (ABS(T3-T4) < 4*SPACING(T3)) THEN
    S = (T-T2)*(T-T3)*S1/((T1-T2)*(T1-T3)) &
      + (T-T1)*(T-T3)*S2/((T2-T1)*(T2-T3)) &
      + (T-T1)*(T-T2)*S3/((T3-T1)*(T3-T2))
  ELSE
    S = (T-T2)*(T-T3)*(T-T4)*S1/((T1-T2)*(T1-T3)*(T1-T4)) &
      + (T-T1)*(T-T3)*(T-T4)*S2/((T2-T1)*(T2-T3)*(T2-T4)) &
      + (T-T1)*(T-T2)*(T-T4)*S3/((T3-T1)*(T3-T2)*(T3-T4)) &
      + (T-T1)*(T-T2)*(T-T3)*S4/((T4-T1)*(T4-T2)*(T4-T3))
  ENDIF
END SUBROUTINE INTERP_SOURCE



SUBROUTINE VARIABLE_BRDF_SURFACE (NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                                  MU2, PHI2, SRCTYPE, SOLARMU, DIRFLUXSFC, &
                                  SFCTYPE, SFCPARMS, SFCPLANCK, &
                                  SFCRADIANCE, BOTRAD)
 ! Computes the upwelling radiance at the bottom boundary for one 
 ! outgoing direction using the specified bidirectional reflectance 
 ! distribution function.  The upwelling radiance includes the reflection of
 ! the incident radiance, the thermal emission (emissivity*Planck function),
 ! and the reflected direct solar flux (if applicable).  The upwelling 
 ! radiance is the integral over all incident directions of the BRDF times 
 ! the downwelling radiance, so a discrete sum is done and the integral 
 ! weights (WTDO) are included.  The general BRDF function is called 
 ! to compute the reflectance for the particular BRDF type (SFCTYPE) 
 ! with parameters (SFCPARMS) for the incident and outgoing directions.  
 ! The emissivity is computed implicitly from the integral of the BRDF.  
 ! The outgoing direction is specified with (MU2,PHI2), and the BRDF is 
 ! computed for all incident directions (loop over JMU,JPHI).  The incident 
 ! downwelling radiances are input in SFCRADIANCE(:,:) and the outgoing 
 ! upwelling radiance is output in BOTRAD.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
  REAL,    INTENT(IN) :: MU(NMU), PHI(NPHI0MAX,NMU), WTDO(NMU,NPHI0MAX)
  REAL,    INTENT(IN) :: MU2, PHI2, SOLARMU, DIRFLUXSFC, SFCPLANCK
  REAL,    INTENT(IN) :: SFCPARMS(*), SFCRADIANCE(NPHI0MAX,NMU)
  REAL,    INTENT(OUT):: BOTRAD
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE, SFCTYPE
  INTEGER :: JMU, JPHI
  REAL    :: OPI, SOLARAZ, REFLECT, W

  OPI = 1.0/ACOS(-1.0)
  SOLARAZ = 0.0

  ! Initialize the upwelling boundary radiances to zero (for thermal)
  !   or to the reflected direct solar flux.
  IF (SRCTYPE == 'T') THEN
    BOTRAD = 0.0
  ELSE
    CALL SURFACE_BRDF (SFCTYPE, SFCPARMS, MU2, PHI2, SOLARMU, SOLARAZ, &
                       REFLECT)
    BOTRAD = OPI*REFLECT*DIRFLUXSFC
  ENDIF

  ! Integrate over the incident discrete ordinate directions (JMU,JPHI)
  DO JMU = 1, NMU/2
    DO JPHI = 1, NPHI0(JMU)
      CALL SURFACE_BRDF (SFCTYPE, SFCPARMS, MU2, PHI2, &
                         MU(JMU),PHI(JPHI,JMU), REFLECT)
      W = OPI*ABS(MU(JMU))*WTDO(JMU,JPHI)
      BOTRAD = BOTRAD + W*REFLECT*SFCRADIANCE(JPHI,JMU) &
                      + W*(1-REFLECT)*SFCPLANCK
    ENDDO
  ENDDO

END SUBROUTINE VARIABLE_BRDF_SURFACE




SUBROUTINE SURFACE_BRDF (SFCTYPE, REFPARMS, &
                         MU2, PHI2, MU1, PHI1, REFLECT)
 !  Returns the reflection coefficient for the general bidirectional
 ! reflection distribution function of the specified type (SFCTYPE).
 ! The incident direction is (MU1,PHI1), and the outgoing direction
 ! is (MU2,PHI2) (MU is cosine of zenith angle, and PHI is the azimuthal
 ! angle in radians).  The incident directions have mu<0 (downward), 
 ! while the outgoing directions have mu>0 (upward). The reflection 
 ! function is normalized so that for a Lambertian surface (uniform 
 ! reflection) the returned value (REFLECT) is simply the albedo.
 ! This routine calls the desired BRDF function, passing the 
 ! appropriate parameters.  More BRDF surface types may be added easily 
 ! by putting in the appropriate function calls.
 !        Type        Parameters
 !    R  Rahman et al  rho0, k, Theta
  IMPLICIT NONE
  REAL, INTENT(IN) :: REFPARMS(*), MU1, PHI1, MU2, PHI2
  REAL, INTENT(OUT):: REFLECT
  CHARACTER(LEN=1), INTENT(IN) :: SFCTYPE
  REAL  :: PI
  REAL  :: RPV_REFLECTION

  PI = ACOS(-1.0)
  IF (SFCTYPE == 'R') THEN
    ! R: Rahman, Pinty, and Verstraete
    REFLECT = RPV_REFLECTION (REFPARMS(1),REFPARMS(2),REFPARMS(3), &
                              -MU1, MU2, PHI1-PHI2-PI)
  ELSE
    STOP 'SURFACE_BRDF: Unknown BRDF type'
  ENDIF
END SUBROUTINE SURFACE_BRDF




FUNCTION RPV_REFLECTION (RHO0, K, THETA, MU1, MU2, PHI)
 !  Computes the Rahman, Pinty, Verstraete BRDF.  The incident
 ! and outgoing cosine zenith angles are MU1 and MU2, respectively,
 ! and the relative azimuthal angle is PHI.  In this case the incident
 ! direction is where the radiation is coming from (i.e. opposite of the
 ! discrete ordinate), so MU1>0 and the hot spot is MU2=MU1 and PHI=0.
 ! The reference is:
 !    Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere 
 !    Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable 
 !    With NOAA Advanced Very High Resolution Radiometer Data,
 !    J. Geophys. Res., 98, 20791-20801.
  IMPLICIT NONE
  REAL, INTENT(IN) :: RHO0, K, THETA, MU1, MU2, PHI
  REAL :: RPV_REFLECTION
  REAL M, F, H, COSPHI, SIN1, SIN2, COSG, TAN1, TAN2, CAPG

  M = MU1**(K-1) * MU2**(K-1) / (MU1 + MU2)**(1-K)
  COSPHI = COS(PHI)
  SIN1 = SQRT(1.0-MU1**2)
  SIN2 = SQRT(1.0-MU2**2)
  COSG = MU1*MU2 + SIN1*SIN2*COSPHI
  F = (1-THETA**2) / (1 + 2*THETA*COSG + THETA**2)**1.5
  TAN1 = SIN1/MU1
  TAN2 = SIN2/MU2
  CAPG = SQRT( TAN1**2 + TAN2**2 - 2*TAN1*TAN2*COSPHI )
  H = 1 + (1-RHO0)/(1+CAPG)
  RPV_REFLECTION = RHO0 * M * F * H
END FUNCTION RPV_REFLECTION




SUBROUTINE COMPUTE_RADIANCE (NLAY, ALBEDOSC, LEGENSC, &
                             MAXLEG, NLEGP, LEGENP,  DELTAM, &
                             NPTS, HEIGHTG, TAUG, TAUGSC, IXP, IXG, &
                             ML, MM, NLM, NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                             SRCTYPE, SOLARFLUX, SOLARMU, &
                             SKYRAD, WAVELEN, WAVENO, SFCTEMP, &
                             SFCTYPE, SFCPARMS, SOURCE, RADIANCE, &
                             NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                             RADOUT)
 ! Computes the radiances for the specified tau's, mu's, and phi's.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, MAXLEG, NLEGP(NLAY)
  INTEGER, INTENT(IN) :: NPTS, IXP(2,NLAY), IXG(NPTS)
  INTEGER, INTENT(IN) :: ML, MM, NLM, NMU, NPHI0MAX, NPHI0(NMU)
  INTEGER, INTENT(IN) :: NZOUT, NMUOUT, NPHIOUT
  LOGICAL, INTENT(IN) :: DELTAM
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE, SFCTYPE
  REAL,    INTENT(IN) :: ALBEDOSC(NLAY), LEGENSC(0:Ml,NLAY)
  REAL,    INTENT(IN) :: LEGENP(0:MAXLEG,NLAY)
  REAL,    INTENT(IN) :: HEIGHTG(NPTS), TAUG(NPTS), TAUGSC(NPTS)
  REAL,    INTENT(IN) :: MU(NMU), PHI(NPHI0MAX,NMU), WTDO(NMU,NPHI0MAX)
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, SKYRAD
  REAL,    INTENT(IN) :: SFCPARMS(*), SFCTEMP, WAVELEN, WAVENO(2)
  REAL,    INTENT(IN) :: SOURCE(NLM,NPTS), RADIANCE(NPHI0MAX,NMU,NPTS)
  REAL,    INTENT(INOUT) :: ZOUT(NZOUT), MUOUT(NMUOUT), PHIOUT(NPHIOUT)
  REAL,    INTENT(OUT) :: RADOUT(NZOUT,NMUOUT,NPHIOUT)

  INTEGER :: I, J, K, L
  REAL    :: ZDIF, RADDEG, PHIRAD
  REAL    :: DIRFLUXSFC, SFCPLANCK, BNDRAD, PLANCKFUNC
  INTEGER, ALLOCATABLE :: LZOUT(:)
  REAL, ALLOCATABLE :: DIRFLUX(:), SOURCE1(:), RADVEC(:)

  RADDEG = ACOS(-1.0)/180.0
  ALLOCATE (DIRFLUX(NPTS), SOURCE1(NPTS), RADVEC(NPTS), LZOUT(NZOUT))

  ! Find the closest level in height to the output heights
  DO I = 1, NZOUT
    ZDIF = 1.0E10
    DO L = 1, NPTS
      IF (ABS(HEIGHTG(L)-ZOUT(I)) < ZDIF) THEN
        LZOUT(I) = L
        ZDIF = ABS(HEIGHTG(L)-ZOUT(I))
      ENDIF
    ENDDO
    ZOUT(I) = HEIGHTG(LZOUT(I))
  ENDDO

  DO L = 1, NPTS
    ! Make the perpendicular direct beam flux for all the levels
    DIRFLUX(IXG(L)) = SOLARFLUX/ABS(SOLARMU) *EXP(-TAUGSC(L)/ABS(SOLARMU))
  ENDDO

  IF (SFCTYPE /= 'L') THEN
    IF (SRCTYPE /= 'T') THEN
      DIRFLUXSFC = SOLARFLUX*EXP(-TAUGSC(IXP(2,NLAY))/ABS(SOLARMU))
      SFCPLANCK = 0.0
    ELSE
      SFCPLANCK = PLANCKFUNC(WAVELEN,WAVENO,SFCTEMP)
    ENDIF
  ENDIF


   ! Loop over the output radiance directions
  DO J = 1, NMUOUT
    DO K = 1, NPHIOUT
      PHIRAD = PHIOUT(K)*RADDEG
      IF (MUOUT(J) == 0.0 .OR. ABS(MUOUT(J)) > 1.0) THEN
        WRITE (*,*) 'COMPUTE_RADIANCE: Bad mu for radiance',MUOUT(J)
        RADOUT(:,J,K) = -999.0
        CYCLE
      ENDIF      

      ! Compute the source function throughout grid for this angle
      CALL COMPUTE_ONE_SOURCE (NLAY, ALBEDOSC, LEGENSC, &
                               MAXLEG, NLEGP, LEGENP, &
                               NPTS, IXP, IXG, ML, MM, NLM, &
                               MUOUT(J), PHIRAD, SRCTYPE, DELTAM, &
                               SOLARMU, DIRFLUX, SOURCE, SOURCE1)

      ! Get boundary radiances: either top or bottom
      !  Isotropic top boundary or Lambertian bottom boundary can use
      !  the previously computed boundary radiances in RADIANCE,
      !  otherwise, compute the radiance for this angle by integrating 
      !  over the stored downwelling radiances.
      IF (MUOUT(J) < 0.0) THEN
        BNDRAD = RADIANCE(1,1,IXG(IXP(1,1)))
      ELSE
        IF (SFCTYPE == 'L') THEN
          BNDRAD = RADIANCE(1,NMU/2+1,IXG(IXP(2,NLAY)))
        ELSE
          CALL VARIABLE_BRDF_SURFACE (NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                                      MUOUT(J), PHIRAD, &
                                      SRCTYPE, SOLARMU, DIRFLUXSFC, &
                                      SFCTYPE, SFCPARMS, SFCPLANCK, &
                                      RADIANCE(:,:,IXG(IXP(2,NLAY))), BNDRAD)
        ENDIF
      ENDIF

      ! Integrate the source function for this direction to get the radiances
      !   at all the levels (RADVEC), then pull out the desired one.
      CALL INTEGRATE_SOURCE (NLAY, NPTS, IXP, IXG, TAUGSC, &
                             MUOUT(J), BNDRAD, SOURCE1, RADVEC)
      DO I = 1, NZOUT
        RADOUT(I,J,K) = RADVEC(LZOUT(I))
      ENDDO
    ENDDO
  ENDDO

  DEALLOCATE (RADVEC, SOURCE1, DIRFLUX, LZOUT)
END SUBROUTINE COMPUTE_RADIANCE



SUBROUTINE COMPUTE_ONE_SOURCE (NLAY, ALBEDOSC, LEGENSC, &
                               MAXLEG, NLEGP, LEGENP, &
                               NPTS, IXP, IXG, ML, MM, NLM, &
                               MUO, PHIO, SRCTYPE, DELTAM, &
                               SOLARMU, DIRFLUX, SOURCE, SOURCE1)
 ! Computes the source function (SOURCE1) in the direction (MUO,PHIO)
 ! for all levels (NPTS).  The spherical harmonic source function
 ! series is input in SOURCE.  DIRFLUX is the attenuated direct beam 
 ! perpendicular flux.
 ! For a solar source if delta-M then use Nakajima and Tanaka TMS
 ! procedure, replacing delta-M single scattering with single scattering
 ! for unscaled untruncated phase function.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, MAXLEG, NLEGP(NLAY), NPTS
  INTEGER, INTENT(IN) :: ML, MM, NLM, IXP(2,NLAY), IXG(NPTS)
  LOGICAL, INTENT(IN) :: DELTAM
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE
  REAL,    INTENT(IN) :: ALBEDOSC(NLAY), LEGENSC(0:ML,NLAY)
  REAL,    INTENT(IN) :: LEGENP(0:MAXLEG,NLAY)
  REAL,    INTENT(IN) :: MUO, PHIO, SOLARMU, DIRFLUX(NPTS)
  REAL,    INTENT(IN) :: SOURCE(NLM,NPTS)
  REAL,    INTENT(OUT) :: SOURCE1(NPTS)
  INTEGER :: I, J, K, L, LAY, M, ME
  DOUBLE PRECISION  :: F, A, COSSCAT
  REAL, ALLOCATABLE :: YLMDIR(:), YLMSUN(:)
  DOUBLE PRECISION, ALLOCATABLE :: SUNDIRLEG(:), SINGSCAT(:)

  ALLOCATE (YLMDIR(NLM), YLMSUN(NLM), SUNDIRLEG(0:MAXLEG), SINGSCAT(NLAY))

  ! Precompute Ylm's for output direction
  CALL YLMALL (MUO, PHIO, ML, MM, 1, YLMDIR)

  ! Special case for solar source and Delta-M
  IF (SRCTYPE /= 'T' .AND. DELTAM) THEN
    ! Compute the Legendre polynomials for the scattering angle 
    !   for the untruncated solar single scattering computation.
    COSSCAT = SOLARMU*MUO + SQRT((1.0D0-SOLARMU**2)*(1.0D0-MUO**2))*COS(PHIO)
    CALL LEGENDRE_ALL (COSSCAT, MAXLEG, SUNDIRLEG)
    SUNDIRLEG(:) = SUNDIRLEG(:)/(4.0D0*ACOS(-1.0D0))
    CALL YLMALL (SOLARMU, 0.0, ML, MM, 1, YLMSUN)

    ! Precompute the single scattering contribution for the original 
    ! unscaled phase function (except for the direct beam transmission).
    ! Need to put the inverse of the tau scaling in the source function 
    ! because tau is still scaled.
    DO LAY = 1, NLAY
      SINGSCAT(LAY) = 0.0
      IF (ML+1 > NLEGP(LAY)) THEN
        F = 0
      ELSE
        F = LEGENP(ML+1,LAY)/(2*(ML+1)+1)
      ENDIF
      DO L = 0, NLEGP(LAY)
        SINGSCAT(LAY) = SINGSCAT(LAY) & 
                       + SUNDIRLEG(L) *ALBEDOSC(LAY)*LEGENP(L,LAY)/(1-F)
      ENDDO
      IF (SINGSCAT(LAY) < 0.0) THEN
        WRITE (*,*) 'COMPUTE_ONE_SOURCE: negative source function for layer: ',LAY
        STOP
      ENDIF
    ENDDO

    ! Loop over all the grid points, computing the source function 
    !   at the viewing angle from the spherical harmonic source function.
    DO LAY = 1, NLAY
      DO I = IXP(1,LAY), IXP(2,LAY)
        K = IXG(I)
        ! Sum over the spherical harmonic series of the source function
        SOURCE1(K) = SUM(SOURCE(:,K)*YLMDIR(:))
        ! First subtract off the truncated single scattering 
        J = 1
        DO L = 0, ML
          ME = MIN(L,MM)
          A = DIRFLUX(K)*ALBEDOSC(LAY)*LEGENSC(L,LAY)/(2*L+1)
          DO M = 0, ME
            SOURCE1(K) = SOURCE1(K) - A*YLMDIR(J)*YLMSUN(J)
            J = J + 1
          ENDDO
        ENDDO
        ! Then add in the single scattering contribution
        SOURCE1(K) = SOURCE1(K) + DIRFLUX(K)*SINGSCAT(LAY)
      ENDDO
    ENDDO

  ELSE  ! for thermal sources or no delta-M it is very simple:
    ! Sum over the spherical harmonic series of the source function
    DO K = 1, NPTS
      SOURCE1(K) = SUM(SOURCE(:,K)*YLMDIR(:))
    ENDDO
  ENDIF

  DEALLOCATE (YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT)
END SUBROUTINE COMPUTE_ONE_SOURCE
 



SUBROUTINE INTEGRATE_SOURCE (NLAY, NPTS, IXP, IXG, TAUGSC, &
                             MUO, BNDRAD, SOURCE1, RADVEC)
 ! Integrates the radiative transfer equation through the medium specified 
 ! by the optical depth (TAUGSC) and source function (SOURCE1) in direction
 ! cosine MUO.  The radiance at the approprate boundary (top or bottom
 ! depending on the sign of MUO) is input (BNDRAD).  The radiance for
 ! the direction the source function was made at all the adaptive grid 
 ! levels is output (RADVEC with indexing corresponding to TAUGSC).
 ! Uses the same quadratic/cubic polynomial source function interpolation
 ! as in PATH_INTEGRATION.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, NPTS, IXP(2,NLAY), IXG(NPTS)
  REAL,    INTENT(IN) :: MUO, BNDRAD
  REAL,    INTENT(IN) :: TAUGSC(NPTS), SOURCE1(NPTS)
  REAL,    INTENT(OUT) :: RADVEC(NPTS)
  INTEGER :: LAY, I, I1, I2
  REAL, ALLOCATABLE :: SRCINT(:,:)
  DOUBLE PRECISION :: A(4), B(4), C(4), D(4), E, T, T1, T3, T4, TR

  ! Calculate and store the source function integration coefficients
  ALLOCATE (SRCINT(0:4,NPTS))
  DO LAY = 1, NLAY
    I1 = IXP(1,LAY)
    I2 = IXP(2,LAY)-1
    IF (I1 == I2) THEN   ! special case for layer with just two points
      T = (TAUGSC(I1+1) - TAUGSC(I1))/ABS(MUO)
      TR = EXP(-T)        ! sublayer transmission along ordinate
      IF (T < 1.0E-4) THEN
        D(1) = 0.5*T - (1.0D0/3)*T**2
        D(2) = 0.5*T - (1.0D0/6)*T**2
      ELSE
        D(1) = (1-TR*(1+T))/T
        D(2) = (T-(1-TR))/T
      ENDIF
      SRCINT(0,I1) = TR
      IF (MUO < 0.0) THEN
        SRCINT(1,I1) = D(2)
        SRCINT(2,I1) = D(1)
      ELSE
        SRCINT(1,I1) = D(1)
        SRCINT(2,I1) = D(2)
      ENDIF
    ELSE
      DO I = I1, I2
        A(:)=0.0 ; B(:)=0.0 ; C(:)=0.0 ; D(:)=0.0
        T3 = (TAUGSC(I+1)-TAUGSC(I))/ABS(MUO)
        IF (I == I1) THEN
          T4 = (TAUGSC(I+2)-TAUGSC(I))/ABS(MUO)
          E=T3*T4      ; B(2)=1/E ; C(2)=-(T3+T4)/E ; D(2)=1.0
          E=T3*(T3-T4) ; B(3)=1/E ; C(3)=-T4/E      ; D(3)=0.0
          E=T4*(T4-T3) ; B(4)=1/E ; C(4)=-T3/E      ; D(4)=0.0
        ELSE IF (I < I2) THEN
          T1 = (TAUGSC(I-1)-TAUGSC(I))/ABS(MUO)
          T4 = (TAUGSC(I+2)-TAUGSC(I))/ABS(MUO)
          E=T1*(T1-T3)*(T1-T4) ; A(1)=1/E
          B(1)=-(T3+T4)/E    ; C(1)=T3*T4/E ; D(1)=0.0
          E=-T1*T3*T4 ; A(2)=1/E 
          B(2)=-(T1+T3+T4)/E ; C(2)=(T1*T3+T1*T4+T3*T4)/E ; D(2)=1.0
          E=T3*(T3-T1)*(T3-T4) ; A(3)=1/E
          B(3)=-(T1+T4)/E    ; C(3)=T1*T4/E ; D(3)=0.0
          E=T4*(T4-T1)*(T4-T3) ; A(4)=1/E
          B(4)=-(T1+T3)/E    ; C(4)=T1*T3/E ; D(4)=0.0
        ELSE
          T1 = (TAUGSC(I-1)-TAUGSC(I))/ABS(MUO)
          E=T1*(T1-T3) ; B(1)=1/E ; C(1)=-T3/E      ; D(1)=0.0
          E=T1*T3      ; B(2)=1/E ; C(2)=-(T1+T3)/E ; D(2)=1.0
          E=T3*(T3-T1) ; B(3)=1/E ; C(3)=-T1/E      ; D(3)=0.0
        ENDIF
        T = T3
        TR = EXP(-T)
        SRCINT(0,I) = TR
        IF (MUO < 0.0) THEN
          SRCINT(1:4,I) = A(:)*(6*TR+T**3-3*T**2+6*T-6) + D(:)*(1-TR) &
                        + B(:)*(-2*TR+T**2-2*T+2) + C(:)*(TR+T-1)
        ELSE        
          SRCINT(1:4,I) = A(:)*(6-TR*(T**3+3*T**2+6*T+6)) + D(:)*(1-TR) &
                        + B(:)*(2-TR*(T**2+2*T+2)) + C(:)*(1-TR*(T+1)) 
        ENDIF
      ENDDO
    ENDIF
  ENDDO


  IF (MUO < 0.0) THEN
    ! For downwelling direction integrate source function from the top down
    RADVEC(IXP(1,1)) = BNDRAD
    DO LAY = 1, NLAY
      I1 = IXP(1,LAY)
      I2 = IXP(2,LAY)-1
      IF (I1 == I2) THEN
        RADVEC(I1+1) = SRCINT(0,I1)*RADVEC(I1) &
             + SRCINT(1,I1)*SOURCE1(IXG(I1)) + SRCINT(2,I1)*SOURCE1(IXG(I1+1))
      ELSE
        RADVEC(I1+1) = SRCINT(0,I1)*RADVEC(I1) + SRCINT(2,I1)*SOURCE1(IXG(I1)) &
             + SRCINT(3,I1)*SOURCE1(IXG(I1+1)) + SRCINT(4,I1)*SOURCE1(IXG(I1+2))
        DO I = I1+1, I2-1
          RADVEC(I+1) = SRCINT(0,I)*RADVEC(I) + SRCINT(1,I)*SOURCE1(IXG(I-1)) &
               + SRCINT(2,I)*SOURCE1(IXG(I)) + SRCINT(3,I)*SOURCE1(IXG(I+1)) &
               + SRCINT(4,I)*SOURCE1(IXG(I+2))
        ENDDO
        RADVEC(I2+1) = SRCINT(0,I2)*RADVEC(I2) + SRCINT(1,I2)*SOURCE1(IXG(I2-1)) &
             + SRCINT(2,I2)*SOURCE1(IXG(I2)) + SRCINT(3,I2)*SOURCE1(IXG(I2+1))
      ENDIF
      IF (LAY<NLAY) RADVEC(IXP(1,LAY+1)) = RADVEC(IXP(2,LAY))
    ENDDO
  ELSE
    ! For upwelling direction integrate source function from the bottom up
    RADVEC(IXP(2,NLAY)) = BNDRAD
    DO LAY = NLAY, 1, -1
      I1 = IXP(1,LAY)
      I2 = IXP(2,LAY)-1
      IF (I1 == I2) THEN
        RADVEC(I1) = SRCINT(0,I1)*RADVEC(I1+1) &
             + SRCINT(1,I1)*SOURCE1(IXG(I1)) + SRCINT(2,I1)*SOURCE1(IXG(I1+1))
      ELSE
        RADVEC(I2) = SRCINT(0,I2)*RADVEC(I2+1) + SRCINT(1,I2)*SOURCE1(IXG(I2-1)) &
             + SRCINT(2,I2)*SOURCE1(IXG(I2)) + SRCINT(3,I2)*SOURCE1(IXG(I2+1))
        DO I = I2-1, I1+1, -1
          RADVEC(I) = SRCINT(0,I)*RADVEC(I+1) + SRCINT(1,I)*SOURCE1(IXG(I-1)) &
             + SRCINT(2,I)*SOURCE1(IXG(I)) + SRCINT(3,I)*SOURCE1(IXG(I+1)) &
             + SRCINT(4,I)*SOURCE1(IXG(I+2))
        ENDDO
        RADVEC(I1) = SRCINT(0,I1)*RADVEC(I1+1) + SRCINT(2,I1)*SOURCE1(IXG(I1)) &
            + SRCINT(3,I1)*SOURCE1(IXG(I1+1)) + SRCINT(4,I1)*SOURCE1(IXG(I1+2))
      ENDIF
      IF (LAY>1) RADVEC(IXP(2,LAY-1)) = RADVEC(IXP(1,LAY))
    ENDDO
  ENDIF
  DEALLOCATE (SRCINT)
END SUBROUTINE INTEGRATE_SOURCE
 



SUBROUTINE INIT_RADIANCE (NLAY, ML, TAUPSC, ALBEDOSC, LEGENSC, TEMPP, &
                          SRCTYPE, SOLARFLUX, SOLARMU, GNDALB, SFCTEMP, &
                          SKYRAD, WAVELEN, WAVENO,  NPTS, TAUGSC, IXP, NLM, &
                          SHRADIANCE)
 ! Initializes the spherical harmonic radiance field by solving the
 ! plane-parallel Eddington system to give the first two terms of the
 ! spherical harmonic series at the NPTS sublayer interfaces.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, ML, NPTS, IXP(2,NLAY), NLM
  REAL,    INTENT(IN) :: TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML,NLAY)
  REAL,    INTENT(IN) :: TEMPP(NLAY+1)
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, GNDALB, SFCTEMP
  REAL,    INTENT(IN) :: SKYRAD, WAVELEN, WAVENO(2), TAUGSC(NPTS)
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE
  REAL,    INTENT(OUT):: SHRADIANCE(NLM,NPTS)

  INTEGER :: NLAYER, LAY, I, L
  LOGICAL, PARAMETER :: DELTAM = .FALSE.
  REAL    :: PI, C0, C1, GNDEMIS, TAU1, U
  REAL, ALLOCATABLE :: OPTDEPTHS(:), ALBEDOS(:), ASYMMETRIES(:)
  REAL, ALLOCATABLE :: TEMPS(:), FLUXES(:,:)
  CHARACTER(LEN=1) :: SRCTYPE2

  NLAYER = SUM(IXP(2,:)-IXP(1,:))
  ALLOCATE (OPTDEPTHS(NLAYER), ALBEDOS(NLAYER), ASYMMETRIES(NLAYER))
  ALLOCATE (TEMPS(NLAYER+1), FLUXES(3,NLAYER+1))

  SRCTYPE2 = SRCTYPE
  IF (SRCTYPE == 'B') SRCTYPE2='T'
  PI = ACOS(-1.0)
  C0 = SQRT(1.0/PI)
  C1 = SQRT(3.0/(4*PI))

  ! Make layer properties for the Eddington routine
  L = 0
  DO LAY = 1, NLAY
    TAU1 = TAUGSC(IXP(1,LAY))
    DO I = IXP(1,LAY), IXP(2,LAY)-1
      L = L + 1
      OPTDEPTHS(L) = TAUGSC(I+1)-TAUGSC(I)
      ALBEDOS(L) = ALBEDOSC(LAY)
      ASYMMETRIES(L) = LEGENSC(1,LAY)/3.0
      IF (TAUPSC(LAY) > 0.0) THEN
        U = (TAUGSC(I)-TAU1)/TAUPSC(LAY)
      ELSE
        U = (I-IXP(1,LAY))/FLOAT(IXP(2,LAY)-IXP(1,LAY))
      ENDIF
      TEMPS(L) = TEMPP(LAY) + U*(TEMPP(LAY+1)-TEMPP(LAY))
    ENDDO
  ENDDO
  IF (L /= NLAYER) STOP 'INIT_RADIANCE: initial grid error'
  TEMPS(NLAYER+1) = TEMPP(NLAY+1)

     ! Call the Eddington flux routine
  GNDEMIS = 1.0-GNDALB
  CALL EDDRTF (NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES, &
               TEMPS, DELTAM, SRCTYPE2, SOLARFLUX, SOLARMU, &
               SFCTEMP, GNDEMIS, SKYRAD, WAVELEN, WAVENO, FLUXES)


  ! Convert fluxes to first two moments of spherical harmonics
  SHRADIANCE(:,:) = 0.0
  L = 0
  DO LAY = 1, NLAY
    DO I = IXP(1,LAY), IXP(2,LAY)
      L = L + 1
      SHRADIANCE(1,I) = C0*(FLUXES(1,L)+FLUXES(2,L))
      SHRADIANCE(2,I) = C1*(FLUXES(1,L)-FLUXES(2,L))
    ENDDO
    L = L - 1
  ENDDO

  DEALLOCATE (OPTDEPTHS, ALBEDOS, ASYMMETRIES, TEMPS, FLUXES)
END SUBROUTINE INIT_RADIANCE


 
 
SUBROUTINE EDDRTF (NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES, &
                   TEMPS, DELTAM, SRCTYPE, SOLARFLUX, SOLARMU, &
                   SFCTEMP, GNDEMIS, SKYRAD, WAVELEN, WAVENO,  FLUXES)
 ! EDDRTF computes the layer interface fluxes for a plane-parallel
 ! atmosphere with either solar or thermal source of radiation using 
 ! the Eddington approximation.  The medium is specified by a number 
 ! of homogeneous layers.  For a thermal source the Planck function is
 ! linear with optical depth, while for a solar source it is exponential.
 ! The temperatures, optical depth, single scattering albedo, and 
 ! asymmetry parameter are specified for each layer.  The boundary
 ! conditions such as the solar flux, and reflection and/or emission from 
 ! ground surface are also specified. Delta Eddington scaling may be 
 ! used.  The diffuse Eddington fluxes and the solar direct flux at 
 ! each level are returned.
 ! The model works by calculating the reflection, transmission, and
 ! source terms for each layer from the input properties.  A
 ! tri-diagonal matrix solver is then used to compute the diffuse fluxes 
 ! at each layer interface from the applied boundary conditions.
 !
 ! Parameters:
 !   Input:
 ! NLAYER         integer      Number of homogenous layers
 !                              (layers are specified from the top down)
 ! OPTDEPTHS      real array   Optical thickness of layers
 ! ALBEDOS        real array   Single scattering albedos
 ! ASYMMETRIES    real array   Asymmetry parameters
 ! TEMPS          real array   Temperatures (K) at layer interfaces
 !                              (e.g. TEMPS(1) is at top of top layer, 
 !                               TEMPS(2) is at bottom of top layer).
 ! DELTAM         logical      True for delta-Eddington scaling
 ! SRCTYPE        character    'S' for solar source, 'T' for thermal source
 ! SOLARFLUX      real         Incident solar flux on horizonal plane
 ! SOLARMU        real         Cosine of the solar zenith angle
 ! SFCTEMP        real         Surface temperature (Kelvin)
 ! GNDEMIS        real         Ground emissivity (1-albedo)
 ! SKYRAD         real         Radiance (for solar) or brightness 
 !                              temperature (for thermal) of isotropic 
 !                              incident radiation from above
 ! WAVELEN        real         Wavelength (micron)
 ! WAVENO(2)      real         Wavenumber range (cm^-1) (if WAVELEN=0)
 !
 !   Output:
 ! FLUXES         real         Eddington fluxes at layer interfaces.
 !                               FLUXES(1,L) is upwelling diffuse, 
 !                               FLUXES(2,L) is downwelling diffuse,
 !                               FLUXES(3,L) is downwelling direct,
 !                               L=1 is top, L=NLAYER+1 is bottom
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAYER
  LOGICAL, INTENT(IN) :: DELTAM
  REAL,    INTENT(IN) :: TEMPS(NLAYER+1)
  REAL,    INTENT(IN) :: OPTDEPTHS(NLAYER), ALBEDOS(NLAYER), ASYMMETRIES(NLAYER)
  REAL,    INTENT(IN) :: SFCTEMP, GNDEMIS, SKYRAD, SOLARFLUX, SOLARMU
  REAL,    INTENT(IN) :: WAVELEN, WAVENO(2)
  CHARACTER(LEN=1), INTENT(IN) ::  SRCTYPE
  REAL,    INTENT(OUT) :: FLUXES(3,NLAYER+1)
 
  INTEGER :: N, L, I
  DOUBLE PRECISION :: DELTAU, G, OMEGA, F
  DOUBLE PRECISION :: LAMBDA, R, T, D, CP, CM, A, B, X1, X2
  DOUBLE PRECISION :: REFLECT, TRANS, SOURCEP, SOURCEM
  DOUBLE PRECISION :: RADP1P, RADP1M, RADP2P, RADP2M
  DOUBLE PRECISION :: PI, MU0, SKYFLUX,GNDFLUX, PLANCK1,PLANCK2,C,TAU
  DOUBLE PRECISION :: EXLP, EXLM, V, DS, B1, B2, SOLPP, SOLPM
  DOUBLE PRECISION, ALLOCATABLE :: LOWER(:), UPPER(:), DIAG(:), RHS(:)
  PARAMETER (PI=3.1415926535)
  REAL :: PLANCKFUNC

  ! Compute the reflection, transmission, and source coefficients 
  !   for each layer for the diffuse Eddington two stream problem.
  N = 2*NLAYER+2
  ALLOCATE (LOWER(N), UPPER(N), DIAG(N), RHS(N))
  IF (SRCTYPE == 'T') THEN
    PLANCK1 = PI*PLANCKFUNC(WAVELEN,WAVENO,TEMPS(1))
  ENDIF
  MU0 = ABS(SOLARMU)
  TAU = 0.0
  I = 2
  DO L = 1, NLAYER
    DELTAU = OPTDEPTHS(L)
    IF (DELTAU < 0.0) STOP 'EDDRTF: TAU<0'
    ! Special case for zero optical depth
    IF (DELTAU == 0.0) THEN
      TRANS = 1.0
      REFLECT = 0.0
      SOURCEP = 0.0
      SOURCEM = 0.0
    ELSE
      OMEGA = ALBEDOS(L)
      G = ASYMMETRIES(L)
      IF (DELTAM) THEN
        F = G**2
        DELTAU = (1-OMEGA*F)*DELTAU
        OMEGA = (1-F)*OMEGA/(1-OMEGA*F)
        G = (G-F)/(1-F)
      ENDIF
      R = ( 1.0 - OMEGA*(4.0-3.0*G) )/4.0
      T = ( 7.0 - OMEGA*(4.0+3.0*G) )/4.0
      LAMBDA = SQRT( 3.0*(1.0-OMEGA)*(1.0-OMEGA*G) )
      ! Special case for conservative scattering (lambda=0)
      IF (LAMBDA == 0.0) THEN
        D = 1.0/(1.0+T*DELTAU)
        TRANS = D
        REFLECT = -R*DELTAU*D
      ELSE
        X1 = -R
        X2 = LAMBDA + T
        EXLP = DEXP(MIN(LAMBDA*DELTAU,75.D0))
        EXLM = 1.0/EXLP
        TRANS = 2.*LAMBDA/(X2*EXLP + (LAMBDA-T)*EXLM)
        REFLECT = X1*(EXLP - EXLM) *TRANS /(2.*LAMBDA)
        D = 1.0/(X2**2 *EXLP - X1**2 *EXLM)
      ENDIF

      IF (SRCTYPE == 'T') THEN
        ! Calculate thermal source terms
        PLANCK2 = PI*PLANCKFUNC(WAVELEN,WAVENO,TEMPS(L+1))
        V = 2.0*(PLANCK2-PLANCK1)/(3.0*(1.-OMEGA*G)*DELTAU)
        RADP1P = -V + PLANCK1
        RADP2M =  V + PLANCK2
        RADP2P = -V + PLANCK2
        RADP1M =  V + PLANCK1
        IF (LAMBDA .EQ. 0.0) THEN
          A =  (R*DELTAU*RADP1P - RADP2M) *D
          B = -(R*RADP1P + T*RADP2M) *D
          SOURCEP = (B - T*(A+B*DELTAU))/R + RADP2P
          SOURCEM = A + RADP1M
        ELSE
          CP  =  (X1*EXLM*RADP1P - X2*RADP2M) *D
          CM = (-X2*EXLP*RADP1P + X1*RADP2M) *D
          SOURCEP = X1*CP*EXLP + X2*CM*EXLM + RADP2P
          SOURCEM = X2*CP + X1*CM + RADP1M
        ENDIF
        PLANCK1 = PLANCK2
        FLUXES(3,L) = 0.0
      ELSE
        ! Calculate solar source terms
        FLUXES(3,L) = SOLARFLUX*EXP(-TAU/MU0)
        DS = 1.0/(LAMBDA**2-1.0/MU0**2)
        B1 = 0.5*OMEGA*(SOLARFLUX/MU0)*EXP(-TAU/MU0) *DS
        B2 = 0.5*OMEGA*(SOLARFLUX/MU0)*EXP(-(TAU+DELTAU)/MU0) *DS
        SOLPP =  1.0 + 1.5*G*MU0
        SOLPM = -1.0 + 1.5*G*MU0
        RADP1P = ( (T+1.0/MU0)*SOLPP + R*SOLPM )*B1
        RADP2M = ((-T+1.0/MU0)*SOLPM - R*SOLPP )*B2
        RADP2P = ( (T+1.0/MU0)*SOLPP + R*SOLPM )*B2
        RADP1M = ((-T+1.0/MU0)*SOLPM - R*SOLPP )*B1
        IF (LAMBDA .EQ. 0.0) THEN
          A =  (R*DELTAU*RADP1P - RADP2M) *D
          B = -(R*RADP1P + T*RADP2M) *D
          SOURCEP = (B - T*(A+B*DELTAU))/R + RADP2P
          SOURCEM = A + RADP1M
        ELSE
          CP  =  (X1*EXLM*RADP1P - X2*RADP2M) *D
          CM = (-X2*EXLP*RADP1P + X1*RADP2M) *D
          SOURCEP = X1*CP*EXLP + X2*CM*EXLM + RADP2P
          SOURCEM = X2*CP + X1*CM + RADP1M
        ENDIF
        TAU = TAU + DELTAU
      ENDIF
    ENDIF
    DIAG(I) = -REFLECT
    DIAG(I+1) = -REFLECT
    LOWER(I) = 1.0
    LOWER(I+1) = -TRANS
    UPPER(I) = -TRANS
    UPPER(I+1) = 1.0
    RHS(I) = SOURCEM
    RHS(I+1) = SOURCEP
    I = I + 2
  ENDDO

  ! Set up boundary radiances
  IF (SRCTYPE == 'S') THEN
    FLUXES(3,NLAYER+1) = SOLARFLUX*EXP(-TAU/MU0)
    GNDFLUX = (1.0-GNDEMIS)*SOLARFLUX*EXP(-TAU/MU0)
    SKYFLUX = PI*SKYRAD
  ELSE
    FLUXES(3,NLAYER+1) = 0.0
    GNDFLUX = PI*GNDEMIS*PLANCKFUNC(WAVELEN,WAVENO,SFCTEMP)
    SKYFLUX = PI*PLANCKFUNC(WAVELEN,WAVENO,SKYRAD)
  ENDIF
  ! Setup for and call the tri-diagonal matrix solver
  RHS(1) = SKYFLUX
  DIAG(1) = 0.0
  UPPER(1) = 1.0
  DIAG(N) = -(1.0-GNDEMIS)
  LOWER(N) = 1.0
  RHS(N) = GNDFLUX
  CALL TRIDIAG (N, LOWER, DIAG, UPPER, RHS)
  ! Put the fluxes in the output array
  I = 1
  DO L = 1, NLAYER+1 
    FLUXES(1,L) = RHS(I)
    FLUXES(2,L) = RHS(I+1)
    I = I + 2
  ENDDO
 
  DEALLOCATE (LOWER, UPPER, DIAG, RHS)
END SUBROUTINE EDDRTF
 


SUBROUTINE TRIDIAG (N, LOWER, DIAG, UPPER, RHS)
 ! Computes the solution to a tridiagonal system. 
 ! N is order of the matrix.  LOWER(2..N) is the subdiagonal,
 ! DIAG(1..N) is the diagonal, and UPPER(1..N-1) is the 
 ! superdiagonal.  On input RHS is the right hand side, while
 ! on output it is the solution vector.  Everything is destroyed.
 ! Hacked from Linpack DGTSL.
  IMPLICIT NONE
  INTEGER :: N 
  DOUBLE PRECISION :: LOWER(*), DIAG(*), UPPER(*), RHS(*)
  INTEGER :: K, KB
  DOUBLE PRECISION :: T

  IF (N .EQ. 1) THEN
    IF (DIAG(1) .EQ. 0.0) GOTO 990
    RHS(1) = RHS(1)/DIAG(1)
  ENDIF
  LOWER(1) = DIAG(1)
  DIAG(1) = UPPER(1)
  UPPER(1) = 0.0
  UPPER(N) = 0.0
  DO K = 1, N-1
    ! Interchange this and next row to the get the largest pivot.
    IF (ABS(LOWER(K+1)) .GE. ABS(LOWER(K))) THEN
      T = LOWER(K+1)
      LOWER(K+1) = LOWER(K)
      LOWER(K) = T
      T = DIAG(K+1)
      DIAG(K+1) = DIAG(K)
      DIAG(K) = T
      T = UPPER(K+1)
      UPPER(K+1) = UPPER(K)
      UPPER(K) = T
      T = RHS(K+1)
      RHS(K+1) = RHS(K)
      RHS(K) = T
    ENDIF
    IF (LOWER(K) .EQ. 0.0) GOTO 990
    T = -LOWER(K+1)/LOWER(K)
    LOWER(K+1) = DIAG(K+1) + T*DIAG(K)
    DIAG(K+1) = UPPER(K+1) + T*UPPER(K)
    UPPER(K+1) = 0.0
    RHS(K+1) = RHS(K+1) + T*RHS(K)
  ENDDO
  IF (LOWER(N) .EQ. 0.0) GOTO 990
    ! Back substitute
  RHS(N) = RHS(N)/LOWER(N)
  RHS(N-1) = (RHS(N-1) - DIAG(N-1)*RHS(N))/LOWER(N-1)
  DO KB = 1, N-2
    K = N - 2 - KB + 1
    RHS(K) = (RHS(K) -DIAG(K)*RHS(K+1) -UPPER(K)*RHS(K+2))/LOWER(K)
  ENDDO
  RETURN
990 CONTINUE
    STOP 'Singular matrix in TRIDIAG'
END SUBROUTINE TRIDIAG
 

 



SUBROUTINE LEGENDRE_ALL (COSSCAT, NLEG, P)
 ! This subroutine computes a set of Legendre polynomials for
 ! a particular scattering angle COSSCAT.  NLEG is the maximum term.
 ! The Legendre functions evaluated at COSSCAT are returned in 
 ! P, starting at l=0 and ending with l=NLEG  (NLEG+1 terms).
  IMPLICIT NONE
      INTEGER NLEG
      DOUBLE PRECISION COSSCAT, P(0:NLEG)
      INTEGER L
      DOUBLE PRECISION X, PL, PL1, PL2

      X = DBLE(COSSCAT)
      IF (X*X .GT. 1.) STOP 'LEGENDRE_ALL: |COSSCAT| larger than 1'
      ! Use the stable upward recursion on l, starting from P_0
      PL2 = 1.0D0
      P(0) = PL2
      IF (NLEG .GT. 1) THEN
        PL1 = X
        P(1) = X
      ENDIF
      DO L = 2, NLEG
        PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
        P(L) = PL
        PL2 = PL1
        PL1 = PL
      ENDDO
END SUBROUTINE LEGENDRE_ALL




SUBROUTINE YLMALL (MU, PHI, ML, MM, NCS, P)
 ! This subroutine computes a set of normalized spherical harmonic 
 ! functions, P(J), for a particular direction mu,phi. 
 ! ML is the maximum meridional mode, MM is the maximum azimuthal mode,
 ! and NCS is the azimuthal mode flag (|NCS|=1 for cosine only, |NCS|=2 for 
 ! sines and cosines).  Returns normalized associated Legendre functions 
 ! only if NCS<0. The set is returned for triangular truncation: 
 ! J = NCS*(L*(L+1))/2 + M+1  for L<=MM
 ! J = (NCS*MM+1)*L-MM*(2+NCS*(MM-1))/2 + M+1  for L>MM
  IMPLICIT NONE
      INTEGER ML, MM, NCS
      REAL    MU, PHI, P(*)
      INTEGER J, L, M, C
      DOUBLE PRECISION X, Y, A, PMM, PL, PL1, PL2, PHI8

      C = ABS(NCS)
      IF (C .NE. 1 .AND. C .NE. 2)  STOP 'YLMALL: bad NCS'
      IF (MM .GT. ML)  STOP 'YLMALL: MM greater than LM'
      IF (MU*MU .GT. 1.) STOP 'YLMALL: |MU| larger than 1'
      X = DBLE(MU)
      Y = SQRT(1.0D0-X*X)
      ! Use the stable upward recursion on l, starting from P^m_m
      ! Put in the spherical harmonic normalization as it goes
      PMM = 1.0D0/SQRT(4.0D0*ACOS(-1.0D0))
      DO M = 0, MM
        IF (M .GT. 0)  PMM = -PMM*Y*SQRT((2*M+1.0D0)/(2.0D0*M))
        J = C*(M*(M+1))/2 + M+1
        P(J) = PMM
        PL2 = PMM
        IF (M .LT. ML) THEN
          IF (M+1.LE.MM) J=C*((M+1)*(M+2))/2 +M+1
          IF (M+1.GT.MM) J=(C*MM+1)*(M+1)-(MM*(2+C*(MM-1)))/2+M+1
          PL1 = SQRT(2*M+3.0D0)*X*PMM
          P(J) = PL1
        ENDIF
        DO L = M+1, ML-1
          IF (L+1.LE.MM) J=C*((L+1)*(L+2))/2 +M+1
          IF (L+1.GT.MM) J=(C*MM+1)*(L+1)-(MM*(2+C*(MM-1)))/2+M+1
          A = 1.0D0/((L+M+1.D0)*(L-M+1.D0))
          PL = SQRT((2*L+1)*A*(2*L+3)) *X*PL1 &
                  - SQRT((2*L+3)*A*(L+M)*(L-M)/(2*L-1.)) *PL2
          P(J) = PL
          PL2 = PL1
          PL1 = PL
        ENDDO
        IF (M .EQ. 0) PMM = PMM*SQRT(2.0D0)
      ENDDO
      ! If there are M<0 terms then fill them in
      IF (C .EQ. 2) THEN
        DO L = 0, ML
          DO M = 1, MIN(L,MM)
            IF (L .LE. MM) J = L*(L+1) +M+1
            IF (L .GT. MM) J = MM*(2*L-MM) +L+M+1
            P(J-2*M) = P(J)
          ENDDO
        ENDDO
      ENDIF
      ! Put in the azimuthal dependence
      PHI8 = DBLE(PHI)
      IF (NCS .GT. 0) THEN
        DO M = (1-NCS)*MM, MM
          IF (M .LT. 0) THEN
            A = SIN(-M*PHI8)
          ELSE IF (M .GT. 0) THEN
            A = COS(M*PHI8)
          ELSE
            A = 1.0D0
          ENDIF
          DO L = ABS(M), ML
            IF (L .LE. MM) THEN
              J = C*(L*(L+1))/2 +M+1
            ELSE
              J = (C*MM+1)*L-(MM*(2+C*(MM-1)))/2 + M+1
            ENDIF
            P(J) = A*P(J)
          ENDDO
        ENDDO
      ENDIF
END SUBROUTINE YLMALL


 




SUBROUTINE GAUSQUADS (N, XA, WT)
 ! Generates the abscissas (X) and weights (W) for an N point
 ! Gauss-Legendre quadrature.  The XA are returned in this order: 
 ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK, 
 ! K=N/2, N must be even).
  IMPLICIT NONE
      INTEGER  N
      REAL     XA(*), WT(*)
      INTEGER  K, I, J, L
      DOUBLE PRECISION  X, XP, PL, PL1, PL2, DPL, TINY
      PARAMETER (TINY=3.0D-7)

      K = N/2
      IF (2*K .NE. N) STOP 'GAUSQUADS: N must be even'
      DO J = 1, K
        X = COS(3.141592654*(J-.25)/(N+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
          ENDDO
          DPL = N*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        XA(J)     = -X
        XA(J+K) = X
        WT(J  )   = 2.0D0/((1.0D0-X*X)*DPL*DPL)
        WT(J+K) = WT(J)
      ENDDO
END SUBROUTINE GAUSQUADS




SUBROUTINE DGAUSQUADS (N, XA, WT)
 ! Generates the abscissas (X) and weights (W) for an N point
 ! Double-Gauss-Legendre quadrature.  The XA are returned in this order: 
 ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK, 
 ! K=N/2, N must be even).
  IMPLICIT NONE
      INTEGER  N
      REAL     XA(*), WT(*)
      INTEGER  K, I, J, L, N2
      DOUBLE PRECISION  X, XP, PL, PL1, PL2, DPL, TINY
      PARAMETER (TINY=3.0D-7)

      N2 = N/2
      IF (2*N2 .NE. N) STOP 'DGAUSQUADS: N must be even'
      K = (N2+1)/2
      DO J = 1, K
        X = COS(3.141592654*(J-.25)/(N2+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO L = 2, N2
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
          ENDDO
          DPL = N2*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        XA(J+N2) = (1+X)/2
        XA(N2+1-J+N2) = (1-X)/2
        WT(J+N2) = 1/((1-X*X)*DPL*DPL)
        WT(N2+1-J+N2) = WT(J+N2)
        XA(J) = -XA(J+N2)
        XA(N2+1-J) = -XA(N2+1-J+N2)
        WT(J) = WT(J+N2)
        WT(N2+1-J) = WT(J+N2)
      ENDDO
END SUBROUTINE DGAUSQUADS




SUBROUTINE INTEGRATE_PLANCK (X1, X2, F)
 ! Returns integral of x**3/(exp(x)-1) from x1 to x2.  Accurate to better 
 ! than 1 part in 10**9 for integral from 0 to infinity (pi**4/15).
  IMPLICIT NONE
  DOUBLE PRECISION :: X1, X2, F
  DOUBLE PRECISION :: INTLOW, INTHIGH, C

  C = 1.0D0
  IF (X1 < C .AND. X2 < C) THEN
    F = INTLOW(X2) - INTLOW(X1)
  ELSE IF (X1 < C .AND. X2 >= C) THEN
    F = INTLOW(C) - INTLOW(X1) + INTHIGH(C) - INTHIGH(X2)
  ELSE IF (X1 >= C .AND. X2 >= C) THEN
    F = INTHIGH(X1) - INTHIGH(X2)
  ELSE
    STOP 'X1 and X2 out of order'
  ENDIF
END SUBROUTINE INTEGRATE_PLANCK


FUNCTION INTLOW (X)
 ! Integral of x**3/(exp(x)-1) from 0 to x.  Accurate for x less than about 1.
 ! Uses Taylor series expansion around x=0.
  IMPLICIT NONE 
  DOUBLE PRECISION :: X, INTLOW
  INTEGER :: N
  DOUBLE PRECISION :: SUM, F, A(29)
  DATA    A/0.0D0, 0.0D0,  0.333333333333333333333D0, &
    -0.125D0,  0.016666666666666666667D0, 0.0D0, -0.0001984126984126984127D0, &
     0.0D0,  0.36743092298647854203D-5,   0.0D0, -0.75156325156325156325D-7, &
     0.0D0,  0.16059043836821614599D-8,   0.0D0, -0.35227934257916621232D-10, &
     0.0D0,  0.78720803121674581370D-12,  0.0D0, -0.17840422612224120352D-13, &
     0.0D0,  0.40886009791799259829D-15,  0.0D0, -0.94559508632959211871D-17, &
     0.0D0,  0.22036011313440918061D-18,  0.0D0, -0.51683202540046382743D-20, &
     0.0D0,  0.12188644964239543006D-20/

  SUM = A(4)*X**4
  F = X
  DO N = 3, 29, 2
    F = X*X*F
    SUM = SUM + A(N)*F
  ENDDO
  INTLOW = SUM
END FUNCTION INTLOW

FUNCTION INTHIGH (X)
 ! Integral of x**3/(exp(x)-1) from x to infinity.
 ! Accurate for x greater than about 1.
  IMPLICIT NONE
  DOUBLE PRECISION :: X, INTHIGH
  INTEGER :: N
  DOUBLE PRECISION :: SUM

  SUM = 0.0D0
  DO N = 1, 15
    SUM = SUM + EXP(-N*X) *(X**3/N + 3*X**2/N**2 + 6*X/N**3 + 6.0D0/N**4)
  ENDDO
  INTHIGH = SUM
END FUNCTION INTHIGH




      SUBROUTINE SSORT (X, Y, N, KFLAG)
! ***BEGIN PROLOGUE  SSORT
! ***PURPOSE  Sort an array and optionally make the same interchanges in
!             an auxiliary array.  The array may be sorted in increasing
!             or decreasing order.  A slightly modified QUICKSORT
!             algorithm is used.
! ***LIBRARY   SLATEC
! ***CATEGORY  N6A2B
! ***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
! ***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
! ***AUTHOR  Jones, R. E., (SNLA)
!            Wisniewski, J. A., (SNLA)
! ***DESCRIPTION
! 
!    SSORT sorts array X and optionally makes the same interchanges in
!    array Y.  The array X may be sorted in increasing order or
!    decreasing order.  A slightly modified quicksort algorithm is used.
! 
!    Description of Parameters
!       X - array of values to be sorted   (usually abscissas)
!       Y - array to be (optionally) carried along
!       N - number of values in array X to be sorted
!       KFLAG - control parameter
!             =  2  means sort X in increasing order and carry Y along.
!             =  1  means sort X in increasing order (ignoring Y)
!             = -1  means sort X in decreasing order (ignoring Y)
!             = -2  means sort X in decreasing order and carry Y along.
! 
! ***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                  for sorting with minimal storage, Communications of
!                  the ACM, 12, 3 (1969), pp. 185-187.
! ***END PROLOGUE  SSORT
!      .. Scalar Arguments ..
      INTEGER KFLAG, N
!      .. Array Arguments ..
!       REAL X(*), Y(*)
      REAL X(*)
      INTEGER Y(*)
!      .. Local Scalars ..
!       REAL R, T, TT, TTY, TY
      REAL R, T, TT
      INTEGER TY, TTY
      INTEGER I, IJ, J, K, KK, L, M, NN
!      .. Local Arrays ..
      INTEGER IL(51), IU(51)
!      .. External Subroutines ..
!      .. Intrinsic Functions ..
      INTRINSIC ABS, INT
! ***First executable statement  SSORT
      NN = N
      IF (NN .LT. 1) THEN
         STOP 'The number of values to be sorted is not positive.'
      ENDIF

      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
        STOP 'The sort control parameter, K, is not 2, 1, -1, or -2.'
      ENDIF

!      Alter array X to get decreasing order if needed
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF

      IF (KK .EQ. 2) GO TO 100

!      Sort X only
      M = 1
      I = 1
      J = NN
      R = 0.375E0

   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF

   30 K = I

!      Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)

!      If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J

!      If last element of array is less than than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)

!         If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF

!      Find an element in the second half of the array which is smaller
!      than T
   40 L = L-1
      IF (X(L) .GT. T) GO TO 40

!      Find an element in the first half of the array which is greater
!      than T
   50 K = K+1
      IF (X(K) .LT. T) GO TO 50

!      Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF

!      Save upper and lower subscripts of the array yet to be sorted
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70

!      Begin again on another portion of the unsorted array
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1

   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = X(I+1)
      IF (X(I) .LE. T) GO TO 80
      K = I

   90 X(K+1) = X(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 90
      X(K+1) = T
      GO TO 80

!      Sort X and carry Y along
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0

  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
! 
  120 K = I

!      Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)

!      If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)
      ENDIF
      L = J

!      If last element of array is less than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)

!         If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
         ENDIF
      ENDIF

!      Find an element in the second half of the array which is smaller
!      than T
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130

!      Find an element in the first half of the array which is greater
!      than T
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140

!      Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY
         GO TO 130
      ENDIF

!      Save upper and lower subscripts of the array yet to be sorted
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160

!      Begin again on another portion of the unsorted array
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1

  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)
      TY = Y(I+1)
      IF (X(I) .LE. T) GO TO 170
      K = I

  180 X(K+1) = X(K)
      Y(K+1) = Y(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y(K+1) = TY
      GO TO 170

!      Clean up
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
      RETURN
      END



SUBROUTINE ADD_IN_KDIST (NLAY, HEIGHTP, TAUP, ALBEDOP, &
                         NZCKD, ZCKD, GASABS,  TAUPM, ALBEDOPM)
 ! Adds the gaseous absorption in with the particle optical properties.
 ! The gaseous absorption optical depth is found for each of the NLAY
 ! optical property layers by interpolated the absorption linearly
 ! between the levels in ZCKD.
  INTEGER, INTENT(IN) :: NLAY, NZCKD
  REAL,    INTENT(IN) :: HEIGHTP(NLAY+1), TAUP(NLAY), ALBEDOP(NLAY)
  REAL,    INTENT(IN) :: ZCKD(NZCKD), GASABS(NZCKD)
  REAL,    INTENT(OUT) :: TAUPM(NLAY), ALBEDOPM(NLAY)
  INTEGER :: LAY, I
  REAL    :: GASTAU, Z, Z1, Z2, F, ABS

  IF (NZCKD <= 0) THEN
    TAUPM(1:NLAY) = TAUP(1:NLAY)
    ALBEDOPM(1:NLAY) = ALBEDOP(1:NLAY)
  ELSE
    IF (HEIGHTP(1) > ZCKD(1) .OR. HEIGHTP(NLAY+1) < ZCKD(NZCKD)) THEN
      WRITE (*,*) 'Layer heights outside of CKD file height range.'
      STOP
    ENDIF
    I = 1
    DO LAY = 1, NLAY
      ! First, use trapezoidal integration to get optical depth of 
      ! gaseous absorption assuming linear interpolation.
      GASTAU = 0.0
      Z1 = HEIGHTP(LAY)
      Z2 = HEIGHTP(LAY+1)
      DO WHILE (ZCKD(I) >= Z1)
        I = I + 1
      ENDDO
      I = I - 1
      Z = Z1
      DO WHILE (ZCKD(I) > Z2)
        F = (Z-ZCKD(I))/(ZCKD(I+1)-ZCKD(I))
        ABS = (1-F)*GASABS(I) + F*GASABS(I+1)
        GASTAU = GASTAU + 0.5*(ABS+GASABS(I+1))*(Z-ZCKD(I+1))
        Z = ZCKD(I+1)
        I = I + 1
      ENDDO
      I = I - 1
      F = (Z2-ZCKD(I))/(ZCKD(I+1)-ZCKD(I))
      ABS = (1-F)*GASABS(I) + F*GASABS(I+1)
      GASTAU = GASTAU - 0.5*(ABS+GASABS(I+1))*(Z2-ZCKD(I+1))
      ! Then, add the gaseous optical depth to the particle properties
      TAUPM(LAY) = TAUP(LAY) + GASTAU
      IF (TAUPM(LAY) > 0.0) THEN
        ALBEDOPM(LAY) = ALBEDOP(LAY)*TAUP(LAY)/TAUPM(LAY)
      ELSE
        ALBEDOPM(LAY) = 0.0
      ENDIF
      ! print *, 'add_kdist: ',LAY,GASTAU,TAUPM(LAY),ALBEDOPM(LAY)
    ENDDO
  ENDIF
END SUBROUTINE ADD_IN_KDIST



SUBROUTINE SUM_KDIST_OUTPUT (IG, WT, NLAY, TAUP, IXP, NPTS, TAUG, &
                             HEIGHTP, HEIGHTG, FLUXUPM, FLUXUP, &
                             FLUXDNM, FLUXDN, MEANRADM, MEANRAD, &
                             NZOUT, NMUOUT, NPHIOUT, RADOUTM, RADOUT)
 ! Sums the output arrays over the k-distribution.  Only sums the output
 ! at the input layer interfaces because the adaptive grid changes for
 ! each k.  Also fixes up the tau grid and height grid accordingly.
  INTEGER, INTENT(IN) :: IG, NLAY, NZOUT, NMUOUT, NPHIOUT
  INTEGER, INTENT(INOUT) :: NPTS, IXP(2,NLAY)
  REAL,    INTENT(IN) :: WT, TAUP(NLAY), HEIGHTP(NLAY+1)
  REAL,    INTENT(INOUT) :: TAUG(*), HEIGHTG(*)
  REAL,    INTENT(IN) :: FLUXUPM(*), FLUXDNM(*), MEANRADM(*)
  REAL,    INTENT(IN) :: RADOUTM(NZOUT,NMUOUT,NPHIOUT)
  REAL,    INTENT(OUT) :: FLUXUP(*), FLUXDN(*), MEANRAD(*)
  REAL,    INTENT(OUT) :: RADOUT(NZOUT,NMUOUT,NPHIOUT)

  IF (IG == 1) THEN
    FLUXUP(1:NLAY+1) = 0.0
    FLUXDN(1:NLAY+1) = 0.0
    MEANRAD(1:NLAY+1) = 0.0
    RADOUT(:,:,:) = 0.0
  ENDIF
  
  DO L = 1, NLAY
    FLUXUP(L) = FLUXUP(L) + WT*FLUXUPM(IXP(1,L))
    FLUXDN(L) = FLUXDN(L) + WT*FLUXDNM(IXP(1,L))
    MEANRAD(L) = MEANRAD(L) + WT*MEANRADM(IXP(1,L))
  ENDDO
  FLUXUP(NLAY+1) = FLUXUP(NLAY+1) + WT*FLUXUPM(IXP(2,NLAY))
  FLUXDN(NLAY+1) = FLUXDN(NLAY+1) + WT*FLUXDNM(IXP(2,NLAY))
  MEANRAD(NLAY+1) = MEANRAD(NLAY+1) + WT*MEANRADM(IXP(2,NLAY))
  RADOUT(:,:,:) = RADOUT(:,:,:) + WT*RADOUTM(:,:,:)

  HEIGHTG(1:NLAY+1) = HEIGHTP(1:NLAY+1)
  TAUG(1) = 0.0
  DO L = 1, NLAY
    TAUG(L+1) = TAUG(L) + TAUP(L)
    IXP(1,L) = L
    IXP(2,L) = L+1
  ENDDO
  NPTS = NLAY+1
END SUBROUTINE SUM_KDIST_OUTPUT




SUBROUTINE USER_INPUT (PROPFILE, CKDFILE, KDIST, NMU, NPHI, &
                       SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                       SFCTEMP, SFCTYPE, SFCPARMS, WAVELEN, WAVENO, &
                       ACCELFLAG, SOLACC, MAXITER, SPLITACC, &
                       FLUXFILE, FLUXPARM, RADFILE, MAXPAR, NZOUT, ZOUT, &
                       NMUOUT, MUOUT, NPHIOUT, PHIOUT)
!       Obtains the input parameters for the program by writing prompts
!     to and reading responses from the terminal.  Echos the inputs to
!     make a useful log file when running non-interactively.  See the overall
!     program documentation for the list of input parameters.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: MAXPAR
  INTEGER, INTENT(OUT) :: NMU, NPHI
  INTEGER, INTENT(OUT) :: MAXITER, FLUXPARM, NZOUT, NMUOUT, NPHIOUT
  LOGICAL, INTENT(OUT) :: KDIST, ACCELFLAG
  REAL,    INTENT(OUT) :: SOLARFLUX, SOLARMU
  REAL,    INTENT(OUT) :: SFCTEMP, SKYRAD, WAVELEN, WAVENO(2)
  REAL,    INTENT(OUT) :: SOLACC, SPLITACC
  REAL,    INTENT(OUT) :: SFCPARMS(MAXPAR)
  REAL,    INTENT(OUT) :: ZOUT(MAXPAR), MUOUT(MAXPAR), PHIOUT(MAXPAR)
  CHARACTER(LEN=1),  INTENT(OUT) :: SRCTYPE, SFCTYPE
  CHARACTER(LEN=80), INTENT(OUT) :: PROPFILE, CKDFILE, FLUXFILE, RADFILE 

  WRITE (*,*)
  WRITE (*,'(1X,A)') 'Atmospheric optical properties file name'
  READ (*,'(A)') PROPFILE
  WRITE (*,*) PROPFILE

  WRITE (*,'(1X,A)') 'Correlated k-distribution file name (or NONE)'
  READ (*,'(A)') CKDFILE
  WRITE (*,*) CKDFILE
  KDIST = CKDFILE(1:2) .NE. 'NO'

  WRITE (*,'(1X,A)') 'Number of discrete ordinates in mu and phi'
  READ (*,*) NMU, NPHI
    WRITE (*,*) NMU, NPHI

  WRITE (*,'(1X,A)') 'Thermal, solar, or both source (T, S, B)'
  READ (*,'(A)') SRCTYPE
    WRITE (*,*) SRCTYPE
  IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
    WRITE (*,'(1X,A)') 'Solar flux and direction (F, mu0)'
    READ (*,*) SOLARFLUX, SOLARMU
    WRITE (*,*) SOLARFLUX, SOLARMU
    WRITE (*,'(1X,A)') 'Isotropic sky radiance'
    READ (*,*) SKYRAD
      WRITE (*,*) SKYRAD
    SOLARMU = -ABS(SOLARMU)
  ENDIF
  IF (SRCTYPE .EQ. 'T' .OR. SRCTYPE .EQ. 'B') THEN
    WRITE (*,'(1X,A)') 'Surface temperature'
    READ (*,*) SFCTEMP
      WRITE (*,*) SFCTEMP
  ENDIF
  IF (SRCTYPE .EQ. 'T') THEN
    WRITE (*,'(1X,A)') 'Sky temperature'
    READ (*,*) SKYRAD
      WRITE (*,*) SKYRAD
  ENDIF

  WRITE (*,'(1X,A)') 'Surface type (L-Lambertian, R-RPV vegetation BRDF)'
  READ (*,'(A1)') SFCTYPE
    WRITE (*,*) SFCTYPE
  IF (SFCTYPE == 'L') THEN
    WRITE (*,'(1X,A)') 'Surface albedo'
    READ (*,*) SFCPARMS(1)
      WRITE (*,*) SFCPARMS(1)
  ELSE IF (SFCTYPE == 'R') THEN
    WRITE (*,'(1X,A)') 'RPV parameters (rho0, k, Theta)'
    READ (*,*) SFCPARMS(1:3)
      WRITE (*,*) SFCPARMS(1:3)
  ELSE
    WRITE (*,*) 'Unknown surface type'
    STOP
  ENDIF

  WRITE (*,'(1X,A)') 'Wavelength (micron) (<=0 for wavenumber range)'
  READ (*,*) WAVELEN
    WRITE (*,*) WAVELEN
  IF (KDIST .AND. WAVELEN > 0.0) THEN
    WRITE (*,*) 'Must use wavenumber range for k-distribution'
    STOP
  ENDIF
  IF (WAVELEN <= 0.0) THEN
    WRITE (*,'(1X,A)') 'Wavenumber range (cm^-1)'
    READ (*,*) WAVENO(1:2)
      WRITE (*,*) WAVENO(1:2)
  ENDIF

  WRITE (*,'(1X,A)') 'Layer splitting accuracy'
  READ (*,*) SPLITACC
    WRITE (*,*) SPLITACC

  WRITE (*,'(1X,A,A)') 'Solution accuracy'
  READ (*,*) SOLACC
    WRITE (*,*) SOLACC
 
  WRITE (*,'(1X,A,A)') 'Acceleration flag (T or F) and max iterations'
  READ (*,*) ACCELFLAG, MAXITER
    WRITE (*,*) ACCELFLAG, MAXITER

  WRITE (*,'(1X,A,A)') 'Flux output file name'
  READ (*,'(A)') FLUXFILE
    WRITE (*,'(A)') FLUXFILE
  WRITE (*,'(1X,A)') 'Fluxes output at input layer boundaries or adaptive levels (1,2)'
  READ (*,*) FLUXPARM
    WRITE (*,*) FLUXPARM

  WRITE (*,'(1X,A,A)') 'Radiance output file name'
  READ (*,'(A)') RADFILE
    WRITE (*,'(A)') RADFILE

  WRITE (*,*) 'Number of heights for output radiance'
  READ (*,*) NZOUT
    WRITE (*,*) NZOUT
  WRITE (*,*) 'Height levels for output radiance'
  READ (*,*) ZOUT(1:NZOUT)
    WRITE (*,*) ZOUT(1:NZOUT)

  WRITE (*,*) 'Number of cosine zenith angles for output radiance'
  READ (*,*) NMUOUT
    WRITE (*,*) NMUOUT
  WRITE (*,*) 'Cosine zenith angles for output radiance'
  READ (*,*) MUOUT(1:NMUOUT)
    WRITE (*,*) MUOUT(1:NMUOUT)

  WRITE (*,*) 'Number of azimuth angles for output radiance'
  READ (*,*) NPHIOUT
    WRITE (*,*) NPHIOUT
  WRITE (*,*) 'Azimuth angles for output radiance (degrees)'
  READ (*,*) PHIOUT(1:NPHIOUT)
    WRITE (*,*) PHIOUT(1:NPHIOUT)

  WRITE (*,*)
END SUBROUTINE USER_INPUT


 

SUBROUTINE READ_PROPERTY_SIZE (PROPFILE, NLAY, MAXLEG)
!       Reads parts of the property file to get the maximum array sizes
!     needed for allocatable arrays.  
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: NLAY, MAXLEG
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE
  INTEGER :: I, NLEG
  REAL    :: TAU, Z, TTOP, TBOT, SSALB
 
  ! Open the file, get the number of layers and the max number Legendre terms
  OPEN (UNIT=1, FILE=PROPFILE, STATUS='OLD')
  READ (1,*) NLAY
  READ (1,*) Z, TTOP
  MAXLEG = 1
  DO I = 1, NLAY
    READ (1,*) Z, TBOT, TAU, SSALB, NLEG
    MAXLEG = MAX(MAXLEG,NLEG)
  ENDDO
  CLOSE (1)
END SUBROUTINE READ_PROPERTY_SIZE



SUBROUTINE READ_PROPERTIES (PROPFILE, NLAY, MAXLEG, &
                            HEIGHTP, TEMPP, TAUP, ALBEDOP, NLEGP, LEGENP)
 ! Reads the plane-parallel medium properties from the file into the
 ! property arrays.  The optical properties are uniform within each layer.
 ! The heights must decrease (top of atmosphere downward).
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NLAY, MAXLEG
  INTEGER, INTENT(OUT) :: NLEGP(NLAY)
  REAL,    INTENT(OUT) :: TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG,NLAY)
  REAL,    INTENT(OUT) :: HEIGHTP(NLAY+1), TEMPP(NLAY+1)
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE
  INTEGER :: NL, I, L

  ! Open the file, get the number of layers and the max number Legendre terms
  OPEN (UNIT=1, FILE=PROPFILE, STATUS='OLD')
  READ (1,*) NL
  READ (1,*) HEIGHTP(1), TEMPP(1)
  DO I = 1, NLAY
    LEGENP(:,I) = 0.0
    READ (1,*) HEIGHTP(I+1), TEMPP(I+1), TAUP(I), ALBEDOP(I), &
               NLEGP(I), (LEGENP(L,I), L=1,NLEGP(I))
    IF (HEIGHTP(I+1) > HEIGHTP(I)) THEN
      WRITE (*,*) 'Heights must decrease in property file.'
      STOP
    ENDIF
    LEGENP(0,I) = 1.0
  ENDDO
  CLOSE (1)
END SUBROUTINE READ_PROPERTIES
 


SUBROUTINE READ_CKD_SIZE (CKDFILE, WAVENO, NG, NZCKD)
 ! Find the size of the k-distribution arrays for this wavenumber band.
  INTEGER, INTENT(OUT) :: NG, NZCKD
  REAL,    INTENT(IN)  :: WAVENO(2)
  CHARACTER(LEN=*), INTENT(IN) :: CKDFILE
  INTEGER :: NB, IB, JB, KB, N
  REAL    :: WAVENUM1, WAVENUM2, SF
      
  OPEN (UNIT=1, FILE=CKDFILE, STATUS='OLD')
  READ (1,*)
  READ (1,*) NB
  READ (1,*)   
  ! Read band information until the right one is found
  JB = 0
  DO IB = 1, NB
    READ (1,*) KB, WAVENUM1, WAVENUM2, SF, N
    IF ( ABS(WAVENUM1-WAVENO(1)) < 1.0 .AND. &
         ABS(WAVENUM2-WAVENO(2)) < 1.0 ) THEN
      BACKSPACE (1)
      NG = N
      READ (1,*) JB, WAVENUM1, WAVENUM2, SF, N
    ENDIF
  ENDDO  
  IF (JB == 0) THEN
    WRITE (*,*) ' Wavenumber range not found in CKD file : ', &
                WAVENO(1), WAVENO(2), CKDFILE
    STOP
  ENDIF 
  READ (1,*) NZCKD
  CLOSE (1)
END SUBROUTINE READ_CKD_SIZE
               
               
               
SUBROUTINE READ_CKD (CKDFILE, WAVENO, NG, NZCKD, SOLFLUX, DELG, ZCKD, KABS)
 ! Reads the information appropriate for one band from a
 ! correlated k-distribution file.   The number of "g"'s are input in NG,
 ! and the number of levels are input in NZCKD. The wavenumber range (cm^-1)
 ! in the file must match the desired in WAVENO.  The solar flux 
 ! is returned in SOLFLUX, the weights or delta g's in DELG, the
 ! levels from the top down in ZCKD, and the absorption coefficients
 ! for each g and level in KABS. 
  INTEGER, INTENT(IN) :: NG, NZCKD
  REAL,    INTENT(INOUT) :: WAVENO(2)
  REAL,    INTENT(OUT) :: SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD,NG)
  CHARACTER(LEN=*), INTENT(IN) :: CKDFILE
  INTEGER :: NB, IB, JB, KB, IG, I, J, N
  REAL    :: WAVENUM1, WAVENUM2, SF
      
  OPEN (UNIT=1, FILE=CKDFILE, STATUS='OLD')
  READ (1,*)
  READ (1,*) NB
  READ (1,*)   
  ! Read band information until the right one is found
  JB = 0
  DO IB = 1, NB
    READ (1,*) KB, WAVENUM1, WAVENUM2, SF, N   
    IF ( ABS(WAVENUM1-WAVENO(1)) < 1.0 .AND. &
         ABS(WAVENUM2-WAVENO(2)) < 1.0 ) THEN
      BACKSPACE (1)
      READ (1,*) JB, WAVENO(1), WAVENO(2), SOLFLUX, &
                 N, (DELG(IG), IG=1, NG)
    ENDIF
  ENDDO
  ! Read the Z level heights
  READ (1,*) N
  READ (1,*)
  READ (1,*)
  DO I = 1, NZCKD
    READ (1,*) ZCKD(I)
  ENDDO
  ! Skip over the irrelevant bands and read the absorption coefficients
  READ (1,*)
  DO I = 1, (JB-1)*NZCKD
    READ (1,*)
  ENDDO
  DO I = 1, NZCKD
    READ (1,*) IB, J, (KABS(I,IG), IG = 1, NG)
    DO IG = 1, NG
      IF (KABS(I,IG) < 0.0) THEN
        WRITE (*,'(1X,A,A,I3,A,I3,A,I2)') 'Negative absorption in CKD file. ', &
            'Band: ',IB, '   Level: ',I, '   k: ',IG
      ENDIF
    ENDDO  
  ENDDO    
  CLOSE (1)
END SUBROUTINE READ_CKD




SUBROUTINE OUTPUT_RESULTS (NLAY, NPTS, NMU, NPHI, NG, &
                           PROPFILE, CKDFILE, SRCTYPE, SOLARFLUX, SOLARMU, &
                           SKYRAD, SFCTYPE, SFCTEMP, SFCPARMS, WAVELEN, WAVENO, & 
                           SPLITACC, SOLACC, MAXITER, ITER, &
                           IXP, HEIGHTG, FLUXUP, FLUXDN, MEANRAD, FLUXPARM, &
                           NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                           RADOUT, OUTTYPE, OUTFILE)
!       Writes the desired type of output file from the output fields
!     with a header giving the input parameters.
!     There are two types (OUTTYPE) of output: 'R' - for radiance,  
!     'F' - for hemispheric/actinic flux.
!     See the overall program documentation for output parameters.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, NPTS, NMU, NPHI, NG
  INTEGER, INTENT(IN) :: NZOUT, NMUOUT, NPHIOUT, FLUXPARM
  INTEGER, INTENT(IN) :: MAXITER, ITER, IXP(2,NLAY)
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, SFCTEMP, SFCPARMS(*), SKYRAD
  REAL,    INTENT(IN) :: SOLACC, SPLITACC, WAVELEN, WAVENO(2)
  REAL,    INTENT(IN) :: HEIGHTG(NPTS)
  REAL,    INTENT(IN) :: FLUXUP(NPTS), FLUXDN(NPTS), MEANRAD(NPTS)
  REAL,    INTENT(IN) :: ZOUT(NZOUT), MUOUT(NMUOUT), PHIOUT(NPHIOUT)
  REAL,    INTENT(IN) :: RADOUT(NZOUT,NMUOUT,NPHIOUT)
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE, SFCTYPE, OUTTYPE
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE, CKDFILE, OUTFILE
  INTEGER :: I, J, K, L, N
  REAL    :: FOURPI
  CHARACTER(LEN=32) :: SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

  FOURPI=4*ACOS(-1.0)
  IF (SRCTYPE == 'S') THEN
    SOURCENAME = 'SOLAR'
  ELSE IF (SRCTYPE == 'T') THEN
    SOURCENAME = 'THERMAL'
  ELSE IF (SRCTYPE == 'B') THEN
    SOURCENAME = 'SOLAR/THERMAL'
  ENDIF

  IF (SFCTYPE == 'L') THEN 
    SFCNAME = 'LAMBERTIAN'
  ELSE IF (SFCTYPE == 'R') THEN 
    SFCNAME = 'RPV-BRDF'
  ENDIF

  IF (WAVELEN > 0.0) THEN
    IF (OUTTYPE == 'R') THEN
      UNITSNAME = 'WATTS/(M^2 MICRON STER)'
    ELSE
      UNITSNAME = 'WATTS/(M^2 MICRON)'
    ENDIF
  ELSE
    IF (OUTTYPE == 'R') THEN
      UNITSNAME = 'WATTS/(M^2 STER)'
    ELSE
      UNITSNAME = 'WATTS/(M^2)'
    ENDIF
  ENDIF
  IF (OUTTYPE == 'R')  OUTNAME = 'RADIANCE'
  IF (OUTTYPE == 'F')  OUTNAME = 'FLUX'


  OPEN (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')
  WRITE (2,'(A)') '! SHDOM Plane-parallel Radiative Transfer Output'
  WRITE (2,'(2(A,I3))') '!  NMU=',NMU, '  NPHI=',NPHI
  WRITE (2,'(2(A,I4))') '!  NLAYER(input)=',NLAY, '    NPTS=',NPTS
  WRITE (2,'(A,A32)')     '!  PROPERTY_FILE=', PROPFILE
  WRITE (2,'(A,A32,A,I2)') '!  CORRELATED_K-DIST_FILE=', CKDFILE, & 
                           '   NUM_G=', NG
  WRITE (2,'(A,A14)')    '!  SOURCE_TYPE=', SOURCENAME
  WRITE (2,'(A,A20)') '!  SURFACE_TYPE=', SFCNAME
  IF (SFCTYPE == 'L') THEN
    WRITE (2,'(A,F9.7)') '!  SURFACE_ALBEDO=',SFCPARMS(1)
  ELSE IF (SFCTYPE == 'R') THEN
    WRITE (2,'(3(A,F8.5))') '!  SURFACE_RHO0=',SFCPARMS(1), &
          '   SURFACE_K=',SFCPARMS(2), '   SURFACE_THETA=',SFCPARMS(3)
  ENDIF
  WRITE (2,'(A,F8.3,A,E12.5)') &
          '!  SURFACE_TEMP=',SFCTEMP, '  SKY_RAD=',SKYRAD

  IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
    WRITE (2,'(A,E13.6,A,F10.7)') &
         '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU
  ENDIF
  IF (WAVELEN > 0.0) THEN
    WRITE (2,'(A,A24,A,F10.2)') '!  UNITS=',UNITSNAME, &
                                '   WAVELENGTH(MICRON)=',WAVELEN
  ELSE
    WRITE (2,'(A,A18,A,2F10.3)') '!  UNITS=',UNITSNAME, &
                                '   WAVENUMBER_RANGE(CM^-1)=',WAVENO(1:2)
  ENDIF
  WRITE (2,'(A,E10.3)') '!  SPLITTING_ACCURACY=', SPLITACC
  WRITE (2,'(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
  WRITE (2,'(2(A,I4))') '!  MAXIMUM_ITERATIONS=', MAXITER, &
                        '   NUMBER_ITERATIONS=', ITER
  WRITE (2,'(A,A20)')   '!  OUTPUT_TYPE=', OUTNAME


  IF (OUTTYPE .EQ. 'R') THEN
    ! Radiance output
    WRITE (2,'(3(A,I3))') '!   NZOUT=',NZOUT, &
                          '  NMUOUT=',NMUOUT, '  NPHIOUT=',NPHIOUT
    WRITE (2,'(A)') '!    Z       MU     PHI     RADIANCE'
    DO I = 1, NZOUT
      DO J = 1, NMUOUT
        DO K = 1, NPHIOUT
          WRITE (2,'(1X,F7.3,1X,F8.5,1X,F6.2,1X,E13.6)') &
              ZOUT(I), MUOUT(J), PHIOUT(K), RADOUT(I,J,K)
        ENDDO
      ENDDO
    ENDDO

  ELSE IF (OUTTYPE .EQ. 'F') THEN
    ! Hemispheric flux output
    !   Format 1 is fluxes at every property grid point
    !   Format 2 is fluxes at every grid point
    WRITE (2,'(A)') '!    Z         F_UP        F_DOWN       F_ACTINIC'
    IF (FLUXPARM == 1) THEN
      DO L = 1, NLAY
        I = IXP(1,L)
        WRITE (2,'(1X,F7.3,3(2X,E12.5))') &
            HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
      ENDDO
      I = IXP(2,NLAY)
      WRITE (2,'(1X,F7.3,3(2X,E12.5))') &
          HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
    ELSE
      DO L = 1, NLAY
        DO I = IXP(1,L), IXP(2,L)-1
          WRITE (2,'(1X,F7.3,3(2X,E12.5))') &
            HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
        ENDDO
      ENDDO
      I = IXP(2,NLAY)
      WRITE (2,'(1X,F7.3,3(2X,E12.5))') &
          HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
    ENDIF

  ENDIF
  CLOSE (2)
END SUBROUTINE OUTPUT_RESULTS
 
