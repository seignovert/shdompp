PROGRAM DISORTSH
 ! Driver program that calls DISORT to run plane-parallel radiative transfer
 ! with an SHDOMPP property file input.

  IMPLICIT NONE
  ! DISORT array sizes
  INTEGER :: MAXCLY, MAXMOM, MAXULV, MAXUMU, MAXPHI
    ! MAXCLY   :  Max. number of computational layers
    ! MAXMOM   :  Max. number of phase function moment coefficients
    ! MAXULV   :  Max. number of user levels
    ! MAXUMU   :  Max. number of user polar angles
    ! MAXPHI   :  Max. number of user azimuth angles
  !    DISORT variables - see below
  LOGICAL :: LAMBER, PLANK, ONLYFL, PRNT(5), USRANG, USRTAU
  INTEGER :: IBCND, NLYR, NMOM, NUMU, NSTR, NPHI, NTAU
  REAL    :: FBEAM, UMU0, PHI0, FISOT, TEMIS, TTEMP, BTEMP, ALBEDO
  REAL    :: WVNMLO, WVNMHI, ACCUR
  REAL, ALLOCATABLE :: DTAUC(:), SSALB(:), PMOM(:,:), TEMPER(:)
  REAL, ALLOCATABLE :: UTAU(:), PHI(:), UMU(:)
  REAL, ALLOCATABLE :: RFLDIR(:), RFLDN(:), FLUP(:), DFDT(:), UAVG(:)
  REAL, ALLOCATABLE :: UU(:,:,:), ALBMED(:), TRNMED(:)
  CHARACTER(LEN=127) ::  HEADER


  ! SHDOMPP variables
  ! Number of output files and inputs parameter per file
  INTEGER, PARAMETER   :: MAXPAR=100  ! Max array size for user input arrays
  REAL                 :: ZOUT(MAXPAR), MUOUT(MAXPAR), PHIOUT(MAXPAR)

  INTEGER :: NMU, NUMOUT, MAXLEG, NZOUT, NMUOUT, NPHIOUT
  INTEGER :: NRAD, IRAD, I, J, K, N, L, NCASES, ICASE, NG, NZCKD, IG
  INTEGER, ALLOCATABLE :: NLEG(:), LOUT(:)
  LOGICAL :: KDIST
  REAL    :: SOLARFLUX, SOLFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO(2)
  REAL    :: ZDIF, MINDIF
  REAL,    ALLOCATABLE :: HEIGHTS(:), TAUP(:), ALBEDOP(:)
  REAL,    ALLOCATABLE :: DELG(:), ZCKD(:), KABS(:,:)
  REAL,    ALLOCATABLE :: FLUXUP(:), FLUXDN(:), MEANRAD(:), RADOUT(:,:,:)
  CHARACTER(LEN=1) :: SRCTYPE
  CHARACTER(LEN=80) :: PROPFILE, CKDFILE, FLUXFILE, RADFILE


  CALL USER_INPUT (PROPFILE, CKDFILE, KDIST, NMU, SRCTYPE, &
                   SOLARFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO, &
                   FLUXFILE, RADFILE, MAXPAR, &
                   NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT)

  NMU = MAX(2, 2*INT((NMU+1)/2) )
  NSTR = NMU

    ! Find out how many radiance directions we will be generating
  ALLOCATE (UMU(NMUOUT), PHI(NPHIOUT), RADOUT(NZOUT,NMUOUT,NPHIOUT))
  NPHI = NPHIOUT
  MAXPHI = NPHIOUT
  NUMU = NMUOUT
  MAXUMU = NMUOUT
  UMU(1:NMUOUT) = MUOUT(1:NMUOUT)
  PHI(1:NPHIOUT) = PHIOUT(1:NPHIOUT)


  CALL READ_PROPERTY_SIZE (PROPFILE, NLYR, MAXLEG)
  MAXLEG = MAX(MAXLEG,NMU)
  ALLOCATE (HEIGHTS(NLYR+1), TEMPER(NLYR+1), TAUP(NLYR), ALBEDOP(NLYR))
  ALLOCATE (DTAUC(NLYR), SSALB(NLYR), NLEG(NLYR))
  ALLOCATE (PMOM(0:MAXLEG,NLYR))
  MAXCLY = NLYR
  MAXULV = NLYR+1
  MAXMOM = MAXLEG
  NMOM = MAXLEG
  ALLOCATE (RFLDIR(MAXULV), RFLDN(MAXULV), FLUP(MAXULV), DFDT(MAXULV))
  ALLOCATE (UTAU(MAXULV), UAVG(MAXULV), UU(MAXUMU,MAXULV,MAXPHI))
  ALLOCATE (ALBMED(MAXUMU), TRNMED(MAXUMU))
  ALLOCATE (LOUT(MAXULV))
  ALLOCATE (FLUXUP(MAXULV), FLUXDN(MAXULV), MEANRAD(MAXULV))

  ! Read in the properties of the medium
  CALL READ_PROPERTIES (PROPFILE, NLYR, MAXLEG, &
                        HEIGHTS, TEMPER, TAUP, ALBEDOP, NLEG, PMOM)
  DO L = 0, MAXLEG   ! Convert Legendre coefficients to DISORT definition
    PMOM(L,:) = PMOM(L,:)/(2*L+1)
  ENDDO

  ! Have disort output radiances
  ONLYFL = .FALSE.

  ! Radiant quantities are to be returned at boundary of every 
  ! computational layer.
  USRTAU = .FALSE.

  ! Radiant quantities are to be returned at user polar angles. 
  USRANG = .TRUE.

  ! General top and bottom boundary conditions 
  IBCND = 0

  IF (SRCTYPE == 'S' .OR. SRCTYPE == 'B') THEN
    FISOT = SKYRAD  ! Intensity of top-boundary isotropic illumination.
    ! Intensity of incident parallel beam at top boundary (relative units)
    FBEAM = SOLARFLUX/SOLARMU
  ELSE
    TTEMP = SKYRAD  ! Temperature of top boundary (K).
    TEMIS = 1.0     ! Emissivity of top boundary.
    FBEAM = 0.0
  ENDIF

  ! Polar angle cosine of incident beam (positive).
  UMU0 = ABS(SOLARMU)
  ! Azimuth angle of incident beam (0 to 360 degrees)
  PHI0 = 0.0

  ! If true, isotropically reflecting (Lambertian) bottom boundary
  LAMBER = .TRUE.
  ! Bottom-boundary albedo
  ALBEDO = SFCALB

  IF (SRCTYPE == 'T' .OR. SRCTYPE == 'B') THEN
    ! Temperature of bottom boundary (K)  (bottom emissivity from ALBEDO)
    BTEMP = SFCTEMP
    PLANK = .TRUE.   ! Include thermal emission for thermal source
                     ! Wavenumbers (inv cm) of spectral interval of interest
    WVNMLO = WAVENO(1)
    WVNMHI = WAVENO(2)
  ELSE
    BTEMP = 0.0
    PLANK = .FALSE.
  ENDIF


  ! Find levels of the desired radiances for output (pick closest)
  DO I = 1, NZOUT
    MINDIF = 1.0E10
    DO L = 1, NLYR+1
      ZDIF = ABS(ZOUT(I)-HEIGHTS(L))
      IF (ZDIF < MINDIF) THEN
        MINDIF = ZDIF
        LOUT(I) = L
      ENDIF
    ENDDO
    ZOUT(I) = HEIGHTS(LOUT(I))
  ENDDO

  !              CONTROL FLAGS

  ! Convergence criterion for azimuthal (Fourier cosine) series.
  ACCUR = 0.0001

  ! Array of LOGICAL print flags causing the following prints:
  !         L        quantities printed
  !        --        ------------------
  !         1        input variables (except PMOM)
  !         2        fluxes
  !         3        intensities at user levels and angles
  !         4        planar transmissivity and planar albedo
  !                  as a function solar zenith angle ( IBCND = 1 )
  !         5        phase function moments PMOM for each layer
  !                  ( only if PRNT(1) = TRUE, and only for layers
  !                  with scattering )
  PRNT(1) = .FALSE.
  PRNT(2) = .FALSE.
  PRNT(3) = .FALSE.
  PRNT(4) = .FALSE.
  PRNT(5) = .FALSE.

  ! A 127 character header for prints, embedded in a DISORT banner;
  ! setting HEADER = '' will eliminate both the banner and the header.
  HEADER = ' '


  ! If doing a k-distribution then get the band info from the CKD file
  IF (KDIST) THEN
    CALL READ_CKD_SIZE (CKDFILE, WAVENO, NG, NZCKD)
    ALLOCATE (DELG(NG), ZCKD(NZCKD), KABS(NZCKD,NG))
    CALL READ_CKD (CKDFILE, WAVENO, NG, NZCKD, SOLFLUX, DELG, ZCKD, KABS)
    FBEAM = SOLFLUX*SOLARFLUX
  ELSE
    NG = 1
    NZCKD = 0
    ALLOCATE (DELG(NG), ZCKD(NZCKD), KABS(NZCKD,NG))
    DELG(1) = 1.0
  ENDIF


!  WRITE (*,*) 'Number of times to run DISORT routine for timing'
!  READ (*,*) NCASES
!    WRITE (*,*) NCASES
!  DO ICASE = 1, NCASES

  ! Loop over the k-distribution g's from low to high absorption
  DO IG = 1, NG
    ! Add in the layer gaseous absorption if doing a k-distribution
    CALL ADD_IN_KDIST (NLYR, HEIGHTS, TAUP, ALBEDOP, &
                       NZCKD, ZCKD, KABS(:,IG),  DTAUC, SSALB)

    CALL  DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, &
                WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, &
                NUMU, UMU, NPHI, PHI, IBCND,  &
                FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, &
                BTEMP, TTEMP, TEMIS, PLANK, ONLYFL, &
                ACCUR, PRNT, HEADER, &
                MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM, &
                RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

    IF (KDIST) THEN
      CALL SUM_KDIST_OUTPUT (IG, DELG(IG), NLYR, FLUP, FLUXUP, &
                             RFLDIR+RFLDN, FLUXDN, UAVG, MEANRAD, &
                             NZOUT, NMUOUT, NPHIOUT, LOUT, UU, RADOUT)
    ELSE
      ! Put the upwelling and downwelling fluxes in the output array
      FLUXUP(:) = FLUP(:)
      FLUXDN(:) = RFLDIR(:) + RFLDN(:)
      MEANRAD(:) = UAVG(:)
      ! Put the radiances at desired levels in the output array
      DO J = 1, NMUOUT
        DO K = 1, NPHIOUT
          DO I = 1, NZOUT
            RADOUT(I,J,K) = UU(J,LOUT(I),K)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  ENDDO  ! End of k-distribution loop

!  ENDDO


  CALL OUTPUT_RESULTS (NLYR, NSTR, NG, HEIGHTS, PROPFILE, CKDFILE, SRCTYPE, &
                       SOLARFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO, &
                       FLUXUP, FLUXDN, MEANRAD,  NZOUT, ZOUT, &
                       NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
                       'F', FLUXFILE)

  IF (NZOUT*NMUOUT*NPHIOUT > 0) THEN
    CALL OUTPUT_RESULTS (NLYR, NSTR, NG, HEIGHTS, PROPFILE, CKDFILE, SRCTYPE, &
                       SOLARFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO, &
                       FLUXUP, FLUXDN, MEANRAD,  NZOUT, ZOUT, &
                       NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
                       'R', RADFILE)
  ENDIF

  DEALLOCATE (RFLDIR, RFLDN, FLUP, DFDT, UAVG, UTAU, UU, ALBMED, TRNMED)
  DEALLOCATE (DTAUC, SSALB, NLEG, PMOM, TEMPER)
  DEALLOCATE (UMU, PHI, RADOUT, FLUXUP, FLUXDN, MEANRAD)
  DEALLOCATE (HEIGHTS, TAUP, ALBEDOP, DELG, ZCKD, KABS)
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



SUBROUTINE SUM_KDIST_OUTPUT (IG, WT, NLAY, FLUP, FLUXUP, &
                             FLDN, FLUXDN, UAVG, MEANRAD, &
                             NZOUT, NMUOUT, NPHIOUT, LOUT, UU, RADOUT)
 ! Sums the output arrays over the k-distribution.  
  INTEGER, INTENT(IN) :: IG, NLAY, NZOUT, NMUOUT, NPHIOUT, LOUT(NLAY+1)
  REAL,    INTENT(IN) :: WT, UU(NMUOUT,NLAY+1,NPHIOUT)
  REAL,    INTENT(IN) :: FLUP(NLAY+1), FLDN(NLAY+1), UAVG(NLAY+1)
  REAL,    INTENT(OUT) :: FLUXUP(NLAY+1), FLUXDN(NLAY+1), MEANRAD(NLAY+1)
  REAL,    INTENT(OUT) :: RADOUT(NZOUT,NMUOUT,NPHIOUT)
  INTEGER :: I, J, K

  IF (IG == 1) THEN
    FLUXUP(:) = 0.0
    FLUXDN(:) = 0.0
    MEANRAD(:) = 0.0
    RADOUT(:,:,:) = 0.0
  ENDIF
  
  FLUXUP(:) = FLUXUP(:) + WT*FLUP(:)
  FLUXDN(:) = FLUXDN(:) + WT*FLDN(:)
  MEANRAD(:) = MEANRAD(:) + WT*UAVG(:)

  ! Add in the radiances at desired levels in the output array
  DO J = 1, NMUOUT
    DO K = 1, NPHIOUT
      DO I = 1, NZOUT
        RADOUT(I,J,K) = RADOUT(I,J,K) + WT*UU(J,LOUT(I),K)
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE SUM_KDIST_OUTPUT





SUBROUTINE USER_INPUT (PROPFILE, CKDFILE, KDIST, NMU, &
                       SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                       SFCTEMP, SFCALB, WAVENO, &
                       FLUXFILE, RADFILE, MAXPAR, &
                       NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT)
!       Obtains the input parameters for the program by writing prompts
!     to and reading responses from the terminal.  Echos the inputs to
!     make a useful log file when running non-interactively.  See the overall
!     program documentation for the list of input parameters.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: MAXPAR
  INTEGER, INTENT(OUT) :: NMU, NZOUT, NMUOUT, NPHIOUT
  LOGICAL, INTENT(OUT) :: KDIST
  REAL,    INTENT(OUT) :: SOLARFLUX, SOLARMU
  REAL,    INTENT(OUT) :: SFCTEMP, SFCALB, SKYRAD, WAVENO(2)
  REAL,    INTENT(OUT) :: ZOUT(MAXPAR), MUOUT(MAXPAR), PHIOUT(MAXPAR)
  CHARACTER(LEN=1),  INTENT(OUT) :: SRCTYPE
  CHARACTER(LEN=80), INTENT(OUT) :: PROPFILE, CKDFILE, FLUXFILE, RADFILE 

  WRITE (*,*)
  WRITE (*,'(1X,A)') 'Atmospheric optical properties file name'
  READ (*,'(A)') PROPFILE
  WRITE (*,*) PROPFILE

  WRITE (*,'(1X,A)') 'Correlated k-distribution file name (or NONE)'
  READ (*,'(A)') CKDFILE
  WRITE (*,*) CKDFILE
  KDIST = CKDFILE(1:2) .NE. 'NO'

  WRITE (*,'(1X,A)') 'Number of discrete ordinates in -1 < mu < 1'
  READ (*,*) NMU
    WRITE (*,*) NMU

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
    SOLARMU = ABS(SOLARMU)
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

  WRITE (*,'(1X,A)') 'Lambertian surface albedo'
  READ (*,*) SFCALB
    WRITE (*,*) SFCALB

  IF (KDIST .OR. SRCTYPE .EQ. 'T' .OR. SRCTYPE .EQ. 'B') THEN
    WRITE (*,'(1X,A)') 'Wavenumber range (cm^-1)'
    READ (*,*) WAVENO(1:2)
      WRITE (*,*) WAVENO(1:2)
  ENDIF

  WRITE (*,'(1X,A,A)') 'Flux output file name'
  READ (*,'(A)') FLUXFILE
    WRITE (*,'(A)') FLUXFILE

  WRITE (*,'(1X,A,A)') 'Radiance output file name'
  READ (*,'(A)') RADFILE
    WRITE (*,'(A)') RADFILE

  WRITE (*,*) 'Number of optical depth levels for output radiance'
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




SUBROUTINE OUTPUT_RESULTS (NLAY, NMU, NG, HEIGHTS, PROPFILE, CKDFILE, SRCTYPE, &
                         SOLARFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO, & 
                         FLUXUP, FLUXDN, MEANRAD,  NZOUT, ZOUT, &
                         NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
                         OUTTYPE, OUTFILE)
!       Writes the desired type of output file from the output fields
!     with a header giving the input parameters.
!     There are two types (OUTTYPE) of output: 'R' - for radiance,  
!     'F' - for hemispheric/actinic flux.
!     See the overall program documentation for output parameters.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY, NMU, NG
  INTEGER, INTENT(IN) :: NZOUT, NMUOUT, NPHIOUT
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, SFCTEMP, SFCALB, SKYRAD
  REAL,    INTENT(IN) :: WAVENO(2)
  REAL,    INTENT(IN) :: HEIGHTS(NLAY+1)
  REAL,    INTENT(IN) :: FLUXUP(NLAY+1), FLUXDN(NLAY+1), MEANRAD(NLAY+1)
  REAL,    INTENT(IN) :: ZOUT(NZOUT), MUOUT(NMUOUT), PHIOUT(NPHIOUT)
  REAL,    INTENT(IN) :: RADOUT(NZOUT,NMUOUT,NPHIOUT)
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE, OUTTYPE
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

  SFCNAME = 'LAMBERTIAN'

  IF (OUTTYPE == 'R') THEN
    UNITSNAME = 'WATTS/(M^2 STER)'
  ELSE
    UNITSNAME = 'WATTS/(M^2)'
  ENDIF
  IF (OUTTYPE == 'R')  OUTNAME = 'RADIANCE'
  IF (OUTTYPE == 'F')  OUTNAME = 'FLUX'


  OPEN (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')
  WRITE (2,'(A)') '! DISORT Plane-parallel Radiative Transfer Output'
  WRITE (2,'(2(A,I3))') '!  N_STREAMS=',NMU, '    N_LAYERS=',NLAY
  WRITE (2,'(A,A32)')   '!  PROPERTY_FILE=', PROPFILE
  WRITE (2,'(A,A32,A,I2)') '!  CORRELATED_K-DIST_FILE=', CKDFILE, & 
                           '   NUM_G=', NG
  WRITE (2,'(A,A14)')   '!  SOURCE_TYPE=', SOURCENAME
  WRITE (2,'(A,A20)')   '!  SURFACE_TYPE=', SFCNAME
  WRITE (2,'(A,F9.7)')  '!  SURFACE_ALBEDO=',SFCALB
  WRITE (2,'(A,F8.3,A,E12.5)') '!  SURFACE_TEMP=',SFCTEMP, '  SKY_RAD=',SKYRAD
  IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
    WRITE (2,'(A,E13.6,A,F10.7)') &
         '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU
  ENDIF
  WRITE (2,'(A,A18,A,2F10.3)') '!  UNITS=',UNITSNAME, &
                                '   WAVENUMBER_RANGE(CM^-1)=',WAVENO(1:2)
  WRITE (2,'(A,A20)')   '!  OUTPUT_TYPE=', OUTNAME


  IF (OUTTYPE .EQ. 'R') THEN
    ! Radiance output
    WRITE (2,'(3(A,I3))') '!   NZOUT=',NZOUT, &
                          '  NMUOUT=',NMUOUT, '  NPHIOUT=',NPHIOUT
    WRITE (2,'(A)') '!    Z       MU     PHI     RADIANCE'
    DO I = 1, NZOUT
      DO J = 1, NMUOUT
        DO K = 1, NPHIOUT
          WRITE (2,'(1X,F7.3,1X,F8.5,1X,F6.2,1X,E12.5)') &
              ZOUT(I), MUOUT(J), PHIOUT(K), RADOUT(I,J,K)
        ENDDO
      ENDDO
    ENDDO

  ELSE IF (OUTTYPE .EQ. 'F') THEN
    ! Hemispheric flux output
    WRITE (2,'(A)') '!    Z         F_UP        F_DOWN       F_ACTINIC'
    DO I = 1, NLAY+1
      WRITE (2,'(1X,F7.3,3(2X,E12.5))') HEIGHTS(I), &
            FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
    ENDDO

  ENDIF
  CLOSE (2)
END SUBROUTINE OUTPUT_RESULTS
 

