PROGRAM PPMIEPRP
 ! Makes an SHDOMPP property file by combining optical properties of
 ! various components from Mie scattering tables.  The Mie table files
 ! are output by cloudprp.f or make_mie_table.f90.  Each has the optical 
 ! properties for one component as a function of particle effective radius 
 ! for an equivalent liquid water content of unity.  For non water or ice 
 ! particles the material density is assumed to be 1 for definition of the
 ! equivalent mass content (called LWC).  The liquid water path and
 ! effective radius are input for each component and each layer from
 ! stdin (if a component is absent from a layer simply set LWP=0).
 ! In addition to particle scattering with optical properties in
 ! Mie tables, molecular Rayleigh scattering, and molecular absorption
 ! are allowed.  Rayleigh scattering is handled by inputting a wavelength
 ! dependent Rayleigh coefficient, and the height and temperature profile.
 ! The density and hence Rayleigh optical depth are calculate from the 
 ! temperature and pressure difference across the layer (from the 
 ! hypsometric equation).  The molecular absorption for each layer is 
 ! assumed to be computed by some other program (e.g MODTRAN) and
 ! simply input from stdin.
 !
 !   Written by Frank Evans, University of Colorado,   October 2002
 !
 !    compile with  pgf90 -fast -o ppmieprp ppmieprp.f90

  IMPLICIT NONE
  INTEGER, PARAMETER :: MAXCOMP=5, MAXNRE=50, MAXLAY=100, MAXLEG=5000
  INTEGER :: NCOMP, NLAY, I
  INTEGER :: NRETAB(MAXCOMP), NLEGTAB(MAXNRE,MAXCOMP)
  REAL    :: SRETAB(MAXCOMP), ERETAB(MAXCOMP)
  REAL    :: EXTINCTAB(MAXNRE,MAXCOMP), SSALBTAB(MAXNRE,MAXCOMP)
  REAL    :: LEGENTAB(0:MAXLEG,MAXNRE,MAXCOMP)
  REAL    :: HEIGHTS(MAXLAY+1), TEMPS(MAXLAY+1)
  REAL    :: LWP(MAXCOMP,MAXLAY), REFF(MAXCOMP,MAXLAY)
  REAL    :: MOLABS(MAXLAY), RAYTAU(MAXLAY), RAYLCOEF
  CHARACTER(LEN=80) :: MIEFILES(MAXCOMP), PROPFILE


  CALL USER_INPUT (MAXCOMP, MAXLAY, NCOMP, NLAY, MIEFILES, &
                   HEIGHTS, TEMPS, LWP, REFF, MOLABS, RAYLCOEF, PROPFILE)

  ! Read in the Mie scattering tables
  DO I = 1, NCOMP
    CALL READ_MIE_TABLE (MIEFILES(I), MAXNRE, NRETAB(I), SRETAB(I), ERETAB(I), &
                         EXTINCTAB(1,I), SSALBTAB(1,I), MAXLEG, &
                         NLEGTAB(1,I), LEGENTAB(0,1,I))
  ENDDO

  ! Compute the molecular Rayleigh scattering from density
  CALL RAYLEIGH_TAU (NLAY, HEIGHTS, TEMPS, RAYLCOEF, RAYTAU)


  ! Combine the component optical properties and output the property file
  CALL OUTPUT_PROPERTIES (MAXCOMP, MAXNRE, MAXLEG, NCOMP, NLAY, &
                          NRETAB, SRETAB, ERETAB, &
                          EXTINCTAB, SSALBTAB, NLEGTAB, LEGENTAB, &
                          HEIGHTS, TEMPS, LWP, REFF, RAYTAU, MOLABS, &
                          PROPFILE)

END





SUBROUTINE USER_INPUT (MAXCOMP, MAXLAY, NCOMP, NLAY, MIEFILES, &
                       HEIGHTS, TEMPS, LWP, REFF, MOLABS, RAYLCOEF, PROPFILE)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MAXCOMP, MAXLAY
  INTEGER, INTENT(OUT) :: NCOMP, NLAY
  REAL,    INTENT(OUT) :: HEIGHTS(MAXLAY+1), TEMPS(MAXLAY+1)
  REAL,    INTENT(OUT) :: LWP(MAXCOMP,MAXLAY), REFF(MAXCOMP,MAXLAY)
  REAL,    INTENT(OUT) :: MOLABS(MAXLAY), RAYLCOEF
  CHARACTER(LEN=*), INTENT(OUT) :: MIEFILES(MAXCOMP), PROPFILE
  INTEGER :: I, L

  WRITE (*,*) 'Number of optical components (Mie tables)'
  READ (*,*) NCOMP
  IF (NCOMP > MAXCOMP) STOP 'USER_INPUT: MAXCOMP exceeded'
    WRITE (*,*) NCOMP

  DO I = 1, NCOMP
    WRITE (*,'(A,I2)') 'Mie table file name',I
    READ (*,'(A)') MIEFILES(I)
      WRITE (*,*) MIEFILES(I)
  ENDDO

  WRITE (*,*) 'Number of layers'
  READ (*,*) NLAY
  IF (NLAY > MAXLAY) STOP 'USER_INPUT: MAXLAY exceeded'
    WRITE (*,*) NLAY

  WRITE (*,*) 'Height of layer boundaries from the top down (km)'
  READ (*,*) HEIGHTS(1:NLAY+1)

  WRITE (*,*) 'Layer boundary temperatures from the top down (K)'
  READ (*,*) TEMPS(1:NLAY+1)

  DO L = 1, NLAY
    WRITE (*,'(A,I2)') 'Equivalent liquid water paths (g/m^2) of each component for layer',L
    READ (*,*) LWP(1:NCOMP,L)
      WRITE (*,*) LWP(1:NCOMP,L)
    WRITE (*,'(A,I2)') 'Effective radii (micron) of each component for layer',L
    READ (*,*) REFF(1:NCOMP,L)
      WRITE (*,*) REFF(1:NCOMP,L)
  ENDDO

  WRITE (*,*) 'Molecular absorption optical depth for each layer'
  READ (*,*) MOLABS(1:NLAY)
    WRITE (*,*) MOLABS(1:NLAY)

  WRITE (*,*) 'Molecular scattering coefficient (K/(km mb))'
  READ (*,*) RAYLCOEF
    WRITE (*,*) RAYLCOEF

  WRITE (*,*) 'Output SHDOMPP property file name'
  READ (*,'(A)') PROPFILE
    WRITE (*,*) PROPFILE

END SUBROUTINE USER_INPUT



SUBROUTINE OUTPUT_PROPERTIES (MAXCOMP, MAXNRE, MAXLEG, NCOMP, NLAY, &
                              NRETAB, SRETAB, ERETAB, &
                              EXTINCTAB, SSALBTAB, NLEGTAB, LEGENTAB, &
                              HEIGHTS, TEMPS, LWP, REFF, RAYTAU, MOLABS, &
                              PROPFILE)
 ! Combines the optical properties of the components in each layer
 ! according to the equivalent liquid water path and effective radius
 ! of the components and outputs the SHDOMPP property file.
 ! Interpolates the Mie scattering property tables in effective radius
 ! for the non zero components in each layer, and adds in the molecular 
 ! Rayleigh scattering optical depth (RAYTAU) and the molecular 
 ! absorption (MOLABS),
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MAXCOMP, MAXNRE, MAXLEG, NCOMP, NLAY
  INTEGER, INTENT(IN) :: NRETAB(NCOMP), NLEGTAB(MAXNRE,NCOMP)
  REAL,    INTENT(IN) :: SRETAB(NCOMP), ERETAB(NCOMP)
  REAL,    INTENT(IN) :: EXTINCTAB(MAXNRE,NCOMP), SSALBTAB(MAXNRE,NCOMP)
  REAL,    INTENT(IN) :: LEGENTAB(0:MAXLEG,MAXNRE,NCOMP)
  REAL,    INTENT(IN) :: HEIGHTS(NLAY+1), TEMPS(NLAY+1)
  REAL,    INTENT(IN) :: LWP(MAXCOMP,NLAY), REFF(MAXCOMP,NLAY)
  REAL,    INTENT(IN) :: RAYTAU(NLAY), MOLABS(NLAY)
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE
  INTEGER :: I, J, L, NLEG
  REAL    :: TAU_TOT, TAU_SCA, LEGEN(0:MAXLEG), RJ, F
  REAL    :: EXT, TAU, SSALB

  ! Open the property file and output the header
  OPEN (UNIT=2, FILE=PROPFILE, STATUS='UNKNOWN')
  WRITE (2,'(I4)') NLAY
  WRITE (2,'(1X,F7.3,1X,F6.2)') HEIGHTS(1), TEMPS(1)

  DO L = 1, NLAY
    TAU_TOT = MOLABS(L) + RAYTAU(L)
    TAU_SCA = RAYTAU(L)
    LEGEN(0:MAXLEG) = 0.0
    LEGEN(0) = RAYTAU(L)
    LEGEN(2) = 0.5*RAYTAU(L)
    NLEG = 0
    IF (TAU_SCA > 0) NLEG = 2
    DO I = 1, NCOMP
      IF (LWP(I,L) > 0.0) THEN
        RJ = 1 + (REFF(I,L)-SRETAB(I))*(NRETAB(I)-1)/(ERETAB(I)-SRETAB(I))
        IF (RJ > NRETAB(I)) THEN
          WRITE (*,'(A,2I4,F9.5,F8.2)') 'Warning: effective radius beyond table :', &
             L, I, LWP(I,L), REFF(I,L)
        ENDIF
        RJ = MAX(1.0,MIN(FLOAT(NRETAB(I)),RJ))
        J = INT(RJ)
        F = RJ - J
        EXT = (1-F)*EXTINCTAB(J,I) + F*EXTINCTAB(J+1,I)
        TAU = LWP(I,L)*EXT/1000.
        SSALB = (1-F)*SSALBTAB(J,I) + F*SSALBTAB(J+1,I)
        TAU_TOT = TAU_TOT + TAU
        TAU_SCA = TAU_SCA + TAU*SSALB
        NLEG = MAX(NLEG, MAX(NLEGTAB(J,I),NLEGTAB(J+1,I)) )
        LEGEN(0:NLEG) = LEGEN(0:NLEG) + TAU*SSALB  &
                  *( (1-F)*LEGENTAB(0:NLEG,J,I) + F*LEGENTAB(0:NLEG,J+1,I) )
      ENDIF
    ENDDO
    IF (TAU_TOT > 0.0) THEN
      SSALB = TAU_SCA/TAU_TOT
    ELSE
      SSALB = 0.0
    ENDIF
    IF (TAU_SCA > 0.0) THEN
      LEGEN(0:NLEG) = LEGEN(0:NLEG)/TAU_SCA
    ELSE
      LEGEN(0:NLEG) = 0.0
    ENDIF
    WRITE (2,'(1X,F7.3,1X,F6.2,1X,F10.5,1X,F8.6,1X,I4,10000F10.5)') &
                 HEIGHTS(L+1), TEMPS(L+1), TAU_TOT, SSALB, NLEG, LEGEN(1:NLEG)
  ENDDO
  CLOSE (2)
END SUBROUTINE OUTPUT_PROPERTIES



SUBROUTINE RAYLEIGH_TAU (NLAY, HEIGHTS, TEMPS, RAYLCOEF, RAYTAU)
 ! Computes the molecular Rayleigh layer optical depth profile RAYTAU
 ! from the temperature profile TEMPS [K] at HEIGHTS [km].  Assumes
 ! a linear lapse rate between levels to compute the pressure at
 ! each level.  The Rayleigh extinction is proportional to air
 ! density, with the coefficient RAYLCOEF in [K/(mb km)].
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLAY
  REAL,    INTENT(IN) :: HEIGHTS(NLAY+1), TEMPS(NLAY+1), RAYLCOEF
  REAL,    INTENT(OUT) :: RAYTAU(NLAY)
  INTEGER :: I
  REAL    :: PRES1, PRES2, LAPSE, TS, DZ

  ! Find surface pressure by integrating hydrostatic relation
  !   for a dry atmosphere up to surface height.
  PRES1 = 1013.
  TS = TEMPS(NLAY+1)
  LAPSE = 6.5*0.001
  PRES1 = PRES1*(TS/(TS+LAPSE*HEIGHTS(NLAY+1)*1000.))**(9.8/(287.*LAPSE))

  ! Use layer mean temperature to compute fractional pressure change.
  DO I = NLAY, 1, -1
    DZ = 1000.*(HEIGHTS(I)-HEIGHTS(I+1))
    LAPSE = (TEMPS(I+1)-TEMPS(I))/DZ
    IF (ABS(LAPSE) > 0.00001) THEN
      PRES2 = PRES1*(TEMPS(I)/TEMPS(I+1))**(9.8/(287.*LAPSE))
    ELSE
      PRES2 = PRES1*EXP(-9.8*DZ/(287.*TEMPS(I+1)))
    ENDIF
    ! RAYTAU(I) = RAYLCOEF*(0.5*PRES1/TEMPS(I+1)+0.5*PRES2/TEMPS(I+1))*DZ/1000.
    RAYTAU(I) = RAYLCOEF*(287./9.8)*(PRES1-PRES2)/1000.
    ! print '(1X,F5.2,1X,F6.2,1X,F6.1,1X,F6.2,1X,E10.3)', &
    !     HEIGHTS(I),TEMPS(I),PRES2,PRES1-PRES2,RAYTAU(I)
    PRES1 = PRES2
  ENDDO  
END SUBROUTINE RAYLEIGH_TAU



SUBROUTINE READ_MIE_TABLE (MIEFILE, MAXNRE, NRETAB, SRETAB, ERETAB, &
                           EXTINCT, SSALB, MAXLEG, NLEG, LEGEN)
 ! Reads a table of Mie scattering properties as a function of 
 ! effective radius.  
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: MIEFILE
  INTEGER, INTENT(IN) :: MAXNRE, MAXLEG
  INTEGER, INTENT(OUT) :: NRETAB, NLEG(MAXNRE)
  REAL,    INTENT(OUT) :: SRETAB, ERETAB
  REAL,    INTENT(OUT) :: EXTINCT(MAXNRE), SSALB(MAXNRE), LEGEN(0:MAXLEG,MAXNRE)
  INTEGER :: I, L
  REAL    :: REFF

  OPEN (UNIT=3, FILE=MIEFILE, STATUS='OLD')
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*) NRETAB, SRETAB, ERETAB
  IF (NRETAB .GT. MAXNRE) STOP 'READ_MIE_TABLE: MAXNRE exceeded'
  DO I = 1, NRETAB
    READ (3,*) REFF, EXTINCT(I), SSALB(I), NLEG(I)
    IF (NLEG(I) > MAXLEG) STOP 'READ_MIE_TABLE: MAXLEG exceeded'
    READ (3,*) (LEGEN(L,I), L=0,NLEG(I)) 
    LEGEN(NLEG(I)+1:MAXLEG,I) = 0.0
  ENDDO
  CLOSE (3)
END SUBROUTINE READ_MIE_TABLE



