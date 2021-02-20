      !#########################################################################
      !>
      !> TEST GTOR AND RTOG FUNCTIONS
      !>
      !> \remark TEST DIRECT AND INVERSE GEODETIC PROGLEM SOLVING
      !> \n      RESULT IS WRITTEN TO 'TEST_GTOR.TXT' FILE
      !>
      PROGRAM TEST_GTOR_RTOG
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! CONSTANTS
      REAL :: C = 300000000.0          ! SPEED OF LIGHT IN VACUUM,[m/s]
      REAL :: RZ = 6371000.0           ! RADIUS OF EARTH,[m]
      REAL(4) :: DGTORD = ASIN(1.)/90. ! DEGREES TO RADIANS = ASIN(1.D0)/90.D0 = 0.017...
      REAL(4) :: RDTODG = 90./ASIN(1.) ! RADIANS TO DEGREES = 90.D0/ASIN(1.D0) = 57.296...
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      REAL TETTA1,PHI1,R,AZ,AZEND,TETTA3,PHI3,AZEND3
      REAL TETTA2(8),PHI2(8)
      CHARACTER(250) SHARPSLINE   ! LINE OF '#' CHARACTERS
      CHARACTER(100) ELIPSOIDNAME ! NAME OF ELLIPSOID
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! COMMON BLOCKS
      COMMON /CONSTS/ C,RZ,DGTORD,RDTODG                 ! CONSTANTS BLOCK
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      OPEN(900,FILE='TEST_GTOR.TXT',STATUS='UNKNOWN')
      DO I=1,250
          SHARPSLINE(I:I)='#'
      ENDDO
      TETTA1=0.
      PHI1=0.
      TETTA2=(/30.,30.,0.,-30.,-30.,-30.,0.,30./) ! LATITUDE,DEG
      PHI2=(/0.,30.,30.,30.,0.,-30.,-30.,-30./)   ! LONGITUDE,DEG
      R=0.
      AZ=0.
      AZEND=0.
      ITESTSCOUNT=0
      ITESTSPASSED=0
      TOL=0.25 ! TOLERANCE OF LAT/LON, DEG
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! TEST 5 ELLIPSOIDS:
      ! 1 : SPHERE R=6371 [km] (AVERAGE RADIUS)
      ! 2 : SPHERE R=6378 [km] (EQUATORIAL RADIUS)
      ! 3 : ELLIPSOID WGS-84 (6378137.0 [m], 1 / 298.257223563)
      ! 4 : ELLIPSOID PZ-90  (6378136.0 [m], 1 / 298.257839303)
      ! 5 : ELLIPSOID GRS-80 (6378137.0 [m], 1 / 298.257222101)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO 12 NELIPSOID=1,5
      SELECTCASE(NELIPSOID)
      CASE(1)
         ELIPSOIDNAME='SPHERE R=6371 [km] (AVERAGE RADIUS)'
      CASE(2)
         ELIPSOIDNAME='SPHERE R=6378 [km] (EQUATORIAL RADIUS)'
      CASE(3)
         ELIPSOIDNAME=
     *   'ELLIPSOID WGS-84 (6378137.0 [m], 1 / 298.257223563)'
      CASE(4)
         ELIPSOIDNAME=
     *   'ELLIPSOID PZ-90  (6378136.0 [m], 1 / 298.257839303)'
      CASE(5)
         ELIPSOIDNAME=
     *   'ELLIPSOID GRS-80 (6378137.0 [m], 1 / 298.257222101)'
      ENDSELECT
      NUM=1 ! NUMBER OF TEST
      DO 11 I=1,8
      DO 10 J=1,2
      ITESTSCOUNT=ITESTSCOUNT+1
      WRITE(900,'(/1X, A250)') SHARPSLINE
      WRITE(900,'(/1X, A100)') ELIPSOIDNAME
      WRITE(900,1) NUM
1     FORMAT(1X, 'Test #', I2)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL GTOR(NELIPSOID,TETTA1,PHI1,TETTA2(I)*J,PHI2(I)*J,R,AZ,AZEND)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      WRITE(900,2) TETTA1,PHI1,TETTA2(I)*J,PHI2(I)*J,R,AZ,AZEND
2     FORMAT(
     * /1X, 'StartPoint:',
     * /1X, ' TETTA1 =',F11.4,
     * /1X, ' PHI1   =',F11.4,
     * //1X,'EndPoint:',
     * /1X, ' TETTA2 =',F11.4,
     * /1X, ' PHI2   =',F11.4,
     * //1X,'SOLVING INVERSE PROBLEM (G TO R)',
     * //1X,' R      =',F11.4,
     * /1X, ' AZ     =',F11.4,
     * /1X, ' AZEND  =',F11.4)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL RTOG(NELIPSOID,TETTA1,PHI1,R,AZ,TETTA3,PHI3,AZEND3)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      WRITE(900,3) TETTA3,PHI3,AZEND3 !
3     FORMAT(
     * //1X,'SOLVING DIRECT PROBLEM (R TO G)',
     * //1X,' TETTA3 =',F11.4,
     * /1X, ' PHI3   =',F11.4,
     * /1X, ' AZEND3 =',F11.4)

      IF(ABS(TETTA2(I)*J-TETTA3).LE.TOL.AND.
     *   ABS(PHI2(I)*J-PHI3).LE.TOL) ITESTSPASSED=ITESTSPASSED+1

      NUM=NUM+1! INCREMENT NUMBER OF TEST
10    CONTINUE ! CYCLE FOR COORDINATES
11    CONTINUE ! CYCLE FOR COORDINATES
12    CONTINUE ! CYCLE FOR ELLIPSOIDS

      IF(ITESTSPASSED.EQ.ITESTSCOUNT)THEN
         WRITE(900,'(/1X, A)') 'TESTS PASSED SUCCESSFULLY'
         PRINT *,'TESTS PASSED SUCCESSFULLY'
      ELSE
         WRITE(900,'(/1X, A)') 'TESTS FAILED!'
         PRINT *,'TESTS FAILED!'
      ENDIF
      PRINT *,'PRESS ANY KEY TO EXIT...'
      READ *
      CLOSE(900)
      END
