      !#########################################################################
      !>
      !> \brief CONVERTS ANGLE TO [0..360) DEGREES OR [0..2*PI) RADIANS
      !>
      !> \param ANGLE   ANGLE TO BE CONVERTED
      !> \param P       D - ANGLE IN [DEGREES], R - ANGLE IN [RADIANS]
      !> \return RETURN UNITS ARE THE SAME AS INPUT DUE TO 'P' PARAMETER
      !>
      REAL FUNCTION TO360(ANGLE,P)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      REAL ANGLE     ! ANGLE TO BE CONVERTED TO [0..360) DEGREES OR [0..2*PI) RADIANS
      CHARACTER(1) P ! D - ANGLE IN [DEGREES], R - ANGLE IN [RADIANS]
      !-------------------------------------------------------------------------
      REAL PI
      REAL K
      K=ANGLE
      IP=0
      IF(P.EQ.'D')IP=1
      IF(P.EQ.'R')IP=2

      SELECT CASE (IP)
         CASE (1)!'D'
            DO WHILE (K.GE.360.)
                K=K-360.
            ENDDO
            DO WHILE (K.LT.0.)
                K=K+360.
            ENDDO
         CASE (2)!'R'
            PI=ACOS(-1.)
            DO WHILE (K.GE.(2.*PI))
                K=REAL(K-(2.*PI))
            ENDDO
            DO WHILE (K.LT.0.)
                K=REAL(K+(2.*PI))
            ENDDO
      END SELECT
      TO360=K
      RETURN
      END

      !#########################################################################
      !>
      !> \brief CONVERTS ANGLE TO [0..360) DEGREES OR [0..2*PI) RADIANS
      !>
      !> \param ANGLE   ANGLE TO BE CONVERTED
      !> \param P       D - ANGLE IN [DEGREES], R - ANGLE IN [RADIANS]
      !> \return RETURN UNITS ARE THE SAME AS INPUT DUE TO 'P' PARAMETER
      !>
      REAL(8) FUNCTION TO360_DOUBLE(ANGLE,P)
      !DOUBLEPRECISION FUNCTION TO360_DOUBLE(ANGLE,P)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      REAL(8) ANGLE  ! ANGLE TO BE CONVERTED TO [0..360) DEGREES OR [0..2*PI) RADIANS
      CHARACTER(1) P ! D - ANGLE IN [DEGREES], R - ANGLE IN [RADIANS]
      !-------------------------------------------------------------------------
      REAL(8) PI
      REAL(8) K
      K=ANGLE
      IP=0
      IF(P.EQ.'D')IP=1
      IF(P.EQ.'R')IP=2
      SELECT CASE (IP)
         CASE (1) !'D'
            DO WHILE (K.GE.360.D0)
                K=K-360.D0
            ENDDO
            DO WHILE (K.LT.0.D0)
                K=K+360.D0
            ENDDO
         CASE (2) !'R'
            PI=ACOS(-1.D0)
            DO WHILE (K.GE.(2.D0*PI))
                K=(K-(2.D0*PI))
            ENDDO
            DO WHILE (K.LT.0.D0)
                K=(K+(2.D0*PI))
            ENDDO
      END SELECT
      TO360_DOUBLE=K
      RETURN
      END

      !#########################################################################
      !>
      !> \brief CONVERTS GEOGRAPHIC COORDINATES TO RADIOLOCATION COORDINATES (INVERSE TASK)
      !>
      !> \param N   NUMBER OF ELLIPSOID:
      !> \n         0 : SPHERE WITH RADIUS RZ FROM COMMON BLOCK
      !> \n         1 : SPHERE R=6371 [km] (AVERAGE RADIUS)
      !> \n         2 : SPHERE R=6378 [km] (EQUATORIAL RADIUS)
      !> \n         3 : ELLIPSOID WGS-84 (6378137.0 [m], 1 / 298.257223563)
      !> \n         4 : ELLIPSOID PZ-90  (6378136.0 [m], 1 / 298.257839303)
      !> \n         5 : ELLIPSOID GRS-80 (6378137.0 [m], 1 / 298.257222101)
      !> \param LAT_START LATITUDE OF START POINT,[deg]
      !> \param LON_START LONGITUDE OF START POINT,[deg]
      !> \param LAT_END   LATITUDE OF END POINT,[deg]
      !> \param LON_END   LONGITUDE OF END POINT,[deg]
      !> \param R         RANGE,[km]
      !> \param AZ        AZIMUTH,[deg]
      !> \param AZ_END    END AZIMUTH,[deg]
      !> \remark RESULT RETURNS IN R, AZ, AZ_END
      !> \warning SUBROUTINE USES COMMON BLOCK:
      !> \n       /CONSTS/
      !> \n       SEE SOURCE FOR DETAILS!
      !>
      SUBROUTINE GTOR
     * (N,LAT_START,LON_START,LAT_END,LON_END,R,AZ,AZ_END)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      INTEGER N      ! NUMBER OF ELLIPSOID
      REAL LAT_START ! LATITUDE OF START POINT,[deg]
      REAL LON_START ! LONGITUDE OF START POINT,[deg]
      REAL LAT_END   ! LATITUDE OF END POINT,[deg]
      REAL LON_END   ! LONGITUDE OF END POINT,[deg]
      REAL R         ! RANGE,[km]
      REAL AZ        ! AZIMUTH,[deg]
      REAL AZ_END    ! END AZIMUTH,[deg]
      !-------------------------------------------------------------------------
      ! LIST OF COMMON BLOCKS
      REAL(4) DGTORD,RDTODG
      COMMON /CONSTS/ C,RZ,DGTORD,RDTODG ! CONSTANTS
      !-------------------------------------------------------------------------
      ! DECLARATIONS
      REAL(8) COS_LAT_END,SIN_LAT_END,COS_LAT_START,SIN_LAT_START
     * ,SIN_LON_END_MINUS_LON_START,COS_LON_END_MINUS_LON_START,F1,F2,F3
     * ,DGTORD_D,RDTODG_D
      REAL(8) REARTH,A,B,F,L,U1,U2,SINU1,COSU1,SINU2,COSU2,
     * LAMBD,LAMBD_NEW,SINSIGMA,COSSIGMA,SIGMA,SINALPHA,
     * COSSQALPHA,COS2SIGMAM,SINLAMBDA,COSLAMBDA,CC,AA,BB,DELTASIGMA,
     * USQ,PI,TMPD,TMPD2,TO360_DOUBLE
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DGTORD_D=ASIN(1.D0)/90.D0
      RDTODG_D=90.D0/ASIN(1.D0)

      SELECTCASE(N)
      CASE(0,1,2)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !
      !                          SUPPOSE EARTH IS A SPHERE
      !
      IF(N.EQ.0) REARTH=DBLE(RZ)
      IF(N.EQ.1) REARTH=6371000.D0
      IF(N.EQ.2) REARTH=6378000.D0

      COS_LAT_END  =COS(DBLE(LAT_END)*DGTORD_D)
      SIN_LAT_END  =SIN(DBLE(LAT_END)*DGTORD_D)
      COS_LAT_START=COS(DBLE(LAT_START)*DGTORD_D)
      SIN_LAT_START=SIN(DBLE(LAT_START)*DGTORD_D)
      SIN_LON_END_MINUS_LON_START=SIN(DBLE(LON_END-LON_START)*DGTORD_D)
      COS_LON_END_MINUS_LON_START=COS(DBLE(LON_END-LON_START)*DGTORD_D)

      ! AZIMUTH
      F1=COS_LAT_END*SIN_LON_END_MINUS_LON_START
      F2=COS_LAT_START*SIN_LAT_END
      F3=SIN_LAT_START*COS_LAT_END*COS_LON_END_MINUS_LON_START
      AZ=TO360(REAL(ATAN2(F1,F2-F3)*RDTODG_D),'D')

      ! END AZIMUTH
      F1=COS_LAT_START*SIN_LON_END_MINUS_LON_START
      F2=COS_LAT_START*SIN_LAT_END*COS_LON_END_MINUS_LON_START
      F3=SIN_LAT_START*COS_LAT_END
      !AZ_END=TO360((ATAN2(F1,F2-F3)*RDTODG+180.),'D')
      AZ_END=TO360(REAL(ATAN2(F1,F2-F3)*RDTODG_D),'D')

      ! DISTANCE
      F1=SIN_LAT_START*SIN_LAT_END
      F2=COS_LAT_START*COS_LAT_END*COS_LON_END_MINUS_LON_START
      F3=F1+F2;
      R=REAL(ACOS(F3)*REARTH/1000.D0)

      CASE(3,4,5)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !
      !            SUPPOSE EARTH IS ELLIPSOID, USE VINCENTY FORMULAS
      !
      IF(N.EQ.3)THEN ! WGS-84
         A=6378137.0D0
         F=1.0D0 / 298.257223563D0
      ENDIF
      IF(N.EQ.4)THEN ! PZ-90
         A=6378136.0D0
         F=1.0D0 / 298.257839303D0
      ENDIF
      IF(N.EQ.5)THEN ! GRS-80
         A=6378137.0D0
         F=1.0D0 / 298.257222101D0
      ENDIF
      B=(1-F)*A
      PI=ACOS(-1.D0)
      ! ELLIPSOID DETERMINED, MAKE CALCULATIONS:
      L=DBLE(LON_END-LON_START)*DGTORD_D
      U1=ATAN((1.0D0-F)*TAN(DBLE(LAT_START)*DGTORD_D))
      U2=ATAN((1.0D0-F)*TAN(DBLE(LAT_END)*DGTORD_D))
      SINU1=SIN(U1)
      COSU1=COS(U1)
      SINU2=SIN(U2)
      COSU2=COS(U2)
      LAMBD=L
      LAMBD_NEW=0.D0
      ITER_LIMIT=100
      SINSIGMA=0.D0
      COSSIGMA=0.D0
      SIGMA=0.D0
      SINALPHA=0.D0
      COSSQALPHA=0.D0
      COS2SIGMAM=0.D0
      SINLAMBDA=0.D0
      COSLAMBDA=0.D0
      CE=0.D0

      !DO WHILE((ABS(LAMBD-LAMBD_NEW)/LAMBD).GT.1E-15)
1     CONTINUE
         SINLAMBDA=SIN(LAMBD)
         COSLAMBDA=COS(LAMBD)

         ! EQ.14
         SINSIGMA=SQRT(((COSU2*SINLAMBDA)*(COSU2*SINLAMBDA)+
     *                 (COSU1*SINU2-SINU1*COSU2*COSLAMBDA)*
     *                 (COSU1*SINU2-SINU1*COSU2*COSLAMBDA)))
         IF(ABS(SINSIGMA).LT.0.00000001D0) RETURN ! CO-INCIDENT POINTS

         ! EQ.15
         COSSIGMA=SINU1*SINU2+COSU1*COSU2*COSLAMBDA

         ! EQ.16
         SIGMA=ATAN2(SINSIGMA, COSSIGMA)

         ! EQ.17 Careful!  sin2sigma might be almost 0
         SINALPHA=COSU1*COSU2*SINLAMBDA/SINSIGMA;
         COSSQALPHA=1.0D0-SINALPHA*SINALPHA;

         ! EQ.18 Careful!  cos2alpha might be almost 0!
         COS2SIGMAM=COSSIGMA-2.0D0*SINU1*SINU2/COSSQALPHA !!!!!!!!!!!

         IF(ISNAN(COS2SIGMAM)) COS2SIGMAM=0.D0 ! equatorial line: cosSqAlpha = 0

         ! EQ.10
         CC=(F/16.D0)*COSSQALPHA*(4.D0+F*(4.D0-3.D0*COSSQALPHA))

         LAMBD_NEW=LAMBD

         ! EQ.11 (modified)
         LAMBD=L+(1.D0-CC)*F*SINALPHA*
     *         (SIGMA+CC*SINSIGMA*(COS2SIGMAM+CC*COSSIGMA*
     *         (-1.D0+2.D0*COS2SIGMAM*COS2SIGMAM)))

         ITER_LIMIT=ITER_LIMIT-1
         !IF(ITER_LIMIT.LT.0)EXIT ! ITERATIONS ARE OUT OF LIMIT !
         IF(ITER_LIMIT.LT.0)GOTO 2 ! ITERATIONS ARE OUT OF LIMIT !

      !ENDDO
      IF((ABS(LAMBD-LAMBD_NEW)/LAMBD).GT.1E-15.AND.
     *    ITER_LIMIT.GT.0) GOTO 1
2     CONTINUE

      USQ=COSSQALPHA*(A*A-B*B)/(B*B)

      ! EQ.3
      AA=1.D0+USQ/16384.D0*
     * (4096.D0+USQ*(-768.D0+USQ*(320.D0-175.D0*USQ)))

      ! EQ.4
      BB=USQ/1024.D0*(256.D0+USQ*(-128.D0+USQ*(74.D0-47.D0*USQ)))

      ! EQ.6
      DELTASIGMA=BB*SINSIGMA*(COS2SIGMAM+(BB/4.D0)*
     *           (COSSIGMA*(-1.D0+2.D0*COS2SIGMAM*COS2SIGMAM)-(BB/6.D0)*
     *           COS2SIGMAM*(-3.D0+4.D0*SINSIGMA*SINSIGMA)*
     *           (-3.D0+4.D0*COS2SIGMAM*COS2SIGMAM)))

      ! EQ.19
      R=REAL((B*AA*(SIGMA-DELTASIGMA))/1000.D0) ! [km]

      ! EQ.20
      TMPD=ATAN2((COSU2*SINLAMBDA),(COSU1*SINU2-SINU1*COSU2*COSLAMBDA))
      AZ=REAL(RDTODG_D*TO360_DOUBLE(TMPD,'R'))

      ! EQ.21
      TMPD2=
     * ATAN2((COSU1*SINLAMBDA),(-SINU1*COSU2+COSU1*SINU2*COSLAMBDA))   ! +PI - NEEDED TO FIND REVERSE
                                                                       !AZIMUTH, NOT DIRECT
      AZ_END=REAL(RDTODG_D*TO360_DOUBLE(TMPD2,'R'))
      ENDSELECT
      RETURN
      END

      !#########################################################################
      !>
      !> \brief CONVERTS RADIOLOCATION COORDINATES TO GEOGRAPHIC COORDINATES (DIRECT TASK)
      !>
      !> \param N   NUMBER OF ELLIPSOID:
      !> \n         0 : SPHERE WITH RADIUS RZ FROM COMMON BLOCK
      !> \n         1 : SPHERE R=6371 [km] (AVERAGE RADIUS)
      !> \n         2 : SPHERE R=6378 [km] (EQUATORIAL RADIUS)
      !> \n         3 : ELLIPSOID WGS-84 (6378137.0 [m], 1 / 298.257223563)
      !> \n         4 : ELLIPSOID PZ-90  (6378136.0 [m], 1 / 298.257839303)
      !> \n         5 : ELLIPSOID GRS-80 (6378137.0 [m], 1 / 298.257222101)
      !> \param LAT_START LATITUDE OF START POINT,[deg]
      !> \param LON_START LONGITUDE OF START POINT,[deg]
      !> \param R         RANGE,[km]
      !> \param AZ        AZIMUTH,[deg]
      !> \param LAT_END   LATITUDE OF END POINT,[deg]
      !> \param LON_END   LONGITUDE OF END POINT,[deg]
      !> \param AZ_END    END AZIMUTH,[deg]
      !> \remark RESULT RETURNS IN LAT_END, LON_END, AZ_END
      !> \warning SUBROUTINE USES COMMON BLOCK:
      !> \n       /CONSTS/
      !> \n       SEE SOURCE FOR DETAILS!
      !>
      SUBROUTINE RTOG
     * (N,LAT_START,LON_START,R,AZ,LAT_END,LON_END,AZ_END)!
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      INTEGER N      ! NUMBER OF ELLIPSOID
      REAL LAT_START ! LATITUDE OF START POINT,[deg]
      REAL LON_START ! LONGITUDE OF START POINT,[deg]
      REAL R         ! RANGE,[km]
      REAL AZ        ! AZIMUTH,[deg]
      REAL LAT_END   ! LATITUDE OF END POINT,[deg]
      REAL LON_END   ! LONGITUDE OF END POINT,[deg]
      REAL AZ_END    ! END AZIMUTH,[deg]
      !-------------------------------------------------------------------------
      ! LIST OF COMMON BLOCKS
      REAL(4) DGTORD,RDTODG
      COMMON /CONSTS/ C,RZ,DGTORD,RDTODG ! CONSTANTS
      !-------------------------------------------------------------------------
      ! DECLARATIONS
      REAL(8) RNG,SIN_LAT_START,COS_LAT_START,COS_RNG,SIN_RNG,COS_AZ
     * ,F1,F2,F3
     * ,DGTORD_D,RDTODG_D
      REAL(8) PI,REARTH,A,B,F,S,ALPHA1,SINALPHA1,COSALPHA1
     * ,TANU1,COSU1,SINU1,SIGMA1,SINALPHA,COSSQALPHA,USQ,AA,BB
     * ,SIGMA,PREVSIGMA,COS2SIGMAM,SINSIGMA,COSSIGMA,DELTASIGMA,TMP
     * ,TETTA,LAMBD,CC,L,PHI,ALPHA2,BZV,TO360_DOUBLE
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DGTORD_D=ASIN(1.D0)/90.D0
      RDTODG_D=90.D0/ASIN(1.D0)
      REARTH=0.D0

      SELECTCASE(N)
      CASE(0,1,2)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !
      !                          SUPPOSE EARTH IS A SPHERE
      !
      IF(N.EQ.0) REARTH=DBLE(RZ)
      IF(N.EQ.1) REARTH=6371000.D0
      IF(N.EQ.2) REARTH=6378000.D0

      RNG=DBLE(R)*1000.D0/REARTH ! CONVERT RANGE TO A SINGLE SPHERE

      SIN_LAT_START=SIN(DBLE(LAT_START)*DGTORD_D)
      COS_LAT_START=COS(DBLE(LAT_START)*DGTORD_D)
      COS_RNG=COS(RNG)
      SIN_RNG=SIN(RNG)
      COS_AZ=COS(DBLE(AZ)*DGTORD_D)

      ! LATITUDE
      F1=SIN_LAT_START*COS_RNG
      F2=COS_LAT_START*SIN_RNG*COS_AZ
      LAT_END=REAL(ASIN(F1+F2)*RDTODG_D)

      ! LONGITUDE
      F1=SIN_RNG*SIN(DBLE(AZ)*DGTORD_D)
      F2=COS_LAT_START*COS_RNG
      F3=SIN_LAT_START*SIN_RNG*COS_AZ
      LON_END=REAL(DBLE(LON_START)+ATAN2(F1,F2-F3)*RDTODG_D)

      ! END AZIMUTH
      F1=COS_LAT_START*SIN(DBLE(AZ)*DGTORD_D)
      F2=COS_LAT_START*COS_RNG*COS_AZ
      F3=SIN_LAT_START*SIN_RNG
      AZ_END=TO360(REAL(ATAN2(F1,F2-F3)*RDTODG_D),'D') ! FINAL DIRECT BEARING

      CASE(3,4,5)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !
      !            SUPPOSE EARTH IS ELLIPSOID, USE VINCENTY FORMULAS
      !
      IF(N.EQ.3)THEN ! WGS-84
         A=6378137.0D0
         F=1.0D0 / 298.257223563D0
      ENDIF
      IF(N.EQ.4)THEN ! PZ-90
         A=6378136.0D0
         F=1.0D0 / 298.257839303D0
      ENDIF
      IF(N.EQ.5)THEN ! GRS-80
         A=6378137.0D0
         F=1.0D0 / 298.257222101D0
      ENDIF
      B=(1-F)*A
      PI=ACOS(-1.D0)
      ! ELLIPSOID DETERMINED, MAKE CALCULATIONS:

      S=R*1000. ! [m]
      ALPHA1=DBLE(AZ)*DGTORD_D

      SINALPHA1=SIN(ALPHA1)
      COSALPHA1=COS(ALPHA1)

      TANU1=(1.D0-F)*TAN(DBLE(LAT_START)*DGTORD_D)
      COSU1=1.D0 / SQRT(1.D0+TANU1*TANU1)
      SINU1=TANU1*COSU1

      ! EQ.1
      SIGMA1=ATAN2(TANU1, COSALPHA1)

      ! EQ.2
      SINALPHA=COSU1*SINALPHA1
      COSSQALPHA=1.D0-SINALPHA*SINALPHA

      USQ=COSSQALPHA*(A*A-B*B)/(B*B)

      ! EQ.3
      AA=1.D0+
     * (USQ/16384.D0)*(4096.D0+USQ*(-768.D0+USQ*(320.D0-175.D0*USQ)))

      ! EQ.4
      BB=(USQ/1024.D0)*(256.D0+USQ*(-128.D0+USQ*(74.D0-47.D0*USQ)))

      ! ITERATE UNTIL THERE IS A NEGLIGIBLE CHANGE IN SIGMA
      SIGMA=S/(B*AA)
      PREVSIGMA=SIGMA

      COS2SIGMAM=0.D0
      SINSIGMA=0.D0
      COSSIGMA=0.D0
      DELTASIGMA=0.D0

      DO WHILE(.TRUE.)

      ! EQ.5
      COS2SIGMAM=COS((2.D0*SIGMA1+SIGMA))
      SINSIGMA=SIN(SIGMA)
      COSSIGMA=COS(SIGMA)

      ! EQ.6
      DELTASIGMA=BB*SINSIGMA*(COS2SIGMAM+(BB/4.D0)*(COSSIGMA*
     * (-1.D0+2.D0*COS2SIGMAM*COS2SIGMAM)-(BB/6.D0)*COS2SIGMAM*
     * (-3.D0+4.D0*SINSIGMA*SINSIGMA)*
     * (-3.D0+4.D0*COS2SIGMAM*COS2SIGMAM)))

      ! EQ.7
      SIGMA=(S/(B*AA))+DELTASIGMA

      ! BREAK AFTER CONVERGING TO TOLERANCE
      IF((ABS(SIGMA-PREVSIGMA).LT.1E-15).OR.
     *   ISNAN(SIGMA))EXIT

      PREVSIGMA=SIGMA
      ENDDO

      COS2SIGMAM=COS(2.D0*SIGMA1+SIGMA)
      SINSIGMA=SIN(SIGMA)
      COSSIGMA=COS(SIGMA)

      TMP=SINU1*SINSIGMA-COSU1*COSSIGMA*COSALPHA1

      ! EQ.8
      TETTA=ATAN2(SINU1*COSSIGMA+COSU1*SINSIGMA*COSALPHA1,
     * (1.D0-F)*SQRT((SINALPHA*SINALPHA+TMP*TMP)))

      ! EQ.9
      LAMBD=ATAN2((SINSIGMA*SINALPHA1),
     * (COSU1*COSSIGMA-SINU1*SINSIGMA*COSALPHA1))

      ! EQ.10
      CC=(F/16.D0)*COSSQALPHA*(4.D0+F*(4.D0-3.D0*COSSQALPHA))

      ! EQ.11
      L=LAMBD-(1.D0-CC)*F*SINALPHA*(SIGMA+CC*SINSIGMA*
     * (COS2SIGMAM+CC*COSSIGMA*(-1.D0+2.D0*COS2SIGMAM*COS2SIGMAM)))

      !PHI=(LON_START*DGTORD+L+3.D0*PI)%(2.D0*PI)-PI; ! TO -180.. 180 original!
      PHI=MOD((LON_START*DGTORD+L+3.D0*PI),(2.D0*PI))-PI

      ! EQ.12
      ALPHA2=ATAN2(SINALPHA,-TMP) ! FINAL BEARING, IF REQUIRED

      BZV=TO360_DOUBLE(ALPHA2,'R') ! DIRECT AZIMUTH AT THE END POINT (FINAL BEARING)

      ! RESULT
      LAT_END=REAL(TETTA*RDTODG_D)
      LON_END=REAL(PHI*RDTODG_D)
      AZ_END=REAL(BZV*RDTODG_D)

      ENDSELECT
      RETURN
      END
