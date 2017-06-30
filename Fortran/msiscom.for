*
*      ----------------------------------------------------------------
*
*                              MSISCOM.FOR
*
*  This file contains common routines between the MSIS models.
*
*                          Companion code for
*             Fundamentals of Astrodynamics and Applications
*                                  2007
*                            by David Vallado
*
*     (W) 719-573-2600, email dvallado@agi.com
*
*     *****************************************************************
*
*  Current :
*            14 Feb 03  David Vallado
*                         Misc updates
*  Changes :
*            28 May 02  David Vallado
*                         Fix densu per L Schwartz correction - tests correct
*            15 Apr 02  David Vallado
*                         Work on gsurf
*            24 Feb 02  David Vallado
*                         Finish double conversion
*            20 Jan 02  David Vallado
*                         Original baseline
*
*     *****************************************************************
*
*     Uses object files:
*       None
*
*      Common to ALL MSIS Models:
*
*      SUBROUTINE TSELEC (SV)
*
*      SUBROUTINE TRETRV (SVV)
*
*      REAL*8 FUNCTION DNET (DD, DM, ZHM, XMM, XM)
*
*      REAL*8 FUNCTION CCOR ( ALT, R, H1, ZH)
*
*      SUBROUTINE METERS(METER)
*
*      Common to 90 and 00 MSIS Models:
*
*      REAL*8 FUNCTION DENSU (ALT, DLB, TINF, TLB, XM, ALPHA, TZ, ZLB,
*     &                       S2, MN1, ZN1, TN1, TGN1)
*
*      REAL*8 FUNCTION DENSM (ALT, D0, XM, TZ, MN3, ZN3, TN3, TGN3, MN2,
*     &                       ZN2, TN2, TGN2)
*
*
*      SUBROUTINE SPLINE(X, Y, N, YP1, YPN, Y2)
*
*      SUBROUTINE SPLINT(XA, YA, Y2A, N, X, Y)
*
*      SUBROUTINE SPLINI(XA, YA, Y2A, N, X, YI)
*
*      SUBROUTINE GLATF (LAT, GV, REFF)
*
*      REAL*8 FUNCTION VTST (IYD, SEC, GLAT, GLONG, STL, F107A, F107,
*     &                      AP, IC)
*
*
*
*     *****************************************************************
*        Set up routines that are common to all the MSIS models
*     *****************************************************************

C----------------------------------------------------------------------
C     SET SWITCHES
C     Output in  COMMON/CSW/SW(25),ISW,SWC(25)
C     SW FOR MAIN TERMS, SWC FOR CROSS TERMS
C
C     TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SV),
C     WHERE SV IS A 25 ELEMENT ARRAY CONTAINING 0.0D0 FOR OFF, 1.0D0
C     FOR ON, OR 2.0D0 FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C
C     To get current values of SW: CALL TRETRV(SW)
C----------------------------------------------------------------------

      SUBROUTINE TSELEC (SV)
        IMPLICIT NONE
        REAL*8 SV(1)

        REAL*8 SAV(25), SVV(1)
        INTEGER I

        REAL*8  SW(25), SWC(25)
        INTEGER ISW
        COMMON/CSW/SW,SWC,ISW

        SAVE   ! may not need this statement
c------------------------------ begin ---------------------------------
        DO I = 1,25
            SAV(I) = SV(I)
            SW(I) = DMOD(SV(I),2.0D0)
            IF (DABS(SV(I)).EQ.1.0D0 .OR. DABS(SV(I)).EQ.2.0D0) THEN
                SWC(I) = 1.0D0
              ELSE
                SWC(I) = 0.0D0
              ENDIF
          ENDDO

        ISW = 64999

      RETURN
      END

C----------------------------------------------------------------------

      SUBROUTINE TRETRV (SVV)
        IMPLICIT NONE
        INTEGER i
        REAL*8 SVV(25), SAV(25)

        Save

        DO I = 1,25
            SVV(I) = SAV(I)
        ENDDO

      RETURN
      END
*
C----------------------------------------------------------------------
C     TURBOPAUSE CORRECTION FOR MSIS MODELS
C     Eq. A12b
C         Root mean density
C          DD   - diffusive density
C          DM   - full mixed density
C          ZHM  - transition scale length
C          XMM  - full mixed molecular weight
C          XM   - species molecular weight
C          DNET - combined density
C       8/20/80
C----------------------------------------------------------------------

      REAL*8 FUNCTION DNET (DD, DM, ZHM, XMM, XM)
        IMPLICIT NONE
        REAL*8 DD, DM, XM, XMM, ZHM

        REAL*8 A, YLOG
C------------------------------ begin ---------------------------------
        A = ZHM/(XMM-XM)

        IF (DM.le.0.0D0 .or. DD.le.0.0D0 ) THEN
            WRITE(*,*) 'DNET LOG ERROR',DM,DD,XM

            IF (DD.EQ.0.0D0 .AND. DM.EQ.0.0D0) DD = 1.0D0
            IF (DM.EQ.0.0D0) DNET = DD
            IF (DD.EQ.0.0D0) DNET = DM

          ELSE
            ! ---- Eq. A12a ----
            YLOG = A*DLOG(DM/DD)
            IF (YLOG.LT.-10.0D0) THEN
                DNET = DD
              ELSE
                IF (YLOG.GT.10.0D0) THEN
                    DNET = DM
                  ELSE
                    DNET = DD*(1.0D0 + DEXP(YLOG))**(1.0D0/A)
                  ENDIF
              ENDIF
          ENDIF

      RETURN
      END

C----------------------------------------------------------------------
C     CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
C        ALT - altitude
C        R   - target ratio
C        H1  - transition scale length
C        ZH  - altitude of 1/2 R
C     Eq. A20a or Eq. A21
C----------------------------------------------------------------------

      REAL*8 FUNCTION CCOR ( ALT, R, H1, ZH)
        IMPLICIT NONE
        REAL*8 ALT, R, ZH

        REAL*8 E, EX, H1
C------------------------------ begin ---------------------------------
        E = (ALT-ZH)/H1
        IF (E.GT.70.0D0) THEN
            CCOR = 0.0D0
          ELSE
            IF (E.LT.-70.0D0) THEN
                CCOR = R
              ELSE
                EX = DEXP(E)
                CCOR = R/(1.0D0 + EX)
              ENDIF
          ENDIF

        CCOR = DEXP(CCOR)
      RETURN
      END

*
C-----------------------------------------------------------------------
C      Convert outputs to Kg & Meters if METER true
C-----------------------------------------------------------------------

      SUBROUTINE METERS (METER)
        IMPLICIT NONE
        LOGICAL METER

        INTEGER IMR
        COMMON/METSEL/IMR

C------------------------------ begin ---------------------------------
        IMR=0
        IF (METER) THEN
            IMR=1
          ENDIF

      RETURN
      END
*
*     *****************************************************************
*     *****************************************************************
*       Set up routines that are common to the MSIS-90 and 00 models
*     *****************************************************************

C----------------------------------------------------------------------
C       Calculate Temperature and Density Profiles for MSIS models
C       New lower thermo polynomial 10/30/89
c
C----------------------------------------------------------------------

      REAL*8 FUNCTION DENSU (ALT, DLB, TINF, TLB, XM, ALPHA, TZ, ZLB,
     &                       S2, MN1, ZN1, TN1, TGN1)
        IMPLICIT NONE
        REAL*8 alt, dlb, tinf, tlb, xm, alpha, tz, zlb, s2, zn1(5),
     &         tn1(5), tgn1(2)
        INTEGER mn1

        REAL*8 RGAS, Zeta, ZZ, ZL, ZA, Z, ZG2, DTA, T1, T2,
     &         ZGDif, ZG, XS(5), Z1, YD1, YD2, Y, Y2out(5), Glb, Gamma,
     &         Expl, Densa, YS(5), Gamm, TT, TA, Z2, X, YI

        INTEGER MN, K

c        DIMENSION ZN1(MN1),TN1(MN1)

        REAL*8 GSurf, RE
        COMMON/PARMB90/GSURF,RE

        REAL*8  QPB(50), DV(60)
        INTEGER MP, II, JG, LT, IErr, IFun, N, J
        COMMON/LSQV90/ MP, II, JG, LT, IERR, IFUN, N, J, QPB, DV

        SAVE

        DATA RGAS/831.4D0/

c       Statement Function
        ZETA(ZZ,ZL) = (ZZ-ZL)*(RE+ZL)/(RE+ZZ)

C------------------------------ begin ---------------------------------
C       WRITE(6,*) 'DB',ALT,DLB,TINF,TLB,XM,ALPHA,ZLB,S2,MN1,ZN1,TN1
        DENSU = 1.0D0

        ! -------- Joining altitude of Bates and spline
        ZA = ZN1(1)
        Z  = AMAX1(ALT,ZA)

        ! -------- Geopotential altitude difference from ZLB
        ZG2 = ZETA(Z,ZLB)

        ! -------- Bates temperature
        TT = TINF-(TINF-TLB)*DEXP(-S2*ZG2)
        TA = TT
        TZ = TT
        DENSU = TZ

      ! --------------- CALCULATE TEMPERATURE BELOW ZA ----------------
      ! -------- Temperature gradient at ZA from Bates profile
      IF (ALT.lt.ZA) THEN
          DTA = (TINF-TA)*S2*((RE+ZLB)/(RE+ZA))**2
          TGN1(1) = DTA
          TN1(1) = TA
          Z  = AMAX1(ALT,ZN1(MN1))
          MN = MN1
          Z1 = ZN1(1)
          Z2 = ZN1(MN)
          T1 = TN1(1)
          T2 = TN1(MN)

          ! -------- Geopotental difference from Z1
          ZG = ZETA(Z,Z1)
          ZGDIF = ZETA(Z2,Z1)

          ! -------- Set up spline nodes
          DO K = 1,MN
              XS(K) = ZETA(ZN1(K),Z1)/ZGDIF
              YS(K) = 1.0D0/TN1(K)
            ENDDO

          ! -------- End node derivatives
          YD1 = -TGN1(1)/(T1*T1)*ZGDIF
          YD2 = -TGN1(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2

          ! -------- Calculate spline coefficients
          CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
          X = ZG/ZGDIF
          CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)

          ! -------- temperature at altitude
          TZ = 1.0D0/Y
          DENSU = TZ
cdav
        ENDIF

          ! ------------ CALCULATE DENSITY ABOVE ZA -------------------
          IF (XM.ne.0.0D0) THEN
              GLB  = GSURF  / (1.0D0 + ZLB/RE)**2
              GAMMA= XM*GLB / (S2*RGAS*TINF)
              EXPL = DEXP(-S2*GAMMA*ZG2)
              IF (EXPL.GT.50 .OR. TT.LE.0.0D0) THEN
                EXPL = 50.0D0
              ENDIF
              ! -------- Density at altitude
              DENSA = DLB*(TLB/TT)**(1.0D0+ALPHA+GAMMA)*EXPL
              DENSU = DENSA

              IF (ALT.lt.ZA) THEN
                  ! ---------- CALCULATE DENSITY BELOW ZA -------------
                  GLB  = GSURF/(1.0D0 + Z1/RE)**2
                  GAMM = XM*GLB*ZGDIF/RGAS

                  ! -------- integrate spline temperatures
                  CALL SPLINI (XS,YS,Y2OUT,MN,X,YI)
                  EXPL = GAMM*YI
                  IF (EXPL.GT.50.0D0.OR.TZ.LE.0.0D0) THEN
                      EXPL = 50.0D0
                    ENDIF

                  ! -------- Density at altitude
                  DENSU = DENSU*(T1/TZ)**(1.0D0+ALPHA)*DEXP(-EXPL)
                ENDIF
            ENDIF

      RETURN
      END
*
C----------------------------------------------------------------------
c
c
c Calculate Temperature and Density Profiles for lower atmos.
c
C----------------------------------------------------------------------

      REAL*8 FUNCTION DENSM (ALT, D0, XM, TZ, MN3, ZN3, TN3, TGN3, MN2,
     &                       ZN2, TN2, TGN2)
        IMPLICIT NONE
        REAL*8 Alt, D0, Xm, Tz, ZN3(5), TN3(5), TGN3(2), ZN2(4), TN2(4),
     &         TGN2(2)
        INTEGER MN3, MN2

        REAL*8 RGAS, Zeta, ZZ, ZL, Z, T1, T2,
     &         ZGDif, ZG, XS(10), Z1, YD1, YD2, Y, Y2out(10), Glb,
     &         Ys(10), Gamm, Z2, X, YI, EXPL
        INTEGER MN, K

c        DIMENSION ZN3(MN3),TN3(MN3)
c        DIMENSION ZN2(MN2),TN2(MN2)

        REAL*8 GSurf, RE
        COMMON/PARMB90/GSURF,RE

        REAL*8 TAF
        COMMON/FIT90/TAF

        REAL*8  QPB(50), DV(60)
        INTEGER MP, II, JG, LT, IErr, IFun, N, J
        COMMON/LSQV90/ MP, II, JG, LT, IERR, IFUN, N, J, QPB, DV

        SAVE

        DATA RGAS/831.4D0/

c       Statement Function
        ZETA(ZZ,ZL) = (ZZ-ZL)*(RE+ZL)/(RE+ZZ)

C------------------------------ begin ---------------------------------
        DENSM = D0

      IF (ALT.le.ZN2(1)) THEN
          ! -------- STRATOSPHERE/MESOSPHERE TEMPERATURE
          Z = AMAX1(ALT,ZN2(MN2))
          MN = MN2
          Z1 = ZN2(1)
          Z2 = ZN2(MN)
          T1 = TN2(1)
          T2 = TN2(MN)
          ZG = ZETA(Z,Z1)
          ZGDIF = ZETA(Z2,Z1)

          ! -------- Set up spline nodes
          DO K = 1,MN
              XS(K) = ZETA(ZN2(K),Z1)/ZGDIF
              YS(K) = 1.0D0/TN2(K)
            ENDDO

          YD1 = -TGN2(1)/(T1*T1)*ZGDIF
          YD2 = -TGN2(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2

          ! -------- Calculate spline coefficients
          CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
          X = ZG/ZGDIF
          CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)

          ! -------- Temperature at altitude
          TZ = 1.0D0/Y

          IF (XM.ne.0.0D0) THEN
              ! ------ CALCULATE STRATOSPHERE/MESOSPHERE DENSITY ------
              GLB = GSURF/(1.0D0+Z1/RE)**2
              GAMM = XM*GLB*ZGDIF/RGAS

              ! -------- Integrate temperature profile
              CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
              EXPL = GAMM*YI
              IF (EXPL.GT.50.0D0) EXPL = 50.0D0

              ! -------- Density at altitude
              DENSM = DENSM*(T1/TZ)*DEXP(-EXPL)
            ENDIF

          IF (ALT.le.ZN3(1)) THEN
              ! -------- TROPOSPHERE/STRATOSPHERE TEMPERATURE ---------
              Z = ALT
              MN = MN3
              Z1 = ZN3(1)
              Z2 = ZN3(MN)
              T1 = TN3(1)
              T2 = TN3(MN)
              ZG = ZETA(Z,Z1)
              ZGDIF = ZETA(Z2,Z1)

              ! -------- Set up spline nodes
              DO K = 1,MN
                  XS(K) = ZETA(ZN3(K),Z1)/ZGDIF
                  YS(K) = 1.0D0/TN3(K)
                ENDDO
              YD1 = -TGN3(1)/(T1*T1)*ZGDIF
              YD2 = -TGN3(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2

              ! -------- Calculate spline coefficients
              CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
              X = ZG/ZGDIF
              CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)

              ! -------- temperature at altitude
              TZ = 1.0D0/Y
              IF (XM.ne.0.0D0) THEN
                  ! --- CALCULATE TROPOSPHERIC/STRATOSPHERE DENSITY ---
                  GLB = GSURF/(1.0D0+Z1/RE)**2
                  GAMM = XM*GLB*ZGDIF/RGAS

                  ! -------- Integrate temperature profile
                  CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
                  EXPL = GAMM*YI
                  IF (EXPL.GT.50.0D0) EXPL = 50.0D0

                  ! -------- Density at altitude
                  DENSM = DENSM*(T1/TZ)*DEXP(-EXPL)
                ENDIF
            ENDIF
        ENDIF

      IF (XM.EQ.0.0D0) THEN
          DENSM = TZ
        ENDIF

      RETURN
      END
*
C----------------------------------------------------------------------
c  CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
c  ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
c  X,Y    : ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
c  N      : SIZE OF ARRAYS X,Y
c  YP1,YPN: SPECIFIED DERIVATIVES AT X(1) AND X(N); VALUES
c           > =  1D30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
c  Y2     : OUTPUT ARRAY OF SECOND DERIVATIVES
C----------------------------------------------------------------------

      SUBROUTINE SPLINE(X, Y, N, YP1, YPN, Y2)
        IMPLICIT NONE
        REAL*8 X(5), Y(5), YP1, YPN, Y2(5), U
        INTEGER NMAX,N, I, K

        REAL*8 Sig, P, Qn, Un
        PARAMETER (NMAX = 100)
c        DIMENSION X(N),Y(N),Y2(N),U(NMAX)
        DIMENSION U(NMAX)

c        SAVE

C------------------------------ begin ---------------------------------
        IF (YP1.GT.0.99D30) THEN
            Y2(1) = 0.0D0
            U(1)  = 0.0D0
          ELSE
            Y2(1) = -0.5D0
            U(1)  = (3.0D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
          ENDIF

        DO I = 2,N-1
            SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
            P   = SIG*Y2(I-1) + 2.0D0
            Y2(I) = (SIG-1.0D0)/P
            U(I)  = (6.0D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     &             /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1)) / P
          ENDDO

        IF (YPN.GT.0.99D30) THEN
            QN = 0.0D0
            UN = 0.0D0
          ELSE
            QN = 0.5D0
            UN = (3.0D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/
     &           (X(N)-X(N-1)))
          ENDIF

        Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.0D0)
        DO K = N-1,1,-1
            Y2(K) = Y2(K)*Y2(K+1)+U(K)
          ENDDO
      RETURN
      END
*
C----------------------------------------------------------------------
c  CALCULATE CUBIC SPLINE INTERP VALUE
c  ADAPTED FROM NUMBERICAL RECIPES BY PRESS ET AL.
c  XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
c  Y2A  : ARRAY OF SECOND DERIVATIVES
c  N    : SIZE OF ARRAYS XA,YA,Y2A
c  X    : ABSCISSA FOR INTERPOLATION
c  Y    : OUTPUT VALUE
C----------------------------------------------------------------------

      SUBROUTINE SPLINT(XA, YA, Y2A, N, X, Y)
        IMPLICIT NONE
        REAL*8 XA(5), YA(5), Y2A(5)
        INTEGER N, K, KHI, KLO

        REAL*8 H, A, B, Y, X
c        DIMENSION XA(N),YA(N),Y2A(N)

c        SAVE

C------------------------------ begin ---------------------------------
        KLO = 1
        KHI = N
        DO WHILE (KHI-KLO.GT.1)
            K = (KHI+KLO) * 0.5D0
            IF (XA(K).GT.X) THEN
                KHI = K
              ELSE
                KLO = K
              ENDIF
          ENDDO

        H = XA(KHI)-XA(KLO)
        IF (H.EQ.0.0D0) THEN
            WRITE(*,*) 'BAD XA INPUT TO SPLINT'
          ENDIF
        A = (XA(KHI)-X)/H
        B = (X-XA(KLO))/H
        Y = A*YA(KLO)+B*YA(KHI)+
     &      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*H*H/6.0D0

      RETURN
      END
*
C----------------------------------------------------------------------
C     INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
c  XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
c  Y2A  : ARRAY OF SECOND DERIVATIVES
c  N    : SIZE OF ARRAYS XA,YA,Y2A
c  X    : ABSCISSA ENDPOINT FOR INTEGRATION
c  Y    : OUTPUT VALUE
C----------------------------------------------------------------------

      SUBROUTINE SPLINI(XA, YA, Y2A, N, X, YI)
        IMPLICIT NONE
        REAL*8 XA(5), YA(5), Y2A(5), X, YI
        INTEGER N, KHI, KLO

        REAL*8 H, A, B, XX, A2, B2
c        DIMENSION XA(N),YA(N),Y2A(N)

c        SAVE

C------------------------------ begin ---------------------------------
        YI = 0
        KLO = 1
        KHI = 2

        DO WHILE (X.GT.XA(KLO) .AND. KHI.LE.N)
            XX = X
            IF (KHI.LT.N) THEN
                XX = AMIN1(X,XA(KHI))
              ENDIF
            H  = XA(KHI)-XA(KLO)
            A  = (XA(KHI)-XX)/H
            B  = (XX-XA(KLO))/H
            A2 = A*A
            B2 = B*B
            YI = YI + ((1.0D0 - A2)*YA(KLO)*0.5D0 + B2*YA(KHI)*0.5D0 +
     &           ((-(1.0D0 + A2*A2)*0.25D0 + A2*0.5D0)*Y2A(KLO) +
     &           (B2*B2/4.0D0 - B2*0.5D0)*Y2A(KHI))*H*H/6.0D0) * H
            KLO = KLO + 1
            KHI = KHI + 1
          ENDDO

      RETURN
      END

C----------------------------------------------------------------------
C      CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
C      RADIUS (REFF)
c
cdav chg to *8
C----------------------------------------------------------------------

      SUBROUTINE GLATF (LAT, GV, REFF)
      IMPLICIT NONE
        REAL*8 LAT, GV, REFF

        REAL*8 LATL, DGTR, C2

        DATA DGTR/1.74533D-2/,LATL/-999.0D0/

C------------------------------ begin ---------------------------------
        IF(LAT.NE.LATL) THEN
            C2 = DCOS(2.0D0*DGTR*LAT)
          ENDIF
        LATL = LAT
        GV   = 980.616D0*(1.0D0-0.0026373D0*C2)
        REFF = 2.0D0*GV/(3.085462D-6 + 2.27D-9*C2)*1.0D-5

      RETURN
      END
*
C----------------------------------------------------------------------
C       Test if geophysical variables or switches changed and save
C       Return 0 if unchanged and 1 if changed
C----------------------------------------------------------------------

      REAL*8 FUNCTION VTST (IYD, SEC, GLAT, GLONG, STL, F107A, F107,
     &                      AP, IC)
        IMPLICIT NONE
        REAL*8 Sec, GLAT, GLong, STL, F107A, F107, Ap(7),
     &         SWCL(25,2), APL(7,2), FL(2), Secl(2), GLatl(2), Gll(2)

        REAL*8 STLL(2), FAL(2), SWL(25,2)
        INTEGER IYD, IC, I, IYDL(2)

        REAL*8  SW(25), SWC(25)
        INTEGER ISW
        COMMON/CSW/SW,SWC,ISW

c        SAVE

        DATA IYDL/2*-999/,SECL/2*-999.0D0/,GLATL/2*-999.0D0/,
     &       GLL/2*-999.0D0/
        DATA STLL/2*-999.0D0/,FAL/2*-999.0D0/,FL/2*-999.0D0/,
     &       APL/14*-999.0D0/
        DATA SWL/50*-999.0D0/,SWCL/50*-999.0D0/

c------------------------------ begin ---------------------------------
        VTST = 0
        IF(  IYD.NE.IYDL(IC))  GOTO 10
        IF(  SEC.NE.SECL(IC))  GOTO 10
        IF( GLAT.NE.GLATL(IC)) GOTO 10
        IF(GLONG.NE.GLL(IC))   GOTO 10
        IF(  STL.NE.STLL(IC))  GOTO 10
        IF(F107A.NE.FAL(IC))   GOTO 10
        IF( F107.NE.FL(IC))    GOTO 10

        DO I = 1,7
            IF(AP(I).NE.APL(I,IC)) GOTO 10
          ENDDO

        DO I = 1,25
            IF(SW(I).NE.SWL(I,IC)) GOTO 10
            IF(SWC(I).NE.SWCL(I,IC)) GOTO 10
          ENDDO

        GOTO 20

   10 CONTINUE

        VTST = 1

        IYDL(IC)  = IYD
        SECL(IC)  = SEC
        GLATL(IC) = GLAT
        GLL(IC)   = GLONG
        STLL(IC)  = STL
        FAL(IC)   = F107A
        FL(IC)    = F107

        DO I = 1,7
            APL(I,IC)=AP(I)
          ENDDO

        DO I = 1,25
            SWL(I,IC) = SW(I)
            SWCL(I,IC) = SWC(I)
          ENDDO

   20   CONTINUE

      RETURN
      END

