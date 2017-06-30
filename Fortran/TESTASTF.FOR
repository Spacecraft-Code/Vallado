*
*      -----------------------------------------------------------------
*
*                              TESTASTF.FOR
*
*   This file contains sample subroutines to test each of the functions
*   in the astrodynamic libraries.
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
*            28 Jan 04  David Vallado
*                         Fix headers
*  Changes :
*            28 Feb 03  David Vallado
*                         New baseline
*            14 May 01  David Vallado
*                         2nd edition baseline
*            23 Nov 87  David Vallado
*                         Original baseline
*
*     *****************************************************************
*
*     Uses all ast, etc object files
*
      PROGRAM TESTASTF
        IMPLICIT NONE
        INTEGER i, NumChose
        CHARACTER*12 IFileName
        CHARACTER*4 Blank4
        CHARACTER*15 Title

        INCLUDE 'astmath.cmn'
        INCLUDE 'astreduc.cmn'
        INCLUDE 'astconst.cmn'

c        Write(*,*) 'Input Filename '
c        Read(*,*) IFileName
        IFileName = 'testastf.dat'
        OPEN( 10,FILE=IFileName,ACCESS='SEQUENTIAL',STATUS='OLD' )
        OPEN( 20,FILE='testastf.out',ACCESS='SEQUENTIAL',
     &        STATUS='UNKNOWN' )

        i = 1
        DO WHILE (i .ge. 0)
            ! -------------------- Tests for Math ----------------------
            Read( 10,'(I4,4X,A15)', END=999 ) NumChose,Title
            IF (NumChose .eq. 0 ) THEN
                Read(10,'(I4,4X,A15)',END=999 ) NumChose,Title
              ENDIF
            Write(*,*) ' ------ xx',Title,' Test Case --'
            Write(20,*) '-- Test Case --','------ xx',NumChose,Title,
     &                           '-----------------------'

            IF (NumChose .eq. 2 ) THEN
                CALL TestGETPART
              ENDIF
            IF (NumChose .eq. 3 ) THEN
                CALL TestGETPARTL
              ENDIF
            IF (NumChose .eq. 4 ) THEN
                CALL TestGetPartR
              ENDIF
            IF (NumChose .eq. 5 ) THEN
                CALL TestFACTORIAL
              ENDIF
            IF (NumChose .eq. 6 ) THEN
                CALL TestBINOMIAL
              ENDIF
            IF (NumChose .eq. 7 ) THEN
                CALL TestMIN
              ENDIF
            IF (NumChose .eq. 8 ) THEN
                CALL TestMAX
              ENDIF
            IF (NumChose .eq. 9 ) THEN
                CALL TestPLANE
              ENDIF
            IF (NumChose .eq. 10 ) THEN
                CALL TestARCCOSH
              ENDIF
            IF (NumChose .eq. 11 ) THEN
                CALL TestSINH
              ENDIF
            IF (NumChose .eq. 12 ) THEN
                CALL TestARCSINH
              ENDIF
            IF (NumChose .eq. 13 ) THEN
                CALL TestARCTANH
              ENDIF
            IF (NumChose .eq. 14 ) THEN
                CALL TestDOT
              ENDIF
            IF (NumChose .eq. 15 ) THEN
                CALL TestCROSS
              ENDIF
            IF (NumChose .eq. 16 ) THEN
                CALL TestMAG
              ENDIF
            IF (NumChose .eq. 17 ) THEN
                CALL TestNORM
              ENDIF
            IF (NumChose .eq. 18 ) THEN
                CALL TestROT1
              ENDIF
            IF (NumChose .eq. 19 ) THEN
                CALL TestROT2
              ENDIF
            IF (NumChose .eq. 20 ) THEN
                CALL TestROT3
              ENDIF
            IF (NumChose .eq. 21 ) THEN
                CALL TestADDVEC
              ENDIF
            IF (NumChose .eq. 22 ) THEN
                CALL TestADD3VEC
              ENDIF
            IF (NumChose .eq. 23 ) THEN
                CALL TestLNCOM1
              ENDIF
            IF (NumChose .eq. 24 ) THEN
                CALL TestLNCOM2
              ENDIF
            IF (NumChose .eq. 25 ) THEN
                CALL TestLNCOM3
              ENDIF
            IF (NumChose .eq. 26 ) THEN
                CALL TestPOLYFIT
              ENDIF
            IF (NumChose .eq. 27 ) THEN
                CALL TestANGLE
              ENDIF
            IF (NumChose .eq. 29 ) THEN
                CALL TestFACTOR
              ENDIF
            IF (NumChose .eq. 30 ) THEN
                CALL TestQUADRATIC
              ENDIF
            IF (NumChose .eq. 31 ) THEN
                CALL TestCUBIC
              ENDIF
            IF (NumChose .eq. 32 ) THEN
                CALL TestQUARTIC
              ENDIF
            IF (NumChose .eq. 33 ) THEN
                CALL TestMATSCALE
              ENDIF
            IF (NumChose .eq. 34 ) THEN
                CALL TestMATMULT
              ENDIF
            IF (NumChose .eq. 35 ) THEN
                CALL TestMATADD
              ENDIF
            IF (NumChose .eq. 36 ) THEN
                CALL TestMATSUB
              ENDIF
            IF (NumChose .eq. 37 ) THEN
                CALL TestMATTRANS
              ENDIF
            IF (NumChose .eq. 38 ) THEN
                CALL TestMAKEMAT
              ENDIF
            IF (NumChose .eq. 41 ) THEN
                CALL TestMATINVERSE
              ENDIF
            IF (NumChose .eq. 42 ) THEN
                CALL TestPRINTMAT
              ENDIF
            IF (NumChose .eq. 43 ) THEN
                CALL TestFILEPRINTMAT
              ENDIF
            IF (NumChose .eq. 44 ) THEN
                CALL TestFILEEXPPRINTMAT
              ENDIF
            IF (NumChose .eq. 45 ) THEN
                CALL TestDETERMINANT
              ENDIF

            ! -------------------- Tests for Time ----------------------
            IF (NumChose .eq. 100 ) THEN
                CALL TestUpCaseSt
              ENDIF
            IF (NumChose .eq. 101 ) THEN
                CALL TestINITTIME
              ENDIF
            IF (NumChose .eq. 102 ) THEN
                CALL TestGETINTMON
              ENDIF
            IF (NumChose .eq. 103 ) THEN
                CALL TestGETINTDAY
              ENDIF
            IF (NumChose .eq. 104 ) THEN
                CALL TestDAYOFWEEK
              ENDIF
            IF (NumChose .eq. 105 ) THEN
                CALL TestDAYLIGHTST
              ENDIF
            IF (NumChose .eq. 106 ) THEN
                CALL TestjDAY
              ENDIF
            IF (NumChose .eq. 107 ) THEN
                CALL TestjDAYALL
              ENDIF
            IF (NumChose .eq. 108 ) THEN
                CALL TestDAYS2MDHMS
              ENDIF
            IF (NumChose .eq. 109 ) THEN
                CALL TestINVjDAY
              ENDIF
            IF (NumChose .eq. 110 ) THEN
                CALL TestFINDDAYS
              ENDIF
            IF (NumChose .eq. 111 ) THEN
                CALL TestLSTIME
              ENDIF
            IF (NumChose .eq. 112 ) THEN
                CALL TestSUNRISESET
              ENDIF
            IF (NumChose .eq. 113 ) THEN
                CALL TestMOONRISESET
              ENDIF
            IF (NumChose .eq. 114 ) THEN
                CALL TestHMS_SEC
              ENDIF
            IF (NumChose .eq. 115 ) THEN
                CALL TestHMS_UT
              ENDIF
            IF (NumChose .eq. 116 ) THEN
                CALL TestHMS_RAD
              ENDIF
            IF (NumChose .eq. 117 ) THEN
                CALL TestDMS_RAD
              ENDIF
            ! -------------------- Tests for 2body ---------------------
            IF (NumChose .eq. 150 ) THEN
                CALL TestELORB
              ENDIF
            IF (NumChose .eq. 151 ) THEN
                CALL TestRANDV
              ENDIF
            IF (NumChose .eq. 152 ) THEN
                CALL TestFINDC2C3
              ENDIF
            IF (NumChose .eq. 153 ) THEN
                CALL TestNEWTONE
              ENDIF
            IF (NumChose .eq. 154 ) THEN
                CALL TestNEWTONM
              ENDIF
            IF (NumChose .eq. 155 ) THEN
                CALL TestNEWTONNU
              ENDIF
            IF (NumChose .eq. 156 ) THEN
                CALL TestKEPLER
              ENDIF
            IF (NumChose .eq. 157 ) THEN
                CALL TestFINDTOF
              ENDIF
            IF (NumChose .eq. 158 ) THEN
                CALL TestIJK2llA
              ENDIF
            IF (NumChose .eq. 159 ) THEN
                CALL TestGEOCGEOD
              ENDIF
            IF (NumChose .eq. 160 ) THEN
                CALL TestSIGHT
              ENDIF
            IF (NumChose .eq. 161 ) THEN
                CALL TestSUN
              ENDIF
            IF (NumChose .eq. 162 ) THEN
                CALL TestMOON
              ENDIF
            IF (NumChose .eq. 163 ) THEN
                CALL TestLIGHT
              ENDIF
            IF (NumChose .eq. 164 ) THEN
                CALL TestCHECKHITEARTH
              ENDIF
            IF (NumChose .eq. 165 ) THEN
                CALL TestSATFOV
              ENDIF
            IF (NumChose .eq. 166 ) THEN
                CALL TestRNGAZ
              ENDIF
            IF (NumChose .eq. 167 ) THEN
                CALL TestPATH
              ENDIF

            ! -------------------- Tests for iod -----------------------
            IF (NumChose .eq. 200 ) THEN
                CALL TestSITE
              ENDIF
            IF (NumChose .eq. 201 ) THEN
                CALL TestANGLESLAPLACE
              ENDIF
            IF (NumChose .eq. 202 ) THEN
                CALL TestANGLESGAUSS
              ENDIF
            IF (NumChose .eq. 203 ) THEN
                CALL TestRV_RADEC
              ENDIF
            IF (NumChose .eq. 204 ) THEN
                CALL TestRV_TRADEC
              ENDIF
            IF (NumChose .eq. 205 ) THEN
                CALL TestRV_RAZEL
              ENDIF
            IF (NumChose .eq. 206 ) THEN
                CALL TestRV_ELATLON
              ENDIF
            IF (NumChose .eq. 207 ) THEN
                CALL TestRVSEZ_RAZEL
              ENDIF
            IF (NumChose .eq. 208 ) THEN
                CALL TestRADEC_ELATLON
              ENDIF
            IF (NumChose .eq. 209 ) THEN
                CALL TestRADEC_AZEL
              ENDIF
            IF (NumChose .eq. 210 ) THEN
                CALL TestGIBBS
              ENDIF
            IF (NumChose .eq. 211 ) THEN
                CALL TestHERRGIBBS
              ENDIF
            IF (NumChose .eq. 212 ) THEN
                CALL TestLAMBERTUNIV
              ENDIF
            IF (NumChose .eq. 213 ) THEN
                CALL TestLAMBERTBATTIN
              ENDIF
            IF (NumChose .eq. 214 ) THEN
                CALL TestTARGET
              ENDIF

            ! -------------------- Tests for reduc ---------------------
            IF (NumChose .eq. 300 ) THEN
                CALL Testinitreduc
              ENDIF
            IF (NumChose .eq. 301 ) THEN
                CALL TestCONVTIME
              ENDIF
            IF (NumChose .eq. 302 ) THEN
                CALL TestPRECESSION
              ENDIF
            IF (NumChose .eq. 303 ) THEN
                CALL TestNUTATION
              ENDIF
            IF (NumChose .eq. 304 ) THEN
                CALL TestSIDEREAL
              ENDIF
            IF (NumChose .eq. 305 ) THEN
                CALL TestPOLARM
              ENDIF
            IF (NumChose .eq. 306 ) THEN
                CALL TestTrueMean
              ENDIF
            IF (NumChose .eq. 307 ) THEN
                CALL TestFK4
              ENDIF
            IF (NumChose .eq. 310 ) THEN
                CALL TestAllReduc
              ENDIF

            ! -------------------- Tests for manv ---------------------
            IF (NumChose .eq. 401 ) THEN
                CALL TestHohmann
              ENDIF
            IF (NumChose .eq. 402 ) THEN
                CALL TestBiElliptic
              ENDIF
            IF (NumChose .eq. 403 ) THEN
                CALL TestOneTangent
              ENDIF
            IF (NumChose .eq. 404 ) THEN
                CALL TestIOnlyChg
              ENDIF
            IF (NumChose .eq. 405 ) THEN
                CALL TestNodeOnlyChg
              ENDIF
            IF (NumChose .eq. 406 ) THEN
                CALL TestIandNodeChg
              ENDIF
            IF (NumChose .eq. 407 ) THEN
                CALL TestMinCombinedPlaneChg
              ENDIF
            IF (NumChose .eq. 408 ) THEN
                CALL TestCombinedPlaneChg
              ENDIF
            IF (NumChose .eq. 409 ) THEN
                CALL TestRendezvous
              ENDIF
            IF (NumChose .eq. 410 ) THEN
                CALL TestNonCoplanarRendz
              ENDIF
            IF (NumChose .eq. 411 ) THEN
                CALL TestHillsR
              ENDIF
            IF (NumChose .eq. 412 ) THEN
                CALL TestHillsV
              ENDIF
            IF (NumChose .eq. 413 ) THEN
                CALL TestIJK_RSW
              ENDIF
            IF (NumChose .eq. 414 ) THEN
                CALL TestCow2Hill
              ENDIF
            IF (NumChose .eq. 415 ) THEN
                CALL TestHill2Cow
              ENDIF

            ! -------------------- Tests for pert ---------------------
            IF (NumChose .eq. 501 ) THEN
                CALL TestPkepler
              ENDIF
            IF (NumChose .eq. 502 ) THEN
                CALL TestJ2DragPert
              ENDIF
            IF (NumChose .eq. 503 ) THEN
                CALL TestPredict
              ENDIF
            IF (NumChose .eq. 504 ) THEN
                CALL TestInitGravityField
              ENDIF
            IF (NumChose .eq. 505 ) THEN
                CALL TestLegPoly
              ENDIF
            IF (NumChose .eq. 506 ) THEN
                CALL TestDeriv
              ENDIF
            IF (NumChose .eq. 507 ) THEN
                CALL TestPertAccel
              ENDIF
            IF (NumChose .eq. 508 ) THEN
                CALL TestPDeriv
              ENDIF
            IF (NumChose .eq. 509 ) THEN
                CALL TestRK4
              ENDIF
            IF (NumChose .eq. 510 ) THEN
                CALL TestRKF45
              ENDIF
            IF (NumChose .eq. 511 ) THEN
                CALL TestCowell
              ENDIF
            IF (NumChose .eq. 512 ) THEN
                CALL TestAtmos
              ENDIF

            ! -------------------- Tests for DC -----------------------
            IF (NumChose .eq. 601 ) THEN
                CALL TestNonlin
              ENDIF
            IF (NumChose .eq. 602 ) THEN
                CALL TestFindAtwaAtwb
              ENDIF
            IF (NumChose .eq. 603 ) THEN
                CALL TestLeastSquares
              ENDIF
            IF (NumChose .eq. 604 ) THEN
*                CALL TestSequential
              ENDIF


           ENDDO

  999   CONTINUE

        CLOSE( 10 )
        CLOSE( 20 )
      STOP
      END
*
* ------------------------------------------------------------------------
*
*                       Individual SUBROUTINES
*
*  These subroutines implement each function or subroutine using a datafile
*    with the number of each routine to identify which subroutine to use.
*
*
* ------------------------------------------------------------------------

      SUBROUTINE TestGETPART
        IMPLICIT NONE
        EXTERNAL GetPart
        INTEGER GetPart, i, LocStart,Leng
        CHARACTER*250 InStr

        Read(10,*) InStr
        Write(20,*)  InStr
        Read(10,*) LocStart,Leng
        Write(20,*)  LocStart,Leng

        i = GETPART   ( InStr, LocStart,Leng )

        Write(20,*) '  Results:'
        Write(20,*) i

      RETURN
      END

      SUBROUTINE TestGETPARTL
        IMPLICIT NONE
        EXTERNAL GetPartL
        INTEGER LocStart,Leng
        INTEGER*4 GL,GetPartL
        CHARACTER*250 InStr

        Read(10,*) InStr
        Write(20,*)  InStr
        Read(10,*) LocStart,Leng
        Write(20,*)  LocStart,Leng

        Gl = GETPARTL( InStr, LocStart,Leng )

        Write(20,*) '  Results:'
        Write(20,*) Gl

      RETURN
      END

      SUBROUTINE TestGetPartR
        IMPLICIT NONE
        EXTERNAL GetPartR
        INTEGER LocStart,Leng
        REAL*8 GetPartR, GR
        CHARACTER*250 InStr

        Read(10,*) InStr
        Write(20,*)  InStr
        Read(10,*) LocStart,Leng
        Write(20,*)  LocStart,Leng

        GR = GetPartR   ( InStr, LocStart,Leng )

        Write(20,*) '  Results:'
        Write(20,*) GR

      RETURN
      END

      SUBROUTINE TestFACTORIAL
        IMPLICIT NONE
        EXTERNAL Factorial
        INTEGER x
        REAL*8 Factorial, Fact

        Read(10,*) x

        Fact = FACTORIAL( x )

        Write(20,*) '  Results:'
        Write(20,*) Fact

      RETURN
      END

      SUBROUTINE TestBINOMIAL
        IMPLICIT NONE
        EXTERNAL Binomial
        INTEGER i,j
        REAL*8 Binomial, Bio

        Read(10,*) i,j
        Write(20,*)  i,j

        Bio = BINOMIAL ( i,j )

        Write(20,*) '  Results:'
        Write(20,*) Bio

      RETURN
      END

      SUBROUTINE TestMIN
        IMPLICIT NONE
        EXTERNAL Min
        REAL*8 Min, Mn, x, y

        Read(10,*) x,y
        Write(20,*)  x,y

        Mn = MIN      ( X, Y )

        Write(20,*) '  Results:'
        Write(20,*) Mn

      RETURN
      END

      SUBROUTINE TestMAX
        IMPLICIT NONE
        EXTERNAL Max
        REAL*8 Max, Mx, x, y

        Read(10,*) x,y
        Write(20,*)  x,y

        Mx = MAX      ( X, Y )

        Write(20,*) '  Results:'
        Write(20,*) Mx

      RETURN
      END

      SUBROUTINE TestPLANE

        Read(10,*) x1,y1,z1
        Write(20,*)  x1,y1,z1
        Read(10,*) x2,y2,z2
        Write(20,*)  x2,y2,z2
        Read(10,*) x3,y3,z3
        Write(20,*)  x3,y3,z3

        CALL PLANE         ( x1,y1,z1,x2,y2,z2,x3,y3,z3, a,b,c,d )

        Write(20,*) '  Results:'
        Write(20,*) a,b,c,d

      RETURN
      END

      SUBROUTINE TestARCCOSH
        IMPLICIT NONE
        EXTERNAL ArcCosh
        REAL*8 ArcCosh, XVal, AC

        Read(10,*) XVal
        Write(20,*)  XVal

        AC = ARCCOSH( XVal )

        Write(20,*) '  Results:'
        Write(20,*) AC

      RETURN
      END

      SUBROUTINE TestSINH
        IMPLICIT NONE
        EXTERNAL Sinh
        REAL*8 Sinh, XVal, S

        Read(10,*) XVal
        Write(20,*)  XVal

        S = SINH   ( XVal )

        Write(20,*) '  Results:'
        Write(20,*) S

      RETURN
      END

      SUBROUTINE TestARCSINH
        IMPLICIT NONE
        EXTERNAL ArcSinh
        REAL*8 ArcSinh, XVal, AS

        Read(10,*) XVal
        Write(20,*)  XVal

        AS = ARCSINH( XVal )

        Write(20,*) '  Results:'
        Write(20,*) AS

      RETURN
      END

      SUBROUTINE TestARCTANH
        IMPLICIT NONE
        EXTERNAL ArcTanh
        REAL*8 ArcTanh, XVal, AT

        Read(10,*) XVal
        Write(20,*)  XVal

        AT = ARCTANH( XVal )

        Write(20,*) '  Results:'
        Write(20,*) AT

      RETURN
      END

      SUBROUTINE TestDOT
        IMPLICIT NONE
        EXTERNAL Dot
        REAL*8 Dot,Dt, Vec1(3),Vec2(3)

        Read(10,*) Vec1(1),Vec1(2),Vec1(3)
        Write(20,*)  Vec1(1),Vec1(2),Vec1(3)
        Read(10,*) Vec2(1),Vec2(2),Vec2(3)
        Write(20,*)  Vec2(1),Vec2(2),Vec2(3)

        Dt = DOT    ( Vec1,Vec2 )

        Write(20,*) '  Results:'
        Write(20,*) Dt

      RETURN
      END

      SUBROUTINE TestCROSS
        IMPLICIT NONE
        REAL*8 Vec1(3),Vec2(3), OutVec(3)

        Read(10,*) Vec1(1),Vec1(2),Vec1(3)
        Write(20,*)  Vec1(1),Vec1(2),Vec1(3)
        Read(10,*) Vec2(1),Vec2(2),Vec2(3)
        Write(20,*)  Vec2(1),Vec2(2),Vec2(3)

        CALL CROSS       ( Vec1,Vec2, OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestMAG
        IMPLICIT NONE
        REAL*8 Vec(3), magr,MAG
        EXTERNAL MAG

        Read(10,*) Vec(1),Vec(2),Vec(3)
        Write(20,*)  Vec(1),Vec(2),Vec(3)

        magr =  MAG ( Vec )
        Write(20,*) '  Results:'
        Write(20,*) magr

      RETURN
      END

      SUBROUTINE TestNORM
        IMPLICIT NONE
        REAL*8 Vec(3), OutVec(3)

        Read(10,*) Vec(1),Vec(2),Vec(3)
        Write(20,*)  Vec(1),Vec(2),Vec(3)

        CALL NORM        ( Vec, OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestROT1
        IMPLICIT NONE
        REAL*8 Vec(3), XVal, OutVec(3)

        Read(10,*) XVal
        Write(20,*)  XVal
        Read(10,*) Vec(1),Vec(2),Vec(3)
        Write(20,*)  Vec(1),Vec(2),Vec(3)

        CALL ROT1        ( Vec, XVal, OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestROT2
        IMPLICIT NONE
        REAL*8 Vec(3), XVal, OutVec(3)

        Read(10,*) XVal
        Write(20,*)  XVal
        Read(10,*) Vec(1),Vec(2),Vec(3)
        Write(20,*)  Vec(1),Vec(2),Vec(3)

        CALL ROT2        ( Vec, XVal, OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestROT3
        IMPLICIT NONE
        REAL*8 Vec(3), XVal, OutVec(3)

        Read(10,*) XVal
        Write(20,*)  XVal
        Read(10,*) Vec(1),Vec(2),Vec(3)
        Write(20,*)  Vec(1),Vec(2),Vec(3)

        CALL ROT3        ( Vec, XVal, OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestADDVEC
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), OutVec(3)

        Read(10,*) Vec1(1),Vec1(2),Vec1(3)
        Write(20,*)  Vec1(1),Vec1(2),Vec1(3)
        Read(10,*) Vec2(1),Vec2(2),Vec2(3)
        Write(20,*)  Vec2(1),Vec2(2),Vec2(3)

        CALL ADDVEC      ( Vec1,Vec2, OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestADD3VEC
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), Vec3(3), OutVec(3)

        Read(10,*) Vec1(1),Vec1(2),Vec1(3)
        Write(20,*)  Vec1(1),Vec1(2),Vec1(3)
        Read(10,*) Vec2(1),Vec2(2),Vec2(3)
        Write(20,*)  Vec2(1),Vec2(2),Vec2(3)
        Read(10,*) Vec3(1),Vec3(2),Vec3(3)
        Write(20,*)  Vec3(1),Vec3(2),Vec3(3)

        CALL ADD3VEC     ( Vec1,Vec2,Vec3, OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestLNCOM1
        IMPLICIT NONE
        REAL*8 Vec(3), A, OutVec(3)

        Read(10,*) a
        Write(20,*)  a
        Read(10,*) Vec(1),Vec(2),Vec(3)
        Write(20,*)  Vec(1),Vec(2),Vec(3)

        CALL LNCOM1      ( A, Vec,  OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestLNCOM2
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), A1, A2, OutVec(3)

        Read(10,*) A1,A2
        Write(20,*)  A1,A2
        Read(10,*) Vec1(1),Vec1(2),Vec1(3)
        Write(20,*)  Vec1(1),Vec1(2),Vec1(3)
        Read(10,*) Vec2(1),Vec2(2),Vec2(3)
        Write(20,*)  Vec2(1),Vec2(2),Vec2(3)

        CALL LNCOM2      ( A1, A2, Vec1, Vec2,  OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestLNCOM3
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), Vec3(3), A1, A2, A3, OutVec(3)

        Read(10,*) A1,A2, A3
        Write(20,*)  A1,A2, A3
        Read(10,*) Vec1(1),Vec1(2),Vec1(3)
        Write(20,*)  Vec1(1),Vec1(2),Vec1(3)
        Read(10,*) Vec2(1),Vec2(2),Vec2(3)
        Write(20,*)  Vec2(1),Vec2(2),Vec2(3)
        Read(10,*) Vec3(1),Vec3(2),Vec3(3)
        Write(20,*)  Vec3(1),Vec3(2),Vec3(3)

        CALL LNCOM3      ( A1, A2, A3, Vec1, Vec2, Vec3,  OutVec )

        Write(20,*) '  Results:'
        Write(20,*) OutVec(1),OutVec(2),OutVec(3)

      RETURN
      END

      SUBROUTINE TestPOLYFIT
        IMPLICIT NONE
        INTEGER Degree, NumPts, i
        REAL*8 DataPoints(50,2),Coeff(50,1), MinX, MinY

        Read(10,*) Degree, Numpts
        Write(20,*)  Degree, Numpts
        Do i = 1,NumPts
            Read(10,*) datapoints(i,1), DataPoints(i,2)
            Write(20,*)  datapoints(i,1), DataPoints(i,2)
          ENDDO

        CALL POLYFIT     ( Degree,NumPts,DataPoints,Coeff,
     &                         50, 50, MinX,MinY )

        Write(20,*) '  Results:'
        DO i=1,Degree
            Write(20,*) Coeff(i,1)
          ENDDO
        Write(20,*) Minx, MinY

      RETURN
      END

      SUBROUTINE TestANGLE
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), Theta
        EXTERNAL Dot

        Read(10,*) Vec1(1),Vec1(2),Vec1(3)
        Write(20,*)  Vec1(1),Vec1(2),Vec1(3)
        Read(10,*) Vec2(1),Vec2(2),Vec2(3)
        Write(20,*)  Vec2(1),Vec2(2),Vec2(3)

        CALL ANGLE       ( Vec1,Vec2, Theta )

        Write(20,*) '  Results:'
        Write(20,*) Theta

      RETURN
      END

      SUBROUTINE TestFACTOR
        IMPLICIT NONE
        REAL*8 Poly(16), Roots(15,2)
        INTEGER NRootS,i

        Read(10,*) NRoots
        Write(20,*)  NRoots
        DO i=1,NRoots+1
            Read(10,*) Poly(i)
            Write(20,*) Poly(i)
          ENDDO

        CALL FACTOR      ( Poly, NRootS, RootS )

        Write(20,*) '  Results:'
        Do i = 1,Nroots
            Write(20,*) Roots(i,1),Roots(i,2)
          ENDDO

      RETURN
      END

      SUBROUTINE TestQUADRATIC
        IMPLICIT NONE
        REAL*8 R1r,R1i,R2r,R2i,a,b,c

        Read(10,*) a,b,c
        Write(20,*)  a,b,c

        CALL QUADRATIC   ( a,b,c, R1r,R1i,R2r,R2i )

        Write(20,*) '  Results:'
        Write(20,*) R1r,R1i,R2r,R2i

      RETURN
      END

      SUBROUTINE TestCUBIC
        IMPLICIT NONE
        REAL*8 R1r,R1i,R2r,R2i,R3r,R3i,a,b,c,d

        Read(10,*) a,b,c,d
        Write(20,*)  a,b,c,d

        CALL CUBIC       ( a,b,c,d, R1r,R1i,R2r,R2i,R3r,R3i )

        Write(20,*) '  Results:'
        Write(20,*) R1r,R1i,R2r,R2i,R3r,R3i

      RETURN
      END

      SUBROUTINE TestQUARTIC
        IMPLICIT NONE
        REAL*8 R1r,R1i,R2r,R2i,R3r,R3i,R4r,R4i,a,b,c,d,e

        Read(10,*) a,b,c,d,e
        Write(20,*)  a,b,c,d,e

        CALL QUARTIC     ( a,b,c,d,e, R1r,R1i,R2r,R2i,
     &                         R3r,R3i,R4r,R4i )

        Write(20,*) '  Results:'
        Write(20,*) R1r,R1i,R2r,R2i,R3r,R3i,R4r,R4i

      RETURN
      END

      SUBROUTINE TestMATSCALE
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1(20,20),Scale,Mat2(20,20)
        CHARACTER*64 Title

        Read(10,*) Mat1r,Mat1c,Scale
        Write(20,*)  Mat1r,Mat1c,Scale
        CALL ReadMat( Mat1,Mat1r,Mat1c,20,20 )

        CALL MATSCALE    ( Mat1,Scale,Mat1r,Mat1c,20,20, Mat2 )

        Write(20,*) '  Results:'
        Title = ' Answer '
        CALL FilePrintMat( Mat2,Mat1r,Mat1c,20,20,Title )

      RETURN
      END

      SUBROUTINE TestMATMULT
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c,Mat2c
        REAL*8 Mat1(20,20),Mat2(20,20),Mat3(20,20)
        CHARACTER*64 Title

        Read(10,*) Mat1r,Mat1c,Mat2c
        Write(20,*)  Mat1r,Mat1c,Mat2c
        CALL ReadMat( Mat1,Mat1r,Mat1c,20,20 )
        CALL ReadMat( Mat2,Mat1c,Mat2c,20,20 )

        CALL MATMULT     ( Mat1,Mat2, Mat1r,Mat1c,Mat2c,20,20,20, Mat3)

        Write(20,*) '  Results:'
        Title = ' Answer '
        CALL FilePrintMat( Mat3,Mat1r,Mat2c,20,20,Title )

      RETURN
      END

      SUBROUTINE TestMATADD
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1(20,20),Mat2(20,20),Mat3(20,20)
        CHARACTER*64 Title

        Read(10,*) Mat1r,Mat1c
        Write(20,*)  Mat1r,Mat1c
        CALL ReadMat( Mat1,Mat1r,Mat1c,20,20 )
        CALL ReadMat( Mat2,Mat1r,Mat1c,20,20 )

        CALL MATADD      ( Mat1,Mat2, Mat1r,Mat1c,20,20, Mat3 )

        Write(20,*) '  Results:'
        Title = ' Answer '
        CALL FilePrintMat( Mat3,Mat1r,Mat1c,20,20,Title )

      RETURN
      END

      SUBROUTINE TestMATSUB
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1(20,20),Mat2(20,20),Mat3(20,20)
        CHARACTER*64 Title

        Read(10,*) Mat1r, Mat1c
        Write(20,*)  Mat1r, Mat1c
        CALL ReadMat( Mat1,Mat1r,Mat1c,20,20 )
        CALL ReadMat( Mat2,Mat1r,Mat1c,20,20 )

        CALL MATSUB      ( Mat1,Mat2, Mat1r,Mat1c,20,20, Mat3 )

        Write(20,*) '  Results:'
        Title = ' Answer '
        CALL FilePrintMat( Mat3,Mat1r,Mat1c,20,20,Title )

      RETURN
      END

      SUBROUTINE TestMATTRANS
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1(20,20),Mat2(20,20)
        CHARACTER*64 Title

        Read(10,*) Mat1r, Mat1c
        Write(20,*)  Mat1r, Mat1c
        CALL ReadMat( Mat1,Mat1r,Mat1c,20,20 )

        CALL MATTRANS    ( Mat1, Mat1r,Mat1c,20,20, Mat2 )

        Write(20,*) '  Results:'
        Title = ' Answer '
        CALL FilePrintMat( Mat2,Mat1c,Mat1r,20,20,Title )

      RETURN
      END

      SUBROUTINE TestMAKEMAT
        IMPLICIT NONE
        REAL*8 Angl, Matr(3,3)
        INTEGER Numbr
        CHARACTER*64 Title

        Read(10,*) Angl,Numbr
        Write(20,*) Angl,Numbr

        CALL MAKEMAT     ( Angl, Numbr, Matr )

        Write(20,*) '  Results:'
        Title = ' Answer '
        CALL FilePrintMat( Matr,3,3,3,3,Title )

      RETURN
      END

      SUBROUTINE TestMATINVERSE
        IMPLICIT NONE
        INTEGER Order
        REAL*8 Mat(20,20), MatInv(20,20)
        CHARACTER*64 Title

        Read(10,*) Order
        Write(20,*)  Order
        CALL ReadMat( Mat,Order,Order,20,20 )

        Title = ' test before '
       CALL filePRINTMAT    ( Mat,order,order,20,20,Title )
        CALL MATINVERSE  ( Mat, Order,20, MatInv )

        Write(20,*) '  Results:'
        Title = ' Answer '
        CALL FilePrintMat( MatInv,Order,Order,20,20,Title )

      RETURN
      END

      SUBROUTINE TestPRINTMAT
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1( 20,20  )
        CHARACTER*64 Title

        Read(10,*) Mat1r,Mat1c,Title
        Write(20,*)  Mat1r,Mat1c
        CALL ReadMat( Mat1,Mat1r,Mat1c,20,20 )

        CALL PRINTMAT    ( Mat1, Mat1r,Mat1c,20,20,Title )

        Write(20,*) '  Results:'
        Write(20,*)

      RETURN
      END

      SUBROUTINE TestFILEPRINTMAT
        IMPLICIT NONE
        INTEGER Mat1r,Mat1c
        REAL*8 Mat1( 20,20 )
        CHARACTER*64 Title

        Read(10,*) Mat1r,Mat1c,Title
        Write(20,*)  Mat1r,Mat1c
        CALL ReadMat( Mat1,Mat1r,Mat1c,20,20 )

        CALL FILEPRINTMAT( Mat1, Mat1r,Mat1c,20,20,Title )

        Write(20,*) '  Results:'
        Write(20,*)

      RETURN
      END

      SUBROUTINE TestFILEEXPPRINTMAT
        IMPLICIT NONE
        INTEGER Mat1r, Mat1c
        REAL*8 Mat1(20,20)
        CHARACTER*64 Title

        Read(10,*) Mat1r,Mat1c,Title
        Write(20,*)  Mat1r,Mat1c
        CALL ReadMat( Mat1,Mat1r,Mat1c,20,20 )

        CALL FILEEXPPRINTMAT ( Mat1,Mat1r,Mat1c,20,20,Title )

        Write(20,*) '  Results:'
        Write(20,*)

      RETURN
      END

      SUBROUTINE TestDETERMINANT
        IMPLICIT NONE
        INTEGER Order
        EXTERNAL Determinant
        REAL*8 Det, Determinant, Mat1( 20,20 )

        Read(10,*) Order
        Write(20,*)  Order
        CALL ReadMat( Mat1,Order,Order,20,20 )

        Det = DETERMINANT( Mat1, Order,20 )

        Write(20,*) '  Results:'
        Write(20,*) Det

      RETURN
      END

   ! -------------------- Time s -----------------------------------

      SUBROUTINE TestUpCaseSt
        IMPLICIT NONE
        CHARACTER*250 S
        EXTERNAL UpCaseSt
        CHARACTER*250 UpCaseSt, UpS

        Read(10,*) S
        Write(20,*)  S

        UpS = UpCaseSt ( S )

        Write(20,*) '  Results:'
        Write(20,*) Ups

      RETURN
      END


      SUBROUTINE TestINITTIME
        IMPLICIT NONE

        Read(10,*)
        Write(20,*) 'In '

        CALL INITTIME

        Write(20,*) '  Results:'
        Write(20,*)

      RETURN
      END

      SUBROUTINE TestGETINTMON
        IMPLICIT NONE
        CHARACTER*3 MonStr
        EXTERNAL GetIntMon
        INTEGER Month, GetintMon

        Read(10,*) monstr
        Write(20,*)  monstr

        Month = GETINTMON  ( MonStr )

        Write(20,*) '  Results:'
        Write(20,*) monstr,month

      RETURN
      END

      SUBROUTINE TestGETINTDAY
        IMPLICIT NONE
        CHARACTER*3 DayStr
        EXTERNAL GETINTDAY
        INTEGER GETINTDAY, DayVal

        Read(10,*) DayStr
        Write(20,*)  DayStr

        DayVal = GETINTDAY  ( DayStr )

        Write(20,*) '  Results:'
        Write(20,*) DayVal

      RETURN
      END

      SUBROUTINE TestDAYOFWEEK
        IMPLICIT NONE
        REAL*8 JD
        EXTERNAL DAYOFWEEK
        INTEGER DayW, DAYOFWEEK

        Read(10,*) JD
        Write(20,*)  JD

        DayW = DAYOFWEEK  ( JD )

        Write(20,*) '  Results:'
        Write(20,*) DayW

      RETURN
      END


      SUBROUTINE TestDAYLIGHTST
        IMPLICIT NONE
        INTEGER Year, StartDay, StopDay
        REAL*8  JDStartDST, JDStopDST

        Read(10,*) Year
        Write(20,*)  Year

        CALL DAYLIGHTST  ( Year, StartDay, StopDay, JDStartDST,
     &                         JDStopDST )

        Write(20,*) '  Results:'
        Write(20,*) StartDay, StopDay, JDStartDST, JDStopDST

      RETURN
      END

      SUBROUTINE TestjDAY
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, Min
        REAL*8  Sec, JD

        Read(10,*) Year, Mon, Day, Hr, Min, Sec
        Write(20,*)  Year, Mon, Day, Hr, Min, Sec

        CALL jDAY   ( Year,Mon,Day,Hr,Min, Sec, JD )

        Write(20,*) '  Results:'
        Write(20,*) JD

      RETURN
      END

      SUBROUTINE TestjDAYALL
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, Min
        CHARACTER WhichType
        REAL*8  Sec, JD

        Read(10,*) Year, Mon, Day, Hr, Min, Sec,  WhichType
        Write(20,*)  Year, Mon, Day, Hr, Min, Sec, WhichType

        CALL jDAYALL( Year,Mon,Day,Hr,Min,Sec, WhichType, JD )

        Write(20,*) '  Results:'
        Write(20,*) JD

      RETURN
      END

      SUBROUTINE TestDAYS2MDHMS
        IMPLICIT NONE
        REAL*8 Days,Sec
        INTEGER Year, Mon, Day, Hr, Min

        Read(10,*) Year, Days
        Write(20,*)  Year, Days

        CALL DAYS2MDHMS   ( Year,Days,  Mon,Day,Hr,Min,Sec )

        Write(20,*) '  Results:'
        Write(20,*) Mon,Day,Hr,Min,Sec

      RETURN
      END

      SUBROUTINE TestINVjDAY
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, Min
        REAL*8  Sec, JD

        Read(10,*) JD
        Write(20,*)  JD

        CALL INVjDAY( JD, Year,Mon,Day,Hr,Min, Sec )

        Write(20,*) '  Results:'
        Write(20,*) Year,Mon,Day,Hr,Min,Sec

      RETURN
      END

      SUBROUTINE TestFINDDAYS
        IMPLICIT NONE
        INTEGER Year, Month, Day, Hr, Min
        REAL*8 Sec, Days

        Read(10,*) year,Month,Day,Hr,Min,Sec
        Write(20,*)  year,Month,Day,Hr,Min,Sec

        CALL FINDDAYS    ( Year,Month,Day,Hr,Min, Sec,  Days )

        Write(20,*) '  Results:'
        Write(20,*) Days

      RETURN
      END

      SUBROUTINE TestLSTIME
        IMPLICIT NONE
        REAL*8 Lon, JD, LST, GST

        INCLUDE 'astmath.cmn'

        Read(10,*) Lon,JD
        Write(20,*)  Lon,JD
        Lon = Lon * deg2Rad

        CALL LSTIME      ( Lon, JD, LST, GST )

        Write(20,*) '  Results:'
        Write(20,*) LST*rad2deg,GST*rad2deg

      RETURN
      END

      SUBROUTINE TestSUNRISESET
        IMPLICIT NONE
        REAL*8 JD, Latgd, Lon, UTSunRise, UTSunSet
        CHARACTER WhichKind
        CHARACTER*12 Error
        INCLUDE 'astmath.cmn'

        Read(10,*) JD, Latgd, Lon, Whichkind
        Write(20,*)  JD, Latgd, Lon,  WhichKind
        Latgd= Latgd * deg2Rad
        Lon = Lon * deg2Rad

        CALL SUNRISESET  ( JD,Latgd,Lon, WhichKind, UTSunRise,
     &                         UTSunSet, Error )

        Write(20,*) '  Results:'
        Write(20,*) UTSunRise, UTSunSet, Error

      RETURN
      END

      SUBROUTINE TestMOONRISESET
        IMPLICIT NONE
        REAL*8 JD, Latgd, Lon, UTMoonRise, UTMoonSet, MoonPhaseAng
        CHARACTER*12 Error
        INCLUDE 'astmath.cmn'

        Read(10,*) JD, Latgd, Lon
        Write(20,*)  JD, Latgd, Lon
        Latgd= Latgd *deg2Rad
        Lon = Lon *deg2Rad

        CALL MOONRISESET ( JD,Latgd,Lon, UTMoonRise, UTMoonSet,
     &                         MoonPhaseAng, Error )

        Write(20,*) '  Results:'
        Write(20,*) UTMoonRise, UTMoonSet, MoonPhaseAng, Error

      RETURN
      END

      SUBROUTINE TestHMS_SEC
        IMPLICIT NONE
        INTEGER Hr, Min
        REAL*8  Sec, UTSec
        CHARACTER*4 Direction

        Read(10,*) Hr, Min, Sec, Direction
        Write(20,*)  Hr, Min, Sec,':', Direction,':'

        CALL HMS_SEC     ( Hr,Min, Sec, Direction, UTSec )

        Write(20,*) '  Results:'
        Write(20,*) UTSec

      RETURN
      END

      SUBROUTINE TestHMS_UT
        IMPLICIT NONE
        INTEGER Hr, Min
        REAL*8  Sec, UT
        CHARACTER*4 Direction

        Read(10,*) Hr, Min, Sec, Direction
        Write(20,*)  Hr, Min, Sec,  Direction

        CALL HMS_UT      ( Hr,Min, Sec, Direction, UT )

        Write(20,*) '  Results:'
        Write(20,*) UT

      RETURN
      END

      SUBROUTINE TestHMS_RAD
        IMPLICIT NONE
        INTEGER Hr, Min
        REAL*8  Sec, HMS
        CHARACTER*4 Direction

        Read(10,*) Hr, Min, Sec, Direction
        Write(20,*)  Hr, Min, Sec,  Direction

        CALL HMS_RAD     ( Hr,Min, Sec, Direction, HMS )

        Write(20,*) '  Results:'
        Write(20,*) HMS

      RETURN
      END

      SUBROUTINE TestDMS_RAD
        IMPLICIT NONE
        INTEGER Deg, Min
        REAL*8  Sec, DMS
        CHARACTER*4 Direction
        INCLUDE 'astmath.cmn'

        Read(10,*) Deg, Min, Sec,  Direction
        Write(20,*)  Deg, Min, Sec,  Direction
        Deg = Deg *deg2Rad

        CALL DMS_RAD     ( Deg,Min, Sec, Direction, DMS )

        Write(20,*) '  Results:'
        Write(20,*) DMS

      RETURN
      END



   ! --------------------  2 body s --------------------------------
      SUBROUTINE TestELORB
        IMPLICIT NONE
        REAL*8 R(3),V(3),P,A,Ecc,Incl,Omega,Argp,Nu,M,ArgLat,TrueLon,
     &         LonPer
        REAL*8 af,ag,meanlon,n,chi,psi
        REAL*8 ttt,jdut1,lod,xp,yp,latgc,lon
        REAL*8 rmag,vmag,rtasc,decl,fpa,fpav,az,rad
        REAL*8 rrac(3),vrac(3),rivc(3),vivc(3),tm(3,3)

        INTEGER Year, Mon, Day, Hr, MIN
        INTEGER TimeZone
        CHARACTER TypeUTIn
        CHARACTER*12  Error
        REAL*8 DUT1, DAT, UT1, TUT1, UTC, TAI, TT,
     &  Sec, JDTT, TDB, TTDB, JDTDB,DDpsi,DDeps

        INCLUDE 'astmath.cmn'

        rad = rad2deg

        Read(10,*) R(1),R(2),R(3)
        Read(10,*) V(1),V(2),V(3)
        Write(20,'(a,6(f14.7))') 'r', R(1),R(2),R(3),V(1),V(2),V(3)

       year = 2000
       mon  =  12
       day  =  15
       hr   =  16
       min  =  58
       sec  =  50.208
       dut1 =  0.10597
       dat  = 32
       xp   =  0.0
       yp   =  0.0
       lod  =  0.0
       timezone= 0
        Write(20,*) 'In',Year,Mon,Day,Hr,Min,Sec,TimeZone
        CALL CONVTIME    ( Year, Mon, Day, Hr, MIN, SEC,
     &                         TimeZone, 'C', DUT1, DAT, xp, yp,
     &                         UT1, TUT1, JDUT1, UTC, TAI, TT, TTT,
     &                         JDTT, TDB, TTDB, JDTDB, DDPSi, DDEps,
     &                         LOD, Error )

        CALL rv2coe( R,V,  P,A,Ecc,Incl,Omega,Argp,Nu,M,ArgLat,
     &              TrueLon,LonPer )

        Write(20,*) '  Results: '
        Write(20,'(3(f14.8))') P,A,Ecc
        Write(20,'(5(f14.7))') Incl*rad,Omega*rad,Argp*rad,Nu*rad,M*rad
        Write(20,'(3(f14.7))') ArgLat*rad, Truelon*rad, LonPer*rad


        CALL coe2rv(p,ecc,incl,omega,argp,nu,arglat,truelon,lonper,r,v)
        Write(20,'(a,6(f14.7))') 'r', R(1),R(2),R(3),V(1),V(2),V(3)

        ! ---- equinoctial elements
        CALL rv2eq ( r,v,  af,ag,meanlon,n,chi,psi )
        write(20,*)'           af           ag            meanlon    ',
     &            '       n           chi        psi'
        write(20,'(6f14.7)') af,ag,meanlon*rad,n,chi,psi
        CALL eq2rv( af, ag, meanlon, n, chi, psi, r, v)
        Write(20,'(a,6(f14.7))') 'r', R(1),R(2),R(3),V(1),V(2),V(3)

            ! ---- flight elements
        CALL rv2flt ( r,v,ttt,jdut1,lod,xp,yp,2, rmag,vmag,latgc,lon,
     &             fpa,az)
        write(20,*)'         rmag km       vmag km/s     latgc deg ',
     &            '      lon deg       fpa deg       az deg'
        write(20,'(6f14.7)') rmag,vmag,latgc*rad,lon*rad,fpa*rad,az*rad
        CALL flt2rv ( rmag,vmag,latgc,lon,fpa,az,ttt,jdut1,lod,xp,yp,
     &             2,r,v )
        Write(20,'(a,6(f14.7))') 'r', R(1),R(2),R(3),V(1),V(2),V(3)

            ! ----  adbarv elements
        CALL rv2adbar ( r,v, rmag,vmag,rtasc,decl,fpav,az )
        write(20,*)'          rmag km      vmag km/s     rtasc deg ',
     &            '      decl deg      fpav deg      az deg'
        write(20,'(6f14.7)') rmag,vmag,rtasc*rad,decl*rad,fpav*rad,
     &             az*rad
        CALL adbar2rv ( rmag,vmag,rtasc,decl,fpav,az, r,v )
        Write(20,'(a,6(f14.7))') 'r', R(1),R(2),R(3),V(1),V(2),V(3)

            ! ---- radial, along-track, cross-track
        CALL rv2rac(r,v, rrac, vrac, tm)
        Write(20,'(a,3(f14.7))') 'r', rRac(1),rRac(2),rRac(3)
        Write(20,'(3(f14.7))')  Vrac(1),Vrac(2),Vrac(3)

            ! ---- in-radial, velocity, cross-track
        CALL rv2ivc(r,v, rivc, vivc, tm)
        Write(20,'(a,3(f14.7))') 'r', Rivc(1),Rivc(2),Rivc(3)
        Write(20,'(3(f14.7))')  Vivc(1),Vivc(2),Vivc(3)

      RETURN
      END

      SUBROUTINE TestRANDV
        IMPLICIT NONE
        REAL*8 R(3),V(3),P,a,Ecc,Incl,Omega,Argp,Nu,ArgLat,TrueLon,
     &         LonPer

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        Read(10,*) p,a,ecc,incl,Omega,Argp,Nu,ArgLat,TrueLon,LonPer
        Write(20,*) p,a,ecc,incl,Omega,Argp,Nu,ArgLat,TrueLon,
     &        LonPer
        p = p * rekm
        incl= incl*deg2rad
        Omega= Omega*deg2rad
        Argp= Argp*deg2rad
        Nu= Nu*deg2rad
        ArgLat= ArgLat*deg2rad
        TrueLon= TrueLon*deg2rad
        LonPer= LonPer*deg2rad

        CALL coe2rv( P,Ecc,Incl,Omega,Argp,Nu,ArgLat,TrueLon,LonPer,R,V)

        Write(20,*) '  Results:'
        Write(20,'(4(f14.8))') R(1),R(2),R(3),R(3)
        Write(20,'(4(f14.8))') V(1),V(2),V(3),V(3)

      RETURN
      END

      SUBROUTINE TestFINDC2C3
        IMPLICIT NONE
        REAL*8 ZNew, C2New, C3New

        Read(10,*) ZNew
        Write(20,*)  ZNew

        CALL FINDC2C3    ( ZNew, C2New,C3New )

        Write(20,*) '  Results:'
        Write(20,*) C2New, C3New

      RETURN
      END

      SUBROUTINE TestNEWTONE
        IMPLICIT NONE
        REAL*8 Ecc, E0, M, Nu
        INCLUDE 'astmath.cmn'

        Read(10,*) Ecc, E0
        Write(20,*)  Ecc, E0
        E0 = E0 *deg2Rad

        CALL NEWTONE     ( Ecc,E0,  M,Nu )

        Write(20,*) '  Results:'
        Write(20,*) M*Rad2deg, Nu*Rad2deg

      RETURN
      END

      SUBROUTINE TestNEWTONM
        IMPLICIT NONE
        REAL*8 Ecc, M, E0, Nu

        INCLUDE 'astmath.cmn'
        Read(10,*) Ecc, M
        Write(20,*)  Ecc, M
        M = M *deg2rad

        CALL NEWTONM     ( Ecc,M, E0,Nu )

        Write(20,*) '  Results:'
        Write(20,*) E0*Rad2deg, Nu*Rad2deg

      RETURN
      END


      SUBROUTINE TestNEWTONNU
        IMPLICIT NONE
        REAL*8 Ecc, Nu, E0, M
        INCLUDE 'astmath.cmn'

        Read(10,*) Ecc, Nu
        Write(20,*)  Ecc, Nu
        Nu = Nu* Deg2Rad

        CALL NEWTONNU    ( Ecc,Nu, E0,M )

        Write(20,*) '  Results:'
        Write(20,*) E0*Rad2deg, M*Rad2deg

      RETURN
      END

      SUBROUTINE TestKEPLER
        IMPLICIT NONE
        REAL*8 Ro(3),Vo(3), dtsec,R(3),V(3)
        INTEGER i
        CHARACTER*12 Error
        REAL*8 vkmpersec

        INCLUDE 'astconst.cmn'

        VKmPerSec = 7.905366149846D0

        Read(10,*) Ro(1),Ro(2),Ro(3)
        Read(10,*) Vo(1),Vo(2),Vo(3)
        Read(10,*) dtsec

        DO i = 1,3
            ro(i) = ro(i) * rekm
            vo(i) = vo(i) * VKmPerSec
          ENDDO

        Write(20,*)  Ro(1),Ro(2),Ro(3)
        Write(20,*)  Vo(1),Vo(2),Vo(3)
        Write(20,*)  dtsec

        CALL KEPLER      ( Ro,Vo, dtsec, R,V, Error )

        Write(20,*) '  Results:'
        Write(20,*) R(1),R(2),R(3)
        Write(20,*) V(1),V(2),V(3)

      RETURN
      END

      SUBROUTINE TestFINDTOF
        IMPLICIT NONE
        REAL*8 Ro(3), R(3), p, Tof
        INTEGER i

        INCLUDE 'astconst.cmn'

        Read(10,*) Ro(1),Ro(2),Ro(3)
        Read(10,*) R(1),R(2),R(3), p

        p = p * rekm
        DO i = 1,3
            ro(i) = ro(i) * rekm
            r(i) = r(i) * rekm
          ENDDO

        Write(20,*)  Ro(1),Ro(2),Ro(3)
        Write(20,*)  R(1),R(2),R(3), p

        CALL FINDTOF     ( Ro,R, p, Tof )

        Write(20,*) '  Results:'
        Write(20,*) tof

      RETURN
      END

      SUBROUTINE TestIJK2llA
        IMPLICIT NONE
        REAL*8 R(3), JD, Latgc,Latgd,Lon,Hellp
        INTEGER i
        EXTERNAL GSTIME
        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        Read(10,*) r(1),r(2),r(3)
        Read(10,*) JD

        DO i = 1,3
            r(i) = r(i) * rekm
          ENDDO
        Write(20,*)  r(1),r(2),r(3)

        CALL IJK2llA( R, JD, Latgc,Latgd,Lon,Hellp )

        Write(20,*) '  Results:'
        Write(20,*) LatGc*rad2deg,LatGd*Rad2deg,Lon*rad2deg,Hellp

      RETURN
      END

      SUBROUTINE TestGEOCGEOD
        IMPLICIT NONE
        REAL*8 Latgc,Latgd
        CHARACTER*4 Direction

        INCLUDE 'astmath.cmn'

        Read(10,*) Latgc, Direction
        Write(20,*)  Latgc, Direction
        Latgc = Latgc * deg2rad

        CALL gc2gd    ( Latgc, Direction, Latgd )

        Write(20,*) '  Results:'
        Write(20,*) Latgd*rad2deg

      RETURN
      END

      SUBROUTINE TestSIGHT
        IMPLICIT NONE
        REAL*8 R1(3),R2(3)
        INTEGER i
        CHARACTER WhichKind
        CHARACTER*3 LOS

        INCLUDE 'astconst.cmn'

        Read(10,*) r1(1),r1(2),r1(3)
        Read(10,*) r2(1),r2(2),r2(3),  WhichKind

        DO i = 1,3
            r1(i) = r1(i) * rekm
            r2(i) = r2(i) * rekm
          ENDDO

        Write(20,*)  r1(1),r1(2),r1(3)
        Write(20,*)  r2(1),r2(2),r2(3),  WhichKind

        CALL SIGHT       ( R1,R2, WhichKind, LOS )

        Write(20,*) '  Results:'
        Write(20,*) LOS

      RETURN
      END

      SUBROUTINE TestSUN
        IMPLICIT NONE
        REAL*8 JD,RSun(3), RtAsc, Decl

        Read(10,*) JD
        Write(20,*)  JD

        CALL SUN         ( JD, RSun, RtAsc,Decl )

        Write(20,*) '  Results:'
        Write(20,*) RSun(1),RSun(2),RSun(3), RSun(3)
        Write(20,*) RtAsc, Decl

      RETURN
      END

      SUBROUTINE TestMOON
        IMPLICIT NONE
        REAL*8 JD, RMoon(3), RtAsc, Decl

        Read(10,*) JD
        Write(20,*)  JD

        CALL MOON        ( JD, RMoon, RtAsc,Decl )

        Write(20,*) '  Results:'
        Write(20,*) RMoon(1),RMoon(2),RMoon(3), RMoon(3)
        Write(20,*) RtAsc, Decl

      RETURN
      END

      SUBROUTINE TestLIGHT
        IMPLICIT NONE
        REAL*8 R(3), JD
        CHARACTER  WhichKind, Lit(3)
        INTEGER i

        INCLUDE 'astconst.cmn'

        Read(10,*) r(1),r(2),r(3)
        Read(10,*) JD,  WhichKind

        DO i = 1,3
            r(i) = r(i) * rekm
          ENDDO
        Write(20,*)  r(1),r(2),r(3)
        Write(20,*)  JD,  WhichKind

        CALL LIGHT       ( R, JD, WhichKind, LIT )

        Write(20,*) '  Results:'
        Write(20,*) LIT

      RETURN
      END

      SUBROUTINE TestCHECKHITEARTH
        IMPLICIT NONE
        REAL*8 RInt(3), V1t(3), RTgt(3),V2t(3)
        INTEGER i
        CHARACTER HitEarth

        INCLUDE 'astconst.cmn'

        Read(10,*) Rint(1)
        Read(10,*) Rtgt(1)

        DO i = 1,3
            rint(i) = rint(i) * rekm
            rtgt(i) = rtgt(i) * rekm
          ENDDO

        Write(20,*)  Rint(1)

        CALL CHECKHITEARTH ( Rint,V1t,Rtgt,V2t, HitEarth )

        Write(20,*) '  Results:'
        Write(20,*) HitEarth

      RETURN
      END

      SUBROUTINE TestSATFOV
        IMPLICIT NONE
        REAL*8 Incl, Az, SLatgd, SLon, SAlt, tFOV, EtaCtr, FovMAx,
     &     TotalRng, RhoMax, RhoMin,TgtLat, TgtLon
        INCLUDE 'astmath.cmn'

        Read(10,*) Incl,Az,SLatgd,SLon,SAlt,tFOV,EtaCtr,FOVMax
        Write(20,*)  Incl,Az,SLatgd,SLon,SAlt,tFOV,EtaCtr,FOVMax
        Incl = Incl *deg2rad
        az = az *deg2rad
        slatgd=slatgd *deg2rad
        slon=slon*deg2rad

        CALL SATFOV      ( Incl,Az, SLatgd, SLon, SAlt,tFOV,EtaCtr,
     &                       FovMax, TotalRng, RhoMax, RhoMin,TgtLat,
     &                       TgtLon )

        Write(20,*) '  Results:'
        Write(20,*) TotalRng, RhoMax, RhoMin, TgtLat,TgtLon

      RETURN
      END

      SUBROUTINE TestRNGAZ
        IMPLICIT NONE
        REAL*8 LLat, LLon, TLat, TLon, Tof, Range, Az
        INCLUDE 'astmath.cmn'

        Read(10,*) LLat,LLon,TLat,TLon
        Write(20,*)  LLat,LLon,TLat,TLon
        LLat = LLat *deg2rad
        LLon = LLon *deg2rad
        TLat = TLat *deg2rad
        TLon = TLon *deg2rad

        CALL RNGAZ       ( LLat,LLon,TLat,TLon,Tof, Range, Az )

        Write(20,*) '  Results:'
        Write(20,*) TOF,Range,Az*Rad2deg


      RETURN
      END

      SUBROUTINE TestPATH
        IMPLICIT NONE
        REAL*8 LLat, LLon, Range, Az, TLat, TLon

        INCLUDE 'astmath.cmn'

        Read(10,*) LLat,LLon,Range,Az
        Write(20,*)  LLat,LLon,Range,Az
        LLat = LLat *deg2rad
        LLon = LLon *deg2rad
        Az = Az *deg2rad

        CALL PATH        ( LLat, LLon, Range, Az, TLat, TLon )

        Write(20,*) '  Results:'
        Write(20,*) TLat*rad2deg,TLon*rad2deg

      RETURN
      END

   ! --------------------  IOD s -----------------------------------
      SUBROUTINE TestSITE
        IMPLICIT NONE
        REAL*8 LatGd,Alt,Lon,RS(3),VS(3), ReFt
        INCLUDE 'astmath.cmn'

        Reft = 20925644.0288713D0

        Read(10,*) Latgd,Alt,Lon
        Write(20,*)  Latgd,Alt,Lon
        Latgd= Latgd *deg2Rad
        Lon= Lon *deg2Rad
        Alt= Alt / REFt

        CALL  SITE        ( Latgd,Alt,Lon, RS,VS )

        Write(20,*) '  Results:'
        Write(20,*) RS(1),RS(2),RS(3)
        Write(20,*) VS(1),VS(2),VS(3)

        Latgd= (-7.0D0 - 54.0D0/60.0D0 - 23.886D0/3600.0D0) *deg2Rad
        Lon= (345.0D0 + 35.0D0/60.0D0 + 51.0D0/3600.0D0) *deg2Rad
        Alt= 0.056D0

        CALL  SITE        ( Latgd,Alt,Lon, RS,VS )

        Write(20,*) '  Results:'
        Write(20,*) RS(1),RS(2),RS(3)
        Write(20,*) VS(1),VS(2),VS(3)

      RETURN
      END

      SUBROUTINE TestANGLESLAPLACE
        IMPLICIT NONE
        REAL*8 Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,JD2,JD3,
     &         RS1(3),RS2(3),RS3(3),r2(3),v2(3)


        Read(10,*) Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,JD2,JD3
        Write(20,*)  Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,
     &         JD2,JD3

        CALL  ANGLESLAPLACE( Delta1,Delta2,Delta3,Alpha1,Alpha2,
     &                          Alpha3,JD1,JD2,JD3,RS1,RS2,RS3, r2,v2)

        Write(20,*) '  Results:'
        Write(20,*) R2(1),R2(2),R2(3)
        Write(20,*) V2(1),V2(2),V2(3)

      RETURN
      END

      SUBROUTINE TestANGLESGAUSS
        IMPLICIT NONE
        REAL*8 Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,JD2,JD3,
     &         RS1(3),RS2(3),RS3(3),r2(3),v2(3)


        Read(10,*) Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,JD2,JD3
        Write(20,*)  Delta1,Delta2,Delta3,Alpha1,Alpha2,Alpha3,JD1,
     &      JD2,JD3

        CALL  ANGLESGAUSS ( Delta1,Delta2,Delta3,Alpha1,Alpha2,
     &                         Alpha3,JD1,JD2,JD3,RS1,RS2,RS3, r2,v2 )

        Write(20,*) '  Results:'
        Write(20,*) R2(1),R2(2),R2(3)
        Write(20,*) V2(1),V2(2),V2(3)

      RETURN
      END

      SUBROUTINE TestRV_RADEC
        IMPLICIT NONE
        REAL*8 Rijk(3),Vijk(3),rr,RtAsc,Decl,DRr,DRtAsc,DDecl,
     &       rekm,vkmpersec
        CHARACTER*4 Direction
        INTEGER i

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        VKmPerSec = 7.905366149846D0

        Read(10,*) rijk(1),rijk(2),rijk(3)
        Read(10,*) vijk(1),vijk(2),vijk(3), Direction
        Write(20,*) vijk(1),vijk(2),vijk(3), Direction

        CALL  RV_RADEC    ( Rijk,Vijk, Direction, rr,RtAsc,Decl,
     &                         DRr,DRtAsc,DDecl )

        Write(20,*) '  Results:'
        Write(20,*) rr,rtasc*rad2deg,decl*rad2deg,drr,drtasc,ddecl

        Direction = 'FROM'
        CALL  RV_RADEC    ( Rijk,Vijk, Direction, rr,RtAsc,Decl,
     &                         DRr,DRtAsc,DDecl )
        Write(20,*) '  Results:'
        Write(20,*) Rijk(1),Rijk(2),Rijk(3)
        Write(20,*) Vijk(1),Vijk(2),Vijk(3)

      RETURN
      END

      SUBROUTINE TestRV_TRADEC
        IMPLICIT NONE
        REAL*8 Rijk(3),VIjk(3),RS(3),Rho,TRtAsc,TDecl,DRho,DTRtAsc,
     &         DTDecl,rekm,vkmpersec
        CHARACTER*4 Direction
        INTEGER i
        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        VKmPerSec = 7.905366149846D0

        Read(10,*) rijk(1),rijk(2),rijk(3)
        Read(10,*) vijk(1),vijk(2),vijk(3)
        Read(10,*) rs(1),rs(2),rs(3), Direction
        Write(20,*) rs(1),rs(2),rs(3), Direction

        CALL  RV_TRADEC   ( Rijk,Vijk,RS, Direction, Rho,TRtAsc,
     &                         TDecl,DRho,DTRtAsc,DTDecl )

        Write(20,*) '  Results:'
        Write(20,*) Rho,TRtasc*rad2deg,TDecl*rad2deg,DRho,DTRtasc,DTDecl

        Direction = 'FROM'
        CALL  RV_TRADEC   ( Rijk,Vijk,RS, Direction, Rho,TRtAsc,
     &                         TDecl,DRho,DTRtAsc,DTDecl )

        Write(20,*) '  Results:'
        Write(20,*) Rijk(1),Rijk(2),Rijk(3)
        Write(20,*) Vijk(1),Vijk(2),Vijk(3)
        Write(20,*) Rs(1),Rs(2),Rs(3)

      RETURN
      END

      SUBROUTINE TestRV_RAZEL
        IMPLICIT NONE
        REAL*8 Reci(3),Veci(3),Latgd,Lon,alt,Rho,Az,El,DRho,DAz,DEl
        INTEGER terms
        CHARACTER*4 Direction
        INTEGER Year, Mon, Day, Hr, MIN,terms
        CHARACTER*12  Error
        REAL*8 DUT1, DAT, UT1, TUT1, UTC, TAI, TT,lod,xp,yp,
     &  Sec, JDTT, TDB, TTDB, JDTDB,DDpsi,DDeps,jdut1,ttt

        INCLUDE 'astmath.cmn'

        Read(10,*) Reci(1),Reci(2),Reci(3)
        Write(20,*) Reci(1),Reci(2),Reci(3)
        Read(10,*) Veci(1),Veci(2),Veci(3)
        Write(20,*)  Veci(1),Veci(2),Veci(3)
        Read(10,*)Latgd,Lon,alt, Direction
        Write(20,*) Latgd,Lon,alt, Direction
        Read(10,*) year,mon,day,hr,min,sec
        Write(20,*) year,mon,day,hr,min,sec
        Read(10,*) dut1,dat,lod,xp,yp,terms
        Write(20,*) dut1,dat,lod,xp,yp,terms
        Latgd= Latgd *deg2Rad
        Lon= Lon *deg2Rad

        CALL CONVTIME    ( Year, Mon, Day, Hr, MIN, SEC,
     &                         0, 'C', DUT1, DAT, xp, yp,
     &                         UT1, TUT1, JDUT1, UTC, TAI, TT, TTT,
     &                         JDTT, TDB, TTDB, JDTDB, DDPSi, DDEps,
     &                         LOD, Error )

        CALL RV_RAZEL    ( Reci,Veci,Latgd,Lon,alt,TTT,jdut1,lod,
     &                     xp,yp,terms, Direction,
     &                     Rho,Az,El,DRho,DAz,DEl )

        Write(20,*) '  Results:'
        Write(20,*) Rho,Az*rad2deg,El*rad2deg,DRho,DAz,Del

        Direction = 'FROM'
        CALL RV_RAZEL    ( Reci,Veci,Latgd,Lon,alt,TTT,jdut1,lod,
     &                     xp,yp,terms, Direction,
     &                     Rho,Az,El,DRho,DAz,DEl )

        Write(20,*) '  Results:'
        Write(20,*) Reci(1),Reci(2),Reci(3)
        Write(20,*) Veci(1),Veci(2),Veci(3)

      RETURN
      END

      SUBROUTINE TestRV_ELATLON
        IMPLICIT NONE
        REAL*8 Rijk(3), Vijk(3), rr,EclLat,EclLon,DRr,DEclLat,DEclLon,
     &         rekm,VKmPerSec
        CHARACTER*4 Direction
        INTEGER i
        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        VKmPerSec = 7.905366149846D0


        Read(10,*) Rijk(1),Rijk(2),Rijk(3)
        Write(20,*)  Rijk(1),Rijk(2),Rijk(3)
        Read(10,*) Vijk(1),Vijk(2),Vijk(3),Direction
        Write(20,*)  Vijk(1),Vijk(2),Vijk(3),Direction
        DO i = 1,3
            rijk(i) = rijk(i) / rekm
            vijk(i) = vijk(i) / VKmPerSec
          ENDDO

        CALL  RV_ELATLON  ( Rijk,Vijk, Direction, rr,EclLat,EclLon,
     &                         DRr,DEclLat,DEclLon )

        Write(20,*) '  Results:'
        Write(20,*) rr,EclLat*rad2deg,EclLon*Rad2deg,Drr,DEclLat,DEclLon

        Direction = 'FROM'
        CALL  RV_ELATLON  ( Rijk,Vijk, Direction, rr,EclLat,EclLon,
     &                         DRr,DEclLat,DEclLon )

        Write(20,*) '  Results:'
        Write(20,*) Rijk(1),Rijk(2),Rijk(3)
        Write(20,*) Vijk(1),Vijk(2),Vijk(3)

      RETURN
      END

      SUBROUTINE TestRVSEZ_RAZEL
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 RhoSez(3), DRhoSez(3),Rho,Az,El,DRho,DAz, DEl
        INCLUDE 'astmath.cmn'

        Read(10,*) Rhosez(1),Rhosez(2),Rhosez(3)
        Write(20,*)  Rhosez(1),Rhosez(2),Rhosez(3)
        Read(10,*) DRhosez(1),DRhosez(2),DRhosez(3),Direction
        Write(20,*)  DRhosez(1),DRhosez(2),DRhosez(3),Direction

        CALL  RVSEZ_RAZEL ( Rhosez,DRhosez,Direction, Rho,Az,El,
     &                         DRho,DAz,DEl )

        Write(20,*) '  Results:'
        Write(20,*) Rho, Az*rad2deg, El*rad2deg, DRho, Daz, Del

        Direction = 'FROM'
        CALL  RVSEZ_RAZEL ( Rhosez,DRhosez,Direction, Rho,Az,El,
     &                         DRho,DAz,DEl )

        Write(20,*) '  Results:'
        Write(20,*) Rhosez(1),Rhosez(2),Rhosez(3)
        Write(20,*) DRhosez(1),DRhosez(2),DRhosez(3)

      RETURN
      END

      SUBROUTINE TestRADEC_ELATLON
        IMPLICIT NONE
        REAL*8 RtAsc,Decl,EclLat,EclLon
        CHARACTER*4 Direction
        INCLUDE 'astmath.cmn'

        Read(10,*) RtAsc,Decl,Direction
        Write(20,*)  RtAsc,Decl,Direction
        RtAsc = RtAsc *deg2rad
        Decl = Decl *deg2rad

        CALL  RADEC_ELATLON ( RtAsc,Decl,Direction, EclLat, EclLon)

        Write(20,*) '  Results:'
        Write(20,*) EclLat*rad2deg, EclLon*rad2deg

        Direction = 'FROM'
        CALL  RADEC_ELATLON ( RtAsc,Decl,Direction, EclLat, EclLon)

        Write(20,*) '  Results:'
        Write(20,*) RtAsc*rad2deg, Decl*rad2deg

      RETURN
      END

      SUBROUTINE TestRADEC_AZEL
        IMPLICIT NONE
        REAL*8 RtAsc,Decl,LST,Latgc,Az,El
        CHARACTER*4 Direction
        INCLUDE 'astmath.cmn'

        Read(10,*) RtAsc,Decl,LST,LatGc,Direction
        Write(20,*)  RtAsc,Decl,LST,LatGc,Direction

        CALL  RADEC_AZEL  ( RtAsc,Decl,LST,Latgc, Direction, Az,El)

        Write(20,*) '  Results:'
        Write(20,*) Az*rad2deg,El*rad2deg

      RETURN
      END

      SUBROUTINE TestGIBBS
        IMPLICIT NONE
        REAL*8 R1(3), R2(3), R3(3), V2(3), Theta, Theta1, Copa
        CHARACTER*12 Error
        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        Read(10,*) R1(1),R1(2),R1(3)
        Read(10,*) R2(1),R2(2),R2(3)
        Read(10,*) R3(1),R3(2),R3(3)

        CALL lncom1(rekm,r1,r1)
        CALL lncom1(rekm,r2,r2)
        CALL lncom1(rekm,r3,r3)

        Write(20,*)  R1(1),R1(2),R1(3)
        Write(20,*)  R2(1),R2(2),R2(3)
        Write(20,*)  R3(1),R3(2),R3(3)

        CALL  GIBBS       ( R1,R2,R3, V2, Theta,Theta1,Copa, Error)

        Write(20,*) '  Results: '
        Write(20,*) V2(1),V2(2),V2(3)
        Write(20,*) 'angle ',theta*rad2deg,Theta1*rad2deg,
     &              Copa*rad2deg,Error

      RETURN
      END

      SUBROUTINE TestHERRGIBBS
        IMPLICIT NONE
        REAL*8 R1(3), R2(3), R3(3), JD1, JD2, JD3, V2(3), Theta,
     &         Theta1, Copa,t1,t2,t3
        CHARACTER*12 Error
        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        Read(10,*) R1(1),R1(2),R1(3),t1
        JD1 = 2448608.0D0 + t1/86400.0D0
        Read(10,*) R2(1),R2(2),R2(3),t2
        JD2 = 2448608.0D0 + t2/86400.0D0
        Read(10,*) R3(1),R3(2),R3(3),t3
        JD3 = 2448608.0D0 + t3/86400.0D0

        CALL lncom1(rekm,r1,r1)
        CALL lncom1(rekm,r2,r2)
        CALL lncom1(rekm,r3,r3)

        Write(20,*)  R1(1),R1(2),R1(3),JD1
        Write(20,*)  R2(1),R2(2),R2(3),JD2
        Write(20,*)  R3(1),R3(2),R3(3),JD3

        CALL  HERRGIBBS   ( R1,R2,R3,JD1,JD2,JD3, V2, Theta,Theta1,
     &                         Copa, Error )

        Write(20,*) '  Results: '
        Write(20,*) V2(1),V2(2),V2(3)
        Write(20,*) 'angle ',theta*rad2deg,Theta1*rad2deg,
     &              Copa*rad2deg,Error

      RETURN
      END

      SUBROUTINE TestLAMBERTUNIV
        IMPLICIT NONE
        REAL*8 Ro(3), r(3), dtsec, Vo(3),v(3)
        INTEGER i
        CHARACTER Dm, OverRev
        CHARACTER*12 Error

        INCLUDE 'astconst.cmn'

        Read(10,*) Ro(1),Ro(2),Ro(3)
        Read(10,*) R(1),R(2),R(3)
        Read(10,*) dtsec,dm,OverRev
        DO i = 1,3
            ro(i) = ro(i) * rekm
            r(i) = r(i) * rekm
          ENDDO
        dtsec = dtsec * 806.810991D0

        Write(20,*)  Ro(1),Ro(2),Ro(3)
        Write(20,*)  R(1),R(2),R(3)
        Write(20,'(5X,A1,A1,F16.6)') dm,OverRev,dtsec

        Write(20,*)  dm,' : ',overrev,' ',dtsec

        CALL  LAMBERTUNIV ( ro,r, dm,OverRev, dtsec, vo,v, Error )

        Write(20,*) '  Results:'
        Write(20,*) Vo(1),Vo(2),Vo(3)
        Write(20,*) V(1),V(2),V(3)

      RETURN
      END

      SUBROUTINE TestLAMBERTBATTIN
        IMPLICIT NONE
        REAL*8 Ro(3), r(3), dtsec, Vo(3),v(3)
        INTEGER i
        CHARACTER Dm, OverRev
        CHARACTER*12 Error

        INCLUDE 'astconst.cmn'

        Read(10,*) Ro(1),Ro(2),Ro(3)
        Read(10,*) R(1),R(2),R(3)
        Read(10,*) dtsec, dm,OverRev

        DO i = 1,3
            ro(i) = ro(i) * rekm
            r(i) = r(i) * rekm
          ENDDO
        dtsec = dtsec * 806.810991D0

        Write(20,*)  dm,OverRev,dtsec
        Write(20,*)  Ro(1),Ro(2),Ro(3)
        Write(20,*)  R(1),R(2),R(3)


        CALL  LAMBERTBATTIN ( ro,r, dm,OverRev, dtsec, vo,v, Error )

        Write(20,*) '  Results:'
        Write(20,*) Vo(1),Vo(2),Vo(3)
        Write(20,*) V(1),V(2),V(3)

      RETURN
      END

      SUBROUTINE TestTARGET
        IMPLICIT NONE
        REAL*8 RInt(3),VInt(3),RTgt(3),VTgt(3),dtsec,V1t(3),V2t(3),
     &         DV1(3),DV2(3)
        CHARACTER Kind, Dm
        CHARACTER*12 Error

        Read(10,*)
        Write(20,*) 'In '

        CALL  TARGET      ( RInt,VInt,RTgt,VTgt, Dm,Kind, dtsec,
     &                         V1t,V2t,DV1,DV2, Error  )

        Write(20,*) '  Results:'
        Write(20,*)

      RETURN
      END

   ! --------------------  REDUC's -----------------------------------

      SUBROUTINE Testinitreduc
       IMPLICIT NONE
        CHARACTER*64 FileN1
c       INTEGER IAr00(5,106), IAR80(5,106)
c       REAL*8  RAr00(4,106), RAr80(4,106)
c       INTEGER IAr96(5,263)
c       REAL*8  RAr96(4,263)
c       INTEGER pIAr96(10,112)
c       REAL*8  pRAr96(4,112)

        INCLUDE 'astreduc.cmn'

        Read(10,*) FileN1
c        Read(10,*) FileN2
        Write(20,*) 'In '

c        CALL initreduc ( FileN1,FileN2,IAr00,RAr00, IAr80,RAr80,
c     &                      IAr96,RAr96, pIAr96,pRAr96 )
        CALL initreduc ( FileN1 )

        Write(20,*) '  Results:'
        Write(20,*) IAr80(1,1),' ',RAr80(1,1)*3600.0D0/0.0001D0
        Write(20,*) IAr80(2,1),' ',RAr80(2,1)*3600.0D0/0.0001D0
        Write(20,*) IAr80(3,1),' ',RAr80(3,1)*3600.0D0/0.0001D0
      RETURN
      END

      SUBROUTINE TestConvTime
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, MIN
        INTEGER TimeZone
        CHARACTER TypeUTIn
        CHARACTER*12  Error
        REAL*8 DUT1, DAT, xp, yp, UT1, TUT1, JDUT1, UTC, TAI, TT, TTT,
     &  Sec, JDTT, TDB, TTDB, JDTDB,DDpsi,DDeps,lod

        Read(10,*) Year,Mon,Day,Hr,Min,Sec,TimeZone,TypeUtIn
        Write(20,*) 'In',Year,Mon,Day,Hr,Min,Sec,TimeZone,TypeUtIn

        CALL CONVTIME    ( Year, Mon, Day, Hr, MIN, SEC,
     &                         TimeZone, TypeUTIn, DUT1, DAT, xp, yp,
     &                         UT1, TUT1, JDUT1, UTC, TAI, TT, TTT,
     &                         JDTT, TDB, TTDB, JDTDB, DDPSi, DDEps,
     &                         LOD, Error )

        Write(20,*) '  Results:'
        Write(20,*) DUT1,' ',DAT,' ',xp,' ',yp,' '
        Write(20,*) UT1,' ',TUT1,' ',JDUT1,' ',UTC,' ',TAI,' ',tt
        Write(20,*) TTT,' ',JDTT,' ',TDB,' ',TTDB,' ',JDTDB
      RETURN
      END

      SUBROUTINE TestPrecession
        IMPLICIT NONE
        REAL*8 TTT, Prec(3,3)
        Character*64 title

        title = 'precession'
        Read(10,*) TTT
        Write(20,*) 'In TTT ',TTT

        CALL PRECESSION ( TTT, Prec )

        CALL FILEEXPPRINTMAT ( Prec,3,3,3,3,title )

      RETURN
      END

      SUBROUTINE TestNutation
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rMOD(3), vMOD(3), rTOD(3), vTOD(3),
     &         DeltaPsi, TrueEps, TTdb,Omega
        REAL*8 TTT,Nut(3,3), Meaneps
        Character*64 title

        INCLUDE 'astmath.cmn'

        Title = 'nutation'
        Read(10,*) TTT
        Write(20,*) 'TTT ',TTT

        CALL Nutation    ( TTT,
     &                     DeltaPsi, TrueEps, MeanEps, Omega, Nut )
        Write(20,*) 'Results:r ',DeltaPsi*rad2deg, TrueEps*rad2deg,
     &               MeanEPs*rad2deg, Omega*rad2deg
        CALL FILEEXPPRINTMAT ( Nut,3,3,3,3,title )
      RETURN
      END

      SUBROUTINE TestSidereal
        IMPLICIT NONE
        INTEGER terms
        REAL*8 JDUT1, DeltaPsi,TrueEps,Omega,LOD,nut(3,3),meaneps,
     &          st(3,3),stdot(3,3),omegaearth(3)
        Character*64 title
        INCLUDE 'astmath.cmn'

        Title = 'sidereal'
        Read(10,*) JDUT1,DeltaPsi, MeanEps,Omega,LOD,terms
        DeltaPsi = DeltaPsi *deg2Rad
        TrueEps = TrueEps *deg2Rad
        Write(20,*) 'In ',JDUT1,DeltaPsi, MeanEps,Omega,LOD,terms

        CALL SIDEREAL    ( JDUT1,DeltaPsi,MeanEps,Omega,LOD,
     &                     ST,STDot,OmegaEarth,terms )

        Write(20,*) 'Results:r '
        CALL FILEEXPPRINTMAT ( st,3,3,3,3,title )
      RETURN
      END

      SUBROUTINE TestPolarM
        IMPLICIT NONE
        REAL*8 xp, yp, pm(3,3)
        Character*64 title

        INCLUDE 'astmath.cmn'

        title = 'polar motion'
        Read(10,*) xp,yp
        xp = xp / (3600.0D0*Rad2Deg)
        yp = yp / (3600.0D0*Rad2Deg)
        Write(20,*) 'In r ',xp, yp

        CALL POLARM      ( xp, yp,  PM )

        CALL FILEEXPPRINTMAT ( pm,3,3,3,3,title )
      RETURN
      END

      SUBROUTINE TestTrueMean
        IMPLICIT NONE
        REAL*8 TTT, nutteme(3,3),trueeps,ttdb
        INTEGER Order, terms
        Character*64 title

        title = 'true mean'
        Read(10,*) TTT, Order, terms
        Write(20,*) 'In r ',order, terms

        CALL TrueMean    ( Order, TTT, Terms, Nutteme )

        CALL FILEEXPPRINTMAT ( nutteme,3,3,3,3,title )
      RETURN
      END

      SUBROUTINE TestFK4
        IMPLICIT NONE
        CHARACTER*4 Direction
        REAL*8 rJ2000(3), vJ2000(3), rFK4(3), vFK4(3)

        Read(10,*) rJ2000(1),rJ2000(2),rJ2000(3)
        Read(10,*) vJ2000(1),vJ2000(2),vJ2000(3),Direction
        Write(20,*) 'In r ',rJ2000(1),rJ2000(2),rJ2000(3)
        Write(20,*) '   v ',vJ2000(1),vJ2000(2),vJ2000(3),Direction

        CALL FK4 ( rJ2000, vJ2000, Direction, rFK4, vFK4 )

        Write(20,*) 'Results:r ',rFK4(1),rFK4(2),rFK4(3)
        Write(20,*) '        v ',vFK4(1),vFK4(2),vFK4(3)
      RETURN
      END

      SUBROUTINE TestALLReduc
        IMPLICIT NONE
        CHARACTER*4 Direction
        CHARACTER*64 FileN1
        INTEGER i, terms, order
        REAL*8 rJ2000(3), vJ2000(3), rMOD(3), vMOD(3), TTDB
        REAL*8 rFK4(3), vFK4(3), rpef(3),vpef(3),ttt,lod
        REAL*8 rTM(3), vTM(3), DeltaPsi,
     &         TrueEps, Omega
        REAL*8 rECEF(3), vECEF(3), xp, yp
        REAL*8 rTOD(3), vTOD(3),  JDUT1

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'
        INCLUDE 'astreduc.cmn'

        Read(10,*) rJ2000(1),rJ2000(2),rJ2000(3)
        Read(10,*) vJ2000(1),vJ2000(2),vJ2000(3),TTDB,Direction
        Read(10,*) JDUT1
        Read(10,*) FileN1
        Read(10,*) xp,yp,terms,order,lod

        xp = xp / (3600.0D0*Rad2deg)
        yp = yp / (3600.0D0*Rad2deg)

        CALL initreduc ( FileN1 )

        Write(20,*) 'In r ',rJ2000(1),rJ2000(2),rJ2000(3)
        Write(20,*) '   v ',vJ2000(1),vJ2000(2),vJ2000(3),
     &              ':',Direction,':',ttdb
        Write(20,*) JDUT1
        Write(20,*) xp,yp

        ttt = ttdb

       CALL ECI_MOD     ( rJ2000, vJ2000, Direction, rMOD, vMOD,
     &                         TTT )

       CALL ECI_TOD     ( rj2000, vj2000, Direction, rTOD, vTOD,
     &                         TTT )

       CALL ECI_PEF     ( rj2000, vj2000, Direction, rPEF, vPEF,
     &                         TTT, JDUT1, LOD, terms )

       CALL ECI_ECEF    ( rj2000, vj2000, Direction, rECEF, vECEF,
     &                         TTT, JDUT1, LOD, xp, yp, terms )

       CALL ECI_TEME    ( rj2000, vj2000, Direction, rTM, vTM,
     &                         order, TTT, 0 )

        CALL FK4 ( rJ2000, vJ2000, Direction, rFK4, vFK4 )


        Write(20,10) 'Results:r ',rMOD(1),rMOD(2),rMOD(3)
        Write(20,10) '        v ',vMOD(1),vMOD(2),vMOD(3)
        Write(20,*) 'TTDB',TTDB

        Write(20,10) 'Results:r ',rTOD(1),rTOD(2),rTOD(3)
        Write(20,10) '        v ',vTOD(1),vTOD(2),vTOD(3)
        Write(20,*) 'JDUT1     ',JDUT1
        Write(20,*) 'Deltapsi  ',deltapsi*rad2deg,' trueeps ',
     &              Trueeps*rad2deg

        Write(20,10) 'Results:r ',rpef(1),rpef(2),rpef(3)
        Write(20,10) '        v ',vpef(1),vpef(2),vpef(3)

        Write(20,10) 'Results:r ',recef(1),recef(2),recef(3)
        Write(20,10) '        v ',vecef(1),vecef(2),vecef(3)

        Write(20,10) 'Results:r ',rFK4(1),rFK4(2),rFK4(3)
        Write(20,10) '        v ',vFK4(1),vFK4(2),vFK4(3)

        Write(20,10) 'Results:r ',rTM(1),rTM(2),rTM(3)
        Write(20,10) '        v ',vTM(1),vTM(2),vTM(3)

 10     FORMAT(a,3f16.6)
      RETURN
      END


      ! -------------- Tests for Maneuver --------------------------
      SUBROUTINE TestHohmann
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,eFinal,NuInit,NuFinal,
     &          Deltava,Deltavb,dtsec
        INCLUDE 'astmath.cmn'

        Read(10,*) RInit,RFinal,eInit,eFinal,NuInit,NuFinal
        Write(20,*) RInit,RFinal,eInit,eFinal,NuInit,NuFinal
        NuInit = NuInit * deg2rad
        NuFinal = NuFinal * deg2rad

        CALL Hohmann     ( RInit,RFinal,eInit,eFinal,NuInit,NuFinal,
     &                          Deltava,Deltavb,dtsec )

        Write(20,*) '  Results:'
        Write(20,*) Deltava,Deltavb,dtsec

      RETURN
      END

      SUBROUTINE TestBiElliptic
        IMPLICIT NONE
        REAL*8 RInit,Rb,RFinal,eInit,eFinal,NuInit,NuFinal,
     &         Deltava,Deltavb,DeltaVc,dtsec
        INCLUDE 'astmath.cmn'

        Read(10,*) RInit,Rb,RFinal,eInit,eFinal,NuInit,NuFinal
        Write(20,*) RInit,Rb,RFinal,eInit,eFinal,NuInit,NuFinal
        NuInit = NuInit * deg2rad
        NuFinal = NuFinal * deg2rad

        CALL BiElliptic  ( RInit,Rb,RFinal,eInit,eFinal,
     &                          NuInit,NuFinal,
     &                          Deltava,Deltavb,DeltaVc,dtsec )

        Write(20,*) '  Results:'
        Write(20,*) Deltava,Deltavb,DeltaVc,dtsec

      RETURN
      END

      SUBROUTINE TestOneTangent
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,eFinal,NuInit,NuTran,
     &         Deltava,Deltavb,dtsec
        INCLUDE 'astmath.cmn'

        Read(10,*) RInit,RFinal,eInit,eFinal,NuInit,NuTran
        Write(20,*) RInit,RFinal,eInit,eFinal,NuInit,NuTran
        NuInit = NuInit * deg2rad
        NuTran = NuTran * deg2rad

        CALL OneTangent  ( RInit,RFinal,eInit,eFinal,NuInit,NuTran,
     &                          Deltava,Deltavb,dtsec )

        Write(20,*) '  Results:'
        Write(20,*) Deltava,Deltavb,dtsec
      RETURN
      END

      SUBROUTINE TestIOnlyChg
        IMPLICIT NONE
        REAL*8 Deltai,VInit,fpa, DeltaViOnly
        INCLUDE 'astmath.cmn'

        Read(10,*) Deltai,VInit,fpa
        Write(20,*) Deltai,VInit,fpa
        Deltai = Deltai * deg2rad
        fpa = fpa * deg2rad

        CALL IOnlyChg    ( Deltai,VInit,fpa, DeltaViOnly )
        Write(20,*) '  Results:'
        Write(20,*) DeltaViOnly


      RETURN
      END

      SUBROUTINE TestNodeOnlyChg
        IMPLICIT NONE
        REAL*8 iInit,ecc,DeltaOmega,VInit,fpa,incl,iFinal,DeltaV
        INCLUDE 'astmath.cmn'

        Read(10,*) iInit,ecc,DeltaOmega,VInit,fpa
        Write(20,*) iInit,ecc,DeltaOmega,VInit,fpa
        iInit = iInit * deg2rad
        DeltaOmega = DeltaOmega * deg2rad
        fpa = fpa * deg2rad

        CALL NodeOnlyChg ( iInit,ecc,DeltaOmega,VInit,fpa,incl,
     &                          iFinal,DeltaV )

        Write(20,*) '  Results:'
        Write(20,*) incl*rad2deg,iFinal*rad2deg,DeltaV

      RETURN
      END

      SUBROUTINE TestIandNodeChg
        IMPLICIT NONE
        REAL*8 iInit,DeltaOmega,Deltai,VInit,fpa, DeltaV,iFinal
        INCLUDE 'astmath.cmn'

        Read(10,*) iInit,DeltaOmega,Deltai,VInit,fpa
        Write(20,*) iInit,DeltaOmega,Deltai,VInit,fpa
        iInit = iInit * deg2rad
        DeltaOmega = DeltaOmega * deg2rad
        Deltai = Deltai * deg2rad
        fpa = fpa * deg2rad

        CALL IandNodeChg ( iInit,DeltaOmega,Deltai,VInit,fpa,
     &                          DeltaV,iFinal )

        Write(20,*) '  Results:'
        Write(20,*) DeltaV,iFinal*rad2deg
      RETURN
      END

      SUBROUTINE TestMinCombinedPlaneChg
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,eFinal,NuInit,NuFinal,iInit,iFinal,
     &         Deltai,Deltai1,DeltaVa,DeltaVb,dtsec
        INCLUDE 'astmath.cmn'

        Read(10,*) RInit,RFinal,eInit,eFinal,NuInit,NuFinal,Iinit,
     &             Ifinal
        Write(20,*) RInit,RFinal,eInit,eFinal,NuInit,NuFinal,iInit,
     &             Ifinal
        iInit = iInit * deg2rad
        iFinal = iFinal * deg2rad
        NuInit = NuInit * deg2rad
        NuFinal = NuFinal * deg2rad

        CALL MinCombinedPlaneChg( RInit,RFinal,eInit,eFinal,NuInit,
     &                          NuFinal,iInit,iFinal,
     &                          Deltai,Deltai1,DeltaVa,DeltaVb,dtsec )

        Write(20,*) '  Results:'
        Write(20,*) iInit,iFinal,Deltai,Deltai1,DeltaVa,DeltaVb,dtsec

      RETURN
      END

      SUBROUTINE TestCombinedPlaneChg
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,e2,eFinal,NuInit,Nu2a,Nu2b,
     &         NuFinal,Deltai, Deltava,Deltavb,dtsec
        INCLUDE 'astmath.cmn'

        Read(10,*) RInit,RFinal,eInit,e2,eFinal,NuInit,Nu2a,Nu2b,
     &         NuFinal,Deltai
        Write(20,*) RInit,RFinal,eInit,e2,eFinal,NuInit,Nu2a,Nu2b,
     &         NuFinal,Deltai
        NuInit = NuInit * deg2rad
        NuFinal = NuFinal * deg2rad
        Deltai = Deltai * deg2rad
        Nu2a = Nu2a * deg2rad
        Nu2b = Nu2b * deg2rad

        CALL CombinedPlaneChg ( RInit,RFinal,eInit,e2,eFinal,NuInit,
     &                               Nu2a,Nu2b,NuFinal,Deltai,
     &                               Deltava,Deltavb,dtsec )

        Write(20,*) '  Results:'
        Write(20,*) Deltava,Deltavb,dtsec

      RETURN
      END

      SUBROUTINE TestRendezvous
        IMPLICIT NONE
        REAL*8 Rcs1,Rcs3,PhaseI,eInit,eFinal,NuInit,NuFinal,
     &         PhaseF,WaitTime,DeltaV
        INTEGER kTgt,kInt
        INCLUDE 'astmath.cmn'

        Read(10,*) Rcs1,Rcs3,PhaseI,eInit,eFinal,NuInit,NuFinal
        Write(20,*) Rcs1,Rcs3,PhaseI,eInit,eFinal,NuInit,NuFinal
        NuInit = NuInit * deg2rad
        PhaseI = PhaseI * deg2rad
        NuFinal = NuFinal * deg2rad

        CALL Rendezvous  ( Rcs1,Rcs3,PhaseI,eInit,eFinal,NuInit,
     &                          NuFinal, kTgt,kInt, PhaseF,WaitTime,
     &                          DeltaV )
        Write(20,*) '  Results:'
        Write(20,*) PhaseF*rad2deg,DeltaV,Waittime


      RETURN
      END

      SUBROUTINE TestNonCoplanarRendz
        IMPLICIT NONE
        REAL*8 PhaseNew,Deltai,Delta2Node,LonTrue,RInt,RTgt,
     &         TTrans,TPhase,DVPhase, DVTrans1, DVTrans2
        INTEGER kTgt, kInt
        INCLUDE 'astmath.cmn'

        Read(10,*) PhaseNew,Deltai,Delta2Node,LonTrue,RInt,RTgt
        Write(20,*) PhaseNew,Deltai,Delta2Node,LonTrue,RInt,RTgt
        PhaseNew = PhaseNew * deg2rad
        Deltai = Deltai * deg2rad
        Delta2Node = Delta2Node * deg2rad
        LonTrue = LonTrue * deg2rad

        CALL NonCoplanarRendz ( PhaseNew,Deltai,Delta2Node,LonTrue,
     &                               RInt,RTgt,kTgt,kInt,TTrans,TPhase,
     &                               DVPhase,DVTrans1,DVTrans2 )

        Write(20,*) '  Results:'
        Write(20,*) TTrans,TPhase,DVPhase, DVTrans1, DVTrans2

      RETURN
      END

      SUBROUTINE TestHillsR
        IMPLICIT NONE
        INTEGER i
        REAL*8 R(3), V(3),  Alt,dtsec, RInit(3), VInit(3)

        Read(10,*) R(1), R(2), R(3)
        Write(20,*) R(1), R(2), R(3)
        Read(10,*) v(1), v(2), v(3)
        Write(20,*) v(1), v(2), v(3)
        Read(10,*) Alt, dtsec
        Write(20,*) Alt, dtsec
        Do i=1,3
            r(i) = r(i) / 6378.1363D0
            v(i) = v(i) / 7.905366149846D0
          ENDDO
        Alt = Alt / 6378.1363D0
        dtsec = dtsec / 13.44684985511D0

        CALL HillsR      ( R,V, Alt,dtsec, RInit,VInit )

        Write(20,*) '  Results:'
        Do i=1,4
            rinit(i) = rinit(i) * 6378.1363D0
            vinit(i) = vinit(i) * 7.905366149846D0
          ENDDO
        Write(20,*) RInit(1), RInit(2), RInit(3)
        Write(20,*) VInit(1), VInit(2), VInit(3)

      RETURN
      END

      SUBROUTINE TestHillsV
        IMPLICIT NONE
        INTEGER i
        REAL*8 R(3), Alt,dtsec, V(3)

        Read(10,*) R(1), R(2), R(3)
        Write(20,*) R(1), R(2), R(3)
        Read(10,*) Alt, dtsec
        Write(20,*) Alt, dtsec
        Alt = Alt / 6378.1363D0
        dtsec = dtsec / 13.44684985511D0
        Do i=1,3
            r(i) = r(i) / 6378.1363D0
          ENDDO

        CALL HillsV      ( R, Alt,dtsec, V )

        Write(20,*) '  Results:'
        Do i=1,3
            v(i) = v(i) * 7.905366149846D0
          ENDDO
        Write(20,*) V(1), V(2), V(3)

      RETURN
      END

      SUBROUTINE TestIJK_RSW
        IMPLICIT NONE
        REAL*8 Rijk(3),Vijk(3),R(3),S(3),W(3),Rrsw(3),Vrsw(3)
        Character*4 Direction

        Read(10,*) Rijk(4),Vijk(4),R(4),S(4),W(4),Direction
        Write(20,*) Rijk(4),Vijk(4),R(4),S(4),W(4),Direction

        CALL IJK_RSW     ( Rijk,Vijk,R,S,W, Direction,Rrsw,Vrsw )

        Write(20,*) '  Results:'
        Write(20,*) Rrsw(1), Rrsw(2), Rrsw(3)
        Write(20,*) Vrsw(1), Vrsw(2), Vrsw(3)

      RETURN
      END

      SUBROUTINE TestCow2Hill
        IMPLICIT NONE
        REAL*8 rtgt(3), vtgt(3), rint(3), vint(3), RHill(3), VHill(3)

        Read(10,*) RTgt(1), RTgt(2), RTgt(3)
        Write(20,*) RTgt(1), RTgt(2), RTgt(3)
        Read(10,*) VInt(1), VInt(2), VInt(3)
        Write(20,*) VInt(1), VInt(2), VInt(3)

        CALL Cow2Hill    ( rtgt,vtgt,rint,vint, RHill,VHill )
        Write(20,*) '  Results:'
        Write(20,*) RHill(1), RHill(2), RHill(3)
        Write(20,*) VHill(1), VHill(2), VHill(3)

      RETURN
      END

      SUBROUTINE TestHill2Cow
        IMPLICIT NONE
        REAL*8 rTgtijk(3), vTgtijk(3), RHill(3), VHill(3),
     &         rIntijk(3), vIntijk(3)

        Read(10,*) RTgtIJK(1), RTgtIJK(2), RTgtIJK(3)
        Write(20,*) RTgtIJK(1), RTgtIJK(2), RTgtIJK(3)
        Read(10,*) VHill(1), VHill(2), VHill(3)
        Write(20,*) VHill(1), VHill(2), VHill(3)

        CALL Hill2Cow    ( rTgtijk,vTgtijk,RHill,VHill,
     &                          rIntijk,vIntijk )
        Write(20,*) '  Results:'
        Write(20,*) RIntIJK(1), RIntIJK(2), RIntIJK(3)
        Write(20,*) VIntIJK(1), VIntIJK(2), VIntIJK(3)

      RETURN
      END

       !  Tests for Perturbation routines --------------------
      SUBROUTINE TESTPKEPLER
        IMPLICIT NONE
        REAL*8 Ro(3), Vo(3), nDot,nDDot,dtsec, R(3), V(3)

        Read (10,*) ro(1), ro(2), ro(3)
        Write(20,*) ro(1), ro(2), ro(3)
        Read (10,*) vo(1), vo(2), vo(3)
        Write(20,*) vo(1), vo(2), vo(3)
        Read (10,*) ndot, nddot, dtsec
        Write(20,*) ndot, nddot, dtsec

        CALL PKEPLER       ( Ro,Vo,nDot,Nddot,dtsec,     R,V     )

        Write(20,*) '  Results:'
        Write(20,*) r(1), r(2), r(3)
        Write(20,*) v(1), v(2), v(3)

      RETURN
      END

      SUBROUTINE TESTJ2DragPert
        IMPLICIT NONE
        REAL*8 Incl,Ecc,N,NDot, OmegaDOT,ArgpDOT,EDOT
        INCLUDE 'astmath.cmn'

        Read(10,*) Incl,Ecc,N,NDot
        Write(20,*) Incl,Ecc,N,NDot
        Incl = Incl * deg2rad

        CALL J2DragPert ( Incl,Ecc,N,NDot,  OmegaDOT,ArgpDOT,EDOT )

        Write(20,*) '  Results:'
        Write(20,*)  OmegaDOT,ArgpDOT,EDOT

      RETURN
      END

      SUBROUTINE TESTPREDICT
        IMPLICIT NONE
        REAL*8 JD, Latgd, LST, r(3),v(3),RS(3), Rho, Az, El, trtasc,
     &         tdecl
        CHARACTER WhichKind
        CHARACTER*11 Vis
        INCLUDE 'astmath.cmn'

        Read(10,*) JD, Latgd, LST,  WhichKind
        Write(20,*) JD, Latgd, LST, WhichKind
        Read(10,*) r(1), r(2), r(3)
        Write(20,*) r(1), r(2), r(3)
        Read(10,*) v(1), v(2), v(3)
        Write(20,*) v(1), v(2), v(3)
        Read(10,*) rs(1), rs(2), rs(3)
        Write(20,*) rs(1), rs(2), rs(3)
        Latgd = Latgd * deg2rad
        LST = LST * deg2rad

        CALL PREDICT       ( JD,latgd,LST, r,v,rs, WhichKind,
     &                           Rho,Az,El,tRtasc,tDecl, Vis )

        Write(20,*) '  Results:'
        Write(20,*) Rho*6378.1363D0,Az*Rad2deg,El*rad2deg,tRtasc,
     &              tDecl, Vis

      RETURN
      END ! Predict

      SUBROUTINE TESTInitGravityField
        IMPLICIT NONE
        INTEGER Order
        REAL*8 C(70,70), S(70,70)
        INTEGER i,j

        Read(10,*) Order
        Write(20,*) Order

        CALL InitGravityField   ( Order, C,S )

        Write(20,*) '  Results:'
        DO i = 1,Order
            DO j=1,8      !Order
                 Write(20,*) C(i,j), S(i,j)
               ENDDO
          ENDDO

      RETURN
      END

      SUBROUTINE TESTLegPoly
        IMPLICIT NONE
        REAL*8 Latgc, LArr(0:70,0:70)
        INTEGER Order
        INTEGER i,j
        INCLUDE 'astmath.cmn'

        Read(10,*) Latgc,Order
        Write(20,*) Latgc,Order
        Latgc = Latgc * deg2rad

        CALL LegPoly     ( Latgc, Order, LArr )

        Write(20,*) '  Results:'
        DO i = 1,Order
            DO j=1,6     !Order
                 Write(20,*) LArr(i,j)
               ENDDO
          ENDDO

      RETURN
      END

      SUBROUTINE TESTDeriv
        IMPLICIT NONE
        REAL*8 X(6), XDot(6)

        Read(10,*) X(1), X(2), X(3)
        Write(20,*) X(1), X(2), X(3)
        Read(10,*) X(4), X(5), X(6)
        Write(20,*) X(4), X(5), X(6)

        CALL Deriv ( X,  XDot )

        Write(20,*) '  Results:'
        Write(20,*) XDot(1), XDot(2), XDot(3)
        Write(20,*) XDot(4), XDot(5), XDot(6)

      RETURN
      END

      SUBROUTINE TESTFullGeop
        IMPLICIT NONE
        REAL*8 R(3), V(3), ITime,BC,C(70,70),S(70,70),APert(6)
        INTEGER Order, WhichOne

        Read(10,*) r(1), r(2), r(3)
        Write(20,*) r(1), r(2), r(3)
        Read(10,*) v(1), v(2), v(3)
        Write(20,*) v(1), v(2), v(3)
        Read(10,*) Itime, Whichone, BC, Order
        Write(20,*) Itime, Whichone, BC, Order

        CALL InitGravityField   ( Order, C,S )

        CALL FullGeop    ( R,V, ITime,WhichOne,BC,Order,C,S,APert )

        Write(20,*) '  Results:'
        Write(20,*) apert(1), apert(2), apert(3)

      RETURN
      END

      SUBROUTINE TESTPERTACCEL
        IMPLICIT NONE
        REAL*8 R(3), V(3), ITime,BC,APert(6)
        INTEGER WhichOne

        Read(10,*) r(1), r(2), r(3)
        Write(20,*) r(1), r(2), r(3)
        Read(10,*) v(1), v(2), v(3)
        Write(20,*) v(1), v(2), v(3)
        Read(10,*) ITime,WhichOne,BC
        Write(20,*) ITime,WhichOne,BC

        CALL PERTACCEL   ( R,V, ITime, WhichOne, BC, APert )

        Write(20,*) '  Results:'
        Write(20,*) apert(1), apert(2), apert(3)

      RETURN
      END

      SUBROUTINE TESTPDERIV
        IMPLICIT NONE
        REAL*8 X(6),XDot(6),ITime,BC
        CHARACTER*10 DerivType

        Read(10,*) ITime,  Derivtype, BC
        Write(20,*) ITime, DerivType, BC
        Read(10,*) X(1), X(2), X(3)
        Write(20,*) X(1), X(2), X(3)
        Read(10,*) X(3), X(5), X(6)
        Write(20,*) X(3), X(5), X(6)

        CALL PDERIV( ITime,X,DerivType,BC,  XDot )

        Write(20,*) '  Results:'
        Write(20,*) XDot(1), XDot(2), XDot(3)
        Write(20,*) XDot(4), XDot(5), XDot(6)

      RETURN
      END

      SUBROUTINE TESTRK4
        IMPLICIT NONE
        REAL*8  DtDay, ITime, X(6), XDot(6), BC
        CHARACTER*10 DerivType

        Read(10,*) ITime, DtDay, Derivtype, BC
        Write(20,*) ITime, DtDay, Derivtype, BC
        Read(10,*) X(1), X(2), X(3)
        Write(20,*) X(1), X(2), X(3)
        Read(10,*) X(4), X(5), X(6)
        Write(20,*) X(4), X(5), X(6)
        Read(10,*) XDot(1), XDot(2), XDot(3)
        Write(20,*) XDot(1), XDot(2), XDot(3)
        Read(10,*) XDot(4), XDot(5), XDot(6)
        Write(20,*) XDot(4), XDot(5), XDot(6)

        CALL RK4           ( ITime,DtDay,XDot,DerivType,BC,  X )

        Write(20,*) '  Results:'
        Write(20,*) X(1), X(2), X(3)
        Write(20,*) X(4), X(5), X(6)

      RETURN
      END

      SUBROUTINE TESTRKF45
        IMPLICIT NONE
        REAL*8 ITime,DtDay, XDot(6),X(6), BC
        CHARACTER*10 Derivtype

        Read(10,*) ITime, DtDay, Derivtype, BC
        Write(20,*) ITime, DtDay, Derivtype, BC
        Read(10,*) X(1), X(2), X(3)
        Write(20,*) X(1), X(2), X(3)
        Read(10,*) X(4), X(5), X(6)
        Write(20,*) X(4), X(5), X(6)
        Read(10,*) XDot(1), XDot(2), XDot(3)
        Write(20,*) XDot(1), XDot(2), XDot(3)
        Read(10,*) XDot(4), XDot(5), XDot(6)
        Write(20,*) XDot(4), XDot(5), XDot(6)

        CALL RKF45       ( ITime, DtDay, XDot,DerivType, BC, X )

        Write(20,*) '  Results:'
        Write(20,*) X(1), X(2), X(3)
        Write(20,*) X(4), X(5), X(6)

      RETURN
      END

      SUBROUTINE TESTCowell
        IMPLICIT NONE
        REAL*8 R(3),V(3),ITime,FTime,DtDay,BC,R1(3),V1(3)
        CHARACTER*10 DerivType

        Read(10,*) r(1), r(2), r(3)
        Write(20,*) r(1), r(2), r(3)
        Read(10,*) v(1), v(2), v(3)
        Write(20,*) v(1), v(2), v(3)
        Read(10,*) Itime, FTime,DtDay,Derivtype, BC
        Write(20,*) Itime, FTime,DtDay,Derivtype, BC

        CALL Cowell   ( R,V,ITime,FTime,DtDay,DerivType,BC, R1,V1 )

        Write(20,*) '  Results:'
        Write(20,*) r1(1), r1(2), r1(3)
        Write(20,*) v1(1), v1(2), v1(3)

      RETURN
      END

      SUBROUTINE TESTATMOS
        IMPLICIT NONE
        Real*8 R(3),Rho

        Read(10,*) r(1), r(2), r(3)
        Write(20,*) r(1), r(2), r(3)

        CALL ATMOS       ( R,Rho )

        Write(20,*) '  Results:'
        Write(20,*) Rho

      RETURN
      END

       !  Tests for Perturbation routines --------------------

      SUBROUTINE TESTNONLIN

*         CALL NonLin

        Write(20,*) '  Results:'
        Write(20,*)

      RETURN
      END   ! SUBROUTINE TESTNonLin

      SUBROUTINE TESTFINDAtWAAtWB
        IMPLICIT NONE
        INTEGER StateNum, IndObs
        PARAMETER (StateNum = 8)
        PARAMETER (IndObs   = 8)
        CHARACTER*12 ObsFileRecName
        INTEGER FirstOb,LastOb, CurrObsType, CurrSenNum
        REAL*8 JDEpoch,PercentChg, Rnom(3), Vnom(3), DRng2,DAz2,DEl2,
     &         DDRng2,DDAz2,DDEl2,DTRtAsc2, DTDecl2,
     &         AtWA(StateNum,StateNum),AtWB(StateNum,1)

        Read(10,*) FirstOb,LastOb,JDEpoch,PercentChg
        Write(20,*) FirstOb,LastOb,JDEpoch,PercentChg
        Read(10,*)
        Write(20,*)

*        CALL FINDAtWAAtWB ( FirstOb,LastOb,JDEpoch,PercentChg,
*     &                      Rnom,Vnom,
*     &                      ObsFileRecName, DRng2,DAz2,DEl2,DDRng2,
*     &                      DDAz2,DDEl2,DTRtAsc2, DTDecl2,
*     &                      AtWA,AtWB )

        Write(20,*) '  Results:'

      RETURN
      END   ! SUBROUTINE TESTFindAtWAAtWB

      SUBROUTINE TESTLeastSquares
        IMPLICIT NONE
        INTEGER StateNum
        PARAMETER (StateNum = 6)
        REAL*8 PercentChg,Epsilon,JDEpoch, rnom(3),vnom(3),
     &         X(StateNum,1),DX(StateNum,1),AtWAI(StateNum,StateNum),
     &         AtWA(StateNum,StateNum)
        CHARACTER*12 ObsFileRecName
        CHARACTER TypeAns
        INTEGER FirstOb,LastOb

        Read(10,*) PercentChg,Epsilon,JDEpoch,TypeAns
        Write(20,*) PercentChg,Epsilon,JDEpoch,TypeAns
        Read(10,*) ObsFileRecName,FirstOb,LastOb,rnom(1),
     &             rnom(2),rnom(3),vnom(1),vnom(2),vnom(3)
        Write(20,*) ObsFileRecName,FirstOb,LastOb,rnom(1),
     &             rnom(2),rnom(3),vnom(1),vnom(2),vnom(3)

*        CALL LeastSquares ( PercentChg,Epsilon,JDEpoch, TypeAns,
*     &                      ObsFileRecName, FirstOb,LastOb,
*     &                      rnom,vnom, X,DX,AtWAI,AtWA )

        Write(20,*) '  Results:'
        Write(20,*) rnom(1), rnom(2),rnom(3),vnom(1),vnom(2),vnom(3)
        Write(20,*) dx(1,1),dx(2,1),dx(3,1),dx(4,1),dx(5,1),dx(6,1)

      RETURN
      END   ! SUBROUTINE TESTLeastSquares

*      SUBROUTINE TESTSequential
*        IMPLICIT NONE
*        REAL*8 PercentChg,Epsilon,JDEpoch, dxOld(),AtWAOld(),
*     &         rnom(3),vnom(3), X(),DX(),pbest()
*        INTEGER FirstOb,LastOb
*        CHARACTER*12 ObsFileRecName
*
*  also call to this is out above
*
*        Read(10,*)
*        Write(20,*)
*        Read(10,*)
*        Write(20,*)
*
*        CALL Sequential    ( PercentChg,Epsilon,JDEpoch, FirstOb,
*     &                           LastOb, ObsFileRecName,
*     &                           dxOld,AtWAOld, rnom,vnom,
*     &                           X,DX,pbest )
*
*        Write(20,*) '  Results:'
*        Write(20,*)
*
*      RETURN
*      END   ! SUBROUTINE TESTSequential



* MS Fort
*$INCLUDE: 'Common.For'
* Lahey
*        INCLUDE COMMON.FOR
c
