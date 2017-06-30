*   -------------------------------------------------------------------
*
*                              ASTUTIL.FOR
*
*   This file contains some utility routines for character operations.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                   2007
*                             by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*            30 may 07  david vallado
*                           3rd edition baseline
*    changes :
*            28 Jan 04  David Vallado
*                         Fix headers
*            28 Feb 03  David Vallado
*                         New baseline
*            14 May 01  David Vallado
*                         2nd edition baseline
*            10 Nov 98  David Vallado
*                         Original baseline
*
*     *****************************************************************
*
*
*      -------------------- Misc functions for parsing ----------------
*
*      CHARACTER*255 FUNCTION COPY     ( InStr, LocStart,Leng )
*
*      CHARACTER*255 FUNCTION RmSpcs   ( InStr, LengIn, LenOut )
*
*      INTEGER FUNCTION GETPART        ( InStr, LocStart,Leng )
*
*      INTEGER*4 FUNCTION GETPARTL     ( InStr, LocStart,Leng )
*
*      REAL*8 FUNCTION GetPartR        ( InStr, LocStart,Leng )
*
*      CHARACTER*250 FUNCTION UpCaseSt ( S )
*
*
      CHARACTER*255 FUNCTION COPY      ( InStr, LocStart,Leng )
        IMPLICIT NONE
        INTEGER LocStart, Leng
        CHARACTER*255 InStr

* -----------------------------  Locals  ------------------------------
        CHARACTER*255 TempStr

        ! --------------------  Implementation   ----------------------
        TempStr = InStr(LocStart:LocStart+Leng)
        COPY= TempStr
      RETURN
      END


      CHARACTER*255 FUNCTION RMSPCS     ( InStr,LengIn, LengOut )
        IMPLICIT NONE
        CHARACTER*255 InStr
        INTEGER LengIn, LengOut

* -----------------------------  Locals  ------------------------------
        CHARACTER*255 TempStr
        INTEGER iStart,iEnd, i

        ! --------------------  Implementation   ----------------------

        i = 1
        DO WHILE ((i.lt.LengIn).and.(Instr(i:i).eq.' '))
            i = i + 1
          ENDDO
        iStart = i

c       Go from the end back to first non blk
        i = LengIn
        DO WHILE ((i.gt.iStart).and.(Instr(i:i).eq.' '))
            i = i - 1
          ENDDO
        iEnd = i

        LengOut = iEnd - iStart + 1
        TempStr = InStr(iStart:iEnd)

        Rmspcs = TempStr
      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           FUNCTION GETPART
*
*  These functions parse a section of a string and return various types of
*    values - Real, INTEGER, Byte, etc., depending on the last letter in the
*    name.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    InStr       - Input String
*    LocStart    - Location to start parsing
*    Leng        - Length of dersired variable
*
*  Outputs       :
*    FUNCTION    - variable of correct type
*
*  Locals        :
*    TempStr     - Temporary string
*    Temp        - Temporary variable
*    Code        - Error variable
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------

      INTEGER FUNCTION GETPART         ( InStr, LocStart,Leng )
        IMPLICIT NONE
        INTEGER LocStart, Leng
        EXTERNAL RmSpcs
        CHARACTER*255 InStr,RmSpcs
* -----------------------------  Locals  ------------------------------
        INTEGER Temp,i,j,k,LengOut, FirstPlc
        CHARACTER*255 TempStr
        LOGICAL Neg

        ! --------------------  Implementation   ----------------------
        Neg = .FALSE.
        TempStr= InStr(LocStart:LocStart+Leng-1)
        TempStr = RmSpcs( TempStr,Leng,LengOut )

        IF (Index(InStr,'-').ne.0) THEN
            FirstPlc = 1
            DO WHILE ((FirstPlc.lt.Leng).and.
     &                (InStr(FirstPlc:FirstPlc).eq.' '))
                FirstPlc = FirstPlc + 1
              ENDDO
            IF (FirstPlc.eq.Index(Instr,'-')) THEN
                Neg = .TRUE.
                TempStr= InStr(INDEX(InStr,'-')+1:LocStart+Leng-1)
                k = Leng - (Index(Instr,'-') - LocStart + 1)
              ELSE
                k = LengOut
              ENDIF
          ELSE
            k = LengOut
          ENDIF

        Temp = 0
        DO i = 1,k
            j = ICHAR(TempStr(i:i)) - 48
            Temp = Temp + j*(10**(k-i))
          ENDDO
        IF (Neg.eqv..TRUE.) Temp = -Temp

        GETPART= Temp
      RETURN
      END
*
      INTEGER*4 FUNCTION GETPARTL      ( InStr, LocStart,Leng )
        IMPLICIT NONE
        INTEGER LocStart, Leng
        EXTERNAL RmSpcs
        CHARACTER*255 InStr, RmSpcs

* -----------------------------  Locals  ------------------------------
        INTEGER*4 Temp
        INTEGER i,j,k,LengOut,FirstPlc
        CHARACTER*255 TempStr
        LOGICAL Neg

        ! --------------------  Implementation   ----------------------
        Neg = .FALSE.
        TempStr= InStr(LocStart:LocStart+Leng-1)
        TempStr = RmSpcs( TempStr,Leng,LengOut )

        IF (Index(InStr,'-').ne.0) THEN
            FirstPlc = 1
            DO WHILE ((FirstPlc.lt.Leng).and.
     &                (InStr(FirstPlc:FirstPlc).eq.' '))
                FirstPlc = FirstPlc + 1
              ENDDO
            IF (FirstPlc.eq.Index(Instr,'-')) THEN
                Neg = .TRUE.
                TempStr= InStr(INDEX(InStr,'-')+1:LocStart+Leng-1)
                k = Leng - (Index(Instr,'-') - LocStart + 1)
              ELSE
                k = LengOut
              ENDIF  
          ELSE
            k = LengOut
          ENDIF

        Temp = 0
        DO i = 1,k
            j = ICHAR(TempStr(i:i)) - 48
            Temp = Temp + j*(10**(k-i))
          ENDDO
        IF (Neg.eqv..TRUE.) Temp = -Temp

        GETPARTL= Temp
      RETURN
      END
*
      REAL*8 FUNCTION GetPartR         ( InStr, LocStart,Leng )
        IMPLICIT NONE
        INTEGER LocStart, Leng
        CHARACTER*255 InStr
* -----------------------------  Locals  ------------------------------
        INTEGER i,j
        REAL*8 TempR1,TempR2
        LOGICAL Fraction
        CHARACTER*255 TempStr

        ! --------------------  Implementation   ----------------------
        TempStr= InStr(LocStart:LocStart+Leng)
c do negatives!!

        TempR1 = 0.0D0
        TempR2 = 0.0D0
        Fraction = .FALSE.
        DO i = 1,Leng
            j = ICHAR(TempStr(i:i)) - 48
            IF (Fraction .eqv..FALSE.) THEN
                TempR1 = TempR1 + j*(10.0D0**(Leng-i))
              ELSE
                TempR2 = TempR2 + j*(10.0D0**(-Leng+i))
              ENDIF
          ENDDO

        GetPartR= TempR1 + TempR2
      RETURN
      END

      CHARACTER*250 FUNCTION UpCaseSt  ( S,lens )
        IMPLICIT NONE
        CHARACTER*250 S
        INTEGER Lens

* ----------------------------  Locals  -------------------------------
        INTEGER i,ChrVal

        ! --------------------  Implementation   ----------------------

        DO i = 1, LenS
             ChrVal = ICHAR( S(i:i) )
             IF ( (ChrVal.ge.97).and.(ChrVal.le.122) ) THEN
                 S(i:i) = CHAR( ChrVal - 32 )
               ENDIF
          ENDDO
        UpCaseSt = S
      RETURN
      END

*
