*   -------------------------------------------------------------------
*
*                              ASTDC.FOR
*
*   This file contains fundamental Astrodynamic procedures and functions
*   relating to Differential Correction functions.  See chap 10 for a
*   detailed discussion and derivation of the fundamental routines.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                   2007
*                             by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              15 mar 07  david vallado
*                           3rd edition baseline
*    changes :
*            28 Jan 04  David Vallado
*                         Fix headers
*            28 Feb 03  David Vallado
*                         New baseline
*            14 May 01  David Vallado
*                         2nd edition baseline
*            23 Nov 87  David Vallado
*                         Original baseline
*
*     *****************************************************************
*
*  Uses object files:
*      astMath,
*      asttime,
*      ast2body,
*      astiod,
*      senset, ! because str12 is in senset
*
*     TYPE
*         Individual types are as follows:
*                        0 DRng
*                        1 Az, El
*                        2 Rng, Az, El
*                        3 Rng, Az, El, DRng
*                        4 Rng, Az, El, DRng, DAz, DEl
*                        5 TRtAsc, TDecl
*                        6 Rng
*           ObsRec = RECORD
*                      SenNum   : INTEGER
*                      JD       : REAL*8
*                      RSVec(3) : REAL*8
*                      ObsType  : BYTE
*                      Rng,Az,El,DRng,DAz,DEl,TRtAsc,TDecl : REAL*8
*                    ENDIF
*           ObsRecType = File of ObsRec

*      SUBROUTINE NONLIN
*
*      SUBROUTINE FINDAtWAAtWB ( FirstOb,LastOb,JDEpoch,PercentChg,
*     &                          Rnom,Vnom,NDot,NDDot,
*     &                          DRng2,DAz2,DEl2,DDRng2,
*     &                          DDAz2,DDEl2,DTRtAsc2, DTDecl2,
*     &                          AtWA,AtWB )
*
*      SUBROUTINE LeastSquares ( PercentChg,Epsilon,JDEpoch, TypeAns,
*     &                          FirstOb,LastOb,
*     &                          Rnom,Vnom,NDot,NDDot,
*     &                          X,DX,AtWAI,AtWA )
*
*      SUBROUTINE Sequential    ( PercentChg,Epsilon,JDEpoch, FirstOb,
*     &                           LastOb,
*     &                           dxOld,AtWAOld, rnom,vnom,
*     &                           X,DX,pbest )
*
*      SUBROUTINE FindPhi       ( Rnom,Vnom,Dtsec,PercentChg,Phi)
*
*      SUBROUTINE FINDH         ( XBest, CurrObsRec, ObsRecFile,
*     &                           PercentChg,Dtsec,CurrSenPtr,LST,H )
*
*      SUBROUTINE LKF           ( PercentChg,Epsilon, JDEpoch,Q,R,
*     &                           FirstOb,LastOb, ObsRecFile, XOld,
*     &                           rnom,vnom, Xbest,dxbest,pbest )
*
*      SUBROUTINE EKF           ( PercentChg,Epsilon, JDEpoch,Q,R,
*     &                           FirstOb,LastOb,ObsRecFile,rnom,
*     &                           vnom,Xbest,dxbest,pbest )
*
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE NONLIN
*
* this subroutine finds the A and b matrices for the differential correction
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    NumObs      - Number of observations
*    Epoch       - Epoch to which the DC places all answers
*    PercentChg  - Amount to change the vectors by in finite differencing
*
*
*
*  OutPuts       :                                    *    A
*
*
*  Locals        :
*    z23
*
*
*
*
*  Coupling      :
*    KEPLER      - Propagation
*
*  References    :
*    Vallado       1997, 675-682
*
* -----------------------------------------------------------------------------

      SUBROUTINE NONLIN
* ----------------------------  Locals  -------------------------------
        INTEGER NumObs, StateNum
        PARAMETER (NumObs   = 4)
        PARAMETER (StateNum = 2)

        REAL*8 DX(StateNum,1),a(NumObs,StateNum),b(NumObs,StateNum),
     &         atw(StateNum,NumObs),X(StateNum,1),at(StateNum,NumObs),
     &         AtWA(StateNum,StateNum),AtWAi(StateNum,StateNum),
     &         AtWB(NumObs,1),y(NumObs)
        REAL*8 xx(20)
        REAL*8 sigmaold,sigmanew,w1,w2, AlpNom,BetNom,epsilon
        INTEGER i,j,Loops
        CHARACTER*64 Title

        ! --------------------  Implementation   ----------------------
        ! -------- Get nominal estimate for problem -------------------
        betnom= DLOG( 50.0D0/19.0D0 ) / DLOG (4.0D0/3.0D0)
        alpnom= 19.0D0/(3.0D0**betnom)   ! use 3rd point now
        write(*,'(A,F11.7,F11.7)') 'alpnom 7 betnom ',alpnom,betnom
* 
        ! -- Calculate corrections based on Least Squares Analysis ----
        Epsilon= 0.01D0 
        Write(*,*)  ' Enter RMS tolerance to stop '
        Read(*,*) Epsilon

        SigmaOld=   1.0D0
        SigmaNew= 100.0D0 
        Loops   = 0 

        ! ----------------------- Find W matrix -----------------------
*        W1= 1.0D0 / (SigmaX*SigmaX)
*        W2= 1.0D0 / (SigmaY*SigmaY)
        W1= 1.0D0 
        W2= 1.0D0 

        Write(20,'(A,F6.4)') ' Epsilon DO convergence ',epsilon
*
        DO WHILE ((DABS( (SigmaNew-SigmaOld)/SigmaOld ) .gt. Epsilon)
     &           .and. (Loops .lt. 10))
            SigmaOld= SigmaNew
            IF ( Loops .eq. 0 ) THEN
*                INITMATRIX( StateNum,1,X )
              ENDIF
            X(1,1) = alpNom
            X(2,1) = betNom
            Write(20,'(A,F12.8,f12.8)') 'StateVect  ',  X(1,1), X(2,1)

            DO i= 1, NumObs

                ! --------- Find B matrix -----------------------------
                B(i ,1)= Y(i) - alpnom*xx(i)**betnom

                ! --------- Find A matrix -----------------------------
                A(i,1)= xx(i) ** Betnom
                A(i,2)= Alpnom* xx(i) ** Betnom *DLOG(xx(i))
              ENDDO

            Title = ' B Matrix '
c            CALL FILEwriteMAT( B,Numobs, 1,StateNum,1, Title )
*            CALL MATTRANS  ( A         , NumObs,StateNum,8,8  , At  )

            IF ( Loops .eq. 0 ) THEN
*               CALL writeMAT( at,4,'At ' )
              ENDIF

            DO i= 1, StateNum
                j= 1
                DO WHILE (j .le. NumObs)
                    AtW(i,j) =  At(i,j )*W1
                    AtW(i+1,j) =  At(i+1,j )*W2
                    j= j + 1 
                  ENDDO
              ENDDO
            IF ( Loops .eq. 0 ) THEN
*                CALL writeMAT( atw,4,'Atw ' )
              ENDIF

*            CALL MATMULT   ( AtW,A , StateNum,NumObs,StateNum,8,8,AtWA)
            IF ( Loops .eq. 0 ) THEN
*                CALL writeMAT( AtWA,4,'Ata ' )
              ENDIF

*            CALL FILEwriteMAT( A,Statenum,1,' A Matrix ' )

            CALL MATINVERSE( AtWA      , StateNum ,8, AtWAI )
*            CALL FILEEXPwriteMAT( AtWAi,Statenum,StateNum,
*     &                            ' AtWAI - COVARIANCE Matrix ' )

*            CALL MATMULT   ( AtW,B     , StateNum,NumObs,1     , AtWB )
            IF ( Loops .eq. 0 ) THEN
*                CALL writeMAT( AtWB,4,'Atb ' )
              ENDIF
*            CALL MATMULT   ( AtWAI,AtWB, StateNum,StateNum,1     , DX )

            alpNom= alpNom +   DX(1,1 )
            betNom= betNom +   DX(2,1 )

            Write(*,'(A,2(F12.8))') 'Corrections    ', DX(1,1), DX(2,1)
            Write(20,'(A,2(F12.8))') 'Corrections    ', DX(1,1), DX(2,1)
            Write(20,*)

            SigmaNew= DSQRT( (b(1,1)**2 + b(2,1)**2 +
     &                        b(3,1)**2 + b(4,1)**2 ) / (NumObs-1) )
            Write(*,'(A,I2,2(A,F12.6),F14.8)')  'Iteration # ',loops,
     &            ' SigmaOld ',sigmaold,' SigmaNew ',SigmaNew,
     &            DABS( (SigmaNew-SigmaOld)/SigmaOld)
            Write(20,'(A,I2,2(A,F12.6),F14.8)')  'Iteration # ',loops,
     &            ' SigmaOld ',sigmaold,' SigmaNew ',SigmaNew,
     &            DABS( (SigmaNew-SigmaOld)/SigmaOld)
            Write(20,*) ' ------------------------------------------- '

            Loops= Loops + 1 

          ENDDO   ! DO WHILE processing loop

      RETURN
      END   ! Subroutine NonLin
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE FINDAtWBAtWA
*
* this subroutine finds the A .and. b matrices DO the differential correction
*   problem.  Remember that it isn't critical DO the propagations to use
*   the highest fidelity techniques because we're only trying to find the
*   "slope". k is an index that allows us to do multiple rows at once. It's
*   used for both the b .and. A matrix calculations.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    NumObs      - Number of observations
*    JDEpoch     - Julian date DO solution       Days from 4713 BC
*    PercentChg  - Amount to modify the vectors
*                  by in finite differencing
*    RNom        - Nom position vector at epoch   ER
*    VNom        - Nom velocity vector at epoch   ER/TU
*    ObsRecFile    - Array of RECORDs containing:
*                  SenNum, JD, RSVec, ObsType,
*                  Rng, Az, El, DRng, DAz, DEl,
*                  TRtasc, TDecl
*
*  OutPuts       :
*    DRng2       - Range residual squared
*    DAz2        - Range residual squared
*    DEl2        - Range residual squared
*    DTRtAsc2    -
*    DTDecl2     -
*    A           - A matrix
*    B           - B matrix
*    AtWA        - AtWA matrix
*    AtWB        - AtWB matrix
*
*  Locals        :
*    IndObs      -
*    At          -
*    AtW         -
*    W1          -
*    W2          -
*    W3          -
*    LST         -
*    GST         -
*    Dtsec        -
*    DeltaAmt    -
*    RngPert     -
*    AzPert      - Modified azimuth               -2Pi to 2Pi
*    ElPert      - Modified azimuth               -Pi/2 to Pi/2
*    DRng        -
*    DAz         -
*    DEl         -
*    Error       -
*    i           -
*    j           -
*    k           -
*
*  Coupling      :
*    KEPLER      - Propagation
*    LSTIME      - Local SIDEREAL Time
*    Findsenptr  - Find sensor data
*    RV_RAZEL    - Find r .and. v given range, az, el, .and. rates
*    RV_TRadec   - Find r .and. v given topocentric rtasc .and. decl
*
*
*  References    :
*    Vallado       1997, 682-686
*
* -----------------------------------------------------------------------------
*
      SUBROUTINE FINDAtWAAtWB ( FirstOb,LastOb,JDEpoch,PercentChg,
     &                          Rnom,Vnom,NDot,NDDot,
     &                          DRng2,DAz2,DEl2,DDRng2,
     &                          DDAz2,DDEl2,DTRtAsc2, DTDecl2,
     &                          AtWA,AtWB )
        IMPLICIT NONE
        INTEGER StateNum, IndObs
        PARAMETER (StateNum = 8)
        PARAMETER (IndObs   = 8)
        INTEGER FirstOb,LastOb
        REAL*8 JDEpoch,PercentChg, Rnom(3), Vnom(3), DRng2,DAz2,DEl2,
     &         DDRng2,DDAz2,DDEl2,DTRtAsc2, DTDecl2,NDot,NDDot,
     &         AtWA(StateNum,StateNum),AtWB(StateNum,1)
* ----------------------------  Locals  -------------------------------
        INTEGER rowc,colc, r,c, i,j, NObs
        REAL*8 r3(3),v3(3),rPert(3),Vpert(3),A(StateNum,StateNum),
     &        B(IndObs,1),AtWAacc(StateNum,StateNum),r1(3),v1(3),
     &        AtWBacc(StateNum,1),At(StateNum,IndObs),
     &        Atw(StateNum,IndObs), W1,W2,W3,W4,W5,W6, LST,GST,Dtsec,
     &        DeltaAmt,RngPert,AzPert,ElPert,DRngPert,DAzPert,DElPert,
     &        TRtAscNom,TDeclNom,DTRtAscNom,DTDeclNom,TRtAscPert,
     &        TDeclPert,DTRtAscPert,DTDeclPert,Weight, RngNom,AzNom,
     &        ElNom,DRngNom,DAzNom,DElNom
        CHARACTER*12 error
        CHARACTER*64 Title

         REAL*8 latgd,lon,alt,sec
         INTEGER Year,month,day,hr,minute

* ----------------------------  Commons  ------------------------------
        COMMON /SensorRec/
     &         SenNum, SenName, SenLat, SenLon, SenAlt,
     &         RngMin, RngMax, AzMin, AzMax, ElMin, ElMax,
     &         BiasRng, BiasAz, BiasEl, BiasDRng, BiasDAz, BiasDEl,
     &         BiasTRtAsc, BiasTDecl,
     &         NoiseRng, NoiseAz, NoiseEl, NoiseDRng, NoiseDAz,
     &         NoiseDEl, NoiseTRtAsc, NoiseTDecl
        INTEGER SenNum
        CHARACTER*36 SenName
        REAL*8 SenLat, SenLon, SenAlt,
     &         RngMin, RngMax, AzMin, AzMax, ElMin, ElMax,
     &         BiasRng, BiasAz, BiasEl, BiasDRng, BiasDAz, BiasDEl,
     &         BiasTRtAsc, BiasTDecl,
     &         NoiseRng, NoiseAz, NoiseEl, NoiseDRng, NoiseDAz,
     &         NoiseDEl, NoiseTRtAsc, NoiseTDecl
        COMMON /ObsRec/
     &         JD, RSVec, ObsType, SensNum,
     &         Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl
        INTEGER ObsType, SensNum
        REAL*8 JD, RSVec(3), Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl

         REAL*8 MFME, DUT1, DAT, xp, yp, LOD, DDPsi, DDEps,Conv1
         REAL*8 MJD, TTT,
     &          jdut1, deltapsi,trueeps,meaneps,omega,
     &          pm(3,3),nutteme(3,3),rteme(3),vteme(3)
         REAL*8  UT1, TUT1, UTC, TAI, TDT, JDTDT, TDB, TTDB, JDTDB

         INTEGER RecNum,timezone,testnum,
     &           terms, order

         INCLUDE 'astreduc.cmn'

        ! --------------------  Implementation   ----------------------

        ! -------  Figure out the obstype and assign the NObs ---------
        Read(15,REC=1) SensNum, JD, (RSVec(i),i=1,3),
     &                 ObsType, Rng, Az, El, DRng,
     &                 Daz, Del, TRtasc, TDecl
        IF ( (ObsType.eq.0).or.(ObsType.eq.6) ) THEN
            NObs= 1
          ELSEIF ( (ObsType.eq.1).or.(ObsType.eq.5) ) THEN
            NObs= 2
          ELSEIF ( ObsType.eq.2 ) THEN
            NObs= 3
          ELSEIF ( ObsType.eq.3 ) THEN
            NObs= 4
          ELSEIF ( ObsType.eq.4 ) THEN
            NObs= 6
          ENDIF   ! Case
        ! Zero out the matricies
        DO i = 1, StateNum
            IF (i.le.Indobs) THEN
                b(i,1) = 0.0D0
              ENDIF
            ATWbacc(i,1) = 0.0D0
            ATWb(i,1) = 0.0D0
            DO j = 1, StateNum
                a(i,j) = 0.0D0
                ATWA(i,j) = 0.0D0
                ATWAacc(i,j) = 0.0D0
              ENDDO
          ENDDO

        ! -------------- Loop through all the observations ------------
        DO i= FirstOb, LastOb
            Read(15,REC=i) SensNum, JD, (RSVec(j),j=1,3),
     &                 ObsType, Rng, Az, El, DRng,
     &                 Daz, Del, TRtasc, TDecl

c           write(20,*) jd,'jd in side',jdepoch
           write(20,*) i,' rng',rng,' jd',jd
cdav fix jd diff, seems off
                       Dtsec= (JD - JDEpoch) *86400.0D0 - 0.00024140D0

            ! ---- Propagate the nominal vector to the epoch time -----
*            CALL KEPLER( RNom,VNom,Dtsec, r3,v3,error )
            CALL PKEPLER( RNom,VNom,ndot,nddot,Dtsec, r1,v1 )  ! this is in eci----
           write(*,*) sensnum,'r1',r1(1),dtsec
          
            ! ------------ Determine sensor characteristics -----------
            CALL FindSenPtr( SensNum )

cdav            CALL LSTIME( SenLon,JD,  LST,GST )

            ! ------------------ Find B matrix ------------------------
            IF ( ObsType .eq. 5 ) THEN
                CALL RV_TRADEC( R1,V1,RSVec,'TOO',
     &                     RngNom,TRtAscNom,TDeclNom,DRngNom,DTRtAscNom,
     &                     DTDeclNom )
              ELSE
         lod  =  0.0D0
         timezone=0
         order = 106
         terms = 2
         dut1 =  0.3261068D0
         dat  = 29
         xp   =  -0.11554D0* 3.1415926535879D0 / (3600.0D0*180.0D0)
         yp   =  0.48187D0* 3.1415926535879D0 / (3600.0D0*180.0D0)

         CALL InvJDAY( jd, year,month,day,hr,minute,sec )

         CALL CONVTIME ( Year, Month, Day, Hr, minute, SEC,
     &                  TimeZone, 'C', DUT1, DAT, xp, yp, UT1,
     &                  TUT1, JDUT1, UTC, TAI, TDT, TTT, JDTDT,TDB,
     &                  TTDB, JDTDB, DDPsi,DDeps,Lod, Error )

             CALL RV_RAZEL( R1,V1,senLat,senLon,senalt,TTT,jdut1,lod,
     &                      xp,yp,terms, 'TOO',
     &                      RngNom,AzNom,ElNom,DRngNom,DAzNom,DElNom )
cdav                CALL RV_RAZEL( R1,V1,RSVec,SenLat,LST,
cdav     &                  'TOO',RngNom,AzNom,ElNom,DRngNom,DAzNom,DElNom )
              ENDIF
*          write(20,*) 'r1',r1(1),r1(2),r1(3)

            IF ( ObsType .eq. 0 ) THEN
                B(1,1) = DRng - DRngNom
              ELSEIF ( ObsType .eq. 1 ) THEN
                B(1,1) = Az  - AzNom
                B(2,1) = El  - ElNom
              ELSEIF (( ObsType .ge. 2 ).and.
     &                ( ObsType .le. 4 )) THEN
                B(1,1) = Rng - RngNom
                B(2,1) = Az  - AzNom
                B(3,1) = El  - ElNom
                IF ( ObsType .gt. 2 ) THEN
                    B(4,1) = DRng - DRngNom
                    IF ( ObsType .eq. 4 ) THEN
                        B(5,1) = DAz  - DAzNom
                        B(6,1) = DEl  - DElNom
                      ENDIF
                  ENDIF
              ELSEIF ( ObsType .eq. 5 ) THEN
                B(1,1) = TRtAsc - TRtAscNom
                B(2,1) = TDecl  - TDeclNom
              ELSEIF ( ObsType .eq. 6 ) THEN
                B(1,1) = Rng - RngNom
              ENDIF  ! Case construct
*
            ! ------------------- Find A matrix -----------------------
            ! ---- Find modified obs values after propagating state ---
            DO j= 1,3
                RPert(j)= RNom(j)  ! this is in eci----
                VPert(j)= VNom(j)
              ENDDO

            DO j= 1, 6
                IF ( j .le. 3 ) THEN
                    DeltaAmt= RPert(j)*PercentChg
                    RPert(j)= RPert(j) + DeltaAmt
                  ELSE
                    DeltaAmt= VPert(j-3)*PercentChg
                    VPert(j-3)= VPert(j-3) + DeltaAmt
                  ENDIF 
  
*                CALL KEPLER( RPert,VPert,Dtsec,R3,V3,error )
                CALL PKEPLER( RPert,VPert,NDot,NDDot,Dtsec,R3,V3 )

                IF ( ObsType .eq. 5 ) THEN
                    CALL RV_TRADEC( R3,V3,RSVec,'TOO',
     &                      RngPert,TRtAscPert,TDeclPert,DRngPert,
     &                      DTRtAscPert,DTDeclPert )
                  ELSE
         CALL CONVTIME ( Year, Month, Day, Hr, minute, SEC,
     &                  TimeZone, 'C', DUT1, DAT, xp, yp, UT1,
     &                  TUT1, JDUT1, UTC, TAI, TDT, TTT, JDTDT,TDB,
     &                  TTDB, JDTDB, DDPsi,DDeps,Lod, Error )

             CALL RV_RAZEL( R3,V3,senLat,senLon,senalt,TTT,jdut1,lod,
     &               xp,yp,terms, 'TOO',
     &               Rngpert,Azpert,Elpert,DRngpert,DAzpert,DElpert )
cdav                    CALL RV_RAZEL( R3,V3,RSVec,SenLat,
cdav     &                      LST,'TOO',RngPert,AzPert,ElPert,DRngPert,
cdav     &                      DAzPert,DElPert )
                  ENDIF
                  
                  
c          write(20,'(a,8f18.8)') 'rnom ',r1(1),r1(2),r1(3), RngNom,
c     &              AzNom,ElNom ,percentchg  , DeltaAmt     
c          write(20,'(a,7f18.8)') 'rpert',r3(1),r3(2),r3(3), Rngpert,
c     &              Azpert,Elpert,dtsec
                  IF ( ObsType .eq. 0 ) THEN
                    A(1,j) = (DRngPert - DRngNom) / DeltaAmt
                  ELSEIF ( ObsType .eq. 1 ) THEN
                    A(1,j) = (AzPert  - AzNom)  / DeltaAmt
                    A(2,j) = (ElPert  - ElNom)  / DeltaAmt
                  ELSEIF (( ObsType .ge. 2 ).and.
     &                    ( ObsType .le. 4 )) THEN
                    A(1,j) = (RngPert - RngNom) / DeltaAmt
                    A(2,j) = (AzPert  - AzNom)  / DeltaAmt
                    A(3,j) = (ElPert  - ElNom)  / DeltaAmt
                    IF ( ObsType .gt. 2 ) THEN
                        A(4,j) = DRngPert - DRngNom
                        IF ( ObsType .eq. 4 ) THEN
                            A(5,j) = DAzPert  - DAzNom
                            A(6,j) = DElPert  - DElNom
                          ENDIF
                      ENDIF
                  ELSEIF ( ObsType .eq. 5 ) THEN
                    A(1,j) = (TRtAscPert  - TRtAscNom) / DeltaAmt
                    A(2,j) = (TDeclPert   - TDeclNom)  / DeltaAmt
                  ELSEIF ( ObsType .eq. 6 ) THEN
                    A(1,j) = (RngPert - RngNom) / DeltaAmt
                  ENDIF  ! Case construct

                ! --------------- Reset the modified vector -----------
                IF ( j .le. 3 ) THEN
                    RPert(j)= RNom(j) ! this is in eci----
                  ELSE
                    VPert(j-3)= VNom(j-3)
                  ENDIF
              ENDDO   ! DO j = 1 to 6

            ! ------------ Now form the matrix combinations -----------
            CALL MATTRANS  ( A,  NObs,6, 8,8,  At )
*
            ! -------------------- Assign weights ---------------------
            IF ( ObsType .eq. 0 ) THEN
                W1= 1.0D0 / (NoiseDRng*NoiseDRng)
                DRng2= DRng2 +   B(1,1)*  B(1,1)*W1
              ELSEIF ( ObsType .eq. 1 ) THEN
                W1= 1.0D0 / (NoiseAz*NoiseAz)
                W2= 1.0D0 / (NoiseEl*NoiseEl)
                DAz2 = DAz2 +   B(1,1)*  B(1,1)*W1
                DEl2 = DEl2 +   B(2,1)*  B(2,1)*W2
              ELSEIF (( ObsType .ge. 2 ).and.
     &                ( ObsType .le. 4 )) THEN
                W1= 1.0D0 / (NoiseRng*NoiseRng)
                W2= 1.0D0 / (NoiseAz*NoiseAz)
                W3= 1.0D0 / (NoiseEl*NoiseEl)
                DRng2= DRng2 +   B(1,1)*  B(1,1)*W1
                DAz2 = DAz2  +   B(2,1)*  B(2,1)*W2
                DEl2 = DEl2  +   B(3,1)*  B(3,1)*W3
                IF ( ObsType .gt. 2 ) THEN
                    W4= 1.0D0 / (NoiseDRng*NoiseDRng)
                    DDRng2= DDRng2 +   B(4,1)*  B(4,1)*W4
                    IF ( ObsType .eq. 4 ) THEN
                        W5= 1.0D0 / (NoiseDAz*NoiseDAz)
                        W6= 1.0D0 / (NoiseDEl*NoiseDEl)
                        DDAz2 = DDAz2 +   B(5,1)*  B(5,1)*W5
                        DDEl2 = DDEl2 +   B(6,1)*  B(6,1)*W6
                      ENDIF
                  ENDIF
              ELSEIF ( ObsType .eq. 5 ) THEN
                W1= 1.0D0 / (NoiseTRtAsc*NoiseTRtAsc)
                W2= 1.0D0 / (NoiseTDecl*NoiseTDecl)
                DTRtAsc2= DTRtAsc2 +   B(1,1)*  B(1,1)*W1
                DTDecl2 = DTDecl2  +   B(2,1)*  B(2,1)*W2
              ELSEIF ( ObsType .eq. 6 ) THEN
                W1= 1.0D0 / (NoiseRng*NoiseRng)
              ENDIF  ! Case construct

            DO Rowc= 1, 6
                DO Colc= 1, NObs
                    IF (Colc .eq. 1) THEN
                        Weight= W1
                      ELSEIF (Colc .eq. 2) THEN
                        Weight= W2
                      ELSEIF (Colc .eq. 3) THEN
                        Weight= W3
                      ELSEIF (Colc .eq. 4) THEN
                        Weight= W4
                      ELSEIF (Colc .eq. 5) THEN
                        Weight= W5
                      ELSEIF (Colc .eq. 6) THEN
                        Weight= W6
                      ENDIF   ! Case
                    AtW(Rowc,ColC) = At(Rowc,Colc)*Weight
                  ENDDO   ! DO Colc
              ENDDO  ! DO Rowc
 
            ! ----------------- Find the final answer -----------------
            CALL MATMULT( AtW,A, 6,NObs,6, 8,8,8  , AtWAacc )
            CALL MATMULT( AtW,B, 6,NObs,1, 8,8,1  , AtWBacc )

            ! --------------- Accumulate the matricies ----------------
            DO r=1, 6
                DO c= 1, 6
                     AtWA(r,c) = AtWAAcc(r,c)+ AtWA(r,c)
                  ENDDO
              ENDDO

            c= 1
            DO r=1, 6
                 AtWB(r,c) = AtWBAcc(r,c)+ AtWB(r,c)
              ENDDO
c        Title = ' AtWA Matrix '
c        CALL FILEEXPwriteMAT( AtWA,6,6, 8,8,Title )

          ENDDO  ! DO i through the observations

        Title = ' AtWA Matrix '
        CALL FILEEXPwriteMAT( AtWA,6,6, 8,8,Title )
        Title = ' AtWB Matrix '
        CALL FILEwriteMAT( AtWB,6,1, 8,1,Title )
        Title = ' A Matrix '
        CALL FILEEXPwriteMAT( A,6,6, 8,8,Title )
        Title = ' B Matrix '
        CALL FILEwriteMAT( B,6,1, 8,1,Title )

      RETURN
      END   ! Subroutine FindAtWAAtWB
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE LeastSquares
*
* this subroutine performs orbit determination using least squares differential
*   correction method. A variety of observation types are possible.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    PercentChg  - Amount to modify the vectors
*                  by in finite differencing
*    Epsilon     - Tolerance DO calculations
*    JDEpoch     - Julian date DO solution       Days from 4713 BC
*    TypeAns     - Type of DC (IOD .or. Nom)        'I','N'
*    FirstOb     - First obs record to use
*    LastOb      - Last obs record to use
*    ObsRecFile  - File of RECORD containing:
*                  SeNum, JD, RSVec, ObsType,
*                  Rng, Az, El, DRng, DAz, DEl,
*                  TRtasc, TDecl data
*
*  OutPuts       :
*    RNom        - Nom position vector at epoch   ER
*    VNom        - Nom velocity vector at epoch   ER/TU
*    X           - State vector                               ER, ER/TU
*    DX          - Last iteration of state corr   ER, ER/TU
*    AtWAI       - Covariance Matrix Inverse
*    AtWA        - Covariance Matrix
*
*  Locals        :
*    Loops       - Number of iterations performed
*    NumWork     - Number of observations working with
*    AtWB        - Accumulated matrix of partial derivatives
*    SigmaOld    -
*    SigmaNew    -
*    DRng2       -
*    DAz2        -
*    DEl2        -
*    DTRtAsc2    -
*    DTDecl2     -
*
*  Coupling      :
*    FindAtWAAtWB- Find combination matrices DO use in final equations
*    MATINVERSE  - Matrix inverse routine
*    CALL MATMULT     - Matrix multiply routine
*
*  References    :
*    Vallado       1997, 682-702
*
* -----------------------------------------------------------------------------

      SUBROUTINE LeastSquares ( PercentChg,Epsilon,JDEpoch, TypeAns,
     &                          FirstOb,LastOb,
     &                          Rnom,Vnom,NDot,NDDot,
     &                          X,DX,AtWAI,AtWA,atwb )
        IMPLICIT NONE
        INTEGER StateNum
        PARAMETER (StateNum = 8)
        REAL*8 PercentChg,Epsilon,JDEpoch, rnom(3),vnom(3),NDot,NDDot,
     &         X(StateNum,1),DX(StateNum,1),AtWAI(StateNum,StateNum),
     &         AtWA(StateNum,StateNum), AtWB(StateNum,1)
        CHARACTER TypeAns
        INTEGER FirstOb,LastOb

* ----------------------------  Locals  -------------------------------
        INTEGER i,j,loops,NumWork
        REAL*8 SigmaNew,SigmaOld, DRng2,DAz2,DEl2,
     &         DDRng2,DDAz2,DDel2,DTRtAsc2,DTDecl2

        COMMON /ObsRec/
     &         JD, RSVec, ObsType, SensNum,
     &         Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl
        INTEGER ObsType, SensNum
        REAL*8 JD, RSVec(3), Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl
        CHARACTER*64 Title

        INCLUDE 'astconst.cmn'
        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Write(20,*) 'Input data observations are: '
        DO i= FirstOb, LastOb
            Read(15,REC=i) SensNum, JD, (RSVec(j),j=1,3),
     &                 ObsType, Rng, Az, El, DRng,
     &                 Daz, Del, TRtasc, TDecl
            Write(20,*) i,' ',SensNum,JD,(RSVec(j),j=1,3),
     &                 ObsType, Rng,'~',Az*rad2deg,
     &                 El*rad2deg,TRtasc*rad2deg,TDecl*rad2deg
          ENDDO

        ! --------------------- Initialize parameters -----------------
        SigmaNew= 100.0D0 
        SigmaOld=   1.0D0  ! initialize to let work for repeat until
        NumWork = LastOb - FirstOb + 2
        IF ( TypeAns .eq. 'I' ) THEN
            NumWork = NumWork - 2   ! 2 obs produce only position vector
          ENDIF
*
        Loops= 0 
        DO WHILE ( (DABS((SigmaNew-SigmaOld)/SigmaOld) .ge. Epsilon)
     &         .and. (Loops .lt. 10) )
            SigmaOld= SigmaNew 
            DRng2   = 0.0D0
            DAz2    = 0.0D0
            DEl2    = 0.0D0
            DDRng2  = 0.0D0
            DDAz2   = 0.0D0
            DDEl2   = 0.0D0
            DTRtAsc2= 0.0D0
            DTDecl2 = 0.0D0

            ! --- Place the nominal value into the state vector ----
            DO i= 1, 3
                X(i,1)  = rnom(i)  ! this is in eci----
                X(i+3,1)= vnom(i)
              ENDDO

            CALL FindAtWAAtWB( FirstOb,LastOb,JDEpoch,PercentChg,RNom,
     &                       VNom,NDot,NDdot,DRng2,DAz2,DEl2,DDRng2,
     &                       DDAz2,DDel2,DTRtAsc2,DTDecl2,AtWA,AtWB )

            CALL MATINVERSE( AtWA      , 6,8         , AtWAI )
            CALL MATMULT   ( AtWAI,AtWB, 6,6,1, 8,8,1   , DX   )

            ! --- Update state vector ----
            DO i= 1, 3
                RNom(i) = rNom(i) + DX(i,1)  ! this is in eci----
                VNom(i) = vNom(i) + DX(i+3,1)
                X(i,1)  = rnom(i)
                X(i+3,1)= vnom(i)
              ENDDO

            IF ( ObsType .eq. 0 ) THEN
                SigmaNew= DSQRT( DDRng2 / NumWork )
              ELSEIF ( ObsType .eq. 1 ) THEN
                SigmaNew= DSQRT( (DAz2 + DEl2) / NumWork )
              ELSEIF ( ObsType .eq. 2 ) THEN
                SigmaNew= DSQRT( (DRng2 + DAz2 + DEl2) / NumWork )
              ELSEIF ( ObsType .eq. 3 ) THEN
                SigmaNew= DSQRT( (DRng2 + DAz2 + DEl2 + DDRng2) /
     &                    NumWork )
              ELSEIF ( ObsType .eq. 4 ) THEN
                SigmaNew= DSQRT( (DRng2 + DAz2 + DEl2 + DDRng2 +
     &                    DDAz2 + DDEl2) / NumWork )
              ELSEIF ( ObsType .eq. 5 ) THEN
                SigmaNew= DSQRT( (DTRtAsc2 + DTDecl2) / NumWork )
              ELSEIF ( ObsType .eq. 6 ) THEN
                SigmaNew= DSQRT( DRng2 / NumWork )
              ENDIF

            ! ------------------- Write out Data ------------------------
*            Write(*,*)  FileOut,'Perturb Delta % ',PercentChg:6:4,
*    &       ' Epsilon for convergence ',epsilon:6:4 )

            Title = ' AtWAi Matrix '
            CALL FILEwriteMAT( AtWAi,6,6, 8,8,Title )
            CALL writeMAT( AtWAi,6,6, 8,8,Title )
            Write(*,'(A,6f14.6)')  'Corrections    ',DX(1,1),DX(2,1),
     &               DX(3,1),DX(4,1),DX(5,1),DX(6,1)
            Write(*,'(A,I2,A,6f12.6)') 'Loop ',Loops,' Nominal',RNom(1),
     &               RNom(2),RNom(3),VNom(1),VNom(2),VNom(3)
            Write(*,*)
            Write(20,'(A,6f14.6)')  'Corrections    ',DX(1,1),DX(2,1),
     &               DX(3,1),DX(4,1),DX(5,1),DX(6,1)
            Write(20,'(A,I2,A,6f12.6)') 'Loop ',Loops,' Nominal',
     &               RNom(1),RNom(2),RNom(3),VNom(1),VNom(2),VNom(3)
            Write(20,*)
            Write(20,'(A,I2,A,3f14.7,3f12.6)') 'Loop ',Loops,' Nominal',
     &            RNom(1)/rekm,RNom(2)/rekm,
     &            RNom(3)/rekm,VNom(1)/7.905366149D0,
     &            VNom(2)/7.905366149D0,VNom(3)/7.905366149D0
            Write(*,'(A,I2,A,3f14.7,3f12.8)') 'Loop ',Loops,' Nominal',
     &            RNom(1)/rekm,RNom(2)/rekm,
     &            RNom(3)/rekm,VNom(1)/7.905366149D0,
     &            VNom(2)/7.905366149D0,VNom(3)/7.905366149D0
            Write(*,'(A,f12.6,A,f12.6)')  'SigmaOld ',sigmaold,
     &               ' SigmaNew ',SigmaNew
            Write(20,'(A,f12.6,A,f12.6)')  'SigmaOld ',sigmaold,
     &              ' SigmaNew ',SigmaNew
            Write(20,'(A,6f12.6)')  'r v sigmas ',dsqrt(atwai(1,1)),
     &              dsqrt(atwai(2,2)),dsqrt(atwai(3,3)),
     &              dsqrt(atwai(4,4)),dsqrt(atwai(5,5)),
     &              dsqrt(atwai(6,6))
            Write(*,*)  ' -------------------------------------------'
            Write(20,*)  ' -------------------------------------------'
            Write(*,*)  'k for loop ',loops
            Write(20,*)  'k for loop ',loops

            Loops= Loops + 1
          ENDDO

      RETURN
      END   ! Subroutine LeastSquares
*
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE Sequential
*
* this subroutine calculates the sequential least squares differential
*   correction. A variety of observation types are posible. Remember
*   that we must run a least squares technique to obtain the starting values.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    PercentChg  - Amount to modify the vectors
*                  by in finite differencing
*    ObsRecFile    - Array of RECORDs containing:
*                  SeNum, JD, RSVec, ObsType,
*                  Rng, Az, El, DRng, DAz, DEl,
*                  TRtasc, TDecl data
*    Epsilon     - Tolerance DO calculations
*    JDEpoch     - Julian date DO solution       Days from 4713 BC
*    UseObs      - Number of observations to use
*                  ( a subset of the array )
*    dxOld       - Satrting iteration of the state corrcetions    ER, ER/TU
*    AtWAOld     - Starting Covariance Matrix
*
*  OutPuts       :
*    RNom        - Nom position vector at epoch   ER
*    VNom        - Nom velocity vector at epoch   ER/TU
*    X           - State vector                               ER, ER/TU
*    X           - State vector                               ER, ER/TU
*    AtWA        - Covariance Matrix
*
*  Locals        :
*    Loops       - Number of iterations performed
*    NumWork     - Number of observations working with
*    A           - A matrix of partial derivatives
*    b           - b matrix of observation residuals
*
*
*
*  Coupling      :
*    MATINVERSE  - Matrix inverse routine
*    MATMULT     - Matrix multiply routine
*    MATADD      - Matrix add routine
*    FindAtWAAtWB- Find combination matrices DO use in final equations
*
*  References    :
*    Vallado       1997, 702-707
*
* -----------------------------------------------------------------------------
*
      SUBROUTINE Sequential    ( PercentChg,Epsilon,JDEpoch, FirstOb,
     &                           LastOb,
     &                           Rnom,Vnom,
     &                           AtWAold,atwbold,
     &                           X,DX,pbest )
        IMPLICIT NONE
        INTEGER StateNum
        PARAMETER (StateNum = 8)
        REAL*8 PercentChg,Epsilon,JDEpoch, rnom(3),vnom(3),NDot,NDDot,
     &         X(StateNum,1),DX(StateNum,1),PBest(StateNum,StateNum),
     &         AtWAOld(StateNum,StateNum), DXOld(Statenum,1),
     &         Atwbold(statenum,1)
        CHARACTER TypeAns
        INTEGER FirstOb,LastOb

* ----------------------------  Locals  -------------------------------
        INTEGER i,j,loops,NumWork
        REAL*8 PDXOld(Statenum,1), AtWB(StateNum,1),
     &         SigmaNew,SigmaOld, DRng2,DAz2,DEl2, AtWBP(Statenum,1),
     &         DDRng2,DDAz2,DDel2,DTRtAsc2,DTDecl2,Atwap(Statenum,1),
     &         ATWA(Statenum,Statenum)

        COMMON /ObsRec/
     &         JD, RSVec, ObsType, SensNum,
     &         Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl
        INTEGER ObsType, SensNum
        REAL*8 JD, RSVec(3), Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl
        CHARACTER*64 Title

        INCLUDE 'astconst.cmn'
        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Write(20,*) 'Input data observations are: '
        DO i= FirstOb, LastOb
            Read(15,REC=i) SensNum, JD, (RSVec(j),j=1,3),
     &                 ObsType, Rng, Az, El, DRng,
     &                 Daz, Del, TRtasc, TDecl
            Write(20,*) i,' ',SensNum,JD,(RSVec(j),j=1,3),
     &                 ObsType, Rng,'~',Az*rad2deg,
     &                 El*rad2deg,TRtasc*rad2deg,TDecl*rad2deg
          ENDDO

        SigmaNew= 100.0D0
        SigmaOld=   1.0D0  ! initialize to let work for repeat until
        NumWork = LastOb - FirstOb + 1

        Loops= 0
        DO WHILE ( (DABS((SigmaNew-SigmaOld)/SigmaOld) .ge. Epsilon)
     &         .and. (Loops .lt. 10) )
            SigmaOld= SigmaNew 
            DRng2   = 0.0D0
            DAz2    = 0.0D0
            DEl2    = 0.0D0
            DDRng2  = 0.0D0
            DDAz2   = 0.0D0
            DDEl2   = 0.0D0
            DTRtAsc2= 0.0D0
            DTDecl2 = 0.0D0

            ! --- Place the nominal value into the state vector ----
            DO i= 1, 3
                X(i,1)  = rnom(i)  ! this is in eci----
                X(i+3,1)= vnom(i)
              ENDDO

           ndot = 0.0D0
           nddot=0.0d0

            CALL FindAtWAAtWB( FirstOb,LastOb,JDEpoch,PercentChg,RNom,
     &                       VNom,NDot,NDdot,DRng2,DAz2,DEl2,DDRng2,
     &                       DDAz2,DDel2,DTRtAsc2,DTDecl2,AtWA,AtWB )

            Title = ' AtWA  Matrix '
            CALL FILEwriteMAT( AtWA,6,6, 8,8,Title )
            Title = ' AtWB  Matrix '
            CALL FILEwriteMAT( AtWB,6,1, 8,1,Title )

            ! --- Implement Sequential equations ----
            CALL MATADD ( AtWB,AtwbOld,  6,1, 8,1,  AtWBP )

            CALL MATADD ( AtWA,AtWAold,  6,6, 8,8,  AtWAP )

            CALL MATINVERSE( AtWAp,      6,8,     PBest   )    ! P = (PnewI + PoldI)I

            CALL MATMULT( PBest,AtWBP, 6,6,1, 8,8,1, DX      )    ! DX = DX +

            ! --- Update state vector ----
            DO i= 1, 3
                RNom(i) = rNom(i) + DX(i,1)  ! this is in eci----
                VNom(i) = vNom(i) + DX(i+3,1)
                X(i,1)  = rnom(i)
                X(i+3,1)= vnom(i)
              ENDDO

            IF ( ObsType .eq. 0 ) THEN
                SigmaNew= DSQRT( DDRng2 / NumWork )
              ELSEIF ( ObsType .eq. 1 ) THEN
                SigmaNew= DSQRT( (DAz2 + DEl2) / NumWork )
              ELSEIF ( ObsType .eq. 2 ) THEN
                SigmaNew= DSQRT( (DRng2 + DAz2 + DEl2) / NumWork )
              ELSEIF ( ObsType .eq. 3 ) THEN
                SigmaNew= DSQRT( (DRng2 + DAz2 + DEl2 + DDRng2) /
     &                    NumWork )
              ELSEIF ( ObsType .eq. 4 ) THEN
                SigmaNew= DSQRT( (DRng2 + DAz2 + DEl2 + DDRng2 +
     &                    DDAz2 + DDEl2) / NumWork )
              ELSEIF ( ObsType .eq. 5 ) THEN
                SigmaNew= DSQRT( (DTRtAsc2 + DTDecl2) / NumWork )
              ELSEIF ( ObsType .eq. 6 ) THEN
                SigmaNew= DSQRT( DRng2 / NumWork )
              ENDIF

            ! ------------------- Write out Data ------------------------
*            Write(*,*)  FileOut,'Perturb Delta % ',PercentChg:6:4,
*    &       ' Epsilon for convergence ',epsilon:6:4 )

            Title = ' PBest Matrix '
            CALL FILEwriteMAT( PBest,6,6, 8,8,Title )
            CALL writeMAT( PBest,6,6, 8,8,Title )
            Write(*,'(A,6f14.6)')  'Corrections    ',DX(1,1),DX(2,1),
     &               DX(3,1),DX(4,1),DX(5,1),DX(6,1)
c            Write(*,'(A,I2,A,6f12.6)') 'Loop ',Loops,' Nominal',RNom(1),
            Write(*,*) 'Loop ',Loops,' Nominal',RNom(1),
     &               RNom(2),RNom(3),VNom(1),VNom(2),VNom(3)
c            Write(*,*)
c            Write(20,'(A,6f14.6)')  'Corrections    ',DX(1,1),DX(2,1),
            Write(20,*)  'Corrections    ',DX(1,1),DX(2,1),
     &               DX(3,1),DX(4,1),DX(5,1),DX(6,1)
c            Write(20,'(A,I2,A,6f12.6)') 'Loop ',Loops,' Nominal',
            Write(20,*) 'Loop ',Loops,' Nominal',
     &               RNom(1),RNom(2),RNom(3),VNom(1),VNom(2),VNom(3)
            Write(20,*)  'Corrections    ',DX(1,1),DX(2,1),
     &               DX(3,1),DX(4,1),DX(5,1),DX(6,1)
            Write(20,*) 'Loop ',Loops,' Nominal',
     &               RNom(1),RNom(2),RNom(3),VNom(1),VNom(2),VNom(3)
c            Write(20,*)
c            Write(20,'(A,I2,A,3f14.7,3f12.6)') 'Loop ',Loops,' Nominal',
c     &            RNom(1)*rekm,RNom(2)*rekm,
c     &            RNom(3)*rekm,VNom(1)*7.905366149D0,
c     &            VNom(2)*7.905366149D0,VNom(3)*7.905366149D0
c            Write(*,'(A,I2,A,3f14.7,3f12.8)') 'Loop ',Loops,' Nominal',
c     &            RNom(1)*rekm,RNom(2)*rekm,
c     &            RNom(3)*rekm,VNom(1)*7.905366149D0,
c     &            VNom(2)*7.905366149D0,VNom(3)*7.905366149D0
c            Write(*,11)  'SigmaOld ',sigmaold,
c     &               ' SigmaNew ',SigmaNew
c            Write(20,'(A,f12.6,A,f12.6)')  'SigmaOld ',sigmaold,
c     &              ' SigmaNew ',SigmaNew
c            Write(20,'(A,6f12.6)')  'r v sigmas ',dsqrt(atwai(1,1)),
c            Write(20,'(A,6f12.6)')  'r v sigmas ',dsqrt(Pbest(1,1)),
            Write(20,*)  'r v sigmas ',dsqrt(Pbest(1,1)),
     &              dsqrt(pbest(2,2)),dsqrt(pbest(3,3)),
     &              dsqrt(pbest(4,4)),dsqrt(pbest(5,5)),
     &              dsqrt(pbest(6,6))
            Write(*,*)  ' -------------------------------------------'
            Write(20,*)  ' -------------------------------------------'
            Write(*,*)  'k for loop ',loops
            Write(20,*)  'k for loop ',loops

            Loops= Loops + 1
          ENDDO

 11     FORMAT(A,f12.6,A,f12.6)

      RETURN
      END   ! Subroutine Sequential
*
** -----------------------------------------------------------------------------
**
**                           SUBROUTINE FINDPHI
**
** this subroutine finds the Phi matrix via numerical techniques.
**
**
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RNom        - Nom position vector at epoch   ER
*    VNom        - Nom velocity vector at epoch   ER/TU
*    PercentChg  - Amount to modify the vectors
*                  by in finite differencing
*
*
*
*  OutPuts       :
*    A           - A matrix
*
*
*  Locals        :
*    z23
*
*
*
*
*
*
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    MATMULT     - Matrix Multiplication
*    MATADD      - Matrix Addition
*    MATSCALE    - Matrix Multiplication by a constant
*
*  References    :
*    Vallado       1997, 689-695, 709-712
*
* -----------------------------------------------------------------------------

*      SUBROUTINE FindPhi       ( Rnom,Vnom,Dtsec,PercentChg,Phi)
*        IMPLICIT NONE
*        REAL*8 Rnom(3),Vnom(3), Dtsec,PercentChg, Phi()
** ----------------------------  Locals  -------------------------------
        ! --------------------  Implementation   ----------------------
*
*        REAL*8 RPert(3),VPert(3),VDot(3),VPDot(3), DeltaAmt,RCubed
*        F,IM,Fdt,Fdt2,IFdt,tFdt : Matrix
*        INTEGER i,j
*         ! -- local vars ---
*        INITMATRIX( 6,6,F )
*        INITMATRIX( 6,6,IM )
*          IM(1,1) = 1.0D0
*          IM(2,2) = 1.0D0
*          IM(3,3) = 1.0D0
*          IM(4,4) = 1.0D0
*          IM(5,5) = 1.0D0
*          IM(6,6) = 1.0D0
*
*        ! -------- Find state derivative from nominal vector ----------
*        RCubed= Rnom(3)*Rnom(3)*Rnom(3)
*        VDot(1)= -Rnom(1)/RCubed
*        VDot(2)= -Rnom(2)/RCubed
*        VDot(3)= -Rnom(3)/RCubed
*
*        DO i =1,3
*            RPert(i)= RNom(i)
*            VPert(i)= VNom(i)
*          ENDDO
*
*        ! -------- Find state derivative from perturbed vectors -------
*        DO j= 1, 6
*            IF ( j .le. 3 ) THEN
*                DeltaAmt= RPert(j)*PercentChg
*                RPert(j)= RPert(j) + DeltaAmt
*              ELSE
*                DeltaAmt= VPert(j-3)*PercentChg
*                VPert(j-3)= VPert(j-3) + DeltaAmt
*              ENDIF
*            CALL MAG(RPert)
*
*            RCubed= Rpert(3)*Rpert(3)*Rpert(3)
*            VPDot(1)= -Rpert(1)/RCubed
*            VPDot(2)= -Rpert(2)/RCubed
*            VPDot(3)= -Rpert(3)/RCubed
*
*            F(1,j)= (VPert(1)-VNom(1)) / DeltaAmt )
*            F(2,j)= (VPert(2)-VNom(2)) / DeltaAmt )
*            F(3,j)= (VPert(3)-VNom(3)) / DeltaAmt )
*            F(4,j)= (VPDot(1)-VDot(1)) / DeltaAmt )
*            F(5,j)= (VPDot(2)-VDot(2)) / DeltaAmt )
*            F(6,j)= (VPDot(3)-VDot(3)) / DeltaAmt )
*
*            ! --- Reset perturbed element ----
*            IF ( j .le. 3 ) THEN
*                RPert(j)= RNom(j)
*              ELSE
*                VPert(j-3)= VNom(j-3)
*              ENDIF
*          ENDDO
*
*        ! ---------------------- Find Phi matrix ----------------------
*        CALL MATSCALE( F,Dtsec,  Fdt  )
*        CALL MATADD( IM,FDt,6,6, IFdt )
*
*        CALL MATMULT( Fdt,Fdt,6,6,6, tFdt )
*        CALL MATSCALE( tFdt,0.25D0, Fdt2 )
*        CALL MATADD( IFdt,Fdt2,6,6, Phi )
*
*        Write(*,*)  'Done with Phi ' )
*      RETURN
*      END   ! Subroutine FindPhi
**
** -----------------------------------------------------------------------------
**
**                           SUBROUTINE FINDH
**
** this subroutine finds the H matrix via numerical techniques.
**     Partial of obs at ent time
**     LST is passed back DO b matrix calculations
**
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
**
**  Inputs          Description                    Range / Units
**    RNom        - Nom position vector at epoch   ER
**    VNom        - Nom velocity vector at epoch   ER/TU
**    PercentChg  - Amount to modify the vectors
**                  by in finite differencing
**
**
**
**  OutPuts       :
**    A           - A matrix
**    LST         - Local SIDEREAL Time
**
**  Locals        :
**    RNom
**
**
**
**
**
**
**
**
*
*
*
*  Coupling      :
*    KEPLER      - Propagation
*    RAZEL       -
*    MAG         -
*    LSTIME      -
*
*  References    :
*    Vallado       1997, 688-689, 714
*
* -----------------------------------------------------------------------------
*
*      SUBROUTINE FINDH         ( XBest, Rec, ObsRecFile,
*     &                           PercentChg,Dtsec,SenPtr,LST,H )
*        IMPLICIT NONE
*        REAL*8 XBest(), PercentChg,Dtsec, LST, H()
*        Rec                            : ObsRec
*        VAR ObsRecFile                        : ObsRecType
*        SenPtr                            : SenPtr
** ----------------------------  Locals  -------------------------------
*        RS,RNom,VNom,RPert,VPert,RNew,VNew : Vector
*        RngNom,AzNom,ElNom,GST,DeltaAmt,RngPert,AzPert,ElPert,DRng,DAz,DEl : REAL*8
*        Loops,i,j      : INTEGER
*        Error:Str12
*
        ! --------------------  Implementation   ----------------------
*        write(*,*) 'insideh'
*        INITMATRIX( 3,6,H )
*        DO i= 1, 3
*            RNom(i)=   XBest(i,1 )
*            VNom(i)=   XBest(i+3,1 )
*          ENDDO
*        CALL MAG( RNom )
*        CALL MAG( VNom )
*
*        ! ---------------------- Get nominal values -------------------
*        CALL LSTIME( SenPtr^.SenLon,JD,  LST,GST )
*
*        CALL KEPLER( RNom,VNom,Dtsec, rnew,vnew,error )
*
*
**Write(*,*)  rsvec(1):14:7 )
**readln
*        CALL RV_RAZEL( Rnew,Vnew,RSvec,SenPtr^.SenLat,LST,'TOO',
*                  RngNom,AzNom,ElNom,DRng,DAz,DEl )
*        IF ( ElNom .gt. Pi ) THEN
*            ElNom= ElNom - 2.0D0*Pi
*          ENDIF
*
*        DO i =1,3
*            RPert(i)= RNom(i)
*            VPert(i)= VNom(i)
*          ENDDO
*
*      ! ---------------------- Evaluate H at the time ----------------
*        DO j= 1, 6
*            IF ( j .le. 3 ) THEN
*                DeltaAmt= RPert(j)*PercentChg
*                RPert(j)= RPert(j) + DeltaAmt
*              ELSE
*                DeltaAmt= VPert(j-3)*PercentChg
*                VPert(j-3)= VPert(j-3) + DeltaAmt
*              ENDIF
*
**         CASE xx
*            CALL KEPLER( RPert,VPert,Dtsec, rnew,vnew,error )
*          'P' : PKEPLER()
*             ENDIF    Case
*
*            CALL RV_RAZEL( RNew,VNew,RSVec,SenPtr^.Senlat,LST,
*     &                    'TOO',RngPert,AzPert,ElPert,DRng,DAz,DEl )
*            IF ( ElPert .gt. Pi ) THEN
*                ElPert= ElPert - 2.0D0*Pi
*              ENDIF
*
*            H(1,j) = (RngPert - RngNom) / DeltaAmt
*            H(2,j) = (AzPert  - AzNom)  / DeltaAmt
*            H(3,j) = (ElPert  - ElNom)  / DeltaAmt
*
*            IF ( j .le. 3 ) THEN
*                RPert(j)= RNom(j)
*              ELSE
*                VPert(j-3)= VNom(j-3)
*              ENDIF
*          ENDDO
*
*      RETURN
*      END   ! Subroutine FindH
**
* -----------------------------------------------------------------------------
*
*                           SUBROUTINE LKF
*
* this subroutine completes the linearized Kalman Filter differential
*   correction.  Note that a prior batch version must be run to initiate the
*   process.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RNom        - Nom position vector at epoch   ER
*    VNom        - Nom velocity vector at epoch   ER/TU
*    PercentChg  - Amount to modify the vectors
*                  by in finite differencing
*    ObsRecFile    - Array of RECORDs containing:
*                  SeNum, JD, RSVec, ObsType,
*                  Rng, Az, El, DRng, DAz, DEl,
*                  TRtasc, TDecl data
*    Epsilon     - Tolerance DO calculations
*    JDEpoch     - Julian date DO solution       Days from 4713 BC
*    UseObs      - Number of observations to use
*                  ( a subset of the array )
*
*  OutPuts       :
*    A           - A matrix
*
*
*  Locals        :
*    z23
*
*
*
*
*  Coupling      :
*    KEPLER      - Propagation
*
*  References    :
*    Vallado       1997, 722-724
*
* -----------------------------------------------------------------------------

*      SUBROUTINE LKF            ( PercentChg,Epsilon, JDEpoch,Q,R,
*     &                            FirstOb,LastOb, ObsRecFile, XOld,
*     &                            rnom,vnom, Xbest,dxbest,pbest )
*        IMPLICIT NONE
*        REAL*8 PercentChg,Epsilon, JDEpoch, Q(),R(), XOld(),
*     &         rnom(3),vnom(3), Xbest(),dxbest(),pbest()
*        INTEGER FirstOb,LastOb
*        VAR ObsRecFile                        : ObsRecType
** ----------------------------  Locals  -------------------------------
*
*        REAL*8 SigmaOld,SigmaNew, Drng2,Daz2,DEl2, rnew(3), vnew(3)
*        Rec : ObsRec
*        INTEGER  NumObs,  Loops, i,j
*
*        REAL*8 XPred,PPred,H,Ht,phi,DxPred,Phit,Ppt,PPPt,Pht,hpht,
*     &   a,b,AtWA,AtWB, hphtr,hphtri,Kg,hdx,bhx,kbhx,HP,KHP : Matrix
*
*        CHARACTER*12 Str12
*        REAL*8 Dtsec, tuday,JDateOld,JDate,LST,RngPred,Azpred,ElPred,
*     &         Drng,Daz,Del : REAL*8
*        First : BOOLEAN
*
        ! --------------------  Implementation   ----------------------
*        tuday= 0.0093380913806D0
*
*write('inside lkf')
*        SigmaNew= 100.0D0
*        INITMATRIX( (LastOb-FirstOb)*3, 6,A )
*        INITMATRIX( (LastOb-FirstOb)*3, 1,B )
*write('inside lkf')
*        INITMATRIX( 6, 1, XPred )
*        INITMATRIX( 3,3,R )
**       R(1,1) = SigmaRho*SigmaRho )
*          R,2,2, SigmaAz*SigmaAz )
*          R,3,3, SigmaEl*SigmaEl )
**      default settings
*          R(1,1) = 1.0D0 )
*          R,2,2, 1.0D0 )
*          R,3,3, 1.0D0 )
*
*        INITMATRIX( 3, 6, H )
*        INITMATRIX( 6, 6, PBest )
*write('after ')
*        INITMATRIX( 6, 6, Q )
*          Q(1,1) = 0.001D0  )
*          Q,2,2, 0.001D0 )
*          Q,3,3, 0.001D0 )
*          Q,4,4, 0.00001D0  )
*          Q,5,5, 0.00001D0 )
*          Q,6,6, 0.00001D0 )
*        writeMAT( Q,5,' Q Matrix ' )
**   EnterMatrix( Q,4,'Q Matrix ' )
*        writeMAT( Q,5,' Q Matrix ' )
*
*        Say( 1,27,'  Time     Residual        Delta X     Covariance  ' )
**
*        ! ---------- Loop through new data one at a time --------------
*        First   = TRUE
*
*        DO Loops= FirstOb, LastOb
*
*            Rec= Rec
*Write(*,*)  sennum:4,jd:20:12 )
*            Dtsec= (JD - JDEpoch) / TUDay
*
*            FindSenPtr( SenNum,SHead,STail,SenPtr )
*            FindPhi( RNom,VNom,Dtsec,PercentChg,  Phi )
*
*            ! --------------- Prediction ----------------
*            DO i= 1, 3
*                RNom(i)=   XOld(i,1 )     ! LKF use old value of a
*                VNom(i)=   XOld(i+3,1 )
*              ENDDO
*            CALL MAG( RNom )
*            CALL MAG( VNom )
*
*            CALL KEPLER( RNom,VNom,Dtsec,  RNew,VNew,error )
*
*            DO i= 1, 3
*                XPred(i,1)  = RNew(i)
*                XPred(i+3,1)= VNew(i)
*              ENDDO
*
*            CALL MATMULT   ( Phi,DXBest,6,6,1, DXPred )  ! DX = ÌDX
*
*            CALL MATTRANS  ( Phi       ,6,6  , Phit   )
*            CALL MATMULT   ( PBest,Phit,6,6,6, PPt    )
*            CALL MATMULT   ( Phi,PPt   ,6,6,6, PPPt   )
*            CALL MATADD    ( PPPt,Q    ,6,6,   PPred  )  ! P = ÌPÌt + Q
*
*            CALL FindH( XOld,Rec,ObsRecFile,PercentChg,Dtsec,SenPtr,  LST,H )
**    xpred .or. rnom(xold)??
*
*               ! --------------- Correction ----------------
*            CALL MATTRANS  ( H        ,3,6  , Ht    )
*            CALL MATMULT   ( PPred,Ht ,6,6,3, PHt   )
*            CALL MATMULT   ( H,PHt    ,3,6,3, HPHt  )
*            CALL MATADD    ( HPHt,R   ,3,3,   HPHtR )
*            CALL MATINVERSE( HPHtR    ,3    , HPHtRI)
*            CALL MATMULT   ( PHt,HPHtRI,6,3,3, KG   )     ! K = PHt( HPHt + R )-1
*
*            CALL MATMULT   ( H,DXPred ,3,6,1, HdX   )
*
*            CALL RV_RAZEL( RNew,VNew,RSVec,SenPtr^.SenLat,LST,'TOO',
*                      RngPred,AzPred,ElPred,DRng,DAz,DEl )
*
*            B(1,1)= Rng - RngPred )
*            B(2,1)= Az  - AzPred )
*            B(3,1)= El  - ElPred )
*
*            CALL MATSUB    ( B,Hdx    ,3,1  , bHX   )
*            CALL MATMULT   ( KG,bHX   ,6,3,1, KbHX  )
*            CALL MATADD    ( DXPred,KbHX, 6,1,DXBest)     ! DX = DX + K(b - Hdx)
*
*            CALL MATMULT   ( H,PPred  ,3,6,6 ,HP    )
*            CALL MATMULT   ( KG,HP    ,6,6,6, KHP   )
*            CALL MATSUB    ( PPred,KHP,  6,6, PBest )     ! P = P - KHP
*
**}!}!}       CALL MATADD    ( XPred,DXBest,6,1, XBest )    ! X = X + DX
*
*            Write(*,*)  ' ---- Predictions ---- ' )
*            writeMAT( XPred,5,' X Pred - State Vector ' )
*            FILEwriteMAT( xpred,' xpred Matrix ',1.0D0,6,FileOut )
*            writeMAT( dXPred,7,' DX Pred ' )
*            FILEwriteMAT( dxpred,' dxpred Matrix ',1.0D0,6,FileOut )
*            writeMAT( PPred,5,' P Pred - Pred COVARIANCE Matrix ' )
*            FILEwriteMAT( ppred,' ppred Matrix ',1.0D0,6,FileOut )
*
*            Write(*,*)  ' ---- Corrections ---- ' )
*            writeMAT( KG,5,   ' K      - Kalman Gain Matrix ' )
*            FILEwriteMAT( kg,' kg Matrix ',1.0D0,6,FileOut )
*
*            writeMAT( XBest,5,' X Best - State Vector ' )
*            FILEwriteMAT( xBest,' xBest Matrix ',1.0D0,6,FileOut )
*            writeMAT( dXBest,7,' DX Best ' )
*            FILEwriteMAT( dxBest,' dxBest Matrix ',1.0D0,6,FileOut )
*            writeMAT( PBest,5,' P Best - Best COVARIANCE Matrix ' )
*            FILEwriteMAT( pBest,' pBest Matrix ',1.0D0,6,FileOut )
*            Readln
*
*           ENDIF   ! DO WHILE not EOF
*
*        DO i= 1, 3
*            RNom(i)= XBest(i,1 )
*            VNom(i)= XBest(i+3,1 )
*          ENDDO
*        CALL MAG( RNom )
*        CALL MAG( VNom )
*      RETURN
*      END   ! Subroutine LKF
**
** -----------------------------------------------------------------------------
**
**                           SUBROUTINE EKF
**
** this subroutine completes the REAL*8 Kalman Filter differential
**   correction.  Note that a prior batch version must be run to initiate the
**   process.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RNom        - Nom position vector at epoch   ER
*    VNom        - Nom velocity vector at epoch   ER/TU
*    PercentChg  - Amount to modify the vectors
*                  by in finite differencing
*    ObsRecFile    - Array of RECORDs containing:
*                  SeNum, JD, RSVec, ObsType,
*                  Rng, Az, El, DRng, DAz, DEl,
*                  TRtasc, TDecl data
*    Epsilon     - Tolerance DO calculations
*    JDEpoch     - Julian date DO solution       Days from 4713 BC
*    UseObs      - Number of observations to use
*                  ( a subset of the array )
*
*  OutPuts       :
*    A           - A matrix
*
*
*  Locals        :
*    z23
*
*
*
*
*  Coupling      :
*    KEPLER      - Propagation
*
*  References    :
*    Vallado       1997, 724-727
*
* -----------------------------------------------------------------------------

*      SUBROUTINE EKF             ( PercentChg,Epsilon, JDEpoch,Q,R,
*     &                             FirstOb,LastOb,ObsRecFile,rnom,
*     &                             vnom,Xbest,dxbest,pbest )
*        IMPLICIT NONE
*        REAL*8 PercentChg,Epsilon, JDEpoch, Q(), R(), rnom(3),vnom(3),
*     &         Xbest(),dxbest(),pbest()
*        INTEGER FirstOb,LastOb
*        VAR ObsRecFile                        : ObsRecType
** ----------------------------  Locals  -------------------------------
*        REAL*8 SigmaOld,SigmaNew, Drng2,Daz2,DEl2, rnew(3),vnew(3)
*        Rec : ObsRec
*        INTEGER NumObs, Loops, i,j
*
*        REAL*8 xold,XPred,PPred,H,Ht,phi,DxPred,Phit,Ppt,PPPt,Pht,hpht,
*        a,b,AtWA,AtWB, hphtr,hphtri,Kg,kb,HP,KHP : Matrix
*
*        REAL*8 Dtsec,tuday,JDateOld,JDate,LST,RngPred,Azpred,ElPred,Drng,Daz,Del : REAL*8
*        error:str12
*        First : BOOLEAN
*
        ! --------------------  Implementation   ----------------------
*        tuday= 0.0093380913806D0
*
*        SigmaNew= 100.0D0
*        INITMATRIX( (LastOb-FirstOb)*3, 6,A )
*        INITMATRIX( (LastOb-FirstOb)*3, 1,B )
*
*        INITMATRIX( 6, 1, XOld )
*        XOld= XBest
*
*        INITMATRIX( 6, 1, XPred )
*        INITMATRIX( 3,3,R )
**       R(1,1) = SigmaRho*SigmaRho )
*          R,2,2, SigmaAz*SigmaAz )
*          R,3,3, SigmaEl*SigmaEl )
**      default settings
*          R(1,1) = 1.0D0 )
*          R,2,2, 1.0D0 )
*          R,3,3, 1.0D0 )
*
*        INITMATRIX( 3, 6, H )
*        INITMATRIX( 6, 6, PBest )
*
*        INITMATRIX( 6, 6, Q )
*          Q(1,1) = 0.001D0  )
*          Q,2,2, 0.001D0 )
*          Q,3,3, 0.001D0 )
*          Q,4,4, 0.00001D0  )
*          Q,5,5, 0.00001D0 )
*          Q,6,6, 0.00001D0 )
*        writeMAT( Q,5,' Q Matrix ' )
**   EnterMatrix( Q,4,'Q Matrix ' )
*        writeMAT( Q,5,' Q Matrix ' )
*
*        Say( 1,27,'  Time     Residual        Delta X     Covariance  ' )
**
*         ! ---------- Loop through new data one at a time --------------
*        First   = TRUE
*
*        DO Loops= FirstOb, LastOb
*
*            Rec= Rec
*Write(*,*)  sennum:4,jd:20:12 )
*
*            Dtsec= (JD - JDEpoch) / TUDay
*
*            CALL FindSenPtr( SenNum,SHead,STail,SenPtr )
*
*            CALL FindPhi( RNom,VNom,Dtsec,PercentChg,  Phi )
*
*               ! --------------- Prediction ----------------
*            DO i= 1, 3
*                RNom(i)=   XBest(i,1 )     ! EKF use new value of a
*                VNom(i)=   XBest(i+3,1 )
*              ENDDO
*            CALL MAG( RNom )
*            CALL MAG( VNom )
*            CALL KEPLER( RNom,VNom,Dtsec,  RNew,VNew,error )
*
*            DO i= 1, 3
*                  XPred,i,1,   RNew(i) )
*                  XPred,i+3,1, VNew(i) )
*              ENDDO
*
**         dxpred=0
*
*            CALL MATTRANS  ( Phi       ,6,6  , Phit   )
*            CALL MATMULT   ( PBest,Phit,6,6,6, PPt    )
*            CALL MATMULT   ( Phi,PPt   ,6,6,6, PPPt   )
*            CALL MATADD    ( PPPt,Q    ,6,6,   PPred  )  ! P = ÌPÌt + Q
*
*            CALL FindH( XOld,Rec,ObsRecFile,PercentChg,Dtsec,SenPtr,  LST,H )
**    xpred .or. rnom(xold)??
*
*               ! --------------- Correction ----------------
*            CALL MATTRANS  ( H        ,3,6  , Ht    )
*            CALL MATMULT   ( PPred,Ht ,6,6,3, PHt   )
*            CALL MATMULT   ( H,PHt    ,3,6,3, HPHt  )
*            CALL MATADD    ( HPHt,R   ,3,3,   HPHtR )
*            CALL MATINVERSE( HPHtR    ,3    , HPHtRI)
*            CALL MATMULT   ( PHt,HPHtRI,6,3,3, KG   )     ! K = PHt( HPHt + R )-1
*
*            CALL RV_RAZEL( RNew,VNew,RSVec,SenPtr^.SenLat,LST,'TOO',
*                      RngPred,AzPred,ElPred,DRng,DAz,DEl )
*
*            B(1,1)= Rng - RngPred
*            B(2,1)= Az  - AzPred
*            B(3,1)= El  - ElPred
*
*            CALL MATMULT   ( KG,b     ,6,3,1, Kb  )
*            CALL MATADD    ( DXPred,Kb, 6,1,DXBest)     ! DX = DX + Kb
*
*            CALL MATMULT   ( H,PPred  ,3,6,6 ,HP    )
*            CALL MATMULT   ( KG,HP    ,6,6,6, KHP   )
*            CALL MATSUB    ( PPred,KHP,  6,6, PBest )     ! P = P - KHP
*
*            CALL MATADD    ( XPred,DXBest,6,1, XBest )    ! X = X + DX
*
*            Write(*,*)  ' ---- Predictions ---- ' )
*            writeMAT( XPred,5,' X Pred - State Vector ' )
*            FILEwriteMAT( xpred,' xpred Matrix ',1.0D0,6,FileOut )
*            writeMAT( PPred,5,' P Pred - Pred COVARIANCE Matrix ' )
*            FILEwriteMAT( ppred,' ppred Matrix ',1.0D0,6,FileOut )
*
*            Write(*,*)  ' ---- Corrections ---- ' )
*            writeMAT( KG,5,   ' K      - Kalman Gain Matrix ' )
*            FILEwriteMAT( kg,' kg Matrix ',1.0D0,6,FileOut )
*
*            writeMAT( XBest,5,' X Best - State Vector ' )
*            FILEwriteMAT( xBest,' xBest Matrix ',1.0D0,6,FileOut )
*            writeMAT( dXBest,7,' DX Best ' )
*            FILEwriteMAT( dxBest,' dxBest Matrix ',1.0D0,6,FileOut )
*            writeMAT( PBest,5,' P Best - Best COVARIANCE Matrix ' )
*            FILEwriteMAT( pBest,' pBest Matrix ',1.0D0,6,FileOut )
*            Readln
*
*           ENDIF   ! DO WHILE not EOF
*
*        DO i= 1, 3
*            RNom(i)= XBest(i,1 )
*            VNom(i)= XBest(i+3,1 )
*          ENDDO
*
*        CALL MAG( RNom )
*        CALL MAG( VNom )
*      RETURN
*      END   ! Subroutine EKF
*
