       PROGRAM DCF
        IMPLICIT NONE

*  This program tests the DC capability for radar and optical data
*  There are fixes that are needed to run for other cases.
*
*                          Companion code for
*             Fundamentals of Astrodynamics and Applications
*                                  2001
*                            by David Vallado
*
*     (H)               email valladodl@worldnet.att.net
*     (W) 719-573-2600, email dvallado@stk.com
*
*     *****************************************************************
*
*  Current :
*            28 Jan 04  David Vallado
*                         Fix headers
*  Changes :
*            14 May 01  David Vallado
*                         2nd edition baseline
*            23 Nov 87  David Vallado
*                         Original baseline
*
*     *****************************************************************
*
*
*
* Uses objects dcf,astutil,astmath,asttime,astiod,ast2body,astreduc,
*              astmanv,astpert,astdc

        REAL*8 PercentChg,Epsilon,JDEpoch, RNom(3),VNom(3),
     &         NDot,NDDot
        INTEGER NumObs, FirstOb, LastOb,i
        CHARACTER TypeAns, ans

        CHARACTER*12 FileName
        INTEGER StateNum
        PARAMETER (StateNum = 8)
        REAL*8 X(StateNum,1),DX(StateNum,1),AtWAI(StateNum,StateNum),
     &         AtWA(StateNum,StateNum),Atwb(Statenum,1),
     &         PBest(StateNum,StateNum),
     &         AtWAOld(StateNum,StateNum), atwbOld(Statenum,1)

        CHARACTER*64 FileN1
        REAL*8 Mag,magrnom,magvnom
        EXTERNAL MAG



* -------------------------- Initialize Values ------------------------
        INCLUDE 'astconst.cmn'

        ! -------- Open obs file
c        Write(*,*) 'Input filename  GEOS3.dat '
c        Read(*,*) FileName
        filename = 'geos3a.dat'
        filename = 'geos6.inp'
        OPEN( 12, FILE=FileName,STATUS='UNKNOWN' )
        READ(12,*)  ! blank line if a header is included

        ! -------- Output file
        OPEN( 20, FILE='dcf.out',STATUS='UNKNOWN' )

        ! ----- Open up and form file of record of obs
        OPEN(15,FILE='ObsRec.bak',ACCESS='DIRECT', FORM='UNFORMATTED',
     &       RECL=350,STATUS='UNKNOWN')

        OPEN( 44, FILE='timerec.rec',FORM='UNFORMATTED',
     &       ACCESS='DIRECT',RECL=80,STATUS='UNKNOWN')
        FileN1 = 'nutation.dat'
        CALL InitReduc   ( FileN1 )


        ! -------- Do preliminary orbit determination -----------------
        Write(*,*) 'Input type of run - Initial OD (I) or Nominal (N) '
        Read(*,*) TypeAns
c        TypeAns = 'i'

        NDot = 0.0D0  ! these can be changed
        NDDot = 0.0D0

        IF ((TypeAns .eq. 'I').or.(TypeAns .eq. 'i')) THEN
            CALL Prelim ( TypeAns,NumObs,JDEpoch,NDot,NDDot,RNom,VNom )
          ELSE
            CALL Prelim ( TypeAns,NumObs,JDEpoch,NDot,NDDot,RNom,VNom )
*            Write(*,*) 'Input rnom eci (1-3) in m'
*            Read(*,*) rnom(1),rnom(2),rnom(3)
*            Write(*,*) 'Input vnom (1-3) in m/s'
*            Read(*,*) vnom(1),vnom(2),vnom(3)

            ! ----- Dan's test case with more obs
            rnom(1)=  -6851129.39D0          ! m
            rnom(2)=   -740891.64D0
            rnom(3)=   2913372.63D0
            vnom(1)=      2532.43400D0
            vnom(2)=     -4817.24380D0
            vnom(3)=      4755.05040D0

            rnom(1)=   6440618.63D0
            rnom(2)=  947717.14D0
            rnom(3)=  -3044197.18D0
            vnom(1)= -2998.14726D0
            vnom(2)=  4974.36356D0
            vnom(3)= -4782.23403D0
            DO i=1,3
                Rnom(i) = Rnom(i)*0.001D0
                vnom(i) = vnom(i)*0.001D0
              ENDDO

*            ! ---- simply assign the nominal value for the test runs
*            rnom(1)=  0.902412858729549D0    ! ER
*            rnom(2)=  0.419216277668340D0
*            rnom(3)=  0.539314634745115D0
*            vnom(1)=  0.545048068050250D0
*            vnom(2)= -0.243555229697641D0
*            vnom(3)= -0.725903949888685D0

            magrnom= MAG(rnom)
            magvnom= MAG(vnom)
          ENDIF
            rnom(1)=  5975.2904D0    ! km
            rnom(2)=  2568.6400D0           
            rnom(3)=  3120.5845D0
            vnom(1)=  3.983846D0
            vnom(2)= -2.071159D0
            vnom(3)= -5.917095D0

        Write(20,*)  ' Nominal '
        Write(20,*)  ' R Nominal ',RNom(1),RNom(2),RNom(3)
        Write(20,*)  ' V Nominal ',VNom(1),VNom(2),VNom(3)
 10     FORMAT( A11,3F12.6 )

        ! --- Calculate corrections based on Least Squares Analysis ---
        PercentChg   =   0.01D0
        Epsilon =   0.01D0
        FirstOb = 1
        write(*,*) numobs,'numobs'
        LastOb = NumObs
       write(*,*) 'input lastob',numobs,' are the total'
       read(*,*) lastob

        CALL LeastSquares ( PercentChg,Epsilon,JDEpoch, TypeAns,
     &                      FirstOb,LastOb,
     &                      rnom,vnom,NDot,NDDot,
     &                      X,DX,AtWAI,AtWA,atwb )

        Write(*,*) 'Beginning sequential estimate'
        Write(20,*) '---------Beginning sequential estimate'
        Write(*,*) 'Continue with sequential? y,n'
        read(*,*) ans

        if (ans.eq.'y') THEN
            Firstob= 11
            lastob = 18
            Write(*,'(A,6f12.6)') ' Nominal',RNom(1),
     &               RNom(2),RNom(3),VNom(1),VNom(2),VNom(3)

            CALL Sequential    ( PercentChg,Epsilon,JDEpoch,
     &               FirstOb,LastOb,
     &               Rnom,Vnom,
     &               AtWA,atwb,
     &               X,DX,pbest )
          endif

        CLOSE( 20 )
        CLOSE( 12 )
      STOP
      END
*
* ---------------------------------------------------------------------
*      This subroutine finds the parameters for the sensor. The variables are
*      passed back through a common block, although a record would be better.

      SUBROUTINE FindSenPtr ( SensorNum )
        IMPLICIT NONE
        INTEGER SensorNum

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

        REAL*8 DegPerSec

* ------------------------  Implementation   --------------------------
        INCLUDE 'astconst.cmn'
        INCLUDE 'astmath.cmn'

        DegPerSec     =   1.0D0/0.0710151110315201D0

        OPEN(14, FILE='sensor.txt',ACCESS='SEQUENTIAL',
     &           FORM='FORMATTED', STATUS='UNKNOWN' )

        SenNum = -1 ! must be less than the sensor number used - allow 1
        DO WHILE (SenNum .lt. SensorNum )
            READ(14,*)
     &           SenNum, SenName, SenLat, SenLon, SenAlt,
     &           RngMin, RngMax, AzMin, AzMax, ElMin, ElMax,
     &           BiasRng, BiasAz, BiasEl, BiasDRng, BiasDAz, BiasDEl,
     &           BiasTRtAsc, BiasTDecl,
     &           NoiseRng, NoiseAz, NoiseEl, NoiseDRng, NoiseDAz,
     &           NoiseDEl, NoiseTRtAsc, NoiseTDecl
          ENDDO

        SenLat = SenLat*Deg2Rad
        SenLon = SenLon*Deg2Rad
        SenAlt = SenAlt/(1000.0D0)

        NoiseRng    = NoiseRng * 0.001D0
        NoiseAz     = NoiseAz * Deg2Rad
        NoiseEl     = NoiseEl * Deg2Rad
        NoiseDRng   = NoiseDRng * 0.001D0
        NoiseDAz    = NoiseDAz * DegPerSec
        NoiseDEl    = NoiseDEl * DegPerSec
        NoiseTRtAsc = NoiseTRtAsc * Deg2Rad
        NoiseTDecl  = NoiseTDecl * Deg2Rad

        BiasRng    = BiasRng * 0.001D0
        BiasAz     = BiasAz * Deg2Rad
        BiasEl     = BiasEl * Deg2Rad
        BiasDRng   = BiasDRng * 0.001D0
        BiasDAz    = BiasDAz * DegPerSec
        BiasDEl    = BiasDEl * DegPerSec
        BiasTRtAsc = BiasTRtAsc * Deg2Rad
        BiasTDecl  = BiasTDecl * Deg2Rad

        IF (SenNum .ne. SensorNum) THEN
            Write(*,*) 'Sensor not found, using defaults '
          ENDIF

        CLOSE( 14 )

      RETURN
      END   ! FindSenPtr
*
* ---------------------------------------------------------------------
* This procedure reads data for a pass from a file, and forms the initial
*  r and v vector if its an initial OD run.
      SUBROUTINE GetRV ( TypeAns, NumObs, R,V )
        IMPLICIT NONE
        CHARACTER TypeAns
        INTEGER NumObs
        REAL*8 R(3),V(3)

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

* ----------------------------  Locals  -------------------------------
        REAL*8 VS(3)
        REAL*8 Sec,LST,GST,Tr,Ta,Te,
     &         DegPerSec, latgd,lon,alt
        INTEGER Year,Month,Day,Hr,Min,i
        INTEGER*4 SatNum

         REAL*8 MFME, DUT1, DAT, xp, yp, LOD, DDPsi, DDEps,Conv1
         REAL*8 MJD, TTT,
     &          jdut1, deltapsi,trueeps,meaneps,omega,
     &          pm(3,3),nutteme(3,3),rteme(3),vteme(3)
         REAL*8  UT1, TUT1, UTC, TAI, TDT, JDTDT, TDB, TTDB, JDTDB

         INTEGER RecNum,timezone,testnum,
     &           terms, order
         CHARACTER*12 Error

         INCLUDE 'astmath.cmn'
         INCLUDE 'astconst.cmn'
         INCLUDE 'astreduc.cmn'

* ------------------------  Implementation   --------------------------
        DegPerSec     =   1.0D0/0.0710151110315201D0

        Rng = 0.0D0
        Az  = 0.0D0
        El  = 0.0D0
        DRng= 0.0D0
        DAz = 0.0D0
        DEl = 0.0D0
        TRtAsc= 0.0D0
        TDecl = 0.0D0

        ! ------------------------
        !  Notice that this section reads a duplicate line. This is because
        !  FORTRAN doesn't seem to be able to read a portion of a line, and then
        !  the rest of the line. In PAS, or C, the option would be to read the
        !  beginning of the line, and then use a CASE statement to read
        !  the correct quantities for the remaining data.
        ! ------------------------

cdav        READ(12,*) ObsType,SatNum,SensNum
cdav fix this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        obstype = 2
        satnum=7734
        sensnum = 932

        CALL FindSenPtr( SensNum )  ! CurrSen COMMON contains all the data

        ! ---------- Add in bias values and convert units -------------
        IF (ObsType.eq.0) THEN
             READ(12,*) ObsType,SatNum,SensNum,Year,Month,Day,Hr,Min,
     &                    Sec,DRng
             DRng= DRng - BiasDRng
           ELSEIF (ObsType.eq.1) THEN
             READ(12,*) ObsType,SatNum,SensNum,Year,Month,Day,Hr,Min,
     &                    Sec,Az,El
             Az  = Az * Deg2Rad - BiasAz
             El  = El * Deg2Rad - BiasEl
           ELSEIF (ObsType.eq.2) THEN
             READ(12,*) ObsType,SatNum,SensNum,Year,Month,Day,Hr,Min,
     &                    Sec,Rng,Az,El
             Rng = Rng - BiasRng
             Az  = Az * Deg2Rad - BiasAz
             El  = El * Deg2Rad - BiasEl
           ELSEIF (ObsType.eq.3) THEN
             READ(12,*) ObsType,SatNum,SensNum,Year,Month,Day,Hr,Min,
     &                    Sec,Rng,Az,El,DRng
             Rng = Rng - BiasRng
             Az  = Az * Deg2Rad - BiasAz
             El  = El * Deg2Rad - BiasEl
             DRng= DRng - BiasDRng
           ELSEIF (ObsType.eq.4) THEN
             READ(12,*) ObsType,SatNum,SensNum,Year,Month,Day,Hr,Min,
     &                    Sec,Rng,Az,El,DRng,DAz,DEl
             Rng = Rng - BiasRng
             Az  = Az * Deg2Rad - BiasAz
             El  = El * Deg2Rad - BiasEl
             DRng= DRng - BiasDRng
             DAz = (DAz - BiasDAz)* DegPerSec
             DEl = (DEl - BiasDEl)* DegPerSec
           ELSEIF (ObsType.eq.5) THEN
             READ(12,*) ObsType,SatNum,SensNum,Year,Month,Day,Hr,Min,
     &                    Sec,TRtAsc,TDecl
             TRtAsc= TRtAsc * Deg2Rad - BiasTRtAsc
             TDecl = TDecl  * Deg2Rad - BiasTDecl
           ELSEIF (ObsType.eq.6) THEN
             READ(12,*) ObsType,SatNum,SensNum,Year,Month,Day,Hr,Min,
     &                    Sec,Rng
             Rng= Rng - BiasRng
           ENDIF  ! Case

         CALL JDay( Year,Month,Day,Hr,Min,Sec,   JD )
cdav         CALL LSTIME( SenLon,JD,        LST,Gst )
cdav         CALL SITE( SenLat,SenAlt,LST,  RSVec,VS )
cdav         CALL SITE( Latgd,Alt,Lon, RSecef,VSecef )
         CALL SITE( SenLat,SenAlt,SenLon, RSVec,VS ) !ecef

*       Temp code to "create" rtasc dec data from the az el data
*         CALL RADEC_AZEL( tRtAsc,tDecl,lst,senlat,'FROM',Az,El )
*        Write(20,*) NumObs,'rtasc ',TRtasc/deg2rad,' ',TDecl/deg2rad,
*     &               el/deg2rad


cdav         Write(*,*) NumObs,' ',lst*deg2rad,' ',JD,' ',lst*deg2rad
         Write(*,*) NumObs,' ',JD
c
cdav at this point, all the rsvec are in ecef
c
         Write(15,Rec=NumObs)
     &         SensNum, JD, (RSVec(i),i=1,3), ObsType,
     &         Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl

         ! ------ Write Header information on first pass through ------
         IF (NumObs .eq. 1) THEN
             Write( 20,*) ' LSDC RESULTS FOR Satellite ',SatNum,
     &                    ' sensor #',SenNum,SenName
             Write( 20,*)
             Write( 20,'(3(A,F11.6),A)') ' Lat ',SenLat/Deg2Rad,' Lon ',
     &                    SenLon/Deg2Rad,' Alt ',
     &                    SenAlt*1000.0,'m'
             Write( 20,'(3(A,F9.6),A)') ' Sensor Biases:  Rho ',
     &            BiasRng*1000,
     &            ' m Az ',BiasAz*rad2deg,' ø El ',BiasEl*rad2deg,' ø'
             Write( 20,'(3(A,F9.6),A)') ' Sensor Noises:  Rho ',
     &            NoiseRng*1000,' m Az ',
     &            NoiseAz*rad2deg,' ø El ',NoiseEl*rad2Deg,' ø'
           ENDIF

         ! -------------- Perform Initial Orbit Determination ---------
         IF ((TypeAns .eq. 'I').or.(TypeAns .eq. 'i')) THEN
             Tr= 0.0
             Ta= 0.0
             Te= 0.0
c             CALL RV_RAZEL( R,V,RSVec,SenLat,LST, 'FROM',
c     &                      Rng,Az,El,Tr,Ta,Te )  ! this is in eci----

         dut1 =  0.3261068D0
         dat  = 29
         xp   =  -0.11554D0* 3.1415926535879D0 / (3600.0D0*180.0D0)
         yp   =  0.48187D0* 3.1415926535879D0 / (3600.0D0*180.0D0)
         lod  =  0.0D0
         timezone=0
         order = 106
         terms = 2


         CALL CONVTIME ( Year, Month, Day, Hr, MIN, SEC,
     &                  TimeZone, 'C', DUT1, DAT, xp, yp, UT1,
     &                  TUT1, JDUT1, UTC, TAI, TDT, TTT, JDTDT,TDB,
     &                  TTDB, JDTDB, DDPsi,DDeps,Lod, Error )


             CALL RV_RAZEL( R,V,senLat,senLon,senalt,TTT,jdut1,lod,
     &                         xp,yp,terms, 'FROM',
     &                         Rng,Az,El,tr,ta,te )   ! this is in eci----
!             CALL AnglesGauss( TRtasc1, .. TDecl, (order!),
!                               JD1,JD2,JD3, RS1,RS2,RS3,r2,v2 )

             Write( 20,*)  ttt,senlon*57.29557,r(1)
     

           ENDIF

      RETURN
      END  ! Procedure GetRV
!
* ---------------------------------------------------------------------
*  This Subroutine reads in the raw data, forms position and velocity vectors,
*    and averages the middle vector calculations to form the nominal vector.
* --
      SUBROUTINE Prelim ( TypeAns,NumObs, JDEpoch,NDot,NDDot,RNom,VNom )
        IMPLICIT NONE
        CHARACTER TypeAns
        INTEGER NumObs
        REAL*8 JDEpoch, RNom(3),VNom(3), NDot, NDDot

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

* ----------------------------  Locals  -------------------------------
        REAL*8 R2V2Rec(1:60,1:2,1:3),velkmps
        REAL*8 Theta, CoplAngl, Dt,Jd1,Jd2,Jd3,p,Theta1,
     &        temp,a,ecc,inc,omega,argp,m,nu,u,l,cappi,
     &        R1(3), V1(3), R2(3), V2(3), R3(3), V3(3)
        INTEGER UseObs, i, j,NObs
        CHARACTER*12 Error
        REAL*8 Mag,magrnom,magvnom
        EXTERNAL MAG

* ------------------------  Implementation   --------------------------
        INCLUDE 'astconst.cmn'
        INCLUDE 'astmath.cmn'

        Velkmps = 7.905365719014D0
cdav        TUDay   =     0.00933809017716D0

        ! ------ Find initial orbit values based on given data --------
        CALL GetRV( TypeAns,1,R1,V1 )    ! this is in eci----
        JD1= JD

        Write( 20,* )
        Write( 20,*) ' r and v from initial data '
        Write( 20,* )
        Write( 20,* ) '       rx         ry         rz          vx   ',
     &                '      vy         vz      Coplnr     angle btwn'
        Write( 20,'(A,3(f14.7))') 'R1',R1(1),R1(2),R1(3)

        CALL GetRV( TypeAns,2,R2,V2 )
        JD2= JD

        ! ------------ Find the number of obs, then reset --------------
        NObs = 2
   98   READ( 12,*,END=99 )
cdav        READ( 12,*,END=99 ) ! do twice because the obs occur in pairs
        NObs = NObs + 1
        GOTO 98
   99   CONTINUE
        REWIND (12)
        READ( 12,* ) ! get back to current file position
        READ( 12,* ) 
        READ( 12,* ) 
        NObs = NObs - 1

        ! --------------- Read in all data and form vectors ------------
        NumObs= 2
        DO WHILE (NumObs.le.NObs)
           NumObs= NumObs + 1
           CALL GetRV( TypeAns,NumObs,R3,V3 )
           JD3= JD

           IF ((TypeAns .eq. 'I').or.(TypeAns .eq. 'i')) THEN
               CALL HerrGibbs( R1,R2,R3,jd1,jd2,jd3, V2 ,Theta,theta1,
     &                         Coplangl,Error )
               IF (Theta .gt. 1.0* Deg2Rad) THEN
                   CALL GIBBS( r1,r2,r3,  v2,theta,theta1,Coplangl,
     &                         Error )     ! this is in eci----
                 ENDIF

c               Write( 20,'(A,3(F14.7))')'R2',R2(1),R2(2),R2(3)
               Write( 20,'(A,3(F14.7))')'R2',R2(1)/rekm,
     &               R2(2)/rekm,R2(3)/rekm
               Write( 20,'(A,3(F14.7),3(F8.4),A)')' ',
     &               V2(1)/velkmps,V2(2)/velkmps,V2(3)/velkmps,
     &                 Theta/Deg2Rad,theta1/Deg2Rad,coplangl/Deg2Rad,'ø'

               IF (NumObs .lt. 60) THEN ! Limit to 60 for an initial IOD guess
                   DO i = 1,3
                       R2V2Rec(NumObs-2,1,i)= R2(i)  ! this is in eci----
                       R2V2Rec(NumObs-2,2,i)= V2(i)
                     ENDDO
                 ENDIF
               DO i = 1, 3
                   r1(i) = r2(i)
                   r2(i) = r3(i)
                 ENDDO
               JD1= JD2
               JD2= JD3
             ENDIF  ! If TypeAns = I

         ENDDO ! While

       ! -------------- Find r2 v2 vectors at epoch time -------------
       Write(20,*)
       Read(15,Rec=1)
     &         SensNum, JD, (RSVec(i),i=1,3), ObsType,
     &         Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl
       JDEpoch= JD ! get first time for entry option
       Write( 20,'(A,F18.10,A)') 'Data ranges from ',JDEpoch,' to '
       Read(15,Rec=NumObs)
     &         SensNum, JD, (RSVec(i),i=1,3), ObsType,
     &         Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl
       Temp= JD
       Write(20,'(f18.10)') Temp
       Write(20,*) ' Move r2, v2 vectors to EPOCH time '
       Write(20,*)
       Write(20,*)' Move r2, v2 vectors to EPOCH time '

       ! -------------------- Set the epoch time ---------------------
       Write(*,*) 'Input Epoch JD time 2449746.61015046D0 '
cdav       READ(*,*) JDEpoch
       JDEpoch = 2449746.61015046D0
      write(*,*) jdepoch
*      !  xxxxxxx   fix this date. The Epoch here is only for the test case

       IF ((TypeAns .eq. 'I').or.(TypeAns .eq.'i')) THEN
           DO i= 1, 3
               RNom(i)= 0.0D0
               VNom(i)= 0.0D0
             ENDDO

           IF (NumObs .lt. 60) THEN
               UseObs= NumObs - 2
             ELSE
               UseObs= 60
             ENDIF

           Write(  *,'(A,f18.8,A,I5)')'Epoch = ',JDEpoch,
     &               ' Number of obs ',UseObs
           Write( 20,'(A,f18.8,A,I5)')'Epoch = ',JDEpoch,
     &               ' Number of obs ',UseObs
           Write( 20,* )
           Write( 20,* ) '       rx         ry         rz          vx ',
     &                   '        vy         vz         dt'

           DO i= 1, UseObs
               Read(15,Rec=i+1)
     &             SensNum, JD, (RSVec(j),j=1,3), ObsType,
     &             Rng, Az, El, DRng, DAz, DEl, TRtAsc, TDecl
               Dt= (JDEpoch - JD)*86400.0D0
               DO j = 1,3
                   R2(j)= R2V2Rec(i,1,j)
                   V2(j)= R2V2Rec(i,2,j)
                 ENDDO

*               CALL KEPLER( r2,v2,dt, r3,v3,error )
               CALL PKEPLER( r2,v2,NDot,NDDot,dt, r3,v3 )  ! this is in eci----

               Write( 20,'(A,6(f14.5),f8.3,A,f18.10)')'R ',R3(1),R3(2),
     &                    R3(3),V3(1),V3(2),V3(3),dt*1440.0,' min ',JD

               ! - This is where you COULD throw out obs, ie 6,17 453 -
               DO j= 1, 3
                   RNom(j)= RNom(j) + r3(j)
                   VNom(j)= VNom(j) + v3(j)
                 ENDDO
             ENDDO

           ! -- Complete preliminary orbit determination and find the -
           ! ------------- average, or nominal state vector -----------
          DO i= 1, 3
               RNom(i)= RNom(i) / UseObs
               VNom(i)= VNom(i) / UseObs
             ENDDO
           magrnom= MAG(rnom)
           magvnom= MAG(vnom)

           Write(20,*)
           Write( *,'(A,3(f12.7))') 'Reci Noml ',RNom(1),RNom(2),RNom(3)
           Write( *,'(A,3(f12.7))') 'V Nominal ',VNom(1),VNom(2),VNom(3)
           Write( 20,* )
           Write( 20,'(A,3(f12.7))') 'Reci Nominal ',RNom(1),RNom(2),
     &                RNom(3)
           Write( 20,'(A,3(f12.7))') 'V Nominal ',VNom(1),VNom(2),
     &                VNom(3)

           CALL RV2COE( RNom,VNom, p,a,ecc,inc,omega,argp,nu,m,
     &                  u,l,cappi)
           Write( 20,* )
           Write( 20,'(A,F14.5,A)') ' a      ',a,' km'
           Write( 20,'(A,F12.8,A)') ' e      ',ecc
           Write( 20,'(A,F12.8,A)') ' i      ',inc/Deg2Rad,' ø'
           Write( 20,'(A,F12.8,A)') ' ê      ',Omega/Deg2Rad,' ø'
           Write( 20,'(A,F12.8,A)') ' w      ',Argp/Deg2Rad,' ø'
           Write( 20,'(A,F12.8,A)') ' v      ',Nu/Deg2Rad,' ø'
           Write( 20,'(A,F12.8,A)') ' M      ',M/Deg2Rad,' ø'
           Write( *,'(A,F14.5,A)') ' a      ',a,' km'
           Write( *,'(A,F12.8,A)') ' e      ',ecc
           Write( *,'(A,F12.8,A)') ' i      ',inc/Deg2Rad,' ø'
           Write( *,'(A,F12.8,A)') ' ê      ',Omega/Deg2Rad,' ø'
           Write( *,'(A,F12.8,A)') ' w      ',Argp/Deg2Rad,' ø'
           Write( *,'(A,F12.8,A)') ' v      ',Nu/Deg2Rad,' ø'
           Write( *,'(A,F12.8,A)') ' M      ',M/Deg2Rad,' ø'
         ENDIF  ! If TypeAns = 'I'

      RETURN
      END  ! Procedure Prelim



