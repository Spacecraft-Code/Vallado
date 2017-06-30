*   -------------------------------------------------------------------
*
*                              ASTMANV.FOR
*
*  this file contains fundamental astrodynamic procedures and functions
*  relating to orbit transfer calculations. ch 6 describes each of these
*  routines.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                   2007
*                             by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              30 may 07  david vallado
*                           3rd edition baseline
*    changes :
*              21 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
*
*     *****************************************************************
*
*     Uses Object files:
*         Astmath,
*         AstTime
*     Uses Common files:
*         Astmath.cmn
*
*
*      ----------- Routines for Orbit Transfer calculations ------------
*
*      SUBROUTINE Hohmann     ( RInit,RFinal,eInit,eFinal,NuInit,NuFinal,
*    &                          Deltava,Deltavb,DtTU )
*
*      SUBROUTINE BiElliptic  ( RInit,Rb,RFinal,eInit,eFinal,
*    &                          NuInit,NuFinal,
*    &                          Deltava,Deltavb,DeltaVc,DtTU )
*
*      SUBROUTINE OneTangent  ( RInit,RFinal,eInit,eFinal,NuInit,NuTran,
*    &                          Deltava,Deltavb,DtTU )
*
*      SUBROUTINE IOnlyChg    ( Deltai,VInit,fpa, DeltaViOnly )
*
*      SUBROUTINE NodeOnlyChg ( iInit,ecc,DeltaOmega,VInit,fpa,incl,
*    &                          iFinal,DeltaV )
*
*      SUBROUTINE IandNodeChg ( iInit,DeltaOmega,Deltai,VInit,fpa,
*    &                          DeltaV,iFinal )
*
*      SUBROUTINE MinCombinedPlaneChg( RInit,RFinal,eInit,eFinal,NuInit,
*    &                          NuFinal,iInit,iFinal,
*    &                          Deltai,Deltai1,DeltaVa,DeltaVb,DtTU )
*
*      SUBROUTINE CombinedPlaneChg ( RInit,RFinal,eInit,e2,eFinal,NuInit,
*    &                               Nu2a,Nu2b,NuFinal,Deltai,
*    &                               Deltava,Deltavb,DtTU )
*
*      SUBROUTINE Rendezvous  ( Rcs1,Rcs3,PhaseI,eInit,eFinal,NuInit,
*    &                          NuFinal, kTgt,kInt, PhaseF,WaitTime,
*    &                          DeltaV )
*
*      SUBROUTINE NonCoplanarRendz ( PhaseNew,Deltai,Delta2Node,LonTrue,
*    &                               RInt,RTgt,kTgt,kInt,TTrans,TPhase,
*    &                               DVPhase,DVTrans1,DVTrans2 )
*
*      SUBROUTINE HillsR      ( R,V, Alt,DtTU, RInit,VInit )
*
*      SUBROUTINE HillsV      ( R, Alt,DtTU, V )
*
*      SUBROUTINE Cow2Hill    ( rtgt,vtgt,rint,vint, RHill,VHill )
*
*      SUBROUTINE Hill2Cow    ( rTgtijk,vTgtijk,RHill,VHill,
*    &                          rIntijk,vIntijk )
*
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HOHMANN
*
*  this subroutine calculates the delta v's DO a Hohmann transfer DO either
*    circle to circle, or ellipse to ellipse.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInit       - Initial position magnitude     ER
*    RFinal      - Final position magnitude       ER
*    eInit       - Eccentricity of first orbit
*    eFinal      - Eccentricity of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad
*
*  OutPuts       :
*    DeltaVa     - Change in velocity at point a  ER / TU
*    DeltaVb     - Change in velocity at point b  ER / TU
*    DtTU        - Time of Flight DO the trans   TU
*
*  Locals        :
*    SME1        - Mech Energy of first orbit     ER2 / TU
*    SME2        - Mech Energy of transfer orbit  ER2 / TU
*    SME3        - Mech Energy of final orbit     ER2 / TU
*    VInit       - Velocity of first orbit at a   ER / TU
*    vTransa     - Velocity of trans orbit at a   ER / TU
*    vTransb     - Velocity of trans orbit at b   ER / TU
*    VFinal      - Velocity of final orbit at b   ER / TU
*    aInit       - Semimajor axis of first orbit  ER
*    aTrans      - Semimajor axis of Trans orbit  ER
*    aFinal      - Semimajor axis of final orbit  ER
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 327, Alg 36, Ex 6-1
*
* ------------------------------------------------------------------------------

      SUBROUTINE Hohmann     ( RInit,RFinal,eInit,eFinal,NuInit,
     &                          NuFinal,Deltava,Deltavb,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,eFinal,NuInit,NuFinal,
     &          Deltava,Deltavb,DtTU
* ----------------------------  Locals  -------------------------------
        REAL*8 VInit,VTrana,VTranb,VFinal, aInit,aTran,aFinal

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        aInit  = (rInit*(1.0D0+eInit*DCOS(NuInit)))
     &           / (1.0D0 - eInit*eInit )
        aTran  = ( RInit + RFinal ) / 2.0D0
        aFinal = (rFinal*(1.0D0+eFinal*DCOS(NuFinal)))
     &           / (1.0D0 - eFinal*eFinal )
        DeltaVa= 0.0D0
        DeltaVb= 0.0D0
        DtTU   = 0.0D0

        IF ( ( eInit .lt. 1.0D0 ) .or. ( eFinal .lt. 1.0D0 ) ) THEN
            ! ----------------  Find Delta v at point a  --------------
            VInit  = DSQRT( 2.0D0/rInit - (1.0D0/aInit) )
            VTrana = DSQRT( 2.0D0/rInit - (1.0D0/aTran) )
            DeltaVa= DABS( VTrana - VInit )

            ! ----------------  Find Delta v at point b  --------------
            VFinal = DSQRT( 2.0D0/rFinal - (1.0D0/aFinal) )
            VTranb = DSQRT( 2.0D0/rFinal - (1.0D0/aTran) )
            DeltaVb= DABS( VFinal - VTranb )

            ! ---------------  Find Transfer Time of Flight  ----------
            DtTU= Pi * DSQRT( aTran*aTran*aTran ) 
          ENDIF

      RETURN
      END  ! SUBROUTINE Hohmann
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE BIELLIPTIC
*
*  this subroutine calculates the delta v's DO a Bi-elliptic transfer DO either
*    circle to circle, .or. ellipse to ellipse.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInit       - Initial position magnitude     ER
*    R2          - Interim orbit magnitude        ER
*    RFinal      - Final position magnitude       ER
*    eInit       - Eccentricity of first orbit
*    eFinal      - Eccentricity of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad, Opp of NuInit
*
*  OutPuts       :
*    DeltaVa     - Change in velocity at point a  ER / TU
*    DeltaVb     - Change in velocity at point b  ER / TU
*    DtTU        - Time of Flight DO the trans   TU
*
*  Locals        :
*    SME1        - Mech Energy of first orbit     ER2 / TU
*    SME2        - Mech Energy of transfer orbit  ER2 / TU
*    SME3        - Mech Energy of final orbit     ER2 / TU
*    VInit       - Velocity of first orbit at a   ER / TU
*    vTransa     - Velocity of trans orbit at a   ER / TU
*    vTransb     - Velocity of trans orbit at b   ER / TU
*    VFinal      - Velocity of final orbit at b   ER / TU
*    aInit       - Semimajor axis of first orbit  ER
*    aTrans      - Semimajor axis of Trans orbit  ER
*    aFinal      - Semimajor axis of final orbit  ER
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 327, Alg 37, Ex 6-2
*
* ------------------------------------------------------------------------------  

      SUBROUTINE BiElliptic  ( RInit,Rb,RFinal,eInit,eFinal,
     &                         NuInit,NuFinal,
     &                         Deltava,Deltavb,DeltaVc,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,Rb,RFinal,eInit,eFinal,NuInit,NuFinal,
     &         Deltava,Deltavb,DeltaVc,DtTU
* ----------------------------  Locals  -------------------------------
        REAL*8 VInit,VTran1a,VTran1b,VTran2b,VTran2c,VFinal,
     &         aInit,aTran1,aTran2,aFinal

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        aInit  = (rInit*(1.0D0+eInit*DCOS(NuInit)))
     &           / (1.0D0 - eInit*eInit )
        aTran1 = (RInit + Rb) * 0.5D0
        aTran2 = (Rb + RFinal) * 0.5D0
        aFinal = (rFinal*(1.0D0+eFinal*DCOS(NuFinal)))
     &           / (1.0D0 - eFinal*eFinal )

        DeltaVa= 0.0D0
        DeltaVb= 0.0D0
        DeltaVc= 0.0D0
        DtTU   = 0.0D0

        IF ( ( eInit .lt. 1.0D0 ) .and. ( eFinal .lt. 1.0D0 ) ) THEN
            ! ----------------  Find Delta v at point a  --------------
            VInit  = DSQRT( 2.0D0/rInit - (1.0D0/aInit) )
            VTran1a= DSQRT( 2.0D0/rInit - (1.0D0/aTran1) )
            DeltaVa= DABS( VTran1a - VInit )

            ! ----------------  Find Delta v at point b  --------------
            VTran1b= DSQRT( 2.0D0/rb - (1.0D0/aTran1) )
            VTran2b= DSQRT( 2.0D0/rb - (1.0D0/aTran2) )
            DeltaVb= DABS( VTran1b - VTran2b )

            ! ----------------  Find Delta v at point c  --------------
            VTran2c= DSQRT( 2.0D0/rFinal - (1.0D0/aTran2) )
            VFinal = DSQRT( 2.0D0/rFinal - (1.0D0/aFinal) )
            DeltaVc= DABS( VFinal - VTran2c )

            ! ---------------  Find Transfer Time of Flight  ----------
            DtTU= Pi * DSQRT( aTran1*aTran1*aTran1 ) +
     &              Pi * DSQRT( aTran2*aTran2*aTran2 )
          ENDIF
      RETURN
      END   ! SUBROUTINE BiElliptic
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ONETANGENT
*
*  this subroutine calculates the delta v's DO a One Tangent transfer DO either
*    circle to circle, .or. ellipse to ellipse.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInit       - Initial position magnitude     ER
*    RFinal      - Final position magnitude       ER
*    eInit       - Eccentricity of first orbit
*    eFinal      - Eccentricity of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    Nu2         - True Anomaly of second orbit   Same quad as NuInit, rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad
*
*  OutPuts       :
*    DeltaVa     - Change in velocity at point a  ER / TU
*    DeltaVb     - Change in velocity at point b  ER / TU
*    DtTU        - Time of Flight DO the transf  TU
*
*  Locals        :
*    SME1        - Mech Energy of first orbit     ER2 / TU
*    SME2        - Mech Energy of transfer orbit  ER2 / TU
*    SME3        - Mech Energy of final orbit     ER2 / TU
*    VInit       - Velocity of first orbit at a   ER / TU
*    vTransa     - Velocity of trans orbit at a   ER / TU
*    vTransb     - Velocity of trans orbit at b   ER / TU
*    VFinal      - Velocity of final orbit at b   ER / TU
*    aInit       - Semimajor axis of first orbit  ER
*    aTrans      - Semimajor axis of Trans orbit  ER
*    aFinal      - Semimajor axis of final orbit  ER
*    E           - Ecc anomaly of trans at b      rad
*    Ratio       - Ratio of initial to final
*                    orbit radii
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 335, Alg 38, Ex 6-3
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE OneTangent  ( RInit,RFinal,eInit,eFinal,NuInit,NuTran,
     &                         Deltava,Deltavb,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,eFinal,NuInit,NuTran,
     &         Deltava,Deltavb,DtTU
* ----------------------------  Locals  -------------------------------
        REAL*8 EAInit,VInit,VTrana,VTranb,VFinal, eTran,aInit,aTran,
     &         aFinal, fpaTranb,fpaFinal, E, Sinv,Cosv,Ratio

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        DeltaVa= 0.0D0 
        DeltaVb= 0.0D0
        DtTU   = 0.0D0
        Ratio  = rInit/rFinal 
        IF ( DABS(NuInit) .lt. 0.01D0 ) THEN ! check 0 .or. 180
            eTran  = ( Ratio-1.0D0 ) / ( DCOS(NuTran)-Ratio )   ! init at perigee
            EAInit= 0.0D0 
          ELSE
            eTran  = ( Ratio-1.0D0 ) / ( DCOS(NuTran)+Ratio )  ! init at apogee
            EAInit= Pi 
          ENDIF 
        IF ( eTran .ge. 0.0D0 ) THEN
            aInit = (rInit*(1.0D0+eInit*DCOS(NuInit)))
     &              / (1.0D0 - eInit*eInit )
            aFinal= (rFinal*(1.0D0+eFinal*DCOS(NuTran)))
     &              / (1.0D0 - eFinal*eFinal )
*                       nutran is used since it = nufinal!!  
*aInit = rinit
*afinal= rfinal
            IF ( DABS( eTran-1.0D0 ) .gt. 0.000001D0 ) THEN
                IF ( DABS(NuInit) .lt. 0.01D0 ) THEN ! check 0 .or. 180
                    aTran = (rInit*(1.0D0+eTran*DCOS(NuInit)))
     &                      / (1.0D0 - eTran*eTran ) ! per
                  ELSE
*                 aTran = (rInit*(1.0D0+eTran*DCOS(NuInit))) / (1.0D0 + eTran*eTran )  apo  
                    aTran= RInit/(1.0D0 + eTran)
                  ENDIF
              ELSE
                aTran = 999999.9D0   ! Infinite DO Parabolic orbit
              ENDIF 

            ! ----------------  Find Delta v at point a  --------------
            VInit  = DSQRT( 2.0D0/rInit - (1.0D0/aInit) )
            VTrana = DSQRT( 2.0D0/rInit - (1.0D0/aTran) )
            DeltaVa= DABS( VTrana - VInit )

            ! ----------------  Find Delta v at point b  --------------
            VFinal  = DSQRT( 2.0D0/rFinal - (1.0D0/aFinal) )
            VTranb  = DSQRT( 2.0D0/rFinal - (1.0D0/aTran) )
            fpaTranb= DATAN( ( eTran*DSIN(NuTran) )
     &                / ( 1.0D0 + eTran*DCOS(NuTran) ) )
            fpaFinal= DATAN( ( eFinal*DSIN(NuTran) )
     &                / ( 1.0D0 + eFinal*DCOS(NuTran) ) )
            DeltaVb = DSQRT( VTranb*VTranb + VFinal*VFinal
     &                       - 2.0D0*VTranb*VFinal
     &                       *DCOS( fpaTranb-fpaFinal ) )

            ! ---------------  Find Transfer Time of Flight  ----------
            IF ( eTran .lt. 0.99999D0 ) THEN
                Sinv= ( DSQRT( 1.0D0-eTran*eTran )*DSIN(NuTran) )
     &                 / ( 1.0D0 + eTran*DCOS(NuTran) )
                Cosv= (eTran+DCOS(NuTran))/(1.0D0+eTran*DCOS(NuTran) )
                E   = DATAN2( Sinv,Cosv ) 
                DtTU= DSQRT( aTran*aTran*aTran ) *
     &                    ( E - eTran*DSIN(E)
     &                    - (EAInit - ETran*DSIN(EAInit)) )
              ELSE
                IF ( DABS( eTran-1.0D0 ) .lt. 0.000001D0 ) THEN
*                  Parabolic DtTU  
                  ELSE
*                  Hyperbolic DtTU  
                  ENDIF
              ENDIF

          ELSE
            Write(*,*) 'one tangent burn is not possible DO this case '
          ENDIF
      RETURN
      END   ! SUBROUTINE OneTangent
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE IONLYCHG
*
*  this subroutine calculates the delta v's DO a change in inclination only.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    DeltaI      - Change in inclination          rad
*    VInit       - Initial velocity vector        ER/TU
*    fpa         - Flight path angle              rad
*
*  OutPuts       :
*    DeltaVionly - answer
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 346, Alg 39, Ex 6-4
*
* ------------------------------------------------------------------------------  

      SUBROUTINE IOnlyChg    ( Deltai,VInit,fpa, DeltaViOnly )
        IMPLICIT NONE
        REAL*8 Deltai,VInit,fpa, DeltaViOnly

        ! --------------------  Implementation   ----------------------
        DeltaViOnly = 2.0D0*VInit*DCOS(fpa)*DSIN(0.5D0*Deltai) 
      RETURN
      END   ! SUBROUTINE IOnlyChg
    
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NODEONLYCHG
*
*  this subroutine calculates the delta v's for a change in longitude of
*    ascending node only.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    DeltaOmega  - Change in Node                 Rad
*    ecc         - Ecc of first orbit
*    VInit       - Initial velocity vector        ER/TU
*    fpa         - Flight path angle              rad
*    Incl        - Inclination                    rad
*
*
*  OutPuts       :
*    iFinal      - Final inclination              rad
*    Deltav      - Change in velocity             ER/TU
*
*  Locals        :
*    VFinal      - Final velocity vector          ER/TU
*    ArgLat      - Argument of latitude           rad
*    ArgLat1     - Final Argument of latitude     rad
*    NuInit      - Initial true anomaly           rad
*    Theta       -
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 349, Alg 40, Ex 6-5
*
* ------------------------------------------------------------------------------  

      SUBROUTINE NodeOnlyChg ( iInit,ecc,DeltaOmega,VInit,fpa,incl,
     &                         iFinal,DeltaV )
        IMPLICIT NONE
        REAL*8 iInit,ecc,DeltaOmega,VInit,fpa,incl,iFinal,DeltaV
* ----------------------------  Locals  -------------------------------
        REAL*8 ArgLat,ArgLat1,Theta

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        IF ( DABS(ecc) .gt. 0.00000001D0 ) THEN
            ! ------------------------ Elliptical ---------------------
            Theta = DATAN( DSIN(iInit)*DTAN(DeltaOmega) ) 
            iFinal= DASIN( DSIN(Theta)/DSIN(DeltaOmega) ) 
            DeltaV= 2.0D0*VInit*DCOS(fpa)*DSIN(0.5D0*Theta)

            ArgLat = Pi*0.5D0  ! set at 90 deg
            ArgLat1= DACOS( DCOS(incl)*DSIN(incl)*
     &               (1.0D0-DCOS(DeltaOmega))/ DSIN(Theta) )
          ELSE
            ! ------------------------- Circular ----------------------
            theta = DACOS( DCOS(iinit)**2
     &                +DSIN(iinit)**2*DCOS(DeltaOmega) )
            DeltaV= 2.0D0*VInit*DSIN(0.5D0*Theta)

            ArgLat = DACOS( DTAN(iInit)*(DCOS(DeltaOmega)-DCOS(Theta))
     &                      / DSIN(Theta) )
            ArgLat1= DACOS( DCOS(incl)*DSIN(incl)*
     &               (1.0D0-DCOS(DeltaOmega))/ DSIN(Theta) )

          ENDIF

      RETURN
      END   ! SUBROUTINE NodeOnlyChg
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE IandNodeChg
*
*  this subroutine calculates the delta v's for a change in inclination .and.
*    longitude of ascending node.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    VInit       - Initial velocity vector        ER/TU
*    iInit       - Initial inclination            rad
*    fpa         - Flight path angle              rad
*    DeltaOmega  - Change in Node                 Rad
*    DeltaI      - Change in inclination          Rad
*    RFinal      - Final position magnitude       ER
*
*  OutPuts       :
*    iFinal      - Final inclination              rad
*    Deltav      - Change in velocity             ER/TU
*
*  Locals        :
*    ArgLat      - Argument of latitude           rad
*    ArgLat1     - Final Argument of latitude     rad
*    Theta       -
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 350, Alg 41, Ex 6-6
*
* ------------------------------------------------------------------------------

      SUBROUTINE IandNodeChg ( iInit,DeltaOmega,Deltai,VInit,fpa,
     &                         DeltaV,iFinal )
        IMPLICIT NONE
        REAL*8 iInit,DeltaOmega,Deltai,VInit,fpa, DeltaV,iFinal
* ----------------------------  Locals  -------------------------------
        REAL*8 ArgLat, ArgLat1, Theta

        ! --------------------  Implementation   ----------------------
        iFinal= iInit - Deltai
        theta = DACOS( DCOS(iinit)*DCOS(ifinal) +
     &                    DSIN(iinit)*DSIN(ifinal)*DCOS(DeltaOmega) )
        DeltaV= 2.0D0*VInit*DCOS(fpa)*DSIN(0.5D0*Theta)

        ArgLat = DACOS( (DSIN(ifinal)*DCOS(DeltaOmega) -
     &                  DCOS(Theta)*DSIN(iinit))
     &                  / (DSIN(Theta)*DCOS(iinit)) )
        ArgLat1= DACOS( (DCOS(iInit)*DSIN(iFinal) -
     &                  DSIN(IInit)*DCOS(iFinal)*DCOS(DeltaOmega))
     &                  / DSIN(Theta) )

      RETURN
      END   ! SUBROUTINE IandNodeChg
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MINCOMBINEDPLANECHG
*
*  this subroutine calculates the delta v's .and. the change in inclination
*    necessary DO the minimum change in velocity when traveling between two
*    non-coplanar orbits.  The notation used is from the initial orbit (1) at
*    point a, transfer is made to the transfer orbit (2), .and. to the final
*    orbit (3) at point b.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    RInit       - Initial position magnitude     ER
*    RFinal      - Final position magnitude       ER
*    eInit       - Ecc of first orbit
*    e2          - Ecc of trans orbit
*    eFinal      - Ecc of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad
*    iInit       - Incl of the first orbit        rad
*    iFinal      - Incl of the second orbit       rad
*
*  OutPuts       :
*    Deltai1     - Amount of incl chg req at a    rad
*    DeltaVa     - Change in velocity at point a  ER / TU
*    DeltaVb     - Change in velocity at point b  ER / TU
*    DtTU        - Time of Flight DO the trans   TU
*    NumIter     - Number of iterations
*
*  Locals        :
*    SME1        - Mech Energy of first orbit     ER2 / TU
*    SME2        - Mech Energy of transfer orbit  ER2 / TU
*    SME3        - Mech Energy of final orbit     ER2 / TU
*    VInit       - Velocity of first orbit at a   ER / TU
*    vTransa     - Velocity of trans orbit at a   ER / TU
*    vTransb     - Velocity of trans orbit at b   ER / TU
*    VFinal      - Velocity of final orbit at b   ER / TU
*    aInit       - Semimajor axis of first orbit  ER
*    aTrans      - Semimajor axis of Trans orbit  ER
*    aFinal      - Semimajor axis of final orbit  ER
*    e2          - Eccentricity of second orbit
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 355, Alg 42, Table 6-3
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE MinCombinedPlaneChg( RInit,RFinal,eInit,eFinal,NuInit,
     &                         NuFinal,iInit,iFinal,
     &                         Deltai,Deltai1,DeltaVa,DeltaVb,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,eFinal,NuInit,NuFinal,iInit,iFinal,
     &         Deltai,Deltai1,DeltaVa,DeltaVb,DtTU
* ----------------------------  Locals  -------------------------------
        INTEGER numiter
        REAL*8 deltainew,temp,TDi, SME1,SME2,SME3, VInit,
     &         V1t,V3t,VFinal, a1,a2,a3

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        ! -------------------  Initialize values   --------------------
        a1  = (RInit*(1.0D0+eInit*DCOS(NuInit))) / (1.0D0 - eInit**2 )
        a2  = 0.5D0 * (RInit+RFinal) 
        a3  = (RFinal*(1.0D0+eFinal*DCOS(NuFinal))) / (1.0D0-eFinal**2 )
        SME1= -1.0D0 / (2.0D0*a1) 
        SME2= -1.0D0 / (2.0D0*a2) 
        SME3= -1.0D0 / (2.0D0*a3) 

        ! ----------------------- Find velocities ---------------------
        VInit = DSQRT( 2.0D0*( (1.0D0/RInit) + SME1 ) ) 
        V1t   = DSQRT( 2.0D0*( (1.0D0/RInit) + SME2 ) ) 

        VFinal= DSQRT( 2.0D0*( (1.0D0/RFinal) + SME3 ) ) 
        V3t   = DSQRT( 2.0D0*( (1.0D0/RFinal) + SME2 ) ) 

        ! ---------- Find the optimum change of inclination -----------
        TDi = iFinal-iInit 

        Temp= (1.0D0/TDi) * DATAN( (RFinal/RInit**1.5D0 - DCOS(TDi))
     &        / DSIN(TDi) )
        Temp= (1.0D0/TDi) * DATAN( DSIN(TDi)
     &        / (RFinal/RInit**1.5D0 + DCOS(TDi)) )

        DeltaVa= DSQRT( V1t*V1t + VInit*VInit
     &                 - 2.0D0*V1t*VInit*DCOS(Temp*Tdi) )
        DeltaVb= DSQRT( V3t*V3t + VFinal*VFinal
     &                 - 2.0D0*V3t*VFinal*DCOS(TDi*(1.0D0-Temp)) )

        Deltai = Temp*TDi 
        Deltai1= TDi*(1.0D0-Temp) 

        ! ---------------  Find Transfer Time of Flight  --------------
        DtTU= Pi * DSQRT( A2*A2*A2 ) 

        ! ---- Iterate to find the optimum change of inclination ------
        DeltaiNew  = Deltai    ! 1st guess, 0.01D0 to 0.025D0 seems good
        Deltai1    = 100.0D0   ! IF ( going to smaller orbit, should be
        NumIter    = 0         ! 1.0D0 - 0.025D0!

        DO WHILE (DABS(DeltaiNew-Deltai1) .gt. 0.000001D0)
            Deltai1= DeltaiNew
            DeltaVa= DSQRT( V1t*V1t + VInit*VInit
     &               - 2.0D0*V1t*VInit* DCOS(Deltai1) )

            DeltaVb= DSQRT( V3t*V3t + VFinal*VFinal
     &               - 2.0D0*V3t*VFinal* DCOS(TDi-Deltai1) )

            DeltaiNew= DASIN( (DeltaVa*VFinal*V3t*DSIN(TDi-Deltai1))
     &                 / (VInit*V1t*DeltaVb) )
            NumIter = NumIter + 1
          ENDDO   ! DO WHILE DABS()

      RETURN
      END   ! SUBROUTINE MinCombinedPlaneChg
*
*    Vallado       2007, 355, Ex 6-7

      SUBROUTINE CombinedPlaneChg ( RInit,RFinal,eInit,e2,eFinal,NuInit,
     &                              Nu2a,Nu2b,NuFinal,Deltai,
     &                              Deltava,Deltavb,DtTU )
        IMPLICIT NONE
        REAL*8 RInit,RFinal,eInit,e2,eFinal,NuInit,Nu2a,Nu2b,
     &         NuFinal,Deltai, Deltava,Deltavb,DtTU
* ----------------------------  Locals  -------------------------------
        REAL*8 SME1,SME2,SME3, VInit,vTransa,vTransb,VFinal, a1,a2,a3,
     &         fpa1,fpa2a,fpa2b,fpa3,E,Eo,Sinv,Cosv

        ! --------------------  Implementation   ----------------------
        a1 = (RInit*(1.0D0+eInit*DCOS(NuInit))) / (1.0D0-eInit*eInit )
        IF ( DABS( e2-1.0D0 ) .gt. 0.000001D0 ) THEN
            a2  = (RInit*(1.0D0+e2*DCOS(Nu2a))) / (1.0D0 - e2*e2 )
            SME2= -1.0D0 / (2.0D0*a2)
          ELSE
            a2  = 999999.9D0   ! Undefined DO Parabolic orbit
            SME2= 0.0D0 
          ENDIF 
        a3 = (RFinal*(1.0D0+eFinal*DCOS(NuFinal)))
     &        / (1.0D0 - eFinal*eFinal )
        SME1= -1.0D0 / (2.0D0*a1) 
        SME3= -1.0D0 / (2.0D0*a3)

        ! ----------------  Find Delta v at point a  ------------------
        VInit = DSQRT( 2.0D0*( (1.0D0/RInit) + SME1 ) ) 
        vTransa= DSQRT( 2.0D0*( (1.0D0/RInit) + SME2 ) )
        fpa2a= DATAN( ( e2*DSIN(Nu2a) ) / ( 1.0D0 + e2*DCOS(Nu2a) ) ) 
        fpa1 = DATAN( ( eInit*DSIN(NuInit) )
     &                 / ( 1.0D0 + eInit*DCOS(NuInit) ) )
        DeltaVa= DSQRT( vTransa*vTransa + VInit*VInit
     &           - 2.0D0*vTransa*VInit*( DSIN(fpa2a)*DSIN(fpa1)
     &           + DCOS(fpa2a)*DCOS(fpa1)*DCOS(Deltai)) )

        ! ----------------  Find Delta v at point b  ------------------
        VFinal = DSQRT( 2.0D0*( (1.0D0/RFinal) + SME3 ) ) 
        vTransb= DSQRT( 2.0D0*( (1.0D0/RFinal) + SME2 ) ) 
        fpa2b= DATAN( ( e2*DSIN(Nu2b) ) / ( 1.0D0 + e2*DCOS(Nu2b) ) ) 
        fpa3 = DATAN( ( eFinal*DSIN(NuFinal) )
     &                / ( 1.0D0 + eFinal*DCOS(NuFinal) ) )
        DeltaVb= DSQRT( vTransb*vTransb + VFinal*VFinal
     &           - 2.0D0*vTransb*VFinal*( DSIN(fpa2b)*DSIN(fpa3)
     &           + DCOS(fpa2b)*DCOS(fpa3)*DCOS(Deltai)) )

        ! ---------------  Find Transfer Time of Flight  --------------
        Sinv= ( DSQRT( 1.0D0-e2*e2 )*DSIN(Nu2b) )
     &        / ( 1.0D0 + e2*DCOS(Nu2b) )
        Cosv= ( e2+DCOS(Nu2b) ) / ( 1.0D0 + e2*DCOS(Nu2b) )
        E= DATAN2( Sinv,Cosv ) 
        Sinv= ( DSQRT( 1.0D0-e2*e2 )*DSIN(Nu2a) )
     &        / ( 1.0D0 + e2*DCOS(Nu2a) )
        Cosv= ( e2+DCOS(Nu2a) ) / ( 1.0D0 + e2*DCOS(Nu2a) ) 
        Eo= DATAN2( Sinv,Cosv ) 
        DtTU= DSQRT( A2**3 ) * ( (E - e2*DSIN(E)) - (Eo-e2*DSIN(Eo)) )

      RETURN
      END   ! SUBROUTINE CombinedPlaneChg
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE RENDEZVOUS
*
*  this subroutine calculates parameters for a Hohmann transfer rendezvous.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rcs1        - Radius of circular orbit int   ER
*    Rcs2        - Radius of circular orbit tgt   ER
*    eInit       - Ecc of first orbit
*    eFinal      - Ecc of final orbit
*    NuInit      - True Anomaly of first orbit    0 .or. Pi rad
*    NuFinal     - True Anomaly of final orbit    0 .or. Pi rad
*    PhaseI      - Initial phase angle (Tgt-Int)  +(ahead) .or. -(behind) rad
*    NumRevs     - Number of revs to wait
*    kTgt        -
*    kInt        -
*
*  OutPuts       :
*    PhaseF      - Final Phase Angle              rad
*    WaitTime    - Wait before next intercept opp TU
*    Deltav      - Change in velocity             ER/TU
*
*  Locals        :
*    DtTUTrans   - Time of flight of trans orbit  TU
*    ATrans      - Semimajor axis of trans orbit  ER
*    AngVelTgt   - Angular velocity of target     rad / TU
*    AngVelInt   - Angular velocity of int        rad / TU
*    LeadAng     - Lead Angle                     rad
*
*  Coupling      :
*    None
*
*  References    :
*    Vallado       2007, 364, Alg 44, Alg 45, Ex 6-8, Ex 6-9
*
* ------------------------------------------------------------------------------  
*
      SUBROUTINE Rendezvous  ( Rcs1,Rcs3,PhaseI,eInit,eFinal,NuInit,
     &                         NuFinal, kTgt,kInt, PhaseF,WaitTime,
     &                         DeltaV )
        IMPLICIT NONE
        REAL*8 Rcs1,Rcs3,PhaseI,eInit,eFinal,NuInit,NuFinal,
     &         PhaseF,WaitTime,DeltaV
        INTEGER kTgt,kInt
* ----------------------------  Locals  -------------------------------
        REAL*8 PeriodTrans,Rp,DtTUTrans,LeadAng,aTrans,
     &         AngVelTgt,AngVelInt, a1,a2,a3,VInit,vTransa,VFinal,
     &         vTransb,SME1,SME2,SME3,DeltaVa,DeltaVb,VInt,VTrans

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        ATrans    = (Rcs1 + Rcs3) / 2.0D0 
        DtTUTrans = Pi*DSQRT( ATrans*ATrans*ATrans ) 
        AngVelInt = 1.0D0 / ( DSQRT(Rcs1*Rcs1*Rcs1) ) 
        AngVelTgt = 1.0D0 / ( DSQRT(Rcs3*Rcs3*Rcs3) ) 
        VInt      = DSQRT( 1.0D0/Rcs1 ) 

        ! --------- Check DO satellites in the same orbits ------------
        IF ( DABS( AngVelInt - AngVelTgt ) .lt. 0.000001D0 ) THEN
            PeriodTrans= ( kTgt*TwoPi + PhaseI ) / AngVelTgt
            aTrans     = (PeriodTrans/(TwoPi*kInt))**(2.0D0/3.0D0)
            Rp         = 2.0D0*aTrans - Rcs1 
            IF ( Rp .lt. 1.0D0 ) THEN
                Write(*,*) 'Error - transfer orbit intersects Earth'
              ENDIF
            VTrans  = DSQRT( (2.0D0/Rcs1) - (1.0D0/aTrans) ) 
            DeltaV  = 2.0D0*(VTrans-VInt) 
            WaitTime= 0.0D0 
      leadang= 0.0D0
          ELSE
            LeadAng = AngVelTgt * DtTUTRans
            PhaseF  = LeadAng - Pi 
            WaitTime= ( PhaseF - PhaseI + 2.0D0*Pi*kTgt )
     &            / ( AngVelInt - AngVelTgt )

            a1  = (rcs1*(1.0D0+eInit*DCOS(NuInit)))
     &            / (1.0D0 - eInit*eInit )
            a2  = ( Rcs1 + Rcs3 ) / 2.0D0 
            a3  = (rcs3*(1.0D0+eFinal*DCOS(NuFinal)))
     &            / (1.0D0 - eFinal*eFinal )
            SME1= -1.0D0 / (2.0D0*a1) 
            SME2= -1.0D0 / (2.0D0*a2) 
            SME3= -1.0D0 / (2.0D0*a3) 
            ! ----------------  Find Delta v at point a  --------------
            VInit  = DSQRT( 2.0D0*( (1.0D0/Rcs1) + SME1 ) )
            vTransa= DSQRT( 2.0D0*( (1.0D0/Rcs1) + SME2 ) ) 
            DeltaVa= DABS( vTransa - VInit ) 

            ! ----------------  Find Delta v at point b  --------------
            VFinal = DSQRT( 2.0D0*( (1.0D0/Rcs3) + SME3 ) ) 
            vTransb= DSQRT( 2.0D0*( (1.0D0/Rcs3) + SME2 ) ) 
            DeltaVb= DABS( VFinal - vTransb )
            DeltaV = DeltaVa + DeltaVb
          ENDIF 

      RETURN
      END   ! SUBROUTINE Rendezvous
*
* -----
*    Vallado       2007, 370, Alg 46, Ex 6-10
*
* ------
      SUBROUTINE NonCoplanarRendz ( PhaseNew,Deltai,Delta2Node,LonTrue,
     &                              RInt,RTgt,kTgt,kInt,TTrans,TPhase,
     &                              DVPhase,DVTrans1,DVTrans2 )
        IMPLICIT NONE
        REAL*8 PhaseNew,Deltai,Delta2Node,LonTrue,RInt,RTgt,
     &         TTrans,TPhase,DVPhase, DVTrans1, DVTrans2
        INTEGER kTgt, kInt
* ----------------------------  Locals  -------------------------------
        REAL*8 AngVelInt,AngVelTgt,atrans,aphase,lead,
     &         leadNew,TNode,LonTrueNew,VInt,VTgt,VPhase,VTrans1,
     &         VTrans2

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        AngVelInt= DSQRT( 1.0D0/(RInt*RInt*RInt) ) 
        AngVelTgt= DSQRT( 1.0D0/(RTgt*RTgt*RTgt) ) 
        ATrans   = (RInt + RTgt) * 0.5D0
        TTrans = Pi*DSQRT( ATrans*ATrans*ATrans ) 

        Lead = AngVelTgt * TTRans 

        TNode= Delta2Node/AngVelInt 

        LonTrueNew= LonTrue + AngVelTgt*TNode 
*fix     PhaseNew= 13.5D0/Rad  
        LeadNew= Pi + PhaseNew

        TPhase= (LeadNew - Lead + TwoPi*kTgt) / AngVelTgt 

        aPhase = (TPhase/(TwoPi*kInt))**(2.0D0/3.0D0)

        ! ----------------  Find Deltav's  -----------------
        VInt= DSQRT(1.0D0/RInt) 
        VPhase= DSQRT(2.0D0/RInt - 1.0D0/aPhase)
        DVPhase= VPhase - VInt 

        VTrans1= DSQRT(2.0D0/RInt - 1.0D0/aTrans) 
        DVTrans1= VTrans1 - VPhase 

        VTrans2= DSQRT(2.0D0/RTgt - 1.0D0/aTrans) 
        VTgt= DSQRT(1.0D0/RTgt) 
        DVTrans2= DSQRT(VTgt*VTgt + VTrans2*VTrans2
     &            - 2.0D0*VTgt*VTrans2*DCOS(Deltai))

      RETURN
      END   ! SUBROUTINE NonCoplanarRendz
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HILLSR
*
*  this subroutine calculates various position information for Hills equations.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Initial Position vector of INT ER
*    V           - Initial Velocity Vector of INT ER / TU
*    Alt         - Altitude of TGT satellite      ER
*    DtTU        - Desired Time                   TU
*
*  Outputs       :
*    RInit       - Final Position vector of INT   ER
*    VInit       - Final Velocity Vector of INT   ER / TU
*
*  Locals        :
*    nt          - Angular velocity times time    rad
*    Omega       -
*    Sinnt       - Sine of nt
*    Cosnt       - Cosine of nt
*    Radius      - Magnitude of range vector      ER
*
*  Coupling      :
*
*
*  References    :
*    Vallado       2007, 397, Alg 47, Ex 6-14
*
* ------------------------------------------------------------------------------

      SUBROUTINE HillsR      ( R,V, Alt,DtTU, RInit,VInit )
        IMPLICIT NONE
        REAL*8 R(3), V(3),  Alt,DtTU, RInit(3), VInit(3)
* ----------------------------  Locals  -------------------------------
        REAL*8 SinNt,CosNt,Omega,nt,Radius

        ! --------------------  Implementation   ----------------------
        Radius= 1.0D0 + Alt
        Omega = DSQRT( 1.0D0 / (Radius*Radius*Radius) ) 
        nt    = Omega*DtTU 
        CosNt = DCOS( nt ) 
        SinNt = DSIN( nt ) 

        ! --------------- Determine new positions  --------------------
        RInit(1)= ( V(1)/Omega ) * SinNt -
     &           ( (2.0D0*V(2)/Omega) + 3.0D0*R(1) ) * CosNt +
     &           ( (2.0D0*V(2)/Omega) + 4.0D0*R(1) )
        RInit(2)= ( 2.0D0*V(1)/Omega ) * CosNt +
     &           ( (4.0D0*V(2)/Omega) + 6.0D0*R(1) ) * SinNt +
     &           ( R(2) - (2.0D0*V(1)/Omega) ) -
     &           ( 3.0D0*V(2) + 6.0D0*Omega*R(1) )*DtTU
        RInit(3)= R(3)*CosNt + (V(3)/Omega)*SinNt

        ! --------------- Determine new velocities  -------------------
        VInit(1)= V(1)*CosNt + (2.0D0*V(2)+3.0D0*Omega*R(1))*SinNt
        VInit(2)= -2.0D0*V(1)*SinNt + (4.0D0*V(2)
     &         +6.0D0*Omega*R(1))*CosNt - (3.0D0*V(2)+6.0D0*Omega*R(1))
        VInit(3)= -R(3)*Omega*SinNt + V(3)*CosNt

      RETURN
      END   ! SUBROUTINE HillsR
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HILLSV
*
*  this subroutine calculates initial velocity DO Hills equations.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Initial Position vector of INT ER
*    Alt         - Altitude of TGT satellite      ER
*    DtTU        - Desired Time                   TU
*
*  Outputs       :
*    V           - Initial Velocity Vector of INT ER / TU
*
*  Locals        :
*    Numer       -
*    Denom       -
*    nt          - Angular velocity times time    rad
*    Omega       -
*    Sinnt       - Sine of nt
*    Cosnt       - Cosine of nt
*    Radius      - Magnitude of range vector      ER
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 410, Eq 6-60, Ex 6-15
*
* ------------------------------------------------------------------------------  

      SUBROUTINE HillsV      ( R, Alt,DtTU, V )
        IMPLICIT NONE
        REAL*8 R(3), Alt,DtTU, V(3)
* ----------------------------  Locals  -------------------------------
        REAL*8 Numer,Denom,SinNt,CosNt,Omega,nt,Radius, DCot
        EXTERNAL DCOT

        ! --------------------  Implementation   ----------------------
        Radius= 1.0D0 + Alt
        Omega = DSQRT( 1.0D0 / (Radius*Radius*Radius) ) 
        nt    = Omega*DtTU 
        CosNt = DCOS( nt ) 
        SinNt = DSIN( nt )

        ! --------------- Determine initial Velocity ------------------
        Numer= ( (6.0D0*r(1)*(nt-SinNt)-r(2))*Omega*SinNt
     &          - 2.0D0*Omega*r(1)*(4.0D0-3.0D0*CosNt)*(1.0D0-CosNt) )
        Denom= (4.0D0*SinNt-3.0D0*nt)*SinNt + 4.0D0*( 1.0D0-CosNt )
     &          *( 1.0D0-CosNt )

        IF ( DABS( Denom ) .gt. 0.000001D0 ) THEN
            V(2)= Numer / Denom
          ELSE
            V(2)= 0.0D0
          ENDIF
        IF ( DABS( SinNt ) .gt. 0.000001D0 ) THEN
            V(1)= -( Omega*r(1)*(4.0D0-3.0D0*CosNt)
     &            +2.0D0*(1.0D0-CosNt)*v(2) ) / ( SinNt )
          ELSE
            V(1)= 0.0D0
          ENDIF
        V(3)= -R(3)*Omega*DCot(nt)

      RETURN
      END   ! SUBROUTINE HillsV
*
      SUBROUTINE IJK_RSW     ( Rijk,Vijk,R,S,W, Direction,Rrsw,Vrsw )
        IMPLICIT NONE
        REAL*8 Rijk(3),Vijk(3),R(3),S(3),W(3),Rrsw(3),Vrsw(3)
        Character*4 Direction

        REAL*8 ROTMat(3,3)
        INTEGER i

        IF ( Direction .eq. 'FROM' ) THEN
            ! ---------------- Form Rotation Matrix -------------------
            DO i= 1, 3
                ROTMat(i,1)= R(i)
                ROTMat(i,2)= S(i)
                ROTMat(i,3)= W(i)
              ENDDO

            ! ----------------- Do multiplication ---------------------
            DO i= 1, 4
                rijk(i)= 0.0D0
                vijk(i)= 0.0D0
              ENDDO
            CALL MatVecMult( ROTMat,Rrsw,3,3,1, 3,3,1, Rijk )
            CALL MatVecMult( ROTMat,Vrsw,3,3,1, 3,3,1, Vijk )
          ELSE
            ! --------------- Form Rotation Matrix --------------------
            DO i= 1, 3
                ROTMat(1,i)= R(i)
                ROTMat(2,i)= S(i)
                ROTMat(3,i)= W(i)
              ENDDO

            ! ----------------- Do multiplication ---------------------
            DO i= 1, 4
                Rrsw(i)= 0.0D0
                Vrsw(i)= 0.0D0
              ENDDO
            CALL MatVecMult( ROTMat,Rijk,3,3,1, 3,3,1, Rrsw )
            CALL MatVecMult( ROTMat,Vijk,3,3,1, 3,3,1, Vrsw )
          ENDIF 

      RETURN
      END   ! SUBROUTINE IJKRSW
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE COW2Hill
*
*  this subroutine finds the equivalent relative motion vector given a geocentric
*    vector.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rtgt        - Position vector of tgt         ER
*    Vtgt        - Velocity Vector of tgt         ER / TU
*    Rint        - Position vector of int         ER
*    Vint        - Velocity Vector of int         ER / TU
*
*  Outputs       :
*    RHill       - Position vector of int rel to
*                  target                         ER
*    VHill       - Velocity Vector of int rel to
*                  target                         ER / TU
*
*  Locals        :
*    None.
*
*  Coupling      :
*    CROSS       - Cross product of two vectors
*    NORM        - Unit vector
*    IJK_RSW     -
*    LNCOM2      - Linear combination of two scalars .and. two vectors
*    ROT3        - Rotation about the 3rd axis
*
*  References    :
*    Vallado       2007, 413
*
* ------------------------------------------------------------------------------

      SUBROUTINE Cow2Hill    ( rtgt,vtgt,rint,vint, RHill,VHill )
        IMPLICIT NONE
        REAL*8 rtgt(3), vtgt(3), rint(3), vint(3), RHill(3), VHill,
     &         MAG

        EXTERNAL MAG

* ----------------------------  Locals  -------------------------------
        REAL*8 r(3), s(3), w(3), rv(3),  rt(3), magrt, magri,
     &         vt(3), ri(3), vi(3), angly, anglz

        ! --------------------  Implementation   ----------------------
        ! ---------- Form RSW unit vectors DO transformation ----------
        CALL NORM ( rtgt,  R )
        CALL CROSS( rtgt,vtgt, RV )
        CALL NORM ( RV,  W )
        CALL CROSS( W,R, S )

        CALL IJK_RSW( RTgt,VTgt,R,S,W,'TOO', Rt,Vt )   ! IJK to RSW
        magrt = MAG(rt)
        CALL IJK_RSW( RInt,VInt,R,S,W,'TOO', Ri,Vi )   ! IJK to RSW

       ! --- Determine z offset to correct vector ----
        IF ( DABS(ri(3)) .gt. 0.0000001D0 ) THEN
            anglz= DATAN(ri(3)/magrt)  ! coord sys based at tgt
            CALL ROT2( Ri, -anglz, Ri )
            CALL ROT2( Vi, -anglz, Vi )  ! should be ROT2(a), but opp
          ELSE
            anglz= 0.0D0 
          ENDIF

        ! -------------- Determine y offset to correct vector ---------
        IF ( DABS(ri(2)) .gt. 0.0000001D0 ) THEN
            angly= DATAN(ri(2)/magrt)  ! should be ROT3(-a), but opp, but sign DO later
            CALL ROT3( Ri, angly, Ri )
            CALL ROT3( Vi, angly, Vi )
          ELSE
            angly= 0.0D0
          ENDIF

        ! ---------------------------- Do all 3 here ------------------
*     LNCOM2( 1.0D0,-1.0D0,ri,rt, rHill )
        CALL LNCOM2( 1.0D0,-1.0D0,vi,vt, VHill )

        ! ------------------- Now add in corrections ------------------
        magri = MAG(ri)
        RHill(1)= ri(1) - magrt  !4if not do rotri or 1
        RHill(2)= angly*magri
        RHill(3)= anglz*magri

      RETURN
      END   ! SUBROUTINE Cow2Hill
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE HILL2COW
*
*  this subroutine finds the equivalent geocentric vector given the target
*    and relative motion vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Rtgt        - Position vector of tgt         ER
*    Vtgt        - Velocity Vector of tgt         ER / TU
*    RHill       - Position vector of int rel to
*                  target                         ER
*    VHill       - Velocity Vector of int rel to
*                  target                         ER / TU
*
*  Outputs       :
*    Rint        - Position vector of int         ER
*    Vint        - Velocity Vector of int         ER / TU
*
*  Locals        :
*    None.
*
*  Coupling      :
*    IJKRSW      - Form translation matrix given position and velocity vectors
*    ROT2        - Rotation about the 2nd axis
*    ROT3        - Rotation about the 3rd axis
*    MAG         - Magnitude of a vector
*    NORM        - Unit vector
*    CROSS       - Cross product of two vectors
*
*  References    :
*    Vallado       2007, 414
*
* ------------------------------------------------------------------------------

      SUBROUTINE Hill2Cow    ( rTgtijk,vTgtijk,RHill,VHill,
     &                         rIntijk,vIntijk )
        IMPLICIT NONE
        REAL*8 rTgtijk(3), vTgtijk(3), RHill(3), VHill(3),
     &         rIntijk(3), vIntijk(3), MAG
        EXTERNAL MAG

* ----------------------------  Locals  -------------------------------
        REAL*8 angly, anglz, rtem(3), vtem(3),  RTgtrsw(3), VTgtrsw(3),
     &         rv(3), R(3), S(3), W(3), magrtem

        ! --------------------  Implementation   ----------------------
        ! --- Form RSW unit vectors DO transformation ----
        CALL NORM ( rTgtijk,  R )
        CALL CROSS( rTgtijk,vTgtijk, RV )
        CALL NORM ( RV,  W )
        CALL CROSS( W,R, S )

        ! IJK to RSW
        CALL IJK_RSW( RTgtijk,VTgtijk,R,S,W,'TOO', RTgtrsw,VTgtrsw )

        RTem(1)= RTgtrsw(1)+RHill(1)   ! in RSW
        RTem(2)= RTgtrsw(2)
        RTem(3)= RTgtrsw(3)
        vTem(1)= vTgtrsw(1)+VHill(1)
        vTem(2)= vTgtrsw(2)+VHill(2)
        vTem(3)= vTgtrsw(3)+VHill(3)
        magrtem = MAG( rTem)

        ! --- Now perform rotation to fix y ----
        IF ( DABS(RHill(2)) .gt. 0.0000001D0 ) THEN
            angly= DATAN(RHill(2)/magrTem)  ! rtgt, but IF ( x non-zero, needs extra
            CALL ROT3( rTem,-angly,  rTem )  ! should be ROT3(a) but opp
            CALL ROT3( vTem,-angly,  vTem )
          ELSE
            angly= 0.0D0
          ENDIF

        ! --- Now perform rotation to fix z ----
        IF ( DABS(RHill(3)) .gt. 0.0000001D0 ) THEN
            anglz= DATAN(RHill(3)/magrTem)  ! should be ROT2(-a), but opp
            CALL ROT2( rTem,anglz,  rTem )
            CALL ROT2( vTem,anglz,  vTem )
          ELSE
            anglz= 0.0D0
          ENDIF

        ! --- Now RSW to IJK via MatMult!! ----
        CALL IJK_RSW( RIntijk,VIntijk,R,S,W, 'FROM', RTem,VTem )
      RETURN
      END   ! SUBROUTINE Hill2Cow
*
