*   -------------------------------------------------------------------
*
*                              ASTPERT.FOR
*
*   this file contains fundamental astrodynamic procedures and functions
*   allowing the calcualtions of perturbations. most of these routines
*   are discussed in ch 8.
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
*     Uses object files:
*          Astmath
*          AstTime
*          Ast2body
*          AstIod
*
*      SUBROUTINE PKEPLER     ( Ro,Vo,nDot,Nddot,DtSec,     R,V          )
*
*      SUBROUTINE J2DragPert  ( Incl,Ecc,N,NDot,  OmegaDOT,ArgpDOT,EDOT )
*
*      SUBROUTINE PREDICT     ( JD,latgd,LST, r,v,rs, WhichKind,
*     &                           Rho,Az,El,tRtasc,tDecl, Vis )
*
*      SUBROUTINE Deriv       ( X,  XDot )
*
*      SUBROUTINE InitGravityField   ( Order, C,S )
*
*      SUBROUTINE LegPoly     ( Latgc, Order, LArr )
*
*      SUBROUTINE FullGeop    ( R,V, ITime,WhichOne,BC,Order,C,S,APert  )
*
*      SUBROUTINE PERTACCEL   ( R,V, ITime, WhichOne, BC, APert )
*
*      SUBROUTINE PDERIV      ( ITime,X,DerivType,BC,  XDot )
*
*      SUBROUTINE RK4         ( ITime,DtDay,XDot,DerivType,BC,  X )
*
*      SUBROUTINE RKF45       ( ITime, DtDay, XDot,DerivType, BC, X )
*
*      SUBROUTINE Cowell      ( R,V,ITime,FTime,DtDay,DerivType,BC, R1,V1 )
*
*      SUBROUTINE ATMOS       ( R, Rho )
*
*
* ------------------- Constants used in this Library ------------------
*
*     J2         : REAL*8 =    0.00108263D0
*     J3         : REAL*8 =   -0.00000254D0
*     J4         : REAL*8 =   -0.00000161D0
*     GMS        : REAL*8 =    3.329529364D0E05
*     GMM        : REAL*8 =    0.01229997D0
*
* ---------------------------------------------------------------------

*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PKEPLER
*
*  this subroutine propagates a satellite's position and velocity vector over
*    a given time period accounting for perturbations caused by J2.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ro          - original position vector       km
*    Vo          - original velocity vector       km/s
*    NDot        - Time rate of change of n       rad/s
*    NDDot       - Time accel of change of n      rad/s2
*    DtSec        - Change in time                 s
*
*  Outputs       :
*    R           - updated position vector        km
*    V           - updated velocity vector        km/s
*
*  Locals        :
*    P           - Semi-paramter                  km
*    A           - semior axis                    km
*    Ecc         - eccentricity
*    incl        - inclination                    rad
*    Argp        - argument of periapsis          rad
*    ArgpDot     - change in argument of perigee  rad/s
*    Omega       - longitude of the asc node      rad
*    OmegaDot    - change in Omega                rad
*    E0          - eccentric anomaly              rad
*    E1          - eccentric anomaly              rad
*    M           - mean anomaly                   rad/s
*    MDot        - change in mean anomaly         rad/s
*    ArgLat      - argument of latitude           rad
*    ArgLatDot   - change in argument of latitude rad/s
*    TrueLon     - true longitude of vehicle      rad
*    TrueLonDot  - change in the true longitude   rad/s
*    LonPerg     - longitude of periapsis         rad
*    LonPeroDot  - longitude of periapsis change  rad/s
*    N           - mean angular motion            rad/s
*    NUo         - true anomaly                   rad
*    J2oP2       - J2 over p sqyared
*    Sinv,Cosv   - Sine .and. Cosine of Nu
*
*  Coupling:
*    rv2coe       - Orbit Elements from position .and. Velocity vectors
*    coe2rv       - Position .and. Velocity Vectors from orbit elements
*    NEWTONM     - Newton Rhapson to find Nu .and. Eccentric anomaly
*
*  References    :
*    Vallado       2007, 687, Alg 64
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE PKEPLER       ( Ro,Vo,nDot,Nddot,DtSec,     R,V     )
        IMPLICIT NONE
        REAL*8 Ro(3), Vo(3), nDot,nDDot,DtSec, R(3), V(3)
* -----------------------------  Locals  ------------------------------
        REAL*8 Tndto3,p, a, Ecc, Incl, Omega, Argp, Nu, M, ArgLat,
     &         TrueLon, LonPerg, OmegaDot, E0, ArgpDot, MDot,ArgLatDot,
     &         TrueLonDot, LonPerDot, n, J2oP2, J2

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        J2    =  0.00108263D0

        CALL rv2coe( Ro,Vo,  p,a,Ecc,incl,Omega,Argp,Nu,M,ArgLat,
     &              TrueLon,LonPerg )
        n= DSQRT(mu/(A*A*A))

        ! ------------ Find the value of J2 perturbations -------------
        J2oP2   = (n*rekm*rekm*1.5D0*J2) / (p*p)
*     NBar    = n*( 1.0D0 + J2oP2*DSQRT(1.0D0-Ecc*Ecc)* (1.0D0 - 1.5D0*DSIN(Incl)*DSIN(Incl)) )
        OmegaDot= -J2oP2 * DCOS(Incl)
        ArgpDot =  J2oP2 * (2.0D0-2.5D0*DSIN(Incl)*DSIN(Incl))
        MDot    =  N

        Tndto3= 2.0D0*NDot*DtSec / (3.0D0*n)
        a     = a - Tndto3 * a
*     edot  = -Tndto3 * (1.0D0-Ecc)/DtSec
        Ecc   = Ecc - Tndto3 * (1.0D0-Ecc)
        p     = a*(1.0D0 - Ecc*Ecc)

        ! ---- Update the orbital elements DO each orbit type ---------
        IF ( Ecc .lt. Small ) THEN
           ! --------------  Circular Equatorial  ---------------------
           IF ( (incl .lt. Small) .or. (DABS(Incl-Pi).lt.Small) ) THEN
               TrueLonDot= OmegaDot + ArgpDot + MDot
               TrueLon   = TrueLon  + TrueLonDot * DtSec
               TrueLon   = DMOD(TrueLon, TwoPi)
             ELSE
               ! ---------------  Circular Inclined    ----------------
               Omega    = Omega + OmegaDot * DtSec
               Omega    = DMOD(Omega, TwoPi)
               ArgLatDot= ArgpDot + MDot
               ArgLat   = ArgLat + ArgLatDot * DtSec
               ArgLat   = DMOD(ArgLat, TwoPi)
             ENDIF
           ELSE
             ! ----- Elliptical, Parabolic, Hyperbolic Equatorial -----
             IF ( ( incl .lt. Small ) .or.
     &            ( DABS(Incl-Pi) .lt. Small ) ) THEN
                  LonPerDot= OmegaDot + ArgpDot
                  LonPerg  = LonPerg + LonPerDot * DtSec
                  LonPerg  = DMOD(LonPerg, TwoPi)
                  M        = M + MDOT*DtSec +NDot*DtSec**2 +
     &                       NDdot*DtSec**3
                  M        = DMOD(M, TwoPi)
                  CALL NEWTONM( Ecc,M,  e0,Nu )
                ELSE
                  ! ---- Elliptical, Parabolic, Hyperbolic Inclined ---
                  Omega= Omega + OmegaDot * DtSec
                  Omega= DMOD(Omega, TwoPi)
                  Argp = Argp  + ArgpDot  * DtSec
                  Argp = DMOD(Argp, TwoPi)
                  M    = M + MDOT*DtSec + NDot*DtSec**2 + NDdot*DtSec**3
                  M    = DMOD(M, TwoPi)
                  CALL NEWTONM( Ecc,M,  e0,Nu )
                ENDIF
           ENDIF

         ! ------------ Use coe2rv to find new vectors -----------------
         CALL coe2rv(P,Ecc,Incl,Omega,Argp,Nu,ArgLat,TrueLon,LonPerg,
     &               R,V)
      RETURN
      END  ! SUBROUTINE PKEPLER
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE J2DRAGPERT
*
*  this subroutine calculates the perturbations for the PREDICT problem
*    involving secular rates of change resulting from J2 .and. Drag only.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    incl        - Inclination                    rad
*    Ecc         - Eccentricity
*    N           - Mean Motion                    rad/s
*    NDot        - Mean Motion rate               rad / 2TU2
*
*  Outputs       :
*    OmegaDot    - Long of Asc Node rate          rad / s
*    ArgpDot     - Argument of perigee rate       rad / s
*    EDot        - Eccentricity rate              / s
*
*  Locals        :
*    P           - Semiparameter                  km
*    A           - Semimajor axis                 km
*    NBar        - Mean Mean motion               rad / s
*
*  Coupling      :
*    None
*
*  References    :
*    Vallado       2007, 645, Eq 9-37, 646, Eq 9-39, 648, Eq 9-50
*
* ------------------------------------------------------------------------------

      SUBROUTINE J2DragPert ( Incl,Ecc,N,NDot,  OmegaDOT,ArgpDOT,EDOT )
        IMPLICIT NONE
        REAL*8 Incl,Ecc,N,NDot, OmegaDOT,ArgpDOT,EDOT

* -----------------------------  Locals  ------------------------------
        REAL*8 P,A,J2,NBar

        ! --------------------  Implementation   ----------------------
        J2  =  0.00108263D0

        a   = (1.0D0/n) ** (2.0D0/3.0D0)
        p   = a*(1.0D0 - ecc**2)
        NBar= n*( 1.0D0+1.5D0*J2*(DSQRT(1.0D0-Ecc*Ecc)/(p*p))*
     &          ( 1.0D0-1.5D0*DSIN(Incl)**2 ))

* ------------------------- Find dot Terms  ---------------------------
        OmegaDot = -1.5D0*( J2/(p*p) ) * DCOS(Incl) * NBar
        ArgpDot  =  1.5D0*( J2/(p*p) ) * (2.0D0-2.5D0*DSIN(Incl)**2) *
     &                    Nbar
        EDot     = -(4.0D0/3.0D0) * (1.0D0-Ecc) * (NDot/Nbar)

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PREDICT
*
*  this subroutine determines the azimuth .and. elevation DO the viewing
*    of a satellite from a known ground SITE.  Notice the Julian Date is left
*    in it's usual DAYS format because the DOT terms are input as radians per
*    day, thus no extra need DO conversion. Setup with vectors to simplify the
*    use with any propagator.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date of desired obs     Day
*    Latgd       - Geodetic Latitude of SITE      -Pi to Pi rad
*    LST         - Local SIDEREAL Time            -2Pi to 2Pi rad
*    R           - updated position vector        km
*    V           - updated velocity vector        km/s
*    RS          - IJK SITE Vector                km
*    WhichKind   - Type of Sunrise                'S''C''N''A'
*
*  OutPuts       :
*    Rho         - Range from SITE to satellite   km
*    Az          - Azimuth                        rad
*    El          - Elevation                      rad
*    TRtAsc      - Topo Right ascension           rad
*    TDecl       - Topo Declination               rad
*    Vis         - Visibility
*                  'Radar SUN'   - both in SUN
*                  'Eye'  - SITE dark, sat in SUN
*                  'Radar Nite'  - both dark
*                  'Not Visible' - sat below horz
*  Locals        :
*    Temp        - Temporary Real value
*    SRtAsc      - Suns Right ascension           rad
*    SDecl       - Suns Declination               rad
*    SatAngle    - ANGLE between IJK SUN .and. Sat  rad
*    Dist        - Ppdculr dist of sat from RSun  km
*    rr          - Range rate
*    Drr         - Range acceleration
*    Dtrtasc     - Topocentric rtasc rate
*    DRho        - Slant range rate
*    DAz         - Azimuth rate
*    Del         - Elevation rate
*    SunAngle    - ANGLE between SUN .and. SITE     rad
*    AngleLimit  - ANGLE DO twilight conditions  rad
*    RhoVec      - SITE to sat vector in SEZ      km
*    TempVec     - Temporary vector
*    RHoV        - SITE to sat vector in IJK      km
*    RSun        - SUN vector                     AU
*    C           - Temporary Vector
*
*  Coupling      :
*    SUN         - Position vector of SUN
*    CROSS       - CROSS Product of two vectors
*    ROT2,ROT3   - Rotations about 2nd .and. 3rd axis
*    LNCOM1      - Combination of a vector .and. a scalar
*    LNCOM2      - Combination of two vectors .and. two scalars
*    RV_RAZEL    - Conversion with vectors .and. range azimuth elevation
*    RV_TRADEC   - Conversion with topocentric right ascension declination
*    ANGLE       - ANGLE between two vectors
*
*  References    :
*    Vallado       2007, 900, Alg 73, Ex 11-6
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE PREDICT       ( JD,latgd,LST, r,v,rs, WhichKind,
     &                           Rho,Az,El,tRtasc,tDecl, Vis )
        IMPLICIT NONE
        REAL*8 JD, Latgd, LST, r(3),v(3),RS(3), Rho, Az, El, trtasc,
     &         tdecl
        CHARACTER WhichKind
        CHARACTER*11 Vis
        EXTERNAL MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 RhoVec(3), TempVec(3), RhoV(3), RSun(3), C(3), MAG,
     &         rr, drr, dtrtasc, dtdecl, SRtAsc, SDecl, Dist,
     &         drho, daz, del, SunAngle, SatAngle, AngleLimit,
     &         magc, magrsun,magr

        INCLUDE 'astmath.cmn'
        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        Az    =  0.0D0
        El    =  0.0D0
        Rho   =  0.0D0
        TRtAsc=  0.0D0
        TDecl =  0.0D0

        ! ------ Find IJK range vector from SITE to satellite ---------
        CALL LNCOM2( 1.0D0,-1.0D0,R,RS,  RhoV )
        Rho= MAG(RhoV)

        ! ------- Calculate Topocentric Rt Asc .and. Declination ------
        CALL RV_TRADEC(r,v,rs,'TOO',rr,trtasc,tdecl,Drr,Dtrtasc,Dtdecl)

        ! ---------------------- Rotate to SEZ ------------------------
        CALL ROT3( RhoV,       LST   ,  TempVec )
        CALL ROT2( TempVec,HalfPi-Latgd,   RhoVec )

        ! --------------- Check visibility constraints ----------------
        ! ------------------ Is it above the Horizon ------------------
        IF ( RhoVec(3) .gt. 0.0D0 ) THEN
            ! --------- Is the SITE in the LIGHT, or the dark? --------
            CALL SUN( JD,RSun,SRtAsc,SDecl )
            CALL LNCOM1( AUER,RSun, RSun )
            CALL ANGLE( RSun,RS, SunAngle )
            IF ( WHICHKind .eq.'S') THEN
               AngleLimit= (90.0D0 + 50.0D0/60.0D0)*Deg2Rad
              ELSE
                IF ( WHICHKind .eq.'C') THEN
                    AngleLimit=  96.0D0*Deg2Rad
                  ELSE
                    IF ( WHICHKind .eq.'N') THEN
                        AngleLimit= 102.0D0*Deg2Rad
                      ELSE
                        IF ( WHICHKind .eq.'A') THEN
                            AngleLimit= 108.0D0*Deg2Rad
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF

            IF ( SunAngle .lt. AngleLimit ) THEN
                Vis= 'Day        '
              ELSE
                ! ---------- This assumes a conical shadow ------------
                ! ----- Is the satellite in the shadow .or. not? ------
                CALL CROSS( RSun, R, C )
                Magc = MAG(c)
                Magr = MAG(r)
                Magrsun = MAG(rsun)
                SatAngle= DASIN( magc/ (magrsun*magr) )
                Dist= magr*DCOS( SatAngle - HalfPi )
                IF ( Dist .gt. 1.0D0 ) THEN
                    Vis= 'Terminator '
                  ELSE
                    Vis= 'Night      '
                  ENDIF
              ENDIF
          ELSE
            Vis= 'not visible'
          ENDIF      

        ! -----------  Calculate Azimuth .and. Elevation  -------------
c        CALL RV_RAZEL( Reci,Veci,Latgd,Lon,alt,TTT,jdut1,lod,
c     &                         xp,yp,terms, 'TOO',
c     &                         Rho,Az,El,DRho,DAz,DEl )
c need to define transformation coordinates.


c old way        CALL RV_RAZEL( r,v,rs,latgd,LST,'TOO', rho,az,el,drho,daz,del )
      RETURN
      END ! Predict
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE DERIV
*
*  this subroutine calculates the derivative of the two-body state vector
*    use with the Runge-Kutta algorithm.  Note time is not needed.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    X           - State Vector                   km, km/s
*
*  Outputs       :
*    XDot        - Derivative of State Vector     km/s,  km/TU2
*
*  Locals        :
*    RCubed      - Cube of R
*
*  Coupling      :
*    None.
*
*  References    :
*    None.
*
* ------------------------------------------------------------------------------

      SUBROUTINE Deriv ( X,  XDot )
        IMPLICIT NONE
        REAL*8 X(6), XDot(6)

* -----------------------------  Locals  ------------------------------
        Real*8 RCubed

        ! --------------------  Implementation   ----------------------
        RCubed= ( DSQRT( X(1)**2 + X(2)**2 + X(3)**2 ) )**3

        write(*,*) rcubed,'   rcubed',x(4)
        write(*,*) x(5),'  ',x(6)
        ! -----------------  Velocity Terms  --------------------------
        XDot(1)= X(4)
        XDot(2)= X(5)
        XDot(3)= X(6)

        ! --------------  Acceleration Terms   ------------------------
        XDot(4)= -X(1) / RCubed
        XDot(5)= -X(2) / RCubed
        XDot(6)= -X(3) / RCubed

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE InitGravityField
*
*  this subroutine reads .and. stores the gravity field DO use in the program.
*    coefficients. The routine can be configured DO either normalized .or.
*    unnormalized values.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Order       - Size of gravity field          1..70D0
*
*  Outputs       :
*    C           - Gravitational Coefficients
*    S           - Gravitational Coefficients
*
*  Locals        :
*    RCubed      - Cube of R
*
*  Coupling      :
*    None.
*
*  References    :
*    None.
*
* ------------------------------------------------------------------------------

      SUBROUTINE InitGravityField   ( Order, C,S )
        IMPLICIT NONE
        INTEGER Order
        REAL*8 C(70,70), S(70,70)

        INTEGER l, m, cexp, sexp
        REAL*8 cnor, snor, Cval, Sval
        CHARACTER*8 NOTEOF

        ! --------------------  Implementation   ----------------------
        OPEN( UNIT=25,File='c:\tplib\wgs84.dat',STATUS='OLD' )

* --------------- Set up Loop to READ through Input File --------------
        NOTEOF = 'TRUE'
        DO WHILE (NOTEOF.eq.'TRUE')
            Read(25,*,END=999) l,m,Cnor,Snor
            C(l,m)= Cnor  ! unnormalized values
            S(l,m)= Snor
          ENDDO  ! While Not EOF ( END=999 in Read )

 999    CLOSE( 25 )

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE LegPoly
*
*  this subroutine finds the Legendre polynomials DO the gravity field.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Latgc       - Geocentric Latitude of SITE    -Pi to Pi rad
*    Order       - Size of gravity field          1..70D0
*
*  Outputs       :
*    LArr        - Array of Legendre Polynomials
*
*  Locals        :
*    L,m         - Indices of gravitational potential
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 593, Eq 8-56
*
* ------------------------------------------------------------------------------

      SUBROUTINE LegPoly     ( Latgc, Order, LArr )
        IMPLICIT NONE
        REAL*8 Latgc, LArr(0:70,0:70)
        INTEGER Order

        INTEGER L,m

        ! --------------------  Implementation   ----------------------
        LArr(0,1)= 0.0D0
        LArr(0,0)= 1.0D0
        LArr(1,0)= DSIN(Latgc)
        LArr(1,1)= DCOS(Latgc)

        ! ------------------- Perform Recursions ----------------------
        DO L= 2,Order
            LArr(0,L-1)= 0.0D0
            DO m= 0,L
                IF ( m .eq. 0 ) THEN
                      LArr(L,0)= ( (2*L-1)* LArr(1,0) * LArr(L-1,0)
     &                           - (L-1)* LArr(L-2,0) )/L
                  ELSE
                    IF ( m .eq. L ) THEN
                        LArr(L,m)= (2*L-1) * LArr(1,1) * LArr(L-1,m-1)
                      ELSE
                        LArr(L,m)= LArr(L-2,m)
     &                           + (2*L-1) * LArr(1,1) * LArr(L-1,m-1)
                      ENDIF
                  ENDIF  
              ENDDO   ! DO m
          ENDDO   ! DO L

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE FullGeop
*
*  this subroutine finds the Legendre polynomial value DO the gravity field
*    DO a given order.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Order       - Size of gravity field          1..70D0
*
*  Outputs       :
*    C           - Gravitational Coefficients
*    S           - Gravitational Coefficients
*
*  Locals        :
*    RCubed      - Cube of R
*
*  Coupling      :
*    IJKTOLATLONA- Find sub satellite point
*
*  References    :
*    Vallado       2007,
*
* ------------------------------------------------------------------------------

      SUBROUTINE FullGeop    ( R,V, ITime,WhichOne,BC,Order,C,S,APert )
        IMPLICIT NONE
        REAL*8 R(3), V(3), ITime,BC,C(70,70),S(70,70),APert(6)
        INTEGER Order, WhichOne
        EXTERNAL MAG
        INTEGER L, m
        REAL*8 LArr(0:70,0:70), OORDelta, Temp, OOr,  SumM1, SumM2,MAG,
     &         SumM3, DistPartr, DistPartPhi, DistPartLon, RDelta,
     &         Latgc, latgd, hellp, Lon, LastOOr, SumL1,SumL2,magr
        ! --------------------  Implementation   ----------------------
        CALL IJK2llA( R,ITime, latgc,latgd,lon,hellp )

        ! -------------------- Find Legendre polynomials --------------
        CALL LegPoly( Latgc,Order, LArr )

        ! --------- Partial derivatives of disturbing potential -------
        magr = MAG(r)
        OOr= 1.0D0/magr
        LastOOr= 1.0D0/magr
        SumM1= 0.0D0
        SumM2= 0.0D0
        SumM3= 0.0D0
        DO L= 2,Order
            DO m= 0,L
                SumM1= SumM1 + LArr(L,m) * C(L,m)*DCOS(m*Lon)
     &                       + S(L,m)*DSIN(m*Lon)
                SumM2= SumM2 + C(L,m)*DCOS(m*Lon)
     &                       + S(L,m)*DSIN(m*Lon) *
     &                       ( LArr(L,m+1) - LArr(L,m)*m*DTAN(Latgc) )
                SumM3= SumM3 + m*LArr(L,m) * (S(L,m)*DCOS(m*Lon)
     &                       - C(L,m)*DSIN(m*Lon))
                SumL1 = 0.0D0 ! fixnjnjnjjjjnjnjnjn
                SumL2 = 0.0D0 ! fix
              ENDDO
          ENDDO

        DistPartR  = -OOr*OOr*SumL1 * SumM1
        DistPartPhi=  OOr*SumL2     * SumM2
        DistPartLon=  OOr*SumL2     * SumM3

        ! --------- Non-spherical pertubative acceleration ------------
        RDelta  = DSQRT( r(1)*r(1) + r(2)*r(2) )
        OORdelta= 1.0D0/RDelta
        Temp    = OOr*DistPartR - r(3)*OOr*OOr*OORDelta*DistPartPhi

        APert(1)= Temp*r(1) - OORDelta*DistPartLon*r(2)
        APert(2)= Temp*r(2) + OORDelta*DistPartLon*r(1)
        APert(3)= OOr*DistPartR*r(3) + OOR*OOr*RDelta*DistPartPhi

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PERTACCEL
*
*  this subroutine calculates the actual value of the perturbing acceleration.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Radius vector                  km
*    V           - Velocity vector                km/s
*    Time        - Initial time (Julian Date)     Days from 4713 BC
*    WhichOne    - Which perturbation to calc     1 2 3 4 5 ...
*    BC          - Ballistic Coefficient          kg/m2
*
*  Outputs       :
*    APert       - Perturbing acceleration        km/TU2
*
*  Locals        :
*    rs2         - SUN radius vector **2
*    rs3         - SUN radius vector **3
*    rm2         - MOON radius vector **2
*    rm3         - MOON radius vector **3
*    r32         - "z" component of Radius vec **2
*    r33         - "z" component of Radius vec **3
*    r34         - "z" component of Radius vec **4
*    r2          - Radius vector **2
*    r3          - Radius vector **3
*    r4          - Radius vector **4
*    r5          - Radius vector **5
*    r7          - Radius vector **7
*    Beta        -
*    Temp        - Temporary Real Value
*    rho         - Atmospheric Density
*    Va          - Relative Velocity Vector       km / s
*    RSun        - Radius Vector to SUN           AU
*    RMoon       - Radius Vector to MOON          km
*    RtAsc       - Right Ascension                rad
*    Decl        - Declination                    rad
*    AUER        - Conversion from AU to km
*    Temp1       -
*    Temp2       -
*    i           - Index
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    SUN         - SUN vector
*    MOON        - MOON vector
*    ATMOS       - Atmospheric density
*
*  References    :
*
* ------------------------------------------------------------------------------

      SUBROUTINE PERTACCEL   ( R,V, ITime, WhichOne, BC, APert )
        IMPLICIT NONE
        REAL*8 R(3), V(3), ITime,BC,APert(6)
        INTEGER WhichOne
        EXTERNAL MAG
* -----------------------------  Locals  ------------------------------
        INTEGER i
* fix the c and s vars
        REAL*8 C(70), S(70)

        REAL*8 Va(3), RSun(3), RMoon(3), rs2, rm2, rs3, rm3, r32, r33,
     &         r34, r2, r3, r4, r5, r7, Beta, Temp, rho, srtasc, magva,
     &         sdecl, mrtasc, mdecl, Temp1, Temp2, J2, J3, magapert,
     &         J4, GMS, GMM, DOT, Mag, magr, magv, magrsun, magrmoon
        EXTERNAL DOT

        INCLUDE 'astconst.cmn'

        ! --------------------  Implementation   ----------------------
        J2         =    0.00108263D0
        J3         =   -0.00000254D0
        J4         =   -0.00000161D0
        GMS        =    3.329529364D5
        GMM        =    0.01229997D0
        magr = MAG( R )
        magv = MAG( V )
        R2 = magr*magr
        R3 = R2*magr
        R4 = R2*R2
        R5 = R2*R3
        R7 = R5*R2
        R32= r(3)*r(3)
        R33= R32*r(3)
        R34= R32*R32

        ! -----------------   J2 Acceleration   -----------------------
        IF ( WhichOne .eq. 1 ) THEN
                Temp1=  (-1.5D0*J2) / R5
                Temp2=  1.0D0 - (5.0D0*R32) / R2
                APert(1)= Temp1*r(1) * Temp2  ! recheck with formulae
                APert(2)= Temp1*r(2) * Temp2
                APert(3)= Temp1*r(3) * ( 3.0D0-(5.0D0*R32) / R2 )
              ENDIF

        ! ------------------   J3 Acceleration   ----------------------
        IF ( WhichOne .eq. 2 ) THEN
                Temp1=  (-2.5D0*J3) / R7
                Temp2=  3.0D0*r(3)-(7.0D0*R33) / R2
                APert(1)= Temp1*r(1) * Temp2
                APert(2)= Temp1*r(2) * Temp2
                IF ( DABS( r(3) ) .gt. 0.0000001D0 ) THEN
                    APert(3)= Temp1*r(3) * ((6.0D0*r(3))-((7.0D0*R33)
     &                           / R2) - ( (3.0D0*r2) / (5.0D0*r(3)) ))
                  ELSE
                    APert(3)= 0.0D0
                  ENDIF
              ENDIF

        ! -----------------    J4 Acceleration   ----------------------
        IF ( WhichOne .eq. 3 ) THEN
                Temp1=  (-1.875D0*J4) / R7
                Temp2=  1.0D0-((14.0D0*R32)/R2)+((21.0D0*R34) / R4)
                APert(1)= Temp1*r(1) * Temp2
                APert(2)= Temp1*r(2) * Temp2
                APert(3)= Temp1*r(3) * (5.0D0-((70.0D0*R32)/(3.0D0*R2))
     &                    +((21.0D0*R34) / R4 ))
              ENDIF

        ! -----------------   SUN Acceleration   ----------------------
        IF ( WhichOne .eq. 4 ) THEN
                CALL SUN( ITime,RSun,SRtAsc,SDecl )
                DO i= 1,3
                    RSun(i)= RSun(i)*AuER    ! chg AU's to km's
                  ENDDO

                magrsun = MAG(rsun)

                RS2= magrsun*magrsun
                RS3= RS2*magrsun
                Temp= DOT( R,RSun )
                Temp1= -GMS/RS3
                Temp2= 3.0D0*Temp/RS2
                APert(1)= Temp1 * (r(1) - Temp2*RSun(1))
                APert(2)= Temp1 * (r(2) - Temp2*RSun(2))
                APert(3)= Temp1 * (r(3) - Temp2*RSun(3))
              ENDIF

        ! -----------------  MOON Acceleration   ----------------------
        IF ( WhichOne .eq. 5 ) THEN
                CALL MOON( ITime,RMoon,MRtAsc,MDecl )
                magrmoon = MAG(rmoon)
                RM2= magRMoon**2
                RM3= RM2*magRMoon
                Temp= DOT( R,RMoon )
                Temp1= -GMM/RM3
                Temp2= 3.0D0*Temp/RM2
                APert(1)= Temp1 * (r(1) - Temp2*RMoon(1))
                APert(2)= Temp1 * (r(2) - Temp2*RMoon(2))
                APert(3)= Temp1 * (r(3) - Temp2*RMoon(3))
              ENDIF

        ! -----------------  Drag Acceleration   ----------------------
        IF ( WhichOne .eq. 6 ) THEN
                Va(1)= V(1) + (OmegaEarth*r(2))   ! km/s
                Va(2)= V(2) - (OmegaEarth*r(1))
                Va(3)= V(3)
                magva = MAG( Va )

                CALL ATMOS( R, Rho )

                Temp= -1000.0D0 * magVa * 0.5D0*Rho* ( 1.0D0/BC )*
     &                 6378137.0D0
                APert(1)= Temp*Va(1)
                APert(2)= Temp*Va(2)
                APert(3)= Temp*Va(3)
              ENDIF

        ! ----------------- Solar Acceleration   ----------------------
        IF ( WhichOne .eq. 7 ) THEN
                CALL SUN( ITime,RSun,SRtAsc,SDecl )

                Beta= 0.4D0                          ! reflectivity
                magAPert= (Beta*2.0D0*4.51D-06)/BC   ! assume Csr = 2.0D0
                Temp= -magAPert/magrsun
                APert(1)= Temp*RSun(1)
                APert(2)= Temp*RSun(2)
                APert(3)= Temp*RSun(3)
              ENDIF

        ! ------------------- Square Gravity Field --------------------
        IF ( WhichOne .eq. 10 ) THEN

            CALL FullGeop( R,V,ITime,WhichOne,BC,50,C,S,APert )
          ENDIF

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE PDERIV
*
*  this subroutine calculates the derivative of the state vector DO use with
*    the Runge-Kutta algorithm.  The DerivType string is used to determine
*    which perturbation equations are used.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    ITime       - Initial Time (Julian Date)     Days from 4713 BC
*    X           - State Vector                   km  ,  km/s
*    DerivType   - String of which perts to incl  'Y' .and. 'N'
*                  Options are in order : J2, J3,
*                  J4, SUN, MOON, Drag, Solarrad
*    BC          - Ballistic Coefficient          kg/m2
*
*  Outputs       :
*    XDot        - Derivative of State Vector     km/s, km/TU2
*
*  Locals        :
*    RCubed      - Radius vector cubed            ER3
*    Ro          - Radius vector                  km
*    Vo          - Velocity vector                km/s
*    APert       - Perturbing acceleration        km/TU2
*    TempPert    - Temporary acceleration         km/TU2
*    i           - Index
*
*  Coupling      :
*    PERTACCEL   - Calculates the actual values of each perturbing acceleration
*    ADDVEC      - Adds two vectors together
*    MAG         - Magnitude of a vector
*
*  References    :
*
* ------------------------------------------------------------------------------

      SUBROUTINE PDERIV( ITime,X,DerivType,BC,  XDot )
        IMPLICIT NONE
        REAL*8 X(6),XDot(6),ITime,BC
        CHARACTER*10 DerivType

* -----------------------------  Locals  ------------------------------
        REAL*8 RCubed,Ro(3),Vo(3),APert(3),TempPert(3), magr, magv, MAG
        INTEGER i

        EXTERNAL MAG
        
        ! --------------------  Implementation   ----------------------
        DO i= 1, 3
            APert(i)= 0.0D0
            Ro(i)   = X(i)
            Vo(i)   = X(i+3)
        ENDDO
        magr = MAG( Ro )
        magv = MAG( Vo )
*        APert(4)= 0.0D0
        RCubed = magr**3

* -------------------------  Velocity Terms  --------------------------
        XDot(1)= X(4)
        XDot(2)= X(5)
        XDot(3)= X(6)

* ----------------------  Acceleration Terms  -------------------------
        IF ( DerivType(1:1).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,1,BC, APert )
          ENDIF
        IF ( DerivType(2:2).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,2,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(3:3).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,3,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(4:4).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,4,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(5:5).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,5,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(6:6).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,6,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        IF ( DerivType(7:7).eq.'Y' ) THEN
            CALL PertAccel( Ro,Vo,ITime,7,BC, TempPert )
            CALL AddVec( TempPert,APert,APert )
          ENDIF
        ! -------------------- new full gravity field -----------------
        IF ( DerivType(10:10) .eq. 'Y' ) THEN
            CALL PERTACCEL( Ro,Vo,ITime,10,BC, TempPert )
            CALL ADDVEC( TempPert,APert,APert )
          ENDIF

        XDot(4)= (-X(1) / RCubed) + APert(1)
        XDot(5)= (-X(2) / RCubed) + APert(2)
        XDot(6)= (-X(3) / RCubed) + APert(3)

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                                SUBROUTINE RK4
*
*  this subroutine is a fourth order Runge-Kutta integrator DO a 6 dimension
*    First Order differential equation.  The intended use is DO a satellite
*    equation of motion.  The user must provide an external SUBROUTINE containing
*    the system Equations of Motion.  Notice time is included since some
*    applications in PDERIV may need this.  The LAST position in DerivType is a
*    flag DO two-body motion.  Two-Body motion is used IF ( the 10th element is
*    set to '2', otherwise the Yes .and. No values determine which perturbations
*    to use. Be careful with the units. The ITime parameter comes as a JD because
*    it's used DO SUN/MOON calcs later on. The stepsize is changed to s
*    operations involving the state vector because it is canonical.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*                                 12 Nov 1993 - fix DO s, DT, etc.
*  Inputs          Description                    Range / Units
*    ITime       - Initial Time (Julian Date)     Days from 4713 BC
*    DtDay       - Step size                      Day
*    XDot        - Derivative of State Vector
*    DerivType   - String of which perts to incl  'Y' .and. 'N'
*                  Options are in order : J2, J3,
*                  J4, SUN, MOON, Drag, Solarrad
*    BC          - Ballistic Coefficient          kg/m2
*    X           - State vector at initial time   km, km/s
*
*  Outputs       :
*    X           - State vector at new time       km, km/s
*
*  Locals        :
*    K           - Storage DO values of state
*                   vector at different times
*    Temp        - Storage DO state vector
*    TempTime    - Temporary time storage half
*                   way between DtDay             Day
*    J           - Index
*    DtSec        - Step size                      s
*
*  Coupling      :
*    DERIV       - SUBROUTINE for Derivatives of EOM
*    PDeriv      - SUBROUTINE for Perturbed Derivatives of EOM
*
*  References    :
*    Vallado       2007, 526
*
* ------------------------------------------------------------------------------
*
      SUBROUTINE RK4           ( ITime,DtDay,XDot,DerivType,BC,  X )
        IMPLICIT NONE
        REAL*8  DtDay, ITime, X(6), XDot(6), BC
        CHARACTER*10 DerivType

* -----------------------------  Locals  ------------------------------
        REAL*8 K(6,3), TEMP(6), TempTime, DtSec, TUDay
        INTEGER J

        TUDay = 0.00933809017716D0
          DtSec= DtDay/TUDay
          ! --------- Evaluate 1st Taylor Series Term -----------------
          IF (DerivType(10:10).eq.'2') THEN
              CALL DERIV( X,XDot )
            ELSE
              CALL PDERIV( ITime,X,DerivType,BC,XDot )
            ENDIF

          TempTime = ITime + DtDay*0.5D0

          ! -------- Evaluate 2nd Taylor Series Term ------------------
          DO J = 1,6
             K(J,1)  = DtSec * XDot(J)
             TEMP(J) = X(J) + 0.5D0*K(J,1)
          ENDDO

          IF (DerivType(10:10).eq.'2') THEN
              CALL DERIV( Temp,XDot )
            ELSE
              CALL PDERIV( TempTime,Temp,DerivType,BC,XDot )
            ENDIF

          ! ---------- Evaluate 3rd Taylor Series Term ----------------
          DO J = 1,6
             K(J,2)  = DtSec * XDot(J)
             TEMP(J) = X(J) + 0.5D0*K(J,2)
          ENDDO

          IF (DerivType(10:10).eq.'2') THEN
              CALL DERIV( Temp,XDot )
            ELSE
              CALL PDERIV( TempTime,Temp,DerivType,BC,XDot )
            ENDIF

          ! ---------- Evaluate 4th Taylor Series Term ----------------
          DO J = 1,6
             K(J,3)  = DtSec * XDot(J)
             TEMP(J) = X(J) + K(J,3)
          ENDDO

          IF (DerivType(10:10).eq.'2') THEN
              CALL DERIV( Temp,XDot )
            ELSE
              CALL PDERIV( ITime+DtDay,Temp,DerivType,BC,XDot )
            ENDIF

          ! ---- Update the state vector, perform integration  --------
          DO J = 1,6
             X(J) = X(J) + ( K(J,1) + 2.0D0*(K(J,2) + K(J,3)) +
     &                       DtSec*XDot(J) ) / 6.0D0
          ENDDO

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                                SUBROUTINE RKF45
*
*  this subroutine is a fourth order Runge-Kutta-Fehlberg integrator DO a 6-D
*    First Order differential equation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    ITime       - Initial Time (Julian Date)     Days from 4713 BC
*    DtDay       - Step size                      Day
*    XDot        - Derivative of State Vector
*    DerivType   - String of which perts to incl  'Y' .and. 'N'
*                  Options are in order : J2, J3,
*                  J4, SUN, MOON, Drag, Solarrad
*    BC          - Ballistic Coefficient          kg/m2
*    X           - State vector at initial time   km, km/s
*
*  Outputs       :
*    X           - State vector at new time       km, km/s
*
*  Locals        :
*    K           - Storage DO values of state
*                    vector at different times
*    Temp        - Storage DO state vector
*    TempTime    - Temporary time storage half
*                    way between DtDay            Day
*    J           - Index
*    DtSec        - Step size                      s
*
*  Coupling      :
*    DERIV       - SUBROUTINE DO Derivatives of EOM
*    PDeriv      - SUBROUTINE DO Perturbed Derivatives of EOM
*
*  References    :
*    Vallado       2007, 526
*
* ------------------------------------------------------------------------------  

      SUBROUTINE RKF45       ( ITime, DtDay, XDot,DerivType, BC, X )
        IMPLICIT NONE
        REAL*8 ITime,DtDay, XDot(6),X(6), BC
        CHARACTER*10 Derivtype

* -----------------------------  Locals  ------------------------------
        INTEGER Ktr, J
        REAL*8 K(6,6), Temp(6,1), DtSec, HMin, HMax, TStop, Time, Err,
     &         S, TempTime, Small, TUDay

        Small =     0.000001D0  ! this is pretty sensitive for RF45
        TUDay =     0.0093380913806D0
        HMin = DtDay/64.0D0
        HMax = DtDay*64.0D0
        Time = ITime
        TStop= ITime + DtDay
        DtSec = DtDay/TUDay

        Ktr= 1
        DO WHILE (Time .lt. TStop)
            IF ( Time + DtDay .gt. TStop ) THEN  ! Make sure you END exactly on the step
                DtDay= TStop - Time
              ENDIF

            ! ------------- Evaluate 1st Taylor Series Term -----------
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( X, XDot )
              ELSE
                CALL PDeriv( Time,X,DerivType,BC, XDot )
              ENDIF

            TempTime= Time + DtDay*0.25D0
            ! ------------- Evaluate 2nd Taylor Series Term -----------
            DO j= 1,6
                K(J,1)   = DtSec * XDot(J)   ! set # 1
                Temp(J,1)= X(J) + 0.25D0 * K(J,1)   !get ready DO 2
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( TempTime,Temp,DerivType,BC,  XDot )
              ENDIF

            TempTime= Time + DtDay*0.375D0
            ! ------------- Evaluate 3rd Taylor Series Term -----------
            DO j= 1,6
                K(J,2)   = DtSec * XDot(J)
                Temp(J,1)= X(J) + 0.09375D0 * K(J,1)
     &                     + 0.28125D0 * K(J,2)
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( TempTime,Temp,DerivType,BC,  XDot )
              ENDIF

            TempTime= Time + DtDay*12.0D0/13.0D0

            ! ------------- Evaluate 4th Taylor Series Term -----------
            DO j= 1,6
                  K(J,3)    = DtSec * XDot(J)
                  Temp(J,1) = X(J) + K(J,1) * 1932.0D0/2197.0D0
     &                        - K(J,2)*7200.0D0/2197.0D0
     &                        + K(J,3)*7296.0D0/2197.0D0
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( TempTime,Temp,DerivType,BC,  XDot )
              ENDIF

            ! ------------- Evaluate 5th Taylor Series Term -----------
            DO j= 1,6
                  K(J,4)    = DtSec * XDot(J)
                  Temp(J,1) = X(J) + K(J,1)* 439.0D0/ 216.0D0
     &                      - K(J,2) * 8.0D0 + K(J,3)*3680.0D0/ 513.0D0
     &                      - K(J,4) * 845.0D0/4104.0D0
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( Time+DtDay,Temp,DerivType,BC,  XDot )
              ENDIF 

            TempTime= Time + DtDay*0.5D0 

            ! ------------- Evaluate 6th Taylor Series Term -----------
            DO j= 1,6
                  K(J,5)    = DtSec * XDot(J)
                  Temp(J,1) = X(J) - K(J,1)*8.0D0/27.0D0
     &                      + K(J,2)* 2.0D0 - K(J,3)*3544.0D0/2565.0D0
     &                      + K(J,4)*1859.0D0/4104.0D0 - K(J,5)*0.275D0
              ENDDO
            IF ( DerivType(10:10) .eq. '2' ) THEN
                CALL DERIV( Temp, XDot )
              ELSE
                CALL PDeriv( TempTime,Temp,DerivType,BC,  XDot )
              ENDIF 

            DO j= 1,6
                  K(J,6)=  DtSec * XDot(J)
              ENDDO

            ! ------------------- Check DO convergence ---------------
            Err= 0.0D0 
            DO j= 1,6
                Err= DABS( K(J,1)*1.0D0/360.0D0
     &                 - K(J,3)*128.0D0/4275.0D0
     &                 - K(J,4)*2197.0D0/75240.0D0
     &                 + K(J,5)*0.02D0 + K(J,6)*2.0D0/55.0D0 )
              ENDDO

            ! ----- Update the State vector, perform integration ------
            IF ( ( Err .lt. Small ) .or.
     &           ( DtDay .le. 2.0D0*HMin+Small ) ) THEN
                DO j= 1,6
                    X(J)= X(J) + K(J,1)*25.0D0/216.0D0
     &                      + K(J,3)*1408.0D0/2565.0D0
     &                      + K(J,4)*2197.0D0/4104.0D0 - K(J,5)*0.2D0
                  ENDDO
                Time= Time + DtDay
                s   = 0.0D0
                Ktr = 1
              ELSE
                S= 0.84D0* (Small*DtDay/Err)**0.25D0
                IF ( ( S.lt.0.75D0 ).and.(DtDay.gt.2.0D0*HMin ) ) THEN  ! Reduce  Step  Size
                    DtDay= DtDay * 0.5D0
                  ENDIF
                IF ( ( S.gt.1.5D0 ).and.(2.0D0*DtDay.lt.HMax ) ) THEN    ! Increase Step Size
                    DtDay= DtDay * 2.0D0
                  ENDIF
                Ktr = Ktr + 1
              ENDIF

          ENDDO   ! WHILE

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                                SUBROUTINE COWELL
*
*  this subroutine uses a fourth order Runge-Kutta integrator on a 6 dimension
*    First Order differential equation.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - Initial position vector        km
*    V           - Initial velocity vector        km/s
*    ITime       - Initial Time (Julian Date)     Days from 4713 BC
*    FTime       - Final Time (Julian Date)       Days from 4713 BC
*    DtDay       - Step size                      Day
*    DerivType   - String of which perts to incl  'Y' .and. 'N'
*                  Options are in order : J2, J3,
*                  J4, SUN, MOON, Drag, Solarrad
*    BC          - Ballistic Coefficient          kg/m2
*
*  Outputs       :
*    R1          - Final position vector          km
*    V1          - Final velocity vector          km/s
*
*  Locals        :
*    Time        - Current time during the loop   Days from 4713 BC
*    X           - State vector at each time      km, km/s
*
*  Coupling      :
*    RK4         - Runge-Kutta algorithm
*    MAG         - Magnitude of a vector
*
*  References    :
*
* ------------------------------------------------------------------------------

      SUBROUTINE Cowell   ( R,V,ITime,FTime,DtDay,DerivType,BC, R1,V1 )
        IMPLICIT NONE
        REAL*8 R(3),V(3),ITime,FTime,DtDay,BC,R1(3),V1(3)
        CHARACTER*10 DerivType

* -----------------------------  Locals  ------------------------------
        REAL*8 Time, X(6), XDot(6),FT
        INTEGER i

        DO i= 1, 6
            IF ( i .le. 3 ) THEN
                X(i)= r(i)
              ELSE
                X(i)= v(i-3)
              ENDIF
          ENDDO

        ! --------- Loop through the time interval desired ------------
        Time= ITime
        Ft = FTime
        DO WHILE (Time .le. Ft)
            IF ( Time+DtDay .gt. Ft ) THEN
                DtDay = Ft - Time
                Ft = FTime - 1.0D0
              write(*,*) 'asifgy'
              ENDIF

            CALL RK4( Time,DtDay,XDot,DerivType,BC, X )
         write(*,'(6(f8.4))') (x(i), i=1,6)
*            CALL RKF45( Time,DtDay,XDot,DerivType,BC, X )

            Time = Time + DtDay
         write(*,*) Ft,' ',time,' ',dtday
          ENDDO

        ! ----------------- Update the state vector -------------------
        DO i= 1, 6
            IF ( i .le. 3 ) THEN
                r1(i)  =   X(i)
              ELSE
                v1(i-3)=   X(i)
              ENDIF
          ENDDO
*        CALL MAG( R1 )
*        CALL MAG( V1 )

      RETURN
      END
*
* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ATMOS
*
*  this subroutine finds the atmospheric density at an altitude above an
*    oblate earth given the position vector in the Geocentric Equatorial
*    frame.  The position vector is in km's .and. the density is in gm/cm**3.
*    DO certain applications, it may not be necessary to find the Hellp
*    exact height difference.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - IJK Position vector            km
*
*  Outputs       :
*    Rho         - Density                        kg/m**3
*
*  Locals        :
*    Hellp       - Height above ellipsoid         km
*    OldDelta    - Previous value of DeltaLat     rad
*    Latgd       - Geodetic Latitude              -Pi/2 to Pi/2 rad
*    SinTemp     - Sine of Temp
*    RhoNom      - Nominal density at particular alt      gm/cm**3
*    NextBaseAlt - Next Base Altitude
*    LastBaseAlt - Last Base Altitude
*    H           - Scale Height                   km
*    i           - index
*    AtmosFile   - File of data DO the
*                    exponential atmosphere
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*
*  References    :
*    Vallado       2007, 562, Ex 8-4
*
* ------------------------------------------------------------------------------

      SUBROUTINE ATMOS       ( R,Rho )
        IMPLICIT NONE
        Real*8 R(3),Rho
        EXTERNAL MAG
        INTEGER i
        REAL*8 Hellp, OldDelta, Latgd, SinTemp, c, Decl, Temp, H,
     &         RhoNom, BaseAlt, LastBaseAlt, Small,
     &         LastH, LastRhoNom, MAG, magr

        INCLUDE 'astconst.cmn'

        ! -------------------  Initialize values   --------------------
        Small      =     0.0000001D0
        OPEN( UNIT=25,File='atmosexp.dat',STATUS='OLD' )

        magr = MAG( R )
        Decl = DASIN( R(3) / magr )
        Latgd= Decl

        ! ---- Iterate to find Geocentric .and. Geodetic Latitude  ----
        Temp = DSQRT( R(1)*R(1) + R(2)*R(2) )
        i= 1
        OldDelta = Latgd*2.0D0
        DO WHILE ( ( DABS(OldDelta-Latgd).ge.Small ) .and. (i.lt.10) )
            OldDelta= Latgd
            SinTemp = DSIN( Latgd )
            c       = 1.0D0 / (DSQRT( 1.0D0-eeSqrd*SinTemp*SinTemp ))
            Latgd   = DATAN( (r(3)+c*eeSqrd*SinTemp)/Temp )
            i = i + 1
          ENDDO
        Hellp = ( (Temp/DCOS(Latgd)) - c ) * rekm

        IF ( i .ge. 10 ) THEN
            Write(*,*)  'IJKtoLatLon did NOT converge '
          ENDIF

        ! ---------- Determine density based on altitude --------------

        ! ---------- Set up Loop to READ through Input File -----------
        READ( 25,*,END=999) LastBaseAlt, LastRhoNom, LastH

        DO WHILE (LastBaseAlt .lt. Hellp )
            READ( 25,*,END=999) BaseAlt, RhoNom, H
            IF (BaseAlt .lt. Hellp ) THEN
                LastBaseAlt= BaseAlt
                LastRhoNom = RhoNom
                LastH      = H
              ENDIF
          ENDDO

        RHO   = LastRHONOM * EXP((LastBaseAlt-Hellp)/LastH)

 999    CLOSE( 25 )
      RETURN
      END
*
