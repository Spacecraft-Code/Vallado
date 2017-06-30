% ----------------------------------------------------------------------------
%
%                           function truemean
%
%  this function forms the transformation matrix to go between the
%    norad true equator mean equinox of date and the mean equator mean equinox
%    of date (eci).  the results approximate the effects of nutation and
%    precession.
%
%  author        : david vallado                  719-573-2600   25 jun 2002
%
%  revisions
%    vallado     - fixes to order                                29 sep 2002
%    vallado     - fixes to all options                           6 may 2003
%
%  inputs          description                    range / units
%    ttt         - julian centuries of tt
%    order       - number of terms for nutation   4, 50, 106, ...
%    eqeterms    - number of terms for eqe        0, 2
%    opt         - option for processing          a - complete nutation
%                                                 b - truncated nutation
%                                                 c - truncated transf matrix
%
%  outputs       :
%    nutteme     - matrix for mod - teme - an approximation for nutation
%
%  locals        :
%    prec        - matrix for mod - j2000
%    tm          - combined matrix for teme
%    ttt2        - ttt squared
%    ttt3        - ttt cubed
%    l           -                                rad
%    ll          -                                rad
%    f           -                                rad
%    d           -                                rad
%    omega       -                                rad
%    deltapsi    - nutation angle                 rad
%    deltaeps    - change in obliquity            rad
%    eps         - mean obliquity of the ecliptic rad
%    trueeps     - true obliquity of the ecliptic rad
%    meaneps     - mean obliquity of the ecliptic rad
%
%  coupling      :
%
%
%  references    :
%    vallado       2004, 230
%
% [deltapsi,trueeps,meaneps,omega,eqe,nutteme] = truemean ( ttt,order,eqeterms,opt );
% ----------------------------------------------------------------------------

function [deltapsi,trueeps,meaneps,omega,eqe,nutteme] = truemean ( ttt,order,eqeterms,opt );

        deg2rad = pi/180.0;

        [iar80,rar80] = iau80in;

        % ---- determine coefficients for iau 1980 nutation theory ----
        ttt2= ttt*ttt;
        ttt3= ttt2*ttt;
        ttt4= ttt2*ttt2;

        meaneps = -46.8150 *ttt - 0.00059 *ttt2 + 0.001813 *ttt3 + 84381.448;
        meaneps = rem( meaneps/3600.0 ,360.0  );
        meaneps = meaneps * deg2rad;

        l    =  134.96340251  + ( 1717915923.2178 *ttt + ...
                31.8792 *ttt2 + 0.051635 *ttt3 - 0.00024470 *ttt4 ) / 3600.0;
        l1   =  357.52910918  + (  129596581.0481 *ttt - ...
                 0.5532 *ttt2 - 0.000136 *ttt3 - 0.00001149*ttt4 )  / 3600.0;
        f    =   93.27209062  + ( 1739527262.8478 *ttt - ...
                12.7512 *ttt2 + 0.001037 *ttt3 + 0.00000417*ttt4 )  / 3600.0;
        d    =  297.85019547  + ( 1602961601.2090 *ttt - ...
                 6.3706 *ttt2 + 0.006593 *ttt3 - 0.00003169*ttt4 )  / 3600.0;
        omega=  125.04455501  + (   -6962890.2665 *ttt + ...
                 7.4722 *ttt2 + 0.007702 *ttt3 - 0.00005939*ttt4 )  / 3600.0;

        l    = rem( l,360.0  )     * deg2rad;
        l1   = rem( l1,360.0  )    * deg2rad;
        f    = rem( f,360.0  )     * deg2rad;
        d    = rem( d,360.0  )     * deg2rad;
        omega= rem( omega,360.0  ) * deg2rad;

        deltapsi= 0.0;
        deltaeps= 0.0;

        for i= 1:order   % the eqeterms in nut80.dat are already sorted
            tempval= iar80(i,1)*l + iar80(i,2)*l1 + iar80(i,3)*f + ...
                     iar80(i,4)*d + iar80(i,5)*omega;
            deltapsi= deltapsi + (rar80(i,1)+rar80(i,2)*ttt) * sin( tempval );
            deltaeps= deltaeps + (rar80(i,3)+rar80(i,4)*ttt) * cos( tempval );
          end

        % --------------- find nutation parameters --------------------
        deltapsi = rem( deltapsi,360.0  ) * deg2rad;
        deltaeps = rem( deltaeps,360.0  ) * deg2rad;
        trueeps  = meaneps + deltaeps;
fprintf(1,'dpsi %16.9f   deps %16.9f  trueps %16.8f meaneps %16.8f degpre \n',deltapsi/deg2rad,deltaeps/deg2rad,trueeps/deg2rad, (trueeps-deltaeps)/deg2rad);
        cospsi  = cos(deltapsi);
        sinpsi  = sin(deltapsi);
        coseps  = cos(meaneps);
        sineps  = sin(meaneps);
        costrueeps = cos(trueeps);
        sintrueeps = sin(trueeps);

        jdttt = ttt*36525.0 + 2451545.0;
        % small disconnect with ttt instead of ut1
        if (jdttt > 2450449.5 ) & (eqeterms > 0)
            eqe= deltapsi* cos(meaneps) ...
                + 0.00264*pi /(3600*180)*sin(omega) ...
                + 0.000063*pi /(3600*180)*sin(2.0 *omega);
          else
            eqe= deltapsi* cos(meaneps);
          end

        nut(1,1) =  cospsi;
        nut(1,2) =  costrueeps * sinpsi;
        if (opt=='b')
            nut(1,2) = 0.0;
          end;
        nut(1,3) =  sintrueeps * sinpsi;
        nut(2,1) = -coseps * sinpsi;
        if (opt=='b')
            nut(2,1) = 0.0;
          end;
        nut(2,2) =  costrueeps * coseps * cospsi + sintrueeps * sineps;
        nut(2,3) =  sintrueeps * coseps * cospsi - sineps * costrueeps;
        nut(3,1) = -sineps * sinpsi;
        nut(3,2) =  costrueeps * sineps * cospsi - sintrueeps * coseps;
        nut(3,3) =  sintrueeps * sineps * cospsi + costrueeps * coseps;

        st(1,1) =  cos(eqe);
        st(1,2) = -sin(eqe);
        st(1,3) =  0.0;
        st(2,1) =  sin(eqe);
        st(2,2) =  cos(eqe);
        st(2,3) =  0.0;
        st(3,1) =  0.0;
        st(3,2) =  0.0;
        st(3,3) =  1.0;

        nutteme = nut*st;

        if (opt=='c')
            nutteme(1,1) =  1.0;
            nutteme(1,2) =  0.0;
            nutteme(1,3) =  deltapsi * sineps;
            nutteme(2,1) =  0.0;
            nutteme(2,2) =  1.0;
            nutteme(2,3) =  deltaeps;
            nutteme(3,1) = -deltapsi * sineps;
            nutteme(3,2) = -deltaeps;
            nutteme(3,3) =  1.0;
          end;

%        tm = nutteme * prec;

