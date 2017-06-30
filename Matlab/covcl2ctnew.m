% ----------------------------------------------------------------------------
%
%                           function covcl2ct
%
%  this function transforms a six by six covariance matrix expressed in classical elements
%    into one expressed in cartesian elements
%
%  author        : david vallado
%
%  revisions
%    vallado     - major update                                  26 aug 2015
%
%  inputs          description                    range / units
%    classcov    - 6x6 classical covariance matrix
%    classstate  - 6x1 classical orbit state      (a e i O w nu/m)
%    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
%
%  outputs       :
%    cartcov     - 6x6 cartesian covariance matrix
%    tm          - transformation matrix
%
%  locals        :
%    r           - matrix of partial derivatives
%    a           - semimajor axis                 km
%    ecc         - eccentricity
%    incl        - inclination                    0.0  to pi rad
%    omaga       - longitude of ascending node    0.0  to 2pi rad
%    argp        - argument of perigee            0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%    m           - mean anomaly                   0.0  to 2pi rad
%    p1,p2,p3,p4 - denominator terms for the partials
%    e0          - eccentric anomaly              0.0  to 2pi rad
%    true1, true2- temp true anomaly              0.0  to 2pi rad
%
%  coupling      :
%    newtonm     - newton iteration for m and ecc to nu
%    newtonnu    - newton iteration for nu and ecc to m
%    constastro
%
%  references    :
%    Vallado and Alfano 2015
%
%   [cartcov,tm] = covcl2ct( classcov,classstate,anom )
% ----------------------------------------------------------------------------

function [cartcov,tm] = covcl2ctnew( classcov, classstate, anom )

       % -------- define gravitational constant
        constastro;

        % --------- determine which set of variables is in use ---------
        % ---- parse the input vector into the classical elements -----
        a = classstate(1);
        n = sqrt(mum/a^3);
        ecc     = classstate(2);
        incl    = classstate(3);
        omega   = classstate(4);
        argp    = classstate(5);
        % -------- if mean anomaly is used, convert to true anomaly
        % -------- eccentric anomaly (e) is needed for both
        if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
            mean = classstate(6);
           [e, nu] = newtonm(ecc, mean);
        else
            if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            % note that mean is not used in the partials, but nu is! 
               nu = classstate(6);
               [e, mean] = newtonnu(ecc, nu);
            end
        end
        
        p = a*(1-ecc^2)/1000;  % needs to be in km
        [r,v] = coe2rv ( p,ecc,incl,omega,argp,nu,0.0,0.0,0.0 );
        rx = r(1)*1000;
        ry = r(2)*1000;
        rz = r(3)*1000;
        vx = v(1)*1000;
        vy = v(2)*1000;
        vz = v(3)*1000;
        sin_inc = sin(incl);
        cos_inc = cos(incl);
        sin_anode = sin(omega);
        cos_anode = cos(omega);
        sin_w = sin(argp);
        cos_w = cos(argp);
        sin_nu = sin(nu);
        cos_nu = cos(nu);
       
        p5 = (1.0 - ecc^2)^1.5 / ( (1.0 + ecc*cos(nu))^2 );  % dm/dv
        p6 = -sin(nu)*((ecc*cos(nu) + 1)*(ecc+cos(nu))/sqrt((ecc + cos(nu))^2) + 1.0 - 2.0*ecc^2 - ecc^3*cos(nu)) / ( (ecc*cos(nu) + 1.0)^2 * sqrt(1-ecc^2) );  % dm/de   
        
        % ---------------- calculate matrix elements ------------------
        % ---- partials of rx wrt (a e i O w m)
        drx_da = cos_anode*sin_w*sin_nu - cos_w*cos_anode*cos_nu + cos_inc*cos_w*sin_anode*sin_nu + cos_inc*cos_nu*sin_w*sin_anode;
        tm(1,1) = drx_da*(ecc*ecc - 1)/(1 + ecc*cos_nu);
        tm(1,1) = rx/a; 
        drx_decc = cos_anode*sin_w*sin_nu - cos_w*cos_anode*cos_nu + cos_inc*cos_w*sin_anode*sin_nu + cos_inc*cos_nu*sin_w*sin_anode;
        tm(1,2) = drx_decc*a*(cos_nu*ecc*ecc + 2*ecc + cos_nu)/((1 + ecc*cos_nu)*(1 + ecc*cos_nu));
        drx_dinc = sin_inc*cos_w*sin_anode*sin_nu + sin_inc*cos_nu*sin_w*sin_anode;
        tm(1,3) = drx_dinc*a*(1 - ecc*ecc)/(1 + ecc*cos_nu);
        drx_danode = cos_w*cos_nu*sin_anode - sin_w*sin_anode*sin_nu + cos_inc*cos_w*cos_anode*sin_nu + cos_inc*cos_anode*cos_nu*sin_w;
        tm(1,4) = drx_danode*a*(ecc*ecc - 1)/(1 + ecc*cos_nu);
        drx_dw = cos_w*cos_anode*sin_nu + cos_anode*cos_nu*sin_w + cos_inc*cos_w*cos_nu*sin_anode - cos_inc*sin_w*sin_anode*sin_nu;
        tm(1,5) = drx_dw*a*(ecc*ecc - 1)/(1 + ecc*cos_nu);
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            drx_dnu = ecc*cos_anode*sin_w + cos_w*cos_anode*sin_nu + cos_anode*cos_nu*sin_w + ecc*cos_inc*cos_w*sin_anode;
            drx_dnu = drx_dnu + cos_inc*cos_w*cos_nu*sin_anode - cos_inc*sin_w*sin_anode*sin_nu;
            tm(1,6) = drx_dnu*a*(ecc*ecc - 1)/((1 + ecc*cos_nu)*(1 + ecc*cos_nu));
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p0 = a * (ecc^2 - 1.0) / (ecc*cos(nu) + 1.0)^2;
                tm(1,6) = p0 * (ecc*cos(omega)*sin(argp)+cos(omega)*cos(argp)*sin(nu)+cos(omega)*sin(argp)*cos(nu)+ecc*cos(incl)*sin(omega)*cos(argp)+cos(incl)*sin(omega)*cos(argp)*cos(nu)-cos(incl)*sin(omega)*sin(argp)*sin(nu));
                tm(1,6) = tm(1,6) / p5;
                tm(1,2) = tm(1,2) - tm(1,6) * p6;
            end
        end
        
        % ---- partials of ry wrt (a e i O w nu/m)
        dry_da = cos_w*cos_nu*sin_anode - sin_w*sin_anode*sin_nu + cos_inc*cos_w*cos_anode*sin_nu + cos_inc*cos_anode*cos_nu*sin_w;
        tm(2,1) =  - dry_da*(ecc*ecc - 1)/(1 + ecc*cos_nu);
        tm(2,1) = ry/a; 
        dry_decc = cos_w*cos_nu*sin_anode - sin_w*sin_anode*sin_nu + cos_inc*cos_w*cos_anode*sin_nu + cos_inc*cos_anode*cos_nu*sin_w;
        tm(2,2) =  - dry_decc*a*(cos_nu*ecc*ecc + 2*ecc + cos_nu)/((1 + ecc*cos_nu)*(1 + ecc*cos_nu));
        dry_dinc = sin_inc*cos_w*cos_anode*sin_nu + sin_inc*cos_nu*sin_w*cos_anode;
        tm(2,3) = dry_dinc*a*(ecc*ecc - 1)/(1 + ecc*cos_nu);
        dry_danode = sin_w*sin_nu*cos_anode - cos_w*cos_anode*cos_nu + cos_inc*cos_w*sin_anode*sin_nu + cos_inc*sin_anode*cos_nu*sin_w;
        tm(2,4) = dry_danode*a*(ecc*ecc - 1)/(1 + ecc*cos_nu);
        dry_dw = cos_w*sin_anode*sin_nu + sin_anode*cos_nu*sin_w - cos_inc*cos_w*cos_nu*cos_anode + cos_inc*sin_w*cos_anode*sin_nu;
        tm(2,5) = dry_dw*a*(ecc*ecc - 1)/(1 + ecc*cos_nu);
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            dry_dnu = ecc*sin_anode*sin_w + cos_w*sin_anode*sin_nu + sin_anode*cos_nu*sin_w - ecc*cos_inc*cos_w*cos_anode;
            dry_dnu = dry_dnu - cos_inc*cos_w*cos_nu*cos_anode + cos_inc*sin_w*cos_anode*sin_nu;
            tm(2,6) = dry_dnu*a*(ecc*ecc - 1)/((1 + ecc*cos_nu)*(1 + ecc*cos_nu));
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p0 = a * (ecc^2 - 1.0) / (ecc*cos(nu) + 1.0)^2;
                tm(2,6) = p0 * (ecc*sin(omega)*sin(argp)+ sin(omega)*cos(argp)*sin(nu)+ sin(omega)*sin(argp)*cos(nu)- ecc*cos(incl)*cos(omega)*cos(argp)- cos(incl)*cos(omega)*cos(argp)*cos(nu)+ cos(incl)*cos(omega)*sin(argp)*sin(nu));
                tm(2,6) = tm(2,6) / p5;
                tm(2,2) = tm(2,2) - tm(2,6) * p6;
           end
        end
               
        % ---- partials of rz wrt (a e i O w nu/m)
        drz_da = sin_inc*cos_w*sin_nu + sin_inc*cos_nu*sin_w;
        tm(3,1) = drz_da*(1 - ecc*ecc)/(1 + ecc*cos_nu);
        tm(3,1) = rz/a; 
        drz_decc = cos_w*sin_nu + cos_nu*sin_w;
        tm(3,2) =  - drz_decc*a*sin_inc*(cos_nu*ecc*ecc + 2*ecc + cos_nu)/((1 + ecc*cos_nu)*(1 + ecc*cos_nu));
        drz_dinc = cos_inc*cos_w*sin_nu + cos_inc*cos_nu*sin_w;
        tm(3,3) =  - drz_dinc*a*(ecc*ecc - 1)/(1 + ecc*cos_nu);
        tm(3,4) = 0;
        drz_dw = sin_inc*cos_w*cos_nu - sin_inc*sin_w*sin_nu;
        tm(3,5) =  - drz_dw*a*(ecc*ecc - 1)/(1 + ecc*cos_nu);
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            drz_dnu = cos_w*cos_nu - sin_w*sin_nu + ecc*cos_w;
            tm(3,6) = drz_dnu*a*sin_inc*(1 - ecc*ecc)/((1 + ecc*cos_nu)*(1 + ecc*cos_nu));
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p0 = -a * (ecc^2 - 1.0) / (ecc*cos(nu) + 1.0)^2;
                tm(3,6) = p0 * sin(incl)*(cos(argp+nu)+ecc*cos(argp));
                tm(3,6) = tm(3,6) / p5;
                tm(3,2) = tm(3,2) - tm(3,6) * p6;
            end
        end
                
        % ---- partials of vx wrt (a e i O w nu/m)
        sqrt_term = sqrt(mum/(a*(1 - ecc*ecc)));
        dvx_da = ecc*cos_anode*sin_w + cos_w*cos_anode*sin_nu + cos_anode*cos_nu*sin_w + ecc*cos_inc*cos_w*sin_anode;
        dvx_da = dvx_da + cos_inc*cos_w*cos_nu*sin_anode - cos_inc*sin_w*sin_anode*sin_nu;
        tm(4,1) = dvx_da*mum/(2*a*a*(1 - ecc*ecc)*sqrt_term);
        tm(4,1) = -vx/(2.0*a);  
        dvx_decc = ecc*cos_anode*sin_w + cos_w*cos_anode*sin_nu + cos_anode*cos_nu*sin_w + ecc*cos_inc*cos_w*sin_anode;
        dvx_decc = dvx_decc + cos_inc*cos_w*cos_nu*sin_anode - cos_inc*sin_w*sin_anode*sin_nu;
        dvx_decc =  - dvx_decc*mum*ecc/(a*(1 - ecc*ecc)*(1 - ecc*ecc)*sqrt_term);
        tm(4,2) = dvx_decc - (cos_anode*sin_w + cos_inc*cos_w*sin_anode)*sqrt_term;
        dvx_dinc = cos_w*cos_nu - sin_w*sin_nu + ecc*cos_w;
        tm(4,3) = dvx_dinc*sin_inc*sin_anode*sqrt_term;
        dvx_danode = ecc*sin_w*sin_anode + cos_w*sin_anode*sin_nu + cos_nu*sin_w*sin_anode - ecc*cos_inc*cos_w*cos_anode;
        dvx_danode = dvx_danode - cos_inc*cos_w*cos_anode*cos_nu + cos_inc*cos_anode*sin_w*sin_nu;
        tm(4,4) = dvx_danode*sqrt_term;
        dvx_dw = cos_anode*sin_w*sin_nu - cos_w*cos_anode*cos_nu - ecc*cos_w*cos_anode + ecc*cos_inc*sin_w*sin_anode;
        dvx_dw = dvx_dw + cos_inc*cos_w*sin_anode*sin_nu + cos_inc*cos_nu*sin_w*sin_anode;
        tm(4,5) = dvx_dw*sqrt_term;
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            dvx_dnu = cos_anode*sin_w*sin_nu - cos_w*cos_anode*cos_nu + cos_inc*cos_w*sin_anode*sin_nu + cos_inc*cos_nu*sin_w*sin_anode;
            tm(4,6) = dvx_dnu*sqrt_term;
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p0 = sqrt(-mum/(a*(ecc^2-1.0)));
                tm(4,6) = p0 * (cos(omega)*sin(argp)*sin(nu)-cos(omega)*cos(argp)*cos(nu)+cos(incl)*sin(omega)*cos(argp)*sin(nu)+cos(incl)*sin(omega)*sin(argp)*cos(nu));
                tm(4,6) = tm(4,6) / p5;
                tm(4,2) = tm(4,2) - tm(4,6) * p6;
            end
        end
        
        % ---- partials of vy wrt (a e i O w nu/m)
        dvy_da = ecc*sin_anode*sin_w + cos_w*sin_anode*sin_nu + sin_anode*cos_nu*sin_w - ecc*cos_inc*cos_w*cos_anode;
        dvy_da = dvy_da - cos_inc*cos_w*cos_nu*cos_anode + cos_inc*sin_w*cos_anode*sin_nu;
        tm(5,1) = dvy_da*mum/(2*a*a*(1 - ecc*ecc)*sqrt_term);
        tm(5,1) = -vy/(2.0*a);  
        dvy_decc = ecc*sin_anode*sin_w + cos_w*sin_anode*sin_nu + sin_anode*cos_nu*sin_w - ecc*cos_inc*cos_w*cos_anode;
        dvy_decc = dvy_decc - cos_inc*cos_w*cos_nu*cos_anode + cos_inc*sin_w*cos_anode*sin_nu;
        dvy_decc =  - dvy_decc*mum*ecc/(a*(1 - ecc*ecc)*(1 - ecc*ecc)*sqrt_term);
        tm(5,2) = dvy_decc - (sin_anode*sin_w - cos_inc*cos_w*cos_anode)*sqrt_term;
        dvy_dinc = ecc*sin_inc*cos_w*cos_anode + sin_inc*cos_w*cos_anode*cos_nu - sin_inc*cos_anode*sin_w*sin_nu;
        tm(5,3) =  - dvy_dinc*sqrt_term;
        dvy_danode = ecc*cos_anode*sin_w + cos_w*cos_anode*sin_nu + cos_anode*cos_nu*sin_w + ecc*cos_inc*cos_w*sin_anode;
        dvy_danode = dvy_danode + cos_inc*cos_w*cos_nu*sin_anode - cos_inc*sin_w*sin_anode*sin_nu;
        tm(5,4) =  - dvy_danode*sqrt_term;
        dvy_dw = ecc*sin_anode*cos_w + cos_w*sin_anode*cos_nu - sin_anode*sin_nu*sin_w + ecc*cos_inc*sin_w*cos_anode;
        dvy_dw = dvy_dw + cos_inc*cos_w*sin_nu*cos_anode + cos_inc*sin_w*cos_anode*cos_nu;
        tm(5,5) =  - dvy_dw*sqrt_term;
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            dvy_dnu = sin_anode*cos_w*cos_nu - sin_w*sin_anode*sin_nu + cos_inc*cos_w*cos_anode*sin_nu + cos_inc*cos_nu*sin_w*cos_anode;
            tm(5,6) =  - dvy_dnu*sqrt_term;
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p0 = sqrt(-mum/(a*(ecc^2-1.0)));
                tm(5,6) = -p0 * (sin(omega)*cos(argp)*cos(nu)-sin(omega)*sin(argp)*sin(nu)+cos(incl)*cos(omega)*cos(argp)*sin(nu)+cos(incl)*cos(omega)*sin(argp)*cos(nu));
                tm(5,6) = tm(5,6) / p5;
                tm(5,2) = tm(5,2) - tm(5,6) * p6;
            end
        end
        
        % ---- partials of vz wrt (a e i O w nu/m)
        dvz_da = cos_w*cos_nu - sin_w*sin_nu + ecc*cos_w;
        tm(6,1) = dvz_da*mum*sin_inc/(2*a*a*(ecc*ecc - 1)*sqrt_term);
        tm(6,1) = -vz/(2.0*a);  
        dvz_decc = cos_w*cos_nu - sin_w*sin_nu + ecc*cos_w;
        dvz_decc = dvz_decc*mum*ecc*sin_inc/(a*(1 - ecc*ecc)*(1 - ecc*ecc)*sqrt_term);
        tm(6,2) = dvz_decc + sin_inc*cos_w*sqrt_term;
        tm(6,3) = cos_inc*(cos_w*cos_nu - sin_w*sin_nu + ecc*cos_w)*sqrt_term;
        tm(6,4) = 0;
        tm(6,5) =  -sin_inc*(cos_w*sin_nu + cos_nu*sin_w + ecc*sin_w)*sqrt_term;
        if strcmp(anom,'truea') == 1 || strcmp(anom,'truen') == 1  % 1 is true
            tm(6,6) =  -sin_inc*(cos_w*sin_nu + cos_nu*sin_w)*sqrt_term;
        else
            if strcmp(anom,'meana') == 1 || strcmp(anom,'meann') == 1  % 1 is true
                p0 = sqrt(-mum/(a*(ecc^2-1.0)));
                tm(6,6) = p0 * (-sin(incl)*sin(argp+nu));
                tm(6,6) = tm(6,6) / p5;
                tm(6,2) = tm(6,2) - tm(6,6) * p6;
            end
        end
        
        % ---------- calculate the output covariance matrix -----------
        cartcov =  tm*classcov*tm';

