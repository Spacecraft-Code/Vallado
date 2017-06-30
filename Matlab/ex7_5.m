%     -----------------------------------------------------------------
%
%                              Ex7_5.m
%
%  this file demonstrates example 7-5.
%
%                          companion code for
%             fundamentals of astrodynamics and applications
%                                 2007
%                            by david vallado
%
%     (w) 719-573-2600, email dvallado@agi.com
%
%     *****************************************************************
%
%  current :
%            30 mar 07  david vallado
%                         original
%  changes :
%            13 feb 07  david vallado
%                         original baseline
%
%     *****************************************************************
    fid = 1;
constastro;
    
  % targeting example - target fixed
    ro = [-6518.1083, -2403.8479, -22.1722];
    vo = [ 2.604057, -7.105717, -0.263218];
    r1 = [6697.4756, 1794.5832, 0.0];
    v1 = [-1.962373, 7.323674, 0.0];
    dtsec = (4010)*60.0;
% dtsec = (13)*60.0;  % in sec from min,  4565
% ro = [0.5, 0.6, 0.7]*re;
% r1 = [0.0, 1.0, 0.0]*re;
    dm = 's';
 nrev = 0;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r1 );
    [vodv,v1dv,errorl] = lambertu ( ro,r1, dm,nrev, dtsec,fid );
    fprintf(1,' vodv %16.8f%16.8f%16.8f\n',vodv );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f\n',v1dv );
    
%     % original orbits
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vo, re, mu);
%     fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
%         p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'period %16.8f min \n',period/60);
%     fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',ro,vo,a );
% 
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1, re, mu);
%     fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
%         p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'period %16.8f min \n',period/60);
%     fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',r1,v1,a );
% 
%     % transfer orbit
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vodv, re, mu);
%     fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
%         p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'period %16.8f min \n',period/60);
%     fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',ro,vodv,a );    
    
    
    fprintf(1,'\n-------- lambertu test book pg 467 \n' );
    ro = [ 2.500000,    0.000000 ,   0.000000]*6378.137;
    r  = [ 1.9151111,   1.6069690,   0.000000]*6378.137;
    dtsec = 76.0*60.0;
    dm = 's';
    nrev = 0;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r );
    [vo,v,errorl] = lambertu ( ro,r, dm,nrev, dtsec,fid );
    fprintf(1,' vo %16.8f%16.8f%16.8f\n',vo );
    fprintf(1,' v  %16.8f%16.8f%16.8f\n',v );
    fprintf(1,'\n-------- lambertu test \n' );
    [vodv,v1dv,errorl] = lambertu ( ro,r, dm,nrev, dtsec,fid );
    fprintf(1,' vodv %16.8f%16.8f%16.8f\n',vodv );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f\n',v1dv );
    % transfer orbit
%     [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vodv, re, mu);
%     fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
%         p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
%     period = 2.0 * pi*sqrt(a^3/mu);
%     fprintf(1,'period %16.8f min \n',period/60);
%     fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',ro,vodv,a );    

    
    fprintf(1,'\n-------- lambertu test book pg 467 long way \n' );
    ro = [ 2.500000,    0.000000 ,   0.000000]*6378.137;
    r  = [ 1.9151111,   1.6069690,   0.000000]*6378.137;
    dtsec = 76.0*60.0;
    dm = 's';
    nrev = 0;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r );
    fprintf(1,'\n-------- lambertu test \n' );
    [vodv,v1dv,errorl] = lambertu ( ro,r, 'l',nrev, dtsec,fid );
    fprintf(1,' vodv %16.8f%16.8f%16.8f\n',vodv );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f\n',v1dv );
    % original orbit, assume circular
    vo = [0 0 0];
    vo(2) = sqrt(mu/ro(1));
   
    fprintf(1,' vo %16.8f%16.8f%16.8f\n',vo );

    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vo, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
    fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',ro,vo,a );    
    % transfer orbit
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vodv, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
    fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',ro,vodv,a );    

    [vodv,v1dv,errorl] = lambertu ( ro,r, 's',nrev, dtsec,fid );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vodv, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
    fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',ro,vodv,a );    

    pause;

    fprintf(1,'\n-------- lambertb test book pg 467 \n' );
    ro = [ 2.500000    0.000000    0.000000]*6378.137;
    r  = [ 1.9151111   1.6069690   0.000000]*6378.137;
    dtsec = 76.0*60.0;
    dm = 's';
    nrev = 0;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r );
    [vo,v,errorl] = lambertb ( ro,r, dm,nrev, dtsec );
    fprintf(1,' vo %16.8f%16.8f%16.8f\n',vo );
    fprintf(1,' v  %16.8f%16.8f%16.8f\n',v );
    fprintf(1,'\n-------- lambertu test \n' );
    [vodv,v1dv,errorl] = lambertu ( ro,r, dm,nrev, dtsec,fid );
    fprintf(1,' vodv %16.8f%16.8f%16.8f\n',vodv );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f\n',v1dv );

    fprintf(1,'\n-------- lambertb test book pg 467 long way \n' );
    ro = [ 2.500000    0.000000    0.000000]*6378.137;
    r  = [ 1.9151111   1.6069690   0.000000]*6378.137;
    dtsec = 76.0*60.0;
    dm = 's';
    nrev = 0;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r );
    [vo,v,errorl] = lambertb ( ro,r, 'l',nrev, dtsec );
    fprintf(1,' vo %16.8f%16.8f%16.8f\n',vo );
    fprintf(1,' v  %16.8f%16.8f%16.8f\n',v );
    fprintf(1,'\n-------- lambertu test \n' );
    [vodv,v1dv,errorl] = lambertu ( ro,r, 'l',nrev, dtsec,fid );
    fprintf(1,' vodv %16.8f%16.8f%16.8f\n',vodv );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f\n',v1dv );

    pause;
  
    dm = 's';
    nrev = 0;
    ro = [-8518.1083, -2403.8479, -2132.1722];
  %  ro = [23518.1083, 12403.8479, -232.1722];
    vo = [ 2.604057, -6.105717, -1.263218];
    r1 = [8697.4756, 1794.5832, 800.0];
    r1 = [2697.4756, 8794.5832, 1800.0];  % new
    r1 = [-8518.4756, -2400.0000, 1800.0];  % new
    v1 = [-1.962373, 6.323674, 1.02];
    
%     ro = [-18518.1083, -2403.8479, -2132.1722];
%     vo = [ 2.604057, -4.105717, -1.263218];
%     r1 = [8697.4756, 21794.5832, 800.0];
%     v1 = [-1.962373, 4.323674, 1.02];
%     
%     ro = [-18518.1083, -22403.8479, -2132.1722];
%     vo = [ 2.604057, -2.105717, -1.263218];
%     r1 = [28697.4756, 21794.5832, 800.0];
%     v1 = [-1.962373, 2.323674, 1.02];
   % [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vo, re, mu);
   % plotorb(ro,vo,a);
   % hold on;
   % [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1, re, mu);
   % plotorb(r1,v1,a);

    directory = 'd:\codes\library\matlab\';
    outfile = fopen(strcat(directory,'tlambert.out'), 'wt');
  
    % test what the iteratinos looklike on this case
%     [vodv,v1dv,errorl] = lambertu ( ro,r, dm,1, 22000,outfile );
%    
%     [vodv,v1dv,errorl] = lambertu ( ro,r, dm,1, 20000,outfile );
     [vodv,v1dv,errorl] = lambertu ( ro,r, dm,1, 21182,outfile );
     [vodv,v1dv,errorl] = lambertu ( ro,r, dm,1, 21183,outfile );
     [vodv,v1dv,errorl] = lambertu ( ro,r, dm,1, 21184,outfile );
     [vodv,v1dv,errorl] = lambertu ( ro,r, dm,1, 21185,outfile );
%     [vodv,v1dv,errorl] = lambertu ( ro,r, dm,1, 21700,outfile );

psi0 = (39+157)/2;
psiL = ((39+157)/2)*0.4;

[c2,c3] = findc2c3( psi0 );
c2dot0 = 1.0/(2*psi0) * (1.0 - psi0*c3 - 2.0*c2);
c2ddot0 = -2.0/psi0^2 + (4.0 * c2)/psi0^2;

[c2,c3] = findc2c3( psiL );
c2dotL = 1.0/(2*psiL) * (1.0 - psiL*c3 - 2.0*c2);
c2ddotL = -2.0/psiL^2 + (4.0 * c2)/psiL^2;

for dd = 1: 10
    psi = (psi0+psiL)*0.5;    
    % next guess
    if c2dot0 < c2dotL
    end;
        
    if c2dotL < c2dot0
    end;

    psi0 = psi - c2dot/c2ddot;
    fprintf(1,'%12.5f %12.5f %12.5f %12.5f %12.5f \n', psi, psi0, c2, c2dot, c2ddot);
end
pause;


    fprintf(1,'now do findtbi calcs \n');
    tbi = findtbi(ro, r, dm)
    %tbi = findlambertmins( ro, r, dm ); 
    %tbi = findlambertminsb( ro, r, dm ); 
    pause;
% lps, nrev dm,dtnew,           magro,      magr,         vara,          y,             x
%C 11   1   s 21200.0000000 9103.9997818 15945.3423626 4059.6664884 24754.6525924 1038.4818222   86.9010942 
%C 11   1   s 22000.0000000 9103.9997818 15945.3423626 4059.6664884 25069.7389723 1055.9930749   88.9604236 
%  50   1   s 762069.3486263 9103.9997818 15945.3423626 4059.6664884 30393.8496620 3988.5688838  139.6673768 
%
%fprintf( 1,'   conv     #s nrev  dm      y         xold         dtnew        upper      lower      psinew   slope  \n');
%    fprintf( outfile,'   conv     #s nrev  dm      y         xold         dtnew         vx       vy     vz     v1x    v1y     v1z     upper      lower      psinew   slope  \n');

    for zz = 20: 20:1000
        dtsec = zz * 60.0;
        [vodv,v1dv,errorl] = lambertu ( ro,r, dm,nrev, dtsec,outfile );
    end
    fprintf(1,'  \n');
    dm = 's';
    nrev = 0;
    for zz = 270: 5:1000  % 7000
        dtsec = zz * 60.0;
        if dtsec > 20000  % these are the min points
            nrev = 1;
        end
        if dtsec > 33900
            nrev = 2;
        end
        if dtsec > 131610
            nrev = 3;
        end
        if dtsec > 139853
            nrev = 4;
        end
        [vodv,v1dv,errorl] = lambertu ( ro,r, dm,nrev, dtsec,outfile );

        %[p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1, re, mu);
        %plotorb(ro,vodv,a);
        %pause;
    end
    fprintf(1,'  \n');

    dm = 'l';
    nrev = 0;
    for zz = 270: 5:1000  % 7000
        dtsec = zz * 60.0;
        if dtsec > 20000
            nrev = 1;
        end
        if dtsec > 33900
            nrev = 2;
        end
        if dtsec > 131610
            nrev = 3;
        end
        if dtsec > 139853
            nrev = 4;
        end
        [vodv,v1dv,errorl] = lambertu ( ro,r, dm,nrev, dtsec,outfile );

        %[p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1, re, mu);
        %plotorb(ro,vodv,a);
        %pause;
    end

    pause;

    fprintf(1,'\n-------- lambertb test ben joseph \n' );
    ro = [ 6822.88933   -5147.86167    -454.39488];
    r  = [ -4960.67860   10585.2504   927.1937739];
    dtsec = 4976.002;
    dm = 's';
    nrev = 0;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r );
    [vo,v,errorl] = lambertb ( ro,r, dm,nrev, dtsec );
    fprintf(1,' vo %16.8f%16.8f%16.8f\n',vo );
    fprintf(1,' v  %16.8f%16.8f%16.8f\n',v );
    fprintf(1,'\n-------- lambertu test \n' );
    [vodv,v1dv,errorl] = lambertu ( ro,r, dm,nrev, dtsec,fid );
    fprintf(1,' vodv %16.8f%16.8f%16.8f\n',vodv );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f\n',v1dv );

    pause;

    fprintf(1,'\n-------- lambertu test giuseppe \n' );
    mu         = 398600.4418;      % km3/s2
    re = 6378.137;
    velkmps = sqrt(mu / 6378.137);
    rad = 180.0/pi;

    ro = [ 0.500000,   0.600000,   0.700000]*re;
    r  = [ 0.000000,   1.000000,   0.000000]*re;
    dtsec = 779.999;
    dm = 's';
    nrev = 0;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r );
    [vo,v,errorl] = lambertu ( ro,r, dm,nrev, dtsec,fid );
    fprintf(1,' vo %16.8f%16.8f%16.8f\n',vo );
    fprintf(1,' v  %16.8f%16.8f%16.8f\n',v );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vo, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r,v, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);

    [vo,v,errorl] = lambertb ( ro,r, 's',nrev, dtsec );
    fprintf(1,' vo %16.8f%16.8f%16.8f\n',vo );
    fprintf(1,' v  %16.8f%16.8f%16.8f\n',v );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vo, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r,v, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
   
    fprintf(1,'\n-------- lambertu test multi-rev \n' );
    ro = [ 2.500000,    0.000000 ,   0.000000]*6378.137;
    r  = [ 1.9151111,   1.6069690,   0.000000]*6378.137;
   dtsec = 50000; % put in between, on short way side
    dm = 's';
    nrev = 1;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r );
    [vo,v,errorl] = lambertu ( ro,r, dm,nrev, dtsec,fid );
    fprintf(1,' vo %16.8f%16.8f%16.8f\n',vo );
    fprintf(1,' v  %16.8f%16.8f%16.8f\n',v );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vo, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r,v, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
 
    fprintf(1,'\n-------- lambertu test multi-rev \n' );
    ro = [ 22323.479500,         0.000000,     0.000000];
    r  = [ -7820.961704,     26780.193936,     0.000000];
    dtsec = 90000.0;
    dm = 's';
    nrev = 2;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r );
    [vo,v,errorl] = lambertu ( ro,r, dm,nrev, dtsec,fid );
    fprintf(1,' vo %16.8f%16.8f%16.8f\n',vo );
    fprintf(1,' v  %16.8f%16.8f%16.8f\n',v );
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper] = rv2coeh (ro,vo, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper] = rv2coeh (r,v, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);


    % targeting example - target fixed
    ro = [-6518.1083, -2403.8479, -22.1722];
    vo = [ 2.604057, -7.105717, -0.263218];
    r1 = [6697.4756, 1794.5832, 0.0];
    v1 = [-1.962373, 7.323674, 0.0];
    dtsec = 9000;  % 4010
    dm = 's';
    nrev = 1;
    fprintf(1,' ro %16.8f%16.8f%16.8f\n',ro );
    fprintf(1,' r  %16.8f%16.8f%16.8f\n',r1 );
    [vodv,v1dv,errorl] = lambertu ( ro,r1, dm,nrev, dtsec,fid );
    fprintf(1,' vodv %16.8f%16.8f%16.8f\n',vodv );
    fprintf(1,' v1dv %16.8f%16.8f%16.8f\n',v1dv );
    
    % original orbits
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vo, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
    fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',ro,vo,a );

    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (r1,v1, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
    fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',r1,v1,a );

    % transfer orbit
    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coeh (ro,vodv, re, mu);
    fprintf(1,'ans coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n',...
        p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );
    period = 2.0 * pi*sqrt(a^3/mu);
    fprintf(1,'period %16.8f min \n',period/60);
    fprintf(1,'plotorb([%16.8f,%16.8f,%16.8f], [%16.8f,%16.8f,%16.8f], %16.8f)\n',ro,vodv,a );
    
    
    
    
    