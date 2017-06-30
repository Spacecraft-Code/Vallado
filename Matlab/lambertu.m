% ------------------------------------------------------------------------------
%
%                           function lambertu
%
%  this function solves the lambert problem for orbit determination and returns
%    the velocity vectors at each of two given position vectors.  the solution
%    uses universal variables for calculation and a bissection technique
%    updating psi.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    r1          - ijk position vector 1          km
%    r2          - ijk position vector 2          km
%    dm          - direction of motion            'l','s'
%    dtsec       - time between r1 and r2         s
%    nrev        - multiple revoluions            0, 1, ...
%
%  outputs       :
%    v1          - ijk velocity vector            km / s
%    v2          - ijk velocity vector            km / s
%    error       - error flag                     'ok', ...
%
%  locals        :
%    vara        - variable of the iteration,
%                  not the semi-axis
%    y           - area between position vectors
%    upper       - upper bound for z
%    lower       - lower bound for z
%    cosdeltanu  - cosine of true anomaly change  rad
%    f           - f expression
%    g           - g expression
%    gdot        - g dot expression
%    x        - old universal variable x
%    xcubed   - x cubed
%    zold        - old value of z
%    znew        - new value of z
%    c2       - c2(z) function
%    c3       - c3(z) function
%    timenew     - new time                       s
%    small       - tolerance for roundoff errors
%    i, j        - index
%
%  coupling      :
%    mag         - magnitude of a vector
%    dot         - dot product of two vectors
%    findc2c3    - find c2 and c3 functions
%
%  references    :
%    vallado       2001, 459-464, alg 55, ex 7-5
%
% [vo,v,errorl] = lambertu ( ro,r, dm, nrev, dtsec );
% ------------------------------------------------------------------------------

%function [vo,v,errorl] = lambertu ( ro,r, dm, nrev, dtsec );
function [vo,v,errorl] = lambertu ( ro,r, dm, nrev, dtsec, outfile )

% -------------------------  implementation   -------------------------
        constmath;
        constastro;
small = 0.00001; % can affect cases where znew is multiples of 2pi^2
        numiter= 40;
        errorl  = '      ok';
        for i= 1 : 3
            vo(i)= 0.0;
            v(i) = 0.0;
        end

 %fprintf(1,'%11.7f %11.7f nrev %3i %1c \n',cosdeltanu*rad, vara , nrev, dm);

        % --------- set up initial bounds for the bissection ----------
        if ( nrev == 0 )  
            lower = -4.0*pi*pi;  % allow hyperbolic and parabolic solutions
            upper =  4.0*pi*pi;  % could be negative infinity for all cases
        else
            % set absolute limits for multi-rev case
            lower = 4.0*nrev^2*pi*pi;
            upper = 4.0*(nrev + 1.0)^2*pi*pi;       
            % adjust based on long or short way
            if dm == 'l'
              upper = lower + (upper - lower)*0.6;             
            else
              lower = lower + (upper - lower)*0.4;             
            end          
        end

        % ---------------  form initial guesses   ---------------------
        dtdpsi = 0.0;
        x   = 0.0;
        psinew = 0.0;
        if (nrev == 0)
            % use log to get initial guess
            psiold = (log(dtsec/806.811874)-1.2357)/0.118;
            if psiold > upper
                psiold = upper - pi;
            end
        else
            if (dm == 's')
                psiold = lower + (upper - lower)*0.75;
            else
                psiold = lower + (upper - lower)*0.25;
            end
        end
        
        [c2,c3] = findc2c3( psiold );

        magro = mag(ro);
        magr  = mag(r);

        % find all the nrev cases for tbi etc
        %[tbi] = findlambertmins( ro, r1, dm );

        cosdeltanu= dot(ro,r)/(magro*magr);  % nrev +
        %cosdeltanu*180/pi
        if ( dm == 'l' )  
            vara = -sqrt( magro*magr*(1.0+cosdeltanu) );
        else
            vara =  sqrt( magro*magr*(1.0+cosdeltanu) );
        end

%        chord = sqrt( magro^2 + magr^2 - 2.0*magro*magr*cosdeltanu );
%            nrev = 1;
%        chord = sqrt( magro^2 + magr^2 - 2.0*magro*magr*cosdeltanu );
%        s     = ( magro + magr + chord )*0.5;
%        betam = 2.0* asin( sqrt((s-chord)/chord) );  % comes out imaginary?? just for hyperbolic??
%        tmin  = ((2.0*nrev+1.0)*pi-betam + sin(betam))/sqrt(mu);

        % -------  determine if  the orbit is possible at all ---------
        if ( abs( vara ) > small )  
            loops  = 0;
            ynegktr= 1;  % y neg ktr
            dtnew = -10.0;
            while ((abs(dtnew-dtsec) >= small) && (loops < numiter) && (ynegktr <= 10))
                % fprintf(1,'%3i  dtnew-dtsec %11.7f yneg %3i \n',loops,dtnew-dtsec,ynegktr );
                if ( abs(c2) > small )
                    y= magro + magr - ( vara*(1.0 - psiold*c3)/sqrt(c2) );
                else
                    y= magro + magr;
                end
                % ----------- check for negative values of y ----------
                if ( (vara > 0.0) && (y < 0.0) )  % ( vara > 0.0 ) &
                    ynegktr= 1;
                    while (( y < 0.0 ) && ( ynegktr < 10 ))
                        psinew = 0.8*(1.0 / c3)*( 1.0 - (magro + magr)*sqrt(c2)/vara  );  
                        % -------- find c2 and c3 functions -----------
                        [c2,c3] = findc2c3( psinew );
                        psiold = psinew;
                        lower  = psiold;
                        if ( abs(c2) > small )
                            y= magro + magr - ( vara*(1.0-psiold*c3) / sqrt(c2) );
                        else
                            y= magro + magr;
                        end
       %             fprintf(outfile,'yneg %3i  y %11.7f lower %11.7f c2 %11.7f psinew %11.7f yneg %3i \n',loops,y,lower,c2,psinew,ynegktr );
                        ynegktr = ynegktr + 1;
                    end % while
                end  % if  y neg

                if ( ynegktr < 10 )  
                    if ( abs(c2) > small )  
                        x= sqrt( y / c2 );
                    else
                        x= 0.0;
                    end
                    xcubed= x^3;
                   
                    dtnew    = (xcubed*c3 + vara*sqrt(y)) / sqrt(mu);
                    % try newton rhapson iteration to update psi
                    if abs(psiold) > 1e-5 
                        c2dot = 0.5/psiold * (1.0 - psiold*c3 - 2.0*c2);
                        c3dot = 0.5/psiold * (c2 - 3.0*c3);
                    else
                        c2dot = -1.0/factorial(4) + 2.0*psiold/factorial(6) - 3.0*psiold^2/factorial(8) + 4.0*psiold^3/factorial(10) - 5.0*psiold^4/factorial(12);
                        c3dot = -1.0/factorial(5) + 2.0*psiold/factorial(7) - 3.0*psiold^2/factorial(9) + 4.0*psiold^3/factorial(11) - 5.0*psiold^4/factorial(13);
                    end    
                    dtdpsi = (xcubed*(c3dot - 3.0*c3*c2dot/(2.0*c2)) + vara/8.0 * (3.0*c3*sqrt(y)/c2 + vara/x)) / sqrt(mu);
                    psinew =  psiold - (dtnew - dtsec)/dtdpsi;  
                    
                    % check if newton guess for psi is outside bounds (too steep a slope)
                    if abs(psinew) > upper || psinew < lower 
                        % --------  readjust upper and lower bounds -------
                        if ( dtnew < dtsec )
                            if psiold > lower
                                lower= psiold; %upper - pi; %psiold;
                            end
                        end
                        if ( dtnew > dtsec )
                            if psiold < upper
                                upper= psiold;
                            end
                        end
                        psinew= (upper+lower) * 0.5;
                    end

                    % ------------- find c2 and c3 functions ----------
                    [c2,c3] = findc2c3( psinew );
%if nrev > 0
%                    fprintf(outfile,'%3i  y %11.7f x %11.7f %11.7f dtnew %11.7f %11.7f %11.7f psinew %11.7f %11.7f \n', ...
%                            loops,y,x,dtsec, dtnew, lower, upper, psinew, (dtnew - dtsec)/dtdpsi );  % c2dot, c3dot
%end
                    psiold = psinew;
                    loops = loops + 1;

                    % --- make sure the first guess isn't too close ---
                    if ( (abs(dtnew - dtsec) < small) && (loops == 1) );
                        dtnew= dtsec - 1.0;
                    end
                end  % if  ynegktr < 10
                
  %            fprintf(1,'%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n',loops,y,x,dtnew,psinew );
%              fprintf(1,'%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,x/sqrt(re),dtnew/tusec,psinew );
%              fprintf(1,'%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n',loops,y/re,x/sqrt(re),dtnew/60.0,psinew );
            end % while loop

            if ( (loops >= numiter) || (ynegktr >= 10) )
                errorl= strcat('gnotconv',num2str(abs(dtnew - dtsec)));
                if ( ynegktr >= 10 )
                    errorl= 'y negati';
                end
            else
                % --- use f and g series to find velocity vectors -----
                f   = 1.0 - y/magro;
                gdot= 1.0 - y/magr;
                %if (nrev > 0)
                %    g   = 1.0/(-2.0*pi*nrev) + 1.0 / (vara*sqrt( y/mu ));  % 1 over g
                %    fdot = (magro*magr + y^2 - magr*y - magro*y - 1.0)*sqrt(mu)/(-2.0*pi*nrev*sqrt(mu) + vara*sqrt(y));
                %else
                    g   = 1.0 / (vara*sqrt( y/mu ));  % 1 over g
                %    fdot = (magro*magr + y^2 - magr*y - magro*y - 1.0)*sqrt(mu)/(-2.0*pi*nrev*sqrt(mu) + vara*sqrt(y));
                %end
                %fdot = sqrt(mu*y)*(-magr-magro+y)/(magro*magr*vara);
                %f*gdot - fdot/g
                for i= 1 : 3
                    vo(i)= ( r(i) - f*ro(i) )*g;
                    v(i) = ( gdot*r(i) - ro(i) )*g;
                end
            end   % if  the answer has converged
        else
            errorl= 'impos180';
        end  % if  var a > 0.0
          
        % fprintf( fid,'psinew %11.5f  %11.5f %11.5f  \n',psinew, dtsec/60.0, x/rad);
        % write out final iteration results
        %fprintf( 1,'%3i  y %11.7f x %11.7f dtnew %11.7f %11.7f %11.7f psinew %11.7f %11.7f \n', ...
        %         loops,y,x,dtnew, lower, upper, psinew, (dtnew - dtsec)/dtdpsi );  % c2dot, c3dot
        %fprintf( 1,'   conv     #s nrev  dm      y         x         dtnew        upper      lower      psinew   slope  \n');
        if (strcmp(errorl, '      ok') ~= 0)
            [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (ro,vo);
            fprintf(outfile,'%10s %3i %3i  %2s %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.8f \n', ...
                     errorl, loops, nrev, dm, dtnew, y, x, vo(1), vo(2), vo(3), v(1), v(2), v(3), lower, upper, psinew, (dtnew - dtsec)/dtdpsi, ecc );  % c2dot, c3dot
            %fprintf(1,'%10s %3i %3i  %2s %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.8f \n', ...
            %         errorl, loops, nrev, dm, dtnew, y, x, vo(1), vo(2), vo(3), v(1), v(2), v(3), lower, upper, psinew, (dtnew - dtsec)/dtdpsi, ecc );  % c2dot, c3dot
            fprintf(1,'C%3i %3i  %2s %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f  %11.7f \n', ...
                   loops, nrev, dm, dtnew, magro, magr, vara, y, x, psinew)
        else
            fprintf(outfile,'%s \n',errorl);
            fprintf(1,'%s \n',errorl);
        end
    
end         
% ------------------------------------------------------------------------------

%C  4   0   s 21000.0000000 9103.9997818 15945.3423626 4059.6664884 27881.0982880 566.4217904 
%  50   1   s 928398.5940708 9103.9997818 9031.3025942 12518.5302475 34615.6582754 4256.2996647 
       
       
        
        
       
       
