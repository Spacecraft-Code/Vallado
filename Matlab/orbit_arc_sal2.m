
% ------------------------------------------------------------------------------
%
%                           function orbit_arc_sal2
%
%
% Computes the length of an arc on a closed orbit
% between true anomaly angles nu1 and nu2 (in radians). 
% Arc length is determined by dividing the arc in small straight segments.
%
%  author        : sal alfano         719-573-2600   24 sep 2012
%

function [arc] = orbit_arc_sal2(a,ecc,nu1,nu2, ang_step)

    % determine eccentric anomalies from true
%                sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu1) ) / ( 1.0 +ecc*cos(nu1) );
%                cose= ( ecc + cos(nu1) ) / ( 1.0  + ecc*cos(nu1) );
%                E1  = atan2( sine,cose );
    E1 = 2*atan(sqrt((1-ecc)/(1+ecc))*tan(nu1*0.5));
%                sine= ( sqrt( 1.0 -ecc*ecc ) * sin(nu2) ) / ( 1.0 +ecc*cos(nu2) );
%                cose= ( ecc + cos(nu2) ) / ( 1.0  + ecc*cos(nu2) );
%                E2  = atan2( sine,cose );
    E2 = 2*atan(sqrt((1-ecc)/(1+ecc))*tan(nu2*0.5));
    E_diff = E2-E1;

    % test for very small orbital arc, if small then ignore
    if abs(E_diff) > 1e-10   % old is -10

        % correct direction by changing E2 to be within half rev of E1
        if E_diff < -pi
            E_diff = E_diff+2*pi;
        end;
        if E_diff > pi
            E_diff = E_diff-2*pi;
        end;
        E2 = E1+E_diff;

        % set angle step size (adjust ang_step as needed)
        % this tolerance should be identical to find_nu_from_orbit_arc_sal.m
%        ang_step = 0.0000001;  % ' old is -8
        num_steps = ceil(abs(E_diff)/ang_step)+2;
        dE = ang_step*sign(E_diff);

        % loop through and sum parts
        arc_unit = 0;

        E = E1-dE*0.5;
        for stp = 1:num_steps
            E = E+dE;
            ecc_cos = ecc*cos(E);
            d_arc = sqrt(1-ecc_cos*ecc_cos);
            arc_unit = arc_unit + d_arc;
            if abs(E-E1) >=  abs(E_diff)
                arc_unit = arc_unit-d_arc * (E-E1-E_diff)/dE;
                break;
            end;
        end;
        arc = arc_unit*a*dE;
    else
        % if E_diff is less than a small amount
        arc = a*E_diff;
    end;  






