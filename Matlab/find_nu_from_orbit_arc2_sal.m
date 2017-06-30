% ------------------------------------------------------------------------------
%
%                           function find_nu_from_orbit_arc_sal
%
%
% Finds true anomaly associated with arc length on a closed orbit
% starts at nu1 and walks to nu2 (in radians) to tolerance dE. 
% Arc length is determined by dividing the arc in small straight segments.
%
%  author        : sal alfano         719-573-2600   13 aug 2011
%
% ------------------------------------------------------------------------------

function [nu2] = find_nu_from_orbit_arc2_sal(a,ecc,nu1,arc, ang_step)

% test for very small orbital arc, if small then ignore   
    if abs(arc) < 1e-10
      nu2 = nu1;
      return
    end;
    
% determine eccentric anomaly from true
    E1 = 2*atan(sqrt((1-ecc)/(1+ecc))*tan(nu1*0.5));
    
% set angle step size (adjust ang_step as needed)
% this tolerance should be identical to orbit_arc_sal.m
%    ang_step = 0.0000001; % old -8 this tolerance should be identical to 
    num_steps = pi/ang_step;
    dE = ang_step*sign(arc);
    
% loop through and sum parts until just past arc, then go back a bit
    arc_unit = 0;
    arc_end = abs(arc/(a*dE));
    E = E1-dE*0.5;
    for stp = 1:num_steps
      E = E+dE;
      ecc_cos = ecc*cos(E);
      d_arc = sqrt(1-ecc_cos*ecc_cos);
      arc_unit = arc_unit+d_arc;
      if arc_unit >=  arc_end
        E = E+dE*(arc_end-arc_unit)/d_arc;
        break;
      end;
    end;
 
% find nu2
    nu2 = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E*0.5));

    

