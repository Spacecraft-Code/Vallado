% ------------------------------------------------------------------------------
%
%                           procedure iandnode 
%
%  this procedure calculates the delta v's for a change in inclination and
%    longitude of ascending node.
%
%  author        : david vallado                  719-573-2600    1 mar 2001
%
%  inputs          description                    range / units
%    vinit       - initial velocity vector        er/tu
%    iinit       - initial inclination            rad
%    fpa         - flight path angle              rad
%    deltaomega  - change in node                 rad
%    deltai      - change in inclination          rad
%    rfinal      - final position magnitude       er
%
%  outputs       :
%    ifinal      - final inclination              rad
%    deltav      - change in velocity             er/tu
%
%  locals        :
%    arglat      - argument of latitude           rad
%    arglat1     - final argument of latitude     rad
%    theta       -
%
%  coupling      :
%    acos      - arc cosine function
%
%  references    :
%    vallado       2007, 350, alg 41, ex 6-6
%function [deltav] = iandnode(iinit,deltaomega,ifinal,vinit,fpa);
% ----------------------------------------------------------------------------- }

function [deltav] = iandnode(iinit,deltaomega,ifinal,vinit,fpa);
     rad  = 57.29577951308230;
     deltai= iinit - ifinal;
     theta = acos( cos(iinit) * cos(ifinal) + ...
                      sin(iinit) * sin(ifinal) * cos(deltaomega) );
     deltav= 2.0 * vinit * cos(fpa) * sin(0.5 * theta);

     arglat = acos( (sin(ifinal) * cos(deltaomega) - ...
                       cos(theta) * sin(iinit)) / (sin(theta) * cos(iinit)) );
     arglat1= acos( (cos(iinit) * sin(ifinal) - ...
                       sin(iinit) * cos(ifinal) * cos(deltaomega)) / sin(theta) );

     fprintf(1,' theta   %11.7f  \n',theta*rad );
     fprintf(1,' arglat   %11.7f  %11.7f  \n',arglat*rad, arglat1*rad );

