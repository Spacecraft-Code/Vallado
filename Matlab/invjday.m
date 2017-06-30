% ------------------------------------------------------------------------------
%
%                           function invjday
%
%  this function finds the year, month, day, hour, minute and second
%    given the julian date. tu can be ut1, tdt, tdb, etc.
%
%  author        : david vallado                  719-573-2600   27 may 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    jd          - julian date                    days from 4713 bc
%
%  outputs       :
%    year        - year                           1900 .. 2100
%    mon         - month                          1 .. 12
%    day         - day                            1 .. 28,29,30,31
%    hr          - hour                           0 .. 23
%    min         - minute                         0 .. 59
%    sec         - second                         0.0 .. 59.999
%
%  locals        :
%    days        - day of year plus fractional
%                  portion of a day               days
%    tu          - julian centuries from 0 h
%                  jan 0, 1900
%    temp        - temporary real values
%    leapyrs     - number of leap years from 1900
%
%  coupling      :
%    days2mdhms  - finds month, day, hour, minute and second given days and year
%
%  references    :
%    vallado       2007, 208, alg 22, ex 3-13
%
% [year,mon,day,hr,min,sec] = invjday ( jd, jdfrac );
% -----------------------------------------------------------------------------

function [year,mon,day,hr,min,sec] = invjday ( jd, jdfrac );

     % check jdfrac for multiple days
     if abs(jdfrac) >= 1.0 
         jd = jd + floor(jdfrac);
         jdfrac = jdfrac - floor(jdfrac);
     end

     % check for fraction of a day included in the jd
	 dt = jd - floor(jd);
	 if (abs(dt - 0.5) > 0.00000001)
	 	jdfrac = jdfrac + dt - 0.5;
     end

     % ----------------- find year and days of the year ---------------
     temp   = jd + jdfrac - 2415019.5;
     tu     = temp / 365.25;
     year   = 1900 + floor( tu );
     leapyrs= floor( ( year-1901 )*0.25 );
     days   = floor(temp - ((year-1900)*365.0 + leapyrs ));

     % ------------ check for case of beginning of a year -------------
     if days < 1.0
         year   = year - 1;
         leapyrs= floor( ( year-1901 )*0.25 );
         days   = floor(temp - ((year-1900)*365.0 + leapyrs ));
     end

     % ------------------- find remaining data  -----------------------
     % now add the daily time in to preserve accuracy
     [mon,day,hr,min,sec] = days2mdh( year, days + jdfrac );

     
