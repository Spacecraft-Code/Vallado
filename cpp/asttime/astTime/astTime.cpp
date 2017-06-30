/*     ----------------------------------------------------------------
*
*                                asttime.cpp
*
*   This file contains fundamental Astrodynamic procedures and functions
*   relating to the time functions. These routines are discussed in Ch 3
*   and Ch 5.
*
*                          Companion code for
*             Fundamentals of Astrodynamics and Applications
*                                  2013
*                            by David Vallado
*
*       (w) 719-573-2600, email dvallado@agi.com, davallado@gmail.com
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*               4 may 09  david vallado
*                           misc updates
*              13 feb 08  david vallado
*                           add getmon
*              15 mar 07  david vallado
*                           3rd edition baseline
*              21 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */

// stdafx.h must be the first in cpp files
//#include "stdafx.h"

#include "astTime.h"

namespace astTime {

	int getmon
		(
		char instr[3]
		)
	{
		int ktr;
		typedef char str3[4];
		str3 monstr[13];

		// ------------------------  implementation   --------------------------
		strcpy_s(monstr[1], "Jan");
		strcpy_s(monstr[2], "Feb");
		strcpy_s(monstr[3], "Mar");
		strcpy_s(monstr[4], "Apr");
		strcpy_s(monstr[5], "May");
		strcpy_s(monstr[6], "Jun");
		strcpy_s(monstr[7], "Jul");
		strcpy_s(monstr[8], "Aug");
		strcpy_s(monstr[9], "Sep");
		strcpy_s(monstr[10], "Oct");
		strcpy_s(monstr[11], "Nov");
		strcpy_s(monstr[12], "Dec");

		ktr = 1;
		while ((strcmp(instr, monstr[ktr]) != 0) && (ktr <= 12))
			ktr = ktr + 1;

		return (ktr);
	} // getmon


	/*       ----------------------------------------------------------------      */
	int getday
		(
		char instr[3]
		)
	{
		int ktr;
		typedef char str3[4];
		str3 monstr[8];

		// ------------------------  implementation   --------------------------
		strcpy_s(monstr[1], "Sun");
		strcpy_s(monstr[2], "Mon");
		strcpy_s(monstr[3], "Tue");
		strcpy_s(monstr[4], "Wed");
		strcpy_s(monstr[5], "Thr");
		strcpy_s(monstr[6], "Fri");
		strcpy_s(monstr[7], "Sat");

		ktr = 1;
		while ((strcmp(instr, monstr[ktr]) != 0) && (ktr <= 7))
			ktr = ktr + 1;

		return (ktr);
	} // getday


	/* -----------------------------------------------------------------------------
	*
	*                           function dayofweek
	*
	*  this function finds the day of the week. integers are used for the days,
	*    1 = 'Sun', 2 = 'Mon', ... 7 = 'Sat'.
	*
	*  Author        : David Vallado                  719-573-2600    1 Mar 2001
	*
	*  Inputs          Description                    Range / Units
	*    JD          - Julian date of interest        days from 4713 BC
	*
	*  OutPuts       :
	*    dayofweek   - answer                         1 to 7
	*
	*  Locals        :
	*    None.
	*
	*  Coupling      :
	*    None.
	*
	*  References    :
	* --------------------------------------------------------------------------- */

	int dayofweek
		(
		double jd
		)
	{
		int temp;
		// ----- Be sure jd is at 0.0 h on the day of interest -----
		jd = floor(jd + 0.5);

		temp = int(floor(jd - 7 * floor((jd + 1) / 7) + 2));
		return temp;
	}  // function dayofweek


	/* -----------------------------------------------------------------------------
	*
	*                           procedure jday
	*
	*  this procedure finds the julian date given the year, month, day, and time.
	*    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
	*
	*  algorithm     : calculate the answer in one step for efficiency
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    year        - year                           1900 .. 2100
	*    mon         - month                          1 .. 12
	*    day         - day                            1 .. 28,29,30,31
	*    hr          - universal time hour            0 .. 23
	*    min         - universal time min             0 .. 59
	*    sec         - universal time sec             0.0 .. 59.999
	*
	*  outputs       :
	*    jd          - julian date (days only)           days from 4713 bc
	*    jdFrac      - julian date (fraction of a day)   days from 0 hr of the day
	*
	*  locals        :
	*    none.
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 183, alg 14, ex 3-4
	*
	* --------------------------------------------------------------------------- */

	void jday
		(
		int year, int mon, int day, int hr, int minute, double sec,
		double& jd, double& jdFrac
		)
	{
		jd = 367.0 * year -
			floor((7 * (year + floor((mon + 9) / 12.0))) * 0.25) +
			floor(275 * mon / 9.0) +
			day + 1721013.5;  // use - 678987.0 to go to mjd directly
		jdFrac = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;

		// check that the day and fractional day are correct
		if (fabs(jdFrac) >= 1.0)
		{
			double dtt = floor(jdFrac);
			jd = jd + dtt;
			jdFrac = jdFrac - dtt;
		}

		// - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure jdayall
	*
	*  this procedure finds the julian date given the year, month, day, and time.
	*    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
	*
	*  algorithm     : calculate the answer in one step for efficiency
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    year        - year                           all, 1900 .. 2100
	*    mon         - month                          1 .. 12
	*    day         - day                            1 .. 28,29,30,31
	*    hr          - universal time hour            0 .. 23
	*    min         - universal time min             0 .. 59
	*    sec         - universal time sec             0.0 .. 59.999
	*    whichtype   - julian or gregorian calender   'j' or 'g'
	*
	*  outputs       :
	*    jd          - julian date                    days from 4713 bc
	*
	*  locals        :
	*    b           - var to aid gregorian dates
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 183, alg 14, ex 3-4
	* --------------------------------------------------------------------------- */

	void    jdayall
		(
		int year, int mon, int day, int hr, int minute, double sec,
		char whichtype, double& jd
		)
	{
		double b;

		if (mon <= 2)
		{
			year = year - 1;
			mon = mon + 12;
		}
		/* --------- use for julian calender, every 4 years --------- */
		if (whichtype == 'j')
			b = 0.0;
		else
		{
			/* ---------------------- use for gregorian ----------------- */
			b = 2 - floor(year * 0.01) + floor(floor(year * 0.01) * 0.25);
			jd = floor(365.25 * (year + 4716)) +
				floor(30.6001 * (mon + 1)) +
				day + b - 1524.5 +
				((sec / 60.0 + minute) / 60.0 + hr) / 24.0;  // ut in days
		}
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure days2mdhms
	*
	*  this procedure converts the day of the year, days, to the equivalent month
	*    day, hour, minute and second.
	*
	*  algorithm     : set up array for the number of days per month
	*                  find leap year - use 1900 because 2000 is a leap year
	*                  loop through a temp value while the value is < the days
	*                  perform int conversions to the correct day and month
	*                  convert remainder into h m s using type conversions
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    year        - year                           1900 .. 2100
	*    days        - julian day of the year         1.0  .. 366.0
	*
	*  outputs       :
	*    mon         - month                          1 .. 12
	*    day         - day                            1 .. 28,29,30,31
	*    hr          - hour                           0 .. 23
	*    min         - minute                         0 .. 59
	*    sec         - second                         0.0 .. 59.999
	*
	*  locals        :
	*    dayofyr     - day of year
	*    temp        - temporary extended values
	*    inttemp     - temporary int value
	*    i           - index
	*    lmonth[12]  - int array containing the number of days per month
	*
	*  coupling      :
	*    none.
	* --------------------------------------------------------------------------- */

	void    days2mdhms
		(
		int year, double days,
		int& mon, int& day, int& hr, int& minute, double& sec
		)
	{
		int i, inttemp, dayofyr;
		double    temp;
		// start indicies at 1
		int lmonth[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

		dayofyr = (int)floor(days);
		/* ----------------- find month and day of month ---------------- */
		if ((year % 4) == 0)
			lmonth[2] = 29;

		i = 1;
		inttemp = 0;
		while ((dayofyr > inttemp + lmonth[i]) && (i < 12))
		{
			inttemp = inttemp + lmonth[i];
			i++;
		}
		mon = i;
		day = dayofyr - inttemp;

		/* ----------------- find hours minutes and seconds ------------- */
		temp = (days - dayofyr) * 24.0;
		hr = floor(temp);
		temp = (temp - hr) * 60.0;
		minute = floor(temp);
		sec = (temp - minute) * 60.0;
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure invjday
	*
	*  this procedure finds the year, month, day, hour, minute and second
	*  given the julian date. tu can be ut1, tdt, tdb, etc.
	*
	*  algorithm     : set up starting values
	*                  find leap year - use 1900 because 2000 is a leap year
	*                  find the elapsed days through the year in a loop
	*                  call routine to find each individual value
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    jd          - julian date (days only)           days from 4713 bc
	*    jdFrac      - julian date (fraction of a day)   days from 0 hr of the day
	*
	*  outputs       :
	*    year        - year                           1900 .. 2100
	*    mon         - month                          1 .. 12
	*    day         - day                            1 .. 28,29,30,31
	*    hr          - hour                           0 .. 23
	*    min         - minute                         0 .. 59
	*    sec         - second                         0.0 .. 59.999
	*
	*  locals        :
	*    days        - day of year plus fractional
	*                  portion of a day               days
	*    tu          - julian centuries from 0 h
	*                  jan 0, 1900
	*    temp        - temporary double values
	*    leapyrs     - number of leap years from 1900
	*
	*  coupling      :
	*    days2mdhms  - finds month, day, hour, minute and second given days and year
	*
	*  references    :
	*    vallado       2013, 203, alg 22, ex 3-13
	* --------------------------------------------------------------------------- */

	void    invjday
		(
		double jd, double jdFrac,
		int& year, int& mon, int& day,
		int& hr, int& minute, double& sec
		)
	{
		int leapyrs;
		double    days, tu, temp;

		// check jdfrac for multiple days
		if (fabs(jdFrac) >= 1.000000)
		{
			double dtt = floor(jdFrac);
			jd = jd + dtt;
			jdFrac = jdFrac - dtt;
		}

		// check for fraction of a day included in the jd
		double dt = jd - floor(jd);
		if (fabs(dt) - 0.5 > 0.00000001)
			jdFrac = jdFrac + dt - 0.5;

		/* --------------- find year and days of the year --------------- */
		temp = jd - 2415019.5;
		tu = temp / 365.25;
		year = 1900 + (int)floor(tu);
		leapyrs = (int)floor((year - 1901) * 0.25);

		days = floor(temp - ((year - 1900) * 365.0 + leapyrs));

		/* ------------ check for case of beginning of a year ----------- */
		if (days < 1.0)
		{
			year = year - 1;
			leapyrs = (int)floor((year - 1901) * 0.25);
			days = floor(temp - ((year - 1900) * 365.0 + leapyrs));
		}

		/* ----------------- find remaing data  ------------------------- */
		// now add the daily time in to preserve accuracy
		days2mdhms(year, days + jdFrac, mon, day, hr, minute, sec);
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure finddays
	*
	*  this procedure finds the fractional days through a year given the year,
	*    month, day, hour, minute and second.
	*
	*  algorithm     : set up array for the number of days per month
	*                  find leap year - use 1900 because 2000 is a leap year
	*                  check for a leap year
	*                  loop to find the elapsed days in the year
	*
	*  author        : david vallado                  719-573-2600
	*
	*  inputs          description                    range / units
	*    year        - year                           1900 .. 2100
	*    mon         - month                          1 .. 12
	*    day         - day                            1 .. 28,29,30,31
	*    hr          - hour                           0 .. 23
	*    min         - minute                         0 .. 59
	*    sec         - second                         0.0 .. 59.999
	*
	*  outputs       :
	*    days        - day of year plus fraction of a
	*                    day                          days
	*
	*  locals        :
	*    lmonth      - length of months of year
	*    i           - index
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 201, ex 3-12
	* --------------------------------------------------------------------------- */

	void    finddays
		(
		int year, int month, int day, int hr, int minute,
		double sec, double& days
		)
	{
		// shift index values to be 1-12
		int lmonth[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
		int i;

		if (((year - 1900) % 4) == 0)
			lmonth[2] = 29;

		i = 1;
		days = 0.0;
		while ((i < month) && (i < 12))
		{
			days = days + lmonth[i];
			i = i + 1;
		}
		days = days + day + hr / 24.0 + minute / 1440.0 + sec / 86400.0;
	}

	/* -----------------------------------------------------------------------------
	*
	*                           function gstime
	*
	*  this function finds the greenwich sidereal time (iau-82).
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    jdut1       - julian date in ut1             days from 4713 bc
	*
	*  outputs       :
	*    gstime      - greenwich sidereal time        0 to 2pi rad
	*
	*  locals        :
	*    temp        - temporary variable for doubles   rad
	*    tut1        - julian centuries from the
	*                  jan 1, 2000 12 h epoch (ut1)
	*
	*  coupling      :
	*    none
	*
	*  references    :
	*    vallado       2013, 187, eq 3-45
	* --------------------------------------------------------------------------- */

	double  gstime
		(
		double jdut1
		)
	{
		const double twopi = 2.0 * pi;
		const double deg2rad = pi / 180.0;
		double       temp, tut1;

		tut1 = (jdut1 - 2451545.0) / 36525.0;
		temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
			(876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841;  // sec
		temp = fmod(temp * deg2rad / 240.0, twopi); //360/86400 = 1/240, to deg, to rad

		// ------------------------ check quadrants ---------------------
		if (temp < 0.0)
			temp += twopi;

		return temp;
	}

	/* -----------------------------------------------------------------------------
	*                           procedure lstime
	*
	*  this procedure finds the local sidereal time at a given location.
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    lon         - site longitude (west -)        -2pi to 2pi rad
	*    jdut1       - julian date in ut1             days from 4713 bc
	*
	*  outputs       :
	*    lst         - local sidereal time            0.0 to 2pi rad
	*    gst         - greenwich sidereal time        0.0 to 2pi rad
	*
	*  locals        :
	*    none.
	*
	*  coupling      :
	*    gstime        finds the greenwich sidereal time
	*
	*  references    :
	*    vallado       2013, 188, alg 15, ex 3-5
	* --------------------------------------------------------------------------- */

	void    lstime
		(
		double lon, double jdut1, double& lst, double& gst
		)
	{
		const double twopi = 2.0 * pi;

		gst = gstime(jdut1);
		lst = lon + gst;

		/* ------------------------ check quadrants --------------------- */
		lst = fmod(lst, twopi);
		if (lst < 0.0)
			lst = lst + twopi;
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure hms_sec
	*
	*  this procedure converts hours, minutes and seconds to seconds from the
	*  beginning of the day.
	*
	*  author        : david vallado                  719-573-2600   25 jun 2002
	*
	*  revisions
	*                -
	*
	*  inputs          description                    range / units
	*    utsec       - seconds                        0.0 .. 86400.0
	*
	*  outputs       :
	*    hr          - hours                          0 .. 24
	*    min         - minutes                        0 .. 59
	*    sec         - seconds                        0.0 .. 59.99
	*
	*  locals        :
	*    temp        - temporary variable
	*
	*  coupling      :
	*    none.
	* --------------------------------------------------------------------------- */

	void    hms_sec
		(
		int& hr, int& min, double& sec, edirection direct, double& utsec
		)
	{
		double temp;

		// ------------------------  implementation   ------------------
		if (direct == eTo)
			utsec = hr * 3600.0 + min * 60.0 + sec;
		else
		{
			temp = utsec / 3600.0;
			hr = (int)floor(temp);
			min = (int)floor((temp - hr)* 60.0);
			sec = (temp - hr - min / 60.0) * 3600.0;
		}
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure hms_ut
	*
	*  this procedure converts hours, minutes and seconds into universal time.
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    hr          - hours                          0 .. 24
	*    min         - minutes                        0 .. 59
	*    sec         - seconds                        0.0 .. 59.99
	*    direction   - which set of vars to output    from  too
	*
	*  outputs       :
	*    ut          - universal time                 hrmin.sec
	*
	*  locals        :
	*    none.
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 199, alg 21, ex 3-10
	* --------------------------------------------------------------------------- */

	void    hms_ut
		(
		int& hr, int& min, double& sec, edirection direct, double& ut
		)
	{
		// ------------------------  implementation   ------------------
		if (direct == eTo)
			ut = hr * 100.0 + min + sec * 0.01;
		else
		{
			hr = (int)floor(ut * 0.01);
			min = (int)floor(ut - hr*100.0);
			sec = (ut - hr*100.0 - min) * 100.0;
		}
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure hms_rad
	*
	*  this procedure converts hours, minutes and seconds into radians.
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    hr          - hours                          0 .. 24
	*    min         - minutes                        0 .. 59
	*    sec         - seconds                        0.0 .. 59.99
	*    direction   - which set of vars to output    from  too
	*
	*  outputs       :
	*    hms         - result                         rad
	*
	*  locals        :
	*    temp        - conversion from hours to rad   0.261799
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 198, alg 19 alg 20, ex 3-9
	* --------------------------------------------------------------------------- */

	void    hms_rad
		(
		int& hr, int& min, double& sec, edirection direct, double& hms
		)
	{
		const double rad2deg = 57.29577951308230;
		double temp;

		// ------------------------  implementation   ------------------
		temp = 15.0 / rad2deg;
		if (direct == eTo)
			hms = hr + min / 60.0 + sec / 3600.0;
		else
		{
			temp = hms / temp;
			hr = int(temp);
			min = int((temp - hr)* 60.0);
			sec = (temp - hr - min / 60.0) * 3600.0;
		}
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure dms_rad
	*
	*  this procedure converts degrees, minutes and seconds into radians.
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    deg         - degrees                        0 .. 360
	*    min         - minutes                        0 .. 59
	*    sec         - seconds                        0.0 .. 59.99
	*    direction   - which set of vars to output    from  too
	*
	*  outputs       :
	*    dms         - result                         rad
	*
	*  locals        :
	*    temp        - temporary variable
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 197, alg 17 alg 18, ex 3-8
	* --------------------------------------------------------------------------- */

	void    dms_rad
		(
		int& deg, int& min, double& sec, edirection direct, double& dms
		)
	{
		const double rad2deg = 57.29577951308230;
		double temp;

		// ------------------------  implementation   ------------------
		if (direct == eTo)
			dms = (deg + min / 60.0 + sec / 3600.0) / rad2deg;
		else
		{
			temp = dms * rad2deg;
			deg = (int)floor(temp);
			min = (int)floor((temp - deg)* 60.0);
			sec = (temp - deg - min / 60.0) * 3600.0;
		}
	}

	/* -----------------------------------------------------------------------------
	*
	*                           function jd2sse
	*
	*  this function finds the seconds since epoch (1 Jan 2000) given the julian date
	*
	*  author        : david vallado                  719-573-2600   12 dec 2002
	*
	*  revisions
	*                -
	*
	*  inputs          description                    range / units
	*    jd          - julian date                    days from 4713 bc
	*
	*  outputs       :
	*    sse         - seconds since epoch 1 jan 2000
	*
	*  locals        :
	*    none.
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*
	* sse = jd2sse( jd );
	* ----------------------------------------------------------------------------- */

	double jd2sse
		(
		double jd
		)
	{
		double temp;
		// ------------------------  implementation   ------------------
		temp = (jd - 2451544.5) * 86400.0;
		return temp;
	}  // function jd2sse


	/* -----------------------------------------------------------------------------
	*
	*                           procedure convtime
	*
	*  this procedure finds the time parameters and julian century values for inputs
	*    of utc or ut1. numerous outputs are found as shown in the local variables.
	*    because calucations are in utc, you must include timezone if ( you enter a
	*    local time, otherwise it should be zero.
	*
	*  author        : david vallado                  719-573-2600    4 jun 2002
	*
	*  revisions
	*    vallado     - fix documentation for dut1                     8 oct 2002
	*
	*  inputs          description                    range / units
	*    year        - year                           1900 .. 2100
	*    mon         - month                          1 .. 12
	*    day         - day                            1 .. 28,29,30,31
	*    hr          - universal time hour            0 .. 23
	*    min         - universal time min             0 .. 59
	*    sec         - universal time sec (utc)            0.0  .. 59.999
	*    timezone    - offset to utc from local site  0 .. 23 hr
	*    dut1        - delta of ut1 - utc             sec
	*    dat         - delta of utc - tai             sec
	*
	*  outputs       :
	*    ut1         - universal time                 sec
	*    tut1        - julian centuries of ut1
	*    jdut1       - julian date (days only)           days from 4713 bc
	*    jdut1Frac   - julian date (fraction of a day)   days from 0 hr of the day
	*    utc         - coordinated universal time     sec
	*    tai         - atomic time                    sec
	*    tdt         - terrestrial dynamical time     sec
	*    ttdt        - julian centuries of tdt
	*    jdtt        - julian date (days only)           days from 4713 bc
	*    jdttFrac    - julian date (fraction of a day)   days from 0 hr of the day
	*    tcg         - geocentric coordinate time     sec
	*    tdb         - terrestrial barycentric time   sec
	*    ttdb        - julian centuries of tdb
	*    jdtdb       - julian date (days only)           days from 4713 bc
	*    jdtdbFrac   - julian date (fraction of a day)   days from 0 hr of the day
	*
	*  locals        :
	*    hrtemp      - temporary hours                hr
	*    mintemp     - temporary minutes              min
	*    sectemp     - temporary seconds              sec
	*    localhr     - difference to local time       hr
	*    jd          - julian date of request         days from 4713 bc
	*    me          - mean anomaly of the earth      rad
	*
	*  coupling      :
	*    hms_2_sec   - conversion between hr-min-sec .and. seconds
	*    jday        - find the julian date
	*
	*  references    :
	*    vallado       2013, 195, alg 16, ex 3-7
	* --------------------------------------------------------------------------- */

	void    convtime
		(
		int year, int mon, int day, int hr, int min, double sec, int timezone,
		double dut1, int dat,
		double& ut1, double& tut1, double& jdut1, double& jdut1Frac, double& utc, double& tai,
		double& tt, double& ttt, double& jdtt, double& jdttFrac, double& tcg, double& tdb,
		double& ttdb, double& jdtdb, double& jdtdbFrac, double& tcb
		)
	{
		double deg2rad, jd, jdFrac, sectemp, me;
		int localhr, hrtemp, mintemp;

		deg2rad = pi / 180.0;

		// ------------------------  implementation   ------------------
		jday(year, mon, day, 0, 0, 0.0, jd, jdFrac);
		//     mjd  = jd - 2400000.5;
		//     mfme = hr*60.0 + min + sec/60.0;

		// ------------------ start if ( ut1 is known ------------------
		localhr = timezone + hr;

		hms_sec(localhr, min, sec, eTo, utc);
		ut1 = utc + dut1;
		hms_sec(hrtemp, mintemp, sectemp, eFrom, ut1);
		jday(year, mon, day, hrtemp, mintemp, sectemp, jdut1, jdut1Frac);
		tut1 = (jdut1 - 2451545.0) / 36525.0;

		tai = utc + dat;

		tt = tai + 32.184;   // sec
		hms_sec(hrtemp, mintemp, sectemp, eFrom, tt);
		jday(year, mon, day, hrtemp, mintemp, sectemp, jdtt, jdttFrac);
		ttt = (jdtt - 2451545.0) / 36525.0;

		tcg = tt + 6.969290134e-10*(jdut1 - 2443144.5)*86400.0; // sec

		me = 357.5277233 + 35999.05034 *ttt;
		me = fmod(me, 360.0);
		me = me * deg2rad;
		tdb = tt + 0.001657  * sin(me) + 0.00001385 *sin(2.0 *me);
		hms_sec(hrtemp, mintemp, sectemp, eFrom, tdb);
		jday(year, mon, day, hrtemp, mintemp, sectemp, jdtdb, jdtdbFrac);
		ttdb = (jdtdb - 2451545.0) / 36525.0;

		tcb = tdb + 1.55051976772e-8*(jdtt - 2443144.5)*86400.0; // sec

		//     fprintf( 'time %14.6f%14.6f%9.6f%18.5f\n',utc,ut1,ttt,jdut1 );
	}

}  // namespace
