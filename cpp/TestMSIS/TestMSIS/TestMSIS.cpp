// TestMSIS.cpp : main project file.
/*
 * test driver for msis00 model
 * dav 12 nov 2004 fix ap array set (to 8) and include sw tests for each
 * dav 4 oct 2004
 *
 * The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
 * Doug Drob. They also wrote a NRLMSISE-00 distribution package in
 * FORTRAN which is available at
 * http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm

 
 remove .h from include do not do it for other libraries,
 remove clrscr(); everywhere and replace by system("cls");
 finally after writing all #include stuff add this-
 using namespace std; it allows you to do things like cin>> cout<< etc.
 
 */

/* ------------------------------ INCLUDES --------------------------- */

//#include "stdafx.h"

#include <stdio.h>

#include "astTime.h"
#include "astMath.h"
#include "EopSpw.h"
#include "MSIS_Vers.h"



typedef struct msise_inputrec
{
	int year;      /* year, currently ignored */
	int doy;       /* day of year */
	double sec;    /* seconds in day (UT) */
	double alt;    /* altitude in kilometes */
	double glat;  /* geodetic latitude */
	double glong; /* geodetic longitude */
	double stl;    /* local apparent solar time (hours), see note below */
	double f107a;  /* 81 day average of F10.7 flux (centered on doy) */
	double f107;   /* daily F10.7 flux for previous day */
} msisinput;

        int sv[26]  = {0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
        int svt[26] = {0,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

        FILE *outfile;


// local functions
// test the routine against the nrl test cases
void main();
//void int_tmain();

void testjb06();

// interface for plugin to hpop
void mainplugin( double[], double[], double );



//using namespace System;

//void int_tmain(array<System::String ^> ^args)
void main()
{
 	msisinput input[17];
	int i, mass;
	msistype msis00r;
	msistype msis90r;
        msistype msis86r;
        lpolytype lpoly;
        fittype fit;
        lsqvtype lsqv;
        double d[10],t[3], aph[8], apar[8], apar9[8];

	/* input values */
  	for (i=0; i<8; i++)
          {
            aph[i]   = 100;
            apar[i]  =   4;
            apar9[i] =  40;
          }

        mass = 48;

	for (i=0; i<17; i++)
        {
	    input[i].doy=172;
   	    input[i].year=0; /* without effect */
  	    input[i].sec=29000;
	    input[i].alt=400;
	    input[i].glat=60;
	    input[i].glong=-70;
	    input[i].stl=16;
	    input[i].f107a=150;
	    input[i].f107=150;
	}

        // setup differences for each test
	input[1].doy   = 81;
	input[2].sec   = 75000;
	input[2].alt   = 1000;
	input[3].alt   = 100;
	input[4].glat  = 0;
	input[5].glong = 0;
	input[6].stl   = 4;
	input[7].f107a = 70;
	input[8].f107  = 180;
	input[10].alt  = 0;
	input[11].alt  = 10;
	input[12].alt  = 30;
	input[13].alt  = 50;
	input[14].alt  = 70;
	input[16].alt  = 100;

	printf("now do the msis-86 test cases --------------------------------------\n");
        MSIS_Vers::msis86init(msis86r);

	/* evaluate std test cases */
  	for (i = 0; i < 17; i++)
          {
            printf("------case number %3i\n",i);
            if (i == 2)  // it's different for the msis00 test case
                input[i].alt = 400.0;
            if (i == 9)
				MSIS_Vers::gts5(msis86r, lpoly, lsqv,
                     input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                     input[i].stl, input[i].f107a, input[i].f107, apar9, mass, d, t);
              else
              {
                if (i >= 15)
                  {
// test switch settings
					MSIS_Vers::tselec(msis86r.csw, svt);
					MSIS_Vers::gts5(msis86r, lpoly, lsqv,
                         input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                         input[i].stl, input[i].f107a, input[i].f107, aph, mass, d, t);
                   }
                  else
					  MSIS_Vers::gts5(msis86r, lpoly, lsqv,
                         input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                         input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
              }
            printf("%3i %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n",
                   input[i].doy,input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                   input[i].stl, input[i].f107a, input[i].f107 );
            printf("%11.7f %11.7f %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g \n",
                     t[1],t[2],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]);
          }

	printf("now do the msis-90 test cases --------------------------------------\n");
	MSIS_Vers::msis90init(msis90r);

	/* evaluate std test cases */
  	for (i = 0; i < 17; i++)
          {
            printf("------case number %3i\n",i);
            if (i == 2)  // it's different for the msis00 test case
                input[i].alt = 400.0;
            if (i == 9)
				MSIS_Vers::gtd6(msis90r, lpoly, fit, lsqv,
                     input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                     input[i].stl, input[i].f107a, input[i].f107, apar9, mass, d, t);
              else
              {
                if (i >= 15)
                  {
// test switch settings
					MSIS_Vers::tselec(msis90r.csw, svt);
					MSIS_Vers::gtd6(msis90r, lpoly, fit, lsqv,
                         input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                         input[i].stl, input[i].f107a, input[i].f107, aph, mass, d, t);
                  }
                  else
					  MSIS_Vers::gtd6(msis90r, lpoly, fit, lsqv,
                         input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                         input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
              }
            printf("%3i %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n",
                   input[i].doy,input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                   input[i].stl, input[i].f107a, input[i].f107 );
            printf("%11.7f %11.7f %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g \n",
                     t[1],t[2],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]);
          }

	printf("now do the msis-00 test cases --------------------------------------\n");
	MSIS_Vers::msis00init(msis00r);

	/* evaluate std test cases */
  	for (i = 0; i < 17; i++)
          {
            printf("------case number %3i\n",i);
            if (i ==2)   // alt above 500 km
              {
                input[i].alt = 1000.0;
                printf("-- high alt, gtd7 results \n");
				MSIS_Vers::gtd7(msis00r, lpoly, fit, lsqv,
                   input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                   input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
                printf("%3i %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n",
                       input[i].doy,input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                       input[i].stl, input[i].f107a, input[i].f107 );
                printf("%11.7f %11.7f %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g \n",
                         t[1],t[2],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]);
                printf("-- high alt, gtd7d results \n");
				MSIS_Vers::gtd7d(msis00r, lpoly, fit, lsqv,
                   input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                   input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
              }
              else
                if (i ==9)
					MSIS_Vers::gtd7(msis00r, lpoly, fit, lsqv,
                       input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                       input[i].stl, input[i].f107a, input[i].f107, apar9, mass, d, t);
                   else
                   {
                    if (i >=15)
                      {
// test switch settings
						MSIS_Vers::tselec(msis00r.csw, svt);
						MSIS_Vers::gtd7(msis00r, lpoly, fit, lsqv,
                           input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                           input[i].stl, input[i].f107a, input[i].f107, aph, mass, d, t);
                        }
                      else
						  MSIS_Vers::gtd7(msis00r, lpoly, fit, lsqv,
                           input[i].doy, input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                           input[i].stl, input[i].f107a, input[i].f107, apar, mass, d, t);
                   }
            printf("%3i %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n",
                   input[i].doy,input[i].sec, input[i].alt, input[i].glat, input[i].glong,
                   input[i].stl, input[i].f107a, input[i].f107 );
            printf("%11.7f %11.7f %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g %11.7g \n",
                     t[1],t[2],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8]);
          }

	printf("\n");
        }


void mainplugin
       (
         double r[3], double v[3], double jd
       )
{
	int i, mass;
	msistype msis00r;
        double apar[8];
        double f107, f107a, mfme;
        char interp, inputtype, fluxtype, f81type;

        double f107ctr81;
        double ap;
        double apavg;
        double aparray[8];
        double kp;
        double sumkp;
        double kparray[8];
        spwdata spwarr[2500];
        double jdspwstart;

        mass = 48;

//	doy   =   172;
//   	year  =  2005;
//  	sec   = 29000;
//	alt   =   400;
//	glat  =    60;
//	glong =   -70;
//	stl   =     6;


	/* get space weather data values */
  	for (i=0; i<8; i++)
            apar[i]  =   4;
      	f107a =   150;
	f107  =   150;

        //    interp      - interpolation       n-none, a-ap only, f-f10.7 only, b-both
        interp = 'n';
        //    inputtype   - input type          a-actual, u - user   c - constant
        inputtype = 'a'; // actual

        fluxtype = 'a';  // adjusted
        f81type = 'c';   // centered

//        initeop( eoparr, jdeopsatrt);
//		EopSpw::initspw(spwarr, "D:/Codes/LIBRARY/CPP/TestMSIS/Data/SpaceWeather-All-v1.2_11-03-2014.txt",  jdspwstart);

        jd = 2452363.5;
        r[0] = 1;
        r[1] = 2;
        r[2] = 3;

        // initialize the msis routine - probably do in the initialize part of the code for plug-in
		MSIS_Vers::msis00init(msis00r);

		mfme = 0.0;
//        EopSpw::findatmosparam( jd, mfme, interp, fluxtype, f81type, inputtype, spwarr, jdspwstart, f107, f107ctr81,
//                        ap, apavg, aparray, kp, sumkp, kparray );

 //       interfaceatmos( jd, mfme, r, interp, fluxtype, f81type, inputtype, msis00r, spwarr, jdspwstart );
   }


/*     test case for running jb2006
c     input to jb2006:
c           amjd   : date and time, in modified julian days
c                    and fraction (mjd =3d jd-2400000.5)
c           sun(1) : right ascension of sun (radians)
c           sun(2) : declination of sun (radians)
c           sat(1) : right ascension of position (radians)
c           sat(2) : geocentric latitude of position (radians)
c           sat(3) : height of position (km)
c           geo(1) : 10.7-cm solar flux (1.0e-22*watt/(m**2*hertz))
c                    (tabular time 1.0 day earlier)
c           geo(2) : 10.7-cm solar flux, ave.
c                    81-day centered on the input time
c           geo(3) : geomagnetic planetary 3-hour index
c                    a-sub-p for a tabular time 0.279 days earlier
c                    (6.7 hours earlier)
c           s10    : euv index (26-34 nm) scaled to f10
c                    (tabular time 1.0 day earlier)
c           s10b   : euv 81-day ave. centered index
c           xm10   : mg2 index scaled to f10
c           xm10b  : mg2 81-day ave. centered index
c                    (tabular time 5.0 days earlier)
c
c     output from jb2006:
c           temp(1): exospheric temperature above input position (deg k)
c           temp(2): temperature at input position (deg k)
c           rho    : total mass-desnity at input position (kg/m**3)
*/

void testjb06()
   {
      //const double pi       = 3.1415926535879323846;
      double sun[3],sat[4],geo[4],temp1[4];
//      double  pi = 3.1415927;
      double s10, s10b, f10, f10b, xm10, xm10b, ap, d1950, amjd, rho;
      int iyr, iyy, iday, id1950;

//    set solar indices
//     use 1 day lag for euv and f10 influence
      s10  = 140;
      s10b = 100;
      f10  = 135;
      f10b = 95;

//     use 5 day lag for mg fuv influence
      xm10  = 130;
      xm10b = 95;

//     use 6.7 hr lag for ap influence
      ap = 30;
      geo[1] = f10;
      geo[2] = f10b;
      geo[3] = ap;

//     set point of interest location (radians and km)
      sat[1] = 90.0*pi/180.0;
      sat[2] = 45.0*pi/180.0;
      sat[3] = 400.0;

//     set sun location (radians)
      sun[1] = 90.0*pi/180.0;
      sun[2] = 20.0*pi/180.0;

//    set time of interest
      iyr  = 1;
      iday = 200;
      if (iyr < 50) iyr = iyr + 100;
      iyy  = (iyr-50)*365 + ((iyr-1)/4-12);
      id1950 = iyy + iday;
      d1950  = id1950;
      amjd   = d1950 + 33281.0;

//    compute density kg/m3 rho
//      jb2006 (amjd,sun,sat,geo,s10,s10b,xm10,xm10b,temp1, rho);

	  rho = 0.0;
      printf("%8.0f  %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %10.1f %10.1f %14.8g  \n",
               d1950,s10,s10b,f10,f10b,xm10,xm10b,ap,temp1[1],temp1[2],rho);
//     output results:
// 18828.  140. 100. 135.  95. 130.  95.  30.  1145.8  1137.7 0.4066d-11

}

