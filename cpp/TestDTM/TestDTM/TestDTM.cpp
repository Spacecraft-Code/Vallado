// TestDTM.cpp : main project file.

/*
 * test driver for DTM-2012 model
 * dav 6 aug 2013
 *
*/

// system includes
//#include "stdafx.h" // important - put this one first!!!

#include <stdio.h>
#include <string.h>
#include <vector>

#pragma once


// dtm includes
#include "astTime.h"
#include "astMath.h"
#include "ast2Body.h"
#include "EopSpw.h"
#include "DTM_12.h"


// local functions
// test the routine against the dtm test cases

//#define pi 3.14159265358979323846

void main
	(
	);

// example for plugin to HPOP
void mainplugin
	( 
        double recef[3], double vecef[3], double jdut1, double dut1
    );

//void ijk2ll
//     (
//       double recef[3], double jdut1,
//       double& latgc, double& latgd, double& lon, double& hellp
//     );
//
//void gc_gd
//     (
//       double&    latgc ,
//       edirection direct,
//       double&    latgd
//     );
//
//double  mag
//        (
//          double x[3]
//        );
//
//double  gstime
//        (
//          double jdut1
//        );


using namespace std; // it allows you to do things like cin>> cout<< etc.
//using namespace System;

//----------------------------------------------------------------------
//                                                                            //
//> @brief           Parameters for the DTM wrapper and DTM 2012 subroutines
//                                                                            //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                                                                            //
// Software context            : ATMOP                                        //
// Library Responsible         : Noelia Sanchez                               //
//>Subroutine Author           : @author Raul Dominguez
// Company                     : DEIMOS Space S.L.                            //
// Programming Language        : Fortran 90                                   //
// Associated File Name        : dtm_wrapper.f90                              //
// Development Compiler and OS : Windows 7 32 bit, cygwin                     //
// Compiler Version            : GNU gfortran 4.5.3                           //
// Compiling Options           : standard                                     //
//                                                                            //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                                                                            //
//  Method                                                                    //
//  ======
//>  @details
//>   This module allows the user to set certain parameters required for the
//>   dtm2012 and dtm_wrapper modules (namely, unit numbers and path to files)
//>   To do so, certain variables are exported.
//>   The user can change the values stored in these variables to choose
//>   custom units and file paths.
//>   In case the user does not set these values, default values are provided
//                                                                            //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                                                                            //
// Function history                                                           //
//                                                                            //
//> @version 1.0
//> @date 17/09/2012
//  Initial coding
//                                                                            //
//----------------------------------------------------------------------


//----------------------------------------------------------------------
//                                                                            //
//> @brief           Wrapper for the dtm 2012 subroutine
//                                                                            //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                                                                            //
// Software context            : ATMOP                                        //
// Library Responsible         : Noelia Sanchez                               //
//>Subroutine Author           : @author Raul Dominguez
// Company                     : DEIMOS Space S.L.                            //
// Programming Language        : Fortran 90                                   //
// Associated File Name        : dtm_wrapper.f90                              //
// Development Compiler and OS : Windows 7 32 bit, cygwin                     //
// Compiler Version            : GNU gfortran 4.5.3                           //
// Compiling Options           : standard                                     //
//                                                                            //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                                                                            //
//  Method                                                                    //
//  ======
//>  @details
//>   This module is a wrapper for the DTM2012 subroutine. Its main purpose is
//>   to allow the user to call the DTM2012 subroutine giving only the date,
//>   the longitude and latitude and the altitude. The computation of local
//>   hour, day of year, and the computations related to the sun and
//>   geomagnetic proxies are computed automatically within it
//                                                                            //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                                                                            //
//                                                                            //
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//                                                                            //
// Function history                                                           //
// ================                                                           //
//                                                                            //
//   @version 1.0
//   @date    31/07/12
//   First coding
//
//   @version 1.1
//   @date    10/08/12
//   Added "test mode" to allow testing the wrapper with user-fed f,fbar, akp
//         Added a custom type "date" to simplify input to the wrapper
//
//   @version 1.2
//   @date    14/09/12
//   Changed module so it only works with real variables (no double to real
//         conversions anymore)
//
//>  @version 1.3
//>  @date    03/10/12
//   Deprecated test mode                                                                            //
//----------------------------------------------------------------------


//public components. Only these variables and subroutines can be seen from outside the dtm_wrap module
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//public :: dtm_wrapper                                                    //wrap the DTM2012 routine
//public :: local_solar_time                                               //allows to compute the local solar time anywhere
//public :: day_of_year                                                    //allows to compute the day of year anywhere
//public :: dtm_date                                                       //export the type(date)
//public :: process_input_date                                             //export the routine which fills the missing fields of the type(date)
//public :: f_out, fbar_out, akp_out                                       //export interesting intermediate results
//public :: hl_out,dayofyear_out                                           //export interesting intermediate results
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//




/* -----------------------------------------------------------------------------
*
*                           function main
*
*  this function does all the claculations necessary to solve the dtm 
*  atmospoheric model. Specific parametersare read, space weather data is input,
*  the specific test case is entered,values are interpolated, and the answer 
*  (and uncertainty) are found. 
*
*  author        : sean bruinsma                                       2012
*
*  revisions
*    conversion form for95 to c++
*    note that all arrays are larger by +1 from the for95 code to permit
*    the indices to be the same in both codes
*                : david vallado                  719-573-2600   5 aug 2013
*
*  inputs          description                    range / units
*    latgd        - geodetic latitude              deg
*    hellp        - altitude                       km
*    lon        - longitude                      deg
*    today       - structure whholding YMD HMS time information
*
*  outputs       :
*    ro          - density at the given altitude  kg/m3
*    ro_unc      - density uncertainty            %
*    tinf        - exospheric temperature
*    tz          - temperature at the given altitude
*    tp120       - vertical temperature gradient at 120 km
*    wmm         - mean molecular mass
*    d[7]        - concentrations of
*                   1 atomic hydrogen
*                   2 helium
*                   3 atomic oxygen
*                   4 molecular nitrogen
*                   5 molecular oxygen
*                   6 atomic nitrogen (currently unused)
*    eoparr      - array of eop data records
*    jdeopstart  - julian date of the start of the eoparr data
*
*  locals        :
*                -
*
*  coupling      :
*
*  references    :
*
*  -------------------------------------------------------------------------- */

void main()
{
    // Math parameters
    double deg2rad; //< Convert from degrees to radians
    double rad2deg; //< Convert from radians to degress
    double Pi;      //< Constant Pi value

	// initialize all structures - not all seem to be used...
	dtm_daterectype today;
	plgdtmtype plgdtm;
	hlocaltype hlocal;
	eclipttype eclipt;
	datmotype datmo;
	constype cons;
	pardtmtype pardtm;

	// ---------------------- density uncertainty files and parameters
    char PATH_TO_DEN_UNC_1[100];  //< Full Path (absolute or relative) to this file
                                            //< (Default: data/smoothrelvar_350km_kpmax3)
    char PATH_TO_DEN_UNC_2[100];  //< Path (absolute or relative) to this file
                                            //< (Default: data/smoothratio_kpmin4-to-max3)
    char PATH_TO_DEN_UNC_3[100];  //< Path (absolute or relative) to this file
                                            //< (Default: data/smoothratio_500to350km)
    static double nomgrid[25][19];
    static double alt_scale[25][19];
    static double kp_scale[25][19];
    
	// ---------------------- Space Weather geomagnetic and solar flux data
    char PATH_TO_DTM2012_IN_FILE[100];  //< Full Path (absolute or relative) to this file
                                                   //< (Default: data/dtm_2012_NF.dat)
    char PATH_TO_DTM2012_OUT_FILE[100]; //< Path (absolute or relative) to this file
                                                   //< (Default: ZVERIF_xxx.DAT)
    char path_to_a_file[100]; //< Path (absolute or relative) to this file
                                             //< (Default: data/am_file_spider.dat)
    char path_to_f_file[100]; //< Path (absolute or relative) to this file
                                            //< (Default: data/proxies_unadjusted.dat)
	a_structtype a_indexes[3000];
	f_structtype f_indexes[3000];
	int number_of_lines_in_a_file, number_of_lines_in_f_file;

	double hellp, latgd, lat_rads, lon_rads, lon; 
	double s_soltime_hours;
	int i;
    double ro, tinf, tz, tp120, wmm, ro_unc, d[7];   
	
	static int IVERIF = 0;  //// IVERIF = 0/1 >> print model in u_Out

    // Auxiliary variables. We use them to report the results of the proxy file
    // interpolation outside the module
    // cdav increment arrays by 1 to have the same indices the same between languages
    double f_out[3], fbar_out[3];
	double f[3], fbar[3], akp[5];
	double akp_out[5];
    double hl_out, dayofyear_out;
	double dayofyear, hl;


    // ---------------------- begin code ------------------------
    // math constants
	deg2rad = 1.74532925199e-2;
    rad2deg = 57.2957795131;
    Pi      = 3.14159265358979323846;

	// initialize variables for the density uncertainty and dtm
	strcpy_s(PATH_TO_DTM2012_IN_FILE, sizeof(char), "");
	strcpy_s(PATH_TO_DTM2012_IN_FILE, sizeof(PATH_TO_DTM2012_IN_FILE), "D:/Codes/LIBRARY/CPP/TestDTM/Data/dtm_2012_NF.dat");
	strcpy_s(PATH_TO_DTM2012_OUT_FILE, sizeof(char), "");
	strcpy_s(PATH_TO_DTM2012_OUT_FILE, sizeof(PATH_TO_DTM2012_OUT_FILE), "D:/Codes/LIBRARY/CPP/TestDTM/Data/dtm_2012_NF.out");
	strcpy_s(PATH_TO_DEN_UNC_1, sizeof(char), "");
	strcpy_s(PATH_TO_DEN_UNC_1, sizeof(PATH_TO_DEN_UNC_1), "D:/Codes/LIBRARY/CPP/TestDTM/Data/smoothrelvar_350km_kpmax3");
	strcpy_s(PATH_TO_DEN_UNC_2, sizeof(char), "");
	strcpy_s(PATH_TO_DEN_UNC_2, sizeof(PATH_TO_DEN_UNC_2), "D:/Codes/LIBRARY/CPP/TestDTM/Data/smoothratio_kpmin4-to-max3");
	strcpy_s(PATH_TO_DEN_UNC_3, sizeof(char), "");
	strcpy_s(PATH_TO_DEN_UNC_3, sizeof(PATH_TO_DEN_UNC_3), "D:/Codes/LIBRARY/CPP/TestDTM/Data/smoothratio_500to350km");

	DTM_12::DTMInit(PATH_TO_DEN_UNC_1, PATH_TO_DEN_UNC_2, PATH_TO_DEN_UNC_3, PATH_TO_DTM2012_IN_FILE, PATH_TO_DTM2012_OUT_FILE, 
		     IVERIF, nomgrid, kp_scale, alt_scale,  pardtm );

    // initialize Space Weather geomagnetic and solar flux data
#ifdef _MSC_VER
	strcpy_s(path_to_a_file, sizeof(path_to_a_file), "");
	strcpy_s(path_to_a_file, sizeof(path_to_a_file), "D:/Codes/LIBRARY/CPP/TestDTM/Data/am_file_spiderx.dat");  // geomagnetic indices
	strcpy_s(path_to_f_file, sizeof(path_to_f_file), "");
	strcpy_s(path_to_f_file, sizeof(path_to_f_file), "D:/Codes/LIBRARY/CPP/TestDTM/Data/proxies_unadjustedx.dat");  // solar flux values
#else
	strcpy(path_to_a_file, "");
	strcpy(path_to_a_file, "D:/Codes/LIBRARY/CPP/TestDTM/Data/am_file_spiderx.dat");  // geomagnetic indices
	strcpy(path_to_f_file, "");
	strcpy(path_to_f_file, "D:/Codes/LIBRARY/CPP/TestDTM/Data/proxies_unadjustedx.dat");  // solar flux values
#endif
	DTM_12::InitSPW(path_to_a_file, path_to_f_file,
		     number_of_lines_in_a_file,  number_of_lines_in_f_file, a_indexes, f_indexes );

	// ----------------- setup specific test case
    hellp = 200.0; //real
    lon =  -4.0;  // deg
    latgd =  40.0;
	//longitude and latitude in radians
    lat_rads = latgd * deg2rad;
    lon_rads = lon * deg2rad;

    today.type_flag=2; // Set this flag to 2 to input a calendar date
    today.day    = 3;
    today.month  = 9;
    today.year   = 2011;
    today.hour   =  8;
    today.minute = 32;
    today.second = 34.23; // warning// real*8

	// get space weather parameters for specific test case date
	DTM_12::dtm_wrapper(today, lon, number_of_lines_in_a_file, number_of_lines_in_f_file,
	             a_indexes, f_indexes, f, fbar, akp, dayofyear, hl );

    // ------------------------------- dtm2012 subroutine  
	DTM_12::dtm2012(plgdtm, hlocal, eclipt, datmo, cons, pardtm,
		            dayofyear, f, fbar, akp, hellp, hl, lat_rads, lon_rads, tz, tinf, tp120, ro, d, wmm);


    // get the density uncertainty
    // note that the local solar time must be in radians
    s_soltime_hours = hl*24.0/(2.0*Pi);
    if (s_soltime_hours >= 23.999) 
          s_soltime_hours = 0.0;
	DTM_12::density_uncertainty(nomgrid, alt_scale, kp_scale, hellp, lat_rads, s_soltime_hours,
		                 f[1], akp[1], ro_unc);

    //write the results to public global variables
    for (i = 1; i < 5; i++)
	  {
          if (i <=3)
		  {
    		  f_out[i] = f[i];
              fbar_out[i] = fbar[i];
		  }
          akp_out[i] = akp[i];
	  }
    dayofyear_out = dayofyear;
    hl_out = s_soltime_hours; //required in hours instead of radians
	
	// output the results
    printf("\n -----------CALL TO DTM2012----------- \n Inputs \n" );
    printf("    Date: %12.7f ( %2i/%2i/%4i ) at %2i:%2i:%12.9f )  \n", today.mjd2000, 
                                         today.day, today.month,today.year, today.hour, today.minute, today.second );
    printf("   day of year            %10.6f \n", dayofyear_out );
    printf("   local time (hours)     %12.7f \n", hl_out );
    printf("   latitude               %12.7f \n", latgd );
    printf("   longitude              %12.7f \n", lon );
    printf("   altitude               %12.7f \n", hellp );
    printf("   f          %12.7lf  %12.7f \n", f_out[1], f_out[2] );   //This variable is provided by dtm wrapper
    printf("   fbar       %12.7lf  %12.7f \n", fbar_out[1], fbar_out[2] );  //This variable is provided by dtm wrapper
    printf("   akp        %12.7lf  %12.7f %12.7f %12.7f \n", akp_out[1],akp_out[2],akp_out[3],akp_out[4] );   //This variable is provided by dtm wrapper
    printf(" \n outputs \n");
    printf("   Temp at altitude       %12.7f \n", tz );
    printf("   exospheric tmp         %12.7f \n", tinf );
    printf("   atomic hydrogen        %16.7g \n", d[1] );
    printf("   helium                 %16.7g \n", d[2] );
    printf("   atomic oxygen          %16.7g \n", d[3] );
    printf("   molecular nitro        %16.7g \n", d[4] );
    printf("   molecular oxygen       %16.7g \n", d[5] );
    printf("   density (g/cm^3)       %16.7g    %6.2f pct \n", ro, ro_unc );
    //write(*,'("   atomic nitrogen : ", ES16.7)') d[6] Ignore this value
    printf("   mean molec mass        %16.7f \n", wmm );

	//--------------------
	// test hpop like implementation
	//--------------------
	double recef[3], vecef[3], jdutc, dut1;
	recef[0] = -605.790430800;
	recef[1] = -5870.23040700;
	recef[2] =  3493.05200400;
	vecef[0] = -1.568251615;
	vecef[1] = -3.702348353;
	vecef[2] = -6.479484915;
	DTM_12::JD2000(jdutc, today.year, today.month, today.day, today.hour, today.minute, today.second);
	dut1 = 0.6453; // need to update this for each date selected, but it's small
//	mainplugin(recef, vecef, jdutc, dut1);

}   // end main



// test a simulated HPOP interface using position vector and julian date
void mainplugin
       (
         double recef[3], double vecef[3], double jdutc, double dut1
       )
{
    double f107, mfme;
    char interp, inputtype, fluxtype, f81type;

    double f107ctr81;
    double ap;
    double apavg;
    double aparray[8];
    double kp;
    double sumkp;
    double kparray[8];
	std::vector< spwdata > spwarr(spwsize);
	double jdspwstart, jdspwstartF;

    // Math parameters
    double deg2rad; //< Convert from degrees to radians
    double rad2deg; //< Convert from radians to degress
    double Pi;      //< Constant Pi value

	// initialize all structures - not all seem to be used...
	dtm_daterectype today;
	plgdtmtype plgdtm;
	hlocaltype hlocal;
	eclipttype eclipt;
	datmotype datmo;
	constype cons;
	pardtmtype pardtm;

	// ---------------------- density uncertainty files and parameters
    char PATH_TO_DEN_UNC_1[100];  //< Full Path (absolute or relative) to this file
                                            //< (Default: data/smoothrelvar_350km_kpmax3)
    char PATH_TO_DEN_UNC_2[100];  //< Path (absolute or relative) to this file
                                            //< (Default: data/smoothratio_kpmin4-to-max3)
    char PATH_TO_DEN_UNC_3[100];  //< Path (absolute or relative) to this file
                                            //< (Default: data/smoothratio_500to350km)
    static double nomgrid[25][19];
    static double alt_scale[25][19];
    static double kp_scale[25][19];
    
	// ---------------------- Space Weather geomagnetic and solar flux data
    char PATH_TO_DTM2012_IN_FILE[100];  //< Full Path (absolute or relative) to this file
                                                   //< (Default: data/dtm_2012_NF.dat)
    char PATH_TO_DTM2012_OUT_FILE[100]; //< Path (absolute or relative) to this file
                                                   //< (Default: ZVERIF_xxx.DAT)
//    char path_to_a_file[100]; //< Path (absolute or relative) to this file
                                             //< (Default: data/am_file_spider.dat)
//    char path_to_f_file[100]; //< Path (absolute or relative) to this file
                                            //< (Default: data/proxies_unadjusted.dat)
//	a_structtype a_indexes[3000];
//	f_structtype f_indexes[3000];

	double hellp, latgd, latgc, lat_rads, lon_rads, lon; 
	double s_soltime_hours;
    double ro, tinf, tz, tp120, wmm, ro_unc, d[7];   
	
	static int IVERIF = 0;  //// IVERIF = 0/1 >> print model in u_Out

    // Auxiliary variables. We use them to report the results of the proxy file
    // interpolation outside the module
    // cdav increment arrays by 1 to have the same indices the same between languages
	double f[3], fbar[3], akp[5];
	double dayofyear, hl;


    // ---------------------- begin code ------------------------
    // math constants
	deg2rad = 1.74532925199e-2;
    rad2deg = 57.2957795131;
    Pi      = 3.14159265358979323846;

	// initialize variables for the density uncertainty and dtm
    strcpy_s( PATH_TO_DTM2012_IN_FILE,"");
    strcpy_s( PATH_TO_DTM2012_IN_FILE,"D:/Codes/LIBRARY/CPP/TestDTM/Data/dtm_2012_NF.dat");
    strcpy_s( PATH_TO_DTM2012_OUT_FILE,"");
    strcpy_s( PATH_TO_DTM2012_OUT_FILE,"D:/Codes/LIBRARY/CPP/TestDTM/Data/dtm_2012_NF.out");
    strcpy_s( PATH_TO_DEN_UNC_1,"");
    strcpy_s( PATH_TO_DEN_UNC_1,"D:/Codes/LIBRARY/CPP/TestDTM/Data/smoothrelvar_350km_kpmax3");
    strcpy_s( PATH_TO_DEN_UNC_2,"");
    strcpy_s( PATH_TO_DEN_UNC_2,"D:/Codes/LIBRARY/CPP/TestDTM/Data/smoothratio_kpmin4-to-max3");
    strcpy_s( PATH_TO_DEN_UNC_3,"");
    strcpy_s( PATH_TO_DEN_UNC_3,"D:/Codes/LIBRARY/CPP/TestDTM/Data/smoothratio_500to350km");
    
	DTM_12::DTMInit(PATH_TO_DEN_UNC_1, PATH_TO_DEN_UNC_2, PATH_TO_DEN_UNC_3,
		     PATH_TO_DTM2012_IN_FILE, PATH_TO_DTM2012_OUT_FILE, IVERIF,
			 nomgrid, kp_scale, alt_scale,  pardtm );

    // initialize Space Weather geomagnetic and solar flux data
//	strcpy_s( path_to_a_file,"");
//    strcpy_s( path_to_a_file,"D:/Codes/LIBRARY/CPP/TestDTM/Data/am_file_spider.dat");  // geomagnetic indices
//    strcpy_s( path_to_f_file,"");
//    strcpy_s( path_to_f_file,"D:/Codes/LIBRARY/CPP/TestDTM/Data/proxies_unadjusted.dat");  // solar flux values
//	InitSPW( path_to_a_file, path_to_f_file, 
//		     number_of_lines_in_a_file,  number_of_lines_in_f_file, a_indexes, f_indexes );
	EopSpw::initspw(spwarr, "D:/Codes/LIBRARY/CPP/TestDTM/Data/SpaceWeather-All-v1.2_11-03-2014.txt", jdspwstart, jdspwstartF);
    
	// ----------------- setup specific test case
	// convert pos/vel to lat,lon,alt
	ast2Body::ijk2ll( recef, jdutc + dut1/86400.0, latgc, latgd, lon, hellp );

	//longitude and latitude already in radians
    lat_rads = latgd;
    lon_rads = lon;

    DTM_12::DJ2000( jdutc, today.year, today.month, today.day, today.hour, today.minute, today.second );
    today.type_flag=2; // Set this flag to 2 to input a calendar date

	// get space weather parameters for specific test case date
//    dtm_wrapper( today, lon, number_of_lines_in_a_file, number_of_lines_in_f_file, 
//	             a_indexes, f_indexes, f, fbar, akp, dayofyear, hl );
	interp = 'l';  // linear
	f81type = 'c';  // centered
	fluxtype = 'o';  // observed
	inputtype = 'a';  // actual data
	mfme = 0;

	EopSpw::findatmosparam(jdutc, mfme, interp, fluxtype, f81type, inputtype, spwarr, jdspwstart, f107, f107ctr81,
                        ap, apavg, aparray, kp, sumkp, kparray );
	
    // ------------------------------- dtm2012 subroutine  
	DTM_12::dtm2012(plgdtm, hlocal, eclipt, datmo, cons, pardtm,
		      dayofyear, f, fbar, akp, hellp, hl, lat_rads, lon_rads, tz, tinf, tp120, ro, d, wmm);

    // get the density uncertainty
    // note that the local solar time must be in radians
    s_soltime_hours = hl*24.0/(2.0*pi);
    if (s_soltime_hours >= 23.999) 
          s_soltime_hours = 0.0;
	DTM_12::density_uncertainty(nomgrid, alt_scale, kp_scale, hellp, lat_rads, s_soltime_hours,
		                 f[1], akp[1], ro_unc);


	// output the results
    printf("\n -----------CALL TO DTM2012----------- \n Inputs \n" );
    printf("    Date: %12.7f ( %2i/%2i/%4i ) at %2i:%2i:%12.9f )  \n", today.mjd2000, 
                                         today.day, today.month,today.year, today.hour, today.minute, today.second );
    printf("   day of year            %10.6f \n", dayofyear );
    printf("   local time (hours)     %12.7f \n", s_soltime_hours );
    printf("   latitude               %12.7f \n", latgd );
    printf("   longitude              %12.7f \n", lon );
    printf("   altitude               %12.7f \n", hellp );
    printf("   f          %12.7lf  %12.7f \n", f[1], f[2] );   //This variable is provided by dtm wrapper
    printf("   fbar       %12.7lf  %12.7f \n", fbar[1], fbar[2] );  //This variable is provided by dtm wrapper
    printf("   akp        %12.7lf  %12.7f %12.7f %12.7f \n", akp[1],akp[2],akp[3],akp[4] );   //This variable is provided by dtm wrapper
    printf(" \n outputs \n");
    printf("   Temp at altitude       %12.7f \n", tz );
    printf("   exospheric tmp         %12.7f \n", tinf );
    printf("   atomic hydrogen        %16.7g \n", d[1] );
    printf("   helium                 %16.7g \n", d[2] );
    printf("   atomic oxygen          %16.7g \n", d[3] );
    printf("   molecular nitro        %16.7g \n", d[4] );
    printf("   molecular oxygen       %16.7g \n", d[5] );
    printf("   density (g/cm^3)       %16.7g    %6.2f pct \n", ro, ro_unc );
    //write(*,'("   atomic nitrogen : ", ES16.7)') d[6] Ignore this value
    printf("   mean molec mass        %16.7f \n", wmm );

} // end mainplugin

//
///* -----------------------------------------------------------------------------
//*
//*                           function ijk2ll
//*
//*  these subroutines convert a geocentric equatorial position vector into
//*    latitude and longitude.  geodetic and geocentric latitude are found. the
//*    inputs must be ecef.
//*
//*  author        : david vallado                  719-573-2600    6 dec 2005
//*
//*  revisions
//*
//*  inputs          description                    range / units
//*    recef       - ecef position vector           km
//*    jdut1       - julian date (ut1)              days from 4713 bc
//*
//*  outputs       :
//*    latgc       - geocentric latitude            -Pi to Pi rad
//*    latgd       - geodetic latitude              -Pi to Pi rad
//*    lon         - longitude (west -)             -2pi to 2pi rad
//*    hellp       - height above the ellipsoid     km
//*
//*  locals        :
//*    temp        - diff between geocentric/
//*                  geodetic lat                   rad
//*    gst         - greenwich sidereal time        rad
//*    sintemp     - sine of temp                   rad
//*    olddelta    - previous value of deltalat     rad
//*    rtasc       - right ascension                rad
//*    decl        - declination                    rad
//*    i           - index
//*
//*  coupling      :
//*    mag         - magnitude of a vector
//*    gstime      - greenwich sidereal time
//*    gcgd        - converts between geocentric and geodetic latitude
//*
//*  references    :
//*    vallado       2007, 179-180, alg 12 and alg 13, ex 3-3
//* --------------------------------------------------------------------------- */
//
//void ijk2ll
//     (
//       double recef[3], double jdut1,
//       double& latgc, double& latgd, double& lon, double& hellp
//     )
//     {
//       const double twopi      =    2.0 * Pi;
//       const double small      =    0.00000001;         // small value for tolerances
//       const double re         = 6378.137;
//       const double eesqrd     =    0.006694385000;     // eccentricity of earth sqrd
//       double magr, gst, decl, rtasc, olddelta, temp, sintemp, s, c = 0.0;
//       int i;
//
//        // ---------------------------  implementation   -----------------------
//        magr = mag( recef );
//
//        // ---------------------- find longitude value  ------------------------
//        temp = sqrt( recef[0]*recef[0] + recef[1]*recef[1] );
//        if ( fabs( temp ) < small )
//            rtasc= sgn(recef[2])*Pi*0.5;
//          else
//            rtasc= atan2( recef[1], recef[0] );
//
//        gst  = gstime( jdut1 );
//        lon  = rtasc - gst;
//        if ( fabs(lon) >= Pi )   // mod it ?
//          {
//            if ( lon < 0.0  )
//                lon= twopi + lon;
//              else
//                lon= lon - twopi;
//
//          }
//        decl = asin( recef[2] / magr );
//        latgd= decl;
//
//        // ----------------- iterate to find geodetic latitude -----------------
//        i = 1;
//        olddelta = latgd + 10.0;
//
//        while ((fabs(olddelta - latgd) >= small) && (i<10))
//          {
//            olddelta= latgd;
//            sintemp = sin( latgd );
//            c       = re  / (sqrt( 1.0 - eesqrd*sintemp*sintemp ));
//            latgd   = atan( (recef[2] + c*eesqrd*sintemp)/temp );
//            i = i + 1;
//          }
//
//        if ((Pi*0.5 - fabs(latgd)) > Pi/180.0)  // 1 deg
//           hellp   = (temp/cos(latgd)) - c;
//         else
//         {
//           s = c * (1.0 - eesqrd);
//           hellp   = recef[2]/sin(latgd) - s;
//         }
//
//        gc_gd(latgc, eFrom, latgd);
//   }   // function ijk2ll
//
//
///* -----------------------------------------------------------------------------
//*
//*                           function gc_gd
//*
//*  this function converts from geodetic to geocentric latitude for positions
//*    on the surface of the earth.  notice that (1-f) squared = 1-esqrd.
//*
//*  author        : david vallado                  719-573-2600    6 dec 2005
//*
//*  revisions
//*
//*  inputs          description                    range / units
//*    latgd       - geodetic latitude              -Pi to Pi rad
//*
//*  outputs       :
//*    latgc       - geocentric latitude            -Pi to Pi rad
//*
//*  locals        :
//*    none.
//*
//*  coupling      :
//*    none.
//*
//*  references    :
//*    vallado       2007, 148, eq 3-11
//* --------------------------------------------------------------------------- */
//
//void gc_gd
//     (
//       double&    latgc ,
//       edirection direct,
//       double&    latgd
//     )
//     {
//       const double eesqrd = 0.006694385000;     // eccentricity of earth sqrd
//
//       if (direct == eTo)
//           latgd= atan( tan(latgc)/(1.0  - eesqrd) );
//         else
//           latgc= atan( (1.0  - eesqrd)*tan(latgd) );
//     }   // function gc_gd
//
///* -----------------------------------------------------------------------------
//*
//*                           function mag
//*
//*  this procedure finds the magnitude of a vector.  the tolerance is set to
//*    0.000001, thus the 1.0e-12 for the squared test of underflows.
//*
//*  author        : david vallado                  719-573-2600    1 mar 2001
//*
//*  inputs          description                    range / units
//*    vec       - vector
//*
//*  outputs       :
//*    vec       - answer stored in function return
//*
//*  locals        :
//*    none.
//*
//*  coupling      :
//*    none.
//*
//* --------------------------------------------------------------------------- */
//
//double  mag
//        (
//          double x[3]
//        )
//   {
//     return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
//   }  // end mag
//
//
///* -----------------------------------------------------------------------------
//*
//*                           function gstime
//*
//*  this function finds the greenwich sidereal time (iau-82).
//*
//*  author        : david vallado                  719-573-2600    1 mar 2001
//*
//*  inputs          description                    range / units
//*    jdut1       - julian date in ut1             days from 4713 bc
//*
//*  outputs       :
//*    gstime      - greenwich sidereal time        0 to 2pi rad
//*
//*  locals        :
//*    temp        - temporary variable for doubles   rad
//*    tut1        - julian centuries from the
//*                  jan 1, 2000 12 h epoch (ut1)
//*
//*  coupling      :
//*    none
//*
//*  references    :
//*    vallado       2007, 193, eq 3-43
//*
//* --------------------------------------------------------------------------- */
//
//double  gstime
//        (
//          double jdut1
//        )
//   {
//     const double twopi = 2.0 * Pi;
//     const double deg2rad = Pi / 180.0;
//     double       temp, tut1;
//
//     tut1 = (jdut1 - 2451545.0) / 36525.0;
//     temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
//             (876600.0*3600 + 8640184.812866) * tut1 + 67310.54841;  // sec
//     temp = fmod(temp * deg2rad / 240.0, twopi); //360/86400 = 1/240, to deg, to rad
//
//     // ------------------------ check quadrants ---------------------
//     if (temp < 0.0)
//         temp += twopi;
//
//     return temp;
//   }
//






