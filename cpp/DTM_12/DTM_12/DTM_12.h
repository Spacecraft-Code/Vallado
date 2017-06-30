#ifndef _DTM_12_h_
#define _DTM_12_h_
/*     ----------------------------------------------------------------
*
*                              DTM_12.h
*
*  this file contains routines for the dtm-2012 atmospheric model
*
*     (w) 719-573-2600, email dvallado@agi.com
*
*     *****************************************************************
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*               6 aug 13  david vallado
*                           conversion to c ++
*              10 may 13  sean bruinsma
*                           original baseline
*     *****************************************************************       */

#include <math.h>
#include <io.h>      
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// cdav increment all by one to start at 1 instead of 0
// create a struct for the numebrs since they don't seem to work as globals in c++
typedef struct dtmindextype
{
	static const int Nb_lat        = 7;
    static const int Nb_flux       = 13;
    static const int Nb_kp         = 16;
    static const int Nb_SLat       = 3;
    static const int Nb_SASLat     = 3;
    static const int Nb_NSLat      = 4;
    static const int Nb_SANSLat    = 4;
    static const int Nb_DiAn       = 13;
    static const int Nb_SDiAn      = 13;
    static const int Nb_TDi        = 3;
    static const int Nb_AMg        = 10;
    static const int Nb_Lon        = 11;  // old 11
    static const int Nb_dPhas      = 5;
} dtmindextype;

   static const int Nb_lat          = 7;
   static const int Nb_flux         = 13;
   static const int Nb_kp           = 16;
   static const int Nb_SLat         = 3;
   static const int Nb_SASLat       = 3;
   static const int Nb_NSLat        = 4;
   static const int Nb_SANSLat      = 4;
   static const int Nb_DiAn         = 13;
   static const int Nb_SDiAn        = 13;
   static const int Nb_TDi          = 3;
   static const int Nb_AMg          = 10;
   static const int Nb_Lon          = 11;  // only difference
   static const int Nb_dPhas        = 5;


//   double  DTM_12_T_one;                  //  SUM(T_x * T_dx) except T_dPhas
//   double  DTM_12_T_lat[Nb_lat];          //  in latitude
//   double  DTM_12_T_flux[Nb_flux];        //  in flux
//   double  DTM_12_T_kp[Nb_kp];            //  in kp
//   double  DTM_12_T_SLat[Nb_SLat];        //  annual & symetric in latitude
//   double  DTM_12_T_SASLat[Nb_SASLat];    //  semi-annual & symetric in latitude
//   double  DTM_12_T_NSLat[Nb_NSLat];      //  annual & non-symetric in latitude
//   double  DTM_12_T_SANSLat[Nb_SANSLat];  //  semi-annual  non-symetric in latitude
//   double  DTM_12_T_DiAn[Nb_DiAn];        //  diurnal (& annual coupled)
//   double  DTM_12_T_SDiAn[Nb_SDiAn];      //  semi-diurnal (& annual coupled)
//   double  DTM_12_T_TDi[Nb_TDi];          //  ter-diurnal
//   double  DTM_12_T_AMg[Nb_AMg];          //  activity (Magnetic)
//   double  DTM_12_T_Lon[Nb_Lon];          //  in longitude
//   double  DTM_12_T_dPhas[Nb_dPhas];      //  derivative of phases of the annual and semi-annual terms

typedef struct DTM_12type
  {
	double T_one, T_lat[Nb_lat], T_flux[Nb_flux], T_kp[Nb_kp], T_SLat[Nb_SLat], T_SASLat[Nb_SASLat], 
		   T_NSLat[Nb_NSLat], T_SANSLat[Nb_SANSLat], T_DiAn[Nb_DiAn], T_SDiAn[Nb_SDiAn], T_TDi[Nb_TDi], 
		   T_AMg[Nb_AMg], T_Lon[Nb_Lon], T_dPhas[Nb_dPhas]; 
  } DTM_12type;

//   double  ddtm_12_T_one;                  //  SUM(T_x * T_dx) except T_dPhas
//   double  ddtm_12_T_lat[Nb_lat];          //  in latitude
//   double  ddtm_12_T_flux[Nb_flux];        //  in flux
//   double  ddtm_12_T_kp[Nb_kp];            //  in kp
//   double  ddtm_12_T_SLat[Nb_SLat];        //  annual & symetric in latitude
//   double  ddtm_12_T_SASLat[Nb_SASLat];    //  semi-annual & symetric in latitude
//   double  ddtm_12_T_NSLat[Nb_NSLat];      //  annual & non-symetric in latitude
//   double  ddtm_12_T_SANSLat[Nb_SANSLat];  //  semi-annual  non-symetric in latitude
//   double  ddtm_12_T_DiAn[Nb_DiAn];        //  diurnal (& annual coupled)
//   double  ddtm_12_T_SDiAn[Nb_SDiAn];      //  semi-diurnal (& annual coupled)
//   double  ddtm_12_T_TDi[Nb_TDi];          //  ter-diurnal
//   double  ddtm_12_T_AMg[Nb_AMg];          //  activity (Magnetic)
//   double  ddtm_12_T_Lon[Nb_Lon];          //  in longitude
//   double  ddtm_12_T_dPhas[Nb_dPhas];      //  derivative of phases of the annual and semi-annual terms

typedef struct ddtm_12type
  {
	double T_one, T_lat[Nb_lat], T_flux[Nb_flux], T_kp[Nb_kp], T_SLat[Nb_SLat], T_SASLat[Nb_SASLat], 
		   T_NSLat[Nb_NSLat], T_SANSLat[Nb_SANSLat], T_DiAn[Nb_DiAn], T_SDiAn[Nb_SDiAn], T_TDi[Nb_TDi], 
		   T_AMg[Nb_AMg], T_Lon[Nb_Lon], T_dPhas[Nb_dPhas]; 
  } ddtm_12type;

//.. Common Blocks .. 
typedef struct plgdtmtype
  {
	double p10, p20, p30, p40, p50, p60, p11, p21, p31, p41, p51, p22, p32, p42, p52, 
           p62, p33, p10mg, p20mg, p40mg;
  } plgdtmtype;

typedef struct constype
  {
	double deupi, cdr, sard;
  } constype;

typedef struct datmotype
  {
	int npara, itype, ilin;
  } datmotype;

typedef struct eclipttype
  {
	double cecl, secl, c2ecl, s2ecl, c3ecl, s3ecl, p10ecl, p20ecl, p30ecl, 
           p40ecl, p50ecl, p11ecl, p21ecl, p31ecl, p41ecl, p51ecl, p22ecl, 
           p32ecl, p42ecl, p52ecl, p33ecl;
  } eclipttype;

typedef struct hlocaltype
  {
	double hl0, ch, sh, c2h, s2h, c3h, s3h;
  } hlocaltype;

typedef struct pardtmtype
  {
	 DTM_12type tt, h, he, ox, az2, o2, az, t0, tp;
	 ddtm_12type dtt, dh, dhe, dox, daz2, do2, daz, dt0, dtp;
  }  pardtmtype;


    // additional type definitions
    //
    //Custom type definitions
    //These types contain all the data which is read directly from input files
    // cdav don't need todaytype as the dtm_daterecordtype does this?
	//typedef struct todaytype
    //{
  	//   int day, month, year, hour, minute, type_flag;
	//   double second;
    //} todaytype;

    //>@brief Structure to store data from solar flux "f" file
    typedef struct f_structtype
	{
        int year;      //< Year
        int month;     //< Month
        int day;       //< Day
        double f;      //< F from file
        double fbar;   //< FBAR from file
	} f_structtype;

    //>@brief Structure to store data from geomagnetic "a" file
    // increment by 1 to have the same indices the same between languages
    typedef struct a_structtype
	{
        int year;      //< Year
        int month;     //< Month
        int day;       //< Day
        int dayofyear; //< Day of year (1-366)
        int mean_am;   //< Mean value of a_m
        int am[9];     //< The eight a_m values in each line of the file
    } a_structtype;
		
    //> @brief Structure to hold the date input by the user
    //>This type allows the user to enter either a MJD2000 date or a day/month/year hh:mm:sec
    //>date, without the need of passing a large number of parameters.
    //>The user must set a type_flag and the required variables
    //>Unused variables do not need to be initialized
    typedef struct dtm_daterectype
	{
       int type_flag;  //<1 for MJD2000 date, 2 for calendar date
       double mjd2000; //< Date in MJD2000
       int day;        //< Day
       int month;      //< Month
       int year;       //< Year
       int hour;       //< Hour
       int minute;     //< Minute
       double second;  //< Seconds
	} dtm_daterectype;

	// store the uncertainty parameters
	typedef struct dtm_unctype
    {
       double nomgrid[25][19], kp_scale[25][19], alt_scale[25][19];
    } dtm_unctype;


#pragma once

// version :  2012
//   #define pi 3.14159265358979323846 

	namespace DTM_12 {

		//	public ref class DTM_12Cl
		//	{

		// make sure they are all visible
		//	public:

		// ---------------- routines ----------------
		void dtm2012
			(
			plgdtmtype& plgdtm, hlocaltype& hlocal, eclipttype& eclipt, datmotype& datmo, constype& cons, pardtmtype& pardtm,
			double day, double f[3], double fbar[3], double akp[5], double hellp, double hl, double latgd, double lon,
			double& tz, double& tinf, double& tp120, double& ro, double d[7], double& wmm
			);


		void density_uncertainty
			(
			double nomgrid[25][19], double alt_scale[25][19], double kp_scale[25][19],
			double hellp, double latgd, double lst, double flux, double kp,
			double& unc
			);

		void gldtm_XX
			(
			plgdtmtype& plgdtm, hlocaltype& hlocal, datmotype& datmo,
			double f[3], double fbar[3], double akp[5], double day, DTM_12type& DTM_12, ddtm_12type& ddtm_12,
			double& gdel, double ff0, double lon
			);


		// additional routines used by DTM

		double a2K
			(
			double a
			);

		double local_solar_time
			(
			double longitude, int hour, int minute, double sec
			);

		int leap_day
			(
			int year
			);

		double linear_interpolation
			(
			double x1, double y1, double x2, double y2, double x
			);

		double  hms2hr
			(
			int hr, int mins, double sec
			);

		int day_of_year
			(
			int day, int month, int year
			);

		bool mdyhms_ok
			(
			int day, int month, int year, int hour, int minute,
			double second
			);

		void DJ2000
			(
			double MJDay, int& year, int& mon, int& day, int& hr, int& minute, double& sec
			);

		void JD2000
			(
			double& MJDay, int year, int month, int day, int hr, int minute, double sec
			);

		void fatal_error
			(
			int code
			);

		void warning
			(
			int code
			);

		void process_input_date
			(
			dtm_daterectype& in_date
			);


		void InitSPW
			(
			char path_to_a_file[100], char path_to_f_file[100], int& number_of_lines_in_a_file,
			int& number_of_lines_in_f_file, a_structtype a_indexes[3000], f_structtype f_indexes[3000]
			);

		void DTMInit
			(
			char PATH_TO_DEN_UNC_1[100], char PATH_TO_DEN_UNC_2[100], char PATH_TO_DEN_UNC_3[100],
			char PATH_TO_DTM2012_IN_FILE[100], char PATH_TO_DTM2012_OUT_FILE[100], int IVERIF,
			double nomgrid[25][19], double kp_scale[25][19], double alt_scale[25][19],
			pardtmtype& pardtm
			);

		void dtm_wrapper
			(
			dtm_daterectype& in_date, double lon, int number_of_lines_in_a_file, int number_of_lines_in_f_file,
			a_structtype a_indexes[3000], f_structtype f_indexes[3000],
			double f[3], double fbar[3], double akp[5], double& dayofyear, double& hl
			);


		//	}; // class

	};  // namespace

#endif
