// MSIS_Vers.h
/*      ----------------------------------------------------------------
*
*                              MSIS_Vers.h
*
*  this file contains routines for the msis-86, msis-90, and nrlmsise-00 atmospheric models.
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*              31 mar 08  david vallado
*                           misc updates
*              15 mar 07  david vallado
*                           3rd edition baseline
*               6 aug 04  david vallado
*                           convert to c++
*               6 sep 03  david vallado
*                           fix low alt test cases (long in lpoly)
*              14 feb 03  david vallado
*                           misc updates
*              28 feb 02  david vallado
*                           finish double conversion
*               8 oct 88  nrl
*                           original baseline
*     http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
*
*     *****************************************************************       */

#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

//using namespace System;

	typedef struct gts3crecord
	{
		double tlb, s, db04, db16, db28, db32, db40, db48, db01,
		za, t0, z0, g0, rl, dd, db14, tr12;
	} gts3ctype;

	typedef struct datimerecord
	{
		char isdate[10], istime[10], name[10];
	} datimetype;

	typedef struct cswrecord
	{
		int sw[26], swc[26];
		int isw;
	} cswtype;

	typedef struct dmixrecord
	{
		double dm04, dm16, dm28, dm32, dm40, dm01, dm14;
	} dmixtype;

	typedef struct parmbrecord
	{
		double  gsurf, re;
	} parmbtype;

	typedef struct metselrecord
	{
		int imr;
	} metseltype;

	typedef struct lsqvrecord
	{
		double  qpb[51], dv[61];
		int mp, ii, jg, lt, ierr, ifun, n, j;
	} lsqvtype;

	typedef struct fitrecord
	{
		double taf;
	} fittype;

	typedef struct lpolyrecord
	{
		double plg[10][5], ctloc, stloc, c2tloc, s2tloc, c3tloc, s3tloc,
		day, df, dfa, apd, apdf, apt[5], xlong,
		clong, slong;
		int iyr;
	} lpolytype;

	typedef struct mesorecord
	{
		double tn1[6], tn2[5], tn3[6], tgn1[3], tgn2[3], tgn3[3];
	} mesotype;

	typedef struct lowerrecord
	{
		double ptm[11], pdm[11][9];
	} lowertype;

	typedef struct parmrecord
	{
		double pt[151], pd[151][10], ps[151], pdl[26][3], ptl[101][5],
		pma[101][11], sam[101],
		// these added for function calls when a subarray needed.
		// there is probably an easier way to do this!!!
		pddav[10][151], ptldav[5][101], pmadav[11][101];
	} parmtype;

	// cdav i deleted this because it was just a string manipulation for titles
	//      typedef struct datim7record
	//       {
	//          int isd[3], ist[3], nam[3];
	//        } datim7type;

	typedef struct mavgrecord
	{
		double pavgm[11];
	} mavgtype;

	typedef struct msisrecord
	{
		gts3ctype   gts3c;
		mesotype    meso;
		lowertype   lower;
		parmtype    parm;
		datimetype  datime;
		cswtype     csw;
		mavgtype    mavg;
		dmixtype    dmix;
		parmbtype   parmb;
		metseltype  metsel;
	} msistype;

	//
	//      this common does not get used in the fortran code.
	//      gb and rout are nebver used at all.
	//      tinfg is assigned in globe7, but the value is passed back through the
	//          function value.
	//      tt is set in glob7s, but again, it is passed back through the function value.
	//      typedef struct ttestrecord
	//        {
	//          double tinfg, gb, rout, t[16];
	//        } ttesttype;
	//


	namespace MSIS_Vers {

//		public ref class MSISClass
//		{

			// ---------------- MSIS-com ----------------
			// ------- common to all MSIS models --------

			void tselec
			(
			cswtype&, int[]
			);

			double dnet
				(
				double, double, double, double, double
				);

			double  ccor
				(
				double, double, double, double
				);

			double  ccor2
				(
				double, double, double, double, double
				);

			void meters
				(
				metseltype&
				);

			// ------- common to MSIS-86 and MSIS-00 models --------

			double g0
				(double a, double p[151]);

			double sumex
				(double ex);

			double sg0
				(double ex, double p[151], double ap[8]);

			// ------- common to MSIS-90 and MSIS-00 models --------

			double zeta
				(
				double zz, double zl, double re
				);

			double densu
				(
				parmbtype&, lsqvtype&,
				double, double, double, double, double, double, double&, double,
				double, int, double[], double[], double[]
				);

			double densm
				(
				parmbtype&, fittype&, lsqvtype&,
				double, double, double, double&, int, double[], double[], double[], int,
				double[], double[], double[]
				);

			void spline
				(
				double[], double[], int, double, double, double[]
				);

			void splint
				(
				double[], double[], double[], int, double, double&
				);

			void splini
				(
				double[], double[], double[], int, double, double&
				);

			void glatf
				(
				double, double&, double&
				);

			double  vtst
				(
				cswtype&,
				int, double, double, double, double, double, double,
				double ap[], int
				);

			// ---------------- MSIS-86 model ----------------
			void gts5
				(
				msistype& msis86r,
				lpolytype& lpoly,
				lsqvtype& lsqv,
				int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
				double f107, double ap[8], int mass, double d[10], double t[3]
				);

			void msis86init
				(
				msistype& msis86r
				);

			// ---------------- MSIS-90 model ----------------
			void gtd6
				(
				msistype& msis90r,
				lpolytype& lpoly,
				fittype& fit,
				lsqvtype& lsqv,
				int iyd, double sec, double alt, double glat, double glong, double stl,
				double f107a, double f107, double ap[8], int mass,
				double d[9], double t[3]
				);

			void gts6
				(
				msistype& msis90r,
				lpolytype& lpoly,
				lsqvtype& lsqv,
				int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
				double f107, double ap[8], int mass, double d[10], double t[3]
				);

			void msis90init
				(
				msistype& msis90r
				);

			// -------------- NRLMSISE-00 model --------------
			void gtd7
				(
				msistype& msis00r,
				lpolytype& lpoly,
				fittype& fit,
				lsqvtype& lsqv,
				int iyd, double sec, double alt, double glat, double glong, double stl,
				double f107a, double f107, double ap[8], int mass,
				double d[10], double t[3]
				);

			void gtd7d
				(
				msistype& msis00r,
				lpolytype& lpoly,
				fittype& fit,
				lsqvtype& lsqv,
				int iyd, double sec, double alt, double glat, double glong,
				double stl, double f107a, double f107, double ap[8], int mass,
				double d[10], double t[3]
				);

			void msis00init
				(
				msistype& msis00r
				);

//		};  // MSISClass

	}  // namespace MSISFuncs


