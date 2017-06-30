#ifndef _coordFK5_h_
#define _coordFK5_h_
/*      ----------------------------------------------------------------
*
*                              coordFK5.h
*
*  this file contains routines for iau-76/fk5 conversions.
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*              22 may 09  david vallado
*                           add all transformation
*              21 jan 08  david vallado
*                           fix matrix operations
*              30 may 07  david vallado
*                           3rd edition baseline
*               1 dec 05  david vallado
*                           add frame bias, restructure
*              24 mar 05  david vallado
*                           conversion to c++
*              21 feb 05  david vallado
*                           original baseline
*     *****************************************************************       */
//#include "stdafx.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "astMath.h"
#include "astTime.h"
#include "ast2Body.h"
//#include "EopSpw.h"

#define pi 3.14159265358979323846


#pragma once

//using namespace System;


/*    *****************************************************************
*     type definitions
*     ***************************************************************** */

typedef enum eOpt
  {
    e80,
    e96,
    e00a,
    e00b
  } eOpt;


typedef struct iau80data
  {
     int    iar80[107][6];
     double rar80[107][5];
  } iau80data;


namespace coordFK5 {

	//	public ref class coordFK5Cl
	//	{


	/*    *************************** routines ********************************** */
	void sethelp
		(
		char& iauhelp, char iauopt
		);

	void ddpsiddeps_dxdy
		(
		double ttt, double& ddpsi, double& ddeps,
		edirection direct,
		double& dx, double& dy
		);

	void iau76fk5_itrf_gcrf
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		iau80data& iau80rec,
		char   nutopt,
		double ttt, double jdut1, double lod, double xp,
		double yp, int eqeterms, double ddpsi, double ddeps,
		double deltapsi, double deltaeps,
		std::vector< std::vector<double> > trans
		);

	void iau76fk5all_itrf_gcrf
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		double rpef[3], double vpef[3], double apef[3],
		double rtod[3], double vtod[3], double atod[3],
		double rmod[3], double vmod[3], double amod[3],
		iau80data& iau80rec,
		char   nutopt,
		double ttt, double jdut1, double lod, double xp,
		double yp, int eqeterms, double ddpsi, double ddeps,
		double deltapsi, double deltaeps,
		std::vector< std::vector<double> > trans
		);

	void iau76fk5_itrf_j2k
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rj2k[3], double vj2k[3], double aj2k[3],
		iau80data& iau80rec,
		double ttt, double jdut1, double lod,
		double xp, double yp, int eqeterms,
		std::vector< std::vector<double> > &trans
		);

	void iau76fk5_itrf_mod
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rmod[3], double vmod[3], double amod[3],
		iau80data& iau80rec,
		double ttt, double jdut1, double lod, double xp,
		double yp, int eqeterms, double ddpsi, double ddeps,
		std::vector< std::vector<double> > &trans
		);

	void iau76fk5_itrf_teme
		(
		double ritrf[3], double vitrf[3], double aitrf[3],
		edirection direct,
		double rteme[3], double vteme[3], double ateme[3],
		double ttt, double xp, double yp,
		double jdut1, double lod,
		std::vector< std::vector<double> > &trans
		);

	void pef_gcrf
		(
		double rpef[3], double vpef[3], double apef[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		iau80data& iau80rec,
		double ttt, double jdut1, double lod, int eqeterms,
		double ddpsi, double ddeps
		);

	void tod_gcrf
		(
		double rtod[3], double vtod[3], double atod[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		iau80data& iau80rec,
		double ttt, double ddpsi, double ddeps
		);

	void mod_gcrf
		(
		double rmod[3], double vmod[3], double amod[3],
		edirection direct,
		double rgcrf[3], double vgcrf[3], double agcrf[3],
		double ttt
		);

	//void teme_j2k
	//	(
	//	double rteme[3], double vteme[3], double ateme[3],
	//	edirection direct,
	//	double rj2k[3], double vj2k[3], double aj2k[3],
	//	iau80data& iau80rec,
	//	double ttt, int order, int eqeterms, char optteme
	//	);

	void teme_pef
		(
		double rteme[3], double vteme[3], double ateme[3],
		edirection direct,
		double rpef[3], double vpef[3], double apef[3],
		double jdut1
		);

	void iau80in
		(
		iau80data& iau80rec,
		char EopLoc[85]
		);

	void fundarg
		(
		double ttt, eOpt opt,
		double& l, double& l1, double& f, double& d, double& omega,
		double& lonmer, double& lonven, double& lonear, double& lonmar,
		double& lonjup, double& lonsat, double& lonurn, double& lonnep,
		double& precrate
		);

	void precess
		(
		double ttt, eOpt opt,
		double& psia, double& wa, double& epsa, double& chia,
		std::vector< std::vector<double> > &prec
		);

	void nutation
		(
		double ttt, double ddpsi, double ddeps,
		iau80data& iau80rec, char nutopt,
		double& deltapsi, double& deltaeps, double& trueeps, double& meaneps, double& omega,
		std::vector< std::vector<double> > &nut
		);

	void sidereal
		(
		double jdut1, double deltapsi, double meaneps, double omega,
		double lod, int eqeterms,
		std::vector< std::vector<double> > &st,
		std::vector< std::vector<double> > &stdot
		);

	void polarm
		(
		double xp, double yp, double ttt, eOpt opt, std::vector< std::vector<double> > &pm
		);

	void framebias
		(
		char opt,
		double& term1, double& term2, double& term3, std::vector< std::vector<double> > &fb
		);

	void truemean
		(
		double ttt, int order, int eqeterms, char opt,
		iau80data& iau80rec,
		std::vector< std::vector<double> > &nutteme
		);

	//	};  // Class

}  // namespace coordFK5

#endif
