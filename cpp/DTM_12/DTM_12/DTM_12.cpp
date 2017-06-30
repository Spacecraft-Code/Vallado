/*        ----------------------------------------------------------------
 * 
 *                               DTM_12.cpp
 * 
 *   this file contains common routines for the DTM-2012 atmospheric model.
 *   note that the indices are increased by one so they will appear the same as
 *   original fortran-90 source.
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
/*  ---------------------------------------------------------------------- */

//#include "stdafx.h"

#include "DTM_12.h"

/*  -------------------------- local rountines --------------------------- */

namespace DTM_12 {

 //---------------------------------------------------------------------------
 //
 // ROUTINE: dtm2012
 //
 //> @author Sean Bruinsma
 //>
 //> @brief  calculation of temperature and density with DTM2012
 //
 // PROTOTYPE:
 //                     dtm2012 (
 //                              double intent(in)   day
 //                              double intent(in)   f[2]
 //                              double intent(in)   fbar[2]
 //                              double intent(in)   akp[4]
 //                              double intent(in)   hellp
 //                              double intent(in)   hl
 //                              double intent(in)   alat
 //                              double intent(in)   lon
 //                              double intent(out)   tz
 //                              double intent(out)   tinf
 //                              double intent(out)   tp120
 //                              double intent(out)   ro
 //                              double intent(out)   d[6]
 //                              double intent(out)   wmm
 //                             )
 //
 //
 // INPUT ARGUMENTS:
 //>                    @param[in] day day-of-year [1-366]
 //>                    @param[in] f[2] f[1] = instantaneous flux at (day - 1) / f[2] = 0
 //>                    @param[in] fbar[2] fbar[1] = mean flux at t / fbar[2] = 0.
 //>                    @param[in] akp[4] kp delayed by 3 hours, akp[3] = mean of last 24 hours / akp[2] & akp[4] = 0
 //>                    @param[in] hellp altitude (in km) greater than 120 km
 //>                    @param[in] hl local time (in radian)
 //>                    @param[in] alat latitude (in radian)
 //>                    @param[in] lon longitude (in radian)
 //
 //
 // OUTPUT ARGUMENTS:
 //>                    @param[out] ro           density (g/cm^3) at the given position
 //>                    @param[out] tinf         exosphere temperature
 //>                    @param[out] tz           temperature at the given height
 //>                    @param[out] wmm          mean molecular mass
 //>                    @param[out] tp120	temperature gradient at 120 km
 //>                    @param[out] d[6]         partial density atomic hydrogen [1]
 //>                                             partial density helium [2]
 //>                                             partial density atomic oxygen [3]
 //>                                             partial density molecular nitrogen [4]
 //>                                             partial density molecular oxygen [5]
 //>                                             partial density atomic nitrogen [6] (unused)
 //
 //> @date 06/2012
 //
 //---------------------------------------------------------------------------
 // * aut SB         Mod PS     
 // * ver 06/10/2009   06/2012
 //
 // * rol calculation of temperature and density with DTM2012  
 // * par * *  INPUT * * 
 //     day      = day-of-year [1-366]
 //     f        = f[1] = instantaneous flux at (day - 1) / f[2] = 0.
 //     fbar     = fbar[1] = average flux at t  / fbar[2] = 0.
 //     akp      = akp[1] =  kp delayed by 3 hours, akp[3] = mean of last 24 hours / akp[2] & akp[4] = 0.
 //     hellp     = altitude (in km) greater than 120 km
 //     hl       = local time (in radian)
 //     alat     = latitude (in radian)
 //     lon     = longitude (in radian)
 // * par * *  OUTPUT * * 
 //     tz        = temperature at altitude -> hellp
 //     tinf      = exosphere temperature
 //     d[1]      = partial density atomic hydrogen
 //     d[2]      = partial density helium
 //     d[3]      = partial density atomic oxygen
 //     d[4]      = partial density molecular nitrogen
 //     d[5]      = partial density molecular oxygen
 //     d[6]      = partial density atomic nitrogen
 //     ro        = density (total) in g/cm3
 //     wmm       = mean molecular mass in g
 // cdav increment indicies by one for compatibility to fortran code 
/*  ---------------------------------------------------------------------- */

void dtm2012
    (
	  plgdtmtype& plgdtm, hlocaltype& hlocal, eclipttype& eclipt, datmotype& datmo, constype& cons, pardtmtype& pardtm,
	  double day, double f[3], double fbar[3], double akp[5], double hellp, double hl, double alat, double lon,
      double& tz, double& tinf, double& tp120, double& ro, double d[7], double& wmm
	)
    {
     //.. Local Scalars .. 
     int i, ialt;
     int ityp = 0;
     double c,c2,c4,cmg,dt120,dtinf, 
            dtp120,gamma,gdelaz,gdelh,gdelo,gdelt0,gdeltp,s2,tinftz,zeta;
     double cpmg = 0.19081, re = 6356.77, rgas = 831.4;
 
     double clmlmg,cmg2,cmg4,dzeta, 
            dzeta2,expsz,gdelaz2,gdelhe,gdelo2,gdelt,glb,s,sigma,sigzeta, 
            sp,t120tz,upapg;

     double xlmg = -1.2392;
     double t120 = 0.0, xlog = 0.0;
     double cose = 0.9175, gsurf = 980.665, sine = 0.3978, spmg = 0.98163, zero = 0.0;
     double zlb;
     double zlb0 = 120.0;
   
     //.. Local Arrays 
     double cc[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
     double dbase[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
     double fz[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

     //.. Data Declarations .. 
	 // cdav increment all arrays so they can begin at 1
     static double alefa[7] = {0.0, -0.40, -0.38, 0.0, 0.0, 0.0, 0.0};
     static int    ma[7]    = {0, 1, 4, 16, 28, 32, 14};
     static double vma[7]   = {0.0, 1.6606e-24, 6.6423e-24, 26.569e-24, 46.4958e-24, 53.1381e-24, 23.2479e-24};

     // ... Executable Statements ...
     zlb = zlb0;

     ialt = int(hellp);
     ro = 0.0;

     dtinf = 0.0;
     dt120 = 0.0;
     dtp120 = 0.0;
     fz[1] = 0.0;
     fz[2] = 0.0;
     fz[3] = 0.0;
     fz[4] = 0.0;
     fz[5] = 0.0;
     fz[6] = 0.0;
     //
	 pardtm.dtt.T_one  = 0.0;
     pardtm.dh.T_one   = 0.0;
     pardtm.dhe.T_one  = 0.0;
     pardtm.dox.T_one  = 0.0;
     pardtm.daz2.T_one = 0.0;
     pardtm.do2.T_one  = 0.0;
     pardtm.daz.T_one  = 0.0;
     pardtm.dt0.T_one  = 0.0;
     pardtm.dtp.T_one  = 0.0;
  
     //   calcul des polynomes de legendre
     c   = sin(alat);
     c2  = c * c;
     c4  = c2 * c2;
     s   = cos(alat);
     s2  = s * s;
     plgdtm.p10 = c;
     plgdtm.p20 = 1.5 * c2 - 0.5;
     plgdtm.p30 = c * (2.5 * c2 - 1.5);
     plgdtm.p40 = 4.375 * c4 - 3.75 * c2 + 0.375;
     plgdtm.p50 = c * (7.875 * c4 - 8.75 * c2 + 1.875);
     plgdtm.p60 = (5.5 * c * plgdtm.p50 - 2.5 * plgdtm.p40) / 3.0;
     plgdtm.p11 = s;
     plgdtm.p21 = 3.0 * c * s;
     plgdtm.p31 = s * (7.5 * c2 - 1.5);
     plgdtm.p41 = c * s * (17.5 * c2 - 7.5);
     plgdtm.p51 = s * (39.375 * c4 - 26.25 * c2 + 1.875);
     plgdtm.p22 = 3.0 * s2;
     plgdtm.p32 = 15.0 * c * s2;
     plgdtm.p42 = s2 * (52.5 * c2 - 7.5);
     plgdtm.p52 = 3.0 * c * plgdtm.p42 - 2.0 * plgdtm.p32;
     plgdtm.p62 = 2.75 * c * plgdtm.p52 - 1.75 * plgdtm.p42;
     plgdtm.p33 = 15.0 * s * s2;

     //   calcul des polynomes de legendre / pole magnetique (79n,71w)
     clmlmg = cos(lon - xlmg);
     sp     = s * cpmg * clmlmg + c * spmg;
     
     cmg    = sp; // pole magnetique
     cmg2   = cmg * cmg;
     cmg4   = cmg2 * cmg2;
     plgdtm.p10mg  = cmg;
     plgdtm.p20mg  = 1.5 * cmg2 - 0.5;
     plgdtm.p40mg  = 4.375 * cmg4 - 3.75 * cmg2 + 0.375;
     
     //   heure locale 
     hlocal.hl0 = hl;
     hlocal.ch  = cos(hlocal.hl0);
     hlocal.sh  = sin(hlocal.hl0);
     hlocal.c2h = hlocal.ch * hlocal.ch - hlocal.sh * hlocal.sh;
     hlocal.s2h = 2.0 * hlocal.ch * hlocal.sh;
     hlocal.c3h = hlocal.c2h * hlocal.ch - hlocal.s2h * hlocal.sh;
     hlocal.s3h = hlocal.s2h * hlocal.ch + hlocal.c2h * hlocal.sh;
     
     //   calcul of fonction g(l) / tinf, t120, tp120
     //// new version declaration of tt,dtt en struct dtm
     gldtm_XX(plgdtm, hlocal, datmo, f, fbar, akp, day, pardtm.tt, pardtm.dtt, gdelt, 1.0, lon);
     pardtm.dtt.T_one = 1.0 + gdelt;
     tinf = pardtm.tt.T_one * pardtm.dtt.T_one;

     gldtm_XX(plgdtm, hlocal, datmo, f, fbar, akp, day, pardtm.t0, pardtm.dt0, gdelt0, 1.0, lon);
     pardtm.dt0.T_one = 1.0 + gdelt0;
     t120 = pardtm.t0.T_one * pardtm.dt0.T_one;

     gldtm_XX(plgdtm, hlocal, datmo, f, fbar, akp, day, pardtm.tp, pardtm.dtp, gdeltp, 1.0, lon);
     pardtm.dtp.T_one = 1.0 + gdeltp;
     tp120 = pardtm.tp.T_one * pardtm.dtp.T_one;

     //-----------------------------------------------------------------------------
     //   calcul of concentrations n(z): H, HE, O, N2, O2, N
     sigma = tp120 / (tinf - t120);
     dzeta = (re + zlb) / (re + hellp);
     zeta  = (hellp - zlb) * dzeta;
     dzeta2 = dzeta * dzeta;
     sigzeta =  sigma * zeta;
     expsz  = exp(-sigzeta);
     tz    = tinf - (tinf - t120) * expsz;

     gldtm_XX(plgdtm, hlocal, datmo, f, fbar, akp, day, pardtm.h, pardtm.dh, gdelh, 0.0, lon);
     pardtm.dh.T_one = exp(gdelh);
     dbase[1] = pardtm.h.T_one * pardtm.dh.T_one;

     gldtm_XX(plgdtm, hlocal, datmo, f, fbar, akp, day, pardtm.he, pardtm.dhe, gdelhe, 0.0, lon);
     pardtm.dhe.T_one = exp(gdelhe);
     dbase[2] =  pardtm.he.T_one * pardtm.dhe.T_one;

     gldtm_XX(plgdtm, hlocal, datmo, f, fbar, akp, day, pardtm.ox, pardtm.dox, gdelo, 1.0, lon);
     pardtm.dox.T_one = exp(gdelo);
     dbase[3] = pardtm.ox.T_one * pardtm.dox.T_one;
 
     gldtm_XX(plgdtm, hlocal, datmo, f, fbar, akp, day, pardtm.az2, pardtm.daz2, gdelaz2, 1.0, lon);
     pardtm.daz2.T_one = exp(gdelaz2);
     dbase[4] = pardtm.az2.T_one * pardtm.daz2.T_one;

     gldtm_XX(plgdtm, hlocal, datmo, f, fbar, akp, day, pardtm.o2, pardtm.do2, gdelo2, 1.0, lon);
     pardtm.do2.T_one = exp(gdelo2);
     dbase[5] =  pardtm.o2.T_one * pardtm.do2.T_one;

     gldtm_XX(plgdtm, hlocal, datmo, f, fbar, akp, day, pardtm.az, pardtm.daz, gdelaz, 1.0, lon);
     pardtm.daz.T_one = exp(gdelaz);
     dbase[6] = pardtm.az.T_one * pardtm.daz.T_one;
     //
     glb = gsurf / pow(1.0 + zlb/re, 2);
     glb = glb / (sigma * rgas * tinf);
     t120tz = t120 / tz;
     tinftz = tinf / tz;
     for (i = 1; i <= 6; i++)
	 {
         gamma = ma[i] * glb;
         upapg = 1.0 + alefa[i] + gamma;
         fz[i] = pow( t120tz, upapg ) * exp(-sigzeta * gamma);
         //   concentrations en H, HE, O, N2, O2, N
         cc[i] = dbase[i] * fz[i];
         //   densites en H, HE, O, N2, O2, N
         d[i] = cc[i] * vma[i];
         //   densite totale
         ro = ro + d[i];
	 }

     //  average of atomic mass                              
     wmm = ro/(vma[1] * (cc[1] + cc[2] + cc[3] + cc[4] + cc[5] + cc[6]));
    } // function dtm2012


//---------------------------------------------------------------------------
//
// ROUTINE: densisty_uncertainty
//
//> @author Sean Bruinsma
//>
//> @brief  calculation of density uncertainty
//
// PROTOTYPE:
//        density_uncertainty (
//                              double intent(in)   alt
//                              double intent(in)   latgd
//                              double intent(in)   lst
//                              double intent(in)   flux
//                              double intent(in)   kp
//                              double intent(out)   unc
//                             )
//
//
// INPUT ARGUMENTS:
//>                    @param[in] alt altitude in km (NB: uncertainties NOT correct - too big - for altitudes<300 km)
//>                    @param[in] latgd latitude in degrees
//>                    @param[in] lst local solar time in hours (0. - 23.999)
//>                    @param[in] flux mean solar flux (sfu)
//>                    @param[in] kp geomagnetic activity (0. - 9.0)
//
//
// OUTPUT ARGUMENTS:
//>                    @param[out] unc uncertainty (1-sigma) due to density variations below model resolution, in _ of density
//
//> @date
//
//
// alt  = altitude in km (NB: uncertainties NOT correct - too big - for altitudes<300 km)
// latgd  = latitude in degrees
// lst  = local solar time in hours (0. - 23.999)
// flux = mean solar flux (sfu)
// kp   = geomagnetic activity (0. - 9.0)
// unc  = uncertainty (1-sigma) due to density variations below model resolution, in _ of density
//
//NB: based on relative density variations with scales < 2400 km, at 350 and 500 km altitude only,
//    below and above that altitude uncertainties are extrapolated or constants
//
//  Raul Dominguez: Modification for ATMOP DTM2012 wrapper:
//                  Modified so the input files are read only once
//
//  Raul Dominguez: Modification for ATMOP DTM2012 wrapper:
//                  modified so arguments alt and flux are not changed within the routine
//                  I see no reason for such behaviour
//
//  Raul Dominguez: Modified to allow setting file units in a external module
//                  This improves portability
//---------------------------------------------------------------------------

void density_uncertainty
   ( 
	 double nomgrid[25][19], double alt_scale[25][19], double kp_scale[25][19], 
	 double alt, double latgd, double lst, double flux, double kp, double& unc
   )
   {
     double steplat, unc_nom, unc_alt, unc_kp, unc_f, alti_factor, kp_factor, flux_factor, flux2;
	 double rad2deg;
     int  ind_lat, ind_lst;

     //radg modification
     double  alt_aux, flux_aux;
// cdav move definitions to test

	 rad2deg = 57.2957795131;
     //radg modification
     alt_aux = alt;
     flux_aux = flux;

     //radg modification to allow setting file units and names externally
     unc = -999.0;

// cdav move if firstrun section to main
     steplat = 10.0;
	 // cdav change latgd to deg for the index calculation * rad2deg
     ind_lat = int((latgd * rad2deg + 90.0)/steplat) + 1;
     ind_lst = int(lst) + 1;
     // ---
     unc_nom = nomgrid[ind_lst][ind_lat];     //nominal relative variability: @350 km, flux = 65, kp<3

     if(alt_aux > 650.0) 
 		 alt_aux = 650.0;                     //in fact no information above 500 km
     alti_factor = (alt_aux - 350.0)/150.0;   //altitude effect: larger with increasing altitude
     if(alt_aux < 350.0) 
		 alti_factor = 0.0;                   //no information yet; will be computed with GOCE
     unc_alt = unc_nom * (1.0 + alti_factor * alt_scale[ind_lst][ind_lat]);

     kp_factor = (kp - 4.0)/2.0;              //geomagnetic activity effect
     if(kp < 4.0) 
 		 kp_factor = 0.0;
     unc_kp = unc_alt * (1.0 + kp_factor * kp_scale[ind_lst][ind_lat]);

     if(flux_aux > 200.0) 
	 	 flux_aux = 200.0;
     flux2 = flux_aux * flux_aux;             //81-day mean solar activity effect
     flux_factor = 2.5506 - 3.7403e-2 * flux_aux + 2.4455e-4 * flux2 - (5.5582e-7 * flux_aux * flux2);
     unc_f = flux_factor * unc_kp;

     // --- uncertainty in _
     unc  = unc_f * 100.0;
    } // function density_uncertainty



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
// SB 06/03/2009 MOD PS 04/2012
// Corrections & new variables
// * rol calculation of function g(l)  for dtm2009 & dtm2012
//
//     DTM =  table of coefficients for g(l) or for each component
//     ddtm = table of derivatives dg(l)/da
//     ff0 = 1 : for oxygen, nitrogen, helium, temperature
//     ff0 = 0 : hydrogen
//
//     gdel = result for g(l)
//     lon = longitude (not used) 
//       double intent(in)   day,ff0,lon;
//       double intent(out)   gdel;
//       double dimension[2], intent(in)   f,fbar;
//       double dimension[4], intent(in)   akp;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

void gldtm_XX
	(
	  plgdtmtype& plgdtm, hlocaltype& hlocal, datmotype& datmo, 
	  double f[3], double fbar[3], double akp[5], double day, DTM_12type& DTM_12, ddtm_12type& ddtm_12, 
	  double& gdel, double ff0, double lon
    )
	{
     int   kle_eq = 0;
     //.. Local Scalars .. 
     int   i;
     int   ikp = 0;
     double  f1f,fp;
     double  c2fi = 0.0;
     double  clfl,cos2te,coste,dakp,dakpm,dkp,dkpm,f0,slfl;
     double  rot = 0.017214206,
		     rot2 = 0.034428412,
			 roth = 0.261799387, 
             rots = 7.27220e-05;
     // 
     //.. Local Arrays .. 
	 // cdav increment array size by 1 so indicies match fortran
     double fbm150[3],fmfb[3];

	 fbm150[1] = 0.0;
	 fbm150[2] = 0.0;
	 fmfb[1] = 0.0;
	 fmfb[2] = 0.0;

     // ... Executable Statements ...
     //   termes in latitude
     ddtm_12.T_lat[1] = plgdtm.p20;
     ddtm_12.T_lat[2] = plgdtm.p40;
     ddtm_12.T_lat[3] = plgdtm.p10;
     ddtm_12.T_lat[4] = plgdtm.p30;
     ddtm_12.T_lat[5] = plgdtm.p50;
     ddtm_12.T_lat[6] = plgdtm.p60;

     //   termes of flux
     fmfb[1] = f[1] - fbar[1];
     fmfb[2] = f[2] - fbar[2];
     fbm150[1] = fbar[1] - 150.0;
     fbm150[2] = fbar[2];

     ddtm_12.T_flux[1] = fmfb[1];
     ddtm_12.T_flux[3] = fbm150[1];
       
     ddtm_12.T_flux[1] = ddtm_12.T_flux[1] + DTM_12.T_flux[5] * fmfb[2];
     ddtm_12.T_flux[3] = ddtm_12.T_flux[3] + DTM_12.T_flux[6] * fbm150[2];
    
     ddtm_12.T_flux[5] = fmfb[2] * (DTM_12.T_flux[1] + 2.0 * DTM_12.T_flux[2] * ddtm_12.T_flux[1] 
                         + DTM_12.T_flux[7] * plgdtm.p10 + DTM_12.T_flux[8] * plgdtm.p20 + DTM_12.T_flux[9] * plgdtm.p30);

     ddtm_12.T_flux[6] = fbm150[2] * (DTM_12.T_flux[3] + 2.0 * DTM_12.T_flux[4] * ddtm_12.T_flux[3] 
                         + DTM_12.T_flux[10] * plgdtm.p10 + DTM_12.T_flux[11] * plgdtm.p20 + DTM_12.T_flux[12] * plgdtm.p30);

     ddtm_12.T_flux[2]  = ddtm_12.T_flux[1] * ddtm_12.T_flux[1];
     ddtm_12.T_flux[4]  = ddtm_12.T_flux[3] * ddtm_12.T_flux[3];
     ddtm_12.T_flux[7]  = ddtm_12.T_flux[1] * plgdtm.p10;
     ddtm_12.T_flux[8]  = ddtm_12.T_flux[1] * plgdtm.p20;
     ddtm_12.T_flux[9]  = ddtm_12.T_flux[1] * plgdtm.p30;
     ddtm_12.T_flux[10] = ddtm_12.T_flux[3] * plgdtm.p20;
     ddtm_12.T_flux[11] = ddtm_12.T_flux[3] * plgdtm.p30;
     ddtm_12.T_flux[12] = ddtm_12.T_flux[3] * plgdtm.p40;

     //   termes in kp
     c2fi = 1.0 - plgdtm.p10mg * plgdtm.p10mg;

     dkp = akp[1] + (DTM_12.T_kp[5] + c2fi * DTM_12.T_kp[6]) * akp[2];
     dakp = DTM_12.T_kp[1] + DTM_12.T_kp[2] * plgdtm.p20mg + DTM_12.T_kp[11] * plgdtm.p40mg + 
            2.0 * dkp * (DTM_12.T_kp[3] + DTM_12.T_kp[4] * plgdtm.p20mg + DTM_12.T_kp[14] * 2.0 * dkp * dkp);

     ddtm_12.T_kp[5] = dakp * akp[2];
     ddtm_12.T_kp[6] =  ddtm_12.T_kp[5] * c2fi;

     dkpm = akp[3] + DTM_12.T_kp[10] * akp[4];

     dakpm = DTM_12.T_kp[7] + DTM_12.T_kp[8] * plgdtm.p20mg + DTM_12.T_kp[12] * plgdtm.p40mg + 
             2.0 * dkpm * (DTM_12.T_kp[9] + DTM_12.T_kp[13] * plgdtm.p20mg + 
             DTM_12.T_kp[15] * 2.0 * dkpm * dkpm);

     ddtm_12.T_kp[10] = dakpm * akp[4];

     ddtm_12.T_kp[1]  = dkp;
     ddtm_12.T_kp[2]  = plgdtm.p20mg * dkp;
     ddtm_12.T_kp[11] = plgdtm.p40mg * dkp;
     ddtm_12.T_kp[3]  = dkp * dkp;
     ddtm_12.T_kp[4]  = plgdtm.p20mg * ddtm_12.T_kp[3];
     ddtm_12.T_kp[14] = ddtm_12.T_kp[3] * ddtm_12.T_kp[3];
     ddtm_12.T_kp[7]  = dkpm;
     ddtm_12.T_kp[8]  = plgdtm.p20mg * dkpm;
     ddtm_12.T_kp[12] = plgdtm.p40mg * dkpm;
     ddtm_12.T_kp[9]  = dkpm * dkpm;
     ddtm_12.T_kp[13] = plgdtm.p20mg * ddtm_12.T_kp[9];
     ddtm_12.T_kp[15] = ddtm_12.T_kp[9] * ddtm_12.T_kp[9];

     //   function: g(l) non periodic
     f0 = DTM_12.T_flux[1] * ddtm_12.T_flux[1] + DTM_12.T_flux[2] * ddtm_12.T_flux[2]  
         + DTM_12.T_flux[3] * ddtm_12.T_flux[3] + DTM_12.T_flux[4] * ddtm_12.T_flux[4] 
         + DTM_12.T_flux[7] * ddtm_12.T_flux[7] + DTM_12.T_flux[8] * ddtm_12.T_flux[8]  
         + DTM_12.T_flux[9] * ddtm_12.T_flux[9]  + DTM_12.T_flux[10] * ddtm_12.T_flux[10]  
         + DTM_12.T_flux[11] * ddtm_12.T_flux[11] + DTM_12.T_flux[12] * ddtm_12.T_flux[12];

     f1f = 1.0 + f0 * ff0;
     //     write(6, * ) ' ff0,f0,f1f ',ff0,f0,f1f

     f0 = f0 + DTM_12.T_lat[1] * ddtm_12.T_lat[1] + DTM_12.T_lat[2] * ddtm_12.T_lat[2] 
            + DTM_12.T_lat[3] * ddtm_12.T_lat[3] + DTM_12.T_lat[4] * ddtm_12.T_lat[4] 
            + DTM_12.T_kp[1] * ddtm_12.T_kp[1] + DTM_12.T_kp[2] * ddtm_12.T_kp[2] 
            + DTM_12.T_kp[3] * ddtm_12.T_kp[3] + DTM_12.T_kp[4] * ddtm_12.T_kp[4] 
            + DTM_12.T_kp[11] * ddtm_12.T_kp[11] + DTM_12.T_kp[7] * ddtm_12.T_kp[7] 
            + DTM_12.T_kp[8] * ddtm_12.T_kp[8] + DTM_12.T_kp[9] * ddtm_12.T_kp[9] 
            + DTM_12.T_kp[12] * ddtm_12.T_kp[12] + DTM_12.T_kp[13] * ddtm_12.T_kp[13] 
            + DTM_12.T_kp[14] * ddtm_12.T_kp[14] + DTM_12.T_kp[15]  * ddtm_12.T_kp[15]  
            + DTM_12.T_lat[5] * ddtm_12.T_lat[5] + DTM_12.T_lat[6] * ddtm_12.T_lat[6];

     //   terms annual & symetric in latitude 
     ddtm_12.T_SLat[1] = cos(rot * (day - DTM_12.T_dPhas[1]));
     ddtm_12.T_SLat[2] = plgdtm.p20 * ddtm_12.T_SLat[1];

     //   terms semi-annual & symetric in latitude
     ddtm_12.T_SASLat[1] = cos(rot2 * (day - DTM_12.T_dPhas[2]));
     ddtm_12.T_SASLat[2] = plgdtm.p20 * ddtm_12.T_SASLat[1];

     //   terms annual & non-symetric in latitude 
     coste = cos(rot * (day - DTM_12.T_dPhas[3]));
     ddtm_12.T_NSLat[1] = plgdtm.p10 * coste;
     ddtm_12.T_NSLat[2] = plgdtm.p30 * coste;
     ddtm_12.T_NSLat[3] = ddtm_12.T_flux[3] * ddtm_12.T_NSLat[1];
     //      ddtm_12.T_NSLat[3] = p50 * coste

     //   terms semi-annual  non-symetric in latitude 
     cos2te = cos(rot2 * (day - DTM_12.T_dPhas[4]));
     ddtm_12.T_SANSLat[1] = plgdtm.p10 * cos2te;
     ddtm_12.T_SANSLat[2] = plgdtm.p30 * cos2te;
     ddtm_12.T_SANSLat[3] = ddtm_12.T_flux[3] * ddtm_12.T_SANSLat[1];
     //      ddtm_12.T_SANSLat[3] = p50 * cos2te

     //   terms diurnal (& annual coupled)
     ddtm_12.T_DiAn[1]  = plgdtm.p11 * hlocal.ch;
     ddtm_12.T_DiAn[2]  = plgdtm.p31 * hlocal.ch;
     ddtm_12.T_DiAn[3]  = ddtm_12.T_flux[3] * ddtm_12.T_DiAn[1] ;
     ddtm_12.T_DiAn[4]  = ddtm_12.T_DiAn[1] * coste;
     ddtm_12.T_DiAn[5]  = plgdtm.p21 * hlocal.ch * coste; 
     ddtm_12.T_DiAn[6]  = plgdtm.p11 * hlocal.sh;
     ddtm_12.T_DiAn[7]  = plgdtm.p31 * hlocal.sh;
     ddtm_12.T_DiAn[8]  = ddtm_12.T_flux[3] * ddtm_12.T_DiAn[6];
     ddtm_12.T_DiAn[9]  = ddtm_12.T_DiAn[6] * coste;
     ddtm_12.T_DiAn[10] = plgdtm.p21 * hlocal.sh * coste;
     ddtm_12.T_DiAn[11] = plgdtm.p51 * hlocal.ch;
     ddtm_12.T_DiAn[12] = plgdtm.p51 * hlocal.sh;  

     //   terms semi-diurnes (& annual coupled)
     ddtm_12.T_SDiAn[1]  = plgdtm.p22 * hlocal.c2h;
     ddtm_12.T_SDiAn[2]  = plgdtm.p42 * hlocal.c2h;
     ddtm_12.T_SDiAn[3]  = plgdtm.p32 * hlocal.c2h * coste;
     ddtm_12.T_SDiAn[4]  = plgdtm.p22 * hlocal.s2h;
     ddtm_12.T_SDiAn[5]  = plgdtm.p42 * hlocal.s2h;
     ddtm_12.T_SDiAn[6]  = plgdtm.p32 * hlocal.s2h * coste;
     ddtm_12.T_SDiAn[7]  = plgdtm.p32 * hlocal.c2h; //coeff. rajoute pour tp120/t120 (slb)
     ddtm_12.T_SDiAn[8]  = plgdtm.p32 * hlocal.s2h;
     ddtm_12.T_SDiAn[9]  = ddtm_12.T_flux[3] * ddtm_12.T_SDiAn[1];
     ddtm_12.T_SDiAn[10] = ddtm_12.T_flux[3] * ddtm_12.T_SDiAn[4];
     ddtm_12.T_SDiAn[11] = plgdtm.p62 * hlocal.c2h;
     ddtm_12.T_SDiAn[12] = plgdtm.p62 * hlocal.s2h;
  
     //   terms ter-diurnes
     ddtm_12.T_TDi[1] = plgdtm.p33 * hlocal.c3h;
     ddtm_12.T_TDi[2] = plgdtm.p33 * hlocal.s3h;
  
     //   function periodic -> g(l) 
     fp = DTM_12.T_SLat[1]    * ddtm_12.T_SLat[1]   + DTM_12.T_SLat[2]   * ddtm_12.T_SLat[2] 
         + DTM_12.T_SASLat[1] * ddtm_12.T_SASLat[1] + DTM_12.T_SASLat[2] * ddtm_12.T_SASLat[2]
         + DTM_12.T_NSLat[1]  * ddtm_12.T_NSLat[1]  + DTM_12.T_NSLat[2]  * ddtm_12.T_NSLat[2]  
         + DTM_12.T_NSLat[3]  * ddtm_12.T_NSLat[3]  + DTM_12.T_SANSLat[1]* ddtm_12.T_SANSLat[1] 
         + DTM_12.T_DiAn[1]   * ddtm_12.T_DiAn[1]   + DTM_12.T_DiAn[2]   * ddtm_12.T_DiAn[2] 
         + DTM_12.T_DiAn[3]   * ddtm_12.T_DiAn[3]   + DTM_12.T_DiAn[4]   * ddtm_12.T_DiAn[4] 
         + DTM_12.T_DiAn[5]   * ddtm_12.T_DiAn[5]   + DTM_12.T_DiAn[6]   * ddtm_12.T_DiAn[6] 
         + DTM_12.T_DiAn[7]   * ddtm_12.T_DiAn[7]   + DTM_12.T_DiAn[8]   * ddtm_12.T_DiAn[8] 
         + DTM_12.T_DiAn[9]   * ddtm_12.T_DiAn[9]   + DTM_12.T_DiAn[10]  * ddtm_12.T_DiAn[10] 
         + DTM_12.T_SDiAn[1]  * ddtm_12.T_SDiAn[1]  + DTM_12.T_SDiAn[3]  * ddtm_12.T_SDiAn[3] 
         + DTM_12.T_SDiAn[4]  * ddtm_12.T_SDiAn[4]  + DTM_12.T_SDiAn[6]  * ddtm_12.T_SDiAn[6] 
         + DTM_12.T_TDi[1]    * ddtm_12.T_TDi[1]    + DTM_12.T_TDi[2]    * ddtm_12.T_TDi[2] 
         + DTM_12.T_SDiAn[2]  * ddtm_12.T_SDiAn[2]  + DTM_12.T_SDiAn[5]  * ddtm_12.T_SDiAn[5] 
         + DTM_12.T_SANSLat[2]* ddtm_12.T_SANSLat[2]+ DTM_12.T_SANSLat[3]* ddtm_12.T_SANSLat[3] 
         + DTM_12.T_SDiAn[7]  * ddtm_12.T_SDiAn[7]  + DTM_12.T_SDiAn[8]  * ddtm_12.T_SDiAn[8]  
         + DTM_12.T_SDiAn[9]  * ddtm_12.T_SDiAn[9]  + DTM_12.T_SDiAn[10] * ddtm_12.T_SDiAn[10] 
         + DTM_12.T_SDiAn[11] * ddtm_12.T_SDiAn[11] + DTM_12.T_SDiAn[12] * ddtm_12.T_SDiAn[12]  
         + DTM_12.T_DiAn[11]  * ddtm_12.T_DiAn[11]  + DTM_12.T_DiAn[12]  * ddtm_12.T_DiAn[12]; 
  
     //   terms magnetic activity 
     ddtm_12.T_AMg[1] = plgdtm.p10 * coste * dkp;
     ddtm_12.T_AMg[2] = plgdtm.p30 * coste * dkp;
     ddtm_12.T_AMg[3] = plgdtm.p50 * coste * dkp;
     ddtm_12.T_AMg[4] = plgdtm.p11 * hlocal.ch * dkp;
     ddtm_12.T_AMg[5] = plgdtm.p31 * hlocal.ch * dkp;
     ddtm_12.T_AMg[6] = plgdtm.p51 * hlocal.ch * dkp;
     ddtm_12.T_AMg[7] = plgdtm.p11 * hlocal.sh * dkp;
     ddtm_12.T_AMg[8] = plgdtm.p31 * hlocal.sh * dkp;
     ddtm_12.T_AMg[9] = plgdtm.p51 * hlocal.sh * dkp;
   
     //   function g(l) (additional periodic)
     fp = fp + DTM_12.T_AMg[1] * ddtm_12.T_AMg[1] + DTM_12.T_AMg[2] * ddtm_12.T_AMg[2] 
             + DTM_12.T_AMg[3] * ddtm_12.T_AMg[3] + DTM_12.T_AMg[4] * ddtm_12.T_AMg[4] 
             + DTM_12.T_AMg[5] * ddtm_12.T_AMg[5] + DTM_12.T_AMg[6] * ddtm_12.T_AMg[6] 
             + DTM_12.T_AMg[7] * ddtm_12.T_AMg[7] + DTM_12.T_AMg[8] * ddtm_12.T_AMg[8] 
             + DTM_12.T_AMg[9] * ddtm_12.T_AMg[9];
     //
     dakp = (DTM_12.T_AMg[1] * plgdtm.p10 + DTM_12.T_AMg[2] * plgdtm.p30 + DTM_12.T_AMg[3] * plgdtm.p50) * coste +
            (DTM_12.T_AMg[4] * plgdtm.p11 + DTM_12.T_AMg[5] * plgdtm.p31 + DTM_12.T_AMg[6] * plgdtm.p51) * hlocal.ch + 
            (DTM_12.T_AMg[7] * plgdtm.p11 + DTM_12.T_AMg[8] * plgdtm.p31 + DTM_12.T_AMg[9] * plgdtm.p51) * hlocal.sh;

     ddtm_12.T_kp[5] = ddtm_12.T_kp[5] + dakp * akp[2];
     ddtm_12.T_kp[6] = ddtm_12.T_kp[5] + dakp * c2fi * akp[2];

     //   terms in longitude
     clfl = cos(lon);
     ddtm_12.T_Lon[1] = plgdtm.p11 * clfl;
     ddtm_12.T_Lon[2] = plgdtm.p21 * clfl;
     ddtm_12.T_Lon[3] = plgdtm.p31 * clfl;
     ddtm_12.T_Lon[4] = plgdtm.p41 * clfl;
     ddtm_12.T_Lon[5] = plgdtm.p51 * clfl;
     slfl = sin(lon);
     ddtm_12.T_Lon[6] = plgdtm.p11 * slfl;
     ddtm_12.T_Lon[7] = plgdtm.p21 * slfl;
     ddtm_12.T_Lon[8] = plgdtm.p31 * slfl;
     ddtm_12.T_Lon[9] = plgdtm.p41 * slfl;
     ddtm_12.T_Lon[10] = plgdtm.p51 * slfl;
   
     //   function g(l) periodic (additional)
     for (i = 1; i <=  10; i++ )  //Nb_Lon cdav Nb_Lon is not recording??
         fp = fp + DTM_12.T_Lon[i] * ddtm_12.T_Lon[i];

     //   function g(l) sum (coupled with flux)
     gdel = f0 + fp * f1f;

    } // function gldtm_XX


	// additional routines

 //---------------------------------------------------------------------------
 //
 // FUNCTION: a2K
 //
 //> @brief   This function gets K from a, based on the "tramkm.f" subroutine
 //>
 //> @authors Sean Bruinsma
 //>          Raul Dominguez
 //
 // INPUT ARGUMENTS:   int a            => "a" index
 //>                   @param[in] a a index value
 //
 // OUTPUT ARGUMENTS:  double a2K           <= "K" index corresponding to "a"
 //
 // REVISION HISTORY:
 //
 //  @version 1.0
 //  @date 09/08/2012
 //  Upgraded from FORTRAN 77 code
 //
 //  @version 1.1
 //  @date 14/09/2011
 //  Change output from double to real
 //
 //  @version 1.2
 //  @date 20/09/2011
 //  Check input "a" value. If it is less than 0, fatal error
 //
 //> @version 1.3
 //> @date 01/10/2011
 //  Changed input type. Now the routine accepts a real number (instead of
 //  an int)
 //---------------------------------------------------------------------------
 
double a2K
	(  double a  )
	{
        //Conversion tables
        int table_length;
        table_length = 28;
		// cdav
		double temp;

        double AM[29] = { 0.0,   0.0,   1.4,   3.4,   5.4,   7.4,  10.4,  13.4,  
                         16.4,  20.4,  26.4,  33.4,  40.4,  50.4,  60.4,  70.4,   
                         86.4, 103.4, 120.4, 146.4, 173.4, 200.4, 243.4, 286.4, 
						330.4, 386.4, 443.4, 500.4, 611.4};

        double AKM[29] = {0.0,   0.0,   0.33,  0.66,  1.0,   1.33,  1.66,  2.0,
			              2.33,  2.66,  3.0,   3.33,  3.66,  4.0,   4.33,  4.66,  
						  5.0,   5.33,  5.66,  6.0,   6.33,  6.66,  7.0,   7.33,  
						  7.66,  8.0,   8.33,  8.66,  9.0};
        //Auxiliaries
        int i, im;

		//if a is less than 0, return an error
        if (a < 0.0) 
    		fatal_error(11);

        temp = -1.0;
        for (i=1; i <=table_length; i++)
		{
              if(a == AM[i]) 
                  return temp = AKM[i];
              else if (a <= AM[i]) 
			  {
				  im = i - 1;
                  return temp = AKM[im] + (AKM[i] - AKM[im]) * (a - AM[im]) / (AM[i] - AM[im]);
			  }
		}  // for

         if (temp < 0.0) 
            fatal_error(12);
         return temp;
	 }  // function a2k


 //---------------------------------------------------------------------------
 //
 // FUNCTION: local_solar_time
 //
 //> @brief  find the local solar time in radians given the universal time and
 //>              the longitude
 //
 //> @details This function computes the local solar time in radians from the
 //>          universal time and the longitude, by means of the formula
 //>         @f$LST =2\cdot \pi ( nsec/3600 + \lambda \cdot (24/3600) )@f$
 //>           this function assumes then that local solar time
 //>           ranges from 0 to 2pi.
 //>           It also sanitizes output, so angles out of the 0 to 2pi range
 //>           are converted to the 0-2pi interval
 //
 //> @author Raul Dominguez
 //
 // INPUT ARGUMENTS:
 //>                    @param[in]  longitude  longitude in degress
 //>                    @param[in]   hour      universal time (hour)
 //>                    @param[in]   minute    universal time (minute)
 //>                    @param[in]   sec       universal time (seconds)
 //>
 // OUTPUT ARGUMENTS:
 //>                    @param[out]  local_solar_time  local solar time in radians
 //
 // REVISION HISTORY:
 //
 // version 1.0
 // date 31/07/2012
 // Initial version
 //
 // version 1.1
 // date 11/09/2012
 // Sanitize output, so that it is between 0 and two pi
 //
 //> @version 1.2
 //> @date 14/09/2012
 //> Change in/out from double to real
 //
 //---------------------------------------------------------------------------
 
double local_solar_time 
	  (
	     double longitude, int hour, int minute, double sec
	  )
 {
   //auxiliaries
   double nsec, temp, pi;

   pi = 3.14159265358979323846;
 
   //compute the number of seconds elapsed since start of the day
   //note that seconds are double precision
   nsec = hour * 3600.0 + minute * 60.0 + sec;

   //compute the local solar time
   temp = (2.0 * pi/(24.0)) * ((nsec/3600.0) + (longitude/15.0));

   //sanitize the output
   if (temp > 2.0 * pi) 
   {
       while (temp > 2.0 * pi)
          temp = temp-2.0 * pi;
   }
   else if (temp < 0.0) 
   {
       while (temp < 0.0)
          temp = temp + 2.0 * pi;
   }
   return temp;

} // function localsolartime





 //---------------------------------------------------------------------------
 //
 // FUNCTION: leap day
 //
 //> @brief returns a 1 if the input year is a leap year and 0 otherwise
 //
 //> @author Raul Dominguez
 //
 // INPUT ARGUMENTS:   int year         => year
 //
 // OUTPUT ARGUMENTS:  int leap_day     <= 1 if the year is leap year, 0 otherwise
 //
 // REVISION HISTORY: 31/07/2012 Initial version
 //
 //---------------------------------------------------------------------------
 int leap_day
	 (int year)
   {
	 if ((year % 400) == 0) 
         return 1;
     else if ((year % 100) == 0) 
         return 0;
     else if ((year % 4) == 0) 
         return 1;
     else
         return 0;
   }    // function leap_day

   

 //---------------------------------------------------------------------------
 //
 // FUNCTION: linear_interpolation
 //
 //> @brief perform a linear interpolation between (x1,y1) and (x2,y2)
 //
 //> @author Raul Dominguez
 //
 // INPUT ARGUMENTS:   real x1,y1,x2,y2,x => data vectors (x1,y1), (x2,y2) and
 //                                            point x in which we interpolate
 //
 // OUTPUT ARGUMENTS:  real linear_interpolation <= value of the interpolated
 //                                                   line at x
 //
 // REVISION HISTORY: 31/07/2012 Initial version
 //                   14/09/2012 Changed i/o from double to real
 //
 //---------------------------------------------------------------------------
 double linear_interpolation 
      (
	    double x1, double y1, double x2, double y2, double x 
      )
	{
		double a, b;

        a = (y1 - y2) / (x1 - x2);
        b = y1 - (a * x1);

        return a * x + b;

      //write(*,*) "interpolation:"
      //write(*,*) "X              Y"
      //write(*,*) x1,y1
      //write(*,*) x2,y2
      //write(*,*) "punto", x
      //write(*,*) "resultado", linear_interpolation
      //write(*,*)  " "
	}  // function linear_interpolation


 //---------------------------------------------------------------------------
 //
 // FUNCTION: hms2hr
 //
 //> @brief  convert a time from hours, minutes and seconds to hours (decimal)
 //
 //> @author Raul Dominguez
 //
 // INPUT ARGUMENTS:   int hr            => hours
 //                    int mins          => minutes
 //                    double  sec           => seconds
 //
 // OUTPUT ARGUMENTS:  double  hms2hr        <= decimal number of hours corresponding
 //                                             to the input
 //
 // REVISION HISTORY: 31/07/2012 Initial version
 //                   14/09/2012 Convert output from double to real
 //
 //---------------------------------------------------------------------------

 double hms2hr 
     (
	   int hr, int mins, double sec 
	 )
	{
		return hr + mins/60.0 + sec/3600.0;
	}  // function hms2hr



 //---------------------------------------------------------------------------
 //
 // FUNCTION: day_of_year
 //
 //> @brief this function returns the number of days elapsed since the beginning of the year
 //>              for a given date. the output ranges from 0 to 365 (or 366 days in leap years)
 //
 //> @author Raul Dominguez
 //
 //
 // INPUT ARGUMENTS:   int day         => day
 //                    int month       => month
 //                    int year        => year
 //
 // OUTPUT ARGUMENTS:  int day_of_year <= the day of year corresponding to the input date
 //
 // REVISION HISTORY: 31/07/2012 Initial version
 //
 //---------------------------------------------------------------------------
 
int day_of_year
	(
	  int day, int month, int year 
    )
{
   int month_length[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
   int i, temp;

   month_length[2] = 28 + leap_day(year); //count the leap day

   temp=0;
   for (i=1;i <= month-1; i++)
      temp = temp + month_length[i];
   
   return temp + day;

}  // function day_of_year


 //---------------------------------------------------------------------------
 //
 // FUNCTION: mdyhms_ok
 //
 //> @brief This function checks wheter the user has input a valid date or not
 //>         for example, a date like 34/01/2012 is invalid.
 //>         It returns .TRUE. if the input date is ok, and .FALSE. otherwise
 //
 //> @author Raul Dominguez
 //
 // INPUT ARGUMENTS:   int day         => day
 //                    int month       => month
 //                    int year        => year
 //                    int hour        => hour
 //                    int minute      => minute
 //                    double  second      => second
 //
 // OUTPUT ARGUMENTS:  logical mdyhms_ok <= .TRUE. if the user-entered parameters
 //                                                are a proper date. .FALSE.
 //                                                otherwise
 //
 // REVISION HISTORY: 06/09/2012 Initial version
 //
 //---------------------------------------------------------------------------

 bool mdyhms_ok
	 (
	    int day, int month, int year, int hour, int minute,
		double second 
     )
{
   //table with the duration of months
   int month_length[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

   month_length[2] = 28 + leap_day(year); //count the leap day

    //check day and month
    if ((month < 1) || (month > 12)) 
	{
     return false;
	}

    if ((day < 1) || (day > month_length[month])) 
	{
     return false;
	}

    if ((hour < 0) || (hour > 24)) 
	{
     return false;
	}

    if ((minute < 0) || (minute > 59)) 
	{
     return false;
	}

    if ((second < 0.0) || (second >= 60.0)) 
	{
     return false;
	}

    return true;

} // function mdyhms_ok

	 
 //---------------------------------------------------------------------------
 //
 // SUBROUTINE: DJ2000
 //
 //> @brief COMPUTES CALENDER DATE FROM MODIFIED JULIAN DAY 2000
 //>              VALID FOR DATES BETWEEN 1950/JAN/1 AND 2099/DEC/31.
 //>              MJD(2000) = MJD(1950) - 18262.0 IS = 0 ON 2000/01/01 AT 00:00:00.
 //
 //> @author legacy code
 //
 // INPUT ARGUMENTS:   (double) DAY    =>    MOD. JULIAN DAY, REFERRED TO 2000 (MAY BE NEGATIVE).
 //
 // OUTPUT ARGUMENTS:  (intS)      <=    I=YEAR, J=MONTH, K=DAY, JHR=HOUR, MI=MINUTE
 //                    (double)        <=    SEC=SECOND.
 //
 // REVISION HISTORY: legacy code
 //
 //---------------------------------------------------------------------------
 void DJ2000
	  ( 
	    double MJDay, int& year, int& mon, int& day, int& hr, int& minute, double& sec
	  )
 {
	 double JDAY;
	 int L, M, N, JJ;
	 //  MAKE SURE TO ROUND-OFF ONLY DOWN, ALSO FOR NEGATIVE MJD:
	 JDAY = MJDay + 18262.0;
	 L = floor((4000.0 * (JDAY + 18204)) / 1461001);
	 N = JDAY - floor((1461.0 * L) / 4) + 18234;
	 M = floor((80.0 * N) / 2447);
	 day = N - floor((2447.0 * M) / 80);
	 JJ = (int)floor(M / 11.0);
	 mon = M + 2 - 12 * JJ;
	 year = 1900 + L + JJ;
	 sec = (MJDay - (JDAY - 18262)) * 24.0;
	 hr = (int)floor(sec);
	 sec = (sec - (hr)) * 60.0;
	 minute = (int)floor(sec);
	 sec = (sec - (minute)) * 60.0;
 } // function DJ2000




 //---------------------------------------------------------------------------
 //
 // SUBROUTINE: JD2000
 //
 //> @brief GIVES THE NEW MOD. JULIAN DAY (MJD=0.0 ON 2000/JAN/1 AT 0:00:00)
 //>        FOR INPUT CALENDAR DATES BETWEEN 1950/JAN/1 AND 2099/DEC/31..
 //>
 //>    MJD(2000) = MJD(1950) - 18262.0
 //
 //> @author legacy code
 //
 // OUTPUT ARGUMENTS:   (double) DAY    =>    MOD. JULIAN DAY, REFERRED TO 2000 (MAY BE NEGATIVE).
 //
 // INPUT ARGUMENTS:  (intS)      <=    I=YEAR, J=MONTH, K=DAY, JHR=HOUR, MI=MINUTE
 //                    (double)        <=    SEC=SECOND.
 //
 // REVISION HISTORY: legacy code
 //
 //---------------------------------------------------------------------------

 void JD2000
	 (
	   double& MJDay, int year, int mon, int day, int hr, int minute, double sec
     )
{
	 double JJ;
	 int L;
	 // note floor returns double rounded towards zero, need to recast if int is needed, 
	 // and args must have a double in them
	 JJ = floor((14.0 - mon) / 12);
	 L = year - JJ - 1900 * (int)floor(year / 1900.0) + 100 * (int)floor(2000.0 / (year + 1951));
	 MJDay = day - 36496 + floor((1461.0 * L) / 4) + floor((367.0 * (mon - 2 + JJ * 12)) / 12);
	 MJDay = MJDay + (((hr * 60 + minute) * 60) + sec) / 86400.0;
 } // function JD2000



 

 //---------------------------------------------------------------------------
 //
 // SUBROUTINE: warning
 //
 //> @brief this routine comes into focus when a non-fatal error happens in
 //>              dtm_wrap module. These errors, altough non fatal, should
 //>              be reported to the final user. The user should modify this routine
 //>              to fit into his/her error handling scheme
 //
 //> @author  Raul Dominguez
 //
 // INPUT ARGUMENTS:   int code         => error code
 //
 // REVISION HISTORY: 31/07/2012 Initial version
 //
 //---------------------------------------------------------------------------
void warning 
    ( 
	    int code
    )
{
    typedef char ls[100];
	ls warning_messages[20];

    //these messages tell the user what happened
#ifdef _MSC_VER
	strcpy_s(warning_messages[1], " "); // Unused
    strcpy_s(warning_messages[2], "Unable to interpolate fbar. Returning a non-interpolated value");
#else
	strcpy(warning_messages[1], " "); // Unused
	strcpy(warning_messages[2], "Unable to interpolate fbar. Returning a non-interpolated value");
#endif

    printf( "\n\n---------------------------------- \n" );
    printf( "WARNING in dtm_wrap module\n" );
    printf( "%i", warning_messages[code] );
    printf( "-----------------------------------" );
} // function warning


 //---------------------------------------------------------------------------
 //
 // SUBROUTINE: process_input_date
 //
 //> @brief this routine fills the data in the provided in_date structure
 //>
 //> @details If the user-input date is in MJD2000, fill the fields containing
 //>          day month year hour minute seconds. If the user provided d-m-yr
 //>          h-m-s, fill the MJD2000 field
 //
 //> @author  Raul Dominguez
 //
 // INPUT ARGUMENTS:
 //                   @param[inout] in_date structure whose fields are to be filled
 //
 // REVISION HISTORY: 14/09/2012 Initial version
 //
 //---------------------------------------------------------------------------

void process_input_date 
	( 
	    dtm_daterectype& in_date
    )
{
    if (in_date.type_flag == 1) 
	{
         //Convert the user entered date to dd/mm/year hh:mm:ss.sssss
         DJ2000(in_date.mjd2000,in_date.year,in_date.month,in_date.day, 
                in_date.hour,in_date.minute,in_date.second);
	}
    else if (in_date.type_flag == 2) 
	{
        //the user gave an universal time in d-m-yr h-m-s format

        //First check that the input date is OK
        if (mdyhms_ok(in_date.day, in_date.month, in_date.year, in_date.hour, 
            in_date.minute, in_date.second) == false) 
		{
            fatal_error(10);
    	} 

        //also convert ymdhms to mjd2000
        JD2000(in_date.mjd2000, in_date.year, in_date.month, in_date.day, 
               in_date.hour, in_date.minute, in_date.second);
	}

} // function process_input_date



 //
 // initialize Space Weather
 // this routine simply opens and reads the geomagnetic and solar flux data
 //
void InitSPW
    (
	  char path_to_a_file[100], char path_to_f_file[100], int& number_of_lines_in_a_file,
      int& number_of_lines_in_f_file, a_structtype a_indexes[3000], f_structtype f_indexes[3000]
    )
{
	  char longstr1[140];
	  FILE *infile;
      //Counter indexes
      int file_io_status;
      int idx;
	  long int max_lines, tmpa, tmpb;
	  char temp8[9];
	  //save first_run,a_indexes,f_indexes, number_of_lines_in_a_file, number_of_lines_in_f_file;

	  //======================================================================
	  max_lines = 3000;

      //This block reads the "f" file
      file_io_status = 0;

#ifdef _MSC_VER
	  errno_t errs6 = fopen_s(&infile, path_to_f_file, "r");
#else
	  infile  = fopen( path_to_f_file, "r");
#endif	  
      // open (f_file_unit,file=path_to_f_file,iostat=file_io_status,status='old')

      if (infile == NULL)
	  {
	      printf("Failed to open file: %s\n", path_to_f_file);
//          exit;
      }

      //skip first 14 lines (header)
      for (idx=1; idx <=14; idx++)
		 fgets( longstr1,140,infile);

      //  read (f_file_unit,'(a1)') dummychar;
      idx=1;
	  while (feof(infile) == 0)
      {
          fgets( longstr1,140,infile);
          //3x,i4,3x,i2,3x,i2,3x,4x,3x,f4.0,3x,f4.0,3x,4x      
		  // %[^_] means to read until an underscore character is encountered
	  	  // use dummy string to bypass NaN's
#ifdef _MSC_VER
		  sscanf_s(longstr1, " %i %i %i %8s %lf %lf ", &f_indexes[idx].year, &f_indexes[idx].month, &f_indexes[idx].day,
				   &temp8, 9 * sizeof(char), &f_indexes[idx].f, &f_indexes[idx].fbar );
#else
		  sscanf(longstr1, " %i %i %i %8s %lf %lf ", &f_indexes[idx].year, &f_indexes[idx].month, &f_indexes[idx].day,
			     &temp8, &f_indexes[idx].f, &f_indexes[idx].fbar);
#endif
		  idx = idx + 1;

          //stop the routine if the file contains more lines than the max_lines parameter
          if (idx > max_lines) 
		  {
		      fatal_error(5);
//		      exit;
		  }
      }    // (strcmp(longstr1, "#")!=0) && 

		  number_of_lines_in_f_file = idx - 1;
          fclose(infile);

          //This block reads the "a" file
          file_io_status = 0;
          
		  // cdav read in geomagnetic indices
#ifdef _MSC_VER
		  errno_t errs5 = fopen_s(&infile, path_to_a_file, "r");
#else
		  infile = fopen(path_to_a_file, "r");
#endif	  

          //	open (a_file_unit,file=path_to_a_file,iostat=file_io_status,status='old')

          if (infile == NULL)
	      {
	          printf("Failed to open file: %s\n", path_to_f_file);
//              exit;
          }
        
          idx=1;
          while (feof(infile) == 0)
  		  {
			  fgets( longstr1, 130, infile);
              //               4i4,4x,9i4   cdav chg to dayofyear
              // use d format for integers with leading 0's
			  // sscanf_s doesn't appear to be able to parse strings "and" values, so extract all values with atoi, then do
			  // substrings
			  //&temp4a, 4 * sizeof(char), &temp4b, 3 * sizeof(char), &temp4c, 4 * sizeof(char),
			  //&a_indexes[idx].month, &a_indexes[idx].day, &a_indexes[idx].year,
			  //   1   12006   12193  12  14   7   7   8  21  16  13  11
			  //   1   22006   22194  15  16  19  16  18   6  12  11  29
			  //   1   32006   32195   6  10   3   3   5   8  15   8   2
#ifdef _MSC_VER
			  // some tricks. reading leading whitespace with string doesn't work. read starts at first non white space. c++ is lousy
			  sscanf_s(longstr1, "%i %i %i %4i %4i %4i %4i %4i %4i %4i %4i %4i ", &a_indexes[idx].month, &tmpa, &tmpb, 
			           &a_indexes[idx].mean_am, &a_indexes[idx].am[1], &a_indexes[idx].am[2], &a_indexes[idx].am[3],
				       &a_indexes[idx].am[4], &a_indexes[idx].am[5], &a_indexes[idx].am[6], &a_indexes[idx].am[7], 
					   &a_indexes[idx].am[8]);
#else
			  sscanf_s(longstr1, "%i %i %i %4i %4i %4i %4i %4i %4i %4i %4i %4i ", &a_indexes[idx].month, &tmpa, &tmpb, 
				  &a_indexes[idx].mean_am, &a_indexes[idx].am[1], &a_indexes[idx].am[2], &a_indexes[idx].am[3],
				  &a_indexes[idx].am[4], &a_indexes[idx].am[5], &a_indexes[idx].am[6], &a_indexes[idx].am[7], 
				  &a_indexes[idx].am[8]);
#endif
			  //strncpy_s(temp4a,&longstr1[0],4);  // get substring value
              //strncpy_s(temp4b,&longstr1[4],4);
              //strncpy_s(temp4c,&longstr1[8],4);
			  //a_indexes[idx].month = atoi(temp4a);
			  //a_indexes[idx].day = atoi(temp4b);
			  //a_indexes[idx].year = atoi(temp4c);
			  a_indexes[idx].day = int(floor(tmpa/10000.0));
			  a_indexes[idx].year = tmpa - a_indexes[idx].day * 10000;
//			  a_indexes[idx].dayofyear =  int(k/10000);

              idx = idx + 1;
		
              //stop the routine if the file contains more lines than the max_lines parameter
              if (idx > max_lines) 
			  { 
                  fatal_error(8);
//				  exit;
			  }
	  	  } //  end while

		  number_of_lines_in_a_file = idx - 1;
		  
          fclose(infile);

}  // function InitSPW


// http://www.daniweb.com/software-development/c/threads/395169/converting-string-to-float-without-atof-for-scientific-notation
// c++ is lousy at strings!!
void ScientificToDouble(char *in_String, double& tepom) 
{
	// Loop Variables
	int        Counter       = 0;
	int        Length        = strlen(in_String) + 1;
	// Flags and signs
	int        NegativeFlag  = 0;
	int        DecimalFlag   = 0;
	int        ExponentSign  = 0;  // -1 = Negative, 0 = None, 1 = Positive
	// Numerical Data
	int        Exponent      = 0;
	int        FinalDivision = 1;
	long       Digits        = 0;
//	double tepom;

	// Loop per each character. Ignore anything weird.
	for (;Counter < Length; Counter++) 
	{
		// Depending on the current character
		switch (in_String[Counter]) 
		{
			// On any digit
			case '0': case '5':
			case '1': case '6':
			case '2': case '7':
			case '3': case '8':
			case '4': case '9':
				// If we haven't reached an exponent yet ("e")
				if (ExponentSign == 0) 
				{
					// Adjust the final division if a decimal was encountered
					if (DecimalFlag) FinalDivision *= 10;
					// Add a digit to our main number
					Digits = (Digits * 10) + (in_String[Counter] - '0');
				// If we passed an "e" at some point
				} else 
				{
					// Add a digit to our exponent
					Exponent = (Exponent * 10) + (in_String[Counter] - '0');
				}
				break;
			// On a negative sign
			case '-':
				// If we passed an 'e'
				if (ExponentSign > 0)
					// The exponent sign will be negative
					ExponentSign = -1;
				// Otherwise we are still dealing with the main number
				else
					// Set the negative flag. We will negate the main number later.
					NegativeFlag = 1;
				break;
			// If we encounter some kind of "e"
			case 'e': case 'E':
				// Set the exponent flag
				ExponentSign = 1;
				break;
			// If we encounter a period
			case '.':
				// Set the decimal flag. We will start tracking decimal depth.
				DecimalFlag = 1;
				break;
			// We gladly accept all sorts of additional garbage.
			default:
				break;
		}
	}
	// If the negative flag is set, negate the main number
	if (NegativeFlag)
		Digits = 0 - Digits;
	
	// If the exponent is supposed to be negative, negate it now
	if (ExponentSign < 0)
		Exponent = 0 - Exponent;
	// Return the calculated result of our observations
	tepom = double(Digits)/double(FinalDivision) * double(pow(10.0, Exponent));
	//return tepom;
}  // function ScientificToDouble


//---------------------------------------------------------------------------
//
// ROUTINE: DTMInit
//
//> @author Sean Bruinsma
//>
//> @brief  read DTM format/version 2012
//
// PROTOTYPE:
//                 P_ReadDTM12 (
//                             )
//
//
// INPUT ARGUMENTS:
//
//
// OUTPUT ARGUMENTS:
//
//> @date 03/2012
//
//
// * rol Programme: read DTM format/version 2012 
//
//      pgf90  -o P_ReadDTM12.x P_ReadDTM12.f90
//
// * vers      : 03/2012
// * aut      : PS
//
//---------------------------------------------------------------------------

void DTMInit
	(
	  char PATH_TO_DEN_UNC_1[100], char PATH_TO_DEN_UNC_2[100], char PATH_TO_DEN_UNC_3[100],
	  char PATH_TO_DTM2012_IN_FILE[100], char PATH_TO_DTM2012_OUT_FILE[100], int IVERIF,
      double nomgrid[25][19], double kp_scale[25][19], double alt_scale[25][19],
	  pardtmtype& pardtm
	)
{ 
    // int IVERIF = 0;  //// IVERIF = 0/1 >> print model in u_Out
    // _._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._._.
    //.. VARIABLES .. 
//	dtmindextype dtmindex;  // needed because numbers weren't passing correctly
	int  i, j, j_test, npdtm;
    int  Iok;
    char titre[110], Com[210];
  
    FILE *infile, *outfile;
    char longstr1[5200];
    FILE *infile1, *infile2, *infile3;
    double tempo;
    char tempstr[12];
    char infilename[100];
  
  // --------------------------------------------------
#ifdef _MSC_VER
	strcpy_s(infilename, sizeof(infilename), PATH_TO_DEN_UNC_1);
	errno_t errs1 = fopen_s(&infile1, infilename, "r");
#else
	strcpy(infilename, PATH_TO_DEN_UNC_1);
	infile1 = fopen(infilename, "r");
#endif

#ifdef _MSC_VER
	strcpy_s(infilename, sizeof(infilename), PATH_TO_DEN_UNC_2);
	errno_t errs2 = fopen_s(&infile2, infilename, "r");
#else
	strcpy(infilename, PATH_TO_DEN_UNC_2);
	infile2 = fopen(infilename, "r");
#endif

#ifdef _MSC_VER
	strcpy_s(infilename, sizeof(infilename), PATH_TO_DEN_UNC_3);
	errno_t errs3 = fopen_s(&infile3, infilename, "r");
#else
	strcpy(infilename, PATH_TO_DEN_UNC_3);
	infile3 = fopen(infilename, "r");
#endif

    // cdav may need to shift each by 1...
    // do a slow but works way. The msvs2010 debugger shows the values incorrect, 
    // but correct AFTER the call
	fgets( longstr1,5200,infile1);  // single line 5185 char
	for (j=0; j <18; j++)
    {
    	for (i=0; i <24; i++)
        {
#ifdef _MSC_VER
			strncpy_s(tempstr, &longstr1[j * 276 + (i + j) * 12], 11); // starts at 0
#else
			strncpy(tempstr, &longstr1[j * 276 + (i + j) * 12]); // starts at 0
#endif
			tempstr[11] = '\0';  // each variable is just 11 chars long. put endofstr marker AT 11
            ScientificToDouble(tempstr, tempo);
            nomgrid[i+1][j+1] = tempo;
	    }  // through i
    }  // through j

	fgets( longstr1,5200,infile2);  // single line 5185 char
	for (j=0; j <18; j++)
    {
    	for (i=0; i <24; i++)
        {
#ifdef _MSC_VER
			strncpy_s(tempstr, &longstr1[j * 276 + (i + j) * 12], 11); // starts at 0
#else
			strncpy(tempstr, &longstr1[j * 276 + (i + j) * 12]); // starts at 0
#endif
			tempstr[11] = '\0';  // each variable is just 11 chars long. put endofstr marker AT 11
            ScientificToDouble(tempstr, tempo);
            kp_scale[i+1][j+1] = tempo;
	    }  // through i
    }  // through j

	fgets( longstr1,5200,infile3);  // single line 5185 char
	for (j=0; j <18; j++)
    {
    	for (i=0; i <24; i++)
        {
#ifdef _MSC_VER
			strncpy_s(tempstr, &longstr1[j * 276 + (i + j) * 12], 11); // starts at 0
#else
			strncpy(tempstr, &longstr1[j * 276 + (i + j) * 12]); // starts at 0
#endif
			tempstr[11] = '\0';  // each variable is just 11 chars long. put endofstr marker AT 11
            ScientificToDouble(tempstr, tempo);
            alt_scale[i+1][j+1] = tempo;
    	}  // through i
    }  // through j

//		fgets( longstr1,5200,infile3);  // single line 5185 char
//		sscanf_s(longstr1,"%{%g} ", &alt_scale );

    fclose (infile1);
    fclose (infile2);
    fclose (infile3);

#ifdef _MSC_VER
	errno_t errs4 = fopen_s(&infile, PATH_TO_DTM2012_IN_FILE, "r");
#else
	infile = fopen(PATH_TO_DTM2012_IN_FILE, "r");
#endif	  
   
	if (IVERIF == 1)
	{
#ifdef _MSC_VER
		errno_t errs5 = fopen_s(&outfile, PATH_TO_DTM2012_OUT_FILE, "w");
#else
		outfile = fopen(PATH_TO_DTM2012_OUT_FILE, "w");
#endif	  
	}

    // --------------------------------------------------
    fgets( longstr1,110,infile);
#ifdef _MSC_VER
	strcpy_s(titre, sizeof(titre), longstr1);
#else
	strcpy(titre, longstr1);
#endif
	//scanf( "%s",&titre );
    fgets( longstr1,210,infile);
    sscanf_s( longstr1,"%i",&npdtm );

    //	  fgets( longstr1,210,infile); re-process the line
#ifdef _MSC_VER
	strcpy_s(Com, sizeof(Com), longstr1);
#else
	strcpy(Com, longstr1);
#endif
	//      scanf( "%s",&Com );
    //      backspace(u_In);

    if (IVERIF == 1) 
    {
        fprintf(outfile, "%110s \n", titre);
        fprintf(outfile, "%210s \n",  Com);
	    printf( "%110s \n",titre);
        printf( "%210s \n", Com);
    }

    // --------------------------------------------------
    Iok = 1;    // test read ok/nok
     
    // ...... termes ONE
    i = 1;
    fgets( longstr1,210,infile);
#ifdef _MSC_VER
	sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
		           &j_test,  &pardtm.tt.T_one  ,&pardtm.dtt.T_one, 
                   &pardtm.h.T_one   ,&pardtm.dh.T_one, &pardtm.he.T_one  ,&pardtm.dhe.T_one, 
                   &pardtm.ox.T_one  ,&pardtm.dox.T_one, &pardtm.az2.T_one ,&pardtm.daz2.T_one, 
                   &pardtm.o2.T_one  ,&pardtm.do2.T_one, &pardtm.az.T_one  ,&pardtm.daz.T_one, 
                   &pardtm.t0.T_one  ,&pardtm.dt0.T_one, &pardtm.tp.T_one  ,&pardtm.dtp.T_one);
#else
	sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
		&j_test, &pardtm.tt.T_one, &pardtm.dtt.T_one,
		&pardtm.h.T_one, &pardtm.dh.T_one, &pardtm.he.T_one, &pardtm.dhe.T_one,
		&pardtm.ox.T_one, &pardtm.dox.T_one, &pardtm.az2.T_one, &pardtm.daz2.T_one,
		&pardtm.o2.T_one, &pardtm.do2.T_one, &pardtm.az.T_one, &pardtm.daz.T_one,
		&pardtm.t0.T_one, &pardtm.dt0.T_one, &pardtm.tp.T_one, &pardtm.dtp.T_one);
#endif
	//// ..:
     if (i != j_test) 
     {
         printf( " * * * WARNING: PB in field T_one " );
         printf( " * * * WARNING: incompatibility between dimension T_lat & DTM file  " );
         printf( " * * *  i = %6i j_test = %6i ", i, j_test);
         Iok = 0; 
         goto fifteen;
	 }


// ...... termes in LAT
	 // cdav lat all these also have to be < and not <= because the array dimension is 1 more
    for (i = 1; i <  Nb_lat; i++ )
 	{
        ////    READ(u_In,640) i,(DTM_12(j).T_lat[i],ddtm_12(j).T_lat[i],j = 1,9)
        fgets( longstr1,210,infile);
#ifdef _MSC_VER
		sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_lat[i]  ,&pardtm.dtt.T_lat[i], 
                  &pardtm.h.T_lat[i]   ,&pardtm.dh.T_lat[i], &pardtm.he.T_lat[i]  ,&pardtm.dhe.T_lat[i], 
                  &pardtm.ox.T_lat[i]  ,&pardtm.dox.T_lat[i], &pardtm.az2.T_lat[i] ,&pardtm.daz2.T_lat[i], 
                  &pardtm.o2.T_lat[i]  ,&pardtm.do2.T_lat[i], &pardtm.az.T_lat[i]  ,&pardtm.daz.T_lat[i], 
                  &pardtm.t0.T_lat[i]  ,&pardtm.dt0.T_lat[i], &pardtm.tp.T_lat[i]  ,&pardtm.dtp.T_lat[i]);
#else
		sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			&j_test, &pardtm.tt.T_lat[i], &pardtm.dtt.T_lat[i],
			&pardtm.h.T_lat[i], &pardtm.dh.T_lat[i], &pardtm.he.T_lat[i], &pardtm.dhe.T_lat[i],
			&pardtm.ox.T_lat[i], &pardtm.dox.T_lat[i], &pardtm.az2.T_lat[i], &pardtm.daz2.T_lat[i],
			&pardtm.o2.T_lat[i], &pardtm.do2.T_lat[i], &pardtm.az.T_lat[i], &pardtm.daz.T_lat[i],
			&pardtm.t0.T_lat[i], &pardtm.dt0.T_lat[i], &pardtm.tp.T_lat[i], &pardtm.dtp.T_lat[i]);
#endif
		if (i != j_test) 
          {
             printf( " * * * WARNING: PB in field T_lat " );
             printf(  " * * * WARNING: incompatibility between dimension T_lat & DTM file  " );
             printf( " * * *  i = %6i j_test = %6i ", i, j_test);
             Iok = 0;
             goto fifteen;
          }
 	}

// --------------------------------------------------
// ...... termes in Flux
	// cdav flux
    for (i = 1; i <  Nb_flux; i++ )
    {
        fgets( longstr1,210,infile);
#ifdef _MSC_VER
		sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_flux[i]  ,&pardtm.dtt.T_flux[i], 
                  &pardtm.h.T_flux[i]   ,&pardtm.dh.T_flux[i], &pardtm.he.T_flux[i]  ,&pardtm.dhe.T_flux[i], 
                  &pardtm.ox.T_flux[i]  ,&pardtm.dox.T_flux[i], &pardtm.az2.T_flux[i] ,&pardtm.daz2.T_flux[i], 
                  &pardtm.o2.T_flux[i]  ,&pardtm.do2.T_flux[i], &pardtm.az.T_flux[i]  ,&pardtm.daz.T_flux[i], 
                  &pardtm.t0.T_flux[i]  ,&pardtm.dt0.T_flux[i], &pardtm.tp.T_flux[i]  ,&pardtm.dtp.T_flux[i] );
#else
		sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			&j_test, &pardtm.tt.T_flux[i], &pardtm.dtt.T_flux[i],
			&pardtm.h.T_flux[i], &pardtm.dh.T_flux[i], &pardtm.he.T_flux[i], &pardtm.dhe.T_flux[i],
			&pardtm.ox.T_flux[i], &pardtm.dox.T_flux[i], &pardtm.az2.T_flux[i], &pardtm.daz2.T_flux[i],
			&pardtm.o2.T_flux[i], &pardtm.do2.T_flux[i], &pardtm.az.T_flux[i], &pardtm.daz.T_flux[i],
			&pardtm.t0.T_flux[i], &pardtm.dt0.T_flux[i], &pardtm.tp.T_flux[i], &pardtm.dtp.T_flux[i]);
#endif
		if (i != j_test) 
        	{
               printf( " * * * WARNING: PB in field T_flux " );
               printf( " * * * WARNING: incompatibility between dimension T_flux & DTM file  " );
               printf( " * * *  i = %6i j_test = %6i ", i, j_test);
               Iok = 0; 
               goto fifteen;
            }
	}

// --------------------------------------------------
// ...... termes in kp
    for (i = 1; i <   Nb_kp; i++ )
	{
         fgets( longstr1,210,infile);
#ifdef _MSC_VER
		 sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_kp[i]  ,&pardtm.dtt.T_kp[i], 
                  &pardtm.h.T_kp[i]   ,&pardtm.dh.T_kp[i], &pardtm.he.T_kp[i]  ,&pardtm.dhe.T_kp[i], 
                  &pardtm.ox.T_kp[i]  ,&pardtm.dox.T_kp[i], &pardtm.az2.T_kp[i] ,&pardtm.daz2.T_kp[i], 
                  &pardtm.o2.T_kp[i]  ,&pardtm.do2.T_kp[i], &pardtm.az.T_kp[i]  ,&pardtm.daz.T_kp[i], 
                  &pardtm.t0.T_kp[i]  ,&pardtm.dt0.T_kp[i], &pardtm.tp.T_kp[i]  ,&pardtm.dtp.T_kp[i] );
#else
		 sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			 &j_test, &pardtm.tt.T_kp[i], &pardtm.dtt.T_kp[i],
			 &pardtm.h.T_kp[i], &pardtm.dh.T_kp[i], &pardtm.he.T_kp[i], &pardtm.dhe.T_kp[i],
			 &pardtm.ox.T_kp[i], &pardtm.dox.T_kp[i], &pardtm.az2.T_kp[i], &pardtm.daz2.T_kp[i],
			 &pardtm.o2.T_kp[i], &pardtm.do2.T_kp[i], &pardtm.az.T_kp[i], &pardtm.daz.T_kp[i],
			 &pardtm.t0.T_kp[i], &pardtm.dt0.T_kp[i], &pardtm.tp.T_kp[i], &pardtm.dtp.T_kp[i]);
#endif
		 if (i != j_test)
	        {
              printf( " * * * WARNING: PB in field T_kp " );
              printf( " * * * WARNING: incompatibility between dimension T_kp & DTM file  " );
              printf( " * * *  i = %6i j_test = %6i ", i, j_test);
              Iok = 0;
              goto fifteen;
            }
	}

// --------------------------------------------------
// ...... termes in SLat
    for (i = 1; i <   Nb_SLat; i++ )
	{   
        fgets( longstr1,210,infile);
#ifdef _MSC_VER
		sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_SLat[i]  ,&pardtm.dtt.T_SLat[i], 
                  &pardtm.h.T_SLat[i]   ,&pardtm.dh.T_SLat[i], &pardtm.he.T_SLat[i]  ,&pardtm.dhe.T_SLat[i], 
                  &pardtm.ox.T_SLat[i]  ,&pardtm.dox.T_SLat[i], &pardtm.az2.T_SLat[i] ,&pardtm.daz2.T_SLat[i], 
                  &pardtm.o2.T_SLat[i]  ,&pardtm.do2.T_SLat[i], &pardtm.az.T_SLat[i]  ,&pardtm.daz.T_SLat[i], 
                  &pardtm.t0.T_SLat[i]  ,&pardtm.dt0.T_SLat[i], &pardtm.tp.T_SLat[i]  ,&pardtm.dtp.T_SLat[i] );
#else
		sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			&j_test, &pardtm.tt.T_SLat[i], &pardtm.dtt.T_SLat[i],
			&pardtm.h.T_SLat[i], &pardtm.dh.T_SLat[i], &pardtm.he.T_SLat[i], &pardtm.dhe.T_SLat[i],
			&pardtm.ox.T_SLat[i], &pardtm.dox.T_SLat[i], &pardtm.az2.T_SLat[i], &pardtm.daz2.T_SLat[i],
			&pardtm.o2.T_SLat[i], &pardtm.do2.T_SLat[i], &pardtm.az.T_SLat[i], &pardtm.daz.T_SLat[i],
			&pardtm.t0.T_SLat[i], &pardtm.dt0.T_SLat[i], &pardtm.tp.T_SLat[i], &pardtm.dtp.T_SLat[i]);
#endif
		if (i != j_test)
	        {
              printf( " * * * WARNING: PB in field T_SLat " );
              printf( " * * * WARNING: incompatibility between dimension T_SLat & DTM file  " );
              printf( " * * *  i = %6i j_test = %6i ", i, j_test);
              Iok = 0; 
              goto fifteen;
            }
	}

// --------------------------------------------------
// ...... termes in SASLat
    for (i = 1; i <   Nb_SASLat; i++ )
	{
        fgets( longstr1,210,infile);
#ifdef _MSC_VER
		sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_SASLat[i]  ,&pardtm.dtt.T_SASLat[i], 
                  &pardtm.h.T_SASLat[i]   ,&pardtm.dh.T_SASLat[i], &pardtm.he.T_SASLat[i]  ,&pardtm.dhe.T_SASLat[i], 
                  &pardtm.ox.T_SASLat[i]  ,&pardtm.dox.T_SASLat[i], &pardtm.az2.T_SASLat[i] ,&pardtm.daz2.T_SASLat[i], 
                  &pardtm.o2.T_SASLat[i]  ,&pardtm.do2.T_SASLat[i], &pardtm.az.T_SASLat[i]  ,&pardtm.daz.T_SASLat[i], 
                  &pardtm.t0.T_SASLat[i]  ,&pardtm.dt0.T_SASLat[i], &pardtm.tp.T_SASLat[i]  ,&pardtm.dtp.T_SASLat[i] );
#else
		sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			&j_test, &pardtm.tt.T_SASLat[i], &pardtm.dtt.T_SASLat[i],
			&pardtm.h.T_SASLat[i], &pardtm.dh.T_SASLat[i], &pardtm.he.T_SASLat[i], &pardtm.dhe.T_SASLat[i],
			&pardtm.ox.T_SASLat[i], &pardtm.dox.T_SASLat[i], &pardtm.az2.T_SASLat[i], &pardtm.daz2.T_SASLat[i],
			&pardtm.o2.T_SASLat[i], &pardtm.do2.T_SASLat[i], &pardtm.az.T_SASLat[i], &pardtm.daz.T_SASLat[i],
			&pardtm.t0.T_SASLat[i], &pardtm.dt0.T_SASLat[i], &pardtm.tp.T_SASLat[i], &pardtm.dtp.T_SASLat[i]);
#endif
		if (i != j_test)
	        {
              printf( " * * * WARNING: PB in field T_SASLat " );
              printf( " * * * WARNING: incompatibility between dimension T_SASLat & DTM file  " );
              printf( " * * *  i = %6i j_test = %6i ", i, j_test);
              Iok = 0;
              goto fifteen;
         }
	}

// --------------------------------------------------
// ...... termes in T_NSLat
    for (i = 1; i <   Nb_NSLat; i++ )   
	{
        fgets( longstr1,210,infile);
#ifdef _MSC_VER
		sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_NSLat[i]  ,&pardtm.dtt.T_NSLat[i], 
                  &pardtm.h.T_NSLat[i]   ,&pardtm.dh.T_NSLat[i], &pardtm.he.T_NSLat[i]  ,&pardtm.dhe.T_NSLat[i], 
                  &pardtm.ox.T_NSLat[i]  ,&pardtm.dox.T_NSLat[i], &pardtm.az2.T_NSLat[i] ,&pardtm.daz2.T_NSLat[i], 
                  &pardtm.o2.T_NSLat[i]  ,&pardtm.do2.T_NSLat[i], &pardtm.az.T_NSLat[i]  ,&pardtm.daz.T_NSLat[i], 
                  &pardtm.t0.T_NSLat[i]  ,&pardtm.dt0.T_NSLat[i], &pardtm.tp.T_NSLat[i]  ,&pardtm.dtp.T_NSLat[i] );
#else
		sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			&j_test, &pardtm.tt.T_NSLat[i], &pardtm.dtt.T_NSLat[i],
			&pardtm.h.T_NSLat[i], &pardtm.dh.T_NSLat[i], &pardtm.he.T_NSLat[i], &pardtm.dhe.T_NSLat[i],
			&pardtm.ox.T_NSLat[i], &pardtm.dox.T_NSLat[i], &pardtm.az2.T_NSLat[i], &pardtm.daz2.T_NSLat[i],
			&pardtm.o2.T_NSLat[i], &pardtm.do2.T_NSLat[i], &pardtm.az.T_NSLat[i], &pardtm.daz.T_NSLat[i],
			&pardtm.t0.T_NSLat[i], &pardtm.dt0.T_NSLat[i], &pardtm.tp.T_NSLat[i], &pardtm.dtp.T_NSLat[i]);
#endif
		if (i != j_test)
	         {
                printf( " * * * WARNING: PB in field T_NSLat " );
                printf( " * * * WARNING: incompatibility between dimension T_NSLat & DTM file  " );
                printf( " * * *  i = %6i j_test = %6i ", i, j_test);
                Iok = 0; 
                goto fifteen;
             }
	}

// --------------------------------------------------
// ...... termes in T_SANSLat 
    for (i = 1; i <   Nb_SANSLat; i++ )   
	{
         fgets( longstr1,210,infile);
#ifdef _MSC_VER
		 sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_SANSLat[i]  ,&pardtm.dtt.T_SANSLat[i], 
                  &pardtm.h.T_SANSLat[i]   ,&pardtm.dh.T_SANSLat[i], &pardtm.he.T_SANSLat[i]  ,&pardtm.dhe.T_SANSLat[i], 
                  &pardtm.ox.T_SANSLat[i]  ,&pardtm.dox.T_SANSLat[i], &pardtm.az2.T_SANSLat[i] ,&pardtm.daz2.T_SANSLat[i], 
                  &pardtm.o2.T_SANSLat[i]  ,&pardtm.do2.T_SANSLat[i], &pardtm.az.T_SANSLat[i]  ,&pardtm.daz.T_SANSLat[i], 
                  &pardtm.t0.T_SANSLat[i]  ,&pardtm.dt0.T_SANSLat[i], &pardtm.tp.T_SANSLat[i]  ,&pardtm.dtp.T_SANSLat[i] );
#else
		 sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			 &j_test, &pardtm.tt.T_SANSLat[i], &pardtm.dtt.T_SANSLat[i],
			 &pardtm.h.T_SANSLat[i], &pardtm.dh.T_SANSLat[i], &pardtm.he.T_SANSLat[i], &pardtm.dhe.T_SANSLat[i],
			 &pardtm.ox.T_SANSLat[i], &pardtm.dox.T_SANSLat[i], &pardtm.az2.T_SANSLat[i], &pardtm.daz2.T_SANSLat[i],
			 &pardtm.o2.T_SANSLat[i], &pardtm.do2.T_SANSLat[i], &pardtm.az.T_SANSLat[i], &pardtm.daz.T_SANSLat[i],
			 &pardtm.t0.T_SANSLat[i], &pardtm.dt0.T_SANSLat[i], &pardtm.tp.T_SANSLat[i], &pardtm.dtp.T_SANSLat[i]);
#endif
		 if (i != j_test)
	         {
                printf( " * * * WARNING: PB in field T_SANSLat " );
                printf( " * * * WARNING: incompatibility between dimension T_SANSLat & DTM file  " );
                printf( " * * *  i = %6i j_test = %6i ", i, j_test);
                Iok = 0; 
                goto fifteen;
             }
	}

// --------------------------------------------------
// ...... termes in DiAn
    for (i = 1; i <   Nb_DiAn; i++ )   
	{
         fgets( longstr1,210,infile);
#ifdef _MSC_VER
		 sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_DiAn[i]  ,&pardtm.dtt.T_DiAn[i], 
                  &pardtm.h.T_DiAn[i]   ,&pardtm.dh.T_DiAn[i], &pardtm.he.T_DiAn[i]  ,&pardtm.dhe.T_DiAn[i], 
                  &pardtm.ox.T_DiAn[i]  ,&pardtm.dox.T_DiAn[i], &pardtm.az2.T_DiAn[i] ,&pardtm.daz2.T_DiAn[i], 
                  &pardtm.o2.T_DiAn[i]  ,&pardtm.do2.T_DiAn[i], &pardtm.az.T_DiAn[i]  ,&pardtm.daz.T_DiAn[i], 
                  &pardtm.t0.T_DiAn[i]  ,&pardtm.dt0.T_DiAn[i], &pardtm.tp.T_DiAn[i]  ,&pardtm.dtp.T_DiAn[i] );
#else
		 sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			 &j_test, &pardtm.tt.T_DiAn[i], &pardtm.dtt.T_DiAn[i],
			 &pardtm.h.T_DiAn[i], &pardtm.dh.T_DiAn[i], &pardtm.he.T_DiAn[i], &pardtm.dhe.T_DiAn[i],
			 &pardtm.ox.T_DiAn[i], &pardtm.dox.T_DiAn[i], &pardtm.az2.T_DiAn[i], &pardtm.daz2.T_DiAn[i],
			 &pardtm.o2.T_DiAn[i], &pardtm.do2.T_DiAn[i], &pardtm.az.T_DiAn[i], &pardtm.daz.T_DiAn[i],
			 &pardtm.t0.T_DiAn[i], &pardtm.dt0.T_DiAn[i], &pardtm.tp.T_DiAn[i], &pardtm.dtp.T_DiAn[i] );
#endif
		 if (i != j_test)
	         {
                printf( " * * * WARNING: PB in field T_DiAn " );
                printf( " * * * WARNING: incompatibility between dimension T_DiAn & DTM file  " );
                printf( " * * *  i = %6i j_test = %6i ", i, j_test);
                Iok = 0; 
                goto fifteen;
               }
	}

// --------------------------------------------------
// ...... termes  SDiAn
    for (i = 1; i <   Nb_SDiAn; i++ )
	{
         fgets( longstr1,210,infile);
#ifdef _MSC_VER
		 sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_SDiAn[i]  ,&pardtm.dtt.T_SDiAn[i], 
                  &pardtm.h.T_SDiAn[i]   ,&pardtm.dh.T_SDiAn[i], &pardtm.he.T_SDiAn[i]  ,&pardtm.dhe.T_SDiAn[i], 
                  &pardtm.ox.T_SDiAn[i]  ,&pardtm.dox.T_SDiAn[i], &pardtm.az2.T_SDiAn[i] ,&pardtm.daz2.T_SDiAn[i], 
                  &pardtm.o2.T_SDiAn[i]  ,&pardtm.do2.T_SDiAn[i], &pardtm.az.T_SDiAn[i]  ,&pardtm.daz.T_SDiAn[i], 
                  &pardtm.t0.T_SDiAn[i]  ,&pardtm.dt0.T_SDiAn[i], &pardtm.tp.T_SDiAn[i]  ,&pardtm.dtp.T_SDiAn[i] );
#else
		 sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			 &j_test, &pardtm.tt.T_SDiAn[i], &pardtm.dtt.T_SDiAn[i],
			 &pardtm.h.T_SDiAn[i], &pardtm.dh.T_SDiAn[i], &pardtm.he.T_SDiAn[i], &pardtm.dhe.T_SDiAn[i],
			 &pardtm.ox.T_SDiAn[i], &pardtm.dox.T_SDiAn[i], &pardtm.az2.T_SDiAn[i], &pardtm.daz2.T_SDiAn[i],
			 &pardtm.o2.T_SDiAn[i], &pardtm.do2.T_SDiAn[i], &pardtm.az.T_SDiAn[i], &pardtm.daz.T_SDiAn[i],
			 &pardtm.t0.T_SDiAn[i], &pardtm.dt0.T_SDiAn[i], &pardtm.tp.T_SDiAn[i], &pardtm.dtp.T_SDiAn[i]);
#endif
		 if (i != j_test)
	         {
                printf(" * * * WARNING: PB in field T_SDiAn " );
                printf(" * * * WARNING: incompatibility between dimension T_SDiAn & DTM file  " );
                printf( " * * *  i = %6i j_test = %6i ", i, j_test);
                Iok = 0; 
                goto fifteen;
               }
	}

// ------------- -------------------------------------
// ...... termes in  TDi
    for (i = 1; i <   Nb_TDi; i++ )
	{
         fgets( longstr1,210,infile);
#ifdef _MSC_VER
		 sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_TDi[i]  ,&pardtm.dtt.T_TDi[i], 
                  &pardtm.h.T_TDi[i]   ,&pardtm.dh.T_TDi[i], &pardtm.he.T_TDi[i]  ,&pardtm.dhe.T_TDi[i], 
                  &pardtm.ox.T_TDi[i]  ,&pardtm.dox.T_TDi[i], &pardtm.az2.T_TDi[i] ,&pardtm.daz2.T_TDi[i], 
                  &pardtm.o2.T_TDi[i]  ,&pardtm.do2.T_TDi[i], &pardtm.az.T_TDi[i]  ,&pardtm.daz.T_TDi[i], 
                  &pardtm.t0.T_TDi[i]  ,&pardtm.dt0.T_TDi[i], &pardtm.tp.T_TDi[i]  ,&pardtm.dtp.T_TDi[i] );
#else
		 sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			 &j_test, &pardtm.tt.T_TDi[i], &pardtm.dtt.T_TDi[i],
			 &pardtm.h.T_TDi[i], &pardtm.dh.T_TDi[i], &pardtm.he.T_TDi[i], &pardtm.dhe.T_TDi[i],
			 &pardtm.ox.T_TDi[i], &pardtm.dox.T_TDi[i], &pardtm.az2.T_TDi[i], &pardtm.daz2.T_TDi[i],
			 &pardtm.o2.T_TDi[i], &pardtm.do2.T_TDi[i], &pardtm.az.T_TDi[i], &pardtm.daz.T_TDi[i],
			 &pardtm.t0.T_TDi[i], &pardtm.dt0.T_TDi[i], &pardtm.tp.T_TDi[i], &pardtm.dtp.T_TDi[i]);
#endif
		 if (i != j_test)
	         {
                printf(" * * * WARNING: PB in field T_TDi " );
                printf( " * * * WARNING: incompatibility between dimension T_TDi & DTM file  " );
                printf( " * * *  i = %6i j_test = %6i ", i, j_test);
                Iok = 0; 
                goto fifteen;
              }
	}

// --------------------------------------------------
// ...... termes in  AMg
    for (i = 1; i <   Nb_AMg; i++ )
	{
        fgets( longstr1,210,infile);
#ifdef _MSC_VER
		sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_AMg[i]  ,&pardtm.dtt.T_AMg[i], 
                  &pardtm.h.T_AMg[i]   ,&pardtm.dh.T_AMg[i], &pardtm.he.T_AMg[i]  ,&pardtm.dhe.T_AMg[i], 
                  &pardtm.ox.T_AMg[i]  ,&pardtm.dox.T_AMg[i], &pardtm.az2.T_AMg[i] ,&pardtm.daz2.T_AMg[i], 
                  &pardtm.o2.T_AMg[i]  ,&pardtm.do2.T_AMg[i], &pardtm.az.T_AMg[i]  ,&pardtm.daz.T_AMg[i], 
                  &pardtm.t0.T_AMg[i]  ,&pardtm.dt0.T_AMg[i], &pardtm.tp.T_AMg[i]  ,&pardtm.dtp.T_AMg[i] );
#else
		sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			&j_test, &pardtm.tt.T_AMg[i], &pardtm.dtt.T_AMg[i],
			&pardtm.h.T_AMg[i], &pardtm.dh.T_AMg[i], &pardtm.he.T_AMg[i], &pardtm.dhe.T_AMg[i],
			&pardtm.ox.T_AMg[i], &pardtm.dox.T_AMg[i], &pardtm.az2.T_AMg[i], &pardtm.daz2.T_AMg[i],
			&pardtm.o2.T_AMg[i], &pardtm.do2.T_AMg[i], &pardtm.az.T_AMg[i], &pardtm.daz.T_AMg[i],
			&pardtm.t0.T_AMg[i], &pardtm.dt0.T_AMg[i], &pardtm.tp.T_AMg[i], &pardtm.dtp.T_AMg[i]);
#endif
		if (i != j_test)
	         {
                printf(" * * * WARNING: PB in field T_AMg \n" );
                printf(" * * * WARNING: incompatibility beetween dimension T_AMg & DTM file  \n" );
                printf( " * * *  i = %6i j_test = %6i \n", i, j_test);
                Iok = 0; 
                goto fifteen;
               }
	}

// --------------------------------------------------
// ...... termes in  Lon
    for (i = 1; i <   Nb_Lon; i++ )
	{
        fgets( longstr1,210,infile);
#ifdef _MSC_VER
		sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_Lon[i]  ,&pardtm.dtt.T_Lon[i], 
                  &pardtm.h.T_Lon[i]   ,&pardtm.dh.T_Lon[i], &pardtm.he.T_Lon[i]  ,&pardtm.dhe.T_Lon[i], 
                  &pardtm.ox.T_Lon[i]  ,&pardtm.dox.T_Lon[i], &pardtm.az2.T_Lon[i] ,&pardtm.daz2.T_Lon[i], 
                  &pardtm.o2.T_Lon[i]  ,&pardtm.do2.T_Lon[i], &pardtm.az.T_Lon[i]  ,&pardtm.daz.T_Lon[i], 
                  &pardtm.t0.T_Lon[i]  ,&pardtm.dt0.T_Lon[i], &pardtm.tp.T_Lon[i]  ,&pardtm.dtp.T_Lon[i] );
#else
		sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			&j_test, &pardtm.tt.T_Lon[i], &pardtm.dtt.T_Lon[i],
			&pardtm.h.T_Lon[i], &pardtm.dh.T_Lon[i], &pardtm.he.T_Lon[i], &pardtm.dhe.T_Lon[i],
			&pardtm.ox.T_Lon[i], &pardtm.dox.T_Lon[i], &pardtm.az2.T_Lon[i], &pardtm.daz2.T_Lon[i],
			&pardtm.o2.T_Lon[i], &pardtm.do2.T_Lon[i], &pardtm.az.T_Lon[i], &pardtm.daz.T_Lon[i],
			&pardtm.t0.T_Lon[i], &pardtm.dt0.T_Lon[i], &pardtm.tp.T_Lon[i], &pardtm.dtp.T_Lon[i]);
#endif
		if (i != j_test)
	         {
                printf(" *** WARNING: PB in field T_Lon \n" );
                printf(" *** WARNING: incompatibility beetween dimension T_Lon & DTM file  \n" );
                printf( " ***  i = %6i j_test = %6i \n", i, j_test);
                Iok = 0; 
                goto fifteen;
               }
	}
//
// ...... termes in dPhas 
    for (i = 1; i <   Nb_dPhas; i++ )
	{
        fgets( longstr1,210,infile);
#ifdef _MSC_VER
		sscanf_s(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			      &j_test,  &pardtm.tt.T_dPhas[i]  ,&pardtm.dtt.T_dPhas[i], 
                  &pardtm.h.T_dPhas[i]   ,&pardtm.dh.T_dPhas[i], &pardtm.he.T_dPhas[i]  ,&pardtm.dhe.T_dPhas[i], 
                  &pardtm.ox.T_dPhas[i]  ,&pardtm.dox.T_dPhas[i], &pardtm.az2.T_dPhas[i] ,&pardtm.daz2.T_dPhas[i], 
                  &pardtm.o2.T_dPhas[i]  ,&pardtm.do2.T_dPhas[i], &pardtm.az.T_dPhas[i]  ,&pardtm.daz.T_dPhas[i], 
                  &pardtm.t0.T_dPhas[i]  ,&pardtm.dt0.T_dPhas[i], &pardtm.tp.T_dPhas[i]  ,&pardtm.dtp.T_dPhas[i] );
#else
		sscanf(longstr1, " %4d %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf %13lf %9lf",
			&j_test, &pardtm.tt.T_dPhas[i], &pardtm.dtt.T_dPhas[i],
			&pardtm.h.T_dPhas[i], &pardtm.dh.T_dPhas[i], &pardtm.he.T_dPhas[i], &pardtm.dhe.T_dPhas[i],
			&pardtm.ox.T_dPhas[i], &pardtm.dox.T_dPhas[i], &pardtm.az2.T_dPhas[i], &pardtm.daz2.T_dPhas[i],
			&pardtm.o2.T_dPhas[i], &pardtm.do2.T_dPhas[i], &pardtm.az.T_dPhas[i], &pardtm.daz.T_dPhas[i],
			&pardtm.t0.T_dPhas[i], &pardtm.dt0.T_dPhas[i], &pardtm.tp.T_dPhas[i], &pardtm.dtp.T_dPhas[i]);
#endif
		if (i != j_test)
	        {
              printf( " * * * WARNING: PB in field T_dPhas " );
              printf( " * * * WARNING: incompatibility between dimension T_dPhas & DTM file  " );
              printf( " * * *  i = %6i j_test = %6i ", i, j_test);
              Iok = 0; 
              goto fifteen;
            }
	}
// --------------------------------------------------

 
// --------------------------------------------------

// * * * * * * * * * * * * Write / Verify * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * 

 if (IVERIF == 1) 
 {
// ...... termes ONE

     i = 1;
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
			      i, pardtm.tt.T_one  ,pardtm.dtt.T_one, 
                  pardtm.h.T_one   ,pardtm.dh.T_one, pardtm.he.T_one  , pardtm.dhe.T_one, 
                  pardtm.ox.T_one  ,pardtm.dox.T_one, pardtm.az2.T_one , pardtm.daz2.T_one, 
                  pardtm.o2.T_one  ,pardtm.do2.T_one, pardtm.az.T_one  , pardtm.daz.T_one, 
                  pardtm.t0.T_one  ,pardtm.dt0.T_one, pardtm.tp.T_one  , pardtm.dtp.T_one );

// ...... termes in LAT
	 // cdav lat
    for (i = 1; i <   Nb_lat; i++ )
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
			 	  i,   pardtm.tt.T_lat[i]  ,pardtm.dtt.T_lat[i], pardtm.h.T_lat[i]   ,pardtm.dh.T_lat[i], 
                  pardtm.he.T_lat[i]  ,pardtm.dhe.T_lat[i],  pardtm.ox.T_lat[i]  ,pardtm.dox.T_lat[i], 
                  pardtm.az2.T_lat[i] ,pardtm.daz2.T_lat[i], pardtm.o2.T_lat[i]  ,pardtm.do2.T_lat[i], 
                  pardtm.az.T_lat[i]  ,pardtm.daz.T_lat[i], pardtm.t0.T_lat[i]  ,pardtm.dt0.T_lat[i], 
                  pardtm.tp.T_lat[i]  ,pardtm.dtp.T_lat[i] );
	}
// ...... termes in Flux
    for (i = 1; i <   Nb_flux; i++ )   
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
			  	  i,   pardtm.tt.T_flux[i]  ,pardtm.dtt.T_flux[i], pardtm.h.T_flux[i]   ,pardtm.dh.T_flux[i], 
                  pardtm.he.T_flux[i]  ,pardtm.dhe.T_flux[i], pardtm.ox.T_flux[i]  ,pardtm.dox.T_flux[i], 
                  pardtm.az2.T_flux[i] ,pardtm.daz2.T_flux[i], pardtm.o2.T_flux[i]  ,pardtm.do2.T_flux[i], 
                  pardtm.az.T_flux[i]  ,pardtm.daz.T_flux[i], pardtm.t0.T_flux[i]  ,pardtm.dt0.T_flux[i], 
                  pardtm.tp.T_flux[i]  ,pardtm.dtp.T_flux[i] );
	}

// ...... termes in kp
    for (i = 1; i <   Nb_kp; i++ )   
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
				 i,  pardtm.tt.T_kp[i]  ,pardtm.dtt.T_kp[i], pardtm.h.T_kp[i]   ,pardtm.dh.T_kp[i], 
                  pardtm.he.T_kp[i]  ,pardtm.dhe.T_kp[i], pardtm.ox.T_kp[i]  ,pardtm.dox.T_kp[i], 
                  pardtm.az2.T_kp[i] ,pardtm.daz2.T_kp[i], pardtm.o2.T_kp[i]  ,pardtm.do2.T_kp[i], 
                  pardtm.az.T_kp[i]  ,pardtm.daz.T_kp[i], pardtm.t0.T_kp[i]  ,pardtm.dt0.T_kp[i], 
                  pardtm.tp.T_kp[i]  ,pardtm.dtp.T_kp[i] );
	}
// ...... termes in 
    for (i = 1; i <   Nb_SLat; i++ )   
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
				 i,  pardtm.tt.T_SLat[i]  ,pardtm.dtt.T_SLat[i], pardtm.h.T_SLat[i]   ,pardtm.dh.T_SLat[i], 
                  pardtm.he.T_SLat[i]  ,pardtm.dhe.T_SLat[i], pardtm.ox.T_SLat[i]  ,pardtm.dox.T_SLat[i], 
                  pardtm.az2.T_SLat[i] ,pardtm.daz2.T_SLat[i], pardtm.o2.T_SLat[i]  ,pardtm.do2.T_SLat[i], 
                  pardtm.az.T_SLat[i]  ,pardtm.daz.T_SLat[i], pardtm.t0.T_SLat[i]  ,pardtm.dt0.T_SLat[i], 
                  pardtm.tp.T_SLat[i]  ,pardtm.dtp.T_SLat[i] );
	}
// ...... termes in 
    for (i = 1; i <   Nb_SASLat; i++ )
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
				  i,  pardtm.tt.T_SASLat[i]  ,pardtm.dtt.T_SASLat[i], pardtm.h.T_SASLat[i]   ,pardtm.dh.T_SASLat[i], 
                  pardtm.he.T_SASLat[i]  ,pardtm.dhe.T_SASLat[i], pardtm.ox.T_SASLat[i]  ,pardtm.dox.T_SASLat[i], 
                  pardtm.az2.T_SASLat[i] ,pardtm.daz2.T_SASLat[i], pardtm.o2.T_SASLat[i]  ,pardtm.do2.T_SASLat[i], 
                  pardtm.az.T_SASLat[i]  ,pardtm.daz.T_SASLat[i], pardtm.t0.T_SASLat[i]  ,pardtm.dt0.T_SASLat[i], 
                  pardtm.tp.T_SASLat[i]  ,pardtm.dtp.T_SASLat[i] );
	}

// ...... termes in 
    for (i = 1; i <   Nb_NSLat; i++ )
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
         		  i,  pardtm.tt.T_NSLat[i]  ,pardtm.dtt.T_NSLat[i], pardtm.h.T_NSLat[i]   ,pardtm.dh.T_NSLat[i], 
                  pardtm.he.T_NSLat[i]  ,pardtm.dhe.T_NSLat[i], pardtm.ox.T_NSLat[i]  ,pardtm.dox.T_NSLat[i], 
                  pardtm.az2.T_NSLat[i] ,pardtm.daz2.T_NSLat[i], pardtm.o2.T_NSLat[i]  ,pardtm.do2.T_NSLat[i], 
                  pardtm.az.T_NSLat[i]  ,pardtm.daz.T_NSLat[i], pardtm.t0.T_NSLat[i]  ,pardtm.dt0.T_NSLat[i], 
                  pardtm.tp.T_NSLat[i]  ,pardtm.dtp.T_NSLat[i] );
	}

// ...... termes in 
    for (i = 1; i <   Nb_SANSLat; i++ )   
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
		 	 	  i,  pardtm.tt.T_SANSLat[i]  ,pardtm.dtt.T_SANSLat[i], pardtm.h.T_SANSLat[i]   ,pardtm.dh.T_SANSLat[i], 
                  pardtm.he.T_SANSLat[i]  ,pardtm.dhe.T_SANSLat[i], pardtm.ox.T_SANSLat[i]  ,pardtm.dox.T_SANSLat[i], 
                  pardtm.az2.T_SANSLat[i] ,pardtm.daz2.T_SANSLat[i], pardtm.o2.T_SANSLat[i]  ,pardtm.do2.T_SANSLat[i], 
                  pardtm.az.T_SANSLat[i]  ,pardtm.daz.T_SANSLat[i], pardtm.t0.T_SANSLat[i]  ,pardtm.dt0.T_SANSLat[i], 
                  pardtm.tp.T_SANSLat[i]  ,pardtm.dtp.T_SANSLat[i] );
	}
// ...... termes in 
    for (i = 1; i <   Nb_DiAn; i++ )   
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
		 		  i,  pardtm.tt.T_DiAn[i]  ,pardtm.dtt.T_DiAn[i], pardtm.h.T_DiAn[i]   ,pardtm.dh.T_DiAn[i], 
                  pardtm.he.T_DiAn[i]  ,pardtm.dhe.T_DiAn[i], pardtm.ox.T_DiAn[i]  ,pardtm.dox.T_DiAn[i], 
                  pardtm.az2.T_DiAn[i] ,pardtm.daz2.T_DiAn[i], pardtm.o2.T_DiAn[i]  ,pardtm.do2.T_DiAn[i], 
                  pardtm.az.T_DiAn[i]  ,pardtm.daz.T_DiAn[i], pardtm.t0.T_DiAn[i]  ,pardtm.dt0.T_DiAn[i], 
                  pardtm.tp.T_DiAn[i]  ,pardtm.dtp.T_DiAn[i] );
	}
// ...... termes in 
    for (i = 1; i <   Nb_SDiAn; i++ )
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
		 		  i,  pardtm.tt.T_SDiAn[i]  ,pardtm.dtt.T_SDiAn[i], pardtm.h.T_SDiAn[i]   ,pardtm.dh.T_SDiAn[i], 
                  pardtm.he.T_SDiAn[i]  ,pardtm.dhe.T_SDiAn[i], pardtm.ox.T_SDiAn[i]  ,pardtm.dox.T_SDiAn[i], 
                  pardtm.az2.T_SDiAn[i] ,pardtm.daz2.T_SDiAn[i], pardtm.o2.T_SDiAn[i]  ,pardtm.do2.T_SDiAn[i], 
                  pardtm.az.T_SDiAn[i]  ,pardtm.daz.T_SDiAn[i], pardtm.t0.T_SDiAn[i]  ,pardtm.dt0.T_SDiAn[i], 
                  pardtm.tp.T_SDiAn[i]  ,pardtm.dtp.T_SDiAn[i] );
	}
// ...... termes in 
    for (i = 1; i <   Nb_TDi; i++ )
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
		 		  i,  pardtm.tt.T_TDi[i]  ,pardtm.dtt.T_TDi[i], pardtm.h.T_TDi[i]   ,pardtm.dh.T_TDi[i], 
                  pardtm.he.T_TDi[i]  ,pardtm.dhe.T_TDi[i], pardtm.ox.T_TDi[i]  ,pardtm.dox.T_TDi[i], 
                  pardtm.az2.T_TDi[i] ,pardtm.daz2.T_TDi[i], pardtm.o2.T_TDi[i]  ,pardtm.do2.T_TDi[i], 
                  pardtm.az.T_TDi[i]  ,pardtm.daz.T_TDi[i], pardtm.t0.T_TDi[i]  ,pardtm.dt0.T_TDi[i], 
                  pardtm.tp.T_TDi[i]  ,pardtm.dtp.T_TDi[i] );
   }

// ...... termes in 
    for (i = 1; i <   Nb_AMg; i++ )
	{
     fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
		   		  i,  pardtm.tt.T_AMg[i]  ,pardtm.dtt.T_AMg[i], pardtm.h.T_AMg[i]   ,pardtm.dh.T_AMg[i], 
                  pardtm.he.T_AMg[i]  ,pardtm.dhe.T_AMg[i], pardtm.ox.T_AMg[i]  ,pardtm.dox.T_AMg[i], 
                  pardtm.az2.T_AMg[i] ,pardtm.daz2.T_AMg[i], pardtm.o2.T_AMg[i]  ,pardtm.do2.T_AMg[i], 
                  pardtm.az.T_AMg[i]  ,pardtm.daz.T_AMg[i], pardtm.t0.T_AMg[i]  ,pardtm.dt0.T_AMg[i], 
                  pardtm.tp.T_AMg[i]  ,pardtm.dtp.T_AMg[i] );
   }

// ...... termes in 
    for (i = 1; i <   Nb_Lon; i++ )   
	{
    fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
				  i,  pardtm.tt.T_Lon[i]  ,pardtm.dtt.T_Lon[i], pardtm.h.T_Lon[i]   ,pardtm.dh.T_Lon[i], 
                  pardtm.he.T_Lon[i]  ,pardtm.dhe.T_Lon[i], pardtm.ox.T_Lon[i]  ,pardtm.dox.T_Lon[i], 
                  pardtm.az2.T_Lon[i] ,pardtm.daz2.T_Lon[i], pardtm.o2.T_Lon[i]  ,pardtm.do2.T_Lon[i], 
                  pardtm.az.T_Lon[i]  ,pardtm.daz.T_Lon[i], pardtm.t0.T_Lon[i]  ,pardtm.dt0.T_Lon[i], 
                  pardtm.tp.T_Lon[i]  ,pardtm.dtp.T_Lon[i]);
   }
// ...... termes in 
	// cdav phas
    for (i = 1; i <   Nb_dPhas; i++ )
	{
		 fprintf(outfile,"%4i %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e %13.6e %9.2e \n", 
			      i,  pardtm.tt.T_dPhas[i]  ,pardtm.dtt.T_dPhas[i], pardtm.h.T_dPhas[i]   ,pardtm.dh.T_dPhas[i], 
                  pardtm.he.T_dPhas[i]  ,pardtm.dhe.T_dPhas[i], pardtm.ox.T_dPhas[i]  ,pardtm.dox.T_dPhas[i], 
                  pardtm.az2.T_dPhas[i] ,pardtm.daz2.T_dPhas[i], pardtm.o2.T_dPhas[i]  ,pardtm.do2.T_dPhas[i], 
                  pardtm.az.T_dPhas[i]  ,pardtm.daz.T_dPhas[i], pardtm.t0.T_dPhas[i]  ,pardtm.dt0.T_dPhas[i], 
                  pardtm.tp.T_dPhas[i]  ,pardtm.dtp.T_dPhas[i]);
   }

 }

fclose(infile);
if (IVERIF == 1)
    if (outfile != NULL)
        fclose(outfile);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * 
// -- First Verif /  model  >>>> ok
//       Do i = 1,Nb_lat
//      write( * , * ) 'ox(',i,') : ',ox.T_Lat[i]
//      write( * , * ) '   dox(',i,') : ',&pardtm.dox.T_Lat[i]
//       EndDo
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * 

fifteen:   

    if (Iok == 0) 
      {
//       stop;
       printf( " exit: LecDTM" );
      }

} // function InitDTM



 //---------------------------------------------------------------------------
 //
 // ROUTINE: dtm_wrapper
 //
 //> @author Raul Dominguez
 //>
 //> @brief  compute the required arguments for dtm2012 subroutine, call it
 //>              and return the results.
 //
 // INPUT ARGUMENTS:
 //>                    @param[in] latgd         latitude in degrees
 //>                    @param[in] lon          longitude in degrees
 //>                    @param[in] hellp        altitude in km (must be greater than 120)
 //>                    @param[in] in_date      UT date in MJD2000 or calendar date
 //
 // OUTPUT ARGUMENTS:
 //>                    @param[out] ro           density (g/cm^3) at the given position
 //>                    @param[out] ro_unc       density uncertainty
 //>                    @param[out] tinf         exospheric temperature
 //>                    @param[out] tp120
 //>                    @param[out] tz           temperature at the given height
 //>                    @param[out] xmm          mean molecular mass
 //>                    @param[out] d[1]         concentration in atomic hydrogen
 //>                    @param[out] d[2]         concentration in helium
 //>                    @param[out] d[3]         concentration in atomic oxygen
 //>                    @param[out] d[4]         concentration in molecular nitrogen
 //>                    @param[out] d(5)         concentration in molecular oxygen
 //>                    @param[out] d(6)         concentration in atomic nitrogen
 //
 // REVISION HISTORY:
 //
 // version 1.0
 // date 31/07/2012
 // Initial version
 //
 // version 1.1
 // date 10/08/2012
 // Added density uncertainty computation
 //
 // version 1.2
 // date 10/08/2012
 // Changed internal variable dayofyear from int to real to match its kind
 // in DTM2012 subroutine
 // Changed date input, to allow the user to use either MJD2000 or calendar date
 //
 // version 1.3
 // date 11/09/2012
 // Bugfix. The density uncertainty was being called with local solar time in
 // radians, and should be hours
 //
 // version 1.4
 // date 14/09/2012
 // Change. Input parameters are now real (before they were double)
 //
 // version 1.6
 // date 19/09/2012
 // Change in the solar proxy file format (email from Sean Bruinsma at 19/09/2012)
 // Fixed a bug: The am values where not being read properly
 //
 // @version 1.7
 // @date 28/09/2012
 // Bugfix. The interpolated akp value was not being computed properly. In order to
 // interpolate, the minutes and seconds of the user input hour were not being
 // considered
 //---------------------------------------------------------------------------

void dtm_wrapper
    (
      dtm_daterectype& in_date, double lon, int number_of_lines_in_a_file, int number_of_lines_in_f_file, 
	  a_structtype a_indexes[3000], f_structtype f_indexes[3000],
      double f[3], double fbar[3], double akp[5], double& dayofyear, double& hl
    )
{
      //User-entered date is stored in these vars
      int usr_day, usr_month, usr_year, usr_hour, usr_min;
      double usr_sec;
      double usr_hour_decimal;

      //Auxiliaries for akp computation
      // cdav increment by 1 to have the same indices the same between languages
      int a_row_current[9], a_row_former[9], a_row_formerformer[9];
      int a_row_aux[11];
      int a_row_aux_idx[11];
      int a_before, a_after;
      int day_index;

	  int idx, idx2;
	  int const max_lines = 3000;  // 3000 with smaller file 15818 19784
	  // cdav add counter i
	  long int i;
	  //save first_run,a_indexes,f_indexes, number_of_lines_in_a_file, number_of_lines_in_f_file;

	  //======================================================================
      //Process the user input date
      process_input_date (in_date);
      usr_day=in_date.day;
      usr_month=in_date.month;
      usr_year=in_date.year;
      usr_hour=in_date.hour;
      usr_min=in_date.minute;
      usr_sec=in_date.second;

      //Get the solar local time and the day of year
      hl = local_solar_time(lon,usr_hour,usr_min,usr_sec);

	  dayofyear = day_of_year(usr_day,usr_month,usr_year) + (usr_hour/24.0) +                                  
                             (usr_min/(24.0*60.0)) + (usr_sec/(24.0*3600.0));

      // Accrding to dtm2012 documentation: f[2]=0 fbar[2]=0
      // f[1] y fbar[1] are read from f_indexes structure
 //
 // cdav this is somewhat inefficient to find the solar flux and geomagnetic values
 // better to use a technique to find the recnum directly
 //          jd          = jd + mfme/1440.0;
 //          jdspwstarto = floor( jd - jdspwstart);
 //          recnum      = int(jdspwstarto);
 //

      f[2] = -1000.0;
      for (idx=1; idx < number_of_lines_in_f_file; idx++)  // < or <=?? cdav
		{
            if (f_indexes[idx].year == usr_year) 
		      {
                if (f_indexes[idx].month == usr_month) 
			      {
                    if (f_indexes[idx].day == usr_day) 
				      {
                        // f[1] must be the f index corresponding to a day before the user input date. Therefore,
                        // if the first line of the "f" file corresponds to the user input date, there is an error
						  
                        if (idx == 1) 
                            fatal_error(1);
                          else if (idx == number_of_lines_in_f_file) 
                            fatal_error(1);
						
						f[2] = 0.0; 
                        fbar[2] = 0.0;

                        //get f[1] by linear interpolation between f(day-1) and f(day)
                        //get fbar[1] by linear interpolation between f(day) and f(day+1)
						// cdav should be 2086, 118.7 98.0 fbar
                        f[1] = linear_interpolation(0.0, f_indexes[idx-1].f, 24.0, f_indexes[idx].f, hms2hr(usr_hour, usr_min, usr_sec));

                        //if the current line is the last of the file, no interpolation is posible. If this happens, do not
                        //interpolate and issue a warning
                        if (idx != number_of_lines_in_f_file) 
                            fbar[1] = linear_interpolation(0.0, f_indexes[idx].fbar, 24.0, f_indexes[idx+1].fbar, hms2hr(usr_hour, usr_min, usr_sec));
                          else
					      {
                             warning(2);
                             fbar[1] = f_indexes[idx].fbar;
					      }
                 //   exit;
				   }  // if
			   } // if
		   } // if
	    }   //  for

      // Check that the loop above has actually read something. If not, the user input date is not in the data file
      if (f[2] == -1000.0) 
           fatal_error(1);

	  //Find AKP from "a" structure
      //
      //According to the comments in DTM2012:
      // akp[2] = akp[4] = 0
      // akp[1] = kp delayed by 3 hours
      // akp[3] = mean of last 24 hours

      //Also take the following into account
      //    - Kp only takes discrete values: 0, 1/3, 2/3, 1, 4/3....
      //    - Averages must be computed from "a", not from Kp
      //    - For the same reason, interpolations are performed on "a"

      akp[2] = -1000.0;
      for (idx=1; idx <= number_of_lines_in_a_file; idx++ )
	  {
         if (a_indexes[idx].year == usr_year) 
		 {
             if (a_indexes[idx].month == usr_month) 
			 {
                 if (a_indexes[idx].day == usr_day) 
				 {
                    usr_hour_decimal = hms2hr(usr_hour,usr_min,usr_sec);

                    //Each row contains 8 values for a day: from 0h to 3h, from 3h to 6h, and so..
                    //So once the day matches, find the "a" corresponding to the input hour

                    if (usr_hour_decimal <= 3.0) 
                        day_index=1;
                    else if (usr_hour_decimal <= 6.0) 
                        day_index=2;
                    else if (usr_hour_decimal <= 9.0) 
                        day_index=3;
                    else if (usr_hour_decimal <= 12.0) 
                        day_index=4;
                    else if (usr_hour_decimal <= 15.0) 
                        day_index=5;
                    else if (usr_hour_decimal <= 18.0) 
                        day_index=6;
                    else if (usr_hour_decimal <= 21.0) 
                        day_index=7;
                    else if (usr_hour_decimal <= 24.0) 
                        day_index=8;
                    else
                        day_index=-1; //Prevent a warning.
					
                    // Get the row corresponding to the current day and the two days before

                    // If the user input date matches the first or second lines, it is not possible
                    // to compute akp[3]
                    if (idx == 1) 
					    fatal_error(9);
                    else if (idx == 2) 
                        fatal_error(9);
					else
					{
                        for (i = 1; i < 9; i++ )  // cdav set each of the geomagnetic values
						{
					    	a_row_formerformer[i] = a_indexes[idx-2].am[i];
                            a_row_former[i] = a_indexes[idx-1].am[i];
                            a_row_current[i] = a_indexes[idx].am[i];
						}
					}

                    // Write akp values
                    akp[2]=0.0;
                    akp[4]=0.0;

                    //for akp[1], linear interpolation between the nearest "a" values, which are
                    //akp_before and akp_after

                    if (day_index-1 == 0) 
					{
                      a_before = (a_row_former[8]);
                      a_after = (a_row_current[1]);
					}
                    else
					{
                      a_before=(a_row_current[day_index-1]);
                      a_after=(a_row_current[day_index]);
					} 

                    //This line interpolates a between akp_before and akp_after, and then finds akp[1] from it
                    //Note that the time is delayed 3 hours

                    akp[1] = a2K(linear_interpolation(0.0, a_before, 3.0, a_after, usr_hour_decimal - 3*(day_index-1) ));

                    //akp[3] is the average of the last 24 hours. Usually, I will need values from both
                    // a_row_current and a_row_former and even a_row_formerformer

                    //Load all the required am values required for akp[3]
                    //into an auxiliary array. The times (3,6,9...) matching
                    //those values are stores in another auxiliary array

                    //Then compute the mean am value using those arrays

                    //To compute mean am, we assume that am is a piecewise function,
                    //each of whose pieces is a straight line

                    //the auxiliary arrays are filled backwards
                    for (idx2=1; idx2 <= 10; idx2++)
					{

                        if ((day_index+1-idx2) >= 1)  //get the value from the current line
			    		{
                           a_row_aux[11-idx2] = a_row_current[day_index+1-idx2];
                           a_row_aux_idx[11-idx2] = 3*(day_index+1-idx2);
			    		}
                        else if ((day_index+1-idx2) >= -7)  //get the value from the previous line
			    		{
                           a_row_aux[11-idx2] = a_row_former[(day_index+1-idx2)+8];
                           a_row_aux_idx[11-idx2] = 3*(day_index+1-idx2);
			    		}
                        else if ((day_index+1-idx2) >= -8)  //get the value from the previous line
				    	{
                           a_row_aux[11-idx2] = a_row_formerformer[(day_index+1-idx2)+16];
                           a_row_aux_idx[11-idx2] = 3*(day_index+1-idx2);
			    		}

					}  // for

                    akp[3] = 0.0;

                    //Integrate from T-24 to the next known point, then between known
                    //points and finally from the last known point before T up to T

                    //First integral
                    akp[3] = akp[3] + (a_row_aux_idx[2]-(usr_hour_decimal-24.0))*0.5* 
                                    (a_row_aux[2]+linear_interpolation(  0.0,            
                                                                     a_row_aux[1], 
                                                                         3.0,            
                                                                     a_row_aux[2], 
                                      (usr_hour_decimal-24.0) -  a_row_aux_idx[1]  
                                    ));

                    //Seven integrals
                    for (idx2 = 1; idx2 <= 7; idx2++)
					{
                     akp[3] = akp[3]+((3.0/2.0)*(a_row_aux[idx2+1] + a_row_aux[idx2+2]));
					}

                    //Last integral
                    akp[3] = akp[3] + (((usr_hour_decimal-a_row_aux_idx[9])*0.5 * 
                             (a_row_aux[9] + 
        				     linear_interpolation (0.0, a_row_aux[9], 3.0, a_row_aux[10], usr_hour_decimal - a_row_aux_idx[9]))));


                    //The akp[3] value is the sum found previously divided by 24, and
                    // passed to the a2K function
                    akp[3] = a2K(akp[3]/24.0);

                  //  exit;

				 }
			 }
		 }
      
      } // for

      //Check that the loop above has actually read something. If not, the user input date is not in the data file
      if (akp[2] == -1000.0) 
         fatal_error(2);

}  // function dtm_wrapper



 //---------------------------------------------------------------------------
 //
 // SUBROUTINE: fatal error
 //
 //> @brief this routine comes into focus when a fatal error happens in
 //>              dtm_wrap module. The user should modify this routine
 //>              to fit into his/her error handling scheme
 //
 //> @author Raul Dominguez
 //
 // INPUT ARGUMENTS:   int code         => error code
 //
 // REVISION HISTORY: 31/07/2012 Initial version
 //
 //---------------------------------------------------------------------------
 
void fatal_error 
    (
	    int code
	)
{
    typedef char ls[120];
	ls error_messages[20];

    //these messages tell the user what happened
#ifdef _MSC_VER
	strcpy_s(error_messages[1], "Cannot read 'f' from fluxes array. Dates in file do not match user input date" );
    strcpy_s(error_messages[2] ,"Cannot read 'a' from geomagnetic data array. Dates in file do not match user input date");
    strcpy_s(error_messages[3] ,"Cannot open solar proxies file");
    strcpy_s(error_messages[4] ,"I/O error while reading solar proxies file");
    strcpy_s(error_messages[5] ,"Solar proxies file longer than max_lines parameter");
    strcpy_s(error_messages[6] ,"Cannot open geomagnetic proxies file");
    strcpy_s(error_messages[7] ,"I/O error while reading geomagnetic proxies file");
    strcpy_s(error_messages[8] ,"Geomagnetic proxies file longer than max_lines parameter");
    strcpy_s(error_messages[9] ,"Geomagnetic proxies file too new. The first record of this file must be 48 hrs or more before the date");
    strcpy_s(error_messages[10],"The provided calendar date is invalid. Please check it");
    strcpy_s(error_messages[11],"Geomagnetic proxy (am) file contains an invalid measurement. Check geomagnetic proxies file");
    strcpy_s(error_messages[12],"Geomagnetic proxy (am) has a strange value (too big). Check geomagnetic proxies file");
#else
	strcpy(error_messages[1], "Cannot read 'f' from fluxes array. Dates in file do not match user input date");
	strcpy(error_messages[2], "Cannot read 'a' from geomagnetic data array. Dates in file do not match user input date");
	strcpy(error_messages[3], "Cannot open solar proxies file");
	strcpy(error_messages[4], "I/O error while reading solar proxies file");
	strcpy(error_messages[5], "Solar proxies file longer than max_lines parameter");
	strcpy(error_messages[6], "Cannot open geomagnetic proxies file");
	strcpy(error_messages[7], "I/O error while reading geomagnetic proxies file");
	strcpy(error_messages[8], "Geomagnetic proxies file longer than max_lines parameter");
	strcpy(error_messages[9], "Geomagnetic proxies file too new. The first record of this file must be 48 hrs or more before the date");
	strcpy(error_messages[10], "The provided calendar date is invalid. Please check it");
	strcpy(error_messages[11], "Geomagnetic proxy (am) file contains an invalid measurement. Check geomagnetic proxies file");
	strcpy(error_messages[12], "Geomagnetic proxy (am) has a strange value (too big). Check geomagnetic proxies file");
#endif
    printf( "FATAL ERROR in dtm_wrap module %i \n", error_messages[code] );
} // end fatal_error


}   // namespace
 
