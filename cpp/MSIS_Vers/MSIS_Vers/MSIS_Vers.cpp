/*      ----------------------------------------------------------------
*
*                              MSIS_Vers.cpp
*
*  this file contains routines for the msis-86, msis-90, and nrlmsise-00 atmospheric models.
*
*                          companion code for
*             fundamentals of astrodynamics and applications
*                                  2013
*                            by david vallado
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

//#include "stdafx.h"

#include "MSIS_Vers.h"


namespace MSIS_Vers {
	// ---------------- MSIS-com ----------------
	/*    *****************************************************************
	*        set up routines that are common to all the msis models
	*     *****************************************************************       */

	//----------------------------------------------------------------------
	//     set switches
	//     output in  common/csw/sw(25),isw,swc(25)
	//     sw for main terms, swc for cross terms
	//
	//     to turn on and off particular variations call tselec(sv),
	//     where sv is a 25 element array containing 0.0 for off, 1.0
	//     for on, or 2.0 for main effects off but cross terms on
	//
	//  cdav this was deleted since one can simply examine the csw structure
	//     to get current values of sw: call tretrv(sw)
	//----------------------------------------------------------------------

	void tselec
		(
		cswtype& csw, int sv[26]
		)
	{
		//        int sav[26];
		int i;

		//------------------------------ begin ---------------------------------
		for (i = 1; i <= 25; i++)
		{
			//            sav[i] = sv[i];
			csw.sw[i] = fmodf(sv[i], 2);
			if ((fabsf(sv[i]) == 1) | (fabsf(sv[i]) == 2))
				csw.swc[i] = 1;
			else
				csw.swc[i] = 0;
		}

		csw.isw = 64999;
	}

	//----------------------------------------------------------------------
	//     turbopause correction for msis models
	//     eq. a12b
	//         root mean density
	//          dd   - diffusive density
	//          dm   - full mixed density
	//          zhm  - transition scale length
	//          xmm  - full mixed molecular weight
	//          xm   - species molecular weight
	//          dnet - combined density
	//       8/20/80
	//----------------------------------------------------------------------

	double dnet
		(
		double dd, double dm, double zhm, double xmm, double xm
		)
	{
		double a, ylog;

		//------------------------------ begin ---------------------------------
		a = zhm / (xmm - xm);

		if ((dm <= 0.0) | (dd <= 0.0))
		{
			// cdav this loop is exercised with test case # 11 and msis86. i changed the return
			// from the first if, and then used else's for the remaining ones.
			// essentially, it appears that msis86 does not properly work down to 0km alt.
			// the test cases match for 100 km, but no others were given. I expanded the
			// original test cases and this cpp version matches 30, 50, and 70 km. the msis86
			//  model gets nan's for 0 and 10 km. the answers below 100 km are all suspect.
			printf("dnet log error %11.7g %11.7g %11.7f %11.7f \n", dm, dd, xm, a);
			printf("dd and dm may need to be output vars for this case\n");

			if (dd == 0.0  &&  dm == 0.0)
				return 1.0;    // dd = 1.0
			else
			{
				if (dm == 0.0)
					return dd;
				else
					//                     if (dd == 0.0)
					return dm;
			}
		}
		else
		{
			// ---- eq. a12a ----
			ylog = a*log(dm / dd);
			if (ylog < -10.0)
				return dd;
			else
			{
				if (ylog > 10.0)
					return dm;
				else
					return dd*pow(1.0 + exp(ylog), (1.0 / a));
			}
		}

	}

	//----------------------------------------------------------------------
	//     chemistry/dissociation correction for msis models
	//        alt - altitude
	//        r   - target ratio
	//        h1  - transition scale length
	//        zh  - altitude of 1/2 r
	//     eq. a20a or eq. a21
	//----------------------------------------------------------------------

	double ccor
		(
		double alt, double r, double h1, double zh
		)
	{
		double e, ex, ccor1;
		//------------------------------ begin ---------------------------------
		e = (alt - zh) / h1;
		if (e > 70.0)
			ccor1 = 0.0;
		else
		{
			if (e < -70.0)
				ccor1 = r;
			else
			{
				ex = exp(e);
				ccor1 = r / (1.0 + ex);
			}
		}
		return exp(ccor1);
	}

	double ccor2
		(
		double alt, double r, double h1, double zh, double h2
		)
	{
		double e1, e2, ccor1;
		//------------------------------ begin ---------------------------------
		e1 = (alt - zh) / h1;
		e2 = (alt - zh) / h2;
		if ((e1 > 70.0) | (e2 > 70.0))
			ccor1 = 0.0;
		else
		{
			if (e1 < -70.0 && e2 < -70.0)
				ccor1 = r;
			else
				ccor1 = r / (1.0 + 0.5*(exp(e1) + exp(e2)));
		}
		return exp(ccor1);
	}

	//-----------------------------------------------------------------------
	//      control option of converting to kg/m3 or lb/ft3.
	//-----------------------------------------------------------------------

	void meters
		(
		metseltype& metsel
		)
	{
		if (metsel.imr == 0)
			metsel.imr = 1;
		else
			metsel.imr = 0;
	}


	/*     *****************************************************************
	*       set up routines that are common to the msis-86, 90, and 00 models
	*     **************************************************************** */

	// cdav these next few functions could be done as inline functions if that speeds
	// things up
	// ---- 3hr magnetic activity functions
	//      eq. a24d
	double g0
		(double a, double p[151])
	{
		return (a - 4.0 + (p[26] - 1.0) *
			(a - 4.0 + (exp(-fabs(p[25])*(a - 4.0)) - 1.0) / fabs(p[25])));
	}
	// ---- eq. a24c
	double sumex
		(double ex)
	{
		return 1.0 + (1.0 - pow(ex, 19)) / (1.0 - ex)*sqrt(ex);
	}
	// ---- eq. a24a
	double sg0
		(double ex, double p[151], double ap[8])
	{
		return (g0(ap[2], p) +
			(g0(ap[3], p)*ex + g0(ap[4], p)*ex*ex + g0(ap[5], p)* pow(ex, 3) +
			(g0(ap[6], p)*pow(ex, 4) + g0(ap[7], p)*pow(ex, 12)) * (1.0 - pow(ex, 8))
			/ (1.0 - ex)
			)
			) / sumex(ex);
	}

	/*     *****************************************************************
	*       set up routines that are common to the msis-90 and 00 models
	*     **************************************************************** */

	//----------------------------------------------------------------------
	//       calculate temperature and density profiles for msis models
	//       new lower thermo polynomial 10/30/89
	//
	//----------------------------------------------------------------------

	// do as in anline function/////////
	double zeta
		(
		double zz, double zl, double re
		)
	{
		return (zz - zl)*(re + zl) / (re + zz);
	}

	double densu
		(
		parmbtype& parmb, lsqvtype& lsqv,
		double alt, double dlb, double tinf, double tlb, double xm, double alpha, double& tz,
		double zlb, double s2, int mn1, double zn1[6], double tn1[6], double tgn1[3]
		)
	{
		double z, za, zg2, dta, t1, t2, densu1,
			zgdif, zg, z1, yd1, yd2, y, y2out[6], glb, gamma,
			expl, densa, gamm, tt, ta, z2, x;
		static double  yi, xs[6], ys[6];
		int k;
		static int mn;

		//        dimension zn1(mn1),tn1(mn1)

		const double rgas = 831.4;

		//       statement function
		//cdav        zeta(zz,zl) = (zz-zl)*(re+zl)/(re+zz);

		//------------------------------ begin ---------------------------------
		//       printf(6,*) 'db',alt,dlb,tinf,tlb,xm,alpha,zlb,s2,mn1,zn1,tn1
		densu1 = 1.0;

		// -------- joining altitude of bates and spline
		za = zn1[1];
		if (alt >= za)
			z = alt;
		else
			z = za;
		//cdav        z = max(alt,za);

		// -------- geopotential altitude difference from zlb
		zg2 = zeta(z, zlb, parmb.re);

		// -------- bates temperature
		tt = tinf - (tinf - tlb)*exp(-s2*zg2);
		ta = tt;
		tz = tt;
		densu1 = tz;

		// --------------- calculate temperature below za ----------------
		// -------- temperature gradient at za from bates profile
		if (alt < za)
		{
			dta = (tinf - ta)*s2*(pow((parmb.re + zlb) / (parmb.re + za), 2));
			tgn1[1] = dta;
			tn1[1] = ta;
			if (alt >= zn1[mn1])
				z = alt;
			else
				z = zn1[mn1];
			//cdav          z  = max(alt,zn1[mn1]);

			mn = mn1;
			z1 = zn1[1];
			z2 = zn1[mn];
			t1 = tn1[1];
			t2 = tn1[mn];

			// -------- geopotental difference from z1
			zg = zeta(z, z1, parmb.re);
			zgdif = zeta(z2, z1, parmb.re);

			// -------- set up spline nodes
			for (k = 1; k <= mn; k++)
			{
				xs[k] = zeta(zn1[k], z1, parmb.re) / zgdif;
				ys[k] = 1.0 / tn1[k];
			}

			// -------- end node derivatives
			yd1 = -tgn1[1] / (t1*t1)*zgdif;
			yd2 = -tgn1[2] / (t2*t2)*zgdif*pow(((parmb.re + z2) / (parmb.re + z1)), 2);

			// -------- calculate spline coefficients
			spline(xs, ys, mn, yd1, yd2, y2out);
			x = zg / zgdif;
			splint(xs, ys, y2out, mn, x, y);

			// -------- temperature at altitude
			tz = 1.0 / y;
			densu1 = tz;
			//cdav
		}

		// ------------ calculate density above za -------------------
		if (xm != 0.0)
		{
			glb = parmb.gsurf / pow((1.0 + zlb / parmb.re), 2);
			gamma = xm*glb / (s2*rgas*tinf);
			expl = exp(-s2*gamma*zg2);
			if ((expl > 50) | (tt <= 0.0))
				expl = 50.0;

			// -------- density at altitude
			densa = dlb*pow((tlb / tt), (1.0 + alpha + gamma))*expl;
			densu1 = densa;

			if (alt < za)
			{
				// ---------- calculate density below za -------------
				glb = parmb.gsurf / pow((1.0 + z1 / parmb.re), 2);
				gamm = xm*glb*zgdif / rgas;

				// -------- integrate spline temperatures
				splini(xs, ys, y2out, mn, x, yi);
				expl = gamm*yi;
				if ((expl > 50.0) | (tz <= 0.0))
					expl = 50.0;

				// -------- density at altitude
				densu1 = densu1*pow((t1 / tz), (1.0 + alpha)) *exp(-expl);
			}
		}

		return densu1;
	}

	//----------------------------------------------------------------------
	//
	//
	// calculate temperature and density profiles for lower atmos.
	//
	//----------------------------------------------------------------------

	double densm
		(
		parmbtype& parmb, fittype& fit, lsqvtype& lsqv,
		double alt, double d0, double xm, double& tz, int mn3, double zn3[6],
		double tn3[6], double tgn3[3], int mn2,
		double zn2[5], double tn2[5], double tgn2[3]
		)
	{
		double z, t1, t2, densm1,
			zgdif, zg, xs[11], z1, yd1, yd2, y, y2out[11], glb,
			ys[11], gamm, z2, x, yi, expl;
		int mn, k;

		//        dimension zn3(mn3),tn3(mn3)
		//        dimension zn2(mn2),tn2(mn2)

		const double rgas = 831.4;

		//       statement function
		//cdav        zeta(zz,zl) = (zz-zl)*(re+zl)/(re+zz);

		//------------------------------ begin ---------------------------------
		densm1 = d0;

		if (alt <= zn2[1])
		{
			// -------- stratosphere/mesosphere temperature
			if (alt >= zn2[mn2])
				z = alt;
			else
				z = zn2[mn2];
			mn = mn2;
			z1 = zn2[1];
			z2 = zn2[mn];
			t1 = tn2[1];
			t2 = tn2[mn];
			zg = zeta(z, z1, parmb.re);
			zgdif = zeta(z2, z1, parmb.re);

			// -------- set up spline nodes
			for (k = 1; k <= mn; k++)
			{
				xs[k] = zeta(zn2[k], z1, parmb.re) / zgdif;
				ys[k] = 1.0 / tn2[k];
			}

			yd1 = -tgn2[1] / (t1*t1)*zgdif;
			yd2 = -tgn2[2] / (t2*t2)*zgdif*pow((parmb.re + z2) / (parmb.re + z1), 2);

			// -------- calculate spline coefficients
			spline(xs, ys, mn, yd1, yd2, y2out);
			x = zg / zgdif;
			splint(xs, ys, y2out, mn, x, y);

			// -------- temperature at altitude
			tz = 1.0 / y;

			if (xm != 0.0)
			{
				// ------ calculate stratosphere/mesosphere density ------
				glb = parmb.gsurf / pow((1.0 + z1 / parmb.re), 2);
				gamm = xm*glb*zgdif / rgas;

				// -------- integrate temperature profile
				splini(xs, ys, y2out, mn, x, yi);
				expl = gamm*yi;
				if (expl > 50.0)
					expl = 50.0;

				// -------- density at altitude
				densm1 = densm1*(t1 / tz)*exp(-expl);
			}

			if (alt <= zn3[1])
			{
				// -------- troposphere/stratosphere temperature ---------
				z = alt;
				mn = mn3;
				z1 = zn3[1];
				z2 = zn3[mn];
				t1 = tn3[1];
				t2 = tn3[mn];
				zg = zeta(z, z1, parmb.re);
				zgdif = zeta(z2, z1, parmb.re);

				// -------- set up spline nodes
				for (k = 1; k <= mn; k++)
				{
					xs[k] = zeta(zn3[k], z1, parmb.re) / zgdif;
					ys[k] = 1.0 / tn3[k];
				}
				yd1 = -tgn3[1] / (t1*t1)*zgdif;
				yd2 = -tgn3[2] / (t2*t2)*zgdif*pow((parmb.re + z2) / (parmb.re + z1), 2);

				// -------- calculate spline coefficients
				spline(xs, ys, mn, yd1, yd2, y2out);
				x = zg / zgdif;
				splint(xs, ys, y2out, mn, x, y);

				// -------- temperature at altitude
				tz = 1.0 / y;
				if (xm != 0.0)
				{
					// --- calculate tropospheric/stratosphere density ---
					glb = parmb.gsurf / pow((1.0 + z1 / parmb.re), 2);
					gamm = xm*glb*zgdif / rgas;

					// -------- integrate temperature profile
					splini(xs, ys, y2out, mn, x, yi);
					expl = gamm*yi;
					if (expl > 50.0)
						expl = 50.0;

					// -------- density at altitude
					densm1 = densm1*(t1 / tz)*exp(-expl);
				}
			}
		}

		if (xm == 0.0)
			densm1 = tz;

		return densm1;
	}

	//----------------------------------------------------------------------
	//  calculate 2nd derivatives of cubic spline interp function
	//  adapted from numerical recipes by press et al
	//  x,y    : arrays of tabulated function in ascending order by x
	//  n      : size of arrays x,y
	//  yp1,ypn: specified derivatives at x[1] and x[n]; values
	//           > =  1d30 signal signal second derivative zero
	//  y2     : output array of second derivatives
	//----------------------------------------------------------------------

	void spline
		(
		double x[6], double y[6], int n, double yp1, double ypn, double y2[6]
		)
	{
		int i, k;

		double sig, p, qn, un, u[101];
		//        parameter (nmax = 101)
		//        dimension x[n],y[n],y2[n],u(nmax)

		//------------------------------ begin ---------------------------------
		if (yp1 > 0.99e30)
		{
			y2[1] = 0.0;
			u[1] = 0.0;
		}
		else
		{
			y2[1] = -0.5;
			u[1] = (3.0 / (x[2] - x[1]))*((y[2] - y[1]) / (x[2] - x[1]) - yp1);
		}

		for (i = 2; i <= n - 1; i++)
		{
			sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
			p = sig*y2[i - 1] + 2.0;
			y2[i] = (sig - 1.0) / p;
			u[i] = (6.0*((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1])
				/ (x[i] - x[i - 1])) / (x[i + 1] - x[i - 1]) - sig*u[i - 1]) / p;
		}

		if (ypn > 0.99e30)
		{
			qn = 0.0;
			un = 0.0;
		}
		else
		{
			qn = 0.5;
			un = (3.0 / (x[n] - x[n - 1]))*(ypn - (y[n] - y[n - 1]) / (x[n] - x[n - 1]));
		}

		y2[n] = (un - qn*u[n - 1]) / (qn*y2[n - 1] + 1.0);
		for (k = n - 1; k >= 1; k--)
			y2[k] = y2[k] * y2[k + 1] + u[k];
	}

	//----------------------------------------------------------------------
	//  calculate cubic spline interp value
	//  adapted from numberical recipes by press et al.
	//  xa,ya: arrays of tabulated function in ascending order by x
	//  y2a  : array of second derivatives
	//  n    : size of arrays xa,ya,y2a
	//  x    : abscissa for interpolation
	//  y    : output value
	//----------------------------------------------------------------------

	void splint
		(
		double xa[6], double ya[6], double y2a[6], int n, double x, double& y
		)
	{
		int k, khi, klo;
		double h, a, b;
		//        dimension xa[n],ya[n],y2a[n]

		//------------------------------ begin ---------------------------------
		klo = 1;
		khi = n;
		while (khi - klo > 1)
		{
			k = (khi + klo) * 0.5;
			if (xa[k] > x)
				khi = k;
			else
				klo = k;
		}

		h = xa[khi] - xa[klo];
		if (h == 0.0)
			printf("bad xa input to splint\n");
		a = (xa[khi] - x) / h;
		b = (x - xa[klo]) / h;
		y = a*ya[klo] + b*ya[khi] +
			((pow(a, 3) - a)*y2a[klo] + (pow(b, 3) - b)*y2a[khi])*h*h / 6.0;

	}

	//----------------------------------------------------------------------
	//     integrate cubic spline function from xa[1] to x
	//  xa,ya: arrays of tabulated function in ascending order by x
	//  y2a  : array of second derivatives
	//  n    : size of arrays xa,ya,y2a
	//  x    : abscissa endpoint for integration
	//  y    : output value
	//----------------------------------------------------------------------

	void splini
		(
		double xa[6], double ya[6], double y2a[6], int n, double x, double& yi
		)
	{
		int khi, klo;
		double h, a, b, xx, a2, b2;

		//------------------------------ begin ---------------------------------
		yi = 0;
		klo = 1;
		khi = 2;

		while (x > xa[klo] && khi <= n)
		{
			xx = x;
			if (khi < n)
			{
				if (x <= xa[khi])
					xx = x;
				else
					xx = xa[khi];
				//                xx = min(x,xa[khi]);
			}

			h = xa[khi] - xa[klo];
			a = (xa[khi] - xx) / h;
			b = (xx - xa[klo]) / h;
			a2 = a*a;
			b2 = b*b;
			yi = yi + ((1.0 - a2)*ya[klo] * 0.5 + b2*ya[khi] * 0.5 +
				((-(1.0 + a2*a2)*0.25 + a2*0.5)*y2a[klo] +
				(b2*b2 / 4.0 - b2*0.5)*y2a[khi])*h*h / 6.0) * h;
			klo = klo + 1;
			khi = khi + 1;
		}

	}

	//----------------------------------------------------------------------
	//      calculate latitude variable gravity (gv) and effective
	//      radius (reff)
	//
	//dav chg to *8
	//----------------------------------------------------------------------

	void glatf
		(
		double lat, double& gv, double& reff
		)
	{
		double c2;

		const double  dgtr = 1.74533e-2;
		double  latl = -999.0;

		//------------------------------ begin ---------------------------------
		if (lat != latl)
			c2 = cos(2.0*dgtr*lat);

		latl = lat;
		gv = 980.616*(1.0 - 0.0026373*c2);
		reff = 2.0*gv / (3.085462e-6 + 2.27e-9*c2)*1.0e-5;

	}

	//----------------------------------------------------------------------
	//       test if geophysical variables or switches changed and save
	//----------------------------------------------------------------------

	double vtst
		(
		cswtype& csw,
		int iyd, double sec, double glat, double glong, double stl,
		double f107a, double f107, double ap[8], int ic
		)
	{

		int i, iydl[3];
		double secl[3], glatl[3], gll[3], stll[3], fal[3], fl[3], apl[8][3],
			swl[26][3], swcl[26][3], vtst1;
		//cdav
		i = 1;
		for (i = 1; i< 3; i++)
		{
			iydl[i] = -999;
			secl[i] = -999.0;
			glatl[i] = -999.0;
			gll[i] = -999.0;
			stll[i] = -999.0;
			fal[i] = -999.0;
			fl[i] = -999.0;
		}
		for (i = 1; i< 8; i++)
		{
			apl[i][1] = -999.0;
			apl[i][2] = -999.0;
		}
		for (i = 1; i< 26; i++)
		{
			swl[i][1] = -999.0;
			swl[i][2] = -999.0;
			swcl[i][1] = -999.0;
			swcl[i][2] = -999.0;
		}

		//------------------------------ begin ---------------------------------
		vtst1 = 0;
		if (iyd != iydl[ic])  goto ten;
		if (sec != secl[ic])  goto ten;
		if (glat != glatl[ic]) goto ten;
		if (glong != gll[ic])   goto ten;
		if (stl != stll[ic])  goto ten;
		if (f107a != fal[ic])   goto ten;
		if (f107 != fl[ic])    goto ten;

		for (i = 1; i <= 7; i++)
		if (ap[i] != apl[i][ic]) goto ten;

		for (i = 1; i <= 25; i++)
		{
			if (csw.sw[i] != swl[i][ic]) goto ten;
			if (csw.swc[i] != swcl[i][ic]) goto ten;
		}

		goto twenty;

	ten:

		vtst1 = 1;

		iydl[ic] = iyd;
		secl[ic] = sec;
		glatl[ic] = glat;
		gll[ic] = glong;
		stll[ic] = stl;
		fal[ic] = f107a;
		fl[ic] = f107;

		for (i = 1; i <= 7; i++)
			apl[i][ic] = ap[i];

		for (i = 1; i <= 25; i++)
		{
			swl[i][ic] = csw.sw[i];
			swcl[i][ic] = csw.swc[i];
		}

	twenty:

		return vtst1;

	}

	// ---------------- MSIS-86 model ----------------
	/*    ******************* local routines ******************************       */
	double denss
		(
		parmbtype& parmb, lsqvtype& lsqv,
		double alt, double dlb, double tinf, double tlb, double xm, double alpha, double& tz,
		double zlb, double s2, double t0, double za, double z0, double tr12, double& taf
		);

	double globe5
		(
		cswtype& csw,
		lpolytype& lpoly,
		double yrd, double sec, double lat, double llong, double tloc, double f107a,
		double f107, double ap[8], double p[151]
		);

	double glob5l
		(
		lpolytype& lpoly,
		cswtype& csw,
		double p[151]
		);



	/* *********************************************************************
	*        msis-86/cira 1986 neutral thermosphere model
	*         a.e.hedin 3/15/85;2/26/87 [variable \names shortened]
	*         10/14/87 increase altitude limit of o mixing calculation
	*             altl[2] from 300.0 to 400.0 km .
	*     input:
	*        iyd - year and day as yyyyddd
	*        sec - ut[sec]
	*        alt - altitude[km] [greater than 85 km]
	*        glat - geodetic latitude[deg]
	*        glong - geodetic longitude[deg]
	*        stl - local apparent solar time[hrs]
	*        f107a - 3 month [81-day] average of f10.7 flux
	*        f107 - daily f10.7 flux for previous day
	*        ap - magnetic index[daily] or when sw[9] = -1.0 :
	*           - array containing:
	*             [1] daily ap
	*             [2] 3 hr ap index for current time
	*             [3] 3 hr ap index for 3 hrs before current time
	*             [4] 3 hr ap index for 6 hrs before current time
	*             [5] 3 hr ap index for 9 hrs before current time
	*             [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
	*                    to current time
	*             [7] average of eight 3 hr ap indicies from 36 to 59 hrs prior
	*                    to current time
	*        mass - mass number [only density for selected gas is
	*                 calculated.  mass 0 is temperature.  mass 48 for all.
	*     output:
	*        d[1] - he number density[cm-3]
	*        d[2] - o number density[cm-3]
	*        d[3] - n2 number density[cm-3]
	*        d[4] - o2 number density[cm-3]
	*        d[5] - ar number density[cm-3]
	*        d[6] - total mass density[gm/cm3]
	*        d[7] - h number density[cm-3]
	*        d[8] - n number density[cm-3]
	*        t[1] - exospheric temperature
	*        t[2] - temperature at alt
	*
	*       to get output in m-3 and kg/m3:   call meters[.true.]
	*
	*          additional comments
	*           [1] lower bound quantities in common/gts3c/
	*           [2] to turn on and off particular variations call tselec[sw]
	*               where sw is a 25 element array containing 0. for off, 1.0
	*               for on, or 2.0 for main effects off but cross terms on
	*               for the following variations
	*               1 - f10.7 effect on mean  2 - time independent
	*               3 - symmetrical annual    4 - symmetrical semiannual
	*               5 - asymmetrical annual   6 - asymmetrical semiannual
	*               7 - diurnal               8 - semidiurnal
	*               9 - daily ap             10 - all ut/long effects
	*              11 - longitudinal         12 - ut and mixed ut/long
	*              13 - mixed ap/ut/long     14 - terdiurnal
	*              15 - departures from diffusive equilibrium
	*              16 - all tinf var         17 - all tlb var
	*              18 - all t0 var           19 - all s var
	*              20 - all z0 var           21 - all nlb var
	*              22 - all tr12 var         23 - turbo scale height var
	*
	*              to get current values of sw: call tretrv[sw]
	*
	*  changes
	*
	*  - name,isd,ist,isdate, and istime were changed to character variables
	*    in gts5 and prmsg5
	*
	*  - the variable dimension of p and ap in globe5 and globe5l was
	*    indicted by *, rather than 1; if this does not work on your system
	*    you may want to use p[150] and ap[7].
	*
	*  - the large data statement in prmsg5 is now read in from file
	*    msis86.dat; some compilers do not allow named commons to be
	*    initialized in a data statement.
	*
	*  - the first call to globe5 should come before the common array sw[25]
	*    is used in gts5.
	*
	* dieter bilitza march 87
	* *********************************************************************       */

	void gts5
		(
		msistype& msis86r,
		lpolytype& lpoly,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
		double f107, double ap[8], int mass, double d[10], double t[3]
		)
	{
		double yrd, tinf, tr12, gggg, taf,
			b01, b04, b14, b16, b28, b32, b40,
			g1, g14, g16, g28, g32, g4, g40,
			hc01, hc04, hc14, hc16, hc32, hc40,
			hcc01, hcc14, hcc16,
			rc01, rc14, rc16, tz, xmd, xmm,
			zc01, zc04, zc14, zc16, zc32, zc40,
			zcc01, zcc14, zcc16,
			zh01, zh04, zh14, zh16, zh28, zh32, zh40, 
			zhm01, zhm04, zhm14, zhm16, zhm28, zhm32, zhm40;
		int i, j;
		double temppd[151];
		double mt[11] = { 0, 48, 0, 4, 16, 28, 32, 40, 1, 49, 14 };
		double altl[9] = { 0, 200.0, 400.0, 150.0, 200.0, 240.0, 450.0, 320.0, 450.0 };

		//------------------------------ begin ---------------------------------
		yrd = iyd;
		// ---- eq. a7
		//////old//// tinf = msis86r.lower.ptm[1]*[1.0+msis86r.csw.sw[16]*globe5[yrd,sec,glat,glong,stl,f107a,f107,
		//////old////  ap,msis86r.parm.pt]]*msis86r.parm.pt[1]

		gggg = globe5(msis86r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap, msis86r.parm.pt);
		tinf = msis86r.lower.ptm[1] * (1.0 + msis86r.csw.sw[16] * gggg)*msis86r.parm.pt[1];
		msis86r.gts3c.za = msis86r.lower.ptm[5] * msis86r.parm.pdl[16][2];

		// ---- eq. a9
		for (j = 0; j <= 26; j++)
			temppd[j] = msis86r.parm.pddav[3][j + 75];
		msis86r.gts3c.t0 = msis86r.lower.ptm[3] * msis86r.parm.pd[76][3] *
			(1.0 + msis86r.csw.sw[18] * glob5l(lpoly, msis86r.csw, temppd));

		// ---- eq. a8
		for (j = 0; j <= 26; j++)
			temppd[j] = msis86r.parm.pddav[3][j + 25];
		msis86r.gts3c.tlb = msis86r.lower.ptm[2] * (1.0 + msis86r.csw.sw[17] *
			glob5l(lpoly, msis86r.csw, temppd))*msis86r.parm.pd[26][3];

		// ---- eq. a10
		for (j = 0; j <= 26; j++)
			temppd[j] = msis86r.parm.pddav[3][j + 50];
		msis86r.gts3c.z0 = msis86r.lower.ptm[7] * (1.0 + msis86r.csw.sw[20] *
			glob5l(lpoly, msis86r.csw, temppd))*msis86r.parm.pd[51][3];

		// ---- eq. a6
		msis86r.gts3c.g0 = msis86r.lower.ptm[4] * msis86r.parm.ps[1] * (1.0 + msis86r.csw.sw[19] *
			globe5(msis86r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap, msis86r.parm.ps));

		// ---- eq. a5
		msis86r.gts3c.s = msis86r.gts3c.g0 / (tinf - msis86r.gts3c.tlb);

		// ---- eq. a11
		for (j = 0; j <= 26; j++)
			temppd[j] = msis86r.parm.pddav[3][j + 100];
		tr12 = msis86r.parm.pd[101][3] * (1.0 + msis86r.csw.sw[22] *
			glob5l(lpoly, msis86r.csw, temppd));
		t[1] = tinf;

		if (mass == 0)
			goto fifty;

		// ---- eq. a18  n2
		for (j = 0; j <= 26; j++)
			temppd[j] = msis86r.parm.pddav[3][j + 0];
		g28 = msis86r.csw.sw[21] * glob5l(lpoly, msis86r.csw, temppd);
		yrd = iyd;
		t[1] = tinf;
		xmm = msis86r.lower.pdm[5][3];

		for (i = 1; i <= 10; i++)
		if (mass == mt[i])
			goto fifteen;
		printf("mass %i5 not valid\n", mass);
		goto ninety;

	fifteen:
		if ((alt <= altl[6]) | (mass == 28) | (mass == 48))
		{
			// ---- **** n2 density ****
			// ---- eq. a18
			msis86r.gts3c.db28 = msis86r.lower.pdm[1][3] * exp(g28)*msis86r.parm.pd[1][3];
			// ---- eq. a13 - a17
			d[3] = denss(msis86r.parmb, lsqv, alt, msis86r.gts3c.db28, tinf, msis86r.gts3c.tlb, 28.0, 0.0,
				t[2], msis86r.lower.ptm[6],
				msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);
			msis86r.gts3c.dd = d[3];
			// ---- eq. a19
			zh28 = msis86r.lower.pdm[3][3];
			zhm28 = msis86r.lower.pdm[4][3] * msis86r.parm.pdl[6][2];
			xmd = 28.0 - xmm;
			b28 = denss(msis86r.parmb, lsqv, zh28, msis86r.gts3c.db28, tinf, msis86r.gts3c.tlb, xmd, -1.0,
				tz, msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za,
				msis86r.gts3c.z0, tr12, taf);
			if (alt <= altl[3] && msis86r.csw.sw[15] != 0)
			{
				msis86r.dmix.dm28 = denss(msis86r.parmb, lsqv, alt, b28, tinf, msis86r.gts3c.tlb, xmm, 0.0,
					tz, msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0,
					msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);
				// ---- eq. a12
				d[3] = dnet(d[3], msis86r.dmix.dm28, zhm28, xmm, 28.0);
			}
		}

		switch (i)
		{
		case 1: goto twenty;
		case 2: goto fifty;
		case 3: goto twenty;
		case 4: goto twentyfive;
		case 5: goto ninety;
		case 6: goto thirtyfive;
		case 7: goto fourty;
		case 8: goto fourtyfive;
		case 9: goto twentyfive;
		case 10: goto fourtyeight;
		}
	twenty:
		// ---- **** he density ****
		// ---- eq. a18
		g4 = msis86r.csw.sw[21] * globe5(msis86r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis86r.parm.pddav[1]);
		msis86r.gts3c.db04 = msis86r.lower.pdm[1][1] * exp(g4)*msis86r.parm.pd[1][1];

		// ---- eq. a13 - a17
		d[1] = denss(msis86r.parmb, lsqv, alt, msis86r.gts3c.db04, tinf, msis86r.gts3c.tlb, 4.0, -0.4, t[2],
			msis86r.lower.ptm[6],
			msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);
		msis86r.gts3c.dd = d[1];
		if (alt <= altl[1] && msis86r.csw.sw[15] != 0)
		{
			// ---- eq. a19
			zh04 = msis86r.lower.pdm[3][1];
			b04 = denss(msis86r.parmb, lsqv, zh04, msis86r.gts3c.db04, tinf, msis86r.gts3c.tlb, 4.0 - xmm, -1.4,
				t[2], msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za,
				msis86r.gts3c.z0, tr12, taf);
			msis86r.dmix.dm04 = denss(msis86r.parmb, lsqv, alt, b04, tinf, msis86r.gts3c.tlb, xmm, 0.0, t[2],
				msis86r.lower.ptm[6], msis86r.gts3c.s,
				msis86r.gts3c.t0, msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);
			// ---- eq. a12
			zhm04 = zhm28;
			d[1] = dnet(d[1], msis86r.dmix.dm04, zhm04, xmm, 4.0);

			// ---- eq. a20b
			msis86r.gts3c.rl = log(b28*msis86r.lower.pdm[2][1] / b04);

			// ---- eq. a20a
			zc04 = msis86r.lower.pdm[5][1] * msis86r.parm.pdl[1][2];
			hc04 = msis86r.lower.pdm[6][1] * msis86r.parm.pdl[2][2];
			d[1] = d[1] * ccor(alt, msis86r.gts3c.rl, hc04, zc04);
		}

		if (mass != 48)
			goto ninety;

	twentyfive:
		// ----**** o density ****
		// ---- eq. a18
		g16 = msis86r.csw.sw[21] * globe5(msis86r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis86r.parm.pddav[2]);
		msis86r.gts3c.db16 = msis86r.lower.pdm[1][2] * exp(g16)*msis86r.parm.pd[1][2];

		// ---- eq. a13 - a17
		d[2] = denss(msis86r.parmb, lsqv, alt, msis86r.gts3c.db16, tinf, msis86r.gts3c.tlb, 16.0, 0.0, t[2],
			msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0,
			msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);
		msis86r.gts3c.dd = d[2];
		if (alt <= altl[2] && msis86r.csw.sw[15] != 0)
		{
			//  corrected from msis86r.lower.pdm[3][1] to msis86r.lower.pdm[3][2]  12/2/85
			// ---- eq. a19
			zh16 = msis86r.lower.pdm[3][2];
			b16 = denss(msis86r.parmb, lsqv, zh16, msis86r.gts3c.db16, tinf, msis86r.gts3c.tlb, 16.0 - xmm, -1.0,
				t[2], msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za,
				msis86r.gts3c.z0, tr12, taf);
			msis86r.dmix.dm16 = denss(msis86r.parmb, lsqv, alt, b16, tinf, msis86r.gts3c.tlb, xmm, 0.0, t[2],
				msis86r.lower.ptm[6], msis86r.gts3c.s,
				msis86r.gts3c.t0, msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);

			// ---- eq. a12
			zhm16 = zhm28;
			d[2] = dnet(d[2], msis86r.dmix.dm16, zhm16, xmm, 16.0);

			// ---- eq. a20b
			msis86r.gts3c.rl = log(b28*msis86r.lower.pdm[2][2] * fabs(msis86r.parm.pdl[17][2]) / b16);

			// ---- eq. a20a
			hc16 = msis86r.lower.pdm[6][2] * msis86r.parm.pdl[4][2];
			zc16 = msis86r.lower.pdm[5][2] * msis86r.parm.pdl[3][2];
			d[2] = d[2] * ccor(alt, msis86r.gts3c.rl, hc16, zc16);

			// ---- eq. a21
			hcc16 = msis86r.lower.pdm[8][2] * msis86r.parm.pdl[14][2];
			zcc16 = msis86r.lower.pdm[7][2] * msis86r.parm.pdl[13][2];
			rc16 = msis86r.lower.pdm[4][2] * msis86r.parm.pdl[15][2];
			d[2] = d[2] * ccor(alt, rc16, hcc16, zcc16);
		}

		if (mass != 48 && mass != 49)
			goto ninety;

	thirtyfive:
		// ---- **** o2 density ****
		// ---- eq. a18
		g32 = msis86r.csw.sw[21] * globe5(msis86r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis86r.parm.pddav[4]);
		msis86r.gts3c.db32 = msis86r.lower.pdm[1][4] * exp(g32)*msis86r.parm.pd[1][4];

		// ---- eq. a13 - a17
		d[4] = denss(msis86r.parmb, lsqv, alt, msis86r.gts3c.db32, tinf, msis86r.gts3c.tlb, 32.0, 0.0, t[2],
			msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0,
			msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);
		if (mass == 49)
			msis86r.gts3c.dd = msis86r.gts3c.dd + 2.0*d[4];
		else
			msis86r.gts3c.dd = d[4];

		if (alt <= altl[4] && msis86r.csw.sw[15] != 0)
		{
			// ---- eq. a19
			zh32 = msis86r.lower.pdm[3][4];
			b32 = denss(msis86r.parmb, lsqv, zh32, msis86r.gts3c.db32, tinf, msis86r.gts3c.tlb, 32.0 - xmm, -1.0,
				t[2], msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za,
				msis86r.gts3c.z0, tr12, taf);
			msis86r.dmix.dm32 = denss(msis86r.parmb, lsqv, alt, b32, tinf, msis86r.gts3c.tlb, xmm, 0.0, t[2],
				msis86r.lower.ptm[6], msis86r.gts3c.s,
				msis86r.gts3c.t0, msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);

			// ---- eq. a12
			zhm32 = zhm28;
			d[4] = dnet(d[4], msis86r.dmix.dm32, zhm32, xmm, 32.0);

			// ---- eq. a20b
			msis86r.gts3c.rl = log(b28*msis86r.lower.pdm[2][4] / b32);

			// ---- eq. a20a
			hc32 = msis86r.lower.pdm[6][4] * msis86r.parm.pdl[8][2];
			zc32 = msis86r.lower.pdm[5][4] * msis86r.parm.pdl[7][2];
			d[4] = d[4] * ccor(alt, msis86r.gts3c.rl, hc32, zc32);
		}

		if (mass != 48)
			goto ninety;

	fourty:
		// ---- **** ar density ****
		// ---- eq. a18
		g40 = msis86r.csw.sw[21] * globe5(msis86r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis86r.parm.pddav[5]);
		msis86r.gts3c.db40 = msis86r.lower.pdm[1][5] * exp(g40)*msis86r.parm.pd[1][5];
		// ---- eq. a13 - a17
		d[5] = denss(msis86r.parmb, lsqv, alt, msis86r.gts3c.db40, tinf, msis86r.gts3c.tlb, 40.0, 0.0, t[2],
			msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0,
			msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);
		msis86r.gts3c.dd = d[5];

		if (alt <= altl[5] && msis86r.csw.sw[15] != 0)
		{
			// ---- eq. a19
			zh40 = msis86r.lower.pdm[3][5];
			b40 = denss(msis86r.parmb, lsqv, zh40, msis86r.gts3c.db40, tinf, msis86r.gts3c.tlb, 40.0 - xmm, -1.0,
				t[2], msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za,
				msis86r.gts3c.z0, tr12, taf);
			msis86r.dmix.dm40 = denss(msis86r.parmb, lsqv, alt, b40, tinf, msis86r.gts3c.tlb, xmm, 0.0, t[2],
				msis86r.lower.ptm[6], msis86r.gts3c.s,
				msis86r.gts3c.t0, msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);

			// ---- eq. a12
			zhm40 = zhm28;
			d[5] = dnet(d[5], msis86r.dmix.dm40, zhm40, xmm, 40.0);

			// ---- eq. a20b
			msis86r.gts3c.rl = log(b28*msis86r.lower.pdm[2][5] / b40);

			// ---- eq. a20a
			hc40 = msis86r.lower.pdm[6][5] * msis86r.parm.pdl[10][2];
			zc40 = msis86r.lower.pdm[5][5] * msis86r.parm.pdl[9][2];
			d[5] = d[5] * ccor(alt, msis86r.gts3c.rl, hc40, zc40);
		}
		if (mass != 48)
			goto ninety;

	fourtyfive:
		// ----  **** hydrogen density ****
		// ---- eq. a18
		g1 = msis86r.csw.sw[21] * globe5(msis86r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107,
			ap, msis86r.parm.pddav[6]);
		msis86r.gts3c.db01 = msis86r.lower.pdm[1][6] * exp(g1)*msis86r.parm.pd[1][6];

		// ---- eq. a13 - a17
		d[7] = denss(msis86r.parmb, lsqv, alt, msis86r.gts3c.db01, tinf, msis86r.gts3c.tlb, 1.0, -0.4, t[2],
			msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za,
			msis86r.gts3c.z0, tr12, taf);
		msis86r.gts3c.dd = d[7];
		if (alt <= altl[7] && msis86r.csw.sw[15] != 0)
		{
			// ---- eq. a19
			zh01 = msis86r.lower.pdm[3][6];
			b01 = denss(msis86r.parmb, lsqv, zh01, msis86r.gts3c.db01, tinf, msis86r.gts3c.tlb, 1.0 - xmm, -1.4,
				t[2], msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za,
				msis86r.gts3c.z0, tr12, taf);
			msis86r.dmix.dm01 = denss(msis86r.parmb, lsqv, alt, b01, tinf, msis86r.gts3c.tlb, xmm, 0.0, t[2],
				msis86r.lower.ptm[6], msis86r.gts3c.s,
				msis86r.gts3c.t0, msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);

			// ---- eq. a12
			zhm01 = zhm28;
			d[7] = dnet(d[7], msis86r.dmix.dm01, zhm01, xmm, 1.0);

			// ---- eq. a20b
			msis86r.gts3c.rl = log(b28*msis86r.lower.pdm[2][6] * fabs(msis86r.parm.pdl[18][2]) / b01);

			// ---- eq. a20a
			hc01 = msis86r.lower.pdm[6][6] * msis86r.parm.pdl[12][2];
			zc01 = msis86r.lower.pdm[5][6] * msis86r.parm.pdl[11][2];
			d[7] = d[7] * ccor(alt, msis86r.gts3c.rl, hc01, zc01);

			// ---- eq. a21
			hcc01 = msis86r.lower.pdm[8][6] * msis86r.parm.pdl[20][2];
			zcc01 = msis86r.lower.pdm[7][6] * msis86r.parm.pdl[19][2];
			rc01 = msis86r.lower.pdm[4][6] * msis86r.parm.pdl[21][2];
			d[7] = d[7] * ccor(alt, rc01, hcc01, zcc01);
		}

	fourtyeight:
		// ----  **** atomic nitrogen density ****
		// ---- eq. a18
		g14 = msis86r.csw.sw[21] * globe5(msis86r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107,
			ap, msis86r.parm.pddav[7]);
		msis86r.gts3c.db14 = msis86r.lower.pdm[1][7] * exp(g14)*msis86r.parm.pd[1][7];

		// ---- eq. a13 - a17
		d[8] = denss(msis86r.parmb, lsqv, alt, msis86r.gts3c.db14, tinf, msis86r.gts3c.tlb, 14.0, 0.0, t[2],
			msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0,
			msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);
		msis86r.gts3c.dd = d[8];
		if (alt <= altl[8] && msis86r.csw.sw[15] != 0)
		{
			// ---- eq. a19
			zh14 = msis86r.lower.pdm[3][7];
			b14 = denss(msis86r.parmb, lsqv, zh14, msis86r.gts3c.db14, tinf, msis86r.gts3c.tlb, 14.0 - xmm, -1.0,
				t[2], msis86r.lower.ptm[6], msis86r.gts3c.s, msis86r.gts3c.t0, msis86r.gts3c.za,
				msis86r.gts3c.z0, tr12, taf);
			msis86r.dmix.dm14 = denss(msis86r.parmb, lsqv, alt, b14, tinf, msis86r.gts3c.tlb, xmm, 0.0, t[2],
				msis86r.lower.ptm[6], msis86r.gts3c.s,
				msis86r.gts3c.t0, msis86r.gts3c.za, msis86r.gts3c.z0, tr12, taf);

			// ---- eq. a12
			zhm14 = zhm28;
			d[8] = dnet(d[8], msis86r.dmix.dm14, zhm14, xmm, 14.0);

			// ---- eq. a20b
			msis86r.gts3c.rl = log(b28*msis86r.lower.pdm[2][7] * fabs(msis86r.parm.pdl[3][1]) / b14);

			// ---- eq. a20a
			hc14 = msis86r.lower.pdm[6][7] * msis86r.parm.pdl[2][1];
			zc14 = msis86r.lower.pdm[5][7] * msis86r.parm.pdl[1][1];
			d[8] = d[8] * ccor(alt, msis86r.gts3c.rl, hc14, zc14);

			// ---- eq. a21
			hcc14 = msis86r.lower.pdm[8][7] * msis86r.parm.pdl[5][1];
			zcc14 = msis86r.lower.pdm[7][7] * msis86r.parm.pdl[4][1];
			rc14 = msis86r.lower.pdm[4][7] * msis86r.parm.pdl[6][1];
			d[8] = d[8] * ccor(alt, rc14, hcc14, zcc14);
		}

		if (mass != 48)
			goto ninety;

		// ---- total mass density
		d[6] = 1.66e-24 * (4.0*d[1] + 16.0*d[2] + 28.0*d[3] +
			32.0 * d[4] + 40.0 * d[5] + d[7] + 14.0*d[8]);
		msis86r.gts3c.db48 = 1.66e-24 * (4.0*msis86r.gts3c.db04 + 16.0*msis86r.gts3c.db16 +
			28.0*msis86r.gts3c.db28 + 32.0 * msis86r.gts3c.db32 +
			40.0 * msis86r.gts3c.db40 + msis86r.gts3c.db01 +
			14.0*msis86r.gts3c.db14);
		goto ninety;

	fifty:
		// cdav ddum never gets used
		//        msis86r.gts3c.ddum = denss(msis86r.parmb, lsqv, alt,1.0, tinf,msis86r.gts3c.tlb,0.0,0.0,t[2],msis86r.lower.ptm[6],
		//                                   msis86r.gts3c.s,msis86r.gts3c.t0,
		//                                   msis86r.gts3c.za,msis86r.gts3c.z0,tr12, taf);
	ninety :
		if (msis86r.metsel.imr == 1)
		{
			for (i = 1; i <= 8; i++)
				d[i] = d[i] * 1.0e6;
			d[6] = d[6] / 1000.0;
		}
	}
	/*
	//--------------------------------------------------------------------
	// ---- calculate temperature and density profiles for msis models
	//--------------------------------------------------------------------
	*/
	// cdav do as in anline function/////////
	//double zeta
	//     (
	//       double zz, double zl, double re
	//     )
	//     {
	//     return (zz-zl)*(re+zl)/(re+zz);
	//     }

	double denss
		(
		parmbtype& parmb, lsqvtype& lsqv,
		double alt, double dlb, double tinf, double tlb, double xm, double alpha, double& tz,
		double zlb, double s2, double t0, double za, double z0, double tr12, double& taf
		)
	{
		double z, zg2, dta, denss1,
			bb, cc, densa, dd, t12, x2, zg0, zg1,
			glb, gamma,
			gamm, tt, ta, x;

		const double rgas = 831.4;

		//------------------------------ begin ---------------------------------
		denss1 = 1.0;

		if (alt >= za)
			z = alt;
		else
			z = za;

		// ---- eq. a4a
		zg2 = zeta(z, zlb, parmb.re);

		// ---- eq. a1a
		tt = tinf - (tinf - tlb)*exp(-s2*zg2);
		ta = tt;
		tz = tt;
		denss1 = tz;

		if (alt < za)
		{
			// ---- eq. a4b
			zg0 = zeta(z0, za, parmb.re);

			// ---- eq. a2b
			dta = (tinf - ta)*s2*pow((parmb.re + zlb) / (parmb.re + za), 2);

			// ---- eq. a3e
			t12 = t0 + tr12*(ta - t0);

			// ---- eq. a4b
			zg1 = zeta(alt, za, parmb.re);

			// ---- calculate temperature below za
			// ---- eq. a3a
			dd = 0.666666*zg0*dta / pow(ta, 2) - 3.11111*(1.0 / ta - 1.0 / t0) +
				7.11111*(1.0 / t12 - 1.0 / t0);

			// ---- eq. a3b
			cc = zg0*dta / (2.0*ta*ta) - (1.0 / ta - 1.0 / t0) - 2.0*dd;

			// ---- eq. a3c
			bb = (1.0 / ta - 1.0 / t0) - cc - dd;

			// ---- eq. a3d
			x = (-(zg1 - zg0) / zg0);

			// ---- eq. a1b
			x2 = x*x;
			tz = 1.0 / (1.0 / t0 + bb*x2 + cc*pow(x2, 2) + dd*pow(x2, 3));
			denss1 = tz;
			taf = (t12 - t0) / (ta - t0);
		}
		if (xm != 0.0)
		{
			if ((ta <= 0.0) | (tz <= 0.0))
			{
				//                printf("  \n",
				//                   alt,xm,tinf,tlb,t0,ta,ii,jg,n,lsqv.dv[j],ifun,s2,zg0,tz);
				tt = tlb;
				ta = tlb;
				tz = tlb;
			}

			// ---- calculate density above za
			// ---- eq. a17a
			glb = parmb.gsurf / pow(1.0 + zlb / parmb.re, 2);

			// ---- eq. a16a
			gamma = xm*glb / (s2*rgas*tinf);

			// ---- eq. a13, a14a,   a15
			densa = dlb*pow(tlb / tt, (1.0 + alpha + gamma)) * exp(-s2*gamma*zg2);
			denss1 = densa;
			if (alt < za)
			{
				// ---- calculate density below za
				// ---- eq. a17b
				glb = parmb.gsurf / pow(1.0 + za / parmb.re, 2);

				// ---- eq. a16b
				gamm = xm*glb*zg0 / rgas;

				// ---- eq. a13, a14b,   a15
				denss1 = densa * pow(ta / tz, (1.0 + alpha)) *
					exp(gamm * ((x - 1.0) / t0 + bb * (x*x2 - 1.0) / 3.0 +
					cc*(x2*x2*x - 1.0) / 5.0 + dd * (pow(x2, 3) * x - 1.0) / 7.0));
			}
		}

		return denss1;
	}

	//// cdav these next few functions could be done as inline functions if that speeds
	//// things up
	//// ---- 3hr magnetic activity functions
	////      eq. a24d
	//      double g0
	//          ( double a, double p[151])
	//          {
	//          return (a - 4.0 + (p[26]-1.0) *
	//                 (a - 4.0 + (exp(-fabs(p[25])*(a-4.0)) -1.0) / fabs(p[25]) ));
	//          }
	//// ---- eq. a24c
	//      double sumex
	//          ( double ex )
	//          {
	//          return 1.0 + (1.0-pow(ex,19)) / (1.0-ex)*sqrt(ex);
	//          }
	//// ---- eq. a24a
	//      double sg0
	//          ( double ex, double p[151], double ap[8] )
	//          {
	//          return ( g0(ap[2],p) +
	//                  ( g0(ap[3],p)*ex + g0(ap[4],p)*ex*ex + g0(ap[5],p)* pow(ex,3) +
	//                   ( g0(ap[6],p)*pow(ex,4) + g0(ap[7],p)*pow(ex,12)) * (1.0-pow(ex,8))
	//                   / (1.0-ex)
	//                  )
	//                 ) / sumex(ex);
	//          }

	/*
	//----------------------------------------------------------------------
	//       calculate g[l] function for msis-86/cira 1986
	//       upper thermosphere parameters
	//----------------------------------------------------------------------
	*/
	double globe5
		(
		cswtype& csw,
		lpolytype& lpoly,
		double yrd, double sec, double lat, double llong, double tloc, double f107a,
		double f107, double ap[8], double p[151]
		)
	{
		int i;
		double c, s, c2, c4, s2, cd14, cd18, cd32, cd39, f1, f2, t71, t72,
			t81, t82, p44, p45, exp1;
		double t[16], tinf;

		double dgtr = 1.74533e-2;
		double dr = 1.72142e-2;
		double xl = 1000.0;
		static double tll = 1000.0;
		//        static double sw9  = 1.0;
		double dayl = -1.0;
		double p14 = -1000.0;
		double p18 = -1000.0;
		double p32 = -1000.0;
		double p39 = -1000.0;
		double hr = 0.2618;
		double sr = 7.2722e-5;
		int sv[26] = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
		double nsw = 14;
		double exp2;

		//------------------------------ begin ---------------------------------
		if (csw.isw != 64999)
			tselec(csw, sv);

		t[10] = 0.0;
		t[11] = 0.0;
		t[12] = 0.0;
		t[13] = 0.0;

		lpoly.iyr = yrd / 1000.0;
		lpoly.day = yrd - lpoly.iyr*1000.0;

		// eq. a22 [remainder of code]
		if (xl != lat)
		{
			// calculate legendre polynomials
			c = sin(lat*dgtr);
			s = cos(lat*dgtr);
			c2 = c*c;
			c4 = c2*c2;
			s2 = s*s;
			lpoly.plg[2][1] = c;
			lpoly.plg[3][1] = 0.5* (3.0*c2 - 1.0);
			lpoly.plg[4][1] = 0.5* (5.0*c*c2 - 3.0*c);
			lpoly.plg[5][1] = (35.0*c4 - 30.0*c2 + 3.0) / 8.0;
			lpoly.plg[6][1] = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c) / 8.0;
			lpoly.plg[7][1] = (11.0*c*lpoly.plg[6][1] - 5.0*lpoly.plg[5][1]) / 6.0;
			lpoly.plg[2][2] = s;
			lpoly.plg[3][2] = 3.0*c*s;
			lpoly.plg[4][2] = 1.5* (5.0*c2 - 1.0)*s;
			lpoly.plg[5][2] = 2.5* (7.0*c2*c - 3.0*c)*s;
			lpoly.plg[6][2] = 1.875* (21.0*c4 - 14.0*c2 + 1.0)*s;
			lpoly.plg[7][2] = (11.0*c*lpoly.plg[6][2] - 6.0*lpoly.plg[5][2]) / 5.0;
			lpoly.plg[3][3] = 3.0*s2;
			lpoly.plg[4][3] = 15.0*s2*c;
			lpoly.plg[5][3] = 7.5*(7.0*c2 - 1.0)*s2;
			lpoly.plg[6][3] = 3.0*c*lpoly.plg[5][3] - 2.0*lpoly.plg[4][3];
			lpoly.plg[7][3] = (11.0*c*lpoly.plg[6][3] - 7.0*lpoly.plg[5][3]) / 4.0;
			lpoly.plg[8][3] = (13.0*c*lpoly.plg[7][3] - 8.0*lpoly.plg[6][3]) / 5.0;
			lpoly.plg[4][4] = 15.0*s2*s;
			lpoly.plg[5][4] = 105.0*s2*s*c;
			lpoly.plg[6][4] = (9.0*c*lpoly.plg[5][4] - 7.0*lpoly.plg[4][4]) / 2.0;
			lpoly.plg[7][4] = (11.0*c*lpoly.plg[6][4] - 8.0*lpoly.plg[5][4]) / 3.0;
			xl = lat;
		}

		if (tll != tloc)
		{
			lpoly.stloc = sin(hr*tloc);
			lpoly.ctloc = cos(hr*tloc);
			lpoly.s2tloc = sin(2.0*hr*tloc);
			lpoly.c2tloc = cos(2.0*hr*tloc);
			lpoly.s3tloc = sin(3.0*hr*tloc);
			lpoly.c3tloc = cos(3.0*hr*tloc);
			tll = tloc;
		}
		if ((lpoly.day != dayl) | (p[14] != p14))
			cd14 = cos(dr*(lpoly.day - p[14]));
		// cdav c2d14 is never used
		//        if (lpoly.day != dayl | p[14] != p14)
		//            c2d14 = cos(dr*2*(lpoly.day-p[14]));
		if ((lpoly.day != dayl) | (p[18] != p18))
			cd18 = cos(2.0*dr*(lpoly.day - p[18]));
		if ((lpoly.day != dayl) | (p[32] != p32))
			cd32 = cos(dr*(lpoly.day - p[32]));
		if ((lpoly.day != dayl) | (p[39] != p39))
			cd39 = cos(2.0*dr*(lpoly.day - p[39]));
		dayl = lpoly.day;
		p14 = p[14];
		p18 = p[18];
		p32 = p[32];
		p39 = p[39];

		// ---- f10.7 effect ----
		lpoly.df = f107 - f107a;
		lpoly.dfa = f107a - 150.0;
		t[1] = p[20] * lpoly.df + p[21] * lpoly.df*lpoly.df + p[22] * lpoly.dfa + p[30] * pow(lpoly.dfa, 2);
		f1 = 1.0 + (p[48] * lpoly.dfa + p[20] * lpoly.df + p[21] * lpoly.df*lpoly.df)*csw.swc[1];
		f2 = 1.0 + (p[50] * lpoly.dfa + p[20] * lpoly.df + p[21] * lpoly.df*lpoly.df)*csw.swc[1];

		// ---- time independent ----
		t[2] = (p[2] * lpoly.plg[3][1] + p[3] * lpoly.plg[5][1] + p[23] * lpoly.plg[7][1])
			+ (p[15] * lpoly.plg[3][1])*lpoly.dfa*csw.swc[1] + p[27] * lpoly.plg[2][1];

		// ---- symmetrical annual ----
		t[3] = p[19] * cd32;

		// ---- symmetrical semiannual ----
		t[4] = (p[16] + p[17] * lpoly.plg[3][1])*cd18;

		// ---- asymmetrical annual ----
		t[5] = f1 * (p[10] * lpoly.plg[2][1] + p[11] * lpoly.plg[4][1])*cd14;

		// ---- asymmetrical semiannual ----
		t[6] = p[38] * lpoly.plg[2][1] * cd39;

		// ---- diurnal ----
		t71 = (p[12] * lpoly.plg[3][2] + p[36] * lpoly.plg[2][2])*cd14*csw.swc[5];
		t72 = (p[13] * lpoly.plg[3][2] + p[37] * lpoly.plg[2][2])*cd14*csw.swc[5];
		t[7] = f2*
			((p[4] * lpoly.plg[2][2] + p[5] * lpoly.plg[4][2] + p[28] * lpoly.plg[6][2]
			+ t71)*lpoly.ctloc
			+ (p[7] * lpoly.plg[2][2] + p[8] * lpoly.plg[4][2] + p[29] * lpoly.plg[6][2]
			+ t72)*lpoly.stloc);

		// ----  semidiurnal ----
		t81 = (p[24] * lpoly.plg[4][3])*cd14*csw.swc[5];
		t82 = (p[34] * lpoly.plg[4][3])*cd14*csw.swc[5];
		t[8] = f2*
			((p[6] * lpoly.plg[3][3] + p[42] * lpoly.plg[5][3] + t81)*lpoly.c2tloc
			+ (p[9] * lpoly.plg[3][3] + p[43] * lpoly.plg[5][3] + t82)*lpoly.s2tloc);

		// ---- terdiurnal ----
		t[14] = f2*
			((p[40] * lpoly.plg[4][4] + (p[94] * lpoly.plg[5][4] + p[47] * lpoly.plg[7][4])*cd14*
			csw.swc[5])*lpoly.s3tloc
			+ (p[41] * lpoly.plg[4][4] + (p[95] * lpoly.plg[5][4] + p[49] * lpoly.plg[7][4])*cd14*
			csw.swc[5])*lpoly.c3tloc);

		// ----    magnetic activity based on daily ap
		if ((csw.sw[9] != -1.0) | (p[52] == 0.0))
		{
			lpoly.apd = ap[1] - 4.0;
			p44 = p[44];
			p45 = p[45];
			//cdav
			if (p44 <= 0.0)
				p44 = 1.0e-5;
			lpoly.apdf = (lpoly.apd + (p45 - 1.0)*(lpoly.apd + (exp(-p44 *lpoly.apd) - 1.0) / p44));
			t[9] = lpoly.apdf*(p[33] + p[46] * lpoly.plg[3][1] + p[35] * lpoly.plg[5][1] +
				(p[101] * lpoly.plg[2][1] + p[102] * lpoly.plg[4][1] + p[103] * lpoly.plg[6][1])*cd14*
				csw.swc[5] + (p[122] * lpoly.plg[2][2] + p[123] * lpoly.plg[4][2] + p[124] * lpoly.plg[6][2])
				*csw.swc[7] * cos(hr*(tloc - p[125])));
		}
		else
		{
			exp1 = exp(-10800.0*fabs(p[52])) / (1.0 + p[139] * (45.0 - fabs(lat)));
			if (exp1 > 0.99999)
				exp1 = 0.99999;
			exp2 = exp(-10800.0*fabs(p[54]));
			if (exp2 > 0.99999)
				exp2 = 0.99999;

			if (p[25] < 1.0e-4)
				p[25] = 1.0e-4;

			lpoly.apt[1] = sg0(exp1, p, ap);
			lpoly.apt[3] = sg0(exp2, p, ap);
			t[9] = lpoly.apt[1] * (p[51] + p[97] * lpoly.plg[3][1] + p[55] * lpoly.plg[5][1] +
				(p[126] * lpoly.plg[2][1] + p[127] * lpoly.plg[4][1] + p[128] * lpoly.plg[6][1])*
				cd14*csw.swc[5] +
				(p[129] * lpoly.plg[2][2] + p[130] * lpoly.plg[4][2] + p[131] * lpoly.plg[6][2])*
				csw.swc[7] * cos(hr*(tloc - p[132])));
		}

		if ((csw.sw[10] == 0) | (llong <= -1000.0)) goto fourtynine;
		// ----  longitudinal
		t[11] = (1.0 + p[90] * lpoly.plg[2][1])*(1.0 + p[81] * lpoly.dfa*csw.swc[1])*
			((p[65] * lpoly.plg[3][2] + p[66] * lpoly.plg[5][2] + p[67] * lpoly.plg[7][2]
			+ p[104] * lpoly.plg[2][2] + p[105] * lpoly.plg[4][2] + p[106] * lpoly.plg[6][2]
			+ csw.swc[5] * (p[110] * lpoly.plg[2][2] + p[111] * lpoly.plg[4][2] + p[112] *
			lpoly.plg[6][2])*cd14) * cos(dgtr*llong)
			+ (p[91] * lpoly.plg[3][2] + p[92] * lpoly.plg[5][2] + p[93] * lpoly.plg[7][2]
			+ p[107] * lpoly.plg[2][2] + p[108] * lpoly.plg[4][2] + p[109] * lpoly.plg[6][2]
			+ csw.swc[5] * (p[113] * lpoly.plg[2][2] + p[114] * lpoly.plg[4][2] + p[115] *
			lpoly.plg[6][2])*cd14) * sin(dgtr*llong));
		// ----  ut and mixed ut,longitude
		t[12] = (1.0 + p[96] * lpoly.plg[2][1])*(1.0 + p[82] * lpoly.dfa*csw.swc[1])*
			(1.0 + p[120] * lpoly.plg[2][1] * csw.swc[5] * cd14)*
			((p[69] * lpoly.plg[2][1] + p[70] * lpoly.plg[4][1] + p[71] * lpoly.plg[6][1])*
			cos(sr*(sec - p[72])));
		t[12] = t[12] + csw.swc[11] *
			(p[77] * lpoly.plg[4][3] + p[78] * lpoly.plg[6][3] + p[79] * lpoly.plg[8][3])*
			cos(sr*(sec - p[80]) + 2.0*dgtr*llong)*(1.0 + p[138] *
			lpoly.dfa*csw.swc[1]);

		// ----  ut,longitude magnetic activity
		if ((csw.sw[9] != -1.0) | (p[52] == 0.0))
		{
			t[13] = lpoly.apdf*csw.swc[11] * (1.0 + p[121] * lpoly.plg[2][1])*
				((p[61] * lpoly.plg[3][2] + p[62] * lpoly.plg[5][2] + p[63] * lpoly.plg[7][2])*
				cos(dgtr*(llong - p[64])))
				+ lpoly.apdf*csw.swc[11] * csw.swc[5] *
				(p[116] * lpoly.plg[2][2] + p[117] * lpoly.plg[4][2] + p[118] * lpoly.plg[6][2])*
				cd14*cos(dgtr*(llong - p[119]))
				+ lpoly.apdf*csw.swc[12] *
				(p[84] * lpoly.plg[2][1] + p[85] * lpoly.plg[4][1] + p[86] * lpoly.plg[6][1])*
				cos(sr*(sec - p[76]));
		}
		else
		{
			t[13] = lpoly.apt[1] * csw.swc[11] * (1.0 + p[133] * lpoly.plg[2][1])*
				((p[53] * lpoly.plg[3][2] + p[99] * lpoly.plg[5][2] + p[68] * lpoly.plg[7][2])*
				cos(dgtr*(llong - p[98])))
				+ lpoly.apt[1] * csw.swc[11] * csw.swc[5] *
				(p[134] * lpoly.plg[2][2] + p[135] * lpoly.plg[4][2] + p[136] * lpoly.plg[6][2])*
				cd14*cos(dgtr*(llong - p[137]))
				+ lpoly.apt[1] * csw.swc[12] *
				(p[56] * lpoly.plg[2][1] + p[57] * lpoly.plg[4][1] + p[58] * lpoly.plg[6][1])*
				cos(sr*(sec - p[59]));
		}

		// ----  parms not used: 60,83,100,140-150
	fourtynine:
		tinf = 0.0;
		if (csw.sw[9] == -1)
			tinf = p[31];
		for (i = 1; i <= nsw; i++)
			tinf = tinf + fabsf(csw.sw[i]) * t[i];
		return tinf;
	}

	/*
	//----------------------------------------------------------------------
	//       limited parameter version of globe 9/2/82
	//       calculate g[l] function for msis-86/cira 1986
	//       lower thermosphere parameters
	//----------------------------------------------------------------------
	*/
	double glob5l
		(
		lpolytype& lpoly,
		cswtype& csw,
		double p[151]
		)
	{
		double cd7, cd9, cd11, t[15], tt;
		int i;
		double dr = 1.72142e-2;
		double dayl = -1.0;
		double p7 = -1000.0;
		double p9 = -1000.0;
		double p11 = -1000.0;

		//------------------------------ begin ---------------------------------
		for (i = 1; i <= 14; i++)
			t[i] = 0.0;

		if ((lpoly.day != dayl) | (p7 != p[7]))
			cd7 = cos(dr*(lpoly.day - p[7]));

		if ((lpoly.day != dayl) | (p9 != p[9]))
			cd9 = cos(2.0*dr*(lpoly.day - p[9]));

		if ((lpoly.day != dayl) | (p11 != p[11]))
			cd11 = cos(dr*(lpoly.day - p[11]));

		dayl = lpoly.day;
		p7 = p[7];
		p9 = p[9];
		p11 = p[11];

		t[1] = p[2] * lpoly.dfa;
		t[2] = p[4] * lpoly.plg[3][1];
		t[3] = p[6] * cd7;
		t[4] = p[8] * cd9;
		t[5] = (p[10] * lpoly.plg[2][1] + p[22] * lpoly.plg[4][1])*cd11;
		t[6] = 0.0;
		t[7] = p[14] * lpoly.plg[2][2] * lpoly.ctloc + p[15] * lpoly.plg[2][2] * lpoly.stloc;
		t[8] = (p[16] * lpoly.plg[3][3] + p[18] * lpoly.plg[5][3]
			+ (p[20] * lpoly.plg[6][3])*cd11*csw.swc[5])*lpoly.c2tloc
			+ (p[17] * lpoly.plg[3][3] + p[19] * lpoly.plg[5][3]
			+ (p[21] * lpoly.plg[6][3])*cd11*csw.swc[5]
			)*lpoly.s2tloc;
		t[14] = p[12] * lpoly.plg[4][4] * lpoly.c3tloc + p[25] * lpoly.plg[4][4] * lpoly.s3tloc;
		if (csw.sw[9] == 1)
			t[9] = lpoly.apdf*(p[23] + p[24] * lpoly.plg[3][1] * csw.swc[2]);

		if (csw.sw[9] == -1)
			t[9] = (p[3] * lpoly.apt[3] + p[5] * lpoly.plg[3][1] * lpoly.apt[3] * csw.swc[2]);

		// ---- parms not used: 13
		tt = 0.0;
		for (i = 1; i <= 14; i++)
			tt = tt + fabsf(csw.sw[i]) * t[i];

		return tt;
	}

	/*
	//----------------------------------------------------------------------
	// calculates geomagnetic longitude [xlm] and latitude [bm]
	// from geografic longitude [xlg] and latitude [bg] for art = 0
	// and reverse for art = 1.0 all angles in degree.
	// latitude:-90 to 90. longitude:0 to 360 east.
	//----------------------------------------------------------------------
	*/
	void ggm
		(
		int art, double& xlg, double& bg, double& xlm, double& bm
		)
	{
		double cbg, cbm, ci, clg, clm, sbg, sbm, si, slg, slm, ylg, zpi;
		double faktor = 0.0174532925;

		//------------------------------ begin ---------------------------------
		zpi = faktor*360.0;
		cbg = 11.4*faktor;
		ci = cos(cbg);
		si = sin(cbg);

		if (art != 0)
		{
			cbm = cos(bm*faktor);
			sbm = sin(bm*faktor);
			clm = cos(xlm*faktor);
			slm = sin(xlm*faktor);
			sbg = sbm*ci - cbm*clm*si;
			bg = asin(sbg);
			cbg = cos(bg);
			slg = (cbm*slm) / cbg;
			clg = (sbm*si + cbm*clm*ci) / cbg;
			if (fabs(clg) > 1.0)
			{
				if (clg > 0.0)
					clg = 1.0;
				else
					clg = -1.0;
			}
			xlg = acos(clg);
			if (slg < 0.0)
				xlg = zpi - acos(clg);

			bg = bg / faktor;
			xlg = xlg / faktor;
			xlg = xlg - 69.8;
			if (xlg < 0.0)
				xlg = xlg + 360.0;
		}
		else
		{
			ylg = xlg + 69.8;
			cbg = cos(bg*faktor);
			sbg = sin(bg*faktor);
			clg = cos(ylg*faktor);
			slg = sin(ylg*faktor);
			sbm = sbg*ci + cbg*clg*si;
			bm = asin(sbm);
			cbm = cos(bm);
			slm = (cbg*slg) / cbm;
			clm = (-sbg*si + cbg*clg*ci) / cbm;
			xlm = acos(clm);
			if (slm < 0.0)
				xlm = zpi - acos(clm);

			bm = bm / faktor;
			xlm = xlm / faktor;
		}
	}

	/*
	//----------------------------------------------------------------------
	//      cira     11-feb-86
	//----------------------------------------------------------------------
	*/
	void msis86init
		(
		msistype& msis86r
		)
	{
		int i, j;

		// cdav i'm not string literate yet in c++ :-)
		//        msis86r.datime.isdate = "11-feb-86";
		//        msis86r.datime.istime = "18:23:31 ";
		//        msis86r.datime.name   = "msis-86  ";

		msis86r.parmb.gsurf = 980.665;
		msis86r.parmb.re = 6356.77;
		//        iiee  = 0

		msis86r.csw.isw = 0;
		//
		// cdav   set this to output in
		msis86r.metsel.imr = 1;

		// ----   temperature
		static double pt[151] =
		{ 0,
		0.996040e+00, 0.385528e-01, 0.303445e-02, -0.105531e+00, -0.607134e-02,
		-0.516278e-03, -0.115622e+00, 0.202240e-02, 0.990156e-02, -0.127371e+00,
		-0.302449e-01, 0.123512e-01, -0.526277e-02, -0.845398e+01, 0.000000e+00,
		0.142370e-01, 0.000000e+00, 0.125818e+03, 0.805486e-02, 0.164419e-02,
		-0.621452e-05, 0.311701e-02, 0.000000e+00, 0.386578e-02, 0.132397e+00,
		0.213315e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -0.641110e-05,
		0.000000e+00, 0.300150e+02, 0.533297e-02, 0.389146e-02, 0.204725e-02,
		0.000000e+00, 0.000000e+00, -0.192645e-01, 0.275905e+01, 0.147284e-02,
		0.341345e-03, -0.117388e-02, -0.354589e-03, 0.113139e+00, 0.169134e+00,
		0.508295e-02, 0.365016e-04, 0.426385e-02, 0.115102e-03, 0.511819e-02,
		0.609108e-02, 0.404995e-04, 0.153049e-02, 0.241470e-04, 0.230764e-02,
		0.155267e-02, 0.133722e-02, -0.182318e-02, -0.263007e+03, 0.000000e+00,
		0.137337e-02, 0.995774e-03, 0.000000e+00, -0.108983e+03, 0.562606e-02,
		0.594053e-02, 0.109358e-02, 0.000000e+00, -0.133410e-01, -0.243409e-01,
		-0.135688e-01, 0.311370e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		-0.283023e+04, 0.845583e-03, 0.538706e-03, 0.000000e+00, 0.247956e+03,
		0.292246e-02, 0.000000e+00, 0.000000e+00, 0.747703e-04, 0.887993e-03,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		-0.116540e-01, -0.449173e-02, -0.353189e-03, -0.173933e-03, -0.153218e-03,
		-0.565411e+00, 0.777272e-02, -0.911784e+02, 0.645187e-03, 0.000000e+00,
		-0.837685e-03, 0.242318e-02, 0.473796e-02, -0.301801e-02, -0.423564e-02,
		-0.248289e-02, 0.919286e-03, 0.216372e-02, 0.863968e-03, 0.189689e-02,
		0.415654e-02, 0.000000e+00, 0.118068e-01, 0.331190e-02, 0.000000e+00,
		0.120222e-02, 0.000000e+00, 0.000000e+00, -0.307246e+01, 0.000000e+00,
		0.000000e+00, 0.672403e-03, 0.108930e-02, 0.972278e-03, 0.468242e+01,
		-0.315034e-03, 0.400059e-02, 0.515036e-02, 0.162989e-02, 0.108824e-02,
		0.995261e-03, 0.418955e+01, -0.364059e+00, 0.170182e-02, 0.000000e+00,
		0.000000e+00, -0.320120e+01, 0.000000e+00, 0.580206e-02, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		};

		static double pd[8][151] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			},
			// ----   he density
			{ 0,
			0.104934e+01, -0.288362e-01, -0.207095e+00, -0.103314e+00, -0.702373e-02,
			0.129664e-01, 0.408853e+00, -0.919895e-02, -0.188660e-01, 0.140927e+01,
			0.175033e+00, 0.187351e-01, 0.110979e+00, -0.742871e+01, 0.000000e+00,
			0.267143e+00, -0.595979e-01, 0.105038e+03, -0.840963e-01, -0.697632e-03,
			0.206521e-05, 0.765306e-03, 0.000000e+00, 0.000000e+00, 0.126762e+00,
			0.128876e+00, -0.504479e-01, -0.130735e-01, -0.224348e-01, 0.000000e+00,
			0.000000e+00, -0.150832e+03, -0.629928e-02, 0.000000e+00, -0.407760e-02,
			0.000000e+00, 0.000000e+00, 0.525725e-01, -0.311486e+02, -0.313351e-02,
			0.275838e-02, 0.000000e+00, 0.000000e+00, 0.111247e+00, 0.108815e+00,
			-0.466713e-01, 0.000000e+00, -0.329329e-02, 0.000000e+00, 0.167838e-02,
			-0.916691e-02, 0.345044e-04, -0.971806e-02, 0.000000e+00, -0.204672e-02,
			-0.786899e-02, -0.798285e-02, 0.536515e-02, -0.531172e+04, 0.000000e+00,
			-0.642781e-02, -0.171690e-02, 0.000000e+00, -0.679131e+02, -0.179912e-01,
			-0.158305e-01, -0.712313e-02, 0.000000e+00, 0.253477e-01, 0.852960e-01,
			0.102163e+00, 0.295009e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			-0.684625e+04, -0.619098e-02, -0.269289e-02, 0.000000e+00, -0.520231e+03,
			-0.633463e-02, 0.000000e+00, 0.000000e+00, -0.602428e-02, -0.407077e-02,
			0.542264e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.407560e-01, 0.282288e-01, 0.908088e-02, 0.000000e+00, 0.000000e+00,
			-0.405204e+00, -0.597931e-01, -0.731823e+02, -0.206620e-02, 0.000000e+00,
			-0.372723e-02, -0.188146e-01, -0.101794e-01, 0.804633e-02, 0.101090e-01,
			0.873253e-02, 0.238268e-01, 0.480444e-02, 0.171088e-02, 0.396369e-01,
			-0.213809e-01, 0.000000e+00, -0.102588e+00, -0.591702e-02, 0.000000e+00,
			0.270923e-02, 0.000000e+00, 0.000000e+00, -0.175043e+03, 0.603489e+00,
			-0.617589e+00, 0.838098e-02, 0.183871e-02, -0.705329e-03, -0.406644e+01,
			-0.509347e-02, -0.284344e-01, -0.124160e-01, 0.133665e-01, 0.393410e-02,
			-0.503723e-03, -0.457683e+01, -0.529542e+00, -0.425812e-02, 0.000000e+00,
			0.000000e+00, 0.191541e+02, 0.000000e+00, 0.323247e-02, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
			},
			// ----   o density
			{ 0,
			0.931113e+00, -0.138721e+00, -0.133457e+00, -0.529542e-01, -0.444983e-02,
			0.135264e-01, 0.598075e-01, -0.362880e-01, -0.312798e-01, 0.372068e+00,
			0.295974e-01, 0.120509e-01, 0.521995e-01, -0.778888e+01, 0.000000e+00,
			0.118634e+00, -0.204495e-01, 0.103280e+03, 0.982432e-01, 0.477694e-03,
			0.000000e+00, 0.274372e-02, 0.000000e+00, 0.000000e+00, 0.757809e-01,
			0.171403e+00, -0.105205e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, -0.873348e+01, -0.581094e-02, 0.000000e+00, -0.814944e-02,
			0.000000e+00, 0.000000e+00, 0.517255e-01, -0.153028e+02, -0.348932e-02,
			0.961771e-03, 0.557732e-02, -0.454180e-03, 0.988213e-01, 0.940456e-01,
			-0.318797e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.232122e-02,
			-0.600220e-02, 0.277654e-04, -0.322019e-02, 0.000000e+00, -0.378551e-02,
			-0.334809e-02, -0.170668e-02, 0.000000e+00, 0.636184e+04, 0.000000e+00,
			0.159986e-02, -0.388204e-02, -0.164825e-02, -0.747955e+02, -0.105360e-01,
			-0.945723e-02, -0.159824e-02, -0.706730e-03, -0.168513e-01, -0.113023e+00,
			-0.636637e-01, -0.137709e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			-0.152368e+05, -0.586061e-02, -0.253108e-02, 0.000000e+00, -0.254837e+04,
			-0.328988e-02, 0.000000e+00, 0.000000e+00, -0.276364e-02, 0.967923e-02,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.434255e-01, 0.114020e-01, -0.618447e-02, 0.000000e+00, 0.000000e+00,
			-0.302568e+00, -0.327694e-01, -0.671589e+02, -0.228340e-02, 0.000000e+00,
			0.306230e-02, -0.465113e-02, -0.973421e-02, 0.128326e-01, 0.788553e-02,
			0.797197e-02, -0.120760e-01, -0.767547e-02, -0.120755e-02, -0.298523e-01,
			-0.126560e-01, 0.000000e+00, -0.568350e-01, -0.153039e-01, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.242911e-02, -0.401347e-02, -0.219074e-02, 0.311281e+01,
			0.323251e-02, -0.639523e-02, -0.663069e-02, -0.304403e-03, -0.401920e-02,
			-0.118708e-02, 0.415211e+01, -0.201896e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
			},
			// ----   n2 density
			{ 0,
			0.106903e+01, 0.377113e-03, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.898481e-01, -0.236325e+02, 0.208180e-01, 0.139638e+03, -0.119444e+00,
			-0.845398e+01, -0.399776e-05, 0.000000e+00, 0.366210e-02, -0.178929e-02,
			0.190412e-01, -0.392257e-01, 0.632343e-02, 0.548144e-02, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -0.243022e-02,
			0.976619e+00, 0.568478e-03, 0.582026e-02, 0.000000e+00, 0.621998e-02,
			0.000000e+00, 0.000000e+00, 0.107674e-01, 0.893820e+02, -0.192414e-01,
			-0.845398e+01, 0.000000e+00, 0.000000e+00, -0.200200e-01, -0.195833e-02,
			-0.938391e-02, 0.131480e-01, -0.260147e-02, -0.808556e-03, 0.511651e-04,
			0.255717e-02, 0.000000e+00, 0.466814e-02, 0.664196e-02, 0.000000e+00,
			0.998594e+00, 0.190038e-03, 0.000000e+00, -0.243825e-01, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.522105e-01,
			-0.845398e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.767271e-02, 0.564539e-02, -0.270623e-02, -0.526454e-03, 0.137075e-02,
			0.133060e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.949197e+00, 0.000000e+00, 0.000000e+00, -0.768008e-01, 0.000000e+00,
			0.000000e+00, 0.000000e+00, -0.137993e-01, -0.140136e+01, 0.120481e+00,
			-0.845398e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.987746e-02, 0.175330e-02, -0.688835e-03, 0.287022e-02, 0.000000e+00,
			0.000000e+00, 0.744513e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.152840e+00, 0.000000e+00, 0.000000e+00, 0.116252e+01, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -0.649190e+00,
			-0.845398e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			-0.584949e-01, -0.102105e+00, 0.299153e-01, -0.486227e-01, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
			},
			// ----   msis00r.gts3c.tlb
			{ 0,
			0.931402e+00, 0.137976e+00, 0.000000e+00, 0.323736e-03, 0.000000e+00,
			-0.910906e-02, 0.707506e-01, 0.000000e+00, -0.516650e-01, 0.689755e-01,
			0.000000e+00, 0.000000e+00, 0.000000e+00, -0.845398e+01, 0.000000e+00,
			0.281140e-01, 0.000000e+00, 0.736009e+02, 0.596604e-01, 0.000000e+00,
			0.000000e+00, -0.151792e-02, 0.000000e+00, 0.000000e+00, 0.132397e+00,
			0.213315e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.948758e+01, 0.884541e-02, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.113139e+00, 0.169134e+00,
			0.145192e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.107906e-01, 0.299942e-04, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -0.148930e-01,
			-0.787184e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			-0.683420e-01, -0.441778e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.229730e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
			},
			// ----   o2 density
			{ 0,
			0.868053e+00, 0.236364e+00, 0.134306e+00, 0.103086e-01, 0.000000e+00,
			-0.379164e-02, -0.157806e+00, 0.000000e+00, -0.587644e-01, -0.312508e+00,
			0.000000e+00, 0.437387e-01, -0.354091e-01, -0.223636e+02, 0.000000e+00,
			-0.533976e-01, 0.000000e+00, 0.114091e+03, 0.517497e-01, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.132397e+00,
			0.213315e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.342702e+03, 0.157033e-01, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -0.366278e-02,
			-0.116193e-02, 0.000000e+00, 0.000000e+00, 0.113139e+00, 0.169134e+00,
			0.178431e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.162864e-01, 0.316963e-04, 0.127968e-01, 0.000000e+00, 0.000000e+00,
			-0.704599e-02, 0.207921e-02, 0.636660e-02, 0.229940e+05, 0.000000e+00,
			0.127833e-01, -0.208036e-02, -0.461820e-02, -0.629391e+02, -0.120745e-01,
			0.136675e-01, 0.136011e-01, -0.537162e-02, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.192509e+05, 0.835522e-02, 0.419439e-02, 0.000000e+00, 0.120366e+05,
			0.000000e+00, 0.000000e+00, 0.000000e+00, -0.100034e-01, -0.233267e-02,
			0.972374e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			-0.265079e-01, -0.209125e-01, -0.109465e-01, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.217252e-01, -0.712385e+02, -0.189428e-02, 0.000000e+00,
			-0.602006e-02, 0.169058e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.290646e-01,
			0.348971e-02, 0.000000e+00, 0.501174e-01, 0.550595e-01, 0.000000e+00,
			-0.955897e-02, 0.000000e+00, 0.000000e+00, -0.151693e+04, 0.000000e+00,
			0.000000e+00, 0.129306e-01, 0.269567e-02, 0.000000e+00, 0.392243e+01,
			-0.847690e-02, 0.116896e-01, 0.000000e+00, 0.148967e-01, 0.544521e-02,
			0.000000e+00, 0.564918e+01, 0.000000e+00, -0.772178e-02, 0.000000e+00,
			0.000000e+00, -0.734042e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
			},
			// ----   ar density
			{ 0,
			0.127515e+01, -0.210472e+00, -0.177924e+00, 0.218900e+00, 0.288436e-01,
			0.190077e-01, 0.291001e+00, 0.217437e-01, -0.105186e-01, 0.436141e+00,
			0.107605e+00, 0.330755e-01, 0.400581e-01, -0.958051e+01, 0.000000e+00,
			0.154028e-01, 0.000000e+00, 0.734194e+02, 0.496540e-01, -0.595906e-02,
			0.384512e-04, -0.136000e-01, 0.000000e+00, 0.000000e+00, 0.132397e+00,
			0.213315e+00, -0.416610e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.146276e+03, -0.198408e-01, 0.000000e+00, 0.132530e-01,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, -0.104687e-03,
			-0.147562e-02, 0.000000e+00, 0.000000e+00, 0.113139e+00, 0.169134e+00,
			-0.126913e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, -0.608370e-02,
			-0.257587e-01, 0.319022e-04, 0.000000e+00, 0.000000e+00, 0.156644e-01,
			0.103640e-01, 0.105771e-02, 0.000000e+00, 0.357949e+04, 0.000000e+00,
			-0.125672e-02, 0.152783e-02, 0.130518e-02, 0.755558e+01, -0.920341e-02,
			-0.209142e-01, -0.134106e-01, 0.000000e+00, -0.483312e-01, 0.830900e-01,
			0.988009e-01, -0.141148e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			-0.105513e+04, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.673442e-02, 0.201691e-02,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.598019e-01, 0.633298e-02, -0.112871e-02, 0.000000e+00, 0.000000e+00,
			0.000000e+00, -0.128604e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			-0.494960e-02, -0.136415e-01, -0.115039e-01, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, -0.586860e-02, -0.141732e-02, 0.213697e-02, 0.263845e+01,
			-0.834186e-02, -0.187336e-01, -0.190870e-01, -0.803810e-02, -0.284279e-02,
			0.256722e-02, 0.171429e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
			},
			{ 0,
			0.573587e+02, -0.398747e+00, 0.000000e+00, -0.529554e+00, -0.582186e-02,
			0.714177e-01, -0.679279e+00, -0.167715e+00, -0.642434e-01, -0.211569e+00,
			-0.159922e+00, -0.171024e-03, -0.115885e+00, 0.651603e+01, 0.000000e+00,
			-0.176683e+00, 0.650395e-01, 0.143504e+01, 0.928208e-01, 0.511662e-02,
			0.000000e+00, 0.995121e-02, 0.000000e+00, 0.000000e+00, 0.132397e+00,
			0.213315e+00, 0.101451e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.567667e+02, 0.238192e-02, 0.000000e+00, -0.188240e-01,
			0.000000e+00, 0.000000e+00, 0.476218e-01, 0.235206e+02, 0.475901e-02,
			0.576162e-02, 0.151815e-01, -0.192730e-01, 0.113139e+00, 0.169134e+00,
			-0.288771e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.118418e-02,
			-0.368927e-02, 0.314704e-04, 0.882198e-02, 0.000000e+00, -0.192562e-01,
			-0.258674e-02, -0.219913e-01, 0.000000e+00, 0.438655e+04, 0.000000e+00,
			0.760126e-02, 0.259438e-02, 0.172310e-02, 0.779204e+02, 0.797786e-03,
			-0.770510e-02, 0.190982e-02, 0.272707e-02, 0.101016e-01, 0.116537e+00,
			-0.312236e-02, 0.139783e+05, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			-0.130712e+04, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, -0.320544e-02, -0.206970e-01,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.159010e-01, -0.191427e-02, -0.342829e-01, 0.000000e+00, 0.000000e+00,
			0.000000e+00, -0.345379e-01, 0.894518e+02, 0.171556e-02, 0.000000e+00,
			-0.765278e-02, -0.208987e-03, -0.157393e-01, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, -0.860673e-02, -0.119922e-01, -0.646356e-02, -0.300107e+01,
			-0.932511e-02, -0.150205e-01, -0.867835e-02, -0.764801e-02, -0.131495e-01,
			-0.676720e-02, -0.182396e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
			}
		};
		// ----    msis86r.gts3c.ps param
		static double ps[151] =
		{ 0,
		0.951363e+00, -0.467542e-01, 0.120260e+00, 0.000000e+00, 0.000000e+00,
		0.191357e-01, 0.000000e+00, 0.000000e+00, 0.125429e-02, -0.133240e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, -0.845398e+01, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.252317e-02, 0.000000e+00, -0.973404e-02, 0.132397e+00,
		0.213315e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, -0.718482e-03, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.787683e-02, -0.233698e-02, 0.113139e+00, 0.169134e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
		0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
		};
		// ----    turbo
		static double pdl[3][26] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0
			},
			{ 0,
			0.933804e+00, 0.547446e+01, 0.153263e+00, 0.919303e+00, 0.164109e+02,
			0.427083e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
			0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00
			},
			{ 0,
			0.115897e+01, 0.471094e+00, 0.109459e+01, 0.525012e+01, 0.100000e+01,
			0.100000e+01, 0.103999e+01, 0.767132e+00, 0.110514e+01, 0.175636e+01,
			0.110845e+01, 0.233439e+01, 0.796532e+00, 0.431520e+01, 0.407300e+01,
			0.101885e+01, 0.239547e+00, 0.253791e-05, 0.842931e+00, 0.104192e+01,
			0.200202e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01
			}
		};
		// ----   lower boundary
		static double ptm[9] =
		{ 0,
		0.104130e+04, 0.386000e+03, 0.190000e+03, 0.166728e+02, 0.115000e+03,
		0.120000e+03, 0.945537e+02, 0.000000e+00
		};
		static double pdm[8][9] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0
			},
			{ 0,
			0.245600e+08, 0.671072e-05, 0.100000e+03, 0.000000e+00, 0.110000e+03,
			0.100000e+02, 0.000000e+00, 0.000000e+00
			},
			{ 0,
			0.859400e+11, 0.540000e+00, 0.105000e+03, -0.800000e+01, 0.110000e+03,
			0.100000e+02, 0.900000e+02, 0.200000e+01
			},
			{ 0,
			0.281000e+12, 0.000000e+00, 0.105000e+03, 0.280000e+02, 0.289500e+02,
			0.000000e+00, 0.000000e+00, 0.000000e+00
			},
			{ 0,
			0.330000e+11, 0.268270e+00, 0.105000e+03, 0.000000e+00, 0.110000e+03,
			0.100000e+02, 0.000000e+00, 0.000000e+00
			},
			{ 0,
			0.133000e+10, 0.119615e-01, 0.105000e+03, 0.000000e+00, 0.110000e+03,
			0.100000e+02, 0.000000e+00, 0.000000e+00
			},
			{ 0,
			0.176100e+06, 0.100000e+01, 0.950000e+02, -0.800000e+01, 0.110000e+03,
			0.100000e+02, 0.900000e+02, 0.200000e+01
			},
			{ 0,
			0.100000e+08, 0.100000e+01, 0.105000e+03, -0.800000e+01, 0.110000e+03,
			0.100000e+02, 0.900000e+02, 0.200000e+01
			}
		};

		// now assign the record values
		memcpy(msis86r.parm.pt, pt, 151 * sizeof(double));
		memcpy(msis86r.parm.ps, ps, 151 * sizeof(double));
		memcpy(msis86r.lower.ptm, ptm, 9 * sizeof(double));

		// since c reads in rows first, while the fortran code assumed columns,
		// we do a loop

		for (i = 0; i <= 7; i++)
		for (j = 0; j <= 8; j++)
			msis86r.lower.pdm[j][i] = pdm[i][j];
		for (i = 0; i <= 7; i++)
		for (j = 0; j <= 150; j++)
			msis86r.parm.pd[j][i] = pd[i][j];
		for (i = 0; i <= 2; i++)
		for (j = 0; j <= 25; j++)
			msis86r.parm.pdl[j][i] = pdl[i][j];

		// cdav
		// these are the extra arrays needed so a whole array can be accessed, but in
		// reverse order from the standard [row,col] approach in the fortran.
		//    extras to aid function calls for a subarry
		memcpy(msis86r.parm.pddav, pd, 8 * 151 * sizeof(double));

	}
	
	// ---------------- MSIS-90 model ----------------
	/*    ******************* local routines ******************************       */
	void ghp6
		(
		msistype& msis90r,
		lpolytype& lpoly,
		fittype& fit,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
		double f107, double ap[8],
		double d[9], double t[3], double press
		);

	double globe6
		(
		cswtype& csw,
		lpolytype& lpoly,
		double yrd, double sec, double lat, double llong, double tloc, double f107a,
		double f107, double ap[8], double p[151]
		);

	double glob6s
		(
		lpolytype& lpoly,
		cswtype& csw,
		double p[101]
		);



	/* ----------------------------------------------------------------------
	*          neutral atmosphere empirical model from the surface to lower
	*            exosphere  msise90 [jgr, 96, 1159-1172, 1991]
	*           a.e.hedin 4/24/90;6/3/91[add save]
	*           2/11/93 correct switch initialization and mks calculation
	*           2/11/97 [aeh] cprrect error in ghp6 when using meter6(.true.)
	*             see void ghp6 to specify a pressure rather than
	*             altitude.
	*     input:
	*          iyd - year and day as yyyyddd or just ddd [day of year from 1 to 365]
	*          sec - ut[sec]
	*          alt - altitude[km]
	*          glat - geodetic latitude[deg]
	*          glong - geodetic longitude[deg]
	*          stl - local apparent solar time[hrs]
	*          f107a - 3 month average of f10.7 flux
	*          f107 - daily f10.7 flux for previous day
	*          ap - magnetic index[daily] or when sw[9] = -1. :
	*             - array containing:
	*               [1] daily ap
	*               [2] 3 hr ap index for current time
	*               [3] 3 hr ap index for 3 hrs before current time
	*               [4] 3 hr ap index for 6 hrs before current time
	*               [5] 3 hr ap index for 9 hrs before current time
	*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
	*                      to current time
	*               [7] average of eight 3 hr ap indicies from 36 to 59 hrs prior
	*                      to current time
	*          mass - mass number [only density for selected gas is
	*                   calculated.  mass 0 is temperature.  mass 48 for all.
	*     note:  ut, local time, and longitude are used independently in the
	*              model and are not of equal importance for every situation.
	*              for the most physically realistic calculation these three
	*              variables should be consistent [stl = sec/3600+glong/15].
	*              f107, f107a, and ap effects are not large below 80 km
	*              and these can be set to 150., 150., and 4. respectively.
	*     output:
	*          d[1] - he number density[cm-3]
	*          d[2] - o number density[cm-3]
	*          d[3] - n2 number density[cm-3]
	*          d[4] - o2 number density[cm-3]
	*          d[5] - ar number density[cm-3]
	*          d[6] - total mass density[gm/cm3]
	*          d[7] - h number density[cm-3]
	*          d[8] - n number density[cm-3]
	*          t[1] - exospheric temperature
	*          t[2] - temperature at alt
	*
	*      to get output in m-3 and kg/m3:   call meters[.true.]
	*
	*      o, h, and n set to zero below 72.5 km
	*      exospheric temperature set to average for altitudes below 120 km.
	*
	*             the following is for test and special purposes:
	*              to turn on and off particular variations call tselec[sw]
	*                 where sw is a 25 element array containing 0. for off, 1.
	*                 for on, or 2. for main effects off but cross terms on
	*                 for the following variations
	*                 1 - f10.7 effect on mean  2 - time independent
	*                 3 - symmetrical annual    4 - symmetrical semiannual
	*                 5 - asymmetrical annual   6 - asymmetrical semiannual
	*                 7 - diurnal               8 - semidiurnal
	*                 9 - daily ap             10 - all ut/long effects
	*                11 - longitudinal         12 - ut and mixed ut/long
	*                13 - mixed ap/ut/long     14 - terdiurnal
	*                15 - departures from diffusive equilibrium
	*                16 - all tinf var         17 - all tlb var
	*                18 - all tn1 var           19 - all s var
	*                20 - all tn2 var           21 - all nlb var
	*                22 - all tn3 var           23 - turbo scale height var
	*
	*                to get current values of sw: call tretrv[sw]
	*
	* ----------------------------------------------------------------------      */

	void gtd6
		(
		msistype& msis90r,
		lpolytype& lpoly,
		fittype& fit,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong, double stl,
		double f107a, double f107, double ap[8], int mass,
		double d[9], double t[3]
		)
	{
		double d6[9], ts[3], v1, xlat, xmm, altt, dm28m, dmc, dz28, dmr, tz;
		int i, mss;

		// make sure and initilize the data through msis00init at the start

		int mn3 = 5;
		double zn3[6] = { 0, 32.5, 20.0, 15.0, 10.0, 0.0 };
		int mn2 = 4;
		double zn2[5] = { 0, 72.5, 55.0, 45.0, 32.5 };
		double zmix = 62.5;

		double mssl = -999;
		int sv[26] = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

		// treat as local, but save in between calls
		static double alast = 99999.0;

		//------------------------------ begin ---------------------------------
		if (msis90r.csw.isw != 64999)
			tselec(msis90r.csw, sv);

		// ----  test for changed input
		v1 = vtst(msis90r.csw, iyd, sec, glat, glong, stl, f107a, f107, ap, 1);

		// ---- latitude variation of gravity (none for msis90r.csw.sw[2]=0)
		xlat = glat;
		if (msis90r.csw.sw[2] == 0)
			xlat = 45.0;
		glatf(xlat, msis90r.parmb.gsurf, msis90r.parmb.re);

		xmm = msis90r.lower.pdm[5][3];

		// ---- thermosphere/mesosphere (above zn2[1])
		if (alt >= zn2[1])
			altt = alt;
		else
			altt = zn2[1];

		mss = mass;
		// ---- only calculate n2 in thermosphere if alt in mixed region
		if (alt < zmix && mass > 0)
			mss = 28;

		// ---- only calculate thermosphere if input parameters changed
		// ----   or altitude above zn2[1] in mesosphere
		if ((v1 == 1.0) | (alt > zn2[1]) | (alast > zn2[1]) | (mss != mssl))
		{
			gts6(msis90r, lpoly, lsqv, iyd, sec, altt, glat, glong, stl,
				f107a, f107, ap, mss, d6, ts);
			dm28m = msis90r.dmix.dm28;
			// ----   metric adjustment
			if (msis90r.metsel.imr == 1)
				dm28m = msis90r.dmix.dm28*1.0e6;
			mssl = mss;
		}
		t[1] = ts[1];
		t[2] = ts[2];
		if (alt >= zn2[1])
		{
			for (i = 1; i <= 9; i++)
				d[i] = d6[i];
			goto ten;
		}

		// ---- lower mesosphere/upper stratosphere [between zn3[1] and zn2[1]]
		// ----   temperature at nodes and gradients at end nodes
		// ----   inverse temperature a linear function of spherical harmonics
		// ----   only calculate nodes if input changed
		if ((v1 == 1.0) | (alast >= zn2[1]))
		{
			msis90r.meso.tgn2[1] = msis90r.meso.tgn1[2];
			msis90r.meso.tn2[1] = msis90r.meso.tn1[5];
			msis90r.meso.tn2[2] = msis90r.parm.pma[1][1] * msis90r.mavg.pavgm[1] /
				(1.0 - msis90r.csw.sw[20] * glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[1]));
			msis90r.meso.tn2[3] = msis90r.parm.pma[1][2] * msis90r.mavg.pavgm[2] /
				(1.0 - msis90r.csw.sw[20] * glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[2]));
			msis90r.meso.tn2[4] = msis90r.parm.pma[1][3] * msis90r.mavg.pavgm[3] /
				(1.0 - msis90r.csw.sw[20] * msis90r.csw.sw[22] *
				glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[3]));
			msis90r.meso.tgn2[2] = msis90r.mavg.pavgm[9] * msis90r.parm.pma[1][10] *
				(1.0 + msis90r.csw.sw[20] * msis90r.csw.sw[22] *
				glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[10]))
				*msis90r.meso.tn2[4] * msis90r.meso.tn2[4] /
				pow((msis90r.parm.pma[1][3] * msis90r.mavg.pavgm[3]), 2);
			msis90r.meso.tn3[1] = msis90r.meso.tn2[4];
		}
		if (alt < zn3[1])
		{
			// ---- lower stratosphere and troposphere [below zn3[1]]
			// ----   temperature at nodes and gradients at end nodes
			// ----   inverse temperature a linear function of spherical harmonics
			// ----   only calculate nodes if input changed
			if ((v1 == 1.0) | (alast >= zn3[1]))
			{
				msis90r.meso.tgn3[1] = msis90r.meso.tgn2[2];
				msis90r.meso.tn3[2] = msis90r.parm.pma[1][4] * msis90r.mavg.pavgm[4]
					/ (1.0 - msis90r.csw.sw[22] * glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[4]));
				msis90r.meso.tn3[3] = msis90r.parm.pma[1][5] * msis90r.mavg.pavgm[5]
					/ (1.0 - msis90r.csw.sw[22] * glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[5]));
				msis90r.meso.tn3[4] = msis90r.parm.pma[1][6] * msis90r.mavg.pavgm[6]
					/ (1.0 - msis90r.csw.sw[22] * glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[6]));
				msis90r.meso.tn3[5] = msis90r.parm.pma[1][7] * msis90r.mavg.pavgm[7]
					/ (1.0 - msis90r.csw.sw[22] * glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[7]));
				msis90r.meso.tgn3[2] = msis90r.parm.pma[1][8] * msis90r.mavg.pavgm[8] *
					(1.0 + msis90r.csw.sw[22] * glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[8]))
					*msis90r.meso.tn3[5] * msis90r.meso.tn3[5] /
					pow((msis90r.parm.pma[1][7] * msis90r.mavg.pavgm[7]), 2);
			}
		}

		if (mass == 0)
			goto fifty;
		// ----    linear transition to full mixing below zn2[1]
		dmc = 0;
		if (alt > zmix)
			dmc = 1.0 - (zn2[1] - alt) / (zn2[1] - zmix);
		dz28 = d6[3];

		// ----***** n2 density ****
		dmr = d6[3] / dm28m - 1.0;
		d[3] = densm(msis90r.parmb, fit, lsqv, alt, dm28m, xmm, tz, mn3, zn3,
			msis90r.meso.tn3, msis90r.meso.tgn3, mn2, zn2, msis90r.meso.tn2,
			msis90r.meso.tgn2);
		d[3] = d[3] * (1.0 + dmr*dmc);

		// ----***** he density ****
		d[1] = 0.0;
		if ((mass == 4) | (mass == 48))
		{
			dmr = d6[1] / (dz28*msis90r.lower.pdm[2][1]) - 1.0;
			d[1] = d[3] * msis90r.lower.pdm[2][1] * (1.0 + dmr*dmc);
		}

		// ----**** o density ****
		d[2] = 0.0;

		// ----***** o2 density ****
		d[4] = 0.0;
		if ((mass == 32) | (mass == 48))
		{
			dmr = d6[4] / (dz28*msis90r.lower.pdm[2][4]) - 1.0;
			d[4] = d[3] * msis90r.lower.pdm[2][4] * (1.0 + dmr*dmc);
		}

		// ----***** ar density ****
		d[5] = 0.0;
		if ((mass == 40) | (mass == 48))
		{
			dmr = d6[5] / (dz28*msis90r.lower.pdm[2][5]) - 1.0;
			d[5] = d[3] * msis90r.lower.pdm[2][5] * (1.0 + dmr*dmc);
		}

		// ----***** hydrogen density ****
		d[7] = 0.0;

		// ----***** atomic nitrogen density ****
		d[8] = 0.0;

		// ---- total mass density

		if (mass == 48)
		{
			d[6] = 1.66e-24*(4.0*d[1] + 16.0*d[2] + 28.0*d[3] + 32.0*d[4] + 40.0*d[5] + d[7] + 14.0*d[8]);
			if (msis90r.metsel.imr == 1)
				d[6] = d[6] / 1000.0;
		}
		t[2] = tz;
	ten:
		goto ninety;
	fifty:
		msis90r.gts3c.dd = densm(msis90r.parmb, fit, lsqv,
			alt, 1.0, 0.0, tz, mn3, zn3, msis90r.meso.tn3, msis90r.meso.tgn3, mn2, zn2,
			msis90r.meso.tn2, msis90r.meso.tgn2);
		t[2] = tz;
	ninety:
		alast = alt;
	}

	/* ----------------------------------------------------------------------
	*         find altitude of pressure surface [press] from gtd6
	*     input:
	*          iyd - year and day as yyyyddd
	*          sec - ut[sec]
	*          glat - geodetic latitude[deg]
	*          glong - geodetic longitude[deg]
	*          stl - local apparent solar time[hrs]
	*          f107a - 3 month average of f10.7 flux
	*          f107 - daily f10.7 flux for previous day
	*          ap - magnetic index[daily] or when sw[9] = -1.0 :
	*             - array containing:
	*               [1] daily ap
	*               [2] 3 hr ap index for current time
	*               [3] 3 hr ap index for 3 hrs before current time
	*               [4] 3 hr ap index for 6 hrs before current time
	*               [5] 3 hr ap index for 9 hrs before current time
	*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
	*                      to current time
	*               [7] average of eight 3 hr ap indicies from 36 to 59 hrs prior
	*                      to current time
	*          press - pressure level[mb]
	*     output:
	*          alt - altitude[km]
	*          d[1] - he number density[cm-3]
	*          d[2] - o number density[cm-3]
	*          d[3] - n2 number density[cm-3]
	*          d[4] - o2 number density[cm-3]
	*          d[5] - ar number density[cm-3]
	*          d[6] - total mass density[gm/cm3]
	*          d[7] - h number density[cm-3]
	*          d[8] - n number density[cm-3]
	*          t[1] - exospheric temperature
	*          t[2] - temperature at alt
	*
	* ----------------------------------------------------------------------      */

	void ghp6
		(
		msistype& msis90r,
		lpolytype& lpoly,
		fittype& fit,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
		double f107, double ap[8],
		double d[9], double t[3], double press
		)
	{
		double pl, zi, cl, cl2, cd, ca, z, l, xm, g, xn, p, sh, diff;
		int iday;

		double bm = 1.3806e-19;
		double rgas = 831.4;
		double test = 0.00043;

		//------------------------------ begin ---------------------------------
		pl = log10(press);

		//       initial altitude estimate
		if (pl >= -5.0)
		{
			if (pl > 2.5)   zi = 18.06*(3.00 - pl);
			if (pl > 0.75 &&  pl <= 2.5) zi = 14.98*(3.08 - pl);
			if (pl > -1.0 &&  pl <= 0.75) zi = 17.8*(2.72 - pl);
			if (pl > -2.0 &&  pl <= -1.0) zi = 14.28*(3.64 - pl);
			if (pl > -4.0 &&  pl <= -2.0) zi = 12.72*(4.32 - pl);
			if (pl <= -4.0) zi = 25.3*(0.11 - pl);
			iday = iyd - int(iyd / 1000.0)*1000.0;
			printf("not iday %12i \n", iday);
			iday = iyd % 1000;
			printf("not iday %12i \n", iday);
			cl = glat / 90.0;
			cl2 = cl*cl;
			if (iday < 182)
				cd = 1.0 - iday / 91.25;
			if (iday >= 182)
				cd = iday / 91.25 - 3.0;
			ca = 0;
			if (pl > -1.11 && pl <= -.23)
				ca = 1.0;
			if (pl > -0.23)
				ca = (2.79 - pl) / (2.79 + 0.23);
			if (pl <= -1.11 && pl > -3.0)
				ca = (-2.93 - pl) / (-2.93 + 1.11);
			z = zi - 4.87*cl*cd*ca - 1.64*cl2*ca + 0.31*ca*cl;
		}
		if (pl < -5.0)
			z = 22.0*pow(pl + 4.0, 2) + 110.0;

		//       iteration loop
		l = 0;
	ten:
		l = l + 1;
		gtd6(msis90r, lpoly, fit, lsqv,
			iyd, sec, z, glat, glong, stl, f107a, f107, ap, 48, d, t);
		xn = d[1] + d[2] + d[3] + d[4] + d[5] + d[7] + d[8];
		p = bm*xn*t[2];
		if (msis90r.metsel.imr == 1)
			p = p*1.0e-6;
		diff = pl - log10(p);
		if ((fabs(diff) < test) | (l == 6))
			goto twenty;
		xm = d[6] / xn / 1.66e-24;
		if (msis90r.metsel.imr == 1)
			xm = xm*1.0e3;
		g = msis90r.parmb.gsurf / pow((1.0 + z / msis90r.parmb.re), 2);
		sh = rgas*t[2] / (xm*g);

		// ----   new altitude estimate using scale height
		if (l < 6)
			z = z - sh*diff*2.302;
		else
			z = z - sh*diff;
		goto ten;
	twenty:
		if (l == 6)
			printf("not converging for press %12.2f %12.2f \n", press, diff);
		alt = z;
	}

	/* ----------------------------------------------------------------------
	*          neutral thermosphere model above 72.5 km for msise-90
	*           a.e.hedin 3/9/90
	*           coefficients not changed for 120km and above, but results may differ
	*          by a few percent from msis-86 [gts5] with introduction of a
	*          latitude dependent accel. of gravity.
	*           lower thermosphere reformulated for better continuation into
	*          lower atmosphere.
	*          for efficiency:
	*           exospheric temperature left at average value for alt below 120km;
	*           120 km gradient left at average value for alt below 72 km;
	*     input:
	*          iyd - year and day as yyyyddd
	*          sec - ut[sec]
	*          alt - altitude[km] [greater than 72.5 km]
	*          glat - geodetic latitude[deg]
	*          glong - geodetic longitude[deg]
	*          stl - local apparent solar time[hrs]
	*          f107a - 3 month average of f10.7 flux
	*          f107 - daily f10.7 flux for previous day
	*          ap - magnetic index[daily] or when sw[9] = -1.0 :
	*             - array containing:
	*               [1] daily ap
	*               [2] 3 hr ap index for current time
	*               [3] 3 hr ap index for 3 hrs before current time
	*               [4] 3 hr ap index for 6 hrs before current time
	*               [5] 3 hr ap index for 9 hrs before current time
	*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
	*                      to current time
	*               [7] average of eight 3 hr ap indicies from 36 to 59 hrs prior
	*                      to current time
	*          mass - mass number [only density for selected gas is
	*                   calculated.  mass 0 is temperature.  mass 48 for all.
	*     note:  ut, local time, and longitude are used independently in the
	*              model and are not of equal importance for every situation.
	*              for the most physically realistic calculation these three
	*              variables should be consistent [stl = sec/3600+glong/15].
	*     output:
	*          d[1] - he number density[cm-3]
	*          d[2] - o number density[cm-3]
	*          d[3] - n2 number density[cm-3]
	*          d[4] - o2 number density[cm-3]
	*          d[5] - ar number density[cm-3]
	*          d[6] - total mass density[gm/cm3]
	*          d[7] - h number density[cm-3]
	*          d[8] - n number density[cm-3]
	*          t[1] - exospheric temperature
	*          t[2] - temperature at alt
	*
	*             the following is for test and special purposes:
	*             [1] lower bound quantities in common/gts3c/
	*             [2] to turn on and off particular variations tselec[sw]
	*                 where sw is a 25 element array containing 0.0 for off, 1.0
	*                 for on, or 2.0 for main effects off but cross terms on
	*                 for the following variations
	*                 1 - f10.7 effect on mean  2 - time independent
	*                 3 - symmetrical annual    4 - symmetrical semiannual
	*                 5 - asymmetrical annual   6 - asymmetrical semiannual
	*                 7 - diurnal               8 - semidiurnal
	*                 9 - daily ap             10 - all ut/long effects
	*                11 - longitudinal         12 - ut and mixed ut/long
	*                13 - mixed ap/ut/long     14 - terdiurnal
	*                15 - departures from diffusive equilibrium
	*                16 - all tinf var         17 - all tlb var
	*                18 - all tn1 var           19 - all s var
	*                20 - all tn2 var           21 - all nlb var
	*                22 - all tn3 var           23 - turbo scale height var
	*
	*                to get current values of sw: tretrv[sw]
	*
	* ----------------------------------------------------------------------      */

	void gts6
		(
		msistype& msis90r,
		lpolytype& lpoly,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
		double f107, double ap[8], int mass, double d[10], double t[3]
		)
	{
		double v2, yrd, tinf, zlb,
			b01, b04, b14, b16, b28, b32, b40,
			g1, g14, g16, g28, g32, g4, g40,
			hc01, hc04, hc14, hc16, hc32, hc40,
			hcc01, hcc14, hcc16,
			rc01, rc14, rc16, tz, xmd, xmm,
			zc01, zc04, zc14, zc16, zc32, zc40,
			zcc01, zcc14, zcc16,
			zh01, zh04, zh14, zh16, zh28, zh32, zh40, zhf,
			zhm01, zhm04, zhm14, zhm16, zhm28, zhm32, zhm40;
		int i;
		double mt[11] = { 0, 48, 0, 4, 16, 28, 32, 40, 1, 49, 14 };
		int mn1 = 5;
		double zn1[6] = { 0, 120.0, 110.0, 100.0, 90.0, 72.5 };
		double dgtr = 1.74533e-2;
		double dr = 1.72142e-2;
		double altl[9] = { 0, 200.0, 400.0, 160.0, 200.0, 240.0, 450.0, 320.0, 450.0 };

		// treat as local, but save in between calls
		static double alast = -999.0;

		//------------------------------ begin ---------------------------------
		// ----  test for changed input
		v2 = vtst(msis90r.csw, iyd, sec, glat, glong, stl, f107a, f107, ap, 2);

		yrd = iyd;
		msis90r.gts3c.za = msis90r.parm.pdl[16][2];
		zn1[1] = msis90r.gts3c.za;
		for (i = 1; i <= 8; i++)
			d[i] = 0.0;

		// ----  tinf variations not important below za or zn1[1]
		if (alt > zn1[1])
		{
			if ((v2 == 1.0) | (alast <= zn1[1]))
				tinf = msis90r.lower.ptm[1] * msis90r.parm.pt[1] *
				(1.0 + msis90r.csw.sw[16] *
				globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107,
				ap, msis90r.parm.pt));
		}
		else
			tinf = msis90r.lower.ptm[1] * msis90r.parm.pt[1];

		t[1] = tinf;
		// ----    gradient variations not important below zn1[5]
		if (alt > zn1[5])
		{
			if ((v2 == 1) | (alast <= zn1[5]))
				msis90r.gts3c.g0 = msis90r.lower.ptm[4] * msis90r.parm.ps[1] * (1.0 + msis90r.csw.sw[19] *
				globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107,
				ap, msis90r.parm.ps));
		}
		else
			msis90r.gts3c.g0 = msis90r.lower.ptm[4] * msis90r.parm.ps[1];

		//  calculate these temperatures only if input changed
		if (v2 == 1.0)
			msis90r.gts3c.tlb = msis90r.lower.ptm[2] *
			(1.0 + msis90r.csw.sw[17] *
			globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107,
			ap, msis90r.parm.pddav[4])) *
			msis90r.parm.pd[1][4];

		msis90r.gts3c.s = msis90r.gts3c.g0 / (tinf - msis90r.gts3c.tlb);
		// ---- lower thermosphere temp variations not significant for
		// ----  density above 300 km
		if (alt < 300.0)
		{
			if ((v2 == 1.0) | (alast >= 300.0))
			{
				msis90r.meso.tn1[2] = msis90r.lower.ptm[7] * msis90r.parm.ptl[1][1] /
					(1.0 - msis90r.csw.sw[18] * glob6s(lpoly, msis90r.csw, msis90r.parm.ptldav[1]));
				msis90r.meso.tn1[3] = msis90r.lower.ptm[3] * msis90r.parm.ptl[1][2] /
					(1.0 - msis90r.csw.sw[18] * glob6s(lpoly, msis90r.csw, msis90r.parm.ptldav[2]));
				msis90r.meso.tn1[4] = msis90r.lower.ptm[8] * msis90r.parm.ptl[1][3] /
					(1.0 - msis90r.csw.sw[18] * glob6s(lpoly, msis90r.csw, msis90r.parm.ptldav[3]));
				msis90r.meso.tn1[5] = msis90r.lower.ptm[5] * msis90r.parm.ptl[1][4] /
					(1.0 - msis90r.csw.sw[18] * msis90r.csw.sw[20] *
					glob6s(lpoly, msis90r.csw, msis90r.parm.ptldav[4]));
				msis90r.meso.tgn1[2] = msis90r.lower.ptm[9] * msis90r.parm.pma[1][9] *
					(1.0 + msis90r.csw.sw[18] * msis90r.csw.sw[20] *
					glob6s(lpoly, msis90r.csw, msis90r.parm.pmadav[9])) *
					msis90r.meso.tn1[5] * msis90r.meso.tn1[5] /
					pow((msis90r.lower.ptm[5] * msis90r.parm.ptl[1][4]), 2);
			}
		}
		else
		{
			msis90r.meso.tn1[2] = msis90r.lower.ptm[7] * msis90r.parm.ptl[1][1];
			msis90r.meso.tn1[3] = msis90r.lower.ptm[3] * msis90r.parm.ptl[1][2];
			msis90r.meso.tn1[4] = msis90r.lower.ptm[8] * msis90r.parm.ptl[1][3];
			msis90r.meso.tn1[5] = msis90r.lower.ptm[5] * msis90r.parm.ptl[1][4];
			msis90r.meso.tgn1[2] = msis90r.lower.ptm[9] * msis90r.parm.pma[1][9] * msis90r.meso.tn1[5] *
				msis90r.meso.tn1[5] /
				pow((msis90r.lower.ptm[5] * msis90r.parm.ptl[1][4]), 2);
		}

		msis90r.gts3c.z0 = zn1[4];
		msis90r.gts3c.t0 = msis90r.meso.tn1[4];
		zlb = msis90r.lower.ptm[6];
		msis90r.gts3c.tr12 = 1.0;

		if (mass == 0)
			goto fifty;

		// ---- n2 variation factor at zlb
		g28 = msis90r.csw.sw[21] *
			globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap, msis90r.parm.pddav[3]);
		lpoly.day = int(yrd) % 1000;

		// ----  variation of turbopause height
		zhf = msis90r.parm.pdl[25][2] * (1.0 + msis90r.csw.sw[5] * msis90r.parm.pdl[25][1] * sin(dgtr*glat)*
			cos(dr*(lpoly.day - msis90r.parm.pt[14])));
		yrd = iyd;
		t[1] = tinf;
		xmm = msis90r.lower.pdm[5][3];

		for (i = 1; i <= 10; i++)
		if (mass == mt[i])   goto fifteen;

		printf("mass %i5 not valid\n", mass);
		goto ninety;

	fifteen:
		if ((alt <= altl[6]) | (mass == 28) | (mass == 48))
		{
			// ---- **** n2 density ****
			//       diffusive density at zlb
			//dav note dexp may be different in f90 and f77
			msis90r.gts3c.db28 = msis90r.lower.pdm[1][3] * exp(g28)*msis90r.parm.pd[1][3];

			//       diffusive density at alt
			d[3] = densu(msis90r.parmb, lsqv, alt, msis90r.gts3c.db28, tinf,
				msis90r.gts3c.tlb, 28.0, 0.0, t[2], zlb, msis90r.gts3c.s, mn1,
				zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
			msis90r.gts3c.dd = d[3];

			//       turbopause
			zh28 = msis90r.lower.pdm[3][3] * zhf;
			zhm28 = msis90r.lower.pdm[4][3] * msis90r.parm.pdl[6][2];
			xmd = 28.0 - xmm;

			//       mixed density at zlb
			b28 = densu(msis90r.parmb, lsqv, zh28, msis90r.gts3c.db28, tinf,
				msis90r.gts3c.tlb, xmd, -1.0, tz, zlb, msis90r.gts3c.s,
				mn1, zn1, msis90r.meso.tn1, msis90r.meso.tgn1);

			if ((alt <= altl[3]) && (msis90r.csw.sw[15] != 0))
			{
				//       mixed density at alt
				msis90r.dmix.dm28 = densu(msis90r.parmb, lsqv, alt, b28, tinf,
					msis90r.gts3c.tlb, xmm, 0.0, tz, zlb, msis90r.gts3c.s,
					mn1, zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
				//       net density at alt
				d[3] = dnet(d[3], msis90r.dmix.dm28, zhm28, xmm, 28.0);
			}
		}

		switch (i)
		{
		case 1: goto twenty;
		case 2: goto fifty;
		case 3: goto twenty;
		case 4: goto twentyfive;
		case 5: goto ninety;
		case 6: goto thirtyfive;
		case 7: goto fourty;
		case 8: goto fourtyfive;
		case 9: goto twentyfive;
		case 10: goto fourtyeight;
		}
	twenty:
		// ---- **** he density ****
		// ---- density variation factor at zlb
		g4 = msis90r.csw.sw[21] * globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis90r.parm.pddav[1]);

		//       diffusive density at zlb
		msis90r.gts3c.db04 = msis90r.lower.pdm[1][1] * exp(g4)*msis90r.parm.pd[1][1];

		//       diffusive density at alt
		d[1] = densu(msis90r.parmb, lsqv, alt, msis90r.gts3c.db04, tinf,
			msis90r.gts3c.tlb, 4.0, -0.4, t[2], zlb, msis90r.gts3c.s, mn1,
			zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
		msis90r.gts3c.dd = d[1];

		if (alt <= altl[1] && msis90r.csw.sw[15] != 0)
		{
			//       turbopause
			zh04 = msis90r.lower.pdm[3][1];
			zhm04 = zhm28;

			//       mixed density at zlb
			b04 = densu(msis90r.parmb, lsqv, zh04, msis90r.gts3c.db04, tinf,
				msis90r.gts3c.tlb, 4.0 - xmm, -1.4, t[2], zlb, msis90r.gts3c.s,
				mn1, zn1, msis90r.meso.tn1, msis90r.meso.tgn1);

			//       mixed density at alt
			msis90r.dmix.dm04 = densu(msis90r.parmb, lsqv, alt, b04, tinf, msis90r.gts3c.tlb,
				xmm, 0.0, t[2], zlb, msis90r.gts3c.s, mn1, zn1,
				msis90r.meso.tn1, msis90r.meso.tgn1);

			//       net density at alt
			d[1] = dnet(d[1], msis90r.dmix.dm04, zhm04, xmm, 4.0);

			//       correction to specified mixing ratio at ground
			msis90r.gts3c.rl = log(b28*msis90r.lower.pdm[2][1] / b04);
			zc04 = msis90r.lower.pdm[5][1] * msis90r.parm.pdl[1][2];
			hc04 = msis90r.lower.pdm[6][1] * msis90r.parm.pdl[2][2];

			//       net density corrected at alt
			d[1] = d[1] * ccor(alt, msis90r.gts3c.rl, hc04, zc04);
		}

		if (mass != 48)   goto ninety;

	twentyfive:
		// ----**** o density ****
		// ---- density variation factor at zlb
		g16 = msis90r.csw.sw[21] * globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107,
			ap, msis90r.parm.pddav[2]);

		//       diffusive density at zlb
		msis90r.gts3c.db16 = msis90r.lower.pdm[1][2] * exp(g16)*msis90r.parm.pd[1][2];

		// ---- diffusive density at alt
		d[2] = densu(msis90r.parmb, lsqv, alt, msis90r.gts3c.db16, tinf, msis90r.gts3c.tlb,
			16.0, 0.0, t[2], zlb, msis90r.gts3c.s,
			mn1, zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
		msis90r.gts3c.dd = d[2];
		if (alt <= altl[2] && msis90r.csw.sw[15] != 0)
		{
			//       corrected from msis90r.lower.pdm[3][1] to msis90r.lower.pdm[3][2]  12/2/85
			// ---- turbopause
			zh16 = msis90r.lower.pdm[3][2];
			zhm16 = zhm28;

			//       mixed density at zlb
			b16 = densu(msis90r.parmb, lsqv, zh16, msis90r.gts3c.db16, tinf, msis90r.gts3c.tlb, 16.0 - xmm, -1.0, t[2], zlb, msis90r.gts3c.s,
				mn1, zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
			//       mixed density at alt
			msis90r.dmix.dm16 = densu(msis90r.parmb, lsqv, alt, b16, tinf, msis90r.gts3c.tlb, xmm, 0.0, t[2], zlb, msis90r.gts3c.s, mn1,
				zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
			//       net density at alt
			d[2] = dnet(d[2], msis90r.dmix.dm16, zhm16, xmm, 16.0);

			// ---- correction to specified mixing ratio at ground
			msis90r.gts3c.rl = log(b28*msis90r.lower.pdm[2][2] * fabs(msis90r.parm.pdl[17][2]) / b16);
			hc16 = msis90r.lower.pdm[6][2] * msis90r.parm.pdl[4][2];
			zc16 = msis90r.lower.pdm[5][2] * msis90r.parm.pdl[3][2];
			d[2] = d[2] * ccor(alt, msis90r.gts3c.rl, hc16, zc16);

			// ---- chemistry correction
			hcc16 = msis90r.lower.pdm[8][2] * msis90r.parm.pdl[14][2];
			zcc16 = msis90r.lower.pdm[7][2] * msis90r.parm.pdl[13][2];
			rc16 = msis90r.lower.pdm[4][2] * msis90r.parm.pdl[15][2];

			//       net density corrected at alt
			d[2] = d[2] * ccor(alt, rc16, hcc16, zcc16);
		}
		if (mass != 48 && mass != 49) goto ninety;

	thirtyfive:
		// ---- **** o2 density ****
		// ---- density variation factor at zlb
		g32 = msis90r.csw.sw[21] * globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis90r.parm.pddav[5]);
		//       diffusive density at zlb
		msis90r.gts3c.db32 = msis90r.lower.pdm[1][4] * exp(g32)*msis90r.parm.pd[1][5];

		// ---- diffusive density at alt
		d[4] = densu(msis90r.parmb, lsqv, alt, msis90r.gts3c.db32, tinf, msis90r.gts3c.tlb, 32.0, 0.0, t[2], zlb, msis90r.gts3c.s, mn1, zn1,
			msis90r.meso.tn1, msis90r.meso.tgn1);
		if (mass == 49)
			msis90r.gts3c.dd = msis90r.gts3c.dd + 2.0*d[4];
		else
			msis90r.gts3c.dd = d[4];

		if (alt <= altl[4] && msis90r.csw.sw[15] != 0)
		{
			// ---- turbopause
			zh32 = msis90r.lower.pdm[3][4];
			zhm32 = zhm28;

			//       mixed density at zlb
			b32 = densu(msis90r.parmb, lsqv, zh32, msis90r.gts3c.db32, tinf, msis90r.gts3c.tlb, 32.0 - xmm, -1.0,
				t[2], zlb, msis90r.gts3c.s, mn1, zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
			//       mixed density at alt
			msis90r.dmix.dm32 = densu(msis90r.parmb, lsqv, alt, b32, tinf, msis90r.gts3c.tlb, xmm, 0.0, t[2], zlb, msis90r.gts3c.s, mn1, zn1,
				msis90r.meso.tn1, msis90r.meso.tgn1);
			//       net density at alt
			d[4] = dnet(d[4], msis90r.dmix.dm32, zhm32, xmm, 32.0);

			// ---- correction to specified mixing ratio at ground
			msis90r.gts3c.rl = log(b28*msis90r.lower.pdm[2][4] / b32);
			hc32 = msis90r.lower.pdm[6][4] * msis90r.parm.pdl[8][2];
			zc32 = msis90r.lower.pdm[5][4] * msis90r.parm.pdl[7][2];

			//       net density corrected at alt
			d[4] = d[4] * ccor(alt, msis90r.gts3c.rl, hc32, zc32);
		}
		if (mass != 48)   goto ninety;

	fourty:
		// ---- **** ar density ****
		// ---- density variation factor at zlb
		g40 = msis90r.csw.sw[21] * globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis90r.parm.pddav[6]);
		//       diffusive density at zlb
		msis90r.gts3c.db40 = msis90r.lower.pdm[1][5] * exp(g40)*msis90r.parm.pd[1][6];
		// ---- diffusive density at alt
		d[5] = densu(msis90r.parmb, lsqv, alt, msis90r.gts3c.db40, tinf, msis90r.gts3c.tlb, 40.0, 0.0, t[2], zlb, msis90r.gts3c.s, mn1,
			zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
		msis90r.gts3c.dd = d[5];

		if (alt <= altl[5] && msis90r.csw.sw[15] != 0)
		{
			// ---- turbopause
			zh40 = msis90r.lower.pdm[3][5];
			zhm40 = zhm28;

			//       mixed density at zlb
			b40 = densu(msis90r.parmb, lsqv, zh40, msis90r.gts3c.db40, tinf, msis90r.gts3c.tlb, 40.0 - xmm, -1.0,
				t[2], zlb, msis90r.gts3c.s, mn1, zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
			//       mixed density at alt
			msis90r.dmix.dm40 = densu(msis90r.parmb, lsqv, alt, b40, tinf, msis90r.gts3c.tlb, xmm, 0.0, t[2], zlb, msis90r.gts3c.s, mn1, zn1,
				msis90r.meso.tn1, msis90r.meso.tgn1);

			//       net density at alt
			d[5] = dnet(d[5], msis90r.dmix.dm40, zhm40, xmm, 40.0);

			// ---- correction to specified mixing ratio at ground
			msis90r.gts3c.rl = log(b28*msis90r.lower.pdm[2][5] / b40);
			hc40 = msis90r.lower.pdm[6][5] * msis90r.parm.pdl[10][2];
			zc40 = msis90r.lower.pdm[5][5] * msis90r.parm.pdl[9][2];

			//       net density corrected at alt
			d[5] = d[5] * ccor(alt, msis90r.gts3c.rl, hc40, zc40);
		}
		if (mass != 48)   goto ninety;

	fourtyfive:
		// ----  **** hydrogen density ****
		// ---- density variation factor at zlb
		g1 = msis90r.csw.sw[21] * globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis90r.parm.pddav[7]);

		//       diffusive density at zlb
		msis90r.gts3c.db01 = msis90r.lower.pdm[1][6] * exp(g1)*msis90r.parm.pd[1][7];
		// ---- diffusive density at alt
		d[7] = densu(msis90r.parmb, lsqv, alt, msis90r.gts3c.db01, tinf, msis90r.gts3c.tlb, 1.0, -0.4, t[2], zlb, msis90r.gts3c.s, mn1, zn1,
			msis90r.meso.tn1, msis90r.meso.tgn1);
		msis90r.gts3c.dd = d[7];

		if (alt <= altl[7] && msis90r.csw.sw[15] != 0)
		{
			// ---- turbopause
			zh01 = msis90r.lower.pdm[3][6];
			zhm01 = zhm28;

			//       mixed density at zlb
			b01 = densu(msis90r.parmb, lsqv, zh01, msis90r.gts3c.db01, tinf, msis90r.gts3c.tlb, 1.0 - xmm, -1.4,
				t[2], zlb, msis90r.gts3c.s, mn1, zn1, msis90r.meso.tn1, msis90r.meso.tgn1);

			//       mixed density at alt
			msis90r.dmix.dm01 = densu(msis90r.parmb, lsqv, alt, b01, tinf, msis90r.gts3c.tlb, xmm, 0.0, t[2], zlb, msis90r.gts3c.s, mn1, zn1,
				msis90r.meso.tn1, msis90r.meso.tgn1);

			//       net density at alt
			d[7] = dnet(d[7], msis90r.dmix.dm01, zhm01, xmm, 1.0);

			// ---- correction to specified mixing ratio at ground
			msis90r.gts3c.rl = log(b28*msis90r.lower.pdm[2][6] * fabs(msis90r.parm.pdl[18][2]) / b01);
			hc01 = msis90r.lower.pdm[6][6] * msis90r.parm.pdl[12][2];
			zc01 = msis90r.lower.pdm[5][6] * msis90r.parm.pdl[11][2];
			d[7] = d[7] * ccor(alt, msis90r.gts3c.rl, hc01, zc01);

			// ---- chemistry correction
			hcc01 = msis90r.lower.pdm[8][6] * msis90r.parm.pdl[20][2];
			zcc01 = msis90r.lower.pdm[7][6] * msis90r.parm.pdl[19][2];
			rc01 = msis90r.lower.pdm[4][6] * msis90r.parm.pdl[21][2];

			//       net density corrected at alt
			d[7] = d[7] * ccor(alt, rc01, hcc01, zcc01);
		}

	fourtyeight:
		// ----  **** atomic nitrogen density ****
		// ---- density variation factor at zlb
		g14 = msis90r.csw.sw[21] * globe6(msis90r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis90r.parm.pddav[8]);
		//       diffusive density at zlb
		msis90r.gts3c.db14 = msis90r.lower.pdm[1][7] * exp(g14)*msis90r.parm.pd[1][8];
		// ---- diffusive density at alt
		d[8] = densu(msis90r.parmb, lsqv, alt, msis90r.gts3c.db14, tinf, msis90r.gts3c.tlb, 14.0, 0.0, t[2], zlb, msis90r.gts3c.s, mn1,
			zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
		msis90r.gts3c.dd = d[8];

		if (alt <= altl[8] && msis90r.csw.sw[15] != 0)
		{
			// ---- turbopause
			zh14 = msis90r.lower.pdm[3][7];
			zhm14 = zhm28;

			//        mixed density at zlb
			b14 = densu(msis90r.parmb, lsqv, zh14, msis90r.gts3c.db14, tinf, msis90r.gts3c.tlb, 14.0 - xmm, -1.0,
				t[2], zlb, msis90r.gts3c.s, mn1, zn1, msis90r.meso.tn1, msis90r.meso.tgn1);

			//       mixed density at alt
			msis90r.dmix.dm14 = densu(msis90r.parmb, lsqv, alt, b14, tinf, msis90r.gts3c.tlb, xmm, 0.0, t[2], zlb, msis90r.gts3c.s, mn1,
				zn1, msis90r.meso.tn1, msis90r.meso.tgn1);
			//       net density at alt
			d[8] = dnet(d[8], msis90r.dmix.dm14, zhm14, xmm, 14.0);

			// ---- correction to specified mixing ratio at ground
			msis90r.gts3c.rl = log(b28*msis90r.lower.pdm[2][7] * fabs(msis90r.parm.pdl[3][1]) / b14);
			hc14 = msis90r.lower.pdm[6][7] * msis90r.parm.pdl[2][1];
			zc14 = msis90r.lower.pdm[5][7] * msis90r.parm.pdl[1][1];
			d[8] = d[8] * ccor(alt, msis90r.gts3c.rl, hc14, zc14);

			// ---- chemistry correction
			hcc14 = msis90r.lower.pdm[8][7] * msis90r.parm.pdl[5][1];
			zcc14 = msis90r.lower.pdm[7][7] * msis90r.parm.pdl[4][1];
			rc14 = msis90r.lower.pdm[4][7] * msis90r.parm.pdl[6][1];

			//       net density corrected at alt
			d[8] = d[8] * ccor(alt, rc14, hcc14, zcc14);
		}
		if (mass != 48) goto ninety;

		// ---- total mass density
		d[6] = 1.66e-24*(4.0*d[1] + 16.0*d[2] + 28.0*d[3] + 32.0*d[4] +
			40.0*d[5] + d[7] + 14.0*d[8]);
		msis90r.gts3c.db48 = 1.66e-24*(4.0*msis90r.gts3c.db04 + 16.0*msis90r.gts3c.db16 + 28.0*msis90r.gts3c.db28 + 32.0*msis90r.gts3c.db32 +
			40.0*msis90r.gts3c.db40 + msis90r.gts3c.db01 + 14.0*msis90r.gts3c.db14);
		goto ninety;

		// ---- temperature at altitude
	fifty:
		//cdav
		//        ddum is never used so this call is unnecessary
		//        ddum = densu(msis90r.parmb, lsqv, alt,1.0,tinf,msis90r.gts3c.tlb,0.0,0.0,t[2],zlb,msis90r.gts3c.s,mn1,
		//               zn1,msis90r.meso.tn1,msis90r.meso.tgn1);
		goto ninety;
	ninety:

		// ---- adjust densities from cgs to kgm
		if (msis90r.metsel.imr == 1)
		{
			for (i = 1; i <= 8; i++)
				d[i] = d[i] * 1.0e6;
			d[6] = d[6] / 1000.0;
		}
		alast = alt;
	}

	//// cdav these next few functions could be done as inline functions if that speeds
	//// things up
	//// ---- 3hr magnetic activity functions
	////      eq. a24d
	//      double g0
	//          ( double a, double p[151])
	//          {
	//          return (a - 4.0 + (p[26]-1.0) *
	//                 (a - 4.0 + (exp(-fabs(p[25])*(a-4.0)) -1.0) / fabs(p[25]) ));
	//          }
	//// ---- eq. a24c
	//      double sumex
	//          ( double ex )
	//          {
	//          return 1.0 + (1.0-pow(ex,19)) / (1.0-ex)*sqrt(ex);
	//          }
	//// ---- eq. a24a
	//      double sg0
	//          ( double ex, double p[151], double ap[8] )
	//          {
	//          return ( g0(ap[2],p) +
	//                  ( g0(ap[3],p)*ex + g0(ap[4],p)*ex*ex + g0(ap[5],p)* pow(ex,3) +
	//                   ( g0(ap[6],p)*pow(ex,4) + g0(ap[7],p)*pow(ex,12)) * (1.0-pow(ex,8))
	//                   / (1.0-ex)
	//                  )
	//                 ) / sumex(ex);
	//          }
	/* ----------------------------------------------------------------------
	*         calculate g[l] function
	*         upper thermosphere parameters
	* ----------------------------------------------------------------------      */

	double globe6
		(
		cswtype& csw,
		lpolytype& lpoly,
		double yrd, double sec, double lat, double llong, double tloc, double f107a,
		double f107, double ap[8], double p[151]
		)
	{
		int i, iyr;
		double c, s, c2, c4, s2, cd14, cd18, cd32, cd39, f1, f2, t71, t72,
			t81, t82, p44, p45, exp1;
		double t[16], tinfg;

		double dgtr = 1.74533e-2;
		double dr = 1.72142e-2;
		double xl = 1000.0;
		static double tll = 1000.0;
		static double sw9 = 1.0;
		double dayl = -1.0;
		double p14 = -1000.0;
		double p18 = -1000.0;
		double p32 = -1000.0;
		double hr = 0.2618;
		double sr = 7.2722e-5;
		int sv[26] = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
		double nsw = 14;
		double p39 = -1000.0;

		// are these needed? used from m90?
		double exp2;

		static double longl = -999.0;

		//------------------------------ begin ---------------------------------
		if (csw.isw != 64999)
			tselec(csw, sv);

		for (i = 1; i <= 15; i++)
			t[i] = 0.0;

		if (csw.sw[9] > 0)
			sw9 = 1.0;
		if (csw.sw[9] < 0)
			sw9 = -1.0;
		iyr = yrd / 1000.0;
		lpoly.day = yrd - iyr*1000.0;
		lpoly.xlong = llong;
		//  eq. a22 (remainder of code)
		if (xl != lat)
		{
			// ----    calculate legendre polynomials
			c = sin(lat*dgtr);
			s = cos(lat*dgtr);
			c2 = c*c;
			c4 = c2*c2;
			s2 = s*s;
			lpoly.plg[2][1] = c;
			lpoly.plg[3][1] = 0.5*(3.0*c2 - 1.0);
			lpoly.plg[4][1] = 0.5*(5.0*c*c2 - 3.0*c);
			lpoly.plg[5][1] = (35.0*c4 - 30.0*c2 + 3.0) / 8.0;
			lpoly.plg[6][1] = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c) / 8.0;
			lpoly.plg[7][1] = (11.0*c*lpoly.plg[6][1] - 5.0*lpoly.plg[5][1]) / 6.0;
			//        lpoly.plg[8][1] = (13.0*c*lpoly.plg[7][1] - 6.0*lpoly.plg[6][1])/7.0;
			lpoly.plg[2][2] = s;
			lpoly.plg[3][2] = 3.0*c*s;
			lpoly.plg[4][2] = 1.5*(5.0*c2 - 1.0)*s;
			lpoly.plg[5][2] = 2.5*(7.0*c2*c - 3.0*c)*s;
			lpoly.plg[6][2] = 1.875*(21.0*c4 - 14.0*c2 + 1.0)*s;
			lpoly.plg[7][2] = (11.0*c*lpoly.plg[6][2] - 6.0*lpoly.plg[5][2]) / 5.0;
			//        lpoly.plg[8][2] = (13.0*c*lpoly.plg[7][2]-7.0*lpoly.plg[6][2])/6.0;
			//        lpoly.plg[9][2] = (15.0*c*lpoly.plg[8][2]-8.0*lpoly.plg[7][2])/7.0;
			lpoly.plg[3][3] = 3.0*s2;
			lpoly.plg[4][3] = 15.0*s2*c;
			lpoly.plg[5][3] = 7.5*(7.0*c2 - 1.0)*s2;
			lpoly.plg[6][3] = 3.0*c*lpoly.plg[5][3] - 2.0*lpoly.plg[4][3];
			lpoly.plg[7][3] = (11.0*c*lpoly.plg[6][3] - 7.0*lpoly.plg[5][3]) / 4.0;
			lpoly.plg[8][3] = (13.0*c*lpoly.plg[7][3] - 8.0*lpoly.plg[6][3]) / 5.0;
			lpoly.plg[4][4] = 15.0*s2*s;
			lpoly.plg[5][4] = 105.0*s2*s*c;
			lpoly.plg[6][4] = (9.0*c*lpoly.plg[5][4] - 7.0*lpoly.plg[4][4]) / 2.0;
			lpoly.plg[7][4] = (11.0*c*lpoly.plg[6][4] - 8.0*lpoly.plg[5][4]) / 3.0;

			xl = lat;
		}
		if (tll == tloc) goto sixteen;
		if (csw.sw[7] == 0 && csw.sw[8] == 0 && csw.sw[14] == 0) goto sixteen;
		lpoly.stloc = sin(hr*tloc);
		lpoly.ctloc = cos(hr*tloc);
		lpoly.s2tloc = sin(2.0*hr*tloc);
		lpoly.c2tloc = cos(2.0*hr*tloc);
		lpoly.s3tloc = sin(3.0*hr*tloc);
		lpoly.c3tloc = cos(3.0*hr*tloc);
		tll = tloc;
	sixteen:
		if (llong != longl)
		{
			lpoly.clong = cos(dgtr*llong);
			lpoly.slong = sin(dgtr*llong);
		}
		longl = llong;
		if ((lpoly.day != dayl) | (p[14] != p14))
			cd14 = cos(dr*(lpoly.day - p[14]));
		if ((lpoly.day != dayl) | (p[18] != p18))
			cd18 = cos(2.0*dr*(lpoly.day - p[18]));

		if ((lpoly.day != dayl) | (p[32] != p32))
			cd32 = cos(dr*(lpoly.day - p[32]));
		if ((lpoly.day != dayl) | (p[39] != p39))
			cd39 = cos(2.0*dr*(lpoly.day - p[39]));

		dayl = lpoly.day;
		p14 = p[14];
		p18 = p[18];
		p32 = p[32];
		p39 = p[39];

		// ----   f10.7 effect
		lpoly.df = f107 - f107a;
		lpoly.dfa = f107a - 150.0;
		t[1] = p[20] * lpoly.df + p[21] * lpoly.df*lpoly.df +
			p[22] * lpoly.dfa + p[30] * lpoly.dfa*lpoly.dfa;
		f1 = 1.0 + (p[48] * lpoly.dfa + p[20] * lpoly.df + p[21] * lpoly.df*lpoly.df)*csw.swc[1];
		f2 = 1.0 + (p[50] * lpoly.dfa + p[20] * lpoly.df + p[21] * lpoly.df*lpoly.df)*csw.swc[1];

		// ----  time independent
		t[2] = (p[2] * lpoly.plg[3][1] + p[3] * lpoly.plg[5][1] + p[23] * lpoly.plg[7][1])
			+ (p[15] * lpoly.plg[3][1])*lpoly.dfa*csw.swc[1] + p[27] * lpoly.plg[2][1];

		// ----  symmetrical annual
		t[3] = (p[19])*cd32;

		// ----  symmetrical semiannual
		t[4] = (p[16] + p[17] * lpoly.plg[3][1])*cd18;

		// ----  asymmetrical annual
		t[5] = f1* (p[10] * lpoly.plg[2][1] + p[11] * lpoly.plg[4][1])*cd14;

		// ----   asymmetrical semiannual
		t[6] = p[38] * lpoly.plg[2][1] * cd39;

		// ----  diurnal
		if (csw.sw[7] != 0)
		{
			t71 = (p[12] * lpoly.plg[3][2])*cd14*csw.swc[5];
			t72 = (p[13] * lpoly.plg[3][2])*cd14*csw.swc[5];
			t[7] = f2* ((p[4] * lpoly.plg[2][2] + p[5] * lpoly.plg[4][2] + p[28] * lpoly.plg[6][2]
				+ t71)*lpoly.ctloc + (p[7] * lpoly.plg[2][2] + p[8] * lpoly.plg[4][2]
				+ p[29] * lpoly.plg[6][2] + t72)*lpoly.stloc);
		}
		// ----  semidiurnal
		if (csw.sw[8] != 0)
		{
			t81 = (p[24] * lpoly.plg[4][3] + p[36] * lpoly.plg[6][3])*cd14*csw.swc[5];
			t82 = (p[34] * lpoly.plg[4][3] + p[37] * lpoly.plg[6][3])*cd14*csw.swc[5];
			t[8] = f2*((p[6] * lpoly.plg[3][3] + p[42] * lpoly.plg[5][3] + t81)*lpoly.c2tloc
				+ (p[9] * lpoly.plg[3][3] + p[43] * lpoly.plg[5][3] + t82)*lpoly.s2tloc);
		}
		// ----  terdiurnal
		if (csw.sw[14] != 0)
		{
			t[14] = f2*
				((p[40] * lpoly.plg[4][4] + (p[94] * lpoly.plg[5][4]
				+ p[47] * lpoly.plg[7][4])*cd14*csw.swc[5])*lpoly.s3tloc
				+ (p[41] * lpoly.plg[4][4] + (p[95] * lpoly.plg[5][4]
				+ p[49] * lpoly.plg[7][4])*cd14*csw.swc[5])*lpoly.c3tloc);
		}

		// ----    magnetic activity based on daily ap
		if (sw9 != -1.0)
		{
			lpoly.apd = (ap[1] - 4.0);
			p44 = p[44];
			p45 = p[45];
			// cdav at one time i saw an error where p44 was less than zero - unfortunately, i didn't
			//      write exactly how i got there. it may have been during debugging when i had other
			//      things wrong
			if (p44 <= 0.0)
				p44 = 1.0e-5;
			lpoly.apdf = lpoly.apd + (p45 - 1.0)*(lpoly.apd + (exp(-p44  *lpoly.apd) - 1.0) / p44);
			if (csw.sw[9] == 0)
				goto fourty;
			t[9] = lpoly.apdf*(p[33] + p[46] * lpoly.plg[3][1] + p[35] * lpoly.plg[5][1] +
				(p[101] * lpoly.plg[2][1] + p[102] * lpoly.plg[4][1] + p[103] * lpoly.plg[6][1])*cd14*csw.swc[5] +
				(p[122] * lpoly.plg[2][2] + p[123] * lpoly.plg[4][2] + p[124] * lpoly.plg[6][2])*csw.swc[7] *
				cos(hr*(tloc - p[125])));
			goto fourty;
		}

		if (p[52] == 0)
			goto fourty;
		exp1 = exp(-10800.0*fabs(p[52]) / (1.0 + p[139] * (45.0 - fabs(lat))));
		exp2 = exp(-10800.0*fabs(p[54]));
		if (exp1 > .99999)
			exp1 = 0.99999;
		if (p[25] < 1.0e-4)
			p[25] = 1.0e-4;
		lpoly.apt[1] = sg0(exp1, p, ap);
		//        apt[2] = sg2[exp1]
		lpoly.apt[3] = sg0(exp2, p, ap);
		//        apt[4] = sg2[exp2]

		if (csw.sw[9] == 0) goto fourty;
		t[9] = lpoly.apt[1] * (p[51] + p[97] * lpoly.plg[3][1] + p[55] * lpoly.plg[5][1] +
			(p[126] * lpoly.plg[2][1] + p[127] * lpoly.plg[4][1] + p[128] * lpoly.plg[6][1])*cd14*csw.swc[5] +
			(p[129] * lpoly.plg[2][2] + p[130] * lpoly.plg[4][2] + p[131] * lpoly.plg[6][2])*csw.swc[7] *
			cos(hr*(tloc - p[132])));
	fourty:

		if ((csw.sw[10] == 0) | (llong <= -1000.0))
			goto fourtynine;

		// ----  longitudinal
		if (csw.sw[11] != 0)
		{
			t[11] = (1.0 + p[81] * lpoly.dfa*csw.swc[1])*
				((p[65] * lpoly.plg[3][2] + p[66] * lpoly.plg[5][2] + p[67] * lpoly.plg[7][2]
				+ p[104] * lpoly.plg[2][2] + p[105] * lpoly.plg[4][2] + p[106] * lpoly.plg[6][2]
				+ csw.swc[5] * (p[110] * lpoly.plg[2][2] + p[111] * lpoly.plg[4][2] + p[112] * lpoly.plg[6][2])*cd14)*
				cos(dgtr*llong)
				+ (p[91] * lpoly.plg[3][2] + p[92] * lpoly.plg[5][2] + p[93] * lpoly.plg[7][2]
				+ p[107] * lpoly.plg[2][2] + p[108] * lpoly.plg[4][2] + p[109] * lpoly.plg[6][2]
				+ csw.swc[5] * (p[113] * lpoly.plg[2][2] + p[114] * lpoly.plg[4][2] + p[115] * lpoly.plg[6][2])*cd14)*
				sin(dgtr*llong));
		}

		// ----  ut and mixed ut,longitude
		if (csw.sw[12] != 0)
		{
			t[12] = (1.0 + p[96] * lpoly.plg[2][1])*(1.0 + p[82] * lpoly.dfa*csw.swc[1])*
				(1.0 + p[120] * lpoly.plg[2][1] * csw.swc[5] * cd14)*
				((p[69] * lpoly.plg[2][1] + p[70] * lpoly.plg[4][1] + p[71] * lpoly.plg[6][1])*
				cos(sr*(sec - p[72])));
			t[12] = t[12] + csw.swc[11] *
				(p[77] * lpoly.plg[4][3] + p[78] * lpoly.plg[6][3] + p[79] * lpoly.plg[8][3])*
				cos(sr*(sec - p[80]) + 2.0*dgtr*llong)*(1.0 + p[138] * lpoly.dfa*
				csw.swc[1]);
		}

		// ----  ut,longitude magnetic activity
		if (csw.sw[13] == 0) goto fourtyeight;

		if (sw9 != -1.0)
		{
			t[13] = lpoly.apdf*csw.swc[11] * (1.0 + p[121] * lpoly.plg[2][1])*
				((p[61] * lpoly.plg[3][2] + p[62] * lpoly.plg[5][2] + p[63] * lpoly.plg[7][2])*
				cos(dgtr*(llong - p[64])))
				+ lpoly.apdf*csw.swc[11] * csw.swc[5] *
				(p[116] * lpoly.plg[2][2] + p[117] * lpoly.plg[4][2] + p[118] * lpoly.plg[6][2])*
				cd14*cos(dgtr*(llong - p[119]))
				+ lpoly.apdf*csw.swc[12] *
				(p[84] * lpoly.plg[2][1] + p[85] * lpoly.plg[4][1] + p[86] * lpoly.plg[6][1])*
				cos(sr*(sec - p[76]));
			goto fourtyeight;
		}

		if (p[52] == 0) goto fourtyeight;
		t[13] = lpoly.apt[1] * csw.swc[11] * (1.0 + p[133] * lpoly.plg[2][1])*
			((p[53] * lpoly.plg[3][2] + p[99] * lpoly.plg[5][2] + p[68] * lpoly.plg[7][2])*
			cos(dgtr*(llong - p[98])))
			+ lpoly.apt[1] * csw.swc[11] * csw.swc[5] *
			(p[134] * lpoly.plg[2][2] + p[135] * lpoly.plg[4][2] + p[136] * lpoly.plg[6][2])*
			cd14*cos(dgtr*(llong - p[137]))
			+ lpoly.apt[1] * csw.swc[12] *
			(p[56] * lpoly.plg[2][1] + p[57] * lpoly.plg[4][1] + p[58] * lpoly.plg[6][1])*
			cos(sr*(sec - p[59]));
	fourtyeight:
	fourtynine :
		tinfg = p[31];
			   for (i = 1; i <= nsw; i++)
				   tinfg = tinfg + fabsf(csw.sw[i])*t[i];

			   return tinfg;
	}

	/*
	//----------------------------------------------------------------------
	//      version of globe for lower atmosphere 1/17/90
	//dav chg to *8
	//----------------------------------------------------------------------
	*/
	double glob6s
		(
		lpolytype& lpoly,
		cswtype& csw,
		double p[101]
		)
	{
		double cd32, cd18, cd14, cd39, t[15], t71, t72, t81, t82, tt;

		int i;

		double dr = 1.72142e-2;
		double dayl = -1.0;
		double p32 = -1000.0;
		double p18 = -1000.0;
		double p14 = -1000.0;
		double p39 = -1000.0;

		//------------------------------ begin ---------------------------------
		for (i = 1; i <= 14; i++)
			t[i] = 0.0;

		if ((lpoly.day != dayl) | (p32 != p[32]))
			cd32 = cos(dr*(lpoly.day - p[32]));
		if ((lpoly.day != dayl) | (p18 != p[18]))
			cd18 = cos(2.0*dr*(lpoly.day - p[18]));
		if ((lpoly.day != dayl) | (p14 != p[14]))
			cd14 = cos(dr*(lpoly.day - p[14]));
		if ((lpoly.day != dayl) | (p39 != p[39]))
			cd39 = cos(2.0*dr*(lpoly.day - p[39]));

		dayl = lpoly.day;
		p32 = p[32];
		p18 = p[18];
		p14 = p[14];
		p39 = p[39];

		// ---- f10.7
		t[1] = p[22] * lpoly.dfa;

		// ---- time independent
		t[2] = p[2] * lpoly.plg[3][1] + p[3] * lpoly.plg[5][1] + p[23] * lpoly.plg[7][1]
			+ p[27] * lpoly.plg[2][1] + p[28] * lpoly.plg[4][1] + p[29] * lpoly.plg[6][1];

		// ---- symmetrical annual
		t[3] = (p[19] + p[48] * lpoly.plg[3][1] + p[30] * lpoly.plg[5][1])*cd32;

		// ---- symmetrical semiannual
		t[4] = (p[16] + p[17] * lpoly.plg[3][1] + p[31] * lpoly.plg[5][1])*cd18;

		// ---- asymmetrical annual
		t[5] = (p[10] * lpoly.plg[2][1] + p[11] * lpoly.plg[4][1] + p[36] * lpoly.plg[6][1])*cd14;

		// ---- asymmetrical semiannual
		t[6] = (p[38] * lpoly.plg[2][1])*cd39;

		// ----  diurnal
		if (csw.sw[7] != 0)
		{
			t71 = p[12] * lpoly.plg[3][2] * cd14*csw.swc[5];
			t72 = p[13] * lpoly.plg[3][2] * cd14*csw.swc[5];
			t[7] = ((p[4] * lpoly.plg[2][2] + p[5] * lpoly.plg[4][2] + t71)*lpoly.ctloc
				+ (p[7] * lpoly.plg[2][2] + p[8] * lpoly.plg[4][2] + t72)*lpoly.stloc);
		}

		// ----  semidiurnal
		if (csw.sw[8] != 0)
		{
			t81 = (p[24] * lpoly.plg[4][3] + p[47] * lpoly.plg[6][3])*cd14*csw.swc[5];
			t82 = (p[34] * lpoly.plg[4][3] + p[49] * lpoly.plg[6][3])*cd14*csw.swc[5];
			t[8] = ((p[6] * lpoly.plg[3][3] + p[42] * lpoly.plg[5][3] + t81)*lpoly.c2tloc
				+ (p[9] * lpoly.plg[3][3] + p[43] * lpoly.plg[5][3] + t82)*lpoly.s2tloc);
		}

		// ----  terdiurnal
		if (csw.sw[14] != 0)
			t[14] = p[40] * lpoly.plg[4][4] * lpoly.s3tloc + p[41] * lpoly.plg[4][4] * lpoly.c3tloc;

		// ---- magnetic activity
		if (csw.sw[9] != 0)
		{
			if (csw.sw[9] == 1)
				t[9] = lpoly.apdf*(p[33] + p[46] * lpoly.plg[3][1] * csw.swc[2]);
			if (csw.sw[9] == -1)
				t[9] = (p[51] * lpoly.apt[3] + p[97] * lpoly.plg[3][1] * lpoly.apt[3] * csw.swc[2]);
		}
		if ((csw.sw[10] == 0) | (csw.sw[11] == 0) | (lpoly.xlong <= -1000.0))
			goto fourtynine;

		// ----  longitudinal
		t[11] = (1.0 + lpoly.plg[2][1] * (p[81] * csw.swc[5] * cos(dr*(lpoly.day - p[82]))
			+ p[86] * csw.swc[6] * cos(2.0*dr*(lpoly.day - p[87])))
			+ p[84] * csw.swc[3] * cos(dr*(lpoly.day - p[85]))
			+ p[88] * csw.swc[4] * cos(2.0*dr*(lpoly.day - p[89])))
			*((p[65] * lpoly.plg[3][2] + p[66] * lpoly.plg[5][2] + p[67] * lpoly.plg[7][2]
			+ p[75] * lpoly.plg[2][2] + p[76] * lpoly.plg[4][2] + p[77] * lpoly.plg[6][2])*lpoly.clong
			+ (p[91] * lpoly.plg[3][2] + p[92] * lpoly.plg[5][2] + p[93] * lpoly.plg[7][2]
			+ p[78] * lpoly.plg[2][2] + p[79] * lpoly.plg[4][2] + p[80] * lpoly.plg[6][2])*lpoly.slong);
	fourtynine:
		tt = 0.0;
		for (i = 1; i <= 14; i++)
			tt = tt + fabsf(csw.sw[i])*t[i];

		return tt;
	}

	/*
	//----------------------------------------------------------------------
	//-    msise 90 12-mar-90
	//----------------------------------------------------------------------
	*/
	void msis90init
		(
		msistype& msis90r
		)
	{
		int i, j;

		// cdav i'm not string literate yet in c++ :-)
		//        msis90r.datime.isdate = "02-nov-97";
		//        msis90r.datime.istime = "15:49:27 ";
		//        msis90r.datime.name   = "msise-90 ";

		msis90r.csw.isw = 0;
		//
		// cdav   set this to output in
		msis90r.metsel.imr = 1;

		// ----   temperature
		static double pt[151] =
		{ 0,
		9.96040e-01, 3.85528e-02, 3.03445e-03, -1.05531e-01, -6.07134e-03,
		-5.16278e-04, -1.15622e-01, 2.02240e-03, 9.90156e-03, -1.27371e-01,
		-3.02449e-02, 1.23512e-02, -5.26277e-03, -8.45398e+00, 0.00000e+00,
		1.42370e-02, 0.00000e+00, 1.25818e+02, 8.05486e-03, 1.64419e-03,
		-6.21452e-06, 3.11701e-03, 0.00000e+00, 3.86578e-03, 1.32397e-01,
		2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, -6.41110e-06,
		0.00000e+00, 3.00150e+01, 5.33297e-03, 3.89146e-03, 2.04725e-03,
		0.00000e+00, 0.00000e+00, -1.92645e-02, 2.75905e+00, 1.47284e-03,
		3.41345e-04, -1.17388e-03, -3.54589e-04, 1.13139e-01, 1.69134e-01,
		5.08295e-03, 3.65016e-05, 4.26385e-03, 1.15102e-04, 5.11819e-03,
		6.09108e-03, 4.04995e-05, 1.53049e-03, 2.41470e-05, 2.30764e-03,
		1.55267e-03, 1.33722e-03, -1.82318e-03, -2.63007e+02, 0.00000e+00,
		1.37337e-03, 9.95774e-04, 0.00000e+00, -1.08983e+02, 5.62606e-03,
		5.94053e-03, 1.09358e-03, 0.00000e+00, -1.33410e-02, -2.43409e-02,
		-1.35688e-02, 3.11370e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-2.83023e+03, 8.45583e-04, 5.38706e-04, 0.00000e+00, 2.47956e+02,
		2.92246e-03, 0.00000e+00, 0.00000e+00, 7.47703e-05, 8.87993e-04,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-1.16540e-02, -4.49173e-03, -3.53189e-04, -1.73933e-04, -1.53218e-04,
		-5.65411e-01, 7.77272e-03, -9.11784e+01, 6.45187e-04, 0.00000e+00,
		-8.37685e-04, 2.42318e-03, 4.73796e-03, -3.01801e-03, -4.23564e-03,
		-2.48289e-03, 9.19286e-04, 2.16372e-03, 8.63968e-04, 1.89689e-03,
		4.15654e-03, 0.00000e+00, 1.18068e-02, 3.31190e-03, 0.00000e+00,
		1.20222e-03, 0.00000e+00, 0.00000e+00, -3.07246e+00, 0.00000e+00,
		0.00000e+00, 6.72403e-04, 1.08930e-03, 9.72278e-04, 4.68242e+00,
		-3.15034e-04, 4.00059e-03, 5.15036e-03, 1.62989e-03, 1.08824e-03,
		9.95261e-04, 4.18955e+00, -3.64059e-01, 1.70182e-03, 0.00000e+00,
		0.00000e+00, -3.20120e+00, 0.00000e+00, 5.80206e-03, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		};
		// ----   he density
		static double pd[10][151] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			},
			{ 0,
			1.04934e+00, -2.88362e-02, -2.07095e-01, -1.03314e-01, -7.02373e-03,
			1.29664e-02, 4.08853e-01, -9.19895e-03, -1.88660e-02, 1.40927e+00,
			1.75033e-01, 1.87351e-02, 1.10979e-01, -7.42871e+00, 0.00000e+00,
			2.67143e-01, -5.95979e-02, 1.05038e+02, -8.40963e-02, -6.97632e-04,
			2.06521e-06, 7.65306e-04, 0.00000e+00, 0.00000e+00, 1.26762e-01,
			1.28876e-01, -5.04479e-02, -1.30735e-02, -2.24348e-02, 0.00000e+00,
			0.00000e+00, -1.50832e+02, -6.29928e-03, 0.00000e+00, -4.07760e-03,
			0.00000e+00, 0.00000e+00, 5.25725e-02, -3.11486e+01, -3.13351e-03,
			2.75838e-03, 0.00000e+00, 0.00000e+00, 1.11247e-01, 1.08815e-01,
			-4.66713e-02, 0.00000e+00, -3.29329e-03, 0.00000e+00, 1.67838e-03,
			-9.16691e-03, 3.45044e-05, -9.71806e-03, 0.00000e+00, -2.04672e-03,
			-7.86899e-03, -7.98285e-03, 5.36515e-03, -5.31172e+03, 0.00000e+00,
			-6.42781e-03, -1.71690e-03, 0.00000e+00, -6.79131e+01, -1.79912e-02,
			-1.58305e-02, -7.12313e-03, 0.00000e+00, 2.53477e-02, 8.52960e-02,
			1.02163e-01, 2.95009e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-6.84625e+03, -6.19098e-03, -2.69289e-03, 0.00000e+00, -5.20231e+02,
			-6.33463e-03, 0.00000e+00, 0.00000e+00, -6.02428e-03, -4.07077e-03,
			5.42264e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			4.07560e-02, 2.82288e-02, 9.08088e-03, 0.00000e+00, 0.00000e+00,
			-4.05204e-01, -5.97931e-02, -7.31823e+01, -2.06620e-03, 0.00000e+00,
			-3.72723e-03, -1.88146e-02, -1.01794e-02, 8.04633e-03, 1.01090e-02,
			8.73253e-03, 2.38268e-02, 4.80444e-03, 1.71088e-03, 3.96369e-02,
			-2.13809e-02, 0.00000e+00, -1.02588e-01, -5.91702e-03, 0.00000e+00,
			2.70923e-03, 0.00000e+00, 0.00000e+00, -1.75043e+02, 6.03489e-01,
			-6.17589e-01, 8.38098e-03, 1.83871e-03, -7.05329e-04, -4.06644e+00,
			-5.09347e-03, -2.84344e-02, -1.24160e-02, 1.33665e-02, 3.93410e-03,
			-5.03723e-04, -4.57683e+00, -5.29542e-01, -4.25812e-03, 0.00000e+00,
			0.00000e+00, 1.91541e+01, 0.00000e+00, 3.23247e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   o density
			{ 0,
			9.31113e-01, -1.38721e-01, -1.33457e-01, -5.29542e-02, -4.44983e-03,
			1.35264e-02, 5.98075e-02, -3.62880e-02, -3.12798e-02, 3.72068e-01,
			2.95974e-02, 1.20509e-02, 5.21995e-02, -7.78888e+00, 0.00000e+00,
			1.18634e-01, -2.04495e-02, 1.03280e+02, 9.82432e-02, 4.77694e-04,
			0.00000e+00, 2.74372e-03, 0.00000e+00, 0.00000e+00, 7.57809e-02,
			1.71403e-01, -1.05205e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -8.73348e+00, -5.81094e-03, 0.00000e+00, -8.14944e-03,
			0.00000e+00, 0.00000e+00, 5.17255e-02, -1.53028e+01, -3.48932e-03,
			9.61771e-04, 5.57732e-03, -4.54180e-04, 9.88213e-02, 9.40456e-02,
			-3.18797e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.32122e-03,
			-6.00220e-03, 2.77654e-05, -3.22019e-03, 0.00000e+00, -3.78551e-03,
			-3.34809e-03, -1.70668e-03, 0.00000e+00, 6.36184e+03, 0.00000e+00,
			1.59986e-03, -3.88204e-03, -1.64825e-03, -7.47955e+01, -1.05360e-02,
			-9.45723e-03, -1.59824e-03, -7.06730e-04, -1.68513e-02, -1.13023e-01,
			-6.36637e-02, -1.37709e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-1.52368e+04, -5.86061e-03, -2.53108e-03, 0.00000e+00, -2.54837e+03,
			-3.28988e-03, 0.00000e+00, 0.00000e+00, -2.76364e-03, 9.67923e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			4.34255e-02, 1.14020e-02, -6.18447e-03, 0.00000e+00, 0.00000e+00,
			-3.02568e-01, -3.27694e-02, -6.71589e+01, -2.28340e-03, 0.00000e+00,
			3.06230e-03, -4.65113e-03, -9.73421e-03, 1.28326e-02, 7.88553e-03,
			7.97197e-03, -1.20760e-02, -7.67547e-03, -1.20755e-03, -2.98523e-02,
			-1.26560e-02, 0.00000e+00, -5.68350e-02, -1.53039e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.42911e-03, -4.01347e-03, -2.19074e-03, 3.11281e+00,
			3.23251e-03, -6.39523e-03, -6.63069e-03, -3.04403e-04, -4.01920e-03,
			-1.18708e-03, 4.15211e+00, -2.01896e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   n2 density
			{ 0,
			1.06903e+00, 0.00000e+00, 0.00000e+00, 3.66210e-03, 0.00000e+00,
			1.90412e-02, -1.78929e-03, 0.00000e+00, -3.92257e-02, -1.19444e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -8.45398e+00, 0.00000e+00,
			2.08180e-02, 0.00000e+00, 1.39638e+02, 8.98481e-02, 0.00000e+00,
			0.00000e+00, 3.77113e-04, 0.00000e+00, 0.00000e+00, 1.32397e-01,
			2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -2.36325e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.43022e-03,
			-3.99776e-06, 6.32343e-03, 5.48144e-03, 1.13139e-01, 1.69134e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.41470e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   msis00r.gts3c.tlb
			{ 0,
			9.76619e-01, 0.00000e+00, 0.00000e+00, -2.00200e-02, 0.00000e+00,
			-9.38391e-03, -1.95833e-03, 0.00000e+00, 1.31480e-02, -1.92414e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -8.45398e+00, 0.00000e+00,
			1.07674e-02, 0.00000e+00, 8.93820e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 5.68478e-04, 0.00000e+00, 0.00000e+00, 1.32397e-01,
			2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 4.66814e-03, 0.00000e+00, 0.00000e+00,
			5.11651e-05, 2.55717e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -2.60147e-03, -8.08556e-04, 1.13139e-01, 1.69134e-01,
			6.64196e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			5.82026e-03, 2.41470e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 6.21998e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   o2 density
			{ 0,
			9.31402e-01, 1.37976e-01, 0.00000e+00, 3.23736e-04, 0.00000e+00,
			-9.10906e-03, 7.07506e-02, 0.00000e+00, -5.16650e-02, 6.89755e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -8.45398e+00, 0.00000e+00,
			2.81140e-02, 0.00000e+00, 7.36009e+01, 5.96604e-02, 0.00000e+00,
			0.00000e+00, -1.51792e-03, 0.00000e+00, 0.00000e+00, 1.32397e-01,
			2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 9.48758e+00, 8.84541e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.13139e-01, 1.69134e-01,
			1.45192e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.07906e-02, 2.99942e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.48930e-02,
			-7.87184e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-6.83420e-02, -4.41778e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.29730e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   ar density
			{ 0,
			8.68053e-01, 2.36364e-01, 1.34306e-01, 1.03086e-02, 0.00000e+00,
			-3.79164e-03, -1.57806e-01, 0.00000e+00, -5.87644e-02, -3.12508e-01,
			0.00000e+00, 4.37387e-02, -3.54091e-02, -2.23636e+01, 0.00000e+00,
			-5.33976e-02, 0.00000e+00, 1.14091e+02, 5.17497e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.32397e-01,
			2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 3.42702e+02, 1.57033e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.66278e-03,
			-1.16193e-03, 0.00000e+00, 0.00000e+00, 1.13139e-01, 1.69134e-01,
			1.78431e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.62864e-02, 3.16963e-05, 1.27968e-02, 0.00000e+00, 0.00000e+00,
			-7.04599e-03, 2.07921e-03, 6.36660e-03, 2.29940e+04, 0.00000e+00,
			1.27833e-02, -2.08036e-03, -4.61820e-03, -6.29391e+01, -1.20745e-02,
			1.36675e-02, 1.36011e-02, -5.37162e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.92509e+04, 8.35522e-03, 4.19439e-03, 0.00000e+00, 1.20366e+04,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -1.00034e-02, -2.33267e-03,
			9.72374e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-2.65079e-02, -2.09125e-02, -1.09465e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.17252e-02, -7.12385e+01, -1.89428e-03, 0.00000e+00,
			-6.02006e-03, 1.69058e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.90646e-02,
			3.48971e-03, 0.00000e+00, 5.01174e-02, 5.50595e-02, 0.00000e+00,
			-9.55897e-03, 0.00000e+00, 0.00000e+00, -1.51693e+03, 0.00000e+00,
			0.00000e+00, 1.29306e-02, 2.69567e-03, 0.00000e+00, 3.92243e+00,
			-8.47690e-03, 1.16896e-02, 0.00000e+00, 1.48967e-02, 5.44521e-03,
			0.00000e+00, 5.64918e+00, 0.00000e+00, -7.72178e-03, 0.00000e+00,
			0.00000e+00, -7.34042e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    h density
			{ 0,
			1.27515e+00, -2.10472e-01, -1.77924e-01, 2.18900e-01, 2.88436e-02,
			1.90077e-02, 2.91001e-01, 2.17437e-02, -1.05186e-02, 4.36141e-01,
			1.07605e-01, 3.30755e-02, 4.00581e-02, -9.58051e+00, 0.00000e+00,
			1.54028e-02, 0.00000e+00, 7.34194e+01, 4.96540e-02, -5.95906e-03,
			3.84512e-05, -1.36000e-02, 0.00000e+00, 0.00000e+00, 1.32397e-01,
			2.13315e-01, -4.16610e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.46276e+02, -1.98408e-02, 0.00000e+00, 1.32530e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.04687e-04,
			-1.47562e-03, 0.00000e+00, 0.00000e+00, 1.13139e-01, 1.69134e-01,
			-1.26913e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, -6.08370e-03,
			-2.57587e-02, 3.19022e-05, 0.00000e+00, 0.00000e+00, 1.56644e-02,
			1.03640e-02, 1.05771e-03, 0.00000e+00, 3.57949e+03, 0.00000e+00,
			-1.25672e-03, 1.52783e-03, 1.30518e-03, 7.55558e+00, -9.20341e-03,
			-2.09142e-02, -1.34106e-02, 0.00000e+00, -4.83312e-02, 8.30900e-02,
			9.88009e-02, -1.41148e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-1.05513e+03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 6.73442e-03, 2.01691e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			5.98019e-02, 6.33298e-03, -1.12871e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.28604e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-4.94960e-03, -1.36415e-02, -1.15039e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -5.86860e-03, -1.41732e-03, 2.13697e-03, 2.63845e+00,
			-8.34186e-03, -1.87336e-02, -1.90870e-02, -8.03810e-03, -2.84279e-03,
			2.56722e-03, 1.71429e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    n density
			{ 0,
			5.73587e+01, -3.98747e-01, 0.00000e+00, -5.29554e-01, -5.82186e-03,
			7.14177e-02, -6.79279e-01, -1.67715e-01, -6.42434e-02, -2.11569e-01,
			-1.59922e-01, -1.71024e-04, -1.15885e-01, 6.51603e+00, 0.00000e+00,
			-1.76683e-01, 6.50395e-02, 1.43504e+00, 9.28208e-02, 5.11662e-03,
			0.00000e+00, 9.95121e-03, 0.00000e+00, 0.00000e+00, 1.32397e-01,
			2.13315e-01, 1.01451e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 5.67667e+01, 2.38192e-03, 0.00000e+00, -1.88240e-02,
			0.00000e+00, 0.00000e+00, 4.76218e-02, 2.35206e+01, 4.75901e-03,
			5.76162e-03, 1.51815e-02, -1.92730e-02, 1.13139e-01, 1.69134e-01,
			-2.88771e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.18418e-03,
			-3.68927e-03, 3.14704e-05, 8.82198e-03, 0.00000e+00, -1.92562e-02,
			-2.58674e-03, -2.19913e-02, 0.00000e+00, 4.38655e+03, 0.00000e+00,
			7.60126e-03, 2.59438e-03, 1.72310e-03, 7.79204e+01, 7.97786e-04,
			-7.70510e-03, 1.90982e-03, 2.72707e-03, 1.01016e-02, 1.16537e-01,
			-3.12236e-03, 1.39783e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-1.30712e+03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -3.20544e-03, -2.06970e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.59010e-02, -1.91427e-03, -3.42829e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -3.45379e-02, 8.94518e+01, 1.71556e-03, 0.00000e+00,
			-7.65278e-03, -2.08987e-04, -1.57393e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -8.60673e-03, -1.19922e-02, -6.46356e-03, -3.00107e+00,
			-9.32511e-03, -1.50205e-02, -8.67835e-03, -7.64801e-03, -1.31495e-02,
			-6.76720e-03, -1.82396e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----  spare
			{ 0,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -8.45398e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.32397e-01,
			2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.13139e-01, 1.69134e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			}
		};
		// ----    msis90r.gts3c.ps param
		static double ps[151] =
		{ 0,
		9.51363e-01, -4.67542e-02, 1.20260e-01, 0.00000e+00, 0.00000e+00,
		1.91357e-02, 0.00000e+00, 0.00000e+00, 1.25429e-03, -1.33240e-01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -8.45398e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 2.52317e-03, 0.00000e+00, -9.73404e-03, 1.32397e-01,
		2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -7.18482e-04, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 7.87683e-03, -2.33698e-03, 1.13139e-01, 1.69134e-01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		};
		// ----    turbo
		static double pdl[3][26] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0
			},
			{ 0,
			9.33804e-01, 5.47446e+00, 1.53263e-01, 9.19303e-01, 1.64109e+01,
			4.27083e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.40925e-01,
			},
			{ 0,
			1.15897e+00, 4.71094e-01, 1.09459e+00, 5.25012e+00, 1.00000e+00,
			1.00000e+00, 1.03999e+00, 7.67132e-01, 1.10514e+00, 1.75636e+00,
			1.10845e+00, 2.33439e+00, 7.96532e-01, 4.31520e+00, 4.07300e+00,
			1.22807e+02, 2.39547e-01, 2.53791e-06, 8.42931e-01, 1.04192e+00,
			2.00202e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 9.62736e-01
			}
		};
		// ----   lower bouneary
		static double ptm[11] =
		{ 0,
		1.04130e+03, 3.86000e+02, 1.95000e+02, 1.66728e+01, 2.13000e+02,
		1.20000e+02, 2.40000e+02, 1.87000e+02, -2.00000e+00, 0.00000e+00
		};
		static double pdm[9][11] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			},
			{ 0,
			2.45600e+07, 6.71072e-06, 1.00000e+02, 0.00000e+00, 1.10000e+02,
			1.00000e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			8.59400e+10, 5.40000e-01, 1.05000e+02, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			2.81000e+11, 0.00000e+00, 1.05000e+02, 2.80000e+01, 2.89500e+01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			3.30000e+10, 2.68270e-01, 1.05000e+02, 0.00000e+00, 1.10000e+02,
			1.00000e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			1.33000e+09, 1.19615e-02, 1.05000e+02, 0.00000e+00, 1.10000e+02,
			1.00000e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			1.76100e+05, 1.00000e+00, 9.50000e+01, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			1.00000e+07, 1.00000e+00, 1.05000e+02, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			1.00000e+07, 1.00000e+00, 1.05000e+02, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
			}
		};
		// ----   msis90r.meso.tn1[2]
		static double ptl[5][101] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			},
			{ 0,
			1.02083e+00, 4.08449e-02, -2.34582e-02, 4.38274e-04, -1.52380e-02,
			-2.09089e-02, 4.46355e-03, -3.41250e-03, -1.12961e-02, -7.03277e-02,
			-4.82724e-02, 0.00000e+00, 0.00000e+00, -6.20496e+00, 0.00000e+00,
			-9.80197e-03, -1.45065e-02, -1.13226e+02, 2.28455e-02, 0.00000e+00,
			0.00000e+00, 4.93658e-04, 0.00000e+00, 3.79078e-03, 1.32397e-01,
			2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -8.89051e+03, 2.25900e-03, 1.76142e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.55015e-04,
			2.21388e-03, -5.99073e-04, -3.52331e-03, 1.13139e-01, 1.69134e-01,
			7.79156e-03, -1.93458e-03, -1.08596e-02, -4.39285e-04, 0.00000e+00,
			3.83994e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 6.76608e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   msis90r.meso.tn1[3]
			{ 0,
			9.24880e-01, 7.41986e-02, -6.37629e-03, 6.00575e-03, 1.29382e-03,
			6.97550e-03, -1.70782e-03, 2.80584e-03, -8.87214e-03, -4.35703e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 4.31515e+00, 0.00000e+00,
			-1.81474e-02, -6.06627e-02, -8.43503e+01, 8.46944e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -2.17081e-02, -2.19500e-03, 1.32397e-01,
			2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.47580e+02, 4.41585e-03, 7.80466e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 6.44155e-04,
			-2.49166e-03, 2.90482e-03, -3.40501e-04, 1.13139e-01, 1.69134e-01,
			-6.01460e-03, -1.63368e-03, 0.00000e+00, -4.31340e-03, 0.00000e+00,
			4.53979e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -5.43660e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   msis90r.meso.tn1[4]
			{ 0,
			9.72669e-01, -4.26748e-02, 1.12876e-02, -8.44951e-03, 7.04114e-03,
			1.26036e-02, -3.88164e-03, -5.20509e-04, -6.09710e-04, 1.31603e-01,
			1.13804e-01, 0.00000e+00, 0.00000e+00, -6.15970e+00, 0.00000e+00,
			-2.14214e-02, -6.62913e-02, -2.02884e-01, 2.35350e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.13573e-02, -1.84905e-03, 1.32397e-01,
			2.13315e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.42645e+00, -2.64405e-03, -5.57771e-04, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -2.20621e+01, -1.10313e-03,
			3.97063e-05, 5.47632e-05, 3.57577e-03, 1.13139e-01, 1.69134e-01,
			0.00000e+00, 1.18897e-03, 0.00000e+00, 7.62305e-04, 0.00000e+00,
			-3.52015e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -9.52550e-04,
			8.56253e-04, 4.33114e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.21223e-03,
			2.38694e-04, 9.15245e-04, 1.28385e-03, 8.67668e-04, -5.61425e-06,
			1.04445e+00, 3.41112e+01, 0.00000e+00, -8.40704e-01, -2.39639e+02,
			7.06668e-01, -2.05873e+01, -3.63696e-01, 2.39245e+01, 1.00000e+01,
			-1.06657e-03, -7.67292e-04, 1.54534e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   msis90r.meso.tn1[5] msis90r.meso.tn2[1]
			{ 0,
			9.99368e-01, 4.33893e-02, -2.07009e-03, 1.09617e-03, 1.05440e-03,
			4.83408e-04, 9.77040e-04, 9.24791e-04, 4.80247e-04, 4.94737e-02,
			1.05985e-03, 0.00000e+00, 0.00000e+00, 2.74409e+00, 0.00000e+00,
			-4.96656e-03, -1.51684e-02, 4.65158e+01, -7.51133e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 6.63808e-04, 1.32397e-01,
			2.13315e-01, -2.06652e-03, -6.32046e-03, 0.00000e+00, 0.00000e+00,
			5.94545e-03, -1.90958e+02, 0.00000e+00, -4.16892e-03, 0.00000e+00,
			-1.67499e-02, 0.00000e+00, 2.58987e-03, 5.97781e+02, 0.00000e+00,
			0.00000e+00, 4.44890e-04, 4.66444e-04, 1.13139e-01, 1.69134e-01,
			0.00000e+00, 7.11360e-04, 1.32186e-02, 2.23948e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.60571e-03,
			6.28078e-04, 5.05469e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.57829e-03,
			-4.00855e-04, 5.04077e-05, -1.39001e-03, -2.33406e-03, -4.81197e-04,
			1.46758e+00, 6.20332e+00, 0.00000e+00, 3.66476e-01, -6.19760e+01,
			3.09198e-01, -1.98999e+01, 0.00000e+00, -3.29933e+02, 0.00000e+00,
			-1.10080e-03, -9.39310e-05, 1.39638e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			}
		};
		// ----    msis90r.meso.tn2[2]
		static double pma[11][101] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			},
			{ 0,
			9.81637e-01, -1.41317e-03, 3.87323e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.58707e-02,
			-8.63658e-03, 0.00000e+00, 0.00000e+00, -2.02226e+00, 0.00000e+00,
			-8.69424e-03, -1.91397e-02, 8.76779e+01, 4.52188e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -7.07572e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-4.11210e-03, 3.50060e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			2.23760e-02, 0.00000e+00, -8.36657e-03, 1.61347e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.45130e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.24152e-03,
			6.43365e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.33255e-03,
			2.42657e-03, 1.60666e-03, -1.85728e-03, -1.46874e-03, -4.79163e-06,
			1.22464e+00, 3.53510e+01, 0.00000e+00, 4.49223e-01, -4.77466e+01,
			4.70681e-01, 8.41861e+00, -2.88198e-01, 1.67854e+02, 0.00000e+00,
			7.11493e-04, 6.05601e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    msis90r.meso.tn2[3]
			{ 0,
			1.00422e+00, -7.11212e-03, 5.24480e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.28914e-02,
			-2.41301e-02, 0.00000e+00, 0.00000e+00, -2.12219e+01, 0.00000e+00,
			-3.28077e-03, 1.65727e-02, 1.68564e+00, -6.68154e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 8.42365e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -4.34645e-03, -1.03830e-02, -8.08279e-03, 2.16780e-02,
			0.00000e+00, -1.38459e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.45155e-02, 0.00000e+00, 7.04573e-03, -4.73204e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.08767e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.21769e-04,
			-2.27387e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.26769e-03,
			3.16901e-03, 4.60316e-04, -1.01431e-04, 1.02131e-03, 9.96601e-04,
			1.25707e+00, 2.50114e+01, 0.00000e+00, 4.24472e-01, -2.77655e+01,
			3.44625e-01, 2.75412e+01, 0.00000e+00, 7.94251e+02, 0.00000e+00,
			2.45835e-03, 1.38871e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    msis90r.meso.tn2[4] tn3[1]
			{ 0,
			1.01890e+00, -2.46603e-02, 1.00078e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -6.70977e-02,
			-4.02286e-02, 0.00000e+00, 0.00000e+00, -2.29466e+01, 0.00000e+00,
			2.26580e-03, 2.63931e-02, 3.72625e+01, -6.39041e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.85291e-03, -7.47019e-03, -7.07265e-03, 0.00000e+00,
			0.00000e+00, 1.39717e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			9.58383e-03, 0.00000e+00, 9.19771e-03, -3.69121e+02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.57067e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.92953e-03,
			-2.77739e-03, -4.40092e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.47280e-03,
			2.95035e-04, -1.81246e-03, 2.81945e-03, 4.27296e-03, 9.78863e-04,
			1.40545e+00, -6.19173e+00, 0.00000e+00, 0.00000e+00, -7.93632e+01,
			4.44643e-01, -4.03085e+02, 0.00000e+00, 1.15603e+01, 0.00000e+00,
			2.25068e-03, 8.48557e-04, -2.98493e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    msis90r.meso.tn3[2]
			{ 0,
			9.75801e-01, 3.80680e-02, -3.05198e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.85575e-02,
			5.04057e-02, 0.00000e+00, 0.00000e+00, -1.76046e+02, 0.00000e+00,
			-1.48297e-03, -3.68560e-03, 3.02185e+01, -3.23338e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.15558e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 4.89620e-03, 1.44594e-02, 9.91215e-03, -1.00616e-02,
			-8.21324e-03, -1.57757e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.53569e-02, 0.00000e+00, 6.63564e-03, 4.58410e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -2.51280e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -8.73148e-04,
			-1.29648e-03, -7.32026e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.68110e-03,
			-4.66003e-03, -1.31567e-03, -7.39390e-04, 6.32499e-04, -4.65588e-04,
			-1.29785e+00, -1.57139e+02, 0.00000e+00, 2.58350e-01, -3.69453e+01,
			4.10672e-01, 9.78196e+00, -1.52064e-01, -3.85084e+03, 0.00000e+00,
			-8.52706e-04, -1.40945e-03, -7.26786e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    msis90r.meso.tn3[3]
			{ 0,
			9.60722e-01, 7.03757e-02, -3.00266e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.22671e-02,
			4.10423e-02, 0.00000e+00, 0.00000e+00, -1.63070e+02, 0.00000e+00,
			5.40747e-04, 7.79481e-03, 1.44908e+02, 1.51484e-04, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.41844e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 5.77884e-03, 1.06073e-02, 5.36685e-03, 9.74319e-03,
			0.00000e+00, -2.88015e+03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.97547e-02, 0.00000e+00, -4.44902e-03, -2.92760e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 2.34419e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.65325e-04,
			-5.50628e-04, 3.31465e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.06179e-03,
			-3.08575e-03, -7.93589e-04, -1.08629e-04, 5.95511e-04, -9.05050e-04,
			1.18997e+00, 4.15924e+01, 0.00000e+00, -4.72064e-01, -9.47150e+02,
			3.98723e-01, 1.98304e+01, 0.00000e+00, 3.73219e+03, 0.00000e+00,
			-1.50040e-03, -1.14933e-03, -1.56769e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    msis90r.meso.tn3[4]
			{ 0,
			1.03123e+00, -7.05124e-02, 8.71615e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.82621e-02,
			-9.80975e-03, 0.00000e+00, 0.00000e+00, 2.89286e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 8.66153e+01, 7.91938e-04, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 4.68917e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 7.86638e-03, 9.57341e-03, 5.72268e-03, 9.90827e-03,
			0.00000e+00, 6.55573e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -4.00200e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 7.07457e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.04970e-04,
			1.21560e-03, -8.05579e-06, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.49941e-03,
			-4.57256e-04, -1.59311e-04, 2.96481e-04, -1.77318e-03, -6.37918e-04,
			1.02395e+00, 1.28172e+01, 0.00000e+00, 1.49903e-01, -2.63818e+01,
			0.00000e+00, 4.70628e+01, -2.22139e-01, 4.82292e-02, 0.00000e+00,
			-8.67075e-04, -5.86479e-04, 5.32462e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    tn3[5] surface temp tsl
			{ 0,
			1.00828e+00, -9.10404e-02, -2.26549e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.32420e-02,
			-9.08925e-03, 0.00000e+00, 0.00000e+00, 3.36105e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.24957e+01, -5.87939e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.79765e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 2.01237e+03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.75553e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.29699e-03,
			1.26659e-03, 2.68402e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.17894e-03,
			1.48746e-03, 1.06478e-04, 1.34743e-04, -2.20939e-03, -6.23523e-04,
			6.36539e-01, 1.13621e+01, 0.00000e+00, -3.93777e-01, 2.38687e+03,
			0.00000e+00, 6.61865e+02, -1.21434e-01, 9.27608e+00, 0.00000e+00,
			1.68478e-04, 1.24892e-03, 1.71345e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    tgn3[2] surface grae tslg
			{ 0,
			1.57293e+00, -6.78400e-01, 6.47500e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -7.62974e-02,
			-3.60423e-01, 0.00000e+00, 0.00000e+00, 1.28358e+02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 4.68038e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.67898e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.90994e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 3.15706e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    tgn2[1] tgn1[2]
			{ 0,
			8.66492e-01, 3.55807e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.12111e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.82458e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.01024e+02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 6.54251e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.56959e-02,
			1.91001e-02, 3.15971e-02, 1.00982e-02, -6.71565e-03, 2.57693e-03,
			1.38692e+00, 2.82132e-01, 0.00000e+00, 0.00000e+00, 3.81511e+02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    tgn3[1] tgn2[2]
			{ 0,
			1.06029e+00, -5.25231e-02, 3.73034e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.31072e-02,
			-3.88409e-01, 0.00000e+00, 0.00000e+00, -1.65295e+02, 0.00000e+00,
			-4.38916e-02, -3.22716e-01, -8.82393e+01, 1.18458e-01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.19782e-01, -2.13801e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.62229e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-4.35863e-01, 0.00000e+00, 0.00000e+00, -5.37443e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -4.55788e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.84009e-02,
			3.96733e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.05494e-02,
			7.39617e-02, 1.92200e-02, -8.46151e-03, -1.34244e-02, 1.96338e-02,
			1.50421e+00, 1.88368e+01, 0.00000e+00, 0.00000e+00, -5.13114e+01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			5.11923e-02, 3.61225e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			}
		};
		// ----   middle atmosphere averages
		static double pavgm[11] =
		{ 0,
		2.61000e+02, 2.64000e+02, 2.29000e+02, 2.17000e+02, 2.17000e+02,
		2.23000e+02, 2.86760e+02, -2.93940e+00, 2.50000e+00, 0.00000e+00
		};

		// now assign the record values
		memcpy(msis90r.parm.pt, pt, 151 * sizeof(double));
		memcpy(msis90r.parm.ps, ps, 151 * sizeof(double));
		memcpy(msis90r.lower.ptm, ptm, 11 * sizeof(double));
		memcpy(msis90r.mavg.pavgm, pavgm, 11 * sizeof(double));

		// since c reads in rows first, while the fortran code assumed columns,
		// we do a loop

		for (i = 0; i <= 8; i++)
		for (j = 0; j <= 10; j++)
			msis90r.lower.pdm[j][i] = pdm[i][j];
		for (i = 0; i <= 4; i++)
		for (j = 0; j <= 100; j++)
			msis90r.parm.ptl[j][i] = ptl[i][j];
		for (i = 0; i <= 9; i++)
		for (j = 0; j <= 150; j++)
			msis90r.parm.pd[j][i] = pd[i][j];

		for (i = 0; i <= 2; i++)
		for (j = 0; j <= 25; j++)
			msis90r.parm.pdl[j][i] = pdl[i][j];
		for (i = 0; i <= 10; i++)
		for (j = 0; j <= 100; j++)
			msis90r.parm.pma[j][i] = pma[i][j];

		// cdav
		// these are the extra arrays needed so a whole array can be accessed, but in
		// reverse order from the standard [row,col] approach in the fortran.
		//    extras to aid function calls for a subarry
		memcpy(msis90r.parm.pddav, pd, 10 * 151 * sizeof(double));
		memcpy(msis90r.parm.ptldav, ptl, 5 * 101 * sizeof(double));
		memcpy(msis90r.parm.pmadav, pma, 11 * 101 * sizeof(double));

	}

	
	// ---------------- NRLMSISE-00 model ----------------
	/*    ******************* local routines ******************************       */
	void gts7
		(
		msistype& msis00r,
		lpolytype& lpoly,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
		double f107, double ap[8], int mass, double d[10], double t[3]
		);

	double scalh
		(
		parmbtype& parmb,
		double alt, double xm, double temp
		);

	double g0
		(
		double a, double p[151]
		);

	double sumex
		(
		double ex
		);
	double sg0
		(
		double ex, double p[151], double ap[8]
		);

	double globe7
		(
		cswtype& csw,
		lpolytype& lpoly,
		double yrd, double sec, double lat, double llong, double tloc, double f107a,
		double f107, double ap[8], double p[151]
		);

	double glob7s
		(
		lpolytype& lpoly,
		cswtype& csw,
		double p[101]
		);



	/*  ----------------------------------------------------------------------
	*
	*     nrlmsise-00
	*     -----------
	*          neutral atmosphere empirical model from the surface to lower
	*          exosphere
	*
	*          new features:
	*            *extensive satellite drag database used in model generation
	*            *revised o2 (and o) in lower thermosphere
	*            *additional nonlinear solar activity term
	*            *"anomalous oxygen" number density, output d[9]
	*             at high altitudes (> 500 km), hot atomic oxygen or ionized
	*             oxygen can become appreciable for some ranges of subroutine
	*             inputs, thereby affecting drag on satellites and debris. we
	*             group these species under the term "anomalous oxygen," since
	*             their individual variations are not presently separable with
	*             the drag data used to define this model component.
	*
	*          subroutines for special outputs:
	*
	*          high altitude drag: effective total mass density
	*          (subroutine gtd7d, output d[6])
	*             for atmospheric drag calculations at altitudes above 500 km,
	*             call subroutine gtd7d to compute the "effective total mass
	*             density" by including contributions from "anomalous oxygen."
	*             see "notes on output variables" below on d[6].
	*
	*          pressure grid (subroutine ghp7)
	*            see subroutine ghp7 to specify outputs at a pressure level
	*            rather than at an altitude.
	*
	*          output in m-3 and kg/m3:   call meters(.true.)
	*
	*     input variables:
	*          iyd - year and day as yyddd (day of year from 1 to 365 (or 366))
	*                (year ignored in current model)
	*          sec - ut(sec)
	*          alt - altitude(km)
	*          glat - geodetic latitude(deg)
	*          glong - geodetic longitude(deg)
	*          stl - local apparent solar time(hrs; see note below)
	*          f107a - 81 day average of f10.7 flux (centered on day ddd)
	*          f107 - daily f10.7 flux for previous day
	*          ap - magnetic index(daily) or when msis00r.csw.sw[9]=-1. :
	*             - array containing:
	*               [1] daily ap
	*               [2] 3 hr ap index for current time
	*               [3] 3 hr ap index for 3 hrs before current time
	*               [4] 3 hr ap index for 6 hrs before current time
	*               [5] 3 hr ap index for 9 hrs before current time
	*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
	*                      to current time
	*               [7] average of eight 3 hr ap indicies from 36 to 57 hrs prior
	*                      to current time
	*          mass - mass number (only density for selected gas is
	*                   calculated.  mass 0 is temperature.  mass 48 for all.
	*                   mass 17 is anomalous o only.)
	*
	*     notes on input variables:
	*          ut, local time, and longitude are used independently in the
	*          model and are not of equal importance for every situation.
	*          for the most physically realistic calculation these three
	*          variables should be consistent (stl=sec/3600+glong/15).
	*          the equation of time departures from the above formula
	*          for apparent local time can be included if available but
	*          are of minor importance.
	*
	*          f107 and f107a values used to generate the model correspond
	*          to the 10.7 cm radio flux at the actual distance of the earth
	*          from the sun rather than the radio flux at 1 au. the following
	*          site provides both classes of values:
	*          ftp://ftp.ngdc.noaa.gov/stp/solar_data/solar_radio/flux/
	*
	*          f107, f107a, and ap effects are neither large nor well
	*          established below 80 km and these parameters should be set to
	*          150., 150., and 4. respectively.
	*
	*     output variables:
	*          d[1] - he number density(cm-3)
	*          d[2] - o number density(cm-3)
	*          d[3] - n2 number density(cm-3)
	*          d[4] - o2 number density(cm-3)
	*          d[5] - ar number density(cm-3]
	*          d[6] - total mass density(gm/cm3]
	*          d[7] - h number density(cm-3]
	*          d[8] - n number density(cm-3]
	*          d[9] - anomalous oxygen number density(cm-3]
	*          t[1] - exospheric temperature
	*          t[2] - temperature at alt
	*
	*     notes on output variables:
	*          to get output in m-3 and kg/m3:   call meters(.true.)
	*
	*          o, h, and n are set to zero below 72.5 km
	*
	*          t[1], exospheric temperature, is set to global average for
	*          altitudes below 120 km. the 120 km gradient is left at global
	*          average value for altitudes below 72 km.
	*
	*          d[6], total mass density, is not the same for subroutines gtd7
	*          and gtd7d
	*
	*            subroutine gtd7 -- d[6] is the sum of the mass densities of the
	*            species labeled by indices 1-5 and 7-8 in output variable d.
	*            this includes he, o, n2, o2, ar, h, and n but does not include
	*            anomalous oxygen (species index 9).
	*
	*            subroutine gtd7d -- d[6] is the "effective total mass density
	*            for drag" and is the sum of the mass densities of all species
	*            in this model, including anomalous oxygen.
	*
	*     msis00r.csw.switches: the following is for test and special purposes:
	*
	*          to turn on and off particular variations call tselec(msis00r.csw.sw),
	*          where msis00r.csw.sw is a 25 element array containing 0. for off, 1.
	*          for on, or 2. for main effects off but cross terms on
	*          for the following variations
	*                 1 - f10.7 effect on mean  2 - time independent
	*                 3 - symmetrical annual    4 - symmetrical semiannual
	*                 5 - asymmetrical annual   6 - asymmetrical semiannual
	*                 7 - diurnal               8 - semidiurnal
	*                 9 - daily ap             10 - all ut/long effects
	*                11 - longitudinal         12 - ut and mixed ut/long
	*                13 - mixed ap/ut/long     14 - terdiurnal
	*                15 - departures from diffusive equilibrium
	*                16 - all tinf var         17 - all tlb var
	*                18 - all tn1 var           19 - all s var
	*                20 - all tn2 var           21 - all nlb var
	*                22 - all tn3 var           23 - turbo scale height var
	*
	*          to get current values of msis00r.csw.sw: call tretrv(msis00r.csw.sw)
	*
	* --------------------------------------------------------------------------- */

	void gtd7
		(
		msistype& msis00r,
		lpolytype& lpoly,
		fittype& fit,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong, double stl,
		double f107a, double f107, double ap[8], int mass,
		double d[10], double t[3]
		)
	{
		double ds[10], ts[3], v1, xlat, xmm, altt, dm28m, dmc, dz28, dmr, tz;
		int i, mss;

		// make sure and initilize the data through msis00init at the start

		int mn3 = 5;
		double zn3[6] = { 0, 32.5, 20.0, 15.0, 10.0, 0.0 };
		int mn2 = 4;
		double zn2[5] = { 0, 72.5, 55.0, 45.0, 32.5 };
		double zmix = 62.5;

		double mssl = -999;
		int sv[26] = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

		// treat as local, but save in between calls
		static double alast = 99999.0;

		//------------------------------ begin ---------------------------------
		if (msis00r.csw.isw != 64999)
			tselec(msis00r.csw, sv);

		// ----  test for changed input
		v1 = vtst(msis00r.csw, iyd, sec, glat, glong, stl, f107a, f107, ap, 1);

		// ---- latitude variation of gravity (none for msis00r.csw.sw[2]=0)
		xlat = glat;
		if (msis00r.csw.sw[2] == 0)
			xlat = 45.0;
		glatf(xlat, msis00r.parmb.gsurf, msis00r.parmb.re);

		xmm = msis00r.lower.pdm[5][3];

		// ---- thermosphere/mesosphere (above zn2[1])
		if (alt >= zn2[1])
			altt = alt;
		else
			altt = zn2[1];

		mss = mass;
		// ---- only calculate n2 in thermosphere if alt in mixed region
		if (alt < zmix && mass > 0)
			mss = 28;

		// ---- only calculate thermosphere if input parameters changed
		// ----   or altitude above zn2[1] in mesosphere
		if ((v1 == 1.0) | (alt > zn2[1]) | (alast > zn2[1]) | (mss != mssl))
		{
			gts7(msis00r, lpoly, lsqv, iyd, sec, altt, glat, glong, stl,
				f107a, f107, ap, mss, ds, ts);
			dm28m = msis00r.dmix.dm28;
			// ----   metric adjustment
			if (msis00r.metsel.imr == 1)
				dm28m = msis00r.dmix.dm28*1.0e6;
			mssl = mss;
		}
		t[1] = ts[1];
		t[2] = ts[2];
		if (alt >= zn2[1])
		{
			for (i = 1; i <= 9; i++)
				d[i] = ds[i];
			goto ten;
		}

		// ---- lower mesosphere/upper stratosphere [between zn3[1] and zn2[1]]
		// ----   temperature at nodes and gradients at end nodes
		// ----   inverse temperature a linear function of spherical harmonics
		// ----   only calculate nodes if input changed
		if ((v1 == 1.0) | (alast >= zn2[1]))
		{
			msis00r.meso.tgn2[1] = msis00r.meso.tgn1[2];
			msis00r.meso.tn2[1] = msis00r.meso.tn1[5];
			msis00r.meso.tn2[2] = msis00r.parm.pma[1][1] * msis00r.mavg.pavgm[1] /
				(1.0 - msis00r.csw.sw[20] * glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[1]));
			msis00r.meso.tn2[3] = msis00r.parm.pma[1][2] * msis00r.mavg.pavgm[2] /
				(1.0 - msis00r.csw.sw[20] * glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[2]));
			msis00r.meso.tn2[4] = msis00r.parm.pma[1][3] * msis00r.mavg.pavgm[3] /
				(1.0 - msis00r.csw.sw[20] * msis00r.csw.sw[22] *
				glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[3]));
			msis00r.meso.tgn2[2] = msis00r.mavg.pavgm[9] * msis00r.parm.pma[1][10] *
				(1.0 + msis00r.csw.sw[20] * msis00r.csw.sw[22] *
				glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[10]))
				*msis00r.meso.tn2[4] * msis00r.meso.tn2[4] /
				pow((msis00r.parm.pma[1][3] * msis00r.mavg.pavgm[3]), 2);
			msis00r.meso.tn3[1] = msis00r.meso.tn2[4];
		}
		if (alt < zn3[1])
		{
			// ---- lower stratosphere and troposphere [below zn3[1]]
			// ----   temperature at nodes and gradients at end nodes
			// ----   inverse temperature a linear function of spherical harmonics
			// ----   only calculate nodes if input changed
			if ((v1 == 1.0) | (alast >= zn3[1]))
			{
				msis00r.meso.tgn3[1] = msis00r.meso.tgn2[2];
				msis00r.meso.tn3[2] = msis00r.parm.pma[1][4] * msis00r.mavg.pavgm[4]
					/ (1.0 - msis00r.csw.sw[22] * glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[4]));
				msis00r.meso.tn3[3] = msis00r.parm.pma[1][5] * msis00r.mavg.pavgm[5]
					/ (1.0 - msis00r.csw.sw[22] * glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[5]));
				msis00r.meso.tn3[4] = msis00r.parm.pma[1][6] * msis00r.mavg.pavgm[6]
					/ (1.0 - msis00r.csw.sw[22] * glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[6]));
				msis00r.meso.tn3[5] = msis00r.parm.pma[1][7] * msis00r.mavg.pavgm[7]
					/ (1.0 - msis00r.csw.sw[22] * glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[7]));
				msis00r.meso.tgn3[2] = msis00r.parm.pma[1][8] * msis00r.mavg.pavgm[8] *
					(1.0 + msis00r.csw.sw[22] * glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[8]))
					*msis00r.meso.tn3[5] * msis00r.meso.tn3[5] /
					pow((msis00r.parm.pma[1][7] * msis00r.mavg.pavgm[7]), 2);
			}
		}

		if (mass == 0)
			goto fifty;
		// ----    linear transition to full mixing below zn2[1]
		dmc = 0;
		if (alt > zmix)
			dmc = 1.0 - (zn2[1] - alt) / (zn2[1] - zmix);
		dz28 = ds[3];

		// ----***** n2 density ****
		dmr = ds[3] / dm28m - 1.0;
		d[3] = densm(msis00r.parmb, fit, lsqv, alt, dm28m, xmm, tz, mn3, zn3,
			msis00r.meso.tn3, msis00r.meso.tgn3, mn2, zn2, msis00r.meso.tn2,
			msis00r.meso.tgn2);
		d[3] = d[3] * (1.0 + dmr*dmc);

		// ----***** he density ****
		d[1] = 0.0;
		if ((mass == 4) | (mass == 48))
		{
			dmr = ds[1] / (dz28*msis00r.lower.pdm[2][1]) - 1.0;
			d[1] = d[3] * msis00r.lower.pdm[2][1] * (1.0 + dmr*dmc);
		}

		// ----**** o density ****
		d[2] = 0.0;
		d[9] = 0.0;

		// ----***** o2 density ****
		d[4] = 0.0;
		if ((mass == 32) | (mass == 48))
		{
			dmr = ds[4] / (dz28*msis00r.lower.pdm[2][4]) - 1.0;
			d[4] = d[3] * msis00r.lower.pdm[2][4] * (1.0 + dmr*dmc);
		}

		// ----***** ar density ****
		d[5] = 0.0;
		if ((mass == 40) | (mass == 48))
		{
			dmr = ds[5] / (dz28*msis00r.lower.pdm[2][5]) - 1.0;
			d[5] = d[3] * msis00r.lower.pdm[2][5] * (1.0 + dmr*dmc);
		}

		// ----***** hydrogen density ****
		d[7] = 0.0;

		// ----***** atomic nitrogen density ****
		d[8] = 0.0;

		// ---- total mass density

		if (mass == 48)
		{
			d[6] = 1.66e-24*(4.0*d[1] + 16.0*d[2] + 28.0*d[3] + 32.0*d[4] + 40.0*d[5] + d[7] + 14.0*d[8]);
			if (msis00r.metsel.imr == 1)
				d[6] = d[6] / 1000.0;
		}
		t[2] = tz;
	ten:
		goto ninety;
	fifty:
		msis00r.gts3c.dd = densm(msis00r.parmb, fit, lsqv,
			alt, 1.0, 0.0, tz, mn3, zn3, msis00r.meso.tn3, msis00r.meso.tgn3, mn2, zn2,
			msis00r.meso.tn2, msis00r.meso.tgn2);
		t[2] = tz;
	ninety:
		alast = alt;
	}

	/* -----------------------------------------------------------------------
	*
	*     nrlmsise-00
	*     -----------
	*          this void provides effective total mass density for
	*          output d[6] which includes contributions from "anomalous
	*          oxygen" which can affect satellite drag above 500 km.  this
	*          void is part of the distribution package for the
	*          neutral atmosphere empirical model from the surface to lower
	*          exosphere.  see void gtd7 for more extensive comments.
	*
	*     input variables:
	*          iyd - year and day as yyddd (day of year from 1 to 365 (or 366])
	*                (year ignored in current model)
	*          sec - ut(sec)
	*          alt - altitude(km)
	*          glat - geodetic latitude(deg)
	*          glong - geodetic longitude(deg)
	*          stl - local apparent solar time(hrs; see note below)
	*          f107a - 81 day average of f10.7 flux (centered on day ddd)
	*          f107 - daily f10.7 flux for previous day
	*          ap - magnetic index(daily) or when msis00r.csw.sw[9] = -1.0 :
	*             - array containing:
	*               [1] daily ap
	*               [2] 3 hr ap index for current time
	*               [3] 3 hr ap index for 3 hrs before current time
	*               [4] 3 hr ap index for 6 hrs before current time
	*               [5] 3 hr ap index for 9 hrs before current time
	*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
	*                      to current time
	*               [7] average of eight 3 hr ap indicies from 36 to 57 hrs prior
	*                      to current time
	*          mass - mass number (only density for selected gas is
	*                   calculated.  mass 0 is temperature.  mass 48 for all.
	*                   mass 17 is anomalous o only.)
	*
	*     notes on input variables:
	*          ut, local time, and longitude are used independently in the
	*          model and are not of equal importance for every situation.
	*          for the most physically realistic calculation these three
	*          variables should be consistent (stl = sec/3600+glong/15).
	*          the equation of time departures from the above formula
	*          for apparent local time can be included if available but
	*          are of minor importance.
	*
	*          f107 and f107a values used to generate the model correspond
	*          to the 10.7 cm radio flux at the actual distance of the earth
	*          from the sun rather than the radio flux at 1 au.
	*
	*     output variables:
	*          d[1] - he number density(cm-3]
	*          d[2] - o number density(cm-3]
	*          d[3] - n2 number density(cm-3]
	*          d[4] - o2 number density(cm-3]
	*          d[5] - ar number density(cm-3]
	*          d[6] - total mass density(gm/cm3] [includes anomalous oxygen]
	*          d[7] - h number density(cm-3]
	*          d[8] - n number density(cm-3]
	*          d[9] - anomalous oxygen number density(cm-3]
	*          t[1] - exospheric temperature
	*          t[2] - temperature at alt
	* ----------------------------------------------------------------------      */

	void gtd7d
		(
		msistype& msis00r,
		lpolytype& lpoly,
		fittype& fit,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong,
		double stl, double f107a, double f107, double ap[8], int mass,
		double d[10], double t[3]
		)
	{
		//------------------------------ begin ---------------------------------
		gtd7(msis00r, lpoly, fit, lsqv,
			iyd, sec, alt, glat, glong, stl, f107a, f107, ap, mass, d, t);

		// ---- total mass density
		if (mass == 48)
		{
			d[6] = 1.66e-24 * (4.0*d[1] + 16.0*d[2] + 28.0*d[3] + 32.0*d[4] +
				40.0*d[5] + d[7] + 14.0*d[8] + 16.0*d[9]);
			if (msis00r.metsel.imr == 1)
				d[6] = d[6] / 1000.0;
		}
	}

	/* ----------------------------------------------------------------------
	*         find altitude of pressure surface (press) from gtd7
	*     input:
	*          iyd - year and day as yyddd
	*          sec - ut(sec)
	*          glat - geodetic latitude(deg)
	*          glong - geodetic longitude(deg)
	*          stl - local apparent solar time(hrs)
	*          f107a - 3 month average of f10.7 flux
	*          f107 - daily f10.7 flux for previous day
	*          ap - magnetic index(daily) or when msis00r.csw.sw[9] = -1.0 :
	*             - array containing:
	*               [1] daily ap
	*               [2] 3 hr ap index for current time
	*               [3] 3 hr ap index for 3 hrs before current time
	*               [4] 3 hr ap index for 6 hrs before current time
	*               [5] 3 hr ap index for 9 hrs before current time
	*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
	*                      to current time
	*               [7] average of eight 3 hr ap indicies from 36 to 59 hrs prior
	*                      to current time
	*          press - pressure level(mb)
	*     output:
	*          alt - altitude(km)
	*          d[1] - he number density(cm-3]
	*          d[2] - o number density(cm-3]
	*          d[3] - n2 number density(cm-3]
	*          d[4] - o2 number density(cm-3]
	*          d[5] - ar number density(cm-3]
	*          d[6] - total mass density(gm/cm3]
	*          d[7] - h number density(cm-3]
	*          d[8] - n number density(cm-3]
	*          d[9] - hot o number density(cm-3]
	*          t[1] - exospheric temperature
	*          t[2] - temperature at alt
	* ----------------------------------------------------------------------      */

	void ghp7
		(
		msistype& msis00r,
		lpolytype& lpoly,
		fittype& fit,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
		double f107, double ap[8],
		double d[10], double t[3], double press
		)
	{
		double pl, zi, cl, cl2, cd, ca, z, l, xm, g, xn, p, sh, diff;
		int iday;

		double bm = 1.3806e-19;
		double rgas = 831.4;
		double test = 0.00043;
		int   ltest = 12;

		//------------------------------ begin ---------------------------------
		pl = log10(press);

		//       initial altitude estimate
		if (pl >= -5.0)
		{
			if (pl > 2.5)   zi = 18.06*(3.00 - pl);
			if (pl > 0.75 &&  pl <= 2.5) zi = 14.98*(3.08 - pl);
			if (pl > -1.0 &&  pl <= 0.75) zi = 17.8*(2.72 - pl);
			if (pl > -2.0 &&  pl <= -1.0) zi = 14.28*(3.64 - pl);
			if (pl > -4.0 &&  pl <= -2.0) zi = 12.72*(4.32 - pl);
			if (pl <= -4.0) zi = 25.3*(0.11 - pl);
			iday = iyd - int(iyd / 1000.0)*1000.0;
			printf("not iday %12i \n", iday);
			iday = iyd % 1000;
			printf("not iday %12i \n", iday);
			cl = glat / 90.0;
			cl2 = cl*cl;
			if (iday < 182)
				cd = 1.0 - iday / 91.25;
			if (iday >= 182)
				cd = iday / 91.25 - 3.0;
			ca = 0;
			if (pl > -1.11 && pl <= -.23)
				ca = 1.0;
			if (pl > -0.23)
				ca = (2.79 - pl) / (2.79 + 0.23);
			if (pl <= -1.11 && pl > -3.0)
				ca = (-2.93 - pl) / (-2.93 + 1.11);
			z = zi - 4.87*cl*cd*ca - 1.64*cl2*ca + 0.31*ca*cl;
		}
		if (pl < -5.0)
			z = 22.0*pow(pl + 4.0, 2) + 110.0;

		//       iteration loop
		l = 0;
	ten:
		l = l + 1;
		gtd7(msis00r, lpoly, fit, lsqv,
			iyd, sec, z, glat, glong, stl, f107a, f107, ap, 48, d, t);
		xn = d[1] + d[2] + d[3] + d[4] + d[5] + d[7] + d[8];
		p = bm*xn*t[2];
		if (msis00r.metsel.imr == 1)
			p = p*1.0e-6;
		diff = pl - log10(p);
		if ((fabs(diff) < test) | (l == ltest))
			goto twenty;
		xm = d[6] / xn / 1.66e-24;
		if (msis00r.metsel.imr == 1)
			xm = xm*1.0e3;
		g = msis00r.parmb.gsurf / pow((1.0 + z / msis00r.parmb.re), 2);
		sh = rgas*t[2] / (xm*g);

		// ----   new altitude estimate using scale height
		if (l < 6)
			z = z - sh*diff*2.302;
		else
			z = z - sh*diff;
		goto ten;
	twenty:
		if (l == ltest)
			printf("not converging for press %12.2f %12.2f \n", press, diff);
		alt = z;
	}

	/*  ----------------------------------------------------------------------
	*     thermospheric portion of nrlmsise-00
	*     see gtd7 for more extensive comments
	*
	*          output in m-3 and kg/m3:   call meters(.true.)
	*
	*     input variables:
	*          iyd - year and day as yyddd (day of year from 1 to 365 (or 366])
	*                (year ignored in current model)
	*          sec - ut(sec)
	*          alt - altitude(km) (>72.5 km)
	*          glat - geodetic latitude(deg)
	*          glong - geodetic longitude(deg)
	*          stl - local apparent solar time(hrs; see note below)
	*          f107a - 81 day average of f10.7 flux (centered on day ddd)
	*          f107 - daily f10.7 flux for previous day
	*          ap - magnetic index(daily) or when msis00r.csw.sw[9] = -1.0 :
	*             - array containing:
	*               [1] daily ap
	*               [2] 3 hr ap index for current time
	*               [3] 3 hr ap index for 3 hrs before current time
	*               [4] 3 hr ap index for 6 hrs before current time
	*               [5] 3 hr ap index for 9 hrs before current time
	*               [6] average of eight 3 hr ap indicies from 12 to 33 hrs prior
	*                      to current time
	*               [7] average of eight 3 hr ap indicies from 36 to 57 hrs prior
	*                      to current time
	*          mass - mass number (only density for selected gas is
	*                   calculated.  mass 0 is temperature.  mass 48 for all.
	*                   mass 17 is anomalous o only.)
	*
	*     notes on input variables:
	*          ut, local time, and longitude are used independently in the
	*          model and are not of equal importance for every situation.
	*          for the most physically realistic calculation these three
	*          variables should be consistent (stl = sec/3600+glong/15).
	*          the equation of time departures from the above formula
	*          for apparent local time can be included if available but
	*          are of minor importance.
	*
	*          f107 and f107a values used to generate the model correspond
	*          to the 10.7 cm radio flux at the actual distance of the earth
	*          from the sun rather than the radio flux at 1 au. the following
	*          site provides both classes of values:
	*          ftp:*ftp.ngdc.noaa.gov/stp/solar_data/solar_radio/flux/
	*
	*          f107, f107a, and ap effects are neither large nor well
	*          established below 80 km and these parameters should be set to
	*          150.0, 150.0, and 4.0 respectively.
	*
	*     output variables:
	*          d[1] - he number density(cm-3]
	*          d[2] - o number density(cm-3]
	*          d[3] - n2 number density(cm-3]
	*          d[4] - o2 number density(cm-3]
	*          d[5] - ar number density(cm-3]
	*          d[6] - total mass density(gm/cm3] [anomalous o not included]
	*          d[7] - h number density(cm-3]
	*          d[8] - n number density(cm-3]
	*          d[9] - anomalous oxygen number density(cm-3]
	*          t[1] - exospheric temperature
	*          t[2] - temperature at alt
	* -----------------------------------------------------------------------     */

	void gts7
		(
		msistype& msis00r,
		lpolytype& lpoly,
		lsqvtype& lsqv,
		int iyd, double sec, double alt, double glat, double glong, double stl, double f107a,
		double f107, double ap[8], int mass, double d[10], double t[3]
		)
	{
		double v2, yrd, tinf, 
			b01, b04, b14, b16, b28, b32, b40,
			db16h, g1, g14, g16, g16h, g28, g32, g4, g40,
			hc01, hc04, hc14, hc16, hc216, hc32, hc40,
			hcc01, hcc14, hcc16, hcc232, hcc32,
			rc01, rc14, rc16, rc32, t2, tho, tz, xmd, xmm,
			z, zc01, zc04, zc14, zc16, zc32, zc40,
			zcc01, zcc14, zcc16, zcc32,
			zh01, zh04, zh14, zh16, zh28, zh32, zh40, zhf,
			zhm01, zhm04, zhm14, zhm16, zhm28, zhm32, zhm40, zmho, zsho, zsht;
		int i;

		double mt[12] = { 0, 48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17 };
		int mn1 = 5;
		double zn1[6] = { 0, 120.0, 110.0, 100.0, 90.0, 72.5 };

		double dgtr = 1.74533e-2;
		double dr = 1.72142e-2;
		double alpha[10] = { 0, -0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0 };
		double altl[9] = { 0, 200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0 };

		// treat as local, but save in between calls
		static double alast = -999.0;

		//------------------------------ begin ---------------------------------
		// ----  test for changed input
		v2 = vtst(msis00r.csw, iyd, sec, glat, glong, stl, f107a, f107, ap, 2);

		yrd = iyd;
		msis00r.gts3c.za = msis00r.parm.pdl[16][2];
		zn1[1] = msis00r.gts3c.za;
		for (i = 1; i <= 9; i++)
			d[i] = 0.0;

		// ----  tinf variations not important below za or zn1[1]
		if (alt > zn1[1])
		{
			if ((v2 == 1.0) | (alast <= zn1[1]))
				tinf = msis00r.lower.ptm[1] * msis00r.parm.pt[1] *
				(1.0 + msis00r.csw.sw[16] *
				globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107,
				ap, msis00r.parm.pt));
		}
		else
			tinf = msis00r.lower.ptm[1] * msis00r.parm.pt[1];

		t[1] = tinf;
		// ----    gradient variations not important below zn1[5]
		if (alt > zn1[5])
		{
			if ((v2 == 1) | (alast <= zn1[5]))
				msis00r.gts3c.g0 = msis00r.lower.ptm[4] * msis00r.parm.ps[1] * (1.0 + msis00r.csw.sw[19] *
				globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107,
				ap, msis00r.parm.ps));
		}
		else
			msis00r.gts3c.g0 = msis00r.lower.ptm[4] * msis00r.parm.ps[1];

		//  calculate these temperatures only if input changed
		if ((v2 == 1.0) | (alt < 300.0))
			msis00r.gts3c.tlb = msis00r.lower.ptm[2] *
			(1.0 + msis00r.csw.sw[17] *
			globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107,
			ap, msis00r.parm.pddav[4])) * msis00r.parm.pd[1][4];

		msis00r.gts3c.s = msis00r.gts3c.g0 / (tinf - msis00r.gts3c.tlb);
		// ---- lower thermosphere temp variations not significant for
		// ----  density above 300 km
		if (alt < 300.0)
		{
			if ((v2 == 1.0) | (alast >= 300.0))
			{
				msis00r.meso.tn1[2] = msis00r.lower.ptm[7] * msis00r.parm.ptl[1][1] /
					(1.0 - msis00r.csw.sw[18] * glob7s(lpoly, msis00r.csw, msis00r.parm.ptldav[1]));
				msis00r.meso.tn1[3] = msis00r.lower.ptm[3] * msis00r.parm.ptl[1][2] /
					(1.0 - msis00r.csw.sw[18] * glob7s(lpoly, msis00r.csw, msis00r.parm.ptldav[2]));
				msis00r.meso.tn1[4] = msis00r.lower.ptm[8] * msis00r.parm.ptl[1][3] /
					(1.0 - msis00r.csw.sw[18] * glob7s(lpoly, msis00r.csw, msis00r.parm.ptldav[3]));
				msis00r.meso.tn1[5] = msis00r.lower.ptm[5] * msis00r.parm.ptl[1][4] /
					(1.0 - msis00r.csw.sw[18] * msis00r.csw.sw[20] *
					glob7s(lpoly, msis00r.csw, msis00r.parm.ptldav[4]));
				msis00r.meso.tgn1[2] = msis00r.lower.ptm[9] * msis00r.parm.pma[1][9] *
					(1.0 + msis00r.csw.sw[18] * msis00r.csw.sw[20] *
					glob7s(lpoly, msis00r.csw, msis00r.parm.pmadav[9])) *
					msis00r.meso.tn1[5] * msis00r.meso.tn1[5] /
					pow((msis00r.lower.ptm[5] * msis00r.parm.ptl[1][4]), 2);
			}
		}
		else
		{
			msis00r.meso.tn1[2] = msis00r.lower.ptm[7] * msis00r.parm.ptl[1][1];
			msis00r.meso.tn1[3] = msis00r.lower.ptm[3] * msis00r.parm.ptl[1][2];
			msis00r.meso.tn1[4] = msis00r.lower.ptm[8] * msis00r.parm.ptl[1][3];
			msis00r.meso.tn1[5] = msis00r.lower.ptm[5] * msis00r.parm.ptl[1][4];
			msis00r.meso.tgn1[2] = msis00r.lower.ptm[9] * msis00r.parm.pma[1][9] * msis00r.meso.tn1[5] *
				msis00r.meso.tn1[5] /
				pow((msis00r.lower.ptm[5] * msis00r.parm.ptl[1][4]), 2);
		}

		msis00r.gts3c.z0 = zn1[4];
		msis00r.gts3c.t0 = msis00r.meso.tn1[4];
		msis00r.gts3c.tr12 = 1.0;

		if (mass == 0)
			goto fifty;

		// ---- n2 variation factor at zlb
		g28 = msis00r.csw.sw[21] *
			globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap, msis00r.parm.pddav[3]);
		lpoly.day = int(yrd) % 1000;

		// ----  variation of turbopause height
		zhf = msis00r.parm.pdl[25][2] * (1.0 + msis00r.csw.sw[5] * msis00r.parm.pdl[25][1] * sin(dgtr*glat)*
			cos(dr*(lpoly.day - msis00r.parm.pt[14])));
		yrd = iyd;
		t[1] = tinf;
		xmm = msis00r.lower.pdm[5][3];
		z = alt;

		for (i = 1; i <= 11; i++)
		if (mass == mt[i])   goto fifteen;

		printf("mass %i5 not valid\n", mass);
		goto ninety;
	fifteen:
		if ((z <= altl[6]) | (mass == 28) | (mass == 48))
		{
			// ---- **** n2 density ****
			//       diffusive density at zlb
			msis00r.gts3c.db28 = msis00r.lower.pdm[1][3] * exp(g28)*msis00r.parm.pd[1][3];
			//       diffusive density at alt
			d[3] = densu(msis00r.parmb, lsqv, z, msis00r.gts3c.db28, tinf, msis00r.gts3c.tlb,
				28.0, alpha[3], t[2], msis00r.lower.ptm[6], msis00r.gts3c.s,
				mn1, zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
			msis00r.gts3c.dd = d[3];
			//       turbopause
			zh28 = msis00r.lower.pdm[3][3] * zhf;
			zhm28 = msis00r.lower.pdm[4][3] * msis00r.parm.pdl[6][2];
			xmd = 28.0 - xmm;
			//       mixed density at zlb
			b28 = densu(msis00r.parmb, lsqv, zh28, msis00r.gts3c.db28, tinf, msis00r.gts3c.tlb,
				xmd, alpha[3] - 1.0, tz, msis00r.lower.ptm[6],
				msis00r.gts3c.s, mn1, zn1, msis00r.meso.tn1, msis00r.meso.tgn1);

			if ((z <= altl[3]) && (msis00r.csw.sw[15] != 0))
			{
				//       mixed density at alt
				msis00r.dmix.dm28 = densu(msis00r.parmb, lsqv, z, b28, tinf, msis00r.gts3c.tlb,
					xmm, alpha[3], tz, msis00r.lower.ptm[6], msis00r.gts3c.s, mn1,
					zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
				//       net density at alt
				d[3] = dnet(d[3], msis00r.dmix.dm28, zhm28, xmm, 28.0);
			}
		}
		switch (i)
		{
		case 1: goto twenty;
		case 2: goto fifty;
		case 3: goto twenty;
		case 4: goto twentyfive;
		case 5: goto ninety;
		case 6: goto thirtyfive;
		case 7: goto fourty;
		case 8: goto fourtyfive;
		case 9: goto twentyfive;
		case 10: goto fourtyeight;
		case 11: goto fourtysix;
		}
	twenty:
		// ---- **** he density ****
		// ---- density variation factor at zlb
		g4 = msis00r.csw.sw[21] *
			globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap, msis00r.parm.pddav[1]);

		//       diffusive density at zlb
		msis00r.gts3c.db04 = msis00r.lower.pdm[1][1] * exp(g4)*msis00r.parm.pd[1][1];

		//      diffusive density at alt
		d[1] = densu(msis00r.parmb, lsqv, z, msis00r.gts3c.db04, tinf, msis00r.gts3c.tlb, 4.0, alpha[1],
			t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1,
			zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
		msis00r.gts3c.dd = d[1];
		if (z <= altl[1] && msis00r.csw.sw[15] != 0)
		{
			//      turbopause
			zh04 = msis00r.lower.pdm[3][1];

			//      mixed density at zlb
			b04 = densu(msis00r.parmb, lsqv, zh04, msis00r.gts3c.db04, tinf, msis00r.gts3c.tlb,
				4.0 - xmm, alpha[1] - 1.0, t[2], msis00r.lower.ptm[6], msis00r.gts3c.s,
				mn1, zn1, msis00r.meso.tn1, msis00r.meso.tgn1);

			//      mixed density at alt
			msis00r.dmix.dm04 = densu(msis00r.parmb, lsqv, z, b04, tinf, msis00r.gts3c.tlb,
				xmm, 0.0, t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1,
				zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
			zhm04 = zhm28;

			//      net density at alt
			d[1] = dnet(d[1], msis00r.dmix.dm04, zhm04, xmm, 4.0);

			//      correction to specified mixing ratio at ground
			msis00r.gts3c.rl = log(b28*msis00r.lower.pdm[2][1] / b04);
			zc04 = msis00r.lower.pdm[5][1] * msis00r.parm.pdl[1][2];
			hc04 = msis00r.lower.pdm[6][1] * msis00r.parm.pdl[2][2];

			//      net density corrected at alt
			d[1] = d[1] * ccor(z, msis00r.gts3c.rl, hc04, zc04);
		}

		if (mass != 48)   goto ninety;

	twentyfive:
		// ----**** o density ****
		// ---- density variation factor at zlb
		g16 = msis00r.csw.sw[21] *
			globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis00r.parm.pddav[2]);
		//       diffusive density at zlb
		msis00r.gts3c.db16 = msis00r.lower.pdm[1][2] * exp(g16)*msis00r.parm.pd[1][2];
		// ---- diffusive density at alt
		d[2] = densu(msis00r.parmb, lsqv, z, msis00r.gts3c.db16, tinf, msis00r.gts3c.tlb,
			16.0, alpha[2], t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1,
			zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
		msis00r.gts3c.dd = d[2];

		if (z <= altl[2] && msis00r.csw.sw[15] != 0)
		{
			//  corrected from msis00r.lower.pdm[3,1) to msis00r.lower.pdm[3][2]  12/2/85
			// ---- turbopause
			zh16 = msis00r.lower.pdm[3][2];
			//       mixed density at zlb
			b16 = densu(msis00r.parmb, lsqv, zh16, msis00r.gts3c.db16, tinf, msis00r.gts3c.tlb,
				16.0 - xmm, alpha[2] - 1.0, t[2], msis00r.lower.ptm[6], msis00r.gts3c.s,
				mn1, zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
			//       mixed density at alt
			msis00r.dmix.dm16 = densu(msis00r.parmb, lsqv, z, b16, tinf, msis00r.gts3c.tlb, xmm, 0.0,
				t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1, zn1, msis00r.meso.tn1,
				msis00r.meso.tgn1);
			zhm16 = zhm28;
			//       net density at alt
			d[2] = dnet(d[2], msis00r.dmix.dm16, zhm16, xmm, 16.0);
			//   3/16/99 change form to match o2 departure from diff equil near 150
			//   km and add dependence on f10.7
			//       rl = dlog(b28*msis00r.lower.pdm[2][2]*fabs(msis00r.parm.pdl[17][2])/b16]
			msis00r.gts3c.rl = msis00r.lower.pdm[2][2] * msis00r.parm.pdl[17][2] *
				(1.0 + msis00r.csw.sw[1] * msis00r.parm.pdl[24][1] * (f107a - 150.0));
			hc16 = msis00r.lower.pdm[6][2] * msis00r.parm.pdl[4][2];
			zc16 = msis00r.lower.pdm[5][2] * msis00r.parm.pdl[3][2];
			hc216 = msis00r.lower.pdm[6][2] * msis00r.parm.pdl[5][2];
			d[2] = d[2] * ccor2(z, msis00r.gts3c.rl, hc16, zc16, hc216);

			// ---- chemistry correction
			hcc16 = msis00r.lower.pdm[8][2] * msis00r.parm.pdl[14][2];
			zcc16 = msis00r.lower.pdm[7][2] * msis00r.parm.pdl[13][2];
			rc16 = msis00r.lower.pdm[4][2] * msis00r.parm.pdl[15][2];
			//       net density corrected at alt
			d[2] = d[2] * ccor(z, rc16, hcc16, zcc16);
		}
		if (mass != 48 && mass != 49) goto ninety;

	thirtyfive:
		// ---- **** o2 density ****
		// ---- density variation factor at zlb
		g32 = msis00r.csw.sw[21] *
			globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap,
			msis00r.parm.pddav[5]);

		//       diffusive density at zlb
		msis00r.gts3c.db32 = msis00r.lower.pdm[1][4] * exp(g32)*msis00r.parm.pd[1][5];

		// ---- diffusive density at alt
		d[4] = densu(msis00r.parmb, lsqv, z, msis00r.gts3c.db32, tinf, msis00r.gts3c.tlb,
			32.0, alpha[4], t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1,
			zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
		if (mass == 49)
			msis00r.gts3c.dd = msis00r.gts3c.dd + 2.0*d[4];
		else
			msis00r.gts3c.dd = d[4];

		if (msis00r.csw.sw[15] != 0)
		{
			if (z <= altl[4])
			{
				// ---- turbopause
				zh32 = msis00r.lower.pdm[3][4];
				//       mixed density at zlb
				b32 = densu(msis00r.parmb, lsqv, zh32, msis00r.gts3c.db32, tinf, msis00r.gts3c.tlb,
					32.0 - xmm, alpha[4] - 1.0, t[2], msis00r.lower.ptm[6], msis00r.gts3c.s,
					mn1, zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
				//       mixed density at alt
				msis00r.dmix.dm32 = densu(msis00r.parmb, lsqv, z, b32, tinf, msis00r.gts3c.tlb, xmm, 0.0, t[2],
					msis00r.lower.ptm[6], msis00r.gts3c.s, mn1, zn1, msis00r.meso.tn1,
					msis00r.meso.tgn1);
				zhm32 = zhm28;
				//       net density at alt
				d[4] = dnet(d[4], msis00r.dmix.dm32, zhm32, xmm, 32.0);

				// ---- correction to specified mixing ratio at ground
				msis00r.gts3c.rl = log(b28*msis00r.lower.pdm[2][4] / b32);
				hc32 = msis00r.lower.pdm[6][4] * msis00r.parm.pdl[8][2];
				zc32 = msis00r.lower.pdm[5][4] * msis00r.parm.pdl[7][2];
				d[4] = d[4] * ccor(z, msis00r.gts3c.rl, hc32, zc32);
			}

			//       correction for general departure from diffusive equilibrium above zlb
			hcc32 = msis00r.lower.pdm[8][4] * msis00r.parm.pdl[23][2];
			hcc232 = msis00r.lower.pdm[8][4] * msis00r.parm.pdl[23][1];
			zcc32 = msis00r.lower.pdm[7][4] * msis00r.parm.pdl[22][2];
			rc32 = msis00r.lower.pdm[4][4] * msis00r.parm.pdl[24][2] * (
				1.0 + msis00r.csw.sw[1] * msis00r.parm.pdl[24][1] * (f107a - 150.0));

			//       net density corrected at alt
			d[4] = d[4] * ccor2(z, rc32, hcc32, zcc32, hcc232);
		}

		if (mass != 48)   goto ninety;

	fourty:
		// ---- **** ar density ****
		// ---- density variation factor at zlb
		g40 = msis00r.csw.sw[21] *
			globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap, msis00r.parm.pddav[6]);
		//      diffusive density at zlb
		msis00r.gts3c.db40 = msis00r.lower.pdm[1][5] * exp(g40)*msis00r.parm.pd[1][6];

		// ---- diffusive density at alt
		d[5] = densu(msis00r.parmb, lsqv, z, msis00r.gts3c.db40, tinf, msis00r.gts3c.tlb,
			40.0, alpha[5], t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1,
			zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
		msis00r.gts3c.dd = d[5];

		if (z <= altl[5] && msis00r.csw.sw[15] != 0)
		{
			// ---- turbopause
			zh40 = msis00r.lower.pdm[3][5];
			//       mixed density at zlb
			b40 = densu(msis00r.parmb, lsqv, zh40, msis00r.gts3c.db40, tinf, msis00r.gts3c.tlb,
				40.0 - xmm, alpha[5] - 1.0, t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1, zn1,
				msis00r.meso.tn1, msis00r.meso.tgn1);
			//       mixed density at alt
			msis00r.dmix.dm40 = densu(msis00r.parmb, lsqv, z, b40, tinf, msis00r.gts3c.tlb, xmm,
				0.0, t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1, zn1,
				msis00r.meso.tn1, msis00r.meso.tgn1);
			zhm40 = zhm28;
			//       net density at alt
			d[5] = dnet(d[5], msis00r.dmix.dm40, zhm40, xmm, 40.0);

			// ---- correction to specified mixing ratio at ground
			msis00r.gts3c.rl = log(b28*msis00r.lower.pdm[2][5] / b40);
			hc40 = msis00r.lower.pdm[6][5] * msis00r.parm.pdl[10][2];
			zc40 = msis00r.lower.pdm[5][5] * msis00r.parm.pdl[9][2];

			//       net density corrected at alt
			d[5] = d[5] * ccor(z, msis00r.gts3c.rl, hc40, zc40);
		}

		if (mass != 48)   goto ninety;

	fourtyfive:
		// ----  **** hydrogen density ****
		// ---- density variation factor at zlb
		g1 = msis00r.csw.sw[21] *
			globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap, msis00r.parm.pddav[7]);
		//      diffusive density at zlb
		msis00r.gts3c.db01 = msis00r.lower.pdm[1][6] * exp(g1)*msis00r.parm.pd[1][7];

		// ---- diffusive density at alt
		d[7] = densu(msis00r.parmb, lsqv, z, msis00r.gts3c.db01, tinf, msis00r.gts3c.tlb,
			1.0, alpha[7], t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1,
			zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
		msis00r.gts3c.dd = d[7];

		if (z <= altl[7] && msis00r.csw.sw[15] != 0)
		{
			// ---- turbopause
			zh01 = msis00r.lower.pdm[3][6];
			//       mixed density at zlb
			b01 = densu(msis00r.parmb, lsqv, zh01, msis00r.gts3c.db01, tinf, msis00r.gts3c.tlb,
				1.0 - xmm, alpha[7] - 1.0, t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1, zn1,
				msis00r.meso.tn1, msis00r.meso.tgn1);
			//       mixed density at alt
			msis00r.dmix.dm01 = densu(msis00r.parmb, lsqv, z, b01, tinf, msis00r.gts3c.tlb, xmm,
				0.0, t[2], msis00r.lower.ptm[6], msis00r.gts3c.s,
				mn1, zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
			zhm01 = zhm28;

			//       net density at alt
			d[7] = dnet(d[7], msis00r.dmix.dm01, zhm01, xmm, 1.0);

			// ---- correction to specified mixing ratio at ground
			msis00r.gts3c.rl = log(b28*msis00r.lower.pdm[2][6] * fabs(msis00r.parm.pdl[18][2]) / b01);
			hc01 = msis00r.lower.pdm[6][6] * msis00r.parm.pdl[12][2];
			zc01 = msis00r.lower.pdm[5][6] * msis00r.parm.pdl[11][2];
			d[7] = d[7] * ccor(z, msis00r.gts3c.rl, hc01, zc01);

			// ---- chemistry correction
			hcc01 = msis00r.lower.pdm[8][6] * msis00r.parm.pdl[20][2];
			zcc01 = msis00r.lower.pdm[7][6] * msis00r.parm.pdl[19][2];
			rc01 = msis00r.lower.pdm[4][6] * msis00r.parm.pdl[21][2];

			//      net density corrected at alt
			d[7] = d[7] * ccor(z, rc01, hcc01, zcc01);
		}

		if (mass != 48)   goto ninety;

	fourtyeight:
		// ----  **** atomic nitrogen density ****
		// ---- density variation factor at zlb
		g14 = msis00r.csw.sw[21] *
			globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap, msis00r.parm.pddav[8]);

		//       diffusive density at zlb
		msis00r.gts3c.db14 = msis00r.lower.pdm[1][7] * exp(g14)*msis00r.parm.pd[1][8];

		// ---- diffusive density at alt
		d[8] = densu(msis00r.parmb, lsqv, z, msis00r.gts3c.db14, tinf, msis00r.gts3c.tlb, 14.0, alpha[8],
			t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1,
			zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
		msis00r.gts3c.dd = d[8];

		if (z <= altl[8] && msis00r.csw.sw[15] != 0)
		{
			// ---- turbopause
			zh14 = msis00r.lower.pdm[3][7];

			//       mixed density at zlb
			b14 = densu(msis00r.parmb, lsqv, zh14, msis00r.gts3c.db14, tinf, msis00r.gts3c.tlb,
				14.0 - xmm, alpha[8] - 1.0, t[2], msis00r.lower.ptm[6], msis00r.gts3c.s, mn1, zn1,
				msis00r.meso.tn1, msis00r.meso.tgn1);

			//       mixed density at alt
			msis00r.dmix.dm14 = densu(msis00r.parmb, lsqv, z, b14, tinf, msis00r.gts3c.tlb, xmm, 0.0,
				t[2], msis00r.lower.ptm[6],
				msis00r.gts3c.s, mn1, zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
			zhm14 = zhm28;
			//       net density at alt
			d[8] = dnet(d[8], msis00r.dmix.dm14, zhm14, xmm, 14.0);

			// ---- correction to specified mixing ratio at ground
			msis00r.gts3c.rl = log(b28*msis00r.lower.pdm[2][7] * fabs(msis00r.parm.pdl[3][1]) / b14);
			hc14 = msis00r.lower.pdm[6][7] * msis00r.parm.pdl[2][1];
			zc14 = msis00r.lower.pdm[5][7] * msis00r.parm.pdl[1][1];
			d[8] = d[8] * ccor(z, msis00r.gts3c.rl, hc14, zc14);

			// ---- chemistry correction
			hcc14 = msis00r.lower.pdm[8][7] * msis00r.parm.pdl[5][1];
			zcc14 = msis00r.lower.pdm[7][7] * msis00r.parm.pdl[4][1];
			rc14 = msis00r.lower.pdm[4][7] * msis00r.parm.pdl[6][1];

			//       net density corrected at alt
			d[8] = d[8] * ccor(z, rc14, hcc14, zcc14);
		}

		if (mass != 48) goto ninety;

	fourtysix:
		// ----  **** anomalous oxygen density ****
		g16h = msis00r.csw.sw[21] *
			globe7(msis00r.csw, lpoly, yrd, sec, glat, glong, stl, f107a, f107, ap, msis00r.parm.pddav[9]);
		db16h = msis00r.lower.pdm[1][8] * exp(g16h)*msis00r.parm.pd[1][9];
		tho = msis00r.lower.pdm[10][8] * msis00r.parm.pdl[7][1];
		msis00r.gts3c.dd = densu(msis00r.parmb, lsqv, z, db16h, tho, tho, 16.0, alpha[9], t2,
			msis00r.lower.ptm[6], msis00r.gts3c.s, mn1,
			zn1, msis00r.meso.tn1, msis00r.meso.tgn1);
		zsht = msis00r.lower.pdm[6][8];
		zmho = msis00r.lower.pdm[5][8];
		zsho = scalh(msis00r.parmb, zmho, 16.0, tho);
		d[9] = msis00r.gts3c.dd * exp(-zsht / zsho*(exp(-(z - zmho) / zsht) - 1.0));

		if (mass != 48) goto ninety;

		// ---- total mass density
		d[6] = 1.66e-24 *
			(4.0*d[1] + 16.0*d[2] + 28.0*d[3] + 32.0*d[4] + 40.0*d[5] + d[7] + 14.0*d[8]);
		msis00r.gts3c.db48 = 1.66e-24*(4.0*msis00r.gts3c.db04 + 16.0*msis00r.gts3c.db16 +
			28.0*msis00r.gts3c.db28 + 32.0*msis00r.gts3c.db32 +
			40.0*msis00r.gts3c.db40 + msis00r.gts3c.db01 + 14.0*msis00r.gts3c.db14);

		goto ninety;
		// ---- temperature at altitude
	fifty:
		z = fabs(alt);
		//cdav
		//        ddum is never used so this call is unnecessary
		//        ddum  = densu(msis00r.parmb, lsqv, z,1.0, tinf,msis00r.gts3c.tlb,0.0,0.0,t[2],
		//                  msis00r.lower.ptm[6],msis00r.gts3c.msis00r.gts3c.s,mn1,
		//                  zn1,msis00r.meso.tn1,msis00r.meso.tgn1);
	ninety:
		// ---- adjust densities from cgs to kgm
		if (msis00r.metsel.imr == 1)
		{
			for (i = 1; i <= 9; i++)
				d[i] = d[i] * 1.0e6;

			d[6] = d[6] / 1000.0;
		}

		alast = alt;
	}

	/* ----------------------------------------------------------------------
	*      calculate scale height (km)
	* ----------------------------------------------------------------------      */

	double scalh
		(
		parmbtype& parmb,
		double alt, double xm, double temp
		)
	{
		double g;
		double rgas = 831.4;

		g = parmb.gsurf / pow((1.0 + alt / parmb.re), 2);
		return rgas*temp / (g*xm);
	}

	//// cdav these next few functions could be done as inline functions if that speeds
	//// things up
	//// ---- 3hr magnetic activity functions
	////      eq. a24d
	//double g0
	//       (double a, double p[151])
	//          {
	//          return (a - 4.0 + (p[26]-1.0) *
	//                 (a - 4.0 + (exp(-fabs(p[25])*(a-4.0)) -1.0) / fabs(p[25]) ));
	//          }
	//// ---- eq. a24c
	//double sumex
	//       (double ex)
	//          {
	//          return 1.0 + (1.0-pow(ex,19)) / (1.0-ex)*sqrt(ex);
	//          }
	//// ---- eq. a24a
	//double sg0
	//       (double ex, double p[151], double ap[8])
	//          {
	//          return ( g0(ap[2],p) +
	//                  ( g0(ap[3],p)*ex + g0(ap[4],p)*ex*ex + g0(ap[5],p)* pow(ex,3) +
	//                   ( g0(ap[6],p)*pow(ex,4) + g0(ap[7],p)*pow(ex,12)) * (1.0-pow(ex,8))
	//                   / (1.0-ex)
	//                  )
	//                 ) / sumex(ex);
	//          }
	/* ----------------------------------------------------------------------
	*         calculate g(l) function
	*         upper thermosphere parameters
	* ----------------------------------------------------------------------      */

	double globe7
		(
		cswtype& csw,
		lpolytype& lpoly,
		double yrd, double sec, double lat, double llong, double tloc, double f107a,
		double f107, double ap[8], double p[151]
		)
	{
		int i, iyr;
		double c, s, c2, c4, s2, cd14, cd18, cd32, cd39, f1, f2, t71, t72,
			t81, t82, p44, p45, exp1;
		double t[16], tinfg;

		double dgtr = 1.74533e-2;
		double dr = 1.72142e-2;
		double xl = 1000.0;
		static double tll = 1000.0;
		static double sw9 = 1.0;
		double dayl = -1.0;
		double p14 = -1000.0;
		double p18 = -1000.0;
		double p32 = -1000.0;
		double hr = 0.2618;
		double sr = 7.2722e-5;
		int sv[26] = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
		double nsw = 14;
		double p39 = -1000.0;

		//------------------------------ begin ---------------------------------
		if (csw.isw != 64999)
			tselec(csw, sv);

		for (i = 1; i <= 15; i++)
			t[i] = 0.0;

		if (csw.sw[9] > 0)
			sw9 = 1.0;
		if (csw.sw[9] < 0)
			sw9 = -1.0;
		iyr = yrd / 1000.0;
		lpoly.day = yrd - iyr*1000.0;
		lpoly.xlong = llong;
		//  eq. a22 (remainder of code)
		if (xl != lat)
		{
			// ----    calculate legendre polynomials
			c = sin(lat*dgtr);
			s = cos(lat*dgtr);
			c2 = c*c;
			c4 = c2*c2;
			s2 = s*s;
			lpoly.plg[2][1] = c;
			lpoly.plg[3][1] = 0.5*(3.0*c2 - 1.0);
			lpoly.plg[4][1] = 0.5*(5.0*c*c2 - 3.0*c);
			lpoly.plg[5][1] = (35.0*c4 - 30.0*c2 + 3.0) / 8.0;
			lpoly.plg[6][1] = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c) / 8.0;
			lpoly.plg[7][1] = (11.0*c*lpoly.plg[6][1] - 5.0*lpoly.plg[5][1]) / 6.0;
			//        lpoly.plg[8][1] = (13.0*c*lpoly.plg[7][1] - 6.0*lpoly.plg[6][1])/7.0;
			lpoly.plg[2][2] = s;
			lpoly.plg[3][2] = 3.0*c*s;
			lpoly.plg[4][2] = 1.5*(5.0*c2 - 1.0)*s;
			lpoly.plg[5][2] = 2.5*(7.0*c2*c - 3.0*c)*s;
			lpoly.plg[6][2] = 1.875*(21.0*c4 - 14.0*c2 + 1.0)*s;
			lpoly.plg[7][2] = (11.0*c*lpoly.plg[6][2] - 6.0*lpoly.plg[5][2]) / 5.0;
			//        lpoly.plg[8][2] = (13.0*c*lpoly.plg[7][2]-7.0*lpoly.plg[6][2])/6.0;
			//        lpoly.plg[9][2] = (15.0*c*lpoly.plg[8][2]-8.0*lpoly.plg[7][2])/7.0;
			lpoly.plg[3][3] = 3.0*s2;
			lpoly.plg[4][3] = 15.0*s2*c;
			lpoly.plg[5][3] = 7.5*(7.0*c2 - 1.0)*s2;
			lpoly.plg[6][3] = 3.0*c*lpoly.plg[5][3] - 2.0*lpoly.plg[4][3];
			lpoly.plg[7][3] = (11.0*c*lpoly.plg[6][3] - 7.0*lpoly.plg[5][3]) / 4.0;
			lpoly.plg[8][3] = (13.0*c*lpoly.plg[7][3] - 8.0*lpoly.plg[6][3]) / 5.0;
			lpoly.plg[4][4] = 15.0*s2*s;
			lpoly.plg[5][4] = 105.0*s2*s*c;
			lpoly.plg[6][4] = (9.0*c*lpoly.plg[5][4] - 7.0*lpoly.plg[4][4]) / 2.0;
			lpoly.plg[7][4] = (11.0*c*lpoly.plg[6][4] - 8.0*lpoly.plg[5][4]) / 3.0;

			xl = lat;
		}
		if (tll == tloc) goto sixteen;
		if (csw.sw[7] == 0 && csw.sw[8] == 0 && csw.sw[14] == 0) goto sixteen;
		lpoly.stloc = sin(hr*tloc);
		lpoly.ctloc = cos(hr*tloc);
		lpoly.s2tloc = sin(2.0*hr*tloc);
		lpoly.c2tloc = cos(2.0*hr*tloc);
		lpoly.s3tloc = sin(3.0*hr*tloc);
		lpoly.c3tloc = cos(3.0*hr*tloc);
		tll = tloc;
	sixteen:
		if ((lpoly.day != dayl) | (p[14] != p14))
			cd14 = cos(dr*(lpoly.day - p[14]));
		if ((lpoly.day != dayl) | (p[18] != p18))
			cd18 = cos(2.0*dr*(lpoly.day - p[18]));

		if ((lpoly.day != dayl) | (p[32] != p32))
			cd32 = cos(dr*(lpoly.day - p[32]));
		if ((lpoly.day != dayl) | (p[39] != p39))
			cd39 = cos(2.0*dr*(lpoly.day - p[39]));

		dayl = lpoly.day;
		p14 = p[14];
		p18 = p[18];
		p32 = p[32];
		p39 = p[39];

		// ----   f10.7 effect
		lpoly.df = f107 - f107a;
		lpoly.dfa = f107a - 150.0;
		t[1] = p[20] * lpoly.df*(1.0 + p[60] * lpoly.dfa) + p[21] * lpoly.df*lpoly.df +
			p[22] * lpoly.dfa + p[30] * lpoly.dfa*lpoly.dfa;
		f1 = 1.0 + (p[48] * lpoly.dfa + p[20] * lpoly.df + p[21] * lpoly.df*lpoly.df)*csw.swc[1];
		f2 = 1.0 + (p[50] * lpoly.dfa + p[20] * lpoly.df + p[21] * lpoly.df*lpoly.df)*csw.swc[1];

		// ----  time independent
		t[2] = (p[2] * lpoly.plg[3][1] + p[3] * lpoly.plg[5][1] + p[23] * lpoly.plg[7][1])
			+ (p[15] * lpoly.plg[3][1])*lpoly.dfa*csw.swc[1] + p[27] * lpoly.plg[2][1];

		// ----  symmetrical annual
		t[3] = (p[19])*cd32;

		// ----  symmetrical semiannual
		t[4] = (p[16] + p[17] * lpoly.plg[3][1])*cd18;

		// ----  asymmetrical annual
		t[5] = f1* (p[10] * lpoly.plg[2][1] + p[11] * lpoly.plg[4][1])*cd14;

		// ----   asymmetrical semiannual
		t[6] = p[38] * lpoly.plg[2][1] * cd39;

		// ----  diurnal
		if (csw.sw[7] != 0)
		{
			t71 = (p[12] * lpoly.plg[3][2])*cd14*csw.swc[5];
			t72 = (p[13] * lpoly.plg[3][2])*cd14*csw.swc[5];
			t[7] = f2* ((p[4] * lpoly.plg[2][2] + p[5] * lpoly.plg[4][2] + p[28] * lpoly.plg[6][2]
				+ t71)*lpoly.ctloc + (p[7] * lpoly.plg[2][2] + p[8] * lpoly.plg[4][2]
				+ p[29] * lpoly.plg[6][2] + t72)*lpoly.stloc);
		}
		// ----  semidiurnal
		if (csw.sw[8] != 0)
		{
			t81 = (p[24] * lpoly.plg[4][3] + p[36] * lpoly.plg[6][3])*cd14*csw.swc[5];
			t82 = (p[34] * lpoly.plg[4][3] + p[37] * lpoly.plg[6][3])*cd14*csw.swc[5];
			t[8] = f2*((p[6] * lpoly.plg[3][3] + p[42] * lpoly.plg[5][3] + t81)*lpoly.c2tloc
				+ (p[9] * lpoly.plg[3][3] + p[43] * lpoly.plg[5][3] + t82)*lpoly.s2tloc);
		}
		// ----  terdiurnal
		if (csw.sw[14] != 0)
		{
			t[14] = f2*
				((p[40] * lpoly.plg[4][4] + (p[94] * lpoly.plg[5][4]
				+ p[47] * lpoly.plg[7][4])*cd14*csw.swc[5])*lpoly.s3tloc
				+ (p[41] * lpoly.plg[4][4] + (p[95] * lpoly.plg[5][4]
				+ p[49] * lpoly.plg[7][4])*cd14*csw.swc[5])*lpoly.c3tloc);
		}

		// ----    magnetic activity based on daily ap
		if (sw9 != -1.0)
		{
			lpoly.apd = (ap[1] - 4.0);
			p44 = p[44];
			p45 = p[45];
			// cdav at one time i saw an error where p44 was less than zero - unfortunately, i didn't
			//      write exactly how i got there. it may have been during debugging when i had other
			//      things wrong
			if (p44 <= 0.0)
				p44 = 1.0e-5;
			lpoly.apdf = lpoly.apd + (p45 - 1.0)*(lpoly.apd + (exp(-p44  *lpoly.apd) - 1.0) / p44);
			if (csw.sw[9] == 0)
				goto fourty;
			t[9] = lpoly.apdf*(p[33] + p[46] * lpoly.plg[3][1] + p[35] * lpoly.plg[5][1] +
				(p[101] * lpoly.plg[2][1] + p[102] * lpoly.plg[4][1] + p[103] * lpoly.plg[6][1])*cd14*csw.swc[5] +
				(p[122] * lpoly.plg[2][2] + p[123] * lpoly.plg[4][2] + p[124] * lpoly.plg[6][2])*csw.swc[7] *
				cos(hr*(tloc - p[125])));
			goto fourty;
		}

		if (p[52] == 0)
			goto fourty;
		exp1 = exp(-10800.0*fabs(p[52]) / (1.0 + p[139] * (45.0 - fabs(lat))));
		if (exp1 > .99999)
			exp1 = 0.99999;
		if (p[25] < 1.0e-4)
			p[25] = 1.0e-4;
		lpoly.apt[1] = sg0(exp1, p, ap);
		//        apt[2] = sg2[exp1]
		//        apt[3] = sg0[exp2]
		//        apt[4] = sg2[exp2]

		if (csw.sw[9] == 0) goto fourty;
		t[9] = lpoly.apt[1] * (p[51] + p[97] * lpoly.plg[3][1] + p[55] * lpoly.plg[5][1] +
			(p[126] * lpoly.plg[2][1] + p[127] * lpoly.plg[4][1] + p[128] * lpoly.plg[6][1])*cd14*csw.swc[5] +
			(p[129] * lpoly.plg[2][2] + p[130] * lpoly.plg[4][2] + p[131] * lpoly.plg[6][2])*csw.swc[7] *
			cos(hr*(tloc - p[132])));
	fourty:

		if ((csw.sw[10] == 0) | (llong <= -1000.0))
			goto fourtynine;

		// ----  longitudinal
		if (csw.sw[11] != 0)
		{
			t[11] = (1.0 + p[81] * lpoly.dfa*csw.swc[1])*
				((p[65] * lpoly.plg[3][2] + p[66] * lpoly.plg[5][2] + p[67] * lpoly.plg[7][2]
				+ p[104] * lpoly.plg[2][2] + p[105] * lpoly.plg[4][2] + p[106] * lpoly.plg[6][2]
				+ csw.swc[5] * (p[110] * lpoly.plg[2][2] + p[111] * lpoly.plg[4][2] + p[112] * lpoly.plg[6][2])*cd14)*
				cos(dgtr*llong)
				+ (p[91] * lpoly.plg[3][2] + p[92] * lpoly.plg[5][2] + p[93] * lpoly.plg[7][2]
				+ p[107] * lpoly.plg[2][2] + p[108] * lpoly.plg[4][2] + p[109] * lpoly.plg[6][2]
				+ csw.swc[5] * (p[113] * lpoly.plg[2][2] + p[114] * lpoly.plg[4][2] + p[115] * lpoly.plg[6][2])*cd14)*
				sin(dgtr*llong));
		}

		// ----  ut and mixed ut,longitude
		if (csw.sw[12] != 0)
		{
			t[12] = (1.0 + p[96] * lpoly.plg[2][1])*(1.0 + p[82] * lpoly.dfa*csw.swc[1])*
				(1.0 + p[120] * lpoly.plg[2][1] * csw.swc[5] * cd14)*
				((p[69] * lpoly.plg[2][1] + p[70] * lpoly.plg[4][1] + p[71] * lpoly.plg[6][1])*
				cos(sr*(sec - p[72])));
			t[12] = t[12] + csw.swc[11] *
				(p[77] * lpoly.plg[4][3] + p[78] * lpoly.plg[6][3] + p[79] * lpoly.plg[8][3])*
				cos(sr*(sec - p[80]) + 2.0*dgtr*llong)*(1.0 + p[138] * lpoly.dfa*
				csw.swc[1]);
		}

		// ----  ut,longitude magnetic activity
		if (csw.sw[13] == 0) goto fourtyeight;

		if (sw9 != -1.0)
		{
			t[13] = lpoly.apdf*csw.swc[11] * (1.0 + p[121] * lpoly.plg[2][1])*
				((p[61] * lpoly.plg[3][2] + p[62] * lpoly.plg[5][2] + p[63] * lpoly.plg[7][2])*
				cos(dgtr*(llong - p[64])))
				+ lpoly.apdf*csw.swc[11] * csw.swc[5] *
				(p[116] * lpoly.plg[2][2] + p[117] * lpoly.plg[4][2] + p[118] * lpoly.plg[6][2])*
				cd14*cos(dgtr*(llong - p[119]))
				+ lpoly.apdf*csw.swc[12] *
				(p[84] * lpoly.plg[2][1] + p[85] * lpoly.plg[4][1] + p[86] * lpoly.plg[6][1])*
				cos(sr*(sec - p[76]));
			goto fourtyeight;
		}

		if (p[52] == 0) goto fourtyeight;
		t[13] = lpoly.apt[1] * csw.swc[11] * (1.0 + p[133] * lpoly.plg[2][1])*
			((p[53] * lpoly.plg[3][2] + p[99] * lpoly.plg[5][2] + p[68] * lpoly.plg[7][2])*
			cos(dgtr*(llong - p[98])))
			+ lpoly.apt[1] * csw.swc[11] * csw.swc[5] *
			(p[134] * lpoly.plg[2][2] + p[135] * lpoly.plg[4][2] + p[136] * lpoly.plg[6][2])*
			cd14*cos(dgtr*(llong - p[137]))
			+ lpoly.apt[1] * csw.swc[12] *
			(p[56] * lpoly.plg[2][1] + p[57] * lpoly.plg[4][1] + p[58] * lpoly.plg[6][1])*
			cos(sr*(sec - p[59]));
	fourtyeight:
		//       parms not used: 83, 90,100,140-150
	fourtynine :
		tinfg = p[31];
			   for (i = 1; i <= nsw; i++)
				   tinfg = tinfg + fabsf(csw.sw[i])*t[i];

			   return tinfg;

	}
	/* ----------------------------------------------------------------------
	*      version of globe for lower atmosphere 10/26/99
	* ----------------------------------------------------------------------      */

	double glob7s
		(
		lpolytype& lpoly,
		cswtype& csw,
		double p[101]
		)
	{
		double cd32, cd18, cd14, cd39, t[15], t71, t72, t81, t82, tt;

		double dr = 1.72142e-2;
		double dgtr = 1.74533e-2;
		double pset = 2.0;
		double dayl = -1.0;
		double p32 = -1000.0;
		double p18 = -1000.0;
		double p14 = -1000.0;
		double p39 = -1000.0;

		int i;
		//------------------------------ begin ---------------------------------
		// ---- confirm parameter set
		if (p[100] == 0)
			p[100] = pset;
		if (p[100] != pset)
		{
			printf("wrong parameter set for glob7s %10.1f %10.1f \n", pset, p[100]);
			return 0;
		}

		for (i = 1; i <= 14; i++)
			t[i] = 0.0;

		if ((lpoly.day != dayl) | (p32 != p[32]))
			cd32 = cos(dr*(lpoly.day - p[32]));
		if ((lpoly.day != dayl) | (p18 != p[18]))
			cd18 = cos(2.0*dr*(lpoly.day - p[18]));

		if ((lpoly.day != dayl) | (p14 != p[14]))
			cd14 = cos(dr*(lpoly.day - p[14]));
		if ((lpoly.day != dayl) | (p39 != p[39]))
			cd39 = cos(2.0*dr*(lpoly.day - p[39]));

		dayl = lpoly.day;
		p32 = p[32];
		p18 = p[18];
		p14 = p[14];
		p39 = p[39];

		// ---- f10.7
		t[1] = p[22] * lpoly.dfa;

		// ---- time independent
		t[2] = p[2] * lpoly.plg[3][1] + p[3] * lpoly.plg[5][1] + p[23] * lpoly.plg[7][1]
			+ p[27] * lpoly.plg[2][1] + p[15] * lpoly.plg[4][1] + p[60] * lpoly.plg[6][1];

		// ---- symmetrical annual
		t[3] = (p[19] + p[48] * lpoly.plg[3][1] + p[30] * lpoly.plg[5][1])*cd32;

		// ---- symmetrical semiannual
		t[4] = (p[16] + p[17] * lpoly.plg[3][1] + p[31] * lpoly.plg[5][1])*cd18;

		// ---- asymmetrical annual
		t[5] = (p[10] * lpoly.plg[2][1] + p[11] * lpoly.plg[4][1] + p[21] * lpoly.plg[6][1])*cd14;

		// ---- asymmetrical semiannual
		t[6] = (p[38] * lpoly.plg[2][1])*cd39;

		// ----  diurnal
		if (csw.sw[7] != 0)
		{
			t71 = p[12] * lpoly.plg[3][2] * cd14*csw.swc[5];
			t72 = p[13] * lpoly.plg[3][2] * cd14*csw.swc[5];
			t[7] = ((p[4] * lpoly.plg[2][2] + p[5] * lpoly.plg[4][2] + t71)*lpoly.ctloc
				+ (p[7] * lpoly.plg[2][2] + p[8] * lpoly.plg[4][2] + t72)*lpoly.stloc);
		}

		// ----  semidiurnal
		if (csw.sw[8] != 0)
		{
			t81 = (p[24] * lpoly.plg[4][3] + p[36] * lpoly.plg[6][3])*cd14*csw.swc[5];
			t82 = (p[34] * lpoly.plg[4][3] + p[37] * lpoly.plg[6][3])*cd14*csw.swc[5];
			t[8] = ((p[6] * lpoly.plg[3][3] + p[42] * lpoly.plg[5][3] + t81)*lpoly.c2tloc
				+ (p[9] * lpoly.plg[3][3] + p[43] * lpoly.plg[5][3] + t82)*lpoly.s2tloc);
		}

		// ----  terdiurnal
		if (csw.sw[14] != 0)
			t[14] = p[40] * lpoly.plg[4][4] * lpoly.s3tloc + p[41] * lpoly.plg[4][4] * lpoly.c3tloc;

		// ---- magnetic activity
		if (csw.sw[9] != 0)
		{
			if (csw.sw[9] == 1)
				t[9] = lpoly.apdf*(p[33] + p[46] * lpoly.plg[3][1] * csw.swc[2]);
			if (csw.sw[9] == -1)
				t[9] = (p[51] * lpoly.apt[1] + p[97] * lpoly.plg[3][1] * lpoly.apt[1] * csw.swc[2]);
		}
		if ((csw.sw[10] == 0) | (csw.sw[11] == 0) | (lpoly.xlong <= -1000.0))
			goto fourtynine;

		// ----  longitudinal
		t[11] = (1.0 + lpoly.plg[2][1] * (p[81] * csw.swc[5] * cos(dr*(lpoly.day - p[82]))
			+ p[86] * csw.swc[6] * cos(2.0*dr*(lpoly.day - p[87])))
			+ p[84] * csw.swc[3] * cos(dr*(lpoly.day - p[85]))
			+ p[88] * csw.swc[4] * cos(2.0*dr*(lpoly.day - p[89])))
			*((p[65] * lpoly.plg[3][2] + p[66] * lpoly.plg[5][2] + p[67] * lpoly.plg[7][2]
			+ p[75] * lpoly.plg[2][2] + p[76] * lpoly.plg[4][2] + p[77] * lpoly.plg[6][2])*cos(dgtr*lpoly.xlong)
			+ (p[91] * lpoly.plg[3][2] + p[92] * lpoly.plg[5][2] + p[93] * lpoly.plg[7][2]
			+ p[78] * lpoly.plg[2][2] + p[79] * lpoly.plg[4][2] + p[80] * lpoly.plg[6][2])*sin(dgtr*lpoly.xlong));
	fourtynine:
		tt = 0.0;
		for (i = 1; i <= 14; i++)
			tt = tt + fabsf(csw.sw[i])*t[i];

		return tt;
	}

	/* ----------------------------------------------------------------------
	* nrlmsise-00 01-feb-02 dist 17
	*
	* this function loads all the iniital data for the msis00 model
	*
	* ----------------------------------------------------------------------      */

	void msis00init
		(
		msistype& msis00r
		)
	{

		int i, j;

		// cdav i'm not string literate yet in c++ :-)
		//        msis00r.datime.isdate = "01-feb-02";
		//        msis00r.datime.istime = "15:49:27 ";
		//        msis00r.datime.name   = "msise-00 ";

		msis00r.csw.isw = 0;
		//
		// cdav   set this to output in
		msis00r.metsel.imr = 1;

		// ----   temperature
		static double pt[151] =
		{ 0,
		9.86573e-01, 1.62228e-02, 1.55270e-02, -1.04323e-01, -3.75801e-03,
		-1.18538e-03, -1.24043e-01, 4.56820e-03, 8.76018e-03, -1.36235e-01,
		-3.52427e-02, 8.84181e-03, -5.92127e-03, -8.61650e+00, 0.00000e+00,
		1.28492e-02, 0.00000e+00, 1.30096e+02, 1.04567e-02, 1.65686e-03,
		-5.53887e-06, 2.97810e-03, 0.00000e+00, 5.13122e-03, 8.66784e-02,
		1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, -7.27026e-06,
		0.00000e+00, 6.74494e+00, 4.93933e-03, 2.21656e-03, 2.50802e-03,
		0.00000e+00, 0.00000e+00, -2.08841e-02, -1.79873e+00, 1.45103e-03,
		2.81769e-04, -1.44703e-03, -5.16394e-05, 8.47001e-02, 1.70147e-01,
		5.72562e-03, 5.07493e-05, 4.36148e-03, 1.17863e-04, 4.74364e-03,
		6.61278e-03, 4.34292e-05, 1.44373e-03, 2.41470e-05, 2.84426e-03,
		8.56560e-04, 2.04028e-03, 0.00000e+00, -3.15994e+03, -2.46423e-03,
		1.13843e-03, 4.20512e-04, 0.00000e+00, -9.77214e+01, 6.77794e-03,
		5.27499e-03, 1.14936e-03, 0.00000e+00, -6.61311e-03, -1.84255e-02,
		-1.96259e-02, 2.98618e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		6.44574e+02, 8.84668e-04, 5.05066e-04, 0.00000e+00, 4.02881e+03,
		-1.89503e-03, 0.00000e+00, 0.00000e+00, 8.21407e-04, 2.06780e-03,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		-1.20410e-02, -3.63963e-03, 9.92070e-05, -1.15284e-04, -6.33059e-05,
		-6.05545e-01, 8.34218e-03, -9.13036e+01, 3.71042e-04, 0.00000e+00,
		4.19000e-04, 2.70928e-03, 3.31507e-03, -4.44508e-03, -4.96334e-03,
		-1.60449e-03, 3.95119e-03, 2.48924e-03, 5.09815e-04, 4.05302e-03,
		2.24076e-03, 0.00000e+00, 6.84256e-03, 4.66354e-04, 0.00000e+00,
		-3.68328e-04, 0.00000e+00, 0.00000e+00, -1.46870e+02, 0.00000e+00,
		0.00000e+00, 1.09501e-03, 4.65156e-04, 5.62583e-04, 3.21596e+00,
		6.43168e-04, 3.14860e-03, 3.40738e-03, 1.78481e-03, 9.62532e-04,
		5.58171e-04, 3.43731e+00, -2.33195e-01, 5.10289e-04, 0.00000e+00,
		0.00000e+00, -9.25347e+04, 0.00000e+00, -1.99639e-03, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		};
		// ----   he density
		static double pd[10][151] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			},
			{ 0,
			1.09979e+00, -4.88060e-02, -1.97501e-01, -9.10280e-02, -6.96558e-03,
			2.42136e-02, 3.91333e-01, -7.20068e-03, -3.22718e-02, 1.41508e+00,
			1.68194e-01, 1.85282e-02, 1.09384e-01, -7.24282e+00, 0.00000e+00,
			2.96377e-01, -4.97210e-02, 1.04114e+02, -8.61108e-02, -7.29177e-04,
			1.48998e-06, 1.08629e-03, 0.00000e+00, 0.00000e+00, 8.31090e-02,
			1.12818e-01, -5.75005e-02, -1.29919e-02, -1.78849e-02, -2.86343e-06,
			0.00000e+00, -1.51187e+02, -6.65902e-03, 0.00000e+00, -2.02069e-03,
			0.00000e+00, 0.00000e+00, 4.32264e-02, -2.80444e+01, -3.26789e-03,
			2.47461e-03, 0.00000e+00, 0.00000e+00, 9.82100e-02, 1.22714e-01,
			-3.96450e-02, 0.00000e+00, -2.76489e-03, 0.00000e+00, 1.87723e-03,
			-8.09813e-03, 4.34428e-05, -7.70932e-03, 0.00000e+00, -2.28894e-03,
			-5.69070e-03, -5.22193e-03, 6.00692e-03, -7.80434e+03, -3.48336e-03,
			-6.38362e-03, -1.82190e-03, 0.00000e+00, -7.58976e+01, -2.17875e-02,
			-1.72524e-02, -9.06287e-03, 0.00000e+00, 2.44725e-02, 8.66040e-02,
			1.05712e-01, 3.02543e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-6.01364e+03, -5.64668e-03, -2.54157e-03, 0.00000e+00, 3.15611e+02,
			-5.69158e-03, 0.00000e+00, 0.00000e+00, -4.47216e-03, -4.49523e-03,
			4.64428e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			4.51236e-02, 2.46520e-02, 6.17794e-03, 0.00000e+00, 0.00000e+00,
			-3.62944e-01, -4.80022e-02, -7.57230e+01, -1.99656e-03, 0.00000e+00,
			-5.18780e-03, -1.73990e-02, -9.03485e-03, 7.48465e-03, 1.53267e-02,
			1.06296e-02, 1.18655e-02, 2.55569e-03, 1.69020e-03, 3.51936e-02,
			-1.81242e-02, 0.00000e+00, -1.00529e-01, -5.10574e-03, 0.00000e+00,
			2.10228e-03, 0.00000e+00, 0.00000e+00, -1.73255e+02, 5.07833e-01,
			-2.41408e-01, 8.75414e-03, 2.77527e-03, -8.90353e-05, -5.25148e+00,
			-5.83899e-03, -2.09122e-02, -9.63530e-03, 9.77164e-03, 4.07051e-03,
			2.53555e-04, -5.52875e+00, -3.55993e-01, -2.49231e-03, 0.00000e+00,
			0.00000e+00, 2.86026e+01, 0.00000e+00, 3.42722e-04, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   o density
			{ 0,
			1.02315e+00, -1.59710e-01, -1.06630e-01, -1.77074e-02, -4.42726e-03,
			3.44803e-02, 4.45613e-02, -3.33751e-02, -5.73598e-02, 3.50360e-01,
			6.33053e-02, 2.16221e-02, 5.42577e-02, -5.74193e+00, 0.00000e+00,
			1.90891e-01, -1.39194e-02, 1.01102e+02, 8.16363e-02, 1.33717e-04,
			6.54403e-06, 3.10295e-03, 0.00000e+00, 0.00000e+00, 5.38205e-02,
			1.23910e-01, -1.39831e-02, 0.00000e+00, 0.00000e+00, -3.95915e-06,
			0.00000e+00, -7.14651e-01, -5.01027e-03, 0.00000e+00, -3.24756e-03,
			0.00000e+00, 0.00000e+00, 4.42173e-02, -1.31598e+01, -3.15626e-03,
			1.24574e-03, -1.47626e-03, -1.55461e-03, 6.40682e-02, 1.34898e-01,
			-2.42415e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 6.13666e-04,
			-5.40373e-03, 2.61635e-05, -3.33012e-03, 0.00000e+00, -3.08101e-03,
			-2.42679e-03, -3.36086e-03, 0.00000e+00, -1.18979e+03, -5.04738e-02,
			-2.61547e-03, -1.03132e-03, 1.91583e-04, -8.38132e+01, -1.40517e-02,
			-1.14167e-02, -4.08012e-03, 1.73522e-04, -1.39644e-02, -6.64128e-02,
			-6.85152e-02, -1.34414e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			6.07916e+02, -4.12220e-03, -2.20996e-03, 0.00000e+00, 1.70277e+03,
			-4.63015e-03, 0.00000e+00, 0.00000e+00, -2.25360e-03, -2.96204e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			3.92786e-02, 1.31186e-02, -1.78086e-03, 0.00000e+00, 0.00000e+00,
			-3.90083e-01, -2.84741e-02, -7.78400e+01, -1.02601e-03, 0.00000e+00,
			-7.26485e-04, -5.42181e-03, -5.59305e-03, 1.22825e-02, 1.23868e-02,
			6.68835e-03, -1.03303e-02, -9.51903e-03, 2.70021e-04, -2.57084e-02,
			-1.32430e-02, 0.00000e+00, -3.81000e-02, -3.16810e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -9.05762e-04, -2.14590e-03, -1.17824e-03, 3.66732e+00,
			-3.79729e-04, -6.13966e-03, -5.09082e-03, -1.96332e-03, -3.08280e-03,
			-9.75222e-04, 4.03315e+00, -2.52710e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   n2 density
			{ 0,
			1.16112e+00, 0.00000e+00, 0.00000e+00, 3.33725e-02, 0.00000e+00,
			3.48637e-02, -5.44368e-03, 0.00000e+00, -6.73940e-02, 1.74754e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.74712e+02, 0.00000e+00,
			1.26733e-01, 0.00000e+00, 1.03154e+02, 5.52075e-02, 0.00000e+00,
			0.00000e+00, 8.13525e-04, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -2.50482e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.48894e-03,
			6.16053e-04, -5.79716e-04, 2.95482e-03, 8.47001e-02, 1.70147e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.47425e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   msis00r.gts3c.tlb
			{ 0,
			9.44846e-01, 0.00000e+00, 0.00000e+00, -3.08617e-02, 0.00000e+00,
			-2.44019e-02, 6.48607e-03, 0.00000e+00, 3.08181e-02, 4.59392e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.74712e+02, 0.00000e+00,
			2.13260e-02, 0.00000e+00, -3.56958e+02, 0.00000e+00, 1.82278e-04,
			0.00000e+00, 3.07472e-04, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 3.83054e-03, 0.00000e+00, 0.00000e+00,
			-1.93065e-03, -1.45090e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.23493e-03, 1.36736e-03, 8.47001e-02, 1.70147e-01,
			3.71469e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			5.10250e-03, 2.47425e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 3.68756e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   o2 density
			{ 0,
			1.35580e+00, 1.44816e-01, 0.00000e+00, 6.07767e-02, 0.00000e+00,
			2.94777e-02, 7.46900e-02, 0.00000e+00, -9.23822e-02, 8.57342e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 2.38636e+01, 0.00000e+00,
			7.71653e-02, 0.00000e+00, 8.18751e+01, 1.87736e-02, 0.00000e+00,
			0.00000e+00, 1.49667e-02, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -3.67874e+02, 5.48158e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
			1.22631e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			8.17187e-03, 3.71617e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.10826e-03,
			-3.13640e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-7.35742e-02, -5.00266e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.94965e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----   ar density
			{ 0,
			1.04761e+00, 2.00165e-01, 2.37697e-01, 3.68552e-02, 0.00000e+00,
			3.57202e-02, -2.14075e-01, 0.00000e+00, -1.08018e-01, -3.73981e-01,
			0.00000e+00, 3.10022e-02, -1.16305e-03, -2.07596e+01, 0.00000e+00,
			8.64502e-02, 0.00000e+00, 9.74908e+01, 5.16707e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 3.46193e+02, 1.34297e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.48509e-03,
			-1.54689e-04, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
			1.47753e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.89320e-02, 3.68181e-05, 1.32570e-02, 0.00000e+00, 0.00000e+00,
			3.59719e-03, 7.44328e-03, -1.00023e-03, -6.50528e+03, 0.00000e+00,
			1.03485e-02, -1.00983e-03, -4.06916e-03, -6.60864e+01, -1.71533e-02,
			1.10605e-02, 1.20300e-02, -5.20034e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-2.62769e+03, 7.13755e-03, 4.17999e-03, 0.00000e+00, 1.25910e+04,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -2.23595e-03, 4.60217e-03,
			5.71794e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-3.18353e-02, -2.35526e-02, -1.36189e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.03522e-02, -6.67837e+01, -1.09724e-03, 0.00000e+00,
			-1.38821e-02, 1.60468e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.51574e-02,
			-5.44470e-04, 0.00000e+00, 7.28224e-02, 6.59413e-02, 0.00000e+00,
			-5.15692e-03, 0.00000e+00, 0.00000e+00, -3.70367e+03, 0.00000e+00,
			0.00000e+00, 1.36131e-02, 5.38153e-03, 0.00000e+00, 4.76285e+00,
			-1.75677e-02, 2.26301e-02, 0.00000e+00, 1.76631e-02, 4.77162e-03,
			0.00000e+00, 5.39354e+00, 0.00000e+00, -7.51710e-03, 0.00000e+00,
			0.00000e+00, -8.82736e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    h density
			{ 0,
			1.26376e+00, -2.14304e-01, -1.49984e-01, 2.30404e-01, 2.98237e-02,
			2.68673e-02, 2.96228e-01, 2.21900e-02, -2.07655e-02, 4.52506e-01,
			1.20105e-01, 3.24420e-02, 4.24816e-02, -9.14313e+00, 0.00000e+00,
			2.47178e-02, -2.88229e-02, 8.12805e+01, 5.10380e-02, -5.80611e-03,
			2.51236e-05, -1.24083e-02, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, -3.48190e-02, 0.00000e+00, 0.00000e+00, 2.89885e-05,
			0.00000e+00, 1.53595e+02, -1.68604e-02, 0.00000e+00, 1.01015e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.84552e-04,
			-1.22181e-03, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
			-1.04927e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.91313e-03,
			-2.30501e-02, 3.14758e-05, 0.00000e+00, 0.00000e+00, 1.26956e-02,
			8.35489e-03, 3.10513e-04, 0.00000e+00, 3.42119e+03, -2.45017e-03,
			-4.27154e-04, 5.45152e-04, 1.89896e-03, 2.89121e+01, -6.49973e-03,
			-1.93855e-02, -1.48492e-02, 0.00000e+00, -5.10576e-02, 7.87306e-02,
			9.51981e-02, -1.49422e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			2.65503e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 6.37110e-03, 3.24789e-04,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			6.14274e-02, 1.00376e-02, -8.41083e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.27099e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-3.94077e-03, -1.28601e-02, -7.97616e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -6.71465e-03, -1.69799e-03, 1.93772e-03, 3.81140e+00,
			-7.79290e-03, -1.82589e-02, -1.25860e-02, -1.04311e-02, -3.02465e-03,
			2.43063e-03, 3.63237e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----    n density
			{ 0,
			7.09557e+01, -3.26740e-01, 0.00000e+00, -5.16829e-01, -1.71664e-03,
			9.09310e-02, -6.71500e-01, -1.47771e-01, -9.27471e-02, -2.30862e-01,
			-1.56410e-01, 1.34455e-02, -1.19717e-01, 2.52151e+00, 0.00000e+00,
			-2.41582e-01, 5.92939e-02, 4.39756e+00, 9.15280e-02, 4.41292e-03,
			0.00000e+00, 8.66807e-03, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 9.74701e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 6.70217e+01, -1.31660e-03, 0.00000e+00, -1.65317e-02,
			0.00000e+00, 0.00000e+00, 8.50247e-02, 2.77428e+01, 4.98658e-03,
			6.15115e-03, 9.50156e-03, -2.12723e-02, 8.47001e-02, 1.70147e-01,
			-2.38645e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.37380e-03,
			-8.41918e-03, 2.80145e-05, 7.12383e-03, 0.00000e+00, -1.66209e-02,
			1.03533e-04, -1.68898e-02, 0.00000e+00, 3.64526e+03, 0.00000e+00,
			6.54077e-03, 3.69130e-04, 9.94419e-04, 8.42803e+01, -1.16124e-02,
			-7.74414e-03, -1.68844e-03, 1.42809e-03, -1.92955e-03, 1.17225e-01,
			-2.41512e-02, 1.50521e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			1.60261e+03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -3.54403e-04, -1.87270e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			2.76439e-02, 6.43207e-03, -3.54300e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -2.80221e-02, 8.11228e+01, -6.75255e-04, 0.00000e+00,
			-1.05162e-02, -3.48292e-03, -6.97321e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.45546e-03, -1.31970e-02, -3.57751e-03, -1.09021e+00,
			-1.50181e-02, -7.12841e-03, -6.64590e-03, -3.52610e-03, -1.87773e-02,
			-2.22432e-03, -3.93895e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			// ----  hot o density
			{ 0,
			6.04050e-02, 1.57034e+00, 2.99387e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.51018e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -8.61650e+00, 1.26454e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 5.50878e-03, 0.00000e+00, 0.00000e+00, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 6.23881e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 8.47001e-02, 1.70147e-01,
			-9.45934e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			}
		};
		// ----    msis00r.gts3c.ps param
		static double ps[151] =
		{ 0,
		9.56827e-01, 6.20637e-02, 3.18433e-02, 0.00000e+00, 0.00000e+00,
		3.94900e-02, 0.00000e+00, 0.00000e+00, -9.24882e-03, -7.94023e-03,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 1.74712e+02, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 2.74677e-03, 0.00000e+00, 1.54951e-02, 8.66784e-02,
		1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, -6.99007e-04, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 1.24362e-02, -5.28756e-03, 8.47001e-02, 1.70147e-01,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 2.47425e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		};
		// ----    turbo
		static double pdl[3][26] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0
			},
			{ 0,
			1.09930e+00, 3.90631e+00, 3.07165e+00, 9.86161e-01, 1.63536e+01,
			4.63830e+00, 1.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.28840e+00, 3.10302e-02, 1.18339e-01
			},
			{ 0,
			1.00000e+00, 7.00000e-01, 1.15020e+00, 3.44689e+00, 1.28840e+00,
			1.00000e+00, 1.08738e+00, 1.22947e+00, 1.10016e+00, 7.34129e-01,
			1.15241e+00, 2.22784e+00, 7.95046e-01, 4.01612e+00, 4.47749e+00,
			1.23435e+02, -7.60535e-02, 1.68986e-06, 7.44294e-01, 1.03604e+00,
			1.72783e+02, 1.15020e+00, 3.44689e+00, -7.46230e-01, 9.49154e-01
			}
		};
		// ----   lower bouneary
		static double ptm[11] =
		{ 0,
		1.04130e+03, 3.86000e+02, 1.95000e+02, 1.66728e+01, 2.13000e+02,
		1.20000e+02, 2.40000e+02, 1.87000e+02, -2.00000e+00, 0.00000e+00
		};
		static double pdm[9][11] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			},
			{ 0,
			2.45600e+07, 6.71072e-06, 1.00000e+02, 0.00000e+00, 1.10000e+02,
			1.00000e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			8.59400e+10, 1.00000e+00, 1.05000e+02, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			2.81000e+11, 0.00000e+00, 1.05000e+02, 2.80000e+01, 2.89500e+01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			3.30000e+10, 2.68270e-01, 1.05000e+02, 1.00000e+00, 1.10000e+02,
			1.00000e+01, 1.10000e+02, -1.00000e+01, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			1.33000e+09, 1.19615e-02, 1.05000e+02, 0.00000e+00, 1.10000e+02,
			1.00000e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			1.76100e+05, 1.00000e+00, 9.50000e+01, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			1.00000e+07, 1.00000e+00, 1.05000e+02, -8.00000e+00, 1.10000e+02,
			1.00000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 0.00000e+00
			},
			{ 0,
			1.00000e+06, 1.00000e+00, 1.05000e+02, -8.00000e+00, 5.50000e+02,
			7.60000e+01, 9.00000e+01, 2.00000e+00, 0.00000e+00, 4.00000e+03
			}
		};
		// ----   msis00r.meso.tn1[2]
		static double ptl[5][101] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			},
			{ 0,
			1.00858e+00, 4.56011e-02, -2.22972e-02, -5.44388e-02, 5.23136e-04,
			-1.88849e-02, 5.23707e-02, -9.43646e-03, 6.31707e-03, -7.80460e-02,
			-4.88430e-02, 0.00000e+00, 0.00000e+00, -7.60250e+00, 0.00000e+00,
			-1.44635e-02, -1.76843e-02, -1.21517e+02, 2.85647e-02, 0.00000e+00,
			0.00000e+00, 6.31792e-04, 0.00000e+00, 5.77197e-03, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -8.90272e+03, 3.30611e-03, 3.02172e-03, 0.00000e+00,
			-2.13673e-03, -3.20910e-04, 0.00000e+00, 0.00000e+00, 2.76034e-03,
			2.82487e-03, -2.97592e-04, -4.21534e-03, 8.47001e-02, 1.70147e-01,
			8.96456e-03, 0.00000e+00, -1.08596e-02, 0.00000e+00, 0.00000e+00,
			5.57917e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 9.65405e-03, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----   msis00r.meso.tn1[3]
			{ 0,
			9.39664e-01, 8.56514e-02, -6.79989e-03, 2.65929e-02, -4.74283e-03,
			1.21855e-02, -2.14905e-02, 6.49651e-03, -2.05477e-02, -4.24952e-02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 1.19148e+01, 0.00000e+00,
			1.18777e-02, -7.28230e-02, -8.15965e+01, 1.73887e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.44691e-02, 2.80259e-04, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.16584e+02, 3.18713e-03, 7.37479e-03, 0.00000e+00,
			-2.55018e-03, -3.92806e-03, 0.00000e+00, 0.00000e+00, -2.89757e-03,
			-1.33549e-03, 1.02661e-03, 3.53775e-04, 8.47001e-02, 1.70147e-01,
			-9.17497e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			3.56082e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.00902e-02, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----   msis00r.meso.tn1[4]
			{ 0,
			9.85982e-01, -4.55435e-02, 1.21106e-02, 2.04127e-02, -2.40836e-03,
			1.11383e-02, -4.51926e-02, 1.35074e-02, -6.54139e-03, 1.15275e-01,
			1.28247e-01, 0.00000e+00, 0.00000e+00, -5.30705e+00, 0.00000e+00,
			-3.79332e-02, -6.24741e-02, 7.71062e-01, 2.96315e-02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 6.81051e-03, -4.34767e-03, 8.66784e-02,
			1.58727e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.07003e+01, -2.76907e-03, 4.32474e-04, 0.00000e+00,
			1.31497e-03, -6.47517e-04, 0.00000e+00, -2.20621e+01, -1.10804e-03,
			-8.09338e-04, 4.18184e-04, 4.29650e-03, 8.47001e-02, 1.70147e-01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-4.04337e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -9.52550e-04,
			8.56253e-04, 4.33114e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.21223e-03,
			2.38694e-04, 9.15245e-04, 1.28385e-03, 8.67668e-04, -5.61425e-06,
			1.04445e+00, 3.41112e+01, 0.00000e+00, -8.40704e-01, -2.39639e+02,
			7.06668e-01, -2.05873e+01, -3.63696e-01, 2.39245e+01, 0.00000e+00,
			-1.06657e-03, -7.67292e-04, 1.54534e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----   msis00r.meso.tn1[5] msis00r.meso.tn2[1]
			{ 0,
			1.00320e+00, 3.83501e-02, -2.38983e-03, 2.83950e-03, 4.20956e-03,
			5.86619e-04, 2.19054e-02, -1.00946e-02, -3.50259e-03, 4.17392e-02,
			-8.44404e-03, 0.00000e+00, 0.00000e+00, 4.96949e+00, 0.00000e+00,
			-7.06478e-03, -1.46494e-02, 3.13258e+01, -1.86493e-03, 0.00000e+00,
			-1.67499e-02, 0.00000e+00, 0.00000e+00, 5.12686e-04, 8.66784e-02,
			1.58727e-01, -4.64167e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			4.37353e-03, -1.99069e+02, 0.00000e+00, -5.34884e-03, 0.00000e+00,
			1.62458e-03, 2.93016e-03, 2.67926e-03, 5.90449e+02, 0.00000e+00,
			0.00000e+00, -1.17266e-03, -3.58890e-04, 8.47001e-02, 1.70147e-01,
			0.00000e+00, 0.00000e+00, 1.38673e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.60571e-03,
			6.28078e-04, 5.05469e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.57829e-03,
			-4.00855e-04, 5.04077e-05, -1.39001e-03, -2.33406e-03, -4.81197e-04,
			1.46758e+00, 6.20332e+00, 0.00000e+00, 3.66476e-01, -6.19760e+01,
			3.09198e-01, -1.98999e+01, 0.00000e+00, -3.29933e+02, 0.00000e+00,
			-1.10080e-03, -9.39310e-05, 1.39638e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			}
		};
		// ----    msis00r.meso.tn2[2]
		static double pma[11][101] = {
			{ 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			},
			{ 0,
			9.81637e-01, -1.41317e-03, 3.87323e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.58707e-02,
			-8.63658e-03, 0.00000e+00, 0.00000e+00, -2.02226e+00, 0.00000e+00,
			-8.69424e-03, -1.91397e-02, 8.76779e+01, 4.52188e-03, 0.00000e+00,
			2.23760e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -7.07572e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			-4.11210e-03, 3.50060e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -8.36657e-03, 1.61347e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.45130e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.24152e-03,
			6.43365e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.33255e-03,
			2.42657e-03, 1.60666e-03, -1.85728e-03, -1.46874e-03, -4.79163e-06,
			1.22464e+00, 3.53510e+01, 0.00000e+00, 4.49223e-01, -4.77466e+01,
			4.70681e-01, 8.41861e+00, -2.88198e-01, 1.67854e+02, 0.00000e+00,
			7.11493e-04, 6.05601e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----    msis00r.meso.tn2[3]
			{ 0,
			1.00422e+00, -7.11212e-03, 5.24480e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -5.28914e-02,
			-2.41301e-02, 0.00000e+00, 0.00000e+00, -2.12219e+01, -1.03830e-02,
			-3.28077e-03, 1.65727e-02, 1.68564e+00, -6.68154e-03, 0.00000e+00,
			1.45155e-02, 0.00000e+00, 8.42365e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -4.34645e-03, 0.00000e+00, 0.00000e+00, 2.16780e-02,
			0.00000e+00, -1.38459e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 7.04573e-03, -4.73204e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.08767e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -8.08279e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.21769e-04,
			-2.27387e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.26769e-03,
			3.16901e-03, 4.60316e-04, -1.01431e-04, 1.02131e-03, 9.96601e-04,
			1.25707e+00, 2.50114e+01, 0.00000e+00, 4.24472e-01, -2.77655e+01,
			3.44625e-01, 2.75412e+01, 0.00000e+00, 7.94251e+02, 0.00000e+00,
			2.45835e-03, 1.38871e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----    msis00r.meso.tn2[4] tn3[1]
			{ 0,
			1.01890e+00, -2.46603e-02, 1.00078e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -6.70977e-02,
			-4.02286e-02, 0.00000e+00, 0.00000e+00, -2.29466e+01, -7.47019e-03,
			2.26580e-03, 2.63931e-02, 3.72625e+01, -6.39041e-03, 0.00000e+00,
			9.58383e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.85291e-03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 1.39717e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 9.19771e-03, -3.69121e+02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.57067e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -7.07265e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.92953e-03,
			-2.77739e-03, -4.40092e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.47280e-03,
			2.95035e-04, -1.81246e-03, 2.81945e-03, 4.27296e-03, 9.78863e-04,
			1.40545e+00, -6.19173e+00, 0.00000e+00, 0.00000e+00, -7.93632e+01,
			4.44643e-01, -4.03085e+02, 0.00000e+00, 1.15603e+01, 0.00000e+00,
			2.25068e-03, 8.48557e-04, -2.98493e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----    tn3[2]
			{ 0,
			9.75801e-01, 3.80680e-02, -3.05198e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.85575e-02,
			5.04057e-02, 0.00000e+00, 0.00000e+00, -1.76046e+02, 1.44594e-02,
			-1.48297e-03, -3.68560e-03, 3.02185e+01, -3.23338e-03, 0.00000e+00,
			1.53569e-02, 0.00000e+00, -1.15558e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 4.89620e-03, 0.00000e+00, 0.00000e+00, -1.00616e-02,
			-8.21324e-03, -1.57757e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 6.63564e-03, 4.58410e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -2.51280e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 9.91215e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -8.73148e-04,
			-1.29648e-03, -7.32026e-05, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.68110e-03,
			-4.66003e-03, -1.31567e-03, -7.39390e-04, 6.32499e-04, -4.65588e-04,
			-1.29785e+00, -1.57139e+02, 0.00000e+00, 2.58350e-01, -3.69453e+01,
			4.10672e-01, 9.78196e+00, -1.52064e-01, -3.85084e+03, 0.00000e+00,
			-8.52706e-04, -1.40945e-03, -7.26786e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----    tn3[3]
			{ 0,
			9.60722e-01, 7.03757e-02, -3.00266e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.22671e-02,
			4.10423e-02, 0.00000e+00, 0.00000e+00, -1.63070e+02, 1.06073e-02,
			5.40747e-04, 7.79481e-03, 1.44908e+02, 1.51484e-04, 0.00000e+00,
			1.97547e-02, 0.00000e+00, -1.41844e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 5.77884e-03, 0.00000e+00, 0.00000e+00, 9.74319e-03,
			0.00000e+00, -2.88015e+03, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -4.44902e-03, -2.92760e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 2.34419e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.36685e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -4.65325e-04,
			-5.50628e-04, 3.31465e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.06179e-03,
			-3.08575e-03, -7.93589e-04, -1.08629e-04, 5.95511e-04, -9.05050e-04,
			1.18997e+00, 4.15924e+01, 0.00000e+00, -4.72064e-01, -9.47150e+02,
			3.98723e-01, 1.98304e+01, 0.00000e+00, 3.73219e+03, 0.00000e+00,
			-1.50040e-03, -1.14933e-03, -1.56769e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----    tn3[4]
			{ 0,
			1.03123e+00, -7.05124e-02, 8.71615e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -3.82621e-02,
			-9.80975e-03, 0.00000e+00, 0.00000e+00, 2.89286e+01, 9.57341e-03,
			0.00000e+00, 0.00000e+00, 8.66153e+01, 7.91938e-04, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 4.68917e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 7.86638e-03, 0.00000e+00, 0.00000e+00, 9.90827e-03,
			0.00000e+00, 6.55573e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -4.00200e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 7.07457e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.72268e-03,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.04970e-04,
			1.21560e-03, -8.05579e-06, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.49941e-03,
			-4.57256e-04, -1.59311e-04, 2.96481e-04, -1.77318e-03, -6.37918e-04,
			1.02395e+00, 1.28172e+01, 0.00000e+00, 1.49903e-01, -2.63818e+01,
			0.00000e+00, 4.70628e+01, -2.22139e-01, 4.82292e-02, 0.00000e+00,
			-8.67075e-04, -5.86479e-04, 5.32462e-04, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----    tn3[5] surface temp tsl
			{ 0,
			1.00828e+00, -9.10404e-02, -2.26549e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -2.32420e-02,
			-9.08925e-03, 0.00000e+00, 0.00000e+00, 3.36105e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.24957e+01, -5.87939e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.79765e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 2.01237e+03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -1.75553e-02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.29699e-03,
			1.26659e-03, 2.68402e-04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 1.17894e-03,
			1.48746e-03, 1.06478e-04, 1.34743e-04, -2.20939e-03, -6.23523e-04,
			6.36539e-01, 1.13621e+01, 0.00000e+00, -3.93777e-01, 2.38687e+03,
			0.00000e+00, 6.61865e+02, -1.21434e-01, 9.27608e+00, 0.00000e+00,
			1.68478e-04, 1.24892e-03, 1.71345e-03, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----    tgn3[2] surface grae tslg
			{ 0,
			1.57293e+00, -6.78400e-01, 6.47500e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -7.62974e-02,
			-3.60423e-01, 0.00000e+00, 0.00000e+00, 1.28358e+02, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 4.68038e+01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.67898e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.90994e+04, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 3.15706e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----    tgn2[1] tgn1[2]
			{ 0,
			8.60028e-01, 3.77052e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.17570e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 7.77757e-03, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 1.01024e+02, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 6.54251e+02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, -1.56959e-02,
			1.91001e-02, 3.15971e-02, 1.00982e-02, -6.71565e-03, 2.57693e-03,
			1.38692e+00, 2.82132e-01, 0.00000e+00, 0.00000e+00, 3.81511e+02,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			},
			// ----    tgn3[1] tgn2[2]
			{ 0,
			1.06029e+00, -5.25231e-02, 3.73034e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.31072e-02,
			-3.88409e-01, 0.00000e+00, 0.00000e+00, -1.65295e+02, -2.13801e-01,
			-4.38916e-02, -3.22716e-01, -8.82393e+01, 1.18458e-01, 0.00000e+00,
			-4.35863e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, -1.19782e-01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 2.62229e+01, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, -5.37443e+01, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, -4.55788e-01, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 3.84009e-02,
			3.96733e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 5.05494e-02,
			7.39617e-02, 1.92200e-02, -8.46151e-03, -1.34244e-02, 1.96338e-02,
			1.50421e+00, 1.88368e+01, 0.00000e+00, 0.00000e+00, -5.13114e+01,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			5.11923e-02, 3.61225e-02, 0.00000e+00, 0.00000e+00, 0.00000e+00,
			0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 2.00000e+00
			}
		};
		// ----    semiannual mult sam
		static double sam[101] =
		{ 0,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
		0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00
		};
		// ----   middle atmosphere averages
		static double pavgm[11] =
		{ 0,
		2.61000e+02, 2.64000e+02, 2.29000e+02, 2.17000e+02, 2.17000e+02,
		2.23000e+02, 2.86760e+02, -2.93940e+00, 2.50000e+00, 0.00000e+00
		};

		// now assign the record values
		memcpy(msis00r.parm.pt, pt, 151 * sizeof(double));
		memcpy(msis00r.parm.ps, ps, 151 * sizeof(double));
		memcpy(msis00r.lower.ptm, ptm, 11 * sizeof(double));
		memcpy(msis00r.parm.sam, sam, 101 * sizeof(double));
		memcpy(msis00r.mavg.pavgm, pavgm, 11 * sizeof(double));

		// since c reads in rows first, while the fortran code assumed columns,
		// we do a loop

		for (i = 0; i <= 8; i++)
		for (j = 0; j <= 10; j++)
			msis00r.lower.pdm[j][i] = pdm[i][j];
		for (i = 0; i <= 4; i++)
		for (j = 0; j <= 100; j++)
			msis00r.parm.ptl[j][i] = ptl[i][j];
		for (i = 0; i <= 9; i++)
		for (j = 0; j <= 150; j++)
			msis00r.parm.pd[j][i] = pd[i][j];

		for (i = 0; i <= 2; i++)
		for (j = 0; j <= 25; j++)
			msis00r.parm.pdl[j][i] = pdl[i][j];
		for (i = 0; i <= 10; i++)
		for (j = 0; j <= 100; j++)
			msis00r.parm.pma[j][i] = pma[i][j];

		// cdav
		// these are the extra arrays needed so a whole array can be accessed, but in
		// reverse order from the standard [row,col] approach in the fortran.
		//    extras to aid function calls for a subarry
		memcpy(msis00r.parm.pddav, pd, 10 * 151 * sizeof(double));
		memcpy(msis00r.parm.ptldav, ptl, 5 * 101 * sizeof(double));
		memcpy(msis00r.parm.pmadav, pma, 11 * 101 * sizeof(double));

	}

}  // namespace MSIS_Vers