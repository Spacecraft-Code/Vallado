/* -----------------------------------------------------------------------------
*
*                                astIOD.cpp
*
*   this file contains fundamental astrodynamic procedures and functions
*   relating to the initial orbit determination techniques. see ch 4, and ch 7
*   for a complete discussion of these routines.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2013
*                             by david vallado
*
*      (w) 719-573-2600, email dvallado@agi.com
*
*      *****************************************************************
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*              31 mar 08  david vallado
*                           misc updates
*              14 feb 08  david vallado
*                           fix razel conversions
*              15 mar 07  david vallado
*                           3rd edition baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
----------------------------------------------------------------------      */

#include "astIOD.h"

namespace astIOD {

/* utility functions for lambertbattin, etc */
static double k(double v);

static double see(double v);

/*---------------------------------------------------------------------------
*
*                           procedure site
*
*  this function finds the position and velocity vectors for a site.  the
*    answer is returned in the geocentric equatorial (ecef) coordinate system.
*    note that the velocity is zero because the coordinate system is fixed to
*    the earth.
*
*  author        : david vallado                  719-573-2600   25 jun 2002
*
*  inputs          description                    range / units
*    latgd       - geodetic latitude              -pi/2 to pi/2 rad
*    lon         - longitude of site              -2pi to 2pi rad
*    alt         - altitude                       km
*
*  outputs       :
*    rsecef      - ecef site position vector      km
*    vsecef      - ecef site velocity vector      km/s
*
*  locals        :
*    sinlat      - variable containing  sin(lat)  rad
*    temp        - temporary real value
*    rdel        - rdel component of site vector  km
*    rk          - rk component of site vector    km
*    cearth      -
*
*  coupling      :
*    none
*
*  references    :
*    vallado       2013, 430, alg 51, ex 7-1
----------------------------------------------------------------------------*/

void site
(
double latgd, double lon, double alt,
double rsecef[3], double vsecef[3]
)
{
	const double rearth = 6378.137;  // km
	const double eesqrd = 0.00669437999013;
	double   sinlat, cearth, rdel, rk;

	/* ---------------------  initialize values   ------------------- */
	sinlat = sin(latgd);

	/* -------  find rdel and rk components of site vector  --------- */
	cearth = rearth / sqrt(1.0 - (eesqrd * sinlat * sinlat));
	rdel = (cearth + alt) * cos(latgd);
	rk = ((1.0 - eesqrd) * cearth + alt) * sinlat;

	/* ----------------  find site position vector  ----------------- */
	rsecef[0] = rdel * cos(lon);
	rsecef[1] = rdel * sin(lon);
	rsecef[2] = rk;

	/* ----------------  find site velocity vector  ----------------- */
	vsecef[0] = 0.0;
	vsecef[1] = 0.0;
	vsecef[2] = 0.0;
}  // procedure site



/* -------------------------- angles-only techniques ------------------------ */


/*------------------------------------------------------------------------------
*
*                           procedure anglesgauss
*
*  this procedure solves the problem of orbit determination using three
*    optical sightings.  the solution procedure uses the gaussian technique.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    trtasc1      - right ascension #1            rad
*    trtasc2      - right ascension #2            rad
*    trtasc3      - right ascension #3            rad
*    tdecl1       - declination #1                rad
*    tdecl2       - declination #2                rad
*    tdecl3       - declination #3                rad
*    jd1          - julian date of 1st sighting   days from 4713 bc
*    jd2          - julian date of 2nd sighting   days from 4713 bc
*    jd3          - julian date of 3rd sighting   days from 4713 bc
*    rsijk        - ijk site position vector      er
*
*  outputs        :
*    r            - ijk position vector at t2     er
*    v            - ijk velocity vector at t2     er / tu
*
*  locals         :
*    l1           - line of sight vector for 1st
*    l2           - line of sight vector for 2nd
*    l3           - line of sight vector for 3rd
*    tau          - taylor expansion series about
*                   tau ( t - to )
*    tausqr       - tau squared
*    t21t23       - (t2-t1) * (t2-t3)
*    t31t32       - (t3-t1) * (t3-t2)
*    i            - index
*    d            -
*    rho          - range from site to sat at t2  er
*    rhodot       -
*    dmat         -
*    rs1          - site vectors
*    rs2          -
*    rs3          -
*    earthrate    - velocity of earth rotation
*    p            -
*    q            -
*    oldr         -
*    oldv         -
*    f1           - f coefficient
*    g1           -
*    f3           -
*    g3           -
*    l2dotrs      -
*
*  coupling       :
*    astMath::mag          - astMath::magnitude of a vector
*    detrminant   - evaluate the determinant of a matrix
*    factor       - find roots of a polynomial
*    matmult      - multiply two matrices together
*    assignval    - assign a value to a matrix
*    getval       - get a value from a matrix
*    initmatrix   - initialize a matrix and fil with 0.0's
*    delmatrix    - delete a matrix
*    gibbs        - gibbs method of orbit determination
*    hgibbs       - herrick gibbs method of orbit determination
*    angle        - angle between two vectors
*
*  references     :
*    vallado       2013, 442, alg 52, ex 7-2
-----------------------------------------------------------------------------*/

void anglesgauss
(
double tdecl1, double tdecl2, double tdecl3,
double trtasc1, double trtasc2, double trtasc3,
double jd1, double jd2, double jd3,
double rs1[3], double rs2[3], double rs3[3], double r2[3], double v2[3]
)
{
	const double tuday = 0.00933809017716;
	const double mu = 1.0;
	const double small = 0.0000001;
	const double rad = 180.0 / pi;
	int i, j, ll;
	char error[12];
	double poly[16];
	double roots[15][2];
	double r1[3], r3[3], l1[3], l2[3], l3[3];
//	std::vector< std::vector<double> > lmatii(3, 3), cmat(3, 3), rhomat(3, 3),
//		lmati(3, 3), rsmat(3, 3), lir(3, 3);
	std::vector< std::vector<double> > lmatii = std::vector< std::vector<double> >(3, std::vector<double>(3));
	std::vector< std::vector<double> > cmat = std::vector< std::vector<double> >(3, std::vector<double>(3));
	std::vector< std::vector<double> > rhomat = std::vector< std::vector<double> >(3, std::vector<double>(3));
	std::vector< std::vector<double> > lmati = std::vector< std::vector<double> >(3, std::vector<double>(3));
	std::vector< std::vector<double> > rsmat = std::vector< std::vector<double> >(3, std::vector<double>(3));
	std::vector< std::vector<double> > lir = std::vector< std::vector<double> >(3, std::vector<double>(3));

	double rdot, tau1, tau3, u, udot, p, f1, g1, f3, g3, a, ecc, incl, omega, argp,
		nu, m, l, argper, bigr2, a1, a1u, a3, a3u, d, d1, d2, c1, c3, l2dotrs,
		rhoold1, rhoold2, rhoold3, temproot, theta, theta1, copa, tausqr;

	/* ----------------------   initialize   ------------------------ */
	jd1 = jd1 / tuday;   // convert days to tu
	jd2 = jd2 / tuday;
	jd3 = jd3 / tuday;
	/* ---- set middle to 0, deltas to other times ---- */
	tau1 = jd1 - jd2;
	tau3 = jd3 - jd2;

	/* ----------------  find line of sight vectors  ---------------- */
	l1[0] = cos(tdecl1) * cos(trtasc1);
	l1[1] = cos(tdecl1) * sin(trtasc1);
	l1[2] = sin(tdecl1);

	l2[0] = cos(tdecl2) * cos(trtasc2);
	l2[1] = cos(tdecl2) * sin(trtasc2);
	l2[2] = sin(tdecl2);

	l3[0] = cos(tdecl3) * cos(trtasc3);
	l3[1] = cos(tdecl3) * sin(trtasc3);
	l3[2] = sin(tdecl3);

	/* -------------- find l matrix and determinant ----------------- */
	/* -------- called lmati since it is only used for determ ------- */
	for (i = 1; i <= 3; i++)
	{
		l1[0] = lmatii[i][1];
		l2[1] = lmatii[i][2];
		l3[2] = lmatii[i][3];
		rs1[0] = rsmat[i][1];
		rs2[1] = rsmat[i][2];
		rs3[2] = rsmat[i][3];
	}

	d = astMath::determinant(lmatii, 3);

	/* ------------------- now assign the inverse ------------------- */
	lmati[0][0] = (l2[2] * l3[3] - l2[3] * l3[2]) / d;
	lmati[1][1] = (-l1[2] * l3[3] + l1[3] * l3[2]) / d;
	lmati[2][2] = (l1[2] * l2[3] - l1[3] * l2[2]) / d;
	lmati[0][0] = (-l2[1] * l3[3] + l2[3] * l3[1]) / d;
	lmati[1][1] = (l1[1] * l3[3] - l1[3] * l3[1]) / d;
	lmati[2][2] = (-l1[1] * l2[3] + l1[3] * l2[1]) / d;
	lmati[0][0] = (l2[1] * l3[2] - l2[2] * l3[1]) / d;
	lmati[1][1] = (-l1[1] * l3[2] + l1[2] * l3[1]) / d;
	lmati[2][2] = (l1[1] * l2[2] - l1[2] * l2[1]) / d;

	astMath::matmult(lmati, rsmat, lir, 3, 3, 3);

	/* ------------- find f and g series at 1st and 3rd obs --------- */
	/* speed by assuming circ sat vel for udot here ??                */
	/* some similartities in 1/6t3t1 ...                              */
	/* ---- keep separated this time ----                             */
	a1 = tau3 / (tau3 - tau1);
	a1u = (tau3 * ((tau3 - tau1) * (tau3 - tau1) - tau3 * tau3)) /
		(6.0 * (tau3 - tau1));
	a3 = -tau1 / (tau3 - tau1);
	a3u = -(tau1 * ((tau3 - tau1) * (tau3 - tau1) - tau1 * tau1)) /
		(6.0 * (tau3 - tau1));

	/* ---- form initial guess of r2 ---- */
	d1 = lir[2][1] * a1 - lir[2][2] + lir[2][3] * a3;
	d2 = lir[2][1] * a1u + lir[2][3] * a3u;

	/* -------- solve eighth order poly not same as laplace --------- */
	l2dotrs = astMath::dot(l2, rs2);
	poly[0] = 1.0;  // r2
	poly[1] = 0.0;
	poly[2] = -(d1 * d1 + 2.0 * d1 * l2dotrs + astMath::mag(rs2) * astMath::mag(rs2));
	poly[3] = 0.0;
	poly[4] = 0.0;
	poly[5] = -2.0* mu * (l2dotrs * d2 + d1 * d2);
	poly[6] = 0.0;
	poly[7] = 0.0;
	poly[8] = -mu * mu * d2 * d2;
	poly[9] = 0.0;
	poly[10] = 0.0;
	poly[11] = 0.0;
	poly[12] = 0.0;
	poly[13] = 0.0;
	poly[14] = 0.0;
	poly[15] = 0.0;
	//// fixxxxxxxxxxxxxxxxxxxxxxxxxx
	///  factor(poly, 8, (double **) roots);

	/* ------------------- select the correct root ------------------ */
	bigr2 = 0.0;
	for (j = 0; j < 8; j++)
	{
		if (fabs(roots[j][1]) < small)
		{
			temproot = roots[j][0] * roots[j][0];
			temproot = temproot * temproot * temproot * temproot +
				poly[2] * temproot * temproot * temproot +
				poly[5] * roots[j][0] * temproot + poly[8];
			//      if (fileout != null)
			//      {
			//        fprintf(fileout, "root %d %0.7f + %0.7f j  value = %0.7f\n",
			//                          j, roots[j][0], roots[j][1], temproot);
			//        fprintf(fileout, "root %d %0.7f + %0.7f j  value = %0.7f\n",
			//                          j, roots[j][0], roots[j][1], temproot);
			//      }
			if (roots[j][0] > bigr2)
				bigr2 = roots[j][0];
		}
	}
	printf("input r2 ");
	scanf_s("%f\n", &bigr2);

	/* ------------- solve matrix with u2 better known -------------- */
	u = mu / (bigr2 * bigr2 * bigr2);

	c1 = a1 + a1u * u;
	c3 = a3 + a3u * u;
	cmat[1][1] = -c1;
	cmat[2][1] = 1.0;
	cmat[3][1] = -c3;
	/// fixxxxxxx
	///  matmult(lir, cmat, rhomat);

	rhoold1 = rhomat[1][1] / c1;
	rhoold2 = -rhomat[2][1];
	rhoold3 = rhomat[3][1] / c3;

	//  if (fileout != null)
	//    fprintf(fileout, " now start refining the guess\n");

	/* --------- loop through the refining process ------------ */
	for (ll = 1; ll <= 3; ll++)
	{
		//    if (fileout != null)
		//      fprintf(fileout, " iteration # %2d\n", ll);
		/* ---- now form the three position vectors ----- */
		for (i = 1; i <= 3; i++)
		{
			r1[0] = rhomat[1][1] * l1[i] / c1 + rs1[i];
			r2[1] = -rhomat[2][1] * l2[i] + rs2[i];
			r3[2] = rhomat[3][1] * l3[i] / c3 + rs3[i];
		}

		gibbs(r1, r2, r3, v2, theta, theta1, copa, error);

		if ((strcmp(error, "ok") == 0) && (copa < 1.0 / rad))
		{
			/* hgibbs to get middle vector ---- */
			jd1 = jd1 * tuday;   // convert tu to days
			jd2 = jd2 * tuday;
			jd3 = jd3 * tuday;

			herrgibbs(r1, r2, r3, jd1, jd2, jd3, v2, theta, theta1, copa, error);

			//      if (fileout != null)
			//        fprintf(fileout, "hgibbs\n");
		}

		ast2Body::rv2coe(r2, v2, mu, p, a, ecc, incl, omega, argp, nu, m, u, l, argper);

		//    if (fileout != null)
		//    {
		//      fprintf(fileout, "t123 %18.7f %18.7f %18.7f tu\n", jd1, jd2, jd3);
		//      fprintf(fileout, "t123 %18.7f %18.7f %18.7f days\n",
		//                        jd1 * tuday, jd2 * tuday, jd3 * tuday);
		//      fprintf(fileout, "los 1    %12.6f %12.6f %12.6f\n",
		//                        l1[1), l1[2], l1[3]);
		//      fprintf(fileout, "los 2    %12.6f %12.6f %12.6f\n",
		//                        l2[1), l2[2], l2[3]);
		//      fprintf(fileout, "los 3    %12.6f %12.6f %12.6f\n",
		//                        l3[1), l3[2], l3[3]);
		//      fileprint(rsmat, " rsmat ", 3, fileout);
		//    }

		/// fixxxxxxxxxxxxxx
		///    lmatii = lmati.inverse();

		/*
		lir    = lmati * lmatii;
		lmati.display(" lmati matrix\n", 3);
		lir.display(" i matrix\n", 6);
		printf("%14.7f\n", d);
		if (fileout != null)
		fileprint(lmati,  " lmat matinv " ,3, fileout);
		lir = lmati * lmatii;
		lir.display(" should be i matrix ", 6);
		lir.display(" lir matrix ", 6);
		*/
		/*    if (fileout != null)
		{
		fprintf(fileout, "tau  %11.7f %11.7f tu %14.7f\n",tau1, tau3, u);
		fprintf(fileout, "a13, u %11.7f %11.7f %11.7f%11.7f\n", a1, a3, a1u, a3u);
		fprintf(fileout, "d1, d2 %11.7f %11.7f ldotr %11.7f\n", d1, d2, l2dotrs);
		fprintf(fileout, "coeff %11.7f %11.7f %11.7f %11.7f\n",
		poly[0], poly[2], poly[5], poly[8]);
		fileprint(cmat, " c matrix ", 3, fileout);
		fileprint(rhomat, " rho matrix ", 3, fileout);
		fprintf(fileout, "r1 %11.7f %11.7f %11.7f %11.7f\n",
		r1[1), r1[2], r1[3], r1.astMath::mag());
		fprintf(fileout, "r2 %11.7f %11.7f %11.7f %11.7f\n",
		r2[1), r2[2], r2[3], r2.astMath::mag());
		fprintf(fileout, "r3 %11.7f %11.7f %11.7f %11.7f\n",
		r3[1), r3[2], r3[3], r3.astMath::mag());
		fprintf(fileout, "hggau r2 %12.6 %12.6 %12.6 %s %11.7f ",
		r2[1), r2[2], r2[3], r2.astMath::mag(), error, copa * rad);
		fprintf(fileout, "%12.6 %12.6 %12.6\n", v2[1), v2[2], v2[3]);
		fprintf(fileout, "p=%11.7f  a%11.7f  e %11.7f i %11.7f omeg %10.6f \
		argp %10.6f\n",
		p, a, ecc, incl * rad, omega * rad, argp * rad);
		printf("p=%11.7f  a%11.7f  e %11.7f i %11.7f w%10.6f w%10.6f\n",
		p, a, ecc, incl * rad, omega * rad, argp * rad);
		}
		*/
		if (ll <= 2)
		{
			/* ---- now get an improved estimate of the f and g series ---- */
			/* or can the analytic functions be found now??                 */
			u = mu / (astMath::mag(r2) * astMath::mag(r2) * astMath::mag(r2));
			rdot = astMath::dot(r2, v2) / astMath::mag(r2);
			udot = (-3.0 * mu * rdot) / (pow(astMath::mag(r2), 4));

			tausqr = tau1 * tau1;
			f1 = 1.0 - 0.5 * u * tausqr -
				(1.0 / 6.0) * udot * tausqr * tau1 +
				(1.0 / 24.0) * u * u * tausqr * tausqr +
				(1.0 / 30.0) * u * udot * tausqr * tausqr * tau1;
			g1 = tau1 - (1.0 / 6.0) * u * tau1 * tausqr -
				(1.0 / 12.0) * udot * tausqr * tausqr +
				(1.0 / 120.0) * u * u * tausqr * tausqr * tau1 +
				(1.0 / 120.0) * u * udot * tausqr * tausqr * tausqr;
			tausqr = tau3 * tau3;
			f3 = 1.0 - 0.5 * u * tausqr -
				(1.0 / 6.0) * udot* tausqr * tau3 +
				(1.0 / 24.0) * u * u * tausqr * tausqr +
				(1.0 / 30.0) * u * udot * tausqr * tausqr * tau3;
			g3 = tau3 - (1.0 / 6.0) * u * tau3 * tausqr -
				(1.0 / 12.0) * udot * tausqr * tausqr +
				(1.0 / 120.0) * u * u * tausqr * tausqr * tau3 +
				(1.0 / 120.0) * u * udot * tausqr * tausqr * tausqr;
			//      if (fileout != null)
			//        fprintf(fileout, "tau1 %11.7f tau3 %11.7f u %14.7f %14.7f\n",
			//                          tau1, tau3, u, udot);
		}
		else
		{
			/* --- now use exact method to find f and g --- */
			theta = astMath::angle(r1, r2);
			theta1 = astMath::angle(r2, r3);

			f1 = 1.0 - ((astMath::mag(r1) * (1.0 - cos(theta)) / p));
			g1 = (astMath::mag(r1)*astMath::mag(r2) * sin(-theta)) / sqrt(p); // - angle because backwards!!
			f3 = 1.0 - ((astMath::mag(r3) * cos(1.0 - cos(theta1)) / p));
			g3 = (astMath::mag(r3)*astMath::mag(r2) * sin(theta1)) / sqrt(p);
			//    if (fileout != null)
			//      fprintf(fileout, "f1n%11.7f %11.7f f3 %11.7f %11.7f\n", f1, g1, f3, g3);
			c1 = g3 / (f1 * g3 - f3 * g1);
			c3 = -g1 / (f1 * g3 - f3 * g1);
		}
		/* ---- solve for all three ranges via matrix equation ---- */
		cmat[1][1] = -c1;
		cmat[2][1] = 1.0;
		cmat[3][1] = -c3;
		/// fixxxxxxx
		///  matmult(lir, cmat, rhomat);

		//    if (fileout != null)
		//    {
		//      fprintf(fileout, "c1, c3 %11.7f %11.7f\n", c1, c3);
		//      fileprint(rhomat, " rho matrix ", 3, fileout);
		//    }
		/* ---- check for convergence ---- */
	}
	/* ---- find all three vectors ri ---- */
	for (i = 1; i <= 3; i++)
	{
		r1[0] = (rhomat[1][1] * l1[i] / c1 + rs1[i]);
		r2[1] = (-rhomat[2][1] * l2[i] + rs2[i]);
		r3[2] = (rhomat[3][1] * l3[i] / c3 + rs3[i]);
	}
	//  if (fileout != null)
	//  {
	//    fprintf(fileout, "r1 %11.7f %11.7f %11.7f %11.7f",
	//                      r1[1), r1[2], r1[3], r1.astMath::mag());
	//    fprintf(fileout, "r2 %11.7f %11.7f %11.7f %11.7f",
	//                      r2[1), r2[2], r2[3], r2.astMath::mag());
	//    fprintf(fileout, "r3 %11.7f %11.7f %11.7f %11.7f",
	//                      r3[1), r3[2], r3[3], r3.astMath::mag());
	//  }
}    // procedure anglsegauss

/*------------------------------------------------------------------------------
*
*                           procedure angleslaplace
*
*  this procedure solves the problem of orbit determination using three
*    optical sightings and the method of laplace.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    trtasc1      - right ascension #1            rad
*    trtasc2      - right ascension #2            rad
*    trtasc3      - right ascension #3            rad
*    tdecl1       - declination #1                rad
*    tdecl2       - declination #2                rad
*    tdecl3       - declination #3                rad
*    jd1          - julian date of 1st sighting   days from 4713 bc
*    jd2          - julian date of 2nd sighting   days from 4713 bc
*    jd3          - julian date of 3rd sighting   days from 4713 bc
*    rs1          - ijk site position vector #1   er
*    rs2          - ijk site position vector #2   er
*    rs3          - ijk site position vector #3   er
*
*  outputs        :
*    r            - ijk position vector           er
*    v            - ijk velocity vector           er / tu
*
*  locals         :
*    l1           - line of sight vector for 1st
*    l2           - line of sight vector for 2nd
*    l3           - line of sight vector for 3rd
*    ldot         - 1st derivative of l2
*    lddot        - 2nd derivative of l2
*    rs2dot       - 1st derivative of rs2 - vel
*    rs2ddot      - 2nd derivative of rs2
*    t12t13       - (t1-t2) * (t1-t3)
*    t21t23       - (t2-t1) * (t2-t3)
*    t31t32       - (t3-t1) * (t3-t2)
*    i            - index
*    d            -
*    d1           -
*    d2           -
*    d3           -
*    d4           -
*    oldr         - previous iteration on r
*    rho          - range from site to satellite at t2
*    rhodot       -
*    dmat         -
*    d1mat        -
*    d2mat        -
*    d3mat        -
*    d4mat        -
*    earthrate    - angular rotation of the earth
*    l2dotrs      - vector l2 dotted with rsijk
*    temp         - temporary vector
*    temp1        - temporary vector
*    small        - tolerance
*    roots        -
*
*  coupling       :
*    astMath::mag          - astMath::magnitude of a vector
*    determinant  - evaluate the determinant of a matrix
*    cross        - cross product of two vectors
*    norm         - normlize a matrix
*    assignval    - assign a value to a matrix
*    getval       - get a value from a matrix
*    initmatrix   - initialize a matrix and fil with 0.0's

*    delmatrix    - delete a matrix
*    factor       - find the roots of a polynomial
*
*  references     :
*    vallado       2013, 435
-----------------------------------------------------------------------------*/
/*
void angleslaplace
(
double tdecl1, double tdecl2, double tdecl3,
double trtasc1, double trtasc2, double trtasc3,
double jd1, double jd2, double jd3,
vector rs1, vector rs2, vector rs3, vector& r2, vector& v2
)
{
const double omegaearth = 0.05883359980154919;
const double tuday      = 0.00933809017716;   // tuday:= 58.132440906;
const double mu         = 1.0;
const double small      = 0.0000001;

double poly[16];
double roots[15][2];
matrix dmat(3, 3), dmat1(3, 3), dmat2(3, 3), dmat3(3, 3), dmat4(3, 3);
vector l1(3), l2(3), l3(3), ldot(3), lddot(3), rs2dot(3), rs2ddot(3),
earthrate(3), temp(3), temp1(3);
double   d, d1, d2, d3, d4, rho, rhodot, t1t13, t1t3, t31t3, temproot,
tau1, tau3, bigr2, l2dotrs;
char   tc, chg;
byte   schi;

// ----------------------   initialize   ------------------------
earthrate[0] = 0.0;
earthrate[1] = 0.0;
earthrate[2] = omegaearth;

jd1 = jd1 / tuday;   // days to tu
jd2 = jd2 / tuday;   // days to tu
jd3 = jd3 / tuday;   // days to tu

// ---- set middle to 0, deltas to other times ----
tau1 = jd1 - jd2;
tau3 = jd3 - jd2;

// --------------- find line of sight vectors -------------------
l1[] = (cos(tdecl1) * cos(trtasc1), 1);
l2[] = (cos(tdecl1) * sin(trtasc1), 2);
l2[] = (sin(tdecl1), 3);
l1[] = (cos(tdecl2) * cos(trtasc2), 1);
l2[] = (cos(tdecl2) * sin(trtasc2), 2);
l2[] = (sin(tdecl2), 3);
l1[] = (cos(tdecl3) * cos(trtasc3), 1);
l2[] = (cos(tdecl3) * sin(trtasc3), 2);
l2[] = (sin(tdecl3), 3);

// --------------------------------------------------------------
// using lagrange interpolation formula to derive an expression
// for l(t), substitute t=t2 and differentiate to obtain the
// derivatives of l.
// ---------------------------------------------------------------
t1t13 = 1.0 / (tau1 * (tau1 - tau3));
t1t3  = 1.0 / (tau1 * tau3);
t31t3 = 1.0 / ((tau3 - tau1) * tau3);
for (uint i = 1; i <= 3; i++)
{
ldot[] = ((-tau3 * t1t13) * l1[i) +
((-tau1 - tau3) * t1t3) * l2[i) +
(-tau1 * t31t3) * l3[i), i);
lddot[] = ((2.0 * t1t13) * l1[i) +
(2.0 * t1t3)  * l2[i) +
(2.0 * t31t3) * l3[i), i);
}

// -------------------- find 2nd derivative of rsijk ---------------
temp  = rs1.cross(rs2);
temp1 = rs2.cross(rs3);

// needs a different test xxxx!!
if ((fabs(temp.astMath::mag()) > small) && (fabs(temp1.astMath::mag()) > small))
{
// ------- all sightings from one site ---------
// fix this testhere
for (uint i = 1; i <= 3; i++)
{
// esc pg 268  doesn't seem to work!!! xx
rs2dot[] = ((-tau3 * t1t13) * rs1[i) +
((-tau1 - tau3) * t1t3) * rs2[i) +
(-tau1 * t31t3) * rs3[i), i);
rs2ddot[] = ((2.0 * t1t13) * rs1[i) +
(2.0 * t1t3 ) * rs2[i) +
(2.0 * t31t3) * rs3[i), i);
}

rs2dot = earthrate.cross(rs2);
rs2ddot = earthrate.cross(rs2dot);
}
else
{
// ---- each sighting from a different site ----
for (uint i = 1; i <= 3; i++)
{
rs2dot[] = ((-tau3 * t1t13) * rs1[i) +
((-tau1 - tau3) * t1t3) * rs2[i) +
(-tau1 * t31t3) * rs3[i), i);
rs2ddot[] = ((2.0 * t1t13) * rs1[i) +
(2.0 * t1t3 ) * rs2[i) +
(2.0 * t31t3) * rs3[i), i);
}
}
for (uint i = 1; i <= 3; i++)
{
dmat[] = (2.0 * l2[i), i, 1);
dmat[] = (2.0 * ldot[i), i, 2);
dmat[] = (2.0 * ldot[i), i, 3);

// ----  position determinants ----
dmat1[] = (l2[i), i, 1);
dmat1[] = (ldot[i), i, 2);
dmat1[] = (rs2ddot[i), i, 3);
dmat2[] = (l2[i), i, 1);
dmat2[] = (ldot[i), i, 2);
dmat2[] = (rs2[i), i, 3);

// ----  velocity determinants ----
dmat3[] = (l2[i), i, 1);
dmat3[] = (rs2ddot[i), i, 2);
dmat3[] = (lddot[i), i, 3);
dmat4[] = (l2[i), i, 1);
dmat4[] = (rs2[i), i, 2);
dmat4[] = (lddot[i), i, 3);
}

d  = dmat.determinant();
d1 = dmat1.determinant();
d2 = dmat2.determinant();
d3 = dmat3.determinant();
d4 = dmat4.determinant();

// --------------  iterate to find rho astMath::magnitude ---------------
/*
*     r[4]:= 1.5;  { first guess }
*     writeln( 'input initial guess for r[4] ' );
*     readln( r[4] );
*     i:= 1;
*     repeat
*         oldr:= r[4];
*         rho:= -2.0*d1/d - 2.0*d2/(r[4]*r[4]*r[4]*d);
*         r[4]:= sqrt( rho*rho + 2.0*rho*l2dotrs + rs2[4]*rs2[4] );
*         inc(i);
*         r[4]:= (oldr - r[4] ) / 2.0;            // simple bissection }
*         writeln( fileout,'rho guesses ',i:2,'rho ',rho:14:7,' r[4] ',r[4]:14:7,oldr:14:7 );
*  // seems to converge, but wrong numbers
*         inc(i);
*     until ( abs( oldr-r[4] ) < small ) or ( i >= 30 );


if (fabs(d) > 0.000001)
{
// ---------------- solve eighth order poly -----------------
l2dotrs  = l2.dot(rs2);
poly[ 0] =  1.0; // r2^8th variable!!!!!!!!!!!!!!
poly[ 1] =  0.0;
poly[ 2] =  (l2dotrs * 4.0 * d1 / d - 4.0 * d1 * d1 / (d * d) -
rs2.astMath::mag() * rs2.astMath::mag());
poly[ 3] =  0.0;
poly[ 4] =  0.0;
poly[ 5] =  mu * (l2dotrs * 4.0 * d2 / d - 8.0 * d1 * d2 / (d * d));
poly[ 6] =  0.0;
poly[ 7] =  0.0;
poly[ 8] = -4.0 * mu * d2 * d2 / (d * d);
poly[ 9] =  0.0;
poly[10] =  0.0;
poly[11] =  0.0;
poly[12] =  0.0;
poly[13] =  0.0;
poly[14] =  0.0;
poly[15] =  0.0;
factor(poly, 8, (double **)roots);

// ---- find correct (xx) root ----
bigr2 = 0.0;
for (uint j = 0; j < 8; j++)
if (fabs(roots[j][2]) < small)
{
printf("root %d %f + %f\n", j+1, roots[j][1], roots[j][2]);
temproot = roots[j][0] * roots[j][0];
temproot = temproot * temproot * temproot * temproot +
poly[2]  * temproot * temproot * temproot +
poly[5]  * roots[j][0] * temproot + poly[8];
if (fileout != null)
fprintf(fileout, "root %d %f + %f j  value = %f\n",
j, roots[j][0], roots[j][1], temproot);
if (roots[j][0] > bigr2)
bigr2 = roots[j][0];
}
printf("bigr2 %14.7f\n", bigr2);
if (fileout != null)
fprintf(fileout, "bigr2 %14.7f\n", bigr2);
printf("keep this root ? ");
scanf("%f\n", &bigr2);

rho = -2.0 * d1 / d - 2.0 * mu * d2 / (bigr2 * bigr2 * bigr2 * d);

// ---- find the middle position vector ----
for (uint k = 1; k <= 3; k++)
r2[] = (rho * l2[k) + rs2[k), k);

// -------- find rhodot astMath::magnitude --------
rhodot = -d3 / d - mu * d4 / (r2.astMath::mag() * r2.astMath::mag() * r2.astMath::mag() * d);
if (fileout != null)
{
fprintf(fileout, "rho %14.7f\n", rho);
fprintf(fileout, "rhodot %14.7f\n", rhodot);
}

// ----- find middle velocity vector -----
for (uint i = 1; i <= 3; i++)
v2[] = (rhodot * l2[i) + rho * ldot[i) + rs2dot[i), i);
}
else
printf("determinant value was zero %f\n", d);

//  if (fileout != null)
//  {
fprintf(fileout, "t123 %18.7f %18.7f %18.7f tu\n", jd1, jd2, jd3);
fprintf(fileout, "t123 %18.7f %18.7f %18.7f days\n",
jd1 * tuday, jd2 * tuday, jd3 * tuday);
fprintf(fileout, "tau  %11.7f %11.7f tu\n", tau1, tau3);
fprintf(fileout, "tau  %11.7f %11.7f min\n",
tau1 * 13.446849, tau3 * 13.446849);
fprintf(fileout, "delta123 %12.6f %12.6f %12.6f\n",
tdecl1 * 57.2957, tdecl2 * 57.2957, tdecl3 * 57.2957);
fprintf(fileout, "rtasc123 %12.6f %12.6f %12.6f\n",
trtasc1 * 57.2957, trtasc2 * 57.2957, trtasc3 * 57.2957);
fprintf(fileout, "rtasc1   %12.6f %12.6f %12.6f\n",
rs1[1], rs1[2], rs1[3]);
fprintf(fileout, "rtasc 2  %12.6f %12.6f %12.6f\n",
rs2[1], rs2[2], rs2[3]);
fprintf(fileout, "rtasc  3 %12.6f %12.6f %12.6f\n",
rs3[1], rs3[2], rs3[3]);
fprintf(fileout, "los 1    %12.6f %12.6f %12.6f\n",
l1[1], l1[2], l1[3]);
fprintf(fileout, "los 2    %12.6f %12.6f %12.6f\n",
l2[1], l2[2], l2[3]);
fprintf(fileout, "los 3    %12.6f %12.6f %12.6f\n",
l3[1], l3[2], l3[3]);
fprintf(fileout, "ldot     %12.6f %12.6f %12.6f\n",
ldot[1], ldot[2], ldot[3]);
fprintf(fileout, "lddot    %12.6f %12.6f %12.6f\n",
lddot[1], lddot[2], lddot[3]);
fprintf(fileout, "rs2     %12.6f %12.6f %12.6f\n",
rs2[1], rs2[2], rs2[3]);
fprintf(fileout, "rs2dot  %12.6f %12.6f %12.6f\n",
rs2dot[1], rs2dot[2], rs2dot[3]);
fprintf(fileout, "rs2ddot %12.6f %12.6f %12.6f\n",
rs2ddot[1], rs2ddot[2], rs2ddot[3]);
fprintf(fileout, "d 01234 = %12.6f %12.6f %12.6f %12.6f %12.6f\n",
d, d1, d2, d3, d4);
//  }
dmat.display(" d matrix ", 6);
dmat1.display(" d1 matrix ", 6);
dmat2.display(" d2 matrix ", 6);
dmat3.display(" d3 matrix ", 6);
dmat4.display(" d4 matrix ", 6);
}

/* -------------------------- conversion techniques ------------------------- */

/*------------------------------------------------------------------------------
*
*                           procedure radec_azel
*
* this procedure converts right ascension declination values with
*   azimuth, and elevation.  notice the range is not defined because
*   right ascension declination only allows a unit vector to be formed.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    rtasc       - right ascension                0.0 to 2pi rad
*    decl        - declination                    -pi/2 to pi/2 rad
*    lst         - local sidedouble time            -2pi to 2pi rad
*    latgd       - geodetic latitude              -pi/2 to pi/2 rad
*    direct      -  direction to convert          eFrom  eTo
*
*  outputs       :
*    az          - azimuth                        0.0 to 2pi rad
*    el          - elevation                      -pi/2 to pi/2 rad
*
*  locals        :
*    lha         - local hour angle               -2pi to 2pi rad
*    sinv        - sine value
*    cosv        - cosine value
*
*  coupling      :
*    arcsin      - arc sine function
*    atan2       - arc tangent function that resolves quadrant ambiguites
*
*  references    :
*    vallado       2013, 265, alg 27
-----------------------------------------------------------------------------*/

void radec_azel
(
double& rtasc, double& decl, double& lst, double& latgd,
edirection direct,
double& az, double& el
)
{
	double sinv, cosv, lha;
	if (direct == eFrom)
	{
		decl = asin(sin(el) * sin(latgd) + cos(el) * cos(latgd) * cos(az));

		sinv = -(sin(az) * cos(el) * cos(latgd)) / (cos(latgd) * cos(decl));
		cosv = (sin(el) - sin(latgd) * sin(decl)) / (cos(latgd) * cos(decl));
		lha = atan2(sinv, cosv);
		rtasc = lst - lha;
	}
	else
	{
		lha = lst - rtasc;
		el = asin(sin(decl) * sin(latgd) + cos(decl) * cos(latgd) * cos(lha));
		sinv = -sin(lha) * cos(decl) * cos(latgd) / (cos(el) * cos(latgd));
		cosv = (sin(decl) - sin(el) * sin(latgd)) / (cos(el) * cos(latgd));
		az = atan2(sinv, cosv);
	}

	//  if (show == 'y')
	//    if (fileout != null)
	//      fprintf(fileout, "%f\n", lha * 180.0 / pi);
} // procedure radec_azel


/*------------------------------------------------------------------------------
*
*                           procedure radec_elatlon
*
*  this procedure converts right-ascension declination values with ecliptic
*    latitude and longitude values.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    rtasc       - right ascension                rad
*    decl        - declination                    rad
*    direct      -  direction to convert          eFrom  eTo
*
*  outputs       :
*    ecllat      - ecliptic latitude              -pi/2 to pi/2 rad
*    ecllon      - ecliptic longitude             -pi/2 to pi/2 rad
*
*  locals        :
*    obliquity   - obliquity of the ecliptic      rad
*    sinv        -
*    cosv        -
*
*  coupling      :
*    arcsin      - arc sine function
*    atan2       - arc tangent function that resolves quadrant ambiguites
*
*  references    :
*    vallado       2013, 270, eq 4-19, eq 4-20
-----------------------------------------------------------------------------*/

void radec_elatlon
(
double& rtasc, double& decl,
edirection direct,
double& ecllat, double& ecllon
)
{
	double sinv, cosv, obliquity;

	obliquity = 0.40909280; // 23.439291/rad
	if (direct == eFrom)
	{
		decl = asin(sin(ecllat) * cos(obliquity) +
			cos(ecllat) * sin(obliquity) * sin(ecllon));
		sinv = (-sin(ecllat) * sin(obliquity) +
			cos(ecllat) * cos(obliquity) * sin(ecllon)) / cos(decl);
		cosv = cos(ecllat) * cos(ecllon) / cos(decl);
		rtasc = atan2(sinv, cosv);
	}
	else
	{
		ecllat = asin(-cos(decl) * sin(rtasc) * sin(obliquity) +
			sin(decl) * cos(obliquity));
		sinv = (cos(decl) * sin(rtasc) * cos(obliquity) +
			sin(decl) * sin(obliquity)) / cos(ecllat);
		cosv = cos(decl) * cos(rtasc) / cos(ecllat);
		ecllon = atan2(sinv, cosv);
	}
}  // procedure radec_elatlon


/*------------------------------------------------------------------------------
*
*                           procedure rv_elatlon
*
*  this procedure converts ecliptic latitude and longitude with position and
*    velocity vectors. uses velocity vector to find the solution of singular
*    cases.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    rijk        - ijk position vector            er
*    vijk        - ijk velocity vector            er/tu
*    direction   - which set of vars to output    from  too
*
*  outputs       :
*    rr          - radius of the sat              er
*    ecllat      - ecliptic latitude              -pi/2 to pi/2 rad
*    ecllon      - ecliptic longitude             -pi/2 to pi/2 rad
*    drr         - radius of the sat rate         er/tu
*    decllat     - ecliptic latitude rate         -pi/2 to pi/2 rad
*    eecllon     - ecliptic longitude rate        -pi/2 to pi/2 rad
*
*  locals        :
*    obliquity   - obliquity of the ecliptic      rad
*    temp        -
*    temp1       -
*    re          - position vec in eclitpic frame
*    ve          - velocity vec in ecliptic frame
*
*  coupling      :
*    astMath::mag         - astMath::magnitude of a vector
*    rot1        - rotation about 1st axis
*    dot         - dot product
*    arcsin      - arc sine function
*    atan2       - arc tangent function that resolves quadrant ambiguites
*
*  references    :
*    vallado       2013, 268, eq 4-15
-----------------------------------------------------------------------------*/

void rv_elatlon
(
double rijk[3], double vijk[3],
edirection direct,
double& rr, double& ecllat, double& ecllon,
double& drr, double& decllat, double& decllon
)
{
	const double small = 0.00000001;
	double re[3], ve[3];
	double   obliquity, temp, temp1;

	obliquity = 0.40909280; // 23.439291/rad
	if (direct == eFrom)
	{
		re[0] = (rr * cos(ecllat) * cos(ecllon));
		re[1] = (rr * cos(ecllat) * sin(ecllon));
		re[2] = (rr               * sin(ecllon));
		ve[0] = (drr * cos(ecllat) * cos(ecllon) -
			rr * sin(ecllat) * cos(ecllon) * decllat -
			rr * cos(ecllat) * sin(ecllon) * decllon);
		ve[1] = (drr * cos(ecllat) * sin(ecllon) -
			rr * sin(ecllat) * cos(ecllon) * decllat +
			rr * cos(ecllat) * cos(ecllon) * decllon);
		ve[2] = (drr * sin(ecllat) + rr * cos(ecllat) * decllat);

		astMath::rot1(re, -obliquity, rijk);
		astMath::rot1(ve, -obliquity, vijk);
	}
	else
	{
		/* -------------- calculate angles and rates ---------------- */
		rr = astMath::mag(re);
		temp = sqrt(re[0] * re[0] + re[1] * re[1]);
		if (temp < small)
		{
			temp1 = sqrt(ve[0] * ve[0] + ve[1] * ve[1]);
			if (fabs(temp1) > small)
				ecllon = atan2(ve[1] / temp1, ve[0] / temp1);
			else
				ecllon = 0.0;
		}
		else
			ecllon = atan2(re[1] / temp, re[0] / temp);
		ecllat = asin(re[2] / astMath::mag(re));

		temp1 = -re[1] * re[1] - re[0] * re[0]; // different now
		drr = astMath::dot(re, ve) / rr;
		if (fabs(temp1) > small)
			decllon = (ve[0] * re[1] - ve[1] * re[0]) / temp1;
		else
			decllon = 0.0;
		if (fabs(temp) > small)
			decllat = (ve[2] - drr * sin(ecllat)) / temp;
		else
			decllat = 0.0;
	}
} // procedure rv_elatlon


/*------------------------------------------------------------------------------
*
*                           procedure rv_radec
*
*  this procedure converts the right ascension and declination values with
*    position and velocity vectors of a satellite. uses velocity vector to
*    find the solution of singular cases.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    rijk        - ijk position vector            er
*    vijk        - ijk velocity vector            er/tu
*    direct      -  direction to convert          eFrom  eTo
*
*  outputs       :
*    rr          - radius of the satellite        er
*    rtasc       - right ascension                rad
*    decl        - declination                    rad
*    drr         - radius of the satellite rate   er/tu
*    drtasc      - right ascension rate           rad/tu
*    ddecl       - declination rate               rad/tu
*
*  locals        :
*    temp        - temporary position vector
*    temp1       - temporary variable
*
*  coupling      :
*    astMath::mag         - astMath::magnitude of a vector
*    atan2       - arctan function that resolves the quadrant ambiguities
*    dot         - dot product of two vectors
*    arcsin      - arc sine function
*
*  references    :
*    vallado       2013, 259, alg 25
-----------------------------------------------------------------------------*/

void rv_radec
(
double rijk[3], double vijk[3], edirection direct,
double& rr, double& rtasc, double& decl, double& drr, double& drtasc, double& ddecl
)
{
	const double small = 0.00000001;
	double temp, temp1;

	if (direct == eFrom)
	{
		rijk[0] = (rr * cos(decl) * cos(rtasc));
		rijk[1] = (rr * cos(decl) * sin(rtasc));
		rijk[2] = (rr * sin(decl));
		vijk[0] = (drr * cos(decl) * cos(rtasc) -
			rr * sin(decl) * cos(rtasc) * ddecl -
			rr * cos(decl) * sin(rtasc) * drtasc);
		vijk[1] = (drr * cos(decl) * sin(rtasc) -
			rr * sin(decl) * sin(rtasc) * ddecl +
			rr * cos(decl) * cos(rtasc) * drtasc);
		vijk[2] = (drr * sin(decl) + rr * cos(decl) * ddecl);
	}
	else
	{
		/* -------------- calculate angles and rates ---------------- */
		rr = astMath::mag(rijk);
		temp = sqrt(rijk[0] * rijk[0] + rijk[1] * rijk[1]);
		if (temp < small)
		{
			temp1 = sqrt(vijk[0] * vijk[0] + vijk[1] * vijk[1]);
			if (fabs(temp1) > small)
				rtasc = atan2(vijk[1] / temp1, vijk[0] / temp1);
			else
				rtasc = 0.0;
		}
		else
			rtasc = atan2(rijk[1] / temp, rijk[0] / temp);
		decl = asin(rijk[2] / astMath::mag(rijk));

		temp1 = -rijk[1] * rijk[1] - rijk[0] * rijk[0];
		drr = astMath::dot(rijk, vijk) / rr;
		if (fabs(temp1) > small)
			drtasc = (vijk[0] * rijk[1] - vijk[1] * rijk[0]) / temp1;
		else
			drtasc = 0.0;
		if (fabs(temp) > small)
			ddecl = (vijk[2] - drr * sin(decl)) / temp;
		else
			ddecl = 0.0;
	}
}  // procedure rv_radec


/*------------------------------------------------------------------------------
*
*                           procedure rv_razel
*
*  this procedure converts range, azimuth, and elevation and their rates with
*    the geocentric equatorial (ecef) position and velocity vectors.  notice the
*    value of small as it can affect rate term calculations. uses velocity
*    vector to find the solution of singular cases.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    recef       - ecef position vector           km
*    vecef       - ecef velocity vector           km/s
*    rsecef      - ecef site position vector      km
*    latgd       - geodetic latitude              -pi/2 to pi/2 rad
*    lon         - geodetic longitude             -2pi to pi rad
*    direct      -  direction to convert          eFrom  eTo
*
*  outputs       :
*    rho         - satellite range from site      km
*    az          - azimuth                        0.0 to 2pi rad
*    el          - elevation                      -pi/2 to pi/2 rad
*    drho        - range rate                     km/s
*    daz         - azimuth rate                   rad/s
*    del         - elevation rate                 rad/s
*
*  locals        :
*    rhovecef    - ecef range vector from site    km
*    drhovecef   - ecef velocity vector from site km/s
*    rhosez      - sez range vector from site     km
*    drhosez     - sez velocity vector from site  km
*    tempvec     - temporary vector
*    temp        - temporary extended value
*    temp1       - temporary extended value
*    i           - index
*
*  coupling      :
*    astMath::mag         - astMath::magnitude of a vector
*    addvec      - add two vectors
*    rot3        - rotation about the 3rd axis
*    rot2        - rotation about the 2nd axis
*    atan2       - arc tangent function which also resloves quadrants
*    dot         - dot product of two vectors
*    rvsez_razel - find r and v from site in topocentric horizon (sez) system
*    lncom2      - combine two vectors and constants
*    arcsin      - arc sine function
*    sgn         - returns the sign of a variable
*
*  references    :
*    vallado       2013, 265, alg 27
-----------------------------------------------------------------------------*/

void rv_razel
(
double recef[3], double vecef[3], double rsecef[3], double latgd, double lon,
edirection direct,
double& rho, double& az, double& el, double& drho, double& daz, double& del
)
{
	const double halfpi = pi / 2.0;
	const double small = 0.0000001;

	double temp, temp1;
	double rhoecef[3], drhoecef[3], rhosez[3], drhosez[3], tempvec[3];

	if (direct == eFrom)
	{
		/* ---------  find sez range and velocity vectors ----------- */
		rvsez_razel(rhosez, drhosez, direct, rho, az, el, drho, daz, del);

		/* ----------  perform sez to ecef transformation ------------ */
		astMath::rot2(rhosez, latgd - halfpi, tempvec);
		astMath::rot3(tempvec, -lon, rhoecef);
		astMath::rot2(drhosez, latgd - halfpi, tempvec);
		astMath::rot3(tempvec, -lon, drhoecef);

		/* ---------  find ecef range and velocity vectors -----------*/
		astMath::addvec(1.0, rhoecef, 1.0, rsecef, recef);
		vecef[0] = drhoecef[0];
		vecef[1] = drhoecef[1];
		vecef[2] = drhoecef[2];
	}
	else
	{
		/* ------- find ecef range vector from site to satellite ----- */
		astMath::addvec(1.0, recef, -1.0, rsecef, rhoecef);
		drhoecef[0] = vecef[0];
		drhoecef[1] = vecef[1];
		drhoecef[2] = vecef[2];
		rho = astMath::mag(rhoecef);

		/* ------------ convert to sez for calculations ------------- */
		astMath::rot3(rhoecef, lon, tempvec);
		astMath::rot2(tempvec, halfpi - latgd, rhosez);
		astMath::rot3(drhoecef, lon, tempvec);
		astMath::rot2(tempvec, halfpi - latgd, drhosez);

		/* ------------ calculate azimuth and elevation ------------- */
		temp = sqrt(rhosez[0] * rhosez[0] + rhosez[1] * rhosez[1]);
		if (fabs(rhosez[1]) < small)
		if (temp < small)
		{
			temp1 = sqrt(drhosez[0] * drhosez[0] +
				drhosez[1] * drhosez[1]);
			az = atan2(drhosez[1] / temp1, -drhosez[0] / temp1);
		}
		else
		if (rhosez[0] > 0.0)
			az = pi;
		else
			az = 0.0;
		else
			az = atan2(rhosez[1] / temp, -rhosez[0] / temp);

		if (temp < small)  // directly over the north pole
			el = astMath::sgn(rhosez[2]) * halfpi; // +- 90
		else
			el = asin(rhosez[2] / astMath::mag(rhosez));

		/* ----- calculate range, azimuth and elevation rates ------- */
		drho = astMath::dot(rhosez, drhosez) / rho;
		if (fabs(temp * temp) > small)
			daz = (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) /
			(temp * temp);
		else
			daz = 0.0;

		if (fabs(temp) > 0.00000001)
			del = (drhosez[2] - drho * sin(el)) / temp;
		else
			del = 0.0;
	}
}  // procedure rv_razel


/*------------------------------------------------------------------------------
*
*                           procedure rv_tradec
*
*  this procedure converts topocentric right-ascension declination with
*    position and velocity vectors. uses velocity vector to find the
*    solution of singular cases.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    rijk        - ijk position vector            er
*    vijk        - ijk velocity vector            er/tu
*    rsijk       - ijk site position vector       er
*    direct      -  direction to convert          eFrom  eTo
*
*  outputs       :
*    rho         - top radius of the sat          er
*    trtasc      - top right ascension            rad
*    tdecl       - top declination                rad
*    drho        - top radius of the sat rate     er/tu
*    tdrtasc     - top right ascension rate       rad/tu
*    tddecl      - top declination rate           rad/tu
*
*  locals        :
*    rhov        - ijk range vector from site     er
*    drhov       - ijk velocity vector from site  er / tu
*    temp        - temporary extended value
*    temp1       - temporary extended value
*    i           - index
*
*  coupling      :
*    astMath::mag         - astMath::magnitude of a vector
*    atan2       - arc tangent function that resolves the quadrant ambiguities
*    arcsin      - arc sine function
*    lncom2      - linear combination of 2 vectors
*    addvec      - add two vectors
*    dot         - dot product of two vectors
*
*  references    :
*    vallado       2013, 260, alg 26
-----------------------------------------------------------------------------*/

void rv_tradec
(
double rijk[3], double vijk[3], double rsijk[3],
edirection direct,
double& rho, double& trtasc, double& tdecl,
double& drho, double& dtrtasc, double& dtdecl
)
{
	const double small = 0.00000001;
	const double omegaearth = 0.05883359221938136;  // earth rot rad/tu

	double earthrate[3], rhov[3], drhov[3], vsijk[3];
	double   latgc, temp, temp1;

	latgc = asin(rsijk[3] / astMath::mag(rsijk));
	earthrate[0] = 0.0;
	earthrate[1] = 0.0;
	earthrate[2] = omegaearth;
	astMath::cross(earthrate, rsijk, vsijk);

	if (direct == eFrom)
	{
		/* --------  calculate topocentric vectors ------------------ */
		rhov[0] = (rho * cos(tdecl) * cos(trtasc));
		rhov[1] = (rho * cos(tdecl) * sin(trtasc));
		rhov[2] = (rho * sin(tdecl));

		drhov[0] = (drho * cos(tdecl) * cos(trtasc) -
			rho * sin(tdecl) * cos(trtasc) * dtdecl -
			rho * cos(tdecl) * sin(trtasc) * dtrtasc);
		drhov[1] = (drho * cos(tdecl) * sin(trtasc) -
			rho * sin(tdecl) * sin(trtasc) * dtdecl +
			rho * cos(tdecl) * cos(trtasc) * dtrtasc);
		drhov[2] = (drho * sin(tdecl) + rho * cos(tdecl) * dtdecl);

		/* ------ find ijk range vector from site to satellite ------ */
		astMath::addvec(1.0, rhov, 1.0, rsijk, rijk);
		astMath::addvec(1.0, drhov, cos(latgc), vsijk, vijk);
		/*
		if (show == 'y')
		if (fileout != null)
		{
		fprintf(fileout, "rtb %18.7f %18.7f %18.7f %18.7f er\n",
		rhov[1], rhov[2], rhov[3], astMath::mag(rhov));
		fprintf(fileout, "vtb %18.7f %18.7f %18.7f %18.7f\n",
		drhov[1], drhov[2], drhov[3], astMath::mag(drhov));
		}
		*/
	}
	else
	{
		/* ------ find ijk range vector from site to satellite ------ */
		astMath::addvec(1.0, rijk, -1.0, rsijk, rhov);
		astMath::addvec(1.0, vijk, -cos(latgc), vsijk, drhov);

		/* -------- calculate topocentric angle and rate values ----- */
		rho = astMath::mag(rhov);
		temp = sqrt(rhov[0] * rhov[0] + rhov[1] * rhov[1]);
		if (temp < small)
		{
			temp1 = sqrt(drhov[0] * drhov[0] + drhov[1] * drhov[1]);
			trtasc = atan2(drhov[1] / temp1, drhov[0] / temp1);
		}
		else
			trtasc = atan2(rhov[1] / temp, rhov[0] / temp);

		tdecl = asin(rhov[2] / astMath::mag(rhov));

		temp1 = -rhov[1] * rhov[1] - rhov[0] * rhov[0];
		drho = astMath::dot(rhov, drhov) / rho;
		if (fabs(temp1) > small)
			dtrtasc = (drhov[0] * rhov[1] - drhov[1] * rhov[0]) / temp1;
		else
			dtrtasc = 0.0;
		if (fabs(temp) > small)
			dtdecl = (drhov[2] - drho * sin(tdecl)) / temp;
		else
			dtdecl = 0.0;
		/*
		if (show == 'y')
		if (fileout != null)
		{
		fprintf(fileout, "rta %18.7f %18.7f %18.7f %18.7f er\n",
		rhov[1], rhov[3], rhov[3], astMath::mag(rhov));
		fprintf(fileout, "vta %18.7f %18.7f %18.7f %18.7f er\n",
		drhov[1], drhov[3], drhov[3], astMath::mag(drhov));
		}
		*/
	}
} // procedure rv_tradec


/*------------------------------------------------------------------------------
*
*                           procedure rvsez_razel
*
*  this procedure converts range, azimuth, and elevation values with slant
*    range and velocity vectors for a satellite from a radar site in the
*    topocentric horizon (sez) system.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    rhovec      - sez satellite range vector     km
*    drhovec     - sez satellite velocity vector  km/s
*    direct      -  direction to convert          eFrom  eTo
*
*  outputs       :
*    rho         - satellite range from site      mk
*    az          - azimuth                        0.0 to 2pi rad
*    el          - elevation                      -pi/2 to pi/2 rad
*    drho        - range rate                     km/s
*    daz         - azimuth rate                   rad/s
*    del         - elevation rate                 rad/s
*
*  locals        :
*    sinel       - variable for sin( el )
*    cosel       - variable for cos( el )
*    sinaz       - variable for sin( az )
*    cosaz       - variable for cos( az )
*    temp        -
*    temp1       -
*
*  coupling      :
*    astMath::mag         - astMath::magnitude of a vector
*    sgn         - returns the sign of a variable
*    dot         - dot product
*    arcsin      - arc sine function
*    atan2       - arc tangent function that resolves quadrant ambiguites
*
*  references    :
*    vallado       2013, 261, eq 4-4, eq 4-5
-----------------------------------------------------------------------------*/

void rvsez_razel
(
double rhosez[3], double drhosez[3],
edirection direct,
double& rho, double& az, double& el, double& drho, double& daz, double& del
)
{
	const double small = 0.00000001;
	const double halfpi = pi / 2.0;

	double temp1, temp, sinel, cosel, sinaz, cosaz;

	if (direct == eFrom)
	{
		sinel = sin(el);
		cosel = cos(el);
		sinaz = sin(az);
		cosaz = cos(az);

		/* ----------------- form sez range vector ------------------ */
		rhosez[0] = (-rho * cosel * cosaz);
		rhosez[1] = (rho * cosel * sinaz);
		rhosez[2] = (rho * sinel);

		/* --------------- form sez velocity vector ----------------- */
		drhosez[0] = (-drho * cosel * cosaz +
			rhosez[2] * del * cosaz + rhosez[1] * daz);
		drhosez[1] = (drho * cosel * sinaz -
			rhosez[2] * del * sinaz - rhosez[0] * daz);
		drhosez[2] = (drho * sinel + rho * del * cosel);
	}
	else
	{
		/* ------------ calculate azimuth and elevation ------------- */
		temp = sqrt(rhosez[0] * rhosez[0] + rhosez[1] * rhosez[1]);
		if (fabs(rhosez[1]) < small)
		if (temp < small)
		{
			temp1 = sqrt(drhosez[0] * drhosez[0] +
				drhosez[1] * drhosez[1]);
			az = atan2(drhosez[1] / temp1, drhosez[0] / temp1);
		}
		else
		if (drhosez[0]  > 0.0)
			az = pi;
		else
			az = 0.0;
		else
			az = atan2(rhosez[1] / temp, rhosez[0] / temp);

		if (temp < small)   // directly over the north pole
			el = astMath::sgn(rhosez[2]) * halfpi;  // +- 90
		else
			el = asin(rhosez[2] / astMath::mag(rhosez));

		/* ------  calculate range, azimuth and elevation rates ----- */
		drho = astMath::dot(rhosez, drhosez) / rho;
		if (fabs(temp * temp) > small)
			daz = (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) /
			(temp * temp);
		else
			daz = 0.0;

		if (fabs(temp) > small)
			del = (drhosez[2] - drho * sin(el)) / temp;
		else
			del = 0.0;
	}
}   // procedure rvsez_razel


/* ------------------------- three vector techniques ------------------------ */

/* -----------------------------------------------------------------------------
*
*                           procedure gibbs
*
*  this procedure performs the gibbs method of orbit determination.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    r1          - ijk position vector #1         km
*    r2          - ijk position vector #2         km
*    r3          - ijk position vector #3         km
*
*  outputs       :
*    v2          - ijk velocity vector for r2     km / s
*    theta       - angle between vectors          rad
*    error       - flag indicating success        'ok',...
*
*  locals        :
*    tover2      -
*    l           -
*    small       - tolerance for roundoff errors
*    r1mr2       - astMath::magnitude of r1 - r2
*    r3mr1       - astMath::magnitude of r3 - r1
*    r2mr3       - astMath::magnitude of r2 - r3
*    p           - p vector     r2 x r3
*    q           - q vector     r3 x r1
*    w           - w vector     r1 x r2
*    d           - d vector     p + q + w
*    n           - n vector (r1)p + (r2)q + (r3)w
*    s           - s vector
*                    (r2-r3)r1+(r3-r1)r2+(r1-r2)r3
*    b           - b vector     d x r2
*    theta1      - temp angle between the vectors rad
*    pn          - p unit vector
*    r1n         - r1 unit vector
*    dn          - d unit vector
*    nn          - n unit vector
*    i           - index
*
*  coupling      :
*    astMath::mag         - astMath::magnitude of a vector
*    cross       - cross product of two vectors
*    dot         - dot product of two vectors
*    add3vec     - add three vectors
*    lncom2      - multiply two vectors by two constants
*    lncom3      - add three vectors each multiplied by a constant
*    norm        - creates a unit vector
*    angle       - angle between two vectors
*
*  references    :
*    vallado       2013, 460, alg 54, ex 7-3
-----------------------------------------------------------------------------*/

void gibbs
(
double r1[3], double r2[3], double r3[3],
double v2[3], double& theta, double& theta1, double& copa, char error[12]
)
{
	const double small = 0.000001;
	double tover2, l, r1mr2, r3mr1, r2mr3, mu;
	double p[3], q[3], w[3], d[3], n[3], s[3], b[3], pn[3], r1n[3], dn[3], nn[3];

	/* --------------------  initialize values   -------------------- */
	mu = 3.986004418e5;
#ifdef _MSC_VER
	strcpy_s(error, sizeof(error), "ok");
#else
	strcpy(error, "ok");
#endif
	
	theta = 0.0;
	theta1 = 0.0;
	copa = 0.0;

	/* ----------------------------------------------------------------
	*  determine if the vectors are coplanar.
	----------------------------------------------------------------- */
	astMath::cross(r2, r3, p);
	astMath::cross(r3, r1, q);
	astMath::cross(r1, r2, w);
	astMath::norm(p, pn);
	astMath::norm(r1, r1n);
	copa = asin(astMath::dot(pn, r1n));

	if (fabs(astMath::dot(r1n, pn)) > 0.017452406)
	{
#ifdef _MSC_VER
		strcpy_s(error, sizeof(error), "not coplanar");
#else
		strcpy(error, "not coplanar");
#endif
	}

	/* ---------------- or don't continue processing ---------------- */
	astMath::addvec3(1.0, p, 1.0, q, 1.0, w, d);
	astMath::addvec3(astMath::mag(r1), p, astMath::mag(r2), q, astMath::mag(r3), w, n);
	astMath::norm(n, nn);
	astMath::norm(d, dn);

	/* ----------------------------------------------------------------
	*  determine if the orbit is possible.  both d and n must be in
	*    the same direction, and non-zero.
	----------------------------------------------------------------- */
	if ((fabs(astMath::mag(d)) < small) || (fabs(astMath::mag(n)) < small) ||
		(fabs(astMath::dot(nn, dn)) < small))
	{
#ifdef _MSC_VER
		strcpy_s(error, sizeof(error), "impossible");
#else
		strcpy(error, "impossible");
#endif
	}
	else
	{
		theta = astMath::angle(r1, r2);
		theta1 = astMath::angle(r2, r3);

		/* ------------ perform gibbs method to find v2 ----------- */
		r1mr2 = astMath::mag(r1) - astMath::mag(r2);
		r3mr1 = astMath::mag(r3) - astMath::mag(r1);
		r2mr3 = astMath::mag(r2) - astMath::mag(r3);
		astMath::addvec3(r1mr2, r3, r3mr1, r2, r2mr3, r1, s);
		astMath::cross(d, r2, b);
		l = sqrt(mu / (astMath::mag(d) * astMath::mag(n)));
		tover2 = l / astMath::mag(r2);
		astMath::addvec(tover2, b, l, s, v2);
	}
	/*
	if (((show == 'y') || (show == 's')) && (strcmp(error, "ok") == 0))
	if (fileout != null)
	{
	fprintf(fileout, "%16s %9.3f %9.3f %9.3f\n",
	"p vector = ", p[1], p[2], p[3]);
	fprintf(fileout, "%16s %9.3f %9.3f %9.3f\n",
	"q vector = ", q[1], q[2], q[3]);
	fprintf(fileout, "%16s %9.3f %9.3f %9.3f\n",
	"w vector = ", w[1], w[2], w[3]);
	fprintf(fileout, "%16s %9.3f %9.3f %9.3f\n",
	"n vector = ", n[1], n[2], n[3]);
	fprintf(fileout, "%16s %9.3f %9.3f %9.3f\n",
	"d vector = ", d[1], d[2], d[3]);
	fprintf(fileout, "%16s %9.3f %9.3f %9.3f\n",
	"s vector = ", s[1], s[2], s[3]);
	fprintf(fileout, "%16s %9.3f %9.3f %9.3f\n",
	"b vector = ", b[1], b[2], b[3]);
	}
	*/
}  // procedure gibbs


/*------------------------------------------------------------------------------
*
*                           procedure herrgibbs
*
*  this procedure implements the herrick-gibbs approximation for orbit
*    determination, and finds the middle velocity vector for the 3 given
*    position vectors.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    r1          - ijk position vector #1         km
*    r2          - ijk position vector #2         km
*    r3          - ijk position vector #3         km
*    jd1         - julian date of 1st sighting    days from 4713 bc
*    jd2         - julian date of 2nd sighting    days from 4713 bc
*    jd3         - julian date of 3rd sighting    days from 4713 bc
*
*  outputs       :
*    v2          - ijk velocity vector for r2     km / s
*    theta       - angle between vectors          rad
*    error       - flag indicating success        'ok',...
*
*  locals        :
*    dt21        - time delta between r1 and r2   s
*    dt31        - time delta between r3 and r1   s
*    dt32        - time delta between r3 and r2   s
*    p           - p vector    r2 x r3
*    pn          - p unit vector
*    r1n         - r1 unit vector
*    theta1      - temporary angle between vec    rad
*    tolangle    - tolerance angle  (1 deg)       rad
*    term1       - 1st term for hgibbs expansion
*    term2       - 2nd term for hgibbs expansion
*    term3       - 3rd term for hgibbs expansion
*    i           - index
*
*  coupling      :
*    astMath::mag         - astMath::magnitude of a vector
*    cross       - cross product of two vectors
*    dot         - dot product of two vectors
*    arcsin      - arc sine function
*    norm        - creates a unit vector
*    lncom3      - combination of three scalars and three vectors
*    angle       - angle between two vectors
*
*  references    :
*    vallado       2013, 466, alg 55, ex 7-4
-----------------------------------------------------------------------------*/

void herrgibbs
(
double r1[3], double r2[3], double r3[3], double jd1, double jd2, double jd3,
double v2[3], double& theta, double& theta1, double& copa, char error[12]
)
{
	double p[3], pn[3], r1n[3], mu;
	double   dt21, dt31, dt32, term1, term2, term3, tolangle;

	/* --------------------  initialize values   -------------------- */
	mu = 3.986004418e5;
	tolangle = 0.017452406;  // (1.0 deg in rad)

#ifdef _MSC_VER
		strcpy_s(error, sizeof(error), "ok");
#else
		strcpy(error, "ok");
#endif

	theta = 0.0;
	theta1 = 0.0;

	dt21 = (jd2 - jd1) * 86400.0;
	dt31 = (jd3 - jd1) * 86400.0;   // differences in times
	dt32 = (jd3 - jd2) * 86400.0;

	/* ----------------------------------------------------------------
	*  determine if the vectors are coplanar.
	---------------------------------------------------------------- */
	astMath::cross(r2, r3, p);
	astMath::norm(p, pn);
	astMath::norm(r1, r1n);
	copa = asin(astMath::dot(pn, r1n));
	if (fabs(astMath::dot(pn, r1n)) > tolangle)
	{
#ifdef _MSC_VER
		strcpy_s(error, sizeof(error), "not coplanar");
#else
		strcpy(error, "not coplanar");
#endif
	}

	/* ----------------------------------------------------------------
	* check the size of the angles between the three position vectors.
	*   herrick gibbs only gives "reasonable" answers when the
	*   position vectors are reasonably close.  1.0 deg is only an estimate.
	---------------------------------------------------------------- */
	theta = astMath::angle(r1, r2);
	theta1 = astMath::angle(r2, r3);
	if ((theta > tolangle) || (theta1 > tolangle))
	{
#ifdef _MSC_VER
		strcpy_s(error, sizeof(error), "angle > 1ø");
#else
		strcpy(error, "angle > 1ø");
#endif
	}

	/* ------------ perform herrick-gibbs method to find v2 --------- */
	term1 = -dt32 *
		(1.0 / (dt21 * dt31) + mu / (12.0 * astMath::mag(r1) * astMath::mag(r1) * astMath::mag(r1)));
	term2 = (dt32 - dt21) *
		(1.0 / (dt21 * dt32) + mu / (12.0 * astMath::mag(r2) * astMath::mag(r2) * astMath::mag(r2)));
	term3 = dt21 *
		(1.0 / (dt32 * dt31) + mu / (12.0 * astMath::mag(r3) * astMath::mag(r3) * astMath::mag(r3)));
	astMath::addvec3(term1, r1, term2, r2, term3, r3, v2);
}  // procedure herrgibbs


/* ----------------------- lambert techniques -------------------- */

/*------------------------------------------------------------------------------
*
*                           procedure lamberbattin
*
*  this procedure solves lambert's problem using battins method. the method is
*    developed in battin (1987).
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    ro          - ijk position vector 1          km
*    r           - ijk position vector 2          km
*   dm          - direction of motion            'l','s'
*    dtsec        - time between r1 and r2         sec
*
*  outputs       :
*    vo          - ijk velocity vector            er / tu
*    v           - ijk velocity vector            er / tu
*    error       - error flag                     1, 2, 3, ... use numbers since c++ is so horrible at strings
*        error = 1;   // a = 0.0
*
*  locals        :
*    i           - index
*    loops       -
*    u           -
*    b           -
*    sinv        -
*    cosv        -
*    rp          -
*    x           -
*    xn          -
*    y           -
*    l           -
*    m           -
*    cosdeltanu  -
*    sindeltanu  -
*    dnu         -
*    a           -
*    tan2w       -
*    ror         -
*    h1          -
*    h2          -
*    tempx       -
*    eps         -
*    denom       -
*    chord       -
*    k2          -
*    s           -
*    f           -
*    g           -
*    fdot        -
*    am          -
*    ae          -
*    be          -
*    tm          -
*    gdot        -
*    arg1        -
*    arg2        -
*    tc          -
*    alpe        -
*    bete        -
*    alph        -
*    beth        -
*    de          -
*    dh          -
*
*  coupling      :
*    arcsin      - arc sine function
*    arccos      - arc cosine function
*    astMath::mag         - astMath::magnitude of a vector
*    arcsinh     - inverse hyperbolic sine
*    arccosh     - inverse hyperbolic cosine
*    sinh        - hyperbolic sine
*    power       - raise a base to a power
*    atan2       - arc tangent function that resolves quadrants
*
*  references    :
*    vallado       2013, 494, Alg 59, ex 7-5
-----------------------------------------------------------------------------*/

void lambertbattin
(
double ro[3], double r[3], char dm, int overrev, double dtsec,
double vo[3], double v[3], int& error
)
{
	const double small = 0.000001;
	const double mu = 398600.4418;  // m3s2
	double rcrossr[3];
	int   i, loops;
	double   u, b, sinv, cosv, rp, x, xn, y, l, m, cosdeltanu, sindeltanu, dnu, a,
		tan2w, ror, h1, h2, tempx, eps, denom, chord, k2, s, f, g, fdot, am,
		ae, be, tm, gdot, arg1, arg2, alpe, bete, alph, beth, de, dh,
		testamt, bigt;

	error = 0;
	cosdeltanu = astMath::dot(ro, r) / (astMath::mag(ro) * astMath::mag(r));
	astMath::cross(ro, r, rcrossr);
	sindeltanu = astMath::mag(rcrossr) / (astMath::mag(ro) * astMath::mag(r));
	dnu = atan2(sindeltanu, cosdeltanu);

	ror = astMath::mag(r) / astMath::mag(ro);
	eps = ror - 1.0;
	tan2w = 0.25 * eps * eps / (sqrt(ror) + ror *(2.0 + sqrt(ror)));
	rp = sqrt(astMath::mag(ro)*astMath::mag(r)) * ((cos(dnu * 0.25) * cos(dnu * 0.25)) + tan2w);

	  if (dnu < pi)
	    l = (pow(sin(dnu * 0.25), 2) + tan2w ) /
	        (pow(sin(dnu * 0.25), 2) + tan2w + cos(dnu * 0.5));
	  else
	    l = (pow(cos(dnu * 0.25), 2) + tan2w - cos(dnu * 0.5)) /
	        (pow(cos(dnu * 0.25), 2) + tan2w);

	m = mu * dtsec * dtsec / (8.0 * rp * rp * rp);
	xn = 0.0;   // 0 for par and hyp
	chord = sqrt(astMath::mag(ro) * astMath::mag(ro) + astMath::mag(r) * astMath::mag(r) -
		2.0 * astMath::mag(ro) * astMath::mag(r) * cos(dnu));
	s = (astMath::mag(ro) + astMath::mag(r) + chord) * 0.5;

	loops = 1;
	while (1 == 1)
	{
		x = xn;
		tempx = see(x);
		denom = 1.0 / ((1.0 + 2.0 * x + l) * (3.0 + x * (1.0 + 4.0 * tempx)));
		h1    = pow(l + x, 2) * ( 1.0 + (1.0 + 3.0 * x) * tempx) * denom;
		h2 = m * (1.0 + (x - l) * tempx) * denom;

		/* ------------------------ evaluate cubic ------------------ */
		b  = 0.25 * 27.0 * h2 / pow(1.0 + h1, 3);
		u = -0.5 * b / (1.0 + sqrt(1.0 + b));
		k2 = k(u);

		y = ((1.0 + h1) / 3.0) * (2.0 + sqrt(1.0 + b) /
			(1.0 - 2.0 * u * k2 * k2));
		xn = sqrt(pow((1.0 - l) * 0.5, 2) + m / (y * y)) - (1.0 + l) * 0.5;

		loops++;

		if ((fabs(xn - x) < small) || (loops > 30))
			break;
	}

	a = mu * dtsec * dtsec / (16.0 * rp * rp * xn * y * y);

	/* -------------------- find eccentric anomalies ---------------- */
	/* -------------------------- hyperbolic ------------------------ */
	if (a < -small)
	{
		arg1 = sqrt(s / (-2.0 * a));
		arg2 = sqrt((s - chord) / (-2.0 * a));
		/* -------- evaluate f and g functions -------- */
		alph = 2.0 * astMath::asinh(arg1);
		beth = 2.0 * astMath::asinh(arg2);
		dh = alph - beth;
		f = 1.0 - (a / astMath::mag(ro)) * (1.0 - cosh(dh));
		gdot = 1.0 - (a / astMath::mag(r)) * (1.0 - cosh(dh));
        g    = dtsec - sqrt(-a*a*a/mu)*(sinh(dh)-dh);
		fdot = -sqrt(-a) * sinh(dh) / (astMath::mag(ro) * astMath::mag(r));
	}
	else
	{
		/* ------------------------- elliptical --------------------- */
		if (a > small)
		{
			arg1 = sqrt(s / (2.0 * a));
			arg2 = sqrt((s - chord) / (2.0 * a));
			sinv = arg2;
			cosv = sqrt(1.0 - (astMath::mag(ro) + astMath::mag(r) - chord) / (4.0 * a));
			bete = 2.0 * acos(cosv);
			bete = 2.0 * asin(sinv);
			if (dnu > pi)
				bete = -bete;

			cosv = sqrt(1.0 - s / (2.0 * a));
			sinv = arg1;

			am = s * 0.5;
			ae = pi;
			be = 2.0 * asin(sqrt((s - chord) / s));
			tm = sqrt(am * am * am/mu) * (ae - (be - sin(be)));
			if (dtsec > tm)
				alpe = 2.0 * pi - 2.0 * asin(sinv);
			else
				alpe = 2.0 * asin(sinv);
			de = alpe - bete;
			f = 1.0 - (a / astMath::mag(ro)) * (1.0 - cos(de));
			gdot = 1.0 - (a / astMath::mag(r)) * (1.0 - cos(de));
			g = dtsec - sqrt(a * a * a/mu) * (de - sin(de));
			fdot = -sqrt(a) * sin(de) / (astMath::mag(ro) * astMath::mag(r));
		}
		else
		{
			/* ------------------------- parabolic -------------------- */
			arg1 = 0.0;
			arg2 = 0.0;
            f = 0.0;
            g = 0.0;
            gdot = 0.0;
            fdot = 0.0;
			error = 1; // a = 0.0
			//      if (fileout != null)
			//        fprintf(fileout, " a parabolic orbit \n");
		}
	}

	for (i = 1; i <= 3; i++)
	{
		vo[i] = ((r[i] - f * ro[i]) / g);
		v[i] = ((gdot * r[i] - ro[i]) / g);
	}

	if (error == 0)
		testamt = f * gdot - fdot * g;
	else
		testamt = 2.0;

	//  if (fileout != null)
	//    fprintf(fileout, "%8.5f %3d\n", testamt, loops);

	bigt = sqrt(8.0 / (s * s * s)) * dtsec;
}  // procedure lambertbattin




/*------------------------------------------------------------------------------
*
*                           procedure lambertuniv
*
*  this procedure solves the lambert problem for orbit determination and returns
*    the velocity vectors at each of two given position vectors.  the solution
*    uses universal variables for calculation and a bissection technique for
*    updating psi.
*
*  algorithm     : setting the initial bounds:
*                  using -8pi and 4pi2 will allow single rev solutions
*                  using -4pi2 and 8pi2 will allow multi-rev solutions
*                  the farther apart the initial guess, the more iterations
*                    because of the iteration
*                  inner loop is for special cases. must be sure to exit both!
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    r1          - ijk position vector 1          er
*    r2          - ijk position vector 2          er
*    dm          - direction of motion            'l','s'
*    dtsec       - time between r1 and r2         sec
*
*  outputs       :
*    v1          - ijk velocity vector            er / tu
*    v2          - ijk velocity vector            er / tu
*    error       - error flag                     1, 2, 3, ... use numbers since c++ is so horrible at strings
*        error = 1;   // g not converged
*        error = 2;   // y negative
*        error = 3;   // impossible 180 transfer
*
*  locals        :
*    vara        - variable of the iteration,
*                  not the semi or axis!
*    y           - area between position vectors
*    upper       - upper bound for z
*    lower       - lower bound for z
*    cosdeltanu  - cosine of true anomaly change  rad
*    f           - f expression
*    g           - g expression
*    gdot        - g dot expression
*    xold        - old universal variable x
*    xoldcubed   - xold cubed
*    zold        - old value of z
*    znew        - new value of z
*    c2new       - c2(z) function
*    c3new       - c3(z) function
*    timenew     - new time                       tu
*    small       - tolerance for roundoff errors
*    i, j        - index
*
*  coupling
*    astMath::mag         - astMath::magnitude of a vector
*    dot         - dot product of two vectors
*    findc2c3    - find c2 and c3 functions
*
*  references    :
*    vallado       2013, 492, alg 58, ex 7-5
-----------------------------------------------------------------------------*/

void lambertuniv
(
double ro[3], double r[3], char dm, int nrev, double dtsec,
double vo[3], double v[3], int& error, FILE *outfile
)
{
	const double twopi = 2.0 * pi;
	const double small = 0.0000001;
	const int numiter = 40;
	const double mu = 398600.4418;  // m3s2

	int loops, ynegktr;
	double vara, y, upper, lower, cosdeltanu, f, g, gdot, xold, xoldcubed, 
		psiold, psinew, c2new, c3new, dtnew, c2dot, c3dot, dtdpsi, psiold2;

	/* --------------------  initialize values   -------------------- */
	error = 0;
	psinew = 0.0;

	/* -------- set up initial bounds for the bissection ------------ */
	if (nrev == 0)
	{
		upper = 4.0 * pi * pi;  // could be negative infinity for all cases
		lower = -4.0 * pi * pi; // allow hyperbolic and parabolic solutions
	}
	else
	{
		lower = 4.0 * nrev * nrev * pi * pi;
		upper = 4.0 * (nrev + 1.0) * (nrev + 1.0) * pi * pi;
	}

	/* ----------------  form initial guesses   --------------------- */
	psinew = 0.0;
	xold = 0.0;
	if (nrev == 0)
		psiold = (log(dtsec / 806.811874) - 1.2357) / 0.118;
	else
	{
		if (dm == 's')
			psiold = lower + (upper - lower)*0.75;
		else
			psiold = lower + (upper - lower)*0.25;
	}
	ast2Body::findc2c3(psiold, c2new, c3new);

	cosdeltanu = astMath::dot(ro, r) / (astMath::mag(ro) * astMath::mag(r));
	if (dm == 'l')
		vara = -sqrt(astMath::mag(ro) * astMath::mag(r) * (1.0 + cosdeltanu));
	else
		vara = sqrt(astMath::mag(ro) * astMath::mag(r) * (1.0 + cosdeltanu));


	/* --------  determine if the orbit is possible at all ---------- */
	if (fabs(vara) > small)
	{
		loops = 0;
		ynegktr = 1; // y neg ktr
		dtnew = -10.0; 
		while ((fabs(dtnew-dtsec) >= small) && (loops < numiter) && (ynegktr <= 10))
		{
			if (fabs(c2new) > small)
				y = astMath::mag(ro) + astMath::mag(r) - (vara * (1.0 - psiold * c3new) / sqrt(c2new));
			else
				y = astMath::mag(ro) + astMath::mag(r);
			/* ------- check for negative values of y ------- */
			if ((vara > 0.0) && (y < 0.0))
			{
				ynegktr = 1;
				while (( y < 0.0 ) && ( ynegktr < 10 ))
				{
					psinew = 0.8 * (1.0 / c3new) *
						(1.0 - (astMath::mag(ro) + astMath::mag(r)) * sqrt(c2new) / vara);

					/* ------ find c2 and c3 functions ------ */
					ast2Body::findc2c3(psinew, c2new, c3new);
					psiold = psinew;
					lower = psiold;
					if (fabs(c2new) > small)
						y = astMath::mag(ro) + astMath::mag(r) -
						(vara * (1.0 - psiold * c3new) / sqrt(c2new));
					else
						y = astMath::mag(ro) + astMath::mag(r);
					//          if (show == 'y')
					//            if (fileout != null)
					//              fprintf(fileout, "%3d %10.5f %10.5f %10.5f %7.3f %9.5f %9.5f\n",
					//                      loops, psiold, y, xold, dtnew, vara, upper, lower);

					ynegktr++;
				}
			}

			if (ynegktr < 10)
			{
				if (fabs(c2new) > small)
					xold = sqrt(y / c2new);
				else
					xold = 0.0;
				xoldcubed = xold * xold * xold;
				dtnew = (xoldcubed * c3new + vara * sqrt(y)) / sqrt(mu);

				// try newton rhapson iteration to update psi
				if (fabs(psiold) > 1e-5)
				{
					c2dot = 0.5 / psiold * (1.0 - psiold * c3new - 2.0 * c2new);
					c3dot = 0.5 / psiold * (c2new - 3.0 * c3new); 
				}
				else
				{
					psiold2 = psiold * psiold;
					c2dot = -1.0 / astMath::factorial(4) + 2.0 * psiold / astMath::factorial(6) - 3.0 * psiold2 / astMath::factorial(8) +
					    	 4.0 * psiold2 * psiold / astMath::factorial(10) - 5.0 * psiold2*psiold2 / astMath::factorial(12);
					c3dot = -1.0 / astMath::factorial(5) + 2.0 * psiold / astMath::factorial(7) - 3.0 * psiold2 / astMath::factorial(9) +
					     	 4.0 * psiold2 * psiold / astMath::factorial(11) - 5.0 * psiold2*psiold2 / astMath::factorial(13);
				}
				dtdpsi = (xoldcubed * (c3dot - 3.0 * c3new * c2dot / (2.0 * c2new)) + vara / 8.0 * (3.0 * c3new * sqrt(y) / c2new + vara / xold)) / sqrt(mu);
				psinew = psiold - (dtnew - dtsec) / dtdpsi;

				// check if newton guess for psi is outside bounds(too steep a slope)
				if (psinew > upper || psinew < lower)
				{
					// --------readjust upper and lower bounds------ -
					if (dtnew < dtsec)
						lower = upper - pi;
					if (dtnew > dtsec)
						upper = psiold;
					psinew = (upper + lower) * 0.5;
				}
					

				/* -------------- find c2 and c3 functions ---------- */
				ast2Body::findc2c3(psinew, c2new, c3new);
				psiold = psinew;
				loops = loops + 1;

				/* ---- make sure the first guess isn't too close --- */
				if ((fabs(dtnew - dtsec) < small) && (loops == 1))
					dtnew = dtsec - 1.0;
			}

		}

		if ((loops >= numiter) || (ynegktr >= 10))
		{
			error = 1; // g not converged

			if (ynegktr >= 10)
			{
				error = 2;  // y negative
			}
		}
		else
		{
			/* ---- use f and g series to find velocity vectors ----- */
			f = 1.0 - y / astMath::mag(ro);
			gdot = 1.0 - y / astMath::mag(r);
			g = 1.0 / (vara * sqrt(y/mu)); // 1 over g
		//	fdot = sqrt(y) * (-astMath::mag(r) - astMath::mag(ro) + y) / (astMath::mag(r) * astMath::mag(ro) * vara);
			for (int i = 0; i <= 3; i++)
			{
				vo[i] = ((r[i] - f * ro[i]) * g);
				v[i] = ((gdot * r[i] - ro[i]) * g);
			}
		}
	}
	else
		error = 3;   // impossible 180 transfer

		if (dtsec < 0.001)
			fprintf(outfile, " \n");
		else
		    fprintf(outfile, "%5i %c %5i %12.5f %12.5f %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n",
		                     nrev, dm, loops, dtsec, dtsec, y, xold, vo[0], vo[1], vo[2], v[0], v[1], v[2], lower, upper, psinew);

	/*
	---- for fig 6-14 dev with testgau.pas ----
	if error = 'ok' then write( fileout,psinew:14:7,dtsec*13.44685:14:7 )
	else write( fileout,' 9999.0 ':14,dtsec*13.44685:14:7 );
	*/
}  // procedure lambertuniv


/*------------------------------------------------------------------------------
*
*                           procedure target
*
*  this procedure accomplishes the targeting problem using kepler/pkepler and
*    lambert.
*
*  author        : david vallado                  719-573-2600   22 jun 2002
*
*  inputs          description                    range / units
*    rint        - initial position vector of int er
*    vint        - initial velocity vector of int er/tu
*    rtgt        - initial position vector of tgt er
*    vtgt        - initial velocity vector of tgt er/tu
*    dm          - direction of motion for gauss  'l','s'
*    kind        - type of propagator             'k','p'
*    dtsec        - time of flight to the int      tu
*
*  outputs       :
*    v1t         - initial transfer velocity vec  er/tu
*    v2t         - final transfer velocity vec    er/tu
*    dv1         - initial change velocity vec    er/tu
*    dv2         - final change velocity vec      er/tu
*    error       - error flag from gauss          'ok', ...
*
*  locals        :
*    transnormal - cross product of trans orbit   er
*    intnormal   - cross product of int orbit     er
*    r1tgt       - position vector after dt, tgt  er
*    v1tgt       - velocity vector after dt, tgt  er/tu
*    rirt        - rint[4] * r1tgt[4]
*    cosdeltanu  - cosine of deltanu              rad
*    sindeltanu  - sine of deltanu                rad
*    deltanu     - angle between position vectors rad
*    i           - index
*
*  coupling      :
*    kepler      - find r and v at future time
*    lambertuniv - find velocity vectors at each end of transfer
*    lncom2      - linear combination of two vectors and constants
*
*  references    :
*    vallado       2013, 503, alg 61
-----------------------------------------------------------------------------*/

void target
(
double rint[3], double vint[3], double rtgt[3], double vtgt[3],
char dm, char kind, double dtsec,
double v1t[3], double v2t[3], double dv1[3], double dv2[3], int error
)
{
	double  r1tgt[3], v1tgt[3];
	FILE *outfile;

	/* ----------- propagate target forward in time ----------------- */
	switch (kind)
	{
	case 'k':
		//      kepler(rtgt, vtgt, dtsec,  r1tgt, v1tgt, error);
		break;
	case 'p':
		//      pkepler(rtgt, vtgt, dtsec,  r1tgt, v1tgt, error);
		break;
	default:
		//      kepler(rtgt, vtgt, dtsec,  r1tgt, v1tgt, error);
		break;
	}

	/* ----------- calculate transfer orbit between r's ------------- */
	if (error == 0)
	{
		lambertuniv(rint, r1tgt, dm, 'n', dtsec, v1t, v2t, error, outfile);

		if (error == 0)
		{
			astMath::addvec(1.0, v1t, -1.0, vint, dv1);
			astMath::addvec(1.0, v1tgt, -1.0, v2t, dv2);
		}
		else
		{
		}
	}
}  // procedure target


/* utility functions for lambertbattin, etc */
/* -------------------------------------------------------------------------- */

static double k
(
double v
)
{
	double d[21] =
	{
		1.0 / 3.0, 4.0 / 27.0,
		8.0 / 27.0, 2.0 / 9.0,
		22.0 / 81.0, 208.0 / 891.0,
		340.0 / 1287.0, 418.0 / 1755.0,
		598.0 / 2295.0, 700.0 / 2907.0,
		928.0 / 3591.0, 1054.0 / 4347.0,
		1330.0 / 5175.0, 1480.0 / 6075.0,
		1804.0 / 7047.0, 1978.0 / 8091.0,
		2350.0 / 9207.0, 2548.0 / 10395.0,
		2968.0 / 11655.0, 3190.0 / 12987.0,
		3658.0 / 14391.0
	};
	double del, delold, term, termold, sum1;
	int i;

	/* ---- process forwards ---- */
	sum1 = d[0];
	delold = 1.0;
	termold = d[0];
	i = 1;
	while ((i <= 20) && (fabs(termold) > 0.000001))
	{
		del = 1.0 / (1.0 - d[i] * v * delold);
		term = termold * (del - 1.0);
		sum1 = sum1 + term;
		i++;
		delold = del;
		termold = term;
	}
	return sum1;
}  // double k


/* -------------------------------------------------------------------------- */

static double see(double v)
{
	double c[21] =
	{
		0.2,
		9.0 / 35.0, 16.0 / 63.0,
		25.0 / 99.0, 36.0 / 143.0,
		49.0 / 195.0, 64.0 / 255.0,
		81.0 / 323.0, 100.0 / 399.0,
		121.0 / 483.0, 144.0 / 575.0,
		169.0 / 675.0, 196.0 / 783.0,
		225.0 / 899.0, 256.0 / 1023.0,
		289.0 / 1155.0, 324.0 / 1295.0,
		361.0 / 1443.0, 400.0 / 1599.0,
		441.0 / 1763.0, 484.0 / 1935.0
	};
	double term, termold, del, delold, sum1, eta, sqrtopv;
	int i;

	sqrtopv = sqrt(1.0 + v);
	eta     = v / pow(1.0 + sqrtopv, 2);

	/* ---- process forwards ---- */
	delold = 1.0;
	termold = c[0];  // * eta
	sum1 = termold;
	i = 1;
	while ((i <= 20) && (fabs(termold) > 0.000001))
	{
		del = 1.0 / (1.0 + c[i] * eta * delold);
		term = termold * (del - 1.0);
		sum1 = sum1 + term;
		i++;
		delold = del;
		termold = term;
	}

	return ((1.0 / (8.0 * (1.0 + sqrtopv))) * (3.0 + sum1 / (1.0 + eta * sum1)));
}  // double see

}   // namespace