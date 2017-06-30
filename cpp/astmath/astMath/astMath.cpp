/*     ----------------------------------------------------------------
*
*                                astMath.cpp
*
*    this file contains miscellaneous math functions.
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
*               7 may 08  david vallado
*                           misc updates, fix sgn, show both options for matrix
*                           multiplications
*              22 jan 08  david vallado
*                           fix some minor errors, fixes to matmult
*              19 oct 07  david vallado
*                           fix sgn baseline comparison
*              30 may 07  david vallado
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

#include "astMath.h"

namespace astMath {

double sgn (
          double x
        )
   {
     if (x < 0.0)
       {
          return -1.0;
       }
       else
       {
          return 1.0;
       }

   }  // end sgn


// round a number to the nearest integer
double round
       (
         double x
       )
   {
       double temp;
       temp = floor( x + 0.5);
       return int(temp);
   }  // end round

/* -----------------------------------------------------------------------------
*
*                           function acosh
*
*  this function evaluates the inverse hyperbolic cosine function.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    xval        - angle value                                  1.0 to infinity
*
*  outputs       :
*    acosh       - result                                       any real
*
*  locals        :
*    temp        - temporary value
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

double  acosh
        (
          double xval
        )
   {
     double temp;
     if ( xval*xval - 1.0 < 0.0 )
       {
         temp = undefined;
         printf("error in arccosh function \n");
       } 
       else
         temp= log( xval + sqrt( xval*xval - 1.0 ) );

     return temp;
   }


/* -----------------------------------------------------------------------------
*
*                           function asinh
*
*  this function evaluates the inverse hyperbolic sine function.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    xval        - angle value                                  any real
*
*  outputs       :
*    asinh       - result                                       any real
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

double  asinh
        (
          double xval
        )
   {
     return log( xval + sqrt( xval*xval + 1.0 ) );
   }  // end asinh


/* ------------------------------------------------------------------------------
*
*                           function cot
*
*  this function finds the cotangent of an input in radians.
*
*  author        : david vallado                  719-573-2600    1 Mar 2001
*
*  inputs          description                    range / units
*    xval        - input to take cotangent of        rad
*
*  outputs       :
*    cot         - result
*
*  locals        :
*    temp        - temporary real variable
 ---------------------------------------------------------------------------- */

double cot
       (
         double xval
       )
   {
      double temp;

      temp = tan( xval );
      if (fabs( temp ) < 0.00000001 )
          return infinite;
        else
          return 1.0 / temp;
   }  // procedure cot


/* -----------------------------------------------------------------------------
*
*                           function dot
*
*  this function finds the dot product of two vectors.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    dot         - result
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

double  dot
        (
          double x[3], double y[3]
        )
   {
     return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
   }  // end dot

/* -----------------------------------------------------------------------------
*
*                           function mag
*
*  this procedure finds the magnitude of a vector.  the tolerance is set to
*    0.000001, thus the 1.0e-12 for the squared test of underflows.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec       - vector
*
*  outputs       :
*    vec       - answer stored in function return
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

double  mag
        (
          double x[3]
        )
   {
     return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
   }  // end mag

/* -----------------------------------------------------------------------------
*
*                           procedure cross
*
*  this procedure crosses two vectors.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    outvec      - vector result of a x b
*
*  locals        :
*    none.
*
*  coupling      :
*    none
 ---------------------------------------------------------------------------- */

void    cross
        (
          double vec1[3], double vec2[3], double outvec[3]
        )
   {
     outvec[0]= vec1[1]*vec2[2] - vec1[2]*vec2[1];
     outvec[1]= vec1[2]*vec2[0] - vec1[0]*vec2[2];
     outvec[2]= vec1[0]*vec2[1] - vec1[1]*vec2[0];
   }  // end cross

/* -----------------------------------------------------------------------------
*
*                           procedure norm
*
*  this procedure calculates a unit vector given the original vector.  if a
*    zero vector is input, the vector is set to zero.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec         - vector
*
*  outputs       :
*    outvec      - unit vector
*
*  locals        :
*    i           - index
*    small       - value defining a small value
*    magv        - magnitude of the vector
*
*  coupling      :
*    mag           magnitude of a vector
* --------------------------------------------------------------------------- */

void    norm
        (
          double vec[3],
          double outvec[3]
        )
   {
     const double small = 0.000001;
     double magv;
     int i;

     magv = mag( vec );
     if (magv > small)
       {
         for (i = 0; i <= 2; i++)
             outvec[i]= vec[i]/magv;
       }
       else
         for (i = 0; i <= 2; i++)
             outvec[i]= 0.0;
   }

/* -----------------------------------------------------------------------------
*
*                           procedure roti
*
*  this procedure performs a rotation about the ith axis. i is specified
*    for each operation.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec         - input vector
*    xval        - angle of rotation              rad
*
*  outputs       :
*    outvec      - vector result
*
*  locals        :
*    c           - cosine of the angle xval
*    s           - sine of the angle xval
*    temp        - temporary extended value
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

void    rot1
        (
          double vec[3],
          double xval,
          double outvec[3]
        )
   {
      double c, s, temp;

      temp= vec[2];
      c= cos( xval );
      s= sin( xval );

      outvec[2]= c*vec[2] - s*vec[1];
      outvec[1]= c*vec[1] + s*temp;
      outvec[0]= vec[0];
   }

void    rot2
        (
          double vec[3],
          double xval,
          double outvec[3]
        )
   {
      double c, s, temp;

      temp= vec[2];
      c= cos( xval );
      s= sin( xval );

      outvec[2]= c*vec[2] + s*vec[0];
      outvec[0]= c*vec[0] - s*temp;
      outvec[1]= vec[1];
   }

void    rot3
        (
          double vec[3],
          double xval,
          double outvec[3]
        )
   {
      double c, s, temp;

      temp= vec[1];
      c= cos( xval );
      s= sin( xval );

      outvec[1]= c*vec[1] - s*vec[0];
      outvec[0]= c*vec[0] + s*temp;
      outvec[2]= vec[2];
   }

/* -----------------------------------------------------------------------------
*
*                                  rotimat
*
*  this function sets up a rotation matrix for an input angle about the first
*    axis.
*
*  author        : david vallado                  719-573-2600   10 jan 2003
*
*  revisions
*                -
*
*  inputs          description                    range / units
*    xval        - angle of rotation              rad
*
*  outputs       :
*    outmat      - matrix result
*
*  locals        :
*    c           - cosine of the angle xval
*    s           - sine of the angle xval
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

void    rot1mat
        (
          double xval,
          std::vector< std::vector<double> > &outmat
//          double& outmat[3][3]
        )
   {
     outmat.resize(3);  // rows
     for (std::vector< std::vector<double> >::iterator it=outmat.begin(); it != outmat.end();++it)
          it->resize(3);
     double c, s;
     c= cos( xval );
     s= sin( xval );

     outmat[0][0]= 1.0;
     outmat[0][1]= 0.0;
     outmat[0][2]= 0.0;

     outmat[1][0]= 0.0;
     outmat[1][1]= c;
     outmat[1][2]= s;

     outmat[2][0]= 0.0;
     outmat[2][1]= -s;
     outmat[2][2]= c;
   }

void    rot2mat
        (
          double xval,
          std::vector< std::vector<double> > &outmat
//          double outmat[3][3]
        )
   {
     outmat.resize(3);  // rows
     for (std::vector< std::vector<double> >::iterator it=outmat.begin(); it != outmat.end();++it)
          it->resize(3);
     double c, s;
     c= cos( xval );
     s= sin( xval );

     outmat[0][0]= c;
     outmat[0][1]= 0.0;
     outmat[0][2]= -s;

     outmat[1][0]= 0.0;
     outmat[1][1]= 1.0;
     outmat[1][2]= 0.0;

     outmat[2][0]= s;
     outmat[2][1]= 0.0;
     outmat[2][2]= c;
   }

void    rot3mat
        (
          double xval,
          std::vector< std::vector<double> > &outmat
//          double outmat[3][3]
        )
   {
     outmat.resize(3);  // rows
     for (std::vector< std::vector<double> >::iterator it=outmat.begin(); it != outmat.end();++it)
          it->resize(3);
     double c, s;
     c= cos( xval );
     s= sin( xval );

     outmat[0][0]= c;
     outmat[0][1]= s;
     outmat[0][2]= 0.0;

     outmat[1][0]= -s;
     outmat[1][1]= c;
     outmat[1][2]= 0.0;

     outmat[2][0]= 0.0;
     outmat[2][1]= 0.0;
     outmat[2][2]= 1.0;
   }


/* -----------------------------------------------------------------------------
*
*                           procedure addvec
*
*  this procedure adds two vectors possibly multiplied by a constant.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    a1          - constant multiplier
*    a2          - constant multiplier
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    outvec      - vector result of a + b
*
*  locals        :
*    row         - index
*
*  coupling      :
*     none
* --------------------------------------------------------------------------- */

void    addvec
        (
          double a1, double vec1[3],
          double a2, double vec2[3],
          double vec3[3]
        )
   {
     int row;

     for (row = 0; row <= 2; row++)
       {
         vec3[row]= 0.0;
         vec3[row] = a1* vec1[row] + a2* vec2[row];
       }
   }

/* -----------------------------------------------------------------------------
*
*                           procedure addvec3
*
*  this procedure adds three vectors possibly multiplied by a constant.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    a1          - constant multiplier
*    a2          - constant multiplier
*    a3          - constant multiplier
*    vec1        - vector number 1
*    vec2        - vector number 2
*    vec3        - vector number 3
*
*  outputs       :
*    outvec      - vector result of a + b + c
*
*  locals        :
*    row         - index
*
*  coupling      :
*     none
* --------------------------------------------------------------------------- */

void    addvec3
        (
          double a1, double vec1[3],
          double a2, double vec2[3],
          double a3, double vec3[3],
          double vec4[3]
        )
   {
     int row;

     for (row = 0; row <= 2; row++)
       {
         vec4[row]= 0.0;
         vec4[row] = a1* vec1[row] + a2* vec2[row]+ a3* vec3[row];
       }
   }


/* -----------------------------------------------------------------------------
*
*                           procedure angle
*
*  this procedure calculates the angle between two vectors.  the output is
*    set to 999999.1 to indicate an undefined value.  be sure to check for
*    this at the output phase.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    theta       - angle between the two vectors  -Pi to Pi
*
*  locals        :
*    temp        - temporary real variable
*    magv1       - magnitude of vec1
*    magv2       - magnitude of vec2
*    small       - value defining a small value
*    undefined   - large number to use in place of a not defined number
*
*  coupling      :
*    dot           dot product of two vectors
*    acos          arc cosine function
*    mag           magnitude of a vector
* --------------------------------------------------------------------------- */

double  angle
        (
          double vec1[3],
          double vec2[3]
        )
   {
     double small, magv1, magv2, temp;
     small     = 0.00000001;

     magv1 = mag(vec1);
     magv2 = mag(vec2);

     if (magv1*magv2 > small*small)
       {
         temp= dot(vec1,vec2) / (magv1*magv2);
         if (fabs( temp ) > 1.0)
             temp= sgn(temp) * 1.0;
         return acos( temp );
       }
       else
         return undefined;
   }  // end angle

// this writes a vector out to the screen
void    writevec
        (
          char title[10],
          double r[3], double v[3], double a[3]
        )
   {
     printf("%10s  %15.8f%15.8f%15.8f",title, r[0],r[1],r[2] );
     printf(" v %15.9f%15.9f%15.9f"   ,v[0],v[1],v[2] );
     printf(" a %14.9f%14.9f%14.9f\n" ,a[0],a[1],a[2] );
   }

/* -----------------------------------------------------------------------------
*
*                           procedure matvecmult
*
*  this procedure multiplies a 3x3 matrix and a 3x1 vector together.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    mat         - 3 x 3 matrix
*    vec         - vector
*
*  outputs       :
*    vecout      - vector result of mat * vec
*
*  locals        :
*    row         - row index
*    col         - column index
*    ktr         - index
*
*  coupling      :
* --------------------------------------------------------------------------- */

void    matvecmult
        (
          std::vector< std::vector<double> > mat,
//          double mat[3][3],
          double vec[3],
          double vecout[3]
        )
   {
     int row, ktr;

     for (row = 0; row <= 2; row++)
       {
         vecout[row]= 0.0;
         for (ktr = 0; ktr <= 2; ktr++)
             vecout[row]= vecout[row] + mat[row][ktr] * vec[ktr];
       }
   }

/* -----------------------------------------------------------------------------
*
*                           procedure matmult
*
*  this procedure multiplies two matricies up to 10x10 together.
*
*  author        : david vallado                  719-573-2600    7 dec 2007
*
*  inputs          description                    range / units
*    mat1        - matrix number 1
*    mat2        - matrix number 2
*    mat1r       - matrix number 1 rows
*    mat1c       - matrix number 1 columns
*    mat2c       - matrix number 2 columns
*
*  outputs       :
*    mat3        - matrix result of mat1 * mat2 of size mat1r x mat2c
*
*  locals        :
*    row         - row index
*    col         - column index
*    ktr         - index
*
*  coupling      :
* --------------------------------------------------------------------------- */

void    matmult
        (
          std::vector< std::vector<double> > mat1,
          std::vector< std::vector<double> > mat2,
          std::vector< std::vector<double> > &mat3,
//          double mat1[3][3],
//          double mat2[3][3],
//          double mat3[3][3],
          int mat1r, int mat1c, int mat2c
        )
   {
     int row,col,ktr;
     // specify the actual sizes
     mat3.resize(mat1r);  // rows
     for (std::vector< std::vector<double> >::iterator it=mat3.begin(); it != mat3.end();++it)
          it->resize(mat2c);


     for (row = 0; row < mat1r; row++)
       {
         for (col = 0; col < mat2c; col++)
           {
             mat3[row][col]= 0.0;
             for (ktr = 0; ktr < mat1c; ktr ++)
                 mat3[row][col]= mat3[row][col] + mat1[row][ktr] * mat2[ktr][col];
           }
       }
   }

/* -----------------------------------------------------------------------------
*
*                           procedure mattrans
*
*  this procedure finds the transpose of a matrix.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    mat1        - matrix number 1
*    mat1r       - matrix number 1 rows
*    mat1c       - matrix number 1 columns
*
*  outputs       :
*    mat2        - matrix result of transpose mat2
*
*  locals        :
*    row         - row index
*    col         - column index
*
*  coupling      :
* --------------------------------------------------------------------------- */

void    mattrans
        (
          std::vector< std::vector<double> > mat1,
          std::vector< std::vector<double> > &mat2,
//          double mat1[3][3],
//          double mat2[3][3],
          int mat1r, int mat1c
        )
   {
     int row,col;

     mat2.resize(mat1c);  // rows
     for (std::vector< std::vector<double> >::iterator it=mat2.begin(); it != mat2.end();++it)
          it->resize(mat1r);

     for (row = 0; row < mat1r; row++)
       {
         for (col = 0; col < mat1c; col++)
             mat2[col][row] = mat1[row][col];
       }
   }

/* ------------------------------------------------------------------------------
*
*                           procedure ludecomp
*
*  this procedure decomposes a matrix into an lu form.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    order       - order of matrix
*
*  outputs       :
*    lu          - lu decomposition matrix
*    index       - index vector for pivoting
*
*  locals        :
*    i           - index
*    j           - index
*    k           - index
*    imax        - pivot row pointer
*    scale       - scale factor vector
*    sum         - temporary variables
*    amax        - temporary variables
*    dum         - temporary variables
*
*  coupling      :
*    none
*
*  references    :
*    numerical recipes - flannery
 ----------------------------------------------------------------------------- */

void ludecomp
     (
       std::vector< std::vector<double> > &lu,
       std::vector< int > &indexx,
//          double lu[3][3],
//          double indexx[3],
       int order
     )
   {
     const double small = 0.000001;
     int i, j, k, imax;
//     std::vector< std::vector<double> > scale(order+1,2);
	 std::vector< std::vector<double> > scale = std::vector< std::vector<double> >(order + 1, std::vector<double>(2, 0.0));

     double sum, amax, dum;

     imax = 0;
     for (i = 1; i <= order; i++)
       {
         amax = 0.0;
         for (j = 1; j <= order; j++)
             if (fabs( lu[i][j] ) > amax)
                 amax = fabs(lu[i][j]);
         if (fabs(amax) < small)
           {
             printf(" singular matrix amax ");
           }
         scale[i][1] =  1.0 / amax;
       }
     for (j = 1; j <= order; j++)
       {
         for (i = 1; i <= j - 1; i++)
           {
             sum = lu[i][j];
             for (k = 1; k <= i - 1; k++)
                 sum = sum - lu[i][k] * lu[k][j];
             lu[i][j] = sum;
           }
         amax = 0.0;
         for (i = j; i <= order; i++)
           {
             sum = lu[i][j];
             for (k = 1; k <= j - 1; k++)
                 sum = sum - lu[i][k] * lu[k][j];
             lu[i][j] = sum;
             dum = scale[i][1] * fabs(sum);
             if (dum >= amax)
               {
                 imax = i;
                 amax = dum;
               }
           }
         if (j != imax)
           {
             for (k = 1; k <= order; k++)
               {
                 dum = lu[imax][k];
                 lu[imax][k] = lu[j][k];
                 lu[j][k] =  dum;
               }
             scale[imax][1] = scale[j][1];
           }
         indexx[j] = imax;
         if (fabs( lu[j][j]) < small)
           {
             printf(" matrix is singular lu ");
           }
         if (j != order)
           {
             dum = 1.0 / lu[j][j];
             for (i = j + 1; i <= order; i++)
                 lu[i][j] = dum * lu[i][j];
           }
       }
   } // procedure ludecmp

/* ------------------------------------------------------------------------------
*
*                           procedure lubksub
*
*  this procedure finds the inverse of a matrix using lu decomposition.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    order       - order of matrix
*    lu          - lu decomposition matrix
*    index       - index vector for pivoting
*
*  outputs       :
*    b           - solution vector
*
*  locals        :
*    i           - index
*    j           - index
*    i0          - pointer to non-zero element
*    iptr        - pivot rwo pointer
*    sum         - temporary variables
*
*  coupling      :
*    none
*
*  references    :
*    numerical recipes - flannery
 ----------------------------------------------------------------------------- */

void lubksub
     (
       std::vector< std::vector<double> > lu,
       std::vector< int > indexx,
//          double lu[3][3],
//          double indexx[3],
       int order,
       std::vector< std::vector<double> > &b
//          double b[3][3]
     )
   {
     int i, j, iptr, i0;
     double sum;

     i0 = 0;
     for (i = 1; i <= order; i++)
       {
         iptr = indexx[i];
         sum  = b[iptr][1];
         b[iptr][1] = b[i][1];
         if (i0 != 0)
             for (j = i0; j <= i - 1; j++)
                 sum = sum - lu[i][j] * b[j][1];
           else
             if (sum != 0.0)
                i0 = i;
         b[i][1] = sum;
       }

     b[order][1] = b[order][1] / lu[order][order];

     for (i = order - 1; i >= 1; i--)
       {
         sum = b[i][1];
         for (j = i + 1; j <= order; j++)
            sum = sum - lu[i][j] * b[j][1];
         b[i][1] = sum / lu[i][i];
       }
  }  //  procedure lubksub

/* ------------------------------------------------------------------------------
*
*                           procedure matinverse
*
*  this procedure finds the inverse of a matrix using lu decomposition.
*  internally, the function uses arrays starting at 1 to be consistent with
*  other procedures.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    mat         - matrix to invert, 0 array starts
*    order       - order of matrix
*
*  outputs       :
*    matinv      - inverted matrix, 0 array starts
*
*  locals        :
*    i           - index
*    j           - index
*    index       - index vector for pivoting
*    lu          - lu decomposition matrix
*    b           - operational vector to form matinv
*
*  coupling      :
*    ludecomp    -
*    lubksub     -
*
*  references    :
*    numerical recipes - flannery
 ----------------------------------------------------------------------------- */

void matinverse
     ( std::vector< std::vector<double> > mat,
//          double mat[3][3],
       int  order,
       std::vector< std::vector<double> > &matinv
//          double matinv[3][3]
     )
   {
      int i, j;
      std::vector< int > indexx(order+1);
      //std::vector< std::vector<double> > lu(order+1,order+1);
	  //std::vector< std::vector<double> >  b(order+1,2);
	  std::vector< std::vector<double> > lu = std::vector< std::vector<double> >(order + 1, std::vector<double>(order + 1, 0.0));
	  std::vector< std::vector<double> >  b = std::vector< std::vector<double> >(order + 1, std::vector<double>(2, 0.0));

      matinv.resize(order);  // rows
      for (std::vector< std::vector<double> >::iterator it=matinv.begin(); it != matinv.end();++it)
          it->resize(order);

      for (i = 1; i <= order; i++)
        {
          indexx[i] = i;
          for (j = 1; j <= order; j++)
              lu[i][j] = mat[i-1][j-1];
        }
      ludecomp(lu, indexx, order);
       for (j = 1; j <= order; j++)
        {
          for (i = 1; i <= order; i++)
            {
             if (i == j)
                 b[i][1] = 1.0;
               else
                 b[i][1] = 0.0;
            }

          lubksub(lu, indexx, order, b);
          for (i = 1; i <= order; i++)
             matinv[i-1][j-1] = b[i][1];
        }
   } // procedure matinverse


/* ---------------------------------------------------------------------------- 
*
*                           procedure determinant
*
*  This function calculates the determinant value using L - U decompisition.
*    The formula must have a non-zero number in the 1, 1 position. if the
*    function senses a non-zero number in row 1, it exchanges row1 for a row
*    with a non-zero number.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    mat1        - matrix to find determinant of
*    order       - order of matrix
*
*  outputs       :
*    determinant - determinant of mat1
*
*  locals        :
*    i, j, k, n  - index
*    index       - index vector for pivoting
*    lu          - lu decomposition matrix
*    b           - operational vector to form matinv
*
*  coupling      :
*    ludecomp    -
*    lubksub     -
*
*  references    :
*    Marion        pg. 168 - 172, 126 - 127
 ---------------------------------------------------------------------------- - */

double determinant
    (
    std::vector< std::vector<double> > mat1,
    int  order
    )
{
	double small = 0.00000001;
	int i, j, k;
	double temp, d, sum;
	std::vector< std::vector<double> > l = std::vector< std::vector<double> >(order, std::vector<double>(order, 0.0));
	std::vector< std::vector<double> > u = std::vector< std::vector<double> >(order, std::vector<double>(order, 0.0));

	sum = 0.0;
	// ----------- Switch a non zero row to the first row---------- 
	if (abs(mat1[1][1]) < small)
	{
		j = 1;
		while (j <= order)
		{
			if (abs(mat1[j][1]) > small)
			{
				for (k = 1; k <= order; k++)
				{
					temp = mat1[1][k];
					mat1[1][k] = mat1[j][k];
					mat1[j][k] = temp;
				}
				j = order + 1;
			}
			j = j + 1;
		} // if abs
	}

	for (i = 1; i <= order; i++)
		l[i][1] = mat1[i][1];
	for (j = 1; j <= order; j++)
		u[1][j] = mat1[1][j] / l[1][1];
	for (j = 2; j <= order; j++)
	{
		for (i = j; i <= order; i++)
		{
			sum = 0.0;
			for (k = 1; k <= j - 1; k++)
				sum = sum + l[i][k] * u[k][j];
			l[i][j] = mat1[i][j] - sum;
		} // for i 
		u[j][j] = 1.0;
		if (j != order)
		{
			for (i = j + 1; i <= order; i++)
			{
				sum = 0.0;
				for (k = 1; k <= j - 1; k++)
					sum = sum + l[j][k] * u[k][i];
				u[j][i] = (mat1[j][i] - sum) / l[j][j];
			} // for i 
		} // if j 
	} //  for j 
	d = 1.0;
	for (i = 1; i <= order; i++)
		d = d * l[i][i];
	return d;
} // determinant


void writemat
     (
       char matname[30],
       std::vector< std::vector<double> > mat,
//          double mat[3][3],
       int row, int col
     )
   {
     int r, c;

     printf("matrix %15s \n",matname );
     for (r = 0; r < row; r++)
       {
         for (c = 0; c < col; c++)
             printf("%16.11f ",mat[r][c] );
         printf(" \n");
       }
   }

void writeexpmat
     (
       char matname[30],
       std::vector< std::vector<double> > mat,
//          double mat[3][3],
       int row, int col
     )
   {
     int r, c;

     printf("matrix %15s \n",matname );
     for (r = 0; r < row; r++)
       {
         for (c = 0; c < col; c++)
             printf("%14g ",mat[r][c] );
         printf(" \n");
       }
   }

/* -----------------------------------------------------------------------------
*
*                           function cubicspl
*
*  this function performs cubic splining of an input zero crossing
*  function in order to find function values.
*
*  author        : david vallado                  719-573-2600     7 aug 2005
*
*  revisions
*                -
*  inputs          description                    range / units
*    p0,p1,p2,p3 - function values used for splining
*    t0,t1,t2,t3 - time values used for splining
*
*  outputs       :
*    acu0..acu3  - splined polynomial coefficients. acu3 t^3, etc
*
*  locals        : none
*
*  coupling      :
*    none
*
*  references    :
*    vallado       2013, 1034
* --------------------------------------------------------------------------- */

void cubicspl
     (
       double p1, double p2, double p3, double p4,
       double& acu0, double& acu1, double& acu2, double& acu3
     )
     {
        acu0 = p2;
        acu1 = -p1/3.0 - 0.5*p2 + p3 -p4/6.0;
        acu2 = 0.5*p1 - p2 + 0.5*p3;
        acu3 = -p1/6.0 + 0.5*p2 - 0.5*p3 + p4/6.0;
     }

/* -----------------------------------------------------------------------------
*
*                           function quadric
*
*  this function solves for the two roots of a quadric equation.  there are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  the general form is y = ax2 + bx + c.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  revisions
*    vallado     - convert to matlab              719-573-2600    3 dec 2002
*
*  inputs          description                    range / units
*    a           - coefficient of x squared term
*    b           - coefficient of x term
*    c           - constant
*    opt         - option for output              I all roots including imaginary
*                                                 R only real roots
*                                                 U only unique real roots (no repeated)
*
*  outputs       :
*    r1r         - real portion of root 1
*    r1i         - imaginary portion of root 1
*    r2r         - real portion of root 2
*    r2i         - imaginary portion of root 2
*
*  locals        :
*    discrim     - discriminate b2 - 4ac
*
*  coupling      :
*    none.
*
*  references    :
*    vallado       2013, 1027
* ----------------------------------------------------------------------------*/

void quadric
     (
       double a, double b, double c, char opt,
       double& r1r, double& r1i, double& r2r, double& r2i
     )
   {
        const double small    = 0.0000001;
        double discrim;
        // --------------------  implementation   ----------------------
        r1r = 0.0;
        r1i = 0.0;
        r2r = 0.0;
        r2i = 0.0;

        discrim = b*b - 4.0 *a*c;

        // ---------------------  real roots  --------------------------
        if ( fabs(discrim) < small  )
          {
            r1r = -b / ( 2.0 *a );
            r2r = r1r;
            if (opt=='U')
                r2r = 99999.9;
          }
          else
           {
            if (fabs(a) < small)
                 r1r = -c/b;
              else
               {
                if ( discrim > 0.0  )
                  {
                    r1r = ( -b + sqrt(discrim) ) / ( 2.0 *a );
                    r2r = ( -b - sqrt(discrim) ) / ( 2.0 *a );
                  }
                  else
                  {
                    // ------------------ complex roots --------------------
                    if (opt=='I')
                      {
                        r1r = -b / ( 2.0 *a );
                        r2r = r1r;
                        r1i =  sqrt(-discrim) / ( 2.0 *a );
                        r2i = -sqrt(-discrim) / ( 2.0 *a );
                      }
                      else
                      {
                        r1r = 99999.9;
                        r2r = 99999.9;
                      }
                  }
                }
             }
   }  // quadric


/* -----------------------------------------------------------------------------
*
*                           function cubic
*
*  this function solves for the three roots of a cubic equation.  there are
*    no restrictions on the coefficients, and imaginary results are passed
*    out as separate values.  the general form is y = a3x3 + b2x2 + c1x + d0.  note
*    that r1i will always be zero since there is always at least one real root.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  revisions
*    vallado     - convert to matlab              719-573-2600   18 dec 2002
*
*  inputs          description                    range / units
*    a3          - coefficient of x cubed term
*    b2          - coefficient of x squared term
*    c1          - coefficient of x term
*    d0          - constant
*    opt         - option for output              I all roots including imaginary
*                                                 R only real roots
*                                                 U only unique real roots (no repeated)
*
*  outputs       :
*    r1r         - real portion of root 1
*    r1i         - imaginary portion of root 1
*    r2r         - real portion of root 2
*    r2i         - imaginary portion of root 2
*    r3r         - real portion of root 3
*    r3i         - imaginary portion of root 3
*
*  locals        :
*    temp1       - temporary value
*    temp2       - temporary value
*    p           - coefficient of x squared term where x cubed term is 1.0
*    q           - coefficient of x term where x cubed term is 1.0
*    r           - coefficient of constant term where x cubed term is 1.0
*    delta       - discriminator for use with cardans formula
*    e0          - angle holder for trigonometric solution
*    phi         - angle used in trigonometric solution
*    cosphi      - cosine of phi
*    sinphi      - sine of phi
*
*  coupling      :
*    quadric     - quadratic roots
*
*  references    :
*    vallado       2013, 1027
* --------------------------------------------------------------------------- */

void cubic
         (
           double a3, double b2, double c1, double d0, char opt,
           double& r1r, double& r1i, double& r2r, double& r2i, double& r3r, double& r3i
         )
   {
        const double rad      = 57.29577951308230;
        const double onethird = 1.0/3.0;
        const double small    = 0.00000001;
        double temp1, temp2, p, q, r, delta, e0, cosphi, sinphi, phi;
        // ------------------------  implementation   --------------------------
        r1r  = 0.0;
        r1i  = 0.0;
        r2r  = 0.0;
        r2i  = 0.0;
        r3r  = 0.0;
        r3i  = 0.0;

        if (fabs(a3) > small)
          {
           // ------------- force coefficients into std form -------------------
            p= b2/a3;
            q= c1/a3;
            r= d0/a3;

            a3= onethird*( 3.0 *q - p*p );
            b2= (1.0 /27.0 )*( 2.0 *p*p*p - 9.0 *p*q + 27.0 *r );

            delta= (a3*a3*a3/27.0 ) + (b2*b2*0.25 );

            // -------------------- use cardans formula ------------------------
            if ( delta > small )
              {
                temp1= (-b2*0.5 )+sqrt(delta);
                temp2= (-b2*0.5 )-sqrt(delta);
                temp1= sgn(temp1)*pow( fabs(temp1),onethird );
                temp2= sgn(temp2)*pow( fabs(temp2),onethird );
                r1r= temp1 + temp2 - p*onethird;

                if (opt=='I')
                  {
                    r2r= -0.5 *(temp1 + temp2) - p*onethird;
                    r2i= -0.5 *sqrt( 3.0  )*(temp1 - temp2);
                    r3r= -0.5 *(temp1 + temp2) - p*onethird;
                    r3i= -r2i;
                  }
                  else
                  {
                    r2r= 99999.9;
                    r3r= 99999.9;
                  }
               }
              else
               {
            // -------------------- evaluate zero point ------------------------
            if ( fabs( delta ) < small  )
              {
                r1r= -2.0*sgn(b2)*pow(fabs(b2*0.5),onethird) - p*onethird;
                r2r=      sgn(b2)*pow(fabs(b2*0.5),onethird) - p*onethird;
                if (opt=='U')
                    r3r= 99999.9;
                  else
                    r3r= r2r;
              }
              else
              {
                // --------------- use trigonometric identities ----------------
                e0     = 2.0 *sqrt(-a3*onethird);
                cosphi = (-b2/(2.0 *sqrt(-a3*a3*a3/27.0 )) );
                sinphi = sqrt( 1.0 -cosphi*cosphi );
                phi    = atan2( sinphi,cosphi );
                if (phi < 0.0)
                    phi = phi + 2.0*pi;
                r1r= e0*cos( phi*onethird ) - p*onethird;
                r2r= e0*cos( phi*onethird + 120.0 /rad ) - p*onethird;
                r3r= e0*cos( phi*onethird + 240.0 /rad ) - p*onethird;
              } // if fabs(delta)...
             }  // if delta > small
          }  // if fabs > small
        else
        {
          quadric( b2, c1, d0, opt, r1r, r1i, r2r, r2i );
          r3r  = 99999.9;
          r3i  = 99999.9;
        }
   }  // cubic

/* -----------------------------------------------------------------------------
*
*                           function cubicinterp
*
*  this function performs a cubic spline. four points are needed.
*
*  author        : david vallado                  719-573-2600   1 dec  2005
*
*  revisions
*
*  inputs          description                    range / units
*    valuein     - kp
*
*  outputs       :
*    out         - ap
*
*  locals        :
*                -
*
*  coupling      :
*    cubicspl
*
*  references    :
*    vallado       2013, 1034
* --------------------------------------------------------------------------- */

double  cubicinterp
        (
          double p1a, double p1b, double p1c, double p1d, double p2a, double p2b,
          double p2c, double p2d, double valuein
        )
   {
        double kc0, kc1, kc2, kc3, ac0, ac1, ac2, ac3,
               r1r, r1i, r2r, r2i, r3r, r3i, value;

           // -------- assign function points ---------
           cubicspl(p1a, p1b, p1c, p1d, ac0, ac1, ac2, ac3);
           cubicspl(p2a, p2b, p2c, p2d, kc0, kc1, kc2, kc3);

           // recover the original function values
           // use the normalized time first, but at an arbitrary interval
           cubic(kc3, kc2, kc1, kc0-valuein, 'R', r1r, r1i, r2r, r2i, r3r, r3i);

           if ((r1r >= -0.000001) && (r1r <= 1.001))
             {
               value = r1r;
             }
             else
             {
               if ((r2r >= -0.000001) && (r2r <= 1.001))
                 {
                   value = r2r;
                 }
                 else
                 {
                   if ((r3r >= -0.000001) && (r3r <= 1.001))
                     {
                       value = r3r;
                     }
                     else
                     {
                       value = 0.0;
                       printf("error in cubicinterp root %17.14f %11.7f %11.7f %11.7f \n",
                               valuein,r1r,r2r,r3r);
                     }
                 }
             }
           return (ac3*pow(value,3)  + ac2*value*value  + ac1*value + ac0);

   } // cubicinterp


/* -----------------------------------------------------------------------------
*
*                           function factorial
*
*  this function performs a factorial. 
*
*  author        : david vallado                  719-573-2600   11 feb 2016
*
*  revisions
*
*  inputs          description                    range / units
*    n           - order in
*
*  outputs       :
*    factorial   - result
*
*  locals        :
*                -
*
*  coupling      :
*    none
* --------------------------------------------------------------------------- */

int factorial(int n)
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
}



} // namespace

