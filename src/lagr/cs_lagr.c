/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Utilitarian functions for the diphasic lagrangian module
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* variable for _orient3D_set_maxvalue */

static double MAX24;               /* bounding box size */
static double SNAP_TO_GRID;        /* for rounding of coordinates */
static double Orient3D_split_2_25; /* split to bit MAX0^2 * 2^25 */
static double C1_ORIENTATION_3D;   /* static orientation filter */
static double C2_ORIENTATION_3D;   /* semi-static orientation filter */
static double C3_ORIENTATION_3D;   /* half degree 3 unit for semi-static filter */
static double C4_PERTURB_OR_3D;    /* static filter for perturbation */

/* variable for _orientation3D */

static const cs_int_t POSITIVE   =  1 ;
static const cs_int_t NEGATIVE   = -1 ;
static const cs_int_t COPLANAR   =  0 ;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if IEEE 754 standard is respected for floating storage for this
 * architecture. Test on the accuracy of the computation.
 *
 * Returns:
 *  integer (EXIT_FAILURE / EXIT_SUCCES)
 *----------------------------------------------------------------------------*/

static cs_int_t
_check_ieee_standard(void)
{
  double  dd;

  double  d = 9007199254740992.0;  /* 2^53 */

  dd = d + 1;
  if (dd != d) return EXIT_FAILURE;  /* compute with too much bits */

  dd = d - 1;
  if (dd == d) return EXIT_FAILURE;  /* compute with not enough bits */

  return EXIT_SUCCESS;
}

/*----------------------------------------------------------------------------
 * The argument is an upperbound on the coordinates.
 * Change the data conditioning.
 *
 * d       --> upper bound on data values.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

static inline void
_orient3D_set_maxvalue(double  d)
{
  float  D;

  /* Coordinates precision */
  double MAX0 = 1.0;                                 /* 2^0 */

  MAX24 = 16777216.0;                                /* 2^24 */
  SNAP_TO_GRID = 6755399441055744.0;                 /* 3 * 2^51 */
  Orient3D_split_2_25 = 453347182355485940514816.0;  /* 3 * 2^77 */
  C1_ORIENTATION_3D = 75497472.0;                    /* 9 * 2^23 */
  C2_ORIENTATION_3D = 3.0/9007199254740992.0;        /* 3 * 2^(-53) */
  C3_ORIENTATION_3D = 0.5;                           /* 1/2 */
  C4_PERTURB_OR_3D  = 5066549580791808.0;            /* 9 * 2^49 */

  D = d * 4503599627370497.0;         /* 2^52 + 1 */

  /* On x86 or compatible, force extended-precision (80-bit) calculation
     using x87 FPU instead of SSE instructions (64-bit). */

#if defined(__GNUC__) && (   defined(__i386) \
                          || defined(__x86_64__) || defined(__amd64__))
  {
    double s = D;
    __asm__ __volatile__("fldl %1\n\t"  \
                         "faddl %2\n\t" \
                         "fsubl %2\n\t" \
                         "fstpl %0"
                         : "=m" (s)
                         : "m" (d), "m" (s));
    D = s;
  }

/* On other architectures, a compiler option such ICC's
   "-fp-model extended" is required. If the architecture does not
   support Intel-type 80-bit precision calculations, the code
   must be adapted. */

#else
  D = (d + D) -D;                     /* get the power of two closer to d */
#endif

  if (D < ((float)d))
    D *= 2;                           /* grid max size 2^b */

  MAX0 = D / 16777216.0;              /* grid step 2^k = 2^b / 2^24 */
  MAX24 = D;
  SNAP_TO_GRID *= MAX0;
  Orient3D_split_2_25 *= MAX0*MAX0;
  C1_ORIENTATION_3D *= MAX0*MAX0*MAX0;
  C3_ORIENTATION_3D *= MAX0*MAX0*MAX0;
  C4_PERTURB_OR_3D *= MAX0*MAX0*MAX0*MAX0;
}

/*----------------------------------------------------------------------------
 * Round a coordinate on the grid
 *
 * value         --> valoe to round
 *
 * Returns:
 *----------------------------------------------------------------------------*/

static inline void
_orient3D_normalize(float  *value)
{
  /* round a value with fixed precision */

  if ( (*value > MAX24) || (*value < -MAX24)) {

    bft_error(__FILE__, __LINE__, 0,
              _("overflow |%g| > %g\nVerify the bounding box"
                " for your data."), *value, MAX24);

    *value=0.0;

  }
  else {

    /* On x86 or compatible, force extended-precision (80-bit) calculation
       using x87 FPU instead of SSE instructions (64-bit). */

#if defined(__GNUC__) && (   defined(__i386) \
                          || defined(__x86_64__) || defined(__amd64__))
    {
      double d = *value, s = SNAP_TO_GRID;
      __asm__ __volatile__("fldl %1\n\t"  \
                           "faddl %2\n\t" \
                           "fsubl %2\n\t" \
                           "fstpl %0"
                           : "=m" (d)
                           : "m" (d), "m" (s));
      *value = d;
    }

    /* On other architectures, a compiler option such ICC's
       "-fp-model extended" is required. If the architecture does not
       support Intel-type 80-bit precision calculations, the code
       must be adapted. */
#else
    *value = ( *value + SNAP_TO_GRID ) - SNAP_TO_GRID;
#endif
  }

}

/*----------------------------------------------------------------------------
 * Split a in two numbers. a=a1+a0, a1 get the most significant digit
 * format specify the range of the input, and how to split if format is
 * _split_i_j it means that a is of degree i (in terms of original coordinates)
 *  the splitting bit is MAX0^i * 2^j
 *
 * inline is not an ansi feature and not recognized everywhere
 * BUT it looks as if -O optimization level lead to inlining anyway
 * on Linux machines
 *
 * a              --> value to treat
 * a0             <-- a = a0 + a1
 * a1             <-- a = a0 + a1
 * format         -->
 *
 * Returns:
 *----------------------------------------------------------------------------*/

static inline void
_orient3D_split(double    a,
                double   *a1,
                double   *a0,
                double    format)
{
  volatile double _a1 = (a + format);
  *a1 = _a1 -format;
  *a0 = a - *a1;
}

/*----------------------------------------------------------------------------
 * Compute the orientation of the four points assumed to be rounded on
 * the grid.
 * The last argument specify if the perturbation technique is used in case of
 * coplanarity
 *
 * ax                   --> X coordinate of the first vertex
 * ay                   --> Y coordinate of the first vertex
 * az                   --> Z coordinate of the first vertex
 * bx                   --> X coordinate of the second vertex
 * by                   --> Y coordinate of the second vertex
 * bz                   --> Z coordinate of the second vertex
 * cx                   --> X coordinate of the third vertex
 * cy                   --> Y coordinate of the third vertex
 * cz                   --> Z coordinate of the third vertex
 * dx                   --> X coordinate of the fourth vertex
 * dy                   --> Y coordinate of the fourth vertex
 * dz                   --> Z coordinate of the fourth vertex
 * pturb                <->
 *
 * Returns:
 *  orientation of the four vertices (a,b,c,d).
 *  0 = COPLANAR, 1 = POSITIVE, -1 = NEGATIVE.
 *----------------------------------------------------------------------------*/

static cs_int_t
_orientation3D(float      ax,
               float      ay,
               float      az,
               float      bx,
               float      by,
               float      bz,
               float      cx,
               float      cy,
               float      cz,
               float      dx,
               float      dy,
               float      dz,
               cs_int_t   PERTURBED)
{
  double  error;
  double  M10, M11, M20, M21, M30, M31;
  double  RA, RB, RC, RD;
  double  RB0, RB1, RC0, RC1, RD0, RD1;

  /* points are assumed to be distincts */

  double  A = (double)bx - (double)ax;    /*    x y z   */
  double  B = (double)by - (double)ay;
  double  C = (double)bz - (double)az;    /*  | A B C | */
  double  D = (double)cx - (double)ax;    /*  | D E F | */
  double  E = (double)cy - (double)ay;    /*  | G H I | */
  double  F = (double)cz - (double)az;
  double  G = (double)dx - (double)ax;
  double  H = (double)dy - (double)ay;
  double  I = (double)dz - (double)az;

  /* minors computation */

  double  M1 = D*H - E*G;
  double  M2 = A*H - B*G;
  double  M3 = A*E - B*D;

  /* determinant computation */

  double  det = (C*M1 - F*M2) + I*M3;

  if (det >   C1_ORIENTATION_3D ) return POSITIVE;
  if (det <  -C1_ORIENTATION_3D ) return NEGATIVE;

  /* det is small */
  /* Orientation 3D, continuing */

  error = (fabs(C*M1) + fabs(F*M2)) + fabs(I*M3);
  error *= C2_ORIENTATION_3D;

  if (error < C3_ORIENTATION_3D) error = 0.0;
  if (det >   error) return POSITIVE;
  if (det <  -error) return NEGATIVE;

  /* det is very small, exact computation must be done */
  /* Orientation 3D, continuing */

  if (error != 0.0) { /* otherwise, exact value 0 already certified */

    _orient3D_split(M1, &M11, &M10, Orient3D_split_2_25);
    _orient3D_split(M2, &M21, &M20, Orient3D_split_2_25);
    _orient3D_split(M3, &M31, &M30, Orient3D_split_2_25);

    det  =  C*M11 - F*M21 + I*M31;
    det += (C*M10 - F*M20 + I*M30);     /* less significant */

    if (det>0) return POSITIVE;
    if (det<0) return NEGATIVE;

  }

  /* points are coplanar */

  if (! PERTURBED)  return COPLANAR;

  /* start perturbation scheme */
  /* (x,y,x)  |-->(x+e^3(x^2+y^2+z^2), y+e^2(x^2+y^2+z^2), z +e(x^2+y^2+z^2)) */

  RA = ax*ax + ay*ay + az*az;
  RB = bx*bx + by*by + bz*bz - RA;
  RC = cx*cx + cy*cy + cz*cz - RA;
  RD = dx*dx + dy*dy + dz*dz - RA;

  /* epsilon term */

  det = (RB*M1 - RC*M2) + RD*M3; /* reuse minors */     /* A B RB */
  if (det >  C4_PERTURB_OR_3D) return POSITIVE;         /* D E RC */
  if (det < -C4_PERTURB_OR_3D) return NEGATIVE;         /* G H RD */

  /* det is small, doing exact computation */

  /* modif EDF deb */
  _orient3D_split(M1, &M11, &M10, Orient3D_split_2_25);
  _orient3D_split(M2, &M21, &M20, Orient3D_split_2_25);
  _orient3D_split(M3, &M31, &M30, Orient3D_split_2_25);
   /* modif EDF fin */

  _orient3D_split(RB, &RB1, &RB0, Orient3D_split_2_25);
  _orient3D_split(RC, &RC1, &RC0, Orient3D_split_2_25);
  _orient3D_split(RD, &RD1, &RD0, Orient3D_split_2_25);

  det  =  RB1*M11 - RC1*M21 + RD1*M31;
  det += (RB1*M10 + RB0*M11) - (RC1*M20 + RC0*M21) + (RD1*M30 + RD0*M31);
  det += (RB0*M10 - RC0*M20 + RD0*M30);      /* less significant */
  if (det>0) return POSITIVE;
  if (det<0) return NEGATIVE;

  /* points coplanar in special direction go to e^2 term */

  M1 = D*I - F*G;                           /* A RB C */
  M2 = A*I - C*G;                           /* D RC F */
  M3 = A*F - C*D;                           /* G RD I */
  det = (-RB*M1 + RC*M2) - RD*M3;
  if (det >  C4_PERTURB_OR_3D) return POSITIVE;
  if (det < -C4_PERTURB_OR_3D) return NEGATIVE;

  /* det is small, doing exact computation */

  _orient3D_split(M1, &M11, &M10, Orient3D_split_2_25);
  _orient3D_split(M2, &M21, &M20, Orient3D_split_2_25);
  _orient3D_split(M3, &M31, &M30, Orient3D_split_2_25);

  det  = -RB1*M11 + RC1*M21 - RD1*M31;
  det += -(RB1*M10 + RB0*M11) + (RC1*M20 + RC0*M21) - (RD1*M30 + RD0*M31);
  det += (-RB0*M10 + RC0*M20 - RD0*M30);      /* less significant */
  if (det>0) return POSITIVE;
  if (det<0) return NEGATIVE;

  /* points coplanar in special direction go to e^3 term */

  M1 = E*I - F*H;                           /* RB B C */
  M2 = B*I - C*H;                           /* RC E F */
  M3 = B*F - C*E;                           /* RD H I */
  det = (RB*M1 - RC*M2) + RD*M3;
  if (det >  C4_PERTURB_OR_3D) return POSITIVE;
  if (det < -C4_PERTURB_OR_3D) return NEGATIVE;

  /* det is small, doing exact computation */

  _orient3D_split (M1, &M11, &M10, Orient3D_split_2_25);
  _orient3D_split (M2, &M21, &M20, Orient3D_split_2_25);
  _orient3D_split (M3, &M31, &M30, Orient3D_split_2_25);
  det  =  RB1*M11 - RC1*M21 + RD1*M31;
  det += (RB1*M10 + RB0*M11) - (RC1*M20 + RC0*M21) + (RD1*M30 + RD0*M31);
  det += (RB0*M10 - RC0*M20 + RD0*M30);      /* less significant */
  if (det>0) return POSITIVE;
  if (det<0) return NEGATIVE;

  /* points are colinear */

  return COPLANAR;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if IEEE 754 standard is respected for floating storage for this
 * architecture. If the standard is not respected the particle trajectography
 * may be wrong.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (csieee,CSIEEE)(void)
{

  if (_check_ieee_standard()) {

    bft_error(__FILE__, __LINE__, 0,
              _("IEEE 754 arithmetic is not supported by the "
                "current architecture."));

  }
}

/*----------------------------------------------------------------------------
 * Check the relative localization of two vertices. We want to know if these
 * two vertices are identical.
 *
 * pvalmax              --> upperbound on coordinates
 * px                   --> X coordinate of the vertex P
 * py                   --> Y coordinate of the vertex P
 * pz                   --> Z coordinate of the vertex P
 * qx                   --> X coordinate of the vertex Q
 * qy                   --> Y coordinate of the vertex Q
 * qz                   --> Z coordinate of the vertex Q
 * sign                 <-> return tag (1 -> identical else 0)
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (coloca,COLOCA)(cs_real_t  *pvalmax,
                         cs_real_t  *px,
                         cs_real_t  *py,
                         cs_real_t  *pz,
                         cs_real_t  *qx,
                         cs_real_t  *qy,
                         cs_real_t  *qz,
                         cs_int_t   *sign)
{
  float  xx0, yy0, zz0;
  float  xx1, yy1, zz1;

  xx0 = (float) *px;
  yy0 = (float) *py;
  zz0 = (float) *pz;

  xx1 = (float) *qx;
  yy1 = (float) *qy;
  zz1 = (float) *qz;

  /* We look for the nearest and upper value to valmax by power of 2 */

  _orient3D_set_maxvalue(*pvalmax);

  _orient3D_normalize(&xx0);
  _orient3D_normalize(&yy0);
  _orient3D_normalize(&zz0);

  _orient3D_normalize(&xx1);
  _orient3D_normalize(&yy1);
  _orient3D_normalize(&zz1);

  /* We check if the two vertices are identical */

  if ( (xx0==xx1) && (yy0==yy1) && (zz0==zz1) )
    *sign = 1;
  else
    *sign = 0;

}

/*----------------------------------------------------------------------------
 * Look for coordinate system orientation to locate particles in relation to
 * faces.
 *
 * pvalmax              --> upper bound on coordinates
 * px                   --> X coordinate of the first vertex
 * py                   --> Y coordinate of the first vertex
 * pz                   --> Z coordinate of the first vertex
 * qx                   --> X coordinate of the second vertex
 * qy                   --> Y coordinate of the second vertex
 * qz                   --> Z coordinate of the second vertex
 * cdgx                 --> X coordinate of the third vertex
 * cdgy                 --> Y coordinate of the third vertex
 * cdgz                 --> Z coordinate of the third vertex
 * crgx                 --> X coordinate of the fourth vertex
 * crgy                 --> Y coordinate of the fourth vertex
 * crgz                 --> Z coordinate of the fourth vertex
 * sign                 <-> orientation of the four vertices.
 * pturb                <->
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (coturn,COTURN)(cs_real_t   *pvalmax,
                         cs_real_t   *px,
                         cs_real_t   *py,
                         cs_real_t   *pz,
                         cs_real_t   *qx,
                         cs_real_t   *qy,
                         cs_real_t   *qz,
                         cs_real_t   *cdgx,
                         cs_real_t   *cdgy,
                         cs_real_t   *cdgz,
                         cs_real_t   *crdx,
                         cs_real_t   *crdy,
                         cs_real_t   *crdz,
                         cs_int_t    *sign,
                         cs_int_t    *pturb)
{
  float xx0, yy0, zz0;
  float xx1, yy1, zz1;
  float xx2, yy2, zz2;
  float xx3, yy3, zz3;

  /* Several "cast" */

  xx0 = (float) *px;
  yy0 = (float) *py;
  zz0 = (float) *pz;

  xx1 = (float) *qx;
  yy1 = (float) *qy;
  zz1 = (float) *qz;

  xx2 = (float) *cdgx;
  yy2 = (float) *cdgy;
  zz2 = (float) *cdgz;

  xx3 = (float) *crdx;
  yy3 = (float) *crdy;
  zz3 = (float) *crdz;

  /* We look for the nearest and upper value to valmax by power of 2 */

  _orient3D_set_maxvalue(*pvalmax);

  _orient3D_normalize(&xx0);
  _orient3D_normalize(&yy0);
  _orient3D_normalize(&zz0);

  _orient3D_normalize(&xx1);
  _orient3D_normalize(&yy1);
  _orient3D_normalize(&zz1);

  _orient3D_normalize(&xx2);
  _orient3D_normalize(&yy2);
  _orient3D_normalize(&zz2);

  _orient3D_normalize(&xx3);
  _orient3D_normalize(&yy3);
  _orient3D_normalize(&zz3);

  /* *sign == 0 => COPLANAR,
     *sign == 1 => POSITIVE orientation,
     *sign == 1 => NEGATIVE orientation */

  *sign = _orientation3D(xx0, yy0, zz0,
                         xx1, yy1, zz1,
                         xx2, yy2, zz2,
                         xx3, yy3, zz3,
                         *pturb);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
