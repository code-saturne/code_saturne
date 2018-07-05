/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_cdo_connect.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdovb_vecteq.h"
#include "cs_equation_param.h"
#include "cs_evaluate.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_param.h"
#include "cs_param_cdo.h"
#include "cs_scheme_geometry.h"
#include "cs_source_term.h"
#include "cs_time_step.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Static global variables
 *============================================================================*/

static cs_cdo_connect_t  *connect = NULL;
static cs_cdo_quantities_t  *quant = NULL;
static cs_time_step_t  *time_step = NULL;

static FILE  *quadrature = NULL;

static const double non_poly_int[] = { 1.1319870772271508,  /* Hexa */
                                       0.0801351697078868}; /* Tetra */
static double non_poly_int_f[2][6] = {
  {0.6587901114231911, 0.6587901114231911, 1.7907771886504419,
   1.7907771886504419, 0.6587901114231911, 1.7907771886504419},
  {0.2231301601484198, 0.2231301601484198, 0.2231301601484198,
   0.5252709594852753, 0., 0.}
};

/*============================================================================
 * Structures definition
 *============================================================================*/

typedef void
(_evaluate_t)(cs_flag_t        dof_flag,
              const cs_xdef_t *def,
              double           retval[]);

typedef enum {

  CONSTANT,
  LINEAR,
  QUADRATIC,
  CUBIC,
  NON_POLYNOMIAL,
  N_FUNCTION_TYPES

} _function_type_t;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/* Test function based on analytic definition
 * Fill array with 1
 */
static void
_unity(cs_real_t         time,
       cs_lnum_t         n_pts,
       const cs_lnum_t  *pt_ids,
       const cs_real_t  *xyz,
       bool              compact,
       void             *input,
       cs_real_t         retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(xyz);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? pt_ids[i] : i;
    const cs_lnum_t  r = compact ? i : p;
    retval[r] = 1.0;
  }
}

/* Test function based on analytic definition
 * Fill array with [1, 1, 1]
 */
static void
_unity_vect(cs_real_t         time,
            cs_lnum_t         n_pts,
            const cs_lnum_t  *pt_ids,
            const cs_real_t  *xyz,
            bool              compact,
            void             *input,
            cs_real_t         retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(xyz);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? pt_ids[i] : i;
    const cs_lnum_t  r = compact ? i : p;
    for (int j = 0; j < 3; j++)
      retval[3*r+j] = 1.0;
  }
}

/* Test function based on analytic definition
 * Fill array with [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
 */
static void
_unity_tens(cs_real_t         time,
            cs_lnum_t         n_pts,
            const cs_lnum_t  *pt_ids,
            const cs_real_t  *xyz,
            bool              compact,
            void             *input,
            cs_real_t         retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(xyz);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? pt_ids[i] : i;
    const cs_lnum_t  r = compact ? i : p;
    for (int j = 0; j < 9; j++)
      retval[9*r+j] = 1.0;
  }
}

/* Test function based on analytic definition
 * Fill array with x+y+z
 */
static void
_linear_xyz(cs_real_t          time,
            cs_lnum_t          n_pts,
            const cs_lnum_t   *pt_ids,
            const cs_real_t   *xyz,
            bool               compact,
            void              *input,
            cs_real_t          retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? pt_ids[i] : i;
    const cs_lnum_t  r = compact ? i : p;
    retval[r] = xyz[3*p] + xyz[3*p+1] + xyz[3*p+2];
  }
}

/* Test function based on analytic definition
 * Fill array with [x, 2*y, 3*z]
 */
static void
_linear_xyz_vect(cs_real_t          time,
                 cs_lnum_t          n_pts,
                 const cs_lnum_t   *pt_ids,
                 const cs_real_t   *xyz,
                 bool               compact,
                 void              *input,
                 cs_real_t          retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? 3*pt_ids[i] : 3*i;
    const cs_lnum_t  r = compact ? 3*i : p;
    for (int j = 0; j < 3; j++)
      retval[r+j] = (j+1) * xyz[p+j];
  }

}

/* Test function based on analytic definition
 * Fill array with [[x+y+z, x+y+z, x+y+z],
 *                  [x+y+z, x+y+z, x+y+z],
 *                  [x+y+z, x+y+z, x+y+z]]
 */
static void
_linear_xyz_tens(cs_real_t          time,
                 cs_lnum_t          n_pts,
                 const cs_lnum_t   *pt_ids,
                 const cs_real_t   *xyz,
                 bool               compact,
                 void              *input,
                 cs_real_t          retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? 3*pt_ids[i] : 3*i;
    const cs_lnum_t  r = compact ? 9*i : 3*p;
    const cs_real_t  xpypz = xyz[p] + xyz[p+1] + xyz[p+2];
    for (int j = 0; j < 9; j++)
      retval[r+j] = xpypz;
  }
}

/* Test function based on analytic definition
 * Fill array with x*x
 */
static void
_quadratic_x2(cs_real_t          time,
              cs_lnum_t          n_pts,
              const cs_lnum_t   *pt_ids,
              const cs_real_t   *xyz,
              bool               compact,
              void              *input,
              cs_real_t          retval[])
{
  CS_UNUSED(input);
  CS_UNUSED(time);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? pt_ids[i] : i;
    const cs_lnum_t  r = compact ? i : p;
    retval[r] = xyz[3*p]*xyz[3*p];
  }
}

/* Test function based on analytic definition
 * Fill array with [x*x, x*x, x*x]
 */
static void
_quadratic_x2_vect(cs_real_t          time,
                   cs_lnum_t          n_pts,
                   const cs_lnum_t   *pt_ids,
                   const cs_real_t   *xyz,
                   bool               compact,
                   void              *input,
                   cs_real_t          retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? 3*pt_ids[i] : 3*i;
    const cs_lnum_t  r = compact ? 3*i : p;
    const cs_real_t  x2 = xyz[p]*xyz[p];
    for (int j = 0; j < 3; j++)
      retval[r+j] = x2;
  }

}

/* Test function based on analytic definition
 * Fill array with [[x*x, x*x, x*x],
                    [x*x, x*x, x*x],
                    [x*x, x*x, x*x]]
 */
static void
_quadratic_x2_tens(cs_real_t          time,
                   cs_lnum_t          n_pts,
                   const cs_lnum_t   *pt_ids,
                   const cs_real_t   *xyz,
                   bool               compact,
                   void              *input,
                   cs_real_t          retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? 3*pt_ids[i] : 3*i;
    const cs_lnum_t  r = compact ? 9*i : 3*p;
    const cs_real_t  x2 = xyz[p]*xyz[p];
    for (int j = 0; j < 9; j++)
      retval[r+j] = x2;
  }

}

/* Test function based on analytic definition
 * Fill array with exp(x+y+z-1.5)
 * Integral over the unitary hexa    = (e-1)^3 / e^(3/2)   = 1.1319870772271508
 * Integral over the reference tetra = (e-2) / (2 e^(3/2)) = 0.0801351697078868
 *
 * I = int_(x=0 to 1) int_(y=0 to 1-x) int_(z=0 to 1-x-y) f
 */
static void
_nonpoly(cs_real_t         time,
         cs_lnum_t         n_pts,
         const cs_lnum_t  *pt_ids,
         const cs_real_t  *xyz,
         bool              compact,
         void             *input,
         cs_real_t         retval[])
{
  CS_UNUSED(input);
  CS_UNUSED(time);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? pt_ids[i] : i;
    const cs_lnum_t  r = compact ? i : p;
    retval[r] = exp(xyz[3*p]+xyz[3*p+1]+xyz[3*p+2]-1.5);
  }
}

/* Test function based on analytic definition
 * Fill array with [exp(x+y+z-1.5), exp(x+y+z-1.5), exp(x+y+z-1.5)]
 */
static void
_nonpoly_vect(cs_real_t          time,
              cs_lnum_t          n_pts,
              const cs_lnum_t   *pt_ids,
              const cs_real_t   *xyz,
              bool               compact,
              void              *input,
              cs_real_t          retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? 3*pt_ids[i] : 3*i;
    const cs_lnum_t  r = compact ? 3*i : p;
    const cs_real_t  eval = exp(xyz[p]+xyz[p+1]+xyz[p+2]-1.5);
    for (int j = 0; j < 3; j++)
      retval[r+j] = eval;
  }

}

/* Test function based on analytic definition
 * Fill array with [[exp(x+y+z-1.5), exp(x+y+z-1.5), exp(x+y+z-1.5)],
 *                  [exp(x+y+z-1.5), exp(x+y+z-1.5), exp(x+y+z-1.5)],
 *                  [exp(x+y+z-1.5), exp(x+y+z-1.5), exp(x+y+z-1.5)]]
 */
static void
_nonpoly_tens(cs_real_t          time,
              cs_lnum_t          n_pts,
              const cs_lnum_t   *pt_ids,
              const cs_real_t   *xyz,
              bool               compact,
              void              *input,
              cs_real_t          retval[])
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  for (cs_lnum_t i = 0; i < n_pts; i++) {
    const cs_lnum_t  p = (pt_ids != NULL) ? 3*pt_ids[i] : 3*i;
    const cs_lnum_t  r = compact ? 9*i : 3*p;
    const cs_real_t  eval = exp(xyz[p]+xyz[p+1]+xyz[p+2]-1.5);
    for (int j = 0; j < 9; j++)
      retval[r+j] = eval;
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the name related to a type of function
 *
 * \param[in]  ftype   type of function to deal with
 *
 * \return a pointer to a string
 */
/*----------------------------------------------------------------------------*/

static inline const char *
_get_ftype_name(const _function_type_t   ftype)
{
  switch (ftype) {
  case CONSTANT:
    return "P0";
  case LINEAR:
    return "P1";
  case QUADRATIC:
    return "P2";
  case NON_POLYNOMIAL:
    return "NP";
  default:
    return "";
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Get the right function pointer according to the requested
 *          dimension and functino type
 *
 * \param[in]  dim     dimension of the evaluation (scalar, vector or tensor)
 * \param[in]  ftype   type of function to deal with
 *
 * \return a pointer to a function
 */
/*----------------------------------------------------------------------------*/

static inline cs_analytic_func_t *
_get_func_to_eval(const int                dim,
                  const _function_type_t   ftype)
{
  cs_analytic_func_t *_func = NULL;

  switch (dim) {

  case 1:

    switch (ftype) {
    case CONSTANT:
      _func = _unity;
      break;
    case LINEAR:
      _func = _linear_xyz;
      break;
    case QUADRATIC:
      _func = _quadratic_x2;
      break;
    case NON_POLYNOMIAL:
      _func = _nonpoly;
      break;

    default:
      break; /* Nothing to do */
    }
    break;

  case 3:

    switch (ftype) {
    case CONSTANT:
      _func = _unity_vect;
      break;
    case LINEAR:
      _func = _linear_xyz_vect;
      break;
    case QUADRATIC:
      _func = _quadratic_x2_vect;
      break;
    case NON_POLYNOMIAL:
      _func = _nonpoly_vect;
      break;

    default:
      break; /* Nothing to do */
    }
    break;

  case 9:

    switch (ftype) {
    case CONSTANT:
      _func = _unity_tens;
      break;
    case LINEAR:
      _func = _linear_xyz_tens;
      break;
    case QUADRATIC:
      _func = _quadratic_x2_tens;
      break;
    case NON_POLYNOMIAL:
      _func = _nonpoly_tens;
      break;

    default:
      break; /* Nothing to do */
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid dimension.\n",
              __func__);

  }  /* Switch */

  return _func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a cs_cell_mesh_t structure for a uniform hexahedral cell
 *          of size a
 *
 * \param[in]    a          length of sides
 * \param[in]    cm         pointer to the cs_cell_mesh_t struct. to build
 */
/*----------------------------------------------------------------------------*/

static void
_define_cm_hexa_unif(double            a,
                     cs_cell_mesh_t   *cm)
{
  short int  _v, _e, _f;
  short int  *ids = NULL, *sgn = NULL;
  cs_quant_t  *q = NULL;

  const double  ah = a/2.;

  cm->c_id = 0;
  cm->type = FVM_CELL_HEXA;

  /* Set all quantities */
  cm->flag = CS_CDO_LOCAL_PV |CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PEQ |
    CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FEQ |
    CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FE | CS_CDO_LOCAL_EFQ |
    CS_CDO_LOCAL_DIAM;
  cm->xc[0] = cm->xc[1] = cm->xc[2] = ah;
  cm->vol_c = a*a*a;

  /* VERTICES */
  cm->n_vc = 8;
  for (int i = 0; i < cm->n_vc; i++) {
    cm->v_ids[i] = i;
    cm->wvc[i] = 1./8.;
  }

  /* Coordinates */
  _v = 0; // V0
  cm->xv[3*_v] = 0, cm->xv[3*_v+1] = 0, cm->xv[3*_v+2] = 0;
  _v = 1; // V1
  cm->xv[3*_v] = a, cm->xv[3*_v+1] = 0, cm->xv[3*_v+2] = 0;
  _v = 2; // V2
  cm->xv[3*_v] = a, cm->xv[3*_v+1] = a, cm->xv[3*_v+2] = 0;
  _v = 3; // V3
  cm->xv[3*_v] = 0, cm->xv[3*_v+1] = a, cm->xv[3*_v+2] = 0;
  _v = 4; // V4
  cm->xv[3*_v] = 0, cm->xv[3*_v+1] = 0, cm->xv[3*_v+2] = a;
  _v = 5; // V5
  cm->xv[3*_v] = a, cm->xv[3*_v+1] = 0, cm->xv[3*_v+2] = a;
  _v = 6; // V6
  cm->xv[3*_v] = a, cm->xv[3*_v+1] = a, cm->xv[3*_v+2] = a;
  _v = 7;
  cm->xv[3*_v] = 0, cm->xv[3*_v+1] = a, cm->xv[3*_v+2] = a;

  /* EDGES */
  cm->n_ec = 12;

  // e0
  _e = 0, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 0, ids[1] = 1, sgn[0] = -1;
  q->center[0] = ah, q->center[1] = 0, q->center[2] = 0;
  q->unitv[0] = 1.0, q->unitv[1] = 0.0, q->unitv[2] = 0.0;

  // e1
  _e = 1, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 0, ids[1] = 3, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = 0;
  q->unitv[1] = 1.0, q->center[1] = ah;
  q->unitv[2] = 0.0, q->center[2] = 0;

  // e2
  _e = 2, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 0, ids[1] = 4, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = 0;
  q->unitv[1] = 0.0, q->center[1] = 0;
  q->unitv[2] = 1.0, q->center[2] = ah;

  // e3
  _e = 3, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 1, ids[1] = 2, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = a;
  q->unitv[1] = 1.0, q->center[1] = ah;
  q->unitv[2] = 0.0, q->center[2] = 0;

  // e4
  _e = 4, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 1, ids[1] = 5, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = a;
  q->unitv[1] = 0.0, q->center[1] = 0;
  q->unitv[2] = 1.0, q->center[2] = ah;

  // e5
  _e = 5, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 2, ids[1] = 6, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = a;
  q->unitv[1] = 0.0, q->center[1] = a;
  q->unitv[2] = 1.0, q->center[2] = ah;

  // e6
  _e = 6, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 2, ids[1] = 3, sgn[0] = -1;
  q->unitv[0] = -1.0, q->center[0] = ah;
  q->unitv[1] =  0.0, q->center[1] = a;
  q->unitv[2] =  0.0, q->center[2] = 0;

  // e7
  _e = 7, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 4, ids[1] = 5, sgn[0] = -1;
  q->unitv[0] = 1.0, q->center[0] = ah;
  q->unitv[1] = 0.0, q->center[1] = 0;
  q->unitv[2] = 0.0, q->center[2] = a;

  // e8
  _e = 8; ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 5, ids[1] = 6, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = a;
  q->unitv[1] = 1.0, q->center[1] = ah;
  q->unitv[2] = 0.0, q->center[2] = a;

  // e9
  _e = 9, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 6, ids[1] = 7, sgn[0] = -1;
  q->unitv[0] = -1.0, q->center[0] = ah;
  q->unitv[1] =  0.0, q->center[1] = a;
  q->unitv[2] =  0.0, q->center[2] = a;

  // e10
  _e = 10; ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge +_e;
  ids[0] = 4, ids[1] = 7, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = 0;
  q->unitv[1] = 1.0, q->center[1] = ah;
  q->unitv[2] = 0.0, q->center[2] = a;

  // e11
  _e = 11, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge +_e;
  ids[0] = 3, ids[1] = 7, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = 0;
  q->unitv[1] = 0.0, q->center[1] = a;
  q->unitv[2] = 1.0, q->center[2] = ah;

  for (short int e = 0; e < cm->n_ec; e++) {
    cm->e_ids[e] = e;
    cm->edge[e].meas = a;
    cm->dface[e].meas = ah*ah;
    for (int k = 0; k < 3; k++) cm->dface[e].unitv[k] = cm->edge[e].unitv[k];
  }

  /* FACES */
  cm->n_fc = 6;
  cm->f2e_idx[0] = 0;
  for (short int f = 0; f < cm->n_fc; f++)
    cm->f2e_idx[f+1] = cm->f2e_idx[f] + 4;

  // f0
  _f = 0, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 0, ids[1] = 3, ids[2] = 6, ids[3] = 1;
  q->unitv[0] =  0.0, q->center[0] = ah;
  q->unitv[1] =  0.0, q->center[1] = ah;
  q->unitv[2] = -1.0, q->center[2] = 0;

  // f1
  _f = 1, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 0, ids[1] = 4, ids[2] = 7, ids[3] = 2;
  q->unitv[0] =  0.0, q->center[0] = ah;
  q->unitv[1] = -1.0, q->center[1] = 0;
  q->unitv[2] =  0.0, q->center[2] = ah;

  // f2
  _f = 2, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 3, ids[1] = 5, ids[2] = 8, ids[3] = 4;
  q->unitv[0] =  1.0, q->center[0] = a;
  q->unitv[1] =  0.0, q->center[1] = ah;
  q->unitv[2] =  0.0, q->center[2] = ah;

  // f3
  _f = 3, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 6, ids[1] = 11, ids[2] = 9, ids[3] = 5;
  q->unitv[0] =  0.0, q->center[0] = ah;
  q->unitv[1] =  1.0, q->center[1] = a;
  q->unitv[2] =  0.0, q->center[2] = ah;

  // f4
  _f = 4, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 1, ids[1] = 11, ids[2] = 10, ids[3] = 2;
  q->unitv[0] = -1.0, q->center[0] = 0;
  q->unitv[1] =  0.0, q->center[1] = ah;
  q->unitv[2] =  0.0, q->center[2] = ah;

  // f5
  _f = 5, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 7, ids[1] = 8, ids[2] = 9, ids[3] = 10;
  q->unitv[0] =  0.0, q->center[0] = ah;
  q->unitv[1] =  0.0, q->center[1] = ah;
  q->unitv[2] =  1.0, q->center[2] = a;

  assert(cm->f2e_idx[cm->n_fc] == 24);

  for (short int f = 0; f < cm->n_fc; f++) {
    cm->f_ids[f] = f;
    cm->f_sgn[f] = 1; // all face are outward-oriented
    cm->hfc[f] = ah;
    cm->face[f].meas = a*a;
    cm->dedge[f].meas = ah;
    cm->f_diam[f] = a * sqrt(2);
    for (int k = 0; k < 3; k++) cm->dedge[f].unitv[k] = cm->face[f].unitv[k];
  }

  for (int i = 0; i < cm->f2e_idx[cm->n_fc]; i++)
    cm->tef[i] = ah*ah;

  cm->diam_c = a * sqrt(3);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a cs_cell_mesh_t structure for a tetrahedral cell
 *          of size a
 *
 * \param[in]    a          length of sides
 * \param[in]    cm         pointer to the cs_cell_mesh_t struct. to build
 */
/*----------------------------------------------------------------------------*/

static void
_define_cm_tetra_ref(double            a,
                     cs_cell_mesh_t   *cm)
{
  short int  _v, _e, _f;
  short int  *ids = NULL, *sgn = NULL;
  cs_quant_t  *q = NULL;

  const double  ah = a/2.;
  const double  sq2 = sqrt(2.), invsq2 = 1./sq2;

  cm->c_id = 0;
  cm->type = FVM_CELL_TETRA;

  /* Set all quantities */
  cm->flag = CS_CDO_LOCAL_PV |CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PEQ |
    CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FEQ |
    CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FE |CS_CDO_LOCAL_EFQ  |
    CS_CDO_LOCAL_DIAM;

  cm->vol_c = cs_math_onesix*a*a*a;
  cm->xc[0] = cm->xc[1] = cm->xc[2] = 0.25*a;

  /* VERTICES */
  cm->n_vc = 4;
  for (int i = 0; i < cm->n_vc; i++) {
    cm->v_ids[i] = i;
    cm->wvc[i] = 0;
  }

  /* Coordinates */
  _v = 0; // V0
  cm->xv[3*_v] = 0, cm->xv[3*_v+1] = 0, cm->xv[3*_v+2] = 0;
  _v = 1; // V1
  cm->xv[3*_v] = a, cm->xv[3*_v+1] = 0, cm->xv[3*_v+2] = 0;
  _v = 2; // V2
  cm->xv[3*_v] = 0, cm->xv[3*_v+1] = a, cm->xv[3*_v+2] = 0;
  _v = 3; // V3
  cm->xv[3*_v] = 0, cm->xv[3*_v+1] = 0, cm->xv[3*_v+2] = a;

  /* EDGES */
  cm->n_ec = 6;
  for (short int e = 0; e < cm->n_ec; e++) cm->e_ids[e] = e;

  // e0
  _e = 0, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 0, ids[1] = 1, sgn[0] = -1;
  q->center[0] = ah, q->center[1] = 0, q->center[2] = 0;
  q->unitv[0] = 1.0, q->unitv[1] = 0.0, q->unitv[2] = 0.0;
  q->meas = a;

  // e1
  _e = 1, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 0, ids[1] = 2, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = 0;
  q->unitv[1] = 1.0, q->center[1] = ah;
  q->unitv[2] = 0.0, q->center[2] = 0;
  q->meas = a;

  // e2
  _e = 2, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 0, ids[1] = 3, sgn[0] = -1;
  q->unitv[0] = 0.0, q->center[0] = 0;
  q->unitv[1] = 0.0, q->center[1] = 0;
  q->unitv[2] = 1.0, q->center[2] = ah;
  q->meas = a;

  // e3
  _e = 3, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 1, ids[1] = 2, sgn[0] = -1;
  q->unitv[0] =-invsq2, q->center[0] = ah;
  q->unitv[1] = invsq2, q->center[1] = ah;
  q->unitv[2] =    0.0, q->center[2] = 0;
  q->meas = a * sq2;

  // e4
  _e = 4, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 1, ids[1] = 3, sgn[0] = -1;
  q->unitv[0] =-invsq2, q->center[0] = ah;
  q->unitv[1] =    0.0, q->center[1] = 0;
  q->unitv[2] = invsq2, q->center[2] = ah;
  q->meas = a * sq2;

  // e5
  _e = 5, ids = cm->e2v_ids + 2*_e; sgn = cm->e2v_sgn + _e, q = cm->edge + _e;
  ids[0] = 2, ids[1] = 3, sgn[0] = -1;
  q->unitv[0] =    0.0, q->center[0] = 0;
  q->unitv[1] =-invsq2, q->center[1] = ah;
  q->unitv[2] = invsq2, q->center[2] = ah;
  q->meas = a * sq2;

  /* FACES */
  cm->n_fc = 4;
  for (short int f = 0; f < cm->n_fc; f++) {
    cm->f_ids[f] = f;
    cm->f_sgn[f] = 1; // all face are outward-oriented
  }

  cm->f2e_idx[0] = 0;
  for (short int f = 0; f < cm->n_fc; f++)
    cm->f2e_idx[f+1] = cm->f2e_idx[f] + 3;

  // f0
  _f = 0, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 0, ids[1] = 3, ids[2] = 1;
  q->unitv[0] =  0.0, q->center[0] = a/3.;
  q->unitv[1] =  0.0, q->center[1] = a/3.;
  q->unitv[2] = -1.0, q->center[2] = 0;
  q->meas = a*ah;

  // f1
  _f = 1, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 0, ids[1] = 4, ids[2] = 2;
  q->unitv[0] =  0.0, q->center[0] = a/3.;
  q->unitv[1] = -1.0, q->center[1] = 0;
  q->unitv[2] =  0.0, q->center[2] = a/3.;
  q->meas = a*ah;

  // f2
  _f = 2, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 1, ids[1] = 5, ids[2] = 2;
  q->unitv[0] = -1.0, q->center[0] = 0;
  q->unitv[1] =  0.0, q->center[1] = a/3.;
  q->unitv[2] =  0.0, q->center[2] = a/3.;
  q->meas = a*ah;

  // f3
  _f = 3, ids = cm->f2e_ids + cm->f2e_idx[_f], q = cm->face + _f;
  ids[0] = 3, ids[1] = 5, ids[2] = 4;
  q->unitv[0] = 1/sqrt(3), q->center[0] = a/3.;
  q->unitv[1] = 1/sqrt(3), q->center[1] = a/3.;
  q->unitv[2] = 1/sqrt(3), q->center[2] = a/3.;
  q->meas = 0.5*sqrt(3)*a*a;

  assert(cm->f2e_idx[cm->n_fc] == 12);

  /* Dual faces, wvc ? */

  /* Compute additional quantities */
  for (short int i = 0; i < 2*cm->n_ec; i++) cm->e2f_ids[i] = -1;

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];

    /* Compute dual edge quantities */
    cs_math_3_length_unitv(cm->xc, pfq.center,
                           &(cm->dedge[f].meas), cm->dedge[f].unitv);

    /* Compute height of the pyramid of basis f */
    cm->hfc[f] = cs_math_3_dot_product(pfq.unitv,
                                       cm->dedge[f].unitv)*cm->dedge[f].meas;
    assert(cm->hfc[f] > 0);

    /* Compute tef */
    for (short int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      cs_nvec3_t  sefc;
      cs_real_3_t  cp_efc, xexf, xexc;

      const short int  e = cm->f2e_ids[i], eshft = 2*e;
      const cs_quant_t  peq = cm->edge[e]; /* Edge quantities */

      cm->tef[i] = cs_compute_area_from_quant(peq, pfq.center);

      /* Compute the vectorial area for the triangle : xc, xf, xe */
      for (int k = 0; k < 3; k++) {
        xexf[k] = pfq.center[k] - peq.center[k];
        xexc[k] = cm->xc[k] - peq.center[k];
      }
      cs_math_3_cross_product(xexf, xexc, cp_efc);
      cs_nvec3(cp_efc, &sefc);

      /* One should have (cp_efc, sefc) > 0 */
      short int  _sgn = 1;
      if (_dp3(sefc.unitv, peq.unitv) < 0) _sgn = -1;

      if (cm->e2f_ids[eshft] == -1) {
        cm->e2f_ids[eshft] = f;
        cm->sefc[eshft].meas = 0.5*sefc.meas;
        for (int k = 0; k < 3; k++)
          cm->sefc[eshft].unitv[k] = _sgn*sefc.unitv[k];
      }
      else {
        assert(cm->e2f_ids[eshft+1] == -1);
        cm->e2f_ids[eshft+1] = f;
        cm->sefc[eshft+1].meas = 0.5*sefc.meas;
        for (int k = 0; k < 3; k++)
          cm->sefc[eshft+1].unitv[k] = _sgn*sefc.unitv[k];
      }

    }

  }  /* Loop on cell faces */

  /* Compute dual face quantities */
  for (short int e = 0; e < cm->n_ec; e++) {

    cs_real_3_t  df;
    const cs_nvec3_t  s1 = cm->sefc[2*e], s2 = cm->sefc[2*e+1];
    for (int k = 0; k < 3; k++)
      df[k] = s1.meas*s1.unitv[k] + s2.meas*s2.unitv[k];
    cs_nvec3(df, &(cm->dface[e]));

  }  /* Loop on cell edges */

  /* Compute dual cell volume */
  for (short int f = 0; f < cm->n_fc; f++) {

    const double  hf_coef = cs_math_onesix * cm->hfc[f];

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];
      const double  half_pef_vol = cm->tef[i]*hf_coef;

      cm->wvc[v1] += half_pef_vol;
      cm->wvc[v2] += half_pef_vol;

    }  /* Loop on face edges */

  }  /* Loop on cell faces */

  /* Reset diam */
  double  dbuf[10];
  short int  vtag[4];
  int  size = cm->n_vc*(cm->n_vc+1)/2;
  int  shift = 0;

  cm->diam_c = -1;
  for (int i = 0; i < size; i++) dbuf[i] = 0.;

  for (short int vi = 0; vi < cm->n_vc; ++vi) {
    shift++;  /* diag entry not taken into account */
    const double *xvi = cm->xv + 3*vi;
    for (short int vj = vi+1; vj < cm->n_vc; vj++, shift++) {
      double  l = dbuf[shift] = cs_math_3_distance(xvi, cm->xv + 3*vj);
      if (l > cm->diam_c) cm->diam_c = l;

    } /* Loop on vj > vi */
  }   /* Loop on vi */

  for (short int f = 0; f < cm->n_fc; ++f) {

    /* Reset vtag */
    for (short int v = 0; v < cm->n_vc; v++) vtag[v] = -1;

    /* Tag face vertices */
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const int  eshft = 2*cm->f2e_ids[i];
      vtag[cm->e2v_ids[eshft  ]] = 1;
      vtag[cm->e2v_ids[eshft+1]] = 1;
    }

    cm->f_diam[f] = -1;
    for (short int vi = 0; vi < cm->n_vc; ++vi) {

      if (vtag[vi] > 0) { /* belong to the current face */
        for (short int vj = vi+1; vj < cm->n_vc; vj++) {
          if (vtag[vj] > 0) { /* belong to the current face */

            shift = vj*(vj+1)/2 + vi;
            const double  l = dbuf[shift];
            if (l > cm->f_diam[f]) cm->f_diam[f] = l;

          }
        } /* Loop on vj > vi */
      }
    }   /* Loop on vi */

  } /* Loop on cell faces */

  const double  invvol = 1/cm->vol_c;
  for (short int v = 0; v < cm->n_vc; v++) cm->wvc[v] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test quadrature functions on scalar, vector or tensor valued
 *          used for computing source terms in case of CDO-Fb schemes
 *
 * \param[in]      out         output file
 * \param[in]      tag         tag
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      dim         dimension of the variable
 * \param[in]      ftype       type of function to deal with
 * \param[in]      val0        first set of values
 * \param[in]      val1        second set of values
 * \param[in]      val2        third set of values
 */
/*----------------------------------------------------------------------------*/

static void
_dump_quad_res(FILE                     *out,
               const char                tag[],
               const cs_cell_mesh_t     *cm,
               const short int           dim,
               const _function_type_t    ftype,
               const cs_real_t          *val0,
               const cs_real_t          *val1,
               const cs_real_t          *val2)
{
  const short int  nf = cm->n_fc;
  const cs_real_t *const  c_val0 = val0 + dim*nf,
                  *const  c_val1 = val1 + dim*nf,
                  *const  c_val2 = val2 + dim*nf;

  fprintf(out, "\n QUADRATURE; %s; %s; DIM = %d\n",
          tag, _get_ftype_name(ftype), dim);
  fprintf(out, " T(ID) ID %12s %12s %12s\n", "LOWEST", "HIGHER", "HIGHEST");

  for (short int f = 0; f < nf; f++) {
    for ( short int i = 0; i < dim; i++)
      fprintf(out, " F(%2d) %2d %10.6e %10.6e %10.6e\n", f, i,
              val0[f*dim+i], val1[f*dim+i], val2[f*dim+i]);
  }
  fprintf(out, "------------------------------------------------\n");

  for ( short int i = 0; i < dim; i++)
    fprintf(out, " C(%2d) %2d %10.6e %10.6e %10.6e\n", 0, i,
            c_val0[i], c_val1[i], c_val2[i]);

  if (ftype == NON_POLYNOMIAL) {

    const short int cell_type = (cm->type == FVM_CELL_TETRA);
    const cs_real_t cex_res = non_poly_int[cell_type] / cm->vol_c;
    cs_real_t fex_res[6];
    for (short int f = 0; f < nf; f++)
      fex_res[f] = non_poly_int_f[cell_type][f] / cm->face[f].meas;

    fprintf(out, " QUADRATURE; %s; %s |ERR|; DIM = %d\n",
            tag, _get_ftype_name(ftype), dim);
    fprintf(out, " T(ID) ID %12s %12s %12s\n", "LOWEST", "HIGHER", "HIGHEST");
    for (short int f = 0; f < nf; f++) {
      for ( short int i = 0; i < dim; i++)
        fprintf(out, " F(%2d) %2d %10.6e %10.6e %10.6e\n", f, i,
                fabs(val0[f*dim+i] - fex_res[f]),
                fabs(val1[f*dim+i] - fex_res[f]),
                fabs(val2[f*dim+i] - fex_res[f]));
      fprintf(out, "------------------------------------------------\n");
    }

    for ( short int i = 0; i < dim; i++)
      fprintf(out, " C(%2d) %2d %10.6e %10.6e %10.6e\n", 0, i,
              fabs(c_val0[i] - cex_res),
              fabs(c_val1[i] - cex_res),
              fabs(c_val2[i] - cex_res));


  } /* Non polynomial function */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test quadrature functions on scalar, vector or tensor valued
 *          used for computing source terms in case of CDO-Fb schemes
 *
 * \param[in]      out         output file
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      dim         dimension of the variable
 * \param[in]      ftype       type of function to deal with
 * \param[in]      n_runs      number of executions to perform
 */
/*----------------------------------------------------------------------------*/

static void
_test_cdofb_source(FILE                     *out,
                   const cs_cell_mesh_t     *cm,
                   const short int           dim,
                   const _function_type_t    ftype,
                   const int                 n_runs)
{
  const short int  totdof =  dim*(cm->n_fc + 1);

  cs_real_t  teval = time_step->t_cur;
  cs_flag_t  state_flag = CS_FLAG_STATE_DENSITY;
  cs_flag_t  meta_flag = cs_source_term_set_default_flag(CS_SPACE_SCHEME_CDOFB);
  cs_cell_builder_t  *cb = NULL;
  cs_cell_sys_t  *csys = NULL;

  if (dim == 1)
    cs_cdofb_scaleq_get(&csys, &cb);
  else
    cs_cdofb_vecteq_get(&csys, &cb);

  cs_real_t  *st0 = NULL, *st1 = NULL, *st2 = NULL, *st3 = NULL;
  BFT_MALLOC(st0, totdof, cs_real_t);
  BFT_MALLOC(st1, totdof, cs_real_t);
  BFT_MALLOC(st2, totdof, cs_real_t);
  BFT_MALLOC(st3, totdof, cs_real_t);

  cs_real_t *const  c_st0 = st0 + dim*cm->n_fc;
  cs_real_t *const  c_st1 = st1 + dim*cm->n_fc;
  cs_real_t *const  c_st2 = st2 + dim*cm->n_fc;
  cs_real_t *const  c_st3 = st3 + dim*cm->n_fc;

  /* Evaluate the performance */
  cs_timer_t  t0, t1, t2, t3, t4;
  cs_timer_counter_t  tc0, tc1, tc2, tc3;
  CS_TIMER_COUNTER_INIT(tc0);
  CS_TIMER_COUNTER_INIT(tc1);
  CS_TIMER_COUNTER_INIT(tc2);
  CS_TIMER_COUNTER_INIT(tc3);

  if (ftype == CONSTANT) { /* definition by value */

    /* Constant value equal to 1 */
    cs_real_t  ones[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    cs_xdef_t  *st =
      cs_xdef_volume_create(CS_XDEF_BY_VALUE, dim, 0,  /* z_id */
                            state_flag | CS_FLAG_STATE_UNIFORM, meta_flag,
                            ones);

    /* Loop on runs to evaluate the performance of each quadrature */
    for (int r = 0; r < n_runs; r++) {

      /* Reset cell values */
      memset(c_st0, 0, dim*sizeof(cs_real_t));

      switch (dim) {

      case 1:
        t0 = cs_timer_time();
        cs_source_term_pcsd_by_value(st, cm, teval, cb, NULL, st0);
        t1 = cs_timer_time();
        break;

      case 3:
        t0 = cs_timer_time();
        cs_source_term_pcvd_by_value(st, cm, teval, cb, NULL, st0);
        t1 = cs_timer_time();
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Need to implement cs_source_term_pctd_by_value.\n",
                  __func__);
        break;
      }

      cs_timer_counter_add_diff(&(tc0), &t0, &t1);

    }

    fprintf(out, "\nCDO.FB; SOURCE_TERM %s; DIM = %d\n",
            _get_ftype_name(ftype), dim);
    fprintf(out, " CS %12s\n", "FBSD");
    for ( short int i = 0; i < dim; i++)
      fprintf(out, "%3d %10.6e\n", i, c_st0[i]);

    st = cs_xdef_free(st);

  }
  else { /* definition by analytic */

    /* In what follows, when computing the exact integral, a unitary cube/hexa
     * was considered. The precision results will be incorrect if another
     * cell-type is considered instead */

    cs_analytic_func_t  *_func = _get_func_to_eval(dim, ftype);
    if (_func == NULL)
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.\n", __func__);

    cs_xdef_analytic_input_t  anai = {.func = _func, .input = NULL };
    cs_xdef_t  *st = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                           dim,
                                           0,  /* z_id */
                                           state_flag, meta_flag,
                                           &anai);

    /* Loop on runs to evaluate the performance of each quadrature */
    for (int r = 0; r < n_runs; r++) {

      /* Reset */
      memset(c_st0, 0, dim*sizeof(cs_real_t));
      memset(c_st1, 0, dim*sizeof(cs_real_t));
      memset(c_st2, 0, dim*sizeof(cs_real_t));
      memset(c_st3, 0, dim*sizeof(cs_real_t));

      switch (dim) {

      case 1:
        t0 = cs_timer_time();
        cs_source_term_pcsd_bary_by_analytic(st, cm, teval, cb, NULL, st0);
        t1 = cs_timer_time();
        cs_xdef_set_quadrature(st, CS_QUADRATURE_BARY_SUBDIV);
        cs_source_term_pcsd_by_analytic(st, cm, teval, cb, NULL, st1);
        t2 = cs_timer_time();
        cs_xdef_set_quadrature(st, CS_QUADRATURE_HIGHER);
        cs_source_term_pcsd_by_analytic(st, cm, teval, cb, NULL, st2);
        t3 = cs_timer_time();
        cs_xdef_set_quadrature(st, CS_QUADRATURE_HIGHEST);
        cs_source_term_pcsd_by_analytic(st, cm, teval, cb, NULL, st3);
        t4 = cs_timer_time();
        break;

      case 3:
        t0 = cs_timer_time();
        cs_source_term_pcvd_bary_by_analytic(st, cm, teval, cb, NULL, st0);
        t1 = cs_timer_time();
        cs_xdef_set_quadrature(st, CS_QUADRATURE_BARY_SUBDIV);
        cs_source_term_pcvd_by_analytic(st, cm, teval, cb, NULL, st1);
        t2 = cs_timer_time();
        cs_xdef_set_quadrature(st, CS_QUADRATURE_HIGHER);
        cs_source_term_pcvd_by_analytic(st, cm, teval, cb, NULL, st2);
        t3 = cs_timer_time();
        cs_xdef_set_quadrature(st, CS_QUADRATURE_HIGHEST);
        cs_source_term_pcvd_by_analytic(st, cm, teval, cb, NULL, st3);
        t4 = cs_timer_time();
        break;

      case 9:
        /*_source_term = cs_source_term_fbta_by_analytic; --> TODO */
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Need to implement cs_source_term_fbta_by_analytic.\n",
                  __func__);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, " %s: Invalid dimension.\n",
                  __func__);

      }  /* Switch */

      cs_timer_counter_add_diff(&(tc0), &t0, &t1);
      cs_timer_counter_add_diff(&(tc1), &t1, &t2);
      cs_timer_counter_add_diff(&(tc2), &t2, &t3);
      cs_timer_counter_add_diff(&(tc3), &t3, &t4);

    } /* Loop on runs */

    fprintf(out, "\nCDO.FB; SOURCE_TERM %s; DIM = %d\n",
            _get_ftype_name(ftype), dim);
    fprintf(out, "  C %12s %12s %12s %12s\n",
            "BARY", "BARY.SUB", "HIGHER", "HIGHEST");
    for ( short int i = 0; i < dim; i++)
      fprintf(out, "%3d %10.6e %10.6e %10.6e %10.6e\n", i,
              c_st0[i], c_st1[i], c_st2[i], c_st3[i]);

    st = cs_xdef_free(st);

  } /* Definition by analytic function */

  fprintf(out, "\nCDO.FB; PERFORMANCE OF SOURCE TERMS %s; DIM = %d\n",
          _get_ftype_name(ftype), dim);
  fprintf(out, " %12s %12s %12s %12s\n",
          "BARY", "BARY.SUB", "HIGHER", "HIGHEST");
  fprintf(out, " %10.6e %10.6e %10.6e %10.6e\n",
          tc0.wall_nsec*1e-9, tc1.wall_nsec*1e-9, tc2.wall_nsec*1e-9,
          tc3.wall_nsec*1e-9);

  BFT_FREE(st0);
  BFT_FREE(st1);
  BFT_FREE(st2);
  BFT_FREE(st3);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test CDO vertex-based schemes
 *
 * \param[in]      out         output file
 * \param[in]      cm          pointer to a cs_cell_mesh_t structure
 * \param[in]      dim         dimension of the variable
 * \param[in]      ftype       type of function to deal with
 * \param[in]      n_runs      number of executions to perform
 */
/*----------------------------------------------------------------------------*/

static void
_test_cdovb_source(FILE                     *out,
                   const cs_cell_mesh_t     *cm,
                   const short int           dim,
                   const _function_type_t    ftype,
                   const int                 n_runs)
{
  const double  tcur = 0.;

  cs_real_t  st0[8], st1[8], st2[8], st3[8];
  cs_flag_t  state_flag = CS_FLAG_STATE_DENSITY;
  cs_flag_t  meta_flag = cs_source_term_set_default_flag(CS_SPACE_SCHEME_CDOVB);

  cs_cell_builder_t  *cb = NULL;
  cs_cell_sys_t  *csys = NULL;

  cs_cdovb_scaleq_get(&csys, &cb);

  /* Evaluate the performance */
  cs_timer_t  t0, t1, t2, t3, t4;
  cs_timer_counter_t  tc0, tc1, tc2, tc3;
  CS_TIMER_COUNTER_INIT(tc0);  /* build system */
  CS_TIMER_COUNTER_INIT(tc1);  /* build system */
  CS_TIMER_COUNTER_INIT(tc2);  /* build system */
  CS_TIMER_COUNTER_INIT(tc3);  /* build system */

  if (ftype == CONSTANT) { /* definition by value */

    /* Constant value equal to 1 */
    cs_real_t  ones[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    cs_xdef_t  *st =
      cs_xdef_volume_create(CS_XDEF_BY_VALUE, dim, 0,  /* z_id */
                            state_flag | CS_FLAG_STATE_UNIFORM,
                            meta_flag,
                            ones);

    /* Loop on runs to evaluate the performance of each quadrature */
    for (int r = 0; r < n_runs; r++) {

      /* Reset */
      for (int v = 0; v < cm->n_vc; v++)
        st0[v] = st1[v] = st2[v] = st3[v] = 0.0;

      switch (dim) {

      case 1:
        t0 = cs_timer_time();
        cs_source_term_dcsd_by_value(st, cm, tcur, cb, NULL, st0);
        t1 = cs_timer_time();
        break;

      case 3:
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Need to implement cs_source_term_fbtd_by_value.\n",
                  __func__);
        break;
      }

      cs_timer_counter_add_diff(&(tc0), &t0, &t1);

    }

    fprintf(out, "\nCDO.VB; SOURCE_TERM %s; DIM = %d\n",
            _get_ftype_name(ftype), dim);
    fprintf(out, " V %12s\n", "DCSD");
    for (int i = 0; i < cm->n_vc; i++)
      fprintf(out, "%2d %10.6e\n", i, st0[i]);

    st = cs_xdef_free(st);

  }
  else { /* Definition by analytic */

    cs_xdef_analytic_input_t  anai = {.func = _get_func_to_eval(dim, ftype),
                                      .input = NULL };

    cs_xdef_t  *st = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                           dim,
                                           0,  /* z_id */
                                           state_flag,
                                           meta_flag,
                                           &anai);

    /* Loop on runs to evaluate the performance of each quadrature */
    for (int r = 0; r < n_runs; r++) {

      /* Reset */
      for (int v = 0; v < cm->n_vc; v++)
        st0[v] = st1[v] = st2[v] = st3[v] = 0.0;

      switch (dim) {

      case 1:
        t0 = cs_timer_time();
        cs_source_term_dcsd_bary_by_analytic(st, cm, tcur, cb, NULL, st0);
        t1 = cs_timer_time();
        cs_source_term_dcsd_q1o1_by_analytic(st, cm, tcur, cb, NULL, st1);
        t2 = cs_timer_time();
        cs_source_term_dcsd_q10o2_by_analytic(st, cm, tcur, cb, NULL, st2);
        t3 = cs_timer_time();
        cs_source_term_dcsd_q5o3_by_analytic(st, cm, tcur, cb, NULL, st3);
        t4 = cs_timer_time();
        break;

      case 3:
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Need to implement cs_source_term_fbtd_by_value.\n",
                  __func__);
        break;
      }

      cs_timer_counter_add_diff(&(tc0), &t0, &t1);
      cs_timer_counter_add_diff(&(tc1), &t1, &t2);
      cs_timer_counter_add_diff(&(tc2), &t2, &t3);
      cs_timer_counter_add_diff(&(tc3), &t3, &t4);

    }

    fprintf(out, "\nCDO.VB; SOURCE_TERM P1\n");
    fprintf(out, " V %12s %12s %12s %12s\n",
            "DCSD_BARY", "DCSD_Q1O1", "DCSD_Q10O2", "DCSD_Q5O3");
    for (int i = 0; i < cm->n_vc; i++)
      fprintf(out, "%2d %10.6e %10.6e %10.6e %10.6e\n",
              i, st0[i], st1[i], st2[i], st3[i]);

    if (ftype == NON_POLYNOMIAL && cm->n_vc == 8) {

      cs_real_t  exact_result[8] = {0.0609162, // V (0.0,0.0,0.0)
                                    0.100434,  // V (1.0,0.0,0.0)
                                    0.165587,  // V (1.0,1.0,0.0)
                                    0.100434,  // V (0.0,1.0,0.0)
                                    0.100434,  // V (0.0,0.0,1.0)
                                    0.165587,  // V (1.0,0.0,1.0)
                                    0.273007,  // V (1.0,1.0,1.0)
                                    0.165587}; // V (0.0,1.0,1.0)

      fprintf(out, " V %12s %12s %12s %12s (ERROR)\n",
              "DCSD_BARY", "DCSD_Q1O1", "DCSD_Q10O2", "DCSD_Q5O3");
      for (int i = 0; i < cm->n_vc; i++)
        fprintf(out, "%2d % 10.6e % 10.6e % 10.6e % 10.6e\n", i,
                st0[i] - exact_result[i],
                st1[i] - exact_result[i],
                st2[i] - exact_result[i],
                st3[i] - exact_result[i]);
    }

    st = cs_xdef_free(st);

  } /* Definition by analytic */

  fprintf(out, "\nCDO.VB; PERFORMANCE OF SOURCE TERMS\n");
  fprintf(out, " %12s %12s %12s %12s\n",
          "DCSD_BARY", "DCSD_Q1O1", "DCSD_Q10O2", "DCSD_Q5O3");
  fprintf(out, " %10.6e %10.6e %10.6e %10.6e\n",
          tc0.wall_nsec*1e-9, tc1.wall_nsec*1e-9,
          tc2.wall_nsec*1e-9, tc3.wall_nsec*1e-9);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test quadrature functions on scalar, vector or tensor valued
 *          analytic functions
 *
 * \param[in]  out     output file
 * \param[in]  cm      pointer to a cs_cell_mesh_t structure
 * \param[in]  dim     dimension of the evaluation (scalar, vector or tensor)
 * \param[in]  ftype   type of function to deal with
 * \param[in]  n_runs  number of execution to perform to assess the performance
 */
/*----------------------------------------------------------------------------*/

static void
_test_quadratures_misc(FILE                     *out,
                       const cs_cell_mesh_t     *cm,
                       const short int           dim,
                       const _function_type_t    ftype,
                       const int                 n_runs)
{
  const short int  nf = cm->n_fc;
  const short int  totdof = dim*(nf + 1);
  const cs_real_t  time = time_step->t_cur;

  cs_analytic_func_t *_func = _get_func_to_eval(dim, ftype);
  if (_func == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.\n", __func__);

  void * f_input = NULL;
  cs_real_t  *st0, *st1, *st2;

  BFT_MALLOC(st0, totdof, cs_real_t);
  BFT_MALLOC(st1, totdof, cs_real_t);
  BFT_MALLOC(st2, totdof, cs_real_t);

  cs_real_t *const c_st0 = st0 + dim*nf,
            *const c_st1 = st1 + dim*nf,
            *const c_st2 = st2 + dim*nf;

  /* Evaluate the performance */
  cs_timer_counter_t  tc0, tc1, tc2;
  CS_TIMER_COUNTER_INIT(tc0);
  CS_TIMER_COUNTER_INIT(tc1);
  CS_TIMER_COUNTER_INIT(tc2);

  /* Set up quadrature related variables */
  cs_quadrature_tetra_integral_t
    *tet_bary = cs_quadrature_get_tetra_integral(dim, CS_QUADRATURE_BARY),
    *tet_i_er = cs_quadrature_get_tetra_integral(dim, CS_QUADRATURE_HIGHER),
    *tet_i_st = cs_quadrature_get_tetra_integral(dim, CS_QUADRATURE_HIGHEST);
  cs_quadrature_tria_integral_t
    *tri_bary = cs_quadrature_get_tria_integral(dim, CS_QUADRATURE_BARY),
    *tri_i_er = cs_quadrature_get_tria_integral(dim, CS_QUADRATURE_HIGHER),
    *tri_i_st = cs_quadrature_get_tria_integral(dim, CS_QUADRATURE_HIGHEST);

  /* Loop on runs to evaluate the performance of each quadrature */
  for (int r = 0; r < n_runs; r++) {

    /* Reset values */
    memset(st0, 0, totdof*sizeof(cs_real_t));
    memset(st1, 0, totdof*sizeof(cs_real_t));
    memset(st2, 0, totdof*sizeof(cs_real_t));

    cs_timer_t  t0 = cs_timer_time();
    cs_xdef_cw_eval_fc_int_by_analytic(cm, time,
                                       _func, f_input,
                                       dim,
                                       tet_bary, tri_bary,
                                       c_st0, st0);
    cs_timer_t  t1 = cs_timer_time();
    cs_xdef_cw_eval_fc_int_by_analytic(cm, time,
                                       _func, f_input,
                                       dim, tet_i_er, tri_i_er,
                                       c_st1, st1);

    cs_timer_t  t2 = cs_timer_time();
    cs_xdef_cw_eval_fc_int_by_analytic(cm, time,
                                       _func, f_input,
                                       dim, tet_i_st, tri_i_st,
                                       c_st2, st2);
    cs_timer_t  t3 = cs_timer_time();

    cs_timer_counter_add_diff(&(tc0), &t0, &t1);
    cs_timer_counter_add_diff(&(tc1), &t1, &t2);
    cs_timer_counter_add_diff(&(tc2), &t2, &t3);

  }

  /* Dump performance and evaluations */
  _dump_quad_res(out, "CELL&FACE", cm, dim, ftype,
                 st0, st1, st2);

  fprintf(out, "\nQUADRATURE; PERFORMANCE OF QUADRATURES; %s; DIM = %d\n",
          _get_ftype_name(ftype), dim);
  fprintf(out, " %12s %12s %12s\n", "LOWEST", "HIGHER", "HIGHEST");
  fprintf(out, " %10.6e %10.6e %10.6e\n",
          tc0.wall_nsec*1e-9, tc1.wall_nsec*1e-9, tc2.wall_nsec*1e-9);

  BFT_FREE(st0);
  BFT_FREE(st1);
  BFT_FREE(st2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test quadrature functions on scalar, vector or tensor valued
 *          analytic functions
 *
 * \param[in]  out     output file
 * \param[in]  cm      pointer to a cs_cell_mesh_t structure
 * \param[in]  dim     dimension of the evaluation (scalar, vector or tensor)
 * \param[in]  ftype   type of function to deal with
 * \param[in]  n_runs  number of execution to perform to assess the performance
 */
/*----------------------------------------------------------------------------*/

static void
_test_quadratures_xdef(FILE                     *out,
                       const cs_cell_mesh_t     *cm,
                       const short int           dim,
                       const _function_type_t    ftype,
                       const int                 n_runs)
{
  const cs_real_t  teval = time_step->t_cur;
  const short int nf = cm->n_fc;
  const short int totdof = dim*(nf + 1);

  cs_analytic_func_t  *_func = _get_func_to_eval(dim, ftype);
  if (_func == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.\n", __func__);

  cs_xdef_analytic_input_t  anai = {.func = _func, .input = NULL };

  cs_real_t  *st0, *st1, *st2;
  BFT_MALLOC(st0, totdof, cs_real_t);
  BFT_MALLOC(st1, totdof, cs_real_t);
  BFT_MALLOC(st2, totdof, cs_real_t);
  cs_real_t *const c_st0 = st0 + dim*nf,
            *const c_st1 = st1 + dim*nf,
            *const c_st2 = st2 + dim*nf;

  /* Evaluate the performance */
  cs_timer_counter_t  tc0, tc1, tc2;
  CS_TIMER_COUNTER_INIT(tc0);
  CS_TIMER_COUNTER_INIT(tc1);
  CS_TIMER_COUNTER_INIT(tc2);

  cs_xdef_cw_eval_int_t  *cell_int = NULL;
  cs_xdef_cw_eval_face_t  *face_int = NULL;

  switch (dim) {
  case 1:
    cell_int = cs_xdef_cw_eval_scalar_avg_by_analytic;
    face_int = cs_xdef_cw_eval_scalar_face_avg_by_analytic;
    break;
  case 3:
    cell_int = cs_xdef_cw_eval_vector_avg_by_analytic;
    face_int = cs_xdef_cw_eval_vector_face_avg_by_analytic;
    break;
  case 9:
    cell_int = cs_xdef_cw_eval_tensor_avg_by_analytic;
    face_int = cs_xdef_cw_eval_tensor_face_avg_by_analytic;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid dimension.\n",
              __func__);

  } /* Switch */

  /* Loop on runs to evaluate the performance of each quadrature */
  for (int r = 0; r < n_runs; r++) {

    /* Reset values */
    memset(st0, 0, totdof*sizeof(cs_real_t));
    memset(st1, 0, totdof*sizeof(cs_real_t));
    memset(st2, 0, totdof*sizeof(cs_real_t));

    cs_timer_t  t0 = cs_timer_time();
    cell_int(cm, teval, (void*)(&anai), CS_QUADRATURE_BARY, c_st0);
    for (short int f = 0; f < nf; f++)
      face_int(cm, f, teval, (void*)(&anai), CS_QUADRATURE_BARY,
               st0 + f*dim);
    cs_timer_t  t1 = cs_timer_time();

    cell_int(cm, teval, (void*)(&anai), CS_QUADRATURE_HIGHER, c_st1);
    for (short int f = 0; f < nf; f++)
      face_int(cm, f, teval, (void*)(&anai), CS_QUADRATURE_HIGHER,
               st1 + f*dim);
    cs_timer_t  t2 = cs_timer_time();

    cell_int(cm, teval, (void*)(&anai), CS_QUADRATURE_HIGHEST, c_st2);
    for (short int f = 0; f < nf; f++)
      face_int(cm, f, teval, (void*)(&anai), CS_QUADRATURE_HIGHEST,
               st2 + f*dim);
    cs_timer_t  t3 = cs_timer_time();

    cs_timer_counter_add_diff(&(tc0), &t0, &t1);
    cs_timer_counter_add_diff(&(tc1), &t1, &t2);
    cs_timer_counter_add_diff(&(tc2), &t2, &t3);

  } /* Loop on runs */

  _dump_quad_res(out, "XDEF AVG", cm, dim, ftype,
                 st0, st1, st2);

  fprintf(out, "\nQUADRATURE; PERFORMANCE OF QUADRATURES; %s; DIM = %d\n",
          _get_ftype_name(ftype), dim);
  fprintf(out, " %12s %12s %12s\n", "LOWEST", "HIGHER", "HIGHEST");
  fprintf(out, " %10.6e %10.6e %10.6e\n",
          tc0.wall_nsec*1e-9, tc1.wall_nsec*1e-9, tc2.wall_nsec*1e-9);

  BFT_FREE(st0);
  BFT_FREE(st1);
  BFT_FREE(st2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test quadrature functions on scalar, vector or tensor valued
 *          analytic functions
 *
 * \param[in]  out     output file
 * \param[in]  cm      pointer to a cs_cell_mesh_t structure
 * \param[in]  dim     dimension of the evaluation (scalar or vector)
 * \param[in]  ftype   type of function to deal with
 */
/*----------------------------------------------------------------------------*/

static void
_test_cdofb_quadatures_avg(FILE                   *out,
                           const cs_cell_mesh_t   *cm,
                           const int               dim,
                           const _function_type_t  ftype)
{
  const cs_real_t  teval = time_step->t_cur;
  const short int  nf = cm->n_fc;
  const short int  totdof = dim*(nf + 1);

  cs_real_t  *st0, *st1, *st2;
  BFT_MALLOC(st0, totdof, cs_real_t);
  BFT_MALLOC(st1, totdof, cs_real_t);
  BFT_MALLOC(st2, totdof, cs_real_t);

  cs_analytic_func_t  *_func =  _get_func_to_eval(dim, ftype);
  if (_func == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.\n", __func__);

  cs_xdef_analytic_input_t  anai = {.func = _func, .input = NULL};

  /* Reset values */
  memset(st0, 0, totdof*sizeof(cs_real_t));
  memset(st1, 0, totdof*sizeof(cs_real_t));
  memset(st2, 0, totdof*sizeof(cs_real_t));

  cs_xdef_cw_eval_int_t  *compute = NULL;
  switch (dim) {

  case 1:
    compute = cs_xdef_cw_eval_scal_avg_reduction_by_analytic;
    break;
  case 3:
    compute = cs_xdef_cw_eval_vect_avg_reduction_by_analytic;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid dimension.\n",
              __func__);

  } /* Switch */

  compute(cm, teval, (void*)(&anai), CS_QUADRATURE_BARY, st0);
  compute(cm, teval, (void*)(&anai), CS_QUADRATURE_HIGHER, st1);
  compute(cm, teval, (void*)(&anai), CS_QUADRATURE_HIGHEST, st2);

  /* Dump performance and evaluations */
  _dump_quad_res(out, "RED AVG", cm, dim, ftype,
                 st0, st1, st2);

  BFT_FREE(st0);
  BFT_FREE(st1);
  BFT_FREE(st2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test quadrature functions on scalar, vector and tensor valued
 *          analytic functions
 *
 * \param[in]    out    output file
 * \param[in]    cm     pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_main_quadratures(FILE             *out,
                  cs_cell_mesh_t   *cm)
{
  const int  n_runs = 1000;

  fprintf(out, "\n **********************\n");
  fprintf(out, " ** WORKING ON ");
  if (cm->type == FVM_CELL_HEXA) fprintf(out, "HEXA  **\n");
  else                           fprintf(out, "TETRA **\n");
  fprintf(out, " **********************\n");

  /* Check computation of source terms */
  {
    /* Scalar-valued functions for CDO-Fb schemes */
    _test_cdofb_source(out, cm, 1, CONSTANT, n_runs);
    _test_cdofb_source(out, cm, 1, LINEAR, n_runs);
    _test_cdofb_source(out, cm, 1, QUADRATIC, n_runs);
    _test_cdofb_source(out, cm, 1, NON_POLYNOMIAL, n_runs);

    /* Vector-valued functions for CDO-Fb schemes */
    _test_cdofb_source(out, cm, 3, CONSTANT, n_runs);
    _test_cdofb_source(out, cm, 3, LINEAR, n_runs);
    _test_cdofb_source(out, cm, 3, QUADRATIC, n_runs);
    _test_cdofb_source(out, cm, 3, NON_POLYNOMIAL, n_runs);

    /* Scalar-valued functions for CDO-Vb schemes */
    _test_cdovb_source(out, cm, 1, CONSTANT, n_runs);
    _test_cdovb_source(out, cm, 1, LINEAR, n_runs);
    _test_cdovb_source(out, cm, 1, QUADRATIC, n_runs);
    _test_cdovb_source(out, cm, 1, NON_POLYNOMIAL, n_runs);
  }

  /* Check the evaluation of integrals */
  {
    /* Scalar-valued functions */
    _test_quadratures_misc(out, cm, 1, CONSTANT, n_runs);
    _test_quadratures_misc(out, cm, 1, LINEAR, n_runs);
    _test_quadratures_misc(out, cm, 1, QUADRATIC, n_runs);
    _test_quadratures_misc(out, cm, 1, NON_POLYNOMIAL, n_runs);

    /* Vector-valued functions */
    _test_quadratures_misc(out, cm, 3, CONSTANT, n_runs);
    _test_quadratures_misc(out, cm, 3, LINEAR, n_runs);
    _test_quadratures_misc(out, cm, 3, QUADRATIC, n_runs);
    _test_quadratures_misc(out, cm, 3, NON_POLYNOMIAL, n_runs);

    /* Tensor-valued functions */
    _test_quadratures_misc(out, cm, 9, CONSTANT, n_runs);
    _test_quadratures_misc(out, cm, 9, LINEAR, n_runs);
    _test_quadratures_misc(out, cm, 9, QUADRATIC, n_runs);
    _test_quadratures_misc(out, cm, 9, NON_POLYNOMIAL, n_runs);
  }

  /* Check the evaluation of integrals with cs_xdef_t structure */
  {
    /* Scalar-valued functions */
    _test_quadratures_xdef(out, cm, 1, CONSTANT, n_runs);
    _test_quadratures_xdef(out, cm, 1, LINEAR, n_runs);
    _test_quadratures_xdef(out, cm, 1, QUADRATIC, n_runs);
    _test_quadratures_xdef(out, cm, 1, NON_POLYNOMIAL, n_runs);

    /* Vector-valued functions */
    _test_quadratures_xdef(out, cm, 3, CONSTANT, n_runs);
    _test_quadratures_xdef(out, cm, 3, LINEAR, n_runs);
    _test_quadratures_xdef(out, cm, 3, QUADRATIC, n_runs);
    _test_quadratures_xdef(out, cm, 3, NON_POLYNOMIAL, n_runs);

    /* Tensor-valued functions */
    _test_quadratures_xdef(out, cm, 9, CONSTANT, n_runs);
    _test_quadratures_xdef(out, cm, 9, LINEAR, n_runs);
    _test_quadratures_xdef(out, cm, 9, QUADRATIC, n_runs);
    _test_quadratures_xdef(out, cm, 9, NON_POLYNOMIAL, n_runs);

  }

  /* Test average reduction */
  {
    /* Scalar-valued functions */
    _test_cdofb_quadatures_avg(out, cm, 1, CONSTANT);
    _test_cdofb_quadatures_avg(out, cm, 1, LINEAR);
    _test_cdofb_quadatures_avg(out, cm, 1, QUADRATIC);
    _test_cdofb_quadatures_avg(out, cm, 1, NON_POLYNOMIAL);


    /* Vector-valued functions */
    _test_cdofb_quadatures_avg(out, cm, 3, CONSTANT);
    _test_cdofb_quadatures_avg(out, cm, 3, LINEAR);
    _test_cdofb_quadatures_avg(out, cm, 3, QUADRATIC);
    _test_cdofb_quadatures_avg(out, cm, 3, NON_POLYNOMIAL);
  }

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Main program to check CDO/HHO algorithms
 *
 * \param[in]    argc
 * \param[in]    argv
 */
/*----------------------------------------------------------------------------*/

int
main(int    argc,
     char  *argv[])
{
  CS_UNUSED(argc);
  CS_UNUSED(argv);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  {
    int t_id;
#pragma omp parallel private(t_id)
    {
      t_id = omp_get_thread_num();
      if (t_id == 0)
        cs_glob_n_threads = omp_get_max_threads();
    }
  }
#endif

  cs_quadrature_setup();

  quadrature = fopen("Quadrature_tests.log", "w");

  /* =========================== */
  /* TEST DISCRETIZATION SCHEMES */
  /* =========================== */

  /* quant, connect and time_step are declared as static */

  /* connectivity */
  BFT_MALLOC(connect, 1, cs_cdo_connect_t);
  connect->n_max_vbyc = 8;
  connect->n_max_ebyc = 12;
  connect->n_max_fbyc = 6;
  connect->n_max_vbyf = 4;
  connect->v_max_cell_range = 8;
  connect->e_max_cell_range = 12;

  /* Nothing to do for quant (work cellwise) */

  /* Time step */
  BFT_MALLOC(time_step, 1, cs_time_step_t);
  time_step->t_cur = 0.; // Useful when analytic function are called

  cs_source_term_set_shared_pointers(quant, connect);

  /* Allocate local structures */
  cs_cell_mesh_t  *cm = cs_cell_mesh_create(connect);

  cs_cdofb_scaleq_init_common(quant, connect, time_step, NULL);
  cs_cdofb_vecteq_init_common(quant, connect, time_step, NULL);
  cs_cdovb_scaleq_init_common(quant, connect, time_step, NULL);
  cs_cdovb_vecteq_init_common(quant, connect, time_step, NULL);

  /* ========= */
  /* TEST HEXA */
  /* ========= */

  _define_cm_hexa_unif(1., cm);

  /* Test quadrature rules */
  _main_quadratures(quadrature, cm);

  /* ========== */
  /* TEST TETRA */
  /* ========== */

  _define_cm_tetra_ref(1., cm);

  /* Test quadrature rules */
  _main_quadratures(quadrature, cm);

  /* Free memory */
  cs_cdofb_scaleq_finalize_common();
  cs_cdofb_vecteq_finalize_common();
  cs_cdovb_scaleq_finalize_common();
  cs_cdovb_vecteq_finalize_common();

  cs_cell_mesh_free(&cm);

  BFT_FREE(connect);
  BFT_FREE(time_step);

  fclose(quadrature);

  printf(" --> Quadrature Tests (Done)\n");
  exit (EXIT_SUCCESS);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
