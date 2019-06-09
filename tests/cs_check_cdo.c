/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_basis_func.h"
#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_local.h"
#include "cs_cdo_quantities.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdovb_scaleq.h"
#include "cs_equation_param.h"
#include "cs_evaluate.h"
#include "cs_hho_builder.h"
#include "cs_hho_scaleq.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_param.h"
#include "cs_param_cdo.h"
#include "cs_scheme_geometry.h"
#include "cs_sdm.h"
#include "cs_source_term.h"
#include "cs_time_step.h"
#include "cs_timer.h"
#include "cs_xdef_cw_eval.h"
#include "cs_xdef_eval.h"


/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Static global variables
 *============================================================================*/

static FILE  *hexa = NULL;
static FILE  *tetra = NULL;
static FILE  *hexa_hho0 = NULL;
static FILE  *hexa_hho1 = NULL;
static FILE  *hexa_hho2 = NULL;
static FILE  *tetra_hho0 = NULL;
static FILE  *tetra_hho1 = NULL;
static FILE  *tetra_hho2 = NULL;

static cs_cdo_connect_t  *connect = NULL;
static cs_cdo_quantities_t  *quant = NULL;
static cs_time_step_t  *time_step = NULL;

static const double non_poly_int[] = { 1.1319870772271508,  /* Hexa */
                                       0.0801351697078868}; /* Tetra */
static double non_poly_int_f[2][6] = {
  {0.6587901114231911, 0.6587901114231911, 1.7907771886504419,
   1.7907771886504419, 0.6587901114231911, 1.7907771886504419},
  {0.2231301601484198, 0.2231301601484198, 0.2231301601484198,
   0.5252709594852753, 0., 0.}
};

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
 * Fill array with x*x + y*y + z*z
 */
static void
_quadratic_xyz2(cs_real_t          time,
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
    retval[r] = cs_math_3_square_norm(xyz + 3*p);
  }

}


/* Test function based on analytic definition
 * Fill array with [x*x*x, y*y*y, z*z*z]
 */
static void
_cubic_xyz3_vect(cs_real_t          time,
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
    for (short int j = 0; j < 3; j++)
      retval[r+j] = cs_math_pow3(xyz[p+j]);
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the orthogonal projection of the point xm on the face f
 *
 * \param[in]      pfq        quantities related to the primal face f
 * \param[in]      xm         coordinates of the point M to project
 * \param[in, out] xp         coordinates of the projection
 */
/*----------------------------------------------------------------------------*/

static inline void
_ortho_proj(const cs_quant_t      pfq,
            const cs_real_3_t     xm,
            cs_real_3_t           xp)
{
  cs_real_3_t  mf = {xm[0] - pfq.center[0],
                     xm[1] - pfq.center[1],
                     xm[2] - pfq.center[2]};

  const cs_real_t  dp = _dp3(mf, pfq.unitv);
  for (int k = 0; k < 3; k++)
    xp[k] = mf[k] - dp * pfq.unitv[k] + pfq.center[k];
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
  cm->flag = CS_FLAG_COMP_PV |CS_FLAG_COMP_PVQ | CS_FLAG_COMP_PEQ |
    CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ |
    CS_FLAG_COMP_DFQ | CS_FLAG_COMP_HFQ | CS_FLAG_COMP_FE | CS_FLAG_COMP_EFQ |
    CS_FLAG_COMP_DIAM;
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
  cm->flag = CS_FLAG_COMP_PV |CS_FLAG_COMP_PVQ | CS_FLAG_COMP_PEQ |
    CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ |
    CS_FLAG_COMP_DFQ | CS_FLAG_COMP_HFQ | CS_FLAG_COMP_FE |CS_FLAG_COMP_EFQ  |
    CS_FLAG_COMP_DIAM;

  cm->vol_c = cs_math_1ov6*a*a*a;
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

  // Dual faces, wvc ?

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

  } // Loop on cell faces

  /* Compute dual face quantities */
  for (short int e = 0; e < cm->n_ec; e++) {

    cs_real_3_t  df;
    const cs_nvec3_t  s1 = cm->sefc[2*e], s2 = cm->sefc[2*e+1];
    for (int k = 0; k < 3; k++)
      df[k] = s1.meas*s1.unitv[k] + s2.meas*s2.unitv[k];
    cs_nvec3(df, &(cm->dface[e]));

  } // Loop on cell edges

  /* Compute dual cell volume */
  for (short int f = 0; f < cm->n_fc; f++) {

    const double  hf_coef = cs_math_1ov6 * cm->hfc[f];

    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];
      const double  half_pef_vol = cm->tef[i]*hf_coef;

      cm->wvc[v1] += half_pef_vol;
      cm->wvc[v2] += half_pef_vol;

    } // Loop on face edges

  } // Loop on cell faces

  /* Reset diam */
  double  dbuf[10];
  short int  vtag[4];
  int  size = cm->n_vc*(cm->n_vc+1)/2;
  int  shift = 0;

  cm->diam_c = -1;
  for (int i = 0; i < size; i++) dbuf[i] = 0.;

  for (short int vi = 0; vi < cm->n_vc; ++vi) {
    shift++; // diag entry not taken into account
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
 * \brief   Dump a local discrete Hodge operator
 *
 * \param[in]    fic        pointer to a FILE structure
 * \param[in]    msg        optional message to print
 * \param[in]    dof_ids    id related to each Degree of Freedom
 * \param[in]    lm         pointer to the cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_locmat_dump(FILE               *fic,
             const char         *msg,
             const int          *dof_ids,
             const cs_sdm_t     *lm)
{
  assert(fic != NULL && lm != NULL);

  if (msg != NULL)
    fprintf(fic, "%s\n", msg);

  /* List sub-entity ids */
  fprintf(fic, "%6s","ID");
  for (int i = 0; i < lm->n_rows; i++) fprintf(fic, " %11d", dof_ids[i]);
  fprintf(fic, "\n");

  for (int i = 0; i < lm->n_rows; i++) {
    fprintf(fic, " %5d", dof_ids[i]);
    for (int j = 0; j < lm->n_rows; j++)
      fprintf(fic, " % 6.4e", lm->val[i*lm->n_rows+j]);
    fprintf(fic, "\n");
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a local discrete Hodge operator
 *
 * \param[in]    fic        pointer to a FILE structure
 * \param[in]    msg        optional message to print
 * \param[in]    csys       pointer to the cs_cell_sys_t  struct.
 */
/*----------------------------------------------------------------------------*/

static void
_locsys_dump(FILE                 *fic,
             const char           *msg,
             const cs_cell_sys_t  *csys)
{
  assert(fic != NULL && csys != NULL);
  const cs_sdm_t  *lm = csys->mat;

  if (msg != NULL)
    fprintf(fic, "%s\n", msg);

  /* List sub-entity ids */
  fprintf(fic, "%6s","ID");
  for (int i = 0; i < lm->n_rows; i++) fprintf(fic, " %11d", csys->dof_ids[i]);
  fprintf(fic, "%11s %11s %11s\n", "RHS", "SOURCE", "VAL_N");
  for (int i = 0; i < lm->n_rows; i++) {
    fprintf(fic, " %5d", csys->dof_ids[i]);
    for (int j = 0; j < lm->n_rows; j++)
      fprintf(fic, " % 6.4e", lm->val[i*lm->n_rows+j]);
    fprintf(fic, " % 6.4e % 6.4e % 6.4e\n",
            csys->rhs[i], csys->source[i], csys->val_n[i]);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Analyse a Hodge operator for Face-based schemes
 *
 * \param[in]    out     output file
 * \param[in]    cm      pointer to a cs_cell_mesh_t structure
 * \param[in]    hdg     pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_test_hodge_fb(FILE               *out,
               cs_cell_mesh_t     *cm,
               cs_sdm_t           *hdg)
{
  fprintf(out, "\n");

  for (short int fi = 0; fi < cm->n_fc; fi++) {
    double  row_sum = 0.;
    double  *row_vals = hdg->val + fi*(cm->n_fc + 1);
    for (short int fj = 0; fj < cm->n_fc + 1; fj++) row_sum += row_vals[fj];
    fprintf(out, "F%d = % 9.6e\n", fi, row_sum);
  }

  double  crow_sum = 0.;
  double  *crow_vals = hdg->val + cm->n_fc*(cm->n_fc + 1);
  for (short int fj = 0; fj < cm->n_fc + 1; fj++) crow_sum += crow_vals[fj];
  fprintf(out, "Cell = % 9.6e\n", crow_sum);

  assert(cm->n_fc + 1 <= 7);
  cs_real_t  a[7], res[7];
  for (int i = 0; i < 7; i++) a[i] = 1.;
  cs_sdm_matvec(hdg, a, res);
  for (int i = 0; i < cm->n_fc + 1; i++)
    fprintf(out, "a[%d] = % -8.5e --> res[%d] = % -8.5e\n",
            i, a[i], i, res[i]);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Analyse a Hodge operator for Vertex-based schemes
 *
 * \param[in]    out     output file
 * \param[in]    cm      pointer to a cs_cell_mesh_t structure
 * \param[in]    hdg     pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_test_hodge_vb(FILE               *out,
               cs_cell_mesh_t     *cm,
               cs_sdm_t           *hdg)
{
  fprintf(out, "\n");

  for (short int vi = 0; vi < cm->n_vc; vi++) {
    double  row_sum = 0.;
    double  *row_vals = hdg->val + vi*cm->n_vc;
    for (short int vj = 0; vj < cm->n_vc; vj++) row_sum += row_vals[vj];
    fprintf(out, "V%d = % 9.6e |delta=% 9.6e\n",
            vi, row_sum, row_sum - cm->wvc[vi]*cm->vol_c);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Analyse a stiffness matrix for Vertex-based schemes
 *
 * \param[in]    out    output file
 * \param[in]    cm     pointer to a cs_cell_mesh_t structure
 * \param[in]    s      pointer to a cs_sdm_t structure (stiffness matrix)
 */
/*----------------------------------------------------------------------------*/

static void
_test_stiffness_vb(FILE               *out,
                   cs_cell_mesh_t     *cm,
                   cs_sdm_t           *s)
{
  fprintf(out, "\nCDO.VB;   %10s %10s\n", "ROW_SUM", "LIN_SUM");
  for (short int vi = 0; vi < cm->n_vc; vi++) {
    double  row_sum = 0., linear_sum = 0.;
    double  *row_vals = s->val + vi*cm->n_vc;
    for (short int vj = 0; vj < cm->n_vc; vj++) {
      row_sum += row_vals[vj];
      linear_sum += row_vals[vj]*cm->xv[3*vj];
    }
    fprintf(out, "  V%d = % 9.6e % 9.6e\n", vi, row_sum, linear_sum);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test CDO vertex-based schemes
 */
/*----------------------------------------------------------------------------*/

static void
_test_cdovb_schemes(FILE                *out,
                    cs_cell_mesh_t      *cm,
                    cs_face_mesh_t      *fm,
                    cs_cell_sys_t       *csys,
                    cs_cell_builder_t   *cb)
{
  /* Initialize a cell view of the algebraic system */
  csys->n_dofs = cm->n_vc;
  cs_sdm_square_init(csys->n_dofs, csys->mat);
  for (short int v = 0; v < cm->n_vc; v++)
    csys->dof_ids[v] = cm->v_ids[v];

  /* Handle anisotropic diffusion */
  cb->dpty_mat[0][0] = 1.0, cb->dpty_mat[0][1] = 0.5, cb->dpty_mat[0][2] = 0.0;
  cb->dpty_mat[1][0] = 0.5, cb->dpty_mat[1][1] = 1.0, cb->dpty_mat[1][2] = 0.5;
  cb->dpty_mat[2][0] = 0.0, cb->dpty_mat[2][1] = 0.5, cb->dpty_mat[2][2] = 1.0;

  /* Useful for a weak enforcement of the BC */
  cs_math_33_eigen((const cs_real_t (*)[3])cb->dpty_mat,
                   &(cb->eig_ratio),
                   &(cb->eig_max));

  /* HODGE */
  /* ===== */

  /* WBS Hodge operator */

  cs_param_hodge_t  hwbs_info = {.is_unity = true,
                                 .is_iso = true,
                                 .inv_pty = false,
                                 .type = CS_PARAM_HODGE_TYPE_VPCD,
                                 .algo = CS_PARAM_HODGE_ALGO_WBS,
                                 .coef = 1.0};

  cs_hodge_vpcd_wbs_get(hwbs_info, cm, cb);
  _locmat_dump(out, "\nCDO.VB; HDG.VPCD.WBS; PERMEABILITY.ISO",
               csys->dof_ids, cb->hdg);
  _test_hodge_vb(out, cm, cb->hdg);

  cs_hodge_compute_wbs_surfacic(fm, cb->hdg);
  _locmat_dump(out, "\nCDO.VB; HDG.VPCD.WBS.FACE; UNITY",
               csys->dof_ids, cb->hdg);

  for (int vi = 0; vi < fm->n_vf; vi++) {
    double  row_sum = 0.0;
    double  *hi = cb->hdg->val + vi*fm->n_vf;
    for (int vj = 0; vj < fm->n_vf; vj++) row_sum += hi[vj];
    fprintf(out, "V%d = %6.4e |delta= %6.4e\n",
            vi, row_sum, row_sum - fm->face.meas/fm->n_vf);
  }

  /* Voronoi Hodge operator */
  cs_param_hodge_t  hvor_info = {.is_unity = true,
                                 .is_iso = true,
                                 .inv_pty = false,
                                 .type = CS_PARAM_HODGE_TYPE_VPCD,
                                 .algo = CS_PARAM_HODGE_ALGO_VORONOI,
                                 .coef = 1.0};
  cs_hodge_vpcd_voro_get(hvor_info, cm, cb);
  _locmat_dump(out, "\nCDO.VB; HDG.VPCD.VORONOI; PERMEABILITY.ISO",
               csys->dof_ids, cb->hdg);
  _test_hodge_vb(out, cm, cb->hdg);

  /* DIFFUSION (Stiffness matrix) */
  /* ============================ */

  cs_equation_param_t  *eqp = cs_equation_create_param("Test.CDOVB",
                                                       CS_EQUATION_TYPE_USER,
                                                       1,
                                                       CS_PARAM_BC_HMG_NEUMANN);

  eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;

  /* Stiffness matrix arising from a Hodge EpFd built with COST algo. */
  eqp->diffusion_hodge.is_unity = true;
  eqp->diffusion_hodge.is_iso = true;
  eqp->diffusion_hodge.inv_pty = false; /* inverse property ? */
  eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
  eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
  eqp->diffusion_hodge.coef = 1./3.; /* DGA algo. */

  { /* Compute the stifness matrix (CO+ST -- DGA) */
    cs_hodge_vb_cost_get_stiffness(eqp->diffusion_hodge, cm, cb);
    _locmat_dump(out,"\nCDO.VB; STIFFNESS WITH HDG.EPFD.DGA; PERMEABILITY.ISO",
                 csys->dof_ids, cb->loc);
    _test_stiffness_vb(out, cm, cb->loc);

    /* Anisotropic case */
    eqp->diffusion_hodge.is_unity = false;
    eqp->diffusion_hodge.is_iso = false;

    cs_hodge_vb_cost_get_stiffness(eqp->diffusion_hodge, cm, cb);
    _locmat_dump(out,
                 "\nCDO.VB; STIFFNESS WITH HDG.EPFD.DGA; PERMEABILITY.ANISO",
                 csys->dof_ids, cb->loc);
    _test_stiffness_vb(out, cm, cb->loc);

    /* Enforce Dirichlet BC */
    cs_equation_set_param(eqp, CS_EQKEY_BC_ENFORCEMENT, "penalization");
    cs_cdo_diffusion_pena_dirichlet(eqp, cm, fm, cb, csys);
    _locsys_dump(out, "\nCDO.VB; PENA.DGA.FLX.COST; PERMEABILITY.ANISO",
                 csys);
    for (int v = 0; v < cm->n_vc; v++) csys->rhs[v] = 0;
    for (int v = 0; v < cm->n_vc*cm->n_vc; v++) csys->mat->val[v] = 0.;

    /* Enforce Dirichlet BC */
    cs_equation_set_param(eqp, CS_EQKEY_BC_ENFORCEMENT, "weak");
    cs_cdo_diffusion_svb_cost_weak_dirichlet(eqp, cm, fm, cb, csys);
    _locsys_dump(out, "\nCDO.VB; WEAK.DGA.FLX.COST; PERMEABILITY.ANISO",
                 csys);

    for (int v = 0; v < cm->n_vc; v++) csys->rhs[v] = 0;
    for (int v = 0; v < cm->n_vc*cm->n_vc; v++) csys->mat->val[v] = 0.;

    /* Enforce Dirichlet BC */
    cs_equation_set_param(eqp, CS_EQKEY_BC_ENFORCEMENT, "weak_sym");
    cs_cdo_diffusion_svb_cost_wsym_dirichlet(eqp, cm, fm, cb, csys);
    _locsys_dump(out, "\nCDO.VB; WSYM.DGA.FLX.COST; PERMEABILITY.ANISO",
                 csys);
    for (int v = 0; v < cm->n_vc; v++) csys->rhs[v] = 0;
    cs_sdm_square_init(cm->n_vc, csys->mat);
  }

  /* ------------------------------------------------------------------------ */

  { /* Compute the stifness matrix (VORONOI) */

    eqp->diffusion_hodge.is_unity = true;
    eqp->diffusion_hodge.is_iso = true;
    eqp->diffusion_hodge.inv_pty = false; /* inverse property ? */
    eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
    eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;

    cs_hodge_vb_voro_get_stiffness(eqp->diffusion_hodge, cm, cb);
    _locmat_dump(out,
                 "\nCDO.VB; STIFFNESS WITH HDG.EPFD.VORO; PERMEABILITY.ISO",
                 csys->dof_ids, cb->loc);
    _test_stiffness_vb(out, cm, cb->loc);

    /* Enforce Dirichlet BC */
    cs_equation_set_param(eqp, CS_EQKEY_BC_ENFORCEMENT, "penalization");
    cs_cdo_diffusion_pena_dirichlet(eqp, cm, fm, cb, csys);
    _locsys_dump(out,
                 "\nCDO.VB; PENA.VORO.FLX.COST; PERMEABILITY.ISO",
                 csys);
    for (int v = 0; v < cm->n_vc; v++) csys->rhs[v] = 0;
    cs_sdm_square_init(cm->n_vc, csys->mat);
  }

  /* ------------------------------------------------------------------------ */

  { /* Stiffness matrix arising from a Hodge EpFd built with WBS algo. */

    eqp->diffusion_hodge.is_unity = true;
    eqp->diffusion_hodge.is_iso = true;
    eqp->diffusion_hodge.inv_pty = false; /* inverse property ? */
    eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
    eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;

    cs_hodge_vb_wbs_get_stiffness(eqp->diffusion_hodge, cm, cb);
    _locmat_dump(out,
                 "\nCDO.VB; STIFFNESS WITH HDG.EPFD.WBS; PERMEABILITY.ISO",
                 csys->dof_ids, cb->loc);
    _test_stiffness_vb(out, cm, cb->loc);

    eqp->diffusion_hodge.is_unity = false;
    eqp->diffusion_hodge.is_iso = false;
    cs_hodge_vb_wbs_get_stiffness(eqp->diffusion_hodge, cm, cb);
    _locmat_dump(out,
                 "\nCDO.VB; STIFFNESS WITH HDG.EPFD.WBS; PERMEABILITY.ANISO",
                 csys->dof_ids, cb->loc);
    _test_stiffness_vb(out, cm, cb->loc);

    /* Enforce Dirichlet BC */
    cs_equation_set_param(eqp, CS_EQKEY_BC_ENFORCEMENT, "penalization");
    cs_cdo_diffusion_pena_dirichlet(eqp, cm, fm, cb, csys);
    _locsys_dump(out, "\nCDO.VB; PENA.WBS.FLX.WBS; PERMEABILITY.ANISO",
                 csys);
    for (int v = 0; v < cm->n_vc; v++) csys->rhs[v] = 0;
    cs_sdm_square_init(cm->n_vc, csys->mat);

    /* Enforce Dirichlet BC */
    cs_equation_set_param(eqp, CS_EQKEY_BC_ENFORCEMENT, "weak");
    cs_cdo_diffusion_svb_wbs_weak_dirichlet(eqp, cm, fm, cb, csys);
    _locsys_dump(out, "\nCDO.VB; WEAK.WBS.FLX.WBS; PERMEABILITY.ANISO", csys);
    for (int v = 0; v < cm->n_vc; v++) csys->rhs[v] = 0;
    cs_sdm_square_init(cm->n_vc, csys->mat);

        /* Enforce Dirichlet BC */
    cs_equation_set_param(eqp, CS_EQKEY_BC_ENFORCEMENT, "weak_sym");
    cs_cdo_diffusion_svb_wbs_wsym_dirichlet(eqp, cm,fm, cb, csys);
    _locsys_dump(out, "\nCDO.VB; WSYM.WBS.FLX.WBS; PERMEABILITY.ANISO", csys);
    for (int v = 0; v < cm->n_vc; v++) csys->rhs[v] = 0;
    cs_sdm_square_init(cm->n_vc, csys->mat);
  }

  /* ADVECTION OPERATOR */
  /* ================== */

  cs_adv_field_t  *beta = cs_advection_field_add_user("Adv.Field");

  /* Numerical settings for the advection scheme */
  eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_CONSERV;
  eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;

  /* Constant advection field */
  cs_real_3_t  vector_field = {1., 0., 0.};
  cs_advection_field_def_by_value(beta, vector_field);
  eqp->adv_field = beta;

  /* Free memory */
  cs_advection_field_destroy_all();
  eqp = cs_equation_free_param(eqp);

  /* ADVECTION: BOUNDARY FLUX OPERATOR */
  /* ================================= */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print the comparison of the evaluation of the basis functions at
 *          different mesh locations
 *
 * \param[in]  out          output file
 * \param[in]  n_cases      number of evaluations * number of basis functions
 * \param[in]  labels       list of labels (n_cases)
 * \param[in]  eval_i       results of the evaluation with inertia basis
 * \param[in]  eval_m       results of the evaluation with monomial bais
 */
/*----------------------------------------------------------------------------*/

static inline void
_dump_eval_cmp(FILE             *out,
               int               n_cases,
               const char       *labels[],
               const cs_real_t   eval_i[],
               const cs_real_t   eval_m[])
{
  fprintf(out, " %8s | %12s | %12s |\n", "Cases", "Inertia", "Monomial");
  for (int i = 0; i < n_cases; i++)
    fprintf(out, " %8s | % .5e | % .5e |\n", labels[i], eval_i[i], eval_m[i]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test basis functions implementation
 *
 * \param[in]  out           output file
 * \param[in]  scheme_order  scheme_order
 * \param[in]  cm            pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_test_basis_functions(FILE               *out,
                      int                 scheme_order,
                      cs_cell_mesh_t     *cm,
                      cs_cell_builder_t  *cb)
{
  switch (scheme_order) {

  case 0:
    fprintf(out,
            "\n ************************************************************\n"
            "          TEST BASIS FUNCTION Oth order\n"
            " ************************************************************\n");
    break;

  case 1:
    fprintf(out,
            "\n ************************************************************\n"
            "          TEST BASIS FUNCTION 1st order\n"
            " ************************************************************\n");

    break;

  case 2:
    fprintf(out,
            "\n ************************************************************\n"
            "           TEST BASIS FUNCTION 2nd order\n"
            " ************************************************************\n");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Scheme order not handled yet.\n",
              __func__);

  }

  /* Test cell basis functions : (flag, order, dim) */
  cs_real_t  c_eval[30], c_eval_mono[30], c_dof[20], c_dof_mono[20];
  cs_real_t  g_eval[180], g_eval_mono[180];
  for (int i = 0; i < 180; i++) g_eval[i] = g_eval_mono[i] = -999.999;

  /* Inertial basis */
  cs_basis_func_t  *cbf = cs_basis_func_create(0, scheme_order, 3);
  cbf->setup(cbf, cm, 0, cm->xc, cb);
  cbf->eval_all_at_point(cbf, cm->xc, c_eval);
  cbf->eval_all_at_point(cbf, cm->xv, c_eval + cbf->size);
  cbf->eval_all_at_point(cbf, cm->xv + 3, c_eval + 2*cbf->size);
  cbf->compute_projector(cbf, cm, 0);
  cbf->compute_factorization(cbf);

  /* Case of a unit constant */
  cbf->project(cbf, cbf->projector->val, c_dof);
  if (scheme_order > 0) /* Case of a (x -xc)/diam function (1st basis func.) */
    cbf->project(cbf, cbf->projector->val + cbf->size, c_dof + cbf->size);

  /* Output inertial basis functions */
  cbf->dump_projector(cbf);
  fprintf(out, "\n >> Inertial cell basis functions\n");
  cs_basis_func_fprintf(out, NULL, cbf);
  cs_sdm_fprintf(out, NULL, 1e-15, cbf->projector);

  /* Define the related gradient basis */
  cs_basis_func_t  *gbf = cs_basis_func_grad_create(cbf);
  cs_basis_func_copy_setup(cbf, gbf);
  const int bsize = 3*gbf->size;

  fprintf(out, "\n >> Inertial cell gradient basis functions\n");
  cs_basis_func_fprintf(out, NULL, gbf);

  gbf->eval_all_at_point(gbf, cm->xc, g_eval);
  gbf->eval_all_at_point(gbf, cm->xv, g_eval + bsize);
  gbf->eval_all_at_point(gbf, cm->xv + 3, g_eval + 2*bsize);

  /* Free basis structures */
  gbf = cs_basis_func_free(gbf);
  cbf = cs_basis_func_free(cbf);

  /* Monomial basis */
  cs_basis_func_t  *cbf_mono = cs_basis_func_create(CS_BASIS_FUNC_MONOMIAL,
                                                    scheme_order, 3);
  cbf_mono->setup(cbf_mono, cm, 0, cm->xc, cb);
  cs_basis_func_dump(cbf_mono);
  cbf_mono->eval_all_at_point(cbf_mono, cm->xc, c_eval_mono);
  cbf_mono->eval_all_at_point(cbf_mono, cm->xv,
                              c_eval_mono + cbf_mono->size);
  cbf_mono->eval_all_at_point(cbf_mono, cm->xv + 3,
                              c_eval_mono + 2*cbf_mono->size);
  cbf_mono->compute_projector(cbf_mono, cm, 0);
  cbf_mono->compute_factorization(cbf_mono);

  /* Case of a unit constant */
  cbf_mono->project(cbf_mono, cbf_mono->projector->val, c_dof_mono);
  if (scheme_order > 0) /* Case of a (x -xc)/diam function (1st basis func.) */
    cbf_mono->project(cbf_mono, cbf_mono->projector->val + cbf_mono->size,
                      c_dof_mono + cbf_mono->size);

  /* Output inertial basis functions */
  cbf_mono->dump_projector(cbf_mono);
  fprintf(out, "\n >> Monomial cell basis functions\n");
  cs_basis_func_fprintf(out, NULL, cbf_mono);
  cs_sdm_fprintf(out, NULL, 1e-15, cbf_mono->projector);

  /* Define the related gradient basis */
  cs_basis_func_t  *gbf_mono = cs_basis_func_grad_create(cbf_mono);
  cs_basis_func_copy_setup(cbf_mono, gbf_mono);

  cs_basis_func_dump(gbf_mono);
  gbf_mono->eval_all_at_point(gbf_mono, cm->xc, g_eval_mono);
  gbf_mono->eval_all_at_point(gbf_mono, cm->xv, g_eval_mono + bsize);
  gbf_mono->eval_all_at_point(gbf_mono, cm->xv + 3, g_eval_mono + 2*bsize);

  fprintf(out, "\n >> Mononial cell gradient basis functions\n");
  cs_basis_func_fprintf(out, NULL, gbf_mono);

  /* Free basis structures */
  gbf_mono = cs_basis_func_free(gbf_mono);
  cbf_mono = cs_basis_func_free(cbf_mono);

  /* Comparison between the two choices of building a basis function */
  printf( "\n Evaluation points for cell:\n xc = (%5.3e, %5.3e, %5.3e)"
          " xv1 = (%5.3e, %5.3e, %5.3e) xv2 = (%5.3e, %5.3e, %5.3e)\n",
          cm->xc[0], cm->xc[1], cm->xc[2], cm->xv[0], cm->xv[1], cm->xv[2],
          cm->xv[3], cm->xv[4], cm->xv[5]);

  const char  *cvv[3] = {"xc", "xv1", "xv2"}, *cdof[2] = {"unit", "px"};
  const char  *c_tags[10] =
    {"0", "x", "y", "z", "xx", "xy", "xz", "yy", "yz", "zz" };
  const char  *g_tags[60] = {"0:x",    "0:y",    "0:z",
                             "gx:x",   "gx:y",   "gx:z",
                             "gy:x",   "gy:y",   "gy:z",
                             "gz:x",   "gz:y",   "gz:z",
                             "gx2:x",  "gx2:y",  "gx2:z",
                             "gxy:x",  "gxy:y",  "gxy:z",
                             "gxz:x",  "gxz:y",  "gxz:z",
                             "gy2:x",  "gy2:y",  "gy2:z",
                             "gyz:x",  "gyz:y",  "gyz:z",
                             "gz2:x",  "gz2:y",  "gz2:z",
                             "gx3:x",  "gx3:y",  "gx3:z",
                             "gx2y:x", "gx2y:y", "gx2y:z",
                             "gx2z:x", "gx2z:y", "gx2z:z",
                             "gxy2:x", "gxy2:y", "gxy2:z",
                             "gxyz:x", "gxyz:y", "gxyz:z",
                             "gxz2:x", "gxz2:y", "gxz2:z",
                             "gy3:x",  "gy3:y",  "gy3:z",
                             "gy2z:x", "gy2z:y", "gy2z:z",
                             "gyz2:x", "gyz2:y", "gyz2:z",
                             "gy3:x",  "gy3:y",  "gy3:z" };

  switch (scheme_order) {

  case 0:
    printf( " --> %s\n", cdof[0]);
    _dump_eval_cmp(stdout, 1, c_tags, c_dof, c_dof_mono);
    for (int i = 0; i < 3; i++) {
      printf( " --> %s\n", cvv[i]);
      _dump_eval_cmp(stdout, 1, c_tags, c_eval + i, c_eval_mono + i);
      _dump_eval_cmp(stdout, 12, g_tags, g_eval + i*bsize, g_eval_mono + i*bsize);
    }
    break;

  case 1:
    for (int i = 0; i < 2; i++) {
      printf( " --> %s\n", cdof[i]);
      _dump_eval_cmp(stdout, 4, c_tags, c_dof + 4*i, c_dof_mono + 4*i);
    }
    for (int i = 0; i < 3; i++) {
      printf( " --> %s\n", cvv[i]);
      _dump_eval_cmp(stdout, 4, c_tags, c_eval + 4*i, c_eval_mono + 4*i);
      _dump_eval_cmp(stdout, 30, g_tags, g_eval + i*bsize, g_eval_mono + i*bsize);
    }
    break;

  case 2:
    for (int i = 0; i < 2; i++) {
      printf( " --> %s\n", cdof[i]);
      _dump_eval_cmp(stdout, 10, c_tags, c_dof + 10*i, c_dof_mono + 10*i);
    }
    for (int i = 0; i < 3; i++) {
      printf( " --> %s\n", cvv[i]);
      _dump_eval_cmp(stdout, 10, c_tags, c_eval + 10*i, c_eval_mono + 10*i);
      _dump_eval_cmp(stdout, 60, g_tags, g_eval + i*bsize, g_eval_mono + i*bsize);
    }
    break;

  default:
    break;
  }

  /* Test face basis functions */
  for (short int f = 0; f < cm->n_fc; f++) {

    short int v = cm->e2v_ids[2*cm->f2e_ids[cm->f2e_idx[f]]];
    cs_real_t  f_eval[12], f_eval_mono[12], f_dof[12], f_dof_mono[12];

    cs_basis_func_t  *fbf = cs_basis_func_create(0, scheme_order, 2);
    fbf->setup(fbf, cm, f, cm->face[f].center, cb);
    cs_basis_func_dump(fbf);
    fbf->eval_all_at_point(fbf, cm->face[f].center, f_eval);
    fbf->eval_all_at_point(fbf, cm->xv + 3*v,
                           f_eval + fbf->size);

    fbf->compute_projector(fbf, cm, f);
    fbf->compute_factorization(fbf);

    /* Case of a unit constant */
    fbf->project(fbf, fbf->projector->val, f_dof);
    if (scheme_order > 0)
      /* Case of a (x - xf)/diam function (1st basis func.) */
      fbf->project(fbf, fbf->projector->val + fbf->size,
                   f_dof + fbf->size);

    fbf = cs_basis_func_free(fbf);

    cs_basis_func_t  *fbf_mono = cs_basis_func_create(CS_BASIS_FUNC_MONOMIAL,
                                                      scheme_order, 2);
    fbf_mono->setup(fbf_mono, cm, f, cm->face[f].center, cb);
    cs_basis_func_dump(fbf_mono);
    fbf_mono->eval_all_at_point(fbf_mono, cm->face[f].center, f_eval_mono);
    fbf_mono->eval_all_at_point(fbf_mono, cm->xv + 3*v,
                                f_eval_mono + fbf_mono->size);

    fbf_mono->compute_projector(fbf_mono, cm, f);
    fbf_mono->compute_factorization(fbf_mono);

    /* Case of a unit constant */
    fbf_mono->project(fbf_mono, fbf_mono->projector->val, f_dof_mono);
    if (scheme_order > 0)
      /* Case of a (x - xf)/diam function (1st basis func.) */
      fbf_mono->project(fbf_mono, fbf_mono->projector->val + fbf_mono->size,
                        f_dof_mono + fbf_mono->size);

    fbf_mono = cs_basis_func_free(fbf_mono);

    printf( "\n Evaluation points for f=%d: xf = (%5.3e, %5.3e, %5.3e)"
            " xv = (%5.3e, %5.3e, %5.3e)\n", f,
            cm->face[f].center[0], cm->face[f].center[1], cm->face[f].center[2],
            cm->xv[3], cm->xv[4], cm->xv[5]);

    /* Comparison between the two choices of building a basis function */
    const char  *fv[2] = {"xf", "xv"};
    const char  *f_tags[6] = {"0", "x", "y", "xx", "xy", "yy"};
    switch (scheme_order) {

    case 0:
      printf( " --> %s (f:%d)\n", cdof[0], f);
      _dump_eval_cmp(stdout, 1, f_tags, f_dof, f_dof_mono);
      for (int i = 0; i < 2; i++) {
        printf( "\n --> %s (f:%d)\n", fv[i], f);
        _dump_eval_cmp(stdout, 1, f_tags, f_eval + i, f_eval_mono + i);
      }
      break;

    case 1:
      for (int i = 0; i < 2; i++) {
        printf( " --> %s (f:%d)\n", cdof[i], f);
        _dump_eval_cmp(stdout, 3, f_tags, f_dof + 3*i, f_dof_mono + 3*i);
        printf( " --> %s (f:%d)\n", fv[i], f);
        _dump_eval_cmp(stdout, 3, f_tags, f_eval + 3*i, f_eval_mono + 3*i);
      }
      break;

    case 2:
      for (int i = 0; i < 2; i++) {
        printf( "\n --> %s (f:%d)\n", cdof[i], f);
        _dump_eval_cmp(stdout, 6, f_tags, f_dof + 6*i, f_dof_mono + 6*i);
        printf( "\n --> %s (f:%d)\n", fv[i], f);
        _dump_eval_cmp(stdout, 6, f_tags, f_eval + 6*i, f_eval_mono + 6*i);
      }
      break;

    default:
      break;
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test HHO (Hybrid High Order) schemes
 *
 * \param[in]      out           output file
 * \param[in]      scheme_order  scheme_order
 * \param[in]      cm            pointer to a cs_cell_mesh_t structure
 * \param[in]      fm            pointer to a cs_face_mesh_t structure
 * \param[in, out] csys          pointer to a cs_cell_sys_t structure
 * \param[in, out] cb            pointer to a cs_cell_bc_t structure
 * \param[in, out] hhob          pointer to a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_test_hho_schemes(FILE                *out,
                  int                  scheme_order,
                  cs_cell_mesh_t      *cm,
                  cs_face_mesh_t      *fm,
                  cs_cell_sys_t       *csys,
                  cs_cell_builder_t   *cb,
                  cs_hho_builder_t    *hhob)
{
  CS_UNUSED(fm);
  CS_UNUSED(csys);
  const double  tcur = 0.;

  switch (scheme_order) {

  case 0:
    fprintf(out,
            "\n ************************************************************\n"
            "                              HHO_P0 scheme\n"
            " ************************************************************\n");
    break;

  case 1:
    fprintf(out,
            "\n ************************************************************\n"
            "                              HHO_P1 scheme\n"
            " ************************************************************\n");

    cs_basis_func_set_hho_flag(0, 0); /* Inertial basis function */
    cs_hho_builder_cellwise_setup(cm, cb, hhob);

    fprintf(out, "\n >> Inertial cell basis functions\n");
    cs_basis_func_fprintf(out, NULL, hhob->cell_basis);
    cs_sdm_fprintf(out, NULL, 1e-15, hhob->cell_basis->projector);
    for (short int f = 0; f < cm->n_fc; f++) {
      fprintf(out, "\n >> Inertial face basis functions (f = %d)\n", f);
      cs_basis_func_fprintf(out, NULL, hhob->face_basis[f]);
      cs_sdm_fprintf(out, NULL, 1e-15, hhob->face_basis[f]->projector);
    }

    cs_hho_builder_compute_grad_reco(cm, cb, hhob);

    fprintf(out,  "\n RHS matrix\n");
    cs_sdm_block_fprintf(out, NULL, 1e-15, cb->aux);

    fprintf(out, "\n Stiffness matrix\n");
    cs_sdm_fprintf(out, NULL, 1e-15, cb->hdg);

    fprintf(out, "\n Gradient Reconstruction matrix\n");
    cs_sdm_block_fprintf(out, NULL, 1e-15, hhob->grad_reco_op);

    cs_hho_builder_diffusion(cm, cb, hhob);
    fprintf(out, "\n Diffusion matrix\n");
    cs_sdm_block_fprintf(out, NULL, 1e-15, cb->loc);

    fprintf(out, "\n Diffusion matrix (Mccgg)\n");
    cs_sdm_block_fprintf(out, NULL, 1e-15, cb->aux);

    fprintf(out, "\n Diffusion matrix (stabilization)\n");
    cs_sdm_block_fprintf(out, NULL, 1e-15, hhob->jstab);

    {
      cs_xdef_analytic_input_t  anai = {.func = _unity, .input = NULL };
      cs_xdef_t  *uni = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                              1,
                                              0, // z_id
                                              0, // state flag
                                              0, // meta flag
                                              &anai);

      anai.func = _linear_xyz;
      cs_xdef_t  *lin = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                              1,
                                              0, // z_id
                                              0, // state flag
                                              0, // meta flag
                                              &anai);

      anai.func = _quadratic_x2;
      cs_xdef_t  *x2 = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                             1,
                                             0, // z_id
                                             0, // state flag
                                             0, // meta flag
                                             &anai);

      cs_real_t  reduction_uni[22], reduction_xyz[22], reduction_x2[22];
      for (int i = 0; i < 3*cm->n_fc+4; i++)
        reduction_uni[i] = reduction_xyz[i] = reduction_x2[i] = 0.0;

      /* cs_hho_builder_reduction_from_analytic(uni, cm, tcur,
                                                cb, hhob, reduction_uni); */
      cs_hho_builder_reduction_from_analytic(lin, cm, tcur, cb, hhob,
                                             reduction_xyz);
      /* cs_hho_builder_reduction_from_analytic(x2 , cm, tcur,
                                                cb, hhob, reduction_x2); */

      fprintf(out, "\n Reduction of polynomial functions.\n"
              "    const   |   linear   | quadratic\n");
      for (int i = 0; i < 3*cm->n_fc+4; i++)
        fprintf(out, " % -5.3e | % -5.3e | % -5.3e\n",
                reduction_uni[i], reduction_xyz[i], reduction_x2[i]);

      /* Evaluate projection at (0.25; 0.25; 0.25) and (1.0, 0.0, 0.0)*/
      cs_real_t  phi_eval[4];
      cs_real_3_t  coord[3] = {{0.25, 0.25, 0.25},
                               {1.00, 0.00, 0.00},
                               {0.10, 0.25, 0.40} };

      for (short int f = 0; f < cm->n_fc; f++) {

        cs_basis_func_t  *fbf = hhob->face_basis[f];
        for (int ic = 0; ic < 3; ic++) {
          fbf->eval_all_at_point(fbf, coord[ic], phi_eval);
          cs_real_t f_eval_at_coord[3] = {0., 0., 0.}, eval_at_coord[3];
          for (int i = 0; i < fbf->size; i++) {
            f_eval_at_coord[0] += reduction_uni[fbf->size*f+i]*phi_eval[i];
            f_eval_at_coord[1] += reduction_xyz[fbf->size*f+i]*phi_eval[i];
            f_eval_at_coord[2] += reduction_x2[fbf->size*f+i]*phi_eval[i];
          }

          cs_real_3_t  a;
          _ortho_proj(cm->face[f], coord[ic], a);

          _unity(0., 1, NULL, a, true, NULL, eval_at_coord);
          _linear_xyz(0., 1, NULL, a, true, NULL, eval_at_coord + 1);
          _quadratic_x2(0., 1, NULL, a, true, NULL, eval_at_coord + 2);

          fprintf(out,
                  "\nface %d (%5.3e, %5.3e, %5.3e) proj(%5.3e, %5.3e, %5.3e)\n",
                  f, coord[ic][0], coord[ic][1], coord[ic][2],
                  a[0], a[1], a[2]);
          fprintf(out, "%10s % -5.3e (exact: % -5.3e)\n",
                  "Unity:", f_eval_at_coord[0], eval_at_coord[0]);
          fprintf(out, "%10s % -5.3e (exact: % -5.3e)\n",
                  "Linear:", f_eval_at_coord[1], eval_at_coord[1]);
          fprintf(out, "%10s % -5.3e (exact: % -5.3e)\n",
                  "QuadX:", f_eval_at_coord[2], eval_at_coord[2]);
        }
      } /* Loop on cell faces */

      int  shift = cm->n_fc*3;
      cs_basis_func_t  *cbf = hhob->cell_basis;
      for (int ic = 0; ic < 3; ic++) {
        cbf->eval_all_at_point(cbf, coord[ic], phi_eval);
        cs_real_t f_eval_at_coord[3] = {0., 0., 0.}, eval_at_coord[3];
        for (int i = 0; i < cbf->size; i++) {
          f_eval_at_coord[0] += reduction_uni[shift+i]*phi_eval[i];
          f_eval_at_coord[1] += reduction_xyz[shift+i]*phi_eval[i];
          f_eval_at_coord[2] += reduction_x2[shift+i]*phi_eval[i];
        }

        _unity(0., 1, NULL, coord[ic], true, NULL, eval_at_coord);
        _linear_xyz(0., 1, NULL, coord[ic], true, NULL, eval_at_coord + 1);
        _quadratic_x2(0., 1, NULL, coord[ic], true, NULL, eval_at_coord + 2);

        fprintf(out, "\ncell for (%5.3e, %5.3e, %5.3e)\n",
                coord[ic][0], coord[ic][1], coord[ic][2]);
        fprintf(out, "%10s % -5.3e (exact: % -5.3e)\n",
                "Unity:", f_eval_at_coord[0], eval_at_coord[0]);
        fprintf(out, "%10s % -5.3e (exact: % -5.3e)\n",
                "Linear:", f_eval_at_coord[1], eval_at_coord[1]);
        fprintf(out, "%10s % -5.3e (exact: % -5.3e)\n",
                "QuadX:", f_eval_at_coord[2], eval_at_coord[2]);
      }

      cs_xdef_free(uni);
      cs_xdef_free(lin);
      cs_xdef_free(x2);
    }

    break;

  case 2:
    fprintf(out,
            "\n ************************************************************\n"
            "                              HHO_P2 scheme\n"
            " ************************************************************\n");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Scheme order not handled yet.\n",
              __func__);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test CDO vertex-based schemes
 *
 * \param[in] out_hho      output file for HHO schemes
 * \param[in] flag         flag storing the order of the polynomial basis
 * \param[in] cm           pointer to a cs_cell_mesh_t structure
 * \param[in] fm           pointer to a cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_main_hho_schemes(FILE             *out_hho,
                  cs_flag_t         flag,
                  cs_cell_mesh_t   *cm,
                  cs_face_mesh_t   *fm)
{
  cs_hho_builder_t  *hhob = NULL;
  cs_cell_builder_t  *cb = NULL;
  cs_cell_sys_t  *csys = NULL;

  cs_hho_scaleq_init_common(flag, quant, connect, time_step, NULL, NULL, NULL);
  cs_hho_scaleq_get(&csys, &cb, &hhob);

  int order = 0;
  if (flag & CS_FLAG_SCHEME_POLY2)
    order = 2;
  else if (flag & CS_FLAG_SCHEME_POLY1)
    order = 1;

  /* Compute the cell tensor inertia */
  cs_real_33_t  inertia;

  cs_compute_inertia_tensor(cm, cm->xc, inertia);
  fprintf(out_hho,
          "                       % .4e % .4e % .4e\n"
          " Cell Inertial tensor  % .4e % .4e % .4e\n"
          "                       % .4e % .4e % .4e\n",
          inertia[0][0], inertia[0][1], inertia[0][2],
          inertia[1][0], inertia[1][1], inertia[1][2],
          inertia[2][0], inertia[2][1], inertia[2][2]);

  _test_basis_functions(out_hho, order, cm, cb);
  _test_hho_schemes(out_hho, order, cm, fm, csys, cb, hhob);

  cs_hho_scaleq_finalize_common();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test CDO vertex-based schemes
 *
 * \param[in]    out       output file
 * \param[in]    case_id   0: hexa; 1: tetra
 * \param[in]    cm        pointer to a cs_cell_mesh_t structure
 * \param[in]    fm        pointer to a cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_main_cdovb_schemes(FILE             *out,
                    int               case_id,
                    cs_cell_mesh_t   *cm,
                    cs_face_mesh_t   *fm)
{
  cs_cell_builder_t  *cb = NULL;
  cs_cell_sys_t  *csys = NULL;

  cs_cdovb_scaleq_init_common(quant, connect, time_step, NULL);
  cs_cdovb_scaleq_get(&csys, &cb);

    /* Initialize a cell view of the BC */
  if (case_id == 0) {

    csys->has_dirichlet = true;
    csys->n_dofs = cm->n_vc;
    for (short int i = 0; i < cm->n_vc; i++) {
      csys->dof_flag[i] = 0;
      csys->dir_values[i] = csys->neu_values[i] = 0;
      csys->rob_values[2*i] = csys->rob_values[2*i+1] = 0.;
    }
    csys->n_bc_faces = 1;
    csys->_f_ids[0] = 4; //f_id = 4
    csys->bf_ids[0] = 4; //f_id = 4
    csys->bf_flag[0] = CS_CDO_BC_DIRICHLET;

    for (short int v = 0; v < fm->n_vf; v++) {
      csys->dir_values[fm->v_ids[v]] = 1.0;
      csys->dof_flag[fm->v_ids[v]] |= CS_CDO_BC_DIRICHLET;
    }

  }
  else if (case_id == 1) {

    csys->has_dirichlet = true;
    csys->n_dofs = cm->n_vc;
    for (short int i = 0; i < cm->n_vc; i++) {
      csys->dof_flag[i] = 0;
      csys->dir_values[i] = csys->neu_values[i] = 0;
      csys->rob_values[2*i] = csys->rob_values[2*i+1] = 0.;
    }
    csys->n_bc_faces = 1;
    csys->_f_ids[0] = 2; //f_id = 2
    csys->bf_ids[0] = 2; //f_id = 2
    csys->bf_flag[0] = CS_CDO_BC_DIRICHLET;

    for (short int v = 0; v < fm->n_vf; v++) {
      csys->dir_values[fm->v_ids[v]] = 1.0;
      csys->dof_flag[fm->v_ids[v]] |= CS_CDO_BC_DIRICHLET;
    }

  }

  _test_cdovb_schemes(out, cm, fm, csys, cb);

  cs_cdovb_scaleq_finalize_common();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence matrix associated to the current cell.
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] div        divergence
 */
/*----------------------------------------------------------------------------*/

static inline void
_compute_divergence(const cs_cell_mesh_t  *cm,
                    cs_real_t              div[])
{
  /* D(\hat{u}) = \frac{1}{|c|} \sum_{f_c} \iota_{fc} u_f.f
   * But, when integrating [[ p, q ]]_{P_c} = |c| p_c q_c
   * Thus, the volume in the divergence drops
   */
  /* It is copied from cdfb_stokes */

  for (short int f = 0; f < cm->n_fc; f++) {

    cs_real_t *_div_f = div + 3 * f;
    const cs_quant_t pfq = cm->face[f];
    const cs_real_t i_f = cm->f_sgn[f] * pfq.meas;

    _div_f[0] = i_f * pfq.unitv[0];
    _div_f[1] = i_f * pfq.unitv[1];
    _div_f[2] = i_f * pfq.unitv[2];

  } /* End for f */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Test quadrature functions on scalar, vector or tensor valued
 *          analytic functions
 *
 * \param[in]    out       output file
 * \param[in]    cm        pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_test_divergence(FILE                     *out,
                 const cs_cell_mesh_t     *cm,
                 cs_cell_builder_t        *cb)
{
  const short int nf = cm->n_fc, totdof = 3*(nf+1);
  const cs_real_t ov = 1. / cm->vol_c;

  cs_real_t *red;
  BFT_MALLOC(red, totdof, cs_real_t);

  cs_real_t *const  div = cb->aux->val;
  _compute_divergence(cm, div);

  /* No cell contribution */
  memset(div + 3* nf, 0, 3*sizeof(cs_real_t));

  /* Error valid only for a unitary cube */
  memset(red, 0, totdof*sizeof(cs_real_t));

  { /* Constant */
    cs_xdef_analytic_input_t  anai = {.func = _unity_vect, .input = NULL };

    cs_xdef_cw_eval_vect_avg_reduction_by_analytic(cm, time_step->t_cur,
                                                   (void*)(&anai),
                                                   CS_QUADRATURE_HIGHEST, red);

    fprintf(out, "\n CDO.FB; DIVERGENCE; P0; |ERR(Div_c)| = %10.6e\n",
            fabs(ov*cs_sdm_dot(totdof, red, div)) );
  }

  { /* Linear */
    cs_xdef_analytic_input_t  anai = {.func = _linear_xyz_vect, .input = NULL };
    cs_real_t  ex_div = 0.0;

    cs_xdef_cw_eval_c_int_by_analytic(cm, time_step->t_cur,
                                         _unity, NULL,
                                         cs_quadrature_tet_1pt_scal,
                                         &ex_div);
    ex_div *= 6. * ov;
    cs_xdef_cw_eval_vect_avg_reduction_by_analytic(cm, time_step->t_cur,
                                                   (void*)(&anai),
                                                   CS_QUADRATURE_HIGHEST, red);

    fprintf(out, " CDO.FB; DIVERGENCE; P1; |ERR(Div_c)| = %10.6e\n",
            fabs(ov*cs_sdm_dot(totdof, red, div) - ex_div) );
  }

  { /* Quadratic */
    cs_xdef_analytic_input_t  anai = {.func = _quadratic_x2_vect,
                                      .input = NULL };
    cs_real_t ex_div = 0.0;

    cs_xdef_cw_eval_c_int_by_analytic(cm, time_step->t_cur,
                                         _linear_xyz, NULL,
                                         cs_quadrature_tet_5pts_scal,
                                         &ex_div);
    ex_div *= 2. *ov;
    cs_xdef_cw_eval_vect_avg_reduction_by_analytic(cm, time_step->t_cur,
                                                   (void*)(&anai),
                                                   CS_QUADRATURE_HIGHEST, red);

    fprintf(out, " CDO.FB; DIVERGENCE; P2; |ERR(Div_c)| = %10.6e\n",
            fabs(ov*cs_sdm_dot(totdof, red, div) - ex_div) );
  }

  { /* Cubic */
    cs_xdef_analytic_input_t  anai = {.func = _cubic_xyz3_vect, .input = NULL };
    cs_real_t ex_div = 0.0;

    cs_xdef_cw_eval_c_int_by_analytic(cm, time_step->t_cur,
                                         _quadratic_xyz2, NULL,
                                         cs_quadrature_tet_5pts_scal, &ex_div);
    ex_div *= 3. * ov;
    cs_xdef_cw_eval_vect_avg_reduction_by_analytic(cm, time_step->t_cur,
                                                   (void*)(&anai),
                                                   CS_QUADRATURE_HIGHEST, red);

    fprintf(out, " CDO.FB; DIVERGENCE; P3; |ERR(Div_c)| = %10.6e\n",
            fabs(ov*cs_sdm_dot(totdof, red, div) - ex_div) );
  }

  { /* Non pol */
    cs_xdef_analytic_input_t  anai = {.func = _nonpoly_vect, .input = NULL };

    memset(red, 0, totdof*sizeof(cs_real_t));
    cs_xdef_cw_eval_vect_avg_reduction_by_analytic(cm, time_step->t_cur,
                                                   (void*)(&anai),
                                                   CS_QUADRATURE_HIGHEST, red);

    const short int  type = (cm->type == FVM_CELL_TETRA);
    const double  inte = non_poly_int[type] * ov;
    const cs_real_t  ex_div = 3*inte;
    for (short int f = 0; f < nf; f++) {
      const double  val = non_poly_int_f[type][f] / cm->face[f].meas;
      red[3*f] = red[3*f+1] = red[3*f+2] = val;
    }
    red[3*nf] = red[3*nf+1] = red[3*nf+2] = inte;

    fprintf(out, " CDO.FB; DIVERGENCE; NP; |ERR(Div_c)| = %10.6e\n",
            fabs(ov*cs_sdm_dot(totdof, red, div) - ex_div) );
  }

  BFT_FREE(red);
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
_main_cdofb_schemes(FILE             *out,
                    cs_cell_mesh_t   *cm)
{
  cs_cell_builder_t  *cb = NULL;
  cs_cell_sys_t  *csys = NULL;

  cs_cdofb_scaleq_init_common(quant, connect, time_step, NULL);
  cs_cdofb_scaleq_get(&csys, &cb);

  /* Initialize a cell view of the algebraic system */
  csys->n_dofs = cm->n_fc + 1;
  cs_sdm_square_init(csys->n_dofs, csys->mat);
  for (short int f = 0; f < cm->n_fc; f++)
    csys->dof_ids[f] = cm->f_ids[f];
  csys->dof_ids[cm->n_fc] = cm->c_id;

  fprintf(out, "\nWORKING ON ");
  if (cm->type == FVM_CELL_HEXA) fprintf(out, "** HEXA **\n");
  else                           fprintf(out, "** TETRA **\n");

  /* Test divergence / gradient */
  _test_divergence(out, cm, cb);

  /* HODGE */
  /* ===== */

  /* COST Hodge operator */
  cs_param_hodge_t  hcost_info = {.is_unity = true,
                                  .is_iso = true,
                                  .inv_pty = false,
                                  .type = CS_PARAM_HODGE_TYPE_FB,
                                  .algo = CS_PARAM_HODGE_ALGO_COST,
                                  .coef = 1.0};

  cs_hodge_fb_get_mass(hcost_info, cm, cb);
  _locmat_dump(out, "\nCDO.FB; HDG.FB.MASS; PERMEABILITY.ISO.UNITY",
               csys->dof_ids, cb->hdg);
  _test_hodge_fb(out, cm, cb->hdg);

  cs_cdofb_scaleq_finalize_common();
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

  hexa = fopen("CDO_tests_Hexa.log", "w");
  tetra = fopen("CDO_tests_Tetra.log", "w");
  hexa_hho0 = fopen("HHO0_Hexa_tests.log", "w");
  hexa_hho1 = fopen("HHO1_Hexa_tests.log", "w");
  hexa_hho2 = fopen("HHO2_Hexa_tests.log", "w");
  tetra_hho0 = fopen("HHO0_Tetra_tests.log", "w");
  tetra_hho1 = fopen("HHO1_Tetra_tests.log", "w");
  tetra_hho2 = fopen("HHO2_Tetra_tests.log", "w");

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

  /* Nothing to do for quant */

  /* Time step */
  BFT_MALLOC(time_step, 1, cs_time_step_t);
  time_step->t_cur = 0.; // Useful when analytic function are called

  cs_source_term_set_shared_pointers(quant, connect);

  /* Allocate local structures */
  cs_cell_mesh_t  *cm = cs_cell_mesh_create(connect);
  cs_face_mesh_t  *fm = cs_face_mesh_create(connect->n_max_vbyf);

  /* ========= */
  /* TEST HEXA */
  /* ========= */

  _define_cm_hexa_unif(1., cm);

  /* Compute the face mesh for the face id 4 */
  cs_face_mesh_build_from_cell_mesh(cm, 4, fm);

  /* Operate several basic tests on CDO-Vb schemes */
  _main_cdovb_schemes(hexa, 0, cm, fm);

  /* Operate several basic tests on CDO-Fb schemes */
  _main_cdofb_schemes(hexa, cm);

  /* _main_hho_schemes(hexa_hho0, CS_FLAG_SCHEME_POLY0, cm, fm); */
  _main_hho_schemes(hexa_hho1, CS_FLAG_SCHEME_POLY1, cm, fm);
  /* _main_hho_schemes(hexa_hho2, CS_FLAG_SCHEME_POLY2, cm, fm); */

  /* ========== */
  /* TEST TETRA */
  /* ========== */

  _define_cm_tetra_ref(1., cm);

  /* Compute the face mesh for the face id 2 */
  cs_face_mesh_build_from_cell_mesh(cm, 2, fm);

  _main_cdovb_schemes(tetra, 1, cm, fm);

  /* Operate several basic tests on CDO-Fb schemes */
  _main_cdofb_schemes(tetra, cm);

  /* _main_hho_schemes(tetra_hho0, CS_FLAG_SCHEME_POLY0, cm, fm); */
  _main_hho_schemes(tetra_hho1, CS_FLAG_SCHEME_POLY1, cm, fm);
  /* _main_hho_schemes(tetra_hho2, CS_FLAG_SCHEME_POLY2, cm, fm); */

  /* Free memory */
  cs_cell_mesh_free(&cm);
  cs_face_mesh_free(&fm);

  BFT_FREE(connect);
  BFT_FREE(time_step);

  fclose(hexa);
  fclose(tetra);
  fclose(hexa_hho0);
  fclose(tetra_hho0);
  fclose(hexa_hho1);
  fclose(tetra_hho1);
  fclose(hexa_hho2);
  fclose(tetra_hho2);

  printf(" --> CDO Tests (Done)\n");
  exit (EXIT_SUCCESS);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
