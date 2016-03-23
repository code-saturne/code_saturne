/*============================================================================
 * Functions and structures to deal with evaluation of quantities
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_math.h"
#include "cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_evaluate.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;
static const cs_time_step_t  *cs_time_step;

static const char _err_empty_array[] =
  " Array storing the evaluation should be allocated before the call"
  " to this function.";
static const char _err_not_handled[] = " This case is not handled yet.";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron of the barycentric subdiv.
 *         using a barycentric quadrature rule
 *
 * \param[in]  tcur   current physical time of the simulation
 * \param[in]  xv     first point of the tetrahedron
 * \param[in]  xe     second point of the tetrahedron
 * \param[in]  xf     third point of the tetrahedron
 * \param[in]  xc     fourth point of the tetrahedron
 * \param[in]  ana    pointer to the analytic function
 *
 * \return the result of the integration
 */
/*----------------------------------------------------------------------------*/

inline static double
_analytic_quad_tet1(double                 tcur,
                    const cs_real_3_t      xv,
                    const cs_real_3_t      xe,
                    const cs_real_3_t      xf,
                    const cs_real_3_t      xc,
                    cs_analytic_func_t    *ana)
{
  int  k;
  cs_real_3_t  xg;
  cs_get_t  evaluation;

  const double  vol_tet = cs_math_voltet(xv, xe, xf, xc);

  for (k = 0; k < 3; k++)
    xg[k] = 0.25*(xv[k] + xe[k] + xf[k] + xc[k]);

  ana(tcur, xg, &evaluation);

  return vol_tet * evaluation.val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron of the barycentric subdiv.
 *         with a quadrature rule using 4 Gauss points and a unique weight
 *
 * \param[in]  tcur   current physical time of the simulation
 * \param[in]  xv     first point of the tetrahedron
 * \param[in]  xe     second point of the tetrahedron
 * \param[in]  xf     third point of the tetrahedron
 * \param[in]  xc     fourth point of the tetrahedron
 * \param[in]  ana    pointer to the analytic function
 *
 * \return the result of the integration
 */
/*----------------------------------------------------------------------------*/

inline static double
_analytic_quad_tet4(double                tcur,
                    const cs_real_3_t     xv,
                    const cs_real_3_t     xe,
                    const cs_real_3_t     xf,
                    const cs_real_3_t     xc,
                    cs_analytic_func_t   *ana)
{
  double  weight;
  cs_real_3_t  gauss_pts[4];
  cs_get_t  evaluation;

  double  result = 0.0;
  const double  vol_tet = cs_math_voltet(xv, xe, xf, xc);

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tet_4pts(xv, xe, xf, xc, vol_tet, gauss_pts, &weight);

  for (int p = 0; p < 4; p++) {
    ana(tcur, gauss_pts[p], &evaluation);
    result += evaluation.val;
  }
  result *= weight;

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a tetrahedron of the barycentric subdiv.
 *         with a quadrature rule using 5 Gauss points and 5 weights
 *
 * \param[in]  tcur   current physical time of the simulation
 * \param[in]  xv     first point of the tetrahedron
 * \param[in]  xe     second point of the tetrahedron
 * \param[in]  xf     third point of the tetrahedron
 * \param[in]  xc     fourth point of the tetrahedron
 * \param[in]  ana    pointer to the analytic function
 *
 * \return the result of the integration
 */
/*----------------------------------------------------------------------------*/

inline static double
_analytic_quad_tet5(double                tcur,
                    const cs_real_3_t     xv,
                    const cs_real_3_t     xe,
                    const cs_real_3_t     xf,
                    const cs_real_3_t     xc,
                    cs_analytic_func_t   *ana)
{
  double  weights[5];
  cs_real_3_t  gauss_pts[5];
  cs_get_t  evaluation;

  double  result = 0.0;
  const double  vol_tet = cs_math_voltet(xv, xe, xf, xc);

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tet_5pts(xv, xe, xf, xc, vol_tet, gauss_pts, weights);

  for (int p = 0; p < 5; p++) {
    ana(tcur, gauss_pts[p], &evaluation);
    result += evaluation.val*weights[p];
  }

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over dual cells of a scalar density field
 *         defined by an analytical function on a selection of (primal) cells
 *
 * \param[in]      ana         pointer to the analytic function
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in]      quad_type   type of quadrature to use
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_dcsd_by_analytic(cs_analytic_func_t       *ana,
                  const cs_lnum_t           n_elts,
                  const cs_lnum_t          *elt_ids,
                  cs_quadra_type_t          quad_type,
                  double                    values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_sla_matrix_t  *c2f = connect->c2f;
  const cs_sla_matrix_t  *f2e = connect->f2e;
  const double  tcur = cs_time_step->t_cur;

  /* Compute dual volumes */
  for (cs_lnum_t  id = 0; id < n_elts; id++) {

    const cs_lnum_t  c_id = (elt_ids == NULL) ? id : elt_ids[id];
    const cs_real_t  *xc = quant->cell_centers + 3*c_id;

    for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

      const cs_lnum_t  f_id = c2f->col_id[i];
      const cs_quant_t  f = quant->face[f_id];

      for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

        const cs_lnum_t  e_id = f2e->col_id[j];
        const cs_quant_t  e = quant->edge[e_id];
        const cs_lnum_t  v1_id = connect->e2v->col_id[2*e_id];
        const cs_lnum_t  v2_id = connect->e2v->col_id[2*e_id+1];
        const cs_real_t  *xv1 = quant->vtx_coord + 3*v1_id;
        const cs_real_t  *xv2 = quant->vtx_coord + 3*v2_id;

        double  add1 = 0.0, add2 = 0.0;

        switch(quad_type) {

        case CS_QUADRATURE_BARY: /* Barycenter of the tetrahedral subdiv. */
          add1 = _analytic_quad_tet1(tcur, xv1, e.center, f.center, xc, ana);
          add2 = _analytic_quad_tet1(tcur, xv2, e.center, f.center, xc, ana);
          break;
        case CS_QUADRATURE_HIGHER: /* Quadrature with a unique weight */
          add1 = _analytic_quad_tet4(tcur, xv1, e.center, f.center, xc, ana);
          add2 = _analytic_quad_tet4(tcur, xv1, e.center, f.center, xc, ana);
          break;
        case CS_QUADRATURE_HIGHEST: /* Most accurate quadrature available */
          add1 = _analytic_quad_tet5(tcur, xv1, e.center, f.center, xc, ana);
          add2 = _analytic_quad_tet5(tcur, xv1, e.center, f.center, xc, ana);
          break;
        default:
          bft_error(__FILE__, __LINE__, 0, _("Invalid quadrature type.\n"));

        } /* Quad rule */

        values[v1_id] += add1;
        values[v2_id] += add2;

      } // Loop on edges

    } // Loop on faces

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over primal cells of a scalar density field
 *         defined by an analytical function on a selection of (primal) cells
 *
 * \param[in]      ana         pointer to the analytic function
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in]      quad_type   type of quadrature to use
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pcsd_by_analytic(cs_analytic_func_t       *ana,
                  const cs_lnum_t           n_elts,
                  const cs_lnum_t          *elt_ids,
                  cs_quadra_type_t          quad_type,
                  double                    values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_sla_matrix_t  *c2f = connect->c2f;
  const cs_sla_matrix_t  *f2e = connect->f2e;
  const double  tcur = cs_time_step->t_cur;

  for (cs_lnum_t  id = 0; id < n_elts; id++) {

    const cs_lnum_t  c_id = (elt_ids == NULL) ? id : elt_ids[id];
    const cs_real_t  *xc = quant->cell_centers + 3*c_id;

    for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

      const cs_lnum_t  f_id = c2f->col_id[i];
      const cs_quant_t  f = quant->face[f_id];

      for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

        const cs_lnum_t  e_id = f2e->col_id[j];
        const cs_lnum_t  v1_id = connect->e2v->col_id[2*e_id];
        const cs_lnum_t  v2_id = connect->e2v->col_id[2*e_id+1];
        const cs_real_t  *xv1 = quant->vtx_coord + 3*v1_id;
        const cs_real_t  *xv2 = quant->vtx_coord + 3*v2_id;

        double  add = 0.0;

        switch(quad_type) {

        case CS_QUADRATURE_BARY: /* Barycenter of the tetrahedral subdiv. */
          add  = _analytic_quad_tet1(tcur, xv1, xv2, f.center, xc, ana);
          break;
        case CS_QUADRATURE_HIGHER: /* Quadrature with a unique weight */
          add = _analytic_quad_tet4(tcur, xv1, xv2, f.center, xc, ana);
          break;
        case CS_QUADRATURE_HIGHEST: /* Most accurate quadrature available */
          add = _analytic_quad_tet5(tcur, xv1, xv2, f.center, xc, ana);
          break;
        default:
          bft_error(__FILE__, __LINE__, 0, _("Invalid quadrature type.\n"));

        } /* Quad rule */

        values[c_id] += add;

      } // Loop on edges

    } // Loop on faces

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a dual cell (or a portion) of a value
 *         defined on a selection of (primal) cells
 *
 * \param[in]      const_val   constant value
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_dcsd_by_value(const double       const_val,
               const cs_lnum_t    n_elts,
               const cs_lnum_t   *elt_ids,
               double             values[])
{
  const cs_connect_index_t  *c2v = cs_cdo_connect->c2v;
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *dual_vol = quant->dcell_vol; /* scan by c2v */

  if (elt_ids == NULL) {

    assert(n_elts == quant->n_cells);
    for (cs_lnum_t c_id = 0; c_id < n_elts; c_id++)
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        values[c2v->ids[j]] += dual_vol[j]*const_val;

  }
  else { /* Loop on selected cells */

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t  c_id = elt_ids[i];
      for (cs_lnum_t  j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        values[c2v->ids[j]] += dual_vol[j]*const_val;
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a (primal) cell of a value related to
 *         scalar density field
 *
 * \param[in]      const_val   constant value
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pcsd_by_value(const double       const_val,
               const cs_lnum_t    n_elts,
               const cs_lnum_t   *elt_ids,
               double             values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  if (elt_ids == NULL) { /* All the support entities are selected */
# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      values[c_id] = quant->cell_vol[c_id]*const_val;
  }

  else { /* Loop on selected cells */
# pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t  c_id = elt_ids[i];
      values[c_id] = quant->cell_vol[c_id]*const_val;
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each primal faces for a scalar potential
 *         defined by an analytical function on a selection of (primal) cells
 *
 * \param[in]      ana         pointer to the analytic function
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pfsp_by_analytic(cs_analytic_func_t    *ana,
                  const cs_lnum_t        n_elts,
                  const cs_lnum_t       *elt_ids,
                  double                 values[])
{
  cs_get_t  result;

  const double  tcur = cs_time_step->t_cur;
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_sla_matrix_t  *c2f = cs_cdo_connect->c2f;

  /* Initialize todo array */
  bool  *todo = NULL;

  BFT_MALLOC(todo, quant->n_vertices, bool);

# pragma omp parallel for if (quant->n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < quant->n_faces; f_id++)
    todo[f_id] = true;

  for (cs_lnum_t i = 0; i < n_elts; i++) { // Loop on selected cells

    cs_lnum_t  c_id = elt_ids[i];

    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      cs_lnum_t  f_id = c2f->col_id[j];
      if (todo[f_id]) {
        ana(tcur, quant->face[f_id].center, &result);
        values[f_id] = result.val;
        todo[f_id] = false;
      }

    } // Loop on cell vertices

  } // Loop on selected cells

  BFT_FREE(todo);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each primal vertices for a scalar potential
 *         defined by an analytical function on a selection of (primal) cells
 *
 * \param[in]      ana         pointer to the analytic function
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pvsp_by_analytic(cs_analytic_func_t    *ana,
                  const cs_lnum_t        n_elts,
                  const cs_lnum_t       *elt_ids,
                  double                 values[])
{
  cs_get_t  result;

  const double  tcur = cs_time_step->t_cur;
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_connect_index_t  *c2v = cs_cdo_connect->c2v;

  /* Initialize todo array */
  bool  *todo = NULL;

  BFT_MALLOC(todo, quant->n_vertices, bool);

# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++)
    todo[v_id] = true;

  for (cs_lnum_t i = 0; i < n_elts; i++) { // Loop on selected cells

    cs_lnum_t  c_id = elt_ids[i];

    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

      cs_lnum_t  v_id = c2v->ids[j];
      if (todo[v_id]) {
        ana(tcur, quant->vtx_coord + 3*v_id, &result);
        values[v_id] = result.val;
        todo[v_id] = false;
      }

    } // Loop on cell vertices

  } // Loop on selected cells

  BFT_FREE(todo);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each primal faces for a scalar potential
 *
 * \param[in]      const_val   constant value
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in, out] values      pointer to the array storing the values
 */
/*----------------------------------------------------------------------------*/

static void
_pfsp_by_value(const double       const_val,
               cs_lnum_t          n_elts,
               const cs_lnum_t   *elt_ids,
               double             values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_sla_matrix_t  *c2f = cs_cdo_connect->c2f;

  /* Initialize todo array */
  bool  *todo = NULL;

  BFT_MALLOC(todo, quant->n_vertices, bool);

# pragma omp parallel for if (quant->n_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < quant->n_faces; f_id++)
    todo[f_id] = true;

  for (cs_lnum_t i = 0; i < n_elts; i++) { // Loop on selected cells

    cs_lnum_t  c_id = elt_ids[i];

    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      cs_lnum_t  f_id = c2f->col_id[j];
      if (todo[f_id])
        values[f_id] = const_val, todo[f_id] = false;

    } // Loop on cell vertices

  } // Loop on selected cells

  BFT_FREE(todo);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a value to each DoF such that a given quantity is put inside
 *         the volume associated to the list of cells
 *
 * \param[in]      quantity_val  amount of quantity to distribute
 * \param[in]      n_loc_elts    number of elements to consider
 * \param[in]      elt_ids       pointer to the list od selected ids
 * \param[in, out] values        pointer to the array storing the values
 */
/*----------------------------------------------------------------------------*/

static void
_pvsp_by_qov(const double       quantity_val,
             cs_lnum_t          n_elts,
             const cs_lnum_t   *elt_ids,
             double             values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *dc_vol = quant->dcell_vol;
  const cs_connect_index_t  *c2v = cs_cdo_connect->c2v;
  const cs_sla_matrix_t  *c2f = cs_cdo_connect->c2f;
  const cs_sla_matrix_t  *f2c = cs_cdo_connect->f2c;
  const cs_sla_matrix_t  *f2e = cs_cdo_connect->f2e;
  const cs_sla_matrix_t  *e2v = cs_cdo_connect->e2v;

  /* Initialize todo array */
  bool  *cell_tag = NULL, *vtx_tag = NULL;

  BFT_MALLOC(cell_tag, quant->n_cells, bool);
  BFT_MALLOC(vtx_tag, quant->n_vertices, bool);

# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++)
    vtx_tag[v_id] = false;
# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
    cell_tag[c_id] = false;

  /* First pass: flag cells and vertices */
# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) { // Loop on selected cells

    const cs_lnum_t  c_id = elt_ids[i];

    cell_tag[c_id] = true;
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
      vtx_tag[c2v->ids[j]] = true;

  } // Loop on selected cells

  /* Second pass: detect cells at the frontier of the selection */
# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) { // Loop on selected cells

    const cs_lnum_t  c_id = elt_ids[i];

    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      const cs_lnum_t  f_id = c2f->col_id[j];

      bool is_ext_face = false;
      for (cs_lnum_t l = f2c->idx[f_id]; l < f2c->idx[f_id+1]; l++)
        if (!cell_tag[f2c->col_id[l]]) is_ext_face = true;

      if (is_ext_face) {
        for (cs_lnum_t l = f2e->idx[f_id]; l < f2e->idx[f_id+1]; l++) {

          const cs_lnum_t  e_id = f2e->col_id[l];
          const cs_lnum_t  shift_e = 2*e_id;

          vtx_tag[e2v->col_id[shift_e]] = false;
          vtx_tag[e2v->col_id[shift_e+1]] = false;

        } // Loop on face edges
      } // This face belongs to the frontier of the selection (only interior)

    } // Loop on cell faces

  } // Loop on selected cells

  /* Third pass: compute the (really) available volume */
  double  volume = 0.;
# pragma omp parallel for reduction(+:volume) if (n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) { // Loop on selected cells

    const cs_lnum_t  c_id = elt_ids[i];

    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
      if (vtx_tag[c2v->ids[j]])
        volume += dc_vol[j];

  } // Loop on selected cells

  double val_to_set = quantity_val;
  if (volume > 0)
    val_to_set /= volume;

# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++)
    if (vtx_tag[v_id])
      values[v_id] = val_to_set;

  BFT_FREE(cell_tag);
  BFT_FREE(vtx_tag);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each primal vertices for a scalar potential
 *
 * \param[in]      const_val   constant value
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in, out] values      pointer to the array storing the values
 */
/*----------------------------------------------------------------------------*/

static void
_pvsp_by_value(const double       const_val,
               cs_lnum_t          n_elts,
               const cs_lnum_t   *elt_ids,
               double             values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_connect_index_t  *c2v = cs_cdo_connect->c2v;

  /* Initialize todo array */
  bool  *todo = NULL;

  BFT_MALLOC(todo, quant->n_vertices, bool);

# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++)
    todo[v_id] = true;

  for (cs_lnum_t i = 0; i < n_elts; i++) { // Loop on selected cells

    cs_lnum_t  c_id = elt_ids[i];

    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

      cs_lnum_t  v_id = c2v->ids[j];
      if (todo[v_id])
        values[v_id] = const_val, todo[v_id] = false;

    } // Loop on cell vertices

  } // Loop on selected cells

  BFT_FREE(todo);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                const cs_cdo_connect_t       *connect,
                                const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
  cs_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value related to each DoF in the case of a density field
 *         The value defined by the analytic function is by unity of volume
 *
 * \param[in]      dof_flag    indicate where the evaluation has to be done
 * \param[in]      ml_id       id related to a cs_mesh_location_t structure
 * \param[in]      ana         accessor to values thanks to a function pointer
 * \param[in]      quad_type   type of quadrature (not always used)
 * \param[in, out] retval      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_density_from_analytic(cs_flag_t              dof_flag,
                                  int                    ml_id,
                                  cs_analytic_func_t    *ana,
                                  cs_quadra_type_t       quad_type,
                                  double                 retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array);

  /* Retrieve information from mesh location structures */
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  /* Sanity checks */
  assert(n_elts != NULL);
  cs_mesh_location_type_t  ml_type = cs_mesh_location_get_type(ml_id);
  if (elt_ids != NULL && ml_type != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0, _err_not_handled);

  /* Perform the evaluation */
  if (dof_flag & CS_FLAG_SCAL) { /* DoF is scalar-valued */

    if (cs_cdo_same_support(dof_flag, cs_cdo_primal_cell))
      _pcsd_by_analytic(ana, n_elts[0], elt_ids, quad_type, retval);

    else if (cs_cdo_same_support(dof_flag, cs_cdo_dual_cell))
      _dcsd_by_analytic(ana, n_elts[0], elt_ids, quad_type, retval);

    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled);

  }
  else
    bft_error(__FILE__, __LINE__, 0, _err_not_handled);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value related to each DoF in the case of a density field
 *         Accessor to the value is by unit of volume
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      ml_id     id related to a cs_mesh_location_t structure
 * \param[in]      get       accessor to the constant value to set
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_density_from_value(cs_flag_t       dof_flag,
                               int             ml_id,
                               cs_get_t        get,
                               double          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array);

  /* Retrieve information from mesh location structures */
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  /* Sanity checks */
  assert(n_elts != NULL);
  cs_mesh_location_type_t  ml_type = cs_mesh_location_get_type(ml_id);
  if (elt_ids != NULL && ml_type != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0, _err_not_handled);

  /* Perform the evaluation */
  if (dof_flag & CS_FLAG_SCAL) { /* DoF is scalar-valued */

    if (cs_cdo_same_support(dof_flag, cs_cdo_primal_cell))
      _pcsd_by_value(get.val, n_elts[0], elt_ids, retval);

    else if (cs_cdo_same_support(dof_flag, cs_cdo_dual_cell))
      _dcsd_by_value(get.val, n_elts[0], elt_ids, retval);

    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled);

  }
  else
    bft_error(__FILE__, __LINE__, 0, _err_not_handled);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a quantity defined by analytic
 *         function for all the degrees of freedom
 *
 * \param[in]      dof_flag    indicate where the evaluation has to be done
 * \param[in]      ml_id       id related to a cs_mesh_location_t structure
 * \param[in]      ana         accessor to values thanks to a function pointer
 * \param[in, out] retval      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_from_analytic(cs_flag_t              dof_flag,
                                    int                    ml_id,
                                    cs_analytic_func_t    *ana,
                                    double                 retval[])
{
  cs_get_t  result;

  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array);

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const double  tcur = cs_time_step->t_cur;

  /* Retrieve information from mesh location structures */
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  /* Sanity checks */
  assert(n_elts != NULL);
  cs_mesh_location_type_t  ml_type = cs_mesh_location_get_type(ml_id);
  if (elt_ids != NULL && ml_type != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0, _err_not_handled);

  /* Perform the evaluation */
  if (dof_flag & CS_FLAG_SCAL) { /* DoF is scalar-valued */

    if (cs_cdo_same_support(dof_flag, cs_cdo_primal_vtx)) {

      if (elt_ids == NULL) { /* All the support entities are selected */
# pragma omp parallel for private(result) if(quant->n_vertices > CS_THR_MIN)
        for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {
          ana(tcur, quant->vtx_coord + 3*v_id, &result);
          retval[v_id] = result.val;
        }
      }
      else
        _pvsp_by_analytic(ana, n_elts[0], elt_ids, retval);

    } /* Located at primal vertices */

    else if (cs_cdo_same_support(dof_flag, cs_cdo_primal_face)) {

      if (elt_ids == NULL) { /* All the support entities are selected */
# pragma omp parallel for  private(result) if(quant->n_faces > CS_THR_MIN)
        for (cs_lnum_t f_id = 0; f_id < quant->n_faces; f_id++) {
          ana(tcur, quant->face[f_id].center, &result);
          retval[f_id] = result.val;
        }
      }
      else
        _pfsp_by_analytic(ana, n_elts[0], elt_ids, retval);

    } /* Located at primal faces */

    else if (cs_cdo_same_support(dof_flag, cs_cdo_primal_cell) ||
             cs_cdo_same_support(dof_flag, cs_cdo_dual_vtx)) {

      if (elt_ids == NULL) { /* All the support entities are selected */
# pragma omp parallel for  private(result) if(quant->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
          ana(tcur, quant->cell_centers + 3*c_id, &result);
          retval[c_id] = result.val;
        }
      }
      else
        for (cs_lnum_t i = 0; i < n_elts[0]; i++) { // Loop on selected cells
          cs_lnum_t  c_id = elt_ids[i];
          ana(tcur, quant->cell_centers + 3*c_id, &result);
          retval[c_id] = result.val;
        }

    } /* Located at primal cells or dual vertices */

    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled);

  }
  else
    bft_error(__FILE__, __LINE__, 0, _err_not_handled);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a value to each DoF in the case of a potential field in order
 *         to put a given quantity inside the volume associated to ml_id
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      ml_id     id related to a cs_mesh_location_t structure
 * \param[in]      get       accessor to the constant value related to the
 *                           quantity to put in the volume spanned by ml_id
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_from_qov(cs_flag_t       dof_flag,
                               int             ml_id,
                               cs_get_t        get,
                               double          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array);

  /* Retrieve information from mesh location structures */
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  /* Sanity checks */
  assert(n_elts != NULL);
  cs_mesh_location_type_t  ml_type = cs_mesh_location_get_type(ml_id);
  if (elt_ids != NULL && ml_type != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0, _err_not_handled);

  /* Perform the evaluation */
  bool check = false;
  if (dof_flag & CS_FLAG_SCAL) { /* DoF is scalar-valued */

    if (cs_cdo_same_support(dof_flag, cs_cdo_primal_vtx))
      if (elt_ids != NULL) {
        _pvsp_by_qov(get.val, n_elts[0], elt_ids, retval);
        check = true;
      }

  } /* Located at primal vertices */

  if (!check)
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop evaluating a potential from 'quantity over volume'.\n"
                " This situation is not handled yet."));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store the value related to each DoF in the case of a potential field
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      ml_id     id related to a cs_mesh_location_t structure
 * \param[in]      get       accessor to the constant value to set
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_from_value(cs_flag_t       dof_flag,
                                 int             ml_id,
                                 cs_get_t        get,
                                 double          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array);

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  /* Retrieve information from mesh location structures */
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  /* Sanity checks */
  assert(n_elts != NULL);
  cs_mesh_location_type_t  ml_type = cs_mesh_location_get_type(ml_id);
  if (elt_ids != NULL && ml_type != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0, _err_not_handled);

  /* Perform the evaluation */
  if (dof_flag & CS_FLAG_SCAL) { /* DoF is scalar-valued */

    if (cs_cdo_same_support(dof_flag, cs_cdo_primal_vtx)) {

      if (elt_ids == NULL) { /* All the support entities are selected */
# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
        for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++)
          retval[v_id] = get.val;
      }
      else
        _pvsp_by_value(get.val, n_elts[0], elt_ids, retval);

    } /* Located at primal vertices */

    else if (cs_cdo_same_support(dof_flag, cs_cdo_primal_face)) {

      if (elt_ids == NULL) { /* All the support entities are selected */
# pragma omp parallel for if (quant->n_faces > CS_THR_MIN)
        for (cs_lnum_t f_id = 0; f_id < quant->n_faces; f_id++)
          retval[f_id] = get.val;
      }
      else
        _pfsp_by_value(get.val, n_elts[0], elt_ids, retval);

    } /* Located at primal faces */

    else if (cs_cdo_same_support(dof_flag, cs_cdo_primal_cell) ||
             cs_cdo_same_support(dof_flag, cs_cdo_dual_vtx)) {

      if (elt_ids == NULL) { /* All the support entities are selected */
# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
          retval[c_id] = get.val;
      }
      else
        for (cs_lnum_t i = 0; i < n_elts[0]; i++) // Loop on selected cells
          retval[elt_ids[i]] = get.val;

    } /* Located at primal cells or dual vertices */

    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled);

  }
  else
    bft_error(__FILE__, __LINE__, 0, _err_not_handled);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
