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

/* Cases currently handled */
static const cs_flag_t  scd_dof =
  CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_DUAL;
static const cs_flag_t  scp_dof =
  CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_PRIMAL;
static const cs_flag_t pvp_dof =
  CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;
static const cs_time_step_t  *cs_time_step;

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
_analytic_quad_tet1(double                tcur,
                    cs_real_3_t           xv,
                    cs_real_3_t           xe,
                    cs_real_3_t           xf,
                    cs_real_3_t           xc,
                    cs_analytic_func_t   *ana)
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
_analytic_quad_tet4(double               tcur,
                    cs_real_3_t          xv,
                    cs_real_3_t          xe,
                    cs_real_3_t          xf,
                    cs_real_3_t          xc,
                    cs_analytic_func_t  *ana)
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
_analytic_quad_tet5(double              tcur,
                    cs_real_3_t         xv,
                    cs_real_3_t         xe,
                    cs_real_3_t         xf,
                    cs_real_3_t         xc,
                    cs_analytic_func_t *ana)
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
 * \brief  Compute the integral over a dual cell of a variable defined by an
 *         analytical function on a selection of (primal) cells
 *
 * \param[in]      ana         pointer to the analytic function
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in]      quad_type   type of quadrature to use
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_scd_dof_from_analytic(cs_analytic_func_t       *ana,
                       const cs_lnum_t          *n_loc_elts,
                       const cs_lnum_t          *elt_ids,
                       cs_quadra_type_t          quad_type,
                       double                    values[])
{
  cs_lnum_t  i, j, k, c_id;
  cs_real_3_t  xc, xv1, xv2;

  double  add1 = 0.0, add2 = 0.0;
  double  tcur = cs_time_step->t_cur;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_sla_matrix_t  *c2f = connect->c2f;
  const cs_sla_matrix_t  *f2e = connect->f2e;

  /* Compute dual volumes */
  for (cs_lnum_t  id = 0; id < n_loc_elts[0]; id++) {

    if (elt_ids == NULL)
      c_id = id;
    else
      c_id = elt_ids[id];

    for (k = 0; k < 3; k++)
      xc[k] = quant->cell_centers[3*c_id+k];

    for (i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) { // Loop on faces

      cs_lnum_t  f_id = c2f->col_id[i];
      cs_quant_t  f = quant->face[f_id];

      for (j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) { // Loop on edges

        cs_lnum_t  e_id = f2e->col_id[j];
        cs_quant_t  e = quant->edge[e_id];
        cs_lnum_t  v1_id = connect->e2v->col_id[2*e_id];
        cs_lnum_t  v2_id = connect->e2v->col_id[2*e_id+1];

        for (k = 0; k < 3; k++) {
          xv1[k] = quant->vtx_coord[3*v1_id+k];
          xv2[k] = quant->vtx_coord[3*v2_id+k];
        }

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

      } // edges
    } // faces

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over primal cells of a variable defined by an
 *         analytical function on a selection of (primal) cells
 *
 * \param[in]      ana         pointer to the analytic function
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in]      quad_type   type of quadrature to use
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_scp_dof_from_analytic(cs_analytic_func_t       *ana,
                       const cs_lnum_t          *n_loc_elts,
                       const cs_lnum_t          *elt_ids,
                       cs_quadra_type_t          quad_type,
                       double                    values[])
{
  cs_lnum_t  i, j, k, c_id;
  cs_real_3_t  xc, xv1, xv2;

  double  add = 0.0;
  double  tcur = cs_time_step->t_cur;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_sla_matrix_t  *c2f = connect->c2f;
  const cs_sla_matrix_t  *f2e = connect->f2e;

  for (cs_lnum_t  id = 0; id < n_loc_elts[0]; id++) {

    if (elt_ids == NULL)
      c_id = id;
    else
      c_id = elt_ids[id];

    for (k = 0; k < 3; k++)
      xc[k] = quant->cell_centers[3*c_id+k];

    for (i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) { // Loop on faces

      cs_lnum_t  f_id = c2f->col_id[i];
      cs_quant_t  f = quant->face[f_id];

      for (j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) { // Loop on edges

        cs_lnum_t  e_id = f2e->col_id[j];
        cs_lnum_t  v1_id = connect->e2v->col_id[2*e_id];
        cs_lnum_t  v2_id = connect->e2v->col_id[2*e_id+1];

        for (k = 0; k < 3; k++) {
          xv1[k] = quant->vtx_coord[3*v1_id+k];
          xv2[k] = quant->vtx_coord[3*v2_id+k];
        }

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

      } // edges
    } // faces

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each primal vertices of a variable defined by an
 *         analytical function on a selection of (primal) cells
 *
 * \param[in]      ana         pointer to the analytic function
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pvp_dof_from_analytic(cs_analytic_func_t    *ana,
                       const cs_lnum_t       *n_loc_elts,
                       const cs_lnum_t       *elt_ids,
                       double                 values[])
{
  cs_lnum_t  i, j, v_id, c_id;
  cs_get_t  evaluation;

  double  tcur = cs_time_step->t_cur;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_connect_index_t  *c2v = cs_cdo_connect->c2v;

  if (elt_ids == NULL) {

    for (v_id = 0; v_id < quant->n_vertices; v_id++) {
      ana(tcur, &(quant->vtx_coord[3*v_id]), &evaluation);
      values[v_id] = evaluation.val;
    }

  }
  else {

    bool  *todo = NULL;

    BFT_MALLOC(todo, quant->n_vertices, bool);
    for (v_id = 0; v_id < quant->n_vertices; v_id++)
      todo[v_id] = true;

    for (i = 0; i < n_loc_elts[0]; i++) { // Loop on selected cells

      c_id = elt_ids[i];

      for (j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

        v_id = c2v->ids[j];
        if (todo[v_id]) {
          ana(tcur, &(quant->vtx_coord[3*v_id]), &evaluation);
          values[v_id] = evaluation.val;
          todo[v_id] = false;
        }

      } // Loop on cell vertices

    } // Loop on selected cells

    BFT_FREE(todo);

  } // elt_ids ?

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
_scd_dof_from_value(const double         const_val,
                    const cs_lnum_t     *n_loc_elts,
                    const cs_lnum_t     *elt_ids,
                    double               values[])
{
  cs_lnum_t  i, j, c_id;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_connect_index_t  *c2v = cs_cdo_connect->c2v;

  if (elt_ids == NULL) { // All cells are selected

    for (c_id = 0; c_id < quant->n_cells; c_id++)
      for (j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        values[c2v->ids[j]] += quant->dcell_vol[j]*const_val;

  }
  else { // selection

    for (i = 0; i < n_loc_elts[0]; i++) {
      c_id = elt_ids[i];
      for (j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        values[c2v->ids[j]] += quant->dcell_vol[j]*const_val;
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a (primal) cell of a value defined on a
 *         selection of (primal) cells
 *
 * \param[in]      const_val   constant value
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_scp_dof_from_value(const double         const_val,
                    const cs_lnum_t     *n_loc_elts,
                    const cs_lnum_t     *elt_ids,
                    double               values[])
{
  cs_lnum_t  i, c_id;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  if (elt_ids == NULL) { // All cells are selected

    for (c_id = 0; c_id < quant->n_cells; c_id++)
      values[c_id] = quant->cell_vol[c_id]*const_val;

  }
  else { // selection

    for (i = 0; i < n_loc_elts[0]; i++) {
      c_id = elt_ids[i];
      values[c_id] = quant->cell_vol[c_id]*const_val;
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each primal vertices of a variable defined by a
 *         constant value
 *
 * \param[in]      const_val   constant value
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list od selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pvp_dof_from_value(const double         const_val,
                    const cs_lnum_t     *n_loc_elts,
                    const cs_lnum_t     *elt_ids,
                    double               values[])
{
  cs_lnum_t  i, j, v_id;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_connect_index_t  *c2v = cs_cdo_connect->c2v;

  if (elt_ids == NULL)
    for (v_id = 0; v_id < quant->n_vertices; v_id++)
      values[v_id] = const_val;

  else {

    bool  *todo = NULL;

    BFT_MALLOC(todo, quant->n_vertices, bool);
    for (v_id = 0; v_id < quant->n_vertices; v_id++)
      todo[v_id] = true;

    for (i = 0; i < n_loc_elts[0]; i++) { // Loop on selected cells

      cs_lnum_t  c_id = elt_ids[i];

      for (j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

        v_id = c2v->ids[j];
        if (todo[v_id]) {
          values[v_id] = const_val;
          todo[v_id] = false;
        }

      } // Loop on cell vertices

    } // Loop on selected cells

    BFT_FREE(todo);

  } // elt_ids ?

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
 * \brief  Compute the contribution related to a quantity defined by analytic
 *         function for all the degrees of freedom
 *
 * \param[in]      dof_flag    indicate where the evaluation has to be done
 * \param[in]      ml_id       id related to a cs_mesh_location_t structure
 * \param[in]      ana         accessor to values thanks to a function pointer
 * \param[in]      quad_type   type of quadrature (not always used)
 * \param[in, out] retval      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_from_analytic(cs_flag_t              dof_flag,
                          int                    ml_id,
                          cs_analytic_func_t    *ana,
                          cs_quadra_type_t       quad_type,
                          double                 retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " array storing the result is not allocated.");

  /* Retrieve information from mesh location structures */
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  /* Sanity checks */
  assert(n_elts != NULL);

  if (dof_flag == scd_dof)
    _scd_dof_from_analytic(ana, n_elts, elt_ids, quad_type, retval);
  else if (dof_flag == scp_dof)
    _scp_dof_from_analytic(ana, n_elts, elt_ids, quad_type, retval);
  else if (dof_flag == pvp_dof)
    _pvp_dof_from_analytic(ana, n_elts, elt_ids, retval);
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of degrees of freedom.\n"
                " This case is not handled yet."));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a quantity defined by analytic
 *         function for all the degrees of freedom
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      ml_id     id related to a cs_mesh_location_t structure
 * \param[in]      value     constant value used for computing the contribution
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_from_value(cs_flag_t       dof_flag,
                       int             ml_id,
                       double          value,
                       double          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " array storing the result is not allocated.");

  /* Retrieve information from mesh location structures */
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  /* Sanity checks */
  assert(n_elts != NULL);

  if (dof_flag == scd_dof)
    _scd_dof_from_value(value, n_elts, elt_ids, retval);
  else if (dof_flag == scp_dof)
    _scp_dof_from_value(value, n_elts, elt_ids, retval);
  else if (dof_flag == pvp_dof)
    _pvp_dof_from_value(value, n_elts, elt_ids, retval);
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of degree of freedom.\n"
                " This case is not handled yet."));

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
