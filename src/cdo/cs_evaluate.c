/*============================================================================
 * Functions and structures to deal with source term computations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

#define CS_EVALUATE_DBG 0

/* Cases currently handled */
static const cs_flag_t scd_tag =
  CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_DUAL;
static const cs_flag_t scp_tag =
  CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_PRIMAL;
static const cs_flag_t pvp_tag =
  CS_PARAM_FLAG_SCAL | CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL;

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

  double vol_tet = cs_voltet(xv, xe, xf, xc);

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
  int  p;
  double  weight;
  cs_real_3_t  gauss_pts[4];
  cs_get_t  evaluation;

  double result = 0.0;
  double vol_tet = cs_voltet(xv, xe, xf, xc);

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tet_4pts(xv, xe, xf, xc, vol_tet, gauss_pts, &weight);

  for (p = 0; p < 4; p++) {
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
  int  p;
  double  weights[5];
  cs_real_3_t  gauss_pts[5];
  cs_get_t  evaluation;

  double result = 0.0;
  double vol_tet = cs_voltet(xv, xe, xf, xc);

  /* Compute Gauss points and its unique weight */
  cs_quadrature_tet_5pts(xv, xe, xf, xc, vol_tet, gauss_pts, weights);

  for (p = 0; p < 5; p++) {
    ana(tcur, gauss_pts[p], &evaluation);
    result += evaluation.val*weights[p];
  }

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a dual cell (or a portion) of a variable
 *         defined by a user function on a selection of (primal) cells
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      ana        pointer to the analytic function
 * \param[in]      tcur       current physical time of the simulation
 * \param[in]      ml_id      id related to a cs_mesh_location_t struct.
 * \param[in]      quad_type  type of quadrature to use
 * \param[in]      use_subdiv consider or not the subdivision into tetrahedra
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_scd_by_analytic_func(const cs_cdo_quantities_t    *quant,
                      const cs_cdo_connect_t       *connect,
                      cs_analytic_func_t           *ana,
                      double                        tcur,
                      int                           ml_id,
                      cs_quadra_type_t              quad_type,
                      bool                          use_subdiv,
                      double                        values[])
{
  cs_lnum_t  i, j, k, id, e_id, f_id, c_id, v1_id, v2_id, n_elts;
  cs_real_3_t  xc, xv1, xv2;
  cs_quant_t  e, f;

  double  add1 = 0.0, add2 = 0.0;

  const cs_lnum_t  *n_loc_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);
  const cs_sla_matrix_t  *c2f = connect->c2f;
  const cs_sla_matrix_t  *f2e = connect->f2e;

  if (elt_ids == NULL)
    n_elts = quant->n_cells;
  else
    n_elts = n_loc_elts[0];

  if (!use_subdiv) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(" Force subdivision into tetrahedra for evaluating the\n"
               " integral over dual cells of an expression defined by\n"
               " analytic function.\n");
  }

  /* Compute dual volumes */
  for (id = 0; id < n_elts; id++) {

    if (elt_ids == NULL)
      c_id = id;
    else
      c_id = elt_ids[id];

    for (k = 0; k < 3; k++)
      xc[k] = quant->cell_centers[3*c_id+k];

    for (i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) { // Loop on faces

      f_id = c2f->col_id[i];
      f = quant->face[f_id]; /* Get quantities related to this edge */

      for (j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) { // Loop on edges

        e_id = f2e->col_id[j];
        e = quant->edge[e_id]; /* Get quantities related to this edge */
        v1_id = connect->e2v->col_id[2*e_id];
        v2_id = connect->e2v->col_id[2*e_id+1];

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
 * \brief  Compute the integral over primal cells (or a portion) of a variable
 *         defined by a analytical function
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      ana        pointer to the analytic function
 * \param[in]      tcur       current physical time of the simulation
 * \param[in]      ml_id      id related to a cs_mesh_location_t struct.
 * \param[in]      quad_type  type of quadrature to use
 * \param[in]      use_subdiv consider or not the subdivision into tetrahedra
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_scp_by_analytic_func(const cs_cdo_quantities_t    *quant,
                      const cs_cdo_connect_t       *connect,
                      cs_analytic_func_t           *ana,
                      double                        tcur,
                      int                           ml_id,
                      cs_quadra_type_t              quad_type,
                      bool                          use_subdiv,
                      double                        values[])
{
  cs_lnum_t  i, j, k, id, e_id, f_id, c_id, v1_id, v2_id, n_elts;
  cs_real_3_t  xc, xv1, xv2;
  cs_quant_t  f, e;

  double  add = 0.0;

  const cs_lnum_t  *n_loc_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);
  const cs_sla_matrix_t  *c2f = connect->c2f;
  const cs_sla_matrix_t  *f2e = connect->f2e;

  if (elt_ids == NULL)
    n_elts = quant->n_cells;
  else
    n_elts = n_loc_elts[0];

  for (id = 0; id < n_elts; id++) {

    if (elt_ids == NULL)
      c_id = id;
    else
      c_id = elt_ids[id];

    for (k = 0; k < 3; k++)
      xc[k] = quant->cell_centers[3*c_id+k];

    if (use_subdiv) {

      for (i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) { // Loop on faces

        f_id = c2f->col_id[i];
        f = quant->face[f_id]; /* Get quantities related to this face */

        for (j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) { // Loop on edges

          e_id = f2e->col_id[j];
          e = quant->edge[e_id];
          v1_id = connect->e2v->col_id[2*e_id];
          v2_id = connect->e2v->col_id[2*e_id+1];

          for (k = 0; k < 3; k++) {
            xv1[k] = quant->vtx_coord[3*v1_id+k];
            xv2[k] = quant->vtx_coord[3*v2_id+k];
          }

          switch(quad_type) {

          case CS_QUADRATURE_BARY: /* Barycenter of the tetrahedral subdiv. */
            add  = _analytic_quad_tet1(tcur, xv1, e.center, f.center, xc, ana);
            add += _analytic_quad_tet1(tcur, xv2, e.center, f.center, xc, ana);
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

    }
    else { // Do not use subdivision into tetrahedra

      cs_get_t  evaluation;

      if (quad_type != CS_QUADRATURE_BARY) {
        cs_base_warn(__FILE__, __LINE__);
        bft_printf(" Subdivision into tetrahedra is not used for evaluating\n"
                   " the integral over cells of an expression defined by\n"
                   " analytic function but the type of quadrature needs one.\n"
                   " Modify the quadrature type and used bary.\n");
      }

      ana(tcur, xc, &evaluation);
      values[c_id] += evaluation.val * quant->cell_vol[c_id];

    } // subdivision or not

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each primal vertices of a variable defined by a
 *         analytical function
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      ana        pointer to the analytic function
 * \param[in]      tcur       current physical time of the simulation
 * \param[in]      ml_id      id related to a cs_mesh_location_t struct.
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pvp_by_analytic_func(const cs_cdo_quantities_t    *quant,
                      const cs_cdo_connect_t       *connect,
                      cs_analytic_func_t           *ana,
                      double                        tcur,
                      int                           ml_id,
                      double                        values[])
{
  cs_lnum_t  i, j, v_id, c_id;
  cs_get_t  evaluation;

  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);
  const cs_connect_index_t  *c2v = connect->c2v;

  if (elt_ids == NULL) {

    for (v_id = 0; v_id < quant->n_vertices; v_id++) {
      ana(tcur, &(quant->vtx_coord[3*v_id]), &evaluation);
      values[v_id] = evaluation.val;
    }

  }
  else {

    bool  *todo = NULL;

    const cs_lnum_t  *n_loc_elts = cs_mesh_location_get_n_elts(ml_id);

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

  } // elt_ids ?

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a dual cell (or a portion) of a value
 *         defined on a selection of (primal) cells
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      const_val  constant value
 * \param[in]      ml_id      id related to a cs_mesh_location_t struct.
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_scd_by_val(const cs_cdo_quantities_t    *quant,
            const cs_cdo_connect_t       *connect,
            const double                  const_val,
            int                           ml_id,
            double                        values[])
{
  cs_lnum_t  i, j, c_id;

  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);
  const cs_connect_index_t  *c2v = connect->c2v;

  /* Compute dual volumes */

  if (elt_ids == NULL) { // All cells are selected

    for (c_id = 0; c_id < quant->n_cells; c_id++)
      for (j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        values[c2v->ids[j]] += quant->dcell_vol[j]*const_val;

  }
  else { // selection

    for (i = 0; i < n_elts[0]; i++) {

      c_id = elt_ids[i];
      for (j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        values[c2v->ids[j]] += quant->dcell_vol[j]*const_val;

    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a (primal) cell of a value defined on a
 *          selection of (primal) cells
 *
 * \param[in]     quant      additional mesh quantities struct.
 * \param[in]     const_val  constant value
 * \param[in]     ml_id      id related to a cs_mesh_location_t struct.
 * \param[in,out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_scp_by_val(const cs_cdo_quantities_t    *quant,
            const double                  const_val,
            int                           ml_id,
            double                        values[])
{
  cs_lnum_t  i, c_id;

  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  if (elt_ids == NULL) { // All cells are selected

    for (c_id = 0; c_id < quant->n_cells; c_id++)
      values[c_id] = quant->cell_vol[c_id]*const_val;

  }
  else { // selection

    for (i = 0; i < n_elts[0]; i++) {
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
 * \param[in]     quant      additional mesh quantities struct.
 * \param[in]     const_val  constant value
 * \param[in]     ml_id      id related to a cs_mesh_location_t struct.
 * \param[in,out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pvp_by_val(const cs_cdo_quantities_t    *quant,
            const double                  const_val,
            int                           ml_id,
            double                        values[])
{
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  if (elt_ids != NULL)
    bft_error(__FILE__, __LINE__, 0,
              " This situation is not handled yet.\n"
              " Please use a mesh location over the full mesh.");

  for (cs_lnum_t  v_id = 0; v_id < quant->n_vertices; v_id++)
    values[v_id] = const_val;

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a source term or a
 *         boundary condition for instance
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      dof_flag   indicate where the evaluation has to be done
 * \param[in]      ml_id      id related to a cs_mesh_location_t struct.
 * \param[in]      def_type   type of definition
 * \param[in]      quad_type  type of quadrature (not always used)
 * \param[in]      use_subdiv consider or not the subdivision into tetrahedra
 * \param[in]      def        access to the definition of the values
 * \param[in, out] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate(const cs_cdo_quantities_t    *quant,
            const cs_cdo_connect_t       *connect,
            const cs_time_step_t         *time_step,
            cs_flag_t                     dof_flag,
            int                           ml_id,
            cs_param_def_type_t           def_type,
            cs_quadra_type_t              quad_type,
            bool                          use_subdiv,
            cs_def_t                      def,
            double                       *p_values[])
{
  int  stride;
  cs_lnum_t  i, n_ent;

  double  *values = *p_values;

  stride = 0, n_ent = 0;
  if (dof_flag == scd_tag || dof_flag == pvp_tag)
    stride = 1, n_ent = quant->n_vertices;
  else if (dof_flag == scp_tag)
    stride = 1, n_ent = quant->n_cells;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid case. Not able to compute an evaluation.\n"));

  /* Initialize values */
  if (values == NULL)
    BFT_MALLOC(values, n_ent*stride, double);
  for (i = 0; i < n_ent*stride; i++)
    values[i] = 0.0;

  switch (def_type) {

  case CS_PARAM_DEF_BY_VALUE: // constant value (simple case)

    if (dof_flag == scd_tag)
      _scd_by_val(quant, connect, def.get.val, ml_id, values);
    else if (dof_flag == scp_tag)
      _scp_by_val(quant, def.get.val, ml_id, values);
    else if (dof_flag == pvp_tag)
      _pvp_by_val(quant, def.get.val, ml_id, values);
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of definition. Stop evaluation.\n"
                  " This case is not handled yet.\n"));
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION: // constant value (simple case)
    {
      const double  tcur = time_step->t_cur;

      if (dof_flag == scd_tag)
        _scd_by_analytic_func(quant, connect,
                              def.analytic,
                              tcur,
                              ml_id, quad_type, use_subdiv,
                              values);

      else if (dof_flag == scp_tag)
        _scp_by_analytic_func(quant, connect,
                              def.analytic,
                              tcur,
                              ml_id, quad_type, use_subdiv,
                              values);

      else if (dof_flag == pvp_tag)
        _pvp_by_analytic_func(quant, connect,
                              def.analytic,
                              tcur,
                              ml_id,
                              values);

      else
        bft_error(__FILE__, __LINE__, 0,
                  _(" Invalid type of definition. Stop evaluation.\n"
                    " This case is not handled yet.\n"));

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _(" Invalid type of definition.\n"));

  } /* switch def_type */

  /* Return values */
  *p_values = values;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
