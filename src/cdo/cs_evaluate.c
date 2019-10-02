/*============================================================================
 * Functions and structures to deal with evaluation of quantities
 *============================================================================*/

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

#include "cs_array_reduce.h"
#include "cs_halo.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_parall.h"
#include "cs_range_set.h"
#include "cs_volume_zone.h"
#include "cs_quadrature.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_evaluate.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

#define _dp3  cs_math_3_dot_product

/*=============================================================================
 * Local static variables
 *============================================================================*/

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;

static const char _err_empty_array[] =
  " %s: Array storing the evaluation should be allocated before the call"
  " to this function.";
static const char _err_not_handled[] = " %s: Case not handled yet.";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Parallel synchronization of the local reduction operations
 *
 * \param[in]      dim     local array dimension (max: 3)
 * \param[in, out] min     resulting min array (size: dim, or 4 if dim = 3)
 * \param[in, out] max     resulting max array (size: dim, or 4 if dim = 3)
 * \param[in, out] wsum    (weighted) sum array (size: dim, or 4 if dim = 3)
 * \param[in, out] asum    (weighted) sum of absolute values (same size as wsum)
 * \param[in, out] ssum    (weighted) sum of squared values (same size as wsum)
 */
/*----------------------------------------------------------------------------*/

static void
_synchronize_reduction(int              dim,
                       cs_real_t       *min,
                       cs_real_t       *max,
                       cs_real_t       *wsum,
                       cs_real_t       *asum,
                       cs_real_t       *ssum)
{
  if (cs_glob_n_ranks < 2)
    return; /* Nothing to do */

  /* Min/Max */
  if (dim == 1) {

    cs_real_t  minmax[2] = {-min[0], max[0]};
    cs_parall_max(2, CS_REAL_TYPE, minmax);

    min[0] = -minmax[0];
    max[0] = minmax[1];

  }
  else {

    assert(dim == 3);
    cs_real_t  minmax[8];
    for (int i = 0; i < 4; i++)
      minmax[i] = -min[i], minmax[4+i] = max[i];

    cs_parall_max(8, CS_REAL_TYPE, minmax);

    for (int i = 0; i < 4; i++)
      min[i] = -minmax[i], max[i] = minmax[4+i];

  }

  /* Sums */
  if (dim == 1) {

    cs_real_t  sums[3] = {wsum[0], asum[0], ssum[0]};
    cs_parall_sum(3, CS_REAL_TYPE, sums);

    wsum[0] = sums[0];
    asum[0] = sums[1];
    ssum[0] = sums[2];

  }
  else {

    assert(dim == 3);
    cs_real_t  sums[12];
    for (int i = 0; i < 4; i++)
      sums[i] = wsum[i], sums[4+i] = asum[i], sums[8+i] = ssum[i];

    cs_parall_sum(12, CS_REAL_TYPE, sums);

    for (int i = 0; i < 4; i++)
      wsum[i] = sums[i], asum[i] = sums[4+i], ssum[i] = sums[8+i];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Unmarked vertices belonging to the frontier of the cell selection
 *
 * \param[in]      c_id       id of the cell to treat
 * \param[in]      c_tags     tag for each cell
 * \param[in, out] v_tags     tag for each vertex
 */
/*----------------------------------------------------------------------------*/

static void
_untag_frontier_vertices(cs_lnum_t         c_id,
                         const cs_lnum_t   c_tags[],
                         cs_lnum_t         v_tags[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_lnum_t  *f2v_idx = m->i_face_vtx_idx;
  const cs_lnum_t  *f2v_lst = m->i_face_vtx_lst;
  const cs_adjacency_t  *c2f = cs_cdo_connect->c2f;

  for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

    const cs_lnum_t  f_id = c2f->ids[j];
    if (f_id < m->n_i_faces) { /* interior face */

      if (c_tags[m->i_face_cells[f_id][0]] == 0 ||
          c_tags[m->i_face_cells[f_id][1]] == 0) {
        for (cs_lnum_t i = f2v_idx[f_id]; i < f2v_idx[f_id+1]; i++)
          v_tags[f2v_lst[i]] = 0;  /* untag */
      }

    } /* This face belongs to the frontier of the selection (only interior) */

  } /* Loop on cell faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Tag vertex and cell entities associated to the selected cells
 *         Perform the parallel synchronization.
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       pointer to the list of selected ids
 * \param[in, out] v_tags        pointer to the array storing the vertex tags
 */
/*----------------------------------------------------------------------------*/

static void
_tag_geometric_entities(cs_lnum_t          n_elts,
                        const cs_lnum_t   *elt_ids,
                        cs_lnum_t          v_tags[])
{
  const cs_mesh_t  *m = cs_glob_mesh;
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;

  cs_lnum_t  *c_tags = NULL;
  BFT_MALLOC(c_tags, m->n_cells_with_ghosts, cs_lnum_t);

  if (n_elts < n_cells) { /* Only some cells are selected */

    memset(v_tags, 0, n_vertices * sizeof(cs_lnum_t));
    memset(c_tags, 0, m->n_cells_with_ghosts * sizeof(cs_lnum_t));

    /* First pass: flag cells and vertices */
#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) { /* Loop on selected cells */

      const cs_lnum_t  c_id = elt_ids[i];
      c_tags[c_id] = 1;
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        v_tags[c2v->ids[j]] = -1; /* activated */

    } /* Loop on selected cells */

  }
  else { /* All cells are selected */

    assert(n_cells == n_elts);

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
      v_tags[v_id] = -1;

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_tags[c_id] = 1;
    for (cs_lnum_t c_id = n_cells; c_id < m->n_cells_with_ghosts; c_id++)
      c_tags[c_id] = 0;

  }

  if (m->halo != NULL)
    cs_halo_sync_num(m->halo, CS_HALO_STANDARD, c_tags);

  /* Second pass: detect cells at the frontier of the selection */
  for (cs_lnum_t i = 0; i < n_elts; i++) {
    const cs_lnum_t  c_id = (n_elts == n_cells) ? i : elt_ids[i];
    _untag_frontier_vertices(c_id, c_tags, v_tags);
  }

  /* Not needed anymore */
  BFT_FREE(c_tags);

  /* Handle parallelism (always the scalar interface) */
  if (cs_glob_n_ranks > 1)
    cs_interface_set_max(cs_cdo_connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         n_vertices,
                         1,           /* stride */
                         true,        /* interlace, not useful here */
                         CS_LNUM_TYPE,
                         (void *)v_tags);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a value to each DoF such that a given quantity is put inside
 *         the volume associated to the list of cells
 *         Case of primal vertices for scalar-valued quantities.
 *
 * \param[in]      quantity_val  amount of quantity to distribute
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       pointer to the list of selected ids
 * \param[in, out] v_vals        pointer to the array storing the v_vals
 */
/*----------------------------------------------------------------------------*/

static void
_pvsp_by_qov(const cs_real_t    quantity_val,
             cs_lnum_t          n_elts,
             const cs_lnum_t   *elt_ids,
             cs_real_t          v_vals[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_real_t  *dc_vol = quant->dcell_vol;
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;

  cs_lnum_t  *v_tags = NULL;

  BFT_MALLOC(v_tags, n_vertices, cs_lnum_t);

  /* Tag selected vertices and cells */
  _tag_geometric_entities(n_elts, elt_ids, v_tags);

  /* Third pass: compute the (really) available volume */
  double  volume_marked = 0.;

#   pragma omp parallel for reduction(+:volume_marked) if (n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) { /* Loop on selected cells */

    const cs_lnum_t  c_id = (n_elts == n_cells) ? i : elt_ids[i];
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
      if (v_tags[c2v->ids[j]] == -1) /* activated */
        volume_marked += dc_vol[j];   /* | dual_cell cap cell | */

  } /* Loop on selected cells */

  /* Handle parallelism */
  if (cs_glob_n_ranks > 1)
    cs_parall_sum(1, CS_DOUBLE, &volume_marked);

  cs_real_t  val_to_set = quantity_val;
  if (volume_marked > 0)
    val_to_set /= volume_marked;

  if (n_elts < n_cells) { /* Only some cells are selected */

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
      if (v_tags[v_id] == -1)
        v_vals[v_id] = val_to_set;

  }
  else { /* All cells are selected */

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
      v_vals[v_id] = val_to_set;

  }

  BFT_FREE(v_tags);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a value to each DoF such that a given quantity is put inside
 *         the volume associated to the list of cells
 *         Case of primal vertices and cells for scalar-valued quantities.
 *
 * \param[in]      quantity_val  amount of quantity to distribute
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       pointer to the list of selected ids
 * \param[in, out] v_vals        pointer to the array storing the vertex values
 * \param[in, out] c_vals        pointer to the array storing the cell values
 */
/*----------------------------------------------------------------------------*/

static void
_pvcsp_by_qov(const cs_real_t    quantity_val,
              cs_lnum_t          n_elts,
              const cs_lnum_t   *elt_ids,
              cs_real_t          v_vals[],
              cs_real_t          c_vals[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_real_t  *dc_vol = quant->dcell_vol;
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;

  cs_lnum_t  *v_tags = NULL;

  BFT_MALLOC(v_tags, n_vertices, cs_lnum_t);

  /* Tag selected vertices and cells */
  _tag_geometric_entities(n_elts, elt_ids, v_tags);

  /* Compute the (really) available volume:
     - 1/4 of the cell volume is associated to the cell unkonwn
     - 3/4 of the dual cell volume is associated to the vertex unknown
  */
  double  volume_marked = 0.;

#   pragma omp parallel for reduction(+:volume_marked) if (n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) { /* Loop on selected cells */

    const cs_lnum_t  c_id = (n_elts == n_cells) ? i : elt_ids[i];

    volume_marked += 0.25 * quant->cell_vol[c_id];
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
      if (v_tags[c2v->ids[j]] == -1) /* activated */
        volume_marked += 0.75 * dc_vol[j];   /* 3/4 * | dual_cell cap cell | */

  } /* Loop on selected cells */

  /* Handle parallelism */
  cs_parall_sum(1, CS_DOUBLE, &volume_marked);

  cs_real_t  val_to_set = quantity_val;
  if (volume_marked > 0)
    val_to_set /= volume_marked;

  if (n_elts < n_cells) { /* Only some cells are selected */

    assert(elt_ids != NULL);

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      c_vals[elt_ids[i]] = val_to_set;

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
      if (v_tags[v_id] == -1)
        v_vals[v_id] = val_to_set;

  }
  else { /* All cells are selected */

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      c_vals[c_id] = val_to_set;

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
      v_vals[v_id] = val_to_set;

  }

  BFT_FREE(v_tags);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over dual cells of a scalar-valued density field
 *         defined by an analytical function on a selection of (primal) cells
 *
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      ana               pointer to the analytic function
 * \param[in]      input             NULL or pointer cast on-the-fly
 * \param[in]      n_elts            number of elements to consider
 * \param[in]      elt_ids           pointer to the list of selected ids
 * \param[in]      compute_integral  function pointer
 * \param[in, out] values            pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_dcsd_by_analytic(cs_real_t                        time_eval,
                  cs_analytic_func_t              *ana,
                  void                            *input,
                  const cs_lnum_t                  n_elts,
                  const cs_lnum_t                 *elt_ids,
                  cs_quadrature_tetra_integral_t  *compute_integral,
                  cs_real_t                        values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;

  /* Computation over dual volumes */
  for (cs_lnum_t  id = 0; id < n_elts; id++) {

    const cs_lnum_t  c_id = (elt_ids == NULL) ? id : elt_ids[id];
    const cs_real_t  *xc = quant->cell_centers + 3*c_id;

    for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

      const cs_lnum_t  f_id = c2f->ids[i];
      const cs_real_t  *xf = cs_quant_get_face_center(f_id, quant);

      for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

        const cs_lnum_t  _2e = 2*f2e->ids[j];
        const cs_lnum_t  v1 = connect->e2v->ids[_2e];
        const cs_lnum_t  v2 = connect->e2v->ids[_2e+1];
        const cs_real_t  *xv1 = quant->vtx_coord + 3*v1;
        const cs_real_t  *xv2 = quant->vtx_coord + 3*v2;

        cs_real_3_t  xe;
        for (int k = 0; k < 3; k++) xe[k] = 0.5 * (xv1[k] + xv2[k]);

        compute_integral(time_eval, xv1, xe, xf, xc, quant->dcell_vol[v1],
                         ana, input, values + v1);
        compute_integral(time_eval, xv2, xe, xf, xc, quant->dcell_vol[v2],
                         ana, input, values + v2);

      } /* Loop on edges */

    } /* Loop on faces */

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over dual cells of a vector-valued density field
 *         defined by an analytical function on a selection of (primal) cells
 *
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      ana               pointer to the analytic function
 * \param[in]      input             NULL or pointer cast on-the-fly
 * \param[in]      n_elts            number of elements to consider
 * \param[in]      elt_ids           pointer to the list of selected ids
 * \param[in]      compute_integral  function pointer
 * \param[in, out] values            pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_dcvd_by_analytic(cs_real_t                        time_eval,
                  cs_analytic_func_t              *ana,
                  void                            *input,
                  const cs_lnum_t                  n_elts,
                  const cs_lnum_t                 *elt_ids,
                  cs_quadrature_tetra_integral_t  *compute_integral,
                  cs_real_t                        values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;
  const int  dim = 3;

  /* Computation over dual volumes */
  for (cs_lnum_t  id = 0; id < n_elts; id++) {

    const cs_lnum_t  c_id = (elt_ids == NULL) ? id : elt_ids[id];
    const cs_real_t  *xc = quant->cell_centers + 3*c_id;

    for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

      const cs_lnum_t  f_id = c2f->ids[i];
      const cs_real_t  *xf = cs_quant_get_face_center(f_id, quant);

      for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

        const cs_lnum_t  _2e = 2*f2e->ids[j];
        const cs_lnum_t  v1 = connect->e2v->ids[_2e];
        const cs_lnum_t  v2 = connect->e2v->ids[_2e+1];
        const cs_real_t  *xv1 = quant->vtx_coord + 3*v1;
        const cs_real_t  *xv2 = quant->vtx_coord + 3*v2;

        cs_real_3_t  xe;
        for (int k = 0; k < 3; k++) xe[k] = 0.5 * (xv1[k] + xv2[k]);

        compute_integral(time_eval, xv1, xe, xf, xc, quant->dcell_vol[v1],
                         ana, input, values + dim*v1);
        compute_integral(time_eval, xv2, xe, xf, xc, quant->dcell_vol[v2],
                         ana, input, values + dim*v2);

      }  /* Loop on edges */

    }  /* Loop on faces */

  }  /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over primal cells of a scalar-valued density
 *         field defined by an analytical function on a selection of (primal)
 *         cells
 *
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      ana               pointer to the analytic function
 * \param[in]      input             NULL or pointer cast on-the-fly
 * \param[in]      n_elts            number of elements to consider
 * \param[in]      elt_ids           pointer to the list of selected ids
 * \param[in]      compute_integral  function pointer
 * \param[in, out] values            pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pcsd_by_analytic(cs_real_t                        time_eval,
                  cs_analytic_func_t              *ana,
                  void                            *input,
                  const cs_lnum_t                  n_elts,
                  const cs_lnum_t                 *elt_ids,
                  cs_quadrature_tetra_integral_t  *compute_integral,
                  cs_real_t                        values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *xv = quant->vtx_coord;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t id = 0; id < n_elts; id++) {

    const cs_lnum_t  c_id = (elt_ids == NULL) ? id : elt_ids[id];
    if (connect->cell_type[c_id] == FVM_CELL_TETRA) {

      const cs_lnum_t  *v = connect->c2v->ids + connect->c2v->idx[c_id];

      compute_integral(time_eval,
                       xv+3*v[0], xv+3*v[1], xv+3*v[2], xv+3*v[3],
                       quant->cell_vol[c_id],
                       ana, input, values + c_id);

    }
    else {

      const cs_real_t  *xc = quant->cell_centers + 3*c_id;

      for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

        const cs_lnum_t  f_id = c2f->ids[i];
        const cs_quant_t  pfq = cs_quant_set_face(f_id, quant);
        const double  hfco =
          cs_math_1ov3 * cs_math_3_dot_product(pfq.unitv,
                                                   quant->dedge_vector+3*i);
        const cs_lnum_t  start = f2e->idx[f_id], end = f2e->idx[f_id+1];

        if (end - start == 3) {

          cs_lnum_t v0, v1, v2;
          cs_connect_get_next_3_vertices(connect->f2e->ids, connect->e2v->ids,
                                         start, &v0, &v1, &v2);
          compute_integral(time_eval, xv + 3*v0, xv + 3*v1, xv + 3*v2, xc,
                           hfco * pfq.meas,
                           ana, input, values + c_id);
        }
        else {

          for (cs_lnum_t j = start; j < end; j++) {

            const cs_lnum_t  _2e = 2*f2e->ids[j];
            const cs_lnum_t  v1 = connect->e2v->ids[_2e];
            const cs_lnum_t  v2 = connect->e2v->ids[_2e+1];

            compute_integral(time_eval, xv + 3*v1, xv + 3*v2, pfq.center, xc,
                             hfco*cs_math_surftri(xv+3*v1, xv+3*v2, pfq.center),
                             ana, input, values + c_id);

          } /* Loop on edges */

        } /* Current face is triangle or not ? */

      } /* Loop on faces */

    } /* Not a tetrahedron */

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over primal cells of a vector-valued density
 *         field defined by an analytical function on a selection of (primal)
 *         cells
 *
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      ana               pointer to the analytic function
 * \param[in]      input             NULL or pointer cast on-the-fly
 * \param[in]      n_elts            number of elements to consider
 * \param[in]      elt_ids           pointer to the list of selected ids
 * \param[in]      compute_integral  function pointer
 * \param[in, out] values            pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pcvd_by_analytic(cs_real_t                        time_eval,
                  cs_analytic_func_t              *ana,
                  void                            *input,
                  const cs_lnum_t                  n_elts,
                  const cs_lnum_t                 *elt_ids,
                  cs_quadrature_tetra_integral_t  *compute_integral,
                  cs_real_t                        values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *xv = quant->vtx_coord;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;
  const int  dim = 3;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t  id = 0; id < n_elts; id++) {

    const cs_lnum_t  c_id = (elt_ids == NULL) ? id : elt_ids[id];
    if (connect->cell_type[c_id] == FVM_CELL_TETRA) {

      const cs_lnum_t  *v = connect->c2v->ids + connect->c2v->idx[c_id];

      compute_integral(time_eval,
                       xv+3*v[0], xv+3*v[1], xv+3*v[2], xv+3*v[3],
                       quant->cell_vol[c_id],
                       ana, input, values + dim*c_id);

    }
    else {

      const cs_real_t  *xc = quant->cell_centers + 3*c_id;

      for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

        const cs_lnum_t  f_id = c2f->ids[i];
        const cs_quant_t  pfq = cs_quant_set_face(f_id, quant);
        const double hfc = cs_math_1ov3 *
                cs_math_3_dot_product(pfq.unitv, quant->dedge_vector+3*i);
        const cs_lnum_t start = f2e->idx[f_id], end = f2e->idx[f_id+1];

        if (end - start == 3) {

          cs_lnum_t v0, v1, v2;
          cs_connect_get_next_3_vertices(connect->f2e->ids, connect->e2v->ids,
                                         start, &v0, &v1, &v2);
          compute_integral(time_eval,xv + 3*v0, xv + 3*v1, xv + 3*v2, xc,
                           hfc * pfq.meas,
                           ana, input, values + 3*c_id);
        }
        else {

          for (cs_lnum_t j = start; j < end; j++) {

            const cs_lnum_t  _2e = 2*f2e->ids[j];
            const cs_lnum_t  v1 = connect->e2v->ids[_2e];
            const cs_lnum_t  v2 = connect->e2v->ids[_2e+1];

            compute_integral(time_eval, xv + 3*v1, xv + 3*v2, pfq.center, xc,
                             hfc*cs_math_surftri(xv+3*v1, xv+3*v2, pfq.center),
                             ana, input, values + 3*c_id);

          } /* Loop on edges */

        } /* Current face is triangle or not ? */

      } /* Loop on faces */

    } /* Not a tetrahedron */

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the average over primal cells of a scalar field defined
 *         by an analytical function on a selection of (primal) cells
 *
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      ana               pointer to the analytic function
 * \param[in]      input             NULL or pointer cast on-the-fly
 * \param[in]      n_loc_elts        number of elements to consider
 * \param[in]      elt_ids           pointer to the list of selected ids
 * \param[in]      compute_integral  function pointer
 * \param[in, out] values            pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pcsa_by_analytic(cs_real_t                        time_eval,
                  cs_analytic_func_t              *ana,
                  void                            *input,
                  const cs_lnum_t                  n_elts,
                  const cs_lnum_t                 *elt_ids,
                  cs_quadrature_tetra_integral_t  *compute_integral,
                  cs_real_t                        values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *xv = quant->vtx_coord;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t id = 0; id < n_elts; id++) {

    const cs_lnum_t  c_id = (elt_ids == NULL) ? id : elt_ids[id];
    if (connect->cell_type[c_id] == FVM_CELL_TETRA) {

      const cs_lnum_t  *v = connect->c2v->ids + connect->c2v->idx[c_id];

      compute_integral(time_eval,
                       xv+3*v[0], xv+3*v[1], xv+3*v[2], xv+3*v[3],
                       quant->cell_vol[c_id],
                       ana, input, values + c_id);

    }
    else {

      const cs_real_t  *xc = quant->cell_centers + 3*c_id;

      for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

        const cs_lnum_t  f_id = c2f->ids[i];
        const cs_quant_t  pfq = cs_quant_set_face(f_id, quant);
        const double  hfco =
          cs_math_1ov3 * cs_math_3_dot_product(pfq.unitv,
                                                   quant->dedge_vector+3*i);
        const cs_lnum_t  start = f2e->idx[f_id], end = f2e->idx[f_id+1];

        if (end - start == 3) {

          cs_lnum_t v0, v1, v2;
          cs_connect_get_next_3_vertices(connect->f2e->ids, connect->e2v->ids,
                                         start, &v0, &v1, &v2);
          compute_integral(time_eval, xv + 3*v0, xv + 3*v1, xv + 3*v2, xc,
                           hfco * pfq.meas,
                           ana, input, values + c_id);
        }
        else {

          for (cs_lnum_t j = start; j < end; j++) {

            const cs_lnum_t  _2e = 2*f2e->ids[j];
            const cs_lnum_t  v1 = connect->e2v->ids[_2e];
            const cs_lnum_t  v2 = connect->e2v->ids[_2e+1];

            compute_integral(time_eval, xv + 3*v1, xv + 3*v2, pfq.center, xc,
                             hfco*cs_math_surftri(xv+3*v1, xv+3*v2, pfq.center),
                             ana, input, values + c_id);

          } /* Loop on edges */

        } /* Current face is triangle or not ? */

      } /* Loop on faces */

    } /* Not a tetrahedron */

    /* Average */
    values[c_id] /= quant->cell_vol[c_id];

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the average over primal cells of a vector field defined
 *         by an analytical function on a selection of (primal) cells
 *
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      ana               pointer to the analytic function
 * \param[in]      input             NULL or pointer cast on-the-fly
 * \param[in]      n_loc_elts        number of elements to consider
 * \param[in]      elt_ids           pointer to the list of selected ids
 * \param[in]      compute_integral  function pointer
 * \param[in, out] values            pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

/* Note: the only difference from the scalar version is that there's 3*c_id in
 * the values. Consider merging the two
 */

static void
_pcva_by_analytic(cs_real_t                        time_eval,
                  cs_analytic_func_t              *ana,
                  void                            *input,
                  const cs_lnum_t                  n_elts,
                  const cs_lnum_t                 *elt_ids,
                  cs_quadrature_tetra_integral_t  *compute_integral,
                  cs_real_t                        values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *xv = quant->vtx_coord;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t id = 0; id < n_elts; id++) {

    const cs_lnum_t  c_id = (elt_ids == NULL) ? id : elt_ids[id];
    cs_real_t *val_i = values + 3*c_id;

    if (connect->cell_type[c_id] == FVM_CELL_TETRA) {

      const cs_lnum_t  *v = connect->c2v->ids + connect->c2v->idx[c_id];

      compute_integral(time_eval,
                       xv+3*v[0], xv+3*v[1], xv+3*v[2], xv+3*v[3],
                       quant->cell_vol[c_id],
                       ana, input, val_i);

    }
    else {

      const cs_real_t  *xc = quant->cell_centers + 3*c_id;

      for (cs_lnum_t i = c2f->idx[c_id]; i < c2f->idx[c_id+1]; i++) {

        const cs_lnum_t  f_id = c2f->ids[i];
        const cs_quant_t  pfq = cs_quant_set_face(f_id, quant);
        const double  hfco =
          cs_math_1ov3 * cs_math_3_dot_product(pfq.unitv,
                                               quant->dedge_vector+3*i);
        const cs_lnum_t  start = f2e->idx[f_id], end = f2e->idx[f_id+1];

        if (end - start == 3) {

          cs_lnum_t v0, v1, v2;
          cs_connect_get_next_3_vertices(connect->f2e->ids, connect->e2v->ids,
                                         start, &v0, &v1, &v2);
          compute_integral(time_eval, xv + 3*v0, xv + 3*v1, xv + 3*v2, xc,
                           hfco * pfq.meas,
                           ana, input, val_i);

        }
        else {

          for (cs_lnum_t j = start; j < end; j++) {

            const cs_lnum_t  _2e = 2*f2e->ids[j];
            const cs_lnum_t  v1 = connect->e2v->ids[_2e];
            const cs_lnum_t  v2 = connect->e2v->ids[_2e+1];

            compute_integral(time_eval, xv + 3*v1, xv + 3*v2, pfq.center, xc,
                             hfco*cs_math_surftri(xv+3*v1, xv+3*v2, pfq.center),
                             ana, input, val_i);

          } /* Loop on edges */

        } /* Current face is triangle or not ? */

      } /* Loop on faces */

    } /* Not a tetrahedron */

    const double _overvol = 1./quant->cell_vol[c_id];
    for (int k = 0; k < 3; k++) val_i[k] *= _overvol;

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a dual cell (or a portion) of a value
 *         defined on a selection of (primal) cells
 *
 * \param[in]      const_val   constant value
 * \param[in]      n_elts      number of elements to consider
 * \param[in]      elt_ids     pointer to the list of selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_dcsd_by_value(const cs_real_t    const_val,
               const cs_lnum_t    n_elts,
               const cs_lnum_t   *elt_ids,
               cs_real_t          values[])
{
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
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
 * \brief  Compute the integral over a dual cell (or a portion) of a
 *         vector-valued density field defined on a selection of (primal) cells
 *
 * \param[in]      const_vec   constant vector
 * \param[in]      n_elts      number of elements to consider
 * \param[in]      elt_ids     pointer to the list of selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_dcvd_by_value(const cs_real_t    const_vec[3],
               const cs_lnum_t    n_elts,
               const cs_lnum_t   *elt_ids,
               cs_real_t          values[])
{
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
  const cs_real_t  *dual_vol = cs_cdo_quant->dcell_vol; /* scan by c2v */

  if (elt_ids == NULL) {

    for (cs_lnum_t c_id = 0; c_id < n_elts; c_id++) {
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {
        const cs_lnum_t  v_id = c2v->ids[j];
        const cs_real_t  vol_vc = dual_vol[j];

        values[3*v_id   ] += vol_vc*const_vec[0];
        values[3*v_id +1] += vol_vc*const_vec[1];
        values[3*v_id +2] += vol_vc*const_vec[2];

      }
    }

  }
  else { /* Loop on selected cells */

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  c_id = elt_ids[i];
      for (cs_lnum_t  j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {
        const cs_lnum_t  v_id = c2v->ids[j];
        const cs_real_t  vol_vc = dual_vol[j];

        values[3*v_id   ] += vol_vc*const_vec[0];
        values[3*v_id +1] += vol_vc*const_vec[1];
        values[3*v_id +2] += vol_vc*const_vec[2];
      }
    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a (primal) cell of a value related to
 *         scalar density field
 *
 * \param[in]      const_val   constant value
 * \param[in]      n_elts      number of elements to consider
 * \param[in]      elt_ids     pointer to the list of selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pcsd_by_value(const cs_real_t    const_val,
               const cs_lnum_t    n_elts,
               const cs_lnum_t   *elt_ids,
               cs_real_t          values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  if (elt_ids == NULL) { /* All the support entities are selected */
#   pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      values[c_id] = quant->cell_vol[c_id]*const_val;
  }

  else { /* Loop on selected cells */
#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t  c_id = elt_ids[i];
      values[c_id] = quant->cell_vol[c_id]*const_val;
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the average over a (primal) cell of a scalar field
 *
 * \param[in]      const_val   constant value
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list of selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static inline void
_pcsa_by_value(const cs_real_t    const_val,
               const cs_lnum_t    n_elts,
               const cs_lnum_t   *elt_ids,
               cs_real_t          values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  if (elt_ids == NULL) { /* All the support entities are selected */
#   pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      values[c_id] = const_val;
  }

  else { /* Loop on selected cells */
#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t  c_id = elt_ids[i];
      values[c_id] = const_val;
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the integral over a (primal) cell of a vector-valued
 *         density field
 *
 * \param[in]      const_vec   constant values
 * \param[in]      n_elts      number of elements to consider
 * \param[in]      elt_ids     pointer to the list of selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pcvd_by_value(const cs_real_t     const_vec[3],
               const cs_lnum_t     n_elts,
               const cs_lnum_t    *elt_ids,
               cs_real_t           values[])
{
  const cs_real_t  *vol = cs_cdo_quant->cell_vol;

  if (elt_ids == NULL) { /* All the support entities are selected */
#   pragma omp parallel for if (cs_cdo_quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < cs_cdo_quant->n_cells; c_id++) {
      const cs_real_t  vol_c = vol[c_id];
      cs_real_t  *val_c = values + 3*c_id;

      val_c[0] = vol_c * const_vec[0];
      val_c[1] = vol_c * const_vec[1];
      val_c[2] = vol_c * const_vec[2];
    }
  }

  else { /* Loop on selected cells */
#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  c_id = elt_ids[i];
      const cs_real_t  vol_c = vol[c_id];
      cs_real_t  *val_c = values + 3*c_id;

      val_c[0] = vol_c * const_vec[0];
      val_c[1] = vol_c * const_vec[1];
      val_c[2] = vol_c * const_vec[2];
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the average over a (primal) cell of a vector-valued field
 *
 * \param[in]      const_vec   constant values
 * \param[in]      n_loc_elts  number of elements to consider
 * \param[in]      elt_ids     pointer to the list of selected ids
 * \param[in, out] values      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static inline void
_pcva_by_value(const cs_real_t     const_vec[3],
               const cs_lnum_t     n_elts,
               const cs_lnum_t    *elt_ids,
               cs_real_t           values[])
{
  if (elt_ids == NULL) { /* All the support entities are selected */
#   pragma omp parallel for if (cs_cdo_quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t c_id = 0; c_id < cs_cdo_quant->n_cells; c_id++) {
      memcpy(values + 3*c_id, const_vec, 3*sizeof(cs_real_t));
    }
  }

  else { /* Loop on selected cells */
#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  c_id = elt_ids[i];
      memcpy(values+3*c_id, const_vec, 3*sizeof(cs_real_t));
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the average at each primal faces for a scalar potential
 *         defined by an analytical function on a selection of (primal) cells
 *
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      ana               pointer to the analytic function
 * \param[in]      input             NULL or pointer cast on-the-fly
 * \param[in]      n_loc_elts        number of elements to consider
 * \param[in]      elt_ids           pointer to the list of selected ids
 * \param[in]      compute_integral  function pointer
 * \param[in, out] values            pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pfsa_by_analytic(cs_real_t                        time_eval,
                  cs_analytic_func_t              *ana,
                  void                            *input,
                  const cs_lnum_t                  n_elts,
                  const cs_lnum_t                 *elt_ids,
                  cs_quadrature_tria_integral_t   *compute_integral,
                  cs_real_t                        values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_adjacency_t  *f2e = cs_cdo_connect->f2e;
  const cs_adjacency_t  *e2v = cs_cdo_connect->e2v;
  const cs_real_t  *xv = quant->vtx_coord;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t f = 0; f < n_elts; f++) {

    const cs_lnum_t  f_id = (elt_ids == NULL) ? f : elt_ids[f];
    const cs_quant_t pfq = cs_quant_set_face(f_id, quant);
    const cs_lnum_t  start_idx = f2e->idx[f_id], end_idx = f2e->idx[f_id+1];
    cs_real_t *val_i = values + f_id;

    switch (end_idx - start_idx) {

    case CS_TRIANGLE_CASE: /* Triangle: one-shot computation */
      {
        cs_lnum_t  v1, v2, v3;

        cs_connect_get_next_3_vertices(f2e->ids, e2v->ids, start_idx,
                                       &v1, &v2, &v3);
        compute_integral(time_eval,
                         xv + 3*v1, xv + 3*v2, xv + 3*v3, pfq.meas,
                         ana, input, val_i);
      }
      break;

    default:
      for (cs_lnum_t j = start_idx; j < end_idx; j++) {

        const cs_lnum_t  *_v = e2v->ids + 2*f2e->ids[j];
        const cs_real_t  *xv1 = xv + 3*_v[0], *xv2 = xv + 3*_v[1];

        compute_integral(time_eval,
                         xv1, xv2, pfq.center,
                         cs_math_surftri(xv1, xv2, pfq.center),
                         ana, input, val_i);

      } /* Loop on edges */
      break;

    } /* End of switch */

    /* Average */
    val_i[0] /= pfq.meas;

  } /* Loop on faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the average at each primal faces for a vector potential
 *         defined by an analytical function on a selection of (primal) cells
 *
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in]      ana                 pointer to the analytic function
 * \param[in]      input               NULL or pointer cast on-the-fly
 * \param[in]      n_loc_elts          number of elements to consider
 * \param[in]      elt_ids             pointer to the list of selected ids
 * \param[in]      compute_integral    function pointer
 * \param[in, out] values              pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_pfva_by_analytic(cs_real_t                       time_eval,
                  cs_analytic_func_t             *ana,
                  void                           *input,
                  const cs_lnum_t                 n_elts,
                  const cs_lnum_t                *elt_ids,
                  cs_quadrature_tria_integral_t  *compute_integral,
                  cs_real_t                       values[])
{
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_adjacency_t  *f2e = cs_cdo_connect->f2e;
  const cs_adjacency_t  *e2v = cs_cdo_connect->e2v;
  const cs_real_t  *xv = quant->vtx_coord;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t f = 0; f < n_elts; f++) {

    const cs_lnum_t  f_id = (elt_ids == NULL) ? f : elt_ids[f];
    const cs_quant_t pfq = cs_quant_set_face(f_id, quant);
    const cs_lnum_t  start_idx = f2e->idx[f_id], end_idx = f2e->idx[f_id+1];
    cs_real_t *val_i = values + 3*f_id;

    switch (end_idx - start_idx) {

    case CS_TRIANGLE_CASE: /* Triangle: one-shot computation */
      {
        cs_lnum_t  v1, v2, v3;

        cs_connect_get_next_3_vertices(f2e->ids, e2v->ids, start_idx,
                                       &v1, &v2, &v3);
        compute_integral(time_eval,
                         xv + 3*v1, xv + 3*v2, xv + 3*v3, pfq.meas,
                         ana, input, val_i);
      }
      break;

    default:
      for (cs_lnum_t j = start_idx; j < end_idx; j++) {

        const cs_lnum_t  *_v = e2v->ids + 2*f2e->ids[j];
        const cs_real_t  *xv1 = xv + 3*_v[0], *xv2 = xv + 3*_v[1];

        compute_integral(time_eval,
                         xv1, xv2, pfq.center,
                         cs_math_surftri(xv1, xv2, pfq.center),
                         ana, input, val_i);

      } /* Loop on edges */
      break;

    } /* End of switch */

    /* Average */
    const double _oversurf = 1./pfq.meas;
    for (short int xyz = 0; xyz < 3; xyz++)
      val_i[xyz] *= _oversurf;

  } /* Loop on faces */

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                const cs_cdo_connect_t       *connect)
{
  /* Assign static const pointers */
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute reduced quantities for an array of size equal to dim * n_x
 *         The quantities computed are synchronized in parallel.
 *
 * \param[in]      dim     local array dimension (max: 3)
 * \param[in]      n_x     number of elements
 * \param[in]      array   array to analyze
 * \param[in]      w_x     weight to apply (may be set to  NULL)
 * \param[in, out] min     resulting min array (size: dim, or 4 if dim = 3)
 * \param[in, out] max     resulting max array (size: dim, or 4 if dim = 3)
 * \param[in, out] wsum    (weighted) sum array (size: dim, or 4 if dim = 3)
 * \param[in, out] asum    (weighted) sum of absolute values (same size as wsum)
 * \param[in, out] ssum    (weighted) sum of squared values (same size as wsum)
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_array_reduction(int                     dim,
                            cs_lnum_t               n_x,
                            const cs_real_t        *array,
                            const cs_real_t        *w_x,
                            cs_real_t              *min,
                            cs_real_t              *max,
                            cs_real_t              *wsum,
                            cs_real_t              *asum,
                            cs_real_t              *ssum)
{
  assert(cs_cdo_quant != NULL && cs_cdo_connect != NULL);

  /* Get reduced quantities for this array and for this MPI rank */
  cs_real_t  dummy[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  cs_array_reduce_simple_norms_l(n_x,
                                 dim,
                                 NULL, NULL, /* elt_list, weight_list */
                                 array,
                                 w_x,
                                 min, max,
                                 dummy,     /* Not useful in this context */
                                 wsum,
                                 asum,
                                 dummy + 4, /* Not useful in this context */
                                 ssum);

  /* Parallel treatment */
  _synchronize_reduction(dim, min, max, wsum, asum, ssum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute reduced quantities for an array attached to either vertex,
 *         face or edge DoFs
 *         The weight to apply to each entity x is scanned using the adjacency
 *         structure. array size is equal to dim * n_x
 *         The quantities computed are synchronized in parallel.
 *
 * \param[in]      dim     local array dimension (max: 3)
 * \param[in]      n_x     number of elements
 * \param[in]      array   array to analyze
 * \param[in]      w_x     weight to apply (may be set to  NULL)
 * \param[in, out] min     resulting min array (size: dim, or 4 if dim = 3)
 * \param[in, out] max     resulting max array (size: dim, or 4 if dim = 3)
 * \param[in, out] wsum    (weighted) sum array (size: dim, or 4 if dim = 3)
 * \param[in, out] asum    (weighted) sum of absolute values (same size as wsum)
 * \param[in, out] ssum    (weighted) sum of squared values (same size as wsum)
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_scatter_array_reduction(int                     dim,
                                    cs_lnum_t               n_x,
                                    const cs_real_t        *array,
                                    const cs_adjacency_t   *c2x,
                                    const cs_real_t        *w_x,
                                    cs_real_t              *min,
                                    cs_real_t              *max,
                                    cs_real_t              *wsum,
                                    cs_real_t              *asum,
                                    cs_real_t              *ssum)
{
  assert(cs_cdo_quant != NULL && cs_cdo_connect != NULL);

  if (c2x == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " %s: One needs an adjacency.\n", __func__);

  /* Get the min/max for this MPI rank */
  cs_array_reduce_minmax_l(n_x, dim, NULL, array, min, max);

  /* Get reduced quantities for this array and for this MPI rank */
  cs_array_scatter_reduce_norms_l(cs_cdo_quant->n_cells,
                                  c2x->idx, c2x->ids,
                                  NULL, /* filter list */
                                  dim,
                                  n_x, array, w_x,
                                  wsum, asum, ssum); /* results */

  /* Parallel treatment */
  _synchronize_reduction(dim, min, max, wsum, asum, ssum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value related to each DoF in the case of a density field
 *         The value defined by the analytic function is by unity of volume
 *
 * \param[in]      dof_flag    indicate where the evaluation has to be done
 * \param[in]      def         pointer to a cs_xdef_t structure
 * \param[in]      time_eval   physical time at which one evaluates the term
 * \param[in, out] retval      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_density_by_analytic(cs_flag_t           dof_flag,
                                const cs_xdef_t    *def,
                                cs_real_t           time_eval,
                                cs_real_t           retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);

    /* Retrieve information from mesh location structures */
  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_lnum_t  *elt_ids
    = (cs_cdo_quant->n_cells == z->n_elts) ? NULL : z->elt_ids;

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)def->input;
  cs_quadrature_tetra_integral_t *qfunc
    = cs_quadrature_get_tetra_integral(def->dim, def->qtype);

  /* Perform the evaluation */
  if (dof_flag & CS_FLAG_SCALAR) { /* DoF is scalar-valued */

    if (cs_flag_test(dof_flag, cs_flag_primal_cell))
      _pcsd_by_analytic(time_eval, anai->func, anai->input,
                        z->n_elts, elt_ids, qfunc,
                        retval);
    else if (cs_flag_test(dof_flag, cs_flag_dual_cell))
      _dcsd_by_analytic(time_eval, anai->func, anai->input,
                        z->n_elts, elt_ids, qfunc,
                        retval);
    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);

  }
  else if (dof_flag & CS_FLAG_VECTOR) { /* DoF is vector-valued */

    if (cs_flag_test(dof_flag, cs_flag_primal_cell))
      _pcvd_by_analytic(time_eval, anai->func, anai->input,
                        z->n_elts, elt_ids, qfunc,
                        retval);
    else if (cs_flag_test(dof_flag, cs_flag_dual_cell))
      _dcvd_by_analytic(time_eval, anai->func, anai->input,
                        z->n_elts, elt_ids, qfunc,
                        retval);
    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);

  }
  else
    bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity defined by a value in the case of a density
 *         field for all the degrees of freedom
 *         Accessor to the value is by unit of volume
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      def       pointer to a cs_xdef_t structure
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_density_by_value(cs_flag_t          dof_flag,
                             const cs_xdef_t   *def,
                             cs_real_t          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);

  /* Retrieve information from mesh location structures */
  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

  /* Perform the evaluation */
  if (dof_flag & CS_FLAG_SCALAR) { /* DoF is scalar-valued */

    const cs_real_t  *constant_val = (const cs_real_t *)def->input;

    if (cs_flag_test(dof_flag, cs_flag_primal_cell))
      _pcsd_by_value(constant_val[0], z->n_elts, z->elt_ids, retval);
    else if (cs_flag_test(dof_flag, cs_flag_dual_cell))
      _dcsd_by_value(constant_val[0], z->n_elts, z->elt_ids, retval);
    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);

  }
  else if (dof_flag & CS_FLAG_VECTOR) { /* DoF is vector-valued */

    const cs_real_t  *constant_vec = (const cs_real_t *)def->input;

    if (cs_flag_test(dof_flag, cs_flag_primal_cell))
      _pcvd_by_value(constant_vec, z->n_elts, z->elt_ids, retval);
    else if (cs_flag_test(dof_flag, cs_flag_dual_cell))
      _dcvd_by_value(constant_vec, z->n_elts, z->elt_ids, retval);
    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);

  }
  else
    bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity attached to a potential field at vertices
 *         when the definition relies on an analytic expression
 *
 * \param[in]      def           pointer to a cs_xdef_t pointer
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      n_v_selected  number of selected vertices
 * \param[in]      selected_lst  list of selected vertices
 * \param[in, out] retval        pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_vertices_by_analytic(const cs_xdef_t   *def,
                                              const cs_real_t    time_eval,
                                              const cs_lnum_t    n_v_selected,
                                              const cs_lnum_t   *selected_lst,
                                              cs_real_t          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)def->input;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;

  /* Perform the evaluation */
  if (n_vertices == n_v_selected)
    anai->func(time_eval,
               n_vertices, NULL, quant->vtx_coord,
               false,  /* compacted output ? */
               anai->input,
               retval);
  else
    anai->func(time_eval,
               n_v_selected, selected_lst, quant->vtx_coord,
               false,  /* compacted output ? */
               anai->input,
               retval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity attached to a potential field at face centers
 *         when the definition relies on an analytic expression
 *
 * \param[in]      def           pointer to a cs_xdef_t pointer
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in]      n_f_selected  number of selected faces
 * \param[in]      selected_lst  list of selected faces
 * \param[in, out] retval        pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_faces_by_analytic(const cs_xdef_t   *def,
                                           const cs_real_t    time_eval,
                                           const cs_lnum_t    n_f_selected,
                                           const cs_lnum_t   *selected_lst,
                                           cs_real_t          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)def->input;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  /* Perform the evaluation */
  if (n_faces == n_f_selected) {

    /* All the support entities are selected:
       - First pass: interior faces
       - Second pass: border faces
    */
    anai->func(time_eval,
               quant->n_i_faces, NULL, quant->i_face_center,
               true, /* Output is compacted ? */
               anai->input,
               retval);
    anai->func(time_eval,
               quant->n_b_faces, NULL, quant->b_face_center,
               true, /* Output is compacted ? */
               anai->input,
               retval + def->dim*quant->n_i_faces);

  }
  else {

    assert(selected_lst != NULL);

    /* Selected faces are stored by increasing number */
    cs_lnum_t  n_i_faces = 0;
    for (cs_lnum_t  i = 0; i < n_f_selected; i++) {
      if (selected_lst[i] < quant->n_i_faces)
        n_i_faces++;
      else
        break;
    }

    /* Interior faces */
    anai->func(time_eval,
               n_i_faces, selected_lst, quant->i_face_center,
               false, /* Output is compacted ? */
               anai->input,
               retval);

    /* Border faces */
    cs_lnum_t n_b_faces = n_f_selected - n_i_faces;
    assert(n_b_faces > -1);
    anai->func(time_eval,
               n_b_faces, selected_lst + n_i_faces, quant->b_face_center,
               false, /* Output is compacted ? */
               anai->input,
               retval);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity attached to a potential field at cell centers
 *         when the definition relies on an analytic expression
 *
 * \param[in]      def         pointer to a cs_xdef_t pointer
 * \param[in]      time_eval   physical time at which one evaluates the term
 * \param[in, out] retval      pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_cells_by_analytic(const cs_xdef_t    *def,
                                           const cs_real_t     time_eval,
                                           cs_real_t           retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  cs_xdef_analytic_input_t  *anai = (cs_xdef_analytic_input_t *)def->input;

  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  if (def->meta & CS_FLAG_FULL_LOC) /* All cells are selected */
    anai->func(time_eval,
               quant->n_cells, NULL, quant->cell_centers,
               false,  /* compacted output */
               anai->input,
               retval);
  else
    anai->func(time_eval,
               z->n_elts, z->elt_ids, quant->cell_centers,
               false,  /* compacted output */
               anai->input,
               retval);

  /* No sync since theses values are computed by only one rank */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a value to each DoF in the case of a potential field in order
 *         to put a given quantity inside the volume associated to the zone
 *         related to the given definition
 *         wvals may be NULL.
 *
 * \param[in]      dof_flag  indicate where the evaluation has to be done
 * \param[in]      def       pointer to a cs_xdef_t pointer
 * \param[in, out] vvals     pointer to the first array of computed values
 * \param[in, out] wvals     pointer to the second array of computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_by_qov(cs_flag_t          dof_flag,
                             const cs_xdef_t   *def,
                             cs_real_t          vvals[],
                             cs_real_t          wvals[])
{
  /* Sanity check */
  if (vvals == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);
  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);

  const cs_real_t  *input = (cs_real_t *)def->input;
  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

  /* Perform the evaluation */
  bool check = false;
  if (dof_flag & CS_FLAG_SCALAR) { /* DoF is scalar-valued */

    const cs_real_t  const_val = input[0];

    if (cs_flag_test(dof_flag, cs_flag_primal_vtx | cs_flag_primal_cell)) {
      if (wvals == NULL)
        bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);
      _pvcsp_by_qov(const_val, z->n_elts, z->elt_ids, vvals, wvals);
      check = true;
    }
    else if (cs_flag_test(dof_flag, cs_flag_primal_vtx)) {
      _pvsp_by_qov(const_val, z->n_elts, z->elt_ids, vvals);
      check = true;
    }

  } /* Located at primal vertices */

  if (!check)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Stop evaluating a potential from 'quantity over volume'."
                "\n This situation is not handled yet."), __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a potential field at vertices from a definition by a
 *         constant value
 *
 * \param[in]      def             pointer to a cs_xdef_t pointer
 * \param[in]      n_v_selected    number of selected vertices
 * \param[in]      selected_lst    list of selected vertices
 * \param[in, out] retval          pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_vertices_by_value(const cs_xdef_t   *def,
                                           const cs_lnum_t    n_v_selected,
                                           const cs_lnum_t   *selected_lst,
                                           cs_real_t          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);
  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  const cs_lnum_t  n_vertices = cs_cdo_quant->n_vertices;
  const cs_real_t  *input = (cs_real_t *)def->input;

  /* Perform the evaluation */
  if (def->dim == 1) { /* DoF is scalar-valued */

    const cs_real_t  const_val = input[0];

    if (n_v_selected == n_vertices) {

#     pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
        retval[v_id] = const_val;

    }
    else {  /* Partial selection */

      assert(selected_lst != NULL);

      /* Loop on selected vertices */
      for (cs_lnum_t i = 0; i < n_vertices; i++)
        retval[selected_lst[i]] = const_val;

    }

  }
  else if (def->dim == 3) {

    const size_t _3real = 3*sizeof(cs_real_t);

    if (n_v_selected == n_vertices) {

#     pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
        memcpy(retval + 3*v_id, input, _3real);

    }
    else { /* Partial selection */

      assert(selected_lst != NULL);

      /* Loop on selected vertices */
      for (cs_lnum_t i = 0; i < n_vertices; i++)
        memcpy(retval + 3*selected_lst[i], input, _3real);

    }

  }
  else
    bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a potential field atface centers from a definition by a
 *         constant value
 *
 * \param[in]      def             pointer to a cs_xdef_t pointer
 * \param[in]      n_f_selected    number of selected faces
 * \param[in]      selected_lst    list of selected faces
 * \param[in, out] retval          pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_faces_by_value(const cs_xdef_t   *def,
                                        const cs_lnum_t    n_f_selected,
                                        const cs_lnum_t   *selected_lst,
                                        cs_real_t          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);
  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  const cs_lnum_t  n_faces = cs_cdo_quant->n_faces;
  const cs_real_t  *input = (cs_real_t *)def->input;

  if (def->dim == 1) { /* DoF is scalar-valued */

    if (n_faces == n_f_selected) {

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++)
        retval[f_id] = input[0];

    }
    else {

      assert(selected_lst != NULL);

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f = 0; f < n_f_selected; f++)
        retval[selected_lst[f]] = input[0];

    }

  }
  else if (def->dim == 3) {

    const size_t _3real = 3*sizeof(cs_real_t);

    if (n_faces == n_f_selected) {

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++)
        memcpy(retval + 3*f_id, input, _3real);

    }
    else {

      assert(selected_lst != NULL);

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f = 0; f < n_f_selected; f++)
        memcpy(retval + 3*selected_lst[f], input, _3real);

    }

  }
  else {

    const size_t s = def->dim*sizeof(cs_real_t);

    if (n_faces == n_f_selected) {

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++)
        memcpy(retval + def->dim*f_id, input, s);

    }
    else {

      assert(selected_lst != NULL);

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f = 0; f < n_f_selected; f++)
        memcpy(retval + def->dim*selected_lst[f], input, s);

    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a potential field at cell centers from a definition by
 *         value
 *
 * \param[in]      def       pointer to a cs_xdef_t pointer
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_cells_by_value(const cs_xdef_t   *def,
                                        cs_real_t          retval[])
{
  /* Sanity check */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);
  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;
  const cs_real_t  *input = (cs_real_t *)def->input;
  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

  if (def->dim == 1) { /* DoF is scalar-valued */

    const cs_real_t  const_val = input[0];

    if (n_cells == z->n_elts) {
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        retval[c_id] = const_val;

    }
    else { /* Partial selection */

#     pragma omp parallel for if (z->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < z->n_elts; i++)
        retval[z->elt_ids[i]] = const_val;

    }

  }
  else if (def->dim == 3) {

    const size_t _3real = 3*sizeof(cs_real_t);

    if (n_cells == z->n_elts) {
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        memcpy(retval + 3*c_id, input, _3real);

    }
    else { /* Partial selection */

#     pragma omp parallel for if (z->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < z->n_elts; i++)
        memcpy(retval + 3*z->elt_ids[i], input, _3real);

    }

  }
  else {

    const size_t s = def->dim*sizeof(cs_real_t);

    if (n_cells == z->n_elts) {
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
        memcpy(retval + def->dim*c_id, input, s);

    }
    else { /* Partial selection */

#     pragma omp parallel for if (z->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < z->n_elts; i++)
        memcpy(retval + def->dim*z->elt_ids[i], input, s);

    }

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the circulation along a selection of (primal) edges.
 *         Circulation is defined thanks to a constant vector field (by value)
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      n_e_selected   number of selected edges
 * \param[in]      selected_lst   list of selected edges
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_circulation_along_edges_by_value(const cs_xdef_t   *def,
                                             const cs_lnum_t    n_e_selected,
                                             const cs_lnum_t   *selected_lst,
                                             cs_real_t          retval[])
{
  /* Sanity checks */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->type == CS_XDEF_BY_VALUE);

  const cs_lnum_t  n_edges = cs_cdo_quant->n_edges;
  const cs_real_t  *edge_vector = cs_cdo_quant->edge_vector;
  const cs_real_t  *input = (cs_real_t *)def->input;

  /* DoF is scalar-valued since this is a circulation but the definition is
   * either scalar-valued meaning that one only gives the tangential part or
   * vector-valued (in this case, one needs to extract the tangential part) */

  switch (def->dim) {

  case 1: /* Scalar-valued integral */
    if (n_edges == n_e_selected) {

#     pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++)
        retval[e_id] = input[0];

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != NULL);

#     pragma omp parallel for if (n_e_selected > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_e_selected; e++) {
        const cs_lnum_t e_id = selected_lst[e];
        retval[e_id] = input[0];
      }

    }
    break;

  case 3:
    if (n_edges == n_e_selected) {

#     pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++)
        retval[e_id] = _dp3(input, edge_vector + 3*e_id);

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != NULL);

#     pragma omp parallel for if (n_e_selected > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_e_selected; e++) {
        const cs_lnum_t e_id = selected_lst[e];
        retval[e_id] = _dp3(input, edge_vector + 3*e_id);
      }

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid dimension value %d. Only 1 and 3 are valid.\n",
              __func__, def->dim);

  } /* End of switch on dimension */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the circulation along a selection of (primal) edges.
 *         Circulation is defined thanks to an array
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      n_e_selected   number of selected edges
 * \param[in]      selected_lst   list of selected edges
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_circulation_along_edges_by_array(const cs_xdef_t   *def,
                                             const cs_lnum_t    n_e_selected,
                                             const cs_lnum_t   *selected_lst,
                                             cs_real_t          retval[])
{
  /* Sanity checks */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->type == CS_XDEF_BY_ARRAY);

  const cs_lnum_t  n_edges = cs_cdo_quant->n_edges;
  const cs_real_t  *edge_vector = cs_cdo_quant->edge_vector;

  cs_xdef_array_input_t  *ainput = (cs_xdef_array_input_t *)def->input;
  assert(cs_flag_test(ainput->loc, cs_flag_primal_edge));

  /* DoF is scalar-valued since this is a circulation but the definition is
   * either scalar-valued meaning that one only gives the tangential part or
   * vector-valued (in this case, one needs to extract the tangential part) */

  switch (def->dim) {

  case 1: /* Scalar-valued integral */
    assert(ainput->stride == 1);
    if (n_edges == n_e_selected) {

#     pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++)
        retval[e_id] = ainput->values[e_id];

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != NULL);

#     pragma omp parallel for if (n_e_selected > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_e_selected; e++) {
        const cs_lnum_t e_id = selected_lst[e];
        retval[e_id] = ainput->values[e_id];
      }

    }
    break;

  case 3:
    assert(ainput->stride == 3);
    if (n_edges == n_e_selected) {

#     pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++)
        retval[e_id] = _dp3(ainput->values + 3*e_id, edge_vector + 3*e_id);

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != NULL);

#     pragma omp parallel for if (n_e_selected > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_e_selected; e++) {
        const cs_lnum_t e_id = selected_lst[e];
        retval[e_id] = _dp3(ainput->values + 3*e_id, edge_vector + 3*e_id);
      }

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid dimension value %d. Only 1 and 3 are valid.\n",
              __func__, def->dim);

  } /* End of switch on dimension */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the circulation along a selection of (primal) edges.
 *         Circulation is defined by an analytical function.
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      n_e_selected   number of selected edges
 * \param[in]      selected_lst   list of selected edges
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_circulation_along_edges_by_analytic(const cs_xdef_t   *def,
                                                const cs_real_t    time_eval,
                                                const cs_lnum_t    n_e_selected,
                                                const cs_lnum_t   *selected_lst,
                                                cs_real_t          retval[])
{
  /* Sanity checks */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  const cs_lnum_t  n_edges = cs_cdo_quant->n_edges;
  const cs_real_t  *edge_vector = cs_cdo_quant->edge_vector;
  const cs_real_t  *xv = cs_cdo_quant->vtx_coord;
  const cs_adjacency_t  *e2v = cs_cdo_connect->e2v;

  cs_xdef_analytic_input_t *anai = (cs_xdef_analytic_input_t *)def->input;
  cs_quadrature_edge_integral_t
    *qfunc = cs_quadrature_get_edge_integral(def->dim, def->qtype);

  /* DoF is scalar-valued since this is a circulation but the definition is
   * either scalar-valued meaning that one only gives the tangential part or
   * vector-valued (in this case, one needs to extract the tangential part) */
  switch (def->dim) {

  case 1: /* Scalar-valued integral */
    if (n_edges == n_e_selected) {

#     pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {

        const cs_lnum_t  *_v = e2v->ids + 2*e_id;

        cs_real_t  e_len = cs_math_3_norm(edge_vector + 3*e_id);
        cs_real_t  integral = 0.;
        qfunc(time_eval, xv + 3*_v[0], xv + 3*_v[1], e_len,
              anai->func, anai->input, &integral);

        retval[e_id] = integral;

      }

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != NULL);

#     pragma omp parallel for if (n_e_selected > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_e_selected; e++) {

        const cs_lnum_t e_id = selected_lst[e];
        const cs_lnum_t  *_v = e2v->ids + 2*e_id;

        cs_real_t  e_len = cs_math_3_norm(edge_vector + 3*e_id);
        cs_real_t  integral = 0.;
        qfunc(time_eval, xv + 3*_v[0], xv + 3*_v[1], e_len,
              anai->func, anai->input, &integral);

        retval[e_id] = integral;

      }

    }
    break;

  case 3: /* Vector-valued case */
    if (n_edges == n_e_selected) {

#   pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++) {

        const cs_lnum_t  *_v = e2v->ids + 2*e_id;

        cs_nvec3_t  e_vec;
        cs_nvec3(edge_vector + 3*e_id, &e_vec);

        cs_real_3_t  integral = {0., 0., 0.};
        qfunc(time_eval, xv + 3*_v[0], xv + 3*_v[1], e_vec.meas,
              anai->func, anai->input, integral);

        retval[e_id] = _dp3(integral, e_vec.unitv);

      }

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != NULL);

#     pragma omp parallel for if (n_e_selected > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_e_selected; e++) {

        const cs_lnum_t e_id = selected_lst[e];
        const cs_lnum_t  *_v = e2v->ids + 2*e_id;

        cs_nvec3_t  e_vec;
        cs_nvec3(edge_vector + 3*e_id, &e_vec);

        cs_real_3_t  integral = {0., 0., 0.};
        qfunc(time_eval, xv + 3*_v[0], xv + 3*_v[1], e_vec.meas,
              anai->func, anai->input, integral);

        retval[e_id] = _dp3(integral, e_vec.unitv);

      }

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid dimension value %d. Only 1 and 3 are valid.\n",
              __func__, def->dim);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the faces
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      n_f_selected   number of selected faces
 * \param[in]      selected_lst   list of selected faces
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_faces_by_value(const cs_xdef_t   *def,
                                      const cs_lnum_t    n_f_selected,
                                      const cs_lnum_t   *selected_lst,
                                      cs_real_t          retval[])
{
  /* Sanity checks */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  const cs_lnum_t  n_faces = cs_cdo_quant->n_faces;
  const cs_real_t  *input = (cs_real_t *)def->input;

  if (n_faces == n_f_selected) {

    if (def->dim == 1) { /* DoF is scalar-valued */

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++)
        retval[f_id] = input[0];

    }
    else { /* Multi-valued case */

      const size_t  s = def->dim*sizeof(cs_real_t);
#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f_id = 0; f_id < n_faces; f_id++)
        memcpy(retval + def->dim*f_id, input, s);

    }

  }
  else { /* Definition does not apply to all entities */

    assert(selected_lst != NULL);

    if (def->dim == 1) {

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f = 0; f < n_f_selected; f++)
        retval[selected_lst[f]] = input[0];

    }
    else { /* Multi-valued case */

      const size_t  s = def->dim*sizeof(cs_real_t);
#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f = 0; f < n_f_selected; f++)
        memcpy(retval + def->dim*selected_lst[f], input, s);

    }

  } /* Deal with a selection of cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the faces.
 *         Warning: retval has to be initialize before calling this function.
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      n_f_selected   number of selected faces
 * \param[in]      selected_lst   list of selected faces
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_faces_by_analytic(const cs_xdef_t    *def,
                                         const cs_real_t     time_eval,
                                         const cs_lnum_t     n_f_selected,
                                         const cs_lnum_t    *selected_lst,
                                         cs_real_t           retval[])
{
  /* Sanity checks */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  cs_quadrature_tria_integral_t
    *qfunc = cs_quadrature_get_tria_integral(def->dim, def->qtype);
  cs_xdef_analytic_input_t *anai = (cs_xdef_analytic_input_t *)def->input;

  switch (def->dim) {

  case 1: /* Scalar-valued */
    _pfsa_by_analytic(time_eval,
                      anai->func, anai->input,
                      n_f_selected, selected_lst,
                      qfunc,
                      retval);
    break;

  case 3: /* Vector-valued */
    _pfva_by_analytic(time_eval,
                      anai->func, anai->input,
                      n_f_selected, selected_lst,
                      qfunc,
                      retval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid dimension of analytical function.\n"), __func__);

  } /* End of switch on dimension */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the cells
 *
 * \param[in]      def       pointer to a cs_xdef_t pointer
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_cells_by_value(const cs_xdef_t   *def,
                                      cs_real_t          retval[])
{
  /* Sanity checks */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);

  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_real_t  *input = (cs_real_t *)def->input;

  switch (def->dim) {

  case 1: /* Scalar-valued */
    _pcsa_by_value(input[0], z->n_elts, z->elt_ids, retval);
    break;

  case 3: /* Vector-valued */
    _pcva_by_value(input, z->n_elts, z->elt_ids, retval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid dimension of analytical function.\n"), __func__);
    break;

  } /* End of switch on dimension */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the cells
 *
 * \param[in]      def       pointer to a cs_xdef_t pointer
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_cells_by_array(const cs_xdef_t   *def,
                                      cs_real_t          retval[])
{
  /* Sanity checks */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);

  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_xdef_array_input_t  *input = (cs_xdef_array_input_t *)def->input;
  const int  stride = input->stride;
  const cs_real_t  *val = input->values;

  if (cs_flag_test(input->loc, cs_flag_primal_cell) == false)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case. Not implemented yet.",
              __func__);

  if (def->meta & CS_FLAG_FULL_LOC)
    memcpy(retval, val, stride*sizeof(cs_real_t)*cs_cdo_quant->n_cells);

  else {

    assert(z->elt_ids != NULL);
    if (stride == 1) {

#     pragma omp parallel for if (z->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < z->n_elts; i++) {
        const cs_lnum_t  c_id = z->elt_ids[i];
        retval[c_id] = val[c_id];
      }

    }
    else {

#     pragma omp parallel for if (z->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < z->n_elts; i++) {
        const cs_lnum_t  c_id = z->elt_ids[i];
        memcpy(retval + stride*c_id, val + stride*c_id,
               stride*sizeof(cs_real_t));
      }

    }

  } /* deal with a selection of cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the average of a function on the cells.
 *         Warning: retval has to be initialize before calling this function.
 *
 * \param[in]      def        pointer to a cs_xdef_t pointer
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] retval     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_cells_by_analytic(const cs_xdef_t   *def,
                                         cs_real_t          time_eval,
                                         cs_real_t          retval[])
{
  /* Sanity checks */
  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def != NULL);
  assert(def->support == CS_XDEF_SUPPORT_VOLUME);

  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_lnum_t  *elt_ids
    = (cs_cdo_quant->n_cells == z->n_elts) ? NULL : z->elt_ids;

  cs_quadrature_tetra_integral_t
    *qfunc = cs_quadrature_get_tetra_integral(def->dim, def->qtype);
  cs_xdef_analytic_input_t *anai = (cs_xdef_analytic_input_t *)def->input;

  switch (def->dim) {

  case 1: /* Scalar-valued */
    if (elt_ids == NULL)
      memset(retval, 0, cs_cdo_quant->n_cells*sizeof(cs_real_t));
    else {
      for (cs_lnum_t i = 0; i < z->n_elts; i++)
        retval[z->elt_ids[i]] = 0;
    }

    _pcsa_by_analytic(time_eval,
                      anai->func, anai->input, z->n_elts, elt_ids, qfunc,
                      retval);
    break;

  case 3: /* Vector-valued */
    if (elt_ids == NULL)
      memset(retval, 0, 3*cs_cdo_quant->n_cells*sizeof(cs_real_t));
    else {
      for (cs_lnum_t i = 0; i < z->n_elts; i++)
        for (int k = 0; k < 3; k++)
          retval[3*z->elt_ids[i]+k] = 0;
    }

    _pcva_by_analytic(time_eval,
                      anai->func, anai->input, z->n_elts, elt_ids, qfunc,
                      retval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid dimension of analytical function.\n"), __func__);
    break;

  } /* End of switch on the dimension */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the integral over the full computational domain of a
 *         quantity defined by an array
 *
 * \param[in]      array_loc  flag indicating where are located values
 * \param[in]      array_val  array of values
 *
 * \return the value of the integration
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_evaluate_scal_domain_integral_by_array(cs_flag_t         array_loc,
                                          const cs_real_t  *array_val)
{
  cs_real_t  result = 0.;

  /* Sanity checks */
  if (array_val == NULL)
    return result;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  if (cs_flag_test(array_loc, cs_flag_primal_cell)) {

#   pragma omp parallel for reduction(+:result)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      result += array_val[c_id] * quant->cell_vol[c_id];

  }
  else if (cs_flag_test(array_loc, cs_flag_primal_vtx)) {

    const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
    const cs_real_t  *dc_vol = quant->dcell_vol;

#   pragma omp parallel for reduction(+:result)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        result += dc_vol[j] * array_val[c2v->ids[j]];

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid array location. Stop evaluation.", __func__);

  if (cs_glob_n_ranks > 1)      /* MPI synchronization */
    cs_parall_sum(1, CS_REAL_TYPE, &result);

  return result;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
