/*============================================================================
 * Functions and structures to deal with evaluation of quantities
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

#include "cs_array.h"
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

/* Block size for superblock algorithm */

#define CS_SBLOCK_BLOCK_SIZE 60

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/*=============================================================================
 * Local static variables
 *============================================================================*/

/* Pointer to shared structures (owned by a cs_domain_t structure) */

static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;
static const cs_mesh_t  *cs_shared_mesh;

static const char _err_empty_array[] =
  " %s: Array storing the evaluation should be allocated before the call"
  " to this function.";
static const char _err_not_handled[] = " %s: Case not handled yet.";

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Sanity checkings before computing norms
 *
 * \param[in]  func_name   name of the calling function
 * \param[in]  c2x         first pointer to check
 * \param[in]  w_c2x       second pointer to check
 */
/*----------------------------------------------------------------------------*/

static inline void
_sanity_checks(const char              func_name[],
               const cs_adjacency_t   *c2x,
               const cs_real_t        *w_c2x)
{
  assert(cs_cdo_quant != nullptr && cs_cdo_connect != nullptr);

  if (c2x == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              " %s: The cs_adjacency_t structure is not allocated.\n",
              func_name);

  if (w_c2x == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              " %s: The array storing weights is not allocated.\n",
              func_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *        When called inside an OpenMP parallel section, this will return the
 *        start an past-the-end indexes for the array range assigned to that
 *        thread. In other cases, the start index is 1, and the past-the-end
 *        index is n;
 *
 * \param[in]       n     size of array
 * \param[in, out]  s_id  start index for the current thread
 * \param[in, out]  e_id  past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

static inline void
_thread_range(cs_lnum_t    n,
              cs_lnum_t   *s_id,
              cs_lnum_t   *e_id)
{
#if defined(HAVE_OPENMP)
  const int t_id = omp_get_thread_num();
  const int n_t = omp_get_num_threads();
  const cs_lnum_t t_n = (n + n_t - 1) / n_t;

  *s_id =  t_id    * t_n;
  *e_id = (t_id+1) * t_n;
  *s_id = cs_align(*s_id, CS_CL);
  *e_id = cs_align(*e_id, CS_CL);
  if (*e_id > n) *e_id = n;
#else
  *s_id = 0;
  *e_id = n;
#endif
}

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

  if (dim == 1) {

    /* Min/Max */

    cs_real_t  minmax[2] = {-min[0], max[0]};
    cs_parall_max(2, CS_REAL_TYPE, minmax);

    min[0] = -minmax[0];
    max[0] = minmax[1];

    /* Sums */

    cs_real_t  sums[3] = {wsum[0], asum[0], ssum[0]};
    cs_parall_sum(3, CS_REAL_TYPE, sums);

    wsum[0] = sums[0];
    asum[0] = sums[1];
    ssum[0] = sums[2];

  }
  else if (dim == 3) {

    /* Min/Max */

    cs_real_t  minmax[8];
    for (int i = 0; i < 4; i++)
      minmax[i] = -min[i], minmax[4+i] = max[i];

    cs_parall_max(8, CS_REAL_TYPE, minmax);

    for (int i = 0; i < 4; i++)
      min[i] = -minmax[i], max[i] = minmax[4+i];

    cs_real_t  sums[12];
    for (int i = 0; i < 4; i++)
      sums[i] = wsum[i], sums[4+i] = asum[i], sums[8+i] = ssum[i];

    cs_parall_sum(12, CS_REAL_TYPE, sums);

    for (int i = 0; i < 4; i++)
      wsum[i] = sums[i], asum[i] = sums[4+i], ssum[i] = sums[8+i];

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid dimension (=%d). Expected 1 or 3.\n",
              __func__, dim);
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
  const cs_mesh_t  *m = cs_shared_mesh;
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
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
  const cs_mesh_t  *mesh = cs_shared_mesh;
  const cs_lnum_t  n_cells_with_ghosts = mesh->n_cells_with_ghosts;

  cs_lnum_t *c_tags = nullptr;
  BFT_MALLOC(c_tags, n_cells_with_ghosts, cs_lnum_t);

  if (n_elts < n_cells) { /* Only some cells are selected */

    cs_array_lnum_fill_zero(n_vertices, v_tags);
    cs_array_lnum_fill_zero(n_cells_with_ghosts, c_tags);

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

    cs_array_lnum_set_value(n_vertices, -1, v_tags);
    cs_array_lnum_set_value(n_cells, 1, c_tags);
    cs_array_lnum_fill_zero(mesh->n_ghost_cells, c_tags + n_cells);

  }

  cs_halo_sync_num(mesh->halo, CS_HALO_STANDARD, c_tags);

  /* Second pass: detect cells at the frontier of the selection */

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    const cs_lnum_t  c_id = (n_elts == n_cells) ? i : elt_ids[i];
    _untag_frontier_vertices(c_id, c_tags, v_tags);
  }

  /* Not needed anymore */

  BFT_FREE(c_tags);

  /* Handle parallelism (always the scalar interface) */

  if (cs_cdo_connect->vtx_ifs != nullptr)
    cs_interface_set_max(cs_cdo_connect->vtx_ifs,
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
  const cs_real_t  *dc_vol = quant->pvol_vc;
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;

  cs_lnum_t *v_tags = nullptr;

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

    cs_array_real_set_scalar(n_vertices, val_to_set, v_vals);

  }

  BFT_FREE(v_tags);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a value to each DoF such that a given quantity is put inside
 *         the volume associated to the list of cells
 *         Case of primal cells for scalar-valued quantities.
 *
 * \param[in]      quantity_val  amount of quantity to distribute
 * \param[in]      n_elts        number of cells to consider
 * \param[in]      elt_ids       pointer to the list of selected cell ids
 * \param[in, out] c_vals        pointer to the array storing the values
 */
/*----------------------------------------------------------------------------*/

static void
_pcsp_by_qov(const cs_real_t    quantity_val,
             cs_lnum_t          n_elts,
             const cs_lnum_t   *elt_ids,
             cs_real_t          c_vals[])
{
  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;
  const cs_lnum_t  n_cells = cdoq->n_cells;

  double  volume_marked = 0.;

  if (n_elts == n_cells) {

#   pragma omp parallel for reduction(+:volume_marked) if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      volume_marked += cdoq->cell_vol[i];

  }
  else {

    assert(elt_ids != nullptr);
#   pragma omp parallel for reduction(+:volume_marked) if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++)
      volume_marked += cdoq->cell_vol[elt_ids[i]];

  }

  /* Handle parallelism */

  cs_parall_sum(1, CS_DOUBLE, &volume_marked);

  cs_real_t  val_to_set = quantity_val;
  if (volume_marked > 0)
    val_to_set /= volume_marked;

  if (n_elts == n_cells)
    cs_array_real_set_scalar(n_cells, val_to_set, c_vals);
  else
    cs_array_real_set_scalar_on_subset(n_elts, elt_ids, val_to_set, c_vals);
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
  const cs_real_t  *dc_vol = quant->pvol_vc;
  const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;

  cs_lnum_t *v_tags = nullptr;

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

    assert(elt_ids != nullptr);

    cs_array_real_set_scalar_on_subset(n_elts, elt_ids, val_to_set, c_vals);

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
      if (v_tags[v_id] == -1)
        v_vals[v_id] = val_to_set;

  }
  else { /* All cells are selected */

    cs_array_real_set_scalar(n_cells, val_to_set, c_vals);
    cs_array_real_set_scalar(n_vertices, val_to_set, v_vals);

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
 * \param[in]      input             nullptr or pointer cast on-the-fly
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

    const cs_lnum_t   c_id = (elt_ids == nullptr) ? id : elt_ids[id];
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

        compute_integral(time_eval, xv1, xe, xf, xc, quant->pvol_vc[v1],
                         ana, input, values + v1);
        compute_integral(time_eval, xv2, xe, xf, xc, quant->pvol_vc[v2],
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
 * \param[in]      input             nullptr or pointer cast on-the-fly
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

    const cs_lnum_t   c_id = (elt_ids == nullptr) ? id : elt_ids[id];
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

        compute_integral(time_eval, xv1, xe, xf, xc, quant->pvol_vc[v1],
                         ana, input, values + dim*v1);
        compute_integral(time_eval, xv2, xe, xf, xc, quant->pvol_vc[v2],
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
 * \param[in]      input             nullptr or pointer cast on-the-fly
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
  constexpr cs_real_t c_1ov3 = 1./3.;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *xv = quant->vtx_coord;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t id = 0; id < n_elts; id++) {

    const cs_lnum_t c_id = (elt_ids == nullptr) ? id : elt_ids[id];
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
          c_1ov3 * cs_math_3_dot_product(pfq.unitv,
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
 * \param[in]      input             nullptr or pointer cast on-the-fly
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
  constexpr cs_real_t c_1ov3 = 1./3.;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *xv = quant->vtx_coord;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;
  const int  dim = 3;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t  id = 0; id < n_elts; id++) {

    const cs_lnum_t c_id = (elt_ids == nullptr) ? id : elt_ids[id];
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
        const double hfc = c_1ov3 *
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
 * \param[in]      input             nullptr or pointer cast on-the-fly
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
  constexpr cs_real_t c_1ov3 = 1./3.;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *xv = quant->vtx_coord;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t id = 0; id < n_elts; id++) {

    const cs_lnum_t c_id = (elt_ids == nullptr) ? id : elt_ids[id];
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
          c_1ov3 * cs_math_3_dot_product(pfq.unitv,
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
 * \param[in]      input             nullptr or pointer cast on-the-fly
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
  constexpr cs_real_t c_1ov3 = 1./3.;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_real_t  *xv = quant->vtx_coord;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *c2f = connect->c2f;
  const cs_adjacency_t  *f2e = connect->f2e;

# pragma omp parallel for if (n_elts > CS_THR_MIN)
  for (cs_lnum_t id = 0; id < n_elts; id++) {

    const cs_lnum_t c_id  = (elt_ids == nullptr) ? id : elt_ids[id];
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
          c_1ov3 * cs_math_3_dot_product(pfq.unitv,
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
  const cs_real_t  *dual_vol = quant->pvol_vc; /* scan by c2v */

  if (elt_ids == nullptr) {

    assert(n_elts == quant->n_cells);
    for (cs_lnum_t c_id = 0; c_id < n_elts; c_id++)
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id + 1]; j++)
        values[c2v->ids[j]] += dual_vol[j] * const_val;
  }
  else { /* Loop on selected cells */

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t c_id = elt_ids[i];
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id + 1]; j++)
        values[c2v->ids[j]] += dual_vol[j] * const_val;
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
_dcvd_by_value(const cs_real_t  const_vec[3],
               const cs_lnum_t  n_elts,
               const cs_lnum_t *elt_ids,
               cs_real_t        values[])
{
  const cs_adjacency_t *c2v      = cs_cdo_connect->c2v;
  const cs_real_t      *dual_vol = cs_cdo_quant->pvol_vc; /* scan by c2v */

  if (elt_ids == nullptr) {

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
 * \brief  Get the average at each primal faces for a scalar potential
 *         defined by an analytical function on a selection of (primal) cells
 *
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      ana               pointer to the analytic function
 * \param[in]      input             nullptr or pointer cast on-the-fly
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

    const cs_lnum_t  f_id      = (elt_ids == nullptr) ? f : elt_ids[f];
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
 * \param[in]      input               nullptr or pointer cast on-the-fly
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

    const cs_lnum_t  f_id      = (elt_ids == nullptr) ? f : elt_ids[f];
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set an array by a constant definition by value
 *
 * \param[in]      n_elts     number of elements
 * \param[in]      stride     number of values for each element
 * \param[in]      elt_ids    nullptr or list of element ids
 * \param[in]      ref_val    value(s) to set
 * \param[in, out] retval     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

static void
_eval_by_value(cs_lnum_t          n_elts,
               int                stride,
               const cs_lnum_t   *elt_ids,
               const cs_real_t    ref_val[],
               cs_real_t          retval[])
{
  switch (stride) {

  case 1:  /* DoF is scalar-valued */
    cs_array_real_set_scalar_on_subset(n_elts, elt_ids, ref_val[0], retval);
    break;

  case 3: /* DoF is vector-valued */
    cs_array_real_set_vector_on_subset(n_elts, elt_ids, ref_val, retval);
    break;

  case 9: /* DoF is tensor-valued */
    {
      const cs_real_t  tensor[3][3] = {{ref_val[0], ref_val[1], ref_val[2]},
                                       {ref_val[3], ref_val[4], ref_val[5]},
                                       {ref_val[6], ref_val[7], ref_val[8]}};

      cs_array_real_set_tensor_on_subset(n_elts, elt_ids, tensor, retval);
    }
    break;

  default:
    cs_array_real_set_value_on_subset(n_elts, stride, elt_ids, ref_val, retval);
    break;
  }
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set shared pointers to main domain members
 *
 * \param[in] quant    pointer to additional mesh quantities for CDO schemes
 * \param[in] connect  pointer to additional mesh connectivities for CDO schemes
 * \param[in] mesh     pointer to the shared mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_init_sharing(const cs_cdo_quantities_t    *quant,
                         const cs_cdo_connect_t       *connect,
                         const cs_mesh_t              *mesh)
{
  /* Assign static const pointers */

  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
  cs_shared_mesh = mesh;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute reduced quantities for an array of size equal to dim * n_x
 *         The computed quantities are synchronized in parallel.
 *
 * \param[in]      dim     local array dimension (max: 3)
 * \param[in]      n_x     number of elements
 * \param[in]      array   array to analyze
 * \param[in]      w_x     weight to apply (may be set to  nullptr)
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
  if (array == nullptr)
    return;

  assert(cs_cdo_quant != nullptr && cs_cdo_connect != nullptr);

  /* Get reduced quantities for this array and for this MPI rank */

  cs_real_t  dummy[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  cs_array_reduce_simple_norms_l(n_x,
                                 dim,
                                 nullptr,
                                 nullptr, /* elt_list, weight_list */
                                 array,
                                 w_x,
                                 min,
                                 max,
                                 dummy, /* Not useful in this context */
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
 *         The computed quantities are synchronized in parallel.
 *
 * \param[in]      dim     local array dimension (max: 3)
 * \param[in]      n_x     number of elements
 * \param[in]      array   array to analyze
 * \param[in]      c2x     pointer to the associated cs_adjacency_t structure
 * \param[in]      w_x     weight to apply (may be set to  nullptr)
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
  if (array == nullptr)
    return;
  if (c2x == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              " %s: One needs an adjacency.\n", __func__);

  assert(cs_cdo_quant != nullptr && cs_cdo_connect != nullptr);

  /* Get the min/max for this MPI rank */

  cs_array_reduce_minmax_l(n_x, dim, nullptr, array, min, max);

  /* Get reduced quantities for this array and for this MPI rank */

  cs_array_scatter_reduce_norms_l(cs_cdo_quant->n_cells,
                                  c2x->idx,
                                  c2x->ids,
                                  nullptr, /* filter list */
                                  dim,
                                  n_x,
                                  array,
                                  w_x,
                                  wsum,
                                  asum,
                                  ssum); /* results */

  /* Parallel treatment */

  _synchronize_reduction(dim, min, max, wsum, asum, ssum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity defined by a value in the case of a density
 *         field for all the degrees of freedom
 *         Accessor to the value is by unit of volume and the return values are
 *         integrated over a volume
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  /* Retrieve information from mesh location structures */

  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;
  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_real_t  *ref_val = (const cs_real_t *)def->context;
  const cs_lnum_t            *elt_ids
    = (cdoq->n_cells == z->n_elts) ? nullptr : z->elt_ids;

  /* Perform the evaluation */

  if (cs_flag_test(dof_flag, cs_flag_primal_cell)) {

    if (dof_flag & CS_FLAG_SCALAR)
      cs_array_real_set_wscalar_on_subset(z->n_elts, elt_ids, ref_val[0],
                                          cdoq->cell_vol, /* weights */
                                          retval);
    else if (dof_flag & CS_FLAG_VECTOR)
      cs_array_real_set_wvector_on_subset(z->n_elts, elt_ids, ref_val,
                                          cdoq->cell_vol, /* weights */
                                          retval);
    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);

  }
  else if (cs_flag_test(dof_flag, cs_flag_dual_cell)) {

    if (dof_flag & CS_FLAG_SCALAR)
      _dcsd_by_value(ref_val[0], z->n_elts, z->elt_ids, retval);
    else if (dof_flag & CS_FLAG_VECTOR)
      _dcvd_by_value(ref_val, z->n_elts, z->elt_ids, retval);
    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);

  }
  else
    bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value related to each DoF in the case of a density field
 *         The value defined by the analytic function is by unity of volume and
 *         the return values are integrated over a volume
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  /* Retrieve information from mesh location structures */

  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_lnum_t  *elt_ids
    = (cs_cdo_quant->n_cells == z->n_elts) ? nullptr : z->elt_ids;

  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)def->context;
  cs_quadrature_tetra_integral_t *qfunc
    = cs_quadrature_get_tetra_integral(def->dim, def->qtype);

  /* Perform the evaluation */

  if (dof_flag & CS_FLAG_SCALAR) { /* DoF is scalar-valued */

    if (cs_flag_test(dof_flag, cs_flag_primal_cell))
      _pcsd_by_analytic(time_eval, ac->func, ac->input,
                        z->n_elts, elt_ids, qfunc,
                        retval);
    else if (cs_flag_test(dof_flag, cs_flag_dual_cell))
      _dcsd_by_analytic(time_eval, ac->func, ac->input,
                        z->n_elts, elt_ids, qfunc,
                        retval);
    else
      bft_error(__FILE__, __LINE__, 0, _err_not_handled, __func__);

  }
  else if (dof_flag & CS_FLAG_VECTOR) { /* DoF is vector-valued */

    if (cs_flag_test(dof_flag, cs_flag_primal_cell))
      _pcvd_by_analytic(time_eval, ac->func, ac->input,
                        z->n_elts, elt_ids, qfunc,
                        retval);
    else if (cs_flag_test(dof_flag, cs_flag_dual_cell))
      _dcvd_by_analytic(time_eval, ac->func, ac->input,
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  /* Perform the evaluation */

  const cs_lnum_t  n_vertices = cs_cdo_quant->n_vertices;
  const cs_real_t  *input = (cs_real_t *)def->context;
  const cs_lnum_t  *elt_ids
    = (n_v_selected == n_vertices) ? nullptr : selected_lst;

  _eval_by_value(n_v_selected, def->dim, elt_ids, input, retval);
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)def->context;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;

  /* Perform the evaluation */

  if (n_vertices == n_v_selected)
    ac->func(time_eval,
             n_vertices,
             nullptr,
             quant->vtx_coord,
             false, /* dense output ? */
             ac->input,
             retval);
  else
    ac->func(time_eval,
             n_v_selected, selected_lst, quant->vtx_coord,
             false,  /* dense output ? */
             ac->input,
             retval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity attached to a potential field at vertices
 *         when the definition relies on a DoF function (Degrees of freedom)
 *
 * \param[in]      def           pointer to a cs_xdef_t pointer
 * \param[in]      n_v_selected  number of selected vertices
 * \param[in]      selected_lst  list of selected vertices
 * \param[in, out] retval        pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_vertices_by_dof_func(const cs_xdef_t   *def,
                                              const cs_lnum_t    n_v_selected,
                                              const cs_lnum_t   *selected_lst,
                                              cs_real_t          retval[])
{
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_DOF_FUNCTION);

  cs_xdef_dof_context_t *dcx
    = static_cast<cs_xdef_dof_context_t *>(def->context);

  /* Perform the evaluation */

  if (cs_cdo_quant->n_vertices == n_v_selected)
    dcx->func(n_v_selected,
              nullptr, /* elt_ids */
              false,   /* dense output ? */
              dcx->input,
              retval);
  else
    dcx->func(n_v_selected,
              selected_lst,   /* elt_ids */
              false,          /* dense output ? */
              dcx->input,
              retval);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a potential field at face centers from a definition by a
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  /* Perform the evaluation */

  const cs_lnum_t  n_faces = cs_cdo_quant->n_faces;
  const cs_real_t  *input = (cs_real_t *)def->context;
  const cs_lnum_t *elt_ids = (n_f_selected == n_faces) ? nullptr : selected_lst;

  _eval_by_value(n_f_selected, def->dim, elt_ids, input, retval);
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)def->context;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  /* Perform the evaluation */

  if (n_faces == n_f_selected) {

    /* All the support entities are selected:
       - First pass: interior faces
       - Second pass: border faces
    */

    ac->func(time_eval,
             quant->n_i_faces,
             nullptr,
             quant->i_face_center,
             true, /* Output is dense ? */
             ac->input,
             retval);
    ac->func(time_eval,
             quant->n_b_faces,
             nullptr,
             quant->b_face_center,
             true, /* Output is dense ? */
             ac->input,
             retval + def->dim * quant->n_i_faces);
  }
  else {

    assert(selected_lst != nullptr);

    /* Selected faces are stored by increasing number */

    cs_lnum_t  n_i_faces = 0;
    for (cs_lnum_t  i = 0; i < n_f_selected; i++) {
      if (selected_lst[i] < quant->n_i_faces)
        n_i_faces++;
      else
        break;
    }

    /* Interior faces */

    ac->func(time_eval,
             n_i_faces, selected_lst, quant->i_face_center,
             false, /* Output is dense ? */
             ac->input,
             retval);

    /* Border faces */

    cs_lnum_t n_b_faces = n_f_selected - n_i_faces;
    assert(n_b_faces > -1);
    ac->func(time_eval,
             n_b_faces, selected_lst + n_i_faces, quant->b_face_center,
             false, /* Output is dense ? */
             ac->input,
             retval);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate a potential field at cell centers from a definition by
 *         array
 *
 * \param[in]      def       pointer to a cs_xdef_t pointer
 * \param[in, out] retval    pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_cells_by_array(const cs_xdef_t   *def,
                                        cs_real_t          retval[])
{
  cs_evaluate_average_on_cells_by_array(def, retval);
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;
  const cs_real_t  *input = (cs_real_t *)def->context;
  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_lnum_t  *elt_ids = (z->n_elts == n_cells) ? nullptr : z->elt_ids;

  _eval_by_value(z->n_elts, def->dim, elt_ids, input, retval);
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  cs_xdef_analytic_context_t  *ac = (cs_xdef_analytic_context_t *)def->context;

  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  if (def->meta & CS_FLAG_FULL_LOC) /* All cells are selected */
    ac->func(time_eval,
             quant->n_cells,
             nullptr,
             quant->cell_centers,
             false, /* dense output */
             ac->input,
             retval);
  else
    ac->func(time_eval,
             z->n_elts, z->elt_ids, quant->cell_centers,
             false,  /* dense output */
             ac->input,
             retval);

  /* No sync since theses values are computed by only one rank */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the quantity attached to a potential field at cells
 *         when the definition relies on a DoF function (Degrees of freedom)
 *
 * \param[in]      def           pointer to a cs_xdef_t pointer
 * \param[in, out] retval        pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_potential_at_cells_by_dof_func(const cs_xdef_t   *def,
                                           cs_real_t          retval[])
{
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_DOF_FUNCTION);

  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

  cs_xdef_dof_context_t *dcx
    = static_cast<cs_xdef_dof_context_t *>(def->context);

  /* Perform the evaluation */

  if (cs_cdo_quant->n_cells == z->n_elts)
    dcx->func(z->n_elts,
              nullptr, /* elt_ids */
              false,   /* dense output ? */
              dcx->input,
              retval);
  else
    dcx->func(z->n_elts,
              z->elt_ids,   /* elt_ids */
              false,        /* dense output ? */
              dcx->input,
              retval);

  /* No sync since theses values are computed by only one rank */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a value to each DoF in the case of a potential field in order
 *         to put a given quantity inside the volume associated to the zone
 *         related to the given definition
 *         wvals may be nullptr.
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
  if (def == nullptr)
    return;
  if (vvals == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_QOV);

  const cs_real_t  *input = (cs_real_t *)def->context;
  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

  /* Perform the evaluation */

  bool check = false;
  if (dof_flag & CS_FLAG_SCALAR) { /* DoF is scalar-valued */

    const cs_real_t  const_val = input[0];

    if (cs_flag_test(dof_flag, cs_flag_primal_vtx | cs_flag_primal_cell)) {

      if (wvals == nullptr)
        bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

      _pvcsp_by_qov(const_val, z->n_elts, z->elt_ids, vvals, wvals);
      check = true;

    }
    else if (cs_flag_test(dof_flag, cs_flag_primal_vtx)) {

      _pvsp_by_qov(const_val, z->n_elts, z->elt_ids, vvals);
      check = true;

    }
    else if (cs_flag_test(dof_flag, cs_flag_primal_cell)) {

      _pcsp_by_qov(const_val, z->n_elts, z->elt_ids, vvals);
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  const cs_lnum_t  n_edges = cs_cdo_quant->n_edges;
  const cs_real_t  *edge_vector = cs_cdo_quant->edge_vector;
  const cs_real_t  *input = (cs_real_t *)def->context;

  /* DoF is scalar-valued since this is a circulation but the definition is
   * either scalar-valued meaning that one only gives the tangential part or
   * vector-valued (in this case, one needs to extract the tangential part) */

  switch (def->dim) {

  case 1: /* Scalar-valued integral */
    if (n_edges == n_e_selected)
      cs_array_real_set_scalar(n_edges, input[0], retval);
    else
      cs_array_real_set_scalar_on_subset(n_e_selected, selected_lst, input[0],
                                         retval);
    break;

  case 3:
    if (n_edges == n_e_selected) {

#     pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++)
        retval[e_id] = _dp3(input, edge_vector + 3*e_id);

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != nullptr);

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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ARRAY);

  const cs_lnum_t  n_edges = cs_cdo_quant->n_edges;
  const cs_real_t  *edge_vector = cs_cdo_quant->edge_vector;

  cs_xdef_array_context_t  *ac = (cs_xdef_array_context_t *)def->context;
  assert(cs_flag_test(ac->value_location, cs_flag_primal_edge));

  /* DoF is scalar-valued since this is a circulation but the definition is
   * either scalar-valued meaning that one only gives the tangential part or
   * vector-valued (in this case, one needs to extract the tangential part) */

  switch (def->dim) {

  case 1: /* Scalar-valued integral */
    assert(ac->stride == 1);
    if (n_edges == n_e_selected)
      cs_array_real_copy(n_edges, ac->values, retval);
    else {
      if (ac->full_length)
        cs_array_real_copy_subset(n_e_selected, 1, selected_lst,
                                  CS_ARRAY_SUBSET_INOUT,
                                  ac->values,
                                  retval);
      else
        cs_array_real_copy_subset(n_e_selected, 1, selected_lst,
                                  CS_ARRAY_SUBSET_OUT,
                                  ac->values,
                                  retval);
    }
    break;

  case 3:
    assert(ac->stride == 3);
    if (n_edges == n_e_selected) {

#     pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e_id = 0; e_id < n_edges; e_id++)
        retval[e_id] = _dp3(ac->values + 3*e_id, edge_vector + 3*e_id);

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != nullptr);

#     pragma omp parallel for if (n_e_selected > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_e_selected; e++) {
        const cs_lnum_t e_id = selected_lst[e];
        retval[e_id] = _dp3(ac->values + 3*e_id, edge_vector + 3*e_id);
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  const cs_lnum_t  n_edges = cs_cdo_quant->n_edges;
  const cs_real_t  *edge_vector = cs_cdo_quant->edge_vector;
  const cs_real_t  *xv = cs_cdo_quant->vtx_coord;
  const cs_adjacency_t  *e2v = cs_cdo_connect->e2v;

  cs_xdef_analytic_context_t *ac = (cs_xdef_analytic_context_t *)def->context;
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
              ac->func, ac->input, &integral);

        retval[e_id] = integral;

      }

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != nullptr);

#     pragma omp parallel for if (n_e_selected > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_e_selected; e++) {

        const cs_lnum_t e_id = selected_lst[e];
        const cs_lnum_t  *_v = e2v->ids + 2*e_id;

        cs_real_t  e_len = cs_math_3_norm(edge_vector + 3*e_id);
        cs_real_t  integral = 0.;
        qfunc(time_eval, xv + 3*_v[0], xv + 3*_v[1], e_len,
              ac->func, ac->input, &integral);

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
              ac->func, ac->input, integral);

        retval[e_id] = _dp3(integral, e_vec.unitv);

      }

    }
    else { /* A selection of edges is selected */

      assert(selected_lst != nullptr);

#     pragma omp parallel for if (n_e_selected > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_e_selected; e++) {

        const cs_lnum_t e_id = selected_lst[e];
        const cs_lnum_t  *_v = e2v->ids + 2*e_id;

        cs_nvec3_t  e_vec;
        cs_nvec3(edge_vector + 3*e_id, &e_vec);

        cs_real_3_t  integral = {0., 0., 0.};
        qfunc(time_eval, xv + 3*_v[0], xv + 3*_v[1], e_vec.meas,
              ac->func, ac->input, integral);

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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  /* Perform the evaluation */

  const cs_lnum_t  n_faces = cs_cdo_quant->n_faces;
  const cs_lnum_t *elt_ids = (n_f_selected == n_faces) ? nullptr : selected_lst;
  const cs_real_t  *values = (cs_real_t *)def->context;

  _eval_by_value(n_f_selected, def->dim, elt_ids, values, retval);
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  cs_quadrature_tria_integral_t
    *qfunc = cs_quadrature_get_tria_integral(def->dim, def->qtype);
  cs_xdef_analytic_context_t *ac = (cs_xdef_analytic_context_t *)def->context;

  switch (def->dim) {

  case 1: /* Scalar-valued */
    _pfsa_by_analytic(time_eval,
                      ac->func, ac->input,
                      n_f_selected, selected_lst,
                      qfunc,
                      retval);
    break;

  case 3: /* Vector-valued */
    _pfva_by_analytic(time_eval,
                      ac->func, ac->input,
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
 * \brief  Evaluate the average value on faces following the given definition
 *
 * \param[in]      def            pointer to a cs_xdef_t pointer
 * \param[in]      time_eval      physical time at which one evaluates the term
 * \param[in]      n_f_selected   number of selected faces
 * \param[in]      selected_lst   list of selected faces
 * \param[in, out] retval         pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_faces(const cs_xdef_t   *def,
                             cs_real_t          time_eval,
                             const cs_lnum_t    n_f_selected,
                             const cs_lnum_t   *selected_lst,
                             cs_real_t          retval[])
{
  if (def == nullptr)
    return;

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    cs_evaluate_average_on_faces_by_value(def,
                                          n_f_selected, selected_lst,
                                          retval);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_evaluate_average_on_faces_by_analytic(def,
                                             time_eval,
                                             n_f_selected, selected_lst,
                                             retval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Case not handled yet.", __func__);

  }
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_VALUE);

  const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;
  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_lnum_t  *elt_ids = (z->n_elts == n_cells) ? nullptr : z->elt_ids;
  const cs_real_t  *values = (cs_real_t *)def->context;

  _eval_by_value(z->n_elts, def->dim, elt_ids, values, retval);
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ARRAY);

  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_xdef_array_context_t  *ac = (cs_xdef_array_context_t *)def->context;
  const int  stride = ac->stride;
  const cs_real_t  *val = ac->values;

  if (cs_flag_test(ac->value_location, cs_flag_primal_cell) == false)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case. Not implemented yet.",
              __func__);
  assert(stride == def->dim);

  if (def->meta & CS_FLAG_FULL_LOC)
    cs_array_real_copy(stride*cs_cdo_quant->n_cells, val, retval);

  else {

    if (ac->full_length)
      cs_array_real_copy_subset(z->n_elts, stride, z->elt_ids,
                                CS_ARRAY_SUBSET_INOUT,
                                val,
                                retval);
    else
      cs_array_real_copy_subset(z->n_elts, stride, z->elt_ids,
                                CS_ARRAY_SUBSET_OUT,
                                val,
                                retval);

  } /* Deal with a selection of cells */
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
  if (def == nullptr)
    return;
  if (retval == nullptr)
    bft_error(__FILE__, __LINE__, 0, _err_empty_array, __func__);

  assert(def->support == CS_XDEF_SUPPORT_VOLUME);
  assert(def->type == CS_XDEF_BY_ANALYTIC_FUNCTION);

  const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;
  const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);
  const cs_lnum_t  *elt_ids = (n_cells == z->n_elts) ? nullptr : z->elt_ids;

  cs_quadrature_tetra_integral_t
    *qfunc = cs_quadrature_get_tetra_integral(def->dim, def->qtype);
  cs_xdef_analytic_context_t *ac = (cs_xdef_analytic_context_t *)def->context;

  switch (def->dim) {

  case 1: /* Scalar-valued */
    if (elt_ids == nullptr)
      cs_array_real_fill_zero(n_cells, retval);
    else
      cs_array_real_set_scalar_on_subset(z->n_elts, elt_ids, 0., retval);

    _pcsa_by_analytic(time_eval,
                      ac->func, ac->input, z->n_elts, elt_ids, qfunc,
                      retval);
    break;

  case 3: /* Vector-valued */
    if (elt_ids == nullptr)
      cs_array_real_fill_zero(3*n_cells, retval);
    else {
      cs_real_t  zero[3] = {0, 0, 0};
      cs_array_real_set_vector_on_subset(z->n_elts, elt_ids, zero, retval);
    }

    _pcva_by_analytic(time_eval,
                      ac->func, ac->input, z->n_elts, elt_ids, qfunc,
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
 * \brief  Evaluate the average value on cells following the given definition
 *         The cells associated to this definition (through the related zone)
 *         are all considered.
 *
 * \param[in]      def        pointer to a cs_xdef_t pointer
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] retval     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_evaluate_average_on_cells(const cs_xdef_t   *def,
                             cs_real_t          time_eval,
                             cs_real_t          retval[])
{
  if (def == nullptr)
    return;

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    cs_evaluate_average_on_cells_by_value(def, retval);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_evaluate_average_on_cells_by_analytic(def, time_eval, retval);
    break;

  case CS_XDEF_BY_ARRAY:
    cs_evaluate_average_on_cells_by_array(def, retval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Case not handled yet.", __func__);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Evaluate the integral over the full computational domain of a
 *         quantity defined by an array. The parallel sum reduction is
 *         performed inside this function.
 *
 * \param[in]  array_loc   flag indicating where are located values
 * \param[in]  array_val   array of values
 *
 * \return the value of the integration (parallel sum reduction done)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_evaluate_scal_domain_integral_by_array(cs_flag_t         array_loc,
                                          const cs_real_t  *array_val)
{
  cs_real_t  result = 0.;

  if (array_val == nullptr)
    return result;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  if (cs_flag_test(array_loc, cs_flag_primal_cell)) {

#   pragma omp parallel for reduction(+:result)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      result += array_val[c_id] * quant->cell_vol[c_id];

  }
  else if (cs_flag_test(array_loc, cs_flag_primal_vtx)) {

    const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;
    const cs_real_t  *wvc = quant->pvol_vc;
    assert(wvc != nullptr);

#   pragma omp parallel for reduction(+:result)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
        result += wvc[j] * array_val[c2v->ids[j]];

  }
  else if (cs_flag_test(array_loc, cs_flag_primal_edge)) {

    const cs_adjacency_t  *c2e = cs_cdo_connect->c2e;
    const cs_real_t  *wec = quant->pvol_ec;
    assert(wec != nullptr);

#   pragma omp parallel for reduction(+:result)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++)
        result += wec[j] * array_val[c2e->ids[j]];

  }
  else if (cs_flag_test(array_loc, cs_flag_primal_face)) {

    const cs_adjacency_t  *c2f = cs_cdo_connect->c2f;
    const cs_real_t  *wfc = quant->pvol_fc;
    assert(wfc != nullptr);

#   pragma omp parallel for reduction(+:result)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
      for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
        result += wfc[j] * array_val[c2f->ids[j]];

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid array location. Stop evaluation.", __func__);

  cs_parall_sum(1, CS_REAL_TYPE, &result);

  return result;
}

/*----------------------------------------------------------------------------*/

#undef _dp3
END_C_DECLS
