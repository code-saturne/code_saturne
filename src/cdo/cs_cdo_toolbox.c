/*============================================================================
 * Set of toolbox functions: shared buffer, balance, synchronization
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_cdo_local.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_xdef_eval.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_toolbox.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_CDO_TOOLBOX_DBG        0 /* Debug level */

/*============================================================================
 * Local private variables
 *============================================================================*/

/* Temporary buffers useful during the building of all algebraic systems */

static size_t  cs_cdo_toolbox_work_buffer_size = 0;
static cs_real_t  *cs_cdo_toolbox_work_buffer = NULL;

/* Shared structure */

static const cs_cdo_connect_t  *cs_shared_connect;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *         Set also shared pointers from the main domain members
 *
 * \param[in]  connect      pointer to a cs_cdo_connect_t structure
 * \param[in]  eb_flag      metadata for Edge-based schemes
 * \param[in]  fb_flag      metadata for Face-based schemes
 * \param[in]  vb_flag      metadata for Vertex-based schemes
 * \param[in]  vcb_flag     metadata for Vertex+Cell-basde schemes
 * \param[in]  hho_flag     metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_toolbox_init(const cs_cdo_connect_t       *connect,
                    cs_flag_t                     eb_flag,
                    cs_flag_t                     fb_flag,
                    cs_flag_t                     vb_flag,
                    cs_flag_t                     vcb_flag,
                    cs_flag_t                     hho_flag)
{
  assert(connect != NULL); /* Sanity check */

  /* Allocate cell-wise and face-wise view of a mesh */

  cs_cdo_local_initialize(connect);

  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];
  const cs_lnum_t  n_vertices = connect->n_vertices;
  const cs_lnum_t  n_edges = connect->n_edges;

  /* Allocate shared buffer and initialize shared structures */

  size_t  cwb_size = n_cells; /* initial cell-wise buffer size */

  /* Allocate and initialize matrix assembler and matrix structures */

  if (vb_flag > 0 || vcb_flag > 0) {

    if (vb_flag & CS_FLAG_SCHEME_SCALAR || vcb_flag & CS_FLAG_SCHEME_SCALAR) {

      if (vb_flag & CS_FLAG_SCHEME_SCALAR)
        cwb_size = CS_MAX(cwb_size, (size_t)n_vertices);

      if (vcb_flag & CS_FLAG_SCHEME_SCALAR)
        cwb_size = CS_MAX(cwb_size, (size_t)(n_vertices + n_cells));

    } /* scalar-valued equations */

    if (vb_flag & CS_FLAG_SCHEME_VECTOR || vcb_flag & CS_FLAG_SCHEME_VECTOR) {

      cwb_size = CS_MAX(cwb_size, (size_t)3*n_cells);
      if (vb_flag & CS_FLAG_SCHEME_VECTOR)
        cwb_size = CS_MAX(cwb_size, (size_t)3*n_vertices);

      if (vcb_flag & CS_FLAG_SCHEME_VECTOR)
        cwb_size = CS_MAX(cwb_size, (size_t)3*(n_vertices + n_cells));

    } /* vector-valued equations */

  } /* Vertex-based schemes and related ones */

  if (eb_flag > 0) {

    if (eb_flag & CS_FLAG_SCHEME_SCALAR) {

      /* This is a vector-valued equation but the DoF is scalar-valued since
       * it is a circulation associated to each edge */

      cwb_size = CS_MAX(cwb_size, (size_t)3*n_cells);
      cwb_size = CS_MAX(cwb_size, (size_t)n_edges);

    } /* vector-valued equations with scalar-valued DoFs */

  } /* Edge-based schemes */

  if (fb_flag > 0 || hho_flag > 0) {

    if (cs_flag_test(fb_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_SCALAR) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_SCALAR)) {

      assert(n_faces > n_cells);
      if (fb_flag & CS_FLAG_SCHEME_SCALAR)
        cwb_size = CS_MAX(cwb_size, (size_t)n_faces);

      if (hho_flag & CS_FLAG_SCHEME_SCALAR)
        cwb_size = CS_MAX(cwb_size, (size_t)n_faces);

    } /* Scalar-valued CDO-Fb or HHO-P0 */

    if (cs_flag_test(fb_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_VECTOR) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_POLY1 | CS_FLAG_SCHEME_SCALAR) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_POLY0 | CS_FLAG_SCHEME_VECTOR)) {

      assert((CS_DOF_FACE_SCAP1 == CS_DOF_FACE_VECT) &&
             (CS_DOF_FACE_SCAP1 == CS_DOF_FACE_VECP0));

      cwb_size = CS_MAX(cwb_size, (size_t)CS_N_DOFS_FACE_1ST * n_faces);

    } /* Vector CDO-Fb or HHO-P1 or vector HHO-P0 */

    if (cs_flag_test(hho_flag,
                     CS_FLAG_SCHEME_POLY2 | CS_FLAG_SCHEME_SCALAR))
      cwb_size = CS_MAX(cwb_size, (size_t)CS_N_DOFS_FACE_2ND * n_faces);

    /* For vector equations and HHO */

    if (cs_flag_test(hho_flag, CS_FLAG_SCHEME_VECTOR | CS_FLAG_SCHEME_POLY1) ||
        cs_flag_test(hho_flag, CS_FLAG_SCHEME_VECTOR | CS_FLAG_SCHEME_POLY2)) {

      if  (hho_flag & CS_FLAG_SCHEME_POLY1)
        cwb_size = CS_MAX(cwb_size, (size_t)3*CS_N_DOFS_FACE_1ST*n_faces);

      else if  (hho_flag & CS_FLAG_SCHEME_POLY2)
        cwb_size = CS_MAX(cwb_size, (size_t)3*CS_N_DOFS_FACE_2ND*n_faces);

    }

  } /* Face-based schemes (CDO or HHO) */

  /* Assign static const pointers: shared pointers with a cs_domain_t */

  cs_shared_connect = connect;

  /* Common buffer for temporary usage */

  cs_cdo_toolbox_work_buffer_size = cwb_size;
  BFT_MALLOC(cs_cdo_toolbox_work_buffer, cwb_size, double);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free buffers shared among the equations solved with CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_toolbox_finalize(void)
{
  /* Free cell-wise and face-wise view of a mesh */

  cs_cdo_local_finalize();

  /* Free common buffer */

  BFT_FREE(cs_cdo_toolbox_work_buffer);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to a buffer of size at least the 2*n_cells
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdo_toolbox_get_tmpbuf(void)
{
  return cs_cdo_toolbox_work_buffer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the allocation size of the temporary buffer
 *
 * \return  the size of the temporary buffer
 */
/*----------------------------------------------------------------------------*/

size_t
cs_cdo_toolbox_get_tmpbuf_size(void)
{
  return cs_cdo_toolbox_work_buffer_size;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cdo_balance_t structure
 *
 * \param[in]  location   where the balance is performed
 * \param[in]  size       size of arrays in the structure
 *
 * \return  a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_balance_t *
cs_cdo_balance_create(cs_flag_t    location,
                      cs_lnum_t    size)
{
  cs_cdo_balance_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_cdo_balance_t);

  b->size = size;
  b->location = location;
  if (cs_flag_test(location, cs_flag_primal_cell) == false &&
      cs_flag_test(location, cs_flag_primal_vtx) == false)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid location", __func__);

  BFT_MALLOC(b->balance, 7*size, cs_real_t);
  b->unsteady_term  = b->balance +   size;
  b->reaction_term  = b->balance + 2*size;
  b->diffusion_term = b->balance + 3*size;
  b->advection_term = b->balance + 4*size;
  b->source_term    = b->balance + 5*size;
  b->boundary_term  = b->balance + 6*size;

  /* Set to zero all members */

  cs_cdo_balance_reset(b);

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset a cs_cdo_balance_t structure
 *
 * \param[in, out] b     pointer to a cs_cdo_balance_t to reset
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_balance_reset(cs_cdo_balance_t   *b)
{
  if (b == NULL)
    return;
  if (b->size < 1)
    return;

  if (b->balance == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: array is not allocated.", __func__);

  size_t  bufsize = b->size *sizeof(cs_real_t);

  memset(b->balance, 0, 7*bufsize);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize balance terms if this is a parallel computation
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in, out] b         pointer to a cs_cdo_balance_t to rsync
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_balance_sync(const cs_cdo_connect_t    *connect,
                    cs_cdo_balance_t          *b)
{
  if (b == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: structure not allocated", __func__);

  if (cs_flag_test(b->location, cs_flag_primal_vtx)) {

    assert(b->size == connect->n_vertices);

    if (connect->vtx_ifs != NULL)
      cs_interface_set_sum(connect->vtx_ifs,
                           b->size,
                           7,   /* stride: 1 for each kind of balance */
                           false,
                           CS_REAL_TYPE,
                           b->balance);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cdo_balance_t structure
 *
 * \param[in, out]  p_balance  pointer to the pointer to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_balance_destroy(cs_cdo_balance_t   **p_balance)
{
  cs_cdo_balance_t *b = *p_balance;

  if (b == NULL)
    return;

  BFT_FREE(b->balance);

  BFT_FREE(b);
  *p_balance = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the volumetric definitions to consider at each vertex
 *
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2v_idx   index array  to define
 * \param[in, out]  def2v_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_sync_vol_def_at_vertices(int                      n_defs,
                                cs_xdef_t              **defs,
                                cs_lnum_t                def2v_idx[],
                                cs_lnum_t                def2v_ids[])
{
  if (n_defs == 0)
    return;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->n_vertices;
  const cs_adjacency_t  *c2v = connect->c2v;

  int  *v2def_ids = NULL;
  BFT_MALLOC(v2def_ids, n_vertices, int);
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t v = 0; v < n_vertices; v++)
    v2def_ids[v] = -1;          /* default */

  for (int def_id = 0; def_id < n_defs; def_id++) {

    /* Get and then set the definition of the initial condition */

    const cs_xdef_t  *def = defs[def_id];
    assert(def->support == CS_XDEF_SUPPORT_VOLUME);

    if (def->meta & CS_FLAG_FULL_LOC) {

#     pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t v = 0; v < n_vertices; v++)
        v2def_ids[v] = def_id;

    }
    else {

      const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

      for (cs_lnum_t i = 0; i < z->n_elts; i++) { /* Loop on selected cells */
        const cs_lnum_t  c_id = z->elt_ids[i];
        for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
          v2def_ids[c2v->ids[j]] = def_id;
      }

    }

  } /* Loop on definitions */

  if (connect->vtx_ifs != NULL) {

    /* Last definition is used in case of conflict */

    cs_interface_set_max(connect->vtx_ifs,
                         n_vertices,
                         1,             /* stride */
                         false,         /* interlace (not useful here) */
                         CS_INT_TYPE,   /* int */
                         v2def_ids);

  }

  /* 0. Initialization */

  cs_lnum_t  *count = NULL;
  BFT_MALLOC(count, n_defs, cs_lnum_t);
  memset(count, 0, n_defs*sizeof(cs_lnum_t));
  memset(def2v_idx, 0, (n_defs+1)*sizeof(cs_lnum_t));

  /* 1. Count number of vertices related to each definition */

  for (cs_lnum_t v = 0; v < n_vertices; v++)
    if (v2def_ids[v] > -1)
      def2v_idx[v2def_ids[v]+1] += 1;

  /* 2. Build index */

  for (int def_id = 0; def_id < n_defs; def_id++)
    def2v_idx[def_id+1] += def2v_idx[def_id];

  /* 3. Build list */

  for (cs_lnum_t v = 0; v < n_vertices; v++) {
    const int def_id = v2def_ids[v];
    if (def_id > -1) {
      def2v_ids[def2v_idx[def_id] + count[def_id]] = v;
      count[def_id] += 1;
    }
  }

  BFT_FREE(v2def_ids);
  BFT_FREE(count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the volumetric definitions to consider at each edge
 *
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2e_idx   index array  to define
 * \param[in, out]  def2e_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_sync_vol_def_at_edges(int                      n_defs,
                             cs_xdef_t              **defs,
                             cs_lnum_t                def2e_idx[],
                             cs_lnum_t                def2e_ids[])
{
  if (n_defs == 0)
    return;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_edges = connect->n_edges;
  const cs_adjacency_t  *c2e = connect->c2e;

  int  *e2def_ids = NULL;
  BFT_MALLOC(e2def_ids, n_edges, int);
# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t e = 0; e < n_edges; e++)
    e2def_ids[e] = -1; /* default: not associated to a definition */

  for (int def_id = 0; def_id < n_defs; def_id++) {

    /* Get and then set the definition of the initial condition */

    const cs_xdef_t  *def = defs[def_id];
    assert(def->support == CS_XDEF_SUPPORT_VOLUME);

    if (def->meta & CS_FLAG_FULL_LOC) {

#     pragma omp parallel for if (n_edges > CS_THR_MIN)
      for (cs_lnum_t e = 0; e < n_edges; e++)
        e2def_ids[e] = def_id;

    }
    else {

      const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

      for (cs_lnum_t i = 0; i < z->n_elts; i++) { /* Loop on selected cells */
        const cs_lnum_t  c_id = z->elt_ids[i];
        for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++)
          e2def_ids[c2e->ids[j]] = def_id;
      }

    }

  } /* Loop on definitions */

  if (connect->edge_ifs != NULL) {

    /* Last definition is used in case of conflict */

    cs_interface_set_max(connect->edge_ifs,
                         n_edges,
                         1,             /* stride */
                         false,         /* interlace (not useful here) */
                         CS_INT_TYPE,   /* int */
                         e2def_ids);

  }

  /* 0. Initialization */

  cs_lnum_t  *count = NULL;
  BFT_MALLOC(count, n_defs, cs_lnum_t);
  memset(count, 0, n_defs*sizeof(cs_lnum_t));
  memset(def2e_idx, 0, (n_defs+1)*sizeof(cs_lnum_t));

  /* 1. Count the number of edges related to each definition */

  for (cs_lnum_t e = 0; e < n_edges; e++)
    if (e2def_ids[e] > -1)
      def2e_idx[e2def_ids[e]+1] += 1;

  /* 2. Build the index */

  for (int def_id = 0; def_id < n_defs; def_id++)
    def2e_idx[def_id+1] += def2e_idx[def_id];

  /* 3. Build the list */

  for (cs_lnum_t e = 0; e < n_edges; e++) {
    const int def_id = e2def_ids[e];
    if (def_id > -1) {
      def2e_ids[def2e_idx[def_id] + count[def_id]] = e;
      count[def_id] += 1;
    }
  }

  BFT_FREE(e2def_ids);
  BFT_FREE(count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the volumetric definitions to consider at each face
 *
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2f_idx   index array  to define
 * \param[in, out]  def2f_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_sync_vol_def_at_faces(int                        n_defs,
                             cs_xdef_t                **defs,
                             cs_lnum_t                  def2f_idx[],
                             cs_lnum_t                  def2f_ids[])
{
  if (n_defs == 0)
    return;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];
  const cs_adjacency_t  *c2f = connect->c2f;

  int  *f2def_ids = NULL;
  BFT_MALLOC(f2def_ids, n_faces, int);
# pragma omp parallel for if (n_faces > CS_THR_MIN)
  for (cs_lnum_t f = 0; f < n_faces; f++)
    f2def_ids[f] = -1;          /* default */

  for (int def_id = 0; def_id < n_defs; def_id++) {

    /* Get and then set the definition of the initial condition */

    const cs_xdef_t  *def = defs[def_id];
    assert(def->support == CS_XDEF_SUPPORT_VOLUME);

    if (def->meta & CS_FLAG_FULL_LOC) {

#     pragma omp parallel for if (n_faces > CS_THR_MIN)
      for (cs_lnum_t f = 0; f < n_faces; f++)
        f2def_ids[f] = def_id;

    }
    else {

      const cs_zone_t  *z = cs_volume_zone_by_id(def->z_id);

      for (cs_lnum_t i = 0; i < z->n_elts; i++) { /* Loop on selected cells */
        const cs_lnum_t  c_id = z->elt_ids[i];
        for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
          f2def_ids[c2f->ids[j]] = def_id;
      }

    }

  } /* Loop on definitions */

  if (connect->face_ifs != NULL) {

    /* Last definition is used in case of conflict */

    cs_interface_set_max(connect->face_ifs,
                         n_faces,
                         1,             /* stride */
                         false,         /* interlace (not useful here) */
                         CS_INT_TYPE,   /* int */
                         f2def_ids);

  }

  /* 0. Initialization */

  cs_lnum_t  *count = NULL;
  BFT_MALLOC(count, n_defs, cs_lnum_t);
  memset(count, 0, n_defs*sizeof(cs_lnum_t));
  memset(def2f_idx, 0, (n_defs+1)*sizeof(cs_lnum_t));

  /* 1. Count number of faces related to each definition */

  for (cs_lnum_t f = 0; f < n_faces; f++)
    if (f2def_ids[f] > -1)
      def2f_idx[f2def_ids[f]+1] += 1;

  /* 2. Build index */

  for (int def_id = 0; def_id < n_defs; def_id++)
    def2f_idx[def_id+1] += def2f_idx[def_id];

  /* 3. Build list */

  for (cs_lnum_t f = 0; f < n_faces; f++) {
    const int def_id = f2def_ids[f];
    if (def_id > -1) {
      def2f_ids[def2f_idx[def_id] + count[def_id]] = f;
      count[def_id] += 1;
    }
  }

  BFT_FREE(f2def_ids);
  BFT_FREE(count);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the mean-value across ranks at each vertex
 *
 * \param[in]       dim         number of entries for each vertex
 * \param[in]       counter     number of occurences on this rank
 * \param[in, out]  values      array to update
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_sync_vertex_mean_values(int                         dim,
                               int                        *counter,
                               cs_real_t                  *values)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->n_vertices;

  if (connect->vtx_ifs != NULL) {

    cs_interface_set_sum(connect->vtx_ifs,
                         n_vertices,
                         1,           /* stride */
                         false,       /* interlace (not useful here) */
                         CS_INT_TYPE, /* int */
                         counter);

    cs_interface_set_sum(connect->vtx_ifs,
                         n_vertices,
                         dim,         /* stride */
                         true,        /* interlace */
                         CS_REAL_TYPE,
                         values);

  }

  if (dim == 1) {

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++)
      if (counter[v_id] > 1)
        values[v_id] /= counter[v_id];

  }
  else { /* dim > 1 */

#   pragma omp parallel for if (n_vertices > CS_THR_MIN)
    for (cs_lnum_t v_id = 0; v_id < n_vertices; v_id++) {
      if (counter[v_id] > 1) {
        const cs_real_t  inv_count = 1./counter[v_id];
        for (int k = 0; k < dim; k++)
          values[dim*v_id + k] *= inv_count;
      }
    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
