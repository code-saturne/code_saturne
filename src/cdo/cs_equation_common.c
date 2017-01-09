/*============================================================================
 * Routines to handle common features for building algebraic system in CDO
 * schemes
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_local.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_cdofb_scaleq.h"
#include "cs_hho_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_common.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

/* Main categories to consider for high-level matrix structures */
#define CS_EQ_COMMON_VERTEX   0
#define CS_EQ_COMMON_FACE     1
#define CS_EQ_N_COMMONS       2

/*============================================================================
 * Local private variables
 *============================================================================*/

/* Temporary buffers useful during the building of all algebraic systems */
static size_t  cs_equation_common_work_buffer_size = 0;
static cs_real_t  *cs_equation_common_work_buffer = NULL;

/* Store the matrix structure and its assembler structures for each family
   of space discretizations */
static cs_matrix_assembler_t  **cs_equation_common_ma = NULL;
static cs_matrix_structure_t  **cs_equation_common_ms = NULL;

/* Structure related to the index of a matrix for vertex-based schemes
   vertex --> vertices through cell connectivity */
cs_connect_index_t  *cs_connect_v2v = NULL;

/* Structure related to the index of a matrix for face-based schemes
   face --> faces through cell connectivity */
cs_connect_index_t  *cs_connect_f2f = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a vertex -> vertices connectivity index
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_connect_index_t *
_get_v2v(const cs_cdo_connect_t     *connect)
{
  /* Build a (sorted) v2v connectivity index */
  const cs_lnum_t  n_vertices = connect->v_info->n_elts;
  const cs_connect_index_t  *c2v = connect->c2v;

  cs_connect_index_t  *v2c = cs_index_transpose(n_vertices, c2v);
  cs_connect_index_t  *v2v = cs_index_compose(n_vertices, v2c, c2v);

  cs_index_sort(v2v);

  /* Update index (v2v has a diagonal entry. We remove it since we have in
     mind an index structure for a  matrix stored using the MSR format */
  cs_lnum_t  shift = 0;
  cs_lnum_t  prev_start = v2v->idx[0];
  cs_lnum_t  prev_end = v2v->idx[1];

  for (cs_lnum_t i = 0; i < n_vertices; i++) {

    for (cs_lnum_t j = prev_start; j < prev_end; j++)
      if (v2v->ids[j] != i)
        v2v->ids[shift++] = v2v->ids[j];

    if (i != n_vertices - 1) { // Update prev_start and prev_end
      prev_start = v2v->idx[i+1];
      prev_end = v2v->idx[i+2];
    }
    v2v->idx[i+1] = shift;

  } // Loop on vertices

  /* Free temporary buffers */
  cs_index_free(&v2c);

  return v2v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a face -> faces connectivity index
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_connect_index_t *
_get_f2f(const cs_cdo_connect_t     *connect)
{
  cs_connect_index_t  *c2f = NULL, *f2c = NULL, *f2f = NULL;

  const cs_lnum_t  n_faces = connect->f_info->n_elts;
  const cs_sla_matrix_t *mc2f = connect->c2f;
  const cs_sla_matrix_t *mf2c = connect->f2c;

  /* Build a face -> face connectivity */
  f2c = cs_index_map(mf2c->n_rows, mf2c->idx, mf2c->col_id);
  c2f = cs_index_map(mc2f->n_rows, mc2f->idx, mc2f->col_id);
  f2f = cs_index_compose(n_faces, f2c, c2f);
  cs_index_sort(f2f);

  /* Free temporary memory */
  cs_index_free(&f2c);
  cs_index_free(&c2f);

  /* Update index (v2v has a diagonal entry. We remove it since we have in
     mind an index structure for a  matrix stored using the MSR format */
  cs_lnum_t  shift = 0;
  cs_lnum_t  prev_start = f2f->idx[0];
  cs_lnum_t  prev_end = f2f->idx[1];

  for (cs_lnum_t i = 0; i < n_faces; i++) {

    for (cs_lnum_t j = prev_start; j < prev_end; j++)
      if (f2f->ids[j] != i)
        f2f->ids[shift++] = f2f->ids[j];

    if (i != n_faces - 1) { // Update prev_start and prev_end
      prev_start = f2f->idx[i+1];
      prev_end = f2f->idx[i+2];
    }
    f2f->idx[i+1] = shift;

  } // Loop on faces

  return f2f;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_matrix_assembler_t structure
 *
 * \param[in]      n_elts   number of elements
 * \param[in]      x2x      pointer to a cs_connect_index_t struct.
 * \param[in]      rs       pointer to a range set or NULL if sequential
 * \param[in, out] ma       pointer to the cs_matrix_assembler_t to update
 */
/*----------------------------------------------------------------------------*/

static void
_build_matrix_assembler(cs_lnum_t                   n_elts,
                        const cs_connect_index_t   *x2x,
                        const cs_range_set_t       *rs,
                        cs_matrix_assembler_t      *ma)
{
  cs_gnum_t  *grows = NULL, *gcols = NULL;

  /* First loop to count max size of the buffer */
  cs_lnum_t  max_size = 0;
  for (cs_lnum_t id = 0; id < n_elts; id++)
    max_size = CS_MAX(max_size, x2x->idx[id+1] - x2x->idx[id]);
  BFT_MALLOC(grows, max_size + 1, cs_gnum_t); // +1 for the diagonal entry
  BFT_MALLOC(gcols, max_size + 1, cs_gnum_t);

  for (cs_lnum_t id = 0; id < n_elts; id++) {

    cs_lnum_t  start = x2x->idx[id];
    cs_lnum_t  end = x2x->idx[id+1];
    cs_gnum_t  grow_id = rs->g_id[id];

    /* Diagonal term is excluded in this connectivity. Add it "manually" */
    grows[0] = grow_id, gcols[0] = grow_id;
    for (cs_lnum_t j = start, i = 1; j < end; j++, i++) {
      grows[i] = grow_id;
      gcols[i] = rs->g_id[x2x->ids[j]];
    }

    cs_matrix_assembler_add_g_ids(ma, end-start+1, grows, gcols);

  } // Loop on entities

  /* Now compute structure */
  cs_matrix_assembler_compute(ma);

  /* Free temporary buffers */
  BFT_FREE(grows);
  BFT_FREE(gcols);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *         Set also shared pointers from the main domain members
 *
 * \param[in]  connect       pointer to a cs_cdo_connect_t structure
 * \param[in]  quant         pointer to additional mesh quantities struct.
 * \param[in]  time_step     pointer to a time step structure
 * \param[in]  scheme_flag   flag to identify which kind of numerical scheme is
 *                           requested to solve the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_allocate_common_structures(const cs_cdo_connect_t     *connect,
                                       const cs_cdo_quantities_t  *quant,
                                       const cs_time_step_t       *time_step,
                                       cs_flag_t                   scheme_flag)
{
  assert(connect != NULL); // Sanity check

  /* Two types of mat. ass. are considered:
     - The one related to matrix based on vertices
     - The one related to matrix based on faces
  */
  BFT_MALLOC(cs_equation_common_ma, CS_EQ_N_COMMONS, cs_matrix_assembler_t *);
  for (int i = 0; i < CS_EQ_N_COMMONS; i++)
    cs_equation_common_ma[i] = NULL;

  BFT_MALLOC(cs_equation_common_ms, CS_EQ_N_COMMONS, cs_matrix_structure_t *);
  for (int i = 0; i < CS_EQ_N_COMMONS; i++)
    cs_equation_common_ms[i] = NULL;

  /* Allocate cell-wise and face-wise view of a mesh */
  cs_cdo_local_initialize(connect);

  const cs_lnum_t  n_cells = connect->c_info->n_elts;
  const cs_lnum_t  n_faces = connect->f_info->n_elts;
  const cs_lnum_t  n_vertices = connect->v_info->n_elts;

  /* Allocate and initialize matrix assembler and matrix structures */
  if (scheme_flag & CS_SCHEME_FLAG_CDOVB ||
      scheme_flag & CS_SCHEME_FLAG_CDOVCB) {

    /* Build the "v2v" connectivity index */
    cs_connect_v2v = _get_v2v(connect);

    cs_matrix_assembler_t  *ma =
      cs_matrix_assembler_create(connect->v_rs->l_range, true); // sep_diag

    _build_matrix_assembler(n_vertices, cs_connect_v2v, connect->v_rs, ma);

    cs_matrix_structure_t  *ms =
      cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

    cs_equation_common_ma[CS_EQ_COMMON_VERTEX] = ma;
    cs_equation_common_ms[CS_EQ_COMMON_VERTEX] = ms;

  } // Vertex-based schemes and related ones

  if (scheme_flag & CS_SCHEME_FLAG_CDOFB || scheme_flag & CS_SCHEME_FLAG_HHO) {

    /* Build the "f2f" connectivity index */
    cs_connect_f2f = _get_f2f(connect);

    cs_matrix_assembler_t  *ma =
      cs_matrix_assembler_create(connect->f_rs->l_range, true); // sep_diag

    _build_matrix_assembler(n_faces, cs_connect_f2f, connect->f_rs, ma);

    cs_matrix_structure_t  *ms =
      cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR, ma);

    cs_equation_common_ma[CS_EQ_COMMON_FACE] = ma;
    cs_equation_common_ms[CS_EQ_COMMON_FACE] = ms;

  } // Face-based schemes and related ones

  /* Allocate shared buffer and initialize shared structures */
  size_t  cwb_size = 2*n_cells; // initial cell-wise buffer size

  if (scheme_flag & CS_SCHEME_FLAG_CDOVB &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    cwb_size = CS_MAX(cwb_size, (size_t)3*n_vertices);

    /* Initialize additional structures */
    cs_cdovb_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_cdovb_scaleq_initialize();

  }

  if (scheme_flag & CS_SCHEME_FLAG_CDOVCB &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    cwb_size = CS_MAX(cwb_size, (size_t)2*(n_vertices + n_cells));

    /* Initialize additional structures */
    cs_cdovcb_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_cdovcb_scaleq_initialize();
  }

  if (scheme_flag & CS_SCHEME_FLAG_CDOFB &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    /* Initialize additional structures */
    cs_cdofb_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_cdofb_scaleq_initialize();

    cwb_size = CS_MAX(cwb_size, (size_t)3*n_faces);

  }

  if (scheme_flag & CS_SCHEME_FLAG_HHO &&
      scheme_flag & CS_SCHEME_FLAG_SCALAR) {

    /* Initialize additional structures */
    cs_hho_scaleq_set_shared_pointers(quant, connect, time_step);
    cs_hho_scaleq_initialize();

    // TODO: Update this value accordingly to HHO needs
    cwb_size = CS_MAX(cwb_size, (size_t)3*n_faces);

  }

  cs_equation_common_work_buffer_size = cwb_size;
  BFT_MALLOC(cs_equation_common_work_buffer, cwb_size, double);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the 2*n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         Call specific structure allocation related to a numerical scheme
 *         according the scheme flag
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \param[in]  scheme_flag   flag to identify which kind of numerical scheme is
 *                           requested to solve the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_free_common_structures(cs_flag_t   scheme_flag)
{
  /* Free cell-wise and face-wise view of a mesh */
  cs_cdo_local_finalize();

  if (scheme_flag & CS_SCHEME_FLAG_CDOVB || scheme_flag & CS_SCHEME_FLAG_CDOVCB)
    cs_index_free(&(cs_connect_v2v));

  if (scheme_flag & CS_SCHEME_FLAG_CDOFB || scheme_flag & CS_SCHEME_FLAG_HHO)
    cs_index_free(&(cs_connect_f2f));

  /* Free common structures specific to a numerical scheme */
  if ((scheme_flag & CS_SCHEME_FLAG_CDOVB) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_cdovb_scaleq_finalize();

  if ((scheme_flag & CS_SCHEME_FLAG_CDOVCB) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_cdovcb_scaleq_finalize();

  if ((scheme_flag & CS_SCHEME_FLAG_CDOFB) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_cdofb_scaleq_finalize();

  if ((scheme_flag & CS_SCHEME_FLAG_HHO) &&
      (scheme_flag & CS_SCHEME_FLAG_SCALAR))
    cs_hho_scaleq_finalize();

  BFT_FREE(cs_equation_common_work_buffer);

  /* matrix assemblers and structures */
  for (int i = 0; i < CS_EQ_N_COMMONS; i++) {
    cs_matrix_structure_destroy(&(cs_equation_common_ms[i]));
    cs_matrix_assembler_destroy(&(cs_equation_common_ma[i]));
  }
  BFT_FREE(cs_equation_common_ms);
  BFT_FREE(cs_equation_common_ma);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assemble a cellwise system related to cell vertices into the global
 *         algebraic system
 *
 * \param[in]       csys      cellwise view of the algebraic system
 * \param[in]       rset      pointer to a cs_range_set_t structure on vertices
 * \param[in]       sys_flag  flag associated to the current system builder
 * \param[in, out]  rhs       array storing the right-hand side
 * \param[in, out]  sources   array storing the contribution of source terms
 * \param[in, out]  mav       pointer to a matrix assembler structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_assemble_v(const cs_cell_sys_t            *csys,
                       const cs_range_set_t           *rset,
                       cs_flag_t                       sys_flag,
                       cs_real_t                      *rhs,
                       cs_real_t                      *sources,
                       cs_matrix_assembler_values_t   *mav)
{
  const short int  n_vc = csys->mat->n_ent;
  const cs_lnum_t  *v_ids = csys->mat->ids;

  cs_gnum_t  grows[CS_CDO_ASSEMBLE_BUF_SIZE], gcols[CS_CDO_ASSEMBLE_BUF_SIZE];
  cs_real_t  vals[CS_CDO_ASSEMBLE_BUF_SIZE];

  /* Assemble the matrix related to the advection/diffusion/reaction terms
     If advection is activated, the resulting system is not symmetric
     Otherwise, the system is symmetric with extra-diagonal terms. */
  /* TODO: Add a symmetric version for optimization */
  int  block_size = 0;
  for (short int i = 0; i < n_vc; i++) {

    const double  *mval_i = csys->mat->val + i*n_vc;
    const cs_lnum_t  vi_id = v_ids[i];

    rhs[vi_id] += csys->rhs[i];
    if (sys_flag & CS_FLAG_SYS_SOURCETERM)
      sources[vi_id] += csys->source[i];

    cs_gnum_t  grow_id = rset->g_id[vi_id];

    /* Diagonal term is excluded in this connectivity. Add it "manually" */
    for (short int j = 0; j < n_vc; j++) {

      grows[block_size] = grow_id;
      gcols[block_size] = rset->g_id[v_ids[j]];
      vals[block_size] = mval_i[j];
      block_size += 1;

      if (block_size == CS_CDO_ASSEMBLE_BUF_SIZE) {
        cs_matrix_assembler_values_add_g(mav,
                                         CS_CDO_ASSEMBLE_BUF_SIZE,
                                         grows, gcols, vals);
        block_size = 0;
      }

    } /* Loop on cell vertices (local cols) */

  } /* Loop on cell vertices (local rows) */

  if (block_size > 0)
    cs_matrix_assembler_values_add_g(mav, block_size, grows, gcols, vals);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the associated cs_matrix_structure_t according
 *         to the space scheme
 *
 * \param[in]  scheme       enum on the discretization scheme used
 *
 * \return  a pointer on a cs_matrix_structure_t *
 */
/*----------------------------------------------------------------------------*/

const cs_matrix_structure_t *
cs_equation_get_matrix_structure(cs_space_scheme_t   scheme)
{
  cs_matrix_structure_t  *ms = NULL;

  switch (scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    ms = cs_equation_common_ms[CS_EQ_COMMON_VERTEX];
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO:
    ms = cs_equation_common_ms[CS_EQ_COMMON_FACE];
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "Invalid space scheme.");
    break;

  } // Switch on space scheme

  return ms;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to the associated cs_matrix_assembler_t according
 *         to the space scheme
 *
 * \param[in]  scheme       enum on the discretization scheme used
 *
 * \return  a pointer on a cs_matrix_assembler_t *
 */
/*----------------------------------------------------------------------------*/

const cs_matrix_assembler_t *
cs_equation_get_matrix_assembler(cs_space_scheme_t   scheme)
{
  cs_matrix_assembler_t  *mav = NULL;

  switch (scheme) {

  case CS_SPACE_SCHEME_CDOVB:
  case CS_SPACE_SCHEME_CDOVCB:
    mav = cs_equation_common_ma[CS_EQ_COMMON_VERTEX];
    break;

  case CS_SPACE_SCHEME_CDOFB:
  case CS_SPACE_SCHEME_HHO:
    mav = cs_equation_common_ma[CS_EQ_COMMON_FACE];
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "Invalid space scheme.");
    break;

  } // Switch on space scheme

  return mav;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the connectivity vertex->vertices for the local rank
 *
 * \return  a pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_connect_index_t *
cs_equation_get_v2v_index(void)
{
  return cs_connect_v2v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the connectivity face->faces for the local rank
 *
 * \return  a pointer to a cs_connect_index_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_connect_index_t *
cs_equation_get_f2f_index(void)
{
  return cs_connect_f2f;
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
cs_equation_get_tmpbuf(void)
{
  return cs_equation_common_work_buffer;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the allocation size of the temporary buffer
 *
 * \return  the size of the temporary buffer
 */
/*----------------------------------------------------------------------------*/

size_t
cs_equation_get_tmpbuf_size(void)
{
  return cs_equation_common_work_buffer_size;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
