/* ===========================================================================
 * Routines to handle low-level routines related to CDO local quantities:
 * - local matrices (stored in dense format),
 * - local mesh structure related to a cell or to a couple cell/face
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <limits.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_math.h"
#include "cs_param.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_local.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_cdo_local.c

  \brief  Routines to handle low-level actions related to CDO local quantities:
  - local mesh structure related to a cell or to a couple cell/face
  - local systems

*/

/*=============================================================================
 * Local type definitions
 *============================================================================*/

#define CS_CDO_LOCAL_DBG       0

/*============================================================================
 * Global variables
 *============================================================================*/

cs_cell_mesh_t        **cs_cdo_local_cell_meshes = NULL;
cs_face_mesh_t        **cs_cdo_local_face_meshes = NULL;
cs_face_mesh_light_t  **cs_cdo_local_face_meshes_light = NULL;

/*============================================================================
 * Local static variables
 *============================================================================*/

static const int  n_robin_parameters = 3;
static int  cs_cdo_local_n_structures = 0;

/* Store predefined flags */
static const cs_flag_t  cs_cdo_local_flag_v =
  CS_FLAG_COMP_PV | CS_FLAG_COMP_PVQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_FV;
static const cs_flag_t  cs_cdo_local_flag_e =
  CS_FLAG_COMP_PE | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ | CS_FLAG_COMP_EV |
  CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EF  | CS_FLAG_COMP_EFQ;
static const cs_flag_t  cs_cdo_local_flag_peq =
  CS_FLAG_COMP_PEQ | CS_FLAG_COMP_FEQ;
static const cs_flag_t  cs_cdo_local_flag_f =
  CS_FLAG_COMP_PF  | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_FE  |
  CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EF  | CS_FLAG_COMP_EFQ | CS_FLAG_COMP_HFQ |
  CS_FLAG_COMP_FV;
static const cs_flag_t  cs_cdo_local_flag_pfq =
  CS_FLAG_COMP_PFQ | CS_FLAG_COMP_HFQ | CS_FLAG_COMP_FEQ;
static const cs_flag_t  cs_cdo_local_flag_deq =
  CS_FLAG_COMP_HFQ | CS_FLAG_COMP_DEQ;
static const cs_flag_t  cs_cdo_local_flag_fe =
  CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EF | CS_FLAG_COMP_EFQ;
static const cs_flag_t  cs_cdo_local_flag_ef =
  CS_FLAG_COMP_EF | CS_FLAG_COMP_EFQ;

/* Auxiliary buffers for computing quantities related to a cs_cell_mesh_t */
static double     **cs_cdo_local_dbuf = NULL;
static short int  **cs_cdo_local_kbuf = NULL;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate global structures related to a cs_cell_mesh_t and
 *         cs_face_mesh_t structures
 *
 * \param[in]   connect   pointer to a \ref cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_local_initialize(const cs_cdo_connect_t     *connect)
{
  /* Sanity check */
  assert(cs_glob_n_threads > 0);

  int  n_vc = connect->n_max_vbyc;
  int  size = cs_glob_n_threads;

  cs_cdo_local_n_structures = size;
  BFT_MALLOC(cs_cdo_local_cell_meshes, size, cs_cell_mesh_t *);
  BFT_MALLOC(cs_cdo_local_face_meshes, size, cs_face_mesh_t *);
  BFT_MALLOC(cs_cdo_local_face_meshes_light, size, cs_face_mesh_light_t *);
  BFT_MALLOC(cs_cdo_local_dbuf, size, double *);
  BFT_MALLOC(cs_cdo_local_kbuf, size, short int *);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cdo_local_cell_meshes[t_id] = cs_cell_mesh_create(connect);
    cs_cdo_local_face_meshes[t_id] = cs_face_mesh_create(connect->n_max_vbyf);
    cs_cdo_local_face_meshes_light[t_id] =
      cs_face_mesh_light_create(connect->n_max_vbyf, connect->n_max_vbyc);

    BFT_MALLOC(cs_cdo_local_dbuf[t_id], n_vc*(n_vc+1)/2, double);
    BFT_MALLOC(cs_cdo_local_kbuf[t_id], CS_MAX(connect->v_max_cell_range,
                                               connect->e_max_cell_range) + 1,
               short int);
  }
#else

  assert(cs_glob_n_threads == 1);

  cs_cdo_local_cell_meshes[0] = cs_cell_mesh_create(connect);
  cs_cdo_local_face_meshes[0] = cs_face_mesh_create(connect->n_max_vbyf);
  cs_cdo_local_face_meshes_light[0] =
    cs_face_mesh_light_create(connect->n_max_vbyf, connect->n_max_vbyc);

  BFT_MALLOC(cs_cdo_local_dbuf[0], n_vc*(n_vc+1)/2, double);
  BFT_MALLOC(cs_cdo_local_kbuf[0], CS_MAX(connect->v_max_cell_range,
                                          connect->e_max_cell_range) + 1,
             short int);

#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free global structures related to \ref cs_cell_mesh_t and
 *         \ref cs_face_mesh_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_local_finalize(void)
{
  if (cs_cdo_local_n_structures < 1)
    return;

  assert(cs_cdo_local_n_structures == cs_glob_n_threads);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_mesh_free(&(cs_cdo_local_cell_meshes[t_id]));
    cs_face_mesh_free(&(cs_cdo_local_face_meshes[t_id]));
    cs_face_mesh_light_free(&(cs_cdo_local_face_meshes_light[t_id]));
    BFT_FREE(cs_cdo_local_dbuf[t_id]);
    BFT_FREE(cs_cdo_local_kbuf[t_id]);

  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_mesh_free(&(cs_cdo_local_cell_meshes[0]));
  cs_face_mesh_free(&(cs_cdo_local_face_meshes[0]));
  cs_face_mesh_light_free(&(cs_cdo_local_face_meshes_light[0]));
  BFT_FREE(cs_cdo_local_dbuf[0]);
  BFT_FREE(cs_cdo_local_kbuf[0]);
#endif /* openMP */

  BFT_FREE(cs_cdo_local_cell_meshes);
  BFT_FREE(cs_cdo_local_face_meshes);
  BFT_FREE(cs_cdo_local_face_meshes_light);
  BFT_FREE(cs_cdo_local_dbuf);
  BFT_FREE(cs_cdo_local_kbuf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a \ref cs_cell_sys_t structure
 *
 * \param[in]   n_max_dofbyc    max number of entries
 * \param[in]   n_max_fbyc      max number of faces in a cell
 * \param[in]   n_blocks        number of blocks in a row/column
 * \param[in]   block_sizes     size of each block or NULL.
 *                              Specific treatment n_blocks = 1.
 *
 * \return a pointer to a new allocated \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_sys_t *
cs_cell_sys_create(int      n_max_dofbyc,
                   int      n_max_fbyc,
                   int      n_blocks,
                   int     *block_sizes)
{
  cs_cell_sys_t  *csys = NULL;

  BFT_MALLOC(csys, 1, cs_cell_sys_t);

  /* Metadata about DoFs */
  csys->cell_flag = 0;
  csys->c_id = -1;
  csys->n_dofs = 0;
  csys->dof_ids = NULL;
  csys->dof_flag = NULL;

  /* System and previous values */
  csys->mat = NULL;
  csys->rhs = NULL;
  csys->source = NULL;
  csys->val_n = NULL;

  /* Internal enforcement */
  csys->has_internal_enforcement = false;
  csys->intern_forced_ids = NULL;

  if (n_max_dofbyc > 0)
    BFT_MALLOC(csys->intern_forced_ids, n_max_dofbyc, cs_lnum_t);

  /* Boundary conditions */
  csys->n_bc_faces = 0;
  csys->_f_ids = NULL;
  csys->bf_ids = NULL;
  csys->bf_flag = NULL;

  csys->has_dirichlet = false;
  csys->dir_values = NULL;

  csys->has_nhmg_neumann = false;
  csys->neu_values = NULL;

  csys->has_robin = false;
  csys->rob_values = NULL;

  csys->has_sliding = false;

  if (n_max_fbyc > 0) {

    BFT_MALLOC(csys->bf_flag, n_max_fbyc, cs_flag_t);
    memset(csys->bf_flag, 0, sizeof(cs_flag_t)*n_max_fbyc);

    BFT_MALLOC(csys->_f_ids, n_max_fbyc, short int);
    memset(csys->_f_ids, 0, sizeof(short int)*n_max_fbyc);

    BFT_MALLOC(csys->bf_ids, n_max_fbyc, cs_lnum_t);
    memset(csys->bf_ids, 0, sizeof(cs_lnum_t)*n_max_fbyc);

  }

  if (n_max_dofbyc > 0) {

    BFT_MALLOC(csys->dof_flag, n_max_dofbyc, cs_flag_t);
    memset(csys->dof_flag, 0, sizeof(cs_flag_t)*n_max_dofbyc);

    BFT_MALLOC(csys->dof_ids, n_max_dofbyc, cs_lnum_t);
    memset(csys->dof_ids, 0, sizeof(cs_lnum_t)*n_max_dofbyc);

    if (block_sizes == NULL)
      csys->mat = cs_sdm_square_create(n_max_dofbyc);

    else if (n_blocks == 1) {

      assert(block_sizes != NULL);
      if (block_sizes[0] == 3) {
        int  n_row_blocks = n_max_dofbyc/3;
        assert(n_max_dofbyc % 3 == 0);
        csys->mat = cs_sdm_block33_create(n_row_blocks, n_row_blocks);
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid initialization of the cellwise block matrix\n",
                  __func__);

    }
    else
      csys->mat = cs_sdm_block_create(n_blocks, n_blocks,
                                      block_sizes,
                                      block_sizes);

    BFT_MALLOC(csys->rhs       , n_max_dofbyc, double);
    BFT_MALLOC(csys->source    , n_max_dofbyc, double);
    BFT_MALLOC(csys->val_n     , n_max_dofbyc, double);
    BFT_MALLOC(csys->dir_values, n_max_dofbyc, double);
    BFT_MALLOC(csys->neu_values, n_max_dofbyc, double);

    const size_t  s = n_max_dofbyc * sizeof(double);

    memset(csys->rhs       , 0, s);
    memset(csys->source    , 0, s);
    memset(csys->val_n     , 0, s);
    memset(csys->dir_values, 0, s);
    memset(csys->neu_values, 0, s);
  }

  int  n_rob_size = n_robin_parameters*CS_MAX(n_max_dofbyc, n_max_fbyc);
  BFT_MALLOC(csys->rob_values, n_rob_size, double);
  memset(csys->rob_values, 0, n_rob_size*sizeof(cs_real_t));

  return csys;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset all members related to BC and some other ones in a
 *         \ref cs_cell_sys_t structure
 *
 * \param[in]      n_fbyc     number of faces in a cell
 * \param[in, out] csys       pointer to the \ref cs_cell_sys_t struct to reset
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_sys_reset(int              n_fbyc,
                  cs_cell_sys_t   *csys)
{
  if (n_fbyc < 1 || csys->n_dofs < 1)
    return;

  const size_t  s = csys->n_dofs * sizeof(double);

  memset(csys->rhs, 0, s);
  memset(csys->source, 0, s);

  csys->n_bc_faces = 0;
  csys->has_dirichlet = csys->has_nhmg_neumann = csys->has_robin = false;

  csys->has_internal_enforcement = false;
  for (int i = 0; i < csys->n_dofs; i++)
    csys->intern_forced_ids[i] = -1; /* Not selected */

#if defined(DEBUG) && !defined(NDEBUG)
  memset(csys->bf_flag , 0, sizeof(cs_flag_t)*n_fbyc);
  memset(csys->_f_ids  , 0, sizeof(short int)*n_fbyc);
  memset(csys->bf_ids  , 0, sizeof(cs_lnum_t)*n_fbyc);
  memset(csys->dof_flag, 0, sizeof(cs_flag_t)*csys->n_dofs);

  memset(csys->dir_values, 0, s);
  memset(csys->neu_values, 0, s);
  memset(csys->rob_values, 0,
         CS_MAX(n_fbyc, csys->n_dofs)*sizeof(double)*n_robin_parameters);
#else
  if ((csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_VERTEX) ||
      (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE)) {

    memset(csys->bf_flag , 0, sizeof(cs_flag_t)*n_fbyc);
    memset(csys->_f_ids  , 0, sizeof(short int)*n_fbyc);
    memset(csys->bf_ids  , 0, sizeof(cs_lnum_t)*n_fbyc);
    memset(csys->dof_flag, 0, sizeof(cs_flag_t)*csys->n_dofs);

    memset(csys->dir_values, 0, s);
    memset(csys->neu_values, 0, s);
    memset(csys->rob_values, 0,
           CS_MAX(n_fbyc, csys->n_dofs)*sizeof(double)*n_robin_parameters);

  } /* Boundary cell -> reset BC-related members */
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_cell_sys_t structure
 *
 * \param[in, out]  p_csys   pointer of pointer to a \ref cs_cell_sys_t struct
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_sys_free(cs_cell_sys_t     **p_csys)
{
  cs_cell_sys_t  *csys = *p_csys;

  if (csys == NULL)
    return;

  BFT_FREE(csys->dof_ids);
  BFT_FREE(csys->dof_flag);

  csys->mat = cs_sdm_free(csys->mat);

  BFT_FREE(csys->rhs);
  BFT_FREE(csys->source);
  BFT_FREE(csys->val_n);

  BFT_FREE(csys->_f_ids);
  BFT_FREE(csys->bf_ids);
  BFT_FREE(csys->bf_flag);
  BFT_FREE(csys->dir_values);
  BFT_FREE(csys->neu_values);
  BFT_FREE(csys->rob_values);

  BFT_FREE(csys->intern_forced_ids);

  BFT_FREE(csys);
  *p_csys= NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a local system for debugging purpose
 *
 * \param[in]       msg     associated message to print
 * \param[in]       csys    pointer to a \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_sys_dump(const char             msg[],
                 const cs_cell_sys_t   *csys)
{
# pragma omp critical
  {
    bft_printf( "[rank:%d] %s\n", cs_glob_rank_id, msg);

    if (csys->cell_flag > 0) {
      bft_printf(">> dirichlet: %s, nhmg_neumann: %s, robin: %s,"
                 " sliding: %s\n",
                 cs_base_strtf(csys->has_dirichlet),
                 cs_base_strtf(csys->has_nhmg_neumann),
                 cs_base_strtf(csys->has_robin),
                 cs_base_strtf(csys->has_sliding));
      bft_printf(">> Boundary faces\n"
                 ">> %-8s | %-8s | %-6s\n", "_ID", "ID", "FLAG");
      for (int i = 0; i < csys->n_bc_faces; i++) {
        short int f = csys->_f_ids[i];
        bft_printf(">> %8d | %8d | %6d\n",
                   f, csys->bf_ids[f], csys->bf_flag[f]);
      }
    }

    if (csys->mat->flag & CS_SDM_BY_BLOCK)
      cs_sdm_block_dump(csys->c_id, csys->mat);
    else
      cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, csys->mat);

    bft_printf(">> %-8s | %-10s | %-10s | %-10s | %-8s |"
               " %-6s | %-10s\n",
               "IDS", "RHS", "TS", "VAL_PREV", "ENFORCED", "FLAG",
               "DIR_VALS");
    for (int i = 0; i < csys->n_dofs; i++)
      bft_printf(">> %8d | % -.3e | % -.3e | % -.3e |"
                 " %8d | %6d | % -.3e\n",
                 csys->dof_ids[i], csys->rhs[i], csys->source[i],
                 csys->val_n[i], csys->intern_forced_ids[i],
                 csys->dof_flag[i], csys->dir_values[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate \ref cs_cell_builder_t structure
 *
 * \return a pointer to the new allocated \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_builder_t *
cs_cell_builder_create(void)
{
  cs_cell_builder_t  *cb = NULL;

  /* Common part to all discretization */
  BFT_MALLOC(cb, 1, cs_cell_builder_t);

  /* Consider a unitary diffusion property by default */
  cb->eig_ratio = 1;
  cb->eig_max = 1;

  cb->dpty_mat[0][0] = cb->dpty_mat[1][1] = cb->dpty_mat[2][2] = 1;
  cb->dpty_mat[0][1] = cb->dpty_mat[1][0] = cb->dpty_mat[2][0] = 0;
  cb->dpty_mat[0][2] = cb->dpty_mat[1][2] = cb->dpty_mat[2][1] = 0;
  cb->dpty_val = 1;
  cb->cpty_mat[0][0] = cb->cpty_mat[1][1] = cb->cpty_mat[2][2] = 1;
  cb->cpty_mat[0][1] = cb->cpty_mat[1][0] = cb->cpty_mat[2][0] = 0;
  cb->cpty_mat[0][2] = cb->cpty_mat[1][2] = cb->cpty_mat[2][1] = 0;
  cb->cpty_val = 1;
  cb->gpty_val = 1;
  cb->tpty_val = 1;
  cb->rpty_val = 1;

  for (int r = 0; r < CS_CDO_N_MAX_REACTIONS; r++) cb->rpty_vals[r] = 1;

  cb->adv_fluxes = NULL;
  cb->ids = NULL;
  cb->values = NULL;
  cb->vectors = NULL;

  cb->hdg = cb->loc = cb->aux = NULL;

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_cell_builder_t structure
 *
 * \param[in, out]  p_cb  pointer of pointer to a \ref cs_cell_builder_t struct
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_builder_free(cs_cell_builder_t     **p_cb)
{
  cs_cell_builder_t  *cb = *p_cb;

  if (cb == NULL)
    return;

  BFT_FREE(cb->adv_fluxes);
  BFT_FREE(cb->ids);
  BFT_FREE(cb->values);
  BFT_FREE(cb->vectors);

  cb->hdg = cs_sdm_free(cb->hdg);
  cb->loc = cs_sdm_free(cb->loc);
  cb->aux = cs_sdm_free(cb->aux);

  BFT_FREE(cb);
  *p_cb = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a cs_cell_mesh_t structure
 *
 * \param[in]  connect        pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_mesh_t *
cs_cell_mesh_create(const cs_cdo_connect_t   *connect)
{
  cs_cell_mesh_t  *cm = NULL;

  BFT_MALLOC(cm, 1, cs_cell_mesh_t);

  /* Sizes used to allocate buffers */
  cm->n_max_vbyc = connect->n_max_vbyc;
  cm->n_max_ebyc = connect->n_max_ebyc;
  cm->n_max_fbyc = connect->n_max_fbyc;

  cm->flag = 0;
  cm->n_vc = 0;
  cm->n_ec = 0;
  cm->n_fc = 0;

  /* Vertex information */
  BFT_MALLOC(cm->v_ids, cm->n_max_vbyc, cs_lnum_t);
  BFT_MALLOC(cm->wvc, cm->n_max_vbyc, double);
  BFT_MALLOC(cm->xv, 3*cm->n_max_vbyc, double);

  /* Edge information */
  BFT_MALLOC(cm->e_ids, cm->n_max_ebyc, cs_lnum_t);
  BFT_MALLOC(cm->edge, cm->n_max_ebyc, cs_quant_t);
  BFT_MALLOC(cm->dface, cm->n_max_ebyc, cs_nvec3_t);
  BFT_MALLOC(cm->e2v_sgn, cm->n_max_ebyc, short int);

  /* Face information */
  BFT_MALLOC(cm->f_ids, cm->n_max_fbyc, cs_lnum_t);
  BFT_MALLOC(cm->f_sgn, cm->n_max_fbyc, short int);
  BFT_MALLOC(cm->f_diam, cm->n_max_fbyc, double);
  BFT_MALLOC(cm->hfc, cm->n_max_fbyc, double);
  BFT_MALLOC(cm->pfc, cm->n_max_fbyc, double);
  BFT_MALLOC(cm->face, cm->n_max_fbyc, cs_quant_t);
  BFT_MALLOC(cm->dedge, cm->n_max_fbyc, cs_nvec3_t);

  /* face --> vertices connectivity */
  BFT_MALLOC(cm->f2v_idx, cm->n_max_fbyc + 1, short int);
  BFT_MALLOC(cm->f2v_ids, 2*cm->n_max_ebyc, short int);

  /* face --> edges connectivity */
  BFT_MALLOC(cm->f2e_idx, cm->n_max_fbyc + 1, short int);
  BFT_MALLOC(cm->f2e_ids, 2*cm->n_max_ebyc, short int);
  BFT_MALLOC(cm->tef, 2*cm->n_max_ebyc, double);

  /* edge --> vertices connectivity */
  BFT_MALLOC(cm->e2v_ids, 2*cm->n_max_ebyc, short int);

  /* edge --> face connectivity */
  BFT_MALLOC(cm->e2f_ids, 2*cm->n_max_ebyc, short int);
  BFT_MALLOC(cm->sefc, 2*cm->n_max_ebyc, cs_nvec3_t);

  cs_cell_mesh_reset(cm);

  return cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_cell_mesh_t structure corresponding to mesh id
 *
 * \param[in]   mesh_id   id in the array of pointer to cs_cell_mesh_t struct.
 *
 * \return a pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_mesh_t *
cs_cdo_local_get_cell_mesh(int    mesh_id)
{
  if (mesh_id < 0 || mesh_id >= cs_glob_n_threads)
    return NULL;

  cs_cell_mesh_t  *cm = cs_cdo_local_cell_meshes[mesh_id];

#if defined(DEBUG) && !defined(NDEBUG)
  /* This is to check that the mesh flag is correctly set */
  cs_cell_mesh_reset(cm);
#endif

  return cm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize to invalid values a cs_cell_mesh_t structure
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_reset(cs_cell_mesh_t   *cm)
{
  cm->n_vc = -1;
  cm->n_ec = -1;
  cm->n_fc = -1;

  /* Cell information */
  cm->c_id = -1;
  cm->xc[0] = cm->xc[1] = cm->xc[2] = -DBL_MAX;
  cm->vol_c = -DBL_MAX;
  cm->diam_c = -DBL_MAX;

  /* Vertex information */
  for (short int v = 0; v < cm->n_max_vbyc; v++) {
    cm->v_ids[v] = -1;
    cm->wvc[v] = -DBL_MAX;
    cm->xv[3*v] = cm->xv[3*v+1] = cm->xv[3*v+2] = -DBL_MAX;
  }

  /* Edge information */
  for (short int e = 0; e < cm->n_max_ebyc; e++) {
    cm->e_ids[e] = -1;
    cm->e2v_sgn[e] = 0;
    cm->edge[e].meas = cm->dface[e].meas = -DBL_MAX;
    cm->edge[e].unitv[0] = cm->dface[e].unitv[0] = -DBL_MAX;
    cm->edge[e].unitv[1] = cm->dface[e].unitv[1] = -DBL_MAX;
    cm->edge[e].unitv[2] = cm->dface[e].unitv[2] = -DBL_MAX;
    cm->edge[e].center[0] = -DBL_MAX;
    cm->edge[e].center[1] = -DBL_MAX;
    cm->edge[e].center[2] = -DBL_MAX;
  }

  /* Face information */
  for (short int f = 0; f < cm->n_max_fbyc; f++) {
    cm->f_ids[f] = -1;
    cm->f_sgn[f] = 0;
    cm->f_diam[f] = -DBL_MAX;
    cm->hfc[f] = -DBL_MAX;
    cm->pfc[f] = -DBL_MAX;
    cm->face[f].meas = cm->dedge[f].meas = -DBL_MAX;
    cm->face[f].unitv[0] = cm->dedge[f].unitv[0] = -DBL_MAX;
    cm->face[f].unitv[1] = cm->dedge[f].unitv[1] = -DBL_MAX;
    cm->face[f].unitv[2] = cm->dedge[f].unitv[2] = -DBL_MAX;
    cm->face[f].center[0] = -DBL_MAX;
    cm->face[f].center[1] = -DBL_MAX;
    cm->face[f].center[2] = -DBL_MAX;
  }

  /* face --> edges connectivity */
  for (short int f = 0; f < cm->n_max_fbyc + 1; f++)
    cm->f2e_idx[f] = cm->f2v_idx[f] = -1;

  for (int i = 0; i < 2*cm->n_max_ebyc; i++) {
    cm->e2v_ids[i] = cm->e2f_ids[i] = cm->f2e_ids[i] = cm->f2v_ids[i] = -1;
    cm->tef[i] = cm->sefc[i].meas = -DBL_MAX;
    cm->sefc[i].unitv[0]=cm->sefc[i].unitv[1]=cm->sefc[i].unitv[2] = -DBL_MAX;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_cell_mesh_t structure
 *
 * \param[in]    cm    pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_dump(const cs_cell_mesh_t     *cm)
{
  if (cm == NULL) {
    bft_printf("\n>> Dump cs_cell_mesh_t %p\n", (const void *)cm);
    return;
  }

  bft_printf("\n>> [rank: %d] Dump cs_cell_mesh_t %p; %s; flag: %d\n"
             " c_id:%d; vol: %9.6e; xc (% .4e % .4e % .4e); diam: % .4e\n",
             cs_glob_rank_id, (const void *)cm, fvm_element_type_name[cm->type],
             cm->flag, cm->c_id, cm->vol_c, cm->xc[0], cm->xc[1], cm->xc[2],
             cm->diam_c);

  /* Information related to primal vertices */
  if (cm->flag & cs_cdo_local_flag_v) {

    bft_printf(" %s | %6s | %35s | %10s\n",
               "v", "id", "coord", "wvc");
    for (short int v = 0; v < cm->n_vc; v++)
      bft_printf("%2d | %6d | % .4e % .4e % .4e | %.4e\n",
                 v, cm->v_ids[v], cm->xv[3*v], cm->xv[3*v+1], cm->xv[3*v+2],
                 cm->wvc[v]);

  } /* Vertex quantities */

  /* Information related to primal edges */
  if (cm->flag & cs_cdo_local_flag_e) {

    bft_printf(" %s | %6s | %3s | %2s | %2s | %9s |"
               " %35s | %35s | %10s | %35s\n",
               "e", "id", "sgn", "v1", "v2", "length", "unit", "coords",
               "df.meas", "df.unit");
    for (short int e = 0; e < cm->n_ec; e++) {

      cs_quant_t  peq = cm->edge[e];
      cs_nvec3_t  dfq = cm->dface[e];
      bft_printf("%2d | %6d | %3d | %2d | %2d | %.3e |"
                 " % .4e % .4e % .4e | % .4e % .4e % .4e | %.4e |"
                 " % .4e % .4e % .4e\n",
                 e, cm->e_ids[e], cm->e2v_sgn[e], cm->e2v_ids[2*e],
                 cm->e2v_ids[2*e+1], peq.meas, peq.unitv[0], peq.unitv[1],
                 peq.unitv[2], peq.center[0], peq.center[1], peq.center[2],
                 dfq.meas, dfq.unitv[0], dfq.unitv[1], dfq.unitv[2]);

    } /* Loop on edges */

  } /* Edge quantities */

  /* Information related to primal faces */
  if (cm->flag & cs_cdo_local_flag_f) {

    bft_printf(" %s | %6s | %9s | %3s | %35s | %35s |"
               " %10s | %35s | %11s  %11s  %11s\n",
               "f", "id", "surf", "sgn", "unit", "coords", "dlen", "dunitv",
               "pfc",  "hfc", "diam");
    for (short int f = 0; f < cm->n_fc; f++) {
      cs_quant_t  pfq = cm->face[f];
      cs_nvec3_t  deq = cm->dedge[f];
      bft_printf("%2d | %6d | %.3e | %3d | % .4e % .4e % .4e |"
                 " % .4e % .4e % .4e | %.4e | % .4e % .4e % .4e | %.3e |"
                 " %.3e | %.3e\n",
                 f, cm->f_ids[f], pfq.meas, cm->f_sgn[f],
                 pfq.unitv[0], pfq.unitv[1], pfq.unitv[2], pfq.center[0],
                 pfq.center[1], pfq.center[2], deq.meas, deq.unitv[0],
                 deq.unitv[1], deq.unitv[2], cm->pfc[f], cm->hfc[f],
                 cm->f_diam[f]);
    }

  } /* Face quantities */

  if (cm->flag & cs_cdo_local_flag_fe) {

    bft_printf(" n_ef | f: pef\n");
    for (short int f = 0; f < cm->n_fc; f++) {
      bft_printf(" %4d |",
                 cm->f2e_idx[f+1] - cm->f2e_idx[f]);
      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++)
        bft_printf(" %2d:%.4e|", cm->f2e_ids[i], cm->tef[i]);
      bft_printf("\n");
    }

  }

  if (cm->flag & cs_cdo_local_flag_ef) {

    bft_printf("%-4s | f0 | %-53s | f1 | %-53s\n",
               "e", "sef0c: meas, unitv", "sef1c: meas, unitv");
    for (short int e = 0; e < cm->n_ec; e++)
      bft_printf(" %3d | %2d | % .4e (% .4e % .4e % .4e) |"
                 " %2d | % .4e (% .4e % .4e % .4e)\n",
                 e, cm->e2f_ids[2*e], cm->sefc[2*e].meas,
                 cm->sefc[2*e].unitv[0], cm->sefc[2*e].unitv[1],
                 cm->sefc[2*e].unitv[2], cm->e2f_ids[2*e+1],
                 cm->sefc[2*e+1].meas, cm->sefc[2*e+1].unitv[0],
                 cm->sefc[2*e+1].unitv[1], cm->sefc[2*e+1].unitv[2]);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cell_mesh_t structure
 *
 * \param[in, out]  p_cm   pointer of pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_free(cs_cell_mesh_t     **p_cm)
{
  cs_cell_mesh_t  *cm = *p_cm;

  if (cm == NULL)
    return;

  BFT_FREE(cm->v_ids);
  BFT_FREE(cm->wvc);
  BFT_FREE(cm->xv);

  BFT_FREE(cm->e_ids);
  BFT_FREE(cm->edge);
  BFT_FREE(cm->dface);

  BFT_FREE(cm->f_ids);
  BFT_FREE(cm->f_sgn);
  BFT_FREE(cm->f_diam);
  BFT_FREE(cm->hfc);
  BFT_FREE(cm->pfc);
  BFT_FREE(cm->face);
  BFT_FREE(cm->dedge);

  BFT_FREE(cm->e2v_ids);
  BFT_FREE(cm->e2v_sgn);

  BFT_FREE(cm->f2v_idx);
  BFT_FREE(cm->f2v_ids);

  BFT_FREE(cm->f2e_idx);
  BFT_FREE(cm->f2e_ids);
  BFT_FREE(cm->tef);

  BFT_FREE(cm->e2f_ids);
  BFT_FREE(cm->sefc);

  BFT_FREE(cm);
  *p_cm = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_cell_mesh_t structure for a given cell id. According
 *         to the requested level, some quantities may not be defined;
 *
 * \param[in]       c_id        cell id
 * \param[in]       build_flag  indicate which members are really built
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       quant       pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  cm          pointer to a cs_cell_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_build(cs_lnum_t                    c_id,
                   cs_flag_t                    build_flag,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_cell_mesh_t              *cm)
{
  if (cm == NULL)
    return;

  cm->flag = build_flag;
  cm->type = connect->cell_type[c_id];

  /* Store information related to cell */
  cm->c_id = c_id;
  cm->vol_c = quant->cell_vol[c_id];
  for (int k = 0; k < 3; k++)
    cm->xc[k] = quant->cell_centers[3*c_id+k];

  /* Store the number of cell faces (useful to allocated boundary quantities) */
  const cs_lnum_t  *c2f_idx = connect->c2f->idx + c_id;

  cm->n_fc = c2f_idx[1] - c2f_idx[0];
  cm->bface_shift = quant->n_i_faces; /* Border faces come after interior
                                       * faces */

  if (build_flag == 0)
    return;

  /* Information related to primal vertices */
  if (build_flag & cs_cdo_local_flag_v) {

    const cs_lnum_t  *c2v_idx = connect->c2v->idx + c_id;
    const cs_lnum_t  *c2v_ids = connect->c2v->ids + c2v_idx[0];

    cm->n_vc = c2v_idx[1] - c2v_idx[0];

    for (short int v = 0; v < cm->n_vc; v++) {
      const cs_lnum_t  v_id = c2v_ids[v];
      cm->v_ids[v] = v_id;
      for (int k = 0; k < 3; k++)
        cm->xv[3*v+k] = quant->vtx_coord[3*v_id+k];
    }

    /* Primal vertices quantities */
    if (build_flag & CS_FLAG_COMP_PVQ) {

      const double  *wvc = quant->dcell_vol + c2v_idx[0];
      const double  invvol = 1/cm->vol_c;
      for (short int v = 0; v < cm->n_vc; v++)
        cm->wvc[v] = invvol * wvc[v];

    }

  } /* vertices */

  /* Information related to primal edges */
  if (build_flag & cs_cdo_local_flag_e) {

    const cs_lnum_t  *c2e_idx = connect->c2e->idx + c_id;
    const cs_lnum_t  *c2e_ids = connect->c2e->ids + c2e_idx[0];

    cm->n_ec = c2e_idx[1] - c2e_idx[0];

    for (short int e = 0; e < cm->n_ec; e++)
      cm->e_ids[e] = c2e_ids[e];

    if (build_flag & cs_cdo_local_flag_peq) {

      assert(build_flag & CS_FLAG_COMP_PV);

      /* Primal edge quantities */
      for (short int e = 0; e < cm->n_ec; e++) {

        const cs_lnum_t  e_id = cm->e_ids[e];
        const cs_nvec3_t  nv = cs_quant_set_edge_nvec(e_id, quant);
        const cs_lnum_t  *v_id = connect->e2v->ids + 2*e_id;
        const cs_real_t  *xv1 = quant->vtx_coord + 3*v_id[0];
        const cs_real_t  *xv2 = quant->vtx_coord + 3*v_id[1];

        cm->edge[e].meas = nv.meas;
        for (int k = 0; k < 3; k++) {
          cm->edge[e].center[k] = 0.5 * (xv1[k] + xv2[k]);
          cm->edge[e].unitv[k] = nv.unitv[k];
        }

      }

    } /* Primal edge quantities */

    /* Dual face quantities related to each edge */
    if (build_flag & CS_FLAG_COMP_DFQ) {

      for (short int e = 0; e < cm->n_ec; e++) {

        const cs_real_t  *sface0 = quant->sface_normal + 6*(c2e_idx[0] + e);
        const cs_real_t  *sface1 = sface0 + 3;

        cs_real_3_t  df_vect;
        for (int k = 0; k < 3; k++)
          df_vect[k] = sface0[k] + sface1[k];

        cs_nvec3_t  df_nvect;
        cs_nvec3(df_vect, &df_nvect);

        cm->dface[e].meas = df_nvect.meas;
        for (int k = 0; k < 3; k++)
          cm->dface[e].unitv[k] = df_nvect.unitv[k];

      }

    } /* Dual face quantities */

  } /* Edge-related quantities */

  /* Information related to primal faces */
  if (build_flag & cs_cdo_local_flag_f) {

    const cs_lnum_t  *c2f_ids = connect->c2f->ids + c2f_idx[0];
    const short int  *c2f_sgn = connect->c2f->sgn + c2f_idx[0];

    for (short int f = 0; f < cm->n_fc; f++) {
      cm->f_ids[f] = c2f_ids[f];
      cm->f_sgn[f] = c2f_sgn[f];
    }

    /* Face related quantities */
    if (build_flag & cs_cdo_local_flag_pfq) {

      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_quant_t  pfq = cs_quant_set_face(cm->f_ids[f], quant);

        cm->face[f].meas = pfq.meas;
        for (int k = 0; k < 3; k++) {
          cm->face[f].center[k] = pfq.center[k];
          cm->face[f].unitv[k] = pfq.unitv[k];
        }

      }

    } /* Primal face quantities */

    if (build_flag & cs_cdo_local_flag_deq) {

      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_nvec3_t  de_nvect = cs_quant_set_dedge_nvec(c2f_idx[0]+f,
                                                             quant);

        /* Copy cs_nvec3_t structure */
        cm->dedge[f].meas = de_nvect.meas;
        for (int k = 0; k < 3; k++)
          cm->dedge[f].unitv[k] = de_nvect.unitv[k];

      }

    } /* Dual edge quantities */

    if (build_flag & CS_FLAG_COMP_HFQ) {

      /* Compute the height of the pyramid of base f whose apex is
         the cell center */
      for (short int f = 0; f < cm->n_fc; f++) {
        cm->hfc[f] = cs_math_3_dot_product(cm->face[f].unitv,
                                           cm->dedge[f].unitv);
        cm->hfc[f] *= cm->dedge[f].meas;
#if defined(DEBUG) && !defined(NDEBUG) && CS_CDO_LOCAL_DBG > 0
        if (cm->hfc[f] <= 0)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Invalid result; hfc = %5.3e < 0 !\n",
                    __func__, cm->hfc[f]);
#endif

        /* Volume of the pyramid of base f and apex x_c */
        cm->pfc[f] = cs_math_1ov3 * cm->hfc[f] * cm->face[f].meas;

      } /* Loop on cell faces */

    } /* Quantities related to the pyramid of base f */

  } /* Face information */

  if (build_flag & CS_FLAG_COMP_EV || build_flag & CS_FLAG_COMP_FV) {

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);
#else
    int t_id = 0;
#endif /* openMP */

    short int  *kbuf = cs_cdo_local_kbuf[t_id];

    /* Store in a compact way ids for vertices: mesh --> cell mesh */
    cs_lnum_t  v_shift = cm->v_ids[0];
    for (short int v = 1; v < cm->n_vc; v++)
      if (cm->v_ids[v] < v_shift)
        v_shift = cm->v_ids[v];
    for (short int v = 0; v < cm->n_vc; v++)
      kbuf[cm->v_ids[v]-v_shift] = v;

    if (build_flag & CS_FLAG_COMP_EV) {
      for (short int e = 0; e < cm->n_ec; e++) {

        const cs_lnum_t  e_id = cm->e_ids[e];

        /* Store only the sign related to the first vertex since the sign
           related to the second one is minus the first one */
        cm->e2v_sgn[e] = connect->e2v->sgn[2*e_id];
        cm->e2v_ids[2*e]   = kbuf[connect->e2v->ids[2*e_id] - v_shift];
        cm->e2v_ids[2*e+1] = kbuf[connect->e2v->ids[2*e_id+1] - v_shift];

      } /* Loop on cell edges */
    } /* edge-vertices information */

    if (build_flag & CS_FLAG_COMP_FV) {

      /* Build the index and the list of face vertices in the cellwise
       * numbering */
      cm->f2v_idx[0] = 0;
      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_lnum_t  f_id = cm->f_ids[f];
        if (f_id < cm->bface_shift) { /* Interior face */

          const cs_lnum_t  *idx = connect->if2v->idx + f_id;
          const cs_lnum_t  n_ent = idx[1] - idx[0];
          const cs_lnum_t  *v_ids = connect->if2v->ids + idx[0];

          cm->f2v_idx[f+1] = cm->f2v_idx[f] + n_ent;

          short int  *_ids = cm->f2v_ids + cm->f2v_idx[f];
          for (cs_lnum_t i = 0; i < n_ent; i++)
            _ids[i] = kbuf[v_ids[i] - v_shift];

        }
        else { /* Border face */

          const cs_lnum_t  bf_id = f_id - cm->bface_shift;
          const cs_lnum_t  *idx = connect->bf2v->idx + bf_id;
          const cs_lnum_t  n_ent = idx[1] - idx[0];
          const cs_lnum_t  *v_ids = connect->bf2v->ids + idx[0];

          cm->f2v_idx[f+1] = cm->f2v_idx[f] + n_ent;

          short int  *_ids = cm->f2v_ids + cm->f2v_idx[f];
          for (cs_lnum_t i = 0; i < n_ent; i++)
            _ids[i] = kbuf[v_ids[i] - v_shift];

        }

      } /* Loop on cell faces */

      assert(cm->f2v_idx[cm->n_fc] == 2*cm->n_ec);

    } /* face --> vertices connectivity */

  } /* EV or FV flag */

  if (build_flag & cs_cdo_local_flag_fe) {

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);
#else
    int t_id = 0;
#endif /* openMP */

    short int  *kbuf = cs_cdo_local_kbuf[t_id];

    /* Store in compact way: mesh --> cell mesh ids for edges */
    cs_lnum_t  shift = cm->e_ids[0];
    for (short int e = 1; e < cm->n_ec; e++)
      if (cm->e_ids[e] < shift)
        shift = cm->e_ids[e];
    for (short int e = 0; e < cm->n_ec; e++)
      kbuf[cm->e_ids[e]-shift] = e;

    const cs_adjacency_t  *f2e = connect->f2e;

    cm->f2e_idx[0] = 0;
    for (short int f = 0; f < cm->n_fc; f++) {

      const cs_lnum_t  f_id = cm->f_ids[f];
      const cs_lnum_t  *idx = f2e->idx + f_id;
      const cs_lnum_t  *ids = f2e->ids + idx[0];
      const cs_lnum_t  n_ef = idx[1] - idx[0];

      cm->f2e_idx[f+1] = cm->f2e_idx[f] + n_ef;
      short int  *_ids = cm->f2e_ids + cm->f2e_idx[f];
      for (cs_lnum_t i = 0; i < n_ef; i++)
        _ids[i] = kbuf[ids[i] - shift];

    } /* Loop on cell faces */

    /* Sanity check */
    assert(cm->f2e_idx[cm->n_fc] == 2*cm->n_ec);

    if (build_flag & CS_FLAG_COMP_FEQ) {

      for (short int f = 0; f < cm->n_fc; f++) {
        for (int ie = cm->f2e_idx[f]; ie < cm->f2e_idx[f+1]; ie++)
          cm->tef[ie] = cs_compute_area_from_quant(cm->edge[cm->f2e_ids[ie]],
                                                   cm->face[f].center);
      }

    } /* face --> edges quantities */

  } /* face --> edges connectivity */

  if (build_flag & cs_cdo_local_flag_ef) {

    /* Build the e2f connectivity */
    for (short int i = 0; i < 2*cm->n_ec; i++) cm->e2f_ids[i] = -1;

    for (short int f = 0; f < cm->n_fc; f++) {
      for (short int ie = cm->f2e_idx[f]; ie < cm->f2e_idx[f+1]; ie++) {

        short int  *ids = cm->e2f_ids + 2*cm->f2e_ids[ie];
        if (ids[0] == -1)
          ids[0] = f;
        else {
          assert(ids[1] == -1);
          ids[1] = f;
        }

      } /* Loop on face edges */
    } /* Loop on cell faces */

    if (build_flag & CS_FLAG_COMP_EFQ) { /* Build cm->sefc */

      cs_nvec3_t  nv;
      const cs_lnum_t  *c2e_idx = connect->c2e->idx + c_id;

      for (short int e = 0; e < cm->n_ec; e++) {

        const cs_real_t  *sface0 = quant->sface_normal + 6*(c2e_idx[0] + e);
        const cs_real_t  *sface1 = sface0 + 3;
        cs_nvec3_t  *_sefc = cm->sefc + 2*e;

        cs_nvec3(sface0, &nv);
        _sefc[0].meas = nv.meas;
        _sefc[0].unitv[0] = nv.unitv[0];
        _sefc[0].unitv[1] = nv.unitv[1];
        _sefc[0].unitv[2] = nv.unitv[2];
        cs_nvec3(sface1, &nv);
        _sefc[1].meas = nv.meas;
        _sefc[1].unitv[0] = nv.unitv[0];
        _sefc[1].unitv[1] = nv.unitv[1];
        _sefc[1].unitv[2] = nv.unitv[2];

      } /* Loop on face edges */

    } /* (edga,face) quantities */

  } /* edge-->faces */

  if (build_flag & CS_FLAG_COMP_DIAM) {

    assert(cs_flag_test(build_flag, CS_FLAG_COMP_EV | CS_FLAG_COMP_FE));

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);
#else
    int t_id = 0;
#endif /* openMP */

    double  *dbuf = cs_cdo_local_dbuf[t_id];
    short int  *vtag = cs_cdo_local_kbuf[t_id];
    int  size = cm->n_vc*(cm->n_vc-1)/2;
    int  shift = 0;

    /* Reset diam */
    cm->diam_c = -1;
    memset(dbuf, 0, sizeof(cs_real_t)*size);

    for (short int vi = 0; vi < cm->n_vc; vi++) {
      const double *xvi = cm->xv + 3*vi;
      for (short int vj = vi+1; vj < cm->n_vc; vj++) {
        double  l = cs_math_3_distance(xvi, cm->xv + 3*vj);
        dbuf[shift++] = l;
        if (l > cm->diam_c) cm->diam_c = l;

      } /* Loop on vj > vi */
    }   /* Loop on vi */

    for (short int f = 0; f < cm->n_fc; f++) {

      /* Reset vtag */
      for (short int v = 0; v < cm->n_vc; v++) vtag[v] = -1;

      /* Tag face vertices */
      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

        const int  eshft = 2*cm->f2e_ids[i];
        vtag[cm->e2v_ids[eshft  ]] = 1;
        vtag[cm->e2v_ids[eshft+1]] = 1;
      }

      cm->f_diam[f] = -1;
      shift = 0;
      for (short int vi = 0; vi < cm->n_vc; vi++) {
        for (short int vj = vi+1; vj < cm->n_vc; vj++) {

          if (vtag[vi] > 0) /* belong to the current face */
            if (vtag[vj] > 0) /* belong to the current face */
              if (dbuf[shift] > cm->f_diam[f]) cm->f_diam[f] = dbuf[shift];
          shift++;

        } /* Loop on vj > vi */
      }   /* Loop on vi */

    } /* Loop on cell faces */

  } /* Compute diameters */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_face_mesh_t structure
 *
 * \param[in]  n_max_vbyf    max. number of vertices fir a face
 *
 * \return a pointer to a new allocated cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_t *
cs_face_mesh_create(short int   n_max_vbyf)
{
  cs_face_mesh_t  *fm = NULL;

  BFT_MALLOC(fm, 1, cs_face_mesh_t);

  fm->n_max_vbyf = n_max_vbyf;

  fm->c_id = -1;
  fm->xc[0] = fm->xc[1] = fm->xc[2] = 0.;

  /* Face-related quantities */
  fm->f_id = -1;
  fm->f_sgn = 0;

  /* Vertex-related quantities */
  fm->n_vf = 0;
  BFT_MALLOC(fm->v_ids, fm->n_max_vbyf, cs_lnum_t);
  BFT_MALLOC(fm->xv, 3*fm->n_max_vbyf, double);
  BFT_MALLOC(fm->wvf, fm->n_max_vbyf, double);

  /* Edge-related quantities */
  fm->n_ef = 0;
  BFT_MALLOC(fm->e_ids, fm->n_max_vbyf, cs_lnum_t);
  BFT_MALLOC(fm->edge,  fm->n_max_vbyf, cs_quant_t);
  BFT_MALLOC(fm->e2v_ids, 2*fm->n_max_vbyf, short int);
  BFT_MALLOC(fm->tef, fm->n_max_vbyf, double);

  return fm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_face_mesh_t structure corresponding to mesh id
 *
 * \param[in]   mesh_id   id in the array of pointer to cs_face_mesh_t struct.
 *
 * \return a pointer to a cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_t *
cs_cdo_local_get_face_mesh(int    mesh_id)
{
  if (mesh_id < 0 || mesh_id >= cs_glob_n_threads)
    return NULL;

  return cs_cdo_local_face_meshes[mesh_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_face_mesh_t structure
 *
 * \param[in, out]  p_fm   pointer of pointer to a cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_free(cs_face_mesh_t     **p_fm)
{
  cs_face_mesh_t  *fm = *p_fm;

  if (fm == NULL)
    return;

  BFT_FREE(fm->v_ids);
  BFT_FREE(fm->xv);
  BFT_FREE(fm->wvf);

  BFT_FREE(fm->e_ids);
  BFT_FREE(fm->edge);
  BFT_FREE(fm->e2v_ids);
  BFT_FREE(fm->tef);

  BFT_FREE(fm);
  *p_fm = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_face_mesh_t structure for a given face/cell id.
 *
 * \param[in]       c_id      cell id
 * \param[in]       f_id      face id in the mesh structure
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_build(cs_lnum_t                    c_id,
                   cs_lnum_t                    f_id,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_face_mesh_t              *fm)
{
  if (fm == NULL)
    return;

  /* Sanity checks */
  assert(c_id > -1);
  assert(f_id > -1);

  fm->c_id = c_id;
  const cs_real_t  *xc = quant->cell_centers + 3*c_id;
  for (int k = 0; k < 3; k++) fm->xc[k] = xc[k];

  /* Face-related quantities */
  const cs_quant_t  pfq = cs_quant_set_face(f_id, quant);

  fm->f_id = f_id;
  fm->face.meas = pfq.meas;
  for (int k = 0; k < 3; k++) {
    fm->face.center[k] = pfq.center[k];
    fm->face.unitv[k] = pfq.unitv[k];
  }

  const cs_lnum_t  *c2f_idx = connect->c2f->idx + c_id;
  const cs_lnum_t  *c2f_ids = connect->c2f->ids + c2f_idx[0];
  const int  n_fc = c2f_idx[1] - c2f_idx[0];

  short int _f = n_fc;
  for (short int f = 0; f < n_fc; f++) {
    if (c2f_ids[f] == f_id) {

      const cs_nvec3_t  de_nvect = cs_quant_set_dedge_nvec(c2f_idx[0]+f,
                                                           quant);
      const short int  *f_sgn = connect->c2f->sgn + c2f_idx[0];

      _f = f;
      fm->dedge.meas = de_nvect.meas;
      for (int k = 0; k < 3; k++) fm->dedge.unitv[k] = de_nvect.unitv[k];
      fm->f_sgn = f_sgn[f];
      break;
    }
  }

  if (_f == n_fc) /* Sanity check */
    bft_error(__FILE__, __LINE__, 0,
              _(" Face %d not found.\n Stop build a face mesh."), f_id);

  const cs_lnum_t  *f2e_idx = connect->f2e->idx + f_id;
  const cs_lnum_t  *f2e_lst = connect->f2e->ids + f2e_idx[0];

  fm->n_vf = fm->n_ef = f2e_idx[1] - f2e_idx[0];
  short int nv = 0;
  for (int i = 0; i < fm->n_vf; i++)
    fm->v_ids[i] = -1;

  for (short int e = 0; e < fm->n_ef; e++) {

    const cs_lnum_t  e_id = f2e_lst[e];
    const cs_nvec3_t  e_nvect = cs_quant_set_edge_nvec(e_id, quant);

    fm->e_ids[e] = e_id;
    fm->edge[e].meas = e_nvect.meas;
    for (int k = 0; k < 3; k++)
      fm->edge[e].unitv[k] = e_nvect.unitv[k];
    /* Still to handle the edge barycenter */

    const cs_lnum_t  *e2v_ids = connect->e2v->ids + 2*e_id;
    short int  v1 = -1, v2 = -1;
    for (int v = 0; v < fm->n_vf && fm->v_ids[v] != -1; v++) {
      if (fm->v_ids[v] == e2v_ids[0])
        v1 = v;
      else if (fm->v_ids[v] == e2v_ids[1])
        v2 = v;
    }

    /* Add vertices if not already identified */
    if (v1 == -1) /* Not found -> Add v1 */
      fm->v_ids[nv] = e2v_ids[0], v1 = nv++;
    if (v2 == -1) /* Not found -> Add v2 */
      fm->v_ids[nv] = e2v_ids[1], v2 = nv++;

    /* Update e2v_ids */
    const int _eshft = 2*e;
    fm->e2v_ids[_eshft]   = v1;
    fm->e2v_ids[_eshft+1] = v2;

  } /* Loop on face edges */

  assert(nv == fm->n_vf); /* Sanity check */

  /* Update vertex coordinates */
  int  shift = 0;
  for (short int v = 0; v < fm->n_vf; v++) {
    const cs_real_t *xv = quant->vtx_coord + 3*fm->v_ids[v];
    for (int k = 0; k < 3; k++)
      fm->xv[shift++] = xv[k];
  }

  /* Update the edge center. Define wvf and tef */
  for (int i = 0; i < fm->n_vf; i++)
    fm->wvf[i] = 0;

  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];
    const cs_real_t  *xv1 = fm->xv + 3*v1;
    const cs_real_t  *xv2 = fm->xv + 3*v2;

    /* Update the edge center */
    for (int k = 0; k < 3; k++)
      fm->edge[e].center[k] = 0.5 * (xv1[k] + xv2[k]);

    /* tef = ||(xe -xf) x e||/2 = s(v1,e,f) + s(v2, e, f) */
    const double  tef = cs_compute_area_from_quant(fm->edge[e], pfq.center);

    fm->wvf[v1] += tef;
    fm->wvf[v2] += tef;
    fm->tef[e] = tef;

  } /* Loop on face edges */

  const double  invf = 0.5/pfq.meas;
  for (short int v = 0; v < fm->n_vf; v++) fm->wvf[v] *= invf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_face_mesh_t structure for a given cell from a
 *         cs_cell_mesh_t structure
 *
 * \param[in]       cm     pointer to the reference cs_cell_mesh_t structure
 * \param[in]       f      face id in the cs_cell_mesh_t structure
 * \param[in, out]  fm     pointer to a cs_face_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_build_from_cell_mesh(const cs_cell_mesh_t    *cm,
                                  short int                f,
                                  cs_face_mesh_t          *fm)
{
  if (fm == NULL || cm == NULL)
    return;

  /* Sanity checks */
  assert(f > -1 && f < cm->n_fc);
  assert(cs_flag_test(cm->flag,
                      CS_FLAG_COMP_PV  | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ |
                      CS_FLAG_COMP_PEQ | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  fm->c_id = cm->c_id;
  for (int k = 0; k < 3; k++) fm->xc[k] = cm->xc[k];

  /* Face-related quantities */
  fm->f_id = f;
  fm->f_sgn = cm->f_sgn[f];

  const cs_quant_t  pfq = cm->face[f];
  fm->face.meas = pfq.meas;
  for (int k = 0; k < 3; k++) {
    fm->face.center[k] = pfq.center[k];
    fm->face.unitv[k] = pfq.unitv[k];
  }

  const cs_nvec3_t  deq = cm->dedge[f];
  fm->dedge.meas = deq.meas;
  for (int k = 0; k < 3; k++)
    fm->dedge.unitv[k] = deq.unitv[k];

  const short int  *f2e_idx = cm->f2e_idx + f;
  const short int  *f2e_ids = cm->f2e_ids + f2e_idx[0];
  const double *_tef = cm->tef + f2e_idx[0];

  fm->n_vf = fm->n_ef = f2e_idx[1] - f2e_idx[0];
  short int nv = 0;
  for (int i = 0; i < fm->n_vf; i++)  fm->v_ids[i] = -1;

  for (short int ef = 0; ef < fm->n_ef; ef++) {

    const short int  ec = f2e_ids[ef];

    fm->e_ids[ef] = ec;
    fm->tef[ef] = _tef[ef];

    const cs_quant_t  peq = cm->edge[ec];
    fm->edge[ef].meas = peq.meas;
    for (int k = 0; k < 3; k++) {
      fm->edge[ef].center[k] = peq.center[k];
      fm->edge[ef].unitv[k] = peq.unitv[k];
    }

    const int  eshft = 2*ec;
    short int  v1c_id = cm->e2v_ids[eshft];
    short int  v2c_id = cm->e2v_ids[eshft+1];

    /* Compact vertex numbering to this face */
    short int  v1 = -1, v2 = -1;
    for (int v = 0; v < fm->n_vf && fm->v_ids[v] != -1; v++) {
      if (fm->v_ids[v] == v1c_id)
        v1 = v;
      else if (fm->v_ids[v] == v2c_id)
        v2 = v;
    }

    /* Add vertices if not already identified */
    if (v1 == -1) /* Not found -> Add v1 */
      fm->v_ids[nv] = v1c_id, v1 = nv++;
    if (v2 == -1) /* Not found -> Add v2 */
      fm->v_ids[nv] = v2c_id, v2 = nv++;

    /* Update e2v_ids */
    const int _eshft = 2*ef;
    fm->e2v_ids[_eshft]   = v1;
    fm->e2v_ids[_eshft+1] = v2;

  } /* Loop on face edges */

  assert(nv == fm->n_vf); /* Sanity check */

  /* Update vertex coordinates */
  int  shift = 0;
  for (short int v = 0; v < fm->n_vf; v++) {
    const cs_real_t *xv = cm->xv + 3*fm->v_ids[v];
    for (int k = 0; k < 3; k++)
      fm->xv[shift++] = xv[k];
  }

  /* Define wvf and tef */
  for (int i = 0; i < fm->n_vf; i++)
    fm->wvf[i] = 0;

  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];

    fm->wvf[v1] += _tef[e];
    fm->wvf[v2] += _tef[e];

  }

  const double  invf = 0.5/pfq.meas;
  for (short int v = 0; v < fm->n_vf; v++) fm->wvf[v] *= invf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_face_mesh_light_t structure
 *
 * \param[in]  n_max_vbyf    max. number of vertices for a face
 * \param[in]  n_max_vbyc    max. number of vertices for a cell
 *
 * \return a pointer to a new allocated cs_face_mesh_light_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_light_t *
cs_face_mesh_light_create(short int   n_max_vbyf,
                          short int   n_max_vbyc)
{
  cs_face_mesh_light_t  *fm = NULL;

  BFT_MALLOC(fm, 1, cs_face_mesh_light_t);

  fm->n_max_vbyf = n_max_vbyf;
  fm->c_id = -1;
  fm->f = -1;

  /* Vertex-related quantities */
  fm->n_vf = 0;
  BFT_MALLOC(fm->v_ids, n_max_vbyc, short int);
  BFT_MALLOC(fm->wvf, n_max_vbyc, double);

  /* Edge-related quantities */
  fm->n_ef = 0;
  BFT_MALLOC(fm->e_ids, fm->n_max_vbyf, short int);
  BFT_MALLOC(fm->tef, fm->n_max_vbyf, double);

  return fm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_face_mesh_light_t structure corresponding to
 *         mesh id
 *
 * \param[in]   mesh_id   id in the cs_face_mesh_light_t array
 *
 * \return a pointer to a cs_face_mesh_light_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_light_t *
cs_cdo_local_get_face_mesh_light(int    mesh_id)
{
  if (mesh_id < 0 || mesh_id >= cs_glob_n_threads)
    return NULL;

  return cs_cdo_local_face_meshes_light[mesh_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_face_mesh_light_t structure
 *
 * \param[in, out]  p_fm   pointer of pointer to a cs_face_mesh_light_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_light_free(cs_face_mesh_light_t     **p_fm)
{
  cs_face_mesh_light_t  *fm = *p_fm;

  if (fm == NULL)
    return;

  BFT_FREE(fm->v_ids);
  BFT_FREE(fm->wvf);
  BFT_FREE(fm->e_ids);
  BFT_FREE(fm->tef);

  BFT_FREE(fm);
  *p_fm = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_face_mesh_light_t structure starting from a
 *         cs_cell_mesh_t structure.
 *
 * \param[in]       cm     pointer to the reference cs_cell_mesh_t structure
 * \param[in]       f      face id in the cs_cell_mesh_t structure
 * \param[in, out]  fm     pointer to a cs_face_mesh_light_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_light_build(const cs_cell_mesh_t    *cm,
                         short int                f,
                         cs_face_mesh_light_t    *fm)
{
  if (fm == NULL || cm == NULL)
    return;

  /* Sanity checks */
  assert(f > -1 && f < cm->n_fc);
  assert(cs_flag_test(cm->flag,
                      CS_FLAG_COMP_PV | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV));

  fm->c_id = cm->c_id;
  fm->f = f;

  const short int  *f2e_idx = cm->f2e_idx + f;
  const short int  *f2e_ids = cm->f2e_ids + f2e_idx[0];

  /* Initialization */
  fm->n_vf = fm->n_ef = f2e_idx[1] - f2e_idx[0];
  for (int i = 0; i < cm->n_vc; i++) {
    fm->v_ids[i] = -1;
    fm->wvf[i] = 0;
  }

  /* Define wvf from the knowledge of tef */
  const double *_tef = cm->tef + f2e_idx[0];

  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  e_cellwise = f2e_ids[e];
    const int  eshft = 2*e_cellwise;
    const short int  v1_cellwise = cm->e2v_ids[eshft];
    const short int  v2_cellwise = cm->e2v_ids[eshft+1];

    fm->e_ids[e] = e_cellwise;
    fm->tef[e] = _tef[e];
    fm->v_ids[v1_cellwise] = 1;
    fm->v_ids[v2_cellwise] = 1;

    /* Build wvf */
    fm->wvf[v1_cellwise] += _tef[e];
    fm->wvf[v2_cellwise] += _tef[e];

  } /* Loop on face edges */

  /* Compact vertex numbering to this face */
  short int nv = 0; /* current vertex id in the face numbering */
  for (short int v = 0; v < cm->n_vc; v++) {
    if (fm->v_ids[v] > 0) {
      fm->v_ids[nv] = v;
      fm->wvf[nv] = fm->wvf[v];
      nv++;
    }
  }

  assert(nv == fm->n_vf); /* Sanity check */
  const double  invf = 0.5/cm->face[f].meas;
  for (short int v = 0; v < fm->n_vf; v++) fm->wvf[v] *= invf;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
