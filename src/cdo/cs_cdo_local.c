/* ===========================================================================
 * Routines to handle low-level routines related to CDO local quantities:
 * - local matrices (stored in dense format),
 * - local mesh structure related to a cell or to a couple cell/face
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
#include <float.h>
#include <limits.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_log.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_local.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

cs_cell_mesh_t   **cs_cdo_local_cell_meshes = NULL;
cs_face_mesh_t   **cs_cdo_local_face_meshes = NULL;

static int  cs_cdo_local_n_structures = 0;

/* Store predefined flags */
static const cs_flag_t  cs_cdo_local_flag_v =
  CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_EV;
static const cs_flag_t  cs_cdo_local_flag_e =
  CS_CDO_LOCAL_PE | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_EV |
  CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EF  | CS_CDO_LOCAL_EFQ;
static const cs_flag_t  cs_cdo_local_flag_peq =
  CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_FEQ;
static const cs_flag_t  cs_cdo_local_flag_f =
  CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_FE |
  CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EF | CS_CDO_LOCAL_EFQ | CS_CDO_LOCAL_HFQ;
static const cs_flag_t  cs_cdo_local_flag_pfq =
  CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FEQ;
static const cs_flag_t  cs_cdo_local_flag_deq =
  CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_DEQ;
static const cs_flag_t  cs_cdo_local_flag_fe =
  CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EF | CS_CDO_LOCAL_EFQ;
static const cs_flag_t  cs_cdo_local_flag_ef =
  CS_CDO_LOCAL_EF | CS_CDO_LOCAL_EFQ;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cell_sys_t structure
 *
 * \param[in]   n_max_ent    max number of entries
 *
 * \return a pointer to a new allocated cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_sys_t *
cs_cell_sys_create(int    n_max_ent)
{
  cs_cell_sys_t  *csys = NULL;

  BFT_MALLOC(csys, 1, cs_cell_sys_t);

  csys->n_dofs = 0;
  csys->mat = NULL;
  csys->rhs = NULL;
  csys->source = NULL;
  csys->val_n = NULL;

  if (n_max_ent > 0) {
    csys->mat = cs_locmat_create(n_max_ent);
    BFT_MALLOC(csys->rhs, n_max_ent, double);
    BFT_MALLOC(csys->source, n_max_ent, double);
    BFT_MALLOC(csys->val_n, n_max_ent, double);
    for (int i = 0; i < n_max_ent; i++)
      csys->rhs[i] = csys->source[i] = csys->val_n[i] = 0.0;
  }

  return csys;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cell_sys_t structure
 *
 * \param[in, out]  p_csys   pointer of pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_sys_free(cs_cell_sys_t     **p_csys)
{
  cs_cell_sys_t  *csys = *p_csys;

  if (csys == NULL)
    return;

  csys->mat = cs_locmat_free(csys->mat);
  BFT_FREE(csys->rhs);
  BFT_FREE(csys->source);
  BFT_FREE(csys->val_n);

  BFT_FREE(csys);
  *p_csys= NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a local system for debugging purpose
 *
 * \param[in]       msg     associated message to print
 * \param[in]       c_id    id related to the cell
 * \param[in]       csys    pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_sys_dump(const char             msg[],
                 const cs_lnum_t        c_id,
                 const cs_cell_sys_t   *csys)
{
#pragma omp critical
  {
    cs_log_printf(CS_LOG_DEFAULT, "%s", msg);

    cs_locmat_dump(c_id, csys->mat);
    cs_log_printf(CS_LOG_DEFAULT,
                  "\n>> %-10s | %-10s | %-10s\n",
                  "RHS", "TS", "VAL_PREV");
    for (int i = 0; i < csys->n_dofs; i++)
      cs_log_printf(CS_LOG_DEFAULT,
                    ">> % .3e | % .3e | % .3e\n",
                    csys->rhs[i], csys->source[i], csys->val_n[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cell_bc_t structure
 *
 * \param[in]   n_max_dofbyc    max. number of entries in a cell
 * \param[in]   n_max_fbyc      max. number of faces in a cell
 *
 * \return a pointer to a new allocated cs_cell_bc_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_bc_t *
cs_cell_bc_create(int    n_max_dofbyc,
                  int    n_max_fbyc)
{
  cs_cell_bc_t  *cbc = NULL;

  BFT_MALLOC(cbc, 1, cs_cell_bc_t);

  cbc->n_bc_faces = 0;
  cbc->bf_ids = NULL;
  cbc->face_flag = NULL;
  if (n_max_fbyc > 0) {
    BFT_MALLOC(cbc->bf_ids, n_max_fbyc, short int);
    BFT_MALLOC(cbc->face_flag, n_max_fbyc, cs_flag_t);
  }

  cbc->n_dofs = 0;
  cbc->n_dirichlet = 0;
  cbc->n_nhmg_neuman = 0;
  cbc->n_robin = 0;

  if (n_max_dofbyc > 0) { // Max. number of DoFs in a cell

    BFT_MALLOC(cbc->dof_flag, n_max_dofbyc, cs_flag_t);
    BFT_MALLOC(cbc->dir_values, n_max_dofbyc, double);
    BFT_MALLOC(cbc->neu_values, n_max_dofbyc, double);
    BFT_MALLOC(cbc->rob_values, 2*n_max_dofbyc, double);

    for (short int i = 0; i < n_max_dofbyc; i++) {
      cbc->dof_flag[i] = 0;
      cbc->dir_values[i] = 0;
      cbc->neu_values[i] = 0;
      cbc->rob_values[2*i] = cbc->rob_values[2*i+1] = 0;
    }

  }

  return cbc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cell_bc_t structure
 *
 * \param[in, out]  p_cbc   pointer of pointer to a cs_cell_bc_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_bc_free(cs_cell_bc_t     **p_cbc)
{
  cs_cell_bc_t  *cbc = *p_cbc;

  if (cbc == NULL)
    return;

  BFT_FREE(cbc->bf_ids);
  BFT_FREE(cbc->face_flag);
  BFT_FREE(cbc->dof_flag);
  BFT_FREE(cbc->dir_values);
  BFT_FREE(cbc->neu_values);
  BFT_FREE(cbc->rob_values);

  BFT_FREE(cbc);
  *p_cbc = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a cs_cell_builder_t structure according to
 *         to the type of discretization which is requested.
 *
 * \param[in]  scheme    type of discretization
 * \param[in]  connect   pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to the new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_builder_t *
cs_cell_builder_create(cs_space_scheme_t         scheme,
                       const cs_cdo_connect_t   *connect)
{
  int  size;

  cs_cell_builder_t  *cb = NULL;

  /* Common part to all discretization */
  BFT_MALLOC(cb, 1, cs_cell_builder_t);

  cb->eig_ratio = -DBL_MAX;
  cb->eig_max = -DBL_MAX;

  cb->pty_mat[0][0] = cb->pty_mat[1][1] = cb->pty_mat[2][2] = 1;
  cb->pty_mat[0][1] = cb->pty_mat[1][0] = cb->pty_mat[2][0] = 0;
  cb->pty_mat[0][2] = cb->pty_mat[1][2] = cb->pty_mat[2][1] = 0;
  cb->pty_val = 1;

  cb->ids = NULL;
  cb->values = NULL;
  cb->vectors = NULL;

  cb->hdg = cb->loc = cb->aux = NULL;

  /* Max number of entities in a cell */
  const int  n_vc = connect->n_max_vbyc;
  const int  n_ec = connect->n_max_ebyc;
  const int  n_fc = connect->n_max_fbyc;

  /* Temporary buffers used during the cellwise construction of operators.
     The size of the allocation corresponds to the max. possible size
     according to any numerical options. This means that one could optimize
     the following allocation if one takes into account the numerical
     settings */

  switch (scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    {
      BFT_MALLOC(cb->ids, n_vc, short int);
      for (int i = 0; i < n_vc; i++) cb->ids[i] = 0;

      size = n_ec*(n_ec+1);
      size = CS_MAX(4*n_ec + 3*n_vc, size);
      BFT_MALLOC(cb->values, size, double);
      for (int i = 0; i < size; i++) cb->values[i] = 0;

      size = 2*n_ec;
      BFT_MALLOC(cb->vectors, size, cs_real_3_t);
      for (int i = 0; i < size; i++)
        cb->vectors[i][0] = cb->vectors[i][1] = cb->vectors[i][2] = 0;

      /* Local square dense matrices used during the construction of
         operators */
      cb->hdg = cs_locmat_create(n_ec);
      cb->loc = cs_locmat_create(n_vc);
      cb->aux = cs_locmat_create(n_vc);
    }
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    {
      BFT_MALLOC(cb->ids, n_vc + 1, short int);
      for (int i = 0; i < n_vc + 1; i++) cb->ids[i] = 0;

      size = 2*n_vc + 3*n_ec + n_fc;
      BFT_MALLOC(cb->values, size, double);
      for (int i = 0; i < size; i++) cb->values[i] = 0;

      size = 2*n_ec + n_vc;
      BFT_MALLOC(cb->vectors, size, cs_real_3_t);
      for (int i = 0; i < size; i++)
        cb->vectors[i][0] = cb->vectors[i][1] = cb->vectors[i][2] = 0;

      /* Local square dense matrices used during the construction of
         operators */
      cb->hdg = cs_locmat_create(n_vc + 1);
      cb->loc = cs_locmat_create(n_vc + 1);
      cb->aux = cs_locmat_create(n_vc + 1);

    }
    break;

  case CS_SPACE_SCHEME_CDOFB:

    size = 4*(n_fc + 1);
    BFT_MALLOC(cb->values, size, double);
    for (int i = 0; i < size; i++) cb->values[i] = 0;

    size = n_fc;
    BFT_MALLOC(cb->vectors, size, cs_real_3_t);
    for (int i = 0; i < size; i++)
      cb->vectors[i][0] = cb->vectors[i][1] = cb->vectors[i][2] = 0;

    /* Local square dense matrices used during the construction of
       operators */
    cb->hdg = cs_locmat_create(n_fc + 1);
    cb->loc = cs_locmat_create(n_fc + 1);
    cb->aux = cs_locmat_create(n_fc + 1);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _("Invalid space scheme."));

  } // End of switch on space scheme

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cell_builder_t structure
 *
 * \param[in, out]  p_cb   pointer of pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_builder_free(cs_cell_builder_t     **p_cb)
{
  cs_cell_builder_t  *cb = *p_cb;

  if (cb == NULL)
    return;

  BFT_FREE(cb->ids);
  BFT_FREE(cb->values);
  BFT_FREE(cb->vectors);

  cb->hdg = cs_locmat_free(cb->hdg);
  cb->loc = cs_locmat_free(cb->loc);
  cb->aux = cs_locmat_free(cb->aux);

  BFT_FREE(cb);
  *p_cb = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate global structures related to cs_cell_mesh_t and
 *         cs_face_mesh_t structures
 *
 * \param[in]   connect   pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_local_initialize(const cs_cdo_connect_t     *connect)
{
  /* Sanity check */
  assert(cs_glob_n_threads > 0);

  int  size = cs_glob_n_threads;
  cs_cdo_local_n_structures = size;
  BFT_MALLOC(cs_cdo_local_cell_meshes, size, cs_cell_mesh_t *);
  BFT_MALLOC(cs_cdo_local_face_meshes, size, cs_face_mesh_t *);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cdo_local_cell_meshes[t_id] = cs_cell_mesh_create(connect);
    cs_cdo_local_face_meshes[t_id] = cs_face_mesh_create(connect->n_max_vbyf);
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cdo_local_cell_meshes[0] = cs_cell_mesh_create(connect);
  cs_cdo_local_face_meshes[0] = cs_face_mesh_create(connect->n_max_vbyf);
#endif /* openMP */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free global structures related to cs_cell_mesh_t and cs_face_mesh_t
 *         structures
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
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_mesh_free(&(cs_cdo_local_cell_meshes[0]));
  cs_face_mesh_free(&(cs_cdo_local_face_meshes[0]));
#endif /* openMP */

  BFT_FREE(cs_cdo_local_cell_meshes);
  BFT_FREE(cs_cdo_local_face_meshes);
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

  return cs_cdo_local_cell_meshes[mesh_id];
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
  BFT_MALLOC(cm->kbuf, CS_MAX(connect->v_max_cell_range,
                              connect->e_max_cell_range) + 1, short int);
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
  BFT_MALLOC(cm->hfc, cm->n_max_fbyc, double);
  BFT_MALLOC(cm->face, cm->n_max_fbyc, cs_quant_t);
  BFT_MALLOC(cm->dedge, cm->n_max_fbyc, cs_nvec3_t);

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
 * \brief  Initialize to invalid values a cs_cell_mesh_t structure
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_reset(cs_cell_mesh_t   *cm)
{
  cm->n_vc = SHRT_MAX;
  cm->n_ec = SHRT_MAX;
  cm->n_fc = SHRT_MAX;

  /* Cell information */
  cm->c_id = SHRT_MIN;
  cm->xc[0] = cm->xc[1] = cm->xc[2] = -DBL_MAX;
  cm->vol_c = -DBL_MAX;

  /* Vertex information */
  for (short int v = 0; v < cm->n_max_vbyc; v++) {
    cm->v_ids[v] = SHRT_MIN;
    cm->wvc[v] = -DBL_MAX;
    cm->xv[3*v] = cm->xv[3*v+1] = cm->xv[3*v+2] = -DBL_MAX;
  }

  /* Edge information */
  for (short int e = 0; e < cm->n_max_ebyc; e++) {
    cm->e_ids[e] = SHRT_MIN;
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
    cm->f_ids[f] = SHRT_MIN;
    cm->f_sgn[f] = 0;
    cm->hfc[f] = -DBL_MAX;
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
    cm->f2e_idx[f] = SHRT_MIN;

  for (int i = 0; i < 2*cm->n_max_ebyc; i++) {
    cm->e2v_ids[i] = cm->e2f_ids[i] = cm->f2e_ids[i] = SHRT_MIN;
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
cs_cell_mesh_dump(cs_cell_mesh_t     *cm)
{
  if (cm == NULL) {
    cs_log_printf(CS_LOG_DEFAULT, "\n>> Dump cs_cell_mesh_t %p\n", (void *)cm);
    return;
  }

  cs_log_printf(CS_LOG_DEFAULT, "\n>> Dump cs_cell_mesh_t %p; flag: %d\n"
                " c_id:%d; vol: %9.6e; xc (% .5e % .5e % .5e)\n",
                (void *)cm, cm->flag, cm->c_id, cm->vol_c,
                cm->xc[0], cm->xc[1], cm->xc[2]);

  /* Information related to primal vertices */
  if (cm->flag & cs_cdo_local_flag_v) {

    cs_log_printf(CS_LOG_DEFAULT, "%-3s %-9s %-38s %-9s\n",
                  "v", "id", "coord", "wvc");
    for (short int v = 0; v < cm->n_vc; v++)
      cs_log_printf(CS_LOG_DEFAULT, "%2d |%8d |% .5e % .5e % .5e| %.5e\n",
                    v, cm->v_ids[v], cm->xv[3*v], cm->xv[3*v+1], cm->xv[3*v+2],
                    cm->wvc[v]);

  } /* Vertex quantities */

  /* Information related to primal edges */
  if (cm->flag & cs_cdo_local_flag_e) {

    cs_log_printf(CS_LOG_DEFAULT, "%-3s %-9s %-9s %-38s %-38s %-11s %-38s\n",
                  "e", "id", "length", "unit", "coords", "df.meas", "df.unit");
    for (short int e = 0; e < cm->n_ec; e++) {
      cs_quant_t  peq = cm->edge[e];
      cs_nvec3_t  dfq = cm->dface[e];
      cs_log_printf(CS_LOG_DEFAULT, "%2d |%8d |%.3e|% .5e % .5e % .5e|"
                    "% .5e % .5e % .5e|%.5e|% .5e % .5e % .5e\n",
                    e, cm->e_ids[e], peq.meas, peq.unitv[0], peq.unitv[1],
                    peq.unitv[2], peq.center[0], peq.center[1], peq.center[2],
                    dfq.meas, dfq.unitv[0], dfq.unitv[1], dfq.unitv[2]);
    }

  } /* Edge quantities */

  /* Information related to primal faces */
  if (cm->flag & cs_cdo_local_flag_f) {

    cs_log_printf(CS_LOG_DEFAULT, "%-3s %-9s %-9s %-4s %-38s %-38s %-11s"
                  "%-11s %-38s\n",
                  "f", "id", "surf", "sgn", "unit", "coords", "hfc", "dlen",
                  "dunitv");
    for (short int f = 0; f < cm->n_fc; f++) {
      cs_quant_t  fq = cm->face[f];
      cs_nvec3_t  eq = cm->dedge[f];
      cs_log_printf(CS_LOG_DEFAULT, "%2d |%8d |%.3e| %2d|% .5e % .5e % .5e|"
                    "% .5e % .5e % .5e|%.5e|%.5e|% .5e % .5e % .5e\n",
                    f, cm->f_ids[f], fq.meas, cm->f_sgn[f], fq.unitv[0],
                    fq.unitv[1], fq.unitv[2], fq.center[0], fq.center[1],
                    fq.center[2], cm->hfc[f], eq.meas, eq.unitv[0], eq.unitv[1],
                    eq.unitv[2]);
    }

  } /* Face quantities */

  if (cm->flag & CS_CDO_LOCAL_EV) {

    cs_log_printf(CS_LOG_DEFAULT, "%-2s (v1, v2) sgn\n", "e");
    for (short int e = 0; e < cm->n_ec; e++)
      cs_log_printf(CS_LOG_DEFAULT, "%2d (%2d, %2d) %2d\n",
                    e, cm->e2v_ids[2*e], cm->e2v_ids[2*e+1], cm->e2v_sgn[e]);

  }

  if (cm->flag & cs_cdo_local_flag_fe) {

    cs_log_printf(CS_LOG_DEFAULT, " n_ef | f: pef\n");
    for (short int f = 0; f < cm->n_fc; f++) {
      cs_log_printf(CS_LOG_DEFAULT, " %4d |",
                    cm->f2e_idx[f+1] - cm->f2e_idx[f]);
      for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++)
        cs_log_printf(CS_LOG_DEFAULT, " %2d:%.5e|", cm->f2e_ids[i], cm->tef[i]);
      cs_log_printf(CS_LOG_DEFAULT, "\n");
    }

  }

  if (cm->flag & cs_cdo_local_flag_ef) {

    cs_log_printf(CS_LOG_DEFAULT, "%-4s | f0 | %-53s | f1 | %-53s\n",
                  "e", "sef0c: meas, unitv", "sef1c: meas, unitv");
    for (short int e = 0; e < cm->n_ec; e++)
      cs_log_printf(CS_LOG_DEFAULT,
                    " %3d | %2d | % .5e (% .5e % .5e % .5e) |"
                    " %2d | % .5e (% .5e % .5e % .5e)\n",
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

  BFT_FREE(cm->kbuf);

  BFT_FREE(cm->v_ids);
  BFT_FREE(cm->wvc);
  BFT_FREE(cm->xv);

  BFT_FREE(cm->e_ids);
  BFT_FREE(cm->edge);
  BFT_FREE(cm->dface);

  BFT_FREE(cm->f_ids);
  BFT_FREE(cm->f_sgn);
  BFT_FREE(cm->hfc);
  BFT_FREE(cm->face);
  BFT_FREE(cm->dedge);

  BFT_FREE(cm->e2v_ids);
  BFT_FREE(cm->e2v_sgn);

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
 * \param[in]       c_id      cell id
 * \param[in]       flag      indicate which members are really defined
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  cm        pointer to a cs_cell_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_build(cs_lnum_t                    c_id,
                   cs_flag_t                    flag,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_cell_mesh_t              *cm)
{
  if (cm == NULL)
    return;

  cm->flag = flag;

  /* Store information related to cell */
  cm->c_id = c_id;
  cm->vol_c = quant->cell_vol[c_id];
  for (int k = 0; k < 3; k++)
    cm->xc[k] = quant->cell_centers[3*c_id+k];

  const cs_lnum_t  *c2v_idx = connect->c2v->idx + c_id;
  const cs_lnum_t  *c2e_idx = connect->c2e->idx + c_id;
  const cs_lnum_t  *c2f_idx = connect->c2f->idx + c_id;

  cm->n_vc = c2v_idx[1] - c2v_idx[0];
  cm->n_ec = c2e_idx[1] - c2e_idx[0];
  cm->n_fc = c2f_idx[1] - c2f_idx[0];

  if (flag == 0)
    return;

  /* Information related to primal vertices */
  if (flag & cs_cdo_local_flag_v) {

    const cs_lnum_t  *c2v_ids = connect->c2v->ids + c2v_idx[0];

    for (short int v = 0; v < cm->n_vc; v++) {
      const cs_lnum_t  v_id = c2v_ids[v];
      cm->v_ids[v] = v_id;
      for (int k = 0; k < 3; k++)
        cm->xv[3*v+k] = quant->vtx_coord[3*v_id+k];
    }

    /* Primal vertices quantities */
    if (flag & CS_CDO_LOCAL_PVQ) {

      const double  *wvc = quant->dcell_vol + c2v_idx[0];
      const double  invvol = 1/cm->vol_c;
      for (short int v = 0; v < cm->n_vc; v++)
        cm->wvc[v] = invvol * wvc[v];

    }

  } // vertices

  /* Information related to primal edges */
  if (flag & cs_cdo_local_flag_e) {

    const cs_lnum_t  *c2e_ids = connect->c2e->ids + c2e_idx[0];

    for (short int e = 0; e < cm->n_ec; e++)
      cm->e_ids[e] = c2e_ids[e];

    if (flag & cs_cdo_local_flag_peq) {

      /* Primal edge quantities */
      for (short int e = 0; e < cm->n_ec; e++) {

        const cs_quant_t  qe = quant->edge[cm->e_ids[e]];
        cm->edge[e].meas = qe.meas;
        for (int k = 0; k < 3; k++) {
          cm->edge[e].center[k] = qe.center[k];
          cm->edge[e].unitv[k] = qe.unitv[k];
        }

      }

    } // Primal edge quantities

    /* Dual face quantities related to each edge */
    if (flag & CS_CDO_LOCAL_DFQ) {

      const cs_dface_t  *qdf = quant->dface + c2e_idx[0];

      for (short int e = 0; e < cm->n_ec; e++) {
        cm->dface[e].meas = qdf[e].sface[0].meas + qdf[e].sface[1].meas;
        const double  inv_meas = 1/cm->dface[e].meas;
        for (int k = 0; k < 3; k++)
          cm->dface[e].unitv[k] = inv_meas*qdf[e].vect[k];
      }

    } // Dual face quantities

  } // Edges

  if (flag & CS_CDO_LOCAL_EV) {

    /* Store in compact way: mesh --> cell mesh ids for vertices */
    cs_lnum_t  shift = cm->v_ids[0];
    for (short int v = 1; v < cm->n_vc; v++)
      if (cm->v_ids[v] < shift)
        shift = cm->v_ids[v];
    for (short int v = 0; v < cm->n_vc; v++)
      cm->kbuf[cm->v_ids[v]-shift] = v;

    for (short int e = 0; e < cm->n_ec; e++) {

      const cs_lnum_t  e_id = cm->e_ids[e];

      /* Store only the sign related to the first vertex since the sign
         related to the second one is minus the first one */
      cm->e2v_sgn[e] = connect->e2v->sgn[2*e_id];
      cm->e2v_ids[2*e]   = cm->kbuf[connect->e2v->col_id[2*e_id] - shift];
      cm->e2v_ids[2*e+1] = cm->kbuf[connect->e2v->col_id[2*e_id+1] - shift];

    } // Loop on cell edges

  } // edge-vertices information

  /* Information related to primal faces */
  if (flag & cs_cdo_local_flag_f) {

    const cs_lnum_t  *c2f_lst = connect->c2f->col_id + c2f_idx[0];
    const short int  *c2f_sgn = connect->c2f->sgn + c2f_idx[0];

    for (short int f = 0; f < cm->n_fc; f++) {
      cm->f_ids[f] = c2f_lst[f];
      cm->f_sgn[f] = c2f_sgn[f];
    } // Loop on cell faces

    /* Face related quantities */
    if (flag & cs_cdo_local_flag_pfq) {

      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_quant_t  pfq = quant->face[cm->f_ids[f]];
        cm->face[f].meas = pfq.meas;
        for (int k = 0; k < 3; k++) {
          cm->face[f].center[k] = pfq.center[k];
          cm->face[f].unitv[k] = pfq.unitv[k];
        }

      }

    } /* Primal face quantities */

    if (flag & cs_cdo_local_flag_deq) {

      const cs_nvec3_t  *deq = quant->dedge + c2f_idx[0];
      for (short int f = 0; f < cm->n_fc; f++) {
        cm->dedge[f].meas = deq[f].meas;
        for (int k = 0; k < 3; k++)
          cm->dedge[f].unitv[k] = deq[f].unitv[k];
      }

    } /* Dual edge quantities */

    if (flag & CS_CDO_LOCAL_HFQ) {

      /* Compute the height of the pyramid of base f whose apex is
         the cell center */
      for (short int f = 0; f < cm->n_fc; f++) {
        cm->hfc[f] = cs_math_3_dot_product(cm->face[f].unitv,
                                           cm->dedge[f].unitv);
        cm->hfc[f] *= cm->dedge[f].meas;
        assert(cm->hfc[f] > 0);
      }

    } /* Quantities related to the pyramid of base f */

  } /* Face information */

  if (flag & cs_cdo_local_flag_fe) {

    /* Store in compact way: mesh --> cell mesh ids for edges */
    cs_lnum_t  shift = cm->e_ids[0];
    for (short int e = 1; e < cm->n_ec; e++)
      if (cm->e_ids[e] < shift)
        shift = cm->e_ids[e];
    for (short int e = 0; e < cm->n_ec; e++)
      cm->kbuf[cm->e_ids[e]-shift] = e;

    const cs_lnum_t  *f2e_idx = connect->f2e->idx;
    const cs_lnum_t  *f2e_ids = connect->f2e->col_id;

    cm->f2e_idx[0] = 0;
    int shift_idx = 0;
    for (short int f = 0; f < cm->n_fc; f++) {

      const cs_lnum_t  f_id = cm->f_ids[f];

      /* Build index */
      const cs_lnum_t  f2e_start = f2e_idx[f_id];
      const cs_lnum_t  f2e_end = f2e_idx[f_id+1];

      cm->f2e_idx[f+1] = cm->f2e_idx[f] + f2e_start - f2e_end;

      for (cs_lnum_t i = f2e_start; i < f2e_end; i++)
        cm->f2e_ids[shift_idx++] = cm->kbuf[f2e_ids[i] - shift];
      cm->f2e_idx[f+1] = shift_idx;

    } // Loop on cell faces

    /* Sanity check */
    assert(cm->f2e_idx[cm->n_fc] == 2*cm->n_ec);

    if (flag & CS_CDO_LOCAL_FEQ) {

      for (short int f = 0; f < cm->n_fc; f++) {
        for (int ie = cm->f2e_idx[f]; ie < cm->f2e_idx[f+1]; ie++)
          cm->tef[ie] = cs_compute_area_from_quant(cm->edge[cm->f2e_ids[ie]],
                                                   cm->face[f].center);
      }

    } // face --> edges quantities

  } // face --> edges connectivity

  if (flag & cs_cdo_local_flag_ef) {

    /* Build the e2f connectivity */
    for (short int i = 0; i < 2*cm->n_ec; i++) cm->e2f_ids[i] = -1;

    for (short int f = 0; f < cm->n_fc; f++) {
      for (short int jf = cm->f2e_idx[f]; jf < cm->f2e_idx[f+1]; jf++) {

        const short int  e = cm->f2e_ids[jf];
        if (cm->e2f_ids[2*e] == -1)
          cm->e2f_ids[2*e] = f;
        else {
          assert(cm->e2f_ids[2*e+1] == -1);
          cm->e2f_ids[2*e+1] = f;
        }

      } /* Loop on face edges */
    } // Loop on cell faces

    if (flag & CS_CDO_LOCAL_EFQ) { /* Build cm->sefc */

      cs_dface_t  *qdf = quant->dface + connect->c2e->idx[c_id];

      for (short int e = 0; e < cm->n_ec; e++) {

        short int  ee = 2*e;
        /* Sanity checks */
        assert(cm->f_ids[cm->e2f_ids[ee]] == qdf[e].parent_id[0]);
        assert(cm->f_ids[cm->e2f_ids[ee+1]] == qdf[e].parent_id[1]);

        cm->sefc[ee].meas = qdf[e].sface[0].meas;
        cm->sefc[ee].unitv[0] = qdf[e].sface[0].unitv[0];
        cm->sefc[ee].unitv[1] = qdf[e].sface[0].unitv[1];
        cm->sefc[ee].unitv[2] = qdf[e].sface[0].unitv[2];
        cm->sefc[ee+1].meas = qdf[e].sface[1].meas;
        cm->sefc[ee+1].unitv[0] = qdf[e].sface[1].unitv[0];
        cm->sefc[ee+1].unitv[1] = qdf[e].sface[1].unitv[1];
        cm->sefc[ee+1].unitv[2] = qdf[e].sface[1].unitv[2];

      } /* Loop on face edges */

    } /* (edga,face) quantities */

  } /* edge-->faces */

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
  const cs_quant_t  pfq = quant->face[f_id];

  fm->f_id = f_id;
  fm->face.meas = pfq.meas;
  for (int k = 0; k < 3; k++) {
    fm->face.center[k] = pfq.center[k];
    fm->face.unitv[k] = pfq.unitv[k];
  }

  const cs_lnum_t  *c2f_idx = connect->c2f->idx + c_id;
  const cs_lnum_t  *c2f_ids = connect->c2f->col_id + c2f_idx[0];
  const int  n_fc = c2f_idx[1] - c2f_idx[0];

  short int _f = n_fc;
  for (short int f = 0; f < n_fc; f++) {
    if (c2f_ids[f] == f_id) {

      const cs_nvec3_t  *deq = quant->dedge + c2f_idx[0];
      const short int  *f_sgn = connect->c2f->sgn + c2f_idx[0];

      _f = f;
      fm->dedge.meas = deq[f].meas;
      for (int k = 0; k < 3; k++) fm->dedge.unitv[k] = deq[f].unitv[k];
      fm->f_sgn = f_sgn[f];
      break;
    }
  }

  if (_f == n_fc) // Sanity check
    bft_error(__FILE__, __LINE__, 0,
              _(" Face %d not found.\n Stop build a face mesh."), f_id);

  const cs_lnum_t  *f2e_idx = connect->f2e->idx + f_id;
  const cs_lnum_t  *f2e_lst = connect->f2e->col_id + f2e_idx[0];

  fm->n_vf = fm->n_ef = f2e_idx[1] - f2e_idx[0];
  short int nv = 0;
  for (int i = 0; i < fm->n_vf; i++)
    fm->v_ids[i] = -1;

  for (short int e = 0; e < fm->n_ef; e++) {

    const cs_lnum_t  e_id = f2e_lst[e];
    const cs_quant_t  qe = quant->edge[e_id];

    fm->e_ids[e] = e_id;
    fm->edge[e].meas = qe.meas;
    for (int k = 0; k < 3; k++) {
      fm->edge[e].center[k] = qe.center[k];
      fm->edge[e].unitv[k] = qe.unitv[k];
    }

    const cs_lnum_t  *e2v_ids = connect->e2v->col_id + 2*e_id;
    short int  v1 = -1, v2 = -1;
    for (int v = 0; v < fm->n_vf && fm->v_ids[v] != -1; v++) {
      if (fm->v_ids[v] == e2v_ids[0])
        v1 = v;
      else if (fm->v_ids[v] == e2v_ids[1])
        v2 = v;
    }

    /* Add vertices if not already identified */
    if (v1 == -1) // Not found -> Add v1
      fm->v_ids[nv] = e2v_ids[0], v1 = nv++;
    if (v2 == -1) // Not found -> Add v2
      fm->v_ids[nv] = e2v_ids[1], v2 = nv++;

    /* Update e2v_ids */
    const int _eshft = 2*e;
    fm->e2v_ids[_eshft]   = v1;
    fm->e2v_ids[_eshft+1] = v2;

  } // Loop on face edges

  assert(nv == fm->n_vf); // Sanity check

  /* Update vertex coordinates */
  int  shift = 0;
  for (short int v = 0; v < fm->n_vf; v++) {
    const cs_real_t *xv = quant->vtx_coord + 3*fm->v_ids[v];
    for (int k = 0; k < 3; k++)
      fm->xv[shift++] = xv[k];
  }

  /* Define wvf and tef */
  for (int i = 0; i < fm->n_vf; i++)
    fm->wvf[i] = 0;

  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];

    /* tef = ||(xe -xf) x e||/2 = s(v1,e,f) + s(v2, e, f) */
    const double  tef = cs_compute_area_from_quant(fm->edge[e], pfq.center);

    fm->wvf[v1] += tef;
    fm->wvf[v2] += tef;
    fm->tef[e] = tef;

  }

  const double  invf = 0.5/pfq.meas;
  for (short int v = 0; v < fm->n_vf; v++) fm->wvf[v] *= invf;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_face_mesh_t structure for a given cell from a
 *         cs_cell_mesh_t structure
 *
 * \param[in]       cm        pointer to the reference cs_cell_mesh_t structure
 * \param[in]       f_id      face id in the cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure to set
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
  assert(cs_test_flag(cm->flag,
                      CS_CDO_LOCAL_PV  | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
                      CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV));

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

    const cs_lnum_t  eshft = 2*ec;
    cs_lnum_t  v1c_id = cm->e2v_ids[eshft];
    cs_lnum_t  v2c_id = cm->e2v_ids[eshft+1];

    /* Compact vertex numbering to this face */
    short int  v1 = -1, v2 = -1;
    for (int v = 0; v < fm->n_vf && fm->v_ids[v] != -1; v++) {
      if (fm->v_ids[v] == v1c_id)
        v1 = v;
      else if (fm->v_ids[v] == v2c_id)
        v2 = v;
    }

    /* Add vertices if not already identified */
    if (v1 == -1) // Not found -> Add v1
      fm->v_ids[nv] = v1c_id, v1 = nv++;
    if (v2 == -1) // Not found -> Add v2
      fm->v_ids[nv] = v2c_id, v2 = nv++;

    /* Update e2v_ids */
    const int _eshft = 2*ef;
    fm->e2v_ids[_eshft]   = v1;
    fm->e2v_ids[_eshft+1] = v2;

  } // Loop on face edges

  assert(nv == fm->n_vf); // Sanity check

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

END_C_DECLS
