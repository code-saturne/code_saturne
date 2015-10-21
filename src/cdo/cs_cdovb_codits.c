/*============================================================================
 * Build an algebraic CDO vertex-based system for convection/diffusion equation
 * with source terms
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

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_log.h"
#include "cs_search.h"
#include "cs_quadrature.h"
#include "cs_evaluate.h"
#include "cs_cdo_bc.h"
#include "cs_hodge.h"
#include "cs_cdo_convection.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovb_codits.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CDOVB_CODITS_DBG 0

/* Algebraic system for CDO vertex-based discretization */

struct _cs_cdovb_codits_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure. This structure is not owner. */

  const cs_equation_param_t  *eqp;

  /* System size (known boundary entities may be removed if BCs are strongly
     enforced) */

  cs_lnum_t  n_vertices;
  cs_lnum_t  n_dof_vertices; /* n_rows = n_cols = n_vertices - dir. vertices */

  /* Index of the system matrix (avoid to build it at each iteration)
     v2v connectivity through cell neighboorhood */
  cs_connect_index_t  *v2v;

  /* Boundary conditions:

     face_bc should not change during the simulation.
     The case of a definition of the BCs which changes of type during the
     simulation is possible but not implemented.
     You just have to call the initialization step each time the type of BCs
     is modified to define an updated cs_cdo_bc_t structure.

     We translate Dirichlet BCs to border vertices
     The values can be modified at each time step in case of transient
     simulation.
     For Neumann and Robin BCs, the treatment is more complex since the
     contributions of these BCs to a dual cell related to a border vertex is
     computed from the contribution of different portions of primal faces
     (these primal border faces form a closure of the dual cell).
     These contributions are computed on the fly.

   */

  cs_param_bc_enforce_t  enforce; // type of enforcement of BCs
  cs_cdo_bc_t           *face_bc; // list of faces sorted by type of BCs
  cs_cdo_bc_list_t      *vtx_dir; // list of vertices attached to a Dirichlet BC
  double                *dir_val; /* size = vtx_dir->n_nhmg_elts */

  /* Indirection between zipped numbering (without BC) and initial numbering
     Allocated only if the boundary conditions are strongly enforced.
  */
  cs_lnum_t          *v_z2i_ids;  // Mapping n_dof_vertices -> n_vertices
  cs_lnum_t          *v_i2z_ids;  // Mapping n_vertices     -> n_dof_vertices

  /* Work buffer */
  cs_lnum_t          *vtag;       /* size: n_vertices, store the local vertex id
                                     or -1 if not activated */
  cs_real_t          *source_terms;
  cs_real_t          *work;

};

typedef struct {

  cs_hodge_builder_t   *hb;  /* Pointer to a hodge builder structure */

  /* compact way to stored the edge --> vertices connect. related to a cell */
  int          n_bits;   // number of bits in a block
  int          n_blocks; // number of blocks in a mask
  cs_flag_t   *emsk;     // list of masks to store the connectivity

} _cdovb_diff_t;

/*============================================================================
 * Private variables
 *============================================================================*/

// Advanceed developper parameters
static const cs_real_t  cs_weak_nitsche_pena_coef = 500;
static const cs_real_t  cs_weak_penalization_weight = 0.01;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Encode E_c cap E_v sets for all vertices of a cell using a mask
 *          of bits
 *
 * \param[in]       c_id      cell id
 * \param[in]       connect   point to a cs_cdo_connect_t struct.
 * \param[in]       vtag      tag on vertices to store local ids
 * \param[in, out]  diff      pointer to a _cdovb_diff_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_encode_edge_masks(cs_lnum_t                  c_id,
                   const cs_cdo_connect_t    *connect,
                   const cs_lnum_t            vtag[],
                   _cdovb_diff_t             *diff)
{
  short int  ie;
  cs_lnum_t  i;

  cs_flag_t  *masks = diff->emsk;

  const int  n_bits = diff->n_bits;
  const cs_connect_index_t  *c2e = connect->c2e;

  /* Reset masks */
  for (i = 0; i < diff->n_blocks*connect->n_max_vbyc; i++)
    masks[i] = 0;

  /* Loop on cell edges */
  for (ie = 0, i = c2e->idx[c_id]; i < c2e->idx[c_id+1]; ie++, i++) {

    cs_lnum_t  eshft = 2*c2e->ids[i];
    short int  vi = vtag[connect->e2v->col_id[eshft]];
    short int  vj = vtag[connect->e2v->col_id[eshft+1]];

    /* Sanity checks */
    assert(vi > -1 && vi < connect->n_max_vbyc);
    assert(vj > -1 && vj < connect->n_max_vbyc);

    /* Encode this edge in the set E_v1,c and E_v2,c */
    short int  b = ie/n_bits, s = ie%n_bits;
    masks[vi*diff->n_blocks + b] |= (1 << s);
    masks[vj*diff->n_blocks + b] |= (1 << s);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a _cdovb_diff_t structure
 *
 * \param[in, out ] diff   pointer to a _cdovb_diff_t structure
 *
 * \return  NULL
 */
/*----------------------------------------------------------------------------*/

static _cdovb_diff_t *
_free_diffusion(_cdovb_diff_t   *diff)
{
  if (diff == NULL)
    return diff;

  BFT_FREE(diff->emsk);

  diff->hb = cs_hodge_builder_free(diff->hb);

  BFT_FREE(diff);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a structure used to build the stiffness matrix
 *
 * \param[in] connect     pointer to a cs_cdo_connect_t struct.
 * \param[in] time_step   pointer to a time step structure
 * \param[in] builder     pointer to a cs_cdovb_codits_t struct.
 *
 * \return a pointer to a new allocated _cdovb_diff_t structure
 */
/*----------------------------------------------------------------------------*/

static _cdovb_diff_t *
_init_diffusion(const cs_cdo_connect_t     *connect,
                const cs_time_step_t       *time_step,
                const cs_cdovb_codits_t    *builder)
{
  _cdovb_diff_t  *diff = NULL;

  const cs_equation_param_t  *eqp = builder->eqp;
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;

  /* Sanity check */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EPFD);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);

  BFT_MALLOC(diff, 1, _cdovb_diff_t);

  /* Initialize mask features */
  diff->n_bits = sizeof(cs_flag_t)*CHAR_BIT;
  diff->n_blocks = connect->n_max_ebyc/diff->n_bits;
  if (connect->n_max_ebyc % diff->n_bits != 0)
    diff->n_blocks += 1;

  BFT_MALLOC(diff->emsk, diff->n_blocks*connect->n_max_vbyc, cs_flag_t);
  for (int  i = 0; i < diff->n_blocks*connect->n_max_vbyc; i++)
    diff->emsk[i] = 0;

  /* Define a builder for the related discrete Hodge operator */
  diff->hb = cs_hodge_builder_init(connect, time_step, h_info);

  return diff;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the full stiffness matrix from a cellwise assembly process
 *
 * \param[in]      c_id        current cell id
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in]      builder     pointer to a cs_cdovb_codits_t struct.
 * \param[in, out] diff        auxiliary structure used to build the diff. term
 * \param[in, out] sloc        pointer to the local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static void
_build_local_stiffness_matrix(cs_lnum_t                   c_id,
                              const cs_cdo_connect_t     *connect,
                              const cs_cdo_quantities_t  *quant,
                              const cs_cdovb_codits_t    *builder,
                              _cdovb_diff_t              *diff,
                              cs_locmat_t                *sloc)
{
  short int  ek, el, vi, sgn_ik, sgn_jl;

  const cs_lnum_t  *vtag = builder->vtag;

  /* Build a local discrete Hodge operator and return a local dense matrix */
  const cs_locmat_t  *hloc = cs_hodge_build_local(c_id,
                                                  connect, quant,
                                                  diff->hb);

  /* Encode edge sets related to each vertex belonging to the current cell */
  _encode_edge_masks(c_id, connect, vtag, diff);

  /* Build local stiffness matrix */
  for (vi = 0; vi < sloc->n_ent; vi++) {

    const short int pos_i = vi*sloc->n_ent;

    for (ek = 0; ek < hloc->n_ent; ek++) {  /* Find edges attached to vi */

      const short int b = ek/diff->n_bits, r = ek % diff->n_bits;

      if (diff->emsk[vi*diff->n_blocks+b] & (1 << r)) { // ek in E_i

        const cs_lnum_t  ek_shft = 2*hloc->ids[ek];
        const int  pos_ek = ek*hloc->n_ent;

        if (connect->e2v->col_id[ek_shft] == sloc->ids[vi])
          sgn_ik = connect->e2v->sgn[ek_shft];
        else {
          assert(connect->e2v->col_id[ek_shft+1] == sloc->ids[vi]);
          sgn_ik = connect->e2v->sgn[ek_shft+1];
        }

        for (el = 0; el < hloc->n_ent; el++) { /* Loop on cell edges */

          cs_lnum_t  el_shft = 2*hloc->ids[el];
          double  val = hloc->mat[pos_ek+el]*sgn_ik;
          cs_lnum_t  v1_id = connect->e2v->col_id[el_shft];
          short int  vj1 = vtag[v1_id];

          sgn_jl = connect->e2v->sgn[el_shft];
          sloc->mat[pos_i+vj1] += val*sgn_jl;

          cs_lnum_t  v2_id = connect->e2v->col_id[el_shft+1];
          short int  vj2 = vtag[v2_id];

          sgn_jl = connect->e2v->sgn[el_shft+1];
          sloc->mat[pos_i+vj2] += val*sgn_jl;

        } /* Loop on cell edges */

      } /* ek in E_i */

    } /* Find edges attached to _vi */

  } /* Loop on vertices vi */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contribution of source terms to the rhs for this
 *          time step
 *
 * \param[in]      m            pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      builder      pointer to a cs_cdovb_codits_t structure
 * \param[in]      time_step    pointer to a time step structure
 * \param[in, out] full_rhs     right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_add_source_terms(const cs_mesh_t            *m,
                  const cs_cdo_connect_t     *connect,
                  const cs_cdo_quantities_t  *quant,
                  const cs_time_step_t       *time_step,
                  cs_cdovb_codits_t          *builder,
                  cs_real_t                   full_rhs[])
{
  cs_lnum_t  i;

  const cs_equation_param_t  *eqp = builder->eqp;

  if (eqp->flag & CS_EQUATION_UNSTEADY) {

    const cs_param_time_t  t_info = eqp->time_info;

    /* Previous values are stored inside builder->source_terms i.e.
       values of the source terms related to t_prev */
    if (t_info.scheme == CS_TIME_SCHEME_EXPLICIT)
      for (i = 0; i < builder->n_vertices; i++)
        full_rhs[i] += builder->source_terms[i];

    else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
             t_info.scheme == CS_TIME_SCHEME_THETA) {

      const double  tcoef = 1 - t_info.theta;

      for (i = 0; i < builder->n_vertices; i++)
        full_rhs[i] += tcoef * builder->source_terms[i];

    }

    /* Update builder->source_term with the value attached to t_cur */
    cs_cdovb_codits_compute_source(m, connect, quant, time_step, builder);

    if (t_info.scheme == CS_TIME_SCHEME_IMPLICIT)
      for (i = 0; i < builder->n_vertices; i++)
        full_rhs[i] += builder->source_terms[i];

    else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
             t_info.scheme == CS_TIME_SCHEME_THETA) {

      for (i = 0; i < builder->n_vertices; i++)
        full_rhs[i] += t_info.theta * builder->source_terms[i];

    }

  }
  else { /* Steady case: source terms have already been computed during
            the initialization step */

    for (i = 0; i < builder->n_vertices; i++)
      full_rhs[i] += builder->source_terms[i];

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the matrix related to the unsteady term
 *
 * \param[in]      builder    pointer to a cs_cdovb_codits_t structure
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_time_matrix(cs_cdovb_codits_t          *builder)
{
  cs_sla_matrix_t  *time_mat = NULL;

  const cs_equation_param_t  *eqp = builder->eqp;

  if (!eqp->flag & CS_EQUATION_UNSTEADY) // Steady-state eq. => Nothing to do
    return time_mat;

  const cs_param_hodge_t  h_info = eqp->time_hodge;

  /* Sanity check */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_VPCD);

  if (h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI)
    time_mat = cs_sla_matrix_create(builder->n_vertices, // n_rows
                                    builder->n_vertices, // n_cols
                                    1,                   // stride
                                    CS_SLA_MAT_MSR,      // type
                                    true);               // symmetric ?
  else
    time_mat = cs_sla_matrix_create_from_index(builder->v2v,
                                               CS_SLA_MAT_MSR,   // type = MSR
                                               true,             // sorted
                                               1);               // stride

  /* Set matrix flag */
  time_mat->flag |= CS_SLA_MATRIX_SYM;
  time_mat->flag |= CS_SLA_MATRIX_SORTED;

  return time_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the time discretization to all the terms of the equation
 *          excepted for the unsteady term and the source term
 *
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      time_mat   pointer to a cs_sla_matrix_t structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to a cs_cdovb_codits_t structure
 * \param[in, out] rhs        pointer to the right-hand side array
 * \param[in, out] matrix     pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_apply_time_scheme(const cs_real_t          *field_val,
                   const cs_sla_matrix_t    *time_mat,
                   double                    dt_cur,
                   cs_cdovb_codits_t        *builder,
                   cs_real_t                *rhs,
                   cs_sla_matrix_t          *sys_mat)

{
  cs_lnum_t  i;

  cs_real_t  *mv_time = builder->work;
  cs_real_t  *mv_sys = builder->work + builder->n_vertices;
  size_t  time_nnz = time_mat->idx[time_mat->n_rows];
  size_t  sys_nnz = sys_mat->idx[sys_mat->n_rows];

  const double inv_dt = 1/dt_cur;
  const cs_equation_param_t  *eqp = builder->eqp;
  const cs_param_time_t  t_info = eqp->time_info;

  /* Sanity checks */
  assert(eqp->flag & CS_EQUATION_UNSTEADY);
  assert(sys_mat->type == CS_SLA_MAT_MSR && time_mat->type == CS_SLA_MAT_MSR);
  assert(sys_mat->n_rows == builder->n_vertices);

  switch (t_info.scheme) {

  case CS_TIME_SCHEME_EXPLICIT:

    /* A_{diff,conv,reac} * p^(n) --> RHS */
    cs_sla_matvec(sys_mat, field_val, &mv_sys, true);
    cs_sla_matvec(time_mat, field_val, &mv_time, true);
    for (i = 0; i < builder->n_vertices; i++)
      rhs[i] += mv_sys[i] + inv_dt * mv_time[i];

    /* Implicit part is only the time_mat */
    for (i = 0; i < sys_mat->n_rows; i++)
      sys_mat->diag[i] = inv_dt * time_mat->diag[i];

    if (sys_nnz != time_nnz) {
      BFT_REALLOC(sys_mat->col_id, time_nnz, cs_lnum_t);
      BFT_REALLOC(sys_mat->val, time_nnz, cs_real_t);
    }
    memcpy(sys_mat->col_id, time_mat->col_id, time_nnz*sizeof(cs_lnum_t));

    for (i = 0; i < time_mat->n_rows; i++)
      sys_mat->idx[i+1] = time_mat->idx[i+1];
    for (size_t ii = 0; ii < time_nnz; ii++)
      sys_mat->val[ii] = inv_dt*time_mat->val[ii];

    break;

  case CS_TIME_SCHEME_IMPLICIT:

    /* Update rhs */
    cs_sla_matvec(time_mat, field_val, &mv_time, true);
    for (i = 0; i < builder->n_vertices; i++)
      rhs[i] += inv_dt * mv_time[i];

    /* Update the system matrix */
    for (i = 0; i < builder->n_vertices; i++)
      sys_mat->diag[i] += inv_dt * time_mat->diag[i];

    if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_WBS ||
        eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_COST) {

      /* sys_mat and time_hodge share the same index */
      assert(sys_nnz == time_nnz);
      for (size_t ii = 0; ii < sys_nnz; ii++)
        sys_mat->val[ii] += inv_dt * time_mat->val[ii];

    }
    break;

  case CS_TIME_SCHEME_CRANKNICO:
  case CS_TIME_SCHEME_THETA:
    {
      const double  tcoef = t_info.theta - 1;

      /* Update rhs += (tcoef*sys_mat + over_dt*time_hodge)*field_val */
      cs_sla_matvec(sys_mat, field_val, &mv_sys, true);
      cs_sla_matvec(time_mat, field_val, &mv_time, true);
      for (i = 0; i < builder->n_vertices; i++)
        rhs[i] += tcoef*mv_sys[i] + inv_dt*mv_time[i];

      /* Update the diag. terms of the system matrix */
      for (i = 0; i < builder->n_vertices; i++) {
        sys_mat->diag[i] *= t_info.theta;
        sys_mat->diag[i] += inv_dt * time_mat->diag[i];
      }

      /* Update the extra-diag. terms of the system matrix: sys_mat *= theta */
      for (size_t ii = 0; ii < sys_nnz; ii++)
        sys_mat->val[ii] *= t_info.theta;

      if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_WBS ||
          eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_COST) {

        /* sys_mat and time_hodge share the same index */
        assert(sys_nnz == time_nnz);

        for (size_t ii = 0; ii < sys_nnz; ii++)
          sys_mat->val[ii] += inv_dt * time_mat->val[ii];

      }

    }
    break;

  default:
    break;

  } // End of switch

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs
 *
 * \param[in]      m            pointer to a cs_mesh_t structure
 * \param[in]      time_step    pointer to a time step structure
 * \param[in]      field_val    pointer to the current value of the field
 * \param[in, out] builder      pointer to a cs_cdovb_codits_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_dir_values(const cs_mesh_t            *m,
                    const cs_time_step_t       *time_step,
                    const cs_real_t            *field_val,
                    const cs_cdovb_codits_t    *builder)
{
  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;
  const cs_equation_param_t  *eqp = builder->eqp;

  if (vtx_dir->n_nhmg_elts == 0)
    return; // Nothing to do

  cs_flag_t  dof_flag =
    CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL | CS_PARAM_FLAG_SCAL;

  /* Get the value of the Dirichlet for the current time */
  cs_cdo_bc_dirichlet_set(dof_flag,
                          time_step,
                          m,
                          eqp->bc,
                          vtx_dir,
                          builder->dir_val);

  if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
      builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
    if (eqp->flag & CS_EQUATION_UNSTEADY) {

      cs_lnum_t  i;

      const cs_param_time_t  t_info = eqp->time_info;

      /* Previous values of the unknown are stored inside field_val (iter n) */
      if (t_info.scheme == CS_TIME_SCHEME_EXPLICIT)
        for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
          builder->dir_val[i] = field_val[vtx_dir->elt_ids[i]];

      else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
               t_info.scheme == CS_TIME_SCHEME_THETA) {

        const double  tcoef = 1 - t_info.theta;

        for (i = 0; i < vtx_dir->n_nhmg_elts; i++) {
          builder->dir_val[i] *= t_info.theta;
          builder->dir_val[i] += tcoef * field_val[vtx_dir->elt_ids[i]];
        }

      }

    } /* Unsteady */

  } /* Enforcement is not strong or penalized */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Update the matrix of the linear system to take into account the
 *          boundary contribution of the convection term
 *
 * \param[in]      m            pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in, out] conv         pointer to a convection builder structure
 * \param[in, out] sys_builder  pointer to a cs_cdovb_codits_t structure
 * \param[in, out] rhs          pointer of pointer to the right-hand side
 * \param[in, out] matrix       pointer to a matrix structure
 */
/*----------------------------------------------------------------------------*/

static void
_add_advection_bc(const cs_mesh_t            *m,
                  const cs_cdo_connect_t     *connect,
                  const cs_cdo_quantities_t  *quant,
                  cs_convection_builder_t    *conv,
                  cs_cdovb_codits_t          *builder,
                  cs_real_t                  *rhs,
                  cs_sla_matrix_t            *matrix)
{
  cs_lnum_t  i;

  /* temporary buffer */
  cs_real_t  *dir_vals = builder->work;
  cs_real_t  *rhs_contrib = builder->work + builder->n_vertices;
  cs_real_t  *diag_contrib = builder->work + builder->n_vertices;

  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;

  /* Initialize arrays */
  for (i = 0; i < builder->n_vertices; i++)
    dir_vals[i] = rhs_contrib[i] = diag_contrib[i] = 0;

  /* Store the value of the Dirichlet condition into a full (i.e. n_vertices)
     array */
  for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
    dir_vals[vtx_dir->elt_ids[i]] = builder->dir_val[i];

  /* Compute the contribution */
  cs_convection_get_bc_contrib(m, connect, quant, dir_vals,
                               conv,
                               rhs_contrib,
                               diag_contrib);

  /* Add boundary contribution */
  for (i = 0; i < builder->n_vertices; i++) {
    rhs[i] += rhs_contrib[i];
    matrix->diag[i] += diag_contrib[i];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply Dirichlet boundary conditions arising from the diffusion op.
 *          Only useful if boundary conditions are weakly imposed using a
 *          Nitsche technique.
 *          Right-hand side and the system matrix are updated.
 *
 * \param[in]      m            pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step    pointer to a time step structure
 * \param[in, out] builder      pointer to a cs_cdovb_codits_t structure
 * \param[in, out] full_rhs     right-hand side
 * \param[in, out] full_matrix  matrix of the linear system
 */
/*----------------------------------------------------------------------------*/

static void
_add_diffusion_bc(const cs_mesh_t            *m,
                  const cs_cdo_connect_t     *connect,
                  const cs_cdo_quantities_t  *quant,
                  const cs_time_step_t       *time_step,
                  cs_cdovb_codits_t          *builder,
                  cs_real_t                   full_rhs[],
                  cs_sla_matrix_t            *full_matrix)
{
  if (builder->enforce !=  CS_PARAM_BC_ENFORCE_WEAK_SYM ||
      builder->enforce !=  CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
    return; // Nothing to do

  const double  tcur = time_step->t_cur;
  const cs_equation_param_t  *eqp = builder->eqp;
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;

  /* Sanity check */
  assert(full_matrix->type == CS_SLA_MAT_MSR);

  /* Sanity check */
  if (h_info.algo != CS_PARAM_HODGE_ALGO_COST)
    bft_error(__FILE__, __LINE__, 0,
              " Only COST algorithm for hodge operators is possible when\n"
              " a weak enforcement of the boundary conditions is requested.");

  cs_lnum_t  i, ie, j, k, _id, v_id;
  cs_real_t  surf, eig_ratio, eig_max, pena_coef;
  cs_real_3_t  xyz, mn, reco_val;
  cs_real_33_t  matpty;

  short int  *loc_e_id = NULL;
  cs_real_t  *over_pec_vol = NULL, *_vec = NULL, *_matvec = NULL;
  cs_real_3_t  *dfv = NULL, *pev = NULL; // dual face and primal edge vectors
  cs_locmat_t  *transp = NULL;

  cs_locmat_t  *ntrgrd = cs_locmat_create(connect->n_max_vbyc);

  /* Temporary buffers stored using builder */
  cs_lnum_t  *loc_v_id = builder->vtag;
  cs_real_t  *dir_vals = builder->work;
  cs_real_t  *v_coef = builder->work + builder->n_vertices;

  /* Initialize v_coef and loc_v_ids */
  for (i = 0; i < builder->n_vertices; i++) {
    v_coef[i] = 0.0;
    loc_v_id[i] = -1;
  }

  BFT_MALLOC(loc_e_id, quant->n_edges, short int);
  for (i = 0; i < quant->n_edges; i++)
    loc_e_id[i] = -1;

  /* Store locally dual face vector for a quicker access */
  BFT_MALLOC(dfv, connect->n_max_ebyc, cs_real_3_t);
  BFT_MALLOC(pev, connect->n_max_ebyc, cs_real_3_t);
  BFT_MALLOC(over_pec_vol, connect->n_max_ebyc, cs_real_t);

  if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

    for (i = 0; i < builder->n_vertices; i++)
      dir_vals[i] = 0.;
    for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
      dir_vals[vtx_dir->elt_ids[i]] = builder->dir_val[i];

    BFT_MALLOC(_vec, connect->n_max_vbyc, cs_real_t);
    BFT_MALLOC(_matvec, connect->n_max_vbyc, cs_real_t);

    transp = cs_locmat_create(connect->n_max_vbyc);
  }

  /* Parameters linked to the discrete Hodge operator used for diffusion */
  bool  is_uniform = cs_param_pty_is_uniform(h_info.pty_id);

  /* Get the anisotropic ratio and the max. eigenvalue (if uniform) */
  if (is_uniform) {
    xyz[0] = xyz[1] = xyz[2] = 0.;
    cs_evaluate_pty(h_info.pty_id,
                    tcur,       // when ?
                    xyz,        // where ?
                    h_info.inv_pty,
                    &matpty);

    cs_eigen_mat33(matpty, &eig_ratio, &eig_max);
  }

  const double  beta = h_info.coef;
  const cs_sla_matrix_t  *e2v = connect->e2v;
  const cs_sla_matrix_t  *f2e = connect->f2e;
  const cs_cdo_bc_list_t  *face_dir = builder->face_bc->dir;
  const cs_connect_index_t  *c2e = connect->c2e;
  const cs_connect_index_t  *c2v = connect->c2v;

  /* Loop on Dirichlet faces */
  for (i = 0; i < face_dir->n_elts; i++) {

    cs_lnum_t  f_id = face_dir->elt_ids[i] + quant->n_i_faces;
    cs_quant_t  qf = quant->face[f_id];
    cs_lnum_t  c_id = connect->f2c->col_id[connect->f2c->idx[f_id]];
    cs_real_t  over_c_vol = 1/quant->cell_vol[c_id];

    /* Sanity check (this is a border face) */
    assert(connect->f2c->idx[f_id+1] - connect->f2c->idx[f_id] == 1);

    if (!is_uniform) {
      cs_evaluate_pty(h_info.pty_id,
                      tcur,                           // when ?
                      &(quant->cell_centers[3*c_id]), // where ?
                      h_info.inv_pty,
                      &matpty);

      cs_eigen_mat33(matpty, &eig_ratio, &eig_max);
    }

    cs_real_t  f_coef = pow(qf.meas, -0.5) * eig_ratio * eig_max;

    /* Compute the product: matpty*face unit normal */
    _mv3(matpty, qf.unitv, &mn);

    /* Define an id local to this cell for each vertex */
    for (j = c2v->idx[c_id], _id = 0; j < c2v->idx[c_id+1]; j++, _id++) {
      v_id = c2v->ids[j];
      loc_v_id[v_id] = _id;
      ntrgrd->ids[_id] = v_id;
    }

    /* Initialize buffers related to vertices */
    ntrgrd->n_ent = _id;

    /* Reset local trace normal*gradient operator */
    for (j = 0; j < ntrgrd->n_ent*ntrgrd->n_ent; j++)
      ntrgrd->mat[j] = 0.0;

    /* Store the local id and useful quantities for each edge */
    for (j = c2e->idx[c_id], _id = 0; j < c2e->idx[c_id+1]; j++, _id++) {

      cs_lnum_t  e_id = c2e->ids[j];
      cs_quant_t qpe = quant->edge[e_id];
      cs_dface_t qdf = quant->dface[j];

      loc_e_id[e_id] = _id;
      for (k = 0; k < 3; k++) {
        pev[_id][k] = qpe.unitv[k]*qpe.meas;
        dfv[_id][k] = qdf.vect[k];
      }
      over_pec_vol[_id] = 3./_dp3(pev[_id], dfv[_id]);

    } // Loop on cell edges

    /* Loop on border face edges */
    for (j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

      cs_lnum_t  e_id = f2e->col_id[j];
      cs_quant_t  qe = quant->edge[e_id];
      cs_lnum_t  e_shft = e2v->idx[e_id];
      cs_lnum_t  v1_id = e2v->col_id[e_shft];
      cs_lnum_t  v2_id = e2v->col_id[e_shft+1];

      short int  _e = loc_e_id[e_id];
      short int  _v1 = loc_v_id[v1_id];
      short int  _v2 = loc_v_id[v2_id];

      /* Sanity checks */
      assert(_e != -1 && _v1 != -1 && _v2 != -1);

      surf = cs_surftri(&(m->vtx_coord[3*v1_id]), qe.center, qf.center);
      v_coef[v1_id] += surf*f_coef;
      v_coef[v2_id] += surf*f_coef;

      for (ie = c2e->idx[c_id]; ie < c2e->idx[c_id+1]; ie++) {

        cs_lnum_t  ek_id = c2e->ids[ie];
        cs_lnum_t  ek_shft = e2v->idx[ek_id];
        cs_lnum_t  vj1_id = e2v->col_id[ek_shft];
        short int  sgn_j1 = e2v->sgn[ek_shft];
        cs_lnum_t  vj2_id = e2v->col_id[ek_shft+1];
        short int  sgn_j2 = e2v->sgn[ek_shft+1];

        short int  _ek = loc_e_id[ek_id];
        short int  _vj1 = loc_v_id[vj1_id];
        short int  _vj2 = loc_v_id[vj2_id];

        /* Compute L_Ec(GRAD(p_j)) on t_{ek,f} */
        if (_ek == _e)
          for (k = 0; k < 3; k++)
            reco_val[k] = dfv[_ek][k]*beta*over_pec_vol[_ek];
        else
          for (k = 0; k < 3; k++)
            reco_val[k] = 0.0;

        cs_real_t  dp = beta*over_pec_vol[_e]*_dp3(dfv[_ek], pev[_e]);

        for (k = 0; k < 3; k++)
          reco_val[k] += over_c_vol*(dfv[_ek][k] - dp*dfv[_e][k]);

        cs_real_t  contrib = _dp3(mn, reco_val)*surf;

        /* Update the coefficients of the local normal trace operator */
        ntrgrd->mat[_v1*ntrgrd->n_ent + _vj1] += sgn_j1 * contrib;
        ntrgrd->mat[_v1*ntrgrd->n_ent + _vj2] += sgn_j2 * contrib;
        ntrgrd->mat[_v2*ntrgrd->n_ent + _vj1] += sgn_j1 * contrib;
        ntrgrd->mat[_v2*ntrgrd->n_ent + _vj2] += sgn_j2 * contrib;

      } // Loop on cell edges

    } // border face edges

    /* Assemble contributions coming from local normal trace operator */
    if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

      /* ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) */
      cs_locmat_add_transpose(ntrgrd, transp);

      cs_sla_assemble_msr_sym(ntrgrd, full_matrix, false); // Not only diagonal

      /* Modify RHS according to the add of transp */
      for (j = 0; j < ntrgrd->n_ent; j++)
        _vec[j] = dir_vals[ntrgrd->ids[j]];

      cs_locmat_matvec(transp, _vec, _matvec);

      for (j = 0; j < ntrgrd->n_ent; j++)
        full_rhs[ntrgrd->ids[j]] += _matvec[j];

    }
    else
      cs_sla_assemble_msr(ntrgrd, full_matrix);


  } // Dirichlet faces

  /* Loop on Dirichlet vertices (first, nhmg vertices then hmg vertices)
     to assemble the contribution coming from the diagonal Hodge penalization */
  for (i = 0; i < vtx_dir->n_elts; i++) {

    v_id = vtx_dir->elt_ids[i];
    pena_coef = cs_weak_nitsche_pena_coef*v_coef[v_id];
    full_matrix->diag[v_id] += pena_coef;
    if (i < vtx_dir->n_nhmg_elts)
      full_rhs[v_id] += pena_coef*builder->dir_val[i];

  }

  /* Free memory */
  if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
    transp = cs_locmat_free(transp);
    BFT_FREE(_vec);
    BFT_FREE(_matvec);
  }

  ntrgrd = cs_locmat_free(ntrgrd);

  BFT_FREE(loc_e_id);
  BFT_FREE(dfv);
  BFT_FREE(pev);
  BFT_FREE(over_pec_vol);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply boundary conditions. Update right-hand side and the system
 *          matrix
 *
 * \param[in, out] builder      pointer to a cs_cdovb_codits_t structure
 * \param[in, out] full_matrix  matrix of the linear system
 * \param[in, out] full_rhs     right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_strong_bc_enforcement(cs_cdovb_codits_t       *builder,
                       cs_real_t              **rhs,
                       cs_sla_matrix_t        **matrix)
{
  cs_lnum_t  i;

  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;

  if (vtx_dir->n_nhmg_elts == 0)
    return;

  /* Sanity check */
  assert(builder->n_vertices > builder->n_dof_vertices);

  cs_sla_matrix_t  *full_matrix = *matrix, *final_matrix = NULL;
  double  *full_rhs = *rhs, *final_rhs = NULL;
  double  *tmp_rhs = builder->work;
  double  *x_bc = builder->work + builder->n_vertices;
  double  *contrib = builder->work + 2*builder->n_vertices;

  for (i = 0; i < builder->n_vertices; i++)
    x_bc[i] = 0.0;
  for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
    x_bc[vtx_dir->elt_ids[i]] = builder->dir_val[i];

  /* Compute full_matrix*Tbc: rhs = rhs - full_matrix*Tbc */
  cs_sla_matvec(full_matrix, x_bc, &contrib, true);
  for (i = 0; i < builder->n_vertices; i++)
    full_rhs[i] -= contrib[i];

  /* Reduce the rhs size. Remove vertices with Dirichlet BC */
  memcpy(tmp_rhs, full_rhs, builder->n_vertices*sizeof(double));
  final_rhs = *rhs;
  for (i = 0; i < builder->n_dof_vertices; i++)
    final_rhs[i] = tmp_rhs[builder->v_z2i_ids[i]];

  /* Reduce the system size from n_vertices to n_dof_vertices.
     Vertices attached to a Dirichlet BC are removed.
     Extract block with degrees of freedom */
  final_matrix = cs_sla_matrix_pack(builder->n_dof_vertices,
                                    builder->n_dof_vertices,
                                    full_matrix,
                                    builder->v_z2i_ids,
                                    builder->v_i2z_ids,
                                    true); // keep sym.

  /* Free buffers */
  full_matrix = cs_sla_matrix_free(full_matrix);

  /* Return pointers */
  *matrix = final_matrix;
  *rhs = final_rhs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Modify the matrix of the linear system and its right hand side to
 *          take into account a strong enforcement or a large penalization of
 *          the boundary conditions.
 *          Nothing to do in case of weak enforcement.
 *
 * \param[in, out] builder      pointer to a cs_cdovb_codits_t structure
 * \param[in, out] full_matrix  matrix of the linear system
 * \param[in, out] full_rhs     right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_enforce_bc(cs_cdovb_codits_t          *builder,
            cs_real_t                 **rhs,
            cs_sla_matrix_t           **matrix)
{
  int  i;

  cs_sla_matrix_t  *full_matrix = *matrix, *final_matrix = NULL;
  double  *full_rhs = *rhs, *final_rhs = NULL;

  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;

  /* Sanity check */
  if (builder->enforce != CS_PARAM_BC_ENFORCE_STRONG &&
      builder->n_vertices != builder->n_dof_vertices)
    bft_error(__FILE__, __LINE__, 0,
              " Error detected: Boundary conditions are not strongly enforced"
              " but there are some removed vertices.");

  /* Treatment differs according to the way to enforce BCs.
     In vertex-based scheme, Dirichlet BC are essential and Neuman BC natural */
  switch (builder->enforce) {

  case CS_PARAM_BC_ENFORCE_STRONG:
    _strong_bc_enforcement(builder, rhs, matrix);
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_PENA:
    {
      cs_real_t  pena_coef = cs_weak_penalization_weight/cs_get_eps_machine();

      for (i = 0; i < builder->vtx_dir->n_elts; i++)
        full_matrix->diag[builder->vtx_dir->elt_ids[i]] += pena_coef;

      for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
        full_rhs[vtx_dir->elt_ids[i]] += pena_coef * builder->dir_val[i];

    }
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    // Nothing to do
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " This kind of BC enforcement is not implemented yet.\n"
              " Please modify your settings.");

  } // End of switch (enforcement)

  /* TODO: Add contribution for Neumann BC (if homogeneous nothing to do)
     and Robin BC */

  if (builder->n_vertices == builder->n_dof_vertices) { // Keep the full system
    final_matrix = full_matrix;
    final_rhs = full_rhs;
  }

  /* Return pointers */
  *matrix = final_matrix;
  *rhs = final_rhs;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_codits_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 * \param[in] connect   pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cdovb_codits_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_codits_init(const cs_equation_param_t  *eqp,
                     const cs_mesh_t            *mesh,
                     const cs_cdo_connect_t     *connect)
{
  cs_lnum_t  i;

  /* Sanity checks */
  assert(eqp != NULL);
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(eqp->type == CS_EQUATION_TYPE_SCAL);

  const cs_lnum_t  n_vertices = connect->v_info->n_ent;
  const cs_lnum_t  n_b_faces = connect->f_info->n_ent_bd;

  cs_cdovb_codits_t  *builder = NULL;

  BFT_MALLOC(builder, 1, cs_cdovb_codits_t);

  builder->eqp = eqp;

  /* Dimensions:
     By default, we set number of DoFs as if there is a weak enforcement of
     the boundary conditions */
  builder->n_vertices = n_vertices;
  builder->n_dof_vertices = n_vertices;

  /* Build a (sorted) v2v connectivity index */
  const cs_connect_index_t  *c2v = connect->c2v;
  cs_connect_index_t  *v2c = cs_index_transpose(n_vertices, c2v);

  builder->v2v = cs_index_compose(n_vertices, v2c, c2v);
  cs_index_sort(builder->v2v);
  cs_index_free(&v2c);

  /* Update index (v2v has a diagonal entry. We remove it since we consider a
     matrix stored using the MSR format */
  cs_lnum_t  shift = 0;
  cs_lnum_t  *tmp_idx = NULL;

  BFT_MALLOC(tmp_idx, n_vertices + 1, cs_lnum_t);
  memcpy(tmp_idx, builder->v2v->idx, sizeof(cs_lnum_t)*(n_vertices+1));

  for (i = 0; i < n_vertices; i++) {

    cs_lnum_t  start = tmp_idx[i], end = tmp_idx[i+1];

    for (cs_lnum_t  j = start; j < end; j++)
      if (builder->v2v->ids[j] != i)
        builder->v2v->ids[shift++] = builder->v2v->ids[j];

    builder->v2v->idx[i+1] = builder->v2v->idx[i] + end-start-1;

  } // Loop on vertices

  assert(shift == builder->v2v->idx[n_vertices]); // Sanity check
  BFT_FREE(tmp_idx);

  /* Set members and structures related to the management of the BCs */

  const cs_param_bc_t  *bc_param = eqp->bc;

  /* Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.
     We also compute also the list of Dirichlet vertices along with their
     related definition.
  */
  builder->face_bc = cs_cdo_bc_init(bc_param, n_b_faces);
  builder->vtx_dir = cs_cdo_bc_vtx_dir_create(mesh, builder->face_bc);

  /* Allocate dir_val */
  BFT_MALLOC(builder->dir_val, builder->vtx_dir->n_nhmg_elts, double);
  for (i = 0; i < builder->vtx_dir->n_nhmg_elts; i++)
    builder->dir_val[i] = 0.0;

  /* Strong enforcement means that we need an indirection list between the
     compress (or zip) and initial numbering of vertices */
  builder->enforce = bc_param->enforcement;
  builder->v_z2i_ids = NULL; // zipped --> initial ids
  builder->v_i2z_ids = NULL; // initial --> zipped ids

  if (builder->enforce == CS_PARAM_BC_ENFORCE_STRONG &&
      builder->vtx_dir->n_elts > 0) {

    const bool  advection = (eqp->flag & CS_EQUATION_CONVECTION) ? true : false;
    const bool  unsteady  = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;

    if (advection || unsteady)
      bft_error(__FILE__, __LINE__, 0,
                " Invalid choice of enforcement of the boundary conditions.\n"
                " Strong enforcement is not implemented when convection or"
                " unsteady terms are activated.\n"
                " Please modify your settings.");

    cs_lnum_t  cur_id = 0;
    bool  *is_kept = NULL;

    builder->n_dof_vertices = n_vertices - builder->vtx_dir->n_elts;

    BFT_MALLOC(is_kept, n_vertices, bool);
    for (i = 0; i < n_vertices; i++)
      is_kept[i] = true;
    for (i = 0; i < builder->vtx_dir->n_elts; i++)
      is_kept[builder->vtx_dir->elt_ids[i]] = false;

    /* Build builder->v_z2i_ids and builder->i2i_ids */
    BFT_MALLOC(builder->v_z2i_ids, builder->n_dof_vertices, cs_lnum_t);
    BFT_MALLOC(builder->v_i2z_ids, builder->n_vertices, cs_lnum_t);

    for (i = 0; i < builder->n_vertices; i++) {
      /* by default, we consider that it's removed */
      builder->v_i2z_ids[i] = -1;
      if (is_kept[i]) {
        builder->v_i2z_ids[i] = cur_id;
        builder->v_z2i_ids[cur_id++] = i;
      }
    }
    assert(cur_id == builder->n_dof_vertices);

    BFT_FREE(is_kept);

  } /* Strong enforcement of BCs */

  /* Work buffer */
  BFT_MALLOC(builder->source_terms, builder->n_vertices, cs_real_t);

  /* Initialize tags */
  BFT_MALLOC(builder->vtag, n_vertices, cs_lnum_t);
  for (i = 0; i < n_vertices; i++)
    builder->vtag[i] = -1;

  BFT_MALLOC(builder->work, 3*n_vertices, cs_real_t);

  return builder;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovb_codits_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovb_codits_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_codits_free(void   *builder)
{
  cs_cdovb_codits_t  *_builder = (cs_cdovb_codits_t *)builder;

  if (_builder == NULL)
    return _builder;

  /* Free BC structure */
  if (_builder->vtx_dir->n_nhmg_elts > 0)
    BFT_FREE(_builder->dir_val);

  _builder->face_bc = cs_cdo_bc_free(_builder->face_bc);
  _builder->vtx_dir = cs_cdo_bc_list_free(_builder->vtx_dir);

  /* Free connectivity index */
  cs_index_free(&(_builder->v2v));

  /* Renumbering (if strong enforcement of BCs for instance) */
  if (_builder->n_vertices > _builder->n_dof_vertices) {
    BFT_FREE(_builder->v_z2i_ids);
    BFT_FREE(_builder->v_i2z_ids);
  }

  BFT_FREE(_builder->vtag);
  BFT_FREE(_builder->source_terms);
  BFT_FREE(_builder->work);
  BFT_FREE(_builder);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in]      m           pointer to a cs_mesh_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step   pointer to a time step structure
 * \param[in, out] builder     pointer to a cs_cdovb_codits_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_compute_source(const cs_mesh_t            *m,
                               const cs_cdo_connect_t     *connect,
                               const cs_cdo_quantities_t  *quant,
                               const cs_time_step_t       *time_step,
                               void                       *builder)
{
  cs_lnum_t  i;
  cs_cdovb_codits_t  *bld = (cs_cdovb_codits_t *)builder;

  const cs_equation_param_t  *eqp = bld->eqp;

  if (eqp->n_source_terms > 0) { /* Add contribution from source terms */

    for (i = 0; i < bld->n_vertices; i++)
      bld->source_terms[i] = 0;

    for (int  st_id = 0; st_id < eqp->n_source_terms; st_id++) {

      const cs_param_source_term_t  st = eqp->source_terms[st_id];

      double  *contrib = bld->work + bld->n_vertices;
      cs_flag_t  dof_flag =
        CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_DUAL | CS_PARAM_FLAG_SCAL;

      /* Sanity check */
      assert(st.var_type == CS_PARAM_VAR_SCAL);

      cs_evaluate(m, quant, connect,  // geometrical and topological info.
                  time_step,
                  dof_flag,
                  st.ml_id,
                  st.def_type,
                  st.quad_type,
                  st.use_subdiv,
                  st.def,             // definition of the explicit part
                  &contrib);          // updated inside this function

      /* Update source term array */
      for (i = 0; i < bld->n_vertices; i++)
        bld->source_terms[i] += contrib[i];

    } // Loop on source terms

  } /* There is at least one source term which is defined */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO vertex-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      m          pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to cs_cdovb_codits_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] sla_mat    pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_build_system(const cs_mesh_t             *m,
                             const cs_cdo_connect_t      *connect,
                             const cs_cdo_quantities_t   *quant,
                             const cs_real_t             *field_val,
                             const cs_time_step_t        *time_step,
                             double                       dt_cur,
                             void                        *builder,
                             cs_real_t                  **rhs,
                             cs_sla_matrix_t            **sla_mat)
{
  cs_lnum_t  i, j;

  cs_cdovb_codits_t  *sys_builder = (cs_cdovb_codits_t *)builder;

  const cs_connect_index_t  *c2v = connect->c2v;
  const cs_equation_param_t  *eqp = sys_builder->eqp;
  const bool  do_diffusion = (eqp->flag & CS_EQUATION_DIFFUSION) ? true : false;
  const bool  do_advection = (eqp->flag & CS_EQUATION_CONVECTION) ? true : false;
  const bool  do_unsteady  = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;

  /* Allocate and initialize a matrix with the larger stencil (that related
     to diffusion => all vertices of a cell are potentially in interaction)
     adr = advection/diffusion/reaction
  */
  cs_sla_matrix_t  *sys_mat =
    cs_sla_matrix_create_from_index(sys_builder->v2v,
                                    CS_SLA_MAT_MSR,   // type = MSR
                                    true,             // sorted
                                    1);               // stride

  if (!do_advection && sys_builder->enforce != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
    sys_mat->flag |= CS_SLA_MATRIX_SYM;

  cs_locmat_t  *adr_loc = cs_locmat_create(connect->n_max_vbyc);

  /* Preparatory step for diffusion term */
  _cdovb_diff_t  *diff = NULL;
  if (do_diffusion)
    diff = _init_diffusion(connect, time_step, sys_builder);

  /* Preparatory step for advection term */
  cs_convection_builder_t  *conv = NULL;
  if (do_advection)
    conv = cs_convection_builder_init(connect, time_step,
                                      eqp->advection,
                                      do_diffusion,
                                      eqp->diffusion_hodge);

  /* Preparatory step for unsteady term */
  cs_sla_matrix_t  *time_mat = _init_time_matrix(sys_builder);
  cs_hodge_builder_t  *time_builder = NULL;
  if (do_unsteady)
    time_builder = cs_hodge_builder_init(connect, time_step, eqp->time_hodge);

  const bool  only_time_diag =
    (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_VORONOI) ? true : false;

  /* Main loop on cells to build the system matrix */
  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Initialize local system matrix */
    adr_loc->n_ent = c2v->idx[c_id+1] - c2v->idx[c_id];
    for (i = c2v->idx[c_id], j = 0; i < c2v->idx[c_id+1]; i++, j++) {
      adr_loc->ids[j] = c2v->ids[i];
      sys_builder->vtag[adr_loc->ids[j]] = j;
    }

    for (i = 0; i < adr_loc->n_ent*adr_loc->n_ent; i++)
      adr_loc->mat[i] = 0;

    /* Define the local matrix related to diffusion/convection/reaction */
    if (do_diffusion)
      _build_local_stiffness_matrix(c_id,
                                    connect, quant, sys_builder,
                                    diff, adr_loc);

    if (do_advection) {
      cs_locmat_t  *adv_loc = cs_convection_build_local_op(c_id,
                                                           connect, quant,
                                                           sys_builder->vtag,
                                                           conv);

      cs_locmat_add(adr_loc, adv_loc);
    }

    if (do_unsteady) { /* Build mass matrix to take into account time effect */

      cs_locmat_t  *time_loc = cs_hodge_build_local(c_id,
                                                    connect, quant,
                                                    time_builder);

      /* Assemble the local matrix into the system matrix */
      cs_sla_assemble_msr_sym(time_loc, time_mat, only_time_diag);

    }

    /* Assemble the matrix related to the advcetion/diffusion/reaction terms
       If advection is activated, the resulting system is not symmetric
       Otherwise, the system is symmetric with extra-diagonal terms. */
    if (sys_mat->flag & CS_SLA_MATRIX_SYM)
      cs_sla_assemble_msr_sym(adr_loc, sys_mat, false);
    else
      cs_sla_assemble_msr(adr_loc, sys_mat);

  } // Main loop on cells

  /* Free temporary buffers and structures */
  adr_loc = cs_locmat_free(adr_loc);

  if (do_unsteady)
    time_builder = cs_hodge_builder_free(time_builder);

  /* Initialize full rhs */
  cs_real_t  *full_rhs = *rhs;
  if (full_rhs == NULL)
    BFT_MALLOC(full_rhs, sys_builder->n_vertices, cs_real_t);
  for (i = 0; i < sys_builder->n_vertices; i++)
    full_rhs[i] = 0.0;

  /* Compute the contribution of source terms to the full rhs for this
     time step */
  _add_source_terms(m, connect, quant, time_step, sys_builder, full_rhs);

  /* Compute the values of the Dirichlet BC.
     TODO: do the analogy for Neumann BC */
  _compute_dir_values(m, time_step, field_val, sys_builder);

  /* Apply boundary conditions */
  if (do_diffusion)
    _add_diffusion_bc(m, connect, quant, time_step,
                      sys_builder, full_rhs, sys_mat);

  if (do_advection)
    _add_advection_bc(m, connect, quant,
                      conv, sys_builder, full_rhs, sys_mat);

  /* Unsteady terms have to be considered at the end in order to deal with
     the system matrix fullfilled with diffusion, convection and reaction terms
  */
  if (do_unsteady) {

    _apply_time_scheme(field_val, time_mat, dt_cur,
                       sys_builder, full_rhs, sys_mat);

    /* sys_mat becomes the system matrix since it now includes the time
       contribution */
    time_mat = cs_sla_matrix_free(time_mat);
  }

  /* Final step in BC management.
     Apply the strong or penalized enforcement.
     Must be call after the application of the time scheme */
  _enforce_bc(sys_builder, &full_rhs, &sys_mat);

  bool do_cleaning = false; // Advanced option
  if (do_cleaning)
    cs_sla_matrix_clean(sys_mat, 10*cs_get_eps_machine());

  /* Free remaining buffers and structures */
  if (do_diffusion)
    diff = _free_diffusion(diff);
  if (do_advection)
    conv = cs_convection_builder_free(conv);

  /* Return pointers */
  *rhs = full_rhs;
  *sla_mat = sys_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO vertex-based scheme.
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      solu       solution array
 * \param[in, out] builder    pointer to cs_cdovb_codits_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_update_field(const cs_cdo_connect_t     *connect,
                             const cs_cdo_quantities_t  *quant,
                             const cs_time_step_t       *time_step,
                             const cs_real_t            *solu,
                             void                       *builder,
                             cs_real_t                  *field_val)
{
  int  i;

  cs_cdovb_codits_t  *_builder = (cs_cdovb_codits_t  *)builder;
  const cs_cdo_bc_list_t  *v_dir = _builder->vtx_dir;

  /* Sanity checks (avoid a warning in debug mode) */
  assert(connect != NULL);
  assert(quant != NULL);
  assert(time_step != NULL);

  /* Set computed solution in field array */
  if (_builder->n_dof_vertices < _builder->n_vertices) {
    for (i = 0; i < _builder->n_vertices; i++)
      field_val[i] = 0;
    for (i = 0; i < _builder->n_dof_vertices; i++)
      field_val[_builder->v_z2i_ids[i]] = solu[i];
  }
  else
    for (i = 0; i < _builder->n_vertices; i++)
      field_val[i] = solu[i];

  /* Set BC in field array if we have this knowledge */
  if (_builder->enforce == CS_PARAM_BC_ENFORCE_STRONG)
    for (i = 0; i < v_dir->n_nhmg_elts; i++)
      field_val[v_dir->elt_ids[i]] = _builder->dir_val[i];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
