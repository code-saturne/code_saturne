/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
 * reaction scalar equations with source terms
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_log.h"
#include "cs_math.h"
#include "cs_search.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_source_term.h"
#include "cs_cdo_bc.h"
#include "cs_hodge.h"
#include "cs_cdovb_advection.h"
#include "cs_cdovb_diffusion.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVB_SCALEQ_DBG  0

/* Algebraic system for CDO vertex-based discretization */

struct _cs_cdovb_scaleq_t {

  /* Pointer to shared structures (i.e. not owned by this structure):
     - Pointer to a cs_equation_param_t structure shared with a cs_equation_t
       structure.
     - Pointer to a cs_cdo_quantities_t  structure
     - Pointer to a cs_cdo_connect_t  structure
     - Pointer to a cs_time_step_t struture
  */

  const cs_equation_param_t  *eqp;
  const cs_cdo_quantities_t  *quant;
  const cs_cdo_connect_t     *connect;
  const cs_time_step_t       *time_step;

  /* System size (known boundary entities may be removed if BCs are strongly
     enforced) */
  cs_lnum_t  n_vertices;
  cs_lnum_t  n_dof_vertices; /* n_rows = n_cols = n_vertices - dir. vertices */

  /* Indirection between zipped numbering (without BC) and initial numbering
     Allocated only if the boundary conditions are strongly enforced.
  */
  cs_lnum_t          *v_z2i_ids;  // Mapping n_dof_vertices -> n_vertices
  cs_lnum_t          *v_i2z_ids;  // Mapping n_vertices     -> n_dof_vertices

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

  /* Source terms */
  cs_real_t          *source_terms;

  /* Hodge^{VpCd,Conf} : only if reaction with the same algo for the discrete
     Hodge is used (in all cases, the matrix index is shared) */
  bool                build_hvpcd_conf;
  cs_sla_matrix_t    *hvpcd_conf;

  /* Work buffer */
  cs_lnum_t          *vtag;       /* size: n_vertices, store the local vertex id
                                     or -1 if not activated */

};

/*============================================================================
 * Private variables
 *============================================================================*/

// Advanced developper parameters (weakly enforced the boundary conditions)
static const cs_real_t  cs_weak_nitsche_pena_coef = 500;
static const cs_real_t  cs_weak_penalization_weight = 0.01;

static size_t  _vbscal_work_size = 0;
static cs_real_t  *_vbscal_work = NULL;

/* Index of a CDO vertex-based matrix (avoid to build it at each iteration and
   for each equation).
   v2v connectivity through cell neighboorhood */
static cs_connect_index_t  *cs_cdovb_v2v = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contribution of source terms to the rhs for this
 *          time step
 *
 * \param[in]      builder      pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] full_rhs     right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_add_source_terms(cs_cdovb_scaleq_t          *builder,
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
    cs_cdovb_scaleq_compute_source(builder);

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
 * \brief   Compute a discrete Hodge op. Vp-->Cd using conforming reco. op.
 *
 * \param[in, out] builder     pointer to a cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_hvpcd_conf(cs_cdovb_scaleq_t    *builder)
{
  const cs_cdo_connect_t  *connect = builder->connect;
  const cs_cdo_quantities_t  *quant = builder->quant;

  cs_param_hodge_t  h_info = {.inv_pty = false,
                              .type = CS_PARAM_HODGE_TYPE_VPCD,
                              .algo = CS_PARAM_HODGE_ALGO_WBS,
                              .coef = 1}; // not useful in this context
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect, h_info);

  builder->build_hvpcd_conf = true;

  /* Initialize matrix structure */
  builder->hvpcd_conf =
    cs_sla_matrix_create_msr_from_index(cs_cdovb_v2v, true, true, 1);

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_locmat_t  *hloc = cs_hodge_build_local(c_id, connect, quant, hb);

    cs_sla_assemble_msr_sym(hloc, builder->hvpcd_conf, false);

  }

  /* Free memory */
  hb = cs_hodge_builder_free(hb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the matrix related to the unsteady term
 *
 * \param[in]      builder    pointer to a cs_cdovb_scaleq_t structure
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_time_matrix(cs_cdovb_scaleq_t          *builder)
{
  cs_sla_matrix_t  *time_mat = NULL;

  const cs_equation_param_t  *eqp = builder->eqp;
  const bool  is_unsteady = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;

  if (!is_unsteady) // Steady-state eq. => Nothing to do
    return time_mat;

  const cs_param_hodge_t  h_info = eqp->time_hodge;

  /* Sanity check */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_VPCD);

  if (h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI) {
    time_mat = cs_sla_matrix_create(builder->n_vertices, // n_rows
                                    builder->n_vertices, // n_cols
                                    1,                   // stride
                                    CS_SLA_MAT_MSR,      // type
                                    true);               // symmetric ?

    /* Set matrix flag */
    time_mat->flag |= CS_SLA_MATRIX_SYM;
    time_mat->flag |= CS_SLA_MATRIX_SORTED;

  }
  else
    time_mat =
      cs_sla_matrix_create_msr_from_index(cs_cdovb_v2v, true, true, 1);

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
 * \param[in, out] builder    pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] rhs        pointer to the right-hand side array
 * \param[in, out] sys_mat    pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_apply_time_scheme(const cs_real_t          *field_val,
                   const cs_sla_matrix_t    *time_mat,
                   double                    dt_cur,
                   cs_cdovb_scaleq_t        *builder,
                   cs_real_t                *rhs,
                   cs_sla_matrix_t          *sys_mat)

{
  cs_lnum_t  i;

  cs_real_t  *mv_time = _vbscal_work;
  cs_real_t  *mv_sys = _vbscal_work + builder->n_vertices;
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
 * \param[in]      mesh         pointer to a cs_mesh_t structure
 * \param[in]      field_val    pointer to the current value of the field
 * \param[in, out] builder      pointer to a cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_dir_values(const cs_mesh_t            *mesh,
                    const cs_real_t            *field_val,
                    const cs_cdovb_scaleq_t    *builder)
{
  cs_lnum_t  i;

  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;
  const cs_equation_param_t  *eqp = builder->eqp;
  const cs_time_step_t  *time_step = builder->time_step;

  if (vtx_dir->n_nhmg_elts == 0)
    return; // Nothing to do

  cs_flag_t  dof_flag = CS_FLAG_VERTEX | CS_FLAG_PRIMAL | CS_FLAG_SCAL;

  /* Get the value of the Dirichlet for the current time */
  cs_cdo_bc_dirichlet_set(dof_flag,
                          time_step,
                          mesh,
                          eqp->bc,
                          vtx_dir,
                          builder->dir_val);

  if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
      builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
    if (eqp->flag & CS_EQUATION_UNSTEADY) {

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
 * \param[in, out] adv          pointer to an advection builder structure
 * \param[in, out] rhs          pointer of pointer to the right-hand side
 * \param[in, out] matrix       pointer to a matrix structure
 */
/*----------------------------------------------------------------------------*/

static void
_add_advection_bc(cs_cdovb_scaleq_t           *builder,
                  cs_cdovb_adv_t              *adv,
                  cs_real_t                   *rhs,
                  cs_sla_matrix_t             *matrix)
{
  cs_lnum_t  i;

  /* temporary buffer */
  cs_real_t  *dir_vals = _vbscal_work;
  cs_real_t  *rhs_contrib = _vbscal_work + builder->n_vertices;
  cs_real_t  *diag_contrib = _vbscal_work + 2*builder->n_vertices;

  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;
  const cs_cdo_connect_t  *connect = builder->connect;
  const cs_cdo_quantities_t  *quant = builder->quant;

  /* Initialize arrays */
  for (i = 0; i < builder->n_vertices; i++)
    dir_vals[i] = rhs_contrib[i] = diag_contrib[i] = 0;

  /* Store the value of the Dirichlet condition into a full (i.e. n_vertices)
     array */
  for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
    dir_vals[vtx_dir->elt_ids[i]] = builder->dir_val[i];

  /* Compute the contribution */
  cs_cdovb_advection_add_bc(connect, quant,
                            dir_vals,
                            adv,
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
 * \param[in, out] builder      pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] diff_builder pointer to a cs_cdovb_diff_t structure
 * \param[in, out] full_rhs     right-hand side
 * \param[in, out] full_matrix  matrix of the linear system
 */
/*----------------------------------------------------------------------------*/

static void
_weak_bc_enforcement(cs_cdovb_scaleq_t     *builder,
                     cs_cdovb_diff_t       *diff_builder,
                     cs_real_t              full_rhs[],
                     cs_sla_matrix_t       *full_matrix)
{
  if (builder->enforce !=  CS_PARAM_BC_ENFORCE_WEAK_SYM &&
      builder->enforce !=  CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
    return; // Nothing to do

  /* Sanity check */
  assert(full_matrix->type == CS_SLA_MAT_MSR);

  const cs_cdo_connect_t  *connect = builder->connect;
  const cs_cdo_quantities_t  *quant = builder->quant;
  const cs_equation_param_t  *eqp = builder->eqp;
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;
  const cs_property_t  *pty = eqp->diffusion_property;

  cs_lnum_t  i, j, v_id;
  cs_real_t  eig_ratio, eig_max, pena_coef;
  cs_real_t  matpty[3][3];

  cs_real_t  *_vec = NULL, *_matvec = NULL;
  cs_locmat_t  *transp = NULL, *ntrgrd = NULL;

  /* Temporary buffers stored using builder */
  cs_lnum_t  *loc_v_id = builder->vtag;
  cs_real_t  *dir_vals = _vbscal_work;
  cs_real_t  *v_coef = _vbscal_work + builder->n_vertices;

  /* Initialize v_coef and loc_v_ids */
  for (i = 0; i < builder->n_vertices; i++) {
    v_coef[i] = 0.0;
    loc_v_id[i] = -1;
  }

  if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

    for (i = 0; i < builder->n_vertices; i++)
      dir_vals[i] = 0.;
    for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
      dir_vals[vtx_dir->elt_ids[i]] = builder->dir_val[i];

  }

  /* Parameters linked to the discrete Hodge operator used for diffusion */
  bool  is_uniform = cs_property_is_uniform(pty);

  /* Get the anisotropic ratio and the max. eigenvalue (if uniform) */
  if (is_uniform) {
    cs_property_get_cell_tensor(0, pty, h_info.inv_pty, matpty);
    cs_math_33_eigen((const cs_real_t (*)[3])matpty, &eig_ratio, &eig_max);
  }

  const cs_cdo_bc_list_t  *face_dir = builder->face_bc->dir;

  /* Loop on Dirichlet faces */
  for (i = 0; i < face_dir->n_elts; i++) {

    cs_lnum_t  f_id = face_dir->elt_ids[i] + quant->n_i_faces;
    cs_lnum_t  c_id = connect->f2c->col_id[connect->f2c->idx[f_id]];

    /* Sanity check (this is a border face) */
    assert(connect->f2c->idx[f_id+1] - connect->f2c->idx[f_id] == 1);

    if (!is_uniform) {
      cs_property_get_cell_tensor(c_id, pty, h_info.inv_pty, matpty);
      cs_math_33_eigen((const cs_real_t (*)[3])matpty, &eig_ratio, &eig_max);
    }

    ntrgrd = cs_cdovb_diffusion_ntrgrd_build(c_id, f_id,
                                             connect, quant,
                                             (const cs_real_t (*)[3])matpty,
                                             eig_ratio, eig_max,
                                             loc_v_id, v_coef,
                                             diff_builder);

    /* Assemble contributions coming from local normal trace operator */
    if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

      /* Additional auxiliary buffers */
      BFT_MALLOC(_vec, connect->n_max_vbyc, cs_real_t);
      BFT_MALLOC(_matvec, connect->n_max_vbyc, cs_real_t);
      transp = cs_locmat_create(connect->n_max_vbyc);

      /* ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) */
      cs_locmat_add_transpose(ntrgrd, transp);

#if CS_CDOVB_SCALEQ_DBG > 1
      bft_printf(">> Local weak sym bc matrix (f_id: %d)", f_id);
      cs_locmat_dump(c_id, ntrgrd);
#endif

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

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply boundary conditions. Update right-hand side and the system
 *          matrix
 *
 * \param[in, out] builder      pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] full_rhs     right-hand side
 * \param[in, out] full_matrix  matrix of the linear system
 */
/*----------------------------------------------------------------------------*/

static void
_strong_bc_enforcement(cs_cdovb_scaleq_t       *builder,
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
  double  *tmp_rhs = _vbscal_work;
  double  *x_bc = _vbscal_work + builder->n_vertices;
  double  *contrib = _vbscal_work + 2*builder->n_vertices;

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
 * \param[in, out] builder      pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] rhs          right-hand side
 * \param[in, out] matrix       matrix of the linear system
 */
/*----------------------------------------------------------------------------*/

static void
_enforce_bc(cs_cdovb_scaleq_t          *builder,
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
      cs_real_t  pena_coef =
        cs_weak_penalization_weight / cs_math_get_machine_epsilon();

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
 * \brief  Allocate work buffer and general structures related to CDO
 *         vertex-based schemes
 *
 * \param[in] connect   pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_initialize(const cs_cdo_connect_t      *connect)
{
  /* Sanity check */
  assert(_vbscal_work == NULL && _vbscal_work_size == 0);

  const cs_lnum_t  n_vertices = connect->v_info->n_ent;
  const cs_lnum_t  n_cells = connect->c_info->n_ent;

  /* Work buffers */
  _vbscal_work_size = CS_MAX(3*n_vertices, n_cells);
  BFT_MALLOC(_vbscal_work, _vbscal_work_size, cs_real_t);

  /* Build a (sorted) v2v connectivity index */
  const cs_connect_index_t  *c2v = connect->c2v;
  cs_connect_index_t  *v2c = cs_index_transpose(n_vertices, c2v);

  cs_cdovb_v2v = cs_index_compose(n_vertices, v2c, c2v);
  cs_index_sort(cs_cdovb_v2v);

  /* Update index (v2v has a diagonal entry. We remove it since we consider a
     matrix stored using the MSR format */
  cs_lnum_t  shift = 0;
  cs_lnum_t  *tmp_idx = NULL;

  BFT_MALLOC(tmp_idx, n_vertices + 1, cs_lnum_t);
  memcpy(tmp_idx, cs_cdovb_v2v->idx, sizeof(cs_lnum_t)*(n_vertices+1));

  for (cs_lnum_t  i = 0; i < n_vertices; i++) {

    cs_lnum_t  start = tmp_idx[i], end = tmp_idx[i+1];

    for (cs_lnum_t  j = start; j < end; j++)
      if (cs_cdovb_v2v->ids[j] != i)
        cs_cdovb_v2v->ids[shift++] = cs_cdovb_v2v->ids[j];

    cs_cdovb_v2v->idx[i+1] = cs_cdovb_v2v->idx[i] + end-start-1;

  } // Loop on vertices

  assert(shift == cs_cdovb_v2v->idx[n_vertices]); // Sanity check

  /* Free temporary buffers */
  cs_index_free(&v2c);
  BFT_FREE(tmp_idx);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO vertex-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_finalize(void)
{
  if (_vbscal_work != NULL) {
    _vbscal_work_size = 0;
    BFT_FREE(_vbscal_work );
  }

  cs_index_free(&(cs_cdovb_v2v));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to a temporary buffer related to scalar equations
 *         discretized with CDO vertex-based schemes
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_scaleq_get_tmpbuf(void)
{
  /* Sanity check */
  assert(_vbscal_work != NULL);

  return _vbscal_work;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_scaleq_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 * \param[in] connect   pointer to a cs_cdo_connect_t structure
 * \param[in] quant      pointer to a cs_cdo_quantities_t structure
 * \param[in] time_step  time_step structure
 *
 * \return a pointer to a new allocated cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_scaleq_init(const cs_equation_param_t   *eqp,
                     const cs_mesh_t             *mesh,
                     const cs_cdo_connect_t      *connect,
                     const cs_cdo_quantities_t   *quant,
                     const cs_time_step_t        *time_step)
{
  cs_lnum_t  i;

  /* Sanity checks */
  assert(eqp != NULL);
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(eqp->var_type == CS_PARAM_VAR_SCAL);

  const cs_lnum_t  n_vertices = connect->v_info->n_ent;
  const cs_lnum_t  n_b_faces = connect->f_info->n_ent_bd;

  cs_cdovb_scaleq_t  *builder = NULL;

  BFT_MALLOC(builder, 1, cs_cdovb_scaleq_t);

  /* Shared pointers */
  builder->eqp = eqp;
  builder->connect = connect;
  builder->quant = quant;
  builder->time_step = time_step;

  /* Dimensions:
     By default, we set number of DoFs as if there is a weak enforcement of
     the boundary conditions */
  builder->n_vertices = n_vertices;
  builder->n_dof_vertices = n_vertices;

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

  /* Members of the structure related to source terms */
  BFT_MALLOC(builder->source_terms, builder->n_vertices, cs_real_t);
  builder->build_hvpcd_conf = false;
  builder->hvpcd_conf = NULL;

  /* Initialize tags */
  BFT_MALLOC(builder->vtag, n_vertices, cs_lnum_t);
  for (i = 0; i < n_vertices; i++)
    builder->vtag[i] = -1;

  return builder;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovb_scaleq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_scaleq_free(void   *builder)
{
  cs_cdovb_scaleq_t  *_builder = (cs_cdovb_scaleq_t *)builder;

  if (_builder == NULL)
    return _builder;

  /* eqp, connect, quant, and time_step are only shared.
     These quantities are freed later. */

  /* Free BC structure */
  if (_builder->vtx_dir->n_nhmg_elts > 0)
    BFT_FREE(_builder->dir_val);

  _builder->face_bc = cs_cdo_bc_free(_builder->face_bc);
  _builder->vtx_dir = cs_cdo_bc_list_free(_builder->vtx_dir);

  /* Free Hodge operator defined from conforming reconstruction op. */
  _builder->hvpcd_conf = cs_sla_matrix_free(_builder->hvpcd_conf);

  /* Renumbering (if strong enforcement of BCs for instance) */
  if (_builder->n_vertices > _builder->n_dof_vertices) {
    BFT_FREE(_builder->v_z2i_ids);
    BFT_FREE(_builder->v_i2z_ids);
  }

  BFT_FREE(_builder->vtag);
  BFT_FREE(_builder->source_terms);
  BFT_FREE(_builder);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in, out] builder     pointer to a cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_compute_source(void            *builder)
{
  cs_cdovb_scaleq_t  *bld = (cs_cdovb_scaleq_t *)builder;
  double  *st_eval = _vbscal_work;
  
  for (cs_lnum_t i = 0; i < bld->n_vertices; i++)
    bld->source_terms[i] = 0;

  const cs_equation_param_t  *eqp = bld->eqp;

  if (eqp->n_source_terms == 0)
    return;

  cs_desc_t  desc;

  if (eqp->flag & CS_EQUATION_HCONF_ST) {
    desc.location = CS_FLAG_SCAL | cs_cdo_primal_vtx;
    desc.state = CS_FLAG_STATE_POTENTIAL;
    
    if (!bld->build_hvpcd_conf)
      _build_hvpcd_conf(bld);

  }
  else {
    desc.location = CS_FLAG_SCAL | cs_cdo_dual_cell;
    desc.state = CS_FLAG_STATE_DENSITY;
  }
  
  for (int  st_id = 0; st_id < eqp->n_source_terms; st_id++) {

    const cs_source_term_t  *st = eqp->source_terms[st_id];

    cs_source_term_compute(desc, st, &st_eval); // updated inside this function

    /* Update source term array */
    if (eqp->flag & CS_EQUATION_HCONF_ST) {

      double  *mv = _vbscal_work + bld->n_vertices;

      cs_sla_matvec(bld->hvpcd_conf, st_eval, &mv, true);
      for (cs_lnum_t i = 0; i < bld->n_vertices; i++)
        bld->source_terms[i] += mv[i];

    }
    else
      for (cs_lnum_t i = 0; i < bld->n_vertices; i++)
        bld->source_terms[i] += st_eval[i];

  } // Loop on source terms

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO vertex-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to cs_cdovb_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] sla_mat    pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_build_system(const cs_mesh_t             *mesh,
                             const cs_real_t             *field_val,
                             double                       dt_cur,
                             void                        *builder,
                             cs_real_t                  **rhs,
                             cs_sla_matrix_t            **sla_mat)
{
  cs_lnum_t  i, j;
  cs_real_t  ptyval;

  bool  diff_pty_uniform = true, time_pty_uniform = true;
  cs_real_33_t  diff_tensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  cs_cdovb_scaleq_t  *sys_builder = (cs_cdovb_scaleq_t *)builder;

  const cs_equation_param_t  *eqp = sys_builder->eqp;
  const cs_cdo_connect_t  *connect = sys_builder->connect;
  const cs_cdo_quantities_t  *quant = sys_builder->quant;
  const cs_connect_index_t  *c2v = connect->c2v;
  const bool  do_diffusion = (eqp->flag & CS_EQUATION_DIFFUSION) ? true:false;
  const bool  do_advection = (eqp->flag & CS_EQUATION_CONVECTION) ? true:false;
  const bool  do_unsteady = (eqp->flag & CS_EQUATION_UNSTEADY) ? true:false;
  const bool  do_reaction = (eqp->flag & CS_EQUATION_REACTION) ? true:false;

  /* Structures used to build the linear system */
  cs_cdovb_diff_t  *diff = NULL;
  cs_cdovb_adv_t  *adv = NULL;
  cs_hodge_builder_t  **reac_builder = NULL;
  cs_sla_matrix_t  *time_mat = NULL;
  cs_hodge_builder_t  *time_builder = NULL;
  cs_locmat_t  *adr_loc = cs_locmat_create(connect->n_max_vbyc);

  /* Allocate and initialize a matrix with the larger stencil (that related
     to diffusion => all vertices of a cell are potentially in interaction)
     adr = advection/diffusion/reaction
  */
  cs_sla_matrix_t  *sys_mat =
    cs_sla_matrix_create_msr_from_index(cs_cdovb_v2v,
                                        false,  // symmetric
                                        true,   // sorted
                                        1);     // stride

  if (!do_advection && sys_builder->enforce != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
    sys_mat->flag |= CS_SLA_MATRIX_SYM;

  /* Preparatory step for diffusion term */
  if (do_diffusion) {
    diff_pty_uniform = cs_property_is_uniform(eqp->diffusion_property);
    diff = cs_cdovb_diffusion_builder_init(connect,
                                           diff_pty_uniform,
                                           eqp->diffusion_hodge,
                                           sys_builder->enforce);
  }

  /* Preparatory step for advection term */
  if (do_advection)
    adv = cs_cdovb_advection_builder_init(connect,
                                          eqp->advection_field,
                                          eqp->advection_info,
                                          do_diffusion);

  /* Preparatory step for reaction term */
  if (do_reaction) {
    BFT_MALLOC(reac_builder, eqp->n_reaction_terms, cs_hodge_builder_t *);
    for (i = 0; i < eqp->n_reaction_terms; i++)
      reac_builder[i] = cs_hodge_builder_init(connect,
                                              eqp->reaction_terms[i].hodge);
  }

  /* Preparatory step for unsteady term */
  if (do_unsteady) {
    time_mat = _init_time_matrix(sys_builder);
    time_builder = cs_hodge_builder_init(connect, eqp->time_hodge);
    time_pty_uniform = cs_property_is_uniform(eqp->diffusion_property);
  }

  const bool  only_time_diag =
    (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_VORONOI) ? true : false;

  /* Main loop on cells to build the system matrix */
  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Initialize local system matrix related to advection/diffusion/reaction */
    adr_loc->n_ent = c2v->idx[c_id+1] - c2v->idx[c_id];
    for (i = c2v->idx[c_id], j = 0; i < c2v->idx[c_id+1]; i++, j++) {
      adr_loc->ids[j] = c2v->ids[i];
      sys_builder->vtag[adr_loc->ids[j]] = j;
    }

    for (i = 0; i < adr_loc->n_ent*adr_loc->n_ent; i++)
      adr_loc->mat[i] = 0;

    /* DIFFUSION TERM */
    if (do_diffusion) { /* Define the local matrix related to diffusion
                           ==> stiffness matrix */

      if (c_id == 0 || diff_pty_uniform == false)
        cs_property_get_cell_tensor(c_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    diff_tensor);

      cs_locmat_t  *diff_loc =
        cs_cdovb_diffusion_build_local(c_id, connect, quant,
                                       sys_builder->vtag,
                                       (const cs_real_3_t (*))diff_tensor,
                                       diff);

      cs_locmat_add(adr_loc, diff_loc);

#if CS_CDOVB_SCALEQ_DBG > 1
      bft_printf(">> Local diffusion matrix");
      cs_locmat_dump(c_id, diff_loc);
#endif
    }

    /* ADVECTION TERM */
    if (do_advection) { /* Define the local matrix related to advection */

      cs_locmat_t  *adv_loc =
        cs_cdovb_advection_build_local(c_id, connect, quant,
                                       sys_builder->vtag,
                                       (const cs_real_3_t (*))diff_tensor,
                                       adv);

#if CS_CDOVB_SCALEQ_DBG > 1
      bft_printf(">> Local advection matrix");
      cs_locmat_dump(c_id, adv_loc);
#endif

      cs_locmat_add(adr_loc, adv_loc);
    }

    /* REACTION TERM */
    if (do_reaction) { /* Define the local matrix related to reaction
                          ==> mass matrix */

      for (i = 0; i < eqp->n_reaction_terms; i++) {

        cs_hodge_builder_t  *hb = reac_builder[i];
        cs_locmat_t  *rea_loc = cs_hodge_build_local(c_id, connect, quant, hb);

#if CS_CDOVB_SCALEQ_DBG > 1
        bft_printf(">> Local reaction matrix for term %d", i);
        cs_locmat_dump(c_id, rea_loc);
#endif
        cs_locmat_add(adr_loc, rea_loc);

      }

    }

    /* TIME TERM */
    if (do_unsteady) { /* Build mass matrix to take into account time effect */

      if (c_id == 0 || time_pty_uniform == false) {
        ptyval = cs_property_get_cell_value(c_id, eqp->time_property);
        cs_hodge_builder_set_val(time_builder, ptyval);
      }

      cs_locmat_t  *time_loc = cs_hodge_build_local(c_id, connect, quant,
                                                    time_builder);

#if CS_CDOVB_SCALEQ_DBG > 1
      bft_printf(">> Local time matrix");
      cs_locmat_dump(c_id, time_loc);
#endif

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

#if CS_CDOVB_SCALEQ_DBG > 0
    bft_printf(">> Local system matrix");
    cs_locmat_dump(c_id, adr_loc);
#endif

  } // Main loop on cells

  /* Free temporary buffers and structures */
  adr_loc = cs_locmat_free(adr_loc);

  if (do_reaction) {
    for (i = 0; i < eqp->n_reaction_terms; i++)
      reac_builder[i] = cs_hodge_builder_free(reac_builder[i]);
    BFT_FREE(reac_builder);
  }

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
  _add_source_terms(sys_builder, full_rhs);

  /* Compute the values of the Dirichlet BC.
     TODO: do the analogy for Neumann BC */
  _compute_dir_values(mesh, field_val, sys_builder);

  if (do_diffusion) { /* Last treatment for the diffusion term */

    /* Weakly enforce Dirichlet boundary conditions */
    _weak_bc_enforcement(sys_builder, diff, full_rhs, sys_mat);

    /* Free associated builder structure */
    diff = cs_cdovb_diffusion_builder_free(diff);

  }

  if (do_advection) { /* Last treatment for the advection term */

    /* Apply boundary conditions */
    _add_advection_bc(sys_builder, adv, full_rhs, sys_mat);

    /* Free associated builder structure */
    adv = cs_cdovb_advection_builder_free(adv);

  }

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
     Apply the strong or penalized enforcement. In case of Nitsche enforcement,
     there is nothing to do (already done).
     Must be call after the application of the time scheme */
  _enforce_bc(sys_builder, &full_rhs, &sys_mat);

  bool do_cleaning = false; // Advanced option
  if (do_cleaning)
    cs_sla_matrix_clean(sys_mat, cs_math_get_machine_epsilon());
  else
    cs_sla_matrix_rmzeros(sys_mat);

#if CS_CDOVB_SCALEQ_DBG > 1
  cs_sla_system_dump("system.log", NULL, sys_mat, full_rhs);
#endif

  /* Return pointers */
  *rhs = full_rhs;
  *sla_mat = sys_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO vertex-based scheme.
 *
 * \param[in]      solu       solution array
 * \param[in, out] builder    pointer to cs_cdovb_scaleq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_update_field(const cs_real_t     *solu,
                             void                *builder,
                             cs_real_t           *field_val)
{
  int  i;

  cs_cdovb_scaleq_t  *b = (cs_cdovb_scaleq_t  *)builder;

  /* Set computed solution in field array */
  if (b->n_dof_vertices < b->n_vertices) {

    const cs_cdo_bc_list_t  *v_dir = b->vtx_dir;

    /* Sanity check */
    assert(b->enforce == CS_PARAM_BC_ENFORCE_STRONG);

    for (i = 0; i < b->n_vertices; i++)
      field_val[i] = 0; // To tag unset values
    for (i = 0; i < b->n_dof_vertices; i++)
      field_val[b->v_z2i_ids[i]] = solu[i];
    for (i = 0; i < v_dir->n_nhmg_elts; i++)
      field_val[v_dir->elt_ids[i]] = b->dir_val[i];

  }
  else
    for (i = 0; i < b->n_vertices; i++)
      field_val[i] = solu[i];

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field strufcture
 * \param[in, out]  builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_extra_op(const char            *eqname,
                         const cs_field_t      *field,
                         void                  *builder)
{
  int  len;

  char *postlabel = NULL;
  cs_cdovb_scaleq_t  *b = (cs_cdovb_scaleq_t  *)builder;

  const cs_equation_param_t  *eqp = b->eqp;

  bool  do_adv = eqp->flag & CS_EQUATION_CONVECTION;
  bool  do_diff = eqp->flag & CS_EQUATION_DIFFUSION;
  bool  do_peclet_post = (eqp->process_flag & CS_EQUATION_POST_PECLET ||
                          eqp->process_flag & CS_EQUATION_POST_UPWIND_COEF);

  if (do_adv && do_diff && do_peclet_post) {

    cs_real_t  *work_c = _vbscal_work;
    cs_real_3_t  base_vect;

    len = strlen(eqname) + 9 + 1;
    BFT_MALLOC(postlabel, len, char);

    for (int k = 0; k < 3; k++) {

      if (k == 0) {
        sprintf(postlabel, "%s.Peclet.X", eqname);
        base_vect[0] = 1, base_vect[1] = base_vect[2] = 0;
      }
      else if (k == 1) {
        sprintf(postlabel, "%s.Peclet.Y", eqname);
        base_vect[1] = 1, base_vect[0] = base_vect[2] = 0;
      }
      else {
        sprintf(postlabel, "%s.Peclet.Z", eqname);
        base_vect[2] = 1, base_vect[1] = base_vect[0] = 0;
      }

      cs_cdovb_advection_get_peclet_cell(b->quant,
                                         eqp->advection_field,
                                         eqp->diffusion_property,
                                         base_vect,
                                         &work_c);

      if (eqp->process_flag & CS_EQUATION_POST_PECLET)
        cs_post_write_var(-1,             // id du maillage de post
                          postlabel,
                          1,
                          true,           // interlace
                          true,           // true = original mesh
                          CS_POST_TYPE_cs_real_t,
                          work_c,         // values on cells
                          NULL,           // values at internal faces
                          NULL,           // values at border faces
                          b->time_step);  // time step management structure

      if (eqp->process_flag & CS_EQUATION_POST_UPWIND_COEF) {

        if (k == 0)
          sprintf(postlabel, "%s.UpwCoefX", eqname);
        else if (k == 1)
          sprintf(postlabel, "%s.UpwCoefY", eqname);
        else
          sprintf(postlabel, "%s.UpwCoefZ", eqname);

        /* Compute in each cell an evaluation of upwind weight value */
        cs_cdovb_advection_get_upwind_coef_cell(b->quant,
                                                eqp->advection_info,
                                                work_c);

        cs_post_write_var(-1,             // id du maillage de post
                          postlabel,
                          1,
                          true,           // interlace
                          true,           // true = original mesh
                          CS_POST_TYPE_cs_real_t,
                          work_c,         // values on cells
                          NULL,           // values at internal faces
                          NULL,           // values at border faces
                          b->time_step);  // time step management structure

      } /* Post upwinding coefficient */

    } // Loop on space dimension

    BFT_FREE(postlabel);

  } // Post a Peclet attached to cells

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
