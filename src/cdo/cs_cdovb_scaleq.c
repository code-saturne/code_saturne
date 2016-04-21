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

#include "cs_cdo_bc.h"
#include "cs_cdovb_advection.h"
#include "cs_cdovb_diffusion.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_search.h"
#include "cs_source_term.h"

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

#define CDO_DIFFUSION  0
#define CDO_ADVECTION  1
#define CDO_REACTION   2
#define CDO_TIME       3
#define N_CDO_TERMS    4

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/* Algebraic system for CDO vertex-based discretization */

struct _cs_cdovb_scaleq_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure.  */

  const cs_equation_param_t  *eqp;

  bool       todo[N_CDO_TERMS];

  /* System size (known boundary entities may be removed if BCs are strongly
     enforced) */
  cs_lnum_t  n_vertices;
  cs_lnum_t  n_dof_vertices; /* n_rows = n_cols = n_vertices - dir. vertices */

  /* Indirection between zipped numbering (without BC) and initial numbering
     Allocated only if the boundary conditions are strongly enforced.
  */
  cs_lnum_t     *v_z2i_ids;  // Mapping n_dof_vertices -> n_vertices
  cs_lnum_t     *v_i2z_ids;  // Mapping n_vertices     -> n_dof_vertices

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
  double                *dir_val; // size = vtx_dir->n_nhmg_elts

  /* In case of a weak enforcement of the Dirichlet BCs */
  cs_lnum_t             *c2bf_bc_idx;  // size: n_cells + 1
  cs_lnum_t             *c2bf_bc_ids;  // cell --> border faces ids

  /* Source terms */
  cs_real_t             *source_terms;

  /* Hodge^{VpCd,Conf} : only if reaction with the same algo for the discrete
     Hodge is used (in all cases, the matrix index is shared) */
  bool                   build_hvpcd_conf;
  cs_sla_matrix_t       *hvpcd_conf;

  /* Builder sub-structures */
  double                *loc_rhs; // local contribution to the RHS
  cs_locmat_t           *adr_mat; // local dense matrix
  cs_locmat_t           *tmp_mat; // local dense matrix

  cs_cdovb_diff_t       *diff;     // builder for diffusion
  cs_cdovb_adv_t        *adv;      // builder for advection
  cs_hodge_builder_t   **hb_reac;  // hodge builders for reaction terms
  cs_hodge_builder_t    *hb_time;  // hodge builder for unsteady term

};

/*============================================================================
 * Private variables
 *============================================================================*/

static size_t  cs_cdovb_scal_work_size = 0;
static cs_real_t  *cs_cdovb_scal_work = NULL;
static cs_cdo_locmesh_t  *cs_cell_mesh = NULL;
static cs_cdo_locsys_t  *cs_cell_sys = NULL;
static cs_cdo_locsys_t  *cs_tmp_sys = NULL;

/* Index of a CDO vertex-based matrix (avoid to build it at each iteration and
   for each equation).
   v2v connectivity through cell neighboorhood */
static cs_connect_index_t  *cs_cdovb_v2v = NULL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdovb_quant;
static const cs_cdo_connect_t  *cs_cdovb_connect;
static const cs_time_step_t  *cs_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contribution of source terms to the rhs for this
 *          time step
 *
 * \param[in, out] b          pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] full_rhs   right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_add_source_terms(cs_cdovb_scaleq_t     *b,
                  cs_real_t              full_rhs[])
{
  const cs_equation_param_t  *eqp = b->eqp;

  if (b->todo[CDO_TIME]) {

    const cs_param_time_t  t_info = eqp->time_info;

    /* Previous values are stored inside b->source_terms i.e.
       values of the source terms related to t_prev */
    if (t_info.scheme == CS_TIME_SCHEME_EXPLICIT)
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        full_rhs[i] += b->source_terms[i];

    else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
             t_info.scheme == CS_TIME_SCHEME_THETA) {

      const double  tcoef = 1 - t_info.theta;

# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        full_rhs[i] += tcoef * b->source_terms[i];

    }

    /* Update b->source_term with the value attached to t_cur */
    cs_cdovb_scaleq_compute_source(b);

    if (t_info.scheme == CS_TIME_SCHEME_IMPLICIT)
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        full_rhs[i] += b->source_terms[i];

    else if (t_info.scheme == CS_TIME_SCHEME_CRANKNICO ||
             t_info.scheme == CS_TIME_SCHEME_THETA) {
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < b->n_vertices; i++)
        full_rhs[i] += t_info.theta * b->source_terms[i];

    }

  }
  else { /* Steady case: source terms have already been computed during
            the initialization step */

# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < b->n_vertices; i++)
      full_rhs[i] += b->source_terms[i];

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a discrete Hodge op. Vp-->Cd using conforming reco. op.
 *
 * \param[in, out] b     pointer to a cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_hvpcd_conf(cs_cdovb_scaleq_t    *b)
{
  const cs_cdo_connect_t  *connect = cs_cdovb_connect;
  const cs_cdo_quantities_t  *quant = cs_cdovb_quant;

  cs_param_hodge_t  h_info = {.inv_pty = false,
                              .type = CS_PARAM_HODGE_TYPE_VPCD,
                              .algo = CS_PARAM_HODGE_ALGO_WBS,
                              .coef = 1}; // not useful in this context
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect, h_info);

  b->build_hvpcd_conf = true;

  /* Initialize matrix structure */
  b->hvpcd_conf =
    cs_sla_matrix_create_msr_from_index(cs_cdovb_v2v, true, true, 1);

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_locmat_t  *hloc = cs_hodge_build_local(c_id, connect, quant, hb);

    cs_sla_assemble_msr_sym(hloc, b->hvpcd_conf, false);

  }

  /* Free memory */
  hb = cs_hodge_builder_free(hb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the matrix related to the unsteady term
 *
 * \param[in]      b    pointer to a cs_cdovb_scaleq_t structure
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_time_matrix(cs_cdovb_scaleq_t    *b)
{
  cs_sla_matrix_t  *time_mat = NULL;

  if (!b->todo[CDO_TIME]) // Steady-state eq. => Nothing to do
    return time_mat;

  const cs_equation_param_t  *eqp = b->eqp;
  const cs_param_hodge_t  h_info = eqp->time_hodge;

  /* Sanity check */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_VPCD);

  if (h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI) {
    time_mat = cs_sla_matrix_create(b->n_vertices,  // n_rows
                                    b->n_vertices,  // n_cols
                                    1,              // stride
                                    CS_SLA_MAT_MSR, // type
                                    true);          // symmetric ?

    /* Set matrix flag */
    time_mat->flag |= CS_SLA_MATRIX_SYM;
    time_mat->flag |= CS_SLA_MATRIX_SORTED;

  }
  else
    time_mat = cs_sla_matrix_create_msr_from_index(cs_cdovb_v2v, true, true, 1);

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
  cs_real_t  *mv_time = cs_cdovb_scal_work;
  cs_real_t  *mv_sys = cs_cdovb_scal_work + builder->n_vertices;
  size_t  time_nnz = time_mat->idx[time_mat->n_rows];
  size_t  sys_nnz = sys_mat->idx[sys_mat->n_rows];

  const double inv_dt = 1/dt_cur;
  const cs_equation_param_t  *eqp = builder->eqp;
  const cs_param_time_t  t_info = eqp->time_info;

  /* Sanity checks */
  assert(eqp->flag & CS_EQUATION_UNSTEADY);
  assert(sys_mat->type == CS_SLA_MAT_MSR && time_mat->type == CS_SLA_MAT_MSR);
  assert(sys_mat->n_rows == builder->n_vertices);
  assert(sys_mat->n_rows == time_mat->n_rows);

  switch (t_info.scheme) {

  case CS_TIME_SCHEME_EXPLICIT:

    /* A_{diff,conv,reac} * p^(n) --> RHS */
    cs_sla_matvec(sys_mat, field_val, &mv_sys, true);
    cs_sla_matvec(time_mat, field_val, &mv_time, true);

# pragma omp parallel for if (builder->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < builder->n_vertices; i++)
      rhs[i] += mv_sys[i] + inv_dt * mv_time[i];

    /* Implicit part is only the time_mat */
# pragma omp parallel for if (sys_mat->n_rows > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < sys_mat->n_rows; i++)
      sys_mat->diag[i] = inv_dt * time_mat->diag[i];

    if (sys_nnz != time_nnz) {
      BFT_REALLOC(sys_mat->col_id, time_nnz, cs_lnum_t);
      BFT_REALLOC(sys_mat->val, time_nnz, cs_real_t);
      sys_nnz = time_nnz;
    }
    memcpy(sys_mat->col_id, time_mat->col_id, time_nnz*sizeof(cs_lnum_t));

# pragma omp parallel for if (time_mat->n_rows > CS_THR_MIN)
    for (cs_lnum_t i = 1; i < time_mat->n_rows + 1; i++)
      sys_mat->idx[i] = time_mat->idx[i];
# pragma omp parallel for if (time_nnz > CS_THR_MIN)
    for (size_t ii = 0; ii < time_nnz; ii++)
      sys_mat->val[ii] = inv_dt*time_mat->val[ii];

    break;

  case CS_TIME_SCHEME_IMPLICIT:

    /* Update rhs */
    cs_sla_matvec(time_mat, field_val, &mv_time, true);
# pragma omp parallel for if (builder->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < builder->n_vertices; i++)
      rhs[i] += inv_dt * mv_time[i];

    /* Update the system matrix */
# pragma omp parallel for if (builder->n_vertices > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < builder->n_vertices; i++)
      sys_mat->diag[i] += inv_dt * time_mat->diag[i];

    if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_WBS ||
        eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_COST) {

      /* sys_mat and time_hodge share the same index */
      assert(sys_nnz == time_nnz);
# pragma omp parallel for if (sys_nnz > CS_THR_MIN)
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
# pragma omp parallel for if (builder->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < builder->n_vertices; i++)
        rhs[i] += tcoef*mv_sys[i] + inv_dt*mv_time[i];

      /* Update the diag. terms of the system matrix */
# pragma omp parallel for if (builder->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < builder->n_vertices; i++) {
        sys_mat->diag[i] *= t_info.theta;
        sys_mat->diag[i] += inv_dt * time_mat->diag[i];
      }

      /* Update the extra-diag. terms of the system matrix: sys_mat *= theta */
      if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_WBS ||
          eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_COST) {

        /* sys_mat and time_hodge share the same index */
        assert(sys_nnz == time_nnz);

# pragma omp parallel for if (sys_nnz > CS_THR_MIN)
        for (size_t ii = 0; ii < sys_nnz; ii++) {
          sys_mat->val[ii] *= t_info.theta;
          sys_mat->val[ii] += inv_dt * time_mat->val[ii];
        }

      } // WBS or COST algo

    } // Crank-Nicolson or theta time scheme
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

  if (vtx_dir->n_nhmg_elts == 0)
    return; // Nothing to do

  cs_flag_t  dof_flag = CS_FLAG_VERTEX | CS_FLAG_PRIMAL | CS_FLAG_SCAL;

  /* Get the value of the Dirichlet for the current time */
  cs_cdo_bc_dirichlet_set(dof_flag,
                          cs_time_step,
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
 * \brief   Apply boundary conditions. Update right-hand side and the system
 *          matrix
 *
 * \param[in, out] bld          pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out] full_rhs     right-hand side
 * \param[in, out] full_matrix  matrix of the linear system
 */
/*----------------------------------------------------------------------------*/

static void
_strong_bc_enforcement(cs_cdovb_scaleq_t       *bld,
                       cs_real_t              **rhs,
                       cs_sla_matrix_t        **matrix)
{
  const cs_cdo_bc_list_t  *vtx_dir = bld->vtx_dir;

  if (vtx_dir->n_nhmg_elts == 0)
    return;

  /* Sanity check */
  assert(bld->n_vertices > bld->n_dof_vertices);

  cs_sla_matrix_t  *full_matrix = *matrix, *final_matrix = NULL;
  double  *full_rhs = *rhs, *final_rhs = NULL;
  double  *tmp_rhs = cs_cdovb_scal_work;
  double  *x_bc = cs_cdovb_scal_work + bld->n_vertices;
  double  *contrib = cs_cdovb_scal_work + 2*bld->n_vertices;

# pragma omp parallel for if (bld->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < bld->n_vertices; i++)
    x_bc[i] = 0.0;
# pragma omp parallel for if (vtx_dir->n_nhmg_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < vtx_dir->n_nhmg_elts; i++)
    x_bc[vtx_dir->elt_ids[i]] = bld->dir_val[i];

  /* Compute full_matrix*Tbc: rhs = rhs - full_matrix*Tbc */
  cs_sla_matvec(full_matrix, x_bc, &contrib, true);
# pragma omp parallel for if (bld->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < bld->n_vertices; i++)
    full_rhs[i] -= contrib[i];

  /* Reduce the rhs size. Remove vertices with Dirichlet BC */
  memcpy(tmp_rhs, full_rhs, bld->n_vertices*sizeof(double));
  final_rhs = *rhs;
# pragma omp parallel for if (bld->n_dof_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < bld->n_dof_vertices; i++)
    final_rhs[i] = tmp_rhs[bld->v_z2i_ids[i]];

  /* Reduce the system size from n_vertices to n_dof_vertices.
     Vertices attached to a Dirichlet BC are removed.
     Extract block with degrees of freedom */
  final_matrix = cs_sla_matrix_pack(bld->n_dof_vertices,
                                    bld->n_dof_vertices,
                                    full_matrix,
                                    bld->v_z2i_ids,
                                    bld->v_i2z_ids,
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
 * \param[in, out]  bld       pointer to a cs_cdovb_scaleq_t structure
 * \param[in, out]  rhs       right-hand side
 * \param[in, out]  matrix    matrix of the linear system
 */
/*----------------------------------------------------------------------------*/

static void
_enforce_bc(cs_cdovb_scaleq_t          *bld,
            cs_real_t                 **rhs,
            cs_sla_matrix_t           **matrix)
{
  /* Sanity check */
  if (bld->enforce != CS_PARAM_BC_ENFORCE_STRONG &&
      bld->n_vertices != bld->n_dof_vertices)
    bft_error(__FILE__, __LINE__, 0,
              " Error detected: Boundary conditions are not strongly enforced"
              " but there are some removed vertices.");

  cs_sla_matrix_t  *full_matrix = *matrix;
  double  *full_rhs = *rhs;

  const cs_cdo_bc_list_t  *vtx_dir = bld->vtx_dir;

  /* Treatment differs according to the way of enforcing BCs.
     In vertex-based scheme, Dirichlet BC are essential and Neuman BC natural */
  switch (bld->enforce) {

  case CS_PARAM_BC_ENFORCE_STRONG:
    _strong_bc_enforcement(bld, rhs, matrix);
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_PENA:
    {
      // Advanced parameters
      const cs_real_t  penalization_coef = 1e-2/cs_math_get_machine_epsilon();

# pragma omp parallel for if (vtx_dir->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < vtx_dir->n_elts; i++)
        full_matrix->diag[vtx_dir->elt_ids[i]] += penalization_coef;

# pragma omp parallel for if (vtx_dir->n_nhmg_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < vtx_dir->n_nhmg_elts; i++)
        full_rhs[vtx_dir->elt_ids[i]] += penalization_coef * bld->dir_val[i];

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

  if (bld->n_vertices == bld->n_dof_vertices) { // Keep the full system
    *matrix = full_matrix;
    *rhs = full_rhs;
  }

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                    const cs_cdo_connect_t       *connect,
                                    const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */
  cs_cdovb_quant = quant;
  cs_cdovb_connect = connect;
  cs_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to CDO
 *         vertex-based schemes
 *
 * \param[in] connect   pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_initialize(void)
{
  /* Sanity check */
  assert(cs_cdovb_scal_work == NULL && cs_cdovb_scal_work_size == 0);

  const cs_cdo_connect_t  *connect = cs_cdovb_connect;
  const cs_lnum_t  n_vertices = cs_cdovb_quant->n_vertices;
  const cs_lnum_t  n_cells = cs_cdovb_quant->n_cells;

  /* Work buffers */
  cs_cdovb_scal_work_size = CS_MAX(3*n_vertices, n_cells);
  BFT_MALLOC(cs_cdovb_scal_work, cs_cdovb_scal_work_size, cs_real_t);

  /* Structure to map local mesh connectivities and quantities */
  cs_cell_mesh = cs_cdo_locmesh_create(connect);

  /* Structure used to build the final system by a cell-wise process */
  cs_cell_sys = cs_cdo_locsys_create(connect->n_max_vbyc);
  cs_tmp_sys = cs_cdo_locsys_create(connect->n_max_vbyc);

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
  if (cs_cdovb_scal_work != NULL) {
    cs_cdovb_scal_work_size = 0;
    BFT_FREE(cs_cdovb_scal_work );
  }

  cs_index_free(&cs_cdovb_v2v);

  /* Free local structures */
  cs_cdo_locsys_free(&cs_cell_sys);
  cs_cdo_locsys_free(&cs_tmp_sys);
  cs_cdo_locmesh_free(&cs_cell_mesh);
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
  assert(cs_cdovb_scal_work != NULL);

  return cs_cdovb_scal_work;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_scaleq_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_scaleq_init(const cs_equation_param_t   *eqp,
                     const cs_mesh_t             *mesh)
{
  /* Sanity checks */
  assert(eqp != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVB &&
      eqp->var_type != CS_PARAM_VAR_SCAL)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued CDO vertex-based equation.");

  const cs_cdo_connect_t  *connect = cs_cdovb_connect;
  const cs_lnum_t  n_vertices = connect->v_info->n_elts;
  const cs_lnum_t  n_b_faces = connect->f_info->n_b_elts;
  const cs_lnum_t  n_i_faces = connect->f_info->n_i_elts;
  const cs_lnum_t  n_cells = connect->c_info->n_elts;

  cs_cdovb_scaleq_t  *bld = NULL;

  BFT_MALLOC(bld, 1, cs_cdovb_scaleq_t);

  /* Shared pointers */
  bld->eqp = eqp;

  /* Store a direct access to which term one has to compute */
  bld->todo[CDO_DIFFUSION] = (eqp->flag & CS_EQUATION_DIFFUSION) ? true:false;
  bld->todo[CDO_ADVECTION] = (eqp->flag & CS_EQUATION_CONVECTION) ? true:false;
  bld->todo[CDO_REACTION] = (eqp->flag & CS_EQUATION_REACTION) ? true:false;
  bld->todo[CDO_TIME] = (eqp->flag & CS_EQUATION_UNSTEADY) ? true:false;

  /* Dimensions: By default, we set number of DoFs as if there is a weak
     enforcement of the boundary conditions */
  bld->n_vertices = n_vertices;
  bld->n_dof_vertices = n_vertices;

  /* Set members and structures related to the management of the BCs */

  const cs_param_bc_t  *bc_param = eqp->bc;

  /* Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.
     We also compute also the list of Dirichlet vertices along with their
     related definition.
  */
  bld->face_bc = cs_cdo_bc_init(bc_param, n_b_faces);
  bld->vtx_dir = cs_cdo_bc_vtx_dir_create(mesh, bld->face_bc);

  /* Allocate and initialize dir_val */
  BFT_MALLOC(bld->dir_val, bld->vtx_dir->n_nhmg_elts, double);
# pragma omp parallel for if (bld->vtx_dir->n_nhmg_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < bld->vtx_dir->n_nhmg_elts; i++)
    bld->dir_val[i] = 0.0;

  /* Strong enforcement means that we need an indirection list between the
     compress (or zip) and initial numbering of vertices */
  bld->enforce = bc_param->enforcement;

  bld->v_z2i_ids = NULL; // zipped --> initial ids
  bld->v_i2z_ids = NULL; // initial --> zipped ids
  bld->c2bf_bc_idx = NULL;
  bld->c2bf_bc_ids = NULL;

  switch (bld->enforce) {
  case CS_PARAM_BC_ENFORCE_STRONG:
    if (bld->vtx_dir->n_elts > 0) {

      if (bld->todo[CDO_ADVECTION] || bld->todo[CDO_TIME])
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid choice of enforcement of the boundary conditions.\n"
                  " Strong enforcement is not implemented when convection or"
                  " unsteady terms are activated.\n"
                  " Please modify your settings.");

      bool  *is_kept = NULL;

      bld->n_dof_vertices = n_vertices - bld->vtx_dir->n_elts;

      /* Build bld->v_z2i_ids and bld->i2i_ids */
      BFT_MALLOC(bld->v_z2i_ids, bld->n_dof_vertices, cs_lnum_t);
      BFT_MALLOC(bld->v_i2z_ids, bld->n_vertices, cs_lnum_t);
      BFT_MALLOC(is_kept, n_vertices, bool);

# pragma omp parallel for if (n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_vertices; i++)
        is_kept[i] = true, bld->v_i2z_ids[i] = -1; // by default, set to remove

# pragma omp parallel for if (bld->vtx_dir->n_elts > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < bld->vtx_dir->n_elts; i++)
        is_kept[bld->vtx_dir->elt_ids[i]] = false;

      cs_lnum_t  cur_id = 0;
      for (cs_lnum_t i = 0; i < bld->n_vertices; i++) {
        if (is_kept[i]) {
          bld->v_i2z_ids[i] = cur_id;
          bld->v_z2i_ids[cur_id++] = i;
        }
      }
      assert(cur_id == bld->n_dof_vertices);

      BFT_FREE(is_kept);

    }
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    if (bld->todo[CDO_DIFFUSION]) {

      short int  *count = NULL;

      const cs_sla_matrix_t  *f2c = connect->f2c;
      const cs_cdo_bc_list_t  *face_dir = bld->face_bc->dir;

      /* Allocation and initialization */
      BFT_MALLOC(count, n_cells, short int);
      BFT_MALLOC(bld->c2bf_bc_idx, n_cells + 1, cs_lnum_t);
      bld->c2bf_bc_idx[n_cells] = 0;
# pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++)
        bld->c2bf_bc_idx[i] = 0, count[i] = 0;

      /* First pass: Loop on Dirichlet faces to build index */
      for (cs_lnum_t i = 0; i < face_dir->n_elts; i++) {

        cs_lnum_t  f_id = face_dir->elt_ids[i] + n_i_faces;
        cs_lnum_t  c_id = f2c->col_id[f2c->idx[f_id]];

        assert(f2c->idx[f_id+1] - f2c->idx[f_id] == 1); // check if border
        count[c_id] += 1;

      }

# pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_cells; i++) {
        bld->c2bf_bc_idx[i+1] = bld->c2bf_bc_idx[i] + count[i];
        count[i] = 0;
      }

      /* Second pass: Loop on Dirichlet faces to build list of ids */
      BFT_MALLOC(bld->c2bf_bc_ids, bld->c2bf_bc_idx[n_cells], cs_lnum_t);
      for (cs_lnum_t i = 0; i < face_dir->n_elts; i++) {

        cs_lnum_t  f_id = face_dir->elt_ids[i] + n_i_faces;
        cs_lnum_t  c_id = connect->f2c->col_id[connect->f2c->idx[f_id]];
        cs_lnum_t  shft = bld->c2bf_bc_idx[c_id] + count[c_id];

        bld->c2bf_bc_ids[shft] = f_id;
        count[c_id] += 1;

      }

      BFT_FREE(count);

    } // Diffusion part to do
    break;

  default: // Nothing to do
    break;

  } /* Strong enforcement of BCs */

  /* Diffusion part */
  bld->diff = NULL;
  if (bld->todo[CDO_DIFFUSION]) {
    bool is_uniform = cs_property_is_uniform(eqp->diffusion_property);
    bld->diff = cs_cdovb_diffusion_builder_init(connect,
                                                is_uniform,
                                                eqp->diffusion_hodge,
                                                bld->enforce);
  }

  /* Advection part */
  bld->adv = NULL;
  if (bld->todo[CDO_ADVECTION])
    bld->adv = cs_cdovb_advection_builder_init(connect,
                                               eqp->advection_field,
                                               eqp->advection_info,
                                               bld->todo[CDO_DIFFUSION]);

  /* Reaction part */
  bld->hb_reac = NULL;
  if (bld->todo[CDO_REACTION]) {

    BFT_MALLOC(bld->hb_reac, eqp->n_reaction_terms, cs_hodge_builder_t *);

    for (int r = 0; r < eqp->n_reaction_terms; r++)
      bld->hb_reac[r] = cs_hodge_builder_init(connect,
                                              eqp->reaction_terms[r].hodge);

  }

  /* Time-dependent part */
  bld->hb_time = NULL;
  if (bld->todo[CDO_TIME])
    bld->hb_time = cs_hodge_builder_init(connect, eqp->time_hodge);

  /* Source term part */
  BFT_MALLOC(bld->source_terms, bld->n_vertices, cs_real_t);
  bld->build_hvpcd_conf = false;
  bld->hvpcd_conf = NULL;

  return bld;
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
  cs_cdovb_scaleq_t  *bld = (cs_cdovb_scaleq_t *)builder;

  if (bld == NULL)
    return bld;

  /* eqp is only shared. Thies structure is freed later. */
  const cs_equation_param_t  *eqp = bld->eqp;

  /* Free BC structure */
  if (bld->vtx_dir->n_nhmg_elts > 0)
    BFT_FREE(bld->dir_val);

  bld->face_bc = cs_cdo_bc_free(bld->face_bc);
  bld->vtx_dir = cs_cdo_bc_list_free(bld->vtx_dir);

  /* Free Hodge operator defined from conforming reconstruction op. */
  bld->hvpcd_conf = cs_sla_matrix_free(bld->hvpcd_conf);

  /* Renumbering (if strong enforcement of BCs for instance) */
  if (bld->n_vertices > bld->n_dof_vertices) {
    BFT_FREE(bld->v_z2i_ids);
    BFT_FREE(bld->v_i2z_ids);
  }

  BFT_FREE(bld->source_terms);

  /* Free builder sub-structures */
  if (bld->todo[CDO_DIFFUSION]) {
    bld->diff = cs_cdovb_diffusion_builder_free(bld->diff);

    if (bld->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM ||
        bld->enforce ==  CS_PARAM_BC_ENFORCE_WEAK_NITSCHE) {
      BFT_FREE(bld->c2bf_bc_idx);
      BFT_FREE(bld->c2bf_bc_ids);
    }

  }

  if (bld->todo[CDO_ADVECTION])
    bld->adv = cs_cdovb_advection_builder_free(bld->adv);

  if (bld->todo[CDO_REACTION]) {
    for (int r = 0; r < eqp->n_reaction_terms; r++)
      bld->hb_reac[r] = cs_hodge_builder_free(bld->hb_reac[r]);
    BFT_FREE(bld->hb_reac);
  }

  if (bld->todo[CDO_TIME])
    bld->hb_time = cs_hodge_builder_free(bld->hb_time);

  /* Last free */
  BFT_FREE(bld);
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
  double  *st_eval = cs_cdovb_scal_work;

# pragma omp parallel for if (bld->n_vertices > CS_THR_MIN)
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

      double  *mv = cs_cdovb_scal_work + bld->n_vertices;

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
  cs_real_t  ptyval;

  cs_cdovb_scaleq_t  *b = (cs_cdovb_scaleq_t *)builder;

  const cs_cdo_quantities_t  *quant = cs_cdovb_quant;
  const cs_cdo_connect_t  *connect = cs_cdovb_connect;
  const cs_equation_param_t  *eqp = b->eqp;

  /* Default flag value for vertex-based scalar equations */
  cs_flag_t  lm_flag = CS_CDO_LOCAL_V | CS_CDO_LOCAL_E | CS_CDO_LOCAL_EV;

  /* Allocate and initialize a matrix with the larger stencil (that related
     to diffusion => all vertices of a cell are potentially in interaction)
     "adr" means Advection/Diffusion/Reaction
  */
  cs_sla_matrix_t  *sys_mat =
    cs_sla_matrix_create_msr_from_index(cs_cdovb_v2v,
                                        false,  // symmetric
                                        true,   // sorted
                                        1);     // stride

  if (!b->todo[CDO_ADVECTION] && b->enforce != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
    sys_mat->flag |= CS_SLA_MATRIX_SYM;

  /* Preparatory step for diffusion term */
  bool  diff_pty_uniform = true;
  cs_real_33_t  diff_tensor = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  if (b->todo[CDO_DIFFUSION]) {

    diff_pty_uniform = cs_property_is_uniform(eqp->diffusion_property);
    if (diff_pty_uniform)
      cs_property_get_cell_tensor(0, // cell_id
                                  eqp->diffusion_property,
                                  eqp->diffusion_hodge.inv_pty,
                                  diff_tensor);

    if (eqp->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_WBS)
      lm_flag |= CS_CDO_LOCAL_F | CS_CDO_LOCAL_FE;

    if (b->enforce == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
        b->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM)
      lm_flag |=  CS_CDO_LOCAL_F | CS_CDO_LOCAL_FE;

  } // DIFFUSION

  /* Preparatory step for advection term */
  if (b->todo[CDO_ADVECTION])
    lm_flag |=  CS_CDO_LOCAL_F | CS_CDO_LOCAL_FE;

  /* Preparatory step for reaction term */
  bool  *reac_pty_uniform = NULL;
  if (b->todo[CDO_REACTION]) {

    BFT_MALLOC(reac_pty_uniform, eqp->n_reaction_terms, bool);

    for (int r = 0; r < eqp->n_reaction_terms; r++) {

      cs_property_t  *r_pty = eqp->reaction_properties[r];

      reac_pty_uniform[r] = cs_property_is_uniform(r_pty);
      if (reac_pty_uniform[r])
        cs_hodge_builder_set_val(b->hb_reac[r],
                                 cs_property_get_cell_value(0, r_pty));

    }

  } // REACTION

  /* Preparatory step for unsteady term */
  cs_sla_matrix_t  *time_mat = NULL;
  bool time_pty_uniform = true;
  const bool  time_is_diag =
    (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_VORONOI) ? true : false;
  if (b->todo[CDO_TIME]) {

    time_mat = _init_time_matrix(b);
    time_pty_uniform = cs_property_is_uniform(eqp->time_property);
    if (time_pty_uniform) {
      ptyval = cs_property_get_cell_value(0, eqp->time_property);
      cs_hodge_builder_set_val(b->hb_time, ptyval);
    }

  }

  /* Initialize full rhs */
  cs_real_t  *full_rhs = *rhs;
  if (full_rhs == NULL)
    BFT_MALLOC(full_rhs, b->n_vertices, cs_real_t);

# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_vertices; i++)
    full_rhs[i] = 0.0;

  /* Add the contribution of source terms to the full rhs for this time step */
  _add_source_terms(b, full_rhs);

  /* Compute the values of the Dirichlet BC.
     TODO: do the analogy for Neumann BC */
  _compute_dir_values(mesh, field_val, b);

  /* Temporary pre-allocated buffers (the following buffer is used in
     _add_source_terms (be careful if the order of calls is changed) */
  cs_real_t  *dir_bc_vals = cs_cdovb_scal_work;
  cs_flag_t  *cell_flag = connect->c_info->flag;

   /* Initialize arrays */
# pragma omp parallel for if (b->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_vertices; i++)
    dir_bc_vals[i] = 0;

  /* Store the Dirichlet values into an array of size n_vertices */
# pragma omp parallel for if (b->vtx_dir->n_nhmg_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->vtx_dir->n_nhmg_elts; i++)
    dir_bc_vals[b->vtx_dir->elt_ids[i]] = b->dir_val[i];

  /* Main loop on cells to build the linear system */
  /* --------------------------------------------- */

  for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

    /* Set the local mesh structure for the current cell */
    cs_cdo_locmesh_build(c_id, lm_flag, connect, quant, cs_cell_mesh);

    /* Cell-wise view of the linear system to build */
    cs_locmat_t  *adr_mat = cs_cell_sys->mat;

    /* Store the local values attached to Dirichlet values if the current cell
       has at least one border face */
    if (cell_flag[c_id] & CS_CDO_CONNECT_BD)
      for (short int v = 0; v < cs_cell_mesh->n_vc; v++)
        cs_cell_sys->dir_bc[v] = dir_bc_vals[cs_cell_mesh->v_ids[v]];

    /* Initialize the local matrix storing advection/diffusion/reaction terms */
    adr_mat->n_ent = cs_cell_mesh->n_vc;
    for (short int v = 0; v < adr_mat->n_ent; v++)
      adr_mat->ids[v] = cs_cell_mesh->v_ids[v];
    for (short int v = 0; v < adr_mat->n_ent*adr_mat->n_ent; v++)
      adr_mat->val[v] = 0;

    /* DIFFUSION TERM */
    if (b->todo[CDO_DIFFUSION]) { /* Define the local stiffness matrix */

      if (diff_pty_uniform == false)
        cs_property_get_cell_tensor(c_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    diff_tensor);

      cs_locmat_t  *diff_mat = // local matrix owned by the diffusion builder
        cs_cdovb_diffusion_build_local(quant,
                                       cs_cell_mesh,
                                       (const cs_real_3_t (*))diff_tensor,
                                       b->diff);

      cs_locmat_add(adr_mat, diff_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      bft_printf(">> Local diffusion matrix");
      cs_locmat_dump(c_id, diff_mat);
#endif

      /* Weakly enforced Dirichlet BCs for cells attached to the boundary */
      if (b->c2bf_bc_idx != NULL && cell_flag[c_id] & CS_CDO_CONNECT_BD) {

        for (cs_lnum_t j = b->c2bf_bc_idx[c_id];
             j < b->c2bf_bc_idx[c_id+1]; j++) {

          cs_cdovb_diffusion_weak_bc(b->c2bf_bc_ids[j], // border face id
                                     quant,
                                     cs_cell_mesh,
                                     (const cs_real_3_t (*))diff_tensor,
                                     b->diff,
                                     cs_cell_sys); // adr_mat is updated inside

          for (short int v = 0; v < cs_cell_mesh->n_vc; v++)
            full_rhs[cs_cell_mesh->v_ids[v]] += cs_cell_sys->rhs[v];

        } // Loop on border faces

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
      bft_printf(">> Local diffusion matrix after weak enforcement");
      cs_locmat_dump(c_id, diff_mat);
#endif
      } /* Weak enforcement of Dirichlets BCs */

    } /* DIFFUSION */

    /* ADVECTION TERM */
    if (b->todo[CDO_ADVECTION]) { /* Define the local advection matrix */

      cs_locmat_t  *adv_mat =
        cs_cdovb_advection_build_local(cs_cell_mesh,
                                       (const cs_real_3_t (*))diff_tensor,
                                       b->adv);

      cs_locmat_add(adr_mat, adv_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      bft_printf(">> Local advection matrix");
      cs_locmat_dump(c_id, adv_mat);
#endif

      /* Last treatment for the advection term: Apply boundary conditions */
      if (cell_flag[c_id] & CS_CDO_CONNECT_BD) {

        cs_cdovb_advection_add_bc(quant, cs_cell_mesh, b->adv, cs_cell_sys);

        for (short int v = 0; v < cs_cell_mesh->n_vc; v++)
          full_rhs[cs_cell_mesh->v_ids[v]] += cs_cell_sys->rhs[v];

      } // Apply BC

    } /* ADVECTION */

    /* REACTION TERM */
    if (b->todo[CDO_REACTION]) { /* Define the local reaction matrix */

      for (int r = 0; r < eqp->n_reaction_terms; r++) {

        cs_hodge_builder_t  *hb = b->hb_reac[r];
        cs_property_t  *r_pty = eqp->reaction_properties[r];

        if (reac_pty_uniform[r] == false) {
          ptyval = cs_property_get_cell_value(c_id, r_pty);
          cs_hodge_builder_set_val(hb, ptyval);
        }

        cs_locmat_t  *rea_mat = cs_hodge_build_cellwise(quant,
                                                        cs_cell_mesh,
                                                        hb);

        cs_locmat_add(adr_mat, rea_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
        bft_printf(">> Local reaction matrix for term %d", i);
        cs_locmat_dump(c_id, rea_mat);
#endif

      } /* Loop on reaction terms */

    } /* REACTION */

    /* UNSTEADY TERM */
    if (b->todo[CDO_TIME]) { /* Build the mass matrix related to time */

      if (time_pty_uniform == false) {
        ptyval = cs_property_get_cell_value(c_id, eqp->time_property);
        cs_hodge_builder_set_val(b->hb_time, ptyval);
      }

      cs_locmat_t  *time_loc = cs_hodge_build_cellwise(quant,
                                                       cs_cell_mesh,
                                                       b->hb_time);

      /* Assemble the local matrix into the system matrix */
      cs_sla_assemble_msr_sym(time_loc, time_mat, time_is_diag);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      bft_printf(">> Local time matrix");
      cs_locmat_dump(c_id, time_loc);
#endif

    } /* Unsteady term */

    /* Assemble the matrix related to the advcetion/diffusion/reaction terms
       If advection is activated, the resulting system is not symmetric
       Otherwise, the system is symmetric with extra-diagonal terms. */
    if (sys_mat->flag & CS_SLA_MATRIX_SYM)
      cs_sla_assemble_msr_sym(adr_mat, sys_mat, false);
    else
      cs_sla_assemble_msr(adr_mat, sys_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
    bft_printf(">> (FINAL) Local system matrix");
    cs_locmat_dump(c_id, adr_mat);
#endif

  } // Main loop on cells

  /* Unsteady terms have to be considered at the end in order to deal with the
     system matrix fullfilled with diffusion, advection and reaction terms */
  if (b->todo[CDO_TIME]) {

    _apply_time_scheme(field_val, time_mat, dt_cur, b, full_rhs, sys_mat);

    /* sys_mat is now the full system matrix since it now includes the time
       contribution */
    time_mat = cs_sla_matrix_free(time_mat);

  }

  /* Final step in BC management.
     Apply the strong or penalized enforcement. In case of Nitsche enforcement,
     there is nothing to do (already done).
     Must be call after the application of the time scheme */
  _enforce_bc(b, &full_rhs, &sys_mat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
  cs_sla_system_dump("system.log", NULL, sys_mat, full_rhs);
#endif

  if (b->todo[CDO_REACTION])
    BFT_FREE(reac_pty_uniform);

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
 * \brief  Compute the diffusive and convective flux across a list of faces
 *
 * \param[in]       builder    pointer to a builder structure
 * \param[in]       pdi        pointer to an array of field values
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in]       direction  indicate in which direction flux is > 0
 * \param[in, out]  diff_flux  pointer to the value of the diffusive flux
 * \param[in, out]  conv_flux  pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_compute_flux_across_plane(const void          *builder,
                                          const cs_real_t     *pdi,
                                          int                  ml_id,
                                          cs_real_3_t          direction,
                                          double              *diff_flux,
                                          double              *conv_flux)
{
  const cs_cdovb_scaleq_t  *b = (const cs_cdovb_scaleq_t  *)builder;
  const cs_equation_param_t  *eqp = b->eqp;
  const bool  do_adv = eqp->flag & CS_EQUATION_CONVECTION;
  const bool  do_diff = eqp->flag & CS_EQUATION_DIFFUSION;

  cs_mesh_location_type_t  ml_t = cs_mesh_location_get_type(ml_id);

  *diff_flux = 0.;
  *conv_flux = 0.;

  if (pdi == NULL)
    return;

  if (ml_t != CS_MESH_LOCATION_INTERIOR_FACES &&
      ml_t != CS_MESH_LOCATION_BOUNDARY_FACES) {
    cs_base_warn(__FILE__, __LINE__);
    bft_printf(_(" Mesh location type is incompatible with the computation\n"
                 " of the flux across faces.\n"));
    return;
  }

  const cs_cdo_connect_t  *connect = cs_cdovb_connect;
  const cs_sla_matrix_t  *f2c = connect->f2c;
  const cs_cdo_quantities_t  *quant = cs_cdovb_quant;
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  if (n_elts[0] > 0 && elt_ids == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Computing the flux across all interior or border faces is not"
                " managed yet."));

  double  pf;
  cs_real_3_t  gc, pty_gc;
  cs_real_33_t  pty_tens;
  cs_nvec3_t  adv_c;

  if (ml_t == CS_MESH_LOCATION_BOUNDARY_FACES) {

    const cs_lnum_t  n_i_faces = connect->f_info->n_i_elts;
    const cs_lnum_t  shift_if = 2*n_i_faces;

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  bf_id = elt_ids[i];
      const cs_lnum_t  f_id = n_i_faces + bf_id;
      const cs_lnum_t  c_id = f2c->col_id[shift_if + bf_id];
      const cs_quant_t  f = quant->face[f_id];
      const short int  sgn = (_dp3(f.unitv, direction) < 0) ? -1 : 1;
      const double  coef = sgn * f.meas;

      if (do_diff) { /* Compute the local diffusive flux */

        cs_reco_grd_cell_from_pv(c_id, connect, quant, pdi, gc);
        cs_property_get_cell_tensor(c_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);
        cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);
        *diff_flux += -coef * _dp3(f.unitv, pty_gc);

      }

      if (do_adv) { /* Compute the local advective flux */

        cs_advection_field_get_cell_vector(c_id, eqp->advection_field, &adv_c);
        cs_reco_pv_at_face_center(f_id, connect, quant, pdi, &pf);
        *conv_flux += coef * adv_c.meas * _dp3(adv_c.unitv, f.unitv) * pf;

      }

    } // Loop on selected border faces

  }
  else if (ml_t == CS_MESH_LOCATION_INTERIOR_FACES) {

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  f_id = elt_ids[i];
      const cs_lnum_t  shift_f = 2*f_id;
      const cs_lnum_t  c1_id = f2c->col_id[shift_f];
      const cs_lnum_t  c2_id = f2c->col_id[shift_f+1];
      const cs_quant_t  f = quant->face[f_id];
      const short int  sgn = (_dp3(f.unitv, direction) < 0) ? -1 : 1;
      const double  coef = 0.5 * sgn * f.meas;

      if (do_diff) { /* Compute the local diffusive flux */

        cs_reco_grd_cell_from_pv(c1_id, connect, quant, pdi, gc);
        cs_property_get_cell_tensor(c1_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);
        cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);
        *diff_flux += -coef * _dp3(f.unitv, pty_gc);

        cs_reco_grd_cell_from_pv(c2_id, connect, quant, pdi, gc);
        cs_property_get_cell_tensor(c2_id,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);
        cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);
        *diff_flux += -coef * _dp3(f.unitv, pty_gc);

      }

      if (do_adv) { /* Compute the local advective flux */

        cs_reco_pv_at_face_center(f_id, connect, quant, pdi, &pf);

        cs_advection_field_get_cell_vector(c1_id, eqp->advection_field, &adv_c);
        *conv_flux += coef * adv_c.meas * _dp3(adv_c.unitv, f.unitv) * pf;

        cs_advection_field_get_cell_vector(c2_id, eqp->advection_field, &adv_c);
        *conv_flux += coef * adv_c.meas * _dp3(adv_c.unitv, f.unitv) * pf;

      }

    } // Loop on selected interior faces

  } // Set of interior or border faces
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

    cs_real_t  *work_c = cs_cdovb_scal_work;
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

      cs_cdovb_advection_get_peclet_cell(cs_cdovb_quant,
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
                          cs_time_step);  // time step management structure

      if (eqp->process_flag & CS_EQUATION_POST_UPWIND_COEF) {

        if (k == 0)
          sprintf(postlabel, "%s.UpwCoefX", eqname);
        else if (k == 1)
          sprintf(postlabel, "%s.UpwCoefY", eqname);
        else
          sprintf(postlabel, "%s.UpwCoefZ", eqname);

        /* Compute in each cell an evaluation of upwind weight value */
        cs_cdovb_advection_get_upwind_coef_cell(cs_cdovb_quant,
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
                          cs_time_step);  // time step management structure

      } /* Post upwinding coefficient */

    } // Loop on space dimension

    BFT_FREE(postlabel);

  } // Post a Peclet attached to cells

}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
