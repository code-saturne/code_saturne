/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
 * reaction of scalar-valued equations with source terms
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_local.h"
#include "cs_cdo_time.h"
#include "cs_cdovb_priv.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_evaluate.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
#include "cs_param.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_scheme_geometry.h"
#include "cs_search.h"
#include "cs_sles.h"
#include "cs_source_term.h"
#include "cs_timer.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_cdovb_scaleq.c

  \brief Build an algebraic CDO vertex-based system for unsteady
         convection-diffusion-reaction of scalar-valued equations with
         source terms

*/

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVB_SCALEQ_DBG     0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Algebraic system for CDO vertex-based discretization */
typedef struct _cs_cdovb_t cs_cdovb_scaleq_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Structure to enable a full cellwise strategy during the system building */
static cs_cell_sys_t      **_vbs_cell_system = NULL;
static cs_cell_builder_t  **_vbs_cell_builder = NULL;

/* Pointer to shared structures */
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;
static const cs_matrix_structure_t  *cs_shared_ms;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local builder structure used for building the system
 *          cellwise
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_cell_builder_t *
_vbs_create_cell_builder(const cs_cdo_connect_t   *connect)
{
  const int  n_vc = connect->n_max_vbyc;
  const int  n_ec = connect->n_max_ebyc;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  BFT_MALLOC(cb->ids, n_vc, short int);
  memset(cb->ids, 0, n_vc*sizeof(short int));

  int  size = n_ec*(n_ec+1);
  size = CS_MAX(4*n_ec + 3*n_vc, size);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(cs_real_t));

  size = 2*n_ec;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of
     operators */
  cb->hdg = cs_sdm_square_create(n_ec);
  cb->loc = cs_sdm_square_create(n_vc);
  cb->aux = cs_sdm_square_create(n_vc);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *         Define an indirection array for the enforcement of internal DoFs
 *         only if needed.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]      t_eval          time at which one evaluates BCs
 * \param[in]      mesh            pointer to a cs_mesh_t structure
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in]      eqb             pointer to a cs_equation_builder_t structure
 * \param[in, out] vtx_bc_flag      pointer to an array of BC flag for each vtx
 * \param[in, out] p_dir_values    pointer to the Dirichlet values to set
 * \param[in, out] p_enforced_ids  pointer to the list of enforced vertices
 */
/*----------------------------------------------------------------------------*/

static void
_vbs_setup(cs_real_t                      t_eval,
           const cs_mesh_t               *mesh,
           const cs_equation_param_t     *eqp,
           const cs_equation_builder_t   *eqb,
           cs_flag_t                      vtx_bc_flag[],
           cs_real_t                     *p_dir_values[],
           cs_lnum_t                     *p_enforced_ids[])
{
  assert(vtx_bc_flag != NULL);  /* Sanity check */
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_real_t  *dir_values = NULL;

  /* Compute the values of the Dirichlet BC */
  BFT_MALLOC(dir_values, quant->n_vertices, cs_real_t);

  /* Compute the values of the Dirichlet BC */
  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   cs_shared_connect,
                                   eqp,
                                   eqb->face_bc,
                                   _vbs_cell_builder[0], /* static variable */
                                   vtx_bc_flag,
                                   dir_values);


  *p_dir_values = dir_values;

  /* Internal enforcement of DoFs  */
  if (cs_equation_param_has_internal_enforcement(eqp)) {

    cs_lnum_t  *enforced_ids = NULL;
    BFT_MALLOC(enforced_ids, quant->n_vertices, cs_lnum_t);
    for (cs_lnum_t i = 0; i < quant->n_vertices; i++)
      enforced_ids[i] = -1;     /* Not selected */

    for (cs_lnum_t i = 0; i < eqp->n_enforced_dofs; i++) {
      cs_lnum_t  id = eqp->enforced_dof_ids[i];
      enforced_ids[id] = i;
    }

    *p_enforced_ids = enforced_ids;
  }
  else
    *p_enforced_ids = NULL;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the local structure for the current cell
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in]      cell_flag    flag related to the current cell
 * \param[in]      cm           pointer to a cellwise view of the mesh
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in]      dir_values   Dirichlet values associated to each vertex
 * \param[in]      vtx_bc_flag  Flag related to BC associated to each vertex
 * \param[in]      forced_ids   indirection in case of internal enforcement
 * \param[in]      field_tn     values of the field at the last computed time
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vbs_init_cell_system(cs_real_t                      t_eval,
                      const cs_flag_t                cell_flag,
                      const cs_cell_mesh_t          *cm,
                      const cs_equation_param_t     *eqp,
                      const cs_equation_builder_t   *eqb,
                      const cs_real_t                dir_values[],
                      const cs_flag_t                vtx_bc_flag[],
                      const cs_lnum_t                forced_ids[],
                      const cs_real_t                field_tn[],
                      cs_cell_sys_t                 *csys,
                      cs_cell_builder_t             *cb)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Cell-wise view of the linear system to build */
  csys->c_id = cm->c_id;
  csys->cell_flag = cell_flag;
  csys->n_dofs = cm->n_vc;

  /* Initialize the local system */
  cs_cell_sys_reset(cm->n_fc, csys); /* Generic part */

  cs_sdm_square_init(cm->n_vc, csys->mat);

  for (short int v = 0; v < cm->n_vc; v++) {
    csys->dof_ids[v] = cm->v_ids[v];
    csys->val_n[v] = field_tn[cm->v_ids[v]];
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Set the bc (specific part) */
    cs_equation_vb_set_cell_bc(cm,
                               connect,
                               cs_shared_quant,
                               eqp,
                               eqb->face_bc,
                               vtx_bc_flag,
                               dir_values,
                               t_eval,
                               csys,
                               cb);

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    cs_dbg_check_hmg_dirichlet_cw(__func__, csys);
#endif
  } /* Border cell */

  /* Special case to handle if enforcement by penalization or algebraic
   * This situation may happen with a tetrahedron with one vertex or an edge
   * lying on the boundary (but no face)
   */
  if (cell_flag == CS_FLAG_BOUNDARY_CELL_BY_VERTEX) {

    assert(vtx_bc_flag != NULL);

    for (short int v = 0; v < cm->n_vc; v++) {
      csys->dof_flag[v] = vtx_bc_flag[cm->v_ids[v]];
      if (cs_cdo_bc_is_dirichlet(csys->dof_flag[v])) {
        csys->has_dirichlet = true;
        csys->dir_values[v] = dir_values[cm->v_ids[v]];
      }
    }

  }

  /* Internal enforcement of DoFs  */
  if (cs_equation_param_has_internal_enforcement(eqp)) {

    assert(forced_ids != NULL);
    for (short int v = 0; v < cm->n_vc; v++) {

      const cs_lnum_t  id = forced_ids[cm->v_ids[v]];

      /* In case of a Dirichlet BC, this BC is applied and the enforcement
         is ignored */
      if (cs_cdo_bc_is_dirichlet(csys->dof_flag[v]))
        csys->intern_forced_ids[v] = -1;
      else {
        csys->intern_forced_ids[v] = id;
        if (id > -1)
          csys->has_internal_enforcement = true;
      }

    } /* Loop on cell vertices */

  }

  /* Set the properties for this cell if not uniform */
  cs_equation_init_properties_cw(eqp, eqb, t_eval, cell_flag, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the local matrices arising from the diffusion, advection,
 *         reaction terms. If asked, a mass matrix is also computed and stored
 *         in cb->hdg.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]      time_eval   time at which analytic function are evaluated
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vbs_advection_diffusion_reaction(double                         time_eval,
                                  const cs_equation_param_t     *eqp,
                                  const cs_equation_builder_t   *eqb,
                                  const cs_cdovb_scaleq_t       *eqc,
                                  const cs_cell_mesh_t          *cm,
                                  cs_face_mesh_t                *fm,
                                  cs_cell_sys_t                 *csys,
                                  cs_cell_builder_t             *cb)
{
  if (cs_equation_param_has_diffusion(eqp)) {   /* DIFFUSION TERM
                                                 * ============== */

    /* Define the local stiffness matrix: local matrix owned by the cellwise
       builder (store in cb->loc) */
    eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

    /* Add the local diffusion operator to the local system */
    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after adding diffusion", csys);
#endif
  }

  if (cs_equation_param_has_convection(eqp)) {  /* ADVECTION TERM
                                                 * ============== */

    /* Define the local advection matrix */
    eqc->get_advection_matrix(eqp, cm, time_eval, fm, cb);

    /* Add it to the local system */
    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after adding advection", csys);
#endif
  }

  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                  * =========== */

    /* Build the mass matrix and store it in cb->hdg */
    eqc->get_mass_matrix(eqc->hdg_mass, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys)) {
      cs_log_printf(CS_LOG_DEFAULT, ">> Cell mass matrix");
      cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, cb->hdg);
    }
#endif
  }

  if (cs_equation_param_has_reaction(eqp)) { /* REACTION TERM
                                              * ============= */

    if (eqb->sys_flag & CS_FLAG_SYS_REAC_DIAG) {

      /* |c|*wvc = |dual_cell(v) cap c| */
      assert(cs_flag_test(eqb->msh_flag, CS_CDO_LOCAL_PVQ));
      const double  ptyc = cb->rpty_val * cm->vol_c;
      for (short int i = 0; i < cm->n_vc; i++)
        csys->mat->val[i*(cm->n_vc + 1)] += cm->wvc[i] * ptyc;

    }
    else {

      assert(cs_flag_test(eqb->sys_flag, CS_FLAG_SYS_MASS_MATRIX));

      /* Update local system matrix with the reaction term
         cb->hdg corresponds to the current mass matrix */
      cs_sdm_add_mult(csys->mat, cb->rpty_val, cb->hdg);

    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after adding reaction", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  First pass to apply boundary conditions enforced weakly in CDO-Vb
 *         schemes. Update the local system before applying the time scheme.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]      time_eval   time at which analytical function are evaluated
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vbs_apply_weak_bc(cs_real_t                      time_eval,
                   const cs_equation_param_t     *eqp,
                   const cs_cdovb_scaleq_t       *eqc,
                   const cs_cell_mesh_t          *cm,
                   cs_face_mesh_t                *fm,
                   cs_cell_sys_t                 *csys,
                   cs_cell_builder_t             *cb)
{
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Neumann boundary conditions */
    if (csys->has_nhmg_neumann) {
      for (short int v  = 0; v < cm->n_vc; v++) {
        if (cs_cdo_bc_is_dirichlet(csys->dof_flag[v]) == false)
          csys->rhs[v] += csys->neu_values[v];
      }
    }

    /* Contribution for the advection term: csys is updated inside
       (matrix and rhs) and Dirichlet BCs are handled inside */
    if (cs_equation_param_has_convection(eqp))
      eqc->add_advection_bc(cm, eqp, time_eval, fm, cb, csys);

    /* The enforcement of the Dirichlet has to be done after all
       other contributions */
    if (cs_equation_param_has_diffusion(eqp)) {

      if (csys->has_robin) {
        assert(eqc->enforce_robin_bc != NULL);
        eqc->enforce_robin_bc(eqp, cm, fm, cb, csys);
      }

      if (csys->has_dirichlet) /* csys is updated inside (matrix and rhs) */
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
          eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after BC treatment", csys);
#endif
  } /* Cell with at least one boundary face */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Second pass to apply boundary conditions. Only Dirichlet BCs which
 *         are enforced strongly. Apply also the enforcement of internal DoFs.
 *         Update the local system after applying the time scheme.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vbs_enforce_values(const cs_equation_param_t     *eqp,
                    const cs_cdovb_scaleq_t       *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_face_mesh_t                *fm,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  if (csys->cell_flag > 0 && csys->has_dirichlet) {

    /* Boundary element (through either vertices or faces) */

    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED) {

      /* csys is updated inside (matrix and rhs) */
      eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after strong BC treatment", csys);
#endif
    }
  }

  if (cs_equation_param_has_internal_enforcement(eqp) == false)
    return;

  /* Internal enforcement of DoFs: Update csys (matrix and rhs) */
  if (csys->has_internal_enforcement) {

    cs_equation_enforced_internal_dofs(eqp, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement",
                       csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the renormalization coefficient for the
 *         the residual norm of the linear system
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cellwise view of the mesh/quantities
 * \param[in]      csys      pointer to a cellwise view of the system
 * \param[in, out] rhs_norm  quantity used for the RHS normalization
 */
/*----------------------------------------------------------------------------*/

static inline void
_vbs_compute_cw_sles_normalization(const cs_equation_param_t    *eqp,
                                   const cs_cell_mesh_t         *cm,
                                   const cs_cell_sys_t          *csys,
                                   cs_real_t                    *rhs_norm)
{
  if (eqp->sles_param.resnorm_type == CS_PARAM_RESNORM_WEIGHTED_RHS) {

    cs_real_t  _rhs_norm = 0;
    for (short int v = 0; v < cm->n_vc; v++)
      _rhs_norm += cm->wvc[v] * csys->rhs[v]*csys->rhs[v];

    *rhs_norm += cm->vol_c * _rhs_norm;

  }
  else if (eqp->sles_param.resnorm_type == CS_PARAM_RESNORM_MAT_DIAG) {

    cs_real_t  _rhs_norm = 0;
    for (short int v = 0; v < cm->n_vc; v++) {
      const double  d_val = csys->mat->val[v*(cm->n_vc+1)];
      _rhs_norm += cm->wvc[v] * d_val *d_val;
    }

    *rhs_norm += cm->vol_c * _rhs_norm;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last stage to compute of the renormalization coefficient for the
 *         the residual norm of the linear system
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in, out] rhs_norm  quantity used for the RHS normalization
 */
/*----------------------------------------------------------------------------*/

static inline void
_vbs_sync_sles_normalization(const cs_equation_param_t    *eqp,
                             cs_real_t                    *rhs_norm)
{
  cs_parall_sum(1, CS_DOUBLE, rhs_norm);

  switch (eqp->sles_param.resnorm_type) {

  case CS_PARAM_RESNORM_WEIGHTED_RHS:
  case CS_PARAM_RESNORM_MAT_DIAG:
    *rhs_norm = sqrt(1/cs_shared_quant->vol_tot*(*rhs_norm));
    if (*rhs_norm < 10*FLT_MIN)
      *rhs_norm = cs_shared_quant->vol_tot/cs_shared_quant->n_g_cells;
    break;

  case CS_PARAM_RESNORM_VOLTOT:
    *rhs_norm = cs_shared_quant->vol_tot/cs_shared_quant->n_g_cells;
    break;

  default:
    *rhs_norm = 1.0;
    break;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from a scalar-valued CDO-Vb scheme
 *
 * \param[in, out] sles      pointer to a cs_sles_t structure
 * \param[in]      matrix    pointer to a cs_matrix_t structure
 * \param[in]      field_id id related to the variable field of this equation
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      rhs_norm  quantity used for the RHS normalization
 * \param[in, out] x         solution of the linear system (in: initial guess)
 * \param[in, out] b         right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

static int
_vbs_solve_system(cs_sles_t                    *sles,
                  const cs_matrix_t            *matrix,
                  const int                     field_id,
                  const cs_equation_param_t    *eqp,
                  const cs_real_t               rhs_norm,
                  cs_real_t                    *x,
                  cs_real_t                    *b)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;

  /* solving info */
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_solving_info_t sinfo;
  cs_field_get_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  sinfo.n_it = 0;
  sinfo.res_norm = DBL_MAX;
  cs_range_set_t  *rset = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  cs_real_t  *xsol = NULL;

  const cs_lnum_t  n_scatter_elts = n_vertices;
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);

  if (n_cols > n_scatter_elts) {
    assert(cs_glob_n_ranks > 1);
    BFT_MALLOC(xsol, n_cols, cs_real_t);
    memcpy(xsol, x, n_scatter_elts*sizeof(cs_real_t));
  }
  else
    xsol = x;

  /* Prepare solving (handle parallelism) */
  cs_gnum_t  nnz = cs_equation_prepare_system(1,            /* stride */
                                              n_scatter_elts,
                                              matrix,
                                              rset,
                                              xsol, b);

  /* Solve the linear solver */
  sinfo.rhs_norm = rhs_norm;
  const cs_param_sles_t  sles_param = eqp->sles_param;

  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    sles_param.eps,
                                                    sinfo.rhs_norm,
                                                    &(sinfo.n_it),
                                                    &(sinfo.res_norm),
                                                    b,
                                                    xsol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  if (sles_param.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%s/sles_cvg> code %-d n_iters %d"
                  " residual % -8.4e nnz %lu\n",
                  eqp->name, code, sinfo.n_it, sinfo.res_norm, nnz);

  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         xsol, x);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         b, b);

  cs_dbg_fprintf_system(eqp->name, cs_shared_time_step->nt_cur,
                        CS_CDOVB_SCALEQ_DBG,
                        x, b, n_vertices);
#endif

  /* Free what can be freed at this stage */
  cs_sles_free(sles);

  if (n_cols > n_scatter_elts)
    BFT_FREE(xsol);

  cs_field_set_key_struct(fld, cs_field_key_id("solving_info"), &sinfo);

  return (sinfo.n_it);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-Vb scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdovb_scaleq_is_initialized(void)
{
  if (_vbs_cell_system == NULL || _vbs_cell_builder == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Allocate work buffer and general structures related to CDO
 *           vertex-based schemes
 *           Set shared pointers.
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_init_common(const cs_cdo_quantities_t    *quant,
                            const cs_cdo_connect_t       *connect,
                            const cs_time_step_t         *time_step,
                            const cs_matrix_structure_t  *ms)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
  cs_shared_ms = ms;

  /* Structure used to build the final system by a cell-wise process */
  assert(cs_glob_n_threads > 0);  /* Sanity check */

  BFT_MALLOC(_vbs_cell_system, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(_vbs_cell_builder, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    _vbs_cell_system[i] = NULL;
    _vbs_cell_builder[i] = NULL;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    _vbs_cell_system[t_id] = cs_cell_sys_create(connect->n_max_vbyc,
                                                 connect->n_max_fbyc,
                                                 1, NULL);
    _vbs_cell_builder[t_id] = _vbs_create_cell_builder(connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  _vbs_cell_system[0] = cs_cell_sys_create(connect->n_max_vbyc,
                                            connect->n_max_fbyc,
                                            1, NULL);
  _vbs_cell_builder[0] = _vbs_create_cell_builder(connect);

#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  *csys = _vbs_cell_system[t_id];
  *cb = _vbs_cell_builder[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO vertex-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_finalize_common(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(_vbs_cell_system[t_id]));
    cs_cell_builder_free(&(_vbs_cell_builder[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(_vbs_cell_system[0]));
  cs_cell_builder_free(&(_vbs_cell_builder[0]));
#endif /* openMP */

  BFT_FREE(_vbs_cell_system);
  BFT_FREE(_vbs_cell_builder);
  _vbs_cell_system = NULL;
  _vbs_cell_builder = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_scaleq_t structure storing data useful
 *         for building and managing such a scheme
 *
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id      id of the variable field
 * \param[in]      bflux_id    id of the boundary flux field
 * \param[in, out] eqb         pointer to a \ref cs_equation_builder_t struct.
 *
 * \return a pointer to a new allocated cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_scaleq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVB && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of equation.\n"
              " Expected: scalar-valued CDO vertex-based equation.", __func__);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->n_vertices;

  cs_cdovb_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdovb_scaleq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  eqc->n_dofs = n_vertices;

  /* Flag to indicate the minimal set of quantities to build in a cell mesh
     According to the situation, additional flags have to be set */
  eqb->msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PE |
    CS_CDO_LOCAL_EV;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
    CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_FV;

  /* Diffusion */
  eqc->get_stiffness_matrix = NULL;
  eqc->enforce_robin_bc = NULL;
  if (cs_equation_param_has_diffusion(eqp)) {

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_cost_get_stiffness;

      eqb->bd_msh_flag |= CS_CDO_LOCAL_DEQ;
      eqc->enforce_robin_bc = cs_cdo_diffusion_svb_cost_robin;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_voro_get_stiffness;

      eqb->bd_msh_flag |= CS_CDO_LOCAL_DEQ;
      eqc->enforce_robin_bc = cs_cdo_diffusion_svb_cost_robin;
      break;

    case CS_PARAM_HODGE_ALGO_WBS:
      eqb->msh_flag |= CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_PEQ
        | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_wbs_get_stiffness;
      eqc->enforce_robin_bc = cs_cdo_diffusion_svb_wbs_robin;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                (" %s: Invalid type of algorithm to build the diffusion term."),
                __func__);

    } /* Switch on Hodge algo. */

  } /* Diffusion */

  /* Boundary conditions */
  BFT_MALLOC(eqc->vtx_bc_flag, n_vertices, cs_flag_t);
  cs_equation_set_vertex_bc_flag(connect, eqb->face_bc, eqc->vtx_bc_flag);

  eqc->enforce_dirichlet = NULL;
  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    eqb->bd_msh_flag |= CS_CDO_LOCAL_DEQ;
    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqc->enforce_dirichlet = cs_cdo_diffusion_svb_cost_weak_dirichlet;
      break;
    case CS_PARAM_HODGE_ALGO_WBS:
      eqc->enforce_dirichlet = cs_cdo_diffusion_svb_wbs_weak_dirichlet;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to enforce the Dirichlet BC.",
                __func__);

    } /* Switch on Hodge algo. */
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    eqb->bd_msh_flag |= CS_CDO_LOCAL_DEQ;
    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqc->enforce_dirichlet = cs_cdo_diffusion_svb_cost_wsym_dirichlet;
      break;
    case CS_PARAM_HODGE_ALGO_WBS:
      eqc->enforce_dirichlet = cs_cdo_diffusion_svb_wbs_wsym_dirichlet;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to enforce the Dirichlet BC.",
                __func__);

    } /* Switch on Hodge algo. */
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  /* Advection */
  eqc->get_advection_matrix = NULL;
  eqc->add_advection_bc = NULL;

  if (cs_equation_param_has_convection(eqp)) {

    cs_xdef_type_t  adv_deftype =
      cs_advection_field_get_deftype(eqp->adv_field);

    if (adv_deftype == CS_XDEF_BY_VALUE)
      eqb->msh_flag |= CS_CDO_LOCAL_DFQ;
    else if (adv_deftype == CS_XDEF_BY_ARRAY)
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ;
    else if (adv_deftype == CS_XDEF_BY_ANALYTIC_FUNCTION)
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_EFQ | CS_CDO_LOCAL_PFQ;

    switch (eqp->adv_formulation) {

    case CS_PARAM_ADVECTION_FORM_CONSERV:

      switch (eqp->adv_scheme) {

      case CS_PARAM_ADVECTION_SCHEME_CENTERED:
        eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
        eqc->get_advection_matrix = cs_cdo_advection_vb_cencsv;
        break;

      case CS_PARAM_ADVECTION_SCHEME_MIX_CENTERED_UPWIND:
        eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
        eqc->get_advection_matrix = cs_cdo_advection_vb_mcucsv;
        break;

      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
      case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
      case CS_PARAM_ADVECTION_SCHEME_SG:
        eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
        if (cs_equation_param_has_diffusion(eqp))
          eqc->get_advection_matrix = cs_cdo_advection_vb_upwcsv_di;
        else
          eqc->get_advection_matrix = cs_cdo_advection_vb_upwcsv;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid advection scheme for vertex-based discretization");
      } /* Scheme */
      break; /* Formulation */

    case CS_PARAM_ADVECTION_FORM_NONCONS:

      switch (eqp->adv_scheme) {
      case CS_PARAM_ADVECTION_SCHEME_CENTERED:
        eqc->get_advection_matrix = cs_cdo_advection_vb_cennoc;
        break;

      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
      case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
      case CS_PARAM_ADVECTION_SCHEME_SG:
        eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
        if (cs_equation_param_has_diffusion(eqp))
          eqc->get_advection_matrix = cs_cdo_advection_vb_upwnoc_di;
        else
          eqc->get_advection_matrix = cs_cdo_advection_vb_upwnoc;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid advection scheme for vertex-based discretization");
      } /* Scheme */
      break; /* Formulation */

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid type of formulation for the advection term");
    }

    /* Boundary conditions for advection */
    eqb->bd_msh_flag |= CS_CDO_LOCAL_PEQ;
    eqc->add_advection_bc = cs_cdo_advection_vb_bc;

  }
  else {

    if (eqp->default_enforcement != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
      eqb->sys_flag |= CS_FLAG_SYS_SYM; /* Algebraic system is symmetric */

  }

  /* Reaction */
  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->do_lumping)
      eqb->sys_flag |= CS_FLAG_SYS_REAC_DIAG;
    else {

      switch (eqp->reaction_hodge.algo) {

      case CS_PARAM_HODGE_ALGO_VORONOI:
        eqb->sys_flag |= CS_FLAG_SYS_REAC_DIAG;
        break;
      case CS_PARAM_HODGE_ALGO_WBS:
        eqb->msh_flag |= CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_PEQ
          | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid choice of algorithm for the reaction term.",
                  __func__);
        break;
      }

    } /* Lumping or not lumping */

  } /* Reaction */

  /* Time */
  if (cs_equation_param_has_time(eqp)) {

    if (eqp->do_lumping)
      eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
    else {

      switch (eqp->time_hodge.algo) {

      case CS_PARAM_HODGE_ALGO_VORONOI:
        eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
        break;
      case CS_PARAM_HODGE_ALGO_WBS:
        eqb->msh_flag |= CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_PEQ
          | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid choice of algorithm for the time term.",
                  __func__);
        break;
      }

    } /* Lumping or not lumping */

  } /* Time part */

  /* Source term */
  eqc->source_terms = NULL;

  if (cs_equation_param_has_sourceterm(eqp)) {

    /* When the deprecated mode will be removed. This array needs to be
       allocated only if a theta scheme is used */
    BFT_MALLOC(eqc->source_terms, eqc->n_dofs, cs_real_t);
#   pragma omp parallel for if (eqc->n_dofs > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < eqc->n_dofs; i++)
      eqc->source_terms[i] = 0;

  } /* There is at least one source term */

  /* Pre-defined a cs_hodge_builder_t structure */
  eqc->hdg_mass.is_unity = true;
  eqc->hdg_mass.is_iso   = true;
  eqc->hdg_mass.inv_pty  = false;
  eqc->hdg_mass.type = CS_PARAM_HODGE_TYPE_VPCD;
  eqc->hdg_mass.algo = CS_PARAM_HODGE_ALGO_WBS;
  eqc->hdg_mass.coef = 1.0; /* not useful in this case */

  eqc->get_mass_matrix = cs_hodge_vpcd_wbs_get;

  /* Assembly process */
  eqc->assemble = cs_equation_assemble_matrix;

  /* Array used for extra-operations */
  eqc->cell_values = NULL;

  return eqc;
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
cs_cdovb_scaleq_free_context(void   *builder)
{
  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)builder;

  if (eqc == NULL)
    return eqc;

  /* These arrays may have not been allocated */
  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->cell_values);
  BFT_FREE(eqc->vtx_bc_flag);

  /* Last free */
  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-Vb schemes.
 *
 * \param[in]      t_eval     time at which one evaluates BCs
 * \param[in]      field_id   id related to the variable field of this equation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to the scheme context (cast on-the-fly)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_real_t  *v_vals = fld->val;

  /* By default, 0 is set as initial condition for the computational domain.

     Warning: This operation has to be done after the settings of the
     Dirichlet boundary conditions where an interface sum is performed
     for vertex-based schemes
  */

  memset(v_vals, 0, quant->n_vertices*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    /* Initialize values at mesh vertices */
    cs_flag_t  dof_flag = CS_FLAG_SCALAR | cs_flag_primal_vtx;

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_defs[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(dof_flag, def, v_vals);
        break;

      case CS_XDEF_BY_QOV:
        cs_evaluate_potential_by_qov(dof_flag, def, v_vals, NULL);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        assert(eqp->dof_reduction == CS_PARAM_REDUCTION_DERHAM);
        cs_evaluate_potential_by_analytic(dof_flag, def, t_eval, v_vals);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid way to initialize field values for eq. %s.\n",
                  __func__, eqp->name);

      } /* Switch on possible type of definition */

    } /* Loop on definitions */

  } /* Initial values to set */

  /* Set the boundary values as initial values: Compute the values of the
     Dirichlet BC */
  cs_real_t  *work_v = cs_equation_get_tmpbuf();

  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   _vbs_cell_builder[0], /* static variable */
                                   eqc->vtx_bc_flag,
                                   work_v);


  for (cs_lnum_t v = 0; v < quant->n_vertices; v++)
    if (cs_cdo_bc_is_dirichlet(eqc->vtx_bc_flag[v]))
      v_vals[v] = work_v[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_steady_state(const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Build an array storing the Dirichlet values at vertices and another one
     with a tags to detect vertices related to a Neumann BC */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *forced_ids = NULL;

  /* First argument is set to t_cur even if this is a steady computation since
   * one can call this function to compute a steady-state solution at each time
   * step of an unsteady computation. */
  _vbs_setup(time_eval, mesh, eqp, eqb, eqc->vtx_bc_flag,
             &dir_values, &forced_ids);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;
  cs_real_t  rhs_norm = 0.0;

  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_equation_get_mav(matrix, eqp->omp_assembly_choice, 1);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, dir_values,   \
         forced_ids, fld, rs, _vbs_cell_system, _vbs_cell_builder, rhs_norm)
  {
    /* Set variables and structures inside the OMP section so that each thread
       has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_matrix_assembler_buf_t  *mab = cs_equation_get_assembly_buffers(t_id);
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = _vbs_cell_system[t_id];
    cs_cell_builder_t  *cb = _vbs_cell_builder[t_id];

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _vbs_init_cell_system(time_eval, cell_flag, cm, eqp, eqb,
                            dir_values, eqc->vtx_bc_flag, forced_ids, fld->val,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vbs_advection_diffusion_reaction(time_eval,
                                        eqp, eqb, eqc, cm, fm, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */
        /* Reset the local contribution */
        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system */
        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        time_eval,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        /* Update the RHS */
        for (short int v = 0; v < cm->n_vc; v++)
          csys->rhs[v] += csys->source[v];

      } /* End of term source */

      /* Compute a norm of the RHS for the normalization of the SLES */
      _vbs_compute_cw_sles_normalization(eqp, cm, csys, &rhs_norm);

      /* Apply boundary conditions (those which are weakly enforced) */
      _vbs_apply_weak_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

      /* Enforce values if needed (internal or Dirichlet) */
      _vbs_enforce_values(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      eqc->assemble(csys, rs, mab, mav);       /* Matrix assembly */

      for (short int v = 0; v < cm->n_vc; v++) /* RHS assembly*/
#       pragma omp atomic
        rhs[cm->v_ids[v]] += csys->rhs[v];

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(forced_ids);

  cs_matrix_assembler_values_finalize(&mav);

  /* Last step in the computation of the renormalization coefficient */
  _vbs_sync_sles_normalization(eqp, &rhs_norm);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Now solve the system */
  _vbs_solve_system(cs_sles_find_or_add(field_id, NULL),
                    matrix,
                    field_id,
                    eqp,
                    rhs_norm,
                    fld->val,
                    rhs);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         Implicit time scheme is used to progress in time.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_implicit(const cs_mesh_t            *mesh,
                               const int                   field_id,
                               const cs_equation_param_t  *eqp,
                               cs_equation_builder_t      *eqb,
                               void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  dt_cur = cs_shared_time_step->dt[0];
  const cs_real_t  time_eval = t_cur + dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Build an array storing the Dirichlet values at vertices and another one
     with a tags to detect vertices related to a Neumann BC */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *forced_ids = NULL;

  _vbs_setup(t_cur + dt_cur, mesh, eqp, eqb, eqc->vtx_bc_flag,
             &dir_values,  &forced_ids);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;
  cs_real_t  rhs_norm = 0.;

  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_equation_get_mav(matrix, eqp->omp_assembly_choice, 1);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, dir_values,   \
         forced_ids, fld, rs, _vbs_cell_system, _vbs_cell_builder, rhs_norm)
  {
    /* Set variables and structures inside the OMP section so that each thread
       has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_matrix_assembler_buf_t  *mab = cs_equation_get_assembly_buffers(t_id);
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = _vbs_cell_system[t_id];
    cs_cell_builder_t  *cb = _vbs_cell_builder[t_id];

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _vbs_init_cell_system(time_eval, cell_flag, cm, eqp, eqb,
                           dir_values, eqc->vtx_bc_flag, forced_ids, fld->val,
                           csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vbs_advection_diffusion_reaction(time_eval,
                                        eqp, eqb, eqc, cm, fm, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */

        /* Reset the local contribution */
        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */
        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        time_eval,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        for (short int v = 0; v < cm->n_vc; v++)
          csys->rhs[v] += csys->source[v];

      } /* End of term source */

      /* Apply boundary conditions (those which are weakly enforced) */
      _vbs_apply_weak_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        /* |c|*wvc = |dual_cell(v) cap c| */
        assert(cs_flag_test(eqb->msh_flag, CS_CDO_LOCAL_PVQ));
        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
         *       >> Update the cellwise system with the time matrix */
        for (short int i = 0; i < cm->n_vc; i++) {

          const double  dval =  ptyc * cm->wvc[i];

          /* Update the RHS with values at time t_n */
          csys->rhs[i] += dval * csys->val_n[i];

          /* Add the diagonal contribution from time matrix */
          csys->mat->val[i*(cm->n_vc + 1)] += dval;

        }

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val * inv_dtcur;
        const cs_sdm_t  *mass_mat = cb->hdg;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
         *       >> Update the cellwise system with the time matrix */

        /* Update rhs with csys->mat*p^n */
        double  *time_pn = cb->values;
        cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);
        for (short int i = 0; i < csys->n_dofs; i++)
          csys->rhs[i] += tpty_coef*time_pn[i];

        /* Update the cellwise system with the time matrix */
        cs_sdm_add_mult(csys->mat, tpty_coef, mass_mat);

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after time", csys);
#endif

      /* Compute a norm of the RHS for the normalization of the SLES */
      _vbs_compute_cw_sles_normalization(eqp, cm, csys, &rhs_norm);

      /* Enforce values if needed (internal or Dirichlet) */
      _vbs_enforce_values(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      eqc->assemble(csys, rs, mab, mav);       /* Matrix assembly */

      for (short int v = 0; v < cm->n_vc; v++) /* RHS assembly */
#       pragma omp atomic
        rhs[cm->v_ids[v]] += csys->rhs[v];

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(forced_ids);
  cs_matrix_assembler_values_finalize(&mav);

  /* Last step in the computation of the renormalization coefficient */
  _vbs_sync_sles_normalization(eqp, &rhs_norm);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Now solve the system */
  _vbs_solve_system(cs_sles_find_or_add(field_id, NULL),
                    matrix,
                    field_id,
                    eqp,
                    rhs_norm,
                    fld->val, rhs);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         Theta time scheme is used to progress in time.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_theta(const cs_mesh_t            *mesh,
                            const int                   field_id,
                            const cs_equation_param_t  *eqp,
                            cs_equation_builder_t      *eqb,
                            void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  inv_dtcur = 1./dt_cur;
  const double  tcoef = 1 - eqp->theta;

  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_THETA ||
         eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO);

  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  cs_real_t  rhs_norm = 0.;

  /* Build an array storing the Dirichlet values at vertices
     and another one with a tags to detect vertices related to a
     Neumann BC */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *forced_ids = NULL;

  _vbs_setup(t_cur + dt_cur, mesh, eqp, eqb, eqc->vtx_bc_flag,
             &dir_values, &forced_ids);

  /* Initialize the local system: rhs */
  cs_real_t  *rhs = NULL;
  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Detect the first call (in this case, we compute the initial source term)*/
  _Bool  compute_initial_source = false;
  if (ts->nt_cur == ts->nt_prev || ts->nt_prev == 0) {

    compute_initial_source = true;

  }
  else { /* Add contribution of the previous computed source term */

    if (eqc->source_terms != NULL) {

      assert(cs_equation_param_has_sourceterm(eqp));
      for (cs_lnum_t v = 0; v < n_vertices; v++)
        rhs[v] += tcoef * eqc->source_terms[v];
      memset(eqc->source_terms, 0, n_vertices * sizeof(cs_real_t));

      if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC ||
          eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED) {

        assert(eqc->vtx_bc_flag != NULL);
        for (cs_lnum_t v = 0; v < n_vertices; v++) {
          if (cs_cdo_bc_is_dirichlet(eqc->vtx_bc_flag[v]))
            rhs[v] = 0.;
        }

      } /* Algebraic or penalized enforcement is set */
    } /* At least one source term is defined */

  }

  /* Initialize the local system: matrix */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_equation_get_mav(matrix, eqp->omp_assembly_choice, 1);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav,               \
         dir_values, fld, forced_ids, rs, compute_initial_source,       \
         _vbs_cell_system, _vbs_cell_builder, rhs_norm)
  {
    /* Set variables and structures inside the OMP section so that each thread
       has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_real_t  time_eval = t_cur + 0.5*dt_cur;

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_matrix_assembler_buf_t  *mab = cs_equation_get_assembly_buffers(t_id);
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = _vbs_cell_system[t_id];
    cs_cell_builder_t  *cb = _vbs_cell_builder[t_id];

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _vbs_init_cell_system(time_eval, cell_flag, cm, eqp, eqb,
                            dir_values, eqc->vtx_bc_flag, forced_ids, fld->val,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vbs_advection_diffusion_reaction(time_eval,
                                        eqp, eqb, eqc, cm, fm, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */
        if (compute_initial_source) {

          /* Reset the local contribution */
          memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

          cs_source_term_compute_cellwise(eqp->n_source_terms,
                      (cs_xdef_t *const *)eqp->source_terms,
                                          cm,
                                          eqb->source_mask,
                                          eqb->compute_source,
                                          t_cur,
                                          NULL,  /* No input structure */
                                          cb,    /* mass matrix is cb->hdg */
                                          csys->source);

          for (short int v = 0; v < cm->n_vc; v++)
            csys->rhs[v] += tcoef * csys->source[v];

        }

        /* Reset the local contribution */
        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */
        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        t_cur + dt_cur,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        for (short int v = 0; v < cm->n_vc; v++)
          csys->rhs[v] += eqp->theta * csys->source[v];

      } /* End of term source */

      /* Apply boundary conditions (those which are weakly enforced) */
      _vbs_apply_weak_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      /* STEP.1 >> Compute the contribution of the "adr" to the RHS:
       *           tcoef*adr_pn where adr_pn = csys->mat * p_n */
      double  *adr_pn = cb->values;
      cs_sdm_square_matvec(csys->mat, csys->val_n, adr_pn);
      for (short int i = 0; i < csys->n_dofs; i++) /* n_dofs = n_vc */
        csys->rhs[i] -= tcoef * adr_pn[i];

      /* STEP.2 >> Multiply csys->mat by theta */
      for (int i = 0; i < csys->n_dofs*csys->n_dofs; i++)
        csys->mat->val[i] *= eqp->theta;

      /* STEP.3 >> Handle the mass matrix
       * Two contributions for the mass matrix
       *  a) add to csys->mat
       *  b) add to rhs mass_mat * p_n */
      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        /* |c|*wvc = |dual_cell(v) cap c| */
        assert(cs_flag_test(eqb->msh_flag, CS_CDO_LOCAL_PVQ));
        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
         *       >> Update the cellwise system with the time matrix */
        for (short int i = 0; i < cm->n_vc; i++) {

          const double  dval = ptyc * cm->wvc[i];

          /* Update the RHS with mass_mat * values at time t_n */
          csys->rhs[i] += dval * csys->val_n[i];

          /* Add the diagonal contribution from time matrix to the local
             system */
          csys->mat->val[i*(cm->n_vc + 1)] += dval;

        }

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val * inv_dtcur;
        const cs_sdm_t  *mass_mat = cb->hdg;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix */

        /* Update rhs with mass_mat*p^n */
        double  *time_pn = cb->values;
        cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);
        for (short int i = 0; i < csys->n_dofs; i++)
          csys->rhs[i] += tpty_coef*time_pn[i];

        /* Update the cellwise system with the time matrix */
        cs_sdm_add_mult(csys->mat, tpty_coef, mass_mat);

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after adding time", csys);
#endif

      /* Compute a norm of the RHS for the normalization of the SLES */
      _vbs_compute_cw_sles_normalization(eqp, cm, csys, &rhs_norm);

      /* Enforce values if needed (internal or Dirichlet) */
      _vbs_enforce_values(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      eqc->assemble(csys, rs, mab, mav);         /* Matrix assembly */

      for (short int v = 0; v < cm->n_vc; v++)   /* RHS assembly */
#       pragma omp atomic
        rhs[cm->v_ids[v]] += csys->rhs[v];

      if (eqc->source_terms != NULL) {
        for (short int v = 0; v < cm->n_vc; v++) /* Source term assembly */
#         pragma omp atomic
          eqc->source_terms[cm->v_ids[v]] += csys->source[v];
      }

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(forced_ids);
  cs_matrix_assembler_values_finalize(&mav);

  /* Last step in the computation of the renormalization coefficient */
  _vbs_sync_sles_normalization(eqp, &rhs_norm);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Now solve the system */
  _vbs_solve_system(cs_sles_find_or_add(field_id, NULL),
                    matrix,
                    field_id,
                    eqp,
                    rhs_norm,
                    fld->val, rhs);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh vertices for the variable field
 *         associated to the given context
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_scaleq_get_vertex_values(void      *context)
{
  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  return pot->val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an array of values at mesh cells by interpolating the
 *         variable field associated to the given context located at mesh
 *         vertices
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_scaleq_get_cell_values(void      *context)
{
  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Reset buffer of values */
  if (eqc->cell_values == NULL)
    BFT_MALLOC(eqc->cell_values, quant->n_cells, cs_real_t);
  memset(eqc->cell_values, 0, quant->n_cells*sizeof(cs_real_t));

  /* Compute the values at cell centers from an interpolation of the field
     values defined at vertices */
  cs_reco_pv_at_cell_centers(connect->c2v, quant, pot->val,
                             eqc->cell_values);

  return eqc->cell_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain between time t_cur and t_cur + dt_cur
 *         Case of scalar-valued CDO vertex-based scheme
 *
 * \param[in]      eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb      pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context  pointer to a scheme builder structure
 *
 * \return a pointer to a \ref cs_equation_balance_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_balance_t *
cs_cdovb_scaleq_balance(const cs_equation_param_t     *eqp,
                        cs_equation_builder_t         *eqb,
                        void                          *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  dt_cur = cs_shared_time_step->dt[0];
  const cs_real_t  time_eval = t_cur + 0.5*dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  /* Allocate and initialize the structure storing the balance evaluation */
  cs_equation_balance_t  *eb = cs_equation_balance_create(cs_flag_primal_vtx,
                                                          quant->n_vertices);

  /* OpenMP block */
#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, eqc, pot, eb, _vbs_cell_builder)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Set inside the OMP section so that each thread has its own value
     * Each thread get back its related structures:
     * Get the cell-wise view of the mesh and the algebraic system */
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = _vbs_cell_builder[t_id];

    /* Set inside the OMP section so that each thread has its own value */
    cs_real_t  _p_cur[10], _p_prev[10], _p_theta[10];
    cs_real_t  *p_cur = NULL, *p_prev = NULL, *p_theta = NULL;

    if (connect->n_max_vbyc > 10) {
      BFT_MALLOC(p_cur, connect->n_max_vbyc, cs_real_t);
      BFT_MALLOC(p_prev, connect->n_max_vbyc, cs_real_t);
      BFT_MALLOC(p_theta, connect->n_max_vbyc, cs_real_t);
    }
    else {
      p_cur = _p_cur;
      p_prev = _p_prev;
      p_theta = _p_theta;
    }

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag = cs_equation_cell_mesh_flag(cell_flag, eqb);

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Initialize the properties at a cellwise level */
      cs_equation_init_properties_cw(eqp, eqb, time_eval, cell_flag, cm, cb);

      /* Set the value of the current potential */
      for (short int v = 0; v < cm->n_vc; v++)
        p_cur[v] = pot->val[cm->v_ids[v]];

      if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX)
        eqc->get_mass_matrix(eqc->hdg_mass, cm, cb); /* stored in cb->hdg */

      /* Unsteady term */
      if (cs_equation_param_has_time(eqp)) {

        /* Set the value of the previous potential */
        for (short int v = 0; v < cm->n_vc; v++)
          p_prev[v] = pot->val_pre[cm->v_ids[v]];

        if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) {

          assert(cs_flag_test(eqb->msh_flag, CS_CDO_LOCAL_PVQ));
          /* |c|*wvc = |dual_cell(v) cap c| */
          const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;
          for (short int v = 0; v < cm->n_vc; v++) {
            cs_real_t  dp = p_cur[v] - p_prev[v];
#           pragma omp atomic
            eb->unsteady_term[cm->v_ids[v]] += ptyc * cm->wvc[v] * dp;
          }

        }
        else {

          const double  ptyc = cb->tpty_val * inv_dtcur;
          cs_real_t  *dp = cb->values;
          cs_real_t  *res = cb->values + cm->n_vc;
          for (short int v = 0; v < cm->n_vc; v++) {
            res[v] = 0.;
            dp[v] = p_cur[v] - p_prev[v];
          }
          cs_sdm_square_matvec(cb->hdg, dp, res);

          for (short int v = 0; v < cm->n_vc; v++) {
#           pragma omp atomic
            eb->unsteady_term[cm->v_ids[v]] += ptyc*res[v];
          }

        } /* Add unsteady contribution */

      } /* TIME */

      /* Set p_theta */
      switch (eqp->time_scheme) {
      case CS_TIME_SCHEME_EULER_EXPLICIT:
        for (short int v = 0; v < cm->n_vc; v++)
          p_theta[v] = p_prev[v];
        break;

      case CS_TIME_SCHEME_CRANKNICO:
        for (short int v = 0; v < cm->n_vc; v++)
          p_theta[v] = 0.5*(p_cur[v] + p_prev[v]);
        break;

      case CS_TIME_SCHEME_THETA:
        for (short int v = 0; v < cm->n_vc; v++)
          p_theta[v] = eqp->theta*p_cur[v] + (1-eqp->theta)*p_prev[v];
        break;

      default:
        for (short int v = 0; v < cm->n_vc; v++)
          p_theta[v] = p_cur[v];
        break;

      } /* Switch on time scheme */

      /* Reaction term */
      if (cs_equation_param_has_reaction(eqp)) {

        /* Define the local reaction property */
        const double  rpty_val = cb->rpty_val;

        cs_real_t  *res = cb->values;
        memset(res, 0, cm->n_vc*sizeof(cs_real_t));
        cs_sdm_square_matvec(cb->hdg, p_theta, res);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->reaction_term[cm->v_ids[v]] += rpty_val * res[v];
        }

      } /* Reaction */

      /* Diffusion */
      if (cs_equation_param_has_diffusion(eqp)) {

        /* local matrix owned by the cellwise builder (store in cb->loc) */
        eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        cs_real_t  *res = cb->values;
        memset(res, 0, cm->n_vc*sizeof(cs_real_t));
        cs_sdm_square_matvec(cb->loc, p_theta, res);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->diffusion_term[cm->v_ids[v]] += res[v];
        }

      } /* Diffusion */

      /* Advection term */
      if (cs_equation_param_has_convection(eqp)) {

        /* Define the local advection matrix */
        eqc->get_advection_matrix(eqp, cm, time_eval, fm, cb);

        cs_real_t  *res = cb->values;
        memset(res, 0, cm->n_vc*sizeof(cs_real_t));
        cs_sdm_square_matvec(cb->loc, p_theta, res);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->advection_term[cm->v_ids[v]] += res[v];
        }

      } /* End of advection */

      /* Source term */
      if (cs_equation_param_has_sourceterm(eqp)) {

        cs_real_t  *src = cb->values;
        memset(src, 0, cm->n_vc*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */
        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        time_eval,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        src);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->source_term[cm->v_ids[v]] += src[v];
        }

      } /* End of term source */

      /* Boundary conditions */
      if (cell_flag &  CS_FLAG_BOUNDARY_CELL_BY_FACE) {

        const cs_cdo_bc_face_t  *face_bc = eqb->face_bc;

        /* Identify which face is a boundary face */
        for (short int f = 0; f < cm->n_fc; f++) {
          const cs_lnum_t  bf_id = cm->f_ids[f] - quant->n_i_faces;
          if (bf_id > -1) { /* Border face */

            /* Advective flux */
            if (cs_equation_param_has_convection(eqp)) {
              cs_advection_field_cw_boundary_f2v_flux(cm,
                                                      eqp->adv_field,
                                                      f,
                                                      time_eval,
                                                      cb->values);

              for (short int v = 0; v < cm->n_vc; v++) {
                const cs_real_t  adv_flux = cb->values[v] * p_cur[v];
#               pragma omp atomic
                eb->boundary_term[bf_id] += adv_flux;
#               pragma omp atomic
                eb->advection_term[cm->v_ids[v]] += adv_flux;
              }
            }

            /* Diffusive flux */
            if (cs_equation_param_has_diffusion(eqp)) {
              if (cs_cdo_bc_is_dirichlet(face_bc->flag[bf_id])) {
                cs_cdovb_diffusion_p0_face_flux(f, cm,
                        (const cs_real_t (*)[3])cb->dpty_mat,
                                                p_cur,
                                                cb->values);

                for (short int v = 0; v < cm->n_vc; v++) {
#                 pragma omp atomic
                  eb->boundary_term[bf_id] += cb->values[v];
#                 pragma omp atomic
                  eb->diffusion_term[cm->v_ids[v]] += cb->values[v];
                }

              }
            }

          } /* Is a boundary face */
        } /* Loop on cell faces */

      } /* Boundary conditions */

    } /* Main loop on cells */

    if (p_cur != _p_cur) {
      BFT_FREE(p_cur);
      BFT_FREE(p_prev);
      BFT_FREE(p_theta);
    }

  } /* OpenMP Block */

  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++)
    eb->balance[v_id] =
      eb->unsteady_term[v_id] + eb->reaction_term[v_id] +
      eb->diffusion_term[v_id] + eb->advection_term[v_id] +
      eb->source_term[v_id];

  /* Parallel synchronisation */
  cs_equation_balance_sync(connect, eb);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);

  return eb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each vertex of a boundary face, the portion of diffusive
 *         flux across the boundary face. The surface attached to each vertex
 *         corresponds to the intersection of its dual cell (associated to
 *         a vertex of the face) with the face.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]       t_eval     time at which one performs the evaluation
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in]       pdi        pointer to an array of field values
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  vf_flux    pointer to the values of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_boundary_diff_flux(const cs_real_t              t_eval,
                                   const cs_equation_param_t   *eqp,
                                   const cs_real_t             *pdi,
                                   cs_equation_builder_t       *eqb,
                                   cs_real_t                   *vf_flux)
{
  if (vf_flux == NULL)
    return;

  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(vf_flux, 0, connect->bf2v->idx[quant->n_b_faces]*sizeof(cs_real_t));

    cs_timer_t  t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
    return;
  }

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, vf_flux, pdi, _vbs_cell_builder)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_cdo_bc_face_t  *face_bc = eqb->face_bc;
    const cs_adjacency_t  *bf2v = connect->bf2v;
    const cs_adjacency_t  *f2c = connect->f2c;
    const cs_lnum_t  fidx_shift = f2c->idx[quant->n_i_faces];

    cs_real_t  *pot = NULL, *flux = NULL;
    BFT_MALLOC(pot, connect->n_max_vbyc + 1, cs_real_t); /* +1 for WBS */
    BFT_MALLOC(flux, connect->n_max_vbyc , cs_real_t);

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_builder_t  *cb = _vbs_cell_builder[t_id];
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);

    /* msh_flag for Neumann and Robin BCs. Add add_flag for the other cases
       when one has to reconstruct a flux */
    cs_flag_t  msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_FV;
    cs_flag_t  add_flag = CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE | CS_CDO_LOCAL_PEQ |
      CS_CDO_LOCAL_PFQ;

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
    case CS_PARAM_HODGE_ALGO_VORONOI:
      add_flag |= CS_CDO_LOCAL_DFQ;
      break;

    case CS_PARAM_HODGE_ALGO_WBS:
      add_flag |= CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_FEQ;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid Hodge algorithm", __func__);

    } /* Switch hodge algo. */

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_equation_set_diffusion_property(eqp, 0, t_eval, 0, cb);

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t bf_id = 0; bf_id < quant->n_b_faces; bf_id++) {

      const cs_lnum_t  f_id = bf_id + quant->n_i_faces;
      const cs_lnum_t  c_id = f2c->ids[bf_id + fidx_shift];
      const cs_lnum_t  *idx  = bf2v->idx + bf_id;

      cs_real_t  *_flx = vf_flux + idx[0];

      switch (face_bc->flag[bf_id]) {

      case CS_CDO_BC_HMG_NEUMANN:
        memset(_flx, 0, (idx[1]-idx[0])*sizeof(cs_real_t));
        break;

      case CS_CDO_BC_NEUMANN:
        {
          cs_real_t  *neu_values = cb->values;

          /* Set the local mesh structure for the current cell */
          cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

          const short int  f = cs_cell_mesh_get_f(f_id, cm);

          cs_equation_compute_neumann_sv(t_eval,
                                         face_bc->def_ids[bf_id],
                                         f,
                                         quant,
                                         eqp,
                                         cm,
                                         neu_values);

          short int n_vf = 0;
          for (int i = cm->f2v_idx[f]; i < cm->f2v_idx[f+1]; i++)
            _flx[n_vf++] = neu_values[cm->f2v_ids[i]];

        }
        break;

      case CS_CDO_BC_ROBIN:
        {
          cs_real_t  *robin_values = cb->values;
          cs_real_t  *wvf = cb->values + 3;

          /* Set the local mesh structure for the current cell */
          cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

          const short int  f = cs_cell_mesh_get_f(f_id, cm);
          const cs_real_t  f_area = quant->b_face_surf[bf_id];

          /* Robin BC expression: K du/dn + alpha*(p - p0) = g */
          cs_equation_compute_robin(t_eval,
                                    face_bc->def_ids[bf_id],
                                    f,
                                    quant,
                                    eqp,
                                    cm,
                                    robin_values);

          const cs_real_t  alpha = robin_values[0];
          const cs_real_t  p0 = robin_values[1];
          const cs_real_t  g = robin_values[2];

          cs_cdo_quantities_compute_b_wvf(connect, quant, bf_id, wvf);

          short int n_vf = 0;
          for (int i = cm->f2v_idx[f]; i < cm->f2v_idx[f+1]; i++) {
            const cs_real_t  pv = pdi[cm->v_ids[cm->f2v_ids[i]]];
            _flx[n_vf] = f_area*wvf[n_vf]*(alpha*(p0 - pv) + g);
            n_vf++;
          }
        }
        break;

      default:
        { /* Reconstruct a normal flux at the boundary face */

          /* Set the local mesh structure for the current cell */
          cs_cell_mesh_build(c_id, msh_flag | add_flag, connect, quant, cm);

          const short int  f = cs_cell_mesh_get_f(f_id, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
          if (cs_dbg_cw_test(eqp, cm, NULL)) cs_cell_mesh_dump(cm);
#endif
          if (!eqb->diff_pty_uniform)
            cs_property_tensor_in_cell(cm,
                                       eqp->diffusion_property,
                                       t_eval,
                                       eqp->diffusion_hodge.inv_pty,
                                       cb->dpty_mat);

          /* Define a local buffer keeping the value of the discrete potential
             for the current cell */
          for (short int v = 0; v < cm->n_vc; v++)
            pot[v] = pdi[cm->v_ids[v]];

          if (eqp->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_WBS) {

            /* Interpolate also the value of the potential at the cell center */
            pot[cm->n_vc] = cs_reco_cw_scalar_pv_at_cell_center(cm, pot);

            cs_cdo_diffusion_wbs_vbyf_flux(f, eqp, cm, pot, cb, flux);

          }
          else
            cs_cdo_diffusion_svb_cost_vbyf_flux(f, eqp, cm, pot, cb, flux);

          /* Fill the global flux array */
          short int n_vf = 0;
          for (int i = cm->f2v_idx[f]; i < cm->f2v_idx[f+1]; i++)
            _flx[n_vf++] = flux[cm->f2v_ids[i]];

        }
        break;

      } /* End of switch */

    } /* End of loop on boundary faces */

    BFT_FREE(pot);
    BFT_FREE(flux);

  } /* End of Open block */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]       normal     indicate in which direction flux is > 0
 * \param[in]       pdi        pointer to an array of field values
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to data specific for this scheme
 * \param[in, out]  d_flux     pointer to the value of the diffusive flux
 * \param[in, out]  c_flux     pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_flux_across_plane(const cs_real_t             normal[],
                                  const cs_real_t            *pdi,
                                  const cs_equation_param_t  *eqp,
                                  int                         ml_id,
                                  cs_equation_builder_t      *eqb,
                                  void                       *context,
                                  double                     *d_flux,
                                  double                     *c_flux)
{
  CS_UNUSED(context);

  *d_flux = 0.,  *c_flux = 0.;

  if (pdi == NULL)
    return;

  cs_mesh_location_type_t  ml_t = cs_mesh_location_get_type(ml_id);

  if (ml_t != CS_MESH_LOCATION_INTERIOR_FACES &&
      ml_t != CS_MESH_LOCATION_BOUNDARY_FACES) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" Mesh location type is incompatible with the computation\n"
                    " of the flux across faces.\n"));
    return;
  }

  const cs_timer_t  t0 = cs_timer_time();
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_adjacency_t  *f2c = connect->f2c;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;

  double  pf;
  cs_real_3_t  gc, pty_gc;
  cs_real_33_t  pty_tens;

  if (ml_t == CS_MESH_LOCATION_BOUNDARY_FACES) { /* Belongs to only one cell */

    const cs_lnum_t  n_i_faces = connect->n_faces[2];
    const cs_lnum_t  *cell_ids = f2c->ids + f2c->idx[n_i_faces];

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  bf_id = (elt_ids != NULL) ? elt_ids[i] : i;
      const cs_lnum_t  f_id = n_i_faces + bf_id;
      const cs_lnum_t  c_id = cell_ids[bf_id];
      const cs_quant_t  f = cs_quant_set_face(f_id, quant);
      const short int  sgn = (_dp3(f.unitv, normal) < 0) ? -1 : 1;
      const double  coef = sgn * f.meas;

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Compute the local diffusive flux */
        cs_reco_grad_cell_from_pv(c_id, connect, quant, pdi, gc);
        cs_property_get_cell_tensor(c_id, t_cur,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);
        cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);

        /* Update the diffusive flux */
        *d_flux += -coef * _dp3(f.unitv, pty_gc);

      }

      if (cs_equation_param_has_convection(eqp)) {

        cs_nvec3_t  adv_c;

        /* Compute the local advective flux */
        cs_advection_field_get_cell_vector(c_id, eqp->adv_field, &adv_c);
        cs_reco_pf_from_pv(f_id, connect, quant, pdi, &pf);

        /* Update the convective flux */
        *c_flux += coef * adv_c.meas * _dp3(adv_c.unitv, f.unitv) * pf;

      }

    } /* Loop on selected border faces */

  }
  else if (ml_t == CS_MESH_LOCATION_INTERIOR_FACES) {

    if (n_elts[0] > 0 && elt_ids == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _(" Computing the flux across all interior faces is"
                  " not managed yet."));

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  f_id = elt_ids[i];
      const cs_quant_t  f = cs_quant_set_face(f_id, quant);
      const short int  sgn = (_dp3(f.unitv, normal) < 0) ? -1 : 1;

      for (cs_lnum_t j = f2c->idx[f_id]; j < f2c->idx[f_id+1]; j++) {

        const cs_lnum_t  c_id = f2c->ids[j];

        if (cs_equation_param_has_diffusion(eqp)) {

          const double  coef = 0.5 * sgn * f.meas; /* mean value at the face */

          /* Compute the local diffusive flux */
          cs_reco_grad_cell_from_pv(c_id, connect, quant, pdi, gc);
          cs_property_get_cell_tensor(c_id, t_cur,
                                      eqp->diffusion_property,
                                      eqp->diffusion_hodge.inv_pty,
                                      pty_tens);
          cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);
          *d_flux += -coef * _dp3(f.unitv, pty_gc);

        }

        if (cs_equation_param_has_convection(eqp)) {

          cs_nvec3_t  adv_c;

          /* Compute the local advective flux */
          cs_reco_pf_from_pv(f_id, connect, quant, pdi, &pf);

          /* Evaluate the advection field at the face */
          cs_advection_field_get_cell_vector(c_id, eqp->adv_field, &adv_c);

          /* Update the convective flux (upwinding according to adv_f) */
          const double  dpc = _dp3(adv_c.unitv, f.unitv);

          cs_real_t  fconv_flux = 0;
          if (dpc > 0) {
            if (f2c->sgn[j] > 0) /* nf points outward c; adv.nf > 0 */
              fconv_flux = adv_c.meas * dpc * sgn * f.meas * pf;
          }
          else if (dpc < 0) {
            if (f2c->sgn[j] < 0) /* nf points inward c; adv.nf < 0 */
              fconv_flux = adv_c.meas * dpc * sgn * f.meas * pf;
          }
          else  /* centered approach */
            fconv_flux = 0.5 * adv_c.meas * dpc * sgn * f.meas * pf;

          *c_flux += fconv_flux;

        }

      }

    }  /* Loop on selected interior faces */

  } /* Set of interior or border faces */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of an approximation of a constant diffusive
 *         flux (a vector) in each cell.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to cs_cdovb_scaleq_t structure
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_diff_flux_in_cells(const cs_real_t             *values,
                                   const cs_equation_param_t   *eqp,
                                   cs_real_t                    t_eval,
                                   cs_equation_builder_t       *eqb,
                                   void                        *context,
                                   cs_real_t                   *diff_flux)
{
  CS_UNUSED(context);

  if (diff_flux == NULL)
    return ;

  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* If no diffusion, return after resetting */
  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(diff_flux, 0, 3*quant->n_cells*sizeof(cs_real_t));
    return;
  }

  cs_timer_t  t0 = cs_timer_time();

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)        \
  shared(t_eval, quant, connect, eqp, eqb, diff_flux, values, _vbs_cell_builder)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    double  *pot = NULL;
    BFT_MALLOC(pot, connect->n_max_vbyc + 1, double); /* +1 for WBS */

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cdo_diffusion_cw_flux_t  *compute_flux = NULL;
    cs_cell_builder_t  *cb = _vbs_cell_builder[t_id];
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_flag_t  msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_EV;

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      msh_flag |= CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_PVQ;
      compute_flux = cs_cdo_diffusion_svb_cost_get_cell_flux;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      msh_flag |= CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_PVQ;
      compute_flux = cs_cdo_diffusion_svb_cost_get_cell_flux;
      break;

    case CS_PARAM_HODGE_ALGO_WBS:
      msh_flag |= CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
        CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
      compute_flux = cs_cdo_diffusion_wbs_get_cell_flux;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid Hodge algorithm", __func__);

    } /* Switch hodge algo. */

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_equation_set_diffusion_property(eqp, 0, t_eval, 0, cb);

    /* Define the flux by cellwise contributions */
#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 3
      if (cs_dbg_cw_test(eqp, cm, NULL)) cs_cell_mesh_dump(cm);
#endif

      if (!eqb->diff_pty_uniform) {
        cs_property_tensor_in_cell(cm,
                                   eqp->diffusion_property,
                                   t_eval,
                                   eqp->diffusion_hodge.inv_pty,
                                   cb->dpty_mat);
        if (eqp->diffusion_hodge.is_iso)
          cb->dpty_val = cb->dpty_mat[0][0];
      }

      /* Define a local buffer keeping the value of the discrete potential
         for the current cell */
      for (short int v = 0; v < cm->n_vc; v++)
        pot[v] = values[cm->v_ids[v]];

      if (eqp->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_WBS) {
        /* Interpolate also the value of the potential at the cell center */
        pot[cm->n_vc] = 0.;
        for (short int v = 0; v < cm->n_vc; v++)
          pot[cm->n_vc] += cm->wvc[v]*pot[v];
      }

      compute_flux(cm, pot, cb, diff_flux + 3*c_id);

    } /* Loop on cells */

    BFT_FREE(pot);

  } /* OMP Section */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux accross dual faces
 *         (a scalar) in each cell.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to cs_cdovb_scaleq_t structure
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_diff_flux_dfaces(const cs_real_t             *values,
                                 const cs_equation_param_t   *eqp,
                                 cs_real_t                    t_eval,
                                 cs_equation_builder_t       *eqb,
                                 void                        *context,
                                 cs_real_t                   *diff_flux)
{
  CS_UNUSED(context);

  if (diff_flux == NULL)
    return ;

  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* If no diffusion, return after resetting */
  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(diff_flux, 0, connect->c2e->idx[quant->n_cells]*sizeof(cs_real_t));
    return;
  }

  cs_timer_t  t0 = cs_timer_time();

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(t_eval, quant, connect, eqp, eqb, diff_flux, values,           \
         _vbs_cell_builder)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif
    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_hodge_t  *get_diffusion_hodge = NULL;
    cs_cdo_diffusion_cw_flux_t  *compute_flux = NULL;
    cs_cell_builder_t  *cb = _vbs_cell_builder[t_id];
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);

    double  *pot = NULL;
    BFT_MALLOC(pot, connect->n_max_vbyc + 1, double); /* +1for WBS algo. */

    cs_flag_t  msh_flag = CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PV |
      CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_DFQ;

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      get_diffusion_hodge = cs_hodge_epfd_cost_get;
      compute_flux = cs_cdo_diffusion_svb_cost_get_dfbyc_flux;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      msh_flag |= CS_CDO_LOCAL_EFQ;
      get_diffusion_hodge = cs_hodge_epfd_voro_get;
      compute_flux = cs_cdo_diffusion_svb_cost_get_dfbyc_flux;
      break;

    case CS_PARAM_HODGE_ALGO_WBS:
      msh_flag |= CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
        CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EFQ;
      compute_flux = cs_cdo_diffusion_wbs_get_dfbyc_flux;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid Hodge algorithm", __func__);
      break;

    } /* Switch hodge algo. */

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_equation_set_diffusion_property(eqp, 0, t_eval, 0, cb);

    /* Define the flux by cellwise contributions */
#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
      if (cs_dbg_cw_test(eqp, cm, NULL)) cs_cell_mesh_dump(cm);
#endif

      if (!eqb->diff_pty_uniform) {
        cs_property_tensor_in_cell(cm,
                                   eqp->diffusion_property,
                                   t_eval,
                                   eqp->diffusion_hodge.inv_pty,
                                   cb->dpty_mat);
        if (eqp->diffusion_hodge.is_iso)
          cb->dpty_val = cb->dpty_mat[0][0];
      }

      /* Define a local buffer keeping the value of the discrete potential
         for the current cell */
      for (short int v = 0; v < cm->n_vc; v++)
        pot[v] = values[cm->v_ids[v]];

      switch (eqp->diffusion_hodge.algo) {

      case CS_PARAM_HODGE_ALGO_COST:
      case CS_PARAM_HODGE_ALGO_VORONOI:
        /* Compute the local Hodge operator and store it in cb->hdg */
        get_diffusion_hodge(eqp->diffusion_hodge, cm, cb);
        break;

      case CS_PARAM_HODGE_ALGO_WBS:
        pot[cm->n_vc] = 0.;
        for (short int v = 0; v < cm->n_vc; v++)
          pot[cm->n_vc] += cm->wvc[v]*pot[v];
        break;

      default:
        break; /* Nothing else to do */

      } /* End of switch */

      /* Compute and store the flux */
      compute_flux(cm, pot, cb, diff_flux + connect->c2e->idx[c_id]);

    } /* Loop on cells */

    BFT_FREE(pot);

  } /* OMP Section */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *context)
{
  CS_UNUSED(field);
  CS_UNUSED(context);

  const cs_timer_t  t0 = cs_timer_time();

  if (cs_equation_param_has_convection(eqp)) {
    if (eqp->process_flag & CS_EQUATION_POST_UPWIND_COEF) {

      int  len = strlen(eqname) + 8 + 1;
      char *postlabel = NULL;
      BFT_MALLOC(postlabel, len, char);
      sprintf(postlabel, "%s.UpwCoef", eqname);

      /* Compute in each cell an evaluation of upwind weight value */
      cs_real_t  *work_c = cs_equation_get_tmpbuf();
      cs_cdo_advection_cell_upwind_coef(cs_shared_quant,
                                        eqp->adv_scheme,
                                        work_c);

      cs_post_write_var(CS_POST_MESH_VOLUME,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        postlabel,
                        1,
                        true,                 /* interlace */
                        true,                 /* true = original mesh */
                        CS_POST_TYPE_cs_real_t,
                        work_c,               /* values on cells */
                        NULL,                 /* values at internal faces */
                        NULL,                 /* values at border faces */
                        cs_shared_time_step); /* time step management struct. */

      BFT_FREE(postlabel);

    }
  } /* Post the upwind coefficient attached to cells */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
