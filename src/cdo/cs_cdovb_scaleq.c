/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
 * reaction of scalar-valued equations with source terms
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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
#include "cs_cdovb_priv.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_evaluate.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
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
static cs_cell_sys_t      **_svb_cell_system = NULL;
static cs_cell_builder_t  **_svb_cell_builder = NULL;

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
_svb_create_cell_builder(const cs_cdo_connect_t   *connect)
{
  const int  n_vc = connect->n_max_vbyc;
  const int  n_ec = connect->n_max_ebyc;
  assert(n_ec > n_vc);
  cs_cell_builder_t *cb = cs_cell_builder_create();

  BFT_MALLOC(cb->ids, n_ec, int);
  memset(cb->ids, 0, n_ec*sizeof(int));

  int  size = n_ec*(n_ec+1);
  size = CS_MAX(4*n_ec + 3*n_vc, size);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(cs_real_t));

  size = 2*n_ec;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of
     operators */
  cb->loc = cs_sdm_square_create(n_vc);
  cb->aux = cs_sdm_square_create(n_ec);

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
 * \param[in, out] vtx_bc_flag     pointer to an array of BC flag for each vtx
 * \param[in, out] p_dir_values    pointer to the Dirichlet values to set
 * \param[in, out] p_enforced_ids  pointer to the list of enforced vertices
 */
/*----------------------------------------------------------------------------*/

static void
_svb_setup(cs_real_t                      t_eval,
           const cs_mesh_t               *mesh,
           const cs_equation_param_t     *eqp,
           const cs_equation_builder_t   *eqb,
           cs_flag_t                      vtx_bc_flag[],
           cs_real_t                     *p_dir_values[],
           cs_lnum_t                     *p_enforced_ids[])
{
  assert(vtx_bc_flag != NULL);  /* Sanity check */
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_real_t  *dir_values = NULL;

  /* Compute the values of the Dirichlet BC */
  BFT_MALLOC(dir_values, quant->n_vertices, cs_real_t);

  /* Compute the values of the Dirichlet BC */
  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   _svb_cell_builder[0], /* static variable */
                                   vtx_bc_flag,
                                   dir_values);

  *p_dir_values = dir_values;

  /* Internal enforcement of DoFs  */
  if (cs_equation_param_has_internal_enforcement(eqp))
    cs_equation_build_dof_enforcement(quant->n_vertices,
                                      connect->c2v,
                                      eqp,
                                      p_enforced_ids);
  else
    *p_enforced_ids = NULL;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the local structure for the current cell
 *         Case of scalar-valued CDO-Vb schemes
 *
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
_svb_init_cell_system(const cs_cell_mesh_t          *cm,
                      const cs_equation_param_t     *eqp,
                      const cs_equation_builder_t   *eqb,
                      const cs_real_t                dir_values[],
                      const cs_flag_t                vtx_bc_flag[],
                      const cs_lnum_t                forced_ids[],
                      const cs_real_t                field_tn[],
                      cs_cell_sys_t                 *csys,
                      cs_cell_builder_t             *cb)
{
  /* Cell-wise view of the linear system to build */
  csys->c_id = cm->c_id;
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
  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Set the bc (specific part) */
    cs_equation_vb_set_cell_bc(cm,
                               eqp,
                               eqb->face_bc,
                               vtx_bc_flag,
                               dir_values,
                               cb->t_bc_eval,
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
  if (cb->cell_flag == CS_FLAG_BOUNDARY_CELL_BY_VERTEX) {

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the local matrices arising from the diffusion, advection,
 *         reaction terms.
 *         mass_hodge could be set to NULL if a Voronoi algo. is used.
 *         Otherwise, the mass matrix is computed.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] mass_hodge  pointer to a cs_hodge_t structure (mass matrix)
 * \param[in, out] diff_hodge  pointer to a cs_hodge_t structure (diffusion)
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_svb_conv_diff_reac(const cs_equation_param_t     *eqp,
                    const cs_equation_builder_t   *eqb,
                    const cs_cdovb_scaleq_t       *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_face_mesh_t                *fm,
                    cs_hodge_t                    *mass_hodge,
                    cs_hodge_t                    *diff_hodge,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                  * =========== */
    assert(mass_hodge != NULL);

    /* Build the mass matrix and store it in mass_hodge->matrix */
    eqc->get_mass_matrix(cm, mass_hodge, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys)) {
      cs_log_printf(CS_LOG_DEFAULT, ">> Celll mass matrix");
      cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids,
                  mass_hodge->matrix);
    }
#endif
  }

  if (cs_equation_param_has_diffusion(eqp)) {   /* DIFFUSION TERM
                                                 * ============== */
    assert(diff_hodge != NULL);

    /* Set the diffusion property */
    if (!(eqb->diff_pty_uniform))
      cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, cb->cell_flag,
                                     diff_hodge);

    /* Define the local stiffness matrix: local matrix owned by the cellwise
       builder (store in cb->loc) */
    eqc->get_stiffness_matrix(cm, diff_hodge, cb);

    /* Add the local diffusion operator to the local system */
    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after diffusion", csys);
#endif
  }

  if (cs_equation_param_has_convection(eqp) &&
      ((cb->cell_flag & CS_FLAG_SOLID_CELL) == 0)) {  /* ADVECTION TERM
                                                       * ============== */

    /* Define the local advection matrix (the diffusion property
     * is given as parameter since some schemes introduce a portion of
     * upwinding weigthed with respect to the ratio advection/diffusion
     */
    cs_property_data_t  *diff_pty = (diff_hodge == NULL) ?
      NULL : diff_hodge->pty_data;
    eqc->get_advection_matrix(eqp, cm, diff_pty, fm, cb);

    /* Add it to the local system */
    if (eqp->adv_scaling_property == NULL)
      cs_sdm_add(csys->mat, cb->loc);

    else {

      if (cs_property_is_uniform(eqp->adv_scaling_property))
        cs_sdm_add_mult(csys->mat,
                        eqp->adv_scaling_property->ref_value, cb->loc);
      else {
        cs_real_t scaling = cs_property_value_in_cell(cm,
                                                      eqp->adv_scaling_property,
                                                      cb->t_pty_eval);
        cs_sdm_add_mult(csys->mat, scaling, cb->loc);
      }

    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after advection", csys);
#endif
  }

  if (cs_equation_param_has_reaction(eqp)) { /* REACTION TERM
                                              * ============= */

    /* Update the value of the reaction property(ies) if needed */
    cs_equation_set_reaction_properties_cw(eqp, eqb, cm, cb);

    if (eqb->sys_flag & CS_FLAG_SYS_REAC_DIAG) {

      /* |c|*wvc = |dual_cell(v) cap c| */
      assert(cs_eflag_test(eqb->msh_flag, CS_FLAG_COMP_PVQ));
      const double  ptyc = cb->rpty_val * cm->vol_c;
      for (short int i = 0; i < cm->n_vc; i++)
        csys->mat->val[i*(cm->n_vc + 1)] += cm->wvc[i] * ptyc;

    }
    else {

      assert(cs_flag_test(eqb->sys_flag, CS_FLAG_SYS_MASS_MATRIX));

      /* Update local system matrix with the reaction term */
      cs_sdm_add_mult(csys->mat, cb->rpty_val, mass_hodge->matrix);

    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after reaction", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  First pass to apply boundary conditions enforced weakly in CDO-Vb
 *         schemes. Update the local system before applying the time scheme.
 *         Case of scalar-valued CDO-Vb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] diff_hodge  pointer to a discrete Hodge op. for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_svb_apply_weak_bc(const cs_equation_param_t     *eqp,
                   const cs_cdovb_scaleq_t       *eqc,
                   const cs_cell_mesh_t          *cm,
                   cs_face_mesh_t                *fm,
                   cs_hodge_t                    *diff_hodge,
                   cs_cell_sys_t                 *csys,
                   cs_cell_builder_t             *cb)
{
  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Neumann boundary conditions */
    if (csys->has_nhmg_neumann) {
      for (short int v  = 0; v < cm->n_vc; v++) {
          csys->rhs[v] += csys->neu_values[v];
      }
    }

    /* Contribution for the advection term: csys is updated inside
       (matrix and rhs) and Dirichlet BCs are handled inside */
    if (cs_equation_param_has_convection(eqp) &&
        ((cb->cell_flag & CS_FLAG_SOLID_CELL) == 0))
      eqc->add_advection_bc(cm, eqp, cb->t_bc_eval, fm, cb, csys);

    /* The enforcement of the Dirichlet has to be done after all
       other contributions */
    if (cs_equation_param_has_diffusion(eqp)) {

      if (csys->has_robin) {
        assert(eqc->enforce_robin_bc != NULL);
        eqc->enforce_robin_bc(eqp, cm, fm, diff_hodge, cb, csys);
      }

      if (csys->has_dirichlet) /* csys is updated inside (matrix and rhs) */
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
          eqc->enforce_dirichlet(eqp, cm, fm, diff_hodge, cb, csys);

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
 * \param[in, out] diff_hodge  pointer to a discrete Hodge op. for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_svb_enforce_values(const cs_equation_param_t     *eqp,
                    const cs_cdovb_scaleq_t       *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_face_mesh_t                *fm,
                    cs_hodge_t                    *diff_hodge,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  /* Internal enforcement of DoFs: Update csys (matrix and rhs) */
  if (csys->has_internal_enforcement) {

    cs_equation_enforced_internal_dofs(eqp, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement",
                       csys);
#endif
  }

  if (cs_cell_has_boundary_elements(cb) && csys->has_dirichlet) {

    /* Boundary element (through either vertices or faces) */

    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED) {

      /* csys is updated inside (matrix and rhs) */
      eqc->enforce_dirichlet(eqp, cm, fm, diff_hodge, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after strong BC treatment", csys);
#endif
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the residual normalization at the cellwise level according
 *         to the requested type of renormalization
 *         Case of scalar-valued CDO vertex-based schemes.
 *
 * \param[in]  type       type of renormalization
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 * \param[in]  csys       pointer to a cs_cell_sys_t structure
 *
 * \return the value of the cellwise contribution to the normalization of
 *         the residual
 */
/*----------------------------------------------------------------------------*/

static double
_svb_cw_rhs_normalization(cs_param_resnorm_type_t     type,
                          const cs_cell_mesh_t       *cm,
                          const cs_cell_sys_t        *csys)
{
  double  _rhs_norm = 0;

  if (type == CS_PARAM_RESNORM_WEIGHTED_RHS) {

    for (short int i = 0; i < csys->n_dofs; i++)
      _rhs_norm += cm->wvc[i] * csys->rhs[i]*csys->rhs[i];

    _rhs_norm = cm->vol_c * _rhs_norm;

  }
  else if (type == CS_PARAM_RESNORM_FILTERED_RHS) {

    if (csys->has_dirichlet || csys->has_internal_enforcement) {

      for (short int i = 0; i < csys->n_dofs; i++) {
        if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET)
          continue;
        else if (csys->intern_forced_ids[i] > -1)
          continue;
        else
          _rhs_norm += csys->rhs[i]*csys->rhs[i];
      }

    }
    else { /* No need to apply a filter */

      for (short int i = 0; i < csys->n_dofs; i++)
        _rhs_norm += csys->rhs[i]*csys->rhs[i];

    }

  } /* Type of residual normalization */

  return _rhs_norm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Perform the assembly step for scalar-valued CDO Vb schemes
 *
 * \param[in]      eqc    context for this kind of discretization
 * \param[in]      cm     pointer to a cellwise view of the mesh
 * \param[in]      csys   pointer to a cellwise view of the system
 * \param[in]      rs     pointer to a cs_range_set_t structure
 * \param[in, out] eqa    pointer to a cs_equation_assemble_t structure
 * \param[in, out] mav    pointer to a cs_matrix_assembler_values_t structure
 * \param[in, out] rhs    right-hand side array
 */
/*----------------------------------------------------------------------------*/

inline static void
_svb_assemble(const cs_cdovb_scaleq_t           *eqc,
              const cs_cell_mesh_t              *cm,
              const cs_cell_sys_t               *csys,
              const cs_range_set_t              *rs,
              cs_equation_assemble_t            *eqa,
              cs_matrix_assembler_values_t      *mav,
              cs_real_t                         *rhs)
{
  /* Matrix assembly */
  eqc->assemble(csys->mat, csys->dof_ids, rs, eqa, mav);

  /* RHS assembly */
#if CS_CDO_OMP_SYNC_SECTIONS > 0
# pragma omp critical
  {
    for (int v = 0; v < cm->n_vc; v++)
      rhs[cm->v_ids[v]] += csys->rhs[v];
  }

  if (eqc->source_terms != NULL) {
#   pragma omp critical
    {
      for (int v = 0; v < cm->n_vc; v++) /* Source term assembly */
        eqc->source_terms[cm->v_ids[v]] += csys->source[v];
    }
  }
#else  /* Use atomic barrier */

  for (int v = 0; v < cm->n_vc; v++)
#   pragma omp atomic
    rhs[cm->v_ids[v]] += csys->rhs[v];

  if (eqc->source_terms != NULL) {
    for (int v = 0; v < cm->n_vc; v++) /* Source term assembly */
#     pragma omp atomic
      eqc->source_terms[cm->v_ids[v]] += csys->source[v];
  }
#endif
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
  if (_svb_cell_system == NULL || _svb_cell_builder == NULL)
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

  BFT_MALLOC(_svb_cell_system, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(_svb_cell_builder, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    _svb_cell_system[i] = NULL;
    _svb_cell_builder[i] = NULL;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    _svb_cell_system[t_id] = cs_cell_sys_create(connect->n_max_vbyc,
                                                 connect->n_max_fbyc,
                                                 1, NULL);
    _svb_cell_builder[t_id] = _svb_create_cell_builder(connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  _svb_cell_system[0] = cs_cell_sys_create(connect->n_max_vbyc,
                                            connect->n_max_fbyc,
                                            1, NULL);
  _svb_cell_builder[0] = _svb_create_cell_builder(connect);

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

  *csys = _svb_cell_system[t_id];
  *cb = _svb_cell_builder[t_id];
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
    cs_cell_sys_free(&(_svb_cell_system[t_id]));
    cs_cell_builder_free(&(_svb_cell_builder[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(_svb_cell_system[0]));
  cs_cell_builder_free(&(_svb_cell_builder[0]));
#endif /* openMP */

  BFT_FREE(_svb_cell_system);
  BFT_FREE(_svb_cell_builder);
  _svb_cell_system = NULL;
  _svb_cell_builder = NULL;
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
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid type of equation.\n"
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
  eqb->msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PVQ | CS_FLAG_COMP_PE |
    CS_FLAG_COMP_EV;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
    CS_FLAG_COMP_FEQ | CS_FLAG_COMP_FV;

  bool  need_eigen =
    (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
     eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) ? true : false;

  /* Diffusion term */
  eqc->diffusion_hodge = NULL;
  eqc->get_stiffness_matrix = NULL;

  if (cs_equation_param_has_diffusion(eqp)) {

    eqc->diffusion_hodge = cs_hodge_init_context(connect,
                                                 eqp->diffusion_property,
                                                 &(eqp->diffusion_hodgep),
                                                 true,        /* tensor ? */
                                                 need_eigen); /* eigen ? */

    const cs_property_data_t  *diff_pty = eqc->diffusion_hodge[0]->pty_data;
    switch (eqp->diffusion_hodgep.algo) {

    case CS_HODGE_ALGO_COST:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
      eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ;
      if (diff_pty->is_iso || diff_pty->is_unity)
        eqc->get_stiffness_matrix = cs_hodge_vb_cost_get_iso_stiffness;
      else
        eqc->get_stiffness_matrix = cs_hodge_vb_cost_get_aniso_stiffness;
      break;

    case CS_HODGE_ALGO_BUBBLE:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
      eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ;
      if (diff_pty->is_iso || diff_pty->is_unity)
        eqc->get_stiffness_matrix = cs_hodge_vb_bubble_get_iso_stiffness;
      else
        eqc->get_stiffness_matrix = cs_hodge_vb_bubble_get_aniso_stiffness;
      break;

    case CS_HODGE_ALGO_OCS2:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ | CS_FLAG_COMP_SEF;
      eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_ocs2_get_aniso_stiffness;
      break;

    case CS_HODGE_ALGO_VORONOI:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
      eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_voro_get_stiffness;
      break;

    case CS_HODGE_ALGO_WBS:
      eqb->msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ
        | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ | CS_FLAG_COMP_PFC;
      eqc->get_stiffness_matrix = cs_hodge_vb_wbs_get_stiffness;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                (" %s: Invalid type of algorithm to build the diffusion term."),
                __func__);

    } /* Switch on Hodge algo. */

  } /* Diffusion term is requested */

  /* Boundary conditions */
  BFT_MALLOC(eqc->vtx_bc_flag, n_vertices, cs_flag_t);
  cs_equation_set_vertex_bc_flag(connect, eqb->face_bc, eqc->vtx_bc_flag);

  eqc->enforce_sliding = NULL;  /* Only useful for vector-valued eq. */
  eqc->enforce_robin_bc = NULL;
  if (cs_equation_param_has_robin_bc(eqp)) {

    assert(cs_equation_param_has_diffusion(eqp)); /* sanity check */

    switch (eqp->diffusion_hodgep.algo) {

    case CS_HODGE_ALGO_WBS:
      eqc->enforce_robin_bc = cs_cdo_diffusion_svb_wbs_robin;
      break;

    case CS_HODGE_ALGO_COST:
    case CS_HODGE_ALGO_BUBBLE:
    case CS_HODGE_ALGO_OCS2:
    case CS_HODGE_ALGO_VORONOI:
      eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_HFQ;
      eqc->enforce_robin_bc = cs_cdo_diffusion_svb_cost_robin;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                (" %s: Invalid type of algorithm with Robin boundaries."),
                __func__);

    } /* Switch on Hodge algo. */

  } /* Robin boundary conditions */

  /* Treatment of the Dirichlet boundary conditions */
  eqc->enforce_dirichlet = NULL;
  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_HFQ;
    switch (eqp->diffusion_hodgep.algo) {

    case CS_HODGE_ALGO_COST:
    case CS_HODGE_ALGO_OCS2:
    case CS_HODGE_ALGO_VORONOI:
    case CS_HODGE_ALGO_BUBBLE:
      eqc->enforce_dirichlet = cs_cdo_diffusion_svb_ocs_weak_dirichlet;
      break;
    case CS_HODGE_ALGO_WBS:
      eqc->enforce_dirichlet = cs_cdo_diffusion_svb_wbs_weak_dirichlet;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to enforce the Dirichlet BC.",
                __func__);

    } /* Switch on Hodge algo. */
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_HFQ;
    switch (eqp->diffusion_hodgep.algo) {

    case CS_HODGE_ALGO_COST:
    case CS_HODGE_ALGO_OCS2:
    case CS_HODGE_ALGO_VORONOI:
    case CS_HODGE_ALGO_BUBBLE:
      eqc->enforce_dirichlet = cs_cdo_diffusion_svb_ocs_wsym_dirichlet;
      break;
    case CS_HODGE_ALGO_WBS:
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

  /* Advection term */
  eqc->get_advection_matrix = NULL;
  eqc->add_advection_bc = NULL;

  if (cs_equation_param_has_convection(eqp)) {

    cs_xdef_type_t  adv_deftype =
      cs_advection_field_get_deftype(eqp->adv_field);

    switch (adv_deftype) {

    case CS_XDEF_BY_VALUE:
      eqb->msh_flag |= CS_FLAG_COMP_DFQ;
      break;
    case CS_XDEF_BY_ARRAY:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ;
      break;
    case CS_XDEF_BY_ANALYTIC_FUNCTION:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_SEF | CS_FLAG_COMP_PFQ;
      break;
    case CS_XDEF_BY_FIELD:
      if (eqp->adv_field->status & CS_ADVECTION_FIELD_LEGACY_FV)
        eqb->msh_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ;
      break;

    default: /* Nothing to add */
      break;
    }

    switch (eqp->adv_formulation) {

    case CS_PARAM_ADVECTION_FORM_CONSERV:

      switch (eqp->adv_scheme) {

      case CS_PARAM_ADVECTION_SCHEME_CENTERED:
        eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
        eqc->get_advection_matrix = cs_cdo_advection_vb_cencsv;
        break;

      case CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND:
        eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
        eqc->get_advection_matrix = cs_cdo_advection_vb_mcucsv;
        break;

      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
      case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
      case CS_PARAM_ADVECTION_SCHEME_SG:
        eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
        if (cs_equation_param_has_diffusion(eqp))
          eqc->get_advection_matrix = cs_cdo_advection_vb_upwcsv_wpty;
        else
          eqc->get_advection_matrix = cs_cdo_advection_vb_upwcsv;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid advection scheme for vertex-based schemes",
                  __func__);
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
        eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
        if (cs_equation_param_has_diffusion(eqp))
          eqc->get_advection_matrix = cs_cdo_advection_vb_upwnoc_wpty;
        else
          eqc->get_advection_matrix = cs_cdo_advection_vb_upwnoc;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid advection scheme for vertex-based scheme",
                  __func__);
      } /* Scheme */
      break; /* Formulation */

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of formulation for the advection term",
                __func__);
    }

    /* Boundary conditions for advection */
    eqb->bd_msh_flag |= CS_FLAG_COMP_PEQ;
    eqc->add_advection_bc = cs_cdo_advection_vb_bc;

  }
  else { /* No advection term is requested */

    if (eqp->default_enforcement != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
      eqb->sys_flag |= CS_FLAG_SYS_SYM; /* Algebraic system is symmetric */

  }

  /* A mass matrix can be requested either for the reaction term, the unsteady
     term or for the source term */
  cs_hodge_algo_t  mass_matrix_algo = CS_HODGE_ALGO_VORONOI;

  /* Reaction term */
  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->do_lumping)
      eqb->sys_flag |= CS_FLAG_SYS_REAC_DIAG;
    else {

      switch (eqp->reaction_hodgep.algo) {

      case CS_HODGE_ALGO_VORONOI:
        eqb->sys_flag |= CS_FLAG_SYS_REAC_DIAG;
        break;
      case CS_HODGE_ALGO_WBS:
        mass_matrix_algo = CS_HODGE_ALGO_WBS;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid choice of algorithm for the reaction term.",
                  __func__);
        break;
      }

    } /* Lumping or not lumping */

  } /* Reaction term is requested */

  /* Unsteady term */
  if (cs_equation_param_has_time(eqp)) {

    if (eqp->do_lumping)
      eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
    else {

      switch (eqp->time_hodgep.algo) {

      case CS_HODGE_ALGO_VORONOI:
        eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
        break;
      case CS_HODGE_ALGO_WBS:
        mass_matrix_algo = CS_HODGE_ALGO_WBS;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid choice of algorithm for the unsteady term.",
                  __func__);
        break;
      }

    } /* Lumping or not lumping */

  } /* Unsteady term is requested */

  /* Source term */
  eqc->source_terms = NULL;

  if (cs_equation_param_has_sourceterm(eqp)) {

    if (cs_equation_param_has_time(eqp)) {

      if (eqp->time_scheme == CS_TIME_SCHEME_THETA ||
          eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO) {

        BFT_MALLOC(eqc->source_terms, eqc->n_dofs, cs_real_t);
#       pragma omp parallel for if (eqc->n_dofs > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < eqc->n_dofs; i++)
          eqc->source_terms[i] = 0;

      } /* Theta scheme */

      /* Check the coherency of the settings --> Display a warning if something
         not consistent is found */
      for (int st_id = 0; st_id < eqp->n_source_terms; st_id++) {

        cs_xdef_t  *st = eqp->source_terms[st_id];

        if ((eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) &&
            (st->meta & CS_FLAG_DUAL)) {
          cs_base_warn(__FILE__, __LINE__);
          cs_log_printf(CS_LOG_DEFAULT,
                        "%s: A better choice for the reduction of the source"
                        " term is on primal entities.", __func__);
        }

      } /* Loop on the definitions of source terms */

    } /* Time-dependent equation */

    /* Need a mass matrix */
    for (int st_id = 0; st_id < eqp->n_source_terms; st_id++) {
      cs_xdef_t  *st = eqp->source_terms[st_id];
      if (st->meta & CS_FLAG_PRIMAL) {
        mass_matrix_algo = CS_HODGE_ALGO_WBS;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
      }
    }

  } /* There is at least one source term */

  /* Pre-defined a cs_hodge_param_t structure */
  eqc->mass_hodgep.inv_pty  = false;
  eqc->mass_hodgep.type = CS_HODGE_TYPE_VPCD;
  eqc->mass_hodgep.algo = mass_matrix_algo;
  eqc->mass_hodgep.coef = 1.0;  /* not useful in this case */

  if (mass_matrix_algo == CS_HODGE_ALGO_WBS)
    eqb->msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ
      | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_PFC;

  /* Array of hodge structure for the mass matrix */
  eqc->mass_hodge = cs_hodge_init_context(connect,
                                          NULL,
                                          &(eqc->mass_hodgep),
                                          false,  /* tensor ? */
                                          false); /* eigen ? */

  if (eqp->verbosity > 1 && eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) {
    cs_log_printf(CS_LOG_SETUP,
                  "#### Parameters of the mass matrix of the equation %s\n",
                  eqp->name);
    cs_hodge_param_log("Mass matrix", NULL, eqc->mass_hodgep);
  }

  /* Set the function pointer */
  eqc->get_mass_matrix = cs_hodge_get_func(__func__, eqc->mass_hodgep);

  /* Assembly process */
  eqc->assemble = cs_equation_assemble_set(CS_SPACE_SCHEME_CDOVB,
                                           CS_CDO_CONNECT_VTX_SCAL);

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

  cs_hodge_free_context(&(eqc->diffusion_hodge));
  cs_hodge_free_context(&(eqc->mass_hodge));

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

    cs_lnum_t  *def2v_ids = (cs_lnum_t *)cs_equation_get_tmpbuf();
    cs_lnum_t  *def2v_idx = NULL;
    BFT_MALLOC(def2v_idx, eqp->n_ic_defs + 1, cs_lnum_t);

    cs_equation_sync_vol_def_at_vertices(connect,
                                         eqp->n_ic_defs,
                                         eqp->ic_defs,
                                         def2v_idx,
                                         def2v_ids);

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_defs[def_id];
      const cs_lnum_t  n_v_selected = def2v_idx[def_id+1] - def2v_idx[def_id];
      const cs_lnum_t  *selected_lst = def2v_ids + def2v_idx[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_at_vertices_by_value(def,
                                                   n_v_selected,
                                                   selected_lst,
                                                   v_vals);
        break;

      case CS_XDEF_BY_QOV:
        /* Synchronization is performed inside */
        cs_evaluate_potential_by_qov(CS_FLAG_SCALAR | cs_flag_primal_vtx,
                                     def,
                                     v_vals, NULL);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        assert(eqp->dof_reduction == CS_PARAM_REDUCTION_DERHAM);
        cs_evaluate_potential_at_vertices_by_analytic(def,
                                                      t_eval,
                                                      n_v_selected,
                                                      selected_lst,
                                                      v_vals);
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
  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   _svb_cell_builder[0], /* static variable */
                                   eqc->vtx_bc_flag,
                                   v_vals);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_steady_state(bool                        cur2prev,
                                   const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Build an array storing the Dirichlet values at vertices and another one
     to detect vertices with an enforcement */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *forced_ids = NULL;

  /* First argument is set to t_cur even if this is a steady computation since
   * one can call this function to compute a steady-state solution at each time
   * step of an unsteady computation. */
  _svb_setup(time_eval, mesh, eqp, eqb, eqc->vtx_bc_flag,
             &dir_values, &forced_ids);

  if (eqb->init_step)
    eqb->init_step = false;

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;
  double  rhs_norm = 0.0;

  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
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
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = _svb_cell_system[t_id];
    cs_cell_builder_t  *cb = _svb_cell_builder[t_id];
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    /* Set times at which one evaluates quantities when needed */
    cb->t_pty_eval = time_eval; /* Dummy parameter if really steady */
    cb->t_bc_eval = time_eval;  /* Dummy parameter if really steady */
    cb->t_st_eval = time_eval;  /* Dummy parameter if really steady */

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, diff_hodge, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */
      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _svb_init_cell_system(cm, eqp, eqb,
                            dir_values, eqc->vtx_bc_flag, forced_ids, fld->val,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction terms into the local
       * system.
       * A mass matrix is also built if needed (stored in mass_hodge->matrix)
       */
      _svb_conv_diff_reac(eqp, eqb, eqc, cm,
                          fm, mass_hodge, diff_hodge, csys, cb);

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
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        csys->source);

        /* Update the RHS */
        for (short int v = 0; v < cm->n_vc; v++)
          csys->rhs[v] += csys->source[v];

      } /* End of term source */

      /* Compute a cellwise norm of the RHS for the normalization of the
         residual during the resolution of the linear system */
      rhs_norm += _svb_cw_rhs_normalization(eqp->sles_param->resnorm_type,
                                            cm, csys);

      /* Apply boundary conditions (those which are weakly enforced) */
      _svb_apply_weak_bc(eqp, eqc, cm, fm, diff_hodge, csys, cb);

      /* Enforce values if needed (internal or Dirichlet) */
      _svb_enforce_values(eqp, eqc, cm, fm, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* Assembly process
       * ================ */

      _svb_assemble(eqc, cm, csys, rs, eqa, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(forced_ids);

  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  if (cur2prev)
    cs_field_current_to_previous(fld);

  /* Solve the linear system */
  /* ======================= */

  /* Last step in the computation of the renormalization coefficient */
  cs_equation_sync_rhs_normalization(eqp->sles_param->resnorm_type,
                                     eqc->n_dofs,
                                     rhs,
                                     &rhs_norm);

  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);

  cs_equation_solve_scalar_system(eqc->n_dofs,
                                  eqp->sles_param,
                                  matrix,
                                  rs,
                                  rhs_norm,
                                  true, /* rhs_redux */
                                  sles,
                                  fld->val,
                                  rhs);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t1, &t2);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_sles_free(sles);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         Implicit time scheme is used to progress in time.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_implicit(bool                        cur2prev,
                               const cs_mesh_t            *mesh,
                               const int                   field_id,
                               const cs_equation_param_t  *eqp,
                               cs_equation_builder_t      *eqb,
                               void                       *context)
{
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_lnum_t  n_vertices = quant->n_vertices;

  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Build an array storing the Dirichlet values at vertices and another one
     to detect vertices with an enforcement */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *forced_ids = NULL;

  _svb_setup(ts->t_cur + ts->dt[0], mesh, eqp, eqb, eqc->vtx_bc_flag,
             &dir_values,  &forced_ids);

  if (eqb->init_step)
    eqb->init_step = false;

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;
  double  rhs_norm = 0.;

  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN)
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
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = _svb_cell_system[t_id];
    cs_cell_builder_t  *cb = _svb_cell_builder[t_id];
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    const cs_real_t  time_eval = ts->t_cur + ts->dt[0];
    const cs_real_t  inv_dtcur = 1./ts->dt[0];

    /* Set times at which one evaluates quantities if needed */
    cb->t_pty_eval = time_eval;
    cb->t_bc_eval = time_eval;
    cb->t_st_eval = time_eval;

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, diff_hodge, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */
      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _svb_init_cell_system(cm, eqp, eqb, dir_values, eqc->vtx_bc_flag,
                            forced_ids, fld->val,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed */
      _svb_conv_diff_reac(eqp, eqb, eqc, cm,
                          fm, mass_hodge, diff_hodge, csys, cb);

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
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        csys->source);

        for (short int v = 0; v < cm->n_vc; v++)
          csys->rhs[v] += csys->source[v];

      } /* End of term source */

      /* Apply boundary conditions (those which are weakly enforced) */
      _svb_apply_weak_bc(eqp, eqc, cm, fm, diff_hodge, csys, cb);

      /* Unsteady term + time scheme
       * =========================== */

      if (!(eqb->time_pty_uniform))
        cb->tpty_val = cs_property_value_in_cell(cm, eqp->time_property,
                                                 cb->t_pty_eval);

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        /* |c|*wvc = |dual_cell(v) cap c| */
        CS_CDO_OMP_ASSERT(cs_eflag_test(eqb->msh_flag, CS_FLAG_COMP_PVQ));
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
        const cs_sdm_t  *mass_mat = mass_hodge->matrix;

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

      /* Compute a norm of the RHS for the normalization of the residual
         of the linear system to solve */
      rhs_norm += _svb_cw_rhs_normalization(eqp->sles_param->resnorm_type,
                                            cm, csys);

      /* Enforce values if needed (internal or Dirichlet) */
      _svb_enforce_values(eqp, eqc, cm, fm, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* Assembly process
       * ================ */
      _svb_assemble(eqc, cm, csys, rs, eqa, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(forced_ids);
  cs_matrix_assembler_values_finalize(&mav);

  /* Copy current field values to previous values */
  if (cur2prev)
    cs_field_current_to_previous(fld);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Solve the linear system */
  /* ======================= */

  /* Last step in the computation of the renormalization coefficient */
  cs_equation_sync_rhs_normalization(eqp->sles_param->resnorm_type,
                                     eqc->n_dofs,
                                     rhs,
                                     &rhs_norm);

  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);

  cs_equation_solve_scalar_system(eqc->n_dofs,
                                  eqp->sles_param,
                                  matrix,
                                  rs,
                                  rhs_norm,
                                  true, /* rhs_redux */
                                  sles,
                                  fld->val,
                                  rhs);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t1, &t2);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_sles_free(sles);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         Theta time scheme is used to progress in time.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_theta(bool                        cur2prev,
                            const cs_mesh_t            *mesh,
                            const int                   field_id,
                            const cs_equation_param_t  *eqp,
                            cs_equation_builder_t      *eqb,
                            void                       *context)
{
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_time_step_t  *ts = cs_shared_time_step;

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Sanity checks */
  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_THETA ||
         eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO);

  /* Build an array storing the Dirichlet values at vertices and another one
     to detect vertices with an enforcement */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *forced_ids = NULL;

  _svb_setup(ts->t_cur + ts->dt[0], mesh, eqp, eqb, eqc->vtx_bc_flag,
             &dir_values, &forced_ids);

  /* Initialize the local system: rhs */
  cs_real_t  *rhs = NULL;
  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  double  rhs_norm = 0.;

  /* Initialize the local system: matrix */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_matrix_assembler_values_init(matrix, NULL, NULL);

  const double  tcoef = 1 - eqp->theta;

  /* Detect the first call (in this case, we compute the initial source term)*/
  bool  compute_initial_source = false;
  if (eqb->init_step) {

    compute_initial_source = true;
    eqb->init_step = false;

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

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
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
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = _svb_cell_system[t_id];
    cs_cell_builder_t  *cb = _svb_cell_builder[t_id];
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    const cs_real_t  t_cur = ts->t_cur;
    const cs_real_t  dt_cur = ts->dt[0];
    const cs_real_t  inv_dtcur = 1./dt_cur;

    /* Set times at which one evaluates quantities when needed
     * t_pty_eval = (1-theta).t^n + theta.t^(n+1) = t^n + theta.dt
     * since t^(n+1) = t^n + dt
     */
    cb->t_pty_eval = t_cur + eqp->theta*dt_cur;
    cb->t_bc_eval = t_cur + dt_cur;
    cb->t_st_eval = t_cur + dt_cur;

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, diff_hodge, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */
      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _svb_init_cell_system(cm, eqp, eqb, dir_values, eqc->vtx_bc_flag,
                            forced_ids, fld->val,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (mass_hodge->matrix) */
      _svb_conv_diff_reac(eqp, eqb, eqc, cm,
                          fm, mass_hodge, diff_hodge, csys, cb);

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
                                          mass_hodge,
                                          cb,
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
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        csys->source);

        for (short int v = 0; v < cm->n_vc; v++)
          csys->rhs[v] += eqp->theta * csys->source[v];

      } /* End of term source */

      /* Apply boundary conditions (those which are weakly enforced) */
      _svb_apply_weak_bc(eqp, eqc, cm, fm, diff_hodge, csys, cb);

      /* Unsteady term + time scheme
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
      if (!(eqb->time_pty_uniform))
        cb->tpty_val = cs_property_value_in_cell(cm, eqp->time_property,
                                                 cb->t_pty_eval);

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        /* |c|*wvc = |dual_cell(v) cap c| */
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
        const cs_sdm_t  *mass_mat = mass_hodge->matrix;

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

      /* Compute a norm of the RHS for the normalization of the residual
         of the linear system to solve */
      rhs_norm += _svb_cw_rhs_normalization(eqp->sles_param->resnorm_type,
                                            cm, csys);

      /* Enforce values if needed (internal or Dirichlet) */
      _svb_enforce_values(eqp, eqc, cm, fm, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* Assembly process
       * ================ */
      _svb_assemble(eqc, cm, csys, rs, eqa, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(forced_ids);
  cs_matrix_assembler_values_finalize(&mav);

  /* Copy current field values to previous values */
  if (cur2prev)
    cs_field_current_to_previous(fld);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Solve the linear system */
  /* ======================= */

  /* Last step in the computation of the renormalization coefficient */
  cs_equation_sync_rhs_normalization(eqp->sles_param->resnorm_type,
                                     eqc->n_dofs,
                                     rhs,
                                     &rhs_norm);

  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);

  cs_equation_solve_scalar_system(eqc->n_dofs,
                                  eqp->sles_param,
                                  matrix,
                                  rs,
                                  rhs_norm,
                                  true, /* rhs_redux */
                                  sles,
                                  fld->val,
                                  rhs);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t1, &t2);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_sles_free(sles);
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
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size: n_vertices)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_scaleq_get_vertex_values(void        *context,
                                  bool         previous)
{
  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;

  if (eqc == NULL)
    return NULL;

  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  if (previous)
    return pot->val_pre;
  else
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
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size: n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_scaleq_get_cell_values(void      *context,
                                bool       previous)
{
  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;

  if (eqc == NULL)
    return NULL;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  cs_real_t  *vtx_values = NULL;
  if (previous)
    vtx_values = pot->val_pre;
  else
    vtx_values = pot->val;
  assert(vtx_values != NULL);

  /* Reset buffer of values */
  if (eqc->cell_values == NULL)
    BFT_MALLOC(eqc->cell_values, quant->n_cells, cs_real_t);
  memset(eqc->cell_values, 0, quant->n_cells*sizeof(cs_real_t));

  /* Compute the values at cell centers from an interpolation of the field
     values defined at vertices */
  cs_reco_pv_at_cell_centers(connect->c2v, quant, vtx_values, eqc->cell_values);

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
# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(quant, connect, eqp, eqb, eqc, pot, eb, _svb_cell_builder)     \
  firstprivate(time_eval, inv_dtcur)
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
    cs_cell_builder_t  *cb = _svb_cell_builder[t_id];
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    /* Set times at which one evaluates quantities if needed */
    cb->t_pty_eval = time_eval;
    cb->t_bc_eval = time_eval;
    cb->t_st_eval = time_eval;

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, diff_hodge, cb);

    /* Set inside the OMP section so that each thread has its own value */
    cs_real_t  *p_cur = cs_cdo_local_get_d_buffer(t_id);
    cs_real_t  *p_prev = p_cur + connect->n_max_vbyc;
    cs_real_t  *p_theta = p_prev + connect->n_max_vbyc;

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */
      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the value of the current potential */
      for (short int v = 0; v < cm->n_vc; v++)
        p_cur[v] = pot->val[cm->v_ids[v]];

      if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX)
        /* stored in mass_hodge->matrix */
        eqc->get_mass_matrix(cm, mass_hodge, cb);

      /* Unsteady term */
      if (cs_equation_param_has_time(eqp)) {

        if (!(eqb->time_pty_uniform))
          cb->tpty_val = cs_property_value_in_cell(cm, eqp->time_property,
                                                   cb->t_pty_eval);

        /* Set the value of the previous potential */
        for (short int v = 0; v < cm->n_vc; v++)
          p_prev[v] = pot->val_pre[cm->v_ids[v]];

        if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) {

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
          cs_sdm_square_matvec(mass_hodge->matrix, dp, res);

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

      default: /* Implicit (Euler or BDF2) */
        for (short int v = 0; v < cm->n_vc; v++)
          p_theta[v] = p_cur[v];
        break;

      } /* Switch on time scheme */

      /* Reaction term */
      if (cs_equation_param_has_reaction(eqp)) {

        cs_equation_set_reaction_properties_cw(eqp, eqb, cm, cb);

        /* Define the local reaction property */
        const double  rpty_val = cb->rpty_val;

        cs_real_t  *res = cb->values;
        memset(res, 0, cm->n_vc*sizeof(cs_real_t));
        cs_sdm_square_matvec(mass_hodge->matrix, p_theta, res);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->reaction_term[cm->v_ids[v]] += rpty_val * res[v];
        }

      } /* Reaction */

      /* Diffusion term */
      if (cs_equation_param_has_diffusion(eqp)) {

        /* Set the diffusion property */
        if (!(eqb->diff_pty_uniform))
          cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, cb->cell_flag,
                                         diff_hodge);

        /* Define the local stiffness matrix: local matrix owned by the cellwise
           builder (store in cb->loc) */
        eqc->get_stiffness_matrix(cm, diff_hodge, cb);

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
        cs_property_data_t  *diff_pty =
          (diff_hodge == NULL) ? NULL : diff_hodge->pty_data;
        eqc->get_advection_matrix(eqp, cm, diff_pty, fm, cb);

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
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        src);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->source_term[cm->v_ids[v]] += src[v];
        }

      } /* End of term source */

      /* Boundary conditions */
      if (cb->cell_flag &  CS_FLAG_BOUNDARY_CELL_BY_FACE) {

        const cs_cdo_bc_face_t  *face_bc = eqb->face_bc;

        /* Identify which face is a boundary face */
        for (short int f = 0; f < cm->n_fc; f++) {
          const cs_lnum_t  bf_id = cm->f_ids[f] - cm->bface_shift;
          if (bf_id > -1) { /* Border face */

            /* Advective flux */
            if (cs_equation_param_has_convection(eqp)) {
              cs_advection_field_cw_boundary_f2v_flux(cm,
                                                      eqp->adv_field,
                                                      f,
                                                      cb->t_bc_eval,
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
                        (const cs_real_t (*)[3])diff_hodge->pty_data->tensor,
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

  } /* OpenMP Block */

  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++)
    eb->balance[v_id] =   eb->unsteady_term[v_id]
                        + eb->reaction_term[v_id]
                        + eb->diffusion_term[v_id]
                        + eb->advection_term[v_id]
                        + eb->source_term[v_id];

  /* Parallel or periodic synchronisation */
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
 * \param[in, out]  context    pointer to a scheme builder structure
 * \param[in, out]  vf_flux    pointer to the values of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_boundary_diff_flux(const cs_real_t              t_eval,
                                   const cs_equation_param_t   *eqp,
                                   const cs_real_t             *pdi,
                                   cs_equation_builder_t       *eqb,
                                   void                        *context,
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

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;

  assert(eqc->diffusion_hodge != NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(quant, connect, eqp, eqb, eqc, vf_flux, pdi,                   \
         _svb_cell_builder)                                             \
  firstprivate(t_eval)
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

    /* Set inside the OMP section so that each thread has its own value */
    double  *tmp = cs_cdo_local_get_d_buffer(t_id);
    cs_real_t  *pot = tmp, /* +1 for WBS */
               *flux = tmp + connect->n_max_vbyc + 1;

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_builder_t  *cb = _svb_cell_builder[t_id];
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_hodge_t  *hodge = eqc->diffusion_hodge[t_id];

    assert(hodge != NULL);

    /* Set times at which one evaluates quantities if needed */
    cb->t_pty_eval = cb->t_bc_eval = cb->t_st_eval = t_eval;

    /* msh_flag for Neumann and Robin BCs. Add add_flag for the other cases
       when one has to reconstruct a flux */
    cs_eflag_t  msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_FV;
    cs_eflag_t  add_flag = CS_FLAG_COMP_EV | CS_FLAG_COMP_FE |
      CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ;

    switch (eqp->diffusion_hodgep.algo) {

    case CS_HODGE_ALGO_COST:
    case CS_HODGE_ALGO_BUBBLE:
    case CS_HODGE_ALGO_VORONOI:
      add_flag |= CS_FLAG_COMP_DFQ;
      break;

    case CS_HODGE_ALGO_WBS:
      add_flag |= CS_FLAG_COMP_PVQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_FEQ |
        CS_FLAG_COMP_HFQ;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "Invalid Hodge algorithm");

    } /* Switch hodge algo. */

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_hodge_set_property_value(0, cb->t_pty_eval, 0, hodge);

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

          cs_equation_compute_neumann_sv(cb->t_bc_eval,
                                         face_bc->def_ids[bf_id],
                                         f,
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
          cs_equation_compute_robin(cb->t_bc_eval,
                                    face_bc->def_ids[bf_id],
                                    f,
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
            cs_hodge_set_property_value_cw(cm, t_eval, 0, hodge);

          /* Define a local buffer keeping the value of the discrete potential
             for the current cell */
          for (short int v = 0; v < cm->n_vc; v++)
            pot[v] = pdi[cm->v_ids[v]];

          if (eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_WBS) {

            /* Interpolate also the value of the potential at the cell center */
            pot[cm->n_vc] = cs_reco_cw_scalar_pv_at_cell_center(cm, pot);

            cs_cdo_diffusion_wbs_vbyf_flux(f, cm, pot, hodge, cb, flux);

          }
          else
            cs_cdo_diffusion_svb_vbyf_flux(f, cm, pot, hodge, cb, flux);

          /* Fill the global flux array */
          short int  n_vf = 0;
          for (int i = cm->f2v_idx[f]; i < cm->f2v_idx[f+1]; i++)
            _flx[n_vf++] = flux[cm->f2v_ids[i]];

        }
        break;

      } /* End of switch */

    } /* End of loop on boundary faces */

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
  cs_real_t  gc[3], pty_gc[3], pty_tens[3][3];

  if (ml_t == CS_MESH_LOCATION_BOUNDARY_FACES) { /* Belongs to only one cell */

    const cs_lnum_t  n_i_faces = connect->n_faces[CS_INT_FACES];
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
                                    eqp->diffusion_hodgep.inv_pty,
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
                                      eqp->diffusion_hodgep.inv_pty,
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
  if (diff_flux == NULL)
    return ;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* If no diffusion, return after resetting */
  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(diff_flux, 0, 3*quant->n_cells*sizeof(cs_real_t));
    return;
  }

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;

  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL && eqc!= NULL);
  assert(eqc->diffusion_hodge != NULL);

  cs_timer_t  t0 = cs_timer_time();

#pragma omp parallel if (quant->n_cells > CS_THR_MIN)                   \
  shared(t_eval, quant, connect, eqp, eqb, diff_flux, values,           \
         _svb_cell_builder)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_builder_t  *cb = _svb_cell_builder[t_id];
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_hodge_t  *diff_hodge = eqc->diffusion_hodge[t_id];
    cs_eflag_t  msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_EV;
    cs_cdo_diffusion_cw_flux_t  *compute_flux =
      (eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_WBS) ?
      cs_cdo_diffusion_wbs_get_cell_flux : cs_cdo_diffusion_svb_get_cell_flux;

    /* Set inside the OMP section so that each thread has its own value */
    double  *pot = cs_cdo_local_get_d_buffer(t_id);

    if (eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_WBS)
      msh_flag |= CS_FLAG_COMP_PVQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ |
        CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ;
    else
      msh_flag |= CS_FLAG_COMP_DFQ | CS_FLAG_COMP_PVQ;

    /* Set times at which one evaluates quantities if needed */
    cb->t_pty_eval = cb->t_bc_eval = cb->t_st_eval = t_eval;

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_hodge_set_property_value(0, cb->t_pty_eval, 0, diff_hodge);

    /* Define the flux by cellwise contributions */
#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 3
      if (cs_dbg_cw_test(eqp, cm, NULL)) cs_cell_mesh_dump(cm);
#endif

      if (!eqb->diff_pty_uniform)
        cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, 0, diff_hodge);

      /* Define a local buffer keeping the value of the discrete potential
         for the current cell */
      for (short int v = 0; v < cm->n_vc; v++)
        pot[v] = values[cm->v_ids[v]];

      if (eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_WBS) {
        /* Interpolate also the value of the potential at the cell center */
        pot[cm->n_vc] = 0.;
        for (short int v = 0; v < cm->n_vc; v++)
          pot[cm->n_vc] += cm->wvc[v]*pot[v];
      }

      compute_flux(cm, pot, diff_hodge, cb, diff_flux + 3*c_id);

    } /* Loop on cells */

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
  if (diff_flux == NULL)
    return ;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* If no diffusion, return after resetting */
  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(diff_flux, 0, connect->c2e->idx[quant->n_cells]*sizeof(cs_real_t));
    return;
  }

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;

  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL && eqc!= NULL);
  assert(eqc->diffusion_hodge != NULL);

  cs_timer_t  t0 = cs_timer_time();

  cs_hodge_compute_t  *get_diffusion_hodge =
    cs_hodge_get_func(__func__, eqp->diffusion_hodgep);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(t_eval, quant, connect, eqp, eqb, diff_flux, values,           \
         get_diffusion_hodge, _svb_cell_builder)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_builder_t  *cb = _svb_cell_builder[t_id];
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_hodge_t  *diff_hodge = eqc->diffusion_hodge[t_id];
    cs_eflag_t  msh_flag = CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PV |
      CS_FLAG_COMP_PVQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_DFQ;
    cs_cdo_diffusion_cw_flux_t  *compute_flux =
      cs_cdo_diffusion_svb_get_dfbyc_flux;

    if (eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_WBS) {
      msh_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_SEF |
        CS_FLAG_COMP_HFQ | CS_FLAG_COMP_PFC | CS_FLAG_COMP_FEQ;
      compute_flux = cs_cdo_diffusion_wbs_get_dfbyc_flux;
    }

    /* Set inside the OMP section so that each thread has its own value */
    double  *pot = cs_cdo_local_get_d_buffer(t_id);

    /* Set times at which one evaluates quantities if needed */
    cb->t_pty_eval = cb->t_bc_eval = cb->t_st_eval = t_eval;

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_hodge_set_property_value(0, cb->t_pty_eval, 0, diff_hodge);

    /* Define the flux by cellwise contributions */
#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
      if (cs_dbg_cw_test(eqp, cm, NULL)) cs_cell_mesh_dump(cm);
#endif

      if (!eqb->diff_pty_uniform)
        cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, 0, diff_hodge);

      /* Define a local buffer keeping the value of the discrete potential
         for the current cell */
      for (short int v = 0; v < cm->n_vc; v++)
        pot[v] = values[cm->v_ids[v]];

      if (eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_WBS) {
        /* Interpolate also the value of the potential at the cell center */
        pot[cm->n_vc] = 0.;
        for (short int v = 0; v < cm->n_vc; v++)
          pot[cm->n_vc] += cm->wvc[v]*pot[v];
      }
      else {

        /* Compute the local Hodge operator */
        get_diffusion_hodge(cm, diff_hodge, cb);

      }

      /* Compute and store the flux */
      compute_flux(cm, pot, diff_hodge,
                   cb, diff_flux + connect->c2e->idx[c_id]);

    } /* Loop on cells */

  } /* OMP Section */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate a current to previous operation for the field associated to
 *         this equation and potentially for related fields/arrays.
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_current_to_previous(const cs_equation_param_t  *eqp,
                                    cs_equation_builder_t      *eqb,
                                    void                       *context)
{
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(eqc->var_field_id);

  cs_field_current_to_previous(fld);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_extra_post(const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *context)
{
  CS_UNUSED(context);

  const cs_timer_t  t0 = cs_timer_time();

  if (cs_equation_param_has_convection(eqp)) {
    if (eqp->process_flag & CS_EQUATION_POST_UPWIND_COEF) {

      int  len = strlen(eqp->name) + 8 + 1;
      char *postlabel = NULL;
      BFT_MALLOC(postlabel, len, char);
      sprintf(postlabel, "%s.UpwCoef", eqp->name);

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
