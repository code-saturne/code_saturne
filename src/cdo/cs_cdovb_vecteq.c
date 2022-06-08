/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
 * reaction of vector-valued equations with source terms
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
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_local.h"
#include "cs_cdo_solve.h"
#include "cs_cdo_system.h"
#include "cs_cdo_toolbox.h"
#include "cs_cdovb_priv.h"
#include "cs_equation_bc.h"
#include "cs_equation_builder.h"
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

#include "cs_cdovb_vecteq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 *  \file cs_cdovb_vecteq.c
 *
 * \brief Build an algebraic CDO vertex-based system for unsteady
 *        convection-diffusion-reaction of vector-valued equations with
 *        source terms
 *
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVB_VECTEQ_DBG     0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private variables
 *============================================================================*/

/* Structure to enable a full cellwise strategy during the system building */

static cs_cell_sys_t      **_vvb_cell_system = NULL;
static cs_cell_builder_t  **_vvb_cell_builder = NULL;

/* Pointer to shared structures */

static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the local builder structure used for building the system
 *         cellwise
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_cell_builder_t *
_vvb_create_cell_builder(const cs_cdo_connect_t   *connect)
{
  const int  n_vc = connect->n_max_vbyc;
  const int  n_ec = connect->n_max_ebyc;

  cs_cell_builder_t  *cb = cs_cell_builder_create();

  BFT_MALLOC(cb->ids, n_vc, int);
  memset(cb->ids, 0, n_vc*sizeof(int));

  int  size = n_ec*(n_ec+1);
  size = CS_MAX(4*n_ec + 3*n_vc, size);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(cs_real_t));

  size = 2*n_ec;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of
     operators */

  cb->aux = cs_sdm_square_create(n_ec);
  cb->loc = cs_sdm_block33_create(n_vc, n_vc);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *         Define an indirection array for the enforcement of internal DoFs
 *         only if needed
 *
 * \param[in]      t_eval        time at which one evaluates BCs
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      eqp           pointer to a cs_equation_param_t structure
 * \param[in, out] eqb           pointer to a cs_equation_builder_t structure
 * \param[in, out] vtx_bc_flag   pointer to an array of BC flag for each vertex
 */
/*----------------------------------------------------------------------------*/

static void
_vvb_setup(cs_real_t                      t_eval,
           const cs_mesh_t               *mesh,
           const cs_equation_param_t     *eqp,
           cs_equation_builder_t         *eqb,
           cs_flag_t                      vtx_bc_flag[])
{
  assert(vtx_bc_flag != NULL);  /* Sanity check */
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Compute the values of the Dirichlet BC */

  BFT_MALLOC(eqb->dir_values, 3*quant->n_vertices, cs_real_t);

  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   _vvb_cell_builder[0], /* static variable */
                                   vtx_bc_flag,
                                   eqb->dir_values);

  /* Internal enforcement of DoFs */

  if (cs_equation_param_has_internal_enforcement(eqp))
    eqb->enforced_values =
      cs_enforcement_define_at_vertices(connect,
                                        eqp->n_enforcements,
                                        eqp->enforcement_params);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell.
 *          Case of vector-valued CDO-Vb schemes
 *
 * \param[in]      cm           pointer to a cellwise view of the mesh
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in]      vtx_bc_flag  Flag related to BC associated to each vertex
 * \param[in]      field_tn     values of the field at the last computed time
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vvb_init_cell_system(const cs_cell_mesh_t           *cm,
                      const cs_equation_param_t      *eqp,
                      const cs_equation_builder_t    *eqb,
                      const cs_flag_t                 vtx_bc_flag[],
                      const cs_real_t                 field_tn[],
                      cs_cell_sys_t                  *csys,
                      cs_cell_builder_t              *cb)
{
  csys->c_id = cm->c_id;
  csys->n_dofs = 3*cm->n_vc;

  /* Cell-wise view of the linear system to build:
     Initialize the local system */

  cs_cell_sys_reset(cm->n_fc, csys); /* Generic part */

  cs_sdm_block33_init(csys->mat, cm->n_vc, cm->n_vc);

  for (int v = 0; v < cm->n_vc; v++) {
    const cs_lnum_t  v_id = cm->v_ids[v];
    const cs_real_t  *val_p = field_tn + 3*v_id;
    for (int k = 0; k < 3; k++) {
      csys->dof_ids[3*v + k] = 3*v_id + k;
      csys->val_n[3*v + k] = val_p[k];
    }
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Set the bc (specific part) */

    cs_equation_bc_set_cw_vb(cm,
                             eqp,
                             eqb->face_bc,
                             vtx_bc_flag,
                             eqb->dir_values,
                             cb->t_bc_eval,
                             csys,
                             cb);

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    cs_dbg_check_hmg_dirichlet_cw(__func__, csys);
#endif
  } /* Border cell */

  /* Special case to handle if there is an enforcement by penalization or
   * algebraic.  This situation may happen with a tetrahedron with one vertex
   * or an edge lying on the boundary (but no face)
   */

  if (cb->cell_flag == CS_FLAG_BOUNDARY_CELL_BY_VERTEX) {

    assert(vtx_bc_flag != NULL);

    for (int v = 0; v < cm->n_vc; v++) {

      for (int k = 0; k < 3; k++)
        csys->dof_flag[3*v+k] = vtx_bc_flag[cm->v_ids[v]];
      if (cs_cdo_bc_is_dirichlet(vtx_bc_flag[cm->v_ids[v]])) {
        csys->has_dirichlet = true;
        const cs_real_t  *bc_val = eqb->dir_values + 3*cm->v_ids[v];
        for (int k = 0; k < 3; k++)
          csys->dir_values[3*v+k] = bc_val[k];
      }

    }

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion, advection,
 *          reaction terms. If asked, a mass matrix is also computed and
 *          stored in mass_hodge->matrix
 *          Case of vector-valued CDO-Vb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] mass_hodge  pointer to a cs_hodge_t structure (mass matrix)
 * \param[in, out] diff_hodge  pointer to a cs_hodge_t structure (diffusion)
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vvb_conv_diff_reac(const cs_equation_param_t     *eqp,
                    const cs_equation_builder_t   *eqb,
                    const cs_cdovb_vecteq_t       *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_hodge_t                    *mass_hodge,
                    cs_hodge_t                    *diff_hodge,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
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

    const cs_real_t  *sval = cb->loc->val;
    for (int bi = 0; bi < cm->n_vc; bi++) {
      for (int bj = 0; bj < cm->n_vc; bj++) {

        /* Retrieve the 3x3 matrix */

        cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bi, bj);
        assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

        const cs_real_t  _val = sval[cm->n_vc*bi+bj];
        bij->val[0] += _val;
        bij->val[4] += _val;
        bij->val[8] += _val;

      }
    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after adding diffusion", csys);
#endif
  }

  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                  * =========== */
    assert(mass_hodge != NULL);

    /* Build the mass matrix and store it in mass_hodge->matrix */

    eqc->get_mass_matrix(cm, mass_hodge, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys)) {
      cs_log_printf(CS_LOG_DEFAULT, ">> Cell mass matrix");
      cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, mass_hodge->matrix);
    }
#endif
  }

  if (cs_equation_param_has_reaction(eqp)) { /* REACTION TERM
                                              * ============= */

    /* Update the value of the reaction property(ies) if needed */

    cs_equation_builder_set_reaction_pty_cw(eqp, eqb, cm, cb);

    if (eqb->sys_flag & CS_FLAG_SYS_REAC_DIAG) {

      /* |c|*wvc = |dual_cell(v) cap c| */

      assert(cs_eflag_test(eqb->msh_flag, CS_FLAG_COMP_PVQ));
      const double  ptyc = cb->rpty_val * cm->vol_c;

      /* Only the diagonal block and its diagonal entries are modified */

      for (int bi = 0; bi < cm->n_vc; bi++) {

        const double  vpty = cm->wvc[bi] * ptyc;

        /* Retrieve the 3x3 matrix */

        cs_sdm_t  *bii = cs_sdm_get_block(csys->mat, bi, bi);
        assert(bii->n_rows == 3 && bii->n_cols == 3);

        bii->val[0] += vpty;
        bii->val[4] += vpty;
        bii->val[8] += vpty;

      }

    }
    else {

      assert(cs_flag_test(eqb->sys_flag, CS_FLAG_SYS_MASS_MATRIX));

      /* Add the local mass matrix to the local system */

      const cs_real_t  *mval = mass_hodge->matrix->val;
      for (int bi = 0; bi < cm->n_vc; bi++) {
        for (int bj = 0; bj < cm->n_vc; bj++) {

          /* Retrieve the 3x3 matrix */

          cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bi, bj);
          assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

          const cs_real_t  _val = cb->rpty_val * mval[cm->n_vc*bi+bj];
          bij->val[0] += _val;
          bij->val[4] += _val;
          bij->val[8] += _val;

        }
      }

    } /* Lumping or not lumping */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after adding reaction", csys);
#endif
  } /* Reaction term */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   First pass to apply boundary conditions enforced weakly. Update
 *          the local system before applying the time scheme.
 *          Case of vector-valued CDO-Vb schemes
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
_vvb_apply_weak_bc(const cs_equation_param_t     *eqp,
                   const cs_cdovb_vecteq_t       *eqc,
                   const cs_cell_mesh_t          *cm,
                   cs_face_mesh_t                *fm,
                   cs_hodge_t                    *diff_hodge,
                   cs_cell_sys_t                 *csys,
                   cs_cell_builder_t             *cb)
{
  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Neumann boundary conditions:
     * The common practice is to define Phi_neu = - lambda * \grad(u) . n_fc
     * An outward flux is a positive flux whereas an inward flux is negative
     * The minus just above implies the minus just below */

    if (csys->has_nhmg_neumann) {
      for (short int i  = 0; i < 3*cm->n_vc; i++)
        csys->rhs[i] -= csys->neu_values[i];
    }

    /* The enforcement of the Dirichlet has to be done after all
       other contributions */

    if (cs_equation_param_has_diffusion(eqp)) {

      if (csys->has_dirichlet) /* csys is updated inside (matrix and rhs) */
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
          eqc->enforce_dirichlet(eqp, cm, fm, diff_hodge, cb, csys);

    }

    if (csys->has_sliding)
      eqc->enforce_sliding(eqp, cm, fm, diff_hodge, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after BC weak treatment", csys);
#endif
  } /* Cell with at least one boundary face */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Second pass to apply boundary conditions. Only Dirichlet BCs which
 *          are enforced strongly. Apply also the enforcement of internal DoFs.
 *          Update the local system after applying the time scheme.
 *          Case of vector-valued CDO-Vb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] diff_hodge  pointer to a discrete Hodge op. for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vvb_enforce_values(const cs_equation_param_t     *eqp,
                    const cs_equation_builder_t   *eqb,
                    const cs_cdovb_vecteq_t       *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_face_mesh_t                *fm,
                    cs_hodge_t                    *diff_hodge,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  if (cs_cell_has_boundary_elements(cb) && csys->has_dirichlet) {

    /* Boundary element (through either vertices or faces) */

    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED) {

      /* csys is updated inside (matrix and rhs) */

      eqc->enforce_dirichlet(eqp, cm, fm, diff_hodge, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 2
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after strong BC treatment", csys);
#endif
    }
  }

  if (cs_equation_param_has_internal_enforcement(eqp)) {

    /* Internal enforcement of DoFs: Update csys (matrix and rhs) */

    cs_equation_builder_enforce_block_dofs(eqb, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 2
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement",
                       csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb scheme
 *
 * \param[in]      csys              pointer to a cs_cell_sys_t structure
 * \param[in]      block             pointer to a block structure
 * \param[in, out] rhs               array of values for the rhs
 * \param[in, out] eqc               context structure for a vector-valued Fb
 * \param[in, out] asb               pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_vvb_assembly(const cs_cell_sys_t            *csys,
              const cs_cdo_system_block_t    *block,
              cs_real_t                      *rhs,
              cs_cdovb_vecteq_t              *eqc,
              cs_cdo_assembly_t              *asb)
{
  assert(asb != NULL && block != NULL); /* Sanity check */
  assert(block->type == CS_CDO_SYSTEM_BLOCK_DEFAULT);
  cs_cdo_system_dblock_t  *db = block->block_pointer;

  /* Matrix assembly */

  db->assembly_func(csys->mat, csys->dof_ids, db->range_set, asb, db->mav);

  /* RHS and if needed source term assembly */

  if (eqc->source_terms != NULL) {

#   pragma omp critical
    for (short int v = 0; v < csys->n_dofs; v++) {

      cs_lnum_t  v_id = csys->dof_ids[v];

      rhs[v_id] += csys->rhs[v];
      eqc->source_terms[v_id] += csys->source[v];

    } /* Loop on DoFs */

  }
  else { /* No source term */

#   pragma omp critical
    for (short int v = 0; v < csys->n_dofs; v++)
      rhs[csys->dof_ids[v]] += csys->rhs[v];

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the residual normalization at the cellwise level according
 *         to the requested type of renormalization
 *         Case of vector-valued system.
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
_vvb_cw_rhs_normalization(cs_param_resnorm_type_t     type,
                          const cs_cell_mesh_t       *cm,
                          const cs_cell_sys_t        *csys)
{
  double  _rhs_norm = 0;

  switch (type) {

  case CS_PARAM_RESNORM_WEIGHTED_RHS:
    for (short int i = 0; i < cm->n_vc; i++) {
      const cs_real_t  w = cm->wvc[i];
      const cs_real_t  *_rhs = csys->rhs + 3*i;
      _rhs_norm += w * (_rhs[0]*_rhs[0] + _rhs[1]*_rhs[1] + _rhs[2]*_rhs[2]);
    }
    _rhs_norm = cm->vol_c * _rhs_norm;
    break;

  case CS_PARAM_RESNORM_FILTERED_RHS:
    for (short int i = 0; i < csys->n_dofs; i++) {
      if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET)
        continue;
      else if (csys->dof_is_forced[i])
        continue;
      else
        _rhs_norm += csys->rhs[i]*csys->rhs[i];
    }
    break;

  case CS_PARAM_RESNORM_NORM2_RHS:
    for (short int i = 0; i < csys->n_dofs; i++)
      _rhs_norm += csys->rhs[i]*csys->rhs[i];
    break;

  default:
    break; /* Nothing to do */

  } /* Type of residual normalization */

  return _rhs_norm;
}

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

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
cs_cdovb_vecteq_is_initialized(void)
{
  if (_vvb_cell_system == NULL || _vvb_cell_builder == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Allocate work buffer and general structures related to CDO
 *           vector-valued vertex-based schemes
 *           Set shared pointers.
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_init_sharing(const cs_cdo_quantities_t    *quant,
                             const cs_cdo_connect_t       *connect,
                             const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */

  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  /* Structure used to build the final system by a cell-wise process */

  assert(cs_glob_n_threads > 0);
  BFT_MALLOC(_vvb_cell_system, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(_vvb_cell_builder, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    _vvb_cell_system[i] = NULL;
    _vvb_cell_builder[i] = NULL;
  }

  const int  n_max_dofs = 3*connect->n_max_vbyc;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_builder_t  *cb = _vvb_create_cell_builder(connect);
    _vvb_cell_builder[t_id] = cb;

    int  block_size = 3;
    _vvb_cell_system[t_id] = cs_cell_sys_create(n_max_dofs,
                                                connect->n_max_fbyc,
                                                1,
                                                &block_size);
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_builder_t  *cb = _vvb_create_cell_builder(connect);
  _vvb_cell_builder[0] = cb;

  int  block_size = 3;
  _vvb_cell_system[0] =  cs_cell_sys_create(n_max_dofs,
                                            connect->n_max_fbyc,
                                            1,
                                            &block_size);
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
cs_cdovb_vecteq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  *csys = _vvb_cell_system[t_id];
  *cb = _vvb_cell_builder[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO vertex-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_finalize_sharing(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(_vvb_cell_system[t_id]));
    cs_cell_builder_free(&(_vvb_cell_builder[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(_vvb_cell_system[0]));
  cs_cell_builder_free(&(_vvb_cell_builder[0]));
#endif /* openMP */

  BFT_FREE(_vvb_cell_system);
  BFT_FREE(_vvb_cell_builder);
  _vvb_cell_builder = NULL;
  _vvb_cell_system = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_vecteq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_vecteq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVB || eqp->dim != 3)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of equation.\n"
              " Expected: vector-valued CDO vertex-based equation.", __func__);

  eqb->sys_flag = CS_FLAG_SYS_VECTOR;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->n_vertices;

  cs_cdovb_vecteq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdovb_vecteq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  eqc->n_dofs = 3*n_vertices;

  /* Flag to indicate the minimal set of quantities to build in a cell mesh
     According to the situation, additional flags have to be set */

  eqb->msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PVQ | CS_FLAG_COMP_PE |
    CS_FLAG_COMP_EV;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */

  eqb->bd_msh_flag = CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
    CS_FLAG_COMP_FEQ;

  bool  need_eigen =
    (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
     eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) ? true : false;

  /* Diffusion term */

  eqc->get_stiffness_matrix = NULL;
  eqc->get_stiffness_matrix = NULL;

  if (cs_equation_param_has_diffusion(eqp)) {

    eqc->diffusion_hodge = cs_hodge_init_context(connect,
                                                 eqp->diffusion_property,
                                                 &(eqp->diffusion_hodgep),
                                                 true,        /* tensor ? */
                                                 need_eigen); /* eigen ? */

    const cs_property_data_t  *diff_pty = eqc->diffusion_hodge[0]->pty_data;

    if (diff_pty->is_iso == false)
      bft_error(__FILE__, __LINE__, 0, " %s: Case not handle yet\n", __func__);

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
      eqb->msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ |
        CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ | CS_FLAG_COMP_PFC;
      eqc->get_stiffness_matrix = cs_hodge_vb_wbs_get_stiffness;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to build the diffusion term.",
                __func__);

    } /* Switch on Hodge algo. */

  } /* Has diffusion */

  /* Boundary conditions */

  BFT_MALLOC(eqc->vtx_bc_flag, n_vertices, cs_flag_t);
  cs_equation_bc_set_vertex_flag(connect, eqb->face_bc, eqc->vtx_bc_flag);

  eqc->enforce_robin_bc = NULL;
  if (cs_equation_param_has_robin_bc(eqp))
    bft_error(__FILE__, __LINE__, 0,
              (" %s: Robin boundary conditions are not handled yet."),
              __func__);

  eqc->enforce_dirichlet = NULL;
  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_block_dirichlet;
    break;
  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_block_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_HFQ;
    eqc->enforce_dirichlet = cs_cdo_diffusion_vvb_ocs_weak_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  eqc->enforce_sliding = NULL;
  if (eqb->face_bc->n_sliding_faces > 0) {

    /* There is at least one face with a sliding condition to handle */

    eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_HFQ;
    eqc->enforce_sliding = cs_cdo_diffusion_vvb_ocs_sliding;

  }

  /* Advection term */

  eqc->get_advection_matrix = NULL;
  eqc->add_advection_bc = NULL;

  /* A mass matrix can be requested either for the reaction term, the unsteady
     term or for the source term */

  cs_hodge_algo_t  reac_hodge_algo = CS_HODGE_N_ALGOS;
  cs_hodge_algo_t  time_hodge_algo = CS_HODGE_N_ALGOS;
  cs_hodge_algo_t  srct_hodge_algo = CS_HODGE_N_ALGOS;

  /* Reaction term */

  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->do_lumping) {

      eqb->sys_flag |= CS_FLAG_SYS_REAC_DIAG;
      reac_hodge_algo = CS_HODGE_ALGO_VORONOI;

    }
    else {

      switch (eqp->reaction_hodgep.algo) {

      case CS_HODGE_ALGO_VORONOI:
        eqb->sys_flag |= CS_FLAG_SYS_REAC_DIAG;
        reac_hodge_algo = CS_HODGE_ALGO_VORONOI;
        break;
      case CS_HODGE_ALGO_WBS:
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        reac_hodge_algo = CS_HODGE_ALGO_WBS;
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

    if (eqp->do_lumping) {

      eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
      time_hodge_algo = CS_HODGE_ALGO_VORONOI;

    }
    else {

      switch (eqp->time_hodgep.algo) {

      case CS_HODGE_ALGO_VORONOI:
        eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
        time_hodge_algo = CS_HODGE_ALGO_VORONOI;
        break;
      case CS_HODGE_ALGO_WBS:
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        time_hodge_algo = CS_HODGE_ALGO_WBS;
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
        memset(eqc->source_terms, 0, eqc->n_dofs*sizeof(cs_real_t));

      } /* Theta scheme */
    } /* Time-dependent equation */

    /* Need a mass matrix */

    for (int st_id = 0; st_id < eqp->n_source_terms; st_id++) {

      cs_xdef_t  *st = eqp->source_terms[st_id];
      if (st->meta & CS_FLAG_PRIMAL) {
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        srct_hodge_algo = CS_HODGE_ALGO_WBS;
      }

    }

  } /* There is at least one source term */

  /* Pre-defined a cs_hodge_builder_t structure */

  eqc->mass_hodgep.inv_pty  = false;
  eqc->mass_hodgep.coef = 1.0;  /* not useful in this case */
  eqc->mass_hodgep.type = CS_HODGE_TYPE_VPCD;

  eqc->mass_hodgep.algo = cs_hodge_set_mass_algo(eqp->name,
                                                 reac_hodge_algo,
                                                 time_hodge_algo,
                                                 srct_hodge_algo);

  if (eqc->mass_hodgep.algo == CS_HODGE_ALGO_WBS)
    eqb->msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ
      | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_PFC;

  /* Initialize the hodge structure for the mass matrix */

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

  /* Helper structures (range set, interface set, matrix structure and all the
     assembly process) */

  cs_cdo_system_helper_t  *sh = NULL;
  cs_lnum_t  col_block_size = 3*n_vertices;

  sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_DEFAULT,
                                   1,                /* n_col_blocks */
                                   &col_block_size,  /* col_block_size */
                                   1);               /* n_blocks */

  /* Choose the right class of matrix to avoid copy.
   * The way to perform the assembly may change if an external librairy is used
   * for solving the linear system */

  cs_cdo_system_matrix_class_t  matclass;
  bool is_unrolled;

  switch (eqp->sles_param->solver_class) {

  case CS_PARAM_SLES_CLASS_CS:
    is_unrolled = false;
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
    break;

  case CS_PARAM_SLES_CLASS_HYPRE:
#if defined(HAVE_HYPRE)
    matclass = CS_CDO_SYSTEM_MATRIX_HYPRE;
#else
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
#endif
    break;

  default:
    is_unrolled = true;
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
    break;

  }

  cs_cdo_system_add_dblock(sh, 0,         /* block_id */
                           matclass,
                           cs_flag_primal_vtx,
                           n_vertices,
                           3,             /* stride */
                           true,          /* interlaced */
                           is_unrolled);  /* unrolled */

  cs_cdo_system_build_block(sh, 0); /* block_id */

  eqb->system_helper = sh;

  /* Array used for extra-operations */

  eqc->cell_values = NULL;

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovb_vecteq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovb_vecteq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_vecteq_free_context(void   *builder)
{
  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)builder;

  if (eqc == NULL)
    return eqc;

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
 *         Case of vector-valued CDO-Vb schemes.
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
cs_cdovb_vecteq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_real_t  *v_vals = fld->val;

  /* By default, 0 is set as initial condition for the computational domain.

     Warning: This operation has to be done after the settings of the
     Dirichlet boundary conditions where an interface sum is performed
     for vertex-based schemes
  */

  memset(v_vals, 0, 3*quant->n_vertices*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    cs_lnum_t  *def2v_ids = (cs_lnum_t *)cs_cdo_toolbox_get_tmpbuf();
    cs_lnum_t  *def2v_idx = NULL;

    BFT_MALLOC(def2v_idx, eqp->n_ic_defs + 1, cs_lnum_t);

    cs_cdo_sync_vol_def_at_vertices(eqp->n_ic_defs, eqp->ic_defs,
                                    def2v_idx, def2v_ids);

    /* Initialize values at mesh vertices */

    cs_flag_t  dof_flag = CS_FLAG_VECTOR | cs_flag_primal_vtx;

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
        cs_evaluate_potential_by_qov(dof_flag, def, v_vals, NULL);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        assert(eqp->dof_reduction == CS_PARAM_REDUCTION_DERHAM);
        cs_evaluate_potential_at_vertices_by_analytic(def,
                                                      t_eval,
                                                      n_v_selected,
                                                      selected_lst,
                                                      v_vals);
        break;

      case CS_XDEF_BY_DOF_FUNCTION:
        cs_evaluate_potential_at_vertices_by_dof_func(def,
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

    BFT_FREE(def2v_idx);

  } /* Initial values to set */

  /* Set the boundary values as initial values: Compute the values of the
     Dirichlet BC */

  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   _vvb_cell_builder[0], /* static variable */
                                   eqc->vtx_bc_flag,
                                   v_vals);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector steady-state
 *         convection/diffusion/reaction equation with a CDO-Vb scheme.
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_solve_steady_state(bool                        cur2prev,
                                   const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)context;
  cs_cdo_system_helper_t  *sh = eqb->system_helper;
  cs_field_t  *fld = cs_field_by_id(field_id);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  /* Build an array storing the Dirichlet values at vertices
   * First argument is set to t_cur + dt_cur even if this is a steady
   * computation since one can call this function to compute a steady-state
   * solution at each time step of an unsteady computation.
   */

  _vvb_setup(time_eval, mesh, eqp, eqb, eqc->vtx_bc_flag);

  /* Initialize the local system: rhs, matrix and assembler values */

  double  rhs_norm = 0.;
  cs_real_t  *rhs = NULL;

  cs_cdo_system_helper_init_system(sh, &rhs);

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
    cs_cell_sys_t  *csys = _vvb_cell_system[t_id];
    cs_cell_builder_t  *cb = _vvb_cell_builder[t_id];
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *diff_hodge =
      (eqc->diffusion_hodge == NULL) ? NULL : eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    /* Set times at which one evaluates quantities if needed */

    cb->t_pty_eval = time_eval;
    cb->t_bc_eval = time_eval;
    cb->t_st_eval = time_eval;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, diff_hodge, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */

      _vvb_init_cell_system(cm, eqp, eqb, eqc->vtx_bc_flag, fld->val,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (mass_hodge->matrix) */

      _vvb_conv_diff_reac(eqp, eqb, eqc, cm, mass_hodge, diff_hodge, csys, cb);

      const bool  has_sourceterm = cs_equation_param_has_sourceterm(eqp);

      if (has_sourceterm) { /* SOURCE TERM
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

        for (int v = 0; v < csys->n_dofs; v++)
          csys->rhs[v] += csys->source[v];

      } /* End of source term */

      /* Apply boundary conditions (those which are weakly enforced) */

      _vvb_apply_weak_bc(eqp, eqc, cm, fm, diff_hodge, csys, cb);

      /* Enforce values if needed (internal or Dirichlet) */

      _vvb_enforce_values(eqp, eqb, eqc, cm, fm, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* Compute a norm of the RHS for the normalization of the residual
         of the linear system to solve */

      rhs_norm += _vvb_cw_rhs_normalization(eqp->sles_param->resnorm_type,
                                            cm, csys);

      /* ASSEMBLY PROCESS
       * ================ */

      _vvb_assembly(csys, sh->blocks[0], rhs, eqc, asb);

    } /* Main loop on cells */

  } /* OPENMP Block */

  /* Free temporary buffers and structures */

  cs_cdo_system_helper_finalize_assembly(sh);
  cs_equation_builder_reset(eqb);

  /* Last step in the computation of the renormalization coefficient */

  cs_cdo_solve_sync_rhs_norm(eqp->sles_param->resnorm_type,
                             quant->vol_tot,
                             eqc->n_dofs, /* 3*n_vertices */
                             rhs,
                             &rhs_norm);

  /* Copy current field values to previous values */

  if (cur2prev)
    cs_field_current_to_previous(fld);

  /* End of the system building */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Solve the linear system (treated as a scalar-valued system
     with 3 times more DoFs) */

  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);
  cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);

  if (sh->blocks[0]->info.unrolled) {

    cs_cdo_solve_scalar_system(eqc->n_dofs, /* 3*n_vertices */
                               eqp->sles_param,
                               matrix,
                               range_set,
                               rhs_norm,
                               true, /* rhs_redux */
                               sles,
                               fld->val,
                               rhs);
  }
  else {

    cs_cdo_solve_vector_system(eqc->n_dofs/3, /* = n_vertices */
                               true,          /* interlaced ? */
                               eqp->sles_param,
                               matrix,
                               range_set,
                               rhs_norm,
                               true, /* rhs_redux */
                               sles,
                               fld->val,
                               rhs);

  }

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t1, &t2);

  cs_sles_free(sles);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdovb_vecteq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_update_field(const cs_real_t            *solu,
                             const cs_real_t            *rhs,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data,
                             cs_real_t                  *field_val)
{
  CS_UNUSED(rhs);
  CS_UNUSED(eqp);

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t  *)data;
  cs_timer_t  t0 = cs_timer_time();

  /* Set the computed solution in field array */

# pragma omp parallel for if (eqc->n_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < eqc->n_dofs; i++)
    field_val[i] = solu[i];

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
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
 * \return  a pointer to an array of cs_real_t (size: 3*n_vertices)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_vecteq_get_vertex_values(void      *context,
                                  bool       previous)
{
  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)context;

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
 * \return  a pointer to an array of cs_real_t (size: 3*n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_vecteq_get_cell_values(void      *context,
                                bool       previous)
{
  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)context;

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
    BFT_MALLOC(eqc->cell_values, 3*quant->n_cells, cs_real_t);
  memset(eqc->cell_values, 0, 3*quant->n_cells*sizeof(cs_real_t));

  /* Compute the values at cell centers from an interpolation of the field
     values defined at vertices */

  cs_reco_vect_pv_at_cell_centers(connect->c2v, quant, vtx_values,
                                  eqc->cell_values);

  return eqc->cell_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate a current to previous operation for the field associated to
 *         this equation and potentially for related fields/arrays.
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_current_to_previous(const cs_equation_param_t  *eqp,
                                    cs_equation_builder_t      *eqb,
                                    void                       *context)
{
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)context;
  cs_field_t  *fld = cs_field_by_id(eqc->var_field_id);

  cs_field_current_to_previous(fld);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_extra_post(const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *context)
{
  CS_UNUSED(context);
  CS_UNUSED(eqp);

  const cs_timer_t  t0 = cs_timer_time();

  /* TODO */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
