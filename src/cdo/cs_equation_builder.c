/*============================================================================
 * Functions to handle the equation builder structure
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
#include <float.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_log.h"
#include "cs_parall.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_builder.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions and macros
 *============================================================================*/

#define CS_EQUATION_BUILDER_DBG        0 /* Debug level */

/*============================================================================
 * Local private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a new structure to handle the building of algebraic system
 *         related to a cs_equation_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_builder_t *
cs_equation_builder_create(const cs_equation_param_t   *eqp,
                           const cs_mesh_t             *mesh)
{
  cs_equation_builder_t  *eqb = NULL;

  BFT_MALLOC(eqb, 1, cs_equation_builder_t);

  eqb->init_step = true;

  /* Initialize flags used to knows what kind of cell quantities to build */

  eqb->msh_flag = 0;
  eqb->bd_msh_flag = 0;
  eqb->st_msh_flag = 0;
  if (eqp->dim > 1)
    eqb->sys_flag = CS_FLAG_SYS_VECTOR;
  else
    eqb->sys_flag = 0;

  /* Handle properties */

  eqb->diff_pty_uniform = true;
  if (cs_equation_param_has_diffusion(eqp))
    eqb->diff_pty_uniform = cs_property_is_uniform(eqp->diffusion_property);

  eqb->curlcurl_pty_uniform = true;
  if (cs_equation_param_has_curlcurl(eqp))
    eqb->curlcurl_pty_uniform = cs_property_is_uniform(eqp->curlcurl_property);

  eqb->graddiv_pty_uniform = true;
  if (cs_equation_param_has_graddiv(eqp))
    eqb->graddiv_pty_uniform = cs_property_is_uniform(eqp->graddiv_property);

  eqb->time_pty_uniform = true;
  if (cs_equation_param_has_time(eqp))
    eqb->time_pty_uniform = cs_property_is_uniform(eqp->time_property);

  if (eqp->n_reaction_terms > CS_CDO_N_MAX_REACTIONS)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Number of reaction terms for an equation is too high.\n"
              " Current value: %d (max: %d)\n"
              " Change the value of CS_CDO_N_MAX_REACTIONS in the code or\n"
              " modify your settings or contact the developpement team.",
              __func__, eqp->n_reaction_terms, CS_CDO_N_MAX_REACTIONS);

  for (int i = 0; i < eqp->n_reaction_terms; i++)
    eqb->reac_pty_uniform[i]
      = cs_property_is_uniform(eqp->reaction_properties[i]);

  /* Handle source terms */

  eqb->source_mask = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    /* Default initialization */

    eqb->st_msh_flag = cs_source_term_init(eqp->space_scheme,
                                           eqp->n_source_terms,
                       (cs_xdef_t *const *)eqp->source_terms,
                                           eqb->compute_source,
                                           &(eqb->sys_flag),
                                           &(eqb->source_mask));

  } /* There is at least one source term */

  /* The helper structure is allocated during the initialization of the context
     structure */

  eqb->system_helper = NULL;

  /* Enforcement of DoFs */

  eqb->enforced_values = NULL;

  /* Incremental algorithm (if a context associated has to be define for the
     incremental algo then this done during the initialization of the equation
     context) */

  eqb->increment = NULL;
  if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE)
    eqb->incremental_algo = cs_iter_algo_create(eqp->incremental_algo_param);
  else
    eqb->incremental_algo = NULL;

  /* Set members and structures related to the management of the BCs
     Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.  */

  eqb->dir_values = NULL;
  eqb->face_bc = cs_cdo_bc_face_define(eqp->default_bc,
                                       true, /* Steady BC up to now */
                                       eqp->dim,
                                       eqp->n_bc_defs,
                                       eqp->bc_defs,
                                       mesh->n_b_faces);

  /* User hook function */

  eqb->hook_context = NULL;
  eqb->hook_function = NULL;

  /* Monitoring */

  CS_TIMER_COUNTER_INIT(eqb->tcb); /* build system */
  CS_TIMER_COUNTER_INIT(eqb->tcs); /* solve system */
  CS_TIMER_COUNTER_INIT(eqb->tce); /* extra operations */

  return eqb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the range set structure associated to a builder structure
 *        for the block defined in block_id in the system helper structure
 *
 * \param[in, out]  builder      pointer to a cs_equation_builder_t
 * \param[in]       block_id     id of the block to consider
 *
 * \return a pointer to a cs_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_matrix_t *
cs_equation_builder_get_matrix(const cs_equation_builder_t  *builder,
                               int                           block_id)
{
  if (builder == NULL)
    return NULL;

  cs_cdo_system_helper_t  *sh = builder->system_helper;

  if (sh == NULL)
    return NULL;

  if (block_id > sh->n_blocks - 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid block_id \"%d\". Number of blocks: %d\n",
              __func__, block_id, sh->n_blocks);

  return cs_cdo_system_get_matrix(sh, block_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the range set structure associated to a builder structure
 *        for the block defined in block_id in the system helper structure
 *
 * \param[in, out]  builder      pointer to a cs_equation_builder_t
 * \param[in]       block_id     id of the block to consider
 *
 * \return a pointer to a cs_range_set structure
 */
/*----------------------------------------------------------------------------*/

const cs_range_set_t *
cs_equation_builder_get_range_set(const cs_equation_builder_t  *builder,
                                  int                           block_id)
{
  if (builder == NULL)
    return NULL;

  cs_cdo_system_helper_t  *sh = builder->system_helper;

  if (sh == NULL)
    return NULL;

  if (block_id > sh->n_blocks - 1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid block_id \"%d\". Number of blocks: %d\n",
              __func__, block_id, sh->n_blocks);

  return cs_cdo_system_get_range_set(sh, block_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_builder_t structure
 *
 * \param[in, out]  p_builder  pointer of pointer to the cs_equation_builder_t
 *                             structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_free(cs_equation_builder_t  **p_builder)
{
  if (p_builder == NULL)
    return;
  if (*p_builder == NULL)
    return;

  cs_equation_builder_t  *eqb = *p_builder;

  cs_equation_builder_reset(eqb);

  if (eqb->source_mask != NULL)
    BFT_FREE(eqb->source_mask);

  cs_cdo_system_helper_free(&eqb->system_helper);

  /* Quantities related to the incremental resolution (may be NULL) */

  BFT_FREE(eqb->increment);

  /* If the context is not NULL, this means that an Anderson algorithm has been
     activated otherwise nothing to do */

  cs_iter_algo_aa_free(eqb->incremental_algo);

  BFT_FREE(eqb->incremental_algo);

  /* Free BC structure */

  eqb->face_bc = cs_cdo_bc_free(eqb->face_bc);

  BFT_FREE(eqb);

  *p_builder = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free some members of a cs_equation_builder_t structure
 *
 * \param[in, out]  eqb   pointer to the cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_reset(cs_equation_builder_t  *eqb)
{
  if (eqb == NULL)
    return;

  eqb->init_step = true;
  BFT_FREE(eqb->enforced_values);
  BFT_FREE(eqb->dir_values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Print a message in the performance output file related to the
 *          monitoring of equation
 *
 * \param[in]  eqp    pointer to a set of equation parameters
 * \param[in]  eqb    pointer to an equation builder  structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_log_performance(const cs_equation_param_t     *eqp,
                                    const cs_equation_builder_t   *eqb)
{
  if (eqb == NULL)
    return;
  if (eqp == NULL)
    return;
  if (eqp->flag & CS_EQUATION_INSIDE_SYSTEM)
    return;

  double t[3] = {eqb->tcb.nsec, eqb->tcs.nsec, eqb->tce.nsec};
  for (int i = 0; i < 3; i++) t[i] *= 1e-9;

  if (eqp->name == NULL)
    cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f %9.3f %9.3f seconds\n",
                  "<CDO/Equation> Runtime", t[0], t[1], t[2]);

  else {

    char *msg = NULL;
    int len = 1 + strlen("<CDO/> Runtime") + strlen(eqp->name);

    BFT_MALLOC(msg, len, char);
    sprintf(msg, "<CDO/%s> Runtime", eqp->name);
    cs_log_printf(CS_LOG_PERFORMANCE, " %-35s %9.3f %9.3f %9.3f seconds\n",
                  msg, t[0], t[1], t[2]);
    BFT_FREE(msg);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize all reaction properties.
 *        This function is shared across all CDO schemes. The cs_cell_builder_t
 *        structure stores the computed property values. If the property is
 *        uniform, a first call to the function
 *        cs_equation_builder_init_properties has to be done before the loop on
 *        cells
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      eqb      pointer to a cs_equation_builder_t structure
 * \param[in]      cm       pointer to a \ref cs_cell_mesh_t structure
 * \param[in, out] cb       pointer to a \ref cs_cell_builder_t structure
 *
 * \return true if the reaction property is not equal to zero
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_builder_set_reaction_pty_cw(const cs_equation_param_t     *eqp,
                                        const cs_equation_builder_t   *eqb,
                                        const cs_cell_mesh_t          *cm,
                                        cs_cell_builder_t             *cb)
{
  assert(cs_equation_param_has_reaction(eqp));

  /* Set the (linear) reaction property */

  cb->rpty_val = 0;
  for (int r = 0; r < eqp->n_reaction_terms; r++)
    if (eqb->reac_pty_uniform[r])
      cb->rpty_val += cb->rpty_vals[r];
    else
      cb->rpty_val += cs_property_value_in_cell(cm,
                                                eqp->reaction_properties[r],
                                                cb->t_pty_eval);

  if (fabs(cb->rpty_val) > 0)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize all properties potentially useful to build the algebraic
 *         system. This function is shared across all CDO schemes.
 *         The \ref cs_cell_builder_t structure stores property values related
 *         to the reaction term, unsteady term and grad-div term.
 *
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in, out] diff_hodge   pointer to the diffusion hodge structure
 * \param[in, out] cb           pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_init_properties(const cs_equation_param_t     *eqp,
                                    const cs_equation_builder_t   *eqb,
                                    cs_hodge_t                    *diff_hodge,
                                    cs_cell_builder_t             *cb)
{
  /* Preparatory step for diffusion term
   * One calls this function with the boundary tag to examine all tests */

  if (diff_hodge != NULL && eqb->diff_pty_uniform)
    cs_hodge_set_property_value(0, /* cell_id */
                                cb->t_pty_eval,
                                CS_FLAG_BOUNDARY_CELL_BY_FACE,
                                diff_hodge);

  /* Preparatory step for the grad-div term */

  if (cs_equation_param_has_graddiv(eqp) && eqb->graddiv_pty_uniform)
    cb->gpty_val = cs_property_get_cell_value(0, cb->t_pty_eval,
                                              eqp->graddiv_property);

  /* Preparatory step for the unsteady term */

  if (cs_equation_param_has_time(eqp) && eqb->time_pty_uniform)
    cb->tpty_val = cs_property_get_cell_value(0, cb->t_pty_eval,
                                              eqp->time_property);

  /* Preparatory step for the reaction term(s) */

  if (cs_equation_param_has_reaction(eqp)) {

    for (int i = 0; i < CS_CDO_N_MAX_REACTIONS; i++) cb->rpty_vals[i] = 1.0;

    for (int r = 0; r < eqp->n_reaction_terms; r++) {
      if (eqb->reac_pty_uniform[r]) {
        cb->rpty_vals[r] =
          cs_property_get_cell_value(0, cb->t_pty_eval,
                                     eqp->reaction_properties[r]);
      }
    } /* Loop on reaction properties */

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Take into account the enforcement of internal DoFs. Apply an
 *          algebraic manipulation. Update members of the cs_cell_sys_t
 *          structure related to the internal enforcement.
 *
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aii  | Aie |     | Aii  |  0  |     |bi|     |bi -Aid.x_enf|
 *          |------------| --> |------------| and |--| --> |-------------|
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aei  | Aee |     |  0   |  Id |     |be|     |   x_enf     |
 *
 * where x_enf is the value of the enforcement for the selected internal DoFs
 *
 * \param[in]       eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_enforce_dofs(const cs_equation_builder_t     *eqb,
                                 cs_cell_builder_t               *cb,
                                 cs_cell_sys_t                   *csys)
{
  /* Enforcement of internal DoFs */

  double  *x_vals = cb->values; /* define with cs_enforcement_dofs_cw() */
  double  *ax = cb->values + csys->n_dofs;

  memset(cb->values, 0, 2*csys->n_dofs*sizeof(double));

  bool do_enforcement = cs_enforcement_dofs_cw(eqb->enforced_values,
                                               csys,
                                               cb->values);

  csys->has_internal_enforcement = do_enforcement;

  if (!do_enforcement)
    return;

  /* Contribution of the DoFs which are enforced */

  cs_sdm_matvec(csys->mat, x_vals, ax);

  /* Second pass: Replace the block of enforced DoFs by a diagonal block */

  for (int i = 0; i < csys->n_dofs; i++) {

    if (csys->dof_is_forced[i]) {

      /* Reset row i */

      memset(csys->mat->val + csys->n_dofs*i, 0, csys->n_dofs*sizeof(double));

      /* Reset column i */

      for (int j = 0; j < csys->n_dofs; j++)
        csys->mat->val[i + csys->n_dofs*j] = 0;
      csys->mat->val[i*(1 + csys->n_dofs)] = 1;

      /* Set the RHS */

      csys->rhs[i] = x_vals[i];

    } /* DoF associated to a Dirichlet BC */
    else
      csys->rhs[i] -= ax[i];  /* Update RHS */

  } /* Loop on degrees of freedom */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Take into account the enforcement of internal DoFs. Case of matrices
 *        defined by blocks. Apply an algebraic manipulation. Update members
 *        of the cs_cell_sys_t structure related to the internal enforcement.
 *
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aii  | Aie |     | Aii  |  0  |     |bi|     |bi -Aid.x_enf|
 *          |------------| --> |------------| and |--| --> |-------------|
 *          |      |     |     |      |     |     |  |     |             |
 *          | Aei  | Aee |     |  0   |  Id |     |be|     |   x_enf     |
 *
 * where x_enf is the value of the enforcement for the selected internal DoFs
 *
 * \param[in]       eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_builder_enforce_block_dofs(const cs_equation_builder_t   *eqb,
                                       cs_cell_builder_t             *cb,
                                       cs_cell_sys_t                 *csys)
{
  /* Enforcement of internal DoFs */

  double  *x_vals = cb->values; /* define with cs_enforcement_dofs_cw() */
  double  *ax = cb->values + csys->n_dofs;

  memset(cb->values, 0, 2*csys->n_dofs*sizeof(double));

  bool do_enforcement = cs_enforcement_dofs_cw(eqb->enforced_values,
                                               csys,
                                               cb->values);

  csys->has_internal_enforcement = do_enforcement;

  if (!do_enforcement)
    return;

  /* Contribution of the DoFs which are enforced */

  cs_sdm_block_matvec(csys->mat, x_vals, ax);

  /* Define the new right-hand side (rhs) */

  for (int i = 0; i < csys->n_dofs; i++) {
    if (csys->dof_is_forced[i])
      csys->rhs[i] = x_vals[i];
    else
      csys->rhs[i] -= ax[i];  /* Update RHS */
  }

  const cs_sdm_block_t  *bd = csys->mat->block_desc;

  /* Second pass: Replace the block of enforced DoFs by a diagonal block */

  int s = 0;
  for (int ii = 0; ii < bd->n_row_blocks; ii++) {

    cs_sdm_t  *db = cs_sdm_get_block(csys->mat, ii, ii);
    const int  bsize = db->n_rows*db->n_cols;

    if (csys->dof_is_forced[s]) {

      /* Identity for the diagonal block */

      memset(db->val, 0, sizeof(cs_real_t)*bsize);
      for (int i = 0; i < db->n_rows; i++) {
        db->val[i*(1+db->n_rows)] = 1;
        assert(csys->dof_is_forced[s+i]);
      }

      /* Reset column and row block jj < ii */

      for (int jj = 0; jj < ii; jj++) {

        cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, ii, jj);
        memset(bij->val, 0, sizeof(cs_real_t)*bsize);

        cs_sdm_t  *bji = cs_sdm_get_block(csys->mat, jj, ii);
        memset(bji->val, 0, sizeof(cs_real_t)*bsize);

      }

      /* Reset column and row block jj < ii */

      for (int jj = ii+1; jj < db->n_rows; jj++) {

        cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, ii, jj);
        memset(bij->val, 0, sizeof(cs_real_t)*bsize);

        cs_sdm_t  *bji = cs_sdm_get_block(csys->mat, jj, ii);
        memset(bji->val, 0, sizeof(cs_real_t)*bsize);

      }

    } /* DoF associated to an enforcement of their values*/

    s += db->n_rows;

  } /* Loop on degrees of freedom */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
