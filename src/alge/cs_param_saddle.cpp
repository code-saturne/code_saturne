/*============================================================================
 * Routines/structure to handle the settings to solve a saddle-point problem
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

#include "cs_base.h"
#include "cs_log.h"
#include "cs_sles.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_param_saddle.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_param_saddle.c

 * \brief Handle the settings of saddle-point systems.
 *        These systems arise from the monolithic coupling of the Navier-Stokes
 *        equations or in mixed formulation of scalar-valued equations.
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local private variables
 *============================================================================*/

static const char
cs_param_saddle_schur_approx_name[CS_PARAM_SADDLE_N_SCHUR_APPROX]
                                 [CS_BASE_STRING_LEN] = {
  N_("None"),
  N_("Based on the diagonal"),
  N_("Identity matrix"),
  N_("Lumped inverse"),
  N_("Scaled mass matrix"),
  N_("Based on the diagonal + scaled mass matrix"),
  N_("Lumped inverse + scaled mass scaling")
};

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a \ref cs_param_sles_t structure for an extra system used in
 *        the construction of the approximation of the Schur complement
 *
 * \param[in] saddlep   set of parameters managing a saddle-point problem
 *
 * \return the pointer to the newly defined set of parameters
 */
/*----------------------------------------------------------------------------*/

static cs_param_sles_t *
_init_xtra_slesp(const cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return nullptr;

  const cs_param_sles_t  *b11_slesp = saddlep->block11_sles_param;
  const char  *basename = cs_param_saddle_get_name(saddlep);

  int  len = strlen(basename) + strlen("_b11_xtra");
  char  *name = nullptr;
  BFT_MALLOC(name, len + 1, char);
  sprintf(name, "%s_b11_xtra", basename);

  /* One starts from a copy of the (1, 1) block settings */

  cs_param_sles_t  *xtra_slesp = cs_param_sles_create(-1, name);

  cs_param_sles_copy_from(b11_slesp, xtra_slesp);

  /* By default, a coarse approximation should be sufficient */

  xtra_slesp->cvg_param.rtol = 1e-3;
  xtra_slesp->cvg_param.n_max_iter = 50;

  BFT_FREE(name);

  return xtra_slesp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy a \ref cs_param_sles_t structure from an existing one. Case of
 *        the set of parameters associated to the additional SLES system.
 *
 * \param[in] slesp_ref     set of parameters managing a SLES (the ref. one)
 * \param[in] saddlep_dest  set of parameters managing a saddle-point problem
 *
 * \return the pointer to the newly defined set of parameters
 */
/*----------------------------------------------------------------------------*/

static cs_param_sles_t *
_copy_xtra_slesp(const cs_param_sles_t    *slesp_ref,
                 const cs_param_saddle_t  *saddlep_dest)
{
  if (slesp_ref == nullptr)
    return nullptr;

  cs_param_sles_t  *slesp_dest = _init_xtra_slesp(saddlep_dest);

  cs_param_sles_copy_from(slesp_ref, slesp_dest);

  return slesp_dest;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a \ref cs_param_sles_t structure for a system solved during
 *        the first step (transformation or initial resolution). The starting
 *        point is the settings for the (1,1) block.
 *
 * \param[in] saddlep  set of parameters managing a saddle-point problem
 *
 * \return the pointer to the newly defined set of parameters
 */
/*----------------------------------------------------------------------------*/

static cs_param_sles_t *
_init_init_slesp(const cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return nullptr;

  const cs_param_sles_t  *b11_slesp = saddlep->block11_sles_param;
  const char  *basename = saddlep->block11_sles_param->name;

  int  len = strlen(basename) + strlen("_init");
  char  *name = nullptr;
  BFT_MALLOC(name, len + 1, char);
  sprintf(name, "%s_init", basename);

  /* One starts from a copy of the (1, 1) block settings */

  cs_param_sles_t  *init_slesp = cs_param_sles_create(-1, name);

  cs_param_sles_copy_from(b11_slesp, init_slesp);

  /* One needs a more accurate approximation to solve the initial system or to
     transform the system at the initial step */

  double  min_tol_threshold = 1e-14;
  double  tol = fmin(0.1 * b11_slesp->cvg_param.rtol,
                     0.1 * saddlep->cvg_param.rtol);

  tol = fmin(tol, 10*saddlep->cvg_param.atol);

  /* Avoid a too small tolerance if algo->cvg_param.atol is very small */

  tol = fmax(tol, min_tol_threshold);

  init_slesp->cvg_param.rtol = tol;
  init_slesp->cvg_param.atol = fmin(saddlep->cvg_param.atol,
                                    b11_slesp->cvg_param.atol);

  BFT_FREE(name);

  return init_slesp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy a \ref cs_param_sles_t structure from an existing one. Case of
 *        the set of parameters associated to the initial system.
 *
 * \param[in] slesp_ref     set of parameters managing a SLES (the ref. one)
 * \param[in] saddlep_dest  set of parameters managing a saddle-point problem
 *
 * \return the pointer to the newly defined set of parameters
 */
/*----------------------------------------------------------------------------*/

static cs_param_sles_t *
_copy_init_slesp(const cs_param_sles_t    *slesp_ref,
                 const cs_param_saddle_t  *saddlep_dest)
{
  if (slesp_ref == nullptr)
    return nullptr;

  cs_param_sles_t  *slesp_dest = _init_init_slesp(saddlep_dest);

  cs_param_sles_copy_from(slesp_ref, slesp_dest);

  return slesp_dest;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define a \ref cs_param_sles_t structure for the Schur complement
 *        system
 *
 * \param[in] saddlep   set of parameters managing a saddle-point problem
 */
/*----------------------------------------------------------------------------*/

static void
_init_schur_slesp(cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return;

  if (saddlep->schur_sles_param != nullptr)
    cs_param_sles_free(&(saddlep->schur_sles_param));

  const char  *basename = cs_param_saddle_get_name(saddlep);

  int  len = strlen(basename) + strlen("_schur_approx");
  char  *name = nullptr;
  BFT_MALLOC(name, len + 1, char);
  sprintf(name, "%s_schur_approx", basename);

  cs_param_sles_t  *schurp = cs_param_sles_create(-1, name);

  int ierr = cs_param_sles_set_solver("fcg", schurp);
  if (ierr != EXIT_SUCCESS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Error raised during the setup of the solver dedicated\n"
              "to the Schur complement approximation (%s).",
              __func__, schurp->name);

  ierr = cs_param_sles_set_precond("amg", schurp);
  assert(ierr == EXIT_SUCCESS);

  schurp->cvg_param.rtol = 1e-4; /* relative tolerance to stop the iterative
                                    solver */

  saddlep->schur_sles_param = schurp;

  /* Extra solve needed to compute the Schur approximation with some
     settings */

  cs_param_sles_t  *xtra_slesp = nullptr;
  if (saddlep->schur_approx == CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE ||
      saddlep->schur_approx == CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE)
    xtra_slesp = _init_xtra_slesp(saddlep);

  /* Associate the parameter structure to the context structure */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_GCR:
    {
    cs_param_saddle_context_block_krylov_t *ctxp =
      static_cast<cs_param_saddle_context_block_krylov_t *>(saddlep->context);

    if (ctxp->xtra_sles_param != nullptr)
      cs_param_sles_free(&(ctxp->xtra_sles_param));

    ctxp->xtra_sles_param = xtra_slesp;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
    cs_param_saddle_context_uzacg_t *ctxp =
      static_cast<cs_param_saddle_context_uzacg_t *>(saddlep->context);

    if (ctxp->xtra_sles_param != nullptr)
      cs_param_sles_free(&(ctxp->xtra_sles_param));

    ctxp->xtra_sles_param = xtra_slesp;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_SIMPLE:
    {
    cs_param_saddle_context_simple_t *ctxp =
      static_cast<cs_param_saddle_context_simple_t *>(saddlep->context);

    if (ctxp->xtra_sles_param != nullptr)
      cs_param_sles_free(&(ctxp->xtra_sles_param));

    ctxp->xtra_sles_param = xtra_slesp;
    }
    break;

  default:
    assert(xtra_slesp == nullptr); /* There should be no xtra_slesp */
    break;
  }

  BFT_FREE(name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Release structure associated to a context
 *
 * \param[in, out] saddlep        set of parameters for solving a saddle-point
 */
/*----------------------------------------------------------------------------*/

static void
_free_context(cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return;
  if (saddlep->context == nullptr)
    return;

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
    cs_param_saddle_context_alu_t *ctxp =
      static_cast<cs_param_saddle_context_alu_t *>(saddlep->context);

    cs_param_sles_free(&(ctxp->init_sles_param));
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
    {
    cs_param_saddle_context_block_krylov_t *ctxp =
      static_cast<cs_param_saddle_context_block_krylov_t *>(saddlep->context);

    cs_param_sles_free(&(ctxp->xtra_sles_param));
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    {
    cs_param_saddle_context_gkb_t *ctxp =
      static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);

    cs_param_sles_free(&(ctxp->init_sles_param));
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
    cs_param_saddle_context_uzacg_t *ctxp =
      static_cast<cs_param_saddle_context_uzacg_t *>(saddlep->context);

    cs_param_sles_free(&(ctxp->init_sles_param));
    cs_param_sles_free(&(ctxp->xtra_sles_param));
    }
    break;

  default:
    /* Nothing else to free */
    break;

  }

  BFT_FREE(saddlep->context);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the number of iterations to store before starting a Krylov solver
 *
 * \param[in, out] saddlep        set of parameters for solving a saddle-point
 * \param[in]      restart_range  number of directions
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_restart_range(cs_param_saddle_t  *saddlep,
                                  int                 restart_range)
{
  if (saddlep == nullptr)
    return;

  switch (saddlep->solver) {
  case CS_PARAM_SADDLE_SOLVER_FGMRES:
  case CS_PARAM_SADDLE_SOLVER_GCR:
    {
    cs_param_saddle_context_block_krylov_t *ctxp =
      static_cast<cs_param_saddle_context_block_krylov_t *>(saddlep->context);

    ctxp->n_stored_directions = restart_range;
    }
    break;

  default: /* Not useful */
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS,
                  "%s: Restart range not taken into account.\n"
                  "%s: Change the saddle-point solver to do that.",
                  __func__, __func__);
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the scaling coefficient used in the Notay's transformation
 *        devised in
 *        "Algebraic multigrid for Stokes equations" SIAM J. Sci. Comput.
 *        Vol. 39 (5), 2017
 *        In this article, this scaling is denoted by alpha
 *
 * \param[in, out] saddlep       set of parameters for solving a saddle-point
 * \param[in]      scaling_coef  value of the scaling coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_notay_scaling(cs_param_saddle_t  *saddlep,
                                  double              scaling_coef)
{
  if (saddlep == nullptr)
    return;

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM)
    return;

  cs_param_saddle_context_notay_t *ctx =
    static_cast<cs_param_saddle_context_notay_t *>(saddlep->context);

  ctx->scaling_coef = scaling_coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the scaling in front of the augmentation term when an ALU or a
 *        GKB algorithm is considered
 *
 * \param[in, out] saddlep  set of parameters for solving a saddle-point
 * \param[in]      coef     value of the scaling coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_augmentation_coef(cs_param_saddle_t  *saddlep,
                                      double              coef)
{
  if (saddlep == nullptr)
    return;

  switch (saddlep->solver) {
  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
    cs_param_saddle_context_alu_t *ctx =
      static_cast<cs_param_saddle_context_alu_t *>(saddlep->context);

    ctx->augmentation_scaling = coef;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    {
    cs_param_saddle_context_gkb_t *ctx =
      static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);

    ctx->augmentation_scaling = coef;
    }
    break;

  default: /* Not useful */
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS,
                  "%s: Augmentation coef. not taken into account.\n"
                  "%s: Saddle-point solver not relevant.",
                  __func__, __func__);
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the scaling coefficient in front of the augmentation term when an
 *        ALU or a GKB algorithm is considered.
 *
 * \param[in] saddlep  set of parameters for solving a saddle-point
 *
 * \return 0 if not relevant or the value of the augmentation coefficient
 */
/*----------------------------------------------------------------------------*/

double
cs_param_saddle_get_augmentation_coef(const cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return 0;

  switch (saddlep->solver) {
  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
    cs_param_saddle_context_alu_t *ctxp =
      static_cast<cs_param_saddle_context_alu_t *>(saddlep->context);

    return ctxp->augmentation_scaling;
    }
  case CS_PARAM_SADDLE_SOLVER_GKB:
    {
    cs_param_saddle_context_gkb_t *ctxp =
      static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);

    return ctxp->augmentation_scaling;
    }
  default: /* Not useful */
    return 0.;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the name of the type of saddle-point solver
 *
 * \param[in] type  type of saddle-point solver
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_saddle_get_type_name(cs_param_saddle_solver_t  type)
{
  switch (type) {
  case CS_PARAM_SADDLE_SOLVER_NONE:
    return "None";

  case CS_PARAM_SADDLE_SOLVER_ALU:
    return "Augmented-Lagrangian Uzawa";

  case CS_PARAM_SADDLE_SOLVER_FGMRES:
    return "Flexible GMRES";

  case CS_PARAM_SADDLE_SOLVER_GCR:
    return "GCR";

  case CS_PARAM_SADDLE_SOLVER_GKB:
    return "Golub-Kahan Bidiagonalization (GKB)";

  case CS_PARAM_SADDLE_SOLVER_MINRES:
    return "MinRES";

  case CS_PARAM_SADDLE_SOLVER_MUMPS:
    return "MUMPS";

  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM:
    return "Notay's algebraic transformation";

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    return "Uzawa-CG";

  case CS_PARAM_SADDLE_SOLVER_SIMPLE:
    return "Simple";

  case CS_PARAM_SADDLE_N_SOLVERS:
    return "Undefined";
  }

  return "Undefined";
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a cs_param_saddle_t structure
 *        No solver is set by default.
 *
 * \return a pointer to the new cs_param_saddle_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_saddle_t *
cs_param_saddle_create(void)
{
  cs_param_saddle_t  *saddlep = nullptr;

  BFT_MALLOC(saddlep, 1, cs_param_saddle_t);

  saddlep->verbosity = 0;
  saddlep->name = nullptr;

  saddlep->solver_class = CS_PARAM_SOLVER_CLASS_CS;
  saddlep->solver = CS_PARAM_SADDLE_SOLVER_NONE;  /* Not used by default */
  saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;

  saddlep->cvg_param = (cs_param_convergence_t){
    .atol       = 1e-12, /* absolute tolerance */
    .rtol       = 1e-6,  /* relative tolerance */
    .dtol       = 1e3,   /* divergence tolerance */
    .n_max_iter = 100,
  };
  /* saddlep->block11_sles_param is shared and thus is only set if a
     saddle-point problem is solved */

  saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_NONE;
  saddlep->schur_sles_param = nullptr;

  /* By default, no context is set */

  saddlep->context = nullptr;

  return saddlep;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the structure storing the parameter settings for a saddle-point
 *        system
 *
 * \param[in, out] p_saddlep    double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_free(cs_param_saddle_t  **p_saddlep)
{
  if (p_saddlep == nullptr)
    return;

  cs_param_saddle_t  *saddlep = *p_saddlep;

  if (saddlep == nullptr)
    return;

  BFT_FREE(saddlep->name);

  cs_param_sles_free(&(saddlep->schur_sles_param));

  /* block11_sles_param is shared */

  _free_context(saddlep);

  BFT_FREE(saddlep);
  *p_saddlep = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the pointer to the set of parameters to handle a SLES. This SLES
 *        is associated to the approximation of the Schur complement. This is
 *        only useful for solving a saddle-point problem relying on an
 *        elaborated approximation of the Schur complement.
 *
 * \param[in] saddlep  pointer to a \ref cs_param_saddle_t structure
 *
 * \return a pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_param_saddle_get_schur_sles_param(const cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return nullptr;

  return saddlep->schur_sles_param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the pointer to the set of parameters to handle a SLES. This SLES
 *        is associated to an extra-operation specific to a saddle-point solver
 *        It returns a non nullptr pointer only for some sadlle-point solver
 *        relying on a more elaborated Schur complement approximation.
 *
 * \param[in] saddlep  pointer to a \ref cs_param_saddle_t structure
 *
 * \return a pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_param_saddle_get_xtra_sles_param(const cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return nullptr;

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_GCR:
    {
    cs_param_saddle_context_block_krylov_t *ctxp =
      static_cast<cs_param_saddle_context_block_krylov_t *>(saddlep->context);

    return ctxp->xtra_sles_param;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
    cs_param_saddle_context_uzacg_t *ctxp =
      static_cast<cs_param_saddle_context_uzacg_t *>(saddlep->context);

    return ctxp->xtra_sles_param;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_SIMPLE:
    {
    cs_param_saddle_context_simple_t *ctxp =
      static_cast<cs_param_saddle_context_simple_t *>(saddlep->context);

    return ctxp->xtra_sles_param;
    }
    break;

  default:
    return nullptr;
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the pointer to the set of parameters to handle a SLES. This SLES
 *        is associated to the initial saddle-point problem. It returns a non
 *        nullptr pointer only for some sadlle-point solver.
 *
 * \param[in] saddlep  pointer to a \ref cs_param_saddle_t structure
 *
 * \return a pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_param_saddle_get_init_sles_param(const cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return nullptr;

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
    cs_param_saddle_context_alu_t *ctxp =
      static_cast<cs_param_saddle_context_alu_t *>(saddlep->context);

    return ctxp->init_sles_param;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    {
    cs_param_saddle_context_gkb_t *ctxp =
      static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);

    return ctxp->init_sles_param;
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
    cs_param_saddle_context_uzacg_t *ctxp =
      static_cast<cs_param_saddle_context_uzacg_t *>(saddlep->context);

    return ctxp->init_sles_param;
    }
    break;

  default:
    return nullptr;
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Retrieve the name of the saddle-point solver
 *
 * \param[in] saddlep  pointer to a set of saddle-point parameters
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

const char *
cs_param_saddle_get_name(const cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return "Undefined";

  if (saddlep->name != nullptr)
    return saddlep->name;

  else {

    if (saddlep->block11_sles_param != nullptr)
      return saddlep->block11_sles_param->name;
    else
      return "Undefined";

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the name of the saddle-point system.
 *
 * \param[in]      name     name associated to this saddle-point system
 * \param[in, out] saddlep  pointer to the structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_name(const char         *name,
                         cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return;

  size_t  len = strlen(name);
  BFT_MALLOC(saddlep->name, len + 1, char);
  strncpy(saddlep->name, name, len + 1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign the \ref cs_param_sles_t structure (shared) related to the
 *        (1,1)-block to the structure managing the resolution of the
 *        saddle-point problems.
 *
 * \param[in, out] saddlep         pointer to the structure to update
 * \param[in]      block11_slesp   set of parameters for the (1,1) block
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_set_block11_sles_param(cs_param_saddle_t      *saddlep,
                                       const cs_param_sles_t  *block11_slesp)
{
  if (saddlep == nullptr)
    return;

  saddlep->block11_sles_param = block11_slesp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the type of preconditioning to apply for this saddle-point system
 *
 * \param[in]      keyval     value of the key for the preconditioner
 * \param[in, out] saddlep    pointer to the structure to update
 *
 * \return 0 if no error detected, > 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_param_saddle_set_precond(const char          *keyval,
                            cs_param_saddle_t   *saddlep)
{
  int  ierr = 1;
  if (keyval == nullptr)
    return ierr;
  if (saddlep == nullptr)
    return ierr;

  if (strcmp(keyval, "none") == 0)
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;
  else if (strcmp(keyval, "diag") == 0)
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_DIAG;
  else if (strcmp(keyval, "lower") == 0)
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_LOWER;
  else if (strcmp(keyval, "sgs") == 0)
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_SGS;
  else if (strcmp(keyval, "upper") == 0)
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_UPPER;
  else if (strcmp(keyval, "uzawa") == 0) {

    saddlep->precond = CS_PARAM_SADDLE_PRECOND_UZAWA;

    /* One expects a Schur approximation in this case */

    if (saddlep->schur_approx == CS_PARAM_SADDLE_SCHUR_NONE)
      saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_MASS_SCALED;

  }
  else
    return ierr;

  return 0; /* No error detected */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the type of Schur approximation to apply to this saddle-point
 *        system
 *
 * \param[in]      keyval     value of the key for the schur approx.
 * \param[in, out] saddlep    pointer to the structure to update
 *
 * \return 0 if no error detected, > 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_param_saddle_set_schur_approx(const char          *keyval,
                                 cs_param_saddle_t   *saddlep)
{
  int  ierr = 1;
  if (keyval == nullptr)
    return ierr;
  if (saddlep == nullptr)
    return ierr;

  if (strcmp(keyval, "none") == 0)
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_NONE;

  else if (strcmp(keyval, "diag_inv") == 0) {

    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE;
    _init_schur_slesp(saddlep);

  }
  else if (strcmp(keyval, "identity") == 0)
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_IDENTITY;

  else if (strcmp(keyval, "lumped_inv") == 0) {

    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE;
    _init_schur_slesp(saddlep);

  }
  else if (strcmp(keyval, "mass") == 0 ||
           strcmp(keyval, "mass_scaled") == 0)

    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_MASS_SCALED;

  else if (strcmp(keyval, "mass_scaled_diag_inv") == 0) {

    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE;
    _init_schur_slesp(saddlep);

  }
  else if (strcmp(keyval, "mass_scaled_lumped_inv") == 0) {

    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE;
    _init_schur_slesp(saddlep);

  }
  else
    return ierr;

  return 0; /* No error detected */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the class of solver to apply for this saddle-point system
 *
 * \param[in]      keyval     value of the key for the preconditioner
 * \param[in, out] saddlep    pointer to the structure to update
 *
 * \return 0 if no error detected, > 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_param_saddle_set_solver_class(const char          *keyval,
                                 cs_param_saddle_t   *saddlep)
{
  int  ierr = 1;
  if (keyval == nullptr)
    return ierr;
  if (saddlep == nullptr)
    return ierr;

  if (strcmp(keyval, "cs") == 0 ||
      strcmp(keyval, "saturne") == 0)
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_CS;

  else if (strcmp(keyval, "petsc") == 0) {

    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_PETSC;

    cs_param_solver_class_t  ret_class =
      cs_param_sles_check_class(CS_PARAM_SOLVER_CLASS_PETSC);
    if (ret_class != CS_PARAM_SOLVER_CLASS_PETSC)
      return 2;

  }
  else if (strcmp(keyval, "mumps") == 0) {

    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_MUMPS;

    cs_param_solver_class_t  ret_class =
      cs_param_sles_check_class(CS_PARAM_SOLVER_CLASS_MUMPS);

    if (ret_class == CS_PARAM_N_SOLVER_CLASSES)
      return 3;

  }
  else
    return ierr;

  return 0; /* No error detected */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the type of solver to apply for this saddle-point system
 *
 * \param[in]      keyval   value of the key for the preconditioner
 * \param[in, out] saddlep  pointer to the structure to update
 *
 * \return 0 if no error detected, > 0 otherwise
 */
/*----------------------------------------------------------------------------*/

int
cs_param_saddle_set_solver(const char          *keyval,
                           cs_param_saddle_t   *saddlep)
{
  int  ierr = 1;
  if (keyval == nullptr)
    return ierr;
  if (saddlep == nullptr)
    return ierr;

  _free_context(saddlep); /* If one switches from one solver to another one */

  if (strcmp(keyval, "none") == 0)
    saddlep->solver = CS_PARAM_SADDLE_SOLVER_NONE;

  else if (strcmp(keyval, "alu") == 0) {

    saddlep->solver = CS_PARAM_SADDLE_SOLVER_ALU;
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_CS;
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_NONE;

    /* Context structure dedicated to this algorithm */

    cs_param_saddle_context_alu_t  *ctxp = nullptr;
    BFT_MALLOC(ctxp, 1, cs_param_saddle_context_alu_t);

    /* Default value for this context */

    ctxp->augmentation_scaling = 100;
    ctxp->dedicated_init_sles = false;

    /* Initialize an additional set of SLES parameters for the transformation
       of the system (this is different from defining a dedicated cs_sles_t
       structure). */

    ctxp->init_sles_param = _init_init_slesp(saddlep);

    saddlep->context = ctxp;

  }
  else if (strcmp(keyval, "fgmres") == 0) {

    saddlep->solver = CS_PARAM_SADDLE_SOLVER_FGMRES;
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_PETSC;

    cs_param_solver_class_t  ret_class =
      cs_param_sles_check_class(CS_PARAM_SOLVER_CLASS_PETSC);
    if (ret_class != CS_PARAM_SOLVER_CLASS_PETSC)
      return 2;

    cs_param_saddle_context_block_krylov_t  *ctxp = nullptr;
    BFT_MALLOC(ctxp, 1, cs_param_saddle_context_block_krylov_t);

    ctxp->n_stored_directions = 30;    /* default value */
    ctxp->xtra_sles_param = nullptr;      /* It should remain to nullptr */

    saddlep->context = ctxp;

  }
  else if (strcmp(keyval, "gcr") == 0) {

    saddlep->solver = CS_PARAM_SADDLE_SOLVER_GCR;
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_CS;

    /* Additional default settings with this choice of solver */

    saddlep->precond = CS_PARAM_SADDLE_PRECOND_UPPER;
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_MASS_SCALED;

    cs_param_saddle_context_block_krylov_t  *ctxp = nullptr;
    BFT_MALLOC(ctxp, 1, cs_param_saddle_context_block_krylov_t);

    ctxp->n_stored_directions = 30;    /* default value */
    ctxp->xtra_sles_param = nullptr;      /* It depends on the type of Schur
                                          approximation */

    saddlep->context = ctxp;

  }
  else if (strcmp(keyval, "gkb") == 0) {

    saddlep->solver = CS_PARAM_SADDLE_SOLVER_GKB;
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_CS; /* by default */
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_NONE;

    /* Context structure dedicated to this algorithm */

    cs_param_saddle_context_gkb_t  *ctxp = nullptr;
    BFT_MALLOC(ctxp, 1, cs_param_saddle_context_gkb_t);

    ctxp->augmentation_scaling = 0;  /* default value */
    ctxp->truncation_threshold = 5;  /* default value */
    ctxp->dedicated_init_sles = false;

    /* Initialize an additional set of SLES parameters for the transformation
       of the system (this is different from defining a dedicated cs_sles_t
       structure). */

    ctxp->init_sles_param = _init_init_slesp(saddlep);

    saddlep->context = ctxp;

  }
  else if (strcmp(keyval, "minres") == 0) {

    saddlep->solver = CS_PARAM_SADDLE_SOLVER_MINRES;
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_CS;

    /* Additional default settings with this choice of solver */

    saddlep->precond = CS_PARAM_SADDLE_PRECOND_DIAG;
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_MASS_SCALED;

    /* This should be a used for solving a Stokes problem. Thus, there is no
       need for an elaborated Schur approximation (a scaled mass matrix is
       sufficient) */

  }
  else if (strcmp(keyval, "mumps") == 0) {

    saddlep->solver = CS_PARAM_SADDLE_SOLVER_MUMPS;
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_MUMPS;
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_NONE;

    cs_param_solver_class_t  ret_class =
      cs_param_sles_check_class(CS_PARAM_SOLVER_CLASS_MUMPS);

    if (ret_class == CS_PARAM_N_SOLVER_CLASSES)
      return 3;

  }
  else if (strcmp(keyval, "notay") == 0) {

    saddlep->solver = CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM;
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_PETSC;
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_NONE;

    /* Context structure dedicated to this algorithm */

    cs_param_saddle_context_notay_t  *ctxp = nullptr;
    BFT_MALLOC(ctxp, 1, cs_param_saddle_context_notay_t);

    ctxp->scaling_coef = 1.0;  /* default value */

    saddlep->context = ctxp;

  }
  else if (strcmp(keyval, "uzawa_cg") == 0) {

    saddlep->solver = CS_PARAM_SADDLE_SOLVER_UZAWA_CG;
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_CS;
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_MASS_SCALED;

    /* Context structure dedicated to this algorithm */

    cs_param_saddle_context_uzacg_t  *ctxp = nullptr;
    BFT_MALLOC(ctxp, 1, cs_param_saddle_context_uzacg_t);

    ctxp->xtra_sles_param = nullptr;  /* It depends on the type of Schur
                                      approximation used */

    ctxp->dedicated_init_sles = false;

    /* Initialize an additional set of SLES parameters for the initial
       transformation of the system (this is different from defining a
       dedicated cs_sles_t structure). The same SLES can be shared but with
       different settings w.r.t. the stopping convergence criteria. */

    ctxp->init_sles_param = _init_init_slesp(saddlep);

    cs_sles_set_epzero(1e-15);  /* Avoid a too early exit */

    saddlep->context = ctxp;

  }
  else if (strcmp(keyval, "simple") == 0) {

    saddlep->solver = CS_PARAM_SADDLE_SOLVER_SIMPLE;
    saddlep->solver_class = CS_PARAM_SOLVER_CLASS_CS;
    saddlep->precond = CS_PARAM_SADDLE_PRECOND_NONE;
    saddlep->schur_approx = CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE;

    /* Context structure dedicated to this algorithm */

    cs_param_saddle_context_simple_t *ctxp = nullptr;
    BFT_MALLOC(ctxp, 1, cs_param_saddle_context_simple_t);

    ctxp->xtra_sles_param = nullptr;  /* It depends on the type of Schur
                                      approximation used */

    ctxp->dedicated_init_sles = false;

    /* Initialize an additional set of SLES parameters for the initial
       transformation of the system (this is different from defining a
       dedicated cs_sles_t structure). The same SLES can be shared but with
       different settings w.r.t. the stopping convergence criteria. */

    ctxp->init_sles_param = _init_init_slesp(saddlep);

    cs_sles_set_epzero(1e-15);  /* Avoid a too early exit */

    saddlep->context = ctxp;
  }

  else
    return ierr;

  return 0; /* No error detected */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize a \ref cs_param_sles_t structure for the Schur
 *        approximation nested inside a \ref cs_param_saddle_t structure. By
 *        default, this member is not allocated. Do nothing if the related
 *        structure is already allocated.
 *
 * \param[in, out] saddlep    pointer to the structure to update
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_try_init_schur_sles_param(cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return;

  if (saddlep->schur_sles_param != nullptr)
    return; /* Initialization has already been performed */

  _init_schur_slesp(saddlep);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy a cs_param_saddle_t structure from ref to dest
 *
 * \param[in]      ref     reference structure to be copied
 * \param[in, out] dest    destination structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_copy(const cs_param_saddle_t  *ref,
                     cs_param_saddle_t        *dest)
{
  if (ref == nullptr)
    return;

  dest->solver_class = ref->solver_class;
  dest->solver = ref->solver;
  dest->precond = ref->precond;
  dest->schur_approx = ref->schur_approx;

  dest->cvg_param.rtol = ref->cvg_param.rtol;
  dest->cvg_param.atol = ref->cvg_param.atol;
  dest->cvg_param.dtol = ref->cvg_param.dtol;
  dest->cvg_param.n_max_iter = ref->cvg_param.n_max_iter;

  /* Shared pointer */

  dest->block11_sles_param = ref->block11_sles_param;

  if (ref->schur_sles_param != nullptr) {

    if (dest->name == nullptr)
      cs_param_saddle_set_name("auto_schur", dest); /* Avoid using the same
                                                       name */

    cs_param_saddle_try_init_schur_sles_param(dest);

    cs_param_sles_copy_from(ref->schur_sles_param, dest->schur_sles_param);

  }

  /* Copy the associated context */

  switch (ref->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
    cs_param_saddle_context_alu_t *ctxp_ref =
      static_cast<cs_param_saddle_context_alu_t *>(ref->context);
    cs_param_saddle_context_alu_t *ctxp_dest = nullptr;

    BFT_MALLOC(ctxp_dest, 1, cs_param_saddle_context_alu_t);

    ctxp_dest->augmentation_scaling = ctxp_ref->augmentation_scaling;
    ctxp_dest->dedicated_init_sles  = ctxp_ref->dedicated_init_sles;

    ctxp_dest->init_sles_param =
      _copy_init_slesp(ctxp_ref->init_sles_param, dest);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_FGMRES:
  case CS_PARAM_SADDLE_SOLVER_GCR:
    {
    cs_param_saddle_context_block_krylov_t *ctxp_ref =
      static_cast<cs_param_saddle_context_block_krylov_t *>(ref->context);
    cs_param_saddle_context_block_krylov_t *ctxp_dest = nullptr;

    BFT_MALLOC(ctxp_dest, 1, cs_param_saddle_context_block_krylov_t);

    ctxp_ref->n_stored_directions = ctxp_ref->n_stored_directions;

    ctxp_dest->xtra_sles_param =
      _copy_xtra_slesp(ctxp_ref->xtra_sles_param, dest);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    {
    cs_param_saddle_context_gkb_t *ctxp_ref =
      static_cast<cs_param_saddle_context_gkb_t *>(ref->context);
    cs_param_saddle_context_gkb_t *ctxp_dest = nullptr;

    BFT_MALLOC(ctxp_dest, 1, cs_param_saddle_context_gkb_t);

    ctxp_dest->augmentation_scaling = ctxp_ref->augmentation_scaling;
    ctxp_dest->truncation_threshold = ctxp_ref->truncation_threshold;
    ctxp_dest->dedicated_init_sles  = ctxp_ref->dedicated_init_sles;

    ctxp_dest->init_sles_param =
      _copy_init_slesp(ctxp_ref->init_sles_param, dest);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
    cs_param_saddle_context_uzacg_t *ctxp_ref =
      static_cast<cs_param_saddle_context_uzacg_t *>(ref->context);
    cs_param_saddle_context_uzacg_t *ctxp_dest = nullptr;

    ctxp_dest->dedicated_init_sles = ctxp_ref->dedicated_init_sles;

    ctxp_dest->init_sles_param =
      _copy_init_slesp(ctxp_ref->init_sles_param, dest);

    ctxp_dest->xtra_sles_param =
      _copy_xtra_slesp(ctxp_ref->xtra_sles_param, dest);
    }
    break;

  default:
    /* Nothing else to free */
    break;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log the setup information for the given cs_param_saddle_t structure
 *
 * \param[in] saddlep     pointer to the structure
 */
/*----------------------------------------------------------------------------*/

void
cs_param_saddle_log(const cs_param_saddle_t  *saddlep)
{
  if (saddlep == nullptr)
    return;

  if (saddlep->solver == CS_PARAM_SADDLE_SOLVER_NONE)
    return;

  const char  *basename = cs_param_saddle_get_name(saddlep);

  char  *prefix = nullptr;
  int  len = strlen(basename) + strlen("  *  |") + 1;
  BFT_MALLOC(prefix, len, char);
  sprintf(prefix, "  * %s |", basename);

  /* Start the logging */

  cs_log_printf(CS_LOG_SETUP,
                "\n### Setup for the saddle-point system: \"%s\"\n", basename);
  cs_log_printf(CS_LOG_SETUP, "%s Verbosity: %d\n", prefix, saddlep->verbosity);

  /* Solver */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
    cs_param_saddle_context_alu_t *ctxp =
      static_cast<cs_param_saddle_context_alu_t *>(saddlep->context);

    cs_log_printf(CS_LOG_SETUP,
                  "%s Solver: Augmented Lagrangian-Uzawa (ALU)\n",
                  prefix);
    cs_log_printf(CS_LOG_SETUP,
                  "%s ALU parameters: gamma=%5.2e\n",
                  prefix,
                  ctxp->augmentation_scaling);
    cs_log_printf(CS_LOG_SETUP,
                  "%s ALU parameters: dedicated_init_sles=%s\n",
                  prefix,
                  cs_base_strtf(ctxp->dedicated_init_sles));
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_FGMRES:
    {
    cs_param_saddle_context_block_krylov_t *ctxp =
      static_cast<cs_param_saddle_context_block_krylov_t *>(saddlep->context);

    cs_log_printf(CS_LOG_SETUP, "%s Solver: Flexible GMRES (FGMRES)\n", prefix);
    cs_log_printf(CS_LOG_SETUP,
                  "%s FGMRES parameters: n_stored_directions=%d\n",
                  prefix,
                  ctxp->n_stored_directions);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
    {
    cs_param_saddle_context_block_krylov_t *ctxp =
      static_cast<cs_param_saddle_context_block_krylov_t *>(saddlep->context);

    cs_log_printf(CS_LOG_SETUP,
                  "%s Solver: Generalized Conjugate Residual (GCR)\n",
                  prefix);
    cs_log_printf(CS_LOG_SETUP,
                  "%s GCR parameters: n_stored_directions=%d\n",
                  prefix,
                  ctxp->n_stored_directions);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    {
    cs_param_saddle_context_gkb_t *ctxp =
      static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);

    cs_log_printf(CS_LOG_SETUP,
                  "%s Solver: Golub-Kahan Bidiagonalization (GKB)\n",
                  prefix);
    cs_log_printf(CS_LOG_SETUP,
                  "%s GKB parameters:"
                  " gamma=%5.2e; trunctation_threshold=%d\n",
                  prefix,
                  ctxp->augmentation_scaling,
                  ctxp->truncation_threshold);
    cs_log_printf(CS_LOG_SETUP,
                  "%s GKB parameters: dedicated_init_sles=%s\n",
                  prefix,
                  cs_base_strtf(ctxp->dedicated_init_sles));
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_MINRES:
    cs_log_printf(CS_LOG_SETUP, "%s Solver: MINRES\n", prefix);
    break;

  case CS_PARAM_SADDLE_SOLVER_MUMPS:
    cs_log_printf(CS_LOG_SETUP, "%s Solver: MUMPS\n", prefix);
    break;

  case CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM:
    {
    cs_param_saddle_context_notay_t *ctxp =
      static_cast<cs_param_saddle_context_notay_t *>(saddlep->context);

    cs_log_printf(CS_LOG_SETUP, "%s Solver: Notay's transformation\n", prefix);
    cs_log_printf(CS_LOG_SETUP,
                  "%s Notay parameters: alpha=%5.3e\n",
                  prefix,
                  ctxp->scaling_coef);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
    cs_param_saddle_context_uzacg_t *ctxp =
      static_cast<cs_param_saddle_context_uzacg_t *>(saddlep->context);

    cs_log_printf(CS_LOG_SETUP, "%s Solver: Uzawa-CG\n", prefix);

    cs_log_printf(CS_LOG_SETUP,
                  "%s Uzawa-CG parameters: dedicated_init_sles=%s\n",
                  prefix,
                  cs_base_strtf(ctxp->dedicated_init_sles));
    }
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "%s Solver: Undefined\n", prefix);
    break;

  } /* Solver */

  /* Preconditioner */

  switch (saddlep->precond) {

  case CS_PARAM_SADDLE_PRECOND_NONE:
    cs_log_printf(CS_LOG_SETUP, "%s Precond: None\n", prefix);
    break;

  case CS_PARAM_SADDLE_PRECOND_DIAG:
    cs_log_printf(CS_LOG_SETUP, "%s Precond: Diagonal blocks\n", prefix);
    break;

  case CS_PARAM_SADDLE_PRECOND_LOWER:
    cs_log_printf(CS_LOG_SETUP, "%s Precond: Lower triangular blocks\n",
                  prefix);
    break;

  case CS_PARAM_SADDLE_PRECOND_SGS:
    cs_log_printf(CS_LOG_SETUP, "%s Precond: Symm. Gauss-Seidel by block\n",
                  prefix);
    break;

  case CS_PARAM_SADDLE_PRECOND_UPPER:
    cs_log_printf(CS_LOG_SETUP, "%s Precond: Upper triangular blocks\n",
                  prefix);
    break;

  case CS_PARAM_SADDLE_PRECOND_UZAWA:
    cs_log_printf(CS_LOG_SETUP, "%s Precond: Uzawa-like\n", prefix);
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "%s Precond: Undefined\n", prefix);
    break;

  } /* Preconditioner */

  /* Convergence criteria */

  if (saddlep->solver != CS_PARAM_SADDLE_SOLVER_MUMPS &&
      saddlep->solver != CS_PARAM_SADDLE_SOLVER_NOTAY_TRANSFORM) {

    cs_log_printf(CS_LOG_SETUP, "%s Convergence.max_iter:  %d\n",
                  prefix, saddlep->cvg_param.n_max_iter);
    cs_log_printf(CS_LOG_SETUP, "%s Convergence.rtol:     % -10.6e\n",
                  prefix, saddlep->cvg_param.rtol);
    cs_log_printf(CS_LOG_SETUP, "%s Convergence.atol:     % -10.6e\n",
                  prefix, saddlep->cvg_param.atol);

  }

  /* Schur complement */
  /* ---------------- */

  /* Type of approximation */

  switch (saddlep->schur_approx) {

  case CS_PARAM_SADDLE_SCHUR_NONE:
  case CS_PARAM_SADDLE_SCHUR_IDENTITY:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED:
    cs_log_printf(CS_LOG_SETUP, "%s Schur approx.: %s.\n", prefix,
                  cs_param_saddle_schur_approx_name[saddlep->schur_approx]);
    break;

  case CS_PARAM_SADDLE_SCHUR_DIAG_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_LUMPED_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_DIAG_INVERSE:
  case CS_PARAM_SADDLE_SCHUR_MASS_SCALED_LUMPED_INVERSE:
    cs_log_printf(CS_LOG_SETUP, "%s Schur approx.: %s.\n", prefix,
                  cs_param_saddle_schur_approx_name[saddlep->schur_approx]);

    /* SLES parameters associated to the Schur complement */

    cs_param_sles_log(saddlep->schur_sles_param);
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "%s Schur approx.: Undefined\n", prefix);
    break;

  }

  /* Additional SLES parameters for some settings */
  /* -------------------------------------------- */

  switch (saddlep->solver) {

  case CS_PARAM_SADDLE_SOLVER_ALU:
    {
    cs_param_saddle_context_alu_t *ctxp =
      static_cast<cs_param_saddle_context_alu_t *>(saddlep->context);

    if (ctxp->dedicated_init_sles)
      cs_param_sles_log(ctxp->init_sles_param);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GCR:
    {
    cs_param_saddle_context_block_krylov_t *ctxp =
      static_cast<cs_param_saddle_context_block_krylov_t *>(saddlep->context);

    cs_param_sles_log(ctxp->xtra_sles_param);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_GKB:
    {
    cs_param_saddle_context_gkb_t *ctxp =
      static_cast<cs_param_saddle_context_gkb_t *>(saddlep->context);

    if (ctxp->dedicated_init_sles)
      cs_param_sles_log(ctxp->init_sles_param);
    }
    break;

  case CS_PARAM_SADDLE_SOLVER_UZAWA_CG:
    {
    cs_param_saddle_context_uzacg_t *ctxp =
      static_cast<cs_param_saddle_context_uzacg_t *>(saddlep->context);

    if (ctxp->dedicated_init_sles)
      cs_param_sles_log(ctxp->init_sles_param);

    cs_param_sles_log(ctxp->xtra_sles_param);
    }
    break;

  default:
    break; /* Do nothing */

  }

  BFT_FREE(prefix);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
