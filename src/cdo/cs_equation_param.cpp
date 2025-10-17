/*============================================================================
 * Functions related to the structure cs_equation_param_t storing the settings
 * related to an equation.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>  /* tolower() function */
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_mem.h"
#include "bft/bft_printf.h"

#include "alge/cs_saddle_solver_setup.h"
#include "alge/cs_param_sles_setup.h"
#include "base/cs_boundary_zone.h"
#include "base/cs_log.h"
#include "base/cs_volume_zone.h"
#include "cdo/cs_cdo_advection.h"
#include "cdo/cs_cdo_bc.h"
#include "cdo/cs_cdocb_priv.h"
#include "cdo/cs_hodge.h"
#include "cdo/cs_source_term.h"
#include "mesh/cs_mesh_location.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cdo/cs_equation_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_equation_param.cpp

  \brief Structure and functions handling the specific settings related
         to a cs_equation_t structure
*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Local private variables
 *============================================================================*/

static const cs_real_t  _weak_pena_bc_coef_by_default = 100.;
static const cs_real_t  _strong_pena_bc_coef_by_default = 1e12;

static const char _err_empty_eqp[] =
  N_(" Stop setting an empty cs_equation_param_t structure.\n"
     " Please check your settings.\n");

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find an existing "initial condition" definition related to the given
 *        zone id or add a new one.
 *
 * \param[in, out] eqp   pointer to a cs_equation_param_t structure
 * \param[in]      z_id  id related to a zone
 *
 * \return the definition id
 */
/*----------------------------------------------------------------------------*/

static int
_find_or_add_ic(cs_equation_param_t *eqp,
                int                  z_id)
{
  int def_id = -1;

  for (int i = 0; i < eqp->n_ic_defs; i++) {

    cs_xdef_t *d = eqp->ic_defs[i];

    if (d->z_id == z_id) {
      def_id = i;
      eqp->ic_defs[i] = cs_xdef_free(d);
      break;
    }

  }

  if (def_id == -1) { // create a new definition

    def_id = eqp->n_ic_defs;
    eqp->n_ic_defs += 1;
    BFT_REALLOC(eqp->ic_defs, eqp->n_ic_defs, cs_xdef_t *);
    eqp->ic_defs[def_id] = nullptr;

  }

  return def_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a parameter attached to a keyname in a \ref cs_equation_param_t
 *        structure
 *
 * \param[in, out] eqp     pointer to a \ref cs_equation_param_t structure
 * \param[in]      key     key related to the member of eq to set
 * \param[in]      keyval  accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

static void
_set_key(cs_equation_param_t *eqp,
         cs_equation_key_t    key,
         const char          *keyval)
{
  int  ierr = 0;

  const char  *eqname = eqp->name;
  const char  emsg[] = " %s: \"%s\" equation\n --> Invalid key value \"%s\" for"
    " the keyword %s.\n";
  const char  wkmsg[] = " %s: \"%s\" equation\n --> Deprecated key \"%s\"\n"
    " --> Please replace it with key \"%s\".\n";

  switch(key) {

  case CS_EQKEY_ADV_CIP_COEF:
    eqp->cip_scaling_coef = atof(keyval);

    break;

  case CS_EQKEY_ADV_EXTRAPOL:
    if (strcmp(keyval, "none") == 0)
      eqp->adv_extrapol = CS_PARAM_ADVECTION_EXTRAPOL_NONE;
    else if (strcmp(keyval, "taylor") == 0)
      eqp->adv_extrapol = CS_PARAM_ADVECTION_EXTRAPOL_TAYLOR_2;
    else if (strcmp(keyval, "adams_bashforth") == 0)
      eqp->adv_extrapol = CS_PARAM_ADVECTION_EXTRAPOL_ADAMS_BASHFORTH_2;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_ADV_EXTRAPOL");
    }
    break;

  case CS_EQKEY_ADV_FORMULATION:
    if (strcmp(keyval, "conservative") == 0)
      eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_CONSERV;
    else if (strcmp(keyval, "non_conservative") == 0)
      eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_NONCONS;
    else if (strcmp(keyval, "skew_symmetric") == 0)
      eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_SKEWSYM;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_ADV_FORMULATION");
    }
    break;

  case CS_EQKEY_ADV_SCHEME:
    if (strcmp(keyval, "upwind") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
    else if (strcmp(keyval, "samarskii") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_SAMARSKII;
    else if (strcmp(keyval, "sg") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_SG;
    else if (strcmp(keyval, "centered") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CENTERED;
    else if (strcmp(keyval, "L0") == 0 || strcmp(keyval, "l0") == 0) {
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_L0;
    }
    else if (strcmp(keyval, "centered_v8") == 0) {
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CENTERED_V8;
    }
    else if (strcmp(keyval, "upwind_v8") == 0) {
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND_V8;
    }
    else if (strcmp(keyval, "mix_centered_upwind") == 0 ||
             strcmp(keyval, "hybrid_centered_upwind") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND;
    else if (strcmp(keyval, "cip") == 0) {

      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CIP;

      /* Automatically switch to a non-conservative formulation */

      eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_NONCONS;

    }
    else if (strcmp(keyval, "cip_cw") == 0) {

      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CIP_CW;

      /* Automatically switch to a non-conservative formulation */

      eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_NONCONS;

    }
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_ADV_SCHEME");
    }
    break;

  case CS_EQKEY_ADV_STRATEGY:
    if (strcmp(keyval, "fully_implicit") == 0 ||
        strcmp(keyval, "implicit") == 0)
      eqp->adv_strategy = CS_PARAM_ADVECTION_IMPLICIT_FULL;
    else if (strcmp(keyval, "implicit_linear") == 0 ||
             strcmp(keyval, "linearized") == 0)
      eqp->adv_strategy = CS_PARAM_ADVECTION_IMPLICIT_LINEARIZED;
    else if (strcmp(keyval, "explicit") == 0)
      eqp->adv_strategy = CS_PARAM_ADVECTION_EXPLICIT;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_ADV_STRATEGY");
    }
    break;

  case CS_EQKEY_ADV_UPWIND_PORTION:
    eqp->upwind_portion = atof(keyval);

    /* Automatic switch to a hybrid upwind/centered scheme for advection */

    eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND;
    break;

  case CS_EQKEY_AMG_TYPE:
    ierr = cs_param_sles_set_amg_type(keyval, eqp->sles_param);

    if (ierr > 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_AMG_TYPE");
    }
    break;

  case CS_EQKEY_BC_ENFORCEMENT:
    if (strcmp(keyval, "algebraic") == 0)
      eqp->default_enforcement = CS_PARAM_BC_ENFORCE_ALGEBRAIC;
    else if (strcmp(keyval, "penalization") == 0)
      eqp->default_enforcement = CS_PARAM_BC_ENFORCE_PENALIZED;
    else if (strcmp(keyval, "weak_sym") == 0)
      eqp->default_enforcement = CS_PARAM_BC_ENFORCE_WEAK_SYM;
    else if (strcmp(keyval, "weak") == 0)
      eqp->default_enforcement = CS_PARAM_BC_ENFORCE_WEAK_NITSCHE;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_BC_ENFORCEMENT");
    }
    break;

  case CS_EQKEY_BC_STRONG_PENA_COEFF:
    eqp->strong_pena_bc_coeff = atof(keyval);
    if (eqp->strong_pena_bc_coeff < 1.)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value of the penalization coefficient %5.3e\n"
                " This should be positive and large.\n"
                " Equation: %s\n",
                __func__, eqp->strong_pena_bc_coeff, eqname);
    break;

  case CS_EQKEY_BC_WEAK_PENA_COEFF:
    eqp->weak_pena_bc_coeff = atof(keyval);
    if (eqp->weak_pena_bc_coeff < 0.)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value of the penalization coefficient %5.3e\n"
                " This should be positive.\n"
                " Equation: %s\n",
                __func__, eqp->weak_pena_bc_coeff, eqname);
    break;

  case CS_EQKEY_DO_LUMPING:
    if (strcmp(keyval, "true") == 0 || strcmp(keyval, "1") == 0)
      eqp->do_lumping = true;
    else
      eqp->do_lumping = false;  /* Should be the default behavior */
    break;

  case CS_EQKEY_DOF_REDUCTION:
    if (strcmp(keyval, "derham") == 0 || strcmp(keyval, "de_rham") == 0)
      eqp->dof_reduction = CS_PARAM_REDUCTION_DERHAM;
    else if (strcmp(keyval, "average") == 0)
      eqp->dof_reduction = CS_PARAM_REDUCTION_AVERAGE;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_DOF_REDUCTION");
    }
    break;

  case CS_EQKEY_EXTRA_OP:
    if (strcmp(keyval, "balance") == 0)
      eqp->post_flag |= CS_EQUATION_POST_BALANCE;
    else if (strcmp(keyval, "peclet") == 0)
      eqp->post_flag |= CS_EQUATION_POST_PECLET;
    else if (strcmp(keyval, "upwind_coef") == 0)
      eqp->post_flag |= CS_EQUATION_POST_UPWIND_COEF;
    else if (strcmp(keyval, "normal_flux") == 0)
      eqp->post_flag |= CS_EQUATION_POST_NORMAL_FLUX;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_EXTRA_OP");
    }
    break;

  case CS_EQKEY_HODGE_DIFF_ALGO:
    if (strcmp(keyval,"cost") == 0 || strcmp(keyval,"ocs") == 0)
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_COST;
    else if (strcmp(keyval, "bubble") == 0) {
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_BUBBLE;
      eqp->diffusion_hodgep.coef = 2./3.;
    }
    else if (strcmp(keyval, "voronoi") == 0)
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_VORONOI;
    else if (strcmp(keyval, "wbs") == 0)
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_WBS;
    else if (strcmp(keyval, "auto") == 0)
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_AUTO;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_HODGE_DIFF_ALGO");
    }
    break;

  case CS_EQKEY_HODGE_DIFF_COEF:
    if (strcmp(keyval, "dga") == 0) {
      eqp->diffusion_hodgep.coef = 1./3.;
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_COST;
    }
    else if (strcmp(keyval, "sushi") == 0) {
      eqp->diffusion_hodgep.coef = 1./sqrt(3.);
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_COST;
    }
    else if (strcmp(keyval, "gcr") == 0) {
      eqp->diffusion_hodgep.coef = 1.0;
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_COST;
    }
    else if (strcmp(keyval, "frac23") == 0 || strcmp(keyval, "2/3") == 0) {
      eqp->diffusion_hodgep.coef = 2./3.;
    }
    else {
      eqp->diffusion_hodgep.coef = atof(keyval);
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_COST;
    }
    break;

  case CS_EQKEY_HODGE_TIME_ALGO:
    if (strcmp(keyval, "voronoi") == 0)
      eqp->time_hodgep.algo = CS_HODGE_ALGO_VORONOI;
    else if (strcmp(keyval,"bubble") == 0) {
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_BUBBLE;
      eqp->diffusion_hodgep.coef = 2./3.;
    }
    else if (strcmp(keyval,"cost") == 0 || strcmp(keyval,"ocs") == 0)
      eqp->time_hodgep.algo = CS_HODGE_ALGO_COST;
    else if (strcmp(keyval, "wbs") == 0)
      eqp->time_hodgep.algo = CS_HODGE_ALGO_WBS;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_HODGE_TIME_ALGO");
    }
    break;

  case CS_EQKEY_HODGE_REAC_ALGO:
    if (strcmp(keyval, "voronoi") == 0)
      eqp->reaction_hodgep.algo = CS_HODGE_ALGO_VORONOI;
    else if (strcmp(keyval,"bubble") == 0) {
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_BUBBLE;
      eqp->diffusion_hodgep.coef = 2./3.;
    }
    else if (strcmp(keyval,"cost") == 0 || strcmp(keyval,"ocs") == 0)
      eqp->time_hodgep.algo = CS_HODGE_ALGO_COST;
    else if (strcmp(keyval, "wbs") == 0)
      eqp->reaction_hodgep.algo = CS_HODGE_ALGO_WBS;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_HODGE_REAC_ALGO");
    }
    break;

  case CS_EQKEY_ITSOL:  /* kept for backward compatibility */
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS, wkmsg, __func__, eqname,
                  "CS_EQKEY_ITSOL", "CS_EQKEY_SOLVER");
    cs_log_printf_flush(CS_LOG_WARNINGS);
    [[fallthrough]]; /* No break wanted */
  case CS_EQKEY_SOLVER:
    ierr = cs_param_sles_set_solver(keyval, eqp->sles_param);

    if (ierr > 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_SOLVER");
    }
    break;

  case CS_EQKEY_ITSOL_ATOL:  /* kept for backward compatibility */
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS, wkmsg, __func__, eqname,
                  "CS_EQKEY_ITSOL_ATOL", "CS_EQKEY_SOLVER_ATOL");
    cs_log_printf_flush(CS_LOG_WARNINGS);
    [[fallthrough]]; /* No break wanted */
  case CS_EQKEY_SOLVER_ATOL:
    eqp->sles_param->cvg_param.atol = atof(keyval);
    break;

  case CS_EQKEY_ITSOL_DTOL:  /* kept for backward compatibility */
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS, wkmsg, __func__, eqname,
                  "CS_EQKEY_ITSOL_DTOL", "CS_EQKEY_SOLVER_DTOL");
    cs_log_printf_flush(CS_LOG_WARNINGS);
    [[fallthrough]]; /* No break wanted */
  case CS_EQKEY_SOLVER_DTOL:
    eqp->sles_param->cvg_param.dtol = atof(keyval);
    break;

  case CS_EQKEY_ITSOL_MAX_ITER:  /* kept for backward compatibility */
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS, wkmsg, __func__, eqname,
                  "CS_EQKEY_ITSOL_MAX_ITER", "CS_EQKEY_SOLVER_MAX_ITER");
    cs_log_printf_flush(CS_LOG_WARNINGS);
    [[fallthrough]]; /* No break wanted */
  case CS_EQKEY_SOLVER_MAX_ITER:
    eqp->sles_param->cvg_param.n_max_iter = atoi(keyval);
    break;

  case CS_EQKEY_SOLVER_NO_OP:
    if (strcmp(keyval, "true") == 0)
      eqp->sles_param->allow_no_op = true;
    else
      eqp->sles_param->allow_no_op = false;
    break;

  case CS_EQKEY_ITSOL_EPS:  /* kept for backward compatibility */
  case CS_EQKEY_ITSOL_RTOL:
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS,
                  " %s: Key \"%s\" or \"%s\" are deprecated.\n"
                  "  Please replace it with \"%s\"\n",
                  __func__, "CS_EQKEY_ITSOL_EPS", "CS_EQKEY_ITSOL_RTOL",
                  "CS_EQKEY_SOLVER_RTOL");
    cs_log_printf_flush(CS_LOG_WARNINGS);
    [[fallthrough]]; /* No break wanted */
  case CS_EQKEY_SOLVER_RTOL:
    eqp->sles_param->cvg_param.rtol = atof(keyval);
    if (eqp->sles_param->cvg_param.atol > eqp->sles_param->cvg_param.rtol)
      eqp->sles_param->cvg_param.atol = eqp->sles_param->cvg_param.rtol;
    break;

  case CS_EQKEY_ITSOL_RESNORM_TYPE: /* kept for backward compatibility */
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS, wkmsg, __func__, eqname,
                  "CS_EQKEY_ITSOL_RESNORM_TYPE",
                  "CS_EQKEY_SOLVER_RESNORM_TYPE");
    cs_log_printf_flush(CS_LOG_WARNINGS);
    [[fallthrough]]; /* No break wanted */
  case CS_EQKEY_SOLVER_RESNORM_TYPE:
    if (strcmp(keyval, "none") == 0 || strcmp(keyval, "false") == 0 ||
        strcmp(keyval, "") == 0)
      eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_NONE;
    else if (strcmp(keyval, "rhs") == 0)
      eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_NORM2_RHS;
    else if (strcmp(keyval, "weighted_rhs") == 0 ||
             strcmp(keyval, "weighted") == 0)
      eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_WEIGHTED_RHS;
    else if (strcmp(keyval, "filtered_rhs") == 0 ||
             strcmp(keyval, "filtered") == 0)
      eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_FILTERED_RHS;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_SOLVER_RESNORM_TYPE");
    }
    break;

  case CS_EQKEY_ITSOL_RESTART:  /* kept for backward compatibility */
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_WARNINGS, wkmsg, __func__, eqname,
                  "CS_EQKEY_ITSOL_RESTART", "CS_EQKEY_SOLVER_RESTART");
    cs_log_printf_flush(CS_LOG_WARNINGS);
    [[fallthrough]]; /* No break wanted */
  case CS_EQKEY_SOLVER_RESTART:
    eqp->sles_param->restart = atoi(keyval);
    break;

  case CS_EQKEY_PRECOND:
    ierr = cs_param_sles_set_precond(keyval, eqp->sles_param);

    if (ierr > 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_PRECOND");
    }
    break;

  case CS_EQKEY_PRECOND_BLOCK_TYPE:
    ierr = cs_param_sles_set_precond_block_type(keyval, eqp->sles_param);

    if (ierr > 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_PRECOND_BLOCK_TYPE");
    }
    break;

  case CS_EQKEY_SADDLE_ATOL:
    eqp->saddle_param->cvg_param.atol = atof(keyval);
    if (eqp->saddle_param->cvg_param.atol < 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_SADDLE_ATOL");
    }
    break;

  case CS_EQKEY_SADDLE_AUGMENT_SCALING:
    {
      double  gamma = atof(keyval);

      if (gamma < 0) {
        const char *_val = keyval;
        bft_error(__FILE__, __LINE__, 0, emsg,
                  __func__, eqname, _val, "CS_EQKEY_SADDLE_AUGMENT_SCALING");
      }

      cs_param_saddle_set_augmentation_coef(eqp->saddle_param, gamma);
    }
    break;

  case CS_EQKEY_SADDLE_DTOL:
    eqp->saddle_param->cvg_param.dtol = atof(keyval);
    if (eqp->saddle_param->cvg_param.dtol < 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_SADDLE_DTOL");
    }
    break;

  case CS_EQKEY_SADDLE_MAX_ITER:
    eqp->saddle_param->cvg_param.n_max_iter = atoi(keyval);
    break;

  case CS_EQKEY_SADDLE_PRECOND:
    ierr = cs_param_saddle_set_precond(keyval, eqp->saddle_param);
    if (ierr != 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_SADDLE_PRECOND");
    }
    break;

  case CS_EQKEY_SADDLE_RTOL:
    eqp->saddle_param->cvg_param.rtol = atof(keyval);
    if (eqp->saddle_param->cvg_param.rtol < 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_SADDLE_RTOL");
    }
    break;

  case CS_EQKEY_SADDLE_SCHUR_APPROX:
    ierr = cs_param_saddle_set_schur_approx(keyval, eqp->saddle_param);
    if (ierr != 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_SADDLE_SCHUR_APPROX");
    }
    break;

  case CS_EQKEY_SADDLE_SOLVER:
    ierr = cs_param_saddle_set_solver(keyval, eqp->saddle_param);

    switch (ierr) {
    case 1:
      {
        const char *_val = keyval;
        bft_error(__FILE__, __LINE__, 0,
                  emsg, __func__, eqname, _val, "CS_EQKEY_SADDLE_SOLVER");
      }
      break;
    case 2:
      bft_error(__FILE__, __LINE__, 0,
                " %s(): Eq. %s Error detected while setting \"%s\" key.\n"
                " PETSc is not available.\n"
                " Please check your installation settings.\n",
                __func__, eqname, "CS_EQKEY_SADDLE_SOLVER");

      break;
    case 3:
      bft_error(__FILE__, __LINE__, 0,
                " %s(): Eq. %s Error detected while setting \"%s\" key.\n"
                " MUMPS is not available.\n"
                " Please check your installation settings.\n",
                __func__, eqname, "CS_EQKEY_SADDLE_SOLVER");
      break;
    default:
      break;
    }
    break;

  case CS_EQKEY_SADDLE_SOLVER_CLASS:
    ierr = cs_param_saddle_set_solver_class(keyval, eqp->saddle_param);

    switch (ierr) {
    case 1:
      {
        const char *_val = keyval;
        bft_error(__FILE__, __LINE__, 0,
                  emsg, __func__, eqname, _val, "CS_EQKEY_SADDLE_SOLVER_CLASS");
      }
      break;
    case 2:
      bft_error(__FILE__, __LINE__, 0,
                " %s(): Eq. %s Error detected while setting \"%s\" key.\n"
                " PETSc is not available.\n"
                " Please check your installation settings.\n",
                __func__, eqname, "CS_EQKEY_SADDLE_SOLVER_CLASS");

      break;
    case 3:
      bft_error(__FILE__, __LINE__, 0,
                " %s(): Eq. %s Error detected while setting \"%s\" key.\n"
                " MUMPS is not available.\n"
                " Please check your installation settings.\n",
                __func__, eqname, "CS_EQKEY_SADDLE_SOLVER_CLASS");
      break;
    default:
      break;
    }
    break;

  case CS_EQKEY_SADDLE_SOLVER_RESTART:
    {
      int  restart_range = atoi(keyval);

      cs_param_saddle_set_restart_range(eqp->saddle_param, restart_range);

      if (restart_range < 2) {
        const char *_val = keyval;
        bft_error(__FILE__, __LINE__, 0,
                  emsg, __func__, eqname, _val,
                  "CS_EQKEY_SADDLE_SOLVER_RESTART");
      }
    }
    break;

  case CS_EQKEY_SADDLE_VERBOSITY:
    eqp->saddle_param->verbosity = atoi(keyval);
    break;

  case CS_EQKEY_SLES_VERBOSITY: /* "verbosity" for SLES structures */
    eqp->sles_param->verbosity = atoi(keyval);
    break;

  case CS_EQKEY_SOLVER_FAMILY:
    ierr = cs_param_sles_set_solver_class(keyval, eqp->sles_param);

    if (ierr > 0) {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_SOLVER_FAMILY");
    }
    break;

  case CS_EQKEY_SPACE_SCHEME:
    if (strcmp(keyval, "cdo_vb") == 0 ||
        strcmp(keyval, "cdovb") == 0) {

      eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;
      eqp->space_poly_degree = 0;

      /* Set the corresponding default settings */

      eqp->time_hodgep.type = CS_HODGE_TYPE_VPCD;
      eqp->time_hodgep.algo = CS_HODGE_ALGO_VORONOI;

      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_EPFD;
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_BUBBLE;
      eqp->diffusion_hodgep.coef = 2./3.;

      eqp->reaction_hodgep.type = CS_HODGE_TYPE_VPCD;
      eqp->reaction_hodgep.algo = CS_HODGE_ALGO_VORONOI;

    }
    else if (strcmp(keyval, "cdo_vcb") == 0 ||
             strcmp(keyval, "cdovcb") == 0) {

      eqp->space_scheme = CS_SPACE_SCHEME_CDOVCB;
      eqp->space_poly_degree = 0;

      eqp->time_hodgep.type = CS_HODGE_TYPE_VPCD;
      eqp->time_hodgep.algo = CS_HODGE_ALGO_WBS;

      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_WBS;
      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_VC;

      eqp->reaction_hodgep.type = CS_HODGE_TYPE_VPCD;
      eqp->reaction_hodgep.algo = CS_HODGE_ALGO_WBS;

    }
    else if (strcmp(keyval, "cdo_fb") == 0 ||
             strcmp(keyval, "cdofb") == 0) {

      eqp->space_scheme = CS_SPACE_SCHEME_CDOFB;
      eqp->space_poly_degree = 0;

      eqp->time_hodgep.type = CS_HODGE_TYPE_FB;
      eqp->time_hodgep.algo = CS_HODGE_ALGO_VORONOI;

      eqp->reaction_hodgep.type = CS_HODGE_TYPE_FB;
      eqp->reaction_hodgep.algo = CS_HODGE_ALGO_VORONOI;

      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_EDFP;
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_COST;
      eqp->diffusion_hodgep.coef = 2./3.;

    }
    else if (strcmp(keyval, "cdo_cb") == 0 ||
             strcmp(keyval, "cdocb") == 0)
      cs_cdocb_init_default_param(eqp);


    else if (strcmp(keyval, "cdo_eb") == 0 ||
             strcmp(keyval, "cdoeb") == 0) {

      eqp->space_scheme = CS_SPACE_SCHEME_CDOEB;
      eqp->space_poly_degree = 0;

      eqp->time_hodgep.type = CS_HODGE_TYPE_EPFD;

      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_FPED;

      eqp->reaction_hodgep.type = CS_HODGE_TYPE_EPFD;
    }

    /* Only diffusion is implemented for HHO schemes up to now */

    else if (strcmp(keyval, "hho_p0") == 0) {

      eqp->space_scheme = CS_SPACE_SCHEME_HHO_P0;
      eqp->space_poly_degree = 0;

      eqp->time_hodgep.type = CS_HODGE_TYPE_CPVD;

      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_EDFP;

    }
    else if (strcmp(keyval, "hho_p1") == 0) {

      eqp->space_scheme = CS_SPACE_SCHEME_HHO_P1;
      eqp->space_poly_degree = 1;

      eqp->time_hodgep.type = CS_HODGE_TYPE_CPVD;

      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_EDFP;

    }
    else if (strcmp(keyval, "hho_p2") == 0) {

      eqp->space_scheme = CS_SPACE_SCHEME_HHO_P2;
      eqp->space_poly_degree = 2;

      eqp->time_hodgep.type = CS_HODGE_TYPE_CPVD;

      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_EDFP;
    }
    else if (strcmp(keyval, "mac_fb") == 0) {

      eqp->space_scheme      = CS_SPACE_SCHEME_MACFB;
      eqp->space_poly_degree = 0;

      eqp->time_hodgep.type = CS_HODGE_TYPE_CPVD;

      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_EDFP;
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_VORONOI;
      eqp->diffusion_hodgep.coef = 0.0;
    }
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_SPACE_SCHEME");
    }
    break;

  case CS_EQKEY_TIME_SCHEME:

    if (strcmp(keyval, "no") == 0 || strcmp(keyval, "steady") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_STEADY;
    }
    else if (strcmp(keyval, "euler_implicit") == 0 ||
             strcmp(keyval, "forward_euler") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_EULER_IMPLICIT;
      eqp->theta = 1.;
    }
    else if (strcmp(keyval, "euler_explicit") == 0 ||
             strcmp(keyval, "backward_euler") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_EULER_EXPLICIT;
      eqp->theta = 0.;
    }
    else if (strcmp(keyval, "crank_nicolson") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_CRANKNICO;
      eqp->theta = 0.5;
    }
    else if (strcmp(keyval, "theta_scheme") == 0)
      eqp->time_scheme = CS_TIME_SCHEME_THETA;
    else if (strcmp(keyval, "bdf2") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_BDF2;
      bft_error(__FILE__, __LINE__, 0, " Eq. %s\n"
                " Soon available...", eqname);
    }
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_TIME_SCHEME");
    }
    break;

  case CS_EQKEY_TIME_THETA:
    eqp->theta = atof(keyval);
    break;

  case CS_EQKEY_VERBOSITY: /* "verbosity" */
    eqp->verbosity = atoi(keyval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid key for setting the equation %s."),
              __func__, eqname);

  } /* Switch on keys */
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a \ref cs_equation_param_t structure with the given
 *        parameters. The remaining parameters are set with default values;
 *
 * \param[in] name        name of the equation
 * \param[in] type        type of equation
 * \param[in] dim         dim of the variable associated to this equation
 * \param[in] default_bc  type of boundary condition set by default
 *
 * \return a pointer to a new allocated \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_create(const char         *name,
                         cs_equation_type_t  type,
                         int                 dim,
                         cs_param_bc_type_t  default_bc)
{
  cs_equation_param_t *eqp = nullptr;

  BFT_MALLOC(eqp, 1, cs_equation_param_t);

  /* Store the name of the equation */

  size_t  len = strlen(name);
  BFT_MALLOC(eqp->name, len + 1, char);
  strncpy(eqp->name, name, len + 1);

  /* Set additional members */

  eqp->type = type;
  eqp->dim = dim;

  /* Other default settings */

  eqp->verbosity = 1;
  eqp->flag = 0;
  eqp->post_flag = 0;

  /* Vertex-based schemes imply specific discrete Hodge operators for
     diffusion, time and reaction terms.
     Default initialization is made in accordance with this choice */

  eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;
  eqp->dof_reduction = CS_PARAM_REDUCTION_DERHAM;
  eqp->space_poly_degree = 0;

  /* Default initialization for the legacy var_col_opt structure which is now
   * shared inside the cs_equation_param_t structure The default value used
   * here should be the same as the one considered in src/base/cs_parameters.c
   * (static variable called _equation_param_default) */

  eqp->iconv  = 1;
  eqp->istat  = 1;
  eqp->idircl = 1;
  eqp->ndircl = 0;
  eqp->idiff  = 1;
  eqp->idifft = 1;
  eqp->idften = CS_ISOTROPIC_DIFFUSION;
  eqp->iswdyn = -1;
  eqp->ischcv = 1;
  eqp->ibdtso = 1;
  eqp->isstpc = 1;
  eqp->nswrgr = 100;
  eqp->nswrsm = 1;
  eqp->imvisf = 0;
  eqp->imrgra = -1;
  eqp->imligr = -1;
  eqp->ircflu = 1;
  eqp->iwgrec = 0;
  eqp->icoupl = -1;
  eqp->theta = 1.;
  eqp->blencv = 1.;
  eqp->blend_st = 0.;
  eqp->epsilo = 1.e-5;
  eqp->epsrsm = 1.e-4;
  eqp->epsrgr = 1.e-4;
  eqp->climgr = 1.5;
  eqp->relaxv = 1.;
  eqp->b_gradient_r = 2;
  eqp->b_diff_flux_rc = 1;

  /* Boundary conditions structure.
     One assigns a boundary condition by default */

  eqp->default_bc = default_bc;
  eqp->n_bc_defs = 0;
  eqp->bc_defs              = nullptr;
  eqp->default_enforcement = CS_PARAM_BC_ENFORCE_ALGEBRAIC;
  eqp->strong_pena_bc_coeff = _strong_pena_bc_coef_by_default;
  eqp->weak_pena_bc_coeff = _weak_pena_bc_coef_by_default;

  /* Initial condition (zero value by default) */

  eqp->n_ic_defs = 0;
  eqp->ic_defs   = nullptr;

  /* Description of the time discretization (default values) */

  eqp->time_property = nullptr;
  eqp->time_scheme = CS_TIME_SCHEME_EULER_IMPLICIT;
  eqp->do_lumping = false;
  eqp->time_hodgep = (cs_hodge_param_t) {
    .inv_pty = false,
    .type = CS_HODGE_TYPE_VPCD,
    .algo = CS_HODGE_ALGO_VORONOI,
    .coef = 1., // not used by default with Voronoi
  };

  /* Description of the discretization of the diffusion term */

  eqp->diffusion_property = nullptr;
  eqp->diffusion_hodgep = (cs_hodge_param_t) {
    .inv_pty = false,
    .type = CS_HODGE_TYPE_EPFD,
    .algo = CS_HODGE_ALGO_BUBBLE,
    .coef = 2/3.,
  };

  /* Description of the discretization of the curl-curl term */

  eqp->curlcurl_property = nullptr;
  eqp->curlcurl_hodgep = (cs_hodge_param_t) {
    .inv_pty = true,
    .type = CS_HODGE_TYPE_FPED,
    .algo = CS_HODGE_ALGO_COST,
    .coef = 1./3.,
  };

  /* Description of the discretization of the grad-div term */

  eqp->graddiv_property = nullptr;
  eqp->graddiv_hodgep = (cs_hodge_param_t) {
    .inv_pty = false,
    .type = CS_HODGE_TYPE_EPFD,
    .algo = CS_HODGE_ALGO_VORONOI,
    .coef = 1., // not used by default with Voronoi
  };

  /* Advection term */

  eqp->adv_field            = nullptr;
  eqp->adv_scaling_property = nullptr;
  eqp->adv_extrapol = CS_PARAM_ADVECTION_EXTRAPOL_NONE;
  eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_CONSERV;
  eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
  eqp->adv_strategy = CS_PARAM_ADVECTION_IMPLICIT_FULL;
  eqp->upwind_portion = 0.15;
  eqp->cip_scaling_coef = -1.0;  // Automatic settings


  /* Description of the discretization of the reaction term.
     No reaction term by default */

  eqp->n_reaction_terms = 0;
  eqp->reaction_properties = nullptr;
  eqp->reaction_hodgep = (cs_hodge_param_t) {
    .inv_pty = false,
    .type = CS_HODGE_TYPE_VPCD,
    .algo = CS_HODGE_ALGO_VORONOI,
    .coef = 1., // not used by default with Voronoi
  };

  /* Source term (always in the right-hand side)
     No source term by default */

  eqp->n_source_terms = 0;
  eqp->source_terms   = nullptr;

  /* Mass injection in the volume term (always in the right-hand side)
     No volume mass injection term by default */

  eqp->n_volume_mass_injections = 0;
  eqp->volume_mass_injections   = nullptr;

  /* Members of the structure handling the enforcement of (internal) DoFs */

  eqp->n_enforcements = 0;
  eqp->enforcement_params = nullptr;

  /* Settings for driving the linear algebra */

  eqp->sles_param = cs_param_sles_create(-1, name); /* field_id, system_name */

  eqp->saddle_param = cs_param_saddle_create(); /* No saddle-point solver by
                                                   default */

  /* By default, there is no incremental solving */

  eqp->incremental_algo_type = CS_PARAM_NL_ALGO_NONE;
  eqp->incremental_algo_cvg = (cs_param_convergence_t) {
    .atol = 1e-6,           /* absolute tolerance */
    .rtol = 1e-2,           /* relative tolerance */
    .dtol = 1e3,            /* divergence tolerance */
    .n_max_iter = 15 };     /* n_max iter. */

  eqp->incremental_relax_factor = 1.0; /* No relaxation by default */

  eqp->incremental_anderson_param = (cs_iter_algo_param_aac_t) {
    .n_max_dir = 6,
    .starting_iter = 3,
    .max_cond = -1, /* No test by default */
    .beta = 1.0,    /* No damping by default */
    .dp_type = CS_PARAM_DOTPROD_EUCLIDEAN };

  eqp->time_control_owner = false;
  eqp->time_control = nullptr; /* nullptr by default (allways active) */

  return eqp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy the settings from one \ref cs_equation_param_t structure to
 *        another one. The name is not copied.
 *
 * \param[in]      ref          pointer to the reference cs_equation_param_t
 * \param[in, out] dst          pointer to the cs_equation_param_t to update
 * \param[in]      copy_fld_id  copy also the field id or not
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_copy_from(const cs_equation_param_t *ref,
                            cs_equation_param_t       *dst,
                            bool                       copy_fld_id)
{
  /* Generic members */

  dst->type = ref->type;
  dst->dim = ref->dim;
  dst->verbosity = ref->verbosity;
  dst->post_flag = ref->post_flag;
  dst->flag = ref->flag;
  dst->space_scheme = ref->space_scheme;
  dst->dof_reduction = ref->dof_reduction;
  dst->space_poly_degree = ref->space_poly_degree;

  /* Members originally located in the cs_var_cal_opt_t structure */

  dst->iconv  = ref->iconv;
  dst->istat  = ref->istat;
  dst->idircl = ref->idircl;
  dst->ndircl = ref->ndircl;
  dst->idiff  = ref->idiff;
  dst->idifft = ref->idifft;
  dst->idften = ref->idften;
  dst->iswdyn = ref->iswdyn;
  dst->ischcv = ref->ischcv;
  dst->ibdtso = ref->ibdtso;
  dst->isstpc = ref->isstpc;
  dst->nswrgr = ref->nswrgr;
  dst->nswrsm = ref->nswrsm;
  dst->imvisf = ref->imvisf;
  dst->imrgra = ref->imrgra;
  dst->imligr = ref->imligr;
  dst->ircflu = ref->ircflu;
  dst->iwgrec = ref->iwgrec;
  dst->icoupl = ref->icoupl;
  dst->theta = ref->theta;
  dst->blencv = ref->blencv;
  dst->blend_st = ref->blend_st;
  dst->epsilo = ref->epsilo;
  dst->epsrsm = ref->epsrsm;
  dst->epsrgr = ref->epsrgr;
  dst->climgr = ref->climgr;
  dst->relaxv = ref->relaxv;

  /* Boundary conditions structure */

  dst->default_bc = ref->default_bc;
  dst->default_enforcement = ref->default_enforcement;
  dst->strong_pena_bc_coeff = ref->strong_pena_bc_coeff;
  dst->n_bc_defs = ref->n_bc_defs;
  BFT_REALLOC(dst->bc_defs, dst->n_bc_defs, cs_xdef_t *);
  for (int i = 0; i < ref->n_bc_defs; i++)
    dst->bc_defs[i] = cs_xdef_copy(ref->bc_defs[i]);

  /* Description of the time discretization */

  dst->time_scheme = ref->time_scheme;
  dst->do_lumping = ref->do_lumping;
  dst->time_property = ref->time_property;

  cs_hodge_copy_parameters(&(ref->time_hodgep), &(dst->time_hodgep));

  /* Initial condition (zero value by default) */

  dst->n_ic_defs = ref->n_ic_defs;
  BFT_REALLOC(dst->ic_defs, dst->n_ic_defs, cs_xdef_t *);
  for (int i = 0; i < ref->n_ic_defs; i++)
    dst->ic_defs[i] = cs_xdef_copy(ref->ic_defs[i]);

  /* Diffusion term */

  dst->diffusion_property = ref->diffusion_property;

  cs_hodge_copy_parameters(&(ref->diffusion_hodgep), &(dst->diffusion_hodgep));

  /* Curl-curl term */

  dst->curlcurl_property = ref->curlcurl_property;

  cs_hodge_copy_parameters(&(ref->curlcurl_hodgep), &(dst->curlcurl_hodgep));

  /* Grad-div term */

  dst->graddiv_property = ref->graddiv_property;

  cs_hodge_copy_parameters(&(ref->graddiv_hodgep), &(dst->graddiv_hodgep));

  /* Advection term */

  dst->adv_formulation = ref->adv_formulation;
  dst->adv_scheme = ref->adv_scheme;
  dst->adv_strategy = ref->adv_strategy;
  dst->adv_extrapol = ref->adv_extrapol;
  dst->upwind_portion = ref->upwind_portion;
  dst->adv_field = ref->adv_field;
  dst->adv_scaling_property = ref->adv_scaling_property;

  /* Reaction term */

  dst->n_reaction_terms = ref->n_reaction_terms;
  BFT_REALLOC(dst->reaction_properties, dst->n_reaction_terms, cs_property_t *);
  for (int i = 0; i < ref->n_reaction_terms; i++)
    dst->reaction_properties[i] = ref->reaction_properties[i];

  cs_hodge_copy_parameters(&(ref->reaction_hodgep), &(dst->reaction_hodgep));

  /* Source term */

  dst->n_source_terms = ref->n_source_terms;
  BFT_MALLOC(dst->source_terms, dst->n_source_terms, cs_xdef_t *);
  for (int i = 0; i < dst->n_source_terms; i++)
    dst->source_terms[i] = cs_xdef_copy(ref->source_terms[i]);

  /* Mass injection term */

  dst->n_volume_mass_injections = ref->n_volume_mass_injections;
  BFT_REALLOC(dst->volume_mass_injections,
              dst->n_volume_mass_injections,
              cs_xdef_t *);
  for (int i = 0; i < dst->n_volume_mass_injections; i++)
    dst->volume_mass_injections[i]
      = cs_xdef_copy(ref->volume_mass_injections[i]);

  /* Enforcement of internal DoFs */

  if (ref->n_enforcements > 0) {

    dst->n_enforcements = ref->n_enforcements;
    BFT_REALLOC(dst->enforcement_params, dst->n_enforcements,
                cs_enforcement_param_t *);

    for (int i = 0; i < ref->n_enforcements; i++)
      dst->enforcement_params[i] =
        cs_enforcement_param_copy(ref->enforcement_params[i]);

  }
  else {

    if (dst->n_enforcements > 0)
      BFT_FREE(dst->enforcement_params);

    dst->enforcement_params = nullptr;
    dst->n_enforcements = 0;

  }

  /* Copy the settings driving the linear algebra algorithms */

  if (copy_fld_id)
    cs_param_sles_copy_from(ref->sles_param, dst->sles_param);

  else {

    int save_field_id = dst->sles_param->field_id;
    cs_param_sles_copy_from(ref->sles_param, dst->sles_param);
    dst->sles_param->field_id = save_field_id;

  }

  cs_param_saddle_copy(ref->saddle_param, dst->saddle_param);

  /* Non-linear algorithm */

  dst->incremental_algo_type = ref->incremental_algo_type;
  dst->incremental_relax_factor = ref->incremental_relax_factor;

  {
    const cs_param_convergence_t  cvgp = ref->incremental_algo_cvg;

    dst->incremental_algo_cvg.rtol = cvgp.rtol;
    dst->incremental_algo_cvg.atol = cvgp.atol;
    dst->incremental_algo_cvg.dtol = cvgp.dtol;
    dst->incremental_algo_cvg.n_max_iter = cvgp.n_max_iter;
  }

  {
    const cs_iter_algo_param_aac_t  aacp = ref->incremental_anderson_param;

    dst->incremental_anderson_param.n_max_dir = aacp.n_max_dir;
    dst->incremental_anderson_param.starting_iter = aacp.starting_iter;
    dst->incremental_anderson_param.max_cond = aacp.max_cond;
    dst->incremental_anderson_param.beta = aacp.beta;
    dst->incremental_anderson_param.dp_type = aacp.dp_type;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Copy only the part dedicated to the boundary conditions and the DoF
 *        (degrees of freedom) enforcement from one \ref cs_equation_param_t
 *        structure to another one.
 *
 * \param[in]      ref  pointer to the reference \ref cs_equation_param_t
 * \param[in, out] dst  pointer to the \ref cs_equation_param_t to update
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_copy_bc(const cs_equation_param_t *ref,
                          cs_equation_param_t       *dst)
{
  if (ref == nullptr)
    return;

  if (dst == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Structure is not allocated.\n"
              "%s: Stop copying a cs_equation_param_t structure.\n",
              __func__, __func__);

  /* Boundary conditions structure */

  dst->default_bc = ref->default_bc;
  dst->default_enforcement = ref->default_enforcement;
  dst->strong_pena_bc_coeff = ref->strong_pena_bc_coeff;
  dst->n_bc_defs = ref->n_bc_defs;
  BFT_REALLOC(dst->bc_defs, dst->n_bc_defs, cs_xdef_t *);
  for (int i = 0; i < ref->n_bc_defs; i++)
    dst->bc_defs[i] = cs_xdef_copy(ref->bc_defs[i]);

  /* Enforcement of internal DoFs */

  if (ref->n_enforcements > 0) {

    dst->n_enforcements = ref->n_enforcements;
    BFT_REALLOC(dst->enforcement_params, dst->n_enforcements,
                cs_enforcement_param_t *);

    for (int i = 0; i < ref->n_enforcements; i++)
      dst->enforcement_params[i] =
        cs_enforcement_param_copy(ref->enforcement_params[i]);

  }
  else {

    if (dst->n_enforcements > 0)
      BFT_FREE(dst->enforcement_params);

    dst->enforcement_params = nullptr;
    dst->n_enforcements = 0;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the contents of a \ref cs_equation_param_t
 *
 * The cs_equation_param_t structure itself is not freed, but all its
 * sub-structures are freed.
 *
 * This is useful for equation parameters which are accessed through
 * field keywords.
 *
 * \param[in, out] eqp  pointer to a \ref cs_equation_param_t
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_clear(cs_equation_param_t *eqp)
{
  if (eqp == nullptr)
    return;

  /* Information related to the definition of the boundary conditions */

  if (eqp->n_bc_defs > 0) {

    for (int i = 0; i < eqp->n_bc_defs; i++)
      eqp->bc_defs[i] = cs_xdef_free(eqp->bc_defs[i]);
    BFT_FREE(eqp->bc_defs);

  }

  /* Information related to the definition of reaction terms */

  if (eqp->n_reaction_terms > 0) {

    BFT_FREE(eqp->reaction_properties);

    /* Remark: properties are freed when the global cs_domain_t structure is
       freed thanks to a call to cs_property_free() */
  }

  /* Information related to the definition of source terms */

  if (eqp->n_source_terms > 0) {

    for (int i = 0; i < eqp->n_source_terms; i++)
      eqp->source_terms[i] = cs_xdef_free(eqp->source_terms[i]);
    BFT_FREE(eqp->source_terms);

  }

  /* Information related to the definition of mass injection terms */

  if (eqp->n_volume_mass_injections > 0) {

    for (int i = 0; i < eqp->n_volume_mass_injections; i++)
      eqp->volume_mass_injections[i]
        = cs_xdef_free(eqp->volume_mass_injections[i]);
    BFT_FREE(eqp->volume_mass_injections);

  }

  /* Information related to the enforcement of internal DoFs */

  if (eqp->n_enforcements > 0) {

    for (int i = 0; i < eqp->n_enforcements; i++)
      cs_enforcement_param_free(&(eqp->enforcement_params[i]));
    BFT_FREE(eqp->enforcement_params);

  }

  /* Information related to the definition of initial conditions */

  if (eqp->n_ic_defs > 0) {

    for (int i = 0; i < eqp->n_ic_defs; i++)
      eqp->ic_defs[i] = cs_xdef_free(eqp->ic_defs[i]);
    BFT_FREE(eqp->ic_defs);

  }

  /* Information related to the linear algebra settings */

  cs_param_sles_free(&(eqp->sles_param));
  cs_param_saddle_free(&(eqp->saddle_param));

  /* Information related to the time_control structure */
  if (eqp->time_control_owner)
    BFT_FREE(eqp->time_control);

  BFT_FREE(eqp->name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free a \ref cs_equation_param_t
 *
 * \param[in, out] eqp  pointer to a \ref cs_equation_param_t
 *
 * \return a null pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_free(cs_equation_param_t *eqp)
{
  if (eqp == nullptr)
    return nullptr;

  cs_equation_param_clear(eqp);

  BFT_FREE(eqp);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a parameter attached to a keyname in a \ref cs_equation_param_t
 *        structure
 *
 * \param[in, out] eqp     pointer to a \ref cs_equation_param_t structure
 * \param[in]      key     key related to the member of eq to set
 * \param[in]      keyval  accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set(cs_equation_param_t *eqp,
                      cs_equation_key_t    key,
                      const char          *keyval)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (keyval == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: Eq: %s: Key value is empty",
              __func__, eqp->name);
  if (eqp->flag & CS_EQUATION_LOCKED)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Equation %s is not modifiable anymore.\n"
                " Please check your settings."), eqp->name, __func__);

  /* Conversion of the string to lower case */

  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  /* Set the couple (key,keyval) */

  _set_key(eqp, key, val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the pointer to the set of parameters to handle the SLES solver
 *        associated to this set of equation parameters
 *
 * \param[in] eqp  pointer to a \ref cs_equation_param_t structure
 *
 * \return a pointer to a cs_param_sles_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_sles_t *
cs_equation_param_get_sles_param(cs_equation_param_t *eqp)
{
  if (eqp == nullptr)
    return nullptr;

  return eqp->sles_param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the pointer to the set of parameters to handle a saddle-point
 *        solver This is only useful if some cases (the momemtum equation when
 *        a monolithic coupling is considered or the solution of a
 *        scalar-valued CDO cell-based schemes)
 *
 * \param[in] eqp  pointer to a \ref cs_equation_param_t structure
 *
 * \return a pointer to a cs_param_saddle_t structure
 */
/*----------------------------------------------------------------------------*/

cs_param_saddle_t *
cs_equation_param_get_saddle_param(cs_equation_param_t *eqp)
{
  if (eqp == nullptr)
    return nullptr;

  return eqp->saddle_param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set parameters for initializing SLES structures used for the
 *        resolution of the linear system.
 *        Settings are related to this equation.
 *
 * \param[in, out] eqp  pointer to a cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_sles(cs_equation_param_t *eqp)
{
  cs_param_saddle_t  *saddlep = eqp->saddle_param;

  /* Perform the setup only if the corresponding system is not a saddle-point
     system. For saddle-point system, the setup is performed in
     cs_saddle_solver_setup_sles() */

  if (saddlep->solver == CS_PARAM_SADDLE_SOLVER_NONE) {

    /* Define a cs_sles_t structure using the field_id related to the variable
     * field associated to this equation */

    int  ierr = cs_param_sles_setup(true, eqp->sles_param);

    if (ierr < 0)
      bft_error(__FILE__, __LINE__, 0,
                "%s: The requested solver setting is not available"
                " for the equation %s\n"
                " Error %d found. Please modify your settings.",
                __func__, eqp->name, ierr);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply the given quadrature rule to all existing definitions under the
 *        cs_equation_param_t structure.
 *
 *        If the default quadrature has been modified by the code, this
 *        function set the value of the quadrture to the given parameter
 *        whatever was the previous value.
 *
 *        To get a more detailed control of the quadrature rule, please
 *        consider the function \ref cs_xdef_set_quadrature
 *
 * \param[in, out] eqp    pointer to a \ref cs_equation_param_t structure
 * \param[in]      qtype  type of quadrature to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_quadrature_to_all(cs_equation_param_t  *eqp,
                                        cs_quadrature_type_t  qtype)
{
  if (eqp == nullptr)
    return;

  /* Apply the quadrature rule to all BC definitions */

  for (int i = 0; i < eqp->n_bc_defs; i++)
    cs_xdef_set_quadrature(eqp->bc_defs[i], qtype);

  /* Apply the quadrature rule to all source term definitions */

  for (int i = 0; i < eqp->n_source_terms; i++)
    cs_xdef_set_quadrature(eqp->source_terms[i], qtype);

  /* Apply the quadrature rule to all IC definitions */

  for (int i = 0; i < eqp->n_ic_defs; i++)
    cs_xdef_set_quadrature(eqp->ic_defs[i], qtype);

  /* Apply the quadrature rule to all volume_mass_injections */

  for (int i = 0; i < eqp->n_volume_mass_injections; i++)
    cs_xdef_set_quadrature(eqp->volume_mass_injections[i], qtype);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Lock settings to prevent from unwanted modifications.
 *
 * \param[in, out] eqp  pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_lock_settings(cs_equation_param_t *eqp)
{
  if (eqp == nullptr)
    return;

  eqp->flag |= CS_EQUATION_LOCKED;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Unlock settings. Be sure that is really wanted (inconsistency between
 *        the setup logging and what is used may appear)
 *
 * \param[in, out] eqp  pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_unlock_settings(cs_equation_param_t *eqp)
{
  if (eqp == nullptr)
    return;

  eqp->flag -= CS_EQUATION_LOCKED;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief At this stage, the numerical settings should not be modified anymore
 *        by the user. One makes a last set of modifications if needed to
 *        ensure a consistent numerical settings.
 *
 * \param[in, out] eqp  pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_ensure_consistent_settings(cs_equation_param_t *eqp)
{
  if (eqp == nullptr)
    return;

  if (eqp->flag & CS_EQUATION_LOCKED)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Equation %s is not modifiable anymore.\n"
                " Please check your settings."), eqp->name, __func__);

  /* Check the consistency of the settings related to the diffusion term */

  switch (eqp->space_scheme) {
  case CS_SPACE_SCHEME_CDOVCB:
    if (eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_COST   ||
        eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_BUBBLE ||
        eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_VORONOI) {

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_WARNINGS,
                    "%s: Incompatible Hodge algo. for the diffusion term.\n"
                    "%s: Equation \"%s\": Switch to a WBS algo.\n"
                    "%s: Please check your settings.\n",
                    __func__, __func__, eqp->name, __func__);
      cs_log_printf_flush(CS_LOG_WARNINGS);

      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_WBS;

    }
    break;

  case CS_SPACE_SCHEME_CDOFB:
    if (eqp->diffusion_hodgep.algo == CS_HODGE_ALGO_WBS) {

      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_WARNINGS,
                    "%s: Incompatible Hodge algo. for the diffusion term.\n"
                    "%s: Equation \"%s\": Switch to a COST algo.\n"
                    "%s: Please check your settings.\n",
                    __func__, __func__, eqp->name, __func__);
      cs_log_printf_flush(CS_LOG_WARNINGS);

      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_COST;

    }
    break;

  default:
    break;  /* Do nothing */
  }

  /* Activate a set of options in order to be consistent with the lumping */

  if (eqp->do_lumping) {

    eqp->reaction_hodgep.algo = CS_HODGE_ALGO_VORONOI;
    eqp->time_hodgep.algo = CS_HODGE_ALGO_VORONOI;

    /* A simple barycentric quadrature rule is applied to all source terms */

    for (int i = 0; i < eqp->n_source_terms; i++)
      cs_xdef_set_quadrature(eqp->source_terms[i], CS_QUADRATURE_BARY);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the detail of a \ref cs_equation_param_t structure
 *
 * \param[in] eqp  pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_log(const cs_equation_param_t *eqp)
{
  if (eqp == nullptr)
    return;

  const char *eqname = eqp->name;

  char  prefix[256];
  assert(strlen(eqname) < 200); /* Check that prefix is large enough */

  /* High-level settings */

  cs_log_printf(CS_LOG_SETUP, "\n### %s | High-level settings\n", eqname);
  cs_log_printf(CS_LOG_SETUP, "  * %s | Type: ", eqname);
  switch (eqp->type) {
  case CS_EQUATION_TYPE_GROUNDWATER:
    cs_log_printf(CS_LOG_SETUP, "Associated to groundwater flows\n");
    break;
  case CS_EQUATION_TYPE_MAXWELL:
    cs_log_printf(CS_LOG_SETUP, "Associated to the Maxwell module\n");
    break;
  case CS_EQUATION_TYPE_NAVSTO:
    cs_log_printf(CS_LOG_SETUP, "Associated to the Navier-Stokes system\n");
    break;
  case CS_EQUATION_TYPE_PREDEFINED:
    cs_log_printf(CS_LOG_SETUP, "Predefined\n");
    break;
  case CS_EQUATION_TYPE_SOLIDIFICATION:
    cs_log_printf(CS_LOG_SETUP, "Associated to the solidification module\n");
    break;
  case CS_EQUATION_TYPE_THERMAL:
    cs_log_printf(CS_LOG_SETUP, "Associated to the thermal module\n");
    break;
  case CS_EQUATION_TYPE_USER:
    cs_log_printf(CS_LOG_SETUP, "User-defined\n");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Eq. %s has no type.\n Please check your settings.", eqname);
  }

  bool  unsteady = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;
  bool  convection = (eqp->flag & CS_EQUATION_CONVECTION) ? true : false;
  bool  diffusion = (eqp->flag & CS_EQUATION_DIFFUSION) ? true : false;
  bool  curlcurl = (eqp->flag & CS_EQUATION_CURLCURL) ? true : false;
  bool  graddiv = (eqp->flag & CS_EQUATION_GRADDIV) ? true : false;
  bool  reaction = (eqp->flag & CS_EQUATION_REACTION) ? true : false;
  bool  source_term = (eqp->n_source_terms > 0) ? true : false;
  bool  force_values = (eqp->flag & CS_EQUATION_FORCE_VALUES) ? true : false;
  bool  inside_system = (eqp->flag & CS_EQUATION_INSIDE_SYSTEM) ? true : false;

  // General information
  // ===================

  cs_log_printf(CS_LOG_SETUP,
                "  * %s | Terms: unsteady:%s, convection:%s, diffusion:%s\n",
                eqname, cs_base_strtf(unsteady), cs_base_strtf(convection),
                cs_base_strtf(diffusion));
  cs_log_printf(CS_LOG_SETUP,
                "  * %s | Terms: curl-curl:%s, grad-div:%s\n",
                eqname, cs_base_strtf(curlcurl), cs_base_strtf(graddiv));
  cs_log_printf(CS_LOG_SETUP,
                "  * %s | Terms: reaction:%s, source term:%s,"
                " force internal values: %s\n",
                eqname,cs_base_strtf(reaction), cs_base_strtf(source_term),
                cs_base_strtf(force_values));

  if (inside_system)
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Attached to a coupled system of equations\n",
                  eqname);

  if (eqp->space_scheme < CS_SPACE_N_SCHEMES)
    cs_log_printf(CS_LOG_SETUP, "  * %s | Space scheme:       %s\n",
                  eqname, cs_param_get_space_scheme_name(eqp->space_scheme));
  else
    bft_error(__FILE__, __LINE__, 0,
              " Undefined space scheme for eq. %s", eqname);

  cs_log_printf(CS_LOG_SETUP, "  * %s | Space poly degree:  %d\n",
                eqname, eqp->space_poly_degree);
  cs_log_printf(CS_LOG_SETUP, "  * %s | Verbosity:          %d\n",
                eqname, eqp->verbosity);

  // Boundary conditions
  // ===================

  cs_log_printf(CS_LOG_SETUP, "\n### %s | Boundary condition settings\n",
                eqname);
  cs_log_printf(CS_LOG_SETUP,
                "  * %s | Boundary conditions | Default: %s\n",
                eqname, cs_param_get_bc_name(eqp->default_bc));
  cs_log_printf(CS_LOG_SETUP,
                "  * %s | Boundary conditions | Enforcement: %s\n",
                eqname,
                cs_param_get_bc_enforcement_name(eqp->default_enforcement));

  if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED)
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Boundary conditions | Strong penalization coeff.:"
                  " %5.3e\n", eqname, eqp->strong_pena_bc_coeff);
  else if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
           eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Boundary conditions | Weak penalization coeff.:"
                  " %5.3e\n", eqname, eqp->weak_pena_bc_coeff);

  cs_log_printf(CS_LOG_SETUP, "  * %s | Boundary conditions |"
                " Number of definitions: %d\n", eqname, eqp->n_bc_defs);

  if (eqp->verbosity > 0) {
    char  desc[128];
    for (int id = 0; id < eqp->n_bc_defs; id++) {
      const cs_xdef_t  *d = eqp->bc_defs[id];

      cs_cdo_bc_get_desc(d->meta, desc);
      sprintf(prefix, "        Definition %3d", id);
      cs_log_printf(CS_LOG_SETUP, "\n%s | Type: %s\n", prefix, desc);
      cs_xdef_log_setup(prefix, d);
    }
  }

  // Time settings
  // =============

  if (unsteady) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Time settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Initial conditions | Number of definitions: %d",
                  eqname, eqp->n_ic_defs);
    if (eqp->n_ic_defs > 0)
      cs_log_printf(CS_LOG_SETUP, "\n\n");
    for (int i = 0; i < eqp->n_ic_defs; i++) {
      sprintf(prefix, "        Definition %3d", i);
      cs_xdef_log_setup(prefix, eqp->ic_defs[i]);
    }

    const char  *time_scheme = cs_param_get_time_scheme_name(eqp->time_scheme);
    if (time_scheme != nullptr) {
      cs_log_printf(CS_LOG_SETUP, "\n  * %s | Time scheme: %s",
                    eqname, time_scheme);
      if (eqp->time_scheme == CS_TIME_SCHEME_THETA)
        cs_log_printf(CS_LOG_SETUP, " with value %f\n", eqp->theta);
      else
        cs_log_printf(CS_LOG_SETUP, "\n");
    }
    else
      bft_error(__FILE__, __LINE__, 0, " Invalid time scheme.");

    cs_log_printf(CS_LOG_SETUP, "  * %s | Mass.Lumping: %s\n",
                  eqname, cs_base_strtf(eqp->do_lumping));
    cs_log_printf(CS_LOG_SETUP, "  * %s | Time property: %s\n\n",
                  eqname, cs_property_get_name(eqp->time_property));

    sprintf(prefix, "        Time Hodge op. ");
    cs_hodge_param_log(prefix, eqp->time_property, eqp->time_hodgep);

  } // Unsteady term

  // Diffusion term
  // ==============

  if (diffusion) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Diffusion term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Diffusion property: %s\n\n",
                  eqname, cs_property_get_name(eqp->diffusion_property));

    sprintf(prefix, "        Diffusion Hodge op. ");
    cs_hodge_param_log(prefix, eqp->diffusion_property, eqp->diffusion_hodgep);

  }

  // curl-curl term
  // ==============

  if (curlcurl) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Curl-Curl term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Curl-Curl property: %s\n\n",
                  eqname, cs_property_get_name(eqp->curlcurl_property));

    sprintf(prefix, "        Curl-curl Hodge op. ");
    cs_hodge_param_log(prefix, eqp->curlcurl_property, eqp->curlcurl_hodgep);

  }

  // grad-div term
  // =============

  if (graddiv) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Grad-Div term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Grad-Div property: %s\n\n",
                  eqname, cs_property_get_name(eqp->graddiv_property));

    sprintf(prefix, "        Grad-Div Hodge op. ");
    cs_hodge_param_log(prefix, eqp->graddiv_property, eqp->graddiv_hodgep);

  }

  // Advection term
  // ==============

  if (convection) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Advection term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Field: \"%s\"\n",
                  eqname, cs_advection_field_get_name(eqp->adv_field));
    if (eqp->adv_scaling_property != nullptr)
      cs_log_printf(CS_LOG_SETUP, "  * %s | Scaling.Property: %s\n",
                    eqname, cs_property_get_name(eqp->adv_scaling_property));

    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Formulation: %s\n",
                  eqname,
                  cs_param_get_advection_form_name(eqp->adv_formulation));

    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Scheme: %s\n",
                  eqname,
                  cs_param_get_advection_scheme_name(eqp->adv_scheme));

    /* Piece of information specific to a scheme */

    if (eqp->adv_scheme == CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND)
      cs_log_printf(CS_LOG_SETUP, "  * %s | Upwind.Portion: %3.2f %%\n",
                    eqname, 100*eqp->upwind_portion);
    else if (eqp->adv_scheme == CS_PARAM_ADVECTION_SCHEME_CIP ||
             eqp->adv_scheme == CS_PARAM_ADVECTION_SCHEME_CIP_CW)
      cs_log_printf(CS_LOG_SETUP, "  * %s | CIP.coef: %.2e\n",
                    eqname, eqp->cip_scaling_coef);

    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Strategy: %s\n",
                  eqname,
                  cs_param_get_advection_strategy_name(eqp->adv_strategy));
    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Extrapolation: %s\n",
                  eqname,
                  cs_param_get_advection_extrapol_name(eqp->adv_extrapol));

  }

  // Reaction term(s)
  // ================

  if (reaction) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Reaction settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Reaction | Number of terms: %d\n",
                  eqname, eqp->n_reaction_terms);

    sprintf(prefix, "        Reaction Hodge op. ");
    cs_hodge_param_log(prefix, nullptr, eqp->reaction_hodgep);
  }

  // Source term(s)
  // ==============

  if (source_term) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Source term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Source terms | Number of terms: %d\n",
                  eqname,
                  eqp->n_source_terms);

    for (int s_id = 0; s_id < eqp->n_source_terms; s_id++) {
      sprintf(prefix, "        Definition %3d", s_id);
      cs_xdef_log_setup(prefix, eqp->source_terms[s_id]);
    }
  }

  // Interior enforcement(s)
  // =======================

  cs_log_printf(CS_LOG_SETUP,
                "\n### %s | Number of interior enforcements: %d\n",
                eqname,
                eqp->n_enforcements);

  if (eqp->n_enforcements > 0) {

    for (int i = 0; i < eqp->n_enforcements; i++)
      cs_enforcement_param_log(eqname, eqp->enforcement_params[i]);
  }

  // Non-linear algorithm
  // ====================

  if (eqp->incremental_algo_type != CS_PARAM_NL_ALGO_NONE) {

    cs_log_printf(CS_LOG_SETUP,
                  "### %s | Incremental algo: %s\n",
                  eqname,
                  cs_param_get_nl_algo_name(eqp->incremental_algo_type));
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Tolerances of the incremental algo:"
                  " rtol: %5.3e; atol: %5.3e; dtol: %5.3e\n",
                  eqname,
                  eqp->incremental_algo_cvg.rtol,
                  eqp->incremental_algo_cvg.atol,
                  eqp->incremental_algo_cvg.dtol);
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Max of non-linear iterations: %d\n",
                  eqname,
                  eqp->incremental_algo_cvg.n_max_iter);
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Relaxation factor: %.3f\n",
                  eqname,
                  eqp->incremental_relax_factor);

    if (eqp->incremental_algo_type == CS_PARAM_NL_ALGO_ANDERSON) {

      const cs_iter_algo_param_aac_t aap = eqp->incremental_anderson_param;

      cs_log_printf(CS_LOG_SETUP,
                    "  * %s | Anderson param: max. dir: %d; "
                    " start: %d; drop. tol: %5.3e; relax: %5.3e\n",
                    eqname,
                    aap.n_max_dir,
                    aap.starting_iter,
                    aap.max_cond,
                    aap.beta);
      cs_log_printf(CS_LOG_SETUP,
                    "  * %s | Anderson param: Dot product type: %s\n",
                    eqname,
                    cs_param_get_dotprod_type_name(aap.dp_type));
    }
  }

  // SLES information
  // ================

  cs_param_saddle_log(eqp->saddle_param);
  cs_param_sles_log(eqp->sles_param);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Ask if the parameter settings of the equation has requested the
 *        treatment of Robin boundary conditions
 *
 * \param[in] eqp  pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_param_has_robin_bc(const cs_equation_param_t *eqp)
{
  if (eqp == nullptr)
    return false;

  for (int i = 0; i < eqp->n_bc_defs; i++) {
    cs_xdef_t *def = eqp->bc_defs[i];
    if (def->meta & CS_CDO_BC_ROBIN)
      return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the initial condition for the unknown related to this
 *        equation. This definition applies to a volume zone.
 *        By default, the unknown is set to zero everywhere.
 *        Here a constant value is set to all the unknows belonging to the
 *        given zone with name z_name
 *        If a definition already exists for the given zone name, then the
 *        previous definition is replaced.
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 * \param[in]      z_name  name of the associated zone (if null or "" all
 *                         cells are considered)
 * \param[in]      val     pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_value(cs_equation_param_t *eqp,
                            const char          *z_name,
                            cs_real_t           *val)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);
  int def_id = _find_or_add_ic(eqp, z_id);

  cs_flag_t meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                       eqp->dim,
                                       z_id,
                                       CS_FLAG_STATE_UNIFORM, /* state flag */
                                       meta_flag,
                                       val);

  /* Update the list of initial conditions for the equation parameters */

  eqp->ic_defs[def_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the initial condition for the unknown related to this
 *        equation. This definition applies to a volume zone.
 *        By default, the unknown is set to zero everywhere.
 *        Here the value set to each unknown belonging to the given zone with
 *        name z_name is such that the integral over the cells of the zone
 *        returns the requested quantity
 *        If a definition already exists for the given zone name, then the
 *        previous definition is replaced.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if null or
 *                           "" all cells are considered)
 * \param[in]      quantity  quantity to distribute over the mesh location
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_qov(cs_equation_param_t *eqp,
                          const char          *z_name,
                          double               quantity)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);
  int def_id = _find_or_add_ic(eqp, z_id);

  cs_flag_t meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t *d = cs_xdef_volume_create(CS_XDEF_BY_QOV,
                                       eqp->dim,
                                       z_id,
                                       0, /* state flag */
                                       meta_flag,
                                       &quantity);

  /* Update the list of initial conditions for the equation parameters */

  eqp->ic_defs[def_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the initial condition for the unknown related to this
 *        equation. This definition applies to a volume zone.
 *        By default, the unknown is set to zero everywhere.
 *        Here the initial value for each unknown associated to the zone with
 *        name z_name is set according to an analytical function
 *        If a definition already exists for the given zone name, then the
 *        previous definition is replaced.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if null or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     null or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_analytic(cs_equation_param_t *eqp,
                               const char          *z_name,
                               cs_analytic_func_t  *analytic,
                               void                *input)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);
  int def_id = _find_or_add_ic(eqp, z_id);

  cs_flag_t meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_context_t ac
    = { .z_id = z_id, .func = analytic, .input = input, .free_input = nullptr};

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim, z_id,
                                        0, /* state flag */
                                        meta_flag,
                                        &ac);

  /* Update the list of initial conditions for the equation parameters */

  eqp->ic_defs[def_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define the initial condition for the unknown related to this
 *        equation. This definition applies to a volume zone.
 *        By default, the unknown is set to zero everywhere.
 *        Case of a definition by a DoF function.
 *        If a definition already exists for the given zone name, then the
 *        previous definition is replaced.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  where information is computed
 * \param[in]      func      pointer to a DoF function
 * \param[in]      input     nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_dof_func(cs_equation_param_t *eqp,
                               const char          *z_name,
                               cs_flag_t            loc_flag,
                               cs_dof_func_t       *func,
                               void                *input)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);
  int def_id = _find_or_add_ic(eqp, z_id);

  /* Define a flag according to the kind of space discretization */

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_dof_context_t context = {
    .z_id         = z_id,
    .dof_location = loc_flag,
    .func         = func,
    .input        = input,
    .free_input   = nullptr,
  };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_DOF_FUNCTION,
                                        eqp->dim, z_id,
                                        0, /* state flag */
                                        meta_flag,
                                        &context);

  /* Update the list of initial conditions for the equation parameters */

  eqp->ic_defs[def_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set a boundary condition from an existing \ref cs_xdef_t structure
 *        The lifecycle of the cs_xdef_t structure is now managed by the
 *        current \ref cs_equation_param_t structure.
 *
 * \param[in, out] eqp   pointer to a cs_equation_param_t structure
 * \param[in]      xdef  pointer to the \ref cs_xdef_t structure to transfer
*/
/*----------------------------------------------------------------------------*/

void
cs_equation_add_xdef_bc(cs_equation_param_t *eqp,
                        cs_xdef_t           *xdef)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = xdef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and initialize a new structure to set a boundary condition
 *        related to the given equation structure
 *        z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out] eqp      pointer to a cs_equation_param_t structure
 * \param[in]      bc_type  type of boundary condition to add
 * \param[in]      z_name   name of the related boundary zone
 * \param[in]      values   pointer to a array storing the values
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_value(cs_equation_param_t      *eqp,
                            const cs_param_bc_type_t  bc_type,
                            const char               *z_name,
                            cs_real_t                *values)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  if (   eqp->space_scheme != CS_SPACE_SCHEME_LEGACY
      && bc_type == CS_BC_WALL_MODELLED)
    bft_error(__FILE__, __LINE__, 0, "%s: To be done.\n", __func__);

  /* Add a new cs_xdef_t structure */

  int  dim = eqp->dim;
  if (bc_type == CS_BC_NEUMANN_FULL)
    dim *= 3;  /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_BC_ROBIN) {

    /* FluxDiff = alpha * (u - u0) + beta => Set (alpha, u0, beta) */

    if (eqp->dim == 1)
      dim = 3;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled yet.\n", __func__);

  }

  cs_flag_t  meta_flag = (eqp-> space_scheme == CS_SPACE_SCHEME_LEGACY) ?
    (cs_flag_t)bc_type : cs_cdo_bc_get_flag(bc_type);

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          dim,
                                          cs_boundary_zone_id_by_name(z_name),
                                          CS_FLAG_STATE_UNIFORM, /* state */
                                          meta_flag,
                                          (void *)values);

  /* Before incrementing the list of definitions, first verify that there
   * isn't an existing one on the same zone. If so, remove it.
   * This is done at the end of the function in order to ensure that it was
   * possible to create the new definition.
   */

  cs_equation_remove_bc(eqp, z_name);

  /* Increment list of boundary conditions for the equation parameters */

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and initialize a new structure to set a boundary condition
 *        related to the given equation structure
 *        z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out] eqp          pointer to a cs_equation_param_t structure
 * \param[in]      bc_type      type of boundary condition to add
 * \param[in]      z_name       name of the related boundary zone
 * \param[in]      loc          information to know where are located values
 * \param[in]      array        pointer to an array
 * \param[in]      is_owner     transfer the lifecycle to the cs_xdef_t struct.
 *                              (true or false)
 * \param[in]      full_length  if true, size of "array" should be allocated
 *                              to the total numbers of entities related to the
 *                              given location. If false, a new list is
 *                              allocated and filled with the related subset
 *                              indirection.
 *
 * \return a pointer to the new allocated \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_array(cs_equation_param_t      *eqp,
                            const cs_param_bc_type_t  bc_type,
                            const char               *z_name,
                            cs_flag_t                 loc,
                            cs_real_t                *array,
                            bool                      is_owner,
                            bool                      full_length)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  if (   eqp->space_scheme != CS_SPACE_SCHEME_LEGACY
      && bc_type == CS_BC_WALL_MODELLED)
    bft_error(__FILE__, __LINE__, 0, "%s: To be done.\n", __func__);

  assert(cs_flag_test(loc, cs_flag_primal_face)   ||
         cs_flag_test(loc, cs_flag_boundary_face) ||
         cs_flag_test(loc, cs_flag_primal_vtx)    ||
         cs_flag_test(loc, cs_flag_primal_edge)); /* for circulation */

  int  z_id = cs_boundary_zone_id_by_name(z_name);
  int  dim = eqp->dim;
  cs_flag_t  state_flag = 0;

  if (loc == cs_flag_primal_face || loc == cs_flag_boundary_face)
    state_flag = CS_FLAG_STATE_FACEWISE;

  if (bc_type == CS_BC_NEUMANN_FULL)
    dim *= 3;  /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_BC_ROBIN) {

    /* FluxDiff = alpha * (u - u0) + beta => Set (alpha, u0, beta) */

    if (eqp->dim == 1)
      dim = 3;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled yet.\n", __func__);

  }

  cs_flag_t  meta_flag = (eqp-> space_scheme == CS_SPACE_SCHEME_LEGACY) ?
    (cs_flag_t)bc_type : cs_cdo_bc_get_flag(bc_type);

  /* Add a new cs_xdef_t structure */

  cs_xdef_array_context_t input = {
    .z_id           = z_id,
    .stride         = dim,
    .value_location = loc,
    .is_owner       = is_owner,
    .full_length    = full_length,
    .values         = array,
    /* Optional parameters */
    .full2subset = nullptr,
    .n_list_elts = 0,
    .elt_ids     = nullptr,
    .adjacency   = nullptr,
  };

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ARRAY,
                                          dim,
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          (void *)&input);

  /* Build the indirection array if only a subset is used */

  if (!full_length)
    cs_xdef_array_build_full2subset(d);

  /* Before incrementing the list of definitions, first verify that there
   * isn't an existing one on the same zone. If so, remove it.
   * This is done at the end of the function in order to ensure that it was
   * possible to create the new definition.
   */

  cs_equation_remove_bc(eqp, z_name);

  /* Increment list of boundary conditions for the equation parameters */
  int  new_id = eqp->n_bc_defs;

  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and initialize a new structure to set a boundary condition
 *        related to the given equation structure
 *        z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out] eqp      pointer to a cs_equation_param_t structure
 * \param[in]      bc_type  type of boundary condition to add
 * \param[in]      z_name   name of the related boundary zone
 * \param[in]      field    pointer to a cs_field_t structure
 *
 * \return a pointer to the new allocated \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_field(cs_equation_param_t      *eqp,
                            const cs_param_bc_type_t  bc_type,
                            const char               *z_name,
                            cs_field_t               *field)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  int  z_id = cs_boundary_zone_id_by_name(z_name);

  int dim = eqp->dim;
  if (bc_type == CS_BC_NEUMANN_FULL)
    dim *= 3;  /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_BC_ROBIN) {

    /* FluxDiff = alpha * (u - u0) + beta => Set (alpha, u0, beta) */

    if (eqp->dim == 1)
      dim = 3;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled yet.\n", __func__);

  }

  if (   eqp->space_scheme != CS_SPACE_SCHEME_LEGACY
      && bc_type == CS_BC_WALL_MODELLED)
    bft_error(__FILE__, __LINE__, 0, "%s: To be done.\n", __func__);

  assert(field != nullptr);
  if (dim != field->dim)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid dimension for field %s\n", __func__, field->name);

  cs_flag_t  state_flag = CS_FLAG_STATE_FACEWISE;
  cs_flag_t  meta_flag = (eqp-> space_scheme == CS_SPACE_SCHEME_LEGACY) ?
    (cs_flag_t)bc_type : cs_cdo_bc_get_flag(bc_type);

  /* Add a new cs_xdef_t structure */

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_FIELD,
                                          dim,
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          field);

  /* Before incrementing the list of definitions, first verify that there
   * isn't an existing one on the same zone. If so, remove it.
   * This is done at the end of the function in order to ensure that it was
   * possible to create the new definition.
   */

  cs_equation_remove_bc(eqp, z_name);

  /* Increment list of boundary conditions for the equation parameters */

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and initialize a new structure to set a boundary condition
 *        related to the given equation param structure
 *        ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      bc_type   type of boundary condition to add
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function defining the value
 * \param[in]      input     nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_analytic(cs_equation_param_t      *eqp,
                               const cs_param_bc_type_t  bc_type,
                               const char               *z_name,
                               cs_analytic_func_t       *analytic,
                               void                     *input)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Set the value for dim */

  int dim = eqp->dim;

  if (bc_type == CS_BC_NEUMANN_FULL)
    dim *= 3; /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_BC_CIRCULATION) {

    /* This is a vector-valued equation but the DoF is scalar-valued since
     * it is a circulation associated to each edge */

    if (eqp->dim == 3)
      dim = 1;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled.\n", __func__);
  }

  if (bc_type == CS_BC_ROBIN) {

    /* FluxDiff = alpha * (u - u0) + beta => Set (alpha, u0, beta) */

    if (eqp->dim == 1)
      dim = 3;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled yet.\n", __func__);
  }

  if (   eqp->space_scheme != CS_SPACE_SCHEME_LEGACY
      && bc_type == CS_BC_WALL_MODELLED)
    bft_error(__FILE__, __LINE__, 0, "%s: To be done.\n", __func__);

  int z_id = cs_boundary_zone_id_by_name(z_name);

  /* Add a new cs_xdef_t structure */

  cs_xdef_analytic_context_t ac
    = { .z_id = z_id, .func = analytic, .input = input, .free_input = nullptr };

  cs_flag_t  meta_flag = (eqp-> space_scheme == CS_SPACE_SCHEME_LEGACY) ?
    (cs_flag_t)bc_type : cs_cdo_bc_get_flag(bc_type);

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          dim,
                                          z_id,
                                          0, /* state */
                                          meta_flag,
                                          &ac);

  /* Before incrementing the list of definitions, first verify that there
   * isn't an existing one on the same zone. If so, remove it.
   * This is done at the end of the function in order to ensure that it was
   * possible to create the new definition.
   */

  cs_equation_remove_bc(eqp, z_name);

  /* Increment list of boundary conditions for the equation parameters */

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and initialize a new structure to set a boundary condition
 *        related to the given equation param structure
 *        ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *        Definition relying on a \ref cs_time_func_t function pointer
 *
 * \param[in, out] eqp      pointer to a cs_equation_param_t structure
 * \param[in]      bc_type  type of boundary condition to add
 * \param[in]      z_name   name of the associated zone (if nullptr or "" if
 *                          all cells are considered)
 * \param[in]      t_func   pointer to an analytic function defining the value
 * \param[in]      input    nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_time_func(cs_equation_param_t      *eqp,
                                const cs_param_bc_type_t  bc_type,
                                const char               *z_name,
                                cs_time_func_t           *t_func,
                                void                     *input)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Set the value for dim */

  int dim = eqp->dim;

  if (bc_type == CS_BC_NEUMANN_FULL)
    dim *= 3; /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_BC_CIRCULATION) {

    /* This is a vector-valued equation but the DoF is scalar-valued since
     * it is a circulation associated to each edge */

    if (eqp->dim == 3)
      dim = 1;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled.\n", __func__);
  }

  if (bc_type == CS_BC_ROBIN) {

    /* FluxDiff = alpha * (u - u0) + beta => Set (alpha, u0, beta) */

    if (eqp->dim == 1)
      dim = 3;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled yet.\n", __func__);
  }

  if (   eqp->space_scheme != CS_SPACE_SCHEME_LEGACY
      && bc_type == CS_BC_WALL_MODELLED)
    bft_error(__FILE__, __LINE__, 0, "%s: To be done.\n", __func__);

  int z_id = cs_boundary_zone_id_by_name(z_name);

  /* Add a new cs_xdef_t structure */

  cs_xdef_time_func_context_t tfc
    = { .z_id = z_id, .func = t_func, .input = input, .free_input = nullptr };

  cs_flag_t  meta_flag = (eqp-> space_scheme == CS_SPACE_SCHEME_LEGACY) ?
    (cs_flag_t)bc_type : cs_cdo_bc_get_flag(bc_type);

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_TIME_FUNCTION,
                                          dim,
                                          z_id,
                                          0, /* state */
                                          meta_flag,
                                          &tfc);

  /* Before incrementing the list of definitions, first verify that there
   * isn't an existing one on the same zone. If so, remove it.
   * This is done at the end of the function in order to ensure that it was
   * possible to create the new definition.
   */

  cs_equation_remove_bc(eqp, z_name);

  /* Increment list of boundary conditions for the equation parameters */

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and initialize a new structure to set a boundary condition
 *        related to the given cs_equation_param_t structure
 *        ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      bc_type   type of boundary condition to add
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  location where values are computed
 * \param[in]      func      pointer to cs_dof_func_t function
 * \param[in]      input     nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_dof_func(cs_equation_param_t      *eqp,
                               const cs_param_bc_type_t  bc_type,
                               const char               *z_name,
                               cs_flag_t                 loc_flag,
                               cs_dof_func_t            *func,
                               void                     *input)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Set the value for dim */

  int dim = eqp->dim;

  if (bc_type == CS_BC_NEUMANN_FULL)
    dim *= 3; /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_BC_CIRCULATION) {

    /* This is a vector-valued equation but the DoF is scalar-valued since
     * it is a circulation associated to each edge */

    if (eqp->dim == 3)
      dim = 1;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled.\n", __func__);
  }

  if (bc_type == CS_BC_ROBIN) {

    /* FluxDiff = alpha * (u - u0) + beta => Set (alpha, u0, beta) */

    if (eqp->dim == 1)
      dim = 3;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled yet.\n", __func__);
  }

  if (   eqp->space_scheme != CS_SPACE_SCHEME_LEGACY
      && bc_type == CS_BC_WALL_MODELLED)
    bft_error(__FILE__, __LINE__, 0, "%s: To be done.\n", __func__);

  int z_id = cs_boundary_zone_id_by_name(z_name);

  /* Add a new cs_xdef_t structure */

  cs_xdef_dof_context_t ctx = { .z_id         = z_id,
                                .dof_location = loc_flag,
                                .func         = func,
                                .input        = input,
                                .free_input   = nullptr };

  cs_flag_t  meta_flag = (eqp-> space_scheme == CS_SPACE_SCHEME_LEGACY) ?
    (cs_flag_t)bc_type : cs_cdo_bc_get_flag(bc_type);

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_DOF_FUNCTION,
                                          dim,
                                          z_id,
                                          0, /* state */
                                          meta_flag,
                                          &ctx);

  /* Before incrementing the list of definitions, first verify that there
   * isn't an existing one on the same zone. If so, remove it.
   * This is done at the end of the function in order to ensure that it was
   * possible to create the new definition.
   */

  cs_equation_remove_bc(eqp, z_name);

  /* Increment list of boundary conditions for the equation parameters */

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to existing boundary condition definition structure
 *        for the given equation param structure and zone.
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 * \param[in]      z_name  name of the associated zone (if nullptr or "" if
 *                         all cells are considered)
 *
 * \return a pointer to the \ref cs_xdef_t structure if present, or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_find_bc(cs_equation_param_t *eqp,
                    const char          *z_name)
{
  if (eqp == nullptr)
    return nullptr;

  int z_id = -2;

  if (z_name != nullptr) {
    const cs_zone_t *z = cs_boundary_zone_by_name_try(z_name);
    if (z != nullptr)
      z_id = z->id;
  }

  /* Search for given BC */

  for (int i = 0; i < eqp->n_bc_defs; i++) {
    if (eqp->bc_defs[i]->z_id == z_id) {
      return eqp->bc_defs[i];
    }
  }

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove boundary condition from the given equation param structure
 *        for a given zone.
 *
 * If no matching boundary condition is found, the function returns
 * silently.
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 * \param[in]      z_name  name of the associated zone (if nullptr or "" if
 *                         all cells are considered)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_remove_bc(cs_equation_param_t *eqp,
                      const char          *z_name)
{
  if (eqp == nullptr)
    return;

  int z_id = -2;

  if (z_name != nullptr) {
    const cs_zone_t *z = cs_boundary_zone_by_name_try(z_name);
    if (z != nullptr)
      z_id = z->id;
  }

  /* Search for given BC */

  int j = -1;
  for (int i = 0; i < eqp->n_bc_defs; i++) {
    if (eqp->bc_defs[i]->z_id == z_id) {
      j = i;
      break;
    }
  }

  /* Remove it if found */

  if (j > -1) {
    eqp->bc_defs[j] = cs_xdef_free(eqp->bc_defs[j]);
    for (int i = j + 1; i < eqp->n_bc_defs; i++) {
      eqp->bc_defs[i - 1] = eqp->bc_defs[i];
    }
    eqp->n_bc_defs -= 1;
    BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Remove initial condition from the given equation param structure
 *        for a given zone.
 *
 * If no matching boundary condition is found, the function returns
 * silently.
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 * \param[in]      z_name  name of the associated zone (if nullptr or "" if
 *                         all cells are considered)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_remove_ic(cs_equation_param_t *eqp,
                      const char          *z_name)
{
  if (eqp == nullptr)
    return;

  int z_id = -2;

  if (z_name != nullptr) {
    const cs_zone_t *z = cs_volume_zone_by_name_try(z_name);
    if (z != nullptr)
      z_id = z->id;
  }

  /* Search for given IC */

  int j = -1;
  for (int i = 0; i < eqp->n_ic_defs; i++) {
    if (eqp->ic_defs[i]->z_id == z_id) {
      j = i;
      break;
    }
  }

  /* Remove it if found */

  if (j > -1) {
    eqp->ic_defs[j] = cs_xdef_free(eqp->ic_defs[j]);
    for (int i = j + 1; i < eqp->n_ic_defs; i++) {
      eqp->ic_defs[i - 1] = eqp->ic_defs[i];
    }
    eqp->n_ic_defs -= 1;
    BFT_REALLOC(eqp->ic_defs, eqp->n_ic_defs, cs_xdef_t *);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define and initialize a new structure to set a sliding boundary
 *        condition related to the given equation structure
 *        z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 * \param[in]      z_name  name of the related boundary zone
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_sliding_condition(cs_equation_param_t *eqp,
                                  const char          *z_name)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (eqp->dim < 3)
    bft_error(
      __FILE__, __LINE__, 0, "%s: Invalid dimension of equation\n", __func__);

  /* Add two definitions: one for the normal component and one for the
     tangential component */

  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs + 1, cs_xdef_t *);

  cs_xdef_t *d   = nullptr;
  cs_real_t  val = 0;

  /* Add the homogeneous Dirichlet on the normal component */

  d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                              1,
                              cs_boundary_zone_id_by_name(z_name),
                              CS_FLAG_STATE_UNIFORM,  /* state flag */
                              CS_CDO_BC_SYMMETRY,      /* meta */
                              (void *)&val);

  eqp->bc_defs[eqp->n_bc_defs] = d;
  eqp->n_bc_defs += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a new term related to the Laplacian operator for the
 *        equation associated to the given \ref cs_equation_param_t structure
 *        Laplacian means div-grad (either for vector-valued or scalar-valued
 *        equations)
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      property  pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_diffusion(cs_equation_param_t *eqp,
                          cs_property_t       *property)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  eqp->flag |= CS_EQUATION_DIFFUSION;
  eqp->diffusion_property = property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a new term related to the curl-curl operator for the
 *        equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 * \param[in]      inversion  > 0: true, false otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_curlcurl(cs_equation_param_t *eqp,
                         cs_property_t       *property,
                         int                  inversion)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  eqp->flag |= CS_EQUATION_CURLCURL;
  eqp->curlcurl_property = property;

  if (inversion > 0)
    eqp->curlcurl_hodgep.inv_pty = true;
  else
    eqp->curlcurl_hodgep.inv_pty = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a new term related to the grad-div operator for the
 *        equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      property  pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_graddiv(cs_equation_param_t *eqp,
                        cs_property_t       *property)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  eqp->flag |= CS_EQUATION_GRADDIV;
  eqp->graddiv_property = property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a new term related to the time derivative operator for
 *        the equation associated to the given \ref cs_equation_param_t
 *        structure
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      property  pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_time(cs_equation_param_t *eqp,
                     cs_property_t       *property)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  eqp->flag |= CS_EQUATION_UNSTEADY;
  eqp->time_property = property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a new term related to the advection operator for the
 *        equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      adv_field  pointer to a cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_advection(cs_equation_param_t *eqp,
                          cs_adv_field_t      *adv_field)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (adv_field == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: Stop adding an empty advection field.",
              __func__, eqp->name);

  eqp->flag |= CS_EQUATION_CONVECTION;
  eqp->adv_field = adv_field;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a scaling property to the advection
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      property  pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_advection_scaling_property(cs_equation_param_t *eqp,
                                           cs_property_t       *property)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  eqp->adv_scaling_property = property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associate a new term related to the reaction operator for the
 *        equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      property  pointer to a cs_property_t structure
 *
 * \return the id related to the reaction term
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_add_reaction(cs_equation_param_t *eqp,
                         cs_property_t       *property)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  /* Only this kind of reaction term is available up to now.
     Add a new reaction term */

  int  new_id = eqp->n_reaction_terms;
  eqp->n_reaction_terms += 1;
  BFT_REALLOC(eqp->reaction_properties, eqp->n_reaction_terms, cs_property_t *);
  eqp->reaction_properties[new_id] = property;

  /* Flag the equation with "reaction" */

  eqp->flag |= CS_EQUATION_REACTION;

  return new_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new source term by initializing a cs_xdef_t structure.
 *        Case of a definition by a constant value
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 * \param[in]      z_name  name of the associated zone (if nullptr or "" all
 *                         cells are considered)
 * \param[in]      val     pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_val(cs_equation_param_t *eqp,
                                   const char          *z_name,
                                   cs_real_t           *val)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);

  /* Define a flag according to the kind of space discretization */

  cs_flag_t  state_flag = 0, meta_flag = 0;
  cs_source_term_set_default_flag(eqp->space_scheme, &state_flag, &meta_flag);
  state_flag |= CS_FLAG_STATE_UNIFORM;

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        (void *)val);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new source term by initializing a cs_xdef_t structure.
 *        Case of a definition by an analytical function
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 * \param[in]      z_name  name of the associated zone (if nullptr or "" if
 *                         all cells are considered)
 * \param[in]      func    pointer to an analytical function
 * \param[in]      input   nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_analytic(cs_equation_param_t *eqp,
                                        const char          *z_name,
                                        cs_analytic_func_t  *func,
                                        void                *input)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Define a flag according to the kind of space discretization */

  cs_flag_t  state_flag = 0, meta_flag = 0;
  cs_source_term_set_default_flag(eqp->space_scheme, &state_flag, &meta_flag);

  int z_id = cs_volume_zone_id_by_name(z_name);
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_context_t ac
    = { .z_id = z_id, .func = func, .input = input, .free_input = nullptr };

  /* Add a new cs_xdef_t structure */

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &ac);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new source term by initializing a cs_xdef_t structure.
 *        Case of a definition by a DoF function
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  location of element ids given as parameter
 * \param[in]      func      pointer to a DoF function
 * \param[in]      input     nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_dof_func(cs_equation_param_t *eqp,
                                        const char          *z_name,
                                        cs_flag_t            loc_flag,
                                        cs_dof_func_t       *func,
                                        void                *input)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);

  /* Define a flag according to the kind of space discretization */

  cs_flag_t  state_flag = 0, meta_flag = 0;
  cs_source_term_set_default_flag(eqp->space_scheme, &state_flag, &meta_flag);

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_dof_context_t context = { .z_id         = z_id,
                                    .dof_location = loc_flag,
                                    .func         = func,
                                    .input        = input,
                                    .free_input   = nullptr };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_DOF_FUNCTION,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &context);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new source term by initializing a cs_xdef_t structure.
 *        Case of a definition by an array.
 *
 * \param[in, out] eqp          pointer to a cs_equation_param_t structure
 * \param[in]      z_name       name of the associated zone (if nullptr or "" if
 *                              all cells are considered)
 * \param[in]      loc          information to know where are located values
 * \param[in]      array        pointer to an array
 * \param[in]      is_owner     transfer the lifecycle to the cs_xdef_t struct.
 *                              (true or false)
 * \param[in]      full_length  if true, the size of "array" should be allocated
 *                              to the total numbers of entities related to the
 *                              given location. If false, a new list is
 *                              allocated and filled with the related subset
 *                              indirection.
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_array(cs_equation_param_t *eqp,
                                     const char          *z_name,
                                     cs_flag_t            loc,
                                     cs_real_t           *array,
                                     bool                 is_owner,
                                     bool                 full_length)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);

  /* Define a flag according to the kind of space discretization */

  cs_flag_t state_flag = 0, meta_flag = 0;
  cs_source_term_set_default_flag(eqp->space_scheme, &state_flag, &meta_flag);

  if (cs_flag_test(loc, cs_flag_primal_vtx) == true)
    state_flag = CS_FLAG_STATE_POTENTIAL; /* erase predefined settings */
  else if (cs_flag_test(loc, cs_flag_primal_cell) == true)
    state_flag |= CS_FLAG_STATE_CELLWISE;

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_array_context_t ctx = {
    .z_id           = z_id,
    .stride         = eqp->dim,
    .value_location = loc,
    .is_owner       = is_owner,
    .full_length    = full_length,
    .values         = array,
    /* Optional parameters */
    .full2subset = nullptr,
    .n_list_elts = 0,
    .elt_ids     = nullptr,
    .adjacency   = nullptr,
  };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ARRAY,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        (void *)&ctx);

  /* Build the indirection array if only a subset is used */

  if (!full_length)
    cs_xdef_array_build_full2subset(d);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new volume mass injection definition source term by
 *        initializing a cs_xdef_t structure, using a constant value.
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 * \param[in]      z_name  name of the associated zone (if nullptr or "" if
 *                         all cells are considered)
 * \param[in]      val     pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_value(cs_equation_param_t *eqp,
                                               const char          *z_name,
                                               double              *val)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);

  cs_flag_t state_flag = 0, meta_flag = 0;

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        val);

  int  new_id = eqp->n_volume_mass_injections;
  eqp->n_volume_mass_injections += 1;
  BFT_REALLOC(eqp->volume_mass_injections,
              eqp->n_volume_mass_injections,
              cs_xdef_t *);
  eqp->volume_mass_injections[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new volume mass injection definition source term by
 *        initializing a cs_xdef_t structure, using a constant quantity
 *        distributed over the associated zone's volume.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
 *                           all cells are considered)
 * \param[in]      quantity  pointer to quantity to distribute over the zone
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_qov(cs_equation_param_t *eqp,
                                             const char          *z_name,
                                             double              *quantity)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);

  cs_flag_t state_flag = 0, meta_flag = 0;

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_QOV,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        quantity);

  int  new_id = eqp->n_volume_mass_injections;
  eqp->n_volume_mass_injections += 1;
  BFT_REALLOC(eqp->volume_mass_injections,
              eqp->n_volume_mass_injections,
              cs_xdef_t *);
  eqp->volume_mass_injections[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new volume mass injection definition source term by
 *        initializing a cs_xdef_t structure, using an analytical function.
 *
 * \param[in, out] eqp     pointer to a cs_equation_param_t structure
 * \param[in]      z_name  name of the associated zone (if nullptr or "" if
 *                         all cells are considered)
 * \param[in]      func    pointer to an analytical function
 * \param[in]      input   nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_analytic(cs_equation_param_t *eqp,
                                                  const char          *z_name,
                                                  cs_analytic_func_t  *func,
                                                  void                *input)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);

  cs_flag_t  state_flag = 0, meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_context_t ctx
    = { .z_id = z_id, .func = func, .input = input, .free_input = nullptr };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &ctx);

  int  new_id = eqp->n_volume_mass_injections;
  eqp->n_volume_mass_injections += 1;
  BFT_REALLOC(eqp->volume_mass_injections,
              eqp->n_volume_mass_injections,
              cs_xdef_t *);
  eqp->volume_mass_injections[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new volume mass injection definition source term by
 *        initializing a cs_xdef_t structure, using a DoF function.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if nullptr or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  where information is computed
 * \param[in]      func      pointer to an analytical function
 * \param[in]      input     nullptr or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_dof_func(cs_equation_param_t *eqp,
                                                  const char          *z_name,
                                                  cs_flag_t            loc_flag,
                                                  cs_dof_func_t       *func,
                                                  void                *input)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */

  int z_id = cs_volume_zone_id_by_name(z_name);

  cs_flag_t  state_flag = 0, meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_dof_context_t ctx = {
    .z_id         = z_id,
    .dof_location = loc_flag,
    .func         = func,
    .input        = input,
    .free_input   = nullptr,
  };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_DOF_FUNCTION,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &ctx);

  int  new_id = eqp->n_volume_mass_injections;
  eqp->n_volume_mass_injections += 1;
  BFT_REALLOC(eqp->volume_mass_injections,
              eqp->n_volume_mass_injections,
              cs_xdef_t *);
  eqp->volume_mass_injections[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new volume mass injection definition source term by
 *        initializing a cs_xdef_t structure, using an array.
 *
 * \param[in, out] eqp          pointer to a cs_equation_param_t structure
 * \param[in]      bc_type      type of boundary condition to add
 * \param[in]      z_name       name of the related boundary zone
 * \param[in]      loc          information to know where are located values
 * \param[in]      array        pointer to an array
 * \param[in]      is_owner     transfer the lifecycle to the cs_xdef_t struct.
 *                              (true or false)
 * \param[in]      full_length  if true, size of "array" should be allocated
 *                              to the total numbers of entities related to the
 *                              given location. If false, a new list is
 *                              allocated and filled with the related subset
 *                              indirection.
 *
 * \return a pointer to the new allocated \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_array(cs_equation_param_t *eqp,
                                               const char          *z_name,
                                               cs_flag_t            loc_flag,
                                               cs_real_t           *array,
                                               bool                 is_owner,
                                               bool                 full_length)
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  int  z_id = cs_volume_zone_id_by_name(z_name);
  int  dim = eqp->dim;

  cs_flag_t  state_flag = 0, meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  /* Add a new cs_xdef_t structure */

  cs_xdef_array_context_t input = {
    .z_id           = z_id,
    .stride         = dim,
    .value_location = loc_flag,
    .is_owner       = is_owner,
    .full_length    = full_length,
    .values         = array,
    /* Optional parameters */
    .full2subset = nullptr,
    .n_list_elts = 0,
    .elt_ids     = nullptr,
    .adjacency   = nullptr,
  };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ARRAY,
                                        dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        (void *)&input);

  /* Build the indirection array if only a subset is used */

  if (!full_length)
    cs_xdef_array_build_full2subset(d);

  int  new_id = eqp->n_volume_mass_injections;
  eqp->n_volume_mass_injections += 1;
  BFT_REALLOC(eqp->volume_mass_injections,
              eqp->n_volume_mass_injections,
              cs_xdef_t *);
  eqp->volume_mass_injections[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an enforcement of the value of degrees of freedom located at
 *        the mesh vertices.
 *        The spatial discretization scheme for the given equation has to be
 *        CDO vertex-based or CDO vertex+cell-based schemes.
 *
 *        One assumes that values are interlaced if eqp->dim > 1
 *        ref_value or elt_values has to be defined. If both parameters are
 *        defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp         pointer to a cs_equation_param_t structure
 * \param[in]      n_vertices  number of vertices to enforce
 * \param[in]      vertex_ids  list of vertices
 * \param[in]      ref_value   default values or ignored (may be nullptr)
 * \param[in]      vtx_values  list of associated values, ignored if nullptr
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_vertex_dof_enforcement(cs_equation_param_t *eqp,
                                       cs_lnum_t            n_vertices,
                                       const cs_lnum_t      vertex_ids[],
                                       const cs_real_t      ref_value[],
                                       const cs_real_t      vtx_values[])
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVB &&
      eqp->space_scheme != CS_SPACE_SCHEME_CDOVCB)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: Invalid space scheme.\n"
              "This should be a vertex-based one.", __func__, eqp->name);
  if (ref_value == nullptr && vtx_values == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: No enforcement value.\n", __func__, eqp->name);

  cs_enforcement_param_t *efp            = nullptr;
  int  enforcement_id = eqp->n_enforcements;

  eqp->n_enforcements++;

  if (vtx_values == nullptr) {

    assert(ref_value != nullptr);
    efp = cs_enforcement_param_create(CS_ENFORCEMENT_SELECTION_VERTICES,
                                      CS_ENFORCEMENT_BY_CONSTANT,
                                      eqp->dim,
                                      n_vertices,
                                      vertex_ids,
                                      ref_value);
  }
  else
    efp = cs_enforcement_param_create(CS_ENFORCEMENT_SELECTION_VERTICES,
                                      CS_ENFORCEMENT_BY_DOF_VALUES,
                                      eqp->dim,
                                      n_vertices,
                                      vertex_ids,
                                      vtx_values);

  BFT_REALLOC(
    eqp->enforcement_params, eqp->n_enforcements, cs_enforcement_param_t *);
  eqp->enforcement_params[enforcement_id] = efp;
  eqp->flag |= CS_EQUATION_FORCE_VALUES;

  return efp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an enforcement of the value of degrees of freedom located at
 *        the mesh edges.
 *        The spatial discretization scheme for the given equation has to be
 *        CDO edge-based schemes.
 *
 *        One assumes that values are interlaced if eqp->dim > 1
 *        ref_value or elt_values has to be defined. If both parameters are
 *        defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp          pointer to a cs_equation_param_t structure
 * \param[in]      n_edges      number of edges to enforce
 * \param[in]      edge_ids     list of edges
 * \param[in]      ref_value    default values or ignored (may be nullptr)
 * \param[in]      edge_values  list of associated values, ignored if nullptr
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_edge_dof_enforcement(cs_equation_param_t *eqp,
                                     cs_lnum_t            n_edges,
                                     const cs_lnum_t      edge_ids[],
                                     const cs_real_t      ref_value[],
                                     const cs_real_t      edge_values[])
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOEB)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: Invalid space scheme.\n"
              "This should be a edge-based one.",
              __func__, eqp->name);
  if (ref_value == nullptr && edge_values == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: No enforcement value.\n",
              __func__, eqp->name);

  cs_enforcement_param_t *efp            = nullptr;
  int  enforcement_id = eqp->n_enforcements;

  eqp->n_enforcements++;

  /* Edge-based schemes are related to a vector-valued equation but DoF are
     scalar-valued. They are circulation. */

  if (edge_values == nullptr) {

    assert(ref_value != nullptr);
    efp = cs_enforcement_param_create(CS_ENFORCEMENT_SELECTION_EDGES,
                                      CS_ENFORCEMENT_BY_CONSTANT,
                                      1,
                                      n_edges,
                                      edge_ids,
                                      ref_value);
  }
  else
    efp = cs_enforcement_param_create(CS_ENFORCEMENT_SELECTION_EDGES,
                                      CS_ENFORCEMENT_BY_DOF_VALUES,
                                      1,
                                      n_edges,
                                      edge_ids,
                                      edge_values);

  BFT_REALLOC(
    eqp->enforcement_params, eqp->n_enforcements, cs_enforcement_param_t *);
  eqp->enforcement_params[enforcement_id] = efp;
  eqp->flag |= CS_EQUATION_FORCE_VALUES;

  return efp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an enforcement of the value of degrees of freedom located at
 *        the mesh faces.
 *        The spatial discretization scheme for the given equation has to be
 *        CDO face-based schemes.
 *
 *        One assumes that values are interlaced if eqp->dim > 1
 *        ref_value or elt_values has to be defined. If both parameters are
 *        defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp          pointer to a cs_equation_param_t structure
 * \param[in]      n_faces      number of faces to enforce
 * \param[in]      face_ids     list of faces
 * \param[in]      ref_value    default values or ignored (may be nullptr)
 * \param[in]      face_values  list of associated values, ignored if nullptr
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_face_dof_enforcement(cs_equation_param_t *eqp,
                                     cs_lnum_t            n_faces,
                                     const cs_lnum_t      face_ids[],
                                     const cs_real_t      ref_value[],
                                     const cs_real_t      face_values[])
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOFB
      && eqp->space_scheme != CS_SPACE_SCHEME_HHO_P0)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: Invalid space scheme.\n"
              "This should be a face-based one.",
              __func__, eqp->name);
  if (ref_value == nullptr && face_values == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: No enforcement value.\n",
              __func__, eqp->name);

  cs_enforcement_param_t *efp            = nullptr;
  int  enforcement_id = eqp->n_enforcements;

  eqp->n_enforcements++;

  if (face_values == nullptr) {

    assert(ref_value != nullptr);
    efp = cs_enforcement_param_create(CS_ENFORCEMENT_SELECTION_FACES,
                                      CS_ENFORCEMENT_BY_CONSTANT,
                                      eqp->dim,
                                      n_faces,
                                      face_ids,
                                      ref_value);

  }
  else
    efp = cs_enforcement_param_create(CS_ENFORCEMENT_SELECTION_FACES,
                                      CS_ENFORCEMENT_BY_DOF_VALUES,
                                      eqp->dim,
                                      n_faces,
                                      face_ids,
                                      face_values);

  BFT_REALLOC(
    eqp->enforcement_params, eqp->n_enforcements, cs_enforcement_param_t *);
  eqp->enforcement_params[enforcement_id] = efp;
  eqp->flag |= CS_EQUATION_FORCE_VALUES;

  return efp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add an enforcement of the value related to the degrees of freedom
 *        associated to the list of selected cells.
 *
 *        One assumes that values are interlaced if eqp->dim > 1
 *        ref_value or elt_values has to be defined. If both parameters are
 *        defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp          pointer to a cs_equation_param_t structure
 * \param[in]      n_cells      number of selected cells
 * \param[in]      cell_ids     list of cell ids
 * \param[in]      ref_value    ignored if nullptr
 * \param[in]      cell_values  list of associated values, ignored if nullptr
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_cell_enforcement(cs_equation_param_t *eqp,
                                 cs_lnum_t            n_cells,
                                 const cs_lnum_t      cell_ids[],
                                 const cs_real_t      ref_value[],
                                 const cs_real_t      cell_values[])
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (ref_value == nullptr && cell_values == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: No enforcement value.\n", __func__, eqp->name);

  cs_enforcement_param_t *efp            = nullptr;
  int  enforcement_id = eqp->n_enforcements;

  eqp->n_enforcements++;

  if (cell_values == nullptr) {

    assert(ref_value != nullptr);
    efp = cs_enforcement_param_create(CS_ENFORCEMENT_SELECTION_CELLS,
                                      CS_ENFORCEMENT_BY_CONSTANT,
                                      eqp->dim,
                                      n_cells,
                                      cell_ids,
                                      ref_value);
  }
  else
    efp = cs_enforcement_param_create(CS_ENFORCEMENT_SELECTION_CELLS,
                                      CS_ENFORCEMENT_BY_DOF_VALUES,
                                      eqp->dim,
                                      n_cells,
                                      cell_ids,
                                      cell_values);

  BFT_REALLOC(
    eqp->enforcement_params, eqp->n_enforcements, cs_enforcement_param_t *);
  eqp->enforcement_params[enforcement_id] = efp;
  eqp->flag |= CS_EQUATION_FORCE_VALUES;

  return efp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a new enforcement if enforcement_id does not exist or replace it
 *        otherwise. Enforcement of the value related to the degrees of freedom
 *        associated to the list of selected cells.
 *
 *        One assumes that values are interlaced if eqp->dim > 1
 *        ref_value or elt_values has to be defined. If both parameters are
 *        defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp             pointer to a cs_equation_param_t structure
 * \param[in]      enforcement_id  id of the enforcement to handle
 * \param[in]      n_cells         number of selected cells
 * \param[in]      cell_ids        list of cell ids
 * \param[in]      ref_value       ignored if nullptr
 * \param[in]      cell_values     list of associated values, ignored if nullptr
 *
 * \return a pointer to a cs_enforcement_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_enforcement_param_t *
cs_equation_add_or_replace_cell_enforcement(cs_equation_param_t *eqp,
                                            int                  enforcement_id,
                                            cs_lnum_t            n_cells,
                                            const cs_lnum_t      cell_ids[],
                                            const cs_real_t      ref_value[],
                                            const cs_real_t      cell_values[])
{
  if (eqp == nullptr)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (ref_value == nullptr && cell_values == nullptr)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: No enforcement value.\n", __func__, eqp->name);

  cs_enforcement_param_t *efp = nullptr;

  if (enforcement_id > eqp->n_enforcements)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid enforcement id.\n", __func__);

  if (enforcement_id == eqp->n_enforcements) { /* Add a new enforcement */

    efp = cs_equation_add_cell_enforcement(eqp,
                                           n_cells, cell_ids,
                                           ref_value, cell_values);

  }
  else { /* Replace an existing parameter structure */

    assert(enforcement_id < eqp->n_enforcements);
    efp = eqp->enforcement_params[enforcement_id];

    if (cell_values == nullptr) {

      assert(ref_value != nullptr);
      cs_enforcement_param_reset(efp,
                                 CS_ENFORCEMENT_SELECTION_CELLS,
                                 CS_ENFORCEMENT_BY_CONSTANT,
                                 eqp->dim,
                                 n_cells,
                                 cell_ids,
                                 ref_value);
    }
    else
      cs_enforcement_param_reset(efp,
                                 CS_ENFORCEMENT_SELECTION_CELLS,
                                 CS_ENFORCEMENT_BY_DOF_VALUES,
                                 eqp->dim,
                                 n_cells,
                                 cell_ids,
                                 cell_values);

  } /* Reset an existant enforcement parameter */

  return efp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Clear time_control data from the given equation param.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_time_control_clear
(
  cs_equation_param_t *eqp /*!<[in,out] Equation param instance */
)
{
  if (eqp->time_control != nullptr) {
    if (eqp->time_control_owner)
      BFT_FREE(eqp->time_control);
    else
      eqp->time_control = nullptr;
  }
  eqp->time_control_owner = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a time control instance to an equation param structure.
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_time_control_add
(
  cs_equation_param_t *eqp,          /*!<[in] pointer to equation param structure */
  cs_time_control_t   *time_control, /*!<[in] pointer to cs_time_control to add. */
  bool                 shallow_copy  /*!<[in] should the addition be only a
                                              shallow copy (pointer) or a real
                                              copy (data copy). */
)
{
  /* Clear existing time control */
  cs_equation_time_control_clear(eqp);

  /* Copy if needed */
  eqp->time_control_owner = (shallow_copy) ? false : true;
  if (shallow_copy)
    eqp->time_control = time_control;
  else {
    BFT_MALLOC(eqp->time_control, 1, cs_time_control_t);
    cs_time_control_copy(time_control, eqp->time_control);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a time control instance to an equation param structure based
 * on default values (hard copy).
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_time_control_add_default
(
  cs_equation_param_t *eqp /*!<[in] pointer to equation param structure */
)
{
  /* Clear existing time control */
  cs_equation_time_control_clear(eqp);

  /* Hard copy from the default values */
  eqp->time_control_owner = true;
  BFT_MALLOC(eqp->time_control, 1, cs_time_control_t);
  cs_time_control_copy_from_default(eqp->time_control);
}

/*----------------------------------------------------------------------------*/
END_C_DECLS
