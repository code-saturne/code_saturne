/*============================================================================
 * Routines to handle specific settings related to a cs_equation_t structure
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
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_boundary_zone.h"
#include "cs_cdo_bc.h"
#include "cs_fp_exception.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_multigrid.h"
#include "cs_sles.h"
#include "cs_source_term.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

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
 * \brief  Set a parameter attached to a keyname in a \ref cs_equation_param_t
 *         structure
 *
 * \param[in, out]  eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in]       key      key related to the member of eq to set
 * \param[in]       keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

static void
_set_key(cs_equation_param_t   *eqp,
         cs_equation_key_t      key,
         const char            *keyval)
{
  const char  *eqname = eqp->name;
  const char  emsg[] = " %s: %s equation --> Invalid key value %s for"
    " keyword %s.\n";

  switch(key) {

  case CS_EQKEY_ADV_EXTRAPOL:
    if (strcmp(keyval, "none") == 0)
      eqp->adv_extrapol = CS_PARAM_ADVECTION_EXTRAPOL_NONE;
    else if (strcmp(keyval, "taylor") == 0)
      eqp->adv_formulation = CS_PARAM_ADVECTION_EXTRAPOL_TAYLOR_2;
    else if (strcmp(keyval, "adams_bashforth") == 0)
      eqp->adv_formulation = CS_PARAM_ADVECTION_EXTRAPOL_ADAMS_BASHFORTH_2;
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
    /* Automatic witch to a hybrid upwind/centered scheme for advection */
    eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND;
    break;

  case CS_EQKEY_AMG_TYPE:
    if (strcmp(keyval, "none") == 0 || strcmp(keyval, "") == 0)
      eqp->sles_param->amg_type = CS_PARAM_AMG_NONE;
    else if (strcmp(keyval, "v_cycle") == 0) {
      eqp->sles_param->amg_type = CS_PARAM_AMG_HOUSE_V;
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_CS;
    }
    else if (strcmp(keyval, "k_cycle") == 0) {
      eqp->sles_param->amg_type = CS_PARAM_AMG_HOUSE_K;
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_CS;
    }
    else if (strcmp(keyval, "boomer") == 0) {
      eqp->sles_param->amg_type = CS_PARAM_AMG_HYPRE_BOOMER;
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_HYPRE;
    }
    else if (strcmp(keyval, "gamg") == 0) {
      eqp->sles_param->amg_type = CS_PARAM_AMG_PETSC_GAMG;
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_PETSC;
    }
    else if (strcmp(keyval, "pcmg") == 0) {
      eqp->sles_param->amg_type = CS_PARAM_AMG_PETSC_PCMG;
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_PETSC;
    }
    else {
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

  case CS_EQKEY_BC_QUADRATURE:
    {
      cs_quadrature_type_t  qtype = CS_QUADRATURE_NONE;

      if (strcmp(keyval, "bary") == 0)
        qtype = CS_QUADRATURE_BARY;
      else if (strcmp(keyval, "bary_subdiv") == 0)
        qtype = CS_QUADRATURE_BARY_SUBDIV;
      else if (strcmp(keyval, "higher") == 0)
        qtype = CS_QUADRATURE_HIGHER;
      else if (strcmp(keyval, "highest") == 0)
        qtype = CS_QUADRATURE_HIGHEST;
      else {
        const char *_val = keyval;
        bft_error(__FILE__, __LINE__, 0,
                  emsg, __func__, eqname, _val, "CS_EQKEY_BC_QUADRATURE");
      }

      for (int i = 0; i < eqp->n_bc_defs; i++)
        cs_xdef_set_quadrature(eqp->bc_defs[i], qtype);

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
    if (strcmp(keyval, "derham") == 0)
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
      eqp->process_flag |= CS_EQUATION_POST_BALANCE;
    else if (strcmp(keyval, "peclet") == 0)
      eqp->process_flag |= CS_EQUATION_POST_PECLET;
    else if (strcmp(keyval, "upwind_coef") == 0)
      eqp->process_flag |= CS_EQUATION_POST_UPWIND_COEF;
    else if (strcmp(keyval, "normal_flux") == 0)
      eqp->process_flag |= CS_EQUATION_POST_NORMAL_FLUX;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_EXTRA_OP");
    }
    break;

  case CS_EQKEY_HODGE_DIFF_ALGO:
    if (strcmp(keyval,"cost") == 0 || strcmp(keyval,"ocs") == 0)
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_COST;
    else if (strcmp(keyval, "ocs2") == 0)
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_OCS2;
    else if (strcmp(keyval, "bubble") == 0)
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_BUBBLE;
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
    if (strcmp(keyval, "dga") == 0)
      eqp->diffusion_hodgep.coef = 1./3.;
    else if (strcmp(keyval, "sushi") == 0)
      eqp->diffusion_hodgep.coef = 1./sqrt(3.);
    else if (strcmp(keyval, "gcr") == 0)
      eqp->diffusion_hodgep.coef = 1.0;
    else if (strcmp(keyval, "frac23") == 0 || strcmp(keyval, "2/3") == 0)
      eqp->diffusion_hodgep.coef = 2./3.;
    else
      eqp->diffusion_hodgep.coef = atof(keyval);
    break;

  case CS_EQKEY_HODGE_TIME_ALGO:
    if (strcmp(keyval, "voronoi") == 0)
      eqp->time_hodgep.algo = CS_HODGE_ALGO_VORONOI;
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
    else if (strcmp(keyval, "wbs") == 0)
      eqp->reaction_hodgep.algo = CS_HODGE_ALGO_WBS;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_HODGE_REAC_ALGO");
    }
    break;

  case CS_EQKEY_ITSOL:
    if (strcmp(keyval, "amg") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_AMG;
    else if (strcmp(keyval, "bicg") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_BICG;
    else if (strcmp(keyval, "bicgstab2") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_BICGSTAB2;
    else if (strcmp(keyval, "cg") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_CG;
    else if (strcmp(keyval, "cr3") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_CR3;
    else if (strcmp(keyval, "fcg") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_FCG;
    else if (strcmp(keyval, "gauss_seidel") == 0 ||
             strcmp(keyval, "gs") == 0) {
      eqp->sles_param->solver = CS_PARAM_ITSOL_GAUSS_SEIDEL;
      eqp->sles_param->precond = CS_PARAM_PRECOND_NONE;
    }
    else if (strcmp(keyval, "gcr") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_GCR;
    else if (strcmp(keyval, "gmres") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_GMRES;
    else if (strcmp(keyval, "fgmres") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_FGMRES;
    else if (strcmp(keyval, "jacobi") == 0 || strcmp(keyval, "diag") == 0 ||
             strcmp(keyval, "diagonal") == 0) {
      eqp->sles_param->solver = CS_PARAM_ITSOL_JACOBI;
      eqp->sles_param->precond = CS_PARAM_PRECOND_NONE;
    }
    else if (strcmp(keyval, "minres") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_MINRES;

    else if (strcmp(keyval, "mumps") == 0 ||
             strcmp(keyval, "mumps_float") == 0 ||
             strcmp(keyval, "mumps_float_ldlt") == 0 ||
             strcmp(keyval, "mumps_ldlt") == 0) {

      eqp->sles_param->precond = CS_PARAM_PRECOND_NONE;

      /* Modify the default and check availability of MUMPS solvers */
      if (eqp->sles_param->solver_class == CS_PARAM_SLES_CLASS_CS ||
          eqp->sles_param->solver_class == CS_PARAM_SLES_CLASS_MUMPS) {
#if defined(HAVE_MUMPS)
        eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_MUMPS;
#else
#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_MUMPS)
        eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_PETSC;
#else  /* MUMPS: no, PETSc: yes: PETSc with MUMPS: no */
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Error detected while setting \"%s\" key for eq. %s\n"
                  " MUMPS is not available with your installation.\n"
                  " Please check your installation settings.\n",
                  __func__, "CS_EQKEY_ITSOL", eqname);
#endif  /* PETSC_HAVE_MUMPS */
#else   /* MUMPS: no, PETSc: no */
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Error detected while setting \"%s\" key for eq. %s\n"
                  " MUMPS is not available with your installation.\n"
                  " Please check your installation settings.\n",
                  __func__, "CS_EQKEY_ITSOL", eqname);
#endif  /* HAVE_PETSC */
#endif  /* HAVE_MUMPS */
      }

      /* Sanity check */
      assert(eqp->sles_param->solver_class != CS_PARAM_SLES_CLASS_CS &&
             eqp->sles_param->solver_class != CS_PARAM_SLES_CLASS_HYPRE);

      if (strcmp(keyval, "mumps") == 0)
        eqp->sles_param->solver = CS_PARAM_ITSOL_MUMPS;
      else if (strcmp(keyval, "mumps_float") == 0)
        eqp->sles_param->solver = CS_PARAM_ITSOL_MUMPS_FLOAT;
      else if (strcmp(keyval, "mumps_float_ldlt") == 0)
        eqp->sles_param->solver = CS_PARAM_ITSOL_MUMPS_FLOAT_LDLT;
      else if (strcmp(keyval, "mumps_ldlt") == 0)
        eqp->sles_param->solver = CS_PARAM_ITSOL_MUMPS_LDLT;

    }
    else if (strcmp(keyval, "sym_gauss_seidel") == 0 ||
             strcmp(keyval, "sgs") == 0) {
      eqp->sles_param->solver = CS_PARAM_ITSOL_SYM_GAUSS_SEIDEL;
      eqp->sles_param->precond = CS_PARAM_PRECOND_NONE;
    }
    else if (strcmp(keyval, "user") == 0) {
      eqp->sles_param->solver = CS_PARAM_ITSOL_USER_DEFINED;
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_CS;
    }
    else if (strcmp(keyval, "none") == 0)
      eqp->sles_param->solver = CS_PARAM_ITSOL_NONE;

    else { /* keyval not found among the available keyvals */
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_ITSOL");
    }
    break;

  case CS_EQKEY_ITSOL_MAX_ITER:
    eqp->sles_param->n_max_iter = atoi(keyval);
    break;

  case CS_EQKEY_ITSOL_EPS:
    eqp->sles_param->eps = atof(keyval);
    break;

  case CS_EQKEY_ITSOL_RESNORM_TYPE:
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
                emsg, __func__, eqname, _val, "CS_EQKEY_ITSOL_RESNORM_TYPE");
    }
    break;

  case CS_EQKEY_ITSOL_RESTART:
    eqp->sles_param->restart = atoi(keyval);
    break;

  case CS_EQKEY_OMP_ASSEMBLY_STRATEGY:
    if (strcmp(keyval, "critical") == 0)
      eqp->omp_assembly_choice = CS_PARAM_ASSEMBLE_OMP_CRITICAL;
    else if (strcmp(keyval, "atomic") == 0)
      eqp->omp_assembly_choice = CS_PARAM_ASSEMBLE_OMP_ATOMIC;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_OMP_ASSEMBLY_STRATEGY");
    }
    break;

  case CS_EQKEY_PRECOND:
    if (strcmp(keyval, "none") == 0) {
      eqp->sles_param->precond = CS_PARAM_PRECOND_NONE;
      eqp->sles_param->amg_type = CS_PARAM_AMG_NONE;
    }
    else if (strcmp(keyval, "jacobi") == 0)
      eqp->sles_param->precond = CS_PARAM_PRECOND_DIAG;
    else if (strcmp(keyval, "block_jacobi") == 0 ||
             strcmp(keyval, "block_jacobi_ilu0") == 0 ||
             strcmp(keyval, "jacobi_block") == 0) {
      eqp->sles_param->precond = CS_PARAM_PRECOND_BJACOB_ILU0;
      /* Default when using PETSc */
      eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_NORM2_RHS;
    }
    else if (strcmp(keyval, "block_jacobi_sgs") == 0 ||
             strcmp(keyval, "block_jacobi_ssor") == 0) {
      eqp->sles_param->precond = CS_PARAM_PRECOND_BJACOB_SGS;
      /* Default when using PETSc */
      eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_NORM2_RHS;
    }
    else if (strcmp(keyval, "poly1") == 0)
      eqp->sles_param->precond = CS_PARAM_PRECOND_POLY1;
    else if (strcmp(keyval, "poly2") == 0)
      eqp->sles_param->precond = CS_PARAM_PRECOND_POLY2;
    else if (strcmp(keyval, "ssor") == 0)
      eqp->sles_param->precond = CS_PARAM_PRECOND_SSOR;
    else if (strcmp(keyval, "ilu0") == 0)
      eqp->sles_param->precond = CS_PARAM_PRECOND_ILU0;
    else if (strcmp(keyval, "icc0") == 0)
      eqp->sles_param->precond = CS_PARAM_PRECOND_ICC0;
    else if (strcmp(keyval, "amg") == 0) {

      eqp->sles_param->precond = CS_PARAM_PRECOND_AMG;

      /* Set the default choice */
      if (eqp->sles_param->solver_class == CS_PARAM_SLES_CLASS_CS)
        eqp->sles_param->amg_type = CS_PARAM_AMG_HOUSE_K;
      else if (eqp->sles_param->solver_class == CS_PARAM_SLES_CLASS_PETSC) {
        eqp->sles_param->amg_type = CS_PARAM_AMG_PETSC_GAMG;
        /* Default when using PETSc */
        eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_NORM2_RHS;
      }
#if defined(PETSC_HAVE_HYPRE)
      /* Up to now HYPRE is available only through the PETSc interface */
      else if (eqp->sles_param->solver_class == CS_PARAM_SLES_CLASS_HYPRE) {
        eqp->sles_param->amg_type = CS_PARAM_AMG_HYPRE_BOOMER;
        /* Default when using PETSc */
        eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_NORM2_RHS;
      }
#endif
      else
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid choice of AMG type for equation %s.\n"
                  " Please modify your settings",
                  __func__, eqname);

    }
    else if (strcmp(keyval, "amg_block") == 0 ||
             strcmp(keyval, "block_amg") == 0) {

      eqp->sles_param->precond = CS_PARAM_PRECOND_AMG;
      eqp->sles_param->pcd_block_type = CS_PARAM_PRECOND_BLOCK_DIAG;

      /* Set the default choice */
      switch (eqp->sles_param->solver_class) {

      case CS_PARAM_SLES_CLASS_CS:
        eqp->sles_param->amg_type = CS_PARAM_AMG_HOUSE_K;
        break;
      case CS_PARAM_SLES_CLASS_PETSC:
        eqp->sles_param->amg_type = CS_PARAM_AMG_PETSC_GAMG;
        /* Default when using PETSc */
        eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_NORM2_RHS;
        break;
      case CS_PARAM_SLES_CLASS_HYPRE:
        /* Set the default choice */
#if defined(PETSC_HAVE_HYPRE)
        eqp->sles_param->amg_type = CS_PARAM_AMG_HYPRE_BOOMER;
        eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_HYPRE;
#else
        cs_base_warn(__FILE__, __LINE__);
        bft_printf(" Switch to PETSc multigrid since Hypre is not available"
                   " for equation %s.\n", eqname);
        eqp->sles_param->amg_type = CS_PARAM_AMG_PETSC_GAMG;
        eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_PETSC;
#endif
        /* Default when using PETSc */
        eqp->sles_param->resnorm_type = CS_PARAM_RESNORM_NORM2_RHS;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid choice for block AMG for equation %s\n"
                  " Please check your settings.",
                  __func__, eqname);
        break;
      }

    }
    else if (strcmp(keyval, "as") == 0)
      eqp->sles_param->precond = CS_PARAM_PRECOND_AS;
    else {
      const char *_val = keyval;
      bft_error(__FILE__, __LINE__, 0,
                emsg, __func__, eqname, _val, "CS_EQKEY_PRECOND");
    }
    break;

  case CS_EQKEY_SLES_VERBOSITY: /* "verbosity" for SLES structures */
    eqp->sles_param->verbosity = atoi(keyval);
    break;

  case CS_EQKEY_SOLVER_FAMILY:
    if (strcmp(keyval, "cs") == 0)
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_CS;
    else if (strcmp(keyval, "hypre") == 0)
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_HYPRE;
    else if (strcmp(keyval, "mumps") == 0)
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_MUMPS;
    else if (strcmp(keyval, "petsc") == 0)
      eqp->sles_param->solver_class = CS_PARAM_SLES_CLASS_PETSC;
    else {
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
      eqp->time_hodgep.type = CS_HODGE_TYPE_VPCD;
      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_EPFD;
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_BUBBLE;
      eqp->diffusion_hodgep.coef = 2*cs_math_1ov3;
      eqp->reaction_hodgep.type = CS_HODGE_TYPE_VPCD;
      eqp->reaction_hodgep.algo = CS_HODGE_ALGO_WBS;
    }
    else if (strcmp(keyval, "cdo_vcb") == 0 ||
             strcmp(keyval, "cdovcb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOVCB;
      eqp->space_poly_degree = 0;
      eqp->time_hodgep.type = CS_HODGE_TYPE_VPCD;
      eqp->diffusion_hodgep.algo = CS_HODGE_ALGO_WBS;
      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_VC;
      eqp->reaction_hodgep.type = CS_HODGE_TYPE_VPCD;
      eqp->reaction_hodgep.algo = CS_HODGE_ALGO_WBS;
    }
    else if (strcmp(keyval, "cdo_fb") == 0 ||
             strcmp(keyval, "cdofb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOFB;
      eqp->space_poly_degree = 0;
      eqp->time_hodgep.type = CS_HODGE_TYPE_CPVD;
      eqp->time_hodgep.algo = CS_HODGE_ALGO_VORONOI;
      eqp->reaction_hodgep.algo = CS_HODGE_ALGO_VORONOI;
      eqp->diffusion_hodgep.type = CS_HODGE_TYPE_EDFP;
    }
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
 * \brief  Create a \ref cs_equation_param_t structure with the given
 *         parameters. The remaining parameters are set with default values;
 *
 * \param[in] name          name of the equation
 * \param[in] type          type of equation
 * \param[in] dim           dim of the variable associated to this equation
 * \param[in] default_bc    type of boundary condition set by default
 *
 * \return a pointer to a new allocated \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_create_param(const char            *name,
                         cs_equation_type_t     type,
                         int                    dim,
                         cs_param_bc_type_t     default_bc)
{
  cs_equation_param_t  *eqp = NULL;

  BFT_MALLOC(eqp, 1, cs_equation_param_t);

  /* Store the name of the equation */
  int  len = strlen(name)+1;
  BFT_MALLOC(eqp->name, len, char);
  strncpy(eqp->name, name, len);

  /* Set additional members */
  eqp->type = type;
  eqp->dim = dim;

  /* Other default settings */
  eqp->verbosity = 2;
  eqp->flag = 0;
  eqp->process_flag = 0;

  /* Vertex-based schemes imply specific discrete Hodge operators for
     diffusion, time and reaction terms.
     Default initialization is made in accordance with this choice */
  eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;
  eqp->dof_reduction = CS_PARAM_REDUCTION_DERHAM;
  eqp->space_poly_degree = 0;

  /* Default initialization for the legacy var_col_opt structure which is now
   * shared inside the cs_equation_param_t structure
   * The default value used here should be the same as the one considered in
   * src/base/cs_parameters.c (static variable called _equation_param_default) */

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
  eqp->thetav = 1.;
  eqp->blencv = 1.;
  eqp->blend_st = 0.;
  eqp->epsilo = 1.e-5;
  eqp->epsrsm = 1.e-4;
  eqp->epsrgr = 1.e-4;
  eqp->climgr = 1.5;
  eqp->extrag = 0.;
  eqp->relaxv = 1.;

  /* Boundary conditions structure.
     One assigns a boundary condition by default */
  eqp->default_bc = default_bc;
  eqp->n_bc_defs = 0;
  eqp->bc_defs = NULL;
  eqp->default_enforcement = CS_PARAM_BC_ENFORCE_ALGEBRAIC;
  eqp->strong_pena_bc_coeff = _strong_pena_bc_coef_by_default;
  eqp->weak_pena_bc_coeff = _weak_pena_bc_coef_by_default;

  /* Initial condition (zero value by default) */
  eqp->n_ic_defs = 0;
  eqp->ic_defs = NULL;

  /* Description of the time discretization (default values) */
  eqp->time_property = NULL;
  eqp->time_scheme = CS_TIME_SCHEME_EULER_IMPLICIT;
  eqp->theta = 1.0;
  eqp->do_lumping = false;
  eqp->time_hodgep = (cs_hodge_param_t) {
    .inv_pty = false,
    .algo = CS_HODGE_ALGO_VORONOI,
    .type = CS_HODGE_TYPE_VPCD,
    .coef = 1.,
  };

  /* Description of the discetization of the diffusion term */
  eqp->diffusion_property = NULL;
  eqp->diffusion_hodgep = (cs_hodge_param_t) {
    .inv_pty = false,
    .algo = CS_HODGE_ALGO_COST,
    .type = CS_HODGE_TYPE_EPFD,
    .coef = 1./3.,
  };

  /* Description of the discetization of the curl-curl term */
  eqp->curlcurl_property = NULL;
  eqp->curlcurl_hodgep = (cs_hodge_param_t) {
    .inv_pty = true,
    .algo = CS_HODGE_ALGO_COST,
    .type = CS_HODGE_TYPE_FPED,
    .coef = 1./3.,
  };

  /* Description of the discetization of the grad-div term */
  eqp->graddiv_property = NULL;
  eqp->graddiv_hodgep = (cs_hodge_param_t) {
    .inv_pty = false,
    .algo = CS_HODGE_ALGO_VORONOI,
    .type = CS_HODGE_TYPE_EPFD,
    .coef = 1./3.,
  };

  /* Advection term */
  eqp->adv_field = NULL;
  eqp->adv_scaling_property = NULL;
  eqp->adv_extrapol = CS_PARAM_ADVECTION_EXTRAPOL_NONE;
  eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_CONSERV;
  eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
  eqp->adv_strategy = CS_PARAM_ADVECTION_IMPLICIT_FULL;
  eqp->upwind_portion = 0.15;

  /* Description of the discretization of the reaction term.
     No reaction term by default */
  eqp->n_reaction_terms = 0;
  eqp->reaction_properties = NULL;
  eqp->reaction_hodgep = (cs_hodge_param_t) {
    .inv_pty = false,
    .algo = CS_HODGE_ALGO_WBS,
    .type = CS_HODGE_TYPE_VPCD,
  };

  /* Source term (always in the right-hand side)
     No source term by default */
  eqp->n_source_terms = 0;
  eqp->source_terms = NULL;

  /* Mass injection in the volume term (always in the right-hand side)
     No volume mass injection term by default */
  eqp->n_volume_mass_injections = 0;
  eqp->volume_mass_injections = NULL;

  /* Members of the structure to handle the enforcement of (internal) DoFs */
  eqp->enforcement_type = 0;
  BFT_MALLOC(eqp->enforcement_ref_value, eqp->dim, cs_real_t);
  for (int i = 0; i < eqp->dim; i++)
    eqp->enforcement_ref_value[i] = 0.;

  eqp->n_enforced_cells = 0;
  eqp->enforced_cell_ids = NULL;
  eqp->enforced_cell_values = NULL;

  eqp->n_enforced_dofs = 0;
  eqp->enforced_dof_ids = NULL;
  eqp->enforced_dof_values = NULL;

  /* Settings for driving the linear algebra */
  eqp->sles_param = cs_param_sles_create(-1, name); /* field_id, system_name */

  /* Settings for the OpenMP strategy */
  eqp->omp_assembly_choice = CS_PARAM_ASSEMBLE_OMP_CRITICAL;

  return eqp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy the settings from one \ref cs_equation_param_t structure to
 *         another one. The name is not copied.
 *
 * \param[in]      ref       pointer to the reference \ref cs_equation_param_t
 * \param[in, out] dst       pointer to the \ref cs_equation_param_t to update
 * \param[in]      copy_fid  copy also the field id or not
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_copy_param_from(const cs_equation_param_t   *ref,
                            cs_equation_param_t         *dst,
                            bool                         copy_fid)
{
  /* Generic members */
  dst->type = ref->type;
  dst->dim = ref->dim;
  dst->verbosity = ref->verbosity;
  dst->process_flag = ref->process_flag;
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
  dst->thetav = ref->thetav;
  dst->blencv = ref->blencv;
  dst->blend_st = ref->blend_st;
  dst->epsilo = ref->epsilo;
  dst->epsrsm = ref->epsrsm;
  dst->epsrgr = ref->epsrgr;
  dst->climgr = ref->climgr;
  dst->extrag = ref->extrag;
  dst->relaxv = ref->relaxv;

  /* Boundary conditions structure */
  dst->default_bc = ref->default_bc;
  dst->default_enforcement = ref->default_enforcement;
  dst->strong_pena_bc_coeff = ref->strong_pena_bc_coeff;
  dst->n_bc_defs = ref->n_bc_defs;
  BFT_MALLOC(dst->bc_defs, dst->n_bc_defs, cs_xdef_t *);
  for (int i = 0; i < ref->n_bc_defs; i++)
    dst->bc_defs[i] = cs_xdef_copy(ref->bc_defs[i]);

  /* Description of the time discretization */
  dst->time_scheme = ref->time_scheme;
  dst->theta = ref->theta;
  dst->do_lumping = ref->do_lumping;
  dst->time_property = ref->time_property;

  cs_hodge_copy_parameters(&(ref->time_hodgep), &(dst->time_hodgep));

  /* Initial condition (zero value by default) */
  dst->n_ic_defs = ref->n_ic_defs;
  BFT_MALLOC(dst->ic_defs, dst->n_ic_defs, cs_xdef_t *);
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
  dst->adv_extrapol = ref->adv_extrapol;
  dst->adv_formulation = ref->adv_formulation;
  dst->adv_scheme = ref->adv_scheme;
  dst->adv_strategy = ref->adv_strategy;
  dst->upwind_portion = ref->upwind_portion;
  dst->adv_field = ref->adv_field;
  dst->adv_scaling_property = ref->adv_scaling_property;

  /* Reaction term */
  dst->n_reaction_terms = ref->n_reaction_terms;
  BFT_MALLOC(dst->reaction_properties, dst->n_reaction_terms, cs_property_t *);
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
  BFT_MALLOC(dst->volume_mass_injections,
             dst->n_volume_mass_injections,
             cs_xdef_t *);
  for (int i = 0; i < dst->n_volume_mass_injections; i++)
    dst->volume_mass_injections[i]
      = cs_xdef_copy(ref->volume_mass_injections[i]);

  /* No enforcement of internal DoFs */
  dst->enforcement_type = ref->enforcement_type;
  BFT_MALLOC(dst->enforcement_ref_value, dst->dim, cs_real_t);
  memcpy(dst->enforcement_ref_value, ref->enforcement_ref_value,
         dst->dim*sizeof(cs_real_t));

  dst->n_enforced_cells = ref->n_enforced_cells;
  if (ref->n_enforced_cells > 0) {

    BFT_MALLOC(dst->enforced_cell_ids, dst->n_enforced_cells, cs_lnum_t);
    memcpy(dst->enforced_cell_ids, ref->enforced_cell_ids,
           dst->n_enforced_cells*sizeof(cs_lnum_t));

    dst->enforced_cell_values= NULL;
    if (ref->enforced_cell_values != NULL) {
      cs_lnum_t  size = dst->n_enforced_cells*dst->dim;
      BFT_MALLOC(dst->enforced_cell_values, size, cs_real_t);
      memcpy(dst->enforced_dof_values, ref->enforced_dof_values,
             size*sizeof(cs_real_t));
    }

  }

  dst->n_enforced_dofs = ref->n_enforced_dofs;
  if (ref->n_enforced_dofs > 0) {

    BFT_MALLOC(dst->enforced_dof_ids, dst->n_enforced_dofs, cs_lnum_t);
    memcpy(dst->enforced_dof_ids, ref->enforced_dof_ids,
           dst->n_enforced_dofs*sizeof(cs_lnum_t));

    dst->enforced_dof_values = NULL;
    if (ref->enforced_dof_values != NULL) {
      cs_lnum_t  size = dst->n_enforced_dofs*dst->dim;
      BFT_MALLOC(dst->enforced_dof_values, size, cs_real_t);
      memcpy(dst->enforced_dof_values, ref->enforced_dof_values,
             size*sizeof(cs_real_t));
    }

  }

  /* Copy the settings driving the linear algebra algorithms */
  if (copy_fid)
    cs_param_sles_copy_from(ref->sles_param, dst->sles_param);

  else {

    int save_field_id = dst->sles_param->field_id;
    cs_param_sles_copy_from(ref->sles_param, dst->sles_param);
    dst->sles_param->field_id = save_field_id;

  }

  /* Settings related to the performance */
  dst->omp_assembly_choice = ref->omp_assembly_choice;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free the contents of a \ref cs_equation_param_t
 *
 * The cs_equation_param_t structure itself is not freed, but all its
 * sub-structures are freed.
 *
 * This is useful for equation parameters which are accessed through
 * field keywords.
 *
 * \param[in, out]  eqp  pointer to a \ref cs_equation_param_t
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_clear_param(cs_equation_param_t   *eqp)
{
  if (eqp == NULL)
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
  BFT_FREE(eqp->enforcement_ref_value);

  if (eqp->n_enforced_cells > 0) {
    eqp->n_enforced_cells = 0;
    BFT_FREE(eqp->enforced_cell_ids);
    BFT_FREE(eqp->enforced_cell_values);
  }

  if (eqp->n_enforced_dofs > 0) {
    eqp->n_enforced_dofs = 0;
    BFT_FREE(eqp->enforced_dof_ids);
    BFT_FREE(eqp->enforced_dof_values);
  }

  /* Information related to the definition of initial conditions */
  if (eqp->n_ic_defs > 0) {

    for (int i = 0; i < eqp->n_ic_defs; i++)
      eqp->ic_defs[i] = cs_xdef_free(eqp->ic_defs[i]);
    BFT_FREE(eqp->ic_defs);

  }

  /* Information related to the linear algebra settings */
  cs_param_sles_free(&(eqp->sles_param));

  BFT_FREE(eqp->name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_equation_param_t
 *
 * \param[in, out] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_free_param(cs_equation_param_t     *eqp)
{
  if (eqp == NULL)
    return NULL;

  cs_equation_clear_param(eqp);

  BFT_FREE(eqp);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname in a \ref cs_equation_param_t
 *         structure
 *
 * \param[in, out]  eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in]       key      key related to the member of eq to set
 * \param[in]       keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_param(cs_equation_param_t   *eqp,
                      cs_equation_key_t      key,
                      const char            *keyval)
{
  /* Sanity checks */
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (keyval == NULL)
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
 * \brief Set parameters for initializing SLES structures used for the
 *        resolution of the linear system.
 *        Settings are related to this equation.
 *
 * \param[in, out]  eqp           pointer to a cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_sles(cs_equation_param_t      *eqp)
{
  /* Define a cs_sles_t structure using the field_id related to the variable
   * field associated to this equation */
  int  ierr = cs_param_sles_set(true, eqp->sles_param);

  if (ierr == -1)
    bft_error(__FILE__, __LINE__, 0,
              "%s: The requested class of solvers is not available"
              " for the equation %s\n"
              " Please modify your settings.", __func__, eqp->name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last modification of the cs_equation_param_t structure before
 *         launching the computation
 *
 * \param[in, out]  eqp      pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_last_stage(cs_equation_param_t   *eqp)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  if (eqp->flag & CS_EQUATION_LOCKED)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Equation %s is not modifiable anymore.\n"
                " Please check your settings."), eqp->name, __func__);

  if (eqp->do_lumping) { /* Activate a set of options */

    eqp->reaction_hodgep.algo = CS_HODGE_ALGO_VORONOI;
    eqp->time_hodgep.algo = CS_HODGE_ALGO_VORONOI;

    /* A simple barycentric quadrature rule is applied to all source terms */
    for (int i = 0; i < eqp->n_source_terms; i++)
      cs_xdef_set_quadrature(eqp->source_terms[i], CS_QUADRATURE_BARY);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the detail of a \ref cs_equation_param_t structure
 *
 * \param[in]  eqp      pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_log(const cs_equation_param_t   *eqp)
{
  if (eqp == NULL)
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

  if (cs_glob_n_threads > 1) {
    if (eqp->omp_assembly_choice == CS_PARAM_ASSEMBLE_OMP_CRITICAL)
      cs_log_printf(CS_LOG_SETUP, "  * %s | OpenMP.Assembly.Choice:  %s\n",
                    eqname, "critical");
    else if (eqp->omp_assembly_choice == CS_PARAM_ASSEMBLE_OMP_ATOMIC)
      cs_log_printf(CS_LOG_SETUP, "  * %s | OpenMP.Assembly.Choice:  %s\n",
                    eqname, "atomic");
  }

  /* Boundary conditions */
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
                  "  * %s | Boundary conditions | Penalization coefficient:"
                  " %5.3e\n", eqname, eqp->strong_pena_bc_coeff);
  else if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
           eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Boundary conditions | Penalization coefficient:"
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
      cs_xdef_log(prefix, d);
    }
  }

  if (unsteady) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Time settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP,
                  "  * %s | Initial conditions | Number of definitions: %d",
                  eqname, eqp->n_ic_defs);
    if (eqp->n_ic_defs > 0)
      cs_log_printf(CS_LOG_SETUP, "\n\n");
    for (int i = 0; i < eqp->n_ic_defs; i++) {
      sprintf(prefix, "        Definition %3d", i);
      cs_xdef_log(prefix, eqp->ic_defs[i]);
    }

    const char  *time_scheme = cs_param_get_time_scheme_name(eqp->time_scheme);
    if (time_scheme != NULL) {
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

  } /* Unsteady term */

  if (diffusion) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Diffusion term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Diffusion property: %s\n\n",
                  eqname, cs_property_get_name(eqp->diffusion_property));

    sprintf(prefix, "        Diffusion Hodge op. ");
    cs_hodge_param_log(prefix, eqp->diffusion_property, eqp->diffusion_hodgep);

  } /* Diffusion term */

  if (curlcurl) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Curl-Curl term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Curl-Curl property: %s\n\n",
                  eqname, cs_property_get_name(eqp->curlcurl_property));

    sprintf(prefix, "        Curl-curl Hodge op. ");
    cs_hodge_param_log(prefix, eqp->curlcurl_property, eqp->curlcurl_hodgep);

  } /* Curl-curl term */

  if (graddiv) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Grad-Div term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Grad-Div property: %s\n\n",
                  eqname, cs_property_get_name(eqp->graddiv_property));

    sprintf(prefix, "        Grad-Div Hodge op. ");
    cs_hodge_param_log(prefix, eqp->graddiv_property, eqp->graddiv_hodgep);

  } /* Curl-curl term */

  if (convection) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Advection term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Field: %s\n",
                  eqname, cs_advection_field_get_name(eqp->adv_field));
    if (eqp->adv_scaling_property != NULL)
      cs_log_printf(CS_LOG_SETUP, "  * %s | Scaling.Property: %s\n",
                    eqname, cs_property_get_name(eqp->adv_scaling_property));

    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Formulation: %s\n",
                  eqname,
                  cs_param_get_advection_form_name(eqp->adv_formulation));

    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Scheme: %s\n",
                  eqname,
                  cs_param_get_advection_scheme_name(eqp->adv_scheme));
    if (eqp->adv_scheme == CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND)
      cs_log_printf(CS_LOG_SETUP, "  * %s | Upwind.Portion: %3.2f %%\n",
                    eqname, 100*eqp->upwind_portion);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Strategy: %s\n",
                  eqname,
                  cs_param_get_advection_strategy_name(eqp->adv_strategy));
    cs_log_printf(CS_LOG_SETUP, "  * %s | Advection.Extrapolation: %s\n",
                  eqname,
                  cs_param_get_advection_extrapol_name(eqp->adv_extrapol));

  } /* Advection term */

  if (reaction) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Reaction settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Reaction | Number of terms: %d\n",
                  eqname, eqp->n_reaction_terms);

    sprintf(prefix, "        Reaction Hodge op. ");
    cs_hodge_param_log(prefix, NULL, eqp->reaction_hodgep);

  } /* Reaction terms */

  if (source_term) {

    cs_log_printf(CS_LOG_SETUP, "\n### %s | Source term settings\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  * %s | Source terms | Number of terms: %d\n",
                  eqname, eqp->n_source_terms);

    for (int s_id = 0; s_id < eqp->n_source_terms; s_id++) {
      sprintf(prefix, "        Definition %3d", s_id);
      cs_xdef_log(prefix, eqp->source_terms[s_id]);
    }

  } /* Source terms */

  /* Iterative solver information */
  cs_param_sles_log(eqp->sles_param);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Ask if the parameter settings of the equation has requested the
 *         treatment of Robin boundary conditions
 *
 * \param[in] eqp          pointer to a \ref cs_equation_param_t
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_param_has_robin_bc(const cs_equation_param_t     *eqp)
{
  if (eqp == NULL)
    return false;

  for (int i = 0; i < eqp->n_bc_defs; i++) {
    cs_xdef_t  *def = eqp->bc_defs[i];
    if (def->meta & CS_CDO_BC_ROBIN)
      return true;
  }
  return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here a constant value is set to all the entities belonging to the
 *         given mesh location
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       z_name    name of the associated zone (if NULL or
 *                            "" all cells are considered)
 * \param[in]       val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_value(cs_equation_param_t    *eqp,
                            const char             *z_name,
                            cs_real_t              *val)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int  z_id = cs_get_vol_zone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        eqp->dim,
                                        z_id,
                                        CS_FLAG_STATE_UNIFORM, /* state flag */
                                        meta_flag,
                                        val);

  int  new_id = eqp->n_ic_defs;
  eqp->n_ic_defs += 1;
  BFT_REALLOC(eqp->ic_defs, eqp->n_ic_defs, cs_xdef_t *);
  eqp->ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this equation
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the value related to all the entities belonging to the
 *         given mesh location is such that the integral over these cells
 *         returns the requested quantity
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       z_name    name of the associated zone (if NULL or
 *                            "" all cells are considered)
 * \param[in]       quantity  quantity to distribute over the mesh location
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_qov(cs_equation_param_t    *eqp,
                          const char             *z_name,
                          double                  quantity)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = cs_get_vol_zone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_QOV,
                                        eqp->dim,
                                        z_id,
                                        0, /* state flag */
                                        meta_flag,
                                        &quantity);

  int  new_id = eqp->n_ic_defs;
  eqp->n_ic_defs += 1;
  BFT_REALLOC(eqp->ic_defs, eqp->n_ic_defs, cs_xdef_t *);
  eqp->ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the unknown related to this
 *         equation. This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_ic_by_analytic(cs_equation_param_t    *eqp,
                               const char             *z_name,
                               cs_analytic_func_t     *analytic,
                               void                   *input)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = cs_get_vol_zone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_context_t  ac = { .z_id = z_id,
                                     .func = analytic,
                                     .input = input,
                                     .free_input = NULL };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim, z_id,
                                        0, /* state flag */
                                        meta_flag,
                                        &ac);

  int  new_id = eqp->n_ic_defs;
  eqp->n_ic_defs += 1;
  BFT_REALLOC(eqp->ic_defs, eqp->n_ic_defs, cs_xdef_t *);
  eqp->ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a boundary condition from an existing \ref cs_xdef_t structure
 *         The lifecycle of the cs_xdef_t structure is now managed by the
 *         current \ref cs_equation_param_t structure.
 *
 * \param[in, out] eqp    pointer to a cs_equation_param_t structure
 * \param[in]      xdef   pointer to the \ref cs_xdef_t structure to transfer
*/
/*----------------------------------------------------------------------------*/

void
cs_equation_add_xdef_bc(cs_equation_param_t        *eqp,
                        cs_xdef_t                  *xdef)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = xdef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given equation structure
 *         z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       z_name    name of the related boundary zone
 * \param[in]       values    pointer to a array storing the values
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_value(cs_equation_param_t         *eqp,
                            const cs_param_bc_type_t     bc_type,
                            const char                  *z_name,
                            cs_real_t                   *values)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int  dim = eqp->dim;
  if (bc_type == CS_PARAM_BC_NEUMANN_FULL ||
      bc_type == CS_PARAM_BC_HMG_NEUMANN)
    dim *= 3;  /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_PARAM_BC_ROBIN) {
    /* FluxNormal + alpha * (u - u_0) = beta => Set (alpha, u_0, beta) */
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
                                          cs_get_bdy_zone_id(z_name),
                                          CS_FLAG_STATE_UNIFORM, /* state flag */
                                          meta_flag,
                                          (void *)values);

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given equation structure
 *         z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       z_name    name of the related boundary zone
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       is_owner  transfer the lifecycle to the cs_xdef_t structure
 *                            (true or false)
 * \param[in]       index     optional pointer to the array index
 *
 * \return a pointer to the new allocated \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_array(cs_equation_param_t        *eqp,
                            const cs_param_bc_type_t    bc_type,
                            const char                 *z_name,
                            cs_flag_t                   loc,
                            cs_real_t                  *array,
                            bool                        is_owner,
                            cs_lnum_t                  *index)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  assert(cs_flag_test(loc, cs_flag_primal_face) ||
         cs_flag_test(loc, cs_flag_primal_vtx)  ||
         cs_flag_test(loc, cs_flag_primal_edge)); /* for circulation */

  int  z_id = cs_get_bdy_zone_id(z_name);

  /* Add a new cs_xdef_t structure */
  cs_xdef_array_context_t  input = {.z_id = z_id,
                                    .stride = eqp->dim,
                                    .loc = loc,
                                    .values = array,
                                    .index = index,
                                    .is_owner = is_owner};

  cs_flag_t  state_flag = 0;
  if (loc == cs_flag_primal_face)
    state_flag = CS_FLAG_STATE_FACEWISE;

  int dim = eqp->dim;
  if (bc_type == CS_PARAM_BC_NEUMANN_FULL ||
      bc_type == CS_PARAM_BC_HMG_NEUMANN)
    dim *= 3;  /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_PARAM_BC_ROBIN) {
    /* FluxNormal = alpha * (u_0 - u) + beta => Set (alpha, beta, u_0) */
    if (eqp->dim == 1)
      dim = 3;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled yet.\n", __func__);
  }

  cs_flag_t  meta_flag = (eqp-> space_scheme == CS_SPACE_SCHEME_LEGACY) ?
    (cs_flag_t)bc_type : cs_cdo_bc_get_flag(bc_type);

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ARRAY,
                                          dim,
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          (void *)&input);

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given equation param structure
 *         ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      bc_type   type of boundary condition to add
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function defining the value
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
*/
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_analytic(cs_equation_param_t        *eqp,
                               const cs_param_bc_type_t    bc_type,
                               const char                 *z_name,
                               cs_analytic_func_t         *analytic,
                               void                       *input)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Set the value for dim */
  int dim = eqp->dim;

  if (bc_type == CS_PARAM_BC_NEUMANN_FULL ||
      bc_type == CS_PARAM_BC_HMG_NEUMANN)
    dim *= 3;  /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_PARAM_BC_CIRCULATION) {
    /* This is a vector-valued equation but the DoF is scalar-valued since
     * it is a circulation associated to each edge */
    if (eqp->dim == 3)
      dim = 1;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled.\n", __func__);
  }

  if (bc_type == CS_PARAM_BC_ROBIN) {
    /* FluxNormal = alpha * (u_0 - u) + beta => Set (alpha, beta, u_0) */
    if (eqp->dim == 1)
      dim = 3;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled yet.\n", __func__);
  }

  int  z_id = cs_get_bdy_zone_id(z_name);

  /* Add a new cs_xdef_t structure */
  cs_xdef_analytic_context_t  ac = {.z_id = z_id,
                                    .func = analytic,
                                    .input = input,
                                    .free_input = NULL};

  cs_flag_t  meta_flag = (eqp-> space_scheme == CS_SPACE_SCHEME_LEGACY) ?
    (cs_flag_t)bc_type : cs_cdo_bc_get_flag(bc_type);

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          dim,
                                          z_id,
                                          0, /* state */
                                          meta_flag,
                                          &ac);

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given cs_equation_param_t structure
 *         ml_name corresponds to the name of a pre-existing cs_mesh_location_t
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      bc_type   type of boundary condition to add
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  location where values are computed
 * \param[in]      func      pointer to cs_dof_func_t function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
*/
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_bc_by_dof_func(cs_equation_param_t        *eqp,
                               const cs_param_bc_type_t    bc_type,
                               const char                 *z_name,
                               cs_flag_t                   loc_flag,
                               cs_dof_func_t              *func,
                               void                       *input)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Set the value for dim */
  int dim = eqp->dim;

  if (bc_type == CS_PARAM_BC_NEUMANN_FULL ||
      bc_type == CS_PARAM_BC_HMG_NEUMANN)
    dim *= 3;  /* vector if scalar eq, tensor if vector eq. */

  if (bc_type == CS_PARAM_BC_CIRCULATION) {
    /* This is a vector-valued equation but the DoF is scalar-valued since
     * it is a circulation associated to each edge */
    if (eqp->dim == 3)
      dim = 1;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled.\n", __func__);
  }

  if (bc_type == CS_PARAM_BC_ROBIN) {
    /* FluxNormal = alpha * (u_0 - u) + beta => Set (alpha, beta, u_0) */
    if (eqp->dim == 1)
      dim = 3;
    else
      bft_error(__FILE__, __LINE__, 0,
                "%s: This situation is not handled yet.\n", __func__);
  }

  int  z_id = cs_get_bdy_zone_id(z_name);

  /* Add a new cs_xdef_t structure */
  cs_xdef_dof_context_t  cx = { .z_id = z_id,
                                .loc = loc_flag,
                                .func = func,
                                .input = input,
                                .free_input = NULL };

  cs_flag_t  meta_flag = (eqp-> space_scheme == CS_SPACE_SCHEME_LEGACY) ?
    (cs_flag_t)bc_type : cs_cdo_bc_get_flag(bc_type);

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_DOF_FUNCTION,
                                          dim,
                                          z_id,
                                          0, /* state */
                                          meta_flag,
                                          &cx);

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return pointer to existing boundary condition definition structure
 *         for the given equation param structure and zone.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 *
 * \return a pointer to the \ref cs_xdef_t structure if present, or NULL
*/
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_find_bc(cs_equation_param_t   *eqp,
                    const char            *z_name)
{
  if (eqp == NULL)
    return NULL;

  int z_id = -2;

  const cs_zone_t  *z = cs_boundary_zone_by_name_try(z_name);
  if (z != NULL)
    z_id = z->id;

  /* Search for given BC. */
  for (int i = 0; i < eqp->n_bc_defs; i++) {
    if (eqp->bc_defs[i]->z_id == z_id) {
      return eqp->bc_defs[i];
    }
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Remove boundary condition from the given equation param structure
 *         for a given zone.
 *
 * If no matching boundary condition is found, the function returns
 * silently.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
*/
/*----------------------------------------------------------------------------*/

void
cs_equation_remove_bc(cs_equation_param_t   *eqp,
                      const char            *z_name)
{
  if (eqp == NULL)
    return;

  int z_id = -2;

  const cs_zone_t  *z = cs_boundary_zone_by_name_try(z_name);
  if (z != NULL)
    z_id = z->id;

  /* Search for given BC. */
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
    for (int i = j+1; i < eqp->n_bc_defs; i++) {
      eqp->bc_defs[i-1] = eqp->bc_defs[i];
    }
    eqp->n_bc_defs -= 1;
    BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a sliding boundary
 *         condition related to the given equation structure
 *         z_name corresponds to the name of a pre-existing cs_zone_t
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the related boundary zone
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_sliding_condition(cs_equation_param_t     *eqp,
                                  const char              *z_name)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (eqp->dim < 3)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid dimension of equation\n",
              __func__);

  /* Add two definitions: one for the normal component and one for the
     tangential component */
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs + 1, cs_xdef_t *);

  cs_xdef_t  *d = NULL;
  cs_real_t  val = 0;

  /* Add the homogeneous Dirichlet on the normal component */
  d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                              1,
                              cs_get_bdy_zone_id(z_name),
                              CS_FLAG_STATE_UNIFORM,  /* state flag */
                              CS_CDO_BC_SLIDING,      /* meta */
                              (void *)&val);

  eqp->bc_defs[eqp->n_bc_defs] = d;
  eqp->n_bc_defs += 1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the Laplacian operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *         Laplacian means div-grad (either for vector-valued or scalar-valued
 *         equations)
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_diffusion(cs_equation_param_t   *eqp,
                          cs_property_t         *property)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  eqp->flag |= CS_EQUATION_DIFFUSION;
  eqp->diffusion_property = property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the curl-curl operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 * \param[in]      inversion  > 0: true, false otherwise
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_curlcurl(cs_equation_param_t   *eqp,
                         cs_property_t         *property,
                         int                    inversion)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == NULL)
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
 * \brief  Associate a new term related to the grad-div operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_graddiv(cs_equation_param_t   *eqp,
                        cs_property_t         *property)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  eqp->flag |= CS_EQUATION_GRADDIV;
  eqp->graddiv_property = property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the time derivative operator for
 *         the equation associated to the given \ref cs_equation_param_t
 *         structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_time(cs_equation_param_t   *eqp,
                     cs_property_t         *property)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  eqp->flag |= CS_EQUATION_UNSTEADY;
  eqp->time_property = property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the advection operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      adv_field  pointer to a cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_advection(cs_equation_param_t   *eqp,
                          cs_adv_field_t        *adv_field)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (adv_field == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: Stop adding an empty advection field.",
              __func__, eqp->name);

  eqp->flag |= CS_EQUATION_CONVECTION;
  eqp->adv_field = adv_field;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a scaling property to the advection
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_advection_scaling_property(cs_equation_param_t   *eqp,
                                           cs_property_t         *property)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq. %s: Stop adding an empty property.",
              __func__, eqp->name);

  eqp->adv_scaling_property = property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a new term related to the reaction operator for the
 *         equation associated to the given \ref cs_equation_param_t structure
 *
 * \param[in, out] eqp        pointer to a cs_equation_param_t structure
 * \param[in]      property   pointer to a cs_property_t structure
 *
 * \return the id related to the reaction term
 */
/*----------------------------------------------------------------------------*/

int
cs_equation_add_reaction(cs_equation_param_t   *eqp,
                         cs_property_t         *property)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (property == NULL)
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
 * \brief  Add a new source term by initializing a cs_xdef_t structure.
 *         Case of a definition by a constant value
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or
 *                            "" all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_val(cs_equation_param_t    *eqp,
                                   const char             *z_name,
                                   cs_real_t              *val)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = cs_get_vol_zone_id(z_name);

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
 * \brief  Add a new source term by initializing a cs_xdef_t structure.
 *         Case of a definition by an analytical function
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      func      pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_analytic(cs_equation_param_t    *eqp,
                                        const char             *z_name,
                                        cs_analytic_func_t     *func,
                                        void                   *input)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Define a flag according to the kind of space discretization */
  cs_flag_t  state_flag = 0, meta_flag = 0;
  cs_source_term_set_default_flag(eqp->space_scheme, &state_flag, &meta_flag);

  int z_id = cs_get_vol_zone_id(z_name);
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_context_t  ac = { .z_id = z_id,
                                     .func = func,
                                     .input = input,
                                     .free_input = NULL };

  /* Add a new cs_xdef_t structure */
  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &ac);

  /* Default setting for quadrature is different in this case */
  cs_xdef_set_quadrature(d, CS_QUADRATURE_BARY_SUBDIV);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new source term by initializing a cs_xdef_t structure.
 *         Case of a definition by a DoF function
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc_flag  location of element ids given as parameter
 * \param[in]      func      pointer to a DoF function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_dof_func(cs_equation_param_t    *eqp,
                                        const char             *z_name,
                                        cs_flag_t               loc_flag,
                                        cs_dof_func_t          *func,
                                        void                   *input)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = cs_get_vol_zone_id(z_name);

  /* Define a flag according to the kind of space discretization */
  cs_flag_t  state_flag = 0, meta_flag = 0;
  cs_source_term_set_default_flag(eqp->space_scheme, &state_flag, &meta_flag);

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_dof_context_t  context = { .func = func,
                                     .input = input,
                                     .free_input = NULL,
                                     .loc = loc_flag };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_DOF_FUNCTION,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &context);

  /* Default setting for quadrature is different in this case */
  cs_xdef_set_quadrature(d, CS_QUADRATURE_BARY_SUBDIV);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new source term by initializing a cs_xdef_t structure.
 *         Case of a definition by an array.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc       information to know where are located values
 * \param[in]      array     pointer to an array
 * \param[in]      is_owner  transfer the lifecycle to the cs_xdef_t structure
 *                           (true or false)
 * \param[in]      index     optional pointer to the array index
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_array(cs_equation_param_t    *eqp,
                                     const char             *z_name,
                                     cs_flag_t               loc,
                                     cs_real_t              *array,
                                     bool                    is_owner,
                                     cs_lnum_t              *index)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = cs_get_vol_zone_id(z_name);

  /* Define a flag according to the kind of space discretization */
  cs_flag_t  state_flag = 0, meta_flag = 0;
  cs_source_term_set_default_flag(eqp->space_scheme, &state_flag, &meta_flag);

  if (cs_flag_test(loc, cs_flag_primal_vtx) == true)
    state_flag = CS_FLAG_STATE_POTENTIAL; /* erase predefined settings */
  else if (cs_flag_test(loc, cs_flag_primal_cell) == true)
    state_flag |= CS_FLAG_STATE_CELLWISE;

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_array_context_t  ctxt = {.z_id = z_id,
                                   .stride = eqp->dim,
                                   .loc = loc,
                                   .values = array,
                                   .is_owner = is_owner,
                                   .index = index };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ARRAY,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        (void *)&ctxt);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new volume mass injection definition source term by
 *         initializing a cs_xdef_t structure, using a constant value.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_value(cs_equation_param_t  *eqp,
                                               const char           *z_name,
                                               double               *val)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = cs_get_vol_zone_id(z_name);

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
 * \brief  Add a new volume mass injection definition source term by
 *         initializing a cs_xdef_t structure, using a constant quantity
 *         distributed over the associated zone's volume.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      quantity  pointer to quantity to distribute over the zone
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_qov(cs_equation_param_t  *eqp,
                                             const char           *z_name,
                                             double               *quantity)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = cs_get_vol_zone_id(z_name);

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
 * \brief  Add a new volume mass injection definition source term by
 *         initializing a cs_xdef_t structure, using an analytical function.
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      func      pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_volume_mass_injection_by_analytic(cs_equation_param_t   *eqp,
                                                  const char            *z_name,
                                                  cs_analytic_func_t    *func,
                                                  void                  *input)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = cs_get_vol_zone_id(z_name);

  cs_flag_t  state_flag = 0, meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_context_t  ac = { .z_id = z_id,
                                     .func = func,
                                     .input = input,
                                     .free_input = NULL };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &ac);

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
 * \brief  Add an enforcement of the value of degrees of freedom located at
 *         mesh vertices.
 *         The spatial discretization scheme for the given equation has to be
 *         CDO-Vertex based or CDO-Vertex+Cell-based schemes.
 *
 *         One assumes that values are interlaced if eqp->dim > 1
 *         ref_value or elt_values has to be defined. If both parameters are
 *         defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp         pointer to a cs_equation_param_t structure
 * \param[in]      n_elts      number of vertices to enforce
 * \param[in]      elt_ids     list of vertices
 * \param[in]      ref_value   ignored if NULL
 * \param[in]      elt_values  list of associated values, ignored if NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_enforce_vertex_dofs(cs_equation_param_t    *eqp,
                                cs_lnum_t               n_elts,
                                const cs_lnum_t         elt_ids[],
                                const cs_real_t         ref_value[],
                                const cs_real_t         elt_values[])
{
  if (n_elts < 1)
    return; /* Nothing to do */

  /* Sanity checks */
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_CELLS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: Two types of enforcement are requested (by DoFs and"
              " by cells).\n", __func__, eqp->name);
  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVB &&
      eqp->space_scheme != CS_SPACE_SCHEME_CDOVCB)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: Invalid space scheme.\n"
              "This should be a vertex-based one.", __func__, eqp->name);
  if (ref_value == NULL && elt_values == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: No enforcement value.\n",
              __func__, eqp->name);

  if (eqp->n_enforced_dofs > 0) { /* Reset the selection of DoFs */
    eqp->n_enforced_dofs = 0;
    BFT_FREE(eqp->enforced_dof_ids);
    BFT_FREE(eqp->enforced_dof_values);
  }
  if (eqp->n_enforced_cells > 0) { /* Reset the selection of cells */
    eqp->n_enforced_cells = 0;
    BFT_FREE(eqp->enforced_cell_ids);
    BFT_FREE(eqp->enforced_cell_values);
  }

  eqp->flag |= CS_EQUATION_FORCE_VALUES;
  eqp->n_enforced_dofs = n_elts;
  eqp->enforcement_type = CS_EQUATION_ENFORCE_BY_DOFS;

  BFT_MALLOC(eqp->enforced_dof_ids, n_elts, cs_lnum_t);
  memcpy(eqp->enforced_dof_ids, elt_ids, n_elts*sizeof(cs_lnum_t));

  if (elt_values == NULL) { /* Use a reference value for all enforced DoFs */

    eqp->enforcement_type |= CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE;
    for (int i = 0; i < eqp->dim; i++)
      eqp->enforcement_ref_value[i] = ref_value[i];

  }
  else {

    /* Copy user-defined data in the structure */
    const cs_lnum_t  size = eqp->dim*n_elts;
    BFT_MALLOC(eqp->enforced_dof_values, size, cs_real_t);
    memcpy(eqp->enforced_dof_values, elt_values, size*sizeof(cs_real_t));

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add an enforcement of the value related to the degrees of freedom
 *         associated to the list of selected cells.
 *
 *         One assumes that values are interlaced if eqp->dim > 1
 *         ref_value or elt_values has to be defined. If both parameters are
 *         defined, one keeps the definition in elt_values
 *
 * \param[in, out] eqp         pointer to a cs_equation_param_t structure
 * \param[in]      n_elts      number of selected cells
 * \param[in]      elt_ids     list of cell ids
 * \param[in]      ref_value   ignored if NULL
 * \param[in]      elt_values  list of associated values, ignored if NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_enforce_value_on_cell_selection(cs_equation_param_t  *eqp,
                                            cs_lnum_t             n_elts,
                                            const cs_lnum_t       elt_ids[],
                                            const cs_real_t       ref_value[],
                                            const cs_real_t       elt_values[])
{
  if (n_elts < 1)
    return; /* Nothing to do */

  /* Sanity checks */
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);
  if (eqp->enforcement_type & CS_EQUATION_ENFORCE_BY_DOFS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: Two types of enforcement are requested (by DoFs and"
              " by cells).\n", __func__, eqp->name);
  if (ref_value == NULL && elt_values == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Eq: %s: No enforcement value.\n",
              __func__, eqp->name);

  if (eqp->n_enforced_dofs > 0) { /* Reset the selection of DoFs */
    eqp->n_enforced_dofs = 0;
    BFT_FREE(eqp->enforced_dof_ids);
    BFT_FREE(eqp->enforced_dof_values);
  }
  if (eqp->n_enforced_cells > 0) { /* Reset the selection of cells */
    eqp->n_enforced_cells = 0;
    BFT_FREE(eqp->enforced_cell_ids);
    BFT_FREE(eqp->enforced_cell_values);
  }

  eqp->n_enforced_cells = n_elts;
  eqp->enforcement_type = CS_EQUATION_ENFORCE_BY_CELLS;
  eqp->flag |= CS_EQUATION_FORCE_VALUES;

  BFT_MALLOC(eqp->enforced_cell_ids, n_elts, cs_lnum_t);
  memcpy(eqp->enforced_cell_ids, elt_ids, n_elts*sizeof(cs_lnum_t));

  if (elt_values == NULL) { /* Use a reference value for all enforced DoFs */

    eqp->enforcement_type |= CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE;
    for (int i = 0; i < eqp->dim; i++)
      eqp->enforcement_ref_value[i] = ref_value[i];

  }
  else {

    /* Copy user-defined data in the structure */
    const cs_lnum_t  size = eqp->dim*n_elts;
    BFT_MALLOC(eqp->enforced_cell_values, size, cs_real_t);
    memcpy(eqp->enforced_cell_values, elt_values, size*sizeof(cs_real_t));

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
