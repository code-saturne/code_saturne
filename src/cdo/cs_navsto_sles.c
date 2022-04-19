/*============================================================================
 * Functions to handle SLES structures used during the resolution of the
 * Navier-Stokes system of equations
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>

#include "cs_fp_exception.h"
#include "cs_navsto_param.h"

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_navsto_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

 /*!
   \file cs_navsto_sles.c

   \brief  Functions to handle SLES structures used during the resolution of
           the Navier-Stokes system of equations
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private variables
 *============================================================================*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_HYPRE)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced parameters for the AMG related to the velocity field
 *         when BoomerAMG from the HYPRE library is used
 */
/*----------------------------------------------------------------------------*/

static void
_setup_velocity_boomeramg(void)
{
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL, "-pc_Ux_hypre_boomeramg_coarsen_type", "HMIS");
  PetscOptionsSetValue(NULL, "-pc_Ux_hypre_boomeramg_interp_type", "ext+i-cc");
  PetscOptionsSetValue(NULL, "-pc_Ux_hypre_boomeramg_agg_nl", "2");
  PetscOptionsSetValue(NULL, "-pc_Ux_hypre_boomeramg_P_max", "4");
  PetscOptionsSetValue(NULL, "-pc_Ux_hypre_boomeramg_strong_threshold", "0.5");
  PetscOptionsSetValue(NULL, "-pc_Ux_hypre_boomeramg_no_CF", "");

  PetscOptionsSetValue(NULL, "-pc_Uy_hypre_boomeramg_coarsen_type", "HMIS");
  PetscOptionsSetValue(NULL, "-pc_Uy_hypre_boomeramg_interp_type", "ext+i-cc");
  PetscOptionsSetValue(NULL, "-pc_Uy_hypre_boomeramg_agg_nl", "2");
  PetscOptionsSetValue(NULL, "-pc_Uy_hypre_boomeramg_P_max", "4");
  PetscOptionsSetValue(NULL, "-pc_Uy_hypre_boomeramg_strong_threshold", "0.5");
  PetscOptionsSetValue(NULL, "-pc_Uy_hypre_boomeramg_no_CF", "");

  PetscOptionsSetValue(NULL, "-pc_Uz_hypre_boomeramg_coarsen_type", "HMIS");
  PetscOptionsSetValue(NULL, "-pc_Uz_hypre_boomeramg_interp_type", "ext+i-cc");
  PetscOptionsSetValue(NULL, "-pc_Uz_hypre_boomeramg_agg_nl", "2");
  PetscOptionsSetValue(NULL, "-pc_Uz_hypre_boomeramg_P_max", "4");
  PetscOptionsSetValue(NULL, "-pc_Uz_hypre_boomeramg_strong_threshold", "0.5");
  PetscOptionsSetValue(NULL, "-pc_Uz_hypre_boomeramg_no_CF", "");
#else
  PetscOptionsSetValue("-pc_Ux_hypre_boomeramg_coarsen_type","HMIS");
  PetscOptionsSetValue("-pc_Ux_hypre_boomeramg_interp_type","ext+i-cc");
  PetscOptionsSetValue("-pc_Ux_hypre_boomeramg_agg_nl","2");
  PetscOptionsSetValue("-pc_Ux_hypre_boomeramg_P_max","4");
  PetscOptionsSetValue("-pc_Ux_hypre_boomeramg_strong_threshold","0.5");
  PetscOptionsSetValue("-pc_Ux_hypre_boomeramg_no_CF","");

  PetscOptionsSetValue("-pc_Uy_hypre_boomeramg_coarsen_type","HMIS");
  PetscOptionsSetValue("-pc_Uy_hypre_boomeramg_interp_type","ext+i-cc");
  PetscOptionsSetValue("-pc_Uy_hypre_boomeramg_agg_nl","2");
  PetscOptionsSetValue("-pc_Uy_hypre_boomeramg_P_max","4");
  PetscOptionsSetValue("-pc_Uy_hypre_boomeramg_strong_threshold","0.5");
  PetscOptionsSetValue("-pc_Uy_hypre_boomeramg_no_CF","");

  PetscOptionsSetValue("-pc_Uz_hypre_boomeramg_coarsen_type","HMIS");
  PetscOptionsSetValue("-pc_Uz_hypre_boomeramg_interp_type","ext+i-cc");
  PetscOptionsSetValue("-pc_Uz_hypre_boomeramg_agg_nl","2");
  PetscOptionsSetValue("-pc_Uz_hypre_boomeramg_P_max","4");
  PetscOptionsSetValue("-pc_Uz_hypre_boomeramg_strong_threshold","0.5");
  PetscOptionsSetValue("-pc_Uz_hypre_boomeramg_no_CF","");
#endif
}
#endif  /* PETSC_HAVE_HYPRE */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Setup advanced parameters for the AMG related to the velocity field
 *         when GAMG from the PETSc library is used
 */
/*----------------------------------------------------------------------------*/

static void
_setup_velocity_gamg(void)
{
#if PETSC_VERSION_GE(3,7,0)
  //  PetscOptionsSetValue(NULL, "-gamg_est_ksp_type cg"
  PetscOptionsSetValue(NULL, "-Ux_mg_levels_ksp_type", "richardson");
  PetscOptionsSetValue(NULL, "-Ux_mg_levels_pc_type", "sor");
  PetscOptionsSetValue(NULL, "-Ux_mg_levels_ksp_max_it", "1");
  PetscOptionsSetValue(NULL, "-pc_Ux_gamg_threshold", "0.02");
  PetscOptionsSetValue(NULL, "-pc_Ux_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue(NULL, "-pc_Ux_gamg_square_graph", "4");

  PetscOptionsSetValue(NULL, "-Uy_mg_levels_ksp_type", "richardson");
  PetscOptionsSetValue(NULL, "-Uy_mg_levels_pc_type", "sor");
  PetscOptionsSetValue(NULL, "-Uy_mg_levels_ksp_max_it", "1");
  PetscOptionsSetValue(NULL, "-pc_Uy_gamg_threshold", "0.02");
  PetscOptionsSetValue(NULL, "-pc_Uy_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue(NULL, "-pc_Uy_gamg_square_graph", "4");

  PetscOptionsSetValue(NULL, "-Uz_mg_levels_ksp_type", "richardson");
  PetscOptionsSetValue(NULL, "-Uz_mg_levels_pc_type", "sor");
  PetscOptionsSetValue(NULL, "-Uz_mg_levels_ksp_max_it", "1");
  PetscOptionsSetValue(NULL, "-pc_Uz_gamg_threshold", "0.02");
  PetscOptionsSetValue(NULL, "-pc_Uz_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue(NULL, "-pc_Uz_gamg_square_graph", "4");
#else
  PetscOptionsSetValue("-Ux_mg_levels_ksp_type", "richardson");
  PetscOptionsSetValue("-Ux_mg_levels_pc_type", "sor");
  PetscOptionsSetValue("-Ux_mg_levels_ksp_max_it", "1");
  PetscOptionsSetValue("-pc_Ux_gamg_threshold", "0.02");
  PetscOptionsSetValue("-pc_Ux_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue("-pc_Ux_gamg_square_graph", "4");

  PetscOptionsSetValue("-Uy_mg_levels_ksp_type", "richardson");
  PetscOptionsSetValue("-Uy_mg_levels_pc_type", "sor");
  PetscOptionsSetValue("-Uy_mg_levels_ksp_max_it", "1");
  PetscOptionsSetValue("-pc_Uy_gamg_threshold", "0.02");
  PetscOptionsSetValue("-pc_Uy_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue("-pc_Uy_gamg_square_graph", "4");

  PetscOptionsSetValue("-Uz_mg_levels_ksp_type", "richardson");
  PetscOptionsSetValue("-Uz_mg_levels_pc_type", "sor");
  PetscOptionsSetValue("-Uz_mg_levels_ksp_max_it", "1");
  PetscOptionsSetValue("-pc_Uz_gamg_threshold", "0.02");
  PetscOptionsSetValue("-pc_Uz_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue("-pc_Uz_gamg_square_graph", "4");
#endif
}
#endif  /* HAVE_PETSC */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * User function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User-defined algorithm to solve a saddle point problem (the system is
 *        stored in a hybrid way). Please refer to cs_saddle_system_t structure
 *        and cs_saddle_block_precond_t structure definitions.
 *
 * \param[in]      nslesp  pointer to a cs_navsto_param_sles_t structure
 * \param[in]      ssys    pointer to a cs_saddle_system_t structure
 * \param[in, out] sbp     Block-preconditioner for the Saddle-point problem
 * \param[in, out] x1      array for the first part
 * \param[in, out] x2      array for the second part
 * \param[in, out] algo    pointer to a cs_iter_algo_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_navsto_sles_solve(const cs_navsto_param_sles_t    *nslesp,
                          cs_saddle_system_t              *ssys,
                          cs_saddle_block_precond_t       *sbp,
                          cs_real_t                       *x1,
                          cs_real_t                       *x2,
                          cs_iter_algo_t                  *algo)
{
  CS_UNUSED(nslesp);
  CS_UNUSED(ssys);
  CS_UNUSED(sbp);
  CS_UNUSED(x1);
  CS_UNUSED(x2);
  CS_UNUSED(algo);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of multiplicative block preconditioner for a CG
 *
 * \param[in, out] context    pointer to optional (untyped) value or structure
 * \param[in, out] ksp_struct pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

void
cs_navsto_sles_amg_block_hook(void     *context,
                              void     *ksp_struct)
{
  const cs_navsto_param_t  *nsp = (cs_navsto_param_t *)context;

  KSP  ksp = ksp_struct;
  cs_equation_param_t  *eqp = cs_navsto_param_get_velocity_param(nsp);
  cs_param_sles_t  *slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  cs_navsto_param_model_t  model = nsp->model;
  if (model == CS_NAVSTO_MODEL_STOKES)
    KSPSetType(ksp, KSPFCG);
  else
    KSPSetType(ksp, KSPFGMRES);

  /* Set KSP tolerances */

  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp->eps,        /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp->n_max_iter); /* max number of iterations */

  switch (slesp->resnorm_type) {

  case CS_PARAM_RESNORM_NORM2_RHS: /* Try to have "true" norm */
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    break;
  case CS_PARAM_RESNORM_NONE:
    KSPSetNormType(ksp, KSP_NORM_NONE);
    break;
  default:
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    break;

  }

  PC  pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCFIELDSPLIT);
  PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE);

  /* Apply modifications to the KSP structure */

  PetscInt  id, n_split;
  KSP  *uvw_subksp;

  PCFieldSplitSetBlockSize(pc, 3);
  id = 0;
  PCFieldSplitSetFields(pc, "Ux", 1, &id, &id);
  id = 1;
  PCFieldSplitSetFields(pc, "Uy", 1, &id, &id);
  id = 2;
  PCFieldSplitSetFields(pc, "Uz", 1, &id, &id);

  PCSetFromOptions(pc);
  PCSetUp(pc);
  PCFieldSplitGetSubKSP(pc, &n_split, &uvw_subksp);
  assert(n_split == 3);

  switch(slesp->amg_type) {

  case CS_PARAM_AMG_HYPRE_BOOMER_V:
  case CS_PARAM_AMG_HYPRE_BOOMER_W:
#if defined(PETSC_HAVE_HYPRE)
    _setup_velocity_boomeramg();
#else
    _setup_velocity_gamg();
#endif
    break;

  case CS_PARAM_AMG_PETSC_PCMG:
  case CS_PARAM_AMG_PETSC_GAMG_V:
  case CS_PARAM_AMG_PETSC_GAMG_W:
    _setup_velocity_gamg();
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid choice of AMG type.\n", __func__);
    break;
  }

  PC _pc;
  for (id = 0; id < 3; id++) {

    KSP  _ksp = uvw_subksp[id];
    KSPSetType(_ksp, KSPPREONLY);
    KSPGetPC(_ksp, &_pc);

    switch(slesp->amg_type) {

    case CS_PARAM_AMG_HYPRE_BOOMER_V:
#if defined(PETSC_HAVE_HYPRE)
      PCSetType(_pc, PCHYPRE);
      PCHYPRESetType(_pc, "boomeramg");
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_cycle_type","V");
#else
      PCSetType(_pc, PCGAMG);
      PCGAMGSetType(_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(_pc, 1);
      PCMGSetCycleType(_pc, PC_MG_CYCLE_V);
#endif
      break;

    case CS_PARAM_AMG_HYPRE_BOOMER_W:
#if defined(PETSC_HAVE_HYPRE)
      PCSetType(_pc, PCHYPRE);
      PCHYPRESetType(_pc, "boomeramg");
      PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_cycle_type","W");
#else
      PCSetType(_pc, PCGAMG);
      PCGAMGSetType(_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(_pc, 1);
      PCMGSetCycleType(_pc, PC_MG_CYCLE_W);
#endif
      break;

    case CS_PARAM_AMG_PETSC_PCMG:
      PCSetType(_pc, PCMG);
      break;

    case CS_PARAM_AMG_PETSC_GAMG_V:
      PCSetType(_pc, PCGAMG);
      PCGAMGSetType(_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(_pc, 1);
      PCMGSetCycleType(_pc, PC_MG_CYCLE_V);
      break;

    case CS_PARAM_AMG_PETSC_GAMG_W:
      PCSetType(_pc, PCGAMG);
      PCGAMGSetType(_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(_pc, 1);
      PCMGSetCycleType(_pc, PC_MG_CYCLE_W);
      break;

    default:
      break;
    }

  }

  /* User function for additional settings */

  cs_user_sles_petsc_hook(context, ksp);

  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */

  if (!slesp->setup_done) {

    cs_sles_petsc_log_setup(ksp);
    slesp->setup_done = true;

  }

  PetscFree(uvw_subksp);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif  /* HAVE_PETSC */

/*----------------------------------------------------------------------------*/

END_C_DECLS
