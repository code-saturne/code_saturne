/*============================================================================
 * Functions dedicated to the linear algebra settings and operations in case
 * of CDO face-based schemes with a monolithic velocity-pressure coupling
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

#if defined(HAVE_PETSC)
#include <petscversion.h>
#include <petscksp.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation.h"
#include "cs_fp_exception.h"
#include "cs_navsto_coupling.h"
#include "cs_sles.h"
#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_monolithic_sles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_monolithic_sles.c
 *
 * \brief Functions dedicated to to the linear algebra settings and operations
 *        in case of CDO face-based schemes with a monolithic velocity-pressure
 *        coupling
 */

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_MONOLITHIC_SLES_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_range_set_t         *cs_shared_range_set;

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
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_coarsen_type", "HMIS");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_interp_type", "ext+i-cc");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_agg_nl", "2");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_P_max", "4");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_strong_threshold", "0.5");
  PetscOptionsSetValue(NULL,
                       "-pc_velocity_hypre_boomeramg_no_CF", "");
#else
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_coarsen_type","HMIS");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_interp_type","ext+i-cc");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_agg_nl","2");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_P_max","4");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_strong_threshold","0.5");
  PetscOptionsSetValue("-pc_velocity_hypre_boomeramg_no_CF","");
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
  PetscOptionsSetValue(NULL, "-mg_velocity_levels_ksp_type", "richardson");
  PetscOptionsSetValue(NULL, "-mg_velocity_levels_pc_type", "sor");
  PetscOptionsSetValue(NULL, "-mg_velocity_levels_ksp_max_it", "1");
  PetscOptionsSetValue(NULL, "-pc_velocity_gamg_threshold", "0.02");
  PetscOptionsSetValue(NULL, "-pc_velocity_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue(NULL, "-pc_velocity_gamg_square_graph", "4");
#else
  PetscOptionsSetValue("-mg_velocity_levels_ksp_type", "richardson");
  PetscOptionsSetValue("-mg_velocity_levels_pc_type", "sor");
  PetscOptionsSetValue("-mg_velocity_levels_ksp_max_it", "1");
  PetscOptionsSetValue("-pc_velocity_gamg_threshold", "0.02");
  PetscOptionsSetValue("-pc_velocity_gamg_reuse_interpolation", "TRUE");
  PetscOptionsSetValue("-pc_velocity_gamg_square_graph", "4");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generate IndexSet for the PETSc FieldSplit preconditioner
 *
 * \param[in, out]  isp     IndexSet for the pressure DoFs
 * \param[in, out]  isv     IndexSet for the velocity DoFs
 */
/*----------------------------------------------------------------------------*/

static void
_build_is_for_fieldsplit(IS   *isp,
                         IS   *isv)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_range_set_t  *rset = cs_shared_range_set;

  PetscInt  n_faces = quant->n_faces;
  PetscInt  n_cells = quant->n_cells;
  PetscInt  *indices = NULL;

  PetscMalloc1(3*n_faces, &indices);

  /* IndexSet for the velocity DoFs */
  if (rset->n_elts[0] == rset->n_elts[1]) {

    for (PetscInt i = 0; i < 3*n_faces; i++)
      indices[i] = rset->g_id[i];
    ISCreateGeneral(PETSC_COMM_WORLD, 3*n_faces, indices, PETSC_COPY_VALUES,
                    isv);

  }
  else {

    PetscInt  n_velocity_elts = 0;
    for (PetscInt i = 0; i < 3*n_faces; i++) {
      cs_gnum_t  g_id = rset->g_id[i];
      if (g_id >= rset->l_range[0] && g_id < rset->l_range[1])
        indices[n_velocity_elts++] = g_id;
    }
    ISCreateGeneral(PETSC_COMM_WORLD, n_velocity_elts, indices,
                    PETSC_COPY_VALUES, isv);

  }

  /* IndexSet for the velocity DoFs
   * Pressure unknowns are located at cell centers so the treatment should be
   * the same in sequential and parallel computation*/
  for (PetscInt i = 0; i < n_cells; i++)
    indices[i] = rset->g_id[i + 3*n_faces];
  ISCreateGeneral(PETSC_COMM_WORLD, n_cells, indices, PETSC_COPY_VALUES,
                  isp);

  PetscFree(indices);
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of additive block preconditioner for a GMRES
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_additive_amg_gmres_hook(void     *context,
                         Mat       a,
                         KSP       ksp)
{
  IS  isv, isp;

  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  const int  n_max_restart = 30;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  KSPSetType(ksp, KSPFGMRES);
  KSPGMRESSetRestart(ksp, n_max_restart);

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp.eps,         /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  /* Try to have "true" norm */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

#if 0 /* JB: TEST TO PERFORM IN 3D*/
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL,
                  "-fieldsplit_velocity_pc_hypre_boomeramg_strong_threshold",
                  "0.5");
#else
  PetscOptionsSetValue(
                  "-fieldsplit_velocity_pc_hypre_boomeramg_strong_threshold",
                  "0.5");
#endif
#endif

  /* Apply modifications to the KSP structure */
  PC up_pc, u_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_ADDITIVE);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPPREONLY);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCJACOBI);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  KSP  u_ksp = up_subksp[0];
  KSPSetType(u_ksp, KSPPREONLY);
  KSPGetPC(u_ksp, &u_pc);

  switch(slesp.amg_type) {

  case CS_PARAM_AMG_HYPRE_BOOMER:
#if defined(PETSC_HAVE_HYPRE)
    PCSetType(u_pc, PCHYPRE);
    PCHYPRESetType(u_pc, "boomeramg");

    _setup_velocity_boomeramg();
#else
    _setup_velocity_gamg();
#endif
    break;

  case CS_PARAM_AMG_PETSC_PCMG:
  case CS_PARAM_AMG_PETSC_GAMG:
    PCSetType(u_pc, PCGAMG);
    PCGAMGSetType(u_pc, PCGAMGAGG);
    PCGAMGSetNSmooths(u_pc, 1);

    _setup_velocity_gamg();
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid choice of AMG type.\n", __func__);
    break;
  }

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);

  KSPSetFromOptions(u_ksp);
  KSPSetUp(u_ksp);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of diagonal Schur preconditioner by block for a GMRES
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_diag_schur_gmres_hook(void     *context,
                       Mat       a,
                       KSP       ksp)
{
  IS  isv, isp;

  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  const int  n_max_restart = 30;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  KSPSetType(ksp, KSPFGMRES);
  KSPGMRESSetRestart(ksp, n_max_restart);

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp.eps,         /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  /* Try to have "true" norm */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

  /* Apply modifications to the KSP structure */
  PC up_pc, u_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_SCHUR);
  PCFieldSplitSetSchurFactType(up_pc, PC_FIELDSPLIT_SCHUR_FACT_DIAG);
  PCFieldSplitSetSchurPre(up_pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPMINRES);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCNONE);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  KSP  u_ksp = up_subksp[0];
  KSPSetType(u_ksp, KSPCG);
  KSPGetPC(u_ksp, &u_pc);

#if defined(PETSC_HAVE_HYPRE)
  PCSetType(u_pc, PCHYPRE);
  PCHYPRESetType(u_pc, "boomeramg");

  _setup_velocity_boomeramg();
#else
  PCSetType(u_pc, PCGAMG);
  PCGAMGSetType(u_pc, PCGAMGAGG);
  PCGAMGSetNSmooths(u_pc, 1);

  _setup_velocity_gamg();
#endif

  KSPSetTolerances(u_ksp,
                   slesp.eps,   /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   5);          /* max number of iterations */


  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);

  KSPSetFromOptions(u_ksp);
  KSPSetUp(u_ksp);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of upper Schur preconditioner by block for a GMRES
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_upper_schur_gmres_hook(void     *context,
                        Mat       a,
                        KSP       ksp)
{
  IS  isv, isp;

  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  const int  n_max_restart = 30;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  KSPSetType(ksp, KSPFGMRES);
  KSPGMRESSetRestart(ksp, n_max_restart);

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp.eps,         /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  /* Try to have "true" norm */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

  /* Apply modifications to the KSP structure */
  PC up_pc, u_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_SCHUR);
  PCFieldSplitSetSchurFactType(up_pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER);
  PCFieldSplitSetSchurPre(up_pc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPMINRES);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCNONE);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  KSP  u_ksp = up_subksp[0];
  KSPSetType(u_ksp, KSPCG);
  KSPGetPC(u_ksp, &u_pc);
#if defined(PETSC_HAVE_HYPRE)
  PCSetType(u_pc, PCHYPRE);
  PCHYPRESetType(u_pc, "boomeramg");

  _setup_velocity_boomeramg();
#else
  PCSetType(u_pc, PCGAMG);
  PCGAMGSetType(u_pc, PCGAMGAGG);
  PCGAMGSetNSmooths(u_pc, 1);

  _setup_velocity_gamg();
#endif

  KSPSetTolerances(u_ksp,
                   slesp.eps,   /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   5);          /* max number of iterations */

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);

  KSPSetFromOptions(u_ksp);
  KSPSetUp(u_ksp);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

#if PETSC_VERSION_GE(3,11,0)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of GKB as a solver with CG(Boomer) as inner solver for the
 *         velocity block.
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_gkb_hook(void     *context,
          Mat       a,
          KSP       ksp)
{
  IS  isv, isp;

  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  KSPSetType(ksp, KSPPREONLY);

  /* Apply modifications to the KSP structure */
  PC up_pc, u_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_GKB);

  PCFieldSplitSetGKBTol(up_pc, 10*slesp.eps);
  PCFieldSplitSetGKBMaxit(up_pc, slesp.n_max_iter);
  PCFieldSplitSetGKBNu(up_pc, 0);
  PCFieldSplitSetGKBDelay(up_pc, 5);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  u_ksp = up_subksp[0];

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPSetType(u_ksp, KSPFGMRES);
  KSPGetPC(u_ksp, &u_pc);
#if defined(PETSC_HAVE_HYPRE)
  PCSetType(u_pc, PCHYPRE);
  PCHYPRESetType(u_pc, "boomeramg");

  _setup_velocity_boomeramg();
#else
  PCSetType(u_pc, PCGAMG);
  PCGAMGSetType(u_pc, PCGAMGAGG);
  PCGAMGSetNSmooths(u_pc, 1);

  _setup_velocity_gamg();
#endif

  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(u_ksp,
                   slesp.eps,   /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);

  KSPSetFromOptions(u_ksp);
  KSPSetUp(u_ksp);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of GKB preconditioner. by block for a GMRES
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_gkb_gmres_hook(void     *context,
                Mat       a,
                KSP       ksp)
{
  IS  isv, isp;

  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  KSPSetType(ksp, KSPFGMRES);

  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp.eps,         /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  /* Apply modifications to the KSP structure */
  PC up_pc, u_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_GKB);

  PCFieldSplitSetGKBTol(up_pc, 1e-1);
  PCFieldSplitSetGKBMaxit(up_pc, 100);
  PCFieldSplitSetGKBNu(up_pc, 0);
  PCFieldSplitSetGKBDelay(up_pc, 5);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  u_ksp = up_subksp[0];

  /* Set KSP tolerances */
  KSPSetType(u_ksp, KSPFGMRES);
  KSPGetPC(u_ksp, &u_pc);
#if defined(PETSC_HAVE_HYPRE)
  PCSetType(u_pc, PCHYPRE);
  PCHYPRESetType(u_pc, "boomeramg");

  _setup_velocity_boomeramg();
#else
  PCSetType(u_pc, PCGAMG);
  PCGAMGSetType(u_pc, PCGAMGAGG);
  PCGAMGSetNSmooths(u_pc, 1);

  _setup_velocity_gamg();
#endif

  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(u_ksp,
                   1e-2,    /* relative convergence tolerance */
                   abstol,  /* absolute convergence tolerance */
                   dtol,    /* divergence tolerance */
                   50);     /* max number of iterations */

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);

  KSPSetFromOptions(u_ksp);
  KSPSetUp(u_ksp);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif  /* GKB available only if version >= 3.11 */
#endif  /* HAVE_PETSC */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set pointers to shared structures
 *
 * \param[in]  quant    pointer to additional mesh quantities
 * \param[in]  rset     pointer to a \ref cs_range_set_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_set_shared(const cs_cdo_quantities_t     *quant,
                                    const cs_range_set_t          *rset)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_range_set = rset;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a monolithic
 *         algorithm is used to couple the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_set_sles(const cs_navsto_param_t    *nsp,
                             void                       *context)
{
  cs_navsto_monolithic_t  *nsc = (cs_navsto_monolithic_t *)context;

  assert(nsp != NULL && nsc != NULL);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);
  int  field_id = cs_equation_get_field_id(nsc->momentum);

  mom_eqp->sles_param.field_id = field_id;
  if (mom_eqp->sles_param.amg_type == CS_PARAM_AMG_NONE)
    mom_eqp->sles_param.amg_type = CS_PARAM_AMG_HYPRE_BOOMER;

  /* Initialization must be called before setting options;
     it does not need to be called before calling
     cs_sles_petsc_define(), as this is handled automatically. */

  switch (nsp->sles_strategy) {

  case CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK: /* "Classical" way to set SLES */
    cs_equation_param_set_sles(mom_eqp);
    break;

#if defined(HAVE_PETSC)
  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _additive_amg_gmres_hook,
                         (void *)mom_eqp);
    break;

  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _diag_schur_gmres_hook,
                         (void *)mom_eqp);
    break;

  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _upper_schur_gmres_hook,
                         (void *)mom_eqp);
    break;

#if PETSC_VERSION_GE(3,11,0)    /* Golub-Kahan Bi-diagonalization */
  case CS_NAVSTO_SLES_GKB:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _gkb_hook,
                         (void *)mom_eqp);
    break;

  case CS_NAVSTO_SLES_GKB_GMRES:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _gkb_gmres_hook,
                         (void *)mom_eqp);
    break;
#else
  case CS_NAVSTO_SLES_GKB:
  case CS_NAVSTO_SLES_GKB_GMRES:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc 3.11.x or greater is required with this option.\n",
              __func__, mom_eqp->name);
    break;
#endif

#else
  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please build a version of Code_Saturne with the PETSc support.",
              __func__, mom_eqp->name);
    break;
#endif /* HAVE_PETSC */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n",
              __func__, mom_eqp->name);
  }

      /* Define the level of verbosity for SLES structure */
  if (mom_eqp->sles_param.verbosity > 1) {

    cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

    /* Set verbosity */
    cs_sles_set_verbosity(sles, mom_eqp->sles_param.verbosity);

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
