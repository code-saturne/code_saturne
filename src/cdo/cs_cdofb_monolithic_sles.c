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

#include "cs_blas.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation.h"
#include "cs_fp_exception.h"
#include "cs_navsto_coupling.h"
#include "cs_parall.h"
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

/* GKB advanced settings */

#define CS_GKB_TRUNCATION_THRESHOLD       5

/* Block size for superblock algorithm */

#define CS_SBLOCK_BLOCK_SIZE 60

/* Cache line multiple, in cs_real_t units */

#define CS_CL  (CS_CL_SIZE/8)

/* This structure follow notations given in the article entitled
 * "An iterative generalized Golub-Kahan algorithm for problems in structural
 *  mechanics" by M. Arioli, C. Kruse, U. Ruede and N. Tardieu
 *
 * M space is isomorphic to the velocity space (size = 3.n_faces)
 * N space is isomorphic to the pressure space (size = n_cells)
 */

typedef struct {

  /* Value of the grad-div coefficient */
  cs_real_t                gamma;

  /* Size of spaces */
  cs_lnum_t                n_u_dofs; /* Size of the space M */
  cs_lnum_t                n_p_dofs; /* Size of the space N */

  /* Vector transformation */
  cs_real_t               *b_tilda;  /* Modified RHS */
  cs_real_t               *u_tilda;  /* Modified velocity unknown */

  /* Auxiliary vectors */
  cs_real_t               *q;        /* vector iterates in space N */
  cs_real_t               *d;        /* vector iterates in space N */
  cs_real_t               *d__v;     /* buffer in space N */
  cs_real_t               *dt_q;     /* buffer in space M */
  cs_real_t               *m__v;     /* vector iterates in space M */
  cs_real_t               *v;        /* vector iterates in space M */

  /* Orthogonalization coefficients */
  cs_real_t                alpha;
  cs_real_t                beta;
  cs_real_t                zeta;

  /* Store z_size zeta coefficients */
  int                      z_size;
  cs_real_t               *zeta_array;
  cs_real_t                zeta_square_sum;

  cs_navsto_algo_info_t    info;     /* Information related to the convergence
                                        of the algorithm */

} cs_gkb_builder_t;

/* This structure follow notations given in the article entitled
 * "An iterative generalized Golub-Kahan algorithm for problems in structural
 *  mechanics" by M. Arioli, C. Kruse, U. Ruede and N. Tardieu
 *
 * M space is isomorphic to the velocity space (size = 3.n_faces)
 * N space is isomorphic to the pressure space (size = n_cells)
 */

typedef struct {

  /* Value of the grad-div coefficient */
  cs_real_t                gamma;

  /* Size of spaces */
  cs_lnum_t                n_u_dofs; /* Size of the space M */
  cs_lnum_t                n_p_dofs; /* Size of the space N */

  /* Vector transformation */
  cs_real_t               *b_tilda;  /* Modified RHS */

  /* Auxiliary vectors */
  cs_real_t               *inv_mp;   /* reciprocal of the pressure mass matrix */
  cs_real_t               *res_p;    /* buffer in space N */
  cs_real_t               *d__v;     /* buffer in space N */
  cs_real_t               *rhs;      /* buffer in space M */

  cs_navsto_algo_info_t    info;     /* Information related to the convergence
                                        of the algorithm */

} cs_uza_builder_t;

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_range_set_t         *cs_shared_range_set;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute array index bounds for a local thread.
 *
 * When called inside an OpenMP parallel section, this will return the
 * start an past-the-end indexes for the array range assigned to that thread.
 * In other cases, the start index is 1, and the past-the-end index is n;
 *
 * \param[in]   n     size of array
 * \param[out]  s_id  start index for the current thread
 * \param[out]  e_id  past-the-end index for the current thread
 */
/*----------------------------------------------------------------------------*/

static inline void
_thread_range(cs_lnum_t   n,
              cs_lnum_t  *s_id,
              cs_lnum_t  *e_id)
{
#if defined(HAVE_OPENMP)
  int t_id = omp_get_thread_num();
  int n_t = omp_get_num_threads();
  cs_lnum_t t_n = (n + n_t - 1) / n_t;
  *s_id =  t_id    * t_n;
  *e_id = (t_id+1) * t_n;
  *s_id = cs_align(*s_id, CS_CL);
  *e_id = cs_align(*e_id, CS_CL);
  if (*e_id > n) *e_id = n;
#else
  *s_id = 0;
  *e_id = n;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dot product between two arrays on face unknowns.
 *         One assumes that input arrays are in a "scattered" distribution
 *         So the size should be 3*n_faces.
 *
 * \param[in]       size   size of arrays
 * \param[in, out]  x      first array
 * \param[in, out]  y      second array
 *
 * \return the computed value
 */
/*----------------------------------------------------------------------------*/

static inline cs_real_t
_face_gdot(cs_lnum_t    size,
           cs_real_t    x[],
           cs_real_t    y[])
{
  CS_UNUSED(size);       /* Avoid a compilation warning in during compilation */
  const cs_range_set_t  *rset = cs_shared_range_set;

  assert(size == rset->n_elts[1]);
  assert(size == 3*cs_shared_quant->n_faces);

  /* x and y are scattered arrays. One assumes that values are synchronized
     across ranks (for instance by using a cs_interface_set_sum()) */
  if (cs_glob_n_ranks > 1) {

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,/* type */
                        1,           /* stride (treated as scalar up to now) */
                        x,           /* in: size = n_sles_scatter_elts */
                        x);          /* out: size = n_sles_gather_elts */

    cs_range_set_gather(rset,
                        CS_REAL_TYPE,/* type */
                        1,           /* stride (treated as scalar up to now) */
                        y,           /* in: size = n_sles_scatter_elts */
                        y);          /* out: size = n_sles_gather_elts */

  }

  cs_real_t  result = cs_gdot(rset->n_elts[0], x, y);

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    cs_range_set_scatter(rset,
                         CS_REAL_TYPE,
                         1,
                         x,
                         x);
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE,
                         1,
                         y,
                         y);

  }

  return result;
}

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set command line options for PC according to the kind of
 *        preconditionner
 *
 * \param[in]   slesp      set of parameters for the linear algebra
 */
/*----------------------------------------------------------------------------*/

static PCType
_petsc_get_pc_type(const cs_param_sles_t    slesp)
{
  PCType  pc_type = PCNONE;

  switch (slesp.precond) {

  case CS_PARAM_PRECOND_NONE:
    return PCNONE;

  case CS_PARAM_PRECOND_DIAG:
    return PCJACOBI;

  case CS_PARAM_PRECOND_BJACOB:
    return PCBJACOBI;

  case CS_PARAM_PRECOND_SSOR:
    return PCSOR;

  case CS_PARAM_PRECOND_ICC0:
    return PCICC;

  case CS_PARAM_PRECOND_ILU0:
    return PCILU;

  case CS_PARAM_PRECOND_AS:
    return PCASM;

  case CS_PARAM_PRECOND_AMG:
    {
      switch (slesp.amg_type) {

      case CS_PARAM_AMG_PETSC_GAMG:
        return PCGAMG;
        break;

      case CS_PARAM_AMG_PETSC_PCMG:
        return PCMG;
        break;

      case CS_PARAM_AMG_HYPRE_BOOMER:
#if defined(PETSC_HAVE_HYPRE)
        return PCHYPRE;
#else
        cs_base_warn(__FILE__, __LINE__);
        cs_log_printf(CS_LOG_DEFAULT,
                      "%s: Switch to MG since BoomerAMG is not available.\n",
                      __func__);
        return PCMG;
#endif

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid AMG type for the PETSc library.", __func__);
        break;

      } /* End of switch on the AMG type */

    } /* AMG as preconditioner */
    break;

  case CS_PARAM_PRECOND_AMG_BLOCK:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Preconditioner not interfaced with PETSc.\n"
              " Switch to a non block version", __func__);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Preconditioner not interfaced with PETSc.", __func__);
    break;
  }

  return pc_type;
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of additive block preconditioner for a GMRES
 *
 * \param[in]      slesp     pointer to a set of SLES settings
 * \param[in, out] u_ksp     pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_set_velocity_ksp(const cs_param_sles_t    slesp,
                  KSP                      u_ksp)
{
  PC u_pc;
  KSPGetPC(u_ksp, &u_pc);
  PCType  pc_type = _petsc_get_pc_type(slesp);

  /* Try to have "true" norm */
  KSPSetNormType(u_ksp, KSP_NORM_UNPRECONDITIONED);

  /* Set the solver */
  switch (slesp.solver) {

  case CS_PARAM_ITSOL_FCG:
  case CS_PARAM_ITSOL_CG:
    KSPSetType(u_ksp, KSPCG);
    break;
  case CS_PARAM_ITSOL_BICG:      /* Improved Bi-CG stab */
    KSPSetType(u_ksp, KSPIBCGS);
    break;
  case CS_PARAM_ITSOL_BICGSTAB2: /* BiCGstab2 */
    KSPSetType(u_ksp, KSPBCGSL);
    break;
  case CS_PARAM_ITSOL_MUMPS:     /* Direct solver (factorization) */
#if defined(PETSC_HAVE_MUMPS)
    KSPSetType(u_ksp, KSPPREONLY);
    PCSetType(u_pc, PCLU);
    PCFactorSetMatSolverType(u_pc, MATSOLVERMUMPS);
#else
    bft_error(__FILE__, __LINE__, 0,
              " %s: MUMPS not interfaced with this installation of PETSc.",
              __func__);
#endif
    break;
  case CS_PARAM_ITSOL_GMRES:
    KSPSetType(u_ksp, KSPGMRES);
    break;

  case CS_PARAM_ITSOL_MUMPS_LDLT:     /* Direct solver (factorization) */
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver. Try mumps.",
              __func__);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid solver.", __func__);
    break;

  } /* Switch on solver */

  if (slesp.solver != CS_PARAM_ITSOL_MUMPS)
    PCSetType(u_pc, pc_type);

  /* Additional settings for the preconditioner */
  switch (slesp.precond) {

  case CS_PARAM_PRECOND_AMG:
    switch (slesp.amg_type) {

    case CS_PARAM_AMG_PETSC_GAMG:
      PCGAMGSetType(u_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(u_pc, 1);
      _setup_velocity_gamg();
      break;

    case CS_PARAM_AMG_HYPRE_BOOMER:
#if defined(PETSC_HAVE_HYPRE)
      PCHYPRESetType(u_pc, "boomeramg");
      _setup_velocity_boomeramg();

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
#endif  /* if 0 */

#else
      PCGAMGSetType(u_pc, PCGAMGAGG);
      PCGAMGSetNSmooths(u_pc, 1);
      _setup_velocity_gamg();
#endif
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid AMG type.", __func__);
      break;

    } /* AMG type */
    break;

  default:
    break; /* Nothing else to do */

  } /* Switch on preconditioner */

  /* Set tolerance and number of iterations */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(u_ksp, &rtol, &abstol, &dtol, &maxit);
  rtol = fmin(1e-4, fmax(1e4*slesp.eps, 1e-5));
  KSPSetTolerances(u_ksp,
                   rtol,        /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   500);        /* max number of iterations */

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);

  KSPSetFromOptions(u_ksp);
  KSPSetUp(u_ksp);
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

  /* Apply modifications to the KSP structure */
  PC up_pc, p_pc;

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

  /* Set the velocity block */
  _set_velocity_ksp(slesp, up_subksp[0]);

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
 *         Case of multiplicative block preconditioner for a GMRES
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_multiplicative_gmres_hook(void     *context,
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
  PC up_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_MULTIPLICATIVE);

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

  /* Set the velocity block */
  _set_velocity_ksp(slesp, up_subksp[0]);
  KSP  u_ksp = up_subksp[0];

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
  PC up_pc, p_pc;

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

  /* Set the velocity block */
  _set_velocity_ksp(slesp, up_subksp[0]);

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
  PC up_pc, p_pc;

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

  /* Set the velocity block */
  _set_velocity_ksp(slesp, up_subksp[0]);

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

#if defined(PETSC_HAVE_MUMPS)
/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of MUMPS via PETSc
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_mumps_hook(void     *context,
            Mat       a,
            KSP       ksp)
{
  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  cs_fp_exception_disable_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */

  PC  pc;
  KSPSetType(ksp, KSPPREONLY);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCLU);
  PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);

  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp.eps,   /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  cs_fp_exception_restore_trap(); /* Avoid trouble with a too restrictive
                                     SIGFPE detection */
}
#endif  /* PETSC_HAVE_MUMPS */

#endif  /* HAVE_PETSC */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a GKB builder structure
 *
 * \param[in]  gamma      value of the grad-div coefficient
 * \param[in]  n_u_dofs   number of velocity DoFs (degrees of freedom)
 * \param[in]  n_p_dofs   number of pressure DoFs
 *
 * \return a pointer to a new allocated GKB builder
 */
/*----------------------------------------------------------------------------*/

static cs_gkb_builder_t *
_init_gkb_builder(cs_real_t             gamma,
                  cs_lnum_t             n_u_dofs,
                  cs_lnum_t             n_p_dofs)
{
  cs_gkb_builder_t  *gkb = NULL;

  BFT_MALLOC(gkb, 1, cs_gkb_builder_t);

  gkb->gamma = gamma;
  gkb->n_u_dofs = n_u_dofs;
  gkb->n_p_dofs = n_p_dofs;

  /* Vector transformation */
  BFT_MALLOC(gkb->u_tilda, n_u_dofs, cs_real_t);
  /* Rk: b_tilda stores quantities in space M and N alternatively */
  assert(n_u_dofs >= n_p_dofs);
  BFT_MALLOC(gkb->b_tilda, n_u_dofs, cs_real_t);

  /* Auxiliary vectors */
  BFT_MALLOC(gkb->v, n_u_dofs, cs_real_t);
  memset(gkb->v, 0, n_u_dofs*sizeof(cs_real_t));

  BFT_MALLOC(gkb->q, n_p_dofs, cs_real_t);
  BFT_MALLOC(gkb->d, n_p_dofs, cs_real_t);
  BFT_MALLOC(gkb->d__v, n_p_dofs, cs_real_t);
  BFT_MALLOC(gkb->dt_q, n_u_dofs, cs_real_t);
  BFT_MALLOC(gkb->m__v, n_u_dofs, cs_real_t);

  /* Orthogonalization coefficients */
  gkb->alpha = gkb->beta = gkb->zeta = 0.;

  /* Convergence members */
  if (gamma < 1)
    gkb->z_size = CS_GKB_TRUNCATION_THRESHOLD + 1;
  else if (gamma < 10)
    gkb->z_size = CS_GKB_TRUNCATION_THRESHOLD;
  else if (gamma < 100)
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 1);
  else if (gamma < 1e3)
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 2);
  else if (gamma < 1e4)
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 3);
  else
    gkb->z_size = CS_MAX(1, CS_GKB_TRUNCATION_THRESHOLD - 4);

  BFT_MALLOC(gkb->zeta_array, gkb->z_size, cs_real_t);
  memset(gkb->zeta_array, 0, gkb->z_size*sizeof(cs_real_t));

  gkb->zeta_square_sum = 0.;

  cs_navsto_algo_info_init(&(gkb->info));

  return gkb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a GKB builder structure
 *
 * \param[in, out]  p_gkb   double pointer to a GKB builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_free_gkb_builder(cs_gkb_builder_t   **p_gkb)
{
  cs_gkb_builder_t  *gkb = *p_gkb;

  if (gkb == NULL)
    return;

  BFT_FREE(gkb->b_tilda);
  BFT_FREE(gkb->u_tilda);

  BFT_FREE(gkb->q);
  BFT_FREE(gkb->d);
  BFT_FREE(gkb->d__v);
  BFT_FREE(gkb->dt_q);
  BFT_FREE(gkb->m__v);
  BFT_FREE(gkb->v);

  BFT_FREE(gkb->zeta_array);

  BFT_FREE(gkb);
  *p_gkb = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a Uzawa builder structure
 *
 * \param[in]  gamma      value of the grad-div coefficient
 * \param[in]  n_u_dofs   number of velocity DoFs (degrees of freedom)
 * \param[in]  n_p_dofs   number of pressure DoFs
 * \param[in]  quant      pointer to additional mesh quantities
 *
 * \return a pointer to a new allocated Uzawa builder
 */
/*----------------------------------------------------------------------------*/

static cs_uza_builder_t *
_init_uzawa_builder(cs_real_t                     gamma,
                    cs_lnum_t                     n_u_dofs,
                    cs_lnum_t                     n_p_dofs,
                    const cs_cdo_quantities_t    *quant)
{
  cs_uza_builder_t  *uza = NULL;

  BFT_MALLOC(uza, 1, cs_uza_builder_t);

  uza->gamma = gamma;
  uza->n_u_dofs = n_u_dofs;
  uza->n_p_dofs = n_p_dofs;

  BFT_MALLOC(uza->b_tilda, n_u_dofs, cs_real_t);

  /* Auxiliary vectors */
  BFT_MALLOC(uza->inv_mp, n_p_dofs, cs_real_t);
# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    uza->inv_mp[ip] = 1./quant->cell_vol[ip];

  BFT_MALLOC(uza->res_p, n_p_dofs, cs_real_t);
  BFT_MALLOC(uza->d__v, n_p_dofs, cs_real_t);
  BFT_MALLOC(uza->rhs, n_u_dofs, cs_real_t);

  cs_navsto_algo_info_init(&(uza->info));

  return uza;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a Uzawa builder structure
 *
 * \param[in, out]  p_uza   double pointer to a Uzawa builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_free_uza_builder(cs_uza_builder_t   **p_uza)
{
  cs_uza_builder_t  *uza = *p_uza;

  if (uza == NULL)
    return;

  BFT_FREE(uza->b_tilda);

  BFT_FREE(uza->inv_mp);
  BFT_FREE(uza->res_p);
  BFT_FREE(uza->d__v);
  BFT_FREE(uza->rhs);

  BFT_FREE(uza);
  *p_uza = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the divergence operator and store the result in div_v
 *
 * \param[in]      div_op  pointer to the values of divergence operator
 * \param[in]      v       vector to apply in velocity space
 * \param[in, out] div_v   resulting vector in pressure space
 */
/*----------------------------------------------------------------------------*/

static void
_apply_div_op(const cs_real_t   *div_op,
              const cs_real_t   *v,
              cs_real_t         *div_v)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    cs_real_t _div_v = 0;
    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++)
      _div_v += cs_math_3_dot_product(div_op + 3*j, v + 3*c2f->ids[j]);
    div_v[c_id] = _div_v;

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the gradient operator (which is the transpose of the
 *         divergence operator) and store the result in dt_q
 *
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in]      q        vector to apply in pressure space
 * \param[in, out] dt_q     resulting vector in velocity space
 */
/*----------------------------------------------------------------------------*/

static void
_apply_div_op_transpose(const cs_real_t   *div_op,
                        const cs_real_t   *q,
                        cs_real_t         *dt_q)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_adjacency_t  *c2f = cs_shared_connect->c2f;

  memset(dt_q, 0, 3*quant->n_faces*sizeof(cs_real_t));

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const cs_real_t  qc = q[c_id];
    for (cs_lnum_t j = c2f->idx[c_id]; j < c2f->idx[c_id+1]; j++) {

      const cs_real_t  *_div_f = div_op + 3*j;

      cs_real_t  *_dt_q = dt_q + 3*c2f->ids[j];
#     pragma omp critical
      {
        _dt_q[0] += qc * _div_f[0];
        _dt_q[1] += qc * _div_f[1];
        _dt_q[2] += qc * _div_f[2];
      }

    } /* Loop on cell faces */

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Transform the initial saddle-point problem. The velocity unknown
 *         is modified and is stored in u_tilda as well as the RHS related to
 *         the mass equation and stored in b_tilda
 *
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in, out] gkb      pointer to a GKB builder structure
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in, out] u_f      initial velocity on faces
 * \param[in, out] b_f      right-hand side (scatter/gather if needed) on faces
 * \param[in, out] b_c      right_hand side on cells (mass equation)
 */
/*----------------------------------------------------------------------------*/

static void
_transform_gkb_system(const cs_matrix_t             *matrix,
                      const cs_equation_param_t     *eqp,
                      const cs_real_t               *div_op,
                      cs_gkb_builder_t              *gkb,
                      cs_sles_t                     *sles,
                      cs_real_t                     *u_f,
                      cs_real_t                     *b_f,
                      cs_real_t                     *b_c)
{
  assert(gkb != NULL);

  cs_real_t  normalization = 1.0; /* TODO */

  /* Modifiy the tolerance in order to be more accurate on this step */
  cs_equation_param_t  *_eqp = NULL;

  BFT_MALLOC(_eqp, 1, cs_equation_param_t);
  BFT_MALLOC(_eqp->name, strlen(eqp->name) + strlen(":gkb0") + 1, char);
  sprintf(_eqp->name, "%s:gkb0", eqp->name);
  _eqp->sles_param.field_id = eqp->sles_param.field_id;

  cs_equation_param_update_from(eqp, _eqp);
  _eqp->sles_param.eps = fmin(0.1*eqp->sles_param.eps, 1e-10);

  bool  rhs_redux = true;
  if (gkb->gamma > 0) {

    rhs_redux = false;

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++)
      gkb->b_tilda[ip] = gkb->gamma*b_c[ip]/cs_shared_quant->cell_vol[ip];

    /* Solve Dt.b_tilda */
    _apply_div_op_transpose(div_op, gkb->b_tilda, gkb->dt_q);

#   pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++)
      gkb->b_tilda[iu] = b_f[iu] + gkb->dt_q[iu];

    if (cs_glob_n_ranks > 1)
      cs_interface_set_sum(cs_shared_range_set->ifs,
                           gkb->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           gkb->b_tilda);

  }
  else
    memcpy(gkb->b_tilda, b_f, gkb->n_u_dofs*sizeof(cs_real_t));

  /* Compute M^-1.(b_f + gamma. Bt.N^-1.b_c) up to now gamma = 0 */
  gkb->info.n_inner_iter
    += (gkb->info.last_inner_iter
        = cs_equation_solve_scalar_system(gkb->n_u_dofs,
                                          _eqp,
                                          matrix,
                                          cs_shared_range_set,
                                          normalization,
                                          rhs_redux,
                                          sles,
                                          gkb->v,
                                          gkb->b_tilda));

  /* Compute the initial u_tilda := u_f - M^-1.b_f */
# pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++)
    gkb->u_tilda[iu] = u_f[iu] - gkb->v[iu];

  /* Compute b_tilda := b_c - div(M^-1.b_f) */
  _apply_div_op(div_op, gkb->v, gkb->d__v);

# pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++)
    gkb->b_tilda[ip] = b_c[ip] - gkb->d__v[ip];

  cs_equation_free_param(_eqp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the GKB algorithm
 *
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      div_op   pointer to the values of divergence operator
 * \param[in, out] gkb      pointer to a GKB builder structure
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in, out] p_c      right_hand side on cells (mass equation)
 */
/*----------------------------------------------------------------------------*/

static void
_init_gkb_algo(const cs_matrix_t             *matrix,
               const cs_equation_param_t     *eqp,
               const cs_real_t               *div_op,
               cs_gkb_builder_t              *gkb,
               cs_sles_t                     *sles,
               cs_real_t                     *p_c)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  size = quant->n_cells;

  double beta2 = 0.0;

  /* Compute beta := ||b_tilta||_N^-1 and q := N^-1(b_tilda)/beta */
# pragma omp parallel reduction(+:beta2) if (size > CS_THR_MIN)
  {
    cs_lnum_t s_id, e_id;
    _thread_range(size, &s_id, &e_id);

    const cs_lnum_t  n = e_id - s_id;
    const cs_real_t  *_w = quant->cell_vol + s_id;
    const cs_real_t  *_b = gkb->b_tilda + s_id;
    const cs_lnum_t  block_size = CS_SBLOCK_BLOCK_SIZE;
    const cs_lnum_t  n_blocks = (n + block_size - 1) / block_size;
    const cs_lnum_t  n_sblocks = (n_blocks > 3) ? sqrt(n_blocks) : 1;
    const cs_lnum_t  blocks_in_sblocks =
      (n + block_size*n_sblocks - 1) / (block_size*n_sblocks);

    cs_real_t  *_q = gkb->q + s_id;
    cs_lnum_t  shift = 0;

    for (cs_lnum_t s = 0; s < n_sblocks; s++) { /* Loop on slices */

      double  s_beta2 = 0.0;

      for (cs_lnum_t b_id = 0; b_id < blocks_in_sblocks; b_id++) {

        const cs_lnum_t  start_id = shift;
        shift += block_size;
        if (shift > n)
          shift = n, b_id = blocks_in_sblocks;
        const cs_lnum_t  end_id = shift;

        double  _beta2 = 0.0;
        for (cs_lnum_t j = start_id; j < end_id; j++) {

          const  cs_real_t  b_ov_w = _b[j]/_w[j];
          _beta2 += b_ov_w*_b[j];
          _q[j] = b_ov_w;

        } /* Loop on block_size */

        s_beta2 += _beta2;

      } /* Loop on blocks */

      beta2 += s_beta2;

    } /* Loop on super-blocks */

  } /* OpenMP block */

  /* Parallel synchronization */
  cs_parall_sum(1, CS_DOUBLE, &beta2);

  /* Keep the value of beta = ||b||_{N^-1} */
  assert(beta2 > -DBL_MIN);
  gkb->beta = sqrt(beta2);

  /* Store M^-1.(b_f + gamma. Bt.N^-1.b_c) in b_tilda which is not useful
   * anymore */
  memcpy(gkb->b_tilda, gkb->v, gkb->n_u_dofs*sizeof(cs_real_t));

  if (fabs(gkb->beta) > FLT_MIN) {
    const  cs_real_t  inv_beta = 1./gkb->beta;
# pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < size; i++)
      gkb->q[i] *= inv_beta;
  }
  else {
    gkb->info.cvg = CS_SLES_CONVERGED;
    return;
  }

  /* Solve M.w = Dt.q */
  _apply_div_op_transpose(div_op, gkb->q, gkb->dt_q);

  if (cs_glob_n_ranks > 1)
    cs_interface_set_sum(cs_shared_range_set->ifs,
                         gkb->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         gkb->dt_q);

  cs_real_t  normalization = 1.0; /* TODO */

  gkb->info.n_inner_iter
    += (gkb->info.last_inner_iter =
        cs_equation_solve_scalar_system(gkb->n_u_dofs,
                                        eqp,
                                        matrix,
                                        cs_shared_range_set,
                                        normalization,
                                        false, /* rhs_redux */
                                        sles,
                                        gkb->v,
                                        gkb->dt_q));

  gkb->alpha = _face_gdot(gkb->n_u_dofs, gkb->v, gkb->dt_q);
  assert(gkb->alpha > -DBL_MIN);
  gkb->alpha = sqrt(gkb->alpha);

  const double ov_alpha = 1./gkb->alpha;

  gkb->zeta = gkb->beta * ov_alpha;

  /* Initialize auxiliary vectors and first update of the solution vectors */

# pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++) {
    gkb->v[iu] *= ov_alpha;
    gkb->u_tilda[iu] = gkb->zeta * gkb->v[iu];
    gkb->m__v[iu] = ov_alpha * gkb->dt_q[iu];
  }

# pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++) {
    gkb->d[ip] = gkb->q[ip] * ov_alpha;
    p_c[ip] = -gkb->zeta * gkb->d[ip];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more GKB iteration
 *
 * \param[in]      nsp     pointer to a set of parameters for Navier-Stokes
 * \param[in, out] gkb     pointer to a GKB builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_gkb_cvg_test(const cs_navsto_param_t    *nsp,
              cs_gkb_builder_t           *gkb)
{
  const cs_real_t  diverg_factor = 100;

  /* Update the sum of square of zeta values (used for renormalization) */
  cs_real_t  z2 = gkb->zeta*gkb->zeta;

  gkb->zeta_square_sum += z2;
  gkb->zeta_array[gkb->info.n_algo_iter % gkb->z_size] = z2;

  /* Increment the number of Picard iterations */
  gkb->info.n_algo_iter += 1;

  /* Compute the relative energy norm. The normalization arises from an
     iterative estimation of the initial error in the energy norm */
  const cs_real_t  prev_res = gkb->info.res;

  int  n = gkb->z_size;
  if (gkb->info.n_algo_iter < gkb->z_size)
    n = gkb->info.n_algo_iter;

  cs_real_t  err2_energy = 0.;
  for (int i = 0; i < n; i++)
    err2_energy += gkb->zeta_array[i];

  double  tau = (gkb->gamma > 0) ?
    gkb->gamma*nsp->residual_tolerance : nsp->residual_tolerance;

  gkb->info.res = sqrt(err2_energy);

  /* Set the convergence status */
#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                "\nGKB.It%02d-- err2 = %6.4e ?<? tau * square_sum %6.4e\n",
                gkb->info.n_algo_iter, err2_energy, tau * gkb->zeta_square_sum);
#endif

  if (err2_energy < tau * gkb->zeta_square_sum)
    gkb->info.cvg = CS_SLES_CONVERGED;
  else if (gkb->info.n_algo_iter >= nsp->max_algo_iter)
    gkb->info.cvg = CS_SLES_MAX_ITERATION;
  else if (gkb->info.res > diverg_factor * prev_res)
    gkb->info.cvg = CS_SLES_DIVERGED;
  else
    gkb->info.cvg = CS_SLES_ITERATING;

  if (nsp->verbosity > 2)
    cs_log_printf(CS_LOG_DEFAULT,
                  "GKB.It%02d-- %5.3e %5d %6d z2:%6.4e renorm:%6.4e cvg:%d\n",
                  gkb->info.n_algo_iter, gkb->info.res,
                  gkb->info.last_inner_iter, gkb->info.n_inner_iter,
                  z2, sqrt(gkb->zeta_square_sum), gkb->info.cvg);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one needs one more Uzawa iteration
 *
 * \param[in]      nsp     pointer to a set of parameters for Navier-Stokes
 * \param[in, out] uza     pointer to a Uzawa builder structure
 */
/*----------------------------------------------------------------------------*/

static void
_uza_cvg_test(const cs_navsto_param_t    *nsp,
              cs_uza_builder_t           *uza)
{
  const cs_real_t  diverg_factor = 100;

  /* Increment the number of Picard iterations */
  uza->info.n_algo_iter += 1;

  /* Compute the new residual based on the norm of the divergence constraint */
  const cs_real_t  prev_res = uza->info.res;

  cs_real_t  res_square = cs_dot_wxx(uza->n_p_dofs, uza->inv_mp, uza->res_p);
  cs_parall_sum(1, CS_DOUBLE, &res_square);
  assert(res_square > -DBL_MIN);
  uza->info.res = sqrt(res_square);

  double  tau = (uza->gamma > 0) ?
    nsp->residual_tolerance/uza->gamma : nsp->residual_tolerance;

  /* Set the convergence status */
#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 0
  cs_log_printf(CS_LOG_DEFAULT,
                "\nUZA.It%02d-- res = %6.4e ?<? eps %6.4e\n",
                uza->info.n_algo_iter, uza->info.res, nsp->residual_tolerance);
#endif

  if (uza->info.res < tau)
    uza->info.cvg = CS_SLES_CONVERGED;
  else if (uza->info.n_algo_iter >= nsp->max_algo_iter)
    uza->info.cvg = CS_SLES_MAX_ITERATION;
  else if (uza->info.res > diverg_factor * prev_res)
    uza->info.cvg = CS_SLES_DIVERGED;
  else
    uza->info.cvg = CS_SLES_ITERATING;

  if (nsp->verbosity > 2)
    cs_log_printf(CS_LOG_DEFAULT,
                  "UZA.It%02d-- %5.3e %5d %6d cvg:%d\n",
                  uza->info.n_algo_iter, uza->info.res,
                  uza->info.last_inner_iter, uza->info.n_inner_iter,
                  uza->info.cvg);

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create an empty cs_cdofb_monolithic_sles_t structure
 *
 * \return a pointer to a newly allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdofb_monolithic_sles_t *
cs_cdofb_monolithic_sles_create(void)
{
  cs_cdofb_monolithic_sles_t  *msles = NULL;

  BFT_MALLOC(msles, 1, cs_cdofb_monolithic_sles_t);

  msles->matrix = NULL;
  msles->div_op = NULL;

  msles->graddiv_coef = 0.;

  msles->sles = NULL;

  msles->u_f = NULL;
  msles->p_c = NULL;
  msles->b_f = NULL;
  msles->b_c = NULL;

  return msles;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create an empty cs_cdofb_monolithic_sles_t structure
 *
 * \param[in, out]  p_msles  double pointer to the structure to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_free(cs_cdofb_monolithic_sles_t   **p_msles)
{
  cs_cdofb_monolithic_sles_t  *msles = *p_msles;

  if (msles == NULL)
    return;

  BFT_FREE(msles->div_op);
  /* other pointer are shared, thus no free at this stage */

  BFT_FREE(msles);
  *p_msles = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set pointers to shared structures
 *
 * \param[in]  connect  pointer to cdo connectivities
 * \param[in]  quant    pointer to additional mesh quantities
 * \param[in]  rset     pointer to a \ref cs_range_set_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_sles_set_shared(const cs_cdo_connect_t        *connect,
                                    const cs_cdo_quantities_t     *quant,
                                    const cs_range_set_t          *rset)
{
  assert(rset != NULL);

  /* Assign static const pointers */
  cs_shared_connect = connect;
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
  cs_param_sles_t  *mom_slesp = &(mom_eqp->sles_param);
  int  field_id = cs_equation_get_field_id(nsc->momentum);

  mom_slesp->field_id = field_id;
  if (mom_slesp->amg_type == CS_PARAM_AMG_NONE)
    mom_slesp->amg_type = CS_PARAM_AMG_HYPRE_BOOMER;

  /* Initialization must be called before setting options;
     it does not need to be called before calling
     cs_sles_petsc_define(), as this is handled automatically. */

  switch (nsp->sles_strategy) {

  case CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK: /* "Classical" way to set SLES */
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_GKB_SATURNE:
     /* Set solver and preconditioner for solving M = A + zeta * Bt*N^-1*B
      * Notice that zeta can be equal to 0 */
    cs_equation_param_set_sles(mom_eqp);
    break;

  case CS_NAVSTO_SLES_UZAWA_AL:
     /* Set solver and preconditioner for solving M = A + zeta * Bt*N^-1*B
      * Notice that zeta can be equal to 0 */
    cs_equation_param_set_sles(mom_eqp);
    break;

#if defined(HAVE_PETSC)
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
#else  /* HAVE_PETSC */
  case CS_NAVSTO_SLES_GKB:
  case CS_NAVSTO_SLES_GKB_GMRES:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please use a version of Code_Saturne built with PETSc.",
              __func__, mom_eqp->name);
    break;

#endif

#if defined(HAVE_PETSC)
  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _additive_amg_gmres_hook,
                         (void *)mom_eqp);
    break;

  case CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _multiplicative_gmres_hook,
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

  case CS_NAVSTO_SLES_MUMPS:
#if defined(PETSC_HAVE_MUMPS)
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _mumps_hook,
                         (void *)mom_eqp);
#else
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc with MUMPS is required with this option.\n",
              __func__, mom_eqp->name);
#endif
    break;

#else
  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
  case CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK:
  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
  case CS_NAVSTO_SLES_MUMPS:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please use a version of Code_Saturne built with PETSc.",
              __func__, mom_eqp->name);
    break;
#endif /* HAVE_PETSC */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n",
              __func__, mom_eqp->name);
  }

  /* Define the level of verbosity for SLES structure */
  if (mom_slesp->verbosity > 1) {

    cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

    /* Set verbosity */
    cs_sles_set_verbosity(sles, mom_slesp->verbosity);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from a scalar-valued CDO-Fb scheme
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the (cumulated) number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_solve(const cs_navsto_param_t       *nsp,
                          const cs_equation_param_t     *eqp,
                          cs_cdofb_monolithic_sles_t    *msles)
{
  CS_UNUSED(nsp);

  const cs_matrix_t  *matrix = msles->matrix;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_lnum_t  n_scatter_elts = 3*n_faces + n_cells;

  /* De-interlace the velocity array and the rhs for the face DoFs */
  cs_real_t  *xsol = NULL;
  BFT_MALLOC(xsol, n_cols, cs_real_t);

  cs_real_t  *b = NULL;
  BFT_MALLOC(b, n_scatter_elts, cs_real_t);

# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(msles, xsol, b)                             \
  firstprivate(n_faces)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    xsol[f            ] = msles->u_f[3*f];
    xsol[f +   n_faces] = msles->u_f[3*f+1];
    xsol[f + 2*n_faces] = msles->u_f[3*f+2];

    b[f            ] = msles->b_f[3*f];
    b[f +   n_faces] = msles->b_f[3*f+1];
    b[f + 2*n_faces] = msles->b_f[3*f+2];

  }

  /* Add the pressure related elements */
  memcpy(xsol + 3*n_faces, msles->p_c, n_cells*sizeof(cs_real_t));
  memcpy(b + 3*n_faces, msles->b_c, n_cells*sizeof(cs_real_t));

  const cs_range_set_t  *rset = cs_shared_range_set;
  int  n_iters = 0;
  double  residual = DBL_MAX;

  /* Prepare solving (handle parallelism) */
  cs_gnum_t  nnz = cs_equation_prepare_system(1,     /* stride */
                                              n_scatter_elts,
                                              matrix,
                                              rset,
                                              true,  /* rhs_redux */
                                              xsol, b);

  /* Solve the linear solver */
  const double  r_norm = 1.0; /* No renormalization by default (TODO) */
  const cs_param_sles_t  sles_param = eqp->sles_param;

  cs_sles_convergence_state_t  code = cs_sles_solve(msles->sles,
                                                    matrix,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    sles_param.eps,
                                                    r_norm,
                                                    &n_iters,
                                                    &residual,
                                                    b,
                                                    xsol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  if (sles_param.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%s/sles_cvg> code %-d n_iters %d"
                  " residual % -8.4e nnz %lu\n",
                  eqp->name, code, n_iters, residual, nnz);

  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         xsol, xsol);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_SLES_DBG > 1
  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         b, b);

  cs_dbg_fprintf_system(eqp->name, cs_shared_time_step->nt_cur,
                        CS_CDOFB_MONOLITHIC_DBG,
                        xsol, b, 3*n_faces);
#endif

  /* Interlace xsol --> u_f and p_c */
# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(msles, xsol) firstprivate(n_faces)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    msles->u_f[3*f]   = xsol[f];
    msles->u_f[3*f+1] = xsol[f +   n_faces];
    msles->u_f[3*f+2] = xsol[f + 2*n_faces];

  }

  /* Copy the part of the solution array related to the pressure in cells */
  memcpy(msles->p_c, xsol + 3*n_faces, n_cells*sizeof(cs_real_t));

  /* Free what can be freed at this stage */
  BFT_FREE(xsol);
  BFT_FREE(b);

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the GKB algorithm to solve the saddle-point problem arising
 *         from CDO-Fb schemes for Stokes and Navier-Stokes with a monolithic
 *         coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_gkb_solve(const cs_navsto_param_t       *nsp,
                              const cs_equation_param_t     *eqp,
                              cs_cdofb_monolithic_sles_t    *msles)
{
  /* Sanity checks */
  assert(nsp != NULL && nsp->sles_strategy == CS_NAVSTO_SLES_GKB_SATURNE);
  assert(cs_shared_range_set != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  *vol = quant->cell_vol;
  const cs_real_t  gamma = msles->graddiv_coef;
  const cs_real_t  *div_op = msles->div_op;
  const cs_matrix_t  *matrix = msles->matrix;

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = msles->b_f;
  cs_real_t  *b_c = msles->b_c;

  /* Allocate and initialize the GKB builder structure */
  cs_gkb_builder_t  *gkb = _init_gkb_builder(gamma,
                                             3*quant->n_faces,
                                             quant->n_cells);

  /* Transformation of the initial saddle-point system */
  _transform_gkb_system(matrix, eqp, div_op, gkb, msles->sles, u_f, b_f, b_c);

  /* Initialization */
  _init_gkb_algo(matrix, eqp, div_op, gkb, msles->sles, p_c);

  /* Main loop */
  /* ========= */

  while (gkb->info.cvg == CS_SLES_ITERATING) {

    /* Compute g (store as an update of d__v), q */
    _apply_div_op(div_op, gkb->v, gkb->d__v);

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++) {
      gkb->d__v[ip] /= vol[ip];
      gkb->d__v[ip] -= gkb->alpha * gkb->q[ip];
    }

    /* Compute beta */
    gkb->beta = cs_dot_wxx(gkb->n_p_dofs, vol, gkb->d__v);
    cs_parall_sum(1, CS_DOUBLE, &(gkb->beta));
    assert(gkb->beta > -DBL_MIN);
    gkb->beta = sqrt(gkb->beta);

    const double  ov_beta = 1./gkb->beta;

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++)
      gkb->q[ip] = ov_beta*gkb->d__v[ip];

    /* Solve M.w_tilda = Dt.q */
    _apply_div_op_transpose(div_op, gkb->q, gkb->dt_q);

    if (cs_glob_n_ranks > 1)
      cs_interface_set_sum(cs_shared_range_set->ifs,
                           gkb->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           gkb->dt_q);

    /* Prepare update of m__v:
     *  m__v(k+1) = 1/alpha(k+1) * (dt_q - beta*m__v(k)) */
#   pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++) {
      gkb->m__v[iu] *= -gkb->beta;
      gkb->m__v[iu] +=  gkb->dt_q[iu];
    }

    cs_real_t  normalization = gkb->alpha; /* TODO */
    gkb->info.n_inner_iter
      += (gkb->info.last_inner_iter =
          cs_equation_solve_scalar_system(gkb->n_u_dofs,
                                          eqp,
                                          matrix,
                                          cs_shared_range_set,
                                          normalization,
                                          false, /* rhs_redux */
                                          msles->sles,
                                          gkb->v,
                                          gkb->m__v));

    /* Compute alpha */
    gkb->alpha = _face_gdot(gkb->n_u_dofs, gkb->v, gkb->m__v);
    assert(gkb->alpha > -DBL_MIN);
    gkb->alpha = sqrt(gkb->alpha);

    const double ov_alpha = 1./gkb->alpha;

    /* zeta(k+1) = -beta/alpha * zeta(k) */
    gkb->zeta *= -gkb->beta * ov_alpha;

    /* Update vectors and solutions */
#   pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++) {
      gkb->v[iu] *= ov_alpha;
      gkb->u_tilda[iu] += gkb->zeta * gkb->v[iu];
      /* Last step: m__v(k+1) = 1/alpha(k+1) * (dt_q - beta*m__v(k)) */
      gkb->m__v[iu] *= ov_alpha;
    }

#   pragma omp parallel for if (gkb->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < gkb->n_p_dofs; ip++) {
      gkb->d[ip] = ov_alpha * (gkb->q[ip] - gkb->beta*gkb->d[ip]);
      p_c[ip] += -gkb->zeta * gkb->d[ip];
    }

    /* Update error norm and test if one needs one more iteration */
    _gkb_cvg_test(nsp, gkb);

  }

  /* Return to the initial velocity formulation
   * u: = u_tilda + M^-1.(b_f + gamma.N^-1.b_c)
   * where M^-1.(b_f + gamma.N^-1.b_c) is stored in b_tilda */
# pragma omp parallel for if (gkb->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < gkb->n_u_dofs; iu++)
    u_f[iu] = gkb->u_tilda[iu] + gkb->b_tilda[iu];

  int n_inner_iter = gkb->info.n_inner_iter;

  /* Last step: Free temporary memory */
  _free_gkb_builder(&gkb);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Use the Uzawa algorithm with an Augmented Lagrangian technique to
 *         solve the saddle-point problem arisinfg from CDO-Fb schemes for
 *         Stokes and Navier-Stokes with a monolithic coupling
 *
 * \param[in]      nsp      pointer to a cs_navsto_param_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] msles    pointer to a cs_cdofb_monolithic_sles_t structure
 *
 * \return the cumulated number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

int
cs_cdofb_monolithic_uzawa_al_solve(const cs_navsto_param_t       *nsp,
                                   const cs_equation_param_t     *eqp,
                                   cs_cdofb_monolithic_sles_t    *msles)
{
  /* Sanity checks */
  assert(nsp != NULL && nsp->sles_strategy == CS_NAVSTO_SLES_UZAWA_AL);
  assert(cs_shared_range_set != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  gamma = msles->graddiv_coef;
  const cs_real_t  *div_op = msles->div_op;

  cs_real_t  *u_f = msles->u_f;
  cs_real_t  *p_c = msles->p_c;
  cs_real_t  *b_f = msles->b_f;
  cs_real_t  *b_c = msles->b_c;

  /* Allocate and initialize the GKB builder structure */
  cs_uza_builder_t  *uza = _init_uzawa_builder(gamma,
                                               3*quant->n_faces,
                                               quant->n_cells,
                                               quant);

  /* Transformation of the initial right-hand side */
  cs_real_t  *btilda_c = uza->d__v;
# pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
  for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++)
    btilda_c[ip] = uza->inv_mp[ip]*b_c[ip];

  _apply_div_op_transpose(div_op, btilda_c, uza->b_tilda);

  if (cs_glob_n_ranks > 1)
    cs_interface_set_sum(cs_shared_range_set->ifs,
                         uza->n_u_dofs,
                         1, false, CS_REAL_TYPE, /* stride, interlaced */
                         uza->b_tilda);

  /* Update the modify right-hand side: b_tilda = b_f + gamma*Dt.W^-1.b_c */
# pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
  for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++) {
    uza->b_tilda[iu] *= gamma;
    uza->b_tilda[iu] += b_f[iu];
  }

  /* Main loop */
  /* ========= */

  while (uza->info.cvg == CS_SLES_ITERATING) {

    /* Compute the RHS for the Uzawa system: rhs = b_tilda - Dt.p_c */
    _apply_div_op_transpose(div_op, p_c, uza->rhs);

    if (cs_glob_n_ranks > 1)
      cs_interface_set_sum(cs_shared_range_set->ifs,
                           uza->n_u_dofs,
                           1, false, CS_REAL_TYPE, /* stride, interlaced */
                           uza->rhs);

#   pragma omp parallel for if (uza->n_u_dofs > CS_THR_MIN)
    for (cs_lnum_t iu = 0; iu < uza->n_u_dofs; iu++) {
      uza->rhs[iu] *= -1;
      uza->rhs[iu] += uza->b_tilda[iu];
    }

    /* Solve AL.u_f = rhs */
    cs_real_t  normalization = 1.; /* TODO */
    uza->info.n_inner_iter
      += (uza->info.last_inner_iter =
          cs_equation_solve_scalar_system(uza->n_u_dofs,
                                          eqp,
                                          msles->matrix,
                                          cs_shared_range_set,
                                          normalization,
                                          false, /* rhs_redux */
                                          msles->sles,
                                          u_f,
                                          uza->rhs));

    /* Update p_c = p_c - gamma * (D.u_f - b_c) */
    _apply_div_op(div_op, u_f, uza->d__v);

#   pragma omp parallel for if (uza->n_p_dofs > CS_THR_MIN)
    for (cs_lnum_t ip = 0; ip < uza->n_p_dofs; ip++) {
      uza->res_p[ip] = uza->d__v[ip] - b_c[ip];
      p_c[ip] += gamma * uza->inv_mp[ip] * uza->res_p[ip];
    }

    /* Update error norm and test if one needs one more iteration */
    _uza_cvg_test(nsp, uza);

  }

  int n_inner_iter = uza->info.n_inner_iter;

  /* Last step: Free temporary memory */
  _free_uza_builder(&uza);

  return  n_inner_iter;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
