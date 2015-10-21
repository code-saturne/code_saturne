/*============================================================================
 * Routines to handle cs_equation_t structure and its related structures
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include <string.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_PETSC)
#include <petscdraw.h>
#include <petscviewer.h>
#include <petscksp.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_base.h"
#include "cs_cdo.h"
#include "cs_mesh_location.h"
#include "cs_field.h"
#include "cs_post.h"
#include "cs_multigrid.h"
#include "cs_timer_stats.h"
#include "cs_param.h"
#include "cs_cdovb_codits.h"
#include "cs_cdofb_codits.h"
#include "cs_cdo_convection.h"

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_equation.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a builder structure
 *
 * \param[in] eq       pointer to a cs_equation_param_t structure
 * \param[in] m        pointer to a mesh structure
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_init_builder_t)(const cs_equation_param_t  *eqp,
                             const cs_mesh_t            *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of source terms for the current time
 *
 * \param[in]      m          pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  time_step structure
 * \param[in, out] builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_compute_source_t)(const cs_mesh_t            *mesh,
                               const cs_cdo_connect_t     *connect,
                               const cs_cdo_quantities_t  *cdoq,
                               const cs_time_step_t       *time_step,
                               void                       *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system within the CDO framework
 *
 * \param[in]      m          pointer to a cs_mesh_t structure
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in]      quant      pointer to a cs_cdo_quantities_t structure
 * \param[in]      time_step  pointer to a time_step structure
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in, out] builder    pointer to builder structure
 * \param[in, out] rhs        pointer to a right-hand side array pointer
 * \param[in, out] sla_mat    pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_build_system_t)(const cs_mesh_t            *mesh,
                             const cs_cdo_connect_t     *connect,
                             const cs_cdo_quantities_t  *cdoq,
                             const cs_time_step_t       *time_step,
                             double                      dt_cur,
                             const cs_real_t            *field_val,
                             void                       *builder,
                             cs_real_t                 **rhs,
                             cs_sla_matrix_t           **sla_mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant      pointer to a cs_cdo_quantities_t struct.
 * \param[in]      time_step  pointer to a time step structure
 * \param[in]      solu       solution array
 * \param[in, out] builder    pointer to builder structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_update_field_t)(const cs_cdo_connect_t     *connect,
                             const cs_cdo_quantities_t  *cdoq,
                             const cs_time_step_t       *time_step,
                             const cs_real_t            *solu,
                             void                       *builder,
                             cs_real_t                  *field_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at each face
 *
 * \param[in]  builder    pointer to a cs_cdofb_codits_t structure
 * \param[in]  field      pointer to a cs_field_t structure
 *
 * \return  a pointer to an array of double (size n_faces)
 */
/*----------------------------------------------------------------------------*/

typedef const double *
(cs_equation_get_f_values_t)(const void          *builder,
                             const cs_field_t    *field);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a builder structure
 *
 * \param[in, out]  builder   pointer to a builder structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_free_builder_t)(void  *builder);

/*============================================================================
 * Local variables
 *============================================================================*/

/* Default initialization */
static cs_equation_algo_t _algo_info_by_default = {
#if defined(HAVE_PETSC)
  CS_EQUATION_ALGO_PETSC_ITSOL, // Family of iterative solvers
#else
  CS_EQUATION_ALGO_CS_ITSOL,    // Family of iterative solvers
#endif
  0,                            // n_iters
  50,                           // max. number of iterations
  0,                            // n_cumulated_iters
  10000,                        // max. number of cumulated iterations
  1e-6                          // stopping criterion
};

static cs_param_itsol_t _itsol_info_by_default = {
#if defined(HAVE_PETSC)
  CS_PARAM_PRECOND_ILU0,  // preconditioner
  CS_PARAM_ITSOL_BICG,    // iterative solver
#else
  CS_PARAM_PRECOND_DIAG,  // preconditioner
  CS_PARAM_ITSOL_CG,      // iterative solver
#endif
  2500,                   // max. number of iterations
  1e-12,                  // stopping criterion on the accuracy
  150,                    // output frequency
  false                   // normalization of the residual (true or false)
};

/* List of available keys for setting an equation */
typedef enum {

  EQKEY_HODGE_DIFF_ALGO,
  EQKEY_HODGE_DIFF_COEF,
  EQKEY_HODGE_TIME_ALGO,
  EQKEY_HODGE_TIME_COEF,
  EQKEY_ITSOL,
  EQKEY_ITSOL_EPS,
  EQKEY_ITSOL_MAX_ITER,
  EQKEY_ITSOL_RESNORM,
  EQKEY_PRECOND,
  EQKEY_SOLVER_FAMILY,
  EQKEY_SPACE_SCHEME,
  EQKEY_VERBOSITY,
  EQKEY_BC_ENFORCEMENT,
  EQKEY_BC_QUADRATURE,
  EQKEY_OUTPUT_FREQ,
  EQKEY_POST_FREQ,
  EQKEY_POST,
  EQKEY_ADV_OP_TYPE,
  EQKEY_ADV_WEIGHT_ALGO,
  EQKEY_ADV_WEIGHT_CRIT,
  EQKEY_ADV_FLUX_QUADRA,
  EQKEY_TIME_SCHEME,
  EQKEY_TIME_THETA,
  EQKEY_ERROR

} eqkey_t;

/* List of keys for setting a source term */
typedef enum {

  STKEY_POST,
  STKEY_QUADRATURE,
  STKEY_ERROR

} stkey_t;

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

struct _cs_equation_t {

  char *restrict         name;    /* Short description */

  cs_equation_param_t   *param;   /* Set of parameters related to an equation */

  /* Variable atached to this equation is also attached to a cs_field_t
     structure */
  char *restrict         varname;
  int                    field_id;

  /* Timer statistic for a light profiling */
  int     main_ts_id;     /* Id of the main timer states structure related
                             to this equation */
  int     pre_ts_id;      /* Id of the timer stats structure gathering all
                             steps before the resolution of the linear
                             systems */
  int     solve_ts_id;    /* Id of the timer stats structure related
                             to the inversion of the linear system */
  int     post_ts_id;     /* Id of the timer stats structure gathering all
                             steps afterthe resolution of the linear systems
                             (post, balance...) */

  bool    do_build;       /* false => keep the system as it is */

  /* Algebraic system */
  cs_matrix_structure_t    *ms;      /* matrix structure (how are stored
                                        coefficients of the matrix a) */
  cs_matrix_t              *matrix;  // matrix to inverse with cs_sles_solve()
  cs_real_t                *rhs;     // right-hand side

  /* System builder depending on the numerical scheme*/
  void                     *builder;

  /* Pointer to functions */
  cs_equation_init_builder_t    *init_builder;
  cs_equation_compute_source_t  *compute_source;
  cs_equation_build_system_t    *build_system;
  cs_equation_update_field_t    *update_field;
  cs_equation_get_f_values_t    *get_f_values;
  cs_equation_free_builder_t    *free_builder;

};

/*============================================================================
 * Private function prototypes
 *============================================================================*/


#if defined(HAVE_PETSC)

/*----------------------------------------------------------------------------
 * \brief Add visualization of the matrix graph
 *
 * \param[in]  ksp     Krylov SubSpace structure
 *----------------------------------------------------------------------------*/

static void
_add_view(KSP          ksp)
{
  const char *p = getenv("CS_USER_PETSC_MAT_VIEW");

  if (p != NULL) {

    /* Get system and preconditioner matrixes */

    Mat a, pa;
    KSPGetOperators(ksp, &a, &pa);

    /* Output matrix in several ways depending on
       CS_USER_PETSC_MAT_VIEW environment variable */

    if (strcmp(p, "DEFAULT") == 0)
      MatView(a, PETSC_VIEWER_DEFAULT);

    else if (strcmp(p, "DRAW_WORLD") == 0)
      MatView(a, PETSC_VIEWER_DRAW_WORLD);

    else if (strcmp(p, "DRAW") == 0) {

      PetscViewer viewer;
      PetscDraw draw;
      PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, "PETSc View",
                          0, 0, 600, 600, &viewer);
      PetscViewerDrawGetDraw(viewer, 0, &draw);
      PetscViewerDrawSetPause(viewer, -1);
      MatView(a, viewer);
      PetscDrawPause(draw);

      PetscViewerDestroy(&viewer);

    }

  }

}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using CG with Jacobi preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_cg_diag_setup_hook(void   *context,
                    KSP     ksp)
{
  PC pc;

  KSPSetType(ksp, KSPCG);   /* Preconditioned Conjugate Gradient */

  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCJACOBI);  /* Jacobi (diagonal) preconditioning */

  _add_view(ksp);
}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using CG with SSOR preconditioner
 *        Warning: this PETSc implementation is only available in serial mode
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_cg_ssor_setup_hook(void   *context,
                    KSP     ksp)
{
  PC pc;

  KSPSetType(ksp, KSPCG);   /* Preconditioned Conjugate Gradient */

  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCSOR);  /* SSOR preconditioning */
  PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);

  _add_view(ksp);
}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using CG with Additive Schwarz preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_cg_as_setup_hook(void   *context,
                  KSP     ksp)
{
  PC pc;

  KSPSetType(ksp, KSPCG);   /* Preconditioned Conjugate Gradient */

  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCASM);

  _add_view(ksp);
}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using CG with ICC preconditioner
 *        Warning: this PETSc implementation is only available in serial mode
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_cg_icc_setup_hook(void    *context,
                   KSP      ksp)
{
  PC pc;

  KSPSetType(ksp, KSPCG);   /* Preconditioned Conjugate Gradient */

  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCICC);
  PCFactorSetLevels(pc, 0);

  _add_view(ksp);
}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using CG with GAMG preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_cg_gamg_setup_hook(void    *context,
                    KSP      ksp)
{
  PC pc;

  KSPSetType(ksp, KSPCG);   /* Preconditioned Conjugate Gradient */

  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCGAMG);

  _add_view(ksp);
}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using CG with Boomer AMG preconditioner (Hypre library)
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_cg_bamg_setup_hook(void    *context,
                    KSP      ksp)
{
  PC pc;

  KSPSetType(ksp, KSPCG);   /* Preconditioned Conjugate Gradient */

  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCHYPRE);

  _add_view(ksp);
}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using GMRES with ILU0 preconditioner
 *        Warning: this PETSc implementation is only available in serial mode
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_gmres_ilu_setup_hook(void    *context,
                      KSP      ksp)
{
  PC pc;

  const int  n_max_restart = 30;

  KSPSetType(ksp, KSPGMRES);   /* Preconditioned GMRES */

  KSPGMRESSetRestart(ksp, n_max_restart);
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCILU);
  PCFactorSetLevels(pc, 0);

  _add_view(ksp);
}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using GMRES with block Jacobi preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_gmres_bjacobi_setup_hook(void    *context,
                          KSP      ksp)
{
  PC pc;

  const int  n_max_restart = 30;

  KSPSetType(ksp, KSPGMRES);   /* Preconditioned GMRES */

  KSPGMRESSetRestart(ksp, n_max_restart);
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCBJACOBI);

  _add_view(ksp);
}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using BiCGStab with ILU0 preconditioner
 *        Warning: this PETSc implementation is only available in serial mode
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_bicg_ilu_setup_hook(void    *context,
                     KSP      ksp)
{
  PC pc;

  KSPSetType(ksp, KSPBCGS);   /* Preconditioned BiCGStab */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCILU);
  PCFactorSetLevels(pc, 0);

  _add_view(ksp);
}

/*----------------------------------------------------------------------------
 * \brief PETSc solver using BiCGStab with block Jacobi preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_bicg_bjacobi_setup_hook(void    *context,
                         KSP      ksp)
{
  PC pc;

  KSPSetType(ksp, KSPBCGS);   /* Preconditioned BICGStab */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); /* Try to have "true" norm */

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCBJACOBI);

  _add_view(ksp);
}

#endif /* defined(HAVE_PETSC) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Given its name, get the id related to a cs_mesh_location_t structure
 *
 * \param[in]      ml_name    name of the location
 * \param[in, out] p_ml_id    pointer on the id of the related mesh location
 */
/*----------------------------------------------------------------------------*/

static void
_check_ml_name(const char   *ml_name,
               int          *p_ml_id)
{
  *p_ml_id = cs_mesh_location_get_id_by_name(ml_name);

  if (*p_ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid mesh location name %s.\n"
                " This mesh location is not already defined.\n"), ml_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the name of the corresponding equation key
 *
 * \param[in] key        name of the key
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

static const char *
_print_eqkey(eqkey_t  key)
{
  switch (key) {
  case EQKEY_HODGE_DIFF_ALGO:
    return "hodge_diff_algo";
  case EQKEY_HODGE_DIFF_COEF:
    return "hodge_diff_coef";
  case EQKEY_HODGE_TIME_ALGO:
    return "hodge_time_algo";
  case EQKEY_HODGE_TIME_COEF:
    return "hodge_time_coef";
  case EQKEY_ITSOL:
    return "itsol";
  case EQKEY_ITSOL_EPS:
    return "itsol_eps";
  case EQKEY_ITSOL_MAX_ITER:
    return "itsol_max_iter";
  case EQKEY_ITSOL_RESNORM:
    return "itsol_resnorm";
  case EQKEY_PRECOND:
    return "precond";
  case EQKEY_SOLVER_FAMILY:
    return "solver_family";
  case EQKEY_SPACE_SCHEME:
    return "space_scheme";
  case EQKEY_VERBOSITY:
    return "verbosity";
  case EQKEY_BC_ENFORCEMENT:
    return "bc_enforcement";
  case EQKEY_BC_QUADRATURE:
    return "bc_quadrature";
  case EQKEY_OUTPUT_FREQ:
    return "output_freq";
  case EQKEY_POST_FREQ:
    return "post_freq";
  case EQKEY_POST:
    return "post";
  case EQKEY_ADV_OP_TYPE:
    return "adv_formulation";
  case EQKEY_ADV_WEIGHT_ALGO:
    return "adv_weight";
  case EQKEY_ADV_WEIGHT_CRIT:
    return "adv_weight_criterion";
  case EQKEY_ADV_FLUX_QUADRA:
    return "adv_flux_quad";
  case EQKEY_TIME_SCHEME:
    return "time_scheme";
  case EQKEY_TIME_THETA:
    return "time_theta";

  default:
    assert(0);
  }

  return NULL; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the name of the corresponding source term key
 *
 * \param[in] key        name of the key
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

static const char *
_print_stkey(stkey_t  key)
{
  switch (key) {
  case STKEY_POST:
    return "post";
  case STKEY_QUADRATURE:
    return "quadrature";
  default:
    assert(0);
  }

  return NULL; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the corresponding enum from the name of an equation key.
 *         If not found, print an error message
 *
 * \param[in] keyname    name of the key
 *
 * \return a eqkey_t
 */
/*----------------------------------------------------------------------------*/

static eqkey_t
_get_eqkey(const char *keyname)
{
  eqkey_t  key = EQKEY_ERROR;

  if (strncmp(keyname, "hodge", 5) == 0) { /* key begins with hodge */
    if (strcmp(keyname, "hodge_diff_coef") == 0)
      key = EQKEY_HODGE_DIFF_COEF;
    else if (strcmp(keyname, "hodge_diff_algo") == 0)
      key = EQKEY_HODGE_DIFF_ALGO;
    else if (strcmp(keyname, "hodge_time_coef") == 0)
      key = EQKEY_HODGE_TIME_COEF;
    else if (strcmp(keyname, "hodge_time_algo") == 0)
      key = EQKEY_HODGE_TIME_ALGO;
  }

  else if (strncmp(keyname, "itsol", 5) == 0) { /* key begins with itsol */
    if (strcmp(keyname, "itsol") == 0)
      key = EQKEY_ITSOL;
    else if (strcmp(keyname, "itsol_eps") == 0)
      key = EQKEY_ITSOL_MAX_ITER;
    else if (strcmp(keyname, "itsol_max_iter") == 0)
      key = EQKEY_ITSOL_MAX_ITER;
    else if (strcmp(keyname, "itsol_resnorm") == 0)
      key = EQKEY_ITSOL_RESNORM;
  }

  else if (strcmp(keyname, "precond") == 0)
    key = EQKEY_PRECOND;

  else if (strcmp(keyname, "solver_family") == 0)
    key = EQKEY_SOLVER_FAMILY;

  else if (strcmp(keyname, "space_scheme") == 0)
    key = EQKEY_SPACE_SCHEME;

  else if (strcmp(keyname, "verbosity") == 0)
    key = EQKEY_VERBOSITY;

  else if (strncmp(keyname, "bc", 2) == 0) { /* key begins with bc */
    if (strcmp(keyname, "bc_enforcement") == 0)
      key = EQKEY_BC_ENFORCEMENT;
    else if (strcmp(keyname, "bc_quadrature") == 0)
      key = EQKEY_BC_QUADRATURE;
  }

  else if (strcmp(keyname, "post") == 0)
    key = EQKEY_POST;
  else if (strcmp(keyname, "post_freq") == 0)
    key = EQKEY_POST_FREQ;
  else if (strcmp(keyname, "output_freq") == 0)
    key = EQKEY_OUTPUT_FREQ;

  else if (strncmp(keyname, "adv_", 4) == 0) {
    if (strcmp(keyname, "adv_formulation") == 0)
      key = EQKEY_ADV_OP_TYPE;
    else if (strcmp(keyname, "adv_weight_criterion") == 0)
      key = EQKEY_ADV_WEIGHT_CRIT;
    else if (strcmp(keyname, "adv_weight") == 0)
      key = EQKEY_ADV_WEIGHT_ALGO;
    else if (strcmp(keyname, "adv_flux_quad") == 0)
      key = EQKEY_ADV_FLUX_QUADRA;
  }

  else if (strncmp(keyname, "time_", 5) == 0) {
    if (strcmp(keyname, "time_scheme") == 0)
      key = EQKEY_TIME_SCHEME;
    else if (strcmp(keyname, "time_theta") == 0)
      key = EQKEY_TIME_THETA;
  }

  return key;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the corresponding enum from the name of a source term key.
 *         If not found, print an error message
 *
 * \param[in] keyname    name of the key
 *
 * \return a stkey_t
 */
/*----------------------------------------------------------------------------*/

static stkey_t
_get_stkey(const char *keyname)
{
  stkey_t  key = STKEY_ERROR;

  if (strcmp(keyname, "post") == 0)
    key = STKEY_POST;
  else if (strcmp(keyname, "quadrature") == 0)
    key = STKEY_QUADRATURE;

  return key;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a cs_equation_param_t
 *
 * \param[in] status
 * \param[in] type             type of equation (scalar, vector, tensor...)
 * \param[in] is_steady        add an unsteady term or not
 * \param[in] do_convection    add a convection term
 * \param[in] do_diffusion     add a diffusion term
 * \param[in] default_bc       type of boundary condition set by default
 *
 * \return a pointer to a new allocated cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_equation_param_t *
_create_equation_param(cs_equation_status_t   status,
                       cs_equation_type_t     type,
                       bool                   is_steady,
                       bool                   do_convection,
                       bool                   do_diffusion,
                       cs_param_bc_type_t     default_bc)
{
  cs_param_var_type_t  var_type = CS_PARAM_N_VAR_TYPES;
  cs_equation_param_t  *eqp = NULL;

  BFT_MALLOC(eqp, 1, cs_equation_param_t);

  eqp->status = status;
  eqp->type = type;
  eqp->verbosity =  0;
  eqp->output_freq = 1;
  eqp->post_freq = 10;
  eqp->post_flag =  0;

  /* Build the equation flag */
  eqp->flag = 0;
  if (!is_steady)
    eqp->flag |= CS_EQUATION_UNSTEADY;
  if (do_convection)
    eqp->flag |= CS_EQUATION_CONVECTION;
  if (do_diffusion)
    eqp->flag |= CS_EQUATION_DIFFUSION;

  eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;

  /* Vertex-based schemes imply the two following discrete Hodge operators
     Default initialization is made in accordance with this choice */
  eqp->is_multiplied_by_rho = true;
  eqp->time_hodge.pty_id = 0;      // Unity (default property)
  eqp->time_hodge.inv_pty = false; // inverse property ?
  eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
  eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;

  eqp->diffusion_hodge.pty_id = 0;      // Unity (default property)
  eqp->diffusion_hodge.inv_pty = false; // inverse property ?
  eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
  eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
  eqp->diffusion_hodge.coef = 1./3.;

  eqp->advection.adv_id = -1; // No default value
  eqp->advection.form = CS_PARAM_ADVECTION_FORM_CONSERV;
  eqp->advection.weight_algo = CS_PARAM_ADVECTION_WEIGHT_ALGO_UPWIND;
  eqp->advection.weight_criterion = CS_PARAM_ADVECTION_WEIGHT_XEXC;
  eqp->advection.quad_type = CS_QUADRATURE_BARY;

  /* Boundary conditions structure.
     One assigns a boundary condition by default */
  eqp->bc = cs_param_bc_create(default_bc);

  /* Description of the time discretization (default values) */
  eqp->time_info.scheme = CS_TIME_SCHEME_IMPLICIT;
  eqp->time_info.theta = 1.0;
  eqp->time_info.do_lumping = false;

  /* Initial condition (zero value by default) */
  eqp->time_info.ic_def_type = CS_PARAM_DEF_BY_VALUE;

  switch (type) {
  case CS_EQUATION_TYPE_SCAL:
    var_type = CS_PARAM_VAR_SCAL;
    break;
  case CS_EQUATION_TYPE_VECT:
    var_type = CS_PARAM_VAR_VECT;
    break;
  case CS_EQUATION_TYPE_TENS:
    var_type = CS_PARAM_VAR_TENS;
    break;
  default:
    break;
  }

  cs_param_set_def(eqp->time_info.ic_def_type, var_type, NULL,
                   &(eqp->time_info.ic_def));

  /* Settings for driving the linear algebra */
  eqp->algo_info = _algo_info_by_default;
  eqp->itsol_info = _itsol_info_by_default;

  /* Source terms */
  eqp->n_source_terms = 0;
  eqp->source_terms = NULL;

  return eqp;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an equation
 *
 * \param[in] eqname           name of the equation
 * \param[in] varname          name of the variable associated to this equation
 * \param[in] status           user or predefined equations
 * \param[in] type             type of equation (scalar, vector, tensor...)
 * \param[in] is_steady        add an unsteady term or not
 * \param[in] do_convection    add a convection term
 * \param[in] do_diffusion     add a diffusion term
 * \param[in] default_bc       type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_create(const char            *eqname,
                   const char            *varname,
                   cs_equation_status_t   status,
                   cs_equation_type_t     type,
                   bool                   is_steady,
                   bool                   do_convection,
                   bool                   do_diffusion,
                   cs_param_bc_type_t     default_bc)
{
  int  len = strlen(eqname)+1;

  cs_equation_t  *eq = NULL;

  BFT_MALLOC(eq, 1, cs_equation_t);

  /* Sanity checks */
  if (varname == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" No variable name associated to an equation structure.\n"
                " Check your initialization."));

  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" No equation name associated to an equation structure.\n"
                " Check your initialization."));

  /* Store eqname */
  BFT_MALLOC(eq->name, len, char);
  strncpy(eq->name, eqname, len);

  /* Store varname */
  len = strlen(varname)+1;
  BFT_MALLOC(eq->varname, len, char);
  strncpy(eq->varname, varname, len);

  eq->field_id = -1;    // field is created in a second step

  /* Set timer statistic structure to a default value */
  eq->main_ts_id = eq->pre_ts_id = eq->solve_ts_id = eq->post_ts_id = -1;

  eq->do_build = true;  // Force the construction of the algebraic system

  eq->param = _create_equation_param(status,
                                     type,
                                     is_steady,
                                     do_convection,
                                     do_diffusion,
                                     default_bc);

  /* Algebraic system: allocated later */
  eq->ms = NULL;
  eq->matrix = NULL;
  eq->rhs = NULL;

  /* Pointer of function */
  eq->init_builder = NULL;
  eq->compute_source = NULL;
  eq->build_system = NULL;
  eq->update_field = NULL;
  eq->get_f_values = NULL;
  eq->free_builder = NULL;

  eq->builder = NULL;

  return  eq;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_equation_t structure
 *
 * \param[in, out] eq    pointer to a cs_equation_t structure
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_free(cs_equation_t  *eq)
{
  if (eq == NULL)
    return eq;

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

  BFT_FREE(eq->name);
  BFT_FREE(eq->varname);

  cs_equation_param_t  *eqp = eq->param;

  if (eqp->bc != NULL) { // Boundary conditions
    if (eqp->bc->n_defs > 0)
      BFT_FREE(eqp->bc->defs);
    BFT_FREE(eqp->bc);
    eqp->bc = NULL;
  }

  if (eqp->n_source_terms > 0) { // Source terms
    for (int i = 0; i< eqp->n_source_terms; i++)
      BFT_FREE(eqp->source_terms[i].name);
    BFT_FREE(eqp->source_terms);
  }

  BFT_FREE(eq->param);

  cs_matrix_structure_destroy(&(eq->ms));
  cs_matrix_destroy(&(eq->matrix));
  BFT_FREE(eq->rhs);

  eq->builder = eq->free_builder(eq->builder);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

  BFT_FREE(eq);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_summary(const cs_equation_t  *eq)
{
  if (eq == NULL)
    return;

  const cs_equation_param_t  *eqp = eq->param;

  bft_printf("\n%s", lsepline);
  bft_printf("\tSummary of settings for %s eq. (variable %s)\n",
             eq->name, eq->varname);
  bft_printf("%s", lsepline);

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB)
    bft_printf("  <%s/Space scheme>  CDO vertex-based\n", eq->name);
  else if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB)
    bft_printf("  <%s/Space scheme>  CDO face-based\n", eq->name);

  bool  unsteady = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;
  bool  convection = (eqp->flag & CS_EQUATION_CONVECTION) ? true : false;
  bool  diffusion = (eqp->flag & CS_EQUATION_DIFFUSION) ? true : false;
  bool  source_term = (eqp->n_source_terms > 0) ? true : false;

  bft_printf("  <%s/Terms>  unsteady [%s], convection [%s], diffusion [%s],"
             " source term [%s]\n",
             eq->name, cs_base_strtf(unsteady), cs_base_strtf(convection),
             cs_base_strtf(diffusion), cs_base_strtf(source_term));

  /* Boundary conditions */
  if (eqp->verbosity > 0) {
    cs_param_bc_t  *bcp = eqp->bc;

    bft_printf("  <%s/Boundary Conditions>\n", eq->name);
    bft_printf("  <bc> Default BC: %s\n",
               cs_param_get_bc_name(bcp->default_bc));
    if (eqp->verbosity > 1)
      bft_printf("  <bc> Enforcement: %s\n",
                 cs_param_get_bc_enforcement_name(bcp->enforcement));
    bft_printf("  <bc> Number of BCs defined: %d\n", bcp->n_defs);
    if (eqp->verbosity > 1) {
      for (int id = 0; id < bcp->n_defs; id++)
        bft_printf("  <bc> Location: %s; Type: %s; Definition type: %s\n",
                   cs_mesh_location_get_name(bcp->defs[id].loc_id),
                   cs_param_get_bc_name(bcp->defs[id].bc_type),
                   cs_param_get_def_type_name(bcp->defs[id].def_type));
    }
  }

  if (eqp->flag & CS_EQUATION_UNSTEADY) {

    const cs_param_time_t  t_info = eqp->time_info;
    const cs_param_hodge_t  h_info = eqp->time_hodge;

    bft_printf("\n  <%s/Unsteady term>\n", eq->name);
    bft_printf("    <Time> scheme: ");
    switch (t_info.scheme) {
    case CS_TIME_SCHEME_IMPLICIT:
      bft_printf("implicit\n");
      break;
    case CS_TIME_SCHEME_EXPLICIT:
      bft_printf("explicit\n");
      break;
    case CS_TIME_SCHEME_CRANKNICO:
      bft_printf("Crank-Nicolson\n");
      break;
    case CS_TIME_SCHEME_THETA:
      bft_printf("theta scheme with value %f\n", t_info.theta);
      break;
    default:
      bft_error(__FILE__, __LINE__, 0, " Invalid time scheme.");
      break;
    }

    bft_printf("    <Time> Initial condition definition type: %s\n",
               cs_param_get_def_type_name(t_info.ic_def_type));
    bft_printf("    <Unsteady> Mass lumping: %s\n",
               cs_base_strtf(t_info.do_lumping));
    bft_printf("    <Unsteady> Property: %s\n",
               cs_param_pty_get_name(h_info.pty_id));

    if (eqp->verbosity > 0) {
      bft_printf("    <Unsteady> Hodge operator: %s / %s\n",
                 cs_param_hodge_get_type_name(h_info),
                 cs_param_hodge_get_algo_name(h_info));
      bft_printf("    <Unsteady> Inversion of property: %s\n",
                 cs_base_strtf(h_info.inv_pty));
      if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
        bft_printf("    <Unsteady> Value of the Hodge coefficient: %.3e\n",
                   h_info.coef);
    }

  }

  if (eqp->flag & CS_EQUATION_DIFFUSION) {

    const cs_param_hodge_t  h_info = eqp->diffusion_hodge;

    bft_printf("\n  <%s/Diffusion term>\n", eq->name);
    bft_printf("    <Diffusion> Property: %s\n",
               cs_param_pty_get_name(h_info.pty_id));

    if (eqp->verbosity > 0) {
      bft_printf("    <Diffusion> Hodge operator: %s / %s\n",
                 cs_param_hodge_get_type_name(h_info),
                 cs_param_hodge_get_algo_name(h_info));
      bft_printf("    <Diffusion> Inversion of property: %s\n",
                 cs_base_strtf(h_info.inv_pty));
      if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
        bft_printf("    <Diffusion> Value of the Hodge coefficient: %.3e\n",
                   h_info.coef);
    }

  } // Diffusion term

  if (eqp->flag & CS_EQUATION_CONVECTION) {

    const cs_param_advection_t  a_info = eqp->advection;

    bft_printf("\n  <%s/Advection term>\n", eq->name);
    bft_printf("    <Advection field>  %s\n",
               cs_param_adv_field_get_name(a_info.adv_id));

    if (eqp->verbosity > 0) {
      bft_printf("    <Advection operator>");
      switch(a_info.form) {
      case CS_PARAM_ADVECTION_FORM_CONSERV:
        bft_printf(" Conservative formulation");
        break;
      case CS_PARAM_ADVECTION_FORM_NONCONS:
        bft_printf(" Non-conservative formulation");
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid operator type for advection.");
      }

      bft_printf(" & Weight scheme:");
      switch(a_info.weight_algo) {
      case CS_PARAM_ADVECTION_WEIGHT_ALGO_CENTERED:
        bft_printf(" centered\n");
        break;
      case CS_PARAM_ADVECTION_WEIGHT_ALGO_UPWIND:
        bft_printf(" upwind\n");
        break;
      case CS_PARAM_ADVECTION_WEIGHT_ALGO_SAMARSKII:
        bft_printf(" Samarskii\n");
        break;
      case CS_PARAM_ADVECTION_WEIGHT_ALGO_SG:
        bft_printf(" Scharfetter-Gummel\n");
        break;
      case CS_PARAM_ADVECTION_WEIGHT_ALGO_D10G5:
        bft_printf(" Specific with delta=10 and gamma=5\n");
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid weight algorithm for advection.");
      }

      if (eqp->verbosity > 1) {
        bft_printf("    <Evaluation of lambda function>");
        switch(a_info.weight_criterion) {
        case CS_PARAM_ADVECTION_WEIGHT_FLUX:
          bft_printf(" Quadrature on the flux");
          break;
        case CS_PARAM_ADVECTION_WEIGHT_XEXC:
          bft_printf(" Midpoint of [xe, xc]");
          break;
        default:
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid weight criterion for advection.");
        }

      } // verbosity > 1

      bft_printf("\n");

    } // verbosity > 0

  } // Advection term

  if (eqp->n_source_terms > 0) {

    bft_printf("\n  <%s/Source terms>\n", eq->name);
    for (int s_id = 0; s_id < eqp->n_source_terms; s_id++) {

      cs_param_source_term_t  st_info = eqp->source_terms[s_id];

      bft_printf("    <st> name: %s\n", cs_param_source_term_get_name(st_info));
      bft_printf("    <st> mesh location: %s\n",
                 cs_mesh_location_get_name(st_info.ml_id));
      bft_printf("    <st> Type: %s; Variable type: %s; Definition type: %s\n",
                 cs_param_source_term_get_type_name(st_info),
                 cs_param_get_var_type_name(st_info.var_type),
                 cs_param_get_def_type_name(st_info.def_type));
      if (eqp->verbosity > 0)
        bft_printf("    <st> Quadrature type: %s; Use subdivision: %s\n",
                   cs_quadrature_get_type_name(st_info.quad_type),
                   cs_base_strtf(st_info.use_subdiv));

    } // Loop on source terms

  } // Source terms

  /* Iterative solver information */
  const cs_param_itsol_t   itsol = eqp->itsol_info;

  bft_printf("\n  <%s/Sparse Linear Algebra>", eq->name);
  if (eqp->algo_info.type == CS_EQUATION_ALGO_CS_ITSOL)
    bft_printf(" Code_Saturne iterative solvers\n");
  else if (eqp->algo_info.type == CS_EQUATION_ALGO_PETSC_ITSOL)
    bft_printf(" PETSc iterative solvers\n");
  bft_printf("    <sla> Solver.MaxIter     %d\n", itsol.n_max_iter);
  bft_printf("    <sla> Solver.Name        %s\n",
             cs_param_get_solver_name(itsol.solver));
  bft_printf("    <sla> Solver.Precond     %s\n",
             cs_param_get_precond_name(itsol.precond));
  bft_printf("    <sla> Solver.Eps        % -10.6e\n", itsol.eps);
  bft_printf("    <sla> Solver.Normalized  %s\n",
             cs_base_strtf(itsol.resid_normalized));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Assign a set of pointer functions for managing the cs_equation_t
 *         structure during the computation
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_last_setup(cs_equation_t  *eq)
{
  if (eq == NULL)
    return;

  const cs_equation_param_t  *eqp = eq->param;

  /* Set timer statistics */
  if (eqp->verbosity > 0) {

    eq->main_ts_id = cs_timer_stats_create("stages", // parent name
                                           eq->name,
                                           eq->name);

    cs_timer_stats_start(eq->main_ts_id);
    cs_timer_stats_set_plot(eq->main_ts_id, 1);

    if (eqp->verbosity > 1) {

      char *label = NULL;

      int len = strlen("_solve") + strlen(eq->name) + 1;
      BFT_MALLOC(label, len, char);
      sprintf(label, "%s_pre", eq->name);
      eq->pre_ts_id = cs_timer_stats_create(eq->name, label, label);
      cs_timer_stats_set_plot(eq->pre_ts_id, 1);

      sprintf(label, "%s_solve", eq->name);
      eq->solve_ts_id = cs_timer_stats_create(eq->name, label, label);
      cs_timer_stats_set_plot(eq->solve_ts_id, 1);

      sprintf(label, "%s_post", eq->name);
      eq->post_ts_id = cs_timer_stats_create(eq->name, label, label);
      cs_timer_stats_set_plot(eq->post_ts_id, 1);

      BFT_FREE(label);

    } // verbosity > 1

  } // verbosity > 0

  /* Set function pointers */
  switch(eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    eq->init_builder = cs_cdovb_codits_init;
    eq->free_builder = cs_cdovb_codits_free;
    eq->compute_source = cs_cdovb_codits_compute_source;
    eq->build_system = cs_cdovb_codits_build_system;
    eq->update_field = cs_cdovb_codits_update_field;
    eq->get_f_values = NULL;
    break;

  case CS_SPACE_SCHEME_CDOFB:
    eq->init_builder = cs_cdofb_codits_init;
    eq->free_builder = cs_cdofb_codits_free;
    eq->compute_source = cs_cdofb_codits_compute_source;
    eq->build_system = cs_cdofb_codits_build_system;
    eq->update_field = cs_cdofb_codits_update_field;
    eq->get_f_values = cs_cdofb_codits_get_face_values;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid scheme for the space discretization.\n"
                " Please check your settings."));
    break;
  }


  /* Initialize cs_sles_t structure */
  const cs_equation_algo_t  algo = eqp->algo_info;
  const cs_param_itsol_t  itsol = eqp->itsol_info;

  switch (algo.type) {
  case CS_EQUATION_ALGO_CS_ITSOL:
    {
      int  poly_degree = 0; // by default: Jacobi preconditioner

      if (itsol.precond == CS_PARAM_PRECOND_POLY1)
        poly_degree = 1;

      if (itsol.precond != CS_PARAM_PRECOND_POLY1 &&
          itsol.precond != CS_PARAM_PRECOND_DIAG)
        bft_error(__FILE__, __LINE__, 0,
                  " Incompatible preconditioner with Code_Saturne solvers.\n"
                  " Please change your settings (try PETSc ?)");

      switch (itsol.solver) { // Type of iterative solver
      case CS_PARAM_ITSOL_CG:
        cs_sles_it_define(eq->field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_PCG,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_BICG:
        cs_sles_it_define(eq->field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_BICGSTAB2,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_GMRES:
        cs_sles_it_define(eq->field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_GMRES,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_AMG:
        {
          cs_multigrid_t  *mg = cs_multigrid_define(eq->field_id,
                                                    NULL);

          /* Advanced setup (default is specified inside the brackets) */
          cs_multigrid_set_solver_options
            (mg,
             CS_SLES_JACOBI,   // descent smoother type (CS_SLES_PCG)
             CS_SLES_JACOBI,   // ascent smoother type (CS_SLES_PCG)
             CS_SLES_PCG,      // coarse solver type (CS_SLES_PCG)
             itsol.n_max_iter, // n max cycles (100)
             5,                // n max iter for descent (10)
             5,                // n max iter for asscent (10)
             1000,             // n max iter coarse solver (10000)
             0,                // polynomial precond. degree descent (0)
             0,                // polynomial precond. degree ascent (0)
             0,                // polynomial precond. degree coarse (0)
             1.0,    // precision multiplier descent (< 0 forces max iters)
             1.0,    // precision multiplier ascent (< 0 forces max iters)
             1);     // requested precision multiplier coarse (default 1)

        }
      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Undefined iterative solver for solving %s equation.\n"
                    " Please modify your settings."), eq->name);
        break;
      } // end of switch

    } // Solver provided by Code_Saturne
    break;

  case CS_EQUATION_ALGO_PETSC_ITSOL:
    {
#if defined(HAVE_PETSC)

      /* Initialization must be called before setting options;
         it does not need to be called before calling
         cs_sles_petsc_define(), as this is handled automatically. */

      PetscBool is_initialized;
      PetscInitialized(&is_initialized);
      if (is_initialized == PETSC_FALSE) {
        PETSC_COMM_WORLD = cs_glob_mpi_comm;
        PetscInitializeNoArguments();
      }

      switch (eqp->itsol_info.solver) {

      case CS_PARAM_ITSOL_CG:
        switch (eqp->itsol_info.precond) {

        case CS_PARAM_PRECOND_DIAG:
          cs_sles_petsc_define(eq->field_id,
                               NULL,
                               MATMPIAIJ,
                               _cg_diag_setup_hook,
                               NULL);
          break;
        case CS_PARAM_PRECOND_SSOR:
          cs_sles_petsc_define(eq->field_id,
                               NULL,
                               MATSEQAIJ,
                               _cg_ssor_setup_hook,
                               NULL);
          break;
        case CS_PARAM_PRECOND_ICC0:
          cs_sles_petsc_define(eq->field_id,
                               NULL,
                               MATSEQAIJ, // Warning SEQ not MPI
                               _cg_icc_setup_hook,
                               NULL);
          break;

        case CS_PARAM_PRECOND_AMG:
          {
            int  amg_type = 1;

            if (amg_type == 0) { // GAMG

              PetscOptionsSetValue("-pc_gamg_agg_nsmooths", "1");
              PetscOptionsSetValue("-mg_levels_ksp_type", "richardson");
              PetscOptionsSetValue("-mg_levels_pc_type", "sor");
              PetscOptionsSetValue("-mg_levels_ksp_max_it", "1");
              PetscOptionsSetValue("-pc_gamg_threshold", "0.02");
              PetscOptionsSetValue("-pc_gamg_reuse_interpolation", "TRUE");
              PetscOptionsSetValue("-pc_gamg_square_graph", "4");

              cs_sles_petsc_define(eq->field_id,
                                   NULL,
                                   MATMPIAIJ,
                                   _cg_gamg_setup_hook,
                                   NULL);

            }
            else if (amg_type == 1) { // Boomer AMG (hypre)

              PetscOptionsSetValue("-pc_type", "hypre");
              PetscOptionsSetValue("-pc_hypre_type","boomeramg");
              PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type","HMIS");
              PetscOptionsSetValue("-pc_hypre_boomeramg_interp_type","ext+i-cc");
              PetscOptionsSetValue("-pc_hypre_boomeramg_agg_nl","2");
              PetscOptionsSetValue("-pc_hypre_boomeramg_P_max","4");
              PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold","0.5");
              PetscOptionsSetValue("-pc_hypre_boomeramg_no_CF","");

              cs_sles_petsc_define(eq->field_id,
                                   NULL,
                                   MATMPIAIJ,
                                   _cg_bamg_setup_hook,
                                   NULL);

            }

          }
          break;
        case CS_PARAM_PRECOND_AS:
          cs_sles_petsc_define(eq->field_id,
                               NULL,
                               MATMPIAIJ,
                               _cg_as_setup_hook,
                               NULL);
          break;
        default:
          bft_error(__FILE__, __LINE__, 0,
                    " Couple (solver, preconditioner) not handled with PETSc.");
          break;

        } // switch on PETSc preconditionner
        break;

      case CS_PARAM_ITSOL_GMRES:

        switch (eqp->itsol_info.precond) {
        case CS_PARAM_PRECOND_ILU0:
          cs_sles_petsc_define(eq->field_id,
                               NULL,
                               MATSEQAIJ, // Warning SEQ not MPI
                               _gmres_ilu_setup_hook,
                               NULL);
          break;
        case CS_PARAM_PRECOND_DIAG:
          cs_sles_petsc_define(eq->field_id,
                               NULL,
                               MATMPIAIJ,
                               _gmres_bjacobi_setup_hook,
                               NULL);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    " Couple (solver, preconditioner) not handled with PETSc.");
          break;

        } // switch on PETSc preconditionner
        break;

      case CS_PARAM_ITSOL_BICG:

        switch (eqp->itsol_info.precond) {
        case CS_PARAM_PRECOND_ILU0:
          cs_sles_petsc_define(eq->field_id,
                               NULL,
                               MATSEQAIJ, // Warning SEQ not MPI
                               _bicg_ilu_setup_hook,
                               NULL);
          break;
        case CS_PARAM_PRECOND_DIAG:
          cs_sles_petsc_define(eq->field_id,
                               NULL,
                               MATMPIAIJ,
                               _bicg_bjacobi_setup_hook,
                               NULL);
          break;

        default:
          bft_error(__FILE__, __LINE__, 0,
                    " Couple (solver, preconditioner) not handled with PETSc.");
          break;

        } // switch on PETSc preconditionner
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, " Solver not handled.");
        break;

      } // switch on PETSc solver
#else
      bft_error(__FILE__, __LINE__, 0,
                _(" PETSC algorithms used to solve %s are not linked.\n"
                  " Please install Code_Saturne with PETSc."), eq->name);

#endif // HAVE_PETSC
    } // Solver provided by PETSc
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Algorithm requested to solve %s is not implemented yet.\n"
                " Please modify your settings."), eq->name);
    break;

  } // end switch on algorithms

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an equation
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       keyname   name of key related to the member of eq to set
 * \param[in]       val       accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set(cs_equation_t       *eq,
                const char          *keyname,
                const void          *val)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop setting an empty cs_equation_t structure.\n"
                " Please check your settings.\n"));

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_equation_param_t  *eqp = eq->param;
  eqkey_t  key = _get_eqkey(keyname);

  if (key == EQKEY_ERROR) {
    bft_printf("\n\n Current key: %s\n", keyname);
    bft_printf(" Possible keys: ");
    for (int i = 0; i < EQKEY_ERROR; i++) {
      bft_printf("%s ", _print_eqkey(i));
      if (i > 0 && i % 3 == 0)
        bft_printf("\n\t");
    }
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting equation %s.\n"
                " Please read listing for more details and"
                " modify your settings."), eq->name);

  } /* Error message */

  switch(key) {

  case EQKEY_SPACE_SCHEME:
    if (strcmp(val, "cdo_vb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
    }
    else if (strcmp(val, "cdo_fb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOFB;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_CPVD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EDFP;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key %s\n"
                  " Choice between cdo_vb or cdo_fb"), _val, keyname);
    }
    break;

  case EQKEY_HODGE_DIFF_ALGO:
    if (strcmp(val,"cost") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
    else if (strcmp(val, "voronoi") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key %s\n"
                  " Choice between cost or voronoi"), _val, keyname);
    }
    break;

  case EQKEY_HODGE_TIME_ALGO:
    if (strcmp(val,"cost") == 0)
      eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
    else if (strcmp(val, "voronoi") == 0)
      eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key %s\n"
                  " Choice between cost or voronoi"), _val, keyname);
    }
    break;

  case EQKEY_HODGE_DIFF_COEF:
    if (strcmp(val, "dga") == 0)
      eqp->diffusion_hodge.coef = 1./3.;
    else if (strcmp(val, "sushi") == 0)
      eqp->diffusion_hodge.coef = 1./sqrt(3.);
    else if (strcmp(val, "gcr") == 0)
      eqp->diffusion_hodge.coef = 1.0;
    else
      eqp->diffusion_hodge.coef = atof(val);
    break;

  case EQKEY_HODGE_TIME_COEF:
    if (strcmp(val, "") == 0)
      eqp->time_hodge.coef = 1./3.;
    else if (strcmp(val, "sushi") == 0)
      eqp->time_hodge.coef = 1./sqrt(3.);
    else if (strcmp(val, "gcr") == 0)
      eqp->time_hodge.coef = 1.0;
    else
      eqp->time_hodge.coef = atof(val);
    break;

  case EQKEY_SOLVER_FAMILY:
    if (strcmp(val, "cs") == 0)
      eqp->algo_info.type = CS_EQUATION_ALGO_CS_ITSOL;
    else if (strcmp(val, "petsc") == 0)
      eqp->algo_info.type = CS_EQUATION_ALGO_PETSC_ITSOL;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key %s\n"
                  " Choice between cs or petsc"), _val, keyname);
    }
    break;

  case EQKEY_ITSOL:
    if (strcmp(val, "cg") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_CG;
    else if (strcmp(val, "bicg") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_BICG;
    else if (strcmp(val, "gmres") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_GMRES;
    else if (strcmp(val, "amg") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_AMG;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key %s\n"
                  " Choice between cg, bicg, gmres or amg"), _val, keyname);
    }
    break;

  case EQKEY_PRECOND:
    if (strcmp(val, "jacobi") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_DIAG;
    else if (strcmp(val, "poly1") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_POLY1;
    else if (strcmp(val, "ssor") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_SSOR;
    else if (strcmp(val, "ilu0") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_ILU0;
    else if (strcmp(val, "icc0") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_ICC0;
    else if (strcmp(val, "amg") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_AMG;
    else if (strcmp(val, "as") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_AS;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key %s\n"
                  " Choice between jacobi, poly1, ssor, ilu0,\n"
                  " icc0, amg or as"), _val, keyname);
    }
    break;

  case EQKEY_ITSOL_MAX_ITER:
    eqp->itsol_info.n_max_iter = atoi(val);
    break;

  case EQKEY_ITSOL_EPS:
    eqp->itsol_info.eps = atof(val);
    break;

  case EQKEY_ITSOL_RESNORM:
    if (strcmp(val, "true") == 0)
      eqp->itsol_info.resid_normalized = true;
    else if (strcmp(val, "false") == 0)
      eqp->itsol_info.resid_normalized = false;
    break;

  case EQKEY_VERBOSITY: // "verbosity"
    eqp->verbosity = atoi(val);
    break;

  case EQKEY_BC_ENFORCEMENT:
    if (strcmp(val, "strong") == 0)
      eqp->bc->enforcement = CS_PARAM_BC_ENFORCE_STRONG;
    else if (strcmp(val, "penalization") == 0)
      eqp->bc->enforcement = CS_PARAM_BC_ENFORCE_WEAK_PENA;
    else if (strcmp(val, "weak_sym") == 0)
      eqp->bc->enforcement = CS_PARAM_BC_ENFORCE_WEAK_SYM;
    else if (strcmp(val, "weak") == 0)
      eqp->bc->enforcement = CS_PARAM_BC_ENFORCE_WEAK_NITSCHE;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value %s related to key %s\n"
                  " Choice between strong, penalization, weak or\n"
                  " weak_sym."), _val, keyname);
    }
    break;

  case EQKEY_BC_QUADRATURE:
    if (strcmp(val, "subdiv") == 0)
      eqp->bc->use_subdiv = true;
    else if (strcmp(val, "bary") == 0)
      eqp->bc->quad_type = CS_QUADRATURE_BARY;
    else if (strcmp(val, "higher") == 0)
      eqp->bc->quad_type = CS_QUADRATURE_HIGHER;
    else if (strcmp(val, "highest") == 0)
      eqp->bc->quad_type = CS_QUADRATURE_HIGHEST;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid key value %s for setting the quadrature behaviour"
                  " of boundary conditions.\n"
                  " Choices are among subdiv, bary, higher and highest."), _val);
    }
    break;

  case EQKEY_OUTPUT_FREQ:
    eqp->output_freq = atoi(val);
    break;

  case EQKEY_POST_FREQ:
    eqp->post_freq = atoi(val);
    break;

  case EQKEY_POST:
    if (strcmp(val, "peclet") == 0)
      eqp->post_flag |= CS_EQUATION_POST_PECLET;
    else if (strcmp(val, "upwind_coef") == 0)
      eqp->post_flag |= CS_EQUATION_POST_UPWIND_COEF;
    break;

  case EQKEY_ADV_OP_TYPE:
    if (strcmp(val, "conservative") == 0)
      eqp->advection.form = CS_PARAM_ADVECTION_FORM_CONSERV;
    else if (strcmp(val, "non_conservative") == 0)
      eqp->advection.form = CS_PARAM_ADVECTION_FORM_NONCONS;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid key value %s for setting the form of the convection"
                  " term.\n"
                  " Choices are among conservative and non_conservative."),
                _val);
    }
    break;

  case EQKEY_ADV_WEIGHT_ALGO:
    if (strcmp(val, "upwind") == 0)
      eqp->advection.weight_algo = CS_PARAM_ADVECTION_WEIGHT_ALGO_UPWIND;
    else if (strcmp(val, "samarskii") == 0)
      eqp->advection.weight_algo = CS_PARAM_ADVECTION_WEIGHT_ALGO_SAMARSKII;
    else if (strcmp(val, "sg") == 0)
      eqp->advection.weight_algo = CS_PARAM_ADVECTION_WEIGHT_ALGO_SG;
    else if (strcmp(val, "d10g5") == 0)
      eqp->advection.weight_algo = CS_PARAM_ADVECTION_WEIGHT_ALGO_D10G5;
    else if (strcmp(val, "centered") == 0)
      eqp->advection.weight_algo = CS_PARAM_ADVECTION_WEIGHT_ALGO_CENTERED;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid key value %s for setting the algorithm for defining"
                  " the proportion of upwinding.\n"
                  " Choices are among upwind, samarskii, sg and centered."),
                _val);
    }
    break;

  case EQKEY_ADV_WEIGHT_CRIT:
    if (strcmp(val, "xexc") == 0)
      eqp->advection.weight_criterion = CS_PARAM_ADVECTION_WEIGHT_XEXC;
    else if (strcmp(val, "flux") == 0)
      eqp->advection.weight_criterion = CS_PARAM_ADVECTION_WEIGHT_FLUX;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid key value %s for setting the algorithm for"
                  " computing the upwinding weight.\n"
                  " Choices are among flux and xexc."),
                _val);
    }
    break;

  case EQKEY_ADV_FLUX_QUADRA:
    if (strcmp(val, "bary") == 0)
      eqp->advection.quad_type = CS_QUADRATURE_BARY;
    else if (strcmp(val, "higher") == 0)
      eqp->advection.quad_type = CS_QUADRATURE_HIGHER;
    else if (strcmp(val, "highest") == 0)
      eqp->advection.quad_type = CS_QUADRATURE_HIGHEST;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid key value %s for setting the quadrature behaviour"
                  " used for computing the advection flux.\n"
                  " Choices are among bary, higher and highest."), _val);
    }
    break;

  case EQKEY_TIME_SCHEME:
    if (strcmp(val, "implicit") == 0) {
      eqp->time_info.scheme = CS_TIME_SCHEME_IMPLICIT;
      eqp->time_info.theta = 1.;
    }
    else if (strcmp(val, "explicit") == 0) {
      eqp->time_info.scheme = CS_TIME_SCHEME_EXPLICIT;
      eqp->time_info.theta = 0.;
    }
    else if (strcmp(val, "crank_nicolson") == 0) {
      eqp->time_info.scheme = CS_TIME_SCHEME_CRANKNICO;
      eqp->time_info.theta = 0.5;
    }
    else if (strcmp(val, "theta_scheme") == 0)
      eqp->time_info.scheme = CS_TIME_SCHEME_THETA;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid key value %s for setting the time scheme.\n"
                  " Choices are among implicit, explicit, crank_nicolson"
                  " and theta_scheme"), _val);
    }
    break;

  case EQKEY_TIME_THETA:
    eqp->time_info.theta = atof(val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key %s is not implemented yet."), keyname);

  } /* Switch on keys */

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Associate a material property or an advection field with an equation
 *         for a given term (diffusion, time, convection)
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       keyword   "time", "diffusion", "advection"...
 * \param[in]       name      name of the property to associate
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_link(cs_equation_t       *eq,
                 const char          *keyword,
                 const char          *name)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_equation_t structure to set with property %s is NULL\n"),
              name);

  cs_equation_param_t  *eqp = eq->param;

  if (strcmp("diffusion", keyword) == 0)
    eqp->diffusion_hodge.pty_id = cs_param_pty_get_id_by_name(name);
  else if (strcmp("time", keyword) == 0)
    eqp->time_hodge.pty_id = cs_param_pty_get_id_by_name(name);
  else if (strcmp("advection", keyword) == 0)
    eqp->advection.adv_id = cs_param_adv_get_id_by_name(name);
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting a property.\n"
                " Current value: %s\n"
                " Possible choice: diffusion or time\n"), keyword);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition of the unknown related to this equation
 *         def_key is among "value", "analytic", "user"
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       def_key   way of defining the value of the bc
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_ic(cs_equation_t    *eq,
                   const char       *def_key,
                   void             *val)
{
  cs_param_var_type_t  var_type = CS_PARAM_N_VAR_TYPES;

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_equation_t structure is NULL\n"
                " Cannot set the initial condition"));

  cs_equation_param_t  *eqp = eq->param;

  /* Get the type of variable */
  switch (eqp->type) {
  case CS_EQUATION_TYPE_SCAL:
    var_type = CS_PARAM_VAR_SCAL;
    break;
  case CS_EQUATION_TYPE_VECT:
    var_type = CS_PARAM_VAR_VECT;
    break;
  case CS_EQUATION_TYPE_TENS:
    var_type = CS_PARAM_VAR_TENS;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of equation for equation %s.\n"), eq->name);
    break;
  }

  /* Get the type of definition */
  if (strcmp(def_key, "value") == 0)
    eqp->time_info.ic_def_type = CS_PARAM_DEF_BY_VALUE;
  else if (strcmp(def_key, "analytic") == 0)
    eqp->time_info.ic_def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the initial condition.\n"
                " Given key: %s\n"
                " Choice among value and analytic.\n"
                " Please modify your settings."), def_key);

  /* Get the definition */
  cs_param_set_def(eqp->time_info.ic_def_type, var_type, val,
                   &(eqp->time_info.ic_def));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an equation
 *         bc_key among "dirichlet", "neumann" or "robin"
 *         def_key among "value", "analytic", "user"
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       bc_key    type of boundary condition to add
 * \param[in]       def_key   way of defining the value of the bc
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc(cs_equation_t    *eq,
                   const char       *ml_name,
                   const char       *bc_key,
                   const char       *def_key,
                   const void       *val)
{
  int  ml_id;

  cs_param_bc_type_t  bc_type = CS_PARAM_N_BC_TYPES;
  cs_param_def_type_t  def_type = CS_PARAM_N_DEF_TYPES;
  cs_param_var_type_t  var_type = CS_PARAM_N_VAR_TYPES;

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_equation_t structure is NULL\n"
                " Cannot add a boundary condition related to mesh location %s"),
              ml_name);

  cs_equation_param_t  *eqp = eq->param;
  cs_param_bc_t  *bc = eqp->bc;

  /* Sanity checks */
  assert(bc != NULL);

  /* Add a new definition */
  int def_id = bc->n_defs;
  bc->n_defs += 1;
  BFT_REALLOC(bc->defs, bc->n_defs, cs_param_bc_def_t);

  /* Get the mesh location id from its name */
  _check_ml_name(ml_name, &ml_id);

  /* Get the type of variable */
  switch (eqp->type) {
  case CS_EQUATION_TYPE_SCAL:
    var_type = CS_PARAM_VAR_SCAL;
    break;
  case CS_EQUATION_TYPE_VECT:
    var_type = CS_PARAM_VAR_VECT;
    break;
  case CS_EQUATION_TYPE_TENS:
    var_type = CS_PARAM_VAR_TENS;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of equation for equation %s.\n"), eq->name);
    break;
  }

  /* Get the type of definition */
  if (strcmp(def_key, "value") == 0)
    def_type = CS_PARAM_DEF_BY_VALUE;
  else if (strcmp(def_key, "field") == 0)
    def_type = CS_PARAM_DEF_BY_FIELD;
  else if (strcmp(def_key, "evaluator") == 0)
    def_type = CS_PARAM_DEF_BY_EVALUATOR;
  else if (strcmp(def_key, "analytic") == 0)
    def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  else if (strcmp(def_key, "user") == 0)
    def_type = CS_PARAM_DEF_BY_USER_FUNCTION;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the type of definition.\n"
                " Given key: %s\n"
                " Choice among value, field, evaluator, analytic, user, law"
                " or file\n"
                " Please modify your settings."), def_key);

  /* Get the type of boundary condition */
  if (strcmp(bc_key, "dirichlet") == 0)
    bc_type = CS_PARAM_BC_DIRICHLET;
  else if (strcmp(bc_key, "neumann") == 0)
    bc_type = CS_PARAM_BC_NEUMANN;
  else if (strcmp(bc_key, "robin") == 0)
    bc_type = CS_PARAM_BC_ROBIN;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the type of boundary condition.\n"
                " Given key: %s\n"
                " Choice among dirichlet, neumann or robin.\n"
                " Please modify your settings."), bc_key);

  /* Check if this is a homogeneous boundary condition */
  if (def_type == CS_PARAM_DEF_BY_VALUE && var_type == CS_PARAM_VAR_SCAL) {
    cs_real_t  value = atof(val);
    if (fabs(value) < DBL_MIN) {
      if (bc_type == CS_PARAM_BC_DIRICHLET)
        bc_type = CS_PARAM_BC_HMG_DIRICHLET;
      if (bc_type == CS_PARAM_BC_NEUMANN)
        bc_type = CS_PARAM_BC_HMG_NEUMANN;
    }
  }

  cs_param_bc_def_set(bc->defs + def_id,
                      ml_id,
                      bc_type,
                      var_type,
                      def_type,
                      val, NULL); // coef2 is not used up to now
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to a source term
 *         st_key among "implicit", "explicit", "imex"...
 *         def_key among "value", "analytic", "user"...
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       st_name   name of the source term or NULL
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       st_key    type of boundary condition to add
 * \param[in]       def_key   way of defining the value of the bc
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_source_term(cs_equation_t   *eq,
                            const char      *st_name,
                            const char      *ml_name,
                            const char      *st_key,
                            const char      *def_key,
                            const void      *val)
{
  int  ml_id;
  char *_st_name = NULL;

  const char  *name;

  cs_param_source_term_type_t  st_type = CS_PARAM_N_SOURCE_TERM_TYPES;
  cs_param_def_type_t  def_type = CS_PARAM_N_DEF_TYPES;
  cs_param_var_type_t  var_type = CS_PARAM_N_VAR_TYPES;

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_equation_t structure is NULL\n"
                " Can not add a source term related to mesh location %s"),
              ml_name);

  cs_equation_param_t  *eqp = eq->param;

  /* Add a new source term */
  int  st_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_param_source_term_t);

  if (st_name == NULL) { /* Define a name by default */
    assert(st_id < 100);
    int len = strlen("sourceterm_00") + 1;
    BFT_MALLOC(_st_name, len, char);
    sprintf(_st_name, "sourceterm_%2d", st_id);
    name = _st_name;
  }
  else
    name = st_name;

  /* Get the mesh location id from its name */
  _check_ml_name(ml_name, &ml_id);

  /* Get the type of variable */
  switch (eqp->type) {
  case CS_EQUATION_TYPE_SCAL:
    var_type = CS_PARAM_VAR_SCAL;
    break;
  case CS_EQUATION_TYPE_VECT:
    var_type = CS_PARAM_VAR_VECT;
    break;
  case CS_EQUATION_TYPE_TENS:
    var_type = CS_PARAM_VAR_TENS;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of equation for equation %s.\n"), eq->name);
    break;
  }

  /* Get the type of definition */
  if (strcmp(def_key, "value") == 0)
    def_type = CS_PARAM_DEF_BY_VALUE;
  else if (strcmp(def_key, "field") == 0)
    def_type = CS_PARAM_DEF_BY_FIELD;
  else if (strcmp(def_key, "evaluator") == 0)
    def_type = CS_PARAM_DEF_BY_EVALUATOR;
  else if (strcmp(def_key, "analytic") == 0)
    def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  else if (strcmp(def_key, "user") == 0)
    def_type = CS_PARAM_DEF_BY_USER_FUNCTION;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the type of definition.\n"
                " Given key: %s\n"
                " Choice among value, field, evaluator, analytic, user, law"
                " or file\n"
                " Please modify your settings."), def_key);

  /* Get the type of source term */
  if (strcmp(st_key, "implicit") == 0)
    st_type = CS_PARAM_SOURCE_TERM_IMPLICIT;
  else if (strcmp(st_key, "explicit") == 0)
    st_type = CS_PARAM_SOURCE_TERM_EXPLICIT;
  else if (strcmp(st_key, "imex") == 0)
    st_type = CS_PARAM_SOURCE_TERM_IMEX;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the type of boundary condition.\n"
                " Given key: %s\n"
                " Choice among dirichlet, neumann or robin.\n"
                " Please modify your settings."), st_key);

  cs_param_source_term_add(eqp->source_terms + st_id,
                           name,
                           ml_id,
                           st_type,
                           var_type,
                           CS_QUADRATURE_BARY,    // default value
                           def_type,
                           val);

  if (st_name == NULL)
    BFT_FREE(_st_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set advanced parameters which are members defined by default in a
 *         source term structure.
 *         keyname among "quadrature", "post"...

 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       st_name   name of the source term
 * \param[in]       keyname   name of the key
 * \param[in]       keyval    pointer to the value to set to the key
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_source_term_set(cs_equation_t    *eq,
                            const char       *st_name,
                            const char       *keyname,
                            const char       *keyval)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Stop setting an empty cs_equation_t structure.\n"
                " Please check your settings.\n"));

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_equation_param_t  *eqp = eq->param;

  /* Look for the requested source term */
  int  st_id = -1;
  for (int id = 0; id < eqp->n_source_terms; id++) {
    if (strcmp(eqp->source_terms[id].name, st_name) == 0) {
      st_id = id;
      break;
    }
  }

  if (st_id == -1) // Error
    bft_error(__FILE__, __LINE__, 0,
              _(" Cannot find source term %s.\n"
                " Please check your settings.\n"), st_name);

  stkey_t  key = _get_stkey(keyname);

  if (key == STKEY_ERROR) {
    bft_printf("\n\n Current key: %s\n", keyname);
    bft_printf(" Possible keys: ");
    for (int i = 0; i < STKEY_ERROR; i++) {
      bft_printf("%s ", _print_stkey(i));
      if (i > 0 && i % 3 == 0)
        bft_printf("\n\t");
    }
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting source term %s.\n"
                " Please read listing for more details and"
                " modify your settings."), st_name);

  } /* Error message */

  switch(key) {

  case STKEY_POST:
    eqp->source_terms[st_id].post = atoi(keyval);
    break;

  case STKEY_QUADRATURE:
    if (strcmp(keyval, "subdiv") == 0)
      eqp->source_terms[st_id].use_subdiv = true;
    else if (strcmp(keyval, "bary") == 0)
      eqp->source_terms[st_id].quad_type = CS_QUADRATURE_BARY;
    else if (strcmp(keyval, "higher") == 0)
      eqp->source_terms[st_id].quad_type = CS_QUADRATURE_HIGHER;
    else if (strcmp(keyval, "highest") == 0)
      eqp->source_terms[st_id].quad_type = CS_QUADRATURE_HIGHEST;
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid key value for setting the quadrature behaviour"
                  " of a source term.\n"
                  " Choices are among subdiv, bary, higher, highest."));
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key %s is not implemented yet."), keyname);

  } /* Switch on keys */

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a field structure related to this cs_equation_t structure
 *         to an equation
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_create_field(cs_equation_t     *eq)
{
  int  dim = 0, location_id = -1; // initialize values to avoid a warning

  int  field_mask = CS_FIELD_INTENSIVE | CS_FIELD_VARIABLE;

  /* Sanity check */
  assert(eq != NULL);

  const cs_equation_param_t  *eqp = eq->param;

  _Bool has_previous = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  /* Define dim */
  switch (eqp->type) {
  case CS_EQUATION_TYPE_SCAL:
    dim = 1;
    break;
  case CS_EQUATION_TYPE_VECT:
    dim = 3;
    break;
  case CS_EQUATION_TYPE_TENS:
    dim = 9;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Type of equation for eq. %s is incompatible with the"
                " creation of a field structure.\n"), eq->name);
    break;
  }

  /* Associate a predefined mesh_location_id to this field */
  switch (eqp->space_scheme) {
  case CS_SPACE_SCHEME_CDOVB:
    location_id = cs_mesh_location_get_id_by_name(N_("vertices"));
    break;
  case CS_SPACE_SCHEME_CDOFB:
    location_id = cs_mesh_location_get_id_by_name(N_("cells"));
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Space scheme for eq. %s is incompatible with a field.\n"
                " Stop adding a cs_field_t structure.\n"), eq->name);
    break;
  }

  if (location_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid mesh location id (= -1) for the current field\n"));

  cs_field_t  *fld = cs_field_create(eq->varname,
                                     field_mask,
                                     location_id,
                                     dim,
                                     true,          // interleave
                                     has_previous);

  /* Store the related field id */
  eq->field_id = cs_field_id_by_name(eq->varname);

  /* Allocate and initialize values */
  cs_field_allocate_values(fld);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the values of a field according to the initial condition
 *         related to its equation
 *
 * \param[in]       mesh       pointer to the mesh structure
 * \param[in]       connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]       cdoq       pointer to a cs_cdo_quantities_t struct.
 * \param[in]       time_step  pointer to a time step structure
 * \param[in, out]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_system(const cs_mesh_t            *mesh,
                        const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *cdoq,
                        const cs_time_step_t       *time_step,
                        cs_equation_t              *eq)
{
  if (eq == NULL)
    return;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  const double t_ini = 0;
  const cs_equation_param_t  *eqp = eq->param;

  /* Allocate and initialize a system builder */
  eq->builder = eq->init_builder(eqp, mesh);

  /* Compute the (initial) source term */
  eq->compute_source(mesh, connect, cdoq, time_step, eq->builder);

  /* Initialize the associated field to the initial conditino */
  if (!(eqp->flag & CS_EQUATION_UNSTEADY)) // Steady equation do not need this
    return;

  int  k, l, stride;
  cs_lnum_t  i, n_elts;
  cs_get_t  get;

  /* Retrieve the associated field */
  cs_field_t  *field = cs_field_by_id(eq->field_id);

  /* Define dim */
  switch (eqp->type) {
  case CS_EQUATION_TYPE_SCAL:
    stride = 1;
    break;
  case CS_EQUATION_TYPE_VECT:
    stride = 3;
    break;
  case CS_EQUATION_TYPE_TENS:
    stride = 9;
    break;
  default:
    stride = 0; // avoid a warning
    bft_error(__FILE__, __LINE__, 0,
              _(" Type of equation for eq. %s is incompatible with the"
                " creation of a field structure.\n"), eq->name);
    break;
  }

  /* Associate a predefined mesh_location_id to this field */
  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB)
    n_elts = cdoq->n_vertices;
  else {
    assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOFB);
    n_elts = cdoq->n_faces;
  }

  /* Get the definition of the initial condition */
  cs_def_t  def = eqp->time_info.ic_def;

  if (eqp->time_info.ic_def_type == CS_PARAM_DEF_BY_VALUE) {

    if (stride == 1)
      for (i = 0; i < n_elts; i++)
        field->val[i] = def.get.val;
    else if (stride == 3) { // Interleave by construction
      for (i = 0; i < n_elts; i++)
        for (k = 0; k < 3; k++)
          field->val[3*i+k] = def.get.vect[k];
    }
    else if (stride == 9) {
      for (i = 0; i < n_elts; i++)
        for (k = 0; k < 3; k++)
          for (l = 0; l < 3; l++)
          field->val[9*i+3*k+l] = def.get.tens[k][l];
    }

  }
  else if (eqp->time_info.ic_def_type == CS_PARAM_DEF_BY_ANALYTIC_FUNCTION) {

    if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB) {

      for (i = 0; i < n_elts; i++)  {
        def.analytic(t_ini, &(mesh->vtx_coord[3*i]), &get);
        if (stride == 1)
          field->val[i] = def.get.val;
        else if (stride == 3) { // Interleave by construction
          for (k = 0; k < 3; k++)
            field->val[3*i+k] = get.vect[k];
        }
        else { // stride = 9
          for (k = 0; k < 3; k++)
            for (l = 0; l < 3; l++)
              field->val[9*i+3*k+l] = get.tens[k][l];
        }
      } // Loop on vertices

    }
    else { // Face-based schemes

      for (i = 0; i < n_elts; i++)  {
        def.analytic(t_ini, cdoq->face[i].center, &get);
        if (stride == 1)
          field->val[i] = get.val;
        else if (stride == 3) { // Interleave by construction
          for (k = 0; k < 3; k++)
            field->val[3*i+k] = get.vect[k];
        }
        else {
          for (k = 0; k < 3; k++)
            for (l = 0; l < 3; l++)
              field->val[9*i+3*k+l] = get.tens[k][l];
        }

      } // Loop on faces

    } // Test on discretization scheme

  } /* Definition using an analytical function */

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if one has to build the linear system
 *
 * \param[in]  eq        pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_needs_build(const cs_equation_t    *eq)
{
  return eq->do_build;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system for this equation
 *
 * \param[in]       m          pointer to a cs_mesh_t structure
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]       time_step  pointer to a time step structure
 * \param[in]       dt_cur     value of the current time step
 * \param[in, out]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_system(const cs_mesh_t            *m,
                         const cs_cdo_connect_t     *connect,
                         const cs_cdo_quantities_t  *cdoq,
                         const cs_time_step_t       *time_step,
                         double                      dt_cur,
                         cs_equation_t              *eq)
{
  cs_sla_matrix_t  *sla_mat = NULL;

  const char *eqn = eq->name;
  const cs_equation_param_t  *eqp = eq->param;
  const cs_field_t  *fld = cs_field_by_id(eq->field_id);

  if (eq->pre_ts_id > -1)
    cs_timer_stats_start(eq->pre_ts_id);

  eq->build_system(m, connect, cdoq, time_step, dt_cur, fld->val,
                   eq->builder,
                   &(eq->rhs),
                   &(sla_mat));

  /* Get information on the matrix related to this linear system */
  cs_sla_matrix_info_t  minfo = cs_sla_matrix_analyse(sla_mat);

  if (eqp->verbosity > 1) {
    if (time_step->nt_cur % eqp->output_freq == 0) {
      bft_printf("\n Sparse Linear Algebra (SLA) sumup:\n");
      bft_printf("  <%s/sla> A.size         %d\n", eqn, sla_mat->n_rows);
      bft_printf("  <%s/sla> A.nnz          %lu\n", eqn, minfo.nnz);
      bft_printf("  <%s/sla> A.FillIn       %5.2e %%\n", eqn, minfo.fillin);
      bft_printf("  <%s/sla> A.StencilMin   %d\n", eqn, minfo.stencil_min);
      bft_printf("  <%s/sla> A.StencilMax   %d\n", eqn, minfo.stencil_max);
      bft_printf("  <%s/sla> A.StencilMean  %5.2e\n", eqn, minfo.stencil_mean);
    }
  }

  /* Map a cs_sla_matrix_t structure into a cs_matrix_t structure */
  assert(sla_mat->type == CS_SLA_MAT_MSR);

  /* First step: create a matrix structure */
  if (eq->ms == NULL)
    eq->ms = cs_matrix_structure_create_msr(CS_MATRIX_MSR,      // type
                                            true,               // transfer
                                            true,               // have_diag
                                            sla_mat->n_rows,    // n_rows
                                            sla_mat->n_cols,    // n_cols_ext
                                            &(sla_mat->idx),    // row_index
                                            &(sla_mat->col_id), // col_id
                                            NULL,               // halo
                                            NULL);              // numbering

  if (eq->matrix == NULL)
    eq->matrix = cs_matrix_create(eq->ms); // ms is also stored inside matrix

  const cs_lnum_t  *row_index, *col_id;
  cs_matrix_get_msr_arrays(eq->matrix, &row_index, &col_id, NULL, NULL);

  /* Second step: associate coefficients to a matrix structure */
  cs_matrix_transfer_coefficients_msr(eq->matrix,
                                      false,             // symmetric values ?
                                      NULL,              // diag. block
                                      NULL,              // extra-diag. block
                                      row_index,         // row_index
                                      col_id,            // col_id
                                      &(sla_mat->diag),  // diag. values
                                      &(sla_mat->val));  // extra-diag. values

  /* Free non-transferred parts of sla_mat */
  sla_mat = cs_sla_matrix_free(sla_mat);

  eq->do_build = false;
  if (eqp->flag & CS_EQUATION_UNSTEADY)
    eq-> do_build = true; /* Improvement: exhibit cases where a new build
                             is not needed */

  if (eq->pre_ts_id > -1)
    cs_timer_stats_stop(eq->pre_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system for this equation
 *
 * \param[in]       connect    pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]       time_step  pointer to a time step structure
 * \param[in, out]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(const cs_cdo_connect_t     *connect,
                  const cs_cdo_quantities_t  *cdoq,
                  const cs_time_step_t       *time_step,
                  cs_equation_t              *eq)
{
  double  r_norm;
  cs_sles_convergence_state_t  cvg;
  cs_sla_sumup_t  ret;

  if (eq->solve_ts_id > -1)
    cs_timer_stats_start(eq->solve_ts_id);

  cs_halo_rotation_t  halo_rota = CS_HALO_ROTATION_IGNORE;

  cs_real_t  *x = NULL;
  cs_sles_t  *sles = cs_sles_find_or_add(eq->field_id, NULL);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(eq->matrix);
  const cs_param_itsol_t  itsol_info = eq->param->itsol_info;

  printf("\n# <<ITER %5d >> Solve Ax = b for %s with %s\n"
         "# System size: %8d ; eps: % -8.5e ;\n",
         time_step->nt_cur, eq->name,
         cs_param_get_solver_name(itsol_info.solver), n_rows, itsol_info.eps);

  if (itsol_info.resid_normalized)
    r_norm = cs_euclidean_norm(n_rows, eq->rhs) / n_rows;
  else
    r_norm = 1.0;

  BFT_MALLOC(x, n_rows, cs_real_t);
  for (int i = 0; i < n_rows; i++)
    x[i] = 0.;

  cvg = cs_sles_solve(sles,
                      eq->matrix,
                      halo_rota,
                      itsol_info.eps,
                      r_norm,
                      &(ret.iter),
                      &(ret.residual),
                      eq->rhs,
                      x,
                      0,      // aux. size
                      NULL);  // aux. buffers

  if (eq->param->verbosity > 1) {
    if (time_step->nt_cur % eq->param->output_freq == 0) {
      bft_printf("  <%s/sla> code           %d\n", eq->name, cvg);
      bft_printf("  <%s/sla> n_iters        %d\n", eq->name, ret.iter);
      bft_printf("  <%s/sla> residual      % -8.4e\n", eq->name, ret.residual);
    }
  }

  printf("# n_iters = %d with a residual norm = %8.5e for %s\n",
         ret.iter, ret.residual, eq->name);

  if (eq->solve_ts_id > -1)
    cs_timer_stats_stop(eq->solve_ts_id);

  /* Store the solution in the related field structure */
  if (eq->post_ts_id > -1)
    cs_timer_stats_start(eq->post_ts_id);

  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Define the new field value for the current time */
  eq->update_field(connect, cdoq, time_step, x, eq->builder, fld->val);

  if (eq->post_ts_id > -1)
    cs_timer_stats_stop(eq->post_ts_id);

  /* Free memory */
  cs_sles_free(sles);
  BFT_FREE(x);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-processing related to this equation
 *
 * \param[in]  mesh       pointer to the mesh structure
 * \param[in]  cdoq       pointer to a cs_cdo_quantities_t struct.
 * \param[in]  time_step  pointer to a time step structure
 * \param[in]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_post(const cs_mesh_t            *mesh,
                 const cs_cdo_quantities_t  *cdoq,
                 const cs_time_step_t       *time_step,
                 const cs_equation_t        *eq)
{
  int  len;

  const int  nt_cur = time_step->nt_cur;
  char *postlabel = NULL;

  const cs_field_t  *field = cs_field_by_id(eq->field_id);
  const cs_equation_param_t  *eqp = eq->param;

  /* Cases where a post-processing is not required */
  if (eqp->post_freq == -1)
    return;
  if (nt_cur == 0) {
    if (eqp->flag & CS_EQUATION_UNSTEADY)
      return;
  }
  else { /* nt_cur > 0 */
    if (!(eqp->flag & CS_EQUATION_UNSTEADY))
      return;
    if (eqp->post_freq == 0)
      return;
    if (nt_cur % eqp->post_freq > 0)
      return;
  }
  bft_printf(" <post/var> %s\n", field->name);

  /* Perform the post-processing */
  if (eq->post_ts_id > -1)
    cs_timer_stats_start(eq->post_ts_id);

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    cs_post_write_vertex_var(-1,              // id du maillage de post
                             field->name,
                             field->dim,
                             true,            // interlace
                             true,            // true = original mesh
                             CS_POST_TYPE_cs_real_t,
                             field->val,      // values on vertices
                             time_step);      // time step management structure
    break;

  case CS_SPACE_SCHEME_CDOFB:
    {
      const cs_lnum_t  n_i_faces = mesh->n_i_faces;
      const cs_real_t  *face_pdi = cs_equation_get_face_values(eq);

      cs_post_write_var(-1,              // id du maillage de post
                        field->name,
                        field->dim,
                        field->interleaved, // interlace
                        true,               // true = original mesh
                        CS_POST_TYPE_cs_real_t,
                        field->val,         // values on cells
                        NULL,               // values at internal faces
                        NULL,               // values at border faces
                        time_step);         // time step management structure


      len = strlen(field->name) + 8 + 1;
      BFT_MALLOC(postlabel, len, char);
      sprintf(postlabel, "%s.Border", field->name);
      cs_post_write_var(-2,                    // id du maillage de post
                        postlabel,
                        field->dim,
                        field->interleaved,
                        true,                  // true = original mesh
                        CS_POST_TYPE_cs_real_t,
                        NULL,                  // values on cells
                        NULL,                  // values at internal faces
                        face_pdi + n_i_faces,  // values at border faces
                        time_step);            // time step management structure


      BFT_FREE(postlabel);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Space scheme for eq. %s is incompatible with a field.\n"
                " Stop adding a cs_field_t structure.\n"), eq->name);
    break;

  } // Switch on space_scheme

  if ( (eqp->flag & CS_EQUATION_CONVECTION) &&
       (eqp->flag & CS_EQUATION_DIFFUSION) &&
       (eqp->post_flag & CS_EQUATION_POST_PECLET ||
        eqp->post_flag & CS_EQUATION_POST_UPWIND_COEF) ) {

    cs_real_t  *work_c = NULL;
    cs_real_3_t  base_vect;

    const cs_param_advection_t  a_info = eqp->advection;
    const cs_param_hodge_t  d_info = eqp->diffusion_hodge;

    len = strlen(eq->name) + 9 + 1;
    BFT_MALLOC(postlabel, len, char);

    for (int k = 0; k < 3; k++) {

      if (k == 0) {
        sprintf(postlabel, "%s.Peclet.X", eq->name);
        base_vect[0] = 1, base_vect[1] = base_vect[2] = 0;
      }
      else if (k == 1) {
        sprintf(postlabel, "%s.Peclet.Y", eq->name);
        base_vect[1] = 1, base_vect[0] = base_vect[2] = 0;
      }
      else {
        sprintf(postlabel, "%s.Peclet.Z", eq->name);
        base_vect[2] = 1, base_vect[1] = base_vect[0] = 0;
      }

      cs_convection_get_peclet_cell(cdoq,
                                    a_info,
                                    d_info,
                                    base_vect,
                                    time_step->t_cur,
                                    &work_c);

      if (eqp->post_flag & CS_EQUATION_POST_PECLET)
        cs_post_write_var(-1,           // id du maillage de post
                          postlabel,
                          1,
                          true,         // interlace
                          true,         // true = original mesh
                          CS_POST_TYPE_cs_real_t,
                          work_c,       // values on cells
                          NULL,         // values at internal faces
                          NULL,         // values at border faces
                          time_step);   // time step management structure

      if (eqp->post_flag & CS_EQUATION_POST_UPWIND_COEF) {

        if (k == 0)
          sprintf(postlabel, "%s.UpwCoefX", eq->name);
        else if (k == 1)
          sprintf(postlabel, "%s.UpwCoefY", eq->name);
        else
          sprintf(postlabel, "%s.UpwCoefZ", eq->name);

        cs_convection_get_upwind_coef_cell(cdoq,
                                           a_info,
                                           work_c);

        cs_post_write_var(-1,           // id du maillage de post
                          postlabel,
                          1,
                          true,         // interlace
                          true,         // true = original mesh
                          CS_POST_TYPE_cs_real_t,
                          work_c,       // values on cells
                          NULL,         // values at internal faces
                          NULL,         // values at border faces
                          time_step);   // time step management structure

      } /* Post upwinding coefficient */

    } // Loop on space dimension

    BFT_FREE(postlabel);
    BFT_FREE(work_c);

  } // Post a Peclet attached to cells

  if (eq->post_ts_id > -1)
    cs_timer_stats_stop(eq->post_ts_id);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return true is the given equation is steady otherwise false
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_is_steady(const cs_equation_t    *eq)
{
  bool  is_steady = true;

  if (eq->param->flag & CS_EQUATION_UNSTEADY)
    is_steady = false;

  return is_steady;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the values of the associated field at each face of the mesh
 *         If the pointer storing the values is NULL, it is alloacted inside the
 *         function
 *
 * \param[in]       eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to the values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_face_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;

  return eq->get_f_values(eq->builder, cs_field_by_id(eq->field_id));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the name related to the given cs_equation_t structure
 *         to an equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a name or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const char *
cs_equation_get_name(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the field structure associated to a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_field_t *
cs_equation_get_field(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return cs_field_by_id(eq->field_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the cs_equation_param_t structure associated to a
 *         cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a cs_equation_param_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

const cs_equation_param_t *
cs_equation_get_param(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of numerical scheme used for the discretization in
 *         space
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  a cs_space_scheme_t variable
 */
/*----------------------------------------------------------------------------*/

cs_space_scheme_t
cs_equation_get_space_scheme(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return CS_SPACE_N_SCHEMES;
  else
    return eq->param->space_scheme;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
