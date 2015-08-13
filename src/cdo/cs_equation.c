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
#include "cs_multigrid.h"
#include "cs_timer_stats.h"
#include "cs_param.h"
#include "cs_cdovb_codits.h"
#include "cs_cdofb_codits.h"

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
 * \brief  Build a linear system within the CDO framework
 *
 * \param[in]      m        pointer to a cs_mesh_t structure
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a cs_cdo_quantities_t structure
 * \param[in]      tcur     current physical time of the simulation
 * \param[in, out] builder  pointer to builder structure
 * \param[in, out] rhs      pointer to a right-hand side array pointer
 * \param[in, out] sla_mat  pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_build_system_t)(const cs_mesh_t            *mesh,
                             const cs_cdo_connect_t     *connect,
                             const cs_cdo_quantities_t  *cdoq,
                             double                      tcur,
                             void                       *builder,
                             cs_real_t                 **rhs,
                             cs_sla_matrix_t           **sla_mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *
 * \param[in]      connect  pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant    pointer to a cs_cdo_quantities_t struct.
 * \param[in]      solu     solution array
 * \param[in, out] builder  pointer to builder structure
 * \param[in, out] field    pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_update_field_t)(const cs_cdo_connect_t       *connect,
                             const cs_cdo_quantities_t    *cdoq,
                             const cs_real_t              *solu,
                             void                         *builder,
                             cs_field_t                   *field);

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
  CS_PARAM_PRECOND_ICC0,  // preconditioner
  CS_PARAM_ITSOL_CG,      // iterative solver
#else
  CS_PARAM_PRECOND_DIAG,  // preconditioner
  CS_PARAM_ITSOL_CG,      // iterative solver
#endif
  2500,                   // max. number of iterations
  1e-12,                  // stopping criterion on the accuracy
  150,                    // output frequency
  false                   // normalization of the residual (true or false)
};

/* List of key stored within a enum structure */
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
  EQKEY_ERROR

} eqkey_t;


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

  bool   do_build;       /* false => keep the system as it is */

  /* Algebraic system */
  cs_matrix_structure_t    *ms;      /* matrix structure (how are stored
                                        coefficients of the matrix a) */
  cs_matrix_t              *matrix;  // matrix to inverse with cs_sles_solve()
  cs_real_t                *rhs;     // right-hand side

  /* System builder depending on the numerical scheme*/
  void                     *builder;

  /* Pointer to functions */
  cs_equation_init_builder_t   *init_builder;
  cs_equation_build_system_t   *build_system;
  cs_equation_update_field_t   *update_field;
  cs_equation_get_f_values_t   *get_f_values;
  cs_equation_free_builder_t   *free_builder;

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
 * PETSc solver using CG with Jacobi preconditioner
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp     <-> pointer to PETSc KSP context
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
 * PETSc solver using CG with SSOR preconditioner
 * Warning: PETSc implementation is only available in serial mode computation
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp     <-> pointer to PETSc KSP context
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
 * PETSc solver using CG with Additive Schwarz preconditioner
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp     <-> pointer to PETSc KSP context
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
 * PETSc solver using CG with ICC preconditioner
 * Warning: PETSc implementation is only available in serial mode computation
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp     <-> pointer to PETSc KSP context
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
 * PETSc solver using CG with GAMG preconditioner
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp     <-> pointer to PETSc KSP context
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
 * PETSc solver using CG with Boomer AMG preconditioner (Hypre library)
 *
 * parameters:
 *   context <-> pointer to optional (untyped) value or structure
 *   ksp     <-> pointer to PETSc KSP context
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
 * \brief  Print the name of the corresponding key
 *
 * \param[in] key        name of the key
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

static const char *
_print_key(eqkey_t  key)
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
  default:
    assert(0);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the corresponding enum from the name of a key. If not found,
 *         print an error message
 *
 * \param[in] keyname    name of the key
 *
 * \return a cs_equation_key_t
 */
/*----------------------------------------------------------------------------*/

static eqkey_t
_get_key(const char *keyname)
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

  cs_equation_param_t  *eqp = NULL;

  BFT_MALLOC(eqp, 1, cs_equation_param_t);

  eqp->status = status;
  eqp->type = type;
  eqp->verbosity = 0;

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
  eqp->unsteady_hodge.pty_id = 0;      // Unity (default property)
  eqp->unsteady_hodge.inv_pty = false; // inverse property ?
  eqp->unsteady_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
  eqp->unsteady_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;

  eqp->diffusion_hodge.pty_id = 0;      // Unity (default property)
  eqp->diffusion_hodge.inv_pty = false; // inverse property ?
  eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
  eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
  eqp->diffusion_hodge.coef = 1./3.;

  /* Boundary conditions structure.
     One assigns a boundary condition by default */
  eqp->bc = cs_param_bc_create(default_bc);

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
 * \brief  Resume a cs_equation_t structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_resume(const cs_equation_t  *eq)
{
  if (eq == NULL)
    return;

  const cs_equation_param_t  *eqp = eq->param;

  bft_printf("\n");
  bft_printf(lsepline);
  bft_printf("\tResume settings for %s eq. (variable %s)\n",
             eq->name, eq->varname);
  bft_printf(lsepline);

  bool  unsteady = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;
  bool  convection = (eqp->flag & CS_EQUATION_CONVECTION) ? true : false;
  bool  diffusion = (eqp->flag & CS_EQUATION_DIFFUSION) ? true : false;
  bool  source_term = (eqp->n_source_terms > 0) ? true : false;

  bft_printf("  <Equation>  unsteady [%s], convection [%s], diffusion [%s],"
             " source term [%s]\n",
             cs_base_strtf(unsteady), cs_base_strtf(convection),
             cs_base_strtf(diffusion), cs_base_strtf(source_term));

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB)
    bft_printf("  <Equation>  space discretization CDO vertex-based\n");
  else if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB)
    bft_printf("  <Equation>  space discretization CDO face-based\n");

  if (eqp->flag & CS_EQUATION_DIFFUSION) {

    const cs_param_hodge_t  h_info = eqp->diffusion_hodge;

    bft_printf("\n  <Equation/Diffusion term>\n");
    bft_printf("\t--> Property related to the diffusion term: %s\n",
               cs_param_pty_get_name(h_info.pty_id));

    if (eqp->verbosity > 0) {
      bft_printf("\t--> Hodge operator: %s / %s\n",
                 cs_param_hodge_get_type_name(h_info),
                 cs_param_hodge_get_algo_name(h_info));
      bft_printf("\t--> Inversion of the material property: %s\n",
                 cs_base_strtf(h_info.inv_pty));
      if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
        bft_printf("\t--> Coefficient value for COST algo: %.3e\n",
                   h_info.coef);
    }

  } // Diffusion term

  if (eqp->n_source_terms > 0) {

    bft_printf("\n  <Equation/Source terms>\n");
    for (int s_id = 0; s_id < eqp->n_source_terms; s_id++) {

      cs_param_source_term_t  st_info = eqp->source_terms[s_id];

      bft_printf("\t--> Source term name: %s\n",
                 cs_param_source_term_get_name(st_info));
      bft_printf("\t--> Related mesh location: %s\n",
                 cs_mesh_location_get_name(st_info.ml_id));
      bft_printf("\t--> Type: %s; Variable type: %s; Definition type: %s\n",
                 cs_param_source_term_get_type_name(st_info),
                 cs_param_get_var_type_name(st_info.var_type),
                 cs_param_get_def_type_name(st_info.def_type));
      if (eqp->verbosity > 0)
        bft_printf("\t--> Quadrature type: %s\n",
                   cs_quadrature_get_type_name(st_info.quad_type));

    } // Loop on source terms

  } // Source terms

  /* Iterative solver information */
  const cs_param_itsol_t   itsol = eqp->itsol_info;

  bft_printf("\n  <Equation/Iterative Solver Parameters>");
  if (eqp->algo_info.type == CS_EQUATION_ALGO_CS_ITSOL)
    bft_printf(" Code_Saturne iterative solvers\n");
  else if (eqp->algo_info.type == CS_EQUATION_ALGO_PETSC_ITSOL)
    bft_printf(" PETSc iterative solvers\n");
  bft_printf("\t-sla- Solver.MaxIter     %d\n", itsol.n_max_iter);
  bft_printf("\t-sla- Solver.Name        %s\n",
             cs_param_get_solver_name(itsol.solver));
  bft_printf("\t-sla- Solver.Precond     %s\n",
             cs_param_get_precond_name(itsol.precond));
  bft_printf("\t-sla- Solver.Eps        % -10.6e\n", itsol.eps);
  bft_printf("\t-sla- Solver.Normalized  %s\n",
             cs_base_strtf(itsol.resid_normalized));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a builder structure
 *
 * \param[in]       mesh     pointer to a cs_mesh_t structure
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_create_builder(const cs_mesh_t   *mesh,
                           cs_equation_t     *eq)
{
  if (eq == NULL)
    return;

  eq->builder = eq->init_builder(eq->param, mesh);
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
cs_equation_last_init(cs_equation_t  *eq)
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
    eq->build_system = cs_cdovb_codits_build_system;
    eq->update_field = cs_cdovb_codits_update_field;
    eq->get_f_values = NULL;
    break;

  case CS_SPACE_SCHEME_CDOFB:
    eq->init_builder = cs_cdofb_codits_init;
    eq->free_builder = cs_cdofb_codits_free;
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
                    " Couple (solver, preconditioner) not handled.");
          break;

        } // switch on PETSc reconditionner
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
  eqkey_t  key = _get_key(keyname);

  if (key == EQKEY_ERROR) {
    bft_printf("\n\n Current key: %s\n", keyname);
    bft_printf(" Possible keys: ");
    for (int i = 0; i < EQKEY_ERROR; i++) {
      bft_printf("%s ", _print_key(i));
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
      eqp->unsteady_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
    }
    else if (strcmp(val, "cdo_fb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOFB;
      eqp->unsteady_hodge.type = CS_PARAM_HODGE_TYPE_CPVD;
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
      eqp->unsteady_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
    else if (strcmp(val, "voronoi") == 0)
      eqp->unsteady_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
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
      eqp->unsteady_hodge.coef = 1./3.;
    else if (strcmp(val, "sushi") == 0)
      eqp->unsteady_hodge.coef = 1./sqrt(3.);
    else if (strcmp(val, "gcr") == 0)
      eqp->unsteady_hodge.coef = 1.0;
    else
      eqp->unsteady_hodge.coef = atof(val);
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

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key %s is not implemented yet."), keyname);

  } /* Switch on keys */

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an equation
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       pty_key   "time", "diffusion"...
 * \param[in]       pty_name  name of the material property to associate
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_pty(cs_equation_t       *eq,
                    const char          *pty_key,
                    const char          *pty_name)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" cs_equation_t structure to set with property %s is NULL\n"),
              pty_name);

  cs_equation_param_t  *eqp = eq->param;

  if (strcmp("diffusion", pty_key) == 0)
    eqp->diffusion_hodge.pty_id = cs_param_pty_get_id_by_name(pty_name);
  else if (strcmp("time", pty_key) == 0)
    eqp->unsteady_hodge.pty_id = cs_param_pty_get_id_by_name(pty_name);
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting a property.\n"
                " Current value: %s\n"
                " Possible choice: diffusion or time\n"), pty_key);
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
                " Can not add a boundary condition related to mesh location %s"),
              ml_name);

  cs_equation_param_t  *eqp = eq->param;
  cs_param_bc_t  *bc = eqp->bc;

  /* Sanity checks */
  assert(bc != NULL);

  /* Add a new definition */
  int  def_id = bc->n_defs;
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
cs_equation_add_source_term(cs_equation_t    *eq,
                            const char       *st_name,
                            const char       *ml_name,
                            const char       *st_key,
                            const char       *def_key,
                            const void       *val)
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
  // TODO

  return;
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
cs_equation_create_field(cs_equation_t    *eq)
{
  int  dim, location_id;

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
 * \brief  Check if one has to solve this equation now
 *
 * \param[in]  eq          pointer to a cs_equation_t structure
 * \param[in]  time_iter   id of the time iteration
 * \param[in]  tcur        current physical time
 *
 * \return true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_equation_needs_solve(const cs_equation_t   *eq,
                        int                    time_iter,
                        double                 tcur)
{
  if (time_iter == -1) // pre-stage after the loop in time
    if (eq->param->flag & CS_EQUATION_UNSTEADY)
      return false;

  return true;
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
 * \param[in]       m        pointer to a cs_mesh_t structure
 * \param[in]       connect  pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq     pointer to a cs_cdo_quantities_t structure
 * \param[in]       tcur     current physical time of the simulation
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_system(const cs_mesh_t            *m,
                         const cs_cdo_connect_t     *connect,
                         const cs_cdo_quantities_t  *cdoq,
                         double                      tcur,
                         cs_equation_t              *eq)
{
  cs_sla_matrix_t  *sla_mat = NULL;

  if (eq->pre_ts_id > -1)
    cs_timer_stats_start(eq->pre_ts_id);

  eq->build_system(m, connect, cdoq, tcur,
                   eq->builder,
                   &(eq->rhs),
                   &(sla_mat));

  /* Get information on the matrix related to this linear system */
  cs_sla_matrix_info_t  minfo = cs_sla_matrix_analyse(sla_mat);

  bft_printf("\n  <Linear system for equation %s>\n", eq->name);
  bft_printf(" -sla- A.size         %d\n", sla_mat->n_rows);
  bft_printf(" -sla- A.nnz          %lu\n", minfo.nnz);
  bft_printf(" -sla- A.FillIn       %5.2e %%\n", minfo.fillin);
  bft_printf(" -sla- A.StencilMin   %d\n", minfo.stencil_min);
  bft_printf(" -sla- A.StencilMax   %d\n", minfo.stencil_max);
  bft_printf(" -sla- A.StencilMean  %5.2e\n", minfo.stencil_mean);

  /* Map a cs_sla_matrix_t structure into a cs_matrix_t structure */
  assert(sla_mat->type == CS_SLA_MAT_MSR);

    /* First step: create a matrix structure */
  eq->ms =  cs_matrix_structure_create_msr(CS_MATRIX_MSR,      // type
                                           true,               // transfer
                                           true,               // have_diag
                                           sla_mat->n_rows,    // n_rows
                                           sla_mat->n_cols,    // n_cols_ext
                                           &(sla_mat->idx),    // row_index
                                           &(sla_mat->col_id), // col_id
                                           NULL,               // halo
                                           NULL);              // numbering

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

  if (eq->pre_ts_id > -1)
    cs_timer_stats_stop(eq->pre_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system for this equation
 *
 * \param[in]       connect  pointer to a cs_cdo_connect_t structure
 * \param[in]       cdoq     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(const cs_cdo_connect_t     *connect,
                  const cs_cdo_quantities_t  *cdoq,
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

  printf("\n# Solve Ax = b for %s with %s\n"
         "# System size: %8d ; eps: % -8.5e ;\n",
         eq->name, cs_param_get_solver_name(itsol_info.solver),
         n_rows, itsol_info.eps);

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

  bft_printf("\n <iterative solver convergence sumup>\n");
  bft_printf(" -sla- code        %d\n", cvg);
  bft_printf(" -sla- n_iters     %d\n", ret.iter);
  bft_printf(" -sla- residual    % -8.4e\n", ret.residual);
  printf("# n_iters = %d with a residual norm = %8.5e\n",
         ret.iter, ret.residual);

  if (eq->solve_ts_id > -1)
    cs_timer_stats_stop(eq->solve_ts_id);

  /* Store the solution in the related field structure */
  if (eq->post_ts_id > -1)
    cs_timer_stats_start(eq->post_ts_id);

  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  eq->update_field(connect, cdoq, x, eq->builder, fld);

  if (eq->post_ts_id > -1)
    cs_timer_stats_stop(eq->post_ts_id);

  /* Free memory */
  cs_sles_free(sles);
  BFT_FREE(x);
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
