/*============================================================================
 * Routines to handle specific settings related to a cs_equation_t structure
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#if defined(HAVE_PETSC)
#include <petscversion.h>
#include <petscdraw.h>
#include <petscviewer.h>
#include <petscksp.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_cdo_bc.h"
#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_multigrid.h"
#include "cs_source_term.h"
#include "cs_volume_zone.h"

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

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

static cs_param_itsol_t _itsol_info_by_default = {

  CS_PARAM_PRECOND_DIAG,  // preconditioner
  CS_PARAM_ITSOL_BICG,    // iterative solver
  2500,                   // max. number of iterations
  1e-12,                  // stopping criterion on the accuracy
  false                   // normalization of the residual (true or false)

};

static const char _err_empty_eqp[] =
  N_(" Stop setting an empty cs_equation_param_t structure.\n"
     " Please check your settings.\n");

/*============================================================================
 * Private function prototypes
 *============================================================================*/

#if defined(HAVE_PETSC)

/*----------------------------------------------------------------------------
 * \brief Add visualization of the matrix graph
 *
 * \param[in]  ksp     Krylov SubSpace structure
 *----------------------------------------------------------------------------*/

static inline void
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
 * \brief Set PETSc solver and preconditioner
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_petsc_setup_hook(void   *context,
                  KSP     ksp)
{
  PC pc;

  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_itsol_t  info = eqp->itsol_info;
  cs_param_precond_type_t  precond = info.precond;

  /* Set the solver */
  switch (info.solver) {

  case CS_PARAM_ITSOL_CG:  /* Preconditioned Conjugate Gradient */
    KSPSetType(ksp, KSPCG);
    break;

  case CS_PARAM_ITSOL_GMRES:  /* Preconditioned GMRES */
    {
      const int  n_max_restart = 30;

      KSPSetType(ksp, KSPGMRES);
      KSPGMRESSetRestart(ksp, n_max_restart);
    }
    break;

  case CS_PARAM_ITSOL_BICG: /* Preconditioned Bi-CG */
    KSPSetType(ksp, KSPBICG);
    break;

  case CS_PARAM_ITSOL_BICGSTAB2: /* Preconditioned BiCGstab2 */
    KSPSetType(ksp, KSPBCGS);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Iterative solver not interfaced with PETSc.");
  }

  /* Set the preconditioner */
  KSPGetPC(ksp, &pc);

  if (cs_glob_n_ranks > 1) {
    if (precond == CS_PARAM_PRECOND_SSOR ||
        precond == CS_PARAM_PRECOND_ILU0) {
      precond = CS_PARAM_PRECOND_BJACOB;
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    " Modify the requested preconditioner to able a parallel"
                    " computation with PETSC.\n"
                    " Switch to a block jacobi preconditioner.\n"
                    " Please check your settings.");
    }
  }


  switch (precond) {

  case CS_PARAM_PRECOND_DIAG:
    PCSetType(pc, PCJACOBI);  /* Jacobi (diagonal) preconditioning */
    break;
  case CS_PARAM_PRECOND_BJACOB:
    PCSetType(pc, PCBJACOBI);  /* Block-Jacobi (diagonal) preconditioning */
#if PETSC_VERSION_GE(3,7,0)
    PetscOptionsSetValue(NULL, "-sub_pc_factor_levels", "2");
#else
    PetscOptionsSetValue("-sub_pc_factor_levels", "2");
#endif
    break;
  case CS_PARAM_PRECOND_SSOR:
    PCSetType(pc, PCSOR);
    PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP);
    break;
  case CS_PARAM_PRECOND_ICC0:
    PCSetType(pc, PCICC);
    PCFactorSetLevels(pc, 0);
    break;
  case CS_PARAM_PRECOND_ILU0:
    PCSetType(pc, PCILU);
    PCFactorSetLevels(pc, 0);
    break;

  case CS_PARAM_PRECOND_AS:
    PCSetType(pc, PCASM);
    break;

  case CS_PARAM_PRECOND_AMG:
    {
      int  amg_type = 1;

      if (amg_type == 0) { // GAMG
#if PETSC_VERSION_GE(3,7,0)
        PetscOptionsSetValue(NULL, "-pc_gamg_agg_nsmooths", "1");
        PetscOptionsSetValue(NULL, "-mg_levels_ksp_type", "richardson");
        PetscOptionsSetValue(NULL, "-mg_levels_pc_type", "sor");
        PetscOptionsSetValue(NULL, "-mg_levels_ksp_max_it", "1");
        PetscOptionsSetValue(NULL, "-pc_gamg_threshold", "0.02");
        PetscOptionsSetValue(NULL, "-pc_gamg_reuse_interpolation", "TRUE");
        PetscOptionsSetValue(NULL, "-pc_gamg_square_graph", "4");
#else
        PetscOptionsSetValue("-pc_gamg_agg_nsmooths", "1");
        PetscOptionsSetValue("-mg_levels_ksp_type", "richardson");
        PetscOptionsSetValue("-mg_levels_pc_type", "sor");
        PetscOptionsSetValue("-mg_levels_ksp_max_it", "1");
        PetscOptionsSetValue("-pc_gamg_threshold", "0.02");
        PetscOptionsSetValue("-pc_gamg_reuse_interpolation", "TRUE");
        PetscOptionsSetValue("-pc_gamg_square_graph", "4");
#endif
        PCSetType(pc, PCGAMG);
      }
      else if (amg_type == 1) { // Boomer AMG (hypre)
#if PETSC_VERSION_GE(3,7,0)
        PetscOptionsSetValue(NULL,"-pc_type", "hypre");
        PetscOptionsSetValue(NULL,"-pc_hypre_type","boomeramg");
        PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_coarsen_type","HMIS");
        PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_interp_type","ext+i-cc");
        PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_agg_nl","2");
        PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_P_max","4");
        PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_strong_threshold","0.5");
        PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_no_CF","");
#else
        PetscOptionsSetValue("-pc_type", "hypre");
        PetscOptionsSetValue("-pc_hypre_type","boomeramg");
        PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type","HMIS");
        PetscOptionsSetValue("-pc_hypre_boomeramg_interp_type","ext+i-cc");
        PetscOptionsSetValue("-pc_hypre_boomeramg_agg_nl","2");
        PetscOptionsSetValue("-pc_hypre_boomeramg_P_max","4");
        PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold","0.5");
        PetscOptionsSetValue("-pc_hypre_boomeramg_no_CF","");
#endif

        PCSetType(pc, PCHYPRE);
      }

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Preconditioner not interfaced with PETSc.");
  }

  /* Try to have "true" norm */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

  _add_view(ksp);
}

#endif /* defined(HAVE_PETSC) */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the volume zone if from the zone name (If name = NULL or
 *         has an empty length, all entities are selected)
 *
 * \param[in] z_name            name of the zone
 *
 * \return the id of the related zone
 */
/*----------------------------------------------------------------------------*/

static inline int
_get_vzone_id(const char   *z_name)
{
  int z_id = 0;
  if (z_name != NULL) {
    if (strlen(z_name) > 0) {
      const cs_volume_zone_t  *z = cs_volume_zone_by_name(z_name);
      z_id = z->id;
    }
  }
  return z_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the boundary zone if from the zone name (If name = NULL or
 *         has an empty length, all entities are selected)
 *
 * \param[in] z_name            name of the zone
 *
 * \return the id of the related zone
 */
/*----------------------------------------------------------------------------*/

static inline int
_get_bzone_id(const char   *z_name)
{
  int z_id = 0;
  if (z_name != NULL) {
    if (strlen(z_name) > 0) {
      const cs_boundary_zone_t  *z = cs_boundary_zone_by_name(z_name);
      z_id = z->id;
    }
  }
  return z_id;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a cs_equation_param_t
 *
 * \param[in] type             type of equation
 * \param[in] dim              dim of the variable associated to this equation
 * \param[in] default_bc       type of boundary condition set by default
 *
 * \return a pointer to a new allocated cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_create_param(cs_equation_type_t     type,
                         int                    dim,
                         cs_param_bc_type_t     default_bc)
{
  cs_equation_param_t  *eqp = NULL;

  BFT_MALLOC(eqp, 1, cs_equation_param_t);

  eqp->type = type;
  eqp->dim = dim;
  eqp->verbosity = 0;
  eqp->sles_verbosity = 0;
  eqp->process_flag = 0;

  /* Build the equation flag */
  eqp->flag = 0;
  eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;
  eqp->dof_reduction = CS_PARAM_REDUCTION_DERHAM;
  eqp->space_poly_degree = 0;

  /* Vertex-based schemes imply the two following discrete Hodge operators
     Default initialization is made in accordance with this choice */
  eqp->time_hodge.is_unity = true;
  eqp->time_hodge.is_iso = true;
  eqp->time_hodge.inv_pty = false; // inverse property ?
  eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
  eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
  eqp->time_property = NULL;

  /* Description of the time discretization (default values) */
  eqp->time_scheme = CS_TIME_SCHEME_IMPLICIT;
  eqp->theta = 1.0;
  eqp->do_lumping = false;

  /* Initial condition (zero value by default) */
  eqp->n_ic_defs = 0;
  eqp->ic_defs = NULL;

  /* Diffusion term */
  eqp->diffusion_property = NULL;
  eqp->diffusion_hodge.is_unity = false;
  eqp->diffusion_hodge.is_iso = true;
  eqp->diffusion_hodge.inv_pty = false; // inverse property ?
  eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
  eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
  eqp->diffusion_hodge.coef = 1./3.; // DGA algo.

  /* Advection term */
  eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_CONSERV;
  eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
  eqp->adv_field = NULL;

  /* No reaction term by default */
  eqp->reaction_hodge.is_unity = false;
  eqp->reaction_hodge.is_iso = true;
  eqp->reaction_hodge.inv_pty = false;
  eqp->reaction_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
  eqp->reaction_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;

  eqp->n_reaction_terms = 0;
  eqp->reaction_properties = NULL;

  /* No source term by default (always in the right-hand side) */
  eqp->n_source_terms = 0;
  eqp->source_terms = NULL;

  /* Boundary conditions structure.
     One assigns a boundary condition by default */
  eqp->default_bc = default_bc;
  eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_PENA;
  eqp->n_bc_defs = 0;
  eqp->bc_defs = NULL;

  /* Settings for driving the linear algebra */
  eqp->solver_class = CS_EQUATION_SOLVER_CLASS_CS;
  eqp->itsol_info = _itsol_info_by_default;

  return eqp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_equation_param_t
 *
 * \param[in, out] eqp          pointer to a cs_equation_param_t
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_free_param(cs_equation_param_t     *eqp)
{
  if (eqp == NULL)
    return NULL;

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

  /* Information related to the definition of initial conditions */
  if (eqp->n_ic_defs > 0) {

    for (int i = 0; i < eqp->n_ic_defs; i++)
      eqp->ic_defs[i] = cs_xdef_free(eqp->ic_defs[i]);
    BFT_FREE(eqp->ic_defs);

  }

  BFT_FREE(eqp);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname in a cs_equation_param_t
 *         structure
 *
 * \param[in, out]  eqp      pointer to a cs_equation_param_t structure
 * \param[in]       key      key related to the member of eq to set
 * \param[in]       keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_param(cs_equation_param_t   *eqp,
                      cs_equation_key_t      key,
                      const char            *keyval)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  if (eqp->flag & CS_EQUATION_LOCKED)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: The current equation is not modifiable anymore.\n"
                " Please check your settings."), __func__);

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {

  case CS_EQKEY_SPACE_SCHEME:
    if (strcmp(val, "cdo_vb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;
      eqp->space_poly_degree = 0;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
    }
    else if (strcmp(val, "cdo_vcb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOVCB;
      eqp->space_poly_degree = 0;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_VC;
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CIP;
    }
    else if (strcmp(val, "cdo_fb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOFB;
      eqp->space_poly_degree = 0;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_CPVD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EDFP;
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_PENA;
    }
    else if (strcmp(val, "hho_p0") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_HHO_P0;
      eqp->space_poly_degree = 0;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_CPVD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EDFP;
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_PENA;
    }
    else if (strcmp(val, "hho_p1") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_HHO_P1;
      eqp->space_poly_degree = 1;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_CPVD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EDFP;
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_PENA;
    }
    else if (strcmp(val, "hho_p2") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_HHO_P2;
      eqp->space_poly_degree = 2;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_CPVD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EDFP;
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_PENA;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_SPACE_SCHEME\n"),
                _val);
    }
    break;

  case CS_EQKEY_DOF_REDUCTION:
    if (strcmp(val, "derham") == 0)
      eqp->dof_reduction = CS_PARAM_REDUCTION_DERHAM;
    else if (strcmp(val, "average") == 0)
      eqp->dof_reduction = CS_PARAM_REDUCTION_AVERAGE;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_DOF_REDUCTION\n"
                  " Choice between \"derham\" or \"average\"."), _val);
    }
    break;

  case CS_EQKEY_HODGE_DIFF_ALGO:
    if (strcmp(val,"cost") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
    else if (strcmp(val, "voronoi") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
    else if (strcmp(val, "wbs") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
    else if (strcmp(val, "auto") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_AUTO;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_HODGE_DIFF_ALGO\n"
                  " Choice between cost, wbs, auto or voronoi"), _val);
    }
    break;

  case CS_EQKEY_HODGE_TIME_ALGO:
    if (strcmp(val,"cost") == 0)
      eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
    else if (strcmp(val, "voronoi") == 0)
      eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
    else if (strcmp(val, "wbs") == 0)
      eqp->time_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_HODGE_TIME_ALGO\n"
                  " Choice between cost, wbs, voronoi"), _val);
    }
    break;

  case CS_EQKEY_HODGE_DIFF_COEF:
    if (strcmp(val, "dga") == 0)
      eqp->diffusion_hodge.coef = 1./3.;
    else if (strcmp(val, "sushi") == 0)
      eqp->diffusion_hodge.coef = 1./sqrt(3.);
    else if (strcmp(val, "gcr") == 0)
      eqp->diffusion_hodge.coef = 1.0;
    else
      eqp->diffusion_hodge.coef = atof(val);
    break;

  case CS_EQKEY_HODGE_TIME_COEF:
    if (strcmp(val, "dga") == 0)
      eqp->time_hodge.coef = 1./3.;
    else if (strcmp(val, "sushi") == 0)
      eqp->time_hodge.coef = 1./sqrt(3.);
    else if (strcmp(val, "gcr") == 0)
      eqp->time_hodge.coef = 1.0;
    else
      eqp->time_hodge.coef = atof(val);
    break;

  case CS_EQKEY_SOLVER_FAMILY:
    if (strcmp(val, "cs") == 0)
      eqp->solver_class = CS_EQUATION_SOLVER_CLASS_CS;
    else if (strcmp(val, "petsc") == 0)
      eqp->solver_class = CS_EQUATION_SOLVER_CLASS_PETSC;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_SOLVER_FAMILY\n"
                  " Choice between cs or petsc"), _val);
    }
    break;

  case CS_EQKEY_PRECOND:
    if (strcmp(val, "none") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_NONE;
    else if (strcmp(val, "jacobi") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_DIAG;
    else if (strcmp(val, "block_jacobi") == 0)
      eqp->itsol_info.precond = CS_PARAM_PRECOND_BJACOB;
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
                _(" Invalid val %s related to key CS_EQKEY_PRECOND\n"
                  " Choice between jacobi, block_jacobi, poly1, ssor, ilu0,\n"
                  " icc0, amg or as"), _val);
    }
    break;

  case CS_EQKEY_ITSOL:

    if (strcmp(val, "jacobi") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_JACOBI;
    else if (strcmp(val, "cg") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_CG;
    else if (strcmp(val, "bicg") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_BICG;
    else if (strcmp(val, "bicgstab2") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_BICGSTAB2;
    else if (strcmp(val, "cr3") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_CR3;
    else if (strcmp(val, "gmres") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_GMRES;
    else if (strcmp(val, "amg") == 0)
      eqp->itsol_info.solver = CS_PARAM_ITSOL_AMG;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_ITSOL\n"
                  " Choice between cg, bicg, bicgstab2, cr3, gmres or amg"),
                _val);
    }
    break;

  case CS_EQKEY_ITSOL_MAX_ITER:
    eqp->itsol_info.n_max_iter = atoi(val);
    break;

  case CS_EQKEY_ITSOL_EPS:
    eqp->itsol_info.eps = atof(val);
    break;

  case CS_EQKEY_ITSOL_RESNORM:
    if (strcmp(val, "true") == 0)
      eqp->itsol_info.resid_normalized = true;
    else
      eqp->itsol_info.resid_normalized = false;
    break;

  case CS_EQKEY_VERBOSITY: // "verbosity"
    eqp->verbosity = atoi(val);
    break;

  case CS_EQKEY_SLES_VERBOSITY: // "verbosity" for SLES structures
    eqp->sles_verbosity = atoi(val);
    break;

  case CS_EQKEY_BC_ENFORCEMENT:
    if (strcmp(val, "strong") == 0)
      eqp->enforcement = CS_PARAM_BC_ENFORCE_STRONG;
    else if (strcmp(val, "penalization") == 0)
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_PENA;
    else if (strcmp(val, "weak_sym") == 0)
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_SYM;
    else if (strcmp(val, "weak") == 0)
      eqp->enforcement = CS_PARAM_BC_ENFORCE_WEAK_NITSCHE;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value %s related to key CS_EQKEY_BC_ENFORCEMENT\n"
                  " Choice between strong, penalization, weak or weak_sym."),
                _val);
    }
    break;

  case CS_EQKEY_BC_QUADRATURE:
    {
      cs_quadrature_type_t  qtype = CS_QUADRATURE_NONE;

      if (strcmp(val, "bary") == 0)
        qtype = CS_QUADRATURE_BARY;
      else if (strcmp(val, "bary_subdiv") == 0)
        qtype = CS_QUADRATURE_BARY_SUBDIV;
      else if (strcmp(val, "higher") == 0)
        qtype = CS_QUADRATURE_HIGHER;
      else if (strcmp(val, "highest") == 0)
        qtype = CS_QUADRATURE_HIGHEST;
      else {
        const char *_val = val;
        bft_error(__FILE__, __LINE__, 0,
                  _(" Invalid value \"%s\" for key CS_EQKEY_BC_QUADRATURE\n"
                    " Valid choices are \"bary\", \"bary_subdiv\", \"higher\""
                    " and \"highest\"."), _val);
      }

      for (int i = 0; i < eqp->n_bc_defs; i++)
        cs_xdef_set_quadrature(eqp->bc_defs[i], qtype);

    }
    break;

  case CS_EQKEY_EXTRA_OP:
    if (strcmp(val, "peclet") == 0)
      eqp->process_flag |= CS_EQUATION_POST_PECLET;
    else if (strcmp(val, "upwind_coef") == 0)
      eqp->process_flag |= CS_EQUATION_POST_UPWIND_COEF;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                (" Invalid value \"%s\" for CS_EQKEY_EXTRA_OP\n"
                 " Valid keys are \"peclet\", or \"upwind_coef\"."), _val);
    }
    break;

  case CS_EQKEY_ADV_FORMULATION:
    if (strcmp(val, "conservative") == 0)
      eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_CONSERV;
    else if (strcmp(val, "non_conservative") == 0)
      eqp->adv_formulation = CS_PARAM_ADVECTION_FORM_NONCONS;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value \"%s\" for CS_EQKEY_ADV_FORMULATION\n"
                  " Valid keys are \"conservative\" or \"non_conservative\"."),
                _val);
    }
    break;

  case CS_EQKEY_ADV_SCHEME:
    if (strcmp(val, "upwind") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
    else if (strcmp(val, "samarskii") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_SAMARSKII;
    else if (strcmp(val, "sg") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_SG;
    else if (strcmp(val, "centered") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CENTERED;
    else if (strcmp(val, "cip") == 0)
      eqp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CIP;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value \"%s\" for CS_EQKEY_ADV_SCHEME\n"
                  " Valid choices are \"upwind\", \"samarskii\", \"sg\" or"
                  " \"centered\"."), _val);
    }
    break;

  case CS_EQKEY_TIME_SCHEME:
    if (strcmp(val, "implicit") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_IMPLICIT;
      eqp->theta = 1.;
    }
    else if (strcmp(val, "explicit") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_EXPLICIT;
      eqp->theta = 0.;
    }
    else if (strcmp(val, "crank_nicolson") == 0) {
      eqp->time_scheme = CS_TIME_SCHEME_CRANKNICO;
      eqp->theta = 0.5;
    }
    else if (strcmp(val, "theta_scheme") == 0)
      eqp->time_scheme = CS_TIME_SCHEME_THETA;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value \"%s\" for CS_EQKEY_TIME_SCHEME\n"
                  " Valid choices are \"implicit\", \"explicit\","
                  " \"crank_nicolson\", and \"theta_scheme\"."), _val);
    }
    break;

  case CS_EQKEY_TIME_THETA:
    eqp->theta = atof(val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting an equation."));

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set parameters for initializing SLES structures used for the
 *        resolution of the linear system.
 *        Settings are related to this equation.
 *
 * \param[in]   eqname       pointer to an cs_equation_t structure
 * \param[in]   eqp          pointer to a cs_equation_param_t struct.
 * \param[in]   field_id     id of the cs_field_t struct. for this equation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_set_sles(const char               *eqname,
                           cs_equation_param_t      *eqp,
                           int                       field_id)
{
  const cs_param_itsol_t  itsol = eqp->itsol_info;

  switch (eqp->solver_class) {
  case CS_EQUATION_SOLVER_CLASS_CS:
    {
      int  poly_degree = 0; // by default: Jacobi preconditioner

      if (itsol.precond == CS_PARAM_PRECOND_POLY1)
        poly_degree = 1;
      if (itsol.precond == CS_PARAM_PRECOND_NONE)
        poly_degree = -1;

      if (itsol.precond != CS_PARAM_PRECOND_POLY1 &&
          itsol.precond != CS_PARAM_PRECOND_DIAG &&
          itsol.precond != CS_PARAM_PRECOND_NONE)
        bft_error(__FILE__, __LINE__, 0,
                  " Incompatible preconditioner with Code_Saturne solvers.\n"
                  " Please change your settings (try PETSc ?)");

      switch (itsol.solver) { // Type of iterative solver

      case CS_PARAM_ITSOL_JACOBI:
        assert(poly_degree == -1);
        cs_sles_it_define(field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_JACOBI,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_CG:
        cs_sles_it_define(field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_PCG,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_BICG:
        cs_sles_it_define(field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_BICGSTAB,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_BICGSTAB2:
        cs_sles_it_define(field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_BICGSTAB2,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_CR3:
        cs_sles_it_define(field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_PCR3,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_GMRES:
        cs_sles_it_define(field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_GMRES,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_AMG:
        {
          cs_multigrid_t  *mg = cs_multigrid_define(field_id, NULL);

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
             -1,               // polynomial precond. degree coarse (0)
             1.0,    // precision multiplier descent (< 0 forces max iters)
             1.0,    // precision multiplier ascent (< 0 forces max iters)
             1);     // requested precision multiplier coarse (default 1)

        }
      default:
        bft_error(__FILE__, __LINE__, 0,
                  _(" Undefined iterative solver for solving %s equation.\n"
                    " Please modify your settings."), eqname);
        break;
      } // end of switch

      /* Define the level of verbosity for SLES structure */
      if (eqp->sles_verbosity > 3) {

        cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);
        cs_sles_it_t  *sles_it = (cs_sles_it_t *)cs_sles_get_context(sles);

        cs_sles_it_set_plot_options(sles_it, eqname,
                                    true);    /* use_iteration instead of
                                                 wall clock time */

      }

    } // Solver provided by Code_Saturne
    break;

  case CS_EQUATION_SOLVER_CLASS_PETSC:
    {
#if defined(HAVE_PETSC)

      /* Initialization must be called before setting options;
         it does not need to be called before calling
         cs_sles_petsc_define(), as this is handled automatically. */

      PetscBool is_initialized;
      PetscInitialized(&is_initialized);
      if (is_initialized == PETSC_FALSE) {
#if defined(HAVE_MPI)
        PETSC_COMM_WORLD = cs_glob_mpi_comm;
#endif
        PetscInitializeNoArguments();
      }

      if (eqp->itsol_info.precond == CS_PARAM_PRECOND_SSOR ||
          eqp->itsol_info.precond == CS_PARAM_PRECOND_ILU0 ||
          eqp->itsol_info.precond == CS_PARAM_PRECOND_ICC0) {

        if (cs_glob_n_ranks > 1)
          bft_error(__FILE__, __LINE__, 0,
                    " Incompatible PETSc settings for parallel run.\n");

        cs_sles_petsc_define(field_id,
                             NULL,
                             MATSEQAIJ, // Warning SEQ not MPI
                             _petsc_setup_hook,
                             (void *)eqp);

      }
      else
        cs_sles_petsc_define(field_id,
                             NULL,
                             MATMPIAIJ,
                             _petsc_setup_hook,
                             (void *)eqp);
#else
      bft_error(__FILE__, __LINE__, 0,
                _(" PETSC algorithms used to solve %s are not linked.\n"
                  " Please install Code_Saturne with PETSc."), eqname);

#endif // HAVE_PETSC
    } // Solver provided by PETSc
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Algorithm requested to solve %s is not implemented yet.\n"
                " Please modify your settings."), eqname);
    break;

  } // end switch on algorithms

  /* Define the level of verbosity for SLES structure */
  if (eqp->sles_verbosity > 1) {

    cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

    /* Set verbosity */
    cs_sles_set_verbosity(sles, eqp->sles_verbosity);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of a cs_equation_param_t structure
 *
 * \param[in]  eqname   name of the related equation
 * \param[in]  eqp      pointer to a cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_summary_param(const char                  *eqname,
                          const cs_equation_param_t   *eqp)
{
  if (eqp == NULL)
    return;

  switch (eqp->type) {
  case CS_EQUATION_TYPE_USER:
    cs_log_printf(CS_LOG_SETUP, "  <%s/type> User-defined\n", eqname);
    break;
  case CS_EQUATION_TYPE_PREDEFINED:
    cs_log_printf(CS_LOG_SETUP, "  <%s/type> Predefined\n", eqname);
    break;
  case CS_EQUATION_TYPE_GROUNDWATER:
    cs_log_printf(CS_LOG_SETUP, "  <%s/type> Associated to groundwater flows\n",
                  eqname);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              " Eq. %s has no type.\n Please check your settings.", eqname);
  }

  const char *space_scheme = cs_param_get_space_scheme_name(eqp->space_scheme);
  if (eqp->space_scheme != CS_SPACE_N_SCHEMES)
    cs_log_printf(CS_LOG_SETUP, "  <%s/space scheme> %s\n",
                  eqname, space_scheme);
  else
    bft_error(__FILE__, __LINE__, 0,
              " Undefined space scheme for eq. %s", eqname);

  cs_log_printf(CS_LOG_SETUP, "  <%s/space poly degree>  %d\n",
                eqname, eqp->space_poly_degree);

  bool  unsteady = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;
  bool  convection = (eqp->flag & CS_EQUATION_CONVECTION) ? true : false;
  bool  diffusion = (eqp->flag & CS_EQUATION_DIFFUSION) ? true : false;
  bool  reaction = (eqp->flag & CS_EQUATION_REACTION) ? true : false;
  bool  source_term = (eqp->n_source_terms > 0) ? true : false;

  cs_log_printf(CS_LOG_SETUP,
                "  <%s/Terms>  unsteady:%s, convection:%s, diffusion:%s,"
                " reaction:%s, source term:%s\n",
                eqname, cs_base_strtf(unsteady), cs_base_strtf(convection),
                cs_base_strtf(diffusion), cs_base_strtf(reaction),
                cs_base_strtf(source_term));

  /* Boundary conditions */
  if (eqp->verbosity > 0) {

    cs_log_printf(CS_LOG_SETUP,
                  "  <%s/Boundary Conditions> default: %s, enforcement: %s\n",
                  eqname, cs_param_get_bc_name(eqp->default_bc),
                  cs_param_get_bc_enforcement_name(eqp->enforcement));
    cs_log_printf(CS_LOG_SETUP, "    <%s/n_bc_definitions> %d\n",
                  eqname, eqp->n_bc_defs);
    if (eqp->verbosity > 1) {
      for (int id = 0; id < eqp->n_bc_defs; id++)
        cs_xdef_log(eqp->bc_defs[id]);
    }
  }

  if (unsteady) {

    const cs_param_hodge_t  h_info = eqp->time_hodge;

    cs_log_printf(CS_LOG_SETUP, "\n  <%s/Unsteady term>\n", eqname);
    cs_log_printf(CS_LOG_SETUP,
                  "  <%s/Initial.Condition> number of definitions %d\n",
                  eqname, eqp->n_ic_defs);

    for (int i = 0; i < eqp->n_ic_defs; i++)
      cs_xdef_log(eqp->ic_defs[i]);

    const char  *time_scheme = cs_param_get_time_scheme_name(eqp->time_scheme);
    if (time_scheme != NULL) {
      cs_log_printf(CS_LOG_SETUP, "  <%s/Time.Scheme> %s", eqname, time_scheme);
      if (eqp->time_scheme == CS_TIME_SCHEME_THETA)
        cs_log_printf(CS_LOG_SETUP, " with value %f\n", eqp->theta);
      else
        cs_log_printf(CS_LOG_SETUP, "\n");
    }
    else
      bft_error(__FILE__, __LINE__, 0, " Invalid time scheme.");

    cs_log_printf(CS_LOG_SETUP, "  <%s/Mass.Lumping> %s\n",
                  eqname, cs_base_strtf(eqp->do_lumping));
    cs_log_printf(CS_LOG_SETUP, "  <%s/Time.Property> %s\n",
                  eqname, cs_property_get_name(eqp->time_property));

    if (eqp->verbosity > 0) {
      cs_log_printf(CS_LOG_SETUP, "  <%s/Time.Hodge> %s - %s\n",
                    eqname, cs_param_hodge_get_type_name(h_info),
                    cs_param_hodge_get_algo_name(h_info));
      cs_log_printf(CS_LOG_SETUP,
                    "    <%s/Time.Hodge.Inv> Inversion of property  %s\n",
                    eqname, cs_base_strtf(h_info.inv_pty));
      if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
        cs_log_printf(CS_LOG_SETUP, "    <%s/Time.Hodge.Coef> %.3e\n",
                      eqname, h_info.coef);
    }

  } /* Unsteady term */

  if (diffusion) {

    const cs_param_hodge_t  h_info = eqp->diffusion_hodge;

    cs_log_printf(CS_LOG_SETUP, "\n  <%s/Diffusion term>\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  <%s/Diffusion.Property> %s\n",
                  eqname, cs_property_get_name(eqp->diffusion_property));

    if (eqp->verbosity > 0) {
      cs_log_printf(CS_LOG_SETUP, "  <%s/Diffusion.Hodge> %s - %s\n",
                    eqname, cs_param_hodge_get_type_name(h_info),
                    cs_param_hodge_get_algo_name(h_info));
      cs_log_printf(CS_LOG_SETUP, "    <%s/Diffusion.Hodge.Inv>", eqname);
      cs_log_printf(CS_LOG_SETUP, " Inversion of property: %s\n",
                    cs_base_strtf(h_info.inv_pty));
      if (h_info.algo == CS_PARAM_HODGE_ALGO_COST ||
          h_info.algo == CS_PARAM_HODGE_ALGO_AUTO) {
        cs_log_printf(CS_LOG_SETUP, "    <%s/Diffusion.Hodge.Coef> %.3e\n",
                      eqname, h_info.coef);
      }
    }

  } /* Diffusion term */

  if (convection) {

    cs_log_printf(CS_LOG_SETUP, "\n  <%s/Advection term>\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  <Advection field>  %s\n",
                  cs_advection_field_get_name(eqp->adv_field));

    if (eqp->verbosity > 1) {
      cs_log_printf(CS_LOG_SETUP, "  <%s/Advection.Formulation>", eqname);
      switch(eqp->adv_formulation) {
      case CS_PARAM_ADVECTION_FORM_CONSERV:
        cs_log_printf(CS_LOG_SETUP, " Conservative\n");
        break;
      case CS_PARAM_ADVECTION_FORM_NONCONS:
        cs_log_printf(CS_LOG_SETUP, " Non-conservative\n");
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid operator type for advection.");
      }

      cs_log_printf(CS_LOG_SETUP, "  <%s/Advection.Scheme> ", eqname);
      switch(eqp->adv_scheme) {
      case CS_PARAM_ADVECTION_SCHEME_CENTERED:
        cs_log_printf(CS_LOG_SETUP, " centered\n");
        break;
      case CS_PARAM_ADVECTION_SCHEME_CIP:
        cs_log_printf(CS_LOG_SETUP, " continuous interior penalty\n");
        break;
      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
        cs_log_printf(CS_LOG_SETUP, " upwind\n");
        break;
      case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
        cs_log_printf(CS_LOG_SETUP,
                      " upwind weighted with Samarskii function\n");
        break;
      case CS_PARAM_ADVECTION_SCHEME_SG:
        cs_log_printf(CS_LOG_SETUP,
                      " upwind weighted with Scharfetter-Gummel function\n");
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid weight algorithm for advection.");
      }

    } // verbosity > 0

  } // Advection term

  if (reaction) {

    cs_log_printf(CS_LOG_SETUP, "\n  <%s/Number of reaction terms> %d\n",
                  eqname, eqp->n_reaction_terms);

    if (eqp->verbosity > 0) {

      const cs_param_hodge_t  h_info = eqp->reaction_hodge;

      cs_log_printf(CS_LOG_SETUP, "  <%s/Reaction.Hodge> %s - %s\n",
                    eqname, cs_param_hodge_get_type_name(h_info),
                    cs_param_hodge_get_algo_name(h_info));
      if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
        cs_log_printf(CS_LOG_SETUP,
                      "    <%s/Reaction.Hodge.Coefficient> %.3e\n",
                      eqname, h_info.coef);

    }

  } // Reaction terms

  if (source_term) {

    cs_log_printf(CS_LOG_SETUP, "\n  <%s/Source terms>\n", eqname);
    for (int s_id = 0; s_id < eqp->n_source_terms; s_id++)
      cs_xdef_log(eqp->source_terms[s_id]);

  } // Source terms

  /* Iterative solver information */
  const cs_param_itsol_t   itsol = eqp->itsol_info;

  cs_log_printf(CS_LOG_SETUP, "\n  <%s/Sparse.Linear.Algebra>", eqname);

  if (eqp->solver_class == CS_EQUATION_SOLVER_CLASS_CS)
    cs_log_printf(CS_LOG_SETUP, " Code_Saturne iterative solvers\n");
  else if (eqp->solver_class == CS_EQUATION_SOLVER_CLASS_PETSC)
    cs_log_printf(CS_LOG_SETUP, " PETSc iterative solvers\n");

  cs_log_printf(CS_LOG_SETUP, "    <%s/sla> Solver.MaxIter     %d\n",
                eqname, itsol.n_max_iter);
  cs_log_printf(CS_LOG_SETUP, "    <%s/sla> Solver.Name        %s\n",
                eqname, cs_param_get_solver_name(itsol.solver));
  cs_log_printf(CS_LOG_SETUP, "    <%s/sla> Solver.Precond     %s\n",
                eqname, cs_param_get_precond_name(itsol.precond));
  cs_log_printf(CS_LOG_SETUP, "    <%s/sla> Solver.Eps        % -10.6e\n",
                eqname, itsol.eps);
  cs_log_printf(CS_LOG_SETUP, "    <%s/sla> Solver.Normalized  %s\n",
                eqname, cs_base_strtf(itsol.resid_normalized));

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
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_ic_by_value(cs_equation_param_t    *eqp,
                            const char             *z_name,
                            cs_real_t              *val)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int  z_id = _get_vzone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        eqp->dim,
                                        z_id,
                                        CS_FLAG_STATE_UNIFORM, // state flag
                                        meta_flag,
                                        val);

  int  new_id = eqp->n_ic_defs;
  eqp->n_ic_defs += 1;
  BFT_REALLOC(eqp->ic_defs, eqp->n_ic_defs, cs_xdef_t *);
  eqp->ic_defs[new_id] = d;
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
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_ic_by_qov(cs_equation_param_t    *eqp,
                          const char             *z_name,
                          double                  quantity)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = _get_vzone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_QOV,
                                        eqp->dim,
                                        z_id,
                                        0, // state flag
                                        meta_flag,
                                        &quantity);

  int  new_id = eqp->n_ic_defs;
  eqp->n_ic_defs += 1;
  BFT_REALLOC(eqp->ic_defs, eqp->n_ic_defs, cs_xdef_t *);
  eqp->ic_defs[new_id] = d;
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
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_ic_by_analytic(cs_equation_param_t    *eqp,
                               const char             *z_name,
                               cs_analytic_func_t     *analytic,
                               void                   *input)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = _get_vzone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_input_t  anai = {.func = analytic,
                                    .input = input };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim, z_id,
                                        0, // state flag
                                        meta_flag,
                                        &anai);

  int  new_id = eqp->n_ic_defs;
  eqp->n_ic_defs += 1;
  BFT_REALLOC(eqp->ic_defs, eqp->n_ic_defs, cs_xdef_t *);
  eqp->ic_defs[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given equation structure
 *         z_name corresponds to the name of a pre-existing cs_boundary_zone_t
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       z_name    name of the related boundary zone
 * \param[in]       values    pointer to a array storing the values
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc_by_value(cs_equation_param_t         *eqp,
                            const cs_param_bc_type_t     bc_type,
                            const char                  *z_name,
                            cs_real_t                   *values)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int  dim = eqp->dim;
  if (bc_type == CS_PARAM_BC_NEUMANN||
      bc_type == CS_PARAM_BC_HMG_NEUMANN)
    dim *= 3; // vector if scalar eq, tensor if vector eq.
  else if (bc_type == CS_PARAM_BC_ROBIN)
    dim *= 4;

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          dim,
                                          _get_bzone_id(z_name),
                                          CS_FLAG_STATE_UNIFORM, // state flag
                                          cs_cdo_bc_get_flag(bc_type), // meta
                                          (void *)values);

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to set a boundary condition
 *         related to the given equation structure
 *         z_name corresponds to the name of a pre-existing cs_boundary_zone_t
 *
 * \param[in, out]  eqp       pointer to a cs_equation_param_t structure
 * \param[in]       bc_type   type of boundary condition to add
 * \param[in]       z_name    name of the related boundary zone
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       index     optional pointer to the array index
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc_by_array(cs_equation_param_t        *eqp,
                            const cs_param_bc_type_t    bc_type,
                            const char                 *z_name,
                            cs_flag_t                   loc,
                            cs_real_t                  *array,
                            cs_lnum_t                  *index)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  assert(cs_flag_test(loc, cs_flag_primal_face) ||
         cs_flag_test(loc, cs_flag_primal_vtx));

  /* Add a new cs_xdef_t structure */
  cs_xdef_array_input_t  input = {.stride = eqp->dim,
                                  .loc = loc,
                                  .values = array,
                                  .index = index };

  cs_flag_t  state_flag = 0;
  if (loc == cs_flag_primal_face)
    state_flag = CS_FLAG_STATE_FACEWISE;

  int dim = eqp->dim;
  if (bc_type == CS_PARAM_BC_NEUMANN||
      bc_type == CS_PARAM_BC_HMG_NEUMANN)
    dim *= 3; // vector if scalar eq, tensor if vector eq.
  else if (bc_type == CS_PARAM_BC_ROBIN)
    dim *= 4;

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ARRAY,
                                          dim,
                                          _get_bzone_id(z_name),
                                          state_flag,
                                          cs_cdo_bc_get_flag(bc_type), // meta
                                          (void *)&input);

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;
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
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_bc_by_analytic(cs_equation_param_t        *eqp,
                               const cs_param_bc_type_t    bc_type,
                               const char                 *z_name,
                               cs_analytic_func_t         *analytic,
                               void                       *input)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  cs_xdef_analytic_input_t  anai = {.func = analytic,
                                    .input = input };

  int dim = eqp->dim;
  if (bc_type == CS_PARAM_BC_NEUMANN||
      bc_type == CS_PARAM_BC_HMG_NEUMANN)
    dim *= 3; // vector if scalar eq, tensor if vector eq.
  else if (bc_type == CS_PARAM_BC_ROBIN)
    dim *= 4;

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          dim,
                                          _get_bzone_id(z_name),
                                          0, // state
                                          cs_cdo_bc_get_flag(bc_type), // meta
                                          &anai);

  int  new_id = eqp->n_bc_defs;
  eqp->n_bc_defs += 1;
  BFT_REALLOC(eqp->bc_defs, eqp->n_bc_defs, cs_xdef_t *);
  eqp->bc_defs[new_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to a diffusion term
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
  assert(property != NULL);

  eqp->flag |= CS_EQUATION_DIFFUSION;
  eqp->diffusion_property = property;
  cs_property_type_t  type = cs_property_get_type(eqp->diffusion_property);
  if (type == CS_PROPERTY_ISO)
    eqp->diffusion_hodge.is_iso = true;
  else
    eqp->diffusion_hodge.is_iso = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an unsteady term
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
  assert(property != NULL);

  eqp->flag |= CS_EQUATION_UNSTEADY;
  eqp->time_property = property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to an advection term
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
  assert(adv_field != NULL);

  eqp->flag |= CS_EQUATION_CONVECTION;
  eqp->adv_field = adv_field;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to a reaction term
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
 * \brief  Define a new source term structure and initialize it by value
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or
 *                            "" all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new cs_xdef_t structure
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
  int z_id = _get_vzone_id(z_name);

  /* Define a flag according to the kind of space discretization */
  cs_flag_t  state_flag = CS_FLAG_STATE_DENSITY | CS_FLAG_STATE_UNIFORM;
  cs_flag_t  meta_flag = cs_source_term_set_default_flag(eqp->space_scheme);

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
 * \brief  Define a new source term structure and initialize it by an analytical
 *         function
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_analytic(cs_equation_param_t    *eqp,
                                        const char             *z_name,
                                        cs_analytic_func_t     *ana,
                                        void                   *input)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = _get_vzone_id(z_name);

  /* Define a flag according to the kind of space discretization */
  cs_flag_t  state_flag = CS_FLAG_STATE_DENSITY;
  cs_flag_t  meta_flag = cs_source_term_set_default_flag(eqp->space_scheme);

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_input_t  anai = {.func = ana,
                                    .input = input };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &anai);

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
 * \brief  Define a new source term defined by an array
 *
 * \param[in, out] eqp       pointer to a cs_equation_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      loc       information to know where are located values
 * \param[in]      array     pointer to an array
 * \param[in]      index     optional pointer to the array index
 *
 * \return a pointer to the new cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_equation_add_source_term_by_array(cs_equation_param_t    *eqp,
                                     const char             *z_name,
                                     cs_flag_t               loc,
                                     cs_real_t              *array,
                                     cs_lnum_t              *index)
{
  if (eqp == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: %s\n", __func__, _err_empty_eqp);

  /* Add a new cs_xdef_t structure */
  int z_id = _get_vzone_id(z_name);

  /* Define a flag according to the kind of space discretization */
  cs_flag_t  state_flag = CS_FLAG_STATE_DENSITY;
  if (cs_flag_test(loc, cs_flag_primal_cell) == true)
    state_flag |= CS_FLAG_STATE_CELLWISE;

  cs_flag_t  meta_flag = cs_source_term_set_default_flag(eqp->space_scheme);

  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_array_input_t  input = {.stride = eqp->dim,
                                  .loc = loc,
                                  .values = array,
                                  .index = index };


  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ARRAY,
                                        eqp->dim,
                                        z_id,
                                        state_flag,
                                        meta_flag,
                                        &input);

  int  new_id = eqp->n_source_terms;
  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_xdef_t *);
  eqp->source_terms[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
