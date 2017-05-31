/*============================================================================
 * Routines to handle specific settings related to a cs_equation_t structure
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_log.h"
#include "cs_mesh_location.h"
#include "cs_multigrid.h"

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
  CS_PARAM_PRECOND_DIAG,  // preconditioner
  CS_PARAM_ITSOL_BICG,      // iterative solver
  2500,                   // max. number of iterations
  1e-12,                  // stopping criterion on the accuracy
  150,                    // output frequency
  false                   // normalization of the residual (true or false)
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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a cs_equation_param_t
 *
 * \param[in] type             type of equation
 * \param[in] var_type         type of variable (scalar, vector, tensor...)
 * \param[in] default_bc       type of boundary condition set by default
 *
 * \return a pointer to a new allocated cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_equation_param_create(cs_equation_type_t     type,
                         cs_param_var_type_t    var_type,
                         cs_param_bc_type_t     default_bc)
{
  cs_equation_param_t  *eqp = NULL;

  BFT_MALLOC(eqp, 1, cs_equation_param_t);

  eqp->type = type;
  eqp->var_type = var_type;
  eqp->verbosity = 0;
  eqp->sles_verbosity = 0;
  eqp->process_flag = 0;

  /* Build the equation flag */
  eqp->flag = 0;
  eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;
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
  eqp->time_info.scheme = CS_TIME_SCHEME_IMPLICIT;
  eqp->time_info.theta = 1.0;
  eqp->time_info.do_lumping = false;

  /* Initial condition (zero value by default) */
  eqp->time_info.n_ic_definitions = 0;
  eqp->time_info.ic_definitions = NULL;

  /* Diffusion term */
  eqp->diffusion_property = NULL;
  eqp->diffusion_hodge.is_unity = false;
  eqp->diffusion_hodge.is_iso = true;
  eqp->diffusion_hodge.inv_pty = false; // inverse property ?
  eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
  eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
  eqp->diffusion_hodge.coef = 1./3.; // DGA algo.

  /* Advection term */
  eqp->advection_info.formulation = CS_PARAM_ADVECTION_FORM_CONSERV;
  eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
  eqp->advection_info.weight_criterion = CS_PARAM_ADVECTION_WEIGHT_XEXC;
  eqp->advection_info.quad_type = CS_QUADRATURE_BARY;
  eqp->advection_field = NULL;

  /* No reaction term by default */
  eqp->reaction_hodge.is_unity = false;
  eqp->reaction_hodge.is_iso = true;
  eqp->reaction_hodge.inv_pty = false;
  eqp->reaction_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
  eqp->reaction_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;

  eqp->n_reaction_terms = 0;
  eqp->reaction_info = NULL;
  eqp->reaction_properties = NULL;

  /* No source term by default (always in the right-hand side) */
  eqp->n_source_terms = 0;
  eqp->source_terms = NULL;

  /* Boundary conditions structure.
     One assigns a boundary condition by default */
  eqp->bc = cs_param_bc_create(default_bc);

  /* Settings for driving the linear algebra */
  eqp->algo_info = _algo_info_by_default;
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
cs_equation_param_free(cs_equation_param_t     *eqp)
{
  if (eqp == NULL)
    return NULL;

  if (eqp->bc != NULL) { // Boundary conditions
    cs_param_bc_t  *bc = eqp->bc;

    BFT_FREE(bc->ml_ids);
    BFT_FREE(bc->def_types);
    BFT_FREE(bc->types);
    BFT_FREE(bc->defs);

    BFT_FREE(bc);
    eqp->bc = NULL;
  }

  if (eqp->n_reaction_terms > 0) { // reaction terms

    for (int i = 0; i< eqp->n_reaction_terms; i++)
      BFT_FREE(eqp->reaction_info[i].name);

    BFT_FREE(eqp->reaction_info);
    BFT_FREE(eqp->reaction_properties);

    /* Remark: properties are freed when the global cs_domain_t structure is
       freed thanks to a call to cs_property_free() */
  }

  if (eqp->n_source_terms > 0)
    eqp->source_terms = cs_source_term_destroy(eqp->n_source_terms,
                                               eqp->source_terms);

  cs_param_time_t  t_info = eqp->time_info;
  if (t_info.n_ic_definitions > 0)
    BFT_FREE(t_info.ic_definitions);

  BFT_FREE(eqp);

  return NULL;
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
cs_equation_param_summary(const char                  *eqname,
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

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB)
    cs_log_printf(CS_LOG_SETUP,
                  "  <%s/space scheme>  CDO vertex-based\n", eqname);
  else if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB)
    cs_log_printf(CS_LOG_SETUP,
                  "  <%s/space scheme>  CDO vertex+cell-based\n", eqname);
  else if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB)
    cs_log_printf(CS_LOG_SETUP,
                  "  <%s/space scheme>  CDO face-based\n", eqname);
  else if (eqp->space_scheme == CS_SPACE_SCHEME_HHO)
    cs_log_printf(CS_LOG_SETUP,
                  "  <%s/space scheme>  HHO\n", eqname);
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
    cs_param_bc_t  *bcp = eqp->bc;

    cs_log_printf(CS_LOG_SETUP, "  <%s/Boundary Conditions>\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "    <%s/Default.BC> %s\n",
                  eqname, cs_param_get_bc_name(bcp->default_bc));
    if (eqp->verbosity > 1)
      cs_log_printf(CS_LOG_SETUP, "    <%s/BC.Enforcement> %s\n",
                    eqname, cs_param_get_bc_enforcement_name(bcp->enforcement));
    cs_log_printf(CS_LOG_SETUP, "    <%s/n_bc_definitions> %d\n",
                  eqname, bcp->n_defs);
    if (eqp->verbosity > 1) {
      for (int id = 0; id < bcp->n_defs; id++)
        cs_log_printf(CS_LOG_SETUP, "      <%s/BC.Def> Location: %s; Type: %s;"
                      " Definition type: %s\n",
                      eqname, cs_mesh_location_get_name(bcp->ml_ids[id]),
                      cs_param_get_bc_name(bcp->types[id]),
                      cs_param_get_def_type_name(bcp->def_types[id]));
    }
  }

  if (unsteady) {

    const cs_param_time_t  t_info = eqp->time_info;
    const cs_param_hodge_t  h_info = eqp->time_hodge;

    cs_log_printf(CS_LOG_SETUP, "\n  <%s/Unsteady term>\n", eqname);
    cs_log_printf(CS_LOG_SETUP,
                  "  <%s/Initial.Condition> number of definitions %d\n",
                  eqname, t_info.n_ic_definitions);
    for (int i = 0; i < t_info.n_ic_definitions; i++) {
      const cs_param_def_t  *ic = t_info.ic_definitions + i;
      cs_log_printf(CS_LOG_SETUP,
                    "    <%s/Initial.Condition> Location %s; Definition %s\n",
                    eqname, cs_mesh_location_get_name(ic->ml_id),
                    cs_param_get_def_type_name(ic->def_type));
    }
    cs_log_printf(CS_LOG_SETUP, "  <%s/Time.Scheme> ", eqname);
    switch (t_info.scheme) {
    case CS_TIME_SCHEME_IMPLICIT:
      cs_log_printf(CS_LOG_SETUP, "implicit\n");
      break;
    case CS_TIME_SCHEME_EXPLICIT:
      cs_log_printf(CS_LOG_SETUP, "explicit\n");
      break;
    case CS_TIME_SCHEME_CRANKNICO:
      cs_log_printf(CS_LOG_SETUP, "Crank-Nicolson\n");
      break;
    case CS_TIME_SCHEME_THETA:
      cs_log_printf(CS_LOG_SETUP, "theta scheme with value %f\n", t_info.theta);
      break;
    default:
      bft_error(__FILE__, __LINE__, 0, " Invalid time scheme.");
      break;
    }
    cs_log_printf(CS_LOG_SETUP, "  <%s/Mass.Lumping> %s\n",
                  eqname, cs_base_strtf(t_info.do_lumping));
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

    const cs_param_advection_t  a_info = eqp->advection_info;

    cs_log_printf(CS_LOG_SETUP, "\n  <%s/Advection term>\n", eqname);
    cs_log_printf(CS_LOG_SETUP, "  <Advection field>  %s\n",
                  cs_advection_field_get_name(eqp->advection_field));

    if (eqp->verbosity > 0) {
      cs_log_printf(CS_LOG_SETUP, "  <%s/Advection.Formulation>", eqname);
      switch(a_info.formulation) {
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
      switch(a_info.scheme) {
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

    for (int r_id = 0; r_id < eqp->n_reaction_terms; r_id++) {
      cs_param_reaction_t  r_info = eqp->reaction_info[r_id];
      cs_log_printf(CS_LOG_SETUP,
                    "    <%s/Reaction.%02d> Property: %s; Operator type: %s\n",
                    eqname, r_id,
                    cs_property_get_name(eqp->reaction_properties[r_id]),
                    cs_param_reaction_get_type_name(r_info.type));
    }

  } // Reaction terms

  if (source_term) {

    cs_log_printf(CS_LOG_SETUP, "\n  <%s/Source terms>\n", eqname);
    for (int s_id = 0; s_id < eqp->n_source_terms; s_id++)
      cs_source_term_summary(eqname, eqp->source_terms + s_id);

  } // Source terms

  /* Iterative solver information */
  const cs_param_itsol_t   itsol = eqp->itsol_info;

  cs_log_printf(CS_LOG_SETUP, "\n  <%s/Sparse.Linear.Algebra>", eqname);
  if (eqp->algo_info.type == CS_EQUATION_ALGO_CS_ITSOL)
    cs_log_printf(CS_LOG_SETUP, " Code_Saturne iterative solvers\n");
  else if (eqp->algo_info.type == CS_EQUATION_ALGO_PETSC_ITSOL)
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
 * \brief Initialize SLES structure for the resolution of the linear system
 *        according to the settings related to this equation
 *
 * \param[in]   eqname       pointer to an cs_equation_t structure
 * \param[in]   eqp          pointer to a cs_equation_param_t struct.
 * \param[in]   field_id     id of the cs_field_t struct. for this equation
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_param_init_sles(const char                 *eqname,
                            const cs_equation_param_t  *eqp,
                            int                         field_id)
{
  const cs_equation_algo_t  algo = eqp->algo_info;
  const cs_param_itsol_t  itsol = eqp->itsol_info;

  switch (algo.type) {
  case CS_EQUATION_ALGO_CS_ITSOL:
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
             0,                // polynomial precond. degree coarse (0)
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

  case CS_EQUATION_ALGO_PETSC_ITSOL:
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
    cs_sles_it_t  *sles_it = (cs_sles_it_t *)cs_sles_get_context(sles);

    /* Set verbosity */
    cs_sles_set_verbosity(sles, eqp->sles_verbosity);

  }

}
