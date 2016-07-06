/*============================================================================
 * Routines to handle cs_equation_t structure and its related structures
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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
#include "cs_cdovb_scaleq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_cdofb_scaleq.h"
#include "cs_evaluate.h"
#include "cs_mesh_location.h"
#include "cs_multigrid.h"
#include "cs_timer_stats.h"

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
 * \param[in] eq         pointer to a cs_equation_param_t structure
 * \param[in] mesh       pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated builder structure
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_init_builder_t)(const cs_equation_param_t  *eqp,
                             const cs_mesh_t            *mesh);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of source terms for the current time
 *
 * \param[in, out] builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_compute_source_t)(void          *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system within the CDO framework
 *
 * \param[in]      m          pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in, out] builder    pointer to builder structure
 * \param[in, out] rhs        pointer to a right-hand side array pointer
 * \param[in, out] sla_mat    pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_build_system_t)(const cs_mesh_t            *mesh,
                             const cs_real_t            *field_val,
                             double                      dt_cur,
                             void                       *builder,
                             cs_real_t                 **rhs,
                             cs_sla_matrix_t           **sla_mat);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in, out] builder    pointer to builder structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_update_field_t)(const cs_real_t            *solu,
                             const cs_real_t            *rhs,
                             void                       *builder,
                             cs_real_t                  *field_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Extra-operation related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in, out]  builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_extra_op_t)(const char                 *eqname,
                         const cs_field_t           *field,
                         void                       *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at a different location than that of the
 *         field associated to this equation
 *
 * \param[in]  builder    pointer to a builder structure
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

typedef double *
(cs_equation_get_extra_values_t)(const void          *builder);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *
 * \param[in]       builder    pointer to a builder structure
 * \param[in]       f_vals     pointer to an array of field values
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in]       direction  indicate in which direction flux is > 0
 * \param[in, out]  diff_flux  pointer to the value of the diffusive flux
 * \param[in, out]  conv_flux  pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_get_fap_t)(const void            *builder,
                        const cs_real_t       *f_vals,
                        int                    ml_id,
                        const cs_real_t        direction[],
                        double                *diff_flux,
                        double                *conv_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all faces.
 *         Primal or dual faces are considered according to the space scheme
 *
 * \param[in]       builder    pointer to a builder structure
 * \param[in]       f_vals     pointer to an array of field values
 * \param[in, out]  diff_flux  pointer to the value of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_get_cw_dflux_t)(const cs_real_t    *f_vals,
                             void               *builder,
                             cs_real_t          *d_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to a buffer of size at least the number of unknows
 *
 * \param[in, out]  builder    pointer to a builder structure
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

typedef cs_real_t *
(cs_equation_get_tmpbuf_t)(void);

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_sla_matrix_t related to the system to solve
 *
 * \param[in, out]  builder   pointer to a builder structure
 * \param[in, out]  matrix    pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_free_sysmat_t)(void              *builder,
                            cs_sla_matrix_t   *matrix);

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

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

struct _cs_equation_t {

  char *restrict         name;    /* Short description */

  cs_equation_param_t   *param;   /* Set of parameters related to an equation */

  /* Variable attached to this equation is defined as a cs_field_t structure */
  char *restrict         varname;
  int                    field_id;

  /* Timer statistic for a "light" profiling */
  int     main_ts_id;   /* Id of the main timer states structure related
                           to this equation */
  int     pre_ts_id;    /* Id of the timer stats structure gathering all
                           steps before the resolution of the linear systems */
  int     solve_ts_id;  /* Id of the timer stats structure related to the
                           inversion of the linear system */
  int     post_ts_id;   /* Id of the timer stats structure gathering all steps
                           after the resolution of the linear systems
                           (post, balance...) */

  bool    do_build;     /* false => keep the system as it is */

  /* Algebraic system */
  cs_matrix_structure_t    *ms;      /* matrix structure (how are stored
                                        coefficients of the matrix) */
  cs_matrix_t              *matrix;  // matrix to inverse with cs_sles_solve()
  cs_real_t                *rhs;     // right-hand side

  /* System builder depending on the numerical scheme*/
  void                     *builder;

  /* Pointer to functions (see prototypes just above) */
  cs_equation_init_builder_t      *init_builder;
  cs_equation_free_builder_t      *free_builder;
  cs_equation_build_system_t      *build_system;
  cs_equation_free_sysmat_t       *free_system_matrix;
  cs_equation_update_field_t      *update_field;
  cs_equation_compute_source_t    *compute_source;
  cs_equation_get_fap_t           *compute_flux_across_plane;
  cs_equation_get_cw_dflux_t      *compute_cellwise_diff_flux;
  cs_equation_extra_op_t          *postprocess;
  cs_equation_get_extra_values_t  *get_extra_values;
  cs_equation_get_tmpbuf_t        *get_tmpbuf;

};

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_eq[] =
  N_(" Stop setting an empty cs_equation_t structure.\n"
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
      bft_printf(" Modify the requested preconditioner to able a parallel"
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
    PetscOptionsSetValue("-sub_pc_factor_levels", "2");
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
        PetscOptionsSetValue("-pc_gamg_agg_nsmooths", "1");
        PetscOptionsSetValue("-mg_levels_ksp_type", "richardson");
        PetscOptionsSetValue("-mg_levels_pc_type", "sor");
        PetscOptionsSetValue("-mg_levels_ksp_max_it", "1");
        PetscOptionsSetValue("-pc_gamg_threshold", "0.02");
        PetscOptionsSetValue("-pc_gamg_reuse_interpolation", "TRUE");
        PetscOptionsSetValue("-pc_gamg_square_graph", "4");

        PCSetType(pc, PCGAMG);
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

/*----------------------------------------------------------------------------
 * \brief Initialize SLES structure for the resolution of the linear system
 *
 * \param[in]      eq       pointer to an cs_equation_t structure
 *----------------------------------------------------------------------------*/

static void
_sles_initialization(const cs_equation_t  *eq)
{
  const cs_equation_param_t  *eqp = eq->param;
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
                          CS_SLES_BICGSTAB,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_BICGSTAB2:
        cs_sles_it_define(eq->field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_BICGSTAB2,
                          poly_degree,
                          itsol.n_max_iter);
        break;
      case CS_PARAM_ITSOL_CR3:
        cs_sles_it_define(eq->field_id,  // give the field id (future: eq_id ?)
                          NULL,
                          CS_SLES_PCR3,
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

      /* Define the level of verbosity for SLES structure */
      int  sles_verbosity = eq->param->sles_verbosity;
      if (sles_verbosity > 1) {

        cs_sles_t  *sles = cs_sles_find_or_add(eq->field_id, NULL);
        cs_sles_it_t  *sles_it = (cs_sles_it_t *)cs_sles_get_context(sles);

        /* Set verbosity */
        cs_sles_set_verbosity(sles, sles_verbosity);

        if (sles_verbosity > 2) /* Add plot */
          cs_sles_it_set_plot_options(sles_it,
                                      eq->name,
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
          eqp->itsol_info.precond == CS_PARAM_PRECOND_ICC0)
        cs_sles_petsc_define(eq->field_id,
                             NULL,
                             MATSEQAIJ, // Warning SEQ not MPI
                             _petsc_setup_hook,
                             (void *)eqp);
      else
        cs_sles_petsc_define(eq->field_id,
                             NULL,
                             MATMPIAIJ,
                             _petsc_setup_hook,
                             (void *)eqp);
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

}

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
  if (ml_name == NULL)
    bft_error(__FILE__, __LINE__, 0, _(" Mesh location name is NULL."));

  *p_ml_id = cs_mesh_location_get_id_by_name(ml_name);

  if (*p_ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid mesh location name %s.\n"
                " This mesh location is not already defined.\n"), ml_name);
}

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

static cs_equation_param_t *
_create_equation_param(cs_equation_type_t     type,
                       cs_param_var_type_t    var_type,
                       cs_param_bc_type_t     default_bc)
{
  cs_equation_param_t  *eqp = NULL;

  BFT_MALLOC(eqp, 1, cs_equation_param_t);

  eqp->type = type;
  eqp->var_type = var_type;
  eqp->verbosity =  0;
  eqp->sles_verbosity =  0;
  eqp->process_flag =  0;

  /* Build the equation flag */
  eqp->flag = 0;
  eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;

  /* Vertex-based schemes imply the two following discrete Hodge operators
     Default initialization is made in accordance with this choice */
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
 * \brief  Set the initial values for the variable related to an equation
 *
 * \param[in, out]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_initialize_field_from_ic(cs_equation_t     *eq)
{
  int  ml_id = -1;
  cs_flag_t  dof_flag = 0;
  cs_equation_param_t  *eqp = eq->param;

  switch (eqp->var_type) {
  case CS_PARAM_VAR_SCAL:
    dof_flag |= CS_FLAG_SCAL;
    break;
  case CS_PARAM_VAR_VECT:
    dof_flag |= CS_FLAG_VECT;
    break;
  case CS_PARAM_VAR_TENS:
    dof_flag |= CS_FLAG_TENS;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Incompatible type of variable for equation %s."), eq->name);
    break;
  }

  /* Retrieve the associated field */
  cs_field_t  *field = cs_field_by_id(eq->field_id);
  cs_param_time_t  t_info = eqp->time_info;

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB ||
      eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB) {

    dof_flag |= cs_cdo_primal_vtx;

    for (int def_id = 0; def_id < t_info.n_ic_definitions; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_param_def_t  *ic = t_info.ic_definitions + def_id;

      if (strlen(ic->ml_name) > 0)
        ml_id = cs_mesh_location_get_id_by_name(ic->ml_name);
      else
        ml_id = cs_mesh_location_get_id_by_name(N_("vertices"));

      if (ic->def_type == CS_PARAM_DEF_BY_VALUE)
        cs_evaluate_potential_from_value(dof_flag, ml_id, ic->def.get,
                                         field->val);
      else if (ic->def_type == CS_PARAM_DEF_BY_QOV)
        cs_evaluate_potential_from_qov(dof_flag, ml_id, ic->def.get,
                                       field->val);
      else if (ic->def_type == CS_PARAM_DEF_BY_ANALYTIC_FUNCTION)
        cs_evaluate_potential_from_analytic(dof_flag, ml_id, ic->def.analytic,
                                            field->val);

    } // Loop on definitions

  } // VB or VCB schemes --> initialize on vertices

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB) {

    cs_real_t  *face_values = eq->get_extra_values(eq->builder);
    assert(face_values != NULL);

    for (int def_id = 0; def_id < t_info.n_ic_definitions; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_param_def_t  *ic = t_info.ic_definitions + def_id;

      if (strlen(ic->ml_name) > 0)
        ml_id = cs_mesh_location_get_id_by_name(ic->ml_name);
      else
        ml_id = cs_mesh_location_get_id_by_name(N_("cells"));

      /* Initialize cell-based array */
      cs_flag_t  face_flag = dof_flag | cs_cdo_primal_face;

      if (ic->def_type == CS_PARAM_DEF_BY_VALUE)
        cs_evaluate_potential_from_value(face_flag, ml_id, ic->def.get,
                                         face_values);

      else if (ic->def_type == CS_PARAM_DEF_BY_ANALYTIC_FUNCTION)
        cs_evaluate_potential_from_analytic(face_flag, ml_id, ic->def.analytic,
                                            face_values);

    } // Loop on definitions

  } // FB schemes --> initialize on faces

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB ||
      eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB) {

    /* Initialize cell-based array */
    cs_flag_t  cell_flag = dof_flag | cs_cdo_primal_cell;

    cs_real_t  *cell_values = field->val;
    if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB)
      cell_values = eq->get_extra_values(eq->builder);
    assert(cell_values != NULL);

    for (int def_id = 0; def_id < t_info.n_ic_definitions; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_param_def_t  *ic = t_info.ic_definitions + def_id;

      if (strlen(ic->ml_name) > 0)
        ml_id = cs_mesh_location_get_id_by_name(ic->ml_name);
      else
        ml_id = cs_mesh_location_get_id_by_name(N_("cells"));

      if (ic->def_type == CS_PARAM_DEF_BY_VALUE)
        cs_evaluate_potential_from_value(cell_flag, ml_id, ic->def.get,
                                         cell_values);

      else if (ic->def_type == CS_PARAM_DEF_BY_ANALYTIC_FUNCTION)
        cs_evaluate_potential_from_analytic(cell_flag, ml_id, ic->def.analytic,
                                            cell_values);

    } // Loop on definitions

  } // FB or VCB schemes --> initialize on cells

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
 * \param[in] eqtype           type of equation (user, predefined...)
 * \param[in] vartype          type of variable (scalar, vector, tensor...)
 * \param[in] default_bc       type of boundary condition set by default
 *
 * \return  a pointer to the new allocated cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_t *
cs_equation_create(const char            *eqname,
                   const char            *varname,
                   cs_equation_type_t     eqtype,
                   cs_param_var_type_t    vartype,
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

  eq->param = _create_equation_param(eqtype, vartype, default_bc);

  eq->field_id = -1;    // field is created in a second step
  eq->do_build = true;  // Force the construction of the algebraic system

  /* Set timer statistic structure to a default value */
  eq->main_ts_id = eq->pre_ts_id = eq->solve_ts_id = eq->post_ts_id = -1;

  /* Algebraic system: allocated later */
  eq->ms = NULL;
  eq->matrix = NULL;
  eq->rhs = NULL;

  /* Builder structure for this equation */
  eq->builder = NULL;

  /* Pointer of function */
  eq->init_builder = NULL;
  eq->compute_source = NULL;
  eq->build_system = NULL;
  eq->update_field = NULL;
  eq->postprocess = NULL;
  eq->get_extra_values = NULL;
  eq->get_tmpbuf = NULL;
  eq->free_builder = NULL;

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
    cs_timer_stats_start(eq->main_ts_id);

  BFT_FREE(eq->name);
  BFT_FREE(eq->varname);

  cs_equation_param_t  *eqp = eq->param;

  if (eqp->bc != NULL) { // Boundary conditions
    if (eqp->bc->n_defs > 0)
      BFT_FREE(eqp->bc->defs);
    BFT_FREE(eqp->bc);
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

  if (eqp->n_source_terms > 0) { // Source terms

    for (int i = 0; i< eqp->n_source_terms; i++)
      eqp->source_terms[i] = cs_source_term_free(eqp->source_terms[i]);
    BFT_FREE(eqp->source_terms);

  }

  cs_param_time_t  t_info = eqp->time_info;
  if (t_info.n_ic_definitions > 0) {
    for (int i = 0; i < t_info.n_ic_definitions; i++) {
      cs_param_def_t  *ic = t_info.ic_definitions + i;
      BFT_FREE(ic->ml_name);
    }
    BFT_FREE(t_info.ic_definitions);
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

  switch (eqp->type) {
  case CS_EQUATION_TYPE_USER:
    bft_printf("  <%s/type> User-defined\n", eq->name);
    break;
  case CS_EQUATION_TYPE_PREDEFINED:
    bft_printf("  <%s/type> Predefined\n", eq->name);
    break;
  case CS_EQUATION_TYPE_GROUNDWATER:
    bft_printf("  <%s/type> Associated to groundwater flows\n", eq->name);
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              " Eq. %s has no type.\n"
              " Please check your settings.", eq->name);
  }

  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB)
    bft_printf("  <%s/space scheme>  CDO vertex-based\n", eq->name);
  else if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVCB)
    bft_printf("  <%s/space scheme>  CDO vertex+cell-based\n", eq->name);
  else if (eqp->space_scheme == CS_SPACE_SCHEME_CDOFB)
    bft_printf("  <%s/space scheme>  CDO face-based\n", eq->name);

  bool  unsteady = (eqp->flag & CS_EQUATION_UNSTEADY) ? true : false;
  bool  convection = (eqp->flag & CS_EQUATION_CONVECTION) ? true : false;
  bool  diffusion = (eqp->flag & CS_EQUATION_DIFFUSION) ? true : false;
  bool  reaction = (eqp->flag & CS_EQUATION_REACTION) ? true : false;
  bool  source_term = (eqp->n_source_terms > 0) ? true : false;

  bft_printf("  <%s/Terms>  unsteady:%s, convection:%s, diffusion:%s,"
             " reaction:%s, source term:%s\n",
             eq->name, cs_base_strtf(unsteady), cs_base_strtf(convection),
             cs_base_strtf(diffusion), cs_base_strtf(reaction),
             cs_base_strtf(source_term));

  /* Boundary conditions */
  if (eqp->verbosity > 0) {
    cs_param_bc_t  *bcp = eqp->bc;

    bft_printf("  <%s/Boundary Conditions>\n", eq->name);
    bft_printf("\t<BC/Default> %s\n", cs_param_get_bc_name(bcp->default_bc));
    if (eqp->verbosity > 1)
      bft_printf("\t<BC/Enforcement> %s\n",
                 cs_param_get_bc_enforcement_name(bcp->enforcement));
    bft_printf("\t<BC/N_Definitions> %d\n", bcp->n_defs);
    if (eqp->verbosity > 1) {
      for (int id = 0; id < bcp->n_defs; id++)
        bft_printf("\t\t<BC> Location: %s; Type: %s; Definition type: %s\n",
                   cs_mesh_location_get_name(bcp->defs[id].loc_id),
                   cs_param_get_bc_name(bcp->defs[id].bc_type),
                   cs_param_get_def_type_name(bcp->defs[id].def_type));
    }
  }

  if (unsteady) {

    const cs_param_time_t  t_info = eqp->time_info;
    const cs_param_hodge_t  h_info = eqp->time_hodge;

    bft_printf("\n  <%s/Unsteady term>\n", eq->name);
    bft_printf("  <Time/Initial condition> number of definitions %d\n",
               t_info.n_ic_definitions);
    for (int i = 0; i < t_info.n_ic_definitions; i++) {
      const cs_param_def_t  *ic = t_info.ic_definitions + i;
      bft_printf("\t<Time/Initial condition> Location %s; Definition %s\n",
                 ic->ml_name, cs_param_get_def_type_name(ic->def_type));
    }
    bft_printf("  <Time/Scheme> ");
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
    bft_printf("  <Time/Mass lumping> %s\n", cs_base_strtf(t_info.do_lumping));
    bft_printf("  <Time/Property> %s\n",
               cs_property_get_name(eqp->time_property));

    if (eqp->verbosity > 0) {
      bft_printf("  <Time/Hodge> %s - %s\n",
                 cs_param_hodge_get_type_name(h_info),
                 cs_param_hodge_get_algo_name(h_info));
      bft_printf("\t<Time/Hodge> Inversion of property: %s\n",
                 cs_base_strtf(h_info.inv_pty));
      if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
        bft_printf("\t<Time/Hodge> Value of the coercivity coef.: %.3e\n",
                   h_info.coef);
    }

  } /* Unsteady term */

  if (diffusion) {

    const cs_param_hodge_t  h_info = eqp->diffusion_hodge;

    bft_printf("\n  <%s/Diffusion term>\n", eq->name);
    bft_printf("  <Diffusion> Property: %s\n",
               cs_property_get_name(eqp->diffusion_property));

    if (eqp->verbosity > 0) {
      bft_printf("  <Diffusion/Hodge> %s - %s\n",
                 cs_param_hodge_get_type_name(h_info),
                 cs_param_hodge_get_algo_name(h_info));
      bft_printf("\t<Diffusion/Hodge> Inversion of property: %s\n",
                 cs_base_strtf(h_info.inv_pty));
      if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
        bft_printf("\t<Diffusion/Hodge> Value of the coercivity coef.: %.3e\n",
                   h_info.coef);
    }

  } /* Diffusion term */

  if (convection) {

    const cs_param_advection_t  a_info = eqp->advection_info;

    bft_printf("\n  <%s/Advection term>\n", eq->name);
    bft_printf("  <Advection field>  %s\n",
               cs_advection_field_get_name(eqp->advection_field));

    if (eqp->verbosity > 0) {
      bft_printf("  <Advection/Formulation>");
      switch(a_info.formulation) {
      case CS_PARAM_ADVECTION_FORM_CONSERV:
        bft_printf(" Conservative\n");
        break;
      case CS_PARAM_ADVECTION_FORM_NONCONS:
        bft_printf(" Non-conservative\n");
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid operator type for advection.");
      }

      bft_printf("  <Advection/Scheme> ");
      switch(a_info.scheme) {
      case CS_PARAM_ADVECTION_SCHEME_CENTERED:
        bft_printf(" centered\n");
        break;
      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
        bft_printf(" upwind\n");
        break;
      case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
        bft_printf(" upwind weighted with Samarskii function\n");
        break;
      case CS_PARAM_ADVECTION_SCHEME_SG:
        bft_printf(" upwind weighted with Scharfetter-Gummel function\n");
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid weight algorithm for advection.");
      }

    } // verbosity > 0

  } // Advection term

  if (reaction) {

    bft_printf("\n  <%s/Number of reaction terms> %d\n",
               eq->name, eqp->n_reaction_terms);

    if (eqp->verbosity > 0) {

      const cs_param_hodge_t  h_info = eqp->reaction_hodge;

      bft_printf("  <Reaction/Hodge> %s - %s\n",
                 cs_param_hodge_get_type_name(h_info),
                 cs_param_hodge_get_algo_name(h_info));
      if (h_info.algo == CS_PARAM_HODGE_ALGO_COST)
        bft_printf("\t<Reaction/Hodge> Value of the coercivity coef.: %.3e\n",
                   h_info.coef);

    }

    for (int r_id = 0; r_id < eqp->n_reaction_terms; r_id++) {
      cs_param_reaction_t  r_info = eqp->reaction_info[r_id];
      bft_printf("  <Reaction tem %02d> Property: %s; Operator type: %s\n",
                 r_id, cs_property_get_name(eqp->reaction_properties[r_id]),
                 cs_param_reaction_get_type_name(r_info.type));

    }

  } // Reaction terms

  if (source_term) {

    bft_printf("\n  <%s/Source terms>\n", eq->name);
    for (int s_id = 0; s_id < eqp->n_source_terms; s_id++)
      cs_source_term_summary(eq->name, eqp->source_terms[s_id]);

  } // Source terms

  /* Iterative solver information */
  const cs_param_itsol_t   itsol = eqp->itsol_info;

  bft_printf("\n  <%s/Sparse Linear Algebra>", eq->name);
  if (eqp->algo_info.type == CS_EQUATION_ALGO_CS_ITSOL)
    bft_printf(" Code_Saturne iterative solvers\n");
  else if (eqp->algo_info.type == CS_EQUATION_ALGO_PETSC_ITSOL)
    bft_printf(" PETSc iterative solvers\n");
  bft_printf("\t<sla> Solver.MaxIter     %d\n", itsol.n_max_iter);
  bft_printf("\t<sla> Solver.Name        %s\n",
             cs_param_get_solver_name(itsol.solver));
  bft_printf("\t<sla> Solver.Precond     %s\n",
             cs_param_get_precond_name(itsol.precond));
  bft_printf("\t<sla> Solver.Eps        % -10.6e\n", itsol.eps);
  bft_printf("\t<sla> Solver.Normalized  %s\n",
             cs_base_strtf(itsol.resid_normalized));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create timer statistics structures to enable a "home-made" profiling
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_timer_stats(cs_equation_t  *eq)
{
  if (eq == NULL)
    return;

  cs_equation_param_t  *eqp = eq->param;

  /* Set timer statistics */
  if (eqp->verbosity > 0) {

    eq->main_ts_id = cs_timer_stats_create(NULL, // new root
                                           eq->name,
                                           eq->name);

    cs_timer_stats_start(eq->main_ts_id);

    if (eqp->verbosity > 1) {

      char *label = NULL;

      int  len = strlen("_solve") + strlen(eq->name) + 1;
      BFT_MALLOC(label, len, char);
      sprintf(label, "%s_pre", eq->name);
      eq->pre_ts_id = cs_timer_stats_create(eq->name, label, label);

      sprintf(label, "%s_solve", eq->name);
      eq->solve_ts_id = cs_timer_stats_create(eq->name, label, label);

      sprintf(label, "%s_post", eq->name);
      eq->post_ts_id = cs_timer_stats_create(eq->name, label, label);

      BFT_FREE(label);

    } // verbosity > 1

  } // verbosity > 0

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

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_equation_param_t  *eqp = eq->param;

  /* Set function pointers */
  switch(eqp->space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    eq->init_builder = cs_cdovb_scaleq_init;
    eq->free_builder = cs_cdovb_scaleq_free;
    eq->build_system = cs_cdovb_scaleq_build_system;
    eq->free_system_matrix = cs_cdovb_scaleq_free_sysmat;
    eq->compute_source = cs_cdovb_scaleq_compute_source;
    eq->update_field = cs_cdovb_scaleq_update_field;
    eq->postprocess = cs_cdovb_scaleq_extra_op;
    eq->get_tmpbuf = cs_cdovb_scaleq_get_tmpbuf;
    eq->get_extra_values = NULL;
    eq->compute_flux_across_plane = cs_cdovb_scaleq_compute_flux_across_plane;
    eq->compute_cellwise_diff_flux = cs_cdovb_scaleq_compute_cw_diff_flux;
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    eq->init_builder = cs_cdovcb_scaleq_init;
    eq->free_builder = cs_cdovcb_scaleq_free;
    eq->build_system = cs_cdovcb_scaleq_build_system;
    eq->free_system_matrix = cs_cdovcb_scaleq_free_sysmat;
    eq->compute_source = cs_cdovcb_scaleq_compute_source;
    eq->update_field = cs_cdovcb_scaleq_update_field;
    eq->postprocess = cs_cdovcb_scaleq_extra_op;
    eq->get_tmpbuf = cs_cdovcb_scaleq_get_tmpbuf;
    eq->get_extra_values = cs_cdovcb_scaleq_get_cell_values;
    eq->compute_flux_across_plane = cs_cdovcb_scaleq_compute_flux_across_plane;
    eq->compute_cellwise_diff_flux = cs_cdovcb_scaleq_compute_cw_diff_flux;
    break;

  case CS_SPACE_SCHEME_CDOFB:
    eq->init_builder = cs_cdofb_scaleq_init;
    eq->free_builder = cs_cdofb_scaleq_free;
    eq->build_system = cs_cdofb_scaleq_build_system;
    eq->free_system_matrix = cs_cdofb_scaleq_free_sysmat;
    eq->compute_source = cs_cdofb_scaleq_compute_source;
    eq->update_field = cs_cdofb_scaleq_update_field;
    eq->postprocess = cs_cdofb_scaleq_extra_op;
    eq->get_tmpbuf = cs_cdofb_scaleq_get_tmpbuf;
    eq->get_extra_values = cs_cdofb_scaleq_get_face_values;
    eq->compute_flux_across_plane = NULL;
    eq->compute_cellwise_diff_flux = NULL;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid scheme for the space discretization.\n"
                " Please check your settings."));
    break;
  }

  /* Initialize cs_sles_t structure */
  _sles_initialization(eq);

  /* Flag this equation such that parametrization is not modifiable anymore */
  eqp->flag |= CS_EQUATION_LOCKED;

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter in a cs_equation_t structure attached to keyname
 *
 * \param[in, out]  eq       pointer to a cs_equation_t structure
 * \param[in]       key      key related to the member of eq to set
 * \param[in]       keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_param(cs_equation_t       *eq,
                      cs_equation_key_t    key,
                      const char          *keyval)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_equation_param_t  *eqp = eq->param;

  if (eqp->flag & CS_EQUATION_LOCKED)
    bft_error(__FILE__, __LINE__, 0,
              _(" Equation %s is not modifiable anymore.\n"
                " Please check your settings."), eq->name);

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {

  case CS_EQKEY_SPACE_SCHEME:
    if (strcmp(val, "cdo_vb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOVB;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
      eqp->diffusion_hodge.type = CS_PARAM_HODGE_TYPE_EPFD;
    }
    if (strcmp(val, "cdo_vcb") == 0) {
      eqp->space_scheme = CS_SPACE_SCHEME_CDOVCB;
      eqp->time_hodge.type = CS_PARAM_HODGE_TYPE_VPCD;
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
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
                _(" Invalid val %s related to key CS_EQKEY_SPACE_SCHEME\n"
                  " Choice between cdo_vb or cdo_fb"), _val);
    }
    break;

  case CS_EQKEY_HODGE_DIFF_ALGO:
    if (strcmp(val,"cost") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_COST;
    else if (strcmp(val, "voronoi") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_VORONOI;
    else if (strcmp(val, "wbs") == 0)
      eqp->diffusion_hodge.algo = CS_PARAM_HODGE_ALGO_WBS;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_HODGE_DIFF_ALGO\n"
                  " Choice between cost, wbs or voronoi"), _val);
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
                  " Choice between cost, wbs or voronoi"), _val);
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
      eqp->algo_info.type = CS_EQUATION_ALGO_CS_ITSOL;
    else if (strcmp(val, "petsc") == 0)
      eqp->algo_info.type = CS_EQUATION_ALGO_PETSC_ITSOL;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid val %s related to key CS_EQKEY_SOLVER_FAMILY\n"
                  " Choice between cs or petsc"), _val);
    }
    break;

  case CS_EQKEY_PRECOND:
    if (strcmp(val, "jacobi") == 0)
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
    if (strcmp(val, "cg") == 0)
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

  case CS_EQKEY_ITSOL_VERBOSITY: // "verbosity" for SLES structures
    eqp->sles_verbosity = atoi(val);
    break;

  case CS_EQKEY_BC_ENFORCEMENT:
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
                _(" Invalid value %s related to key CS_EQKEY_BC_ENFORCEMENT\n"
                  " Choice between strong, penalization, weak or weak_sym."),
                _val);
    }
    break;

  case CS_EQKEY_BC_QUADRATURE:
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
                _(" Invalid value \"%s\" for key CS_EQKEY_BC_QUADRATURE\n"
                  " Valid choices are \"subdiv\", \"bary\", \"higher\" and"
                  " \"highest\"."), _val);
    }
    break;

  case CS_EQKEY_EXTRA_OP:
    if (strcmp(val, "peclet") == 0)
      eqp->process_flag |= CS_EQUATION_POST_PECLET;
    else if (strcmp(val, "upwind_coef") == 0)
      eqp->process_flag |= CS_EQUATION_POST_UPWIND_COEF;
    else
      eqp->process_flag |= CS_EQUATION_POST_NONE;
    break;

  case CS_EQKEY_ADV_FORMULATION:
    if (strcmp(val, "conservative") == 0)
      eqp->advection_info.formulation = CS_PARAM_ADVECTION_FORM_CONSERV;
    else if (strcmp(val, "non_conservative") == 0)
      eqp->advection_info.formulation = CS_PARAM_ADVECTION_FORM_NONCONS;
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
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
    else if (strcmp(val, "samarskii") == 0)
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_SAMARSKII;
    else if (strcmp(val, "sg") == 0)
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_SG;
    else if (strcmp(val, "centered") == 0)
      eqp->advection_info.scheme = CS_PARAM_ADVECTION_SCHEME_CENTERED;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value \"%s\" for CS_EQKEY_ADV_SCHEME\n"
                  " Valid choices are \"upwind\", \"samarskii\", \"sg\" or"
                  " \"centered\"."), _val);
    }
    break;

  case CS_EQKEY_ADV_FLUX_QUADRA:
    if (strcmp(val, "bary") == 0)
      eqp->advection_info.quad_type = CS_QUADRATURE_BARY;
    else if (strcmp(val, "higher") == 0)
      eqp->advection_info.quad_type = CS_QUADRATURE_HIGHER;
    else if (strcmp(val, "highest") == 0)
      eqp->advection_info.quad_type = CS_QUADRATURE_HIGHEST;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid value \"%s\" for CS_EQKEY_ADV_FLUX_QUADRA\n"
                  " Valid choices are \"bary\", \"higher\" and \"highest\"."),
                _val);
    }
    break;

  case CS_EQKEY_TIME_SCHEME:
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
                _(" Invalid value \"%s\" for CS_EQKEY_TIME_SCHEME\n"
                  " Valid choices are \"implicit\", \"explicit\","
                  " \"crank_nicolson\", and \"theta_scheme\"."), _val);
    }
    break;

  case CS_EQKEY_TIME_THETA:
    eqp->time_info.theta = atof(val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting an equation."));

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
 * \param[in]       keyword   "time", "diffusion", "advection"
 * \param[in]       pointer   pointer to a given structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_link(cs_equation_t       *eq,
                 const char          *keyword,
                 void                *pointer)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  cs_equation_param_t  *eqp = eq->param;

  if (strcmp("diffusion", keyword) == 0) {

    eqp->flag |= CS_EQUATION_DIFFUSION;
    eqp->diffusion_property = (cs_property_t *)pointer;

  }
  else if (strcmp("time", keyword) == 0) {

    eqp->flag |= CS_EQUATION_UNSTEADY;
    eqp->time_property = (cs_property_t *)pointer;

  }
  else if (strcmp("advection", keyword) == 0) {

    eqp->flag |= CS_EQUATION_CONVECTION;
    eqp->advection_field = (cs_adv_field_t *)pointer;

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid keyword for linking an equation.\n"
                " Current value: \"%s\"\n"
                " Valid choices: \"diffusion\", \"time\", \"advection\"."),
              keyword);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition of the unknown related to this equation
 *         This definition can be done mesh location by mesh location
 *         Available types of definition are: "value" and "analytic"
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       ml_name   name of the associated mesh location (if NULL or
 *                            "" all entities are considered)
 * \param[in]       def_key   way of defining the value of the BC
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_ic(cs_equation_t    *eq,
                   const char       *ml_name,
                   const char       *def_key,
                   void             *val)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  cs_equation_param_t  *eqp = eq->param;
  cs_param_time_t  t_info = eqp->time_info;
  int  id = t_info.n_ic_definitions;

  BFT_REALLOC(t_info.ic_definitions, id+1, cs_param_def_t);

  cs_param_def_t  *ic = t_info.ic_definitions + id;

  /* Get the type of definition */
  if (strcmp(def_key, "value") == 0)
    ic->def_type = CS_PARAM_DEF_BY_VALUE;
  else if (strcmp(def_key, "quantity_over_volume") == 0)
    ic->def_type = CS_PARAM_DEF_BY_QOV;
  else if (strcmp(def_key, "analytic") == 0)
    ic->def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the initial condition.\n"
                " Given key: \"%s\"\n"
                " Valid choices are \"value\" and \"analytic\".\n"
                " Please modify your settings."), def_key);

  /* Handle the name of the mesh location */
  if (ml_name == NULL) {
    BFT_MALLOC(ic->ml_name, 1, char);
    strcpy(ic->ml_name, "");
  }
  else {
    BFT_MALLOC(ic->ml_name, strlen(ml_name) + 1, char);
    strcpy(ic->ml_name, ml_name);
  }

  /* Set the definition */
  cs_param_set_def(ic->def_type, eqp->var_type, val, &(ic->def));

  /* Update the number of definitions */
  t_info.n_ic_definitions += 1;
  eqp->time_info = t_info;
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

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  cs_equation_param_t  *eqp = eq->param;
  cs_param_bc_t  *bc = eqp->bc;
  cs_param_bc_type_t  bc_type = CS_PARAM_N_BC_TYPES;
  cs_param_def_type_t  def_type = CS_PARAM_N_DEF_TYPES;

  /* Sanity checks */
  assert(bc != NULL);

  /* Add a new definition */
  int def_id = bc->n_defs;
  bc->n_defs += 1;
  BFT_REALLOC(bc->defs, bc->n_defs, cs_param_bc_def_t);

  /* Get the mesh location id from its name */
  _check_ml_name(ml_name, &ml_id);

  /* Get the type of definition */
  if (strcmp(def_key, "value") == 0)
    def_type = CS_PARAM_DEF_BY_VALUE;
  else if (strcmp(def_key, "array") == 0)
    def_type = CS_PARAM_DEF_BY_ARRAY;
  else if (strcmp(def_key, "analytic") == 0)
    def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  else if (strcmp(def_key, "user") == 0)
    def_type = CS_PARAM_DEF_BY_USER_FUNCTION;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the type of definition.\n"
                " Given key: \"%s\"\n"
                " Valid keys: \"analytic\", \"user\", \"law\" or \"array\"\n"
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
                " Given key: \"%s\"\n"
                " Valid keys: \"dirichlet\", \"neumann\" or \"robin\".\n"
                " Please modify your settings."), bc_key);

  /* Check if this is a homogeneous boundary condition */
  if (def_type == CS_PARAM_DEF_BY_VALUE && eqp->var_type == CS_PARAM_VAR_SCAL) {
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
                      eqp->var_type,
                      def_type,
                      val, NULL); // coef2 is not used up to now
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to a reaction term
 *
 * \param[in, out] eq         pointer to a cs_equation_t structure
 * \param[in]      property   pointer to a cs_property_t struct.
 * \param[in]      r_name     name of the reaction term (optional, i.e. NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_linear_reaction(cs_equation_t   *eq,
                                cs_property_t   *property,
                                const char      *r_name)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  /* Only this kind of reaction term is available up to now */
  cs_equation_param_t  *eqp = eq->param;

  int  r_id = eqp->n_reaction_terms;

  /* Add a new source term */
  eqp->n_reaction_terms += 1;
  BFT_REALLOC(eqp->reaction_info, eqp->n_reaction_terms, cs_param_reaction_t);

  /* Set the type of reaction term to consider */
  eqp->reaction_info[r_id].type = CS_PARAM_REACTION_TYPE_LINEAR;

  /* Associate a name to this reaction term */
  if (r_name == NULL) { /* Define a name by default */

    assert(r_id < 100);
    int  len = strlen("reaction_00") + 1;
    BFT_MALLOC(eqp->reaction_info[r_id].name, len, char);
    sprintf(eqp->reaction_info[r_id].name, "reaction_%02d", r_id);

  }
  else { /* Copy the given name */

    int  len = strlen(r_name) + 1;
    BFT_MALLOC(eqp->reaction_info[r_id].name, len, char);
    strncpy(eqp->reaction_info[r_id].name, r_name, len);

  }

  /* Associate a property to this reaction term */
  BFT_REALLOC(eqp->reaction_properties, eqp->n_reaction_terms, cs_property_t *);
  eqp->reaction_properties[r_id] = property;

  /* Flag the equation with "reaction" */
  eqp->flag |= CS_EQUATION_REACTION;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize a new structure to store parameters related
 *         to a source term
 *         def_key among "value", "analytic", "user"...
 *
 * \param[in, out] eq            pointer to a cs_equation_t structure
 * \param[in]      ml_id         id related to a mesh location
 * \param[in]      array_desc    short description of this array (mask of bits)
 * \param[in]      array_values  pointer to the array values
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_gravity_source_term(cs_equation_t   *eq,
                                    int              ml_id,
                                    cs_desc_t        array_desc,
                                    cs_real_t       *array_values)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  cs_source_term_type_t  st_type = CS_SOURCE_TERM_GRAVITY;
  cs_equation_param_t  *eqp = eq->param;

  /* Add a new source term */
  int  st_id = eqp->n_source_terms;

  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_source_term_t *);

  /* Create and set new source term structure */
  switch (eqp->space_scheme) {
  case CS_SPACE_SCHEME_CDOVB:
    eqp->source_terms[st_id] = cs_source_term_create("gravity_source",
                                                     ml_id,
                                                     st_type,
                                                     CS_SOURCE_TERM_REDUC_DUAL,
                                                     eqp->var_type);
    break;

  case CS_SPACE_SCHEME_CDOFB:
    eqp->source_terms[st_id] = cs_source_term_create("gravity_source",
                                                     ml_id,
                                                     st_type,
                                                     CS_SOURCE_TERM_REDUC_PRIM,
                                                     eqp->var_type);
    break;

  case CS_SPACE_SCHEME_CDOVCB:
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid numerical scheme to set a source term."));

  }

  cs_source_term_def_by_array(eqp->source_terms[st_id],
                              array_desc,
                              array_values);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize by value a new structure to store parameters
 *         related to a source term defined by a user
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       st_name   name of the source term or NULL
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       val       pointer to the value
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_source_term_by_val(cs_equation_t   *eq,
                                   const char      *st_name,
                                   const char      *ml_name,
                                   const void      *val)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  int  ml_id;
  char *_st_name = NULL;
  cs_source_term_type_t  st_type = CS_SOURCE_TERM_USER;
  cs_equation_param_t  *eqp = eq->param;

  const char  *name;

  /* Add a new source term */
  int  st_id = eqp->n_source_terms;

  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_source_term_t *);

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

  /* Create and set new source term structure */
  switch (eqp->space_scheme) {
  case CS_SPACE_SCHEME_CDOVB:
    eqp->source_terms[st_id] = cs_source_term_create(name,
                                                     ml_id,
                                                     st_type,
                                                     CS_SOURCE_TERM_REDUC_DUAL,
                                                     eqp->var_type);
    break;

  case CS_SPACE_SCHEME_CDOFB:
    eqp->source_terms[st_id] = cs_source_term_create(name,
                                                     ml_id,
                                                     st_type,
                                                     CS_SOURCE_TERM_REDUC_PRIM,
                                                     eqp->var_type);
    break;

  case CS_SPACE_SCHEME_CDOVCB: // TODO
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid numerical scheme to set a source term."));

  }

  cs_source_term_def_by_value(eqp->source_terms[st_id], val);

  if (st_name == NULL)
    BFT_FREE(_st_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define and initialize by an analytical function a new structure
 *         related to a source term defined by a user
 *
 * \param[in, out]  eq        pointer to a cs_equation_t structure
 * \param[in]       st_name   name of the source term or NULL
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       ana       pointer to an analytical function
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_add_source_term_by_analytic(cs_equation_t        *eq,
                                        const char           *st_name,
                                        const char           *ml_name,
                                        cs_analytic_func_t   *ana)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  int  ml_id;
  char *_st_name = NULL;
  cs_source_term_type_t  st_type = CS_SOURCE_TERM_USER;
  cs_equation_param_t  *eqp = eq->param;

  const char  *name;

  /* Add a new source term */
  int  st_id = eqp->n_source_terms;

  eqp->n_source_terms += 1;
  BFT_REALLOC(eqp->source_terms, eqp->n_source_terms, cs_source_term_t *);

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

  cs_source_term_reduction_t  default_red = CS_SOURCE_TERM_REDUC_PRIM;
  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB)
    default_red = CS_SOURCE_TERM_REDUC_DUAL;

  /* Create and set new source term structure */
  eqp->source_terms[st_id] = cs_source_term_create(name,
                                                   ml_id,
                                                   st_type,
                                                   default_red,
                                                   eqp->var_type);

  cs_source_term_def_by_analytic(eqp->source_terms[st_id], ana);

  if (st_name == NULL)
    BFT_FREE(_st_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of quadrature to use for computing a source term
 *         If st_name is NULL, all source terms of the given equation are set
 *
 * \param[in, out]  eq         pointer to a cs_equation_t structure
 * \param[in]       st_name    name of the source term
 * \param[in]       quad_type  type of quadrature to use
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_source_term_quadrature(cs_equation_t      *eq,
                                       const char         *st_name,
                                       cs_quadra_type_t    quad_type)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_equation_param_t  *eqp = eq->param;

  /* Look for the requested source term */
  int  st_id = -1;
  if (st_name != NULL) { // Look for the related source term structure

    for (int i = 0; i < eqp->n_source_terms; i++) {
      if (strcmp(cs_source_term_get_name(eqp->source_terms[i]), st_name) == 0) {
        st_id = i;
        break;
      }
    }

    if (st_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _(" Cannot find source term %s among defined source terms.\n"
                  " Please check your settings for equation %s.\n"),
                st_name, eq->name);

  } // st_name != NULL

  if (st_id != -1)
    cs_source_term_set_quadrature(eqp->source_terms[st_id], quad_type);
  else
    for (int i = 0; i < eqp->n_source_terms; i++)
      cs_source_term_set_quadrature(eqp->source_terms[i], quad_type);

  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of quadrature to use for computing a source term
 *         If st_name is NULL, all source terms of the given equation are set
 *
 * \param[in, out]  eq         pointer to a cs_equation_t structure
 * \param[in]       st_name    name of the source term
 * \param[in]       type       type of reduction to apply
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_set_source_term_reduction(cs_equation_t               *eq,
                                      const char                  *st_name,
                                      cs_source_term_reduction_t   type)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);

  cs_equation_param_t  *eqp = eq->param;

  /* Look for the requested source term */
  int  st_id = -1;
  if (st_name != NULL) { // Look for the related source term structure

    for (int i = 0; i < eqp->n_source_terms; i++) {
      if (strcmp(cs_source_term_get_name(eqp->source_terms[i]), st_name) == 0) {
        st_id = i;
        break;
      }
    }

    if (st_id == -1)
      bft_error(__FILE__, __LINE__, 0,
                _(" Cannot find source term %s among defined source terms.\n"
                  " Please check your settings for equation %s.\n"),
                st_name, eq->name);

  } // st_name != NULL

  if (st_id != -1)
    cs_source_term_set_reduction(eqp->source_terms[st_id], type);
  else
    for (int i = 0; i < eqp->n_source_terms; i++)
      cs_source_term_set_reduction(eqp->source_terms[i], type);

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
  switch (eqp->var_type) {
  case CS_PARAM_VAR_SCAL:
    dim = 1;
    break;
  case CS_PARAM_VAR_VECT:
    dim = 3;
    break;
  case CS_PARAM_VAR_TENS:
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
  case CS_SPACE_SCHEME_CDOVCB:
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
                                     has_previous);

  /* Set default value for default keys */
  cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
  cs_field_set_key_int(fld, cs_field_key_id("post_vis"), 1);

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
 * \param[in, out]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_init_system(const cs_mesh_t            *mesh,
                        cs_equation_t              *eq)
{
  if (eq == NULL)
    return;

  if (eq->main_ts_id > -1) {
    cs_timer_stats_start(eq->main_ts_id);
    if (eq->pre_ts_id > -1)
      cs_timer_stats_start(eq->pre_ts_id);
  }

  const cs_equation_param_t  *eqp = eq->param;

  /* Allocate and initialize a system builder */
  eq->builder = eq->init_builder(eqp, mesh);

  /* Compute the (initial) source term */
  eq->compute_source(eq->builder);

  /* Initialize the associated field to the initial condition if unsteady */
  if (!(eqp->flag & CS_EQUATION_UNSTEADY)) {
    if (eq->main_ts_id > -1) {
      if (eq->pre_ts_id > -1)
        cs_timer_stats_stop(eq->pre_ts_id);
      cs_timer_stats_stop(eq->main_ts_id);
    }
    return;
  }

  cs_param_time_t  t_info = eqp->time_info;

  if (t_info.n_ic_definitions == 0) {
    if (eq->main_ts_id > -1) {
      if (eq->pre_ts_id > -1)
        cs_timer_stats_stop(eq->pre_ts_id);
      cs_timer_stats_stop(eq->main_ts_id);
    }
    return; // By default, 0 is set
  }

  _initialize_field_from_ic(eq);

  if (eq->main_ts_id > -1) {
    cs_timer_stats_stop(eq->main_ts_id);
    if (eq->pre_ts_id > -1)
      cs_timer_stats_stop(eq->pre_ts_id);
  }
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
 * \param[in]       m           pointer to a cs_mesh_t structure
 * \param[in]       time_step   pointer to a time step structure
 * \param[in]       dt_cur      value of the current time step
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_build_system(const cs_mesh_t            *mesh,
                         const cs_time_step_t       *time_step,
                         double                      dt_cur,
                         cs_equation_t              *eq)
{
  cs_sla_matrix_t  *sla_mat = NULL;

  const char *eqn = eq->name;
  const cs_equation_param_t  *eqp = eq->param;
  const cs_field_t  *fld = cs_field_by_id(eq->field_id);

  if (eq->main_ts_id > -1) {
    cs_timer_stats_start(eq->main_ts_id);
    if (eq->pre_ts_id > -1)
      cs_timer_stats_start(eq->pre_ts_id);
  }

  eq->build_system(mesh, fld->val, dt_cur,
                   eq->builder,
                   &(eq->rhs),
                   &(sla_mat));

  /* Get information on the matrix related to this linear system */
  if (eqp->verbosity > 1 && time_step->nt_cur == 0) {

    cs_sla_matrix_set_info(sla_mat);

    cs_sla_matrix_info_t  minfo = sla_mat->info;

    bft_printf("\n Sparse Linear Algebra (SLA) sumup:\n");
    bft_printf("  <%s/sla> A.size         %d\n", eqn, sla_mat->n_rows);
    bft_printf("  <%s/sla> A.nnz          %lu\n", eqn, minfo.nnz);
    bft_printf("  <%s/sla> A.FillIn       %5.2e %%\n", eqn, minfo.fillin);
    bft_printf("  <%s/sla> A.StencilMin   %d\n", eqn, minfo.stencil_min);
    bft_printf("  <%s/sla> A.StencilMax   %d\n", eqn, minfo.stencil_max);
    bft_printf("  <%s/sla> A.StencilMean  %5.2e\n", eqn, minfo.stencil_mean);
  }

  /* Map a cs_sla_matrix_t structure into a cs_matrix_t structure */
  assert(sla_mat->type == CS_SLA_MAT_MSR);

  bool  do_idx_transfer = false;
  if (eqp->space_scheme == CS_SPACE_SCHEME_CDOVB)
    if (eqp->bc->enforcement == CS_PARAM_BC_ENFORCE_STRONG)
      do_idx_transfer = true;

  /* First step: create a matrix structure */
  if (eq->ms == NULL)
    eq->ms = cs_matrix_structure_create_msr(CS_MATRIX_MSR,      // type
                                            do_idx_transfer,    // transfer
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
  eq->free_system_matrix(eq->builder, sla_mat);
  eq->do_build = false;

  if (eq->main_ts_id > -1) {
    if (eq->pre_ts_id > -1)
      cs_timer_stats_stop(eq->pre_ts_id);
    cs_timer_stats_stop(eq->main_ts_id);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system for this equation
 *
 * \param[in, out]  eq          pointer to a cs_equation_t structure
 * \param[in]       do_logcvg   output information on convergence or not
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_solve(cs_equation_t   *eq,
                  bool             do_logcvg)
{
  int  n_iters = 0;
  double  r_norm = DBL_MAX, residual = DBL_MAX;

  if (eq->main_ts_id > -1)
    cs_timer_stats_start(eq->main_ts_id);
  if (eq->solve_ts_id > -1)
    cs_timer_stats_start(eq->solve_ts_id);

  cs_real_t  *x = eq->get_tmpbuf();
  cs_sles_t  *sles = cs_sles_find_or_add(eq->field_id, NULL);
  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(eq->matrix);
  const cs_param_itsol_t  itsol_info = eq->param->itsol_info;

  if (eq->param->sles_verbosity > 0)
    printf("\n# %s >> Solve Ax = b with %s as solver and %s as precond.\n"
           "# System size: %8d ; eps: % -8.5e ;\n",
           eq->name, cs_param_get_solver_name(itsol_info.solver),
           cs_param_get_precond_name(itsol_info.precond),
           n_rows, itsol_info.eps);

  if (itsol_info.resid_normalized)
    r_norm = cs_euclidean_norm(n_rows, eq->rhs) / n_rows;
  else
    r_norm = 1.0;

  /* Sanity check (up to now, only scalar field are handled) */
  assert(fld->dim == 1);
  for (cs_lnum_t  i = 0; i < n_rows; i++)
    x[i] = fld->val[i];

  cs_sles_convergence_state_t code = cs_sles_solve(sles,
                                                   eq->matrix,
                                                   CS_HALO_ROTATION_IGNORE,
                                                   itsol_info.eps,
                                                   r_norm,
                                                   &n_iters,
                                                   &residual,
                                                   eq->rhs,
                                                   x,
                                                   0,      // aux. size
                                                   NULL);  // aux. buffers

  if (do_logcvg)
    bft_printf("  <%s/sles_cvg> code  %d n_iters  %d residual  % -8.4e\n",
               eq->name, code, n_iters, residual);

  if (eq->param->sles_verbosity > 0)
    printf("# %s >> n_iters = %d with a residual norm = %8.5e\n",
           eq->name, n_iters, residual);

  if (eq->solve_ts_id > -1)
    cs_timer_stats_stop(eq->solve_ts_id);
  if (eq->post_ts_id > -1)
    cs_timer_stats_start(eq->post_ts_id);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Define the new field value for the current time */
  eq->update_field(x, eq->rhs, eq->builder, fld->val);

  if (eq->post_ts_id > -1)
    cs_timer_stats_stop(eq->post_ts_id);
  if (eq->main_ts_id > -1)
    cs_timer_stats_stop(eq->main_ts_id);

  if (eq->param->flag & CS_EQUATION_UNSTEADY)
    eq->do_build = true; /* Improvement: exhibit cases where a new build
                            is not needed */
  /* Free memory */
  cs_sles_free(sles);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]  eq         pointer to a cs_equation_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_equation_extra_op(const cs_equation_t        *eq)
{
  if (eq == NULL)
    return;

  const cs_field_t  *field = cs_field_by_id(eq->field_id);
  const cs_equation_param_t  *eqp = eq->param;

  /* Cases where a post-processing is not required */
  if (eqp->process_flag & CS_EQUATION_POST_NONE)
    return;

  if (eq->main_ts_id > -1) {
    cs_timer_stats_start(eq->main_ts_id);
    if (eq->solve_ts_id > -1)
      cs_timer_stats_start(eq->solve_ts_id);
  }

  /* Perform the post-processing */
  eq->postprocess(eq->name, field, eq->builder);

  if (eq->main_ts_id > -1) {
    if (eq->post_ts_id > -1)
      cs_timer_stats_stop(eq->post_ts_id);
    cs_timer_stats_stop(eq->main_ts_id);
  }
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
 * \brief  Return the flag associated to an equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a flag (cs_flag_t type)
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_equation_get_flag(const cs_equation_t    *eq)
{
  cs_flag_t  ret_flag = 0;

  if (eq == NULL)
    return ret_flag;

  ret_flag = eq->param->flag;

  return ret_flag;
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
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         diffusion term for this equation (NULL if not activated).
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_diffusion_property(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->param->diffusion_property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         unsteady term for this equation (NULL if not activated).
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_time_property(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  else
    return eq->param->time_property;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a pointer to the cs_property_t structure associated to the
 *         reaction term called r_name and related to this equation
 *
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return a pointer to a cs_property_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_equation_get_reaction_property(const cs_equation_t    *eq,
                                  const char             *r_name)
{
  if (eq == NULL)
    return NULL;

  if (r_name == NULL)
    return NULL;

  const cs_equation_param_t  *eqp = eq->param;

  /* Look for the requested reaction term */
  int  r_id = -1;
  for (int i = 0; i < eqp->n_reaction_terms; i++) {
    if (strcmp(eqp->reaction_info[i].name, r_name) == 0) {
      r_id = i;
      break;
    }
  }

  if (r_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" Cannot find the reaction term %s in equation %s.\n"
                " Please check your settings.\n"), r_name, eq->name);

  return eqp->reaction_properties[r_id];
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
/*!
 * \brief  Return the type of variable solved by this equation
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the type of variable (sclar, vector...) associated to this equation
 */
/*----------------------------------------------------------------------------*/

cs_param_var_type_t
cs_equation_get_var_type(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return CS_PARAM_N_VAR_TYPES;
  else
    return eq->param->var_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return the type of equation for the given equation structure
 *
 * \param[in]  eq       pointer to a cs_equation_t structure
 *
 * \return  the type of the given equation
 */
/*----------------------------------------------------------------------------*/

cs_equation_type_t
cs_equation_get_type(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return CS_EQUATION_N_TYPES;
  else
    return eq->param->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each face of the mesh for the field unknowns
 *         related to this equation.
 *
 * \param[in]   eq        pointer to a cs_equation_t structure
 *
 * \return a pointer to the face values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_face_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  if (eq->get_extra_values == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" No function defined for getting the face values in eq. %s"),
              eq->name);

  if (eq->param->space_scheme == CS_SPACE_SCHEME_CDOFB)
    return eq->get_extra_values(eq->builder);
  else
    return NULL; // Not implemented
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at each cell centers for the field unknowns
 *         related to this equation.
 *
 * \param[in]   eq        pointer to a cs_equation_t structure
 *
 * \return a pointer to the cell values
 */
/*----------------------------------------------------------------------------*/

const cs_real_t *
cs_equation_get_cell_values(const cs_equation_t    *eq)
{
  if (eq == NULL)
    return NULL;
  if (eq->get_extra_values == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" No function defined for getting the cell values in eq. %s"),
              eq->name);

  switch (eq->param->space_scheme) {
  case CS_SPACE_SCHEME_CDOFB:
    {
      cs_field_t  *fld = cs_field_by_id(eq->field_id);

      return fld->val;
    }

  case CS_SPACE_SCHEME_CDOVCB:
    return eq->get_extra_values(eq->builder);

  default:
    return NULL; // Not implemented
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux accross a plane defined
 *         by a mesh location structure attached to the name ml_name.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in]      ml_name     name of the related mesh location
 * \param[in]      direction   vector indicating in which direction flux is > 0
 * \param[in, out] diff_flux   value of the diffusive part of the flux
 * \param[in, out] conv_flux   value of the convective part of the flux
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_flux_across_plane(const cs_equation_t   *eq,
                                      const char            *ml_name,
                                      const cs_real_3_t      direction,
                                      cs_real_t             *diff_flux,
                                      cs_real_t             *conv_flux)
{
  int  ml_id;

  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);
  if (eq->compute_flux_across_plane == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Computation of the diffusive and convective flux across\n"
                " a plane is not available for equation %s\n"), eq->name);

  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  /* Get the mesh location id from its name */
  _check_ml_name(ml_name, &ml_id);

  /* Perform the computation */
  eq->compute_flux_across_plane(eq->builder,
                                fld->val,
                                ml_id,
                                direction,
                                diff_flux, conv_flux);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all faces.
 *         Primal or dual faces are considered according to the space scheme.
 *
 * \param[in]      eq          pointer to a cs_equation_t structure
 * \param[in, out] diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_equation_compute_diff_flux(const cs_equation_t   *eq,
                              cs_real_t             *diff_flux)
{
  if (eq == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_eq);
  if (eq->compute_cellwise_diff_flux == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Cellwise computation of the diffusive flux is not\n"
                " available for equation %s\n"), eq->name);

  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  /* Perform the computation */
  eq->compute_cellwise_diff_flux(fld->val,
                                 eq->builder,
                                 diff_flux);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
