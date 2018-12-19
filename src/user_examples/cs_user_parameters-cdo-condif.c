/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

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
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_advection_field.h"
#include "cs_base.h"
#include "cs_domain_setup.h"
#include "cs_equation.h"
#include "cs_equation_param.h"
#include "cs_field.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_multigrid.h"
#include "cs_halo.h"
#include "cs_param.h"
#include "cs_property.h"
#include "cs_sles.h"
#include "cs_time_step.h"
#include "cs_walldistance.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-base.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \subpage parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Give the explicit definition of the advection field.
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time       when ?
 * \param[in]      n_elts     number of elements to consider
 * \param[in]      pt_ids     list of elements ids (to access coords and fill)
 * \param[in]      coords     where ?
 * \param[in]      compact    true:no indirection, false:indirection for filling
 * \param[in]      input      NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval     result of the function
 */
/*----------------------------------------------------------------------------*/

static void
_define_adv_field(cs_real_t           time,
                  cs_lnum_t           n_pts,
                  const cs_lnum_t    *pt_ids,
                  const cs_real_t    *xyz,
                  bool                compact,
                  void               *input,
                  cs_real_t          *res)
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const cs_lnum_t  id = (pt_ids == NULL) ? p : pt_ids[p];
    const cs_lnum_t  ii = compact ? p : id;
    const cs_real_t  *pxyz = xyz + 3*id;
    cs_real_t  *pres = res + 3*ii;

    pres[0] = pxyz[1] - 0.5;
    pres[1] = 0.5 - pxyz[0];
    pres[2] = pxyz[2] + 1;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Give the explicit definition of the dirichlet boundary conditions
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time      when ?
 * \param[in]      n_elts    number of elements to consider
 * \param[in]      pt_ids    list of elements ids (to access coords and fill)
 * \param[in]      coords    where ?
 * \param[in]      compact   true:no indirection, false:indirection for filling
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 * \param[in, out] res       result of the function
 */
/*----------------------------------------------------------------------------*/

static void
_define_bcs(cs_real_t           time,
            cs_lnum_t           n_pts,
            const cs_lnum_t    *pt_ids,
            const cs_real_t    *xyz,
            bool                compact,
            void               *input,
            cs_real_t          *res)
{
  CS_UNUSED(time);
  CS_UNUSED(input);

  const double  pi = 4.0*atan(1.0);
  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const cs_lnum_t  id = (pt_ids == NULL) ? p : pt_ids[p];
    const cs_lnum_t  ii = compact ? p : id;
    const cs_real_t  *_xyz = xyz + 3*id;
    const double  x = _xyz[0], y = _xyz[1], z = _xyz[2];

    res[ii] = 1 + sin(pi*x)*sin(pi*(y+0.5))*sin(pi*(z+0.25));

  }

}

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{
  /* Activate CDO/HHO module so that main additional structure are built */

  /*! [param_cdo_activation] */
  {
    cs_domain_t  *domain = cs_glob_domain;

    cs_domain_set_cdo_mode(domain, CS_DOMAIN_CDO_MODE_ONLY);
  }
  /*! [param_cdo_activation] */

  /* ======================
     Boundary of the domain
     ====================== */

  /*! [param_cdo_domain_boundary] */
  {
    cs_domain_t  *domain = cs_glob_domain;
    cs_boundary_t  *bdy = domain->boundaries;

    /* Choose a boundary by default */
    cs_boundary_set_default(bdy, CS_BOUNDARY_WALL);

    /* Add a new boundary
     *  >> cs_domain_add_boundary(boundary, type_of_boundary, zone_name);
     *
     * zone_name is either a predefined one or user-defined one
     * type_of_boundary is one of the following keywords:
     *   CS_BOUNDARY_WALL,
     *   CS_BOUNDARY_INLET,
     *   CS_BOUNDARY_OUTLET,
     *   CS_BOUNDARY_SYMMETRY
     */

    cs_boundary_add(bdy, CS_BOUNDARY_INLET, "in");
    cs_boundary_add(bdy, CS_BOUNDARY_OUTLET, "out");
  }
  /*! [param_cdo_domain_boundary] */

  /* =========================
     Generic output management
     ========================= */

  /*! [param_cdo_domain_output] */
  {
    cs_domain_t  *domain = cs_glob_domain;

    cs_domain_set_output_param(domain,
                               -1,     // restart frequency
                               10,     // output log frequency
                               2);     // verbosity (-1: no, 0, ...)

  }
  /*! [param_cdo_domain_output] */

    /* ====================
       Time step management
       ==================== */

  /*! [param_cdo_time_step] */
  {
    cs_domain_t  *domain = cs_glob_domain;

    /*
       If there is an inconsistency between the max. number of iteration in
       time and the final physical time, the first condition encountered stops
       the calculation.
    */

    cs_domain_set_time_param(domain,
                             100,     // nt_max or -1 (automatic)
                             -1.);    // t_max or < 0. (automatic)

    /* Define the value of the time step
       >> cs_domain_def_time_step_by_value(domain, dt_val);
       >> cs_domain_def_time_step_by_func(domain, dt_func);

       The second way to define the time step enable complex definitions.
    */

    cs_domain_def_time_step_by_value(domain, 1.0);

  }
  /*! [param_cdo_time_step] */

  /*! [param_cdo_wall_distance] */
  {
    cs_walldistance_activate();
  }
  /*! [param_cdo_wall_distance] */

  /*! [param_cdo_add_user_equation] */
  {
    /* Add a new user equation:
       Set the default boundary condition among:
       CS_PARAM_BC_HMG_DIRICHLET or
       CS_PARAM_BC_HMG_NEUMANN

       By default, initial values are set to zero (or the value given by the
       restart file in case of restart).
    */

    cs_equation_add_user("AdvDiff.Upw", // equation name
                         "Pot.Upw",     // associated variable field name
                         1,             // dimension of the unknown
                         CS_PARAM_BC_HMG_DIRICHLET); // default boundary

    cs_equation_add_user("AdvDiff.SG",  // equation name
                         "Pot.SG",      // associated variable field name
                         1,             // dimension of the unknown
                         CS_PARAM_BC_HMG_DIRICHLET); // default boundary
  }
  /*! [param_cdo_add_user_equation] */

  /* ========================================
     Add material properties
     ======================================== */

  /*! [param_cdo_add_user_properties] */
  {
    cs_property_add("conductivity",      /* property name */
                    CS_PROPERTY_ANISO);  /* type of material property */
    cs_property_add("rho.cp",            /* property name */
                    CS_PROPERTY_ISO);    /* type of material property */

  }
  /*! [param_cdo_add_user_properties] */

  /*! [param_cdo_add_user_properties_opt] */
  {
    /* Retrieve a property named "conductivity"  */
    cs_property_t  *pty = cs_property_by_name("conductivity");

    /* Activate the computation of the Fourier number for this property */
    cs_property_set_option(pty, CS_PTYKEY_POST_FOURIER);
  }
  /*! [param_cdo_add_user_properties_opt] */

  /* ========================================
     Add advection fields
     ======================================== */

  /*! [param_cdo_add_user_adv_field] */
  {
    /* Add a user-defined advection field named "adv_field"  */
    cs_adv_field_t  *adv = cs_advection_field_add_user("adv_field");
  }
  /*! [param_cdo_add_user_adv_field] */

  /*! [param_cdo_add_user_adv_field_opt] */
  {
    /* Retrieve an advection field named "adv_field"  */
    cs_adv_field_t  *adv = cs_advection_field_by_name("adv_field");

    /* Compute the Courant number (if unsteady simulation) */
    cs_advection_field_set_option(adv, CS_ADVKEY_POST_COURANT);

    /* Set other advanced options: for instance, define an interpolation of the
       advection field at vertices */
    cs_advection_field_set_option(adv, CS_ADVKEY_DEFINE_AT_VERTICES);

    /* Both options in one call */
    cs_advection_field_set_option(adv,
                                  CS_ADVKEY_POST_COURANT |
                                  CS_ADVKEY_DEFINE_AT_VERTICES);
  }
  /*! [param_cdo_add_user_adv_field_opt] */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(void)
{
  /* ─────────────────────────────────────────────────────────────────────────
   *  Modify the setting of an equation using a generic process
   *  cf. the DOXYGEN documentation (website or in your installation path)
   *
   *         cs_equation_set_param(eqp, CS_EQKEY_*, "key_val")
   *
   * ─────────────────────────────────────────────────────────────────────────*/

  /*! [param_cdo_numerics] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff.Upw");

    /* The modification of the space discretization should be apply first */
    cs_equation_set_param(eqp, CS_EQKEY_SPACE_SCHEME, "cdo_vb");
    cs_equation_set_param(eqp, CS_EQKEY_ADV_SCHEME, "upwind");

    /* Modify other parameters than the space discretization */
    cs_equation_set_param(eqp, CS_EQKEY_VERBOSITY, "2");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_ALGO, "cost");
    cs_equation_set_param(eqp, CS_EQKEY_HODGE_DIFF_COEF, "dga");

    /* Linear algebra settings */
#if defined(HAVE_PETSC)
    cs_equation_set_param(eqp, CS_EQKEY_SOLVER_FAMILY, "petsc");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "amg");
#else
    cs_equation_set_param(eqp, CS_EQKEY_SOLVER_FAMILY, "cs");
    cs_equation_set_param(eqp, CS_EQKEY_PRECOND, "jacobi");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL, "cg");
#endif
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL_MAX_ITER, "2500");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL_EPS, "1e-12");
    cs_equation_set_param(eqp, CS_EQKEY_ITSOL_RESNORM, "false");

  }
  /*! [param_cdo_numerics] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Advanced user-defined settings for the linear algebra related
 *         to CDO equations
 *         This is closed to cs_user_linear_solvers() but called once the fields
 *         and equations have been created (this happens at a different stage
 *         in the CDO framework)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_linear_solvers(void)
{
/*! [param_cdo_mg_aggreg] */
  {
    cs_equation_t  *eq = cs_equation_by_name("AdvDiff.Upw");
    cs_equation_param_t  *eqp = cs_equation_get_param(eq);
    cs_field_t  *fld = cs_equation_get_field(eq);

    /* In case of a in-house K-cylcle multigrid as a preconditioner of a
       linear iterative solver */
    if (eqp->sles_param.precond == CS_PARAM_PRECOND_AMG) {
      /* If multigrid is the chosen preconditioner */
      if (eqp->sles_param.amg_type == CS_PARAM_AMG_HOUSE_K) {
        /* If this is a K-cycle multigrid */

        /* Retrieve the different context structures to apply additional
           settings */
        cs_sles_t  *sles = cs_sles_find_or_add(fld->id, NULL);
        cs_sles_it_t  *itsol = cs_sles_get_context(sles);
        cs_sles_pc_t  *pc = cs_sles_it_get_pc(itsol);
        cs_multigrid_t  *mg = cs_sles_pc_get_context(pc);

        /* Available settings:
         * - max. number of elements in an aggregation
         * - type of algorithm to perform the aggregation
         * - max. number of levels (i.e. grids)
         * - max globalnumber of rows at the coarsest level
         * - type of relaxation (weighting between a P_0 and P_1). For K-cycle,
         * this should be equal to 0.
         * - Activation of the postprocessing for the aggregation if > 0.
         * Aggregation set is numbered by its coarse row number modulo this
         * value
         */
        cs_multigrid_set_coarsening_options(mg,
                                            8,   /* aggregation_limit*/
                                            CS_GRID_COARSENING_SPD_PW,
                                            10,  /* n_max_levels */
                                            30,  /* min_g_cells (default 30) */
                                            0.,  /* P0P1 relaxation */
                                            12); /* postprocess (default 0) */

      } /* K-cycle */
    }   /* Multigrid as preconditioner */
  }
  /*! [param_cdo_mg_aggreg] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  - Specify the elements such as properties, advection fields,
 *           user-defined equations and modules which have been previously
 *           added.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_cdo_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* =======================
     User-defined properties
     ======================= */

  /*! [param_cdo_setup_property] */
  {
    cs_property_t  *conductivity = cs_property_by_name("conductivity");
    cs_real_33_t  tensor = {{1.0,  0.5, 0.0}, {0.5, 1.0, 0.5}, {0.0, 0.5, 1.0}};

    cs_property_def_aniso_by_value(conductivity, // property structure
                                   "cells",      // name of the volume zone
                                   tensor);      // values of the property

    cs_property_t  *rhocp = cs_property_by_name("rho.cp");
    cs_real_t  iso_val = 2.0;

    cs_property_def_iso_by_value(rhocp,    // property structure
                                 "cells",  // name of the volume zone
                                 iso_val); // value of the property

  }
  /*! [param_cdo_setup_property] */

  /* =============================
     User-defined advection fields
     ============================= */

  /*! [param_cdo_setup_advfield] */
  {
    cs_adv_field_t  *adv = cs_advection_field_by_name("adv_field");

    cs_advection_field_def_by_analytic(adv, _define_adv_field, NULL);

    /* Enable also the defintion of the advection field at mesh vertices */
    cs_advection_field_set_option(adv, CS_ADVKEY_DEFINE_AT_VERTICES);

    /* Activate the post-processing of the related Courant number */
    cs_advection_field_set_option(adv, CS_ADVKEY_POST_COURANT);
  }
  /*! [param_cdo_setup_advfield] */

  /* ======================
     User-defined equations
     ====================== */

  /* Define the boundary conditions
     >> cs_equation_add_bc_by_analytic(eqp,
                                       bc_type,
                                       "zone_name",
                                       analytic_function);

     -> eq is the structure related to the equation to set
     -> type of boundary condition:
        CS_PARAM_BC_DIRICHLET, CS_PARAM_BC_HMG_DIRICHLET,
        CS_PARAM_BC_NEUMANN, CS_PARAM_BC_HMG_NEUMANN, CS_PARAM_BC_ROBIN

     >> cs_equation_add_bc_by_value(eqp,
                                    bc_type,
                                    "mesh_location_name",
                                    values); // pointer
  */

  /*! [param_cdo_setup_bcs] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff.Upw");

    cs_equation_add_bc_by_analytic(eqp,
                                   CS_PARAM_BC_DIRICHLET,
                                   "boundary_faces",  // zone name
                                   _define_bcs,       // pointer to the function
                                   NULL);             // input structure

  }
  /*! [param_cdo_setup_bcs] */

  /*! [param_cdo_add_terms] */
  {
    cs_equation_param_t  *eqp = cs_equation_param_by_name("AdvDiff.Upw");
    cs_property_t  *rhocp = cs_property_by_name("rho.cp");
    cs_property_t  *conductivity = cs_property_by_name("conductivity");
    cs_adv_field_t  *adv = cs_advection_field_by_name("adv_field");

    /* Activate the unsteady term */
    cs_equation_add_time(eqp, rhocp);

    /* Activate the diffusion term */
    cs_equation_add_diffusion(eqp, conductivity);

    /* Activate advection effect */
    cs_equation_add_advection(eqp, adv);

    /* Simple definition with cs_equation_add_source_term_by_val
       where the value of the source term is given by m^3
    */
    cs_real_t  st_val = -0.1;
    cs_xdef_t  *st = cs_equation_add_source_term_by_val(eqp,
                                                        "cells",
                                                        &st_val);

  }
  /*! [param_cdo_add_terms] */

  /*! [param_cdo_copy_settings] */
  {
    /* Copy the settings for AdvDiff.Upw */
    cs_equation_param_t  *eqp_ref = cs_equation_param_by_name("AdvDiff.Upw");
    cs_equation_t *eq = cs_equation_by_name("AdvDiff.SG");
    cs_equation_param_t  *eqp = cs_equation_get_param(eq);

    /* Copy the settings */
    cs_equation_param_update_from(eqp_ref, eqp);

    /* Keep all the settings from "AdvDiff.Upw and then only change the
       advection scheme for the second equation */
    cs_equation_set_param(eqp, CS_EQKEY_ADV_SCHEME, "sg");

    /* Call this function to be sure that the linear solver is set to what
       one wants */
    cs_equation_param_set_sles(eqp, cs_equation_get_field_id(eq));

  }
  /*! [param_cdo_copy_settings] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
