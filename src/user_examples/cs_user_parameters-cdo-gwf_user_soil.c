/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

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
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \subpage parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Local structure definition
 *============================================================================*/

/*! [param_cdo_gwf_tracy_struct] */
/* Parameters defining the unsaturated soil model devised by Tracy */

typedef struct {

  double    L;                  /* column length */
  double    h_s;                /* saturated head reference */
  double    h_r;                /* residual head reference */
  double    theta_r;            /* residual moisture */
  double    theta_s;            /* saturated moisture */

} cs_tracy_param_t;
/*! [param_cdo_gwf_tracy_struct] */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  The prototype of this function should fit the one of the generic
 *         function pointer \ref cs_gwf_soil_update_t
 *         The model follows in this function has been defined by Tracy
 *         for verification purposes in the case of unsaturated single phase
 *         flows.
 *         The purpose of this function is to update soil properties such as
 *         the moisture content (liquid saturation), (full) permeability and
 *         the soil capacity.
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 * \param[in]      quant       pointer to a cs_cdo_quantities_t structure
 * \param[in]      zone        pointer to a cs_zone_t
 * \param[in, out] soil        pointer to the soil structure to update
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_gwf_set_user_update_soil] */
static void
tracy_update(const cs_real_t              t_eval,
             const cs_mesh_t             *mesh,
             const cs_cdo_connect_t      *connect,
             const cs_cdo_quantities_t   *quant,
             const cs_zone_t             *zone,
             cs_gwf_soil_t               *soil)
{
  CS_UNUSED(mesh);
  CS_UNUSED(connect);
  CS_UNUSED(quant);
  CS_UNUSED(t_eval);

  /* Retrieve the soil parameters */

  const cs_tracy_param_t  *sp = (cs_tracy_param_t *)soil->model_param;

  /* Retrieve the hydraulic context */

  cs_gwf_unsaturated_single_phase_t  *hc = soil->hydraulic_context;

  /* Additional parameters */

  const cs_real_t  *head_values = hc->head_in_law;
  const double  delta_m = soil->porosity - sp->theta_r;
  const double  k_s = soil->abs_permeability[0][0];

  /* Retrieve field values associated to properties to update */

  cs_real_t  *permeability = hc->permeability_field->val;
  cs_real_t  *moisture = hc->moisture_field->val;
  cs_real_t  *capacity = hc->capacity_field->val;

  for (cs_lnum_t i = 0; i < zone->n_elts; i++) {

    const cs_lnum_t  c_id = zone->elt_ids[i];
    const cs_real_t  h = head_values[c_id];
    const cs_real_t  k_r = (h - sp->h_r)/(sp->h_s - sp->h_r);

    /* Set the permeability value */

    permeability[c_id] = k_s * k_r;

    /* Set the moisture content (Se = 1 in this case)*/

    moisture[c_id] = sp->theta_r + k_r * delta_m;

    /* Set the capacity values */

    capacity[c_id] = delta_m /(sp->h_s - sp->h_r);

  } /* Loop on selected cells */

}
/*! [param_cdo_gwf_set_user_update_soil] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  The prototype of this function should fit the one of the generic
 *         function pointer \ref cs_gwf_soil_free_param_t
 *         The purpose of this function is to free the context structure
 *         associated to the user-defined model.
 *
 * \param[in, out] p_soil_param   double pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_gwf_set_user_free_soil] */
static void
tracy_free_param(void         **p_soil_param)
{
  cs_tracy_param_t  *sp = (cs_tracy_param_t *)(*p_soil_param);

  BFT_FREE(sp);
  *p_soil_param = NULL;
}
/*! [param_cdo_gwf_set_user_free_soil] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the boundary condition values for the Richards equation.
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time         when ?
 * \param[in]      n_pts        number of elements to consider
 * \param[in]      pt_ids       list of elements ids (e.g. to access xyz)
 * \param[in]      xyz          where ?
 * \param[in]      dense_output true:no indirection, false:apply pt_ids
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval       result of the function
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_gwf_set_bc_analytic] */
static void
get_bc(cs_real_t           time,
       cs_lnum_t           n_pts,
       const cs_lnum_t    *pt_ids,
       const cs_real_t    *xyz,
       bool                dense_output,
       void               *input,
       cs_real_t          *retval)
{
  const cs_gwf_soil_t  *soil = input;
  assert(soil != NULL);
  const cs_tracy_param_t  *tp = soil->model_param;
  assert(tp != NULL);

  /* Physical parameters */

  const double  overL = 1./tp->L;
  const double  dtheta = soil->porosity - tp->theta_r;
  const double  td =
    -5 * tp->L * tp->L * dtheta /( 6 * tp->h_r * soil->abs_permeability[0][0] );
  const double  alpha = 6 - 5*time/td;

  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const cs_lnum_t  id = (pt_ids == NULL) ? p : pt_ids[p];
    const cs_lnum_t  r_id = dense_output ? p : id;
    const double  xll = (xyz[3*id] - tp->L)*overL, beta = xll*xll;

    retval[r_id] = tp->h_r*(1 - beta/alpha);

  } /* Loop on selected points */
}
/*! [param_cdo_gwf_set_bc_analytic] */

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the explicit definition of the initial solution for the Richards
 *         equation. Same as get_sol but optimize for time=0.
 *         pt_ids is optional. If not NULL, it enables to access in coords
 *         at the right location and the same thing to fill retval if compact
 *         is set to false
 *         Rely on a generic function pointer for an analytic function
 *
 * \param[in]      time         when ?
 * \param[in]      n_pts        number of elements to consider
 * \param[in]      pt_ids       list of elements ids (e.g. to access xyz)
 * \param[in]      xyz          where ?
 * \param[in]      dense_output true:no indirection, false:apply pt_ids
 * \param[in]      input        NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval       result of the function
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_set_ic_by_analytic] */
static void
get_ic(cs_real_t           time,
       cs_lnum_t           n_pts,
       const cs_lnum_t    *pt_ids,
       const cs_real_t    *xyz,
       bool                dense_output,
       void               *input,
       cs_real_t          *retval)
{
  CS_UNUSED(time);

  cs_tracy_param_t  *tp = input;
  assert(input != NULL);

  const double  one6 = 1./6;

  for (cs_lnum_t p = 0; p < n_pts; p++) {

    const cs_lnum_t  id = (pt_ids == NULL) ? p : pt_ids[p];
    const cs_lnum_t  r_id = (dense_output) ? p : id;
    const double  x = xyz[3*id], xll = (x - tp->L)/tp->L;

    retval[r_id] = tp->h_r*(1-one6*xll*xll);

  } /* Loop on selected points */
}
/*! [param_cdo_set_ic_by_analytic] */

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
  cs_domain_t  *domain = cs_glob_domain;

  /* ===================
     Activate CDO module
     =================== */

  cs_domain_set_cdo_mode(domain, CS_DOMAIN_CDO_MODE_ONLY);

  /* ======================================
     Boundaries of the computational domain
     ====================================== */

  /* Choice of the default type of boundary.
     A valid choice is either CS_BOUNDARY_WALL or CS_BOUNDARY_SYMMETRY */

  cs_boundary_set_default(domain->boundaries, CS_BOUNDARY_SYMMETRY);

  cs_boundary_add(domain->boundaries, CS_BOUNDARY_INLET, "left");
  cs_boundary_add(domain->boundaries, CS_BOUNDARY_OUTLET, "right");

  /* 1. Activate the groundwater flow module */

  /*! [param_cdo_activate_gwf_uspf] */

  cs_gwf_activate(CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE,
                  0,                          /* No option to set */
                  CS_GWF_POST_PERMEABILITY);  /* Post-processing options */

  /*! [param_cdo_activate_gwf_uspf] */

  /* 2. Add and define soils */

  /*! [param_cdo_gwf_add_user_soil] */

  cs_gwf_soil_t  *s =
    cs_gwf_add_iso_soil("cells",
                        1.0,        /* bulk mass density (useless here) */
                        1.15741e-4, /* absolute permeability */
                        0.45,       /* porosity */
                        CS_GWF_SOIL_USER);

  /* 2.a Create and define a structure of parameters to manage the soil */

  cs_tracy_param_t  *tp = NULL;

  BFT_MALLOC(tp, 1, cs_tracy_param_t);

  tp->L = 200;
  tp->h_s = 0.;
  tp->h_r = -100;
  tp->theta_r = 0.15;

  /* 2.b Associate the parameter structure and the user-defined functions
   *     to manage the soil */

  cs_gwf_soil_set_user(s,                  /* soil structure */
                       tp,                 /* soil parameter structure */
                       tracy_update,       /* function to update the soil */
                       tracy_free_param);  /* function to free the structure */

  /*! [param_cdo_gwf_add_user_soil] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t     *domain)
{
  /* =========================
     Generic output management
     ========================= */

  cs_domain_set_output_param(domain,
                             CS_RESTART_INTERVAL_ONLY_AT_END,
                             1,  /* output log frequency */
                             2); /* verbosity (-1: no, 0, ...) */

  /* ====================
     Time step management
     ====================

     If there is an inconsistency between the max. number of time iterations
     and the final physical time, the first condition encountered stops
     the calculation.
  */

  cs_domain_set_time_param(domain,
                           200,       /* nt_max or -1 (automatic) */
                           864000.);  /* t_max or < 0. (automatic) */

  /* Define the value of the time step */

  cs_domain_def_time_step_by_value(domain, 4320);

  /* ====================
     Numerical parameters
     ==================== */

  /* Modify the setting of an equation using a generic process
   * See Doxygen HTML documentation for more details
   */

  /* Richards equation
     ================= */

  cs_equation_param_t  *eqp = cs_equation_param_by_name("Richards");

  /* Linear algebra settings. Switch solver w.r.t. to which librairy is
     available */

#if defined(HAVE_PETSC)
  cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "petsc");
  cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "cg");
  cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "amg");
#if defined(PETSC_HAVE_HYPRE)
  cs_equation_param_set(eqp, CS_EQKEY_AMG_TYPE, "boomer");
#else
  cs_equation_param_set(eqp, CS_EQKEY_AMG_TYPE, "gamg");
#endif
#else
  cs_equation_param_set(eqp, CS_EQKEY_SOLVER_FAMILY, "cs");
  cs_equation_param_set(eqp, CS_EQKEY_PRECOND, "jacobi");
  cs_equation_param_set(eqp, CS_EQKEY_ITSOL, "cg");
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify output user parameters.
 *
 * For CDO schemes, this function concludes the setup of properties,
 * equations, source terms...
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /*! [param_cdo_gwf_get_user_soil] */
  /* Final settings for the Richards equation */

  cs_equation_param_t  *eqp = cs_equation_param_by_name("Richards");


  /* 1. Retrieve the soil by its name and then its parameter structure */

  cs_gwf_soil_t  *soil = cs_gwf_soil_by_name("cells");
  cs_tracy_param_t  *tp = soil->model_param;

  /* Define the boundary conditions  */

  cs_equation_add_bc_by_analytic(eqp,
                                 CS_PARAM_BC_DIRICHLET,
                                 "left", // boundary zone name
                                 get_bc,
                                 soil);

  cs_equation_add_bc_by_value(eqp,
                              CS_PARAM_BC_DIRICHLET,
                              "right", // boundary zone name
                              &(tp->h_r));

  /* Set the initial condition */

  cs_equation_add_ic_by_analytic(eqp,
                                 NULL,
                                 get_ic,
                                 tp);

  /*! [param_cdo_gwf_get_user_soil] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
