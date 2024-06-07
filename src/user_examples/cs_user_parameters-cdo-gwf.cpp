/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-cdo-gwf.c
 *
 * \brief User functions for setting a calculation using the groundwater flow
 *        module with CDO schemes
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/* Permeability in each subdomain */

static const double k1 = 1e5, k2 = 1;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the initial condition for the tracer attached to the GWF problem
 *         Generic function pointer for an evaluation relying on an analytic
 *         function.
 *         pt_ids is optional. If not NULL, it enables to access to the xyz
 *         array with an indirection. The same indirection can be applied to
 *         fill retval if dense_output is set to false.
 *
 * \param[in]      time          when ?
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      pt_ids        list of elements ids (in coords and retval)
 * \param[in]      xyz           where ? Coordinates array
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         NULL or pointer to a structure cast on-the-fly
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

/*! [param_cdo_gwf_get_tracer_ic] */
static inline void
get_tracer_ic(cs_real_t          time,
              cs_lnum_t          n_elts,
              const cs_lnum_t   *pt_ids,
              const cs_real_t   *xyz,
              bool               dense_output,
              void              *input,
              cs_real_t         *retval)
{
  CS_UNUSED(input);

  /* Physical parameters */

  const double  magnitude = 2*k1/(k1 + k2);
  const double  x_front = magnitude * time;

  /* Fill retval array */

  for (cs_lnum_t  i = 0; i < n_elts; i++) {

    const cs_lnum_t  id = (pt_ids == NULL) ? i : pt_ids[i];
    const cs_lnum_t  r_id = (dense_output) ? i : id;
    const double  x = xyz[3*id];

    if (x <= x_front)
      retval[r_id] = 1;
    else
      retval[r_id] = 0;

  } /* Loop on selected elements */

}
/*! [param_cdo_gwf_get_tracer_ic] */

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

  /* ======================
     Boundary of the domain
     ====================== */

  cs_boundary_t  *bdy = domain->boundaries;

  /* Choose a boundary by default */

  cs_boundary_set_default(bdy, CS_BOUNDARY_SYMMETRY);

  /* Add new boundaries */

  cs_boundary_add(bdy, CS_BOUNDARY_INLET, "left");
  cs_boundary_add(bdy, CS_BOUNDARY_OUTLET, "right");

  /* Activate CDO/HHO module so that main additional structure are built */

  cs_param_cdo_mode_set(CS_PARAM_CDO_MODE_ONLY);

  /* =================================
     Define the groundwater flow model
     ================================= */

  /* 1. Activate groundwater flow module */
  /* ----------------------------------- */

  /*! [param_cdo_activate_gwf] */
  {
    /* For the groundwater flow module:
       cs_gwf_activate(model_type, option_flag, post_flag);

       If option_flag or post_flag = 0, then there is nothing to set. */

    cs_flag_t  option_flag = 0, post_flag = 0;

    cs_gwf_activate(CS_GWF_MODEL_SATURATED_SINGLE_PHASE,
                    option_flag, post_flag);
  }
  /*! [param_cdo_activate_gwf] */

  /*! [param_cdo_activate_gwf_b] */
  {
    cs_gwf_activate(CS_GWF_MODEL_UNSATURATED_SINGLE_PHASE,
                    /* Physical modelling option or numerical option */
                    CS_GWF_GRAVITATION,/* Take into account the gravity */
                    /* Automatic postprocessing options */
                    CS_GWF_POST_PERMEABILITY |
                    CS_GWF_POST_DARCY_FLUX_BALANCE);

    /* In this case, the gravity vector has to be defined (either using the GUI
       or in cs_user_parameters() function */
  }
  /*! [param_cdo_activate_gwf_b] */

  /*! [param_cdo_post_gwf] */
  /* Specify post-processing options */

  cs_gwf_set_post_options(CS_GWF_POST_PERMEABILITY |
                          CS_GWF_POST_DARCY_FLUX_BALANCE |
                          CS_GWF_POST_DARCY_FLUX_DIVERGENCE |
                          CS_GWF_POST_DARCY_FLUX_AT_BOUNDARY,
                          false); /* No reset */

  /*! [param_cdo_post_gwf] */

  /* 2. Add and define soils (must be done before adding tracers) */
  /* ----------------------- */

  /*! [param_cdo_gwf_add_define_iso_saturated_soil] */
  {
    /* Example 1: Simplest case. A "saturated" isotropic soils */

    /* saturated isotropic permeability */

    const double  iso_val = 1e-10;
    const double  theta_s = 0.85;
    const cs_real_t  bulk_density = 2500; /* useless if no tracer is
                                             considered */

    cs_gwf_add_iso_soil("cells", bulk_density, iso_val, theta_s,
                        CS_GWF_SOIL_SATURATED);
  }
  /*! [param_cdo_gwf_add_define_iso_saturated_soil] */

  /*! [param_cdo_gwf_add_define_aniso_saturated_soil] */
  {
    /* Example 2: Two "saturated" and anisotropic soils */

    /* saturated anisotropic permeability */

    double  aniso_val1[3][3] = {{1e-3, 0,    0},
                                {   0, 1,    0},
                                {   0, 0, 1e-1}};
    double  aniso_val2[3][3] = {{1e-5, 0,    0},
                                {   0, 1,    0},
                                {   0, 0, 1e-2}};
    const double  theta_s = 1;
    const double  bulk_density = 1.0; /* useless if no tracer is considered */

    /* One assumes that two (volume) zones have be defined called "soil1" and
       "soil2". */

    cs_gwf_add_aniso_soil("soil1", bulk_density, aniso_val1, theta_s,
                          CS_GWF_SOIL_SATURATED);

    cs_gwf_add_aniso_soil("soil2", bulk_density, aniso_val2, theta_s,
                          CS_GWF_SOIL_SATURATED);
  }
  /*! [param_cdo_gwf_add_define_aniso_saturated_soil] */

  /*! [param_cdo_gwf_add_define_genuchten_soil] */
  {
    /* Example 3: Add a new soil based on the Van Genuchten-Mualem model for
       all the mesh cells */

    cs_gwf_soil_model_t  model = CS_GWF_SOIL_VGM_SINGLE_PHASE;
    cs_gwf_soil_t  *s = cs_gwf_add_iso_soil("cells", /* volume zone name */
                                            1800,    /* bulk mass density */
                                            3e-1,    /* absolute permeability */
                                            0.9,     /* porosity */
                                            model);

    /* Set now the parameters associated to this model of soil */

    cs_gwf_soil_set_vgm_spf_param(s,
                                  0.078,      /* residual moisture */
                                  0.036,      /* scaling parameter */
                                  1.56,       /* (n) shape parameter */
                                  0.5);       /* (L) tortuosity */
  }
  /*! [param_cdo_gwf_add_define_genuchten_soil] */

  /*! [param_cdo_gwf_add_user_soil] */
  {
    /* Example 4: Add a new user-defined soil for all the cells */

    cs_gwf_add_iso_soil("cells", /* zone name associated to the soil */
                        1800,    /* bulk mass density (useless without tracer */
                        3e-1,    /* absolute permeability */
                        0.9,     /* porosity */
                        CS_GWF_SOIL_USER);
  }
  /*! [param_cdo_gwf_add_user_soil] */

  /*! [param_cdo_gwf_get_soil] */
  {
    /* If a soil structure has been added but not already defined, one can
     * retrieve it thanks to \ref cs_gwf_soil_by_name
     *
     * The name of the soil is the same as the name of the volume zone used at
     * the creation of the soil
     */

    cs_gwf_soil_t  *s1 = cs_gwf_soil_by_name("soil1");

    /* s1 pointer can be used to set futher options in advanced settings */

    assert(s1 != NULL);
  }
  /*! [param_cdo_gwf_get_soil] */

  /* 3. Add and define tracer equations */
  /* ---------------------------------- */

  /*! [param_cdo_gwf_add_tracer] */
  {
    /*
      Add a tracer equation which is unsteady and convected by the darcean
      flux. The created equation is called "eqname" along with a new field
      called "varname".
    */

    cs_gwf_tracer_model_t  model = 0; /* Default model without precipitation
                                         effect */

    cs_gwf_tracer_t  *tr = cs_gwf_add_tracer(/* Tracer model */
                                             model,
                                             /* Equation name */
                                             "Tracer_01",
                                             /* Variable name */
                                             "C");

    /* For a "default" tracer, the definition can be made as follows.
     *
     * The parameters defining the tracer behavior can be set soil by soil
     * (thanks to the soil name as the second argument) or to all soils in one
     * call (using NULL as the second argument)
     */

    cs_gwf_tracer_set_soil_param(tr,
                                 NULL,     /* soil name or NULL for all */
                                 0.,       /* water molecular diffusivity */
                                 1., 0.,   /* alpha (longi. and transvesal) */
                                 1e-4);    /* distribution coef. */
  }
  /*! [param_cdo_gwf_add_tracer] */

  /*! [param_cdo_gwf_add_rtracer] */
  {
    /*
      Add a radioactive tracer equation which is unsteady and convected by the
      darcean flux. A reaction term related to the first order decay
      coefficient is automatically added to the tracer equation. The created
      equation is called "eqname" along with a new field called "varname".
    */

    /* Default tracer model with precipitation effects */

    cs_gwf_tracer_model_t  model = CS_GWF_TRACER_PRECIPITATION;

    cs_gwf_tracer_t  *rtr = cs_gwf_add_radioactive_tracer(/* Tracer model */
                                                          model,
                                                          /* Equation name */
                                                          "Tracer",
                                                          /* Variable name */
                                                          "C",
                                                          /* Decay coeff. */
                                                          0.01);

    /* For a "default" tracer, the definition can be made as follows.
     *
     * The parameters defining the tracer behavior can be set soil by soil
     * (thanks to the soil name as the second argument) or to all soils in one
     * call (using NULL as the second argument)
     */

    cs_gwf_tracer_set_soil_param(rtr,
                                 NULL,     /* soil name or NULL for all */
                                 0.,       /* water molecular diffusivity */
                                 1., 0.,   /* alpha (longi. and transvesal) */
                                 1e-4);    /* distribution coef. */

    /* Set the parameters for the precipitation model */

    cs_gwf_tracer_set_precip_param(rtr,
                                   NULL,    /* soil name/NULL for all */
                                   1e-4);   /* liquid concentration threshold */
  }
  /*! [param_cdo_gwf_add_rtracer] */

  /*! [param_cdo_gwf_add_decay_chain] */
  {
    /*
      Add a radioactive decay chain defining a set of linked tracer equations
    */

    /* Example of a chain with 3 tracers with the last one taking into account
       the precipitation effects. Each tracer is associated to a quantity in
       mole */

    int  n_tracers = 3;
    const char  *species_names[3] = {"grandfather", "father", "son"};
    cs_gwf_tracer_model_t  models[3] = {0, 0, CS_GWF_TRACER_PRECIPITATION};
    double  decay_rates[3] = {1e-3, 1e-5, 1e-4};

    cs_gwf_tracer_decay_chain_t  *tdc =
      cs_gwf_add_decay_chain(n_tracers,
                             CS_GWF_TRACER_UNIT_MOLE,
                             "my_decay_chain",
                             species_names,
                             models,
                             decay_rates);

    /* Now set the parameter associated to each couple (soil, tracer).
     *
     * Example with two soils named: "soil1" and "soil2"
     */

    const char  *soil_names[2] = {"soil1", "soil2"};
    const double  kd_values[2][3] = {{1e-3, 5e-3, 4e-5},
                                     {1e-3, 1e-3, 1e-5}};

    for (int i = 0; i < tdc->n_tracers; i++) {

      /* Retrieve the tracer in the chain */

      cs_gwf_tracer_t  *rtr = tdc->tracers[i];

      /* Set now the tracer for each soil */

      for (int is = 0; is < 2; is++)
        cs_gwf_tracer_set_soil_param(rtr,
                                     soil_names[is],
                                     0.,       /* water molecular diffusivity */
                                     1., 0.,   /* alpha (longi. and trans.) */
                                     kd_values[is][i]); /* distribution coef. */

      /* Set the parameters for the precipitation model. One assumes the
         parameter is the same for the two soils */

      if (i == tdc->n_tracers - 1)
        cs_gwf_tracer_set_precip_param(rtr,
                                       NULL,  /* soil name/NULL for all */
                                       1e-4); /* liquid conc. threshold */

    } /* Loop on tracers */

  }
  /*! [param_cdo_gwf_add_decay_chain] */

  /*! [param_cdo_gwf_get_tracer] */
  {
    /* If the tracer structure has been added but not totally defined, one can
     * retrieve it thanks to cs_gwf_tracer_by_name
     *
     * The name of the tracer is the same as the name of the associated
     * equation given at the creation of the tracer
     */

    cs_gwf_tracer_t  *tr = cs_gwf_tracer_by_name("Tracer_01");

    CS_NO_WARN_IF_UNUSED(tr);
  }
  /*! [param_cdo_gwf_get_tracer] */
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
cs_user_parameters(cs_domain_t *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);

  /*! [param_cdo_gwf_get_equation_param_from_decay_chain] */
  {
    /* If the decay chain structure has been added and one wants to access
     * members of this structure, one can retrieve it thanks to
     * cs_gwf_tracer_decay_chain_by_name
     */

    cs_gwf_tracer_decay_chain_t  *tdc =
      cs_gwf_tracer_decay_chain_by_name("my_decay_chain");
    assert(tdc != NULL);

    /* Set the same numerical option for all the tracer equations related to a
       decay chain */

    for (int it = 0; it < tdc->n_tracers; it++) {

      cs_equation_param_t  *tr_eqp =
        cs_gwf_tracer_decay_chain_get_equation_param(tdc, it);

      cs_equation_set_param(tr_eqp, CS_EQKEY_SLES_VERBOSITY, "2");

    }

  }
  /*! [param_cdo_gwf_get_equation_param_from_decay_chain] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Specify the elements such as properties, advection fields,
 *           user-defined equations and modules which have been previously
 *           added.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
*/
/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

  /* 1. Set the Richards equation */
  /* ---------------------------- */

  cs_equation_param_t  *r_eqp = cs_equation_param_by_name("Richards");

  /* Define the boundary conditions  */

  cs_real_t  right_val = 0.0, left_val = 1.0;

  cs_equation_add_bc_by_value(r_eqp,
                              CS_BC_DIRICHLET,
                              "left",       /* boundary zone name */
                              &left_val);   /* value to set */

  cs_equation_add_bc_by_value(r_eqp,
                              CS_BC_DIRICHLET,
                              "right",       /* boundary zone name */
                              &right_val);   /* value to set */

  /* 2. Set the tracer equation */
  /* ---------------------------- */

  cs_equation_param_t  *t_eqp = cs_equation_param_by_name("Tracer_01");

  cs_equation_add_bc_by_value(t_eqp,
                              CS_BC_DIRICHLET,
                              "left",      /* boundary zone name */
                              &left_val);  /* value to set */

  /*! [param_cdo_gwf_set_ic] */
  {
    /* Define the initial condition with an analytic function (if nothing is
       done, the default initialization is zero) */

    cs_equation_add_ic_by_analytic(t_eqp, "cells", get_tracer_ic, NULL);
  }
  /*! [param_cdo_gwf_set_ic] */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
