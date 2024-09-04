/*============================================================================
 * Compressible models data
 *============================================================================*/

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "cs_array.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_parameters_check.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_physical_properties.h"
#include "cs_prototypes.h"   // for cs_add_model_thermal_field_indexes
#include "cs_restart_default.h"
#include "cs_velocity_pressure.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cf_model.h"
#include "cs_cf_thermo.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cf_model.c
        Compressible models data.
*/
/*----------------------------------------------------------------------------*/

/*!
  \ defgroup compressible Compressible models options

  @addtogroup compressible
  @{

  \struct cs_cf_model_t

  \brief Compressible model general options descriptor

  Members of these fluid properties are publicly accessible, to allow for
  concise syntax, as they are expected to be used in many places.

  \var  cs_cf_model_t::ieos
        indicator of equation of state
        -  CS_EOS_IDEAL_GAS: ideal gas with a constant adiabatic coefficient
        -  CS_EOS_STIFFENED_GAS: stiffened gas
        -  CS_EOS_GAS_MIX: mix of ideal gas
        -  CS_EOS_HOMOGENEOUS_TWO_PHASE: two-phase homogeneous model only,
           each phase follows a stiffened gas law.
        -  CS_EOS_MOIST_AIR: moist air equation of state with condensable
           (mixture of two ideal gas)

  \var  cs_cf_model_t::ithvar
        indicator for thermodynamic variables initialization

  \var  cs_cf_model_t::icfgrp
        indicator for hydrostatic balance in boundary conditions

        In the cases where gravity is predominant, taking into account
        the hydrostatic pressure allows to get rid of perturbations which
        may appear near the horizontal walls when the flow is weakly convective.

        Otherwise, when \ref icfgrp=0, the pressure condition is calculated
        from the solution of the unidimensional Euler equations for a perfect
        gas near a wall, for the variables "normal velocity", "density" and
        "pressure":

        Case of an expansion (M <= 0):
        \f{align*}{
           P_p &= 0 & \textrm{if } 1 + \displaystyle\frac{\gamma-1}{2}M<0
        \\ P_p &= P_i \left(1 + \displaystyle\frac{\gamma-1}{2}M\right)
           ^{\frac{2\gamma}{\gamma-1}} & \textrm{otherwise}
           \f}

        Case of a schock (M > 0):
        \f{eqnarray*}{
        P_p = P_i \left(1 + \displaystyle\frac{\gamma(\gamma+1)}{4}M^2
         +\gamma M \displaystyle\sqrt{1+\displaystyle\frac{(\gamma+1)^2}{16}M^2}\right)
        \f}

        with \f$M = \displaystyle\frac{\vect{u}_i \cdot \vect{n}}{c_i}\f$,
        internal Mach number calculated with the variables taken in the cell.

  \var  cs_cf_model_t::psginf
        stiffened gas limit pressure (zero in perfect gas) for single phase
        model in Pa

  \var  cs_cf_model_t::gammasg
        stiffened gas polytropic coefficient (dimensionless) for single phase
        model

  \defgroup comp_homogeneous Homogeneous two-phase compressible model options

  @addtogroup comp_homogeneous
  @{

  \var  cs_cf_model_t::hgn_relax_eq_st
        source term step indicator for two-phase homogeneous model:
        - -1 disabled
        -  0 enabled

  @}

  @}

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Type definitions
 *============================================================================*/

/* main compressible model structure */

static cs_cf_model_t  _cf_model =
{
  .ieos            = -1,
  .ithvar          = 10000,
  .icfgrp          = 1,
  .psginf          = 0.,
  .gammasg         = 1.4,
  .hgn_relax_eq_st = -1
};

const cs_cf_model_t  *cs_glob_cf_model = &_cf_model;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_cf_model_get_pointers(int    **ithvar);

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get pointers to members of the global compressible model structure.
 *
 * This function is intended for use by Fortran wrappers, and
 * enables mapping to Fortran global pointers.
 *
 * parameters:
 *   ithvar           --> pointer to cs_glob_cf_model->ithvar
 *----------------------------------------------------------------------------*/

void
cs_f_cf_model_get_pointers(int    **ithvar)
{
  *ithvar = &(_cf_model.ithvar);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to compressible model global structure cs_glob_cf_model
 */
/*----------------------------------------------------------------------------*/

cs_cf_model_t *
cs_get_glob_cf_model(void)
{
  return &_cf_model;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Variable field definitions for the compressible module,
 *        according to calculation type selected by the user.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_add_variable_fields(void)
{
  const int kivisl = cs_field_key_id("diffusivity_id");
  const int kvisl0 = cs_field_key_id("diffusivity_ref");

  const cs_real_t epzero = 1e-12;

  /* Set thermal model */
  {
    cs_thermal_model_t *thermal_model = cs_get_glob_thermal_model();
    thermal_model->thermal_variable = CS_THERMAL_MODEL_TOTAL_ENERGY;
    thermal_model->temperature_scale = CS_TEMPERATURE_SCALE_KELVIN;
  }

  /* Total energy */
  {
    cs_field_t *f
      = cs_field_by_id(cs_variable_field_create("total_energy",
                                                "TotEner",
                                                CS_MESH_LOCATION_CELLS,
                                                1));
    cs_add_model_thermal_field_indexes(f->id);

    cs_field_pointer_map(CS_ENUMF_(e_tot), f);

    /* Reference value for diffusivity */
    cs_field_set_key_int (f, kivisl, -1);
    cs_field_set_key_double(f, kvisl0, epzero);
  }

  /* Temperature (postprocessing);
     TODO: should be a property, not a variable */
  {
    cs_field_t *f
      = cs_field_by_id(cs_variable_field_create("temperature",
                                                "TempK",
                                                CS_MESH_LOCATION_CELLS,
                                                1));
    cs_add_model_field_indexes(f->id);

    /* Map to both temperature and secondary t_kelvin pointers */
    cs_field_pointer_map(CS_ENUMF_(t), f);
    cs_field_pointer_map(CS_ENUMF_(t_kelvin),
                       cs_field_by_name_try("temperature"));


    /* Reference value for conductivity */
    cs_field_set_key_int (f, kivisl, -1);
    cs_field_set_key_double(f, kvisl0, epzero);
  }

  /* Mixture fractions (two-phase homogeneous flows) */

  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] == 2) {

    const int keyrf = cs_field_key_id("restart_file");

    const char *f_names[] = {"volume_fraction",
                             "mass_fraction",
                             "energy_fraction"};
    const char *f_labels[] = {"Volume Fraction",
                              "Mass Fraction",
                              "Energy Fraction"};
    const cs_field_pointer_id_t f_pointers[] = {CS_ENUMF_(volume_f),
                                                CS_ENUMF_(mass_f),
                                                CS_ENUMF_(energy_f)};

    /* Volume fraction of phase 1 (with respect to the EOS parameters),
       Mass fraction of phase 1, and
       Energy fraction of phase 1 */

    for (int idx = 0; idx < 3; idx++) {
      cs_field_t *f
        = cs_field_by_id(cs_variable_field_create(f_names[idx],
                                                  f_labels[idx],
                                                  CS_MESH_LOCATION_CELLS,
                                                  1));
      cs_add_model_field_indexes(f->id);

      cs_field_pointer_map(f_pointers[idx], f);

      /* Reference value for diffusivity */
      cs_field_set_key_int (f, kivisl, -1);
      cs_field_set_key_double(f, kvisl0, epzero);

      /* Pure convection equation */
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      eqp->idifft= 0;

      /* Set restart file for fractions */
      cs_field_set_key_int(f, keyrf, CS_RESTART_MAIN);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Property field definitions for the compressible module,
 *        according to calculation type selected by the user.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_add_property_fields(void)
{
  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

  /* Dynamic viscosity of reference of the scalar total energy */

  const int kivisl = cs_field_key_id("diffusivity_id");

  cs_field_t *f_e_tot = cs_field_by_name("total_energy");
  int ifcvsl = cs_field_get_key_int(f_e_tot, kivisl);

  if (ifcvsl < 0 && fp->icv >= 0)
    cs_field_set_key_int(f_e_tot, kivisl, 0);

  /* Property field definitions according to their variability */

  const int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
  const int klbl   = cs_field_key_id("label");
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");

  if (fp->icv >= 0) {
    cs_field_t *f = cs_field_create("specific_heat_const_vol",
                                    field_type,
                                    CS_MESH_LOCATION_CELLS,
                                    1,       /* dim */
                                    false);  /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 0);
    cs_field_set_key_str(f, klbl, "Cv");
    cs_physical_property_define_from_field(f->name, f->type,
                                           f->location_id, f->dim, false);
    fp->icv = f->id;

    cs_field_pointer_map(CS_ENUMF_(cv), f);
  }

  if (fp->iviscv >= 0) {
    cs_field_t *f = cs_field_create("volume_viscosity",
                                    field_type,
                                    CS_MESH_LOCATION_CELLS,
                                    1,       /* dim */
                                    false);  /* has_previous */
    cs_field_set_key_int(f, keyvis, 0);
    cs_field_set_key_int(f, keylog, 0);
    cs_field_set_key_str(f, klbl, "Volume_Viscosity");
    cs_physical_property_define_from_field(f->name, f->type,
                                           f->location_id, f->dim, false);
    fp->iviscv = f->id;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup options specific to the compressible model.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_setup(void)
{
  /* Tranported variables
     -------------------- */

  // Does the temperature scalar behave like a solved temperature
  // (regarding the handling of Cp) ?
  // TODO check this; should be 1 for temperature unless handled in
  // another manner, which migh be the case using Cv instead of Cp...

  const int kscacp  = cs_field_key_id("is_temperature");
  cs_field_set_key_int(cs_field_by_name("temperature"), kscacp, 0);

  // Set upwind convection scheme fo all fields

  const int n_fields = cs_field_n_fields();
  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (   f->type & CS_FIELD_VARIABLE
        && !(f->type & CS_FIELD_CDO)) {
      cs_equation_param_t *eqp = cs_field_get_equation_param(f);
      if (eqp != NULL) {
        eqp->blencv = 0;
      }
    }
  }

  /* Default computation options
     --------------------------- */

  // Variable density

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  fp->irovar = 1;

  /* Parameter checks
     ---------------- */

  cs_parameters_is_equal_int(CS_ABORT_IMMEDIATE,
                             _("Compressible model not compatible "
                               "with pseudo coupled pressure-velocity solver."),
                               "cs_glob_velocity_pressure_param->ipucou",
                               cs_glob_velocity_pressure_param->ipucou,
                               0);

  cs_parameters_is_in_range_int(CS_ABORT_IMMEDIATE,
                             _("Compressible model setup"),
                               "cs_glob_cf_model->icfgrp",
                               cs_glob_cf_model->icfgrp,
                               0, 2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print the compressible module options to setup.log.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_model_log_setup(void)
{
  if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] < 0)
    return;

  cs_log_printf(CS_LOG_SETUP,
                _("\n"
                  "Compressible model\n"
                  "------------------\n\n"));

  const char *icfgrp_value_str[] = {N_("0 (ignored)"),
                                    N_("1 (taken into account)")};

  cs_log_printf(CS_LOG_SETUP,
                _("  Pressure BC with dominant hydrostatic effect:\n"));
  cs_log_printf(CS_LOG_SETUP,
                _("    icfgrp:        %s\n"),
                _(icfgrp_value_str[cs_glob_cf_model->icfgrp]));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize variables of the compressible flow model.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_initialize(void)
{
  /* Compute variable Cv in order to have a correct initialization
     of the total energy (computed in inivar by a call to a thermodynamic
     function), now that initial gas mixture composition is known.
     Note that the only eos with a variable Cv is the ideal gas mix (ieos=3). */

  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  if (fluid_props->icv > -1) {
    cs_real_t *cpro_cp = cs_field_by_id(fluid_props->icp)->val;
    cs_real_t *cpro_cv = cs_field_by_id(fluid_props->icv)->val;
    cs_real_t *mix_mol_mas = cs_field_by_name("mix_mol_mas")->val;

    cs_cf_thermo_cv(cpro_cp, mix_mol_mas, cpro_cv, cs_glob_mesh->n_cells);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute variable physical properties for the  compressible module.
 */
/*----------------------------------------------------------------------------*/

void
cs_cf_physical_properties(void)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  /* Update Lambda/Cv
     ---------------- */

  // It has been checked before this subroutine that cv0 was non zero.
  // If Cv is variable and zero, it is an error due to the user.
  // Here a test is performed at each call (not optimal).
  // If the diffusivity of the total energy is constant, then the thermal
  // conductivity and the isochoric specific heat should be constant.

  const cs_fluid_properties_t *fluid_props = cs_glob_fluid_properties;
  const int kivisl  = cs_field_key_id("diffusivity_id");
  int ifcven = cs_field_get_key_int(CS_F_(e_tot), kivisl);

  if (ifcven >= 0) {

    cs_real_t *cpro_venerg = cs_field_by_id(ifcven)->val;

    int ifclam = cs_field_get_key_int(CS_F_(t), kivisl);
    if (ifclam >= 0) {
      const cs_real_t *cpro_lambda = cs_field_by_id(ifclam)->val;
      cs_array_real_copy(n_cells, cpro_lambda, cpro_venerg);
    }
    else {
      const int kvisl0 = cs_field_key_id("diffusivity_ref");
      double visls_0 = cs_field_get_key_double(CS_F_(t), kvisl0);
      cs_array_real_set_scalar(n_cells, visls_0, cpro_venerg);
    }

    if (fluid_props->icv > -1) {
      cs_real_t *cpro_cp = cs_field_by_id(fluid_props->icp)->val;
      cs_real_t *cpro_cv = cs_field_by_id(fluid_props->icv)->val;
      cs_real_t *mix_mol_mas = cs_field_by_name("mix_mol_mas")->val;

      cs_cf_thermo_cv(cpro_cp, mix_mol_mas, cpro_cv, n_cells);

#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        if (cpro_cv[c_id] <= 0)
          bft_error(__FILE__, __LINE__, 0,
                    _("The isochoric specific heat has at least one\n"
                      " negative or zero value: %g."), cpro_cv[c_id]);
        cpro_venerg[c_id] /= cpro_cv[c_id];
      }
    }
    else {
      cs_real_t cv0 = fluid_props->cv0;
#     pragma omp parallel for if (n_cells > CS_THR_MIN)
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_venerg[c_id] /= cv0;
      }
    }
  }

  else {
    // TODO: this part should be done at setup time,
    // instead of modifying field keywards in the time loop
    // (i.e. after setup logging), which is ugly and risky.

    const int kvisl0 = cs_field_key_id("diffusivity_ref");
    double visls_0 = cs_field_get_key_double(CS_F_(t), kvisl0);
    visls_0 /= fluid_props->cv0;
    cs_field_set_key_double(CS_F_(e_tot), kvisl0, visls_0);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
