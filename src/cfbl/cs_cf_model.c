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

#include "bft_mem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cf_model.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_mesh_location.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_physical_properties.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Prototypes for Fortran functions and variables.
 *============================================================================*/

extern int *cs_glob_cf_icvfli = NULL;
extern int *cs_glob_cf_ifbet = NULL;

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
cs_f_cf_model_get_pointers(int    **ieos,
                           int    **ithvar,
                           int    **icfgrp,
                           double **psginf,
                           double **gammasg,
                           int    **hgn_relax_eq_st);

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
 *   ieos             --> pointer to cs_glob_cf_model->ieos
 *   ithvar           --> pointer to cs_glob_cf_model->ithvar
 *   psginf           --> pointer to cs_glob_cf_model->psginf
 *   gammasg          --> pointer to cs_glob_cf_model->gammasg
 *   hgn_relax_eq_st  --> pointer to cs_glob_cf_model->hgn_relax_eq_st
 *----------------------------------------------------------------------------*/

void
cs_f_cf_model_get_pointers(int    **ieos,
                           int    **ithvar,
                           int    **icfgrp,
                           double **psginf,
                           double **gammasg,
                           int    **hgn_relax_eq_st)
{
  *ieos             = &(_cf_model.ieos);
  *ithvar           = &(_cf_model.ithvar);
  *icfgrp           = &(_cf_model.icfgrp);
  *psginf           = &(_cf_model.psginf);
  *gammasg          = &(_cf_model.gammasg);
  *hgn_relax_eq_st  = &(_cf_model.hgn_relax_eq_st);
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
 * \brief Provide access to boundary face indicator array of convection flux
 *        - 0 upwind scheme
 *        - 1 imposed flux
 */
/*----------------------------------------------------------------------------*/

int *
cs_cf_get_icvfli(void)
{
  return cs_glob_cf_icvfli;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Provide access to imposed thermal flux indicator at the boundary
 *        (some boundary contributions of the total energy eq. have to be
 *         cancelled)
 */
/*----------------------------------------------------------------------------*/

int *
cs_cf_get_ifbet(void)
{
  return cs_glob_cf_ifbet;
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

  /* MAP to field pointers */
  cs_field_pointer_map_compressible();
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

END_C_DECLS
