/*============================================================================
 * Unit tests for the solidification module (CDO)
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bft/bft_error.h"
#include "bft/bft_mem_usage.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_mem.h"
#include "base/cs_timer_stats.h"

#include "cdo/cs_advection_field.h"
#include "cdo/cs_equation.h"
#include "cdo/cs_equation_param.h"
#include "cdo/cs_property.h"
#include "cdo/cs_solidification.h"
#include "cdo/cs_thermal_system.h"

/*============================================================================
 * Unit Test Cases
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test 1: Activation of CS_SOLIDIFICATION_MODEL_VELOCITY_GIVEN and
 *        advection field linkage to the thermal equation.
 */
/*----------------------------------------------------------------------------*/

static void
_test_advection_field_linkage(void)
{
  bft_printf("  -> Running Test 1:"
             " Model activation & Advection field linkage...\n");

  /* 1. Create a user advection field */
  cs_adv_field_t *adv
    = cs_advection_field_add("test_adv_flux",
                             CS_ADVECTION_FIELD_USER |
                             CS_ADVECTION_FIELD_TYPE_SCALAR_FLUX);
  assert(adv != nullptr);

  /* 2. Activate solidification with frozen/given velocity */
  cs_solidification_t *solid
    = cs_solidification_activate(CS_SOLIDIFICATION_MODEL_VELOCITY_GIVEN,
                                 0,
                                 CS_SOLIDIFICATION_POST_CELL_STATE);
  assert(solid != nullptr);
  assert(solid->compute_navsto == false);

  /* 3. Retrieve thermal equation */
  cs_equation_t *th_eq = cs_thermal_system_get_equation();
  assert(th_eq != nullptr);
  cs_equation_param_t *th_eqp = cs_equation_get_param(th_eq);
  assert(th_eqp != nullptr);

  /* 4. Associate advection to thermal equation */
  cs_equation_add_advection(th_eqp, adv);

  /* Check that the equation convection flag is enabled and field is correctly
     linked */
  assert(cs_equation_param_has_convection(th_eqp) == true);
  assert(th_eqp->adv_field == adv);

  /* 5. Set physical properties */
  cs_solidification_set_voller_model(273.05, 273.25, 3.33e5);


  /* Cleanup */
  cs_solidification_destroy_all();
  cs_thermal_system_destroy();
  cs_equation_destroy_all();
  cs_advection_field_destroy_all();
  cs_property_destroy_all();

  bft_printf("     [OK] Advection field linkage.\n");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Test 2: Memory lifecycle and non-NavSto safety
 *        Verifies that cs_solidification_destroy_all cleans up completely
 *        without crashing on unallocated Navier-Stokes structures.
 */
/*----------------------------------------------------------------------------*/

static void
_test_lifecycle_and_cleanup(void)
{
  bft_printf("  -> Running Test 2: Memory lifecycle & NavSto decoupling...\n");

  /* Model without velocity field */
  cs_solidification_t *solid
    = cs_solidification_activate(CS_SOLIDIFICATION_MODEL_VOLLER_PRAKASH_87,
                                 CS_SOLIDIFICATION_NO_VELOCITY_FIELD,
                                 CS_SOLIDIFICATION_POST_CELL_STATE |
                                 CS_SOLIDIFICATION_POST_ENTHALPY);
  assert(solid != nullptr);
  assert(solid->compute_navsto == false);
  assert(solid->viscosity == nullptr); // Must be nullptr without NavSto

  /* Example of settings */
  cs_solidification_set_voller_model(273.05, 273.25, 3.33e5);

  /* Teardown */
  cs_solidification_destroy_all();
  assert(cs_solidification_is_activated() == false);
  cs_thermal_system_destroy();
  cs_equation_destroy_all();
  cs_property_destroy_all();

  bft_printf("     [OK] Decoupled lifecycle completed successfully.\n");
}

/*============================================================================
 * Main entry point
 *============================================================================*/

int
main([[maybe_unused]] int   argc,
     [[maybe_unused]] char *argv[])
{
  /* Initialize BFT and memory tracking */
  bft_mem_usage_init();
  cs_mem_init("cs_check_solidification_mem.log");
  cs_timer_stats_initialize();

  bft_printf("\n======================================================\n");
  bft_printf("  CDO Solidification Unit Test Suite\n");
  bft_printf("======================================================\n");

  _test_advection_field_linkage();
  _test_lifecycle_and_cleanup();

  bft_printf("======================================================\n");
  bft_printf("  All solidification unit tests PASSED.\n");
  bft_printf("======================================================\n\n");

  cs_timer_stats_finalize();
  cs_mem_end();
  exit(EXIT_SUCCESS);
}
