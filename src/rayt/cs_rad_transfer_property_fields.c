/*============================================================================
 * Radiation solver operations.
 *============================================================================*/

/* This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_parameters.h"
#include "cs_post.h"

#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_rad_transfer_property_fields.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*  \file cs_rad_transfer_property_field.c    */
/*  \brief Properties definitions for radiative model.     */


/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create property fields for radiative solver
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_prp(void)
{
  const int keylbl = cs_field_key_id("label");
  const int keyvis = cs_field_key_id("post_vis");
  const int keylog = cs_field_key_id("log");

  if (cs_glob_rad_transfer_params->type > CS_RAD_TRANSFER_NONE) {

    cs_field_t *f = NULL;

    int field_type = CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY;
    int location_id = CS_MESH_LOCATION_CELLS;

    {
      f = cs_field_create("luminance",
                          field_type,
                          location_id,
                          1,
                          false);

      /* hide property */
      cs_field_set_key_int(f, keyvis, 0);
      cs_field_set_key_int(f, keylog, 0);

      /* set label */
      cs_field_set_key_str(f, keylbl, "Luminance");

      cs_field_pointer_map(CS_ENUMF_(rad_lumin), f);
    }

    {
      f = cs_field_create("radiative_flux",
                          field_type,
                          location_id,
                          3,
                          false);

      /* hide property */
      cs_field_set_key_int(f, keyvis, 0);
      cs_field_set_key_int(f, keylog, 0);

      /* set label */
      cs_field_set_key_str(f, keylbl, "Qrad");

      cs_field_pointer_map(CS_ENUMF_(rad_q), f);
    }

    for (int irphas = 0;
         irphas < cs_glob_rad_transfer_params->nrphas;
         irphas++) {

      char suffix[16];

      if (irphas > 0)
        snprintf(suffix,  15, "_%02d", irphas+1);
      else
        suffix[0]  = '\0';

      suffix[15] = '\0';

      char f_name[64], f_label[64];
      {
        snprintf(f_name, 63, "rad_st%s", suffix); f_name[63] = '\0';
        snprintf(f_label, 63, "Srad%s", suffix); f_label[63]  = '\0';
        f = cs_field_create(f_name,
                            field_type,
                            location_id,
                            1,
                            false);

        /* hide property */
        cs_field_set_key_int(f, keyvis, 0);
        cs_field_set_key_int(f, keylog, 0);

        /* set label */
        cs_field_set_key_str(f, keylbl, f_label);

        cs_field_pointer_map_indexed(CS_ENUMF_(rad_est), irphas, f);
      }

      {
        snprintf(f_name, 63, "rad_st_implicit%s", suffix); f_name[63] = '\0';
        snprintf(f_label, 63, "ITSRI%s", suffix); f_label[63]  = '\0';
        f = cs_field_create(f_name,
                            field_type,
                            location_id,
                            1,
                            false);

        /* hide property */
        cs_field_set_key_int(f, keyvis, 0);
        cs_field_set_key_int(f, keylog, 0);

        /* set label */
        cs_field_set_key_str(f, keylbl, f_label);

        cs_field_pointer_map_indexed(CS_ENUMF_(rad_ist), irphas, f);
      }

      {
        snprintf(f_name, 63, "rad_absorption%s", suffix); f_name[63] = '\0';
        snprintf(f_label, 63, "Absorp%s", suffix); f_label[63]  = '\0';
        f = cs_field_create(f_name,
                            field_type,
                            location_id,
                            1,
                            false);

        /* hide property */
        cs_field_set_key_int(f, keyvis, 0);
        cs_field_set_key_int(f, keylog, 0);

        /* set label */
        cs_field_set_key_str(f, keylbl, f_label);

        cs_field_pointer_map_indexed(CS_ENUMF_(rad_abs), irphas, f);
      }

      {
        snprintf(f_name, 63, "rad_emission%s", suffix); f_name[63] = '\0';
        snprintf(f_label, 63, "Emiss%s", suffix); f_label[63]  = '\0';
        f = cs_field_create(f_name,
                            field_type,
                            location_id,
                            1,
                            false);

        /* hide property */
        cs_field_set_key_int(f, keyvis, 0);
        cs_field_set_key_int(f, keylog, 0);

        /* set label */
        cs_field_set_key_str(f, keylbl, f_label);

        cs_field_pointer_map_indexed(CS_ENUMF_(rad_emi), irphas, f);
      }

      {
        snprintf(f_name, 63, "rad_absorption_coeff%s", suffix); f_name[63] = '\0';
        snprintf(f_label, 63, "CoefAb%s", suffix); f_label[63]  = '\0';
        f = cs_field_create(f_name,
                            field_type,
                            location_id,
                            1,
                            false);
        /* hide property */
        cs_field_set_key_int(f, keyvis, 0);
        cs_field_set_key_int(f, keylog, 0);

        /* set label */
        cs_field_set_key_str(f, keylbl, f_label);

        cs_field_pointer_map_indexed(CS_ENUMF_(rad_cak), irphas, f);
      }

      /* Fields for atmospheric infrared absorption model */

      if (cs_glob_rad_transfer_params->atmo_ir_absorption) {

        const char *name[4] = {"rad_flux_up",
                               "rad_flux_down",
                               "rad_absorption_coeff_up",
                               "rad_absorption_coeff_down"};
        const char *label[4] = {"Upwards radiative flux",
                                "Downwards radiative flux",
                                "",
                                ""};
        const bool hide[4] = {false, false, true, true};

        for (int i = 0; i < 4; i++) {

          f = cs_field_create(name[i],
                              field_type,
                              location_id,
                              1,
                              false);

          /* show or hide property */
          if (!hide[i])
            cs_field_set_key_int(f, keyvis, 1);
          cs_field_set_key_int(f, keylog, 1);

          /* set label */
          if (strlen(label[i]) > 0)
            cs_field_set_key_str(f, keylbl, label[i]);

        }

      }

    }

    /* Boundary face fields */
    location_id = CS_MESH_LOCATION_BOUNDARY_FACES;

    {
      f = cs_field_by_name_try("boundary_temperature");
      if (f == NULL)
        f = cs_parameters_add_boundary_temperature();

      if (!cs_field_is_key_set(f, keylog))
        cs_field_set_key_int(f, keylog, 1);
      if (!cs_field_is_key_set(f, keyvis))
        cs_field_set_key_int(f, keyvis, CS_POST_ON_LOCATION);
    }

    {
      f = cs_field_create("rad_incident_flux",
                          field_type,
                          location_id,
                          1,
                          false);
      cs_field_set_key_str(f, keylbl, "Incident_flux");
      cs_field_pointer_map(CS_ENUMF_(qinci), f);

    }

    {
      if (   cs_glob_rad_transfer_params->imoadf >= 1
          || cs_glob_rad_transfer_params->imfsck == 1) {
        f = cs_field_create("spectral_rad_incident_flux",
                            field_type,
                            location_id,
                            cs_glob_rad_transfer_params->nwsgg,
                            false);
        cs_field_set_key_str(f, keylbl, "Spectral_incident_flux");
        cs_field_pointer_map(CS_ENUMF_(qinsp), f);
      }
    }

    {
      f = cs_field_create("wall_thermal_conductivity",
                          field_type,
                          location_id,
                          1,
                          false);
      cs_field_set_key_str(f, keylbl, "Th_conductivity");
      cs_field_pointer_map(CS_ENUMF_(xlam), f);
    }

    {
      f = cs_field_create("wall_thickness",
                          field_type,
                          location_id,
                          1,
                          false);
      cs_field_set_key_str(f, keylbl, "Thickness");
      cs_field_pointer_map(CS_ENUMF_(epa), f);
    }

    {
      f = cs_field_create("emissivity",
                          field_type,
                          location_id,
                          1,
                          false);
      cs_field_set_key_str(f, keylbl, "Emissivity");
      cs_field_pointer_map(CS_ENUMF_(emissivity), f);
    }

    {
      f = cs_field_create("rad_net_flux",
                          field_type,
                          location_id,
                          1,
                          false);
      cs_field_set_key_str(f, keylbl, "Net_flux");
      cs_field_pointer_map(CS_ENUMF_(fnet), f);
    }

    {
      f = cs_field_create("rad_convective_flux",
                          field_type,
                          location_id,
                          1,
                          false);
      cs_field_set_key_str(f, keylbl, "Convective_flux");
      cs_field_pointer_map(CS_ENUMF_(fconv), f);
    }

    {
      f = cs_field_create("rad_exchange_coefficient",
                          field_type,
                          location_id,
                          1,
                          false);
      cs_field_set_key_str(f, keylbl, "Convective_exch_coef");
      cs_field_pointer_map(CS_ENUMF_(hconv), f);
    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
