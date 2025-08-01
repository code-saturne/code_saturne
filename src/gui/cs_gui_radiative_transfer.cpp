/*============================================================================
 * Management of the GUI parameters file: radiative transfer
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_mem.h"
#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "fvm/fvm_selector.h"

#include "base/cs_base.h"
#include "base/cs_boundary_zone.h"
#include "gui/cs_gui_util.h"
#include "gui/cs_gui_boundary_conditions.h"
#include "gui/cs_gui_specific_physics.h"
#include "gui/cs_gui.h"
#include "mesh/cs_mesh.h"
#include "base/cs_post.h"
#include "base/cs_parameters.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_restart.h"
#include "base/cs_thermal_model.h"
#include "base/cs_field.h"
#include "base/cs_field_pointer.h"

#include "rayt/cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "gui/cs_gui_radiative_transfer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure associated to boundary conditions definition
 *----------------------------------------------------------------------------*/

typedef struct {
  int       n_zones;
  char     **label;                /* label for each boundary zone            */
  char     **nature;               /* nature for each boundary zone           */
  int      *type;
  double   *emissivity;
  double   *conductivity;
  double   *thickness;
  double   *thermal_conductivity;
  double   *external_temp;
  double   *internal_temp;
  double   *conduction_flux;
} cs_radiative_transfer_boundary_t;

/*----------------------------------------------------------------------------
 * Private global variables for boundary conditions
 *----------------------------------------------------------------------------*/

static cs_radiative_transfer_boundary_t  *_boundary = nullptr;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return status for logging and postprocessing of radiation
 *
 * parameters:
 *   tn_rt      <-- pointer to tree node associated with radiative transfer
 *   name       <-- name of property
 *   f_log      --> value for logging
 *   f_post_vis --> value for postprocessing
 *
 * return:
 *   label for output
 *----------------------------------------------------------------------------*/

static const char *
_radiative_transfer_char_post(cs_tree_node_t  *tn_rt,
                              const char      *name,
                              int             *f_log,
                              int             *f_post_vis)
{
  cs_tree_node_t *tn = cs_tree_get_node(tn_rt, "property");
  tn = cs_tree_node_get_sibling_with_tag(tn, "name", name);

  const char *label = cs_tree_node_get_tag(tn, "label");

  cs_gui_node_get_child_status_int(tn, "listing_printing", f_log);
  cs_gui_node_get_child_status_int(tn, "postprocessing_recording", f_post_vis);

  if (tn != nullptr && *f_post_vis == -1) /* align with defaults in GUI */
    *f_post_vis = 1;

  return label;
}

/*-----------------------------------------------------------------------------
 * Return value of the type of absorption coefficient for radiation
 *
 * parameters:
 *   tn_rt   <-- tree node associated to radiative transfers
 *   param   <-- name of parameter "absorption coefficient"
 *   keyword --> value of the type of the coefficient
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_type(cs_tree_node_t  *tn_rt,
                         const char      *param,
                         int             *keyword)
{
  cs_tree_node_t *tn = cs_tree_get_node(tn_rt, param);
  const char *type = (tn != nullptr) ? cs_gui_node_get_tag(tn, "type") : nullptr;

  if (type != nullptr) {
    if (cs_gui_strcmp(type, "constant"))
      *keyword = 0;
    else if (cs_gui_strcmp(type, "variable"))
      *keyword = 1;
    else if (cs_gui_strcmp(type, "formula"))
      *keyword = 2;
    else if (cs_gui_strcmp(type, "modak"))
      *keyword = 3;
    else if (cs_gui_strcmp(type, "new_grey_body"))
      *keyword = 4;
    else if (cs_gui_strcmp(type, "spectral_model"))
      *keyword = 5;
    else {
      bft_error (__FILE__, __LINE__, 0,
                 _("unknow type %s\n"), type);
    }
  }
}

/*----------------------------------------------------------------------------
 *  Return int value of the type of radiative condition
 *
 *  parameters:
 *    tn <--tree node associated with this boundary's radiative data
 *
 *  return:
 *    type of radiation condition
 *----------------------------------------------------------------------------*/

static int
_radiative_boundary_type(cs_tree_node_t  *tn_bc)
{
  int result = -999;

  const char *type = cs_tree_node_get_child_value_str(tn_bc, "choice");

  if (cs_gui_strcmp(type, "itpimp"))
    result = CS_BOUNDARY_RAD_WALL_GRAY;
  else if (cs_gui_strcmp(type, "ipgrno"))
    result = CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T;
  else if (cs_gui_strcmp(type, "iprefl"))
    result = CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T;
  else if (cs_gui_strcmp(type, "ifgrno"))
    result = CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX;
  else if (cs_gui_strcmp(type, "ifrefl"))
    result = CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX;

  if (result == -999)
    bft_error (__FILE__, __LINE__, 0,
               _("node request failed \n"));

  return result;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Free GUI-defined radiative transfer parameters
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfers_finalize(void)
{
  if (_boundary != nullptr) {
    int zones = _boundary->n_zones;
    for (int i = 0; i < zones; i++) {
      CS_FREE(_boundary->label[i]);
      CS_FREE(_boundary->nature[i]);
    }
    CS_FREE(_boundary->label);
    CS_FREE(_boundary->nature);
    CS_FREE(_boundary->type);
    CS_FREE(_boundary->emissivity);
    CS_FREE(_boundary->thickness);
    CS_FREE(_boundary->thermal_conductivity);
    CS_FREE(_boundary->external_temp);
    CS_FREE(_boundary->internal_temp);
    CS_FREE(_boundary->conduction_flux);
    CS_FREE(_boundary);
  }
}

/*----------------------------------------------------------------------------
 * Read GUI-defined radiative transfer parameters
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfer_parameters(void)
{
  const char *model = cs_gui_get_thermophysical_model("radiative_transfer");

  if (cs_gui_strcmp(model, "off"))
    cs_glob_rad_transfer_params->type = CS_RAD_TRANSFER_NONE;
  else if (cs_gui_strcmp(model, "dom"))
    cs_glob_rad_transfer_params->type = CS_RAD_TRANSFER_DOM;
  else if (cs_gui_strcmp(model, "p-1"))
    cs_glob_rad_transfer_params->type = CS_RAD_TRANSFER_P1;

  if (cs_glob_rad_transfer_params->type > CS_RAD_TRANSFER_NONE) {

    cs_tree_node_t *tn0
      = cs_tree_get_node(cs_glob_tree,
                         "thermophysical_models/radiative_transfer");

    int isuird = -1;
    cs_gui_node_get_child_status_int(tn0, "restart", &isuird);
    if (! cs_restart_present() || isuird == 0)
      cs_glob_rad_transfer_params->restart = false;
    else if (isuird == 1)
      cs_glob_rad_transfer_params->restart = true;

    cs_gui_node_get_child_int(tn0, "quadrature",
                              &cs_glob_rad_transfer_params->i_quadrature);
    cs_gui_node_get_child_int(tn0, "directions_number",
                              &cs_glob_rad_transfer_params->ndirec);
    cs_gui_node_get_child_int(tn0, "thermal_radiative_source_term",
                              &cs_glob_rad_transfer_params->idiver);
    cs_gui_node_get_child_int(tn0, "temperature_listing_printing",
                              &cs_glob_rad_transfer_params->iimpar);
    cs_gui_node_get_child_int(tn0, "intensity_resolution_listing_printing",
                              &cs_glob_rad_transfer_params->verbosity);
    {
      int ac_type = 0;
      _radiative_transfer_type(tn0, "absorption_coefficient", &ac_type);
      if (ac_type == 3)
        cs_glob_rad_transfer_params->imgrey = 1;
      else if (ac_type == 4)
        cs_glob_rad_transfer_params->imgrey = 2;
      else if (ac_type == 5) {
        if (   (cs_glob_physical_model_flag[CS_COMBUSTION_3PT]  >= 0)
            || (cs_glob_physical_model_flag[CS_COMBUSTION_SLFM] >= 0)) {
          cs_glob_rad_transfer_params->imrcfsk = 1;
          cs_glob_rad_transfer_params->imfsck  = 0;
        } else {
          cs_glob_rad_transfer_params->imrcfsk = 0;
          cs_glob_rad_transfer_params->imfsck  = 1;
        }
      }
    }
    cs_gui_node_get_child_int
      (tn0, "frequency",
       &(cs_glob_rad_transfer_params->time_control.interval_nt));

  }
#if _XML_DEBUG_
  bft_printf("==> cs_gui_radiative_transfer_parameters\n");
  bft_printf("--radiative model: %s  (type = %i)\n",
             model, cs_glob_rad_transfer_params->type);
  if (cs_glob_rad_transfer_params->type > 0) {
    bft_printf("--isuird = %d\n", isuird);
    bft_printf("--quadra = %d\n", cs_glob_rad_transfer_params->i_quadrature);
    bft_printf("--ndirec = %d\n", cs_glob_rad_transfer_params->ndirec);
    bft_printf("--interval_nt = %d\n",
               cs_glob_rad_transfer_params->time_control.interval_nt);
    bft_printf("--idiver = %d\n", cs_glob_rad_transfer_params->idiver);
    bft_printf("--iimpar = %d\n", cs_glob_rad_transfer_params->iimpar);
    bft_printf("--verbosity = %d\n", cs_glob_rad_transfer_params->verbosity);
    bft_printf("--absorption coefficient type: %d\n", ac_type);
    bft_printf("--absorption coefficient by grey body model: %i\n",
               cs_glob_rad_transfer_params->imgrey);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Set the radiative transfer absorption coefficient
 *
 * parameters:
 *   ck --> absorption coefficient at cells
 *----------------------------------------------------------------------------*/

void
cs_gui_rad_transfer_absorption(cs_real_t  ck[])
{

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_tree_node_t *tn
    = cs_tree_get_node(cs_glob_tree,
                       "thermophysical_models/radiative_transfer");
  int ac_type = 0;
  _radiative_transfer_type(tn, "absorption_coefficient", &ac_type);
  if (ac_type == 0) {
    double value = 0.;
    cs_gui_node_get_child_real(tn, "absorption_coefficient", &value);
    for(cs_lnum_t i = 0; i < n_cells; i++)
      ck[i] = value;
  }
#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  bft_printf("--absorption coefficient type: %d\n", ac_type);
  if (ac_type == 0)
    bft_printf("--absorption coefficient value = %f\n", value);
#endif
}

/*----------------------------------------------------------------------------
 * Postprocessing settings for radiative transfer
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfer_postprocess(void)
{
  const int n_rad_b_f = 8;

  const char  *b_rad_names[8] = {
    "rad_incident_flux",
    "spectral_rad_incident_flux",
    "wall_thermal_conductivity",
    "wall_thickness",
    "emissivity",
    "rad_net_flux",
    "rad_convective_flux",
    "rad_exchange_coefficient"};

  cs_field_t * b_rad_f[8] = {
    CS_F_(qinci),
    CS_F_(qinsp),
    CS_F_(xlam),
    CS_F_(epa),
    CS_F_(emissivity),
    CS_F_(fnet),
    CS_F_(fconv),
    CS_F_(hconv)
  };

#if _XML_DEBUG_
  bft_printf("%s\n", __func__);
#endif

  if (cs_glob_rad_transfer_params->type > CS_RAD_TRANSFER_NONE) {

    int k_lbl = cs_field_key_id("label");
    int k_vis = cs_field_key_id("post_vis");
    int k_log = cs_field_key_id("log");

    cs_tree_node_t *tn0
      = cs_tree_get_node(cs_glob_tree,
                         "thermophysical_models/radiative_transfer");

    for (int i = 0; i < n_rad_b_f; i++) {

      cs_field_t *f = b_rad_f[i];
      if (f == nullptr)
        continue;

      int f_post_vis =  0;
      int f_log =  1;
      if (i == 0)
        f_post_vis = CS_POST_ON_LOCATION;
      else
        f_post_vis = -1; /* default not handled as desired in GUI,
                            so handle here */
      const char *label = _radiative_transfer_char_post(tn0,
                                                        b_rad_names[i],
                                                        &f_log,
                                                        &f_post_vis);

#if _XML_DEBUG_
      bft_printf(_("--output boundary faces: %s log %d, postprocess %d\n"),
                 b_rad_names[i], f_log, f_post_vis);
#endif

      if (f_post_vis >= 0)
        cs_field_set_key_int(f, k_vis, f_post_vis);
      if (f_log >= 0)
        cs_field_set_key_int(f, k_log, f_log);
      if (label)
        cs_field_set_key_str(f, k_lbl, label);
    }
  }
}

/*----------------------------------------------------------------------------
 *  Radiative transfer boundary conditions
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfer_bcs(const    int   itypfb[],
                              int           *isothp,
                              double        *epsp,
                              double        *epap,
                              double        *textp,
                              double        *xlamp)
{
  cs_tree_node_t *tn_b0 = cs_tree_get_node(cs_glob_tree,
                                           "boundary_conditions/boundary");
  cs_tree_node_t *tn_w0 = cs_tree_get_node(cs_glob_tree,
                                           "boundary_conditions/wall");

  /* First call only: memory allocation */

  if (_boundary == nullptr) {

    int n_zones = cs_tree_get_node_count(cs_glob_tree,
                                         "boundary_conditions/boundary");

    CS_MALLOC(_boundary, 1, cs_radiative_transfer_boundary_t);
    _boundary->n_zones = n_zones;

    CS_MALLOC(_boundary->label,                n_zones, char *);
    CS_MALLOC(_boundary->nature,               n_zones, char *);
    CS_MALLOC(_boundary->type,                 n_zones, int);
    CS_MALLOC(_boundary->emissivity,           n_zones, double);
    CS_MALLOC(_boundary->thickness,            n_zones, double);
    CS_MALLOC(_boundary->thermal_conductivity, n_zones, double);
    CS_MALLOC(_boundary->external_temp,        n_zones, double);
    CS_MALLOC(_boundary->internal_temp,        n_zones, double);
    CS_MALLOC(_boundary->conduction_flux,      n_zones, double);

    /* loop on boundary zones */

    int z_id = 0;

    for (cs_tree_node_t *tn = tn_b0;
         tn != nullptr;
         tn = cs_tree_node_get_next_of_name(tn), z_id++) {

      /* nature, label and description (color or group)
         of the ith initialization zone */

      const char *nature = cs_tree_node_get_tag(tn, "nature");
      const char *label = cs_tree_node_get_tag(tn, "label");

      CS_MALLOC(_boundary->label[z_id], strlen(label)+1, char);
      strcpy(_boundary->label[z_id], label);

      CS_MALLOC(_boundary->nature[z_id], strlen(nature)+1, char);
      strcpy(_boundary->nature[z_id], nature);

      /* Default initialization */

      _boundary->type[z_id] = -1;
      _boundary->emissivity[z_id] = -1.e12;
      _boundary->thickness[z_id] = -1.e12;
      _boundary->thermal_conductivity[z_id] = -1.e12;
      _boundary->external_temp[z_id] = -1.e12;
      _boundary->internal_temp[z_id] = -1.e12;
      _boundary->conduction_flux[z_id] = 1.e30;

      if (cs_gui_strcmp(nature, "wall")) {

        cs_tree_node_t *tn_w
          = cs_tree_node_get_sibling_with_tag(tn_w0, "label", label);
        cs_tree_node_t *tn_rd = cs_tree_node_get_child(tn_w, "radiative_data");

        _boundary->type[z_id] = _radiative_boundary_type(tn_rd);
        cs_gui_node_get_child_real(tn_rd, "emissivity",
                                   &_boundary->emissivity[z_id]);
        cs_gui_node_get_child_real(tn_rd, "thickness",
                                   &_boundary->thickness[z_id]);
        cs_gui_node_get_child_real(tn_rd, "wall_thermal_conductivity",
                                   &_boundary->thermal_conductivity[z_id]);
        cs_gui_node_get_child_real(tn_rd, "external_temperature_profile",
                                   &_boundary->external_temp[z_id]);
        cs_gui_node_get_child_real(tn_rd, "internal_temperature_profile",
                                   &_boundary->internal_temp[z_id]);
        cs_gui_node_get_child_real(tn_rd, "flux",
                                   &_boundary->conduction_flux[z_id]);

      } /* if (cs_gui_strcmp(nature, "wall")) */

    }  /* for izones */

  }  /* if (_boundary == nullptr)*/

  int z_id = 0;

  for (cs_tree_node_t *tn = tn_b0;
       tn != nullptr;
       tn = cs_tree_node_get_next_of_name(tn), z_id++) {

    const char *label = cs_tree_node_get_tag(tn, "label");

    const cs_zone_t *z = cs_boundary_zone_by_name(label);

    cs_lnum_t n_faces = z->n_elts;
    const cs_lnum_t *faces_list = z->elt_ids;

    if (cs_gui_strcmp(_boundary->nature[z_id], "wall")) {

      cs_field_t *fth = cs_thermal_model_field();
      cs_real_t *th_rcodcl3 = fth->bc_coeffs->rcodcl3;

      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t face_id = faces_list[i];

        if (itypfb[face_id] != CS_SMOOTHWALL && itypfb[face_id] != CS_ROUGHWALL)
          bft_error
            (__FILE__, __LINE__, 0,
             _("Definition of radiative boundary conditions on a boundary "
               "which is not a wall.\n"
               "The definition of the boundary natures given in the GUI"
               " (wall, inlet, outlet, ...)\n"
               "has been modified in a user function"
               " (such as cs_user_boundary_conditions).\n"
               "The radiative boundary conditions defined in the GUI"
               " must be coherent\n"
               "with these new natures.\n"));

        isothp[face_id] = _boundary->type[z_id];
        if (isothp[face_id] == CS_BOUNDARY_RAD_WALL_GRAY) {
          if (_boundary->emissivity[z_id] >= 0.)
            epsp[face_id] = _boundary->emissivity[z_id];
        }
        else if (isothp[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T) {
          if (_boundary->thermal_conductivity[z_id] >= 0.)
            xlamp[face_id] = _boundary->thermal_conductivity[z_id];
          if (_boundary->thickness[z_id] >= 0.)
            epap[face_id] = _boundary->thickness[z_id];

          textp[face_id] = _boundary->external_temp[z_id];
          if (_boundary->emissivity[z_id] >= 0.)
            epsp[face_id] = _boundary->emissivity[z_id];
          if (cs_gui_is_equal_real(_boundary->emissivity[z_id], 0.))
            isothp[face_id] = CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T;
        }
        else if (isothp[face_id] == CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX) {
          th_rcodcl3[face_id]  = _boundary->conduction_flux[z_id];
          if (_boundary->emissivity[z_id] >= 0.)
            epsp[face_id] = _boundary->emissivity[z_id];
          if (cs_gui_is_equal_real(_boundary->emissivity[z_id], 0.))
            isothp[face_id] = CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX;
        }
      }

    } /* if nature == "wall" */

  } /* for z_id */

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  int zones = z_id;
  for (z_id = 0; z_id < zones; z_id++) {
     bft_printf("--label zone = %s\n", _boundary->label[z_id]);
     if (cs_gui_strcmp(_boundary->nature[z_id], "wall")) {
       bft_printf("----type = %i\n", _boundary->type[z_id]);
       bft_printf("----emissivity = %f\n", _boundary->emissivity[z_id]);
       bft_printf("----thickness= %f\n", _boundary->thickness[z_id]);
       bft_printf("----wall_thermal_conductivity = %f\n",
                  _boundary->thermal_conductivity[z_id]);
       bft_printf("----external_temp = %f\n", _boundary->external_temp[z_id]);
       bft_printf("----internal_temp = %f\n", _boundary->internal_temp[z_id]);
       bft_printf("----conduction_flux= %f\n",
                  _boundary->conduction_flux[z_id]);
    }
  }
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
