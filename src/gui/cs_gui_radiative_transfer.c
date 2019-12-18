/*============================================================================
 * Management of the GUI parameters file: radiative transfer
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

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

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_selector.h"

#include "cs_base.h"
#include "cs_boundary_zone.h"
#include "cs_gui_variables.h"
#include "cs_gui_util.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_gui_specific_physics.h"
#include "cs_gui.h"
#include "cs_mesh.h"
#include "cs_post.h"
#include "cs_parameters.h"
#include "cs_restart.h"
#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_rad_transfer.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_radiative_transfer.h"

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

static cs_radiative_transfer_boundary_t  *_boundary = NULL;

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

  if (*f_post_vis == -1) /* align with defaults in GUI */
    *f_post_vis = 1;

  return label;
}

/*-----------------------------------------------------------------------------
 * Return value of the type of absorption coefficient for radiation
 *
 * parameters:
 *   tn_rt   <-- tree node associated to radiative transfers
 *   param   <-- name of parameter "absorption coefficient"
 *   keyword --> value of the type of the coefficent
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_type(cs_tree_node_t  *tn_rt,
                         const char      *param,
                         int             *keyword)
{
  cs_tree_node_t *tn = cs_tree_get_node(tn_rt, param);
  const char *type = cs_gui_node_get_tag(tn, "type");

  if (type != NULL) {
    if (cs_gui_strcmp(type, "constant"))
      *keyword = 0;
    else if (cs_gui_strcmp(type, "variable"))
      *keyword = 1;
    else if (cs_gui_strcmp(type, "formula"))
      *keyword = 2;
    else if (cs_gui_strcmp(type, "modak"))
      *keyword = 3;
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
    result = cs_glob_rad_transfer_params->itpimp;
  else if (cs_gui_strcmp(type, "ipgrno"))
    result = cs_glob_rad_transfer_params->ipgrno;
  else if (cs_gui_strcmp(type, "iprefl"))
    result = cs_glob_rad_transfer_params->iprefl;
  else if (cs_gui_strcmp(type, "ifgrno"))
    result = cs_glob_rad_transfer_params->ifgrno;
  else if (cs_gui_strcmp(type, "ifrefl"))
    result = cs_glob_rad_transfer_params->ifrefl;

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
  if (_boundary != NULL) {
    int zones = _boundary->n_zones;
    for (int i = 0; i < zones; i++) {
      BFT_FREE(_boundary->label[i]);
      BFT_FREE(_boundary->nature[i]);
    }
    BFT_FREE(_boundary->label);
    BFT_FREE(_boundary->nature);
    BFT_FREE(_boundary->type);
    BFT_FREE(_boundary->emissivity);
    BFT_FREE(_boundary->thickness);
    BFT_FREE(_boundary->thermal_conductivity);
    BFT_FREE(_boundary->external_temp);
    BFT_FREE(_boundary->internal_temp);
    BFT_FREE(_boundary->conduction_flux);
    BFT_FREE(_boundary);
  }
}

/*----------------------------------------------------------------------------
 * Read GUI-defined radiative transfer parameters
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfer_parameters(void)
{
  const char *model = cs_gui_get_thermophysical_model("radiative_transfer");

  int ac_type = 0;

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
    cs_gui_node_get_child_int(tn0, "frequency",
                              &cs_glob_rad_transfer_params->nfreqr);
    cs_gui_node_get_child_int(tn0, "thermal_radiative_source_term",
                              &cs_glob_rad_transfer_params->idiver);
    cs_gui_node_get_child_int(tn0, "temperature_listing_printing",
                              &cs_glob_rad_transfer_params->iimpar);
    cs_gui_node_get_child_int(tn0, "intensity_resolution_listing_printing",
                              &cs_glob_rad_transfer_params->iimlum);
    if (!cs_gui_get_activ_thermophysical_model()) {
      _radiative_transfer_type(tn0, "absorption_coefficient", &ac_type);
      if (ac_type == 3)
        cs_glob_rad_transfer_params->imodak = 1;
    }
  }
#if _XML_DEBUG_
  bft_printf("==> cs_gui_radiative_transfer_parameters\n");
  bft_printf("--radiative model: %s  (type = %i)\n",
             model, cs_glob_rad_transfer_params->type);
  if (cs_glob_rad_transfer_params->type > 0) {
    bft_printf("--isuird = %d\n", isuird);
    bft_printf("--quadra = %d\n", cs_glob_rad_transfer_params->i_quadrature);
    bft_printf("--ndirec = %d\n", cs_glob_rad_transfer_params->ndirec);
    bft_printf("--nfreqr = %d\n", cs_glob_rad_transfer_params->nfreqr);
    bft_printf("--idiver = %d\n", cs_glob_rad_transfer_params->idiver);
    bft_printf("--iimpar = %d\n", cs_glob_rad_transfer_params->iimpar);
    bft_printf("--iimlum = %d\n", cs_glob_rad_transfer_params->iimlum);
    bft_printf("--absorption coefficient type: %d\n", ac_type);
    bft_printf("--absorption coefficient by modak: %i\n",
               cs_glob_rad_transfer_params->imodak);
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
  double value = 0.;
  int ac_type = 0;

  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  if (!cs_gui_get_activ_thermophysical_model()) {

    cs_tree_node_t *tn
      = cs_tree_get_node(cs_glob_tree,
                         "thermophysical_models/radiative_transfer");

    _radiative_transfer_type(tn, "absorption_coefficient", &ac_type);
    if (ac_type == 0) {
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
      if (f == NULL)
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
                              int            nvar,
                              int            ivart,
                              int           *isothp,
                              double        *epsp,
                              double        *epap,
                              double        *tintp,
                              double        *textp,
                              double        *xlamp,
                              double        *rcodcl)
{
  const cs_lnum_t  n_b_faces = cs_glob_mesh->n_b_faces;

  cs_tree_node_t *tn_b0 = cs_tree_get_node(cs_glob_tree,
                                           "boundary_conditions/boundary");
  cs_tree_node_t *tn_w0 = cs_tree_get_node(cs_glob_tree,
                                           "boundary_conditions/wall");

  /* First call only: memory allocation */

  if (_boundary == NULL) {

    int n_zones = cs_tree_get_node_count(cs_glob_tree,
                                         "boundary_conditions/boundary");

    BFT_MALLOC(_boundary, 1, cs_radiative_transfer_boundary_t);
    _boundary->n_zones = n_zones;

    BFT_MALLOC(_boundary->label,                n_zones, char *);
    BFT_MALLOC(_boundary->nature,               n_zones, char *);
    BFT_MALLOC(_boundary->type,                 n_zones, int);
    BFT_MALLOC(_boundary->emissivity,           n_zones, double);
    BFT_MALLOC(_boundary->thickness,            n_zones, double);
    BFT_MALLOC(_boundary->thermal_conductivity, n_zones, double);
    BFT_MALLOC(_boundary->external_temp,        n_zones, double);
    BFT_MALLOC(_boundary->internal_temp,        n_zones, double);
    BFT_MALLOC(_boundary->conduction_flux,      n_zones, double);

    /* loop on boundary zones */

    int izone = 0;

    for (cs_tree_node_t *tn = tn_b0;
         tn != NULL;
         tn = cs_tree_node_get_next_of_name(tn), izone++) {

      /* nature, label and description (color or group)
         of the ith initialization zone */

      const char *nature = cs_tree_node_get_tag(tn, "nature");
      const char *label = cs_tree_node_get_tag(tn, "label");

      BFT_MALLOC(_boundary->label[izone], strlen(label)+1, char);
      strcpy(_boundary->label[izone], label);

      BFT_MALLOC(_boundary->nature[izone], strlen(nature)+1, char);
      strcpy(_boundary->nature[izone], nature);

      /* Default initialization */

      _boundary->type[izone] = -1;
      _boundary->emissivity[izone] = -1.e12;
      _boundary->thickness[izone] = -1.e12;
      _boundary->thermal_conductivity[izone] = -1.e12;
      _boundary->external_temp[izone] = -1.e12;
      _boundary->internal_temp[izone] = -1.e12;
      _boundary->conduction_flux[izone] = 1.e30;

      if (cs_gui_strcmp(nature, "wall")) {

        cs_tree_node_t *tn_w
          = cs_tree_node_get_sibling_with_tag(tn_w0, "label", label);
        cs_tree_node_t *tn_rd = cs_tree_node_get_child(tn_w, "radiative_data");

        _boundary->type[izone] = _radiative_boundary_type(tn_rd);
        cs_gui_node_get_child_real(tn_rd, "emissivity",
                                   &_boundary->emissivity[izone]);
        cs_gui_node_get_child_real(tn_rd, "thickness",
                                   &_boundary->thickness[izone]);
        cs_gui_node_get_child_real(tn_rd, "wall_thermal_conductivity",
                                   &_boundary->thermal_conductivity[izone]);
        cs_gui_node_get_child_real(tn_rd, "external_temperature_profile",
                                   &_boundary->external_temp[izone]);
        cs_gui_node_get_child_real(tn_rd, "internal_temperature_profile",
                                   &_boundary->internal_temp[izone]);
        cs_gui_node_get_child_real(tn_rd, "flux",
                                   &_boundary->conduction_flux[izone]);

      } /* if (cs_gui_strcmp(nature, "wall")) */

    }  /* for izones */

  }  /* if (_boundary == NULL)*/

  int izone = 0;

  for (cs_tree_node_t *tn = tn_b0;
       tn != NULL;
       tn = cs_tree_node_get_next_of_name(tn), izone++) {

    const char *label = cs_tree_node_get_tag(tn, "label");

    const cs_zone_t *z = cs_boundary_zone_by_name(label);

    cs_lnum_t n_faces = z->n_elts;
    const cs_lnum_t *faces_list = z->elt_ids;

    if (cs_gui_strcmp(_boundary->nature[izone], "wall")) {

      cs_lnum_t _nvar = nvar; /* ensure lnum type for multiplication */
      cs_lnum_t _ivart = ivart;

      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t ifbr = faces_list[i];

        if (itypfb[ifbr] != CS_SMOOTHWALL && itypfb[ifbr] != CS_ROUGHWALL)
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

        isothp[ifbr] = _boundary->type[izone];
        if (isothp[ifbr] == cs_glob_rad_transfer_params->itpimp) {
          epsp[ifbr] = _boundary->emissivity[izone];
          tintp[ifbr] = _boundary->internal_temp[izone];
        }
        else if (isothp[ifbr] == cs_glob_rad_transfer_params->ipgrno) {
          xlamp[ifbr] = _boundary->thermal_conductivity[izone];
          epap[ifbr] = _boundary->thickness[izone];
          textp[ifbr] = _boundary->external_temp[izone];
          tintp[ifbr] = _boundary->internal_temp[izone];
          epsp[ifbr] = _boundary->emissivity[izone];
          if (cs_gui_is_equal_real(_boundary->emissivity[izone], 0.))
            isothp[ifbr] = cs_glob_rad_transfer_params->iprefl;
        }
        else if (isothp[ifbr] == cs_glob_rad_transfer_params->ifgrno) {
          rcodcl[2 * n_b_faces*_nvar + _ivart*n_b_faces + ifbr]
            = _boundary->conduction_flux[izone];
          tintp[ifbr] = _boundary->internal_temp[izone];
          if (cs_gui_is_equal_real(_boundary->emissivity[izone], 0.))
            isothp[ifbr] = cs_glob_rad_transfer_params->ifrefl;
          else
            epsp[ifbr] = _boundary->emissivity[izone];
        }
      }

    } /* if nature == "wall" */

  } /* for izone */

#if _XML_DEBUG_
  bft_printf("==> %s\n", __func__);
  int zones = izone;
  for (izone = 0; izone < zones; izone++) {
     bft_printf("--label zone = %s\n", _boundary->label[izone]);
     if (cs_gui_strcmp(_boundary->nature[izone], "wall")) {
       bft_printf("----type = %i\n", _boundary->type[izone]);
       bft_printf("----emissivity = %f\n", _boundary->emissivity[izone]);
       bft_printf("----thickness= %f\n", _boundary->thickness[izone]);
       bft_printf("----wall_thermal_conductivity = %f\n",
                  _boundary->thermal_conductivity[izone]);
       bft_printf("----external_temp = %f\n", _boundary->external_temp[izone]);
       bft_printf("----internal_temp = %f\n", _boundary->internal_temp[izone]);
       bft_printf("----conduction_flux= %f\n",
                  _boundary->conduction_flux[izone]);
    }
  }
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
