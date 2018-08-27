/*============================================================================
 * Management of the GUI parameters file: radiative transfer
 *============================================================================*/

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
 * Return integer parameters for radiation
 *
 *   parameters:
 *    param    <--   name of parameter
 *    keyword  -->   value of parameter
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_int(const char  *param,
                              int  *keyword)
{
  char *path;
  int value = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "radiative_transfer",
                        param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_int(path, &value)) *keyword = value;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return float parameters for radiation
 *
 *   parameters:
 *    param    <--   name of parameter
 *    keyword  -->   value of parameter
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_double(const char  *param,
                           double      *keyword)
{
  char *path;
  double value;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "radiative_transfer",
                        param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &value)) *keyword = value;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return value of the parameter of the character type for radiation
 *
 *   parameters:
 *    param    <--   name of parameter
 *    keyword  -->   value of parameter
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_char(const char  *param,
                               int   *keyword)
{
  char *path;
  int result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "radiative_transfer",
                        param);
  cs_xpath_add_attribute(&path, "status");

  if(cs_gui_get_status(path, &result)) *keyword = result;

  BFT_FREE(path);
}


/*-----------------------------------------------------------------------------
 * Return status and label of the property for post treatment of radiation
 *
 * parameters:
 *   name  <-- name of property
 *   value --> value of status
 *----------------------------------------------------------------------------*/

static char *
_radiative_transfer_char_post(const char  *name,
                              int         *list_value,
                              int         *record_value)
{
  char *path = NULL;
  char *path1 = NULL;
  char *path2 = NULL;
  char *label = NULL;
  int result;

  path = cs_xpath_init_path();

  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "radiative_transfer",
                        "property");
  cs_xpath_add_test_attribute(&path, "name", name);

  BFT_MALLOC(path1, strlen(path)+1, char);
  strcpy(path1, path);
  BFT_MALLOC(path2, strlen(path)+1, char);
  strcpy(path2, path);

  cs_xpath_add_attribute(&path, "label");
  label = cs_gui_get_attribute_value(path);

  cs_xpath_add_element(&path1, "listing_printing");
  cs_xpath_add_attribute(&path1, "status");
  if (cs_gui_get_status(path1, &result)) {
    *list_value = 1;
  }

  cs_xpath_add_element(&path2, "postprocessing_recording");
  cs_xpath_add_attribute(&path2, "status");
  if (cs_gui_get_status(path2, &result)) {
    *record_value = -1;
  }

  BFT_FREE(path);
  BFT_FREE(path1);
  BFT_FREE(path2);

  return label;
}

/*-----------------------------------------------------------------------------
 * Return value of the type of absorption coefficient for radiation
 *
 *   parameters:
 *    param    <--   name of parameter "absorption coefficient"
 *    keyword  -->   value of the type of the coefficent
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_type(const char  *param,
                         int         *keyword)
{
  char *path;
  char *type;

  path = cs_xpath_init_path();

  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "radiative_transfer",
                        param);

  cs_xpath_add_attribute(&path, "type");

  type = cs_gui_get_attribute_value(path);

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
    BFT_FREE(type);
  }
  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 *  Return value of radiative variable
 *
 *  parameters:
 *    label <-- label of boundary nature
 *    param <-- name of the  variable
 *    value --> value of the variable
 *----------------------------------------------------------------------------*/

static void
_radiative_boundary(const char  *label,
                    const char  *param,
                    double      *value)
{
  char *path = NULL;
  double res = 0.0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,
                        "boundary_conditions",
                        "wall");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_elements(&path, 2,
                        "radiative_data",
                        param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &res)) {
    if (!cs_gui_is_equal_real(res, *value))
      *value = res;
  }

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 *  Return int value of the type of radiative condition
 *
 *  parameters:
 *    label    <--   label of boundary "wall"
 *    itpimp   <--   if wall faces with imposed temperature
 *    ipgrno   <--   if grey or black wall faces
 *    iprefl   <--   if reflecting wall faces
 *    ifgrno   <--   if grey or black wall faces and conduction flux imposed
 *    ifrefl   <--   if refecting wall faces and conduction flux imposed
 *----------------------------------------------------------------------------*/

static int
_radiative_boundary_type(const char  *label)
{
  char *path = NULL;
  char *type = NULL;
  int result = -999;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2,
                        "boundary_conditions",
                        "wall");
  cs_xpath_add_test_attribute(&path, "label", label);

  cs_xpath_add_element(&path, "radiative_data");
  cs_xpath_add_attribute(&path,"choice");
  type = cs_gui_get_attribute_value(path);

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
               _("Xpath request failed %s \n"), path);

  BFT_FREE(path);
  BFT_FREE(type);

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
    int zones = cs_gui_boundary_zones_number();
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
  if (!cs_gui_file_is_loaded())
    return;

  const char *model = cs_gui_get_thermophysical_model("radiative_transfer");

  int isuird = 0;
  int ac_type = 0;

  if (cs_gui_strcmp(model, "off"))
    cs_glob_rad_transfer_params->type = CS_RAD_TRANSFER_NONE;
  else if (cs_gui_strcmp(model, "dom"))
    cs_glob_rad_transfer_params->type = CS_RAD_TRANSFER_DOM;
  else if (cs_gui_strcmp(model, "p-1"))
    cs_glob_rad_transfer_params->type = CS_RAD_TRANSFER_P1;

  if (cs_glob_rad_transfer_params->type > CS_RAD_TRANSFER_NONE) {
    _radiative_transfer_char("restart", &isuird);
    if (isuird) {
      if (cs_restart_present())
        cs_glob_rad_transfer_params->restart = true;
    }
    _radiative_transfer_int("quadrature",
                            &cs_glob_rad_transfer_params->i_quadrature);
    _radiative_transfer_int("directions_number",
                            &cs_glob_rad_transfer_params->ndirec);
    _radiative_transfer_int("frequency", &cs_glob_rad_transfer_params->nfreqr);
    _radiative_transfer_int("thermal_radiative_source_term",
                            &cs_glob_rad_transfer_params->idiver);
    _radiative_transfer_int("temperature_listing_printing",
                            &cs_glob_rad_transfer_params->iimpar);
    _radiative_transfer_int("intensity_resolution_listing_printing",
                            &cs_glob_rad_transfer_params->iimlum);
    if (!cs_gui_get_activ_thermophysical_model()) {
      _radiative_transfer_type("absorption_coefficient", &ac_type);
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
    _radiative_transfer_type("absorption_coefficient", &ac_type);
    if (ac_type == 0) {
      _radiative_transfer_double("absorption_coefficient", &value);
      for(cs_lnum_t i = 0; i < n_cells; i++)
        ck[i] = value;
    }
#if _XML_DEBUG_
    bft_printf("==>UIRAY3\n");
    bft_printf("--absorption coefficient type: %d\n", ac_type);
    if (type == 0)
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
  if (!cs_gui_file_is_loaded())
    return;

  const int n_rad_b_f = 7;

  const char  *b_rad_names[7] = {
    "flux_incident",
    "thickness",
    "wall_thermal_conductivity",
    "emissivity",
    "flux_net",
    "flux_convectif",
    "coeff_ech_conv"};

  cs_field_t * b_rad_f[7] = {
    CS_F_(qinci),
    CS_F_(epa),
    CS_F_(xlam),
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
    for (int i = 0; i < n_rad_b_f; i++) {
      int f_post_vis =  0;
      int f_log =  1;
      if (i == 0)
        f_post_vis = CS_POST_ON_LOCATION;
      char *label = _radiative_transfer_char_post(b_rad_names[i],
                                                  &f_log,
                                                  &f_post_vis);
#if _XML_DEBUG_
      bft_printf(_("--output boundary faces: %s log %d, postprocess %d\n"),
                 b_rad_names[i], log, post_vis);
#endif
      cs_field_t *f = b_rad_f[i];
      if (f != NULL) {
        if (f_post_vis >= 0)
          cs_field_set_key_int(f, k_vis, f_post_vis);
        if (f_log >= 0)
          cs_field_set_key_int(f, k_log, f_log);
        if (label)
          cs_field_set_key_str(f, k_lbl, label);
      }
      BFT_FREE(label);
    }
  }
}

/*----------------------------------------------------------------------------
 *  Radiative transfer boundary conditions
 *----------------------------------------------------------------------------*/

void
cs_gui_radiative_transfer_bcs(const    int   itypfb[],
                              cs_lnum_t      nvar,
                              cs_lnum_t      ivart,
                              int           *isothp,
                              double        *epsp,
                              double        *epap,
                              double        *tintp,
                              double        *textp,
                              double        *xlamp,
                              double        *rcodcl)
{
  if (!cs_gui_file_is_loaded())
    return;

  int izone;
  int ith_zone;
  int ifbr;
  cs_lnum_t faces = 0;
  char *nature = NULL;
  char *label = NULL;

  const cs_lnum_t  n_b_faces = cs_glob_mesh->n_b_faces;

  int zones = cs_gui_boundary_zones_number();

  /* First call only: memory allocation */

  if (_boundary == NULL) {

    BFT_MALLOC(_boundary, 1, cs_radiative_transfer_boundary_t);
    BFT_MALLOC(_boundary->label,                zones, char *);
    BFT_MALLOC(_boundary->nature,               zones, char *);
    BFT_MALLOC(_boundary->type,                 zones, int);
    BFT_MALLOC(_boundary->emissivity,           zones, double);
    BFT_MALLOC(_boundary->thickness,            zones, double);
    BFT_MALLOC(_boundary->thermal_conductivity, zones, double);
    BFT_MALLOC(_boundary->external_temp,        zones, double);
    BFT_MALLOC(_boundary->internal_temp,        zones, double);
    BFT_MALLOC(_boundary->conduction_flux,      zones, double);

    for (izone = 0; izone < zones; izone++) {

    /* nature, label and description (color or group)
       of the ith initialization zone */

      ith_zone = izone + 1;

      nature = cs_gui_boundary_zone_nature(ith_zone);

      label = cs_gui_boundary_zone_label(ith_zone);

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
        _boundary->type[izone] = _radiative_boundary_type(label);
        _radiative_boundary(label, "emissivity", &_boundary->emissivity[izone]);
        _radiative_boundary(label, "thickness", &_boundary->thickness[izone]);
        _radiative_boundary(label, "wall_thermal_conductivity",
                            &_boundary->thermal_conductivity[izone]);
        _radiative_boundary(label, "external_temperature_profile",
                            &_boundary->external_temp[izone]);
        _radiative_boundary(label, "internal_temperature_profile",
                            &_boundary->internal_temp[izone]);
        _radiative_boundary(label, "flux",
                            &_boundary->conduction_flux[izone]);

      } /* if (cs_gui_strcmp(nature, "wall")) */

      BFT_FREE(nature);
      BFT_FREE(label);

    }  /* for izones */

  }  /* if (boundaries == NULL)*/

  for (izone = 0; izone < zones; izone++) {

    /* list of faces building */

    /* Test for when boundaries are not initialized when only the radiation module
     * is called */
    const cs_lnum_t *faces_list;
    if (boundaries != NULL) {
      faces_list
        = cs_gui_get_boundary_faces(boundaries->label[izone], &faces);
    } else {
      faces_list
        = cs_gui_get_boundary_faces(cs_gui_boundary_zone_label(izone+1), &faces);
    }

    if (cs_gui_strcmp(_boundary->nature[izone], "wall")) {
      for (cs_lnum_t i = 0; i < faces; i++) {
        ifbr = faces_list[i];

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
          rcodcl[2 * n_b_faces * nvar + ivart * n_b_faces + ifbr]
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
  bft_printf("==>UIRAY2\n");
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
