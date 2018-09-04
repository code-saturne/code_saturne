/*============================================================================
 * Management of the GUI parameters file: radiative transfer
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include "cs_field.h"
#include "cs_field_pointer.h"

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
  int      *output_zone;
  int      *type;
  double   *emissivity;
  double   *conductivity;
  double   *thickness;
  double   *thermal_conductivity;
  double   *external_temp;
  double   *internal_temp;
  double   *conduction_flux;
} cs_radiative_boundary_t;

/*----------------------------------------------------------------------------
 * Private global variables for boundary conditions
 *----------------------------------------------------------------------------*/

static cs_radiative_boundary_t *boundary = NULL;

/*----------------------------------------------------------------------------
 * Private global variables for the treatment
 * of NOMVAR. NOMVAR is a characters fortran array
 *----------------------------------------------------------------------------*/

static int      _cs_gui_max_vars = 0;
static char  ** _cs_gui_var_rayt = NULL;


/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return integer parameters for radiation
 *
 *   parameters:
 *    param    -->   name of parameter
 *    keyword  <--   value of parameter
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_int(const char *const param,
                              int *const keyword)
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
 *    param    -->   name of parameter
 *    keyword  <--   value of parameter
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_double(const char   *const param,
                                 double *const keyword)
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
 *    param    -->   name of parameter
 *    keyword  <--   value of parameter
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_char(const char *const param,
                               int  *const keyword)
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
 *   name  --> name of property
 *   value <-- value of status
 *----------------------------------------------------------------------------*/

static char *
_radiative_transfer_char_post(const char *const name,
                              int        *const list_value,
                              int        *const record_value)
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
 *    param    -->   name of parameter "absorption coefficient"
 *    keyword  <--   value of the type of the coefficent
 *----------------------------------------------------------------------------*/

static void
_radiative_transfer_type(const char *const param,
                                     int  *const keyword)
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
 *   parameters:
 *    label    -->   label of boundary nature
 *    param    -->   name of the  variable
 *    value    <--   value of the variable
 *----------------------------------------------------------------------------*/

static void
_radiative_boundary(const   char *const label,
                    const   char *const param,
                          double *const value)
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
 *   parameters:
 *    label    -->   label of boundary "wall"
 *    itpimp   <--   if wall faces with imposed temperature
 *    ipgrno   <--   if grey or black wall faces
 *    iprefl   <--   if reflecting wall faces
 *    ifgrno   <--   if grey or black wall faces and conduction flux imposed
 *    ifrefl   <--   if refecting wall faces and conduction flux imposed
 *----------------------------------------------------------------------------*/

static int
_radiative_boundary_type(const char *const label,
                                 const int itpimp,
                                 const int ipgrno,
                                 const int iprefl,
                                 const int ifgrno,
                                 const int ifrefl)
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
    result = itpimp;
  else if (cs_gui_strcmp(type, "ipgrno"))
    result = ipgrno;
  else if (cs_gui_strcmp(type, "iprefl"))
    result = iprefl;
  else if (cs_gui_strcmp(type, "ifgrno"))
    result = ifgrno;
  else if (cs_gui_strcmp(type, "ifrefl"))
    result = ifrefl;

  if (result == -999)
    bft_error (__FILE__, __LINE__, 0,
               _("Xpath request failed %s \n"), path);

  BFT_FREE(path);
  BFT_FREE(type);

  return result;
}

/*----------------------------------------------------------------------------
 *  Return maximum value of output zone
 *----------------------------------------------------------------------------*/

static int
_radiative_boundary_output_zone_max(void)
{
  char *path;
  int nb_zone, zone_max = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4,
                        "boundary_conditions",
                        "wall",
                        "radiative_data",
                        "output_zone" );

  nb_zone = cs_gui_get_nb_element(path);

  if (nb_zone > 0) {
    cs_xpath_add_function_text(&path);
    zone_max = cs_gui_get_max_value(path);
  }

  BFT_FREE(path);

  return zone_max;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiray1, UIRAY1) (int *const iirayo,
                                int *const isuird,
                                int *const i_quad,
                                int *const ndirec,
                                int *const nfreqr,
                                int *const idiver,
                                int *const iimpar,
                                int *const iimlum,
                                int *const imodak)
{
    char *model = NULL;
    int ac_type = 0;

    model = cs_gui_get_thermophysical_model("radiative_transfer");

    if (cs_gui_strcmp(model, "off"))
        *iirayo = 0;
    else if (cs_gui_strcmp(model, "dom"))
        *iirayo = 1;
    else if (cs_gui_strcmp(model, "p-1"))
        *iirayo = 2;

    if (*iirayo)
    {
        _radiative_transfer_char("restart", isuird);
        _radiative_transfer_int("quadrature", i_quad);
        _radiative_transfer_int("directions_number", ndirec);
        _radiative_transfer_int("frequency", nfreqr);
        _radiative_transfer_int("thermal_radiative_source_term", idiver);
        _radiative_transfer_int("temperature_listing_printing", iimpar);
        _radiative_transfer_int("intensity_resolution_listing_printing", iimlum);
        if (!cs_gui_get_activ_thermophysical_model()) {
            _radiative_transfer_type("absorption_coefficient", &ac_type);
            if (ac_type == 3) {
                *imodak = 1;
            }
         }
    }
#if _XML_DEBUG_
    bft_printf("==>UIRAY1\n");
    bft_printf("--rayonnement : %s  (iirayo = %i)\n", model, *iirayo);
    if (*iirayo)
    {
        bft_printf("--isuird = %d\n", *isuird);
        bft_printf("--quadra = %d\n", *i_quad);
        bft_printf("--ndirec = %d\n", *ndirec);
        bft_printf("--nfreqr = %d\n", *nfreqr);
        bft_printf("--idiver = %d\n", *idiver);
        bft_printf("--iimpar = %d\n", *iimpar);
        bft_printf("--iimlum = %d\n", *iimlum);
        bft_printf("--absorption coefficient type: %d\n", ac_type);
        bft_printf("--absorption coefficient by modak: %i\n", *imodak);
    }
#endif
    BFT_FREE(model);
}

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiray4, UIRAY4) (int *const iirayo)
{
  const int n_rad_b_f = 8;

  const char *const b_rad_names[8] = {
    "wall_temp",
    "flux_incident",
    "thickness",
    "wall_thermal_conductivity",
    "emissivity",
    "flux_net",
    "flux_convectif",
    "coeff_ech_conv"};

  cs_field_t * b_rad_f[8] = {
    CS_F_(tparo),
    CS_F_(qinci),
    CS_F_(epa),
    CS_F_(xlam),
    CS_F_(emissivity),
    CS_F_(fnet),
    CS_F_(fconv),
    CS_F_(hconv)
  };

#if _XML_DEBUG_
  bft_printf("==> UIRAY4\n");
#endif

  if (*iirayo) {
    int k_lbl = cs_field_key_id("label");
    int k_vis = cs_field_key_id("post_vis");
    int k_log = cs_field_key_id("log");
    for (int i = 0; i < n_rad_b_f; i++) {
      int f_post_vis =  1;
      int f_log =  1;
      char *label = _radiative_transfer_char_post(b_rad_names[i],
                                                  &f_log,
                                                  &f_post_vis);
#if _XML_DEBUG_
      bft_printf(_("--output boundary faces: %s log %d, postprocess %d\n"),
                 b_rad_names[i], log, post_vis);
#endif
      cs_field_t *f = b_rad_f[i];
      if (f != NULL) {
        cs_field_set_key_int(f, k_vis, f_post_vis);
        cs_field_set_key_int(f, k_log, f_log);
        if (label)
          cs_field_set_key_str(f, k_lbl, label);
      }
      BFT_FREE(label);
    }
  }
}


/*----------------------------------------------------------------------------
 *  Radiative transfer model usray2.F
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiray2, UIRAY2)
(
 const    int *const itypfb,
 const    int *const iparoi,
 const    int *const iparug,
 const    int *const ivart,
          int *const izfrdp,
          int *const isothp,
 const    int *const itpimp,
 const    int *const ipgrno,
 const    int *const iprefl,
 const    int *const ifgrno,
 const    int *const ifrefl,
 const    int *const nozppm,
 const    int *const nfabor,
 const    int *const nvar,
       double *const epsp,
       double *const epap,
       double *const tintp,
       double *const textp,
       double *const xlamp,
       double *const rcodcl
)
{
  int zones = 0;
  int output_zone_max = 0;
  int izone;
  int ith_zone;
  int ifbr;
  int j, n;
  cs_lnum_t *faces_list = NULL;
  cs_lnum_t faces = 0;
  int iok = 0;
  double tmp = 0.;
  char *nature = NULL;
  char *label = NULL;

  zones   = cs_gui_boundary_zones_number();
  output_zone_max = _radiative_boundary_output_zone_max();

 /* Fisrt iteration only : memory allocation */
  if (boundary == NULL) {

    BFT_MALLOC(boundary,                           1, cs_radiative_boundary_t);
    BFT_MALLOC(boundary->label,                zones, char*                  );
    BFT_MALLOC(boundary->nature,               zones, char*                  );
    BFT_MALLOC(boundary->output_zone,          zones, int                    );
    BFT_MALLOC(boundary->type,                 zones, int                    );
    BFT_MALLOC(boundary->emissivity,           zones, double                 );
    BFT_MALLOC(boundary->thickness,            zones, double                 );
    BFT_MALLOC(boundary->thermal_conductivity, zones, double                 );
    BFT_MALLOC(boundary->external_temp,        zones, double                 );
    BFT_MALLOC(boundary->internal_temp,        zones, double                 );
    BFT_MALLOC(boundary->conduction_flux,      zones, double                 );

    for (izone = 0; izone < zones; izone++) {

    /* nature, label and description (color or group)
       of the ith initialization zone */

        ith_zone = izone + 1;

        nature = cs_gui_boundary_zone_nature(ith_zone);

        label = cs_gui_boundary_zone_label(ith_zone);

        BFT_MALLOC(boundary->label[izone], strlen(label)+1, char);
        strcpy(boundary->label[izone], label);

        BFT_MALLOC(boundary->nature[izone], strlen(nature)+1, char);
        strcpy(boundary->nature[izone], nature);

        /* Default initialization: these values are the same that in raycli
           but given on each face in raycli whereas here one does not
           necessarily have boundary faces (parallism) -> duplication */
        boundary->type[izone] = -1;
        boundary->output_zone[izone] = -1;
        boundary->emissivity[izone] = -1.e12;
        boundary->thickness[izone] = -1.e12;
        boundary->thermal_conductivity[izone] = -1.e12;
        boundary->external_temp[izone] = -1.e12;
        boundary->internal_temp[izone] = -1.e12;
        boundary->conduction_flux[izone] = 1.e30;

        if (cs_gui_strcmp(nature, "wall")) {
          boundary->type[izone] = _radiative_boundary_type(label,
                                                           *itpimp, *ipgrno, *iprefl,
                                                           *ifgrno, *ifrefl);
          tmp = (double) boundary->output_zone[izone];
          _radiative_boundary(label, "output_zone", &tmp);
          boundary->output_zone[izone] = (int) tmp;
          _radiative_boundary(label, "emissivity", &boundary->emissivity[izone]);
          _radiative_boundary(label, "thickness", &boundary->thickness[izone]);
          _radiative_boundary(label, "wall_thermal_conductivity", &boundary->thermal_conductivity[izone]);
          _radiative_boundary(label, "external_temperature_profile", &boundary->external_temp[izone]);
          _radiative_boundary(label, "internal_temperature_profile", &boundary->internal_temp[izone]);
          _radiative_boundary(label, "flux", &boundary->conduction_flux[izone]);

        } /* if (cs_gui_strcmp(nature, "wall")) */

        BFT_FREE(nature);
        BFT_FREE(label);

    }  /* for izones */

  }  /* if (boundaries == NULL)*/

  for (izone = 0; izone < zones; izone++) {

    /* list of faces building */

    faces_list = cs_gui_get_faces_list(izone,
                                       boundaries->label[izone],
                                       *nfabor, *nozppm, &faces);

    if (cs_gui_strcmp(boundary->nature[izone], "wall"))
    {
      for (n = 0; n < faces; n++)
      {
        ifbr = faces_list[n];

        if (itypfb[ifbr] != *iparoi && itypfb[ifbr] != *iparug)
          bft_error(__FILE__, __LINE__, 0,
                    _("One tries to define radiative boundary conditions on boundary which is not a wall.\n"
                      "The definition of the boundaries natures given in GUI (wall, inlet, outlet,...) \n"
                      "is modified in a users subroutine (like cs_user_boundary_conditions, USCPCL,...). \n"
                      "The radiative boundary conditions given in GUI must be coherent \n"
                      "with these new natures.\n"));

        izfrdp[ifbr] = boundary->output_zone[izone];
        isothp[ifbr] = boundary->type[izone];
        if (isothp[ifbr] == *itpimp)
        {
          epsp[ifbr] = boundary->emissivity[izone];
          tintp[ifbr] = boundary->internal_temp[izone];
        }
        else if (isothp[ifbr] == *ipgrno)
        {
          xlamp[ifbr] = boundary->thermal_conductivity[izone];
          epap[ifbr] = boundary->thickness[izone];
          textp[ifbr] = boundary->external_temp[izone];
          tintp[ifbr] = boundary->internal_temp[izone];
          epsp[ifbr] = boundary->emissivity[izone];
          if (cs_gui_is_equal_real(boundary->emissivity[izone], 0.))
            isothp[ifbr] = *iprefl;
        }
        else if (isothp[ifbr] == *ifgrno)
        {
          rcodcl[2 * (*nfabor) * (*nvar) + (*ivart - 1) * (*nfabor) + ifbr]
            = boundary->conduction_flux[izone];
          tintp[ifbr] = boundary->internal_temp[izone];
          if (cs_gui_is_equal_real(boundary->emissivity[izone], 0.))
            isothp[ifbr] = *ifrefl;
          else
            epsp[ifbr] = boundary->emissivity[izone];
        }
      }

    } else {
      j = output_zone_max++;
      for (n = 0; n < faces; n++) {
        ifbr = faces_list[n];
        izfrdp[ifbr] = j;
      }
    } /* if nature == "wall" */

    BFT_FREE(faces_list);

  } /* for izone */

  iok = 0;
  for (n = 0; n < *nfabor; n++) {
    if (izfrdp[n] == -1) iok = 1;
  }
  if (iok == 1) {
    bft_printf("Warning: radiative boundary conditions in GUI are not totally defined \n");
    if (zones)
      bft_printf("These are radiative boundary conditions defined in GUI: \n");
    for (izone = 0; izone < zones; izone++) {
       bft_printf("  nature: %s label: %s\n", boundary->nature[izone], boundary->label[izone]);
       if (cs_gui_strcmp(boundary->nature[izone], "wall")) {
         bft_printf("    output_zone = %i\n", boundary->output_zone[izone]);
         bft_printf("    type = %i\n", boundary->type[izone]);
         bft_printf("    emissivity = %f\n", boundary->emissivity[izone]);
         bft_printf("    thickness= %f\n", boundary->thickness[izone]);
         bft_printf("    wall_thermal_conductivity = %f\n", boundary->thermal_conductivity[izone]);
         bft_printf("    external_temp = %f\n", boundary->external_temp[izone]);
         bft_printf("    internal_temp = %f\n", boundary->internal_temp[izone]);
         bft_printf("    conduction_flux= %f\n", boundary->conduction_flux[izone]);
       }
    }
  }

#if _XML_DEBUG_
  bft_printf("==>UIRAY2\n");
  for (izone = 0; izone < zones; izone++) {
     bft_printf("--label zone = %s\n", boundary->label[izone]);
     if (cs_gui_strcmp(boundary->nature[izone], "wall")) {
       bft_printf("----output_zone = %i\n", boundary->output_zone[izone]);
       bft_printf("----type = %i\n", boundary->type[izone]);
       bft_printf("----emissivity = %f\n", boundary->emissivity[izone]);
       bft_printf("----thickness= %f\n", boundary->thickness[izone]);
       bft_printf("----wall_thermal_conductivity = %f\n", boundary->thermal_conductivity[izone]);
       bft_printf("----external_temp = %f\n", boundary->external_temp[izone]);
       bft_printf("----internal_temp = %f\n", boundary->internal_temp[izone]);
       bft_printf("----conduction_flux= %f\n", boundary->conduction_flux[izone]);
    }
  }
#endif
}

/*----------------------------------------------------------------------------
 *  Radiative transfer model usray3.F
 *----------------------------------------------------------------------------*/


void CS_PROCF (uiray3, UIRAY3) (double   *ck)
{
    int ac_type = 0;
    double value = 0.;

    if (!cs_gui_get_activ_thermophysical_model())
    {
        _radiative_transfer_type("absorption_coefficient", &ac_type);

        if (ac_type == 0)
        {
            cs_lnum_t i;
            cs_lnum_t n_cells = cs_glob_mesh->n_cells;
            _radiative_transfer_double("absorption_coefficient", &value);
            for(i = 0; i < n_cells; i++)
                ck[i] = value;
        }
#if _XML_DEBUG_
    bft_printf("==>UIRAY3\n");
    bft_printf("--absorption coefficient type: %d\n", ac_type);
    if (ac_type == 0)
      bft_printf("--absorption coefficient value = %f\n", value);
#endif
    }
}

/*-----------------------------------------------------------------------------
 * Free memory: clean global private variables.
 *
 * Fortran Interface:
 *
 * SUBROUTINE MEMUI2
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (memui2, MEMUI2) (void)
{
  int zones = 0;
  int i;

  if (boundary != NULL) {

  /* clean memory for global private structure boundaries */

    zones = cs_gui_boundary_zones_number();
    for (i=0 ; i < zones ; i++) {
      BFT_FREE(boundary->label[i]);
      BFT_FREE(boundary->nature[i]);
    }
    BFT_FREE(boundary->label);
    BFT_FREE(boundary->nature);
    BFT_FREE(boundary->output_zone);
    BFT_FREE(boundary->type);
    BFT_FREE(boundary->emissivity);
    BFT_FREE(boundary->thickness);
    BFT_FREE(boundary->thermal_conductivity);
    BFT_FREE(boundary->external_temp);
    BFT_FREE(boundary->internal_temp);
    BFT_FREE(boundary->conduction_flux);
    BFT_FREE(boundary);
  }

  /* clean memory for fortran name of variables */

  for (i = 0; i < _cs_gui_max_vars; i++)
    BFT_FREE(_cs_gui_var_rayt[i]);
  BFT_FREE(_cs_gui_var_rayt);

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
