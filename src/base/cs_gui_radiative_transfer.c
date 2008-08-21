/*============================================================================
 *
 *                    Code_Saturne version 1.3
 *                    ------------------------
 *
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

/*============================================================================
 * Reader of the parameters file: radiative transfer
 *============================================================================*/

#if defined(_CS_HAVE_XML)

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
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_gui.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "fvm_selector.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_radiative_transfer.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


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
static int      _cs_gui_last_var = 0;
static char  ** _cs_gui_var_rayt = NULL;

/*============================================================================
 * C Private Functions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return integer parameters for radiation
 *
 *   parameters:
 *    param    -->   name of parameter
 *    keyword  <--   value of parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_radiative_transfer_int(const char *const param,
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

static void cs_gui_radiative_transfer_double(const char   *const param,
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

static void cs_gui_radiative_transfer_char(const char *const param,
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
 * Return status and label of the property for post traitment of radiation
 *
 *   parameters:
 *    name     -->   name of property
 *    value    <--   value of status
 *----------------------------------------------------------------------------*/

static char *cs_gui_radiative_transfer_char_post(const char *const name,
                                                       int  *const value)
{
  char *path = NULL;
  char *path1 = NULL;
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

  cs_xpath_add_attribute(&path, "label");
  cs_xpath_add_attribute(&path1, "status");

  label = cs_gui_get_attribute_value(path);

  if (cs_gui_get_status(path1, &result)) {
    result = (result ? 1 : -1 );
    *value = result;
  }

  BFT_FREE(path);
  BFT_FREE(path1);

  return label;
}

/*-----------------------------------------------------------------------------
 * Retourne  la valeur du type de coeff d'absorption (rayonnement)
 * Return value of the type of absorption coefficient for radiation
 *
 *   parameters:
 *    param    -->   name of parameter "absorption coefficient"
 *    keyword  <--   value of the type of the coefficent
 *----------------------------------------------------------------------------*/

static void cs_gui_radiative_transfer_type(const char *const param,
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

static void cs_gui_radiative_boundary(const   char *const label,
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
                        "wall_radiative_condition",
                        param );
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &res)){
    if (res != *value)
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

static int cs_gui_radiative_boundary_type(const char *const label,
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

  cs_xpath_add_element(&path, "wall_radiative_condition");
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

static int cs_gui_radiative_boundary_output_zone_max(void)
{
  char *path;
  int nb_zone, zone_max = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4,
                        "boundary_conditions",
                        "wall",
                        "wall_radiative_condition",
                        "output_zone" );

  nb_zone = cs_gui_get_nb_element(path);

  if (nb_zone > 0) {
    cs_xpath_add_function_text(&path);
    zone_max = cs_gui_get_max_value(path);
  }

  BFT_FREE(path);

  return zone_max;
}

/*-----------------------------------------------------------------------------
 * Copy a variable name to private variable names array
 *
 * parameters:
 *   varname        -->  name or label of the variable/scalar/property
 *   ipp            -->  index from the fortran array associated to varname
 *----------------------------------------------------------------------------*/


static void _cs_gui_copy_varname(const char *varname,
                                       int   ipp)
{
  size_t  l;

  if (ipp < 1 || ipp > _cs_gui_last_var)
    bft_error(__FILE__, __LINE__, 0,
              _("Variable index %d out of bounds (1 to %d)"),
              ipp, _cs_gui_last_var);

  l = strlen(varname);

  if (_cs_gui_var_rayt[ipp-1] == NULL)
    BFT_MALLOC(_cs_gui_var_rayt[ipp-1], l + 1, char);

  else if (strlen(_cs_gui_var_rayt[ipp-1]) != l)
    BFT_REALLOC(_cs_gui_var_rayt[ipp-1], l + 1, char);

  strcpy(_cs_gui_var_rayt[ipp-1], varname);
}

/*============================================================================
 * C API public functions
 *============================================================================*/

/*============================================================================
 * Fortran API public functions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 *----------------------------------------------------------------------------*/


void CS_PROCF (uiray1, UIRAY1) (int *const nbrayb,
                                int *const nbrayf,
                                int *const nphas,
                                int *const irayon,
                                int *const isuird,
                                int *const ndirec,
                                int *const nfreqr,
                                int *const idiver,
                                int *const iimpar,
                                int *const iimlum,
                                int *const irayvp,
                                int *const irayvf)
{
  int i, ind, iphas = 0;
  char *model = NULL;
  char *label = NULL;

  const char *const _cs_properties_name1[5] = {
    "srad",
    "qrad",
    "absorp",
    "emiss",
    "coefAb" };

  const char *const _cs_properties_name2[8] = {
     "wall_temp",
     "flux_incident",
     "thickness",
     "thermal_conductivity",
     "emissivity",
     "flux_net",
     "flux_convectif",
     "coeff_ech_conv"};

  if (!cs_gui_get_activ_thermophysical_model()) {

    model = cs_gui_get_thermophysical_model("radiative_transfer");

    if (cs_gui_strcmp(model, "off"))
      irayon[iphas] = 0;
    else if (cs_gui_strcmp(model, "dom"))
      irayon[iphas] = 1;
    else if (cs_gui_strcmp(model, "p-1"))
      irayon[iphas] = 2;
  }

  if (irayon[iphas] != 0) {
    cs_gui_radiative_transfer_char("restart", isuird);
    cs_gui_radiative_transfer_int("directions_number", ndirec);
    cs_gui_radiative_transfer_int("frequency", nfreqr);
    cs_gui_radiative_transfer_int("thermal_radiative_source_term", idiver);
    cs_gui_radiative_transfer_int("temperature_listing_printing", iimpar);
    cs_gui_radiative_transfer_int("intensity_resolution_listing_printing", iimlum);

    for (i=0 ; i < *nbrayb; i++) {
      ind = -1;
      label = cs_gui_radiative_transfer_char_post(_cs_properties_name1[i], &ind);
      for (iphas=0 ; iphas < *nphas ; iphas++) {
        irayvp[(*nbrayb)*iphas + i] = ind;
        if (label) _cs_gui_copy_varname(label, i + 1 + (*nbrayb)*iphas);
      }
      BFT_FREE(label);
    }

    for (i=0 ; i < *nbrayf ; i++) {
      ind = -1;
      label = cs_gui_radiative_transfer_char_post(_cs_properties_name2[i], &ind);
      for (iphas=0 ; iphas < *nphas ; iphas++) {
        irayvf[(*nbrayf)*iphas + i] = ind;
        if (label) _cs_gui_copy_varname(label, i + 1 + (*nbrayb)*(*nphas) + (*nbrayf)*iphas);
      }
      BFT_FREE(label);
    }
  }

#if _XML_DEBUG_
  bft_printf("==>UIRAY1\n");
  for (iphas=0 ; iphas < *nphas ; iphas++)
    bft_printf("--rayonnement : %s  (irayon = %d)\n", model, irayon[iphas]);
  bft_printf("--isuird = %d\n", *isuird);
  bft_printf("--ndirec = %d\n", *ndirec);
  bft_printf("--nfreqr = %d\n", *nfreqr);
  bft_printf("--idiver = %d\n", *idiver);
  bft_printf("--iimpar = %d\n", *iimpar);
  bft_printf("--iimlum = %d\n", *iimlum);
  for (i=0 ; i < *nbrayb ; i++) {
    for (iphas=0 ; iphas < *nphas ; iphas++) {
      bft_printf("--output cel: %s value %i \n",
                 _cs_gui_var_rayt[i + (*nbrayb)*iphas],
                 irayvp[(*nbrayb)*iphas + i]);
    }
  }
  for (i=0 ; i < *nbrayf ; i++) {
    for (iphas=0 ; iphas < *nphas ; iphas++) {
      bft_printf("--output fb: %s value %i \n",
      _cs_gui_var_rayt[i + (*nbrayb)*(*nphas) + (*nbrayf)*iphas],
      irayvf[(*nbrayf)*iphas + i]);
    }
  }
#endif

  BFT_FREE(model);
}


/*----------------------------------------------------------------------------
 * Copy variable name from Fortran to C
 *----------------------------------------------------------------------------*/

void CS_PROCF(fcnmra, FCNMRA)
(
 const char          *const fstr,    /* --> Fortran string */
 int                 *const len,     /* --> String Length  */
 int                 *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
)
{
  int i, i1, i2, l;
  char *cstr = NULL;

  assert(*var_id > 0);

  /* Resize array if necessary */

  if (*var_id > _cs_gui_max_vars) {

    if (_cs_gui_max_vars == 0)
      _cs_gui_max_vars = 16;

    while (_cs_gui_max_vars <= *var_id)
      _cs_gui_max_vars *= 2;

    BFT_REALLOC(_cs_gui_var_rayt, _cs_gui_max_vars, char *);
    for (i = _cs_gui_last_var; i < _cs_gui_max_vars; i++)
      _cs_gui_var_rayt[i] = NULL;
  }

  /* Compute string length (removing start or end blanks) */

  for (i1 = 0;
       i1 < *len && (fstr[i1] == ' ' || fstr[i1] == '\t');
       i1++);

  for (i2 = *len - 1;
       i2 > i1 && (fstr[i2] == ' ' || fstr[i2] == '\t');
       i2--);

  l = i2 - i1 + 1;

  /* Should be called once per variable only */
  assert(_cs_gui_var_rayt[*var_id - 1] == NULL);

  if (l > 0) {

    /* Allocate and copy */
    BFT_MALLOC(cstr, l + 1, char);

  for (i = 0 ; i < l ; i++, i1++)
    cstr[i] = fstr[i1];

  cstr[l] = '\0';

    _cs_gui_var_rayt[*var_id - 1] = cstr;

  }

  /* Update variable counter */
  _cs_gui_last_var = *var_id;

}


/*----------------------------------------------------------------------------
 * Copy variable name from C to Fortran
 *----------------------------------------------------------------------------*/


void CS_PROCF(cfnmra, CFNMRA)
(
 char          *const fstr,    /* --> Fortran string */
 int           *const len,     /* --> String Length  */
 int           *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
)
{
  int i;
  int l = 0;
  char *cstr = NULL;

  /* Check that variable name was set */

  if (*var_id < 1 || *var_id > _cs_gui_last_var)
    bft_error(__FILE__, __LINE__, 0,
              _("Name of variable %d was never set.\n"), var_id);

  /* Copy string */

  cstr = _cs_gui_var_rayt[*var_id - 1];

  if (cstr != NULL) {

  /* Compute string length (removing start or end blanks) */

  l = strlen(cstr);
  if (l > *len)
    l = *len;

    for (i = 0; i < l; i++)
      fstr[i] = cstr[i];

  }

  /* Pad with blanks if necessary */

  for (i = l; i < *len; i++)
    fstr[i] = ' ';
}

/*----------------------------------------------------------------------------
 *  Modele de rayonnement usray2.F
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiray2, UIRAY2)
(
 const    int *const itypfb,
 const    int *const iparoi,
 const    int *const iparug,
 const    int *const ivart,
 const    int *const iph,
 const    int *const nphast,
          int *const izfrdp,
          int *const isothp,
 const    int *const itpimp,
 const    int *const ipgrno,
 const    int *const iprefl,
 const    int *const ifgrno,
 const    int *const ifrefl,
 const    int *const nfabor,
 const    int *const nfml,
 const    int *const ifmfbr,
 const    int *const iprfml,
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
  int j, k, n;
  int *faces_list;
  int faces = 0;
  int iok = 0;
  double tmp = 0.;
  char *nature = NULL;
  char *label = NULL;
  char *description = NULL;

  k = (*iph -1) * (*nfabor);

  BFT_MALLOC(faces_list, *nfabor, int);

  zones   = cs_gui_boundary_zones_number();
  output_zone_max = cs_gui_radiative_boundary_output_zone_max();


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

        description = cs_gui_boundary_zone_localization(nature, label);


        fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                                    description,
                                    &faces,
                                    faces_list);

        BFT_FREE(description);

        /*Initialisation par defaut : ces valeurs sont les memes que
        dans raycli, mais elles sont donnees par face dans raycli
        et ici on n'a pas forcement de face de bord (parallele)
          -> on duplique */
        boundary->type[izone] = -1;
        boundary->output_zone[izone] = -1;
        boundary->emissivity[izone] = -1.e12;
        boundary->thickness[izone] = -1.e12;
        boundary->thermal_conductivity[izone] = -1.e12;
        boundary->external_temp[izone] = -1.e12;
        boundary->internal_temp[izone] = -1.e12;
        boundary->conduction_flux[izone] = 1.e30;

        if (cs_gui_strcmp(nature, "wall")) {
          boundary->type[izone] = cs_gui_radiative_boundary_type(label,
                                                                                   *itpimp, *ipgrno, *iprefl,
                                                                                   *ifgrno, *ifrefl);
          tmp = (double) boundary->output_zone[izone];
          cs_gui_radiative_boundary(label, "output_zone", &tmp);
          boundary->output_zone[izone] = (int) tmp;
          cs_gui_radiative_boundary(label, "emissivity", &boundary->emissivity[izone]);
          cs_gui_radiative_boundary(label, "thickness", &boundary->thickness[izone]);
          cs_gui_radiative_boundary(label, "thermal_conductivity", &boundary->thermal_conductivity[izone]);
          cs_gui_radiative_boundary(label, "external_temperature_profile", &boundary->external_temp[izone]);
          cs_gui_radiative_boundary(label, "internal_temperature_profile", &boundary->internal_temp[izone]);
          cs_gui_radiative_boundary(label, "flux", &boundary->conduction_flux[izone]);

        } /* if (cs_gui_strcmp(nature, "wall")) */

        BFT_FREE(nature);
        BFT_FREE(label);

    }  /* for izones */

  }  /* if (boundaries == NULL)*/

  for (izone = 0; izone < zones; izone++) {

    /* list of faces building */

    description = cs_gui_boundary_zone_localization(boundary->nature[izone],
                                                    boundary->label[izone]);

    fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                          description,
                          &faces,
                          faces_list);

    BFT_FREE(description);

    if (cs_gui_strcmp(boundary->nature[izone], "wall")) {

      for (n = 0; n < faces; n++) {
        ifbr = faces_list[n]-1;

        if (itypfb[ifbr] != *iparoi && itypfb[ifbr] != *iparug)
          bft_error(__FILE__, __LINE__, 0,
                    _("One tries to define radiative boundary conditions on boundary which is not a wall.\n"
                      "The definition of the boundaries natures given in GUI (wall, inlet, outlet,...) \n"
                      "is modified in a users subroutine (like USCLIM, USCPCL,...). \n"
                      "The radiative boundary conditions given in GUI must be coherent \n"
                      "with these new natures.\n"));

        izfrdp[k + ifbr] = boundary->output_zone[izone];
        isothp[k + ifbr] = boundary->type[izone];
        if (isothp[k + ifbr] == *itpimp) {
          epsp[k + ifbr] = boundary->emissivity[izone];
          tintp[k + ifbr] = boundary->internal_temp[izone];
        } else if (isothp[k + ifbr] == *ipgrno) {
          epsp[k + ifbr] = boundary->emissivity[izone];
          xlamp[k + ifbr] = boundary->thermal_conductivity[izone];
          epap[k + ifbr] = boundary->thickness[izone];
          textp[k + ifbr] = boundary->external_temp[izone];
          tintp[k + ifbr] = boundary->internal_temp[izone];
        } else if (isothp[k + ifbr] == *iprefl) {
          xlamp[k + ifbr] = boundary->thermal_conductivity[izone];
          epap[k + ifbr] = boundary->thickness[izone];
          textp[k + ifbr] = boundary->external_temp[izone];
          tintp[k + ifbr] = boundary->internal_temp[izone];
        } else if (isothp[k + ifbr] == *ifgrno) {
          epsp[k + ifbr] = boundary->emissivity[izone];
          rcodcl[2 * (*nfabor) * (*nvar) + (*ivart) * (*nfabor) + ifbr] = boundary->conduction_flux[izone];
          tintp[k + ifbr] = boundary->internal_temp[izone];
        } else if (isothp[k + ifbr] == *ifrefl) {
          rcodcl[2 * (*nfabor) * (*nvar) + (*ivart) * (*nfabor) + ifbr] = boundary->conduction_flux[izone];
          tintp[k + ifbr] = boundary->internal_temp[izone];
        }
      }

    } else {
      j = output_zone_max++;
      for (n = 0; n < faces; n++) {
        ifbr = faces_list[n]-1;
        izfrdp[k + ifbr] = j;
      }
    } /* if nature == "wall" */

  } /* for izone */

  BFT_FREE(faces_list);

  iok = 0;
  for (n = 0; n < *nfabor; n++) {
    if (izfrdp[k + n] == -1) iok = 1;
  }
  if (iok == 1) {
    bft_printf(_("Warning: radiative boundary conditions in GUI are not totally defined \n"));
    if (zones)
      bft_printf(_("These are radiative boundary conditions defined in GUI: \n"));
    for (izone = 0; izone < zones; izone++) {
       bft_printf(_("  nature: %s label: %s\n"), boundary->nature[izone], boundary->label[izone]);
       if (cs_gui_strcmp(boundary->nature[izone], "wall")) {
         bft_printf(_("    output_zone = %i\n"), boundary->output_zone[izone]);
         bft_printf(_("    type = %i\n"), boundary->type[izone]);
         bft_printf(_("    emissivity = %f\n"), boundary->emissivity[izone]);
         bft_printf(_("    thickness= %f\n"), boundary->thickness[izone]);
         bft_printf(_("    thermal_conductivity = %f\n"), boundary->thermal_conductivity[izone]);
         bft_printf(_("    external_temp = %f\n"), boundary->external_temp[izone]);
         bft_printf(_("    internal_temp = %f\n"), boundary->internal_temp[izone]);
         bft_printf(_("    conduction_flux= %f\n"), boundary->conduction_flux[izone]);
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
       bft_printf("----thermal_conductivity = %f\n", boundary->thermal_conductivity[izone]);
       bft_printf("----external_temp = %f\n", boundary->external_temp[izone]);
       bft_printf("----internal_temp = %f\n", boundary->internal_temp[izone]);
       bft_printf("----conduction_flux= %f\n", boundary->conduction_flux[izone]);
    }
  }
#endif
}

/*----------------------------------------------------------------------------
 *  Modele de rayonnement usray3.F
 *----------------------------------------------------------------------------*/


void CS_PROCF (uiray3, UIRAY3) (      double *const ck,
                                const    int *const iph,
                                const    int *const ncelet,
                                const    int *const ncel)
{
  double value;
  int i, type;

  if (!cs_gui_get_activ_thermophysical_model()) {

    cs_gui_radiative_transfer_type("absorption_coefficient", &type);
    cs_gui_radiative_transfer_double("absorption_coefficient", &value);

    if (type == 0)
      for(i = 0; i < *ncel; i++)
         ck[(*iph-1) * (*ncelet) + i] = value;

#if _XML_DEBUG_
  bft_printf("==>UIRAY3\n");
  bft_printf("--absorption coefficient type: %d\n", type);
  if (type == 0) bft_printf("--absorption coefficient value = %f\n", value);
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

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CS_HAVE_XML */



