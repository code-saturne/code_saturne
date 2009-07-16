/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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
 * Management of the GUI parameters file: particles tracking
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "cs_gui.h"
#include "cs_gui_util.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui_particles.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* debugging switch */
#define _XML_DEBUG_ 0

/*============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structures associated to lagrangian particles definition
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Private global variables for lagrangian particles
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Private global variables for the treatment
 * of NOMLAG, NOMLAV and NOMBRD (characters fortran arrays).
 *----------------------------------------------------------------------------*/

static int      _max_mean_vars = 0;
static int      _last_mean_var = 0;
static char  ** _array_mean_varname = NULL;

static int      _max_variance_vars = 0;
static int      _last_variance_var = 0;
static char  ** _array_variance_varname = NULL;

static int      _max_boundary_vars = 0;
static int      _last_boundary_var = 0;
static char  ** _array_boundary_varname = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return value of the particles model
 *----------------------------------------------------------------------------*/

static void
_get_particles_model(int *const imodel)
{
  char *path;
  char *model;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "lagrangian", "particles_models");
  cs_xpath_add_attribute(&path, "model");
  model = cs_gui_get_attribute_value(path);

  if (model != NULL) {
    if (cs_gui_strcmp(model, "off"))
      *imodel = 0;
    else if (cs_gui_strcmp(model, "one_way"))
      *imodel = 1;
    else if (cs_gui_strcmp(model, "two_way"))
      *imodel = 2;
    else if (cs_gui_strcmp(model, "frozen_fields"))
      *imodel = 3;
    else if (cs_gui_strcmp(model, "thermal"))
      *imodel = 1;
    else if (cs_gui_strcmp(model, "coal"))
      *imodel = 2;
    BFT_FREE(model);
  }
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return value of the parameter of the character type for lagrangian
 *
 *   parameters:
 *   keyword   <--   value of parameter
 *   nbr       -->   size of the labels list
 *   ...       -->   list of labels in the path
 *----------------------------------------------------------------------------*/

static void
_get_status(int *const keyword, const int nbr, ...)
{
  va_list list;

  char *elt = NULL;
  char *path;
  int i;
  int result;

  path = cs_xpath_init_path();

  va_start(list, nbr);

  for(i=0; i<nbr; i++) {

    elt = va_arg(list, char *);

    if (elt != NULL) {

      BFT_REALLOC(path,
                  strlen(path)+ strlen(elt)+ strlen("/") +1,
                  char);

      strcat(path, "/");
      strcat(path, elt);
    }
  }
  va_end(list);

  cs_xpath_add_attribute(&path, "status");
  if(cs_gui_get_status(path, &result))
    *keyword = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return integer parameters for lagrangian
 *
 *   parameters:
 *   keyword   <--   value of parameter
 *   nbr       -->   size of the labels list
 *   ...       -->   list of labels in the path
 *----------------------------------------------------------------------------*/

static void
_get_int(int *const keyword, const int nbr, ...)
{
  va_list list;

  char *elt = NULL;
  char *path;
  int value = 0;
  int i;

  path = cs_xpath_init_path();

  va_start(list, nbr);

  for(i=0; i<nbr; i++) {

    elt = va_arg(list, char *);

    if (elt != NULL) {

      BFT_REALLOC(path,
                  strlen(path)+ strlen(elt)+ strlen("/") +1,
                  char);

      strcat(path, "/");
      strcat(path, elt);
    }
  }
  va_end(list);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_int(path, &value))
    *keyword = value;

  BFT_FREE(path);

}


/*-----------------------------------------------------------------------------
 * Return float parameters for lagrangian
 *
 *   parameters:
 *   keyword   <--   value of parameter
 *   nbr       -->   size of the labels list
 *   ...       -->   list of labels in the path
 *----------------------------------------------------------------------------*/

static void
_get_double(double *const keyword, const int nbr, ...)
{
  va_list list;

  char *elt = NULL;
  char *path;
  double value = 0;
  int i;

  path = cs_xpath_init_path();

  va_start(list, nbr);

  for(i=0; i<nbr; i++) {

    elt = va_arg(list, char *);

    if (elt != NULL) {

      BFT_REALLOC(path,
                  strlen(path)+ strlen(elt)+ strlen("/") +1,
                  char);

      strcat(path, "/");
      strcat(path, elt);
    }
  }
  va_end(list);

  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &value))
    *keyword = value;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return value of the attribute of the character type for larangian
 *
 *   parameters:
 *   param     <--   name of the attribute
 *   nbr       -->   size of the labels list
 *   ...       -->   list of labels in the path
 *----------------------------------------------------------------------------*/

static char*
_get_attr(const char *const param, const int nbr, ...)
{
  va_list list;

  int i;
  char *elt = NULL;
  char *path;
  char *name;

  path = cs_xpath_init_path();

  va_start(list, nbr);

  for(i=0; i<nbr; i++) {

    elt = va_arg(list, char *);

    if (elt != NULL) {

      BFT_REALLOC(path,
                  strlen(path)+ strlen(elt)+ strlen("/") +1,
                  char);

      strcat(path, "/");
      strcat(path, elt);
    }
  }
  va_end(list);

  cs_xpath_add_attribute(&path, param);

  name = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return name;
}

/*-----------------------------------------------------------------------------
 * Return float parameters for coal parameters
 *
 *   parameters:
 *    param         -->   value to modify
 *    name          -->   name of property
 *    icoal         -->   number of coal
 *----------------------------------------------------------------------------*/

static void
_get_coal_double(double *const param, const char *const name, int icoal)
{
  double result = 0;
  char *path = NULL;
  char scoal[2];

  sprintf(scoal, "%i", icoal);

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4, "lagrangian", "particles_models", "coal_fouling", name);
  cs_xpath_add_test_attribute(&path, "coal", scoal);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *param = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return status and label of the property for post treatment
 *
 *   parameters:
 *    type          -->   type of property ('volume' or 'boundary')
 *    name          -->   name of property
 *    list_value    <--   status for listing
 *    record_value  <--   status for post processing
 *----------------------------------------------------------------------------*/

static char*
_get_char_post(const char *const type,
               const char *const name,
               int  *list_value,
               int  *record_value)
{
  char *path, *path1, *path2 = NULL;
  char *label = NULL;
  int result;

  *list_value = 1;
  *record_value = 1;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4, "lagrangian", "statistics", type, "property");
  cs_xpath_add_test_attribute(&path, "name", name);
  BFT_MALLOC(path1, strlen(path)+1, char);
  strcpy(path1, path);
  BFT_MALLOC(path2, strlen(path)+1, char);
  strcpy(path2, path);
  cs_xpath_add_attribute(&path, "label");
  label = cs_gui_get_attribute_value(path);

  if (cs_gui_strcmp(type, "volume")) {

    cs_xpath_add_element(&path1, "monitoring_point");
    cs_xpath_add_attribute(&path1, "status");
    if (cs_gui_get_status(path1, &result))
      *list_value = result;
  }

  else if (cs_gui_strcmp(type, "boundary")) {

    cs_xpath_add_element(&path1, "listing_printing");
    cs_xpath_add_attribute(&path1, "status");
    if (cs_gui_get_status(path1, &result))
      *list_value = result;

    cs_xpath_add_element(&path2, "postprocessing_recording");
    cs_xpath_add_attribute(&path2, "status");
    if (cs_gui_get_status(path2, &result))
      *record_value = -1;
  }

  BFT_FREE(path);
  BFT_FREE(path1);
  BFT_FREE(path2);

  return label;
}

/*-----------------------------------------------------------------------------
 * Copy a variable name to the variable names array
 *
 * parameters:
 *   varname        -->  name or label of the variable/scalar/property
 *   ipp            -->  index from the fortran array associated to varname
 *----------------------------------------------------------------------------*/

static void
_copy_mean_varname(char *varname, int ipp)
{
  int i;
  size_t  l;

  assert(ipp > 0);

  /* Resize array if necessary */

  if (ipp > _max_mean_vars) {

    if (_max_mean_vars == 0)
      _max_mean_vars = 16;

    while (_max_mean_vars <= ipp)
      _max_mean_vars *= 2;

    BFT_REALLOC(_array_mean_varname, _max_mean_vars, char *);
    for (i = _last_mean_var; i < _max_mean_vars; i++)
      _array_mean_varname[i] = NULL;
  }

  if (ipp < 1 || ipp > _last_mean_var+1)
    bft_error(__FILE__, __LINE__, 0,
              _("Variable index %d out of bounds (1 to %d)"),
              ipp, _last_mean_var);

  l = strlen(varname);

  if (_array_mean_varname[ipp-1] == NULL)
    BFT_MALLOC(_array_mean_varname[ipp-1], l + 1, char);

  else if (strlen(_array_mean_varname[ipp-1]) != l)
    BFT_REALLOC(_array_mean_varname[ipp-1], l + 1, char);

  strcpy(_array_mean_varname[ipp-1], varname);

  _last_mean_var = ipp;

}

/*-----------------------------------------------------------------------------
 * Copy a variable name to the variance variable names array
 *
 * parameters:
 *   varname        -->  name or label of the variable/scalar/property
 *   ipp            -->  index from the fortran array associated to varname
 *----------------------------------------------------------------------------*/

static void
_copy_variance_varname(char *varname, int  ipp)
{
  int i;
  size_t  l;

  assert(ipp > 0);

  /* Resize array if necessary */

  if (ipp > _max_variance_vars) {

    if (_max_variance_vars == 0)
      _max_variance_vars = 16;

    while (_max_variance_vars <= ipp)
      _max_variance_vars *= 2;

    BFT_REALLOC(_array_variance_varname, _max_variance_vars, char *);
    for (i = _last_variance_var; i < _max_variance_vars; i++)
      _array_variance_varname[i] = NULL;
  }

  if (ipp < 1 || ipp > _last_variance_var+1)
    bft_error(__FILE__, __LINE__, 0,
              _("Variable index %d out of bounds (1 to %d)"),
              ipp, _last_variance_var);

  l = strlen(varname);

  if (_array_variance_varname[ipp-1] == NULL)
    BFT_MALLOC(_array_variance_varname[ipp-1], l + 1, char);

  else if (strlen(_array_variance_varname[ipp-1]) != l)
    BFT_REALLOC(_array_variance_varname[ipp-1], l + 1, char);

  strcpy(_array_variance_varname[ipp-1], varname);

  _last_variance_var = ipp;

}

/*-----------------------------------------------------------------------------
 * Copy a variable name to the variance variable names array
 *
 * parameters:
 *   varname        -->  name or label of the variable/scalar/property
 *   ipp            -->  index from the fortran array associated to varname
 *----------------------------------------------------------------------------*/

static void
_copy_boundary_varname(char *varname, int  ipp)
{
  int i;
  size_t  l;

  assert(ipp > 0);

  /* Resize array if necessary */

  if (ipp > _max_boundary_vars) {

    if (_max_boundary_vars == 0)
      _max_boundary_vars = 16;

    while (_max_boundary_vars <= ipp)
      _max_boundary_vars *= 2;

    BFT_REALLOC(_array_boundary_varname, _max_boundary_vars, char *);
    for (i = _last_boundary_var; i < _max_boundary_vars; i++)
      _array_boundary_varname[i] = NULL;
  }

  if (ipp < 1 || ipp > _last_boundary_var+1)
    bft_error(__FILE__, __LINE__, 0,
              _("Variable index %d out of bounds (1 to %d)"),
              ipp, _last_boundary_var);

  l = strlen(varname);

  if (_array_boundary_varname[ipp-1] == NULL)
    BFT_MALLOC(_array_boundary_varname[ipp-1], l + 1, char);

  else if (strlen(_array_boundary_varname[ipp-1]) != l)
    BFT_REALLOC(_array_boundary_varname[ipp-1], l + 1, char);

  strcpy(_array_boundary_varname[ipp-1], varname);

  _last_boundary_var = ipp;

}

/*============================================================================
 * Public Fortran function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Copy variable name from Fortran to C
 *----------------------------------------------------------------------------*/

void CS_PROCF(fclag1, FCLAG1)
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

  if (*var_id > _max_mean_vars) {

    if (_max_mean_vars == 0)
      _max_mean_vars = 16;

    while (_max_mean_vars <= *var_id)
      _max_mean_vars *= 2;

    BFT_REALLOC(_array_mean_varname, _max_mean_vars, char *);
    for (i = _last_mean_var; i < _max_mean_vars; i++)
      _array_mean_varname[i] = NULL;
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
  assert(_array_mean_varname[*var_id - 1] == NULL);

  if (l > 0) {

    /* Allocate and copy */
    BFT_MALLOC(cstr, l + 1, char);

  for (i = 0 ; i < l ; i++, i1++)
    cstr[i] = fstr[i1];

  cstr[l] = '\0';

    _array_mean_varname[*var_id - 1] = cstr;

  }

  /* Update variable counter */
  _last_mean_var = *var_id;

}

/*----------------------------------------------------------------------------
 * Copy variable name from Fortran to C
 *----------------------------------------------------------------------------*/

void CS_PROCF(fclag2, FCLAG2)
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

  if (*var_id > _max_variance_vars) {

    if (_max_variance_vars == 0)
      _max_variance_vars = 16;

    while (_max_variance_vars <= *var_id)
      _max_variance_vars *= 2;

    BFT_REALLOC(_array_variance_varname, _max_variance_vars, char *);
    for (i = _last_variance_var; i < _max_variance_vars; i++)
      _array_variance_varname[i] = NULL;
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
  assert(_array_variance_varname[*var_id - 1] == NULL);

  if (l > 0) {

    /* Allocate and copy */
    BFT_MALLOC(cstr, l + 1, char);

  for (i = 0 ; i < l ; i++, i1++)
    cstr[i] = fstr[i1];

  cstr[l] = '\0';

    _array_variance_varname[*var_id - 1] = cstr;

  }

  /* Update variable counter */
  _last_variance_var = *var_id;

}

/*----------------------------------------------------------------------------
 * Copy variable name from Fortran to C
 *----------------------------------------------------------------------------*/

void CS_PROCF(fclag3, FCLAG3)
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

  if (*var_id > _max_boundary_vars) {

    if (_max_boundary_vars == 0)
      _max_boundary_vars = 16;

    while (_max_boundary_vars <= *var_id)
      _max_boundary_vars *= 2;

    BFT_REALLOC(_array_boundary_varname, _max_boundary_vars, char *);
    for (i = _last_boundary_var; i < _max_boundary_vars; i++)
      _array_boundary_varname[i] = NULL;
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
  assert(_array_boundary_varname[*var_id - 1] == NULL);

  if (l > 0) {

    /* Allocate and copy */
    BFT_MALLOC(cstr, l + 1, char);

  for (i = 0 ; i < l ; i++, i1++)
    cstr[i] = fstr[i1];

  cstr[l] = '\0';

    _array_boundary_varname[*var_id - 1] = cstr;

  }

  /* Update variable counter */
  _last_boundary_var = *var_id;

}

/*----------------------------------------------------------------------------
 * Copy variable name from C to Fortran
 *----------------------------------------------------------------------------*/

void CS_PROCF(cfname, CFNAME)
(
 int           *const flag,    /* --> flag for array = 1, 2, or 3 */
 char          *const fstr,    /* --> Fortran string */
 int           *const len,     /* --> String Length  */
 int           *const var_id   /* --> Variable Id (1 to n) */
 CS_ARGF_SUPP_CHAINE
)
{
  int i;
  int l = 0;
  char *cstr = NULL;

  assert( *flag==1 || *flag==2 || *flag==3 );

  /* Check that variable name was set and copy string */

  switch(*flag) {
  case 1:
    if (*var_id < 1 || *var_id > _last_mean_var)
      bft_error(__FILE__, __LINE__, 0,
                _("Name of variable %d was never set.\n"), var_id);
    cstr = _array_mean_varname[*var_id - 1];
    break;
  case 2:
    if (*var_id < 1 || *var_id > _last_variance_var)
      bft_error(__FILE__, __LINE__, 0,
               _("Name of variable %d was never set.\n"), var_id);
    cstr = _array_variance_varname[*var_id - 1];
    break;
  case 3:
    if (*var_id < 1 || *var_id > _last_boundary_var)
      bft_error(__FILE__, __LINE__, 0,
                _("Name of variable %d was never set.\n"), var_id);
    cstr = _array_boundary_varname[*var_id - 1];
    break;
  }

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
 * Fortran Interface:
 *
 * SUBROUTINE UILAG1
 * *****************
 *
 * INTEGER          IILAGR     <--   type of lagrangian model used
 * INTEGER          ISUILA     <--   lagrangian restart
 * INTEGER          ISUIST     <--   lagrangian restart for statistics
 * INTEGER          NBPMAX     <--   maximum number of particles
 * INTEGER          ISTTIO     <--   stationnary calculus
 * INTEGER          INJCON     <--   continuous injection of particles
 * INTEGER          IPHYLA     <--   physical model for particles
 * INTEGER          IDPVAR     <--   equation on diameter if iphyla = 1
 * INTEGER          IMPVAR     <--   equation on mass if iphyla = 1
 * INTEGER          ITPVAR     <--   equation on temperature if iphyla = 1
 * INTEGER          IENCRA     <--   coal fouliing if iphyla = 2
 * DOUBLE           TPRENC     <--   particle coal temperature if iphyla = 2
 * DOUBLE           VISREF     <--   particle critical viscosity if iphyla = 2
 * DOUBLE           ENC1       <--   Watt and Fereday coefficient 1
 * DOUBLE           ENC2       <--   Watt and Fereday coefficient 2
 * INTEGER          NSTITS     <--   iteration number for instationnary
 * INTEGER          LTSDYN     <--   reverse coupling on dynamic
 * INTEGER          LTSMAS     <--   reverse coupling on mass
 * INTEGER          LTSTHE     <--   reverse coupling on temperature
 * INTEGER          NORDRE     <--   stochastic  differential equation order
 * INTEGER          IDISTU     <--   particle turbulent dispersion
 * INTEGER          IDIFFL     <--   particle fluid diffusion
 * INTEGER          MODCPL     <--   complete turbulent dispersion model
 * INTEGER          IDIRLA     <--   direction of the complete model
 * INTEGER          IENSI1     <--   post-processing in trajectory mode
 * INTEGER          IENSI2     <--   post-processing in movement mode
 * INTEGER          NTLAL      <--   listing printing frequency
 * INTEGER          NBVIS      <--   number of particles for display
 * INTEGER          NVISLA     <--   output period for post-processing
 * INTEGER          IVISV1     <--   display of variable 'fluid velocity'
 * INTEGER          IVISV2     <--   display of variable 'particles velocity'
 * INTEGER          IVISTP     <--   display of variable 'resident time'
 * INTEGER          IVISDM     <--   display of variable 'particle diameter'
 * INTEGER          IVISTE     <--   display of variable 'particle temperature'
 * INTEGER          IVISMP     <--   display of variable 'particle mass'
 * INTEGER          IVISHP     <--   display of variable 'coal temp. particle'
 * INTEGER          IVISDK     <--   display of variable 'core diameter of part.'
 * INTEGER          IVISCH     <--   display of variable 'mass of reactive coal'
 * INTEGER          IVISCK     <--   display of variable 'mass of char'
 * INTEGER          ISTALA     <--   calculation of volumic statistics
 * INTEGER          NBCLST     <--   number of particle clusters
 * INTEGER          SEUIL      <--   limit statistical weight value for volumic stat.
 * INTEGER          IDSTNT     <--   iteration number for volumic statistics
 * CHAR             NOMLAG     <--   mean variable name of volumic statistics
 * CHAR             NOMLAV     <--   variance variable name of volumic statistics
 * INTEGER          IHSLAG     <--   output of variable
 * INTEGER          IENSI3     <--   calculation of boundaries statistics
 * INTEGER          SEUILF     <--   limit statistical weight value for boundaries stat.
 * INTEGER          NSTBOR     <--   iteration number for boundaries statistics
 * INTEGER          INBRBD     <--   recording of particle/boundary interactions
 * INTEGER          IFLMBD     <--   recording of mass flow related to interactions
 * INTEGER          IANGBD     <--   recording of angle between particle traj./boundary
 * INTEGER          IVITBD     <--   recording of velocity of particle in an interaction
 * INTEGER          IENCBD     <--   recording of mass of coal particles
 * CHAR             NOMBRD     <--   variable name of boundaries statistics
 * INTEGER          IMOYBR     <--   cumulated value for particule/boundary interaction
 *----------------------------------------------------------------------------*/

void CS_PROCF (uilag1, UILAG1) (int *const iilagr,
                                int *const isuila,
                                int *const isuist,
                                int *const nbpmax,
                                int *const isttio,
                                int *const injcon,
                                int *const iphyla,
                                int *const idpvar,
                                int *const itpvar,
                                int *const impvar,
                                int *const iencra,
                                double tprenc[],
                                double visref[],
                                double enc1[],
                                double enc2[],
                                int *const nstits,
                                int *const ltsdyn,
                                int *const ltsmas,
                                int *const ltsthe,
                                int *const nordre,
                                int *const idistu,
                                int *const idiffl,
                                int *const modcpl,
                                int *const idirla,
                                int *const iensi1,
                                int *const iensi2,
                                int *const ntlal,
                                int *const nbvis,
                                int *const nvisla,
                                int *const ivisv1,
                                int *const ivisv2,
                                int *const ivistp,
                                int *const ivisdm,
                                int *const iviste,
                                int *const ivismp,
                                int *const ivishp,
                                int *const ivisdk,
                                int *const ivisch,
                                int *const ivisck,
                                int *const istala,
                                int *const nbclst,
                                double *const seuil,
                                int *const idstnt,
                                int ihslag[],
                                int *const iensi3,
                                double *const seuilf,
                                int *const nstbor,
                                int *const inbrbd,
                                int *const iflmbd,
                                int *const iangbd,
                                int *const ivitbd,
                                int *const iencbd,
                                int imoybr[])
{
  int i, icoal, ncoals = 0;
  int list_ind, record_ind = 1;
  char *label = NULL;
  char *attr = NULL;
  char *path1 = NULL;
  char *fmt, *opt;

  /* Global settings */

  _get_particles_model(iilagr);
  _get_status(isuila, 2, "lagrangian", "restart");
  _get_status(isttio, 2, "lagrangian", "carrier_field_stationary");
  _get_status(injcon, 2, "lagrangian", "continuous_injection");
  _get_int(nbpmax,    2, "lagrangian", "particles_max_number");

  /* Particles model */

  _get_particles_model(iphyla);

  switch (*iphyla) {
  case 1:
    _get_status(idpvar, 3, "lagrangian", "particles_models", "break_up");
    _get_status(impvar, 3, "lagrangian", "particles_models", "evaporation");
    _get_status(itpvar, 3, "lagrangian", "particles_models", "thermal");
//     if (*itpvar == 1) {
//       _get_double(tpart,  4, "lagrangian", "particles_models", "thermal", "particle_temperature");
//       _get_double(cppart, 4, "lagrangian", "particles_models", "thermal", "particle_specific_heat");
//     }
    break;
  case 2:
    _get_status(iencra, 3, "lagrangian", "particles_models", "coal_fouling");
    path1 = cs_xpath_init_path();
    cs_xpath_add_elements(&path1, 4, "lagrangian", "particles_models", "coal_fouling", "threshold_temperature");
    ncoals = cs_gui_get_nb_element(path1);
    BFT_FREE(path1);

    for (icoal=1; icoal <= ncoals; icoal++)
    {
        _get_coal_double(&tprenc[icoal-1], "threshold_temperature", icoal);
        _get_coal_double(&visref[icoal-1], "critical_viscosity",    icoal);
        _get_coal_double(&enc1[icoal-1], "fouling_coefficient_1", icoal);
        _get_coal_double(&enc2[icoal-1], "fouling_coefficient_2", icoal);
    }
    break;
  }

  /* Two-way coupling */

  if (*iilagr == 2) {
    _get_int(nstits, 3, "lagrangian", "two_way_coupling", "iteration_start");
    _get_status(ltsdyn, 3, "lagrangian", "two_way_coupling", "dynamic");
    _get_status(ltsmas, 3, "lagrangian", "two_way_coupling", "mass");
    _get_status(ltsthe, 3, "lagrangian", "two_way_coupling", "thermal");
  }

  /* Numerical modeling */

  attr = _get_attr("choice", 2, "lagrangian", "scheme_order");
  if (attr) {
    *nordre = atoi(attr);
    BFT_FREE(attr);
  }
  attr = _get_attr("choice", 2, "lagrangian", "complete_model_direction");
  if (attr) {
    *idirla = atoi(attr);
    BFT_FREE(attr);
  }
  _get_status(idistu, 2, "lagrangian", "turbulent_dispersion");
  _get_status(idiffl, 2, "lagrangian", "fluid_particles_turbulent_diffusion");
  _get_int(modcpl, 2, "lagrangian", "complete_model");

  /* Output */

  _get_status(iensi1, 3, "lagrangian", "output", "trajectory");
  _get_status(iensi2, 3, "lagrangian", "output", "particles");
  _get_status(ivisv1, 3, "lagrangian", "output", "velocity_fluid_seen");
  _get_status(ivisv2, 3, "lagrangian", "output", "velocity_particles");
  _get_status(ivistp, 3, "lagrangian", "output", "resident_time");
  _get_status(ivisdm, 3, "lagrangian", "output", "diameter");
  _get_status(iviste, 3, "lagrangian", "output", "temperature");
  _get_status(ivismp, 3, "lagrangian", "output", "mass");

  if (*iphyla == 2) {
    _get_status(ivishp, 3, "lagrangian", "output", "coal_temperature");
    _get_status(ivisdk, 3, "lagrangian", "output", "shrinking_core_diameter");
    _get_status(ivisch, 3, "lagrangian", "output", "raw_coal_mass_fraction");
    _get_status(ivisck, 3, "lagrangian", "output", "char_mass_fraction");
  }

  _get_int(nbvis,  3, "lagrangian", "output", "number_of_particles");
  _get_int(nvisla, 3, "lagrangian", "output", "postprocessing_frequency");
  _get_int(ntlal,  3, "lagrangian", "output", "listing_printing_frequency");
  fmt = _get_attr("choice", 3, "lagrangian", "output", "postprocessing_format");
  opt = _get_attr("choice", 3, "lagrangian", "output", "postprocessing_options");
  BFT_FREE(fmt);
  BFT_FREE(opt);

  /* Statistics */

  _get_int(nbclst, 3, "lagrangian", "statistics", "statistics_groups_of_particles");
  _get_status(isuist, 3, "lagrangian", "statistics", "restart");
  _get_status(istala, 3, "lagrangian", "statistics", "volume");

  if (*istala == 1) {
    _get_double(seuil, 4, "lagrangian", "statistics", "volume", "threshold_volume");
    _get_int(idstnt, 4, "lagrangian", "statistics", "volume", "iteration_start_volume");

    /* labels */

    label = _get_char_post("volume", "mean_velocity_U", &list_ind, &record_ind);
    if (label) _copy_mean_varname(label, 1);
    ihslag[1] = list_ind;

    label = _get_char_post("volume", "variance_velocity_U", &list_ind, &record_ind);
    if (label) _copy_variance_varname(label, 1);

    label = _get_char_post("volume", "mean_velocity_V", &list_ind, &record_ind);
    if (label) _copy_mean_varname(label, 2);
    ihslag[2] = list_ind;

    label = _get_char_post("volume", "variance_velocity_V", &list_ind, &record_ind);
    if (label) _copy_variance_varname(label, 2);

    label = _get_char_post("volume", "mean_velocity_W", &list_ind, &record_ind);
    if (label) _copy_mean_varname(label, 3);
    ihslag[3] = list_ind;

    label = _get_char_post("volume", "variance_velocity_W", &list_ind, &record_ind);
    if (label) _copy_variance_varname(label, 3);

    label = _get_char_post("volume", "mean_mass_fraction", &list_ind, &record_ind);
    if (label) _copy_mean_varname(label, 4);
    ihslag[4] = list_ind;

    label = _get_char_post("volume", "variance_mass_fraction", &list_ind, &record_ind);
    if (label) _copy_variance_varname(label, 4);

    label = _get_char_post("volume", "mean_resident_time", &list_ind, &record_ind);
    if (label) _copy_mean_varname(label, 5);
    ihslag[5] = list_ind;

    label = _get_char_post("volume", "variance_resident_time", &list_ind, &record_ind);
    if (label) _copy_variance_varname(label, 5);

    i = 5;

    if (*iphyla == 1) {

      if (*itpvar == 1) {
        i++;
        label = _get_char_post("volume",  "mean_temperature", &list_ind, &record_ind);
        if (label) _copy_mean_varname(label, i);
        ihslag[i] = list_ind;

        label = _get_char_post("volume", "variance_temperature", &list_ind, &record_ind);
        if (label) _copy_variance_varname(label, i);
      }

      if (*idpvar == 1) {
        i++;
        label = _get_char_post("volume", "mean_diameter", &list_ind, &record_ind);
        if (label) _copy_mean_varname(label, i);
        ihslag[i] = list_ind;

        label = _get_char_post("volume", "variance_diameter", &list_ind, &record_ind);
        if (label) _copy_variance_varname(label, i);
      }
    }

    else if (*iphyla == 2) {
      /*
      i++;
      label = _get_char_post("volume", "coal_temperature", &list_ind, &record_ind);
      if (label) _copy_mean_varname(label, i);

      label = _get_char_post("volume", "coal_temperature", &list_ind, &record_ind);
      if (label) _copy_variance_varname(label, i);
      */
      i++;
      label = _get_char_post("volume", "mean_shrinking_core_diameter", &list_ind, &record_ind);
      if (label) _copy_mean_varname(label, i);
      ihslag[i] = list_ind;

      label = _get_char_post("volume", "variance_shrinking_core_diameter", &list_ind, &record_ind);
      if (label) _copy_variance_varname(label, i);

      i++;
      label = _get_char_post("volume", "mean_raw_coal_mass_fraction", &list_ind, &record_ind);
      if (label) _copy_mean_varname(label, i);
      ihslag[i] = list_ind;

      label = _get_char_post("volume", "variance_raw_coal_mass_fraction", &list_ind, &record_ind);
      if (label) _copy_variance_varname(label, i);

      i++;
      label = _get_char_post("volume", "mean_char_mass_fraction", &list_ind, &record_ind);
      if (label) _copy_mean_varname(label, i);
      ihslag[i] = list_ind;

      label = _get_char_post("volume", "variance_char_mass_fraction", &list_ind, &record_ind);
      if (label) _copy_variance_varname(label, i);
    }

    i++;
    label = _get_char_post("volume", "statistical_weight", &list_ind, &record_ind);
    if (label) _copy_mean_varname(label, i);
    ihslag[i] = list_ind;
  }

  _get_status(iensi3, 3, "lagrangian", "statistics", "boundary");

  if (*iensi3 == 1) {

    _get_double(seuilf, 4, "lagrangian", "statistics", "boundary", "threshold_boundary");
    _get_int(nstbor, 4, "lagrangian", "statistics", "boundary", "iteration_start_boundary");

    i = 0;

    label = _get_char_post("boundary", "impacts", inbrbd, &record_ind);
    if (*inbrbd) {
      i++;
      if (label) _copy_boundary_varname(label, i);
      imoybr[i] = record_ind;
    }

    label = _get_char_post("boundary", "mass_flux", iflmbd, &record_ind);
    if (*iflmbd) {
      i++;
      if (label) _copy_boundary_varname(label, i);
      imoybr[i] = record_ind;
    }

    label = _get_char_post("boundary", "angle", iangbd, &record_ind);
    if (*iangbd) {
      i++;
      if (label) _copy_boundary_varname(label, i);
      imoybr[i] = record_ind;
    }

    label = _get_char_post("boundary", "velocity", ivitbd, &record_ind);
    if (*ivitbd) {
      i++;
      if (label) _copy_boundary_varname(label, i);
      imoybr[i] = record_ind;
    }
    label = _get_char_post("boundary", "coal_fouling", iencbd, &record_ind);
    if (*iencbd) {
      i++;
      if (label) _copy_boundary_varname(label, i);
      imoybr[i] = record_ind;
    }
  }
  BFT_FREE(label);

#if _XML_DEBUG_
  printf("==>UILAG1\n");
  printf("--iilagr = %i\n", *iilagr);
  printf("--isuila = %i\n", *isuila);
  printf("--isttio = %i\n", *isttio);
  printf("--nbpmax = %i\n", *nbpmax);
  printf("--isttio = %i\n", *isttio);
  printf("--injcon = %i\n", *injcon);
  printf("--iphyla = %i\n", *iphyla);
  switch(*iphyla) {
  case 0:
    break;
  case 1:
    printf("--idpvar = %i\n", *idpvar);
    printf("--impvar = %i\n", *impvar);
    printf("--itpvar = %i\n", *itpvar);
    break;
  case 2:
    printf("--iencra = %i\n", *iencra);
    for (icoal=1; icoal<=ncoals; icoal++)
      {
        printf("--tprenc[%i] = %f\n", icoal, tprenc[icoal-1]);
        printf("--visref[%i] = %f\n", icoal, visref[icoal-1]);
        printf("--enc1[%i] = %f\n", icoal, enc1[icoal-1]);
        printf("--enc2[%i] = %f\n", icoal, enc2[icoal-1]);
      }
    break;
  }

  if (*iilagr == 2) {
    printf("--nstits = %i\n", *nstits);
    printf("--ltsdyn = %i\n", *ltsdyn);
    printf("--ltsmas = %i\n", *ltsmas);
    printf("--ltsthe = %i\n", *ltsthe);
  }

  printf("--nordre = %i\n", *nordre);
  printf("--idistu = %i\n", *idistu);
  printf("--idiffl = %i\n", *idiffl);
  printf("--modcpl = %i\n", *modcpl);
  printf("--idirla = %i\n", *idirla);

  printf("--iensi1 = %i\n", *iensi1);
  printf("--iensi2 = %i\n", *iensi2);
  printf("--ivisv1 = %i\n", *ivisv1);
  printf("--ivisv2 = %i\n", *ivisv2);
  printf("--ivistp = %i\n", *ivistp);
  printf("--ivisdm = %i\n", *ivisdm);
  printf("--iviste = %i\n", *iviste);
  printf("--ivismp = %i\n", *ivismp);

  if (*iphyla == 2) {
    printf("--ivishp = %i\n", *ivishp);
    printf("--ivisdk = %i\n", *ivisdk);
    printf("--ivisch = %i\n", *ivisch);
    printf("--ivisck = %i\n", *ivisck);
  }

  printf("--nbvis  = %i\n", *nbvis);
  printf("--nvisla = %i\n", *nvisla);

  printf("--isuist = %i\n", *isuist);
  printf("--nbclst = %i\n", *nbclst);

  printf("--istala = %i\n", *istala);
  if (*istala == 1) {
    printf("--idstnt = %i\n", *idstnt);
    printf("--seuil  = %f\n", *seuil);

    /*
    printf("--i        nomlag             nomlav              ihslag\n");
    for (i=1; i <= 5; i++)
      printf("  %i %30s %30s %5i\n", i, nomlag[i], nomlav[i], ihslag[i]);
    i = 5;
    if (*iphyla == 1) {
      if (*itpvar == 1) {
        i++;
        printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      }
      if (*idpvar == 1) {
        i++;
        printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      }
    }
    else if (*iphyla == 2) {
      //i++;
      //printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      i++;
      printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      i++;
      printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      i++;
      printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
    }
    i++;
    printf("  %i %s \n", i, nomlag[i]);
    */
  }

  printf("--iensi3 = %i\n", *iensi3);
  if (*iensi3 == 1) {
    printf("--nstbor = %i\n", *nstbor);
    printf("--seuilf = %f\n", *seuilf);
    printf("--inbrbd = %i\n", *inbrbd);
    printf("--iflmbd = %i\n", *iflmbd);
    printf("--iangbd = %i\n", *iangbd);
    printf("--ivitbd = %i\n", *ivitbd);
    printf("--iencbd = %i\n", *iencbd);
  }

#endif

}

/*-----------------------------------------------------------------------------
 * Fortran Interface:
 *
 * subroutine uilag2
 * *****************
 *
 * integer          iphyla     -->   physica model associated to the particles
 * integer          nflagm     -->   max number of boundaries
 * integer          iusncl     <--   array for particles class(es) number
 * integer          iusclb     <--   array for particles boundary conditions
 * integer          iuslag     <--   array for integer variables
 * double precision ruslag     <--   array for real variables
 *----------------------------------------------------------------------------*/

void CS_PROCF (uilag2, UILAG2) (const int *const iphyla,
                                const int *const nclagm,
                                const int *const nflagm,
                                int     iusncl[],
                                int     iusclb[],
                                int     iuslag[],
                                double  ruslag[])
{
  int izone, nzones;
  int iclas, nclasses;
  char *interaction = NULL;
  char szone[15], sclass[10];
  char *path, *path1, *path2;
  char *choice, *nature, *zlabel = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "boundary");
  nzones = cs_gui_get_nb_element(path);
  path = NULL;

  for (izone=1; izone<=nzones; izone++)
  {
      path = cs_xpath_init_path();
      sprintf(szone, "boundary[%i]", izone);
      cs_xpath_add_elements(&path, 2, "boundary_conditions", szone);
      zlabel = _get_attr("label",  2, "boundary_conditions", szone);
      nature = _get_attr("nature", 2, "boundary_conditions", szone);

      path2 = cs_xpath_init_path();
      cs_xpath_add_elements(&path2, 2, "boundary_conditions", nature);
      cs_xpath_add_test_attribute(&path2, "label", zlabel);
      cs_xpath_add_element(&path2, "particles");

      BFT_MALLOC(path1, strlen(path2)+1, char);
      strcpy(path1, path2);
      cs_xpath_add_attribute(&path1, "choice");
      interaction = cs_gui_get_attribute_value(path1);
      nclasses = 0;

      if (interaction != NULL) {

        if (cs_gui_strcmp(interaction, "inlet"))
        {
            iusclb[izone] = 1;
            strcpy(path1, path2);
            cs_xpath_add_element(&path1, "class");
            nclasses = cs_gui_get_nb_element(path1);
            iusncl[izone] = nclasses;
            strcpy(path1, path2);

            for (iclas=1; iclas<=nclasses; iclas++)
            {
                sprintf(sclass, "class[%i]", iclas);
                BFT_REALLOC(path2, 20+strlen(nature)+10+strlen(zlabel)+13+strlen(sclass)+1, char);
                strcpy(path2, "");
                sprintf(path2, "boundary_conditions/%s[@label='%s']/particles/%s", nature, zlabel, sclass);

                // Arrays IUSLAG and RUSLAG defined in lagran.h as
                //
                // IUSLAG(NCLAGM, NFLAGM, NDLAIM) with NCLAGM = 20, NFLAGM = 100, NDLAIM = 10
                // RUSLAG(NCLAGM, NFLAGM, NDLAGM) with NCLAGM = 20, NFLAGM = 100, NDLAGM = 50
                //
                // F77: IUSLAG(ICLAS, IZONE, I)
                // C  : &iuslag[i][izone][iclas] = &iuslag[ i*(*nclagm)*(*nflagm) + izone*(*nclagm) + iclas ]

                _get_int(&iuslag[1*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "number");
                _get_int(&iuslag[2*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "frequency");
                _get_int(&iuslag[3*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "statistical_groups");

                choice = _get_attr("choice", 2, path2, "velocity");
                if (cs_gui_strcmp(choice, "fluid"))
                  iuslag[4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas] = -1;
                else if (cs_gui_strcmp(choice, "norm")) {
                  iuslag[4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas] = 0;
                  _get_double(&ruslag[1*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 3, path2, "velocity", "norm");
                }
                else if (cs_gui_strcmp(choice, "components")) {
                  iuslag[4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas] = 1;
                   _get_double(&ruslag[2*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 3, path2, "velocity", "velocity_x");
                   _get_double(&ruslag[3*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 3, path2, "velocity", "velocity_y");
                   _get_double(&ruslag[4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 3, path2, "velocity", "velocity_z");
                }
                else if (cs_gui_strcmp(choice, "subroutine"))
                  iuslag[4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas] = 2;

                choice = _get_attr("choice", 2, path2, "temperature");
                if (cs_gui_strcmp(choice, "prescribed")) {
                  iuslag[5*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas] = 1;
                  _get_double(&ruslag[5*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "temperature");
                }
                else if (cs_gui_strcmp(choice, "subroutine"))
                  iuslag[5*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas] = 2;

                choice = _get_attr("choice", 2, path2, "diameter");
                if (cs_gui_strcmp(choice, "prescribed")) {
                  iuslag[6*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas] = 1;
                  _get_double(&ruslag[6*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "diameter");
                  _get_double(&ruslag[7*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "diameter_standard_deviation");
                }
                else if (cs_gui_strcmp(choice, "subroutine"))
                  iuslag[6*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas] = 2;

                _get_double(&ruslag[8*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "density");

                if (*iphyla == 1) {
                  _get_double(&ruslag[9*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "specific_heat");
                }

                _get_double(&ruslag[10*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "emissivity");
                _get_double(&ruslag[11*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "statitical_weight");
                _get_double(&ruslag[12*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "mass_flow_rate");

                if (*iphyla == 2) {
                  _get_int(&iuslag[7*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "coal_number");
                  _get_double(&ruslag[13*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "coal_temperature");
                  _get_double(&ruslag[14*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "raw_coal_mass_fraction");
                  _get_double(&ruslag[15*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas], 2, path2, "char_mass_fraction");
                }

              }

          }

        else if(cs_gui_strcmp(interaction, "outlet"))
          iusclb[izone] = 2;

        else if(cs_gui_strcmp(interaction, "bounce"))
          iusclb[izone] = 3;

        else if(cs_gui_strcmp(interaction, "deposit1"))
          iusclb[izone] = 4;

        else if(cs_gui_strcmp(interaction, "deposit2"))
          iusclb[izone] = 5;

        else if(cs_gui_strcmp(interaction, "deposit3"))
          iusclb[izone] = 6;

        else if(cs_gui_strcmp(interaction, "fouling"))
          iusclb[izone] = 13;

      }

      BFT_FREE(path1);
      BFT_FREE(path2);

    }
  BFT_FREE(path);

#if _XML_DEBUG_
  printf("==>UILAG2\n");
  for (izone=1; izone<=nzones; izone++)
    {
      printf("--iusclb[%i] = %i has %i class(es) \n", izone, iusclb[izone], iusncl[izone]);

      if ( (iusclb[izone]==1) && (iusncl[izone]!=0) ) {
        nclasses = iusncl[izone];
        for (iclas=1; iclas<=nclasses; iclas++)
          {
            printf("--iusncl[%i] : class number %i \n", iclas-1, iclas);

            printf("--iuslag[%i] = %i \n", 1*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, iuslag[1*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            printf("--iuslag[%i] = %i \n", 2*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, iuslag[2*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            printf("--iuslag[%i] = %i \n", 3*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, iuslag[3*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            printf("--iuslag[%i] = %i velocity (-1:fluid, 0:norm, 1:components)\n", 4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, iuslag[4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            if (iuslag[4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]==0)
              printf("--ruslag[%i] = %f \n", 1*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[1*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            else if (iuslag[4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]==1) {
              printf("--ruslag[%i] = %f \n", 2*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[2*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
              printf("--ruslag[%i] = %f \n", 3*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[3*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
              printf("--ruslag[%i] = %f \n", 4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[4*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            }
            printf("--iuslag[%i] = %i temperature(1:prescribed, 2:subroutine)\n", 5*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, iuslag[5*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            if (iuslag[5*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]==1)
              printf("--ruslag[%i] = %f \n", 5*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[5*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);

            printf("--iuslag[%i] = %i diameter(1:prescribed, 2:subroutine)\n", 6*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, iuslag[6*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            if (iuslag[6*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]==1) {
              printf("--ruslag[%i] = %f \n", 6*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[6*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
              printf("--ruslag[%i] = %f \n", 7*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[7*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            }
            printf("--ruslag[%i] = %f \n", 8*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[8*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            if (*iphyla == 1)
                printf("--ruslag[%i] = %f \n", 9*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[9*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            printf("--ruslag[%i] = %f \n", 10*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[10*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            printf("--ruslag[%i] = %f \n", 11*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[11*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            printf("--ruslag[%i] = %f \n", 12*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[12*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            if (*iphyla == 2) {
              printf("--iuslag[%i] = %i \n", 7*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, iuslag[7*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
              printf("--ruslag[%i] = %f \n", 13*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[13*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
              printf("--ruslag[%i] = %f \n", 14*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[14*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
              printf("--ruslag[%i] = %f \n", 15*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas, ruslag[15*(*nclagm)*(*nflagm)+izone*(*nclagm)+iclas]);
            }
          }
      }
    }
#endif
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
