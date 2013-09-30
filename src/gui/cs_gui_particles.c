/*============================================================================
 * Management of the GUI parameters file: particles tracking
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_gui.h"
#include "cs_gui_util.h"
#include "cs_gui_boundary_conditions.h"
#include "cs_prototypes.h"

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

/*
  rcodcl[ k * dim1 *dim2 + j *dim1 + i]
*/

/*============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structures associated to lagrangian particles definition
 *----------------------------------------------------------------------------*/

typedef struct {
  char       **label;             /* label for each boundary zone                    */
  char       **nature;            /* nature for each boundary zone                   */
  char       **p_nature;          /* specific nature of the boundary for particles   */
  int         *n_classes;         /* number of classes for each zone                 */
  int        **n_particles;       /* number of particles for each class              */
  int        **frequency;         /* frequency of injection for each class           */
  int        **statistical_groups;/* frequency of injection for each class           */
  double     **statistical_weight;/* number of real particles for numerical particles*/
  double     **mass_flow_rate;    /* mass flow rate of particles                     */
  double     **density;           /* density for each class                          */
  double     **diameter;          /* diameter for each class                         */
  double     **standard_deviation;/* standard deviation of diameter for each class   */
  double     **specific_heat;     /* specific heat for each class                    */
  double     **emissivity;        /* emissivity for each class                       */
#if 0
  mei_tree_t **velocity;          /* formula for norm or mass flow rate of velocity  */
  mei_tree_t **direction;         /* formula for direction of velocity               */
#endif
} cs_particles_boundary_t;

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
 * Static global variables
 *============================================================================*/

/* Pointer on the main boundaries structure */

extern cs_boundary_t *boundaries;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return value of the particles model
 *----------------------------------------------------------------------------*/

static void
_get_particles_model(const char *const model, int *const imodel)
{
  char *path;
  char *attr;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "lagrangian", model);
  cs_xpath_add_attribute(&path, "model");
  attr = cs_gui_get_attribute_value(path);

  if (attr != NULL) {
    if (cs_gui_strcmp(attr, "off"))
      *imodel = 0;
    else if (cs_gui_strcmp(attr, "one_way"))
      *imodel = 1;
    else if (cs_gui_strcmp(attr, "two_way"))
      *imodel = 2;
    else if (cs_gui_strcmp(attr, "frozen"))
      *imodel = 3;
    else if (cs_gui_strcmp(attr, "thermal"))
      *imodel = 1;
    else if (cs_gui_strcmp(attr, "coal"))
      *imodel = 2;
    BFT_FREE(attr);
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
               int  *record_value)
{
  char *path, *path1, *path2 = NULL;
  char *label = NULL;
  int result;

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

    cs_xpath_add_element(&path1, "postprocessing_recording");
    cs_xpath_add_attribute(&path1, "status");
    if (cs_gui_get_status(path1, &result))
      *record_value = result;
  }

  else if (cs_gui_strcmp(type, "boundary")) {

    cs_xpath_add_element(&path2, "postprocessing_recording");
    cs_xpath_add_attribute(&path2, "status");
    if (cs_gui_get_status(path2, &result))
      *record_value = result;
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
_copy_mean_varname(const char *varname, int ipp)
{
  size_t  l;
  assert(ipp > 0);

  if (ipp < 1 || ipp > _last_mean_var+1)
    bft_error(__FILE__, __LINE__, 0,
              _("Variable index %i out of bounds (1 to %i)"),
              ipp, _last_mean_var);

  l = strlen(varname);

  if (_array_mean_varname[ipp-1] == NULL)
    BFT_MALLOC(_array_mean_varname[ipp-1], l + 1, char);

  else if (strlen(_array_mean_varname[ipp-1]) != l)
    BFT_REALLOC(_array_mean_varname[ipp-1], l + 1, char);

  strcpy(_array_mean_varname[ipp-1], varname);
}

/*-----------------------------------------------------------------------------
 * Copy a variable name to the variance variable names array
 *
 * parameters:
 *   varname        -->  name or label of the variable/scalar/property
 *   ipp            -->  index from the fortran array associated to varname
 *----------------------------------------------------------------------------*/

static void
_copy_variance_varname(const char *varname, int ipp)
{
  size_t  l;
  assert(ipp > 0);

  if (ipp < 1 || ipp > _last_variance_var+1)
    bft_error(__FILE__, __LINE__, 0,
              _("Variable index %i out of bounds (1 to %i)"),
              ipp, _last_variance_var);

  l = strlen(varname);

  if (_array_variance_varname[ipp-1] == NULL)
    BFT_MALLOC(_array_variance_varname[ipp-1], l + 1, char);

  else if (strlen(_array_variance_varname[ipp-1]) != l)
    BFT_REALLOC(_array_variance_varname[ipp-1], l + 1, char);

  strcpy(_array_variance_varname[ipp-1], varname);
}

/*-----------------------------------------------------------------------------
 * Copy a variable name to the variance variable names array
 *
 * parameters:
 *   varname        -->  name or label of the variable/scalar/property
 *   ipp            -->  index from the fortran array associated to varname
 *----------------------------------------------------------------------------*/

static void
_copy_boundary_varname(const char *varname, int ipp)
{
  size_t  l;
  assert(ipp > 0);

  if (ipp < 1 || ipp > _last_boundary_var+1)
    bft_error(__FILE__, __LINE__, 0,
              _("Variable index %i out of bounds (1 to %i)"),
              ipp, _last_boundary_var);

  l = strlen(varname);

  if (_array_boundary_varname[ipp-1] == NULL)
    BFT_MALLOC(_array_boundary_varname[ipp-1], l + 1, char);

  else if (strlen(_array_boundary_varname[ipp-1]) != l)
    BFT_REALLOC(_array_boundary_varname[ipp-1], l + 1, char);

  strcpy(_array_boundary_varname[ipp-1], varname);
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
                _("Name of variable %i was never set.\n"), *var_id);
    cstr = _array_mean_varname[*var_id - 1];
    break;
  case 2:
    if (*var_id < 1 || *var_id > _last_variance_var)
      bft_error(__FILE__, __LINE__, 0,
               _("Name of variable %i was never set.\n"), *var_id);
    cstr = _array_variance_varname[*var_id - 1];
    break;
  case 3:
    if (*var_id < 1 || *var_id > _last_boundary_var)
      bft_error(__FILE__, __LINE__, 0,
                _("Name of variable %i was never set.\n"), *var_id);
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
 * INTEGER          NLAYER     -->   max number of layer per coal particle
 * INTEGER          IILAGR     <--   type of lagrangian model used
 * INTEGER          ISUILA     <--   lagrangian restart
 * INTEGER          ISUIST     <--   lagrangian restart for statistics
 * INTEGER          NBPMAX     <--   maximum number of particles
 * INTEGER          ISTTIO     <--   stationnary calculus
 * INTEGER          INJCON     <--   continuous injection of particles
 * INTEGER          IDEPST     <--   particle deposition submodel
 * INTEGER          IPHYLA     <--   physical model for particles
 * INTEGER          IDPVAR     <--   equation on diameter if iphyla = 1
 * INTEGER          IMPVAR     <--   equation on mass if iphyla = 1
 * INTEGER          ITPVAR     <--   equation on temperature if iphyla = 1
 * INTEGER          IENCRA     <--   coal fouling if iphyla = 2
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
 * INTEGER          NTLAL      <--   listing printing frequency
 * INTEGER          IVISV1     <--   display of variable 'fluid velocity'
 * INTEGER          IVISV2     <--   display of variable 'particles velocity'
 * INTEGER          IVISTP     <--   display of variable 'resident time'
 * INTEGER          IVISDM     <--   display of variable 'particle diameter'
 * INTEGER          IVISTE     <--   display of variable 'particle temperature'
 * INTEGER          IVISMP     <--   display of variable 'particle mass'
 * INTEGER          IVISDK     <--   display of variable 'core diameter of part.'
 * INTEGER          IVISWAT    <--   display of variable 'mass of moisture'
 * INTEGER          IVISCH     <--   display of variable 'mass of reactive coal'
 * INTEGER          IVISCK     <--   display of variable 'mass of char'
 * INTEGER          ISTALA     <--   calculation of volumic statistics
 * INTEGER          NBCLST     <--   number of particle clusters
 * INTEGER          SEUIL      <--   limit statistical weight value for volumic stat.
 * INTEGER          IDSTNT     <--   iteration number for volumic statistics
 * INTEGER          NSTIST     <--   iteration number for steady-state volumic statistics
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
 * INTEGER          IENCNBBD   <--   recording of particle/boundary interactions with fouling
 * INTEGER          IENCMABD   <--   recording of flux mass of coal particles (fouling)
 * INTEGER          IENCDIBD   <--   recording of coal particles diameter (fouling)
 * INTEGER          IENCCKBD   <--   recording of coal particles coke fraction (fouling)
 * CHAR             NOMBRD     <--   variable name of boundaries statistics
 * INTEGER          IMOYBR     <--   cumulated value for particule/boundary interaction
 *----------------------------------------------------------------------------*/

void CS_PROCF (uilag1, UILAG1) (int *const nlayer,
                                int *const iilagr,
                                int *const isuila,
                                int *const isuist,
                                int *const nbpmax,
                                int *const isttio,
                                int *const injcon,
                                int *const idepst,
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
                                int *const ntlal,
                                int *const ivisv1,
                                int *const ivisv2,
                                int *const ivistp,
                                int *const ivisdm,
                                int *const iviste,
                                int *const ivismp,
                                int *const ivisdk,
                                int *const iviswat,
                                int *const ivisch,
                                int *const ivisck,
                                int *const istala,
                                int *const nbclst,
                                double *const seuil,
                                int *const idstnt,
                                int *const nstist,
                                int ihslag[],
                                int *const iensi3,
                                double *const seuilf,
                                int *const nstbor,
                                int *const inbrbd,
                                int *const iflmbd,
                                int *const iangbd,
                                int *const ivitbd,
                                int *const iencnbbd,
                                int *const iencmabd,
                                int *const iencdibd,
                                int *const iencckbd,
                                int imoybr[],
                                int *const iactfv,
                                int *const iactvx,
                                int *const iactvy,
                                int *const iactvz,
                                int *const iactts)
{
  int i, icoal, ilayer, ncoals = 0;
  int list_ind = 1;
  int record_ind = 1;
  char *name = NULL;
  char *attr = NULL;
  char *path1 = NULL;

  attr = _get_attr("model", 1, "lagrangian");
  if (attr == NULL || cs_gui_strcmp(attr, "off")) {
    *iilagr = 0;
    BFT_FREE(attr);
#if _XML_DEBUG_
    bft_printf("==>UILAG1\n");
    bft_printf("--iilagr = %i\n", *iilagr);
#endif
    return;
  }
  BFT_FREE(attr);

  /* Global settings */

  _get_particles_model("coupling_mode", iilagr);
  _get_status(isuila, 2, "lagrangian", "restart");
  _get_status(isttio, 2, "lagrangian", "carrier_field_stationary");
  _get_status(injcon, 2, "lagrangian", "continuous_injection");
  _get_status(idepst, 2, "lagrangian", "deposition_submodel");

  bft_printf("idepst = %d",*idepst);

  _get_int(nbpmax,    2, "lagrangian", "particles_max_number");

  /* Particles model */

  _get_particles_model("particles_models", iphyla);

  switch (*iphyla) {
  case 1:
    _get_status(idpvar, 3, "lagrangian", "particles_models", "break_up");
    _get_status(impvar, 3, "lagrangian", "particles_models", "evaporation");
    _get_status(itpvar, 3, "lagrangian", "particles_models", "thermal");
    /*
    if (*itpvar == 1) {
      _get_double(tpart,  4, "lagrangian", "particles_models", "thermal", "particle_temperature");
      _get_double(cppart, 4, "lagrangian", "particles_models", "thermal", "particle_specific_heat");
    }
    */
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

  _get_status(ivisv1, 3, "lagrangian", "output", "velocity_fluid_seen");
  _get_status(ivisv2, 3, "lagrangian", "output", "velocity_particles");
  _get_status(ivistp, 3, "lagrangian", "output", "resident_time");
  _get_status(ivisdm, 3, "lagrangian", "output", "diameter");
  _get_status(iviste, 3, "lagrangian", "output", "temperature");
  _get_status(ivismp, 3, "lagrangian", "output", "mass");

  if (*iphyla == 2) {
    _get_status(ivisdk,  3, "lagrangian", "output", "shrinking_core_diameter");
    _get_status(iviswat, 3, "lagrangian", "output", "moisture_mass_fraction");
    _get_status(ivisch,  3, "lagrangian", "output", "raw_coal_mass_fraction");
    _get_status(ivisck,  3, "lagrangian", "output", "char_mass_fraction");
  }

  _get_int(ntlal,  3, "lagrangian", "output", "listing_printing_frequency");

  /* Statistics */

  _get_int(nbclst, 3, "lagrangian", "statistics", "statistics_groups_of_particles");
  _get_status(isuist, 3, "lagrangian", "statistics", "restart");
  _get_status(istala, 3, "lagrangian", "statistics", "volume");

  if (*istala == 1) {
    _get_double(seuil, 4, "lagrangian", "statistics", "volume", "threshold_volume");
    _get_int(idstnt, 4, "lagrangian", "statistics", "volume", "iteration_start_volume");
    _get_int(nstist, 4, "lagrangian", "statistics", "volume", "iteration_steady_start_volume");

    /* labels */

    i  = 0;

    _get_char_post("volume", "Part_vol_frac", iactfv);
    if (*iactfv)
    {
      i++;
      _copy_mean_varname("Part_vol_frac", i);
      _copy_variance_varname("var_Part_vol_frac", i);
      ihslag[i-1] = 1;
    }

    _get_char_post("volume", "Part_velocity_X", iactvx);
    if (*iactvx)
    {
      i++;
      _copy_mean_varname("Part_velocity_X", i);
      _copy_variance_varname("var_Part_velocity_X", i);
      ihslag[i-1] = 1;
    }

    _get_char_post("volume", "Part_velocity_Y", iactvy);
    if (*iactvy)
    {
      i++;
      _copy_mean_varname("Part_velocity_Y", i);
      _copy_variance_varname("var_Part_velocity_X", i);
      ihslag[i-1] = 1;
    }

    _get_char_post("volume", "Part_velocity_Z", iactvz);
    if (*iactvz)
    {
      i++;
      _copy_mean_varname("Part_velocity_Z", i);
      _copy_variance_varname("var_Part_velocity_Z", i);
      ihslag[i-1] = 1;
    }

    _get_char_post("volume", "Part_resid_time", iactts);
    if (*iactts)
    {
      i++;
      _copy_mean_varname("Part_resid_time", i);
      _copy_variance_varname("var_Part_resid_time", i);
      ihslag[i-1] = 1;
    }

    if (*iphyla == 1) {

      if (*itpvar == 1) {
        i++;
        _copy_mean_varname("Part_temperature", i);
        _copy_variance_varname("var_Part_temperature", i);
        ihslag[i-1] = 1;
      }

      if (*idpvar == 1) {
        i++;
        _copy_mean_varname("Part_diameter", i);
        _copy_variance_varname("var_Part_diameter", i);
        ihslag[i-1] = 1;
      }
    }

    else if (*iphyla == 2) {

      /* Particle Mass */
      i++;
      _copy_mean_varname("Part_mass", i);
      _copy_variance_varname("var_Part_mass", i);
      ihslag[i-1] = 1;

      /* Particle Temperature */
      if (*nlayer < 2) {
            i++;
        _copy_mean_varname("Part_temperature", i);
        _copy_variance_varname("var_Part_temperature", i);
        ihslag[i-1] = 1;
      }
      else {
        BFT_MALLOC(name, strlen("Part_temperature_layer_")+1 + 2, char);
        for (ilayer = 0; ilayer < *nlayer; ilayer++){
          i++;
          ihslag[i-1] = 1;
          /* Mean name */
          BFT_REALLOC(name, strlen("Part_temperature_layer_")+1 + 2, char);
          sprintf(name,
                  "%s%2.2i",
                  "Part_temperature_layer_",
                  ilayer+1);
          _copy_mean_varname(name, i);
          /* Variance name */
          BFT_REALLOC(name, strlen("var_Part_temperature_layer_")+1 + 2, char);
          sprintf(name,
                  "%s%2.2i",
                  "var_Part_temperature_layer_",
                  ilayer+1);
          _copy_variance_varname(name, i);
        }
        BFT_FREE(name);
      }

      /* Particle Water Mass */
      i++;
      _copy_mean_varname("Part_wat_mass", i);
      _copy_variance_varname("var_Part_wat_mass", i);
      ihslag[i-1] = 1;

      /* Particle Reactive Coal Mass */
      if (*nlayer < 2) {
        i++;
        _copy_mean_varname("Part_ch_mass", i);
        _copy_variance_varname("var_Part_ch_mass", i);
        ihslag[i-1] = 1;
      }
      else {
        BFT_MALLOC(name, strlen("Part_ch_mass_layer_")+1 + 2, char);
        for (ilayer = 0; ilayer < *nlayer; ilayer++){
          i++;
          ihslag[i-1] = 1;
          /* Mean name */
          BFT_REALLOC(name, strlen("Part_ch_mass_layer_")+1 + 2, char);
          sprintf(name,
                  "%s%2.2i",
                  "Part_ch_mass_layer_",
                  ilayer+1);
          _copy_mean_varname(name, i);
          /* Variance name */
          BFT_REALLOC(name, strlen("var_Part_ch_mass_layer_")+1 + 2, char);
          sprintf(name,
                  "%s%2.2i",
                  "var_Part_ch_mass_layer_",
                  ilayer+1);
          _copy_variance_varname(name, i);
        }
        BFT_FREE(name);
      }

      /* Particle Coke Mass */
      if (*nlayer < 2) {
        i++;
        _copy_mean_varname("Part_ck_mass", i);
        _copy_variance_varname("var_Part_ck_mass", i);
        ihslag[i-1] = 1;
      }
      else {
        BFT_MALLOC(name, strlen("Part_ck_mass_layer_")+1 + 2, char);
        for (ilayer = 0; ilayer < *nlayer; ilayer++){
          i++;
          ihslag[i-1] = 1;
          /* Mean name */
          BFT_REALLOC(name, strlen("Part_ck_mass_layer_")+1 + 2, char);
          sprintf(name,
                  "%s%2.2i",
                  "Part_ck_mass_layer_",
                  ilayer+1);
          _copy_mean_varname(name, i);
          /* Variance name */
          BFT_REALLOC(name, strlen("var_Part_ck_mass_layer_")+1 + 2, char);
          sprintf(name,
                  "%s%2.2i",
                  "var_Part_ck_mass_layer_",
                  ilayer+1);
          _copy_variance_varname(name, i);
        }
        BFT_FREE(name);
      }

      /* Particle shrinking core diameter */
      i++;
      _copy_mean_varname("Part_shrink_core_diam", i);
      _copy_variance_varname("var_Part_shrink_core_diam", i);
      ihslag[i] = 1;
    }

    i++;
    _get_char_post("volume", "Part_statis_weight",  &record_ind);
    _copy_mean_varname("Part_statis_weight", i);
    ihslag[i-1] = 1;
  }

  _get_status(iensi3, 3, "lagrangian", "statistics", "boundary");

  if (*iensi3 == 1) {

    _get_double(seuilf, 4, "lagrangian", "statistics", "boundary", "threshold_boundary");
    _get_int(nstbor, 4, "lagrangian", "statistics", "boundary", "iteration_start_boundary");

    i = 0;

    _get_char_post("boundary", "Part_impact_number", inbrbd);
    if (*inbrbd) {
      imoybr[i] = 0;
      i++;
      _copy_boundary_varname("Part_impact_number", i);
    }

    _get_char_post("boundary", "Part_bndy_mass_flux", iflmbd);
    if (*iflmbd) {
      imoybr[i] = 1;
      i++;
      _copy_boundary_varname("Part_bndy_mass_flux", i);
    }

    _get_char_post("boundary", "Part_impact_angle", iangbd);
    if (*iangbd) {
      imoybr[i] = 2;
      i++;
      _copy_boundary_varname("Part_impact_angle", i);
    }

    _get_char_post("boundary", "Part_impact_velocity", ivitbd);
    if (*ivitbd) {
      imoybr[i] = 2;
      i++;
      _copy_boundary_varname("Part_impact_velocity", i);
    }

    /* Coal fouling statistics*/
    _get_char_post("boundary", "Part_fouled_impact_number", iencnbbd);
    if (*iencnbbd) {
      imoybr[i] = 0;
      i++;
      _copy_boundary_varname("Part_fouled_impact_number", i);
    }
    _get_char_post("boundary", "Part_fouled_mass_flux", iencmabd);
    if (*iencmabd) {
      imoybr[i] = 1;
      i++;
      _copy_boundary_varname("Part_fouled_mass_flux", i);
    }
    _get_char_post("boundary", "Part_fouled_diam", iencdibd);
    if (*iencdibd) {
      imoybr[i] = 3;
      i++;
      _copy_boundary_varname("Part_fouled_diam", i);
    }
    _get_char_post("boundary", "Part_fouled_Xck", iencckbd);
    if (*iencckbd) {
      imoybr[i] = 3;
      i++;
      _copy_boundary_varname("Part_fouled_Xck", i);
    }
  }

#if _XML_DEBUG_
  bft_printf("==>UILAG1\n");
  bft_printf("--iilagr = %i\n", *iilagr);
  bft_printf("--isuila = %i\n", *isuila);
  bft_printf("--isttio = %i\n", *isttio);
  bft_printf("--nbpmax = %i\n", *nbpmax);
  bft_printf("--isttio = %i\n", *isttio);
  bft_printf("--injcon = %i\n", *injcon);
  bft_printf("--idepst = %i\n", *idepst);
  bft_printf("--iphyla = %i\n", *iphyla);
  switch(*iphyla) {
  case 0:
    break;
  case 1:
    bft_printf("--idpvar = %i\n", *idpvar);
    bft_printf("--impvar = %i\n", *impvar);
    bft_printf("--itpvar = %i\n", *itpvar);
    break;
  case 2:
    bft_printf("--iencra = %i\n", *iencra);
    for (icoal=1; icoal<=ncoals; icoal++)
    {
      bft_printf("--tprenc[%i] = %f\n", icoal, tprenc[icoal-1]);
      bft_printf("--visref[%i] = %f\n", icoal, visref[icoal-1]);
      bft_printf("--enc1[%i] = %f\n", icoal, enc1[icoal-1]);
      bft_printf("--enc2[%i] = %f\n", icoal, enc2[icoal-1]);
    }
    break;
  }

  if (*iilagr == 2) {
    bft_printf("--nstits = %i\n", *nstits);
    bft_printf("--ltsdyn = %i\n", *ltsdyn);
    bft_printf("--ltsmas = %i\n", *ltsmas);
    bft_printf("--ltsthe = %i\n", *ltsthe);
  }

  bft_printf("--nordre = %i\n", *nordre);
  bft_printf("--idistu = %i\n", *idistu);
  bft_printf("--idiffl = %i\n", *idiffl);
  bft_printf("--modcpl = %i\n", *modcpl);
  bft_printf("--idirla = %i\n", *idirla);

  bft_printf("--ivisv1 = %i\n", *ivisv1);
  bft_printf("--ivisv2 = %i\n", *ivisv2);
  bft_printf("--ivistp = %i\n", *ivistp);
  bft_printf("--ivisdm = %i\n", *ivisdm);
  bft_printf("--iviste = %i\n", *iviste);
  bft_printf("--ivismp = %i\n", *ivismp);

  if (*iphyla == 2) {
    bft_printf("--ivisdk  = %i\n", *ivisdk);
    bft_printf("--iviswat = %i\n", *iviswat);
    bft_printf("--ivisch  = %i\n", *ivisch);
    bft_printf("--ivisck  = %i\n", *ivisck);
  }

  bft_printf("--isuist = %i\n", *isuist);
  bft_printf("--nbclst = %i\n", *nbclst);

  bft_printf("--istala = %i\n", *istala);
  if (*istala == 1) {
    bft_printf("--idstnt = %i\n", *idstnt);
    bft_printf("--nstist = %i\n", *nstist);
    bft_printf("--seuil  = %f\n", *seuil);

    /*
    bft_printf("--i        nomlag             nomlav              ihslag\n");
    for (i=1; i <= 5; i++)
      bft_printf("  %i %30s %30s %5i\n", i, nomlag[i-1], nomlav[i-1], ihslag[i-1]);
    i = 5;
    if (*iphyla == 1) {
      if (*itpvar == 1) {
        i++;
        bft_printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      }
      if (*idpvar == 1) {
        i++;
        bft_printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      }
    }
    else if (*iphyla == 2) {
      //i++;
      //bft_printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      i++;
      bft_printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      i++;
      bft_printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
      i++;
      bft_printf("  %i %s %s \n", i, nomlag[i], nomlav[i]);
    }
    i++;
    bft_printf("  %i %s \n", i, nomlag[i]);
    */
  }

  bft_printf("--iensi3 = %i\n", *iensi3);
  if (*iensi3 == 1) {
    bft_printf("--nstbor   = %i\n", *nstbor);
    bft_printf("--seuilf   = %f\n", *seuilf);
    bft_printf("--inbrbd   = %i\n", *inbrbd);
    bft_printf("--iflmbd   = %i\n", *iflmbd);
    bft_printf("--iangbd   = %i\n", *iangbd);
    bft_printf("--ivitbd   = %i\n", *ivitbd);
    bft_printf("--iencnbbd = %i\n", *iencnbbd);
    bft_printf("--iencmabd = %i\n", *iencmabd);
    bft_printf("--iencdibd = %i\n", *iencdibd);
    bft_printf("--iencckbd = %i\n", *iencckbd);
  }

#endif

}

/*-----------------------------------------------------------------------------
 * Fortran Interface:
 *
 * subroutine uilag2
 * *****************
 *
 * integer          nfabor  -->  number of boundary faces
 * integer          nozppm  -->  max number of boundary conditions zone
 * integer          nclagm  -->  max number of classes
 * integer          nflagm  -->  max number of boundaries
 * integer          iphyla  -->  physical model associated to the particles
 * ..
 * integer          nlayer  <--  number of layer for coal
 * integer          inuchl  <--  particle coal number
 * integer          irawcl  <--  coal particle composition injection condition
 * integer          ihpt    <--  coal temperature in K (for each layer)
 * integer          ifrlag  -->  type of boundary face
 * integer          iusncl  <--  array for particles class(es) number
 * integer          iusclb  <--  array for particles boundary conditions
 *----------------------------------------------------------------------------*/


void CS_PROCF (uilag2, UILAG2) (const int *const nfabor,
                                const int *const nozppm,
                                const int *const nbclst,
                                const int *const ientrl,
                                const int *const isortl,
                                const int *const idepo1,
                                const int *const idepo2,
                                const int *const idepfa,
                                const int *const iencrl,
                                const int *const irebol,
                                const int *const isymtl,
                                const int *const iphyla,
                                const int *const ijnbp,
                                const int *const ijfre,
                                const int *const iclst,
                                const int *const ijuvw,
                                const int *const iuno,
                                const int *const iupt,
                                const int *const ivpt,
                                const int *const iwpt,
                                const int *const ijprpd,
                                const int *const ipoit,
                                const int *const idebt,
                                const int *const ijprdp,
                                const int *const idpt,
                                const int *const ivdpt,
                                const int *const iropt,
                                const int *const ijprtp,
                                const int *const itpt,
                                const int *const icpt,
                                const int *const iepsi,
                                const int *const nlayer,
                                const int *const inuchl,
                                const int *const irawcl,
                                const int  const ihpt[],
                                int     ifrlag[],
                                int     iusncl[],
                                int     iusclb[])
{
  int izone, zones;
  int iclas, ilayer, icoal;
  int ielt, ifac, nelt = 0;
  char *interaction = NULL;
  char sclass[10];
  char *path1, *path2;
  char *choice;
  int *faces_list = NULL;

  cs_int_t  i_cz_params[20];  /* size: current ndlaim (=10) + margin */
  cs_real_t r_cz_params[100+4* *nlayer]; /* size: current ndlagm (=50+4*nlayer) + margin */

  zones = cs_gui_boundary_zones_number();

#if _XML_DEBUG_
  bft_printf("==>UILAG2\n");
#endif

  /* First iteration only: memory allocation */

  /*
    if (boundaries == NULL)
      _init_boundaries(nfabor, nozppm);
  */

  for (izone=0; izone < zones; izone++) {

    faces_list = cs_gui_get_faces_list(izone,
                                       boundaries->label[izone],
                                       *nfabor, *nozppm, &nelt);

    for ( ielt=0; ielt < nelt; ielt++ ) {
      ifac = faces_list[ielt];
      ifrlag[ifac-1] = izone+1;
    }

    path2 = cs_xpath_init_path();
    cs_xpath_add_elements(&path2, 2, "boundary_conditions", boundaries->nature[izone]);
    cs_xpath_add_test_attribute(&path2, "label", boundaries->label[izone]);
    cs_xpath_add_element(&path2, "particles");

    BFT_MALLOC(path1, strlen(path2)+1, char);
    strcpy(path1, path2);
    cs_xpath_add_attribute(&path1, "choice");
    interaction = cs_gui_get_attribute_value(path1);

    if (interaction != NULL) {

      if (cs_gui_strcmp(interaction, "inlet"))
        iusclb[izone] = *ientrl;

      else if(cs_gui_strcmp(interaction, "outlet"))
        iusclb[izone] = *isortl;

      else if(cs_gui_strcmp(interaction, "bounce"))
        iusclb[izone] = *irebol;

      else if(cs_gui_strcmp(interaction, "part_symmetry"))
        iusclb[izone] = *isymtl;

      else if(cs_gui_strcmp(interaction, "deposit1"))
        iusclb[izone] = *idepo1;

      else if(cs_gui_strcmp(interaction, "deposit2"))
        iusclb[izone] = *idepo2;

      else if(cs_gui_strcmp(interaction, "fouling") && *iphyla == 2)
        iusclb[izone] = *iencrl;

      else if(cs_gui_strcmp(interaction, "fouling") && (*iphyla == 0  || *iphyla == 1))
        iusclb[izone] = *idepfa;

#if _XML_DEBUG_
      bft_printf("--iusclb[%i] = %i has %i class(es) \n", izone, iusclb[izone], iusncl[izone]);

      bft_printf("--zone %i : class number %i \n", izone, iusncl[izone]);
      bft_printf("--        : label    %s \n", boundaries->label[izone]);
      bft_printf("--        : nature   %s \n", boundaries->nature[izone]);
      bft_printf("--        : p_nature %i \n", iusclb[izone]);
#endif

      /* Additional info for inlet */

      if (iusclb[izone] == *ientrl) {

        strcpy(path1, path2);
        cs_xpath_add_element(&path1, "class");
        iusncl[izone] = cs_gui_get_nb_element(path1);
        strcpy(path1, path2);

        for (iclas=0; iclas < iusncl[izone]; iclas++) {

          cs_lagr_init_zone_class_param(i_cz_params, r_cz_params);

          sprintf(sclass, "class[%i]", iclas+1);
          BFT_REALLOC(path2,
                      ( 20+strlen(boundaries->nature[izone])
                       +10+strlen(boundaries->label[izone])
                       +13+strlen(sclass)+1),
                      char);
          strcpy(path2, "");
          sprintf(path2,
                  "boundary_conditions/%s[@label='%s']/particles/%s",
                  boundaries->nature[izone],
                  boundaries->label[izone],
                  sclass);

          _get_int(&(i_cz_params[*ijnbp -1]), 2, path2, "number");
          _get_int(&(i_cz_params[*ijfre -1]), 2, path2, "frequency");
          _get_int(&(i_cz_params[*iclst -1]), 2, path2, "statistical_groups");

          /* velocity */

          choice = _get_attr("choice", 2, path2, "velocity");

          if (cs_gui_strcmp(choice, "fluid"))
            i_cz_params[*ijuvw -1] = -1;

          else if (cs_gui_strcmp(choice, "norm")) {
            i_cz_params[*ijuvw -1] = 0;
            _get_double(&(r_cz_params[*iuno -1]), 3, path2, "velocity", "norm");
          }
          else if (cs_gui_strcmp(choice, "components")) {
            i_cz_params[*ijuvw -1] = 1;
            _get_double(&(r_cz_params[*iupt -1]), 3, path2, "velocity", "velocity_x");
            _get_double(&(r_cz_params[*ivpt -1]), 3, path2, "velocity", "velocity_y");
            _get_double(&(r_cz_params[*iwpt -1]), 3, path2, "velocity", "velocity_z");
          }
          else if (cs_gui_strcmp(choice, "subroutine"))
            i_cz_params[*ijuvw -1] = 2;

          BFT_FREE(choice);

          /* statistical_weight, mass_flow_rate*/

          choice = _get_attr("choice", 2, path2, "statistical_weight");

          if (cs_gui_strcmp(choice, "prescribed")) {
            i_cz_params[*ijprpd -1] = 1;
            _get_double(&(r_cz_params[*ipoit -1]), 2, path2, "statistical_weight");
            r_cz_params[*idebt -1] = 0;
          }
          else if (cs_gui_strcmp(choice, "rate")) {
            i_cz_params[*ijprpd -1] = 1;
            _get_double(&(r_cz_params[*idebt -1]), 2, path2, "mass_flow_rate");
            r_cz_params[*ipoit -1] = 1;
          }
          else if (cs_gui_strcmp(choice, "subroutine")) {
            i_cz_params[*ijprpd -1] = 2;
            _get_double(&(r_cz_params[*ipoit -1]), 2, path2, "statistical_weight");
            r_cz_params[*idebt -1] = 0;
          }

          BFT_FREE(choice);

          /* diameter */

          choice = _get_attr("choice", 2, path2, "diameter");

          if (cs_gui_strcmp(choice, "prescribed")) {
            i_cz_params[*ijprdp -1] = 1;
            _get_double(&(r_cz_params[*idpt -1]), 2, path2, "diameter");
            _get_double(&(r_cz_params[*ivdpt -1]), 2, path2, "diameter_standard_deviation");
          }
          else if (cs_gui_strcmp(choice, "subroutine"))
            i_cz_params[*ijprdp -1] = 2;

          BFT_FREE(choice);

          /* density */
          if (*iphyla != 2)
            _get_double(&(r_cz_params[*iropt -1]), 2, path2, "density");

          if (*iphyla == 1) {

            /* temperature, specific_heat, emissivity */

            choice = _get_attr("choice", 2, path2, "temperature");

            if (cs_gui_strcmp(choice, "prescribed")) {
              i_cz_params[*ijprtp -1] = 1;
              _get_double(&(r_cz_params[*itpt -1]), 2, path2, "temperature");
            }
            else if (cs_gui_strcmp(choice, "subroutine"))
              i_cz_params[*ijprtp -1] = 2;

            _get_double(&(r_cz_params[*icpt -1]), 2, path2, "specific_heat");
            _get_double(&(r_cz_params[*iepsi -1]), 2, path2, "emissivity");

            BFT_FREE(choice);

          }

          /* coal */

          if (*iphyla == 2) {
            /* Read the coal number */
            _get_int(&(i_cz_params[*inuchl -1]), 2, path2, "coal_number");
            icoal = i_cz_params[*inuchl -1];

            // Fresh coal or user defined
            choice  = _get_attr("choice", 2, path2, "coal_composition");

            if (cs_gui_strcmp(choice, "raw_coal_as_received")) {

              /* Data are read in pulverised fuel combustion module */
              i_cz_params[*irawcl -1]  = 1;

              for (ilayer=0; ilayer < *nlayer; ilayer++)
                _get_double(&(r_cz_params[ihpt[ilayer] -1]), 2, path2, "coal_temperature");

            }
            else if (cs_gui_strcmp(choice, "subroutine")) {
              i_cz_params[*irawcl-1]  = 0;
            }

            BFT_FREE(choice);
          }

          /* Complete class paramaters definition */

          cs_lagr_define_zone_class_param(iclas+1, izone+1,
                                          i_cz_params, r_cz_params);

#if _XML_DEBUG_

          bft_printf("---number = %i \n", i_cz_params[*ijnbp -1]);
          bft_printf("---frequency = %i \n", i_cz_params[*ijfre -1]);
          bft_printf("---statistical_groups = %i \n", i_cz_params[*iclst -1]);

          bft_printf("---velocity choice: %i  (-1: fluid, 0: norm, 1: components, 2: subroutine)\n", i_cz_params[*ijuvw -1]);

          if (i_cz_params[*ijuvw -1] == 0)

            bft_printf("----norm = %f \n", r_cz_params[*iuno -1]);

          else if (i_cz_params[*ijuvw -1] == 1) {

            bft_printf("----u = %f \n", r_cz_params[*iupt -1]);
            bft_printf("----v = %f \n", r_cz_params[*ivpt -1]);
            bft_printf("----w = %f \n", r_cz_params[*iwpt -1]);
          }

          bft_printf("---statistical weight choice: %i  (1: prescribed, 2: subroutine)\n", i_cz_params[*ijprpd -1]);

          if (i_cz_params[*ijprpd -1] == 1) {
            bft_printf("----statistical weight = %f \n", r_cz_params[*ipoit -1]);
            bft_printf("----mass flow rate = %f \n", r_cz_params[*idebt -1]);
          }

          bft_printf("---diameter choice = %i (1: prescribed, 2: subroutine)\n", i_cz_params[*ijprdp -1]);

          if (i_cz_params[*ijprdp -1] == 1) {
            bft_printf("----diameter = %f \n", r_cz_params[*idpt -1]);
            bft_printf("----standard deviation = %f \n", r_cz_params[*ivdpt -1]);
          }

          if (*iphyla != 2)
            bft_printf("---density = %f \n", r_cz_params[*iropt -1]);

          if (*iphyla == 1) {

            bft_printf("---temperature choice = %i (1: prescribed, 2: subroutine)\n", i_cz_params[*ijprtp -1]);

            if (i_cz_params[*ijprtp -1] == 1)
              bft_printf("----temperature = %f \n", r_cz_params[*itpt -1]);

            bft_printf("---specific heat = %f \n", r_cz_params[*icpt -1]);
            bft_printf("---emissivity = %f \n", r_cz_params[*iepsi -1]);
          }

          if (*iphyla == 2) {
            bft_printf("---coal number = %i \n",            i_cz_params[*inuchl -1]);
            bft_printf("---coal composition = %i (1: raw coal, 2: user defined)\n", i_cz_params[*rawcl -1]);
          }

#endif /* _XML_DEBUG_ */

        } /* End of loop on class */

      } /* End of test on inlet */

    } /* End of loop on zones */

    BFT_FREE(path1);
    BFT_FREE(path2);
    BFT_FREE(faces_list);
    BFT_FREE(interaction);
  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Free global GUI structures related to particles.
 *----------------------------------------------------------------------------*/

void
cs_gui_particles_free(void)
{
  int i;

  for (i = 0; i < _last_mean_var; i++)
    BFT_FREE(_array_mean_varname[i]);
  BFT_FREE(_array_mean_varname);
  _max_mean_vars = 0;
  _last_mean_var = 0;

  for (i = 0; i < _last_variance_var; i++)
    BFT_FREE(_array_variance_varname[i]);
  BFT_FREE(_array_variance_varname);
  _max_variance_vars = 0;
  _last_variance_var = 0;

  for (i = 0; i < _last_boundary_var; i++)
    BFT_FREE(_array_boundary_varname[i]);
  BFT_FREE(_array_boundary_varname);
  _max_boundary_vars = 0;
  _last_boundary_var = 0;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
