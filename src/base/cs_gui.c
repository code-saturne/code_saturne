/*============================================================================
 *
 *                    Code_Saturne version 1.3
 *                    ------------------------
 *
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2007 EDF S.A., France
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
 * Reader of the parameters file: main parameters, boundary conditions
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
 * libxml2 library headers
 *----------------------------------------------------------------------------*/

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "fvm_selector.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_gui.h"

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
 * Variables and scalars management structure
 *----------------------------------------------------------------------------*/

typedef struct {
  char  *model;            /* particular physical model                       */
  char  *model_value;      /* particular physical model value                 */
  char **head;             /* name of the head                                */
  char **type;             /* type of markup: 'variable' or 'scalar'          */
  char **name;             /* variables and scalars label                     */
  char **label;            /* scalars label                                   */
  int   *rtp;              /* variables position in fortran array RTP         */
  int    nvar;             /* total number of variables and scalars           */
  int    nscaus;           /* user scalar number                              */
  int    nscapp;           /* predifined physics scalar number                */
  int    nprop;            /* proprietes number                               */
  int    nsalpp;           /* particular physical proprietes number           */
  int    ntimaver;         /* time averages number                            */
  char **properties_name;  /* label of properties                             */
  int   *properties_ipp;   /* properties position in fortran array PROPCE     */
} cs_var_t;

/*----------------------------------------------------------------------------
 * Structures associated to boundary conditions definition
 *----------------------------------------------------------------------------*/

typedef struct {
  double val1;             /* fortran array RCODCL(.,.,1) mapping             */
  double val2;             /* fortran array RCODCL(.,.,2) mapping             */
  double val3;             /* fortran array RCODCL(.,.,3) mapping             */
} cs_val_t;

typedef struct {
  char      **label;       /* label for each boundary zone                    */
  char      **nature;      /* nature for each boundary zone                   */
  int        *iqimp;       /* 1 if a flow rate is applied                     */
  int        *ientat;      /* 1 if boundary is an inlet for air (coal)        */
  int        *ientcp;      /* 1 if boundary is an inlet for coal              */
  int        *icalke;      /* automatic boundaries for turbulent variables    */
  double     *qimp;        /* inlet flow rate                                 */
  double     *qimpat;      /* inlet flow rate (coal combustion)               */
  double     *timpat;      /* inlet air temperature (coal combustion)         */
  double    **qimpcp;      /* inlet coal flow rate (coal combustion)          */
  double    **timpcp;      /* inlet coal temperature (coal)                   */
  double     *dh;          /* inlet hydraulic diameter                        */
  double     *xintur;      /* inlet turbulent intensity                       */
  int       **type_code;   /* type of boundary for each variables             */
  cs_val_t  **values;      /* fortran array RCODCL mapping                    */
  double   ***distch;      /* ratio for each coal                             */
  double     *rough;       /* roughness size                                  */
} cs_boundary_t;

/*============================================================================
 *  External global variables
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Gestion du document xml
 *----------------------------------------------------------------------------*/

extern xmlXPathContextPtr xpathCtx;   /* Pointer on the Context       */
extern xmlNodePtr node;               /* Pointer on the root node     */

/*============================================================================
 * Private global variables
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Private global variables for the treatment of NOMVAR.
 * NOMVAR is a characters fortran array
 *----------------------------------------------------------------------------*/

static int      _cs_gui_max_vars = 0;
static int      _cs_gui_last_var = 0;
static char  ** _cs_gui_var_name = NULL;

/*----------------------------------------------------------------------------
 * Private global variables for boundaru conditions
 *----------------------------------------------------------------------------*/

static cs_boundary_t *boundaries = NULL;
static cs_var_t *vars = NULL;

/*============================================================================
 * Private functions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Copy a variable name to private variable names array
 *
 * parameters:
 *   varname        -->  name or label of the variable/scalar/property
 *   ipp            -->  index from the fortran array associated to varname
 *----------------------------------------------------------------------------*/

static void _gui_copy_varname(const char *varname, int ipp)
{
  size_t l;

  if (ipp < 1 || ipp > _cs_gui_last_var)
    bft_error(__FILE__, __LINE__, 0,
              _("Variable index %d out of bounds (1 to %d)"),
                 ipp, _cs_gui_last_var);

  l = strlen(varname);

  if (_cs_gui_var_name[ipp-1] == NULL)
    BFT_MALLOC(_cs_gui_var_name[ipp-1], l + 1, char);

  else if (strlen(_cs_gui_var_name[ipp-1]) != l)
    BFT_REALLOC(_cs_gui_var_name[ipp-1], l + 1, char);

  strcpy(_cs_gui_var_name[ipp-1], varname);
}

/*----------------------------------------------------------------------------
 * Turbulence model parameters.
 *
 * parameters:
 *   param                -->  name of the parameters
 *   keyword             <--   turbulence model parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_advanced_options_turbulence(const char *const param,
                                                      int *const keyword)
{
  char *path = NULL;
  int  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models", "turbulence", param);

  if (cs_gui_strcmp("gravity_terms", param)) {

    cs_xpath_add_attribute(&path, "status");
    if (cs_gui_get_status(path, &result)) *keyword = result;

  } else if (cs_gui_strcmp("scale_model", param)) {

    cs_xpath_add_function_text(&path);
    if (cs_gui_get_int(path, &result)) *keyword = result;

  } else
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s \n"), path);

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return user scalar number.
 *----------------------------------------------------------------------------*/

static int cs_gui_get_number_user_scalar(void)
{
  char *path = NULL;
  int   nb;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "additional_scalars", "scalar");
  nb = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return nb;
}

/*-----------------------------------------------------------------------------
 * Return the activated particular physics scalar number
 *----------------------------------------------------------------------------*/

static int cs_gui_model_scalar_number(const char* model)
{
  char *path = NULL;
  int   nb;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "scalar");

  nb = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return nb;
}

/*-----------------------------------------------------------------------------
 * Return the name of the related scalar if the scalar "num_sca" is a variance
 *
 * parameter:
 *   num_sca           -->  scalar number
 *----------------------------------------------------------------------------*/

static char *cs_gui_scalar_variance(const int num_sca)
{
  char *path = NULL;
  char *variance = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "scalar", num_sca);
  cs_xpath_add_element(&path, "variance");
  cs_xpath_add_function_text(&path);

  variance = cs_gui_get_text_value(path);

  BFT_FREE(path);

  return variance;
}

/*-----------------------------------------------------------------------------
 * Return the user thermal scalar indicator.
 *----------------------------------------------------------------------------*/

static int cs_gui_thermal_scalar(void)
{
  char *model_name = NULL;
  int   test;

  model_name = cs_gui_get_thermophysical_model("thermal_scalar");

  if (cs_gui_strcmp(model_name, "off"))
    test = 0;
  else {
    if (cs_gui_strcmp(model_name, "enthalpy"))
      test =  2 ;
    else if (cs_gui_strcmp(model_name, "temperature_kelvin"))
      test =  1 ;
    else if (cs_gui_strcmp(model_name, "temperature_celsius"))
      test = -1 ;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Invalid thermal model: %s\n"), model_name);
    }

  BFT_FREE(model_name);

  return test;
}

/*----------------------------------------------------------------------------
 * Get thermal user scalar number if it is exist.
 *
 * parameters:
 *   iscalt               <--  thermal scalar number order
 *   iscsth               <--  nature of the thermal scalar (C, K, J/kg)
 *----------------------------------------------------------------------------*/

static void cs_gui_thermal_scalar_number(int *const iscalt,
                                         int *const iscsth)
{
  int ind_thermal;
  int i, index, size;
  int iphas = 0;
  char *path = NULL;
  char **name = NULL;

  ind_thermal = cs_gui_thermal_scalar();

  if (ind_thermal) {
    path = cs_xpath_init_path();
    cs_xpath_add_elements(&path, 2, "additional_scalars", "/@type");
    name = cs_gui_get_attribute_values(path, &size);

    index = -1;
    for (i=0; i < size; i++) {
      if (cs_gui_strcmp(name[i], "thermal")) index = i;
    }
    iscalt[iphas] = index+1;
    iscsth[index] = ind_thermal;

    BFT_FREE(path);
    for (i=0; i < size; i++) BFT_FREE(name[i]);
    BFT_FREE(name);
  }
}

/*-----------------------------------------------------------------------------
 * Return the value of choice for user scalar's property
 *
 * parameters:
 *   scalar_num     --> number of scalar
 *   property_name  --> name of property
 *   choice         --> choice for property
 *----------------------------------------------------------------------------*/

static int cs_gui_scalar_properties_choice(const int         scalar_num,
                                           const char *const property_name,
                                                 int  *const choice)
{
  char *path = NULL;
  char *buff = NULL;
  int   ichoice;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "scalar", scalar_num);
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_attribute(&path, "choice");

  buff = cs_gui_get_attribute_value(path);

  if (buff == NULL) {
    ichoice = 0;

  } else {
    ichoice = 1;

    if (cs_gui_strcmp(buff, "variable"))
      *choice = 1;
    else if (cs_gui_strcmp(buff, "constant"))
      *choice = 0;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  }

  BFT_FREE(path);
  BFT_FREE(buff);

  return ichoice;
}

/*-----------------------------------------------------------------------------
 * Return value of diffusion coefficient for user scalars
 *        return 1 if value exists
 *        return 0 if not
 *
 * parameters:
 *   num_sca  --> number of scalar
 *   value   <--  value of diffusion coefficient
 *----------------------------------------------------------------------------*/

static void cs_gui_scalar_diffusion_value(const int           num_sca,
                                                double *const value)
{
  char  *path = NULL;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "scalar", num_sca);
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_element(&path, "initial_value");
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *value = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return value for iale method
 *
 * parameters:
 *   param               -->  iale parameter
 *   keyword             <--  value of the iale parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_iale_parameter(const char   *const param,
                                        double *const keyword)
{
  char   *path   = NULL;
  char   *type = NULL;
  double  result = 0.0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models", "ale_method", param);

  if (cs_gui_strcmp(param,"mesh_viscosity") ){

    cs_xpath_add_attribute(&path, "type");
    type = cs_gui_get_attribute_value(path);
    if(cs_gui_strcmp(type, "isotrop"))
      *keyword = 0;
    else if (cs_gui_strcmp(type, "orthotrop"))
      *keyword = 1;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  } else {

    cs_xpath_add_function_text(&path);
    if (cs_gui_get_double(path, &result)) *keyword = result;

  }
  BFT_FREE(type);
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Get the status of steady management.
 *
 * parameter:
 *   keyword         <--  if 1 unsteady management else steady management
 *----------------------------------------------------------------------------*/

static void cs_gui_get_steady_status(int *const keyword)
{
  char *path = NULL;
  int   result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "steady_management");
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
    *keyword = result;
  else
    *keyword = 1;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return  parameters for steady management.
 *
 * parameter:
 *   param           -->  steady parameter
 *   keyword         <--  new value for the steady parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_steady_parameters(const char   *const param,
                                           double *const keyword)
{
  char   *path   = NULL;
  double  result = 0.0;
  int     status = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "steady_management", param);

  if (cs_gui_strcmp(param,"zero_iteration")){

    cs_xpath_add_attribute(&path, "status");
    if(cs_gui_get_status(path, &status)) *keyword = status;

  } else {

    cs_xpath_add_function_text(&path);
    if (cs_gui_get_double(path, &result)) *keyword = result;

  }
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify time parameters.
 *
 * parameters:
 *   param              -->  time parameter
 *   keyword            <--  new value of the time parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_time_parameters(const char   *const param,
                                         double *const keyword)
{
  char   *path   = NULL;
  double  result = 0.0;
  int     status = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "time_parameters", param);

  if (cs_gui_strcmp(param,"zero_time_step") ||
      cs_gui_strcmp(param,"thermal_time_step")) {

    cs_xpath_add_attribute(&path, "status");
    if(cs_gui_get_status(path, &status)) *keyword = status;

  } else {

    cs_xpath_add_function_text(&path);
    if (cs_gui_get_double(path, &result)) *keyword = result;

  }
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify restart files format.
 *
 * parameters:
 *   param               -->  restart file name
 *   keyword            <-->  new value of the restart file format
 *----------------------------------------------------------------------------*/

static void cs_gui_restart_parameters_file_format (const char *const param,
                                                          int *const format)
{
  char *path = NULL;
  char *result = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "calcul_management", "start_restart", param);
  cs_xpath_add_attribute(&path, "format");

  result = cs_gui_get_attribute_value(path);

  if (result != NULL) {
    if (cs_gui_strcmp(result, "binary"))
      *format = 0;
    else if (cs_gui_strcmp(result, "ascii"))
      *format = 1;
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Restart file format unknown: %s.\nXpath: %s\n"),
                result, path);

  }

  BFT_FREE(result);
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify restart parameters.
 *
 * parameters:
 *   param               -->  restart parameter
 *   keyword            <-->  new value of the restart parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_restart_parameters_status(const char *const param,
                                                    int *const keyword)
{
  int   result;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "calcul_management", "start_restart", param);
  cs_xpath_add_attribute(&path, "status");

  if(cs_gui_get_status(path, &result))
    *keyword = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return value of numerical parameter markup
 *
 * parameters:
 *   variable_type  --> name of variable
 *   value_type     --> name of numerical parameter parkup
 *   value          <-- value of numerical parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_variable_value(const char   *const variable_type,
                                  const char   *const value_type,
                                        double *const value)
{
  char  *path = NULL;
  double result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", variable_type);
  cs_xpath_add_element(&path, value_type);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result)) *value = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the attribute value from the xpath query.
 *
 * parameters:
 *   path          --> path for xpath query
 *   child         --> child markup
 *   keyword      <--  value of attribute node
 *----------------------------------------------------------------------------*/

static void _attribute_value(      char *      path,
                             const char *const child,
                                   int  *const keyword)
{
  char *choice = NULL;
  int   result;

  assert(path != NULL);
  assert(child != NULL);

  if (cs_gui_strcmp(child, "order_scheme")) {

    /* *keyword = 1; */
    cs_xpath_add_attribute(&path, "choice");
    choice = cs_gui_get_attribute_value(path);

    if (cs_gui_strcmp(choice, "centered"))
      *keyword = 1;
    else if (cs_gui_strcmp(choice, "solu"))
      *keyword = 0;
    BFT_FREE(choice);

  } else {

    cs_xpath_add_attribute(&path, "status");

    if (cs_gui_get_status(path, &result)) {
      *keyword = result;

      if (cs_gui_strcmp(child, "slope_test")) {
        if (result == 1) *keyword = 0;
        if (result == 0) *keyword = 1;
      }
    } else {

      if (cs_gui_strcmp(child, "postprocessing_recording") ||
          cs_gui_strcmp(child, "listing_printing")) *keyword = 1;

    }
  }
  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the attribute value associated to a child markup from a variable.
 *
 * parameters:
 *   name          -->  name of the variable markup
 *   child         -->  child markup
 *   keyword      <--   value of attribute node contained in the child markup
 *----------------------------------------------------------------------------*/

static void cs_gui_variable_attribute(const char *const name,
                                      const char *const child,
                                            int  *const keyword)
{
  char *path = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, child);

  _attribute_value(path, child, keyword);
}

/*----------------------------------------------------------------------------
 * Get the text value associated to a child markup from a scalar.
 *
 * parameters:
 *   label                -->  label of the scalar markup
 *   child                -->  name of the child markup
 *   value               <--   value of text node contained in the child markup
 *----------------------------------------------------------------------------*/

static void cs_gui_scalar_value(const char   *const label,
                                const char   *const child,
                                      double *const value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, child);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result)) *value = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the attribute value associated to a child markup from a scalar.
 *
 * parameters:
 *   label         -->  name of the scalar markup
 *   child         -->  child markup
 *   keyword      <--   value of attribute node contained in the child markup
 *----------------------------------------------------------------------------*/

static void cs_gui_scalar_attribute(const char *const label,
                                    const char *const child,
                                          int  *const keyword)
{
  char *path = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, child);

  _attribute_value(path, child, keyword);
}

/*-----------------------------------------------------------------------------
 * Donne les valeurs liees aux scalaire model : min, max ...
 *       La fonction retourne 1 si la valeur existe
 *                            0 sinon
 * le resultat est stocke dans value
 *----------------------------------------------------------------------------*/


static void cs_gui_model_scalar_value(const   char *const model,
                                      const   char *const name,
                                      const   char *const keyword,
                                            double *const value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, keyword);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path,&result))
    *value = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the attribute value associated to a child markup from a scalar.
 *
 * parameters:
 *   model         -->  model markup
 *   name          -->  name of the scalar markup
 *   child         -->  child markup
 *   keyword      <--   value of attribute node contained in the child markup
 *----------------------------------------------------------------------------*/

static void cs_gui_model_scalar_output_status(const char *const model,
                                              const char *const name,
                                              const char *const child,
                                                    int  *const keyword)
{
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, child);

  _attribute_value(path, child, keyword);
}

/*-----------------------------------------------------------------------------
 * Modify double numerical parameters.
 *
 * parameters:
 *   param               -->  label of the numerical parameter
 *   keyword            <-->  value of the numerical parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_numerical_double_parameters(const char   *const param,
                                                     double *const keyword)
{
  char  *path = NULL;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "numerical_parameters");
  cs_xpath_add_element(&path, param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *keyword = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify integer numerical parameters.
 *
 * parameters:
 *   param               -->  label of the numerical parameter
 *   keyword            <-->  value of the numerical parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_numerical_int_parameters(const char *const param,
                                           int  *const keyword)
{
  char *path = NULL;
  char *choice = NULL;
  int   result;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "numerical_parameters");
  cs_xpath_add_element(&path, param);

  if (cs_gui_strcmp(param, "gradient_reconstruction")){

    cs_xpath_add_attribute(&path, "choice");
    choice = cs_gui_get_attribute_value(path);
    if (choice) *keyword = atoi(choice);
    BFT_FREE(choice);

  } else {

    cs_xpath_add_attribute(&path, "status");
    if (cs_gui_get_status(path, &result)) *keyword = result;

  }
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Modify gravity parameters.
 *
 * parameters:
 *   param               -->  gravity parameter (GX, GY, GZ)
 *   keyword            <-->  new value of the gravity parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_gravity_value(const char   *const param,
                                       double *const value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "physical_properties", "gravity", param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result)) *value = result;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Get initial value from property markup.
 *
 * parameters:
 *   property_name       -->  name of the property
 *   value              <--   new initial value of the property
 *----------------------------------------------------------------------------*/

static void  cs_gui_properties_value(const char   *const property_name,
                                           double *const value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", property_name);
  cs_xpath_add_element(&path, "initial_value");
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result)) *value = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the value of the choice attribute from a property markup.
 * Return 1 if the xpath request has succeeded, 0 otherwise.
 *
 * parameters:
 *   property_name        -->  name of the property
 *   choice              <--   value of the attribute choice
 *----------------------------------------------------------------------------*/

static int cs_gui_properties_choice(const char *const property_name,
                                          int  *      choice)
{
  char *path = NULL;
  char *buff = NULL;
  int   iok;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", property_name);
  cs_xpath_add_attribute(&path, "choice");

  buff = cs_gui_get_attribute_value(path);

  if (buff == NULL)
    iok = 0;

  else {
    iok = 1;

    if (cs_gui_strcmp(buff, "variable"))
      *choice = 1;
    else if (cs_gui_strcmp(buff, "constant"))
      *choice = 0;
    else
      bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);
  }

  BFT_FREE(buff);
  BFT_FREE(path);

  return iok;
}

/*-----------------------------------------------------------------------------
 * Get reference value of pressure
 *
 * parameters:
 *   p0              <--   value of pressure
 *----------------------------------------------------------------------------*/

static void cs_gui_reference_pressure(double *const p0)
{
  char *path = NULL;
  double value;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "reference_pressure");
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &value)) *p0 = value;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Get reference value of temperature
 *
 * parameters:
 *   model           -->   name of activated model
 *   t0              <--   value of temperature
 *----------------------------------------------------------------------------*/

static void cs_gui_reference_temperature(char *const model, double *const t0)
{
  char *path = NULL;
  double value;

  path = cs_xpath_short_path();
  cs_xpath_add_elements(&path, 2, model,"reference_temperature");
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &value)) *t0 = value;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Get reference value of mass molar molecular
 *
 * parameters:
 *   model           -->   name of activated model
 *   m0              <--   value of mass molar molecular
 *----------------------------------------------------------------------------*/

static void cs_gui_reference_mass_molar(char *const model, double *const m0)
{
  char *path = NULL;
  double value;

  path = cs_xpath_short_path();
  cs_xpath_add_elements(&path, 2, model,"reference_mass_molar");
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &value)) *m0 = value;
  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Turbulence model parameters.
 *
 * parameters:
 *   param                -->  name of the parameters
 *   keyword             <--   turbulence model parameter
 *----------------------------------------------------------------------------*/

static void _option_turbulence_double(const char *const param,
                                          double *const keyword)
{
  char *path = NULL;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models", "turbulence", param);

  cs_xpath_add_function_text(&path);
  if (cs_gui_get_double(path, &result)) *keyword = result;

  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Initialization choice of the turbulence variables parameters.
 *
 * parameters:
 *   param                -->  name of the parameters
 *   value               <--   initialization choice
 *----------------------------------------------------------------------------*/

static void cs_gui_turbulence_initialization(const char   *const param,
                                            double *const value)
{
  char   *path = NULL;
  double  result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 4,
                        "thermophysical_models",
                        "turbulence",
                        "initialization",
                        param);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result)) *value = result;
  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Return the initialization choice of the turbulence variables.
 *----------------------------------------------------------------------------*/

static char *cs_gui_turbulence_initialization_choice(void)
{
  char *path = NULL;
  char *initialization_choice;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "thermophysical_models",
                        "turbulence",
                        "initialization");
  cs_xpath_add_attribute(&path, "choice");

  initialization_choice = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return initialization_choice;
}

/*================================
 * Entrees-sorties
 *===============================*/

/*----------------------------------------------------------------------------
 * Get output control value parameters.
 *
 * parameters:
 *   param                -->  name of the parameter
 *   keyword             <--   output control parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_output_value(const char *const param,
                               int  *const keyword)
{
  char *path = NULL;
  char *choice = NULL;
  int   result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "output", param);

  if (cs_gui_strcmp(param, "auxiliary_restart_file_writing") ||
      cs_gui_strcmp(param, "fluid_domain") ||
      cs_gui_strcmp(param, "domain_boundary") ||
      cs_gui_strcmp(param, "syrthes_boundary") ) {

    cs_xpath_add_attribute(&path, "status");
    if(cs_gui_get_status(path, &result)) *keyword = result;

  }else if (cs_gui_strcmp(param, "postprocessing_mesh_options")){
    cs_xpath_add_attribute(&path, "choice");
    choice = cs_gui_get_attribute_value(path);
    if (choice) *keyword = atoi(choice);
  } else {

    cs_xpath_add_function_text(&path);
    if (cs_gui_get_int(path, &result)) *keyword = result;

  }

  BFT_FREE(choice);
  BFT_FREE(path);
}

/*----------------------------------------------------------------------------
 * Get the output format and options for postprocessing.
 *
 * parameters:
 *   param                -->  name of the parameter
 *   keyword             <--   output control parameter
 *   size_key             -->  keyword string size
 *----------------------------------------------------------------------------*/

static void cs_gui_output_choice(const char *const param,
                                char *const keyword,
                          const int  *const size_key)
{
  char *path = NULL;
  char *choice = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "output", param);
  cs_xpath_add_attribute(&path, "choice");
  choice = cs_gui_get_attribute_value(path);

  if (choice != NULL) cs_gui_strcpy_c2f(keyword, choice, *size_key);

  BFT_FREE(choice);
  BFT_FREE(path);
}

/*==================================
 * TRAITMENTS FOR TIME AVERAGES
 *=================================*/

/*-----------------------------------------------------------------------------
 * Return the number of time_averages
 *----------------------------------------------------------------------------*/

static int cs_gui_get_means_number(void)
{
  char *path = NULL;
  int   number = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "time_averages", "time_average");
  number = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return number ;
}

/*----------------------------------------------------------------------------
 * Get list of variables or properties or scalar's names for calculation mean
 *
 * parameters:
 *   id           -->  number of mean (imom)
 *   list         <--  output control parameter
 *----------------------------------------------------------------------------*/

static int cs_gui_get_mean_names_number( int   const id)
{
  char *path = NULL;
  char *str_id = NULL;
  int   number = 0;

  BFT_MALLOC(str_id,
             cs_gui_characters_number(id)+1,
             char);
  sprintf(str_id, "%i", id);

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "time_averages", "time_average");
  cs_xpath_add_test_attribute(&path,"id",str_id);
  cs_xpath_add_element(&path, "var_prop");
  number = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return number;

}

/*----------------------------------------------------------------------------
 * Get mean value parameters.
 *
 * parameters:
 *   id              -->  number of mean (imom)
 *   param           -->  name of the parameter
 *   keyword         <--   output control parameter
 *----------------------------------------------------------------------------*/

static void cs_gui_get_mean_value(      int   const id,
                                  const char *const param,
                                        int  *const keyword)
{
  char *path = NULL;
  char *str_id = NULL;
  int   result = 0;

  BFT_MALLOC(str_id,
             cs_gui_characters_number(id)+1,
             char);
  sprintf(str_id, "%i", id);

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "time_averages", "time_average");
  cs_xpath_add_test_attribute(&path,"id",str_id);
  cs_xpath_add_element(&path, param);

  cs_xpath_add_function_text(&path);
  if (cs_gui_get_int(path, &result)) *keyword = result;

  BFT_FREE(path);
  BFT_FREE(str_id);
}

/*----------------------------------------------------------------------------
 * Get variable or properties or scalar's name for one mean
 *
 * parameters:
 *   id           -->  number of mean (imom)
 *   nb           -->  number of order in list of var_prop of the mean
 *----------------------------------------------------------------------------*/

static char *cs_gui_get_mean_prop( int const id, int const nb)
{
  char *path = NULL;
  char *name = NULL;
  char *str_id = NULL;

  BFT_MALLOC(str_id,
             cs_gui_characters_number(id)+1,
             char);
  sprintf(str_id, "%i", id);

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "time_averages", "time_average");
  cs_xpath_add_test_attribute(&path,"id",str_id);
  cs_xpath_add_element_num(&path, "var_prop", nb);
  cs_xpath_add_attribute(&path, "name");

  name = cs_gui_get_attribute_value(path);

  BFT_FREE(path);
  BFT_FREE(str_id);

  return name;
}

/*----------------------------------------------------------------------------
 * Get label of mean
 *
 * parameters:
 *   nb           -->  number of order in list of mean
 *----------------------------------------------------------------------------*/

static char *cs_gui_get_mean_label( int const nb)
{
  char *path = NULL;
  char *label = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "analysis_control", "time_averages");
  cs_xpath_add_element_num(&path, "time_average", nb);
  cs_xpath_add_attribute(&path,"label");

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return label;
}

/*===================
 * FOR PROBES
 *==================*/

/*-----------------------------------------------------------------------------
 * Return the number of the <probe> markups.
 *----------------------------------------------------------------------------*/

static int cs_gui_probes_number(void)
{
  char *path = NULL;
  int   number = 0;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "output", "probe");
  number = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return number ;
}

/*-----------------------------------------------------------------------------
 * Return a single coordinate of a monitoring probe
 *
 * parameters
 *   num_probe            -->  number aka name of the monitoring probe
 *   probe_coord          -->  one coordinate of the monitoring probe
 *----------------------------------------------------------------------------*/

static double cs_gui_probe_coordinate(const int         num_probe,
                                      const char *const probe_coord)
{
  char  *path = NULL;
  char  *str_num_probe = NULL;
  double result = 0.0;

  assert(num_probe>0);

  BFT_MALLOC(str_num_probe,
             cs_gui_characters_number(num_probe)+1,
             char);
  sprintf(str_num_probe, "%i", num_probe);


  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "analysis_control", "output", "probe");
  cs_xpath_add_test_attribute(&path, "name", str_num_probe);
  cs_xpath_add_element(&path, probe_coord);
  cs_xpath_add_function_text(&path);

  if (!cs_gui_get_double(path, &result))
    bft_error(__FILE__, __LINE__, 0,
              _("Coordinate %s of the monitoring probe number %i "
                "not found.\nXpath: %s\n"), probe_coord, num_probe, path);

  BFT_FREE(str_num_probe);
  BFT_FREE(path);

  return result;
}

/*-----------------------------------------------------------------------------
 * Retourne le nombre de sous-balises "probe recording" situees dans la balise
 * <variable>
 *----------------------------------------------------------------------------*/

static int cs_gui_variable_number_probes (const char *const variable)
{
  char *path = NULL;
  char *choice = NULL;
  int   nb_probes ;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", variable);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_attribute(&path, "choice");
  choice = cs_gui_get_attribute_value(path);

  if (choice) {
    nb_probes = atoi(choice);
  } else
    nb_probes = -1;

  BFT_FREE(choice);
  BFT_FREE(path);

  return nb_probes;
}

/*-----------------------------------------------------------------------------
 * Return probe number for balice "probe_recording" for variable
 *
 * parameters:
 *   variable   -->  name of variable
 *   num_probe  --> number of balise "probe_recording"
 *----------------------------------------------------------------------------*/

static int cs_gui_variable_probe_name (const char *const variable,
                                             int         num_probe)
{
  char *path = NULL;
  char *strvalue = NULL;
  int   intvalue;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", variable);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_element_num(&path, "probe_recording", num_probe);
  cs_xpath_add_attribute(&path, "name");

  strvalue = cs_gui_get_attribute_value(path);

  if (strvalue == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s \n"), path);

  intvalue = atoi(strvalue);

  BFT_FREE(strvalue);
  BFT_FREE(path);

  return intvalue;
}

/*-----------------------------------------------------------------------------
 * Return label of variable
 *
 * parameters:
 *   variable   --> name of variable
 *----------------------------------------------------------------------------*/

static char *cs_gui_variable_label (const char *const variable)
{
  char *path = NULL;
  char *label = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", variable);
  cs_xpath_add_attribute(&path, "label");

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return label;
}

/*-----------------------------------------------------------------------------
 * Prise en compte des options de post-traitement pour les variables
 *  (Vitesse, Pression...)
 *    le tableau "globale" est construit dans CSENSO
 *          globale[num_probe] = num_saturne_probe
 *      ou num_saturne_probe est le numero de la sonde dans le code
 *      num_probe et num_saturne_probe different lorsque des sondes
 *                                    sont desactivees dans le fichier XML
 *----------------------------------------------------------------------------*/

static void cs_gui_thermophysical_post(const char *const variable,
                                       const int         ipp,
                                             int  *const ihisvr,
                                             int  *const ilisvr,
                                             int  *const ichrvr,
                                       const int  *const nvppmx)
{
  int   nb_probes;
  int   iprob;
  char *varname = NULL;
  int   num_probe;

  if (ipp == 1) return;

  cs_gui_variable_attribute(variable,
                            "postprocessing_recording",
                            &ichrvr[ipp-1]);

  cs_gui_variable_attribute(variable,
                            "listing_printing",
                            &ilisvr[ipp-1]);

  nb_probes = cs_gui_variable_number_probes(variable);

  ihisvr[0 + (ipp - 1)] = nb_probes;

  if (nb_probes > 0) {
    for (iprob =0; iprob < nb_probes; iprob++) {
      num_probe = cs_gui_variable_probe_name(variable, iprob+1);
      ihisvr[(iprob+1)*(*nvppmx) + (ipp - 1)] = num_probe;
    }
  }

  varname = cs_gui_variable_label(variable);
  _gui_copy_varname(varname, ipp);

  BFT_FREE(varname);
}

/*-----------------------------------------------------------------------------
 * Nombre de sous-balise "probe_recording" pour les scalaires utilisateurs
 *----------------------------------------------------------------------------*/

static int cs_gui_scalar_number_probes(const int scalar_num)
{
  char *path = NULL;
  char *choice = NULL;
  int   nb_probes ;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "scalar", scalar_num);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_attribute(&path, "choice");
  choice = cs_gui_get_attribute_value(path);

  if (choice) {
    nb_probes = atoi(choice);
    BFT_FREE(choice);
  } else
    nb_probes = -1;

  BFT_FREE(path);

  return nb_probes;
}

/*-----------------------------------------------------------------------------
 * Return probe number for number of balise  "probe_recording"
 *
 * parameters:
 *   scalar_num  --> number of scalar
 *   num_probe   --> number of balise "probe_recording"
 *----------------------------------------------------------------------------*/

static int cs_gui_scalar_probe_name(const int scalar_num,
                                    const int num_probe)
{
  char *path = NULL;
  char *strvalue = NULL;
  int   value;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "additional_scalars");
  cs_xpath_add_element_num(&path, "scalar", scalar_num);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_element_num(&path, "probe_recording", num_probe);
  cs_xpath_add_attribute(&path, "name");

  strvalue = cs_gui_get_attribute_value(path);

  if (strvalue == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  value = atoi(strvalue);

  BFT_FREE(path);
  BFT_FREE(strvalue);

  return value;
}

/*-----------------------------------------------------------------------------
 * Prise en compte des options de post-traitement pour les scalaires
 *    utilisateurs, le tableau "globale" est construit dans CSENSO
 *          globale[num_probe] = num_saturne_probe
 *                 ou num_saturne_probe est le numero de la sonde dans le code
 *                 num_probe et num_saturne_probe different lorsque des sondes
 *                 sont desactivees dans le fichier XML
 *----------------------------------------------------------------------------*/

static void cs_gui_scalar_post(const  int        num_sca,
                                      int *const ihisvr,
                                      int *const ilisvr,
                                      int *const ichrvr,
                               const  int *const ipprtp,
                               const  int *const isca,
                               const  int *const nvppmx)
{
  int ipp;
  int nb_probes;
  int iprob;
  int num_probe;

  assert(vars != NULL);

  ipp = ipprtp[isca[num_sca] -1 ];

  if (ipp == 1) return;

  /* frequence des sorties ensight */
  cs_gui_scalar_attribute(vars->label[num_sca],
                          "postprocessing_recording",
                          &ichrvr[ipp - 1]);

  /* frequence des sorties listing */
  cs_gui_scalar_attribute(vars->label[num_sca],
                          "listing_printing",
                          &ilisvr[ipp - 1]);

  /* sondes actives */
  nb_probes = cs_gui_scalar_number_probes(num_sca+1);
  /*ihisvr[0][ipp - 1] = nb_probes;*/
  ihisvr[0 + (ipp - 1)] = nb_probes;

  if (nb_probes > 0) {
    for (iprob=0; iprob < nb_probes; iprob++){
      num_probe = cs_gui_scalar_probe_name(num_sca+1, iprob+1);
      ihisvr[(iprob+1)*(*nvppmx) + (ipp - 1)] = num_probe;
    }
  }

  _gui_copy_varname(vars->label[num_sca], ipp);
}

/*-----------------------------------------------------------------------------
 * Return number of sub balises "probe_recording" for model scalars
 *
 * parameters:
 *   model      -->  Type of model
 *   name       -->  scalar name
 *----------------------------------------------------------------------------*/

static int cs_gui_model_scalar_number_probes(const char* const model,
                                             const char *const name)
{
  char *path = NULL;
  char *choice = NULL;
  int   nb_probes ;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_attribute(&path, "choice");
  choice = cs_gui_get_attribute_value(path);

  if (choice) {
    nb_probes = atoi(choice);
    BFT_FREE(choice);
  } else
    nb_probes = -1;

  BFT_FREE(path);

  return nb_probes;
}

/*-----------------------------------------------------------------------------
 * Return probe number for sub balise "probe_recording" for model scalar
 *
 * parameters:
 *   model      --> type of model
 *   name       --> scalar name
 *   num_probe  --> number of balise "probe_recording"
 *----------------------------------------------------------------------------*/

static int cs_gui_model_scalar_probe_name (const char *const model,
                                           const char *const name,
                                           const int         num_probe)
{
  char *path = NULL;
  char *strvalue = NULL;
  int   value;

/* Construction de la requete */
  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_element_num(&path, "probe_recording", num_probe);
  cs_xpath_add_attribute(&path, "name");

/* Evaluation de la requete */
  strvalue = cs_gui_get_attribute_value(path);

  if (strvalue == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  value = atoi(strvalue);

  BFT_FREE(path);
  BFT_FREE(strvalue);

  return value;
}

/*-----------------------------------------------------------------------------
 * Prise en compte des options de post-traitement pour
 *   les scalaires thermiques et model
 *    le tableau "globale" est construit dans CSENSO
 *          globale[num_probe] = num_saturne_probe
 *                 ou num_saturne_probe est le numero de la sonde dans le code
 *                 num_probe et num_saturne_probe different lorsque des sondes
 *                 sont desactivees dans le fichier XML
 *----------------------------------------------------------------------------*/


static void cs_gui_model_scalar_post(const char  *const model,
                                     const int          num_sca,
                                           int   *const ihisvr,
                                           int   *const ilisvr,
                                           int   *const ichrvr,
                                     const int   *const ipprtp,
                                     const int   *const isca,
                                     const int   *const nvppmx)
{
  int ipp;
  int nb_probes;
  int iprob;
  int num_probe;

  assert(vars != NULL);

  ipp = ipprtp[isca[num_sca] -1];

  if (ipp == 1) return;

  /* frequence des sorties ensight */
  cs_gui_model_scalar_output_status(model, vars->label[num_sca],
                                    "postprocessing_recording",
                                    &ichrvr[ipp - 1]);

  /* frequence des sorties listing */
  cs_gui_model_scalar_output_status(model, vars->label[num_sca],
                                    "listing_printing",
                                    &ilisvr[ipp - 1]);

  /* sondes actives */
  nb_probes = cs_gui_model_scalar_number_probes(model, vars->label[num_sca]);

  /* ihisvr[0][ipp - 1] = nb_probes; */
  ihisvr[0 + (ipp - 1)] = nb_probes;

  if (nb_probes > 0) {
    for (iprob =0; iprob < nb_probes; iprob++) {
      num_probe = cs_gui_model_scalar_probe_name(model, vars->label[num_sca], iprob+1);
      ihisvr[(iprob+1)*(*nvppmx) + (ipp - 1)] = num_probe;
    }
  }

  _gui_copy_varname(vars->label[num_sca], ipp);
}

/*-----------------------------------------------------------------------------
 * Return number of sub balises  "probe_recording" for propeety of model scalar
 *
 * parameters:
 *   model    -->  Type of model
 *   num_sca  -->  scalar number
 *----------------------------------------------------------------------------*/

static int cs_gui_model_property_number_probes(const char *const model,
                                               const char *const name)
{
  char *path = NULL;
  char *choice = NULL;
  int   nb_probes ;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_attribute(&path, "choice");
  choice = cs_gui_get_attribute_value(path);

  if (choice) {
    nb_probes = atoi(choice);
    BFT_FREE(choice);
  } else
    nb_probes = -1;

  BFT_FREE(path);

  return nb_probes;
}

/*-----------------------------------------------------------------------------
 * Return probe number for sub balise "probe_recording" for physical model's
 * property
 *
 * parameters:
 *   model      --> type of model
 *   num_prop   --> number of property
 *   num_probe  --> number of balise "probe_recording"
 *----------------------------------------------------------------------------*/

static int cs_gui_model_property_probe_name(const char *const model,
                                            const char *const name,
                                            const int   num_probe)
{
  char *path = NULL;
  char *strvalue = NULL;
  int   value;

/* Construction de la requete */
  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_element_num(&path, "probe_recording", num_probe);
  cs_xpath_add_attribute(&path, "name");

/* Evaluation de la requete */
  strvalue = cs_gui_get_attribute_value(path);

  if (strvalue == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  value = atoi(strvalue);

  BFT_FREE(path);
  BFT_FREE(strvalue);

  return value;
}

/*-----------------------------------------------------------------------------
 * Return the label model's property
 *
 * parameters:
 *   model             -->  modele
 *   num_prop          <--  property's number
 *----------------------------------------------------------------------------*/

static char *cs_gui_get_model_property_label(const char *const model,
                                             const char *const name)
{
  char *path = NULL;
  char *label_name = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_attribute(&path, "label");

/* Evaluation de la requete */
  label_name = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return label_name;
}

/*-----------------------------------------------------------------------------
 * Return status of the property for physical model
 *
 * parameters:
 *   model          --> type of model
 *   num_pro        --> property name
 *   value_type     --> type of value (listing_printing, postprocessing ..)
 *   keyword       <--  value for the Fortran array
 *----------------------------------------------------------------------------*/

static void cs_gui_model_property_output_status (const char *const model,
                                                 const char *const name,
                                                 const char *const value_type,
                                                       int  *const keyword)
{
  char *path = NULL;
  int   result;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_element(&path, value_type);
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
    *keyword = result;
  else
    *keyword = 1;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Prise en compte des options de post-traitement pour
 *   les scalaires thermiques et model
 *    le tableau "globale" est construit dans CSENSO
 *          globale[num_probe] = num_saturne_probe
 *                 ou num_saturne_probe est le numero de la sonde dans le code
 *                 num_probe et num_saturne_probe different lorsque des sondes
 *                 sont desactivees dans le fichier XML
 *----------------------------------------------------------------------------*/


static void cs_gui_model_property_post (const char  *const model,
                                        const int          num_prop,
                                              int   *const ihisvr,
                                              int   *const ilisvr,
                                              int   *const ichrvr,
                                        const int   *const ipppro,
                                        const int   *const ipproc,
                                        const int   *const nvppmx)
{
  int ipp;
  int nb_probes;
  int iprob;
  int num_probe;
  char *varname = NULL;

  assert(vars != NULL);

  ipp = vars->properties_ipp[num_prop];

  if (ipp == 1) return;

  /* frequence des sorties ensight */
  cs_gui_model_property_output_status(model,
                                      vars->properties_name[num_prop],
                                      "postprocessing_recording",
                                      &ichrvr[ipp - 1]);

  /* frequence des sorties listing */
  cs_gui_model_property_output_status(model,
                                      vars->properties_name[num_prop],
                                      "listing_printing",
                                      &ilisvr[ipp - 1]);


  /* sondes actives */
  nb_probes = cs_gui_model_property_number_probes(model,
                                                  vars->properties_name[num_prop]);

  ihisvr[0 + (ipp - 1)] = nb_probes;

  if (nb_probes > 0) {
    for (iprob=0; iprob<nb_probes; iprob++){
      num_probe = cs_gui_model_property_probe_name(model,
                                                   vars->properties_name[num_prop],
                                                   iprob+1);
      ihisvr[(iprob+1)*(*nvppmx) + (ipp - 1)] = num_probe;
    }
  }

  /* prise en compte du label */

  varname = cs_gui_get_model_property_label(model, vars->properties_name[num_prop]);
  _gui_copy_varname(varname, ipp);

  BFT_FREE(varname);
}

/*-----------------------------------------------------------------------------
 * Return number of probes for property
 *
 * parameters:
 *   property_name  --> name of property
 *----------------------------------------------------------------------------*/

static int cs_gui_properties_number_probes(const char *const property_name)
{
  char *path = NULL;
  char *choice = NULL;
  int   nb_probes ;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", property_name);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_attribute(&path, "choice");
  choice = cs_gui_get_attribute_value(path);

  if (choice) {
    nb_probes = atoi(choice);
    BFT_FREE(choice);
  } else
    nb_probes = -1;

  BFT_FREE(path);

  return nb_probes;
}

/*-----------------------------------------------------------------------------
 * Return probe number for sub balise "probe_recording" for properties
 *
 * parameters:
 *   property_name   --> name of property
 *   num_probe       --> number of balise "probe_recording"
 *----------------------------------------------------------------------------*/

static int cs_gui_properties_probe_name(const char *const property_name,
                                        const int         num_probe)
{
  char *path = NULL;
  char *strvalue = NULL;
  int   value;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", property_name);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_element_num(&path, "probe_recording", num_probe);
  cs_xpath_add_attribute(&path, "name");

  strvalue = cs_gui_get_attribute_value(path);

  if (strvalue == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  value = atoi(strvalue);

  BFT_FREE(path);
  BFT_FREE(strvalue);

  return value;
}

/*-----------------------------------------------------------------------------
 * Return the label attribute of a property markup.
 *
 * parameters:
 *   property_name        -->  name of the property
 *----------------------------------------------------------------------------*/

static char *cs_gui_properties_label(const char *const property_name)
{
  char *path = NULL;
  char *label = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", property_name);
  cs_xpath_add_attribute(&path, "label");

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return label;
}

/*-----------------------------------------------------------------------------
 * Return status of thr property markup
 *
 * parameters:
 *   property_name  --> name of property
 *   value_type     --> type of balise (listing_printing, postprocessing ..)
 *   keyword        <-- number of balise "probe_recording"
 *----------------------------------------------------------------------------*/

static void cs_gui_properties_status(const char *const property_name,
                                     const char *const value_type,
                                     int        *const keyword)
{
  char *path = NULL;
  int   result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "property");
  cs_xpath_add_test_attribute(&path, "name", property_name);
  cs_xpath_add_element(&path, value_type);
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
    *keyword = result;
  else
    *keyword = 1;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Prise en compte des options de post-traitement pour les proprietes physiques
 *    le tableau "globale" est construit dans CSENSO
 *          globale[num_probe] = num_saturne_probe
 *                 ou num_saturne_probe est le numero de la sonde dans le code
 *                 num_probe et num_saturne_probe different lorsque des sondes
 *                 sont desactivees dans le fichier XML
 *----------------------------------------------------------------------------*/


static void cs_gui_properties_post(const char *const property_name,
                                   const int         ipp,
                                         int  *const ichrvr,
                                         int  *const ilisvr,
                                         int  *const ihisvr,
                                   const int  *const nvppmx)
{
  int nb_probes;
  int iprob;
  char *varname = NULL;
  int num_probe;

  if (ipp == 1) return;

  varname = cs_gui_properties_label(property_name);
  if (varname == NULL) return;

  _gui_copy_varname(varname, ipp);
  BFT_FREE(varname);

  cs_gui_properties_status(property_name,
                           "postprocessing_recording",
                           &ichrvr[ipp - 1]);

  cs_gui_properties_status(property_name,
                           "listing_printing",
                           &ilisvr[ipp - 1]);

  nb_probes = cs_gui_properties_number_probes(property_name);

  ihisvr[0 + (ipp - 1)] = nb_probes;

  if (nb_probes > 0) {
    for (iprob =0; iprob < nb_probes; iprob++){
      num_probe = cs_gui_properties_probe_name(property_name,
                                               iprob+1);

      ihisvr[(iprob+1)*(*nvppmx) + (ipp - 1)] = num_probe;
    }
  }

}

/*-----------------------------------------------------------------------------
 * Return number of probes for time average of property
 *
 * parameters:
 *   property_name -->  label of property
 *----------------------------------------------------------------------------*/

static int cs_gui_time_average_number_probes(const char *const property_name)
{
  char *path = NULL;
  char *choice = NULL;
  int   nb_probes ;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "time_average");
  cs_xpath_add_test_attribute(&path, "label", property_name);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_attribute(&path, "choice");
  choice = cs_gui_get_attribute_value(path);

  if (choice) {
    nb_probes = atoi(choice);
    BFT_FREE(choice);
  } else
    nb_probes = -1;

  BFT_FREE(path);

  return nb_probes;
}

/*-----------------------------------------------------------------------------
 * Return probe number for sub balise "probe_recording" for time average of
 * properties
 *
 * parameters:
 *   property_name    --> label of property
 *   num_probe        --> number of balise "probe_recording"
 *----------------------------------------------------------------------------*/

static int cs_gui_time_average_probe_name(const char *const property_name,
                                          const int         num_probe)
{
  char *path = NULL;
  char *strvalue = NULL;
  int   value;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "time_average");
  cs_xpath_add_test_attribute(&path, "label", property_name);
  cs_xpath_add_element(&path, "probes");
  cs_xpath_add_element_num(&path, "probe_recording", num_probe);
  cs_xpath_add_attribute(&path, "name");

  strvalue = cs_gui_get_attribute_value(path);

  if (strvalue == NULL)
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  value = atoi(strvalue);

  BFT_FREE(path);
  BFT_FREE(strvalue);

  return value;
}

/*-----------------------------------------------------------------------------
 * Return status of time average markup
 *
 * parameters:
 *   property_name  --> label of property
 *   value_type     --> type of balise (listing_printing, postprocessing ..)
 *   keyword        <-- number of balise "probe_recording"
 *----------------------------------------------------------------------------*/


static void cs_gui_time_average_status(const char *const property_name,
                                       const char *const value_type,
                                             int  *const keyword)
{
  char *path = NULL;
  int   result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "time_average");
  cs_xpath_add_test_attribute(&path, "label", property_name);
  cs_xpath_add_element(&path, value_type);
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result))
    *keyword = result;
  else
    *keyword = 1;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Prise en compte des options de post-traitement pour les moyennes temporelles
 *    le tableau "globale" est construit dans CSENSO
 *          globale[num_probe] = num_saturne_probe
 *                 ou num_saturne_probe est le numero de la sonde dans le code
 *                 num_probe et num_saturne_probe different lorsque des sondes
 *                 sont desactivees dans le fichier XML
 *----------------------------------------------------------------------------*/

static void cs_gui_time_average_post (const char *const property_name,
                               const int         ipp,
                                     int  *const ichrvr,
                                     int  *const ilisvr,
                                     int  *const ihisvr,
                               const int  *const nvppmx)
{
  int nb_probes;
  int iprob;
  int num_probe;

  if (ipp == 1) return;

  cs_gui_time_average_status(property_name,
                             "postprocessing_recording",
                             &ichrvr[ipp - 1]);

  cs_gui_time_average_status(property_name,
                             "listing_printing",
                             &ilisvr[ipp - 1]);

  nb_probes = cs_gui_time_average_number_probes(property_name);

  ihisvr[0 + (ipp - 1)] = nb_probes;

  if (nb_probes > 0) {
    for (iprob =0; iprob < nb_probes; iprob++){
      num_probe = cs_gui_time_average_probe_name(property_name,
                                               iprob+1);

      ihisvr[(iprob+1)*(*nvppmx) + (ipp - 1)] = num_probe;
    }
  }
  _gui_copy_varname(property_name, ipp);

}

/*-----------------------------------------------------------------------------
 * Return the label attribute of scalars.
 *
 * parameters:
 *   markup               -->  parent markup of the scalar
 *   scalar_num          <--   number of the searching scalar
 *----------------------------------------------------------------------------*/

static char *cs_gui_scalar_label(const char *const markup,
                                 const int         scalar_num)
{
  char *path = NULL;
  char *strvalue = NULL;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, markup);
  cs_xpath_add_element_num(&path, "scalar", scalar_num);
  cs_xpath_add_attribute(&path, "label");

  strvalue = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return strvalue;
}

/*===============================
 * FOR VOLUMICS ZONES
 *==========================*/

/*-----------------------------------------------------------------------------
 * Return  the number of zones of initialization
 *----------------------------------------------------------------------------*/

static int cs_gui_volumic_zones_number(void)
{
  int zones = 0;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3,
                        "solution_domain",
                        "volumic_conditions",
                        "zone");
  zones = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return zones;
}

/*-----------------------------------------------------------------------------
 * Return the name of the volumic zone
 *
 * parameters:
 *   ith_zone        -->  id of volumic zone
 *----------------------------------------------------------------------------*/

static char *cs_gui_volumic_zone_name(const int ith_zone)
{
  char *path = NULL;
  char *name = NULL;

  /* 1) get the name of the ith initialization zone */
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "solution_domain", "volumic_conditions");
  cs_xpath_add_element_num(&path, "zone", ith_zone);
  cs_xpath_add_attribute(&path, "name");

  name = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return name;
}

/*-----------------------------------------------------------------------------
 * Return the localisation for the volumic zone named name
 *
 * parameters:
 *   name        -->  name of volumic zone
 *----------------------------------------------------------------------------*/

static char *cs_gui_volumic_zone_localization(const char *const name)
{
  char *path = NULL;
  char *description = NULL;

  /* 2) get the description (color and groups) of the ith initialization zone */
  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "solution_domain",
                                  "volumic_conditions",
                                  "zone");
  cs_xpath_add_test_attribute(&path, "name", name);
  cs_xpath_add_function_text(&path);

  description = cs_gui_get_text_value(path);

  BFT_FREE(path);

  return description;
}

/*-----------------------------------------------------------------------------
 * Return the initial value of variable for the volumic zone named name
 *
 * parameters:
 *   variable_name    -->  name of variable
 *   zone_name        -->  name of volumic zone
 *   initial_value    <--  initial value
 *----------------------------------------------------------------------------*/

static void cs_gui_variable_initial_value(const char   *const variable_name,
                                   const char   *const zone_name,
                                         double *const initial_value)
{
  char *path = NULL;
  double result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "variable");
  cs_xpath_add_test_attribute(&path, "name", variable_name);
  cs_xpath_add_element(&path, "initial_value");
  cs_xpath_add_test_attribute(&path, "zone", zone_name);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *initial_value = result;
  else
    *initial_value = 0.0;

  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return the initial value of scalar for the volumic zone named name
 *
 * parameters:
 *   parent           -->  name of balise parent for the scalar
 *   label            -->  label of scalar
 *   zone_name        -->  name of volumic zone
 *   initial_value    <--  initial value
 *----------------------------------------------------------------------------*/

static void cs_gui_scalar_initial_value(const char   *const parent,
                                 const char   *const label,
                                 const char   *const zone_name,
                                       double *const initial_value)
{
  char *path = NULL;
  char *scalar_name = NULL;
  double result;

  path = cs_xpath_short_path();
  cs_xpath_add_elements(&path, 2, parent, "scalar");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, "initial_value");
  cs_xpath_add_test_attribute(&path, "zone", zone_name);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result))
    *initial_value = result;
  else
    *initial_value = 0.0;

  BFT_FREE(scalar_name);
  BFT_FREE(path);
}

/*===============================
 * FOR BOUNDARIES CONDITIONS
 *===============================*/

/*-----------------------------------------------------------------------------
 * Return the choice for the scalar of boundary condition type
 *
 * parameters:
 *   nature      -->  nature of boundary condition (inlet, wall, symmetry ..)
 *   label       -->  label of boundary condition
 *   var_sca     -->  name of variable(velocity_pressure, turbulence ...)
 *----------------------------------------------------------------------------*/

static char *cs_gui_boundary_choice(const char *const nature,
                                    const char *const label,
                                    const char *const var_sca)
{
  char *path = NULL;
  char *choice = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, var_sca);
  cs_xpath_add_attribute(&path, "choice");

  choice = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return choice;
}

/*-----------------------------------------------------------------------------
 * Put value of dirichlet for variable of velocity_pressure input boundaries.
 *
 * parameters:
 *   nature      -->  nature of boundary condition (inlet, wall, symmetry ..)
 *   label       -->  label of boundary condition
 *   izone       -->  number of zone
 *   ivar        -->  number of variable
 *----------------------------------------------------------------------------*/

static void cs_gui_boundary_dirichlet(const char *const nature,
                                      const char *const label,
                                      const int         izone,
                                      const int         ivar)
{
  char *path = NULL;
  double result = 0.0;

  assert(vars != NULL);

  path = cs_xpath_init_path();

  cs_xpath_add_element(&path, "boundary_conditions");
  cs_xpath_add_element(&path, nature);
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, "velocity_pressure");

  if (cs_gui_strcmp(nature, "wall")) {
    cs_xpath_add_test_attribute(&path, "choice", "on");
  } else if (cs_gui_strcmp(nature, "inlet")) {
    cs_xpath_add_test_attribute(&path, "choice", "dirichlet");
  } else {
    bft_error(__FILE__, __LINE__, 0,
              _("Unknown conditions type in this context: %s.\nXpath: %s\n"),
              nature, path);
  }

  cs_xpath_add_element(&path, "dirichlet");
  cs_xpath_add_test_attribute(&path, "name", vars->name[ivar]);
  cs_xpath_add_function_text(&path);

  if (cs_gui_get_double(path, &result)) {
    boundaries->type_code[vars->rtp[ivar]][izone] = DIRICHLET;
    boundaries->values[vars->rtp[ivar]][izone].val1 = result;
  }
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Put value of roughness for wall
 *
 * parameters:
 *   label       -->  label of boundary condition
 *   izone       -->  number of zone
 *----------------------------------------------------------------------------*/

static void cs_gui_boundary_rough(const char *const label,
                                  const int         izone)
{
  char *path = NULL;
  double result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "wall");
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_element(&path, "rough_wall");

  if (cs_gui_get_double(path, &result)) boundaries->rough[izone] = result;
  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Put values of inlet flow'parameters input boundaries.
 *
 * parameters:
 *   label       -->  label of boundary condition
 *   izone       -->  number of zone
 *----------------------------------------------------------------------------*/

static void cs_gui_boundary_flow(const char *const label,
                                           double *qimp,
                                           double *timp)
{
  char  *path1 = NULL;
  char  *path2 = NULL;
  double result;

  path1 = cs_xpath_init_path();
  cs_xpath_add_elements(&path1, 2, "boundary_conditions", "inlet");
  cs_xpath_add_test_attribute(&path1, "label", label);
  cs_xpath_add_element(&path1, "velocity_pressure");

  BFT_MALLOC(path2, strlen(path1)+1, char);
  strcpy(path2, path1);

  /* flow rate */

  cs_xpath_add_element(&path1, "flow1");
  cs_xpath_add_function_text(&path1);

  if (cs_gui_get_double(path1, &result)){
    *qimp = result;
  }
  BFT_FREE(path1);

  /* temperature */

  cs_xpath_add_element(&path2, "temperature");
  cs_xpath_add_function_text(&path2);

  if (cs_gui_get_double(path2, &result)) {
    *timp = result;
  }

  BFT_FREE(path2);
}

/*-----------------------------------------------------------------------------
 * Put values of inlet turbulence'parameters input boundaries.
 *
 * parameters:
 *   choice      -->  type of choice to calculate turbulence
 *   izone       -->  number of zone
 *----------------------------------------------------------------------------*/

static void cs_gui_boundary_turbulence(const char *const choice,
                                const  int        izone)
{
  char *path1 = NULL;
  char *path2 = NULL;
  double result;

  if (cs_gui_strcmp(choice, "hydraulic_diameter")) {
    boundaries->icalke[izone] = 1  ;
  } else if(cs_gui_strcmp(choice, "turbulent_intensity")) {
    boundaries->icalke[izone] = 2  ;
  } else {
    return;
  }

  path1 = cs_xpath_init_path();
  cs_xpath_add_elements(&path1, 2, "boundary_conditions", "inlet");
  cs_xpath_add_test_attribute(&path1, "label", boundaries->label[izone]);
  cs_xpath_add_element(&path1, "turbulence");

  BFT_MALLOC(path2, strlen(path1) + 1, char);
  strcpy(path2, path1);

  cs_xpath_add_element(&path1, "hydraulic_diameter");
  cs_xpath_add_function_text(&path1);

  if (cs_gui_get_double(path1, &result)) {
    boundaries->dh[izone] = result;
  }
  BFT_FREE(path1);

  if(cs_gui_strcmp(choice, "turbulent_intensity")) {

    cs_xpath_add_element(&path2, "turbulent_intensity");
    cs_xpath_add_function_text(&path2);

    if (cs_gui_get_double(path2, &result)) {
      boundaries->xintur[izone] = result * 0.01;
    }
  }
  BFT_FREE(path2);
}

/*-----------------------------------------------------------------------------
 * Put scalar's values input boundaries.
 *
 * parameters:
 *   nature      -->  nature of boundary condition
 *   izone       -->  number of zone
 *   nsca        -->  number of user scalar
 *----------------------------------------------------------------------------*/

static void cs_gui_boundary_value_scalar(const char *const nature,
                                         const int         izone,
                                         const int         nsca)
{
  int numvar;
  char *path = NULL;
  char *path_commun = NULL;
  char *path2 = NULL;
  char *choice = NULL;
  double result;

  assert(vars != NULL);

  numvar  = vars->nvar - vars->nscaus - vars->nscapp;
  numvar = numvar + nsca;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", nature);
  cs_xpath_add_test_attribute(&path, "label", boundaries->label[izone]);
  cs_xpath_add_element(&path, "scalar");
  cs_xpath_add_test_attribute(&path, "label", vars->label[nsca]);

  BFT_MALLOC(path_commun, strlen(path)+1, char);
  strcpy(path_commun, path);

  BFT_MALLOC(path2, strlen(path)+1, char);
  strcpy(path2, path);

  cs_xpath_add_attribute(&path_commun, "choice");
  choice = cs_gui_get_attribute_value(path_commun);

  if (choice != NULL) {

    if (cs_gui_strcmp(choice, "dirichlet") || cs_gui_strcmp(choice, "exchange_coefficient") || cs_gui_strcmp(choice, "wall_function")) {
      cs_xpath_add_element(&path, "dirichlet");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_double(path, &result)) {
        if (cs_gui_strcmp(choice, "wall_function")) {
          boundaries->type_code[vars->rtp[numvar]][izone] = WALL_FUNCTION;
        } else {
          boundaries->type_code[vars->rtp[numvar]][izone] = DIRICHLET;
        }
        boundaries->values[vars->rtp[numvar]][izone].val1 = result;
      }

    } else if(cs_gui_strcmp(choice, "neumann")) {
      cs_xpath_add_element(&path, "neumann");
      cs_xpath_add_function_text(&path);
      if (cs_gui_get_double(path, &result)) {
        boundaries->type_code[vars->rtp[numvar]][izone] = NEUMANN;
        boundaries->values[vars->rtp[numvar]][izone].val3 = result;
      }
    }

    if (cs_gui_strcmp(choice, "exchange_coefficient")) {
      cs_xpath_add_element(&path2, "exchange_coefficient");
      cs_xpath_add_function_text(&path2);
      if (cs_gui_get_double(path2, &result)) {
        boundaries->type_code[vars->rtp[numvar]][izone] = COEF_ECHANGE;
        boundaries->values[vars->rtp[numvar]][izone].val2 = result;
      }
    }

    BFT_FREE(choice);
  }

  BFT_FREE(path);
  BFT_FREE(path2);
  BFT_FREE(path_commun);
}

/*-----------------------------------------------------------------------------
 * Put coal's values input boundaries.
 *
 * parameters:
 *   izone       -->  number of zone
 *   ncharb      -->
 *   nclpch      -->
 *----------------------------------------------------------------------------*/

static void cs_gui_coal_boundary_coalflow(const int         izone,
                                   const int  *const ncharb,
                                   const int  *const nclpch)
{
  int    icharb;
  int    iratio;
  char  *path1 = NULL;
  char  *path2 = NULL;
  char  *path3 = NULL;
  char  *path4 = NULL;
  char  *path5 = NULL;
  char  *coalname = NULL;
  char  *classname = NULL;
  double value;

  path1 = cs_xpath_init_path();
  cs_xpath_add_elements(&path1, 2, "boundary_conditions", "inlet");
  cs_xpath_add_test_attribute(&path1, "label", boundaries->label[izone]);
  cs_xpath_add_elements(&path1, 2, "velocity_pressure", "coal_flow");

  BFT_MALLOC(coalname,  4 + 2 + 1, char);
  BFT_MALLOC(classname, 5 + 2 + 1, char);

  for (icharb = 0; icharb < *ncharb; icharb++) {

    BFT_MALLOC(path2, strlen(path1) + 1, char);
    strcpy(path2, path1);

    sprintf(coalname, "%.4s%2.2i", "coal", icharb+1);
    cs_xpath_add_test_attribute(&path2, "name", coalname);

    BFT_MALLOC(path3, strlen(path2) + 1, char);
    strcpy(path3, path2);

    BFT_MALLOC(path4, strlen(path2) + 1, char);
    strcpy(path4, path2);

    /* flow rate */

    cs_xpath_add_element(&path3, "flow1");
    cs_xpath_add_function_text(&path3);
    if (cs_gui_get_double(path3, &value)) {
      boundaries->qimpcp[izone][icharb] = value;
    }

    /* temperature */

    cs_xpath_add_element(&path4, "temperature");
    cs_xpath_add_function_text(&path4);
    if (cs_gui_get_double(path4, &value)) {
      boundaries->timpcp[izone][icharb] = value;
    }

    /* ratio */

    for (iratio=0; iratio < nclpch[icharb]; iratio++) {

      BFT_MALLOC(path5, strlen(path2) + 1, char);
      strcpy(path5, path2);

      cs_xpath_add_element(&path5, "ratio");
      sprintf(classname, "%.5s%2.2i", "class", iratio+1);
      cs_xpath_add_test_attribute(&path5, "name", classname);
      cs_xpath_add_function_text(&path5);

      if (cs_gui_get_double(path5, &value))
        boundaries->distch[izone][icharb][iratio] = value;

      BFT_FREE(path5);

    }
    BFT_FREE(path2);
    BFT_FREE(path3);
    BFT_FREE(path4);
  }

  BFT_FREE(path1);
  BFT_FREE(coalname);
  BFT_FREE(classname);
}

/*===========================================
 * Fonctions pour les physiques particulieres
 *===========================================*/

/*-----------------------------------------------------------------------------
 * Return number of model's properties for particular physical
 *
 * parameters:
 *   model       -->  name of type of model
 *----------------------------------------------------------------------------*/

static int cs_gui_get_number_model_properties(const char* model)
{
  char *path = NULL;
  int   nb;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element(&path, "property");

  nb = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return nb;
}

/*-----------------------------------------------------------------------------
 *  Return value of attribute name for the scalar number for particular physical
 *
 * parameters:
 *   model        -->  name of type of model
 *   scalar_num   -->  scalar number
 *----------------------------------------------------------------------------*/

static char* cs_gui_get_model_scalar_name (const char * const model,
                                           const int          scalar_num)
{
  char * path = NULL;
  char * name_scalar = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element_num(&path, "scalar", scalar_num);
  cs_xpath_add_attribute(&path, "name");

  name_scalar = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return name_scalar;
}

/*-----------------------------------------------------------------------------
 *  Return scalar number for scalar named name for particular physical
 *
 * parameters:
 *   model       -->  name of type of model
 *   name        -->  name of scalar
 *----------------------------------------------------------------------------*/

static int cs_gui_get_model_scalar_number(const  char * const model,
                                          const  char * const name)
{
  int i;
  int nbsca = 0;
  char * nametmp = NULL;
  int numsca;

  nbsca = cs_gui_model_scalar_number(model);
  numsca = 0;

  for (i = 0; i < nbsca; i++){
    nametmp = cs_gui_get_model_scalar_name(model, i+1);
    if (cs_gui_strcmp(name, nametmp)) {
      numsca = i;
      BFT_FREE(nametmp);
      break;
    }
    BFT_FREE(nametmp);
  }

  if (i == nbsca)
    bft_error(__FILE__, __LINE__, 0, _("Invalid scalar name: %s.\n"), name);

  return numsca;
}

/*-----------------------------------------------------------------------------
 * Return value of attribute name for the property number for particular physical
 *
 * parameters:
 *   model        -->  name of type of model
 *   pp_num       -->  number of property
 *----------------------------------------------------------------------------*/

static char* cs_gui_get_model_property_name (const char * const model,
                                             const int          pp_num)
{
  char * path = NULL;
  char * name_pp = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_element(&path, "thermophysical_models");
  cs_xpath_add_element(&path, model);
  cs_xpath_add_element_num(&path, "property", pp_num);
  cs_xpath_add_attribute(&path, "name");

  name_pp = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return name_pp;
}

/*-----------------------------------------------------------------------------
 *  Return property number for particular physical
 *
 * parameters:
 *   model       -->  name of type of model
 *   name        -->  name of property
 *----------------------------------------------------------------------------*/

static int cs_gui_get_model_property_number (const  char * const model,
                                      const  char * const name)
{
  int i;
  int nbpp = 0;
  char * nametmp = NULL;
  int numpp;

  nbpp = cs_gui_get_number_model_properties(model);
  numpp = 0;

  for (i = 0; i < nbpp; i++){
    nametmp = cs_gui_get_model_property_name(model, i + 1);
    if (cs_gui_strcmp(name, nametmp)) {
      numpp = i;
      BFT_FREE(nametmp);
      break;
    }
    BFT_FREE(nametmp);
  }

  if (i == nbpp)
    bft_error(__FILE__, __LINE__, 0, _("Invalid property name: %s.\n"), name);

  return numpp;
}

/*============================
 * Functions for ALE method
 *============================*/


/*-----------------------------------------------------------------------------
 * Return the status of ALE method
 *
 * parameters:
 *   keyword        <--  status of ale balise
 *----------------------------------------------------------------------------*/

static void cs_gui_get_ale_status(int  *const keyword)
{
  char *path = NULL;
  int   result;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", "ale_method");
  cs_xpath_add_attribute(&path, "status");

  if(cs_gui_get_status(path, &result))
    *keyword = result;
  else
    *keyword = 0;


  BFT_FREE(path);
}

/*-----------------------------------------------------------------------------
 * Return the viscosity's type of ALE method
 *
 * parameters:
 *   type        <--  type of viscosity's type
 *----------------------------------------------------------------------------*/

static void cs_gui_get_ale_viscosity_type(int  * type)
{
  char *path = NULL;
  char *buff = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "thermophysical_models", "ale_method", "mesh_viscosity");
  cs_xpath_add_attribute(&path, "type");

  buff = cs_gui_get_attribute_value(path);

  if (cs_gui_strcmp(buff, "orthotrop"))
    *type = 1;
  else if (cs_gui_strcmp(buff, "isotrop"))
    *type = 0;
  else
    bft_error(__FILE__, __LINE__, 0, _("Invalid xpath: %s\n"), path);

  BFT_FREE(path);
  BFT_FREE(buff);
}

/*----------------------------------------------------------------------------
 * Return integer value for calculation of size of user arrays
 *----------------------------------------------------------------------------*/

static int _user_array(const char *const keyword1,
                       const char *const keyword2)
{
  char *path = NULL;
  int value = 0;

  path = cs_xpath_short_path();
  cs_xpath_add_elements(&path, 2, keyword1, keyword2);
  cs_xpath_add_function_text(&path);
  cs_gui_get_int(path, &value);
  BFT_FREE(path);
  return value;
}

/*============================================================================
 * C API public functions
 *============================================================================*/

/*-----------------------------------------------------------------------------
 * Return the name of a thermophysical model.
 *
 * parameter:
 *   model_thermo          -->  thermophysical model
 *----------------------------------------------------------------------------*/


char *cs_gui_get_thermophysical_model(const char *const model_thermo)
{
  char *model = NULL;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "thermophysical_models", model_thermo);
  cs_xpath_add_attribute(&path, "model");

  model = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return model;
}

/*-----------------------------------------------------------------------------
 * Return 1 if a particular physics model is activated. Store in the global
 * structure vars:
 *   vars->model         <= thermophysical model
 *   vars->model_value   <= related model name
 *----------------------------------------------------------------------------*/

int cs_gui_get_activ_thermophysical_model(void)
{
  int isactiv = 0;
  char *value = NULL;

  assert(vars != NULL);

  if (vars->model != NULL && vars->model_value != NULL) {
    isactiv = 1;
    return isactiv;
  }

  value = cs_gui_get_thermophysical_model("pulverized_coal");

  if (!cs_gui_strcmp(value, "off")) {
    BFT_MALLOC(vars->model, strlen("pulverized_coal")+1, char);
    strcpy(vars->model, "pulverized_coal");

    BFT_MALLOC(vars->model_value, strlen(value)+1, char);
    strcpy(vars->model_value, value);

    isactiv = 1;

  } else {
    vars->model = NULL;
    vars->model_value = NULL;
  }

/* modeltmp = cs_gui_get_thermophysical_model(JOULE_EFFECT); */
/* modeltmp = cs_gui_get_thermophysical_model(GAS_COMBUSTION); */
/* modeltmp = cs_gui_get_thermophysical_model(RADIATIVE_TRANSFER); */

  BFT_FREE(value);

  return isactiv;
}

/*-----------------------------------------------------------------------------
 * Return number of boundary regions definition
 *----------------------------------------------------------------------------*/

int cs_gui_boundary_zones_number(void)
{
  int zones = 0;
  char *path = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "boundary_definition");
  cs_xpath_add_all_elements(&path);
  cs_xpath_add_attribute(&path, "label");

  zones = cs_gui_get_nb_element(path);

  BFT_FREE(path);

  return zones;
}

/*-----------------------------------------------------------------------------
 * Return the nature of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

char *cs_gui_boundary_zone_nature(const int ith_zone)
{
  char *path = NULL;
  char *nature = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "boundary_definition");
  cs_xpath_add_element_num(&path, "*", ith_zone);

  nature = cs_gui_get_node_name(path);

  BFT_FREE(path);

  return nature;
}


/*-----------------------------------------------------------------------------
 * Return the label of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

char *cs_gui_boundary_zone_label(const int ith_zone)
{
  char *path = NULL;
  char *label = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "boundary_definition");
  cs_xpath_add_element_num(&path, "*", ith_zone);
  cs_xpath_add_attribute(&path, "label");

  label = cs_gui_get_attribute_value(path);

  BFT_FREE(path);

  return label;
}


/*-----------------------------------------------------------------------------
 * Return the zone number of boundary condition for the given zone
 *----------------------------------------------------------------------------*/

int cs_gui_boundary_zone_number(const int ith_zone)
{
  char *path = NULL;
  char *czone = NULL;
  int zone;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 2, "boundary_conditions", "boundary_definition");
  cs_xpath_add_element_num(&path, "*", ith_zone);
  cs_xpath_add_attribute(&path, "zone");

  czone = cs_gui_get_attribute_value(path);
  zone = atoi(czone);

  BFT_FREE(path);
  BFT_FREE(czone);

  return zone;
}

/*-----------------------------------------------------------------------------
 * Return the description of a boundary zone
 *
 * parameters:
 *   nature                -->  nature of boundary zone (inlet, wall,...)
 *   label                   -->  label of boundary zone
 *----------------------------------------------------------------------------*/

char *cs_gui_boundary_zone_localization(const char *const nature,
                                        const char *const label)
{
  char *path = NULL;
  char *localization = NULL;

  path = cs_xpath_init_path();
  cs_xpath_add_elements(&path, 3, "boundary_conditions",
                                  "boundary_definition",
                                   nature);
  cs_xpath_add_test_attribute(&path, "label", label);
  cs_xpath_add_function_text(&path);

  localization = cs_gui_get_text_value(path);

  BFT_FREE(path);

  return localization;
}

/*============================================================================
 * Fortran API public functions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Turbulence model.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTURB (ITURB, IDEUCH, IGRAKE, IGRAKI, XLOMLG)
 * *****************
 *
 * INTEGER          ITURB   <--   turbulence model
 * INTEGER          IDEUCH  <--   wall law treatment
 * INTEGER          IGRAKE  <--   k-eps gravity effects
 * INTEGER          IGRAKI  <--   Rij-eps gravity effects
 * DOUBLE PRECISION XLOMLG  <--   mixing_length_scale
 *----------------------------------------------------------------------------*/


void CS_PROCF (csturb, CSTURB) (int    *const iturb,
                                int    *const ideuch,
                                int    *const igrake,
                                int    *const igrari,
                                double *const xlomlg)
{
  char *model = NULL;
  int iphas = 0;

  model = cs_gui_get_thermophysical_model("turbulence");
  if (model == NULL) return;

  if (cs_gui_strcmp(model, "off"))
     iturb[iphas] = 0;
  else if (cs_gui_strcmp(model, "mixing_length")){
     iturb[iphas] = 10;
     _option_turbulence_double("mixing_length_scale", &xlomlg[iphas]);
   }
  else if (cs_gui_strcmp(model, "k-epsilon")){
     iturb[iphas] = 20;
     cs_gui_advanced_options_turbulence("scale_model", &ideuch[iphas]);
     cs_gui_advanced_options_turbulence("gravity_terms", &igrake[iphas]);
   }
  else if (cs_gui_strcmp(model, "k-epsilon-PL")){
     iturb[iphas] = 21;
     cs_gui_advanced_options_turbulence("scale_model", &ideuch[iphas]);
     cs_gui_advanced_options_turbulence("gravity_terms", &igrake[iphas]);
   }
  else if (cs_gui_strcmp(model, "Rij-epsilon")){
     iturb[iphas] = 30;
     cs_gui_advanced_options_turbulence("scale_model", &ideuch[iphas]);
     cs_gui_advanced_options_turbulence("gravity_terms", &igrari[iphas]);
   }
  else if (cs_gui_strcmp(model, "Rij-SSG")){
     iturb[iphas] = 31;
     cs_gui_advanced_options_turbulence("scale_model", &ideuch[iphas]);
     cs_gui_advanced_options_turbulence("gravity_terms", &igrari[iphas]);
   }
  else if (cs_gui_strcmp(model, "LES_Smagorinsky")){
     iturb[iphas] = 40;
   }
  else if (cs_gui_strcmp(model, "LES_dynamique")){
     iturb[iphas] = 41;
   }
  else if (cs_gui_strcmp(model, "v2f-phi")){
     iturb[iphas] = 50;
     cs_gui_advanced_options_turbulence("scale_model", &ideuch[iphas]);
     cs_gui_advanced_options_turbulence("gravity_terms", &igrake[iphas]);
   }
  else if (cs_gui_strcmp(model, "k-omega-SST")){
     iturb[iphas] = 60;
     cs_gui_advanced_options_turbulence("scale_model", &ideuch[iphas]);
     cs_gui_advanced_options_turbulence("gravity_terms", &igrake[iphas]);
   }
  else
     bft_error(__FILE__, __LINE__, 0,
               _("Invalid turbulence model: %s.\n"), model);

#if _XML_DEBUG_
  bft_printf(_("==>CSTURB\n"));
  bft_printf(_("--model: %s\n"), model);
  bft_printf(_("--iturb = %i\n"), iturb[iphas]);
  bft_printf(_("--igrake = %i\n"), igrake[iphas]);
  bft_printf(_("--igrari = %i\n"), igrari[iphas]);
  bft_printf(_("--ideuch = %i\n"), ideuch[iphas]);
  bft_printf(_("--xlomlg = %f\n"), xlomlg[iphas]);
#endif

  BFT_FREE(model);
}

/*----------------------------------------------------------------------------
 * Specific heat variable or constant indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSCPVA (ICP)
 * *****************
 *
 * INTEGER          ICP     <--   specific heat variable or constant indicator
 *----------------------------------------------------------------------------*/

void CS_PROCF (cscpva, CSCPVA) (int *const icp)
{
  int choice;
  int iphas = 0;

  if (cs_gui_properties_choice("specific_heat", &choice)) icp[iphas] = choice;

#if _XML_DEBUG_
  bft_printf(_("==>CSCPVA\n"));
  bft_printf(_("--icp = %i\n"), icp[iphas]);
#endif
}

/*----------------------------------------------------------------------------
 * User scalars number.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNSCA (NSCAUS)
 * *****************
 *
 * INTEGER          NSCAUS     <--   user scalars number
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnsca, CSNSCA) (int *const nscaus)
{
  int   i   = 0;
  char *label = NULL;

  *nscaus = cs_gui_get_number_user_scalar();

  if (vars == NULL) {
    BFT_MALLOC(vars, 1, cs_var_t);
    vars->model = NULL;
    vars->model_value = NULL;
    vars->head = NULL;
    vars->type = NULL;
    vars->name = NULL;
    vars->label = NULL;
    vars->rtp = NULL;
    vars->nvar = 0;
    vars->nscaus = 0;
    vars->nscapp = 0;
    vars->nprop = 0;
    vars->nsalpp = 0;
    vars->ntimaver = 0;
    vars->properties_name = NULL;
    vars->properties_ipp = NULL;
  }  else  {
    bft_error(__FILE__, __LINE__, 0,
     _("Trouble with the allocated memory for the global variable 'vars'.\n"));
  }

  vars->nscaus = *nscaus;

  BFT_MALLOC(vars->label, *nscaus, char*);

  for (i=0; i<vars->nscaus; i++) {
    label = cs_gui_scalar_label("additional_scalars", i+1);
    BFT_MALLOC(vars->label[i], strlen(label)+1, char);
    strcpy(vars->label[i], label);
    BFT_FREE(label);
  }

#if _XML_DEBUG_
  bft_printf(_("==>CSNSCA\n"));
  bft_printf(_("--user scalars number: %i\n"), vars->nscaus);
  for (i=0; i<*nscaus; i++)
    bft_printf(_("--label of scalar[%i]: %s\n"), i, vars->label[i]);
#endif
}

/*-----------------------------------------------------------------------------
 * Predefined physics indicator.
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIPPMO
 * *****************
 *
 * INTEGER          IPPMOD <--  predefined physics indicator array
 * INTEGER          ICOD3P  --> diffusion flame en chimie complete rapide
 * INTEGER          ICODEQ  --> diffusion flame en chimie rapide vers l'equilibre
 * INTEGER          ICOEBU  --> Eddy Break Up premixing flame
 * INTEGER          ICOBML  --> Bray - Moss - Libby premixing flame
 * INTEGER          ICOLWC  --> Libby Williams premixing flame
 * INTEGER          ICP3PL  --> Coal combustion. Combustible moyen local
 * INTEGER          ICPL3C  --> Coal combustion coupled with lagrangien approach
 * INTEGER          ICFUEL  --> Fuel combustion
 * INTEGER          IELJOU  --> Joule effect
 * INTEGER          IELARC  --> electrical arc
 * INTEGER          IELION  --> ionique mobility
 * INTEGER          ICOMPF  --> compressible without shock
 * INTEGER          INDJON  --> INDJON=1: a JANAF enthalpy-temperature
 *                              tabulation is used. INDJON=1: users tabulation
 * INTEGER          IEQCO2  --> CO2 massic fraction transport
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uippmo, UIPPMO)(int *const ippmod,
                               int *const icod3p,
                               int *const icodeq,
                               int *const icoebu,
                               int *const icobml,
                               int *const icolwc,
                               int *const icp3pl,
                               int *const icpl3c,
                               int *const icfuel,
                               int *const ieljou,
                               int *const ielarc,
                               int *const ielion,
                               int *const icompf,
                               int *const indjon,
                               int *const ieqco2)
{
  int isactiv = 0;
  int nscapp = 0;

  assert(vars != NULL);

  /* init */
  ippmod[*icod3p - 1] = -1;
  ippmod[*icodeq - 1] = -1;
  ippmod[*icoebu - 1] = -1;
  ippmod[*icobml - 1] = -1;
  ippmod[*icolwc - 1] = -1;
  ippmod[*icp3pl - 1] = -1;
  ippmod[*icpl3c - 1] = -1;
  ippmod[*icfuel - 1] = -1;
  ippmod[*ieljou - 1] = -1;
  ippmod[*ielarc - 1] = -1;
  ippmod[*ielion - 1] = -1;
  ippmod[*icompf - 1] = -1;

  *indjon = 1;
  *ieqco2 = 0;

  /* cherche la physique particuliere active et donne la valeur de
     l'attribut model associe */
  isactiv = cs_gui_get_activ_thermophysical_model();

  if (isactiv) {

    if (cs_gui_strcmp(vars->model, "pulverized_coal")) {

      if (cs_gui_strcmp(vars->model_value, "coal_homo")) {
        ippmod[*icp3pl - 1] = 0;
      } else if (cs_gui_strcmp(vars->model_value, "coal_homo2")) {
        ippmod[*icp3pl - 1] = 1;
      } else if (cs_gui_strcmp(vars->model_value, "coal_lagr")) {
        ippmod[*icpl3c - 1] = 1;
      } else {
        bft_error(__FILE__, __LINE__, 0,
                  _("Invalid coal model: %s.\n"), vars->model_value);
      }
    }
    /* si le model est actif on prend les scalaires physique particuliere */
    nscapp = cs_gui_model_scalar_number(vars->model);
  }

  vars->nscapp = nscapp;

#if _XML_DEBUG_
  bft_printf(_("==>UIPPMO\n"));
  if (isactiv) {
    bft_printf(_("--thermophysical model: %s\n"), vars->model);
    bft_printf(_("--thermophysical value: %s\n"), vars->model_value);
    bft_printf(_("--model scalars number: %i\n"), vars->nscapp);
  }
#endif

}

/*----------------------------------------------------------------------------
 * User scalars which are variance.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISCA (ISCAVR)
 * *****************
 *
 * INTEGER          ISCAVR     <--   user scalars variance array
 *----------------------------------------------------------------------------*/

void CS_PROCF (csisca, CSISCA) (int *const iscavr)
{
  int i;
  int j;
  char *variance = NULL;

  assert(vars != NULL);

  if (vars->nscaus > 0) {

    for (i=0; i < vars->nscaus; i++) {
      variance = cs_gui_scalar_variance(i+1);

      if (variance != NULL) {

        for (j=0 ; j < vars->nscaus ; j++) {

          if (cs_gui_strcmp(variance, vars->label[j])) {

            if ( i == j )
              bft_error(__FILE__, __LINE__, 0,
                 _("Scalar: %i and its variance: %i are the same.\n"),
                 i, j);
            iscavr[i] = j + 1;
          }
        }
        BFT_FREE(variance);
      }
    }

#if _XML_DEBUG_
    bft_printf(_("==>CSISCA\n"));
    for (i = 0 ; i < vars->nscaus ; i++)
      bft_printf(_("--iscavr[%i] = %i \n"), i, iscavr[i]);
#endif
  }
  return;
}

/*----------------------------------------------------------------------------
 * Constant or variable indicator for the user scalar laminar viscosity.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIVIS (ISCAVR, IVISLS, ISCALT, ISCSTH)
 * *****************
 *
 * INTEGER          ISCAVR  <-->  number of the related variance if any
 * INTEGER          IVISLS  <--   indicator for the user scalar viscosity
 * INTEGER          ISCALT  <-->  number of the user thermal scalar if any
 * INTEGER          ISCSTH  <-->  type of the user thermal scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (csivis, CSIVIS) (int *const iscavr,
                                int *const ivisls,
                                int *const iscalt,
                                int *const iscsth)
{
  int iphas = 0;
  int i;
  int choice1, choice2;
  int test1, test2;

  assert(vars != NULL);

  if (vars->nscaus > 0) {

    if (cs_gui_thermal_scalar()) {
      test1 = cs_gui_properties_choice("thermal_conductivity", &choice1);
      test2 = cs_gui_properties_choice("specific_heat", &choice2);

      if (test1 && test2) {
        cs_gui_thermal_scalar_number(iscalt, iscsth);

        if (choice1 || choice2)
          ivisls[iscalt[iphas]-1] = 1;
        else
          ivisls[iscalt[iphas]-1] = 0;
      }
    }

    for (i=0 ; i < vars->nscaus; i++) {
      if (iscavr[i] <= 0 ) {
        if (cs_gui_scalar_properties_choice(i+1,
                                            "diffusion_coefficient",
                                            &choice1))
        if (iscalt[iphas] != i+1) ivisls[i] = choice1;
      }
    }

#if _XML_DEBUG_
    bft_printf(_("==>CSIVIS\n"));
    for (i=0 ; i < vars->nscaus ; i++)
      bft_printf(_("--ivisls[%i] = %i\n"), i, ivisls[i]);
#endif
    }
}

/*----------------------------------------------------------------------------
 * Time passing parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIDTV (IDTVAR)
 * *****************
 *
 * INTEGER          IDTVAR  <--   fixed or variable time step
 *----------------------------------------------------------------------------*/

void CS_PROCF (csidtv, CSIDTV) (int *const idtvar)
{
  double param;
  int steady = 0;

  cs_gui_get_steady_status(&steady);
  if (steady){
    *idtvar = -1;
  }
  else{
    param = (double) *idtvar;
    cs_gui_time_parameters("time_passing", &param);
    *idtvar = (int) param;
  }

#if _XML_DEBUG_
  bft_printf(_("==>CSIDTV\n"));
  bft_printf(_("--idtvar = %i\n"), *idtvar);
#endif
}

/*----------------------------------------------------------------------------
 * Hydrostatic pressure parameter.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIPHY (IPHYDR)
 * *****************
 *
 * INTEGER          IPHYDR  <--   hydrostatic pressure
 *----------------------------------------------------------------------------*/

void CS_PROCF (csiphy, CSIPHY) (int *const iphydr)
{
  char *path = NULL;
  int   result;

  path = cs_xpath_short_path();
  cs_xpath_add_element(&path, "hydrostatic_pressure");
  cs_xpath_add_attribute(&path, "status");

  if (cs_gui_get_status(path, &result)) *iphydr = result;

  BFT_FREE(path);

#if _XML_DEBUG_
  bft_printf(_("==>CSIPHY\n"));
  bft_printf(_("--iphydr = %i\n"), *iphydr);
#endif
}

/*----------------------------------------------------------------------------
 * ALE related keywords
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIALIN
 * *****************
 *
 * INTEGER          IALE    <--  iale method activation
 * INTEGER          NALINF  <--  number of sub iteration of initialization
 *                               of fluid
 * INTEGER          NALIMX  <--  max number of iterations of implicitation of
 *                               the displacement of the structures
 * DOUBLE           EPALIM  <--  realtive precision of implicitation of
 *                               the displacement of the structures
 * INTEGER          IORTVM  <--  type of viscosity of mesh
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (uialin, UIALIN) (int    *const iale,
                                int    *const nalinf,
                                int    *const nalimx,
                                double *const epalim,
                                int    *const iortvm)
{
  double value;

  cs_gui_get_ale_status(iale);

  if (*iale) {
    value =(double) *nalinf;
    cs_gui_iale_parameter("fluid_initialization_sub_iterations", &value);
    *nalinf = (int) value;

    value =(double) *nalimx;
    cs_gui_iale_parameter("max_iterations_implicitation", &value);
    *nalimx = (int) value;

    cs_gui_iale_parameter("implicitation_precision", epalim);

    value =(double) *iortvm;
    cs_gui_iale_parameter("mesh_viscosity", &value);
    *iortvm = (int) value;
  }

#if _XML_DEBUG_
  bft_printf(_("==>UIALIN\n"));
  bft_printf(_("--iale = %i\n"), *iale);
  if (*iale) {
    bft_printf(_("--nalinf = %i\n"), *nalinf);
    bft_printf(_("--nalimx = %i\n"), *nalimx);
    bft_printf(_("--epalim = %g\n"), *epalim);
    bft_printf(_("--iortvm = %i\n"), *iortvm);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Constructs an indirection between an internal index vars->rtp and
 * the fortran array RTP.
 * This function is called after the third call to VARPOS routine.
 * Warning: in vars->rtp, variables are stored first and then
 * scalars (not like VARPOS because of ALE variables).
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSVNUM
 * *****************
 *
 * INTEGER          IPHYDR  <--   hydrostatic pressure
 *----------------------------------------------------------------------------*/

void CS_PROCF (csvnum, CSVNUM) (const int *const nvar,
                                const int *const iu,
                                const int *const iv,
                                const int *const iw,
                                const int *const ipr,
                                const int *const iturb,
                                const int *const ik,
                                const int *const iep,
                                const int *const ir11,
                                const int *const ir22,
                                const int *const ir33,
                                const int *const ir12,
                                const int *const ir13,
                                const int *const ir23,
                                const int *const iomg,
                                const int *const iphi,
                                const int *const ifb,
                                const int *const iale,
                                const int *const iuma,
                                const int *const ivma,
                                const int *const iwma,
                                const int *const isca,
                                const int *const iscapp)
{
  int iphas = 0;
  int n = 0;
  int i, j, k;

  assert(vars != NULL);

  BFT_MALLOC(vars->rtp,  *nvar, int);
  BFT_MALLOC(vars->head, *nvar, char*);
  BFT_MALLOC(vars->type, *nvar, char*);
  BFT_MALLOC(vars->name, *nvar, char*);

  /* Warning!!  vars->nscaus is already fill in CSNSCA */
  /*            vars->label  is already fill in CSNSCA */
  /*            vars->nscapp is already fill in UIPPMO */

  vars->nvar   = *nvar;

  /* 1) pressure and velocity variables */

  k = n;
  vars->rtp[n] = ipr[iphas] -1;
  BFT_MALLOC(vars->name[n], strlen("pressure")+1, char);
  strcpy(vars->name[n++], "pressure");

  vars->rtp[n] = iu[iphas]  -1;
  BFT_MALLOC(vars->name[n], strlen("velocity_U")+1, char);
  strcpy(vars->name[n++], "velocity_U");

  vars->rtp[n] = iv[iphas]  -1;
  BFT_MALLOC(vars->name[n], strlen("velocity_V")+1, char);
  strcpy(vars->name[n++], "velocity_V");

  vars->rtp[n] = iw[iphas]  -1;
  BFT_MALLOC(vars->name[n], strlen("velocity_W")+1, char);
  strcpy(vars->name[n++], "velocity_W");

  for (i=k; i < n; i++) {
    BFT_MALLOC(vars->head[i], strlen("velocity_pressure")+1, char);
    strcpy(vars->head[i], "velocity_pressure");
  }

  /* 2) turbulence variables */

  k = n;

  if (iturb[iphas] == 20 || iturb[iphas] == 21) {

    vars->rtp[n] = ik[iphas]  -1;
    BFT_MALLOC(vars->name[n], strlen("turb_k")+1, char);
    strcpy(vars->name[n++], "turb_k");

    vars->rtp[n] = iep[iphas] -1;
    BFT_MALLOC(vars->name[n], strlen("turb_eps")+1, char);
    strcpy(vars->name[n++], "turb_eps");

  } else if (iturb[iphas] == 30 || iturb[iphas] == 31) {

    vars->rtp[n] = ir11[iphas] -1;
    BFT_MALLOC(vars->name[n], strlen("component_R11")+1, char);
    strcpy(vars->name[n++], "component_R11");

    vars->rtp[n] = ir22[iphas] -1;
    BFT_MALLOC(vars->name[n], strlen("component_R22")+1, char);
    strcpy(vars->name[n++], "component_R22");

    vars->rtp[n] = ir33[iphas] -1;
    BFT_MALLOC(vars->name[n], strlen("component_R33")+1, char);
    strcpy(vars->name[n++], "component_R33");

    vars->rtp[n] = ir12[iphas] -1;
    BFT_MALLOC(vars->name[n], strlen("component_R12")+1, char);
    strcpy(vars->name[n++], "component_R12");

    vars->rtp[n] = ir13[iphas] -1;
    BFT_MALLOC(vars->name[n], strlen("component_R13")+1, char);
    strcpy(vars->name[n++], "component_R13");

    vars->rtp[n] = ir23[iphas] -1;
    BFT_MALLOC(vars->name[n], strlen("component_R23")+1, char);
    strcpy(vars->name[n++], "component_R23");

    vars->rtp[n] = iep[iphas]  -1;
    BFT_MALLOC(vars->name[n], strlen("turb_eps")+1, char);
    strcpy(vars->name[n++], "turb_eps");

  } else if (iturb[iphas] == 50) {

    vars->rtp[n] = ik[iphas]   -1;
    BFT_MALLOC(vars->name[n], strlen("turb_k")+1, char);
    strcpy(vars->name[n++], "turb_k");

    vars->rtp[n] = iep[iphas]  -1;
    BFT_MALLOC(vars->name[n], strlen("turb_eps")+1, char);
    strcpy(vars->name[n++], "turb_eps");

    vars->rtp[n] = iphi[iphas] -1;
    BFT_MALLOC(vars->name[n], strlen("turb_phi")+1, char);
    strcpy(vars->name[n++], "turb_phi");

    vars->rtp[n] = ifb[iphas]  -1;
    BFT_MALLOC(vars->name[n], strlen("turb_fb")+1, char);
    strcpy(vars->name[n++], "turb_fb");

  } else if (iturb[iphas] == 60) {

    vars->rtp[n] = ik[iphas]   -1;
    BFT_MALLOC(vars->name[n], strlen("turb_k")+1, char);
    strcpy(vars->name[n++], "turb_k");

    vars->rtp[n] = iomg[iphas] -1;
    BFT_MALLOC(vars->name[n], strlen("turb_omega")+1, char);
    strcpy(vars->name[n++], "turb_omega");
  }

  for (i=k; i < n; i++) {
    BFT_MALLOC(vars->head[i], strlen("turbulence")+1, char);
    strcpy(vars->head[i], "turbulence");
  }

  /* 3) ALE variables */

  if (*iale) {
    k = n;

    vars->rtp[n] = *iuma -1;
    BFT_MALLOC(vars->name[n], strlen("mesh_velocity_U")+1, char);
    strcpy(vars->name[n++], "mesh_velocity_U");

    vars->rtp[n] = *ivma -1;
    BFT_MALLOC(vars->name[n], strlen("mesh_velocity_V")+1, char);
    strcpy(vars->name[n++], "mesh_velocity_V");

    vars->rtp[n] = *iwma -1;
    BFT_MALLOC(vars->name[n], strlen("mesh_velocity_W")+1, char);
    strcpy(vars->name[n++], "mesh_velocity_W");

    for (i=k; i < n; i++) {
      BFT_MALLOC(vars->head[i], strlen("ale_method")+1, char);
      strcpy(vars->head[i], "ale_method");
    }
  }

  /* 4) update vars->type for variables */
  k = vars->nvar -vars->nscapp -vars->nscaus;
  for (i=0; i < k; i++) {
    BFT_MALLOC(vars->type[i], strlen("variable")+1, char);
    strcpy(vars->type[i], "variable");
  }

  /* 5) user scalars */

  for (i=0; i < vars->nscaus; i++) {
    vars->rtp[n++] = isca[i] -1;

    BFT_MALLOC(vars->name[k+i], strlen(vars->label[i]) +1, char);
    strcpy(vars->name[k+i], vars->label[i]);

    BFT_MALLOC(vars->type[k+i], strlen("scalar")+1, char);
    strcpy(vars->type[k+i], "scalar");

    BFT_MALLOC(vars->head[k+i], strlen("additional_scalar")+1, char);
    strcpy(vars->head[k+i], "additional_scalar");
  }

  /* 6) model scalars */

  k = vars->nvar -vars->nscapp;
  for (i=0; i < vars->nscapp; i++) {
    j = iscapp[i] -1;
    vars->rtp[n++] = isca[j] -1;

    BFT_MALLOC(vars->name[k+j], strlen(vars->label[j]) +1, char);
    strcpy(vars->name[k+j], vars->label[j]);

    BFT_MALLOC(vars->type[k+j], strlen("scalar")+1, char);
    strcpy(vars->type[k+j], "scalar");

    BFT_MALLOC(vars->head[k+j], strlen(vars->model)+1, char);
    strcpy(vars->head[k+j], vars->model);
  }


  /* 6) check for errors */

  if (n != *nvar)
    bft_error(__FILE__, __LINE__, 0,
              _("The kernel variables number %i and the "
                "calculated one by the GUI %i are not the same.\n"),
                *nvar, n);

#if _XML_DEBUG_
  bft_printf(_("==>CSVNUM\n"));
  bft_printf(_("--variables and scalars name: \n"));
  for (i=0; i < vars->nvar; i++)
    bft_printf(_("---name: %s\n"), vars->name[i]);
  /* for (i=0; i < vars->nscapp+vars->nscaus; i++)
    bft_printf(_("--scalars: %s\n"), vars->label[i]); */
#endif

}

/*----------------------------------------------------------------------------
 * Restart files format.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSIFOA (IFOAVA, IFOAVX)
 * *****************
 *
 * INTEGER          IFOAVA  <--   main restart file format
 * INTEGER          IFOAVX  <--   auxiliary restart file format
 *----------------------------------------------------------------------------*/

void CS_PROCF (csifoa, CSIFOA) (int *const ifoava,
                                int *const ifoavx)
{
  cs_gui_restart_parameters_file_format("main_restart",      ifoava);
  cs_gui_restart_parameters_file_format("auxiliary_restart", ifoavx);

#if _XML_DEBUG_
  bft_printf(_("==>CSIFOA\n"));
  bft_printf(_("--ifoava = %i\n"), *ifoava);
  bft_printf(_("--ifoavx = %i\n"), *ifoavx);
#endif
}

/*----------------------------------------------------------------------------
 * Restart parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSISUI (ISUITE, ILEAUX, ICCVFG)
 * *****************
 *
 * INTEGER          ISUITE  <--   restart
 * INTEGER          ILEAUX  <--   restart with auxiliary
 * INTEGER          ICCFVG  <--   restart with frozen field
 *----------------------------------------------------------------------------*/

void CS_PROCF (csisui, CSISUI) (int *const isuite,
                                int *const ileaux,
                                int *const iccvfg)
{
  cs_gui_restart_parameters_status("restart",                isuite);
  cs_gui_restart_parameters_status("restart_with_auxiliary", ileaux);
  cs_gui_restart_parameters_status("frozen_field",           iccvfg);

#if _XML_DEBUG_
  bft_printf(_("==>CSISUI\n"));
  bft_printf(_("--isuite = %i\n"), *isuite);
  bft_printf(_("--ileaux = %i\n"), *ileaux);
  bft_printf(_("--iccvfg = %i\n"), *iccvfg);
#endif
}

/*----------------------------------------------------------------------------
 * Time passing parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTIME (INPDT0, IPTLTO, NTMABS, DTREF,
 * *****************  DTMIN,  DTMAX,  COUMAX, FOUMAX, VARRDT)
 *
 * INTEGER          INPDT0  <--   zero time step
 * INTEGER          IPTLTO  <--   thermal time step control
 * INTEGER          NTMABS  <--   iterations numbers
 * INTEGER          IDTVAR  <--   time steps'options
 * DOUBLE PRECISION DTREF   <--   time step
 * DOUBLE PRECISION DTMIN   <--   minimal time step
 * DOUBLE PRECISION DTMAX   <--   maximal time step
 * DOUBLE PRECISION COUMAX  <--   maximal courant number
 * DOUBLE PRECISION FOUMAX  <--   maximal fournier number
 * DOUBLE PRECISION VARRDT  <--   max time step variation between 2 iterations
 * DOUBLE PRECISION RELXST  <--   relaxation coefficient id idtvar = -1
 *----------------------------------------------------------------------------*/

void CS_PROCF (cstime, CSTIME) (int    *const inpdt0,
                                int    *const iptlro,
                                int    *const ntmabs,
                                int    *const idtvar,
                                double *const dtref,
                                double *const dtmin,
                                double *const dtmax,
                                double *const coumax,
                                double *const foumax,
                                double *const varrdt,
                                double *const relxst)
{
  double value;

  if (*idtvar == -1){
    cs_gui_steady_parameters("relaxation_coefficient", relxst);

    value =(double) *inpdt0;
    cs_gui_steady_parameters("zero_iteration", &value);
    *inpdt0 = (int) value;

    value =(double) *ntmabs;
    cs_gui_steady_parameters("iterations", &value);
    *ntmabs = (int) value;
  }
  else{
    cs_gui_time_parameters("time_step_ref", dtref);
    cs_gui_time_parameters("time_step_min", dtmin);
    cs_gui_time_parameters("time_step_max", dtmax);
    cs_gui_time_parameters("max_courant_num", coumax);
    cs_gui_time_parameters("max_fourier_num", foumax);
    cs_gui_time_parameters("time_step_var", varrdt);

    value =(double) *ntmabs;
    cs_gui_time_parameters("iterations", &value);
    *ntmabs = (int) value;

    value =(double) *inpdt0;
    cs_gui_time_parameters("zero_time_step", &value);
    *inpdt0 = (int) value;

    value =(double) *iptlro;
    cs_gui_time_parameters("thermal_time_step", &value);
    *iptlro = (int) value;
  }

#if _XML_DEBUG_
  bft_printf(_("==>CSTIME\n"));
  bft_printf(_("--idtvar = %i\n"), *idtvar);
  if (*idtvar == -1){
    bft_printf(_("--inpdt0 = %i\n"), *inpdt0);
    bft_printf(_("--relxst = %i\n"), *relxst);
  }
  else{
    bft_printf(_("--inpdt0 = %i\n"), *inpdt0);
    bft_printf(_("--iptlro = %i\n"), *iptlro);
    bft_printf(_("--ntmabs = %i\n"), *ntmabs);
    bft_printf(_("--dtref = %f\n"),  *dtref);
    bft_printf(_("--dtmin = %f\n"),  *dtmin);
    bft_printf(_("--dtmax = %f\n"),  *dtmax);
    bft_printf(_("--coumax = %f\n"), *coumax);
    bft_printf(_("--foumax = %f\n"), *foumax);
    bft_printf(_("--varrdt = %f\n"), *varrdt);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Check if a users thermal scalar is defined.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSSCA1 (ISCALT, ISCSTH)
 * *****************
 *
 * INTEGER          ISCALT  <--   number of the user thermal scalar if any
 * INTEGER          ISCSTH  <--   type of the user thermal scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca1, CSSCA1) (int *const iscalt,
                                int *const iscsth)
{
#if _XML_DEBUG_
  int i;
  int iphas = 0;
#endif

  assert(vars != NULL);

  cs_gui_thermal_scalar_number(iscalt, iscsth);

#if _XML_DEBUG_
  bft_printf(_("==>CSSCA1\n"));
  bft_printf(_("--iscalt[0]=%i \n"), iscalt[iphas]);
  for (i = 0 ; i < vars->nscaus ; i++)
    bft_printf(_("--iscsth[%i]=%i \n"), i, iscsth[i]);
#endif
}

/*----------------------------------------------------------------------------
 * Traitement des aspects numeriques locaux :
 *     BLENCV, ISCHCV, ISSTPC, IRCFLU, CDTVAR, NITMAX, EPSILO
 *----------------------------------------------------------------------------*/

void CS_PROCF (uinum1, UINUM1) (const    int *const isca,
                                const    int *const iscapp,
                                      double *const blencv,
                                         int *const ischcv,
                                         int *const isstpc,
                                         int *const ircflu,
                                      double *const cdtvar,
                                         int *const nitmax,
                                      double *const epsilo)
{
  int i, j, jj, k;
  double tmp;

  assert(vars != NULL);

  k = vars->nvar - vars->nscaus - vars->nscapp;

  /* 1) variables from velocity_pressure (but not pressure) and turbulence */

  for (i=1; i < k; i++) {
     j = vars->rtp[i];
     cs_gui_variable_value(vars->name[i], "blending_factor", &blencv[j]);
     cs_gui_variable_value(vars->name[i], "solveur_precision", &epsilo[j]);
     tmp = (double) nitmax[j];
     cs_gui_variable_value(vars->name[i], "max_iter_number", &tmp);
     nitmax[j] = (int) tmp;
     cs_gui_variable_attribute(vars->name[i], "order_scheme", &ischcv[j]);
     cs_gui_variable_attribute(vars->name[i], "slope_test", &isstpc[j]);
     cs_gui_variable_attribute(vars->name[i], "flux_reconstruction", &ircflu[j]);
  }

  /* 2) user scalars */

  if (vars->nscaus > 0 ) {
    for (i=0 ; i < vars->nscaus; i++) {
      j = isca[i]-1;
      cs_gui_scalar_value(vars->label[i], "blending_factor", &blencv[j]);
      cs_gui_scalar_value(vars->label[i], "solveur_precision", &epsilo[j]);
      cs_gui_scalar_value(vars->label[i], "time_step_factor", &cdtvar[j]);
      tmp = (double) nitmax[j];
      cs_gui_scalar_value(vars->label[i], "max_iter_number", &tmp);
      nitmax[j] = (int) tmp;
      cs_gui_scalar_attribute(vars->label[i], "order_scheme", &ischcv[j]);
      cs_gui_scalar_attribute(vars->label[i], "slope_test", &isstpc[j]);
      cs_gui_scalar_attribute(vars->label[i], "flux_reconstruction", &ircflu[j]);
    }
  }

  /* 3) model scalars */

  if (vars->nscapp > 0 ) {
    for (i=0 ; i < vars->nscapp ; i++) {
      j = iscapp[i] -1;
      jj = isca[j]-1;
      cs_gui_model_scalar_value(vars->model, vars->label[j], "blending_factor", &blencv[jj]);
      cs_gui_model_scalar_value(vars->model, vars->label[j], "solveur_precision", &epsilo[jj]);
      cs_gui_model_scalar_value(vars->model, vars->label[j], "time_step_factor", &cdtvar[jj]);
      tmp = (double) nitmax[jj];
      cs_gui_model_scalar_value(vars->model, vars->label[j], "max_iter_number", &tmp);
      nitmax[jj] = (int) tmp;
      cs_gui_model_scalar_output_status(vars->model, vars->label[j], "order_scheme", &ischcv[jj]);
      cs_gui_model_scalar_output_status(vars->model, vars->label[j], "slope_test", &isstpc[jj]);
      cs_gui_model_scalar_output_status(vars->model, vars->label[j], "flux_reconstruction", &ircflu[jj]);
    }
  }

#if _XML_DEBUG_
  bft_printf(_("==>UINUM1\n"));
  for (i=0; i < vars->nvar - vars->nscaus - vars->nscapp; i++) {
    bft_printf(_("-->variable[%i] = %s\n"), i, vars->name[i]);
    bft_printf(_("--blencv = %f\n"), blencv[vars->rtp[i]]);
    bft_printf(_("--epsilo = %g\n"), epsilo[vars->rtp[i]]);
    bft_printf(_("--cdtvar = %g\n"), cdtvar[vars->rtp[i]]);
    bft_printf(_("--nitmax = %i\n"), nitmax[vars->rtp[i]]);
    bft_printf(_("--ischcv = %i\n"), ischcv[vars->rtp[i]]);
    bft_printf(_("--isstpc = %i\n"), isstpc[vars->rtp[i]]);
    bft_printf(_("--ircflu = %i\n"), ircflu[vars->rtp[i]]);
  }
  for (i=0 ; i < vars->nscaus + vars->nscapp ; i++) {
    bft_printf(_("-->scalar[%i]: %s\n"), isca[i]-1, vars->label[i]);
    bft_printf(_("--blencv = %f\n"), blencv[isca[i]-1]);
    bft_printf(_("--epsilo = %g\n"), epsilo[isca[i]-1]);
    bft_printf(_("--cdtvar = %g\n"), cdtvar[isca[i]-1]);
    bft_printf(_("--nitmax = %i\n"), nitmax[isca[i]-1]);
    bft_printf(_("--ischcv = %i\n"), ischcv[isca[i]-1]);
    bft_printf(_("--isstpc = %i\n"), isstpc[isca[i]-1]);
    bft_printf(_("--ircflu = %i\n"), ircflu[isca[i]-1]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Global numerical parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSNUM2 (IVISSE, RELAXP, IPUCOU, EXTRAG, IMRGRA)
 * *****************
 * INTEGER          IVISSE  <--   gradient transpose
 * DOUBLE PRECISION RELAXP  <--   pressure relaxation
 * INTEGER          IPUCOU  <--   velocity pressure coupling
 * DOUBLE PRECISION EXTRAG  <--   wall pressure extrapolation
 * INTEGER          IMRGRA  <--   gradient reconstruction
 *----------------------------------------------------------------------------*/

void CS_PROCF (csnum2, CSNUM2)(   int *const ivisse,
                               double *const relaxp,
                                  int *const ipucou,
                               double *const extrag,
                                  int *const imrgra)
{
  cs_gui_numerical_int_parameters("gradient_transpose", ivisse);
  cs_gui_numerical_int_parameters("velocity_pressure_coupling", ipucou);
  cs_gui_numerical_int_parameters("gradient_reconstruction", imrgra);
  cs_gui_numerical_double_parameters("wall_pressure_extrapolation", extrag);
  cs_gui_numerical_double_parameters("pressure_relaxation", relaxp);

#if _XML_DEBUG_
  bft_printf(_("==>CSNUM2\n"));
  bft_printf(_("--ivisse = %i\n"), *ivisse);
  bft_printf(_("--ipucou = %i\n"), *ipucou);
  bft_printf(_("--imrgra = %i\n"), *imrgra);
  bft_printf(_("--extrag = %f\n"), *extrag);
  bft_printf(_("--relaxp = %f\n"), *relaxp);
#endif
}

/*----------------------------------------------------------------------------
 * Traitement de la gravite et des proprietes physiques du fluide
 * Initialisation de la pression de reference et de la temperature si thermique
 *----------------------------------------------------------------------------*/

void CS_PROCF (csphys, CSPHYS)
(
 const    int *const nmodpp,
          int *const irovar,
          int *const ivivar,
       double *const gx,
       double *const gy,
       double *const gz,
       double *const ro0,
       double *const viscl0,
       double *const cp0,
       double *const t0,
       double *const p0
)
{
  int iphas = 0;
  int choice;

  assert(vars != NULL);

  cs_gui_gravity_value("gravity_x", gx);
  cs_gui_gravity_value("gravity_y", gy);
  cs_gui_gravity_value("gravity_z", gz);

  cs_gui_properties_value("density", &ro0[iphas]);
  cs_gui_properties_value("molecular_viscosity", &viscl0[iphas]);
  cs_gui_properties_value("specific_heat", &cp0[iphas]);

  cs_gui_reference_pressure(p0);

  /* rho et viscl variables */
  if (*nmodpp == 0) {
    if (cs_gui_properties_choice("density", &choice))
      irovar[iphas] = choice;

    if (cs_gui_properties_choice("molecular_viscosity", &choice))
      ivivar[iphas] = choice;
  }

  /* T0 si nécessaire */

  if (vars->model != NULL)
    cs_gui_reference_temperature(vars->model, t0);

#if _XML_DEBUG_
  bft_printf(_("==>CSPHYS\n"));
  bft_printf(_("--gx = %f \n"),*gx);
  bft_printf(_("--gy = %f \n"),*gy);
  bft_printf(_("--gz = %f \n"),*gz);
  bft_printf(_("--rho = %g , variable %i\n"), ro0[iphas], irovar[iphas]);
  bft_printf(_("--mu = %g , variable %i \n"), viscl0[iphas], ivivar[iphas]);
  bft_printf(_("--Cp = %g \n"), cp0[0]);
  bft_printf(_("--T0 = %f \n"), *t0);
  bft_printf(_("--P0 = %f \n"), *p0);
#endif
}

/*----------------------------------------------------------------------------
 * User scalar min and max values for clipping.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSSCA2 (ISCAVR, SCAMIN, SCAMAX)
 * *****************
 *
 * INTEGER          ISCAVR   -->  number of the related variance if any
 * DOUBLE PRECISION SCAMIN  <--   user scalar min array
 * DOUBLE PRECISION SCAMAX  <--   user scalar max array
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssca2, CSSCA2) (const    int *const iscavr,
                                      double *const scamin,
                                      double *const scamax)
{
  /* Coal combustion: the min max of the model scalar are not given */

  int i;

  assert(vars != NULL);

  if (vars->nscaus > 0 ) {
    for (i=0 ; i < vars->nscaus; i++) {
      if (iscavr[i] <= 0 ) {
        cs_gui_scalar_value(vars->label[i], "min_value", &scamin[i]);
        cs_gui_scalar_value(vars->label[i], "max_value", &scamax[i]);
      }
    }

#if _XML_DEBUG_
    bft_printf(_("==>CSSCA2\n"));
    for (i=0 ; i < vars->nscaus ; i++) {
      bft_printf(_("--scamin[%i] = %f\n"), i, scamin[i]);
      bft_printf(_("--scamax[%i] = %f\n"), i, scamax[i]);
    }
#endif
  }
}


/*----------------------------------------------------------------------------
 * Lecture de la viscosite dynamique de reference des scalaires utilisateurs
 *----------------------------------------------------------------------------*/


void CS_PROCF (cssca3, CSSCA3) (const    int *const iscalt,
                                const    int *const iscavr,
                                      double *const visls0,
                                      double *const t0,
                                      double *const p0)
{
  int i, iphas = 0;
  double result, coeff, density;

  assert(vars != NULL);

  if (vars->nscaus > 0) {

    if (cs_gui_thermal_scalar()) {
      result = 0;
      cs_gui_properties_value("specific_heat", &result);
      if (!result) {
        bft_error(__FILE__, __LINE__, 0,
                  _("Specific heat value is zero or not found in the xml file.\n"));
      }
      i = iscalt[iphas]-1;
      cs_gui_properties_value("thermal_conductivity", &visls0[i]);
      visls0[i] = visls0[i]/result;
    }

  /* Scalaires utilisateurs */
  /* Dans l'interface, l'utilisateur donne le coefficient de diffusion, alors */
  /* que dans Saturne on parle de la diffusivite, donc il faut multiplier */
  /* ce coefficient par la masse volumique pour etre cohérent */

    for (i=0 ; i < vars->nscaus; i++) {
      if (iscavr[i] <= 0 && i != iscalt[iphas]-1) {

        if (vars->model != NULL) {
          result = 0.028966;
          cs_gui_reference_mass_molar(vars->model, &result);
          if (!result)
            bft_error(__FILE__, __LINE__, 0,
                      _("mass molar value is zero or not found in the xml file.\n"));
          density = *p0 * result / (8.31434 *(*t0));
        }
        else
          cs_gui_properties_value("density", &density);

        if (!density)
          bft_error(__FILE__, __LINE__, 0,
                    _("Density value is zero or not found in the xml file.\n"));

        coeff = visls0[i] / density ;
        cs_gui_scalar_diffusion_value(i+1, &coeff);
        visls0[i] = coeff * density;
      }
    }
#if _XML_DEBUG_
    bft_printf(_("==>CSSCA3\n"));
    for (i=0 ; i < vars->nscaus; i++)
      bft_printf(_("--visls0[%i] = %f\n"), i, visls0[i]);
#endif
  }
}

/*----------------------------------------------------------------------------
 * Turbulence initialization parameters.
 *
 * Fortran Interface:
 *
 * SUBROUTINE CSTINI (UREF, ALMAX)
 * *****************
 *
 * INTEGER          UREF   <--   reference velocity
 * INTEGER          ALMAX  <--   reference length
 *----------------------------------------------------------------------------*/

void CS_PROCF (cstini, CSTINI) (double *const uref,
                                double *const almax)
{
  int iphas = 0;
  char* turb_ini_choice = NULL;

  cs_gui_turbulence_initialization("reference_velocity", &uref[iphas]);

  turb_ini_choice = cs_gui_turbulence_initialization_choice();

  if (cs_gui_strcmp(turb_ini_choice, "reference_velocity_length"))
    cs_gui_turbulence_initialization("reference_length", &almax[iphas]);

  BFT_FREE(turb_ini_choice);

#if _XML_DEBUG_
  bft_printf(_("==>CSTINI\n"));
  bft_printf(_("--almax = %f\n"), almax[iphas]);
  bft_printf(_("--uref  = %f\n"), uref[iphas]);
#endif
}

/*----------------------------------------------------------------------------
 * Tableau des propriétés utilisées dans le calcul
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiprop, UIPROP) (const int *const irom,
                                const int *const iviscl,
                                const int *const ivisct,
                                const int *const ivisls,
                                const int *const icour,
                                const int *const ifour,
                                const int *const ismago,
                                const int *const iale,
                                const int *const icp,
                                const int *const iscalt,
                                const int *const iscavr,
                                const int *const iprtot,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const icmome,
                                const int *const ipptx,
                                const int *const ippty,
                                const int *const ipptz,
                                const int *const ippdt,
                                const int *const ivisma,
                                const int *const iappel)
{
  int iphas = 0;
  int itype = 0;
  int n;
  int i = 0;
  int nbp = 6;
  int ncar = 0;
  char *name;
  char *suf = NULL;

  assert(vars != NULL);

  /* Compute the new size of vars->properties_ipp and vars->properties_name */

  if (ismago[iphas] != -1 ) nbp++;
  if (icp[iphas]>0) nbp++;
  if (vars->nscaus >0) {
    for (i=0; i<vars->nscaus; i++) {
	if (ivisls[i] > 0 && iscavr[i] <= 0 )nbp++;
    }
  }
  if (*iale) {
    cs_gui_get_ale_viscosity_type(&itype);
    if (itype == 1) {
      nbp = nbp + 3;
    } else {
      nbp++;
    }
  }

  n = vars->nprop;

  if (*iappel == 0) {

    vars->nprop = vars->nprop + nbp;

    /* Fisrt step : before the third call of VARPOS in INIUSI */

    BFT_REALLOC(vars->properties_ipp,  vars->nprop, int);
    BFT_REALLOC(vars->properties_name, vars->nprop, char*);

    vars->properties_ipp[n] = ipppro[ ipproc[ irom[iphas]-1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen("density")+1, char);
    strcpy(vars->properties_name[n++], "density");

    vars->properties_ipp[n] = ipppro[ ipproc[ iviscl[iphas]-1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen("molecular_viscosity")+1, char);
    strcpy(vars->properties_name[n++], "molecular_viscosity");

    vars->properties_ipp[n] = ipppro[ ipproc[ ivisct[iphas]-1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen("turb_viscosity")+1, char);
    strcpy(vars->properties_name[n++], "turb_viscosity");

    vars->properties_ipp[n] = ipppro[ ipproc[ icour[iphas]-1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen("courant_number")+1, char);
    strcpy(vars->properties_name[n++], "courant_number");

    vars->properties_ipp[n] = ipppro[ ipproc[ ifour[iphas]-1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen("fourier_number")+1, char);
    strcpy(vars->properties_name[n++], "fourier_number");

    if (ismago[iphas] != -1 ) {
      vars->properties_ipp[n] = ipppro[ ipproc[ ismago[iphas]-1 ]-1 ];
      BFT_MALLOC(vars->properties_name[n], strlen("smagorinsky_constant")+1, char);
      strcpy(vars->properties_name[n++], "smagorinsky_constant");
    }

    if (icp[iphas]>0) {
      vars->properties_ipp[n] = ipppro[ ipproc[ icp[iphas]-1 ]-1 ];
      BFT_MALLOC(vars->properties_name[n], strlen("specific_heat")+1, char);
      strcpy(vars->properties_name[n++], "specific_heat");
    }

    vars->properties_ipp[n] = ipppro[ ipproc[ iprtot[iphas]-1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen("total_pressure")+1, char);
    strcpy(vars->properties_name[n++], "total_pressure");

    if (*iale) {
      vars->properties_ipp[n] = ipppro[ ipproc[ ivisma[0]-1 ]-1 ];
      BFT_MALLOC(vars->properties_name[n], strlen("mesh_viscosity_1")+1, char);
      strcpy(vars->properties_name[n++], "mesh_viscosity_1");

      if (itype == 1) {
        vars->properties_ipp[n] = ipppro[ ipproc[ ivisma[1]-1 ]-1 ];
        BFT_MALLOC(vars->properties_name[n], strlen("mesh_viscosity_2")+1, char);
        strcpy(vars->properties_name[n++], "mesh_viscosity_2");

        vars->properties_ipp[n] = ipppro[ ipproc[ ivisma[2]-1 ]-1 ];
        BFT_MALLOC(vars->properties_name[n], strlen("mesh_viscosity_3")+1, char);
        strcpy(vars->properties_name[n++], "mesh_viscosity_3");
      }

    }

    /* scalar diffusivity */

    if (vars->nscaus > 0) {

      /* search lenght of first character of scalar property's suffixe : '_' */
      for (i=0; i < vars->nscaus; i++) {

        if (iscavr[i] <= 0 && ivisls[i] > 0) {

	    if (iscalt[iphas] == i+1){
		vars->properties_ipp[n] = ipppro[ ipproc[ ivisls[i]-1 ]-1 ];
		BFT_MALLOC(vars->properties_name[n], strlen("thermal_conductivity")+1, char);
		strcpy(vars->properties_name[n++], "thermal_conductivity");
	    } else {
		vars->properties_ipp[n] = ipppro[ ipproc[ ivisls[i]-1 ]-1 ];
		/* search lenght of second character scalar property's suffixe: number i */
		ncar = cs_gui_characters_number(i+1);
		
		BFT_MALLOC(name, strlen("diffusion_coefficient") +2 +ncar, char);
		BFT_MALLOC(suf, 1 + ncar, char);
		sprintf(suf, "%i", i+1);
		strcpy(name, "diffusion_coefficient");
		strcat(name, "_");
		strcat(name, suf);
		BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
		strcpy(vars->properties_name[n++], name);
		BFT_FREE(suf);
	    }
        }
      }
    }

  } else {

    /* Second step : before the fourth call of VARPOS in INIUSI */

    vars->nprop = vars->nprop + 4 + vars->ntimaver;
    BFT_REALLOC(vars->properties_ipp,  vars->nprop, int);
    BFT_REALLOC(vars->properties_name, vars->nprop, char*);

    vars->properties_ipp[n] = *ippdt;
    BFT_MALLOC(vars->properties_name[n], strlen("local_time_step")+1, char);
    strcpy(vars->properties_name[n++], "local_time_step");

    vars->properties_ipp[n] = *ipptx;
    BFT_MALLOC(vars->properties_name[n], strlen("weight_matrix_X")+1, char);
    strcpy(vars->properties_name[n++], "weight_matrix_X");

    vars->properties_ipp[n] = *ippty;
    BFT_MALLOC(vars->properties_name[n], strlen("weight_matrix_Y")+1, char);
    strcpy(vars->properties_name[n++], "weight_matrix_Y");

    vars->properties_ipp[n] = *ipptz;
    BFT_MALLOC(vars->properties_name[n], strlen("weight_matrix_Z")+1, char);
    strcpy(vars->properties_name[n++], "weight_matrix_Z");

    for (i=0; i<vars->ntimaver; i++) {
      vars->properties_ipp[n] = ipppro[ ipproc[ icmome[i]-1 ]-1 ];
      name = cs_gui_get_mean_label(i+1);
      BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
      strcpy(vars->properties_name[n++], name);
    }
  }

  if (n != vars->nprop)
    bft_error(__FILE__, __LINE__, 0,
              _("number of properties is not correct: %i instead of: %i\n"),
              n, vars->nprop);

#if _XML_DEBUG_
  bft_printf(_("==>UIPROP %i\n"),*iappel);
  bft_printf(_("-->nombre de proprietes = %i\n"), vars->nprop);
  for (i=0 ; i<vars->nprop ; i++) {
    bft_printf(_("-->properties_ipp[%i]: %i properties_name[%i]: %s\n"),
                    i, vars->properties_ipp[i], i, vars->properties_name[i]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Traitement des moyennes temporelles
 *----------------------------------------------------------------------------*/

void CS_PROCF (uimoyt, UIMOYT) (const int *const ndgmox,
                                const int *const isca,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const icmome,
                                      int *const ntdmom,
                                      int *const imoold,
                                      int *const idfmom)
{
  int nmean = 0;
  int imom = 0;
  int isuite = 0;
  int i,j,n,nb;
  char *name;

  assert(vars != NULL);

  vars->ntimaver = cs_gui_get_means_number();

  /* for each average */
  for (i=0; i < vars->ntimaver; i++) {

    imom = i + 1;
    cs_gui_get_mean_value(imom, "time_step_start", &ntdmom[i]);

    /* test on isuite */
    cs_gui_restart_parameters_status("restart", &isuite);

    if (isuite != 0) {
      cs_gui_get_mean_value(imom, "restart_from_time_average", &imoold[i]);
      if (imoold[i] == imom) imoold[i] = -2;
    }

    nmean = cs_gui_get_mean_names_number(imom);
    for (n=0; n<nmean; n++) {
      nb = n + 1;
      name = cs_gui_get_mean_prop(imom, nb);
      for (j=0 ; j<(vars->nvar - vars->nscapp - vars->nscaus) ; j++){
        if (cs_gui_strcmp(name,  vars->name[j])) {
          idfmom[(imom-1)*(*ndgmox) + n] = vars->rtp[j] +1;
        }
      }
      if (vars->nscaus > 0 ) {
        for (j=0 ; j<vars->nscaus ; j++) {
          if (cs_gui_strcmp(name,  vars->label[j]) ) {
            idfmom[(imom-1)*(*ndgmox) + n] = isca[j];
          }
        }
      }

      for (j=0 ; j<vars->nprop; j++) {
        if (cs_gui_strcmp(name, vars->properties_name[j]))
          idfmom[(imom-1)*(*ndgmox) + n] = -(vars->properties_ipp[j]);
      }
    }
  }
#if _XML_DEBUG_
  bft_printf(_("==>UIMOYT\n"));
  for (i=0; i<vars->ntimaver; i++) {
    bft_printf(_("-->ntdmom =  %i\n"), ntdmom[i]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Traitement des entrees-sorties
 *----------------------------------------------------------------------------*/

void CS_PROCF (csenso, CSENSO)
(
 const    int *const nvppmx,
          int *const ncapt,
          int *const nthist,
          int *const ntlist,
          int *const ichrvl,
          int *const ichrbo,
          int *const ichrsy,
          int *const ichrmd,
         char *const fmtchr,
          int *const size_fmt,
         char *const optchr,
          int *const size_opt,
          int *const ntchr,
          int *const iecaux,
          int *const ichrvr,
          int *const ilisvr,
          int *const ihisvr,
 const    int *const isca,
 const    int *const iscapp,
 const    int *const ipprtp,
 const    int *const ipppro,
 const    int *const ipproc,
       double *const xyzcap)
{
  int i, j;
  int ipp;

  assert(vars != NULL);

  cs_gui_output_value("fluid_domain", ichrvl);
  cs_gui_output_value("domain_boundary", ichrbo);
  cs_gui_output_value("syrthes_boundary", ichrsy);
  cs_gui_output_value("auxiliary_restart_file_writing", iecaux);
  cs_gui_output_value("listing_printing_frequency", ntlist);
  cs_gui_output_value("postprocessing_frequency", ntchr);
  cs_gui_output_value("probe_recording_frequency", nthist);
  cs_gui_output_value("postprocessing_mesh_options", ichrmd);
  cs_gui_output_choice("postprocessing_format", fmtchr, size_fmt);
  cs_gui_output_choice("postprocessing_options", optchr, size_opt);

  *ncapt = cs_gui_probes_number();
  for (i=0; i < *ncapt; i++) {
    xyzcap[0 + i*3] = cs_gui_probe_coordinate(i+1, "probe_x");
    xyzcap[1 + i*3] = cs_gui_probe_coordinate(i+1, "probe_y");
    xyzcap[2 + i*3] = cs_gui_probe_coordinate(i+1, "probe_z");
  }

  /* Sorties des variables vitesses et turbulence */
  for (i=0; i<vars->nvar - vars->nscaus - vars->nscapp; i++) {
     ipp = ipprtp[vars->rtp[i]];
     cs_gui_thermophysical_post(vars->name[i],
                                ipp,
                                ihisvr, ilisvr, ichrvr,
                                nvppmx);
  }

  /* User scalar */
  if (vars->nscaus > 0 ) {
    for (i=0 ; i < vars->nscaus; i++) {
      cs_gui_scalar_post(i, ihisvr, ilisvr, ichrvr,
                         ipprtp, isca, nvppmx);
    }
  }

  /* Scalaire model */
  if (vars->nscapp > 0) {
    for (i=0 ; i < vars->nscapp; i++) {
      j = iscapp[i]-1 ;
      cs_gui_model_scalar_post(vars->model, j,
                               ihisvr, ilisvr, ichrvr,
                               ipprtp, isca, nvppmx);
    }
  }

  /* Proprietes physiques */

  if (vars->nsalpp > 0) {
    for (i=0 ; i < vars->nsalpp; i++) {
      cs_gui_model_property_post(vars->model, i,
                                 ihisvr, ilisvr, ichrvr,
                                 ipppro, ipproc, nvppmx);
    }
  }

  for (i=vars->nsalpp ; i < vars->nprop ; i++) {
    if (vars->ntimaver != 0 && i>=vars->nprop-vars->ntimaver){
      cs_gui_time_average_post(vars->properties_name[i],
                               vars->properties_ipp[i],
                               ichrvr,
                               ilisvr,
                               ihisvr,
                               nvppmx);
    }
    else
      cs_gui_properties_post(vars->properties_name[i],
                             vars->properties_ipp[i],
                             ichrvr,
                             ilisvr,
                             ihisvr,
                             nvppmx);
  }

#if _XML_DEBUG_
  bft_printf(_("==>CSENSO\n"));
  bft_printf(_("--iecaux = %i\n"), *iecaux);
  bft_printf(_("--ichrvl = %i\n"), *ichrvl);
  bft_printf(_("--ichrbo = %i\n"), *ichrbo);
  bft_printf(_("--ichrsy = %i\n"), *ichrsy);
  bft_printf(_("--fmtchr = %s\n"), "à vérifier en fortran");
  bft_printf(_("--optchr = %s\n"), "à vérifier en fortran");
  bft_printf(_("--ntlist = %i\n"), *ntlist);
  bft_printf(_("--ntchr  = %i\n"), *ntchr);
  bft_printf(_("--nthist = %i\n"), *nthist);
  bft_printf(_("--ncapt  = %i\n"), *ncapt);
  for (i = 0; i < *ncapt; i++) {
    bft_printf(_("--xyzcap[%i][0] = %f\n"), i, xyzcap[0 +i*3]);
    bft_printf(_("--xyzcap[%i][1] = %f\n"), i, xyzcap[1 +i*3]);
    bft_printf(_("--xyzcap[%i][2] = %f\n"), i, xyzcap[2 +i*3]);
  }
  for (i=0; i < vars->nvar - vars->nscaus - vars->nscapp; i++){
    ipp = ipprtp[vars->rtp[i]];
    bft_printf(_("-->variable ipprtp[%i] = %s\n"), ipp, vars->name[i]);
    bft_printf(_("--ichrvr[%i] = %i \n"), ipp, ichrvr[ipp-1]);
    bft_printf(_("--ilisvr[%i] = %i \n"), ipp, ilisvr[ipp-1]);
    bft_printf(_("--ihisvr[0][%i]= %i \n"), ipp, ihisvr[0 + (ipp-1)]);
    if (ihisvr[0 + (ipp-1)]>0)
      for (j=0; j<ihisvr[0 + (ipp-1)]; j++)
        bft_printf(_("--ihisvr[%i][%i]= %i \n"), j+1, ipp,
                      ihisvr[(j+1)*(*nvppmx) + (ipp-1)]);
  }
  for (i=0 ; i < vars->nscaus + vars->nscapp ; i++) {
    ipp = ipprtp[isca[i] -1];
    bft_printf(_("-->scalar ipprtp[%i]: %s\n"), ipp, vars->label[i]);
    bft_printf(_("--ichrvr[%i] = %i \n"), ipp, ichrvr[ipp-1]);
    bft_printf(_("--ilisvr[%i] = %i \n"), ipp, ilisvr[ipp-1]);
    bft_printf(_("--ihisvr[0][%i]= %i \n"), ipp, ihisvr[0 + (ipp-1)]);
    if (ihisvr[0 + (ipp-1)]>0)
      for (j=0; j<ihisvr[0 + (ipp-1)]; j++)
        bft_printf(_("--ihisvr[%i][%i]= %i \n"), j+1, ipp,
                      ihisvr[(j+1)*(*nvppmx) + (ipp-1)]);
  }
  for (i=0 ; i<vars->nprop ; i++) {
    ipp = vars->properties_ipp[i];
    bft_printf(_("-->properties_name[%i]: %s\n"), i, vars->properties_name[i]);
    bft_printf(_("--ichrvr[%i] = %i \n"), ipp, ichrvr[ipp-1]);
    bft_printf(_("--ilisvr[%i] = %i \n"), ipp, ilisvr[ipp-1]);
    bft_printf(_("--ihisvr[0][%i]= %i \n"), ipp, ihisvr[0 + (ipp-1)]);
    if (ihisvr[0 + (ipp-1)]>0)
      for (j=0; j<ihisvr[0 + (ipp-1)]; j++)
        bft_printf(_("--ihisvr[%i][%i]= %i \n"), j+1, ipp,
                      ihisvr[(j+1)*(*nvppmx) + (ipp-1)]);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Users arrays
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIUSAR (ICOFTU)
 * *****************
 *
 * INTEGER          ICOFTU   -->  Dimension coef for user arrays
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiusar, UIUSAR) (int *const icoftu)
{
  icoftu[0] = _user_array("integer_user_array", "ncelet");
  icoftu[1] = _user_array("integer_user_array", "nfac");
  icoftu[2] = _user_array("integer_user_array", "nfabor");
  icoftu[3] = _user_array("integer_user_array", "dimless");

  icoftu[4] =_user_array("real_user_array", "ncelet");
  icoftu[5] = _user_array("real_user_array", "nfac");
  icoftu[6] = _user_array("real_user_array", "nfabor");
  icoftu[7] = _user_array("real_user_array", "dimless");


#if _XML_DEBUG_
  bft_printf(_("==>UIUSAR\n"));
  bft_printf(_("--icoftu = %i %i %i %i\n"),
                icoftu[0],icoftu[1],icoftu[2],icoftu[3]);
  bft_printf(_("           %i %i %i %i\n"),
                icoftu[4],icoftu[5],icoftu[6],icoftu[7]);
#endif
}


/*----------------------------------------------------------------------------
 * Copy variable name from Fortran to C
 *----------------------------------------------------------------------------*/

void CS_PROCF(fcnmva, FCNMVA)
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

    BFT_REALLOC(_cs_gui_var_name, _cs_gui_max_vars, char *);
    for (i = _cs_gui_last_var; i < _cs_gui_max_vars; i++)
      _cs_gui_var_name[i] = NULL;
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
  assert(_cs_gui_var_name[*var_id - 1] == NULL);

  if (l > 0) {

    /* Allocate and copy */
    BFT_MALLOC(cstr, l + 1, char);

  for (i = 0 ; i < l ; i++, i1++)
    cstr[i] = fstr[i1];

  cstr[l] = '\0';

    _cs_gui_var_name[*var_id - 1] = cstr;

  }

  /* Update variable counter */
  _cs_gui_last_var = *var_id;

}


/*----------------------------------------------------------------------------
 * Copy variable name from C to Fortran
 *----------------------------------------------------------------------------*/


void CS_PROCF(cfnmva, CFNMVA)
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

  cstr = _cs_gui_var_name[*var_id - 1];

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
 * Variables and user scalars initialization
 *
 * Fortran Interface:
 *
 * SUBROUTINE UIINIV (NCELET, ISCA, RTP)
 * *****************
 *
 * INTEGER          ISCAVR   -->  number of cells
 * INTEGER          ISCA     -->  indirection array for scalar number
 * DOUBLE PRECISION RTP     <--   variables and scalars array
 *----------------------------------------------------------------------------*/

void CS_PROCF(uiiniv, UIINIV)(const int    *const ncelet,
                              const int    *const isca,
                                    double *const rtp)
{
  /* Coal combustion: the initialization of the model scalar are not given */

  int i, j, icel, iel, c_id;
  int zones = 0;
  int cells = 0;
  int *cells_list = NULL;
  double initial_value = 0;
  char *choice = NULL;
  char *name = NULL;
  char *description = NULL;

  assert(vars != NULL);

  /* number of volumic zone */

  zones = cs_gui_volumic_zones_number();

#if _XML_DEBUG_
  bft_printf(_("==>UIINIV\n"));
  bft_printf(_("--initialization zones number: %i\n"), zones);
#endif

  for (i=1; i < zones+1; i++) {

    /* name and description (color or group) of the ith initialization zone */
    name = cs_gui_volumic_zone_name(i);
    description = cs_gui_volumic_zone_localization(name);

    /* build list of cells */
    BFT_MALLOC(cells_list, *ncelet, int);

    c_id = fvm_selector_get_list(cs_glob_mesh->select_cells,
                                 description,
                                 &cells,
                                 cells_list);

    if (fvm_selector_n_missing(cs_glob_mesh->select_cells, c_id) > 0) {
      const char *missing
        = fvm_selector_get_missing(cs_glob_mesh->select_cells, c_id, 0);
      bft_error(__FILE__, __LINE__, 0,
                _("Le groupe ou attribut \"%s\" figurant dans le\n"
                  "critère de sélection:\n"
                  "\"%s\"\n ne correspond à aucune cellule."),
                missing, description);
    }


    /* Velocity variables initialization */
    for (j=1; j < 4; j++) {

      cs_gui_variable_initial_value(vars->name[j], name, &initial_value);

      for (icel = 0; icel < cells; icel++) {
        iel = cells_list[icel]-1;
        rtp[vars->rtp[j]*(*ncelet) + iel] = initial_value;
      }
    }

    /* Turbulence variables initialization */
    choice = cs_gui_turbulence_initialization_choice();

    if (cs_gui_strcmp(choice, "values")) {
      for (j=4; j < vars->nvar - vars->nscaus - vars->nscapp; j++) {

        cs_gui_variable_initial_value(vars->name[j], name, &initial_value);

        for (icel = 0; icel < cells; icel++) {
          iel = cells_list[icel]-1;
          rtp[vars->rtp[j]*(*ncelet) + iel] = initial_value;
        }
      }
    }
    BFT_FREE(choice);

    /* User Scalars initialization */
    for (j=0; j < vars->nscaus; j++) {

      cs_gui_scalar_initial_value("additional_scalars",
                                  vars->label[j],
                                  name,
                                  &initial_value);

      for (icel = 0; icel < cells; icel++) {
        iel = cells_list[icel]-1;
        rtp[(isca[j]-1)*(*ncelet) + iel] = initial_value;
      }
    }
    BFT_FREE(cells_list);

#if _XML_DEBUG_
    bft_printf(_("--zone name and description: %s, %s\n"), name, description);
    bft_printf(_("--zone's element number: %i\n"), cells);

    for (j=1; j < vars->nvar - vars->nscaus - vars->nscapp; j++){
      cs_gui_variable_initial_value(vars->name[j], name, &initial_value);
      bft_printf(_("--initial value for %s: %f\n"),
        vars->name[j], initial_value);
    }

    for (j=0; j < vars->nscaus; j++) {
      cs_gui_scalar_initial_value("additional_scalars",
                                  vars->label[j],
                                  name,
                                  &initial_value);
      bft_printf(_("--initial value for %s: %f\n"), vars->label[j], initial_value);
    }
#endif
    BFT_FREE(name);
    BFT_FREE(description);
  } /* zones+1 */
}

/*----------------------------------------------------------------------------
 * Boundary conditions treatment
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICLIM
 * *****************
 *
 * INTEGER          NOZPPM  --> max number of boundary conditions zone
 * INTEGER          NFABOR  --> number of boundary faces
 * INTEGER          IINDEF  --> type of boundary: not defined
 * INTEGER          IENTRE  --> type of boundary: inlet
 * INTEGER          IPAROI  --> type of boundary: smooth wall
 * INTEGER          IPARUG  --> type of boundary: rough wall
 * INTEGER          ISYMET  --> type of boundary: symetry
 * INTEGER          ISOLIB  --> type of boundary: outlet
 * INTEGER          IQIMP   --> 1 if flow rate is applied
 * INTEGER          ICALKE  --> 1 for automatic turbulent boundary conditions
 * INTEGER          ITYPFB  --> type of boundary for each face
 * INTEGER          IZFPPP  --> zone number
 * INTEGER          ICODCL  --> boundary conditions array type
 * DOUBLE PRECISION QIMP    --> flow rate value if applied
 * DOUBLE PRECISION DH      --> hydraulic diameter
 * DOUBLE PRECISION XINTUR  --> turbulent intensity
 * DOUBLE PRECISION RCODCL  --> boundary conditions array value
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclim, UICLIM)(const    int *const nozppm,
                               const    int *const nfabor,
                               const    int *const iindef,
                               const    int *const ientre,
                               const    int *const iparoi,
                               const    int *const iparug,
                               const    int *const isymet,
                               const    int *const isolib,
                                        int *const iqimp,
                                        int *const icalke,
                                        int *const itypfb,
                                        int *const izfppp,
                                        int *const icodcl,
                                     double *const qimp,
                                     double *const dh,
                                     double *const xintur,
                                     double *const rcodcl)
{
  int iphas = 0;
  int faces = 0;
  int zones = 0;
  int ifac, izone, ith_zone, zone_nbr;
  int ifbr, i, iwall, c_id;
  int ivar, isca;
  double qimp2 = 0.;
  double timp = 0.;
  char *choice = NULL;
  char *nature = NULL;
  char *label = NULL;
  char *description = NULL;
  int *faces_list = NULL;

  assert(vars != NULL);

  zones   = cs_gui_boundary_zones_number();

  /* First iteration only: memory allocation */

  if (boundaries == NULL) {

    BFT_MALLOC(boundaries,            1,          cs_boundary_t);
    BFT_MALLOC(boundaries->label,     zones,      char*        );
    BFT_MALLOC(boundaries->nature,    zones,      char*        );
    BFT_MALLOC(boundaries->type_code, vars->nvar, int*         );
    BFT_MALLOC(boundaries->values,    vars->nvar, cs_val_t*    );
    BFT_MALLOC(boundaries->iqimp,     zones,      int          );
    BFT_MALLOC(boundaries->qimp,      zones,      double       );
    BFT_MALLOC(boundaries->icalke,    zones,      int          );
    BFT_MALLOC(boundaries->dh,        zones,      double       );
    BFT_MALLOC(boundaries->xintur,    zones,      double       );
    BFT_MALLOC(boundaries->rough,     zones,      double       );

    boundaries->ientat = NULL;
    boundaries->qimpat = NULL;
    boundaries->timpat = NULL;
    boundaries->ientcp = NULL;
    boundaries->qimpcp = NULL;
    boundaries->timpcp = NULL;
    boundaries->distch = NULL;

    for (ivar = 0; ivar < vars->nvar; ivar++) {
      i = vars->rtp[ivar];
      BFT_MALLOC(boundaries->type_code[i], zones, int);
      BFT_MALLOC(boundaries->values[i], zones, cs_val_t);
    }

    for (izone = 0; izone < zones; izone++) {
      boundaries->iqimp[izone]  = 0;
      boundaries->qimp[izone]   = 0;
      boundaries->icalke[izone] = 0;
      boundaries->dh[izone]     = 0;
      boundaries->xintur[izone] = 0;
      boundaries->rough[izone]  = -999;
    }

    /* Initialization of boundary->type_code and boundary->values */

    for (ivar = 0; ivar < vars->nvar; ivar++) {
      i = vars->rtp[ivar];
      for (izone = 0; izone < zones; izone++) {
        boundaries->type_code[i][izone] = -1;
        boundaries->values[i][izone].val1 = 1.e30;
        boundaries->values[i][izone].val2 = 1.e30;
        boundaries->values[i][izone].val3 = 0.;
      }
    }

   /* First iteration only : filling of the "boundaries" structure */

    for (izone = 0; izone < zones; izone++) {

     /* nature, label of the ith initialization zone */

      ith_zone = izone + 1;
      nature = cs_gui_boundary_zone_nature(ith_zone);
      label = cs_gui_boundary_zone_label(ith_zone);

      BFT_MALLOC(boundaries->label[izone], strlen(label)+1, char);
      strcpy(boundaries->label[izone], label);

      BFT_MALLOC(boundaries->nature[izone], strlen(nature)+1, char);
      strcpy(boundaries->nature[izone], nature);

      if (cs_gui_strcmp(nature, "inlet")) {

        /* Inlet: VELOCITY */
        choice = cs_gui_boundary_choice("inlet", label, "velocity_pressure");

        if (cs_gui_strcmp(choice, "dirichlet")) {

          for (ivar = 1; ivar < 4; ivar++) {
            boundaries->iqimp[izone] = 0;
            cs_gui_boundary_dirichlet("inlet", label, izone, ivar);
          }

        } else if (cs_gui_strcmp(choice, "flow1")) {

          boundaries->iqimp[izone] = 1;
            cs_gui_boundary_flow(label,&qimp2,&timp);
            boundaries->qimp[izone] = qimp2;
          /* TODO : remplir la direction normale a la face */
          /* boundaries->values[1][izone].val1 = directionU ; */
          /* boundaries->values[2][izone].val1 = directionV ; */
          /* boundaries->values[3][izone].val1 = directionW ; */
        }

        BFT_FREE(choice);

        /* Inlet: TURBULENCE */
        choice = cs_gui_boundary_choice("inlet", label, "turbulence");
        cs_gui_boundary_turbulence(choice, izone);
        BFT_FREE(choice);

        /* Inlet: USER SCALARS */
        for (isca = 0; isca < vars->nscaus; isca++) {
          cs_gui_boundary_value_scalar("inlet", izone, isca);
        }

      } else if (cs_gui_strcmp(nature, "wall")) {

        /* Wall: VELOCITY */
        choice = cs_gui_boundary_choice("wall", label, "velocity_pressure");

        if (cs_gui_strcmp(choice, "on")) {
          for (ivar = 1; ivar < 4; ivar++)
            cs_gui_boundary_dirichlet("wall", label, izone, ivar);
        }
        BFT_FREE(choice);

        /* Wall: ROUGH */
        cs_gui_boundary_rough(label, izone);

        /* Wall: USER SCALARS */
        for (isca = 0; isca < vars->nscaus; isca++) {
          cs_gui_boundary_value_scalar("wall", izone, isca);
        }

      } else if (cs_gui_strcmp(nature, "outlet")) {

        /* Outlet: USER SCALARS */
        for (isca = 0; isca < vars->nscaus; isca++) {
          cs_gui_boundary_value_scalar("outlet", izone, isca);
        }

      }  /* if (cs_gui_strcmp(nature, "outlet")) */

      BFT_FREE(nature);
      BFT_FREE(label);

    }  /* for izones */

  }  /* if (boundaries == NULL)*/

 /* A chaque iteration, boucle sur les faces de bord :
    on remplit itypfb, rcodcl et icodcl a partir des tableaux
    de la structures conditions.limites definie
    dans la premiere partie de la fonction.
    Remember: rdoccl[k][j][i] = rcodcl[ k * dim1 *dim2 + j *dim1 + i] */

  for (izone=0 ; izone < zones ; izone++) {

    ith_zone = izone + 1;
    zone_nbr = cs_gui_boundary_zone_number(ith_zone);
    if (zone_nbr > *nozppm)
      bft_error(__FILE__, __LINE__, 0,
                _("zone's label number %i is greater than %i,"
                  " the maximum allowed \n"),
           zone_nbr, *nozppm);

    description = cs_gui_boundary_zone_localization(boundaries->nature[izone],
                                                    boundaries->label[izone]);
    /* list of faces building */
    BFT_MALLOC(faces_list, *nfabor, int);

    c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                                 description,
                                 &faces,
                                 faces_list);

    if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
      const char *missing
        = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
      bft_error(__FILE__, __LINE__, 0,
                _("Le groupe ou attribut \"%s\" figurant dans le\n"
                  "critère de sélection:\n"
                  "\"%s\"\n ne correspond à aucune face de bord."),
                missing, description);
    }

    BFT_FREE(description);

    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {

      /* Update the depending zone arrays (iqimp, dh, xintur, icalke, qimp)
         because they are initialized at each time step in PRECLI routine */

      iqimp[zone_nbr-1]  = boundaries->iqimp[izone];
      dh[zone_nbr-1]     = boundaries->dh[izone];
      xintur[zone_nbr-1] = boundaries->xintur[izone];
      icalke[zone_nbr-1] = boundaries->icalke[izone];
      qimp[zone_nbr-1]   = boundaries->qimp[izone];

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = *ientre;
        for (i = 0; i < vars->nvar; i++) {
          ivar = vars->rtp[i];
          rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
          = boundaries->values[ivar][izone].val1;
        }
      }

    } else if (cs_gui_strcmp(boundaries->nature[izone], "wall")) {

      if (boundaries->rough[izone] >= 0.) {
        iwall = *iparug;
        /* roughness value is only stored in Velocity_U */
        ivar = 1;
        for (ifac = 0; ifac < faces; ifac++) {
          ifbr = faces_list[ifac]-1;
          icodcl[ivar *(*nfabor) + ifbr] = 6;
          rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
          = boundaries->rough[izone];
        }
      } else {
        iwall = *iparoi;
      }

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = iwall;
      }

      for (i = 0; i < vars->nvar; i++) {
        ivar = vars->rtp[i];

        switch (boundaries->type_code[ivar][izone]) {

          case NEUMANN :

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]-1;
              icodcl[ivar *(*nfabor) + ifbr] = 3;
              rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val3;
            }
          break;

          case DIRICHLET :

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]-1;
              icodcl[ivar *(*nfabor) + ifbr] = 5;
              /* si wall_function icodcl[ivar *(*nfabor) + ifbr] = 1; */
              rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1;
            }
          break;

          case WALL_FUNCTION :

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]-1;
              icodcl[ivar *(*nfabor) + ifbr] = 5;
              rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1;
            }
          break;

          case COEF_ECHANGE :

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]-1;
              icodcl[ivar *(*nfabor) + ifbr] = 5;
              rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1;
              rcodcl[1 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val2;
            }
          break;
        }
      }

    } else if (cs_gui_strcmp(boundaries->nature[izone], "outlet")) {

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = *isolib;
      }

      for (i = 0; i < vars->nvar; i++) {
        ivar = vars->rtp[i];

        switch (boundaries->type_code[ivar][izone]) {

          case DIRICHLET :

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]-1;
              icodcl[ivar *(*nfabor) + ifbr] = 1;
              rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1;
            }
          break;

        }
      }

    } else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = *isymet;
      }

    } else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = *iindef;
      }

    } else {
        bft_error(__FILE__, __LINE__, 0,
                  _("boundary nature %s is unknown \n"), boundaries->nature[izone]);
    }
    BFT_FREE(faces_list);
  } /*  for izone */

#if _XML_DEBUG_
  bft_printf(_("==>UICLIM\n"));
  bft_printf(_("--boundary zones number: %i\n"), zones);

  for (izone=0 ; izone < zones ; izone++) {

    BFT_MALLOC(faces_list, *nfabor, int);
    description = cs_gui_boundary_zone_localization(boundaries->nature[izone],
                                                    boundaries->label[izone]);

    c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                                 description,
                                 &faces,
                                 faces_list);

    if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
      const char *missing
        = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
      bft_error(__FILE__, __LINE__, 0,
                _("Le groupe ou attribut \"%s\" figurant dans le\n"
                  "critère de sélection:\n"
                  "\"%s\"\n ne correspond à aucune face de bord."),
                missing, description);
    }

    zone_nbr = cs_gui_boundary_zone_number(izone+1);

    bft_printf(_("---zone %i label: %s\n"), zone_nbr, boundaries->label[izone]);
    bft_printf(_("---zone %i nature: %s\n"), zone_nbr, boundaries->nature[izone]);
    bft_printf(_("---zone %i number of faces: %i\n"), zone_nbr, faces);
    bft_printf(_("----localization: %s\n"), description);
    BFT_FREE(description);

    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {
      bft_printf(_("-----iqimp=%i, qimp=%12.5e \n"),
                   iqimp[zone_nbr-1], qimp[zone_nbr-1]);
      bft_printf(_("-----icalke=%i, dh=%12.5e, xintur=%12.5e \n"),
                   icalke[zone_nbr-1], dh[zone_nbr-1], xintur[zone_nbr-1]);
    }

    if (faces>0) {
        ifbr = faces_list[0]-1;
        for (i=0; i<vars->nvar - vars->nscaus - vars->nscapp ; i++) {
          ivar = vars->rtp[i];
          bft_printf(_("-----%s: itypfb=%i, icodcl=%i, "
                           "rcodcl(1)=%12.5e, rcodcl(2)=%12.5e, rcodcl(3)=%12.5e\n"),
                         vars->char2[ivar],
                         itypfb[iphas *(*nfabor) + ifbr],
                         icodcl[ivar *(*nfabor) + ifbr ],
                         rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr],
                         rcodcl[1 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr],
                         rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]);
        }

        for (ivar=0; ivar<vars->nscaus + vars->nscapp ; ivar++) {
          bft_printf(_("-----%s: itypfb=%i, icodcl=%i, "
                           "rcodcl(1)=%12.5e, rcodcl(2)=%12.5e, rcodcl(3)=%12.5e\n"),
                         vars->label[ivar],
                         itypfb[iphas *(*nfabor) + ifbr],
                         icodcl[ivar *(*nfabor) + ifbr ],
                         rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr],
                         rcodcl[1 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr],
                         rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]);
        }
    }
    BFT_FREE(faces_list);
  }
#endif
}

/*----------------------------------------------------------------------------
 * Density under relaxation
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICPI1 (SRROM)
 * *****************
 * DOUBLE PRECISION SRROM   <--   density relaxation
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpi1, UICPI1) (double *const srrom)
{
  cs_gui_numerical_double_parameters("density_relaxation", srrom);

#if _XML_DEBUG_
  bft_printf(_("==>UICPI1\n"));
  bft_printf(_("--srrom  = %f\n"), *srrom);
#endif
}

/*-----------------------------------------------------------------------------
 *  Indirection entre la numérotation noyau et la numérotation XML
 *  des propriétés physiques de la physique particulière active
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicppr, UICPPR) (const int *const nclass,
                                const int *const nsalpp,
                                const int *const nsalto,
                                const int *const ippmod,
                                const int *const icp3pl,
                                const int *const ipppro,
                                const int *const ipproc,
                                const int *const itemp1,
                                const int *const irom1,
                                const int *const ym1,
                                const int *const imel,
                                const int *const itemp2,
                                const int *const ix2,
                                const int *const irom2,
                                const int *const idiam2,
                                const int *const igmdch,
                                const int *const igmdv1,
                                const int *const igmdv2,
                                const int *const igmhet,
                                const int *const igmsec,
                                const int *const ilumi)
{
  int i = 0;
  int n;
  char *name = NULL;
  char *snumpp = NULL;

  assert(vars != NULL);

  n = vars->nprop;
  vars->nprop  = *nsalpp;
  vars->nsalpp = *nsalpp;

  BFT_MALLOC(vars->properties_ipp,  vars->nsalpp, int);
  BFT_MALLOC(vars->properties_name, vars->nsalpp, char*);

 /* ITEMP1 */
  vars->properties_ipp[n] = ipppro[ ipproc[ *itemp1 -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("Temp_GAZ")+1, char);
  strcpy(vars->properties_name[n++], "Temp_GAZ");

 /* IROM1 */
  vars->properties_ipp[n] = ipppro[ ipproc[ *irom1 -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("ROM_GAZ")+1, char);
  strcpy(vars->properties_name[n++], "ROM_GAZ");

 /*  YM_CHX1M */
  vars->properties_ipp[n] = ipppro[ ipproc[ ym1[0] -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CHx1m")+1, char);
  strcpy(vars->properties_name[n++], "YM_CHx1m");

 /*  YM_CHX2M */
  vars->properties_ipp[n] = ipppro[ ipproc[ ym1[1] -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CHx2m")+1, char);
  strcpy(vars->properties_name[n++], "YM_CHx2m");

 /*  YM_CO */
  vars->properties_ipp[n] = ipppro[ ipproc[ ym1[2] -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CO")+1, char);
  strcpy(vars->properties_name[n++], "YM_CO");

 /*  YM_O2 */
  vars->properties_ipp[n] = ipppro[ ipproc[ ym1[3] -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_O2")+1, char);
  strcpy(vars->properties_name[n++], "YM_O2");

 /*  YM_CO2 */
  vars->properties_ipp[n] = ipppro[ ipproc[ ym1[4] -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_CO2")+1, char);
  strcpy(vars->properties_name[n++], "YM_CO2");

 /*  YM_H2O */
  vars->properties_ipp[n] = ipppro[ ipproc[ ym1[5] -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_H2O")+1, char);
  strcpy(vars->properties_name[n++], "YM_H2O");

 /*  YM_N2 */
  vars->properties_ipp[n] = ipppro[ ipproc[ ym1[6] -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("YM_N2")+1, char);
  strcpy(vars->properties_name[n++], "YM_N2");

 /* IMEL */
  vars->properties_ipp[n] = ipppro[ ipproc[ *imel -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("XM")+1, char);
  strcpy(vars->properties_name[n++], "XM");

 /* ITEMP2 boucle sur les classes */
  BFT_MALLOC(name, strlen("Temp_CP")+1 + 2, char);
  BFT_MALLOC(snumpp, 1 + 2, char);
  strcpy(name, "Temp_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ itemp2[i] -1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Temp_CP");
  }

 /* IX2 boucle sur les classes */
  BFT_REALLOC(name, strlen("Frm_CP")+1 + 2, char);
  strcpy(name, "Frm_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ ix2[i] -1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Frm_CP");
  }

 /* IROM2 boucle sur les classes */
  BFT_REALLOC(name, strlen("Rho_CP")+1 + 2, char);
  strcpy(name, "Rho_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ irom2[i] -1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Rho_CP");
  }

 /* IDIAM2 boucle sur les classes */
  BFT_REALLOC(name, strlen("Dia_CK")+1 + 2, char);
  strcpy(name, "Dia_CK");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ idiam2[i] -1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Dia_CK");
  }

 /* IGMDCH boucle sur les classes */
  BFT_REALLOC(name, strlen("Ga_DCH")+1 + 2, char);
  strcpy(name, "Ga_DCH");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ igmdch[i] -1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_DCH");
  }

 /* IGMDV1 boucle sur les classes */
  BFT_REALLOC(name, strlen("Ga_DV1")+1 + 2, char);
  strcpy(name, "Ga_DV1");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ igmdv1[i] -1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_DV1");
  }

 /* IGMDV2 boucle sur les classes */
  BFT_REALLOC(name, strlen("Ga_DV2")+1 + 2, char);
  strcpy(name, "Ga_DV2");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ igmdv2[i] -1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_DV2");
  }

 /* IGMHET boucle sur les classes */
  BFT_REALLOC(name, strlen("Ga_HET")+1 + 2, char);
  strcpy(name, "Ga_HET");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumpp, "%2.2i", i+1);
    strcat(name, snumpp);

    vars->properties_ipp[n] = ipppro[ ipproc[ igmhet[i] -1 ]-1 ];
    BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
    strcpy(vars->properties_name[n++], name);

    strcpy(name, "Ga_HET");
  }

  if (ippmod[*icp3pl -1] == 1) {
   /* IGMSEC boucle sur les classes */
    BFT_REALLOC(name, strlen("Ga_SEC")+1 + 2, char);
    strcpy(name, "Ga_SEC");
    for (i = 0; i < *nclass; i++) {
      sprintf(snumpp, "%2.2i", i+1);
      strcat(name, snumpp);

      vars->properties_ipp[n] = ipppro[ ipproc[ igmsec[i] -1 ]-1 ];
      BFT_MALLOC(vars->properties_name[n], strlen(name)+1, char);
      strcpy(vars->properties_name[n++], name);

      strcpy(name, "Ga_SEC");
    }
  }

 /* ILUMI */
  vars->properties_ipp[n] = ipppro[ ipproc[ *ilumi -1 ]-1 ];
  BFT_MALLOC(vars->properties_name[n], strlen("ntLuminance_4PI")+1, char);
  strcpy(vars->properties_name[n++], "ntLuminance_4PI");

  BFT_FREE(name);
  BFT_FREE(snumpp);

  if (n != vars->nsalpp)
    bft_error(__FILE__, __LINE__, 0,
              _("number of properties is not correct: %i instead of: %i\n"),
                n, vars->nsalpp);

#if _XML_DEBUG_
  bft_printf(_("==>UICPPR\n"));
  bft_printf(_("-->nombre de proprietes = %i\n"), vars->nprop);
  for (i=0 ; i<vars->nprop ; i++)
    bft_printf(_("-->properties_ipp[%i]: %i properties_name[%i]: %s\n"),
                    i, vars->properties_ipp[i], i, vars->properties_name[i]);
#endif
}

/*------------------------------------------------------------------------------
 *  Indirection entre la numérotation noyau et
 *  la numérotation XML des scalaires model
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpsc, UICPSC) (const int *const ncharb,
                                const int *const nclass,
                                const int *const ippmod,
                                const int *const icp3pl,
                                const int *const ieqco2,
                                const int *const ihm,
                                const int *const inp,
                                const int *const ixch,
                                const int *const ixck,
                                const int *const ixwt,
                                const int *const ih2,
                                const int *const if1m,
                                const int *const if2m,
                                const int *const if3m,
                                const int *const if4p2m,
                                const int *const if5m,
                                const int *const iyco2)
{
  int i;
  char *name = NULL;
  char *snumsca = NULL;

  assert(vars != NULL);

  if (vars->nscaus > 0) {
    BFT_REALLOC(vars->label, vars->nscapp + vars->nscaus, char*);
  } else {
    BFT_MALLOC(vars->label, vars->nscapp, char*);
  }

  /* IHM */
  BFT_MALLOC(vars->label[*ihm -1], strlen("Enthalpy")+1, char);
  strcpy(vars->label[*ihm -1], "Enthalpy");

  /* Boucle sur les classes IH2, INP, IXCH, IXCK */
  BFT_MALLOC(snumsca, 1 + 2, char);

  /* IH2 */
  BFT_MALLOC(name, strlen("ENT_CP")+1 + 2, char);
  strcpy(name, "ENT_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[ih2[i] -1], strlen(name)+1, char);
    strcpy(vars->label[ih2[i] -1], name);

    strcpy(name, "ENT_CP");
  }

  /* INP */
  BFT_REALLOC(name, strlen("NP_CP")+1 + 2, char);
  strcpy(name, "NP_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[inp[i] -1], strlen(name)+1, char);
    strcpy(vars->label[inp[i] -1], name);

    strcpy(name, "NP_CP");
  }

  /* IXCH */
  BFT_REALLOC(name, strlen("XCH_CP")+1 + 2, char);
  strcpy(name, "XCH_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[ixch[i] -1], strlen(name)+1, char);
    strcpy(vars->label[ixch[i] -1], name);

    strcpy(name, "XCH_CP");
  }

  /* IXCK */
  BFT_REALLOC(name, strlen("XCK_CP")+1 + 2, char);
  strcpy(name, "XCK_CP");
  for (i = 0; i < *nclass; i++) {
    sprintf(snumsca,"%2.2i", i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[ixck[i] -1], strlen(name)+1, char);
    strcpy(vars->label[ixck[i] -1], name);

    strcpy(name, "XCK_CP");
  }

  /* Boucle sur les charbons  IFM1 IFM2 */

  BFT_REALLOC(name, strlen("Fr_MV1")+1 + 2, char);
  strcpy(name, "Fr_MV1");
  for (i = 0; i < *ncharb; i++) {
    sprintf(snumsca,"%2.2i",i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[if1m[i] -1], strlen(name)+1, char);
    strcpy(vars->label[if1m[i] -1], name);

    strcpy(name, "Fr_MV1");
  }

  BFT_REALLOC(name, strlen("Fr_MV2")+1 + 2, char);
  strcpy(name, "Fr_MV2");
  for (i = 0; i < *ncharb; i++) {
    sprintf(snumsca,"%2.2i",i+1);
    strcat(name, snumsca);

    BFT_MALLOC(vars->label[if2m[i] -1], strlen(name)+1, char);
    strcpy(vars->label[if2m[i] -1], name);

    strcpy(name, "Fr_MV2");
  }

 /* IF3M */
  BFT_MALLOC(vars->label[*if3m -1], strlen("Fr_HET")+1, char);
  strcpy(vars->label[*if3m -1], "Fr_HET");

 /* IF4P2M */
  BFT_MALLOC(vars->label[*if4p2m -1], strlen("Var_AIR")+1, char);
  strcpy(vars->label[*if4p2m -1], "Var_AIR");

  if (ippmod[*icp3pl -1] == 1) {
   /* IXWT */
    BFT_MALLOC(name, strlen("XWT_CP")+1 + 2, char);
    strcpy(name, "XWT_CP");
    for (i = 0; i < *nclass; i++) {
      sprintf(snumsca,"%2.2i", i+1);
      strcat(name, snumsca);

      BFT_MALLOC(vars->label[ixwt[i] -1], strlen(name)+1, char);
      strcpy(vars->label[ixwt[i] -1], name);
      strcpy(name, "XWT_CP");
    }

   /* IF5M */
    BFT_MALLOC(vars->label[*if5m -1], strlen("FR_H20")+1, char);
    strcpy(vars->label[*if5m -1], "FR_H20");
  }

  if (*ieqco2 == 1) {
   /* IYCO2 */
    BFT_MALLOC(vars->label[*iyco2 -1], strlen("FR_CO2")+1, char);
    strcpy(vars->label[*iyco2 -1], "FR_CO2");
  }
  BFT_FREE(name);
  BFT_FREE(snumsca);

#if _XML_DEBUG_
  bft_printf(_("==>UICPSC\n"));
  for (i=0; i< vars->nscaus+vars->nscapp; i++) {
    bft_printf(_("--label of scalar[%i]: %s\n"), i, vars->label[i]);
  }
#endif

}

/*----------------------------------------------------------------------------
 * Traitement des conditions aux limites pour le charbon
 *----------------------------------------------------------------------------*/

void CS_PROCF (uicpcl, UICPCL)(const    int *const nozppm,
                               const    int *const ncharm,
                               const    int *const ncharb,
                               const    int *const nclpch,
                               const    int *const nfabor,
                               const    int *const iindef,
                               const    int *const ientre,
                               const    int *const iparoi,
                               const    int *const iparug,
                               const    int *const isymet,
                               const    int *const isolib,
                                        int *const itypfb,
                                        int *const icodcl,
                                     double *const rcodcl,
                                     double *const surfbo,
                                        int *const ientat,
                                        int *const iqimp,
                                     double *const qimpat,
                                     double *const timpat,
                                        int *const ientcp,
                                     double *const qimpcp,
                                     double *const timpcp,
                                     double *const distch,
                                        int *const icalke,
                                     double *const dh,
                                     double *const xintur,
                                        int *const izfppp)
{
  int iphas =0;
  int ivar, ifbr, iwall, c_id;
  int izone, icharb, isca, ith_zone, zones, ifac, zone_nbr;
  int i, k;
  double qimp = 0.;
  double timp = 0.;
  double norm = 0.;
  char *choice = NULL;
  char *label = NULL;
  char *nature = NULL;
  char *description = NULL;
  int *faces_list = NULL;
  int faces = 0;

  assert(vars != NULL);

  zones = cs_gui_boundary_zones_number();

  /* First iteration only : memory allocation */

  if (boundaries == NULL){

    BFT_MALLOC(boundaries,            1,          cs_boundary_t);
    BFT_MALLOC(boundaries->label,     zones,      char*        );
    BFT_MALLOC(boundaries->nature,    zones,      char*        );
    BFT_MALLOC(boundaries->type_code, vars->nvar, int*         );
    BFT_MALLOC(boundaries->values,    vars->nvar, cs_val_t*    );
    BFT_MALLOC(boundaries->ientat,    zones,      int          );
    BFT_MALLOC(boundaries->iqimp,     zones,      int          );
    BFT_MALLOC(boundaries->qimpat,    zones,      double       );
    BFT_MALLOC(boundaries->timpat,    zones,      double       );
    BFT_MALLOC(boundaries->ientcp,    zones,      int          );
    BFT_MALLOC(boundaries->icalke,    zones,      int          );
    BFT_MALLOC(boundaries->dh,        zones,      double       );
    BFT_MALLOC(boundaries->xintur,    zones,      double       );
    BFT_MALLOC(boundaries->qimpcp,    zones,      double*      );
    BFT_MALLOC(boundaries->timpcp,    zones,      double*      );
    BFT_MALLOC(boundaries->distch,    zones,      double**     );
    BFT_MALLOC(boundaries->rough,     zones,      double       );

    boundaries->qimp = NULL;

    for (izone = 0; izone < zones; izone++) {
      BFT_MALLOC(boundaries->qimpcp[izone], *ncharb, double );
      BFT_MALLOC(boundaries->timpcp[izone], *ncharb, double );
      BFT_MALLOC(boundaries->distch[izone], *ncharb, double*);

      for (icharb = 0; icharb < *ncharb; icharb++) {
        BFT_MALLOC(boundaries->distch[izone][icharb],
                   nclpch[icharb],
                   double);
      }
    }

    for (ivar = 0; ivar < vars->nvar; ivar++) {
      i = vars->rtp[ivar];
      BFT_MALLOC(boundaries->type_code[i], zones, int);
      BFT_MALLOC(boundaries->values[i], zones, cs_val_t);
    }

    for (izone = 0; izone < zones; izone++) {
      boundaries->iqimp[izone]  = 0;
      boundaries->ientat[izone] = 0;
      boundaries->ientcp[izone] = 0;
      boundaries->dh[izone]     = 0;
      boundaries->xintur[izone] = 0;
      boundaries->icalke[izone] = 0;
      boundaries->qimpat[izone] = 0;
      boundaries->timpat[izone] = 0;
      boundaries->rough[izone]  = -999;

      for (icharb = 0; icharb < *ncharb; icharb++) {
        boundaries->qimpcp[izone][icharb] = 0;
        boundaries->timpcp[izone][icharb] = 0;

        for (k = 0; k < nclpch[icharb]; k++)
          boundaries->distch[izone][icharb][k] = 0;
      }
    }

    /* Initialization of boundary->type_code and boundary->values */

    for (ivar = 0; ivar < vars->nvar; ivar++) {
      i = vars->rtp[ivar];
      for (izone = 0; izone < zones; izone++) {
        boundaries->type_code[i][izone] = -1;
        boundaries->values[i][izone].val1 = 1.e30;
        boundaries->values[i][izone].val2 = 1.e30;
        boundaries->values[i][izone].val3 = 0.;
      }
    }

    for (izone = 0; izone < zones; izone++) {

   /* nature, label of the ith initialization zone */

      ith_zone = izone + 1;
      nature = cs_gui_boundary_zone_nature(ith_zone);
      label = cs_gui_boundary_zone_label(ith_zone);

      BFT_MALLOC(boundaries->label[izone], strlen(label)+1, char);
      strcpy(boundaries->label[izone], label);

      BFT_MALLOC(boundaries->nature[izone], strlen(nature)+1, char);
      strcpy(boundaries->nature[izone], nature);

      if (cs_gui_strcmp(nature, "inlet")) {

        /* INLET: VELOCITY */
        choice = cs_gui_boundary_choice("inlet", label, "velocity_pressure");

        if (cs_gui_strcmp(choice, "coal_flow")) {

          boundaries->ientcp[izone] = 1;
          boundaries->iqimp[izone]  = 1;
          cs_gui_coal_boundary_coalflow(izone, ncharb, nclpch);
          cs_gui_boundary_flow(label,&qimp,&timp);
          boundaries->qimpat[izone] = qimp;
          boundaries->timpat[izone] = timp;

        } else if (cs_gui_strcmp(choice, "flow1")) {

          boundaries->ientat[izone] = 1;
          boundaries->iqimp[izone]  = 1;
          cs_gui_boundary_flow(label,&qimp,&timp);
          boundaries->qimpat[izone] = qimp;
          boundaries->timpat[izone] = timp;
          /* TODO : remplir la direction normale a la face */
          /* boundaries->values[1][izone].val1 = directionU ; */
          /* boundaries->values[2][izone].val1 = directionv ; */
          /* boundaries->values[3][izone].val1 = directionw ; */
        }

        BFT_FREE(choice);

        /* INLET: TURBULENCE */
        choice = cs_gui_boundary_choice("inlet", label, "turbulence");
        cs_gui_boundary_turbulence(choice, izone);
        BFT_FREE(choice);

        /* INLET: USER SCALARS */
        for (isca=0 ; isca < vars->nscaus ; isca++) {
          cs_gui_boundary_value_scalar("inlet", izone, isca);
        }

      } else if (cs_gui_strcmp(nature, "wall")) {

        /* WALL: VELOCITY */
        choice = cs_gui_boundary_choice("wall", label, "velocity_pressure");

        if (cs_gui_strcmp(choice, "on")) {
          for (ivar = 1; ivar < 4; ivar++)
            cs_gui_boundary_dirichlet("wall", label, izone, ivar);
        }
        BFT_FREE(choice);

        /* Wall: ROUGH */
        cs_gui_boundary_rough(label, izone);

        /* WALL: USER SCALARS */
        for (isca = 0; isca < vars->nscaus; isca++) {
          cs_gui_boundary_value_scalar("wall", izone, isca);
        }

      } else if (cs_gui_strcmp(nature, "outlet")) {

        /* OUTLET: USER SCALARS */
        for (isca=0; isca < vars->nscaus ; isca++) {
          cs_gui_boundary_value_scalar("outlet", izone, isca);
        }
      } /* if (cs_gui_strcmp(nature, "outlet")) */

      BFT_FREE(nature);
      BFT_FREE(label);

    } /* for izones */

  }  /* if (boundaries == NULL)*/

 /* A chaque iteration, boucle sur les faces de bord :
    on remplit itypfb, rcodcl et icodcl a partir des tableaux
    de la structures conditions.limites definie
    dans la premiere partie de la fonction.
    Remember: rdoccl[k][j][i] = rcodcl[ k * dim1 *dim2 + j *dim1 + i] */

  for (izone=0 ; izone < zones ; izone++) {

    ith_zone = izone + 1;
    zone_nbr = cs_gui_boundary_zone_number(ith_zone);
    if (zone_nbr > *nozppm)
      bft_error(__FILE__, __LINE__, 0,
                _("zone's label number %i is greater than %i the maximum allowed \n"),
                 zone_nbr, *nozppm);

    description = cs_gui_boundary_zone_localization(boundaries->nature[izone],
                                                    boundaries->label[izone]);

    /* build list of faces */
    BFT_MALLOC(faces_list, *nfabor, int);

    c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                                 description,
                                 &faces,
                                 faces_list);

    if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
      const char *missing
        = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
      bft_error(__FILE__, __LINE__, 0,
                _("Le groupe ou attribut \"%s\" figurant dans le\n"
                  "critère de sélection:\n"
                  "\"%s\"\n ne correspond à aucune face de bord."),
                missing, description);
    }

    BFT_FREE(description);

    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {

      /* Update the depending zone arrays (iqimp, dh, xintur, icalke,...)
         because they are initialized at each time step in PPPRCL routine */

      iqimp[zone_nbr-1]  = boundaries->iqimp[izone];
      ientat[zone_nbr-1] = boundaries->ientat[izone];
      ientcp[zone_nbr-1] = boundaries->ientcp[izone];
      dh[zone_nbr-1]     = boundaries->dh[izone];
      xintur[zone_nbr-1] = boundaries->xintur[izone];
      icalke[zone_nbr-1] = boundaries->icalke[izone];
      qimpat[zone_nbr-1] = boundaries->qimpat[izone];
      timpat[zone_nbr-1] = boundaries->timpat[izone];

      for (icharb = 0; icharb < *ncharb; icharb++) {
        qimpcp[icharb *(*nozppm)+zone_nbr-1] =boundaries->qimpcp[izone][icharb];
        timpcp[icharb *(*nozppm)+zone_nbr-1] =boundaries->timpcp[izone][icharb];

        for (k = 0; k < nclpch[icharb]; k++) {
          distch[k * (*nozppm) * (*ncharm) +icharb * (*nozppm) +zone_nbr-1]
          = boundaries->distch[izone][icharb][k];
        }
      }

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = *ientre;
        for (i = 0; i < vars->nvar; i++) {
          ivar = vars->rtp[i];
          rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
          = boundaries->values[ivar][izone].val1 ;
        }
        norm = 1.0 / ( sqrt( surfbo[3 * ifbr + 0] * surfbo[3 * ifbr + 0]
                           + surfbo[3 * ifbr + 1] * surfbo[3 * ifbr + 1]
                           + surfbo[3 * ifbr + 2] * surfbo[3 * ifbr + 2] ) );
        rcodcl[vars->rtp[1] * (*nfabor) + ifbr] = -surfbo[3 * ifbr + 0]*norm;
        rcodcl[vars->rtp[2] * (*nfabor) + ifbr] = -surfbo[3 * ifbr + 1]*norm;
        rcodcl[vars->rtp[3] * (*nfabor) + ifbr] = -surfbo[3 * ifbr + 2]*norm;
      }

    } else if (cs_gui_strcmp(boundaries->nature[izone], "wall")) {

      if (boundaries->rough[izone] >= 0.) {
        iwall = *iparug;
        /* roughness value is only stored in Velocity_U */
        ivar = 1;
        for (ifac = 0; ifac < faces; ifac++) {
          ifbr = faces_list[ifac]-1;
          icodcl[ivar *(*nfabor) + ifbr] = 6;
          rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
          = boundaries->rough[izone];
        }
      } else {
        iwall = *iparoi;
      }

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = iwall;
      }

      for (i = 0; i < vars->nvar; i++) {
        ivar = vars->rtp[i];

        switch (boundaries->type_code[ivar][izone]) {

          case NEUMANN :

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]-1;
              icodcl[ivar *(*nfabor) + ifbr ] = 3 ;
              rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val3 ;
            }
          break;

          case DIRICHLET :

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]-1;
              icodcl[ivar *(*nfabor) + ifbr ] = 5 ;
              rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1 ;
            }
          break;

          case COEF_ECHANGE :

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]-1;
              icodcl[ivar *(*nfabor) + ifbr ] = 5 ;
              rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1 ;
              rcodcl[1 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val2 ;
            }
          break;
        }
      }

    } else if (cs_gui_strcmp(boundaries->nature[izone], "outlet")) {

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = *isolib;
      }

      for (i = 0; i < vars->nvar; i++) {
        ivar = vars->rtp[i];

        switch (boundaries->type_code[ivar][izone]) {

          case DIRICHLET :

            for (ifac = 0; ifac < faces; ifac++) {
              ifbr = faces_list[ifac]-1;
              icodcl[ivar *(*nfabor) + ifbr ] = 5 ;
              rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]
              = boundaries->values[ivar][izone].val1 ;
            }
          break;

        }
      }

    } else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = *isymet;
      }

    } else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {

      for (ifac = 0; ifac < faces; ifac++) {
        ifbr = faces_list[ifac]-1;
        izfppp[ifbr] = zone_nbr;
        itypfb[iphas *(*nfabor) +ifbr] = *iindef;
      }

    } else
        bft_error(__FILE__, __LINE__, 0,
                  _("boundary nature %s is unknown \n"),
                    boundaries->nature[izone]);

    BFT_FREE(faces_list);
  } /*  for izone */

#if _XML_DEBUG_
  bft_printf(_("==>UICPCL\n"));
  bft_printf(_("--boundary zones number: %i\n"), zones);

  for (izone=0 ; izone < zones ; izone++) {

    BFT_MALLOC(faces_list, *nfabor, int);
    description = cs_gui_boundary_zone_localization(boundaries->nature[izone],
                                                    boundaries->label[izone]);
    c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                                 description,
                                 &faces,
                                 faces_list);

    if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
      const char *missing
        = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
      bft_error(__FILE__, __LINE__, 0,
                _("Le groupe ou attribut \"%s\" figurant dans le\n"
                  "critère de sélection:\n"
                  "\"%s\"\n ne correspond à aucune face de bord."),
                missing, description);
    }

    zone_nbr = cs_gui_boundary_zone_number(izone+1);

    bft_printf(_("---zone %i label: %s\n"), zone_nbr, boundaries->label[izone]);
    bft_printf(_("---zone %i nature: %s\n"), zone_nbr, boundaries->nature[izone]);
    bft_printf(_("---zone %i number of faces: %i\n"), zone_nbr, faces);
    bft_printf(_("----localization: %s\n"), description);
    BFT_FREE(description);

    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {
      bft_printf(_("-----iqimp=%i, qimpat=%12.5e \n"),
                   iqimp[zone_nbr-1], qimpat[zone_nbr-1]);
      bft_printf(_("-----icalke=%i, dh=%12.5e, xintur=%12.5e \n"),
                   icalke[zone_nbr-1], dh[zone_nbr-1], xintur[zone_nbr-1]);
      bft_printf(_("-----ientat=%i, ientcp=%i, timpat=%12.5e \n"),
                   ientat[zone_nbr-1], ientcp[zone_nbr-1], timpat[zone_nbr-1]);

      for (icharb = 0; icharb < *ncharb; icharb++) {
        bft_printf(_("-----coal=%i, qimpcp=%12.5e, timpcp=%12.5e \n"),
                      icharb, qimpcp[icharb *(*nozppm)+zone_nbr-1],
                      timpcp[icharb *(*nozppm)+zone_nbr-1]);

        for (k = 0; k < nclpch[icharb]; k++)
          bft_printf(_("-----coal=%i, class=%i, distch=%f \n"),
                       icharb, k,
                       distch[k * (*nozppm) * (*ncharm) +icharb * (*nozppm) +zone_nbr-1]);
      }
    }

    if (faces>0) {
        ifbr = faces_list[0]-1;
        for (i=0; i<vars->nvar - vars->nscaus - vars->nscapp ; i++) {
          ivar = vars->rtp[i];
          bft_printf(_("-----%s: itypfb=%i, icodcl=%i, rcodcl(1)=%12.5e, rcodcl(2)=%12.5e, rcodcl(3)=%12.5e\n"),
                         vars->name[ivar],
                         itypfb[iphas *(*nfabor) + ifbr],
                         icodcl[ivar *(*nfabor) + ifbr ],
                         rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr],
                         rcodcl[1 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr],
                         rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]);
        }

        for (ivar=0; ivar<vars->nscaus + vars->nscapp ; ivar++) {
          bft_printf(_("-----%s: itypfb=%i, icodcl=%i, rcodcl(1)=%12.5e, rcodcl(2)=%12.5e, rcodcl(3)=%12.5e\n"),
                         vars->label[ivar],
                         itypfb[iphas *(*nfabor) + ifbr],
                         icodcl[ivar *(*nfabor) + ifbr ],
                         rcodcl[0 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr],
                         rcodcl[1 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr],
                         rcodcl[2 * (*nfabor * (vars->nvar)) + ivar * (*nfabor) + ifbr]);
        }
    }
    BFT_FREE(faces_list);

  }
#endif
}

/*----------------------------------------------------------------------------
 * Boundary conditions input verification
 *
 * Fortran Interface:
 *
 * SUBROUTINE UICLVE
 * *****************
 *
 * INTEGER          NFABOR  --> number of boundary faces
 * INTEGER          IINDEF  --> type of boundary: not defined
 * INTEGER          IENTRE  --> type of boundary: inlet
 * INTEGER          IPAROI  --> type of boundary: smooth wall
 * INTEGER          IPARUG  --> type of boundary: rough wall
 * INTEGER          ISYMET  --> type of boundary: symetry
 * INTEGER          ISOLIB  --> type of boundary: outlet
 * INTEGER          ITYPFB  --> type of boundary for each face
 * INTEGER          IZFPPP  --> zone number
 *----------------------------------------------------------------------------*/

void CS_PROCF (uiclve, UICLVE)(const int *const nfabor,
                               const int *const iindef,
                               const int *const ientre,
                               const int *const iparoi,
                               const int *const iparug,
                               const int *const isymet,
                               const int *const isolib,
                                     int *const itypfb,
                                     int *const izfppp)
{
  int ifbr, ifac, c_id;
  int izone, zones, zone_nbr;
  int inature, inature2;
  int *faces_list = NULL;
  int faces = 0, iphas = 0;
  char *description = NULL;

  zones   = cs_gui_boundary_zones_number();

  for (izone=0 ; izone < zones ; izone++) {

    zone_nbr =  cs_gui_boundary_zone_number(izone + 1);

    description = cs_gui_boundary_zone_localization(boundaries->nature[izone],
                                                    boundaries->label[izone]);

    /* build list of faces */
    BFT_MALLOC(faces_list, *nfabor, int);

    c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                                 description,
                                 &faces,
                                 faces_list);

    if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
      const char *missing
        = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
      bft_error(__FILE__, __LINE__, 0,
                _("Le groupe ou attribut \"%s\" figurant dans le\n"
                  "critère de sélection:\n"
                  "\"%s\"\n ne correspond à aucune face de bord."),
                missing, description);
    }

    BFT_FREE(description);

    if (cs_gui_strcmp(boundaries->nature[izone], "inlet")) {

      inature = *ientre;

    } else if (cs_gui_strcmp(boundaries->nature[izone], "wall")) {

      inature = *iparug;
      if (boundaries->rough[izone] <0.){
        inature = *iparoi;
      }

    } else if (cs_gui_strcmp(boundaries->nature[izone], "outlet")) {

      inature = *isolib;

    } else if (cs_gui_strcmp(boundaries->nature[izone], "symmetry")) {

      inature = *isymet;

    } else if (cs_gui_strcmp(boundaries->nature[izone], "undefined")) {

      inature = *iindef;

    } else
        bft_error(__FILE__, __LINE__, 0,
                  _("boundary nature %s is unknown \n"),
                     boundaries->nature[izone]);

    for (ifac = 0; ifac < faces; ifac++) {
      ifbr = faces_list[ifac]-1;

      if (izfppp[ifbr] != zone_nbr)
        bft_error(__FILE__, __LINE__, 0,
        _("@                                                            \n"
          "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
          "@                                                            \n"
          "@ @@ WARNING: BOUNDARY CONDITIONS ERROR                      \n"
          "@    *******                                                 \n"
          "@                                                            \n"
          "@    The zone %s does not have the same id number            \n"
          "@    in the GUI and in the user subroutine.                  \n"
          "@                                                            \n"
          "@    GUI zone number:             %i                         \n"
          "@    USER SUBROUTINE zone number: %i                         \n"
          "@                                                            \n"
          "@    The id number given in the GUI cannot be modified       \n"
          "@    in the user subroutine (fortran array IZFPPP).          \n"
          "@                                                            \n"
          "@    The calculation will stop.                              \n"
          "@                                                            \n"
          "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
          "@                                                            \n"),
          boundaries->label[izone], zone_nbr, izfppp[ifbr]);

      inature2 = itypfb[iphas *(*nfabor) +ifbr];

        /* The nature of the boundary can be changed from smooth wall to
           rough wall or vice-versa between the GUI and the FORTRAN */
        if (inature2 == *iparug ) inature2 = *iparoi;
        if (inature == *iparug ) inature = *iparoi;

      if (inature2 != inature)
        bft_error(__FILE__, __LINE__, 0,
        _("@                                                            \n"
          "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
          "@                                                            \n"
          "@ @@ WARNING: BOUNDARY CONDITIONS ERROR                      \n"
          "@    *******                                                 \n"
          "@                                                            \n"
          "@    The zone %s does not have the same nature               \n"
          "@    in the GUI and in the user subroutine.                  \n"
          "@                                                            \n"
          "@    GUI zone nature:             %s                         \n"
          "@    USER SUBROUTINE ITYPFB:      %i                         \n"
          "@                                                            \n"
          "@    The nature given in the GUI cannot be modified          \n"
          "@    in the user subroutine (fortran array ITYPFB).          \n"
          "@                                                            \n"
          "@    The calculation will stop.                              \n"
          "@                                                            \n"
          "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
          "@                                                            \n"),
          boundaries->label[izone], boundaries->nature[izone], inature2);
    }

    BFT_FREE(description);
    BFT_FREE(faces_list);

  } /*  for izone */

}

/*----------------------------------------------------------------------------
 * Free memory: clean global private variables and libxml2 variables
 *
 * Fortran Interface:
 *
 * SUBROUTINE MEMUI1
 * *****************
 *
 * INTEGER          NCHARB  --> number of coal
 *----------------------------------------------------------------------------*/

void CS_PROCF (memui1, MEMUI1) (const int *const ncharb)
{
  int i;
  int ivar;
  int izone;
  int zones;
  int icharb;

  assert(vars != NULL);

  if (boundaries != NULL) {

  /* clean memory for global private structure boundaries */

    zones = cs_gui_boundary_zones_number();
    for (izone=0 ; izone < zones ; izone++) {
      BFT_FREE(boundaries->label[izone]);
      BFT_FREE(boundaries->nature[izone]);
    }

    for (i=0; i < vars->nvar; i++) {
      ivar = vars->rtp[i];
      BFT_FREE(boundaries->type_code[ivar]);
      BFT_FREE(boundaries->values[ivar]);
    }

    if (cs_gui_strcmp(vars->model, "pulverized_coal")) {
      for (izone=0 ; izone < zones ; izone++) {
        BFT_FREE(boundaries->qimpcp[izone]);
        BFT_FREE(boundaries->timpcp[izone]);
        for (icharb=0; icharb < *ncharb; icharb++)
          BFT_FREE(boundaries->distch[izone][icharb]);
        BFT_FREE(boundaries->distch[izone]);
      }
      BFT_FREE(boundaries->ientat);
      BFT_FREE(boundaries->ientcp);
      BFT_FREE(boundaries->qimpat);
      BFT_FREE(boundaries->timpat);
      BFT_FREE(boundaries->qimpcp);
      BFT_FREE(boundaries->timpcp);
      BFT_FREE(boundaries->distch);
    }

    BFT_FREE(boundaries->label);
    BFT_FREE(boundaries->nature);
    BFT_FREE(boundaries->iqimp);
    BFT_FREE(boundaries->icalke);
    BFT_FREE(boundaries->qimp);
    BFT_FREE(boundaries->dh);
    BFT_FREE(boundaries->xintur);
    BFT_FREE(boundaries->type_code);
    BFT_FREE(boundaries->values);
    BFT_FREE(boundaries);
  }

  if (vars != NULL) {

  /* clean memory for global private structure vars */

    for (i=0; i < vars->nvar; i++) {
      BFT_FREE(vars->type[i]);
      BFT_FREE(vars->name[i]);
    }
    for (i=0; i < vars->nscaus+vars->nscapp; i++)
      BFT_FREE(vars->label[i]);
    for (i=0; i < vars->nprop; i++)
      BFT_FREE(vars->properties_name[i]);
    BFT_FREE(vars->label);
    BFT_FREE(vars->model);
    BFT_FREE(vars->model_value);
    BFT_FREE(vars->rtp);
    BFT_FREE(vars->name);
    BFT_FREE(vars->properties_name);
    BFT_FREE(vars->properties_ipp);
    BFT_FREE(vars);
  }

  /* clean memory for fortran name of variables */

  for (i = 0; i < _cs_gui_max_vars; i++)
    BFT_FREE(_cs_gui_var_name[i]);
  BFT_FREE(_cs_gui_var_name);

  /* clean memory for xml document */

  if (xpathCtx != NULL) xmlXPathFreeContext(xpathCtx);
  if (node != NULL) xmlFreeNode(node);

  /* Shutdown libxml */

  xmlCleanupParser();
  xmlMemoryDump();
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _CS_HAVE_XML */
