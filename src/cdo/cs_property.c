/*============================================================================
 * Manage the definition/setting of properties
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include <assert.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>
#include "cs_mesh_location.h"
#include "cs_reco.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_property.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_PROPERTY_DBG  1

typedef struct {

  cs_param_def_type_t    def_type;  // by value, by analytic function...
  cs_def_t               def;       // accessor to the definition

  const void            *struc;     /* if this property related to a subdomain
                                       hinges on a structure. Only shared */

} cs_property_subdef_t;

/* Set of parameters attached to a property */
struct _cs_property_t {

  char  *restrict name;

  cs_flag_t    flag;      /* Short descriptor (mask of bits) */
  int          post_freq; /* -1: no post, 0: at the beginning otherwise at each
                             post_freq iteration.
                             If post_freq > -1, a related cs_field_t structure
                             is created. */

  /* Pointer to the main structures (not owned, only shared) */
  const cs_cdo_quantities_t   *cdoq;
  const cs_cdo_connect_t      *connect;
  const cs_time_step_t        *time_step;

  /* The number of values to set depends on the type of property
       - isotropic   = 1 => CS_PARAM_VAR_SCAL
       - orthotropic = 3 => CS_PARAM_VAR_VECT
       - anisotropic = 9 => CS_PARAM_VAR_TENS
  */
  cs_property_type_t   type;     // isotropic, anistotropic...

  /* How the property is defined */
  cs_param_def_type_t  def_type; // by value, by analytic function...
  cs_def_t             def;      // accessor to the definition

  /* Useful buffers to deal with more complex definitions */
  cs_flag_t    array_flag;  // short description of the related array
  cs_real_t   *array;       // if the property hinges on an array

  const void  *struc;       // if one needs a structure. Only shared

  /* Specific members for a definition by subdomains (up to now, only subdomain
     built from an union of cells is considered) */
  int                    n_subdomains;
  cs_property_subdef_t  *sub_defs;   // list of definitions for each subdomain
  short int             *def_ids;    // id in sub_defs for each cell

};

/* List of available keys for setting a property */
typedef enum {

  PTYKEY_POST_FREQ,
  PTYKEY_ERROR

} ptykey_t;

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_pty[] =
  " Stop setting an empty cs_property_t structure.\n"
  " Please check your settings.\n";

static const cs_flag_t  cs_var_support_pc =
  CS_PARAM_FLAG_PRIMAL | CS_PARAM_FLAG_CELL;
static const cs_flag_t  cs_var_support_pv =
  CS_PARAM_FLAG_PRIMAL | CS_PARAM_FLAG_VERTEX;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the name of the corresponding property key
 *
 * \param[in] key        name of the key
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

static const char *
_print_ptykey(ptykey_t  key)
{
  switch (key) {
  case PTYKEY_POST_FREQ:
    return "post_freq";

  default:
    assert(0);
  }

  return NULL; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the corresponding enum from the name of a property key.
 *         If not found, return a key error.
 *
 * \param[in] keyname    name of the key
 *
 * \return a ptykey_t
 */
/*----------------------------------------------------------------------------*/

static ptykey_t
_get_ptykey(const char  *keyname)
{
  ptykey_t  key = PTYKEY_ERROR;

  if (strcmp(keyname, "post_freq") == 0)
    key = PTYKEY_POST_FREQ;

  return key;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the settings are valid
 *
 * \param[in]  pty       pointer to a cs_property_t structure
 * \param[in]  get       accessor to the tensor values
 */
/*----------------------------------------------------------------------------*/

static void
_check_tensor_symmetry(const cs_property_t    *pty,
                       cs_get_t                get)
{
  if ((get.tens[0][1] - get.tens[1][0]) > cs_get_zero_threshold() ||
      (get.tens[0][2] - get.tens[2][0]) > cs_get_zero_threshold() ||
      (get.tens[1][2] - get.tens[2][1]) > cs_get_zero_threshold())
    bft_error(__FILE__, __LINE__, 0,
              _(" The definition of the tensor related to the"
                " property %s is not symmetric.\n"
                " This case is not handled. Please check your settings.\n"),
              pty->name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the settings are valid
 *
 * \param[in]  pty       pointer to a cs_property_t structure
 * \param[in]  ml_name   name of the related mesh location
 * \param[in]  ml_id     id of the related mesh location
 */
/*----------------------------------------------------------------------------*/

static void
_check_def_by_subdomain(const cs_property_t    *pty,
                        const char             *ml_name,
                        int                     ml_id)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  if (pty->n_subdomains == 0 && pty->def_type != CS_PARAM_N_DEF_TYPES)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid definition by subdomain for property %s.\n"
              " Another type of definition has already been used.", pty->name);

  if (pty->n_subdomains > 0 && pty->def_type != CS_PARAM_DEF_BY_SUBDOMAIN)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid definition by subdomain for property %s.\n"
              " Another type of definition has already been used.", pty->name);

  if (ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" mesh location %s has not been found.\n"
                " Please check your settings."), ml_name);

  if (cs_mesh_location_get_type(ml_id) != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of mesh location for mesh location  %s.\n"
                " The expected type is CS_MESH_LOCATION_CELLS.\n"), ml_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new definition to a cs_property_t structure defined by domain
 *
 * \param[in]  pty       pointer to a cs_property_t structure
 * \param[in]  new_id    id of the new definition to set
 * \param[in]  ml_id     id of the related mesh location
 */
/*----------------------------------------------------------------------------*/

static void
_add_def_by_subdomain(cs_property_t    *pty,
                      int               new_id,
                      int               ml_id)
{
  cs_lnum_t  i;

  const cs_lnum_t  n_cells = pty->cdoq->n_cells;

  pty->def_type = CS_PARAM_DEF_BY_SUBDOMAIN;

  if (new_id == 0) { /* Initialization */
    BFT_MALLOC(pty->def_ids, n_cells, short int);
    for (i = 0; i < n_cells; i++)
      pty->def_ids[i] = -1; // unset by default
  }

  /* Assign new id to the selected cells */
  const cs_lnum_t *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  if (elt_ids == NULL) {
    assert(n_elts[0] == n_cells); /* Sanity check. Up to now only this
                                     case is considered */
    for (i = 0; i < n_elts[0]; i++)
      pty->def_ids[i] = new_id;

  }
  else
    for (i = 0; i < n_elts[0]; i++)
      pty->def_ids[elt_ids[i]] = new_id;

  pty->n_subdomains += 1;
  BFT_REALLOC(pty->sub_defs, pty->n_subdomains, cs_property_subdef_t);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value using a law with one argument
 *
 * \param[in]       pty     pointer to a cs_property_t structure
 * \param[in]       get     accessor to the value
 * \param[in, out]  tensor  result stored in a 3x3 tensor
 */
/*----------------------------------------------------------------------------*/

static void
_get_tensor_by_value(const cs_property_t      *pty,
                     cs_get_t                  get,
                     cs_real_3_t              *tensor)
{
  int  k, l;

  switch (pty->type) {

  case CS_PROPERTY_ISO:
    tensor[0][0] = tensor[1][1] = tensor[2][2] = get.val;
    break;

  case CS_PROPERTY_ORTHO:
    for (k = 0; k < 3; k++)
      tensor[k][k] = get.vect[k];
    break;

  case CS_PROPERTY_ANISO:
    for (k = 0; k < 3; k++)
      for (l = 0; l < 3; l++)
        tensor[k][l] = get.tens[k][l];
    break;

  default:
    assert(0);
    break;

  } // Property type
}
 
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value using a law with one argument
 *
 * \param[in]       c_id    cell id
 * \param[in]       pty     pointer to a cs_property_t structure
 * \param[in]       law     function pointer to the law
 * \param[in]       struc   pointer to a structure
 * \param[in, out]  get     pointer to a union used to retrieve the result
 */
/*----------------------------------------------------------------------------*/

static void
_get_result_by_onevar_law(cs_lnum_t                 c_id,
                          const cs_property_t      *pty,
                          cs_onevar_law_func_t     *law,
                          const void               *struc,
                          cs_get_t                 *get)
{
  assert(pty->array != NULL); /* Sanity check */

  /* Test if flag has at least the pattern of the reference support */
  if ((pty->array_flag & cs_var_support_pc) == cs_var_support_pc)
    law(pty->array[c_id], struc, get);

  /* Test if flag has at least the pattern of the reference support */
  else if ((pty->array_flag & cs_var_support_pv) == cs_var_support_pv) {

    cs_real_t  val_xc;

    /* Reconstruct (or interpolate) value at the current cell center */
    cs_reco_pv_at_cell_center(c_id,
                              pty->connect->c2v,
                              pty->cdoq,
                              pty->array, &val_xc);
    
    law(val_xc, struc, get);

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " Invalid support for evaluating the property %s"
              " by law with one argument.", pty->name);

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new property structure
 *
 * \param[in]  name        name of the property
 * \param[in]  key_type    keyname of the type of property
 * \param[in]  cdoq        pointer to a cs_cdo_quantities_t struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_create(const char                  *name,
                   const char                  *key_type,
                   const cs_cdo_quantities_t   *cdoq,
                   const cs_cdo_connect_t      *connect,
                   const cs_time_step_t        *time_step)
{
  cs_property_t  *pty = NULL;

  BFT_MALLOC(pty, 1, cs_property_t);

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(pty->name, len, char);
  strncpy(pty->name, name, len);

  /* Shared pointers for defining the property */
  pty->cdoq = cdoq;
  pty->connect = connect;
  pty->time_step = time_step;

  /* Assign a type */
  if (strcmp(key_type, "isotropic") == 0)
    pty->type = CS_PROPERTY_ISO;
  else if (strcmp(key_type, "orthotropic") == 0)
    pty->type = CS_PROPERTY_ORTHO;
  else if (strcmp(key_type, "anisotropic") == 0)
    pty->type = CS_PROPERTY_ANISO;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key %s for setting the type of property.\n"
                " Key is one of the following: isotropic, orthotropic or"
                " anisotropic.\n"
                " Please modify your settings."), key_type);

  /* Default initialization */
  pty->post_freq = -1;
  pty->flag = 0;
  pty->def_type = CS_PARAM_N_DEF_TYPES;
  pty->def.get.val = 0;

  /* Specific members for more complex definition */
  pty->array_flag = 0;
  pty->array = NULL;

  /* Specific members for a definition by subdomain */
  pty->n_subdomains = 0;
  pty->sub_defs = NULL;
  pty->def_ids = NULL;

  //  pty->struc is only shared if needed

  return pty;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_property_t structure
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_free(cs_property_t   *pty)
{
  if (pty == NULL)
    return pty;

  BFT_FREE(pty->name);

  if (pty->def_type == CS_PARAM_DEF_BY_SUBDOMAIN) {
    BFT_FREE(pty->sub_defs);
    BFT_FREE(pty->def_ids);
  }

  if (pty->array_flag & CS_PARAM_FLAG_OWNER) {
    if (pty->array != NULL)
      BFT_FREE(pty->array);
  }

  /* All other pointers are shared */
  BFT_FREE(pty);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the given property has the name ref_name
 *
 * \param[in]  pty         pointer to a cs_property_t structure to test
 * \param[in]  ref_name    name of the property to find
 *
 * \return true if the name of the property is ref_name otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_property_check_name(const cs_property_t   *pty,
                       const char            *ref_name)
{
  if (pty == NULL)
    return false;

  int  reflen = strlen(ref_name);
  int  len = strlen(pty->name);

  if (reflen == len) {
    if (strcmp(ref_name, pty->name) == 0)
      return true;
    else
      return false;
  }
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  returns true if the property is uniform, otherwise false
 *
 * \param[in]    pty    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_property_is_uniform(const cs_property_t   *pty)
{
  if (pty == NULL)
    return false;

  if (pty->flag & CS_PARAM_FLAG_UNIFORM)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of a property
 *
 * \param[in]    pty    pointer to a property
 *
 * \return  the name of the related property
 */
/*----------------------------------------------------------------------------*/

const char *
cs_property_get_name(const cs_property_t   *pty)
{
  if (pty == NULL)
    return NULL;

  return pty->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the type of a property
 *
 * \param[in]    pty    pointer to a property
 *
 * \return  the type of the related property
 */
/*----------------------------------------------------------------------------*/

cs_property_type_t
cs_property_get_type(const cs_property_t   *pty)
{
  if (pty == NULL)
    return CS_PROPERTY_N_TYPES;

  return pty->type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a summary of a cs_property_t structure
 *
 * \param[in]  pty      pointer to a cs_property_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_property_summary(const cs_property_t   *pty)
{
  if (pty == NULL)
    return;

  _Bool  is_uniform = false, is_steady = true;

  if (pty->flag & CS_PARAM_FLAG_UNIFORM)  is_uniform = true;
  if (pty->flag & CS_PARAM_FLAG_UNSTEADY) is_steady = false;

  bft_printf(" %s >> uniform [%s], steady [%s], ",
             pty->name, cs_base_strtf(is_uniform), cs_base_strtf(is_steady));

  switch(pty->type) {
  case CS_PROPERTY_ISO:
    bft_printf("type: isotropic\n");
    break;
  case CS_PROPERTY_ORTHO:
    bft_printf("type: orthotropic\n");
    break;
  case CS_PROPERTY_ANISO:
    bft_printf("type: anisotropic\n");
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of property."));
    break;
  }

  switch (pty->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    {
      const cs_get_t  mat = pty->def.get;

      switch(pty->type) {

      case CS_PROPERTY_ISO:
        bft_printf("       definition by value: % 5.3e\n", mat.val);
        break;
      case CS_PROPERTY_ORTHO:
        bft_printf("       definition by value: (% 5.3e, % 5.3e, % 5.3e)\n",
                   mat.vect[0], mat.vect[1], mat.vect[2]);
        break;
      case CS_PROPERTY_ANISO:
        bft_printf("                            |% 5.3e, % 5.3e, % 5.3e|\n"
                   "       definition by value: |% 5.3e, % 5.3e, % 5.3e|\n"
                   "                            |% 5.3e, % 5.3e, % 5.3e|\n",
                   mat.tens[0][0], mat.tens[0][1], mat.tens[0][2],
                   mat.tens[1][0], mat.tens[1][1], mat.tens[1][2],
                   mat.tens[2][0], mat.tens[2][1], mat.tens[2][2]);
        break;
      default:
        break;

      } // pty->type
    }
    break; // BY_VALUE

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    bft_printf("       definition by an analytical function\n");
    break;

  case CS_PARAM_DEF_BY_ONEVAR_LAW:
    bft_printf("       definition by a law depending on one variable\n");
    break;

  case CS_PARAM_DEF_BY_SUBDOMAIN:
    bft_printf("       definition by %d subdomains\n", pty->n_subdomains);

    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of definition for a property."));
    break;

  } /* switch on def_type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the value of an isotropic property
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       val      value of the property
 */
/*----------------------------------------------------------------------------*/

void
cs_property_iso_def_by_value(cs_property_t    *pty,
                             double            val)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if (pty->type != CS_PROPERTY_ISO)
    bft_error(__FILE__, __LINE__, 0,
              " One wants to set an isotropic property but %s property is not"
              " isotropic", pty->name);

  pty->def_type = CS_PARAM_DEF_BY_VALUE;
  pty->flag |= CS_PARAM_FLAG_UNIFORM;
  pty->def.get.val = val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the values of an orthotropic property
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       val      values of the property
 */
/*----------------------------------------------------------------------------*/

void
cs_property_ortho_def_by_value(cs_property_t    *pty,
                               const double      val[])
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if (pty->type != CS_PROPERTY_ORTHO)
    bft_error(__FILE__, __LINE__, 0,
              " One wants to set an orthotropic property but %s property is not"
              " orthotropic", pty->name);

  pty->def_type = CS_PARAM_DEF_BY_VALUE;
  pty->flag |= CS_PARAM_FLAG_UNIFORM;

  pty->def.get.vect[0] = val[0];
  pty->def.get.vect[1] = val[1];
  pty->def.get.vect[2] = val[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the values of an anisotropic property
 *
 * \param[in, out]  pty      pointer to a cs_property_t structure
 * \param[in]       tens     3x3 matrix with the values of the property
 */
/*----------------------------------------------------------------------------*/

void
cs_property_aniso_def_by_value(cs_property_t    *pty,
                               const double      tens[3][3])
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));
  if (pty->type != CS_PROPERTY_ANISO)
    bft_error(__FILE__, __LINE__, 0,
              " One wants to set an anisotropic property but %s property is not"
              " anisotropic", pty->name);

  pty->def_type = CS_PARAM_DEF_BY_VALUE;
  pty->flag |= CS_PARAM_FLAG_UNIFORM;

  pty->def.get.tens[0][0] = tens[0][0];
  pty->def.get.tens[0][1] = tens[0][1];
  pty->def.get.tens[0][2] = tens[0][2];
  pty->def.get.tens[1][0] = tens[1][0];
  pty->def.get.tens[1][1] = tens[1][1];
  pty->def.get.tens[1][2] = tens[1][2];
  pty->def.get.tens[2][0] = tens[2][0];
  pty->def.get.tens[2][1] = tens[2][1];
  pty->def.get.tens[2][2] = tens[2][2];

  /* Check the symmetry */
  _check_tensor_symmetry(pty, pty->def.get);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic way to define the value of a cs_property_t structure
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_value(cs_property_t    *pty,
                         const char       *val)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->def_type = CS_PARAM_DEF_BY_VALUE;
  pty->flag |= CS_PARAM_FLAG_UNIFORM;

  switch (pty->type) {

  case CS_PROPERTY_ISO:
    cs_param_set_get(CS_PARAM_VAR_SCAL, (const void *)val, &(pty->def.get));
    break;

  case CS_PROPERTY_ORTHO:
    cs_param_set_get(CS_PARAM_VAR_VECT, (const void *)val, &(pty->def.get));
    break;

  case CS_PROPERTY_ANISO:
    cs_param_set_get(CS_PARAM_VAR_TENS, (const void *)val, &(pty->def.get));

    /* Check the symmetry */
    _check_tensor_symmetry(pty, pty->def.get);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _(" Invalid type of property."));
    break;

  } /* switch on property type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an analytic function
 *
 * \param[in, out]  pty     pointer to a cs_property_t structure
 * \param[in]       func    pointer to a function
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_analytic(cs_property_t        *pty,
                            cs_analytic_func_t   *func)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  pty->def.analytic = func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to a law function
 *
 * \param[in, out]  pty     pointer to a cs_property_t structure
 * \param[in]       func    pointer to a function
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_law(cs_property_t          *pty,
                       cs_onevar_law_func_t   *func)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->def_type = CS_PARAM_DEF_BY_ONEVAR_LAW;
  pty->def.law1_func = func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an array of values
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       support   flag to know where is defined the values
 * \param[in]       array     pointer to an array
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_array(cs_property_t    *pty,
                         cs_flag_t         support,
                         cs_real_t        *array)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->def_type = CS_PARAM_DEF_BY_ARRAY;
  pty->array_flag = support;
  pty->array = array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure in a subdomain by value
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_subdomain_by_value(cs_property_t    *pty,
                                   const char       *ml_name,
                                   const char       *key_val)
{
  int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  _check_def_by_subdomain(pty, ml_name, ml_id);

  int  new_id = pty->n_subdomains;

  _add_def_by_subdomain(pty, new_id, ml_id);

  cs_property_subdef_t  *new_sub = pty->sub_defs + new_id;

  new_sub->def_type = CS_PARAM_DEF_BY_VALUE;

  switch (pty->type) {

  case CS_PROPERTY_ISO:
    cs_param_set_get(CS_PARAM_VAR_SCAL, (const void *)key_val,
                     &(new_sub->def.get));
    break;

  case CS_PROPERTY_ORTHO:
    cs_param_set_get(CS_PARAM_VAR_VECT, (const void *)key_val,
                     &(new_sub->def.get));
    break;

  case CS_PROPERTY_ANISO:
    cs_param_set_get(CS_PARAM_VAR_TENS, (const void *)key_val,
                     &(new_sub->def.get));

    /* Check the symmetry */
    _check_tensor_symmetry(pty, new_sub->def.get);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _(" Invalid type of property."));
    break;

  } /* switch on property type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an isotropic cs_property_t structure in a subdomain by value
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       val       value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_iso_def_subdomain_by_value(cs_property_t    *pty,
                                       const char       *ml_name,
                                       double            val)
{
  int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  _check_def_by_subdomain(pty, ml_name, ml_id);

  if (pty->type != CS_PROPERTY_ISO)
    bft_error(__FILE__, __LINE__, 0,
              " One wants to set an isotropic property but %s property is not"
              " isotropic", pty->name);

  int  new_id = pty->n_subdomains;

  _add_def_by_subdomain(pty, new_id, ml_id);

  cs_property_subdef_t  *new_sub = pty->sub_defs + new_id;

  new_sub->def_type = CS_PARAM_DEF_BY_VALUE;
  new_sub->def.get.val = val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define orthotropic cs_property_t structure in a subdomain by value
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       val       values to set (vector of size 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_ortho_def_subdomain_by_value(cs_property_t    *pty,
                                         const char       *ml_name,
                                         const double      val[])
{
  int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  _check_def_by_subdomain(pty, ml_name, ml_id);

  if (pty->type != CS_PROPERTY_ORTHO)
    bft_error(__FILE__, __LINE__, 0,
              " One wants to set an orthotropic property but %s property is not"
              " orthotropic", pty->name);

  int  new_id = pty->n_subdomains;

  _add_def_by_subdomain(pty, new_id, ml_id);

  cs_property_subdef_t  *new_sub = pty->sub_defs + new_id;

  new_sub->def_type = CS_PARAM_DEF_BY_VALUE;
  new_sub->def.get.vect[0] = val[0];
  new_sub->def.get.vect[1] = val[1];
  new_sub->def.get.vect[2] = val[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an anisotropic cs_property_t structure in a subdomain by value
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       tens      values to set (3x3 tensor)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_aniso_def_subdomain_by_value(cs_property_t    *pty,
                                         const char       *ml_name,
                                         const double      tens[3][3])
{
  int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  _check_def_by_subdomain(pty, ml_name, ml_id);

  if (pty->type != CS_PROPERTY_ANISO)
    bft_error(__FILE__, __LINE__, 0,
              " One wants to set an anisotropic property but %s property is not"
              " anisotropic", pty->name);

  int  new_id = pty->n_subdomains;

  _add_def_by_subdomain(pty, new_id, ml_id);

  cs_property_subdef_t  *new_sub = pty->sub_defs + new_id;

  new_sub->def_type = CS_PARAM_DEF_BY_VALUE;
  new_sub->def.get.tens[0][0] = tens[0][0];
  new_sub->def.get.tens[0][1] = tens[0][1];
  new_sub->def.get.tens[0][2] = tens[0][2];
  new_sub->def.get.tens[1][0] = tens[1][0];
  new_sub->def.get.tens[1][1] = tens[1][1];
  new_sub->def.get.tens[1][2] = tens[1][2];
  new_sub->def.get.tens[2][0] = tens[2][0];
  new_sub->def.get.tens[2][1] = tens[2][1];
  new_sub->def.get.tens[2][2] = tens[2][2];

  /* Check the symmetry */
  _check_tensor_symmetry(pty, new_sub->def.get);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure in a subdomain by law
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       struc     pointer to a structure (may be NULL)
 * \param[in]       func      pointer to a law function defined by subdomain
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_subdomain_by_law(cs_property_t             *pty,
                                 const char                *ml_name,
                                 const void                *struc,
                                 cs_onevar_law_func_t      *func)
{
  int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  _check_def_by_subdomain(pty, ml_name, ml_id);

  int  new_id = pty->n_subdomains;

  _add_def_by_subdomain(pty, new_id, ml_id);

  cs_property_subdef_t  *new_sub = pty->sub_defs + new_id;

  new_sub->def_type = CS_PARAM_DEF_BY_ONEVAR_LAW;
  new_sub->def.law1_func = func;
  new_sub->struc = struc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set "array" members of a cs_property_t structure
 *
 * \param[in, out]  pty          pointer to a cs_property_t structure
 * \param[in]       array_flag   information on the support of the array
 * \param[in]       array        pointer to an array of values
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_array(cs_property_t    *pty,
                      cs_flag_t         array_flag,
                      cs_real_t        *array)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->array_flag = array_flag;
  pty->array = array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set "array" members of a cs_property_t structure
 *
 * \param[in, out]  pty          pointer to a cs_property_t structure
 * \param[in]       structure    structure to associate to this property
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_struct(cs_property_t    *pty,
                       const void       *structure)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->struc = structure;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set optional parameters related to a cs_property_t structure
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       keyname   name of key related to the member of pty to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_option(cs_property_t    *pty,
                       const char       *keyname,
                       const char       *keyval)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  ptykey_t  key = _get_ptykey(keyname);

  if (key == PTYKEY_ERROR) {

    bft_printf("\n\n Current key: %s\n", keyname);
    bft_printf(" Possible keys: ");
    for (int i = 0; i < PTYKEY_ERROR; i++) {
      bft_printf("%s ", _print_ptykey(i));
      if (i > 0 && i % 3 == 0)
        bft_printf("\n\t");
    }
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the property %s.\n"
                " Please read listing for more details and"
                " modify your settings."), pty->name);

  } /* Error message */

  switch(key) {

  case PTYKEY_POST_FREQ:
    pty->post_freq = atoi(keyval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key %s is not implemented yet."), keyname);

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached a property at the cell
 *         center
 *
 * \param[in]      c_id           id of the current cell
 * \param[in]      pty            pointer to a cs_property_t structure
 * \param[in]      do_inversion   true or false
 * \param[in, out] tensor         3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_cell_tensor(cs_lnum_t             c_id,
                            const cs_property_t  *pty,
                            bool                  do_inversion,
                            cs_real_3_t          *tensor)
{
  int  k, l;
  cs_get_t  get;

  if (pty == NULL)
    return;

  /* Initialize extra-diag. values of the tensor */
  for (k = 0; k < 3; k++)
    for (l = k+1; l < 3; l++)
      tensor[k][l] = 0;

  switch (pty->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    _get_tensor_by_value(pty, pty->def.get, tensor);
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    {
      const cs_real_t  *xc = pty->cdoq->cell_centers + 3*c_id;
      const double  t_cur = pty->time_step->t_cur;

      /* Call the analytic function. result is stored in get and then converted
         into a 3x3 tensor */
      pty->def.analytic(t_cur, xc, &get);
      _get_tensor_by_value(pty, get, tensor);

    }
    break;

  case CS_PARAM_DEF_BY_ONEVAR_LAW:
    _get_result_by_onevar_law(c_id, pty, pty->def.law1_func, pty->struc, &get);
    _get_tensor_by_value(pty, get, tensor);
    break;

 case CS_PARAM_DEF_BY_SUBDOMAIN:
   {
     int  def_id = pty->def_ids[c_id];
     cs_property_subdef_t  *sub = pty->sub_defs + def_id;

     switch (sub->def_type) {

     case CS_PARAM_DEF_BY_VALUE:
       _get_tensor_by_value(pty, sub->def.get, tensor);
       break;

     case CS_PARAM_DEF_BY_ONEVAR_LAW:
       _get_result_by_onevar_law(c_id, pty, sub->def.law1_func, sub->struc,
                                 &get);
       _get_tensor_by_value(pty, get, tensor);
       break;

     default:
       bft_error(__FILE__, __LINE__, 0,
                 " Stop computing the cell value related to property %s.\n"
                 " Type of definition for this subdomain is not handled yet.",
                 pty->name);
       break;

     } // Type of definition on this subdomain

   }
   break; // DEF_BY_SUBDOMAIN

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the cell tensor related to property %s.\n"
              " Type of definition not handled yet.", pty->name);
    break;

  } /* type of definition */

  if (do_inversion) {

#if defined(DEBUG) && !defined(NDEBUG) && CS_PROPERTY_DBG > 0
    /* Sanity check */
    for (k = 0; k < 3; k++)
      if (fabs(tensor[k][k]) < cs_get_zero_threshold())
        bft_error(__FILE__, __LINE__, 0,
                  " Potential problem in the inversion of the tensor attached"
                  " to property %s in cell %d.\n"
                  " Tensor[%d][%d] = %5.3e",
                  pty->name, c_id, k, k, tensor[k][k]);
#endif

    if (pty->type == CS_PROPERTY_ISO || pty->type == CS_PROPERTY_ORTHO)
      for (k = 0; k < 3; k++)
        tensor[k][k] /= 1.0;

    else { /* anisotropic */

      cs_real_33_t  invmat;

      _invmat33((const cs_real_3_t (*))tensor, invmat);
      for (k = 0; k < 3; k++)
        for (l = 0; l < 3; l++)
          tensor[k][l] = invmat[k][l];

    }

  } /* Inversion of the tensor */

#if defined(DEBUG) && !defined(NDEBUG) && CS_PROPERTY_DBG > 1
  bft_printf("\n  Tensor property for cell %d\n"
             "   | % 10.6e  % 10.6e  % 10.6e |\n"
             "   | % 10.6e  % 10.6e  % 10.6e |\n"
             "   | % 10.6e  % 10.6e  % 10.6e |\n", c_id,
             tensor[0][0], tensor[0][1], tensor[0][2],
             tensor[1][0], tensor[1][1], tensor[1][2],
             tensor[2][0], tensor[2][1], tensor[2][2]);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a property at the cell center
 *
 * \param[in]   c_id           id of the current cell
 * \param[in]   pty            pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_get_cell_value(cs_lnum_t              c_id,
                           const cs_property_t   *pty)
{
  cs_real_t  result = 0;
  cs_get_t  get;

 if (pty == NULL)
    return result;

 if (pty->type != CS_PROPERTY_ISO)
   bft_error(__FILE__, __LINE__, 0,
             " Invalid type of property for this function.\n"
             " Property %s has to be isotropic.", pty->name);

 switch (pty->def_type) {

 case CS_PARAM_DEF_BY_VALUE:
   result = pty->def.get.val;
   break;

 case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
   {
     const cs_real_t  *xc = pty->cdoq->cell_centers + 3*c_id;
     const double  t_cur = pty->time_step->t_cur;

     /* Call the analytic function. result is stored in get */
     pty->def.analytic(t_cur, xc, &get);
     result = get.val;
   }
   break;

 case CS_PARAM_DEF_BY_ONEVAR_LAW:
   _get_result_by_onevar_law(c_id, pty, pty->def.law1_func, pty->struc, &get);
   result = get.val;
   break;

 case CS_PARAM_DEF_BY_SUBDOMAIN:
   {
     int  def_id = pty->def_ids[c_id];
     cs_property_subdef_t  *sub = pty->sub_defs + def_id;

     switch (sub->def_type) {

     case CS_PARAM_DEF_BY_VALUE:
       result = sub->def.get.val;
       break;

     case CS_PARAM_DEF_BY_ONEVAR_LAW:
       _get_result_by_onevar_law(c_id, pty, sub->def.law1_func, sub->struc,
                                 &get);
       result = get.val;
       break;

     default:
       bft_error(__FILE__, __LINE__, 0,
                 " Stop computing the cell value related to property %s.\n"
                 " Type of definition for this subdomain is not handled yet.",
                 pty->name);
       break;

     } // Type of definition on this subdomain

   }
   break; // DEF_BY_SUBDOMAIN

 default:
   bft_error(__FILE__, __LINE__, 0,
             " Stop computing the cell value related to property %s.\n"
             " Type of definition not handled yet.", pty->name);
   break;

 } /* type of definition */

 return result;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
