/*============================================================================
 * Manage the definition/setting of properties
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

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_defs.h"
#include "cs_log.h"
#include "cs_math.h"
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

#define CS_PROPERTY_DBG  0

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_pty[] =
  " Stop setting an empty cs_property_t structure.\n"
  " Please check your settings.\n";

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;
static const cs_time_step_t  *cs_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the settings are valid
 *
 * \param[in]  pty       pointer to a cs_property_t structure
 * \param[in]  get       accessor to the tensor values
 */
/*----------------------------------------------------------------------------*/

static inline void
_check_tensor_symmetry(const cs_property_t    *pty,
                       cs_get_t                get)
{
  if ((get.tens[0][1] - get.tens[1][0]) > cs_math_zero_threshold ||
      (get.tens[0][2] - get.tens[2][0]) > cs_math_zero_threshold ||
      (get.tens[1][2] - get.tens[2][1]) > cs_math_zero_threshold)
    bft_error(__FILE__, __LINE__, 0,
              _(" The definition of the tensor related to the"
                " property %s is not symmetric.\n"
                " This case is not handled. Please check your settings.\n"),
              pty->name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a new definition to a cs_property_t structure defined by domain
 *         Sanity checks on the settings related to this definition.

 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 *
 * \return a pointer to a new definition to set
 */
/*----------------------------------------------------------------------------*/

static cs_param_def_t *
_init_new_def(cs_property_t     *pty,
              const char        *ml_name)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  int  new_id = pty->n_subdomains;

  if (new_id == pty->n_max_subdomains)
    bft_error(__FILE__, __LINE__, 0,
              _(" Max. number of subdomains has been reached for property %s.\n"
                " Please check your settings."), pty->name);

  int  ml_id = cs_mesh_location_get_id_by_name(ml_name);

  if (ml_id == -1)
    bft_error(__FILE__, __LINE__, 0,
              _(" mesh location %s has not been found.\n"
                " Please check your settings."), ml_name);

  if (cs_mesh_location_get_type(ml_id) != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of mesh location for mesh location  %s.\n"
                " The expected type is CS_MESH_LOCATION_CELLS.\n"), ml_name);

  pty->n_subdomains += 1;

  cs_param_def_t  *new_def = pty->defs + new_id;

  new_def->ml_id = ml_id;

  return new_def;
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
  switch (pty->type) {

  case CS_PROPERTY_ISO:
    tensor[0][0] = tensor[1][1] = tensor[2][2] = get.val;
    break;

  case CS_PROPERTY_ORTHO:
    for (int k = 0; k < 3; k++)
      tensor[k][k] = get.vect[k];
    break;

  case CS_PROPERTY_ANISO:
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        tensor[k][l] = get.tens[k][l];
    break;

  default:
    assert(0);
    break;

  } // Property type
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the value for a cell from an array
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]  cm        pointer to a cs_cell_mesh_t structure
 * \param[in]  desc      information about the array to handle
 * \param[in]  array     array of values (mesh view not cellwise)
 */
/*----------------------------------------------------------------------------*/

static double
_get_cell_val_from_array_cm(const cs_cell_mesh_t    *cm,
                            const cs_desc_t          desc,
                            const cs_real_t          array[])
{
  double  cell_val = 0.;

  if (!(desc.location & CS_FLAG_SCALAR))
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of variable. Stop computing the cell value.");

  /* Test if flag has the pattern of a reference support */
  if ((desc.location & cs_cdo_primal_cell) == cs_cdo_primal_cell)
    cell_val = array[cm->c_id];

  else if ((desc.location & cs_cdo_primal_vtx) == cs_cdo_primal_vtx) {

    /* Sanity checks */
    assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_PVQ));

    /* Reconstruct (or interpolate) value at the current cell center */
    for (short int v = 0; v < cm->n_vc; v++)
      cell_val += cm->wvc[v] * array[cm->v_ids[v]];

  }
  else
    bft_error(__FILE__, __LINE__, 0,
              " Invalid support for evaluating the value of an array");

  return cell_val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the vector at the cell center for cell c_id from an array
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]       cm        pointer to a cs_cell_mesh_t structure
 * \param[in]       desc      information about the array to handle
 * \param[in]       array     values
 * \param[in, out]  vect_val  vector at the cell center
 */
/*----------------------------------------------------------------------------*/

static void
_get_cell_vec_from_array_cm(const cs_cell_mesh_t   *cm,
                            const cs_desc_t         desc,
                            const cs_real_t         array[],
                            cs_real_3_t             vect_val)
{
  /* Test if flag has the pattern of a reference support */
  if ((desc.location & cs_cdo_primal_cell) == cs_cdo_primal_cell)
    for (int k = 0; k < 3; k++)
      vect_val[k] = array[3*cm->c_id+k];

  else if ((desc.location & cs_cdo_dual_face_byc) == cs_cdo_dual_face_byc)

    /* Reconstruct (or interpolate) value at the current cell center */
    cs_reco_dfbyc_in_cell(cm,
                          array + cs_cdo_connect->c2e->idx[cm->c_id],
                          vect_val);

  else
    bft_error(__FILE__, __LINE__, 0,
              " Invalid support for evaluating the vector from an array.");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the value at the cell center for cell c_id from an array
 *
 * \param[in]  c_id      id of the cell to treat
 * \param[in]  desc      information about the array to handle
 * \param[in]  array     values
 */
/*----------------------------------------------------------------------------*/

static double
_get_cell_val_from_array(cs_lnum_t         c_id,
                         const cs_desc_t   desc,
                         const cs_real_t   array[])
{
  double  cell_val = 0.;

  if (!(desc.location & CS_FLAG_SCALAR))
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of variable. Stop computing the cell value.");

  /* Test if flag has the pattern of a reference support */
  if ((desc.location & cs_cdo_primal_cell) == cs_cdo_primal_cell)
    cell_val = array[c_id];

  else if ((desc.location & cs_cdo_primal_vtx) == cs_cdo_primal_vtx)

    /* Reconstruct (or interpolate) value at the current cell center */
    cs_reco_pv_at_cell_center(c_id,
                              cs_cdo_connect->c2v,
                              cs_cdo_quant,
                              array,
                              &cell_val);

  else
    bft_error(__FILE__, __LINE__, 0,
              " Invalid support for evaluating the value of an array");

  return cell_val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the vector at the cell center for cell c_id from an array
 *
 * \param[in]       c_id      if of the cell to treat
 * \param[in]       desc      information about the array to handle
 * \param[in]       array     values
 * \param[in, out]  vect_val  vector at the cell center
 */
/*----------------------------------------------------------------------------*/

static void
_get_cell_vec_from_array(cs_lnum_t         c_id,
                         const cs_desc_t   desc,
                         const cs_real_t   array[],
                         cs_real_3_t       vect_val)
{
  /* Test if flag has the pattern of a reference support */
  if ((desc.location & cs_cdo_primal_cell) == cs_cdo_primal_cell)
    for (int k = 0; k < 3; k++)
      vect_val[k] = array[3*c_id+k];

  else if ((desc.location & cs_cdo_dual_face_byc) == cs_cdo_dual_face_byc)

    /* Reconstruct (or interpolate) value at the current cell center */
    cs_reco_dfbyc_at_cell_center(c_id,
                                 cs_cdo_connect->c2e, cs_cdo_quant, array,
                                 vect_val);

  else
    bft_error(__FILE__, __LINE__, 0,
              " Invalid support for evaluating the vector from an array.");
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                const cs_cdo_connect_t       *connect,
                                const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
  cs_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new property structure
 *
 * \param[in]  name          name of the property
 * \param[in]  key_type      keyname of the type of property
 * \param[in]  n_subdomains  piecewise definition on n_subdomains
 *
 * \return a pointer to a new allocated cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

cs_property_t *
cs_property_create(const char    *name,
                   const char    *key_type,
                   int            n_subdomains)
{
  cs_property_t  *pty = NULL;

  BFT_MALLOC(pty, 1, cs_property_t);

  pty->post_flag = 0;

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(pty->name, len, char);
  strncpy(pty->name, name, len);

  /* Assign a type */
  if (strcmp(key_type, "isotropic") == 0)
    pty->type = CS_PROPERTY_ISO;
  else if (strcmp(key_type, "orthotropic") == 0)
    pty->type = CS_PROPERTY_ORTHO;
  else if (strcmp(key_type, "anisotropic") == 0)
    pty->type = CS_PROPERTY_ANISO;
  else
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key \"%s\" for setting the type of property.\n"
                " Valid keys: \"isotropic\", \"orthotropic\" and"
                " \"anisotropic\".\n"
                " Please modify your settings."), key_type);

  /* Default initialization */
  pty->flag.location = 0;
  pty->flag.state = 0;

  /* Members useful to define the property by subdomain */
  assert(n_subdomains > 0);
  pty->n_max_subdomains = n_subdomains;
  pty->n_subdomains = 0; // No definition set by default

  BFT_MALLOC(pty->defs, n_subdomains, cs_param_def_t);

  pty->def_ids = NULL;

  /* Specific members for more complex definitions */
  pty->desc1.location = pty->desc1.state = 0;
  pty->array1 = NULL;
  pty->desc2.location = pty->desc2.state = 0;
  pty->array2 = NULL;

  return pty;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last stage of the definition of a property based on several
 *         subdomains
 *
 * \param[in, out] pty       pointer to cs_property_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_property_last_definition_stage(cs_property_t  *pty)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  if (pty->n_subdomains > 1) { /* Initialization of def_ids */

    const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;

    BFT_MALLOC(pty->def_ids, n_cells, short int);

#   pragma omp parallel for if (n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells; i++)
      pty->def_ids[i] = -1; // Unset by default

    for (int sub_id = 0; sub_id < pty->n_subdomains; sub_id++) {

      const cs_param_def_t  *def = pty->defs + sub_id;
      const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(def->ml_id);
      const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(def->ml_id);

      if (elt_ids == NULL) {

        assert(n_elts[0] == n_cells);
#       pragma omp parallel for if (n_cells > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_cells; i++)
          pty->def_ids[i] = sub_id;

      }
      else {

#       pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts[0]; i++)
          pty->def_ids[elt_ids[i]] = sub_id;

      }

    } // Loop on subdomains

  } // n_max_subdomains > 1

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
  BFT_FREE(pty->def_ids);
  BFT_FREE(pty->defs);

  if (pty->desc1.state & CS_FLAG_STATE_OWNER)
    if (pty->array1 != NULL)
      BFT_FREE(pty->array1);

  if (pty->desc2.state & CS_FLAG_STATE_OWNER)
    if (pty->array2 != NULL)
      BFT_FREE(pty->array2);

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

  if (pty->flag.state & CS_FLAG_STATE_UNIFORM)
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
 * \brief  Set optional parameters related to a cs_property_t structure
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       key       key related to the member of pty to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_option(cs_property_t       *pty,
                       cs_property_key_t    key,
                       const char          *keyval)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {
  case CS_PTYKEY_POST:
    if (strcmp(val, "fourier") == 0)
      pty->post_flag |= CS_PROPERTY_POST_FOURIER;
    else
      bft_error(__FILE__, __LINE__, 0,
                N_(" Invalid value %s for setting key CS_PTYKEY_POST\n"
                   " Valid choices are \"fourier\".\n"
                   " Please modify your setting.\n"), val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key not implemented for setting an advection field."));

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure by value for entities attached to
 *         the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_value(cs_property_t    *pty,
                         const char       *ml_name,
                         const char       *key_val)
{
  cs_param_def_t  *new_def = _init_new_def(pty, ml_name);

  new_def->def_type = CS_PARAM_DEF_BY_VALUE;
  if (pty->n_max_subdomains == 1)
    pty->flag.state |= CS_FLAG_STATE_UNIFORM;

  switch (pty->type) {

  case CS_PROPERTY_ISO:
    cs_param_set_get(CS_PARAM_VAR_SCAL, (const void *)key_val,
                     &(new_def->def.get));
    break;

  case CS_PROPERTY_ORTHO:
    cs_param_set_get(CS_PARAM_VAR_VECT, (const void *)key_val,
                     &(new_def->def.get));
    break;

  case CS_PROPERTY_ANISO:
    cs_param_set_get(CS_PARAM_VAR_TENS, (const void *)key_val,
                     &(new_def->def.get));

    /* Check the symmetry */
    _check_tensor_symmetry(pty, new_def->def.get);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _(" Invalid type of property."));
    break;

  } /* switch on property type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an isotropic cs_property_t structure by value for entities
 *         attached to the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       val       value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_property_iso_def_by_value(cs_property_t    *pty,
                             const char       *ml_name,
                             double            val)
{
  cs_param_def_t  *new_def = _init_new_def(pty, ml_name);

  if (pty->type != CS_PROPERTY_ISO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not isotropic.\n"
              " Please check your settings.", pty->name);

  new_def->def_type = CS_PARAM_DEF_BY_VALUE;
  if (pty->n_max_subdomains == 1)
    pty->flag.state |= CS_FLAG_STATE_UNIFORM;
  new_def->def.get.val = val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define orthotropic cs_property_t structure by value for entities
 *         attached to the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       val       values to set (vector of size 3)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_ortho_def_by_value(cs_property_t    *pty,
                               const char       *ml_name,
                               const double      val[])
{
  cs_param_def_t  *new_def = _init_new_def(pty, ml_name);

  if (pty->type != CS_PROPERTY_ORTHO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not orthotropic.\n"
              " Please check your settings.", pty->name);

  new_def->def_type = CS_PARAM_DEF_BY_VALUE;
  if (pty->n_max_subdomains == 1)
    pty->flag.state |= CS_FLAG_STATE_UNIFORM;

  new_def->def.get.vect[0] = val[0];
  new_def->def.get.vect[1] = val[1];
  new_def->def.get.vect[2] = val[2];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define an anisotropic cs_property_t structure by value for entities
 *         attached to the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       tens      values to set (3x3 tensor)
 */
/*----------------------------------------------------------------------------*/

void
cs_property_aniso_def_by_value(cs_property_t    *pty,
                               const char       *ml_name,
                               const double      tens[3][3])
{
  cs_param_def_t  *new_def = _init_new_def(pty, ml_name);

  if (pty->type != CS_PROPERTY_ANISO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid setting: property %s is not anisotropic.\n"
              " Please check your settings.", pty->name);

  new_def->def_type = CS_PARAM_DEF_BY_VALUE;
  if (pty->n_max_subdomains == 1)
    pty->flag.state |= CS_FLAG_STATE_UNIFORM;

  new_def->def.get.tens[0][0] = tens[0][0];
  new_def->def.get.tens[0][1] = tens[0][1];
  new_def->def.get.tens[0][2] = tens[0][2];
  new_def->def.get.tens[1][0] = tens[1][0];
  new_def->def.get.tens[1][1] = tens[1][1];
  new_def->def.get.tens[1][2] = tens[1][2];
  new_def->def.get.tens[2][0] = tens[2][0];
  new_def->def.get.tens[2][1] = tens[2][1];
  new_def->def.get.tens[2][2] = tens[2][2];

  /* Check the symmetry */
  _check_tensor_symmetry(pty, new_def->def.get);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an analytic function in
 *         a subdomain attached to the mesh location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       func      pointer to a cs_analytic_func_t function
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_analytic(cs_property_t        *pty,
                            const char           *ml_name,
                            cs_analytic_func_t   *func)
{
  cs_param_def_t  *new_def = _init_new_def(pty, ml_name);

  new_def->def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  new_def->def.analytic = func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to law depending on one
 *         scalar variable in a subdomain attached to the mesh location named
 *         ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       context   pointer to a structure (may be NULL)
 * \param[in]       func      pointer to a law function defined by subdomain
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_onevar_law(cs_property_t             *pty,
                              const char                *ml_name,
                              const void                *context,
                              cs_onevar_law_func_t      *func)
{
  cs_param_def_t  *new_def = _init_new_def(pty, ml_name);

  new_def->def_type = CS_PARAM_DEF_BY_ONEVAR_LAW;
  new_def->def.law1_func = func;
  new_def->context = context;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to a law depending on
 *         two scalars variables in a subdomain attached to the mesh location
 *         named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       context   pointer to a structure (may be NULL)
 * \param[in]       func      pointer to a function
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_twovar_law(cs_property_t          *pty,
                              const char             *ml_name,
                              const void             *context,
                              cs_twovar_law_func_t   *func)
{
  cs_param_def_t  *new_def = _init_new_def(pty, ml_name);

  new_def->def_type = CS_PARAM_DEF_BY_TWOVAR_LAW;
  new_def->def.law2_func = func;
  new_def->context = context;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to an array of values
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       desc      information about this array
 * \param[in]       array     pointer to an array
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_array(cs_property_t    *pty,
                         cs_desc_t         desc,
                         cs_real_t        *array)
{
  cs_param_def_t  *new_def = _init_new_def(pty, "cells");

  if (pty->n_max_subdomains != 1)
    bft_error(__FILE__, __LINE__, 0,
              " When a definition by array is requested, the max. number"
              " of subdomains to consider should be equal to 1.\n"
              " Current value is %d for property %s.\n"
              " Please modify your settings.",
              pty->n_max_subdomains, pty->name);

  new_def->def_type = CS_PARAM_DEF_BY_ARRAY;
  pty->desc1.location = desc.location;
  pty->desc1.state = desc.state;
  pty->array1 = array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the "array" member of a cs_property_t structure
 *
 * \param[in, out]  pty          pointer to a cs_property_t structure
 * \param[in]       desc         information about this array
 * \param[in]       array        pointer to an array of values
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_array(cs_property_t    *pty,
                      cs_desc_t         desc,
                      cs_real_t        *array)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->desc1.location = desc.location;
  pty->desc1.state = desc.state;
  pty->array1 = array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the second "array" member of a cs_property_t structure
 *
 * \param[in, out]  pty        pointer to a cs_property_t structure
 * \param[in]       desc       information about this array
 * \param[in]       array      pointer to an array of values
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_second_array(cs_property_t    *pty,
                             cs_desc_t         desc,
                             cs_real_t        *array)
{
  if (pty == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_pty));

  pty->desc2.location = desc.location;
  pty->desc2.state = desc.state;
  pty->array2 = array;
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
  if (pty == NULL)
    return;

  /* Initialize extra-diag. values of the tensor */
  tensor[0][1] = tensor[1][0] = tensor[2][0] = 0;
  tensor[0][2] = tensor[1][2] = tensor[2][1] = 0;

  int  def_id = -1;
  if (pty->n_max_subdomains == 1)
    def_id = 0;
  else
    def_id = pty->def_ids[c_id];

  cs_param_def_t  *sub = pty->defs + def_id;

  switch (sub->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    _get_tensor_by_value(pty, sub->def.get, tensor);
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    {
      const cs_real_t  *xc = cs_cdo_quant->cell_centers + 3*c_id;
      const double  t_cur = cs_time_step->t_cur;

      /* Call the analytic function. result is stored in get and then converted
         into a 3x3 tensor */
      switch (pty->type) {

      case CS_PROPERTY_ISO:
        {
          double  eval;

          sub->def.analytic(t_cur, 1, xc, &eval);
          tensor[0][0] = tensor[1][1] = tensor[2][2] = eval;
        }
        break;

      case CS_PROPERTY_ORTHO:
        {
          double  eval[3];

          sub->def.analytic(t_cur, 1, xc, eval);
          for (int k = 0; k < 3; k++)
            tensor[k][k] = eval[k];
        }
        break;

      case CS_PROPERTY_ANISO:
        sub->def.analytic(t_cur, 1, xc, (cs_real_t *)tensor);
        break;

      default:
        assert(0);
        break;
      }

    }
    break;

  case CS_PARAM_DEF_BY_ONEVAR_LAW:
    {
      assert(pty->array1 != NULL); /* Sanity check */

      cs_real_t  val_xc = _get_cell_val_from_array(c_id,
                                                   pty->desc1, pty->array1);

      sub->def.law1_func(1, NULL, &val_xc, sub->context, (cs_real_t *)tensor);

    }
    break;

  case CS_PARAM_DEF_BY_TWOVAR_LAW:
    {
      assert(pty->array1 != NULL && pty->array2 != NULL); /* Sanity check */

      cs_real_t  val1 = _get_cell_val_from_array(c_id, pty->desc1, pty->array1);

      if ((pty->desc2.state & CS_FLAG_STATE_POTENTIAL) &&
          (pty->desc2.location & CS_FLAG_SCALAR)) {

        cs_real_t  val2 = _get_cell_val_from_array(c_id,
                                                   pty->desc2, pty->array2);

        /* Compute the value of the law for this cell */
        sub->def.law2_func(1, NULL, &val1, &val2, sub->context,
                           (cs_real_t *)tensor);

      }
      else if ((pty->desc2.state & CS_FLAG_STATE_FLUX) ||
               (pty->desc2.location & CS_FLAG_VECTOR)) {

        cs_real_t  vect_val[3] = {0, 0, 0};
        _get_cell_vec_from_array(c_id, pty->desc2, pty->array2, vect_val);

        /* Retrieve the result of the law function for this cell */
        sub->def.law2_func(1, NULL, &val1, vect_val, sub->context,
                           (cs_real_t *)tensor);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid case. Stop evaluating the property %s\n",
                  pty->name);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the cell tensor related to property %s.\n"
              " Type of definition not handled yet.", pty->name);
    break;

  } /* type of definition */

  if (do_inversion) {

#if defined(DEBUG) && !defined(NDEBUG) && CS_PROPERTY_DBG > 0 /* Sanity check */
    for (int k = 0; k < 3; k++)
      if (fabs(tensor[k][k]) < cs_math_zero_threshold)
        bft_error(__FILE__, __LINE__, 0,
                  " Potential problem in the inversion of the tensor attached"
                  " to property %s in cell %d.\n"
                  " Tensor[%d][%d] = %5.3e",
                  pty->name, c_id, k, k, tensor[k][k]);
#endif

    if (pty->type == CS_PROPERTY_ISO || pty->type == CS_PROPERTY_ORTHO)
      for (int k = 0; k < 3; k++)
        tensor[k][k] /= 1.0;

    else { /* anisotropic */

      cs_real_33_t  invmat;

      cs_math_33_inv((const cs_real_3_t (*))tensor, invmat);
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          tensor[k][l] = invmat[k][l];

    }

  } /* Inversion of the tensor */

#if defined(DEBUG) && !defined(NDEBUG) && CS_PROPERTY_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT,
                "\n  Tensor property for cell %d\n"
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

  if (pty == NULL)
    return result;

  if (pty->type != CS_PROPERTY_ISO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of property for this function.\n"
              " Property %s has to be isotropic.", pty->name);

  int  def_id = -1;
  if (pty->n_max_subdomains == 1)
    def_id = 0;
  else
    def_id = pty->def_ids[c_id];

  cs_param_def_t  *sub = pty->defs + def_id;

  switch (sub->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    result = sub->def.get.val;
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    sub->def.analytic(cs_time_step->t_cur,
                      1,
                      cs_cdo_quant->cell_centers + 3*c_id,
                      &result);
    break;

  case CS_PARAM_DEF_BY_ONEVAR_LAW:
    {
      assert(pty->array1 != NULL); /* Sanity check */

      cs_real_t  val_xc = _get_cell_val_from_array(c_id,
                                                   pty->desc1, pty->array1);

      sub->def.law1_func(1, NULL, &val_xc, sub->context, &result);
    }
    break;

  case CS_PARAM_DEF_BY_TWOVAR_LAW:
    {
      assert(pty->array1 != NULL && pty->array2 != NULL); /* Sanity check */

      cs_real_t  val1 = _get_cell_val_from_array(c_id, pty->desc1, pty->array1);

      if ((pty->desc2.state & CS_FLAG_STATE_POTENTIAL) &&
          (pty->desc2.location & CS_FLAG_SCALAR)) {

        cs_real_t  val2 = _get_cell_val_from_array(c_id,
                                                   pty->desc2, pty->array2);

        /* Compute the value of the law for this cell */
        sub->def.law2_func(1, NULL, &val1, &val2, sub->context, &result);

      }
      else if ((pty->desc2.state & CS_FLAG_STATE_FLUX) ||
               (pty->desc2.location & CS_FLAG_VECTOR)) {

        cs_real_t  vect_val[3] = {0, 0, 0};
        _get_cell_vec_from_array(c_id, pty->desc2, pty->array2, vect_val);

        /* Retrieve the result of the law function for this cell */
        sub->def.law2_func(1, NULL, &val1, vect_val, sub->context, &result);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid case. Stop evaluating the property %s\n",
                  pty->name);

    }
    break;

  case CS_PARAM_DEF_BY_ARRAY:
    result = _get_cell_val_from_array(c_id, pty->desc1, pty->array1);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the cell value related to property %s.\n"
              " Type of definition not handled yet.", pty->name);
    break;

  } /* type of definition */

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the tensor attached a property at the cell
 *         center
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      pty            pointer to a cs_property_t structure
 * \param[in]      do_inversion   true or false
 * \param[in, out] tensor         3x3 matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_property_tensor_in_cell(const cs_cell_mesh_t   *cm,
                           const cs_property_t    *pty,
                           bool                    do_inversion,
                           cs_real_3_t            *tensor)
{
  if (pty == NULL)
    return;

  /* Initialize extra-diag. values of the tensor */
  tensor[0][1] = tensor[1][0] = tensor[2][0] = 0;
  tensor[0][2] = tensor[1][2] = tensor[2][1] = 0;

  int  def_id = -1;
  if (pty->n_max_subdomains == 1)
    def_id = 0;
  else
    def_id = pty->def_ids[cm->c_id];

  cs_param_def_t  *sub = pty->defs + def_id;
  switch (sub->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    _get_tensor_by_value(pty, sub->def.get, tensor);
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    {
      const double  t_cur = cs_time_step->t_cur;

      /* Call the analytic function. result is stored in get and then converted
         into a 3x3 tensor */
      switch (pty->type) {

      case CS_PROPERTY_ISO:
        {
          double  eval;

          sub->def.analytic(t_cur, 1, cm->xc, &eval);
          tensor[0][0] = tensor[1][1] = tensor[2][2] = eval;
        }
        break;

      case CS_PROPERTY_ORTHO:
        {
          double  eval[3];

          sub->def.analytic(t_cur, 1, cm->xc, eval);
          for (int k = 0; k < 3; k++)
            tensor[k][k] = eval[k];
        }
        break;

      case CS_PROPERTY_ANISO:
        sub->def.analytic(t_cur, 1, cm->xc, (cs_real_t *)tensor);
        break;

      default:
        assert(0);
        break;
      }

    }
    break;

  case CS_PARAM_DEF_BY_ONEVAR_LAW:
    {
      assert(pty->array1 != NULL); /* Sanity check */

      cs_real_t  val_xc =
        _get_cell_val_from_array_cm(cm, pty->desc1, pty->array1);

      sub->def.law1_func(1, NULL, &val_xc, sub->context, (cs_real_t *)tensor);

      if (pty->type == CS_PROPERTY_ISO)
        tensor[2][2] = tensor[1][1] = tensor[0][0];
    }
    break;

  case CS_PARAM_DEF_BY_TWOVAR_LAW:
    {
      assert(pty->array1 != NULL && pty->array2 != NULL); /* Sanity check */

      cs_real_t  val1 =
        _get_cell_val_from_array_cm(cm, pty->desc1, pty->array1);

      if ((pty->desc2.state & CS_FLAG_STATE_POTENTIAL) &&
          (pty->desc2.location & CS_FLAG_SCALAR)) {

        cs_real_t  val2 =
          _get_cell_val_from_array_cm(cm, pty->desc2, pty->array2);

        /* Compute the value of the law for this cell */
        sub->def.law2_func(1, NULL, &val1, &val2, sub->context,
                           (cs_real_t *)tensor);

      }
      else if ((pty->desc2.state & CS_FLAG_STATE_FLUX) ||
               (pty->desc2.location & CS_FLAG_VECTOR)) {

        cs_real_t  vect_val[3] = {0, 0, 0};
        _get_cell_vec_from_array_cm(cm, pty->desc2, pty->array2, vect_val);

        /* Retrieve the result of the law function for this cell */
        sub->def.law2_func(1, NULL, &val1, vect_val, sub->context,
                           (cs_real_t *)tensor);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid case. Stop evaluating the property %s\n",
                  pty->name);

      if (pty->type == CS_PROPERTY_ISO)
        tensor[2][2] = tensor[1][1] = tensor[0][0];

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the cell tensor related to property %s.\n"
              " Type of definition not handled yet.", pty->name);
    break;

  } /* type of definition */

  if (do_inversion) {

#if defined(DEBUG) && !defined(NDEBUG) && CS_PROPERTY_DBG > 0 /* Sanity check */
    for (int k = 0; k < 3; k++)
      if (fabs(tensor[k][k]) < cs_math_zero_threshold)
        bft_error(__FILE__, __LINE__, 0,
                  " Potential problem in the inversion of the tensor attached"
                  " to property %s in cell %d.\n"
                  " Tensor[%d][%d] = %5.3e",
                  pty->name, cm->c_id, k, k, tensor[k][k]);
#endif

    if (pty->type == CS_PROPERTY_ISO || pty->type == CS_PROPERTY_ORTHO)
      for (int k = 0; k < 3; k++)
        tensor[k][k] /= 1.0;

    else { /* anisotropic */

      cs_real_33_t  invmat;

      cs_math_33_inv((const cs_real_3_t (*))tensor, invmat);
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          tensor[k][l] = invmat[k][l];

    }

  } /* Inversion of the tensor */

#if defined(DEBUG) && !defined(NDEBUG) && CS_PROPERTY_DBG > 1
  cs_log_printf(CS_LOG_DEFAULT,
                "\n  Tensor property for cell %d\n"
                "   | % 10.6e  % 10.6e  % 10.6e |\n"
                "   | % 10.6e  % 10.6e  % 10.6e |\n"
                "   | % 10.6e  % 10.6e  % 10.6e |\n", cm->c_id,
                tensor[0][0], tensor[0][1], tensor[0][2],
                tensor[1][0], tensor[1][1], tensor[1][2],
                tensor[2][0], tensor[2][1], tensor[2][2]);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of a property at the cell center
 *         Version using a cs_cell_mesh_t structure
 *
 * \param[in]   cm        pointer to a cs_cell_mesh_t structure
 * \param[in]   pty       pointer to a cs_property_t structure
 *
 * \return the value of the property for the given cell
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_property_value_in_cell(const cs_cell_mesh_t   *cm,
                          const cs_property_t    *pty)
{
  cs_real_t  result = 0;

  if (pty == NULL)
    return result;

  if (pty->type != CS_PROPERTY_ISO)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of property for this function.\n"
              " Property %s has to be isotropic.", pty->name);

  int  def_id = -1;
  if (pty->n_max_subdomains == 1)
    def_id = 0;
  else
    def_id = pty->def_ids[cm->c_id];

  cs_param_def_t  *sub = pty->defs + def_id;

  switch (sub->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    result = sub->def.get.val;
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    sub->def.analytic(cs_time_step->t_cur, 1, cm->xc, &result);
    break;

  case CS_PARAM_DEF_BY_ONEVAR_LAW:
    {
      /* Sanity checks */
      assert((pty->array1 != NULL) && (pty->type == CS_PROPERTY_ISO));

      cs_real_t  val_xc =
        _get_cell_val_from_array_cm(cm, pty->desc1, pty->array1);

      sub->def.law1_func(1, NULL, &val_xc, sub->context, &result);
    }
    break;

  case CS_PARAM_DEF_BY_TWOVAR_LAW:
    {
      /* Sanity checks */
      assert(pty->array1 != NULL && pty->array2 != NULL);
      assert(pty->type == CS_PROPERTY_ISO);

      cs_real_t  val1 = _get_cell_val_from_array_cm(cm,
                                                    pty->desc1, pty->array1);

      if ((pty->desc2.state & CS_FLAG_STATE_POTENTIAL) &&
          (pty->desc2.location & CS_FLAG_SCALAR)) {

        cs_real_t  val2 = _get_cell_val_from_array_cm(cm,
                                                      pty->desc2, pty->array2);

        /* Compute the value of the law for this cell */
        sub->def.law2_func(1, NULL, &val1, &val2, sub->context, &result);

      }
      else if ((pty->desc2.state & CS_FLAG_STATE_FLUX) ||
               (pty->desc2.location & CS_FLAG_VECTOR)) {

        cs_real_t  vect_val[3] = {0, 0, 0};
        _get_cell_vec_from_array_cm(cm, pty->desc2, pty->array2, vect_val);

        /* Retrieve the result of the law function for this cell */
        sub->def.law2_func(1, NULL, &val1, vect_val, sub->context, &result);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid case. Stop evaluating the property %s\n",
                  pty->name);
    }
    break;

  case CS_PARAM_DEF_BY_ARRAY:
    result = _get_cell_val_from_array_cm(cm, pty->desc1, pty->array1);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the cell value related to property %s.\n"
              " Type of definition not handled yet.", pty->name);
    break;

  } /* type of definition */

  return result;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Fourier number in each cell
 *
 * \param[in]      pty        pointer to the diffusive property struct.
 * \param[in]      dt         value of the current time step
 * \param[in, out] fourier    pointer to an array storing Fourier numbers
 */
/*----------------------------------------------------------------------------*/

void
cs_property_get_fourier(const cs_property_t     *pty,
                        double                   dt,
                        cs_real_t                fourier[])
{
  assert(fourier != NULL); // Sanity check
  assert(dt > 0.);

  const bool  pty_uniform = cs_property_is_uniform(pty);
  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;

  if (pty->type == CS_PROPERTY_ISO) {

    cs_real_t  ptyval = 0.;
    if (pty_uniform)
      ptyval = cs_property_get_cell_value(0, pty);

    for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

      if (!pty_uniform)
        ptyval = cs_property_get_cell_value(c_id, pty);

      const cs_real_t  hc2 = pow(cdoq->cell_vol[c_id], 2*cs_math_onethird);

      fourier[c_id] = dt * ptyval / hc2;

    } // Loop on cells

  }
  else { // Property is orthotropic or anisotropic

    cs_real_t  eig_max, eig_ratio;
    cs_real_t  ptymat[3][3];

    /* Get the value of the material property at the first cell center */
    if (pty_uniform) {
      cs_property_get_cell_tensor(0, pty, false, ptymat);
      cs_math_33_eigen((const cs_real_t (*)[3])ptymat, &eig_ratio, &eig_max);
    }

    for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

      /* Get the value of the material property at the cell center */
      if (!pty_uniform) {
        cs_property_get_cell_tensor(c_id, pty, false, ptymat);
        cs_math_33_eigen((const cs_real_t (*)[3])ptymat, &eig_ratio, &eig_max);
      }

      const cs_real_t  hc2 = pow(cdoq->cell_vol[c_id], 2*cs_math_onethird);

      fourier[c_id] = dt * eig_max / hc2;

    } // Loop on cells

  } // Type of property

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

  if (pty->flag.state & CS_FLAG_STATE_UNIFORM)  is_uniform = true;
  if (pty->flag.state & CS_FLAG_STATE_UNSTEADY) is_steady = false;

  cs_log_printf(CS_LOG_SETUP, "  %s >> uniform [%s], steady [%s], ",
                pty->name, cs_base_strtf(is_uniform), cs_base_strtf(is_steady));

  switch(pty->type) {
  case CS_PROPERTY_ISO:
    cs_log_printf(CS_LOG_SETUP, "type: isotropic\n");
    break;
  case CS_PROPERTY_ORTHO:
    cs_log_printf(CS_LOG_SETUP, "type: orthotropic\n");
    break;
  case CS_PROPERTY_ANISO:
    cs_log_printf(CS_LOG_SETUP, "type: anisotropic\n");
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of property."));
    break;
  }

  /* Definition */
  cs_log_printf(CS_LOG_SETUP,
                "  %s >> n_subdomains    %d\n", pty->name, pty->n_subdomains);

  for (int i = 0; i < pty->n_subdomains; i++) {

    cs_param_def_t  *pdef  = pty->defs + i;

    cs_log_printf(CS_LOG_SETUP,
                  "  %s >> location  %s,", pty->name,
                  cs_mesh_location_get_name(pdef->ml_id));

    switch (pdef->def_type) {

    case CS_PARAM_DEF_BY_VALUE:
      {
        const cs_get_t  mat = pdef->def.get;

        switch(pty->type) {

        case CS_PROPERTY_ISO:
          cs_log_printf(CS_LOG_SETUP,
                        " definition by value: % 5.3e\n", mat.val);
          break;
        case CS_PROPERTY_ORTHO:
          cs_log_printf(CS_LOG_SETUP,
                        " definition by value: (% 5.3e, % 5.3e, % 5.3e)\n",
                        mat.vect[0], mat.vect[1], mat.vect[2]);
          break;
        case CS_PROPERTY_ANISO:
          cs_log_printf(CS_LOG_SETUP,
                        "\n                       |% 5.3e, % 5.3e, % 5.3e|\n"
                        "  definition by value: |% 5.3e, % 5.3e, % 5.3e|\n"
                        "                       |% 5.3e, % 5.3e, % 5.3e|\n",
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
      cs_log_printf(CS_LOG_SETUP, "  definition by an analytical function\n");
      break;

    case CS_PARAM_DEF_BY_ONEVAR_LAW:
      cs_log_printf(CS_LOG_SETUP,
                    "  definition by a law based on one variable\n");
      break;

    case CS_PARAM_DEF_BY_TWOVAR_LAW:
      cs_log_printf(CS_LOG_SETUP,
                    "  definition by law based on two variables\n");
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of definition for a property."));
      break;

    } /* switch on def_type */

  } // Loop on subdomains

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
