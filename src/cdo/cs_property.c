/*============================================================================
 * Manage the definition/setting of properties
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
#include <assert.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_defs.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_reco.h"
#include "cs_timer_stats.h"

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

/* Set of parameters attached to a property */
struct _cs_property_t {

  char  *restrict name;

  /* Short descriptor to know where is defined the property and what kind of
     property one considers (mask of bits) */

  cs_desc_t   flag;

  /* The number of values to set depends on the type of property
     - isotropic   = 1 => CS_PARAM_VAR_SCAL
     - orthotropic = 3 => CS_PARAM_VAR_VECT
     - anisotropic = 9 => CS_PARAM_VAR_TENS
  */

  cs_property_type_t   type;  // isotropic, anistotropic...

  /* Members to define the value of the property by subdomains (up to now,
     only subdomains built from an union of cells are considered) */

  int               n_max_subdomains; // requested number of subdomains
  int               n_subdomains;     // current number of subddomains defined
  cs_param_def_t   *defs;             // list of definitions for each subdomain
  short int        *def_ids;          /* id in of the definition related to each
                                         cell.
                                         NULL is only one definition is set */

  /* Useful buffers to deal with more complex definitions */

  cs_real_t   *array1;   // if the property hinges on an array
  cs_desc_t    desc1;    // short description of the related array
  cs_real_t   *array2;   // if the property hinges on a second array
  cs_desc_t    desc2;    // short description of the related array

};

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

/* Id related to cs_timer_stats structure used for monitoring */
static int  property_ts_id = -1;

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

  if (pty->n_max_subdomains > 1) { /* Assign new id to the selected cells */

    const cs_lnum_t  n_cells = cs_cdo_quant->n_cells;
    const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
    const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

    if (elt_ids == NULL) {

      assert(n_elts[0] == n_cells); // Sanity check. Only this case is handled
# pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_elts[0]; i++)
        pty->def_ids[i] = new_id;

    }
    else {

# pragma omp parallel for if (n_elts[0] > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < n_elts[0]; i++)
        pty->def_ids[elt_ids[i]] = new_id;

    }

  } // n_max_subdomains > 1

  pty->n_subdomains += 1;

  cs_param_def_t  *new_def = pty->defs + new_id;

  /* Copy mesh location name */
  int  len = strlen(ml_name) + 1;
  BFT_MALLOC(new_def->ml_name, len, char);
  strncpy(new_def->ml_name, ml_name, len);

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
 * \brief  Retrieve the value at the cell center of cell c_id from an array
 *
 * \param[in]  c_id      if of the cell to treat
 * \param[in]  desc      information about the array to handle
 * \param[in]  array     values
 */
/*----------------------------------------------------------------------------*/

static double
_get_cell_value_from_array(cs_lnum_t         c_id,
                           const cs_desc_t   desc,
                           const cs_real_t   array[])
{
  double  cell_val = 0.;

  if ((desc.location & CS_FLAG_SCAL) != CS_FLAG_SCAL)
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
 * \brief  Retrieve the vector at the cell center of cell c_id from an array
 *
 * \param[in]       c_id      if of the cell to treat
 * \param[in]       desc      information about the array to handle
 * \param[in]       array     values
 * \param[in, out]  vect_val  vector at the cell center
 */
/*----------------------------------------------------------------------------*/

static void
_get_cell_vector_from_array(cs_lnum_t         c_id,
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

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value using a law with one scalar variable
 *
 * \param[in]       c_id     cell id
 * \param[in]       pty      pointer to a cs_property_t structure
 * \param[in]       law      function pointer to the law
 * \param[in]       context  pointer to a structure
 * \param[in, out]  get      pointer to a union used to retrieve the result
 */
/*----------------------------------------------------------------------------*/

static void
_get_result_by_onevar_law(cs_lnum_t                 c_id,
                          const cs_property_t      *pty,
                          cs_onevar_law_func_t     *law,
                          const void               *context,
                          cs_get_t                 *get)
{
  assert(pty->array1 != NULL); /* Sanity check */

  cs_real_t  val_xc = _get_cell_value_from_array(c_id, pty->desc1, pty->array1);

  law(val_xc, context, get);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value using a law with two scalar variables
 *
 * \param[in]       c_id     cell id
 * \param[in]       pty      pointer to a cs_property_t structure
 * \param[in]       law      function pointer to the law
 * \param[in]       context  pointer to a structure
 * \param[in, out]  get      pointer to a union used to retrieve the result
 */
/*----------------------------------------------------------------------------*/

static void
_get_result_by_twovar_law(cs_lnum_t                 c_id,
                          const cs_property_t      *pty,
                          cs_twovar_law_func_t     *law,
                          const void               *context,
                          cs_get_t                 *get)
{
  assert(pty->array1 != NULL && pty->array2 != NULL); /* Sanity check */

  cs_real_t  val1 = _get_cell_value_from_array(c_id, pty->desc1, pty->array1);
  cs_real_t  val2 = _get_cell_value_from_array(c_id, pty->desc2, pty->array2);

  /* Compute the value of the law for this cell */
  law(val1, val2, context, get);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value using a law with two arguments: one scalar and
 *         one vector-valued
 *
 * \param[in]       c_id    cell id
 * \param[in]       pty     pointer to a cs_property_t structure
 * \param[in]       law     function pointer to the law
 * \param[in]       context   pointer to a structure
 * \param[in, out]  get     pointer to a union used to retrieve the result
 */
/*----------------------------------------------------------------------------*/

static void
_get_result_by_scavec_law(cs_lnum_t                 c_id,
                          const cs_property_t      *pty,
                          cs_scavec_law_func_t     *law,
                          const void               *context,
                          cs_get_t                 *get)
{
  /* Sanity checks */
  assert(pty->array1 != NULL && pty->array2 != NULL);

  cs_real_t  scal_val = _get_cell_value_from_array(c_id,
                                                   pty->desc1, pty->array1);
  cs_real_t  vect_val[3] = {0, 0, 0};
  _get_cell_vector_from_array(c_id, pty->desc2, pty->array2, vect_val);

  /* Retrieve the result of the law function */
  law(scal_val, vect_val, context, get);
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
 * \brief   Initialize cs_timer_stats_t structure for monitoring purpose
 *
 * \param[in]  level      level of details requested
 */
/*----------------------------------------------------------------------------*/

void
cs_property_set_timer_stats(int   level)
{
  if (level < 1)
    return;

  /* Timer statistics */
  property_ts_id = cs_timer_stats_create("operations", "property", "property");

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
  if (n_subdomains > 1) { /* Initialization of def_ids */

    BFT_MALLOC(pty->def_ids, cs_cdo_quant->n_cells, short int);

# pragma omp parallel for if (cs_cdo_quant->n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < cs_cdo_quant->n_cells; i++)
      pty->def_ids[i] = -1; // Unset by default

  }

  /* Specific members for more complex definitions */
  pty->desc1.location = pty->desc1.state = 0;
  pty->array1 = NULL;
  pty->desc2.location = pty->desc2.state = 0;
  pty->array2 = NULL;

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
  BFT_FREE(pty->def_ids);
  for (int i = 0; i < pty->n_subdomains; i++)
    BFT_FREE(pty->defs[i].ml_name);
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

  bft_printf("  %s >> uniform [%s], steady [%s], ",
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

  /* Definition */
  bft_printf("  %s >> n_subdomains    %d\n", pty->name, pty->n_subdomains);

  for (int i = 0; i < pty->n_subdomains; i++) {

    cs_param_def_t  *pdef  = pty->defs + i;

    bft_printf("  %s >> location  %s,", pty->name, pdef->ml_name);

    switch (pdef->def_type) {

    case CS_PARAM_DEF_BY_VALUE:
      {
        const cs_get_t  mat = pdef->def.get;

        switch(pty->type) {

        case CS_PROPERTY_ISO:
          bft_printf(" definition by value: % 5.3e\n", mat.val);
          break;
        case CS_PROPERTY_ORTHO:
          bft_printf(" definition by value: (% 5.3e, % 5.3e, % 5.3e)\n",
                     mat.vect[0], mat.vect[1], mat.vect[2]);
          break;
        case CS_PROPERTY_ANISO:
          bft_printf("\n                       |% 5.3e, % 5.3e, % 5.3e|\n"
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
      bft_printf("  definition by an analytical function\n");
      break;

    case CS_PARAM_DEF_BY_LAW_ONESCA:
      bft_printf("  definition by a law based on one scalar\n");
      break;

    case CS_PARAM_DEF_BY_LAW_SCAVEC:
      bft_printf("  definition by law based on one scalar + one vector\n");
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                _(" Invalid type of definition for a property."));
      break;

    } /* switch on def_type */

  } // Loop on subdomains

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
cs_property_def_by_law(cs_property_t             *pty,
                       const char                *ml_name,
                       const void                *context,
                       cs_onevar_law_func_t      *func)
{
  cs_param_def_t  *new_def = _init_new_def(pty, ml_name);

  new_def->def_type = CS_PARAM_DEF_BY_LAW_ONESCA;
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

  new_def->def_type = CS_PARAM_DEF_BY_LAW_TWOSCA;
  new_def->def.law2_func = func;
  new_def->context = context;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_property_t structure thanks to a law depending on
 *         a scalar and a vector variables in a subdomain attached to the mesh
 *         location named ml_name
 *
 * \param[in, out]  pty       pointer to a cs_property_t structure
 * \param[in]       ml_name   name of the related mesh location
 * \param[in]       context   pointer to a structure (may be NULL)
 * \param[in]       func      pointer to a law function defined by subdomain
 */
/*----------------------------------------------------------------------------*/

void
cs_property_def_by_scavec_law(cs_property_t             *pty,
                              const char                *ml_name,
                              const void                *context,
                              cs_scavec_law_func_t      *func)
{
  cs_param_def_t  *new_def = _init_new_def(pty, ml_name);

  new_def->def_type = CS_PARAM_DEF_BY_LAW_SCAVEC;
  new_def->def.law_scavec_func = func;
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
  cs_param_def_t  *new_def = _init_new_def(pty, N_("cells"));

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
  cs_get_t  get;

  if (pty == NULL)
    return;

  if (property_ts_id > -1)
    cs_timer_stats_start(property_ts_id);

  /* Initialize extra-diag. values of the tensor */
  for (int k = 0; k < 3; k++)
    for (int l = k+1; l < 3; l++)
      tensor[k][l] = 0;

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
      sub->def.analytic(t_cur, xc, &get);
      _get_tensor_by_value(pty, get, tensor);

    }
    break;

  case CS_PARAM_DEF_BY_LAW_ONESCA:
    _get_result_by_onevar_law(c_id, pty, sub->def.law1_func, sub->context,
                              &get);
    _get_tensor_by_value(pty, get, tensor);
    break;

  case CS_PARAM_DEF_BY_LAW_TWOSCA:
    _get_result_by_twovar_law(c_id, pty, sub->def.law2_func, sub->context,
                              &get);
    _get_tensor_by_value(pty, get, tensor);
    break;

  case CS_PARAM_DEF_BY_LAW_SCAVEC:
    _get_result_by_scavec_law(c_id, pty, sub->def.law_scavec_func, sub->context,
                              &get);
    _get_tensor_by_value(pty, get, tensor);
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
  bft_printf("\n  Tensor property for cell %d\n"
             "   | % 10.6e  % 10.6e  % 10.6e |\n"
             "   | % 10.6e  % 10.6e  % 10.6e |\n"
             "   | % 10.6e  % 10.6e  % 10.6e |\n", c_id,
             tensor[0][0], tensor[0][1], tensor[0][2],
             tensor[1][0], tensor[1][1], tensor[1][2],
             tensor[2][0], tensor[2][1], tensor[2][2]);
#endif

  if (property_ts_id > -1)
    cs_timer_stats_stop(property_ts_id);
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

  if (property_ts_id > -1)
    cs_timer_stats_start(property_ts_id);

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
    {
      const cs_real_t  *xc = cs_cdo_quant->cell_centers + 3*c_id;
      const double  t_cur = cs_time_step->t_cur;

      /* Call the analytic function. result is stored in get */
      sub->def.analytic(t_cur, xc, &get);
      result = get.val;
    }
    break;

  case CS_PARAM_DEF_BY_LAW_ONESCA:
    _get_result_by_onevar_law(c_id, pty, sub->def.law1_func, sub->context,
                              &get);
    result = get.val;
    break;

  case CS_PARAM_DEF_BY_LAW_TWOSCA:
    _get_result_by_twovar_law(c_id, pty, sub->def.law2_func, sub->context,
                              &get);
    result = get.val;
    break;

  case CS_PARAM_DEF_BY_ARRAY:
    result = _get_cell_value_from_array(c_id, pty->desc1, pty->array1);
    break;

  case CS_PARAM_DEF_BY_LAW_SCAVEC:
    _get_result_by_scavec_law(c_id,
                              pty, sub->def.law_scavec_func, sub->context,
                              &get);
    result = get.val;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the cell value related to property %s.\n"
              " Type of definition not handled yet.", pty->name);
    break;

  } /* type of definition */

  if (property_ts_id > -1)
    cs_timer_stats_stop(property_ts_id);

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

      const cs_real_t  hc = pow(cdoq->cell_vol[c_id], cs_math_onethird);

      fourier[c_id] = dt * ptyval / (hc * hc);

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

      const cs_real_t  hc = pow(cdoq->cell_vol[c_id], cs_math_onethird);

      fourier[c_id] = dt * eig_max / (hc * hc);

    } // Loop on cells

  } // Type of property

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
