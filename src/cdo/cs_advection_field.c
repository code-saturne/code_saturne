/*============================================================================
 * Manage the definition/setting of advection fields
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_post.h"
#include "cs_reco.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_ADVECTION_FIELD_DBG  0

/* Redefined names of function from cs_math to get shorter names */
#define _dp3 cs_math_3_dot_product

#define CS_ADVECTION_FIELD_POST_FIELD   (1 << 0)  // postprocess adv. field
#define CS_ADVECTION_FIELD_POST_COURANT (1 << 1)  // postprocess Courant number

struct _cs_adv_field_t {

  char  *restrict name;

  cs_desc_t  desc;          // Short descriptor (mask of bits)
  cs_flag_t  post_flag;     // Short descriptor dedicated to postprocessing
  int        vtx_field_id;  // id among cs_field_t structures (-1 if not used)
  int        cell_field_id; // id among cs_field_t structures (-1 if not used)

  cs_param_def_type_t  def_type;  // type of definition used
  cs_def_t             def;       // structure used for simple definition

  /* Useful buffers to deal with more complex definitions
     N.B.: array and struc are not owned by this structure.  */
  cs_desc_t    array_desc;  // short description of the related array

  const cs_real_t  *array;  // if the advection field hinges on an array
  const void       *struc;  // if the advection field hinges on a structure

};

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_adv[] =
  " Stop setting an empty cs_adv_field_t structure.\n"
  " Please check your settings.\n";

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;
static const cs_time_step_t  *cs_time_step;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

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
cs_advection_field_set_shared_pointers(const cs_cdo_quantities_t    *quant,
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
 * \brief  Create and initialize a new advection field structure
 *
 * \param[in]  name        name of the advection field
 *
 * \return a pointer to a new allocated cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_create(const char   *name)
{
  cs_adv_field_t  *adv = NULL;

  BFT_MALLOC(adv, 1, cs_adv_field_t);

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(adv->name, len, char);
  strncpy(adv->name, name, len);

  /* Default initialization */
  adv->desc.location = adv->desc.state = 0;
  adv->post_flag = 0;
  adv->vtx_field_id = -1;
  adv->cell_field_id = -1;
  adv->def_type = CS_PARAM_N_DEF_TYPES;
  adv->def.get.val = 0;

  /* If a complex definition is requested */
  adv->array_desc.location = adv->array_desc.state = 0;
  // adv->array and adv->struc

  return adv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_adv_field_t structure
 *
 * \param[in, out]  adv      pointer to a cs_adv_field_t structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_free(cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return adv;

  BFT_FREE(adv->name);
  BFT_FREE(adv);

  /* All other pointers are shared */

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the given advection field has the name ref_name
 *
 * \param[in]  adv         pointer to a cs_adv_field_t structure to test
 * \param[in]  ref_name    name of the advection field to find
 *
 * \return true if the name of the advection field is ref_name otherwise false
 */
/*----------------------------------------------------------------------------*/

bool
cs_advection_field_check_name(const cs_adv_field_t   *adv,
                              const char             *ref_name)
{
  if (adv == NULL)
    return false;

  int  reflen = strlen(ref_name);
  int  len = strlen(adv->name);

  if (reflen == len) {
    if (strcmp(ref_name, adv->name) == 0)
      return true;
    else
      return false;
  }
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  returns true if the advection field is uniform, otherwise false
 *
 * \param[in]    adv    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_advection_field_is_uniform(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return false;

  if (adv->desc.state & CS_FLAG_STATE_UNIFORM)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  returns true if the advection field is uniform in each cell
 *         otherwise false
 *
 * \param[in]    adv    pointer to a property to test
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_advection_field_is_cellwise(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return false;

  if (adv->desc.state & CS_FLAG_STATE_UNIFORM)
    return true;
  if (adv->desc.state & CS_FLAG_STATE_CELLWISE)
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of an advection field
 *
 * \param[in]    adv    pointer to an advection field structure
 *
 * \return  the name of the related advection field
 */
/*----------------------------------------------------------------------------*/

const char *
cs_advection_field_get_name(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return NULL;

  return adv->name;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the type of definition used to set the current advection
 *         field structure
 *
 * \param[in]    adv    pointer to an advection field structure
 *
 * \return  the type of definition
 */
/*----------------------------------------------------------------------------*/

cs_param_def_type_t
cs_advection_field_get_deftype(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return CS_PARAM_N_DEF_TYPES;

  return adv->def_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print a summary of a cs_adv_field_t structure
 *
 * \param[in]  adv      pointer to a cs_adv_field_t structure to summarize
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_summary(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return;

  _Bool  is_uniform = false, is_steady = true;

  if (adv->desc.state & CS_FLAG_STATE_UNIFORM)  is_uniform = true;
  if (adv->desc.state & CS_FLAG_STATE_UNSTEADY) is_steady = false;

  cs_log_printf(CS_LOG_SETUP,
                "  %s >> uniform [%s], steady [%s], ", adv->name,
                cs_base_strtf(is_uniform), cs_base_strtf(is_steady));

  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    {
      const cs_get_t  get = adv->def.get;
      cs_log_printf(CS_LOG_SETUP, "value: (% 5.3e, % 5.3e, % 5.3e)\n",
                    get.vect[0], get.vect[1], get.vect[2]);
    }
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    cs_log_printf(CS_LOG_SETUP, "definition by an analytical function\n");
    break;

  case CS_PARAM_DEF_BY_ARRAY:
    cs_log_printf(CS_LOG_SETUP, "definition by an array \n");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid type of definition for a material property."));
    break;

  } /* switch on def_type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set optional parameters related to a cs_adv_field_t structure
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       key       key related to the member of adv to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_set_option(cs_adv_field_t            *adv,
                              cs_advection_field_key_t   key,
                              const char                *keyval)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {

  case CS_ADVKEY_DEFINE_AT:
    if (strcmp(val,"cells") == 0)
      adv->desc.location |= CS_FLAG_CELL;
    else if (strcmp(val, "vertices") == 0)
      adv->desc.location |= CS_FLAG_VERTEX;
    else
      bft_error(__FILE__, __LINE__, 0,
                N_(" Invalid value %s for setting key CS_ADVKEY_DEFINE_AT\n"
                   " Valid choices are \"cells\" or \"vertices\".\n"
                   " Please modify your setting.\n"), val);
    break;

  case CS_ADVKEY_POST:
    if (strcmp(val, "field") == 0)
      adv->post_flag |= CS_ADVECTION_FIELD_POST_FIELD;
    else if (strcmp(val, "courant") == 0)
      adv->post_flag |= CS_ADVECTION_FIELD_POST_COURANT;
    else
      bft_error(__FILE__, __LINE__, 0,
                N_(" Invalid value %s for setting key CS_ADVKEY_POST\n"
                   " Valid choices are \"field\" or \"courant\".\n"
                   " Please modify your setting.\n"), val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key not implemented for setting an advection field."));

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of a cs_adv_field_t structure
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       val       accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_value(cs_adv_field_t    *adv,
                                const char        *val)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  adv->def_type = CS_PARAM_DEF_BY_VALUE;
  adv->desc.state |= CS_FLAG_STATE_UNIFORM;

  cs_param_set_def(adv->def_type, CS_PARAM_VAR_VECT, (const void *)val,
                   &(adv->def));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an analytic function
 *
 * \param[in, out]  adv     pointer to a cs_adv_field_t structure
 * \param[in]       func    pointer to a function
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_analytic(cs_adv_field_t        *adv,
                                   cs_analytic_func_t    *func)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  adv->def_type = CS_PARAM_DEF_BY_ANALYTIC_FUNCTION;
  adv->def.analytic = func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an array of values
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       desc      information about this array
 * \param[in]       array     pointer to an array
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_array(cs_adv_field_t     *adv,
                                cs_desc_t           desc,
                                const cs_real_t    *array)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  adv->def_type = CS_PARAM_DEF_BY_ARRAY;
  adv->array_desc.location = desc.location;
  adv->array_desc.state = desc.state;
  adv->array = array;

  if (cs_cdo_same_support(desc.location, cs_cdo_dual_face_byc) ||
      cs_cdo_same_support(desc.location, cs_cdo_primal_cell))
    adv->desc.state |= CS_FLAG_STATE_CELLWISE;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a cs_field_t structure related to an advection field
 *
 * \param[in, out] adv     pointer to a cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_create_field(cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return;

  int  len;
  char *field_name = NULL;

  _Bool has_previous = (adv->desc.state & CS_FLAG_STATE_UNSTEADY) ? true:false;
  int  field_mask = CS_FIELD_PROPERTY;

  if (adv->desc.location & CS_FLAG_VERTEX) { // Add a field attached to vertices

    /* Define the name of the field */
    len = strlen(adv->name) + strlen("_vertices") + 1;
    BFT_MALLOC(field_name, len, char);
    sprintf(field_name, "%s_vertices", adv->name);

    cs_field_t  *fld = cs_field_create(field_name,
                                       field_mask,
                                       CS_MESH_LOCATION_VERTICES,
                                       3,    // always a vector-valued field
                                       has_previous);

    adv->vtx_field_id = cs_field_id_by_name(field_name);

    /* Allocate and initialize values */
    cs_field_allocate_values(fld);

    BFT_FREE(field_name);

  } // Add a field attached to vertices

  if (adv->desc.location & CS_FLAG_CELL) { // Add a field attached to cells

    /* Define the name of the field */
    len = strlen(adv->name) + strlen("_cells") + 1;
    BFT_MALLOC(field_name, len, char);
    sprintf(field_name, "%s_cells", adv->name);

    cs_field_t  *fld = cs_field_create(field_name,
                                       field_mask,
                                       CS_MESH_LOCATION_CELLS,
                                       3,    // always a vector-valued field
                                       has_previous);

    adv->cell_field_id = cs_field_id_by_name(field_name);

    /* Allocate and initialize values */
    cs_field_allocate_values(fld);

    BFT_FREE(field_name);

  } // Add a field attached to vertices

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at the cell center
 *
 * \param[in]      c_id    id of the current cell
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 * \param[in, out] vect    pointer to a cs_nvec3_t structure (meas + unitv)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_cell_vector(cs_lnum_t               c_id,
                                   const cs_adv_field_t   *adv,
                                   cs_nvec3_t             *vect)
{
  /* Initialize the vector */
  vect->meas = 0.;
  for (int k = 0; k < 3; k++)
    vect->unitv[k] = 0;

  if (adv == NULL)
    return;

  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    cs_nvec3(adv->def.get.vect, vect);
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    {
      cs_get_t  get;

      const cs_real_t  *xc = cs_cdo_quant->cell_centers + 3*c_id;
      const double  t_cur = cs_time_step->t_cur;

      /* Call the analytic function. result is stored in get */
      adv->def.analytic(t_cur, xc, &get);

      /* Switch to a cs_nvec3_t representation */
      cs_nvec3(get.vect, vect);
    }
    break;

  case CS_PARAM_DEF_BY_ARRAY:
    {
      cs_real_3_t  val;

      /* Test if location has at least the pattern of the reference support */
      if (cs_cdo_same_support(adv->array_desc.location, cs_cdo_dual_face_byc))
        cs_reco_dfbyc_at_cell_center(c_id,
                                     cs_cdo_connect->c2e,
                                     cs_cdo_quant,
                                     adv->array, val);

      else if (cs_cdo_same_support(adv->array_desc.location,
                                   cs_cdo_primal_cell)) {

        for (int k = 0; k < 3; k++)
          val[k] = adv->array[3*c_id+k];

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid support for evaluating the advection field %s"
                  " at the cell center of cell %d.", adv->name, c_id);

      /* Switch to a cs_nvec3_t representation */
      cs_nvec3(val, vect);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the vector field for cell id %d related to"
              " the advection field %s.\n"
              " Type of definition not handled yet.", c_id, adv->name);
    break;

  } /* type of definition */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field for a given face
 *
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 * \param[in]      xyz     coordinates where to evaluate the advection field
 * \param[in, out] vect    pointer to a cs_nvec3_t structure (meas + unitv)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_at_xyz(const cs_adv_field_t   *adv,
                              const cs_real_3_t       xyz,
                              cs_nvec3_t             *vect)
{
  /* Initialize the vector */
  vect->meas = 0.;
  for (int k = 0; k < 3; k++)
    vect->unitv[k] = 0;

  if (adv == NULL)
    return;

  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    cs_nvec3(adv->def.get.vect, vect);
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    {
      cs_get_t  get;

      /* Call the analytic function. result is stored in get */
      adv->def.analytic(cs_time_step->t_cur, xyz, &get);

      /* Switch to a cs_nvec3_t representation */
      cs_nvec3(get.vect, vect);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the vector field at face centers related to"
              " the advection field %s.\n"
              " Type of definition not handled yet.", adv->name);
    break;

  } /* type of definition */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at cell centers
 *
 * \param[in]      adv           pointer to a cs_adv_field_t structure
 * \param[in, out] cell_values   array of values at cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_at_cells(const cs_adv_field_t  *adv,
                            cs_real_t             *cell_values)
{
  if (adv == NULL)
    return;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    { // Uniform value inside the computational domain
      const cs_get_t  get = adv->def.get;

      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        const cs_lnum_t  shift_c = 3*c_id;

        cell_values[shift_c]    = get.vect[0];
        cell_values[shift_c +1] = get.vect[1];
        cell_values[shift_c +2] = get.vect[2];

      } // Loop on cells

    }
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    {
      cs_get_t  get;

      const double  t_cur = cs_time_step->t_cur;

      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        const cs_lnum_t  shift_c = 3*c_id;
        const cs_real_t  *xc = quant->cell_centers + shift_c;

        /* Call the analytic function. result is stored in get */
        adv->def.analytic(t_cur, xc, &get);

        cell_values[shift_c]    = get.vect[0];
        cell_values[shift_c +1] = get.vect[1];
        cell_values[shift_c +2] = get.vect[2];

      } // Loop on cells

    }
    break;

  case CS_PARAM_DEF_BY_ARRAY:
    {
      cs_real_3_t  recoval;

      /* Test if location has at least the pattern of the reference support */
      if (cs_cdo_same_support(adv->array_desc.location, cs_cdo_dual_face_byc)) {

        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

          const cs_lnum_t  shift_c = 3*c_id;

          cs_reco_dfbyc_at_cell_center(c_id,
                                       cs_cdo_connect->c2e,
                                       quant,
                                       adv->array, recoval);

          cell_values[shift_c]    = recoval[0];
          cell_values[shift_c +1] = recoval[1];
          cell_values[shift_c +2] = recoval[2];

        } // Loop on cells

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid support for evaluating the advection field %s"
                  " at cell centers.", adv->name);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the vector field at cell centers related to"
              " the advection field %s.\n"
              " Type of definition not handled yet.", adv->name);
    break;

  } /* type of definition */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at vertices
 *
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 * \param[in, out] vect    pointer to a cs_nvec3_t structure (meas/unitv)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_at_vertices(const cs_adv_field_t  *adv,
                               cs_real_t             *vtx_values)
{
  if (adv == NULL)
    return;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;

  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    { // Uniform value inside the computational domain
      const cs_get_t  get = adv->def.get;

      for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {

        const cs_lnum_t  shift_v = 3*v_id;

        vtx_values[shift_v]     = get.vect[0];
        vtx_values[shift_v + 1] = get.vect[1];
        vtx_values[shift_v + 2] = get.vect[2];

      } // Loop on vertices

    }
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    {
      cs_get_t  get;

      const double  t_cur = cs_time_step->t_cur;

      for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {

        const cs_lnum_t  shift = 3*v_id;
        const cs_real_t  *xv = quant->vtx_coord + shift;

        /* Call the analytic function. result is stored in get */
        adv->def.analytic(t_cur, xv, &get);

        /* Store computed values */
        for (int k = 0; k < 3; k++)
          vtx_values[shift+k] = get.vect[k];

      } // Loop on vertices

    }
    break;

  case CS_PARAM_DEF_BY_ARRAY:
    {
      cs_real_3_t  val;
      double  *dc_vol = NULL;

      BFT_MALLOC(dc_vol, quant->n_vertices, double);
      for (cs_lnum_t j = 0; j < quant->n_vertices; j++) {
        dc_vol[j] = 0;
        vtx_values[3*j] = vtx_values[3*j+1] = vtx_values[3*j+2] = 0.0;
      }

      const cs_connect_index_t  *c2v = cs_cdo_connect->c2v;

      /* Test if flag has at least the pattern of the reference support */
      if (cs_cdo_same_support(adv->array_desc.location, cs_cdo_dual_face_byc)) {

        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

          cs_reco_dfbyc_at_cell_center(c_id,
                                       cs_cdo_connect->c2e,
                                       quant,
                                       adv->array, val);

          for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

            const cs_lnum_t  v_id = c2v->ids[j];
            const cs_lnum_t  shift_v = 3*v_id;
            const double  dcc_vol = quant->dcell_vol[j];

            dc_vol[v_id] += dcc_vol;
            for (int k = 0; k < 3; k++)
              vtx_values[shift_v + k] += dcc_vol * val[k];

          } // Loop on cell vertices

        } // Loop on cells

      }
      else if (cs_cdo_same_support(adv->array_desc.location,
                                   cs_cdo_primal_cell)) {

        for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
          for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

            const cs_lnum_t  v_id = c2v->ids[j];
            const cs_lnum_t  shift_v = 3*v_id;
            const double  dcc_vol = quant->dcell_vol[j];

            dc_vol[v_id] += dcc_vol;
            for (int k = 0; k < 3; k++)
              vtx_values[shift_v + k] += dcc_vol * adv->array[3*c_id+k];

          } // Loop on cell vertices
        } // Loop on cells

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid support for evaluating the advection field %s"
                  " at vertices.", adv->name);

      for (cs_lnum_t j = 0; j < quant->n_vertices; j++) {

        const double  inv_dcvol = 1/dc_vol[j];
        const cs_lnum_t  shift_v = 3*j;

        for (int k = 0; k < 3; k++)
          vtx_values[shift_v+k] *= inv_dcvol;

      } // Loop on vertices

      /* Free temporary buffer */
      BFT_FREE(dc_vol);
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the vector field at vertices related to"
              " the advection field %s.\n"
              " Type of definition not handled yet.", adv->name);
    break;

  } /* type of definition */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the max. value of the advection field among cells
 *
 * \param[in]      adv       pointer to a cs_adv_field_t structure
 *
 * \return the  max. value of the magnitude of the field at cell centers
 */
/*----------------------------------------------------------------------------*/

double
cs_advection_field_get_cell_max(const cs_adv_field_t      *adv)
{
  double  beta_max = 0.0;

  if (adv == NULL)
    return beta_max;

  const cs_cdo_quantities_t  *quant = cs_cdo_quant;
  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    beta_max = cs_math_3_norm(adv->def.get.vect);
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
  case CS_PARAM_DEF_BY_ARRAY:
    {
      cs_nvec3_t  adv_cell;

      for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

        cs_advection_field_get_cell_vector(c_id, adv, &adv_cell);
        beta_max = fmax(beta_max, adv_cell.meas);

      } // Loop on cells

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Stop computing the max. ratio of the vector field inside each"
              " cell for field called %s\n"
              " Type of definition not handled yet.", adv->name);
    break;

  } /* type of definition */

  return beta_max;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         the dual faces of a cell
 *
 * \param[in]      c_id     id of the current cell
 * \param[in]      a_info   set of parameters for the advection operator
 * \param[in]      adv      pointer to a cs_adv_field_t structure
 * \param[in, out] fluxes   array of values attached to dual faces of a cell
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_flux_dfaces(cs_lnum_t                     c_id,
                                   const cs_param_advection_t    a_info,
                                   const cs_adv_field_t         *adv,
                                   cs_real_t                    *fluxes)
{
  if (adv == NULL)
    return;

  if (fluxes == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " fluxes array should be allocated before the call.");

  cs_lnum_t  _i, jf, je, k;
  cs_nvec3_t  adv_vect;

  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_connect_index_t  *c2e = connect->c2e;

  if (adv->desc.state & CS_FLAG_STATE_UNIFORM ||
      adv->desc.state & CS_FLAG_STATE_CELLWISE) {

    cs_advection_field_get_cell_vector(c_id, adv, &adv_vect);

    /* Loop on cell edges */
    for (je = c2e->idx[c_id], _i = 0; je < c2e->idx[c_id+1]; je++, _i++) {

      const cs_dface_t  qdf = cdoq->dface[je];

      fluxes[_i] = adv_vect.meas * _dp3(adv_vect.unitv, qdf.vect);

    } // Loop on cell edges

  } /* Uniform in this cell */
  else {

    switch (adv->def_type) {

    case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
      {
        cs_real_t  w;
        cs_real_3_t  gpts[3], xg;
        cs_get_t  get;

        const double  t_cur = cs_time_step->t_cur;
        const cs_real_t  *xc = cdoq->cell_centers + 3*c_id;

        /* Loop on cell edges */
        for (je = c2e->idx[c_id], _i = 0; je < c2e->idx[c_id+1]; je++, _i++) {

          const cs_dface_t  qdf = cdoq->dface[je];
          const cs_lnum_t  e_id = c2e->ids[je];
          const cs_quant_t  qe = cdoq->edge[e_id];

          fluxes[_i] = 0.;

          for (jf = 0; jf < 2; jf++) {

            const cs_nvec3_t  tef = qdf.sface[jf];
            const cs_lnum_t  f_id = qdf.parent_id[jf];
            const cs_quant_t  qf = cdoq->face[f_id];

            switch (a_info.quad_type) {

            case CS_QUADRATURE_BARY:
            case CS_QUADRATURE_BARY_SUBDIV:
              for (k = 0; k < 3; k++)
                xg[k] = cs_math_onethird *(xc[k] + qe.center[k] + qf.center[k]);
              adv->def.analytic(t_cur, xg, &get);
              fluxes[_i] += tef.meas * _dp3(get.vect, tef.unitv);
              break;

            case CS_QUADRATURE_HIGHER:
              {
                cs_quadrature_tria_3pts(qe.center, qf.center, xc, tef.meas,
                                        gpts, &w);

                cs_real_t  add = 0;
                for (int p = 0; p < 3; p++) {
                  adv->def.analytic(t_cur, gpts[p], &get);
                  add += _dp3(get.vect, tef.unitv);
                }
                fluxes[_i] += add * w;

              }
              break;

            case CS_QUADRATURE_HIGHEST: // Not yet implemented
            default:
              bft_error(__FILE__, __LINE__, 0, " Invalid type of quadrature.");
              break;

            } // switch type of quadrature

          } // Loop on the two triangles composing the dual face inside a cell

        } // Loop on cell edges

      } // definition by analytic function
      break;

    case CS_PARAM_DEF_BY_ARRAY:
      {
        /* Test if location has at least the pattern of the reference support */
        if (cs_cdo_same_support(adv->array_desc.location, cs_cdo_dual_face_byc))
          for (je = c2e->idx[c_id], _i = 0; je < c2e->idx[c_id+1]; je++, _i++)
            fluxes[_i] = adv->array[je];
        else
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid support for evaluating the advection field %s"
                    " at the cell center of cell %d.", adv->name, c_id);
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, "Incompatible type of definition.");
      break;

    } // def_type

  } /* Not uniform in this cell */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         triangle defined by the two vertices of an edge and the barycenter
 *         of a face.
 *
 * \param[in]  adv       pointer to a cs_adv_field_t structure
 * \param[in]  a_info    set of parameters for the advection operator
 * \param[in]  cm        pointer to a cs_cell_mesh_t structure
 * \param[in]  f         id of the face in the current cell
 * \param[in]  e         id of the edge in the current cell
 * \param[in]  v1        id of the first vertex in the current cell
 * \param[in]  v2        id of the second vertex in the current cell
 *
 * \return the value of the flux across tef
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_advection_field_get_flux_tef(const cs_adv_field_t        *adv,
                                const cs_param_advection_t   a_info,
                                const cs_cell_mesh_t        *cm,
                                short int                    f,
                                short int                    e,
                                short int                    v1,
                                short int                    v2)
{
  cs_real_t  adv_flx = 0;

  if (adv == NULL)
    return adv_flx;

  const cs_real_t  *xv1 = cm->xv + 3*v1;
  const cs_real_t  *xv2 = cm->xv + 3*v2;
  const cs_quant_t  pfq = cm->face[f];

  cs_real_t  tef = cs_math_surftri(xv1, xv2, pfq.center);

  /* Compute the flux accros the portion of primal face */
  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    adv_flx = tef * _dp3(adv->def.get.vect, pfq.unitv);
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    {
      cs_get_t  get;

      const double  t_cur = cs_time_step->t_cur;

      switch (a_info.quad_type) {

      case CS_QUADRATURE_BARY:
      case CS_QUADRATURE_BARY_SUBDIV:
        {
          cs_real_3_t  xg;

          for (int k = 0; k < 3; k++)
            xg[k] = cs_math_onethird * (xv1[k] + xv2[k] + pfq.center[k]);

          /* Call the analytic function. result is stored in get */
          adv->def.analytic(t_cur, xg, &get);
          adv_flx = tef * _dp3(get.vect, pfq.unitv);
        }
        break;

      case CS_QUADRATURE_HIGHER:
        {
          cs_real_t  w, add = 0.;
          cs_real_3_t  gpts[3];

          cs_quadrature_tria_3pts(xv1, xv2, pfq.center, tef, gpts, &w);

          for (int p = 0; p < 3; p++) {
            adv->def.analytic(t_cur, gpts[p], &get);
            add += _dp3(get.vect, pfq.unitv);
          }
          adv_flx += add * w;
        }
        break;

      case CS_QUADRATURE_HIGHEST: // Not yet implemented
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid type of quadrature for computing the flux of %s"
                  " across an elementary triangle s(v,e,f).\n"
                  " This functionality is not implemented yet.", adv->name);
        break;
      }

    }
    break; // DEF_ANALYTIC_FUNCTION

    case CS_PARAM_DEF_BY_ARRAY:
      {
        cs_real_3_t  reco;

        /* Test if flag has at least the pattern of the reference support */
        if (cs_cdo_same_support(adv->array_desc.location,
                                cs_cdo_dual_face_byc)) {

          const cs_connect_index_t  *c2e = cs_cdo_connect->c2e;

          /* Compute the reconstruction of the flux in pec */
          cs_reco_dfbyc_in_pec(cm->c_id, cm->e_ids[e], c2e, cs_cdo_quant,
                               adv->array,
                               reco);

          /* The reconstruction yields a constant vector field */
          adv_flx = tef * _dp3(pfq.unitv, reco);

        }
        else
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid support for evaluating the advection field %s"
                    " across tef.", adv->name);

      } // DEF_BY_ARRAY
      break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of definition for computing the flux of %s"
              " across the triangle tef.\n"
              " This functionality is not implemented yet.", adv->name);
    break;

  } // switch def_type

  return adv_flx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         triangle defined by a vertex, the face and edge barycenters
 *
 * \param[in]  v_id     id of the current vertex
 * \param[in]  e_id     id of the current edge
 * \param[in]  f_id     id of the current face
 * \param[in]  a_info   set of parameters for the advection operator
 * \param[in]  adv      pointer to a cs_adv_field_t structure
 *
 * \return the value of the flux across s(v,e,f)
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_advection_field_get_flux_svef(cs_lnum_t                    v_id,
                                 cs_lnum_t                    e_id,
                                 cs_lnum_t                    f_id,
                                 const cs_param_advection_t   a_info,
                                 const cs_adv_field_t        *adv)
{
  int  k;

  cs_real_t  adv_flx = 0;

  if (adv == NULL)
    return adv_flx;

  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;
  const cs_quant_t  pfq = cdoq->face[f_id];
  const cs_quant_t  peq = cdoq->edge[e_id];
  const cs_real_t  *xv = cdoq->vtx_coord + 3*v_id;

  cs_real_t  surf = cs_math_surftri(xv, peq.center, pfq.center);

  /* Compute the flux accros the portion of primal face */
  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    adv_flx = surf * _dp3(adv->def.get.vect, pfq.unitv);
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    {
      cs_real_t  w;
      cs_real_3_t  gpts[3], xg;
      cs_get_t  get;

      const double  t_cur = cs_time_step->t_cur;

      switch (a_info.quad_type) {

      case CS_QUADRATURE_BARY:
      case CS_QUADRATURE_BARY_SUBDIV:
        for (k = 0; k < 3; k++)
          xg[k] = cs_math_onethird * (xv[k] + peq.center[k] + pfq.center[k]);

        /* Call the analytic function. result is stored in get */
        adv->def.analytic(t_cur, xg, &get);
        adv_flx = surf * _dp3(get.vect, pfq.unitv);
        break;

      case CS_QUADRATURE_HIGHER:
        {
          cs_quadrature_tria_3pts(peq.center, pfq.center, xv, surf, gpts, &w);

          cs_real_t  add = 0;
          for (int p = 0; p < 3; p++) {
            adv->def.analytic(t_cur, gpts[p], &get);
            add += _dp3(get.vect, pfq.unitv);
          }
          adv_flx += add * w;
        }
        break;

      case CS_QUADRATURE_HIGHEST: // Not yet implemented
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid type of quadrature for computing the flux of %s"
                  " across an elementary triangle s(v,e,f).\n"
                  " This functionality is not implemented yet.", adv->name);
        break;
      }

    }
    break; // DEF_ANALYTIC_FUNCTION

    case CS_PARAM_DEF_BY_ARRAY:
      {
        cs_real_3_t  reco;

        /* Test if flag has at least the pattern of the reference support */
        if (cs_cdo_same_support(adv->array_desc.location,
                                cs_cdo_dual_face_byc)) {

          const cs_connect_index_t  *c2e = cs_cdo_connect->c2e;
          const cs_sla_matrix_t  *f2c = cs_cdo_connect->f2c;

          cs_lnum_t  c_id = f2c->col_id[f2c->idx[f_id]];
          assert(f2c->idx[f_id+1] - f2c->idx[f_id] == 1);

          /* Compute the reconstruction of the flux in pec */
          cs_reco_dfbyc_in_pec(c_id, e_id, c2e, cdoq, adv->array, reco);

          /* The reconstruction yields a constant vector field */
          adv_flx = surf * _dp3(pfq.unitv, reco);

        }
        else
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid support for evaluating the advection field %s"
                    " across s(v,e,f).", adv->name);

      } // DEF_BY_ARRAY
      break;


  default:
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of definition for computing the flux of %s"
              " across an elementary triangle s(v,e,f).\n"
              " This functionality is not implemented yet.", adv->name);
    break;

  } // switch def_type

  return adv_flx;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the values of the related field(s)
 *
 * \param[in, out]     adv     pointer to a cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_update(cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return;

  if (adv->vtx_field_id > -1) { /* Field stored at vertices */

    cs_field_t  *fld = cs_field_by_id(adv->vtx_field_id);

    /* Copy current field values to previous values */
    cs_field_current_to_previous(fld);

    /* Set new values */
    cs_advection_field_at_vertices(adv, fld->val);

  }

  if (adv->cell_field_id > -1) { /* Field stored at cell centers */

    cs_field_t  *fld = cs_field_by_id(adv->cell_field_id);

    /* Copy current field values to previous values */
    cs_field_current_to_previous(fld);

    /* Set new values */
    cs_advection_field_at_cells(adv, fld->val);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Peclet number in each cell
 *
 * \param[in]      adv        pointer to the advection field struct.
 * \param[in]      diff       pointer to the diffusion property struct.
 * \param[in, out] peclet     pointer to an array storing Peclet number
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_get_peclet(const cs_adv_field_t        *adv,
                        const cs_property_t         *diff,
                        cs_real_t                    peclet[])
{
  cs_real_t  ptymat[3][3];
  cs_real_3_t  ptydir;
  cs_nvec3_t  adv_c;

  const bool  pty_uniform = cs_property_is_uniform(diff);
  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;

  assert(peclet != NULL); // Sanity check

  /* Get the value of the material property at the first cell center */
  if (pty_uniform)
    cs_property_get_cell_tensor(0, diff, false, ptymat);

  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    /* Get the value of the material property at the cell center */
    if (!pty_uniform)
      cs_property_get_cell_tensor(c_id, diff, false, ptymat);

    cs_advection_field_get_cell_vector(c_id, adv, &adv_c);

    cs_real_t  hc = pow(cdoq->cell_vol[c_id], cs_math_onethird);

    cs_math_33_3_product((const cs_real_t (*)[3])ptymat, adv_c.unitv, ptydir);

    peclet[c_id] = hc * adv_c.meas / _dp3(adv_c.unitv, ptydir);

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined post-processing output for advection fields.
 *
 * \param[in]  adv         pointer to a cs_adv_field_t structure
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 * \param[in]  dt_cur      value of the current time step
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_extra_post(const cs_adv_field_t      *adv,
                              const cs_time_step_t      *time_step,
                              double                     dt_cur)
{
  if (adv == NULL)
    return;
  if (adv->post_flag == 0)
    return;

  const bool post_field =
    (adv->post_flag & CS_ADVECTION_FIELD_POST_FIELD) ? true : false;
  const bool post_courant =
    (adv->post_flag & CS_ADVECTION_FIELD_POST_COURANT) ? true : false;

  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;

#if CS_ADVECTION_FIELD_DBG && defined(DEBUG) && !defined(NDEBUG)
  cs_log_printf(CS_LOG_SETUP,
                " <post/advection_field %s> iter: %d\n",
                adv->name, time_step->nt_cur);
#endif

  /* Field is defined at vertices ? */
  cs_field_t  *fld = NULL;
  if (adv->vtx_field_id > -1)
    fld = cs_field_by_id(adv->vtx_field_id);

  if (fld != NULL && post_field)
    cs_post_write_vertex_var(CS_POST_MESH_VOLUME,
                             CS_POST_WRITER_ALL_ASSOCIATED,
                             fld->name,
                             3,               // dim
                             true,            // interlace
                             true,            // true = original mesh
                             CS_POST_TYPE_cs_real_t,
                             fld->val,        // values on vertices
                             time_step);      // time step structure

  /* Field is defined at cells ? */
  fld = NULL;
  if (adv->cell_field_id > -1)
    fld = cs_field_by_id(adv->cell_field_id);

  if (fld != NULL && post_field)
    cs_post_write_var(CS_POST_MESH_VOLUME,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      fld->name,
                      3,               // dim
                      true,            // interlace
                      true,            // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      fld->val,        // values in cells
                      NULL, NULL,      // values at internal/border faces
                      time_step);      // time step structure

  if (post_courant) { /* Compute and postprocess the Courant number */

    double  hc;
    cs_nvec3_t  adv_c;

    char  *label = NULL;
    double  *courant = NULL;
    int  len = strlen(adv->name) + 8 + 1;

    BFT_MALLOC(courant, cdoq->n_cells, double);
    BFT_MALLOC(label, len, char);
    sprintf(label, "%s.Courant", adv->name);

    if (fld != NULL) { /* field is defined at cell centers */
      for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {
        cs_nvec3(fld->val + 3*c_id, &adv_c);
        hc = pow(cdoq->cell_vol[c_id], cs_math_onethird);
        courant[c_id] = dt_cur * adv_c.meas / hc;
      }
    }
    else {
      for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {
        cs_advection_field_get_cell_vector(c_id, adv, &adv_c);
        hc = pow(cdoq->cell_vol[c_id], cs_math_onethird);
        courant[c_id] = dt_cur * adv_c.meas / hc;
      }
    }

    cs_post_write_var(CS_POST_MESH_VOLUME,
                      CS_POST_WRITER_ALL_ASSOCIATED,
                      label,
                      1,
                      true,           // interlace
                      true,           // true = original mesh
                      CS_POST_TYPE_cs_real_t,
                      courant,        // values on cells
                      NULL,           // values at internal faces
                      NULL,           // values at border faces
                      time_step);     // time step management struct.

    BFT_FREE(label);
    BFT_FREE(courant);
  }

}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
