/*============================================================================
 * Manage the definition/setting of advection fields
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_mesh_location.h"
#include "cs_field.h"
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

#define CS_ADVECTION_FIELD_DBG  1

struct _cs_adv_field_t {

  char  *restrict name;

  cs_flag_t  flag;       /* Short descriptor (mask of bits) */
  int     post_freq;     /* -1 (no post-processing), 0 (at the beginning)
                            otherwise every post_freq iteration(s) */
  bool    post_unitv;    /* Postprocess or not the unit vector field */
  int     vtx_field_id;  /* id among cs_field_t structures (-1 if not used) */
  int     cell_field_id; /* id among cs_field_t structures (-1 if not used) */

  cs_param_def_type_t    def_type;
  cs_def_t               def;

  /* Pointer to the main structures (not owned, only shared) */
  const cs_cdo_quantities_t   *cdoq;
  const cs_cdo_connect_t      *connect;
  const cs_time_step_t        *time_step;

  /* Useful buffers to deal with more complex definitions
     var and struc are not owned by this structure.  */
  cs_flag_t    array_flag;  // short description of the related array
  const cs_real_t  *array;  // if the advection field hinges on an array
  const void       *struc;  // if the advection field hinges on a structure

};

/* List of available keys for setting an advection field */
typedef enum {

  ADVKEY_POST_FREQ,
  ADVKEY_POST_UNITV,
  ADVKEY_CELL_FIELD,
  ADVKEY_VERTEX_FIELD,
  ADVKEY_ERROR

} advkey_t;

/*============================================================================
 * Private variables
 *============================================================================*/

static const double  one_third = 1./3.;
static const char _err_empty_adv[] =
  " Stop setting an empty cs_adv_field_t structure.\n"
  " Please check your settings.\n";

/* Array support is on dual faces and scan thanks to the c2e connectivity */
static const cs_flag_t  cs_var_support_dfbyc =
  CS_PARAM_FLAG_FACE | CS_PARAM_FLAG_DUAL | CS_PARAM_FLAG_BY_CELL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print the name of the corresponding advection key
 *
 * \param[in] key        name of the key
 *
 * \return a string
 */
/*----------------------------------------------------------------------------*/

static const char *
_print_advkey(advkey_t  key)
{
  switch (key) {

  case ADVKEY_POST_FREQ:
    return "post_freq";
  case ADVKEY_POST_UNITV:
    return "post_unitv";
  case ADVKEY_CELL_FIELD:
    return "cell_field";
  case ADVKEY_VERTEX_FIELD:
    return "vertex_field";
  default:
    assert(0);

  }

  return NULL; // avoid a warning
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the corresponding enum from the name of an advection key.
 *         If not found, return a key error.
 *
 * \param[in] keyname    name of the key
 *
 * \return a advkey_t
 */
/*----------------------------------------------------------------------------*/

static advkey_t
_get_advkey(const char  *keyname)
{
  advkey_t  key = ADVKEY_ERROR;

  if (strcmp(keyname, "post_freq") == 0)
    key = ADVKEY_POST_FREQ;
  else if (strcmp(keyname, "post_unitv") == 0)
    key = ADVKEY_POST_UNITV;
  else if (strcmp(keyname, "cell_field") == 0)
    key = ADVKEY_CELL_FIELD;
  else if (strcmp(keyname, "vertex_field") == 0)
    key = ADVKEY_VERTEX_FIELD;

  return key;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a new advection field structure
 *
 * \param[in]  name        name of the advection field
 * \param[in]  cdoq        pointer to a cs_cdo_quantities_t struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a cs_time_step_t struct.
 *
 * \return a pointer to a new allocated cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_create(const char                  *name,
                          const cs_cdo_quantities_t   *cdoq,
                          const cs_cdo_connect_t      *connect,
                          const cs_time_step_t        *time_step)
{
  cs_adv_field_t  *adv = NULL;

  BFT_MALLOC(adv, 1, cs_adv_field_t);

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(adv->name, len, char);
  strncpy(adv->name, name, len);

  /* Shared pointers for defining the adv_field */
  adv->cdoq = cdoq;
  adv->connect = connect;
  adv->time_step = time_step;

  /* Default initialization */
  adv->post_unitv = false;
  adv->post_freq = -1;
  adv->vtx_field_id = -1;
  adv->cell_field_id = -1;
  adv->flag = 0;
  adv->def_type = CS_PARAM_N_DEF_TYPES;
  adv->def.get.val = 0;

  /* If a complex definition is requested */
  adv->array_flag = 0;
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

  if (adv->flag & CS_PARAM_FLAG_UNIFORM)
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

  if (adv->flag & CS_PARAM_FLAG_UNIFORM)
    return true;
  if (adv->flag & CS_PARAM_FLAG_CELLWISE)
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

  if (adv->flag & CS_PARAM_FLAG_UNIFORM)  is_uniform = true;
  if (adv->flag & CS_PARAM_FLAG_UNSTEADY) is_steady = false;

  bft_printf(" %s >> uniform [%s], steady [%s], ",
             adv->name, cs_base_strtf(is_uniform), cs_base_strtf(is_steady));

  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    {
      const cs_get_t  get = adv->def.get;
      bft_printf("value: (% 5.3e, % 5.3e, % 5.3e)\n",
                 get.vect[0], get.vect[1], get.vect[2]);
    }
    break;

  case CS_PARAM_DEF_BY_ANALYTIC_FUNCTION:
    bft_printf("definition by an analytical function\n");
    break;

  case CS_PARAM_DEF_BY_ARRAY:
    bft_printf("definition by an array \n");
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
 * \param[in]       keyname   name of key related to the member of adv to set
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_set_option(cs_adv_field_t   *adv,
                              const char       *keyname,
                              const char       *keyval)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  advkey_t  key = _get_advkey(keyname);

  if (key == ADVKEY_ERROR) {

    bft_printf("\n\n Current key: %s\n", keyname);
    bft_printf(" Possible keys: ");
    for (int i = 0; i < ADVKEY_ERROR; i++) {
      bft_printf("%s ", _print_advkey(i));
      if (i > 0 && i % 3 == 0)
        bft_printf("\n\t");
    }
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid key for setting the advection field %s.\n"
                " Please read listing for more details and"
                " modify your settings."), adv->name);

  } /* Error message */

  switch(key) {

  case ADVKEY_POST_FREQ:
    adv->post_freq = atoi(keyval);
    break;
  case ADVKEY_POST_UNITV:
    if (!strcmp(keyval, "false"))
      adv->post_unitv = false;
    else
      adv->post_unitv = true;
    break;
  case ADVKEY_CELL_FIELD:
    adv->flag |= CS_PARAM_FLAG_CELL;
    break;
  case ADVKEY_VERTEX_FIELD:
    adv->flag |= CS_PARAM_FLAG_VERTEX;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key %s is not implemented yet."), keyname);

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of a cs_adv_field_t structure
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       keyval    accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_value(cs_adv_field_t    *adv,
                                const char        *val)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  adv->def_type = CS_PARAM_DEF_BY_VALUE;
  adv->flag |= CS_PARAM_FLAG_UNIFORM;

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
 * \param[in]       support   flag to know where is defined the values
 * \param[in]       array     pointer to an array
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_array(cs_adv_field_t     *adv,
                                cs_flag_t           support,
                                const cs_real_t    *array)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  adv->def_type = CS_PARAM_DEF_BY_ARRAY;
  adv->array_flag = support;
  adv->array = array;
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

  _Bool has_previous = (adv->flag & CS_PARAM_FLAG_UNSTEADY) ? true : false;
  int  field_mask = CS_FIELD_PROPERTY;

  if (adv->flag & CS_PARAM_FLAG_VERTEX) { // Add a field attached to vertices

    /* Define the name of the field */
    len = strlen(adv->name) + strlen("_vertices") + 1;
    BFT_MALLOC(field_name, len, char);
    sprintf(field_name, "%s_vertices", adv->name);

    cs_field_t  *fld = cs_field_create(field_name,
                                       field_mask,
                                       CS_MESH_LOCATION_VERTICES,
                                       3,    // always a vector-valued field
                                       true, // interleave
                                       has_previous);

    adv->vtx_field_id = cs_field_id_by_name(field_name);

    /* Allocate and initialize values */
    cs_field_allocate_values(fld);

    BFT_FREE(field_name);

  } // Add a field attached to vertices

  if (adv->flag & CS_PARAM_FLAG_CELL) { // Add a field attached to cells

    /* Define the name of the field */
    len = strlen(adv->name) + strlen("_cells") + 1;
    BFT_MALLOC(field_name, len, char);
    sprintf(field_name, "%s_cells", adv->name);

    cs_field_t  *fld = cs_field_create(field_name,
                                       field_mask,
                                       CS_MESH_LOCATION_CELLS,
                                       3,    // always a vector-valued field
                                       true, // interleave
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
cs_advection_field_get_cell_vector(cs_lnum_t              c_id,
                                   const cs_adv_field_t  *adv,
                                   cs_nvec3_t            *vect)
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

      const cs_real_t  *xc = adv->cdoq->cell_centers + 3*c_id;
      const double  t_cur = adv->time_step->t_cur;

      /* Call the analytic function. result is stored in get */
      adv->def.analytic(t_cur, xc, &get);

      /* Switch to a cs_nvec3_t representation */
      cs_nvec3(get.vect, vect);
    }
    break;

  case CS_PARAM_DEF_BY_ARRAY:
    {
      cs_real_3_t  recoval;

      /* Test if flag has at least the pattern of the reference support */
      if ((adv->array_flag & cs_var_support_dfbyc) == cs_var_support_dfbyc)
        cs_reco_dfbyc_at_cell_center(c_id, adv->connect->c2e, adv->cdoq,
                                     adv->array, recoval);
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid support for evaluating the advection field %s"
                  " at the cell center of cell %d.", adv->name, c_id);

      /* Switch to a cs_nvec3_t representation */
      cs_nvec3(recoval, vect);
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

  cs_lnum_t  c_id;

  const cs_cdo_quantities_t  *quant = adv->cdoq;

  switch (adv->def_type) {

  case CS_PARAM_DEF_BY_VALUE:
    { // Uniform value inside the computational domain
      const cs_get_t  get = adv->def.get;

      for (c_id = 0; c_id < quant->n_cells; c_id++) {

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

      const double  t_cur = adv->time_step->t_cur;

      for (c_id = 0; c_id < quant->n_cells; c_id++) {

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

      /* Test if flag has at least the pattern of the reference support */
      if ((adv->array_flag & cs_var_support_dfbyc) == cs_var_support_dfbyc) {

        for (c_id = 0; c_id < quant->n_cells; c_id++) {

          const cs_lnum_t  shift_c = 3*c_id;

          cs_reco_dfbyc_at_cell_center(c_id, adv->connect->c2e, quant,
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

  const cs_cdo_quantities_t  *quant = adv->cdoq;

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

      const double  t_cur = adv->time_step->t_cur;

      for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {

        const cs_lnum_t  shift = 3*v_id;
        const cs_real_t  *xv = adv->cdoq->vtx_coord + shift;

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
      cs_lnum_t  j, c_id;
      cs_real_3_t  recoval;

      const cs_cdo_connect_t  *topo = adv->connect;

      /* Test if flag has at least the pattern of the reference support */
      if ((adv->array_flag & cs_var_support_dfbyc) == cs_var_support_dfbyc) {

        double  *dc_vol = NULL;

        BFT_MALLOC(dc_vol, quant->n_vertices, double);
        for (j = 0; j < quant->n_vertices; j++) {
          dc_vol[j] = 0;
          vtx_values[3*j] = vtx_values[3*j+1] = vtx_values[3*j+2] = 0.0;
        }

        for (c_id = 0; c_id < quant->n_cells; c_id++) {

          cs_reco_dfbyc_at_cell_center(c_id, topo->c2e, quant,
                                       adv->array, recoval);

          for (j = topo->c2v->idx[c_id]; j < topo->c2v->idx[c_id+1]; j++) {

            const cs_lnum_t  v_id = topo->c2v->ids[j];
            const cs_lnum_t  shift_v = 3*v_id;
            const double  dcc_vol = quant->dcell_vol[j];

            dc_vol[v_id] += dcc_vol;

            for (int k = 0; k < 3; k++)
              vtx_values[shift_v + k] += dcc_vol * recoval[k];

          } // Loop on cell vertices

        } // Loop on cells

        for (j = 0; j < quant->n_vertices; j++) {

          const double  inv_dcvol = 1/dc_vol[j];
          const cs_lnum_t  shift_v = 3*j;

          for (int k = 0; k < 3; k++)
            vtx_values[shift_v+k] *= inv_dcvol;

        } // Loop on vertices

        /* Free temporary buffer */
        BFT_FREE(dc_vol);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid support for evaluating the advection field %s"
                  " at vertices.", adv->name);

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

  const cs_cdo_quantities_t  *cdoq = adv->cdoq;
  const cs_cdo_connect_t  *connect = adv->connect;
  const cs_connect_index_t  *c2e = connect->c2e;

  if (adv->flag & CS_PARAM_FLAG_UNIFORM ||
      adv->flag & CS_PARAM_FLAG_CELLWISE) {

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

        const double  t_cur = adv->time_step->t_cur;
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
              for (k = 0; k < 3; k++)
                xg[k] = one_third * (xc[k] + qe.center[k] + qf.center[k]);
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

        /* Test if flag has at least the pattern of the reference support */
        if ((adv->array_flag & cs_var_support_dfbyc) == cs_var_support_dfbyc)
          for (je = c2e->idx[c_id], _i = 0; je < c2e->idx[c_id+1]; je++, _i++)
            fluxes[_i] = adv->array[je];

        else
          bft_error(__FILE__, __LINE__, 0,
                    " Invalid support for evaluating the advection field %s"
                    " at the cell center of cell %d.", adv->name, c_id);

      } // def. by array
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

  const cs_cdo_quantities_t  *cdoq = adv->cdoq;
  const cs_quant_t  pfq = cdoq->face[f_id];
  const cs_quant_t  peq = cdoq->edge[e_id];
  const cs_real_t  *xv = cdoq->vtx_coord + 3*v_id;

  cs_real_t  surf = cs_surftri(xv, peq.center, pfq.center);

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

      const double  t_cur = adv->time_step->t_cur;

      switch (a_info.quad_type) {

      case CS_QUADRATURE_BARY:
        for (k = 0; k < 3; k++)
          xg[k] = one_third * (xv[k] + peq.center[k] + pfq.center[k]);

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
        if ((adv->array_flag & cs_var_support_dfbyc) == cs_var_support_dfbyc) {

          const cs_connect_index_t  *c2e = adv->connect->c2e;
          const cs_sla_matrix_t  *f2c = adv->connect->f2c;

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
 * \brief  Perform the postprocessing if needed
 *
 * \param[in]      adv     pointer to a cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_post(const cs_adv_field_t  *adv)
{
  cs_lnum_t  v_id, c_id;
  cs_nvec3_t  advect;

  cs_lnum_t  unitv_size = 0;
  cs_real_t  *unitv = NULL;

  const cs_cdo_quantities_t  *cdoq = adv->cdoq;
  const cs_time_step_t  *time_step = adv->time_step;
  const int  nt_cur = time_step->nt_cur;

  if (adv->post_freq < 0)
    return;

  bool  do_post = true;

  if (nt_cur == 0) {
    if (adv->post_freq > 0)
      do_post = false;
  }
  else { /* nt_cur > 0 */
    if (adv->post_freq == 0)
      do_post = false;
    else if (nt_cur % adv->post_freq > 0)
      do_post = false;
  }

  if (do_post) {

    bft_printf(" <post/advection_field> %s\n", adv->name);

    if (adv->post_unitv) {  /* Compute buf_size */

      cs_lnum_t  n_cells = 0, n_vertices = 0;

      if (adv->cell_field_id > -1)
        n_cells = cdoq->n_cells;
      if (adv->vtx_field_id > -1)
        n_vertices = cdoq->n_vertices;

      unitv_size = 3*CS_MAX(n_cells, n_vertices);
      BFT_MALLOC(unitv, unitv_size, cs_real_t);

    }

    if (adv->vtx_field_id > -1) {

      cs_field_t  *fld = cs_field_by_id(adv->vtx_field_id);

      cs_post_write_vertex_var(-1,              // id du maillage de post
                               fld->name,
                               3,               // dim
                               true,            // interlace
                               true,            // true = original mesh
                               CS_POST_TYPE_cs_real_t,
                               fld->val,        // values on vertices
                               time_step);      // time step structure


      if (adv->post_unitv) {

        /* Evaluate the value of the advection field at each vertex */
        for (v_id = 0; v_id < cdoq->n_vertices; v_id++) {

          const cs_lnum_t  shift_v = 3*v_id;

          cs_nvec3(fld->val + shift_v, &advect);
          for (int k = 0; k < 3; k++)
            unitv[shift_v+k] = advect.unitv[k];

        } // Loop on vertices

        char  *label = NULL;
        int  len = strlen(fld->name) + 1 + 10;

        BFT_MALLOC(label, len, char);
        sprintf(label, "%s.Unit", fld->name);

        cs_post_write_vertex_var(-1,              // id du maillage de post
                                 label,
                                 3,               // dim
                                 true,            // interlace
                                 true,            // true = original mesh
                                 CS_POST_TYPE_cs_real_t,
                                 unitv,           // values on vertices
                                 time_step);      // time step structure

        BFT_FREE(label);

      } /* Postprocessing of the unit vector */

    } /* Postprocessing at vertices */

    if (adv->cell_field_id > -1) {

      cs_field_t  *fld = cs_field_by_id(adv->cell_field_id);

      cs_post_write_var(-1,              // id du maillage de post
                        fld->name,
                        3,               // dim
                        true,            // interlace
                        true,            // true = original mesh
                        CS_POST_TYPE_cs_real_t,
                        fld->val,        // values in cells
                        NULL, NULL,      // values at internal/border faces
                        time_step);      // time step structure

      if (adv->post_unitv) {

        /* Evaluate the value of the advection field at each vertex */
        for (c_id = 0; c_id < cdoq->n_cells; c_id++) {

          const cs_lnum_t  shift_c = 3*c_id;

          cs_nvec3(fld->val + shift_c, &advect);
          for (int k = 0; k < 3; k++)
            unitv[shift_c+k] = advect.unitv[k];

        } // Loop on cells

        char  *label = NULL;
        int  len = strlen(fld->name) + 1 + 10;

        BFT_MALLOC(label, len, char);
        sprintf(label, "%s.Unit", fld->name);

        cs_post_write_var(-1,            // id du maillage de post
                          label,
                          3,             // dim
                          true,          // interlace
                          true,          // true = original mesh
                          CS_POST_TYPE_cs_real_t,
                          unitv,         // values in cells
                          NULL, NULL,    // values at internal/border faces
                          time_step);    // time step structure

        BFT_FREE(label);

      } /* Postprocessing of the unit vector */

    } /* Postprocessing at cells */

      /* Free temporay buffers */
    if (adv->post_unitv && unitv_size > 0)
      BFT_FREE(unitv);

  } // Do post

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
