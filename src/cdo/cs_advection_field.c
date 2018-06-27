/*============================================================================
 * Manage the definition/setting of advection fields
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
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_evaluate.h"
#include "cs_field.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_param_cdo.h"
#include "cs_reco.h"
#include "cs_volume_zone.h"
#include "cs_xdef.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

/* Redefined names of function from cs_math to get shorter names */
#define _dp3 cs_math_3_dot_product

#define CS_ADVECTION_FIELD_DBG  0

/*============================================================================
 * Private variables
 *============================================================================*/

static const char _err_empty_adv[] =
  " Stop setting an empty cs_adv_field_t structure.\n"
  " Please check your settings.\n";

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_cdo_quant;
static const cs_cdo_connect_t  *cs_cdo_connect;

  /* Advection fields attached to the computational domain */
static int  _n_adv_fields = 0;
static cs_adv_field_t  **_adv_fields = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the contribution of the flux
 *
 * \param[in]      cdoq       pointer to a cs_cdo_quantities_t structure
 * \param[in]      f2e        face --> edge connectivity
 * \param[in]      e2v        edge --> vertices connectivity
 * \param[in]      bf_id      boundary face id
 * \param[in]      face_flux  value of the normal flux across the face
 * \param[in, out] fluxes     array of values to update
 */
/*----------------------------------------------------------------------------*/

static void
_fill_uniform_boundary_flux(const cs_cdo_quantities_t  *cdoq,
                            const cs_adjacency_t       *f2e,
                            const cs_adjacency_t       *e2v,
                            cs_lnum_t                   bf_id,
                            cs_real_t                   face_flux,
                            cs_real_t                  *fluxes)
{
  const cs_real_t  face_coef = 0.5*face_flux/cdoq->b_face_surf[bf_id];
  const cs_lnum_t  f_id = cdoq->n_i_faces + bf_id;
  const cs_real_t  *xf = cdoq->b_face_center + 3*bf_id;

  for (cs_lnum_t i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

    const cs_lnum_t  eshift = 2*f2e->ids[i];
    const cs_lnum_t  v0 = e2v->ids[eshift];
    const cs_lnum_t  v1 = e2v->ids[eshift+1];
    const double  tef = cs_math_surftri(cdoq->vtx_coord + 3*v0,
                                        cdoq->vtx_coord + 3*v1,
                                        xf);

    const double  weighted_flux = tef * face_coef;

    fluxes[v0] += weighted_flux;
    fluxes[v1] += weighted_flux;

  } /* Loop on face edges */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the contribution of the flux for each vertex
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      f          face id in the cellwise numbering
 * \param[in]      face_flux  value of the normal flux across the face
 * \param[in, out] fluxes     normal boundary flux for each vertex of the face
 */
/*----------------------------------------------------------------------------*/

static void
_fill_cw_uniform_boundary_flux(const cs_cell_mesh_t   *cm,
                               short int               f,
                               cs_real_t               face_flux,
                               cs_real_t              *fluxes)
{
  const cs_real_t  face_coef = 0.5*face_flux/cm->face[f].meas;

  for (short int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

    const short int  eshift = 2*cm->f2e_ids[i];
    const short int  v0 = cm->e2v_ids[eshift];
    const short int  v1 = cm->e2v_ids[eshift+1];
    const double  weighted_flux = face_coef * cm->tef[i];

    fluxes[v0] += weighted_flux;
    fluxes[v1] += weighted_flux;

  } /* Loop on face edges */
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
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_set_shared_pointers(const cs_cdo_quantities_t  *quant,
                                       const cs_cdo_connect_t     *connect)
{
  /* Assign static const pointers */
  cs_cdo_quant = quant;
  cs_cdo_connect = connect;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the number of allocated cs_adv_field_t structures
 *
 * \return the number of advection fields
 */
/*----------------------------------------------------------------------------*/

int
cs_advection_field_get_n_fields(void)
{
  return _n_adv_fields;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search in the array of advection field structures which one has
 *         the name given in argument
 *
 * \param[in]  name        name of the advection field
 *
 * \return a pointer to a cs_adv_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_by_name(const char   *name)
{
  if (_n_adv_fields <= 0)
    return NULL;

  for (int i = 0; i < _n_adv_fields; i++) {

    cs_adv_field_t  *adv = _adv_fields[i];
    if (cs_advection_field_check_name(adv, name))
      return adv;

  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Search in the array of advection field structures which one has
 *         the id given in argument
 *
 * \param[in]  id        identification number
 *
 * \return a pointer to a cs_adv_field_t structure or NULL if not found
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_by_id(int      id)
{
  if (_n_adv_fields <= 0)
    return NULL;
  if (id >= _n_adv_fields || id < 0)
    return NULL;
  if (_adv_fields == NULL)
    return NULL;

  return _adv_fields[id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add and initialize a new user-defined advection field structure
 *
 * \param[in]  name        name of the advection field
 *
 * \return a pointer to the new allocated cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_add_user(const char  *name)
{
  return cs_advection_field_add(name, CS_ADVECTION_FIELD_USER);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add and initialize a new advection field structure
 *
 * \param[in]  name        name of the advection field
 * \param[in]  type        type of advection field
 *
 * \return a pointer to the new allocated cs_adv_field_t structure
 */
/*----------------------------------------------------------------------------*/

cs_adv_field_t *
cs_advection_field_add(const char                  *name,
                       cs_advection_field_type_t    type)
{
  if (name == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " A non-empty name is mandatory to add a new advection field");

  cs_adv_field_t  *adv = cs_advection_field_by_name(name);
  if (adv != NULL) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" An existing advection field has already the name %s.\n"
                    " Stop adding this advection field.\n"), name);
    return adv;
  }

  int  new_id = _n_adv_fields;
  _n_adv_fields++;
  BFT_REALLOC(_adv_fields, _n_adv_fields, cs_adv_field_t *);
  _adv_fields[new_id] = NULL;

  BFT_MALLOC(adv, 1, cs_adv_field_t);

  adv->id = new_id;
  adv->type = type;

  /* Copy name */
  int  len = strlen(name) + 1;
  BFT_MALLOC(adv->name, len, char);
  strncpy(adv->name, name, len);

  /* Default initialization */
  adv->flag = 0;
  adv->vtx_field_id = -1;
  adv->cell_field_id = -1;
  adv->bdy_field_id = -1;

  adv->definition = NULL;
  adv->n_bdy_flux_defs = 0;
  adv->bdy_flux_defs = NULL;
  adv->bdy_def_ids = NULL;

  /* Store the new advection field */
  _adv_fields[new_id] = adv;

  return adv;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free all alllocated cs_adv_field_t structures and its related array
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_destroy_all(void)
{
  if (_adv_fields == NULL)
    return;

  for (int i = 0; i < _n_adv_fields; i++) {

    cs_adv_field_t  *adv = _adv_fields[i];

    adv->definition = cs_xdef_free(adv->definition);

    for (int id = 0; id < adv->n_bdy_flux_defs; id++)
      adv->bdy_flux_defs[id] = cs_xdef_free(adv->bdy_flux_defs[id]);

    if (adv->n_bdy_flux_defs > 0) BFT_FREE(adv->bdy_flux_defs);
    if (adv->bdy_def_ids != NULL)   BFT_FREE(adv->bdy_def_ids);

    BFT_FREE(adv->name);
    BFT_FREE(adv);

    /* All other pointers are shared */

  }  /* Loop on advection fields */

  BFT_FREE(_adv_fields);
  _n_adv_fields = 0;
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

  if (adv->definition->state & CS_FLAG_STATE_UNIFORM)
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

  cs_flag_t  state = adv->definition->state;

  if (state & CS_FLAG_STATE_UNIFORM)
    return true;
  if (state & CS_FLAG_STATE_CELLWISE)
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

cs_xdef_type_t
cs_advection_field_get_deftype(const cs_adv_field_t   *adv)
{
  if (adv == NULL)
    return CS_N_XDEF_TYPES;

  return cs_xdef_get_type(adv->definition);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Print all setup information related to cs_adv_field_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_log_setup(void)
{
  if (_adv_fields == NULL)
    return;

  cs_log_printf(CS_LOG_SETUP, "\n%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, "\tSummary of the advection field\n");
  cs_log_printf(CS_LOG_SETUP, "%s", lsepline);
  cs_log_printf(CS_LOG_SETUP, " -msg- n_advection_fields       %d\n",
                _n_adv_fields);

  for (int i = 0; i < _n_adv_fields; i++) {

    const cs_adv_field_t  *adv = _adv_fields[i];

    cs_log_printf(CS_LOG_SETUP, " <AdvectionField/%s> id: %d\n",
                  adv->name, adv->id);

    switch (adv->type) {
    case CS_ADVECTION_FIELD_NAVSTO:
      cs_log_printf(CS_LOG_SETUP,
                    " <AdvectionField/%s> Related to Navier-Stokes\n",
                    adv->name);
      break;

    case CS_ADVECTION_FIELD_GWF:
      cs_log_printf(CS_LOG_SETUP,
                    " <AdvectionField/%s> Related to the \"Groundwater Flow\""
                    " module\n", adv->name);
      break;

    case CS_ADVECTION_FIELD_USER:
      cs_log_printf(CS_LOG_SETUP,
                    " <AdvectionField/%s> User-defined\n", adv->name);
      break;

    default:
      break;
    }

    if (adv->cell_field_id > -1)
      cs_log_printf(CS_LOG_SETUP, " <AdvectionField/%s> Defined at cells\n",
                    adv->name);
    if (adv->vtx_field_id > -1)
      cs_log_printf(CS_LOG_SETUP, " <AdvectionField/%s> Defined at vertices\n",
                    adv->name);
    if (adv->flag & CS_ADVECTION_FIELD_POST_COURANT)
      cs_log_printf(CS_LOG_SETUP,
                    " <AdvectionField/%s> Postprocess the Courant number\n",
                    adv->name);
    if (adv->flag & CS_ADVECTION_FIELD_STEADY)
      cs_log_printf(CS_LOG_SETUP,
                    " <AdvectionField/%s> Steady-state\n", adv->name);

    cs_xdef_log(adv->definition);

  }  /* Loop on advection fields */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set optional parameters related to a cs_adv_field_t structure
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       key       key related to the member of adv to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_set_option(cs_adv_field_t            *adv,
                              cs_advection_field_key_t   key)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  switch(key) {

  case CS_ADVKEY_DEFINE_AT_VERTICES:
    adv->vtx_field_id = -2;     /* To be created */
    break;
  case CS_ADVKEY_POST_COURANT:
    adv->flag |= CS_ADVECTION_FIELD_POST_COURANT;
    break;
  case CS_ADVKEY_STATE_STEADY:
    adv->flag |= CS_ADVECTION_FIELD_STEADY;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" Key not implemented for setting an advection field."));
    break;

  } /* Switch on keys */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of a cs_adv_field_t structure
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       vector    values to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_value(cs_adv_field_t    *adv,
                                cs_real_t          vector[3])
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_CELLWISE;
  cs_flag_t  meta_flag = CS_FLAG_FULL_LOC;

  adv->definition = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                          3,  /* dim. */
                                          0,  /* zone_id = 0 => all cells */
                                          state_flag,
                                          meta_flag,
                                          vector);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an analytic function
 *
 * \param[in, out]  adv     pointer to a cs_adv_field_t structure
 * \param[in]       func    pointer to a function
 * \param[in]       input   NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_analytic(cs_adv_field_t        *adv,
                                   cs_analytic_func_t    *func,
                                   void                  *input)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = CS_FLAG_FULL_LOC;
  cs_xdef_analytic_input_t  anai = {.func = func,
                                    .input = input };

  adv->definition = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          3,  /* dim. */
                                          0,  /* zone_id = 0 => all cells */
                                          state_flag,
                                          meta_flag,
                                          &anai);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an array of values
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       index     optional pointer to the array index
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_array(cs_adv_field_t    *adv,
                                cs_flag_t          loc,
                                cs_real_t         *array,
                                cs_lnum_t         *index)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = 0; /* Will be updated during the creation */
  cs_flag_t  meta_flag = CS_FLAG_FULL_LOC;
  cs_xdef_array_input_t  input = {.stride = 3,
                                  .loc = loc,
                                  .values = array,
                                  .index = index };

  adv->definition = cs_xdef_volume_create(CS_XDEF_BY_ARRAY,
                                          3,  /* dim */
                                          0,  /* zone_id */
                                          state_flag,
                                          meta_flag,
                                          &input);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_adv_field_t structure thanks to an array of values
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       field     pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_by_field(cs_adv_field_t    *adv,
                                cs_field_t        *field)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = 0; /* Will be updated during the creation */
  cs_flag_t  meta_flag = CS_FLAG_FULL_LOC;

  assert(field->dim == 3); /* sanity check since fields are either at vertices
                              or at cells in this case */
  adv->definition = cs_xdef_volume_create(CS_XDEF_BY_FIELD,
                                          3,
                                          0,  /* zone_id */
                                          state_flag,
                                          meta_flag,
                                          field);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of the boundary normal flux for the given
 *         cs_adv_field_t structure
 *
 * \param[in, out]  adv           pointer to a cs_adv_field_t structure
 * \param[in]       zname         name of the boundary zone to consider
 * \param[in]       normal_flux   value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_boundary_flux_by_value(cs_adv_field_t    *adv,
                                              const char        *zname,
                                              cs_real_t          normal_flux)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = CS_FLAG_STATE_UNIFORM | CS_FLAG_STATE_FACEWISE;
  cs_flag_t  meta_flag = 0;

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          1,  /* dim. */
                                          cs_get_bdy_zone_id(zname),
                                          state_flag,
                                          meta_flag,
                                          (void *)&normal_flux);

  int  def_id = adv->n_bdy_flux_defs;
  adv->n_bdy_flux_defs += 1;
  BFT_REALLOC(adv->bdy_flux_defs, adv->n_bdy_flux_defs, cs_xdef_t *);
  adv->bdy_flux_defs[def_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of the boundary normal flux for the given
 *         cs_adv_field_t structure using an analytic function
 *
 * \param[in, out]  adv     pointer to a cs_adv_field_t structure
 * \param[in]       zname   name of the boundary zone to consider
 * \param[in]       func    pointer to a function
 * \param[in]       input   NULL or pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_boundary_flux_by_analytic(cs_adv_field_t        *adv,
                                                 const char            *zname,
                                                 cs_analytic_func_t    *func,
                                                 void                  *input)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag = 0;
  cs_flag_t  meta_flag = 0;
  cs_xdef_analytic_input_t  anai = {.func = func,
                                    .input = input };

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          1,  /* dim. */
                                          cs_get_bdy_zone_id(zname),
                                          state_flag,
                                          meta_flag,
                                          &anai);

  int  def_id = adv->n_bdy_flux_defs;
  adv->n_bdy_flux_defs += 1;
  BFT_REALLOC(adv->bdy_flux_defs, adv->n_bdy_flux_defs, cs_xdef_t *);
  adv->bdy_flux_defs[def_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the value of the boundary normal flux for the given
 *         cs_adv_field_t structure using an array of values
 *
 * \param[in, out]  adv       pointer to a cs_adv_field_t structure
 * \param[in]       zname     name of the boundary zone to consider
 * \param[in]       loc       information to know where are located values
 * \param[in]       array     pointer to an array
 * \param[in]       index     optional pointer to the array index
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_def_boundary_flux_by_array(cs_adv_field_t    *adv,
                                              const char        *zname,
                                              cs_flag_t          loc,
                                              cs_real_t         *array,
                                              cs_lnum_t         *index)
{
  if (adv == NULL)
    bft_error(__FILE__, __LINE__, 0, _(_err_empty_adv));

  cs_flag_t  state_flag =  0;
  cs_flag_t  meta_flag = 0;
  cs_xdef_array_input_t  input = {.stride = 1,
                                  .loc = loc,
                                  .values = array,
                                  .index = index };

  int  z_id = cs_get_bdy_zone_id(zname);
  if (z_id == 0)
    meta_flag  |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ARRAY,
                                          1,  /* dim. */
                                          z_id,
                                          state_flag,
                                          meta_flag,
                                          &input);

  int  def_id = adv->n_bdy_flux_defs;
  adv->n_bdy_flux_defs += 1;
  BFT_REALLOC(adv->bdy_flux_defs, adv->n_bdy_flux_defs, cs_xdef_t *);
  adv->bdy_flux_defs[def_id] = d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create all needed cs_field_t structures related to an advection
 *         field
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_create_fields(void)
{
  int  len;
  char *field_name = NULL;

  for (int i = 0; i < _n_adv_fields; i++) {

    cs_adv_field_t  *adv = _adv_fields[i];
    assert(adv != NULL);

    bool  has_previous = (adv->flag & CS_ADVECTION_FIELD_STEADY) ? true:false;
    int  field_mask = CS_FIELD_PROPERTY;

    {  /* (Normal flux) Field at boundary faces:
          Always create a field at the boundary faces for taking into account
          the normal flux used in the treatment of the boundary conditions */

      /* Define the name of the field */
      len = strlen(adv->name) + strlen("_boundary_flux") + 1;
      BFT_MALLOC(field_name, len, char);
      sprintf(field_name, "%s_boundary_flux", adv->name);

      cs_field_t  *fld = cs_field_create(field_name,
                                         field_mask,
                                         CS_MESH_LOCATION_BOUNDARY_FACES,
                                         1,   /* always a scalar-valued field */
                                         has_previous);

      cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
      cs_field_set_key_int(fld, cs_field_key_id("post_vis"), 1);

      adv->bdy_field_id = cs_field_id_by_name(field_name);

      BFT_FREE(field_name);

    } /* Add a field attached to boundary faces */

    { /* Add a field attached to cells (Always created since it's used to
         define the numerical scheme for advection */

      /* Define the name of the field */
      len = strlen(adv->name) + strlen("_cells") + 1;
      BFT_MALLOC(field_name, len, char);
      sprintf(field_name, "%s_cells", adv->name);

      cs_field_t  *fld = cs_field_create(field_name,
                                         field_mask,
                                         CS_MESH_LOCATION_CELLS,
                                         3,  /* always a vector-valued field */
                                         has_previous);

      cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
      cs_field_set_key_int(fld, cs_field_key_id("post_vis"), 1);

      adv->cell_field_id = cs_field_id_by_name(field_name);

      BFT_FREE(field_name);

    } /* Add a field attached to cells */

    if (adv->vtx_field_id == -2) { /* Add a field attached to vertices */

      /* Define the name of the field */
      len = strlen(adv->name) + strlen("_vertices") + 1;
      BFT_MALLOC(field_name, len, char);
      sprintf(field_name, "%s_vertices", adv->name);

      cs_field_t  *fld = cs_field_create(field_name,
                                         field_mask,
                                         CS_MESH_LOCATION_VERTICES,
                                         3,  /* always a vector-valued field */
                                         has_previous);

      cs_field_set_key_int(fld, cs_field_key_id("log"), 1);
      cs_field_set_key_int(fld, cs_field_key_id("post_vis"), 1);

      adv->vtx_field_id = cs_field_id_by_name(field_name);

      BFT_FREE(field_name);

    } /* Add a field attached to vertices */

  } /* Loop on advection fields */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Last stage of the definition of an advection field based on several
 *         definitions (i.e. definition by subdomains on the boundary)
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_finalize_setup(void)
{
  if (_n_adv_fields == 0)
    return;

  for (int i = 0; i < _n_adv_fields; i++) {

    cs_adv_field_t  *adv = _adv_fields[i];
    assert(adv != NULL);

    if (adv->n_bdy_flux_defs > 1) {

      const cs_lnum_t  n_b_faces = cs_cdo_quant->n_b_faces;

      BFT_MALLOC(adv->bdy_def_ids, n_b_faces, short int);
#     pragma omp parallel for if (n_b_faces > CS_THR_MIN)
      for (cs_lnum_t j = 0; j < n_b_faces; j++)
        adv->bdy_def_ids[j] = -1;  /* Unset by default */

      for (short int def_id = 0; def_id < adv->n_bdy_flux_defs; def_id++) {

        const cs_xdef_t  *def = adv->bdy_flux_defs[def_id];
        assert(def->z_id > 0);
        assert(def->support == CS_XDEF_SUPPORT_BOUNDARY);

        const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);
        assert(z != NULL);

#       pragma omp parallel for if (z->n_elts > CS_THR_MIN)
        for (cs_lnum_t j = 0; j < z->n_elts; j++)
          adv->bdy_def_ids[z->elt_ids[j]] = def_id;

      }

    } /* More than one definition */

  } /* Loop on advection fields */

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

  cs_field_t  *f = cs_advection_field_get_field(adv, CS_MESH_LOCATION_CELLS);

  cs_nvec3(f->val + 3*c_id, vect);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the mean-value of the advection field inside each cell
 *
 * \param[in]      adv           pointer to a cs_adv_field_t structure
 * \param[in]      time_eval     physical time at which one evaluates the term
 * \param[in, out] cell_values   array of values at cell centers
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_in_cells(const cs_adv_field_t  *adv,
                            cs_real_t              time_eval,
                            cs_real_t             *cell_values)
{
  if (adv == NULL)
    return;

  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;

  cs_xdef_t  *def = adv->definition;

  /* Sanity checks */
  assert(def != NULL);

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    {
      const cs_real_t  *constant_val = (cs_real_t *)def->input;

#     pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_cells; i++) {
        cell_values[3*i  ] = constant_val[0];
        cell_values[3*i+1] = constant_val[1];
        cell_values[3*i+2] = constant_val[2];
      }

    }
    break; /* definition by value */

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *array_input = (cs_xdef_array_input_t *)def->input;

      const int  stride = array_input->stride;

      if (cs_flag_test(array_input->loc, cs_flag_primal_cell)) {

        assert(stride == 3);
        memcpy(cell_values, array_input->values,
               stride*cdoq->n_cells * sizeof(cs_real_t));

      }
      else if (cs_flag_test(array_input->loc, cs_flag_dual_face_byc)) {

        assert(array_input->index == cs_cdo_connect->c2e->idx);

#       pragma omp parallel for if (cdoq->n_cells > CS_THR_MIN)
        for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++)
          cs_reco_dfbyc_at_cell_center(c_id,
                                       cs_cdo_connect->c2e,
                                       cdoq,
                                       array_input->values,
                                       cell_values + 3*c_id);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid location for array", __func__);

    }
    break; /* definition by array */

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *field = (cs_field_t *)def->input;
      assert(field != NULL);

      if (field->location_id == cs_mesh_location_get_id_by_name("cells")) {

        assert(field->dim == 3);

        /* If not the field associated to the advection field */
        if (field->id != adv->cell_field_id)
          memcpy(cell_values, field->val, 3*cdoq->n_cells*sizeof(cs_real_t));

      }
      else if (field->location_id ==
               cs_mesh_location_get_id_by_name("vertices")) {

        assert(field->dim == 3);

        cs_reco_vect_pv_at_cell_centers(cs_cdo_connect->c2v,
                                        cdoq,
                                        field->val,
                                        cell_values);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid case for the input field", __func__);

    }
    break; /* definition by field */

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_evaluate_average_on_cells_by_analytic(def, time_eval, cell_values);
    break; /* definition by analytic */

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Incompatible type of definition.",
              __func__);
    break;

  } /* Type of definition */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at vertices
 *
 * \param[in]      adv          pointer to a cs_adv_field_t structure
 * \param[in]      time_eval    physical time at which one evaluates the term
 * \param[in, out] vtx_values   array storing the results
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_at_vertices(const cs_adv_field_t  *adv,
                               cs_real_t              time_eval,
                               cs_real_t             *vtx_values)
{
  if (adv == NULL)
    return;

  cs_xdef_t  *def = adv->definition;

  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    {
      const cs_real_t  *constant_val = (cs_real_t *)def->input;

#     pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
      for (cs_lnum_t i = 0; i < cdoq->n_vertices; i++) {
        vtx_values[3*i  ] = constant_val[0];
        vtx_values[3*i+1] = constant_val[1];
        vtx_values[3*i+2] = constant_val[2];
      }

    }
    break; /* definition by value */

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *array_input = (cs_xdef_array_input_t *)def->input;

      const int  stride = array_input->stride;

      if (cs_flag_test(array_input->loc, cs_flag_primal_vtx)) {

        assert(stride == 3);
        memcpy(vtx_values, array_input->values,
               3*cdoq->n_vertices * sizeof(cs_real_t));

      }
      else if (cs_flag_test(array_input->loc, cs_flag_primal_cell)) {

        assert(stride == 3);
        cs_reco_vect_pv_from_pc(cs_cdo_connect->c2v,
                                cdoq,
                                array_input->values,
                                vtx_values);

      }
      else if (cs_flag_test(array_input->loc, cs_flag_dual_face_byc)) {

        assert(array_input->index == cs_cdo_connect->c2e->idx);
        memset(vtx_values, 0, 3*cdoq->n_vertices*sizeof(cs_real_t));

        const cs_adjacency_t  *c2v = cs_cdo_connect->c2v;

        for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

          cs_real_t  cell_vector[3];
          cs_reco_dfbyc_at_cell_center(c_id,
                                       cs_cdo_connect->c2e,
                                       cdoq,
                                       array_input->values,
                                       cell_vector);

          for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

            const cs_real_t  vc_vol = cdoq->dcell_vol[j];
            cs_real_t  *_val = vtx_values + 3*c2v->ids[j];

            _val[0] += vc_vol * cell_vector[0];
            _val[1] += vc_vol * cell_vector[1];
            _val[2] += vc_vol * cell_vector[2];

          }

        } /* Loop on cells */

        cs_real_t  *dual_vol = NULL;
        BFT_MALLOC(dual_vol, cdoq->n_vertices, cs_real_t);
        cs_cdo_quantities_compute_dual_volumes(cdoq, c2v, dual_vol);

#       pragma omp parallel for if (cdoq->n_vertices > CS_THR_MIN)
        for (cs_lnum_t v_id = 0; v_id < cdoq->n_vertices; v_id++) {
          const cs_real_t  invvol = 1./dual_vol[v_id];
          for (int k = 0; k < 3; k++)
            vtx_values[3*v_id+k] *= invvol;
        }

        BFT_FREE(dual_vol);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid location for array", __func__);

    }
    break; /* definition by array */

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *field = (cs_field_t *)def->input;
      assert(field != NULL);

      if (field->location_id == cs_mesh_location_get_id_by_name("cells")) {

        assert(field->dim == 3);
        cs_reco_vect_pv_from_pc(cs_cdo_connect->c2v,
                                cdoq,
                                field->val,
                                vtx_values);

      }
      else if (field->location_id ==
               cs_mesh_location_get_id_by_name("vertices")) {

        assert(field->dim == 3);

        /* If not the field associated to the advection field */
        if (field->id != adv->vtx_field_id)
          memcpy(vtx_values, field->val, 3*cdoq->n_vertices*sizeof(cs_real_t));

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid case for the input field", __func__);

    }
    break; /* definition by field */

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      cs_flag_t  dof_flag = cs_flag_primal_vtx | CS_FLAG_VECTOR;

      cs_evaluate_potential_by_analytic(dof_flag, def, time_eval, vtx_values);
    }
    break; /* definition by analytic */

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Incompatible type of definition.",
              __func__);
    break;

  } /* Type of definition */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the advection field at a specific location
 *         inside a cell
 *
 * \param[in]      adv          pointer to a cs_adv_field_t structure
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      xyz          location where to perform the evaluation
 * \param[in]      time_eval    physical time at which one evaluates the term
 * \param[in, out] eval         pointer to a cs_nvec3_t
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_eval_at_xyz(const cs_adv_field_t  *adv,
                               const cs_cell_mesh_t  *cm,
                               const cs_real_3_t      xyz,
                               cs_real_t              time_eval,
                               cs_nvec3_t            *eval)
{
  if (adv == NULL)
    return;

  cs_xdef_t  *def = adv->definition;
  cs_real_3_t  vector_values = {0, 0, 0};

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    {
      const cs_real_t  *constant_val = (cs_real_t *)def->input;

      cs_nvec3(constant_val, eval);
    }
    break; /* definition by value */

  case CS_XDEF_BY_ARRAY:
    cs_xdef_cw_eval_vector_at_xyz_by_array(cm, 1, xyz, time_eval, def->input,
                                           vector_values);
    cs_nvec3(vector_values, eval);
    break;

  case CS_XDEF_BY_FIELD:
    cs_xdef_cw_eval_vector_at_xyz_by_field(cm, 1, xyz, time_eval, def->input,
                                           vector_values);
    cs_nvec3(vector_values, eval);
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    cs_xdef_cw_eval_at_xyz_by_analytic(cm, 1, xyz, time_eval, def->input,
                                       vector_values);
    cs_nvec3(vector_values, eval);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Incompatible type of definition.",
              __func__);
    break;

  } /* Type of definition */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the normal flux of the advection field
 *         across the boundary faces
 *
 * \param[in]      adv          pointer to a cs_adv_field_t structure
 * \param[in]      time_eval    physical time at which one evaluates the term
 * \param[in, out] flx_values   array storing the results
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_across_boundary(const cs_adv_field_t  *adv,
                                   cs_real_t              time_eval,
                                   cs_real_t             *flx_values)
{
  if (adv == NULL)
    return;

  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;
  const cs_lnum_t  n_b_faces = cdoq->n_b_faces;
  const cs_lnum_t  n_i_faces = cdoq->n_i_faces;

  cs_nvec3_t  nvec;

  if (adv->n_bdy_flux_defs == 0) {

    /* No specific definition of the boundary flux */
    cs_xdef_t  *def = adv->definition;

    assert(def->dim == 3);

    switch (def->type) {

    case CS_XDEF_BY_VALUE:
      {
        const cs_real_t  *constant_val = (cs_real_t *)def->input;

        cs_nvec3(constant_val, &nvec);
        for (cs_lnum_t i = 0; i < n_b_faces; i++) {
          const cs_quant_t  qf = cs_quant_set_face(n_i_faces + i, cdoq);
          flx_values[i] = qf.meas * nvec.meas * _dp3(qf.unitv, nvec.unitv);
        }

      }
      break;

    case CS_XDEF_BY_ANALYTIC_FUNCTION:
      {
        const cs_adjacency_t  *f2e = cs_cdo_connect->f2e;
        const cs_adjacency_t  *e2v = cs_cdo_connect->e2v;
        const cs_real_t  *xv = cdoq->vtx_coord;

        cs_quadrature_tria_integral_t
          *compute_integral = cs_quadrature_get_tria_integral(def->dim,
                                                              def->qtype);
        cs_xdef_analytic_input_t *anai = (cs_xdef_analytic_input_t *)def->input;

        for (cs_lnum_t i = 0; i < n_b_faces; i++) {

          const cs_lnum_t  f_id = n_i_faces + i;
          const cs_quant_t  pfq = cs_quant_set_face(f_id, cdoq);
          const cs_lnum_t  start_idx = f2e->idx[f_id];
          const cs_lnum_t  end_idx = f2e->idx[f_id+1];

          cs_real_t  val[3] = {0, 0, 0};

          switch (end_idx - start_idx) {

          case CS_TRIANGLE_CASE: /* Triangle: one-shot computation */
            {
              cs_lnum_t  v1, v2, v3;

              cs_connect_get_next_3_vertices(f2e->ids, e2v->ids, start_idx,
                                             &v1, &v2, &v3);

              compute_integral(time_eval,
                               xv + 3*v1, xv + 3*v2, xv + 3*v3, pfq.meas,
                               anai->func, anai->input, val);

            }
            break;

          default:
            for (cs_lnum_t j = start_idx; j < end_idx; j++) {

              const cs_lnum_t  _2e = 2*f2e->ids[j];
              const cs_real_t  *xv1 = xv + 3*e2v->ids[_2e];
              const cs_real_t  *xv2 = xv + 3*e2v->ids[_2e+1];
              const cs_real_t  tef_meas = cs_math_surftri(xv1, xv2, pfq.center);

              /* val is updated (+=) */
              compute_integral(time_eval, xv1, xv2, pfq.center, tef_meas,
                               anai->func, anai->input, val);

            } /* Loop on face edges */

          } /* End of switch */

          flx_values[i] = _dp3(pfq.unitv, val);

        }

      }
      break; /* definition by analytic */

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Incompatible type of definition.",
                __func__);
      break;

    } /* Type of definition */

  }
  else {

    for (int def_id = 0; def_id < adv->n_bdy_flux_defs; def_id++) {

      const cs_xdef_t  *def = adv->bdy_flux_defs[def_id];
      const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);
      assert(def->support == CS_XDEF_SUPPORT_BOUNDARY);

      switch (def->type) {

      case CS_XDEF_BY_VALUE:
        {
          const cs_real_t  *constant_val = (cs_real_t *)def->input;

          if (z->elt_ids == NULL) {
            assert(z->n_elts == n_b_faces);
#           pragma omp parallel for if (n_b_faces > CS_THR_MIN)
            for (cs_lnum_t i = 0; i < n_b_faces; i++)
              flx_values[i] = constant_val[0];
          }
          else {
#           pragma omp parallel for if (z->n_elts > CS_THR_MIN)
            for (cs_lnum_t i = 0; i < z->n_elts; i++)
              flx_values[z->elt_ids[i]] = constant_val[0];
          }

        }
        break;

      case CS_XDEF_BY_ARRAY:
        {
          const cs_xdef_array_input_t  *input =
            (cs_xdef_array_input_t *)def->input;
          const cs_real_t  *val = input->values;

          assert(input->stride == 1);
          assert(def->meta & CS_FLAG_FULL_LOC || z->elt_ids == NULL);

          if (cs_flag_test(input->loc, cs_flag_primal_face))
            memcpy(flx_values, val, sizeof(cs_real_t)*n_b_faces);

          else if (cs_flag_test(input->loc, cs_flag_dual_closure_byf)) {

            const cs_lnum_t  *idx = input->index;
            for (cs_lnum_t bf_id = 0; bf_id < n_b_faces; bf_id++) {
              flx_values[bf_id] = 0;
              for (cs_lnum_t i = idx[bf_id]; i < idx[bf_id+1]; i++) {
                flx_values[bf_id] += val[i];
              }
            } /* Loop on border faces */

            }
          else
            bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.", __func__);

        }
        break; /* definition by array */

      case CS_XDEF_BY_FIELD:
        {
          cs_field_t  *field = (cs_field_t *)def->input;

          assert(field->dim == 1);
          assert(def->meta & CS_FLAG_FULL_LOC || z->elt_ids == NULL);

          if (field->location_id ==
              cs_mesh_location_get_id_by_name(N_("boundary faces"))) {
            memcpy(flx_values, field->val, sizeof(cs_real_t)*n_b_faces);
          }
          else
            bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.", __func__);

        }
        break; /* definition by field */

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        {
          cs_xdef_analytic_input_t  *anai =
            (cs_xdef_analytic_input_t *)def->input;

          anai->func(time_eval,
                     z->n_elts, z->elt_ids, cdoq->b_face_center,
                     false,   /* compacted output ? */
                     anai->input,
                     flx_values);

        }
        break; /* definition by analytic */

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Incompatible type of definition.", __func__);
        break;

      } /* Type of definition */

    } /*  Loop on boundary definitions */

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         the closure of the dual cell related to each vertex
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      adv        pointer to a cs_adv_field_t structure
 * \param[in]      f          face id in the cellwise numbering
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] fluxes     normal boundary flux for each vertex of the face
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_f2v_boundary_flux(const cs_cell_mesh_t   *cm,
                                         const cs_adv_field_t   *adv,
                                         short int               f,
                                         cs_real_t               time_eval,
                                         cs_real_t              *fluxes)
{
  if (adv == NULL)
    return;

  if (fluxes == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " fluxes array should be allocated before the call.");

  cs_lnum_t  bf_id = cm->f_ids[f] - cs_cdo_quant->n_i_faces;
  assert(bf_id > -1);
  assert(cs_flag_test(cm->flag,  CS_CDO_LOCAL_PFQ| CS_CDO_LOCAL_FEQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FE));

  for (short int v = 0; v < cm->n_vc; v++) fluxes[v] = 0;

  if (adv->n_bdy_flux_defs == 0) {

    /* Implicit definition of the boundary flux from the definition of the
       advection field inside the computational domain */
    cs_xdef_t  *def = adv->definition;

    assert(adv->bdy_def_ids == NULL);
    switch (def->type) {

    case CS_XDEF_BY_ANALYTIC_FUNCTION:
      {
        cs_quadrature_tria_integral_t
          *compute_integral = cs_quadrature_get_tria_integral(def->dim,
                                                              def->qtype);
        cs_xdef_analytic_input_t *anai = (cs_xdef_analytic_input_t *)def->input;

        const cs_quant_t  pfq = cm->face[f];
        const short int  end_idx = cm->f2e_idx[f+1], start_idx = cm->f2e_idx[f];
        const short int n_ef = end_idx - start_idx; /* #vertices (= #edges) */
        const short int *f2e_ids = cm->f2e_ids + start_idx;
        const double  *tef = cm->tef + start_idx;

        for (short int e = 0; e < n_ef; e++) { /* Loop on face edges */

          const short int  _2e = 2*f2e_ids[e];
          const short int  v1 = cm->e2v_ids[_2e];
          const short int  v2 = cm->e2v_ids[_2e+1];

          cs_real_t  val[3] = {0, 0, 0};

          /* val is updated (+=) */
          compute_integral(time_eval,
                           cm->xv + 3*v1, cm->xv + 3*v2, pfq.center, tef[e],
                           anai->func, anai->input, val);

          const double  e_contrib = 0.5*_dp3(pfq.unitv, val);

          fluxes[v1] += e_contrib;
          fluxes[v2] += e_contrib;

        } /* Loop on face edges */

      }
      break;

    default:
      {
        const cs_field_t  *bflx =
          cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

        _fill_cw_uniform_boundary_flux(cm, f, bflx->val[bf_id], fluxes);
      }
      break;

    } /* End of switch */

  }
  else {

    cs_xdef_t  *def = NULL;
    if (adv->bdy_def_ids == NULL) {
      def = adv->bdy_flux_defs[0];
      assert(adv->n_bdy_flux_defs == 1);
    }
    else
      def = adv->bdy_flux_defs[adv->bdy_def_ids[bf_id]];

    const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);

    switch (def->type) {

    case CS_XDEF_BY_VALUE:
      {
        const cs_real_t  *constant_val = (cs_real_t *)def->input;

        _fill_cw_uniform_boundary_flux(cm, f, constant_val[0], fluxes);

      }
      break; /* by value */

    case CS_XDEF_BY_ARRAY:
      {
        const cs_xdef_array_input_t  *input =
          (cs_xdef_array_input_t *)def->input;
        const cs_real_t  *val = input->values;

        assert(input->stride == 1);
        assert(z->id == 0);

        if (cs_flag_test(input->loc, cs_flag_primal_face))
          _fill_cw_uniform_boundary_flux(cm, f, val[bf_id], fluxes);

        else if (cs_flag_test(input->loc, cs_flag_dual_closure_byf)) {

          const cs_adjacency_t  *bf2v = cs_cdo_connect->bf2v;
          const cs_lnum_t  *idx = bf2v->idx;
          assert(idx == input->index);

          for (cs_lnum_t i = idx[bf_id]; i < idx[bf_id+1]; i++) {
            const cs_lnum_t  v_id = bf2v->ids[i];
            short int v = -1;
            for (v = 0; v < cm->n_vc; v++)
              if (cm->v_ids[v] == v_id)
                break;
            assert(v > -1 && v != cm->n_vc);
            fluxes[v] += val[i];
          }

        }
        else
          bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.", __func__);

      }
      break; /* by_array */

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid case", __func__);
      break;

    } /* End of switch on the type of definition */

  }   /* n_bdy_flux_defs > 0 */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         the (primal) faces of a cell
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      adv        pointer to a cs_adv_field_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] fluxes     array of values attached to dual faces of a cell
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_cw_face_flux(const cs_cell_mesh_t       *cm,
                                    const cs_adv_field_t       *adv,
                                    cs_real_t                   time_eval,
                                    cs_real_t                  *fluxes)
{
  if (adv == NULL)
    return;

  if (fluxes == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " fluxes array should be allocated before the call.");

  cs_xdef_t  *def = adv->definition;

  cs_field_t  *bflux =
    cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);
  const cs_real_t  *bflux_val = bflux->val;

  /* Sanity checks */
  assert(def != NULL);
  assert(bflux_val != NULL);

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    {
      /* Loop on cell faces */
      for (short int f = 0; f < cm->n_fc; f++)
        cs_xdef_cw_eval_flux_by_val(cm, f, time_eval, def->input, fluxes);

    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_PFQ));

      /* Loop on cell faces */
      for (short int f = 0; f < cm->n_fc; f++)
        cs_xdef_cw_eval_flux_by_analytic(cm, f,
                                         time_eval,
                                         def->input, def->qtype,
                                         fluxes);

    } /* definition by analytic function */
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *input = (cs_xdef_array_input_t *)def->input;

      /* Test if location has at least the pattern of the reference support */
      if (cs_flag_test(input->loc, cs_flag_primal_face)) {

        /* Sanity check */
        assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PF));
        for (short int f = 0; f < cm->n_fc; f++)
          fluxes[f] = input->values[cm->f_ids[f]];

      }
      else if (cs_flag_test(input->loc, cs_flag_primal_cell)) {

        /* Retrieve the advection field:
           Switch to a cs_nvec3_t representation */
        cs_nvec3_t  adv_vect;
        cs_nvec3(input->values + 3*cm->c_id, &adv_vect);

        /* Sanity check */
        assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ));

        /* Loop on cell faces */
        for (short int f = 0; f < cm->n_fc; f++)
          fluxes[f] = adv_vect.meas*cm->face[f].meas * _dp3(adv_vect.unitv,
                                                            cm->face[f].unitv);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid support for evaluating the advection field %s"
                  " at the cell center of cell %d.", adv->name, cm->c_id);
    }
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *fld = (cs_field_t *)def->input;

      if (fld->location_id == cs_mesh_location_get_id_by_name(N_("cells"))) {

        /* Retrieve the advection field:
           Switch to a cs_nvec3_t representation */
        cs_nvec3_t  adv_vect;
        cs_nvec3(fld->val + 3*cm->c_id, &adv_vect);

        /* Sanity check */
        assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ));

        /* Loop on cell faces */
        for (short int f = 0; f < cm->n_fc; f++)
          fluxes[f] = adv_vect.meas*cm->face[f].meas * _dp3(adv_vect.unitv,
                                                            cm->face[f].unitv);

      }
      else
        bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "Incompatible type of definition.");
    break;

  } /* def_type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the value of the flux of the advection field across the
 *         the dual faces of a cell
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      adv        pointer to a cs_adv_field_t structure
 * \param[in]      time_eval  physical time at which one evaluates the term
 * \param[in, out] fluxes     array of values attached to dual faces of a cell
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_get_cw_dface_flux(const cs_cell_mesh_t     *cm,
                                     const cs_adv_field_t     *adv,
                                     cs_real_t                 time_eval,
                                     cs_real_t                *fluxes)
{
  if (adv == NULL)
    return;

  if (fluxes == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " fluxes array should be allocated before the call.");

  cs_xdef_t  *def = adv->definition;

  /* Sanity checks */
  assert(def != NULL);

  switch (def->type) {

  case CS_XDEF_BY_VALUE:
    {
      /* Retrieve the advection field: Switch to a cs_nvec3_t representation */
      cs_real_3_t  cell_vector;
      cs_xdef_cw_eval_vector_by_val(cm, time_eval, def->input, cell_vector);
      cs_nvec3_t  adv_vect;
      cs_nvec3(cell_vector, &adv_vect);

      /* Sanity check */
      assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_PEQ));

      /* Loop on cell edges */
      for (short int e = 0; e < cm->n_ec; e++)
        fluxes[e] = adv_vect.meas * cm->dface[e].meas *
          _dp3(adv_vect.unitv, cm->dface[e].unitv);

    }
    break;

  case CS_XDEF_BY_ANALYTIC_FUNCTION:
    {
      assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_EFQ |
                          CS_CDO_LOCAL_PFQ));

      cs_quadrature_type_t  qtype = cs_xdef_get_quadrature(def);

      /* Loop on cell edges */
      for (short int e = 0; e < cm->n_ec; e++) {

        const cs_quant_t  edge = cm->edge[e];

        /* Two triangles composing the dual face inside a cell */
        const short int  f0 = cm->e2f_ids[2*e];
        const cs_nvec3_t  sef0 = cm->sefc[2*e];
        const cs_quant_t  qf0 = cm->face[f0];

        const short int  f1 = cm->e2f_ids[2*e+1];
        const cs_nvec3_t  sef1 = cm->sefc[2*e+1];
        const cs_quant_t  qf1 = cm->face[f1];

        fluxes[e] = 0.;
        switch (qtype) {

        case CS_QUADRATURE_NONE:
        case CS_QUADRATURE_BARY:
        case CS_QUADRATURE_BARY_SUBDIV:
          {
            cs_real_3_t  xg[2], adv_xg[2];

            for (int k = 0; k < 3; k++) {
              const double  xec = cm->xc[k] + edge.center[k];
              xg[0][k] = xec + qf0.center[k];
              xg[0][k] *= cs_math_onethird;
              xg[1][k] = xec + qf1.center[k];
              xg[1][k] *= cs_math_onethird;
            }

            cs_xdef_cw_eval_at_xyz_by_analytic(cm,
                                               2, (const cs_real_t *)xg,
                                               time_eval,
                                               def->input,
                                               (cs_real_t *)adv_xg);

            fluxes[e] = sef0.meas * _dp3(adv_xg[0], sef0.unitv)
              + sef1.meas * _dp3(adv_xg[1], sef1.unitv);
          }
          break;

        case CS_QUADRATURE_HIGHER:
          {
            cs_real_t  w[2];
            cs_real_3_t  gpts[6], eval[6];

            /* Two triangles composing the dual face inside a cell */
            cs_quadrature_tria_3pts(edge.center, qf0.center, cm->xc,
                                    sef0.meas,
                                    gpts, w);

            /* Evaluate the field at the three quadrature points */
            cs_quadrature_tria_3pts(edge.center, qf1.center, cm->xc,
                                    sef1.meas,
                                    gpts + 3, w + 1);

            cs_xdef_cw_eval_at_xyz_by_analytic(cm,
                                               6, (const cs_real_t *)gpts,
                                               time_eval,
                                               def->input,
                                               (cs_real_t *)eval);

            cs_real_t  add0 = 0, add1 = 0;
            for (int p = 0; p < 3; p++) add0 += _dp3(eval[p], sef0.unitv);
            add0 *= w[0];
            for (int p = 0; p < 3; p++) add1 += _dp3(eval[p+3], sef1.unitv);
            add1 *= w[1];

            fluxes[e] = add0 + add1;
          }
          break;

        case CS_QUADRATURE_HIGHEST:  /* Not implemented yet */
        default:
          bft_error(__FILE__, __LINE__, 0, " Invalid type of quadrature.");
          break;

        }  /* switch type of quadrature */

      }  /* Loop on cell edges */

    }  /* definition by analytic function */
    break;

  case CS_XDEF_BY_ARRAY:
    {
      cs_xdef_array_input_t  *input = (cs_xdef_array_input_t *)def->input;

      /* Test if location has at least the pattern of the reference support */
      if (cs_flag_test(input->loc, cs_flag_dual_face_byc)) {

        /* Sanity check */
        assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PE));

        const cs_real_t  *flux_array =
          input->values + cs_cdo_connect->c2e->idx[cm->c_id];
        for (short int e = 0; e < cm->n_ec; e++)
          fluxes[e] = flux_array[e];

      }
      else if (cs_flag_test(input->loc, cs_flag_primal_cell)) {

        /* Retrieve the advection field:
           Switch to a cs_nvec3_t representation */
        cs_nvec3_t  adv_vect;
        cs_nvec3(input->values + 3*cm->c_id, &adv_vect);

        /* Sanity check */
        assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_PEQ));

        /* Loop on cell edges */
        for (short int e = 0; e < cm->n_ec; e++)
          fluxes[e] = adv_vect.meas*cm->dface[e].meas*_dp3(adv_vect.unitv,
                                                           cm->dface[e].unitv);

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid support for evaluating the advection field %s"
                  " at the cell center of cell %d.", adv->name, cm->c_id);
    }
    break;

  case CS_XDEF_BY_FIELD:
    {
      cs_field_t  *f = (cs_field_t *)def->input;

      if (f->location_id == cs_mesh_location_get_id_by_name(N_("cells"))) {

        /* Retrieve the advection field:
           Switch to a cs_nvec3_t representation */
        cs_nvec3_t  adv_vect;
        cs_nvec3(f->val + 3*cm->c_id, &adv_vect);

        /* Sanity check */
        assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_PEQ));

        /* Loop on cell edges */
        for (short int e = 0; e < cm->n_ec; e++)
          fluxes[e] = adv_vect.meas*cm->dface[e].meas*_dp3(adv_vect.unitv,
                                                           cm->dface[e].unitv);

      }
      else
        bft_error(__FILE__, __LINE__, 0, "%s: TODO.", __func__);

    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "Incompatible type of definition.");
    break;

  } /* def_type */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   For each cs_adv_field_t structures, update the values of the
 *          related field(s)
 *
 * \param[in]  t_eval     physical time at which one evaluates the term
 * \param[in]  cur2prev   true or false
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_field_update(cs_real_t    t_eval,
                          bool         cur2prev)
{
  for (int i = 0; i < _n_adv_fields; i++) {

    cs_adv_field_t  *adv = _adv_fields[i];

    /* Sanity checks */
    assert(adv != NULL);
    assert(adv->bdy_field_id > -1 && adv->cell_field_id > -1);

    /* GWF and NAVSTO type advection fields are updated elsewhere
       except if there is a field defined at vertices */

    if (adv->type == CS_ADVECTION_FIELD_USER) {

      /* Field storing the boundary normal flux */
      cs_field_t  *bfld = cs_field_by_id(adv->bdy_field_id);

      /* Copy current field values to previous values */
      if (cur2prev)
        cs_field_current_to_previous(bfld);

      /* Set the new values */
      cs_advection_field_across_boundary(adv, t_eval, bfld->val);

      /* Field stored at cell centers */
      cs_field_t  *cfld = cs_field_by_id(adv->cell_field_id);

      /* Copy current field values to previous values */
      if (cur2prev)
        cs_field_current_to_previous(cfld);

      /* Set the new values */
      cs_advection_field_in_cells(adv, t_eval, cfld->val);

    }

    if (adv->vtx_field_id > -1) { /* Field stored at vertices */

      cs_field_t  *vfld = cs_field_by_id(adv->vtx_field_id);

      /* Copy current field values to previous values */
      if (cur2prev)
        cs_field_current_to_previous(vfld);

      /* Set the new values */
      cs_advection_field_at_vertices(adv, t_eval, vfld->val);

    }

  } /* Loop on advection fields */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the Peclet number in each cell
 *
 * \param[in]      adv        pointer to the advection field struct.
 * \param[in]      diff       pointer to the diffusion property struct.
 * \param[in]      t_eval     time at which one evaluates the advection field
 * \param[in, out] peclet     pointer to an array storing Peclet number
 */
/*----------------------------------------------------------------------------*/

void
cs_advection_get_peclet(const cs_adv_field_t     *adv,
                        const cs_property_t      *diff,
                        cs_real_t                 t_eval,
                        cs_real_t                 peclet[])
{
  cs_real_t  ptymat[3][3];
  cs_real_3_t  ptydir;
  cs_nvec3_t  adv_c;

  const bool  pty_uniform = cs_property_is_uniform(diff);
  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;

  assert(peclet != NULL);  /* Sanity check */

  /* Get the value of the material property at the first cell center */
  if (pty_uniform)
    cs_property_get_cell_tensor(0, t_eval, diff, false, ptymat);

  for (cs_lnum_t c_id = 0; c_id < cdoq->n_cells; c_id++) {

    /* Get the value of the material property at the cell center */
    if (!pty_uniform)
      cs_property_get_cell_tensor(c_id, t_eval, diff, false, ptymat);

    cs_advection_field_get_cell_vector(c_id, adv, &adv_c);

    cs_real_t  hc = pow(cdoq->cell_vol[c_id], cs_math_onethird);

    cs_math_33_3_product((const cs_real_t (*)[3])ptymat, adv_c.unitv, ptydir);

    peclet[c_id] = hc * adv_c.meas / _dp3(adv_c.unitv, ptydir);

  }  /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the divergence of the advection field at vertices
 *          Useful for CDO Vertex-based schemes
 *
 * \param[in]      adv         pointer to the advection field struct.
 * \param[in]      t_eval      time at which one evaluates the advection field
 *
 * \return a pointer to an array storing the result
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_advection_field_divergence_at_vertices(const cs_adv_field_t     *adv,
                                          cs_real_t                 t_eval)
{
  CS_UNUSED(t_eval); /* Useful in case of analytic definition */

  cs_real_t  *divergence = NULL;

  if (adv == NULL)
    return divergence;

  const cs_cdo_quantities_t  *cdoq = cs_cdo_quant;
  const cs_cdo_connect_t  *connect = cs_cdo_connect;
  const cs_adjacency_t  *f2e = connect->f2e;
  const cs_adjacency_t  *e2v = connect->e2v;

  BFT_MALLOC(divergence, cdoq->n_vertices, cs_real_t);
  memset(divergence, 0, sizeof(cs_real_t)*cdoq->n_vertices);

  { /* Volume part */
    const cs_xdef_t  *def = adv->definition;

    switch (def->type) {

    case CS_XDEF_BY_ARRAY:
      {
        cs_xdef_array_input_t  *ai = (cs_xdef_array_input_t *)def->input;

        if (cs_flag_test(ai->loc, cs_flag_dual_face_byc)) {

          const cs_adjacency_t  *c2e = connect->c2e;

          for (cs_lnum_t  c_id = 0; c_id < cdoq->n_cells; c_id++) {
            for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++) {

              const cs_lnum_t  e_id = c2e->ids[j];
              const cs_real_t  flx = ai->values[j];
              const cs_lnum_t  eshift = 2*e_id;
              const cs_lnum_t  v0 = e2v->ids[eshift];
              const cs_lnum_t  v1 = e2v->ids[eshift+1];
              const short int  sgn = e2v->sgn[eshift];

              divergence[v0] += -sgn*flx;
              divergence[v1] +=  sgn*flx;

            } /* Loop on cell edges */
          } /* Loop on cells */

        }
        else
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Invalid location for the array.", __func__);

      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.", __func__);

    } /* End of switch */

  } /* Volume part */

  /* Boundary part */
  if (adv->n_bdy_flux_defs > 0) {

    for (int def_id = 0; def_id < adv->n_bdy_flux_defs; def_id++) {

      const cs_xdef_t  *def = adv->bdy_flux_defs[def_id];
      const cs_zone_t  *z = cs_boundary_zone_by_id(def->z_id);
      assert(def->support == CS_XDEF_SUPPORT_BOUNDARY);

      switch (def->type) {

      case CS_XDEF_BY_VALUE:
        {
          const cs_real_t  *constant_val = (cs_real_t *)def->input;

          for (cs_lnum_t id = 0; id < z->n_elts; id++) {

            const cs_lnum_t  bf_id = (z->elt_ids == NULL) ? id: z->elt_ids[id];

            _fill_uniform_boundary_flux(cdoq, f2e, e2v, bf_id, constant_val[0],
                                        divergence);

          } /* Loop on boundary faces */
        }
        break; /* by value */

      case CS_XDEF_BY_ARRAY:
        {
          const cs_xdef_array_input_t  *input =
            (cs_xdef_array_input_t *)def->input;
          const cs_real_t  *val = input->values;

          assert(input->stride == 1);
          assert(z->id == 0);

          if (cs_flag_test(input->loc, cs_flag_primal_face)) {

            for (cs_lnum_t bf_id = 0; bf_id < cdoq->n_b_faces; bf_id++)
              _fill_uniform_boundary_flux(cdoq, f2e, e2v, bf_id, val[bf_id],
                                          divergence);

          }
          else if (cs_flag_test(input->loc, cs_flag_dual_closure_byf)) {

            const cs_adjacency_t  *bf2v = connect->bf2v;
            const cs_lnum_t  *idx = bf2v->idx;
            assert(idx == input->index);

            for (cs_lnum_t bf_id = 0; bf_id < cdoq->n_b_faces; bf_id++) {
              for (cs_lnum_t i = idx[bf_id]; i < idx[bf_id+1]; i++) {
                const cs_lnum_t  v_id = bf2v->ids[i];
                divergence[v_id] += val[i];
              }
            }

          }
          else
            bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.", __func__);

        }
        break; /* by_array */

      default:
        bft_error(__FILE__, __LINE__, 0, " %s: Invalid case", __func__);
        break;

      } /* End of switch on the type of definition */
    }   /* Loop on definitions */

  }
  else {

    /* Handle the boundary flux */
    const cs_field_t  *bflx =
      cs_advection_field_get_field(adv, CS_MESH_LOCATION_BOUNDARY_FACES);

    for (cs_lnum_t bf_id = 0; bf_id < cdoq->n_b_faces; bf_id++) {

      const cs_real_t  face_flx = bflx->val[bf_id];
      const cs_real_t  invsurf = 1./cdoq->b_face_surf[bf_id];
      const cs_lnum_t  f_id = cdoq->n_i_faces + bf_id;

      for (cs_lnum_t i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

        const cs_lnum_t  e_id = f2e->ids[i];
        const cs_lnum_t  eshift = 2*e_id;
        const cs_lnum_t  v0 = e2v->ids[eshift];
        const cs_lnum_t  v1 = e2v->ids[eshift+1];
        const double  tef = cs_math_surftri(cdoq->vtx_coord + 3*v0,
                                            cdoq->vtx_coord + 3*v1,
                                            cdoq->b_face_center + 3*bf_id);

        const double  weighted_flux = 0.5 * tef * invsurf * face_flx;

        divergence[v0] += weighted_flux;
        divergence[v1] += weighted_flux;

      } /* Loop on face edges */

    } /* Loop on boundary faces */

  } /* Boundary part */

#if defined(HAVE_MPI) /* Parallel synchronisation if needed */
  if (cs_glob_n_ranks > 1)
    cs_interface_set_sum(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                         cdoq->n_vertices,
                         1,             /* stride */
                         false,         /* interlace (not useful here) */
                         CS_REAL_TYPE,
                         divergence);
#endif  /* MPI */

  return divergence;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
