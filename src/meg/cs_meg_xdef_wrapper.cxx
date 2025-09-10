/*============================================================================
 * MEG (Mathematical Expression Generator) functions xdef wrapper
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "base/cs_array.h"
#include "base/cs_base.h"
#include "base/cs_mem.h"
#include "meg/cs_meg_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "meg/cs_meg_xdef_wrapper.h"

BEGIN_C_DECLS

/*============================================================================
 * Local variables
 *============================================================================*/

static int  _n_meg_defs = 0;
static cs_meg_xdef_input_t **_meg_defs = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy MEG xdef wrapper arrays
 */
/*----------------------------------------------------------------------------*/

static void
_meg_xdef_wrapper_finalize(void)
{
  for (int i = 0; i < _n_meg_defs; i++)
    CS_FREE(_meg_defs[i]);

  CS_FREE(_meg_defs);
  _n_meg_defs = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get number of elements based on location type.
 */
/*----------------------------------------------------------------------------*/

static cs_lnum_t
_meg_alloc_size_from_location
(
  cs_mesh_location_type_t location
)
{
  cs_lnum_t retval = 0;

  switch(location) {
  case CS_MESH_LOCATION_CELLS:
    retval = cs_glob_mesh->n_cells;
    break;
  case CS_MESH_LOCATION_BOUNDARY_FACES:
    retval = cs_glob_mesh->n_b_faces;
    break;
  case CS_MESH_LOCATION_VERTICES:
    retval = cs_glob_mesh->n_vertices;
    break;
  default:
    assert(0);
    break;
  }

  return retval;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add a MEG function xdef wrapper input data. Allocated memory is
 *  deleted by cs_meg_xdef_wrapper_finalize
 *
 * \param[in] type              type of meg function linked to this input
 * \param[in] z_id              id of zone on which this function is defined
 * \param[in] stride            stride of data
 * \param[in] name              name related to function
 * \param[in] additional_data   additional data (char *) provided to function,
 *                              such as condition or source type
 *
 * \returns pointer to newly allocated input data structure
 */
/*----------------------------------------------------------------------------*/

cs_meg_xdef_input_t *
cs_meg_xdef_wrapper_add_input(const cs_meg_function_type_t   type,
                              const int                      z_id,
                              const cs_mesh_location_type_t  location,
                              const int                      stride,
                              const char                    *name,
                              const char                    *additional_data)
{
  if (_n_meg_defs == 0)
    cs_base_at_finalize(_meg_xdef_wrapper_finalize);

  int new_id = _n_meg_defs;
  _n_meg_defs += 1;

  CS_REALLOC(_meg_defs, _n_meg_defs, cs_meg_xdef_input_t *);

  cs_meg_xdef_input_t *d = NULL;
  CS_MALLOC(d, 1, cs_meg_xdef_input_t);

  d->type = type;
  d->z_id = z_id;
  d->location = location;
  d->stride = stride;

  if (name == NULL || (name != NULL && strlen(name) == 0))
    bft_error(__FILE__, __LINE__, 0,
              _("%s: empty name provided."), __func__);
  snprintf(d->name, 511, "%s", name);
  d->name[511] = '\0';

  if (additional_data != NULL) {
    snprintf(d->additional_data, 511, "%s", additional_data);
    d->additional_data[511] = '\0';
  }
  else
    d->additional_data[0] = '\0';

  _meg_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper function allowing to call MEG functions by xdef structres.
 * This is done by using the cs_xdef_analytic_function type.
 *
 * \param[in] time          when ?
 * \param[in] n_elts        number of elements to consider
 * \param[in] elt_ids       list of elements ids (in coords and retval)
 * \param[in] coords        where ?
 * \param[in] dense_output  perform an indirection in retval or not
 * \param[in] input         pointer to cs_meg_xdef_input_t
 * \param[in] retval        resultint value(s). Must be allocated
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_xdef_wrapper(cs_real_t         time,
                    cs_lnum_t         n_elts,
                    const cs_lnum_t  *elt_ids,
                    const cs_real_t  *coords,
                    bool              dense_output,
                    void             *input,
                    cs_real_t        *retval)
{
  CS_UNUSED(time);

  const cs_meg_xdef_input_t *_input
    = reinterpret_cast<const cs_meg_xdef_input_t *>(input);

  const cs_real_3_t *_coords = (const cs_real_3_t *)coords;

  cs_real_t *meg_vals = NULL;
  /* Volume function takes as an input arrays over the entire domain */
  if (_input->type == CS_MEG_VOLUME_FUNC) {
    if (dense_output) {
      cs_lnum_t n_elts_alloc = _meg_alloc_size_from_location(_input->location);
      CS_MALLOC(meg_vals, n_elts_alloc * _input->stride, cs_real_t);
    }
    else
      meg_vals = retval;
  }
  else {
    if (dense_output)
      meg_vals = retval;
    else {
      cs_lnum_t n_elts_alloc = _meg_alloc_size_from_location(_input->location);
      CS_MALLOC(meg_vals, n_elts_alloc * _input->stride, cs_real_t);
    }
  }

  switch(_input->type) {
  case CS_MEG_BOUNDARY_FUNC:
    {
      const cs_zone_t *z = cs_boundary_zone_by_id(_input->z_id);
      cs_meg_boundary_function(z->name,
                               n_elts,
                               elt_ids,
                               _coords,
                               _input->name,
                               _input->additional_data,
                               meg_vals);

      /* Copy values to retval, knowing that "meg_vals" is dense! */
      if (!dense_output) {
        cs_array_real_copy_subset(n_elts,
                                  _input->stride,
                                  elt_ids,
                                  CS_ARRAY_SUBSET_OUT,
                                  meg_vals,
                                  retval);
        CS_FREE(meg_vals);
      }
    }
    break;
  case CS_MEG_VOLUME_FUNC:
    {
      const cs_zone_t *z = cs_volume_zone_by_id(_input->z_id);
      cs_meg_volume_function(z->name,
                             n_elts,
                             elt_ids,
                             _coords,
                             _input->name,
                             &(meg_vals));
      if (dense_output) {
        cs_array_real_copy_subset(n_elts,
                                  _input->stride,
                                  elt_ids,
                                  CS_ARRAY_SUBSET_IN,
                                  meg_vals,
                                  retval);
        CS_FREE(meg_vals);
      }
    }
    break;
  case CS_MEG_INITIALIZATION_FUNC:
    {
      const cs_zone_t *z = cs_volume_zone_by_id(_input->z_id);
      cs_meg_initialization(z->name,
                            n_elts,
                            elt_ids,
                            _coords,
                            _input->name,
                            meg_vals);

      /* Copy values to retval, knowing that "meg_vals" is dense! */
      if (!dense_output) {
        cs_array_real_copy_subset(n_elts,
                                  _input->stride,
                                  elt_ids,
                                  CS_ARRAY_SUBSET_OUT,
                                  meg_vals,
                                  retval);
        CS_FREE(meg_vals);
      }
    }
    break;
  case CS_MEG_SOURCE_TERM_FUNC:
    {
      const cs_zone_t *z = cs_volume_zone_by_id(_input->z_id);
      cs_meg_source_terms(z->name,
                          n_elts,
                          elt_ids,
                          _coords,
                          _input->name,
                          _input->additional_data,
                          meg_vals);

      /* Copy values to retval, knowing that "meg_vals" is dense! */
      if (!dense_output) {
        cs_array_real_copy_subset(n_elts,
                                  _input->stride,
                                  elt_ids,
                                  CS_ARRAY_SUBSET_OUT,
                                  meg_vals,
                                  retval);
        CS_FREE(meg_vals);
      }
    }
    break;

  case CS_MEG_CALCULATOR_FUNC:
    {
      cs_meg_post_calculator(_input->name,
                             n_elts,
                             elt_ids,
                             _coords,
                             meg_vals);
      if (!dense_output) {
        cs_array_real_copy_subset(n_elts,
                                  _input->stride,
                                  elt_ids,
                                  CS_ARRAY_SUBSET_OUT,
                                  meg_vals,
                                  retval);
        CS_FREE(meg_vals);
      }
    }
    break;

  default:
    {
      bft_error
        (__FILE__, __LINE__, 0,
         _("\"%s\" was called with an incompatible type of MEG function."),
         __func__);
    }
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Wrapper function allowing to call MEG functions by xdef structres.
 * This is done by using the cs_xdef_function type.
 *
 * \param[in] time          when ?
 * \param[in] n_elts        number of elements to consider
 * \param[in] elt_ids       list of elements ids (in coords and retval)
 * \param[in] coords        where ?
 * \param[in] dense_output  perform an indirection in retval or not
 * \param[in] input         pointer to cs_meg_xdef_input_t
 * \param[in] retval        resultint value(s). Must be allocated
 */
/*----------------------------------------------------------------------------*/

void
cs_meg_xdef_eval_func_wrapper(int              location_id,
                              cs_lnum_t        n_elts,
                              const cs_lnum_t *elt_ids,
                              void            *input,
                              void            *retval)
{
    const cs_mesh_location_type_t loc_type
      = cs_mesh_location_get_type(location_id);

    const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

    const cs_real_3_t *base_coords = nullptr;

    if (loc_type == CS_MESH_LOCATION_CELLS)
      base_coords = mq->cell_cen;
    else if (loc_type == CS_MESH_LOCATION_INTERIOR_FACES)
      base_coords = mq->i_face_cog;
    else if (loc_type == CS_MESH_LOCATION_BOUNDARY_FACES)
      base_coords = mq->b_face_cog;
    else if (loc_type == CS_MESH_LOCATION_VERTICES)
      base_coords = reinterpret_cast<const cs_real_3_t *>(cs_glob_mesh->vtx_coord);

  const cs_meg_xdef_input_t *_input
    = reinterpret_cast<const cs_meg_xdef_input_t *>(input);

  if (_input->type == CS_MEG_CALCULATOR_FUNC) {
    cs_meg_post_calculator(_input->name,
                           n_elts,
                           elt_ids,
                           base_coords,
                           static_cast<cs_real_t *>(retval));
  }
  else
    bft_error(__FILE__, __LINE__, 0,
              _("%s: cannot be called for MEG function of type %d\n"),
              __func__, _input->type);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
