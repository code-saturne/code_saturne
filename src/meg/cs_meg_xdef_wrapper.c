/*============================================================================
 * MEG (Mathematical Expression Generator) functions xdef wrapper
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#include "cs_array.h"
#include "cs_base.h"
#include "cs_meg_prototypes.h"
#include "cs_meg_xdef_wrapper.h"

BEGIN_C_DECLS


/*============================================================================
 * Public functions
 *============================================================================*/

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

  const cs_meg_xdef_input_t *_input = input;

  const cs_real_3_t *_coords = (const cs_real_3_t *)coords;


  cs_real_t *meg_vals = NULL;
  /* Volume function takes as an input arrays over the entire domain */
  if (_input->type == CS_MEG_VOLUME_FUNC) {
    if (dense_output)
      bft_error(__FILE__, __LINE__, 0,
                _("Error: MEG Volume function cannot be used with a dense output.\n"));
    else
      meg_vals = retval;
  }
  else {
    if (dense_output)
      meg_vals = retval;
    else
      BFT_MALLOC(meg_vals, cs_glob_mesh->n_cells * _input->stride, cs_real_t);
  }

  switch(_input->type) {
  case CS_MEG_BOUNDARY_FUNC:
    {
      const cs_zone_t *z = cs_boundary_zone_by_id(_input->z_id);
      cs_meg_boundary_function(z->name,
                               n_elts,
                               elt_ids,
                               _coords,
                               _input->input_name,
                               _input->additional_data,
                               meg_vals);

      /* Copy values to retval, knowing that "meg_vals" is dense! */
      if (!dense_output) {
        cs_array_real_copy_subset(n_elts,
                                  _input->stride * n_elts,
                                  elt_ids,
                                  CS_ARRAY_SUBSET_OUT,
                                  meg_vals,
                                  retval);
        BFT_FREE(meg_vals);
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
                             _input->input_name,
                             &(meg_vals));
    }
    break;
  case CS_MEG_INITIALIZATION_FUNC:
    {
      const cs_zone_t *z = cs_volume_zone_by_id(_input->z_id);
      cs_meg_initialization(z->name,
                            n_elts,
                            elt_ids,
                            _coords,
                            _input->input_name,
                            meg_vals);

      /* Copy values to retval, knowing that "meg_vals" is dense! */
      if (!dense_output) {
        cs_array_real_copy_subset(n_elts,
                                  _input->stride * n_elts,
                                  elt_ids,
                                  CS_ARRAY_SUBSET_OUT,
                                  meg_vals,
                                  retval);
        BFT_FREE(meg_vals);
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
                          _input->input_name,
                          _input->additional_data,
                          meg_vals);

      /* Copy values to retval, knowing that "meg_vals" is dense! */
      if (!dense_output) {
        cs_array_real_copy_subset(n_elts,
                                  _input->stride * n_elts,
                                  elt_ids,
                                  CS_ARRAY_SUBSET_OUT,
                                  meg_vals,
                                  retval);
        BFT_FREE(meg_vals);
      }
    }
    break;
  default:
    {
      bft_error(__FILE__, __LINE__, 0,
                _("Error: \"%s\" was called with an incompatible type of MEG function.\n"),
                __func__);
    }
    break;
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
