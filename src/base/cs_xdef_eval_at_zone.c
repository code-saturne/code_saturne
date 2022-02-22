/*============================================================================
 * Boundary condition handling.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_field.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_quantities.h"
#include "cs_parameters.h"
#include "cs_boundary_zone.h"
#include "cs_volume_zone.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_xdef_eval_at_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_boundary_conditions.c
        Boundary condition handling.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values defined by a cs_xdef_t structure associated to
 *          a given volume or boundary zone.
 *
 * \param[in]      def                 pointer to a cs_xdef_t structure
 * \param[in]      caller_description  information string (for error logging)
 * \param[in]      t_eval              time at which values are evaluated
 * \param[in]      dense               if true, values on zone location;
 *                                     if false, values on parent location.
 * \param[in, out] values              pointer to the array of values
 */
/*----------------------------------------------------------------------------*/

void
cs_xdef_eval_at_zone(const cs_xdef_t   *def,
                     const char        *caller_description,
                     cs_real_t          t_eval,
                     bool               dense,
                     cs_real_t         *values)
{
  const cs_zone_t  *z = NULL;
  int base_location = CS_MESH_LOCATION_NONE;

  const char *description
    = caller_description != NULL ? caller_description : cs_empty_string;

  switch (def->support) {
  case CS_XDEF_SUPPORT_BOUNDARY:
    z = cs_boundary_zone_by_id(def->z_id);
    base_location = CS_MESH_LOCATION_BOUNDARY_FACES;
    break;
  case CS_XDEF_SUPPORT_VOLUME:
    z = cs_volume_zone_by_id(def->z_id);
    base_location = CS_MESH_LOCATION_CELLS;
    break;
  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: %s\n"
                " zone %s, definition type %s;\n"
                " cs_xdef_t %p->support expected as"
                " CS_XDEF_SUPPORT_BOUNDARY or CS_XDEF_SUPPORT_VOLUME,\n"
                " not %d."),
              __func__, description,
              z->name, cs_xdef_type_get_name(def->type), z,
              def->support);
  }

  const cs_lnum_t  *elt_ids = z->elt_ids;
  const cs_lnum_t  n_elts = z->n_elts;
  const cs_lnum_t def_dim = def->dim;

  int use_threads = (n_elts > CS_THR_MIN) ? 1 : 0;

  switch(def->type) {

  case CS_XDEF_BY_ANALYTIC_FUNCTION: /* Analytic function */
    {
      cs_xdef_analytic_context_t  *ac
        = (cs_xdef_analytic_context_t *)def->context;

      const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
      const cs_real_3_t *base_coords = NULL;
      if (def->support == CS_XDEF_SUPPORT_BOUNDARY)
        base_coords = (const cs_real_3_t *)(mq->b_face_cog);
      else if (def->support == CS_XDEF_SUPPORT_VOLUME)
        base_coords = (const cs_real_3_t *)(mq->cell_cen);

      ac->func(t_eval,
               n_elts, elt_ids, (const cs_real_t *)base_coords,
               dense,
               ac->input,
               values);
    }
    break;

  case CS_XDEF_BY_ARRAY: /* array defined on full location or same zone */
    {
      cs_xdef_array_context_t  *ac = (cs_xdef_array_context_t *)def->context;

      if (ac->stride != def_dim)
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: %s:\n"
                    " zone %s, definition type %s;\n"
                    " array stride (%d) does not match dimension (%d)."),
                  __func__, description,
                  z->name, cs_xdef_type_get_name(def->type),
                  ac->stride, def_dim);

      if (ac->z_id == 0) {  /* Full location */
        if (dense) {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[i*def_dim + coo_id]
                = ac->values[elt_id*def_dim + coo_id];
            }
          }
        }
        else {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[elt_id*def_dim + coo_id]
                = ac->values[elt_id*def_dim + coo_id];
            }
          }
        }
      }
      else if (ac->z_id == def->z_id) {  /* Same zone */
        if (dense) {
          cs_lnum_t _n_elts = n_elts * def_dim;
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < _n_elts; i++) {
            values[i] = ac->values[i];
          }
        }
        else {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[elt_id*def_dim + coo_id]
                = ac->values[i*def_dim + coo_id];
            }
          }
        }
      }
      else
        bft_error
          (__FILE__, __LINE__, 0,
           _(" %s: %s:\n"
             " zone %s, type %s;\n"
             " array zone id should be %d for this zone or 0 for all elements,"
             " not %d."),
           __func__, description, z->name, cs_xdef_type_get_name(def->type),
           def->z_id, ac->z_id);
    }
    break;

  case CS_XDEF_BY_DOF_FUNCTION: /* DOF function */
    {
      cs_xdef_dof_context_t  *cx = (cs_xdef_dof_context_t *)def->context;

      cs_flag_t test_flag = 0;
      switch(def->support) {
      case CS_XDEF_SUPPORT_BOUNDARY:
        test_flag = cs_flag_primal_face;
        break;
      case CS_XDEF_SUPPORT_VOLUME:
        test_flag = cs_flag_primal_cell;
        break;
      default:
        break;
      }

      if (cs_flag_test(cx->loc, test_flag) == false)
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid location for definition by DoFs.\n",
                  __func__);

      cx->func(n_elts, elt_ids, dense,
               cx->input,
               values);
    }
    break;

  case CS_XDEF_BY_FIELD: /* array defined by associated field */
    {
      cs_field_t  *df = (cs_field_t *)def->context;

      if (df->dim != def_dim)
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: %s\n"
                    " zone %s, type %s;\n"
                    " associated field %s dimension (%d)\n"
                    " does not match dimension (%d)."),
                  __func__, description, z->name, cs_xdef_type_get_name(def->type),
                  df->name, df->dim, def_dim);

      if (df->location_id == base_location) {  /* Full location */
        if (dense) {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[i*def_dim + coo_id]
                = df->val[elt_id*def_dim + coo_id];
            }
          }
        }
        else {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[elt_id*def_dim + coo_id]
                = df->val[elt_id*def_dim + coo_id];
            }
          }
        }
      }
      else if (df->location_id == z->location_id) {  /* Same zone */
        if (dense) {
          cs_lnum_t _n_elts = n_elts * def_dim;
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < _n_elts; i++) {
            values[i] = df->val[i];
          }
        }
        else {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[elt_id*def_dim + coo_id]
                = df->val[i*def_dim + coo_id];
            }
          }
        }
      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: %s\n"
                    " zone %s, type %s;\n"
                    " associated field %s location id should be %d (%s)\n"
                    " or 0 (%s)."),
                  __func__, description, z->name, cs_xdef_type_get_name(def->type),
                  df->name, df->location_id,
                  cs_mesh_location_get_name(df->location_id),
                  cs_mesh_location_get_name(CS_MESH_LOCATION_BOUNDARY_FACES));

    }
    break;

  case CS_XDEF_BY_QOV:
    {
      const cs_real_t  *constant_val = (cs_real_t *)def->context;
      const cs_real_t factor = 1. / z->f_measure;

      if (def_dim > 1) {
        if (dense) {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[i*def_dim + coo_id] = constant_val[coo_id] * factor;
            }
          }
        }
        else {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[elt_id*def_dim + coo_id] = constant_val[coo_id] * factor;
            }
          }
        }
      }
      else {
        const cs_real_t v = constant_val[0] * factor;
        if (dense) {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            values[i] = v;
          }
        }
        else {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            values[elt_id] = v;
          }
        }
      }
    }
    break;

  case CS_XDEF_BY_VALUE:
    {
      const cs_real_t  *constant_val = (cs_real_t *)def->context;

      if (def_dim > 1) {
        if (dense) {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[i*def_dim + coo_id] = constant_val[coo_id];
            }
          }
        }
        else {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {
              values[elt_id*def_dim + coo_id] = constant_val[coo_id];
            }
          }
        }
      }
      else {
        if (dense) {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            values[i] = constant_val[0];
          }
        }
        else {
#         pragma omp parallel for if (use_threads)
          for (cs_lnum_t i = 0; i < n_elts; i++) {
            const cs_lnum_t  elt_id = elt_ids[i];
            values[elt_id] = constant_val[0];
          }
        }
      }
    }
    break;

  case CS_XDEF_BY_TIME_FUNCTION:
    if (def_dim == 1)
      cs_xdef_eval_scalar_at_cells_by_time_func(z->n_elts,
                                                z->elt_ids,
                                                dense,
                                                NULL, /* m */
                                                NULL, /* connect */
                                                NULL, /* quant */
                                                t_eval,
                                                def->context,
                                                values);
    else if (def_dim == 3)
      cs_xdef_eval_vector_at_cells_by_time_func(z->n_elts,
                                                z->elt_ids,
                                                dense,
                                                NULL, /* m */
                                                NULL, /* connect */
                                                NULL, /* quant */
                                                t_eval,
                                                def->context,
                                                values);
    else if (def_dim == 6)
      cs_xdef_eval_tensor_at_cells_by_time_func(z->n_elts,
                                                z->elt_ids,
                                                dense,
                                                NULL, /* m */
                                                NULL, /* connect */
                                                NULL, /* quant */
                                                t_eval,
                                                def->context,
                                                values);
    else
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: %s\n"
                  " zone %s, type %s;\n"
                  " evaluation by time function not handled for dimension %d."),
                  __func__, description, z->name,
                cs_xdef_type_get_name(def->type),
                def_dim);

    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: %s\n"
                " zone %s, type %s;\n"
                " unsupported type of definition."),
              __func__, description, z->name, cs_xdef_type_get_name(def->type));

  } /* switch on def_type */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
