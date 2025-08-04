/*============================================================================
 * Boundary condition handling.
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

#include "base/cs_defs.h"

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

#include <ple_locator.h>

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_ale.h"
#include "base/cs_base.h"
#include "base/cs_base_accel.h"
#include "alge/cs_blas.h"
#include "base/cs_boundary.h"
#include "cfbl/cs_cf_thermo.h"
#include "base/cs_coupling.h"
#include "cdo/cs_equation.h"
#include "base/cs_function.h"
#include "alge/cs_gradient.h"
#include "gui/cs_gui_util.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_flag_check.h"
#include "base/cs_halo.h"
#include "base/cs_math.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_connect.h"
#include "mesh/cs_mesh_quantities.h"
#include "base/cs_parall.h"
#include "base/cs_parameters.h"
#include "pprt/cs_physical_model.h"
#include "base/cs_prototypes.h"
#include "base/cs_post.h"
#include "base/cs_time_control.h"
#include "turb/cs_turbulence_bc.h"
#include "turb/cs_turbulence_model.h"
#include "base/cs_xdef_eval_at_zone.h"

#include "fvm/fvm_nodal.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_boundary_conditions.h"
#include "base/cs_boundary_conditions_priv.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_boundary_conditions.cpp
        Boundary condition handling.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Mapped inlet */
/*--------------*/

typedef struct {

  int             bc_location_id;      /* location id of boundary zone */
  int             source_location_id;  /* location id of source elements */
  cs_real_t       coord_shift[3];      /* coordinates shift relative to
                                          selected boundary faces */
  double          tolerance;           /* search tolerance */

  ple_locator_t  *locator;             /* associated locator */

} cs_bc_map_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static int *_bc_type;

static int *_bc_pm_face_zone;

static int _n_bc_maps = 0;
static cs_bc_map_t *_bc_maps = nullptr;

static int _n_bc_open = 0;
static cs_boundary_conditions_open_t **_bc_open = nullptr;

static cs_real_t *_b_head_loss = nullptr;

/*============================================================================
 * Global variables
 *============================================================================*/

const int *cs_glob_bc_type = nullptr;

cs_boundary_condition_pm_info_t  *cs_glob_bc_pm_info = nullptr;

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

int *
cs_f_boundary_conditions_get_bc_type(void);

void
cs_f_boundary_conditions_mapped_set(int                        field_id,
                                    ple_locator_t             *locator,
                                    cs_mesh_location_type_t    location_type,
                                    int                        enforce_balance,
                                    int                        interpolate,
                                    cs_lnum_t                  n_faces,
                                    const cs_lnum_t           *faces,
                                    cs_real_t                 *balance_w,
                                    int                        nvar,
                                    cs_real_t                 *rcodcl);

void
cs_f_boundary_conditions_get_pointers(int  **itypfb,
                                      int  **izfppp);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Unset a bit from a flag set if present
 *
 * parameters:
 *   flag <-- flag
 *   mask <-- bit mast to unset
 *
 * return:
 *   updated bit mask
 *----------------------------------------------------------------------------*/

static int
_unset_flag(int  mask,
            int  flag)
{
  int retval = (mask | flag) - flag;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return a pointer to equation parameters based on a field or equation name.
 *
 * parameters:
 *   name <-- field or equation name
 *
 * return:
 *   pointer to matching child string
 *----------------------------------------------------------------------------*/

static cs_equation_param_t *
_get_equation_param(const char  *name)
{
  cs_equation_param_t *eqp = nullptr;

  cs_field_t *f = cs_field_by_name_try(name);
  if (f != nullptr)
    eqp = cs_field_get_equation_param(f);

  else
    eqp = cs_equation_param_by_name(name);

  return eqp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build or update map of shifted boundary face coordinates on
 *        cells or boundary faces for automatic interpolation.
 *
 * \warning
 *
 * This function does not currently update the location in case of moving
 * meshes. It could be modifed to do so.
 *
 * \param[in]  map_id      id of defined map
 */
/*----------------------------------------------------------------------------*/

static void
_update_bc_map(int  map_id)
{
  assert(map_id > -1 && map_id < _n_bc_maps);

  cs_bc_map_t *bc_map = _bc_maps + map_id;

  if (bc_map->locator != nullptr)
    return;

  cs_mesh_location_type_t location_type
    = cs_mesh_location_get_type(bc_map->source_location_id);

  cs_lnum_t n_location_elts
    = cs_mesh_location_get_n_elts(bc_map->source_location_id)[0];
  cs_lnum_t n_faces
    = cs_mesh_location_get_n_elts(bc_map->bc_location_id)[0];

  const cs_lnum_t *location_elts
    = cs_mesh_location_get_elt_ids_try(bc_map->source_location_id);
  const cs_lnum_t *faces
    = cs_mesh_location_get_elt_ids_try(bc_map->bc_location_id);

  bc_map->locator = cs_boundary_conditions_map(location_type,
                                               n_location_elts,
                                               n_faces,
                                               location_elts,
                                               faces,
                                               &(bc_map->coord_shift),
                                               0,
                                               bc_map->tolerance);
}

/*----------------------------------------------------------------------------
 * Compute balance at inlet
 *
 * parameters:
 *   f               <-- associated field
 *   m               <-- associated mesh
 *   mq              <-- mesh quantities
 *   enforce_balance <-- balance handling option:
 *                         0: values are simply mapped
 *                         1: values are mapped, then multiplied
 *                            by a constant factor so that their
 *                            surface integral on selected faces
 *                            is preserved (relative to the
 *                            input values)
 *                         2: as 1, but with a boundary-defined
 *                            weight, defined by balance_w
 *                         3: as 1, but with a cell-defined
 *                            weight, defined by balance_w
 *   n_faces         <-- number of selected boundary faces
 *   faces           <-- list of selected boundary faces (0 to n-1),
 *                       or nullptr if no indirection is needed
 *   balance_w       <-- optional balance weight, or nullptr
 *   inlet_sum       --> inlet sum
 *----------------------------------------------------------------------------*/

static void
_inlet_sum(const cs_field_t            *f,
           const cs_mesh_t             *m,
           const cs_mesh_quantities_t  *mq,
           int                          enforce_balance,
           cs_lnum_t                    n_faces,
           const cs_lnum_t             *faces,
           cs_real_t                   *balance_w,
           cs_real_t                    inlet_sum[])
{
  const int dim = f->dim;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_real_t *f_surf = mq->b_face_surf;

  /* Get field's variable id */

  assert(dim <= 9);

  for (cs_lnum_t j = 0; j < dim; j++) {

    cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1 + j*n_b_faces;

    inlet_sum[j] = 0.;

    if (enforce_balance == 1) {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != nullptr) ? faces[i] : i;
        inlet_sum[j] += rcodcl1[f_id]*f_surf[f_id];
      }
    }
    else if (enforce_balance == 2) {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != nullptr) ? faces[i] : i;
        inlet_sum[j] += rcodcl1[f_id]*f_surf[f_id]*balance_w[f_id];
      }
    }
    else if (enforce_balance == 3) {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != nullptr) ? faces[i] : i;
        const cs_lnum_t c_id = m->b_face_cells[f_id];
        inlet_sum[j] += rcodcl1[f_id]*f_surf[f_id]*balance_w[c_id];
      }
    }

  }

  cs_parall_sum(dim, CS_REAL_TYPE, inlet_sum);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the homogeneous Dirichlet BCs
 *
 * \param[in]       mesh        pointer to mesh structure
 * \param[in]       boundaries  pointer to associated boundaries
 * \param[in]       eqp         pointer to a cs_equation_param_t
 * \param[in]       def         pointer to a boundary condition definition
 * \param[in]       icodcl_m    multiplier for assigned icodcl values (1 or -1)
 * \param[in, out]  icodcl      boundary conditions type array
 * \param[in, out]  rcodcl1     boundary conditions values array
 */
/*----------------------------------------------------------------------------*/

static void
_compute_hmg_dirichlet_bc(const cs_mesh_t            *mesh,
                          const cs_boundary_t        *boundaries,
                          const cs_equation_param_t  *eqp,
                          const cs_xdef_t            *def,
                          int                         icodcl_m,
                          int                        *icodcl,
                          double                     *rcodcl1)
{
  assert(eqp->dim == def->dim);

  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_zone_t *bz = cs_boundary_zone_by_id(def->z_id);
  const cs_lnum_t *face_ids = bz->elt_ids;

  cs_boundary_type_t boundary_type = boundaries->default_type;

  if (boundaries->types != nullptr) {
    assert(boundaries->n_boundaries > 0);
    boundary_type =
      boundaries->types[cs_boundary_id_by_zone_id(boundaries, def->z_id)];
  }

  int bc_code = (boundary_type & CS_BOUNDARY_WALL) ? 5 : 1;
  if (boundary_type & CS_BOUNDARY_ROUGH_WALL)
    bc_code = 6;

  bc_code *= icodcl_m;

  for (cs_lnum_t coo_id = 0; coo_id < def->dim; coo_id++) {

    int        *_icodcl  = icodcl + coo_id*n_b_faces;
    cs_real_t  *_rcodcl1 = rcodcl1 + coo_id*n_b_faces;

#   pragma omp parallel for if (bz->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < bz->n_elts; i++) {
      const cs_lnum_t  elt_id = (face_ids == nullptr) ? i : face_ids[i];
      _icodcl[elt_id]  = bc_code;
      _rcodcl1[elt_id] = 0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs
 *
 * \param[in]       mesh         pointer to mesh structure
 * \param[in]       boundaries   pointer to associated boundaries
 * \param[in]       f            pointer to associated field
 * \param[in]       eqp          pointer to a cs_equation_param_t
 * \param[in]       def          pointer to a boundary condition definition
 * \param[in]       description  description string (for error logging)
 * \param[in]       t_eval       time at which one evaluates the boundary cond.
 * \param[in]       icodcl_m     multiplier for assigned icodcl values (1 or -1)
 * \param[in, out]  icodcl       boundary conditions type array
 * \param[in, out]  rcodcl1      boundary conditions values array
 * \param           eval_buf     evaluation bufferb (work array)
 */
/*----------------------------------------------------------------------------*/

static void
_compute_dirichlet_bc(const cs_mesh_t            *mesh,
                      const cs_boundary_t        *boundaries,
                      const cs_field_t           *f,
                      const cs_equation_param_t  *eqp,
                      const cs_xdef_t            *def,
                      const char                 *description,
                      cs_real_t                   t_eval,
                      int                         icodcl_m,
                      int                        *icodcl,
                      double                     *rcodcl1,
                      cs_real_t                   eval_buf[])
{
  CS_NO_WARN_IF_UNUSED(t_eval);

  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_zone_t *bz = cs_boundary_zone_by_id(def->z_id);
  const cs_lnum_t *elt_ids = bz->elt_ids;
  const cs_lnum_t  n_elts = bz->n_elts;
  const cs_lnum_t  def_dim = def->dim;
  assert(eqp->dim == def->dim);

  cs_boundary_type_t boundary_type = boundaries->default_type;

  if (boundaries->types != nullptr) {
    assert(boundaries->n_boundaries > 0);
    boundary_type =
      boundaries->types[cs_boundary_id_by_zone_id(boundaries, def->z_id)];
  }

  int bc_code = (boundary_type & CS_BOUNDARY_WALL) ? 5 : 1;
  if (boundary_type & CS_BOUNDARY_ROUGH_WALL)
    bc_code = 6;
  if (boundary_type & CS_BOUNDARY_CONVECTIVE_INLET)
    bc_code = 13;

  bc_code *= icodcl_m;

  /* Special case for mesh velocity (ALE) wall BC type. */

  if (f->dim == 3) {
    if (   (strcmp(f->name, "mesh_velocity") == 0)
        && (boundary_type & CS_BOUNDARY_WALL))
      bc_code = 1;
  }

  if (f->dim != def_dim)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Boundary condition definition:\n"
                " field %s, zone %s, type %s;\n"
                " dimension %d does not match field dimension (%d)."),
              __func__, f->name, bz->name, cs_xdef_type_get_name(def->type),
              def_dim, f->dim);

  switch(def->type) {

  case CS_XDEF_BY_VALUE:  /* direct application (optimization)
                             for constant value */
    {
      const cs_real_t  *constant_val = (cs_real_t *)def->context;

      for (cs_lnum_t coo_id = 0; coo_id < def_dim; coo_id++) {

        int        *_icodcl  = icodcl + coo_id*n_b_faces;
        cs_real_t  *_rcodcl1 = rcodcl1 + coo_id*n_b_faces;

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  elt_id = elt_ids[i];
          _icodcl[elt_id]  = bc_code;
          _rcodcl1[elt_id] = constant_val[coo_id];
        }

      }
    }
    break;

  default:
    {
      cs_xdef_eval_at_zone(def,
                           description,
                           t_eval,
                           true,  /* dense */
                           eval_buf);

      const cs_lnum_t _dim = def->dim;

      for (cs_lnum_t coo_id = 0; coo_id < _dim; coo_id++) {

        int        *_icodcl  = icodcl + coo_id*n_b_faces;
        cs_real_t  *_rcodcl1 = rcodcl1 + coo_id*n_b_faces;

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  elt_id = elt_ids[i];
          _icodcl[elt_id]  = bc_code;
          _rcodcl1[elt_id] = eval_buf[_dim*i + coo_id];
        }

      }
    }
    break;

  } /* Switch on the type of definition */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the homogeneous Neumann BCs
 *
 * \param[in]       mesh        pointer to mesh structure
 * \param[in]       def         pointer to a boundary condition definition
 * \param[in, out]  icodcl      boundary conditions type array
 * \param[in, out]  rcodcl3     boundary conditions values array
 */
/*----------------------------------------------------------------------------*/

static void
_compute_hmg_neumann_bc(const cs_mesh_t  *mesh,
                        const cs_xdef_t  *def,
                        int              *icodcl,
                        double           *rcodcl3)
{
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_zone_t *bz = cs_boundary_zone_by_id(def->z_id);
  const cs_lnum_t *elt_ids = bz->elt_ids;
  const cs_lnum_t  n_elts = bz->n_elts;

  for (cs_lnum_t coo_id = 0; coo_id < def->dim; coo_id++) {

    int        *_icodcl  = icodcl + coo_id*n_b_faces;
    cs_real_t  *_rcodcl3 = rcodcl3 + coo_id*n_b_faces;

#   pragma omp parallel for if (n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t  elt_id = elt_ids[i];
      _icodcl[elt_id] = 3;
      _rcodcl3[elt_id] = 0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Neumann BCs
 *
 * \param[in]       mesh         pointer to mesh structure
 * \param[in]       eqp          pointer to a cs_equation_param_t
 * \param[in]       def          pointer to a boundary condition definition
 * \param[in]       description  description string (for error logging)
 * \param[in]       t_eval       time at which one evaluates the boundary cond.
 * \param[in, out]  icodcl       boundary conditions type array
 * \param[in, out]  rcodcl3      boundary conditions values array
 * \param           eval_buf     evaluation bufferb (work array)
 */
/*----------------------------------------------------------------------------*/

static void
_compute_neumann_bc(const cs_mesh_t            *mesh,
                    const cs_equation_param_t  *eqp,
                    const cs_xdef_t            *def,
                    const char                 *description,
                    cs_real_t                   t_eval,
                    int                        *icodcl,
                    double                     *rcodcl3,
                    cs_real_t                   eval_buf[])
{
  const cs_lnum_t  n_b_faces = mesh->n_b_faces;
  const cs_zone_t *bz = cs_boundary_zone_by_id(def->z_id);
  const cs_lnum_t *elt_ids = bz->elt_ids;
  const cs_lnum_t  n_elts = bz->n_elts;

  assert(eqp->dim == def->dim);

  const int bc_code = 3;

  switch(def->type) {

  case CS_XDEF_BY_VALUE:
    {
      const int dim = def->dim;
      const cs_real_t *constant_vals = (cs_real_t *)def->context;

      for (cs_lnum_t coo_id = 0; coo_id < dim; coo_id++) {

        const cs_real_t value = constant_vals[coo_id];

        int        *_icodcl = icodcl + coo_id*n_b_faces;
        cs_real_t  *_rcodcl3 = rcodcl3 + coo_id*n_b_faces;

#       pragma omp parallel for if (bz->n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < bz->n_elts; i++) {
          const cs_lnum_t  elt_id = elt_ids[i];
          _icodcl[elt_id]  = bc_code;
          _rcodcl3[elt_id] = value;
        }

      }
    }
    break;

  default:
    {
      cs_xdef_eval_at_zone(def,
                           description,
                           t_eval,
                           true,  /* dense */
                           eval_buf);

      const cs_lnum_t _dim = def->dim;

      for (cs_lnum_t coo_id = 0; coo_id < _dim; coo_id++) {

        int        *_icodcl = icodcl + coo_id*n_b_faces;
        cs_real_t  *_rcodcl3 = rcodcl3 + coo_id*n_b_faces;

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  elt_id = elt_ids[i];
          _icodcl[elt_id]  = bc_code;
          _rcodcl3[elt_id] = eval_buf[_dim*i + coo_id];
        }

      }
    }
    break;

  } /* Switch on the type of definition */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the values of the Dirichlet BCs
 *
 * \param[in]       mesh         pointer to mesh structure
 * \param[in]       boundaries   pointer to associated boundaries
 * \param[in]       f            pointer to associated field
 * \param[in]       eqp          pointer to a cs_equation_param_t
 * \param[in]       def          pointer to a boundary condition definition
 * \param[in]       description  description string (for error logging)
 * \param[in]       t_eval       time at which one evaluates the boundary cond.
 * \param[in, out]  icodcl       boundary conditions type array
 * \param[in, out]  rcodcl1      boundary conditions values array
 * \param[in, out]  rcodcl2      boundary conditions values array
 * \param           eval_buf     evaluation bufferb (work array)
 */
/*----------------------------------------------------------------------------*/

static void
_compute_robin_bc(const cs_mesh_t            *mesh,
                  const cs_boundary_t        *boundaries,
                  const cs_field_t           *f,
                  const cs_equation_param_t  *eqp,
                  const cs_xdef_t            *def,
                  const char                 *description,
                  cs_real_t                   t_eval,
                  int                        *icodcl,
                  double                     *rcodcl1,
                  double                     *rcodcl2,
                  cs_real_t                   eval_buf[])
{
  CS_UNUSED(t_eval);

  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_zone_t *bz = cs_boundary_zone_by_id(def->z_id);
  const cs_lnum_t *elt_ids = bz->elt_ids;
  const cs_lnum_t  n_elts = bz->n_elts;
  const cs_lnum_t  def_dim = def->dim;

  const cs_lnum_t stride = 1 + eqp->dim + eqp->dim*eqp->dim;

  cs_boundary_type_t boundary_type = boundaries->default_type;

  if (boundaries->types != nullptr) {
    assert(boundaries->n_boundaries > 0);
    boundary_type =
      boundaries->types[cs_boundary_id_by_zone_id(boundaries, def->z_id)];
  }

  int bc_code = (boundary_type & CS_BOUNDARY_WALL) ? 5 : 1;
  if (boundary_type & CS_BOUNDARY_ROUGH_WALL)
    bc_code = 6;

  if (stride != def_dim) {
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Boundary condition definition:\n"
                " field %s, zone %s, type %s;\n"
                " dimension %d does not match expected dimension (%d)."),
              __func__, f->name, bz->name, cs_xdef_type_get_name(def->type),
              def_dim, stride);
  }

  switch(def->type) {

  case CS_XDEF_BY_VALUE:  /* direct application (optimization)
                             for constant value */
    {
      const cs_real_t  *constant_val = (cs_real_t *)def->context;

      for (cs_lnum_t coo_id = 0; coo_id < eqp->dim; coo_id++) {

        int        *_icodcl  = icodcl + coo_id*n_b_faces;
        cs_real_t  *_rcodcl1 = rcodcl1 + coo_id*n_b_faces;
        cs_real_t  *_rcodcl2 = rcodcl2 + coo_id*n_b_faces;

        const cs_real_t alpha = constant_val[0];
        const cs_real_t u0 = constant_val[1 + coo_id];

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  elt_id = elt_ids[i];
          _icodcl[elt_id]  = bc_code;
          _rcodcl1[elt_id] = u0;
          _rcodcl2[elt_id] = -alpha;
        }

      }
    }
    break;

  default:
    {
      cs_xdef_eval_at_zone(def,
                           description,
                           t_eval,
                           true,  /* dense */
                           eval_buf);

      for (cs_lnum_t coo_id = 0; coo_id < eqp->dim; coo_id++) {

        int        *_icodcl  = icodcl + coo_id*n_b_faces;
        cs_real_t  *_rcodcl1 = rcodcl1 + coo_id*n_b_faces;
        cs_real_t  *_rcodcl2 = rcodcl2 + coo_id*n_b_faces;

#       pragma omp parallel for if (n_elts > CS_THR_MIN)
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          const cs_lnum_t  elt_id = elt_ids[i];
          _icodcl[elt_id]  = bc_code;
          _rcodcl1[elt_id] = eval_buf[stride*i + 1 + coo_id];
          _rcodcl2[elt_id] = -eval_buf[stride*i];
        }

      }
    }
    break;

  } /* Switch on the type of definition */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute the velocity at boundary faces
 *        using a uniform norm.
 *
 * For the calling function, elt_ids is optional. If non-null, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces * stride
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to cs_gui_boundary_vel_context_t
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_vel_const_uniform_normal(cs_lnum_t         n_elts,
                              const cs_lnum_t  *elt_ids,
                              bool              dense_output,
                              void             *input,
                              cs_real_t        *retval)
{
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *f_n = (const cs_real_3_t *)mq->b_face_u_normal;

  cs_boundary_conditions_open_t *c = (cs_boundary_conditions_open_t *)input;

  cs_real_t u_norm = c->vel_values[3];

  if (c->vel_rescale == CS_BC_VEL_RESCALE_VOLUME_FLOW_RATE)
    u_norm /= c->zone->f_measure;

  if (elt_ids != nullptr) {
    if (dense_output) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = elt_ids[i];
        for (cs_lnum_t k = 0; k < 3; k++)
          retval[i*3 + k] = -f_n[j][k] * u_norm;
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = elt_ids[i];
        for (cs_lnum_t k = 0; k < 3; k++)
          retval[j*3 + k] = -f_n[j][k] * u_norm;
      }
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      for (cs_lnum_t k = 0; k < 3; k++)
        retval[i*3 + k] = -f_n[i][k] * u_norm;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute the velocity at boundary faces
 *        using a vector per face.
 *
 * For the calling function, elt_ids is optional. If non-null, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces * stride
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to cs_gui_boundary_vel_context_t
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_vel_from_buffer(cs_lnum_t         n_elts,
                     const cs_lnum_t  *elt_ids,
                     bool              dense_output,
                     void             *input,
                     cs_real_t        *retval)
{
  cs_boundary_conditions_open_t *c = (cs_boundary_conditions_open_t *)input;

  const cs_real_3_t *vel_b = (const cs_real_3_t *)c->vel_buffer;

  if (elt_ids != nullptr && dense_output == false) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t j = elt_ids[i];
      for (cs_lnum_t k = 0; k < 3; k++)
        retval[j*3 + k] = vel_b[i][k];
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      for (cs_lnum_t k = 0; k < 3; k++)
        retval[i*3 + k] = vel_b[i][k];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_dof_func_t function to compute the velocity at boundary faces
 *        using a vector per face.
 *
 * For the calling function, elt_ids is optional. If non-null, array(s) should
 * be accessed with an indirection. The same indirection can be applied to fill
 * retval if dense_output is set to false.
 * In the current case, retval is allocated to mesh->n_b_faces * stride
 *
 * \param[in]      n_elts        number of elements to consider
 * \param[in]      elt_ids       list of elements ids
 * \param[in]      dense_output  perform an indirection in retval or not
 * \param[in]      input         pointer to cs_gui_boundary_vel_context_t
 * \param[in, out] retval        resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_dof_vel_from_buffer_uniform(cs_lnum_t         n_elts,
                             const cs_lnum_t  *elt_ids,
                             bool              dense_output,
                             void             *input,
                             cs_real_t        *retval)
{
  cs_boundary_conditions_open_t *c = (cs_boundary_conditions_open_t *)input;

  const cs_real_t *vel_b = c->vel_buffer;

  if (elt_ids != nullptr && dense_output == false) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t j = elt_ids[i];
      for (cs_lnum_t k = 0; k < 3; k++)
        retval[j*3 + k] = vel_b[k];
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      for (cs_lnum_t k = 0; k < 3; k++)
        retval[i*3 + k] = vel_b[k];
    }
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant velocity to an inlet.
 *
 * \param[in]  z       pointer to associated zone
 * \param[in]  u_norm  associated constant normal
 */
/*----------------------------------------------------------------------------*/

static void
_clear_inlet_outlet_vel(cs_boundary_conditions_open_t *c)
{
  CS_FREE(c->vel_buffer);

  for (int i = 0; i < 4; i++)
    c->vel_values[i] = 0;

  c->vel_func = nullptr;
  c->vel_func_input = nullptr;
  c->flow_func = nullptr;
  c->flow_func_input = nullptr;
  c->scale_func = nullptr;
  c->scale_func_input = nullptr;

  c->dof_func = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_eval_at_location_t function to compute the velocity at
 *        boundary faces using a constant uniform vector value.
 *
 * For the calling function, elt_ids is optional. If non-null, array(s) should
 * be accessed with an indirection.
 *
 * \param[in]       location_id   base associated mesh location id
 * \param[in]       n_elts        number of elements to consider
 * \param[in]       elt_ids       list of elements ids
 * \param[in]       input         pointer to cs_boundary_conditions_open_t
 * \param[in, out]  vals          resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_vel_profile_constant_uniform(int               location_id,
                              cs_lnum_t         n_elts,
                              const cs_lnum_t  *elt_ids,
                              void             *input,
                              void             *vals_p)
{
  CS_UNUSED(elt_ids);

  assert(   cs_mesh_location_get_type(location_id)
         == CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_3_t *vals = (cs_real_3_t *)vals_p;

  cs_boundary_conditions_open_t  *c
    = (cs_boundary_conditions_open_t *)input;

  assert(c->vel_flags & CS_BC_OPEN_UNIFORM_DIRECTION);

  const cs_real_t *v = c->vel_values;

  for (cs_lnum_t i = 0; i < n_elts; i++) {
    for (cs_lnum_t k = 0; k < 3; k++)
      vals[i][k] = v[k];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_eval_at_location_t function to compute the velocity at
 *        boundary faces using a constant uniform normal value.
 *
 * For the calling function, elt_ids is optional. If non-null, array(s) should
 * be accessed with an indirection.
 *
 * \param[in]       location_id   base associated mesh location id
 * \param[in]       n_elts        number of elements to consider
 * \param[in]       elt_ids       list of elements ids
 * \param[in]       input         pointer to cs_boundary_conditions_open_t
 * \param[in, out]  vals          resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_vel_profile_constant_uniform_normal(int               location_id,
                                     cs_lnum_t         n_elts,
                                     const cs_lnum_t  *elt_ids,
                                     void             *input,
                                     void             *vals_p)
{
  assert(   cs_mesh_location_get_type(location_id)
         == CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_3_t *vals = (cs_real_3_t *)vals_p;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *f_n = (const cs_real_3_t *)mq->b_face_u_normal;

  cs_boundary_conditions_open_t  *c
    = (cs_boundary_conditions_open_t *)input;

  assert(c->vel_flags & CS_BC_OPEN_UNIFORM_QUANTITY);
  assert(c->vel_flags & CS_BC_OPEN_NORMAL_DIRECTION);

  const cs_real_t v
    = (c->vel_rescale == CS_BC_VEL_RESCALE_NONE) ? c->vel_values[3] : 1.0;

  if (elt_ids != nullptr) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      const cs_lnum_t j = elt_ids[i];
      for (cs_lnum_t k = 0; k < 3; k++)
        vals[i][k] = -f_n[j][k] * v;
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      for (cs_lnum_t k = 0; k < 3; k++)
        vals[i][k] = -f_n[i][k] * v;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the mass flow rate for a given boundary zone.
 *
 * This function should be called only after the boundary condition
 * values have been computed.
 *
 * \param[in]   zone  pointer to queried zone
 *
 * \return  zone mass flow rate.
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_compute_mass_flow_rate(const cs_zone_t  *zone)
{
  cs_real_t qm = 0;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  const cs_real_t  *_rcodcl_v[3];
  for (cs_lnum_t coo_id = 0; coo_id < 3; coo_id++)
    _rcodcl_v[coo_id] = CS_F_(vel)->bc_coeffs->rcodcl1 + coo_id*n_b_faces;

  const cs_lnum_t n_elts = zone->n_elts;
  const cs_lnum_t *elt_ids = zone->elt_ids;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *f_n = (const cs_real_3_t *)mq->b_face_normal;

  cs_real_t *sf;
  CS_MALLOC(sf, n_elts, cs_real_t);

  if (CS_F_(rho_b) != nullptr) {
    const cs_real_t *rho_b = CS_F_(rho_b)->val;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t j = elt_ids[i];
      cs_real_t v[3] = {_rcodcl_v[0][j], _rcodcl_v[1][j], _rcodcl_v[1][j]};
      sf[i] = - cs_math_3_dot_product(v, f_n[j]) * rho_b[j];
    }

    qm = cs_sum(n_elts, sf);
    cs_parall_sum(1, CS_REAL_TYPE, &qm);
  }
  else if (CS_F_(rho) != nullptr) {
    const cs_real_t *rho = CS_F_(rho)->val;
    const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t j = elt_ids[i];
      cs_lnum_t k = b_face_cells[j];
      cs_real_t v[3] = {_rcodcl_v[0][j], _rcodcl_v[1][j], _rcodcl_v[1][j]};
      sf[i] = - cs_math_3_dot_product(v, f_n[j]) * rho[k];
    }

    cs_real_t q_m = cs_sum(n_elts, sf);
    cs_parall_sum(1, CS_REAL_TYPE, &q_m);
  }

  CS_FREE(sf);

  return qm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_eval_at_location_t function to scale the velocity at boundary
 *        faces so as to obtain the requested mass flow rate.
 *
 * For the calling function, elt_ids is optional. If non-null, array(s) should
 * be accessed with an indirection.
 *
 * \param[in]       location_id   base associated mesh location id
 * \param[in]       n_elts        number of elements to consider
 * \param[in]       elt_ids       list of elements ids
 * \param[in]       input         pointer to cs_boundary_conditions_open_t
 * \param[in, out]  vals          resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_scale_vel_mass_flow_rate(int               location_id,
                          cs_lnum_t         n_elts,
                          const cs_lnum_t  *elt_ids,
                          void             *input,
                          void             *vals_p)
{
  CS_UNUSED(elt_ids);

  assert(   cs_mesh_location_get_type(location_id)
         == CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_3_t *vals = (cs_real_3_t *)vals_p;

  cs_boundary_conditions_open_t  *c
    = (cs_boundary_conditions_open_t *)input;

  cs_real_t q_m = c->vel_values[3];
  cs_real_t qf = 0.;

  /* Return zero if mass flow is zero, compute ratio otherwise */

  if (fabs(q_m) > 0) {

    const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
    const cs_real_3_t *f_n = (const cs_real_3_t *)mq->b_face_normal;

    cs_real_t *sf;
    CS_MALLOC(sf, n_elts, cs_real_t);

    /* Compressible prescribed inlet with given pressure and temperature */

    if (c->c_pr >= 0 && c->c_tk > 0) {
      for (cs_lnum_t i = 0; i < n_elts; i++) {
        cs_lnum_t j = (elt_ids != nullptr) ? elt_ids[i] : i;
        sf[i] = -   cs_math_3_dot_product(vals[i], f_n[j])
                  * cs_cf_thermo_b_rho_from_pt(j, c->c_pr, c->c_tk);
      }
    }

    /* Regular case */

    else {

      const cs_field_t *rho_b_f = CS_F_(rho_b);

      if (rho_b_f == nullptr) {
        bft_error
          (__FILE__, __LINE__, 0,
           _("%s: no boundary density field for mass flow scaling\n"
             " for zone \"%s\"."),
           __func__, c->zone->name);
      }

      const cs_real_t *rho_b = rho_b_f->val;

      if (elt_ids != nullptr) {
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          cs_lnum_t j = elt_ids[i];
          sf[i] = - cs_math_3_dot_product(vals[i], f_n[j]) * rho_b[j];
        }
      }
      else {
        for (cs_lnum_t i = 0; i < n_elts; i++) {
          sf[i] = - cs_math_3_dot_product(vals[i], f_n[i]) * rho_b[i];
        }
      }

    }

    cs_real_t q_ini = cs_sum(n_elts, sf);
    cs_parall_sum(1, CS_REAL_TYPE, &q_ini);

    CS_FREE(sf);

    if (fabs(q_ini) < 1e-30)
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: cannot normalize for zone \"%s\"\n"
           "  with zero or quasi-zero initial flow rate)."),
         __func__, c->zone->name);

    qf = c->vel_values[3] / q_ini;
  }

  cs_lnum_t n_vals = n_elts*3;

  for (cs_lnum_t i = 0; i < n_vals; i++) {
    c->vel_buffer[i] *= qf;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief cs_eval_at_location_t function to scale the velocity at boundary
 *        faces so as to obtain the requested volume flow rate.
 *
 * For the calling function, elt_ids is optional. If non-null, array(s) should
 * be accessed with an indirection.
 *
 * \param[in]       location_id   base associated mesh location id
 * \param[in]       n_elts        number of elements to consider
 * \param[in]       elt_ids       list of elements ids
 * \param[in]       input         pointer to cs_boundary_conditions_open_t
 * \param[in, out]  vals          resulting value(s). Must be allocated.
 */
/*----------------------------------------------------------------------------*/

static void
_scale_vel_volume_flow_rate(int               location_id,
                            cs_lnum_t         n_elts,
                            const cs_lnum_t  *elt_ids,
                            void             *input,
                            void             *vals_p)
{
  CS_UNUSED(elt_ids);

  assert(   cs_mesh_location_get_type(location_id)
         == CS_MESH_LOCATION_BOUNDARY_FACES);

  cs_real_t *vals = (cs_real_t *)vals_p;

  cs_boundary_conditions_open_t  *c
    = (cs_boundary_conditions_open_t *)input;

  cs_lnum_t r_step = 3;
  if (c->dof_func == _dof_vel_from_buffer_uniform)
    r_step = 0;

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *f_n = (const cs_real_3_t *)mq->b_face_normal;

  cs_real_t *sf;
  CS_MALLOC(sf, n_elts, cs_real_t);

  if (elt_ids != nullptr) {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      cs_lnum_t j = elt_ids[i];
      sf[i] = - cs_math_3_dot_product(vals + i*r_step, f_n[j]);
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      sf[i] = - cs_math_3_dot_product(vals + i*r_step, f_n[i]);
    }
  }

  cs_real_t q_ini = cs_sum(n_elts, sf);
  cs_parall_sum(1, CS_REAL_TYPE, &q_ini);

  CS_FREE(sf);

  if (fabs(q_ini) < 1e-30)
    bft_error
      (__FILE__, __LINE__, 0,
       _("%s: cannot normalize for zone \"%s\"\n"
         "  with zero or quasi-zero initial flow rate)."),
       __func__, c->zone->name);

  cs_real_t qf = c->vel_values[3] / q_ini;

  cs_lnum_t n_vals = n_elts*3;
  if (r_step == 0)
    n_vals = 3;

  for (cs_lnum_t i = 0; i < n_vals; i++) {
    vals[i] *= qf;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update inlet quantities.
 *
 * This updates values associated with a given zone's inlet/outlet.
 *
 * Depending on the complexity of BC velocity definitions, xdef functions will
 * us simpler contexts directly or use computed values from this context to
 * assign BC values.
 *
 * For example, in most cases where rescaling is necessary, values will be
 * updated here first.
 *
 * \param[in]  c  inlet/outlet structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_inlet_outlet(cs_boundary_conditions_open_t  *c)
{
  const cs_time_step_t *ts = cs_glob_time_step;

  /* At the first time step, always update (we do not assume that in case
     of a restarted computation, data is handled using restart data info
     in general);

     At subsequent time steps, update only if the boundary condition might
     depend on time. */

  if (   (ts->nt_cur > ts->nt_prev + 1)
      && (cs_time_control_is_active(&(c->tc), ts) == false))
    return;

  /* Compute flow rate from function if needed */

  if (c->flow_func != nullptr) {
    cs_real_t q[1] = {0};
    c->flow_func(0, 1, nullptr, c->flow_func_input, (void *)q);
    c->vel_values[3] = q[0];
  }

  if (   c->dof_func == nullptr
      || c->dof_func == _dof_vel_const_uniform_normal)
    return;

  const cs_zone_t *z = c->zone;

  if (z->f_measure <= 0. && c->vel_values[3] >= 0. && c->scale_func != nullptr) {

    if (   c->vel_rescale == CS_BC_VEL_RESCALE_MASS_FLOW_RATE
        || c->vel_rescale == CS_BC_VEL_RESCALE_VOLUME_FLOW_RATE) {
      bft_error
        (__FILE__, __LINE__, 0,
         _("%s: cannot normalize requested flow rate (%g) for zone \"%s\"\n"
           "  with zero surface)."),
         __func__, fabs(c->vel_values[3]), z->name);
    }

  }

  /* For uniform direction with constant volume flow, use vel_buffer
     rather than c_vals as dof_func input, so that if the mesh
     is modified (such as with ALE or turbomachinery), the original
     values can be used to recompute the scaling (without progressive
     numerical loss). */

  if (c->dof_func == _dof_vel_from_buffer_uniform) {

    assert(c->scale_func == _scale_vel_volume_flow_rate);

    c->scale_func(z->location_id,
                  z->n_elts,
                  z->elt_ids,
                  c->scale_func_input,
                  c->vel_buffer);

    return;

  }

  assert(c->dof_func == _dof_vel_from_buffer);

  if (   c->vel_buffer == nullptr
      || cs_glob_mesh->time_dep == CS_MESH_TRANSIENT_CONNECT)
    CS_REALLOC(c->vel_buffer, z->n_elts*3, cs_real_t);

  /* Compute velocity buffer */

  if (c->vel_func != nullptr) {
    c->vel_func(z->location_id,
                z->n_elts,
                z->elt_ids,
                c->vel_func_input,
                c->vel_buffer);
  }
  else {
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Missing velocity evaluation function for zone %s."),
              __func__, z->name);
  }

  /* Scale buffer if needed */

  if (c->scale_func != nullptr) {
    c->scale_func(z->location_id,
                  z->n_elts,
                  z->elt_ids,
                  c->scale_func_input,
                  c->vel_buffer);
  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Return pointer to global boundary conditions type array.
 *----------------------------------------------------------------------------*/

int *
cs_f_boundary_conditions_get_bc_type(void)
{
  return _bc_type;
}

/*----------------------------------------------------------------------------
 * Set mapped boundary conditions for a given field and mapping locator.
 *
 * parameters:
 *   field_id        <-- id of field whose boundary conditions are set
 *   locator         <-- associated mapping locator, as returned
 *                       by cs_boundary_conditions_map().
 *   location_type   <-- matching values location (CS_MESH_LOCATION_CELLS or
 *                       CS_MESH_LOCATION_BOUNDARY_FACES)
 *   enforce_balance <-- balance handling option:
 *                         0: values are simply mapped
 *                         1: values are mapped, then multiplied
 *                            by a constant factor so that their
 *                            surface integral on selected faces
 *                            is preserved (relative to the
 *                            input values)
 *                         2: as 1, but with a boundary-defined
 *                            weight, defined by balance_w
 *                         3: as 1, but with a cell-defined
 *                            weight, defined by balance_w
 *   interpolate     <-- interpolation option:
 *                         0: values are simply based on matching
 *                            cell or face center values
 *                         1: values are based on matching cell
 *                            or face center values, corrected
 *                            by gradient interpolation
 *   n_faces         <-- number of selected boundary faces
 *   faces           <-- list of selected boundary faces (1 to n),
 *                       or nullptr if no indirection is needed
 *   balance_w       <-- optional balance weight, or nullptr
 *   nvar            <-- number of variables requiring BC's
 *   rcodcl          <-> boundary condition values
 *----------------------------------------------------------------------------*/

void
cs_f_boundary_conditions_mapped_set(int                        field_id,
                                    ple_locator_t             *locator,
                                    cs_mesh_location_type_t    location_type,
                                    int                        enforce_balance,
                                    int                        interpolate,
                                    cs_lnum_t                  n_faces,
                                    const cs_lnum_t           *faces,
                                    cs_real_t                 *balance_w,
                                    int                        nvar,
                                    cs_real_t                 *rcodcl)
{
  CS_UNUSED(nvar);
  CS_UNUSED(rcodcl);

  cs_lnum_t *_faces = nullptr;

  if (faces != nullptr) {
    CS_MALLOC(_faces, n_faces, cs_lnum_t);
    for (cs_lnum_t i = 0; i < n_faces; i++)
      _faces[i] = faces[i] - 1;
  }

  cs_boundary_conditions_mapped_set(cs_field_by_id(field_id),
                                    locator,
                                    location_type,
                                    enforce_balance,
                                    interpolate,
                                    n_faces,
                                    _faces,
                                    balance_w);

  CS_FREE(_faces);
}

/*----------------------------------------------------------------------------
 * Get pointers for Fortran bindings
 *----------------------------------------------------------------------------*/

void
cs_f_boundary_conditions_get_pointers(int **itypfb,
                                      int **izfppp)
{
  *itypfb = _bc_type;

  *izfppp = cs_glob_bc_pm_info->izfppp;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Semi-private function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Define automatic turbulence values for open boundaries.
 *----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_turb(void)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_b_faces = m->n_b_faces;
  const cs_lnum_t *b_face_cells = m->b_face_cells;

  const cs_real_t *b_rho = CS_F_(rho_b)->val;
  const cs_real_t *c_mu = CS_F_(mu)->val;

  const cs_real_t  *_rcodcl_v[3];
  for (cs_lnum_t coo_id = 0; coo_id < 3; coo_id++)
    _rcodcl_v[coo_id] = CS_F_(vel)->bc_coeffs->rcodcl1 + coo_id*n_b_faces;

  /* Inlet BC's */

  cs_field_t *w_dist = nullptr;
  cs_halo_type_t halo_type = m->halo_type;

  /* For some models, we will need the gradient of the wall distance
     for the shear direction; prepare things here to avoid excess tests
     in loops. */

  if (cs_glob_turb_model->order == CS_TURB_SECOND_ORDER) {
    w_dist = cs_field_by_name_try("wall_distance");
    if (w_dist->bc_coeffs == nullptr)
      w_dist = nullptr;
    else {
      const cs_equation_param_t *eqp = cs_field_get_equation_param(w_dist);
      cs_gradient_type_t  gradient_type;
      cs_gradient_type_by_imrgra(eqp->imrgra,
                                 &gradient_type,
                                 &halo_type);
    }
  }

  /* Loop on open boundaries */

  for (int io_id = 0; io_id < _n_bc_open; io_id++) {
    cs_boundary_conditions_open_t *c = _bc_open[io_id];

    if (c->turb_compute > CS_BC_TURB_NONE) {

      const cs_real_3_t *f_n = (const cs_real_3_t *)mq->b_face_u_normal;
      const cs_zone_t *z = c->zone;
      const cs_real_t dh = c->hyd_diameter;

      const cs_real_t ant = -sqrt(cs_turb_cmu); /* Note: should be
                                                   model dependant */

      /* Now compute BC values at zone's faces */

      for (cs_lnum_t i = 0; i < z->n_elts; i++) {

        cs_lnum_t face_id = z->elt_ids[i];
        cs_lnum_t cell_id = b_face_cells[face_id];

        cs_real_t vel[3] = {_rcodcl_v[0][face_id],
                            _rcodcl_v[1][face_id],
                            _rcodcl_v[2][face_id]};

        if (cs_math_3_dot_product(vel, f_n[face_id]) > 0) {
          cs_turbulence_bc_set_hmg_neumann(face_id);
          continue;
        }

        if (c->turb_compute == CS_BC_TURB_BY_HYDRAULIC_DIAMETER) {
          cs_real_t uref2 = fmax(cs_math_3_square_norm(vel),
                                 cs_math_epzero);

          cs_real_t *vel_dir = nullptr;
          cs_real_t *shear_dir = nullptr;
          cs_real_t _shear_dir[3];

          /* Compute the wall direction only if needed */
          if (w_dist != nullptr) {
            cs_gradient_scalar_cell(m,
                                    mq,
                                    b_face_cells[face_id],
                                    halo_type,
                                    w_dist->bc_coeffs,
                                    w_dist->val,
                                    nullptr,  /* c_weight */
                                    _shear_dir);
            cs_math_3_normalize(_shear_dir, _shear_dir);
            for (int j = 0; j < 3; j++)
              _shear_dir[j] *= ant;

            vel_dir = vel;
            shear_dir = _shear_dir;
          }

          cs_turbulence_bc_set_uninit_inlet_hyd_diam(face_id,
                                                     vel_dir,
                                                     shear_dir,
                                                     uref2, dh,
                                                     b_rho[face_id],
                                                     c_mu[cell_id]);
        }

        else if (c->turb_compute == CS_BC_TURB_BY_TURBULENT_INTENSITY) {
          cs_real_t uref2 = fmax(cs_math_3_square_norm(vel),
                                 cs_math_epzero);

          cs_turbulence_bc_set_uninit_inlet_turb_intensity(face_id,
                                                           uref2,
                                                           c->turb_intensity,
                                                           dh);
        }

      }

    } /* End for inlet */

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find or add an open boundary context for a given zone.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone
 */
/*----------------------------------------------------------------------------*/

cs_boundary_conditions_open_t *
cs_boundary_conditions_open_find_or_add(const  cs_zone_t   *zone)
{
  assert(zone != nullptr);

  for (int i = 0; i < _n_bc_open; i++) {
    if (_bc_open[i]->zone->id == zone->id)
      return _bc_open[i];
  }

  /* If we have not returned yet, we need to add an inlet/outlet structure;

   We assume we do not have a very high number of open boundaries
   in most cases, so we use a simple allocation scheme.
   We could use more complex schemes to limit reallocation, but just
   as for cs_map_t elements, it is important that the adress of
   a given context does not change through reallocation, so as to
   avoid issues where that context is passed as an input to function
   pointers and its stored adress would be invalid. */

  int zone_num_max = 0;
  for (int i = 0; i < _n_bc_open; i++) {
    if (_bc_open[i]->bc_pm_zone_num > zone_num_max)
      zone_num_max = _bc_open[i]->bc_pm_zone_num;
  }

  CS_REALLOC(_bc_open, _n_bc_open+1, cs_boundary_conditions_open_t *);

  cs_boundary_conditions_open_t *c;
  CS_MALLOC(c, 1, cs_boundary_conditions_open_t);

  _bc_open[_n_bc_open] = c;
  _n_bc_open += 1;

  /* Now initialize structure */

  memset(c, 0, sizeof(cs_boundary_conditions_open_t));

  c->zone = zone;

  cs_time_control_init_by_time_step(&(c->tc),
                                    -1,
                                    -1,
                                    1,
                                    true,
                                    false);

  /* Set some default values */

  c->vel_flags = (  CS_BC_OPEN_CONSTANT
                  | CS_BC_OPEN_NORMAL_DIRECTION
                  | CS_BC_OPEN_UNIFORM_QUANTITY);

  c->vel_rescale = CS_BC_VEL_RESCALE_NONE;
  c->turb_compute = CS_BC_TURB_NONE;

  c->bc_pm_zone_num = zone_num_max + 1;

  c->vel_values[0] = 0.;
  c->vel_values[1] = 0.;
  c->vel_values[2] = 0.;
  c->vel_values[3] = 0.;

  c->vel_buffer = nullptr;

  c->hyd_diameter = -1;
  c->turb_intensity = -1;
  c->c_pr = -1;
  c->c_tk = -1;

  c->vel_func = nullptr;
  c->vel_func_input = nullptr;
  c->flow_func = nullptr;
  c->flow_func_input = nullptr;
  c->scale_func = nullptr;
  c->scale_func_input = nullptr;

  c->dof_func = _dof_vel_const_uniform_normal;

  c->model_inlet = nullptr;
  c->model_inlet_del = nullptr;

  return c;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get an open boundary context structure for a given zone.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone, or nullptr if not found.
 */
/*----------------------------------------------------------------------------*/

cs_boundary_conditions_open_t *
cs_boundary_conditions_open_find(const  cs_zone_t   *zone)
{
  if (zone != nullptr) {

    for (int i = 0; i < _n_bc_open; i++) {
      if (_bc_open[i]->zone->id == zone->id)
        return _bc_open[i];
    }

  }

  return nullptr;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Handling of boundary condition definition errors and
 *        associated output.
 *
 * This function checks for errors, and simply returns if no error is
 * encountered. In case of error, it outputs helpful information so as to
 * make it easier to locate the matching faces.
 *
 * For each boundary face, bc_type defines the boundary condition type.
 * As a convention here, zero values correspond to undefined types,
 * positive values to defined types (with no error), and negative values
 * to defined types with inconsistent or incompatible values, the
 * absolute value indicating the original boundary condition type.
 *
 * An optional label may be used if the error is related to another
 * attribute than the boundary type, for appropriate error reporting.
 *
 * \param[in]  bc_flag     array of BC type ids
 * \param[in]  type_name   name of attribute in error, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_error(const int       *bc_flag,
                             const char      *type_name)
{
  const char type_name_default[] = N_("boundary condition type");
  const char *_type_name = type_name_default;

  if (type_name != nullptr)
    _type_name = type_name;

  /* Check for faces with problems;
     bc_flag[] used to determine if there is an error */

  int have_errors
    = cs_flag_check(_("face with boundary condition definition error"),
                    _type_name,
                    _("BC type"),
                    _("Faces with B.C. error"),
                    _("Faces with valid B.C.'s"),
                    CS_MESH_LOCATION_BOUNDARY_FACES,
                    1, /* min_flag */
                    bc_flag);

  if (have_errors)
    bft_error
      (__FILE__, __LINE__, 0,
       _("\nSome boundary condition definitions are incomplete or incorrect.\n"
         "\n  For details, read the end of the calculation log,\n"
         "  or visualize the error postprocessing output."));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Locate shifted boundary face coordinates on possibly filtered
 *        cells or boundary faces for later interpolation.
 *
 * \param[in]  location_type    matching values location (CS_MESH_LOCATION_CELLS
 *                              or CS_MESH_LOCATION_BOUNDARY_FACES)
 * \param[in]  n_location_elts  number of selected location elements
 * \param[in]  n_faces          number of selected boundary faces
 * \param[in]  location_elts    list of selected location elements (0 to n-1),
 *                              or nullptr if no indirection is needed
 * \param[in]  faces            list of selected boundary faces (0 to n-1),
 *                              or nullptr if no indirection is needed
 * \param[in]  coord_shift      array of coordinates shift relative to selected
 *                              boundary faces
 * \param[in]  coord_stride     access stride in coord_shift: 0 for uniform
 *                              shift, 1 for "per face" shift.
 * \param[in]  tolerance        relative tolerance for point location.
 *
 * \return  associated locator structure
 */
/*----------------------------------------------------------------------------*/

ple_locator_t *
cs_boundary_conditions_map(cs_mesh_location_type_t    location_type,
                           cs_lnum_t                  n_location_elts,
                           cs_lnum_t                  n_faces,
                           const cs_lnum_t           *location_elts,
                           const cs_lnum_t           *faces,
                           cs_real_3_t               *coord_shift,
                           int                        coord_stride,
                           double                     tolerance)
{
  ple_locator_t *locator = nullptr;

  /* Build temporary "donor" location  mesh */
  /*----------------------------------------*/

  fvm_nodal_t *nm = nullptr;

  if (location_type == CS_MESH_LOCATION_CELLS)
    nm = cs_mesh_connect_cells_to_nodal(cs_glob_mesh,
                                        "search mesh",
                                        false, /* include_families */
                                        n_location_elts,
                                        location_elts);
  else if (location_type == CS_MESH_LOCATION_BOUNDARY_FACES)
    nm = cs_mesh_connect_faces_to_nodal(cs_glob_mesh,
                                        "search mesh",
                                        false, /* include_families */
                                        0,
                                        n_location_elts,
                                        nullptr,
                                        location_elts);

  /* Now build locator */
  /*-------------------*/

#if defined(PLE_HAVE_MPI)
  locator = ple_locator_create(cs_glob_mpi_comm, cs_glob_n_ranks, 0);
#else
  locator = ple_locator_create();
#endif

  int options[PLE_LOCATOR_N_OPTIONS];
  for (int i = 0; i < PLE_LOCATOR_N_OPTIONS; i++)
    options[i] = 0;
  options[PLE_LOCATOR_NUMBERING] = 0; /* base 0 numbering */

  /* Build location coordinates */

  ple_coord_t *point_coords;

  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *)cs_glob_mesh_quantities->b_face_cog;

  CS_MALLOC(point_coords, n_faces*3, ple_coord_t);
  if (faces != nullptr) {
    for (cs_lnum_t i = 0; i < n_faces; i++) {
      const cs_lnum_t face_id = faces[i];
      for (cs_lnum_t j = 0; j < 3; j++)
        point_coords[i*3 + j] =   b_face_cog[face_id][j]
                                + coord_shift[i*coord_stride][j];
    }
  }
  else {
    for (cs_lnum_t i = 0; i < n_faces; i++) {
      for (cs_lnum_t j = 0; j < 3; j++) {
        point_coords[i*3 + j] =   b_face_cog[i][j]
                                + coord_shift[i*coord_stride][j];
      }
    }
  }

  ple_locator_set_mesh(locator,
                       nm,
                       options,
                       0., /* tolerance_base */
                       tolerance,
                       3, /* dim */
                       n_faces,
                       nullptr,
                       nullptr, /* point_tag */
                       point_coords,
                       nullptr, /* distance */
                       cs_coupling_mesh_extents,
                       cs_coupling_point_in_mesh_p);

  /* Check that location is OK */

  cs_gnum_t loc_count[2];
  loc_count[0] = ple_locator_get_n_exterior(locator);
  loc_count[1] = ple_locator_get_n_exterior(locator);
  cs_parall_counter(loc_count, 2);

  if (loc_count[1] > 0) {
    bft_error
      (__FILE__, __LINE__, 0,
       _("\nIn function %s,\n"
         "  %llu boundary faces (of %llu selected) were not matched to mesh\n"
         "  elements. Check your coordinate shift definitions."),
       __func__,
       (unsigned long long)loc_count[1],
       (unsigned long long)loc_count[0]);
  }

  CS_FREE(point_coords);

  /* Shift from 1-base to 0-based locations */

  ple_locator_shift_locations(locator, -1);

  /* Nodal mesh is not needed anymore */

  nm = fvm_nodal_destroy(nm);

  /* Return initialized locator */

  return locator;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set mapped boundary conditions for a given field and mapping locator.
 *
 * \param[in]       f                field whose boundary conditions are set
 * \param[in]       locator          associated mapping locator, as returned
 *                                   by \ref cs_boundary_conditions_map.
 * \param[in]       location_type    matching values location
 *                                   (CS_MESH_LOCATION_CELLS or
 *                                   CS_MESH_LOCATION_BOUNDARY_FACES)
 * \param[in]       normalize        normalization option:
 *                                     0: values are simply mapped
 *                                     1: values are mapped, then multiplied
 *                                        by a constant factor so that their
 *                                        surface integral on selected faces
 *                                        is preserved (relative to the
 *                                        input values)
 *                                     2: as 1, but with a boundary-defined
 *                                        weight, defined by balance_w
 *                                     3: as 1, but with a cell-defined
 *                                       weight, defined by balance_w
 * \param[in]       interpolate      interpolation option:
 *                                     0: values are simply based on matching
 *                                        cell or face center values
 *                                     1: values are based on matching cell
 *                                        or face center values, corrected
 *                                        by gradient interpolation
 * \param[in]       n_faces          number of selected boundary faces
 * \param[in]       faces            list of selected boundary faces (0 to n-1),
 *                                   or nullptr if no indirection is needed
 * \param[in]       balance_w        optional balance weight, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_mapped_set(const cs_field_t          *f,
                                  ple_locator_t             *locator,
                                  cs_mesh_location_type_t    location_type,
                                  int                        normalize,
                                  int                        interpolate,
                                  cs_lnum_t                  n_faces,
                                  const cs_lnum_t           *faces,
                                  cs_real_t                 *balance_w)
{
  const int dim = f->dim;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  const ple_lnum_t n_dist = ple_locator_get_n_dist_points(locator);
  const ple_lnum_t *dist_loc = ple_locator_get_dist_locations(locator);
  const ple_coord_t *dist_coords = ple_locator_get_dist_coords(locator);

  cs_field_interpolate_t   interpolation_type = CS_FIELD_INTERPOLATE_MEAN;

  cs_real_t inlet_sum_0[9], inlet_sum_1[9];
  cs_real_t *distant_var, *local_var;

  /* Get field's variable id */

  if (f->bc_coeffs == nullptr)
    return;

  assert(f->location_id == CS_MESH_LOCATION_CELLS);

  /* Compute initial balance if applicable */

  if (normalize > 0) {
    assert(dim <= 9);
    _inlet_sum(f,
               cs_glob_mesh,
               cs_glob_mesh_quantities,
               normalize,
               n_faces,
               faces,
               balance_w,
               inlet_sum_0);
  }

  /* Allocate working array */

  CS_MALLOC(distant_var, n_dist*dim, cs_real_t);
  CS_MALLOC(local_var, n_faces*dim, cs_real_t);

  /* Prepare values to send */
  /*------------------------*/

  if (interpolate)
    interpolation_type = CS_FIELD_INTERPOLATE_GRADIENT;

  assert(sizeof(ple_coord_t) == sizeof(cs_real_t));

  if (location_type == CS_MESH_LOCATION_CELLS || interpolate) {
    /* FIXME: we cheat here with the constedness of the field
       for a possible ghost values update, but having a finer control
       of when syncing is required would be preferable */
    cs_field_t *_f = cs_field_by_id(f->id);
    cs_field_interpolate(_f,
                         interpolation_type,
                         n_dist,
                         dist_loc,
                         (const cs_real_3_t *)dist_coords,
                         distant_var);
  }

  else if (location_type == CS_MESH_LOCATION_BOUNDARY_FACES) {

    const cs_lnum_t *b_face_cells = cs_glob_mesh->b_face_cells;
    const cs_field_bc_coeffs_t   *bc_coeffs = f->bc_coeffs;

    /* If boundary condition coefficients are available */

    if (bc_coeffs != nullptr) {

      if (dim == 1) {
        for (cs_lnum_t i = 0; i < n_dist; i++) {
          cs_lnum_t f_id = dist_loc[i];
          cs_lnum_t c_id = b_face_cells[f_id];
          distant_var[i] =   bc_coeffs->a[f_id]
                           + bc_coeffs->b[f_id]*f->val[c_id];
        }
      }
      else {
        for (cs_lnum_t i = 0; i < n_dist; i++) {
          cs_lnum_t f_id = dist_loc[i];
          cs_lnum_t c_id = b_face_cells[f_id];
          for (cs_lnum_t j = 0; j < dim; j++) {
            distant_var[i*dim+j] = bc_coeffs->a[f_id*dim+j];
            for (cs_lnum_t k = 0; k < dim; k++)
              distant_var[i*dim+j] +=  bc_coeffs->b[(f_id*dim+k)*dim + j]
                                      *f->val[c_id*dim+k];
          }
        }
      }

    }

    /* If no boundary condition coefficients are available */

    else {

      for (cs_lnum_t i = 0; i < n_dist; i++) {
        cs_lnum_t f_id = dist_loc[i];
        cs_lnum_t c_id = b_face_cells[f_id];
        for (cs_lnum_t j = 0; j < dim; j++)
          distant_var[i*dim+j] = f->val[c_id*dim+j];
      }

    }

  }

  ple_locator_exchange_point_var(locator,
                                 distant_var,
                                 local_var,
                                 nullptr,               /* faces indirection */
                                 sizeof(cs_real_t),
                                 f->dim,
                                 0);

  /* Now set boundary condition values */

  for (cs_lnum_t j = 0; j < dim; j++) {

    cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1 + j*n_b_faces;

    for (cs_lnum_t i = 0; i < n_faces; i++) {
      const cs_lnum_t f_id = (faces != nullptr) ? faces[i] : i;
      rcodcl1[f_id] = local_var[i*dim + j];
    }

  }

  CS_FREE(local_var);
  CS_FREE(distant_var);

  /* Compute initial balance if applicable */

  if (normalize > 0) {

    _inlet_sum(f,
               cs_glob_mesh,
               cs_glob_mesh_quantities,
               normalize,
               n_faces,
               faces,
               balance_w,
               inlet_sum_1);

    for (cs_lnum_t j = 0; j < dim; j++) {

      const cs_real_t f_mult = (fabs(inlet_sum_1[j]) > 1.e-24) ?
                               inlet_sum_0[j] / inlet_sum_1[j] : 1.;

      cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1 + j*n_b_faces;

      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != nullptr) ? faces[i] : i;
        rcodcl1[f_id] *= f_mult;
      }

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Add location of locate shifted boundary face coordinates on
 *        cells or boundary faces for automatic interpolation.
 *
 * \note
 * This function is currently restricted to mapping of boundary face
 * locations (usually from boundary zones) to cell of boundary face
 * locations, but could be extended to other location types in the future.
 *
 * \param[in]  bc_location_id      id of selected boundary mesh location;
 *                                 currently restricted to subsets of
 *                                 boundary faces (i.e. boundary zone
 *                                 location ids).
 * \param[in]  source_location_id  id of selected location  mesh location
 *                                 (usually CS_MESH_LOCATION_CELLS but can be
 *                                 a more restricted cell or boundary face zone
 *                                 location location id).
 * \param[in]  coord_shift      coordinates shift relative to selected
 *                              boundary faces
 * \param[in]  tolerance        relative tolerance for point location.
 *
 * \return  id of added map
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_conditions_add_map(int         bc_location_id,
                               int         source_location_id,
                               cs_real_t   coord_shift[3],
                               double      tolerance)
{
  int map_id = _n_bc_maps;

  CS_REALLOC(_bc_maps, _n_bc_maps+1, cs_bc_map_t);

  cs_bc_map_t *bc_map = _bc_maps + _n_bc_maps;

  bc_map->bc_location_id = bc_location_id;
  bc_map->source_location_id = source_location_id;

  for (int i = 0; i < 3; i++)
    bc_map->coord_shift[i] = coord_shift[i];

  bc_map->tolerance = tolerance;

  bc_map->locator = nullptr;

  _n_bc_maps += 1;

  return map_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create the legacy boundary conditions zone data arrays
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_create_legacy_zone_data(void)
{
  if (cs_glob_bc_pm_info != nullptr)
    return;

  /* Allocate legacy zone info
     (deprecated, only for specific physics) */

  CS_MALLOC(cs_glob_bc_pm_info, 1, cs_boundary_condition_pm_info_t);
  cs_glob_bc_pm_info->izfppp = nullptr;

  cs_boundary_condition_pm_info_t *bc_pm_info = cs_glob_bc_pm_info;

  bc_pm_info->iautom = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create the legacy boundary conditions face type and face zone arrays.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_create(void)
{
  if (cs_glob_bc_pm_info == nullptr)
    cs_boundary_conditions_create_legacy_zone_data();

  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  /* Get default boundary type (converting "current" to "legacy" codes) */

  const cs_boundary_t  *boundaries = cs_glob_boundaries;
  int default_type = 0;
  if (boundaries->default_type & CS_BOUNDARY_WALL)
    default_type = CS_SMOOTHWALL;
  else if (boundaries->default_type & CS_BOUNDARY_SYMMETRY)
    default_type = CS_SYMMETRY;

  /* boundary conditions type by boundary face */

  CS_MALLOC_HD(_bc_type, n_b_faces, int, cs_alloc_mode);
  for (cs_lnum_t ii = 0; ii < n_b_faces; ii++) {
    _bc_type[ii] = default_type;
  }
  cs_glob_bc_type = _bc_type;

  /* Allocate legacy zone info
     (deprecated, only for specific physics) */

  cs_boundary_condition_pm_info_t *bc_pm_info = cs_glob_bc_pm_info;
  assert(bc_pm_info != nullptr);

  CS_MALLOC(bc_pm_info->izfppp, n_b_faces, int);
  for (cs_lnum_t ii = 0; ii < n_b_faces; ii++) {
    cs_glob_bc_pm_info->izfppp[ii] = 0;
  }

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
    CS_MALLOC(bc_pm_info->iautom, n_b_faces, int);
    for (cs_lnum_t ii = 0; ii < n_b_faces; ii++) {
      bc_pm_info->iautom[ii] = 0;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free the boundary conditions face type and face zone arrays.
 *
 * This also frees boundary condition mappings which may have been defined.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_free(void)
{
  CS_FREE(_bc_type);
  CS_FREE(_bc_pm_face_zone);

  for (int i = 0; i < _n_bc_maps; i++)
    ple_locator_destroy((_bc_maps + i)->locator);

  CS_FREE(_bc_maps);
  _n_bc_maps = 0;

  for (int i = 0; i < _n_bc_open; i++) {
    cs_boundary_conditions_open_t *c = _bc_open[i];
    CS_FREE(c->vel_buffer);
    if (c->model_inlet != nullptr) {
      if (c->model_inlet_del != nullptr)
        c->model_inlet_del(c->model_inlet);
      else
        CS_FREE(c->model_inlet);
    }
    CS_FREE(c);
    _bc_open[i] = nullptr;
  }
  CS_FREE(_bc_open);
  _n_bc_open = 0;

  if (cs_glob_bc_pm_info != nullptr) {
    CS_FREE(cs_glob_bc_pm_info->iautom);
    CS_FREE(cs_glob_bc_pm_info->izfppp);
    CS_FREE(cs_glob_bc_pm_info);
  }

  CS_FREE(_b_head_loss);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Prepare (reset) condition coefficients for all variable fields.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_reset(void)
{
  const cs_lnum_t n_b_faces
    = cs_mesh_location_get_n_elts(CS_MESH_LOCATION_BOUNDARY_FACES)[0];

  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {

    const cs_field_t  *f = cs_field_by_id(f_id);

    if (f->type & CS_FIELD_VARIABLE && f->bc_coeffs != nullptr) {

      int *icodcl  = f->bc_coeffs->icodcl;
      cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1;
      cs_real_t *rcodcl2 = f->bc_coeffs->rcodcl2;
      cs_real_t *rcodcl3 = f->bc_coeffs->rcodcl3;

      /* Initialize all icodcl and rcodcl values to defaults, to
         to avoid issues in some Fortran initialization loops
         which do not do the appropriate checks for variable type. */

      cs_lnum_t n = n_b_faces * (cs_lnum_t)(f->dim);

      for (cs_lnum_t i = 0; i < n; i++) {
        icodcl[i] = 0;
        rcodcl1[i] = cs_math_infinite_r;
        rcodcl2[i] = cs_math_infinite_r;
        rcodcl3[i] = 0;
      }

    }

  }

  if (cs_glob_ale > CS_ALE_NONE) {
    int *ale_bc_type = cs_glob_ale_data->bc_type;

    for (cs_lnum_t i = 0; i < n_b_faces; i++)
      ale_bc_type[i] = 0;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to boundary conditions BC type array.
 */
/*----------------------------------------------------------------------------*/

int *
cs_boundary_conditions_get_bc_type(void)
{
  return _bc_type;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update per variable boundary condition codes.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_compute([[maybe_unused]] int  bc_type[])
{
  /* Initialization */

  const cs_time_step_t *ts = cs_glob_time_step;
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  for (int map_id = 0; map_id < _n_bc_maps; map_id++)
    _update_bc_map(map_id);

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_boundary_t *boundaries = cs_glob_boundaries;
  const cs_real_t t_eval = ts->t_cur;

  /* Buffer for evaluation of analytical values */

  cs_lnum_t  eval_buf_size = 0;
  cs_real_t *eval_buf = nullptr;

  /* Loop on inlets */

  for (int io_id = 0; io_id < _n_bc_open; io_id++)
    _update_inlet_outlet(_bc_open[io_id]);

  /* Loop on fields */

  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    const cs_equation_param_t *eqp = nullptr;
    if (f->type & CS_FIELD_VARIABLE)
      eqp = cs_field_get_equation_param_const(f);

    if (eqp == nullptr)
      continue;

    /* Only handle legacy discretization here */

    if (eqp->space_scheme != CS_SPACE_SCHEME_LEGACY)
      continue;

    /* Settings in eqp may not be well-defined at this stage. The following
       test should be more appropriate to decide if one skips this field or
       not */
    if (f->type & CS_FIELD_CDO)
      continue;

    if (f->bc_coeffs == nullptr)
      continue;

    char description[128];
    snprintf(description, 127, _("boundary condition definitions for \"%s\""),
             f->name);
    description[127] = '\0';

    /* Conversion flag for enthalpy (temperature given);
       not active yet in GUI and XML. */

    int icodcl_m = 1;

    cs_lnum_t n_max_vals = (cs_lnum_t)(f->dim) * n_b_faces;

    if (n_max_vals > eval_buf_size) {
      eval_buf_size = n_max_vals;
      CS_FREE(eval_buf);
      CS_MALLOC(eval_buf, eval_buf_size, cs_real_t);
    }

    /* Loop on boundary conditions */

    for (int bc_id = 0; bc_id < eqp->n_bc_defs; bc_id++) {
      const cs_xdef_t *def = eqp->bc_defs[bc_id];
      cs_param_bc_type_t bc_type_l = (cs_param_bc_type_t)(def->meta);
      switch (bc_type_l) {

      case CS_BC_HMG_DIRICHLET:
        _compute_hmg_dirichlet_bc(mesh,
                                  boundaries,
                                  eqp,
                                  def,
                                  icodcl_m,
                                  f->bc_coeffs->icodcl,
                                  f->bc_coeffs->rcodcl1);
        break;

      case CS_BC_DIRICHLET:
        _compute_dirichlet_bc(mesh,
                              boundaries,
                              f,
                              eqp,
                              def,
                              description,
                              t_eval,
                              icodcl_m,
                              f->bc_coeffs->icodcl,
                              f->bc_coeffs->rcodcl1,
                              eval_buf);
        break;

      case CS_BC_HMG_NEUMANN:
        _compute_hmg_neumann_bc(mesh,
                                def,
                                f->bc_coeffs->icodcl,
                                f->bc_coeffs->rcodcl3);
        break;

      case CS_BC_NEUMANN:
        _compute_neumann_bc(mesh,
                            eqp,
                            def,
                            description,
                            t_eval,
                            f->bc_coeffs->icodcl,
                            f->bc_coeffs->rcodcl3,
                            eval_buf);
        break;

      case CS_BC_ROBIN:
        {
          cs_lnum_t stride = 1 + f->dim + f->dim*f->dim;
          n_max_vals = stride * n_b_faces;
          if (n_max_vals > eval_buf_size) {
            eval_buf_size = n_max_vals;
            CS_FREE(eval_buf);
            CS_MALLOC(eval_buf, eval_buf_size, cs_real_t);
          }

          _compute_robin_bc(mesh,
                            boundaries,
                            f,
                            eqp,
                            def,
                            description,
                            t_eval,
                            f->bc_coeffs->icodcl,
                            f->bc_coeffs->rcodcl1,
                            f->bc_coeffs->rcodcl2,
                            eval_buf);
        }
        break;

      case CS_BC_WALL_MODELLED:
        _compute_dirichlet_bc(mesh,
                              boundaries,
                              f,
                              eqp,
                              def,
                              description,
                              t_eval,
                              icodcl_m,
                              f->bc_coeffs->icodcl,
                              f->bc_coeffs->rcodcl1,
                              eval_buf);
        break;

      default:
        break;
      }
    }

  }

  CS_FREE(eval_buf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic adjustments for boundary condition codes.
 *
 * Currently handles mapped inlets, after the call to \ref stdtcl.
 * As portions of stdtcl are migrated to C, they should be called here,
 * before mapped inlets.
 *
 * \param[in]  bc_type  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_complete(int  bc_type[])
{
  CS_NO_WARN_IF_UNUSED(bc_type);

  /* Initialization */

  const cs_time_step_t *ts = cs_glob_time_step;

  for (int map_id = 0; map_id < _n_bc_maps; map_id++)
    _update_bc_map(map_id);

  /* Loop on fields */

  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);

    if (! (f->type & CS_FIELD_VARIABLE))
      continue;

    const cs_equation_param_t *eqp
      = cs_field_get_equation_param_const(f);

    /* Only handle legacy discretization here */
    if (eqp != nullptr) {
      if (eqp->space_scheme != CS_SPACE_SCHEME_LEGACY)
        continue;
    }

    /* Settings in eqp may not be well-defined at this stage. The following
       test should be more appropriate to decide if one skips this field or
       not */
    if (f->type & CS_FIELD_CDO)
      continue;

    if (f->bc_coeffs == nullptr)
      continue;

    /* Treatment of mapped inlets */

    for (int map_id = 0; map_id < _n_bc_maps; map_id++) {

      cs_bc_map_t *bc_map = _bc_maps + map_id;

      if (bc_map->locator == nullptr || ts->nt_cur <= 1)
        continue;

      int interpolate = 0;
      int normalize = 0;
      if (f == CS_F_(vel))
        normalize = 1;
      else {
        const int keysca = cs_field_key_id("scalar_id");
        if (cs_field_get_key_int(f, keysca) > 0)
          normalize = 1;
      }

      if (f != CS_F_(p)) {
        cs_mesh_location_type_t location_type
          = cs_mesh_location_get_type(bc_map->source_location_id);
        cs_lnum_t n_faces
          = cs_mesh_location_get_n_elts(bc_map->bc_location_id)[0];
        const cs_lnum_t *faces
          = cs_mesh_location_get_elt_ids_try(bc_map->bc_location_id);

        cs_boundary_conditions_mapped_set(f,
                                          bc_map->locator,
                                          location_type,
                                          normalize,
                                          interpolate,
                                          n_faces,
                                          faces,
                                          nullptr);

        if (f == CS_F_(h)) {
          if (f->bc_coeffs != nullptr) {
            int  *icodcl = f->bc_coeffs->icodcl;

            for (cs_lnum_t i = 0; i < n_faces; i++) {
              const cs_lnum_t j = (faces != nullptr) ? faces[i] : i;
              if (icodcl[j] < 0)
                icodcl[j] *= -1;
            }
          }
        }

      }

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to model inlet structure associated with a given
 *        open (inlet/outlet) boundary.
 *
 * The returned pointer is of type void * as it should be cast to the
 * appropriate (model-dependent) type.

 * If no matching parent open boundary has been created yet, it is created.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone
 */
/*----------------------------------------------------------------------------*/

void *
cs_boundary_conditions_get_model_inlet(const cs_zone_t  *zone)
{
  assert(zone != nullptr);

  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(zone);

  return c->model_inlet;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to model inlet structure associated with a given
 *        boundary zone, if associated, or null otherwise.
 *
 * The returned pointer is of type void * as it should be cast to the
 * appropriate (model-dependent) type.

 * If the zone is not and inlet or no matching parent open boundary has
 * been created yet, it is created, null is returned.
 *
 * \param[in]  zone  pointer to associated zone
 *
 * \return: pointer to structure associated with zone, or null
 */
/*----------------------------------------------------------------------------*/

void *
cs_boundary_conditions_get_model_inlet_try(const cs_zone_t  *zone)
{
  assert(zone != nullptr);

  assert(zone != nullptr);

  for (int i = 0; i < _n_bc_open; i++) {
    if (_bc_open[i]->zone->id == zone->id)
      return _bc_open[i];
  }

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*
 * \brief Return legacy zone number related to a given zone, if available.
 *
 * \param[in]  z  pointer to associated zone
 *
 * \return  number associated with legacy zone, or 0 if unavailable.
 */
/*----------------------------------------------------------------------------*/

int
cs_boundary_conditions_get_legacy_zone_num(const  cs_zone_t  *z)
{
  int zone_num = 0;

 cs_boundary_conditions_open_t *c
   = cs_boundary_conditions_open_find(z);
 if (c != nullptr)
   zone_num = c->bc_pm_zone_num;

  return zone_num;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to boundary head losses array.
 *
 * The array is allocated if not previously available.
 *
 * \param[in]  alloc_if_null   do we need to allocate this if not present ?

 * \return  b_head_loss  pointer to boundary head losses array, or nullptr
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_boundary_conditions_get_b_head_loss(bool  alloc_if_null)
{
  if (_b_head_loss == nullptr && alloc_if_null)
    CS_REALLOC(_b_head_loss, cs_glob_mesh->n_b_faces, cs_real_t);

  return _b_head_loss;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign pointer to model inlet structure associated with a given
 *        open (inlet/outlet) boundary.
 *
 * The returned pointer is of type void * as it should be cast to the
 * appropriate (model-dependent) type.

 * If no matching parent open boundary has been created yet, it is created.
 *
 * \param[in]  zone   pointer to associated zone
 * \param[in]  s_ptr  pointer to associated structure
 * \param[in]  s_del  destructor for associated structure, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_assign_model_inlet(const cs_zone_t  *zone,
                                          void             *s_ptr,
                                          void             *s_del)
{
  assert(zone != nullptr);

  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(zone);

  if (c->model_inlet != nullptr && c->model_inlet != s_ptr)
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: model inlet structure already assigned\n"
                " for zone %s."),
              __func__, zone->name);

  c->model_inlet = s_ptr;
  c->model_inlet_del = (cs_destructor_t *)s_del;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Acess the time control structure of an open (inlet/outlet) boundary.
 *
 * This allows modifying that structure, for example updating the inlet
 * velocity values only in a certain time range, and avoiding
 * uneeded recomputations outside that range.
 *
 * \param[in]  zone  pointer to associated zone
 */
/*----------------------------------------------------------------------------*/

cs_time_control_t *
cs_boundary_conditions_open_get_time_control(const  cs_zone_t  *zone)
{
  cs_boundary_conditions_open_t
    *c = cs_boundary_conditions_open_find_or_add(zone);

  return &(c->tc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant velocity to an open (inlet/outlet) boundary.
 *
 * This function may also be used to define the flow direction if called
 * before one of the \c cs_boundary_conditions_open_set_mass_flow_rate
 * or \c cs_boundary_conditions_open_set_volume_flow_rate functions.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  u  associated velocity value
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_velocity_by_value(const cs_zone_t  *z,
                                                  const cs_real_t   u[3])
{
  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(z);

  /* Temporary copy before clear in case u points to c->vel_flags */
  const cs_real_t u_c[] = {u[0], u[1], u[2]};

  _clear_inlet_outlet_vel(c);

  c->vel_rescale = CS_BC_VEL_RESCALE_NONE;

  c->vel_flags = (  CS_BC_OPEN_CONSTANT
                  | CS_BC_OPEN_UNIFORM_QUANTITY
                  | CS_BC_OPEN_UNIFORM_DIRECTION);

  for (int i = 0; i < 3; i++)
    c->vel_values[i] = u_c[i];

  cs_equation_param_t *eqp = _get_equation_param("velocity");
  cs_equation_remove_bc(eqp, c->zone->name);  /* Replace if already set */

  cs_equation_add_bc_by_value(eqp,
                              CS_BC_DIRICHLET,
                              c->zone->name,
                              c->vel_values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant velocity normal to an inlet.
 *
 * \param[in]  z       pointer to associated zone
 * \param[in]  u_norm  associated constant normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_velocity_by_normal_value(const  cs_zone_t  *z,
                                                         cs_real_t     u_norm)
{
  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(z);

  _clear_inlet_outlet_vel(c);

  c->vel_rescale = CS_BC_VEL_RESCALE_NONE;

  c->vel_flags = (  CS_BC_OPEN_CONSTANT
                  | CS_BC_OPEN_NORMAL_DIRECTION
                  | CS_BC_OPEN_UNIFORM_QUANTITY);

  c->vel_values[3] = u_norm;

  c->dof_func = _dof_vel_const_uniform_normal;

  cs_equation_param_t *eqp = _get_equation_param("velocity");
  cs_equation_remove_bc(eqp, c->zone->name);  /* Replace if already set */

  cs_equation_add_bc_by_dof_func(eqp,
                                 CS_BC_DIRICHLET,
                                 z->name,
                                 cs_flag_boundary_face,  // location flag
                                 c->dof_func,
                                 c);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a normal velocity to an inlet using a provided function.
 *
 * Reminder: if the input pointer is non-null, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated velocity vector evaluation function
 *                    at zone faces
 * \param[in]  input  optional function evaluation input, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_velocity_by_func(const  cs_zone_t       *z,
                                                 cs_eval_at_location_t  *func,
                                                 void                   *input)
{
  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(z);

  _clear_inlet_outlet_vel(c);

  c->vel_rescale = CS_BC_VEL_RESCALE_NONE;

  c->vel_flags = 0;

  c->vel_func = func;
  c->vel_func_input = input;

  c->dof_func = _dof_vel_from_buffer;

  cs_equation_param_t *eqp = _get_equation_param("velocity");
  cs_equation_remove_bc(eqp, z->name);  /* Replace if already set */

  cs_equation_add_bc_by_dof_func(eqp,
                                 CS_BC_DIRICHLET,
                                 z->name,
                                 cs_flag_boundary_face,  // location flag
                                 c->dof_func,
                                 c);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the volume flow rate to an inlet or outlet.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * volume flow rate.
 *
 * \param[in]  z  pointer to associated zone
 *
 * \return  volume flow rate associated with open boundary
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_boundary_conditions_open_get_mass_flow_rate(const  cs_zone_t  *z)
{
  cs_real_t qm = 0;

  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find(z);

  /* If mass flow rate is prescribed, assume it has already been computed */

  if (c != nullptr) {
    if (c->vel_rescale == CS_BC_VEL_RESCALE_MASS_FLOW_RATE)
      qm = c->vel_values[3];
    else
      qm = _compute_mass_flow_rate(z);
  }
  else
    qm = _compute_mass_flow_rate(z);

  return qm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant mass flow rate to an inlet or outlet.
 *
 * By default, the flow direction is considered normal to the boundary.
 * The flow direction may be specified by calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func
 * for the appropriate zone before calling this function.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * mass flow rate.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  q  associated constant mass flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_mass_flow_rate_by_value(const  cs_zone_t  *z,
                                                        cs_real_t          q)
{
  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(z);

  c->vel_rescale = CS_BC_VEL_RESCALE_MASS_FLOW_RATE;

  /* Velocity based on constant flow rate can be constant only if density is;
     if we are not sure of this, remove the flag if present (or add test on
     variable density).
     The other flags (for direction and uniform quantities)
     should remain valid */

  c->vel_flags = _unset_flag(c->vel_flags, CS_BC_OPEN_CONSTANT);

  c->flow_func = nullptr;
  c->flow_func_input = nullptr;

  c->vel_values[3] = q;

  c->dof_func = _dof_vel_from_buffer;

  if (c->vel_func == nullptr) {
    if (c->vel_flags & CS_BC_OPEN_UNIFORM_DIRECTION) {
      c->vel_func = _vel_profile_constant_uniform;
      c->vel_func_input = c;
    }
    else if (   (c->vel_flags & CS_BC_OPEN_UNIFORM_QUANTITY)
             && (c->vel_flags & CS_BC_OPEN_NORMAL_DIRECTION)) {
      c->vel_func = _vel_profile_constant_uniform_normal;
      c->vel_func_input = c;
    }
  }

  /* For most models, apply scaling to convert from mass flow
     to velocity; For some combustion models, scaling is done later. */

  c->scale_func = _scale_vel_mass_flow_rate;
  c->scale_func_input = c;

  if (cs_glob_physical_model_flag[CS_COMBUSTION_LW] >= 0) {
    c->scale_func = nullptr;
    c->scale_func_input = nullptr;
  }

  /* Set equation parameters */

  cs_equation_param_t *eqp = _get_equation_param("velocity");
  cs_equation_remove_bc(eqp, z->name);  /* Replace if already set */

  cs_equation_add_bc_by_dof_func(eqp,
                                 CS_BC_DIRICHLET,
                                 z->name,
                                 cs_flag_boundary_face,  // location flag
                                 c->dof_func,
                                 c);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a mass flow rate to an inlet or outlet  based on
 *        provided function.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * \ref cs_boundary_conditions_open_set_velocity_by_normal_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * mass flow rate.
 *
 * Since the flow rate is a global value, the provided function should
 * be associated with the CS_MESH_LOCATION_NONE location.
 *
 * Note also that during updates, this function will be called before
 * the velocity vector update, so in complex cases where flow rate computation
 * would require feedback from the velocity at this boundary, the user
 * must be aware that values from the previous time step or update will
 * be used, handle this in another manner.
 *
 * Reminder: if the input pointer is non-null, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated scalar (mass flow rate) evaluation function
 * \param[in]  input  optional function evaluation input, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_mass_flow_rate_by_func
  (const  cs_zone_t       *z,
   cs_eval_at_location_t  *func,
   void                   *input)
{
  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(z);

  c->vel_rescale = CS_BC_VEL_RESCALE_MASS_FLOW_RATE;

  c->vel_flags = _unset_flag(c->vel_flags, CS_BC_OPEN_CONSTANT);

  c->flow_func = func;
  c->flow_func_input = input;

  c->dof_func = _dof_vel_from_buffer;

  if (c->vel_func == nullptr) {
    if (c->vel_flags & CS_BC_OPEN_UNIFORM_DIRECTION) {
      c->vel_func = _vel_profile_constant_uniform;
      c->vel_func_input = c;
    }
    else if (   (c->vel_flags & CS_BC_OPEN_UNIFORM_QUANTITY)
             && (c->vel_flags & CS_BC_OPEN_NORMAL_DIRECTION)) {
      c->vel_func = _vel_profile_constant_uniform_normal;
      c->vel_func_input = c;
    }
  }

  /* For most models, apply scaling to convert from mass flow
     to velocity; For combustion models, scaling is done later. */

  c->scale_func = _scale_vel_mass_flow_rate;
  c->scale_func_input = c;

  for (int i = CS_COMBUSTION_SLFM; i < CS_COMBUSTION_COAL; i++) {
    if (cs_glob_physical_model_flag[i] >= 0) {
      c->scale_func = nullptr;
      c->scale_func_input = nullptr;
    }
  }

  /* Set equation parameters */

  cs_equation_param_t *eqp = _get_equation_param("velocity");
  cs_equation_remove_bc(eqp, z->name);  /* Replace if already set */

  cs_equation_add_bc_by_dof_func(eqp,
                                 CS_BC_DIRICHLET,
                                 z->name,
                                 cs_flag_boundary_face,  // location flag
                                 c->dof_func,
                                 c);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a constant volume flow rate to an inlet or outlet.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * volume flow rate.
 *
 * \param[in]  z  pointer to associated zone
 * \param[in]  q  associated constant volume flow rate
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_volume_flow_rate_by_value(const  cs_zone_t  *z,
                                                          cs_real_t          q)
{
  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(z);

  c->vel_rescale = CS_BC_VEL_RESCALE_VOLUME_FLOW_RATE;

  c->flow_func = nullptr;
  c->flow_func_input = nullptr;

  c->vel_values[3] = q;

  if (   c->vel_flags & CS_BC_OPEN_NORMAL_DIRECTION
      && c->vel_flags & CS_BC_OPEN_UNIFORM_QUANTITY) {
    if (c->dof_func == nullptr)
      c->dof_func = _dof_vel_const_uniform_normal;
  }

  else if (   c->vel_flags & CS_BC_OPEN_UNIFORM_QUANTITY
           && c->vel_flags & CS_BC_OPEN_UNIFORM_DIRECTION) {
    if (c->dof_func == nullptr) {
      CS_REALLOC(c->vel_buffer, 3, cs_real_t);
      for (int i = 0; i < 3; i++)
        c->vel_buffer[i] = c->vel_values[i];

      c->dof_func = _dof_vel_from_buffer_uniform;
    }
  }

  if (c->dof_func == nullptr)
    c->dof_func = _dof_vel_from_buffer;

  c->vel_flags = _unset_flag(c->vel_flags, CS_BC_OPEN_CONSTANT);

  c->scale_func = _scale_vel_volume_flow_rate;
  c->scale_func_input = c;

  cs_equation_param_t *eqp = _get_equation_param("velocity");
  cs_equation_remove_bc(eqp, z->name);  /* Replace if already set */

  cs_equation_add_bc_by_dof_func(eqp,
                                 CS_BC_DIRICHLET,
                                 z->name,
                                 cs_flag_boundary_face,  // location flag
                                 c->dof_func,
                                 c);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Assign a volume flow rate to an inlet or outlet based on
 *        provided function.
 *
 * The flow direction may be specified by also calling
 * \ref cs_boundary_conditions_open_set_velocity_by_value,
 * \ref cs_boundary_conditions_open_set_velocity_by_normal_value,
 * or \ref cs_boundary_conditions_open_set_velocity_by_func.
 * In that case, the velocity vector is rescaled so as to obtain the required
 * volume flow rate.
 *
 * Since the flow rate is a global value, the provided function should
 * be associated with the CS_MESH_LOCATION_NONE location.
 *
 * Note also that during updates, this function will be called before
 * the velocity vector update, so in complex cases where flow rate computation
 * would require feedback from the velocity at this boundary, the user
 * must be aware that values from the previous time step or update will
 * be used, handle this in another manner.
 *
 * Reminder: if the input pointer is non-null, it must point to valid data
 * when the selection function is called, so either:
 * - that value or structure should not be temporary (i.e. local);
 * - when a single integer identifier is needed, the input pointer can be
 *   set to that value instead of an actual address;
 *
 * \param[in]  z      pointer to associated zone
 * \param[in]  func   associated scalar (volume flow rate) evaluation function
 * \param[in]  input  optional function evaluation input, or nullptr
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_open_set_volume_flow_rate_by_func
  (const  cs_zone_t       *z,
   cs_eval_at_location_t  *func,
   void                   *input)
{
  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(z);

  c->vel_rescale = CS_BC_VEL_RESCALE_VOLUME_FLOW_RATE;

  c->flow_func = func;
  c->flow_func_input = input;

  if (   c->vel_flags & CS_BC_OPEN_NORMAL_DIRECTION
      && c->vel_flags & CS_BC_OPEN_UNIFORM_QUANTITY) {
    if (c->dof_func == nullptr)
      c->dof_func = _dof_vel_const_uniform_normal;
  }

  else if (   c->vel_flags & CS_BC_OPEN_UNIFORM_QUANTITY
           && c->vel_flags & CS_BC_OPEN_UNIFORM_DIRECTION) {
    if (c->dof_func == nullptr) {
      CS_REALLOC(c->vel_buffer, 3, cs_real_t);
      for (int i = 0; i < 3; i++)
        c->vel_buffer[i] = c->vel_values[i];

      c->dof_func = _dof_vel_from_buffer_uniform;
    }
  }

  if (c->dof_func == nullptr)
    c->dof_func = _dof_vel_from_buffer;

  c->vel_flags = _unset_flag(c->vel_flags, CS_BC_OPEN_CONSTANT);

  c->scale_func = _scale_vel_volume_flow_rate;
  c->scale_func_input = c;

  cs_equation_param_t *eqp = _get_equation_param("velocity");
  cs_equation_remove_bc(eqp, z->name);  /* Replace if already set */

  cs_equation_add_bc_by_dof_func(eqp,
                                 CS_BC_DIRICHLET,
                                 z->name,
                                 cs_flag_boundary_face,  // location flag
                                 c->dof_func,
                                 c);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Base the inlet turbulence values on a a circular duct with smooth
 *        wall (see ref cs_turbulence_bc_ke_hyd_diam).
 *
 * \param[in]  zone  pointer to associated zone
 * \param[in]  hd    associated hydraulic diameter
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_inlet_set_turbulence_hyd_diam(const  cs_zone_t  *zone,
                                                     cs_real_t          hd)
{
  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(zone);

  c->turb_compute = CS_BC_TURB_BY_HYDRAULIC_DIAMETER;
  c->hyd_diameter = hd;
  c->turb_intensity = -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Base the inlet turbulence values on a a circular duct with smooth
 *        wall (see ref cs_turbulence_bc_ke_hyd_diam).
 *
 * \param[in]  zone  pointer to associated zone
 * \param[in]  ti    associated turbulence intensity
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_inlet_set_turbulence_intensity(const  cs_zone_t  *zone,
                                                      cs_real_t          ti)
{
  cs_boundary_conditions_open_t *c
    = cs_boundary_conditions_open_find_or_add(zone);

  c->turb_compute = CS_BC_TURB_BY_TURBULENT_INTENSITY;
  c->turb_intensity = ti;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
