/*============================================================================
 * Boundary condition handling.
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

#include <ple_locator.h>

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_ale.h"
#include "cs_base.h"
#include "cs_boundary.h"
#include "cs_coupling.h"
#include "cs_gradient.h"
#include "cs_gui_util.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_operator.h"
#include "cs_field_pointer.h"
#include "cs_flag_check.h"
#include "cs_halo.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_post.h"
#include "cs_xdef_eval_at_zone.h"
#include "fvm_nodal.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_boundary_conditions.h"

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

typedef struct {

  int             bc_location_id;      /* location id of boundary zone */
  int             source_location_id;  /* location id of source elements */
  cs_real_t       coord_shift[3];      /* coordinates shift relative to
                                          selected boundary faces */
  double          tolerance;           /* search tolerance */

  ple_locator_t  *locator;             /* associated locator */

} cs_bc_map_t;

/*----------------------------------------------------------------------------
 * Local Structure Definitions
 *----------------------------------------------------------------------------*/

/*============================================================================
 *  Global variables
 *============================================================================*/

static int *_bc_type;

static int *_bc_pm_face_zone;

static int _n_bc_maps = 0;
static cs_bc_map_t *_bc_maps = NULL;

const int *cs_glob_bc_type = NULL;

cs_boundary_condition_pm_info_t  *cs_glob_bc_pm_info = NULL;

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
                                      int  **izfppp,
                                      int  **itrifb);

void
cs_f_boundary_conditions_get_ppincl_pointers(int     **iqimp,
                                             int     **icalke,
                                             double  **xintur,
                                             double  **dh);

void
cs_f_boundary_conditions_get_coincl_pointers(int     **ientfu,
                                             int     **ientox,
                                             int     **ientgb,
                                             int     **ientgf,
                                             double  **tkent,
                                             double  **fment,
                                             double  **qimp);

void
cs_f_boundary_conditions_get_cpincl_pointers(int             **ientat,
                                             int             **ientcp,
                                             cs_real_t       **qimpat,
                                             cs_real_t       **timpat,
                                             cs_real_5_t     **qimpcp,
                                             cs_real_5_t     **timpcp,
                                             cs_real_5_20_t  **distch,
                                             int             **inmoxy);

void
cs_f_boundary_conditions_get_atincl_pointers(int **iprofm,
                                             int **iautom);

/*============================================================================
 * Private function definitions
 *============================================================================*/

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

  if (bc_map->locator != NULL)
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
 *                       or NULL if no indirection is needed
 *   balance_w       <-- optional balance weight, or NULL
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
  const cs_real_t *f_surf = mq->b_f_face_surf;

  /* Get field's variable id */

  assert(dim <= 9);

  for (cs_lnum_t j = 0; j < dim; j++) {

    cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1 + j*n_b_faces;

    inlet_sum[j] = 0.;

    if (enforce_balance == 1) {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
        inlet_sum[j] += rcodcl1[f_id]*f_surf[f_id];
      }
    }
    else if (enforce_balance == 2) {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
        inlet_sum[j] += rcodcl1[f_id]*f_surf[f_id]*balance_w[f_id];
      }
    }
    else if (enforce_balance == 3) {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
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

  if (boundaries->types != NULL) {
    assert(boundaries->n_boundaries > 0);
    boundary_type =
      boundaries->types[cs_boundary_id_by_zone_id(boundaries, def->z_id)];
  }

  int bc_type = (boundary_type & CS_BOUNDARY_WALL) ? 5 : 1;
  if (boundary_type & CS_BOUNDARY_ROUGH_WALL)
    bc_type = 6;

  bc_type *= icodcl_m;

  for (cs_lnum_t coo_id = 0; coo_id < def->dim; coo_id++) {

    int        *_icodcl  = icodcl + coo_id*n_b_faces;
    cs_real_t  *_rcodcl1 = rcodcl1 + coo_id*n_b_faces;

#   pragma omp parallel for if (bz->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < bz->n_elts; i++) {
      const cs_lnum_t  elt_id = (face_ids == NULL) ? i : face_ids[i];
      _icodcl[elt_id]  = bc_type;
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

  if (boundaries->types != NULL) {
    assert(boundaries->n_boundaries > 0);
    boundary_type =
      boundaries->types[cs_boundary_id_by_zone_id(boundaries, def->z_id)];
  }

  int bc_type = (boundary_type & CS_BOUNDARY_WALL) ? 5 : 1;
  if (boundary_type & CS_BOUNDARY_ROUGH_WALL)
    bc_type = 6;

  bc_type *= icodcl_m;

  /* Special case for mesh velocity (ALE) wall BC type. */

  if (f->dim == 3) {
    if (   (strcmp(f->name, "mesh_velocity") == 0)
        && (boundary_type & CS_BOUNDARY_WALL))
      bc_type = 1;
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
          _icodcl[elt_id]  = bc_type;
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
          _icodcl[elt_id]  = bc_type;
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

  const int bc_type = 3;

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
          _icodcl[elt_id]  = bc_type;
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
          _icodcl[elt_id]  = bc_type;
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

  if (boundaries->types != NULL) {
    assert(boundaries->n_boundaries > 0);
    boundary_type =
      boundaries->types[cs_boundary_id_by_zone_id(boundaries, def->z_id)];
  }

  int bc_type = (boundary_type & CS_BOUNDARY_WALL) ? 5 : 1;
  if (boundary_type & CS_BOUNDARY_ROUGH_WALL)
    bc_type = 6;

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
          _icodcl[elt_id]  = bc_type;
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
          _icodcl[elt_id]  = bc_type;
          _rcodcl1[elt_id] = eval_buf[stride*i + 1 + coo_id];
          _rcodcl2[elt_id] = -eval_buf[stride*i];
        }

      }
    }
    break;

  } /* Switch on the type of definition */
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
 *                       or NULL if no indirection is needed
 *   balance_w       <-- optional balance weight, or NULL
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

  cs_lnum_t *_faces = NULL;

  if (faces != NULL) {
    BFT_MALLOC(_faces, n_faces, cs_lnum_t);
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

  BFT_FREE(_faces);
}

/*----------------------------------------------------------------------------
 * Get pointers for Fortran bindings
 *----------------------------------------------------------------------------*/

void
cs_f_boundary_conditions_get_pointers(int **itypfb,
                                      int **izfppp,
                                      int **itrifb)
{
  *itypfb = _bc_type;

  *izfppp = cs_glob_bc_pm_info->izfppp;

  *itrifb = cs_glob_bc_pm_info->itrifb;
}

void
cs_f_boundary_conditions_get_ppincl_pointers(int     **iqimp,
                                             int     **icalke,
                                             double  **xintur,
                                             double  **dh)
{
  /* Shift by 1 to compensate for Fortran 1-based access */

  *iqimp  = cs_glob_bc_pm_info->iqimp  + 1;
  *icalke = cs_glob_bc_pm_info->icalke + 1;
  *xintur = cs_glob_bc_pm_info->xintur + 1;
  *dh     = cs_glob_bc_pm_info->dh     + 1;
}

void
cs_f_boundary_conditions_get_coincl_pointers(int     **ientfu,
                                             int     **ientox,
                                             int     **ientgb,
                                             int     **ientgf,
                                             double  **tkent,
                                             double  **fment,
                                             double  **qimp)
{
  /* Shift 1d-arrays by 1 to compensate for Fortran 1-based access */

  *ientfu = cs_glob_bc_pm_info->ientfu + 1;
  *ientox = cs_glob_bc_pm_info->ientox + 1;
  *ientgb = cs_glob_bc_pm_info->ientgb + 1;
  *ientgf = cs_glob_bc_pm_info->ientgf + 1;
  *tkent  = cs_glob_bc_pm_info->tkent  + 1;
  *fment  = cs_glob_bc_pm_info->fment  + 1;
  *qimp   = cs_glob_bc_pm_info->qimp   + 1;
}

void
cs_f_boundary_conditions_get_cpincl_pointers(int             **ientat,
                                             int             **ientcp,
                                             cs_real_t       **qimpat,
                                             cs_real_t       **timpat,
                                             cs_real_5_t     **qimpcp,
                                             cs_real_5_t     **timpcp,
                                             cs_real_5_20_t  **distch,
                                             int             **inmoxy)
{
  /* Shift 1d-arrays by 1 to compensate for Fortran 1-based access */

  *qimpat = cs_glob_bc_pm_info->qimp + 1;

  if (   cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] > -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] > -1) {
    *inmoxy = cs_glob_bc_pm_info->inmoxy + 1;
    *ientat = cs_glob_bc_pm_info->ientat + 1;
    *ientcp = cs_glob_bc_pm_info->ientcp + 1;
    *timpat = cs_glob_bc_pm_info->timpat + 1;
    *qimpcp = cs_glob_bc_pm_info->qimpcp + 1;
    *timpcp = cs_glob_bc_pm_info->timpcp + 1;
    *distch = cs_glob_bc_pm_info->distch + 1;
  }
  else {
    *inmoxy = NULL;
    *ientat = NULL;
    *ientcp = NULL;
    *timpat = NULL;
    *qimpcp = NULL;
    *timpcp = NULL;
    *distch = NULL;
  }
}

void
cs_f_boundary_conditions_get_atincl_pointers(int **iprofm,
                                             int **iautom)
{
  /* Shift 1d-arrays by 1 to compensate for Fortran 1-based access */

  *iprofm = cs_glob_bc_pm_info->iprofm + 1;

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1)
    *iautom = cs_glob_bc_pm_info->iautom;
  else
    *iautom = NULL;
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
 * \param[in]  type_name   name of attribute in error, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_error(const int       *bc_flag,
                             const char      *type_name)
{
  const char type_name_default[] = N_("boundary condition type");
  const char *_type_name = type_name_default;

  if (type_name != NULL)
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
 *                              or NULL if no indirection is needed
 * \param[in]  faces            list of selected boundary faces (0 to n-1),
 *                              or NULL if no indirection is needed
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
  ple_locator_t *locator = NULL;

  /* Build temporary "donor" location  mesh */
  /*----------------------------------------*/

  fvm_nodal_t *nm = NULL;

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
                                        NULL,
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
    = (const cs_real_3_t *restrict)cs_glob_mesh_quantities->b_face_cog;

  BFT_MALLOC(point_coords, n_faces*3, ple_coord_t);
  if (faces != NULL) {
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
                       NULL,
                       NULL, /* point_tag */
                       point_coords,
                       NULL, /* distance */
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

  BFT_FREE(point_coords);

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
 *                                   or NULL if no indirection is needed
 * \param[in]       balance_w        optional balance weight, or NULL
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

  if (f->bc_coeffs == NULL)
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

  BFT_MALLOC(distant_var, n_dist*dim, cs_real_t);
  BFT_MALLOC(local_var, n_faces*dim, cs_real_t);

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

    if (bc_coeffs != NULL) {

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
                                 NULL,               /* faces indirection */
                                 sizeof(cs_real_t),
                                 f->dim,
                                 0);

  /* Now set boundary condition values */

  for (cs_lnum_t j = 0; j < dim; j++) {

    cs_real_t *rcodcl1 = f->bc_coeffs->rcodcl1 + j*n_b_faces;

    for (cs_lnum_t i = 0; i < n_faces; i++) {
      const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
      rcodcl1[f_id] = local_var[i*dim + j];
    }

  }

  BFT_FREE(local_var);
  BFT_FREE(distant_var);

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
        const cs_lnum_t f_id = (faces != NULL) ? faces[i] : i;
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

  BFT_REALLOC(_bc_maps, _n_bc_maps+1, cs_bc_map_t);

  cs_bc_map_t *bc_map = _bc_maps + _n_bc_maps;

  bc_map->bc_location_id = bc_location_id;
  bc_map->source_location_id = source_location_id;

  for (int i = 0; i < 3; i++)
    bc_map->coord_shift[i] = coord_shift[i];

  bc_map->tolerance = tolerance;

  bc_map->locator = NULL;

  _n_bc_maps += 1;

  return map_id;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create the legacy boundary conditions face type and face zone arrays.
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_create(void)
{
  const cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;

  /* Get default boundary type (converting "current" to "legacy" codes) */

  const cs_boundary_t  *boundaries = cs_glob_boundaries;
  int default_type = 0;
  if (boundaries->default_type & CS_BOUNDARY_WALL)
    default_type = CS_SMOOTHWALL;
  else if (boundaries->default_type & CS_BOUNDARY_SYMMETRY)
    default_type = CS_SYMMETRY;

  /* boundary conditions type by boundary face */

  BFT_MALLOC(_bc_type, n_b_faces, int);
  for (cs_lnum_t ii = 0; ii < n_b_faces; ii++) {
    _bc_type[ii] = default_type;
  }
  cs_glob_bc_type = _bc_type;

  /* Allocate legacy zone info
     (deprecated, only for specific physics) */

  BFT_MALLOC(cs_glob_bc_pm_info, 1, cs_boundary_condition_pm_info_t);
  BFT_MALLOC(cs_glob_bc_pm_info->izfppp, n_b_faces, int);
  BFT_MALLOC(cs_glob_bc_pm_info->itrifb, n_b_faces, int);
  for (cs_lnum_t ii = 0; ii < n_b_faces; ii++) {
    cs_glob_bc_pm_info->izfppp[ii] = 0;
    cs_glob_bc_pm_info->itrifb[ii] = 0;
  }

  cs_boundary_condition_pm_info_t *bc_pm_info = cs_glob_bc_pm_info;
  for (int i = 0; i < CS_MAX_BC_PM_ZONE_NUM+1; i++) {
    bc_pm_info->iqimp[i]  = 0;
    bc_pm_info->icalke[i] = 0;
    bc_pm_info->qimp[i]   = 0.;
    bc_pm_info->dh[i]     = 0.;
    bc_pm_info->xintur[i] = 0.;

    //gas combustion
    bc_pm_info->ientfu[i] = 0;
    bc_pm_info->ientox[i] = 0;
    bc_pm_info->ientgb[i] = 0;
    bc_pm_info->ientgf[i] = 0;
    bc_pm_info->tkent[i]  = 0.;
    bc_pm_info->fment[i]  = 0.;

    //atmo
    bc_pm_info->iprofm[i] = 0;
  }

  bc_pm_info->inmoxy = NULL;
  bc_pm_info->ientat = NULL;
  bc_pm_info->ientcp = NULL;
  bc_pm_info->timpat = NULL;
  bc_pm_info->qimpcp = NULL;
  bc_pm_info->timpcp = NULL;
  bc_pm_info->distch = NULL;
  bc_pm_info->iautom = NULL;

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1) {
    BFT_MALLOC(bc_pm_info->iautom, n_b_faces, int);
    for (cs_lnum_t ii = 0; ii < n_b_faces; ii++) {
      bc_pm_info->iautom[ii] = 0;
    }
  }

  /* Arrays present only for coal combustion */

  if (   cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] > -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] > -1) {

    BFT_REALLOC(bc_pm_info->inmoxy, CS_MAX_BC_PM_ZONE_NUM+1, cs_lnum_t);
    BFT_REALLOC(bc_pm_info->ientat, CS_MAX_BC_PM_ZONE_NUM+1, cs_lnum_t);
    BFT_REALLOC(bc_pm_info->ientcp, CS_MAX_BC_PM_ZONE_NUM+1, cs_lnum_t);
    BFT_REALLOC(bc_pm_info->timpat, CS_MAX_BC_PM_ZONE_NUM+1, cs_real_t);
    BFT_REALLOC(bc_pm_info->qimpcp, CS_MAX_BC_PM_ZONE_NUM+1, cs_real_5_t);
    BFT_REALLOC(bc_pm_info->timpcp, CS_MAX_BC_PM_ZONE_NUM+1, cs_real_5_t);
    BFT_REALLOC(bc_pm_info->distch, CS_MAX_BC_PM_ZONE_NUM+1, cs_real_5_20_t);

    for (int i = 0; i < CS_MAX_BC_PM_ZONE_NUM+1; i++) {
      bc_pm_info->inmoxy[i] = 0;
      bc_pm_info->ientat[i] = 0;
      bc_pm_info->ientcp[i] = 0;
      bc_pm_info->timpat[i] = 0;
      for (int j = 0; j < 5; j++) {
        bc_pm_info->qimpcp[i][j] = 0;
        bc_pm_info->timpcp[i][j] = 0;
        for (int k = 0; k < 20; k++)
          bc_pm_info->distch[i][j][k] = 0;
      }
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
  BFT_FREE(_bc_type);
  BFT_FREE(_bc_pm_face_zone);

  for (int i = 0; i < _n_bc_maps; i++)
    ple_locator_destroy((_bc_maps + i)->locator);

  BFT_FREE(_bc_maps);
  _n_bc_maps = 0;

  if (   cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] > -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_COAL] > -1
      || cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] > -1) {

    BFT_FREE(cs_glob_bc_pm_info->inmoxy);
    BFT_FREE(cs_glob_bc_pm_info->ientat);
    BFT_FREE(cs_glob_bc_pm_info->ientcp);
    BFT_FREE(cs_glob_bc_pm_info->timpat);
    BFT_FREE(cs_glob_bc_pm_info->qimpcp);
    BFT_FREE(cs_glob_bc_pm_info->timpcp);
    BFT_FREE(cs_glob_bc_pm_info->distch);
  }

  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] > -1)
    BFT_FREE(cs_glob_bc_pm_info->iautom);

  BFT_FREE(cs_glob_bc_pm_info->izfppp);
  BFT_FREE(cs_glob_bc_pm_info->itrifb);
  BFT_FREE(cs_glob_bc_pm_info);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set convective oulet boundary condition for a scalar.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimp          Flux value to impose
 * \param[in]     cfl           Local Courant number used to convect
 * \param[in]     hint          Internal exchange coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_convective_outlet_scalar(cs_real_t *coefa ,
                                                    cs_real_t *cofaf,
                                                    cs_real_t *coefb,
                                                    cs_real_t *cofbf,
                                                    cs_real_t  pimp,
                                                    cs_real_t  cfl,
                                                    cs_real_t  hint)
{
  /* Gradient BCs */
  *coefb = cfl / (1.0 + cfl);
  *coefa = (1.0 - *coefb) * pimp;

  /* Flux BCs */
  *cofaf = - hint * *coefa;
  *cofbf =   hint * (1.0 - *coefb);
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

    if (f->type & CS_FIELD_VARIABLE && f->bc_coeffs != NULL) {

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
 * \brief Update per variable boundary condition codes.
 *
 * \param[in]  itypfb  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_compute(int  itypfb[])
{
  CS_UNUSED(itypfb);

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
  cs_real_t *eval_buf = NULL;

  /* Loop on fields */

  const int n_fields = cs_field_n_fields();

  for (int f_id = 0; f_id < n_fields; f_id++) {
    cs_field_t  *f = cs_field_by_id(f_id);
    const cs_equation_param_t *eqp = NULL;
    if (f->type & CS_FIELD_VARIABLE)
      eqp = cs_field_get_equation_param_const(f);

    if (eqp == NULL)
      continue;

    /* Only handle legacy discretization here */

    if (eqp->space_scheme != CS_SPACE_SCHEME_LEGACY)
      continue;

    /* Settings in eqp may not be well-defined at this stage. The following
       test should be more appropriate to decide if one skips this field or
       not */
    if (f->type & CS_FIELD_CDO)
      continue;

    if (f->bc_coeffs == NULL)
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
      BFT_FREE(eval_buf);
      BFT_MALLOC(eval_buf, eval_buf_size, cs_real_t);
    }

    /* Loop on boundary conditions */

    for (int bc_id = 0; bc_id < eqp->n_bc_defs; bc_id++) {
      const cs_xdef_t *def = eqp->bc_defs[bc_id];
      cs_param_bc_type_t bc_type = (cs_param_bc_type_t)(def->meta);
      switch (bc_type) {

      case CS_PARAM_BC_HMG_DIRICHLET:
        _compute_hmg_dirichlet_bc(mesh,
                                  boundaries,
                                  eqp,
                                  def,
                                  icodcl_m,
                                  f->bc_coeffs->icodcl,
                                  f->bc_coeffs->rcodcl1);
        break;

      case CS_PARAM_BC_DIRICHLET:
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

      case CS_PARAM_BC_HMG_NEUMANN:
        _compute_hmg_neumann_bc(mesh,
                                def,
                                f->bc_coeffs->icodcl,
                                f->bc_coeffs->rcodcl3);
        break;

      case CS_PARAM_BC_NEUMANN:
        _compute_neumann_bc(mesh,
                            eqp,
                            def,
                            description,
                            t_eval,
                            f->bc_coeffs->icodcl,
                            f->bc_coeffs->rcodcl3,
                            eval_buf);
        break;

      case CS_PARAM_BC_ROBIN:
        {
          cs_lnum_t stride = 1 + f->dim + f->dim*f->dim;
          n_max_vals = stride * n_b_faces;
          if (n_max_vals > eval_buf_size) {
            eval_buf_size = n_max_vals;
            BFT_FREE(eval_buf);
            BFT_MALLOC(eval_buf, eval_buf_size, cs_real_t);
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

      case CS_PARAM_BC_WALL_PRESCRIBED:
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

  BFT_FREE(eval_buf);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Automatic adjustments for boundary condition codes.
 *
 * Currently handles mapped inlets, after the call to \ref stdtcl.
 * As portions of stdtcl are migrated to C, they should be called here,
 * before mapped inlets.
 *
 * \param[in]  itypfb  type of boundary for each face
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_complete(int  itypfb[])
{
  CS_NO_WARN_IF_UNUSED(itypfb);

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
    if (eqp != NULL) {
      if (eqp->space_scheme != CS_SPACE_SCHEME_LEGACY)
        continue;
    }

    /* Settings in eqp may not be well-defined at this stage. The following
       test should be more appropriate to decide if one skips this field or
       not */
    if (f->type & CS_FIELD_CDO)
      continue;

    if (f->bc_coeffs == NULL)
      continue;

    /* Treatment of mapped inlets */

    for (int map_id = 0; map_id < _n_bc_maps; map_id++) {

      cs_bc_map_t *bc_map = _bc_maps + map_id;

      if (bc_map->locator == NULL || ts->nt_cur <= 1)
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
                                          NULL);

        if (f == CS_F_(h)) {
          if (f->bc_coeffs != NULL) {
            int  *icodcl = f->bc_coeffs->icodcl;

            for (cs_lnum_t i = 0; i < n_faces; i++) {
              const cs_lnum_t j = (faces != NULL) ? faces[i] : i;
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
 * \brief  Set generalized BC for an anisotropic symmetric vector for a given
 *         face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     pimpv         Dirichlet value to impose on the normal
 *                              component
 * \param[in]     qimpv         Flux value to impose on the
 *                              tangential components
 * \param[in]     hint          Internal exchange coefficient
 * \param[in]     normal        normal
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_sym_vector_aniso
  (cs_real_t        coefa[3],
   cs_real_t        cofaf[3],
   cs_real_t        coefb[3][3],
   cs_real_t        cofbf[3][3],
   const cs_real_t  hint[6],
   const cs_real_t  normal[3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3])
{
  cs_real_t m[6] = {0., 0., 0., 0., 0., 0.};

  m[0] = hint[1]*hint[2] - hint[4]*hint[4];
  m[1] = hint[0]*hint[2] - hint[5]*hint[5];
  m[2] = hint[0]*hint[1] - hint[3]*hint[3];
  m[3] = hint[4]*hint[5] - hint[3]*hint[2];
  m[4] = hint[3]*hint[5] - hint[0]*hint[4];
  m[5] = hint[3]*hint[4] - hint[1]*hint[5];

  const cs_real_t invdet = 1.0/(hint[0]*m[0] + hint[3]*m[3] + hint[5]*m[5]);

  cs_real_t invh[6] = {0., 0., 0., 0., 0., 0.};
  invh[0] = m[0] * invdet;
  invh[1] = m[1] * invdet;
  invh[2] = m[2] * invdet;
  invh[3] = m[3] * invdet;
  invh[4] = m[4] * invdet;
  invh[5] = m[5] * invdet;

  cs_real_t qshint[3] = {0., 0., 0.};
  cs_real_t hintpv[3] = {0., 0., 0.};
  cs_real_t hintnm[3] = {0., 0., 0.};

  cs_math_sym_33_3_product(invh, qimpv,  qshint);
  cs_math_sym_33_3_product(hint, pimpv,  hintpv);
  cs_math_sym_33_3_product(hint, normal, hintnm);

  for (int isou = 0; isou < 3; isou++) {

    /* Gradient BCs */
    coefa[isou] = - qshint[isou];
    /* "[1 -n(x)n] Qimp / hint" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++) {

      coefa[isou] = coefa[isou] + normal[isou]*normal[jsou]
        * (pimpv[jsou] + qshint[jsou]);

      if (jsou == isou)
        coefb[isou][jsou] = 1.0 - normal[isou]*normal[jsou];
      else
        coefb[isou][jsou] = - normal[isou]*normal[jsou];
    }

    /* Flux BCs */
    cofaf[isou] = qimpv[isou];
    /* "[1 -n(x)n] Qimp" is divided into two */
    for (int jsou = 0; jsou < 3; jsou++){
      cofaf[isou] = cofaf[isou] - normal[isou]*normal[jsou]
                  * (hintpv[jsou] + qimpv[jsou]);

      cofbf[isou][jsou] = hintnm[isou] * normal[jsou];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set generalized Dirichlet BC for an anisotropic vector for a given
 *         face.
 *
 * \param[out]    coefa         explicit BC coefficient for gradients
 * \param[out]    cofaf         explicit BC coefficient for diffusive flux
 * \param[out]    coefb         implicit BC coefficient for gradients
 * \param[out]    cofbf         implicit BC coefficient for diffusive flux
 * \param[in]     hint          Internal exchange coefficient
 * \param[in]     normal        normal
 * \param[in]     pimpv         Dirichlet value to impose on the tangential
 *                              components
 * \param[in]     qimpv         Flux value to impose on the
 *                              normal component
 */
/*----------------------------------------------------------------------------*/

void
cs_boundary_conditions_set_generalized_dirichlet_vector_aniso
  (cs_real_t        coefa[3],
   cs_real_t        cofaf[3],
   cs_real_t        coefb[3][3],
   cs_real_t        cofbf[3][3],
   const cs_real_t  hint[6],
   const cs_real_t  normal[3],
   const cs_real_t  pimpv[3],
   const cs_real_t  qimpv[3])
{
  cs_real_t m[6] = {0., 0., 0., 0., 0., 0.};
  m[0] = hint[1]*hint[2] - hint[4]*hint[4];
  m[1] = hint[0]*hint[2] - hint[5]*hint[5];
  m[2] = hint[0]*hint[1] - hint[3]*hint[3];
  m[3] = hint[4]*hint[5] - hint[3]*hint[2];
  m[4] = hint[3]*hint[5] - hint[0]*hint[4];
  m[5] = hint[3]*hint[4] - hint[1]*hint[5];

  const cs_real_t invdet = 1.0/(hint[0]*m[0] + hint[3]*m[3] + hint[5]*m[5]);

  cs_real_t invh[6] = {0., 0., 0., 0., 0., 0.};
  invh[0] = m[0] * invdet;
  invh[1] = m[1] * invdet;
  invh[2] = m[2] * invdet;
  invh[3] = m[3] * invdet;
  invh[4] = m[4] * invdet;
  invh[5] = m[5] * invdet;

  cs_real_t qshint[3] = {0., 0., 0.};
  cs_real_t hintpv[3] = {0., 0., 0.};
  cs_real_t hintnm[3] = {0., 0., 0.};

  cs_math_sym_33_3_product(invh, qimpv,  qshint);
  cs_math_sym_33_3_product(hint, pimpv,  hintpv);
  cs_math_sym_33_3_product(hint, normal, hintnm);

  for (int isou = 0; isou < 3; isou ++) {

    /* Gradient BCs */
    /* "[1 -n(x)n] Pimp" is divided into two */
    coefa[isou] = pimpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      coefa[isou] = coefa[isou] - normal[isou] * normal[jsou]
                  * (pimpv[jsou] + qshint[jsou]);

      coefb[isou][jsou] = normal[isou] * normal[jsou];
    }

    /* Flux BCs */
    /* "[1 -n(x)n] Pimp" is divided into two */
    cofaf[isou] = -hintpv[isou];
    for (int jsou = 0; jsou < 3; jsou++) {

      cofaf[isou] = cofaf[isou] + normal[isou]*normal[jsou]
        *(qimpv[jsou]+hintpv[jsou]);

      if (jsou == isou)
        cofbf[isou][jsou] = hint[isou]-hintnm[isou]*normal[jsou];
      else
        cofbf[isou][jsou] = -hintnm[isou]*normal[jsou];
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
