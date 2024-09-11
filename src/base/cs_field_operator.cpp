/*============================================================================
 * Field based algebraic operators.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_gradient.h"
#include "cs_gradient_boundary.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_mesh.h"
#include "cs_log.h"
#include "cs_map.h"
#include "cs_parameters.h"
#include "cs_parall.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_internal_coupling.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_field_operator.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_field_operator.c
        Field based algebraic operators.

  \enum cs_field_interpolate_t

  \brief Field interpolation modes

  \var CS_FIELD_INTERPOLATE_MEAN
       Mean element value (P0 interpolation)
  \var CS_FIELD_INTERPOLATE_GRADIENT
       Mean element value + gradient correction (pseudo-P1)
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for functions intended for use only by Fortran wrappers.
 * (descriptions follow, with function bodies).
 *============================================================================*/

void
cs_f_field_gradient_scalar(int                    f_id,
                           int                    use_previous_t,
                           int                    inc,
                           cs_real_3_t  *restrict grad);

void
cs_f_field_gradient_vector(int                     f_id,
                           int                     use_previous_t,
                           int                     inc,
                           cs_real_33_t  *restrict grad);

void
cs_f_field_gradient_tensor(int                     f_id,
                           int                     use_previous_t,
                           int                     inc,
                           cs_real_63_t  *restrict grad);

void
cs_f_field_set_volume_average(int       f_id,
                              cs_real_t va);

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Interpolate field values at a given set of points using P0 interpolation.
 *
 * parameters:
 *   f                  <-- pointer to field
 *   n_points           <-- number of points at which interpolation
 *                          is required
 *   point_location     <-- location of points in mesh elements
 *                          (based on the field location)
 *   val                --> interpolated values
 *----------------------------------------------------------------------------*/

static void
_field_interpolate_by_mean(const cs_field_t   *f,
                           cs_lnum_t           n_points,
                           const cs_lnum_t     point_location[],
                           cs_real_t          *val)
{
  for (cs_lnum_t i = 0; i < n_points; i++) {

    cs_lnum_t cell_id = point_location[i];

    for (cs_lnum_t j = 0; j < f->dim; j++)
      val[i*f->dim + j] =  f->val[cell_id*f->dim + j];

  }
}

/*----------------------------------------------------------------------------
 * Interpolate field values at a given set of points using gradient-corrected
 * interpolation.
 *
 * parameters:
 *   f                  <-- pointer to field
 *   n_points           <-- number of points at which interpolation
 *                          is required
 *   point_location     <-- location of points in mesh elements
 *                          (based on the field location)
 *   point_coords       <-- point coordinates
 *   val                --> interpolated values
 *----------------------------------------------------------------------------*/

static void
_field_interpolate_by_gradient(const cs_field_t   *f,
                               cs_lnum_t           n_points,
                               const cs_lnum_t     point_location[],
                               const cs_real_3_t   point_coords[],
                               cs_real_t          *val)
{
  const cs_lnum_t dim = f->dim;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)(cs_glob_mesh_quantities->cell_cen);

  /* Currently possible only for fields on cell location */

  if (f->location_id != CS_MESH_LOCATION_CELLS)
    bft_error(__FILE__, __LINE__, 0,
              _("Field gradient interpolation for field %s :\n"
                " not implemented for fields on location %s."),
              f->name, cs_mesh_location_type_name[f->location_id]);

  /* Compute field cell gradient */

  cs_real_t *grad;
  BFT_MALLOC(grad, 3*dim*n_cells_ext, cs_real_t);

  if (dim == 1)
    cs_field_gradient_scalar(f,
                             true, /* use_previous_t */
                             1,    /* inc */
                             (cs_real_3_t *)grad);
  else if (dim == 3)
    cs_field_gradient_vector(f,
                             true, /* use_previous_t */
                             1,    /* inc */
                             (cs_real_33_t *)grad);

  else
    bft_error(__FILE__, __LINE__, 0,
              _("Field gradient interpolation for field %s of dimension %d:\n"
                " not implemented."),
              f->name, (int)dim);

  /* Now interpolated values */

  for (cs_lnum_t i = 0; i < n_points; i++) {

    cs_lnum_t cell_id = point_location[i];

    cs_real_3_t d = {point_coords[i][0] - cell_cen[cell_id][0],
                     point_coords[i][1] - cell_cen[cell_id][1],
                     point_coords[i][2] - cell_cen[cell_id][2]};

    for (cs_lnum_t j = 0; j < f->dim; j++) {
      cs_lnum_t k = (cell_id*dim + j)*3;
      val[i*dim + j] =   f->val[cell_id*dim + j]
                       + d[0] * grad[k]
                       + d[1] * grad[k+1]
                       + d[2] * grad[k+2];
    }

  }

  BFT_FREE(grad);
}

/*----------------------------------------------------------------------------
 * For each mesh cell this function finds the local extrema of a
 * scalar field.
 *
 * parameters:
 *   pvar            <-- scalar values
 *   halo_type       <-- halo type
 *   local_max       --> local maximum value
 *   local_min       --> local minimum value
 *----------------------------------------------------------------------------*/

static void
_local_extrema_scalar(const cs_real_t *restrict pvar,
                      cs_halo_type_t            halo_type,
                      cs_real_t                *local_max,
                      cs_real_t                *local_min)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_vertices = m->n_vertices;

  const cs_adjacency_t  *c2v = cs_mesh_adjacencies_cell_vertices();
  const cs_lnum_t *c2v_idx = c2v->idx;
  const cs_lnum_t *c2v_ids = c2v->ids;

# pragma omp parallel for  if (n_cells > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_cells; ii++) {
    local_max[ii] = pvar[ii];
    local_min[ii] = pvar[ii];
  }

  cs_real_t *v_min, *v_max;
  BFT_MALLOC(v_min, n_vertices, cs_real_t);
  BFT_MALLOC(v_max, n_vertices, cs_real_t);

# pragma omp parallel for  if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_vertices; ii++) {
    v_max[ii] = -HUGE_VAL;
    v_min[ii] = HUGE_VAL;
  }

  /* Scatter min/max values to vertices */

  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];
    cs_real_t _c_var = pvar[c_id];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      if (_c_var < v_min[v_id])
        v_min[v_id] = _c_var;
      if (_c_var > v_max[v_id])
        v_max[v_id] = _c_var;
    }
  }

  if (m->vtx_interfaces != NULL) {
    cs_interface_set_min(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         v_min);
    cs_interface_set_max(m->vtx_interfaces,
                         m->n_vertices,
                         1,
                         true,
                         CS_REAL_TYPE,
                         v_max);
  }

  /* Gather min/max values from vertices */

# pragma omp parallel for  if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    cs_lnum_t s_id = c2v_idx[c_id];
    cs_lnum_t e_id = c2v_idx[c_id+1];
    for (cs_lnum_t j = s_id; j < e_id; j++) {
      cs_lnum_t v_id = c2v_ids[j];
      if (v_min[v_id] < local_min[c_id])
        local_min[c_id] = v_min[v_id];
      if (v_max[v_id] > local_max[c_id])
        local_max[c_id] = v_max[v_id];
    }
  }

  /* Free memory */
  BFT_FREE(v_min);
  BFT_FREE(v_max);

  /* Synchronisation */
  if (m->halo != NULL) {
    cs_halo_sync_var(m->halo, halo_type, local_min);
    cs_halo_sync_var(m->halo, halo_type, local_max);
  }
}

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f_id           <-- field id
 *   use_previous_t <-- should we use values from the previous time step ?
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void
cs_f_field_gradient_scalar(int                    f_id,
                           int                    use_previous_t,
                           int                    inc,
                           cs_real_3_t  *restrict grad)
{
  bool _use_previous_t = use_previous_t ? true : false;

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_field_gradient_scalar(f,
                           _use_previous_t,
                           inc,
                           grad);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f_id           <-- field id
 *   use_previous_t <-- should we use values from the previous time step ?
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void
cs_f_field_gradient_vector(int                     f_id,
                           int                     use_previous_t,
                           int                     inc,
                           cs_real_33_t  *restrict grad)
{
  bool _use_previous_t = use_previous_t ? true : false;

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_field_gradient_vector(f,
                           _use_previous_t,
                           inc,
                           grad);
}

/*----------------------------------------------------------------------------
 * Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * parameters:
 *   f_id           <-- field id
 *   use_previous_t <-- should we use values from the previous time step ?
 *   inc            <-- if 0, solve on increment; 1 otherwise
 *   grad           --> gradient
 *----------------------------------------------------------------------------*/

void
cs_f_field_gradient_tensor(int                     f_id,
                           int                     use_previous_t,
                           int                     inc,
                           cs_real_63_t  *restrict grad)
{
  bool _use_previous_t = use_previous_t ? true : false;

  const cs_field_t *f = cs_field_by_id(f_id);

  cs_field_gradient_tensor(f,
                           _use_previous_t,
                           inc,
                           grad);
}

/*----------------------------------------------------------------------------
 * Shift field values in order to set its spatial average to a given value.
 *
 * parameters:
 *   f_id           <-- field id
 *   va             <-- real value of volume average to be set
 *----------------------------------------------------------------------------*/

void
cs_f_field_set_volume_average(int       f_id,
                              cs_real_t va)
{
  cs_field_t *f = cs_field_by_id(f_id);

  cs_field_set_volume_average(f,
                              va);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_field_gradient_scalar(const cs_field_t          *f,
                         bool                       use_previous_t,
                         int                        inc,
                         cs_real_3_t      *restrict grad)
{
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  /* Does the field have a parent (variable) ?
     Field is its own parent if not parent is specified */

  const cs_field_t *parent_f = f;

  const int f_parent_id
    = cs_field_get_key_int(f, cs_field_key_id("parent_field_id"));
  if (f_parent_id > -1)
    parent_f = cs_field_by_id(f_parent_id);

  int imrgra = cs_glob_space_disc->imrgra;
  cs_equation_param_t eqp_default = cs_parameters_equation_param_default();

  /* Get the calculation option from the field */
  const cs_equation_param_t
    *eqp = cs_field_get_equation_param_const(parent_f);

  if (eqp != NULL)
    imrgra = eqp->imrgra;
  else
    eqp = &eqp_default;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  int w_stride = 1;
  cs_real_t *c_weight = NULL;
  cs_internal_coupling_t  *cpl = NULL;

  if (parent_f->type & CS_FIELD_VARIABLE && eqp->idiff > 0) {

    if (eqp->iwgrec == 1) {
      /* Weighted gradient coefficients */
      int key_id = cs_field_key_id("gradient_weighting_id");
      int diff_id = cs_field_get_key_int(parent_f, key_id);
      if (diff_id > -1) {
        cs_field_t *f_weight = cs_field_by_id(diff_id);
        c_weight = f_weight->val;
        w_stride = f_weight->dim;
      }
    }

    /* Internal coupling structure */
    int key_id = cs_field_key_id_try("coupling_entity");
    if (key_id > -1) {
      int coupl_id = cs_field_get_key_int(parent_f, key_id);
      if (coupl_id > -1)
        cpl = cs_internal_coupling_by_id(coupl_id);
    }

  }

  if (f->n_time_vals < 2 && use_previous_t)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: field %s does not maintain previous time step values\n"
                "so \"use_previous_t\" can not be handled."),
              __func__, f->name);

  cs_real_t *var = (use_previous_t) ? f->val_pre : f->val;

  cs_field_bc_coeffs_t *bc_coeffs = NULL;
  if (f->bc_coeffs != NULL) {
    bc_coeffs = f->bc_coeffs;
  }

  cs_gradient_scalar(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     eqp->nswrgr,
                     0, /* hyd_p_flag */
                     w_stride,
                     eqp->verbosity,
                     static_cast<cs_gradient_limit_t>(eqp->imligr),
                     eqp->epsrgr,
                     eqp->climgr,
                     NULL, /* f_ext */
                     bc_coeffs,
                     var,
                     c_weight,
                     cpl, /* internal coupling */
                     grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of scalar array using parameters associated
 *         with a given field.
 *
 * \param[in]       f_id           associated field id
 * \param[in]       inc            if 0, solve on increment; 1 otherwise
 * \param[in]       bc_coeffs      boundary condition structure
 * \param[in, out]  var            gradient's base variable
 * \param[out]      grad           gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_field_gradient_scalar_array(int                         f_id,
                               int                         inc,
                               const cs_field_bc_coeffs_t *bc_coeffs,
                               cs_real_t                   var[],
                               cs_real_3_t                 grad[])
{
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_LSQ;

  const cs_field_t *f = cs_field_by_id(f_id);
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  /* Choose gradient type */

  cs_gradient_type_by_imrgra(eqp->imrgra,
                             &gradient_type,
                             &halo_type);

  /* Check if given field has internal coupling  */
  cs_internal_coupling_t  *cpl = NULL;
  if (f_id > -1) {
    const int key_id = cs_field_key_id_try("coupling_entity");
    if (key_id > -1) {
      int coupl_id = cs_field_get_key_int(f, key_id);
      if (coupl_id > -1)
        cpl = cs_internal_coupling_by_id(coupl_id);
    }
  }

  /* Compute gradient */

  cs_gradient_scalar(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     eqp->nswrgr,
                     0,             /* iphydp */
                     1,             /* w_stride */
                     eqp->verbosity,
                     (cs_gradient_limit_t)(eqp->imligr),
                     eqp->epsrgr,
                     eqp->climgr,
                     NULL,          /* f_ext */
                     bc_coeffs,
                     var,
                     NULL,          /* c_weight */
                     cpl,
                     grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of scalar field or component of vector or
 * tensor field.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[in]       hyd_p_flag      flag for hydrostatic pressure
 * \param[in]       f_ext           exterior force generating
 *                                  the hydrostatic pressure
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_field_gradient_potential(const cs_field_t          *f,
                            bool                       use_previous_t,
                            int                        inc,
                            int                        hyd_p_flag,
                            cs_real_3_t                f_ext[],
                            cs_real_3_t      *restrict grad)
{
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  /* Does the field have a parent (variable) ?
     Field is its own parent if not parent is specified */

  const cs_field_t *parent_f = f;

  const int f_parent_id
    = cs_field_get_key_int(f, cs_field_key_id("parent_field_id"));
  if (f_parent_id > -1)
    parent_f = cs_field_by_id(f_parent_id);

  int imrgra = cs_glob_space_disc->imrgra;
  cs_equation_param_t eqp_default = cs_parameters_equation_param_default();

  /* Get the calculation option from the field */
  const cs_equation_param_t
    *eqp = cs_field_get_equation_param_const(parent_f);

  if (eqp != NULL)
    imrgra = eqp->imrgra;
  else
    eqp = &eqp_default;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f->n_time_vals < 2 && use_previous_t)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: field %s does not maintain previous time step values\n"
                "so \"use_previous_t\" can not be handled."),
              __func__, f->name);

  int w_stride = 1;
  cs_real_t *var = (use_previous_t) ? f->val_pre : f->val;

  cs_real_t *c_weight = NULL;
  cs_internal_coupling_t  *cpl = NULL;

  if (parent_f->type & CS_FIELD_VARIABLE && eqp->idiff > 0) {

    if (eqp->iwgrec == 1) {
      int key_id = cs_field_key_id("gradient_weighting_id");
      int diff_id = cs_field_get_key_int(parent_f, key_id);
      if (diff_id > -1) {
        cs_field_t *f_weight = cs_field_by_id(diff_id);
        c_weight = f_weight->val;
        w_stride = f_weight->dim;
      }
    }

    /* Internal coupling structure */
    int key_id = cs_field_key_id_try("coupling_entity");
    if (key_id > -1) {
      int coupl_id = cs_field_get_key_int(parent_f, key_id);
      if (coupl_id > -1)
        cpl = cs_internal_coupling_by_id(coupl_id);
    }

  }

  const cs_field_bc_coeffs_t *bc_coeffs = NULL;
  if (f->bc_coeffs != NULL) {
    bc_coeffs = f->bc_coeffs;
  }

  if (hyd_p_flag == 2)
    hyd_p_flag = 0;

  cs_gradient_scalar(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     eqp->nswrgr,
                     hyd_p_flag,
                     w_stride,
                     eqp->verbosity,
                     static_cast<cs_gradient_limit_t>(eqp->imligr),
                     eqp->epsrgr,
                     eqp->climgr,
                     f_ext,
                     bc_coeffs,
                     var,
                     c_weight,
                     cpl, /* internal coupling */
                     grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of vector field.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_field_gradient_vector(const cs_field_t          *f,
                         bool                       use_previous_t,
                         int                        inc,
                         cs_real_33_t     *restrict grad)
{
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  int imrgra = cs_glob_space_disc->imrgra;
  cs_equation_param_t eqp_default = cs_parameters_equation_param_default();

  /* Get the calculation option from the field */
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  if (eqp != NULL)
    imrgra = eqp->imrgra;
  else
    eqp = &eqp_default;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  cs_real_t *c_weight = NULL;
  cs_internal_coupling_t  *cpl = NULL;

  if (f->type & CS_FIELD_VARIABLE && eqp->idiff > 0) {

    if (eqp->iwgrec == 1) {
      /* Weighted gradient coefficients */
      int key_id = cs_field_key_id("gradient_weighting_id");
      int diff_id = cs_field_get_key_int(f, key_id);
      if (diff_id > -1) {
        cs_field_t *f_weight = cs_field_by_id(diff_id);
        c_weight = f_weight->val;
      }
    }

    int key_id = cs_field_key_id_try("coupling_entity");
    if (key_id > -1) {
      int coupl_id = cs_field_get_key_int(f, key_id);
      if (coupl_id > -1)
        cpl = cs_internal_coupling_by_id(coupl_id);
    }

  }

  if (f->n_time_vals < 2 && use_previous_t)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: field %s does not maintain previous time step values\n"
                "so \"use_previous_t\" can not be handled."),
              __func__, f->name);

  cs_real_3_t *var = (use_previous_t) ? (cs_real_3_t *)(f->val_pre)
                                      : (cs_real_3_t *)(f->val);

  const cs_field_bc_coeffs_t *bc_coeffs = NULL;

  if (f->bc_coeffs != NULL) {
    /* coupled components */
    int coupled_key_id = cs_field_key_id_try("coupled");
    if (coupled_key_id > 1) {
      if (cs_field_get_key_int(f, coupled_key_id) > 0) {
        bc_coeffs = f->bc_coeffs;
      }
    }
  }

  cs_gradient_vector(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     eqp->nswrgr,
                     eqp->verbosity,
                     static_cast<cs_gradient_limit_t>(eqp->imligr),
                     eqp->epsrgr,
                     eqp->climgr,
                     bc_coeffs,
                     var,
                     c_weight,
                     cpl,
                     grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute cell gradient of tensor field.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_field_gradient_tensor(const cs_field_t          *f,
                         bool                       use_previous_t,
                         int                        inc,
                         cs_real_63_t     *restrict grad)
{
  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  int imrgra = cs_glob_space_disc->imrgra;
  cs_equation_param_t eqp_default = cs_parameters_equation_param_default();

  /* Get the calculation option from the field */
  const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);

  if (eqp != NULL)
    imrgra = eqp->imrgra;
  else
    eqp = &eqp_default;

  cs_gradient_type_by_imrgra(imrgra,
                             &gradient_type,
                             &halo_type);

  if (f->n_time_vals < 2 && use_previous_t)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: field %s does not maintain previous time step values\n"
                "so \"use_previous_t\" can not be handled."),
              __func__, f->name);

  cs_real_6_t *var = (use_previous_t) ? (cs_real_6_t *)(f->val_pre)
                                      : (cs_real_6_t *)(f->val);

  const cs_field_bc_coeffs_t *bc_coeffs_ts;
  if (f->bc_coeffs != NULL) {
    int coupled_key_id = cs_field_key_id_try("coupled");
    if (coupled_key_id > 1) {
      if (cs_field_get_key_int(f, coupled_key_id) > 0) {
        bc_coeffs_ts = f->bc_coeffs;
      }
    }
  }

  cs_gradient_tensor(f->name,
                     gradient_type,
                     halo_type,
                     inc,
                     eqp->nswrgr,
                     eqp->verbosity,
                     static_cast<cs_gradient_limit_t>(eqp->imligr),
                     eqp->epsrgr,
                     eqp->climgr,
                     bc_coeffs_ts,
                     var,
                     grad);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of a scalar field at boundary face I' positions.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_field_gradient_boundary_iprime_scalar(const cs_field_t  *f,
                                         bool               use_previous_t,
                                         cs_lnum_t          n_faces,
                                         const cs_lnum_t   *face_ids,
                                         cs_real_t          var_iprime[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  /* Does the field have a parent (variable) ?
     Field is its own parent if not parent is specified */

  const cs_field_t *parent_f = f;

  const int f_parent_id
    = cs_field_get_key_int(f, cs_field_key_id("parent_field_id"));
  if (f_parent_id > -1)
    parent_f = cs_field_by_id(f_parent_id);

  int b_gradient_r = cs_glob_space_disc->imrgra;
  cs_equation_param_t eqp_default = cs_parameters_equation_param_default();

  /* Get the calculation option from the field */
  const cs_equation_param_t
    *eqp = cs_field_get_equation_param_const(parent_f);

  if (eqp != NULL)
    b_gradient_r = eqp->b_gradient_r;
  else
    eqp = &eqp_default;

  cs_gradient_type_by_imrgra(b_gradient_r,
                             &gradient_type,
                             &halo_type);

  int w_stride = 1;
  cs_real_t *c_weight = NULL;
  cs_internal_coupling_t  *cpl = NULL;

  if (parent_f->type & CS_FIELD_VARIABLE && eqp->idiff > 0) {

    if (eqp->iwgrec == 1) {
      /* Weighted gradient coefficients */
      int key_id = cs_field_key_id("gradient_weighting_id");
      int diff_id = cs_field_get_key_int(parent_f, key_id);
      if (diff_id > -1) {
        cs_field_t *f_weight = cs_field_by_id(diff_id);
        c_weight = f_weight->val;
        w_stride = f_weight->dim;
      }
    }

    /* Internal coupling structure */
    int key_id = cs_field_key_id_try("coupling_entity");
    if (key_id > -1) {
      int coupl_id = cs_field_get_key_int(parent_f, key_id);
      if (coupl_id > -1)
        cpl = cs_internal_coupling_by_id(coupl_id);
    }

  }

  if (f->n_time_vals < 2 && use_previous_t)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: field %s does not maintain previous time step values\n"
                "so \"use_previous_t\" can not be handled."),
              __func__, f->name);

  cs_real_t *var = (use_previous_t) ? f->val_pre : f->val;

  const cs_field_bc_coeffs_t *bc_coeffs = NULL;
  if (f->bc_coeffs != NULL) {
    bc_coeffs = f->bc_coeffs;
  }

  /* With least-squares gradient, we can use a cheaper, boundary-only
     reconstruction */

  if (gradient_type == CS_GRADIENT_LSQ) {

    cs_real_t climgr = (eqp->imligr < 0) ? -1.0 : eqp->climgr;

    cs_gradient_boundary_iprime_lsq_s(m,
                                      fvq,
                                      n_faces,
                                      face_ids,
                                      halo_type,
                                      climgr,
                                      bc_coeffs,
                                      c_weight,
                                      var,
                                      var_iprime);

  }
  else {

    const cs_real_3_t *restrict diipb
      = (const cs_real_3_t *)fvq->diipb;

    cs_real_3_t *grad;
    CS_MALLOC_HD(grad,
                 cs_glob_mesh->n_cells_with_ghosts,
                 cs_real_3_t,
                 cs_alloc_mode);

    cs_gradient_scalar(f->name,
                       gradient_type,
                       halo_type,
                       1, /* inc */
                       eqp->nswrgr,
                       0, /* hyd_p_flag */
                       w_stride,
                       eqp->verbosity,
                       static_cast<cs_gradient_limit_t>(eqp->imligr),
                       eqp->epsrgr,
                       eqp->climgr,
                       NULL, /* f_ext */
                       bc_coeffs,
                       var,
                       c_weight,
                       cpl, /* internal coupling */
                       grad);

    /* Finally, reconstruct value at I' */

    const cs_lnum_t *restrict b_face_cells
      = (const cs_lnum_t *)m->b_face_cells;

    if (face_ids != NULL) {
      for (cs_lnum_t idx = 0; idx < n_faces; idx++) {
        cs_lnum_t i = face_ids[idx];
        cs_lnum_t j = b_face_cells[i];
        var_iprime[idx] = var[j] + (  grad[j][0]*diipb[i][0]
                                    + grad[j][1]*diipb[i][1]
                                    + grad[j][2]*diipb[i][2]);
      }
    }
    else {
      for (cs_lnum_t i = 0; i < n_faces; i++) {
        cs_lnum_t j = b_face_cells[i];
        var_iprime[i] = var[j] + (  grad[j][0]*diipb[i][0]
                                  + grad[j][1]*diipb[i][1]
                                  + grad[j][2]*diipb[i][2]);
      }
    }

    CS_FREE_HD(grad);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of a vector field at boundary face I' positions.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_field_gradient_boundary_iprime_vector(const cs_field_t  *f,
                                         bool               use_previous_t,
                                         cs_lnum_t          n_faces,
                                         const cs_lnum_t   *face_ids,
                                         cs_real_3_t        var_iprime[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  /* Does the field have a parent (variable) ?
     Field is its own parent if not parent is specified */

  const cs_field_t *parent_f = f;

  const int f_parent_id
    = cs_field_get_key_int(f, cs_field_key_id("parent_field_id"));
  if (f_parent_id > -1)
    parent_f = cs_field_by_id(f_parent_id);

  int b_gradient_r = cs_glob_space_disc->imrgra;
  cs_equation_param_t eqp_default = cs_parameters_equation_param_default();

  /* Get the calculation option from the field */
  const cs_equation_param_t
    *eqp = cs_field_get_equation_param_const(parent_f);

  if (eqp != NULL)
    b_gradient_r = eqp->b_gradient_r;
  else
    eqp = &eqp_default;

  cs_gradient_type_by_imrgra(b_gradient_r,
                             &gradient_type,
                             &halo_type);

  cs_real_t *c_weight = NULL;
  cs_internal_coupling_t  *cpl = NULL;

  if (parent_f->type & CS_FIELD_VARIABLE && eqp->idiff > 0) {

    /* Internal coupling structure */
    int key_id = cs_field_key_id_try("coupling_entity");
    if (key_id > -1) {
      int coupl_id = cs_field_get_key_int(parent_f, key_id);
      if (coupl_id > -1)
        cpl = cs_internal_coupling_by_id(coupl_id);
    }

  }

  if (f->n_time_vals < 2 && use_previous_t)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: field %s does not maintain previous time step values\n"
                "so \"use_previous_t\" can not be handled."),
              __func__, f->name);

  const cs_field_bc_coeffs_t *bc_coeffs_v = NULL;

  if (f->bc_coeffs != NULL) {
    /* coupled components */
    int coupled_key_id = cs_field_key_id_try("coupled");
    if (coupled_key_id > 1) {
      if (cs_field_get_key_int(f, coupled_key_id) > 0) {
        bc_coeffs_v = f->bc_coeffs;
      }
    }
  }

  /* With least-squares gradient, we can use a cheaper, boundary-only
     reconstruction */

  if (gradient_type == CS_GRADIENT_LSQ) {

    cs_real_t climgr = (eqp->imligr < 0) ? -1.0 : eqp->climgr;

    const cs_real_3_t *var = (use_previous_t) ?
                               (const cs_real_3_t *)f->val_pre
                             : (const cs_real_3_t *)f->val;

    cs_gradient_boundary_iprime_lsq_v(m,
                                      fvq,
                                      n_faces,
                                      face_ids,
                                      halo_type,
                                      climgr,
                                      bc_coeffs_v,
                                      c_weight,
                                      var,
                                      var_iprime);

  }
  else {

    cs_real_3_t *var = (use_previous_t) ?
                         (cs_real_3_t *)f->val_pre
                       : (cs_real_3_t *)f->val;

    const cs_real_3_t *restrict diipb
      = (const cs_real_3_t *)fvq->diipb;

    cs_real_33_t *grad;
    CS_MALLOC_HD(grad,
                 cs_glob_mesh->n_cells_with_ghosts,
                 cs_real_33_t,
                 cs_alloc_mode);

    cs_gradient_vector(f->name,
                       gradient_type,
                       halo_type,
                       1, /* inc */
                       eqp->nswrgr,
                       eqp->verbosity,
                       static_cast<cs_gradient_limit_t>(eqp->imligr),
                       eqp->epsrgr,
                       eqp->climgr,
                       bc_coeffs_v,
                       var,
                       c_weight,
                       cpl, /* internal coupling */
                       grad);

    /* Finally, reconstruct value at I' */

    const cs_lnum_t *restrict b_face_cells
      = (const cs_lnum_t *)m->b_face_cells;

    for (cs_lnum_t idx = 0; idx < n_faces; idx++) {
      cs_lnum_t i = (face_ids != NULL) ? face_ids[idx] : idx;
      cs_lnum_t j = b_face_cells[i];
      for (cs_lnum_t k = 0; k < 3; k++) {
        var_iprime[idx][k] =   var[j][k]
                             + cs_math_3_dot_product(grad[j][k],
                                                     diipb[i]);
      }
    }

    CS_FREE_HD(grad);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the values of a symmetric tensor field at
 *        boundary face I' positions.
 *
 * \param[in]       f               pointer to field
 * \param[in]       use_previous_t  should we use values from the previous
 *                                  time step ?
 * \param[in]       inc             if 0, solve on increment; 1 otherwise
 * \param[out]      grad            gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_field_gradient_boundary_iprime_tensor(const cs_field_t  *f,
                                         bool               use_previous_t,
                                         cs_lnum_t          n_faces,
                                         const cs_lnum_t   *face_ids,
                                         cs_real_6_t        var_iprime[])
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  cs_halo_type_t halo_type = CS_HALO_STANDARD;
  cs_gradient_type_t gradient_type = CS_GRADIENT_GREEN_ITER;

  /* Does the field have a parent (variable) ?
     Field is its own parent if not parent is specified */

  const cs_field_t *parent_f = f;

  const int f_parent_id
    = cs_field_get_key_int(f, cs_field_key_id("parent_field_id"));
  if (f_parent_id > -1)
    parent_f = cs_field_by_id(f_parent_id);

  int b_gradient_r = cs_glob_space_disc->imrgra;
  cs_equation_param_t eqp_default = cs_parameters_equation_param_default();

  /* Get the calculation option from the field */
  const cs_equation_param_t
    *eqp = cs_field_get_equation_param_const(parent_f);

  if (eqp != NULL)
    b_gradient_r = eqp->b_gradient_r;
  else
    eqp = &eqp_default;

  cs_gradient_type_by_imrgra(b_gradient_r,
                             &gradient_type,
                             &halo_type);

  if (f->n_time_vals < 2 && use_previous_t)
    bft_error(__FILE__, __LINE__, 0,
              _("%s: field %s does not maintain previous time step values\n"
                "so \"use_previous_t\" can not be handled."),
              __func__, f->name);

  const cs_field_bc_coeffs_t *bc_coeffs_ts = NULL;

  if (f->bc_coeffs != NULL) {
    /* coupled components */
    int coupled_key_id = cs_field_key_id_try("coupled");
    if (coupled_key_id > 1) {
      if (cs_field_get_key_int(f, coupled_key_id) > 0) {
        bc_coeffs_ts = f->bc_coeffs;
      }
    }
  }

  /* With least-squares gradient, we can use a cheaper, boundary-only
     reconstruction */

  if (gradient_type == CS_GRADIENT_LSQ) {

    cs_real_t climgr = (eqp->imligr < 0) ? -1.0 : eqp->climgr;

    const cs_real_6_t *var = (use_previous_t) ?
                               (const cs_real_6_t *)f->val_pre
                             : (const cs_real_6_t *)f->val;

    cs_gradient_boundary_iprime_lsq_t(m,
                                      fvq,
                                      n_faces,
                                      face_ids,
                                      halo_type,
                                      climgr,
                                      bc_coeffs_ts,
                                      NULL,
                                      var,
                                      var_iprime);

  }
  else {

    cs_real_6_t *var = (use_previous_t) ?
                         (cs_real_6_t *)f->val_pre
                       : (cs_real_6_t *)f->val;

    const cs_real_3_t *restrict diipb
      = (const cs_real_3_t *)fvq->diipb;

    cs_real_63_t *grad;
    CS_MALLOC_HD(grad,
                 cs_glob_mesh->n_cells_with_ghosts,
                 cs_real_63_t,
                 cs_alloc_mode);

    cs_gradient_tensor(f->name,
                       gradient_type,
                       halo_type,
                       1, /* inc */
                       eqp->nswrgr,
                       eqp->verbosity,
                       static_cast<cs_gradient_limit_t>(eqp->imligr),
                       eqp->epsrgr,
                       eqp->climgr,
                       bc_coeffs_ts,
                       var,
                       grad);

    /* Finally, reconstruct value at I' */

    const cs_lnum_t *restrict b_face_cells
      = (const cs_lnum_t *)m->b_face_cells;

    for (cs_lnum_t idx = 0; idx < n_faces; idx++) {
      cs_lnum_t i = (face_ids != NULL) ? face_ids[idx] : idx;
      cs_lnum_t j = b_face_cells[i];
      for (cs_lnum_t k = 0; k < 6; k++) {
        var_iprime[idx][k] =   var[j][k]
                             + cs_math_3_dot_product(grad[j][k],
                                                     diipb[i]);
      }
    }

    CS_FREE_HD(grad);

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Interpolate field values at a given set of points.
 *
 * \param[in]   f                   pointer to field
 * \param[in]   interpolation_type  interpolation type
 * \param[in]   n_points            number of points at which interpolation
 *                                  is required
 * \param[in]   point_location      location of points in mesh elements
 *                                  (based on the field location)
 * \param[in]   point_coords        point coordinates
 * \param[out]  val                 interpolated values
 */
/*----------------------------------------------------------------------------*/

void
cs_field_interpolate(cs_field_t              *f,
                     cs_field_interpolate_t   interpolation_type,
                     cs_lnum_t                n_points,
                     const cs_lnum_t          point_location[],
                     const cs_real_3_t        point_coords[],
                     cs_real_t               *val)
{
  switch (interpolation_type) {

  case CS_FIELD_INTERPOLATE_MEAN:
    _field_interpolate_by_mean(f,
                               n_points,
                               point_location,
                               val);
    break;

  case CS_FIELD_INTERPOLATE_GRADIENT:
    _field_interpolate_by_gradient(f,
                                   n_points,
                                   point_location,
                                   point_coords,
                                   val);
    break;

  default:
    assert(0);
    break;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Find local extrema of a given scalar field at each cell
 *
 * This assumes the field values have been synchronized.
 *
 * \param[in]       f_id        scalar field id
 * \param[in]       halo_type   halo type
 * \param[in, out]  local_max   local maximum value
 * \param[in, out]  local_min   local minimum value
 */
/*----------------------------------------------------------------------------*/

void
cs_field_local_extrema_scalar(int              f_id,
                              cs_halo_type_t   halo_type,
                              cs_real_t       *local_max,
                              cs_real_t       *local_min)
{
  const cs_mesh_t *m = cs_glob_mesh;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_field_t *f = cs_field_by_id(f_id);
  const cs_real_t *restrict pvar = f->val;

  _local_extrema_scalar(pvar,
                        halo_type,
                        local_max,
                        local_min);

  /* Initialisation of local extrema */

  int key_scamax_id = cs_field_key_id("max_scalar");
  int key_scamin_id = cs_field_key_id("min_scalar");

  cs_real_t scalar_max = cs_field_get_key_double(f, key_scamax_id);
  cs_real_t scalar_min = cs_field_get_key_double(f, key_scamin_id);

  /* Bounded by the global extrema */
# pragma omp parallel for if (n_cells_ext > CS_THR_MIN)
  for (cs_lnum_t ii = 0; ii < n_cells_ext; ii++) {
    local_max[ii] = CS_MIN(local_max[ii], scalar_max);
    local_min[ii] = CS_MAX(local_min[ii], scalar_min);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Shift field values in order to set its spatial average (fluid volume
 * average) to a given value.
 *
 * \param[in]   f   pointer to field
 * \param[in]   va  real value of fluid volume average to be set
 */
/*----------------------------------------------------------------------------*/

void
cs_field_set_volume_average(cs_field_t     *f,
                            const cs_real_t va)
{
  const cs_lnum_t n_cells = cs_glob_mesh->n_cells;

  cs_mesh_quantities_t  *mq = cs_glob_mesh_quantities;
  const cs_real_t *cell_f_vol = mq->cell_f_vol;
  const cs_real_t tot_f_vol = mq->tot_f_vol;

  cs_real_t *restrict val = f->val;
  cs_real_t p_va = 0.;

# pragma omp parallel for  reduction(+:p_va)
  for (cs_lnum_t c_id = 0 ; c_id < n_cells ; c_id++) {
    p_va += cell_f_vol[c_id]*val[c_id];
  }

  cs_parall_sum(1, CS_DOUBLE, &p_va);
  p_va = p_va / tot_f_vol;

  cs_real_t shift = va - p_va;

# pragma omp parallel for  if (n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
    val[c_id] += shift;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize current parallel and periodic field values.
 *
 * This function currently only upates fields based on CS_MESH_LOCATION_CELLS.
 *
 * \param[in, out]   f           pointer to field
 * \param[in]        halo_type   halo type
 */
/*----------------------------------------------------------------------------*/

void
cs_field_synchronize(cs_field_t      *f,
                     cs_halo_type_t   halo_type)
{
  if (f->location_id == CS_MESH_LOCATION_CELLS) {

    const cs_halo_t *halo = cs_glob_mesh->halo;

    if (halo != NULL) {

      if (f->dim == 1)
        cs_halo_sync_var(halo, halo_type, f->val);

      else {

        cs_halo_sync_var_strided(halo, halo_type, f->val, f->dim);

        if (cs_glob_mesh->n_init_perio > 0) {
          switch(f->dim) {
          case 9:
            cs_halo_perio_sync_var_tens(halo, halo_type, f->val);
            break;
          case 6:
            cs_halo_perio_sync_var_sym_tens(halo, halo_type, f->val);
            break;
          case 3:
            cs_halo_perio_sync_var_vect(halo, halo_type, f->val, 3);
            break;
          default:
            break;
          }

        }

      }

    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
