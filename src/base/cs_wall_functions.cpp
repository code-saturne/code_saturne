/*============================================================================
 * Wall functions
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"
#include "base/cs_log.h"
#include "base/cs_mem.h"
#include "base/cs_parall.h"
#include "base/cs_physical_constants.h"
#include "base/cs_thermal_model.h"
#include "base/cs_time_step.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_quantities.h"
#include "turb/cs_turbulence_model.h"
#include "cdo/cs_domain.h"
#include "base/cs_field.h"
#include "base/cs_field_operator.h"
#include "base/cs_field_pointer.h"
#include "base/cs_field_default.h"
#include "pprt/cs_physical_model.h"
#include "atmo/cs_atmo.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_wall_functions.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_wall_functions.cpp
        Wall functions descriptor and computation.
*/

/*! \struct cs_wall_functions_t

  \brief wall functions descriptor.

  Members of this wall functions descriptor are publicly accessible, to allow
  for concise syntax, as it is expected to be used in many places.

  \var  cs_wall_functions_t::iwallf
        Indicates the type of wall function used for the velocity
        boundary conditions on a frictional wall.\n
        - CS_WALL_F_DISABLED: no wall functions
        - CS_WALL_F_1SCALE_POWER: one scale of friction velocities (power law)
        - CS_WALL_F_1SCALE_LOG: one scale of friction velocities (log law)
        - CS_WALL_F_2SCALES_LOG: two scales of friction velocities (log law)
        - CS_WALL_F_SCALABLE_2SCALES_LOG: two scales of friction velocities
          (log law) (scalable wall functions)
        - CS_WALL_F_2SCALES_VDRIEST: two scales of friction velocities
          (mixing length based on V. Driest analysis)
        - CS_WALL_F_2SCALES_SMOOTH_ROUGH: wall function unifying rough and smooth
          friction regimes
        - CS_WALL_F_2SCALES_CONTINUOUS: All \f$ y^+ \f$  for low Reynolds models\n
        \ref iwallf is initialized to CS_WALL_F_1SCALE_LOG for \ref model = 10,
          40, 41 or 70
        (mixing length, LES and Spalart Allmaras).\n
        \ref iwallf is initialized to CS_WALL_F_DISABLED for \ref model = 0, 32,
          50 or 51\n
        \ref iwallf is initialized to CS_WALL_F_2SCALES_LOG for \ref model = 20,
          21, 30, 31 or 60
        (\f$k-\epsilon\f$, \f$R_{ij}-\epsilon\f$ LRR, \f$R_{ij}-\epsilon\f$ SSG
          and \f$k-\omega\f$ SST models).\n
        The v2f model (\ref model=50) is not designed to use wall functions
        (the mesh must be low Reynolds).\n
        The value \ref iwallf = CS_WALL_F_2SCALES_LOG is not compatible with
          \ref model=0, 10, 40
        or 41 (laminar, mixing length and LES).\n
        Concerning the \f$k-\epsilon\f$ and \f$R_{ij}-\epsilon\f$ models, the
        two-scales model is usually at least as satisfactory as the one-scale
        model.\n
        The scalable wall function allows to virtually shift the wall when
        necessary in order to be always in a logarithmic layer. It is used to make up for
        the problems related to the use of High-Reynolds models on very refined
        meshes.\n
        Useful if \ref model is different from 50.
  \var  cs_wall_functions_t::iwalfs
        wall functions for scalar
        - CS_WALL_F_S_ARPACI_LARSEN: three layers (Arpaci and Larsen) or two layers (Prandtl-Taylor) for
             Prandtl number smaller than 0.1
        - CS_WALL_F_S_VDRIEST: consistant with the 2 scales wall function for velocity using Van
             Driest mixing length
        - CS_WALL_F_S_LOUIS: default wall function for atmospheric flows for
             potential temperature. This has influence on the dynamic.
        - CS_WALL_F_S_MONIN_OBUKHOV: Monin Obukhov wall function for atmospheric flows for
             potential temperature. This has influence on the dynamic.

  \var  cs_wall_functions_t::ypluli
        limit value of \f$y^+\f$ for the viscous sublayer

        \ref ypluli depends on the chosen wall function: it is initialized to
        10.88 for the scalable wall function (\ref iwallf=CS_WALL_F_SCALABLE_2SCALES_LOG), otherwise it is
        initialized to \f$1/\kappa\approx 2,38\f$. In LES, \ref ypluli is taken
        by default to be 10.88. Always useful.

*/
/*----------------------------------------------------------------------------*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Wall functions structure and associated pointers

   ypluli is set to -grand*10. If the user has not changed this value,
   its value is modified in cs_parameters_*_complete (10.88 with invariant wall
   laws, 1/kappa otherwise). */

static cs_wall_functions_t  _wall_functions =
{
  .iwallf = CS_WALL_F_UNSET,
  .iwalfs = CS_WALL_F_S_UNSET,
  .ypluli = -1e13
};

const cs_wall_functions_t  * cs_glob_wall_functions = &_wall_functions;

static cs_wall_functions_iwf_options_t  _iwf_options = {
  .ri_min = 1.,
  .y_cross_min = 1.,
  .tolerance = 1.e-4,
  .max_iterations = 1000,
  .first_time_step = 4
};

static const cs_mesh_t                    *_iwf_mesh = nullptr;
static cs_lnum_t                           _iwf_storage_n_b_faces = 0;
static cs_lnum_t                           _iwf_n_configured_faces = 0;
static cs_real_t                          *_iwf_length = nullptr;
static cs_wall_functions_iwf_result_t     *_iwf_result = nullptr;
static bool                                _iwf_finalize_registered = false;

/*----------------------------------------------------------------------------*/
/*!
 * \brief Release all storage owned by the integral wall function extension.
 *
 * This code_saturne integration helper is not a physical equation from the
 * reference below. It is registered once through \ref cs_base_at_finalize and
 * owns the cleanup of both face arrays. It also clears the borrowed mesh
 * reference and resets the associated sizes and counters.
 *
 * \par IWF model implemented by this storage
 * Borel, P., Moulouel, R., Houvin, M., Ferrand, M., Benhamadouche, S., 2026.
 * Indoor turbulent airflow modelling with a new integral buoyant wall
 * function. International Journal of Heat and Fluid Flow 117, 110130.
 * https://doi.org/10.1016/j.ijheatfluidflow.2025.110130
 */
/*----------------------------------------------------------------------------*/

static void
_iwf_finalize(void)
{
  CS_FREE(_iwf_length);
  CS_FREE(_iwf_result);
  _iwf_mesh = nullptr;
  _iwf_storage_n_b_faces = 0;
  _iwf_n_configured_faces = 0;
  _iwf_finalize_registered = false;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Ensure that the per-boundary-face IWF storage matches the mesh.
 *
 * The existing arrays are reused while the mesh address and boundary-face
 * count remain unchanged. Otherwise they are released, reallocated, and fully
 * initialized. The function registers \ref _iwf_finalize only after storage
 * has been allocated successfully and never more than once at a time.
 *
 * This code_saturne lifecycle helper is an integration choice, not a physical
 * equation from the reference below.
 *
 * \par IWF model implemented by this storage
 * Borel, P., Moulouel, R., Houvin, M., Ferrand, M., Benhamadouche, S., 2026.
 * Indoor turbulent airflow modelling with a new integral buoyant wall
 * function. International Journal of Heat and Fluid Flow 117, 110130.
 * https://doi.org/10.1016/j.ijheatfluidflow.2025.110130
 */
/*----------------------------------------------------------------------------*/

static void
_iwf_ensure_storage(void)
{
  const cs_mesh_t *mesh = cs_glob_mesh;

  if (mesh == nullptr)
    return;

  if (   _iwf_mesh == mesh
      && _iwf_storage_n_b_faces == mesh->n_b_faces)
    return;

  CS_FREE(_iwf_length);
  CS_FREE(_iwf_result);

  _iwf_mesh = mesh;
  _iwf_storage_n_b_faces = mesh->n_b_faces;
  _iwf_n_configured_faces = 0;
  CS_MALLOC(_iwf_length, _iwf_storage_n_b_faces, cs_real_t);
  CS_MALLOC(_iwf_result,
            _iwf_storage_n_b_faces,
            cs_wall_functions_iwf_result_t);

  for (cs_lnum_t f_id = 0; f_id < _iwf_storage_n_b_faces; f_id++) {
    _iwf_length[f_id] = 0.;
    _iwf_result[f_id] = {};
  }

  if (!_iwf_finalize_registered) {
    cs_base_at_finalize(_iwf_finalize);
    _iwf_finalize_registered = true;
  }
}

/*----------------------------------------------------------------------------*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *! \brief Provide access to cs_glob_wall_functions
 *
 * needed to initialize structure with GUI
 *----------------------------------------------------------------------------*/

cs_wall_functions_t *
cs_get_glob_wall_functions(void)
{
  return &_wall_functions;
}

/*----------------------------------------------------------------------------*/

void
cs_wall_functions_iwf_set_faces(cs_lnum_t        n_faces,
                                const cs_lnum_t  face_ids[],
                                cs_real_t         characteristic_length)
{
  if (n_faces < 0 || (n_faces > 0 && face_ids == nullptr))
    bft_error(__FILE__, __LINE__, 0, "Invalid IWF face selection.");

  if (n_faces == 0)
    return;

  _iwf_ensure_storage();

  const bool configure = characteristic_length > 0.;

  for (cs_lnum_t i = 0; i < n_faces; i++) {
    const cs_lnum_t f_id = face_ids[i];
    if (f_id < 0 || f_id >= _iwf_storage_n_b_faces)
      bft_error(__FILE__, __LINE__, 0,
                "Invalid IWF boundary face id: %ld.", (long)f_id);

    const bool was_configured = _iwf_length[f_id] > 0.;
    if (configure && !was_configured)
      _iwf_n_configured_faces += 1;
    else if (!configure && was_configured) {
      assert(_iwf_n_configured_faces > 0);
      _iwf_n_configured_faces -= 1;
    }

    _iwf_length[f_id] = characteristic_length;
    _iwf_result[f_id] = {};
  }
}

/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_wall_functions_iwf_n_faces(void)
{
  if (   cs_glob_mesh == nullptr
      || _iwf_mesh != cs_glob_mesh
      || _iwf_storage_n_b_faces != cs_glob_mesh->n_b_faces) {
    _iwf_n_configured_faces = 0;
    return 0;
  }

  return _iwf_n_configured_faces;
}

/*----------------------------------------------------------------------------*/

void
cs_wall_functions_iwf_set_options(
  const cs_wall_functions_iwf_options_t  *options)
{
  if (   options == nullptr
      || options->ri_min < 0.
      || options->y_cross_min < 0.
      || options->tolerance <= 0.
      || options->max_iterations < 1
      || options->first_time_step < 0)
    bft_error(__FILE__, __LINE__, 0, "Invalid integral wall function options.");

  _iwf_options = *options;
}

/*----------------------------------------------------------------------------*/

void
cs_wall_functions_iwf_reset_faces(void)
{
  if (_iwf_length == nullptr) {
    _iwf_n_configured_faces = 0;
    return;
  }

  _iwf_ensure_storage();

  for (cs_lnum_t f_id = 0; f_id < _iwf_storage_n_b_faces; f_id++) {
    _iwf_length[f_id] = 0.;
    _iwf_result[f_id] = {};
  }

  _iwf_n_configured_faces = 0;
}

/*----------------------------------------------------------------------------*/

bool
cs_wall_functions_iwf_get_result(
  cs_lnum_t                       face_id,
  cs_wall_functions_iwf_result_t *result)
{
  if (   result == nullptr
      || cs_glob_mesh == nullptr
      || face_id < 0
      || face_id >= cs_glob_mesh->n_b_faces)
    return false;

  if (cs_wall_functions_iwf_n_faces() == 0) {
    *result = {};
    return true;
  }

  _iwf_ensure_storage();

  *result = _iwf_result[face_id];
  return true;
}

/*----------------------------------------------------------------------------*/

void
cs_wall_functions_iwf_compute(const cs_real_t  wall_temperature[])
{
  _iwf_ensure_storage();

  cs_gnum_t counts[7]
    = {static_cast<cs_gnum_t>(_iwf_n_configured_faces), 0, 0, 0, 0, 0, 0};
  cs_gnum_t invalid_counts[4] = {0, 0, 0, 0};
  cs_real_t iwf_sum[4] = {0., 0., 0., 0.};
  cs_real_t iwf_min[4] = {cs_math_big_r, cs_math_big_r,
                          cs_math_big_r, cs_math_big_r};
  cs_real_t iwf_max[4] = {-cs_math_big_r, -cs_math_big_r,
                          -cs_math_big_r, -cs_math_big_r};
  for (cs_lnum_t f_id = 0; f_id < _iwf_storage_n_b_faces; f_id++)
    _iwf_result[f_id] = {};

  cs_gnum_t n_configured_faces = counts[0];
  cs_parall_counter(&n_configured_faces, 1);
  counts[0] = n_configured_faces;

  if (n_configured_faces == 0)
    return;

  const cs_wall_f_type_t iwallf = cs_glob_wall_functions->iwallf;
  bool supported_wall_function = false;

  switch (iwallf) {
  case CS_WALL_F_1SCALE_POWER:
  case CS_WALL_F_1SCALE_LOG:
  case CS_WALL_F_2SCALES_LOG:
  case CS_WALL_F_SCALABLE_2SCALES_LOG:
  case CS_WALL_F_2SCALES_VDRIEST:
  case CS_WALL_F_2SCALES_CONTINUOUS:
    supported_wall_function = true;
    break;
  case CS_WALL_F_UNSET:
  case CS_WALL_F_DISABLED:
  case CS_WALL_F_2SCALES_SMOOTH_ROUGH:
  default:
    break;
  }

  if (!supported_wall_function)
    bft_error(__FILE__, __LINE__, 0,
              "Integral wall functions require an enabled, non-rough "
              "velocity wall function (supported values: 1, 2, 3, 4, 5, "
              "and 7); current value is %d.", (int)iwallf);

  const cs_lnum_t nt_cur = (cs_glob_time_step != nullptr)
                         ? cs_glob_time_step->nt_cur : 0;
  if (nt_cur < _iwf_options.first_time_step) {
    cs_log_printf(CS_LOG_DEFAULT,
                  "IWF step %ld: configured=%llu, waiting for step %ld.\n",
                  (long)nt_cur,
                  (unsigned long long)n_configured_faces,
                  (long)_iwf_options.first_time_step);
    return;
  }

  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_nreal_3_t *b_face_u_normal = mq->b_face_u_normal;
  const cs_real_t *b_dist = mq->b_dist;

  const cs_field_t *temperature = cs_field_by_name_try("temperature");
  const cs_field_t *b_temperature
    = cs_field_by_name_try("boundary_temperature");
  const cs_field_t *conductivity
    = cs_field_by_name_try("thermal_conductivity");
  const cs_field_t *expansion
    = cs_field_by_name_try("thermal_expansion");
  const cs_field_t *density = cs_field_by_name_try("density");
  const cs_field_t *viscosity
    = cs_field_by_name_try("molecular_viscosity");
  const cs_field_t *velocity = CS_F_(vel);

  const cs_real_t *gravity = cs_get_glob_physical_constants()->gravity;
  const cs_real_t g_norm = cs_math_3_norm(gravity);
  const cs_real_3_t *vel = (velocity != nullptr)
                        ? (const cs_real_3_t *)velocity->val : nullptr;
  const cs_field_t *specific_heat = cs_field_by_name_try("specific_heat");
  const cs_fluid_properties_t *fluid_properties
    = cs_get_glob_fluid_properties();

  for (cs_lnum_t f_id = 0; f_id < _iwf_storage_n_b_faces; f_id++) {
    const cs_real_t length = _iwf_length[f_id];
    if (length <= 0.)
      continue;

    cs_wall_functions_iwf_result_t &r = _iwf_result[f_id];
    const cs_lnum_t c_id = b_face_cells[f_id];

    if (mq->has_disable_flag && mq->c_disable_flag[c_id]) {
      counts[1] += 1;
      continue;
    }

    if (   temperature == nullptr
        || (b_temperature == nullptr && wall_temperature == nullptr)
        || velocity == nullptr
        || g_norm <= 0. || b_dist[f_id] <= 0.) {
      counts[2] += 1;
      invalid_counts[0] += 1;
      continue;
    }

    const cs_real_t t_cell = temperature->val[c_id];
    cs_real_t t_wall = (b_temperature != nullptr)
                     ? b_temperature->val[f_id] : wall_temperature[f_id];
    if (   wall_temperature != nullptr
        && isfinite(wall_temperature[f_id])
        && cs::abs(wall_temperature[f_id]) < 1.e6)
      t_wall = wall_temperature[f_id];
    const cs_real_t delta_cell = cs::abs(t_wall - t_cell);
    const cs_real_t rho = (density != nullptr)
                        ? density->val[c_id] : fluid_properties->ro0;
    const cs_real_t mu = (viscosity != nullptr)
                       ? viscosity->val[c_id] : fluid_properties->viscl0;
    const cs_real_t lambda = (conductivity != nullptr)
                           ? conductivity->val[c_id]
                           : fluid_properties->lambda0;
    const cs_real_t temperature_offset
      = (cs_glob_thermal_model->temperature_scale
         == CS_TEMPERATURE_SCALE_CELSIUS) ? 273.15 : 0.;
    const cs_real_t mean_absolute_temperature
      = 0.5 * (t_cell + t_wall) + temperature_offset;
    cs_real_t beta = -1.;
    if (expansion != nullptr)
      beta = expansion->val[c_id];
    if (beta <= 0. && mean_absolute_temperature > 0.)
      beta = 1. / mean_absolute_temperature;
    const cs_real_t cp_loc = (specific_heat != nullptr)
                           ? specific_heat->val[c_id] : fluid_properties->cp0;

    if (delta_cell <= 1.e-4) {
      counts[2] += 1;
      invalid_counts[1] += 1;
      continue;
    }

    if (   !isfinite(t_cell) || !isfinite(t_wall) || !isfinite(rho)
        || !isfinite(mu) || !isfinite(lambda) || !isfinite(beta)
        || !isfinite(cp_loc) || rho <= 0. || mu <= 0. || lambda <= 0.
        || beta <= 0. || cp_loc <= 0.) {
      counts[2] += 1;
      invalid_counts[2] += 1;
      continue;
    }

    const cs_real_t nu = mu / rho;
    const cs_real_t prandtl = mu * cp_loc / lambda;
    const cs_real_t cos_angle
      = cs_math_3_dot_product(gravity, b_face_u_normal[f_id]) / g_norm;

    const int face_type = _iwf_face_type(cos_angle, t_wall - t_cell);

    cs_real_t t_infinity = t_cell;
    cs_real_t nusselt = 0.;
    cs_real_t y_cross = 0.;
    bool valid = true;

    for (int iter = 1; iter <= _iwf_options.max_iterations; iter++) {
      const cs_real_t delta_t = cs::abs(t_infinity - t_wall);
      const cs_real_t grashof
        = g_norm * beta * delta_t * pow(length, 3.) / cs_math_pow2(nu);
      const cs_real_t rayleigh = grashof * prandtl;

      if (   !isfinite(rayleigh) || rayleigh < 0.
          || !isfinite(prandtl) || prandtl <= 0.) {
        valid = false;
        break;
      }

      nusselt = _iwf_nusselt(rayleigh, prandtl, face_type);
      y_cross = b_dist[f_id] / length * nusselt;

      cs_real_t theta = 1.;
      if (y_cross > 1.e-4) {
        const cs_real_t profile_scale = 0.59 * pow(prandtl, 0.25);
        theta = profile_scale / 0.6
              * (1. - exp(-1.2 * y_cross / profile_scale))
              / (2. * y_cross);
      }

      cs_real_t t_new = t_wall;
      if (theta <= 0.99)
        t_new = (theta * t_wall - t_cell) / (theta - 1.);

      r.n_iterations = iter;
      if (!isfinite(t_new)) {
        valid = false;
        break;
      }
      if (cs::abs(t_new - t_infinity) <= _iwf_options.tolerance) {
        t_infinity = t_new;
        r.converged = true;
        break;
      }
      t_infinity = t_new;
    }

    if (!valid) {
      counts[2] += 1;
      invalid_counts[3] += 1;
      continue;
    }
    if (!r.converged) {
      counts[3] += 1;
      continue;
    }

    const cs_real_t delta_t = cs::abs(t_infinity - t_wall);
    const cs_real_t grashof
      = g_norm * beta * delta_t * pow(length, 3.) / cs_math_pow2(nu);
    const cs_real_t rayleigh = grashof * prandtl;
    nusselt = _iwf_nusselt(rayleigh, prandtl, face_type);
    y_cross = b_dist[f_id] / length * nusselt;

    const cs_real_t speed = cs_math_3_norm(vel[c_id]);
    const cs_real_t ri
      = _iwf_richardson(g_norm, beta, delta_t, length, speed);
    const cs_real_t u_ref = 0.5 * sqrt(g_norm * beta * delta_t * length);
    const cs_real_t reynolds = cs::max(u_ref * length / nu, 1.);
    const cs_real_t cf
      = 0.0776 / cs_math_pow2(log10(reynolds) - 1.88) + 60. / reynolds;
    const cs_real_t u_tau = u_ref * sqrt(0.5 * cf);
    const cs_real_t h_cv = nusselt * lambda / length;
    const cs_real_t h_int
      = h_cv * (t_wall - t_infinity) / (t_wall - t_cell);

    r.h_cv = h_cv;
    r.h_int = h_int;
    r.t_infinity = t_infinity;
    r.ri = ri;
    r.y_cross = y_cross;
    r.u_tau = u_tau;
    r.u_ref = u_ref;

    if (   !isfinite(h_cv) || !isfinite(h_int) || !isfinite(ri)
        || !isfinite(y_cross) || !isfinite(u_tau) || !isfinite(u_ref)
        || h_cv <= 0. || h_int <= 0. || u_tau <= 0.) {
      r.converged = false;
      counts[2] += 1;
      invalid_counts[3] += 1;
      continue;
    }
    if (ri < _iwf_options.ri_min) {
      counts[4] += 1;
      continue;
    }
    if (y_cross < _iwf_options.y_cross_min) {
      counts[5] += 1;
      continue;
    }

    r.active = true;
    counts[6] += 1;
    const cs_real_t values[4] = {h_cv, h_int, u_tau, u_ref};
    for (int i = 0; i < 4; i++) {
      iwf_sum[i] += values[i];
      iwf_min[i] = cs::min(iwf_min[i], values[i]);
      iwf_max[i] = cs::max(iwf_max[i], values[i]);
    }
  }

  cs_parall_counter(counts + 1, 6);
  cs_parall_counter(invalid_counts, 4);
  cs_parall_sum(4, CS_REAL_TYPE, iwf_sum);
  cs_parall_min(4, CS_REAL_TYPE, iwf_min);
  cs_parall_max(4, CS_REAL_TYPE, iwf_max);
  cs_log_printf(CS_LOG_DEFAULT,
                "IWF step %ld: configured=%llu active=%llu solid=%llu "
                "invalid=%llu non_converged=%llu Ri_rejected=%llu "
                "y_cross_rejected=%llu.\n",
                (long)nt_cur,
                (unsigned long long)counts[0],
                (unsigned long long)counts[6],
                (unsigned long long)counts[1],
                (unsigned long long)counts[2],
                (unsigned long long)counts[3],
                (unsigned long long)counts[4],
                (unsigned long long)counts[5]);
  cs_log_printf(CS_LOG_DEFAULT,
                "IWF invalid detail: missing=%llu delta_T=%llu "
                "properties=%llu computed=%llu.\n",
                (unsigned long long)invalid_counts[0],
                (unsigned long long)invalid_counts[1],
                (unsigned long long)invalid_counts[2],
                (unsigned long long)invalid_counts[3]);
  if (counts[6] > 0) {
    for (int i = 0; i < 4; i++)
      iwf_sum[i] /= counts[6];
    cs_log_printf(CS_LOG_DEFAULT,
                  "IWF statistics: h_cv=[%g %g %g] h_int=[%g %g %g] "
                  "u_tau=[%g %g %g] u_ref=[%g %g %g].\n",
                  iwf_min[0], iwf_sum[0], iwf_max[0],
                  iwf_min[1], iwf_sum[1], iwf_max[1],
                  iwf_min[2], iwf_sum[2], iwf_max[2],
                  iwf_min[3], iwf_sum[3], iwf_max[3]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute friction velocity u* and surface sensible heat flux q0
 * for a non neutral atmospheric surface layer using the explicit formula
 * developed for the ECMWF by Louis (1982)

 * \param[in]     utau        tangential mean velocity
 * \param[in]     rough_d     roughness z0
 * \param[in]     distbf      distance between the cell center and
                                  the center of gravity of the border face
 * \param[in]     duplus      1 over dimensionless velocity in neutral
 *                            conditions
 * \param[in]     dtplus      1 over dimensionless temperature in neutral
 *                            conditions
 * \param[in]     prt         turbulent Prandtl number
 * \param[in]     yplus_t     thermal dimensionless wall distance
 * \param[in]     buoyant_param g beta delta Theta
 * \param[in]     delta_t     delta Theta
 * \param[out]    uet         friction velocity
 * \param[out]    cfnns       non neutral correction coefficients for profiles
                              of scalar
 * \param[out]    cfnnk       non neutral correction coefficients
                              for profiles of k
 * \param[out]    cfnne       non neutral correction coefficients
                              for profiles of eps
 * \param[out]    dlmo        inverse Monin Obukhov length (for log only)
 * \param[in]     temp        potential temperature in boundary cell
 * \param[in]     totwt       total water content in boundary cell
 * \param[in]     liqwt       liquid water content in boundary cell
 */
/*----------------------------------------------------------------------------*/

static void
_atmo_louis(const cs_real_t  utau,
            const cs_real_t  rough_d,
            const cs_real_t  distbf,
            const cs_real_t  duplus,
            const cs_real_t  dtplus,
            const cs_real_t  prt,
            const cs_real_t  yplus_t,
            const cs_real_t  buoyant_param,
            const cs_real_t  delta_t,
            cs_real_t       *uet,
            cs_real_t       *cfnns,
            cs_real_t       *cfnnk,
            cs_real_t       *cfnne,
            cs_real_t       *dlmo)
{
  const cs_lnum_t nt_cur = cs_glob_time_step->nt_cur;

  /* Initializations
  --------------- */

  const cs_real_t b = 5.0;
  const cs_real_t c = b; // Already corrected with a equivalent Prantl number equal to 0.71
  const cs_real_t d = b;

  /* Compute layer average Richardson number */

  /* NB: rib = 0 if thermal flux conditions are imposed and tpot1 not defined */
  /* Patch for the initial time step when thermal field is not initalized */
  cs_real_t rib;
  if (cs::abs(utau) < cs_math_epzero || nt_cur == 1)
    rib = 0.0;
  else
    rib = buoyant_param * delta_t * distbf/utau/utau;

  /* Compute correction factors based on ECMWF parametrisation
     Louis (1982) */

  cs_real_t fm, fh, fmden1, fmden2, fhden;
  if (rib >= cs_math_epzero) {
    /* Stable case */
    fm = 1.0 / (1.0 + 2.0*b*rib/sqrt(1.0 + d*rib));
    fh = 1.0 / (1.0 + 3.0*b*rib*sqrt(1.0 + d*rib));
  }
  else {
    /* Unstable case */
    fmden1 = (yplus_t + 1.0) * cs::abs(rib);
    fmden2 = 1.0 + 3.0 * b * c * duplus * prt * dtplus * sqrt(fmden1);
    fm = 1.0 - 2.0 * b * rib / fmden2;
    fhden = 3.0 * b * c * duplus * prt * dtplus * sqrt(yplus_t + 1.0);
    fh = 1.0 - (3.0*b*rib)/(1.0 + fhden * sqrt(cs::abs(rib)));
  }

  if (fm <= cs_math_epzero)
    fm = cs_math_epzero;

  if (cs::abs(fh) <= cs_math_epzero)
    fh = cs_math_epzero;

  /* Stable */
  if ((1.0-rib) > cs_math_epzero) {
    *cfnnk = sqrt(1.0 - rib); /* TODO write it +correction with turbulent Prandtl */
    *cfnne = (1.0 - rib) / sqrt(fm);
  }
  else {
    *cfnnk = 1.0;
    *cfnne = 1.0;
  }

  /* Note: non neutral correction coefficients for profiles of wind
     (Re) compute friction velocity uet (for non neutral)
     uet = U/U^+ = U / U^{+,n} * sqrt(fm) */
  *uet = duplus * utau * sqrt(fm);

  /* Compute surface sensible heat flux q0 (can be useful for post-processing)
     Note: non-consistent with two velocity scales */
  *cfnns = fh / sqrt(fm);
  //Note:  q0 = (tpot1 - tpot2) * (*uet) * dtplus * (*cfnns);

  /* Compute local Obukhov inverse length for log
     1/L =  Ri / (z Phim) */
  *dlmo = rib * sqrt(fm) / (distbf + rough_d);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the friction velocity and \f$y^+\f$ / \f$u^+\f$.
 *
 * \param[in]     iwallf        wall function type
 * \param[in]     l_visc        kinematic viscosity
 * \param[in]     t_visc        turbulent kinematic viscosity
 * \param[in]     vel           wall projected cell center velocity
 * \param[in]     y             wall distance
 * \param[in]     rough_d       roughness length scale
 *                              (not sand grain roughness)
 * \param[in]     rnnb          \f$\vec{n}.(\tens{R}\vec{n})\f$
 * \param[in]     kinetic_en    turbulent kinetic energy (cell center)
 * \param[in]     rough_t       boundary_thermal_roughness
 * \param[in]     buoyant_param Buoyant parameter: g beta delta Theta v
 * \param[in]     delta_t       temperature delta
 * \param[in]     turb_prandtl  turbulent Prandtl
 * \param[in]     icodcl_th_fid  icodcl_th_f[f_id]]
 * \param[in]     flux
 * \param[in,out] iuntur        indicator: 0 in the viscous sublayer
 * \param[in,out] nsubla        counter of cell in the viscous sublayer
 * \param[in,out] nlogla        counter of cell in the log-layer
 * \param[out]    cfbnns        corrective factors
 * \param[out]    cfnnk
 * \param[out]    cfnne
 * \param[out]    dlmo          inverse of Obukhov length
 * \param[out]    ustar         friction velocity
 * \param[out]    uk            friction velocity
 * \param[out]    yplus         dimensionless distance to the wall
 * \param[out]    ypup          yplus projected vel ratio
 * \param[out]    cofimp        \f$\frac{|U_F|}{|U_I^p|}\f$ to ensure a good
 *                              turbulence production
 * \param[out]    dplus         dimensionless shift to the wall for scalable
 *                              wall functions
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_functions_velocity(cs_wall_f_type_t  iwallf,
                           cs_real_t         l_visc,
                           cs_real_t         t_visc,
                           cs_real_t         vel,
                           cs_real_t         y,
                           cs_real_t         rough_d,
                           cs_real_t         rnnb,
                           cs_real_t         kinetic_en,
                           cs_real_t         rough_t,
                           cs_real_t         buoyant_param,
                           cs_real_t         delta_t,
                           cs_real_t         turb_prandtl,
                           int               icodcl_th_fid,
                           cs_real_t         flux,
                           cs_field_t       *f_th,
                           int              *iuntur,
                           cs_gnum_t        *nsubla,
                           cs_gnum_t        *nlogla,
                           cs_real_t        *ustar,
                           cs_real_t        *uk,
                           cs_real_t        *yplus,
                           cs_real_t        *ypup,
                           cs_real_t        *cofimp,
                           cs_real_t        *dplus,
                           cs_real_t        *cfnns,
                           cs_real_t        *cfnnk,
                           cs_real_t        *cfnne,
                           cs_real_t        *dlmo)
{
  cs_real_t lmk;
  bool wf = true;

  /* Pseudo shift of the wall, 0 by default */
  *dplus = 0.;

  /* Activation of wall function by default */
  *iuntur = 1;

  /* Sand Grain roughness */
  cs_real_t sg_rough = rough_d * exp(cs_turb_xkappa*cs_turb_cstlog_rough);

  switch (iwallf) {
  case CS_WALL_F_DISABLED:
    cs_wall_functions_disabled(l_visc,
                               t_visc,
                               vel,
                               y,
                               iuntur,
                               nsubla,
                               nlogla,
                               ustar,
                               uk,
                               yplus,
                               dplus,
                               ypup,
                               cofimp);
    break;
  case CS_WALL_F_1SCALE_POWER:
    cs_wall_functions_1scale_power(l_visc,
                                   vel,
                                   y,
                                   iuntur,
                                   nsubla,
                                   nlogla,
                                   ustar,
                                   uk,
                                   yplus,
                                   ypup,
                                   cofimp);
    break;
  case CS_WALL_F_1SCALE_LOG:
    cs_wall_functions_1scale_log(l_visc,
                                 vel,
                                 y,
                                 iuntur,
                                 nsubla,
                                 nlogla,
                                 ustar,
                                 uk,
                                 yplus,
                                 ypup,
                                 cofimp);
    break;
  case CS_WALL_F_2SCALES_LOG:
    cs_wall_functions_2scales_smooth_rough(l_visc,
                                           t_visc,
                                           vel,
                                           y,
                                           rough_d,
                                           kinetic_en,
                                           iuntur,
                                           nsubla,
                                           nlogla,
                                           ustar,
                                           uk,
                                           yplus,
                                           dplus,
                                           ypup,
                                           cofimp);

    break;
  case CS_WALL_F_SCALABLE_2SCALES_LOG:
    cs_wall_functions_2scales_scalable(l_visc,
                                       t_visc,
                                       vel,
                                       y,
                                       kinetic_en,
                                       iuntur,
                                       nsubla,
                                       nlogla,
                                       ustar,
                                       uk,
                                       yplus,
                                       dplus,
                                       ypup,
                                       cofimp);
    break;
  case CS_WALL_F_2SCALES_VDRIEST:
    cs_wall_functions_2scales_vdriest(rnnb,
                                      l_visc,
                                      vel,
                                      y,
                                      kinetic_en,
                                      iuntur,
                                      nsubla,
                                      nlogla,
                                      ustar,
                                      uk,
                                      yplus,
                                      ypup,
                                      cofimp,
                                      &lmk,
                                      sg_rough,
                                      wf);
    break;
  case CS_WALL_F_2SCALES_SMOOTH_ROUGH:
    cs_wall_functions_2scales_smooth_rough(l_visc,
                                           t_visc,
                                           vel,
                                           y,
                                           rough_d,
                                           kinetic_en,
                                           iuntur,
                                           nsubla,
                                           nlogla,
                                           ustar,
                                           uk,
                                           yplus,
                                           dplus,
                                           ypup,
                                           cofimp);
    break;
  case CS_WALL_F_2SCALES_CONTINUOUS:
    cs_wall_functions_2scales_continuous(rnnb,
                                         l_visc,
                                         t_visc,
                                         vel,
                                         y,
                                         kinetic_en,
                                         iuntur,
                                         nsubla,
                                         nlogla,
                                         ustar,
                                         uk,
                                         yplus,
                                         ypup,
                                         cofimp);
    break;
  default:
    break;
  }

  /* To be coherent with a wall function, clip it to 0 */
  *cofimp = cs::max(*cofimp, 0.);
  /* Louis or Monin Obukhov wall function for atmospheric flows */
  const cs_wall_f_s_type_t iwalfs = cs_glob_wall_functions->iwalfs;

  const cs_real_t xkappa = cs_turb_xkappa;

  // Buoyant correction (atmospheric flows)
  // Louis or Monin-Obukhov
  if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1
      && iwalfs != CS_WALL_F_S_MONIN_OBUKHOV) {
    cs_real_t yk =  *uk * y / l_visc;
    /* 1/U+ for neutral */
    const cs_real_t duplus = *ypup / yk;

    /* 1/T+
     * "y+_t" tends to "y/rough_t" for rough regime and to "y+k"
     * times a shift for smooth regime
     *
     * Rough regime reads:
     *   T+ = Prt/kappa ln(y/rough_t) = Prt * (ln(y/zeta)/kappa + 8.5)
     *      = Prt/kappa ln[y/zeta * exp(8.5 kappa)]
     *
     * Note zeta_t = rough_t * exp(8.5 kappa)
     *
     * Question: is 8.5 really in factor of Prt?
     *
     * Smooth regime reads:
     * T+ = Prt *(ln(y uk/nu)/kappa + 5.2)
     *    = Prt/kappa ln[y uk*exp(5.2 kappa) / nu]
     *
     * Mixed regime reads:
     *   T+ = Prt/kappa ln[y uk*exp(5.2 kappa)/(nu + alpha uk zeta)]
     *      = Prt/kappa ln[  y uk*exp(5.2 kappa)
     *                   / (nu + alpha uk rough_t * exp(8.5 kappa))]
     *      = Prt/kappa ln[  y uk*exp(5.2 kappa)
     *                   / (nu + alpha uk rough_t * exp(8.5 kappa))]
     * with
     *   alpha * exp(8.5 kappa) / exp(5.2 kappa) = 1
     * ie
     *   alpha = exp(-(8.5-5.2) kappa) = 0.25
     * so
     *   T+ = Prt/kappa ln[  y uk*exp(5.2 kappa)
     *                   / (nu + uk rough_t * exp(5.2 kappa))]
     *      = Prt/kappa ln[y+k / (exp(-5.2 kappa) + uk rough_t/nu)]
     */

    /* shifted y+ */
    /* FIXME use log constant */
    const cs_real_t yplus_t
      = yk / (exp(-xkappa * 5.2) + *uk *rough_t / l_visc);
    /* 1/T+ for neutral */
#if 0
    bft_printf("xkappa: %f", xkappa);
    bft_printf("yplus_t: %f", yplus_t);
    bft_printf("turb_prandtl: %f", turb_prandtl);
#endif
    const cs_real_t dtplus = xkappa / log(yplus_t) / turb_prandtl;

    _atmo_louis(vel,
                rough_d,
                y,
                duplus,
                dtplus,
                turb_prandtl,
                yplus_t,
                buoyant_param,
                delta_t,
                ustar,
                cfnns,
                cfnnk,
                cfnne,
                dlmo);

    /* Dimensionless velocity, recomputed and therefore may
       take stability into account */

    cs_real_t uplus = vel / *ustar;

    /* y+/U+ for non neutral is recomputed */
    *ypup = yk / cs::max(uplus, cs_math_epzero);

  }

  else if (iwalfs == CS_WALL_F_S_MONIN_OBUKHOV) {

    /* In case this is called with rough_d = 0  (ie z0 =0)
     * as we want U+ profile to get back to
     * ln(z/z0)/kappa = (ln(z uk/nu) + 5.2)/kappa
     * we clip z0 to be max between rough_d and
     * nu exp( -kappa*5.2) /uk
     */
    cs_real_t z0 = cs::max(rough_d, l_visc * exp(-xkappa * 5.2) / *uk);
    /* Compute local LMO */
    if (cs_glob_physical_model_flag[CS_ATMOSPHERIC] >= 1) {

      cs_mo_compute_from_thermal_diff(y,
                                      z0,
                                      vel,
                                      buoyant_param,
                                      delta_t,
                                      dlmo,
                                      ustar);

      /* Take stability into account for the turbulent velocity scale */
      cs_real_t coef_mom = cs_mo_phim(y + z0, *dlmo);
      const cs_real_t one_minus_ri
        = 1 - (y + z0) * *dlmo / coef_mom;
    }

    else if (icodcl_th_fid == 3) {

      cs_mo_compute_from_thermal_flux(y,
                                      z0,
                                      vel,
                                      buoyant_param,
                                      flux,
                                      dlmo,
                                      ustar);
    }
    else {

      /* No temperature delta: neutral */
      cs_mo_compute_from_thermal_diff(y,
                                      z0,
                                      vel,
                                      buoyant_param,
                                      delta_t,
                                      dlmo,
                                      ustar);
    }

    /* Take stability into account for the turbulent velocity scale */
    cs_real_t coef_mom = cs_mo_phim(y + z0, *dlmo);
    const cs_real_t one_minus_ri
      = 1 - (y + z0) * (*dlmo) / coef_mom;

    if (one_minus_ri > 0) {
      /* Warning: overwritting uk, yplus should be recomputed */
      *uk = *uk / pow(one_minus_ri, 0.25);
      *yplus = y * (*uk) / l_visc;

      /* Epsilon should be modified as well to get
         P+G = P(1-Ri) = epsilon
         P = -R_tn dU/dn = uk^2 uet Phi_m / (kappa z) */
      *cfnne = one_minus_ri * coef_mom;
      /* Nothing done for the moment for really high stability */
    }
    else {
      *cfnne = 1.;
    }
    /* Boundary condition on the velocity to have approximately
        the correct turbulence production */
    coef_mom = cs_mo_phim(y+z0, *dlmo);
    const cs_real_t coef_momm = cs_mo_phim(2 * y + z0, *dlmo);
    cs_real_t rcprod =   2*y*sqrt(  xkappa*(*uk)*coef_mom/t_visc
                                / (y+z0))
                     - coef_momm / (2.0 + z0 / y);

    *iuntur = 1;
    const cs_real_t _uplus = vel / *ustar;
    /* Coupled solving of the velocity components
        The boundary term for velocity gradient is implicit
        modified for non neutral boundary layer (in uplus) */
    *cofimp  = cs::min(cs::max(1-1/(xkappa*_uplus) * rcprod, 0),
                      1);
  } /* End Monin Obukhov */
}

/*----------------------------------------------------------------------------*/
/*!
 *  \brief Compute the correction of the exchange coefficient between the
 *         fluid and the wall for a turbulent flow.
 *
 *  This is function of the dimensionless
 *  distance to the wall \f$ y^+ = \dfrac{\centip \centf u_\star}{\nu}\f$.
 *
 *  Then the return coefficient reads:
 *  \f[
 *  h_{tur} = Pr \dfrac{y^+}{T^+}
 *  \f]
 *
 * \param[in]     iwalfs        type of wall functions for scalar
 * \param[in]     l_visc        kinematic viscosity
 * \param[in]     prl           laminar Prandtl number
 * \param[in]     prt           turbulent Prandtl number
 * \param[in]     rough_t       scalar roughness lenghth scale
 * \param[in]     uk            velocity scale based on TKE
 * \param[in]     yplus         dimensionless distance to the wall
 * \param[in]     dplus         dimensionless distance for scalable
 *                              wall functions
 * \param[out]    htur          corrected exchange coefficient
 * \param[out]    yplim         value of the limit for \f$ y^+ \f$
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_functions_scalar(cs_wall_f_s_type_t  iwalfs,
                         cs_real_t           l_visc,
                         cs_real_t           prl,
                         cs_real_t           prt,
                         cs_real_t           rough_t,
                         cs_real_t           uk,
                         cs_real_t           yplus,
                         cs_real_t           dplus,
                         cs_real_t          *htur,
                         cs_real_t          *yplim)
{
  switch (iwalfs) {
  case CS_WALL_F_S_ARPACI_LARSEN:
    cs_wall_functions_s_arpaci_larsen(l_visc,
                                      prl,
                                      prt,
                                      rough_t,
                                      uk,
                                      yplus,
                                      dplus,
                                      htur,
                                      yplim);
    break;
  case CS_WALL_F_S_VDRIEST:
    cs_wall_functions_s_vdriest(prl,
                                prt,
                                yplus,
                                htur);
    break;
  case CS_WALL_F_S_SMOOTH_ROUGH:
    cs_wall_functions_s_smooth_rough(l_visc,
                                     prl,
                                     prt,
                                     rough_t,
                                     uk,
                                     yplus,
                                     dplus,
                                     htur);
    break;
  default:
    /* TODO Monin Obukhov or Louis atmospheric wall function
     * must be adapted to smooth wall functions.
     * Arpaci and Larsen wall functions are put as in previous versions of
     * code_Saturne.
     * */
    cs_wall_functions_s_arpaci_larsen(l_visc,
                                      prl,
                                      prt,
                                      rough_t,
                                      uk,
                                      yplus,
                                      dplus,
                                      htur,
                                      yplim);
    break;
  }
}

/*----------------------------------------------------------------------------*/
