#ifndef __CS_AT_OPT_INTERP_H__
#define __CS_AT_OPT_INTERP_H__

/*============================================================================
 * Optimal Interpolation.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_field.h"
#include "cs_measures_util.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef enum {

  CS_AT_OPT_INTERP_P0,    /* Interpolation from cell containing the observation */
  CS_AT_OPT_INTERP_P1     /* Interpolation from (partial) extended neighbourhood */

} cs_at_opt_interp_type_t;

typedef struct _cs_at_opt_interp_t {

  const char              *name;                /* Name */
  int                      id;                  /* Id */
  int                      ig_id;
  cs_real_t               *obs_cov;
  bool                     obs_cov_is_diag;
  cs_at_opt_interp_type_t  interp_type;
  cs_real_t               *model_to_obs_proj;
  cs_lnum_t               *model_to_obs_proj_idx;
  cs_lnum_t               *model_to_obs_proj_c_ids;
  cs_real_t               *b_proj;
  cs_real_t                ir[2];
  cs_real_t               *relax;
  int                      nb_times;
  int                     *measures_idx;
  cs_real_t               *times;
  cs_real_t               *times_read;
  int                     *active_time;
  cs_real_t               *time_weights;
  cs_real_t               *time_window;
  int                      n_log_data;
  int                      steady;
  int                      frequency;
  int                      type_nudging;

} cs_at_opt_interp_t;

/*============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an optimal interpolation descriptor.
 *
 * \param[in]  name   optimal interpolation name
 *
 */
/*----------------------------------------------------------------------------*/

cs_at_opt_interp_t *
cs_at_opt_interp_create(const char   *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to an optimal interpolation based on its id.
 *
 * This function requires that an optimal interpolation of the given id is
 * defined.
 *
 * \param[in]  id   optimal interpolation id
 *
 * \return  pointer to the optimal interpolation structure.
 */
/*----------------------------------------------------------------------------*/

cs_at_opt_interp_t  *
cs_at_opt_interp_by_id(int  id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return a pointer to an optimal interpolation based on its name.
 *
 * This function requires that an optimal interpolation of the given name is
 * defined.
 *
 * \param[in]  name   optimal interpolation name
 *
 * \return  pointer to the optimal interpolation structure.
 */
/*----------------------------------------------------------------------------*/

cs_at_opt_interp_t  *
cs_at_opt_interp_by_name(const char  *name);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy all defined optimal interpolations.
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interps_destroy(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read an optimal interpolation file for a given variable
 *        and fill in the matching measures set and optimal interpolation
 *        structures.
 *
 * \param[in]  filename   name of interpolation file
 * \param[in]  ms         measures set structure
 * \param[in]  oi         optimal interpolation structure
 * \param[in]  f_dim      dimension of field
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_read_file(char const           filename[50],
                           cs_measures_set_t   *ms,
                           cs_at_opt_interp_t  *oi,
                           const int            f_dim);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return 1 if a p1 projection has been enabled for at least one
 *        optimal interpolation. This function is used to determine if
 *        extended neighborhood is needed.
 *
 * \return  1 if a p1 proj. is needed, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_at_opt_interp_is_p1_proj_needed(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief (re)Allocate and fill in an optimal interpolation structure from an
 *        optimal interpolation file.
 *
 * \param[in]  oi  pointer to the optimal interpolation
 * \param[in]  ms  pointer to the associated measures set
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_map_values(cs_at_opt_interp_t *oi,
                            cs_measures_set_t  *ms);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute observation operator (H).
 *
 * \param[in]  ms  pointer to measures set
 * \param[in]  oi  pointer to an optimal interpolation
 * \param[in]  ig  pointer to interpol grid
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_obs_operator(cs_measures_set_t  *ms,
                              cs_at_opt_interp_t *oi,
                              cs_interpol_grid_t *ig);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute $\tens{H}\tens{B}\transpose{\tens{H}}$.
 *
 * \param[in]  ms  pointer to measures set
 * \param[in]  oi  pointer to an optimal interpolation
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_project_model_covariance(cs_measures_set_t  *ms,
                                          cs_at_opt_interp_t *oi);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Count active observations and compute time weights in case of
 *        unsteady.
 *
 * \param[in]  ms           pointer to measures set
 * \param[in]  oi           optimal interpolation for field variable
 * \param[in]  f_oia        analysis field of field variable
 * \param[in]  inverse      boolean, true if it necessary to recompute the
 *                          inverse of HB(H)
 * \param[in]  ao_idx       index of active observations
 *
 * \return  number of active observations.
 */
/*----------------------------------------------------------------------------*/

int *
cs_at_opt_interp_get_active_obs(cs_measures_set_t  *ms,
                                cs_at_opt_interp_t *oi,
                                cs_field_t         *f_oia,
                                bool              **inverse,
                                int              ***ao_idx);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute analysis for a given variable.
 *
 * \param[in]  f            field variable of which analysis will be computed
 * \param[in]  oi           optimal interpolation for field variable
 * \param[in]  f_oia        analysis field of field variable
 * \param[in]  n_active_obs number of active observations.
 * \param[in]  ao_idx       index of active observations
 * \param[in]  inverse      boolean, true if it necessary to recompute the
 *                          inverse of HB(H)
 */
/*----------------------------------------------------------------------------*/

void
cs_at_opt_interp_compute_analysis(cs_field_t         *f,
                                  cs_at_opt_interp_t *oi,
                                  cs_field_t         *f_oia,
                                  int                 n_active_obs,
                                  int                *ao_idx,
                                  bool                inverse,
                                  int                 mc_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_AT_OPT_INTERP_H__ */
