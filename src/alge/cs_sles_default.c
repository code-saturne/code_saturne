/*============================================================================
 * Sparse Linear Equation Solvers
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
#include "cs_field.h"
#include "cs_grid.h"
#include "cs_halo.h"
#include "cs_internal_coupling.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_matrix.h"
#include "cs_matrix_default.h"
#include "cs_matrix_util.h"
#include "cs_multigrid.h"
#include "cs_parameters.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_sles_pc.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_sles_default.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_SLES_DEFAULT_N_SETUPS 2  /* Number of concurrent setups allowed */

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

static int           _n_setups = 0;
static cs_sles_t    *_sles_setup[CS_SLES_DEFAULT_N_SETUPS];
static cs_matrix_t  *_matrix_setup[CS_SLES_DEFAULT_N_SETUPS][3];

static const int _poly_degree_default = 0;
static const int _n_max_iter_default = 10000;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Default definition of a sparse linear equation solver

  \param[in]  f_id         associated field id, or < 0
  \param[in]  name         associated name if f_id < 0, or NULL
  \param[in]  matrix_type  matrix type, if available, or CS_MATRIX_N_TYPES
                           if not determined at calling site
  \param[in]  symmetric    indicate if matrix is symmetric
*/
/*----------------------------------------------------------------------------*/

static void
_sles_default_native(int                f_id,
                     const char        *name,
                     cs_matrix_type_t   matrix_type,
                     bool               symmetric)
{
  int multigrid = 0;
  cs_sles_it_type_t sles_it_type = CS_SLES_N_IT_TYPES;
  int n_max_iter = _n_max_iter_default;

  if (name != NULL) {

    if (!strcmp(name, "wall_distance")) { /* distpr.f90 */
      sles_it_type = CS_SLES_PCG;
      multigrid = 1;
    }
    if (!strcmp(name, "yplus_wall")) { /* distyp.f90 */
      sles_it_type = CS_SLES_P_SYM_GAUSS_SEIDEL;
    }
    else if (   !strcmp(name, "hydrostatic_p")  /* calhyd.f90 */
             || !strcmp(name, "potential")) {   /* predfl.f90 */
      /* Copy from pressure if possible */
      cs_field_t *cvar_p = (cs_field_by_name_try("pressure"));
      cs_sles_t *src = NULL;
      if (cvar_p != NULL) {
        if (cvar_p->type & CS_FIELD_VARIABLE)
          src = cs_sles_find_or_add(cvar_p->id, NULL);
      }
      if (src != NULL) {
        cs_sles_t *dest = cs_sles_find_or_add(-1, name);
        if (cs_sles_copy(dest, src) == 0) /* Copy OK, we are done */
          return;
      }
      /* If copying from pressure failed, default to multigrid */
      sles_it_type = CS_SLES_PCG;
      multigrid = 1;
    }
    else if (!strcmp(name, "Prhydro")) { /* prehyd.f90 */
      sles_it_type = CS_SLES_PCG;
    }
    else if (!strcmp(name, "Pr compress")) { /* resopv.f90 */
      sles_it_type = CS_SLES_P_SYM_GAUSS_SEIDEL;
    }
    else if (!strcmp(name, "PoissonL")) { /* lageqp.f90 */
      sles_it_type = CS_SLES_PCG;
      n_max_iter = 1000;
    }
    else if (!strcmp(name, "radiation_p1")) { /* raypun.f90 */
      sles_it_type = CS_SLES_PCG;
      multigrid = 1;
    }
    /* cs_bad_cells_regularisation.c */
    else if (!strcmp(name, "potential_regularisation_scalar")) {
      sles_it_type = CS_SLES_PCG;
    }
    else if (!strcmp(name, "potential_regularisation_vector")) {
      sles_it_type = CS_SLES_PCG;
    }
    else if (!strcmp(name, "potential_regularisation_sym_tensor")) {
      sles_it_type = CS_SLES_PCG;
    }
  }
  else if (f_id > -1) {
    const cs_field_t *f = cs_field_by_id(f_id);
    if (!strcmp(f->name, "hydraulic_head")) {
      sles_it_type = CS_SLES_PCG;
      multigrid = 2;
    }
  }

  /* Final default */

  if (sles_it_type == CS_SLES_N_IT_TYPES) {

    int coupling_id = -1;

    if (f_id > -1) {
      const cs_field_t *f = cs_field_by_id(f_id);
      coupling_id
        = cs_field_get_key_int(f, cs_field_key_id("coupling_entity"));
    }

    if (symmetric) {
      sles_it_type = CS_SLES_PCG;
      if (f_id > -1 && coupling_id < 0)
        multigrid = 1;
    }
    else {
      if (coupling_id < 0)
        sles_it_type = CS_SLES_P_SYM_GAUSS_SEIDEL;
      else
        sles_it_type = CS_SLES_BICGSTAB;
    }
  }

  if (multigrid == 1) {

    /* Multigrid used as preconditioner if possible, as solver otherwise */

    if ((matrix_type == CS_MATRIX_MSR) || (matrix_type == CS_MATRIX_N_TYPES)) {
      if (sles_it_type == CS_SLES_PCG && cs_glob_n_threads > 1)
        sles_it_type = CS_SLES_FCG;
      cs_sles_it_t *c = cs_sles_it_define(f_id,
                                          name,
                                          sles_it_type,
                                          -1, /* poly_degree */
                                          n_max_iter);
      cs_sles_pc_t *pc = cs_multigrid_pc_create(CS_MULTIGRID_V_CYCLE);
      cs_sles_it_transfer_pc(c, &pc);
      cs_sles_t *sc = cs_sles_find(f_id, name);
      cs_sles_set_error_handler(sc, cs_sles_default_error);
    }
    else
      cs_multigrid_define(f_id, name, CS_MULTIGRID_V_CYCLE);

  }
  else if (multigrid == 2)
    cs_multigrid_define(f_id, name, CS_MULTIGRID_V_CYCLE);

  else
    (void)cs_sles_it_define(f_id,
                            name,
                            sles_it_type,
                            _poly_degree_default,
                            n_max_iter);

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Default definition of a sparse linear equation solver

  \param[in]  f_id  associated field id, or < 0
  \param[in]  name  associated name if f_id < 0, or NULL
  \param[in]  a     matrix
*/
/*----------------------------------------------------------------------------*/

void
cs_sles_default(int                 f_id,
                const char         *name,
                const cs_matrix_t  *a)
{
  cs_matrix_type_t type = cs_matrix_get_type(a);
  bool symmetric = cs_matrix_is_symmetric(a);

  _sles_default_native(f_id, name, type, symmetric);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Default setup for sparse linear equation solver API.
 *
 * This includes setup logging.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_default_setup(void)
{
  /* Associate "on the fly default definition" function */

  cs_sles_set_default_define(cs_sles_default);
  cs_sles_set_default_verbosity(cs_sles_default_get_verbosity);

  int key_cal_opt_id = cs_field_key_id("var_cal_opt");

  /* Define for all variable fields */

  const int n_fields = cs_field_n_fields();

  /* Define solver for all variable fields if not already done,
     based on convection/diffusion */

  if (key_cal_opt_id > -1) {

    for (int f_id = 0; f_id < n_fields; f_id++) {
      const cs_field_t *f = cs_field_by_id(f_id);
      if (f->type & CS_FIELD_VARIABLE) {

        void *context = NULL;
        cs_sles_t *sc = cs_sles_find(f->id, NULL);
        if (sc != NULL)
          context = cs_sles_get_context(sc);

        if (context == NULL) {
          /* Get the calculation option from the field */
          cs_var_cal_opt_t var_cal_opt;
          cs_field_get_key_struct(f, key_cal_opt_id, &var_cal_opt);
          bool symmetric = (var_cal_opt.iconv > 0) ? false : true;
          _sles_default_native(f_id, NULL, CS_MATRIX_N_TYPES, symmetric);
        }

      }
    }

  }

  /* Logging */

  cs_log_printf(CS_LOG_SETUP, "\n");
  cs_log_separator(CS_LOG_SETUP);

  if (cs_multigrid_needed() && cs_glob_n_ranks > 1)
    cs_grid_log_merge_options();

  cs_sles_it_log_parallel_options();

  cs_sles_log(CS_LOG_SETUP);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Default finalization for sparse linear equation solver API.
 *
 * This includes performance data logging output.
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_default_finalize(void)
{
  cs_sles_log(CS_LOG_PERFORMANCE);

  cs_multigrid_finalize();
  cs_sles_finalize();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return default verbosity associated to a field id, name couple.
 *
 * \param[in]  f_id  associated field id, or < 0
 * \param[in]  name  associated name if f_id < 0, or NULL
 *
 * \return  verbosity associated with field or name
 */
/*----------------------------------------------------------------------------*/

int
cs_sles_default_get_verbosity(int          f_id,
                              const char  *name)
{
  CS_UNUSED(name);

  int retval = 0;

  static int k_log = -1;
  static int k_cal_opt_id = -1;

  if (k_log < 0)
    k_log = cs_field_key_id("log");
  if (k_cal_opt_id < 0)
    k_cal_opt_id = cs_field_key_id("var_cal_opt");

  if (f_id > -1) {
    const cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      cs_var_cal_opt_t var_cal_opt;
      cs_field_get_key_struct(f, k_cal_opt_id, &var_cal_opt);
      retval = var_cal_opt.iwarni;
    }
    else
      retval = cs_field_get_key_int(f, k_log);
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call sparse linear equation solver setup for convection-diffusion
 *        systems
 *
 * \param[in]       f_id                   associated field id, or < 0
 * \param[in]       name                   associated name if f_id < 0, or NULL
 * \param[in]       diag_block_size        block sizes for diagonal, or NULL
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal,
 *                                         or NULL
 * \param[in]       da                     diagonal values (NULL if zero)
 * \param[in]       xa                     extradiagonal values (NULL if zero)
 * \param[in]       da_conv                diagonal values (NULL if zero)
 * \param[in]       xa_conv                extradiagonal values (NULL if zero)
 * \param[in]       da_diff                diagonal values (NULL if zero)
 * \param[in]       xa_diff                extradiagonal values (NULL if zero)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_setup_native_conv_diff(int                  f_id,
                               const char          *name,
                               const int           *diag_block_size,
                               const int           *extra_diag_block_size,
                               const cs_real_t     *da,
                               const cs_real_t     *xa,
                               const cs_real_t     *da_conv,
                               const cs_real_t     *xa_conv,
                               const cs_real_t     *da_diff,
                               const cs_real_t     *xa_diff)
{
  cs_matrix_t *a = NULL, *a_conv = NULL, *a_diff = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  /* Check if this system has already been setup */

  cs_sles_t *sc = cs_sles_find_or_add(f_id, name);

  int setup_id = 0;
  while (setup_id < _n_setups) {
    if (_sles_setup[setup_id] == sc)
      break;
    else
      setup_id++;
  }

  if (setup_id >= _n_setups) {

    _n_setups += 1;

    if (_n_setups > CS_SLES_DEFAULT_N_SETUPS)
      bft_error
        (__FILE__, __LINE__, 0,
         "Too many linear systems solved without calling cs_sles_free_native\n"
         "  maximum number of systems: %d\n"
         "If this is not an error, increase CS_SLES_DEFAULT_N_SETUPS\n"
         "  in file %s.", CS_SLES_DEFAULT_N_SETUPS, __FILE__);

    /* Set extended contribution for domain coupling */

    if (f_id > -1) {

      const cs_field_t *f = cs_field_by_id(f_id);
      int coupling_id
        = cs_field_get_key_int(f, cs_field_key_id("coupling_entity"));

      if (coupling_id > -1) {
        /* TODO base matrix type on solver type (for external libraries) */
        a = cs_matrix_set_coefficients_coupled(f,
                                               CS_MATRIX_MSR,
                                               false, /* symmetric */
                                               diag_block_size,
                                               extra_diag_block_size,
                                               da,
                                               xa);
      }

    }

    if (a == NULL) {

      a = cs_matrix_msr(false,
                        diag_block_size,
                        extra_diag_block_size);

      cs_matrix_set_coefficients(a,
                                 false,
                                 diag_block_size,
                                 extra_diag_block_size,
                                 m->n_i_faces,
                                 (const cs_lnum_2_t *)(m->i_face_cells),
                                 da,
                                 xa);

      cs_matrix_t *a_ref = cs_matrix_default(false,
                                             diag_block_size,
                                             extra_diag_block_size);
      if (a_conv == NULL)
        a_conv = cs_matrix_create_by_copy(a_ref);

      cs_matrix_set_coefficients(a_conv,
                                 false,
                                 diag_block_size,
                                 extra_diag_block_size,
                                 m->n_i_faces,
                                 (const cs_lnum_2_t *)(m->i_face_cells),
                                 da_conv,
                                 xa_conv);

      if (a_diff == NULL)
        a_diff = cs_matrix_create_by_copy(a_ref);

      cs_matrix_set_coefficients(a_diff,
                                 false,
                                 diag_block_size,
                                 extra_diag_block_size,
                                 m->n_i_faces,
                                 (const cs_lnum_2_t *)(m->i_face_cells),
                                 da_diff,
                                 xa_diff);

    }

    _sles_setup[setup_id] = sc;
    _matrix_setup[setup_id][0] = a;
    _matrix_setup[setup_id][1] = a_conv;
    _matrix_setup[setup_id][2] = a_diff;

  }
  else {
    a      = _matrix_setup[setup_id][0];
    a_conv = _matrix_setup[setup_id][1];
    a_diff = _matrix_setup[setup_id][2];
  }

  /* Setup system */

  if (strcmp(cs_sles_get_type(sc), "cs_multigrid_t") != 0)
    bft_error
      (__FILE__, __LINE__, 0,
       "%s requires a cs_multigrid_t solver type", __func__);

  int verbosity = cs_sles_get_verbosity(sc);

  cs_multigrid_t  *mg = cs_sles_get_context(sc);
  cs_multigrid_setup_conv_diff(mg, name, a, a_conv, a_diff, verbosity);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call sparse linear equation solver setup for systems with internal
 *        coupling.
 *
 * \param[in]  f_id                   associated field id, or < 0
 * \param[in]  name                   associated name if f_id < 0, or NULL
 * \param[in]  symmetric              indicates if matrix coefficients
 *                                    are symmetric
 * \param[in]  diag_block_size        block sizes for diagonal, or NULL
 * \param[in]  extra_diag_block_size  block sizes for extra diagonal,
 *                                    or NULL
 * \param[in]  da                     diagonal values (NULL if zero)
 * \param[in]  xa                     extradiagonal values (NULL if zero)
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_setup_native_coupling(int               f_id,
                              const char       *name,
                              bool              symmetric,
                              const int        *diag_block_size,
                              const int        *extra_diag_block_size,
                              const cs_real_t  *da,
                              const cs_real_t  *xa)
{
  cs_matrix_t *a = NULL;

  /* Check if this system has already been setup */

  cs_sles_t *sc = cs_sles_find_or_add(f_id, name);

  int setup_id = 0;
  while (setup_id < _n_setups) {
    if (_sles_setup[setup_id] == sc)
      break;
    else
      setup_id++;
  }

  if (setup_id >= _n_setups) {

    _n_setups += 1;

    if (_n_setups > CS_SLES_DEFAULT_N_SETUPS)
      bft_error
        (__FILE__, __LINE__, 0,
         "Too many linear systems solved without calling cs_sles_free_native\n"
         "  maximum number of systems: %d\n"
         "If this is not an error, increase CS_SLES_DEFAULT_N_SETUPS\n"
         "  in file %s.", CS_SLES_DEFAULT_N_SETUPS, __FILE__);

    /* Set extended contribution for domain coupling */

    if (f_id > -1) {

      const cs_field_t *f = cs_field_by_id(f_id);
      int coupling_id
        = cs_field_get_key_int(f, cs_field_key_id("coupling_entity"));

      if (coupling_id > -1) {
        /* TODO base matrix type on solver type (for external libraries) */
        a = cs_matrix_set_coefficients_coupled(f,
                                               CS_MATRIX_MSR,
                                               symmetric,
                                               diag_block_size,
                                               extra_diag_block_size,
                                               da,
                                               xa);
      }

    }

    _sles_setup[setup_id] = sc;
    _matrix_setup[setup_id][0] = a;
    _matrix_setup[setup_id][1] = a; /* so it is freed later */

  }
  else {
    a = _matrix_setup[setup_id][0];
  }

  /* Setup system */

  cs_sles_setup(sc, a);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call sparse linear equation solver using native matrix arrays.
 *
 * \param[in]       f_id                   associated field id, or < 0
 * \param[in]       name                   associated name if f_id < 0, or NULL
 * \param[in]       symmetric              indicates if matrix coefficients
 *                                         are symmetric
 * \param[in]       diag_block_size        block sizes for diagonal, or NULL
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal,
 *                                         or NULL
 * \param[in]       da                     diagonal values (NULL if zero)
 * \param[in]       xa                     extradiagonal values (NULL if zero)
 * \param[in]       rotation_mode          halo update option for
 *                                         rotational periodicity
 * \param[in]       precision              solver precision
 * \param[in]       r_norm                 residue normalization
 * \param[out]      n_iter                 number of "equivalent" iterations
 * \param[out]      residue                residue
 * \param[in]       rhs                    right hand side
 * \param[in, out]  vx                     system solution
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_solve_native(int                  f_id,
                     const char          *name,
                     bool                 symmetric,
                     const int           *diag_block_size,
                     const int           *extra_diag_block_size,
                     const cs_real_t     *da,
                     const cs_real_t     *xa,
                     cs_halo_rotation_t   rotation_mode,
                     double               precision,
                     double               r_norm,
                     int                 *n_iter,
                     double              *residue,
                     const cs_real_t     *rhs,
                     cs_real_t           *vx)
{
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;
  cs_matrix_t *a = NULL;

  const cs_mesh_t *m = cs_glob_mesh;

  bool need_msr = false;

  /* Check if this system has already been setup */

  cs_sles_t *sc = cs_sles_find_or_add(f_id, name);

  int setup_id = 0;
  while (setup_id < _n_setups) {
    if (_sles_setup[setup_id] == sc)
      break;
    else
      setup_id++;
  }

  if (setup_id >= _n_setups) {

    _n_setups += 1;

    if (_n_setups > CS_SLES_DEFAULT_N_SETUPS)
      bft_error
        (__FILE__, __LINE__, 0,
         "Too many linear systems solved without calling cs_sles_free_native\n"
         "  maximum number of systems: %d\n"
         "If this is not an error, increase CS_SLES_DEFAULT_N_SETUPS\n"
         "  in file %s.", CS_SLES_DEFAULT_N_SETUPS, __FILE__);

    /* If context has not been defined yet, temporarily set
       matrix coefficients (using native matrix, which has lowest
       overhead as coefficients are provided in that form)
       to define the required context.

       The matrix type might be modified later based on solver
       constraints. */

    if (cs_sles_get_context(sc) == NULL) {
      a = cs_matrix_native(symmetric,
                           diag_block_size,
                           extra_diag_block_size);

      cs_matrix_set_coefficients(a,
                                 symmetric,
                                 diag_block_size,
                                 extra_diag_block_size,
                                 m->n_i_faces,
                                 (const cs_lnum_2_t *)(m->i_face_cells),
                                 da,
                                 xa);

      cs_sles_define_t  *sles_default_func = cs_sles_get_default_define();
      sles_default_func(f_id, name, a);
      cs_matrix_release_coefficients(a);
    }

    assert(cs_sles_get_context(sc) != NULL);

    cs_sles_pc_t  *pc = NULL;
    cs_multigrid_t *mg = NULL;

    if (strcmp(cs_sles_get_type(sc), "cs_sles_it_t") == 0) {
      cs_sles_it_t *c = cs_sles_get_context(sc);
      cs_sles_it_type_t s_type = cs_sles_it_get_type(c);
      if (   s_type >= CS_SLES_P_GAUSS_SEIDEL
          && s_type <= CS_SLES_TS_B_GAUSS_SEIDEL)
        need_msr = true;
      else {
        pc = cs_sles_it_get_pc(c);
        if (pc != NULL) {
          if (strcmp(cs_sles_pc_get_type(pc), "multigrid") == 0)
            mg = cs_sles_pc_get_context(pc);
        }
      }
    }
    else if (strcmp(cs_sles_get_type(sc), "cs_multigrid_t") == 0)
      mg = cs_sles_get_context(sc);

    if (mg != NULL) {
      cs_sles_it_type_t fs_type = cs_multigrid_get_fine_solver_type(mg);
      if (   fs_type >= CS_SLES_P_GAUSS_SEIDEL
          && fs_type <= CS_SLES_TS_B_GAUSS_SEIDEL)
        need_msr = true;
    }

    /* MSR not supported yet for full blocks */
    if (extra_diag_block_size != NULL) {
      if (extra_diag_block_size[0] > 1)
        need_msr = false;
    }

    if (need_msr)
      a = cs_matrix_msr(symmetric,
                        diag_block_size,
                        extra_diag_block_size);
    else
      a = cs_matrix_default(symmetric,
                            diag_block_size,
                            extra_diag_block_size);

    cs_matrix_set_coefficients(a,
                               symmetric,
                               diag_block_size,
                               extra_diag_block_size,
                               m->n_i_faces,
                               (const cs_lnum_2_t *)(m->i_face_cells),
                               da,
                               xa);

    _sles_setup[setup_id] = sc;
    _matrix_setup[setup_id][0] = a;
    _matrix_setup[setup_id][1] = NULL;
    _matrix_setup[setup_id][2] = NULL;

  }
  else
    a = _matrix_setup[setup_id][0];

  /* If system uses specific halo (i.e. when matrix contains more than
     face->cell nonzeroes), allocate specific buffers. */

  cs_real_t *_vx = vx, *_rhs = NULL;
  const cs_real_t *rhs_p = rhs;

  const cs_halo_t *halo = cs_matrix_get_halo(a);
  if (halo != NULL && halo != m->halo) {

    size_t stride = 1;
    if (diag_block_size != NULL)
      stride = diag_block_size[1];
    cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
    cs_lnum_t n_cols_ext = cs_matrix_get_n_columns(a);
    assert(n_rows == m->n_cells);
    cs_lnum_t _n_rows = n_rows*stride;
    BFT_MALLOC(_rhs, n_cols_ext*stride, cs_real_t);
    BFT_MALLOC(_vx, n_cols_ext*stride, cs_real_t);
#   pragma omp parallel for  if(_n_rows > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < _n_rows; i++) {
      _rhs[i] = rhs[i];
      _vx[i] = vx[i];
    }
    cs_matrix_pre_vector_multiply_sync(rotation_mode, a, _rhs);
    rhs_p = _rhs;
  }

  /* Solve system */

  cvg = cs_sles_solve(sc,
                      a,
                      rotation_mode,
                      precision,
                      r_norm,
                      n_iter,
                      residue,
                      rhs_p,
                      _vx,
                      0,
                      NULL);

  BFT_FREE(_rhs);
  if (_vx != vx) {
    size_t stride = 1;
    if (diag_block_size != NULL)
      stride = diag_block_size[1];
    cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
    cs_lnum_t _n_rows = n_rows*stride;
#   pragma omp parallel for  if(_n_rows > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < _n_rows; i++)
      vx[i] = _vx[i];
    BFT_FREE(_vx);
  }

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free sparse linear equation solver setup using native matrix arrays.
 *
 * \param[in]  f_id  associated field id, or < 0
 * \param[in]  name  associated name if f_id < 0, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_free_native(int          f_id,
                    const char  *name)
{
  cs_sles_t *sc = cs_sles_find(f_id, name);

  int setup_id = 0;
  while (setup_id < _n_setups) {
    if (_sles_setup[setup_id] == sc)
      break;
    else
      setup_id++;
  }

  if (setup_id < _n_setups) {

    cs_sles_free(sc);
    if (_matrix_setup[setup_id][0] != NULL)
      cs_matrix_release_coefficients(_matrix_setup[setup_id][0]);
    for (int i = 1; i < 3; i++) { /* Remove "copied" matrixes */
      if (_matrix_setup[setup_id][i] != NULL)
        cs_matrix_destroy(&(_matrix_setup[setup_id][i]));
    }

    _n_setups -= 1;

    if (setup_id < _n_setups) {
      for (int i = 0; i < 3; i++) {
        _matrix_setup[setup_id][i] = _matrix_setup[_n_setups][i];
      _sles_setup[setup_id] = _sles_setup[_n_setups];
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Error handler attempting fallback to alternative solution procedure
 *        for sparse linear equation solver.
 *
 * In case of divergence with an iterative solver, this error handler
 * switches to a default preconditioner, then resets the solution vector.
 *
 * The default error for the solver type handler is then  set, in case
 * the solution fails again.
 *
 * Note that this error handler may rebuild solver contexts, so should not
 * be used in conjunction with shared contexts (such as multigrid
 * ascent/descent contexts), but only for "outer" solvers.
 *
 * \param[in, out]  sles           pointer to solver object
 * \param[in]       state          convergence status
 * \param[in]       a              matrix
 * \param[in]       rotation_mode  halo update option for rotational periodicity
 * \param[in]       rhs            right hand side
 * \param[in, out]  vx             system solution
 *
 * \return  true if fallback solution is possible, false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_sles_default_error(cs_sles_t                    *sles,
                      cs_sles_convergence_state_t   state,
                      const cs_matrix_t            *a,
                      cs_halo_rotation_t            rotation_mode,
                      const cs_real_t               rhs[],
                      cs_real_t                     vx[])
{
  CS_UNUSED(rotation_mode);
  CS_UNUSED(rhs);

  bool alternative = false;

  if (state == CS_SLES_BREAKDOWN)
    return alternative;

  /* Case for iterative solver:
     if multigrid is not adapted to the current system,
     revert to Jacobi preconditioner
     (tested for if CS_SLES_DIVERGED or CS_SLES_MAX_ITERATION) */

  if (strcmp(cs_sles_get_type(sles), "cs_sles_it_t") == 0) {

    cs_sles_it_t  *c_old = cs_sles_get_context(sles);

    cs_sles_pc_t  *pc = cs_sles_it_get_pc(c_old);

    if (pc != NULL) {
      const char *pc_type = cs_sles_pc_get_type(pc);
      if (strcmp(pc_type, "multigrid") == 0)
        alternative = true;
    }

    if (alternative) {

      const cs_sles_it_type_t sles_it_type = cs_sles_it_get_type(c_old);

      const int f_id = cs_sles_get_f_id(sles);
      const char *name = cs_sles_get_name(sles);

      /* Switch to alternative solver if possible */

      bft_printf(_("\n\n"
                   "%s [%s]: divergence\n"
                   "  fallback from %s to Jacobi (diagonal) preconditioning\n"
                   "  for re-try and subsequent solves.\n"),
                 _(cs_sles_it_type_name[sles_it_type]), name,
                 cs_sles_pc_get_type_name(pc));

      cs_sles_free(sles);

      cs_sles_it_t  *c_new = cs_sles_it_define(f_id,
                                               name,
                                               sles_it_type,
                                               0, /* poly_degree */
                                               0);

      cs_sles_it_transfer_parameters(c_old, c_new);

    }

  } /* End of "cs_sles_it_t" case */

  else if (strcmp(cs_sles_get_type(sles), "cs_multigrid_t") == 0) {

    cs_sles_it_t  *c_old = cs_sles_get_context(sles);

    alternative = true;

    if (alternative) {

      const cs_sles_it_type_t sles_it_type = cs_sles_it_get_type(c_old);

      const int f_id = cs_sles_get_f_id(sles);
      const char *name = cs_sles_get_name(sles);

      /* Switch to alternative solver if possible */

      bft_printf(_("\n\n"
                   "%s [%s]: divergence\n"
                   "  fallback from multigrid to %s-preconditionned CG solver\n"
                   "  for re-try and subsequent solves.\n"),
                   "Multigrid", name, "Jacobi");

      cs_sles_free(sles);

      cs_sles_it_t  *c_new = cs_sles_it_define(f_id,
                                               name,
                                               sles_it_type,
                                               0, /* poly_degree */
                                               0);

      cs_sles_it_transfer_parameters(c_old, c_new);

    }

  } /* End of "cs_multigrid_t" case */

  /* Reset solution if new solve is expected */

  if (alternative) {
    const int *db_size = cs_matrix_get_diag_block_size(a);
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * db_size[1];
    for (cs_lnum_t i = 0; i < n_cols; i++)
      vx[i] = 0;
  }

  return alternative;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
