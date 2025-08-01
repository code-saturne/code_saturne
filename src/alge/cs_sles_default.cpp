/*============================================================================
 * Sparse Linear Equation Solvers
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

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_field.h"
#include "base/cs_field_default.h"
#include "alge/cs_grid.h"
#include "base/cs_halo.h"
#include "base/cs_internal_coupling.h"
#include "base/cs_log.h"
#include "base/cs_mem.h"
#include "mesh/cs_mesh.h"
#include "mesh/cs_mesh_adjacencies.h"
#include "mesh/cs_mesh_quantities.h"
#include "alge/cs_matrix.h"
#include "alge/cs_matrix_default.h"
#include "alge/cs_matrix_util.h"
#include "alge/cs_multigrid.h"
#include "base/cs_parameters.h"
#include "alge/cs_sles.h"
#include "alge/cs_sles_it.h"
#include "alge/cs_sles_pc.h"
#include "base/cs_timer.h"

#if defined(HAVE_HYPRE)
#include "alge/cs_sles_hypre.h"
#endif

#if defined(HAVE_PETSC)
#include "alge/cs_sles_petsc.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "alge/cs_sles_default.h"

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
static cs_matrix_t  *_matrix_setup[CS_SLES_DEFAULT_N_SETUPS][2];

static const int _poly_degree_default = 0;
static const int _n_max_iter_default = 10000;
static const int _n_max_iter_default_jacobi = 100;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Default definition of a sparse linear equation solver

  \param[in]  f_id         associated field id, or < 0
  \param[in]  name         associated name if f_id < 0, or nullptr
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

  if (name != nullptr) {

    if (!strcmp(name, "potential")) {   /* _cs_mass_flux_prediction */
      /* Copy from pressure if possible */
      cs_field_t *cvar_p = (cs_field_by_name_try("pressure"));
      cs_sles_t *src = nullptr;
      if (cvar_p != nullptr) {
        if (cvar_p->type & CS_FIELD_VARIABLE)
          src = cs_sles_find_or_add(cvar_p->id, nullptr);
      }
      if (src != nullptr) {
        cs_sles_t *dest = cs_sles_find_or_add(-1, name);
        if (cs_sles_copy(dest, src) == 0) /* Copy OK, we are done */
          return;
      }
      /* If copying from pressure failed, default to multigrid */
      sles_it_type = CS_SLES_FCG;
      multigrid = 1;
    }
    else if (!strcmp(name, "Prhydro")) { /* _hydrostatic_pressure_prediction */
      sles_it_type = CS_SLES_FCG;
    }
    else if (!strcmp(name, "Pr compress")) { /* cs_pressure_correction.cpp */
      sles_it_type = CS_SLES_P_SYM_GAUSS_SEIDEL;
    }
    else if (!strcmp(name, "PoissonL")) { /* _lageqp */
      sles_it_type = CS_SLES_FCG;
      n_max_iter = 1000;
    }
    else if (!strcmp(name, "radiation_p1")) { /* cs_rad_transfer_pun */
      sles_it_type = CS_SLES_FCG;
      multigrid = 1;
    }
    /* cs_bad_cells_regularisation.c */
    else if (!strcmp(name, "potential_regularisation_scalar")) {
      sles_it_type = CS_SLES_FCG;
    }
    else if (!strcmp(name, "potential_regularisation_vector")) {
      sles_it_type = CS_SLES_FCG;
    }
    else if (!strcmp(name, "potential_regularisation_sym_tensor")) {
      sles_it_type = CS_SLES_FCG;
    }
    else if (!strcmp(name, "ITM_diffusion_equation")) { /* cs_vof.cpp */
      sles_it_type = CS_SLES_FCG;
      multigrid = 1;
    }
  }
  else if (f_id > -1) {
    const cs_field_t *f = cs_field_by_id(f_id);
    if (!strcmp(f->name, "hydraulic_head")) {
      sles_it_type = CS_SLES_FCG;
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
      sles_it_type = CS_SLES_FCG;
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

  if (cs_get_device_id() > -1) {
    if (   sles_it_type == CS_SLES_P_GAUSS_SEIDEL
        || sles_it_type == CS_SLES_P_SYM_GAUSS_SEIDEL)
      sles_it_type = CS_SLES_JACOBI;
    else if (sles_it_type == CS_SLES_BICGSTAB)
      sles_it_type = CS_SLES_GCR;
  }

  if (multigrid == 1) {

    /* Multigrid used as preconditioner if possible, as solver otherwise */

    if (   (matrix_type == CS_MATRIX_MSR)
        || (matrix_type >= CS_MATRIX_N_TYPES)) {
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

  else {
    cs_sles_it_t *context = cs_sles_it_define(f_id,
                                              name,
                                              sles_it_type,
                                              _poly_degree_default,
                                              n_max_iter);

    if (   sles_it_type == CS_SLES_JACOBI
        || sles_it_type == CS_SLES_P_GAUSS_SEIDEL
        || sles_it_type == CS_SLES_P_SYM_GAUSS_SEIDEL) {
      cs_sles_it_set_fallback_threshold(context,
                                        CS_SLES_ITERATING,
                                        n_max_iter);
      int n_fallback_iter = _n_max_iter_default_jacobi;
      if (sles_it_type == CS_SLES_JACOBI)
        n_fallback_iter *= 2;
      cs_sles_it_set_n_max_iter(context, n_fallback_iter);
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Setup sparse linear equation solver using native matrix arrays.
 *
 * \param[in]       f_id          associated field id, or < 0
 * \param[in]       name          associated name if f_id < 0, or nullptr
 * \param[in]       setup_id      associated setup id
 * \param[in]       sc            associated solver context
 * \param[in]       symmetric     indicates if matrix coefficients
 *                                are symmetric
 * \param[in]       db_size       block sizes for diagonal
 * \param[in]       eb_size       block sizes for extra diagonal
 * \param[in]       da            diagonal values (nullptr if zero)
 * \param[in]       xa            extradiagonal values (nullptr if zero)
 */
/*----------------------------------------------------------------------------*/

static void
_sles_setup_matrix_native(int                  f_id,
                          const char          *name,
                          int                  setup_id,
                          cs_sles_t           *sc,
                          bool                 symmetric,
                          cs_lnum_t            db_size,
                          cs_lnum_t            eb_size,
                          const cs_real_t     *da,
                          const cs_real_t     *xa)
{
  cs_matrix_t *a = nullptr;

  const cs_mesh_t *m = cs_glob_mesh;

  bool need_msr = false;
  bool need_csr = false;
  bool need_external = false;
  char external_type[32] = "";

  /* If context has not been defined yet, temporarily set
     matrix coefficients (using native matrix, which has lowest
     overhead as coefficients are provided in that form)
     to define the required context.

     The matrix type might be modified later based on solver
     constraints. */

  if (cs_sles_get_context(sc) == nullptr) {
    a = cs_matrix_msr();

    cs_matrix_set_coefficients(a,
                               symmetric,
                               db_size,
                               eb_size,
                               m->n_i_faces,
                               m->i_face_cells,
                               da,
                               xa);

    cs_sles_define_t  *sles_default_func = cs_sles_get_default_define();
    sles_default_func(f_id, name, a);
    cs_matrix_release_coefficients(a);
  }

  assert(cs_sles_get_context(sc) != nullptr);

  cs_sles_pc_t  *pc = nullptr;
  cs_multigrid_t *mg = nullptr;

  if (strcmp(cs_sles_get_type(sc), "cs_sles_it_t") == 0) {
    cs_sles_it_t     *c = static_cast<cs_sles_it_t *>(cs_sles_get_context(sc));
    cs_sles_it_type_t s_type = cs_sles_it_get_type(c);
    if (   s_type >= CS_SLES_P_GAUSS_SEIDEL
        && s_type <= CS_SLES_TS_B_GAUSS_SEIDEL) {
      need_msr = true;
      /* Gauss-Seidel not yet implemented with full blocks */
      if (eb_size > 1)
        need_msr = false;
    }
    pc = cs_sles_it_get_pc(c);
    if (pc != nullptr) {
      if (strcmp(cs_sles_pc_get_type(pc), "multigrid") == 0)
        mg = static_cast<cs_multigrid_t *>(cs_sles_pc_get_context(pc));
    }
  }
  else if (strcmp(cs_sles_get_type(sc), "cs_multigrid_t") == 0)
    mg = static_cast<cs_multigrid_t *>(cs_sles_get_context(sc));

#if defined(HAVE_HYPRE)
  else if (strcmp(cs_sles_get_type(sc), "cs_sles_hypre_t") == 0) {
    cs_sles_hypre_t *c =
      static_cast<cs_sles_hypre_t *>(cs_sles_get_context(sc));
    int use_device = cs_sles_hypre_get_host_device(c);
    if (use_device)
      strncpy(external_type, "HYPRE_ParCSR, device", 31);
    else
      strncpy(external_type, "HYPRE_ParCSR", 31);
    need_external = true;
  }
#endif

#if defined(HAVE_PETSC)
  else if (   strcmp(cs_sles_get_type(sc), "cs_sles_petsc_t") == 0
           && m->have_rotation_perio == 0) {
    void *c = cs_sles_get_context(sc);
    const char *mat_type = cs_sles_petsc_get_mat_type(c);
    if (mat_type == nullptr)
      strncpy(external_type, "PETSc", 31);
    else {
      snprintf(external_type, 31, "PETSc, %s", mat_type);
      external_type[31] = '\0';
    }
    need_external = true;
  }
#endif

  if (mg != nullptr) {
    if (cs_multigrid_need_msr(mg))
      need_msr = true;
  }

  /* For now, simply force MSR when running on accelerated node */
  if (need_external == false && cs_get_device_id() > -1)
    need_msr = true;

#if defined(HAVE_CUDSS)
  if (strcmp(cs_sles_get_type(sc), "cs_sles_cudss_t") == 0) {
    need_msr = false;
    need_csr = true;
  }
#endif

  if (need_msr)
    a = cs_matrix_msr();
  else if (need_csr)
    a = cs_matrix_csr();
  else if (need_external) {
    a = cs_matrix_external(external_type,
                           symmetric,
                           db_size,
                           eb_size);
  }
  else
    a = cs_matrix_default(symmetric,
                          db_size,
                          eb_size);

  cs_matrix_set_coefficients(a,
                             symmetric,
                             db_size,
                             eb_size,
                             m->n_i_faces,
                             m->i_face_cells,
                             da,
                             xa);

  if (mg != nullptr) {
    const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
    const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

    if (ma->cell_i_faces == nullptr)
      cs_mesh_adjacencies_update_cell_i_faces();

    cs_matrix_set_mesh_association(a,
                                   ma->cell_cells_idx,
                                   ma->cell_i_faces,
                                   ma->cell_i_faces_sgn,
                                   mq->cell_cen,
                                   mq->cell_vol,
                                   mq->i_face_u_normal,
                                   mq->i_face_surf);
  }

  cs_matrix_default_set_tuned(a);

  _matrix_setup[setup_id][0] = a;
  _matrix_setup[setup_id][1] = nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call sparse linear equation solver setup for systems requiring
 *        a general matrix assember..
 *
 * \param[in]  f_id                   associated field id, or < 0
 * \param[in]  setup_id               associated setup id
 * \param[in]  sc                     associated solver context
 * \param[in]  symmetric              indicates if matrix coefficients
 *                                    are symmetric
 * \param[in]  db_size                block sizes for diagonal
 * \param[in]  eb_size                block sizes for extra diagonal
 * \param[in]  da                     diagonal values (nullptr if zero)
 * \param[in]  xa                     extradiagonal values (nullptr if zero)
 */
/*----------------------------------------------------------------------------*/

static void
_sles_setup_matrix_by_assembler(int               f_id,
                                int               setup_id,
                                cs_sles_t        *sc,
                                bool              symmetric,
                                const cs_lnum_t   db_size,
                                const cs_lnum_t   eb_size,
                                const cs_real_t  *da,
                                const cs_real_t  *xa)
{
  cs_matrix_t *a = nullptr;

  /* Check field and associated internal coupling if present */

  const cs_field_t *f = nullptr;
  if (f_id > -1)
    f = cs_field_by_id(f_id);

  cs_matrix_type_t m_type = CS_MATRIX_MSR;
  if (strcmp(cs_sles_get_type(sc), "cs_sles_amgx_t") == 0)
    m_type = CS_MATRIX_CSR;

  /* Set coefficients */

  a = cs_matrix_set_coefficients_by_assembler(f,
                                              m_type,
                                              symmetric,
                                              db_size,
                                              eb_size,
                                              da,
                                              xa);

  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  cs_matrix_set_mesh_association(a,
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 mq->cell_cen,
                                 mq->cell_vol,
                                 nullptr,
                                 nullptr);

  cs_matrix_default_set_tuned(a);

  _matrix_setup[setup_id][0] = a;
  _matrix_setup[setup_id][1] = a; /* so it is freed later */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Default definition of a sparse linear equation solver

  \param[in]  f_id  associated field id, or < 0
  \param[in]  name  associated name if f_id < 0, or nullptr
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

  /* Define for all variable fields */

  const int n_fields = cs_field_n_fields();

  /* Define solver for all variable fields if not already done,
     based on convection/diffusion */

  for (int f_id = 0; f_id < n_fields; f_id++) {
    const cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {

      if (f->type & CS_FIELD_CDO) /* Skipped this step for CDO equations */
        continue;                 /* This is done elsewhere. */

      void *context = nullptr;
      cs_sles_t *sc = cs_sles_find(f->id, nullptr);
      if (sc != nullptr)
        context = cs_sles_get_context(sc);

      if (context == nullptr) {
        /* Get the calculation option from the field */
        const cs_equation_param_t *eqp = cs_field_get_equation_param_const(f);
        if (eqp != nullptr) {
          bool symmetric = (eqp->iconv > 0) ? false : true;
          _sles_default_native(f_id, nullptr, CS_MATRIX_N_TYPES, symmetric);
        }
      }

    }

  }

  /* Logging */

  cs_log_printf(CS_LOG_SETUP, "\n");
  cs_log_separator(CS_LOG_SETUP);

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
 * \param[in]  name  associated name if f_id < 0, or nullptr
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

  if (k_log < 0)
    k_log = cs_field_key_id("log");

  if (f_id > -1) {
    const cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      retval = cs_field_get_equation_param_const(f)->verbosity;
    }
    else
      retval = cs_field_get_key_int(f, k_log);
  }

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return pointer to matrix structure matching equation solve.
 *
 * Some matrix properties (such as assembly options and geometric
 * association) are set immediately, but coefficients are not
 * assigned at this stage.
 *
 * \param[in]  f_id       associated field id, or < 0
 * \param[in]  name       associated name if f_id < 0, or nullptr
 * \param[in]  db_size    block sizes for diagonal
 * \param[in]  eb_size    block sizes for extra diagonal
 * \param[in]  symmetric  indicate if matrix is symmetric
 *
 * \return  pointer to matrix structure
 */
/*----------------------------------------------------------------------------*/

cs_matrix_t *
cs_sles_default_get_matrix(int          f_id,
                           const char  *name,
                           cs_lnum_t    db_size,
                           cs_lnum_t    eb_size,
                           bool         symmetric)
{
  cs_matrix_t *a = nullptr;

#if defined(HAVE_PETSC)
  const cs_mesh_t *m = cs_glob_mesh;
#endif
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Check if this system has already been setup */

  const cs_field_t *f = nullptr;
  cs_sles_t *sc = cs_sles_find_or_add(f_id, name);

  /* If context has not been defined yet, force default definition.

     The matrix type might be modified later based on solver
     constraints. */

  if (cs_sles_get_context(sc) == nullptr) {
    _sles_default_native(f_id, name, CS_MATRIX_MSR, symmetric);
  }

  bool need_matrix_assembler = false;
  cs_matrix_type_t mat_type = CS_MATRIX_MSR;

  if (f_id > -1) {
    f = cs_field_by_id(f_id);
    int coupling_id
      = cs_field_get_key_int(f, cs_field_key_id("coupling_entity"));
    if (coupling_id > -1)
      need_matrix_assembler = true;
  }

  const char *sles_type = cs_sles_get_type(sc);
  if (strcmp(sles_type, "cs_sles_amgx_t") == 0) {
    need_matrix_assembler = true;
    mat_type = CS_MATRIX_CSR;
  }

  /* For CSR or MSR internal matrix with assembler, mesh association
     will not contain face information as some coefficients do not
     match faces */

  if (need_matrix_assembler) {
    a = cs_matrix_by_assembler(f, mat_type);

    cs_matrix_set_mesh_association(a,
                                   nullptr,
                                   nullptr,
                                   nullptr,
                                   mq->cell_cen,
                                   mq->cell_vol,
                                   nullptr,
                                   nullptr);

    cs_matrix_set_need_xa(a, false);

    return a;
  }

  /* For general case, the format will depend on the type of solver
     used */

  bool need_msr = false;
  bool need_csr = false;
  bool need_external = false;
  char external_type[32] = "";

  cs_sles_pc_t  *pc = nullptr;
  cs_multigrid_t *mg = nullptr;

  if (strcmp(cs_sles_get_type(sc), "cs_sles_it_t") == 0) {
    cs_sles_it_t     *c = static_cast<cs_sles_it_t *>(cs_sles_get_context(sc));
    cs_sles_it_type_t s_type = cs_sles_it_get_type(c);
    if (   s_type >= CS_SLES_P_GAUSS_SEIDEL
        && s_type <= CS_SLES_TS_B_GAUSS_SEIDEL) {
      need_msr = true;
      /* Gauss-Seidel not yet implemented with full blocks */
      if (eb_size > 1)
        need_msr = false;
    }
    pc = cs_sles_it_get_pc(c);
    if (pc != nullptr) {
      if (strcmp(cs_sles_pc_get_type(pc), "multigrid") == 0)
        mg = static_cast<cs_multigrid_t *>(cs_sles_pc_get_context(pc));
    }
  }
  else if (strcmp(cs_sles_get_type(sc), "cs_multigrid_t") == 0)
    mg = static_cast<cs_multigrid_t *>(cs_sles_get_context(sc));

#if defined(HAVE_HYPRE)
  else if (strcmp(cs_sles_get_type(sc), "cs_sles_hypre_t") == 0) {
    cs_sles_hypre_t *c =
      static_cast<cs_sles_hypre_t *>(cs_sles_get_context(sc));
    int use_device = cs_sles_hypre_get_host_device(c);
    if (use_device)
      strncpy(external_type, "HYPRE_ParCSR, device", 31);
    else
      strncpy(external_type, "HYPRE_ParCSR", 31);
    need_external = true;
  }
#endif

#if defined(HAVE_PETSC)
  else if (   strcmp(cs_sles_get_type(sc), "cs_sles_petsc_t") == 0
           && m->have_rotation_perio == 0) {
    void *c = cs_sles_get_context(sc);
    const char *mat_type_s = cs_sles_petsc_get_mat_type(c);
    if (mat_type_s == nullptr)
      strncpy(external_type, "PETSc", 31);
    else {
      snprintf(external_type, 31, "PETSc, %s", mat_type_s);
      external_type[31] = '\0';
    }
    need_external = true;
  }
#endif

  if (mg != nullptr) {
    if (cs_multigrid_need_msr(mg))
      need_msr = true;
  }

  /* For now, simply force MSR when running on accelerated node */
  if (need_external == false && cs_get_device_id() > -1)
    need_msr = true;

  /* MSR is the current default; DIST could be a future option
     if its performance can be improved, while NATIVE is deprecated. */
  need_msr = true;

#if defined(HAVE_CUDSS)
  if (strcmp(cs_sles_get_type(sc), "cs_sles_cudss_t") == 0) {
    need_csr = true;
    need_msr = false;
  }
#endif

  if (need_external)
    a = cs_matrix_external(external_type,
                           symmetric,
                           db_size,
                           eb_size);

  else if (need_msr)
    a = cs_matrix_msr();
  else if (need_csr)
    a = cs_matrix_csr();

  else
    a = cs_matrix_default(symmetric,
                          db_size,
                          eb_size);

  bool need_xa = false;

  if (mg != nullptr) {
    if (symmetric == false) need_xa = true;

    const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;

    if (ma->cell_i_faces == nullptr)
      cs_mesh_adjacencies_update_cell_i_faces();

    cs_matrix_set_mesh_association(a,
                                   ma->cell_cells_idx,
                                   ma->cell_i_faces,
                                   ma->cell_i_faces_sgn,
                                   mq->cell_cen,
                                   mq->cell_vol,
                                   mq->i_face_u_normal,
                                   mq->i_face_surf);

  }

  cs_matrix_set_need_xa(a, need_xa);

  return a;
}

/*----------------------------------------------------------------------------
 * Release of destroy matrix depending on whether is is cached or not.
 *
 * Matrices built by assembler are destroyed.
 *
 * parameters:
 *   matrix <-> pointer to matrix structure pointer
 *----------------------------------------------------------------------------*/

void
cs_sles_default_release_matrix(cs_matrix_t  **m)
{
  cs_matrix_release(m);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call sparse linear equation solver setup for convection-diffusion
 *        systems
 *
 * \param[in]       f_id                   associated field id, or < 0
 * \param[in]       name                   associated name if f_id < 0, or nullptr
 * \param[in]       diag_block_size        block sizes for diagonal, or nullptr
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal,
 *                                         or nullptr
 * \param[in]       da                     diagonal values (nullptr if zero)
 * \param[in]       xa                     extradiagonal values (nullptr if zero)
 * \param[in]       conv_diff              convection-diffusion mode
 */
/*----------------------------------------------------------------------------*/

void
cs_sles_setup_native_conv_diff(int                  f_id,
                               const char          *name,
                               const cs_lnum_t      diag_block_size,
                               const cs_lnum_t      extra_diag_block_size,
                               const cs_real_t     *da,
                               const cs_real_t     *xa,
                               bool                 conv_diff)
{
  cs_matrix_t *a = nullptr;

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

    if (a == nullptr) {

      a = cs_matrix_msr();

      cs_matrix_set_coefficients(a,
                                 false,
                                 diag_block_size,
                                 extra_diag_block_size,
                                 m->n_i_faces,
                                 m->i_face_cells,
                                 da,
                                 xa);

      const cs_mesh_adjacencies_t *ma = cs_glob_mesh_adjacencies;
      const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

      if (ma->cell_i_faces == nullptr)
        cs_mesh_adjacencies_update_cell_i_faces();

      cs_matrix_set_mesh_association(a,
                                     ma->cell_cells_idx,
                                     ma->cell_i_faces,
                                     ma->cell_i_faces_sgn,
                                     mq->cell_cen,
                                     mq->cell_vol,
                                     mq->i_face_u_normal,
                                     mq->i_face_surf);

    }

    _sles_setup[setup_id] = sc;
    _matrix_setup[setup_id][0] = a;

  }
  else {
    a = _matrix_setup[setup_id][0];
  }

  cs_matrix_default_set_tuned(a);

  /* Setup system */

  if (strcmp(cs_sles_get_type(sc), "cs_multigrid_t") != 0)
    bft_error
      (__FILE__, __LINE__, 0,
       "%s requires a cs_multigrid_t solver type", __func__);

  int verbosity = cs_sles_get_verbosity(sc);

  cs_multigrid_t *mg = static_cast<cs_multigrid_t *>(cs_sles_get_context(sc));
  cs_multigrid_setup_conv_diff(mg, name, a, conv_diff, verbosity);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call sparse linear equation solver for general colocated
 *        cell-centered finite volume scheme.
 *
 * The initial solution is assumed to be 0 (and does not need to
 * be initialized before calling this function).
 *
 * \param[in]       sc                     solver context
 * \param[in]       a                      matrix
 * \param[in]       precision              solver precision
 * \param[in]       r_norm                 residual normalization
 * \param[out]      n_iter                 number of "equivalent" iterations
 * \param[out]      residual               residual
 * \param[in]       rhs                    right hand side
 * \param[out]      vx                     system solution
 *
 * \return  convergence state
 */
/*----------------------------------------------------------------------------*/

cs_sles_convergence_state_t
cs_sles_solve_ccc_fv(cs_sles_t           *sc,
                     cs_matrix_t         *a,
                     double               precision,
                     double               r_norm,
                     int                 *n_iter,
                     double              *residual,
                     const cs_real_t     *rhs,
                     cs_real_t           *vx)
{
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;

  const cs_mesh_t *m = cs_glob_mesh;

  /* If system uses specific halo (i.e. when matrix contains more than
     face->cell nonzeroes), allocate specific buffers and synchronize
     right hand side. */

  cs_real_t *_vx = vx, *_rhs = nullptr;
  const cs_real_t *rhs_p = rhs;

  cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
  cs_lnum_t n_cols_ext = cs_matrix_get_n_columns(a);
  cs_assert(n_rows == m->n_cells);

  if (n_cols_ext > m->n_cells_with_ghosts) {
    cs_alloc_mode_t amode = cs_matrix_get_alloc_mode(a);

    cs_dispatch_context ctx;
    ctx.set_use_gpu(amode >= CS_ALLOC_HOST_DEVICE_SHARED);

    cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
    assert(n_rows == m->n_cells);
    cs_lnum_t _n_rows = n_rows*db_size;
    CS_MALLOC_HD(_vx, n_cols_ext*db_size, cs_real_t, amode);

    if (cs_matrix_get_type(a) == CS_MATRIX_NATIVE) {
      CS_MALLOC_HD(_rhs, n_cols_ext*db_size, cs_real_t, amode);
      rhs_p = _rhs;

      cs_lnum_t _n_ext = (n_cols_ext - n_rows)*db_size;
      ctx.parallel_for(_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
        _rhs[i] = rhs[i];
      });
      ctx.parallel_for(_n_ext, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
        _rhs[_n_rows + i] = 0;
      });
    };

    ctx.wait();
  }

  /* Solve system */

  cvg = cs_sles_solve(sc,
                      a,
                      precision,
                      r_norm,
                      true,
                      n_iter,
                      residual,
                      rhs_p,
                      _vx,
                      0,
                      nullptr);

  CS_FREE(_rhs);
  if (_vx != vx) {

    cs_alloc_mode_t amode = cs_matrix_get_alloc_mode(a);
    cs_dispatch_context ctx;
    ctx.set_use_gpu(amode >= CS_ALLOC_HOST_DEVICE_SHARED);

    cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
    cs_lnum_t _n_rows = n_rows*db_size;
    ctx.parallel_for(_n_rows, [=] CS_F_HOST_DEVICE (cs_lnum_t i) {
      vx[i] = _vx[i];
    });
    CS_FREE(_vx);

  }

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Call sparse linear equation solver using native matrix arrays.
 *
 * \param[in]       f_id                   associated field id, or < 0
 * \param[in]       name                   associated name if f_id < 0, or nullptr
 * \param[in]       symmetric              indicates if matrix coefficients
 *                                         are symmetric
 * \param[in]       diag_block_size        block sizes for diagonal
 * \param[in]       extra_diag_block_size  block sizes for extra diagonal
 * \param[in]       da                     diagonal values (nullptr if zero)
 * \param[in]       xa                     extradiagonal values (nullptr if zero)
 * \param[in]       precision              solver precision
 * \param[in]       r_norm                 residual normalization
 * \param[out]      n_iter                 number of "equivalent" iterations
 * \param[out]      residual               residual
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
                     cs_lnum_t            diag_block_size,
                     cs_lnum_t            extra_diag_block_size,
                     const cs_real_t     *da,
                     const cs_real_t     *xa,
                     double               precision,
                     double               r_norm,
                     int                 *n_iter,
                     double              *residual,
                     const cs_real_t     *rhs,
                     cs_real_t           *vx)
{
  cs_sles_convergence_state_t cvg = CS_SLES_ITERATING;
  cs_matrix_t *a = nullptr;

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

    /* Check if we need to used a matrix assembler */

    bool need_matrix_assembler = false;

    if (f_id > -1) {
      const cs_field_t *f = cs_field_by_id(f_id);
      int coupling_id
        = cs_field_get_key_int(f, cs_field_key_id("coupling_entity"));
      if (coupling_id > -1)
        need_matrix_assembler = true;
    }
    const char *sles_type = cs_sles_get_type(sc);
    if (strcmp(sles_type, "cs_sles_amgx_t") == 0) {
      need_matrix_assembler = true;
    }

    if (need_matrix_assembler)
      _sles_setup_matrix_by_assembler(f_id,
                                      setup_id,
                                      sc,
                                      symmetric,
                                      diag_block_size,
                                      extra_diag_block_size,
                                      da,
                                      xa);

    else
      _sles_setup_matrix_native(f_id,
                                name,
                                setup_id,
                                sc,
                                symmetric,
                                diag_block_size,
                                extra_diag_block_size,
                                da,
                                xa);

    _sles_setup[setup_id] = sc;
  }

  a = _matrix_setup[setup_id][0];

  /* If system uses specific halo (i.e. when matrix contains more than
     face->cell nonzeroes), allocate specific buffers and synchronize
     right hand side. */

  cs_real_t *_vx = vx, *_rhs = nullptr;
  const cs_real_t *rhs_p = rhs;

  const cs_halo_t *halo = cs_matrix_get_halo(a);
  if (halo != nullptr && halo != m->halo) {

    size_t stride = diag_block_size;
    cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
    cs_lnum_t n_cols_ext = cs_matrix_get_n_columns(a);
    assert(n_rows == m->n_cells);
    cs_lnum_t _n_rows = n_rows*stride;
    CS_MALLOC(_rhs, n_cols_ext*stride, cs_real_t);
    CS_MALLOC(_vx, n_cols_ext*stride, cs_real_t);
#   pragma omp parallel for  if(_n_rows > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < _n_rows; i++) {
      _rhs[i] = rhs[i];
      _vx[i] = vx[i];
    }
    cs_matrix_pre_vector_multiply_sync(a, _rhs);
    rhs_p = _rhs;

  }

  /* Solve system */

  cvg = cs_sles_solve(sc,
                      a,
                      precision,
                      r_norm,
                      false,
                      n_iter,
                      residual,
                      rhs_p,
                      _vx,
                      0,
                      nullptr);

  CS_FREE(_rhs);
  if (_vx != vx) {
    size_t stride = diag_block_size;
    cs_lnum_t n_rows = cs_matrix_get_n_rows(a);
    cs_lnum_t _n_rows = n_rows*stride;
#   pragma omp parallel for  if(_n_rows > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < _n_rows; i++)
      vx[i] = _vx[i];
    CS_FREE(_vx);
  }

  return cvg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free sparse linear equation solver setup using native matrix arrays.
 *
 * \param[in]  f_id  associated field id, or < 0
 * \param[in]  name  associated name if f_id < 0, or nullptr
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
    if (_matrix_setup[setup_id][0] != nullptr)
      cs_matrix_release_coefficients(_matrix_setup[setup_id][0]);
    /* Remove "copied" matrices */
    if (_matrix_setup[setup_id][1] != nullptr)
      cs_matrix_destroy(&(_matrix_setup[setup_id][1]));

    _n_setups -= 1;

    if (setup_id < _n_setups) {
      for (int i = 0; i < 2; i++) {
        _matrix_setup[setup_id][i] = _matrix_setup[_n_setups][i];
      }
      _sles_setup[setup_id] = _sles_setup[_n_setups];
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
                      const cs_real_t               rhs[],
                      cs_real_t                     vx[])
{
  CS_UNUSED(rhs);

  bool alternative = false;

  if (state == CS_SLES_BREAKDOWN)
    return alternative;

  /* Case for iterative solver:
     if multigrid is not adapted to the current system,
     revert to Jacobi preconditioner
     (tested for if CS_SLES_DIVERGED or CS_SLES_MAX_ITERATION) */

  if (strcmp(cs_sles_get_type(sles), "cs_sles_it_t") == 0) {
    cs_sles_it_t *c_old =
      static_cast<cs_sles_it_t *>(cs_sles_get_context(sles));

    cs_sles_pc_t  *pc = cs_sles_it_get_pc(c_old);

    if (pc != nullptr) {
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
    cs_sles_it_t *c_old =
      static_cast<cs_sles_it_t *>(cs_sles_get_context(sles));

    alternative = true;

    if (alternative) {

      const cs_sles_it_type_t sles_it_type = cs_sles_it_get_type(c_old);

      const int f_id = cs_sles_get_f_id(sles);
      const char *name = cs_sles_get_name(sles);

      /* Switch to alternative solver if possible */

      bft_printf(_("\n\n"
                   "%s [%s]: divergence\n"
                   "  fallback from multigrid to %s-preconditioned CG solver\n"
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
    const cs_lnum_t db_size = cs_matrix_get_diag_block_size(a);
    const cs_lnum_t n_cols = cs_matrix_get_n_columns(a) * db_size;
    for (cs_lnum_t i = 0; i < n_cols; i++)
      vx[i] = 0;
  }

  return alternative;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
