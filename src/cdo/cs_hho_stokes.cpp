/*============================================================================
 * Build an algebraic system for scalar conv./diff. eq. with Hybrid High Order
 * space discretization
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_mem.h"

#include "base/cs_boundary_zone.h"
#include "cdo/cs_equation_priv.h"
#include "cdo/cs_hho_builder.h"
#include "base/cs_log.h"
#include "base/cs_math.h"
#include "mesh/cs_mesh_location.h"
#include "base/cs_post.h"
#include "cdo/cs_quadrature.h"
#include "cdo/cs_reco.h"
#include "cdo/cs_scheme_geometry.h"
#include "base/cs_search.h"
#include "cdo/cs_sdm.h"
#include "cdo/cs_source_term.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cdo/cs_hho_stokes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_HHO_STOKES_DBG  1
#define CS_HHO_STOKES_MODULO  4

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/* Algebraic system for HHO discretization */

struct _cs_hho_stokes_t {

  /* System size (n_faces * n_face_dofs + n_cells * n_cell_dofs) */
  cs_lnum_t                      n_dofs;
  int                            n_max_loc_dofs;
  int                            n_cell_dofs;
  int                            n_face_dofs;

  /* Solution of the algebraic system at the last computed iteration.
     cell_values is different from the values stored in the associated
     field since here it's the values of polynomial coefficients which
     are stored */
  cs_real_t                     *face_values;  /* DoF unknowns (x) + BCs */
  cs_real_t                     *cell_values;  /* DoF recomputed after the
                                                  static condensation */

  /* Storage of the source term (if one wants to apply a specific time
     discretization) */
  cs_real_t                     *source_terms;

  /* Handle the definition of the BCs */
  short int                     *bf2def_ids;

  /* Static condensation members */
  /* =========================== */

  /* Acc_inv = Inverse of a diagonal matrix (block cell-cell)
     rc_tilda = Acc_inv * rhs_c (perform for each cell) */
  cs_real_t                     *rc_tilda;

  /* Lower-Left block of the full matrix (block cell-vertices).
     Access to the values thanks to the c2f connectivity
     The transposed version of each block is stored for a more efficient
     usage */
  cs_sdm_t                      *acf_tilda;

};

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */
static cs_cell_sys_t     **cs_hho_cell_sys = nullptr;
static cs_cell_builder_t **cs_hho_cell_bld = nullptr;
static cs_hho_builder_t  **cs_hho_builders = nullptr;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;
static const cs_matrix_assembler_t  *cs_shared_ma0;
static const cs_matrix_structure_t  *cs_shared_ms0;
static const cs_matrix_assembler_t  *cs_shared_ma1;
static const cs_matrix_structure_t  *cs_shared_ms1;
static const cs_matrix_assembler_t  *cs_shared_ma2;
static const cs_matrix_structure_t  *cs_shared_ms2;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local builder structure used for building the system
 *          cellwise
 *
 * \param[in]  space_scheme   discretization scheme
 * \param[in]  connect        pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_cell_builder_t *
_cell_builder_create(cs_param_space_scheme_t     space_scheme,
                     const cs_cdo_connect_t     *connect)
{
  int  size;

  const int  n_fc = connect->n_max_fbyc;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  switch (space_scheme) {

  case CS_SPACE_SCHEME_HHO_P0:  /* TODO */
    {
    CS_MALLOC(cb->ids, n_fc + 1, int);
    memset(cb->ids, 0, (n_fc + 1) * sizeof(int));

    /* For post-processing errors = 38 */
    size = cs::max(38, n_fc * (n_fc + 1));
    CS_MALLOC(cb->values, size, double);
    memset(cb->values, 0, size * sizeof(cs_real_t));

    size = cs::max(2 * n_fc, 15);
    CS_MALLOC(cb->vectors, size, cs_real_3_t);
    memset(cb->vectors, 0, size * sizeof(cs_real_3_t));

    /* Local square dense matrices used during the construction of
       operators */
    cb->loc = cs_sdm_square_create(n_fc + 1);
    cb->aux = cs_sdm_square_create(n_fc + 1);
    }
    break;

  case CS_SPACE_SCHEME_HHO_P1:
    {
      /* Store the block size description */
      size = n_fc + 1;
      CS_MALLOC(cb->ids, size, int);
      memset(cb->ids, 0, size*sizeof(int));

      /* Store the face, cell and gradient basis function evaluations and
         the Gauss point weights
         --->  42 = 3 + 4 + 3*10 + 5
         or the factorization of the stiffness matrix used for computing
         the gradient reconstruction operator
         --->  n = 9 (10 - 1) --> facto = n*(n+1)/2 = 45
                              --> tmp buffer for facto = n = 9
                              --> 45 + 9 = 54
         or the factorization of the cell_cell block of size 4 --> 10
      */
      size = cs::max(54, (3*n_fc + 4)*2);
      CS_MALLOC(cb->values, size, double);
      memset(cb->values, 0, size*sizeof(cs_real_t));

      /* Store Gauss points and tensor.n_f products */
      size = cs::max(15, 5 + n_fc);
      CS_MALLOC(cb->vectors, size, cs_real_3_t);
      memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

      /* Local dense matrices used during the construction of operators */
      const int g_size = 9;                   /* basis (P_(k+1)) - 1 */
      for (int i = 0; i < n_fc; i++) cb->ids[i] = 3;
      cb->ids[n_fc] = 4;
      cb->loc = cs_sdm_block_create(n_fc + 1, n_fc + 1, cb->ids, cb->ids);
      cb->aux = cs_sdm_block_create(n_fc + 1, 1, cb->ids, &g_size); /* R_g^T */
    }
    break;

  case CS_SPACE_SCHEME_HHO_P2:
    {
      /* Store the block size description */
      size = n_fc + 1;
      CS_MALLOC(cb->ids, size, int);
      memset(cb->ids, 0, size*sizeof(int));

      /* Store the face, cell and gradient basis function evaluations and
         the Gauss point weights */
      /* Store the face, cell and gradient basis function evaluations and
         the Gauss point weights
         --->  91 = 6 (f) + 10 (c) + 3*20 (g) + 15 (gauss)
         or the factorization of the stiffness matrix used for computing
         the gradient reconstruction operator
         --->  n = 19 (20 - 1) --> facto = n*(n+1)/2 = 190
                               --> tmp buffer for facto = n = 19
                              --> 190 + 19 = 209
      */
      size = cs::max(209, 2*(6*n_fc + 20));
      CS_MALLOC(cb->values, size, double);
      memset(cb->values, 0, size*sizeof(cs_real_t));

      size = 15 + n_fc;  /* Store Gauss points and tensor.n_f products */
      CS_MALLOC(cb->vectors, size, cs_real_3_t);
      memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

      /* Local dense matrices used during the construction of operators */
      const int g_size = 19; /* basis (P_(k+1)) - 1 */
      for (int i = 0; i < n_fc; i++) cb->ids[i] = 6;
      cb->ids[n_fc] = 10;
      cb->loc = cs_sdm_block_create(n_fc + 1, n_fc + 1, cb->ids, cb->ids);
      cb->aux = cs_sdm_block_create(n_fc + 1, 1, cb->ids, &g_size); /* R_g^T */
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _("Invalid space scheme."));

  } // End of switch on space scheme

  return cb;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to HHO schemes
 *         Set shared pointers
 *
 * \param[in]  scheme_flag  flag to identify which kind of numerical scheme is
 *                          requested to solve the computational domain
 * \param[in]  quant        additional mesh quantities struct.
 * \param[in]  connect      pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step    pointer to a time step structure
 * \param[in]  ma0          pointer to a cs_matrix_assembler_t structure (P0)
 * \param[in]  ma1          pointer to a cs_matrix_assembler_t structure (P1)
 * \param[in]  ma2          pointer to a cs_matrix_assembler_t structure (P2)
 * \param[in]  ms0          pointer to a cs_matrix_structure_t structure (P0)
 * \param[in]  ms1          pointer to a cs_matrix_structure_t structure (P1)
 * \param[in]  ms2          pointer to a cs_matrix_structure_t structure (P2)
*/
/*----------------------------------------------------------------------------*/

void
cs_hho_stokes_initialize(cs_flag_t                      scheme_flag,
                         const cs_cdo_quantities_t     *quant,
                         const cs_cdo_connect_t        *connect,
                         const cs_time_step_t          *time_step,
                         const cs_matrix_assembler_t   *ma0,
                         const cs_matrix_assembler_t   *ma1,
                         const cs_matrix_assembler_t   *ma2,
                         const cs_matrix_structure_t   *ms0,
                         const cs_matrix_structure_t   *ms1,
                         const cs_matrix_structure_t   *ms2)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
  cs_shared_ma0 = ma0;
  cs_shared_ms0 = ms0;
  cs_shared_ma1 = ma1;
  cs_shared_ms1 = ms1;
  cs_shared_ma2 = ma2;
  cs_shared_ms2 = ms2;

  const int n_fc = connect->n_max_fbyc;

  int  order = -1, fbs = 0, cbs = 0;
  cs_param_space_scheme_t  space_scheme;

  if (scheme_flag & CS_FLAG_SCHEME_POLY2) {

    space_scheme = CS_SPACE_SCHEME_HHO_P2;
    fbs = CS_N_DOFS_FACE_2ND; // DoF by face
    cbs = CS_N_DOFS_CELL_2ND; // DoF for the cell
    order = 2;

  }
  else if (scheme_flag & CS_FLAG_SCHEME_POLY1) {

    space_scheme = CS_SPACE_SCHEME_HHO_P1;
    fbs = CS_N_DOFS_FACE_1ST; // DoF by face
    cbs = CS_N_DOFS_CELL_1ST;  // DoF for the cell
    order = 1;

  }
  else {

    space_scheme = CS_SPACE_SCHEME_HHO_P0;
    fbs = CS_N_DOFS_FACE_0TH; // DoF by face
    cbs = CS_N_DOFS_CELL_0TH; // DoF for the cell
    order = 0;

  }

  const int n_dofs = n_fc * fbs + cbs;

    /* Structure used to build the final system by a cell-wise process */
  assert(cs_glob_n_threads > 0);  /* Sanity check */

  CS_MALLOC(cs_hho_cell_bld, cs_glob_n_threads, cs_cell_builder_t *);
  CS_MALLOC(cs_hho_cell_sys, cs_glob_n_threads, cs_cell_sys_t *);

  /* Allocate builder structure specific to HHO schemes. This is an additional
     structure with respect to a cs_cell_builder_t structure */
  CS_MALLOC(cs_hho_builders, cs_glob_n_threads, cs_hho_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_hho_cell_bld[i] = nullptr;
    cs_hho_cell_sys[i] = nullptr;
    cs_hho_builders[i] = nullptr;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_builder_t  *cb = _cell_builder_create(space_scheme, connect);
    cs_hho_cell_bld[t_id] = cb;
    cs_hho_builders[t_id] = cs_hho_builder_create(order, n_fc);

    for (int i = 0; i < n_fc; i++) cb->ids[i] = fbs;
    cb->ids[n_fc] = cbs;

    cs_hho_cell_sys[t_id] = cs_cell_sys_create(n_dofs,
                                               fbs*n_fc,
                                               n_fc + 1, cb->ids);
  }
#else
  assert(cs_glob_n_threads == 1);

    cs_cell_builder_t  *cb = _cell_builder_create(space_scheme, connect);
    cs_hho_cell_bld[0] = cb;
    cs_hho_builders[0] = cs_hho_builder_create(order, n_fc);

    for (int i = 0; i < n_fc; i++) cb->ids[i] = fbs;
    cb->ids[n_fc] = cbs;

    cs_hho_cell_sys[0] = cs_cell_sys_create(n_dofs,
                                               fbs*n_fc,
                                               n_fc + 1, cb->ids);
#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys    pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb      pointer to a pointer on a cs_cell_builder_t structure
 * \param[out]  hhob    pointer to a pointer on a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_stokes_get(cs_cell_sys_t       **csys,
                  cs_cell_builder_t   **cb,
                  cs_hho_builder_t    **hhob)
{
  const int  t_id = cs_get_thread_id();

  *csys = cs_hho_cell_sys[t_id];
  *cb = cs_hho_cell_bld[t_id];
  *hhob = cs_hho_builders[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free buffers and generic structures related to HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_stokes_finalize(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();

    cs_cell_sys_free(&(cs_hho_cell_sys[t_id]));
    cs_cell_builder_free(&(cs_hho_cell_bld[t_id]));
    cs_hho_builder_free(&(cs_hho_builders[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_sys_free(&(cs_hho_cell_sys[0]));
  cs_cell_builder_free(&(cs_hho_cell_bld[0]));
  cs_hho_builder_free(&(cs_hho_builders[0]));
#endif /* openMP */

  CS_FREE(cs_hho_cell_sys);
  CS_FREE(cs_hho_cell_bld);
  CS_FREE(cs_hho_builders);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_hho_stokes_t structure storing data useful for
 *         managing such a scheme
 *
 * \param[in] eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb   pointer to a cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_hho_stokes_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_hho_stokes_init_context(const cs_equation_param_t   *eqp,
                           cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != nullptr);
  if (eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Expected: scalar-valued HHO equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];
  const cs_lnum_t  n_cells = connect->n_cells;

  cs_hho_stokes_t *eqc = nullptr;

  CS_MALLOC(eqc, 1, cs_hho_stokes_t);

  /* Mesh flag to know what to build */
  eqb->msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ |
    CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ |
    CS_FLAG_COMP_EV | CS_FLAG_COMP_DIAM;

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_HHO_P0:
    eqc->n_cell_dofs = CS_N_DOFS_CELL_0TH;
    eqc->n_face_dofs = CS_N_DOFS_FACE_0TH;
    break;

  case CS_SPACE_SCHEME_HHO_P1:
    eqc->n_cell_dofs = CS_N_DOFS_CELL_1ST;
    eqc->n_face_dofs = CS_N_DOFS_FACE_1ST;
    break;

  case CS_SPACE_SCHEME_HHO_P2:
    eqc->n_cell_dofs = CS_N_DOFS_CELL_2ND;
    eqc->n_face_dofs = CS_N_DOFS_FACE_2ND;
    break;

    /* TODO: case CS_SPACE_SCHEME_HHO_PK */
  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);

  }

  /* System dimension */
  eqc->n_dofs = eqc->n_face_dofs * n_faces;
  eqc->n_max_loc_dofs = eqc->n_face_dofs*connect->n_max_fbyc + eqc->n_cell_dofs;

  /* TODO: Add a system helper */

  /* Values of each DoF related to the cells */
  const cs_lnum_t  n_cell_dofs = n_cells * eqc->n_cell_dofs;
  CS_MALLOC(eqc->cell_values, n_cell_dofs, cs_real_t);
  memset(eqc->cell_values, 0, sizeof(cs_real_t)*n_cell_dofs);

  /* Values at each face (interior and border) i.e. take into account BCs */
  CS_MALLOC(eqc->face_values, eqc->n_dofs, cs_real_t);
  memset(eqc->face_values, 0, sizeof(cs_real_t)*eqc->n_dofs);

  /* Source term */
  eqc->source_terms = nullptr;
  if (cs_equation_param_has_sourceterm(eqp)) {
    CS_MALLOC(eqc->source_terms, n_cell_dofs, cs_real_t);
    memset(eqc->source_terms, 0, sizeof(cs_real_t)*n_cell_dofs);

  } /* There is at least one source term */

  /* Members related to the static condensation.
     The transposed of acf_tilda is stored to speed-up the computation of
     the static condensation */
  CS_MALLOC(eqc->rc_tilda, n_cell_dofs, cs_real_t);
  memset(eqc->rc_tilda, 0, sizeof(cs_real_t)*n_cell_dofs);

  cs_lnum_t  n_row_blocks = connect->c2f->idx[n_cells];
  int       *row_block_sizes = nullptr;

  CS_MALLOC(row_block_sizes, n_row_blocks, int);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_row_blocks; i++)
    row_block_sizes[i] = eqc->n_face_dofs;

  int  col_block_size = eqc->n_cell_dofs;
  eqc->acf_tilda = cs_sdm_block_create(n_row_blocks, 1,
                                       row_block_sizes, &col_block_size);
  cs_sdm_block_init(eqc->acf_tilda,
                    n_row_blocks, 1,
                    row_block_sizes, &col_block_size);

  CS_FREE(row_block_sizes);

  /* Handle boundary conditions */
  const cs_lnum_t  n_b_faces = connect->n_faces[CS_BND_FACES];
  CS_MALLOC(eqc->bf2def_ids, n_b_faces, short int);

# pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    eqc->bf2def_ids[i] = -1; /* Default BC has no definition --> -1 */

  for (int def_id = 0; def_id < eqp->n_bc_defs; def_id++) {

    const cs_xdef_t  *def = eqp->bc_defs[def_id];
    const cs_zone_t  *bz = cs_boundary_zone_by_id(def->z_id);

#   pragma omp parallel for if (bz->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < bz->n_elts; i++)
      eqc->bf2def_ids[bz->elt_ids[i]] = def_id;

  } /* Loop on BC definitions */

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_hho_stokes_t structure
 *
 * \param[in, out]  data    pointer to a cs_hho_stokes_t structure
 *
 * \return a null pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_hho_stokes_free_context(void   *data)
{
  cs_hho_stokes_t  *eqc = (cs_hho_stokes_t *)data;

  if (eqc == nullptr)
    return eqc;

  /* Free temporary buffers */
  CS_FREE(eqc->cell_values);
  CS_FREE(eqc->face_values);
  CS_FREE(eqc->rc_tilda);
  CS_FREE(eqc->source_terms);
  CS_FREE(eqc->bf2def_ids);

  cs_sdm_free(eqc->acf_tilda);

  /* Last free */
  CS_FREE(eqc);

  return nullptr;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] eqb      pointer to a cs_equation_builder_t structure
 * \param[in, out] data     pointer to a cs_hho_stokes_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_stokes_compute_source(const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data)
{
  cs_hho_stokes_t  *eqc = (cs_hho_stokes_t *)data;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_st_dofs = eqc->n_cell_dofs * quant->n_cells;

  memset(eqc->source_terms, 0, sizeof(cs_real_t)*n_st_dofs);

  if (cs_equation_param_has_sourceterm(eqp) == false)
    return;

  cs_timer_t  t0 = cs_timer_time();

  /* TODO */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a HHO scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_hho_stokes_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_stokes_build_system(const cs_mesh_t            *mesh,
                           const cs_real_t            *field_val,
                           double                      dt_cur,
                           const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *data)
{
  CS_UNUSED(mesh);
  CS_UNUSED(field_val);
  CS_UNUSED(dt_cur);

  /* Sanity checks */
  assert(eqp != nullptr && eqb != nullptr);
  /* The only way to set a Dirichlet up to now */
  assert(eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED);

  /* Test to remove */
  if (cs_equation_param_has_convection(eqp))
    bft_error(__FILE__, __LINE__, 0,
              _(" Convection term is not handled yet.\n"));
  if (cs_equation_param_has_time(eqp))
    bft_error(__FILE__, __LINE__, 0,
              _(" Unsteady terms are not handled yet.\n"));

  cs_hho_stokes_t  *eqc = (cs_hho_stokes_t *)data;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_timer_t  t0 = cs_timer_time();

  /* TODO */
  CS_UNUSED(quant);
  CS_UNUSED(connect);
  CS_UNUSED(eqc);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values required for hybrid discretization
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_hho_stokes_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_stokes_update_field(const cs_real_t            *solu,
                           const cs_real_t            *rhs,
                           const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *data,
                           cs_real_t                  *field_val)
{
  CS_UNUSED(rhs);
  CS_UNUSED(eqp);
  CS_UNUSED(solu);
  CS_UNUSED(field_val);

  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_hho_stokes_t  *eqc = (cs_hho_stokes_t  *)data;

  CS_UNUSED(quant);
  CS_UNUSED(connect);
  CS_UNUSED(eqc);

  /* TODO */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at faces (DoF used in the linear system are
 *         located at primal faces)
 *
 * \param[in]  data    pointer to a data structure
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

double *
cs_hho_stokes_get_face_values(const void          *data)
{
  const cs_hho_stokes_t  *eqc = (const cs_hho_stokes_t  *)data;

  if (eqc == nullptr)
    return nullptr;
  else
    return eqc->face_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at cells (DoF used in the linear system are
 *         located at primal faces)
 *
 * \param[in]  data    pointer to a data structure
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

double *
cs_hho_stokes_get_cell_values(const void          *data)
{
  const cs_hho_stokes_t  *eqc = (const cs_hho_stokes_t  *)data;

  if (eqc == nullptr)
    return nullptr;
  else
    return eqc->cell_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_hho_stokes_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_stokes_extra_op(const char                 *eqname,
                       const cs_field_t           *field,
                       const cs_equation_param_t  *eqp,
                       cs_equation_builder_t      *eqb,
                       void                       *data)
{
  cs_hho_stokes_t  *eqc = (cs_hho_stokes_t  *)data;

  // TODO
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);
  CS_UNUSED(eqc);
  CS_UNUSED(field);
  CS_UNUSED(eqname);
}

/*----------------------------------------------------------------------------*/

#undef _dp3
END_C_DECLS
