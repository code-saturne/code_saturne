/*============================================================================
 * Build an algebraic system for scalar conv./diff. eq. with Hybrid High Order
 * space discretization
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_equation_common.h"
#include "cs_hho_builder.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_scheme_geometry.h"
#include "cs_search.h"
#include "cs_sla.h"
#include "cs_source_term.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hho_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_HHO_SCALEQ_DBG  0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/* Algebraic system for HHO discretization */

struct _cs_hho_scaleq_t {

  /* System size */
  cs_lnum_t  n_faces;
  cs_lnum_t  n_cells;
  cs_lnum_t  n_dofs;

  /* Local system size */
  int        n_max_loc_dofs;
  int        n_dofs_by_cell;
  int        n_dofs_by_face;

  /* Solution of the algebraic system at the last iteration */
  cs_real_t   *face_values;  /* DoF unknowns (x) + BCs */

  cs_real_t   *cell_rhs;  // right-hand side related to cell dofs
  cs_real_t   *acc_inv;   // inverse of a diagonal matrix (block cell-cell)
  cs_real_t   *acf;       /* Lower-Left block of the full matrix
                             (block cell-vertices). Access to the values thanks
                             to the c2f connectivity */

  /* Array storing the value arising from the contribution of all source
     terms */
  cs_real_t   *source_terms;

  /* Boundary conditions */
  double                *dir_val; // size = vtx_dir->n_nhmg_elts

  /* In case of a weak enforcement of the Dirichlet BCs */
  cs_lnum_t             *c2bcbf_idx;  // size: n_cells + 1
  cs_lnum_t             *c2bcbf_ids;  // cell --> border faces ids

};

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */
static cs_cell_sys_t  **cs_hho_cell_sys = NULL;
static cs_cell_builder_t  **cs_hho_cell_bld = NULL;
static cs_hho_builder_t  **cs_hho_builders = NULL;

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
_cell_builder_create(cs_space_scheme_t         space_scheme,
                     const cs_cdo_connect_t   *connect)
{
  int  size;

  const int  n_fc = connect->n_max_fbyc;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  switch (space_scheme) {

  case CS_SPACE_SCHEME_HHO_P0:  /* TODO */
    {
      BFT_MALLOC(cb->ids, n_fc, short int);
      memset(cb->ids, 0, n_fc*sizeof(short int));

      size = CS_MAX(32, n_fc*(n_fc+1));
      BFT_MALLOC(cb->values, size, double);
      memset(cb->values, 0, size*sizeof(cs_real_t));

      size = CS_MAX(2*n_fc, 15);
      BFT_MALLOC(cb->vectors, size, cs_real_3_t);
      memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

      /* Local square dense matrices used during the construction of
         operators */
      cb->hdg = cs_sdm_square_create(n_fc);
      cb->loc = cs_sdm_square_create(n_fc + 1);
      cb->aux = cs_sdm_square_create(n_fc + 1);
    }
    break;

  case CS_SPACE_SCHEME_HHO_P1:
    {
      /* Store the block size description */
      size = n_fc + 1;
      BFT_MALLOC(cb->ids, size, short int);
      memset(cb->ids, 0, size*sizeof(short int));

      /* Store the face, cell and gradient basis function evaluations and
         the Gauss point weights
         --->  42 = 3 + 4 + 3*10 + 5
         or the factorization of the stiffness matrix used for computing
         the gradient reconstruction operator
         --->  n = 9 (10 - 1) --> facto = n*(n+1)/2 = 45
                              --> tmp buffer for facto = n = 9
                              --> 45 + 9 = 54
      */
      size = 54;
      BFT_MALLOC(cb->values, size, double);
      memset(cb->values, 0, size*sizeof(cs_real_t));

      /* Store Gauss points and tensor.n_f products */
      size = CS_MAX(15, 5 + n_fc);
      BFT_MALLOC(cb->vectors, size, cs_real_3_t);
      memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

      /* Local dense matrices used during the construction of operators */
      const short int g_size = 9;                   /* basis (P_(k+1)) - 1 */
      for (int i = 0; i < n_fc; i++) cb->ids[i] = 3;
      cb->ids[n_fc] = 4;

      short int  _sizes[3] = {1, 3, 6}; /* c0, cs-1, cs_kp1 - cs */
      cb->hdg = cs_sdm_block_create(1, 3, &g_size, _sizes);
      cb->loc = cs_sdm_block_create(n_fc + 1, n_fc + 1, cb->ids, cb->ids);
      cb->aux = cs_sdm_block_create(n_fc + 1, 1, cb->ids, &g_size); /* R_g^T */
    }
    break;

  case CS_SPACE_SCHEME_HHO_P2:
    {
      /* Store the block size description */
      size = n_fc + 1;
      BFT_MALLOC(cb->ids, size, short int);
      memset(cb->ids, 0, size*sizeof(short int));

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
      size = 209;
      BFT_MALLOC(cb->values, size, double);
      memset(cb->values, 0, size*sizeof(cs_real_t));

      size = 15 + n_fc;  /* Store Gauss points and tensor.n_f products */
      BFT_MALLOC(cb->vectors, size, cs_real_3_t);
      memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

      /* Local dense matrices used during the construction of operators */
      const short int g_size = 19; /* basis (P_(k+1)) - 1 */
      for (int i = 0; i < n_fc; i++) cb->ids[i] = 6;
      cb->ids[n_fc] = 10;

      short int  _sizes[3] = {1, 9, 10}; /* c0, cs-1, cs_kp1 - cs */
      cb->hdg = cs_sdm_block_create(1, 3, &g_size, _sizes);
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
cs_hho_scaleq_initialize(cs_flag_t                      scheme_flag,
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
  cs_space_scheme_t  space_scheme;

  if (scheme_flag & CS_SCHEME_FLAG_POLY2) {

    space_scheme = CS_SPACE_SCHEME_HHO_P2;
    fbs = CS_N_FACE_DOFS_2ND; // DoF by face
    cbs = CS_N_CELL_DOFS_2ND; // DoF for the cell
    order = 2;
    assert(ma2 != NULL && ms2 != NULL);

  }
  else if (scheme_flag & CS_SCHEME_FLAG_POLY1) {

    space_scheme = CS_SPACE_SCHEME_HHO_P1;
    fbs = CS_N_FACE_DOFS_1ST; // DoF by face
    cbs = CS_N_CELL_DOFS_1ST;  // DoF for the cell
    order = 1;
    assert(ma1 != NULL && ms1 != NULL);

  }
  else {

    space_scheme = CS_SPACE_SCHEME_HHO_P0;
    fbs = CS_N_FACE_DOFS_0TH; // DoF by face
    cbs = CS_N_CELL_DOFS_0TH; // DoF for the cell
    order = 0;
    assert(ma0 != NULL && ms0 != NULL);

  }

  const int n_dofs = n_fc * fbs + cbs;

    /* Structure used to build the final system by a cell-wise process */
  assert(cs_glob_n_threads > 0);  /* Sanity check */

  BFT_MALLOC(cs_hho_cell_bld, cs_glob_n_threads, cs_cell_builder_t *);
  BFT_MALLOC(cs_hho_cell_sys, cs_glob_n_threads, cs_cell_sys_t *);

  /* Allocate builder structure specific to HHO schemes. This is an additional
     structure with respect to a cs_cell_builder_t structure */
  BFT_MALLOC(cs_hho_builders, cs_glob_n_threads, cs_hho_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_hho_cell_bld[i] = NULL;
    cs_hho_cell_sys[i] = NULL;
    cs_hho_builders[i] = NULL;
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
cs_hho_scaleq_get(cs_cell_sys_t       **csys,
                  cs_cell_builder_t   **cb,
                  cs_hho_builder_t    **hhob)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

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
cs_hho_scaleq_finalize(void)
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

  BFT_FREE(cs_hho_cell_sys);
  BFT_FREE(cs_hho_cell_bld);
  BFT_FREE(cs_hho_builders);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_hho_scaleq_t structure storing data useful for
 *         managing such a scheme
 *
 * \param[in] eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb   pointer to a cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_hho_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_hho_scaleq_init_context(const cs_equation_param_t   *eqp,
                           cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL);
  if (eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Expected: scalar-valued HHO equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_faces = connect->n_faces[0];
  const cs_lnum_t  n_cells = connect->n_cells;

  cs_hho_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_hho_scaleq_t);

  /* Mesh flag to know what to build */
  eqb->msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ |
    CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ |
    CS_CDO_LOCAL_EV | CS_CDO_LOCAL_DIAM;

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_HHO_P0:
    eqc->n_dofs_by_cell = CS_N_CELL_DOFS_0TH;
    eqc->n_dofs_by_face = CS_N_FACE_DOFS_0TH;
    break;
  case CS_SPACE_SCHEME_HHO_P1:
    eqc->n_dofs_by_cell = CS_N_CELL_DOFS_1ST;
    eqc->n_dofs_by_face = CS_N_FACE_DOFS_1ST;
    break;
  case CS_SPACE_SCHEME_HHO_P2:
    eqc->n_dofs_by_cell = CS_N_CELL_DOFS_2ND;
    eqc->n_dofs_by_face = CS_N_FACE_DOFS_2ND;
    break;
    /* TODO: case CS_SPACE_SCHEME_HHO_PK */
  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);

  }

  /* System dimension */
  eqc->n_faces = n_faces;
  eqc->n_cells = n_cells;
  eqc->n_dofs = eqc->n_dofs_by_face * n_faces;
  eqc->n_max_loc_dofs =
    eqc->n_dofs_by_face*connect->n_max_fbyc + eqc->n_dofs_by_cell;

  /* Values at each face (interior and border) i.e. take into account BCs */
  BFT_MALLOC(eqc->face_values, eqc->n_dofs, cs_real_t);
# pragma omp parallel for if (eqc->n_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < eqc->n_dofs; i++) eqc->face_values[i] = 0;

  /* Source term */
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, n_cells * eqc->n_dofs_by_cell, cs_real_t);
#pragma omp parallel for if (n_cells * eqc->n_dofs_by_cell > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells * eqc->n_dofs_by_cell; i++)
      eqc->source_terms[i] = 0;

  } /* There is at least one source term */

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_hho_scaleq_t structure
 *
 * \param[in, out]  data    pointer to a cs_hho_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_hho_scaleq_free_context(void   *data)
{
  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t *)data;

  if (eqc == NULL)
    return eqc;

  /* Free temporary buffers */
  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->face_values);

  /* Last free */
  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] eqb      pointer to a cs_equation_builder_t structure
 * \param[in, out] data     pointer to a cs_hho_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_compute_source(const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data)
{
  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t *)data;

  if (eqp->n_source_terms == 0)
    return;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         data structure
 *
 * \param[in]      eqp            pointer to a cs_equation_param_t structure
 * \param[in, out] eqb            pointer to a cs_equation_builder_t structure
 * \param[in, out] data           pointer to generic data structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_initialize_system(const cs_equation_param_t  *eqp,
                                cs_equation_builder_t      *eqb,
                                void                       *data,
                                cs_matrix_t               **system_matrix,
                                cs_real_t                 **system_rhs)
{
  CS_UNUSED(data);
  assert(*system_matrix == NULL && *system_rhs == NULL);

  const cs_matrix_structure_t  *ms;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_timer_t  t0 = cs_timer_time();
  cs_lnum_t  n_elts = 0;

  /* Create the matrix related to the current algebraic system */
  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_HHO_P0:
    ms = cs_shared_ms0;
    n_elts = quant->n_faces;
    break;

  case CS_SPACE_SCHEME_HHO_P1:
    ms = cs_shared_ms1;
    n_elts = CS_N_FACE_DOFS_1ST * quant->n_faces;
    break;

  case CS_SPACE_SCHEME_HHO_P2:
    ms = cs_shared_ms2;
    n_elts = CS_N_FACE_DOFS_2ND * quant->n_faces;
    break;

  default:
    ms = NULL;
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid space discretization.", __func__);
    break;

  }

  *system_matrix = cs_matrix_create(ms);

  /* Allocate and initialize the related right-hand side */
  BFT_MALLOC(*system_rhs, n_elts, cs_real_t);
# pragma omp parallel for if  (quant->n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) (*system_rhs)[i] = 0.0;

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
 * \param[in, out] data       pointer to cs_hho_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_build_system(const cs_mesh_t            *mesh,
                           const cs_real_t            *field_val,
                           double                      dt_cur,
                           const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *data,
                           cs_real_t                  *rhs,
                           cs_matrix_t                *matrix)
{
  /* Sanity checks */
  assert(rhs != NULL && matrix != NULL && eqp != NULL && eqb != NULL);

  /* Test to remove */
  if (cs_equation_param_has_convection(eqp))
    bft_error(__FILE__, __LINE__, 0,
              _(" Convection term is not handled yet.\n"));
  if (cs_equation_param_has_time(eqp))
    bft_error(__FILE__, __LINE__, 0,
              _(" Unsteady terms are not handled yet.\n"));

  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t *)data;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav,        \
         field_val, cs_hho_cell_sys, cs_hho_cell_bld, cs_hho_builders)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_hho_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_hho_cell_bld[t_id];
    cs_hho_builder_t  *hhob = cs_hho_builders[t_id];

    /* Set inside the OMP section so that each thread has its own value */

    /* Initialization of the values of properties */
    double  time_pty_val = 1.0;
    double  reac_pty_vals[CS_CDO_N_MAX_REACTIONS];

    cs_equation_init_properties(eqp, eqb, &time_pty_val, reac_pty_vals, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag =
        cs_equation_get_cell_mesh_flag(cell_flag, eqb);

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Define set of basis functions for cell faces and the current cell */
      cs_hho_builder_cellwise_setup(cm, cb, hhob);

      /* Set the local (i.e. cellwise) structures for the current cell */
      /* _init_cell_structures(cell_flag, cm, b, */
      /*                       dir_values, neu_tags, field_val,  // in */
      /*                       csys, cbc, cb);                   // out */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 2
      if (c_id % CS_HHO_SCALEQ_MODULO == 0) cs_cell_mesh_dump(cm);
#endif

      const short int  face_offset = cm->n_fc + eqc->n_dofs_by_face;

      /* DIFFUSION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ============================================== */

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        if (!(eqb->diff_pty_uniform))
          cs_equation_set_diffusion_property_cw(eqp, cm, cell_flag, cb);

        cs_hho_builder_compute_grad_reco(cm, cb, hhob);

        // local matrix owned by the cellwise builder (store in cb->loc)
        cs_hho_builder_diffusion(cm, cb, hhob);

        // Add the local diffusion operator to the local system
        cs_sdm_block_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 1
        if (c_id % CS_HHO_SCALEQ_MODULO == 0)
          cs_cell_sys_dump("\n>> Local system after diffusion", c_id, csys);
#endif
      } /* END OF DIFFUSION */

      /* SOURCE TERM COMPUTATION */
      /* ======================= */

      if (cs_equation_param_has_sourceterm(eqp)) {

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */
        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (const cs_xdef_t **)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        hhob,
                                        cb,    // mass matrix is cb->hdg
                                        csys); // Fill csys->source

        if (cs_equation_param_has_time(eqp) == false) {

          /* Same strategy as if one applies a implicit scheme */
          cs_real_t  *_rhs = csys->rhs + face_offset;
          const cs_real_t  *_st = csys->source + face_offset;
          for (int i = 0; i < eqc->n_dofs_by_cell; i++)
            _rhs[i] += _st[i];

        }

        /* Reset the value of the source term for the cell DoF
           Source term is only hold by the cell DoF in face-based schemes */
        {
          cs_real_t  *st = eqc->source_terms + c_id * eqc->n_dofs_by_cell;
          const cs_real_t  *_st = csys->source + face_offset;
          for (int i = 0; i < eqc->n_dofs_by_cell; i++)
            st[i] = _st[i];
        }

      } /* End of term source contribution */

      /* TODO: Neumann boundary conditions */
      /* /\* Computing Dirichlet BC *\/ */
      /* if (connect->c_info->flag[c_id] & CS_CDO_CONNECT_BD) */
      /*   if (cs_hho_compute_dir_bc(cm,cb,dt_cur,bdl) > 0) */
      /*     // If there are Dirichlet faces, enforce penalization */
      /*     cs_hho_enforce_penalization(bdl); */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 1
      if (c_id % CS_HHO_SCALEQ_MODULO == 0)
        cs_cell_sys_dump(">> Local system matrix before condensation",
                         c_id, csys);
#endif

      /* Static condensation of the local system matrix of size n_vc + 1 into
         a matrix of size n_vc. Store information in the input structure in
         order to be able to compute the values at cell centers. */
      /* _condense_and_store(connect->c2f, b, cb, csys); */

    } // Main loop on cells

  } // OPENMP Block

  cs_matrix_assembler_values_done(mav); // optional

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 2
  cs_dump_array_to_listing("FINAL RHS_FACE", eqc->n_dofs, rhs,
                           eqc->n_dofs_by_face);
  if (eqc->source_terms != NULL)
    cs_dump_array_to_listing("FINAL RHS_CELL",
                             quant->n_cells,
                             eqc->source_terms, eqc->n_dofs_by_cell);
#endif

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
 * \param[in, out] data       pointer to cs_hho_scaleq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_update_field(const cs_real_t            *solu,
                           const cs_real_t            *rhs,
                           const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *data,
                           cs_real_t                  *field_val)
{
  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t  *)data;
  cs_timer_t  t0 = cs_timer_time();

  /* Set the computed solution in field array */
# pragma omp parallel for if (eqc->n_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < eqc->n_dofs; i++)
    field_val[i] = solu[i];

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
cs_hho_scaleq_get_face_values(const void          *data)
{
  const cs_hho_scaleq_t  *eqc = (const cs_hho_scaleq_t  *)data;

  if (eqc == NULL)
    return NULL;
  else
    return eqc->face_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_hho_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_extra_op(const char                 *eqname,
                       const cs_field_t           *field,
                       const cs_equation_param_t  *eqp,
                       cs_equation_builder_t      *eqb,
                       void                       *data)
{
  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t  *)data;

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
