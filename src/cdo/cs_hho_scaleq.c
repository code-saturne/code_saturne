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

#define CS_HHO_FACE_SIZE_0TH 1
#define CS_HHO_FACE_SIZE_1ST 3
#define CS_HHO_FACE_SIZE_2ND 6

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/* Algebraic system for HHO discretization */

struct _cs_hho_scaleq_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure.  */

  const cs_equation_param_t  *eqp;

  /* System size */
  cs_lnum_t  n_faces;
  cs_lnum_t  n_cells;
  cs_lnum_t  n_dofs;

  /* Local system size */
  int        n_max_loc_dofs;
  int        n_dofs_by_cell;
  int        n_dofs_by_face;

  /* Shortcut to know what to build */
  cs_flag_t    msh_flag;     // Information related to cell mesh
  cs_flag_t    bd_msh_flag;  // Information related to cell mesh (boundary)
  cs_flag_t    st_msh_flag;  // Information related to cell mesh (source term)
  cs_flag_t    sys_flag;     // Information related to the sytem

  /* Solution of the algebraic system at the last iteration */
  cs_real_t   *face_values;  /* DoF unknowns (x) + BCs */

  cs_real_t   *cell_rhs;  // right-hand side related to cell dofs
  cs_real_t   *acc_inv;   // inverse of a diagonal matrix (block cell-cell)
  cs_real_t   *acf;       /* Lower-Left block of the full matrix
                             (block cell-vertices). Access to the values thanks
                             to the c2f connectivity */

  /* Metadata related to associated properties */
  bool         diff_pty_uniform;
  bool         time_pty_uniform;
  bool         reac_pty_uniform[CS_CDO_N_MAX_REACTIONS];

  /* Source terms */
  cs_real_t    *source_terms; /* Array storing the value arising from the
                                 contribution of all source terms */
  cs_mask_t    *source_mask;  /* NULL if at least one source term is not
                                 defined for all cells (size = n_cells) */

  /* Pointer to functions which compute the value of the source term */
  cs_source_term_cellwise_t  *compute_source[CS_N_MAX_SOURCE_TERMS];

  /* Boundary conditions:

     face_bc should not change during the simulation.
     The case of a definition of the BCs which changes of type during the
     simulation is possible but not implemented.
     You just have to call the initialization step each time the type of BCs
     is modified to define an updated cs_cdo_bc_t structure.

     We translate Dirichlet BCs to border vertices
     The values can be modified at each time step in case of transient
     simulation.
     For Neumann and Robin BCs, the treatment is more complex since the
     contributions of these BCs to a dual cell related to a border vertex is
     computed from the contribution of different portions of primal faces
     (these primal border faces form a closure of the dual cell).
     These contributions are computed on the fly.

   */

  cs_cdo_bc_t           *face_bc; // list of faces sorted by type of BCs
  double                *dir_val; // size = vtx_dir->n_nhmg_elts

  /* In case of a weak enforcement of the Dirichlet BCs */
  cs_lnum_t             *c2bcbf_idx;  // size: n_cells + 1
  cs_lnum_t             *c2bcbf_ids;  // cell --> border faces ids

  /* Monitoring the efficiency */
  cs_equation_monitor_t           *monitor;

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


/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                  const cs_cdo_connect_t       *connect,
                                  const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to HHO schemes
 *
 * \param[in]  scheme_flag   flag to identify which kind of numerical scheme is
 *                           requested to solve the computational domain
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_initialize(cs_flag_t   scheme_flag)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  assert(connect != NULL);

  const int n_fc = connect->n_max_fbyc;

  int  order = -1, fbs = 0, cbs = 0;
  cs_space_scheme_t  space_scheme;

  if (scheme_flag & CS_SCHEME_FLAG_POLY2) {
    space_scheme = CS_SPACE_SCHEME_HHO_P2;
    fbs = CS_HHO_FACE_SIZE_2ND; // DoF by face
    cbs = cs_math_binom(5, 3);  // DoF for the cell
    order = 2;
  }
  else if (scheme_flag & CS_SCHEME_FLAG_POLY1) {
    space_scheme = CS_SPACE_SCHEME_HHO_P1;
    fbs = CS_HHO_FACE_SIZE_1ST; // DoF by face
    cbs = cs_math_binom(4, 3);  // DoF for the cell
    order = 1;
  }
  else {
    space_scheme = CS_SPACE_SCHEME_HHO_P0;
    fbs = CS_HHO_FACE_SIZE_0TH; // DoF by face
    cbs = 1;                    // DoF for the cell
    order = 0;
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

    cs_cell_builder_t  *cb = cs_cell_builder_create(space_scheme, connect);
    cs_hho_cell_bld[t_id] = cb;
    cs_hho_builders[t_id] = cs_hho_builder_create(order, n_fc);

    for (int i = 0; i < n_fc; i++) cb->ids[i] = fbs;
    cb->ids[n_fc] = cbs;

    cs_hho_cell_sys[t_id] = cs_cell_sys_create(n_dofs, n_fc + 1, cb->ids);
  }
#else
  assert(cs_glob_n_threads == 1);

    cs_cell_builder_t  *cb = cs_cell_builder_create(space_scheme, connect);
    cs_hho_cell_bld[0] = cb;
    cs_hho_builders[0] = cs_hho_builder_create(order, n_fc);

    for (int i = 0; i < n_fc; i++) cb->ids[i] = fbs;
    cb->ids[n_fc] = cbs;

    cs_hho_cell_sys[0] = cs_cell_sys_create(n_dofs, n_fc + 1, cb->ids);
#endif /* openMP */
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
 * \brief  Initialize a cs_hho_scaleq_t structure
 *
 * \param[in] eqp       pointer to a cs_equation_param_t structure
 * \param[in] mesh      pointer to a cs_mesh_t structure
 *
 * \return a pointer to a new allocated cs_hho_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_hho_scaleq_init(const cs_equation_param_t   *eqp,
                   const cs_mesh_t             *mesh)
{
  /* Sanity checks */
  assert(eqp != NULL);

  if (eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Expected: scalar-valued HHO equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_b_faces = connect->n_faces[1];
  const cs_lnum_t  n_faces = connect->n_faces[0];
  const cs_lnum_t  n_cells = connect->n_cells;

  cs_hho_scaleq_t  *b = NULL;

  BFT_MALLOC(b, 1, cs_hho_scaleq_t);

  /* Shared pointers */
  b->eqp = eqp;

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_HHO_P0:
    b->n_dofs_by_cell = 1;
    b->n_dofs_by_face = CS_HHO_FACE_SIZE_0TH;
    break;
  case CS_SPACE_SCHEME_HHO_P1:
    b->n_dofs_by_cell = 4;
    b->n_dofs_by_face = CS_HHO_FACE_SIZE_1ST;
    break;
  case CS_SPACE_SCHEME_HHO_P2:
    b->n_dofs_by_cell = 10;
    b->n_dofs_by_face = CS_HHO_FACE_SIZE_2ND;
    break;
    /* TODO: case CS_SPACE_SCHEME_HHO_PK */
  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);

  }

  /* System dimension */
  b->n_faces = n_faces;
  b->n_cells = n_cells;
  b->n_dofs = b->n_dofs_by_face * n_faces;
  b->n_max_loc_dofs = b->n_dofs_by_face*connect->n_max_fbyc + b->n_dofs_by_cell;

  /* Flag to indicate what to build in a cell mesh */
  b->msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_PFQ |
    CS_CDO_LOCAL_FE | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_EV |
    CS_CDO_LOCAL_DIAM;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  b->bd_msh_flag = 0;

  /* Set members and structures related to the management of the BCs
     Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.  */
  b->face_bc = cs_cdo_bc_define(eqp->default_bc,
                                eqp->n_bc_desc,
                                eqp->bc_desc,
                                n_b_faces);

  /* Values at each face (interior and border) i.e. take into account BCs */
  BFT_MALLOC(b->face_values, b->n_dofs, cs_real_t);
# pragma omp parallel for if (b->n_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < b->n_dofs; i++) b->face_values[i] = 0;

  /* Default intialization */
  b->st_msh_flag = 0;
  b->source_terms = NULL;

  if (cs_equation_param_has_sourceterm(eqp)) {

    /* Default intialization */
    b->st_msh_flag = cs_source_term_init(eqp->space_scheme,
                                         eqp->n_source_terms,
                                         (const cs_xdef_t **)eqp->source_terms,
                                         b->compute_source,
                                         &(b->sys_flag),
                                         &(b->source_mask));

    BFT_MALLOC(b->source_terms, n_cells * b->n_dofs_by_cell, cs_real_t);
#pragma omp parallel for if (n_cells * b->n_dofs_by_cell > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_cells * b->n_dofs_by_cell; i++)
      b->source_terms[i] = 0;

  } /* There is at least one source term */

  /* Monitoring */
  b->monitor = cs_equation_init_monitoring();

  return b;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_hho_scaleq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_hho_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_hho_scaleq_free(void   *builder)
{
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t *)builder;

  if (b == NULL)
    return b;

  /* eqp is only shared. This structure is freed later. */

  /* Free BC structure */
  b->face_bc = cs_cdo_bc_free(b->face_bc);

  /* Free temporary buffers */
  BFT_FREE(b->source_terms);
  BFT_FREE(b->face_values);

  /* Monitoring structure */
  BFT_FREE(b->monitor);

  /* Last free */
  BFT_FREE(b);
  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in, out] builder     pointer to a cs_hho_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_compute_source(void            *builder)
{
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t *)builder;

  const cs_equation_param_t  *eqp = b->eqp;

  if (eqp->n_source_terms == 0)
    return;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         builder structure
 *
 * \param[in, out] builder        pointer to generic builder structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_initialize_system(void           *builder,
                                cs_matrix_t   **system_matrix,
                                cs_real_t     **system_rhs)
{
  if (builder == NULL)
    return;
  assert(*system_matrix == NULL && *system_rhs == NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_timer_t  t0 = cs_timer_time();
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t *)builder;

  cs_matrix_structure_t  *ms = NULL;
  cs_lnum_t  n_elts = 0;

  /* Create the matrix related to the current algebraic system */
  switch (b->eqp->space_scheme) {

  case CS_SPACE_SCHEME_HHO_P0:
    ms = cs_equation_get_matrix_structure(CS_SPACE_SCHEME_HHO_P0);
    n_elts = quant->n_faces;
    break;

  case CS_SPACE_SCHEME_HHO_P1:
    ms = cs_equation_get_matrix_structure(CS_SPACE_SCHEME_HHO_P1);
    n_elts = CS_HHO_FACE_SIZE_1ST * quant->n_faces;
    break;

  case CS_SPACE_SCHEME_HHO_P2:
    ms = cs_equation_get_matrix_structure(CS_SPACE_SCHEME_HHO_P2);
    n_elts = CS_HHO_FACE_SIZE_2ND * quant->n_faces;
    break;

  default:
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
  cs_timer_counter_add_diff(&(b->monitor->tcb), &t0, &t1);
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
 * \param[in, out] builder    pointer to cs_hho_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_build_system(const cs_mesh_t         *mesh,
                           const cs_real_t         *field_val,
                           double                   dt_cur,
                           void                    *builder,
                           cs_real_t               *rhs,
                           cs_matrix_t             *matrix)
{
  /* Sanity checks */
  assert(rhs != NULL && matrix != NULL);

  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t *)builder;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)         \
  shared(dt_cur, quant, connect, b, rhs, matrix, mav,  \
         field_val, cs_hho_cell_sys, cs_hho_cell_bld, cs_hho_builders)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif
    const cs_equation_param_t  *eqp = b->eqp;

    /* Test to remove */
    if (cs_equation_param_has_convection(eqp))
      bft_error(__FILE__, __LINE__, 0,
                _(" Convection term is not handled yet.\n"));
    if (cs_equation_param_has_time(eqp))
      bft_error(__FILE__, __LINE__, 0,
                _(" Unsteady terms are not handled yet.\n"));

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_hho_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_hho_cell_bld[t_id];
    cs_hho_builder_t  *hhob = cs_hho_builders[t_id];

    /* Set inside the OMP section so that each thread has its own value */

    /* Initialize members of the builder related to the current system
       Preparatory step for diffusion term */
    if (cs_equation_param_has_diffusion(eqp)) {
      if (b->diff_pty_uniform) {

        cs_property_get_cell_tensor(0, // cell_id
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    cb->pty_mat);

        if (eqp->diffusion_hodge.is_iso)
          cb->pty_val = cb->pty_mat[0][0];

      } /* Diffusion property is uniform */
    } /* Diffusion */

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag =
        cs_equation_get_cell_mesh_flag(cell_flag,
                                       b->msh_flag,
                                       b->bd_msh_flag,
                                       b->st_msh_flag);

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

      const short int  face_offset = cm->n_fc + b->n_dofs_by_face;

      /* DIFFUSION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ============================================== */

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        if (!(b->diff_pty_uniform)) {

          cs_property_tensor_in_cell(cm,
                                     eqp->diffusion_property,
                                     eqp->diffusion_hodge.inv_pty,
                                     cb->pty_mat);

          if (eqp->diffusion_hodge.is_iso)
            cb->pty_val = cb->pty_mat[0][0];

        }

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
                                        b->source_mask,
                                        b->compute_source,
                                        hhob,
                                        cb,    // mass matrix is cb->hdg
                                        csys); // Fill csys->source

        if (cs_equation_param_has_time(eqp)) {

          /* Same strategy as if one applies a implicit scheme */
          cs_real_t  *_rhs = csys->rhs + face_offset;
          const cs_real_t  *_st = csys->source + face_offset;
          for (int i = 0; i < b->n_dofs_by_cell; i++)
            _rhs[i] += _st[i];

        }

        /* Reset the value of the source term for the cell DoF
           Source term is only hold by the cell DoF in face-based schemes */
        {
          cs_real_t  *st = b->source_terms + c_id * b->n_dofs_by_cell;
          const cs_real_t  *_st = csys->source + face_offset;
          for (int i = 0; i < b->n_dofs_by_cell; i++)
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
         a matrix of size n_vc. Store information in the builder structure in
         order to be able to compute the values at cell centers. */
      /* _condense_and_store(connect->c2f, b, cb, csys); */

    } // Main loop on cells

  } // OPENMP Block

  cs_matrix_assembler_values_done(mav); // optional

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 2
  cs_dump_array_to_listing("FINAL RHS_FACE", b->n_dofs, rhs, b->n_dofs_by_face);
  if (b->source_terms != NULL)
    cs_dump_array_to_listing("FINAL RHS_CELL",
                             quant->n_cells,
                             b->source_terms, b->n_dofs_by_cell);
#endif

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(b->monitor->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values required for hybrid discretization
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in, out] builder    pointer to builder structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_update_field(const cs_real_t     *solu,
                           const cs_real_t     *rhs,
                           void                *builder,
                           cs_real_t           *field_val)
{
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t  *)builder;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at faces (DoF used in the linear system are
 *         located at primal faces)
 *
 * \param[in]  builder    pointer to a builder structure
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

double *
cs_hho_scaleq_get_face_values(const void          *builder)
{
  const cs_hho_scaleq_t  *b = (const cs_hho_scaleq_t  *)builder;

  if (b == NULL)
    return NULL;
  else
    return b->face_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in, out]  builder    pointer to builder structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_extra_op(const char            *eqname,
                       const cs_field_t      *field,
                       void                  *builder)
{
  cs_hho_scaleq_t  *b = (cs_hho_scaleq_t  *)builder;

  const cs_equation_param_t  *eqp = b->eqp;

  // TODO
  CS_UNUSED(eqp);
  CS_UNUSED(field);
  CS_UNUSED(eqname);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
