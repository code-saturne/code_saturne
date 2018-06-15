/*============================================================================
 * Build an algebraic CDO face-based system for unsteady convection/diffusion
 * reaction of scalar-valued equations with source terms
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_search.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_cdofb_priv.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_SCALEQ_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */
static cs_cell_sys_t      **cs_cdofb_cell_sys = NULL;
static cs_cell_builder_t  **cs_cdofb_cell_bld = NULL;

/* Pointer to shared structures */
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;
static const cs_matrix_structure_t  *cs_shared_ms;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local builder structure used for building the system
 *          cellwise
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_cell_builder_t *
_cell_builder_create(const cs_cdo_connect_t   *connect)
{
  const int  n_fc = connect->n_max_fbyc;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  BFT_MALLOC(cb->ids, n_fc, short int);
  memset(cb->ids, 0, n_fc*sizeof(short int));

  int  size = n_fc*(n_fc+1);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(cs_real_t));

  size = 2*n_fc;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of
     operators */
  cb->hdg = cs_sdm_square_create(n_fc);
  cb->loc = cs_sdm_square_create(n_fc + 1);
  cb->aux = cs_sdm_square_create(n_fc + 1);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      cell_flag   flag related to the current cell
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         pointer to a cs_cdofb_scaleq_t structure
 * \param[in]      dir_values  Dirichlet values associated to each face
 * \param[in]      neu_tags    definition id related to each Neumann face
 * \param[in]      field_tn    values of the field at the last computed time
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_init_cell_system(const cs_flag_t               cell_flag,
                  const cs_cell_mesh_t         *cm,
                  const cs_equation_param_t    *eqp,
                  const cs_equation_builder_t  *eqb,
                  const cs_cdofb_scaleq_t      *eqc,
                  const cs_real_t               dir_values[],
                  const short int               neu_tags[],
                  const cs_real_t               field_tn[],
                  cs_real_t                     t_eval,
                  cs_cell_sys_t                *csys,
                  cs_cell_builder_t            *cb)
{
  CS_UNUSED(cb);

  /* Cell-wise view of the linear system to build */
  const int  n_dofs = cm->n_fc + 1;

  /* Initialize the local system */
  cs_cell_sys_reset(cell_flag, n_dofs, cm->n_fc, csys);

  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;
  cs_sdm_square_init(n_dofs, csys->mat);
  for (short int f = 0; f < cm->n_fc; f++) {
    csys->dof_ids[f] = cm->f_ids[f];
    csys->val_n[f] = eqc->face_values[cm->f_ids[f]];
  }

  csys->dof_ids[cm->n_fc] = cm->c_id;
  csys->val_n[cm->n_fc] = field_tn[cm->c_id];

  /* Update rhs with the previous computation of source term if needed */
  if (cs_equation_param_has_sourceterm(eqp)) {
    if (cs_equation_param_has_time(eqp)) {

      /* Source terms attached to cells: Need to update rhs because the part
         related to cell is used in the static condensation */
      cs_cdo_time_update_rhs(eqp,
                             1, /* stride */
                             1, /* n_dofs */
                             csys->dof_ids + cm->n_fc,
                             eqc->source_terms,
                             csys->rhs + cm->n_fc);


    }
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY) {

    const cs_cdo_connect_t  *connect = cs_shared_connect;

    /* Identify which face is a boundary face */
    for (short int f = 0; f < cm->n_fc; f++) {

      const cs_lnum_t  bf_id = cm->f_ids[f] - connect->n_faces[2]; // n_i_faces
      if (bf_id > -1) // Border face
        cs_equation_fb_set_cell_bc(bf_id, f,
                                   eqb->face_bc->flag[bf_id],
                                   cm,
                                   connect,
                                   cs_shared_quant,
                                   eqp,
                                   dir_values,
                                   neu_tags,
                                   t_eval,
                                   csys,
                                   cb);

    } // Loop on cell faces

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    for (short int f = 0; f < cm->n_fc; f++) {
      if (csys->dof_flag[f] & CS_CDO_BC_HMG_DIRICHLET)
        if (fabs(csys->dir_values[f]) > 10*DBL_MIN)
          bft_error(__FILE__, __LINE__, 0,
                    "Invalid enforcement of Dirichlet BCs on faces");
    }
#endif

  } /* Border cell */

}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to CDO
 *         scalar-valued face-based schemes.
 *         Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_init_common(const cs_cdo_quantities_t     *quant,
                            const cs_cdo_connect_t        *connect,
                            const cs_time_step_t          *time_step,
                            const cs_matrix_structure_t   *ms)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
  cs_shared_ms = ms;

  /* Specific treatment for handling openMP */
  BFT_MALLOC(cs_cdofb_cell_sys, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(cs_cdofb_cell_bld, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_cdofb_cell_sys[i] = NULL;
    cs_cdofb_cell_bld[i] = NULL;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cdofb_cell_sys[t_id] = cs_cell_sys_create(connect->n_max_fbyc + 1,
                                                 connect->n_max_fbyc,
                                                 1, NULL);
    cs_cdofb_cell_bld[t_id] = _cell_builder_create(connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cdofb_cell_sys[0] = cs_cell_sys_create(connect->n_max_fbyc + 1,
                                            connect->n_max_fbyc,
                                            1, NULL);
  cs_cdofb_cell_bld[0] = _cell_builder_create(connect);
#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  *csys = cs_cdofb_cell_sys[t_id];
  *cb = cs_cdofb_cell_bld[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO face-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_finalize_common(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(cs_cdofb_cell_sys[t_id]));
    cs_cell_builder_free(&(cs_cdofb_cell_bld[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(cs_cdofb_cell_sys[0]));
  cs_cell_builder_free(&(cs_cdofb_cell_bld[0]));
#endif /* openMP */

  BFT_FREE(cs_cdofb_cell_sys);
  BFT_FREE(cs_cdofb_cell_bld);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_scaleq_t structure storing data useful for
 *         managing such a scheme
 *
 * \param[in]      eqp    pointer to a cs_equation_param_t structure
 * \param[in, out] eqb    pointer to a cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_scaleq_init_context(const cs_equation_param_t   *eqp,
                             cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOFB && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued CDO face-based equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[0];

  cs_cdofb_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdofb_scaleq_t);

  /* Dimensions of the algebraic system */
  eqc->n_dofs = n_faces + n_cells;

  eqb->msh_flag = CS_CDO_LOCAL_PF | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_CDO_LOCAL_EV | CS_CDO_LOCAL_EF | CS_CDO_LOCAL_EFQ;

  /* Set members and structures related to the management of the BCs
     Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.  */

  /* Values at each face (interior and border) i.e. take into account BCs */
  BFT_MALLOC(eqc->face_values, n_faces, cs_real_t);
# pragma omp parallel for if (n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_faces; i++) eqc->face_values[i] = 0;

  /* Store the last computed values of the field at cell centers and the data
     needed to compute the cell values from the face values.
     No need to synchronize all these quantities since they are only cellwise
     quantities. */
  BFT_MALLOC(eqc->rc_tilda, n_cells, cs_real_t);
  BFT_MALLOC(eqc->acf_tilda, connect->c2f->idx[n_cells], cs_real_t);

  memset(eqc->rc_tilda, 0, sizeof(cs_real_t)*n_cells);
  memset(eqc->acf_tilda, 0, sizeof(cs_real_t)*connect->c2f->idx[n_cells]);

  /* Diffusion part */
  /* -------------- */

  eqc->get_stiffness_matrix = NULL;
  eqc->boundary_flux_op = NULL;
  eqc->enforce_dirichlet = NULL;

  if (cs_equation_param_has_diffusion(eqp)) {

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      eqc->get_stiffness_matrix = cs_hodge_fb_cost_get_stiffness;
      eqc->boundary_flux_op = NULL; //cs_cdovb_diffusion_cost_flux_op;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqc->get_stiffness_matrix = cs_hodge_fb_voro_get_stiffness;
      eqc->boundary_flux_op = NULL; //cs_cdovb_diffusion_cost_flux_op;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                (" Invalid type of algorithm to build the diffusion term."));

    } // Switch on Hodge algo.

    switch (eqp->enforcement) {

    case CS_PARAM_BC_ENFORCE_WEAK_PENA:
      eqc->enforce_dirichlet = cs_cdo_diffusion_pena_dirichlet;
      break;

    case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    default:
      bft_error(__FILE__, __LINE__, 0,
                (" Invalid type of algorithm to enforce Dirichlet BC."));

    }

  } // Diffusion part

  /* Advection part */

  eqc->get_advection_matrix = NULL;
  eqc->add_advection_bc = NULL;
  // TODO

  /* Time part */
  if (cs_equation_param_has_time(eqp))
    eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
  eqc->apply_time_scheme = cs_cdo_time_get_scheme_function(eqb->sys_flag, eqp);

  /* Source term part */
  eqc->source_terms = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, n_cells, cs_real_t);
    memset(eqc->source_terms, 0, sizeof(cs_real_t)*n_cells);

  } /* There is at least one source term */

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_scaleq_t structure
 *
 * \param[in, out]  data   pointer to a cs_cdofb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_scaleq_free_context(void   *data)
{
  cs_cdofb_scaleq_t   *eqc  = (cs_cdofb_scaleq_t *)data;

  if (eqc == NULL)
    return eqc;

  /* Free temporary buffers */
  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->face_values);
  BFT_FREE(eqc->rc_tilda);
  BFT_FREE(eqc->acf_tilda);

  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         data structure
 *
 * \param[in]      eqp            pointer to a cs_equation_param_t structure
 * \param[in, out] eqb            pointer to a cs_equation_builder_t structure
 * \param[in, out] data           pointer to cs_cdofb_scaleq_t structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_initialize_system(const cs_equation_param_t  *eqp,
                                  cs_equation_builder_t      *eqb,
                                  void                       *data,
                                  cs_matrix_t               **system_matrix,
                                  cs_real_t                 **system_rhs)
{
  assert(eqb != NULL);
  assert(*system_matrix == NULL && *system_rhs == NULL);
  CS_UNUSED(data);
  CS_UNUSED(eqp);

  cs_timer_t  t0 = cs_timer_time();

  /* Create the matrix related to the current algebraic system */
  *system_matrix = cs_matrix_create(cs_shared_ms);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Allocate and initialize the related right-hand side */
  BFT_MALLOC(*system_rhs, quant->n_faces, cs_real_t);
# pragma omp parallel for if  (quant->n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < quant->n_faces; i++) (*system_rhs)[i] = 0.0;

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings when the fields
 *         stem from a scalar CDO face-based scheme.
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in, out] eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval      time at which one evaluates BCs
 * \param[in, out] field_val   pointer to the values of the variable field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_set_dir_bc(const cs_mesh_t              *mesh,
                           const cs_equation_param_t    *eqp,
                           cs_equation_builder_t        *eqb,
                           cs_real_t                     t_eval,
                           cs_real_t                     field_val[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_timer_t  t0 = cs_timer_time();

  /* Dirichlet values at boundary faces are first computed */
  cs_equation_compute_dirichlet_fb(mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_eval,
                                   cs_cdofb_cell_bld[0],
                                   field_val);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO face-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdofb_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_build_system(const cs_mesh_t            *mesh,
                             const cs_real_t            *field_val,
                             double                      dt_cur,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data,
                             cs_real_t                  *rhs,
                             cs_matrix_t                *matrix)
{
  CS_UNUSED(dt_cur);
  /* Sanity checks */
  assert(rhs != NULL && matrix != NULL);

  /* Test to remove */
  if (eqp->flag & CS_EQUATION_CONVECTION)
    bft_error(__FILE__, __LINE__, 0,
              _(" Convection term is not handled yet.\n"));

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;

  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)data;

  /* Dirichlet values at boundary faces are first computed */
  cs_real_t  *dir_values = NULL;
  BFT_MALLOC(dir_values, quant->n_b_faces, cs_real_t);
  memset(dir_values, 0, quant->n_b_faces*sizeof(cs_real_t));

  cs_equation_compute_dirichlet_fb(mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_cur + dt_cur,
                                   cs_cdofb_cell_bld[0],
                                   dir_values);

  /* Tag faces with a non-homogeneous Neumann BC */
  short int  *neu_tags = cs_equation_tag_neumann_face(quant, eqp);

  /* No source term related to face DoFs. No update of the RHS. */

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav,        \
         dir_values, neu_tags, field_val,                                \
         cs_cdofb_cell_sys, cs_cdofb_cell_bld)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdofb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];

    /* Set inside the OMP section so that each thread has its own value */

    /* Initialization of the values of properties */
    double  time_pty_val = 1.0;
    double  reac_pty_vals[CS_CDO_N_MAX_REACTIONS];

    const cs_real_t  t_eval_pty = t_cur + 0.5*dt_cur;

    cs_equation_init_properties(eqp, eqb, t_eval_pty,
                                &time_pty_val, reac_pty_vals, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag = cs_equation_cell_mesh_flag(cell_flag, eqb);

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _init_cell_system(cell_flag, cm, eqp, eqb, eqc,
                        dir_values, neu_tags, field_val, t_eval_pty, // in
                        csys, cb);                                   // out

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 2
      if (cs_dbg_cw_test(cm)) cs_cell_mesh_dump(cm);
#endif

      /* DIFFUSION TERM */
      /* ============== */

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        if (!(eqb->diff_pty_uniform))
          cs_equation_set_diffusion_property_cw(eqp, cm, t_eval_pty, cell_flag,
                                                cb);

        /* local matrix owned by the cellwise builder (store in cb->loc) */
        eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        /* Add the local diffusion operator to the local system */
        cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
        if (cs_dbg_cw_test(cm))
          cs_cell_sys_dump("\n>> Local system after diffusion", c_id, csys);
#endif
      } /* END OF DIFFUSION */

      /* SOURCE TERM */
      /* =========== */

      if (cs_equation_param_has_sourceterm(eqp)) {

        /* Reset the local contribution */
        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */
        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (const cs_xdef_t **)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        t_eval_pty,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        csys->rhs[cm->n_fc] += csys->source[cm->n_fc];

        /* Reset the value of the source term for the cell DoF
           Source term is only hold by the cell DoF in face-based schemes */
        eqc->source_terms[c_id] = csys->source[cm->n_fc];

      } /* End of term source contribution */

      /* UNSTEADY TERM + TIME SCHEME */
      /* =========================== */

      if (cs_equation_param_has_time(eqp)) {

        /* Get the value of the time property */
        double  tpty_val = 1/dt_cur;
        if (eqb->time_pty_uniform)
          tpty_val *= time_pty_val;
        else
          tpty_val *= cs_property_get_cell_value(c_id, t_eval_pty,
                                                 eqp->time_property);

        /* Assign local matrix to a mass matrix to define */
        cs_sdm_t  *mass_mat = cb->loc;
        assert(mass_mat->n_rows == mass_mat->n_cols);
        assert(mass_mat->n_rows == cm->n_fc + 1);

        if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) {

          /* Use a vector to deal with the diagonal matrix */
          memset(mass_mat->val, 0, sizeof(cs_real_t)*(cm->n_fc + 1));
          mass_mat->val[cm->n_fc] = cm->vol_c * tpty_val;

        }
        else {

          memset(mass_mat->val, 0,
                 sizeof(cs_real_t)*(cm->n_fc + 1)*(cm->n_fc + 1));

          bft_error(__FILE__, __LINE__, 0,
                    "%s: Not implemented yet.", __func__);

        }

        /* Apply the time discretization to the local system.
           Update csys (matrix and rhs) */
        eqc->apply_time_scheme(eqp, tpty_val, mass_mat, eqb->sys_flag, cb,
                               csys);

      } /* END OF TIME CONTRIBUTION */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> Local system matrix before condensation",
                         c_id, csys);
#endif

      /* Neumann boundary conditions */
      if ((cell_flag & CS_FLAG_BOUNDARY) && csys->has_nhmg_neumann) {
        for (short int f  = 0; f < cm->n_fc; f++)
          csys->rhs[f] += csys->neu_values[f];
      }

      /* Static condensation of the local system matrix of size n_fc + 1 into
         a matrix of size n_fc.
         Store data in rc_tilda and acf_tilda to compute the values at cell
         centers after solving the system */
      cs_static_condensation_scalar_eq(connect->c2f,
                                       eqc->rc_tilda,
                                       eqc->acf_tilda,
                                       cb, csys);

      /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ======================================================= */

      if (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_PENA) {

        /* Weakly enforced Dirichlet BCs for cells attached to the boundary
           csys is updated inside (matrix and rhs) */
        if (cell_flag & CS_FLAG_BOUNDARY)
          eqc->enforce_dirichlet(eqp->diffusion_hodge, cm,   // in
                                 eqc->boundary_flux_op,      // function
                                 fm, cb, csys);              // in/out

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", c_id, csys);
#endif

      /* Assemble the local system (related to vertices only since one applies
         a static condensation) to the global system */
      cs_equation_assemble_f(csys,
                             connect->range_sets[CS_CDO_CONNECT_FACE_SP0],
                             eqp,
                             1, /* n_face_dofs */
                             rhs, mav);

    } // Main loop on cells

  } // OPENMP Block

  cs_matrix_assembler_values_done(mav); // optional

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 2
  cs_dbg_darray_to_listing("FINAL RHS_FACE", quant->n_faces, rhs, 8);
#endif

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(neu_tags);
  cs_matrix_assembler_values_finalize(&mav);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdofb_scaleq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_update_field(const cs_real_t              *solu,
                             const cs_real_t              *rhs,
                             const cs_equation_param_t    *eqp,
                             cs_equation_builder_t        *eqb,
                             void                         *data,
                             cs_real_t                    *field_val)
{
  CS_UNUSED(rhs);
  CS_UNUSED(eqp);

  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)data;
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Set computed solution in builder->face_values */
  memcpy(eqc->face_values, solu, quant->n_faces*sizeof(cs_real_t));

  /* Compute values at cells pc from values at faces pf
     pc = acc^-1*(RHS - Acf*pf) */
  cs_static_condensation_recover_scalar(cs_shared_connect->c2f,
                                        eqc->rc_tilda,
                                        eqc->acf_tilda,
                                        solu,
                                        field_val);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain between time t_cur and t_cur + dt_cur
 *         Case of scalar-valued CDO face-based scheme
 *
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in, out] eqb             pointer to a cs_equation_builder_t structure
 * \param[in, out] context         pointer to a scheme builder structure
 * \param[in]      var_field_id    id of the variable field
 * \param[in]      bflux_field_id  id of the variable field
 * \param[in]      dt_cur          current value of the time step
 *
 * \return a pointer to a cs_equation_balance_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_balance_t *
cs_cdofb_scaleq_balance(const cs_equation_param_t     *eqp,
                        cs_equation_builder_t         *eqb,
                        void                          *context,
                        int                            var_field_id,
                        int                            bflux_field_id,
                        cs_real_t                      dt_cur)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  t_eval_pty = t_cur + 0.5*dt_cur;

  cs_timer_t  t0 = cs_timer_time();

  cs_field_t  *pot = cs_field_by_id(var_field_id);
  cs_field_t  *bflux = cs_field_by_id(bflux_field_id);
  cs_cdofb_scaleq_t  *eqc = (cs_cdofb_scaleq_t *)context;

  /* Allocate and initialize the structure storing the balance evaluation */
  cs_equation_balance_t  *eb = cs_equation_balance_create(cs_flag_primal_cell,
                                                          quant->n_cells);

  cs_equation_balance_reset(eb);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, pot, bflux,             \
         eb, cs_cdofb_cell_bld)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = cs_cdofb_cell_bld[t_id];

    /* Set inside the OMP section so that each thread has its own value */

    cs_real_t  _p_cur[10], _p_prev[10], _p_theta[10];
    cs_real_t  *p_cur = NULL, *p_prev = NULL, *p_theta = NULL;

    if (connect->n_max_vbyc > 10) {
      BFT_MALLOC(p_cur, connect->n_max_vbyc, cs_real_t);
      BFT_MALLOC(p_prev, connect->n_max_vbyc, cs_real_t);
      BFT_MALLOC(p_theta, connect->n_max_vbyc, cs_real_t);
    }
    else {
      p_cur = _p_cur;
      p_prev = _p_prev;
      p_theta = _p_theta;
    }

    /* Initialization of the values of properties */
    double  time_pty_val = 1.0;
    double  reac_pty_vals[CS_CDO_N_MAX_REACTIONS];

    cs_equation_init_properties(eqp, eqb, t_eval_pty,
                                &time_pty_val, reac_pty_vals, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag = cs_equation_cell_mesh_flag(cell_flag, eqb);

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Set the value of the current potential */
      for (short int f = 0; f < cm->n_fc; f++)
        p_cur[f] = eqc->face_values[cm->f_ids[f]];
      p_cur[cm->n_fc] = pot->val[cm->c_id];

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_SCALEQ_DBG > 3
      if (cs_dbg_cw_test(cm)) cs_cell_mesh_dump(cm);
#endif

      /* Set p_theta */
      switch (eqp->time_scheme) {
      case CS_TIME_SCHEME_EXPLICIT:
        for (short int i = 0; i < cm->n_fc + 1; i++)
          p_theta[i] = p_prev[i];
        break;

      case CS_TIME_SCHEME_CRANKNICO:
        for (short int i = 0; i < cm->n_fc + 1; i++)
          p_theta[i] = 0.5*(p_cur[i] + p_prev[i]);
        break;

      case CS_TIME_SCHEME_THETA:
        for (short int i = 0; i < cm->n_fc + 1; i++)
          p_theta[i] = eqp->theta*p_cur[i] + (1-eqp->theta)*p_prev[i];
        break;

      default:
        for (short int i = 0; i < cm->n_fc + 1; i++)
          p_theta[i] = p_cur[i];
        break;

      } /* Switch on time scheme */

      /* DIFFUSION TERM */
      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        if (!(eqb->diff_pty_uniform))
          cs_equation_set_diffusion_property_cw(eqp, cm, t_eval_pty, cell_flag,
                                                cb);

        /* local matrix owned by the cellwise builder (store in cb->loc) */
        eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        cs_real_t  *res = cb->values;
        for (short int v = 0; v < cm->n_vc; v++)
          res[v] = 0.;

        cs_sdm_square_matvec(cb->loc, p_theta, res);

        eb->diffusion_term[cm->c_id] = res[cm->n_fc];

      } /* END OF DIFFUSION */

      /* SOURCE TERM */
      /* =========== */

      if (cs_equation_param_has_sourceterm(eqp)) {

        /* Reset the local contribution */
        cs_real_t  *src = cb->values;
        memset(src, 0, (cm->n_fc + 1)*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */
        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (const cs_xdef_t **)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        t_eval_pty,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        src);

        eb->source_term[cm->c_id] += src[cm->n_fc];

      } /* End of term source contribution */

      /* /\* UNSTEADY TERM + TIME SCHEME *\/ */
      /* /\* =========================== *\/ */

      /* if (cs_equation_param_has_time(eqp)) { */

      /*   /\* Get the value of the time property *\/ */
      /*   double  tpty_val = 1/dt_cur; */
      /*   if (eqb->time_pty_uniform) */
      /*     tpty_val *= time_pty_val; */
      /*   else */
      /*     tpty_val *= cs_property_get_cell_value(c_id, t_eval_pty, */
      /*                                            eqp->time_property); */

      /*   /\* Assign local matrix to a mass matrix to define *\/ */
      /*   cs_sdm_t  *mass_mat = cb->loc; */
      /*   assert(mass_mat->n_rows == mass_mat->n_cols); */
      /*   assert(mass_mat->n_rows == cm->n_fc + 1); */

      /*   if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { */

      /*     /\* Use a vector to deal with the diagonal matrix *\/ */
      /*     memset(mass_mat->val, 0, sizeof(cs_real_t)*(cm->n_fc + 1)); */
      /*     mass_mat->val[cm->n_fc] = cm->vol_c * tpty_val; */

      /*   } */
      /*   else { */

      /*     memset(mass_mat->val, 0, */
      /*            sizeof(cs_real_t)*(cm->n_fc + 1)*(cm->n_fc + 1)); */

      /*     bft_error(__FILE__, __LINE__, 0, */
      /*               "%s: Not implemented yet.", __func__); */

      /*   } */

      /*   /\* Apply the time discretization to the local system. */
      /*      Update csys (matrix and rhs) *\/ */
      /*   eqc->apply_time_scheme(eqp, tpty_val, mass_mat, eqb->sys_flag, cb, */
      /*                          csys); */

      /* } /\* END OF TIME CONTRIBUTION *\/ */


    } // Main loop on cells

    if (p_cur != _p_cur) {
      BFT_FREE(p_cur);
      BFT_FREE(p_prev);
      BFT_FREE(p_theta);
    }

  } // OPENMP Block

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
    eb->balance[c_id] =
      eb->unsteady_term[c_id]  + eb->reaction_term[c_id]  +
      eb->diffusion_term[c_id] + eb->advection_term[c_id] +
      eb->source_term[c_id];

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  return eb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_cdofb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_scaleq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data)
{
  CS_UNUSED(eqname); // avoid a compilation warning
  CS_UNUSED(eqp);

  char *postlabel = NULL;
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_i_faces = connect->n_faces[2];
  const cs_real_t  *face_pdi = cs_cdofb_scaleq_get_face_values(data);

  /* Field post-processing */
  int  len = strlen(field->name) + 8 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.Border", field->name);

  cs_post_write_var(CS_POST_MESH_BOUNDARY,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    postlabel,
                    field->dim,
                    true,
                    true,                  // true = original mesh
                    CS_POST_TYPE_cs_real_t,
                    NULL,                  // values on cells
                    NULL,                  // values at internal faces
                    face_pdi + n_i_faces,  // values at border faces
                    cs_shared_time_step);  // time step management structure


  BFT_FREE(postlabel);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at each face
 *
 * \param[in] data       pointer to cs_cdofb_scaleq_t structure
 *
 * \return  a pointer to an array of double (size n_faces)
 */
/*----------------------------------------------------------------------------*/

double *
cs_cdofb_scaleq_get_face_values(const void    *data)
{
  const cs_cdofb_scaleq_t  *eqc = (const cs_cdofb_scaleq_t *)data;

  if (eqc == NULL)
    return NULL;
  else
    return eqc->face_values;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
