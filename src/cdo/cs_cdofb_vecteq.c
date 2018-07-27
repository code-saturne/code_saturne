/*============================================================================
 * Build an algebraic CDO face-based system for unsteady convection/diffusion
 * reaction of vector-valued equations with source terms
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

#include "cs_cdofb_vecteq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_VECTEQ_DBG      0

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
 * \brief  Initialize the local builder structure used for building the system
 *         cellwise
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

  BFT_MALLOC(cb->ids, n_fc + 1, short int);
  memset(cb->ids, 0, (n_fc + 1)*sizeof(short int));

  int  size = n_fc*(n_fc+1);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(cs_real_t));

  size = 2*n_fc;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  short int  *block_sizes = cb->ids;
  for (int i = 0; i < n_fc + 1; i++)
    block_sizes[i] = 3;

  /* Local square dense matrices used during the construction of
     operators */
  cb->hdg = cs_sdm_square_create(n_fc);
  cb->loc = cs_sdm_block_create(n_fc + 1, n_fc + 1, block_sizes, block_sizes);
  cb->aux = cs_sdm_square_create(n_fc + 1);

  return cb;
}

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      cell_flag   flag related to the current cell
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         pointer to a cs_cdofb_vecteq_t structure
 * \param[in]      dir_values  Dirichlet values associated to each face
 * \param[in]      neu_tags    definition id related to each Neumann face
 * \param[in]      field_tn    values of the field at the last computed time
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_init_cell_system(const cs_flag_t               cell_flag,
                                 const cs_cell_mesh_t         *cm,
                                 const cs_equation_param_t    *eqp,
                                 const cs_equation_builder_t  *eqb,
                                 const cs_cdofb_vecteq_t      *eqc,
                                 const cs_real_t               dir_values[],
                                 const short int               neu_tags[],
                                 const cs_real_t               field_tn[],
                                 cs_real_t                     t_eval,
                                 cs_cell_sys_t                *csys,
                                 cs_cell_builder_t            *cb)
{
  CS_UNUSED(cb);

  /* Cell-wise view of the linear system to build */
  const int  n_blocks = cm->n_fc + 1;
  const int  n_dofs = 3*n_blocks;

  short int  *block_sizes = cb->ids;
  for (int i = 0; i < n_blocks; i++)
    block_sizes[i] = 3;

  /* Initialize the local system */
  cs_cell_sys_reset(cell_flag, n_dofs, cm->n_fc, csys);

  csys->cell_flag = cell_flag;
  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;
  csys->face_shift = cs_shared_connect->n_faces[2]; /* shift = n_i_faces */

  cs_sdm_block_init(csys->mat, n_blocks, n_blocks, block_sizes, block_sizes);

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_lnum_t  f_id = cm->f_ids[f];
    for (int k = 0; k < 3; k++) {
      csys->dof_ids[3*f + k] = 3*f_id + k;
      csys->val_n[3*f + k] = eqc->face_values[3*f_id + k];
    }

  }

  for (int k = 0; k < 3; k++) {

    const cs_lnum_t  dof_id = 3*cm->c_id+k;
    const cs_lnum_t  _shift = 3*cm->n_fc + k;

    csys->dof_ids[_shift] = dof_id;
    csys->val_n[_shift] = field_tn[dof_id];

  }

  /* Update rhs with the previous computation of source term if needed */
  if (cs_equation_param_has_sourceterm(eqp)) {
    if (cs_equation_param_has_time(eqp)) {

      /* Source terms attached to cells: Need to update rhs because the part
         related to cell is used in the static condensation */
      cs_cdo_time_update_rhs(eqp,
                             3, /* stride */
                             1, /* n_dofs */
                             csys->dof_ids + cm->n_fc,
                             eqc->source_terms,
                             csys->rhs + 3*cm->n_fc);

    }
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY) {

    /* Set the generic part */
    cs_equation_init_cell_sys_bc(eqb, cm, csys);

    cs_equation_fb_set_cell_bc(cm,
                               cs_shared_connect,
                               cs_shared_quant,
                               eqp,
                               dir_values,
                               neu_tags,
                               t_eval,
                               csys,
                               cb);

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    cs_dbg_check_hmg_dirichlet_cw(__func__, csys);
#endif
  } /* Border cell */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-Fb scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdofb_vecteq_is_initialized(void)
{
  if (cs_cdofb_cell_sys == NULL || cs_cdofb_cell_bld == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to CDO
 *         vector-valued face-based schemes.
 *         Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_init_common(const cs_cdo_quantities_t     *quant,
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

  const short int  n_blocks = connect->n_max_fbyc + 1;
  const short int  n_max_dofs = 3*n_blocks;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_builder_t  *cb = _cell_builder_create(connect);
    short int  *block_sizes = cb->ids;
    for (int i = 0; i < n_blocks; i++)
      block_sizes[i] = 3;

    cs_cdofb_cell_sys[t_id] = cs_cell_sys_create(n_max_dofs,
                                                 n_blocks - 1,
                                                 n_blocks,
                                                 block_sizes);
    cs_cdofb_cell_bld[t_id] = cb;
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_builder_t  *cb = _cell_builder_create(connect);
  short int  *block_sizes = cb->ids;
  for (int i = 0; i < n_blocks; i++)
    block_sizes[i] = 3;

  cs_cdofb_cell_sys[0] =  cs_cell_sys_create(n_max_dofs,
                                             n_blocks - 1,
                                             n_blocks,
                                             block_sizes);
  cs_cdofb_cell_bld[0] = cb;
#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the pointer to the related cs_matrix_structure_t
 *
 * \return a  pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

const cs_matrix_structure_t *
cs_cdofb_vecteq_matrix_structure(void)
{
  return cs_shared_ms;
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
cs_cdofb_vecteq_get(cs_cell_sys_t       **csys,
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
cs_cdofb_vecteq_finalize_common(void)
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
  cs_cdofb_cell_bld = NULL;
  cs_cdofb_cell_sys = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_vecteq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_vecteq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOFB || eqp->dim != 3)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid type of equation.\n"
              " Expected: vector-valued CDO face-based equation.", __func__);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[0];

  cs_cdofb_vecteq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdofb_vecteq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  /* Dimensions of the algebraic system */
  eqc->n_dofs = 3*(n_faces + n_cells);

  eqb->sys_flag = CS_FLAG_SYS_VECTOR;
  eqb->msh_flag = CS_CDO_LOCAL_PF | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_CDO_LOCAL_EV | CS_CDO_LOCAL_EF | CS_CDO_LOCAL_EFQ;

  BFT_MALLOC(eqc->face_values, 3*n_faces, cs_real_t);
  BFT_MALLOC(eqc->rc_tilda, 3*n_cells, cs_real_t);
# pragma omp parallel if (3*n_cells > CS_THR_MIN)
  {
    /* Values at each face (interior and border) i.e. take into account BCs */
#   pragma omp for nowait
    for (cs_lnum_t i = 0; i < 3*n_faces; i++) eqc->face_values[i] = 0;

    /* Store the last computed values of the field at cell centers and the data
       needed to compute the cell values from the face values.
       No need to synchronize all these quantities since they are only cellwise
       quantities. */
#   pragma omp for
    for (cs_lnum_t i = 0; i < 3*n_cells; i++) eqc->rc_tilda[i] = 0;
  }

  /* Assume the 3x3 matrix is diagonal */
  BFT_MALLOC(eqc->acf_tilda, 3*connect->c2f->idx[n_cells], cs_real_t);
  memset(eqc->acf_tilda, 0, 3*connect->c2f->idx[n_cells]*sizeof(cs_real_t));

  /* Diffusion part */
  /* -------------- */

  eqc->get_stiffness_matrix = NULL;
  eqc->boundary_flux_op = NULL;

  if (cs_equation_param_has_diffusion(eqp)) {

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      eqc->get_stiffness_matrix = cs_hodge_fb_cost_get_stiffness;
      eqc->boundary_flux_op = NULL; //cs_cdofb_diffusion_cost_flux_op;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqc->get_stiffness_matrix = cs_hodge_fb_voro_get_stiffness;
      eqc->boundary_flux_op = NULL; //cs_cdofb_diffusion_cost_flux_op;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to build the diffusion term.",
                __func__);

    } /* Switch on Hodge algo. */

  } /* Diffusion part */

  eqc->enforce_dirichlet = NULL;
  switch (eqp->enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_block_dirichlet;
    break;
  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_block_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

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

    BFT_MALLOC(eqc->source_terms, 3*n_cells, cs_real_t);
#   pragma omp parallel for if (3*n_cells > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < 3*n_cells; i++) eqc->source_terms[i] = 0;

  } /* There is at least one source term */

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdofb_vecteq_t structure
 *
 * \param[in, out]  data   pointer to a cs_cdofb_vecteq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_vecteq_free_context(void   *data)
{
  cs_cdofb_vecteq_t   *eqc  = (cs_cdofb_vecteq_t *)data;

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
 * \param[in, out] data           pointer to cs_cdofb_vecteq_t structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_initialize_system(const cs_equation_param_t  *eqp,
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
  cs_lnum_t  size = 3*quant->n_faces;
  BFT_MALLOC(*system_rhs, size, cs_real_t);
# pragma omp parallel for if  (size > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < size; i++) (*system_rhs)[i] = 0.0;

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings when the fields
 *         stem from a vector CDO face-based scheme.
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in, out] eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval      time at which one evaluates BCs
 * \param[in, out] field_val   pointer to the values of the variable field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_set_dir_bc(const cs_mesh_t              *mesh,
                           const cs_equation_param_t    *eqp,
                           cs_equation_builder_t        *eqb,
                           cs_real_t                     t_eval,
                           cs_real_t                     field_val[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_timer_t  t0 = cs_timer_time();

  /* Compute the values of the Dirichlet BC */
  cs_equation_compute_dirichlet_fb(mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_eval,
                                   cs_cdofb_cell_bld[0], /* static variable */
                                   field_val);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a vector convection/diffusion
 *         equation with a CDO face-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdofb_vecteq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_build_system(const cs_mesh_t            *mesh,
                             const cs_real_t            *field_val,
                             double                      dt_cur,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data,
                             cs_real_t                  *rhs,
                             cs_matrix_t                *matrix)
{
  /* Sanity checks */
  assert(rhs != NULL && matrix != NULL);

  /* Test to remove */
  if (eqp->flag & CS_EQUATION_CONVECTION)
    bft_error(__FILE__, __LINE__, 0,
              _(" Convection term is not handled yet.\n"));
  if (eqp->flag & CS_EQUATION_UNSTEADY)
    bft_error(__FILE__, __LINE__, 0,
              _(" Unsteady terms are not handled yet.\n"));

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;

  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)data;

  /* Dirichlet values at boundary faces are first computed */
  cs_real_t  *dir_values = NULL;
  BFT_MALLOC(dir_values, 3*quant->n_b_faces, cs_real_t);
  memset(dir_values, 0, 3*quant->n_b_faces*sizeof(cs_real_t));

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

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)   \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav,      \
         dir_values, neu_tags, field_val,                              \
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
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, eqp, eqb, eqc,
                                       dir_values, neu_tags,
                                       field_val, t_eval_pty,
                                       csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 2
      if (cs_dbg_cw_test(cm)) cs_cell_mesh_dump(cm);
#endif

      /* DIFFUSION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ============================================== */

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        if (!(eqb->diff_pty_uniform))
          cs_equation_set_diffusion_property_cw(eqp, cm, t_eval_pty, cell_flag,
                                                cb);

        /* local matrix owned by the cellwise builder (store in cb->loc) */
        eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        if (eqp->diffusion_hodge.is_iso == false)
          bft_error(__FILE__, __LINE__, 0, " %s: Case not handle yet\n",
                    __func__);

        /* Add the local diffusion operator to the local system */
        const cs_real_t  *sval = cb->loc->val;
        for (int bi = 0; bi < cm->n_fc + 1; bi++) {
          for (int bj = 0; bj < cm->n_fc + 1; bj++) {

            /* Retrieve the 3x3 matrix */
            cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bi, bj);
            assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

            const cs_real_t  _val = sval[(cm->n_fc+1)*bi+bj];
            bij->val[0] += _val;
            bij->val[4] += _val;
            bij->val[8] += _val;

          }
        }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
        if (cs_dbg_cw_test(cm))
          cs_cell_sys_dump("\n>> Local system after diffusion", c_id, csys);
#endif
      } /* END OF DIFFUSION */

      /* SOURCE TERM COMPUTATION */
      /* ======================= */

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

        for (int k = 0; k < 3; k++)
          csys->rhs[3*cm->n_fc + k] += csys->source[3*cm->n_fc + k];

        /* Reset the value of the source term for the cell DoF
           Source term is only hold by the cell DoF in face-based schemes */
        for (int k = 0; k < 3; k++)
          eqc->source_terms[3*c_id + k] = csys->source[3*cm->n_fc + k];

      } /* End of term source contribution */

      /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
       * Operations that have to be performed BEFORE the static condensation
       */
      if (cell_flag & CS_FLAG_BOUNDARY) {

        /* Neumann boundary conditions */
        if (csys->has_nhmg_neumann) {
          for (short int f = 0; f < 3*cm->n_fc; f++)
            csys->rhs[f] += csys->neu_values[f];
        }

      } /* Boundary cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> Local system matrix before condensation",
                         c_id, csys);
#endif

      /* Static condensation of the local system stored inside a block matrix of
         size n_fc + 1 into a block matrix of size n_fc.
         Store information in the context structure in order to be able to
         compute the values at cell centers. */

      cs_static_condensation_vector_eq(connect->c2f,
                                       eqc->rc_tilda, eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> Local system matrix after condensation",
                         c_id, csys);
#endif

      /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
       * Operations that have to be performed AFTER the static condensation
       */
      if (cell_flag & CS_FLAG_BOUNDARY) {

        if (eqp->enforcement == CS_PARAM_BC_ENFORCE_PENALIZED) {

          /* Weakly enforced Dirichlet BCs for cells attached to the boundary
             csys is updated inside (matrix and rhs) */
          eqc->enforce_dirichlet(eqp->diffusion_hodge, cm,   // in
                                 eqc->boundary_flux_op,      // function
                                 fm, cb, csys);              // in/out

        }

      } /* Boundary cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 0
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", c_id, csys);
#endif

      /* ASSEMBLY */
      /* ======== */

      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];

      /* Matrix assembly */
      cs_equation_assemble_block_matrix(csys, rs, 3, mav);

      /* Assemble RHS */
      for (short int i = 0; i < 3*cm->n_fc; i++) {
#       pragma omp atomic
        rhs[csys->dof_ids[i]] += csys->rhs[i];
      }

    } /* Main loop on cells */

  } /* OpenMP Block */

  cs_matrix_assembler_values_done(mav);    /* optional */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 2
  cs_dbg_darray_to_listing("FINAL RHS_FACE", quant->n_faces, rhs, 9);
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
 * \param[in, out] data       pointer to cs_cdofb_vecteq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_update_field(const cs_real_t              *solu,
                             const cs_real_t              *rhs,
                             const cs_equation_param_t    *eqp,
                             cs_equation_builder_t        *eqb,
                             void                         *data,
                             cs_real_t                    *field_val)
{
  CS_UNUSED(rhs);
  CS_UNUSED(eqp);

  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)data;
  cs_timer_t  t0 = cs_timer_time();

  /* Set computed solution in builder->face_values */
  memcpy(eqc->face_values, solu, 3*cs_shared_quant->n_faces*sizeof(cs_real_t));

  /* Build the field inside each cell */
  cs_static_condensation_recover_vector(cs_shared_connect->c2f,
                                        eqc->rc_tilda,
                                        eqc->acf_tilda,
                                        eqc->face_values,
                                        field_val);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_cdofb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data)
{
  CS_UNUSED(eqname);    /* avoid a compilation warning */
  CS_UNUSED(eqp);

  char *postlabel = NULL;
  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_i_faces = connect->n_faces[2];
  const cs_real_t  *face_pdi = cs_cdofb_vecteq_get_face_values(data);

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
 * \brief  Get the computed values at mesh cells from the inverse operation
 *         w.r.t. the static condensation (DoF used in the linear system are
 *         located at primal faces)
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdofb_vecteq_get_cell_values(void      *context)
{
  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  return  pot->val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh faces for the current context.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of cs_real_t (size n_faces)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdofb_vecteq_get_face_values(void    *context)
{
  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t *)context;

  if (eqc == NULL)
    return NULL;
  else
    return  eqc->face_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read additional arrays (not defined as fields) but useful for the
 *         checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in]       scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_read_restart(cs_restart_t    *restart,
                             const char      *eqname,
                             void            *scheme_context)
{
  /* Only the face values are handled. Cell values are stored in a cs_field_t
     structure and thus are handled automatically. */
  if (restart == NULL)
    return;
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is NULL", __func__);
  if (scheme_context == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Scheme context is NULL", __func__);

  int retcode = CS_RESTART_SUCCESS;
  cs_cdofb_vecteq_t  *eqc = (cs_cdofb_vecteq_t  *)scheme_context;

  char sec_name[128];

  /* Handle interior faces */
  /* ===================== */

  const int  i_ml_id = cs_mesh_location_get_id_by_name(N_("interior_faces"));

  /* Define the section name */
  snprintf(sec_name, 127, "%s::i_face_vals", eqname);

  /* Check section */
  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     i_ml_id,
                                     3, /* vector-valued */
                                     CS_TYPE_cs_real_t);

  /* Read section */
  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      i_ml_id,
                                      3, /* vector-valued */
                                      CS_TYPE_cs_real_t,
                                      eqc->face_values);

  /* Handle boundary faces */
  /* ===================== */

  const int  b_ml_id = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  cs_real_t  *b_values = eqc->face_values + 3*cs_shared_quant->n_i_faces;

  /* Define the section name */
  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Check section */
  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     b_ml_id,
                                     1, /* vector-valued */
                                     CS_TYPE_cs_real_t);

  /* Read section */
  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      b_ml_id,
                                      1, /* vector-valued */
                                      CS_TYPE_cs_real_t,
                                      b_values);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Write additional arrays (not defined as fields) but useful for the
 *         checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in]       scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_vecteq_write_restart(cs_restart_t    *restart,
                              const char      *eqname,
                              void            *scheme_context)
{
  /* Only the face values are handled. Cell values are stored in a cs_field_t
     structure and thus are handled automatically. */
  if (restart == NULL)
    return;
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is NULL", __func__);

  const cs_cdofb_vecteq_t  *eqc = (const cs_cdofb_vecteq_t  *)scheme_context;

  char sec_name[128];

  /* Handle interior faces */
  /* ===================== */

  const int  i_ml_id = cs_mesh_location_get_id_by_name(N_("interior_faces"));

  /* Define the section name */
  snprintf(sec_name, 127, "%s::i_face_vals", eqname);

  /* Write interior face section */
  cs_restart_write_section(restart,
                           sec_name,
                           i_ml_id,
                           3,   /* vector-valued */
                           CS_TYPE_cs_real_t,
                           eqc->face_values);

  /* Handle boundary faces */
  /* ===================== */

  const int  b_ml_id = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  const cs_real_t  *b_values = eqc->face_values + 3*cs_shared_quant->n_i_faces;

  /* Define the section name */
  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Write boundary face section */
  cs_restart_write_section(restart,
                           sec_name,
                           b_ml_id,
                           3,   /* vector-valued */
                           CS_TYPE_cs_real_t,
                           b_values);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
