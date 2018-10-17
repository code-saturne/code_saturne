/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
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
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_local.h"
#include "cs_cdo_time.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_param.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_scheme_geometry.h"
#include "cs_search.h"
#include "cs_source_term.h"
#include "cs_timer.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovb_vecteq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
 *  \file cs_cdovb_vecteq.c
 *
 * \brief Build an algebraic CDO vertex-based system for unsteady
 *        convection-diffusion-reaction of vector-valued equations with
 *        source terms
 *
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVB_VECTEQ_DBG     0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Structure to enable a full cellwise strategy during the system building */
static cs_cell_sys_t      **cs_cdovb_cell_sys = NULL;
static cs_cell_builder_t  **cs_cdovb_cell_bld = NULL;

/* Pointer to shared structures */
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;
static const cs_matrix_structure_t  *cs_shared_ms;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
  const int  n_vc = connect->n_max_vbyc;
  const int  n_ec = connect->n_max_ebyc;

  cs_cell_builder_t  *cb = cs_cell_builder_create();

  BFT_MALLOC(cb->ids, n_vc, short int);
  memset(cb->ids, 0, n_vc*sizeof(short int));

  int  size = n_ec*(n_ec+1);
  size = CS_MAX(4*n_ec + 3*n_vc, size);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(cs_real_t));

  size = 2*n_ec;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of
     operators */
  cb->hdg = cs_sdm_square_create(n_ec);
  cb->aux = cs_sdm_square_create(n_vc);

  short int  *block_sizes = cb->ids;
  for (int i = 0; i < n_vc; i++)
    block_sizes[i] = 3;
  cb->loc = cs_sdm_block_create(n_vc, n_vc, block_sizes, block_sizes);

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
 * \param[in]      dir_values  Dirichlet values associated to each vertex
 * \param[in]      neu_tags    definition id related to each Neumann face
 * \param[in]      field_tn    values of the field at the last computed time
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_init_vb_cell_system(const cs_flag_t               cell_flag,
                     const cs_cell_mesh_t         *cm,
                     const cs_equation_param_t    *eqp,
                     const cs_equation_builder_t  *eqb,
                     const cs_real_t               dir_values[],
                     const short int               neu_tags[],
                     const cs_real_t               field_tn[],
                     cs_real_t                     t_eval,
                     cs_cell_sys_t                *csys,
                     cs_cell_builder_t            *cb)
{
  const int  n_blocks = cm->n_vc;
  const int  n_dofs = 3*n_blocks;

  short int  *block_sizes = cb->ids;
  for (int i = 0; i < n_blocks; i++)
    block_sizes[i] = 3;

  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;
  csys->cell_flag = cell_flag;

  /* Cell-wise view of the linear system to build:
     Initialize the local system */
  cs_cell_sys_reset(cm->n_fc, csys);

  cs_sdm_block_init(csys->mat, n_blocks, n_blocks, block_sizes, block_sizes);

  for (short int v = 0; v < cm->n_vc; v++) {
    const cs_lnum_t  v_id = cm->v_ids[v];
    for (int k = 0; k < 3; k++) {
      csys->dof_ids[3*v + k] = 3*v_id + k;
      csys->val_n[3*v + k] = field_tn[3*v_id + k];
    }
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY) {

    /* Set the generic part */
    cs_equation_init_cell_sys_bc(eqb, cm, csys);

    /* Set the bc (specific part) */
    cs_equation_vb_set_cell_bc(cm,
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

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-Vb scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdovb_vecteq_is_initialized(void)
{
  if (cs_cdovb_cell_sys == NULL || cs_cdovb_cell_bld == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Allocate work buffer and general structures related to CDO
 *           vertex-based schemes
 *           Set shared pointers.
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_init_common(const cs_cdo_quantities_t    *quant,
                            const cs_cdo_connect_t       *connect,
                            const cs_time_step_t         *time_step,
                            const cs_matrix_structure_t  *ms)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
  cs_shared_ms = ms;

  /* Structure used to build the final system by a cell-wise process */
  assert(cs_glob_n_threads > 0);  /* Sanity check */

  BFT_MALLOC(cs_cdovb_cell_sys, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(cs_cdovb_cell_bld, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_cdovb_cell_sys[i] = NULL;
    cs_cdovb_cell_bld[i] = NULL;
  }

  const short int  n_blocks = connect->n_max_vbyc;
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

    cs_cdovb_cell_sys[t_id] = cs_cell_sys_create(n_max_dofs,
                                                 connect->n_max_fbyc,
                                                 n_blocks,
                                                 block_sizes);
    cs_cdovb_cell_bld[t_id] = cb;
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_builder_t  *cb = _cell_builder_create(connect);
  short int  *block_sizes = cb->ids;
  for (int i = 0; i < n_blocks; i++)
    block_sizes[i] = 3;

  cs_cdovb_cell_sys[0] =  cs_cell_sys_create(n_max_dofs,
                                             connect->n_max_fbyc,
                                             n_blocks,
                                             block_sizes);
  cs_cdovb_cell_bld[0] = cb;
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
cs_cdovb_vecteq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  *csys = cs_cdovb_cell_sys[t_id];
  *cb = cs_cdovb_cell_bld[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO vertex-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_finalize_common(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(cs_cdovb_cell_sys[t_id]));
    cs_cell_builder_free(&(cs_cdovb_cell_bld[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(cs_cdovb_cell_sys[0]));
  cs_cell_builder_free(&(cs_cdovb_cell_bld[0]));
#endif /* openMP */

  BFT_FREE(cs_cdovb_cell_sys);
  BFT_FREE(cs_cdovb_cell_bld);
  cs_cdovb_cell_bld = NULL;
  cs_cdovb_cell_sys = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_vecteq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_vecteq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVB || eqp->dim != 3)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of equation.\n"
              " Expected: vector-valued CDO vertex-based equation.", __func__);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->n_vertices;

  cs_cdovb_vecteq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdovb_vecteq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  eqc->n_dofs = 3*n_vertices;

  eqb->sys_flag = CS_FLAG_SYS_VECTOR;

  /* Flag to indicate the minimal set of quantities to build in a cell mesh
     According to the situation, additional flags have to be set */
  eqb->msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PE |
    CS_CDO_LOCAL_EV;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
    CS_CDO_LOCAL_FEQ;

  /* Diffusion part */
  /* -------------- */

  eqc->get_stiffness_matrix = NULL;
  eqc->bdy_flux_op = NULL;
  if (cs_equation_param_has_diffusion(eqp)) {

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_cost_get_stiffness;
      eqc->bdy_flux_op = cs_cdovb_diffusion_cost_flux_op;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_voro_get_stiffness;
      eqc->bdy_flux_op = cs_cdovb_diffusion_cost_flux_op;
      break;

    case CS_PARAM_HODGE_ALGO_WBS:
      eqb->msh_flag |= CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_PEQ |
        CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_wbs_get_stiffness;
      eqc->bdy_flux_op = cs_cdovb_diffusion_wbs_flux_op;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to build the diffusion term.",
                __func__);

    } /* Switch on Hodge algo. */

  } /* Has diffusion */

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
  /* -------------- */

  eqc->get_advection_matrix = NULL;
  eqc->add_advection_bc = NULL;

  /* Reaction part */
  /* ------------- */

  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->reaction_hodge.algo == CS_PARAM_HODGE_ALGO_WBS) {
      eqb->msh_flag |= CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FEQ |
        CS_CDO_LOCAL_HFQ;
      eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid choice of algorithm for the reaction term.",
                __func__);

  } /* Reaction */

  /* Time part */
  /* --------- */

  eqc->apply_time_scheme = NULL;
  if (cs_equation_param_has_time(eqp)) {

    if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_VORONOI) {
      eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
    }
    else if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_WBS) {
      if (eqp->do_lumping)
        eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
      else {
        eqb->msh_flag |= CS_CDO_LOCAL_PVQ|CS_CDO_LOCAL_DEQ|CS_CDO_LOCAL_PFQ |
          CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
      }
    }

    eqc->apply_time_scheme = cs_cdo_time_get_scheme_function(eqb->sys_flag,
                                                             eqp);

  } /* Time part */

  /* Source term part */
  /* ---------------- */

  eqc->source_terms = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, eqc->n_dofs, cs_real_t);
#   pragma omp parallel for if (eqc->n_dofs > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < eqc->n_dofs; i++)
      eqc->source_terms[i] = 0;

  } /* There is at least one source term */

  /* Pre-defined a cs_hodge_builder_t structure */
  eqc->hdg_mass.is_unity = true;
  eqc->hdg_mass.is_iso   = true;
  eqc->hdg_mass.inv_pty  = false;
  eqc->hdg_mass.type = CS_PARAM_HODGE_TYPE_VPCD;
  eqc->hdg_mass.algo = CS_PARAM_HODGE_ALGO_WBS;
  eqc->hdg_mass.coef = 1.0; /* not useful in this case */

  eqc->get_mass_matrix = cs_hodge_vpcd_wbs_get;

  /* Array used for extra-operations */
  eqc->cell_values = NULL;

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovb_vecteq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovb_vecteq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_vecteq_free_context(void   *builder)
{
  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)builder;

  if (eqc == NULL)
    return eqc;

  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->cell_values);

  /* Last free */
  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         builder structure
 *
 * \param[in]      eqp            pointer to a cs_equation_param_t structure
 * \param[in, out] eqb            pointer to a cs_equation_builder_t structure
 * \param[in, out] data           pointer to cs_cdovb_vecteq_t structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_initialize_system(const cs_equation_param_t  *eqp,
                                  cs_equation_builder_t      *eqb,
                                  void                       *data,
                                  cs_matrix_t               **system_matrix,
                                  cs_real_t                 **system_rhs)
{
  CS_UNUSED(eqp);

  if (data == NULL)
    return;
  assert(*system_matrix == NULL && *system_rhs == NULL);

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)data;
  cs_timer_t  t0 = cs_timer_time();

  /* Create the matrix related to the current algebraic system */
  *system_matrix = cs_matrix_create(cs_shared_ms);

  /* Allocate and initialize the related right-hand side */
  BFT_MALLOC(*system_rhs, eqc->n_dofs, cs_real_t);
#pragma omp parallel for if  (eqc->n_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < eqc->n_dofs; i++) (*system_rhs)[i] = 0.0;

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings when the fields
 *         stem from a vector CDO vertex-based scheme.
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in, out] eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval      time at which one evaluates BCs
 * \param[in, out] field_val   pointer to the values of the variable field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_set_dir_bc(const cs_mesh_t              *mesh,
                           const cs_equation_param_t    *eqp,
                           cs_equation_builder_t        *eqb,
                           cs_real_t                     t_eval,
                           cs_real_t                     field_val[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_timer_t  t0 = cs_timer_time();

  /* Compute the values of the Dirichlet BC */
  cs_equation_compute_dirichlet_vb(mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_eval,
                                   cs_cdovb_cell_bld[0], /* static variable */
                                   field_val);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a vector convection/diffusion
 *         equation with a CDO vertex-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the vertex field
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdovcb_vecteq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_build_system(const cs_mesh_t            *mesh,
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

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;

  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)data;

  /* Compute the values of the Dirichlet BC */
  cs_real_t  *dir_values = NULL;
  BFT_MALLOC(dir_values, 3*quant->n_vertices, cs_real_t);
  memset(dir_values, 0, 3*quant->n_vertices*sizeof(cs_real_t));

  cs_cdovb_vecteq_set_dir_bc(mesh, eqp, eqb, t_cur + dt_cur, dir_values);

  /* Tag faces with a non-homogeneous Neumann BC */
  short int  *neu_tags = cs_equation_tag_neumann_face(quant, eqp);

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav,       \
         dir_values, neu_tags, field_val, cs_cdovb_cell_sys, cs_cdovb_cell_bld)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_real_t  time_eval = t_cur + 0.5*dt_cur;

    /* Set inside the OMP section so that each thread has its own value
     * Each thread get back its related structures:
     * Get the cell-wise view of the mesh and the algebraic system */
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdovb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdovb_cell_bld[t_id];

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

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
      _init_vb_cell_system(cell_flag, cm, eqp, eqb,
                           dir_values, neu_tags, field_val, time_eval, // in
                           csys, cb);                                  // out

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 2
      if (cs_dbg_cw_test(cm)) cs_cell_mesh_dump(cm);
#endif

      /* DIFFUSION TERM */
      /* ============== */

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        if (!(eqb->diff_pty_uniform))
          cs_equation_set_diffusion_property_cw(eqp, cm, time_eval, cell_flag,
                                                cb);

        /* local matrix owned by the cellwise builder (store in cb->loc) */
        eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        if (eqp->diffusion_hodge.is_iso == false)
          bft_error(__FILE__, __LINE__, 0, " %s: Case not handle yet\n",
                    __func__);

        /* Add the local diffusion operator to the local system */
        const cs_real_t  *sval = cb->loc->val;
        for (int bi = 0; bi < cm->n_vc; bi++) {
          for (int bj = 0; bj < cm->n_vc; bj++) {

            /* Retrieve the 3x3 matrix */
            cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bi, bj);
            assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

            const cs_real_t  _val = sval[cm->n_vc*bi+bj];
            bij->val[0] += _val;
            bij->val[4] += _val;
            bij->val[8] += _val;

          }
        }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
        if (cs_dbg_cw_test(cm))
          cs_cell_sys_dump("\n>> Cell system after diffusion", csys);
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
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        time_eval,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        for (short int v = 0; v < cm->n_vc; v++)
          csys->rhs[v] += csys->source[v];

      } /* End of term source */

      /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ======================================================= */

      if (cell_flag & CS_FLAG_BOUNDARY) {

        if (csys->has_dirichlet)
          /* Weakly enforced Dirichlet BCs for cells attached to the boundary
             csys is updated inside (matrix and rhs) */
        eqc->enforce_dirichlet(eqp, cm, eqc->bdy_flux_op, fm, cb, csys);

      } /* Boundary cell */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 0
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY */
      /* ======== */

      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_VECT];

      /* Matrix assembly */
      cs_equation_assemble_block_matrix(csys, rs, 3, mav);

      /* Assemble RHS */
      for (short int i = 0; i < 3*cm->n_vc; i++) {
#       pragma omp atomic
        rhs[csys->dof_ids[i]] += csys->rhs[i];
      }

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); // optional

  /* TODO */
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
 * \param[in, out] data       pointer to cs_cdovb_vecteq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_update_field(const cs_real_t            *solu,
                             const cs_real_t            *rhs,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data,
                             cs_real_t                  *field_val)
{
  CS_UNUSED(rhs);
  CS_UNUSED(eqp);

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t  *)data;
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
 * \brief  Retrieve an array of values at mesh vertices for the variable field
 *         associated to the given context
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_vecteq_get_vertex_values(void      *context)
{
  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  return pot->val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute an array of values at mesh cells by interpolating the
 *         variable field associated to the given context located at mesh
 *         vertices
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovb_vecteq_get_cell_values(void      *context)
{
  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Reset buffer of values */
  if (eqc->cell_values == NULL)
    BFT_MALLOC(eqc->cell_values, 3*quant->n_cells, cs_real_t);
  memset(eqc->cell_values, 0, 3*quant->n_cells*sizeof(cs_real_t));

  /* Compute the values at cell centers from an interpolation of the field
     values defined at vertices */
  cs_reco_vect_pv_at_cell_centers(connect->c2v, quant, pot->val,
                                  eqc->cell_values);

  return eqc->cell_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *
 * \param[in]       normal     indicate in which direction flux is > 0
 * \param[in]       pdi        pointer to an array of field values
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to data specific for this scheme
 * \param[in, out]  d_flux     pointer to the value of the diffusive flux
 * \param[in, out]  c_flux     pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_compute_flux_across_plane(const cs_real_t             normal[],
                                          const cs_real_t            *pdi,
                                          int                         ml_id,
                                          const cs_equation_param_t  *eqp,
                                          cs_equation_builder_t      *eqb,
                                          void                       *data,
                                          double                     *d_flux,
                                          double                     *c_flux)
{
  CS_UNUSED(data);

  cs_mesh_location_type_t  ml_t = cs_mesh_location_get_type(ml_id);

  *d_flux = 0.;
  *c_flux = 0.;

  if (pdi == NULL)
    return;

  if (ml_t != CS_MESH_LOCATION_INTERIOR_FACES &&
      ml_t != CS_MESH_LOCATION_BOUNDARY_FACES) {
    cs_base_warn(__FILE__, __LINE__);
    cs_log_printf(CS_LOG_DEFAULT,
                  _(" Mesh location type is incompatible with the computation\n"
                    " of the flux across faces.\n"));
    return;
  }

  const cs_timer_t  t0 = cs_timer_time();
  const cs_lnum_t  *n_elts = cs_mesh_location_get_n_elts(ml_id);
  const cs_lnum_t  *elt_ids = cs_mesh_location_get_elt_list(ml_id);

  if (cs_glob_n_ranks == 1)
    if (n_elts[0] > 0 && elt_ids == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _(" Computing the flux across all interior or border faces is"
                  " not managed yet."));

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_adjacency_t  *f2c = connect->f2c;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* TODO */
  CS_UNUSED(quant);
  CS_UNUSED(f2c);
  CS_UNUSED(normal);
  CS_UNUSED(eqp);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  data        pointer to cs_cdovb_vecteq_t structure
 * \param[in, out]  location    where the flux is defined
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_cellwise_diff_flux(const cs_real_t             *values,
                                   const cs_equation_param_t   *eqp,
                                   cs_real_t                    t_eval,
                                   cs_equation_builder_t       *eqb,
                                   void                        *data,
                                   cs_flag_t                    location,
                                   cs_real_t                   *diff_flux)
{
  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t  *)data;

  /* Sanity checks */
  assert(diff_flux != NULL && eqp != NULL && eqc != NULL && eqb != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  if (!cs_flag_test(location, cs_flag_primal_cell) &&
      !cs_flag_test(location, cs_flag_dual_face_byc))
    bft_error(__FILE__, __LINE__, 0,
              " Incompatible location.\n"
              " Stop computing a cellwise diffusive flux.");

  const cs_timer_t  t0 = cs_timer_time();

  /* TODO */
  CS_UNUSED(quant);
  CS_UNUSED(connect);
  CS_UNUSED(values);
  CS_UNUSED(t_eval);
  CS_UNUSED(eqp);

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
 * \param[in, out]  data       pointer to cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data)
{
  CS_UNUSED(field);
  CS_UNUSED(data);

  const cs_timer_t  t0 = cs_timer_time();

  /* TODO */
  CS_UNUSED(eqname);
  CS_UNUSED(eqp);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
