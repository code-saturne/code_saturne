/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
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
 *        convection-diffusion-reaction of scalar-valued equations with
 *        source terms
 *
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVB_VECTEQ_DBG     0
#define CS_CDOVB_VECTEQ_MODULO  100

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

  /* TODO */
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_vecteq_t structure storing data useful for
 *         managing such a scheme
 *
 * \param[in]      eqp   pointer to a cs_equation_param_t structure
 * \param[in, out] eqb   pointer to a cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_vecteq_init_context(const cs_equation_param_t   *eqp,
                             cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVB && eqp->dim != 3)
    bft_error(__FILE__, __LINE__, 0,
              " Invalid type of equation.\n"
              " Expected: vector-valued CDO vertex-based equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->n_vertices;

  cs_cdovb_vecteq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdovb_vecteq_t);

  eqc->n_dofs = 3*n_vertices;

  /* Flag to indicate the minimal set of quantities to build in a cell mesh
     According to the situation, additional flags have to be set */
  eqb->msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PE |
    CS_CDO_LOCAL_EV;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE;

  /* Diffusion part */
  /* -------------- */

  eqc->get_stiffness_matrix = NULL;
  eqc->boundary_flux_op = NULL;
  eqc->enforce_dirichlet = NULL;

  if (cs_equation_param_has_diffusion(eqp)) {

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_cost_get_stiffness;
      eqc->boundary_flux_op = cs_cdovb_diffusion_cost_flux_op;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_voro_get_stiffness;
      eqc->boundary_flux_op = cs_cdovb_diffusion_cost_flux_op;
      break;

    case CS_PARAM_HODGE_ALGO_WBS:
      eqb->msh_flag |= CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_PEQ |
        CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_wbs_get_stiffness;
      eqc->boundary_flux_op = cs_cdovb_diffusion_wbs_flux_op;
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
      eqb->bd_msh_flag |= CS_CDO_LOCAL_PFQ|CS_CDO_LOCAL_DEQ|CS_CDO_LOCAL_FEQ;
      eqc->enforce_dirichlet = cs_cdovb_diffusion_weak_dirichlet;
      break;

    case CS_PARAM_BC_ENFORCE_WEAK_SYM:
      eqb->bd_msh_flag |= CS_CDO_LOCAL_PFQ|CS_CDO_LOCAL_DEQ|CS_CDO_LOCAL_FEQ;
      eqc->enforce_dirichlet = cs_cdovb_diffusion_wsym_dirichlet;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                (" Invalid type of algorithm to enforce Dirichlet BC."));

    }

  } // Has diffusion

  /* Advection part */
  /* -------------- */

  eqc->get_advection_matrix = NULL;
  eqc->add_advection_bc = NULL;

  if (cs_equation_param_has_convection(eqp)) {

    cs_xdef_type_t  adv_deftype =
      cs_advection_field_get_deftype(eqp->adv_field);

    if (adv_deftype == CS_XDEF_BY_VALUE)
      eqb->msh_flag |= CS_CDO_LOCAL_DFQ;
    else if (adv_deftype == CS_XDEF_BY_ARRAY)
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ;
    else if (adv_deftype == CS_XDEF_BY_ANALYTIC_FUNCTION)
      eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_EFQ | CS_CDO_LOCAL_PFQ;

    switch (eqp->adv_formulation) {

    case CS_PARAM_ADVECTION_FORM_CONSERV:

      switch (eqp->adv_scheme) {

      case CS_PARAM_ADVECTION_SCHEME_CENTERED:
        eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
        eqc->get_advection_matrix = cs_cdo_advection_get_vb_cencsv;
        break;

      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
      case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
      case CS_PARAM_ADVECTION_SCHEME_SG:
        eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
        if (cs_equation_param_has_diffusion(eqp))
          eqc->get_advection_matrix = cs_cdo_advection_get_vb_upwcsvdi;
        else
          eqc->get_advection_matrix = cs_cdo_advection_get_vb_upwcsv;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid advection scheme for vertex-based discretization");
      } // Scheme
      break; // Formulation

    case CS_PARAM_ADVECTION_FORM_NONCONS:

      switch (eqp->adv_scheme) {
      case CS_PARAM_ADVECTION_SCHEME_CENTERED:
        eqc->get_advection_matrix = cs_cdo_advection_get_vb_cennoc;
        break;

      case CS_PARAM_ADVECTION_SCHEME_UPWIND:
      case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
      case CS_PARAM_ADVECTION_SCHEME_SG:
        eqb->msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
        if (cs_equation_param_has_diffusion(eqp))
          eqc->get_advection_matrix = cs_cdo_advection_get_vb_upwnocdi;
        else
          eqc->get_advection_matrix = cs_cdo_advection_get_vb_upwnoc;
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid advection scheme for vertex-based discretization");
      } // Scheme
      break; // Formulation

    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid type of formulation for the advection term");
    }

    /* Boundary conditions for advection */
    eqb->bd_msh_flag |= CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_FEQ;
    eqc->add_advection_bc = cs_cdo_advection_add_vb_bc;

  }
  else {

    if (eqp->enforcement != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
      eqb->sys_flag |= CS_FLAG_SYS_SYM; // Algebraic system is symmetric

  }

  /* Reaction part */
  /* ------------- */

  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->reaction_hodge.algo == CS_PARAM_HODGE_ALGO_WBS) {
      eqb->msh_flag |= CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FEQ |
        CS_CDO_LOCAL_HFQ;
      eqb->sys_flag |= CS_FLAG_SYS_HLOC_CONF;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " Invalid choice of algorithm for the reaction term.");

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
        eqb->sys_flag |= CS_FLAG_SYS_HLOC_CONF;
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
  eqc->hdg_mass.coef = 1.0; // not useful in this case

  eqc->get_mass_matrix = cs_hodge_vpcd_wbs_get;

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

  /* Last free */
  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the contributions of source terms (store inside builder)
 *
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_compute_source(const cs_equation_param_t  *eqp,
                               cs_equation_builder_t      *eqb,
                               void                       *data)
{
  if (data == NULL)
    return;

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)data;

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  if (eqp->n_source_terms == 0)
    return;

  cs_timer_t  t0 = cs_timer_time();

  /* Compute the source term cell by cell */
#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, eqc, cs_cdovb_cell_sys, cs_cdovb_cell_bld)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdovb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdovb_cell_bld[t_id];
    cs_flag_t  msh_flag = eqb->st_msh_flag;

#if defined(DEBUG) && !defined(NDEBUG)
    cs_cell_mesh_reset(cm);
#endif

    if (cs_equation_param_has_sourceterm(eqp)) {
      /* Reset source term array */
#     pragma omp for CS_CDO_OMP_SCHEDULE
      for (cs_lnum_t i = 0; i < eqc->n_dofs; i++) eqc->source_terms[i] = 0;
    }

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t  c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 2
      if (c_id % CS_CDOVB_VECTEQ_MODULO == 0) cs_cell_mesh_dump(cm);
#endif

      /* Build the local dense matrix related to this operator
         Store in cb->hdg inside the cs_cell_builder_t structure */
      if (eqb->sys_flag & CS_FLAG_SYS_SOURCES_HLOC)
        eqc->get_mass_matrix(eqc->hdg_mass, cm, cb);

      /* Initialize the local number of DoFs */
      csys->n_dofs = cm->n_vc;

      /* Compute the contribution of all source terms in each cell */
      cs_source_term_compute_cellwise(eqp->n_source_terms,
                  (const cs_xdef_t **)eqp->source_terms,
                                      cm,
                                      eqb->source_mask,
                                      eqb->compute_source,
                                      NULL,  // No input structure
                                      cb,    // mass matrix is stored in cb->hdg
                                      csys); // Fill csys->source

      /* Assemble the cellwise contribution to the rank contribution */
      for (short int v = 0; v < cm->n_vc; v++)
#       pragma omp atomic
        eqc->source_terms[cm->v_ids[v]] += csys->source[v];

    } // Loop on cells

  } // OpenMP block

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 2
    cs_dbg_darray_to_listing("INIT_SOURCE_TERM_VTX", quant->n_vertices,
                             eqc->source_terms, 8);
#endif

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t0, &t1);
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
 * \brief  Build the linear system arising from a scalar convection/diffusion
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

  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)data;

  /* Compute the values of the Dirichlet BC */
  cs_real_t  *dir_values =
    cs_equation_compute_dirichlet_vb(mesh,
                                     quant,
                                     connect,
                                     cs_shared_time_step,
                                     eqp,
                                     eqb->face_bc->dir,
                                     cs_cdovb_cell_bld[0]);

  /* Tag faces with a non-homogeneous Neumann BC */
  short int  *neu_tags = cs_equation_tag_neumann_face(quant, eqp);

  /* TODO */

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

  /* TODO */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
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

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  data        pointer to cs_cdovb_vecteq_t structure
 * \param[in, out]  location    where the flux is defined
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_cellwise_diff_flux(const cs_real_t             *values,
                                   const cs_equation_param_t   *eqp,
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

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
