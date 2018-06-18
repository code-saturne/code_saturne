/*============================================================================
 * Build an algebraic CDO vertex+cell-based system for unsteady convection
 * diffusion reaction of scalar-valued equations with source terms
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
#include "cs_defs.h"
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
#include "cs_static_condensation.h"
#include "cs_timer.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovcb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_cdovcb_scaleq.c

  \brief Build an algebraic CDO vertex+cell-based system for unsteady
         convection diffusion reaction of scalar-valued equations with source
         terms
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVCB_SCALEQ_DBG       0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/* Algebraic system for CDO vertex+cell-based discretization */

struct _cs_cdovcb_scaleq_t {

  /* System size */
  cs_lnum_t    n_dofs;        // n_vertices + n_cells

  /* Store the values of the field at cell centers and the data needed to
     compute the cell values from the vertex values. No need to synchronize
     all these quantities since they are only cellwise quantities. */
  cs_real_t   *cell_values;
  cs_real_t   *cell_rhs;  // right-hand side related to cell dofs

  /* Members related to the static condensation */
  cs_real_t   *rc_tilda;   /* Acc^-1 * RHS_cell */
  cs_real_t   *acv_tilda;  /* Acc^-1 * Acv
                              Cell-vertices lower-Left block of the full matrix
                              Access to the values thanks to the c2v
                              connectivity */

  /* Array storing the value of the contribution of all source terms */
  cs_real_t   *source_terms;

  /* Pointer of function to build the diffusion term */
  cs_hodge_t                      *get_stiffness_matrix;
  cs_cdo_diffusion_enforce_dir_t  *enforce_dirichlet;
  cs_cdo_diffusion_flux_trace_t   *boundary_flux_op;

  /* Pointer of function to build the advection term */
  cs_cdo_advection_t              *get_advection_matrix;
  cs_cdo_advection_bc_t           *add_advection_bc;

  /* Pointer of function to apply the time scheme */
  cs_cdo_time_scheme_t            *apply_time_scheme;

  /* If one needs to build a local hodge op. for time and reaction */
  cs_param_hodge_t                 hdg_mass;
  cs_hodge_t                      *get_mass_matrix;

};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */
static cs_cell_sys_t      **cs_cdovcb_cell_sys = NULL;
static cs_cell_builder_t  **cs_cdovcb_cell_bld = NULL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;
static const cs_matrix_structure_t  *cs_shared_ms;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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
  const int  n_vc = connect->n_max_vbyc;
  const int  n_ec = connect->n_max_ebyc;
  const int  n_fc = connect->n_max_fbyc;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  int  size = n_vc + 1;
  BFT_MALLOC(cb->ids, size, short int);
  memset(cb->ids, 0, size*sizeof(short int));

  size = 2*n_vc + 3*n_ec + n_fc;
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(cs_real_t));

  size = 2*n_ec + n_vc;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of
     operators */
  cb->hdg = cs_sdm_square_create(n_vc + 1);
  cb->loc = cs_sdm_square_create(n_vc + 1);
  cb->aux = cs_sdm_square_create(n_vc + 1);

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
 * \param[in]      eqc         pointer to a cs_cdovb_scaleq_t structure
 * \param[in]      dir_values  Dirichlet values associated to each vertex
 * \param[in]      neu_tags    Definition id related to each Neumann face
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
                  const cs_cdovcb_scaleq_t     *eqc,
                  const cs_real_t               dir_values[],
                  const short int               neu_tags[],
                  const cs_real_t               field_tn[],
                  cs_real_t                     t_eval,
                  cs_cell_sys_t                *csys,
                  cs_cell_builder_t            *cb)
{
  /* Cell-wise view of the linear system to build */
  const short int  n_dofs = cm->n_vc + 1;

  /* Initialize the local system */
  cs_cell_sys_reset(cell_flag, n_dofs, cm->n_fc, csys);

  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;
  csys->face_shift = cs_shared_connect->n_faces[2]; /* shift = n_i_faces */

  cs_sdm_square_init(n_dofs, csys->mat);

  for (short int v = 0; v < cm->n_vc; v++) {
    csys->dof_ids[v] = cm->v_ids[v];
    csys->val_n[v] = field_tn[cm->v_ids[v]];
  }
  csys->dof_ids[cm->n_vc] = cm->c_id;
  csys->val_n[cm->n_vc] = eqc->cell_values[cm->c_id];

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

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Update the value of stabilization coefficient in given situation
 *
 * \param[in]  eqp    pointer to a cs_equation_param_t  structure
 */
/*----------------------------------------------------------------------------*/

static void
_set_cip_coef(const cs_equation_param_t  *eqp)
{
  const double  gseed = 1e-2;  /* Default value to multiply according to the
                                  problem and the ratio of diameters */

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const double  hc_max = quant->cell_info.h_max;
  const double  hc_min = quant->cell_info.h_min;
  const double  hf_max = quant->face_info.h_max;
  const double  hf_min = quant->face_info.h_min;
  const double  hcMm = hc_max * hc_min;
  const double  hfMm = hf_min * hf_max;
  const double  rho_fc = hcMm / hfMm;

  double  gamma = 10 * gseed * hc_max * hc_max * rho_fc;

  /* If not pure convection */
  if (cs_equation_param_has_diffusion(eqp) ||
      cs_equation_param_has_reaction(eqp) ||
      cs_equation_param_has_time(eqp))
    gamma *= 0.1;

  cs_cdo_advection_set_cip_coef(gamma);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-vertex+cell
 *           based scheme are allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdovcb_scaleq_is_initialized(void)
{
  if (cs_cdovcb_cell_sys == NULL || cs_cdovcb_cell_bld == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Allocate work buffer and general structures related to CDO
 *           vertex+cell-based schemes
 *           Set shared pointers.
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_init_common(const cs_cdo_quantities_t    *quant,
                             const cs_cdo_connect_t       *connect,
                             const cs_time_step_t         *time_step,
                             const cs_matrix_structure_t  *ms)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
  cs_shared_ms = ms;

  /* Specific treatment for handling openMP */
  BFT_MALLOC(cs_cdovcb_cell_sys, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(cs_cdovcb_cell_bld, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_cdovcb_cell_sys[i] = NULL;
    cs_cdovcb_cell_bld[i] = NULL;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cdovcb_cell_sys[t_id] = cs_cell_sys_create(connect->n_max_vbyc + 1,
                                                  connect->n_max_fbyc,
                                                  1, NULL);
    cs_cdovcb_cell_bld[t_id] = _cell_builder_create(connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cdovcb_cell_sys[0] = cs_cell_sys_create(connect->n_max_vbyc + 1,
                                             connect->n_max_fbyc,
                                             1, NULL);
  cs_cdovcb_cell_bld[0] = _cell_builder_create(connect);
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
cs_cdovcb_scaleq_get(cs_cell_sys_t       **csys,
                     cs_cell_builder_t   **cb)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  *csys = cs_cdovcb_cell_sys[t_id];
  *cb = cs_cdovcb_cell_bld[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO vertex-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_finalize_common(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(cs_cdovcb_cell_sys[t_id]));
    cs_cell_builder_free(&(cs_cdovcb_cell_bld[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(cs_cdovcb_cell_sys[0]));
  cs_cell_builder_free(&(cs_cdovcb_cell_bld[0]));
#endif /* openMP */

  BFT_FREE(cs_cdovcb_cell_sys);
  BFT_FREE(cs_cdovcb_cell_bld);
  cs_cdovcb_cell_bld = NULL;
  cs_cdovcb_cell_sys = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_scaleq_t structure storing data useful for
 *         managing such a scheme
 *
 * \param[in] eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb   pointer to a cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovcb_scaleq_init_context(const cs_equation_param_t   *eqp,
                              cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVCB && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Invalid type of equation.\n"
              " Expected: scalar-valued CDO vertex+cell-based equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->n_vertices;
  const cs_lnum_t  n_cells = connect->n_cells;

  cs_cdovcb_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdovcb_scaleq_t);

  /* System dimension */
  eqc->n_dofs = n_vertices + n_cells;

  /* Flag to indicate what to build in a cell mesh */
  eqb->msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_DEQ |
    CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_EV |
    CS_CDO_LOCAL_FE  | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;

  /* Store the last computed values of the field at cell centers and the data
     needed to compute the cell values from the vertex values.
     No need to synchronize all these quantities since they are only cellwise
     quantities. */
  BFT_MALLOC(eqc->cell_values, n_cells, cs_real_t);
  BFT_MALLOC(eqc->rc_tilda, n_cells, cs_real_t);
  BFT_MALLOC(eqc->acv_tilda, connect->c2v->idx[n_cells], cs_real_t);

  memset(eqc->cell_values, 0, sizeof(cs_real_t)*n_cells);
  memset(eqc->rc_tilda, 0, sizeof(cs_real_t)*n_cells);
  memset(eqc->acv_tilda, 0, sizeof(cs_real_t)*connect->c2v->idx[n_cells]);

  /* Diffusion part */
  eqc->enforce_dirichlet = NULL;
  eqc->get_stiffness_matrix = NULL;
  eqc->boundary_flux_op = NULL;
  if (cs_equation_param_has_diffusion(eqp)) {

    eqc->get_stiffness_matrix = cs_hodge_vcb_get_stiffness;
    eqc->boundary_flux_op = cs_cdovcb_diffusion_flux_op;

    switch (eqp->enforcement) {

    case CS_PARAM_BC_ENFORCE_WEAK_PENA:
      eqc->enforce_dirichlet = cs_cdo_diffusion_pena_dirichlet;
      break;
    case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
      eqc->enforce_dirichlet = cs_cdovb_diffusion_weak_dirichlet;
      break;
    case CS_PARAM_BC_ENFORCE_WEAK_SYM:
      eqc->enforce_dirichlet = cs_cdovb_diffusion_wsym_dirichlet;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                (" Invalid type of algorithm to enforce Dirichlet BC."));

    }

  }

  /* Advection part */
  eqc->get_advection_matrix = NULL;
  eqc->add_advection_bc = NULL;

  if (cs_equation_param_has_convection(eqp)) {

    switch (eqp->adv_scheme) {
    case CS_PARAM_ADVECTION_SCHEME_CIP:

      eqb->msh_flag |= CS_CDO_LOCAL_EF;
      _set_cip_coef(eqp);

      eqc->add_advection_bc = cs_cdo_advection_add_vcb_bc;

      if (cs_advection_field_is_cellwise(eqp->adv_field))
        eqc->get_advection_matrix = cs_cdo_advection_get_vcb_cw;
      else
        eqc->get_advection_matrix = cs_cdo_advection_get_vcb;

      break;

    case CS_PARAM_ADVECTION_SCHEME_UPWIND:
    case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
    case CS_PARAM_ADVECTION_SCHEME_SG:
    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid advection scheme for vertex-based discretization");
    } // Scheme

  }
  else {

    if (eqp->enforcement != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
      eqb->sys_flag |= CS_FLAG_SYS_SYM; // Algebraic system is symmetric

  }

  /* Reaction part */
  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->reaction_hodge.algo == CS_PARAM_HODGE_ALGO_WBS)
      eqb->sys_flag |= CS_FLAG_SYS_HLOC_CONF;
    else
      bft_error(__FILE__, __LINE__, 0,
                " Invalid choice of algorithm for the reaction term.");

  } /* Reaction */

  /* Time part */
  if (cs_equation_param_has_time(eqp)) {

    if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_VORONOI)
      eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
    else if (eqp->time_hodge.algo == CS_PARAM_HODGE_ALGO_WBS) {
      if (eqp->do_lumping)
        eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
      else
        eqb->sys_flag |= CS_FLAG_SYS_HLOC_CONF;
    }

  }

  eqc->apply_time_scheme = cs_cdo_time_get_scheme_function(eqb->sys_flag, eqp);

  /* Source term part */
  /* ---------------- */

  eqc->source_terms = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, eqc->n_dofs, cs_real_t);
    memset(eqc->source_terms, 0, sizeof(cs_real_t)*eqc->n_dofs);

  } /* There is at least one source term */

  // Pre-defined a cs_hodge_builder_t struct.
  eqc->hdg_mass.is_unity = true;
  eqc->hdg_mass.is_iso   = true;
  eqc->hdg_mass.inv_pty  = false;
  eqc->hdg_mass.type = CS_PARAM_HODGE_TYPE_VC;
  eqc->hdg_mass.algo = CS_PARAM_HODGE_ALGO_WBS;
  eqc->hdg_mass.coef = 1.0; // not useful in this case

  eqc->get_mass_matrix = cs_hodge_vcb_wbs_get;

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovcb_scaleq_t structure
 *
 * \param[in, out]  data   pointer to a cs_cdovcb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovcb_scaleq_free_context(void   *data)
{
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t *)data;

  if (eqc == NULL)
    return eqc;

  BFT_FREE(eqc->cell_values);
  BFT_FREE(eqc->rc_tilda);
  BFT_FREE(eqc->acv_tilda);

  BFT_FREE(eqc->source_terms);

  /* Last free */
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
 * \param[in, out] context        pointer to cs_cdovcb_scaleq_t structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_initialize_system(const cs_equation_param_t    *eqp,
                                   cs_equation_builder_t        *eqb,
                                   void                         *context,
                                   cs_matrix_t                **system_matrix,
                                   cs_real_t                  **system_rhs)
{
  CS_UNUSED(eqp);

  if (context == NULL)
    return;

  assert(*system_matrix == NULL && *system_rhs == NULL);

  cs_timer_t  t0 = cs_timer_time();

  /* Create the matrix related to the current algebraic system */
  *system_matrix = cs_matrix_create(cs_shared_ms);

  const cs_lnum_t  n_v = cs_shared_quant->n_vertices;

  BFT_MALLOC(*system_rhs, n_v, cs_real_t);
  memset(*system_rhs, 0, sizeof(cs_real_t)*n_v);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings when the fields
 *         stem from a scalar CDO vertex+cell-based scheme.
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in, out] eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval      time at which one evaluates BCs
 * \param[in, out] field_val   pointer to the values of the variable field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_set_dir_bc(const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            cs_real_t                     t_eval,
                            cs_real_t                     field_val[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Compute the values of the Dirichlet BC */
  cs_equation_compute_dirichlet_vb(mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_eval,
                                   cs_cdovcb_cell_bld[0], /* static variable */
                                   field_val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO vertex+cell-based scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the vertex field
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_cdovcb_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_build_system(const cs_mesh_t            *mesh,
                              const cs_real_t            *field_val,
                              double                      dt_cur,
                              const cs_equation_param_t  *eqp,
                              cs_equation_builder_t      *eqb,
                              void                       *data,
                              cs_real_t                  *rhs,
                              cs_matrix_t                *matrix)
{
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t *)data;

  /* Sanity checks */
  assert(rhs != NULL && matrix != NULL);
  assert(eqc != NULL && eqb != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;

  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* Compute the values of the Dirichlet BC. */
  cs_real_t  *dir_values = NULL;
  BFT_MALLOC(dir_values, quant->n_vertices, cs_real_t);
  memset(dir_values, 0, quant->n_vertices*sizeof(cs_real_t));

  cs_cdovcb_scaleq_set_dir_bc(mesh, eqp, eqb, t_cur + dt_cur, dir_values);

  /* Tag faces with a non-homogeneous Neumann BC */
  short int  *neu_tags = cs_equation_tag_neumann_face(quant, eqp);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)          \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav, dir_values, \
         neu_tags, field_val, cs_cdovcb_cell_sys, cs_cdovcb_cell_bld)
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
    cs_cell_sys_t  *csys = cs_cdovcb_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_cdovcb_cell_bld[t_id];
    cs_real_t  *cell_sources = eqc->source_terms + quant->n_vertices;

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
      if (cs_dbg_cw_test(cm))
        cs_cell_mesh_dump(cm);
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

        /* Weakly enforced Dirichlet BCs for cells attached to the boundary
           csys is updated inside (matrix and rhs) */
        if (cell_flag & CS_FLAG_BOUNDARY) {
          if (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
              eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
            eqc->enforce_dirichlet(eqp->diffusion_hodge, cm, // in
                                   eqc->boundary_flux_op,    // function
                                   fm, cb, csys);            // in/out
        } // Border cell

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
        if (cs_dbg_cw_test(cm))
          cs_cell_sys_dump("\n>> Local system after diffusion", c_id, csys);
#endif
      } /* END OF DIFFUSION */

      /* ADVECTION TERM */
      /* ============== */

      if (cs_equation_param_has_convection(eqp)) {

        /* Define the local advection matrix */
        eqc->get_advection_matrix(eqp, cm, t_eval_pty, fm, cb);

        cs_sdm_add(csys->mat, cb->loc);

        /* Last treatment for the advection term: Apply boundary conditions
           csys is updated inside (matrix and rhs) */
        if (cell_flag & CS_FLAG_BOUNDARY)
          eqc->add_advection_bc(cm, eqp, t_eval_pty, fm, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
        if (cs_dbg_cw_test(cm))
          cs_cell_sys_dump("\n>> Local system after advection", c_id, csys);
#endif
      } /* END OF ADVECTION */

      if (eqb->sys_flag & CS_FLAG_SYS_HLOC_CONF)
        eqc->get_mass_matrix(eqc->hdg_mass, cm, cb); // stored in cb->hdg

      /* REACTION TERM */
      /* ============= */

      if (cs_equation_param_has_reaction(eqp)) {

        /* Define the local reaction property */
        double  rpty_val = 0;
        for (int r = 0; r < eqp->n_reaction_terms; r++)
          if (eqb->reac_pty_uniform[r])
            rpty_val += reac_pty_vals[r];
          else
            rpty_val += cs_property_get_cell_value(c_id, t_eval_pty,
                                                   eqp->reaction_properties[r]);

        /* Update local system matrix with the reaction term
           cb->hdg corresponds to the current mass matrix */
        cs_sdm_add_mult(csys->mat, rpty_val, cb->hdg);

      } /* END OF REACTION */

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
                                        NULL,  /* No data structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        for (short int v = 0; v < cm->n_vc; v++)
          csys->rhs[v] += csys->source[v];
        csys->rhs[cm->n_vc] += csys->source[cm->n_vc];

      } /* End of term source */

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

        cs_sdm_t  *mass_mat = cb->hdg;
        if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) {

          /* Switch to cb->loc. Define a diagonal matrix (seen as a vector) */
          mass_mat = cb->loc;

          /* 0.75*|c|*wvc = 0.75*|dual_cell(v) cap c| for vertices
             0.25*|c|*wvc = 0.75*|dual_cell(v) cap c| for the cell */
          const double  ptyc = tpty_val * cm->vol_c;
          for (short int v = 0; v < cm->n_vc; v++)
            mass_mat->val[v] = 0.75 * ptyc * cm->wvc[v];
          mass_mat->val[cm->n_vc] = 0.25 * ptyc;

        }

        /* Apply the time discretization to the local system.
           Update csys (matrix and rhs) */
        eqc->apply_time_scheme(eqp, tpty_val, mass_mat, eqb->sys_flag, cb,
                               csys);

      } /* END OF TIME CONTRIBUTION */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> Local system matrix before condensation",
                         c_id, csys);
#endif

      /* Neumann boundary conditions (Consistent for linear solutions) */
      if ((cell_flag & CS_FLAG_BOUNDARY) && csys->has_nhmg_neumann) {
        for (short int v  = 0; v < cm->n_vc; v++)
          csys->rhs[v] += csys->neu_values[v];
      }

      /* Static condensation of the local system matrix of size n_vc + 1 into
         a matrix of size n_vc.
         Store data in rc_tilda and acv_tilda to compute the values at cell
         centers after solving the system */
      cs_static_condensation_scalar_eq(connect->c2v,
                                       eqc->rc_tilda,
                                       eqc->acv_tilda,
                                       cb, csys);

      /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ======================================================= */

      if (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_PENA) {

        /* Weakly enforced Dirichlet BCs for cells attached to the boundary
           csys is updated inside (matrix and rhs) */
        if (cell_flag & CS_FLAG_BOUNDARY)
          eqc->enforce_dirichlet(eqp->diffusion_hodge, cm,  // in
                                 eqc->boundary_flux_op,     // function
                                 fm, cb, csys);             // in/out

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", c_id, csys);
#endif

      /* ASSEMBLY */
      /* ======== */

      const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];

      /* Matrix assembly */
      cs_equation_assemble_matrix(csys, rs, mav);

      /* Assemble RHS */
      for (short int v = 0; v < cm->n_vc; v++)
#       pragma omp atomic
        rhs[cm->v_ids[v]] += csys->rhs[v];

      if (eqc->source_terms != NULL) { /* Assemble only the part related to
                                          vertices */
        for (short int v = 0; v < cm->n_vc; v++)
#         pragma omp atomic
          eqc->source_terms[cm->v_ids[v]] += csys->source[v];
      }

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); // optional

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 2
  cs_dbg_darray_to_listing("FINAL RHS_VTX", quant->n_vertices, rhs, 8);
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
 * \param[in, out] data       pointer to data structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_update_field(const cs_real_t            *solu,
                              const cs_real_t            *rhs,
                              const cs_equation_param_t  *eqp,
                              cs_equation_builder_t      *eqb,
                              void                       *data,
                              cs_real_t                  *field_val)
{
  CS_UNUSED(eqp);
  CS_UNUSED(rhs);

  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)data;
  cs_timer_t  t0 = cs_timer_time();

  /* Set the values at vertices */
  memcpy(field_val, solu, sizeof(cs_real_t) * cs_shared_quant->n_vertices);

  /* Compute values at cells pc = acc^-1*(RHS - Acv*pv) */
  cs_static_condensation_recover_scalar(cs_shared_connect->c2v,
                                        eqc->rc_tilda,
                                        eqc->acv_tilda,
                                        solu,
                                        eqc->cell_values);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at cell centers (DoF used in the linear
 *         system are located at primal vertices and field related to the
 *         structure equation is also attached to primal vertices
 *
 * \param[in]  data    pointer to a data structure
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

double *
cs_cdovcb_scaleq_get_cell_values(const void          *data)
{
  const cs_cdovcb_scaleq_t  *eqc = (const cs_cdovcb_scaleq_t  *)data;

  if (eqc == NULL)
    return NULL;
  else
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
cs_cdovcb_scaleq_compute_flux_across_plane(const cs_real_t             normal[],
                                           const cs_real_t            *pdi,
                                           int                         ml_id,
                                           const cs_equation_param_t  *eqp,
                                           cs_equation_builder_t      *eqb,
                                           void                       *data,
                                           double                     *d_flux,
                                           double                     *c_flux)
{
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)data;
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

  if (n_elts[0] > 0 && elt_ids == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _(" Computing the flux across all interior or border faces is not"
                " managed yet."));

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_adjacency_t  *f2c = connect->f2c;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;

  double  flx, p_f;
  cs_real_33_t  pty_tens;
  cs_nvec3_t  adv_c;

  // To be modified for an fully integration of openMP
  cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(0);
  cs_cell_builder_t  *cb = cs_cdovcb_cell_bld[0];

  double  *p_v = NULL;
  BFT_MALLOC(p_v, connect->n_max_vbyf, double);

  if (ml_t == CS_MESH_LOCATION_BOUNDARY_FACES) {

    const cs_lnum_t  n_i_faces = connect->n_faces[2];
    const cs_lnum_t  *cell_ids = f2c->ids + f2c->idx[n_i_faces];

    for (cs_lnum_t id = 0; id < n_elts[0]; id++) {

      const cs_lnum_t  bf_id = elt_ids[id];
      const cs_lnum_t  f_id = n_i_faces + bf_id;
      const cs_lnum_t  c_id = cell_ids[bf_id];

      /* Build a face-wise view of the mesh */
      cs_face_mesh_build(c_id, f_id, connect, quant, fm);

      const short int  sgn = (_dp3(fm->face.unitv, normal) < 0) ? -1 : 1;

      /* Store values related to this face */
      for (short int v = 0; v < fm->n_vf; v++)
        p_v[v] = pdi[fm->v_ids[v]];

      p_f = cs_reco_fw_scalar_pv_at_face_center(fm, p_v);

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Compute the local diffusive flux */
        cs_property_get_cell_tensor(c_id, t_cur,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);

        flx = cs_cdo_diffusion_face_wbs_flux(fm,
                                             (const cs_real_3_t (*))pty_tens,
                                             p_v, p_f, eqc->cell_values[c_id],
                                             cb);
        *d_flux += sgn * flx;

      } // Diffusive flux

      if (cs_equation_param_has_convection(eqp)) {

        /* Compute the local advective flux */
        const double  coef = sgn * fm->face.meas * p_f;

        cs_advection_field_get_cell_vector(c_id, eqp->adv_field, &adv_c);
        *c_flux += coef * adv_c.meas * _dp3(adv_c.unitv, fm->face.unitv);

      }

    } // Loop on selected border faces

  }
  else if (ml_t == CS_MESH_LOCATION_INTERIOR_FACES) {

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  f_id = elt_ids[i];

      for (cs_lnum_t j = f2c->idx[f_id]; j < f2c->idx[f_id+1]; j++) {

        const cs_lnum_t  c_id = f2c->ids[j];

        /* Build a face-wise view of the mesh */
        cs_face_mesh_build(c_id, f_id, connect, quant, fm);

        const short int  sgn = (_dp3(fm->face.unitv, normal) < 0) ? -1 : 1;

        /* Store values related to this face */
        for (short int v = 0; v < fm->n_vf; v++)
          p_v[v] = pdi[fm->v_ids[v]];

        p_f = cs_reco_fw_scalar_pv_at_face_center(fm, p_v);

        if (cs_equation_param_has_diffusion(eqp)) {

          /* Compute the diffusive flux seen from cell c1 */
          cs_property_get_cell_tensor(c_id, t_cur,
                                      eqp->diffusion_property,
                                      eqp->diffusion_hodge.inv_pty,
                                      pty_tens);

          flx = cs_cdo_diffusion_face_wbs_flux(fm,
                                               (const cs_real_3_t (*))pty_tens,
                                               p_v, p_f, eqc->cell_values[c_id],
                                               cb);

          *d_flux += sgn * 0.5 * flx;

        } // Diffusive flux

        if (cs_equation_param_has_convection(eqp)){ /* Centered flux */

          /* Compute the local advective flux seen from cell */
          cs_advection_field_get_cell_vector(c_id,
                                             eqp->adv_field,
                                             &adv_c);
          flx = adv_c.meas * _dp3(adv_c.unitv, fm->face.unitv);

          *c_flux += sgn * 0.5 * flx  * p_f * fm->face.meas;

        } // Advective flux

      } // Loop on cells attached to this interior face

    } // Loop on selected interior faces

  } // Set of interior or border faces

  BFT_FREE(p_v);

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
 * \param[in, out]  context     pointer to data specific for this scheme
 * \param[in, out]  location    where the flux is defined
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_cellwise_diff_flux(const cs_real_t             *values,
                                    const cs_equation_param_t   *eqp,
                                    cs_real_t                    t_eval,
                                    cs_equation_builder_t       *eqb,
                                    void                        *context,
                                    cs_flag_t                    location,
                                    cs_real_t                   *diff_flux)
{
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)context;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Sanity checks */
  assert(diff_flux != NULL);
  assert(eqp->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_WBS);
  if (!cs_flag_test(location, cs_flag_primal_cell) &&
      !cs_flag_test(location, cs_flag_dual_face_byc))
    bft_error(__FILE__, __LINE__, 0,
              "Incompatible location.\n"
              " Stop computing a cellwise diffusive flux.");

  if (cs_equation_param_has_diffusion(eqp) == false) {

    size_t  size = 0;
    if (cs_flag_test(location, cs_flag_primal_cell))
      size = 3*quant->n_cells;
    else if (cs_flag_test(location, cs_flag_dual_face_byc))
      size = connect->c2e->idx[quant->n_cells];

#   pragma omp parallel for if (size > CS_THR_MIN)
    for (size_t i = 0; i < size; i++)
      diff_flux[i] = 0;

    return;
  }

  cs_timer_t  t0 = cs_timer_time();

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)   \
  shared(quant, connect, location, eqp, eqb, eqc, diff_flux, values,  \
         t_eval, cs_cdovcb_cell_bld)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif
    double  *pot = NULL;

    BFT_MALLOC(pot, connect->n_max_vbyc + 1, double);

    cs_cdo_cellwise_diffusion_flux_t  *compute_flux = NULL;
    cs_flag_t  msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PFQ |
      CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV;

    if (cs_flag_test(location, cs_flag_primal_cell)) {
      compute_flux = cs_cdo_diffusion_wbs_get_pc_flux;
      msh_flag |= CS_CDO_LOCAL_HFQ;
    }
    else if (cs_flag_test(location, cs_flag_dual_face_byc)) {
      compute_flux = cs_cdo_diffusion_wbs_get_dfbyc_flux;
      msh_flag |= CS_CDO_LOCAL_EFQ;
    }

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = cs_cdovcb_cell_bld[t_id];

    if (eqb->diff_pty_uniform)
      cs_equation_set_diffusion_property(eqp, 0, t_eval, 0, cb);

    /* Define the flux by cellwise contributions */
#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      if (!eqb->diff_pty_uniform)
        cs_equation_set_diffusion_property_cw(eqp, cm, t_eval, 0, cb);

      /* Define a local buffer keeping the value of the discrete potential
         for the current cell */
      for (short int v = 0; v < cm->n_vc; v++)
        pot[v] = values[cm->v_ids[v]];
      pot[cm->n_vc] = eqc->cell_values[c_id];

      if (cs_flag_test(location, cs_flag_primal_cell))
        compute_flux(cm, pot, cb, diff_flux + 3*c_id);
      else if (cs_flag_test(location, cs_flag_dual_face_byc))
        compute_flux(cm, pot, cb, diff_flux + connect->c2e->idx[c_id]);

    } // Loop on cells

    BFT_FREE(pot);

  } // OMP Section

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the discrete gradient at vertices
 *
 * \param[in]       v_values    discrete values for the potential at vertices
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to data structure
 * \param[in, out]  v_gradient  gradient at vertices
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_vtx_gradient(const cs_real_t         *v_values,
                              cs_equation_builder_t   *eqb,
                              void                    *context,
                              cs_real_t               *v_gradient)
{
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)context;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  if (v_gradient == NULL)
    bft_error(__FILE__, __LINE__, 0,
              " Result array has to be allocated prior to the call.");

  cs_real_t  *dualcell_vol = NULL;
  BFT_MALLOC(dualcell_vol, quant->n_vertices, cs_real_t);

# pragma omp parallel for if (3*quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*quant->n_vertices; i++)
    v_gradient[i]  = 0;
# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < quant->n_vertices; i++)
    dualcell_vol[i] = 0;

  cs_timer_t  t0 = cs_timer_time();

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)  \
  shared(quant, connect, eqc, v_gradient, v_values, dualcell_vol, \
         cs_cdovcb_cell_bld, cs_glob_n_ranks)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif
    double  *pot = NULL;
    cs_real_3_t  cgrd;

    BFT_MALLOC(pot, connect->n_max_vbyc + 1, double);

    cs_flag_t  msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PFQ |
      CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_HFQ;

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = cs_cdovcb_cell_bld[t_id];

    /* Define the flux by cellwise contributions */
#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Define a local buffer keeping the value of the discrete potential
         for the current cell */
      for (short int v = 0; v < cm->n_vc; v++)
        pot[v] = v_values[cm->v_ids[v]];
      pot[cm->n_vc] = eqc->cell_values[c_id];

      cs_reco_cw_cgrd_wbs_from_pvc(cm, pot, cb, cgrd);

      for (short int v = 0; v < cm->n_vc; v++) {
        const double dvol = cm->wvc[v] * cm->vol_c;
#       pragma omp atomic
        dualcell_vol[cm->v_ids[v]] += dvol;
        for (int k = 0; k < 3; k++)
#         pragma omp atomic
          v_gradient[3*cm->v_ids[v] + k] += dvol*cgrd[k];
      }

    } // Loop on cells

    if (cs_glob_n_ranks > 1) {

      cs_interface_set_sum(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                           connect->n_vertices,
                           1,
                           true, // interlace
                           CS_REAL_TYPE,
                           dualcell_vol);

      cs_interface_set_sum(connect->interfaces[CS_CDO_CONNECT_VTX_SCAL],
                           connect->n_vertices,
                           3,
                           true, // interlace
                           CS_REAL_TYPE,
                           v_gradient);
    }

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t i = 0; i < quant->n_vertices; i++) {
      cs_real_t  inv_dualcell_vol = 1/dualcell_vol[i];
      for (int k = 0; k < 3; k++)
        v_gradient[3*i + k] *= inv_dualcell_vol;
    }

    BFT_FREE(pot);

  } // OMP Section

  BFT_FREE(dualcell_vol);

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
 * \param[in, out]  context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_extra_op(const char                 *eqname,
                          const cs_field_t           *field,
                          const cs_equation_param_t  *eqp,
                          cs_equation_builder_t      *eqb,
                          void                       *context)
{
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)context;

  // TODO
  CS_UNUSED(field);
  CS_UNUSED(eqname);
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);
  CS_UNUSED(eqc);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
