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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
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
#include "cs_sles.h"
#include "cs_source_term.h"
#include "cs_timer.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*!
  \file cs_cdovb_scaleq.c

  \brief Build an algebraic CDO vertex-based system for unsteady
         convection-diffusion-reaction of scalar-valued equations with
         source terms

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOVB_SCALEQ_DBG     0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

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

  cs_cell_builder_t *cb = cs_cell_builder_create();

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
  cb->loc = cs_sdm_square_create(n_vc);
  cb->aux = cs_sdm_square_create(n_vc);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *
 * \param[in]      t_eval        time at which one evaluates BCs
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      eqp           pointer to a cs_equation_param_t structure
 * \param[in, out] eqb           pointer to a cs_equation_builder_t structure
 * \param[in, out] p_dir_values  pointer to the Dirichlet values to set
 * \param[in, out] p_neu_tags    pointer to the Neumann tags to set
 */
/*----------------------------------------------------------------------------*/

static void
_setup_bc(cs_real_t                     t_eval,
          const cs_mesh_t              *mesh,
          const cs_equation_param_t    *eqp,
          cs_equation_builder_t        *eqb,
          cs_real_t                    *p_dir_values[],
          short int                    *p_neu_tags[])
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_real_t  *dir_values = NULL;

  /* Compute the values of the Dirichlet BC */
  BFT_MALLOC(dir_values, quant->n_vertices, cs_real_t);
  memset(dir_values, 0, quant->n_vertices*sizeof(cs_real_t));

  /* Compute the values of the Dirichlet BC */
  cs_equation_compute_dirichlet_vb(mesh,
                                   quant,
                                   cs_shared_connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_eval,
                                   cs_cdovb_cell_bld[0], /* static variable */
                                   dir_values);
  *p_dir_values = dir_values;

  /* Tag faces with a non-homogeneous Neumann BC */
  *p_neu_tags = cs_equation_tag_neumann_face(quant, eqp);
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
_init_vb_cell_system(const cs_flag_t                cell_flag,
                     const cs_cell_mesh_t          *cm,
                     const cs_equation_param_t     *eqp,
                     const cs_equation_builder_t   *eqb,
                     const cs_real_t                dir_values[],
                     const short int                neu_tags[],
                     const cs_real_t                field_tn[],
                     cs_real_t                      t_eval,
                     cs_cell_sys_t                 *csys,
                     cs_cell_builder_t             *cb)
{
  /* Cell-wise view of the linear system to build */
  csys->c_id = cm->c_id;
  csys->cell_flag = cell_flag;
  csys->n_dofs = cm->n_vc;

  /* Initialize the local system */
  cs_cell_sys_reset(cm->n_fc, csys); /* Generic part */

  cs_sdm_square_init(cm->n_vc, csys->mat);

  for (short int v = 0; v < cm->n_vc; v++) {
    csys->dof_ids[v] = cm->v_ids[v];
    csys->val_n[v] = field_tn[cm->v_ids[v]];
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

  /* Set the diffusion property */
  if (cs_equation_param_has_diffusion(eqp))
    if (!(eqb->diff_pty_uniform))
      cs_equation_set_diffusion_property_cw(eqp, cm, t_eval, cell_flag, cb);

  /* Set the (linear) reaction property */
  if (cs_equation_param_has_reaction(eqp)) {

    /* Define the local reaction property */
    cb->rpty_val = 0;
    for (int r = 0; r < eqp->n_reaction_terms; r++)
      if (eqb->reac_pty_uniform[r])
        cb->rpty_val += cb->rpty_vals[r];
      else
        cb->rpty_val += cs_property_value_in_cell(cm,
                                                  eqp->reaction_properties[r],
                                                  t_eval);

  }

  if (cs_equation_param_has_time(eqp))
    if (!(eqb->time_pty_uniform))
      cb->tpty_val = cs_property_value_in_cell(cm, eqp->time_property, t_eval);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
  if (cs_dbg_cw_test(cm)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion, advection,
 *          reaction terms in CDO-Vb schemes. If asked, a mass matrix is also
 *          computed and stored in cb->hdg
 *
 * \param[in]      time_eval   time at which analytic function are evaluated
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vb_advection_diffusion_reaction(double                         time_eval,
                                 const cs_equation_param_t     *eqp,
                                 const cs_equation_builder_t   *eqb,
                                 const cs_cdovb_scaleq_t       *eqc,
                                 const cs_cell_mesh_t          *cm,
                                 cs_face_mesh_t                *fm,
                                 cs_cell_sys_t                 *csys,
                                 cs_cell_builder_t             *cb)
{
  if (cs_equation_param_has_diffusion(eqp)) {   /* DIFFUSION TERM
                                                 * ============== */

    /* Define the local stiffness matrix: local matrix owned by the cellwise
       builder (store in cb->loc) */
    eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

    /* Add the local diffusion operator to the local system */
    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(cm))
      cs_cell_sys_dump("\n>> Local system after diffusion", csys);
#endif
  }

  if (cs_equation_param_has_convection(eqp)) {  /* ADVECTION TERM
                                                 * ============== */

    /* Define the local advection matrix */
    eqc->get_advection_matrix(eqp, cm, time_eval, fm, cb);

    /* Add it to the local system */
    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(cm))
      cs_cell_sys_dump("\n>> Local system after advection", csys);
#endif
  }

  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                  * =========== */

    /* Build the mass matrix adn store it in cb->hdg */
    eqc->get_mass_matrix(eqc->hdg_mass, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(cm)) {
      cs_log_printf(CS_LOG_DEFAULT, ">> Local mass matrix");
      cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, cb->hdg);
    }
#endif
  }

  if (cs_equation_param_has_reaction(eqp)) { /* REACTION TERM
                                              * ============= */

    /* Update local system matrix with the reaction term
       cb->hdg corresponds to the current mass matrix */
    cs_sdm_add_mult(csys->mat, cb->rpty_val, cb->hdg);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(cm))
      cs_cell_sys_dump("\n>> Local system after reaction", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the local system in CDO-Vb schemes
 *
 * \param[in]      time_eval   time at which analytical function are evaluated
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vb_apply_bc(cs_real_t                      time_eval,
             const cs_equation_param_t     *eqp,
             const cs_cdovb_scaleq_t       *eqc,
             const cs_cell_mesh_t          *cm,
             cs_face_mesh_t                *fm,
             cs_cell_sys_t                 *csys,
             cs_cell_builder_t             *cb)
{
  if ((csys->cell_flag & CS_FLAG_BOUNDARY) == 0)
    return;

  /* Neumann boundary conditions */
  if (csys->has_nhmg_neumann) {
    for (short int v  = 0; v < cm->n_vc; v++)
      csys->rhs[v] += csys->neu_values[v];
  }

  /* Contribution for the advection term: csys is updated inside
     (matrix and rhs) and Dirichlet BCs are handled inside */
  if (cs_equation_param_has_convection(eqp))
    eqc->add_advection_bc(cm, eqp, time_eval, fm, cb, csys);

  /* The enforcement of the Dirichlet has to be done after all
     other contributions */
  if (cs_equation_param_has_diffusion(eqp)) {

    if (csys->has_dirichlet) /* csys is updated inside (matrix and rhs) */
      eqc->enforce_dirichlet(eqp->diffusion_hodge,
                             cm,
                             eqc->boundary_flux_op,
                             fm, cb, csys);

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
  if (cs_dbg_cw_test(cm))
    cs_cell_sys_dump("\n>> Local system after BC treatment", csys);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from a scalar-valued CDO-Vb scheme
 *
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      x0       initial guess for the linear system
 * \param[in]      rhs      pointer to a cs_mesh_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in]      xsol     pointer to an array storing the solution of
 *                          the linear system
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

static int
_solve_system(cs_sles_t                    *sles,
              const cs_matrix_t            *matrix,
              const cs_real_t              *x0,
              const cs_real_t              *rhs,
              const cs_equation_param_t    *eqp,
              cs_real_t                    *p_xsol[])
{
  cs_real_t  *x = NULL, *b = NULL;
  int  n_iters = 0;
  double  residual = DBL_MAX;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  cs_range_set_t  *rset = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];

  /* Prepare solving (handle parallelism) */
  cs_gnum_t  nnz = cs_equation_prepare_system(1,          /* stride */
                                              n_vertices, /* n_scatter_elts */
                                              x0,
                                              rhs,
                                              matrix,
                                              rset,
                                              &x, &b);

  /* Solve the linear solver */
  const double  r_norm = 1.0; /* No renormalization by default (TODO) */
  const cs_param_itsol_t  itsol_info = eqp->itsol_info;

  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    itsol_info.eps,
                                                    r_norm,
                                                    &n_iters,
                                                    &residual,
                                                    b,
                                                    x,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  if (eqp->sles_verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%s/sles_cvg> code %-d n_iters %d"
                  " residual % -8.4e nnz %lu\n",
                  eqp->name, code, n_iters, residual, nnz);

  if (cs_glob_n_ranks > 1) { /* Parallel mode */

    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         x,
                         x);

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
  cs_dbg_array_fprintf(NULL, "sol.log", 1e-16, n_vertices, x, 6);
  cs_dbg_array_fprintf(NULL, "rhs.log", 1e-16, n_vertices, rhs, 6);
#endif

  /* Free what can be freed at this stage */
  if (b != rhs)
    BFT_FREE(b);
  cs_sles_free(sles);

  /* Return pointer to the computed solution */
  *p_xsol = x;

  return n_iters;
}

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
cs_cdovb_scaleq_is_initialized(void)
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
cs_cdovb_scaleq_init_common(const cs_cdo_quantities_t    *quant,
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

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cdovb_cell_sys[t_id] = cs_cell_sys_create(connect->n_max_vbyc,
                                                 connect->n_max_fbyc,
                                                 1, NULL);
    cs_cdovb_cell_bld[t_id] = _cell_builder_create(connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cdovb_cell_sys[0] = cs_cell_sys_create(connect->n_max_vbyc,
                                            connect->n_max_fbyc,
                                            1, NULL);
  cs_cdovb_cell_bld[0] = _cell_builder_create(connect);

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
cs_cdovb_scaleq_get(cs_cell_sys_t       **csys,
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
cs_cdovb_scaleq_finalize_common(void)
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
  cs_cdovb_cell_sys = NULL;
  cs_cdovb_cell_bld = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_scaleq_t structure storing data useful
 *         for building and managing such a scheme
 *
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id      id of the variable field
 * \param[in]      bflux_id    id of the boundary flux field
 * \param[in, out] eqb         pointer to a \ref cs_equation_builder_t struct.
 *
 * \return a pointer to a new allocated cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_scaleq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOVB && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of equation.\n"
              " Expected: scalar-valued CDO vertex-based equation.", __func__);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_vertices = connect->n_vertices;

  cs_cdovb_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdovb_scaleq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  eqc->n_dofs = n_vertices;

  /* Flag to indicate the minimal set of quantities to build in a cell mesh
     According to the situation, additional flags have to be set */
  eqb->msh_flag = CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PE |
    CS_CDO_LOCAL_EV;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FE |
    CS_CDO_LOCAL_FEQ;

  /* DIFFUSION */
  eqc->get_stiffness_matrix = NULL;
  eqc->boundary_flux_op = NULL;
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
                (" %s: Invalid type of algorithm to build the diffusion term."),
                __func__);

    } /* Switch on Hodge algo. */

  } /* DIFFUSION */

  eqc->enforce_dirichlet = NULL;
  eqc->vtx_bc_flag = NULL;
  switch (eqp->enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_dirichlet;
    eqc->vtx_bc_flag = cs_equation_set_vertex_bc_flag(connect, eqb->face_bc);
    break;

  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_dirichlet;
    eqc->vtx_bc_flag = cs_equation_set_vertex_bc_flag(connect, eqb->face_bc);
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    eqb->bd_msh_flag |= CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_FEQ;
    eqc->enforce_dirichlet = cs_cdovb_diffusion_weak_dirichlet;
    if (cs_equation_param_has_diffusion(eqp) == false)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid choice of Dirichlet enforcement.\n"
                " Diffusion term should be active.", __func__);
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    eqb->bd_msh_flag |= CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_FEQ;
    eqc->enforce_dirichlet = cs_cdovb_diffusion_wsym_dirichlet;
    if (cs_equation_param_has_diffusion(eqp) == false)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid choice of Dirichlet enforcement.\n"
                " Diffusion term should be active.", __func__);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  /* ADVECTION */
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
      } /* Scheme */
      break; /* Formulation */

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
      } /* Scheme */
      break; /* Formulation */

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
      eqb->sys_flag |= CS_FLAG_SYS_SYM; /* Algebraic system is symmetric */

  }

  /* REACTION */
  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->reaction_hodge.algo == CS_PARAM_HODGE_ALGO_WBS) {
      eqb->msh_flag |= CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_FEQ |
        CS_CDO_LOCAL_HFQ;
      eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
    }
    else
      bft_error(__FILE__, __LINE__, 0,
                " Invalid choice of algorithm for the reaction term.");

  } /* Reaction */

  /* TIME */
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

  /* SOURCE TERM */
  eqc->source_terms = NULL;

  if (cs_equation_param_has_sourceterm(eqp)) {

    /* When the deprecated mode will be removed. This array needs to be
       allocated only if a theta scheme is used */
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
 * \brief  Destroy a cs_cdovb_scaleq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_scaleq_free_context(void   *builder)
{
  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)builder;

  if (eqc == NULL)
    return eqc;

  /* These arrays may have not been allocated */
  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->cell_values);
  BFT_FREE(eqc->vtx_bc_flag);

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
 * \param[in, out] data           pointer to cs_cdovb_scaleq_t structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_initialize_system(const cs_equation_param_t  *eqp,
                                  cs_equation_builder_t      *eqb,
                                  void                       *data,
                                  cs_matrix_t               **system_matrix,
                                  cs_real_t                 **system_rhs)
{
  CS_UNUSED(eqp);

  if (data == NULL)
    return;
  assert(*system_matrix == NULL && *system_rhs == NULL);

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)data;
  cs_timer_t  t0 = cs_timer_time();

  /* Create the matrix related to the current algebraic system */
  *system_matrix = cs_matrix_create(cs_shared_ms);

  /* Allocate and initialize the related right-hand side */
  BFT_MALLOC(*system_rhs, eqc->n_dofs, cs_real_t);
# pragma omp parallel for if  (eqc->n_dofs > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < eqc->n_dofs; i++) (*system_rhs)[i] = 0.0;

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings when the fields
 *         stem from a scalar CDO vertex-based scheme.
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in, out] eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval      time at which one evaluates BCs
 * \param[in, out] field_val   pointer to the values of the variable field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_set_dir_bc(const cs_mesh_t              *mesh,
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
                                   cs_cdovb_cell_bld[0], /* static variable */
                                   field_val);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_steady_state(double                      dt_cur,
                                   const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  CS_UNUSED(dt_cur);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;

  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Build an array storing the Dirichlet values at vertices and another one
     with a tags to detect vertices related to a Neumann BC (dt_cur is a dummy
     argument) */
  cs_real_t  *dir_values = NULL;
  short int  *neu_tags = NULL;

  _setup_bc(dt_cur, mesh, eqp, eqb, &dir_values, &neu_tags);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav,       \
         dir_values, neu_tags, fld, rs, cs_cdovb_cell_sys, cs_cdovb_cell_bld)
  {
    /* Set variables and structures inside the OMP section so that each thread
       has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_real_t  time_eval = dt_cur; /* dummy variable */

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
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

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _init_vb_cell_system(cell_flag, cm, eqp, eqb,
                           dir_values, neu_tags, fld->val, time_eval,
                           csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vb_advection_diffusion_reaction(time_eval,
                                       eqp, eqb, eqc, cm, fm, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */

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

      /* Apply boundary conditions */
      _vb_apply_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      cs_equation_assemble_matrix(csys, rs, mav); /* Matrix assembly */

      for (short int v = 0; v < cm->n_vc; v++) /* Assemble RHS */
#       pragma omp atomic
        rhs[cm->v_ids[v]] += csys->rhs[v];

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(neu_tags);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Now solve the system */
  cs_real_t  *x_sol = NULL;
  cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

  _solve_system(sles, matrix, fld->val, rhs, eqp, &x_sol);

  /* Update field */
  t0 = cs_timer_time();

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Overwrite the initial guess with the computed solution */
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    fld->val[i] = x_sol[i];

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);

  /* Free remaining buffers */
  BFT_FREE(x_sol);
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
  cs_sles_free(sles);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         Implicit time scheme is used to progress in time.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_implicit(double                      dt_cur,
                               const cs_mesh_t            *mesh,
                               const int                   field_id,
                               const cs_equation_param_t  *eqp,
                               cs_equation_builder_t      *eqb,
                               void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  time_eval = t_cur + dt_cur;

  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_IMPLICIT);

  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Build an array storing the Dirichlet values at vertices and another one
     with a tags to detect vertices related to a Neumann BC */
  cs_real_t  *dir_values = NULL;
  short int  *neu_tags = NULL;

  _setup_bc(t_cur + dt_cur, mesh, eqp, eqb, &dir_values, &neu_tags);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav,       \
         dir_values, neu_tags, fld, rs, cs_cdovb_cell_sys, cs_cdovb_cell_bld)
  {
    /* Set variables and structures inside the OMP section so that each thread
       has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
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

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _init_vb_cell_system(cell_flag, cm, eqp, eqb,
                           dir_values, neu_tags, fld->val, time_eval,
                           csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vb_advection_diffusion_reaction(time_eval,
                                       eqp, eqb, eqc, cm, fm, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */

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

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        /* |c|*wvc = |dual_cell(v) cap c| */
        assert(cs_flag_test(eqb->msh_flag, CS_CDO_LOCAL_PVQ));
        const double  ptyc = cb->tpty_val * cm->vol_c / dt_cur;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix */
        for (short int i = 0; i < cm->n_vc; i++) {

          const double  dval =  ptyc * cm->wvc[i];

          /* Update the RHS with values at time t_n */
          csys->rhs[i] += dval * csys->val_n[i];

          /* Add the diagonal contribution from time matrix */
          csys->mat->val[i*(cm->n_vc + 1)] += dval;

        }

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val/dt_cur;
        const cs_sdm_t  *mass_mat = cb->hdg;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix */

        /* Update rhs with csys->mat*p^n */
        double  *time_pn = cb->values;
        cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);
        for (short int i = 0; i < csys->n_dofs; i++)
          csys->rhs[i] += tpty_coef*time_pn[i];

        /* Update the cellwise system with the time matrix */
        cs_sdm_add_mult(csys->mat, tpty_coef, mass_mat);

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump("\n>> Local system after time", csys);
#endif

      /* Apply boundary conditions */
      _vb_apply_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      cs_equation_assemble_matrix(csys, rs, mav); /* Matrix assembly */

      for (short int v = 0; v < cm->n_vc; v++) /* Assemble RHS */
#       pragma omp atomic
        rhs[cm->v_ids[v]] += csys->rhs[v];

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(neu_tags);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Now solve the system */
  cs_real_t  *x_sol = NULL;
  cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

  _solve_system(sles, matrix, fld->val, rhs, eqp, &x_sol);

  /* Update field */
  t0 = cs_timer_time();

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Overwrite the initial guess with the computed solution */
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    fld->val[i] = x_sol[i];

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);

  /* Free remaining buffers */
  BFT_FREE(x_sol);
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
  cs_sles_free(sles);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-Vb scheme
 *         Theta time scheme is used to progress in time.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_solve_theta(double                      dt_cur,
                            const cs_mesh_t            *mesh,
                            const int                   field_id,
                            const cs_equation_param_t  *eqp,
                            cs_equation_builder_t      *eqb,
                            void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const double  tcoef = 1 - eqp->theta;

  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_THETA ||
         eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO);

  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Build an array storing the Dirichlet values at vertices
     and another one with a tags to detect vertices related to a
     Neumann BC */
  cs_real_t  *dir_values = NULL;
  short int  *neu_tags = NULL;

  _setup_bc(t_cur + dt_cur, mesh, eqp, eqb, &dir_values, &neu_tags);

  /* Initialize the local system: rhs */
  cs_real_t  *rhs = NULL;
  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Detect the first call (in this case, we compute the initial source term)*/
  bool  compute_initial_source = false;
  if (cs_shared_time_step->nt_cur == cs_shared_time_step->nt_prev) {
    compute_initial_source = true;
  }
  else { /* Add contribution of the previous computed source term */

    if (eqc->source_terms != NULL) {

      for (cs_lnum_t v = 0; v < n_vertices; v++)
        rhs[v] += tcoef * eqc->source_terms[v];
      memset(eqc->source_terms, 0, n_vertices * sizeof(cs_real_t));

      if (eqp->enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC ||
          eqp->enforcement == CS_PARAM_BC_ENFORCE_PENALIZED) {

        assert(eqc->vtx_bc_flag != NULL);
        for (cs_lnum_t v = 0; v < n_vertices; v++) {
          if (eqc->vtx_bc_flag[v] & CS_CDO_BC_HMG_DIRICHLET)
            rhs[v] = 0.;
          else if (eqc->vtx_bc_flag[v] & CS_CDO_BC_DIRICHLET)
            rhs[v] = 0.;
        }

      } /* Algebraic or penalized enforcement is set */
    } /* At least one source term is defined */

  }

  /* Initialize the local system: matrix */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav, dir_values, \
         neu_tags, fld, rs, cs_cdovb_cell_sys, cs_cdovb_cell_bld,       \
         compute_initial_source)
  {
    /* Set variables and structures inside the OMP section so that each thread
       has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_real_t  time_eval = t_cur + 0.5*dt_cur;

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
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

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _init_vb_cell_system(cell_flag, cm, eqp, eqb,
                           dir_values, neu_tags, fld->val, time_eval,
                           csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vb_advection_diffusion_reaction(time_eval,
                                       eqp, eqb, eqc, cm, fm, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */
        if (compute_initial_source) {

          /* Reset the local contribution */
          memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

          cs_source_term_compute_cellwise(eqp->n_source_terms,
                      (cs_xdef_t *const *)eqp->source_terms,
                                          cm,
                                          eqb->source_mask,
                                          eqb->compute_source,
                                          t_cur,
                                          NULL,  /* No input structure */
                                          cb,    /* mass matrix is cb->hdg */
                                          csys->source);

          for (short int v = 0; v < cm->n_vc; v++)
            csys->rhs[v] += tcoef * csys->source[v];

        }

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
                                        t_cur + dt_cur,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        for (short int v = 0; v < cm->n_vc; v++)
          csys->rhs[v] += eqp->theta * csys->source[v];

      } /* End of term source */

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      /* STEP.1 >> Compute the contribution of the "adr" to the RHS:
       *           tcoef*adr_pn where adr_pn = csys->mat * p_n */
      double  *adr_pn = cb->values;
      cs_sdm_square_matvec(csys->mat, csys->val_n, adr_pn);
      for (short int i = 0; i < csys->n_dofs; i++) /* n_dofs = n_vc */
        csys->rhs[i] -= tcoef * adr_pn[i];

      /* STEP.2 >> Multiply csys->mat by theta */
      for (int i = 0; i < csys->n_dofs*csys->n_dofs; i++)
        csys->mat->val[i] *= eqp->theta;

      /* STEP.3 >> Handle the mass matrix
       * Two contributions for the mass matrix
       *  a) add to csys->mat
       *  b) add to rhs mass_mat * p_n */
      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        /* |c|*wvc = |dual_cell(v) cap c| */
        assert(cs_flag_test(eqb->msh_flag, CS_CDO_LOCAL_PVQ));
        const double  ptyc = cb->tpty_val * cm->vol_c / dt_cur;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
         *       >> Update the cellwise system with the time matrix */
        for (short int i = 0; i < cm->n_vc; i++) {

          const double  dval = ptyc * cm->wvc[i];

          /* Update the RHS with mass_mat * values at time t_n */
          csys->rhs[i] += dval * csys->val_n[i];

          /* Add the diagonal contribution from time matrix to the local
             system */
          csys->mat->val[i*(cm->n_vc + 1)] += dval;

        }

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val / dt_cur;
        const cs_sdm_t  *mass_mat = cb->hdg;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix */

        /* Update rhs with mass_mat*p^n */
        double  *time_pn = cb->values;
        cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);
        for (short int i = 0; i < csys->n_dofs; i++)
          csys->rhs[i] += tpty_coef*time_pn[i];

        /* Update the cellwise system with the time matrix */
        cs_sdm_add_mult(csys->mat, tpty_coef, mass_mat);

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump("\n>> Local system after time", csys);
#endif

      /* Apply boundary conditions */
      _vb_apply_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      cs_equation_assemble_matrix(csys, rs, mav); /* Matrix assembly */

      for (short int v = 0; v < cm->n_vc; v++) /* Assemble RHS */
#       pragma omp atomic
        rhs[cm->v_ids[v]] += csys->rhs[v];

      if (eqc->source_terms != NULL) {
        for (short int v = 0; v < cm->n_vc; v++) /* Assemble source terms */
#         pragma omp atomic
          eqc->source_terms[cm->v_ids[v]] += csys->source[v];
      }

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(neu_tags);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Now solve the system */
  cs_real_t  *x_sol = NULL;
  cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

  _solve_system(sles, matrix, fld->val, rhs, eqp, &x_sol);

  /* Update field */
  t0 = cs_timer_time();

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Overwrite the initial guess with the computed solution */
# pragma omp parallel for if (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++)
    fld->val[i] = x_sol[i];

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);

  /* Free remaining buffers */
  BFT_FREE(x_sol);
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
  cs_sles_free(sles);
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
 * \param[in, out] context    pointer to cs_cdovb_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_build_system(const cs_mesh_t            *mesh,
                             const cs_real_t            *field_val,
                             double                      dt_cur,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *context,
                             cs_real_t                  *rhs,
                             cs_matrix_t                *matrix)
{
  /* Sanity checks */
  assert(rhs != NULL && matrix != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;

  cs_log_printf(CS_LOG_DEFAULT,
                " %s: Deprecated mode for building the system.\n", __func__);
  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;

  /* Compute the values of the Dirichlet BC */
  cs_real_t  *dir_values = NULL;
  BFT_MALLOC(dir_values, quant->n_vertices, cs_real_t);
  memset(dir_values, 0, quant->n_vertices*sizeof(cs_real_t));

  cs_cdovb_scaleq_set_dir_bc(mesh, eqp, eqb, t_cur + dt_cur, dir_values);

  /* Tag faces with a non-homogeneous Neumann BC */
  short int  *neu_tags = cs_equation_tag_neumann_face(quant, eqp);

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, rhs, matrix, mav,       \
         dir_values, neu_tags, field_val,                               \
         cs_cdovb_cell_sys, cs_cdovb_cell_bld)
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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
      if (cs_dbg_cw_test(cm)) cs_cell_mesh_dump(cm);
#endif

      /* DIFFUSION TERM */
      /* ============== */

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        /* local matrix owned by the cellwise builder (store in cb->loc) */
        eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        /* Add the local diffusion operator to the local system */
        cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
        if (cs_dbg_cw_test(cm))
          cs_cell_sys_dump("\n>> Local system after diffusion", csys);
#endif
      } /* END OF DIFFUSION */

      /* ADVECTION TERM */
      /* ============== */

      if (cs_equation_param_has_convection(eqp)) {

        /* Define the local advection matrix */
        eqc->get_advection_matrix(eqp, cm, time_eval, fm, cb);

        cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
        if (cs_dbg_cw_test(cm))
          cs_cell_sys_dump("\n>> Local system after advection", csys);
#endif
      } /* END OF ADVECTION */

      /* MASS MATRIX */
      /* =========== */

      if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) {
        eqc->get_mass_matrix(eqc->hdg_mass, cm, cb); /* stored in cb->hdg */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
        if (cs_dbg_cw_test(cm)) {
          cs_log_printf(CS_LOG_DEFAULT, ">> Local mass matrix");
          cs_sdm_dump(c_id, csys->dof_ids, csys->dof_ids, cb->hdg);
        }
#endif
      }

      /* REACTION TERM */
      /* ============= */

      if (cs_equation_param_has_reaction(eqp)) {

        /* Update local system matrix with the reaction term
           cb->hdg corresponds to the current mass matrix */
        cs_sdm_add_mult(csys->mat, cb->rpty_val, cb->hdg);

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

      /* UNSTEADY TERM + TIME SCHEME */
      /* =========================== */

      if (cs_equation_param_has_time(eqp)) {

        /* Get the value of the time property */
        double  tpty_val = 1/dt_cur;
        if (eqb->time_pty_uniform)
          tpty_val *= cb->tpty_val;
        else
          tpty_val *= cs_property_value_in_cell(cm,
                                                eqp->time_property,
                                                time_eval);

        cs_sdm_t  *mass_mat = cb->hdg;
        if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) {

          assert(cs_flag_test(eqb->msh_flag, CS_CDO_LOCAL_PVQ));
          /* Switch to cb->loc. Used as a diagonal only */
          mass_mat = cb->loc;

          /* |c|*wvc = |dual_cell(v) cap c| */
          const double  ptyc = tpty_val * cm->vol_c;
          for (short int v = 0; v < cm->n_vc; v++)
            mass_mat->val[v] = ptyc * cm->wvc[v];

        }

        /* Apply the time discretization to the local system.
           Update csys (matrix and rhs) */
        eqc->apply_time_scheme(eqp, tpty_val, mass_mat, eqb->sys_flag, cb,
                               csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
        if (cs_dbg_cw_test(cm))
          cs_cell_sys_dump("\n>> Local system after time", csys);
#endif
      } /* END OF TIME CONTRIBUTION */

      /* BOUNDARY CONDITIONS */
      /* =================== */

      if (cell_flag & CS_FLAG_BOUNDARY) {

        /* Neumann boundary conditions */
        if (csys->has_nhmg_neumann) {
          for (short int v  = 0; v < cm->n_vc; v++)
            csys->rhs[v] += csys->neu_values[v];
        }

        /* Contribution for the advection term: csys is updated inside
           (matrix and rhs) and Dirichlet BCs are handled inside */
        if (cs_equation_param_has_convection(eqp)) {

          eqc->add_advection_bc(cm, eqp, time_eval, fm, cb, csys);

        }

        if (cs_equation_param_has_diffusion(eqp)) {

          /* The enforcement of the Dirichlet has to be done after all
             other contributions */
          if (csys->has_dirichlet) /* csys is updated inside (matrix and rhs) */
            eqc->enforce_dirichlet(eqp->diffusion_hodge,
                                   cm,
                                   eqc->boundary_flux_op,
                                   fm, cb, csys);

        }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
        if (cs_dbg_cw_test(cm))
          cs_cell_sys_dump("\n>> Local system after BC treatment", csys);
#endif
      } /* END OF BOUNDARY CONDITIONS */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(cm))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
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

      if (eqc->source_terms != NULL) {
        for (short int v = 0; v < cm->n_vc; v++)
#         pragma omp atomic
          eqc->source_terms[cm->v_ids[v]] += csys->source[v];
      }

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); // optional

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(neu_tags);
  cs_matrix_assembler_values_finalize(&mav);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
  if (eqc->source_terms != NULL)
    cs_dbg_darray_to_listing("EQ.BUILD >> TS", eqc->n_dofs, eqc->source_terms,
                             8);
#endif

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
 * \param[in, out] data       pointer to cs_cdovb_scaleq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_update_field(const cs_real_t            *solu,
                             const cs_real_t            *rhs,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data,
                             cs_real_t                  *field_val)
{
  CS_UNUSED(rhs);
  CS_UNUSED(eqp);

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t  *)data;
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
cs_cdovb_scaleq_get_vertex_values(void      *context)
{
  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
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
cs_cdovb_scaleq_get_cell_values(void      *context)
{
  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Reset buffer of values */
  if (eqc->cell_values == NULL)
    BFT_MALLOC(eqc->cell_values, quant->n_cells, cs_real_t);
  memset(eqc->cell_values, 0, quant->n_cells*sizeof(cs_real_t));

  /* Compute the values at cell centers from an interpolation of the field
     values defined at vertices */
  cs_reco_pv_at_cell_centers(connect->c2v, quant, pot->val,
                             eqc->cell_values);

  return eqc->cell_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain between time t_cur and t_cur + dt_cur
 *         Case of scalar-valued CDO vertex-based scheme
 *
 * \param[in]      eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb      pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context  pointer to a scheme builder structure
 * \param[in]      dt_cur   current value of the time step
 *
 * \return a pointer to a \ref cs_equation_balance_t structure
 */
/*----------------------------------------------------------------------------*/

cs_equation_balance_t *
cs_cdovb_scaleq_balance(const cs_equation_param_t     *eqp,
                        cs_equation_builder_t         *eqb,
                        void                          *context,
                        cs_real_t                      dt_cur)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  time_eval = t_cur + 0.5*dt_cur;

  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);
  cs_field_t  *bflux = cs_field_by_id(eqc->bflux_field_id);

  /* Assign the boundary flux for faces where Neumann is defined */
  cs_equation_init_boundary_flux_from_bc(t_cur, quant, eqp, bflux->val);

  /* Allocate and initialize the structure storing the balance evaluation */
  cs_equation_balance_t  *eb = cs_equation_balance_create(cs_flag_primal_vtx,
                                                          quant->n_vertices);

  cs_equation_balance_reset(eb);

  /* OpenMP block */
#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(dt_cur, quant, connect, eqp, eqb, eqc, pot, bflux,             \
         eb, cs_cdovb_cell_bld)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Set inside the OMP section so that each thread has its own value
     * Each thread get back its related structures:
     * Get the cell-wise view of the mesh and the algebraic system */
    cs_face_mesh_light_t  *fml = cs_cdo_local_get_face_mesh_light(t_id);
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = cs_cdovb_cell_bld[t_id];

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

      /* Set the value of the current potential */
      for (short int v = 0; v < cm->n_vc; v++)
        p_cur[v] = pot->val[cm->v_ids[v]];

      if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX)
        eqc->get_mass_matrix(eqc->hdg_mass, cm, cb); /* stored in cb->hdg */

      /* UNSTEADY TERM */
      if (cs_equation_param_has_time(eqp)) {

        /* Set the value of the previous potential */
        for (short int v = 0; v < cm->n_vc; v++)
          p_prev[v] = pot->val_pre[cm->v_ids[v]];

        if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) {

          assert(cs_flag_test(eqb->msh_flag, CS_CDO_LOCAL_PVQ));
          /* |c|*wvc = |dual_cell(v) cap c| */
          const double  ptyc = cb->tpty_val * cm->vol_c / dt_cur;
          for (short int v = 0; v < cm->n_vc; v++) {
            cs_real_t  dp = p_cur[v] - p_prev[v];
#           pragma omp atomic
            eb->unsteady_term[cm->v_ids[v]] += ptyc * cm->wvc[v] * dp;
          }

        }
        else {

          const double  ptyc = cb->tpty_val/dt_cur;
          cs_real_t  *dp = cb->values;
          cs_real_t  *res = cb->values + cm->n_vc;
          for (short int v = 0; v < cm->n_vc; v++) {
            res[v] = 0.;
            dp[v] = p_cur[v] - p_prev[v];
          }
          cs_sdm_square_matvec(cb->hdg, dp, res);

          for (short int v = 0; v < cm->n_vc; v++) {
#           pragma omp atomic
            eb->unsteady_term[cm->v_ids[v]] += ptyc*res[v];
          }

        } /* Add unsteady contribution */

      } /* TIME */

      /* Set p_theta */
      switch (eqp->time_scheme) {
      case CS_TIME_SCHEME_EXPLICIT:
        for (short int v = 0; v < cm->n_vc; v++)
          p_theta[v] = p_prev[v];
        break;

      case CS_TIME_SCHEME_CRANKNICO:
        for (short int v = 0; v < cm->n_vc; v++)
          p_theta[v] = 0.5*(p_cur[v] + p_prev[v]);
        break;

      case CS_TIME_SCHEME_THETA:
        for (short int v = 0; v < cm->n_vc; v++)
          p_theta[v] = eqp->theta*p_cur[v] + (1-eqp->theta)*p_prev[v];
        break;

      default:
        for (short int v = 0; v < cm->n_vc; v++)
          p_theta[v] = p_cur[v];
        break;

      } /* Switch on time scheme */

      /* REACTION TERM */
      if (cs_equation_param_has_reaction(eqp)) {

        /* Define the local reaction property */
        double  rpty_val = 0;
        for (int r = 0; r < eqp->n_reaction_terms; r++)
          if (eqb->reac_pty_uniform[r])
            rpty_val += cb->rpty_vals[r];
          else
            rpty_val += cs_property_value_in_cell(cm,
                                                  eqp->reaction_properties[r],
                                                  time_eval);

        cs_real_t  *res = cb->values;
        for (short int v = 0; v < cm->n_vc; v++)
          res[v] = 0.;

        cs_sdm_square_matvec(cb->hdg, p_theta, res);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->reaction_term[cm->v_ids[v]] += rpty_val * res[v];
        }

      } /* REACTION */

      /* DIFFUSION */
      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        if (!(eqb->diff_pty_uniform))
          cs_equation_set_diffusion_property_cw(eqp, cm, time_eval, cell_flag,
                                                cb);

        /* local matrix owned by the cellwise builder (store in cb->loc) */
        eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

        cs_real_t  *res = cb->values;
        for (short int v = 0; v < cm->n_vc; v++)
          res[v] = 0.;

        cs_sdm_square_matvec(cb->loc, p_theta, res);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->diffusion_term[cm->v_ids[v]] += res[v];
        }

      } /* DIFFUSION */

      /* ADVECTION TERM */
      if (cs_equation_param_has_convection(eqp)) {

        /* Define the local advection matrix */
        eqc->get_advection_matrix(eqp, cm, time_eval, fm, cb);

        cs_real_t  *res = cb->values;
        for (short int v = 0; v < cm->n_vc; v++)
          res[v] = 0.;

        cs_sdm_square_matvec(cb->loc, p_theta, res);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->advection_term[cm->v_ids[v]] += res[v];
        }

      } /* END OF ADVECTION */

      /* SOURCE TERM */
      if (cs_equation_param_has_sourceterm(eqp)) {

        cs_real_t  *src = cb->values;
        memset(src, 0, cm->n_vc*sizeof(cs_real_t));

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
                                        src);

        for (short int v = 0; v < cm->n_vc; v++) {
#         pragma omp atomic
          eb->source_term[cm->v_ids[v]] += src[v];
        }

      } /* End of term source */

      /* BOUNDARY CONDITIONS */
      if (cell_flag &  CS_FLAG_BOUNDARY) {

        const cs_cdo_bc_t  *face_bc = eqb->face_bc;

        /* Identify which face is a boundary face */
        for (short int f = 0; f < cm->n_fc; f++) {
          const cs_lnum_t  bf_id = cm->f_ids[f] - quant->n_i_faces;
          if (bf_id > -1) { /* Border face */

            /* Advective flux */
            if (cs_equation_param_has_convection(eqp)) {
              cs_advection_field_get_f2v_boundary_flux(cm,
                                                       eqp->adv_field,
                                                       f,
                                                       time_eval,
                                                       cb->values);

              for (short int v = 0; v < cm->n_vc; v++) {
                const cs_real_t  adv_flux = cb->values[v] * p_cur[v];
#               pragma omp atomic
                eb->boundary_term[bf_id] += adv_flux;
#               pragma omp atomic
                eb->advection_term[cm->v_ids[v]] += adv_flux;
              }
            }

            /* Diffusive flux */
            if (cs_equation_param_has_diffusion(eqp)) {
              if (face_bc->flag[bf_id] & CS_CDO_BC_DIRICHLET ||
                  face_bc->flag[bf_id] & CS_CDO_BC_HMG_DIRICHLET) {
                cs_cdovb_diffusion_face_p0_flux(cm,
                        (const cs_real_t (*)[3])cb->dpty_mat,
                                                p_cur,
                                                f,
                                                time_eval,
                                                cb->values);

                for (short int v = 0; v < cm->n_vc; v++) {
#                 pragma omp atomic
                  eb->boundary_term[bf_id] += cb->values[v];
#                 pragma omp atomic
                  eb->diffusion_term[cm->v_ids[v]] += cb->values[v];
                }

              }
            }

          } /* Is a boundary face */
        } /* Loop on cell faces */

      } /* BOUNDARY CONDITIONS */

    } /* Main loop on cells */

    if (p_cur != _p_cur) {
      BFT_FREE(p_cur);
      BFT_FREE(p_prev);
      BFT_FREE(p_theta);
    }

  } /* OpenMP Block */

  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++)
    eb->balance[v_id] =
      eb->unsteady_term[v_id] + eb->reaction_term[v_id] +
      eb->diffusion_term[v_id] + eb->advection_term[v_id] +
      eb->source_term[v_id];

  /* Parallel synchronisation */
  cs_equation_balance_sync(connect, eb);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);

  return eb;
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
cs_cdovb_scaleq_compute_flux_across_plane(const cs_real_t             normal[],
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
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_adjacency_t  *f2c = connect->f2c;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;

  double  pf;
  cs_real_3_t  gc, pty_gc;
  cs_real_33_t  pty_tens;

  if (ml_t == CS_MESH_LOCATION_BOUNDARY_FACES) { // Belongs to only one cell

    const cs_lnum_t  n_i_faces = connect->n_faces[2];
    const cs_lnum_t  *cell_ids = f2c->ids + f2c->idx[n_i_faces];

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  bf_id = (elt_ids != NULL) ? elt_ids[i] : i;
      const cs_lnum_t  f_id = n_i_faces + bf_id;
      const cs_lnum_t  c_id = cell_ids[bf_id];
      const cs_quant_t  f = cs_quant_set_face(f_id, quant);
      const short int  sgn = (_dp3(f.unitv, normal) < 0) ? -1 : 1;
      const double  coef = sgn * f.meas;

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Compute the local diffusive flux */
        cs_reco_grad_cell_from_pv(c_id, connect, quant, pdi, gc);
        cs_property_get_cell_tensor(c_id, t_cur,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);
        cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);

        /* Update the diffusive flux */
        *d_flux += -coef * _dp3(f.unitv, pty_gc);

      }

      if (cs_equation_param_has_convection(eqp)) {

        cs_nvec3_t  adv_c;

        /* Compute the local advective flux */
        cs_advection_field_get_cell_vector(c_id, eqp->adv_field, &adv_c);
        cs_reco_pf_from_pv(f_id, connect, quant, pdi, &pf);

        /* Update the convective flux */
        *c_flux += coef * adv_c.meas * _dp3(adv_c.unitv, f.unitv) * pf;

      }

    } /* Loop on selected border faces */

  }
  else if (ml_t == CS_MESH_LOCATION_INTERIOR_FACES) {

    if (n_elts[0] > 0 && elt_ids == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _(" Computing the flux across all interior faces is"
                  " not managed yet."));

    for (cs_lnum_t i = 0; i < n_elts[0]; i++) {

      const cs_lnum_t  f_id = elt_ids[i];
      const cs_quant_t  f = cs_quant_set_face(f_id, quant);
      const short int  sgn = (_dp3(f.unitv, normal) < 0) ? -1 : 1;

      for (cs_lnum_t j = f2c->idx[f_id]; j < f2c->idx[f_id+1]; j++) {

        const cs_lnum_t  c_id = f2c->ids[j];

        if (cs_equation_param_has_diffusion(eqp)) {

          const double  coef = 0.5 * sgn * f.meas; // mean value at the face

          /* Compute the local diffusive flux */
          cs_reco_grad_cell_from_pv(c_id, connect, quant, pdi, gc);
          cs_property_get_cell_tensor(c_id, t_cur,
                                      eqp->diffusion_property,
                                      eqp->diffusion_hodge.inv_pty,
                                      pty_tens);
          cs_math_33_3_product((const cs_real_t (*)[3])pty_tens, gc, pty_gc);
          *d_flux += -coef * _dp3(f.unitv, pty_gc);

        }

        if (cs_equation_param_has_convection(eqp)) {

          cs_nvec3_t  adv_c;

          /* Compute the local advective flux */
          cs_reco_pf_from_pv(f_id, connect, quant, pdi, &pf);

          /* Evaluate the advection field at the face */
          cs_advection_field_get_cell_vector(c_id, eqp->adv_field, &adv_c);

          /* Update the convective flux (upwinding according to adv_f) */
          const double  dpc = _dp3(adv_c.unitv, f.unitv);

          cs_real_t  fconv_flux = 0;
          if (dpc > 0) {
            if (f2c->sgn[j] > 0) // nf points outward c; adv.nf > 0
              fconv_flux = adv_c.meas * dpc * sgn * f.meas * pf;
          }
          else if (dpc < 0) {
            if (f2c->sgn[j] < 0) // nf points inward c; adv.nf < 0
              fconv_flux = adv_c.meas * dpc * sgn * f.meas * pf;
          }
          else  // centered approach
            fconv_flux = 0.5 * adv_c.meas * dpc * sgn * f.meas * pf;

          *c_flux += fconv_flux;

        }

      }

    }  /* Loop on selected interior faces */

  } // Set of interior or border faces

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
 * \param[in, out]  data        pointer to cs_cdovb_scaleq_t structure
 * \param[in, out]  location    where the flux is defined
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_scaleq_cellwise_diff_flux(const cs_real_t             *values,
                                   const cs_equation_param_t   *eqp,
                                   cs_real_t                    t_eval,
                                   cs_equation_builder_t       *eqb,
                                   void                        *data,
                                   cs_flag_t                    location,
                                   cs_real_t                   *diff_flux)
{
  cs_cdovb_scaleq_t  *eqc = (cs_cdovb_scaleq_t  *)data;

  /* Sanity checks */
  assert(diff_flux != NULL && eqp != NULL && eqc != NULL && eqb != NULL);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  if (!cs_flag_test(location, cs_flag_primal_cell) &&
      !cs_flag_test(location, cs_flag_dual_face_byc))
    bft_error(__FILE__, __LINE__, 0,
              " Incompatible location.\n"
              " Stop computing a cellwise diffusive flux.");

  /* If no diffusion, return after resetting */
  if (cs_equation_param_has_diffusion(eqp) == false) {

    size_t  size = 0;
    if (cs_flag_test(location, cs_flag_primal_cell))
      size = 3*quant->n_cells;
    else if (cs_flag_test(location, cs_flag_dual_face_byc))
      size = connect->c2e->idx[quant->n_cells];

#   pragma omp parallel for if (size > CS_THR_MIN)
    for (size_t i = 0; i < size; i++) diff_flux[i] = 0;

    return;
  }

  cs_timer_t  t0 = cs_timer_time();

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)          \
  shared(t_eval, quant, connect, location, eqp, eqb, eqc, diff_flux, values, \
         cs_cdovb_cell_bld)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif
    double  *pot = NULL;

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = cs_cdovb_cell_bld[t_id];
    cs_flag_t  msh_flag = CS_CDO_LOCAL_PV;
    cs_hodge_t  *get_diffusion_hodge = NULL;
    cs_cdo_cellwise_diffusion_flux_t  *compute_flux = NULL;

#if defined(DEBUG) && !defined(NDEBUG)
    cs_cell_mesh_reset(cm);
#endif

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      BFT_MALLOC(pot, connect->n_max_vbyc, double);

      msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_EV |
        CS_CDO_LOCAL_PVQ;

      /* Set function pointers */
      if (cs_flag_test(location, cs_flag_primal_cell)) {
        msh_flag |= CS_CDO_LOCAL_EV;
        compute_flux = cs_cdo_diffusion_vcost_get_pc_flux;
      }
      else if (cs_flag_test(location, cs_flag_dual_face_byc)) {
        get_diffusion_hodge = cs_hodge_epfd_cost_get;
        compute_flux = cs_cdo_diffusion_vcost_get_dfbyc_flux;
      }

      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      BFT_MALLOC(pot, connect->n_max_vbyc, double);

      /* Set function pointers */
      get_diffusion_hodge = cs_hodge_epfd_voro_get;
      if (cs_flag_test(location, cs_flag_primal_cell))
        compute_flux = cs_cdo_diffusion_vcost_get_pc_flux;
      else if (cs_flag_test(location, cs_flag_dual_face_byc))
        compute_flux = cs_cdo_diffusion_vcost_get_dfbyc_flux;

      msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_EV |
        CS_CDO_LOCAL_EFQ | CS_CDO_LOCAL_PVQ;
      break;

    case CS_PARAM_HODGE_ALGO_WBS:
      BFT_MALLOC(pot, connect->n_max_vbyc + 1, double);

      msh_flag |= CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PEQ |
        CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_FEQ |
        CS_CDO_LOCAL_EV;

      /* Set function pointers */
      if (cs_flag_test(location, cs_flag_primal_cell)) {
        compute_flux = cs_cdo_diffusion_wbs_get_pc_flux;
        msh_flag |= CS_CDO_LOCAL_HFQ;
      }
      else if (cs_flag_test(location, cs_flag_dual_face_byc)) {
        compute_flux = cs_cdo_diffusion_wbs_get_dfbyc_flux;
        msh_flag |= CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_EFQ;
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0, " Invalid Hodge algorithm");

    } /* Switch hodge algo. */

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_equation_set_diffusion_property(eqp, 0, t_eval, 0, cb);

    /* Define the flux by cellwise contributions */
#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(cm)) cs_cell_mesh_dump(cm);
#endif

      if (!eqb->diff_pty_uniform) {
        cs_property_tensor_in_cell(cm,
                                   eqp->diffusion_property,
                                   t_eval,
                                   eqp->diffusion_hodge.inv_pty,
                                   cb->dpty_mat);
        if (eqp->diffusion_hodge.is_iso)
          cb->dpty_val = cb->dpty_mat[0][0];
      }

      /* Build the local dense matrix related to this operator
         (store in cb->hdg) */
      switch (eqp->diffusion_hodge.algo) {

      case CS_PARAM_HODGE_ALGO_COST:
      case CS_PARAM_HODGE_ALGO_VORONOI:

        if (cs_flag_test(location, cs_flag_dual_face_byc))
          get_diffusion_hodge(eqp->diffusion_hodge, cm, cb);

        /* Define a local buffer keeping the value of the discrete potential
           for the current cell */
        for (short int v = 0; v < cm->n_vc; v++)
          pot[v] = values[cm->v_ids[v]];

        break;

      case CS_PARAM_HODGE_ALGO_WBS:

        /* Define a local buffer keeping the value of the discrete potential
           for the current cell */
        pot[cm->n_vc] = 0.;
        for (short int v = 0; v < cm->n_vc; v++) {
          pot[v] = values[cm->v_ids[v]];
          pot[cm->n_vc] += cm->wvc[v]*pot[v];
        }
        break;

      default:
        bft_error(__FILE__, __LINE__, 0, " Invalid Hodge algorithm");

      } // End of switch

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
cs_cdovb_scaleq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *context)
{
  CS_UNUSED(field);
  CS_UNUSED(context);

  const cs_timer_t  t0 = cs_timer_time();

  if (cs_equation_param_has_convection(eqp)) {
    if (eqp->process_flag & CS_EQUATION_POST_UPWIND_COEF) {

      cs_real_t  *work_c = cs_equation_get_tmpbuf();
      char *postlabel = NULL;
      int  len = strlen(eqname) + 8 + 1;

      BFT_MALLOC(postlabel, len, char);
      sprintf(postlabel, "%s.UpwCoef", eqname);

      /* Compute in each cell an evaluation of upwind weight value */
      cs_cdo_advection_get_upwind_coef_cell(cs_shared_quant,
                                            eqp->adv_scheme,
                                            work_c);

      cs_post_write_var(CS_POST_MESH_VOLUME,
                        CS_POST_WRITER_ALL_ASSOCIATED,
                        postlabel,
                        1,
                        true,                 // interlace
                        true,                 // true = original mesh
                        CS_POST_TYPE_cs_real_t,
                        work_c,               // values on cells
                        NULL,                 // values at internal faces
                        NULL,                 // values at border faces
                        cs_shared_time_step); // time step management struct.

      BFT_FREE(postlabel);

    }
  } // Post a Peclet attached to cells

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
