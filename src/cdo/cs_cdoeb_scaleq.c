/*============================================================================
 * Build an algebraic CDO edge-based system. Degrees of freedom are defined as
 * a circulation. Degrees of freedom are scalar-valued.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdo_diffusion.h"
#include "cs_evaluate.h"
#include "cs_reco.h"

#include "cs_cdoeb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOEB_SCALEQ_DBG      0

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */
static cs_cell_sys_t      **cs_cdoeb_cell_system = NULL;
static cs_cell_builder_t  **cs_cdoeb_cell_builder = NULL;

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
_ebs_create_cell_builder(const cs_cdo_connect_t   *connect)
{
  const int  n_fc = connect->n_max_fbyc;
  const int  n_ec = connect->n_max_ebyc;
  const int  n_max = (n_fc > n_ec) ? n_fc : n_ec;

  cs_cell_builder_t  *cb = cs_cell_builder_create();

  BFT_MALLOC(cb->ids, n_max, int);
  memset(cb->ids, 0, n_max*sizeof(int));

  int  size = n_max*(n_max+1);
  size = CS_MAX(7*n_max, size);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(cs_real_t));

  size = 2*n_max;
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of
     operators */
  cb->hdg = cs_sdm_square_create(n_max);
  cb->aux = cs_sdm_square_create(n_max);
  cb->loc = cs_sdm_square_create(n_ec);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in]      cell_flag   flag related to the current cell
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         pointer to a cs_cdoeb_scaleq_t structure
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_eb_init_cell_system(cs_real_t                            t_eval,
                     const cs_flag_t                      cell_flag,
                     const cs_cell_mesh_t                *cm,
                     const cs_equation_param_t           *eqp,
                     const cs_equation_builder_t         *eqb,
                     const cs_cdoeb_scaleq_t             *eqc,
                     cs_cell_sys_t                       *csys,
                     cs_cell_builder_t                   *cb)
{
  /* Cell-wise view of the linear system to build */
  csys->c_id = cm->c_id;
  csys->cell_flag = cell_flag;
  csys->n_dofs = cm->n_ec;

  /* Initialize the local system */
  cs_cell_sys_reset(cm->n_fc, csys); /* Generic part */

  cs_sdm_square_init(csys->n_dofs, csys->mat);

  for (short int e = 0; e < cm->n_ec; e++) {
    csys->dof_ids[e] = cm->e_ids[e];
    csys->val_n[e] = eqc->edge_values[cm->e_ids[e]];
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Set the bc (specific part) */
    /* cs_equation_eb_set_cell_bc(cm, */
    /*                            eqp, */
    /*                            eqb->face_bc, */
    /*                            edge_bc_flag, */
    /*                            bc_values, */
    /*                            t_eval, */
    /*                            csys, */
    /*                            cb); */

  }

  /* Set the properties for this cell if not uniform */
  cs_equation_init_properties_cw(eqp, eqb, t_eval, cell_flag, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOEB_SCALEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion, advection,
 *          reaction terms in CDO-Eb schemes.
 *
 * \param[in]      time_eval   time at which analytic function are evaluated
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_eb_curlcurl(double                         time_eval,
             const cs_equation_param_t     *eqp,
             const cs_equation_builder_t   *eqb,
             const cs_cdoeb_scaleq_t       *eqc,
             const cs_cell_mesh_t          *cm,
             cs_cell_sys_t                 *csys,
             cs_cell_builder_t             *cb)
{
  CS_UNUSED(eqb);
  CS_UNUSED(time_eval);

  if (cs_equation_param_has_curlcurl(eqp)) {   /* CURL-CURL TERM
                                                * ============== */

    assert(cm->flag & CS_FLAG_COMP_FES);

    /* Define the local stiffness matrix: local matrix owned by the cellwise
       builder (store in cb->hdg) */
    eqc->get_curlcurl_hodge(eqp->curlcurl_hodge, cm, cb);

    /* Build the curl-curl operator in cb->loc */
    cs_sdm_square_init(cm->n_ec, cb->loc);

    for (int fk = 0; fk < cm->n_fc; fk++) {

      const cs_real_t  *h_row = cb->hdg->val + fk*cm->n_fc;

      for (int fl = 0; fl < cm->n_fc; fl++) {

        const cs_real_t  h_kl = h_row[fl];

        for (int ik = cm->f2e_idx[fk]; ik < cm->f2e_idx[fk+1]; ik++) {

          cs_real_t  *loc_row = cb->loc->val + cm->f2e_ids[ik]*cm->n_ec;
          const cs_real_t  ik_kl_coef = cm->f2e_sgn[ik] * h_kl;

          for (int il = cm->f2e_idx[fl]; il < cm->f2e_idx[fl+1]; il++) {

            loc_row[cm->f2e_ids[il]] += ik_kl_coef*cm->f2e_sgn[il];

          } /* Loop on face edges (il) */

        } /* Loop on face edges (ik) */

      } /* Loop on cell faces (l) */
    } /* Loop on cell faces (k)  */

    /* Add the local diffusion operator to the local system */
    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOEB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Local system after curlcurl", csys);
#endif
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Perform the assembly step
 *
 * \param[in]      eqc    context for this kind of discretization
 * \param[in]      cm     pointer to a cellwise view of the mesh
 * \param[in]      csys   pointer to a cellwise view of the system
 * \param[in]      rs     pointer to a cs_range_set_t structure
 * \param[in, out] eqa    pointer to a cs_equation_assemble_t structure
 * \param[in, out] mav    pointer to a cs_matrix_assembler_values_t structure
 * \param[in, out] rhs    right-hand side array
 */
/*----------------------------------------------------------------------------*/

inline static void
_assemble(const cs_cdoeb_scaleq_t           *eqc,
          const cs_cell_mesh_t              *cm,
          const cs_cell_sys_t               *csys,
          const cs_range_set_t              *rs,
          cs_equation_assemble_t            *eqa,
          cs_matrix_assembler_values_t      *mav,
          cs_real_t                         *rhs)
{
  /* Matrix assembly */
  eqc->assemble(csys, rs, eqa, mav);

  /* RHS assembly */
#if CS_CDO_OMP_SYNC_SECTIONS > 0
# pragma omp critical
  {
    for (int e = 0; e < cm->n_ec; e++)
      rhs[cm->e_ids[e]] += csys->rhs[e];
  }

  if (eqc->source_terms != NULL) { /* Source term */

#   pragma omp critical
    {
      for (int e = 0; e < cm->n_ec; e++) /* Source term assembly */
        eqc->source_terms[cm->e_ids[e]] += csys->source[e];
    }

  }
#else

  for (int e = 0; e < cm->n_ec; e++)
#   pragma omp atomic
    rhs[cm->e_ids[e]] += csys->rhs[e];

  if (eqc->source_terms != NULL) { /* Source term */
    for (int e = 0; e < cm->n_ec; e++) /* Source term assembly */
#     pragma omp atomic
      eqc->source_terms[cm->e_ids[e]] += csys->source[e];
  }
#endif
}

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Check if the generic structures for building a CDO-Eb scheme are
 *           allocated
 *
 * \return  true or false
 */
/*----------------------------------------------------------------------------*/

bool
cs_cdoeb_scaleq_is_initialized(void)
{
  if (cs_cdoeb_cell_builder == NULL || cs_cdoeb_cell_system == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Allocate work buffers and general structures related to CDO
 *           edge-based schemes. Set shared pointers.
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a time step structure
 * \param[in]  ms          pointer to a cs_matrix_structure_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_scaleq_init_common(const cs_cdo_quantities_t    *quant,
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

  BFT_MALLOC(cs_cdoeb_cell_system, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(cs_cdoeb_cell_builder, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_cdoeb_cell_system[i] = NULL;
    cs_cdoeb_cell_builder[i] = NULL;
  }

  /* TODO */
  const int  n_max_dofs = connect->n_max_ebyc;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_builder_t  *cb = _ebs_create_cell_builder(connect);
    cs_cdoeb_cell_builder[t_id] = cb;

    int  block_size = 3;
    cs_cdoeb_cell_system[t_id] = cs_cell_sys_create(n_max_dofs,
                                                    connect->n_max_fbyc,
                                                    1, NULL);
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_builder_t  *cb = _ebs_create_cell_builder(connect);
  cs_cdoeb_cell_builder[0] = cb;

  int  block_size = 3;
  cs_cdoeb_cell_system[0] =  cs_cell_sys_create(n_max_dofs,
                                            connect->n_max_fbyc,
                                            1, NULL);
#endif /* openMP */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *         in case of scalar-valued edge-based scheme
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_scaleq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  *csys = cs_cdoeb_cell_system[t_id];
  *cb = cs_cdoeb_cell_builder[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free work buffer and general structure related to CDO edge-based
 *         schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_scaleq_finalize_common(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(cs_cdoeb_cell_system[t_id]));
    cs_cell_builder_free(&(cs_cdoeb_cell_builder[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(cs_cdoeb_cell_system[0]));
  cs_cell_builder_free(&(cs_cdoeb_cell_builder[0]));
#endif /* openMP */

  BFT_FREE(cs_cdoeb_cell_system);
  BFT_FREE(cs_cdoeb_cell_builder);
  cs_cdoeb_cell_builder = NULL;
  cs_cdoeb_cell_system = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdoeb_scaleq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdoeb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdoeb_scaleq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOEB || eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of equation.\n"
              " Expected: scalar-valued CDO edge-based equation.", __func__);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_edges = connect->n_edges;

  cs_cdoeb_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdoeb_scaleq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  /* Dimensions of the algebraic system */
  eqc->n_dofs = n_edges;

  eqb->msh_flag = CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ | CS_FLAG_COMP_EF
    | CS_FLAG_COMP_FES;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_FLAG_COMP_EV | CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ;

  /* Values at each edge (interior and border) i.e. BCs are included */
  BFT_MALLOC(eqc->edge_values, n_edges, cs_real_t);
# pragma omp parallel for if (n_edges > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_edges; i++) eqc->edge_values[i] = 0;

  eqc->edge_values_pre = NULL;
  if (cs_equation_param_has_time(eqp)) {
    BFT_MALLOC(eqc->edge_values_pre, n_edges, cs_real_t);
# pragma omp parallel for if (n_edges > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < n_edges; i++) eqc->edge_values_pre[i] = 0;
  }

  if (cs_equation_param_has_curlcurl(eqp)) {

    eqb->msh_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ;

    switch (eqp->curlcurl_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      eqc->get_curlcurl_hodge = cs_hodge_fped_cost_get;;
      break;

    case CS_PARAM_HODGE_ALGO_BUBBLE:
      eqc->get_curlcurl_hodge = cs_hodge_fped_bubble_get;;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqc->get_curlcurl_hodge = cs_hodge_fped_voro_get;;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to build the diffusion term.",
                __func__);

    } /* Switch on Hodge algo. */

  } /* Diffusion */

  /* Essential boundary condition enforcement. The circulation along boundary
   * edges has the same behavior as enforcing a Dirichlet BC */
  eqc->enforce_essential_bc = NULL;
  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_essential_bc = cs_cdo_diffusion_alge_dirichlet;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  /* Source term */
  eqc->source_terms = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, n_edges, cs_real_t);
    memset(eqc->source_terms, 0, sizeof(cs_real_t)*n_edges);

  } /* There is at least one source term */

  /* Pre-defined a cs_hodge_builder_t structure */
  eqc->hdg_mass.is_unity = true;
  eqc->hdg_mass.is_iso   = true;
  eqc->hdg_mass.inv_pty  = false;
  eqc->hdg_mass.type = CS_PARAM_HODGE_TYPE_EPFD;
  eqc->hdg_mass.coef = 1.0; /* not useful in this case */

  if (eqp->do_lumping ||
      eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG ||
      eqb->sys_flag & CS_FLAG_SYS_REAC_DIAG) {

    eqc->hdg_mass.algo = CS_PARAM_HODGE_ALGO_VORONOI;
    eqc->get_mass_matrix = cs_hodge_epfd_voro_get;

  }
  else { /* WBS algorithm */

    eqc->hdg_mass.algo = CS_PARAM_HODGE_ALGO_WBS;
    eqc->get_mass_matrix = cs_hodge_vpcd_wbs_get;

  }

  /* Assembly process */
  eqc->assemble = cs_equation_assemble_set(CS_SPACE_SCHEME_CDOEB,
                                           CS_CDO_CONNECT_EDGE_SCAL);

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdoeb_scaleq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdoeb_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdoeb_scaleq_free_context(void   *builder)
{
  cs_cdoeb_scaleq_t  *eqc = (cs_cdoeb_scaleq_t *)builder;

  if (eqc == NULL)
    return eqc;

  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->edge_values);
  if (eqc->edge_values_pre != NULL)
    BFT_FREE(eqc->edge_values_pre);

  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-Eb schemes.
 *
 * \param[in]      t_eval     time at which one evaluates BCs
 * \param[in]      field_id   id related to the variable field of this equation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to the scheme context (cast on-the-fly)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_scaleq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdoeb_scaleq_t  *eqc = (cs_cdoeb_scaleq_t *)context;
  cs_real_t  *e_vals = eqc->edge_values;

  /* By default, 0 is set as initial condition for the computational domain */
  memset(e_vals, 0, quant->n_edges*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    /* Initialize values at mesh vertices */
    cs_flag_t  e_dof_flag = CS_FLAG_SCALAR | cs_flag_primal_edge;
    cs_lnum_t  *def2e_ids = (cs_lnum_t *)cs_equation_get_tmpbuf();
    cs_lnum_t  *def2e_idx = NULL;
    BFT_MALLOC(def2e_idx, eqp->n_ic_defs + 1, cs_lnum_t);

    cs_equation_sync_vol_def_at_edges(connect,
                                      eqp->n_ic_defs,
                                      eqp->ic_defs,
                                      def2e_idx,
                                      def2e_ids);

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_defs[def_id];
      const cs_lnum_t  n_e_selected = def2e_idx[def_id+1] - def2e_idx[def_id];
      const cs_lnum_t  *selected_lst = def2e_ids + def2e_idx[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_circulation_along_edges_by_value(def,
                                                     n_e_selected,
                                                     selected_lst,
                                                     e_vals);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_evaluate_circulation_along_edges_by_analytic(def,
                                                        t_eval,
                                                        n_e_selected,
                                                        selected_lst,
                                                        e_vals);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid way to initialize field values for eq. %s.\n",
                  __func__, eqp->name);

      } /* Switch on possible type of definition */

    } /* Loop on definitions */

  } /* Initial values to set */

  /* /\* Set the boundary values as initial values: Compute the values of the */
  /*    circulation where it is known thanks to the BCs *\/ */
  /* cs_equation_compute_circulation_eb(t_eval, */
  /*                                    mesh, */
  /*                                    quant, */
  /*                                    connect, */
  /*                                    eqp, */
  /*                                    eqb->face_bc, */
  /*                                    cs_cdoeb_cell_builder[0], */
  /*                                    e_vals); */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Eb scheme.
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdoeb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_scaleq_solve_steady_state(const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_EDGE_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_edges = quant->n_edges;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  cs_timer_t  t0 = cs_timer_time();

  cs_cdoeb_scaleq_t  *eqc = (cs_cdoeb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id); /* vector-valued cell-based */

  /* Build an array storing the values of the prescribed circulation at
     boundary */
  cs_real_t  *circ_bc_vals = NULL;

  BFT_MALLOC(circ_bc_vals, n_edges, cs_real_t);
  memset(circ_bc_vals, 0, n_edges*sizeof(cs_real_t));

  cs_equation_compute_circulation_eb(time_eval,
                                     mesh,
                                     quant,
                                     connect,
                                     eqp,
                                     circ_bc_vals);

  /* Initialize the local system: matrix and rhs */
  cs_real_t  res_normalization = 0.0;
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_edges, cs_real_t);
  memset(rhs, 0, n_edges*sizeof(cs_real_t));

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, circ_bc_vals, \
         fld, rs, cs_cdoeb_cell_system, cs_cdoeb_cell_builder,          \
         res_normalization)                                             \
  firstprivate(time_eval)
  {
    /* Set variables and structures inside the OMP section so that each thread
       has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdoeb_cell_system[t_id];
    cs_cell_builder_t  *cb = cs_cdoeb_cell_builder[t_id];
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:res_normalization)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _eb_init_cell_system(time_eval, cell_flag, cm, eqp, eqb, eqc,
                           csys, cb);

      /* Build and add the diffusion term to the local system. A mass matrix is
         also built if needed (stored it cb->hdg) */
      _eb_curlcurl(time_eval, eqp, eqb, eqc, cm, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */
        /* Reset the local contribution */
        memset(csys->source, 0, csys->n_dofs*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system */
        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        time_eval,
                                        NULL,  /* No input structure */
                                        cb,    /* mass matrix is cb->hdg */
                                        csys->source);

        /* Update the RHS */
        for (short int i = 0; i < csys->n_dofs; i++)
          csys->rhs[i] += csys->source[i];

      }

      /* Compute a norm of the RHS for the normalization of the residual
         of the linear system to solve */
      cs_equation_cw_scal_res_normalization(eqp->sles_param.resnorm_type,
                                            cm->vol_c, csys, NULL,
                                            &res_normalization);

      /* Boundary conditions */


#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOEB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY PROCESS
       * ================ */

      _assemble(eqc, cm, csys, rs, eqa, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  /* Free temporary buffers and structures */
  BFT_FREE(circ_bc_vals);

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Last step in the computation of the renormalization coefficient */
  cs_equation_sync_res_normalization(eqp->sles_param.resnorm_type,
                                     eqc->n_dofs,
                                     rhs,
                                     &res_normalization);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Solve the linear system */
  cs_equation_solve_scalar_system(eqc->n_dofs,
                                  eqp,
                                  matrix,
                                  rs,
                                  res_normalization,
                                  fld->val,
                                  rhs);

  /* Copy current field values to previous values and update the related
     field values */
  cs_field_current_to_previous(fld);

  /* Update the vector-valued field associated to cells */
  cs_reco_ccen_edge_dofs(connect, quant,
                         eqc->edge_values,
                         &(fld->val));

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_cdoeb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_scaleq_extra_op(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data)
{
  CS_UNUSED(eqname);
  CS_UNUSED(field);
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);
  CS_UNUSED(data);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at mesh edges (the DoFs)
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context     pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdoeb_scaleq_get_edge_values(void      *context)
{
  cs_cdoeb_scaleq_t  *eqc = (cs_cdoeb_scaleq_t *)context;

  return eqc->edge_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at mesh cells from a reconstruction of edge values.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context     pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdoeb_scaleq_get_cell_values(void      *context)
{
  cs_cdoeb_scaleq_t  *eqc = (cs_cdoeb_scaleq_t *)context;
  cs_field_t  *c_field = cs_field_by_id(eqc->var_field_id);

  return c_field->val;
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
cs_cdoeb_scaleq_read_restart(cs_restart_t    *restart,
                             const char      *eqname,
                             void            *scheme_context)
{

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
cs_cdoeb_scaleq_write_restart(cs_restart_t    *restart,
                              const char      *eqname,
                              void            *scheme_context)
{

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
