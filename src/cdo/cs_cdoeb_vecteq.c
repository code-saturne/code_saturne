/*============================================================================
 * Build an algebraic CDO edge-based system. Degrees of freedom are defined as
 * a circulation. Degrees of freedom are scalar-valued but the equation to
 * solve is vector-valued
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include "cs_cdo_toolbox.h"
#include "cs_cdo_solve.h"
#include "cs_evaluate.h"
#include "cs_reco.h"

#include "cs_cdoeb_vecteq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOEB_VECTEQ_DBG      0

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

  /* Local square dense matrices used during the construction of operators */

  cb->aux = cs_sdm_square_create(n_max);
  cb->loc = cs_sdm_square_create(n_ec);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      cm              pointer to a cellwise view of the mesh
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in]      eqb             pointer to a cs_equation_builder_t structure
 * \param[in]      eqc             pointer to a cs_cdoeb_vecteq_t structure
 * \param[in, out] csys            pointer to a cellwise view of the system
 * \param[in, out] cb              pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_eb_init_cell_system(const cs_cell_mesh_t                *cm,
                     const cs_equation_param_t           *eqp,
                     const cs_equation_builder_t         *eqb,
                     const cs_cdoeb_vecteq_t             *eqc,
                     cs_cell_sys_t                       *csys,
                     cs_cell_builder_t                   *cb)
{
  /* Cell-wise view of the linear system to build */

  csys->c_id = cm->c_id;
  csys->n_dofs = cm->n_ec;

  /* Initialize the local system */

  cs_cell_sys_reset(cm->n_fc, csys); /* Generic part */

  cs_sdm_square_init(csys->n_dofs, csys->mat);

  for (short int e = 0; e < cm->n_ec; e++) {
    csys->dof_ids[e] = cm->e_ids[e];
    csys->dof_flag[e] = 0;
    csys->val_n[e] = eqc->edge_values[cm->e_ids[e]];
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE)
    cs_equation_bc_set_cw_eb(cm,
                             eqp,
                             eqb->face_bc,
                             eqb->dir_values, /* circulation on edges */
                             csys,
                             cb);

  /* Special case to handle if enforcement by penalization or algebraic
   * This situation may happen with a tetrahedron with an edge lying on
   * the boundary (but no face)
   */

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_EDGE) {

    for (short int e = 0; e < cm->n_ec; e++) {

      csys->dof_flag[e] = eqc->edge_bc_flag[cm->e_ids[e]];
      if (cs_cdo_bc_is_circulation(csys->dof_flag[e])) {
        csys->has_dirichlet = true;
        csys->dir_values[e] = eqb->dir_values[cm->e_ids[e]];
      }

    }

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOEB_VECTEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion, advection,
 *          reaction terms in CDO-Eb schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] hodge       pointer to a hodge operator structure
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_eb_curlcurl(const cs_equation_param_t     *eqp,
             const cs_equation_builder_t   *eqb,
             const cs_cdoeb_vecteq_t       *eqc,
             const cs_cell_mesh_t          *cm,
             cs_hodge_t                    *hodge,
             cs_cell_sys_t                 *csys,
             cs_cell_builder_t             *cb)
{
  if (!cs_equation_param_has_curlcurl(eqp))
    return;

  /* Add the Curl-Curl term */

  assert(cm->flag & CS_FLAG_COMP_FES);

  /* Set the properties for this cell if not uniform */

  if (!eqb->curlcurl_pty_uniform)
    cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, cb->cell_flag, hodge);

  /* The first step is to build the hodge matrix associated to the curl-curl
     operator (stored in hodge->matrix) */

  eqc->get_curlcurl(cm, hodge, cb);

  /* Build the curl-curl operator in cb->loc */

  cs_sdm_square_init(cm->n_ec, cb->loc);

  for (int fk = 0; fk < cm->n_fc; fk++) {

    const cs_real_t  *h_row = hodge->matrix->val + fk*cm->n_fc;

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

  /* Add the local curl-curl operator to the local system */

  cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOEB_VECTEQ_DBG > 1
  if (cs_dbg_cw_test(eqp, cm, csys))
    cs_cell_sys_dump("\n>> Cell system after curlcurl", csys);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply boundary conditions. Only Dirichlet BCs which
 *         are enforced strongly.
 *         Update the local system after applying the time scheme.
 *         Case of CDO-Eb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] hodge       pointer to a cs_hodge_t structure
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_eb_enforce_values(const cs_equation_param_t     *eqp,
                   const cs_equation_builder_t   *eqb,
                   const cs_cdoeb_vecteq_t       *eqc,
                   const cs_cell_mesh_t          *cm,
                   cs_hodge_t                    *hodge,
                   cs_cell_sys_t                 *csys,
                   cs_cell_builder_t             *cb)
{
  if (cb->cell_flag > 0 && csys->has_dirichlet) {

    /* Boundary element (through either edges or faces) */

    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED) {

      /* csys is updated inside (matrix and rhs) */

      eqc->enforce_essential_bc(eqp, cm, NULL, hodge, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOEB_VECTEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after strong BC treatment", csys);
#endif
    }

  }

  if (cs_equation_param_has_internal_enforcement(eqp)) {

    /* Internal enforcement of DoFs: Update csys (matrix and rhs) */

    cs_equation_builder_enforce_dofs(eqb, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOEB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the residual normalization at the cellwise level according
 *         to the requested type of renormalization
 *         Case of CDO Edge-based vector-valued system.
 *
 * \param[in]  type        type of renormalization
 * \param[in]  cm          pointer to a cs_cell_mesh_t structure
 * \param[in]  csys        pointer to a cs_cell_sys_t structure
 *
 * \return the value of the cellwise contribution to the normalization of
 *         the residual
 */
/*----------------------------------------------------------------------------*/

static double
_eb_cw_rhs_normalization(cs_param_resnorm_type_t     type,
                         const cs_cell_mesh_t       *cm,
                         const cs_cell_sys_t        *csys)
{
  double  _rhs_norm = 0;

  switch (type) {

  case CS_PARAM_RESNORM_WEIGHTED_RHS:
    for (short int i = 0; i < cm->n_ec; i++)
      _rhs_norm += cm->pvol_e[i] * csys->rhs[i]*csys->rhs[i];
    break;

  case CS_PARAM_RESNORM_FILTERED_RHS:
    for (short int i = 0; i < csys->n_dofs; i++) {
      if (csys->dof_flag[i] & CS_CDO_BC_DIRICHLET)
        continue;
      else if (csys->dof_is_forced[i])
        continue;
      else
        _rhs_norm += csys->rhs[i]*csys->rhs[i];
    }
    break;

  case CS_PARAM_RESNORM_NORM2_RHS: /* No need to apply a filter */
    for (short int i = 0; i < csys->n_dofs; i++)
      _rhs_norm += csys->rhs[i]*csys->rhs[i];
    break;

  default:
    break; /* Nothing to do */

  } /* Type of residual normalization */

  return _rhs_norm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Perform the assembly step
 *
 * \param[in]      eqc    context for this kind of discretization
 * \param[in]      cm     pointer to a cellwise view of the mesh
 * \param[in]      csys   pointer to a cellwise view of the system
 * \param[in, out] sh     pointer to a system helper structure
 * \param[in, out] asb    pointer to a cs_cdo_assembly_t structure
 */
/*----------------------------------------------------------------------------*/

inline static void
_eb_assemble(const cs_cdoeb_vecteq_t      *eqc,
             const cs_cell_mesh_t         *cm,
             const cs_cell_sys_t          *csys,
             cs_cdo_system_helper_t       *sh,
             cs_cdo_assembly_t            *asb)
{
  assert(asb != NULL && sh != NULL); /* Sanity check */

  /* Matrix assembly */

  cs_cdo_system_block_t  *b = sh->blocks[0];
  assert(b->type == CS_CDO_SYSTEM_BLOCK_DEFAULT);
  cs_cdo_system_dblock_t  *db = b->block_pointer;

  db->assembly_func(csys->mat, csys->dof_ids, db->range_set, asb, db->mav);

  /* RHS assembly */

#if CS_CDO_OMP_SYNC_SECTIONS > 0
# pragma omp critical
  {
    for (int e = 0; e < cm->n_ec; e++)
      sh->rhs[cm->e_ids[e]] += csys->rhs[e];
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
    sh->rhs[cm->e_ids[e]] += csys->rhs[e];

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
cs_cdoeb_vecteq_is_initialized(void)
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
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_init_sharing(const cs_cdo_quantities_t    *quant,
                             const cs_cdo_connect_t       *connect,
                             const cs_time_step_t         *time_step)
{
  /* Assign static const pointers */

  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  /* Structure used to build the final system by a cell-wise process */

  assert(cs_glob_n_threads > 0);  /* Sanity check */

  BFT_MALLOC(cs_cdoeb_cell_system, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(cs_cdoeb_cell_builder, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_cdoeb_cell_system[i] = NULL;
    cs_cdoeb_cell_builder[i] = NULL;
  }

  const int  n_max_dofs = connect->n_max_ebyc;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_builder_t  *cb = _ebs_create_cell_builder(connect);
    cs_cdoeb_cell_builder[t_id] = cb;
    cs_cdoeb_cell_system[t_id] = cs_cell_sys_create(n_max_dofs,
                                                    connect->n_max_fbyc,
                                                    1, NULL);
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_builder_t  *cb = _ebs_create_cell_builder(connect);
  cs_cdoeb_cell_builder[0] = cb;
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
cs_cdoeb_vecteq_get(cs_cell_sys_t       **csys,
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
cs_cdoeb_vecteq_finalize_sharing(void)
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
 * \brief  Initialize a cs_cdoeb_vecteq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdoeb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdoeb_vecteq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  assert(eqp != NULL && eqb != NULL);

  /* This is a vector-valued equation but the DoF is scalar-valued since
   * it is a circulation associated to each edge */

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOEB || eqp->dim != 3)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of equation.\n"
              " Expected: scalar-valued CDO edge-based equation.", __func__);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_edges = connect->n_edges;

  cs_cdoeb_vecteq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdoeb_vecteq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  /* Dimensions of the algebraic system */

  eqc->n_dofs = n_edges;

  eqb->msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ |
    CS_FLAG_COMP_EF | CS_FLAG_COMP_FES;

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

  bool  need_eigen =
    (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
     eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) ? true : false;

  eqc->curlcurl_hodge = NULL;
  eqc->get_curlcurl = NULL;

  if (cs_equation_param_has_curlcurl(eqp)) {

    eqb->msh_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ;

    eqc->curlcurl_hodge = cs_hodge_init_context(connect,
                                                eqp->curlcurl_property,
                                                &(eqp->curlcurl_hodgep),
                                                true,        /* tensor ? */
                                                need_eigen); /* eigen ? */

    eqc->get_curlcurl = cs_hodge_get_func(__func__, eqp->curlcurl_hodgep);

  } /* curl-curl term */

  /* Essential boundary condition enforcement. The circulation along boundary
   * edges has the same behavior as enforcing a Dirichlet BC */

  BFT_MALLOC(eqc->edge_bc_flag, n_edges, cs_flag_t);
  cs_equation_bc_set_edge_flag(connect, eqb->face_bc, eqc->edge_bc_flag);

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

  /* Pre-defined structures for the discrete Hodge operator playing the role of
     the mass matrix */

  eqc->mass_hodgep.inv_pty  = false;
  eqc->mass_hodgep.algo = CS_HODGE_ALGO_COST;
  eqc->mass_hodgep.type = CS_HODGE_TYPE_EPFD;
  eqc->mass_hodgep.coef = cs_math_1ov3;

  /* Specify the algorithm */

  if (eqp->do_lumping ||
      eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG ||
      eqb->sys_flag & CS_FLAG_SYS_REAC_DIAG)
    eqc->mass_hodgep.algo = CS_HODGE_ALGO_VORONOI;

  eqc->mass_hodge = NULL;
  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) {

    eqc->mass_hodge = cs_hodge_init_context(connect,
                                            NULL,
                                            &(eqc->mass_hodgep),
                                            false,  /* need tensor ? */
                                            false); /* need eigen ? */

    eqc->get_mass = cs_hodge_get_func(__func__, eqc->mass_hodgep);

  } /* Add a mass matrix */

  if (eqp->sles_param->resnorm_type == CS_PARAM_RESNORM_WEIGHTED_RHS)
    eqb->msh_flag |= CS_FLAG_COMP_PEC;

  /* Whatever is the algorithm to solve, the algebraic system is the same */

  /* Helper structures (range set, interface set, matrix structure and all the
     assembly process) */

  cs_cdo_system_helper_t  *sh = NULL;

  sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_DEFAULT,
                                   1,         /* n_col_blocks */
                                   &n_edges,  /* col_block_size */
                                   1);        /* n_blocks */

  /* Choose the right class of matrix to avoid copy.
   * The way to perform the assembly may change if an external librairy is used
   * for solving the linear system */

  cs_cdo_system_matrix_class_t  matclass;

  switch (eqp->sles_param->solver_class) {

  case CS_PARAM_SLES_CLASS_CS:
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
    break;

  case CS_PARAM_SLES_CLASS_HYPRE:
#if defined(HAVE_HYPRE)
    matclass = CS_CDO_SYSTEM_MATRIX_HYPRE;
#else
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
#endif
    break;

  default:
    matclass = CS_CDO_SYSTEM_MATRIX_CS;
    break;

  }

  cs_cdo_system_add_dblock(sh, 0,      /* block_id */
                           matclass,
                           cs_flag_primal_edge,
                           n_edges,
                           1,
                           false,  /* interlaced (useless) */
                           true);  /* unrolled (useless) */

  cs_cdo_system_build_block(sh, 0);

  eqb->system_helper = sh;

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdoeb_vecteq_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdoeb_vecteq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdoeb_vecteq_free_context(void   *builder)
{
  cs_cdoeb_vecteq_t  *eqc = (cs_cdoeb_vecteq_t *)builder;

  if (eqc == NULL)
    return eqc;

  BFT_FREE(eqc->edge_bc_flag);
  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->edge_values);
  if (eqc->edge_values_pre != NULL)
    BFT_FREE(eqc->edge_values_pre);

  cs_hodge_free_context(&(eqc->curlcurl_hodge));
  cs_hodge_free_context(&(eqc->mass_hodge));

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
cs_cdoeb_vecteq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context)
{
  CS_UNUSED(field_id);
  CS_UNUSED(mesh);
  CS_UNUSED(eqb);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdoeb_vecteq_t  *eqc = (cs_cdoeb_vecteq_t *)context;

  /* By default, 0 is set as initial condition for the computational domain */

  memset(eqc->edge_values, 0, quant->n_edges*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    /* Initialize values at mesh vertices */

    cs_lnum_t  *def2e_ids = (cs_lnum_t *)cs_cdo_toolbox_get_tmpbuf();
    cs_lnum_t  *def2e_idx = NULL;
    BFT_MALLOC(def2e_idx, eqp->n_ic_defs + 1, cs_lnum_t);

    cs_cdo_sync_vol_def_at_edges(eqp->n_ic_defs,
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
                                                     eqc->edge_values);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        cs_evaluate_circulation_along_edges_by_analytic(def,
                                                        t_eval,
                                                        n_e_selected,
                                                        selected_lst,
                                                        eqc->edge_values);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid way to initialize field values for eq. %s.\n",
                  __func__, eqp->name);

      } /* Switch on possible type of definition */

    } /* Loop on definitions */

  } /* Initial values to set */

  /* Set the boundary values as initial values: Compute the values of the
     circulation where it is known thanks to the BCs */

  cs_equation_compute_circulation_eb(t_eval,
                                     mesh,
                                     quant,
                                     connect,
                                     eqp,
                                     eqc->edge_values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-Eb scheme.
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdoeb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_solve_steady_state(bool                        cur2prev,
                                   const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  cs_timer_t  t0 = cs_timer_time();

  cs_cdoeb_vecteq_t  *eqc = (cs_cdoeb_vecteq_t *)context;
  cs_cdo_system_helper_t  *sh = eqb->system_helper;
  cs_field_t  *fld = cs_field_by_id(field_id); /* vector-valued cell-based */

  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_edges = quant->n_edges;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];


  /* Build an array storing the values of the prescribed circulation at
     boundary */

  BFT_MALLOC(eqb->dir_values, n_edges, cs_real_t);
  memset(eqb->dir_values, 0, n_edges*sizeof(cs_real_t));

  cs_equation_compute_circulation_eb(time_eval,
                                     mesh,
                                     quant,
                                     connect,
                                     eqp,
                                     eqb->dir_values);

  if (cs_equation_param_has_internal_enforcement(eqp))
    eqb->enforced_values =
      cs_enforcement_define_at_edges(connect,
                                     eqp->n_enforcements,
                                     eqp->enforcement_params);

  /* Initialize the local system: rhs, matrix and assembler values */

  double  rhs_norm = 0.;
  cs_real_t  *rhs = NULL;

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
    /* Set variables and structures inside the OMP section so that each thread
       has its own value */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif
    assert(eqc->curlcurl_hodge != NULL);

    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_cdoeb_cell_system[t_id];
    cs_cell_builder_t  *cb = cs_cdoeb_cell_builder[t_id];
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *curlcurl_hodge = eqc->curlcurl_hodge[t_id];

    cb->t_pty_eval = time_eval; /* Dummy parameter if really steady */
    cb->t_bc_eval = time_eval;  /* Dummy parameter if really steady */
    cb->t_st_eval = time_eval;  /* Dummy parameter if really steady */

    /* Initialization of the curl-curl property value if uniform.
     * One calls this function with the boundary tag to examine all tests */

    cs_hodge_set_property_value(0, cb->t_pty_eval,
                                CS_FLAG_BOUNDARY_CELL_BY_FACE,
                                curlcurl_hodge);

    /* ---------------------------------------------
     * Main loop on cells to build the linear system
     * --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */

      _eb_init_cell_system(cm, eqp, eqb, eqc, csys, cb);

      /* Build and add the diffusion term to the local system. A mass matrix is
         also built if needed (stored it curlcurl_hodge->matrix) */

      _eb_curlcurl(eqp, eqb, eqc, cm, curlcurl_hodge, csys, cb);

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
                                        cb->t_st_eval,
                                        NULL,  /* No input structure */
                                        cb,
                                        csys->source);

        /* Update the RHS */

        for (short int i = 0; i < csys->n_dofs; i++)
          csys->rhs[i] += csys->source[i];

      }

      /* Boundary conditions */

      _eb_enforce_values(eqp, eqb, eqc, cm, curlcurl_hodge, csys, cb);

      /* Compute a norm of the RHS for the normalization of the residual
         of the linear system to solve */

      rhs_norm += _eb_cw_rhs_normalization(eqp->sles_param->resnorm_type,
                                           cm, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOEB_VECTEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY PROCESS
       * ================ */

      _eb_assemble(eqc, cm, csys, sh, asb);

    } /* Main loop on cells */

  } /* OPENMP Block */

  /* Free temporary buffers and structures */

  cs_cdo_system_helper_finalize_assembly(sh);
  cs_equation_builder_reset(eqb);

  /* Last step in the computation of the renormalization coefficient */

  cs_cdo_solve_sync_rhs_norm(eqp->sles_param->resnorm_type,
                             quant->vol_tot,
                             eqc->n_dofs,
                             rhs,
                             &rhs_norm);

  /* End of the system building */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Solve the linear system */

  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);
  cs_matrix_t  *matrix = cs_cdo_system_get_matrix(sh, 0);
  cs_range_set_t  *range_set = cs_cdo_system_get_range_set(sh, 0);

  /* Update edge arrays */

  if (cur2prev && eqc->edge_values_pre != NULL)
    memcpy(eqc->edge_values_pre, eqc->edge_values, sizeof(cs_real_t)*n_edges);

  cs_cdo_solve_scalar_system(eqc->n_dofs,
                             eqp->sles_param,
                             matrix,
                             range_set,
                             rhs_norm,
                             true, /* rhs_redux */
                             sles,
                             eqc->edge_values,
                             rhs);

  cs_timer_t  t2 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t1, &t2);

  /* Update fields associated to cells */

  if (cur2prev)
    cs_field_current_to_previous(fld);

  /* Update the vector-valued field associated to cells */

  cs_reco_ccen_edge_dofs(connect, quant, eqc->edge_values, &(fld->val));

  cs_timer_t  t3 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t2, &t3);

  cs_sles_free(sles);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Operate a current to previous operation for the field associated to
 *         this equation and potentially for related fields/arrays.
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdoeb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_current_to_previous(const cs_equation_param_t  *eqp,
                                    cs_equation_builder_t      *eqb,
                                    void                       *context)
{
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);

  cs_cdoeb_vecteq_t  *eqc = (cs_cdoeb_vecteq_t *)context;
  cs_field_t  *fld = cs_field_by_id(eqc->var_field_id);

  /* Edge values */

  if (eqc->edge_values_pre != NULL)
    memcpy(eqc->edge_values_pre, eqc->edge_values,
           sizeof(cs_real_t)*eqc->n_dofs);

  /* Cell values */

  cs_field_current_to_previous(fld);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdoeb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdoeb_vecteq_extra_post(const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *context)
{
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);
  CS_UNUSED(context);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at mesh edges (the DoFs)
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size: n_edges)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdoeb_vecteq_get_edge_values(void      *context,
                                bool       previous)
{
  cs_cdoeb_vecteq_t  *eqc = (cs_cdoeb_vecteq_t *)context;

  if (eqc == NULL)
    return NULL;

  if (previous)
    return eqc->edge_values_pre;
  else
    return eqc->edge_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the values at mesh cells from a reconstruction of edge values.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return a pointer to an array of cs_real_t (size: n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdoeb_vecteq_get_cell_values(void      *context,
                                bool       previous)
{
  cs_cdoeb_vecteq_t  *eqc = (cs_cdoeb_vecteq_t *)context;

  if (eqc == NULL)
    return NULL;

  cs_field_t  *c_field = cs_field_by_id(eqc->var_field_id);

  if (previous)
    return c_field->val_pre;
  else
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
cs_cdoeb_vecteq_read_restart(cs_restart_t    *restart,
                             const char      *eqname,
                             void            *scheme_context)
{
  CS_UNUSED(restart);
  CS_UNUSED(eqname);
  CS_UNUSED(scheme_context);
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
cs_cdoeb_vecteq_write_restart(cs_restart_t    *restart,
                              const char      *eqname,
                              void            *scheme_context)
{
  CS_UNUSED(restart);
  CS_UNUSED(eqname);
  CS_UNUSED(scheme_context);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
