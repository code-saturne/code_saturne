/*============================================================================
 * Build an algebraic CDO vertex-based system for unsteady convection diffusion
 * reaction of vector-valued equations with source terms
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

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_local.h"
#include "cs_cdo_time.h"
#include "cs_cdovb_priv.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_evaluate.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_parall.h"
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

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Algebraic system for CDO vertex-based discretization */
typedef struct _cs_cdovb_t cs_cdovb_vecteq_t;

/*============================================================================
 * Private variables
 *============================================================================*/

/* Structure to enable a full cellwise strategy during the system building */
static cs_cell_sys_t      **_vbv_cell_system = NULL;
static cs_cell_builder_t  **_vbv_cell_builder = NULL;

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
_vbv_create_cell_builder(const cs_cdo_connect_t   *connect)
{
  const int  n_vc = connect->n_max_vbyc;
  const int  n_ec = connect->n_max_ebyc;

  cs_cell_builder_t  *cb = cs_cell_builder_create();

  BFT_MALLOC(cb->ids, n_vc, int);
  memset(cb->ids, 0, n_vc*sizeof(int));

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
  cb->aux = cs_sdm_square_create(n_ec);
  cb->loc = cs_sdm_block33_create(n_vc, n_vc);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *         Define an indirection array for the enforcement of internal DoFs
 *         only if needed
 *
 * \param[in]      t_eval          time at which one evaluates BCs
 * \param[in]      mesh            pointer to a cs_mesh_t structure
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in]      eqb             pointer to a cs_equation_builder_t structure
 * \param[in, out] vtx_bc_flag      pointer to an array of BC flag for each vtx
 * \param[in, out] p_dir_values    pointer to the Dirichlet values to set
 * \param[in, out] p_enforced_ids  pointer to the list of enforced vertices
 */
/*----------------------------------------------------------------------------*/

static void
_vbv_setup(cs_real_t                      t_eval,
           const cs_mesh_t               *mesh,
           const cs_equation_param_t     *eqp,
           const cs_equation_builder_t   *eqb,
           cs_flag_t                      vtx_bc_flag[],
           cs_real_t                     *p_dir_values[],
           cs_lnum_t                     *p_enforced_ids[])
{
  assert(vtx_bc_flag != NULL);  /* Sanity check */
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_real_t  *dir_values = NULL;

  /* Compute the values of the Dirichlet BC */
  BFT_MALLOC(dir_values, 3*quant->n_vertices, cs_real_t);

  /* Compute the values of the Dirichlet BC */
  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   _vbv_cell_builder[0], /* static variable */
                                   vtx_bc_flag,
                                   dir_values);

  *p_dir_values = dir_values;

  /* Internal enforcement of DoFs  */
  if (cs_equation_param_has_internal_enforcement(eqp))
    cs_equation_build_dof_enforcement(quant->n_vertices,
                                      connect->c2v,
                                      eqp,
                                      p_enforced_ids);
  else
    *p_enforced_ids = NULL;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell.
 *          Case of vector-valued CDO-Vb schemes
 *
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in]      cell_flag    flag related to the current cell
 * \param[in]      cm           pointer to a cellwise view of the mesh
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in]      dir_values   Dirichlet values associated to each vertex
 * \param[in]      vtx_bc_flag  Flag related to BC associated to each vertex
 * \param[in]      forced_ids   indirection in case of internal enforcement
 * \param[in]      field_tn     values of the field at the last computed time
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vbv_init_cell_system(cs_real_t                       t_eval,
                      const cs_flag_t                 cell_flag,
                      const cs_cell_mesh_t           *cm,
                      const cs_equation_param_t      *eqp,
                      const cs_equation_builder_t    *eqb,
                      const cs_real_t                 dir_values[],
                      const cs_flag_t                 vtx_bc_flag[],
                      const cs_lnum_t                 forced_ids[],
                      const cs_real_t                 field_tn[],
                      cs_cell_sys_t                  *csys,
                      cs_cell_builder_t              *cb)
{
  csys->c_id = cm->c_id;
  csys->n_dofs = 3*cm->n_vc;
  csys->cell_flag = cell_flag;

  /* Cell-wise view of the linear system to build:
     Initialize the local system */
  cs_cell_sys_reset(cm->n_fc, csys); /* Generic part */

  cs_sdm_block33_init(csys->mat, cm->n_vc, cm->n_vc);

  for (int v = 0; v < cm->n_vc; v++) {
    const cs_lnum_t  v_id = cm->v_ids[v];
    const cs_real_t  *val_p = field_tn + 3*v_id;
    for (int k = 0; k < 3; k++) {
      csys->dof_ids[3*v + k] = 3*v_id + k;
      csys->val_n[3*v + k] = val_p[k];
    }
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Set the bc (specific part) */
    cs_equation_vb_set_cell_bc(cm,
                               eqp,
                               eqb->face_bc,
                               vtx_bc_flag,
                               dir_values,
                               t_eval,
                               csys,
                               cb);

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    cs_dbg_check_hmg_dirichlet_cw(__func__, csys);
#endif
  } /* Border cell */

  /* Special case to handle if enforcement by penalization or algebraic
   * This situation may happen with a tetrahedron with one vertex or an edge
   * lying on the boundary (but no face)
   */
  if (cell_flag == CS_FLAG_BOUNDARY_CELL_BY_VERTEX) {

    assert(vtx_bc_flag != NULL);

    for (int v = 0; v < cm->n_vc; v++) {
      for (int k = 0; k < 3; k++)
        csys->dof_flag[3*v+k] = vtx_bc_flag[cm->v_ids[v]];
      if (cs_cdo_bc_is_dirichlet(vtx_bc_flag[cm->v_ids[v]])) {
        csys->has_dirichlet = true;
        const cs_real_t  *bc_val = dir_values + 3*cm->v_ids[v];
        for (int k = 0; k < 3; k++)
          csys->dir_values[3*v+k] = bc_val[k];
      }
    }

  }

  /* Internal enforcement of DoFs  */
  if (cs_equation_param_has_internal_enforcement(eqp)) {

    assert(forced_ids != NULL);
    for (int v = 0; v < cm->n_vc; v++) {

      const cs_lnum_t  id = forced_ids[cm->v_ids[v]];

      if (id < 0) { /* No enforcement for this vertex */
        for (int k = 0; k < 3; k++)
          csys->intern_forced_ids[3*v+k] = -1;
      }
      else {

        /* In case of a Dirichlet BC, this BC is applied and the enforcement
           is ignored */
        for (int k = 0; k < 3; k++) {
          int  dof_id = 3*v+k;
          if (cs_cdo_bc_is_dirichlet(csys->dof_flag[dof_id]))
            csys->intern_forced_ids[dof_id] = -1;
          else {
            csys->intern_forced_ids[dof_id] = 3*id+k;
            csys->has_internal_enforcement = true;
          }
        }

      }

    } /* Loop on cell vertices */

  } /* Internal enforcement */

  /* Set the properties for this cell if not uniform */
  cs_equation_init_properties_cw(eqp, eqb, t_eval, cell_flag, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion, advection,
 *          reaction terms. If asked, a mass matrix is also computed and
 *          stored in cb->hdg.
 *          Case of vector-valued CDO-Vb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vbv_advection_diffusion_reaction(const cs_equation_param_t     *eqp,
                                  const cs_equation_builder_t   *eqb,
                                  const cs_cdovb_vecteq_t       *eqc,
                                  const cs_cell_mesh_t          *cm,
                                  cs_cell_sys_t                 *csys,
                                  cs_cell_builder_t             *cb)
{
  if (cs_equation_param_has_diffusion(eqp)) {   /* DIFFUSION TERM
                                                 * ============== */

    /* Define the local stiffness matrix: local matrix owned by the cellwise
       builder (store in cb->loc) */
    eqc->get_stiffness_matrix(eqp->diffusion_hodge, cm, cb);

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
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after adding diffusion", csys);
#endif
  }

  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                  * =========== */

    /* Build the mass matrix and store it in cb->hdg */
    eqc->get_mass_matrix(eqc->hdg_mass, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys)) {
      cs_log_printf(CS_LOG_DEFAULT, ">> Cell mass matrix");
      cs_sdm_dump(csys->c_id, csys->dof_ids, csys->dof_ids, cb->hdg);
    }
#endif
  }

  if (cs_equation_param_has_reaction(eqp)) { /* REACTION TERM
                                              * ============= */

    if (eqb->sys_flag & CS_FLAG_SYS_REAC_DIAG) {

      /* |c|*wvc = |dual_cell(v) cap c| */
      assert(cs_eflag_test(eqb->msh_flag, CS_FLAG_COMP_PVQ));
      const double  ptyc = cb->rpty_val * cm->vol_c;

      /* Only the diagonal block and its diagonal entries are modified */
      for (int bi = 0; bi < cm->n_vc; bi++) {

        const double  vpty = cm->wvc[bi] * ptyc;

        /* Retrieve the 3x3 matrix */
        cs_sdm_t  *bii = cs_sdm_get_block(csys->mat, bi, bi);
        assert(bii->n_rows == 3 && bii->n_cols == 3);

        bii->val[0] += vpty;
        bii->val[4] += vpty;
        bii->val[8] += vpty;

      }

    }
    else {

      assert(cs_flag_test(eqb->sys_flag, CS_FLAG_SYS_MASS_MATRIX));

      /* Add the local mass matrix (stored in cb->hdg) to the local system */
      const cs_real_t  *mval = cb->hdg->val;
      for (int bi = 0; bi < cm->n_vc; bi++) {
        for (int bj = 0; bj < cm->n_vc; bj++) {

          /* Retrieve the 3x3 matrix */
          cs_sdm_t  *bij = cs_sdm_get_block(csys->mat, bi, bj);
          assert(bij->n_rows == bij->n_cols && bij->n_rows == 3);

          const cs_real_t  _val = cb->rpty_val * mval[cm->n_vc*bi+bj];
          bij->val[0] += _val;
          bij->val[4] += _val;
          bij->val[8] += _val;

        }
      }

    } /* Lumping or not lumping */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after adding reaction", csys);
#endif
  } /* Reaction term */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   First pass to apply boundary conditions enforced weakly. Update
 *          the local system before applying the time scheme.
 *          Case of vector-valued CDO-Vb schemes
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
_vbv_apply_weak_bc(cs_real_t                      time_eval,
                   const cs_equation_param_t     *eqp,
                   const cs_cdovb_vecteq_t       *eqc,
                   const cs_cell_mesh_t          *cm,
                   cs_face_mesh_t                *fm,
                   cs_cell_sys_t                 *csys,
                   cs_cell_builder_t             *cb)
{
  CS_UNUSED(time_eval);

  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* The enforcement of the Dirichlet has to be done after all
       other contributions */
    if (cs_equation_param_has_diffusion(eqp)) {

      if (csys->has_dirichlet) /* csys is updated inside (matrix and rhs) */
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
          eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

    }

    if (csys->has_sliding)
      eqc->enforce_sliding(eqp, cm, fm, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after BC weak treatment", csys);
#endif
  } /* Cell with at least one boundary face */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Second pass to apply boundary conditions. Only Dirichlet BCs which
 *          are enforced strongly. Apply also the enforcement of internal DoFs.
 *          Update the local system after applying the time scheme.
 *          Case of vector-valued CDO-Vb schemes
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] fm          pointer to a facewise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_vbv_enforce_values(const cs_equation_param_t     *eqp,
                    const cs_cdovb_vecteq_t       *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_face_mesh_t                *fm,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  if (csys->cell_flag > 0 && csys->has_dirichlet) {

    /* Boundary element (through either vertices or faces) */

    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED) {

      /* csys is updated inside (matrix and rhs) */
      eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after strong BC treatment", csys);
#endif
    }
  }

  if (cs_equation_param_has_internal_enforcement(eqp) == false)
    return;

  /* Internal enforcement of DoFs: Update csys (matrix and rhs) */
  if (csys->has_internal_enforcement) {

    cs_equation_enforced_internal_block_dofs(eqp, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement",
                       csys);
#endif
  }
}

/*! \endcond DOXYGEN_SHOULD_SKIP_THIS */

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
  if (_vbv_cell_system == NULL || _vbv_cell_builder == NULL)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief    Allocate work buffer and general structures related to CDO
 *           vector-valued vertex-based schemes
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

  BFT_MALLOC(_vbv_cell_system, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(_vbv_cell_builder, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    _vbv_cell_system[i] = NULL;
    _vbv_cell_builder[i] = NULL;
  }

  const int  n_max_dofs = 3*connect->n_max_vbyc;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_builder_t  *cb = _vbv_create_cell_builder(connect);
    _vbv_cell_builder[t_id] = cb;

    int  block_size = 3;
    _vbv_cell_system[t_id] = cs_cell_sys_create(n_max_dofs,
                                                connect->n_max_fbyc,
                                                1,
                                                &block_size);
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_builder_t  *cb = _vbv_create_cell_builder(connect);
  _vbv_cell_builder[0] = cb;

  int  block_size = 3;
  _vbv_cell_system[0] =  cs_cell_sys_create(n_max_dofs,
                                            connect->n_max_fbyc,
                                            1,
                                            &block_size);
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

  *csys = _vbv_cell_system[t_id];
  *cb = _vbv_cell_builder[t_id];
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
    cs_cell_sys_free(&(_vbv_cell_system[t_id]));
    cs_cell_builder_free(&(_vbv_cell_builder[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(_vbv_cell_system[0]));
  cs_cell_builder_free(&(_vbv_cell_builder[0]));
#endif /* openMP */

  BFT_FREE(_vbv_cell_system);
  BFT_FREE(_vbv_cell_builder);
  _vbv_cell_builder = NULL;
  _vbv_cell_system = NULL;
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
  eqb->msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PVQ | CS_FLAG_COMP_PE |
    CS_FLAG_COMP_EV;

  /* Store additional flags useful for building boundary operator.
     Only activated on boundary cells */
  eqb->bd_msh_flag = CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_FE |
    CS_FLAG_COMP_FEQ;

  /* Diffusion */
  eqc->get_stiffness_matrix = NULL;
  if (cs_equation_param_has_diffusion(eqp)) {

    if (eqp->diffusion_hodge.is_iso == false)
      bft_error(__FILE__, __LINE__, 0, " %s: Case not handle yet\n", __func__);

    switch (eqp->diffusion_hodge.algo) {

    case CS_PARAM_HODGE_ALGO_COST:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
      if (eqp->diffusion_hodge.is_iso || eqp->diffusion_hodge.is_unity)
        eqc->get_stiffness_matrix = cs_hodge_vb_cost_get_iso_stiffness;
      else
        eqc->get_stiffness_matrix = cs_hodge_vb_cost_get_aniso_stiffness;
      eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ;
      break;

    case CS_PARAM_HODGE_ALGO_BUBBLE:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
      if (eqp->diffusion_hodge.is_iso || eqp->diffusion_hodge.is_unity)
        eqc->get_stiffness_matrix = cs_hodge_vb_bubble_get_iso_stiffness;
      else
        eqc->get_stiffness_matrix = cs_hodge_vb_bubble_get_aniso_stiffness;
      eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ;
      break;

    case CS_PARAM_HODGE_ALGO_OCS2:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ | CS_FLAG_COMP_SEF;
      eqc->get_stiffness_matrix = cs_hodge_vb_ocs2_get_aniso_stiffness;
      eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ;
      break;

    case CS_PARAM_HODGE_ALGO_VORONOI:
      eqb->msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_voro_get_stiffness;
      break;

    case CS_PARAM_HODGE_ALGO_WBS:
      eqb->msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ |
        CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ;
      eqc->get_stiffness_matrix = cs_hodge_vb_wbs_get_stiffness;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of algorithm to build the diffusion term.",
                __func__);

    } /* Switch on Hodge algo. */

  } /* Has diffusion */

  /* Boundary conditions */
  BFT_MALLOC(eqc->vtx_bc_flag, n_vertices, cs_flag_t);
  cs_equation_set_vertex_bc_flag(connect, eqb->face_bc, eqc->vtx_bc_flag);

  eqc->enforce_dirichlet = NULL;
  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_block_dirichlet;
    break;
  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_block_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PEQ |
      CS_FLAG_COMP_HFQ;
    eqc->enforce_dirichlet = cs_cdo_diffusion_vvb_ocs_weak_dirichlet;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  eqc->enforce_sliding = NULL;
  if (eqb->face_bc->n_sliding_faces > 0) {
    /* There is at least one face with a sliding condition to handle */
    eqb->bd_msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_HFQ;
    eqc->enforce_sliding = cs_cdo_diffusion_vvb_ocs_sliding;
  }

  /* Advection */
  eqc->get_advection_matrix = NULL;
  eqc->add_advection_bc = NULL;

  /* Reaction */
  if (cs_equation_param_has_reaction(eqp)) {

    if (eqp->do_lumping)
      eqb->sys_flag |= CS_FLAG_SYS_REAC_DIAG;
    else {

      switch (eqp->reaction_hodge.algo) {

      case CS_PARAM_HODGE_ALGO_VORONOI:
        eqb->sys_flag |= CS_FLAG_SYS_REAC_DIAG;
        break;
      case CS_PARAM_HODGE_ALGO_WBS:
        eqb->msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ
          | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid choice of algorithm for the reaction term.",
                  __func__);
        break;
      }

    } /* Lumping or not lumping */

  } /* Reaction */

  /* Time */
  if (cs_equation_param_has_time(eqp)) {

    if (eqp->do_lumping)
      eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
    else {

      switch (eqp->time_hodge.algo) {

      case CS_PARAM_HODGE_ALGO_VORONOI:
        eqb->sys_flag |= CS_FLAG_SYS_TIME_DIAG;
        break;
      case CS_PARAM_HODGE_ALGO_WBS:
        eqb->msh_flag |= CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ
          | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ;
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid choice of algorithm for the time term.",
                  __func__);
        break;
      }

    } /* Lumping or not lumping */

  } /* Time */

  /* Source term */
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

  /* Assembly process */
  eqc->assemble = cs_equation_assemble_set(CS_SPACE_SCHEME_CDOVB,
                                           CS_CDO_CONNECT_VTX_VECT);

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
  BFT_FREE(eqc->vtx_bc_flag);

  /* Last free */
  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of vector-valued CDO-Vb schemes.
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
cs_cdovb_vecteq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_real_t  *v_vals = fld->val;

  /* By default, 0 is set as initial condition for the computational domain.

     Warning: This operation has to be done after the settings of the
     Dirichlet boundary conditions where an interface sum is performed
     for vertex-based schemes
  */

  memset(v_vals, 0, 3*quant->n_vertices*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    cs_lnum_t  *def2v_ids = (cs_lnum_t *)cs_equation_get_tmpbuf();
    cs_lnum_t  *def2v_idx = NULL;
    BFT_MALLOC(def2v_idx, eqp->n_ic_defs + 1, cs_lnum_t);

    cs_equation_sync_vol_def_at_vertices(connect,
                                         eqp->n_ic_defs,
                                         eqp->ic_defs,
                                         def2v_idx,
                                         def2v_ids);

    /* Initialize values at mesh vertices */
    cs_flag_t  dof_flag = CS_FLAG_VECTOR | cs_flag_primal_vtx;

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_defs[def_id];
      const cs_lnum_t  n_v_selected = def2v_idx[def_id+1] - def2v_idx[def_id];
      const cs_lnum_t  *selected_lst = def2v_ids + def2v_idx[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_at_vertices_by_value(def,
                                                   n_v_selected,
                                                   selected_lst,
                                                   v_vals);
        break;

      case CS_XDEF_BY_QOV:
        /* Synchronization is performed inside */
        cs_evaluate_potential_by_qov(dof_flag, def, v_vals, NULL);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        assert(eqp->dof_reduction == CS_PARAM_REDUCTION_DERHAM);
        cs_evaluate_potential_at_vertices_by_analytic(def,
                                                      t_eval,
                                                      n_v_selected,
                                                      selected_lst,
                                                      v_vals);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid way to initialize field values for eq. %s.\n",
                  __func__, eqp->name);

      } /* Switch on possible type of definition */

    } /* Loop on definitions */

  } /* Initial values to set */

  /* Set the boundary values as initial values: Compute the values of the
     Dirichlet BC */
  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   _vbv_cell_builder[0], /* static variable */
                                   eqc->vtx_bc_flag,
                                   v_vals);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a vector steady-state
 *         convection/diffusion/reaction equation with a CDO-Vb scheme.
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovb_vecteq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_vecteq_solve_steady_state(const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_VECT];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  cs_timer_t  t0 = cs_timer_time();

  cs_cdovb_vecteq_t  *eqc = (cs_cdovb_vecteq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  /* Build an array storing the Dirichlet values at vertices and another one
     to detect vertices with an enforcement */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *forced_ids = NULL;

  /* First argument is set to t_cur + dt_cur even if this is a steady
   * computation since one can call this function to compute a steady-state
   * solution at each time step of an unsteady computation. */
  _vbv_setup(time_eval, mesh, eqp, eqb, eqc->vtx_bc_flag,
             &dir_values, &forced_ids);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;
  cs_real_t  rhs_norm = 0.0;

  assert(3*n_vertices == eqc->n_dofs);
  BFT_MALLOC(rhs, eqc->n_dofs, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < eqc->n_dofs; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN)                   \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav,               \
         dir_values, forced_ids, fld, rs, _vbv_cell_system,             \
         _vbv_cell_builder, rhs_norm)                                   \
  firstprivate(time_eval)
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
    cs_cell_sys_t  *csys = _vbv_cell_system[t_id];
    cs_cell_builder_t  *cb = _vbv_cell_builder[t_id];
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, time_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE reduction(+:rhs_norm)
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, eqb),
                         connect, quant, cm);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _vbv_init_cell_system(time_eval, cell_flag, cm, eqp, eqb,
                            dir_values, eqc->vtx_bc_flag, forced_ids, fld->val,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vbv_advection_diffusion_reaction(eqp, eqb, eqc, cm, csys, cb);

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
        for (int v = 0; v < csys->n_dofs; v++)
          csys->rhs[v] += csys->source[v];

      } /* End of source term */

      /* Compute a norm of the RHS for the normalization of the residual
         of the linear system to solve */
      cs_equation_cw_vect_res_normalization(eqp->sles_param.resnorm_type,
                                            cm->vol_c, csys, cm->wvc,
                                            &rhs_norm);

      /* Apply boundary conditions (those which are weakly enforced) */
      _vbv_apply_weak_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

      /* Enforce values if needed (internal or Dirichlet) */
      _vbv_enforce_values(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_VECTEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY PROCESS
       * ================ */

      eqc->assemble(csys, rs, eqa, mav);               /* Matrix assembly */

#     pragma omp critical
      {
        for (int i = 0; i < csys->n_dofs; i++)   /* RHS Assembly */
          rhs[csys->dof_ids[i]] += csys->rhs[i];
      }

      /* **********************  END OF ASSEMBLY PROCESS  ******************* */

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(forced_ids);

  cs_matrix_assembler_values_finalize(&mav);

  /* Last step in the computation of the renormalization coefficient */
  cs_equation_sync_res_normalization(eqp->sles_param.resnorm_type,
                                     eqc->n_dofs, /* 3*n_vertices */
                                     rhs,
                                     &rhs_norm);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Solve the linear system (treated as a scalar-valued system
     with 3 times more DoFs) */
  cs_sles_t  *sles = cs_sles_find_or_add(eqp->sles_param.field_id, NULL);

  cs_equation_solve_scalar_system(eqc->n_dofs, /* 3*n_vertices */
                                  eqp,
                                  matrix,
                                  rs,
                                  rhs_norm,
                                  true, /* rhs_redux */
                                  sles,
                                  fld->val,
                                  rhs);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_sles_free(sles);
  cs_matrix_destroy(&matrix);
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
