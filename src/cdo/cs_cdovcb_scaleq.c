/*============================================================================
 * Build an algebraic CDO vertex+cell-based system for unsteady convection
 * diffusion reaction of scalar-valued equations with source terms
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
#include <assert.h>
#include <string.h>
#include <float.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_local.h"
#include "cs_defs.h"
#include "cs_equation_assemble.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_evaluate.h"
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

  /* Ids related to the variable field and to the boundary flux field */
  int          var_field_id;
  int          bflux_field_id;

  /* System size */
  cs_lnum_t    n_dofs;     /* n_vertices + n_cells */

  /* Store the values of the field at cell centers and the data needed to
     compute the cell values from the vertex values. No need to synchronize
     all these quantities since they are only cellwise quantities. */
  cs_real_t   *cell_values;
  cs_real_t   *cell_rhs;   /* right-hand side related to cell dofs */

  /* Assembly process */
  cs_equation_assembly_t   *assemble;

  /* Members related to the static condensation */
  cs_real_t   *rc_tilda;   /* Acc^-1 * RHS_cell */
  cs_real_t   *acv_tilda;  /* Acc^-1 * Acv
                              Cell-vertices lower-Left block of the full matrix
                              Access to the values thanks to the c2v
                              connectivity */

  /* Array storing the value of the contribution of all source terms */
  cs_real_t   *source_terms;

  /* Boundary conditions */
  cs_cdo_enforce_bc_t      *enforce_dirichlet;
  cs_cdo_enforce_bc_t      *enforce_robin_bc;
  cs_flag_t                *vtx_bc_flag;

  /* Pointer of function to build the diffusion term */
  cs_hodge_t               *get_stiffness_matrix;

  /* Pointer of function to build the advection term */
  cs_cdovb_advection_t     *get_advection_matrix;
  cs_cdovb_advection_bc_t  *add_advection_bc;

  /* If one needs to build a local hodge op. for time and reaction */
  cs_param_hodge_t          hdg_mass;
  cs_hodge_t               *get_mass_matrix;

};

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Algebraic system for CDO vertex-based discretization */
typedef struct _cs_cdovcb_scaleq_t cs_cdovcb_scaleq_t;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */
static cs_cell_sys_t      **_vcbs_cell_system = NULL;
static cs_cell_builder_t  **_vcbs_cell_builder = NULL;

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
  BFT_MALLOC(cb->ids, size, int);
  memset(cb->ids, 0, size*sizeof(int));

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
 * \brief  Set the boundary conditions known from the settings
 *
 * \param[in]      t_eval        time at which one evaluates BCs
 * \param[in]      mesh          pointer to a cs_mesh_t structure
 * \param[in]      eqp           pointer to a cs_equation_param_t structure
 * \param[in, out] eqb           pointer to a cs_equation_builder_t structure
 * \param[in, out] vtx_bc_flag   pointer to an array of BC flag for each vtx
 * \param[in, out] p_dir_values  pointer to the Dirichlet values to set
 */
/*----------------------------------------------------------------------------*/

static void
_setup_vcb(cs_real_t                     t_eval,
           const cs_mesh_t              *mesh,
           const cs_equation_param_t    *eqp,
           cs_equation_builder_t        *eqb,
           cs_flag_t                    *vtx_bc_flag,
           cs_real_t                    *p_dir_values[])

{
  assert(vtx_bc_flag != NULL);  /* Sanity check */
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_real_t  *dir_values = NULL;

  /* Compute the values of the Dirichlet BC */
  BFT_MALLOC(dir_values, quant->n_vertices, cs_real_t);

  /* Compute the values of the Dirichlet BC */
  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   cs_shared_connect,
                                   eqp,
                                   eqb->face_bc,
                                   _vcbs_cell_builder[0], /* static variable */
                                   vtx_bc_flag,
                                   dir_values);

  *p_dir_values = dir_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      cell_flag    flag related to the current cell
 * \param[in]      cm           pointer to a cellwise view of the mesh
 * \param[in]      eqp          pointer to a cs_equation_param_t structure
 * \param[in]      eqb          pointer to a cs_equation_builder_t structure
 * \param[in]      eqc          pointer to a cs_cdovcb_scaleq_t structure
 * \param[in]      dir_values   Dirichlet values associated to each vertex
 * \param[in]      vtx_bc_flag  Flag related to BC associated to each vertex
 * \param[in]      field_tn     values of the field at the last computed time
 * \param[in]      t_eval       time at which one performs the evaluation
 * \param[in, out] csys         pointer to a cellwise view of the system
 * \param[in, out] cb           pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_init_vcb_cell_system(const cs_flag_t                 cell_flag,
                      const cs_cell_mesh_t           *cm,
                      const cs_equation_param_t      *eqp,
                      const cs_equation_builder_t    *eqb,
                      const cs_cdovcb_scaleq_t       *eqc,
                      const cs_real_t                 dir_values[],
                      const cs_flag_t                 vtx_bc_flag[],
                      const cs_real_t                 field_tn[],
                      cs_real_t                       t_eval,
                      cs_cell_sys_t                  *csys,
                      cs_cell_builder_t              *cb)
{
  /* Cell-wise view of the linear system to build */
  const short int  n_dofs = cm->n_vc + 1;

  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;
  csys->cell_flag = cell_flag;

  /* Initialize the local system */
  cs_cell_sys_reset(cm->n_fc, csys);

  cs_sdm_square_init(n_dofs, csys->mat);

  for (short int v = 0; v < cm->n_vc; v++) {
    csys->dof_ids[v] = cm->v_ids[v];
    csys->val_n[v] = field_tn[cm->v_ids[v]];
  }
  csys->dof_ids[cm->n_vc] = cm->c_id;
  csys->val_n[cm->n_vc] = eqc->cell_values[cm->c_id];

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    cs_equation_vb_set_cell_bc(cm,
                               cs_shared_connect,
                               cs_shared_quant,
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

    for (short int v = 0; v < cm->n_vc; v++) {
      csys->dof_flag[v] = vtx_bc_flag[cm->v_ids[v]];
      if (cs_cdo_bc_is_dirichlet(csys->dof_flag[v])) {
        csys->has_dirichlet = true;
        csys->dir_values[v] = dir_values[cm->v_ids[v]];
      }
    }

  }

  /* Set the properties for this cell if not uniform */
  cs_equation_init_properties_cw(eqp, eqb, t_eval, cell_flag, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion, advection,
 *          reaction terms in CDO-VCb schemes. If asked, a mass matrix is also
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
_vcb_advection_diffusion_reaction(double                         time_eval,
                                  const cs_equation_param_t     *eqp,
                                  const cs_equation_builder_t   *eqb,
                                  const cs_cdovcb_scaleq_t      *eqc,
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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after diffusion", csys);
#endif
  }

  if (cs_equation_param_has_convection(eqp)) {  /* ADVECTION TERM
                                                 * ============== */

    /* Define the local advection matrix */
    eqc->get_advection_matrix(eqp, cm, time_eval, fm, cb);

    /* Add it to the local system */
    cs_sdm_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after advection", csys);
#endif
  }

  if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                  * =========== */

    /* Build the mass matrix adn store it in cb->hdg */
    eqc->get_mass_matrix(eqc->hdg_mass, cm, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
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
      assert(cs_flag_test(eqb->msh_flag, CS_FLAG_COMP_PVQ));
      const double  ptyc = cb->rpty_val * cm->vol_c;
      for (short int i = 0; i < cm->n_vc; i++)
        csys->mat->val[i*(cm->n_vc + 1)] += 0.75 * cm->wvc[i] * ptyc;

      /* Cell DoF */
      csys->mat->val[cm->n_vc*(cm->n_vc + 1)] += 0.25 * ptyc;

    }
    else {

      assert(cs_flag_test(eqb->sys_flag, CS_FLAG_SYS_MASS_MATRIX));

      /* Update local system matrix with the reaction term
         cb->hdg corresponds to the current mass matrix */
      cs_sdm_add_mult(csys->mat, cb->rpty_val, cb->hdg);

    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after reaction", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   First pass to apply boundary conditions enforced weakly in CDO-VCb
 *          schemes. Update the local system before applying the time scheme.
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
_vcb_apply_weak_bc(cs_real_t                      time_eval,
                   const cs_equation_param_t     *eqp,
                   const cs_cdovcb_scaleq_t      *eqc,
                   const cs_cell_mesh_t          *cm,
                   cs_face_mesh_t                *fm,
                   cs_cell_sys_t                 *csys,
                   cs_cell_builder_t             *cb)
{
  /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
   * Operations that have to be performed BEFORE the static condensation */
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    if (cs_equation_param_has_convection(eqp))
      /* Apply boundary conditions related to the advection term
         csys is updated inside (matrix and rhs) */
      eqc->add_advection_bc(cm, eqp, time_eval, fm, cb, csys);

    /* Weakly enforced Dirichlet BCs for cells attached to the boundary
       csys is updated inside (matrix and rhs) */
    if (cs_equation_param_has_diffusion(eqp)) {

      if (csys->has_robin) {
        assert(eqc->enforce_robin_bc != NULL);
        eqc->enforce_robin_bc(eqp, cm, fm, cb, csys);
      }

      if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
          eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM)
        eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

    }

    /* Neumann boundary conditions */
    if (csys->has_nhmg_neumann) {
      for (short int v  = 0; v < cm->n_vc; v++)
        csys->rhs[v] += csys->neu_values[v];
    }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Cell system matrix after weak BC treatment", csys);
#endif
  } /* Cell with at least one boundary face */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Second pass to apply boundary conditions. Only Dirichlet BCs which
 *          are enforced strongly. Apply also the enforcement of internal DoFs.
 *          Update the local system after applying the time scheme and the
 *          static condensation.
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
_vcb_enforce_values(const cs_equation_param_t     *eqp,
                    const cs_cdovcb_scaleq_t      *eqc,
                    const cs_cell_mesh_t          *cm,
                    cs_face_mesh_t                *fm,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb)
{
  /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM
   * Operations that have to be performed AFTER the static condensation */
  if (csys->cell_flag > 0 && csys->has_dirichlet) {

    if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
        eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {

      /* Strongly enforced Dirichlet BCs for cells attached to the boundary
         csys is updated inside (matrix and rhs) */
      eqc->enforce_dirichlet(eqp, cm, fm, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after strong BC treatment", csys);
#endif
    }

  } /* Boundary cell */

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
_assemble(const cs_cdovcb_scaleq_t          *eqc,
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
    for (short int v = 0; v < cm->n_vc; v++)
      rhs[cm->v_ids[v]] += csys->rhs[v];
  }

  if (eqc->source_terms != NULL) {

    /* Source term assembly */
#   pragma omp critical
    {
      for (short int v = 0; v < cm->n_vc; v++)
        eqc->source_terms[cm->v_ids[v]] += csys->source[v];
    }

  }

#else  /* Use atomic barrier */

  for (short int v = 0; v < cm->n_vc; v++)
#   pragma omp atomic
    rhs[cm->v_ids[v]] += csys->rhs[v];

  if (eqc->source_terms != NULL) {

    /* Source term assembly */
    for (int v = 0; v < cm->n_vc; v++)
#     pragma omp atomic
      eqc->source_terms[cm->v_ids[v]] += csys->source[v];

  }

#endif

  if (eqc->source_terms != NULL) {
    cs_real_t  *cell_sources = eqc->source_terms + cs_shared_quant->n_vertices;

    cell_sources[cm->c_id] = csys->source[cm->n_vc];
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from a scalar-valued CDO-Vb scheme
 *
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] x        solution of the linear system (in: initial guess)
 * \param[in, out] b        right-hand side (scatter/gather if needed)
 *
 * \return the number of iterations of the linear solver
 */
/*----------------------------------------------------------------------------*/

static int
_solve_vcb_system(cs_sles_t                    *sles,
                  const cs_matrix_t            *matrix,
                  const cs_equation_param_t    *eqp,
                  cs_real_t                    *x,
                  cs_real_t                    *b)
{
  const cs_lnum_t  n_vertices = cs_shared_quant->n_vertices;

  cs_range_set_t *rset = cs_shared_connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  int  n_iters = 0;
  double  residual = DBL_MAX;
  cs_real_t  *xsol = NULL;

  const cs_lnum_t  n_scatter_elts = n_vertices;
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);

  if (n_cols > n_scatter_elts) {
    assert(cs_glob_n_ranks > 1);
    BFT_MALLOC(xsol, n_cols, cs_real_t);
    memcpy(xsol, x, n_scatter_elts*sizeof(cs_real_t));
  }
  else
    xsol = x;

  /* Prepare solving (handle parallelism) */
  cs_gnum_t  nnz = cs_equation_prepare_system(1,          /* stride */
                                              n_scatter_elts,
                                              matrix,
                                              rset,
                                              xsol, b);

  /* Solve the linear solver */
  const double  r_norm = 1.0; /* No renormalization by default (TODO) */
  const cs_param_sles_t  sles_param = eqp->sles_param;

  cs_sles_convergence_state_t  code = cs_sles_solve(sles,
                                                    matrix,
                                                    CS_HALO_ROTATION_IGNORE,
                                                    sles_param.eps,
                                                    r_norm,
                                                    &n_iters,
                                                    &residual,
                                                    b,
                                                    xsol,
                                                    0,      /* aux. size */
                                                    NULL);  /* aux. buffers */

  /* Output information about the convergence of the resolution */
  if (sles_param.verbosity > 0)
    cs_log_printf(CS_LOG_DEFAULT, "  <%s/sles_cvg> code %-d n_iters %d"
                  " residual % -8.4e nnz %lu\n",
                  eqp->name, code, n_iters, residual, nnz);

  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         xsol, x);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         b, b);

  cs_dbg_fprintf_system(eqp->name, cs_shared_time_step->nt_cur,
                        CS_CDOVCB_SCALEQ_DBG,
                        x, b, n_vertices);
#endif

  /* Free what can be freed at this stage */
  cs_sles_free(sles);

  if (n_cols > n_scatter_elts)
    BFT_FREE(xsol);

  return n_iters;
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
  if (_vcbs_cell_system == NULL || _vcbs_cell_builder == NULL)
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
  BFT_MALLOC(_vcbs_cell_system, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(_vcbs_cell_builder, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    _vcbs_cell_system[i] = NULL;
    _vcbs_cell_builder[i] = NULL;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    _vcbs_cell_system[t_id] = cs_cell_sys_create(connect->n_max_vbyc + 1,
                                                  connect->n_max_fbyc,
                                                  1, NULL);
    _vcbs_cell_builder[t_id] = _cell_builder_create(connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  _vcbs_cell_system[0] = cs_cell_sys_create(connect->n_max_vbyc + 1,
                                             connect->n_max_fbyc,
                                             1, NULL);
  _vcbs_cell_builder[0] = _cell_builder_create(connect);
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

  *csys = _vcbs_cell_system[t_id];
  *cb = _vcbs_cell_builder[t_id];
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
    cs_cell_sys_free(&(_vcbs_cell_system[t_id]));
    cs_cell_builder_free(&(_vcbs_cell_builder[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(_vcbs_cell_system[0]));
  cs_cell_builder_free(&(_vcbs_cell_builder[0]));
#endif /* openMP */

  BFT_FREE(_vcbs_cell_system);
  BFT_FREE(_vcbs_cell_builder);
  _vcbs_cell_builder = NULL;
  _vcbs_cell_system = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovcb_scaleq_t structure storing data useful
 *         for building and  managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovcb_scaleq_init_context(const cs_equation_param_t   *eqp,
                              int                          var_id,
                              int                          bflux_id,
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

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  /* System dimension */
  eqc->n_dofs = n_vertices + n_cells;

  /* Flag to indicate what to build in a cell mesh */
  eqb->msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PVQ | CS_FLAG_COMP_DEQ |
    CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_EV |
    CS_FLAG_COMP_FE  | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ;
  eqb->bd_msh_flag = 0;

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
  eqc->get_stiffness_matrix = NULL;
  eqc->enforce_robin_bc = NULL;
  if (cs_equation_param_has_diffusion(eqp)) {
    eqc->get_stiffness_matrix = cs_hodge_vcb_get_stiffness;
    eqc->enforce_robin_bc = cs_cdo_diffusion_svb_wbs_robin;
  }

  /* Dirichlet boundary condition enforcement */
  BFT_MALLOC(eqc->vtx_bc_flag, n_vertices, cs_flag_t);
  cs_equation_set_vertex_bc_flag(connect, eqb->face_bc, eqc->vtx_bc_flag);

  eqc->enforce_dirichlet = NULL;
  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_dirichlet = cs_cdo_diffusion_alge_dirichlet;
    break;
  case CS_PARAM_BC_ENFORCE_PENALIZED:
    eqc->enforce_dirichlet = cs_cdo_diffusion_pena_dirichlet;
    break;
  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    eqc->enforce_dirichlet = cs_cdo_diffusion_vcb_weak_dirichlet;
    break;
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    eqc->enforce_dirichlet = cs_cdo_diffusion_vcb_wsym_dirichlet;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);
    break;
  }

  /* Non-homogeneous Neumann BCs */
  if (eqb->face_bc->n_nhmg_neu_faces > 0) {
    eqb->bd_msh_flag = CS_FLAG_COMP_FV;
  }

  /* Advection */
  eqc->get_advection_matrix = NULL;
  eqc->add_advection_bc = NULL;

  if (cs_equation_param_has_convection(eqp)) {

    switch (eqp->adv_scheme) {
    case CS_PARAM_ADVECTION_SCHEME_CIP:

      eqb->msh_flag |= CS_FLAG_COMP_EF;
      _set_cip_coef(eqp);

      eqc->add_advection_bc = cs_cdo_advection_vcb_bc;

      if (cs_advection_field_is_cellwise(eqp->adv_field))
        eqc->get_advection_matrix = cs_cdo_advection_vcb_cw_cst;
      else
        eqc->get_advection_matrix = cs_cdo_advection_vcb;

      break;

    case CS_PARAM_ADVECTION_SCHEME_UPWIND:
    case CS_PARAM_ADVECTION_SCHEME_SAMARSKII:
    case CS_PARAM_ADVECTION_SCHEME_SG:
    default:
      bft_error(__FILE__, __LINE__, 0,
                " Invalid advection scheme for vertex-based discretization");
    } /* Scheme */

  }
  else {

    if (eqp->default_enforcement != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
      eqb->sys_flag |= CS_FLAG_SYS_SYM; /* Algebraic system is symmetric */

  }

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
        eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid choice of algorithm for the time term.",
                  __func__);
        break;
      }

    } /* Lumping or not lumping */

  }

  /* Source term */
  bool needs_source_term_array = false;
  if (cs_equation_param_has_time(eqp) &&
      (eqp->time_scheme == CS_TIME_SCHEME_THETA ||
       eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO))
    needs_source_term_array = true;

  eqc->source_terms = NULL;
  if (cs_equation_param_has_sourceterm(eqp) && needs_source_term_array) {

    BFT_MALLOC(eqc->source_terms, eqc->n_dofs, cs_real_t);
    memset(eqc->source_terms, 0, sizeof(cs_real_t)*eqc->n_dofs);

  } /* There is at least one source term */

  /* Pre-defined a cs_hodge_builder_t struct. */
  eqc->hdg_mass.is_unity = true;
  eqc->hdg_mass.is_iso   = true;
  eqc->hdg_mass.inv_pty  = false;
  eqc->hdg_mass.type = CS_PARAM_HODGE_TYPE_VC;
  eqc->hdg_mass.algo = CS_PARAM_HODGE_ALGO_WBS;
  eqc->hdg_mass.coef = 1.0; // not useful in this case

  eqc->get_mass_matrix = cs_hodge_vcb_wbs_get;

  /* Assembly process */
  eqc->assemble = cs_equation_assemble_set(CS_SPACE_SCHEME_CDOVCB,
                                           CS_CDO_CONNECT_VTX_SCAL);

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

  BFT_FREE(eqc->vtx_bc_flag);
  BFT_FREE(eqc->source_terms);

  /* Last free */
  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-VCb schemes.
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
cs_cdovcb_scaleq_init_values(cs_real_t                     t_eval,
                             const int                     field_id,
                             const cs_mesh_t              *mesh,
                             const cs_equation_param_t    *eqp,
                             cs_equation_builder_t        *eqb,
                             void                         *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);
  cs_real_t  *v_vals = fld->val;
  cs_real_t  *c_vals = eqc->cell_values;

  /* By default, 0 is set as initial condition for the computational domain.

     Warning: This operation has to be done after the settings of the
     Dirichlet boundary conditions where an interface sum is performed
     for vertex-based schemes
  */

  memset(v_vals, 0, quant->n_vertices*sizeof(cs_real_t));
  memset(c_vals, 0, quant->n_cells*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    /* Initialize values at mesh vertices */
    cs_flag_t  v_dof_flag = CS_FLAG_SCALAR | cs_flag_primal_vtx;
    cs_flag_t  c_dof_flag = CS_FLAG_SCALAR | cs_flag_primal_cell;

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */
      const cs_xdef_t  *def = eqp->ic_defs[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_by_value(v_dof_flag, def, v_vals);
        cs_evaluate_potential_by_value(c_dof_flag, def, c_vals);
        break;

      case CS_XDEF_BY_QOV:
        cs_evaluate_potential_by_qov(v_dof_flag | cs_flag_primal_cell, def,
                                     v_vals, c_vals);
        break;

      case CS_XDEF_BY_ANALYTIC_FUNCTION:
        if (eqp->dof_reduction != CS_PARAM_REDUCTION_DERHAM)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Incompatible reduction for equation %s.\n",
                    __func__, eqp->name);
        cs_evaluate_potential_by_analytic(v_dof_flag, def, t_eval, v_vals);
        cs_evaluate_potential_by_analytic(c_dof_flag, def, t_eval, c_vals);
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
  cs_real_t  *work_v = cs_equation_get_tmpbuf();

  cs_equation_compute_dirichlet_vb(t_eval,
                                   mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   _vcbs_cell_builder[0], /* static variable */
                                   eqc->vtx_bc_flag,
                                   work_v);


  for (cs_lnum_t v = 0; v < quant->n_vertices; v++)
    if (cs_cdo_bc_is_dirichlet(eqc->vtx_bc_flag[v]))
      v_vals[v] = work_v[v];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar steady-state
 *         convection/diffusion/reaction equation with a CDO-VCb scheme
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_solve_steady_state(const cs_mesh_t            *mesh,
                                    const int                   field_id,
                                    const cs_equation_param_t  *eqp,
                                    cs_equation_builder_t      *eqb,
                                    void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  time_eval = ts->t_cur + ts->dt[0];

  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  cs_timer_t  t0 = cs_timer_time();

  /* Build an array storing the Dirichlet values at vertices and another one
     with a tags to detect vertices related to a Neumann BC */
  cs_real_t  *dir_values = NULL;

  /* First argument is set to t_cur even if this is a steady computation since
   * one can call this function to compute a steady-state solution at each time
   * step of an unsteady computation. */
  _setup_vcb(time_eval, mesh, eqp, eqb, eqc->vtx_bc_flag, &dir_values);

  if (eqb->init_step)
    eqb->init_step = false;

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, dir_values,   \
         fld, rs, _vcbs_cell_system, _vcbs_cell_builder)
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
    cs_cell_sys_t  *csys = _vcbs_cell_system[t_id];
    cs_cell_builder_t  *cb = _vcbs_cell_builder[t_id];
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);

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
      _init_vcb_cell_system(cell_flag, cm, eqp, eqb, eqc,
                            dir_values, eqc->vtx_bc_flag, fld->val, time_eval,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vcb_advection_diffusion_reaction(time_eval,
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
        csys->rhs[cm->n_vc] += csys->source[cm->n_vc];

      } /* End of term source */

      /* Apply boundary conditions (those which are weakly enforced) */
      _vcb_apply_weak_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

      /* STATIC CONDENSATION
       * ===================
       * of the local system matrix of size n_vc + 1 into a matrix of size n_vc.
       * Store data in rc_tilda and acv_tilda to compute the values at cell
       * centers after solving the system */
      cs_static_condensation_scalar_eq(connect->c2v,
                                       eqc->rc_tilda, eqc->acv_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system matrix after condensation", csys);
#endif

      /* Enforce values if needed (internal or Dirichlet) */
      _vcb_enforce_values(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY PROCESS
       * ================ */

      _assemble(eqc, cm, csys, rs, eqa, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Solve the linear system */
  _solve_vcb_system(cs_sles_find_or_add(field_id, NULL),
                    matrix, eqp,
                    fld->val, rhs);

  /* Update field */
  t0 = cs_timer_time();

  /* Compute values at cells pc = acc^-1*(RHS - Acv*pv) */
  cs_static_condensation_recover_scalar(cs_shared_connect->c2v,
                                        eqc->rc_tilda,
                                        eqc->acv_tilda,
                                        fld->val,
                                        eqc->cell_values);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-VCb scheme
 *         Time scheme is an implicit Euler
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_solve_implicit(const cs_mesh_t            *mesh,
                                const int                   field_id,
                                const cs_equation_param_t  *eqp,
                                cs_equation_builder_t      *eqb,
                                void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  time_eval = t_cur + ts->dt[0];
  const cs_real_t  inv_dtcur = 1./ts->dt[0];

  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_EULER_IMPLICIT);

  cs_timer_t  t0 = cs_timer_time();

  /* Build an array storing the Dirichlet values at vertices and another one
     with a tags to detect vertices related to a Neumann BC */
  cs_real_t  *dir_values = NULL;

  _setup_vcb(time_eval, mesh, eqp, eqb, eqc->vtx_bc_flag, &dir_values);

  if (eqb->init_step)
    eqb->init_step = false;

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)     \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, dir_values,   \
         fld, rs, _vcbs_cell_system, _vcbs_cell_builder)
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
    cs_cell_sys_t  *csys = _vcbs_cell_system[t_id];
    cs_cell_builder_t  *cb = _vcbs_cell_builder[t_id];
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);

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
      _init_vcb_cell_system(cell_flag, cm, eqp, eqb, eqc,
                            dir_values, eqc->vtx_bc_flag, fld->val, time_eval,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vcb_advection_diffusion_reaction(time_eval,
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
        csys->rhs[cm->n_vc] += csys->source[cm->n_vc];

      } /* End of term source */

      /* Apply boundary conditions (those which are weakly enforced) */
      _vcb_apply_weak_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

      /* UNSTEADY TERM + TIME SCHEME
       * =========================== */

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        /* |c|*wvc = |dual_cell(v) cap c| */
        assert(cs_flag_test(eqb->msh_flag, CS_FLAG_COMP_PVQ));
        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
           >> Update the cellwise system with the time matrix */
        for (short int i = 0; i < cm->n_vc; i++) {

          const double  dval = 0.75 * ptyc * cm->wvc[i];
          /* Update the RHS with values at time t_n */
          csys->rhs[i] += dval * csys->val_n[i];
          /* Add the diagonal contribution from time matrix */
          csys->mat->val[i*(csys->n_dofs + 1)] += dval;

        }

        /* Cell row */
        const double  dvalc = 0.25 * ptyc;
        /* Update the RHS with values at time t_n */
        csys->rhs[cm->n_vc] += dvalc * csys->val_n[cm->n_vc];
        /* Add the diagonal contribution from time matrix */
        csys->mat->val[cm->n_vc*(csys->n_dofs + 1)] += dvalc;

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val * inv_dtcur;
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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after adding time", csys);
#endif

      /* STATIC CONDENSATION
       * ===================
       * of the local system matrix of size n_vc + 1 into a matrix of size n_vc.
       * Store data in rc_tilda and acv_tilda to compute the values at cell
       * centers after solving the system */
      cs_static_condensation_scalar_eq(connect->c2v,
                                       eqc->rc_tilda, eqc->acv_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system matrix after condensation", csys);
#endif

      /* Enforce values if needed (internal or Dirichlet) */
      _vcb_enforce_values(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY PROCESS
       * ================ */

      _assemble(eqc, cm, csys, rs, eqa, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Solve the linear system */
  _solve_vcb_system(cs_sles_find_or_add(field_id, NULL),
                    matrix, eqp,
                    fld->val, rhs);

  /* Update field */
  t0 = cs_timer_time();

  /* Compute values at cells pc = acc^-1*(RHS - Acv*pv) */
  cs_static_condensation_recover_scalar(cs_shared_connect->c2v,
                                        eqc->rc_tilda,
                                        eqc->acv_tilda,
                                        fld->val,
                                        eqc->cell_values);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve the linear system arising from a scalar unsteady
 *         convection/diffusion/reaction equation with a CDO-VCb scheme
 *         Time scheme is a theta scheme.
 *         One works cellwise and then process to the assembly.
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_id   id of the variable field related to this equation
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to cs_cdovcb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_solve_theta(const cs_mesh_t            *mesh,
                             const int                   field_id,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *context)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = connect->range_sets[CS_CDO_CONNECT_VTX_SCAL];
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_vertices = quant->n_vertices;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  inv_dtcur = 1/dt_cur;
  const cs_real_t  tcoef = 1 - eqp->theta;

  /* Time_eval = (1-theta).t^n + theta.t^(n+1) = t^n + theta.dt
   * since t^(n+1) = t^n + dt
   */
  const cs_real_t  time_eval = t_cur + eqp->theta*dt_cur;

  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t *)context;
  cs_field_t  *fld = cs_field_by_id(field_id);

  assert(cs_equation_param_has_time(eqp) == true);
  assert(eqp->time_scheme == CS_TIME_SCHEME_THETA ||
         eqp->time_scheme == CS_TIME_SCHEME_CRANKNICO);

  cs_timer_t  t0 = cs_timer_time();

  /* Build an array storing the Dirichlet values at vertices and another one
     with a tags to detect vertices related to a Neumann BC */
  cs_real_t  *dir_values = NULL;

  _setup_vcb(t_cur + dt_cur, mesh, eqp, eqb, eqc->vtx_bc_flag, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_ms);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, n_vertices, cs_real_t);
# pragma omp parallel for if  (n_vertices > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_vertices; i++) rhs[i] = 0.0;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_matrix_assembler_values_init(matrix, NULL, NULL);

  /* Detect the first call (in this case, we compute the initial source term)*/
  bool  compute_initial_source = false;
  if (eqb->init_step) {

    eqb->init_step = false;
    if (cs_equation_param_has_sourceterm(eqp))
      compute_initial_source = true; /* First iteration */

  }
  else {

    if (cs_equation_param_has_sourceterm(eqp)) {

      /* Add contribution of the previous computed source term
       * Only vertices (and not cells) since there is an assembly process
       * on vertices (and not on cells) */
      for (cs_lnum_t v = 0; v < n_vertices; v++)
        rhs[v] += tcoef * eqc->source_terms[v];
      memset(eqc->source_terms, 0, n_vertices * sizeof(cs_real_t));

      if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC ||
          eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED) {

        assert(eqc->vtx_bc_flag != NULL);
        for (cs_lnum_t v = 0; v < n_vertices; v++) {
          if (cs_cdo_bc_is_dirichlet(eqc->vtx_bc_flag[v]))
            rhs[v] = 0.;
        }

      } /* Algebraic or penalized enforcement is set */
    } /* At least one source term is defined */

  }   /* Not the first time step */

  /* ------------------------- */
  /* Main OpenMP block on cell */
  /* ------------------------- */

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)   \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav, dir_values, \
         fld, rs, _vcbs_cell_system,                                  \
         _vcbs_cell_builder, compute_initial_source)
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
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);
    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = _vcbs_cell_system[t_id];
    cs_cell_builder_t  *cb = _vcbs_cell_builder[t_id];
    cs_real_t  *cell_sources = eqc->source_terms + n_vertices;

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
      _init_vcb_cell_system(cell_flag, cm, eqp, eqb, eqc,
                            dir_values, eqc->vtx_bc_flag, fld->val, time_eval,
                            csys, cb);

      /* Build and add the diffusion/advection/reaction term to the local
         system. A mass matrix is also built if needed (stored it cb->hdg) */
      _vcb_advection_diffusion_reaction(time_eval,
                                        eqp, eqb, eqc, cm, fm, csys, cb);

      if (cs_equation_param_has_sourceterm(eqp)) { /* SOURCE TERM
                                                    * =========== */
        if (compute_initial_source) { /* First time step */

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
          csys->rhs[cm->n_vc] += tcoef * csys->source[cm->n_vc];

        }
        else { /* Add the contribution of the previous time step */

          /* Contribution at vertices has been done before the loop on cells */
          csys->rhs[cm->n_vc] += tcoef*cell_sources[cm->c_id];

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
        csys->rhs[cm->n_vc] += eqp->theta * csys->source[cm->n_vc];

      } /* End of term source */

      /* Apply boundary conditions (those which are weakly enforced) */
      _vcb_apply_weak_bc(time_eval, eqp, eqc, cm, fm, csys, cb);

      /* UNSTEADY TERM + TIME SCHEME
       * ===========================
       *
       * STEP.1 >> Compute the contribution of the "adr" to the RHS:
       *           tcoef*adr_pn where adr_pn = csys->mat * p_n */
      double  *adr_pn = cb->values;
      cs_sdm_square_matvec(csys->mat, csys->val_n, adr_pn);
      for (short int i = 0; i < csys->n_dofs; i++)
        csys->rhs[i] -= tcoef*adr_pn[i];

      /* STEP.2 >> Multiply csys->mat by theta */
      for (int i = 0; i < csys->n_dofs*csys->n_dofs; i++)
        csys->mat->val[i] *= eqp->theta;

      /* STEP.3 >> Handle the mass matrix
       * Two contributions for the mass matrix
       *  a) add to csys->mat
       *  b) add to rhs: mass_mat * p_n */

      if (eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        /* |c|*wvc = |dual_cell(v) cap c| */
        assert(cs_flag_test(eqb->msh_flag, CS_FLAG_COMP_PVQ));
        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* STEPS >> Compute the time contribution to the RHS: Mtime*pn
         *       >> Update the cellwise system with the time matrix */
        for (short int i = 0; i < cm->n_vc; i++) {

          const double  dval = 0.75 * ptyc * cm->wvc[i];
          /* Update the RHS with values at time t_n */
          csys->rhs[i] += dval * csys->val_n[i];
          /* Add the diagonal contribution from time matrix */
          csys->mat->val[i*(csys->n_dofs + 1)] += dval;

        }

        /* Cell row */
        const double  dvalc = 0.25 * ptyc;
        /* Update the RHS with values at time t_n */
        csys->rhs[cm->n_vc] += dvalc * csys->val_n[cm->n_vc];
        /* Add the diagonal contribution from time matrix */
        csys->mat->val[cm->n_vc*(csys->n_dofs + 1)] += dvalc;

      }
      else { /* Use the mass matrix */

        const double  tpty_coef = cb->tpty_val * inv_dtcur;
        const cs_sdm_t  *mass_mat = cb->hdg;

         /* Update rhs with mass_mat*p^n */
        double  *time_pn = cb->values;
        cs_sdm_square_matvec(mass_mat, csys->val_n, time_pn);
        for (short int i = 0; i < csys->n_dofs; i++)
          csys->rhs[i] += tpty_coef*time_pn[i];

        /* Update the cellwise system with the time matrix */
        cs_sdm_add_mult(csys->mat, tpty_coef, mass_mat);

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after adding time", csys);
#endif

      /* STATIC CONDENSATION
       * ===================
       * of the local system matrix of size n_vc + 1 into a matrix of size n_vc.
       * Store data in rc_tilda and acv_tilda to compute the values at cell
       * centers after solving the system */
      cs_static_condensation_scalar_eq(connect->c2v,
                                       eqc->rc_tilda, eqc->acv_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Cell system matrix after condensation", csys);
#endif

      /* Enforce values if needed (internal or Dirichlet) */
      _vcb_enforce_values(eqp, eqc, cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVCB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY PROCESS
       * ================ */

      _assemble(eqc, cm, csys, rs, eqa, mav, rhs);

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Copy current field values to previous values */
  cs_field_current_to_previous(fld);

  /* Solve the linear system */
  _solve_vcb_system(cs_sles_find_or_add(field_id, NULL),
                    matrix, eqp,
                    fld->val, rhs);

  /* Update field */
  t0 = cs_timer_time();

  /* Compute values at cells pc = acc^-1*(RHS - Acv*pv) */
  cs_static_condensation_recover_scalar(cs_shared_connect->c2v,
                                        eqc->rc_tilda,
                                        eqc->acv_tilda,
                                        fld->val,
                                        eqc->cell_values);

  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);

  /* Free remaining buffers */
  BFT_FREE(rhs);
  cs_matrix_destroy(&matrix);
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
cs_cdovcb_scaleq_get_vertex_values(void      *context)
{
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t *)context;
  cs_field_t  *pot = cs_field_by_id(eqc->var_field_id);

  return pot->val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at mesh cells from the inverse operation
 *         w.r.t. the static condensation (DoF used in the linear system are
 *         located at primal vertices and field related to the structure
 *         equation is also attached to primal vertices)
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdovcb_scaleq_get_cell_values(void     *context)
{
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)context;

  if (eqc == NULL)
    return NULL;
  else
    return eqc->cell_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute for each vertex of a boundary face, the portion of diffusive
 *         flux across the boundary face. The surface attached to each vertex
 *         corresponds to the intersection of its dual cell (associated to
 *         a vertex of the face) with the face.
 *         Case of scalar-valued CDO-VCb schemes
 *
 * \param[in]       t_eval     time at which one performs the evaluation
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in]       pot_v      pointer to an array of field values at vertices
 * \param[in]       pot_c      pointer to an array of field values at cells
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  vf_flux    pointer to the values of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_boundary_diff_flux(const cs_real_t              t_eval,
                                    const cs_equation_param_t   *eqp,
                                    const cs_real_t             *pot_v,
                                    const cs_real_t             *pot_c,
                                    cs_equation_builder_t       *eqb,
                                    cs_real_t                   *vf_flux)
{
  if (vf_flux == NULL)
    return;

  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(vf_flux, 0, connect->bf2v->idx[quant->n_b_faces]*sizeof(cs_real_t));

    cs_timer_t  t1 = cs_timer_time();
    cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
    return;
  }

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)        \
  shared(quant, connect, eqp, eqb, vf_flux, pot_v, pot_c, _vcbs_cell_builder)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_cdo_bc_face_t  *face_bc = eqb->face_bc;
    const cs_adjacency_t  *bf2v = connect->bf2v;
    const cs_adjacency_t  *f2c = connect->f2c;
    const cs_lnum_t  fidx_shift = f2c->idx[quant->n_i_faces];

    cs_real_t  *pot = NULL, *flux = NULL;
    BFT_MALLOC(pot, connect->n_max_vbyc + 1, cs_real_t);
    BFT_MALLOC(flux, connect->n_max_vbyc , cs_real_t);

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_builder_t  *cb = _vcbs_cell_builder[t_id];
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);

    /* msh_flag for Neumann and Robin BCs. Add add_flag for the other cases
       when one has to reconstruct a flux */
    cs_flag_t  msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_FV;
    cs_flag_t  add_flag = CS_FLAG_COMP_EV | CS_FLAG_COMP_FE | CS_FLAG_COMP_PEQ |
      CS_FLAG_COMP_PFQ | CS_FLAG_COMP_PVQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_FEQ;

    if (eqb->diff_pty_uniform)  /* c_id = 0, cell_flag = 0 */
      cs_equation_set_diffusion_property(eqp, 0, t_eval, 0, cb);

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t bf_id = 0; bf_id < quant->n_b_faces; bf_id++) {

      const cs_lnum_t  f_id = bf_id + quant->n_i_faces;
      const cs_lnum_t  c_id = f2c->ids[bf_id + fidx_shift];
      const cs_lnum_t  *idx  = bf2v->idx + bf_id;

      cs_real_t  *_flx = vf_flux + idx[0];

      switch (face_bc->flag[bf_id]) {

      case CS_CDO_BC_HMG_NEUMANN:
        memset(_flx, 0, (idx[1]-idx[0])*sizeof(cs_real_t));
        break;

      case CS_CDO_BC_NEUMANN:
        {
          cs_real_t  *neu_values = cb->values;

          /* Set the local mesh structure for the current cell */
          cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

          const short int  f = cs_cell_mesh_get_f(f_id, cm);

          cs_equation_compute_neumann_sv(t_eval,
                                         face_bc->def_ids[bf_id],
                                         f,
                                         quant,
                                         eqp,
                                         cm,
                                         neu_values);

          short int n_vf = 0;
          for (int i = cm->f2v_idx[f]; i < cm->f2v_idx[f+1]; i++)
            _flx[n_vf++] = neu_values[cm->f2v_ids[i]];

        }
        break;

      case CS_CDO_BC_ROBIN:
        {
          cs_real_t  *robin_values = cb->values;
          cs_real_t  *wvf = cb->values + 3;

          /* Set the local mesh structure for the current cell */
          cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

          const short int  f = cs_cell_mesh_get_f(f_id, cm);
          const cs_real_t  f_area = quant->b_face_surf[bf_id];

          /* Robin BC expression: K du/dn + alpha*(p - p0) = g */
          cs_equation_compute_robin(t_eval,
                                    face_bc->def_ids[bf_id],
                                    f,
                                    quant,
                                    eqp,
                                    cm,
                                    robin_values);

          const cs_real_t  alpha = robin_values[0];
          const cs_real_t  p0 = robin_values[1];
          const cs_real_t  g = robin_values[2];

          cs_cdo_quantities_compute_b_wvf(connect, quant, bf_id, wvf);

          short int n_vf = 0;
          for (int i = cm->f2v_idx[f]; i < cm->f2v_idx[f+1]; i++) {
            const cs_real_t  pv = pot_v[cm->v_ids[cm->f2v_ids[i]]];
            _flx[n_vf] = f_area*wvf[n_vf]*(alpha*(p0 - pv) + g);
            n_vf++;
          }
        }
        break;

      default:
        { /* Reconstruct a normal flux at the boundary face */

          /* Set the local mesh structure for the current cell */
          cs_cell_mesh_build(c_id, msh_flag | add_flag, connect, quant, cm);

          const short int  f = cs_cell_mesh_get_f(f_id, cm);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOVB_SCALEQ_DBG > 2
          if (cs_dbg_cw_test(eqp, cm, NULL)) cs_cell_mesh_dump(cm);
#endif
          if (!eqb->diff_pty_uniform)
            cs_property_tensor_in_cell(cm,
                                       eqp->diffusion_property,
                                       t_eval,
                                       eqp->diffusion_hodge.inv_pty,
                                       cb->dpty_mat);

          /* Define a local buffer keeping the value of the discrete potential
             for the current cell */
          for (short int v = 0; v < cm->n_vc; v++)
            pot[v] = pot_v[cm->v_ids[v]];
          pot[cm->n_vc] = pot_c[c_id];

          cs_cdo_diffusion_wbs_vbyf_flux(f, eqp, cm, pot, cb, flux);

          /* Fill the global flux array */
          short int n_vf = 0;
          for (int i = cm->f2v_idx[f]; i < cm->f2v_idx[f+1]; i++)
            _flx[n_vf++] = flux[cm->f2v_ids[i]];

        }
        break;

      } /* End of switch */

    } /* End of loop on boundary faces */

    BFT_FREE(pot);
    BFT_FREE(flux);

  } /* End of Open block */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *         Case of scalar-valued CDO-VCb schemes
 *
 * \param[in]       normal     indicate in which direction flux is > 0
 * \param[in]       pdi        pointer to an array of field values
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to data specific for this scheme
 * \param[in, out]  d_flux     pointer to the value of the diffusive flux
 * \param[in, out]  c_flux     pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_flux_across_plane(const cs_real_t             normal[],
                                   const cs_real_t            *pdi,
                                   const cs_equation_param_t  *eqp,
                                   int                         ml_id,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context,
                                   double                     *d_flux,
                                   double                     *c_flux)
{
  *d_flux = 0.;
  *c_flux = 0.;

  if (pdi == NULL)
    return;

  cs_mesh_location_type_t  ml_t = cs_mesh_location_get_type(ml_id);

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

  /* To be modified for an fully integration of openMP */
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)context;
  cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(0);
  cs_cell_builder_t  *cb = _vcbs_cell_builder[0];

  double  *p_v = NULL;
  BFT_MALLOC(p_v, connect->n_max_vbyf, double);

  if (ml_t == CS_MESH_LOCATION_BOUNDARY_FACES) {

    const cs_lnum_t  n_i_faces = connect->n_faces[2];
    const cs_lnum_t  *cell_ids = f2c->ids + f2c->idx[n_i_faces];

    for (cs_lnum_t id = 0; id < n_elts[0]; id++) { /* Loop on selected faces */

      const cs_lnum_t  bf_id = elt_ids[id];
      const cs_lnum_t  f_id = n_i_faces + bf_id;
      const cs_lnum_t  c_id = cell_ids[bf_id];

      /* Build a face-wise view of the mesh */
      cs_face_mesh_build(c_id, f_id, connect, quant, fm);

      const short int  sgn = (_dp3(fm->face.unitv, normal) < 0) ? -1 : 1;

      /* Store values of the potential related to this face */
      for (short int v = 0; v < fm->n_vf; v++)
        p_v[v] = pdi[fm->v_ids[v]];

      /* Interpolate a value at the face center */
      p_f = cs_reco_fw_scalar_pv_at_face_center(fm, p_v);

      /* Compute the local diffusive flux */
      if (cs_equation_param_has_diffusion(eqp)) {

        cs_property_get_cell_tensor(c_id, t_cur,
                                    eqp->diffusion_property,
                                    eqp->diffusion_hodge.inv_pty,
                                    pty_tens);

        flx = cs_cdo_diffusion_wbs_face_flux(fm,
                                             (const cs_real_3_t (*))pty_tens,
                                             p_v, p_f, eqc->cell_values[c_id],
                                             cb);
        *d_flux += sgn * flx;

      } /* Diffusive flux */

      /* Compute the local advective flux */
      if (cs_equation_param_has_convection(eqp)) {

        const double  coef = sgn * fm->face.meas * p_f;

        cs_advection_field_get_cell_vector(c_id, eqp->adv_field, &adv_c);
        *c_flux += coef * adv_c.meas * _dp3(adv_c.unitv, fm->face.unitv);

      }

    } /* Loop on selected border faces */

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

        /* Compute the diffusive flux seen from cell c1 */
        if (cs_equation_param_has_diffusion(eqp)) {

          cs_property_get_cell_tensor(c_id, t_cur,
                                      eqp->diffusion_property,
                                      eqp->diffusion_hodge.inv_pty,
                                      pty_tens);

          flx = cs_cdo_diffusion_wbs_face_flux(fm,
                                               (const cs_real_3_t (*))pty_tens,
                                               p_v, p_f, eqc->cell_values[c_id],
                                               cb);

          *d_flux += sgn * 0.5 * flx; /* Centered approximation */

        } /* Diffusive flux */

        if (cs_equation_param_has_convection(eqp)) {

          /* Compute the local advective flux seen from cell */
          cs_advection_field_get_cell_vector(c_id,
                                             eqp->adv_field,
                                             &adv_c);
          flx = adv_c.meas * _dp3(adv_c.unitv, fm->face.unitv);

          /* Centered approximation of the advective flux */
          *c_flux += sgn * 0.5 * flx  * p_f * fm->face.meas;

        } /* Advective flux */

      } /* Loop on cells attached to this interior face */

    } /* Loop on selected interior faces */

  } /* Set of interior or border faces */

  BFT_FREE(p_v);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux in each cells.
 *         Case of scalar-valued CDO-VCb schemes
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to data structure cast on-the-fly
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_diff_flux_in_cells(const cs_real_t             *values,
                                    const cs_equation_param_t   *eqp,
                                    cs_real_t                    t_eval,
                                    cs_equation_builder_t       *eqb,
                                    void                        *context,
                                    cs_real_t                   *diff_flux)
{
  if (diff_flux == NULL)
    return ;

  /* Sanity checks */
  assert(eqp->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_WBS);

  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)context;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(diff_flux, 0, 3*quant->n_cells);
    return;
  }

  cs_timer_t  t0 = cs_timer_time();

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)   \
  shared(quant, connect, eqp, eqb, eqc, diff_flux, values,            \
         t_eval, _vcbs_cell_builder)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    cs_flag_t  msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PFQ |
      CS_FLAG_COMP_HFQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV;

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = _vcbs_cell_builder[t_id];

    double  *pot = NULL;
    BFT_MALLOC(pot, connect->n_max_vbyc + 1, double);

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

      cs_cdo_diffusion_wbs_get_cell_flux(cm, pot, cb, diff_flux + 3*c_id);

    } /* Loop on cells */

    BFT_FREE(pot);

  } /* OMP Section */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across dual faces
 *         Case of scalar-valued CDO-VCb schemes
 *
 * \param[in]       values      discrete values for the potential
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to data structure cast on-the-fly
 * \param[in, out]  diff_flux   value of the diffusive flux
  */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_diff_flux_dfaces(const cs_real_t             *values,
                                  const cs_equation_param_t   *eqp,
                                  cs_real_t                    t_eval,
                                  cs_equation_builder_t       *eqb,
                                  void                        *context,
                                  cs_real_t                   *diff_flux)
{
  if (diff_flux == NULL)
    return ;

  /* Sanity checks */
  assert(eqp->diffusion_hodge.algo == CS_PARAM_HODGE_ALGO_WBS);

  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)context;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  if (cs_equation_param_has_diffusion(eqp) == false) {
    memset(diff_flux, 0, connect->c2e->idx[quant->n_cells]*sizeof(cs_real_t));
    return;
  }

  cs_timer_t  t0 = cs_timer_time();

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)   \
  shared(quant, connect, eqp, eqb, eqc, diff_flux, values,            \
         t_eval, _vcbs_cell_builder)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    cs_flag_t  msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PFQ |
      CS_FLAG_COMP_EFQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV;

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = _vcbs_cell_builder[t_id];
    double  *pot = NULL;
    BFT_MALLOC(pot, connect->n_max_vbyc + 1, double);

    if (eqb->diff_pty_uniform)
      cs_equation_set_diffusion_property(eqp, 0, t_eval, 0, cb);

    /* Define the flux by cellwise contributions */
#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      if (!eqb->diff_pty_uniform) /* cell_flag is set to 0 */
        cs_equation_set_diffusion_property_cw(eqp, cm, t_eval, 0, cb);

      /* Define a local buffer keeping the value of the discrete potential
         for the current cell */
      for (short int v = 0; v < cm->n_vc; v++)
        pot[v] = values[cm->v_ids[v]];
      pot[cm->n_vc] = eqc->cell_values[c_id];

      cs_cdo_diffusion_wbs_get_dfbyc_flux(cm, pot, cb,
                                          diff_flux + connect->c2e->idx[c_id]);

    } /* Loop on cells */

    BFT_FREE(pot);

  } /* OMP Section */

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the discrete gradient at vertices
 *
 * \param[in]       v_values    discrete values for the potential at vertices
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  context     pointer to data structure cast on-the-fly
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
         _vcbs_cell_builder, cs_glob_n_ranks)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif
    double  *pot = NULL;
    cs_real_3_t  cgrd;

    BFT_MALLOC(pot, connect->n_max_vbyc + 1, double);

    cs_flag_t  msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PFQ |
      CS_FLAG_COMP_DEQ | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_HFQ;

    /* Each thread get back its related structures:
       Get the cellwise view of the mesh and the algebraic system */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = _vcbs_cell_builder[t_id];

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
 * \brief  Read additional arrays (not defined as fields) but useful for the
 *         checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in]       scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovcb_scaleq_read_restart(cs_restart_t    *restart,
                              const char      *eqname,
                              void            *scheme_context)
{
  /* Only the cell values are handled. Vertex values are stored in a cs_field_t
     structure and thus are handled automatically. */
  if (restart == NULL)
    return;
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is NULL", __func__);
  if (scheme_context == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Scheme context is NULL", __func__);

  int retcode = CS_RESTART_SUCCESS;
  cs_cdovcb_scaleq_t  *eqc = (cs_cdovcb_scaleq_t  *)scheme_context;

  const int  cell_loc_id = cs_mesh_location_get_id_by_name("cells");

  /* Define the section name */
  char sec_name[128];
  snprintf(sec_name, 127, "%s::cell_vals", eqname);

  /* Check section */
  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     cell_loc_id,
                                     1, /* scalar-valued */
                                     CS_TYPE_cs_real_t);

  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      cell_loc_id,
                                      1, /* scalar-valued */
                                      CS_TYPE_cs_real_t,
                                      eqc->cell_values);
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
cs_cdovcb_scaleq_write_restart(cs_restart_t    *restart,
                               const char      *eqname,
                               void            *scheme_context)
{
  /* Only the cell values are handled. Vertex values are stored in a cs_field_t
     structure and thus are handled automatically. */
  if (restart == NULL)
    return;
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is NULL", __func__);

  const cs_cdovcb_scaleq_t  *eqc = (const cs_cdovcb_scaleq_t  *)scheme_context;
  const int  cell_loc_id = cs_mesh_location_get_id_by_name("cells");

  /* Define the section name */
  char sec_name[128];
  snprintf(sec_name, 127, "%s::cell_vals", eqname);

  /* Write section */
  cs_restart_write_section(restart,
                           sec_name,
                           cell_loc_id,
                           1,   /* scalar-valued */
                           CS_TYPE_cs_real_t,
                           eqc->cell_values);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdovcb_scaleq_t structure
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
