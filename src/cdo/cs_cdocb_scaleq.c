/*============================================================================
 * Build an algebraic CDO cell-based system for the diffusion equation
 * and solved it as one block (monolithic approach of the flux-potential
 * coupling)
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

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_blas.h"
#include "cs_domain.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_blas.h"
#include "cs_cdocb_monolithic_sles.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation_bc.h"
#include "cs_equation_builder.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_field_operator.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_sdm.h"
#include "cs_solid_selection.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdocb_priv.h"
#include "cs_cdocb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdocb_scaleq.c
 *
 * \brief Build an algebraic CDO cell-based system for the diffusion
 *        equation and solved it with a monolithic approach
 */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOCB_SCALEQ_DBG      0

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */

static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;

static cs_cell_sys_t               **_scb_cell_system = NULL;
static cs_cell_builder_t           **_scb_cell_builder = NULL;

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

  BFT_MALLOC(cb->adv_fluxes, n_fc, double);
  memset(cb->adv_fluxes, 0, n_fc*sizeof(double));

  BFT_MALLOC(cb->ids, n_fc, int);
  memset(cb->ids, 0, n_fc*sizeof(int));

  int  size = (n_fc + 1)*(n_fc + 1);
  BFT_MALLOC(cb->values, size, double);
  memset(cb->values, 0, size*sizeof(double));

  size = 2*(n_fc + 1);
  BFT_MALLOC(cb->vectors, size, cs_real_3_t);
  memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

  /* Local square dense matrices used during the construction of operators */

  cb->loc = cs_sdm_square_create(n_fc);
  cb->aux = cs_sdm_square_create(n_fc);

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define several structures such as the cs_range_set_t,
 *        cs_interface_set_t, cs_matrix_assembler_t and cs_matrix_structure_t
 *        in case of a full assembly (i.e. the full saddle-point matrix is
 *        built).
 *        A variant, activated with add_pressure_diag, is available in order to
 *        enforce the pressure.
 *
 * \param[in, out]  block              pointer to a block structure
 * \param[in]       add_pressure_diag  true or false (pressure diagonal block)
 */
/*----------------------------------------------------------------------------*/

static void
_build_shared_full_structures(cs_cdo_system_block_t     *block,
                              bool                       add_pressure_diag)
{
  /* Compute the range set for an array of size n_faces + n_cells
   * velocity is attached to faces (one for each component) and pressure
   * to cells
   *
   * Storage for the global numbering: Vel_X | Vel_Y | Vel_Z | Pressure */

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_lnum_t  size = n_faces + quant->n_cells;

  assert(block->type == CS_CDO_SYSTEM_BLOCK_EXT);
  cs_cdo_system_xblock_t  *xb = block->block_pointer;

  /* 1. Build the interface set and the range set structures */

  cs_interface_set_t *ifs = connect->face_ifs;

  if (ifs != NULL)
    xb->interface_set = cs_interface_set_dup_blocks(ifs, n_faces, 1);

  xb->range_set = cs_range_set_create(xb->interface_set,
                                      NULL,   /* halo */
                                      size,
                                      false,  /* TODO: add balance option */
                                      1,      /* tr_ignore */
                                      0);     /* g_id_base */

  /* 2. Build the matrix assembler structure */

  const cs_adjacency_t  *f2f = connect->f2f;
  const cs_adjacency_t  *f2c = connect->f2c;

  /* The second parameter is set to "true" meaning that the diagonal is stored
   * separately --> MSR storage
   * Create the matrix assembler structure */

  xb->matrix_assembler = cs_matrix_assembler_create(xb->range_set->l_range,
                                                    true);

  /* First loop to count max size of the buffer used to fill the matrix
   * structure. +1 to take into account the diagonal term. */

  int  max_sten = 0;
  for (cs_lnum_t f = 0; f < n_faces; f++) {
    int  sten
      = (f2f->idx[f+1]-f2f->idx[f] + 1) + 2*(f2c->idx[f+1]-f2c->idx[f]);
    max_sten = CS_MAX(max_sten, sten);
  }

  cs_gnum_t  *grows = NULL, *gcols = NULL;
  BFT_MALLOC(grows, max_sten, cs_gnum_t);
  BFT_MALLOC(gcols, max_sten, cs_gnum_t);

  /*
   *
   *     |        |         |
   *     |   A    |    Bt   |  B is the divergence (Bt the gradient)
   *     |        |         |  A is csys->mat in what follows
   *     |--------|---------|  The viscous part arising from the CDO-Cb
   *     |        |         |  schemes for scalar-valued variables and
   *     |   B    |    0    |  additional terms as the linearized
   *     |        |         |  convective term
   *

   *  Block A is n_faces * n_faces
   *  Block B is n_cells * n_faces
   */

  /* Only on faces (B_x is build in the same time as Bt_x for pressure DoFs) */

  for (cs_lnum_t frow_id = 0; frow_id < n_faces; frow_id++) {

    const cs_lnum_t  start = f2f->idx[frow_id];
    const cs_lnum_t  end = f2f->idx[frow_id+1];

    /* A face-face entry + the diagonal which is not
       taken into account in the face --> face connectivity. The B and Bt
       operators have the same sparsity. An entry for the c2f
       connectivity. This is multiply by two since one considers B and Bt. */

    int  n_entries = (end-start + 1)
                   + 2*(f2c->idx[frow_id+1]-f2c->idx[frow_id]);

    const cs_gnum_t  grow_id =  xb->range_set->g_id[frow_id];

    int shift = 0;

    /* Diagonal term is excluded in this connectivity. Add it "manually" */

    grows[shift] = grow_id;
    gcols[shift] = grow_id;
    shift++;

    /* Extra diagonal couples */

    for (cs_lnum_t idx = start; idx < end; idx++) {

      const cs_lnum_t  fcol_id = f2f->ids[idx];
      const cs_gnum_t  gcol_id = xb->range_set->g_id[fcol_id];

      grows[shift] = grow_id;
      gcols[shift] = gcol_id;
      shift++;

    } /* Loop on extra-diag. entries */

    /* Loop on pressure-related  entries */

    for (cs_lnum_t idx = f2c->idx[frow_id]; idx < f2c->idx[frow_id+1]; idx++) {

      const cs_lnum_t  ccol_id = f2c->ids[idx];
      const cs_gnum_t  gcol_id = xb->range_set->g_id[n_faces + ccol_id];

      grows[shift] = grow_id;
      gcols[shift] = gcol_id;
      shift++;

      /* Its transposed */

      grows[shift] = gcol_id;
      gcols[shift] = grow_id;
      shift++;


    } /* Loop on pressure related DoFs */

    cs_matrix_assembler_add_g_ids(xb->matrix_assembler,
                                  n_entries, grows, gcols);
    assert(shift == n_entries);

  } /* Loop on face entities */

  if (add_pressure_diag) {

    const cs_gnum_t  *cell_g_ids = xb->range_set->g_id + n_faces;

    cs_matrix_assembler_add_g_ids(xb->matrix_assembler,
                                  quant->n_cells, cell_g_ids, cell_g_ids);

  }

  /* 3. Build the matrix structure */

  cs_matrix_assembler_compute(xb->matrix_assembler);

  xb->matrix_structure
    = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR,
                                                xb->matrix_assembler);

  /* Free temporary buffers */

  BFT_FREE(grows);
  BFT_FREE(gcols);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the cell system when this should
 *          be done before the static condensation.
 *          Case of CDO-Cb schemes with a monolithic flux-potential coupling
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       eqc              context structure for a scalar-valued Cb
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       fm        pointer to a cs_face_mesh_t structure
 * \param[in]       hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_scb_apply_bc(const cs_equation_param_t     *eqp,
              const cs_cdocb_scaleq_t       *eqc,
              const cs_cell_mesh_t          *cm,
              cs_face_mesh_t                *fm,
              cs_hodge_t                    *hodge,
              cs_cell_sys_t                 *csys,
              cs_cell_builder_t             *cb)
{
  assert(cs_equation_param_has_diffusion(eqp) == true);

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Neumann boundary conditions:
     * The common practice is to define Phi_neu = - lambda * grad(u) . n_fc
     * An outward flux is a positive flux whereas an inward flux is negative
     * The minus just above implies the minus just below */

    if (csys->has_nhmg_neumann)
      eqc->enforce_neumann(eqp, cm, fm, hodge, cb, csys);

    /* The enforcement of the Dirichlet has to be done after all
       other contributions */

    if (csys->has_dirichlet) /* csys is updated inside (matrix and rhs) */
       eqc->enforce_dirichlet(eqp, cm, fm, hodge, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOCB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Cell system matrix after partial BC enforcement",
                       csys);
#endif
  } /* Boundary cell */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a scalar-valued system obtained
 *         with CDO-Cb schemes
 *
 * \param[in]      csys             pointer to a cs_cell_sys_t structure
 * \param[in]      cm               pointer to a cs_cell_mesh_t structure
 * \param[in]      eqb              pointer to a cs_equation_builder_t structure
 * \param[in, out] eqc              context structure for a scalar-valued Cb
 * \param[in, out] asb              pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_full_assembly(const cs_cell_sys_t            *csys,
               const cs_cell_mesh_t           *cm,
               cs_equation_builder_t          *eqb,
               cs_cdocb_scaleq_t              *eqc,
               cs_cdo_assembly_t              *asb)
{
  CS_UNUSED(asb);

  const short int  n_f = cm->n_fc;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_sdm_t  *m = csys->mat;

  cs_cdo_system_helper_t  *sh = eqb->system_helper;
  cs_cdo_system_block_t  *b = sh->blocks[0];
  assert(b->type == CS_CDO_SYSTEM_BLOCK_EXT);
  cs_cdo_system_xblock_t  *xb = b->block_pointer;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

  cs_real_t *div_op = eqc->div_op_cw[t_id];
  cs_real_t *rhs = sh->rhs;

  const cs_range_set_t  *rset = xb->range_set;

  /* 1. Matrix assembly
   * ================== */

  cs_gnum_t  r_gids[CS_CDO_ASSEMBLE_BUF_SIZE];
  cs_gnum_t  c_gids[CS_CDO_ASSEMBLE_BUF_SIZE];
  cs_real_t  values[CS_CDO_ASSEMBLE_BUF_SIZE];

  const cs_gnum_t  p_gid = rset->g_id[n_faces + cm->c_id];

  /* 1.a Add the contribution of velocity DoFs */

  int  bufsize = 0;
  for (int i = 0; i < m->n_rows; i++) {

    /* dof_ids is an interlaced array and global numbering is not interlaced
       that's why one considers f_id */

    const cs_lnum_t  fi_id = cm->f_ids[i];
    const cs_gnum_t  i_gid = rset->g_id[fi_id];

    for (int j = 0; j < m->n_rows; j++) {

      const cs_lnum_t  fj_id = cm->f_ids[j];
      const cs_gnum_t  j_gid = rset->g_id[fj_id];

      /* Add an entry */

      r_gids[bufsize] = i_gid;
      c_gids[bufsize] = j_gid;
      values[bufsize] = m->val[i*m->n_rows + j];
      bufsize += 1;

      if (bufsize == CS_CDO_ASSEMBLE_BUF_SIZE) {
#           pragma omp critical
        cs_matrix_assembler_values_add_g(xb->mav, bufsize,
                                         r_gids, c_gids, values);
        bufsize = 0;
      }

    } /* Loop on column */

    /* 1.b Add the contribution of potential DoFs */

    r_gids[bufsize] = i_gid;
    c_gids[bufsize] = p_gid;
    values[bufsize] = div_op[i];

    bufsize += 1;

    if (bufsize == CS_CDO_ASSEMBLE_BUF_SIZE) {
#       pragma omp critical
      cs_matrix_assembler_values_add_g(xb->mav, bufsize,
                                       r_gids, c_gids, values);
      bufsize = 0;
    }

    /* Its transposed  */

    r_gids[bufsize] = p_gid;
    c_gids[bufsize] = i_gid;
    values[bufsize] = div_op[i];

    bufsize += 1;

    if (bufsize == CS_CDO_ASSEMBLE_BUF_SIZE) {
#       pragma omp critical
      cs_matrix_assembler_values_add_g(xb->mav, bufsize,
                                       r_gids, c_gids, values);
      bufsize = 0;
    }
  } /* Loop on row (i) */

  if (bufsize > 0) {
#   pragma omp critical
    cs_matrix_assembler_values_add_g(xb->mav, bufsize,
                                     r_gids, c_gids, values);
    bufsize = 0;
  }

  /* 2. RHS assembly
   * ============================================== */

  for (short int f = 0; f < n_f; f++)
#   pragma omp atomic
    rhs[csys->dof_ids[f]] += csys->rhs[f];

  rhs[n_faces + cm->c_id] = csys->rhs[n_f];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         pointer to a cs_cdocb_scaleq_t strucutre
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_cdocb_scaleq_init_cell_system(const cs_cell_mesh_t         *cm,
                               const cs_equation_param_t    *eqp,
                               const cs_equation_builder_t  *eqb,
                               const cs_cdocb_scaleq_t      *eqc,
                               cs_cell_sys_t                *csys,
                               cs_cell_builder_t            *cb)
{
  /* Cell-wise view of the linear system to build */

  const int  n_dofs = cm->n_fc;

  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;

  /* Initialize the local system */

  cs_cell_sys_reset(cm->n_fc, csys);

  cs_sdm_square_init(n_dofs, csys->mat);

  memset(csys->rhs, 0., (n_dofs + 1)*sizeof(cs_real_t));

  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  for (short int f = 0; f < cm->n_fc; f++) {

    csys->dof_ids[f] = cm->f_ids[f];
    csys->val_n[f] = eqc->flux[cm->f_ids[f]];

  } /* Loop on cell faces */

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    cs_equation_bc_set_cw_cb(cm,
                             eqp,
                             eqb->face_bc,
                             eqb->dir_values,
                             csys,
                             cb);

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    cs_dbg_check_hmg_dirichlet_cw(__func__, csys);
#endif
  } /* Border cell */

  cs_real_t *div_op = eqc->div_op_cw[t_id];

  for (short int f = 0; f < cm->n_fc; f++)
    div_op[f] = (  csys->bf_flag[f] & CS_CDO_BC_HMG_NEUMANN
                 ||csys->bf_flag[f] & CS_CDO_BC_NEUMANN) ? 0.0: -cm->f_sgn[f];

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOCB_SCALEQ > 2
  if (cs_dbg_cw_test(NULL, cm, csys)) {
#   pragma omp critical
    {
      cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
      for (short int f = 0; f < cm->n_fc; f++)
        cs_log_printf(CS_LOG_DEFAULT, "    f%2d: %- .4e \n",
                      f, div_op[f]);
    } /* Critical section */
  }
#endif

  /* Build local arrays related to the boundary conditions */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOCB_SCALEQ_DBG > 2
  if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system in the steady-state case
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      eqb       pointer to a cs_equation_builder_t structure
 * \param[in,out]  eqc       pointer to a cs_cdocb_scaleq_t strucutre
 */
/*----------------------------------------------------------------------------*/

static void
_steady_build(const cs_equation_param_t      *eqp,
              const cs_equation_builder_t    *eqb,
              cs_cdocb_scaleq_t              *eqc)
{
  /* Retrieve shared structures */

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_shared_time_step;

#if defined(DEBUG) && !defined(NDEBUG)
  if (quant->n_b_faces > 0)
    assert(eqb->dir_values != NULL);
#endif

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_face_mesh_t  *fm = cs_cdo_local_get_face_mesh(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL:eqc->mass_hodge[t_id];

    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdocb_scaleq_get(&csys, &cb);

    /* Set times at which one evaluates quantities when needed */

    cb->t_pty_eval = ts->t_cur; /* Dummy parameter if really steady */
    cb->t_bc_eval = ts->t_cur;  /* Dummy parameter if really steady */
    cb->t_st_eval = ts->t_cur;  /* Dummy parameter if really steady */

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, mass_hodge, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* For the stationary problem, the global system writes:
       *
       *     |        |         |
       *     |   A    |    Bt   |  B is the divergence (Bt the gradient)
       *     |        |         |  A is csys->mat in what follows
       *     |--------|---------|  The viscous part arising from the CDO-Cb
       *     |        |         |  schemes for vector-valued variables and
       *     |   B    |    0    |  additional terms as the linearized
       *     |        |         |  convective term
       *
       * Set the local (i.e. cellwise) structures for the current cell
       */

      _cdocb_scaleq_init_cell_system(cm, eqp, eqb,
                                     eqc, csys, cb);


      /* 2- SCALAR EQUATION */
      /* ================================ */

      cs_cdocb_scaleq_diffusion(eqp, eqb, eqc, cm,
                                mass_hodge, csys, cb);

      /* 3- SOURCE TERM COMPUTATION (for the equation) */
      /* ====================================================== */

      if (cs_equation_param_has_sourceterm(eqp)) {

        /* Reset the local contribution */

        memset(csys->source, 0, (cm->n_fc + 1)*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */

        cs_source_term_compute_cellwise(eqp->n_source_terms,
                    (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        csys->source);

        csys->rhs[cm->n_fc] += csys->source[cm->n_fc];

      }

      /* OTHER RHS CONTRIBUTIONS
       * ===========================
       *
       * BOUNDARY CONDITIONS
       * =================== */

      _scb_apply_bc(eqp, eqc, cm, fm, mass_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOCB_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump("\n>> (FINAL) Cell system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      _full_assembly(csys, cm, eqb, eqc, asb);

    } /* Main loop on cells */

  } /* End of the OpenMP Block */

  cs_cdo_system_helper_finalize_assembly(eqb->system_helper);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current content of the related variables-fields
 *         to previous values
 *
 * \param[in, out] eqc       pointer to a context structure
 */
/*----------------------------------------------------------------------------*/

static void
_cdocb_scaleq_current_to_previous(cs_cdocb_scaleq_t       *eqc)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  if (eqc->flux_pre != NULL)
    memcpy(eqc->flux_pre, eqc->flux,
           quant->n_faces * sizeof(cs_real_t));
  if (eqc->potential_pre != NULL)
    memcpy(eqc->potential_pre, eqc->potential,
           quant->n_cells * sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence in a cell of a flux array defined at
 *         faces (values are defined both at interior and border faces).
 *         Variant based on the usage of \ref cs_cdo_quantities_t structure.
 *
 * \param[in]     cm           pointer to a cellwise view of the mesh
 * \param[in]     f_vals       values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

static double
_cdocb_cell_divergence(const cs_cell_mesh_t         *cm,
                       const cs_real_t              *f_vals)
{
  double  div = 0.0;
  for (cs_lnum_t f = 0; f < cm->n_fc; f++) {

    const cs_real_t  sgn_f = -cm->f_sgn[f];

    cs_lnum_t f_id = cm->f_ids[f];
    div += sgn_f*f_vals[f_id];
  } /* Loop on cell faces */

  div /= cm->vol_c;

  return div;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the variables related to CDO-Cb system after a resolution
 *
 * \param[in, out] tce       pointer to a timer counter
 * \param[in, out] fld       pointer to a cs_field_t structure
 * \param[in, out] eqc       pointer to a context structure
 * \param[in]      cur2prev  true if one performs "current to previous" op.
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_cell_fields(cs_timer_counter_t      *tce,
                    cs_field_t              *fld,
                    cs_cdocb_scaleq_t       *eqc,
                    bool                     cur2prev)
{
  cs_timer_t  t0 = cs_timer_time();

  cs_cdo_quantities_t *quant = cs_shared_quant;

  /* Copy current field values to previous values */

  if (cur2prev)
    cs_field_current_to_previous(fld);

  memcpy(fld->val, eqc->potential, quant->n_cells*sizeof(cs_real_t));

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(tce, &t0, &t1);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  *csys = _scb_cell_system[t_id];
  *cb = _scb_cell_builder[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current content of the related variables-fields
 *         to previous values
 *         Case of the monolithic coupling algorithm.
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out] context   pointer to a scheme context structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_current_to_previous(const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t       *eqb,
                                   void                        *context)
{
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);

  cs_cdocb_scaleq_t *eqc = (cs_cdocb_scaleq_t *)context;
  _cdocb_scaleq_current_to_previous(eqc);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion term in the
 *          scalar-valued CDO-Cb schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] mass_hodge  pointer to a cs_hodge_t structure for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_diffusion(const cs_equation_param_t     *eqp,
                          const cs_equation_builder_t   *eqb,
                          const cs_cdocb_scaleq_t       *eqc,
                          const cs_cell_mesh_t          *cm,
                          cs_hodge_t                    *mass_hodge,
                          cs_cell_sys_t                 *csys,
                          cs_cell_builder_t             *cb)
{
  if (cs_equation_param_has_diffusion(eqp)) {   /* DIFFUSION TERM
                                                 * ============== */

    /* Set the diffusion property */

    assert(mass_hodge != NULL);
    if (!(eqb->diff_pty_uniform))
      cs_hodge_set_property_value_cw(cm, cb->t_pty_eval, cb->cell_flag,
                                     mass_hodge);

    /* Define the local stiffness matrix: local matrix owned by the cellwise
       builder (store in cb->loc) */

    bool computed = eqc->get_mass_matrix(cm, mass_hodge, cb);

    /* Add the local diffusion operator to the local system */
    if(computed)
      cs_sdm_add(csys->mat, mass_hodge->matrix);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOCB_SCALEQ_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after diffusion", csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdocb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_extra_post(const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *context)
{
  CS_UNUSED(eqp);

  cs_timer_t  t0 = cs_timer_time();

  cs_cdocb_scaleq_t  *eqc = (cs_cdocb_scaleq_t *)context;

  const cs_field_t  *field = cs_field_by_id(eqc->var_field_id);
  const cs_lnum_t  n_i_faces = cs_shared_connect->n_faces[CS_INT_FACES];
  const cs_real_t  *bface_values = eqc->flux + n_i_faces;

  /* In case of postprocessing of the border faces, one has to check if there
     is a mesh modification. In particular, a removal of 2D extruded border
     faces*/

  bool  use_parent = (cs_shared_quant->remove_boundary_faces) ? false : true;

  /* Field post-processing */

  char *postlabel = NULL;
  int  len = strlen(field->name) + 8 + 1;
  BFT_MALLOC(postlabel, len, char);
  sprintf(postlabel, "%s.Border", field->name);

  cs_post_write_var(CS_POST_MESH_BOUNDARY,
                    CS_POST_WRITER_ALL_ASSOCIATED,
                    postlabel,
                    field->dim,
                    true,                  /* interlaced arrays */
                    use_parent,
                    CS_POST_TYPE_cs_real_t,
                    NULL,                  /* values on cells */
                    NULL,                  /* values at internal faces */
                    bface_values,          /* values at border faces */
                    cs_shared_time_step);  /* time step management structure */

  BFT_FREE(postlabel);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh faces for the current context.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size n_faces)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdocb_scaleq_get_face_values(void    *context,
                                bool     previous)
{
  cs_cdocb_scaleq_t  *eqc = (cs_cdocb_scaleq_t *)context;

  if (eqc == NULL)
    return NULL;

  if (previous) {
    assert(eqc->flux_pre != NULL);
    return eqc->flux_pre;
  }
  else
    return eqc->flux;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh vertices for the variable field
 *         associated to the given context
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size: n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdocb_scaleq_get_cell_values(void        *context,
                                bool         previous)
{
  cs_cdocb_scaleq_t  *eqc = (cs_cdocb_scaleq_t *)context;
  if (eqc == NULL)
    return NULL;

  cs_field_t *potential_fld = cs_field_by_id(eqc->var_field_id);

  if (previous) {
    assert(potential_fld->val_pre != NULL);
    return potential_fld->val_pre;
  }
  else
    return potential_fld->val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_init_sharing(const cs_cdo_quantities_t     *quant,
                             const cs_cdo_connect_t        *connect,
                             const cs_time_step_t          *time_step)
{
  /* Assign static const pointers */

  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  /* Need to build special range set and interfaces ? */

  int  row_size = connect->n_max_fbyc,
       col_size = connect->n_max_fbyc;

  assert(cs_glob_n_threads > 0);
  BFT_MALLOC(_scb_cell_system, cs_glob_n_threads, cs_cell_sys_t *);
  BFT_MALLOC(_scb_cell_builder, cs_glob_n_threads, cs_cell_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    _scb_cell_system[i] = NULL;
    _scb_cell_builder[i] = NULL;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    _scb_cell_system[t_id] = cs_cell_sys_create(connect->n_max_fbyc + 1,
                                                connect->n_max_fbyc,
                                                1, NULL);
    _scb_cell_builder[t_id] = _cell_builder_create(connect);
  }
#else
  assert(cs_glob_n_threads == 1);
  _scb_cell_system[0] = cs_cell_sys_create(connect->n_max_fbyc + 1,
                                           connect->n_max_fbyc,
                                           1, NULL);
  _scb_cell_builder[0] = _cell_builder_create(connect);

#endif /* openMP */

  /* SLES needs these structures for advanced PETSc hooks */

  cs_cdocb_monolithic_sles_init_sharing(connect, quant);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free shared pointers with lifecycle dedicated to this file
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_finalize_sharing(void)
{
  /* Need to build special range set and interfaces ? */

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_cell_sys_free(&(_scb_cell_system[t_id]));
    cs_cell_builder_free(&(_scb_cell_builder[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);
  cs_cell_sys_free(&(_scb_cell_system[0]));
  cs_cell_builder_free(&(_scb_cell_builder[0]));
#endif /* openMP */

  BFT_FREE(_scb_cell_system);
  BFT_FREE(_scb_cell_builder);
  _scb_cell_system = NULL;
  _scb_cell_builder = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdocb_scaleq_t structure storing data useful
 *         for building and managing such a scheme
 *
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id      id of the variable field
 * \param[in]      bflux_id    id of the boundary flux field
 * \param[in, out] eqb         pointer to a \ref cs_equation_builder_t struct.
 *
 * \return a pointer to a new allocated cs_cdocb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdocb_scaleq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
  assert(eqp != NULL && eqb != NULL);

  if (eqp->space_scheme != CS_SPACE_SCHEME_CDOCB && eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid type of equation.\n"
              " Expected: scalar-valued CDO Cell-based equation.", __func__);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t *quant = cs_shared_quant;
  const cs_lnum_t  n_cells = connect->n_cells;
  const cs_lnum_t  n_faces = connect->n_faces[CS_ALL_FACES];

  const cs_domain_t *domain = cs_glob_domain;

  cs_cdocb_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_cdocb_scaleq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  /* Dimensions of the algebraic system */

  eqc->n_faces = n_faces;
  eqc->n_dofs = n_faces;

  /* Flag to indicate the minimal set of quantities to build in a cell mesh
     According to the situation, additional flags have to be set */
  eqb->msh_flag = CS_FLAG_COMP_PF | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ;

  BFT_MALLOC(eqc->flux, n_faces, cs_real_t);
  BFT_MALLOC(eqc->flux_pre, n_faces, cs_real_t);
  BFT_MALLOC(eqc->potential, n_cells, cs_real_t);
  BFT_MALLOC(eqc->potential_pre, n_cells, cs_real_t);
  BFT_MALLOC(eqc->rc_tilda, n_cells, cs_real_t);

# pragma omp parallel if (n_cells > CS_THR_MIN)
  {
    /* Values at each face (interior and border) i.e. take into account BCs */

#   pragma omp for nowait
    for (cs_lnum_t i = 0; i < n_faces; i++) eqc->flux[i] = 0;

#   pragma omp for nowait
    for (cs_lnum_t i = 0; i < n_faces; i++) eqc->flux_pre[i] = 0;

#   pragma omp for nowait
    for (cs_lnum_t i = 0; i < n_cells; i++) {
      eqc->potential[i] = 0.;
      eqc->potential_pre[i] = 0.;
    }

    /* Store the last computed values of the field at cell centers and the data
       needed to compute the cell values from the face values.
       No need to synchronize all these quantities since they are only cellwise
       quantities. */

#   pragma omp for
    for (cs_lnum_t i = 0; i < n_cells; i++) eqc->rc_tilda[i] = 0;
  }

  /* Assume the matrix is diagonal */

  BFT_MALLOC(eqc->acf_tilda, connect->c2f->idx[n_cells], cs_real_t);
  memset(eqc->acf_tilda, 0, connect->c2f->idx[n_cells]*sizeof(cs_real_t));

  bool  need_eigen =
    (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
     eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) ? true : false;

  /* Diffusion term */

  /* For Cell-based scheme, no diffusion hodge
   * operator is needed, only mass hodge operator
   * is used in a coupled Saddle-Point system. */
  /* ----------------------------------------- */

  eqc->get_stiffness_matrix = NULL;
  eqc->diffusion_hodge = NULL;

  eqc->get_mass_matrix = NULL;
  eqc->mass_hodge = NULL;

  if (cs_equation_param_has_diffusion(eqp)) {

    /* Mass matrix */

    eqc->mass_hodgep.inv_pty = true;
    eqc->mass_hodgep.type = CS_HODGE_TYPE_FPED;
    eqc->mass_hodgep.algo = CS_HODGE_ALGO_COST;
    eqc->mass_hodgep.coef = cs_math_1ov3;

    eqb->sys_flag |= CS_FLAG_SYS_MASS_MATRIX;

    eqc->get_mass_matrix = cs_hodge_fped_cost_get;
    eqc->mass_hodge = cs_hodge_init_context(connect,
                                            eqp->diffusion_property,
                                            &(eqc->mass_hodgep),
                                            true,  /* tensor ? */
                                            need_eigen); /* eigen ? */

    if (eqp->verbosity > 1) {
      cs_log_printf(CS_LOG_SETUP,
                    "#### Parameters of the mass matrix of the equation %s\n",
                    eqp->name);
      cs_hodge_param_log("Mass matrix", NULL, eqc->mass_hodgep);
    }

  } /* Diffusion term is requested */

  /* Boundary conditions */
  /* ------------------- */

  eqc->enforce_robin_bc = NULL;

  switch (eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    eqc->enforce_neumann = cs_cdocb_neumann;
    eqc->enforce_dirichlet = cs_cdocb_dirichlet;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  /* No advection term is requested */

  if (eqp->default_enforcement != CS_PARAM_BC_ENFORCE_WEAK_NITSCHE)
    eqb->sys_flag |= CS_FLAG_SYS_SYM; /* Algebraic system is symmetric */

  /* No Reaction term is requested */

  /* Source term */
  /* ----------- */

  eqc->source_terms = NULL;

  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, n_cells, cs_real_t);
    memset(eqc->source_terms, 0, n_cells*sizeof(cs_real_t));

  } /* There is at least one source term */

  /* Helper structures (range set, interface set, matrix structure and all the
     assembly process) */

  cs_cdo_system_helper_t  *sh = NULL;

  cs_lnum_t  block_size = quant->n_faces + quant->n_cells;

  sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                   1,
                                   &block_size,
                                   1);

  cs_cdocb_monolithic_sles_t  *msles =
    cs_cdocb_monolithic_sles_create(quant->n_faces, quant->n_cells);

  eqc->msles = msles;

  /* Choose the right class of matrix to avoid copy.
   * The way to perform the assembly may change if an external librairy is used
   * for solving the linear system */

  cs_cdo_system_block_t *a =
    cs_cdo_system_add_xblock(sh, 0,        /* block id */
                             block_size);  /* n_dofs */

  _build_shared_full_structures(a, true);

  cs_cdo_system_build_block(sh, 0); /* build/set structures */

  eqb->system_helper = sh;

  /* Renormalization of the residual */

  if (eqp->sles_param->resnorm_type == CS_PARAM_RESNORM_WEIGHTED_RHS)
    eqb->msh_flag |= CS_FLAG_COMP_PFC;

  BFT_MALLOC(eqc->div_op_cw, cs_glob_n_threads, cs_real_t*);

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    BFT_MALLOC(eqc->div_op_cw[t_id], connect->n_max_fbyc, cs_real_t);
  }
#else

  assert(cs_glob_n_threads == 1);

  BFT_MALLOC(eqc->div_op_cw[0], connect->n_max_fbyc, cs_real_t);

#endif /* openMP */

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdocb_scaleq_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdocb_scaleq_free_context(void       *scheme_context)
{
  cs_cdocb_scaleq_t  *eqc = (cs_cdocb_scaleq_t *)scheme_context;

  if (eqc == NULL)
    return eqc;

  /* These arrays may have not been allocated */

  BFT_FREE(eqc->source_terms);

  BFT_FREE(eqc->flux);
  if (eqc->flux_pre != NULL)
    BFT_FREE(eqc->flux_pre);

  BFT_FREE(eqc->potential);
  if (eqc->potential_pre != NULL)
    BFT_FREE(eqc->potential_pre);

  BFT_FREE(eqc->rc_tilda);
  BFT_FREE(eqc->acf_tilda);

  cs_hodge_free_context(&(eqc->diffusion_hodge));
  cs_hodge_free_context(&(eqc->mass_hodge));

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    BFT_FREE(eqc->div_op_cw[t_id]);
  }
#else
  assert(cs_glob_n_threads == 1);
    BFT_FREE(eqc->div_op_cw[0]);
#endif /* openMP */

  BFT_FREE(eqc->div_op_cw);

  /* Free the context structure for solving saddle-point system */

  cs_cdocb_monolithic_sles_free(&(eqc->msles));

  /* Other pointers are only shared (i.e. not owner) */

  BFT_FREE(eqc);

  return;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *         Define an indirection array for the enforcement of internal DoFs
 *         only if needed.
 *         Case of scalar-valued CDO-Cb schemes
 *
 * \param[in]      t_eval          time at which one evaluates BCs
 * \param[in]      mesh            pointer to a cs_mesh_t structure
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in, out] eqb             pointer to a cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_setup(cs_real_t                      t_eval,
                      const cs_mesh_t               *mesh,
                      const cs_equation_param_t     *eqp,
                      cs_equation_builder_t         *eqb)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Compute the values of the Dirichlet BC */

  BFT_MALLOC(eqb->dir_values, quant->n_b_faces, cs_real_t);
  memset(eqb->dir_values, 0, quant->n_b_faces*sizeof(cs_real_t));

  cs_equation_compute_dirichlet_cb(mesh,
                                   quant,
                                   connect,
                                   eqp,
                                   eqb->face_bc,
                                   t_eval,
                                   _scb_cell_builder[0], /* static variable */
                                   eqb->dir_values);

  /* Internal enforcement of DoFs */

  if (cs_equation_param_has_internal_enforcement(eqp))
    eqb->enforced_values =
      cs_enforcement_define_at_faces(connect,
                                     eqp->n_enforcements,
                                     eqp->enforcement_params);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-Cb schemes.
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
cs_cdocb_scaleq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdocb_scaleq_t  *eqc = (cs_cdocb_scaleq_t *)context;
  cs_real_t  *potential_vals = eqc->potential;
  cs_real_t  *flux_vals = eqc->flux;

  /* Check that a face interface has been defined */

  if (eqp->n_ic_defs > 0 && cs_glob_n_ranks > 1 && connect->face_ifs == NULL)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Interface set structure at faces not allocated.\n",
              __func__);

  /* By default, 0 is set as initial condition for the computational domain */

  memset(flux_vals, 0, quant->n_faces*sizeof(cs_real_t));
  memset(potential_vals, 0, quant->n_cells*sizeof(cs_real_t));

  if (eqp->n_ic_defs > 0) {

    cs_lnum_t  *def2f_ids = (cs_lnum_t *)cs_cdo_toolbox_get_tmpbuf();
    cs_lnum_t  *def2f_idx = NULL;

    BFT_MALLOC(def2f_idx, eqp->n_ic_defs + 1, cs_lnum_t);

    cs_cdo_sync_vol_def_at_faces(eqp->n_ic_defs, eqp->ic_defs,
                                 def2f_idx,
                                 def2f_ids);

    for (int def_id = 0; def_id < eqp->n_ic_defs; def_id++) {

      /* Get and then set the definition of the initial condition */

      const cs_xdef_t  *def = eqp->ic_defs[def_id];
      const cs_lnum_t  n_f_selected = def2f_idx[def_id+1] - def2f_idx[def_id];
      const cs_lnum_t  *selected_lst = def2f_ids + def2f_idx[def_id];

      switch(def->type) {

      case CS_XDEF_BY_VALUE:
        cs_evaluate_potential_at_cells_by_value(def, potential_vals);
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid way to initialize field values for eq. %s.\n",
                  __func__, eqp->name);

      } /* Switch on possible type of definition */

    } /* Loop on definitions */

    BFT_FREE(def2f_idx);

  } /* Initial values to set */

  if (eqc->flux_pre != NULL)
    memcpy(eqc->flux_pre, eqc->flux,
           quant->n_faces*sizeof(cs_real_t));

  if (eqc->potential_pre != NULL)
    memcpy(eqc->potential_pre, eqc->potential,
           quant->n_cells*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain
 *         Case of scalar-valued CDO Cell-based scheme
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context   pointer to a scheme builder structure
 *
 * \return a pointer to a \ref cs_cdo_balance_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_balance_t *
cs_cdocb_scaleq_balance(const cs_equation_param_t     *eqp,
                        cs_equation_builder_t         *eqb,
                        void                          *context)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const char *func_name = __func__;

  cs_timer_t  t0 = cs_timer_time();

  cs_cdocb_scaleq_t  *eqc = (cs_cdocb_scaleq_t *)context;

  /* Allocate and initialize the structure storing the balance evaluation */

  cs_cdo_balance_t  *eb = cs_cdo_balance_create(cs_flag_primal_cell,
                                                quant->n_cells);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(quant, connect, ts, eqp, eqb, eqc, eb, _scb_cell_builder, \
         func_name)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = _scb_cell_builder[t_id];
    cs_hodge_t  *mass_hodge =
      (eqc->mass_hodge == NULL) ? NULL : eqc->mass_hodge[t_id];

    /* Set times at which one evaluates quantities when needed */

    cb->t_pty_eval = ts->t_cur; /* Dummy parameter if really steady */
    cb->t_bc_eval = ts->t_cur;  /* Dummy parameter if really steady */
    cb->t_st_eval = ts->t_cur;  /* Dummy parameter if really steady */

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(eqp, eqb, mass_hodge, cb);

    /* Set inside the OMP section so that each thread has its own value */

    cs_real_t *flux, *flux_pre;

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag, eqb),
                         connect, quant, cm);

      /* Set the value of the current potential */

      flux = eqc->flux;

      if (eqb->sys_flag & CS_FLAG_SYS_MASS_MATRIX) { /* MASS MATRIX
                                                      * =========== */
        assert(mass_hodge != NULL);

        /* Build the mass matrix and store it in mass_hodge->matrix */

        eqc->get_mass_matrix(cm, mass_hodge, cb);

      }

      /* Diffusion term */

      if (cs_equation_param_has_diffusion(eqp)) {

        cs_real_t div = 0.;
        div = _cdocb_cell_divergence(cm, flux);

        eb->diffusion_term[cm->c_id] += div;

      } /* End of diffusion */

      /* Source term */

      if (cs_equation_param_has_sourceterm(eqp)) {

        /* Reset the local contribution */

        cs_real_t  *src = cb->values;
        memset(src, 0, (cm->n_fc + 1)*sizeof(cs_real_t));

        /* Source term contribution to the algebraic system
           If the equation is steady, the source term has already been computed
           and is added to the right-hand side during its initialization. */

        cs_source_term_compute_cellwise(eqp->n_source_terms,
                                        (cs_xdef_t *const *)eqp->source_terms,
                                        cm,
                                        eqb->source_mask,
                                        eqb->compute_source,
                                        cb->t_st_eval,
                                        mass_hodge,
                                        cb,
                                        src);

        eb->source_term[cm->c_id] += src[cm->n_fc];

      } /* End of term source contribution */

    } /* Main loop on cells */

  } /* OPENMP Block */

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
 * \brief  Solve the steady Stokes or Oseen system with a CDO cell-based scheme
 *         using a monolithic approach and GKB algorithm
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_solve_steady_state(bool                        cur2prev,
                                   const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
  cs_timer_t  t0 = cs_timer_time();
  /* Retrieve high-level structures */

  cs_cdocb_scaleq_t  *eqc = (cs_cdocb_scaleq_t *) context;
  cs_cdo_system_helper_t  *sh = eqb->system_helper;

  cs_field_t *potential_fld = cs_field_by_id(eqc->var_field_id);

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;

  /* Build an array storing the Dirichlet values at faces */

  cs_cdocb_scaleq_setup(t_cur, mesh, eqp, eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t  *rhs = NULL;  /* Since it is NULL, sh get sthe ownership */

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Main loop on cells to define the linear system to solve */

  _steady_build(eqp, eqb, eqc);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);

  /* Free temporary buffers and structures */

  cs_equation_builder_reset(eqb);

  /* End of the system building */

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */

  _cdocb_scaleq_current_to_previous(eqc);

  /* Solve the linear system */

  cs_timer_t  t_solve_start = cs_timer_time();
  cs_cdocb_monolithic_sles_t  *msles = eqc->msles;

  /* matrix has been already assigned to the msles structure */

  msles->flux = eqc->flux; /* FLux DoFs at faces */
  msles->potential = eqc->potential;    /* pressure DoFs at cells */

  msles->sles = cs_sles_find_or_add(eqp->sles_param->field_id, NULL);

  int  cumulated_inner_iters = cs_cdocb_monolithic_solve(eqp, sh, msles);

  _update_cell_fields(&(eqb->tcs), potential_fld, eqc, false);

  cs_cdo_bc_face_t *face_bc = eqb->face_bc;
  if (face_bc->n_hmg_dir_faces + face_bc->n_nhmg_dir_faces == 0) {
    cs_field_set_volume_average(potential_fld, 0.);
    memcpy(eqc->potential, potential_fld->val, quant->n_cells*sizeof(cs_real_t));
  }

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcs), &t_solve_start, &t_solve_end);

  cs_log_printf(CS_LOG_DEFAULT, " -cvg- CDO Cb: cumulated_inner_iters: %d\n",
      cumulated_inner_iters);
  cs_log_printf_flush(CS_LOG_DEFAULT);

  /* Frees */
  cs_cdocb_monolithic_sles_clean(msles);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
