/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it as one block (monolithic approach of the velocity-pressure
 * coupling)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#if defined(HAVE_PETSC)
#include <petscversion.h>
#include <petscksp.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_blas.h"
#include "cs_cdo_bc.h"
#include "cs_cdofb_priv.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_cdofb_navsto.h"
#include "cs_cdofb_monolithic_sles.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_navsto_coupling.h"
#include "cs_parall.h"
#include "cs_param.h"
#include "cs_post.h"
#include "cs_sdm.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_monolithic_priv.h"
#include "cs_cdofb_monolithic.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdofb_monolithic.c
 *
 * \brief Build an algebraic CDO face-based system for the Navier-Stokes
 * equations and solved it with a monolithic approach
 */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_MONOLITHIC_DBG      0

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_mesh_t              *cs_shared_mesh;
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;

/* If one solves the saddle-point system as it is, then one needs to define
 * new shared structure. Otherwise, one points to shared structures with
 * vector-valued face-based schemes */
static cs_range_set_t         *_shared_range_set = NULL;
static cs_interface_set_t     *_shared_interface_set = NULL;
static cs_matrix_structure_t  *_shared_matrix_structure = NULL;
static cs_matrix_assembler_t  *_shared_matrix_assembler = NULL;

static const cs_range_set_t         *cs_shared_range_set;
static const cs_interface_set_t     *cs_shared_interface_set;
static const cs_matrix_structure_t  *cs_shared_matrix_structure;
static const cs_matrix_assembler_t  *cs_shared_matrix_assembler;

static cs_sdm_t **cs_cdofb_monolithic_cw_mat = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current content of the Navier-Stokes related variables-fields
 *         to previous values
 *
 * \param[in, out]       sc     pointer to \ref cs_cdofb_monolithic_t structure
 */
/*----------------------------------------------------------------------------*/

static inline void
_fields_to_previous(cs_cdofb_monolithic_t   *sc)
{
  assert(sc != NULL);
  cs_field_current_to_previous(sc->velocity);
  cs_field_current_to_previous(sc->pressure);
  cs_field_current_to_previous(sc->divergence);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Test if one has to do one more Picard iteration
 *
 * \param[in]      nsp           pointer to a \ref cs_navsto_param_t structure
 * \param[in]      mom_eqc       scheme context for the momentum equation
 * \param[in]      div_l2_norm   L2 norm of the velocity divergence
 * \param[in, out] ns_info       pointer to a cs_navsto_algo_info_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_picard_cvg_test(const cs_navsto_param_t      *nsp,
                 const cs_cdofb_vecteq_t      *mom_eqc,
                 cs_real_t                     div_l2_norm,
                 cs_navsto_algo_info_t        *ns_info)
{
  const cs_real_t  diverg_factor = 100;
  const cs_navsto_param_sles_t  nslesp = nsp->sles_param;

  cs_real_t  previous_picard_res = ns_info->res;

  /* Increment the number of Picard iterations */
  ns_info->n_algo_iter += 1;

  /* Compute the new residual: L2 norm of the face velocity increment */
  ns_info->res = cs_evaluate_delta_3_square_wc2x_norm(mom_eqc->face_values,
                                                      mom_eqc->face_values_pre,
                                                      cs_shared_connect->c2f,
                                                      cs_shared_quant->pvol_fc);
  assert(ns_info->res > -DBL_MIN);
  ns_info->res = sqrt(ns_info->res);

  /* Set the convergence status */
  if (ns_info->res < nslesp.picard_tolerance &&
      div_l2_norm < nslesp.picard_tolerance) {
    ns_info->cvg = CS_SLES_CONVERGED;
    return;
  }

  if (ns_info->res >  diverg_factor * previous_picard_res) {
    ns_info->cvg = CS_SLES_DIVERGED;
    return;
  }

  if (ns_info->n_algo_iter >= nslesp.picard_n_max_iter) {
    ns_info->cvg = CS_SLES_MAX_ITERATION;
    return;
  }

  ns_info->cvg = CS_SLES_ITERATING;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the cs_range_set_t, cs_interface_set_t, cs_matrix_assembler_t
 *         and cs_matrix_structure_t structures
 */
/*----------------------------------------------------------------------------*/

static void
_build_shared_structures(void)
{
  /* Compute the range set for an array of size 3*n_faces + n_cells
   * velocity is attached to faces (one for each component) and pressure
   * to cells
   *
   * Storage for the global numbering: Vel_X | Vel_Y | Vel_Z | Pressure
   */

  const cs_mesh_t  *m = cs_shared_mesh;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_lnum_t  n_i_faces = m->n_i_faces;
  const cs_lnum_t  n_b_faces = m->n_b_faces;
  const cs_lnum_t  size = 3*n_faces + m->n_cells;
  const cs_gnum_t  n_g_faces = m->n_g_i_faces + m->n_g_b_faces;

  /* 1. Build the global numbering */
  cs_gnum_t  *gnum = NULL;
  BFT_MALLOC(gnum, size, cs_gnum_t);

  if (cs_glob_n_ranks > 1) {

    for (int xyz = 0; xyz < 3; xyz++) {

      cs_gnum_t  *_ignum = gnum + xyz * n_faces;
      cs_gnum_t  *_bgnum = _ignum + n_i_faces;
      const cs_gnum_t  igshift = xyz * n_g_faces;
      const cs_gnum_t  bgshift = igshift + m->n_g_i_faces;

#     pragma omp parallel if (n_i_faces > CS_THR_MIN)
      {
        /* Interior faces (X, Y or Z) */
#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < n_i_faces; i++)
          _ignum[i] = igshift + m->global_i_face_num[i];

        /* Boundary faces (X, Y or Z) */
#       pragma omp for nowait
        for (cs_lnum_t i = 0; i < n_b_faces; i++)
          _bgnum[i] = bgshift + m->global_b_face_num[i];

      } /* End of the OpenMP region */

    } /* Loop on components */

    /* Add pressure DoFs */
    cs_gnum_t  *_pgnum = gnum + 3*n_faces;
    const cs_gnum_t  pgshift = 3*n_g_faces;

#   pragma omp parallel if (n_i_faces > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < m->n_cells; i++)
      _pgnum[i] = m->global_cell_num[i] + pgshift;

  }
  else {

#   pragma omp parallel for if (size > CS_THR_MIN)
    for (cs_gnum_t i = 0; i < (cs_gnum_t)size; i++)
      gnum[i] = i + 1;

  } /* Sequential or parallel run */

  /* 2. Build the interface set and the range set structures */

  /* Do not consider periodicity up to now. Should split the face interface
     into interior and border faces to do this, since only boundary faces
     can be associated to a periodicity */

  _shared_interface_set = cs_interface_set_create(size,
                                                  NULL,
                                                  gnum,
                                                  m->periodicity,
                                                  0, NULL, NULL, NULL);

  _shared_range_set = cs_range_set_create(_shared_interface_set,
                                          NULL,      /* halo */
                                          size,
                                          false,     /* TODO: Ask Yvan */
                                          0);        /* g_id_base */

  /* Free memory */
  BFT_FREE(gnum);

  /* 3. Build the matrix assembler structure */
  const cs_adjacency_t  *f2f = cs_shared_connect->f2f;
  const cs_adjacency_t  *f2c = cs_shared_connect->f2c;

  /* The second paramter is set to "true" meaning that the diagonal is stored
   * separately --> MSR storage
   * Create the matrix assembler structure
   */
  _shared_matrix_assembler =
    cs_matrix_assembler_create(_shared_range_set->l_range, true);

  /* First loop to count max size of the buffer used to fill the matrix
   * structure. +1 to take into account the diagonal term.
   */
  int  max_sten = 0;
  for (cs_lnum_t f = 0; f < n_faces; f++) {
    int  sten
      = 9*(f2f->idx[f+1]-f2f->idx[f] + 1) + 6*(f2c->idx[f+1]-f2c->idx[f]);
    max_sten = CS_MAX(max_sten, sten);
  }

  cs_gnum_t  *grows = NULL, *gcols = NULL;
  BFT_MALLOC(grows, max_sten, cs_gnum_t);
  BFT_MALLOC(gcols, max_sten, cs_gnum_t);

  /*
   *   | A_xx  |       |       | Bt_x  |
   *   |-------|-------|-------|-------|
   *   |       | A_yy  |       | Bt_y  |
   *   |-------|-------|-------|-------|
   *   |       |       | A_zz  | Bt_z  |
   *   |-------|-------|-------|-------|
   *   | B_x   | B_y   | B_z   |  0    |
   *
   *  Each block A_.. is n_faces * n_faces
   *  Each block B_.  is n_cells * n_faces
   */

  /* Only on faces (B_x is build in the same time as Bt_x for pressure DoFs) */
  for (cs_lnum_t frow_id = 0; frow_id < n_faces; frow_id++) {

    const cs_lnum_t  start = f2f->idx[frow_id];
    const cs_lnum_t  end = f2f->idx[frow_id+1];
    const int  n_entries /* for the velocity and the pressure DoFs*/
      = (end-start + 1)*9 + 6*(f2c->idx[frow_id+1]-f2c->idx[frow_id]);

    const cs_gnum_t  grow_ids[3]
      = {_shared_range_set->g_id[frow_id],              /* x-component */
         _shared_range_set->g_id[frow_id +   n_faces],  /* y-component */
         _shared_range_set->g_id[frow_id + 2*n_faces]}; /* z-component */

    int shift = 0;

    /* Diagonal term is excluded in this connectivity. Add it "manually" */
    for (int i = 0; i < 3; i++) {
      const cs_gnum_t  grow_id = grow_ids[i];
      for (int j = 0; j < 3; j++) {
        grows[shift] = grow_id;
        gcols[shift] = grow_ids[j];
        shift++;
      }
    }

    /* Extra diagonal couples */
    for (cs_lnum_t idx = start; idx < end; idx++) {

      const cs_lnum_t  fcol_id = f2f->ids[idx];
      const cs_gnum_t  gcol_ids[3]
        = {_shared_range_set->g_id[fcol_id],              /* x-component */
           _shared_range_set->g_id[fcol_id + n_faces],    /* y-component */
           _shared_range_set->g_id[fcol_id + 2*n_faces]}; /* z-component */

      for (int i = 0; i < 3; i++) {
        const cs_gnum_t  grow_id = grow_ids[i];
        for (int j = 0; j < 3; j++) {
          grows[shift] = grow_id;
          gcols[shift] = gcol_ids[j];
          shift++;
        }
      }

    } /* Loop on extra-diag. entries */

    /* Loop on pressure-related  entries */
    for (cs_lnum_t idx = f2c->idx[frow_id]; idx < f2c->idx[frow_id+1]; idx++) {

      const cs_lnum_t  ccol_id = f2c->ids[idx];
      const cs_gnum_t  gcol_id = _shared_range_set->g_id[3*n_faces + ccol_id];

      for (int i = 0; i < 3; i++) { /* x,y,z-component */

        grows[shift] = grow_ids[i];
        gcols[shift] = gcol_id;
        shift++;

        /* Its transposed B_x, B_y, B_z */
        grows[shift] = gcol_id;
        gcols[shift] = grow_ids[i];
        shift++;

      }

    } /* Loop on pressure related DoFs */

    cs_matrix_assembler_add_g_ids(_shared_matrix_assembler,
                                  n_entries, grows, gcols);
    assert(shift == n_entries);

  } /* Loop on face entities */

  /* 4. Build the matrix structure */
  cs_matrix_assembler_compute(_shared_matrix_assembler);

  _shared_matrix_structure =
    cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR,
                                              _shared_matrix_assembler);

  /* Free temporary buffers */
  BFT_FREE(grows);
  BFT_FREE(gcols);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the local system when this should
 *          be done before the static condensation
 *
 * \param[in]      sc          pointer to a cs_cdofb_monolithic_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      bf_type     type of boundary for the boundary face
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_apply_bc_partly(const cs_cdofb_monolithic_t   *sc,
                 const cs_equation_param_t     *eqp,
                 const cs_cell_mesh_t          *cm,
                 const cs_boundary_type_t      *bf_type,
                 cs_cell_sys_t                 *csys,
                 cs_cell_builder_t             *cb)
{
  assert(cs_equation_param_has_diffusion(eqp) == true);

  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Update the velocity-block and the right-hand side (part related to the
     * momentum equation). */

    /* Neumann boundary conditions */
    if (csys->has_nhmg_neumann)
      for (short int f  = 0; f < 3*cm->n_fc; f++)
        csys->rhs[f] += csys->neu_values[f];

    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */
      const short int  f = csys->_f_ids[i];

      if (bf_type[i] & CS_BOUNDARY_IMPOSED_VEL) {
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          sc->apply_velocity_inlet(f, eqp, cm, cb, csys);
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_WALL) {
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          if (bf_type[i] & CS_BOUNDARY_SLIDING_WALL)
            sc->apply_sliding_wall(f, eqp, cm, cb, csys);
          else
            sc->apply_fixed_wall(f, eqp, cm, cb, csys);
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_SYMMETRY) {
        /* Always weakly enforce the symmetric constraint on the
           velocity-block */
        sc->apply_symmetry(f, eqp, cm, cb, csys);
      }

      /* default: nothing to do (case of a "natural" outlet) */

    } /* Loop on boundary faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump(">> Local system matrix after partial BC enforcement",
                       csys);
#endif
  } /* Boundary cell */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the local system when this should
 *          be done after the static condensation
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] sc          pointer to a cs_cdofb_monolithic_t structure
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 * \param[in, out] nsb         builder structure for the NavSto system
 */
/*----------------------------------------------------------------------------*/

static void
_apply_remaining_bc(const cs_equation_param_t     *eqp,
                    const cs_cell_mesh_t          *cm,
                    cs_cdofb_monolithic_t         *sc,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb,
                    cs_cdofb_navsto_builder_t     *nsb)
{
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    const cs_boundary_type_t  *bf_type = nsb->bf_type;
    cs_real_t  *div_op = nsb->div_op;
    cs_real_t  *mass_rhs = sc->msles->b_c;

    /* Update the divergence operator and the right-hand side related to the
     * mass equation.
     * Enforcement of Dirichlet BC in a stronger way if this is the choice
     */
    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */
      const short int  f = csys->_f_ids[i];

      if (bf_type[i] & CS_BOUNDARY_IMPOSED_VEL) {

        /* Update mass RHS from the knowledge of the mass flux
         * TODO: multiply by rho */
        mass_rhs[cm->c_id] -= cs_math_3_dot_product(csys->dir_values + 3*f,
                                                    div_op + 3*f);

        /* Strong enforcement of u.n (--> dp/dn = 0) on the divergence */
        for (int k = 0; k < 3; k++) div_op[3*f+k] = 0;

        /* Enforcement of the velocity DoFs for the velocity-block
         * Dirichlet BCs on the three components of the velocity field.
         * If a weak enforcement is requested, this was done in a previous
         * step (before the static condensation) */
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC)
          sc->apply_velocity_inlet(f, eqp, cm, cb, csys);

      }
      else if (bf_type[i] & CS_BOUNDARY_IMPOSED_P) {

        /* Close the definition of the pressure gradient for this face */
        for (int k = 0; k < 3; k++)
          csys->rhs[3*f+k] += div_op[3*f+k] * nsb->pressure_bc_val[i];
        break;

      }
      else if (bf_type[i] & CS_BOUNDARY_WALL) {

        /* No need to update the mass RHS since there is no mass flux */

        /* Strong enforcement of u.n (--> dp/dn = 0) on the divergence */
        for (int k = 0; k < 3; k++) div_op[3*f+k] = 0;

        /* Enforcement of the velocity for the velocity-block
         * Dirichlet on the three components of the velocity field */
        if (eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
            eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {
          if (bf_type[i] & CS_BOUNDARY_SLIDING_WALL)
            sc->apply_sliding_wall(f, eqp, cm, cb, csys);
          else
            sc->apply_fixed_wall(f, eqp, cm, cb, csys);
        }

      }
      else if (bf_type[i] & CS_BOUNDARY_SYMMETRY) {

        /* No need to update the mass RHS since there is no mass flux */

        /* Weak-enforcement for the velocity-block (cf. _apply_bc_partly) */

        /* Strong enforcement of u.n (--> dp/dn = 0) on the divergence */
        for (int k = 0; k < 3; k++) div_op[3*f+k] = 0;

      }

      /* default: nothing to do (case of a "natural" outlet) */

    } /* Loop over boundary faces */

  } /* This is a boundary cell */

  /* Internal enforcement of DoFs: Update csys (matrix and rhs) */
  if (csys->has_internal_enforcement) {

    cs_equation_enforced_internal_block_dofs(eqp, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement",
                       csys);
#endif
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a matrix and its related structures needed during the
 *         assembly step.
 *         Default case.
 *
 * \param[in, out]  sc           pointer to the scheme context structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_system_default(cs_cdofb_monolithic_t        *sc)
{
  assert(sc != NULL); /* sanity check */

  /* Initialize the local matrix */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_matrix_structure);

  sc->msles->block_matrices[0] = matrix;

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

  sc->mav_structures[0] = mav;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a matrix and its related structures needed during the
 *         assembly step.
 *         Case of component-wise blocks.
 *
 * \param[in, out]  sc           pointer to the scheme context structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_system_by_blocks(cs_cdofb_monolithic_t        *sc)
{
  assert(sc != NULL); /* sanity check */

  cs_cdofb_monolithic_sles_t  *msles = sc->msles;

  /* Initialize the local matrix */
  for (int i = 0; i < msles->n_row_blocks*msles->n_row_blocks; i++) {

    cs_matrix_t  *matrix = cs_matrix_create(cs_shared_matrix_structure);

    sc->msles->block_matrices[i] = matrix;

    /* Initialize the structure to assemble values */
    cs_matrix_assembler_values_t  *mav =
      cs_matrix_assembler_values_init(matrix, NULL, NULL);

    sc->mav_structures[i] = mav;

  } /* Loop on blocks */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb schemes
 *         Shares similarities with cs_equation_assemble_block_matrix()
 *
 * \param[in]       csys              pointer to a cs_cell_sys_t structure
 * \param[in]       cm                pointer to a cs_cell_mesh_t structure
 * \param[in]       div_op            array with the divergence op. values
 * \param[in]       has_sourceterm    has the equation a source term?
 * \param[in, out]  sc                pointer to scheme context structure
 * \param[in, out]  eqc               context structure for a vector-valued Fb
 * \param[in, out]  eqa               pointer to cs_equation_assemble_t
 */
/*----------------------------------------------------------------------------*/

static void
_assembly_by_blocks(const cs_cell_sys_t            *csys,
                    const cs_cell_mesh_t           *cm,
                    const cs_real_t                *div_op,
                    const bool                      has_sourceterm,
                    cs_cdofb_monolithic_t          *sc,
                    cs_cdofb_vecteq_t              *eqc,
                    cs_equation_assemble_t         *eqa)
{
  const cs_range_set_t  *rs = cs_shared_range_set;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  int  t_id = 0;
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
#endif

  cs_cdofb_monolithic_sles_t  *msles = sc->msles;
  cs_sdm_t  *cw_mat = cs_cdofb_monolithic_cw_mat[t_id];

  /* Convert csys->mat into a block view by component */
  cs_sdm_block_33_to_xyz(csys->mat, cw_mat);

  /* Assembly is performed block by block */
  for (int i = 0; i < msles->n_row_blocks; i++) {
    for (int j = 0; j < msles->n_row_blocks; j++) {

      cs_matrix_assembler_values_t   *mav = sc->mav_structures[3*i+j];
      cs_sdm_t  *m_ij = cs_sdm_get_block(cw_mat, i, j);

      sc->elemental_assembly(m_ij, cm->f_ids, rs, eqa, mav);

    } /* Loop on blocks (j) */
  } /* Loop on blocks (i) */

  /* RHS assembly */
  cs_real_t  *rhs_x = msles->b_f, *rhs_y = rhs_x + msles->n_faces;
  cs_real_t  *rhs_z = rhs_y + msles->n_faces;

# pragma omp critical
  {
    for (short int f = 0; f < cm->n_fc; f++) {
      rhs_x[cm->f_ids[f]] += csys->rhs[3*f];
      rhs_y[cm->f_ids[f]] += csys->rhs[3*f+1];
      rhs_z[cm->f_ids[f]] += csys->rhs[3*f+2];
    }
  }

  /* Reset the value of the source term for the cell DoF
     Source term is only hold by the cell DoF in face-based schemes */
  if (has_sourceterm) {
    cs_real_t  *st = eqc->source_terms + 3*cm->c_id;
    for (int k = 0; k < 3; k++)
      st[k] = csys->source[3*cm->n_fc + k];
  }

  /* 2. Store divergence operator in non assembly
   * ============================================ */

  cs_real_t  *_div = msles->div_op + 3*connect->c2f->idx[cm->c_id];
  memcpy(_div, div_op, 3*cm->n_fc*sizeof(cs_real_t));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb schemes when the GKB or ALU algorithm is used as solver.
 *         Rely on cs_cdofb_vecteq_assembly()
 *
 * \param[in]       csys              pointer to a cs_cell_sys_t structure
 * \param[in]       cm                pointer to a cs_cell_mesh_t structure
 * \param[in]       div_op            array with the divergence op. values
 * \param[in]       has_sourceterm    has the equation a source term?
 * \param[in, out]  sc                pointer to scheme context structure
 * \param[in, out]  eqc               context structure for a vector-valued Fb
 * \param[in, out]  eqa               pointer to cs_equation_assemble_t
 */
/*----------------------------------------------------------------------------*/

static void
_velocity_full_assembly(const cs_cell_sys_t            *csys,
                        const cs_cell_mesh_t           *cm,
                        const cs_real_t                *div_op,
                        const bool                      has_sourceterm,
                        cs_cdofb_monolithic_t          *sc,
                        cs_cdofb_vecteq_t              *eqc,
                        cs_equation_assemble_t         *eqa)
{
  const short int  n_f = cm->n_fc;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_range_set_t  *rs = cs_shared_range_set;

  cs_matrix_assembler_values_t   *mav = sc->mav_structures[0];
  cs_real_t  *rhs = sc->msles->b_f;

  /* 1. Matrix assembly
   * ================== */

  if (sc->msles->graddiv_coef > 0.) {
    cs_real_t  gamma = sc->msles->graddiv_coef / cm->vol_c;
    cs_cdofb_navsto_add_grad_div(cm->n_fc, gamma, div_op, csys->mat);
  }

  cs_cdofb_vecteq_assembly(csys, rs, cm, has_sourceterm, eqc, eqa, mav, rhs);

  /* 2. Store divergence operator in non assembly
   * ============================================ */

  cs_real_t  *_div = sc->msles->div_op + 3*connect->c2f->idx[cm->c_id];
  memcpy(_div, div_op, 3*n_f*sizeof(cs_real_t));

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb schemes
 *         Shares similarities with cs_equation_assemble_block_matrix()
 *
 * \param[in]      csys             pointer to a cs_cell_sys_t structure
 * \param[in]      cm               pointer to a cs_cell_mesh_t structure
 * \param[in]      div_op           array with the divergence op. values
 * \param[in]      has_sourceterm   has the equation a source term ?
 * \param[in, out] sc               pointer to scheme context structure
 * \param[in, out] eqc              context structure for a vector-valued Fb
 * \param[in, out] eqa              pointer to cs_equation_assemble_t
 */
/*----------------------------------------------------------------------------*/

static void
_full_assembly(const cs_cell_sys_t            *csys,
               const cs_cell_mesh_t           *cm,
               const cs_real_t                *div_op,
               const bool                      has_sourceterm,
               cs_cdofb_monolithic_t          *sc,
               cs_cdofb_vecteq_t              *eqc,
               cs_equation_assemble_t         *eqa)
{
  CS_UNUSED(eqa);

  const short int  n_f = cm->n_fc;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_range_set_t  *rset = cs_shared_range_set;
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(bd->n_row_blocks == n_f);

  cs_real_t  *eqc_st = eqc->source_terms;
  cs_matrix_assembler_values_t   *mav = sc->mav_structures[0];
  cs_real_t  *rhs = sc->msles->b_f;

  /* 1. Matrix assembly
   * ================== */

  cs_gnum_t  r_gids[CS_CDO_ASSEMBLE_BUF_SIZE];
  cs_gnum_t  c_gids[CS_CDO_ASSEMBLE_BUF_SIZE];
  cs_real_t  values[CS_CDO_ASSEMBLE_BUF_SIZE];

  const cs_gnum_t  p_gid = rset->g_id[3*n_faces + cm->c_id];

  /* 1.a Add the contribution of velocity DoFs */
  int  bufsize = 0;
  for (int bi = 0; bi < bd->n_row_blocks; bi++) {

    /* dof_ids is an interlaced array and global numbering is not interlaced
       that's why one considers f_id */
    const cs_lnum_t  fi_id = cm->f_ids[bi];
    const cs_gnum_t  bi_gids[3]
      = {rset->g_id[fi_id],              /* x-component */
         rset->g_id[fi_id +   n_faces],  /* y-component */
         rset->g_id[fi_id + 2*n_faces]}; /* z-component */

    for (int bj = 0; bj < bd->n_col_blocks; bj++) {

      const cs_lnum_t  fj_id = cm->f_ids[bj];
      const cs_gnum_t  bj_gids[3]
        = {rset->g_id[fj_id],              /* x-component */
           rset->g_id[fj_id +   n_faces],  /* y-component */
           rset->g_id[fj_id + 2*n_faces]}; /* z-component */

      /* mIJ is a small square matrix of size 3 */
      cs_sdm_t  *mIJ = cs_sdm_get_block(m, bi, bj);

      for (short int ii = 0; ii < 3; ii++) {

        const cs_gnum_t  bi_gid = bi_gids[ii];

        for (short int jj = 0; jj < 3; jj++) {

          /* Add an entry */
          r_gids[bufsize] = bi_gid;
          c_gids[bufsize] = bj_gids[jj];
          values[bufsize] = mIJ->val[3*ii + jj];
          bufsize += 1;

          if (bufsize == CS_CDO_ASSEMBLE_BUF_SIZE) {
#           pragma omp critical
            cs_matrix_assembler_values_add_g(mav, bufsize,
                                             r_gids, c_gids, values);
            bufsize = 0;
          }

        } /* jj */
      } /* ii */

    } /* Loop on column blocks */

    /* 1.b Add the contribution of pressure DoFs */
    for (short int ii = 0; ii < 3; ii++) { /* x,y,z-component */

      r_gids[bufsize] = bi_gids[ii];
      c_gids[bufsize] = p_gid;
      values[bufsize] = div_op[3*bi+ii];
      bufsize += 1;

      if (bufsize == CS_CDO_ASSEMBLE_BUF_SIZE) {
#       pragma omp critical
        cs_matrix_assembler_values_add_g(mav, bufsize, r_gids, c_gids, values);
        bufsize = 0;
      }

      /* Its transposed B_x, B_y, B_z */
      r_gids[bufsize] = p_gid;
      c_gids[bufsize] = bi_gids[ii];
      values[bufsize] = div_op[3*bi+ii];
      bufsize += 1;

      if (bufsize == CS_CDO_ASSEMBLE_BUF_SIZE) {
#       pragma omp critical
        cs_matrix_assembler_values_add_g(mav, bufsize, r_gids, c_gids, values);
        bufsize = 0;
      }

    } /* Loop on components */

  } /* Loop on row blocks (bi) */

  if (bufsize > 0) {
#   pragma omp critical
    cs_matrix_assembler_values_add_g(mav, bufsize, r_gids, c_gids, values);
    bufsize = 0;
  }

  /* 2. RHS assembly (only the part with face DoFs)
   * ============================================== */

  for (short int f = 0; f < 3*n_f; f++)
#   pragma omp atomic
    rhs[csys->dof_ids[f]] += csys->rhs[f];

  /* Reset the value of the source term for the cell DoF
     Source term is only hold by the cell DoF in face-based schemes so
     there is no need to assemble this term. */
  if (has_sourceterm)
    for (int k = 0; k < 3; k++)
      eqc_st[3*cm->c_id + k] = csys->source[3*n_f + k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         steady-state case
 *
 * \param[in]         nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in]         dir_values  array storing the Dirichlet values
 * \param[in]         forced_ids  indirection in case of internal enforcement
 * \param[in, out]    sc          pointer to the scheme context
 */
/*----------------------------------------------------------------------------*/

static void
_steady_build(const cs_navsto_param_t      *nsp,
              const cs_real_t              *dir_values,
              const cs_lnum_t               forced_ids[],
              cs_cdofb_monolithic_t        *sc)
{
  assert(dir_values !=  NULL);  /* Sanity check */

  /* Retrieve shared structures */
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  t_cur = cs_shared_time_step->t_cur;
  const cs_real_t  t_eval = t_cur; /* dummy variable */

  /* Retrieve high-level structures */
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  /* Current values of the velocity field */
  const cs_real_t  *vel_c = sc->velocity->val;

  /* Initialize the matrix and all its related structures needed during
   * the assembly step */
  sc->init_system(sc);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, mom_eq, mom_eqp, mom_eqb, mom_eqc, nsp,        \
         dir_values, forced_ids, vel_c, sc)                             \
  firstprivate(t_eval)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(connect);

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, t_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, cs_equation_cell_mesh_flag(cell_flag, mom_eqb),
                         connect, quant, cm);

      /* For the stationary problem, the global system writes:
       *
       *     |        |         |
       *     |   A    |    Bt   |  B is the divergence (Bt the gradient)
       *     |        |         |  A is csys->mat in what follows
       *     |--------|---------|  The viscous part arising from the CDO-Fb
       *     |        |         |  schemes for vector-valued variables and
       *     |   B    |    0    |  additional terms as the linearized
       *     |        |         |  convective term
       */

      /* Set the local (i.e. cellwise) structures for the current cell */
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, mom_eqp, mom_eqb, mom_eqc,
                                       dir_values, forced_ids, vel_c, t_eval,
                                       csys, cb);

      /* 1- SETUP THE NAVSTO LOCAL BUILDER
       * =================================
       * - Set the type of boundary
       * - Set the pressure boundary conditions (if required)
       * - Define  the divergence operator used in the linear system
       */

      cs_cdofb_navsto_define_builder(t_eval, nsp, cm, csys,
                                     sc->pressure_bc, sc->bf_type,
                                     &nsb);

      /* 2- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_conv_diff_reac(t_eval, mom_eqp, mom_eqc, cm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Local system after advection/diffusion/reaction",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const bool has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm) {

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   t_eval, 1., /* time, scaling */
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* 3b- OTHER RHS CONTRIBUTIONS
       * ===========================
       *
       * First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _apply_bc_partly(sc, mom_eqp, cm, nsb.bf_type, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Local system matrix before condensation", csys);
#endif

      /* 5- STATIC CONDENSATION
       * ======================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */

      cs_static_condensation_vector_eq(connect->c2f,
                                       mom_eqc->rc_tilda,
                                       mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _apply_remaining_bc(mom_eqp, cm, sc, csys, cb, &nsb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      sc->assemble(csys, cm, nsb.div_op, has_sourceterm, sc, mom_eqc, eqa);

    } /* Main loop on cells */

    /* Free temporary buffer */
    cs_cdofb_navsto_free_builder(&nsb);

  } /* End of the OpenMP Block */

  for (int i = 0; i < sc->msles->n_row_blocks*sc->msles->n_row_blocks; i++) {
    cs_matrix_assembler_values_done(sc->mav_structures[i]); /* optional */
    cs_matrix_assembler_values_finalize(&(sc->mav_structures[i]));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         case of an implicit Euler time scheme
 *
 * \param[in]         nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in]         dir_values  array storing the Dirichlet values
 * \param[in]         forced_ids  indirection in case of internal enforcement
 * \param[in, out]    sc          pointer to the scheme context
 */
/*----------------------------------------------------------------------------*/

static void
_implicit_euler_build(const cs_navsto_param_t  *nsp,
                      const cs_real_t          *dir_values,
                      const cs_lnum_t           forced_ids[],
                      cs_cdofb_monolithic_t    *sc)
{
  assert(dir_values !=  NULL);  /* Sanity check */

  /* Retrieve high-level structures */
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  /* Retrieve shared structures */
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  t_eval = t_cur + dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  /* Current values of the velocity field */
  cs_real_t  *vel_c = sc->velocity->val;

  /* Initialize the matrix and all its related structures needed during
   * the assembly step */
  sc->init_system(sc);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, mom_eq, mom_eqp, mom_eqb, mom_eqc, nsp,        \
         dir_values, forced_ids, vel_c, sc)                             \
  firstprivate(t_eval, inv_dtcur)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(connect);

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, t_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, cs_equation_cell_mesh_flag(cell_flag, mom_eqb),
                         connect, quant, cm);

      /* The global system problem writes:
       *
       *     |        |         |
       *     |   A    |    Bt   |  B is the divergence (Bt the gradient)
       *     |        |         |  A is csys->mat in what follows
       *     |--------|---------|  The viscous part arising from the CDO-Fb
       *     |        |         |  schemes for vector-valued variables and
       *     |   B    |    0    |  additional terms as the linearized
       *     |        |         |  convective term + an unsteady term
       */

      /* Set the local (i.e. cellwise) structures for the current cell */
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, mom_eqp, mom_eqb, mom_eqc,
                                       dir_values, forced_ids, vel_c, t_eval,
                                       csys, cb);

      /* 1- SETUP THE NAVSTO LOCAL BUILDER *
       * ================================= *
       * - Set the type of boundary
       * - Set the pressure boundary conditions (if required)
       * - Define  the divergence operator used in the linear system
       */
      cs_cdofb_navsto_define_builder(t_eval, nsp, cm, csys,
                                     sc->pressure_bc, sc->bf_type,
                                     &nsb);

      /* 2- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_conv_diff_reac(t_eval, mom_eqp, mom_eqc, cm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Local system after advection/diffusion/reaction",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm) {

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   t_eval, 1., /* scaling if theta */
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* 3b- OTHER RHS CONTRIBUTIONS *
       * =========================== *
       *
       * First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _apply_bc_partly(sc, mom_eqp, cm, nsb.bf_type, csys, cb);

      /* 4- TIME CONTRIBUTION (mass lumping or vaornoï) */
      /* ==================== */

      if (mom_eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) {

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Get cell-cell block */
        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, cm->n_fc, cm->n_fc);

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*cm->n_fc + k] += ptyc * csys->val_n[3*cm->n_fc+k];
          /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        }

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Only diagonal time treatment available so far.\n");

      /* 5- STATIC CONDENSATION
       * ======================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */

      cs_static_condensation_vector_eq(connect->c2f,
                                       mom_eqc->rc_tilda,
                                       mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _apply_remaining_bc(mom_eqp, cm, sc, csys, cb, &nsb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      sc->assemble(csys, cm, nsb.div_op, has_sourceterm, sc, mom_eqc, eqa);

    } /* Main loop on cells */

    /* Free temporary buffer */
    cs_cdofb_navsto_free_builder(&nsb);

  } /* OPENMP Block */

  for (int i = 0; i < sc->msles->n_row_blocks*sc->msles->n_row_blocks; i++) {
    cs_matrix_assembler_values_done(sc->mav_structures[i]); /* optional */
    cs_matrix_assembler_values_finalize(&(sc->mav_structures[i]));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         case of a theta time scheme
 *
 * \param[in]         nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in]         dir_values  array storing the Dirichlet values
 * \param[in]         forced_ids  indirection in case of internal enforcement
 * \param[in, out]    sc          pointer to the scheme context
 */
/*----------------------------------------------------------------------------*/

static void
_theta_scheme_build(const cs_navsto_param_t  *nsp,
                    const cs_real_t          *dir_values,
                    const cs_lnum_t           forced_ids[],
                    cs_cdofb_monolithic_t    *sc)
{
  assert(dir_values !=  NULL);  /* Sanity check */

  /* Retrieve high-level structures */
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  /* Retrieve shared structures */
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const double  tcoef = 1 - mom_eqp->theta;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  /* Time_eval = (1-theta).t^n + theta.t^(n+1) = t^n + theta.dt
   * since t^(n+1) = t^n + dt
   */
  const cs_real_t  t_eval = t_cur + mom_eqp->theta*dt_cur;

  /* Current values of the velocity field */
  const cs_real_t  *vel_c = sc->velocity->val;

  /* Initialize the matrix and all its related structures needed during
   * the assembly step */
  sc->init_system(sc);

  /* Detect the first call (in this case, we compute the initial source term)*/
  bool  compute_initial_source = false;
  if (ts->nt_cur == ts->nt_prev || ts->nt_prev == 0)
    compute_initial_source = true;

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(quant, connect, mom_eq, mom_eqp, mom_eqb, mom_eqc, nsp,        \
         dir_values, vel_c, sc, compute_initial_source)                 \
  firstprivate(t_eval, t_cur, dt_cur, tcoef, inv_dtcur)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(connect);

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, t_eval, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id,
                         cs_equation_cell_mesh_flag(cell_flag, mom_eqb),
                         connect, quant, cm);

      /* Starts from the stationary Stokes problem where
       *
       *     |        |         |
       *     |   A    |    Bt   |  B is the divergence (Bt the gradient)
       *     |        |         |  A is csys->mat in what follows
       *     |--------|---------|  The viscous part arising from the CDO-Fb
       *     |        |         |  schemes for vector-valued variables
       *     |   B    |    0    |
       *     |        |         |
       */

      /* Set the local (i.e. cellwise) structures for the current cell */
      cs_cdofb_vecteq_init_cell_system(cell_flag, cm, mom_eqp, mom_eqb, mom_eqc,
                                       dir_values, forced_ids, vel_c, t_eval,
                                       csys, cb);

      /* 1- SETUP THE NAVSTO LOCAL BUILDER *
       * ================================= *
       * - Set the type of boundary
       * - Set the pressure boundary conditions (if required)
       * - Define  the divergence operator used in the linear system
       */
      cs_cdofb_navsto_define_builder(t_cur + dt_cur, nsp, cm, csys,
                                     sc->pressure_bc, sc->bf_type,
                                     &nsb);

      /* 2- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_conv_diff_reac(t_eval, mom_eqp, mom_eqc, cm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Local system after advection/diffusion/reaction",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm) {

        if (compute_initial_source) { /* First time step:
                                       * Compute source term at a previous step
                                       */

          cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                     t_cur, tcoef,  /* time, scaling */
                                     cb, mom_eqb, csys);

        }
        else { /* Add the contribution of the previous time step */

          for (short int k = 0; k < 3; k++)
            csys->rhs[3*cm->n_fc + k] += tcoef*mom_eqc->source_terms[3*c_id+k];

        }

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   /* time       , scaling */
                                   t_cur + dt_cur, mom_eqp->theta,
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* 3b- OTHER RHS CONTRIBUTIONS *
       * =========================== *
       *
       * First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _apply_bc_partly(sc, mom_eqp, cm, nsb.bf_type, csys, cb);

      /* 4- UNSTEADY TERM + TIME SCHEME
       * ============================== */

      /* STEP.1 >> Compute the contribution of the "adr" to the RHS:
       *           tcoef*adr_pn where adr_pn = csys->mat * p_n */
      double  *adr_pn = cb->values;
      cs_sdm_block_matvec(csys->mat, csys->val_n, adr_pn);
      for (short int i = 0; i < csys->n_dofs; i++) /* n_dofs = n_vc */
        csys->rhs[i] -= tcoef * adr_pn[i];

      /* STEP.2 >> Multiply csys->mat by theta */
      for (int i = 0; i < csys->n_dofs*csys->n_dofs; i++)
        csys->mat->val[i] *= mom_eqp->theta;

      /* STEP.3 >> Handle the mass matrix
       * Two contributions for the mass matrix
       *  a) add to csys->mat
       *  b) add to rhs mass_mat * p_n */
      if (mom_eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Get cell-cell block */
        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, cm->n_fc, cm->n_fc);

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*cm->n_fc + k] += ptyc * csys->val_n[3*cm->n_fc+k];
          /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        } /* Loop on k */

      }
      else
        bft_error(__FILE__, __LINE__, 0,
                  "Only diagonal time treatment available so far.");

      /* 5- STATIC CONDENSATION
       * ======================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */

      cs_static_condensation_vector_eq(connect->c2f,
                                       mom_eqc->rc_tilda,
                                       mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _apply_remaining_bc(mom_eqp, cm, sc, csys, cb, &nsb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      sc->assemble(csys, cm, nsb.div_op, has_sourceterm, sc, mom_eqc, eqa);

    } /* Main loop on cells */

    /* Free temporary buffer */
    cs_cdofb_navsto_free_builder(&nsb);

  } /* OPENMP Block */

  for (int i = 0; i < sc->msles->n_row_blocks*sc->msles->n_row_blocks; i++) {
    cs_matrix_assembler_values_done(sc->mav_structures[i]); /* optional */
    cs_matrix_assembler_values_finalize(&(sc->mav_structures[i]));
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the divergence of the velocity fields and compute its
 *         weighted L2-norm (squared)
 *
 * \param[in, out]  sc       scheme context
 * \param[in, out]  mom_eqc  context of the momentum equation
 *
 * \return the L2-norm of the divergence
 */
/*----------------------------------------------------------------------------*/

static cs_real_t
_update_divergence(cs_cdofb_monolithic_t         *sc,
                   cs_cdofb_vecteq_t             *mom_eqc)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  *vel_f = mom_eqc->face_values;

  /* Update the divergence of the velocity field */
  cs_real_t  *div = sc->divergence->val;
  cs_real_t  norm2 = 0.;

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN) reduction(+:norm2)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    div[c_id] = cs_cdofb_navsto_cell_divergence(c_id,
                                                quant, cs_shared_connect->c2f,
                                                vel_f);
    norm2 += quant->cell_vol[c_id] * div[c_id] * div[c_id];
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
  cs_dbg_darray_to_listing("VELOCITY_DIV", quant->n_cells, div, 9);
#endif

  /* Parallel treatment */
  if (cs_glob_n_ranks > 1)
    cs_parall_sum(1, CS_REAL_TYPE, &norm2);

  if (norm2 > 0)
    norm2 = sqrt(norm2);

  return norm2;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the velocity fields at cells and rescale the pressure field
 *
 * \param[in, out]  sc       scheme context
 * \param[in, out]  mom_eqc  context of the momentum equation
 */
/*----------------------------------------------------------------------------*/

static void
_update_fields(cs_cdofb_monolithic_t         *sc,
               cs_cdofb_vecteq_t             *mom_eqc)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  *vel_f = mom_eqc->face_values;

  /* Update pressure, velocity and divergence fields */
  cs_field_t  *pr_fld = sc->pressure;
  cs_field_t  *vel_fld = sc->velocity;

  /* Compute values at cells vel_c from values at faces vel_f
     vel_c = acc^-1*(RHS - Acf*vel_f) */
  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda,
                                        mom_eqc->acf_tilda,
                                        vel_f, vel_fld->val);

  /* Rescale pressure */
  cs_cdofb_navsto_set_zero_mean_pressure(quant, pr_fld->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
  cs_dbg_darray_to_listing("VELOCITY", 3*quant->n_faces, vel_f, 9);
  cs_dbg_darray_to_listing("PRESSURE", quant->n_cells, pr_fld->val, 9);
#endif
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  nsp         pointer to NavSto parameter settings
 * \param[in]  mesh        pointer to a cs_mesh_t structure
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_init_common(const cs_navsto_param_t       *nsp,
                                const cs_mesh_t               *mesh,
                                const cs_cdo_quantities_t     *quant,
                                const cs_cdo_connect_t        *connect,
                                const cs_time_step_t          *time_step)
{
  /* Assign static const pointers */
  cs_shared_mesh = mesh;
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  /* Need to build special range set and interfaces ? */
  switch (nsp->sles_param.strategy) {

  case CS_NAVSTO_SLES_BY_BLOCKS:
    {
      cs_shared_range_set = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];
      cs_shared_matrix_structure = cs_cdofb_scaleq_matrix_structure();

      int  block_sizes[3] =
        {connect->n_max_fbyc, connect->n_max_fbyc, connect->n_max_fbyc};

      BFT_MALLOC(cs_cdofb_monolithic_cw_mat, cs_glob_n_threads, cs_sdm_t *);
      for (int i = 0; i < cs_glob_n_threads; i++)
        cs_cdofb_monolithic_cw_mat[i] = NULL;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
# pragma omp parallel
      {
        int t_id = omp_get_thread_num();
        assert(t_id < cs_glob_n_threads);
        cs_cdofb_monolithic_cw_mat[t_id] = cs_sdm_block_create(3, 3,
                                                               block_sizes,
                                                               block_sizes);
      }
#else
      assert(cs_glob_n_threads == 1);
      cs_cdofb_monolithic_cw_mat[0] = cs_sdm_block_create(3, 3,
                                                          block_sizes,
                                                          block_sizes);

#endif /* openMP */
    }
    break;

  case CS_NAVSTO_SLES_GKB_SATURNE:
  case CS_NAVSTO_SLES_UZAWA_AL:
    cs_shared_range_set = connect->range_sets[CS_CDO_CONNECT_FACE_VP0];
    cs_shared_matrix_structure = cs_cdofb_vecteq_matrix_structure();
    break;

  default: /* Build the fully coupled system */
    _build_shared_structures();

    cs_shared_interface_set = _shared_interface_set;
    cs_shared_range_set = _shared_range_set;
    cs_shared_matrix_structure = _shared_matrix_structure;
    cs_shared_matrix_assembler = _shared_matrix_assembler;
    break;

  }

  /* SLES needs these structures for advanced PETSc hooks */
  cs_cdofb_monolithic_sles_set_shared(connect, quant, cs_shared_range_set);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free shared pointers with lifecycle dedicated to this file
 *
 * \param[in]  nsp         pointer to NavSto parameter settings
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_finalize_common(const cs_navsto_param_t       *nsp)
{
  /* Need to build special range set and interfaces ? */
  switch (nsp->sles_param.strategy) {

  case CS_NAVSTO_SLES_BY_BLOCKS:
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#   pragma omp parallel
    {
      int t_id = omp_get_thread_num();
      cs_cdofb_monolithic_cw_mat[t_id] =
        cs_sdm_free(cs_cdofb_monolithic_cw_mat[t_id]);
    }
#else
    assert(cs_glob_n_threads == 1);
    cs_cdofb_monolithic_cw_mat[0] =
      cs_sdm_free(cs_cdofb_monolithic_cw_mat[0]);
#endif /* openMP */

    BFT_FREE(cs_cdofb_monolithic_cw_mat);
    break;

  default:
    break; /* Nothing to do */
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a \ref cs_cdofb_monolithic_t structure
 *
 * \param[in] nsp        pointer to a \ref cs_navsto_param_t structure
 * \param[in] bf_type    type of boundary for each boundary face
 * \param[in] nsc_input  pointer to a \ref cs_navsto_ac_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_monolithic_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_monolithic_init_scheme_context(const cs_navsto_param_t   *nsp,
                                        cs_boundary_type_t        *bf_type,
                                        void                      *nsc_input)
{
  /* Sanity checks */
  assert(nsp != NULL && nsc_input != NULL);
  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n", __func__);

  /* Navier-Stokes scheme context (SC) */
  cs_cdofb_monolithic_t  *sc = NULL;

  BFT_MALLOC(sc, 1, cs_cdofb_monolithic_t);

  /* Cast the coupling context (CC) */
  cs_navsto_monolithic_t  *cc = (cs_navsto_monolithic_t  *)nsc_input;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;

  sc->coupling_context = cc; /* shared with cs_navsto_system_t */

  /* Quick access to the main fields */
  sc->velocity = cs_field_by_name("velocity");
  sc->pressure = cs_field_by_name("pressure");

  if (nsp->post_flag & CS_NAVSTO_POST_VELOCITY_DIVERGENCE)
    sc->divergence = cs_field_by_name("velocity_divergence");
  else
    sc->divergence = NULL;

  /* Boundary treatment */
  sc->bf_type = bf_type;

  /* Processing of the pressure boundary condition */
  sc->pressure_bc = cs_cdo_bc_face_define(CS_CDO_BC_HMG_NEUMANN, /* Default */
                                          true,      /* Steady BC up to now */
                                          1,         /* Dimension */
                                          nsp->n_pressure_bc_defs,
                                          nsp->pressure_bc_defs,
                                          cs_shared_quant->n_b_faces);

  /* Set the way to enforce the Dirichlet BC on the velocity
   * "fixed_wall" means a no-slip BC */
  mom_eqb->bd_msh_flag |= CS_FLAG_COMP_PFC;

  sc->apply_symmetry = cs_cdofb_symmetry;
  sc->apply_sliding_wall = cs_cdofb_block_dirichlet_alge;
  sc->apply_fixed_wall = cs_cdofb_block_dirichlet_alge;

  switch (mom_eqp->default_enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_alge;
    break;

  case CS_PARAM_BC_ENFORCE_PENALIZED:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_pena;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_weak;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_wsym;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of algorithm to enforce Dirichlet BC.",
              __func__);

  }

  /* Set the build function */
  sc->steady_build = _steady_build;

  switch (nsp->time_scheme) {

  case CS_TIME_SCHEME_STEADY:
    sc->build = _steady_build;
    break;

  case CS_TIME_SCHEME_EULER_IMPLICIT:
    sc->build = _implicit_euler_build;
    break;

  case CS_TIME_SCHEME_EULER_EXPLICIT:
  case CS_TIME_SCHEME_CRANKNICO:
  case CS_TIME_SCHEME_THETA:
    sc->build = _theta_scheme_build;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid time scheme.", __func__);

  } /* Switch on time schme */

  /* Handle the resolution of a saddle-point system */
  cs_cdofb_monolithic_sles_t  *msles = cs_cdofb_monolithic_sles_create();

  /* Set the solve and assemble functions */
  switch (nsp->sles_param.strategy) {

  case CS_NAVSTO_SLES_BY_BLOCKS:
    sc->init_system = _init_system_by_blocks;
    sc->solve = cs_cdofb_monolithic_by_blocks_solve;
    sc->assemble = _assembly_by_blocks;
    sc->elemental_assembly = cs_equation_assemble_set(CS_SPACE_SCHEME_CDOFB,
                                                      CS_CDO_CONNECT_FACE_SP0);

    BFT_MALLOC(sc->mav_structures, 9, cs_matrix_assembler_values_t *);

    msles->graddiv_coef = nsp->gd_scale_coef;
    msles->n_row_blocks = 3;
    BFT_MALLOC(msles->block_matrices, 9, cs_matrix_t *);
    BFT_MALLOC(msles->div_op,
               3*cs_shared_connect->c2f->idx[cs_shared_quant->n_cells],
               cs_real_t);
    break;

  case CS_NAVSTO_SLES_GKB_SATURNE:
    sc->init_system = _init_system_default;
    sc->solve = cs_cdofb_monolithic_gkb_solve;
    sc->assemble = _velocity_full_assembly;
    sc->elemental_assembly = cs_equation_assemble_set(CS_SPACE_SCHEME_CDOFB,
                                                      CS_CDO_CONNECT_FACE_VP0);

    BFT_MALLOC(sc->mav_structures, 1, cs_matrix_assembler_values_t *);

    msles->graddiv_coef = nsp->gd_scale_coef;
    msles->n_row_blocks = 1;
    BFT_MALLOC(msles->block_matrices, 1, cs_matrix_t *);
    BFT_MALLOC(msles->div_op,
               3*cs_shared_connect->c2f->idx[cs_shared_quant->n_cells],
               cs_real_t);
    break;

  case CS_NAVSTO_SLES_UZAWA_AL:
    sc->init_system = _init_system_default;
    sc->solve = cs_cdofb_monolithic_uzawa_al_incr_solve;
    sc->assemble = _velocity_full_assembly;
    sc->elemental_assembly = cs_equation_assemble_set(CS_SPACE_SCHEME_CDOFB,
                                                      CS_CDO_CONNECT_FACE_VP0);

    BFT_MALLOC(sc->mav_structures, 1, cs_matrix_assembler_values_t *);

    msles->graddiv_coef = nsp->gd_scale_coef;
    msles->n_row_blocks = 1;
    BFT_MALLOC(msles->block_matrices, 1, cs_matrix_t *);
    BFT_MALLOC(msles->div_op,
               3*cs_shared_connect->c2f->idx[cs_shared_quant->n_cells],
               cs_real_t);
    break;

  default:
    sc->init_system = _init_system_default;
    sc->solve = cs_cdofb_monolithic_solve;
    sc->assemble = _full_assembly;
    sc->elemental_assembly = NULL;

    BFT_MALLOC(sc->mav_structures, 1, cs_matrix_assembler_values_t *);

    msles->n_row_blocks = 1;
    BFT_MALLOC(msles->block_matrices, 1, cs_matrix_t *);
    break;

  }

  /* Set the pointer storing linear algebra features */
  sc->msles = msles;

  /* Monitoring */
  CS_TIMER_COUNTER_INIT(sc->timer);

  return sc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdofb_monolithic_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_monolithic_free_scheme_context(void   *scheme_context)
{
  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t  *)scheme_context;

  if (sc == NULL)
    return sc;

  /* Free BC structure */
  sc->pressure_bc = cs_cdo_bc_free(sc->pressure_bc);

  /* Shared structures which have been allocated only for the function used
   * in the monolithic approach have to be freed */
  if (_shared_interface_set != NULL)
    cs_interface_set_destroy(&_shared_interface_set);
  if (_shared_range_set != NULL)
    cs_range_set_destroy(&_shared_range_set);
  if (_shared_matrix_assembler != NULL)
    cs_matrix_assembler_destroy(&_shared_matrix_assembler);
  if (_shared_matrix_structure != NULL)
    cs_matrix_structure_destroy(&_shared_matrix_structure);

  /* Unset shared pointers */
  cs_shared_range_set = NULL;
  cs_shared_matrix_structure = NULL;
  cs_shared_matrix_assembler = NULL;
  cs_shared_interface_set = NULL;

  BFT_FREE(sc->mav_structures);

  /* Free the context structure for solving saddle-point system */
  cs_cdofb_monolithic_sles_free(&(sc->msles));

  /* Other pointers are only shared (i.e. not owner) */
  BFT_FREE(sc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Stokes or Oseen system with a CDO face-based scheme
 *         using a monolithic approach and GKB algorithm
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_steady(const cs_mesh_t            *mesh,
                           const cs_navsto_param_t    *nsp,
                           void                       *scheme_context)
{
  cs_timer_t  t_start = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t *)scheme_context;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces, n_cells = quant->n_cells;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;

  /* Build an array storing the Dirichlet values at faces and ids of DoFs if
   * an enforcement of (internal) DoFs is requested */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *enforced_ids = NULL;

  cs_cdofb_vecteq_setup(t_cur, mesh, mom_eqp, mom_eqb,
                        &dir_values, &enforced_ids);

  /* Initialize the rhs */
  cs_cdofb_monolithic_sles_init(n_cells, n_faces, sc->msles);

  /* Main loop on cells to define the linear system to solve */
  sc->steady_build(nsp, dir_values, enforced_ids, sc);

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(enforced_ids);

  /* End of the system building */
  cs_timer_t  t_bld_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_bld_end);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */
  _fields_to_previous(sc);

  /* Solve the linear system */
  cs_timer_t  t_solve_start = cs_timer_time();
  cs_cdofb_monolithic_sles_t  *msles = sc->msles;

  /* matrix has been already assigned to the msles structure */
  msles->sles = cs_sles_find_or_add(mom_eq->field_id, NULL);
  msles->u_f = mom_eqc->face_values; /* velocity DoFs at faces */
  msles->p_c = sc->pressure->val;    /* pressure DoFs at cells */

  int  cumulated_inner_iters = sc->solve(nsp, mom_eqp, msles);

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

  /* Compute the velocity divergence and retrieve its L2-norm */
  cs_real_t  div_l2_norm = _update_divergence(sc, mom_eqc);

  if (nsp->verbosity > 1) {
    cs_log_printf(CS_LOG_DEFAULT,
                  " -cvg- cumulated_inner_iters: %d ||div(u)|| = %6.4e\n",
                  cumulated_inner_iters, div_l2_norm);
    cs_log_printf_flush(CS_LOG_DEFAULT);
  }

  /* Now compute/update the velocity and pressure fields */
  _update_fields(sc, mom_eqc);

  /* Frees */
  cs_cdofb_monolithic_sles_clean(msles);

  cs_timer_t  t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Navier-Stokes system with a CDO face-based scheme
 *         using a monolithic approach and Picard iterations to solve the
 *         on-linearities arising from the advection term
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_steady_nl(const cs_mesh_t           *mesh,
                              const cs_navsto_param_t   *nsp,
                              void                      *scheme_context)
{
  cs_timer_t  t_start = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t *)scheme_context;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;

  /*--------------------------------------------------------------------------
   *                    INITIAL BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces, n_cells = quant->n_cells;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;

  /* Build an array storing the Dirichlet values at faces and ids of DoFs if
   * an enforcement of (internal) DoFs is requested */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *enforced_ids = NULL;

  cs_cdofb_vecteq_setup(t_cur, mesh, mom_eqp, mom_eqb,
                        &dir_values, &enforced_ids);

  /* Initialize the rhs */
  cs_cdofb_monolithic_sles_init(n_cells, n_faces, sc->msles);

  /* Main loop on cells to define the linear system to solve */
  sc->steady_build(nsp, dir_values, enforced_ids, sc);

  /* End of the system building */
  cs_timer_t  t_build_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_build_end);

  /*--------------------------------------------------------------------------
   *                   INITIAL BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */
  _fields_to_previous(sc);

  /* Solve the linear system */
  cs_timer_t  t_solve_start = cs_timer_time();
  cs_navsto_algo_info_t  ns_info;
  cs_navsto_algo_info_init(&ns_info);

  cs_cdofb_monolithic_sles_t  *msles = sc->msles;

  msles->sles = cs_sles_find_or_add(mom_eq->field_id, NULL);
  msles->u_f = mom_eqc->face_values; /* velocity DoFs at faces */
  msles->p_c = sc->pressure->val;    /* pressure DoFs at cells */

  ns_info.n_inner_iter =
    (ns_info.last_inner_iter = sc->solve(nsp, mom_eqp, msles));
  ns_info.n_algo_iter += 1;

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

  /* Compute the velocity divergence and retrieve its L2-norm */
  cs_real_t  div_l2_norm = _update_divergence(sc, mom_eqc);

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: START
   *--------------------------------------------------------------------------*/

  if (nsp->verbosity > 1) {
    cs_navsto_algo_info_header("Picard");
    cs_navsto_algo_info_printf("Picard", ns_info, div_l2_norm);
  }

  while (ns_info.cvg == CS_SLES_ITERATING) {

    /* Main loop on cells to define the linear system to solve */
    cs_timer_t  t_build_start = cs_timer_time();

    /* rhs set to zero and cs_sles_t structure is freed in order to do
     * the setup once again since the matrix should be modified */
    cs_cdofb_monolithic_sles_reset(msles);

    sc->steady_build(nsp, dir_values, enforced_ids, sc);

    /* End of the system building */
    t_build_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_build_start, &t_build_end);

    /* Current to previous */
    memcpy(mom_eqc->face_values_pre, mom_eqc->face_values,
           3*n_faces*sizeof(cs_real_t));

    /* Solve the new system */
    t_solve_start = cs_timer_time();

    ns_info.n_inner_iter +=
      (ns_info.last_inner_iter = sc->solve(nsp, mom_eqp, msles));

    t_solve_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

    /* Compute the velocity divergence and retrieve its L2-norm */
    div_l2_norm = _update_divergence(sc, mom_eqc);

    /* Check convergence status and update ns_info following the convergence */
    _picard_cvg_test(nsp, mom_eqc, div_l2_norm, &ns_info);

    if (nsp->verbosity > 1)
      cs_navsto_algo_info_printf("Picard", ns_info, div_l2_norm);

  } /* Loop on Picard iterations */

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: END
   *--------------------------------------------------------------------------*/

  /* Now compute/update the velocity and pressure fields */
  _update_fields(sc, mom_eqc);

  /* Frees */
  cs_cdofb_monolithic_sles_clean(msles);
  BFT_FREE(dir_values);
  BFT_FREE(enforced_ids);

  cs_timer_t  t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a monolithic approach.
 *         According to the settings, this function can handle either an
 *         implicit Euler time scheme or more generally a theta time scheme.
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic(const cs_mesh_t          *mesh,
                    const cs_navsto_param_t  *nsp,
                    void                     *scheme_context)
{
  const cs_timer_t  t_start = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t *)scheme_context;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces, n_cells = quant->n_cells;
  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_eval = ts->t_cur + ts->dt[0];

  /* Build an array storing the Dirichlet values at faces and ids of DoFs if
   * an enforcement of (internal) DoFs is requested */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *enforced_ids = NULL;

  cs_cdofb_vecteq_setup(t_eval, mesh, mom_eqp, mom_eqb,
                        &dir_values, &enforced_ids);

  /* Initialize the rhs */
  cs_cdofb_monolithic_sles_init(n_cells, n_faces, sc->msles);

  /* Main loop on cells to define the linear system to solve */
  sc->build(nsp, dir_values, enforced_ids, sc);

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  BFT_FREE(enforced_ids);

  /* End of the system building */
  cs_timer_t  t_bld_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_bld_end);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */
  _fields_to_previous(sc);

  /* Solve the linear system */
  cs_timer_t  t_solve_start = cs_timer_time();
  cs_cdofb_monolithic_sles_t  *msles = sc->msles;

  msles->sles = cs_sles_find_or_add(mom_eq->field_id, NULL);
  msles->u_f = mom_eqc->face_values; /* velocity DoFs at faces */
  msles->p_c = sc->pressure->val;    /* pressure DoFs at cells */

  sc->solve(nsp, mom_eqp, msles);

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

  /* Compute the velocity divergence and retrieve its L2-norm */
  cs_real_t  div_l2_norm = _update_divergence(sc, mom_eqc);

  if (nsp->verbosity > 1)
    cs_log_printf(CS_LOG_DEFAULT, " ||div(u)|| = %6.4e\n", div_l2_norm);

  /* Now compute/update the velocity and pressure fields */
  _update_fields(sc, mom_eqc);

  /* Frees */
  cs_cdofb_monolithic_sles_clean(msles);

  cs_timer_t  t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a monolithic approach.
 *         According to the settings, this function can handle either an
 *         implicit Euler time scheme or more generally a theta time scheme.
 *         Rely on Picard iterations to solve the on-linearities arising from
 *         the advection term
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_nl(const cs_mesh_t           *mesh,
                       const cs_navsto_param_t   *nsp,
                       void                      *scheme_context)
{
  cs_timer_t  t_start = cs_timer_time();

  /* Retrieve high-level structures */
  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t *)scheme_context;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;

  /*--------------------------------------------------------------------------
   *                    INITIAL BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces, n_cells = quant->n_cells;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_eval = ts->t_cur + ts->dt[0];

  /* Build an array storing the Dirichlet values at faces and ids of DoFs if
   * an enforcement of (internal) DoFs is requested */
  cs_real_t  *dir_values = NULL;
  cs_lnum_t  *enforced_ids = NULL;

  cs_cdofb_vecteq_setup(t_eval, mesh, mom_eqp, mom_eqb,
                        &dir_values, &enforced_ids);

  /* Initialize the rhs */
  cs_cdofb_monolithic_sles_init(n_cells, n_faces, sc->msles);

  /* Main loop on cells to define the linear system to solve */
  sc->build(nsp, dir_values, enforced_ids, sc);

  /* End of the system building */
  cs_timer_t  t_build_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_build_end);

  /*--------------------------------------------------------------------------
   *                   INITIAL BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */
  _fields_to_previous(sc);

  /* Solve the linear system */
  cs_timer_t  t_solve_start = cs_timer_time();
  cs_navsto_algo_info_t  ns_info;
  cs_navsto_algo_info_init(&ns_info);

  cs_cdofb_monolithic_sles_t  *msles = sc->msles;

  msles->sles = cs_sles_find_or_add(mom_eq->field_id, NULL);
  msles->u_f = mom_eqc->face_values; /* velocity DoFs at faces */
  msles->p_c = sc->pressure->val;    /* pressure DoFs at cells */

  ns_info.n_inner_iter =
    (ns_info.last_inner_iter = sc->solve(nsp, mom_eqp, msles));
  ns_info.n_algo_iter += 1;

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

  /* Compute the velocity divergence and retrieve its L2-norm */
  cs_real_t  div_l2_norm = _update_divergence(sc, mom_eqc);

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: START
   *--------------------------------------------------------------------------*/

  if (nsp->verbosity > 1) {
    cs_navsto_algo_info_header("Picard");
    cs_navsto_algo_info_printf("Picard", ns_info, div_l2_norm);
  }

  while (ns_info.cvg == CS_SLES_ITERATING) {

    /* Start of the system building */
    cs_timer_t  t_build_start = cs_timer_time();

    /* rhs set to zero and cs_sles_t structure is freed in order to do
     * the setup once again since the matrix should be modified */
    cs_cdofb_monolithic_sles_reset(msles);

    /* Main loop on cells to define the linear system to solve */
    sc->build(nsp, dir_values, enforced_ids, sc);

    /* End of the system building */
    t_build_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_build_start, &t_build_end);

    /* Current to previous */
    memcpy(mom_eqc->face_values_pre, mom_eqc->face_values,
           3*n_faces*sizeof(cs_real_t));

    /* Solve the new system */
    t_solve_start = cs_timer_time();

    ns_info.n_inner_iter +=
      (ns_info.last_inner_iter = sc->solve(nsp, mom_eqp, msles));

    t_solve_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

    /* Compute the velocity divergence and retrieve its L2-norm */
    div_l2_norm = _update_divergence(sc, mom_eqc);

    /* Check convergence status and update ns_info following the convergence */
    _picard_cvg_test(nsp, mom_eqc, div_l2_norm, &ns_info);

    if (nsp->verbosity > 1)
      cs_navsto_algo_info_printf("Picard", ns_info, div_l2_norm);

  } /* Loop on Picard iterations */

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: END
   *--------------------------------------------------------------------------*/

  /* Now compute/update the velocity and pressure fields */
  _update_fields(sc, mom_eqc);

  /* Frees */
  cs_cdofb_monolithic_sles_clean(msles);
  BFT_FREE(dir_values);

  cs_timer_t  t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
