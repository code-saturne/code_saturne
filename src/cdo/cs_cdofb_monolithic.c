/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it as one block (monolithic approach of the velocity-pressure
 * coupling)
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
#include "cs_math.h"
#include "cs_navsto_coupling.h"
#include "cs_param.h"
#include "cs_post.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*! \struct cs_cdofb_monolithic_t
 *  \brief Context related to CDO face-based discretization when dealing with
 *         vector-valued unknowns
 */

typedef struct {

  /*! \var coupling_context
   *  Pointer to a \ref cs_navsto_monolithic_t (owned by \ref
   *  cs_navsto_system_t) containing the settings related to the monolithic
   *  approach
   */

  cs_navsto_monolithic_t   *coupling_context;

  /*!
   * @name Main field variables
   * Fields for every main variable of the equation. Got from cs_navsto_system_t
   */

  /*! \var velocity
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the velocity
   */

  cs_field_t               *velocity;

  /*! \var pressure
   *  Pointer to \ref cs_field_t (owned by \ref cs_navsto_system_t) containing
   *  the cell DoFs of the pressure
   */

  cs_field_t               *pressure;

  /*! \var divergence
   *  Pointer to \ref cs_real_t containing the values of the divergence on the
   *  cells
   */

  cs_field_t               *divergence;

  /*!
   * @}
   * @name Boundary conditions (BC) management
   * Routines and elements used for enforcing the BCs
   * @{
   *
   *  \var bf_type
   *  Array of boundary type for each boundary face. (Shared)
   */

  const cs_boundary_type_t       *bf_type;

  /*!
   * \var pressure_bc
   * Structure storing the metadata after processing the user-defined boundary
   * conditions related to the pressure field
   */

  cs_cdo_bc_face_t               *pressure_bc;

  /*! \var apply_fixed_wall
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  wall boundary (no slip boundary)
   *
   *  \var apply_sliding_wall
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  wall boundary (a tangential velocity is specified at the wall)
   *
   *  \var apply_velocity_inlet
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  boundary with a fixed velocity at the inlet
   *
   *  \var apply_symmetry
   *  \ref cs_cdo_apply_boundary_t function pointer defining how to apply a
   *  symmetry boundary
   */

  cs_cdo_apply_boundary_t        *apply_fixed_wall;
  cs_cdo_apply_boundary_t        *apply_sliding_wall;
  cs_cdo_apply_boundary_t        *apply_velocity_inlet;
  cs_cdo_apply_boundary_t        *apply_symmetry;

  /*!
   * @}
   * @name Performance monitoring
   * Monitoring the efficiency of the algorithm used to solve the Navier-Stokes
   * system
   * @{
   */

  /*! \var timer
   *  Cumulated elapsed time for building and solving the Navier--Stokes system
   */
  cs_timer_counter_t  timer;

  /*! @} */

} cs_cdofb_monolithic_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */
static const cs_mesh_t              *cs_shared_mesh;
static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;

static cs_range_set_t         *cs_shared_range_set = NULL;
static cs_interface_set_t     *cs_shared_interface_set = NULL;
static cs_matrix_structure_t  *cs_shared_matrix_structure = NULL;
static cs_matrix_assembler_t  *cs_shared_matrix_assembler = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset the two RHS arrays to zero
 *
 * \param[out]  face_rhs   First rhs array related to momentum eq. (3*n_faces)
 * \param[out]  cell_rhs   Second rhs array related to mass eq. (n_cells)
 */
/*----------------------------------------------------------------------------*/

static inline void
_set_rhs_to_zero(cs_real_t   *face_rhs,
                 cs_real_t   *cell_rhs)
{
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_lnum_t  n_cells = cs_shared_quant->n_cells;

#if defined(HAVE_OPENMP)
#pragma omp parallel if (n_faces > CS_THR_MIN)
  {
#   pragma omp for nowait
    for (cs_lnum_t i = 0; i < 3*n_faces; i++) face_rhs[i] = 0.;
#   pragma omp for nowait
    for (cs_lnum_t i = 0; i < n_cells; i++) cell_rhs[i] = 0.;
  }
#else
  memset(face_rhs, 0, 3*n_faces*sizeof(cs_real_t));
  memset(cell_rhs, 0, n_cells*sizeof(cs_real_t));
#endif
}

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

  cs_shared_interface_set = cs_interface_set_create(size,
                                                    NULL,
                                                    gnum,
                                                    m->periodicity,
                                                    0, NULL, NULL, NULL);

  cs_shared_range_set = cs_range_set_create(cs_shared_interface_set,
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
  cs_shared_matrix_assembler
    = cs_matrix_assembler_create(cs_shared_range_set->l_range, true);

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
      = {cs_shared_range_set->g_id[frow_id],              /* x-component */
         cs_shared_range_set->g_id[frow_id +   n_faces],  /* y-component */
         cs_shared_range_set->g_id[frow_id + 2*n_faces]}; /* z-component */

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
        = {cs_shared_range_set->g_id[fcol_id],              /* x-component */
           cs_shared_range_set->g_id[fcol_id + n_faces],    /* y-component */
           cs_shared_range_set->g_id[fcol_id + 2*n_faces]}; /* z-component */

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
      const cs_gnum_t  gcol_id = cs_shared_range_set->g_id[3*n_faces + ccol_id];

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

    cs_matrix_assembler_add_g_ids(cs_shared_matrix_assembler,
                                  n_entries, grows, gcols);
    assert(shift == n_entries);

  } /* Loop on face entities */

  /* 4. Build the matrix structure */
  cs_matrix_assembler_compute(cs_shared_matrix_assembler);

  cs_shared_matrix_structure
    = cs_matrix_structure_create_from_assembler(CS_MATRIX_MSR,
                                                cs_shared_matrix_assembler);

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
 * \param[in]      sc          pointer to a cs_cdofb_monolithic_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 * \param[in, out] nsb         builder structure for the NavSto system
 * \param[in, out] mass_rhs    pointer to the rhs for the mass eq. in this cell
 */
/*----------------------------------------------------------------------------*/

static void
_apply_remaining_bc(const cs_cdofb_monolithic_t   *sc,
                    const cs_equation_param_t     *eqp,
                    const cs_cell_mesh_t          *cm,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb,
                    cs_cdofb_navsto_builder_t     *nsb,
                    cs_real_t                     *mass_rhs)
{
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    const cs_boundary_type_t  *bf_type = nsb->bf_type;
    cs_real_t  *div_op = nsb->div_op;

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
        mass_rhs[0] -= cs_math_3_dot_product(csys->dir_values + 3*f,
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb schemes
 *         Shares similarities with cs_equation_assemble_block_matrix()
 *
 * \param[in]        csys              pointer to a cs_cell_sys_t structure
 * \param[in]        cm                pointer to a cs_cell_mesh_t structure
 * \param[in]        div_op            array with the divergence op. values
 * \param[in]        has_sourceterm    has the equation a source term?
 * \param[in, out]   mav               pointer to cs_matrix_assembler_values_t
 * \param[in, out]   rhs               right-end side of the system
 * \param[in, out]   eqc_st            source term from the context view
 */
/*----------------------------------------------------------------------------*/

static void
_assemble(const cs_cell_sys_t            *csys,
          const cs_cell_mesh_t           *cm,
          const cs_real_t                *div_op,
          const bool                      has_sourceterm,
          cs_matrix_assembler_values_t   *mav,
          cs_real_t                       rhs[],
          cs_real_t                       eqc_st[])
{
  const short int n_f = cm->n_fc;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_range_set_t  *rset = cs_shared_range_set;
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;

  /* Sanity checks */
  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(bd->n_row_blocks == n_f);

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

        const cs_gnum_t  i_gid = bi_gids[ii];

        for (short int jj = 0; jj < 3; jj++) {

          /* Add an entry */
          r_gids[bufsize] = i_gid;
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
     there is no need to assemble this term. Use directly mass_rhs */
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
 * \param[in, out]    sc          pointerr to the scheme context
 * \param[in, out]    matrix      pointer to a \ref cs_matrix_t structure
 * \param[in, out]    mom_rhs     rhs array related to the momentum eq.
 * \param[in, out]    mass_rhs    rhs array related to the mass eq.
 */
/*----------------------------------------------------------------------------*/

static void
_steady_build(const cs_navsto_param_t  *nsp,
              cs_cdofb_monolithic_t    *sc,
              const cs_real_t          *dir_values,
              cs_matrix_t              *matrix,
              cs_real_t                *mom_rhs,
              cs_real_t                *mass_rhs)
{
  /* Sanity check */
  assert(matrix != NULL);

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

  /* Reset RHS arrays */
  _set_rhs_to_zero(mom_rhs, mass_rhs);

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)      \
  shared(quant, connect, mom_eq, mom_eqp, mom_eqb, mom_eqc, matrix, nsp,  \
         mav, mom_rhs, mass_rhs, dir_values, vel_c, sc)                   \
  firstprivate(t_eval)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
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
                                       dir_values, vel_c, t_eval, csys, cb);

      /* 1- SETUP THE NAVSTO LOCAL BUILDER
       * =================================
       * - Set the type of boundary
       * - Set the pressure boundary conditions (if required)
       * - Define  the divergence operator used in the linear system */
      cs_cdofb_navsto_define_builder(t_eval, nsp, cm, csys,
                                     sc->pressure_bc, sc->bf_type,
                                     &nsb);

      /* 2- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */
      cs_cdofb_vecteq_advection_diffusion(t_eval, mom_eqp, mom_eqc, cm,
                                          csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion", csys);
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
        cs_cell_sys_dump(">> Local system matrix before condensation", csys);
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
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _apply_remaining_bc(sc, mom_eqp, cm, csys, cb, &nsb, mass_rhs + c_id);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      _assemble(csys, cm, nsb.div_op,
                has_sourceterm, mav, mom_rhs, mom_eqc->source_terms);

    } /* Main loop on cells */

    /* Free temporary buffer */
    cs_cdofb_navsto_free_builder(&nsb);

  } /* End of the OpenMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */
  cs_matrix_assembler_values_finalize(&mav);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         steady-state case
 *
 * \param[in]         nsp         pointer to a \ref cs_navsto_param_t structure
 * \param[in, out]    sc          pointerr to the scheme context
 * \param[in, out]    matrix      pointer to a \ref cs_matrix_t structure
 * \param[in, out]    mom_rhs     rhs array related to the momentum eq.
 * \param[in, out]    mass_rhs    rhs array related to the mass eq.
 */
/*----------------------------------------------------------------------------*/

static void
_implicit_euler_build(const cs_navsto_param_t  *nsp,
                      cs_cdofb_monolithic_t    *sc,
                      const cs_real_t          *dir_values,
                      cs_matrix_t              *matrix,
                      cs_real_t                *mom_rhs,
                      cs_real_t                *mass_rhs)
{
  /* Sanity check */
  assert(matrix != NULL);

  /* Retrieve shared structures */
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  t_eval = t_cur + dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  /* Retrieve high-level structures */
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  /* Current values of the velocity field */
  const cs_real_t  *vel_c = sc->velocity->val;

  /* Reset RHS arrays */
  _set_rhs_to_zero(mom_rhs, mass_rhs);

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)      \
  shared(quant, connect, mom_eq, mom_eqp, mom_eqb, mom_eqc, matrix, nsp,  \
         mav, mom_rhs, mass_rhs, dir_values, vel_c, sc)                   \
  firstprivate(t_eval, inv_dtcur)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
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
                                       dir_values, vel_c, t_eval,
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

      cs_cdofb_vecteq_advection_diffusion(t_eval, mom_eqp, mom_eqc, cm,
                                          csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after convection/diffusion", csys);
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

      /* 4- TIME CONTRIBUTION */
      /* ==================== */

      if (mom_eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping
                                                          or Hodge-Voronoi */

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
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _apply_remaining_bc(sc, mom_eqp, cm, csys, cb, &nsb, mass_rhs + c_id);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      _assemble(csys, cm, nsb.div_op,
                has_sourceterm, mav, mom_rhs, mom_eqc->source_terms);

    } /* Main loop on cells */

    /* Free temporary buffer */
    cs_cdofb_navsto_free_builder(&nsb);

  } /* OPENMP Block */

  cs_matrix_assembler_values_finalize(&mav);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a linear system arising from a scalar-valued CDO-Fb scheme
 *
 * \param[in, out] sles     pointer to a cs_sles_t structure
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in]      eqp      pointer to a cs_equation_param_t structure
 * \param[in, out] vel_f    initial velocity on faces
 * \param[in, out] pre_c    initial pressure in cells
 * \param[in, out] b_f      right-hand side (scatter/gather if needed) on faces
 * \param[in, out] b_c      right_hand side on cells (mass equation)
 *
 * \return the number of iterations of the solver
 */
/*----------------------------------------------------------------------------*/

static int
_solve_system(cs_sles_t                     *sles,
              const cs_matrix_t             *matrix,
              const cs_equation_param_t     *eqp,
              cs_real_t                     *vel_f,
              cs_real_t                     *pre_c,
              cs_real_t                     *b_f,
              cs_real_t                     *b_c)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_lnum_t  n_cols = cs_matrix_get_n_columns(matrix);
  const cs_lnum_t  n_scatter_elts = 3*n_faces + n_cells;

  cs_real_t  *xsol = NULL;
  cs_real_t  *b = NULL;

  BFT_MALLOC(xsol, n_cols, cs_real_t);
  BFT_MALLOC(b, n_scatter_elts, cs_real_t);

  /* De-interlace the velocity array and the rhs for the face DoFs */
# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(vel_f, b_f, xsol, b)                           \
  firstprivate(n_faces)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    xsol[f            ] = vel_f[3*f];
    xsol[f +   n_faces] = vel_f[3*f+1];
    xsol[f + 2*n_faces] = vel_f[3*f+2];

    b[f            ] = b_f[3*f];
    b[f +   n_faces] = b_f[3*f+1];
    b[f + 2*n_faces] = b_f[3*f+2];

  }

  /* Add the pressure related elements */
  memcpy(xsol + 3*n_faces, pre_c, n_cells*sizeof(cs_real_t));
  memcpy(b + 3*n_faces, b_c, n_cells*sizeof(cs_real_t));

  cs_range_set_t  *rset = cs_shared_range_set;
  int  n_iters = 0;
  double  residual = DBL_MAX;

  /* Prepare solving (handle parallelism) */
  cs_gnum_t  nnz = cs_equation_prepare_system(1,        /* stride */
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
                         xsol, xsol);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         b, b);

  cs_dbg_fprintf_system(eqp->name, cs_shared_time_step->nt_cur,
                        CS_CDOFB_MONOLITHIC_DBG,
                        xsol, b, 3*n_faces);
#endif

  /* Interlace xsol --> vel_f and pre_c */
# pragma omp parallel for if (CS_THR_MIN > n_faces)     \
  shared(vel_f, xsol)                                   \
  firstprivate(n_faces)
  for (cs_lnum_t f = 0; f < n_faces; f++) {

    vel_f[3*f]   = xsol[f];
    vel_f[3*f+1] = xsol[f +   n_faces];
    vel_f[3*f+2] = xsol[f + 2*n_faces];

  }

  /* Add the pressure related elements */
  memcpy(pre_c, xsol + 3*n_faces, n_cells*sizeof(cs_real_t));

  /* Free what can be freed at this stage */
  BFT_FREE(xsol);
  BFT_FREE(b);

  return n_iters;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system and update the velocity fields and the
 *         pressure field
 *
 * \param[in, out]  sc       scheme context
 * \param[in, out]  mom_eq   equation structure related to the momentum
 */
/*----------------------------------------------------------------------------*/

static void
_update_fields(cs_cdofb_monolithic_t         *sc,
               cs_equation_t                 *mom_eq)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Update pressure, velocity and divergence fields */

  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_field_t  *pr_fld = sc->pressure;
  cs_field_t  *vel_fld = sc->velocity;
  cs_real_t  *vel_f = mom_eqc->face_values;

  /* Compute values at cells vel_c from values at faces vel_f
     vel_c = acc^-1*(RHS - Acf*vel_f) */
  cs_static_condensation_recover_vector(connect->c2f,
                                        mom_eqc->rc_tilda,
                                        mom_eqc->acf_tilda,
                                        vel_f, vel_fld->val);

  /* Update the divergence of the velocity field */
  cs_real_t  *div = sc->divergence->val;

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++)
    div[c_id] = cs_cdofb_navsto_cell_divergence(c_id,
                                                quant, connect->c2f, vel_f);

  /* Rescale pressure */
  cs_cdofb_navsto_set_zero_mean_pressure(quant, pr_fld->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
  cs_dbg_darray_to_listing("VELOCITY", 3*quant->n_faces, vel_f, 9);
  cs_dbg_darray_to_listing("PRESSURE", quant->n_cells, pr_fld->val, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV", quant->n_cells, div, 9);
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
 * \param[in]  mesh        pointer to a cs_mesh_t structure
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_init_common(const cs_mesh_t               *mesh,
                                const cs_cdo_quantities_t     *quant,
                                const cs_cdo_connect_t        *connect,
                                const cs_time_step_t          *time_step)
{
  /* Assign static const pointers */
  cs_shared_mesh = mesh;
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;

  _build_shared_structures();

  /* SLES needs these structures for advanced PETSc hooks */
  cs_cdofb_monolithic_sles_set_shared(quant, cs_shared_range_set);
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
  sc->divergence = cs_field_by_name("velocity_divergence");

  /* Boundary treatment */
  sc->bf_type = bf_type;

  /* Processing of the pressure boundary condition */
  sc->pressure_bc = cs_cdo_bc_face_define(CS_CDO_BC_HMG_NEUMANN, /* Default */
                                          true, /* Steady BC up to now */
                                          1,    /* Dimension */
                                          nsp->n_pressure_bc_defs,
                                          nsp->pressure_bc_defs,
                                          cs_shared_quant->n_b_faces);

  /* Set the way to enforce the Dirichlet BC on the velocity
   * "fixed_wall" means a no-slip BC
   */
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

  if (cs_shared_range_set != NULL) { /* Shared structures have to be freed */
    cs_range_set_destroy(&cs_shared_range_set);
    cs_interface_set_destroy(&cs_shared_interface_set);
    cs_matrix_structure_destroy(&cs_shared_matrix_structure);
    cs_matrix_assembler_destroy(&cs_shared_matrix_assembler);
  }

  /* Other pointers are only shared (i.e. not owner) */
  BFT_FREE(sc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Stokes or Oseen system with a CDO face-based scheme
 *         using a monolithic approach
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in]      nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_compute_steady(const cs_mesh_t            *mesh,
                                   const cs_navsto_param_t    *nsp,
                                   void                       *scheme_context)
{
  CS_UNUSED(nsp);

  cs_timer_t  t_start = cs_timer_time();

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces, n_cells = quant->n_cells;

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

  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;

  /* Build an array storing the Dirichlet values at faces. */
  cs_real_t  *dir_values = NULL;
  cs_cdofb_vecteq_setup_bc(t_cur, mesh, mom_eqp, mom_eqb, &dir_values);

  /* Initialize the local matrix */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_matrix_structure);

  /* Initialize the two right-hand side */
  cs_real_t  *mom_rhs = NULL, *mass_rhs = NULL;

  BFT_MALLOC(mom_rhs, 3*n_faces, cs_real_t);
  BFT_MALLOC(mass_rhs, n_cells, cs_real_t);

  /* Main loop on cells to define the linear system to solve */
  _steady_build(nsp, sc, dir_values, matrix, mom_rhs, mass_rhs);

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);

  /* End of the system building */
  cs_timer_t  t_bld_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_bld_end);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */
  _fields_to_previous(sc);

  /* Solve the linear system */
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eq->field_id, NULL);

  _solve_system(sles, matrix, mom_eqp,
                mom_eqc->face_values, /* velocity DoFs at faces */
                sc->pressure->val,    /* pressure DoFs at cells */
                mom_rhs, mass_rhs);

  /* Now compute/update the velocity and pressure fields */
  _update_fields(sc, mom_eq);

  /* Frees */
  cs_sles_free(sles);
  BFT_FREE(mom_rhs);
  BFT_FREE(mass_rhs);
  cs_matrix_destroy(&matrix);

  cs_timer_t  t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a monolithic approach and an Euler time scheme.
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_compute_implicit(const cs_mesh_t          *mesh,
                                     const cs_navsto_param_t  *nsp,
                                     void                     *scheme_context)
{
  CS_UNUSED(nsp);

  cs_timer_t  t_start = cs_timer_time();

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces, n_cells = quant->n_cells;

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

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_eval = ts->t_cur + ts->dt[0];

  /* Build an array storing the Dirichlet values at faces. */
  cs_real_t  *dir_values = NULL;
  cs_cdofb_vecteq_setup_bc(t_eval, mesh, mom_eqp, mom_eqb, &dir_values);

  /* Initialize the local matrix */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_matrix_structure);

  /* Initialize the two right-hand side */
  cs_real_t  *mom_rhs = NULL,  *mass_rhs = NULL;

  BFT_MALLOC(mom_rhs, 3*n_faces, cs_real_t);
  BFT_MALLOC(mass_rhs, n_cells, cs_real_t);

  /* Main loop on cells to define the linear system to solve */
  _implicit_euler_build(nsp, sc, dir_values, matrix, mom_rhs, mass_rhs);

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);

  /* End of the system building */
  cs_timer_t  t_bld_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_bld_end);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */
  _fields_to_previous(sc);

  /* Solve the linear system */
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eq->field_id, NULL);

  _solve_system(sles, matrix, mom_eqp,
                mom_eqc->face_values, /* velocity DoFs at faces */
                sc->pressure->val,    /* pressure DoFs at cells */
                mom_rhs, mass_rhs);

  /* Now compute/update the velocity and pressure fields */
  _update_fields(sc, mom_eq);

  /* Frees */
  cs_sles_free(sles);
  BFT_FREE(mom_rhs);
  BFT_FREE(mass_rhs);
  cs_matrix_destroy(&matrix);

  cs_timer_t  t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a monolithic approach and a theta time scheme.
 *
 * \param[in] mesh            pointer to a \ref cs_mesh_t structure
 * \param[in] nsp             pointer to a \ref cs_navsto_param_t structure
 * \param[in] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_compute_theta(const cs_mesh_t          *mesh,
                                  const cs_navsto_param_t  *nsp,
                                  void                     *scheme_context)
{
  CS_UNUSED(nsp);

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces, n_cells = quant->n_cells;

  /* Retrieve high-level structures */
  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t *)scheme_context;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;
  cs_real_t  *vel_c = sc->velocity->val;

  const cs_timer_t  t_start = cs_timer_time();

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const double  tcoef = 1 - mom_eqp->theta;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  /* Time_eval = (1-theta).t^n + theta.t^(n+1) = t^n + theta.dt
   * since t^(n+1) = t^n + dt
   */
  const cs_real_t  t_eval = t_cur + mom_eqp->theta*dt_cur;

  /* Detect the first call (in this case, we compute the initial source term)*/
  bool  compute_initial_source = false;
  if (ts->nt_cur == ts->nt_prev || ts->nt_prev == 0)
    compute_initial_source = true;

  /* Build an array storing the Dirichlet values at faces. */
  cs_real_t  *dir_values = NULL;
  cs_cdofb_vecteq_setup_bc(t_cur + dt_cur, mesh, mom_eqp, mom_eqb, &dir_values);

  /* Initialize the local matrix */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_matrix_structure);

  /* Initialize the two right-hand side */
  cs_real_t  *mom_rhs = NULL;
  cs_real_t  *mass_rhs = NULL;

  BFT_MALLOC(mom_rhs, 3*n_faces, cs_real_t);
  BFT_MALLOC(mass_rhs, n_cells, cs_real_t);

  /* Reset RHS arrays */
  _set_rhs_to_zero(mom_rhs, mass_rhs);

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(quant, connect, mom_eqp, mom_eqb, mom_eqc, mom_rhs, matrix, nsp,   \
         mass_rhs, mav, dir_values, vel_c, sc, compute_initial_source)  \
  firstprivate(t_eval, t_cur, dt_cur, tcoef, inv_dtcur)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(connect);
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

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
                                       dir_values, vel_c, t_eval,
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
      cs_cdofb_vecteq_advection_diffusion(t_eval, mom_eqp, mom_eqc, cm,
                                          csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion", csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);

      if (has_sourceterm) { /* SOURCE TERM
                             * =========== */
        if (compute_initial_source) { /* First time step */

          cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                     t_cur, tcoef,  /* time, scaling */
                                     cb, mom_eqb, csys);

        }
        else { /* Add the contribution of the previous time step */

          for (short int k = 0; k < 3; k++)
            csys->rhs[3*cm->n_fc + k] += tcoef*mom_eqc->source_terms[3*c_id+k];

        }

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   /* time,        scaling */
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
        cs_cell_sys_dump(">> Local system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _apply_remaining_bc(sc, mom_eqp, cm, csys, cb, &nsb, mass_rhs + c_id);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      _assemble(csys, cm, nsb.div_op,
                has_sourceterm, mav, mom_rhs, mom_eqc->source_terms);

    } /* Main loop on cells */

    /* Free temporary buffer */
    cs_cdofb_navsto_free_builder(&nsb);

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t_bld_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_bld_end);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */
  _fields_to_previous(sc);

  /* Solve the linear system */
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eq->field_id, NULL);

  _solve_system(sles, matrix, mom_eqp,
                mom_eqc->face_values, /* velocity DoFs at faces */
                sc->pressure->val,    /* pressure DoFs at cells */
                mom_rhs, mass_rhs);

  /* Now compute/update the velocity and pressure fields */
  _update_fields(sc, mom_eq);

  /* Frees */
  cs_sles_free(sles);
  BFT_FREE(mom_rhs);
  BFT_FREE(mass_rhs);
  cs_matrix_destroy(&matrix);

  cs_timer_t t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
