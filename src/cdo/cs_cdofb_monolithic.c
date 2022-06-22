/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it as one block (monolithic approach of the velocity-pressure
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
#include "cs_cdo_bc.h"
#include "cs_cdo_blas.h"
#include "cs_cdofb_priv.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_vecteq.h"
#include "cs_cdofb_navsto.h"
#include "cs_cdofb_monolithic_sles.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation_bc.h"
#include "cs_equation_builder.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_navsto_coupling.h"
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
 *        equations and solved it with a monolithic approach
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

static cs_sdm_t **cs_cdofb_monolithic_cw_mat = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current content of the Navier-Stokes related variables-fields
 *         to previous values
 *         Case of the monolithic coupling algorithm.
 *
 * \param[in, out] sc    scheme context (\ref cs_cdofb_monolithic_t struct.)
 * \param[in, out] cc    coupling context (\ref cs_navsto_monolithic_t struct.)
 */
/*----------------------------------------------------------------------------*/

static inline void
_mono_fields_to_previous(cs_cdofb_monolithic_t        *sc,
                         cs_navsto_monolithic_t       *cc)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Cell unknows: velocity, pressure) */

  cs_field_current_to_previous(sc->velocity);
  cs_field_current_to_previous(sc->pressure);

  /* Face unknows: mass flux and face velocity */

  memcpy(sc->mass_flux_array_pre, sc->mass_flux_array,
         quant->n_faces * sizeof(cs_real_t));

  cs_cdofb_vecteq_t  *mom_eqc
    = (cs_cdofb_vecteq_t *)cc->momentum->scheme_context;

  if (mom_eqc->face_values_pre != NULL)
    memcpy(mom_eqc->face_values_pre, mom_eqc->face_values,
           3 * quant->n_faces * sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Make sure that the enforcement is done (remove tolerance that may
 *         arise from the resolution)
 *         Case of a monolithic coupling algorithm.
 *
 * \param[in, out]  vel_f    velocity at faces to enforce
 */
/*----------------------------------------------------------------------------*/

static void
_mono_enforce_solid_face_velocity(cs_real_t           *vel_f)
{
  /* Enforcement of solid cells is always defined as follows for the momentum
   * equation:
   * CS_EQUATION_ENFORCE_BY_CELLS | CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE */

  cs_solid_selection_t  *solid = cs_solid_selection_get();

  if (solid->n_g_cells > 0) {

    for (cs_lnum_t f = 0; f < cs_shared_connect->n_faces[0]; f++) {

      if (solid->face_is_solid[f])
        for (int k = 0; k < 3; k++) vel_f[3*f+k] = 0.;

    } /* Loop on faces and search for faces with a tag */

  } /* There is at least one solid cell globally */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the velocity fields at cells and rescale the pressure field
 *         Case of a monolithic coupling algorithm.
 *
 * \param[in]       nsp      set of parameters for the Navier-Stokes system
 * \param[in, out]  sc       scheme context
 * \param[in, out]  mom_eqc  context of the momentum equation
 */
/*----------------------------------------------------------------------------*/

static void
_mono_update_related_cell_fields(const cs_navsto_param_t       *nsp,
                                 cs_cdofb_monolithic_t         *sc,
                                 cs_cdofb_vecteq_t             *mom_eqc)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_real_t  *vel_f = mom_eqc->face_values;
  cs_real_t  *vel_c = sc->velocity->val;

  /* Update the cell velocity
   * Compute values at cells vel_c from values at faces vel_f
   *     vel_c = acc^-1*(RHS - Acf*vel_f)
   */

  cs_static_condensation_recover_vector(cs_shared_connect->c2f,
                                        mom_eqc->rc_tilda,
                                        mom_eqc->acf_tilda,
                                        vel_f,              /* face values */
                                        vel_c);             /* cell values */

  /* Enforcement of solid cells is always defined as follows for the momentum
   * equation:
   * CS_EQUATION_ENFORCE_BY_CELLS | CS_EQUATION_ENFORCE_BY_REFERENCE_VALUE */

  cs_solid_selection_t  *solid = cs_solid_selection_get();

  for (cs_lnum_t i = 0; i < solid->n_cells; i++) {
    cs_real_t  *_vel = vel_c + 3*solid->cell_ids[i];
    for (int k = 0; k < 3; k++)
      _vel[k] = 0.;
  }

  /* Rescale pressure if needed */

  cs_field_t  *pr_fld = sc->pressure;

  if (sc->pressure_rescaling == CS_BOUNDARY_PRESSURE_RESCALING)
    cs_cdofb_navsto_rescale_pressure_to_ref(nsp, quant, pr_fld->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
  cs_dbg_darray_to_listing("CELL_VELOCITY", 3*quant->n_cells, vel_c, 9);
#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 3
  cs_dbg_darray_to_listing("FACE_VELOCITY", 3*quant->n_faces, vel_f, 9);
#endif
  cs_dbg_darray_to_listing("PRESSURE", quant->n_cells, pr_fld->val, 9);
#endif
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
  /* Compute the range set for an array of size 3*n_faces + n_cells
   * velocity is attached to faces (one for each component) and pressure
   * to cells
   *
   * Storage for the global numbering: Vel_X | Vel_Y | Vel_Z | Pressure */

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_mesh_t  *m = cs_shared_mesh;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_lnum_t  size = 3*n_faces + m->n_cells;

  assert(block->type == CS_CDO_SYSTEM_BLOCK_EXT);
  cs_cdo_system_xblock_t  *xb = block->block_pointer;

  /* 1. Build the interface set and the range set structures */

  cs_interface_set_t *ifs = connect->face_ifs;

  if (ifs != NULL)
    xb->interface_set = cs_interface_set_dup_blocks(ifs, n_faces, 3);

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

    /* A face-face entry corresponds to 3x3 block + the diagonal which is not
       taken into account in the face --> face connectivity. The B and Bt
       operators have the same sparsity. 3x1 entries for the c2f
       connectivity. This is multiply by two since one considers B and Bt. */

    int  n_entries = (end-start + 1)*9
                   + 6*(f2c->idx[frow_id+1]-f2c->idx[frow_id]);

    const cs_gnum_t  grow_ids[3]
      = { xb->range_set->g_id[frow_id],               /* x-component */
          xb->range_set->g_id[frow_id +   n_faces],   /* y-component */
          xb->range_set->g_id[frow_id + 2*n_faces] }; /* z-component */

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
        = { xb->range_set->g_id[fcol_id],               /* x-component */
            xb->range_set->g_id[fcol_id + n_faces],     /* y-component */
            xb->range_set->g_id[fcol_id + 2*n_faces] }; /* z-component */

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
      const cs_gnum_t  gcol_id = xb->range_set->g_id[3*n_faces + ccol_id];

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

    cs_matrix_assembler_add_g_ids(xb->matrix_assembler,
                                  n_entries, grows, gcols);
    assert(shift == n_entries);

  } /* Loop on face entities */

  if (add_pressure_diag) {

    const cs_gnum_t  *cell_g_ids = xb->range_set->g_id + 3*n_faces;

    cs_matrix_assembler_add_g_ids(xb->matrix_assembler,
                                  m->n_cells, cell_g_ids, cell_g_ids);

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
 *          Case of CDO-Fb schemes with a monolithic velocity-pressure coupling
 *
 * \param[in]      sc        pointer to a cs_cdofb_monolithic_t structure
 * \param[in]      mom_eqp   pointer to a cs_equation_param_t structure
 * \param[in]      cm        pointer to a cellwise view of the mesh
 * \param[in]      bf_type   type of boundary for the boundary face
 * \param[in]      diff_pty  pointer to a \cs_property_data_t struct. for diff.
 * \param[in, out] csys      pointer to a cellwise view of the system
 * \param[in, out] cb        pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_mono_apply_bc_partly(const cs_cdofb_monolithic_t   *sc,
                      const cs_equation_param_t     *mom_eqp,
                      const cs_cell_mesh_t          *cm,
                      const cs_boundary_type_t      *bf_type,
                      const cs_property_data_t      *diff_pty,
                      cs_cell_sys_t                 *csys,
                      cs_cell_builder_t             *cb)
{
  assert(cs_equation_param_has_diffusion(mom_eqp) == true);

  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Update the velocity-block and the right-hand side (part related to the
     * momentum equation). */

    /* Neumann boundary conditions:
     * The common practice is to define Phi_neu = - lambda * grad(u) . n_fc
     * An outward flux is a positive flux whereas an inward flux is negative
     * The minus just above implies the minus just below */

    if (csys->has_nhmg_neumann)
      for (short int i  = 0; i < 3*cm->n_fc; i++)
        csys->rhs[i] -= csys->neu_values[i];

    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */

      const short int  f = csys->_f_ids[i];

      if (bf_type[i] & CS_BOUNDARY_IMPOSED_VEL) {
        if (mom_eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            mom_eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          sc->apply_velocity_inlet(f, mom_eqp, cm, diff_pty, cb, csys);
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_WALL) {
        if (mom_eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            mom_eqp->default_enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          if (bf_type[i] & CS_BOUNDARY_SLIDING_WALL)
            sc->apply_sliding_wall(f, mom_eqp, cm, diff_pty, cb, csys);
          else
            sc->apply_fixed_wall(f, mom_eqp, cm, diff_pty, cb, csys);
        }
      }

      else if (bf_type[i] & CS_BOUNDARY_SYMMETRY) {

        /* Always weakly enforce the symmetry constraint on the
           velocity-block */

        sc->apply_symmetry(f, mom_eqp, cm, diff_pty, cb, csys);

      }

      /* default: nothing to do (case of a "natural" outlet) */

    } /* Loop on boundary faces */

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
    if (cs_dbg_cw_test(mom_eqp, cm, csys))
      cs_cell_sys_dump(">> Cell system matrix after partial BC enforcement",
                       csys);
#endif
  } /* Boundary cell */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the boundary conditions to the local system when this should
 *         be done after the static condensation. Apply the internal enforcement
 *         Case of CDO-Fb schemes with a monolithic velocity-pressure coupling
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      eqb       pointer to a cs_equation_builder_t structure
 * \param[in]      cm        pointer to a cellwise view of the mesh
 * \param[in]      diff_pty  pointer to a \cs_property_data_t struct. for diff.
 * \param[in, out] sc        pointer to a cs_cdofb_monolithic_t structure
 * \param[in, out] csys      pointer to a cellwise view of the system
 * \param[in, out] cb        pointer to a cellwise builder
 * \param[in, out] nsb       builder structure for the NavSto system
 */
/*----------------------------------------------------------------------------*/

static void
_mono_apply_remaining_bc(const cs_equation_param_t     *eqp,
                         const cs_equation_builder_t   *eqb,
                         const cs_cell_mesh_t          *cm,
                         const cs_property_data_t      *diff_pty,
                         cs_cdofb_monolithic_t         *sc,
                         cs_cell_sys_t                 *csys,
                         cs_cell_builder_t             *cb,
                         cs_cdofb_navsto_builder_t     *nsb)
{
  if (cb->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    const cs_boundary_type_t  *bf_type = nsb->bf_type;
    cs_real_t  *div_op = nsb->div_op;
    cs_real_t  *mass_rhs = sc->system_helper->rhs + 3*cs_shared_quant->n_faces;

    /* Update the divergence operator and the right-hand side related to the
     * mass equation.
     * Enforcement of Dirichlet BC in a stronger way if this is the choice */

    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */

      const short int  f = csys->_f_ids[i];

      if (bf_type[i] & CS_BOUNDARY_IMPOSED_VEL) {

        /* Update mass RHS (constrain on the velocity divergence) from the
           knowledge of the boundary face velocity */

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
          sc->apply_velocity_inlet(f, eqp, cm, diff_pty, cb, csys);

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
            sc->apply_sliding_wall(f, eqp, cm, diff_pty, cb, csys);
          else
            sc->apply_fixed_wall(f, eqp, cm, diff_pty, cb, csys);
        }

      }
      else if (bf_type[i] & CS_BOUNDARY_SYMMETRY) {

        /* No need to update the mass RHS since there is no mass flux.
         * Weak-enforcement for the velocity-block (cf. _mono_apply_bc_partly)
         * Strong enforcement of u.n (--> dp/dn = 0) on the divergence */

        for (int k = 0; k < 3; k++) div_op[3*f+k] = 0;

      }

      /* default: nothing to do (case of a "natural" outlet) */

    } /* Loop over boundary faces */

  } /* This is a boundary cell */

  if (cs_equation_param_has_internal_enforcement(eqp)) {

    /* Internal enforcement of DoFs: Update csys (matrix and rhs) */

    cs_equation_builder_enforce_block_dofs(eqb, cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
    if (cs_dbg_cw_test(eqp, cm, csys))
      cs_cell_sys_dump("\n>> Cell system after the internal enforcement",
                       csys);
#endif
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb schemes
 *         Shares similarities with cs_cdo_assembly_block_matrix()
 *
 * \param[in]       csys              pointer to a cs_cell_sys_t structure
 * \param[in]       cm                pointer to a cs_cell_mesh_t structure
 * \param[in]       div_op            array with the divergence op. values
 * \param[in, out]  sc                pointer to scheme context structure
 * \param[in, out]  eqc               context structure for a vector-valued Fb
 * \param[in, out]  asb               pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_assembly_by_blocks(const cs_cell_sys_t        *csys,
                    const cs_cell_mesh_t       *cm,
                    const cs_real_t            *div_op,
                    cs_cdofb_monolithic_t      *sc,
                    cs_cdofb_vecteq_t          *eqc,
                    cs_cdo_assembly_t          *asb)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  int  t_id = 0;
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
#endif

  cs_cdo_system_helper_t  *sh = sc->system_helper;
  cs_cdo_system_block_t  *b = sh->blocks[0];
  assert(b->type == CS_CDO_SYSTEM_BLOCK_SPLIT);
  cs_cdo_system_sblock_t  *sb = b->block_pointer;
  cs_cdofb_monolithic_sles_t  *msles = sc->msles;
  cs_sdm_t  *cw_mat = cs_cdofb_monolithic_cw_mat[t_id];

  /* Convert csys->mat into a block view by component */

  cs_sdm_block_33_to_xyz(csys->mat, cw_mat);

  /* Assembly is performed block by block */

  for (int i = 0; i < sh->n_col_blocks; i++) {
    for (int j = 0; j < sh->n_col_blocks; j++) {

      cs_matrix_assembler_values_t   *mav = sb->mav_array[3*i+j];
      cs_sdm_t  *m_ij = cs_sdm_get_block(cw_mat, i, j);

      sb->assembly_func(m_ij, cm->f_ids, sb->range_set, asb, mav);

    } /* Loop on blocks (j) */
  } /* Loop on blocks (i) */

  /* RHS assembly */

# pragma omp critical
  {
    for (short int f = 0; f < cm->n_fc; f++) {
      sh->rhs_array[0][cm->f_ids[f]] += csys->rhs[3*f];
      sh->rhs_array[1][cm->f_ids[f]] += csys->rhs[3*f+1];
      sh->rhs_array[2][cm->f_ids[f]] += csys->rhs[3*f+2];
    }
  }

  /* Reset the value of the source term for the cell DoF
     Source term is only hold by the cell DoF in face-based schemes */

  if (eqc->source_terms != NULL) {
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
 * \param[in, out]  sc                pointer to scheme context structure
 * \param[in, out]  eqc               context structure for a vector-valued Fb
 * \param[in, out]  asb               pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_velocity_full_assembly(const cs_cell_sys_t            *csys,
                        const cs_cell_mesh_t           *cm,
                        const cs_real_t                *div_op,
                        cs_cdofb_monolithic_t          *sc,
                        cs_cdofb_vecteq_t              *eqc,
                        cs_cdo_assembly_t              *asb)
{
  const short int  n_f = cm->n_fc;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_cdo_system_helper_t  *sh = sc->system_helper;
  cs_real_t  *_div = sc->msles->div_op + 3*connect->c2f->idx[cm->c_id];

  /* 1. Store divergence operator in non assembly
   *    Take into account solid zone where DoF is set to zero */
  /* ======================================================== */

  if (csys->has_internal_enforcement) {

    for (int i = 0; i < 3*n_f; i++) {
      if (csys->dof_is_forced[i])
        _div[i] = 0.; /* The velocity-block set the value of this DoF */
      else
        _div[i] = div_op[i];
    }

  }
  else
    memcpy(_div, div_op, 3*n_f*sizeof(cs_real_t));

  /* 1. Matrix assembly
   * ================== */

  if (sc->msles->graddiv_coef > 0.) {

    cs_real_t  gamma = sc->msles->graddiv_coef / cm->vol_c;
    cs_cdofb_navsto_add_grad_div(cm->n_fc, gamma, _div, csys->mat);

  }

  cs_cdofb_vecteq_assembly(csys, sh->blocks[0], sh->rhs, eqc, asb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb schemes
 *         Shares similarities with cs_cdo_assembly_block_matrix()
 *
 * \param[in]      csys             pointer to a cs_cell_sys_t structure
 * \param[in]      cm               pointer to a cs_cell_mesh_t structure
 * \param[in]      div_op           array with the divergence op. values
 * \param[in, out] sc               pointer to scheme context structure
 * \param[in, out] eqc              context structure for a vector-valued Fb
 * \param[in, out] asb              pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_full_assembly(const cs_cell_sys_t            *csys,
               const cs_cell_mesh_t           *cm,
               const cs_real_t                *div_op,
               cs_cdofb_monolithic_t          *sc,
               cs_cdofb_vecteq_t              *eqc,
               cs_cdo_assembly_t              *asb)
{
  CS_UNUSED(asb);

  const short int  n_f = cm->n_fc;
  const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
  const cs_sdm_t  *m = csys->mat;
  const cs_sdm_block_t  *bd = m->block_desc;

  assert(m->flag & CS_SDM_BY_BLOCK);
  assert(m->block_desc != NULL);
  assert(bd->n_row_blocks == bd->n_col_blocks);
  assert(bd->n_row_blocks == n_f);

  cs_cdo_system_helper_t  *sh = sc->system_helper;
  cs_cdo_system_block_t  *b = sh->blocks[0];
  assert(b->type == CS_CDO_SYSTEM_BLOCK_EXT);
  cs_cdo_system_xblock_t  *xb = b->block_pointer;

  cs_real_t  *eqc_st = eqc->source_terms;
  cs_real_t  *rhs = sh->rhs;

  const cs_range_set_t  *rset = xb->range_set;

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
            cs_matrix_assembler_values_add_g(xb->mav, bufsize,
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
        cs_matrix_assembler_values_add_g(xb->mav, bufsize,
                                         r_gids, c_gids, values);
        bufsize = 0;
      }

      /* Its transposed B_x, B_y, B_z */

      r_gids[bufsize] = p_gid;
      c_gids[bufsize] = bi_gids[ii];
      values[bufsize] = div_op[3*bi+ii];
      bufsize += 1;

      if (bufsize == CS_CDO_ASSEMBLE_BUF_SIZE) {
#       pragma omp critical
        cs_matrix_assembler_values_add_g(xb->mav, bufsize,
                                         r_gids, c_gids, values);
        bufsize = 0;
      }

    } /* Loop on components */

  } /* Loop on row blocks (bi) */

  if (bufsize > 0) {
#   pragma omp critical
    cs_matrix_assembler_values_add_g(xb->mav, bufsize,
                                     r_gids, c_gids, values);
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

  if (eqc_st != NULL)
    for (int k = 0; k < 3; k++)
      eqc_st[3*cm->c_id + k] = csys->source[3*n_f + k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb schemes
 *         Shares similarities with cs_cdo_assembly_block_matrix()
 *         This is a _full_assembly() function plus a specific treatment for
 *         the enforcement of pressure to zero inside the solid region(s).
 *
 * \param[in]      csys             pointer to a cs_cell_sys_t structure
 * \param[in]      cm               pointer to a cs_cell_mesh_t structure
 * \param[in]      div_op           array with the divergence op. values
 * \param[in, out] sc               pointer to scheme context structure
 * \param[in, out] eqc              context structure for a vector-valued Fb
 * \param[in, out] asb              pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_solidification_full_assembly(const cs_cell_sys_t            *csys,
                              const cs_cell_mesh_t           *cm,
                              const cs_real_t                *div_op,
                              cs_cdofb_monolithic_t          *sc,
                              cs_cdofb_vecteq_t              *eqc,
                              cs_cdo_assembly_t              *asb)
{
  /* 1. First part shared with the assembly of the full saddle-point problem
   * ======================================================================= */

  _full_assembly(csys, cm, div_op, sc, eqc, asb);

  /* 2. Treatment of the solid zone(s)
   * ================================= */

  cs_solid_selection_t  *solid = cs_solid_selection_get();

  assert(solid != NULL);
  if (solid->cell_is_solid == NULL) /* No solid cell in the whole domain */
    return;

  if (solid->cell_is_solid[cm->c_id]) { /* This cell is solid */

    cs_cdo_system_block_t  *b = sc->system_helper->blocks[0];
    assert(b->type == CS_CDO_SYSTEM_BLOCK_EXT);
    cs_cdo_system_xblock_t  *xb = b->block_pointer;

    const cs_range_set_t  *rset = xb->range_set;

    /* Modify the row related to the cell pressure: Add a +1 to the diagonal
       and reset the divergence operator (rhs is equal to zero in this part by
       default so that there is no need to modify the rhs) */

    const cs_lnum_t  n_faces = cs_shared_quant->n_faces;
    const cs_gnum_t  p_gid = rset->g_id[3*n_faces + cm->c_id];

    cs_gnum_t  r_gids[CS_CDO_ASSEMBLE_BUF_SIZE];
    cs_gnum_t  c_gids[CS_CDO_ASSEMBLE_BUF_SIZE];
    cs_real_t  values[CS_CDO_ASSEMBLE_BUF_SIZE];
    int  bufsize = 0;

    for (int f = 0; f < cm->n_fc; f++) {

      const cs_lnum_t  f_id = cm->f_ids[f];
      const cs_gnum_t  uf_gids[3]
        = { rset->g_id[f_id],              /* x-component */
            rset->g_id[f_id +   n_faces],  /* y-component */
            rset->g_id[f_id + 2*n_faces]}; /* z-component */

      for (int k = 0; k < 3; k++) {

        r_gids[bufsize] = uf_gids[k];
        c_gids[bufsize] = uf_gids[k];
        values[bufsize] = 1e9;            /* Penalized the face velocity */
        bufsize += 1;

        r_gids[bufsize] = p_gid;
        c_gids[bufsize] = uf_gids[k];
        values[bufsize] = -div_op[3*f+k]; /* Reset the value for the div. op. */
        bufsize += 1;

      }

    } /* Loop on cell faces */

    r_gids[bufsize] = p_gid;
    c_gids[bufsize] = p_gid;
    values[bufsize] = 1;                  /* Set 1 on the diagonal */
    bufsize += 1;

    assert(bufsize <= CS_CDO_ASSEMBLE_BUF_SIZE);

#   pragma omp critical
    cs_matrix_assembler_values_add_g(xb->mav, bufsize, r_gids, c_gids, values);

  } /* This cell is tag as solid */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a vector-valued system obtained
 *         with CDO-Fb schemes
 *         Shares similarities with cs_cdo_assembly_block_matrix()
 *         Specify to the Notay's transformation of the saddle-point system.
 *         Rely on the _full_assembly() function
 *
 * \param[in]      csys             pointer to a cs_cell_sys_t structure
 * \param[in]      cm               pointer to a cs_cell_mesh_t structure
 * \param[in]      div_op           array with the divergence op. values
 * \param[in, out] sc               pointer to scheme context structure
 * \param[in, out] eqc              context structure for a vector-valued Fb
 * \param[in, out] asb              pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_notay_full_assembly(const cs_cell_sys_t            *csys,
                     const cs_cell_mesh_t           *cm,
                     const cs_real_t                *div_op,
                     cs_cdofb_monolithic_t          *sc,
                     cs_cdofb_vecteq_t              *eqc,
                     cs_cdo_assembly_t              *asb)
{
  /* 1. First part shared with the assembly of the full saddle-point problem
   * ======================================================================= */

  _full_assembly(csys, cm, div_op, sc, eqc, asb);

  /* 2. Store divergence operator in non assembly
   * ============================================ */

  cs_real_t *_div = sc->msles->div_op + 3*cs_shared_connect->c2f->idx[cm->c_id];

  if (csys->has_internal_enforcement) {

    for (int i = 0; i < 3*cm->n_fc; i++) {
      if (csys->dof_is_forced[i])
        _div[i] = 0.; /* The velocity-block set the value of this DoF */
      else
        _div[i] = div_op[i];
    }

  }
  else
    memcpy(_div, div_op, 3*cm->n_fc*sizeof(cs_real_t));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         steady-state case
 *
 * \param[in]      nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]      vel_f_pre    velocity face DoFs of the previous time step
 * \param[in]      vel_c_pre    velocity cell DoFs of the previous time step
 * \param[in]      vel_f_nm1    NULL (for unsteady computations)
 * \param[in]      vel_c_nm1    NULL (for unsteady computations)
 * \param[in, out] sc           pointer to the scheme context
 */
/*----------------------------------------------------------------------------*/

static void
_steady_build(const cs_navsto_param_t      *nsp,
              const cs_real_t               vel_f_pre[],
              const cs_real_t               vel_c_pre[],
              const cs_real_t               vel_f_nm1[],
              const cs_real_t               vel_c_nm1[],
              cs_cdofb_monolithic_t        *sc)
{
  CS_UNUSED(vel_f_nm1);
  CS_UNUSED(vel_c_nm1);

  /* Retrieve shared structures */

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_shared_time_step;

  /* Retrieve high-level structures */

  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

#if defined(DEBUG) && !defined(NDEBUG)
  if (quant->n_b_faces > 0)
    assert(mom_eqb->dir_values != NULL);
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

    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(nsp,
                                                                    connect);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *diff_hodge =
      (mom_eqc->diffusion_hodge == NULL) ? NULL:mom_eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (mom_eqc->mass_hodge == NULL) ? NULL:mom_eqc->mass_hodge[t_id];

    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    /* Set times at which one evaluates quantities when needed */

    cb->t_pty_eval = ts->t_cur; /* Dummy parameter if really steady */
    cb->t_bc_eval = ts->t_cur;  /* Dummy parameter if really steady */
    cb->t_st_eval = ts->t_cur;  /* Dummy parameter if really steady */

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(mom_eqp, mom_eqb, diff_hodge, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag,
                                                            mom_eqb),
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
       *
       * Set the local (i.e. cellwise) structures for the current cell
       */

      cs_cdofb_vecteq_init_cell_system(cm, mom_eqp, mom_eqb,
                                       vel_f_pre, vel_c_pre,
                                       NULL, NULL, /* no n-1 state is given */
                                       csys, cb);

      /* 1- SETUP THE NAVSTO LOCAL BUILDER
       * =================================
       * - Set the type of boundary
       * - Set the pressure boundary conditions (if required)
       * - Define the divergence operator used in the linear system (div_op is
       *   equal to minus the divergence)
       */
      cs_cdofb_navsto_define_builder(cb->t_bc_eval, nsp, cm, csys,
                                     sc->pressure_bc, sc->bf_type, &nsb);

      /* 2- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_conv_diff_reac(mom_eqp, mom_eqb, mom_eqc, cm,
                                     mass_hodge, diff_hodge, csys, cb);

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      if (cs_equation_param_has_sourceterm(mom_eqp))
        cs_cdofb_vecteq_sourceterm(cm, mom_eqp, cb->t_st_eval,
                                   1., /* time scaling */
                                   mass_hodge,
                                   cb, mom_eqb, csys);

      /* Gravity effects and/or Boussinesq approximation rely on another
         strategy than classical source term. The treatment is more compatible
         with the pressure gradient by doing so. */

      if (sc->add_gravity_term != NULL)
        sc->add_gravity_term(nsp, cm, &nsb, csys);

      /* 3b- OTHER RHS CONTRIBUTIONS
       * ===========================
       *
       * First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _mono_apply_bc_partly(sc, mom_eqp, cm, nsb.bf_type, diff_hodge->pty_data,
                            csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system matrix before condensation", csys);
#endif

      /* 5- STATIC CONDENSATION
       * ======================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */

      cs_static_condensation_vector_eq(connect->c2f,
                                       mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _mono_apply_remaining_bc(mom_eqp, mom_eqb, cm, diff_hodge->pty_data,
                               sc, csys, cb, &nsb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> (FINAL) Cell system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      sc->assemble(csys, cm, nsb.div_op, sc, mom_eqc, asb);

    } /* Main loop on cells */

    /* Free temporary buffer */

    cs_cdofb_navsto_free_builder(&nsb);

  } /* End of the OpenMP Block */

  cs_cdo_system_helper_finalize_assembly(sc->system_helper);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         case of an implicit Euler time scheme
 *
 * \param[in]      nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]      vel_f_n      velocity face DoFs at time step n
 * \param[in]      vel_c_n      velocity cell DoFs at time step n
 * \param[in]      vel_f_nm1    NULL (not needed for this time scheme)
 * \param[in]      vel_c_nm1    NULL (not needed for this time scheme)
 * \param[in, out] sc           pointer to the scheme context
 */
/*----------------------------------------------------------------------------*/

static void
_implicit_euler_build(const cs_navsto_param_t  *nsp,
                      const cs_real_t           vel_f_n[],
                      const cs_real_t           vel_c_n[],
                      const cs_real_t           vel_f_nm1[],
                      const cs_real_t           vel_c_nm1[],
                      cs_cdofb_monolithic_t    *sc)
{
  CS_UNUSED(vel_f_nm1);
  CS_UNUSED(vel_c_nm1);

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

#if defined(DEBUG) && !defined(NDEBUG)
  if (quant->n_b_faces > 0)
    assert(mom_eqb->dir_values != NULL);
#endif

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_real_t  t_eval = ts->t_cur + ts->dt[0];
    const cs_real_t  inv_dtcur = 1./ts->dt[0];

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(nsp,
                                                                    connect);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *diff_hodge =
      (mom_eqc->diffusion_hodge == NULL) ? NULL:mom_eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (mom_eqc->mass_hodge == NULL) ? NULL:mom_eqc->mass_hodge[t_id];

    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    /* Set times at which one evaluates quantities when needed */

    cb->t_pty_eval = t_eval;
    cb->t_bc_eval = t_eval;
    cb->t_st_eval = t_eval;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(mom_eqp, mom_eqb, diff_hodge, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag,
                                                            mom_eqb),
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

      cs_cdofb_vecteq_init_cell_system(cm, mom_eqp, mom_eqb,
                                       vel_f_n, vel_c_n,
                                       NULL, NULL, /* no n-1 state is given */
                                       csys, cb);

      /* 1- SETUP THE NAVSTO LOCAL BUILDER *
       * ================================= *
       * - Set the type of boundary
       * - Set the pressure boundary conditions (if required)
       * - Define the divergence operator used in the linear system (div_op is
       *   equal to minus the divergence)
       */

      cs_cdofb_navsto_define_builder(cb->t_bc_eval, nsp, cm, csys,
                                     sc->pressure_bc, sc->bf_type, &nsb);

      /* 2- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_conv_diff_reac(mom_eqp, mom_eqb, mom_eqc, cm,
                                     mass_hodge, diff_hodge, csys, cb);

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      if (cs_equation_param_has_sourceterm(mom_eqp))
        cs_cdofb_vecteq_sourceterm(cm, mom_eqp, cb->t_st_eval,
                                   1., /* time scaling */
                                   mass_hodge,
                                   cb, mom_eqb, csys);

      /* Gravity effects and/or Boussinesq approximation rely on another
         strategy than classical source term. The treatment is more compatible
         with the pressure gradient by doing so. */

      if (sc->add_gravity_term != NULL)
        sc->add_gravity_term(nsp, cm, &nsb, csys);

      /* 3b- OTHER RHS CONTRIBUTIONS *
       * =========================== *
       *
       * First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _mono_apply_bc_partly(sc, mom_eqp, cm, nsb.bf_type, diff_hodge->pty_data,
                            csys, cb);

      /* 4- TIME CONTRIBUTION (mass lumping or voronoï) */
      /* ==================== */

      if (!(mom_eqb->time_pty_uniform))
        cb->tpty_val = cs_property_value_in_cell(cm,
                                                 mom_eqp->time_property,
                                                 cb->t_pty_eval);

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system matrix before static condensation",
                         csys);
#endif

      /* 5- STATIC CONDENSATION
       * ======================
       * Static condensation of the local system matrix of size n_fc + 1 into
       * a matrix of size n_fc.
       * Store data in rc_tilda and acf_tilda to compute the values at cell
       * centers after solving the system */

      cs_static_condensation_vector_eq(connect->c2f,
                                       mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _mono_apply_remaining_bc(mom_eqp, mom_eqb, cm, diff_hodge->pty_data,
                               sc, csys, cb, &nsb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> (FINAL) Cell system matrix", csys);
#endif

      /* ASSEMBLY PROCESS */
      /* ================ */

      sc->assemble(csys, cm, nsb.div_op, sc, mom_eqc, asb);

    } /* Main loop on cells */

    /* Free temporary buffer */

    cs_cdofb_navsto_free_builder(&nsb);

  } /* OPENMP Block */

  cs_cdo_system_helper_finalize_assembly(sc->system_helper);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system for Stokes, Oseen or Navier-Stokes in the
 *         case of a theta time scheme
 *
 * \param[in]      nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in]      vel_f_n      velocity face DoFs at time step n
 * \param[in]      vel_c_n      velocity cell DoFs at time step n
 * \param[in]      vel_f_nm1    velocity face DoFs at time step n-1 or NULL
 * \param[in]      vel_c_nm1    velocity cell DoFs at time step n-1 or NULL
 * \param[in, out] sc           pointer to the scheme context
 */
/*----------------------------------------------------------------------------*/

static void
_theta_scheme_build(const cs_navsto_param_t  *nsp,
                    const cs_real_t           vel_f_n[],
                    const cs_real_t           vel_c_n[],
                    const cs_real_t           vel_f_nm1[],
                    const cs_real_t           vel_c_nm1[],
                    cs_cdofb_monolithic_t    *sc)
{
  CS_UNUSED(vel_f_nm1);
  CS_UNUSED(vel_c_nm1);

  /* Retrieve high-level structures */

  cs_navsto_monolithic_t  *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;

  /* Retrieve shared structures */

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t *ts = cs_shared_time_step;

#if defined(DEBUG) && !defined(NDEBUG)
  if (quant->n_b_faces > 0)
    assert(mom_eqb->dir_values != NULL);
#endif

  /* Detect the first call (in this case, we compute the initial source term)*/

  bool  compute_initial_source = false;
  if (ts->nt_cur == ts->nt_prev || ts->nt_prev == 0)
    compute_initial_source = true;

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)
  {
#if defined(HAVE_OPENMP) /* Determine the default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_real_t  t_cur = ts->t_cur;
    const cs_real_t  dt_cur = ts->dt[0];
    const double  tcoef = 1 - mom_eqp->theta;
    const cs_real_t  inv_dtcur = 1./dt_cur;
    const cs_real_t  t_eval = t_cur + mom_eqp->theta*dt_cur;

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */

    cs_cdofb_navsto_builder_t  nsb = cs_cdofb_navsto_create_builder(nsp,
                                                                    connect);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cdo_assembly_t  *asb = cs_cdo_assembly_get(t_id);
    cs_hodge_t  *diff_hodge =
      (mom_eqc->diffusion_hodge == NULL) ? NULL:mom_eqc->diffusion_hodge[t_id];
    cs_hodge_t  *mass_hodge =
      (mom_eqc->mass_hodge == NULL) ? NULL : mom_eqc->mass_hodge[t_id];

    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;

    cs_cdofb_vecteq_get(&csys, &cb);

    /* Set times at which one evaluates quantities when needed */
    /* Time_eval = (1-theta).t^n + theta.t^(n+1) = t^n + theta.dt
     * since t^(n+1) = t^n + dt
     */

    cb->t_pty_eval = t_eval;
    cb->t_bc_eval = t_cur + dt_cur;
    cb->t_st_eval = t_cur + dt_cur;

    /* Initialization of the values of properties */

    cs_equation_builder_init_properties(mom_eqp, mom_eqb, diff_hodge, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the current cell flag */

      cb->cell_flag = connect->cell_flag[c_id];

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id,
                         cs_equation_builder_cell_mesh_flag(cb->cell_flag,
                                                            mom_eqb),
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

      cs_cdofb_vecteq_init_cell_system(cm, mom_eqp, mom_eqb,
                                       vel_f_n, vel_c_n,
                                       NULL, NULL, /* no n-1 state is given */
                                       csys, cb);

      /* 1- SETUP THE NAVSTO LOCAL BUILDER *
       * ================================= *
       * - Set the type of boundary
       * - Set the pressure boundary conditions (if required)
       * - Define the divergence operator used in the linear system (div_op is
       *   equal to minus the divergence)
       */

      cs_cdofb_navsto_define_builder(cb->t_bc_eval, nsp, cm, csys,
                                     sc->pressure_bc, sc->bf_type,
                                     &nsb);

      /* 2- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_conv_diff_reac(mom_eqp, mom_eqb, mom_eqc, cm,
                                     mass_hodge, diff_hodge, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system after advection/diffusion/reaction",
                         csys);
#endif

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      if (cs_equation_param_has_sourceterm(mom_eqp)) {

        if (compute_initial_source)
          /* First time step: Compute source term at a previous step */
          cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                     t_cur, tcoef,  /* time, scaling */
                                     mass_hodge,
                                     cb, mom_eqb, csys);

        else /* Add the contribution of the previous time step */
          for (short int k = 0; k < 3; k++)
            csys->rhs[3*cm->n_fc + k] += tcoef*mom_eqc->source_terms[3*c_id+k];

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   /* time       , scaling */
                                   cb->t_st_eval, mom_eqp->theta,
                                   mass_hodge,
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* Gravity effects and/or Boussinesq approximation rely on another
         strategy than classical source term. The treatment is more compatible
         with the pressure gradient by doing so. */

      if (sc->add_gravity_term != NULL)
        sc->add_gravity_term(nsp, cm, &nsb, csys);

      /* 3b- OTHER RHS CONTRIBUTIONS *
       * =========================== *
       *
       * First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _mono_apply_bc_partly(sc, mom_eqp, cm, nsb.bf_type, diff_hodge->pty_data,
                            csys, cb);

      /* 4- UNSTEADY TERM + TIME SCHEME
       * ============================== */

      if (!(mom_eqb->time_pty_uniform))
        cb->tpty_val = cs_property_value_in_cell(cm,
                                                 mom_eqp->time_property,
                                                 cb->t_pty_eval);

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
        }

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
                                       mom_eqc->rc_tilda, mom_eqc->acf_tilda,
                                       cb, csys);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> Cell system matrix after static condensation",
                         csys);
#endif

      /* 6- Remaining part of BOUNDARY CONDITIONS
       * ======================================== */

      _mono_apply_remaining_bc(mom_eqp, mom_eqb, cm, diff_hodge->pty_data,
                               sc, csys, cb, &nsb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump("\n>> (FINAL) Cell system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      sc->assemble(csys, cm, nsb.div_op, sc, mom_eqc, asb);

    } /* Main loop on cells */

    /* Free temporary buffer */

    cs_cdofb_navsto_free_builder(&nsb);

  } /* OPENMP Block */

  cs_cdo_system_helper_finalize_assembly(sc->system_helper);
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
cs_cdofb_monolithic_init_sharing(const cs_navsto_param_t       *nsp,
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

  if (nsp->sles_param->strategy == CS_NAVSTO_SLES_BY_BLOCKS) {

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

  /* SLES needs these structures for advanced PETSc hooks */

  cs_cdofb_monolithic_sles_init_sharing(connect, quant);
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

  switch (nsp->sles_param->strategy) {

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
 * \param[in] nsp          pointer to a \ref cs_navsto_param_t structure
 * \param[in] adv_field    pointer to \ref cs_adv_field_t structure
 * \param[in] mflux        current values of the mass flux across primal faces
 * \param[in] mflux_pre    current values of the mass flux across primal faces
 * \param[in] bf_type      type of boundary for each boundary face
 * \param[in] cc_context   pointer to a \ref cs_navsto_monolithic_t structure
 *
 * \return a pointer to a new allocated \ref cs_cdofb_monolithic_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdofb_monolithic_init_scheme_context(const cs_navsto_param_t  *nsp,
                                        cs_adv_field_t           *adv_field,
                                        cs_real_t                *mflux,
                                        cs_real_t                *mflux_pre,
                                        cs_boundary_type_t       *bf_type,
                                        void                     *cc_context)
{
  assert(nsp != NULL && cc_context != NULL);
  if (nsp->space_scheme != CS_SPACE_SCHEME_CDOFB)
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.\n", __func__);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  /* Navier-Stokes scheme context (SC) */

  cs_cdofb_monolithic_t  *sc = NULL;

  BFT_MALLOC(sc, 1, cs_cdofb_monolithic_t);

  /* Cast the coupling context (CC) */

  cs_navsto_monolithic_t  *cc = (cs_navsto_monolithic_t  *)cc_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;

  /* Quantities shared with the cs_navsto_system_t structure */

  sc->coupling_context = cc;
  sc->adv_field = adv_field;
  sc->mass_flux_array = mflux;
  sc->mass_flux_array_pre = mflux_pre;

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

  sc->pressure_bc = cs_cdo_bc_face_define(CS_PARAM_BC_HMG_NEUMANN, /* Default */
                                          true,      /* Steady BC up to now */
                                          1,         /* Dimension */
                                          nsp->n_pressure_bc_defs,
                                          nsp->pressure_bc_defs,
                                          quant->n_b_faces);

  sc->pressure_rescaling =
    cs_boundary_need_pressure_rescaling(quant->n_b_faces, bf_type);

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

  /* Source terms induced by the gravity effect */

  cs_cdofb_navsto_set_gravity_func(nsp, &(sc->add_gravity_term));

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

  case CS_TIME_SCHEME_BDF2:
  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid time scheme.", __func__);

  } /* Switch on time schme */

  /* Handle the resolution of a saddle-point system */

  cs_cdo_system_helper_t  *sh = NULL;
  cs_cdofb_monolithic_sles_t  *msles =
    cs_cdofb_monolithic_sles_create(quant->n_faces, quant->n_cells);

  /* Define the layout of the system and how to assemble the system. It depends
     on the strategy to solve the saddle-point problem */

  cs_navsto_sles_t  strategy = nsp->sles_param->strategy;

  switch (strategy) {

  case CS_NAVSTO_SLES_BY_BLOCKS:
    {
      sc->assemble = _assembly_by_blocks;

      cs_lnum_t  block_sizes[2];
      block_sizes[0] = 3*quant->n_faces, block_sizes[1] = quant->n_cells;

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       2,
                                       block_sizes,
                                       2);

      /* Choose the right class of matrix to avoid copy.
       * The way to perform the assembly may change if an external librairy is
       * used for solving the linear system */

      cs_cdo_system_matrix_class_t  matclass;

      switch (mom_eqp->sles_param->solver_class) {

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

      /* Add a first block for the (0,0) block and then define the underpinning
         structures */

      cs_cdo_system_add_sblock(sh, 0,                /* block id */
                               matclass,             /* class of matrices */
                               cs_flag_primal_face , /* location */
                               quant->n_faces,       /* n_elements */
                               3);                   /* stride */

      cs_cdo_system_build_block(sh, 0); /* build structures */

      cs_cdo_system_block_t  *bdiv =
        cs_cdo_system_add_ublock(sh, 1,               /* block_id */
                                 connect->c2f,        /* adjacency */
                                 cs_flag_primal_face, /* col.location */
                                 quant->n_faces,      /* n_elements */
                                 3,                   /* stride */
                                 true);               /* interlaced */

      cs_cdo_system_build_block(sh, 1); /* build range_set and interface_set
                                           structures */

      BFT_MALLOC(msles->div_op, 3*connect->c2f->idx[quant->n_cells], cs_real_t);

      cs_cdo_system_ublock_t  *b_ub = bdiv->block_pointer;

      b_ub->values = msles->div_op;              /* shared pointer */
      b_ub->adjacency = connect->c2f;            /* shared pointer */
    }
    break;

  case CS_NAVSTO_SLES_DIAG_SCHUR_GCR:
  case CS_NAVSTO_SLES_DIAG_SCHUR_MINRES:
  case CS_NAVSTO_SLES_GCR:
  case CS_NAVSTO_SLES_GKB_SATURNE:
  case CS_NAVSTO_SLES_LOWER_SCHUR_GCR:
  case CS_NAVSTO_SLES_MINRES:
  case CS_NAVSTO_SLES_SGS_SCHUR_GCR:
  case CS_NAVSTO_SLES_UPPER_SCHUR_GCR:
  case CS_NAVSTO_SLES_UZAWA_AL:
  case CS_NAVSTO_SLES_UZAWA_CG:
  case CS_NAVSTO_SLES_UZAWA_SCHUR_GCR:
    {
      sc->assemble = _velocity_full_assembly;

      cs_lnum_t  block_sizes[2];
      block_sizes[0] = 3*quant->n_faces, block_sizes[1] = quant->n_cells;

      /* Create the system helper */

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       2,
                                       block_sizes,
                                       2);

      /* Add a first block for the (0,0) block and then define the underpinning
         structures */

      cs_cdo_system_matrix_class_t  matclass;

      /* Choose the right class of matrix to avoid copy.
       * The way to perform the assembly may change if an external librairy is
       * used for solving the linear system */

      switch (mom_eqp->sles_param->solver_class) {

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

      cs_cdo_system_block_t  *a =
        cs_cdo_system_add_dblock(sh, 0,                /* block id */
                                 matclass,
                                 cs_flag_primal_face , /* location */
                                 quant->n_faces,       /* n_elements */
                                 3,                    /* stride */
                                 true,                 /* interlaced */
                                 true);                /* unrolled */

      cs_cdo_system_build_block(sh, 0); /* build structures */

      /* Add a second block for the (1,0) and (0,1) blocks and then define the
         underpinning structures. The (0,1) block needs to be transposed before
         using it */

      cs_cdo_system_block_t  *bdiv =
        cs_cdo_system_add_ublock(sh, 1,               /* block_id */
                                 connect->c2f,        /* adjacency */
                                 cs_flag_primal_face, /* column location */
                                 quant->n_faces,      /* n_elements */
                                 3,                   /* stride */
                                 true);               /* interlaced */

      BFT_MALLOC(msles->div_op, 3*connect->c2f->idx[quant->n_cells], cs_real_t);

      cs_cdo_system_dblock_t  *a_db = a->block_pointer;
      cs_cdo_system_ublock_t  *b_ub = bdiv->block_pointer;

      /* Define the bdiv block by hand */

      b_ub->adjacency = connect->c2f;            /* shared pointer */
      b_ub->values = msles->div_op;              /* shared pointer */
      b_ub->shared_structures = true;
      b_ub->range_set = a_db->range_set;         /* shared pointer */
      b_ub->interface_set = a_db->interface_set; /* shared pointer */
    }
    break;

  case CS_NAVSTO_SLES_MUMPS:
    {
      if (nsp->model_flag & CS_NAVSTO_MODEL_WITH_SOLIDIFICATION)
        sc->assemble = _solidification_full_assembly;
      else
        sc->assemble = _full_assembly;

      cs_lnum_t  block_size = 3*quant->n_faces + quant->n_cells;

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       1,
                                       &block_size,
                                       1);

      /* Add the block and then define the underpinning structures */

      cs_cdo_system_block_t  *a =
        cs_cdo_system_add_xblock(sh, 0,        /* block id */
                                 block_size);  /* n_dofs */

      /* Define the xblock (with diagonal pressure block) */

      _build_shared_full_structures(a, true);
    }
    break;

  case CS_NAVSTO_SLES_NOTAY_TRANSFORM: /* Experimental */
    {
      sc->assemble = _notay_full_assembly;

      cs_lnum_t  block_size = 3*quant->n_faces + quant->n_cells;

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       1,
                                       &block_size,
                                       1);

      /* Add the block and then define the underpinning structures */

      cs_cdo_system_block_t  *a =
        cs_cdo_system_add_xblock(sh, 0,             /* block id */
                                 block_size);       /* n_dofs */

      /* Define the xblock (with diagonal pressure block) */

      _build_shared_full_structures(a, false);

      /* Define also the unassembled form of the divergence operator */

      BFT_MALLOC(msles->div_op, 3*connect->c2f->idx[quant->n_cells], cs_real_t);
    }
    break;

  default:
    /* CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK
     * CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG
     * CS_NAVSTO_SLES_DIAG_SCHUR_GMRES
     * CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK
     * CS_NAVSTO_SLES_GKB_PETSC
     * CS_NAVSTO_SLES_GKB_GMRES
     * CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK
     * CS_NAVSTO_SLES_NOTAY_TRANSFORM:
     * CS_NAVSTO_SLES_UPPER_SCHUR_GMRES
     */
    {
      sc->assemble = _full_assembly;

      cs_lnum_t block_size = 3*quant->n_faces + quant->n_cells;

      sh = cs_cdo_system_helper_create(CS_CDO_SYSTEM_SADDLE_POINT,
                                       1,
                                       &block_size,
                                       1);

      /* Add the block and then define the underpinning structures */

      cs_cdo_system_block_t  *a =
        cs_cdo_system_add_xblock(sh, 0,             /* block id */
                                 block_size);       /* n_dofs */

      /* Fill the xblock (without diagonal pressure block) */

      _build_shared_full_structures(a, false);
    }
    break;

  } /* Switch on strategy */

  sc->system_helper = sh;

  /* Set the grad-div scaling coefficient: Enable or not augmentation of the
     linear system */

  switch (strategy) {

  case CS_NAVSTO_SLES_BY_BLOCKS:
  case CS_NAVSTO_SLES_GKB_SATURNE:
  case CS_NAVSTO_SLES_MUMPS:
  case CS_NAVSTO_SLES_UZAWA_AL:
    msles->graddiv_coef = nsp->gd_scale_coef;
    break;

  default:
    /* CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK
     * CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG
     * CS_NAVSTO_SLES_DIAG_SCHUR_GCR
     * CS_NAVSTO_SLES_DIAG_SCHUR_GMRES
     * CS_NAVSTO_SLES_DIAG_SCHUR_MINRES
     * CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK
     * CS_NAVSTO_SLES_GCR
     * CS_NAVSTO_SLES_GKB_GMRES
     * CS_NAVSTO_SLES_GKB_PETSC
     * CS_NAVSTO_SLES_LOWER_SCHUR_GCR
     * CS_NAVSTO_SLES_MINRES
     * CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK
     * CS_NAVSTO_SLES_NOTAY_TRANSFORM
     * CS_NAVSTO_SLES_SGS_SCHUR_GCR
     * CS_NAVSTO_SLES_UPPER_SCHUR_GCR
     * CS_NAVSTO_SLES_UPPER_SCHUR_GMRES
     * CS_NAVSTO_SLES_UZAWA_CG
     * CS_NAVSTO_SLES_UZAWA_SCHUR_GCR
     */
    msles->graddiv_coef = 0;    /* No augmentation */
    break;

  } /* Switch on strategy */

  /* Set the solve function pointer */

  switch (nsp->sles_param->strategy) {

  case CS_NAVSTO_SLES_BY_BLOCKS:
    sc->solve = cs_cdofb_monolithic_krylov_block_precond;
    break;

  case CS_NAVSTO_SLES_GKB_SATURNE:
    sc->solve = cs_cdofb_monolithic_gkb_solve;
    break;

  case CS_NAVSTO_SLES_DIAG_SCHUR_GCR:
  case CS_NAVSTO_SLES_DIAG_SCHUR_MINRES:
  case CS_NAVSTO_SLES_GCR:
  case CS_NAVSTO_SLES_LOWER_SCHUR_GCR:
  case CS_NAVSTO_SLES_MINRES:
  case CS_NAVSTO_SLES_SGS_SCHUR_GCR:
  case CS_NAVSTO_SLES_UPPER_SCHUR_GCR:
  case CS_NAVSTO_SLES_UZAWA_SCHUR_GCR:
    sc->solve = cs_cdofb_monolithic_krylov_block_precond;
    break;

  case CS_NAVSTO_SLES_UZAWA_AL:
    sc->solve = cs_cdofb_monolithic_uzawa_al_incr_solve;
    break;

  case CS_NAVSTO_SLES_UZAWA_CG:
    sc->solve = cs_cdofb_monolithic_uzawa_cg_solve;
    break;

  default:
    /* CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK
     * CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG
     * CS_NAVSTO_SLES_DIAG_SCHUR_GMRES
     * CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK
     * CS_NAVSTO_SLES_GKB_PETSC
     * CS_NAVSTO_SLES_GKB_GMRES
     * CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK
     * CS_NAVSTO_SLES_MUMPS
     * CS_NAVSTO_SLES_NOTAY_TRANSFORM
     * CS_NAVSTO_SLES_UPPER_SCHUR_GMRES
     */
    sc->solve = cs_cdofb_monolithic_solve;
    break;

  }

  /* Set the pointer storing linear algebra features */

  sc->msles = msles;

  /* Iterative algorithm to handle the non-linearity (Picard by default) */

  const cs_navsto_param_sles_t  *nslesp = nsp->sles_param;

  sc->nl_algo = cs_iter_algo_create(nslesp->nl_algo_param);

  if (nslesp->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
    sc->nl_algo->context = cs_iter_algo_aa_create(nslesp->anderson_param,
                                                  quant->n_faces);

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

  /* Free the systemm helper */

  cs_cdo_system_helper_free(&(sc->system_helper));

  /* Free the context structure for solving saddle-point system */

  cs_cdofb_monolithic_sles_free(&(sc->msles));

  /* If the context is not NULL, this means that an Anderson algorithm has been
     activated otherwise nothing to do */

  cs_iter_algo_aa_free(sc->nl_algo);

  BFT_FREE(sc->nl_algo);

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
  cs_cdo_system_helper_t  *sh = sc->system_helper;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;

  /* Build an array storing the Dirichlet values at faces */

  cs_cdofb_vecteq_setup(t_cur, mesh, mom_eqp, mom_eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t  *rhs = NULL;  /* Since it is NULL, sh get sthe ownership */

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Main loop on cells to define the linear system to solve */

  sc->steady_build(nsp,
                   mom_eqc->face_values, sc->velocity->val,
                   NULL, NULL,  /* no value at time step n-1 */
                   sc);

  /* Free temporary buffers and structures */

  cs_equation_builder_reset(mom_eqb);

  /* End of the system building */

  cs_timer_t  t_bld_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_bld_end);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */

  _mono_fields_to_previous(sc, cc);

  /* Solve the linear system */

  cs_timer_t  t_solve_start = cs_timer_time();
  cs_cdofb_monolithic_sles_t  *msles = sc->msles;

  /* matrix has been already assigned to the msles structure */

  msles->sles = cs_sles_find_or_add(mom_eq->field_id, NULL);
  msles->u_f = mom_eqc->face_values; /* velocity DoFs at faces */
  msles->p_c = sc->pressure->val;    /* pressure DoFs at cells */

  int  cumulated_inner_iters = sc->solve(nsp, mom_eqp, sh,
                                         mom_eqp->sles_param, msles);

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

  /* Make sure that the DoFs are correctly enforced after the resolution */

  _mono_enforce_solid_face_velocity(mom_eqc->face_values);

  /* Now update the velocity and pressure fields associated to cells */

  _mono_update_related_cell_fields(nsp, sc, mom_eqc);

  /* Compute the new mass flux used as the advection field */

  cs_cdofb_navsto_mass_flux(nsp, quant, mom_eqc->face_values,
                            sc->mass_flux_array);

  if (nsp->verbosity > 1) {

    cs_log_printf(CS_LOG_DEFAULT, " -cvg- NavSto: cumulated_inner_iters: %d\n",
                  cumulated_inner_iters);
    cs_log_printf_flush(CS_LOG_DEFAULT);

  }

  /* Frees */

  cs_cdofb_monolithic_sles_clean(msles);

  cs_timer_t  t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Navier-Stokes system with a CDO face-based scheme
 *         using a monolithic approach and Picard iterations to solve the
 *         non-linearities arising from the advection term
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
  cs_cdo_system_helper_t  *sh = sc->system_helper;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;
  cs_iter_algo_t  *nl_algo = sc->nl_algo;
  cs_param_nl_algo_t  nl_algo_type = nsp->sles_param->nl_algo_type;

  /*--------------------------------------------------------------------------
   *                    INITIAL BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;

  /* Build an array storing the Dirichlet values at faces */

  cs_cdofb_vecteq_setup(t_cur, mesh, mom_eqp, mom_eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t  *rhs = NULL;  /* Since it is NULL, sh get sthe ownership */

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Main loop on cells to define the linear system to solve */

  sc->steady_build(nsp,
                   mom_eqc->face_values, sc->velocity->val,
                   NULL, NULL,  /* no value at time step n-1 */
                   sc);

  /* End of the system building */

  cs_timer_t  t_build_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_build_end);

  /*--------------------------------------------------------------------------
   *                   INITIAL BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */

  _mono_fields_to_previous(sc, cc);

  /* Solve the linear system */

  cs_timer_t  t_solve_start = cs_timer_time();

  cs_iter_algo_reset_nl(nl_algo_type, nl_algo);

  cs_cdofb_monolithic_sles_t  *msles = sc->msles;

  msles->sles = cs_sles_find_or_add(mom_eq->field_id, NULL);
  msles->u_f = mom_eqc->face_values; /* velocity DoFs at faces */
  msles->p_c = sc->pressure->val;    /* pressure DoFs at cells */

  /* Solve the new system:
   * Update the value of mom_eqc->face_values and sc->pressure->val */

  nl_algo->n_inner_iter =
    (nl_algo->last_inner_iter = sc->solve(nsp, mom_eqp, sh,
                                          mom_eqp->sles_param, msles));

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

  /* Make sure that the DoFs are correctly enforced after the resolution */

  _mono_enforce_solid_face_velocity(mom_eqc->face_values);

  /* Compute the new current mass flux used as the advection field */

  cs_cdofb_navsto_mass_flux(nsp, quant, mom_eqc->face_values,
                            sc->mass_flux_array);

  /* Set the normalization of the non-linear algo to the value of the first
     mass flux norm */

  nl_algo->normalization =
    sqrt(cs_cdo_blas_square_norm_pfsf(sc->mass_flux_array));

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: START
   *--------------------------------------------------------------------------*/

  /* Check the convergence status and update the nl_algo structure related
   * to the convergence monitoring */

  while (cs_cdofb_navsto_nl_algo_cvg(nl_algo_type,
                                     sc->mass_flux_array_pre,
                                     sc->mass_flux_array,
                                     nl_algo) == CS_SLES_ITERATING) {

    /* Main loop on cells to define the linear system to solve */

    cs_timer_t  t_build_start = cs_timer_time();

    /* rhs set to zero and cs_sles_t structure is freed in order to do
     * the setup once again since the matrix should be modified */

    cs_cdo_system_helper_init_system(sh, &rhs);
    cs_cdofb_monolithic_sles_clean(msles);

    sc->steady_build(nsp,
                     /* A current to previous op. has been done */
                     mom_eqc->face_values_pre, sc->velocity->val_pre,
                     NULL, NULL,  /* no value at time step n-1 */
                     sc);

    /* End of the system building */

    t_build_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_build_start, &t_build_end);

    /* Solve the new system:
     * Update the value of mom_eqc->face_values and sc->pressure->val */

    t_solve_start = cs_timer_time();

    nl_algo->n_inner_iter +=
      (nl_algo->last_inner_iter = sc->solve(nsp, mom_eqp, sh,
                                            mom_eqp->sles_param, msles));

    t_solve_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

    /* Make sure that the DoFs are correctly enforced after the resolution */

    _mono_enforce_solid_face_velocity(mom_eqc->face_values);

    /* Compute the new mass flux used as the advection field */

    memcpy(sc->mass_flux_array_pre, sc->mass_flux_array,
           n_faces*sizeof(cs_real_t));

    cs_cdofb_navsto_mass_flux(nsp, quant, mom_eqc->face_values,
                              sc->mass_flux_array);

  } /* Loop on non-linear iterations */

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: END
   *--------------------------------------------------------------------------*/

  if (nsp->verbosity > 1) {

    cs_log_printf(CS_LOG_DEFAULT, " -cvg- NavSto: cumulated_inner_iters: %d\n",
                  nl_algo->n_inner_iter);
    cs_log_printf_flush(CS_LOG_DEFAULT);

  }

  cs_iter_algo_post_check(__func__,
                          mom_eqp->name,
                          cs_param_get_nl_algo_label(nl_algo_type),
                          nl_algo);

  if (nsp->sles_param->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
    cs_iter_algo_aa_free_arrays(nl_algo->context);

  /* Now compute/update the velocity and pressure fields */

  _mono_update_related_cell_fields(nsp, sc, mom_eqc);

  /* Frees */

  cs_equation_builder_reset(mom_eqb);
  cs_cdofb_monolithic_sles_clean(msles);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */

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
cs_cdofb_monolithic(const cs_mesh_t           *mesh,
                    const cs_navsto_param_t   *nsp,
                    void                      *scheme_context)
{
  const cs_timer_t  t_start = cs_timer_time();

  /* Retrieve high-level structures */

  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t *)scheme_context;
  cs_cdo_system_helper_t  *sh = sc->system_helper;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;

  /*--------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_eval = ts->t_cur + ts->dt[0];

  /* Build an array storing the Dirichlet values at faces */

  cs_cdofb_vecteq_setup(t_eval, mesh, mom_eqp, mom_eqb);

  /* Initialize the matrix and all its related structures needed during
   * the assembly step as well as the rhs */

  cs_real_t  *rhs = NULL;  /* Since it is NULL, sh get sthe ownership */

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Main loop on cells to define the linear system to solve */

  sc->build(nsp,
            mom_eqc->face_values, sc->velocity->val,
            mom_eqc->face_values_pre, sc->velocity->val_pre,
            sc);

  /* Free temporary buffers and structures */

  cs_equation_builder_reset(mom_eqb);

  /* End of the system building */

  cs_timer_t  t_bld_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_bld_end);

  /*--------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */

  _mono_fields_to_previous(sc, cc);

  /* Solve the linear system */

  cs_timer_t  t_solve_start = cs_timer_time();
  cs_cdofb_monolithic_sles_t  *msles = sc->msles;

  /* matrix has been already assigned to the msles structure */

  msles->sles = cs_sles_find_or_add(mom_eq->field_id, NULL);
  msles->u_f = mom_eqc->face_values; /* velocity DoFs at faces */
  msles->p_c = sc->pressure->val;    /* pressure DoFs at cells */

  int  cumulated_inner_iters = sc->solve(nsp, mom_eqp, sh,
                                         mom_eqp->sles_param, msles);

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

  /* Make sure that the DoFs are correctly enforced after the resolution */

  _mono_enforce_solid_face_velocity(mom_eqc->face_values);

  /* Now update the velocity and pressure fields associated to cells */

  _mono_update_related_cell_fields(nsp, sc, mom_eqc);

  /* Compute the new mass flux used as the advection field */

  cs_cdofb_navsto_mass_flux(nsp, quant, mom_eqc->face_values,
                            sc->mass_flux_array);

  if (nsp->verbosity > 1) {

    cs_log_printf(CS_LOG_DEFAULT, " -cvg- NavSto: cumulated_inner_iters: %d\n",
                  cumulated_inner_iters);
    cs_log_printf_flush(CS_LOG_DEFAULT);

  }

  /* Frees */

  cs_cdofb_monolithic_sles_clean(msles);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */

  cs_timer_t  t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the unsteady Navier-Stokes system with a CDO face-based scheme
 *         using a monolithic approach.
 *         According to the settings, this function can handle either an
 *         implicit Euler time scheme or more generally a theta time scheme.
 *         Rely on Picard iterations to solve the non-linearities arising from
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

  cs_cdofb_monolithic_t  *sc = scheme_context;
  cs_cdo_system_helper_t  *sh = sc->system_helper;
  cs_navsto_monolithic_t *cc = sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc = mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;
  cs_iter_algo_t  *nl_algo = sc->nl_algo;
  cs_param_nl_algo_t  nl_algo_type = nsp->sles_param->nl_algo_type;

  /*--------------------------------------------------------------------------
   *                    INITIAL BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_eval = ts->t_cur + ts->dt[0];

  /* Build an array storing the Dirichlet values at faces */

  cs_cdofb_vecteq_setup(t_eval, mesh, mom_eqp, mom_eqb);

  /* Initialize the rhs */

  cs_real_t  *rhs = NULL;

  cs_cdo_system_helper_init_system(sh, &rhs);

  /* Main loop on cells to define the linear system to solve */

  sc->build(nsp,
            mom_eqc->face_values, sc->velocity->val,
            mom_eqc->face_values_pre, sc->velocity->val_pre,
            sc);

  /* End of the system building */

  cs_timer_t  t_build_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_start, &t_build_end);

  /*--------------------------------------------------------------------------
   *                   INITIAL BUILD: END
   *--------------------------------------------------------------------------*/

  /* Current to previous for main variable fields */

  _mono_fields_to_previous(sc, cc);

  /* Solve the linear system */

  cs_timer_t  t_solve_start = cs_timer_time();
  cs_cdofb_monolithic_sles_t  *msles = sc->msles;

  msles->sles = cs_sles_find_or_add(mom_eq->field_id, NULL);
  msles->u_f = mom_eqc->face_values; /* velocity DoFs at faces */
  msles->p_c = sc->pressure->val;    /* pressure DoFs at cells */

  cs_iter_algo_reset_nl(nl_algo_type, nl_algo);

  /* Solve the new system:
   * Update the value of mom_eqc->face_values and sc->pressure->val */

  nl_algo->n_inner_iter =
    (nl_algo->last_inner_iter = sc->solve(nsp, mom_eqp, sh,
                                          mom_eqp->sles_param, msles));

  cs_timer_t  t_solve_end = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

  /* Make sure that the DoFs are correctly enforced after the resolution */

  _mono_enforce_solid_face_velocity(mom_eqc->face_values);

  /* Compute the new mass flux used as the advection field */

  cs_cdofb_navsto_mass_flux(nsp, quant, mom_eqc->face_values,
                            sc->mass_flux_array);

  /* Set the normalization of the non-linear algo to the value of the first
     mass flux norm */

  nl_algo->normalization =
    sqrt(cs_cdo_blas_square_norm_pfsf(sc->mass_flux_array));

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: START
   *--------------------------------------------------------------------------*/

  /* Since a current to previous op. has been done:
   *   sc->mass_flux_array_pre -> flux at t^n= t^n,0 (= t^(n-1)
   *   sc->mass_flux_array     -> flux at t^n,1 (call to .._navsto_mass_flux */

  cs_cdofb_navsto_nl_algo_cvg(nl_algo_type,
                              sc->mass_flux_array_pre,
                              sc->mass_flux_array,
                              nl_algo);

  cs_real_t  *mass_flux_array_k = NULL;
  cs_real_t  *mass_flux_array_kp1 = sc->mass_flux_array;

  while (nl_algo->cvg == CS_SLES_ITERATING) {

    /* Start of the system building */

    cs_timer_t  t_build_start = cs_timer_time();

    /* rhs set to zero and cs_sles_t structure is freed in order to do
     * the setup once again since the matrix should be modified */

    cs_cdo_system_helper_init_system(sh, &rhs);
    cs_cdofb_monolithic_sles_clean(msles);

    /* Main loop on cells to define the linear system to solve */

    sc->build(nsp,
              /* A current to previous op. has been done */
              mom_eqc->face_values_pre, sc->velocity->val_pre,
              NULL, NULL, /* no n-1 state is given */
              sc);

    /* End of the system building */

    t_build_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_build_start, &t_build_end);

    /* Solve the new system */

    t_solve_start = cs_timer_time();

    nl_algo->n_inner_iter +=
      (nl_algo->last_inner_iter = sc->solve(nsp, mom_eqp, sh,
                                            mom_eqp->sles_param, msles));

    t_solve_end = cs_timer_time();
    cs_timer_counter_add_diff(&(mom_eqb->tcs), &t_solve_start, &t_solve_end);

    /* Make sure that the DoFs are correctly enforced after the resolution */

    _mono_enforce_solid_face_velocity(mom_eqc->face_values);

    /* mass_flux_array_k <-- mass_flux_array_kp1; update mass_flux_array_kp1 */

    if (mass_flux_array_k == NULL)
      BFT_MALLOC(mass_flux_array_k, n_faces, cs_real_t);
    memcpy(mass_flux_array_k, mass_flux_array_kp1, n_faces*sizeof(cs_real_t));

    cs_cdofb_navsto_mass_flux(nsp, quant, mom_eqc->face_values,
                              mass_flux_array_kp1);

    /* Check the convergence status and update the nl_algo structure related
     * to the convergence monitoring */

    cs_cdofb_navsto_nl_algo_cvg(nl_algo_type,
                                mass_flux_array_k,
                                mass_flux_array_kp1,
                                nl_algo);

  } /* Loop on non-linear iterations */

  /*--------------------------------------------------------------------------
   *                   PICARD ITERATIONS: END
   *--------------------------------------------------------------------------*/

  if (nsp->verbosity > 1) {

    cs_log_printf(CS_LOG_DEFAULT, " -cvg- NavSto: cumulated_inner_iters: %d\n",
                  nl_algo->n_inner_iter);
    cs_log_printf_flush(CS_LOG_DEFAULT);

  }

  cs_iter_algo_post_check(__func__,
                          mom_eqp->name,
                          cs_param_get_nl_algo_label(nl_algo_type),
                          nl_algo);

  if (nsp->sles_param->nl_algo_type == CS_PARAM_NL_ALGO_ANDERSON)
    cs_iter_algo_aa_free_arrays(nl_algo->context);

  /* Now compute/update the velocity and pressure fields */

  _mono_update_related_cell_fields(nsp, sc, mom_eqc);

  /* Frees */

  cs_cdofb_monolithic_sles_clean(msles);
  cs_cdo_system_helper_reset(sh);      /* free rhs and matrix */
  cs_equation_builder_reset(mom_eqb);
  if (mass_flux_array_k != NULL)
    BFT_FREE(mass_flux_array_k);

  cs_timer_t  t_end = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_start, &t_end);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
