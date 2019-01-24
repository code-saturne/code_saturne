/*============================================================================
 * Build an algebraic CDO face-based system for the Navier-Stokes equations
 * and solved it as one block (monolithic approach)
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
#include "cs_dbg.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_navsto_coupling.h"
#include "cs_param.h"
#include "cs_post.h"
#include "cs_sles.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

#if defined(HAVE_PETSC)
#include "cs_sles_petsc.h"
#endif

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
   * @name Performance monitoring
   * Monitoring the efficiency of the algorithm used to solve the Navier-Stokes
   * system
   * @{
   *
   *  \var bf_type
   *  Array of boundary type for each boundary face. (Shared)
   */

  const cs_boundary_type_t       *bf_type;

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

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOFB_MONOLITHIC_DBG      0

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

#if defined(HAVE_PETSC)
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generate IndexSet for the PETSc FieldSplit preconditioner
 *
 * \param[in, out]  isp     IndexSet for the pressure DoFs
 * \param[in, out]  isv     IndexSet for the velocity DoFs
 */
/*----------------------------------------------------------------------------*/

static void
_build_is_for_fieldsplit(IS   *isp,
                         IS   *isv)
{
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_range_set_t  *rset = cs_shared_range_set;

  PetscInt  n_faces = quant->n_faces;
  PetscInt  n_cells = quant->n_cells;
  PetscInt  *indices = NULL;

  PetscMalloc1(3*n_faces, &indices);

  /* IndexSet for the velocity DoFs */
  if (rset->n_elts[0] == rset->n_elts[1]) {

    for (PetscInt i = 0; i < 3*n_faces; i++)
      indices[i] = rset->g_id[i];
    ISCreateGeneral(PETSC_COMM_WORLD, 3*n_faces, indices, PETSC_COPY_VALUES,
                    isv);

  }
  else {

    PetscInt  n_velocity_elts = 0;
    for (PetscInt i = 0; i < 3*n_faces; i++) {
      cs_gnum_t  g_id = rset->g_id[i];
      if (g_id >= rset->l_range[0] && g_id < rset->l_range[1])
        indices[n_velocity_elts++] = g_id;
    }
    ISCreateGeneral(PETSC_COMM_WORLD, n_velocity_elts, indices,
                    PETSC_COPY_VALUES, isv);

  }

  /* IndexSet for the velocity DoFs
   * Pressure unknows are located at cell centers so the treatment should be
   * the same in sequential and parallel computation*/
  for (PetscInt i = 0; i < n_cells; i++)
    indices[i] = rset->g_id[i + 3*n_faces];
  ISCreateGeneral(PETSC_COMM_WORLD, n_cells, indices, PETSC_COPY_VALUES,
                  isp);

  PetscFree(indices);
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of additive block preconditioner for a GMRES
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_additive_amg_gmres_hook(void     *context,
                         Mat       a,
                         KSP       ksp)
{
  IS  isv, isp;

  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  const int  n_max_restart = 30;

  KSPSetType(ksp, KSPFGMRES);
  KSPGMRESSetRestart(ksp, n_max_restart);

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp.eps,         /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  /* Try to have "true" norm */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

#if 0 /* JB: TEST TO PERFORM IN 3D*/
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsSetValue(NULL,
                  "-fieldsplit_velocity_pc_hypre_boomeramg_strong_threshold",
                  "0.5");
#else
  PetscOptionsSetValue(
                  "-fieldsplit_velocity_pc_hypre_boomeramg_strong_threshold",
                  "0.5");
#endif
#endif

  /* Apply modifications to the KSP structure */
  PC up_pc, u_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_ADDITIVE);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPPREONLY);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCJACOBI);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  KSP  u_ksp = up_subksp[0];
  KSPSetType(u_ksp, KSPCG);
  KSPGetPC(u_ksp, &u_pc);
  PCSetType(u_pc, PCHYPRE);
  PCHYPRESetType(u_pc, "boomeramg");
  KSPSetTolerances(u_ksp,
                   slesp.eps,   /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   5);          /* max number of iterations */

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);

  KSPSetFromOptions(u_ksp);
  KSPSetUp(u_ksp);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of diagonal Schur preconditioner by block for a GMRES
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_diag_schur_gmres_hook(void     *context,
                       Mat       a,
                       KSP       ksp)
{
  IS  isv, isp;

  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  const int  n_max_restart = 30;

  KSPSetType(ksp, KSPFGMRES);
  KSPGMRESSetRestart(ksp, n_max_restart);

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp.eps,         /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  /* Try to have "true" norm */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

  /* Apply modifications to the KSP structure */
  PC up_pc, u_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_SCHUR);
  PCFieldSplitSetSchurFactType(up_pc, PC_FIELDSPLIT_SCHUR_FACT_DIAG);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPMINRES);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCNONE);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  KSP  u_ksp = up_subksp[0];
  KSPSetType(u_ksp, KSPCG);
  KSPGetPC(u_ksp, &u_pc);
  PCSetType(u_pc, PCHYPRE);
  PCHYPRESetType(u_pc, "boomeramg");
  KSPSetTolerances(u_ksp,
                   slesp.eps,   /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   5);          /* max number of iterations */

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);
  KSPSetUp(u_ksp);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);
}

/*----------------------------------------------------------------------------
 * \brief  Function pointer: setup hook for setting PETSc solver and
 *         preconditioner.
 *         Case of upper Schur preconditioner by block for a GMRES
 *
 * \param[in, out] context  pointer to optional (untyped) value or structure
 * \param[in, out] a        pointer to PETSc Matrix context
 * \param[in, out] ksp      pointer to PETSc KSP context
 *----------------------------------------------------------------------------*/

static void
_upper_schur_gmres_hook(void     *context,
                        Mat       a,
                        KSP       ksp)
{
  IS  isv, isp;

  cs_equation_param_t  *eqp = (cs_equation_param_t *)context;
  cs_param_sles_t  slesp = eqp->sles_param;

  const int  n_max_restart = 30;

  KSPSetType(ksp, KSPFGMRES);
  KSPGMRESSetRestart(ksp, n_max_restart);

  /* Set KSP tolerances */
  PetscReal rtol, abstol, dtol;
  PetscInt  maxit;
  KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxit);
  KSPSetTolerances(ksp,
                   slesp.eps,         /* relative convergence tolerance */
                   abstol,            /* absolute convergence tolerance */
                   dtol,              /* divergence tolerance */
                   slesp.n_max_iter); /* max number of iterations */

  /* Try to have "true" norm */
  KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);

  /* Apply modifications to the KSP structure */
  PC up_pc, u_pc, p_pc;

  KSPGetPC(ksp, &up_pc);
  PCSetType(up_pc, PCFIELDSPLIT);
  PCFieldSplitSetType(up_pc, PC_COMPOSITE_SCHUR);
  PCFieldSplitSetSchurFactType(up_pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER);

  _build_is_for_fieldsplit(&isp, &isv);

  /* First level Pressure | Velocity (X,Y,Z) */
  PCFieldSplitSetIS(up_pc, "velocity", isv);
  PCFieldSplitSetIS(up_pc, "pressure", isp);

  /* Need to call PCSetUp before configuring the second level (Thanks to
     Natacha Bereux) */
  PCSetFromOptions(up_pc);
  PCSetUp(up_pc);
  KSPSetUp(ksp);

  PetscInt  n_split;
  KSP  *up_subksp;
  PCFieldSplitGetSubKSP(up_pc, &n_split, &up_subksp);
  assert(n_split == 2);

  KSP  p_ksp = up_subksp[1];
  KSPSetType(p_ksp, KSPMINRES);
  KSPGetPC(p_ksp, &p_pc);
  PCSetType(p_pc, PCNONE);

  PCSetFromOptions(p_pc);
  PCSetUp(p_pc);
  KSPSetUp(p_ksp);

  KSP  u_ksp = up_subksp[0];
  KSPSetType(u_ksp, KSPCG);
  KSPGetPC(u_ksp, &u_pc);
  PCSetType(u_pc, PCHYPRE);
  PCHYPRESetType(u_pc, "boomeramg");
  KSPSetTolerances(u_ksp,
                   slesp.eps,   /* relative convergence tolerance */
                   abstol,      /* absolute convergence tolerance */
                   dtol,        /* divergence tolerance */
                   5);          /* max number of iterations */

  PCSetFromOptions(u_pc);
  PCSetUp(u_pc);
  KSPSetUp(u_ksp);

  /* User function for additional settings */
  cs_user_sles_petsc_hook(context, a, ksp);

  /* Apply modifications to the KSP structure */
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  /* Dump the setup related to PETSc in a specific file */
  if (!slesp.setup_done) {
    cs_sles_petsc_log_setup(ksp);
    slesp.setup_done = true;
  }

  PetscFree(up_subksp);
  ISDestroy(&isp);
  ISDestroy(&isv);
}
#endif  /* HAVE_PETSC */

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
  const cs_adjacency_t  *f2f = cs_equation_get_f2f_index();
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
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      bf_type     type of boundary for the boundary face
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_apply_bc_partly(const cs_cdofb_monolithic_t   *sc,
                 const cs_equation_param_t     *eqp,
                 const cs_cdofb_vecteq_t       *eqc,
                 const cs_cell_mesh_t          *cm,
                 const cs_boundary_type_t      *bf_type,
                 cs_cell_sys_t                 *csys,
                 cs_cell_builder_t             *cb)
{
  CS_UNUSED(eqc);
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

      switch (bf_type[i]) {

      case CS_BOUNDARY_INLET:
        if (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          sc->apply_velocity_inlet(f, eqp, cm, cb, csys);
        }
        break;

      case CS_BOUNDARY_SLIDING_WALL:
        if (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          sc->apply_sliding_wall(f, eqp, cm, cb, csys);
        }
        break;

      case CS_BOUNDARY_WALL:
        if (eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_NITSCHE ||
            eqp->enforcement == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
          sc->apply_fixed_wall(f, eqp, cm, cb, csys);
        }
        break;

      case CS_BOUNDARY_SYMMETRY:
        /* Always weakly enforce the symmetric constraint on the
           velocity-block */
        sc->apply_symmetry(f, eqp, cm, cb, csys);
        break;

      default: /* Nothing to do */
        /* Remark: Case of a "natural" outlet */
        break;

      } /* End of switch */

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
 * \param[in]      bf_type     type of boundary for the boundary face
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 * \param[in, out] div_op      array with the divergence op. values
 * \param[in, out] mass_rhs    pointer to the rhs for the mass eq. in this cell
 */
/*----------------------------------------------------------------------------*/

static void
_apply_remaining_bc(const cs_cdofb_monolithic_t   *sc,
                    const cs_equation_param_t     *eqp,
                    const cs_cell_mesh_t          *cm,
                    const cs_boundary_type_t      *bf_type,
                    cs_cell_sys_t                 *csys,
                    cs_cell_builder_t             *cb,
                    cs_real_t                     *div_op,
                    cs_real_t                     *mass_rhs)
{
  if (csys->cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Update the divergence operator and the right-hand side related to the
     * mass equation.
     * Enforcement of Dirichlet BC in a stronger way if this is the choice
     */
    for (short int i = 0; i < csys->n_bc_faces; i++) {

      /* Get the boundary face in the cell numbering */
      const short int  f = csys->_f_ids[i];

      switch (bf_type[i]) {

      case CS_BOUNDARY_INLET:
        /* Update mass RHS from the knowledge of the mass flux
         * TODO: multiply by rho */
        mass_rhs[0] -= cs_math_3_dot_product(csys->dir_values + 3*f,
                                             div_op + 3*f);

        /* Strong enforcement of u.n (--> dp/dn = 0) on the divergence */
        for (int k = 0; k < 3; k++) div_op[3*f+k] = 0;

        /* Enforcement of the velocity for the velocity-block
         * Dirichlet on the three components of the velocity field */
        if (eqp->enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
            eqp->enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {
          sc->apply_velocity_inlet(f, eqp, cm, cb, csys);
        }
        break;

      case CS_BOUNDARY_SLIDING_WALL:
        /* No need to update the mass RHS since there is no mass flux */

        /* Strong enforcement of u.n (--> dp/dn = 0) on the divergence */
        for (int k = 0; k < 3; k++) div_op[3*f+k] = 0;

        /* Enforcement of the velocity for the velocity-block
         * Dirichlet on the three components of the velocity field */
        if (eqp->enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
            eqp->enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {
          sc->apply_sliding_wall(f, eqp, cm, cb, csys);
        }
        break;

      case CS_BOUNDARY_WALL:
        /* No need to update the mass RHS since there is no mass flux */

        /* Strong enforcement of u.n (--> dp/dn = 0) on the divergence */
        for (int k = 0; k < 3; k++) div_op[3*f+k] = 0;

        /* Enforcement of the velocity for the velocity-block
         * Dirichlet on the three components of the velocity field */
        if (eqp->enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
            eqp->enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC) {
          sc->apply_fixed_wall(f, eqp, cm, cb, csys);
        }
        break;

      case CS_BOUNDARY_SYMMETRY:

        /* No need to update the mass RHS since there is no mass flux */

        /* Weak-enforcement for the velocity-block (cf. _apply_bc_partly) */

        /* Strong enforcement of u.n (--> dp/dn = 0) on the divergence */
        for (int k = 0; k < 3; k++) div_op[3*f+k] = 0;

        break;

      default: /* Nothing to do */
        /* Remark: Case of a "natural" outlet */
        break;

      } /* End of switch */

    } /* Loop boundary faces */

  } /* Boundary cell */
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
   * ==================
   */

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

  /* 2. RHS assembly
   * ===============
   */

  for (short int f = 0; f < 3*n_f; f++)
#   pragma omp atomic
    rhs[csys->dof_ids[f]] += csys->rhs[f];

  /* Reset the value of the source term for the cell DoF
     Source term is only hold by the cell DoF in face-based schemes */
  if (has_sourceterm)
    for (int k = 0; k < 3; k++)
      eqc_st[3*cm->c_id + k] = csys->source[3*n_f + k];
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
 */
/*----------------------------------------------------------------------------*/

static void
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

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_VECTEQ_DBG > 1
  if (cs_glob_n_ranks > 1) /* Parallel mode */
    cs_range_set_scatter(rset,
                         CS_REAL_TYPE, 1, /* type and stride */
                         b, b);

  cs_dbg_fprintf_system(eqp->name, cs_shared_time_step->nt_cur,
                        CS_CDOFB_VECTEQ_DBG,
                        xsol, b, 3*n_faces);
#endif

  /* Interlace xsol --> vel_f and pre_c */
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
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the linear system and update the velocity fields and the
 *         pressure field
 *
 * \param[in]      sc       scheme context
 * \param[in]      mom_eq   equation structure related to the momentum
 * \param[in]      matrix   pointer to a cs_matrix_t structure
 * \param[in, out] mom_rhs  right-hand side for the momentum equation
 * \param[in, out] mass_rhs right_hand side for the mass equation
 */
/*----------------------------------------------------------------------------*/

static void
_compute_and_update_fields(const cs_matrix_t             *matrix,
                           cs_cdofb_monolithic_t         *sc,
                           cs_equation_t                 *mom_eq,
                           cs_real_t                     *mom_rhs,
                           cs_real_t                     *mass_rhs)
{
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  /* Copy current field values to previous values */
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t  *mom_eqp = mom_eq->param;
  cs_equation_builder_t  *mom_eqb = mom_eq->builder;

  cs_field_t  *pr_fld = sc->pressure;
  cs_field_t  *vel_fld = sc->velocity;

  cs_timer_t  t_upd = cs_timer_time();

  cs_field_current_to_previous(vel_fld);
  cs_field_current_to_previous(pr_fld);
  cs_field_current_to_previous(sc->divergence);

  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

  /* Now solve the system */
  cs_real_t  *vel_f = mom_eqc->face_values;
  cs_sles_t  *sles = cs_sles_find_or_add(mom_eq->field_id, NULL);

  _solve_system(sles, matrix, mom_eqp, vel_f, pr_fld->val, mom_rhs, mass_rhs);

  /* Update pressure, velocity and divergence fields */
  t_upd = cs_timer_time();

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

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tce), &t_upd, &t_tmp);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
  cs_dbg_darray_to_listing("VELOCITY", quant->n_faces, vel_f, 9);
  cs_dbg_darray_to_listing("PRESSURE", quant->n_cells, pr_fld->val, 9);
  cs_dbg_darray_to_listing("VELOCITY_DIV", quant->n_cells, div, 9);
#endif

  /* Frees */
  cs_sles_free(sles);

}

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
  cs_equation_param_t *mom_eqp = mom_eq->param;

  sc->coupling_context = cc; /* shared with cs_navsto_system_t */

  /* Quick access to the main fields */
  sc->velocity = cs_field_by_name("velocity");
  sc->pressure = cs_field_by_name("pressure");
  sc->divergence = cs_field_by_name("velocity_divergence");

  /* Boundary treatment */
  sc->bf_type = bf_type;

  /* Set the way to enforce the Dirichlet BC on the velocity
   * "fixed_wall" means a no-slip BC
   */
  sc->apply_symmetry = cs_cdofb_symmetry;

  switch (mom_eqp->enforcement) {

  case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_alge;
    sc->apply_sliding_wall = cs_cdofb_block_dirichlet_alge;
    sc->apply_fixed_wall = cs_cdofb_block_dirichlet_alge;
    break;

  case CS_PARAM_BC_ENFORCE_PENALIZED:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_pena;
    sc->apply_sliding_wall = cs_cdofb_block_dirichlet_pena;
    sc->apply_fixed_wall = cs_cdofb_block_dirichlet_pena;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_weak;
    sc->apply_sliding_wall = cs_cdofb_block_dirichlet_weak;
    sc->apply_fixed_wall = cs_cdofb_block_dirichlet_weak;
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    sc->apply_velocity_inlet = cs_cdofb_block_dirichlet_wsym;
    sc->apply_sliding_wall = cs_cdofb_block_dirichlet_wsym;
    sc->apply_fixed_wall = cs_cdofb_block_dirichlet_wsym;
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

  /* Other pointers are only shared (i.e. not owner) */
  BFT_FREE(sc);

  if (cs_shared_range_set != NULL) { /* Shared structures have to be freed */
    cs_range_set_destroy(&cs_shared_range_set);
    cs_interface_set_destroy(&cs_shared_interface_set);
    cs_matrix_structure_destroy(&cs_shared_matrix_structure);
    cs_matrix_assembler_destroy(&cs_shared_matrix_assembler);
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Start setting-up the Navier-Stokes equations when a monolithic
 *         algorithm is used to couple the system.
 *         No mesh information is available at this stage
 *
 * \param[in]      nsp      pointer to a \ref cs_navsto_param_t structure
 * \param[in, out] context  pointer to a context structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_monolithic_set_sles(const cs_navsto_param_t    *nsp,
                             void                       *context)
{
  cs_navsto_monolithic_t  *nsc = (cs_navsto_monolithic_t *)context;

  assert(nsp != NULL && nsc != NULL);

  cs_equation_param_t  *mom_eqp = cs_equation_get_param(nsc->momentum);
  int  field_id = cs_equation_get_field_id(nsc->momentum);

  /* Initialization must be called before setting options;
     it does not need to be called before calling
     cs_sles_petsc_define(), as this is handled automatically. */

  switch (nsp->sles_strategy) {

  case CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK: /* "Classical" way to set SLES */
    cs_equation_param_set_sles(mom_eqp, field_id);
    break;

#if defined(HAVE_PETSC)
  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _additive_amg_gmres_hook,
                         (void *)mom_eqp);
    break;

  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _diag_schur_gmres_hook,
                         (void *)mom_eqp);
    break;

  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
    cs_sles_petsc_init();
    cs_sles_petsc_define(field_id,
                         NULL,
                         MATMPIAIJ,
                         _upper_schur_gmres_hook,
                         (void *)mom_eqp);
    break;
#else
  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n"
              " PETSc is required with this option.\n"
              " Please build a version of Code_Saturne with the PETSc support.",
              __func__, mom_eqp->name);
    break;
#endif /* HAVE_PETSC */

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid strategy for solving the linear system %s\n",
              __func__, mom_eqp->name);
  }

      /* Define the level of verbosity for SLES structure */
  if (mom_eqp->sles_param.verbosity > 1) {

    cs_sles_t  *sles = cs_sles_find_or_add(field_id, NULL);

    /* Set verbosity */
    cs_sles_set_verbosity(sles, mom_eqp->sles_param.verbosity);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Navier-Stokes system with a CDO face-based scheme
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

  cs_timer_t  t_cmpt = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  /* Retrieve high-level structures */
  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t *)scheme_context;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;
  cs_real_t  *vel_c = sc->velocity->val;

  /*----------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur, time_eval = t_cur; /* Dummy arguments */

  cs_timer_t  t_bld = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces. */
  cs_real_t  *dir_values = NULL;
  cs_cdofb_vecteq_setup_bc(t_cur, mesh, mom_eqp, mom_eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_matrix_structure);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if  (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  cs_real_t  *mass_rhs = NULL;
  BFT_MALLOC(mass_rhs, quant->n_cells, cs_real_t);
  /* Initialization done after */

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, mom_eqp, mom_eqb, mom_eqc, rhs, matrix,        \
         mass_rhs, mav, dir_values, vel_c, sc)
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
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;
    cs_real_t  *div_op = NULL;
    cs_boundary_type_t  *bf_type = NULL;

    BFT_MALLOC(div_op, 3*connect->n_max_fbyc, cs_real_t);
    BFT_MALLOC(bf_type, connect->n_max_fbyc, cs_boundary_type_t);
    cs_cdofb_vecteq_get(&csys, &cb);

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, time_eval, cb);

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
                                       dir_values, vel_c, time_eval,
                                       csys, cb);

      const short int n_fc = cm->n_fc;

      /* Define the local bf_type */
      for (short int i = 0; i < csys->n_bc_faces; i++) {

        /* Get the boundary face in the cell numbering */
        const short int  f = csys->_f_ids[i];
        const cs_lnum_t  bf_id = cm->f_ids[f] - csys->face_shift;

        bf_type[i] = sc->bf_type[bf_id];

      }

      /* Initialize the RHS for the mass eq */
      mass_rhs[c_id] = 0.;

      /* 1- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_diffusion(time_eval, mom_eqp, mom_eqb, mom_eqc,
                                cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion", csys);
#endif

      /* 2- PRESSURE (SCALAR) EQUATION */
      /* ============================= */

      cs_cdofb_navsto_divergence_vect(cm, cb->aux->val);

      cs_real_t *_div = cb->aux->val;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
      if (cs_dbg_cw_test(mom_eqp, cm, csys)) {
        const cs_real_t  ovc = 1. / cm->vol_c;
#       pragma omp critical
        {
          cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
          for (short int f = 0; f < n_fc; f++)
            cs_log_printf(CS_LOG_DEFAULT, "    f%2d: %- .4e, %- .4e, %- .4e\n",
                          f, _div[3*f]*ovc, _div[3*f+1]*ovc, _div[3*f+2]*ovc);
        } /* Critical section */
      }
#endif

      /* -divergence is used in the system */
      for (int i = 0; i < 3*n_fc; i++)
        div_op[i] = -_div[i];

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const _Bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm) {

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   time_eval, 1., /* time (dummy), scaling */
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* 3b- OTHER RHS CONTRIBUTIONS *
       * =========================== *
       *
       * First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _apply_bc_partly(sc, mom_eqp, mom_eqc, cm, bf_type, csys, cb);

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

      _apply_remaining_bc(sc, mom_eqp, cm, bf_type, csys, cb, div_op,
                          mass_rhs + c_id);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      _assemble(csys, cm, div_op,
                has_sourceterm, mav, rhs, mom_eqc->source_terms);

    } /* Main loop on cells */

    /* Free temporary buffer */
    BFT_FREE(div_op);
    BFT_FREE(bf_type);

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*----------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Now compute/update the velocity and pressure fields */
  _compute_and_update_fields(matrix, sc, mom_eq, rhs, mass_rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
  cs_dbg_fprintf_system(mom_eqp->name, ts->nt_cur, CS_CDOFB_MONOLITHIC_DBG,
                        vel_f, rhs, 3*n_faces);
#endif

  /* Frees */
  BFT_FREE(rhs);
  BFT_FREE(mass_rhs);
  cs_matrix_destroy(&matrix);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmpt, &t_tmp);
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

  cs_timer_t  t_cmpt = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  /* Retrieve high-level structures */
  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t *)scheme_context;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;
  cs_real_t  *vel_c = sc->velocity->val;

  /*----------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  time_eval = t_cur + dt_cur;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  cs_timer_t  t_bld = cs_timer_time();

  /* Build an array storing the Dirichlet values at faces. */
  cs_real_t  *dir_values = NULL;
  cs_cdofb_vecteq_setup_bc(t_cur + dt_cur, mesh, mom_eqp, mom_eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_matrix_structure);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if  (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  cs_real_t  *mass_rhs = NULL;
  BFT_MALLOC(mass_rhs, quant->n_cells, cs_real_t);
  /* Initialization done after */

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, mom_eqp, mom_eqb, mom_eqc, rhs, matrix,        \
         mass_rhs, mav, dir_values, vel_c, sc)
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
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;
    cs_real_t  *div_op = NULL;
    cs_boundary_type_t  *bf_type = NULL;

    BFT_MALLOC(div_op, 3*connect->n_max_fbyc, cs_real_t);
    BFT_MALLOC(bf_type, connect->n_max_fbyc, cs_boundary_type_t);
    cs_cdofb_vecteq_get(&csys, &cb);

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, time_eval, cb);

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
                                       dir_values, vel_c, time_eval,
                                       csys, cb);

      const short int n_fc = cm->n_fc;

      /* Define the local bf_type */
      for (short int i = 0; i < csys->n_bc_faces; i++) {

        /* Get the boundary face in the cell numbering */
        const short int  f = csys->_f_ids[i];
        const cs_lnum_t  bf_id = cm->f_ids[f] - csys->face_shift;

        bf_type[i] = sc->bf_type[bf_id];

      }

      /* Initialize the RHS for the mass eq */
      mass_rhs[c_id] = 0.;

      /* 1- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_diffusion(time_eval, mom_eqp, mom_eqb, mom_eqc,
                                cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion", csys);
#endif

      /* 2- PRESSURE (SCALAR) EQUATION */
      /* ============================= */

      cs_cdofb_navsto_divergence_vect(cm, cb->aux->val);

      cs_real_t *_div = cb->aux->val;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
      if (cs_dbg_cw_test(mom_eqp, cm, csys)) {
        const cs_real_t  ovc = 1. / cm->vol_c;
#       pragma omp critical
        {
          cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
          for (short int f = 0; f < n_fc; f++)
            cs_log_printf(CS_LOG_DEFAULT, "    f%2d: %- .4e, %- .4e, %- .4e\n",
                          f, _div[3*f]*ovc, _div[3*f+1]*ovc, _div[3*f+2]*ovc);
        } /* Critical section */
      }
#endif

      /* -divergence is used in the system */
      for (int i = 0; i < 3*n_fc; i++)
        div_op[i] = -_div[i];

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const _Bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);
      if (has_sourceterm) {

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   time_eval, 1., /* scaling if theta */
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* 3b- OTHER RHS CONTRIBUTIONS *
       * =========================== *
       *
       * First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _apply_bc_partly(sc, mom_eqp, mom_eqc, cm, bf_type, csys, cb);

      /* 4- TIME CONTRIBUTION */
      /* ==================== */

      if (mom_eqb->sys_flag & CS_FLAG_SYS_TIME_DIAG) { /* Mass lumping
                                                          or Hodge-Voronoi */

        const double  ptyc = cb->tpty_val * cm->vol_c * inv_dtcur;

        /* Get cell-cell block */
        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, n_fc, n_fc);

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*n_fc + k] += ptyc * csys->val_n[3*n_fc+k];
          /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        }

      }
      else
        bft_error(__FILE__, __LINE__, 0, " %s: Only diagonal time treatment "
                  "available so far.\n", __func__);

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

      _apply_remaining_bc(sc, mom_eqp, cm, bf_type, csys, cb, div_op,
                          mass_rhs + c_id);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      _assemble(csys, cm, div_op,
                has_sourceterm, mav, rhs, mom_eqc->source_terms);

    } /* Main loop on cells */

    /* Free temporary buffer */
    BFT_FREE(div_op);
    BFT_FREE(bf_type);

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*----------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Now compute/update the velocity and pressure fields */
  _compute_and_update_fields(matrix, sc, mom_eq, rhs, mass_rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
  cs_dbg_fprintf_system(mom_eqp->name, ts->nt_cur, CS_CDOFB_MONOLITHIC_DBG,
                        vel_f, rhs, 3*n_faces);
#endif

  /* Frees */
  BFT_FREE(rhs);
  BFT_FREE(mass_rhs);
  cs_matrix_destroy(&matrix);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmpt, &t_tmp);
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

  cs_timer_t  t_cmpt = cs_timer_time();

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_lnum_t  n_faces = quant->n_faces;

  /* Retrieve high-level structures */
  cs_cdofb_monolithic_t  *sc = (cs_cdofb_monolithic_t *)scheme_context;
  cs_navsto_monolithic_t *cc = (cs_navsto_monolithic_t *)sc->coupling_context;
  cs_equation_t  *mom_eq = cc->momentum;
  cs_cdofb_vecteq_t  *mom_eqc= (cs_cdofb_vecteq_t *)mom_eq->scheme_context;
  cs_equation_param_t *mom_eqp = mom_eq->param;
  cs_equation_builder_t *mom_eqb = mom_eq->builder;
  cs_real_t  *vel_c = sc->velocity->val;

  /*----------------------------------------------------------------------------
   *                      BUILD: START
   *--------------------------------------------------------------------------*/

  const cs_time_step_t *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];
  const cs_real_t  time_eval = t_cur + 0.5*dt_cur;
  const double  tcoef = 1 - mom_eqp->theta;
  const cs_real_t  inv_dtcur = 1./dt_cur;

  cs_timer_t  t_bld = cs_timer_time();

  /* Detect the first call (in this case, we compute the initial source term)*/
  _Bool  compute_initial_source = false;
  if (ts->nt_cur == ts->nt_prev || ts->nt_prev == 0)
    compute_initial_source = true;

  /* Build an array storing the Dirichlet values at faces. */
  cs_real_t  *dir_values = NULL;
  cs_cdofb_vecteq_setup_bc(t_cur + dt_cur, mesh, mom_eqp, mom_eqb, &dir_values);

  /* Initialize the local system: matrix and rhs */
  cs_matrix_t  *matrix = cs_matrix_create(cs_shared_matrix_structure);
  cs_real_t  *rhs = NULL;

  BFT_MALLOC(rhs, 3*n_faces, cs_real_t);
# pragma omp parallel for if  (3*n_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < 3*n_faces; i++) rhs[i] = 0.0;

  cs_real_t  *mass_rhs = NULL;
  BFT_MALLOC(mass_rhs, quant->n_cells, cs_real_t);
  /* Initialization done after */

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav =
    cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, mom_eqp, mom_eqb, mom_eqc, rhs, matrix,        \
         mass_rhs, mav, dir_values, vel_c, sc, compute_initial_source)
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
    cs_cell_sys_t  *csys = NULL;
    cs_cell_builder_t  *cb = NULL;
    cs_real_t  *div_op = NULL;
    cs_boundary_type_t  *bf_type = NULL;

    BFT_MALLOC(div_op, 3*connect->n_max_fbyc, cs_real_t);
    BFT_MALLOC(bf_type, connect->n_max_fbyc, cs_boundary_type_t);
    cs_cdofb_vecteq_get(&csys, &cb);

    /* Store the shift to access border faces (first interior faces and
       then border faces: shift = n_i_faces */
    csys->face_shift = connect->n_faces[CS_INT_FACES];

    /* Initialization of the values of properties */
    cs_equation_init_properties(mom_eqp, mom_eqb, time_eval, cb);

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
                                       dir_values, vel_c, time_eval,
                                       csys, cb);

      const short int n_fc = cm->n_fc;

      /* Define the local bf_type */
      for (short int i = 0; i < csys->n_bc_faces; i++) {

        /* Get the boundary face in the cell numbering */
        const short int  f = csys->_f_ids[i];
        const cs_lnum_t  bf_id = cm->f_ids[f] - csys->face_shift;

        bf_type[i] = sc->bf_type[bf_id];

      }

      /* Initialize the RHS for the mass eq */
      mass_rhs[c_id] = 0.;

      /* 1- VELOCITY (VECTORIAL) EQUATION */
      /* ================================ */

      cs_cdofb_vecteq_diffusion(time_eval, mom_eqp, mom_eqb, mom_eqc,
                                cm, fm, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 1
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> Local system after diffusion", csys);
#endif

      /* 2- PRESSURE (SCALAR) EQUATION */
      /* ============================= */

      cs_cdofb_navsto_divergence_vect(cm, cb->aux->val);

      cs_real_t *_div = cb->aux->val;

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
      if (cs_dbg_cw_test(mom_eqp, cm, csys)) {
        const cs_real_t  ovc = 1. / cm->vol_c;
#       pragma omp critical
        {
          cs_log_printf(CS_LOG_DEFAULT, ">> Divergence:\n");
          for (short int f = 0; f < n_fc; f++)
            cs_log_printf(CS_LOG_DEFAULT, "    f%2d: %- .4e, %- .4e, %- .4e\n",
                          f, _div[3*f]*ovc, _div[3*f+1]*ovc, _div[3*f+2]*ovc);
        } /* Critical section */
      }
#endif

      /* -divergence is used in the system */
      for (int i = 0; i < 3*n_fc; i++)
        div_op[i] = -_div[i];

      /* 3- SOURCE TERM COMPUTATION (for the momentum equation) */
      /* ====================================================== */

      const _Bool  has_sourceterm = cs_equation_param_has_sourceterm(mom_eqp);

      if (has_sourceterm) { /* SOURCE TERM
                             * =========== */
        if (compute_initial_source) { /* First time step */

          cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                     t_cur, tcoef,  /* time, scaling */
                                     cb, mom_eqb, csys);

        }
        else { /* Add the contribution of the previous time step */

          for (short int k = 0; k < 3; k++)
            csys->rhs[3*n_fc + k] += tcoef * mom_eqc->source_terms[3*c_id + k];

        }

        cs_cdofb_vecteq_sourceterm(cm, mom_eqp,
                                   /* time,      scaling */
                                   t_cur+dt_cur, mom_eqp->theta,
                                   cb, mom_eqb, csys);

      } /* End of term source */

      /* 3b- OTHER RHS CONTRIBUTIONS *
       * =========================== *
       *
       * First part of the BOUNDARY CONDITIONS
       *                   ===================
       * Apply a part of BC before the time scheme */

      _apply_bc_partly(sc, mom_eqp, mom_eqc, cm, bf_type, csys, cb);

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
        cs_sdm_t *acc = cs_sdm_get_block(csys->mat, n_fc, n_fc);

        for (short int k = 0; k < 3; k++) {
          csys->rhs[3*n_fc + k] += ptyc * csys->val_n[3*n_fc+k];
          /* Simply add an entry in mat[cell, cell] */
          acc->val[4*k] += ptyc;
        } /* Loop on k */

      }
      else
        bft_error(__FILE__, __LINE__, 0, " %s: Only diagonal time treatment "
                  "available so far.\n", __func__);

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

      _apply_remaining_bc(sc, mom_eqp, cm, bf_type, csys, cb, div_op,
                          mass_rhs + c_id);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 0
      if (cs_dbg_cw_test(mom_eqp, cm, csys))
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
#endif

      /* ************************* ASSEMBLY PROCESS ************************* */

      _assemble(csys, cm, div_op,
                has_sourceterm, mav, rhs, mom_eqc->source_terms);

    } /* Main loop on cells */

    /* Free temporary buffer */
    BFT_FREE(div_op);
    BFT_FREE(bf_type);

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); /* optional */

  /* Free temporary buffers and structures */
  BFT_FREE(dir_values);
  cs_matrix_assembler_values_finalize(&mav);

  /* End of the system building */
  cs_timer_t  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(mom_eqb->tcb), &t_bld, &t_tmp);

  /*----------------------------------------------------------------------------
   *                      BUILD: END
   *--------------------------------------------------------------------------*/

  /* Now compute/update the velocity and pressure fields */
  _compute_and_update_fields(matrix, sc, mom_eq, rhs, mass_rhs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_CDOFB_MONOLITHIC_DBG > 2
  cs_dbg_fprintf_system(mom_eqp->name, ts->nt_cur, CS_CDOFB_MONOLITHIC_DBG,
                        vel_f, rhs, 3*n_faces);
#endif

  /* Frees */
  BFT_FREE(rhs);
  BFT_FREE(mass_rhs);
  cs_matrix_destroy(&matrix);

  t_tmp = cs_timer_time();
  cs_timer_counter_add_diff(&(sc->timer), &t_cmpt, &t_tmp);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
