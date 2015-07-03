/*============================================================================
 * Build the algebraic system for scalar conv./diff. in CDO hybrid cell-based
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
#include <bft_printf.h>

#include "cs_mesh.h"
#include "cs_timer.h"
#include "cs_log.h"
#include "cs_search.h"
#include "cs_prototypes.h"
#include "cs_field.h"
#include "cs_cdo_bc.h"
#include "cs_quadrature.h"
#include "cs_evaluate.h"
#include "cs_hodge.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdofb_codits.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CDOFB_CODITS_DBG 0

struct  _cdofb_codits_t {

  const cs_param_eq_t  *eq;   /* Set of parameters defining the current
                                 system to solve: BCs, material property...
                                 Not owned by the structure.
                                 This structure is freed later */

  _Bool   build_system;       /* false => keep the system as it was at the
                                 previous time step */

  /* Reduced system (known boundary entities are removed --> Dirichlet) */
  cs_lnum_t  n_cells;
  cs_lnum_t  n_faces;
  cs_lnum_t  n_dof_faces; /* Number of interior faces + Neumann faces */

  /* Algebraic system (size = n_dof_vertices) */

  cs_sla_matrix_t  *A;   // matrix to inverse
  double           *x;   // DoF unknows
  double           *rhs; // right-hand side

  /* Boundary conditions:

     face_bc should not change during the simulation.
     The case of a definition of the BCs which changes of type during the
     simulation is possible but not implemented.
     You just have to call the initialization step each time the type of BCs
     is modified to define an updated cs_cdo_bc_t structure.

     We translate Dirichlet BCs to the center of border faces.
     The values can be modified at each time step in case of transient
     simulation.
     For Neumann and Robin BCs,
     BCs contributions are computed on the fly.

   */

  cs_cdo_bc_t        *face_bc;
  cs_real_t          *dir_val; // size: face_bc->dir->n_nhmg_elts

  /* Indirection between zipped numbering (without BC) and initial numbering
     Allocated only if the boundary conditions are strongly enforced.
  */
  cs_lnum_t          *f_z2i_ids;  // Mapping n_dof_faces -> n_faces
  cs_lnum_t          *f_i2z_ids;  // Mapping n_faces     -> n_dof_faces

  /* Work buffer */
  double    *work;

};

/*============================================================================
 * Private variables
 *============================================================================*/

static  int  cs_cdofb_n_scal_systems = 0;
static  cs_cdofb_codits_t  *cs_cdofb_scal_systems = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize structures related to boundary conditions
 *
 * \param[in]     m     pointer to the mesh structure
 * \param[inout]  sys   pointer to a cs_cdofb_codits_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_bc_structures(const cs_mesh_t         *m,
                    cs_cdofb_codits_t       *sys)
{
  int  i;

  const cs_param_bc_t  *bc_param = sys->eq->bc;

  /* Initialize BC information
     We make the distinction betweenn homogeneous and non-homogeneous BCs
  */
  sys->face_bc = cs_cdo_bc_init(bc_param, m->n_b_faces);

  cs_cdo_bc_list_t  *dir_faces = sys->face_bc->dir;

  /* Strong enforcement => compress (or zip) the numbering of vertices */
  sys->f_z2i_ids = NULL; // zipped --> initial ids
  sys->f_i2z_ids = NULL; // initial --> zipped ids

  if (bc_param->strong_enforcement && dir_faces->n_elts > 0) {

    cs_lnum_t  cur_id = 0;
    _Bool  *is_kept = NULL;

    sys->n_dof_faces = sys->n_faces - dir_faces->n_elts;

    BFT_MALLOC(is_kept, sys->n_faces, _Bool);
    for (i = 0; i < sys->n_faces; i++)
      is_kept[i] = true;
    for (i = 0; i < dir_faces->n_elts; i++)
      is_kept[dir_faces->elt_ids[i]] = false;

    /* Build sys->v_z2i_ids and sys->i2i_ids */
    BFT_MALLOC(sys->f_z2i_ids, sys->n_dof_faces, cs_lnum_t);
    BFT_MALLOC(sys->f_i2z_ids, sys->n_faces, cs_lnum_t);

    for (i = 0; i < sys->n_faces; i++) {
      sys->f_i2z_ids[i] = -1;  /* by default, we consider that it's removed */
      if (is_kept[i]) {
        sys->f_i2z_ids[i] = cur_id;
        sys->f_z2i_ids[cur_id++] = i;
      }
    }
    assert(cur_id == dir_faces->n_elts);

    BFT_FREE(is_kept);

    BFT_MALLOC(sys->dir_val, dir_faces->n_nhmg_elts, cs_real_t);
    for (i = 0; i < dir_faces->n_nhmg_elts; i++)
      sys->dir_val[i] = 0.0;

  } /* Strong enforcement of BCs */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize the matrix related to the diffusion op.
 *          Note: values are filled in a second step
 *
 * \param[in]    connect   pointer to a cs_cdo_connect_t structure
 * \param[in]    quant     pointer to a cs_cdo_quantities_t structure
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_diffusion_matrix(const cs_cdo_connect_t     *connect,
                       const cs_cdo_quantities_t  *quant)
{
  int  i, j, shift;

  cs_connect_index_t  *f2f = NULL, *c2f = NULL, *f2c = NULL;

  const cs_lnum_t  n_faces = quant->n_faces;
  const cs_sla_matrix_t *mc2f = connect->c2f;
  const cs_sla_matrix_t *mf2c = connect->f2c;

  /* Allocate and initialize the matrix */
  cs_sla_matrix_t  *mat = cs_sla_matrix_create(n_faces, n_faces, 1,
                                               CS_SLA_MAT_MSR,
                                               false);

  /* Build a face -> face connectivity */
  f2c = cs_index_map(mf2c->n_rows, mf2c->idx, mf2c->col);
  c2f = cs_index_map(mc2f->n_rows, mc2f->idx, mc2f->col);
  f2f = cs_index_convol(n_faces, f2c, c2f);
  cs_index_sort(f2f);
  mat->flag |= CS_SLA_MATRIX_SORTED;

  /* Update index: f2f has the diagonal entry. Remove it for the Hodge index */
  mat->idx[0] = 0;
  for (i = 0; i < n_faces; i++)
    mat->idx[i+1] = mat->idx[i] + f2f->idx[i+1]-f2f->idx[i]-1;

  /* Fill column num */
  BFT_MALLOC(mat->col, mat->idx[n_faces], int);
  shift = 0;
  for (i = 0; i < n_faces; i++)
    for (j = f2f->idx[i]; j < f2f->idx[i+1]; j++)
      if (f2f->lst[j] != i+1)
        mat->col[shift++] = f2f->lst[j];

  /* Sanity check */
  assert(shift == mat->idx[n_faces]);

  /* Free temporary memory */
  cs_index_free(&f2f);
  cs_index_free(&f2c);
  cs_index_free(&c2f);

  /* Allocate and initialize value array */
  for (i = 0; i < n_faces; i++)
    mat->diag[i] = 0.0;

  BFT_MALLOC(mat->val, mat->idx[n_faces], double);
  for (i = 0; i < mat->idx[n_faces]; i++)
    mat->val[i] = 0.0;

  return mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the final (reduced) matrix for diffusion and its right hand
 *          side (RHS). RHS is the sum of three contributions
 *           - Source terms
 *           - Neumann boundary conditions
 *           - Dirichlet boundary conditions
 *
 * \param[in]    m         pointer to a cs_mesh_t structure
 * \param[in]    connect   pointer to a cs_cdo_connect_t struct.
 * \param[in]    quant     pointer to a cs_cdo_quantities_t struct.
 * \param[in]    tcur      current physical time of the simulation
 * \param[inout] sys       pointer to a cs_cdofb_codits_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_build_diffusion_system(const cs_mesh_t             *m,
                        const cs_cdo_connect_t      *connect,
                        const cs_cdo_quantities_t   *quant,
                        cs_real_t                    tcur,
                        cs_cdofb_codits_t           *sys)
{
  int  i, j, ij, c_id;
  double  dsum, rowsum, invdsum;

  double  *BHCtc = NULL; // local size arrays

  cs_sla_matrix_t  *A = _init_diffusion_matrix(connect, quant);

  cs_toolbox_locmat_t  *_h = cs_toolbox_locmat_create(connect->n_max_fbyc);
  cs_toolbox_locmat_t  *_a = cs_toolbox_locmat_create(connect->n_max_fbyc);
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect->n_max_fbyc);

  const cs_lnum_t  n_cells = quant->n_cells;
  const cs_param_eq_t  *eq = sys->eq;
  const cs_param_hodge_t  h_info = eq->diffusion_hodge;
  const cs_cdo_bc_list_t  *dir_faces = sys->face_bc->dir;

  /* Sanity check */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EDFP);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);

  /* Buffers stored in sys->work
     We assume that n_cells <= n_faces for contrib */
  double  *source_term = sys->work;                            // size: n_cells
  double  *contrib = sys->work + sys->n_cells;                 // size: n_faces
  double  *face_rhs = sys->work + sys->n_cells + sys->n_faces; // size: n_faces
  double  *x_bc = sys->work + sys->n_cells + 2*sys->n_faces;   // size: n_faces

  for (i = 0; i < sys->n_cells; i++)
    source_term[i] = 0;

  /* Compute the contribution from source terms */
  if (eq->n_source_terms) { /* Add contribution from source term */

    for (i = 0; i < eq->n_source_terms; i++) {

      const cs_param_source_term_t  st = eq->source_terms[i];

      cs_flag_t  dof_flag =
        CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_PRIMAL | CS_PARAM_FLAG_SCAL;

      /* Sanity check */
      assert(st.var_type == CS_PARAM_VAR_SCAL);

      if (st.type == CS_PARAM_SOURCE_TERM_BASIC ||
          st.type == CS_PARAM_SOURCE_TERM_IMEX) {
        cs_evaluate(m, quant, connect,  // geometrical and topological info.
                    tcur,
                    dof_flag,
                    st.location_id,
                    st.def_type,
                    st.quad_type,
                    st.exp_def,     // definition of the explicit part
                    &contrib);      // updated inside this function

      } // There is an explicit part of the source term to take into account

      for (i = 0; i < sys->n_cells; i++)
        source_term[i] += contrib[i];

    } // Loop on source terms

  } /* There is at least one source term which is defined */

  /*  Build full-size operators:

           n_cells      n_i_faces
         <--------->   <--------->
        | \            * * * * * * |
        |   \           *  BtHC *  |
        |   DcHcUc     * * * * * * |
        |       \       * * * * *  |
    A = |         \                |
        |* * * * * *               |
        | * CHBt * *         H     |
        |* * * * * *               |
        | * * * * *                |

  */

  /* Allocate local operators */
  BFT_MALLOC(BHCtc, connect->n_max_fbyc, double);

  /* Build the remaining discrete operators */
  for (c_id = 0; c_id < n_cells; c_id++) {

    /* Build a local discrete Hodge operator */
    cs_hodge_cost_build_local(c_id, connect, quant, h_info, _h, hb);

    /* Compute dsum = Dc*_H*Uc where Uc = transpose(Dc) and B*H*Ct */
    dsum = 0;
    for (i = 0; i < _h->n_ent; i++) {
      rowsum = 0;
      for (j = 0; j < _h->n_ent; j++)
        rowsum += _h->mat[i*_h->n_ent+j];
      dsum += rowsum;
      BHCtc[i] = -rowsum;
      for (j = 0; j < _a->n_ent; j++)
        _a->mat[_a->n_ent*j+i] = BHCtc[i];
    }
    invdsum = 1/dsum;

    /* Define local diffusion matrix */
    for (i = 0; i < _a->n_ent; i++) {

      double  coef = invdsum*BHCtc[i];

      for (j = 0; j < _a->n_ent; j++) {
        ij = i*_a->n_ent+j;
        _a->mat[ij] *= coef;
        _a->mat[ij] -= _h->mat[ij];
      }

    }

    /* Assemble local matrices */
    cs_sla_assemble_msr(_a, A);

    /* Assemble RHS (source term contribution) */
    for (i = 0; i < _a->n_ent; i++)
      face_rhs[_a->ids[i]] += BHCtc[i]*invdsum*source_term[c_id];

  } /* End of loop on cells */

  /* Clean entries of the operators */
  cs_sla_matrix_clean(A, cs_get_eps_machine());

  /* Free memory */
  BFT_FREE(BHCtc);
  _h = cs_toolbox_locmat_free(_h);
  _a = cs_toolbox_locmat_free(_a);
  hb = cs_hodge_builder_free(hb);

  /* Take into account Dirichlet BCs to update RHS */
  if (dir_faces->n_nhmg_elts > 0) {

    cs_flag_t  dof_flag =
      CS_PARAM_FLAG_FACE | CS_PARAM_FLAG_PRIMAL | CS_PARAM_FLAG_SCAL;

    cs_cdo_bc_dirichlet_set(dof_flag,
                            tcur,
                            quant,
                            eq->bc,
                            dir_faces,
                            sys->dir_val);


    for (i = 0; i < sys->n_faces; i++)
      x_bc[i] = 0;
    for (i = 0; i < dir_faces->n_nhmg_elts; i++)
      x_bc[dir_faces->elt_ids[i]] = sys->dir_val[i];

    cs_sla_matvec(A, x_bc, &contrib, true);
    for (i = 0; i < sys->n_faces; i++)
      face_rhs[i] -= contrib[i];

  }

  /* Reduce system size if there is a strong enforcement of the BCs */
  if (sys->n_dof_faces < sys->n_faces) {

    sys->A = cs_sla_matrix_pack(sys->n_dof_faces,  /* n_block_rows */
                                sys->n_dof_faces,  /* n_block_cols */
                                A,                 /* full matrix */
                                sys->f_z2i_ids,    /* row_p2f_num */
                                sys->f_i2z_ids,    /* col_f2p_num */
                                true);             /* keep sym. */
    A = cs_sla_matrix_free(A);

    /* Define new rhs */
    for (i = 0; i < sys->n_dof_faces; i++)
      sys->rhs[i] = face_rhs[sys->f_z2i_ids[i]];

  }
  else {
    sys->A = A;
    memcpy(sys->rhs, face_rhs, sizeof(cs_real_t)*sys->n_faces);
  }

  /* Prepare usage of iterative solvers */
  cs_sla_matrix_lu_pattern(sys->A);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the required number of scalar equations based on a face
 *         based discretization
 *
 * \param[in]    n_scal_systems   number of scalar equations
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_create_all(int  n_scal_systems)
{
  assert(n_scal_systems > -1);

  BFT_MALLOC(cs_cdofb_scal_systems, n_scal_systems, cs_cdofb_codits_t);

  cs_cdofb_n_scal_systems = n_scal_systems;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdofb_codits_t
 *
 * \param[in] eq       pointer to a structure storing parameters of an eq.
 * \param[in] m        pointer to a mesh structure
 * \param[in] eq_id    id related to the equation to treat
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_init(const cs_param_eq_t  *eq,
                     const cs_mesh_t      *m,
                     int                   eq_id)
{
  cs_lnum_t  i;

  /* Sanity checks */
  assert(eq != NULL);
  assert(eq->space_scheme == CS_SPACE_SCHEME_CDOFB);
  assert(eq->type == CS_PARAM_EQ_TYPE_SCAL);

  if (eq_id < 0 || eq_id >= cs_cdofb_n_scal_systems)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation id (id = %d)\n"), eq_id);

  cs_cdofb_codits_t  *sys = cs_cdofb_scal_systems + eq_id;

  /* Set of parameters related to this algebraic system */
  sys->eq = eq;
  sys->build_system = true; // call the first build

  /* Dimensions: By default, we set number of DoFs as if there is a
     strong enforcement of the BCs */
  sys->n_cells = m->n_cells;
  sys->n_faces = m->n_i_faces + m->n_b_faces;
  sys->n_dof_faces = sys->n_faces;

  /* Boundary conditions (may modify n_dof_faces)
     Set structures related to the management of the BCs */
  _init_bc_structures(m, sys);

  /* Algebraic system */
  BFT_MALLOC(sys->x, sys->n_dof_faces, double);
  BFT_MALLOC(sys->rhs, sys->n_dof_faces, double);
  for (i = 0; i < sys->n_dof_faces; i++)
    sys->x[i] = sys->rhs[i] = 0.0;

  sys->A = NULL; // allocated later

  /* Work buffer */
  assert(sys->n_cells <= sys->n_faces);
  size_t  work_size = sys->n_cells + 3*sys->n_faces;
  BFT_MALLOC(sys->work, work_size, double);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy all cs_cdofb_codits_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_free_all(void)
{
  int  sys_id;

  for (sys_id = 0; sys_id < cs_cdofb_n_scal_systems; sys_id++) {

    cs_cdofb_codits_t  *sys = cs_cdofb_scal_systems + sys_id;

    /* Do not free eq here (not owned by the structure).
       This is done later. */

    /* Free BC structure */
    if (sys->face_bc->dir->n_nhmg_elts > 0)
      BFT_FREE(sys->dir_val);
    sys->face_bc = cs_cdo_bc_free(sys->face_bc);

    /* Renumbering */
    if (sys->n_faces > sys->n_dof_faces) {
      BFT_FREE(sys->f_z2i_ids);
      BFT_FREE(sys->f_i2z_ids);
    }

    /* Free algebraic system */
    BFT_FREE(sys->rhs);
    BFT_FREE(sys->x);
    sys->A = cs_sla_matrix_free(sys->A);

    BFT_FREE(sys->work);

  } // Loop on algebraic systems

  BFT_FREE(cs_cdofb_scal_systems);
  cs_cdofb_scal_systems = NULL;
  cs_cdofb_n_scal_systems = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a scalar convection/diffusion equation with a CDO
 *         face-based scheme.
 *
 * \param[in] m        pointer to a cs_mesh_t structure
 * \param[in] connect  pointer to a cs_cdo_connect_t structure
 * \param[in] quant    pointer to a cs_cdo_quantities_t structure
 * \param[in] tcur     current physical time of the simulation
 * \param[in] sys_id   pointer to a cs_cdofb_codits_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_solve(const cs_mesh_t            *m,
                      const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant,
                      double                      tcur,
                      int                         sys_id)
{
  cs_timer_t  t0, t1;
  cs_timer_counter_t  time_count;
  cs_sla_sumup_t  ret;

  /* Sanity check */
  if (sys_id < 0 || sys_id >= cs_cdofb_n_scal_systems)
    bft_error(__FILE__, __LINE__, 0, _(" Invalid equation id %d\n"), sys_id);

  cs_cdofb_codits_t  *sys = cs_cdofb_scal_systems + sys_id;

  const cs_param_eq_t  *eq = sys->eq;
  const cs_param_itsol_t  itsol_info = eq->itsol_info;

  /* Test to remove */
  if (eq->flag & CS_PARAM_EQ_CONVECTION)
    bft_error(__FILE__, __LINE__, 0,
              _(" Convection is not handled yet.\n"));
  if (eq->flag & CS_PARAM_EQ_UNSTEADY)
    bft_error(__FILE__, __LINE__, 0,
              _(" Transient simulation is not handled yet.\n"));

  /* Build the algebraic system if required */
  t0 = cs_timer_time();

  if (sys->build_system) {

    /* Build diffusion system: stiffness matrix */
    _build_diffusion_system(m, connect, quant, tcur, sys);

    /* Build convection system */
    // TODO

    /* Get information on the matrix related to this linear system */
    cs_sla_matrix_info_t  minfo = cs_sla_matrix_analyse(sys->A);

    bft_printf("\n  <Information about the linear system for %s equation>\n",
               eq->name);
    bft_printf(" -sla- A.size         %d\n", sys->A->n_rows);
    bft_printf(" -sla- A.nnz          %lu\n", minfo.nnz);
    bft_printf(" -sla- A.FillIn       %5.2e %%\n", minfo.fillin);
    bft_printf(" -sla- A.StencilMin   %d\n", minfo.stencil_min);
    bft_printf(" -sla- A.StencilMax   %d\n", minfo.stencil_max);
    bft_printf(" -sla- A.StencilMean  %5.2e\n", minfo.stencil_mean);
  }

  t1 = cs_timer_time();
  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n  -t- Build system >> %s                    %9.3f s\n"),
                eq->name, time_count.wall_nsec*1e-9);
  t0 = cs_timer_time();

  /* Solve system */
  printf("\n# Solve Ax = b with (%s, %s)\n"
         "# System size: %8d ; eps: % -8.5e ;\n",
         cs_param_get_solver_name(itsol_info.solver),
         cs_param_get_precond_name(itsol_info.precond),
         sys->A->n_rows, itsol_info.eps);

  ret = cs_sla_solve(itsol_info, sys->A, sys->rhs, sys->x);

  bft_printf("\n <iterative solver convergence sumup>\n");
  bft_printf(" -sla- code         %s\n", cs_sla_get_code_name(ret.code));
  bft_printf(" -sla- n_iters     %8d\n", ret.iter);
  bft_printf(" -sla- residual   % -8.4e\n", ret.residual);

  t1 = cs_timer_time();
  time_count = cs_timer_diff(&t0, &t1);
  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\n  -t- Linear solver runtime >> %s           %9.3f s\n"),
                eq->name, time_count.wall_nsec*1e-9);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO face-based scheme
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t struct.
 * \param[in]  quant     pointer to a cs_cdo_quantities_t struct.
 * \param[in]  sys_id    id of the equation/system to treat
 */
/*----------------------------------------------------------------------------*/

void
cs_cdofb_codits_post(const cs_cdo_connect_t     *connect,
                     const cs_cdo_quantities_t  *quant,
                     int                         sys_id)
{
  int  i;

  /* Sanity check */
  if (sys_id < 0 || sys_id >= cs_cdofb_n_scal_systems)
    bft_error(__FILE__, __LINE__, 0, _(" Invalid equation id %d\n"), sys_id);

  const cs_cdofb_codits_t  *sys = cs_cdofb_scal_systems + sys_id;
  const cs_cdo_bc_list_t  *dir_faces = sys->face_bc->dir;
  const cs_param_eq_t  *eq = sys->eq;

  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  /* Set computed solution in field array */
  if (sys->n_dof_faces < sys->n_faces)
    for (i = 0; i < sys->n_dof_faces; i++)
      fld->val[sys->f_z2i_ids[i]] = sys->x[i];
  else
    for (i = 0; i < sys->n_faces; i++)
      fld->val[i] = sys->x[i];

  /* Take into account Dirichlet BCs */
  if (eq->bc->strong_enforcement)
    for (i = 0; i < dir_faces->n_nhmg_elts; i++)
      fld->val[dir_faces->elt_ids[i]] = sys->dir_val[i];

  /* /\* Compute cell_field *\/ */
  /* cs_sla_matvec(sys->BHCt, sys->x, &(sys->cell_field), true); */
  /* for (i = 0; i < sys->n_cells; i++) { */
  /*   tmpval = -sys->cell_field[i] + sys->rhs_cell[i]; */
  /*   sys->cell_field[i] = sys->invdiag[i]*tmpval; */
  /* } */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
