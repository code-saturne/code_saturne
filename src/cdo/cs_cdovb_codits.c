/*============================================================================
 * Build an algebraic CDO vertex-based system for convection/diffusion equation
 * with source terms
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
#include <limits.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_mesh.h"
#include "cs_timer.h"
#include "cs_timer_stats.h"
#include "cs_log.h"
#include "cs_search.h"
#include "cs_quadrature.h"
#include "cs_evaluate.h"
#include "cs_prototypes.h"
#include "cs_field.h"
#include "cs_matrix.h"
#include "cs_sles.h"
#include "cs_cdo_bc.h"
#include "cs_hodge.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdovb_codits.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CDOVB_CODITS_DBG 0

/* Algebraic system for CDO vertex-based discretization */

struct _cdovb_codits_t {

  const cs_param_eq_t  *eq;   /* Set of parameters defining the current
                                 system to solve: BCs, material property...
                                 Not owned by the structure.
                                 This structure is freed later */

  int     main_ts_id;         /* Id of the main timer states structure related
                                 to this equation */
  int     pre_ts_id;          /* Id of the timer stats structure gathering all
                                 steps before the resolution of the linear
                                 systems */
  int     solve_ts_id;       /* Id of the timer stats structure related
                                to the inversion of the linear system */
  int     post_ts_id;        /* Id of the timer stats structure gathering all
                                steps afterthe resolution of the linear systems
                                (post, balance...) */

  _Bool   build_system;       /* false => keep the system as it was at the
                                 previous time step */

  /* System size (known boundary entities may be removed if BCs are strongly
     enforced) */

  cs_lnum_t  n_vertices;
  cs_lnum_t  n_dof_vertices; /* n_rows = n_cols = n_vertices - dir. vertices */

  /* Algebraic system (size = n_dof_vertices) */

  cs_matrix_structure_t       *ms;  /* matrix structure (how are stored
                                       coefficients of the matrix a) */
  cs_matrix_t                 *a;   // matrix to inverse with cs_sles_solve()
  cs_real_t                   *x;   // DoF unknows
  cs_real_t                   *rhs; // right-hand side

  /* Boundary conditions:

     face_bc should not change during the simulation.
     The case of a definition of the BCs which changes of type during the
     simulation is possible but not implemented.
     You just have to call the initialization step each time the type of BCs
     is modified to define an updated cs_cdo_bc_t structure.

     We translate Dirichlet BCs to border vertices
     The values can be modified at each time step in case of transient
     simulation.
     For Neumann and Robin BCs, the treatment is more complex since the
     contributions of these BCs to a dual cell related to a border vertex is
     computed from the contribution of different portions of primal faces
     (these primal border faces form a closure of the dual cell).
     These contributions are computed on the fly.

   */

  cs_cdo_bc_t        *face_bc;
  cs_cdo_bc_list_t   *vtx_dir;
  double             *dir_val; /* size = vtx_dir->n_nhmg_elts
                                  allocated if there is a strong enforcement
                                  of the BCs */

  /* Indirection between zipped numbering (without BC) and initial numbering
     Allocated only if the boundary conditions are strongly enforced.
  */
  cs_lnum_t          *v_z2i_ids;  // Mapping n_dof_vertices -> n_vertices
  cs_lnum_t          *v_i2z_ids;  // Mapping n_vertices     -> n_dof_vertices

  /* Work buffer */
  cs_real_t    *work;

};

/*============================================================================
 * Private variables
 *============================================================================*/

static  int  cs_cdovb_n_scal_systems = 0;
static  cs_cdovb_codits_t  *cs_cdovb_scal_systems = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize structures related to boundary conditions
 *
 * \param[in]     m     pointer to the mesh structure
 * \param[inout]  sys   pointer to a cs_cdovb_codits_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_init_bc_structures(const cs_mesh_t         *m,
                    cs_cdovb_codits_t       *sys)
{
  int  i;

  const cs_param_bc_t  *bc_param = sys->eq->bc;

  /* Initialize BC information
     Compute also the list of Dirichlet vertices and their related definition
     We make the distinction betweenn homogeneous and non-homogeneous BCs
  */
  sys->face_bc = cs_cdo_bc_init(bc_param, m->n_b_faces);
  sys->vtx_dir = cs_cdo_bc_vtx_dir_create(m, sys->face_bc);

  /* Strong enforcement => compress (or zip) the numbering of vertices */
  sys->v_z2i_ids = NULL; // zipped --> initial ids
  sys->v_i2z_ids = NULL; // initial --> zipped ids

  if (bc_param->strong_enforcement && sys->vtx_dir->n_elts > 0) {

    cs_lnum_t  cur_id = 0;
    _Bool  *is_kept = NULL;

    sys->n_dof_vertices = sys->n_vertices - sys->vtx_dir->n_elts;

    BFT_MALLOC(is_kept, sys->n_vertices, _Bool);
    for (i = 0; i < sys->n_vertices; i++)
      is_kept[i] = true;
    for (i = 0; i < sys->vtx_dir->n_elts; i++)
      is_kept[sys->vtx_dir->elt_ids[i]] = false;

    /* Build sys->v_z2i_ids and sys->i2i_ids */
    BFT_MALLOC(sys->v_z2i_ids, sys->n_dof_vertices, cs_lnum_t);
    BFT_MALLOC(sys->v_i2z_ids, sys->n_vertices, cs_lnum_t);

    for (i = 0; i < sys->n_vertices; i++) {
      sys->v_i2z_ids[i] = -1;  /* by default, we consider that it's removed */
      if (is_kept[i]) {
        sys->v_i2z_ids[i] = cur_id;
        sys->v_z2i_ids[cur_id++] = i;
      }
    }
    assert(cur_id == sys->n_dof_vertices);

    BFT_FREE(is_kept);

    /* Allocate dir_val */
    BFT_MALLOC(sys->dir_val, sys->vtx_dir->n_nhmg_elts, double);

  } /* Strong enforcement of BCs */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the right hand side (rhs) for a convection/diffusion
 *          equation.
 *          Take into account several contributions:
 *          --> Dirichlet, Neumann and Robin BCs
 *          --> Source terms
 *
 * \param[in]     m         pointer to a cs_mesh_t structure
 * \param[in]     connect   pointer to a cs_cdo_connect_t structure
 * \param[in]     quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]     ad        pointer to the (full) diffusion matrix
 * \param[in]     tcur      current physical time of the simulation
 * \param[inout]  sys       pointer to a cs_cdovb_codits_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_rhs(const cs_mesh_t            *m,
             const cs_cdo_connect_t     *connect,
             const cs_cdo_quantities_t  *quant,
             const cs_sla_matrix_t      *ad,
             double                      tcur,
             cs_cdovb_codits_t          *sys)
{
  int  i;

  const cs_param_eq_t  *eq = sys->eq;
  const cs_param_bc_t  *bc = eq->bc;
  const cs_cdo_bc_list_t  *vtx_dir = sys->vtx_dir;

  // only scalar equation are handled up to now (TODO)
  assert(eq->type == CS_PARAM_EQ_TYPE_SCAL);

  /* Initialize rhs */
  double *full_rhs = sys->work;

  for (i = 0; i < sys->n_vertices; i++)
    full_rhs[i] = 0.0;

  /* BOUNDARY CONDITIONS */

  if (vtx_dir->n_nhmg_elts > 0) { /* Dirichlet BC */

    cs_flag_t  dof_flag =
      CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL | CS_PARAM_FLAG_SCAL;

    cs_cdo_bc_dirichlet_set(dof_flag,
                            tcur,
                            m,
                            bc,
                            vtx_dir,
                            sys->dir_val);

    if (bc->strong_enforcement) {

      double  *x_bc = sys->work + sys->n_vertices;
      double  *contrib = sys->work + 2*sys->n_vertices;

      for (i = 0; i < sys->n_vertices; i++)
        x_bc[i] = 0.0;
      for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
        x_bc[vtx_dir->elt_ids[i]] = sys->dir_val[i];

      /* Compute ad*Tbc: rhs = rhs - ad*Tbc */
      cs_sla_matvec(ad, x_bc, &contrib, true);
      for (i = 0; i < sys->n_vertices; i++)
        full_rhs[i] -= contrib[i];

    }
    else { /* Take into account Dirichlet BC: Define x_bc */

      for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
        full_rhs[vtx_dir->elt_ids[i]] += bc->penalty_coef * sys->dir_val[i];

    }

  } /* Dirichlet BC */

  /* TODO: Add contribution for Neumann BC (if homogeneous nothing to do)
     and Robin BC */

  /* SOURCE TERMS */

  if (eq->n_source_terms > 0) { /* Add contribution from source term */

    for (i = 0; i < eq->n_source_terms; i++) {

      const cs_param_source_term_t  st = eq->source_terms[i];

      double  *contrib = sys->work + 2*sys->n_vertices;
      cs_flag_t  dof_flag =
        CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_DUAL | CS_PARAM_FLAG_SCAL;

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
                    st.exp_def,         // definition of the explicit part
                    &contrib);          // updated inside this function

        /* Update full rhs */
        for (i = 0; i < sys->n_vertices; i++)
          full_rhs[i] += contrib[i];

      } // There is an explicit part of the source term to take into account

    } // Loop on source terms

  } /* There is at least one source term which is defined */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a stiffness matrix
 *
 * \param[in]  n_vertices  number of vertices
 * \param[in]  connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to an initialized stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_stiffness(cs_lnum_t                n_vertices,
                const cs_cdo_connect_t  *connect)
{
  int  i, j, shift;

  cs_connect_index_t  *v2v = NULL, *v2c = NULL;

  const cs_connect_index_t  *c2v = connect->c2v;

  cs_sla_matrix_t  *s = cs_sla_matrix_create(n_vertices,
                                             n_vertices,
                                             1,     // stride
                                             CS_SLA_MAT_MSR,
                                             true); // symmetry

  /* Initialize index (v2v connectivity) */
  v2c = cs_index_transpose(n_vertices, c2v);
  v2v = cs_index_compose(n_vertices, v2c, c2v);
  cs_index_free(&v2c);

  /* Sort index */
  cs_index_sort(v2v);
  s->flag |= CS_SLA_MATRIX_SORTED;

  /* Update index (v2v has a diagonal entry. We remove it -> MSR) */
  s->idx[0] = 0;
  for (i = 0; i < n_vertices; i++)
    s->idx[i+1] = s->idx[i] + v2v->idx[i+1]-v2v->idx[i]-1;

  /* Fill column id */
  BFT_MALLOC(s->col_id, s->idx[n_vertices], cs_lnum_t);
  shift = 0;
  for (i = 0; i < n_vertices; i++)
    for (j = v2v->idx[i]; j < v2v->idx[i+1]; j++)
      if (v2v->ids[j] != i)
        s->col_id[shift++] = v2v->ids[j];

  /* Sanity check */
  assert(shift == s->idx[n_vertices]);

  /* Free temporary memory */
  cs_index_free(&v2v);

  /* Allocate and initialize value array */
  for (i = 0; i < n_vertices; i++)
    s->diag[i] = 0.0;

  BFT_MALLOC(s->val, s->idx[n_vertices], double);
  for (i = 0; i < s->idx[n_vertices]; i++)
    s->val[i] = 0.0;

  return s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Encode E_c cap E_v sets for all vertices of a cell using a mask
 *          of bits
 *
 * \param[in]     c_id      cell id
 * \param[in]     n_bits    number of bits in a block
 * \param[in]     n_blocks  number of blocks in a mask
 * \param[in]     connect   point to a cs_cdo_connect_t struct.
 * \param[in]     vtag      buffer of tag indicating the local vertex id
 * \param[inout]  masks     list of masks
 */
/*----------------------------------------------------------------------------*/

static void
_encode_edge_masks(cs_lnum_t                c_id,
                   int                      n_bits,
                   int                      n_blocks,
                   const cs_cdo_connect_t  *connect,
                   const short int          vtag[],
                   cs_flag_t                masks[])
{
  short int  ie;
  cs_lnum_t  i;

  const cs_connect_index_t  *c2e = connect->c2e;

  /* Reset masks */
  for (i = 0; i < n_blocks*connect->n_max_vbyc; i++)
    masks[i] = 0;

  /* Loop on cell edges */
  for (ie = 0, i = c2e->idx[c_id]; i < c2e->idx[c_id+1]; ie++, i++) {

    cs_lnum_t  eshft = 2*c2e->ids[i];
    short int  vi = vtag[connect->e2v->col_id[eshft]];
    short int  vj = vtag[connect->e2v->col_id[eshft+1]];

    /* Sanity checks */
    assert(vi > -1 && vi < connect->n_max_vbyc);
    assert(vj > -1 && vj < connect->n_max_vbyc);

    /* Encode this edge in the set E_v1,c and E_v2,c */
    short int  b = ie/n_bits, s = ie%n_bits;
    masks[vi*n_blocks + b] |= (1 << s);
    masks[vj*n_blocks + b] |= (1 << s);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the full stiffness matrix and the matrix related to the
 *          discrete Hodge operator from a cellwise assembly process
 *
 * \param[in]    connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]    quant       pointer to a cs_cdo_quantities_t struct.
 * \param[inout] sys         pointer to a cs_cdovb_codits_t struct.
 *
 * \return a pointer to the full stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_build_stiffness_matrix(const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *quant,
                        cs_cdovb_codits_t          *sys)
{
  int  i, j, n_ent, n_blocks, n_bits, mask_size;
  cs_lnum_t  c_id, v_id;

  short int  *vtag = NULL;
  cs_flag_t  *emsk = NULL;

  cs_toolbox_locmat_t  *hloc = cs_toolbox_locmat_create(connect->n_max_ebyc);
  cs_toolbox_locmat_t  *sloc = cs_toolbox_locmat_create(connect->n_max_vbyc);
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect->n_max_ebyc);
  cs_sla_matrix_t  *s = _init_stiffness(quant->n_vertices, connect);

  const cs_connect_index_t  *c2v = connect->c2v;
  const cs_param_eq_t  *eq = sys->eq;
  const cs_param_hodge_t  h_info = eq->diffusion_hodge;

  /* Sanity check */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EPFD);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);

  /* Initialize mask features */
  n_bits = sizeof(cs_flag_t)*CHAR_BIT;
  n_blocks = connect->n_max_ebyc/n_bits;
  if (connect->n_max_ebyc % n_bits != 0) n_blocks++;
  mask_size = n_blocks*n_bits;

  BFT_MALLOC(emsk, n_blocks*connect->n_max_vbyc, cs_flag_t);
  for (i = 0; i < n_blocks*connect->n_max_vbyc; i++)
    emsk[i] = 0;

  /* Initialize tags */
  BFT_MALLOC(vtag, quant->n_vertices, short int);
  for (i = 0; i < quant->n_vertices; i++)
    vtag[i] = -1;

  for (c_id = 0; c_id < quant->n_cells; c_id++) {   /*** Loop on cells ***/

    short int  ek, el, vi, sgn_ik, sgn_jl;

    /* Build the local Hodge op. */
    cs_hodge_cost_build_local(c_id, connect, quant, h_info, hloc, hb);

    /* Initialize vertex tag and the local stiffness matrix */
    for (i = 0, j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; i++, j++) {
      v_id = c2v->ids[j];
      vtag[v_id] = i;
      sloc->ids[i] = v_id;
    }
    n_ent = i;
    sloc->n_ent = n_ent;

    for (i = 0; i < n_ent*n_ent; i++)
      sloc->mat[i] = 0;

    /* Encode edge sets related to each vertex belonging to the current cell */
    _encode_edge_masks(c_id, n_bits, n_blocks, connect, vtag, emsk);

    /* Build local stiffness matrix */
    for (vi = 0; vi < sloc->n_ent; vi++) {

      short int pos_i = vi*sloc->n_ent;

      for (ek = 0; ek < hloc->n_ent; ek++) {  /* Find edges attached to vi */

        short int b = ek/n_bits, r = ek % n_bits;

        if (emsk[vi*n_blocks+b] & (1 << r)) { // ek in E_i

          cs_lnum_t  ek_shft = 2*hloc->ids[ek];
          int  pos_ek = ek*hloc->n_ent;

          if (connect->e2v->col_id[ek_shft] == sloc->ids[vi])
            sgn_ik = connect->e2v->sgn[ek_shft];
          else {
            assert(connect->e2v->col_id[ek_shft+1] == sloc->ids[vi]);
            sgn_ik = connect->e2v->sgn[ek_shft+1];
          }

          for (el = 0; el < hloc->n_ent; el++) { /* Loop on cell edges */

            cs_lnum_t  el_shft = 2*hloc->ids[el];
            double  val = hloc->mat[pos_ek+el]*sgn_ik;
            cs_lnum_t  v1_id = connect->e2v->col_id[el_shft];
            short int  vj1 = vtag[v1_id];

            sgn_jl = connect->e2v->sgn[el_shft];
            sloc->mat[pos_i+vj1] += val*sgn_jl;

            cs_lnum_t  v2_id = connect->e2v->col_id[el_shft+1];
            short int  vj2 = vtag[v2_id];

            sgn_jl = connect->e2v->sgn[el_shft+1];
            sloc->mat[pos_i+vj2] += val*sgn_jl;

          } /* Loop on cell edges */

        } /* ek in E_i */

      } /* Find edges attached to _vi */

    } /* Loop on vertices vi */

    /* Assemble stiffness matrix */
    cs_sla_assemble_msr(sloc, s);

  } /* Loop on cells */

  /* Free temporary buffers and structures */
  BFT_FREE(vtag);
  BFT_FREE(emsk);

  hloc = cs_toolbox_locmat_free(hloc);
  sloc = cs_toolbox_locmat_free(sloc);
  hb = cs_hodge_builder_free(hb);

  return s;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define the matrix "Ad" for the diffusion term and its right hand
 *          side (rhs).
 *          rhs potentially collects these contributions:
 *           - Source terms
 *           - Neumann boundary conditions
 *           - Robin boundary conditions (TODO) [the explicit part]
 *           - Dirichlet boundary conditions : full Ad*Tbc
 *
 * \param[in]    m         pointer to a cs_mesh_t structure
 * \param[in]    connect   pointer to a cs_cdo_connect_t struct.
 * \param[in]    quant     pointer to a cs_cdo_quantities_t struct.
 * \param[in]    tcur      current physical time of the simulation
 * \param[inout] sys       pointer to a cs_cdovb_codits_t struct.
 *
 * \return a pointer to a cs_sla_matrix_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_build_diffusion_system(const cs_mesh_t            *m,
                        const cs_cdo_connect_t     *connect,
                        const cs_cdo_quantities_t  *quant,
                        double                      tcur,
                        cs_cdovb_codits_t          *sys)
{
  int  i;

  cs_sla_matrix_t  *full_matrix = NULL, *final_matrix = NULL;

  const cs_param_eq_t  *eq = sys->eq;

  /* Sanity check */
  assert(eq->bc->strong_enforcement == true); // TODO

  /* Build the (full) stiffness matrix i.e. without taking into account BCs */
  full_matrix = _build_stiffness_matrix(connect, quant, sys);

  /* Remove entries very small with respect to other coefficients */
  cs_sla_matrix_clean(full_matrix, cs_get_eps_machine());

  /* Compute the full rhs */
  _compute_rhs(m, connect, quant, full_matrix, tcur, sys);

  double  *full_rhs = sys->work; // stored in sys->work

  if (sys->n_vertices == sys->n_dof_vertices) { // Keep the full system

    final_matrix = full_matrix;
    memcpy(sys->rhs, full_rhs, sys->n_vertices*sizeof(double));

  }
  else { /* Reduce size. Extract block with degrees of freedom */

    for (i = 0; i < sys->n_dof_vertices; i++)
      sys->rhs[i] = full_rhs[sys->v_z2i_ids[i]];

    final_matrix = cs_sla_matrix_pack(sys->n_dof_vertices,  // n_block_rows
                                      sys->n_dof_vertices,  // n_block_cols
                                      full_matrix,          // full matrix
                                      sys->v_z2i_ids,       // row_z2i_ids
                                      sys->v_i2z_ids,       // col_i2z_ids
                                      true);                // keep sym.

    /* Free buffers */
    full_matrix = cs_sla_matrix_free(full_matrix);

  }

  //  cs_sla_matrix_dump("Stiffness.log", NULL, final_matrix); // DBG
  return final_matrix;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Switch the matrix represenation from a cs_sla_matrix_t struct. to
 *          a cs_matrix_t struct.
 *          sla_mat is freed inside this routine.
 *
 * \param[inout] sys       pointer to a cs_cdovb_codits_t struct.
 * \param[inout] sla_mat   pointer to a cs_sla_matrix_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_map_to_matrix(cs_cdovb_codits_t          *sys,
               cs_sla_matrix_t            *sla_mat)
{
  /* Sanity check */
  assert(sla_mat->type == CS_SLA_MAT_MSR);

  /* First step: create a matrix structure */
  sys->ms =  cs_matrix_structure_create_msr(CS_MATRIX_MSR,      // type
                                            true,               // transfer
                                            true,               // have_diag
                                            sla_mat->n_rows,    // n_rows
                                            sla_mat->n_cols,    // n_cols_ext
                                            &(sla_mat->idx),    // row_index
                                            &(sla_mat->col_id), // col_id
                                            NULL,               // halo
                                            NULL);              // numbering

  sys->a = cs_matrix_create(sys->ms); // ms is also stored inside a

  const cs_lnum_t  *row_index, *col_id;
  cs_matrix_get_msr_arrays(sys->a, &row_index, &col_id, NULL, NULL);

  /* Second step: associate coefficients to a matrix structure */
  cs_matrix_transfer_coefficients_msr(sys->a,
                                      false,             // symmetric values ?
                                      NULL,              // diag. block
                                      NULL,              // extra-diag. block
                                      row_index,         // row_index
                                      col_id,            // col_id
                                      &(sla_mat->diag),  // diag. values
                                      &(sla_mat->val));  // extra-diag. values

  /* Free non-transferred parts of sla_mat */
  sla_mat = cs_sla_matrix_free(sla_mat);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Solve a linear system Ax = b arising from a CDO vertex-based scheme
 *
 * \param[in]    eq        pointer to a cs_param_eq_t structure
 * \param[in]    a         pointer to a cs_matrix_t struct.
 * \param[in]    rhs       right-hand side
 * \param[inout] x         array of unkowns
 */
/*----------------------------------------------------------------------------*/

static void
_solve_linear_system(const cs_param_eq_t        *eq,
                     const cs_matrix_t          *a,
                     const cs_real_t            *rhs,
                     cs_real_t                  *x)
{
  double  r_norm;
  cs_sles_convergence_state_t  cvg;
  cs_sla_sumup_t  ret;

  cs_halo_rotation_t  halo_rota = CS_HALO_ROTATION_IGNORE;

  cs_sles_t  *sles = cs_sles_find_or_add(eq->field_id, NULL);

  const cs_lnum_t  n_rows = cs_matrix_get_n_rows(a);
  const cs_param_itsol_t  itsol_info = eq->itsol_info;

  printf("\n# Solve Ax = b for %s with %s\n"
         "# System size: %8d ; eps: % -8.5e ;\n",
         eq->name, cs_param_get_solver_name(itsol_info.solver),
         n_rows, itsol_info.eps);

  if (itsol_info.resid_normalized)
    r_norm = cs_euclidean_norm(n_rows, rhs) / n_rows;
  else
    r_norm = 1.0;

  cvg = cs_sles_solve(sles,
                      a,
                      halo_rota,
                      itsol_info.eps,
                      r_norm,
                      &(ret.iter),
                      &(ret.residual),
                      rhs,
                      x,
                      0,      // aux. size
                      NULL);  // aux. buffers

  bft_printf("\n <iterative solver convergence sumup>\n");
  bft_printf(" -sla- code        %d\n", cvg);
  bft_printf(" -sla- n_iters     %d\n", ret.iter);
  bft_printf(" -sla- residual    % -8.4e\n", ret.residual);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate the required number of scalar equations based on a vertex
 *         based discretization
 *
 * \param[in]    n_scal_systems   number of scalar equations
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_create_all(int  n_scal_systems)
{
  assert(n_scal_systems > -1);

  BFT_MALLOC(cs_cdovb_scal_systems, n_scal_systems, cs_cdovb_codits_t);

  cs_cdovb_n_scal_systems = n_scal_systems;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_codits_t
 *
 * \param[in] eq       pointer to a structure storing parameters of an eq.
 * \param[in] m        pointer to a mesh structure
 * \param[in] eq_id    id related to the equation/system to treat
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_init(const cs_param_eq_t  *eq,
                     const cs_mesh_t      *m,
                     int                   eq_id)
{
  cs_lnum_t  i, len;

  /* Sanity checks */
  assert(eq != NULL);
  assert(eq->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(eq->type == CS_PARAM_EQ_TYPE_SCAL);

  if (eq_id < 0 || eq_id >= cs_cdovb_n_scal_systems)
    bft_error(__FILE__, __LINE__, 0,
              _(" Invalid equation id (id = %d)\n"), eq_id);

  cs_cdovb_codits_t  *sys = cs_cdovb_scal_systems + eq_id;

  /* Set of parameters related to this algebraic system */
  sys->eq = eq;

  /* Create timer statistic structure for a simplified profiling */
  sys->main_ts_id = sys->pre_ts_id = sys->solve_ts_id = sys->post_ts_id = -1;

  if (eq->verbosity > 0) {
    sys->main_ts_id = cs_timer_stats_create("stages", // parent name
                                            eq->name,
                                            eq->name);

    cs_timer_stats_start(sys->main_ts_id);
    cs_timer_stats_set_plot(sys->main_ts_id, 1);

    if (eq->verbosity > 1) {

      char *label = NULL;

      len = strlen("_solve") + strlen(eq->name) + 1;
      BFT_MALLOC(label, len, char);

      sprintf(label, "%s_pre", eq->name);
      sys->pre_ts_id = cs_timer_stats_create(eq->name, label, label);
      cs_timer_stats_set_plot(sys->pre_ts_id, 1);

      sprintf(label, "%s_solve", eq->name);
      sys->solve_ts_id = cs_timer_stats_create(eq->name, label, label);
      cs_timer_stats_set_plot(sys->solve_ts_id, 1);

      sprintf(label, "%s_post", eq->name);
      sys->post_ts_id = cs_timer_stats_create(eq->name, label, label);
      cs_timer_stats_set_plot(sys->post_ts_id, 1);

      BFT_FREE(label);

    } // verbosity > 1

  } // verbosity > 0

  sys->build_system = true; // call the first build

  /* Dimensions: By default, we set number of DoFs as if there is a
     weak enforcement of the BCs */
  sys->n_vertices = m->n_vertices;
  sys->n_dof_vertices = sys->n_vertices;

  /* Boundary conditions (may modify n_dof_vertices)
     Set structures related to the management of the BCs */
  _init_bc_structures(m, sys);

  /* Algebraic system */
  BFT_MALLOC(sys->x, sys->n_dof_vertices, cs_real_t);
  BFT_MALLOC(sys->rhs, sys->n_dof_vertices, cs_real_t);
  for (i = 0; i < sys->n_dof_vertices; i++)
    sys->x[i] = 0.0, sys->rhs[i] = 0.;

  /* A matrix is composed of a matrix structure and a coefficient structure
     both are allocated later */
  sys->ms = NULL;
  sys->a = NULL;

  /* Work buffer */
  size_t  work_size = 3*sys->n_vertices;
  BFT_MALLOC(sys->work, work_size, cs_real_t);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy all cs_cdovb_codits_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_free_all(void)
{
  int  sys_id;

  for (sys_id = 0; sys_id < cs_cdovb_n_scal_systems; sys_id++) {

    cs_cdovb_codits_t  *sys = cs_cdovb_scal_systems + sys_id;

    /* Do not free eq here (not owned by the structure).
       This is done later. */

    /* Free BC structure */
    if (sys->eq->bc->strong_enforcement && sys->vtx_dir->n_elts > 0)
      BFT_FREE(sys->dir_val);

    sys->face_bc = cs_cdo_bc_free(sys->face_bc);
    sys->vtx_dir = cs_cdo_bc_list_free(sys->vtx_dir);

    /* Renumbering */
    if (sys->n_vertices > sys->n_dof_vertices) {
    BFT_FREE(sys->v_z2i_ids);
    BFT_FREE(sys->v_i2z_ids);
    }

    /* Free algebraic system */
    BFT_FREE(sys->rhs);
    BFT_FREE(sys->x);

    cs_matrix_structure_destroy(&(sys->ms));
    cs_matrix_destroy(&(sys->a));

    BFT_FREE(sys->work);

    if (sys->main_ts_id > -1)
      cs_timer_stats_stop(sys->main_ts_id);
  }

  BFT_FREE(cs_cdovb_scal_systems);
  cs_cdovb_scal_systems = NULL;
  cs_cdovb_n_scal_systems = 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve a scalar convection/diffusion equation with a CDO
 *         vertex-based scheme.
 *
 * \param[in]  m        pointer to a cs_mesh_t structure
 * \param[in]  connect  pointer to a cs_cdo_connect_t structure
 * \param[in]  quant    pointer to a cs_cdo_quantities_t structure
 * \param[in]  tcur     current physical time of the simulation
 * \param[in]  eq_id    pointer to a cs_cdovb_codits_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_solve(const cs_mesh_t            *m,
                      const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant,
                      double                      tcur,
                      int                         eq_id)
{
  /* Sanity check */
  if (eq_id < 0 || eq_id >= cs_cdovb_n_scal_systems)
    bft_error(__FILE__, __LINE__, 0, _(" Invalid equation id %d\n"), eq_id);

  cs_cdovb_codits_t  *sys = cs_cdovb_scal_systems + eq_id;

  const cs_param_eq_t  *eq = sys->eq;

  /* Test to remove */
  if (eq->flag & CS_PARAM_EQ_CONVECTION)
    bft_error(__FILE__, __LINE__, 0,
              _(" Convection is not handled yet.\n"));
  if (eq->flag & CS_PARAM_EQ_UNSTEADY)
    bft_error(__FILE__, __LINE__, 0,
              _(" Transient simulation is not handled yet.\n"));

  if (sys->build_system) {

    if (sys->pre_ts_id > -1)
      cs_timer_stats_start(sys->pre_ts_id);

    /* Build diffusion system: stiffness matrix */
    cs_sla_matrix_t  *sla_mat =  _build_diffusion_system(m, connect, quant,
                                                         tcur,
                                                         sys);

    /* Build convection system */
    // TODO

    /* Get information on the matrix related to this linear system */
    cs_sla_matrix_info_t  minfo = cs_sla_matrix_analyse(sla_mat);

    bft_printf("\n  <Information about the linear system for %s equation>\n",
               eq->name);
    bft_printf(" -sla- A.size         %d\n", sla_mat->n_rows);
    bft_printf(" -sla- A.nnz          %lu\n", minfo.nnz);
    bft_printf(" -sla- A.FillIn       %5.2e %%\n", minfo.fillin);
    bft_printf(" -sla- A.StencilMin   %d\n", minfo.stencil_min);
    bft_printf(" -sla- A.StencilMax   %d\n", minfo.stencil_max);
    bft_printf(" -sla- A.StencilMean  %5.2e\n", minfo.stencil_mean);

    /* Map sla_mat (cs_sla_matrix_t struct.) into sys->a (cs_matrix_t struct.)
       sla_mat is freed during this call */
    _map_to_matrix(sys, sla_mat);

    if (sys->pre_ts_id > -1)
      cs_timer_stats_stop(sys->pre_ts_id);

  }

  /* Solve system */
  if (sys->solve_ts_id > -1)
    cs_timer_stats_start(sys->solve_ts_id);

  _solve_linear_system(eq, sys->a, sys->rhs, sys->x);

  if (sys->solve_ts_id > -1)
    cs_timer_stats_stop(sys->solve_ts_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO vertex-based scheme.
 *
 * \param[in]  connect   pointer to a cs_cdo_connect_t struct.
 * \param[in]  quant     pointer to a cs_cdo_quantities_t struct.
 * \param[in]  eq_id     id of the equation/system to treat
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_post(const cs_cdo_connect_t     *connect,
                     const cs_cdo_quantities_t  *quant,
                     int                         eq_id)
{
  int  i;

  /* Sanity check */
  if (eq_id < 0 || eq_id >= cs_cdovb_n_scal_systems)
    bft_error(__FILE__, __LINE__, 0, _(" Invalid equation id %d\n"), eq_id);

  const  cs_cdovb_codits_t  *sys = cs_cdovb_scal_systems + eq_id;
  const cs_cdo_bc_list_t  *v_dir = sys->vtx_dir;
  const cs_param_eq_t  *eq = sys->eq;

  cs_field_t  *fld = cs_field_by_id(eq->field_id);

  if (sys->post_ts_id > -1)
    cs_timer_stats_start(sys->post_ts_id);

  /* Set computed solution in field array */
  if (sys->n_dof_vertices < sys->n_vertices)
    for (i = 0; i < sys->n_dof_vertices; i++)
      fld->val[sys->v_z2i_ids[i]] = sys->x[i];
  else
    for (i = 0; i < sys->n_vertices; i++)
      fld->val[i] = sys->x[i];

  /* Set BC in field array if we have this knowledge */
  if (eq->bc->strong_enforcement)
    for (i = 0; i < v_dir->n_nhmg_elts; i++)
      fld->val[v_dir->elt_ids[i]] = sys->dir_val[i];

  if (sys->post_ts_id > -1)
    cs_timer_stats_stop(sys->post_ts_id);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
