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

#include "cs_log.h"
#include "cs_search.h"
#include "cs_quadrature.h"
#include "cs_evaluate.h"
#include "cs_cdo_bc.h"
#include "cs_hodge.h"
#include "cs_cdo_convection.h"

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

struct _cs_cdovb_codits_t {

  /* Pointer to a cs_equation_param_t structure shared with a cs_equation_t
     structure. This structure is not owner. */

  const cs_equation_param_t  *eqp;

  /* System size (known boundary entities may be removed if BCs are strongly
     enforced) */

  cs_lnum_t  n_vertices;
  cs_lnum_t  n_dof_vertices; /* n_rows = n_cols = n_vertices - dir. vertices */

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

  cs_param_bc_enforce_t  enforce; // type of enforcement of BCs
  cs_cdo_bc_t           *face_bc; // list of faces sorted by type of BCs
  cs_cdo_bc_list_t      *vtx_dir; // list of vertices attached to a Dirichlet BC
  double                *dir_val; /* size = vtx_dir->n_nhmg_elts
                                     allocated if there is a strong enforcement
                                     of the BCs */

  /* Indirection between zipped numbering (without BC) and initial numbering
     Allocated only if the boundary conditions are strongly enforced.
  */
  cs_lnum_t          *v_z2i_ids;  // Mapping n_dof_vertices -> n_vertices
  cs_lnum_t          *v_i2z_ids;  // Mapping n_vertices     -> n_dof_vertices

  /* Work buffer */
  cs_real_t          *work;

};

/*============================================================================
 * Private variables
 *============================================================================*/

// Advanceed developper parameters
static const cs_real_t  cs_weak_nitsche_pena_coef = 500;
static const cs_real_t  cs_weak_penalization_weight = 0.01;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a local dense matrix-vector product
 *          matvec has been previously allocated
 *
 * \param[in]      loc    local matrix to use
 * \param[in]      vec    local vector to use
 * \param[in, out] matvec result of the local matrix-vector product
 */
/*----------------------------------------------------------------------------*/

static void
_loc_matvec(const cs_toolbox_locmat_t  *loc,
            const cs_real_t            *vec,
            cs_real_t                  *matvec)
{
  /* Sanity checks */
  assert(loc != NULL && vec != NULL && matvec != NULL);

  cs_real_t  v = vec[0];
  for (int j = 0; j < loc->n_ent; j++)
    matvec[j] = v*loc->mat[j*loc->n_ent];

  for (int i = 1; i < loc->n_ent; i++) {
    v = vec[i];
    for (int j = 0; j < loc->n_ent; j++)
      matvec[j] += v*loc->mat[j*loc->n_ent + i];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a new matrix by adding a local matrix with its transpose.
 *          Keep the transposed matrix for future use/
 *
 * \param[in, out] loc   local matrix to transpose and add
 * \param[in, out] tr    transposed of the local matrix
 */
/*----------------------------------------------------------------------------*/

static void
_add_transpose(cs_toolbox_locmat_t  *loc,
               cs_toolbox_locmat_t  *tr)
{
  /* Sanity check */
  assert(loc != NULL && tr != NULL && tr->n_max_ent == loc->n_max_ent);

  if (loc->n_ent < 1)
    return;

  tr->n_ent = loc->n_ent;

  for (int i = 0; i < loc->n_ent; i++) {

    int  ii = i*loc->n_ent + i;

    tr->ids[i] = loc->ids[i];
    tr->mat[ii] = loc->mat[ii];
    loc->mat[ii] *= 2;

    for (int j = i+1; j < loc->n_ent; j++) {

      int  ij = i*loc->n_ent + j;
      int  ji = j*loc->n_ent + i;

      tr->mat[ji] = loc->mat[ij];
      tr->mat[ij] = loc->mat[ji];
      loc->mat[ij] += tr->mat[ij];
      loc->mat[ji] += tr->mat[ji];

    }
  }

}
/*----------------------------------------------------------------------------*/
/*!
 * \brief   Update the stiffness matrix "Ad" for the diffusion term and
 *          its right hand side (rhs).
 *          rhs potentially collects these contributions:
 *           - Source terms
 *           - Neumann boundary conditions
 *           - Robin boundary conditions (TODO) [the explicit part]
 *           - Dirichlet boundary conditions : full Ad*Tbc
 *
 * \param[in]      m            pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      builder      pointer to a cs_cdovb_codits_t structure
 * \param[in]      tcur         current physical time of the simulation
 * \param[in, out] full_matrix  matrix of the linear system
 * \param[in, out] full_rhs     right-hand side
 */
/*----------------------------------------------------------------------------*/

static void
_enforce_weak_nitsche(const cs_mesh_t            *m,
                      const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant,
                      const cs_cdovb_codits_t    *builder,
                      cs_real_t                   tcur,
                      cs_sla_matrix_t            *full_matrix,
                      cs_real_t                   full_rhs[])
{
  cs_lnum_t  i, ie, j, k, _id;
  cs_real_t  surf, eig_ratio, eig_max, beta;
  cs_real_3_t  xyz, mn, reco_val;
  cs_real_33_t  matpty;

  short int  *loc_v_id = NULL, *loc_e_id = NULL;
  cs_real_3_t  *dfv = NULL, *pev = NULL; // dual face and primal edge vectors
  cs_real_t  *v_coef = NULL, *over_pec_vol = NULL, *dir_vals = NULL;
  cs_real_t  *_vec = NULL, *_matvec = NULL;
  cs_toolbox_locmat_t  *ntrgrd = cs_toolbox_locmat_create(connect->n_max_vbyc);
  cs_toolbox_locmat_t  *transp = NULL;

  const cs_sla_matrix_t  *e2v = connect->e2v;
  const cs_sla_matrix_t  *f2e = connect->f2e;
  const cs_cdo_bc_list_t  *face_dir = builder->face_bc->dir;
  const cs_equation_param_t  *eqp = builder->eqp;
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const cs_connect_index_t  *c2e = connect->c2e;
  const cs_connect_index_t  *c2v = connect->c2v;

  /* Parameters linked to the discrete Hodge operator used for diffusion */
  bool  is_uniform = cs_param_pty_is_uniform(h_info.pty_id);

  if (h_info.algo != CS_PARAM_HODGE_ALGO_COST)
    bft_error(__FILE__, __LINE__, 0,
              " Only COST algorithm for hodge operators is possible when\n"
              " a weak enforcement of the boundary conditions is requested.");

  beta = h_info.coef;

  if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

    const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;

    BFT_MALLOC(dir_vals, quant->n_vertices, cs_real_t);
    for (i = 0; i < quant->n_vertices; i++)
      dir_vals[i] = 0.;
    for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
      dir_vals[vtx_dir->elt_ids[i]] = builder->dir_val[i];

    BFT_MALLOC(_vec, connect->n_max_vbyc, cs_real_t);
    BFT_MALLOC(_matvec, connect->n_max_vbyc, cs_real_t);

    transp = cs_toolbox_locmat_create(connect->n_max_vbyc);
  }

  /* Initialize v_coef */
  BFT_MALLOC(v_coef, quant->n_vertices, cs_real_t);
  BFT_MALLOC(loc_v_id, quant->n_vertices, short int);
  for (i = 0; i < quant->n_vertices; i++) {
    v_coef[i] = 0.0;
    loc_v_id[i] = -1;
  }

  BFT_MALLOC(loc_e_id, quant->n_edges, short int);
  for (i = 0; i < quant->n_edges; i++)
    loc_e_id[i] = -1;

  /* Store locally dual face vector for a quicker access */
  BFT_MALLOC(dfv, connect->n_max_ebyc, cs_real_3_t);
  BFT_MALLOC(pev, connect->n_max_ebyc, cs_real_3_t);
  BFT_MALLOC(over_pec_vol, connect->n_max_ebyc, cs_real_t);

  /* Get the anisotropic ratio and the max. eigenvalue (if uniform) */
  if (is_uniform) {
    xyz[0] = xyz[1] = xyz[2] = 0.;
    cs_evaluate_pty(h_info.pty_id,
                    tcur,       // when ?
                    xyz,        // where ?
                    h_info.inv_pty,
                    &matpty);

    cs_eigen_mat33(matpty, &eig_ratio, &eig_max);
  }

  /* Loop on Dirichlet faces */
  for (i = 0; i < face_dir->n_elts; i++) {

    cs_lnum_t  f_id = face_dir->elt_ids[i] + quant->n_i_faces;
    cs_quant_t  qf = quant->face[f_id];
    cs_lnum_t  c_id = connect->f2c->col_id[connect->f2c->idx[f_id]];
    cs_real_t  over_c_vol = 1/quant->cell_vol[c_id];

    /* Sanity check (this is a border face) */
    assert(connect->f2c->idx[f_id+1] - connect->f2c->idx[f_id] == 1);

    if (!is_uniform) {
      cs_evaluate_pty(h_info.pty_id,
                      tcur,                           // when ?
                      &(quant->cell_centers[3*c_id]), // where ?
                      h_info.inv_pty,
                      &matpty);

      cs_eigen_mat33(matpty, &eig_ratio, &eig_max);
    }

    cs_real_t  f_coef = pow(qf.meas, -0.5) * eig_ratio * eig_max;

    /* Compute the product: matpty*face unit normal */
    _mv3(matpty, qf.unitv, &mn);

    /* Define an id local to this cell for each vertex */
    for (j = c2v->idx[c_id], _id = 0; j < c2v->idx[c_id+1]; j++, _id++) {
      cs_lnum_t  v_id = c2v->ids[j];
      loc_v_id[v_id] = _id;
      ntrgrd->ids[_id] = v_id;
    }

    /* Initialize buffers related to vertices */
    ntrgrd->n_ent = _id;

    /* Reset local trace normal*gradient operator */
    for (j = 0; j < ntrgrd->n_ent*ntrgrd->n_ent; j++)
      ntrgrd->mat[j] = 0.0;

    /* Define an id local to this cell and store useful quantities
       for each edge */
    for (j = c2e->idx[c_id], _id = 0; j < c2e->idx[c_id+1]; j++, _id++) {

      cs_lnum_t  e_id = c2e->ids[j];
      cs_quant_t qpe = quant->edge[e_id];
      cs_dface_t qdf = quant->dface[j];

      loc_e_id[e_id] = _id;
      for (k = 0; k < 3; k++) {
        pev[_id][k] = qpe.unitv[k]*qpe.meas;
        dfv[_id][k] = qdf.vect[k];
      }

      over_pec_vol[_id] = 3./_dp3(pev[_id], dfv[_id]);

    } // Loop on cell edges

    /* Loop on border face edges */
    for (j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

      cs_lnum_t  e_id = f2e->col_id[j];
      cs_quant_t  qe = quant->edge[e_id];
      cs_lnum_t  e_shft = e2v->idx[e_id];
      cs_lnum_t  v1_id = e2v->col_id[e_shft];
      cs_lnum_t  v2_id = e2v->col_id[e_shft+1];

      short int  _e = loc_e_id[e_id];
      short int  _v1 = loc_v_id[v1_id];
      short int  _v2 = loc_v_id[v2_id];

      /* Sanity checks */
      assert(_e != -1 && _v1 != -1 && _v2 != -1);

      surf = cs_surftri(&(m->vtx_coord[3*v1_id]), qe.center, qf.center);
      v_coef[v1_id] += surf*f_coef;
      v_coef[v2_id] += surf*f_coef;

      for (ie = c2e->idx[c_id]; ie < c2e->idx[c_id+1]; ie++) {

        cs_lnum_t  ek_id = c2e->ids[ie];
        cs_lnum_t  ek_shft = e2v->idx[ek_id];
        cs_lnum_t  vj1_id = e2v->col_id[ek_shft];
        short int  sgn_j1 = e2v->sgn[ek_shft];
        cs_lnum_t  vj2_id = e2v->col_id[ek_shft+1];
        short int  sgn_j2 = e2v->sgn[ek_shft+1];

        short int  _ek = loc_e_id[ek_id];
        short int  _vj1 = loc_v_id[vj1_id];
        short int  _vj2 = loc_v_id[vj2_id];

        /* Compute L_Ec(GRAD(p_j)) on t_{ek,f} */
        if (_ek == _e)
          for (k = 0; k < 3; k++)
            reco_val[k] = dfv[_ek][k]*beta*over_pec_vol[_ek];
        else
          for (k = 0; k < 3; k++)
            reco_val[k] = 0.0;

        cs_real_t  dp = beta*over_pec_vol[_e]*_dp3(dfv[_ek], pev[_e]);

        for (k = 0; k < 3; k++)
          reco_val[k] += over_c_vol*(dfv[_ek][k] - dp*dfv[_e][k]);

        cs_real_t  contrib = _dp3(mn, reco_val)*surf;

        /* Update the coefficients of the local normal trace operator */
        ntrgrd->mat[_v1*ntrgrd->n_ent + _vj1] += sgn_j1 * contrib;
        ntrgrd->mat[_v1*ntrgrd->n_ent + _vj2] += sgn_j2 * contrib;
        ntrgrd->mat[_v2*ntrgrd->n_ent + _vj1] += sgn_j1 * contrib;
        ntrgrd->mat[_v2*ntrgrd->n_ent + _vj2] += sgn_j2 * contrib;

      } // Loop on cell edges

    } // border face edges

    /* Assemble contributions coming from local normal trace operator */
    if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {

      /* ntrgrd = ntrgrd + transp and transp = transpose(ntrgrd) */
      _add_transpose(ntrgrd, transp);

      cs_sla_assemble_msr_sym(ntrgrd, full_matrix);

      /* Modify RHS according to the add of transp */
      for (j = 0; j < ntrgrd->n_ent; j++)
        _vec[j] = dir_vals[ntrgrd->ids[j]];

      _loc_matvec(transp, _vec, _matvec);

      for (j = 0; j < ntrgrd->n_ent; j++)
        full_rhs[ntrgrd->ids[j]] += _matvec[j];

    }
    else
      cs_sla_assemble_msr(ntrgrd, full_matrix);


  } // Dirichlet faces

  /* Loop on Dirichlet vertices (first, nhmg vertices then hmg vertices)
     to assemble the contribution coming from the diagonal Hodge penalization */
  for (i = 0; i < builder->vtx_dir->n_elts; i++) {

    cs_lnum_t  v_id = builder->vtx_dir->elt_ids[i];
    cs_real_t  pena_coef = cs_weak_nitsche_pena_coef*v_coef[v_id];

    /* Sanity check */
    assert(full_matrix->type == CS_SLA_MAT_MSR);

    full_matrix->diag[v_id] += pena_coef;
    if (i < builder->vtx_dir->n_nhmg_elts)
      full_rhs[v_id] += pena_coef*builder->dir_val[i];

  }

  /* Free memory */
  if (builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_SYM) {
    transp = cs_toolbox_locmat_free(transp);

    BFT_FREE(_vec);
    BFT_FREE(_matvec);
    BFT_FREE(dir_vals);
  }

  ntrgrd = cs_toolbox_locmat_free(ntrgrd);

  BFT_FREE(v_coef);
  BFT_FREE(loc_v_id);
  BFT_FREE(loc_e_id);
  BFT_FREE(dfv);
  BFT_FREE(pev);
  BFT_FREE(over_pec_vol);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the (full) right hand side (rhs) from the contributions of
 *          source terms
 *
 * \param[in]      m         pointer to a cs_mesh_t structure
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      tcur      current physical time of the simulation
 * \param[in, out] builder   pointer to a cs_cdovb_codits_t structure
 * \param[in, out] rhs       pointer to a pointer of rhs values
 */
/*----------------------------------------------------------------------------*/

static void
_compute_st_rhs(const cs_mesh_t            *m,
                const cs_cdo_connect_t     *connect,
                const cs_cdo_quantities_t  *quant,
                double                      tcur,
                cs_cdovb_codits_t          *builder,
                cs_real_t                 **rhs)
{
  int  i;

  const cs_equation_param_t  *eqp = builder->eqp;

  if (*rhs == NULL)
    BFT_MALLOC(*rhs, builder->n_vertices, cs_real_t);

  /* Initialize rhs */
  double *full_rhs = *rhs;

  for (i = 0; i < builder->n_vertices; i++)
    full_rhs[i] = 0.0;

  if (eqp->n_source_terms > 0) { /* Add contribution from source terms */

    for (i = 0; i < eqp->n_source_terms; i++) {

      const cs_param_source_term_t  st = eqp->source_terms[i];

      double  *contrib = builder->work + builder->n_vertices;
      cs_flag_t  dof_flag =
        CS_PARAM_FLAG_CELL | CS_PARAM_FLAG_DUAL | CS_PARAM_FLAG_SCAL;

      /* Sanity check */
      assert(st.var_type == CS_PARAM_VAR_SCAL);
      assert(st.type == CS_PARAM_SOURCE_TERM_EXPLICIT);

      cs_evaluate(m, quant, connect,  // geometrical and topological info.
                  tcur,
                  dof_flag,
                  st.ml_id,
                  st.def_type,
                  st.quad_type,
                  st.use_subdiv,
                  st.def,             // definition of the explicit part
                  &contrib);          // updated inside this function

      /* Update full rhs */
      for (i = 0; i < builder->n_vertices; i++)
        full_rhs[i] += contrib[i];

    } // Loop on source terms

  } /* There is at least one source term which is defined */

}


/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize a (full diffusion) matrix
 *
 * \param[in]  n_vertices  number of vertices
 * \param[in]  connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to an initialized (stiffness) matrix
 */
/*----------------------------------------------------------------------------*/

static cs_sla_matrix_t *
_init_matrix(cs_lnum_t                n_vertices,
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
 * \param[in]       c_id      cell id
 * \param[in]       n_bits    number of bits in a block
 * \param[in]       n_blocks  number of blocks in a mask
 * \param[in]       connect   point to a cs_cdo_connect_t struct.
 * \param[in]       vtag      buffer of tag indicating the local vertex id
 * \param[in, out]  masks     list of masks
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
 * \brief   Define the full stiffness matrix from a cellwise assembly process
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant       pointer to a cs_cdo_quantities_t struct.
 * \param[in, out] builder     pointer to a cs_cdovb_codits_t struct.
 * \param[in, out] s           pointer to a cs_sla_matrix_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_build_diffusion(const cs_cdo_connect_t     *connect,
                 const cs_cdo_quantities_t  *quant,
                 cs_cdovb_codits_t          *builder,
                 cs_sla_matrix_t            *s)
{
  int  i, j, n_ent, n_blocks, n_bits, mask_size;
  cs_lnum_t  c_id, v_id;

  short int  *vtag = NULL;
  cs_flag_t  *emsk = NULL;

  cs_toolbox_locmat_t  *hloc = cs_toolbox_locmat_create(connect->n_max_ebyc);
  cs_toolbox_locmat_t  *sloc = cs_toolbox_locmat_create(connect->n_max_vbyc);
  cs_hodge_builder_t  *hb = cs_hodge_builder_init(connect->n_max_ebyc);

  const cs_connect_index_t  *c2v = connect->c2v;
  const cs_equation_param_t  *eqp = builder->eqp;
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;

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

  for (c_id = 0; c_id < quant->n_cells; c_id++) {   /* Loop on cells */

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
    cs_sla_assemble_msr_sym(sloc, s);

  } /* Loop on cells */

  /* Free temporary buffers and structures */
  BFT_FREE(vtag);
  BFT_FREE(emsk);

  hloc = cs_toolbox_locmat_free(hloc);
  sloc = cs_toolbox_locmat_free(sloc);
  hb = cs_hodge_builder_free(hb);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Update the matrix of the linear system and the right hand side to
 *          take into account the BCs for the diffusion term
 *
 * \param[in]      m         pointer to a cs_mesh_t structure
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      tcur      current physical time of the simulation
 * \param[in, out] builder   pointer to a cs_cdovb_codits_t structure
 * \param[in, out] rhs       pointer of pointer to the right-hand side
 * \param[in, out] matrix    pointer to the pointer of a matrix structure
 */
/*----------------------------------------------------------------------------*/

static void
_update_diffusion_with_bc(const cs_mesh_t            *m,
                          const cs_cdo_connect_t     *connect,
                          const cs_cdo_quantities_t  *quant,
                          double                      tcur,
                          cs_cdovb_codits_t          *builder,
                          cs_real_t                 **rhs,
                          cs_sla_matrix_t           **matrix)
{
  int  i;

  cs_sla_matrix_t  *full_matrix = *matrix, *final_matrix = NULL;
  double  *full_rhs = *rhs, *final_rhs = NULL;

  const cs_cdo_bc_list_t  *vtx_dir = builder->vtx_dir;

  /* Enforce the Dirchlet boundary conditions on vertices */
  switch (builder->enforce) {

  case CS_PARAM_BC_ENFORCE_STRONG:
    {
      double  *_rhs = builder->work;

      if (vtx_dir->n_nhmg_elts > 0) {

        double  *x_bc = builder->work + builder->n_vertices;
        double  *contrib = builder->work + 2*builder->n_vertices;

        for (i = 0; i < builder->n_vertices; i++)
          x_bc[i] = 0.0;
        for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
          x_bc[vtx_dir->elt_ids[i]] = builder->dir_val[i];

        /* Compute full_matrix*Tbc: rhs = rhs - full_matrix*Tbc */
        cs_sla_matvec(full_matrix, x_bc, &contrib, true);
        for (i = 0; i < builder->n_vertices; i++)
          full_rhs[i] -= contrib[i];

      }

      if (builder->n_vertices > builder->n_dof_vertices) {

        memcpy(_rhs, full_rhs, builder->n_vertices*sizeof(double));
        final_rhs = *rhs;
        for (i = 0; i < builder->n_dof_vertices; i++)
          final_rhs[i] = _rhs[builder->v_z2i_ids[i]];

        /* Reduce size. Extract block with degrees of freedom */
        final_matrix = cs_sla_matrix_pack(builder->n_dof_vertices, // n_pack_rows
                                          builder->n_dof_vertices, // n_pack_cols
                                          full_matrix,             // full matrix
                                          builder->v_z2i_ids,      // row_z2i_ids
                                          builder->v_i2z_ids,      // col_i2z_ids
                                          true);                   // keep sym.

        /* Free buffers */
        full_matrix = cs_sla_matrix_free(full_matrix);

      }
    }
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_PENA:
    {
      /* Sanity check */
      assert(builder->n_vertices == builder->n_dof_vertices);

      cs_real_t  pena_coef = cs_weak_penalization_weight/cs_get_eps_machine();

      for (i = 0; i < builder->vtx_dir->n_elts; i++)
        full_matrix->diag[builder->vtx_dir->elt_ids[i]] += pena_coef;

      for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
        full_rhs[vtx_dir->elt_ids[i]] += pena_coef * builder->dir_val[i];

      // DBG
      printf(" Weak enforcement with penalization coefficient %5.3e\n",
             pena_coef);
    }
    break;

  case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
  case CS_PARAM_BC_ENFORCE_WEAK_SYM:
    assert(builder->n_vertices == builder->n_dof_vertices); /* Sanity check */
    _enforce_weak_nitsche(m, connect, quant, builder, tcur,
                          full_matrix,
                          full_rhs);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " This kind of BC enforcement is not implemented yet.\n"
              " Please modify your settings.");

  } // End of switch (enforcement)

  /* TODO: Add contribution for Neumann BC (if homogeneous nothing to do)
     and Robin BC */

  if (builder->n_vertices == builder->n_dof_vertices) { // Keep the full system
    final_matrix = full_matrix;
    final_rhs = full_rhs;
  }

  /* Return pointers */
  *matrix = final_matrix;
  *rhs = final_rhs;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Update the matrix of the linear system to take into account the
 *          convection term
 *
 * \param[in]      m            pointer to a cs_mesh_t structure
 * \param[in]      connect      pointer to a cs_cdo_connect_t structure
 * \param[in]      quant        pointer to a cs_cdo_quantities_t structure
 * \param[in]      tcur         current physical time of the simulation
 * \param[in, out] sys_builder  pointer to a cs_cdovb_codits_t structure
 * \param[in, out] rhs          pointer of pointer to the right-hand side
 * \param[in, out] matrix       pointer to a matrix structure
 */
/*----------------------------------------------------------------------------*/

static void
_build_convection(const cs_mesh_t            *m,
                  const cs_cdo_connect_t     *connect,
                  const cs_cdo_quantities_t  *quant,
                  double                      tcur,
                  cs_cdovb_codits_t          *sys_builder,
                  cs_real_t                  *rhs,
                  cs_sla_matrix_t            *matrix)
{
  cs_lnum_t  i, j, k, f_id;
  cs_qvect_t  beta;
  cs_real_3_t  beta_uni, beta_g, xg = {0, 0, 0};
  cs_real_t  dp = 0, flux = 0;

  cs_lnum_t  *pos_dir = NULL;

  const cs_real_t  one_third = 1./3;
  const cs_real_t  *xyz = m->vtx_coord;
  const cs_sla_matrix_t  *e2v = connect->e2v;
  const cs_sla_matrix_t  *f2e = connect->f2e;
  const cs_cdo_bc_list_t  *vtx_dir = sys_builder->vtx_dir;
  const cs_equation_param_t  *eqp = sys_builder->eqp;
  const cs_param_hodge_t  h_info = eqp->diffusion_hodge;
  const cs_param_advection_t  a_info = eqp->advection;

  bool  adv_uniform = cs_param_adv_field_is_uniform(a_info.adv_id);
  bool  with_diffusion = (eqp->flag & CS_EQUATION_DIFFUSION) ? true : false;

  cs_convection_builder_t  *adv_builder = cs_convection_builder_init(connect);

  /* Build the first part of the cellwise convection operator and assemble
     it into the system matrix */
  if (with_diffusion)
    cs_convection_with_diffusion(connect, quant,
                                 tcur,
                                 a_info, h_info,
                                 adv_builder,
                                 matrix);
  else
    cs_convection(connect, quant,
                  tcur,
                  a_info,
                  adv_builder,
                  matrix);

  adv_builder = cs_convection_builder_free(adv_builder);

  /* Update now the algebraic according to the boundary conditions */
  if (adv_uniform) {
    cs_evaluate_adv_field(a_info.adv_id, tcur, xg, &beta_uni);
    cs_qvect(beta_uni, &beta);
  }

  /* Store the position in dir_val of vertices with a non-homogeneous
     Dirichlet BC */
  BFT_MALLOC(pos_dir, quant->n_vertices, cs_lnum_t);
  for (i = 0; i < quant->n_vertices; i++)
    pos_dir[i] = -1; // not a vertex with non-homogeneous Dirichlet BC
  for (i = 0; i < vtx_dir->n_nhmg_elts; i++)
    pos_dir[vtx_dir->elt_ids[i]] = i;

  cs_real_t  pcoef = 0;
  if (sys_builder->enforce == CS_PARAM_BC_ENFORCE_WEAK_PENA)
    pcoef = cs_weak_penalization_weight/cs_get_eps_machine();

  /* Add diagonal term for vertices attached to a boundary face where
     the advection field points inward */
  for (f_id = quant->n_i_faces; f_id < quant->n_faces; f_id++) {

    cs_quant_t  qf = quant->face[f_id];

    /* Sanity check (this is a border face) */
    assert(connect->f2c->idx[f_id+1] - connect->f2c->idx[f_id] == 1);

    if (adv_uniform)
      dp = _dp3(beta.unitv, qf.unitv);

    /* Loop on border face edges */
    for (i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

      cs_lnum_t  e_id = f2e->col_id[i];
      cs_quant_t  qe = quant->edge[e_id];
      cs_lnum_t  e_shft = e2v->idx[e_id];

      for (j = 0; j < 2; j++) {

        cs_lnum_t  v_id = e2v->col_id[e_shft+j];

        /* Could be improved with a more accurate computation of the flux */
        if (!adv_uniform) {
          for (k = 0; k < 3; k++)
            xg[k] = one_third * (xyz[3*v_id+k] + qe.center[k] + qf.center[k]);
          cs_evaluate_adv_field(a_info.adv_id, tcur, xg, &beta_g);
          cs_qvect(beta_g, &beta);
          dp = _dp3(beta.unitv, qf.unitv);
        }

        flux = dp * beta.meas * cs_surftri(&(xyz[3*v_id]), qe.center, qf.center);

        if (dp < 0) { // advection field is inward w.r.t. the face normal

          /* Weakly enforcement of Dirichlet BCs */
          assert(pos_dir[v_id] > -1);

          cs_real_t  bc_val = sys_builder->dir_val[pos_dir[v_id]];

          switch (sys_builder->enforce) {

          case CS_PARAM_BC_ENFORCE_WEAK_PENA:
            rhs[v_id] += pcoef * bc_val;
            matrix->diag[v_id] +=  pcoef;
            break;

          case CS_PARAM_BC_ENFORCE_WEAK_NITSCHE:
          case CS_PARAM_BC_ENFORCE_WEAK_SYM:
            rhs[v_id] -= flux * bc_val;
            if (a_info.form == CS_PARAM_ADVECTION_FORM_NONCONS)
              matrix->diag[v_id] -= flux;
            break;

            default:
              bft_error(__FILE__, __LINE__, 0,
                        " Invalid type of BC enforcement when advection is"
                        " activated.");

            } // switch

        }
        else {  // advection is oriented outward

          if (a_info.form == CS_PARAM_ADVECTION_FORM_CONSERV)
            matrix->diag[v_id] += flux;

        }

      } // Loop on edge vertices

    } // Loop on face edges

  } // Loop on border faces

  cs_sla_matrix_clean(matrix, 10*cs_get_eps_machine());

  /* Free memory */
  BFT_FREE(pos_dir);
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdovb_codits_t structure
 *
 * \param[in] eqp      pointer to a cs_equation_param_t structure
 * \param[in] m        pointer to a mesh structure
 *
 * \return a pointer to a new allocated cs_cdovb_codits_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_cdovb_codits_init(const cs_equation_param_t  *eqp,
                     const cs_mesh_t            *m)
{
  cs_lnum_t  i;

  /* Sanity checks */
  assert(eqp != NULL);
  assert(eqp->space_scheme == CS_SPACE_SCHEME_CDOVB);
  assert(eqp->type == CS_EQUATION_TYPE_SCAL);

  cs_cdovb_codits_t  *builder = NULL;

  BFT_MALLOC(builder, 1, cs_cdovb_codits_t);

  builder->eqp = eqp;

  /* Dimensions:
     By default, we set number of DoFs as if there is a weak enforcement of
     the boundary conditions */
  builder->n_vertices = m->n_vertices;
  builder->n_dof_vertices = builder->n_vertices;

  /* Set members and structures related to the management of the BCs */

  const cs_param_bc_t  *bc_param = eqp->bc;

  /* Translate user-defined information about BC into a structure well-suited
     for computation. We make the distinction between homogeneous and
     non-homogeneous BCs.
     We also compute also the list of Dirichlet vertices along with their
     related definition.
  */
  builder->face_bc = cs_cdo_bc_init(bc_param, m->n_b_faces);
  builder->vtx_dir = cs_cdo_bc_vtx_dir_create(m, builder->face_bc);

  /* Allocate dir_val */
  BFT_MALLOC(builder->dir_val, builder->vtx_dir->n_nhmg_elts, double);
  for (i = 0; i < builder->vtx_dir->n_nhmg_elts; i++)
    builder->dir_val[i] = 0.0;

  /* Strong enforcement means that we need an indirection list between the
     compress (or zip) and initial numbering of vertices */
  builder->enforce = bc_param->enforcement;
  builder->v_z2i_ids = NULL; // zipped --> initial ids
  builder->v_i2z_ids = NULL; // initial --> zipped ids

  if (builder->enforce == CS_PARAM_BC_ENFORCE_STRONG &&
      builder->vtx_dir->n_elts > 0) {

    cs_lnum_t  cur_id = 0;
    bool  *is_kept = NULL;

    builder->n_dof_vertices = builder->n_vertices - builder->vtx_dir->n_elts;

    BFT_MALLOC(is_kept, builder->n_vertices, bool);
    for (i = 0; i < builder->n_vertices; i++)
      is_kept[i] = true;
    for (i = 0; i < builder->vtx_dir->n_elts; i++)
      is_kept[builder->vtx_dir->elt_ids[i]] = false;

    /* Build builder->v_z2i_ids and builder->i2i_ids */
    BFT_MALLOC(builder->v_z2i_ids, builder->n_dof_vertices, cs_lnum_t);
    BFT_MALLOC(builder->v_i2z_ids, builder->n_vertices, cs_lnum_t);

    for (i = 0; i < builder->n_vertices; i++) {
      /* by default, we consider that it's removed */
      builder->v_i2z_ids[i] = -1;
      if (is_kept[i]) {
        builder->v_i2z_ids[i] = cur_id;
        builder->v_z2i_ids[cur_id++] = i;
      }
    }
    assert(cur_id == builder->n_dof_vertices);

    BFT_FREE(is_kept);

  } /* Strong enforcement of BCs */

  /* Work buffer */
  if (builder->enforce == CS_PARAM_BC_ENFORCE_STRONG)
    BFT_MALLOC(builder->work, 3*builder->n_vertices, cs_real_t);
  else
    BFT_MALLOC(builder->work, 2*builder->n_vertices, cs_real_t);

  return builder;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_cdovb_codits_t structure
 *
 * \param[in, out]  builder   pointer to a cs_cdovb_codits_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdovb_codits_free(void   *builder)
{
  cs_cdovb_codits_t  *_builder = (cs_cdovb_codits_t *)builder;

  if (_builder == NULL)
    return _builder;

  /* Free BC structure */
  if (_builder->vtx_dir->n_nhmg_elts > 0)
    BFT_FREE(_builder->dir_val);

  _builder->face_bc = cs_cdo_bc_free(_builder->face_bc);
  _builder->vtx_dir = cs_cdo_bc_list_free(_builder->vtx_dir);

  /* Renumbering (if strong enforcement of BCs for instance) */
  if (_builder->n_vertices > _builder->n_dof_vertices) {
    BFT_FREE(_builder->v_z2i_ids);
    BFT_FREE(_builder->v_i2z_ids);
  }

  BFT_FREE(_builder->work);
  BFT_FREE(_builder);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a CDO vertex-based scheme.
 *
 * \param[in]      m        pointer to a cs_mesh_t structure
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to a cs_cdo_quantities_t structure
 * \param[in]      tcur     current physical time of the simulation
 * \param[in, out] builder  pointer to cs_cdovb_codits_t structure
 * \param[in, out] rhs      right-hand side
 * \param[in, out] sla_mat  pointer to cs_sla_matrix_t structure pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_build_system(const cs_mesh_t            *m,
                             const cs_cdo_connect_t     *connect,
                             const cs_cdo_quantities_t  *quant,
                             double                      tcur,
                             void                       *builder,
                             cs_real_t                 **rhs,
                             cs_sla_matrix_t           **sla_mat)
{
  cs_cdovb_codits_t  *sys_builder = (cs_cdovb_codits_t *)builder;

  const cs_cdo_bc_list_t  *vtx_dir = sys_builder->vtx_dir;
  const cs_equation_param_t  *eqp = sys_builder->eqp;

  /* Test to remove */
  if (eqp->flag & CS_EQUATION_UNSTEADY)
    bft_error(__FILE__, __LINE__, 0,
              _(" Unsteady terms are not handled yet.\n"));

  /* Allocate and initialize a matrix with the larger stencil (that related
     to diffusion => all vertices of a cell are potentially in interaction) */
  cs_sla_matrix_t  *system_mat = _init_matrix(quant->n_vertices, connect);

  /* Compute contribution of (explicit) source terms to the rhs */
  _compute_st_rhs(m, connect, quant, tcur, sys_builder, rhs);

  /* Update builder with the current values of Dirchlet BCs */
  if (vtx_dir->n_nhmg_elts > 0) {

    cs_flag_t  dof_flag =
      CS_PARAM_FLAG_VERTEX | CS_PARAM_FLAG_PRIMAL | CS_PARAM_FLAG_SCAL;

    cs_cdo_bc_dirichlet_set(dof_flag,
                            tcur,
                            m,
                            eqp->bc,
                            vtx_dir,
                            sys_builder->dir_val);

  } // Dirichlet BCs */

  /* Build diffusion system: stiffness matrix */
  if (eqp->flag & CS_EQUATION_DIFFUSION) {

    /* Build the (full) stiffness matrix i.e. without taking into account BCs */
    _build_diffusion(connect, quant, sys_builder, system_mat);

    /* Treatment differs according to the enforcement of BCs
       In case of a strong enforcement of the BCs, rhs and matrix are resized */
    _update_diffusion_with_bc(m, connect, quant,
                              tcur,
                              sys_builder,
                              rhs,
                              &system_mat);

  }

  if (eqp->flag & CS_EQUATION_CONVECTION) {

    /* Only possible case: Weak enforcement of BCs */
    assert(sys_builder->n_dof_vertices == sys_builder->n_vertices);

    /* Build convection system (BCs are always weakly enforced) */
    _build_convection(m, connect, quant,
                      tcur,
                      sys_builder,
                      *rhs,
                      system_mat);

    /* Remove zero entries or entries which are very small w.r.t. to other
       entries of the same line */
    //    cs_sla_matrix_clean(system_mat, cs_get_eps_machine());

  }

  *sla_mat = system_mat;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Post-process the solution of a scalar convection/diffusion equation
 *         solved with a CDO vertex-based scheme.
 *
 * \param[in]      connect  pointer to a cs_cdo_connect_t struct.
 * \param[in]      quant    pointer to a cs_cdo_quantities_t struct.
 * \param[in]      solu     solution array
 * \param[in, out] builder  pointer to cs_cdovb_codits_t structure
 * \param[in, out] field    pointer to a cs_field_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdovb_codits_update_field(const cs_cdo_connect_t     *connect,
                             const cs_cdo_quantities_t  *quant,
                             const cs_real_t            *solu,
                             void                       *builder,
                             cs_field_t                 *field)
{
  int  i;

  cs_cdovb_codits_t  *_builder = (cs_cdovb_codits_t  *)builder;
  const cs_cdo_bc_list_t  *v_dir = _builder->vtx_dir;

  /* Set computed solution in field array */
  if (_builder->n_dof_vertices < _builder->n_vertices) {
    for (i = 0; i < _builder->n_vertices; i++)
      field->val[i] = 0;
    for (i = 0; i < _builder->n_dof_vertices; i++)
      field->val[_builder->v_z2i_ids[i]] = solu[i];
  }
  else
    for (i = 0; i < _builder->n_vertices; i++)
      field->val[i] = solu[i];

  /* Set BC in field array if we have this knowledge */
  if (_builder->enforce == CS_PARAM_BC_ENFORCE_STRONG)
    for (i = 0; i < v_dir->n_nhmg_elts; i++)
      field->val[v_dir->elt_ids[i]] = _builder->dir_val[i];
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
