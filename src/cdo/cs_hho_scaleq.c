/*============================================================================
 * Build an algebraic system for scalar conv./diff. eq. with Hybrid High Order
 * space discretization
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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_boundary_zone.h"
#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_diffusion.h"
#include "cs_equation_assemble.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_hho_builder.h"
#include "cs_hodge.h"
#include "cs_log.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_scheme_geometry.h"
#include "cs_search.h"
#include "cs_sdm.h"
#include "cs_source_term.h"

#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hho_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_HHO_SCALEQ_DBG     0

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/* Algebraic system for HHO discretization */

struct _cs_hho_scaleq_t {

  /* Ids related to the variable field and to the boundary flux field */
  int          var_field_id;
  int          bflux_field_id;

  /* System size (n_faces * n_face_dofs + n_cells * n_cell_dofs) */
  cs_lnum_t    n_dofs;
  int          n_max_loc_dofs;
  int          n_cell_dofs;
  int          n_face_dofs;

  /* Structures related to the algebraic sytem construction (shared) */
  const cs_matrix_structure_t    *ms;
  const cs_range_set_t           *rs;

  /* Solution of the algebraic system at the last computed iteration.
     cell_values is different from the values stored in the associated
     field since here it's the values of polynomial coefficients which
     are stored */
  cs_real_t                      *face_values;  /* DoF unknowns (x) + BCs */
  cs_real_t                      *cell_values;  /* DoF recomputed after the
                                                   static condensation */

  /* Storage of the source term (if one wants to apply a specific time
     discretization) */
  cs_real_t                      *source_terms;

  /* Handle the definition of the BCs */
  short int                      *bf2def_ids;

  /* Pointer of function to build the diffusion term */
  cs_cdo_enforce_bc_t            *enforce_dirichlet;

  /* Assembly process */
  /* ================ */

  cs_equation_assembly_t         *assemble;

  /* Static condensation members */
  /* =========================== */

  /* Acc_inv = Inverse of a diagonal matrix (block cell-cell)
     rc_tilda = Acc_inv * rhs_c (perform for each cell) */
  cs_real_t                      *rc_tilda;

  /* Lower-Left block of the full matrix (block cell-faces).
     Access to the values thanks to the c2f connectivity
     The transposed version of each block is stored for a more efficient
     usage */
  cs_sdm_t                       *acf_tilda;

};

/*============================================================================
 * Private variables
 *============================================================================*/

/* Size = 1 if openMP is not used */
static cs_cell_sys_t  **cs_hho_cell_sys = NULL;
static cs_cell_builder_t  **cs_hho_cell_bld = NULL;
static cs_hho_builder_t  **cs_hho_builders = NULL;

/* Pointer to shared structures (owned by a cs_domain_t structure) */
static const cs_cdo_quantities_t  *cs_shared_quant;
static const cs_cdo_connect_t  *cs_shared_connect;
static const cs_time_step_t  *cs_shared_time_step;
static const cs_matrix_structure_t  *cs_shared_ms0;
static const cs_matrix_structure_t  *cs_shared_ms1;
static const cs_matrix_structure_t  *cs_shared_ms2;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local builder structure used for building the system
 *          cellwise
 *
 * \param[in]  space_scheme   discretization scheme
 * \param[in]  connect        pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_cell_builder_t *
_cell_builder_create(cs_param_space_scheme_t     space_scheme,
                     const cs_cdo_connect_t     *connect)
{
  int  size;

  const int  n_fc = connect->n_max_fbyc;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  switch (space_scheme) {

  case CS_SPACE_SCHEME_HHO_P0:  /* TODO */
    {
      BFT_MALLOC(cb->ids, n_fc + 1, int);
      memset(cb->ids, 0, (n_fc + 1)*sizeof(int));

      /* For post-processing errors = 38 */
      size = CS_MAX(38, n_fc*(n_fc+1));
      BFT_MALLOC(cb->values, size, double);
      memset(cb->values, 0, size*sizeof(cs_real_t));

      size = CS_MAX(2*n_fc, 15);
      BFT_MALLOC(cb->vectors, size, cs_real_3_t);
      memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

      /* Local square dense matrices used during the construction of
         operators */
      cb->hdg = cs_sdm_square_create(n_fc);
      cb->loc = cs_sdm_square_create(n_fc + 1);
      cb->aux = cs_sdm_square_create(n_fc + 1);
    }
    break;

  case CS_SPACE_SCHEME_HHO_P1:
    {
      /* Store the block size description */
      size = n_fc + 1;
      BFT_MALLOC(cb->ids, size, int);
      memset(cb->ids, 0, size*sizeof(int));

      /* Store the face, cell and gradient basis function evaluations and
         the Gauss point weights
         --->  42 = 3 + 4 + 3*10 + 5
         or the factorization of the stiffness matrix used for computing
         the gradient reconstruction operator
         --->  n = 9 (10 - 1) --> facto = n*(n+1)/2 = 45
                              --> tmp buffer for facto = n = 9
                              --> 45 + 9 = 54
         or the factorization of the cell_cell block of size 4 --> 10
      */
      size = CS_MAX(54, (3*n_fc + 4)*2);
      BFT_MALLOC(cb->values, size, double);
      memset(cb->values, 0, size*sizeof(cs_real_t));

      /* Store Gauss points and tensor.n_f products */
      size = CS_MAX(15, 5 + n_fc);
      BFT_MALLOC(cb->vectors, size, cs_real_3_t);
      memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

      /* Local dense matrices used during the construction of operators */
      const int g_size = 9;                   /* basis (P_(k+1)) - 1 */
      for (int i = 0; i < n_fc; i++) cb->ids[i] = 3;
      cb->ids[n_fc] = 4;

      int  _sizes[3] = {1, 3, 6}; /* c0, cs-1, cs_kp1 - cs */
      cb->hdg = cs_sdm_block_create(1, 3, &g_size, _sizes);
      cb->loc = cs_sdm_block_create(n_fc + 1, n_fc + 1, cb->ids, cb->ids);
      cb->aux = cs_sdm_block_create(n_fc + 1, 1, cb->ids, &g_size); /* R_g^T */
    }
    break;

  case CS_SPACE_SCHEME_HHO_P2:
    {
      /* Store the block size description */
      size = n_fc + 1;
      BFT_MALLOC(cb->ids, size, int);
      memset(cb->ids, 0, size*sizeof(int));

      /* Store the face, cell and gradient basis function evaluations and
         the Gauss point weights */
      /* Store the face, cell and gradient basis function evaluations and
         the Gauss point weights
         --->  91 = 6 (f) + 10 (c) + 3*20 (g) + 15 (gauss)
         or the factorization of the stiffness matrix used for computing
         the gradient reconstruction operator
         --->  n = 19 (20 - 1) --> facto = n*(n+1)/2 = 190
                               --> tmp buffer for facto = n = 19
                              --> 190 + 19 = 209
      */
      size = CS_MAX(209, 2*(6*n_fc + 20));
      BFT_MALLOC(cb->values, size, double);
      memset(cb->values, 0, size*sizeof(cs_real_t));

      size = 15 + n_fc;  /* Store Gauss points and tensor.n_f products */
      BFT_MALLOC(cb->vectors, size, cs_real_3_t);
      memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

      /* Local dense matrices used during the construction of operators */
      const int g_size = 19; /* basis (P_(k+1)) - 1 */
      for (int i = 0; i < n_fc; i++) cb->ids[i] = 6;
      cb->ids[n_fc] = 10;

      int  _sizes[3] = {1, 9, 10}; /* c0, cs-1, cs_kp1 - cs */
      cb->hdg = cs_sdm_block_create(1, 3, &g_size, _sizes);
      cb->loc = cs_sdm_block_create(n_fc + 1, n_fc + 1, cb->ids, cb->ids);
      cb->aux = cs_sdm_block_create(n_fc + 1, 1, cb->ids, &g_size); /* R_g^T */
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _("Invalid space scheme."));

  } /* End of switch on space scheme */

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structures for the current cell
 *
 * \param[in]      cell_flag   flag related to the current cell
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         pointer to a cs_hho_scaleq_t structure
 * \param[in]      t_eval      time at which one performs the evaluation
 * \param[in, out] hhob        pointer to a cs_hho_builder_t structure
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_init_cell_system(const cs_flag_t               cell_flag,
                  const cs_cell_mesh_t         *cm,
                  const cs_equation_param_t    *eqp,
                  const cs_equation_builder_t  *eqb,
                  const cs_hho_scaleq_t        *eqc,
                  cs_real_t                     t_eval,
                  cs_hho_builder_t             *hhob,
                  cs_cell_sys_t                *csys,
                  cs_cell_builder_t            *cb)
{
  const int  n_dofs = eqc->n_cell_dofs + cm->n_fc*eqc->n_face_dofs;
  const int  n_blocks = cm->n_fc + 1;

  int  *block_sizes = cb->ids;
  for (int i = 0; i < cm->n_fc; i++)
    block_sizes[i] = eqc->n_face_dofs;
  block_sizes[cm->n_fc] = eqc->n_cell_dofs;

  csys->c_id = cm->c_id;
  csys->n_dofs = n_dofs;
  csys->cell_flag = cell_flag;

  /* Initialize the local system */
  cs_cell_sys_reset(cm->n_fc, csys);

  cs_sdm_block_init(csys->mat, n_blocks, n_blocks, block_sizes, block_sizes);

  for (int f = 0; f < cm->n_fc; f++) {

    const int _f_shift = eqc->n_face_dofs*f;
    const cs_lnum_t  f_shift = eqc->n_face_dofs*cm->f_ids[f];

    for (int i = 0; i < eqc->n_face_dofs; i++) {
      csys->dof_ids[_f_shift + i] = f_shift + i;
      csys->val_n[_f_shift + i] = eqc->face_values[f_shift + i];
    }

  }

  const int _c_shift = eqc->n_face_dofs*cm->n_fc;
  const cs_lnum_t  c_shift = eqc->n_cell_dofs*cm->c_id;

  for (int i = 0; i < eqc->n_cell_dofs; i++) {
    csys->dof_ids[_c_shift + i] = c_shift + i;
    csys->val_n[_c_shift + i] = eqc->cell_values[c_shift + i];
  }

  /* Store the local values attached to Dirichlet values if the current cell
     has at least one border face */
  if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

    /* Identify which face is a boundary face */
    for (short int f = 0; f < cm->n_fc; f++) {

      const cs_lnum_t  bf_id = cm->f_ids[f] - cm->bface_shift;
      if (bf_id > -1) {        /*  Border face */

        const cs_flag_t  face_flag = eqb->face_bc->flag[bf_id];
        const int  _f_shift = eqc->n_face_dofs*f;

        csys->bf_flag[csys->n_bc_faces] = face_flag;
        csys->_f_ids[csys->n_bc_faces] = f;
        csys->bf_ids[csys->n_bc_faces] = bf_id;
        csys->n_bc_faces++;

        if (face_flag & CS_CDO_BC_HMG_DIRICHLET) {

          csys->has_dirichlet = true;
          for (int k = 0; k < eqc->n_face_dofs; k++)
            csys->dof_flag[_f_shift + k] |= CS_CDO_BC_HMG_DIRICHLET;

        }
        else if (face_flag & CS_CDO_BC_DIRICHLET) {

          csys->has_dirichlet = true;
          for (int k = 0; k < eqc->n_face_dofs; k++)
            csys->dof_flag[_f_shift + k] |= CS_CDO_BC_DIRICHLET;

          /* Compute the value of the Dirichlet BC */
          const short int  def_id = eqc->bf2def_ids[bf_id];
          assert(def_id > -1);

          cs_real_t  *dir_reduction = csys->dir_values + f*eqc->n_face_dofs;
          cs_hho_builder_compute_dirichlet(eqp->bc_defs[def_id],
                                           f,
                                           cm,
                                           t_eval,
                                           cb,
                                           hhob,
                                           dir_reduction);

        }
        else { /* TODO */
          bft_error(__FILE__, __LINE__, 0, "%s: TODO", __func__);
        }

      } /* Border face id */

    } /* Loop on cell faces */

#if defined(DEBUG) && !defined(NDEBUG) /* Sanity check */
    for (short int f = 0; f < eqc->n_face_dofs*cm->n_fc; f++) {
      if (csys->dof_flag[f] & CS_CDO_BC_HMG_DIRICHLET)
        if (fabs(csys->dir_values[f]) > 10*DBL_MIN)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: Invalid enforcement of Dirichlet BCs on faces",
                    __func__);
    }
#endif

  } /* This is a border cell */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Proceed to a static condensation of the local system and keep
 *          information inside the builder to be able to compute the values
 *          at cell centers
 *
 * \param[in]      c2f       pointer to a cs_adjacency_t structure
 * \param[in, out] eqc       pointer to a cs_hho_scaleq_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 * \param[in, out] csys      pointer to a cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_condense_and_store(const cs_adjacency_t    *c2f,
                    cs_hho_scaleq_t         *eqc,
                    cs_cell_builder_t       *cb,
                    cs_cell_sys_t           *csys)
{
  /* m is the matrix on which the static condensation will be performed. m
   * is divided in blocks as follows:
   *                                          _
   *              ||                |      ||  |
   *              ||                |      ||  |
   *              ||       m_FF     | m_FC ||  | #faces * face_size
   *       m  :=  ||                |      ||  |
   *              ||________________|______|| _|
   *              ||                |      ||  |
   *              ||       m_CF     | m_CC ||  | cell_size
   *              ||                |      || _|
   *
   *               |________________|______|
   *                #faces*face_size cell_size
   *
   */

  cs_sdm_t  *m = csys->mat;
  cs_sdm_block_t  *bd = m->block_desc;

  const int  n_fc = bd->n_row_blocks - 1;
  const int  _f_offset = eqc->n_face_dofs*n_fc;
  const cs_lnum_t  c_offset = eqc->n_cell_dofs*csys->c_id;

  /* Store information to compute the cell values after the resolution */
  const double  *_cell_rhs = csys->rhs + _f_offset;

  /* mCC is a small square matrix:
     Compute its modified Cholesky factorization */
  const cs_sdm_t  *mCC = cs_sdm_get_block(m, n_fc, n_fc);

  /* Store the transposed of m_CF to speed-up the computation of aCF_tilda */
  cs_real_t  *acc_facto = cb->values;
  cs_lnum_t  c2f_shift = c2f->idx[csys->c_id];

  switch (eqc->n_cell_dofs) {
  case 1:                       /* HHO k=0, factor. size = 1 */
    acc_facto[0] = 1./mCC->val[0];
    eqc->rc_tilda[c_offset] = acc_facto[0] * _cell_rhs[0];
    for (int f = 0; f < n_fc; f++) {
      const cs_sdm_t  *mCF = cs_sdm_get_block(m, n_fc, f);
      cs_sdm_t  *_acf = cs_sdm_get_block(eqc->acf_tilda, c2f_shift + f, 0);
      _acf->val[0] = mCF->val[0] * acc_facto[0];
    }
    break;

  case 4:                       /* HHO k=1, facto. size = 10 */
    {
      assert(eqc->n_face_dofs == 3);
      cs_real_t  _rhs[4];

      cs_sdm_44_ldlt_compute(mCC, acc_facto);
      cs_sdm_44_ldlt_solve(acc_facto, _cell_rhs, eqc->rc_tilda + c_offset);

      for (int f = 0; f < n_fc; f++) {
        const cs_sdm_t  *mCF = cs_sdm_get_block(m, n_fc, f);
        cs_sdm_t  *_acf = cs_sdm_get_block(eqc->acf_tilda, c2f_shift + f, 0);
        for (int f_dof = 0; f_dof < 3; f_dof++) {
          for (int c_dof = 0; c_dof < 4; c_dof++)
            _rhs[c_dof] = mCF->val[3*c_dof + f_dof];
          cs_sdm_44_ldlt_solve(acc_facto, _rhs, _acf->val + 4*f_dof);
        }
      }
    }
    break;

  case 10:                      /* HHO k=2, facto. size = 55 + 10 as tmp */
    {
      assert(eqc->n_face_dofs == 6);
      cs_real_t  _rhs[10];

      cs_sdm_ldlt_compute(mCC, acc_facto, acc_facto + 55);
      cs_sdm_ldlt_solve(10, acc_facto, _cell_rhs, eqc->rc_tilda + c_offset);

      for (int f = 0; f < n_fc; f++) {
        const cs_sdm_t  *mCF = cs_sdm_get_block(m, n_fc, f);
        cs_sdm_t  *_acf = cs_sdm_get_block(eqc->acf_tilda, c2f_shift + f, 0);
        for (int f_dof = 0; f_dof < 6; f_dof++) {
          for (int c_dof = 0; c_dof < 10; c_dof++)
            _rhs[c_dof] = mCF->val[6*c_dof + f_dof];
          cs_sdm_ldlt_solve(10, acc_facto, _rhs, _acf->val + 10*f_dof);
        }
      }
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid case.", __func__);
    break;
  }

  /* Update csys:
     Compute m_FF = m_FF - m_FC * m_CC_inv * m_CF */
  cs_real_t  *bf_tilda = cb->values;
  cs_sdm_t  *_aff = cb->aux;

  for (int fi = 0; fi < n_fc; fi++) {

    /* Initial block to update */
    cs_sdm_t  *m_fc = cs_sdm_get_block(m, fi, n_fc);

    cs_sdm_matvec(m_fc, eqc->rc_tilda + c_offset, bf_tilda);

    /* Update RHS: RHS_f = RHS_f - Afc*Acc^-1*s_c */
    for (int k = 0; k < eqc->n_face_dofs; k++)
      csys->rhs[eqc->n_face_dofs*fi + k] -= bf_tilda[k];

    for (int fj = 0; fj < n_fc; fj++) {

      cs_sdm_t  *mFF = cs_sdm_get_block(m, fi, fj);
      cs_sdm_t  *_acf = cs_sdm_get_block(eqc->acf_tilda, c2f_shift + fj, 0);

      cs_sdm_init(eqc->n_face_dofs, eqc->n_face_dofs, _aff);
      cs_sdm_multiply_rowrow(m_fc, _acf, _aff);
      cs_sdm_add_mult(mFF, -1, _aff);

    } /* fj */
  } /* fi */

  /* Reshape matrix */
  int  shift = n_fc;
  for (short int bfi = 1; bfi < n_fc; bfi++) {
    for (short int bfj = 0; bfj < n_fc; bfj++) {

      cs_sdm_t  *mFF_old = cs_sdm_get_block(m, bfi, bfj);

      /* Set the block (i,j) */
      cs_sdm_t  *mFF = bd->blocks + shift;
      cs_sdm_map_array(eqc->n_face_dofs, eqc->n_face_dofs, mFF, mFF_old->val);
      shift++;

    }
  }

  csys->n_dofs = _f_offset;
  m->n_rows = m->n_cols = _f_offset;
  bd->n_row_blocks = n_fc;      /* instead of n_fc + 1 */
  bd->n_col_blocks = n_fc;      /* instead of n_fc + 1 */
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate work buffer and general structures related to HHO schemes
 *         Set shared pointers
 *
 * \param[in]  scheme_flag  flag to identify which kind of numerical scheme is
 *                          requested to solve the computational domain
 * \param[in]  quant        additional mesh quantities struct.
 * \param[in]  connect      pointer to a cs_cdo_connect_t struct.
 * \param[in]  time_step    pointer to a time step structure
 * \param[in]  ms0          pointer to a cs_matrix_structure_t structure (P0)
 * \param[in]  ms1          pointer to a cs_matrix_structure_t structure (P1)
 * \param[in]  ms2          pointer to a cs_matrix_structure_t structure (P2)
*/
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_init_common(cs_flag_t                      scheme_flag,
                          const cs_cdo_quantities_t     *quant,
                          const cs_cdo_connect_t        *connect,
                          const cs_time_step_t          *time_step,
                          const cs_matrix_structure_t   *ms0,
                          const cs_matrix_structure_t   *ms1,
                          const cs_matrix_structure_t   *ms2)
{
  /* Assign static const pointers */
  cs_shared_quant = quant;
  cs_shared_connect = connect;
  cs_shared_time_step = time_step;
  cs_shared_ms0 = ms0;
  cs_shared_ms1 = ms1;
  cs_shared_ms2 = ms2;

  const int n_fc = connect->n_max_fbyc;

  int  order = -1, fbs = 0, cbs = 0;
  cs_param_space_scheme_t  space_scheme;

  if (scheme_flag & CS_FLAG_SCHEME_POLY2) {

    space_scheme = CS_SPACE_SCHEME_HHO_P2;
    fbs = CS_N_FACE_DOFS_2ND; // DoF by face
    cbs = CS_N_CELL_DOFS_2ND; // DoF for the cell
    order = 2;

  }
  else if (scheme_flag & CS_FLAG_SCHEME_POLY1) {

    space_scheme = CS_SPACE_SCHEME_HHO_P1;
    fbs = CS_N_FACE_DOFS_1ST; // DoF by face
    cbs = CS_N_CELL_DOFS_1ST;  // DoF for the cell
    order = 1;

  }
  else {

    space_scheme = CS_SPACE_SCHEME_HHO_P0;
    fbs = CS_N_FACE_DOFS_0TH; // DoF by face
    cbs = CS_N_CELL_DOFS_0TH; // DoF for the cell
    order = 0;

  }

  const int n_dofs = n_fc * fbs + cbs;

    /* Structure used to build the final system by a cell-wise process */
  assert(cs_glob_n_threads > 0);  /* Sanity check */

  BFT_MALLOC(cs_hho_cell_bld, cs_glob_n_threads, cs_cell_builder_t *);
  BFT_MALLOC(cs_hho_cell_sys, cs_glob_n_threads, cs_cell_sys_t *);

  /* Allocate builder structure specific to HHO schemes. This is an additional
     structure with respect to a cs_cell_builder_t structure */
  BFT_MALLOC(cs_hho_builders, cs_glob_n_threads, cs_hho_builder_t *);

  for (int i = 0; i < cs_glob_n_threads; i++) {
    cs_hho_cell_bld[i] = NULL;
    cs_hho_cell_sys[i] = NULL;
    cs_hho_builders[i] = NULL;
  }

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    cs_cell_builder_t  *cb = _cell_builder_create(space_scheme, connect);
    cs_hho_cell_bld[t_id] = cb;
    cs_hho_builders[t_id] = cs_hho_builder_create(order, n_fc);

    for (int i = 0; i < n_fc; i++) cb->ids[i] = fbs;
    cb->ids[n_fc] = cbs;

    cs_hho_cell_sys[t_id] = cs_cell_sys_create(n_dofs,
                                               fbs*n_fc,
                                               n_fc + 1, cb->ids);
  }
#else
  assert(cs_glob_n_threads == 1);

    cs_cell_builder_t  *cb = _cell_builder_create(space_scheme, connect);
    cs_hho_cell_bld[0] = cb;
    cs_hho_builders[0] = cs_hho_builder_create(order, n_fc);

    for (int i = 0; i < n_fc; i++) cb->ids[i] = fbs;
    cb->ids[n_fc] = cbs;

    cs_hho_cell_sys[0] = cs_cell_sys_create(n_dofs,
                                               fbs*n_fc,
                                               n_fc + 1, cb->ids);
#endif /* openMP */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys    pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb      pointer to a pointer on a cs_cell_builder_t structure
 * \param[out]  hhob    pointer to a pointer on a cs_hho_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_get(cs_cell_sys_t       **csys,
                  cs_cell_builder_t   **cb,
                  cs_hho_builder_t    **hhob)
{
  int t_id = 0;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
  t_id = omp_get_thread_num();
  assert(t_id < cs_glob_n_threads);
#endif /* openMP */

  *csys = cs_hho_cell_sys[t_id];
  *cb = cs_hho_cell_bld[t_id];
  *hhob = cs_hho_builders[t_id];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free buffers and generic structures related to HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_finalize_common(void)
{
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();

    cs_cell_sys_free(&(cs_hho_cell_sys[t_id]));
    cs_cell_builder_free(&(cs_hho_cell_bld[t_id]));
    cs_hho_builder_free(&(cs_hho_builders[t_id]));
  }
#else
  assert(cs_glob_n_threads == 1);

  cs_cell_sys_free(&(cs_hho_cell_sys[0]));
  cs_cell_builder_free(&(cs_hho_cell_bld[0]));
  cs_hho_builder_free(&(cs_hho_builders[0]));
#endif /* openMP */

  BFT_FREE(cs_hho_cell_sys);
  BFT_FREE(cs_hho_cell_bld);
  BFT_FREE(cs_hho_builders);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_hho_scaleq_t structure storing data useful for
 *         building and managing such a scheme
 *
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id     id of the variable field
 * \param[in]      bflux_id   id of the boundary flux field
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated cs_hho_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void  *
cs_hho_scaleq_init_context(const cs_equation_param_t   *eqp,
                           int                          var_id,
                           int                          bflux_id,
                           cs_equation_builder_t       *eqb)
{
  /* Sanity checks */
  assert(eqp != NULL);
  if (eqp->dim != 1)
    bft_error(__FILE__, __LINE__, 0, " Expected: scalar-valued HHO equation.");

  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_lnum_t  n_faces = connect->n_faces[0];
  const cs_lnum_t  n_cells = connect->n_cells;

  cs_hho_scaleq_t  *eqc = NULL;

  BFT_MALLOC(eqc, 1, cs_hho_scaleq_t);

  eqc->var_field_id = var_id;
  eqc->bflux_field_id = bflux_id;

  /* Mesh flag to know what to build */
  eqb->msh_flag = CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_PFQ |
    CS_FLAG_COMP_FE | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ |
    CS_FLAG_COMP_EV | CS_FLAG_COMP_DIAM;

  switch (eqp->space_scheme) {

  case CS_SPACE_SCHEME_HHO_P0:
    eqc->n_cell_dofs = CS_N_CELL_DOFS_0TH;
    eqc->n_face_dofs = CS_N_FACE_DOFS_0TH;

    /* Not owner; Only shared */
    eqc->ms = cs_shared_ms0;
    eqc->rs = connect->range_sets[CS_CDO_CONNECT_FACE_SP0];

    /* Assembly process */
    eqc->assemble = cs_equation_assemble_set(CS_SPACE_SCHEME_HHO_P0,
                                             CS_CDO_CONNECT_FACE_SP0);
    break;

  case CS_SPACE_SCHEME_HHO_P1:
    eqc->n_cell_dofs = CS_N_CELL_DOFS_1ST;
    eqc->n_face_dofs = CS_N_FACE_DOFS_1ST;

    /* Not owner; Only shared */
    eqc->ms = cs_shared_ms1;
    eqc->rs = connect->range_sets[CS_CDO_CONNECT_FACE_SP1];

    /* Assembly process */
    eqc->assemble = cs_equation_assemble_set(CS_SPACE_SCHEME_HHO_P1,
                                             CS_CDO_CONNECT_FACE_SP1);
    break;

  case CS_SPACE_SCHEME_HHO_P2:
    eqc->n_cell_dofs = CS_N_CELL_DOFS_2ND;
    eqc->n_face_dofs = CS_N_FACE_DOFS_2ND;

    /* Not owner; Only shared */
    eqc->ms = cs_shared_ms2;
    eqc->rs = connect->range_sets[CS_CDO_CONNECT_FACE_SP2];

    /* Assembly process */
    eqc->assemble = cs_equation_assemble_set(CS_SPACE_SCHEME_HHO_P2,
                                             CS_CDO_CONNECT_FACE_SP2);
    break;

    /* TODO: case CS_SPACE_SCHEME_HHO_PK */
  default:
    bft_error(__FILE__, __LINE__, 0, " %s: Invalid space scheme.", __func__);

  }

  /* System dimension */
  eqc->n_dofs = eqc->n_face_dofs * n_faces;
  eqc->n_max_loc_dofs = eqc->n_face_dofs*connect->n_max_fbyc + eqc->n_cell_dofs;

  /* Values of each DoF related to the cells */
  const cs_lnum_t  n_cell_dofs = n_cells * eqc->n_cell_dofs;
  BFT_MALLOC(eqc->cell_values, n_cell_dofs, cs_real_t);
  memset(eqc->cell_values, 0, sizeof(cs_real_t)*n_cell_dofs);

  /* Values at each face (interior and border) i.e. take into account BCs */
  BFT_MALLOC(eqc->face_values, eqc->n_dofs, cs_real_t);
  memset(eqc->face_values, 0, sizeof(cs_real_t)*eqc->n_dofs);

  /* Source term */
  eqc->source_terms = NULL;
  if (cs_equation_param_has_sourceterm(eqp)) {

    BFT_MALLOC(eqc->source_terms, n_cell_dofs, cs_real_t);
    memset(eqc->source_terms, 0, sizeof(cs_real_t)*n_cell_dofs);

  } /* There is at least one source term */

  /* Members related to the static condensation.
     The transposed of acf_tilda is stored to speed-up the computation of
     the static condensation */
  BFT_MALLOC(eqc->rc_tilda, n_cell_dofs, cs_real_t);
  memset(eqc->rc_tilda, 0, sizeof(cs_real_t)*n_cell_dofs);

  cs_lnum_t  n_row_blocks = connect->c2f->idx[n_cells];
  int  *row_block_sizes = NULL;

  BFT_MALLOC(row_block_sizes, n_row_blocks, int);
# pragma omp parallel for if (n_cells > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_row_blocks; i++)
    row_block_sizes[i] = eqc->n_face_dofs;

  int  col_block_size = eqc->n_cell_dofs;
  eqc->acf_tilda = cs_sdm_block_create(n_row_blocks, 1,
                                       row_block_sizes, &col_block_size);
  cs_sdm_block_init(eqc->acf_tilda,
                    n_row_blocks, 1,
                    row_block_sizes, &col_block_size);

  BFT_FREE(row_block_sizes);

  /* Handle boundary conditions */
  const cs_lnum_t  n_b_faces = connect->n_faces[1];
  BFT_MALLOC(eqc->bf2def_ids, n_b_faces, short int);

# pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_b_faces; i++)
    eqc->bf2def_ids[i] = -1; /* Default BC has no definition --> -1 */

  for (int def_id = 0; def_id < eqp->n_bc_defs; def_id++) {

    const cs_xdef_t  *def = eqp->bc_defs[def_id];
    const cs_zone_t  *bz = cs_boundary_zone_by_id(def->z_id);

#   pragma omp parallel for if (bz->n_elts > CS_THR_MIN)
    for (cs_lnum_t i = 0; i < bz->n_elts; i++)
      eqc->bf2def_ids[bz->elt_ids[i]] = def_id;

  } /* Loop on BC definitions */

  /* Diffusion */
  eqc->enforce_dirichlet = NULL;

  if (cs_equation_param_has_diffusion(eqp)) {

    switch (eqp->default_enforcement) {

    case CS_PARAM_BC_ENFORCE_ALGEBRAIC:
      eqc->enforce_dirichlet = cs_cdo_diffusion_alge_block_dirichlet;
      break;

    case CS_PARAM_BC_ENFORCE_PENALIZED:
      eqc->enforce_dirichlet = cs_cdo_diffusion_pena_block_dirichlet;
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                (" %s Invalid type of algorithm to enforce Dirichlet BC."),
                __func__);

    }

  } /* Has diffusion term to handle */

  return eqc;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a cs_hho_scaleq_t structure
 *
 * \param[in, out]  data    pointer to a cs_hho_scaleq_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_hho_scaleq_free_context(void   *data)
{
  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t *)data;

  if (eqc == NULL)
    return eqc;

  /* Free temporary buffers */
  BFT_FREE(eqc->cell_values);
  BFT_FREE(eqc->face_values);
  BFT_FREE(eqc->rc_tilda);
  BFT_FREE(eqc->source_terms);
  BFT_FREE(eqc->bf2def_ids);

  cs_sdm_free(eqc->acf_tilda);

  /* Last free */
  BFT_FREE(eqc);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued HHO schemes.
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
cs_hho_scaleq_init_values(cs_real_t                     t_eval,
                          const int                     field_id,
                          const cs_mesh_t              *mesh,
                          const cs_equation_param_t    *eqp,
                          cs_equation_builder_t        *eqb,
                          void                         *context)
{
  /* Unused parameters --> generic function pointer */
  CS_UNUSED(field_id);
  CS_UNUSED(eqb);
  CS_UNUSED(t_eval);
  CS_UNUSED(mesh);
  CS_UNUSED(eqp);

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t *)context;
  cs_real_t  *f_vals = eqc->face_values;
  cs_real_t  *c_vals = eqc->cell_values;

  memset(f_vals, 0, quant->n_faces * eqc->n_face_dofs * sizeof(cs_real_t));
  memset(c_vals, 0, quant->n_cells * eqc->n_cell_dofs * sizeof(cs_real_t));

  /* TODO */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         data structure
 *
 * \param[in]      eqp            pointer to a cs_equation_param_t structure
 * \param[in, out] eqb            pointer to a cs_equation_builder_t structure
 * \param[in, out] data           pointer to generic data structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_initialize_system(const cs_equation_param_t  *eqp,
                                cs_equation_builder_t      *eqb,
                                void                       *data,
                                cs_matrix_t               **system_matrix,
                                cs_real_t                 **system_rhs)
{
  CS_UNUSED(eqp);

  assert(*system_matrix == NULL && *system_rhs == NULL);

  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t *)data;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;

  cs_timer_t  t0 = cs_timer_time();
  const cs_lnum_t  n_elts = quant->n_faces * eqc->n_face_dofs;

  *system_matrix = cs_matrix_create(eqc->ms);

  /* Allocate and initialize the related right-hand side */
  BFT_MALLOC(*system_rhs, n_elts, cs_real_t);
# pragma omp parallel for if  (n_elts > CS_THR_MIN)
  for (cs_lnum_t i = 0; i < n_elts; i++) (*system_rhs)[i] = 0.0;

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build the linear system arising from a scalar convection/diffusion
 *         equation with a HHO scheme.
 *         One works cellwise and then process to the assembly
 *
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_hho_scaleq_t structure
 * \param[in, out] rhs        right-hand side
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_build_system(const cs_mesh_t            *mesh,
                           const cs_real_t            *field_val,
                           const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *data,
                           cs_real_t                  *rhs,
                           cs_matrix_t                *matrix)
{
  CS_UNUSED(mesh);
  CS_UNUSED(field_val);

  /* Sanity checks */
  assert(rhs != NULL && matrix != NULL && eqp != NULL && eqb != NULL);
  /* The only way to set a Dirichlet up to now */
  assert(eqp->default_enforcement == CS_PARAM_BC_ENFORCE_PENALIZED ||
         eqp->default_enforcement == CS_PARAM_BC_ENFORCE_ALGEBRAIC);

  /* Test to remove */
  if (cs_equation_param_has_convection(eqp))
    bft_error(__FILE__, __LINE__, 0,
              _(" Convection term is not handled yet.\n"));
  if (cs_equation_param_has_time(eqp))
    bft_error(__FILE__, __LINE__, 0,
              _(" Unsteady terms are not handled yet.\n"));

  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t *)data;

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;
  const cs_time_step_t  *ts = cs_shared_time_step;
  const cs_real_t  t_cur = ts->t_cur;
  const cs_real_t  dt_cur = ts->dt[0];

  cs_timer_t  t0 = cs_timer_time();

  /* Initialize the structure to assemble values */
  cs_matrix_assembler_values_t  *mav
    = cs_matrix_assembler_values_init(matrix, NULL, NULL);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, eqp, eqb, eqc, rhs, matrix, mav,               \
         field_val, cs_hho_cell_sys, cs_hho_cell_bld, cs_hho_builders)  \
  firstprivate(dt_cur, t_cur)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_real_t  t_eval_pty = t_cur + 0.5*dt_cur;

    /* Set inside the OMP section so that each thread has its own value
     * Each thread get back its related structures:
     * Get the cell-wise view of the mesh and the algebraic system */
    cs_equation_assemble_t  *eqa = cs_equation_assemble_get(t_id);
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_sys_t  *csys = cs_hho_cell_sys[t_id];
    cs_cell_builder_t  *cb = cs_hho_cell_bld[t_id];
    cs_hho_builder_t  *hhob = cs_hho_builders[t_id];

    /* Initialization of the values of properties */
    cs_equation_init_properties(eqp, eqb, t_eval_pty, cb);

    /* --------------------------------------------- */
    /* Main loop on cells to build the linear system */
    /* --------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag = cs_equation_cell_mesh_flag(cell_flag, eqb);

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Define set of basis functions for cell faces and the current cell */
      cs_hho_builder_cellwise_setup(cm, cb, hhob);

      /* Set the local (i.e. cellwise) structures for the current cell */
      _init_cell_system(cell_flag, cm, eqp, eqb, eqc,
                        t_cur + dt_cur, /* For Dirichlet up to now */
                        hhob, csys, cb);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 2
      if (cs_dbg_cw_test(eqp, cm, csys)) cs_cell_mesh_dump(cm);
#endif

      const short int  face_offset = cm->n_fc*eqc->n_face_dofs;

      /* DIFFUSION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ============================================== */

      if (cs_equation_param_has_diffusion(eqp)) {

        /* Define the local stiffness matrix */
        if (!(eqb->diff_pty_uniform))
          cs_equation_set_diffusion_property_cw(eqp, cm, t_eval_pty, cell_flag,
                                                cb);

        cs_hho_builder_compute_grad_reco(cm, cb, hhob);

        /* Local matrix owned by the cellwise builder (store in cb->loc) */
        cs_hho_builder_diffusion(cm, cb, hhob);

        /* Add the local diffusion operator to the local system */
        cs_sdm_block_add(csys->mat, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 1
        if (cs_dbg_cw_test(eqp, cm, csys))
          cs_cell_sys_dump("\n>> Local system after diffusion", csys);
#endif
      } /* END OF DIFFUSION */

      /* SOURCE TERM COMPUTATION */
      /* ======================= */

      if (cs_equation_param_has_sourceterm(eqp)) {

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
                                        t_eval_pty,
                                        hhob,   /* input structure */
                                        cb,     /* mass matrix is cb->hdg */
                                        csys->source);

        if (cs_equation_param_has_time(eqp) == false) { /* Steady-case */

          /* Same strategy as if one applies a implicit scheme */
          cs_real_t  *_rhs = csys->rhs + face_offset;
          const cs_real_t  *_st = csys->source + face_offset;
          for (int i = 0; i < eqc->n_cell_dofs; i++)
            _rhs[i] += _st[i];

        }

        /* Reset the value of the source term for the cell DoF
           Source term is only hold by the cell DoF in face-based schemes */
        {
          cs_real_t  *st = eqc->source_terms + c_id * eqc->n_cell_dofs;
          const cs_real_t  *_st = csys->source + face_offset;
          for (int i = 0; i < eqc->n_cell_dofs; i++)
            st[i] = _st[i];
        }

      } /* End of term source contribution */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 1
      if (cs_dbg_cw_test(eqp, cm, csys))
        cs_cell_sys_dump(">> Local system matrix before condensation", csys);
#endif

      /* Static condensation of the local system matrix of n_fc + 1 blocks into
         a matrix of n_fc block. Store information in the input structure in
         order to be able to compute the values at cell centers. */
      _condense_and_store(connect->c2f, eqc, cb, csys);

      /* BOUNDARY CONDITION CONTRIBUTION TO THE ALGEBRAIC SYSTEM */
      /* ======================================================= */

      /* TODO: Neumann boundary conditions */

      if (cell_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {

        if (cs_equation_param_has_diffusion(eqp)) {

          /* Weakly enforced Dirichlet BCs for cells attached to the boundary
             csys is updated inside (matrix and rhs)
             eqp->diffusion_hodge is a dummy parameter (not used)
          */
          eqc->enforce_dirichlet(eqp, cm, NULL, cb, csys);

        } /* diffusion term */

      }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 0
      if (cs_dbg_cw_test(eqp, cm, csys)) {
        cs_cell_sys_dump(">> (FINAL) Local system matrix", csys);
        cs_sdm_block_fprintf(NULL, NULL, 1e-16, csys->mat);
      }
#endif

      /* ASSEMBLY */
      /* ======== */

      /* Matrix assembly */
      eqc->assemble(csys, eqc->rs, eqa, mav);

      /* RHS assembly */
      for (short int i = 0; i < eqc->n_face_dofs*cm->n_fc; i++) {
#       pragma omp atomic
        rhs[csys->dof_ids[i]] += csys->rhs[i];
      }

    } /* Main loop on cells */

  } /* OPENMP Block */

  cs_matrix_assembler_values_done(mav); // optional

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 2
  cs_dbg_darray_to_listing("FINAL RHS_FACE",
                           eqc->n_dofs, rhs, eqc->n_face_dofs);
  if (eqc->source_terms != NULL)
    cs_dbg_darray_to_listing("FINAL RHS_CELL", quant->n_cells,
                             eqc->source_terms, eqc->n_cell_dofs);
#endif

  /* Free temporary buffers and structures */
  cs_matrix_assembler_values_finalize(&mav);

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tcb), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values required for hybrid discretization
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to cs_hho_scaleq_t structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_update_field(const cs_real_t            *solu,
                           const cs_real_t            *rhs,
                           const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *data,
                           cs_real_t                  *field_val)
{
  CS_UNUSED(rhs);
  CS_UNUSED(eqp);

  cs_timer_t  t0 = cs_timer_time();

  const cs_cdo_quantities_t  *quant = cs_shared_quant;
  const cs_cdo_connect_t  *connect = cs_shared_connect;

  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t  *)data;

  /* Reset the solution at cells */
  memset(eqc->cell_values, 0,
         sizeof(cs_real_t) * eqc->n_cell_dofs * quant->n_cells);

# pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)    \
  shared(quant, connect, eqp, eqb, eqc, rhs, solu, field_val,           \
         cs_hho_cell_bld, cs_hho_builders)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    /* Each thread get back its related structures:
       Get the cell-wise view of the mesh and the algebraic system */
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = cs_hho_cell_bld[t_id];
    cs_hho_builder_t  *hhob = cs_hho_builders[t_id];

    /* Set inside the OMP section so that each thread has its own value */

    /* ----------------------------------------------------------- */
    /* Main loop on cells to reconstruct the field at cell centers */
    /* ----------------------------------------------------------- */

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      const cs_lnum_t  c2f_shift = connect->c2f->idx[c_id];
      const cs_flag_t  cell_flag = connect->cell_flag[c_id];
      const cs_flag_t  msh_flag = cs_equation_cell_mesh_flag(cell_flag, eqb);

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Define set of basis functions for cell faces and the current cell */
      cs_hho_builder_cellbasis_setup(cm, cb, hhob);

      /* Evaluate the solution at cell centers */
      cs_real_t  *f_contrib = cb->values + eqc->n_cell_dofs;
      cs_real_t  *c_vals = eqc->cell_values + c_id * eqc->n_cell_dofs;

      memset(f_contrib, 0, sizeof(cs_real_t)*eqc->n_cell_dofs);

      /* Recover the cell DoFs
         x_C = A_CC^-1 b_C (--> rc_tilda) - acf_tilda * x_F  */
      for (short int f = 0; f < cm->n_fc; f++) {

        cs_sdm_t  *_acf = cs_sdm_get_block(eqc->acf_tilda, c2f_shift + f, 0);
        const cs_real_t  *f_vals = solu + cm->f_ids[f] * eqc->n_face_dofs;

        /* Update c_vals c_vals += _acf * f_vals */
        cs_sdm_matvec_transposed(_acf, f_vals, f_contrib);

      }

      const cs_real_t  *_rc = eqc->rc_tilda + c_id * eqc->n_cell_dofs;
      cs_real_t  *cphi_eval = cb->values;

      hhob->cell_basis->eval_all_at_point(hhob->cell_basis, cm->xc, cphi_eval);

      field_val[c_id] = 0.;
      for (short int i = 0; i < eqc->n_cell_dofs; i++) {
        c_vals[i] = _rc[i] - f_contrib[i];
        field_val[c_id] += cphi_eval[i] * c_vals[i];
      }

    } // Main loop on cells

  } // OPENMP Block

  /* Set the computed solution in field array */
  memcpy(eqc->face_values, solu, sizeof(cs_real_t)*eqc->n_dofs);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HHO_SCALEQ_DBG > 2
  cs_dbg_darray_to_listing("FINAL FACE_VALUES",
                           eqc->n_dofs, eqc->face_values, eqc->n_face_dofs);
  cs_dbg_darray_to_listing("FINAL CELL_VALUES", quant->n_cells,
                             eqc->cell_values, eqc->n_cell_dofs);
  cs_dbg_darray_to_listing("FINAL CELL_CENTER_VALUES", quant->n_cells,
                             field_val, 1);
#endif

  cs_timer_t  t1 = cs_timer_time();
  cs_timer_counter_add_diff(&(eqb->tce), &t0, &t1);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at faces (DoF used in the linear system are
 *         located at primal faces)
 *
 * \param[in, out]  data    pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_hho_scaleq_get_face_values(void          *data)
{
  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t  *)data;

  if (eqc == NULL)
    return NULL;
  else
    return eqc->face_values;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at cells (DoF used in the linear system are
 *         located at primal faces)
 *
 * \param[in, out]  data    pointer to a data structure cast on the fly
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_hho_scaleq_get_cell_values(void          *data)
{
  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t  *)data;

  if (eqc == NULL)
    return NULL;
  else
    return eqc->cell_values;
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
cs_hho_scaleq_read_restart(cs_restart_t    *restart,
                           const char      *eqname,
                           void            *scheme_context)
{
  /* Only the face values are handled. Cell values are stored in a cs_field_t
     structure and thus are handled automatically. */
  if (restart == NULL)
    return;
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is NULL", __func__);
  if (scheme_context == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Scheme context is NULL", __func__);

  int retcode = CS_RESTART_SUCCESS;
  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t  *)scheme_context;

  char sec_name[128];

  /* Handle interior faces */
  /* ===================== */

  const int  i_ml_id = cs_mesh_location_get_id_by_name(N_("interior_faces"));

  /* Define the section name */
  snprintf(sec_name, 127, "%s::i_face_vals", eqname);

  /* Check section */
  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     i_ml_id,
                                     eqc->n_face_dofs,
                                     CS_TYPE_cs_real_t);

  /* Read section */
  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      i_ml_id,
                                      eqc->n_face_dofs,
                                      CS_TYPE_cs_real_t,
                                      eqc->face_values);

  /* Handle boundary faces */
  /* ===================== */

  const int  b_ml_id = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  cs_real_t  *b_values =
    eqc->face_values + eqc->n_face_dofs*cs_shared_quant->n_i_faces;

  /* Define the section name */
  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Check section */
  retcode = cs_restart_check_section(restart,
                                     sec_name,
                                     b_ml_id,
                                     eqc->n_face_dofs,
                                     CS_TYPE_cs_real_t);

  /* Read section */
  if (retcode == CS_RESTART_SUCCESS)
    retcode = cs_restart_read_section(restart,
                                      sec_name,
                                      b_ml_id,
                                      eqc->n_face_dofs,
                                      CS_TYPE_cs_real_t,
                                      b_values);

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
cs_hho_scaleq_write_restart(cs_restart_t    *restart,
                            const char      *eqname,
                            void            *scheme_context)
{
  /* Only the face values are handled. Cell values are stored in a cs_field_t
     structure and thus are handled automatically. */
  if (restart == NULL)
    return;
  if (eqname == NULL)
    bft_error(__FILE__, __LINE__, 0, " %s: Name is NULL", __func__);

  const cs_hho_scaleq_t  *eqc = (const cs_hho_scaleq_t  *)scheme_context;

  char sec_name[128];

  /* Handle interior faces */
  /* ===================== */

  const int  i_ml_id = cs_mesh_location_get_id_by_name(N_("interior_faces"));

  /* Define the section name */
  snprintf(sec_name, 127, "%s::i_face_vals", eqname);

  /* Write interior face section */
  cs_restart_write_section(restart,
                           sec_name,
                           i_ml_id,
                           eqc->n_face_dofs,
                           CS_TYPE_cs_real_t,
                           eqc->face_values);

  /* Handle boundary faces */
  /* ===================== */

  const int  b_ml_id = cs_mesh_location_get_id_by_name(N_("boundary_faces"));
  const cs_real_t  *b_values =
    eqc->face_values + eqc->n_face_dofs*cs_shared_quant->n_i_faces;

  /* Define the section name */
  snprintf(sec_name, 127, "%s::b_face_vals", eqname);

  /* Write boundary face section */
  cs_restart_write_section(restart,
                           sec_name,
                           b_ml_id,
                           eqc->n_face_dofs,
                           CS_TYPE_cs_real_t,
                           b_values);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to cs_hho_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hho_scaleq_extra_op(const char                 *eqname,
                       const cs_field_t           *field,
                       const cs_equation_param_t  *eqp,
                       cs_equation_builder_t      *eqb,
                       void                       *data)
{
  cs_hho_scaleq_t  *eqc = (cs_hho_scaleq_t  *)data;

  // TODO
  CS_UNUSED(eqp);
  CS_UNUSED(eqb);
  CS_UNUSED(eqc);
  CS_UNUSED(field);
  CS_UNUSED(eqname);
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
