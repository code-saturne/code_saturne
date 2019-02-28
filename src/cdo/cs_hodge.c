/*============================================================================
 * Build discrete Hodge operators
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
#include <limits.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_log.h"
#include "cs_math.h"
#include "cs_scheme_geometry.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_hodge.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_hodge.c

  \brief Build discrete Hodge operators

*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_HODGE_DBG       0
#define CS_HODGE_MODULO    1

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Local variables
 *============================================================================*/

/*============================================================================
 * Private constant variables
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

static const double  cs_hodge_vc_coef = 3./20;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Check the coherency of the values of a stiffness matrix
 *
 * \param[in] sloc       pointer to a cs_sdm_t struct.
 */
/*----------------------------------------------------------------------------*/

inline static void
_check_stiffness(const cs_sdm_t       *sloc)
{
  assert(sloc != NULL);

  double  print_val = 0.;

  for (int i = 0 ; i < sloc->n_rows; i++) {
    double  rsum = 0.;
    const cs_real_t  *rval = sloc->val + i*sloc->n_rows;
    for (cs_lnum_t j = 0; j < sloc->n_rows; j++)
      rsum += rval[j];

    print_val += fabs(rsum);

    if (rsum > 100*cs_math_get_machine_epsilon()) {
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT, " %s: row %d Row sum = %5.3e > 0 !\n",
                    __func__, i, rsum);
    }

  }
  cs_log_printf(CS_LOG_DEFAULT, " %s: err = % -5.3e\n", __func__, print_val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Check the coherency of the values of a discrete Hodge operator
 *
 * \param[in]      vec     vectors of quantities to test against a hodge
 * \param[in]      res     vectors of quantities to compare with
 * \param[in]      hdg     pointer to a cs_sdm_t structure
 * \param[in]      h_info  parametrization of the discrete Hodge operator
 * \param[in, out] cb      pointer to a cell builder structure
 *                         buffers to store temporary values
 */
/*----------------------------------------------------------------------------*/

inline static void
_check_vector_hodge(const cs_real_3_t       *vec,
                    const cs_real_3_t       *res,
                    const cs_sdm_t       *hdg,
                    const cs_param_hodge_t   h_info,
                    cs_cell_builder_t       *cb)

{
  assert(hdg != NULL);

  cs_real_t  *in = cb->values;
  cs_real_t  *h_in = cb->values + hdg->n_rows;
  cs_real_t  *ref = cb->values + 2*hdg->n_rows;
  double  print_val = 0.;

  if (h_info.is_unity) {
    cb->dpty_mat[0][0] = cb->dpty_mat[1][1] = cb->dpty_mat[2][2] = 1;
    cb->dpty_mat[0][1] = cb->dpty_mat[1][0] = cb->dpty_mat[2][0] = 0;
    cb->dpty_mat[0][2] = cb->dpty_mat[1][2] = cb->dpty_mat[2][1] = 0;
  }
  else if (h_info.is_iso) {
    cb->dpty_mat[0][0] = cb->dpty_mat[1][1] = cb->dpty_mat[2][2] = cb->dpty_val;
    cb->dpty_mat[0][1] = cb->dpty_mat[1][0] = cb->dpty_mat[2][0] = 0;
    cb->dpty_mat[0][2] = cb->dpty_mat[1][2] = cb->dpty_mat[2][1] = 0;
  }

  const cs_real_3_t  a[3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };

  for (int dim = 0; dim < 3; dim++) {

    cs_real_3_t  pty_a;
    cs_math_33_3_product((const cs_real_t (*)[3])cb->dpty_mat, a[dim], pty_a);

    for (int i = 0; i < hdg->n_rows; i++) {
      in[i] = vec[i][dim];
      ref[i] = _dp3(pty_a, res[i]);
    }

    cs_sdm_square_matvec(hdg, in, h_in);

    double  err = 0.;
    for (int i = 0; i < hdg->n_rows; i++)
      err += fabs(ref[i] - h_in[i]);
    print_val += err;
    if (err > 100 * cs_math_get_machine_epsilon()) {
      cs_base_warn(__FILE__, __LINE__);
      cs_log_printf(CS_LOG_DEFAULT,
                    " %s: err = %5.3e, dim = %d\n", __func__, err, dim);
    }

  }

  cs_log_printf(CS_LOG_DEFAULT, "%s: err = % -5.3e\n", __func__, print_val);

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local builder structure used for building Hodge op.
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

  const int  n_vc = connect->n_max_vbyc;
  const int  n_ec = connect->n_max_ebyc;
  const int  n_fc = connect->n_max_fbyc;

  cs_cell_builder_t *cb = cs_cell_builder_create();

  switch (space_scheme) {

  case CS_SPACE_SCHEME_CDOVB:
    size = CS_MAX(4*n_ec + 3*n_vc, n_ec*(n_ec+1));
    BFT_MALLOC(cb->values, size, double);
    memset(cb->values, 0, size*sizeof(cs_real_t));

    size = 2*n_ec;
    BFT_MALLOC(cb->vectors, size, cs_real_3_t);
    memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

    cb->hdg = cs_sdm_square_create(n_ec);
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    size = 2*n_vc + 3*n_ec + n_fc;
    BFT_MALLOC(cb->values, size, double);
    memset(cb->values, 0, size*sizeof(cs_real_t));

    size = 2*n_ec + n_vc;
    BFT_MALLOC(cb->vectors, size, cs_real_3_t);
    memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

    cb->hdg = cs_sdm_square_create(n_vc + 1);
    break;

  case CS_SPACE_SCHEME_CDOFB:
    size = n_fc*(n_fc+1);
    BFT_MALLOC(cb->values, size, double);
    memset(cb->values, 0, size*sizeof(cs_real_t));

    size = 2*n_fc;
    BFT_MALLOC(cb->vectors, size, cs_real_3_t);
    memset(cb->vectors, 0, size*sizeof(cs_real_3_t));

    cb->hdg = cs_sdm_square_create(n_fc);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _("Invalid space scheme."));

  } // End of switch on space scheme

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quantities used for defining the entries of the discrete
 *          Hodge for COST algo. when the property is isotropic
 *          Initialize the local discrete Hodge op. with the consistency part
 *
 * \param[in]      n_ent      number of local entities
 * \param[in]      invcvol    1/|c|
 * \param[in]      ptyval     values of property inside this cell
 * \param[in]      pq         pointer to the first set of quantities
 * \param[in]      dq         pointer to the second set of quantities
 * \param[in, out] alpha      geometrical quantity
 * \param[in, out] kappa      geometrical quantity
 * \param[in, out] hloc       pointer to a cs_sdm_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cost_quant_iso(const int               n_ent,
                        const double            invcvol,
                        const double            ptyval,
                        const cs_real_3_t      *pq,
                        const cs_real_3_t      *dq,
                        double                 *alpha,
                        double                 *kappa,
                        cs_sdm_t            *hloc)
{
  /* Compute several useful quantities
     alpha_ij = delta_ij - pq_j.Consist_i where Consist_i = 1/|c| dq_i
     qmq_ii = dq_i.ptyval.dq_i
     kappa_i = qmq_ii / |subvol_i|
  */

  for (int i = 0; i < n_ent; i++) {

    const double  dsvol_i = _dp3(dq[i], pq[i]);

    double  *alpha_i = alpha + i*n_ent;
    double  *mi = hloc->val + i*n_ent;

    alpha_i[i] = 1 - invcvol * dsvol_i;

    const double  qmq_ii = ptyval * _dp3(dq[i], dq[i]);

    mi[i] = invcvol * qmq_ii;
    kappa[i] = 3. * qmq_ii / dsvol_i;

    for (int j = i+1; j < n_ent; j++) {

      /* Initialize the lower left part of hloc with the consistency part */
      mi[j] = invcvol * ptyval * _dp3(dq[j], dq[i]);

      /* Compute the alpha matrix (not symmetric) */
      alpha_i[j] = -invcvol * _dp3(pq[j], dq[i]);
      alpha[j*n_ent + i] = -invcvol * _dp3(pq[i], dq[j]);

    } /* Loop on entities (J) */

  } /* Loop on entities (I) */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quantities used for defining the entries of the discrete
 *          Hodge for COST algo.
 *          Initialize the local discrete Hodge op. with the consistency part
 *
 * \param[in]      n_ent    number of local entities
 * \param[in]      ovc      1/|c| where |c| is the cell volume
 * \param[in]      pty      values of the tensor related to the material pty
 * \param[in]      pq       pointer to the first set of quantities
 * \param[in]      dq       pointer to the second set of quantities
 * \param[in, out] alpha    geometrical quantity
 * \param[in, out] kappa    geometrical quantity
 * \param[in, out] hloc     pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cost_quant(const int               n_ent,
                    const double            ovc,
                    const cs_real_33_t      pty,
                    const cs_real_3_t      *pq,
                    const cs_real_3_t      *dq,
                    double                 *alpha,
                    double                 *kappa,
                    cs_sdm_t               *hloc)
{
  /* Compute several useful quantities:
       alpha_ij = delta_ij - 1/|c|*pq_j.dq_i
       qmq_ii = dq_i.mat.dq_i
       kappa_i = qmq_ii / |subvol_i| */

  for (int i = 0; i < n_ent; i++) {

    const  double  mdq_i[3]
      = { pty[0][0] * dq[i][0] + pty[0][1] * dq[i][1] + pty[0][2] * dq[i][2],
          pty[1][0] * dq[i][0] + pty[1][1] * dq[i][1] + pty[1][2] * dq[i][2],
          pty[2][0] * dq[i][0] + pty[2][1] * dq[i][1] + pty[2][2] * dq[i][2] };

    double  *h_i = hloc->val + i*n_ent;

    h_i[i] = _dp3(dq[i], mdq_i); /* Add the consistent part */

    double  *alpha_i = alpha + i*n_ent;

    alpha_i[i] = _dp3(dq[i], pq[i]);
    kappa[i] = 3. * h_i[i] / alpha_i[i];
    alpha_i[i] = 1 - ovc*alpha_i[i];
    h_i[i] *= ovc;

    for (int j = i+1; j < n_ent; j++) {

      /* Initialize the lower left part of hloc with the consistency part */
      h_i[j] = ovc * _dp3(dq[j], mdq_i);

      /* Compute the alpha matrix (not symmetric) */
      alpha_i[j] = -ovc * _dp3(pq[j], dq[i]);

    }

    const  double  opq_i[3] = { -ovc*pq[i][0], -ovc*pq[i][1], -ovc*pq[i][2] };
    for (int j = i+1; j < n_ent; j++)
      alpha[j*n_ent + i] = _dp3(opq_i, dq[j]);

  } /* Loop on entities (I) */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local discrete Hodge operator using the generic COST algo.
 *          and a cellwise view of the mesh. Specific for EpFd Hodge operator.
 *          COST means COnsistency + STabilization
 *
 * \param[in]      n_ent    number of local entities
 * \param[in]      beta2    squared value of the stabilization coefficient
 * \param[in]      alpha    precomputed quantities
 * \param[in]      kappa    precomputed quantities
 * \param[in, out] hval     values of the coefficient of the Hodge matrix
 */
/*----------------------------------------------------------------------------*/

static void
_compute_hodge_cost(const int       n_ent,
                    const double    beta2,
                    const double    alpha[],
                    const double    kappa[],
                    double          hval[])
{
  for (int i = 0; i < n_ent; i++) { /* Loop on cell entities I */

    const int  shift_i = i*n_ent;
    const double  *alpha_i = alpha + shift_i;

    double  *mi = hval + shift_i;

    /* Add contribution from the stabilization part for
       each sub-volume related to a primal entity */
    double  stab_part = 0;
    for (int k = 0; k < n_ent; k++) /* Loop over sub-volumes */
      stab_part += kappa[k] * alpha_i[k] * alpha_i[k];

    mi[i] += beta2 * stab_part; // Consistency part has already been computed

    /* Compute extra-diag entries */
    for (int j = i + 1; j < n_ent; j++) { /* Loop on cell entities J */

      const int  shift_j = j*n_ent;
      const double  *alpha_j = alpha + shift_j;
      double  *mj = hval + shift_j;

      /* Add contribution from the stabilization part for
         each sub-volume related to a primal entity */
      stab_part = 0;
      for (int k = 0; k < n_ent; k++) /* Loop over sub-volumes */
        stab_part += kappa[k] * alpha_i[k] * alpha_j[k];

      mi[j] += beta2 * stab_part;
      mj[i] = mi[j]; // Symmetric by construction

    } /* End of loop on J entities */

  } /* End of loop on I entities */

}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc
 *          Case of CDO face-based schemes
 *
 * \param[in]      h_info     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fb_cost_get_stiffness(const cs_param_hodge_t    h_info,
                               const cs_cell_mesh_t     *cm,
                               cs_cell_builder_t        *cb)
{
  /* Sanity checks */
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ));

  /* Initialize the local stiffness matrix */
  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_fc + 1, sloc);

  /* Compute the local discrete Hodge operator */
  cs_hodge_edfp_cost_get(h_info, cm, cb);

  cs_sdm_t  *hloc = cb->hdg;
  double  *sval_crow = sloc->val + cm->n_fc*sloc->n_rows;
  double  full_sum = 0.;

  for (int i = 0; i < hloc->n_rows; i++) {

    const short int  fi_sgn = cm->f_sgn[i];
    const double  *hval_i = hloc->val + i*hloc->n_rows;

    double  *sval_i = sloc->val + i*sloc->n_rows;
    double  row_sum = 0.;
    for (int j = 0; j < hloc->n_rows; j++) {
      const cs_real_t  hsgn_coef = fi_sgn * cm->f_sgn[j] * hval_i[j];
      row_sum += hsgn_coef;
      sval_i[j] = hsgn_coef;
    }

    sval_i[cm->n_fc] = -row_sum;
    sval_crow[i] = -row_sum;
    full_sum += row_sum;

  }

  /* (c, c) diagonal entry */
  sval_crow[cm->n_fc] = full_sum;

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, ">> Local stiffness matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, sloc);
    _check_stiffness(sloc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Voronoi algorithm
 *          The computed matrix is stored in cb->loc
 *          Case of CDO face-based schemes
 *
 * \param[in]      h_info     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fb_voro_get_stiffness(const cs_param_hodge_t    h_info,
                               const cs_cell_mesh_t     *cm,
                               cs_cell_builder_t        *cb)
{
  /* Sanity checks */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EDFP);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PF | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ));

  /* Compute the local discrete Hodge operator */
  cs_hodge_edfp_voro_get(h_info, cm, cb);

  /* Initialize the local stiffness matrix */
  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_fc + 1, sloc);

  double  full_sum = 0.;

  cs_sdm_t  *hloc = cb->hdg;
  double  *sval_crow = sloc->val + cm->n_fc*sloc->n_rows;

  for (int i = 0; i < hloc->n_rows; i++) {

    /* Hodge operator is diagonal */
    const double  *hval_i = hloc->val + i*hloc->n_rows;
    const double  row_sum = hval_i[i];

    double  *sval_i = sloc->val + i*sloc->n_rows;

    sval_i[i] = hval_i[i];
    sval_i[cm->n_fc] = -row_sum;
    sval_crow[i] = -row_sum;
    full_sum += row_sum;

  }

  /* (c, c) diagonal entry */
  sval_crow[cm->n_fc] = full_sum;

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, ">> Local stiffness matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, sloc);
    _check_stiffness(sloc);
  }
#endif

  bft_error(__FILE__, __LINE__, 0, "Under construction");
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          Case of CDO vertex-based schemes
 *
 * \param[in]      h_info     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_cost_get_stiffness(const cs_param_hodge_t    h_info,
                               const cs_cell_mesh_t     *cm,
                               cs_cell_builder_t        *cb)
{
  /* Sanity checks */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EPFD);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ |
                      CS_CDO_LOCAL_EV));

  double  *kappa = cb->values;
  double  *alpha = cb->values + cm->n_ec;
  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_ec;

  /* Initialize the local stiffness matrix */
  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_vc, sloc);

  /* Initialize the hodge matrix */
  cs_sdm_t  *hloc = cb->hdg;
  cs_sdm_square_init(cm->n_ec, hloc);

  /* Set numbering and geometrical quantities Hodge builder */
  for (int ii = 0; ii < cm->n_ec; ii++) {

    cs_nvec3_t  dfq = cm->dface[ii];
    cs_quant_t  peq = cm->edge[ii];

    for (int k = 0; k < 3; k++) {
      dq[ii][k] = dfq.meas * dfq.unitv[k];
      pq[ii][k] = peq.meas * peq.unitv[k];
    }

  } /* Loop on cell edges */

  /* Compute additional geometrical quantities.
     Initial the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  const double  invcvol = 1/cm->vol_c;
  const double  beta2 = h_info.coef*h_info.coef;

  if (h_info.is_iso || h_info.is_unity) {

    double  dpty_val = 0.;
    if (h_info.is_unity)
      dpty_val = 1.0;
    else if (h_info.is_iso)
      dpty_val = cb->dpty_val;

    _compute_cost_quant_iso(cm->n_ec, invcvol, dpty_val,
                            (const cs_real_t (*)[3])pq,
                            (const cs_real_t (*)[3])dq,
                            alpha, kappa, hloc);

  }
  else
    _compute_cost_quant(cm->n_ec, invcvol,
                        (const cs_real_3_t *)cb->dpty_mat,
                        (const cs_real_t (*)[3])pq,
                        (const cs_real_t (*)[3])dq,
                        alpha, kappa, hloc);

  for (int ei = 0; ei < cm->n_ec; ei++) { /* Loop on cell edges I */

    const int  shift_i = ei*cm->n_ec;
    const double  *alpha_i = alpha + shift_i;
    const short int  i1ei = cm->e2v_sgn[ei];
    const short int  i1 = cm->e2v_ids[2*ei];
    const short int  i2 = cm->e2v_ids[2*ei+1];
    const double  *hi = hloc->val + shift_i;

    assert(i1 < i2);

    double  *si1 = sloc->val + i1*sloc->n_rows;
    double  *si2 = sloc->val + i2*sloc->n_rows;

    /* Add contribution from the stabilization part for each sub-volume
       related to a primal entity */
    double  stab_part = 0;
    for (int ek = 0; ek < cm->n_ec; ek++) /* Loop over sub-volumes */
      stab_part += kappa[ek] * alpha_i[ek] * alpha_i[ek];

    /* Diagonal value: consistency part has already been computed */
    const double  dval = hi[ei] + beta2 * stab_part;

    si1[i1] += dval;
    si1[i2] -= dval;
    si2[i2] += dval;

    /* Compute extra-diag entries */
    for (int ej = ei + 1; ej < cm->n_ec; ej++) { /* Loop on cell entities J */

      const int  shift_j = ej*cm->n_ec;
      const double  *alpha_j = alpha + shift_j;
      const short int  j1ej = cm->e2v_sgn[ej];
      const short int  j1 = cm->e2v_ids[2*ej];
      const short int  j2 = cm->e2v_ids[2*ej+1];

      assert(j1 < j2);

      double  *sj1 = sloc->val + j1*sloc->n_rows;
      double  *sj2 = sloc->val + j2*sloc->n_rows;

      /* Add contribution from the stabilization part for each sub-volume
         related to a primal entity */
      stab_part = 0;
      for (int ek = 0; ek < cm->n_ec; ek++) /* Loop over sub-volumes */
        stab_part += kappa[ek] * alpha_i[ek] * alpha_j[ek];

      /* Extra-diagonal value */
      const double xval = (hi[ej] + beta2 * stab_part) * i1ei * j1ej;

      if (i2 < j1) {            /* i1 < i2 < j1 < j2 */
        si1[j1] += xval;
        si1[j2] -= xval;
        si2[j1] -= xval;
        si2[j2] += xval;
      }
      else if (i2 == j1) {      /* i1 < i2 = j1 < j2 */
        si1[j1] += xval;
        si1[j2] -= xval;
        si2[j1] -= 2*xval;
        si2[j2] += xval;
      }
      else if (i2 < j2) {

        assert(i2 > j1);
        if (i1 < j1)            /* i1 < j1 < i2 < j2 */
          si1[j1] += xval;
        else if (i1 == j1)      /* i1 = j1 < i2 < j2 */
          si1[j1] += 2*xval;
        else                    /* j1 < i1 < i2 < j2 */
          sj1[i1] += xval;

        si1[j2] -= xval;
        sj1[i2] -= xval;
        si2[j2] += xval;

      }
      else if (i2 == j2) {

        if (i1 < j1)            /* i1 < j1 < i2 = j2 */
          si1[j1] += xval;
        else if (i1 == j1)      /* i1 = j1 < i2 = j2 */
          si1[j1] += 2*xval;
        else                    /* j1 < i1 < i2 = j2 */
          sj1[i1] += xval;

        si1[j2] -= xval;
        sj1[i2] -= xval;
        si2[j2] += 2*xval;

      }
      else {                    /* i2 > j2 */

        if (i1 < j1) {          /* i1 < j1 < j2 < i2 */
          si1[j1] += xval;
          si1[j2] -= xval;
        }
        else if (i1 == j1) {    /* i1 = j1 < j2 < i2 */
          si1[j1] += 2*xval;
          si1[j2] -= xval;
        }
        else if (i1 < j2) {     /* j1 < i1 < j2 < i2 */
          sj1[i1] += xval;
          si1[j2] -= xval;
        }
        else if (i1 == j2) {    /* j1 < i1 = j2 < i2 */
          sj1[i1] += xval;
          si1[j2] -= 2*xval;
        }
        else {                  /* j1 < j2 < i1 < i2 */
          sj1[i1] += xval;
          sj2[i1] -= xval;
        }

        assert(i2 > j2);
        sj1[i2] -= xval;
        sj2[i2] += xval;

      } /* End of tests */

    } /* End of loop on J entities */

  } /* End of loop on I entities */

  /* Stiffness matrix is symmetric by construction */
  for (int ei = 0; ei < sloc->n_rows; ei++) {
    double *si = sloc->val + ei*sloc->n_rows;
    for (int ej = 0; ej < ei; ej++)
      si[ej] = sloc->val[ej*sloc->n_rows + ei];
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, ">> Local stiffness matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, sloc);
    _check_stiffness(sloc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Voronoi algorithm
 *          The computed matrix is stored in cb->loc
 *          Case of CDO vertex-based schemes
 *
 * \param[in]      h_info     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_voro_get_stiffness(const cs_param_hodge_t    h_info,
                               const cs_cell_mesh_t     *cm,
                               cs_cell_builder_t        *cb)
{
  /* Sanity checks */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EPFD);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ |
                      CS_CDO_LOCAL_EV));

  /* Initialize the local stiffness matrix */
  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_vc, sloc);

  if (h_info.is_iso || h_info.is_unity) {

    double  dpty_val;
    if (h_info.is_unity)
      dpty_val = 1.0;
    else if (h_info.is_iso)
      dpty_val = cb->dpty_val;

    /* Loop on cell edges */
    for (int ii = 0; ii < cm->n_ec; ii++) {

      cs_nvec3_t  dfq = cm->dface[ii];
      cs_quant_t  peq = cm->edge[ii];

      /* Only a diagonal term */
      const double  dval = dpty_val * dfq.meas/peq.meas;
      const short int  vi = cm->e2v_ids[2*ii];
      const short int  vj = cm->e2v_ids[2*ii+1];

      double  *si = sloc->val + vi*sloc->n_rows;
      double  *sj = sloc->val + vj*sloc->n_rows;

      si[vi] += dval;
      sj[vj] += dval;
      si[vj] = sj[vi] = -dval; // sgn_i * sgn_j = -1

    } /* End of loop on cell edges */

  }
  else { /* Diffusion property is anisotropic */

    cs_real_3_t  mv;

    /* Loop on cell edges */
    for (int ii = 0; ii < cm->n_ec; ii++) {

      cs_nvec3_t  dfq = cm->dface[ii];
      cs_quant_t  peq = cm->edge[ii];

      cs_math_33_3_product((const cs_real_3_t *)cb->dpty_mat, dfq.unitv, mv);

      /* Only a diagonal term */
      const double  dval = _dp3(mv, dfq.unitv) * dfq.meas/peq.meas;
      const short int  vi = cm->e2v_ids[2*ii];
      const short int  vj = cm->e2v_ids[2*ii+1];

      double  *si = sloc->val + vi*sloc->n_rows;
      double  *sj = sloc->val + vj*sloc->n_rows;

      si[vi] += dval;
      sj[vj] += dval;
      si[vj] = sj[vi] = -dval; // sgn_j * sgn_i = -1

    } /* End of loop on cell edges */

  } /* Tensor-valued property */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, ">> Local stiffness matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, sloc);
    _check_stiffness(sloc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic WBS algo.
 *          WBS standing for Whitney Barycentric Subdivision (WBS)
 *          algo.
 *
 * \param[in]      h_info     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vb_wbs_get_stiffness(const cs_param_hodge_t    h_info,
                              const cs_cell_mesh_t     *cm,
                              cs_cell_builder_t        *cb)
{
  /* Sanity checks */
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EPFD);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_WBS);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ |
                      CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_HFQ | CS_CDO_LOCAL_FEQ));

  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg;

  cs_real_3_t  *uvc = cb->vectors;
  cs_real_3_t  *glv = cb->vectors + cm->n_vc;
  cs_real_t  *lvc = cb->values;
  cs_real_t  *wvf = cb->values + cm->n_vc;
  cs_real_t  *pefc_vol = cb->values + 2*cm->n_vc;

  /* Set the diffusion tensor */
  cs_real_33_t  tensor = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  if (h_info.is_iso) {
    if (!h_info.is_unity)
      tensor[0][0] = tensor[1][1] = tensor[2][2] = cb->dpty_val;
  }
  else {
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        tensor[k][l] = cb->dpty_mat[k][l];
  }

  /* Initialize the local stiffness matrix */
  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_vc, sloc);

  /* Define the length and unit vector of the segment x_c --> x_v */
  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, lvc + v, uvc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute the gradient of the lagrange function related to a cell
       in each p_{f,c} and the weights for each vertex related to this face */
    cs_compute_fwbs_q2(f, cm, grd_c, wvf, pefc_vol);

    /* Loop on face edges to scan p_{e,f,c} subvolumes */
    for (int i = cm->f2e_idx[f], jj = 0; i < cm->f2e_idx[f+1]; i++, jj++) {

      const double  subvol = pefc_vol[jj];
      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];

      /* Gradient of the lagrange function related to v1 and v2 */
      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])uvc, lvc,
                        grd_v1, grd_v2);

      /* Gradient of the lagrange function related to a face.
         This formula is a consequence of the Partition of the Unity */
      for (int k = 0; k < 3; k++)
        grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

      /* Compute the gradient of the conforming reconstruction functions for
         each vertex of the cell in this subvol (pefc) */
      for (int si = 0; si < sloc->n_rows; si++) {

        for (int k = 0; k < 3; k++)
          glv[si][k] = cm->wvc[si]*grd_c[k];

        if (wvf[si] > 0) // Face contrib.
          for (int k = 0; k < 3; k++)
            glv[si][k] += wvf[si]*grd_f[k];

        if (si == v1) // Vertex 1 contrib
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v1[k];

        if (si == v2) // Vertex 2 contrib
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v2[k];

      } // Loop on cell vertices

      /* Build the upper right part */
      for (int si = 0; si < sloc->n_rows; si++) {

        cs_math_33_3_product((const cs_real_t (*)[3])tensor, glv[si], matg);

        /* Diagonal contribution */
        double  *mi = sloc->val + si*sloc->n_rows;
        mi[si] += subvol * _dp3(matg, glv[si]);

        /* Loop on vertices v_j (j > i) */
        for (int sj = si+1; sj < sloc->n_rows; sj++)
          mi[sj] += subvol * _dp3(matg, glv[sj]);

      } /* Loop on vertices v_i */

    }

  } // Loop on cell faces

  /* Matrix is symmetric by construction */
  for (int si = 0; si < sloc->n_rows; si++) {
    double  *mi = sloc->val + si*sloc->n_rows;
    for (int sj = si+1; sj < sloc->n_rows; sj++)
      sloc->val[sj*sloc->n_rows+si] = mi[sj];
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, ">> Local stiffness matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, sloc);
    _check_stiffness(sloc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic WBS algo.
 *          WBS standing for Whitney Barycentric Subdivision (WBS) algo.
 *          The computed matrix is stored in cb->loc
 *
 * \param[in]      h_info     pointer to a cs_param_hodge_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vcb_get_stiffness(const cs_param_hodge_t    h_info,
                           const cs_cell_mesh_t     *cm,
                           cs_cell_builder_t        *cb)
{
  /* Sanity checks */
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ));

  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg, matg_c;

  cs_real_3_t  *uvc = cb->vectors;
  cs_real_3_t  *glv = cb->vectors + cm->n_vc;
  cs_real_t  *lvc = cb->values;
  cs_real_t  *wvf = cb->values + cm->n_vc;
  cs_real_t  *pefc_vol = cb->values + 2*cm->n_vc;

  /* Set the diffusion tensor */
  cs_real_33_t  tensor = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  if (h_info.is_iso) {
    if (!h_info.is_unity)
      tensor[0][0] = tensor[1][1] = tensor[2][2] = cb->dpty_val;
  }
  else {
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        tensor[k][l] = cb->dpty_mat[k][l];
  }

  /* Initialize the local stiffness matrix */
  const int  nc_dofs = cm->n_vc + 1;
  /* index to the (cell,cell) entry */
  const int  cc = nc_dofs*cm->n_vc + cm->n_vc;

  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(nc_dofs, sloc);

  /* Define the length and unit vector of the segment x_c --> x_v */
  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, lvc + v, uvc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute for the current face:
       - the gradient of the Lagrange function related xc in p_{f,c}
       - weights related to vertices
       - subvolume p_{ef,c} related to edges
    */
    const double  pfc_vol = cs_compute_fwbs_q3(f, cm, grd_c, wvf, pefc_vol);

    /* Compute the contribution to the entry A(c,c) */
    cs_math_33_3_product((const cs_real_t (*)[3])tensor, grd_c, matg_c);
    sloc->val[cc] += pfc_vol * _dp3(grd_c, matg_c);

    /* Loop on face edges to scan p_{e,f,c} subvolumes */
    for (int i = cm->f2e_idx[f], jj = 0; i < cm->f2e_idx[f+1]; i++, jj++) {

      const double  subvol = pefc_vol[jj];
      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];

      /* Gradient of the lagrange function related to v1 and v2 */
      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])uvc, lvc,
                        grd_v1, grd_v2);

      /* Gradient of the Lagrange function related to a face.
         This formula is a consequence of the Partition of the Unity */
      for (int k = 0; k < 3; k++)
        grd_f[k] = -(grd_c[k] + grd_v1[k] + grd_v2[k]);

      /* Compute the gradient of the conforming reconstruction functions for
         each vertex of the cell in this subvol (pefc) */
      for (int si = 0; si < cm->n_vc; si++) {

        for (int k = 0; k < 3; k++)
          glv[si][k] = 0;

        if (wvf[si] > 0) // Face contrib.
          for (int k = 0; k < 3; k++)
            glv[si][k] += wvf[si]*grd_f[k];

        if (si == v1) // Vertex 1 contrib
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v1[k];

        if (si == v2) // Vertex 2 contrib
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v2[k];

      } // Loop on cell vertices

      /* Build the upper right part (v-v and v-c)
         Be careful: sloc->n_rows = cm->n_vc + 1 */
      for (int si = 0; si < cm->n_vc; si++) {

        double  *mi = sloc->val + si*sloc->n_rows;

        /* Add v-c contribution */
        mi[cm->n_vc] += subvol * _dp3(matg_c, glv[si]);

        /* Add v-v contribution */
        cs_math_33_3_product((const cs_real_t (*)[3])tensor, glv[si], matg);

        /* Loop on vertices v_j (j >= i) */
        for (int sj = si; sj < cm->n_vc; sj++)
          mi[sj] += subvol * _dp3(matg, glv[sj]);

      } /* Loop on vertices v_i */

    }

  } // Loop on cell faces

  /* Matrix is symmetric by construction */
  for (int si = 0; si < sloc->n_rows; si++) {
    double  *mi = sloc->val + si*sloc->n_rows;
    for (int sj = si+1; sj < sloc->n_rows; sj++)
      sloc->val[sj*sloc->n_rows+si] = mi[sj];
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, ">> Local stiffness matrix");
    cs_sdm_dump(cm->c_id, NULL, NULL, sloc);
    _check_stiffness(sloc);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a local Hodge operator on a given cell which is equivalent of
 *         a mass matrix. It relies on a CO+ST algo. and is specific to CDO-Fb
 *         schemes.
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fb_get_mass(const cs_param_hodge_t    h_info,
                     const cs_cell_mesh_t     *cm,
                     cs_cell_builder_t        *cb)
{
  CS_UNUSED(h_info);            /* Only in debug mode */

  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_FB);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ));

  const int n_cols = cm->n_fc + 1;
  const cs_real_t  over_cell = 1./(cm->vol_c*cm->vol_c);

  /* Initialize the local matrix related to this discrete Hodge operator */
  cs_sdm_t  *hdg = cb->hdg;
  cs_sdm_square_init(n_cols, hdg);

  /* cell-cell entry (cell-face and face-cell block are null */
  hdg->val[cm->n_fc*(n_cols + 1)] = cm->vol_c;

  /* Compute the inertia (useful for the face-face block */
  cs_real_t  inertia_tensor[3][3];
  cs_compute_inertia_tensor(cm, cm->xc, inertia_tensor);

  for (short int fi = 0; fi < cm->n_fc; fi++) {

    const cs_quant_t  pfq_i = cm->face[fi];
    const short int  sgn_i = cm->f_sgn[fi];

    /* Diagonal entry (a bit more optimized) */
    cs_real_t dval = 0;
    for (int k = 0; k < 3; k++) {
      dval += pfq_i.unitv[k]*pfq_i.unitv[k]*inertia_tensor[k][k];
      for (int l = k+1; l < 3; l++)
        dval += 2*pfq_i.unitv[k]*pfq_i.unitv[l]*inertia_tensor[k][l];
    }

    hdg->val[fi*(n_cols+1)] = over_cell * pfq_i.meas*pfq_i.meas * dval;

    /* Extra-diag entry */
    for (short int fj = fi + 1; fj < cm->n_fc; fj++) {

      const cs_quant_t  pfq_j = cm->face[fj];

      cs_real_t  eval = 0;
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++)
          eval += pfq_i.unitv[k]*pfq_j.unitv[l]*inertia_tensor[k][l];
      }
      eval *= over_cell * sgn_i*pfq_i.meas * cm->f_sgn[fj]*pfq_j.meas;

      /* Symmetric by construction */
      hdg->val[fi*n_cols+fj] = eval;
      hdg->val[fj*n_cols+fi] = eval;

    } /* Loop on cell faces (fj) */

  } /* Loop on cell faces (fi) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the WBS algo.
 *          This function is specific for vertex+cell-based schemes
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vcb_wbs_get(const cs_param_hodge_t    h_info,
                     const cs_cell_mesh_t     *cm,
                     cs_cell_builder_t        *cb)
{
  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_VC);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_WBS);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ |
                      CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ));

  cs_sdm_t  *hdg = cb->hdg;

  /* Initialize the local matrix related to this discrete Hodge operator */
  cs_sdm_square_init(cm->n_vc + 1, hdg);

  double  *wvf = cb->values;
  double  *pefc_vol = cb->values + cm->n_vc;

  const int  msize = cm->n_vc + 1;
  const double  c_coef1 = 0.2*cm->vol_c;
  const double  c_coef2 = cs_hodge_vc_coef * cm->vol_c;

  /* H(c,c) = 0.1*|c| */
  hdg->val[msize*cm->n_vc + cm->n_vc] = 0.1*cm->vol_c;

  /* Initialize the upper part of the local Hodge matrix
     diagonal and cell column entries */
  for (short int vi = 0; vi < cm->n_vc; vi++) {

    double  *mi = hdg->val + vi*msize;

    mi[vi] = c_coef1 * cm->wvc[vi];       // Diagonal entry
    for (short int vj = vi+1; vj < cm->n_vc; vj++)
      mi[vj] = 0.;
    mi[cm->n_vc] = c_coef2 * cm->wvc[vi]; // Cell column

  } // Loop on cell vertices

  /* Loop on each pef and add the contribution */
  for (short int f = 0; f < cm->n_fc; f++) {

    /* Define useful quantities for WBS algo. */
    const double pfc_vol = cs_compute_fwbs_q1(f, cm, wvf, pefc_vol);
    const double f_coef = 0.3 * pfc_vol;

    /* Add face contribution:
       Diagonal entry    H(i,i) += 0.3*wif*wif*pfc_vol
       Extra-diag. entry H(i,j) += 0.3*wjf*wif*pfc_vol */
    for (short int vi = 0; vi < cm->n_vc; vi++) {

      const double  coef_if = f_coef * wvf[vi];
      double  *mi = hdg->val + vi*msize;

      /* Diagonal and Extra-diagonal entries: Add face contribution */
      for (short int vj = vi; vj < cm->n_vc; vj++)
        mi[vj] += coef_if * wvf[vj];

    } // Extra-diag entries

    /* Add edge-face contribution (only extra-diag) = 0.05 * |p_{ef,c}| */
    for (int i = cm->f2e_idx[f], ii = 0; i < cm->f2e_idx[f+1]; i++, ii++) {

      const short int  e = cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[2*e];
      const short int  v2 = cm->e2v_ids[2*e+1];

      /* Sanity check */
      assert(v1 > -1 && v2 > -1);
      if (v1 < v2)
        hdg->val[v1*msize+v2] += 0.05 * pefc_vol[ii];
      else
        hdg->val[v2*msize+v1] += 0.05 * pefc_vol[ii];

    } // Loop on face edges

  } // Loop on cell faces

  /* Take into account the value of the associated property */
  if (!h_info.is_unity) {
    for (short int vi = 0; vi < msize; vi++) {
      double  *mi = hdg->val + vi*msize;
      for (short int vj = vi; vj < msize; vj++)
        mi[vj] *= cb->dpty_val;
    }
  }

  /* Local matrix is symmetric by construction. Set the lower part. */
  for (short int vj = 0; vj < msize; vj++) {
    double  *mj = hdg->val + vj*msize;
    for (short int vi = vj+1; vi < msize; vi++)
      hdg->val[vi*msize + vj] = mj[vi];
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using WBS algo.
 *          Hodge op. from primal vertices to dual cells.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vpcd_wbs_get(const cs_param_hodge_t    h_info,
                      const cs_cell_mesh_t     *cm,
                      cs_cell_builder_t        *cb)
{
  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_VPCD);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_WBS);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PVQ |CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ |
                      CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ));

  double  *wvf = cb->values;
  double  *pefc_vol = cb->values + cm->n_vc;

  cs_sdm_t  *hdg = cb->hdg;
  cs_sdm_square_init(cm->n_vc, hdg);

  const double  c_coef = 0.1*cm->vol_c;

  /* Initialize the upper part of the local Hodge matrix */
  for (short int vi = 0; vi < cm->n_vc; vi++) {

    double  *mi = hdg->val + vi*cm->n_vc;
    const double  vi_coef = 4 * c_coef * cm->wvc[vi];

    // Diag. entry has an additional contrib
    mi[vi] = vi_coef * (0.5 + cm->wvc[vi]);
    for (short int vj = vi + 1; vj < cm->n_vc; vj++)
      mi[vj] = vi_coef * cm->wvc[vj]; // Extra-diagonal entries

  } // Loop on cell vertices

  /* Loop on each pef and add the contribution */
  for (short int f = 0; f < cm->n_fc; f++) {

    /* Define useful quantities for WBS algo. */
    const double pfc_vol = cs_compute_fwbs_q1(f, cm, wvf, pefc_vol);
    const double f_coef = 0.3 * pfc_vol;

    /* Add face contribution:
       Diagonal entry    H(i,i) += 0.3*wif*wif*pfc_vol
       Extra-diag. entry H(i,j) += 0.3*wjf*wif*pfc_vol */
    for (short int vi = 0; vi < cm->n_vc; vi++) {

      double  *mi = hdg->val + vi*cm->n_vc;

      /* Diagonal and Extra-diagonal entries: Add face contribution */
      const double  coef_if = f_coef * wvf[vi];
      for (short int vj = vi; vj < cm->n_vc; vj++)
        mi[vj] += coef_if * wvf[vj];

    } // Face contribution

    /* Add edge-face contribution (only extra-diag) = 0.05 * |p_{ef,c}| */
    for (int i = cm->f2e_idx[f], ii = 0; i < cm->f2e_idx[f+1]; i++, ii++) {

      const short int  eshft = 2*cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[eshft];
      const short int  v2 = cm->e2v_ids[eshft+1];

      /* Sanity check */
      assert(v1 > -1 && v2 > -1);

      if (v1 < v2)
        hdg->val[v1*cm->n_vc+v2] += 0.05 * pefc_vol[ii];
      else
        hdg->val[v2*cm->n_vc+v1] += 0.05 * pefc_vol[ii];

    } // Loop on face edges

  } // Loop on cell faces

  /* Take into account the value of the associated property */
  if (!h_info.is_unity) {
    for (short int vi = 0; vi < cm->n_vc; vi++) {
      double  *mi = hdg->val + vi*cm->n_vc;
      for (short int vj = vi; vj < cm->n_vc; vj++)
        mi[vj] *= cb->dpty_val;
    }
  }

  /* Local matrix is symmetric by construction. Set the lower part. */
  for (short int vj = 0; vj < cm->n_vc; vj++) {
    double  *mj = hdg->val + vj*cm->n_vc;
    for (short int vi = vj+1; vi < cm->n_vc; vi++)
      hdg->val[vi*cm->n_vc + vj] = mj[vi];
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal vertices to dual cells.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_vpcd_voro_get(const cs_param_hodge_t    h_info,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb)
{
  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_VPCD);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PVQ));

  cs_sdm_t  *hdg = cb->hdg;
  cs_sdm_square_init(cm->n_vc, hdg);

  if (h_info.is_unity) {

    for (int v = 0; v < cm->n_vc; v++)
      hdg->val[v*cm->n_vc+v] = cm->wvc[v] * cm->vol_c;

  }
  else {

    for (int v = 0; v < cm->n_vc; v++)
      hdg->val[v*cm->n_vc+v] = cb->dpty_val * cm->wvc[v] * cm->vol_c;

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal edges to dual faces.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_epfd_voro_get(const cs_param_hodge_t    h_info,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb)
{
  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EPFD);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI);
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ | CS_CDO_LOCAL_EFQ));

  /* Initialize the local matrix related to this discrete Hodge operator */
  cs_sdm_t  *hdg = cb->hdg;
  cs_sdm_square_init(cm->n_ec, hdg);

  for (short int e = 0; e < cm->n_ec; e++) {

    if (h_info.is_iso) {
      hdg->val[e*cm->n_ec+e] = cb->dpty_val*cm->dface[e].meas/cm->edge[e].meas;
    }
    else {

      const cs_nvec3_t  sef0c = cm->sefc[2*e], sef1c = cm->sefc[2*e+1];
      const cs_real_3_t *tens = (const cs_real_3_t *)cb->dpty_mat;

      cs_real_3_t  mv;

      /* First sub-triangle contribution */
      cs_math_33_3_product(tens, sef0c.unitv, mv);
      hdg->val[e*cm->n_ec+e] = sef0c.meas * _dp3(mv, sef0c.unitv);
      /* Second sub-triangle contribution */
      cs_math_33_3_product(tens, sef1c.unitv, mv);
      hdg->val[e*cm->n_ec+e] += sef1c.meas * _dp3(mv, sef1c.unitv);

      hdg->val[e*cm->n_ec+e] /= cm->edge[e].meas;

    }

  } // Loop on cell edges

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from primal edges to dual faces.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_epfd_cost_get(const cs_param_hodge_t    h_info,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb)
{
  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EPFD);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ));

  /* Initialize the local matrix related to this discrete Hodge operator */
  cs_sdm_t  *hdg = cb->hdg;
  cs_sdm_square_init(cm->n_ec, hdg);

  double  *kappa = cb->values;
  double  *alpha = cb->values + cm->n_ec;
  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_ec;

  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_nvec3_t  dfq = cm->dface[e];
    const cs_quant_t  peq = cm->edge[e];

    for (int k = 0; k < 3; k++) {
      dq[e][k] = dfq.meas * dfq.unitv[k];
      pq[e][k] = peq.meas * peq.unitv[k];
    }

  } /* Loop on cell edges */

  /* Compute additional geometrical quantities: qmq and T
     Initial the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  if (h_info.is_unity)
    _compute_cost_quant_iso(cm->n_ec, 1/cm->vol_c, 1.0,
                            (const cs_real_t (*)[3])pq,
                            (const cs_real_t (*)[3])dq,
                            alpha, kappa, hdg);
  else if (h_info.is_iso)
    _compute_cost_quant_iso(cm->n_ec, 1/cm->vol_c, cb->dpty_val,
                            (const cs_real_t (*)[3])pq,
                            (const cs_real_t (*)[3])dq,
                            alpha, kappa, hdg);
  else
    _compute_cost_quant(cm->n_ec, 1/cm->vol_c,
                        (const cs_real_3_t *)cb->dpty_mat,
                        (const cs_real_t (*)[3])pq,
                        (const cs_real_t (*)[3])dq,
                        alpha, kappa, hdg);

  _compute_hodge_cost(cm->n_ec, h_info.coef*h_info.coef, alpha, kappa,
                      hdg->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
    _check_vector_hodge((const cs_real_t (*)[3])pq,
                        (const cs_real_t (*)[3])dq,
                        hdg, h_info, cb);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal faces to dual edges.
 *          This function is related to cell-based schemes
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fped_voro_get(const cs_param_hodge_t    h_info,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb)
{
  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_FPED);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ));

  /* Initialize the local matrix related to this discrete Hodge operator */
  cs_sdm_t  *hdg = cb->hdg;
  cs_sdm_square_init(cm->n_fc, hdg);

  for (short int f = 0; f < cm->n_fc; f++) {

    if (h_info.is_iso) {
      hdg->val[f*cm->n_fc+f] = cb->dpty_val*cm->face[f].meas/cm->dedge[f].meas;
    }
    else {

      cs_real_3_t  mv;

      const cs_nvec3_t  deq = cm->dedge[f];
      const cs_real_3_t *tens = (const cs_real_3_t *)cb->dpty_mat;

      cs_math_33_3_product(tens, deq.unitv, mv);
      hdg->val[f*cm->n_fc+f] = deq.meas * _dp3(mv, deq.unitv)/cm->face[f].meas;

    }

  } // Loop on cell faces

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from primal faces to dual edges.
 *          This function is related to cell-based schemes
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_fped_cost_get(const cs_param_hodge_t    h_info,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb)
{
  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_FPED);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ));

  /* Initialize the local matrix related to this discrete Hodge operator */
  cs_sdm_t  *hdg = cb->hdg;
  cs_sdm_square_init(cm->n_fc, hdg);

  double  *kappa = cb->values;
  double  *alpha = cb->values + cm->n_ec;
  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_fc;

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_nvec3_t  deq = cm->dedge[f];
    const cs_quant_t  pfq = cm->face[f];

    for (int k = 0; k < 3; k++) {
      dq[f][k] = deq.meas * deq.unitv[k];
      pq[f][k] = pfq.meas * pfq.unitv[k];
    }

  } /* Loop on cell faces */

  /* Compute additional geometrical quantities:
     Initialize the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  if (h_info.is_unity)
    _compute_cost_quant_iso(cm->n_fc, 1/cm->vol_c, 1.0,
                            (const cs_real_t (*)[3])pq,
                            (const cs_real_t (*)[3])dq,
                            alpha, kappa, hdg);
  else if (h_info.is_iso)
    _compute_cost_quant_iso(cm->n_fc, 1/cm->vol_c, cb->dpty_val,
                            (const cs_real_t (*)[3])pq,
                            (const cs_real_t (*)[3])dq,
                            alpha, kappa, hdg);
  else
    _compute_cost_quant(cm->n_fc, 1/cm->vol_c,
                        (const cs_real_3_t *)cb->dpty_mat,
                        (const cs_real_t (*)[3])pq,
                        (const cs_real_t (*)[3])dq,
                        alpha, kappa, hdg);

  _compute_hodge_cost(cm->n_fc, h_info.coef*h_info.coef, alpha, kappa,
                      hdg->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
    _check_vector_hodge((const cs_real_t (*)[3])pq,
                        (const cs_real_t (*)[3])dq,
                        hdg, h_info, cb);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from dual edges to primal faces.
 *          This function is related to face-based schemes
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_edfp_voro_get(const cs_param_hodge_t    h_info,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb)
{
  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EDFP);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_VORONOI);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ));

  /* Initialize the local matrix related to this discrete Hodge operator */
  cs_sdm_t  *hdg = cb->hdg;
  cs_sdm_square_init(cm->n_fc, hdg);

  for (short int f = 0; f < cm->n_fc; f++) {

    if (h_info.is_iso) {
      hdg->val[f*cm->n_fc+f] = cb->dpty_val*cm->face[f].meas/cm->dedge[f].meas;
    }
    else {

      const cs_quant_t  pfq = cm->face[f];
      const cs_real_3_t *tens = (const cs_real_3_t *)cb->dpty_mat;

      cs_real_3_t  mv;

      cs_math_33_3_product(tens, pfq.unitv, mv);
      hdg->val[f*cm->n_fc+f] = pfq.meas * _dp3(mv, pfq.unitv)/cm->edge[f].meas;

    }

  } // Loop on cell faces

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from dual edges to primal faces.
 *          This function is related to face-based schemes
 *
 * \param[in]      h_info    pointer to a cs_param_hodge_t structure
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_edfp_cost_get(const cs_param_hodge_t    h_info,
                       const cs_cell_mesh_t     *cm,
                       cs_cell_builder_t        *cb)
{
  /* Sanity check */
  assert(cb != NULL && cb->hdg != NULL);
  assert(h_info.type == CS_PARAM_HODGE_TYPE_EDFP);
  assert(h_info.algo == CS_PARAM_HODGE_ALGO_COST);
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ));

  /* Initialize the local matrix related to this discrete Hodge operator */
  cs_sdm_t  *hdg = cb->hdg;
  cs_sdm_square_init(cm->n_fc, hdg);

  double  *kappa = cb->values;
  double  *alpha = cb->values + cm->n_fc;
  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_fc;

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_nvec3_t  deq = cm->dedge[f];
    const cs_quant_t  pfq = cm->face[f];

    for (int k = 0; k < 3; k++) {
      dq[f][k] = deq.meas * deq.unitv[k];
      pq[f][k] = pfq.meas * pfq.unitv[k];
    }

  } /* Loop on cell faces */

  /* Compute additional geometrical quantities:
     Initialize the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  if (h_info.is_unity)
    _compute_cost_quant_iso(cm->n_fc, 1/cm->vol_c, 1.0,
                            (const cs_real_t (*)[3])dq,
                            (const cs_real_t (*)[3])pq,
                            alpha, kappa, hdg);
  else if (h_info.is_iso)
    _compute_cost_quant_iso(cm->n_fc, 1/cm->vol_c, cb->dpty_val,
                            (const cs_real_t (*)[3])dq,
                            (const cs_real_t (*)[3])pq,
                            alpha, kappa, hdg);
  else
    _compute_cost_quant(cm->n_fc, 1/cm->vol_c,
                        (const cs_real_3_t *)cb->dpty_mat,
                        (const cs_real_t (*)[3])dq,
                        (const cs_real_t (*)[3])pq,
                        alpha, kappa, hdg);

  _compute_hodge_cost(cm->n_fc, h_info.coef*h_info.coef, alpha, kappa,
                      hdg->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
    _check_vector_hodge((const cs_real_t (*)[3])dq,
                        (const cs_real_t (*)[3])pq,
                        hdg, h_info, cb);
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute cellwise a discrete hodge operator and multiple it with
 *          a vector
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      h_info    cs_param_hodge_t structure
 * \param[in]      pty       pointer to a cs_property_t structure or NULL
 * \param[in]      in_vals   vector to multiply with the discrete Hodge op.
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in, out] result    array storing the resulting matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_matvec(const cs_cdo_connect_t       *connect,
                const cs_cdo_quantities_t    *quant,
                const cs_param_hodge_t        h_info,
                const cs_property_t          *pty,
                const double                  in_vals[],
                cs_real_t                     t_eval,
                double                        result[])
{
  if (in_vals == NULL)
    return;
  if (result == NULL) {
    bft_error(__FILE__, __LINE__, 0, "Resulting vector must be allocated");
    return; // Avoid a warning
  }
  assert(connect != NULL && quant != NULL); // Sanity checks

#pragma omp parallel if (quant->n_cells > CS_THR_MIN) default(none)        \
  shared(quant, connect, in_vals, t_eval, result, pty)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    cs_flag_t  msh_flag = 0;
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    bool pty_uniform = cs_property_is_uniform(pty);
    cs_hodge_t  *compute = NULL;
    cs_cell_builder_t  *cb = NULL;
    double  *_in = NULL;

    switch (h_info.type) {

    case CS_PARAM_HODGE_TYPE_VPCD:

      msh_flag |= CS_CDO_LOCAL_PVQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOVB, connect);
      BFT_MALLOC(_in, connect->n_max_vbyc, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_vertices; i++) result[i] = 0;

      switch (h_info.algo) {
      case CS_PARAM_HODGE_ALGO_COST:
      case CS_PARAM_HODGE_ALGO_VORONOI:
        msh_flag |= CS_CDO_LOCAL_PVQ;
        compute = cs_hodge_vpcd_voro_get;
        break;
      case CS_PARAM_HODGE_ALGO_WBS:
        msh_flag |= CS_CDO_LOCAL_PVQ |CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_PFQ |
          CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_HFQ;
        compute = cs_hodge_vpcd_wbs_get;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid algorithm to build a VP->CD Hodge operator");
      }
      break;

    case CS_PARAM_HODGE_TYPE_EPFD:

      msh_flag |= CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOVB, connect);
      BFT_MALLOC(_in, connect->n_max_ebyc, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_edges; i++) result[i] = 0;

      switch (h_info.algo) {
      case CS_PARAM_HODGE_ALGO_COST:
        compute = cs_hodge_epfd_cost_get;
        break;
      case CS_PARAM_HODGE_ALGO_VORONOI:
        msh_flag |= CS_CDO_LOCAL_EFQ;
        compute = cs_hodge_epfd_voro_get;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid algorithm to build a EP->FD Hodge operator");
      }
      break;

    case CS_PARAM_HODGE_TYPE_EDFP:

      msh_flag |= CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOFB, connect);
      BFT_MALLOC(_in, connect->n_max_fbyc, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_faces; i++) result[i] = 0;

      switch (h_info.algo) {
      case CS_PARAM_HODGE_ALGO_COST:
        compute = cs_hodge_edfp_cost_get;
        break;
      case CS_PARAM_HODGE_ALGO_VORONOI:
        compute = cs_hodge_edfp_voro_get;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid algorithm to build a EP->FD Hodge operator");
      }
      break;

    case CS_PARAM_HODGE_TYPE_FPED:

      msh_flag |= CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOFB, connect);
      BFT_MALLOC(_in, connect->n_max_fbyc, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_faces; i++) result[i] = 0;

      switch (h_info.algo) {
      case CS_PARAM_HODGE_ALGO_COST:
        compute = cs_hodge_fped_cost_get;
        break;
      case CS_PARAM_HODGE_ALGO_VORONOI:
        compute = cs_hodge_fped_voro_get;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid algorithm to build a EP->FD Hodge operator");
      }
      break;

    case CS_PARAM_HODGE_TYPE_VC:

      msh_flag |= CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_HFQ |
         CS_CDO_LOCAL_DEQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_FEQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOVCB, connect);
      BFT_MALLOC(_in, connect->n_max_vbyc + 1, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_vertices + quant->n_cells; i++)
        result[i] = 0;

      switch (h_info.algo) {
      case CS_PARAM_HODGE_ALGO_WBS:
        compute = cs_hodge_vcb_wbs_get;
        break;
      default:
        bft_error(__FILE__, __LINE__, 0,
                  " Invalid algorithm to build a VP->CD Hodge operator");
      }
      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                "Invalid type of discrete Hodge operator");
    }

    if (pty == NULL) {
      cb->dpty_val = 1;
      cb->dpty_mat[0][0] = cb->dpty_mat[1][1] = cb->dpty_mat[2][2] = 1.0;
      cb->dpty_mat[0][1] = cb->dpty_mat[1][0] = cb->dpty_mat[2][0] = 0.0;
      cb->dpty_mat[0][2] = cb->dpty_mat[1][2] = cb->dpty_mat[2][1] = 0.0;
    }

    if (pty_uniform) { /* Get the value from the first cell */
      cs_property_get_cell_tensor(0, t_eval, pty, h_info.inv_pty, cb->dpty_mat);
      if (h_info.is_iso)
        cb->dpty_val = cb->dpty_mat[0][0];
    }

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */
      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Retrieve the value of the property inside the current cell */
      if (!pty_uniform) {
        cs_property_tensor_in_cell(cm, pty, t_eval, h_info.inv_pty,
                                   cb->dpty_mat);
        if (h_info.is_iso)
          cb->dpty_val = cb->dpty_mat[0][0];
      }

      /* Build the local discrete Hodge operator */
      compute(h_info, cm, cb);

      /* Define a local vector to multiply for the current cell */
      switch (h_info.type) {

      case CS_PARAM_HODGE_TYPE_VPCD:

        for (short int v = 0; v < cm->n_vc; v++)
          _in[v] = in_vals[cm->v_ids[v]];

        /* Local matrix-vector operation */
        cs_sdm_square_matvec(cb->hdg, _in, cb->values);

        /* Assemble the resulting vector */
        for (short int v = 0; v < cm->n_vc; v++)
#         pragma omp atomic
          result[cm->v_ids[v]] += cb->values[v];
        break;

      case CS_PARAM_HODGE_TYPE_EPFD:
        for (short int e = 0; e < cm->n_ec; e++)
          _in[e] = in_vals[cm->e_ids[e]];

        /* Local matrix-vector operation */
        cs_sdm_square_matvec(cb->hdg, _in, cb->values);

        /* Assemble the resulting vector */
        for (short int e = 0; e < cm->n_ec; e++)
#         pragma omp atomic
          result[cm->e_ids[e]] += cb->values[e];
        break;

      case CS_PARAM_HODGE_TYPE_FPED:
        for (short int f = 0; f < cm->n_fc; f++)
          _in[f] = in_vals[cm->f_ids[f]];

        /* Local matrix-vector operation */
        cs_sdm_square_matvec(cb->hdg, _in, cb->values);

        /* Assemble the resulting vector */
        for (short int f = 0; f < cm->n_fc; f++)
#         pragma omp atomic
          result[cm->f_ids[f]] += cb->values[f];
        break;

      case CS_PARAM_HODGE_TYPE_EDFP:
      default:
        bft_error(__FILE__, __LINE__, 0,
                  "Invalid type of discrete Hodge operator");

      } // Hodge type

    } // Main loop on cells

    BFT_FREE(_in);
    cs_cell_builder_free(&cb);

  } // OpenMP Block

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the hodge operator related to a face (i.e. a mass matrix
 *          with unity property) using a Whitney Barycentric Subdivision (WBS)
 *          algorithm
 *
 * \param[in]      fm        pointer to a cs_face_mesh_t structure
 * \param[in, out] hf        pointer to a cs_sdm_t structure to define
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_compute_wbs_surfacic(const cs_face_mesh_t    *fm,
                              cs_sdm_t                *hf)
{
  assert(hf != NULL && fm != NULL);

  /* Reset values */
  cs_sdm_square_init(fm->n_vf, hf);

  for (short int vfi = 0; vfi < fm->n_vf; vfi++) {

    double  *hi = hf->val + vfi*hf->n_rows;

    /* Default contribution */
    const double  default_coef = 0.5 * fm->wvf[vfi] * fm->face.meas;
    for (short int vfj = 0; vfj < fm->n_vf; vfj++)
      hi[vfj] = default_coef * fm->wvf[vfj];

    /* Specific diagonal contribution */
    hi[vfi] += 2 * default_coef * cs_math_1ov3;

  } /* Loop on face vertices */

  /* Specific extra-diag contribution */
  for (short int e = 0; e < fm->n_ef; e++) {

    const short int  v1 = fm->e2v_ids[2*e];
    const short int  v2 = fm->e2v_ids[2*e+1];
    const double  extra_val = cs_math_1ov12 * fm->tef[e];

    hf->val[v1*hf->n_rows + v2] += extra_val;
    hf->val[v2*hf->n_rows + v1] += extra_val;

  } /* Loop on face edges */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 2
  if (fm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Surfacic Hodge op.   ");
    cs_sdm_dump(fm->f_id, NULL, NULL, hf);
  }
#endif
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
