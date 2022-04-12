/*============================================================================
 * Manage discrete Hodge operators and closely related operators
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

static const char
cs_hodge_type_desc[CS_HODGE_N_TYPES][CS_BASE_STRING_LEN] =
  { N_("VpCd"),
    N_("EpFd"),
    N_("FpEd"),
    N_("EdFp"),
    N_("CpVd")  };

static const char
cs_hodge_algo_desc[CS_HODGE_N_ALGOS][CS_BASE_STRING_LEN] =
  { N_("Voronoi"),
    N_("Whitney on the Barycentric Subdivision (WBS)"),
    N_("Orthogonal Consistency/Stabilization (OCS)"),
    N_("Orthogonal Consistency/Sub-Stabilization (OCS2)"),
    N_("Orthogonal Consistency/Bubble-Stabilization (BUBBLE)"),
    N_("Automatic switch") };

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private function prototypes for debugging purpose
 *============================================================================*/

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
/*----------------------------------------------------------------------------*/
/*!
 * \brief   Check the coherency of the values of a stiffness matrix
 *
 * \param[in] c_id       current cell id
 * \param[in] sloc       pointer to a cs_sdm_t struct.
 */
/*----------------------------------------------------------------------------*/

static void
_check_stiffness(cs_lnum_t             c_id,
                 const cs_sdm_t       *sloc)
{
  if (c_id % CS_HODGE_MODULO != 0 || sloc == NULL)
    return;

  cs_log_printf(CS_LOG_DEFAULT, ">> Local stiffness matrix");
  cs_sdm_dump(c_id, NULL, NULL, sloc);

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
 * \param[in]      c_id    current cell id
 * \param[in]      vec     vectors of quantities to test against a hodge
 * \param[in]      res     vectors of quantities to compare with
 * \param[in]      hodge   pointer to a discrete Hodge operator
 * \param[in, out] cb      pointer to a cell builder structure
 *                         buffers to store temporary values
 */
/*----------------------------------------------------------------------------*/

static void
_check_vector_hodge(cs_lnum_t                c_id,
                    const cs_real_3_t       *vec,
                    const cs_real_3_t       *res,
                    cs_hodge_t              *hodge,
                    cs_cell_builder_t       *cb)

{
  if (c_id % CS_HODGE_MODULO != 0 || hmat == NULL)
    return;

  cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
  cs_sdm_dump(c_id, NULL, NULL, hodge->matrix);

  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;
  const cs_sdm_t  *hmat = hodge->matrix;

  cs_real_t  *in = cb->values;
  cs_real_t  *h_in = cb->values + hmat->n_rows;
  cs_real_t  *ref = cb->values + 2*hmat->n_rows;
  double  print_val = 0.;

  if (ptyd->is_unity) {
    ptyd->tensor[0][0] = ptyd->tensor[1][1] = ptyd->tensor[2][2] = 1;
    ptyd->tensor[0][1] = ptyd->tensor[1][0] = ptyd->tensor[2][0] = 0;
    ptyd->tensor[0][2] = ptyd->tensor[1][2] = ptyd->tensor[2][1] = 0;

  }
  else if (ptyd->is_iso) {
    ptyd->tensor[0][0] = ptyd->tensor[1][1] = ptyd->tensor[2][2] = ptyd->value;
    ptyd->tensor[0][1] = ptyd->tensor[1][0] = ptyd->tensor[2][0] = 0;
    ptyd->tensor[0][2] = ptyd->tensor[1][2] = ptyd->tensor[2][1] = 0;
  }

  const cs_real_3_t  a[3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };

  for (int dim = 0; dim < 3; dim++) {

    cs_real_3_t  pty_a;
    cs_math_33_3_product((const cs_real_t (*)[3])ptyd->tensor, a[dim], pty_a);

    for (int i = 0; i < hmat->n_rows; i++) {
      in[i] = vec[i][dim];
      ref[i] = _dp3(pty_a, res[i]);
    }

    cs_sdm_square_matvec(hmat, in, h_in);

    double  err = 0.;
    for (int i = 0; i < hmat->n_rows; i++)
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
#endif  /* DEBUG */

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute a L1-norm on a 3x3 tensor
 *
 * \param[in] tensor    3x3 matrix
 *
 * \return the value of the L1 norm
 */
/*----------------------------------------------------------------------------*/

inline static cs_real_t
_tensor_norm_l1(const cs_real_t   tensor[3][3])
{
  cs_real_t  _norm = 0;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      _norm += fabs(tensor[i][j]);

  return _norm;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Perform a part of a matrix-vector product
 *
 * \param[in]      n_ent   size of the matrix
 * \param[in]      i       starting index
 * \param[in]      dq_pq   matrix values
 * \param[in]      vec     vector
 * \param[in, out] mvec    resulting vector
 */
/*----------------------------------------------------------------------------*/

inline static void
_partial_matvec(const int         i,
                const cs_sdm_t   *dq_pq,
                const double     *restrict vec,
                double           *mvec)

{
  const int n_ent = dq_pq->n_rows;

  for (int irow = i; irow < n_ent; irow++) {
    const double  *restrict m_i = dq_pq->val + irow*n_ent;
    double s = 0;
    for (int j = 0; j < n_ent; j++) s += m_i[j]*vec[j];
    mvec[irow] = s;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the geometrical quantities for EpFd Hodge operators
 *          (Cost and bubble algorithm)
 *
 * \param[in]      cm    pointer to a cs_cell_mesh_t structure
 * \param[in, out] pq    primal vector-valued quantities
 * \param[in, out] dq    dual vector-valued quantities
 */
/*----------------------------------------------------------------------------*/

static inline void
_init_vb_geom_quant(const cs_cell_mesh_t    *cm,
                    cs_real_3_t             *pq,
                    cs_real_3_t             *dq)
{
  for (int ii = 0; ii < cm->n_ec; ii++) {

    cs_nvec3_t  dfq = cm->dface[ii];
    cs_quant_t  peq = cm->edge[ii];

    for (int k = 0; k < 3; k++) {
      dq[ii][k] = dfq.meas * dfq.unitv[k];
      pq[ii][k] = peq.meas * peq.unitv[k];
    }

  } /* Loop on cell edges */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the geometrical quantities for FpEd Hodge operators
 *          (Cost and bubble algorithm)
 *
 * \param[in]      cm    pointer to a cs_cell_mesh_t structure
 * \param[in, out] pq    primal vector-valued quantities
 * \param[in, out] dq    dual vector-valued quantities
 */
/*----------------------------------------------------------------------------*/

static inline void
_init_fb_geom_quant(const cs_cell_mesh_t    *cm,
                    cs_real_3_t             *pq,
                    cs_real_3_t             *dq)
{
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_nvec3_t  deq = cm->dedge[f];
    const cs_quant_t  pfq = cm->face[f];

    for (int k = 0; k < 3; k++) {
      dq[f][k] = deq.meas * deq.unitv[k];
      pq[f][k] = pfq.meas * pfq.unitv[k];
    }

  } /* Loop on cell faces */
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
    break;

  case CS_SPACE_SCHEME_CDOVCB:
    size = 2*n_vc + 3*n_ec + n_fc;
    BFT_MALLOC(cb->values, size, double);
    memset(cb->values, 0, size*sizeof(cs_real_t));

    size = 2*n_ec + n_vc;
    BFT_MALLOC(cb->vectors, size, cs_real_3_t);
    memset(cb->vectors, 0, size*sizeof(cs_real_3_t));
    break;

  case CS_SPACE_SCHEME_CDOFB:
    size = n_fc*(n_fc+1);
    BFT_MALLOC(cb->values, size, double);
    memset(cb->values, 0, size*sizeof(cs_real_t));

    size = 2*n_fc;
    BFT_MALLOC(cb->vectors, size, cs_real_3_t);
    memset(cb->vectors, 0, size*sizeof(cs_real_3_t));
    break;

  case CS_SPACE_SCHEME_CDOEB:
    {
      int  n_ent = CS_MAX(n_fc, n_ec);

      size = n_ent*(n_ent+1);
      BFT_MALLOC(cb->values, size, double);
      memset(cb->values, 0, size*sizeof(cs_real_t));

      size = 2*n_ent;
      BFT_MALLOC(cb->vectors, size, cs_real_3_t);
      memset(cb->vectors, 0, size*sizeof(cs_real_3_t));
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0, _("Invalid space scheme."));

  } /* End of switch on space scheme */

  return cb;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute quantities used for defining the entries of the discrete
 *          Hodge for COST algo. when the property is isotropic
 *          Initialize the local discrete Hodge op. with the consistency part
 *
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in]      hmat       pointer to a local Hodge matrix
 * \param[in, out] sloc       pointer to the local stiffness matrix
 */
/*----------------------------------------------------------------------------*/

static void
_define_vb_stiffness(const cs_cell_mesh_t   *cm,
                     const cs_sdm_t         *hmat,
                     cs_sdm_t               *sloc)
{
  /* Initialize the local stiffness matrix */

  cs_sdm_square_init(cm->n_vc, sloc);

  for (int ei = 0; ei < cm->n_ec; ei++) { /* Loop on cell edges I */

    const short int  *v_sgn = cm->e2v_sgn + ei;
    const short int  *v_ids = cm->e2v_ids + 2*ei;
    const short int  i1 = v_ids[0], i2 = v_ids[1];
    assert(i1 < i2);

    double  *si1 = sloc->val + i1*sloc->n_rows;
    double  *si2 = sloc->val + i2*sloc->n_rows;

    /* Diagonal value: consistency part has already been computed */

    const double  *hii = hmat->val + ei*(1 + cm->n_ec);
    const double  dval = hii[0];

    si1[i1] += dval;
    si1[i2] -= dval;
    si2[i2] += dval;

    /* Compute extra-diag entries */

    for (int _j = 1; _j < cm->n_ec-ei; _j++) { /* Loop on cell entities J */

      const short int  j1 = v_ids[2*_j], j2 = v_ids[2*_j+1];
      assert(j1 < j2);

      double  *sj1 = sloc->val + j1*sloc->n_rows;
      double  *sj2 = sloc->val + j2*sloc->n_rows;

      /* Add contribution from the stabilization part for each sub-volume
         related to a primal entity */

      const double xval = hii[_j] * v_sgn[0]*v_sgn[_j];

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
        else                    /* j1 < i1 < i2 = j2 */
          sj1[i1] += xval;

        /* Remark: the case i1 == j1 is not possible since ei != ej */
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

  cs_sdm_symm_ur(sloc);
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
 * \param[in, out] hmat       pointer to a cs_sdm_t struct.
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
                        cs_sdm_t               *hmat)
{
  /* Compute several useful quantities
     alpha_ij = delta_ij - pq_j.Consist_i where Consist_i = 1/|c| dq_i
     qmq_ii = dq_i.ptyval.dq_i
     kappa_i = qmq_ii / |subvol_i|
  */

  for (int i = 0; i < n_ent; i++) {

    const double  dsvol_i = _dp3(dq[i], pq[i]);

    double  *alpha_i = alpha + i*n_ent;
    double  *mi = hmat->val + i*n_ent;

    alpha_i[i] = 1 - invcvol * dsvol_i;

    const double  qmq_ii = ptyval * _dp3(dq[i], dq[i]);

    mi[i] = invcvol * qmq_ii;
    kappa[i] = 3. * qmq_ii / dsvol_i;

    for (int j = i+1; j < n_ent; j++) {

      /* Initialize the upper right part of hmat with the consistency part */

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
 * \param[in, out] hmat     pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_cost_quant(const int               n_ent,
                    const double            ovc,
                    const cs_real_t         pty[3][3],
                    const cs_real_3_t      *pq,
                    const cs_real_3_t      *dq,
                    double                 *alpha,
                    double                 *kappa,
                    cs_sdm_t               *hmat)
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

    double  *h_i = hmat->val + i*n_ent;

    h_i[i] = _dp3(dq[i], mdq_i); /* Add the consistent part */

    double  *alpha_i = alpha + i*n_ent;

    alpha_i[i] = _dp3(dq[i], pq[i]);
    kappa[i] = 3. * h_i[i] / alpha_i[i];
    alpha_i[i] = 1 - ovc*alpha_i[i];
    h_i[i] *= ovc;

    for (int j = i+1; j < n_ent; j++) {

      /* Initialize the upper right part of hmat with the consistency part */

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
 * \brief   Compute the discrete EpFd Hodge operator (the upper right part).
 *          Co+St algo. in case of isotropic material property.
 *
 * \param[in]      n_ent    number of local entities
 * \param[in]      dbeta2   space dim * squared value of the stabilization coef.
 * \param[in]      ovc      reciprocal of the cell volume
 * \param[in]      pty      values of the the material pty in this cell
 * \param[in]      pq       pointer to the first set of quantities
 * \param[in]      dq       pointer to the second set of quantities
 * \param[in, out] cb       temporary buffers
 * \param[in, out] hmat     pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_iso_hodge_ur(const int               n_ent,
                      const double            dbeta2,
                      const double            ovc,
                      const cs_real_t         pty,
                      const cs_real_3_t      *pq,
                      const cs_real_3_t      *dq,
                      cs_cell_builder_t      *cb,
                      cs_sdm_t               *hmat)
{
  const double  ptyc = pty*ovc;

  double  *kappa = cb->values;                /* size = n_ent */
  double  *kappa_pq_dqi = cb->values + n_ent; /* size = n_ent */
  double  *stab = cb->values + 2*n_ent;       /* size = n_ent */
  double  *dq_pq = cb->aux->val;              /* size = n_ent*n_ent */

  cs_sdm_square_init(n_ent, cb->aux);

  const double  dbetac = dbeta2*ovc;
  const double  dbetac2 = dbetac*ovc;
  const double  beta_coef = ovc * (1 - 2*dbeta2);

  /* Initialize the upper right part of the discrete Hodge op and store useful
     quantities */

  for (int i = 0; i < n_ent; i++) {

    const double  dqi[3] = {dq[i][0], dq[i][1], dq[i][2]};

    double  *dqi_pq = dq_pq + i*n_ent;
    for (int j = 0; j < n_ent; j++)
      dqi_pq[j] = _dp3(dqi, pq[j]);

    const double  dqi_m_dqi = pty * _dp3(dqi, dqi);
    kappa[i] = dqi_m_dqi/dqi_pq[i];

    double  *hi = hmat->val + i*n_ent;

    hi[i] = dqi_m_dqi*beta_coef + dbeta2*kappa[i];
    for (int j = i+1; j < n_ent; j++)
      hi[j] = ptyc * _dp3(dq[j], dqi);

  }

  /* Build the upper right part of the discrete Hodge operator */

  for (int i = 0; i < n_ent; i++) {

    const double  *dqi_pq = dq_pq + n_ent*i;

    for (int k = 0; k < n_ent;k++)
      kappa_pq_dqi[k] = kappa[k] * dqi_pq[k];

    _partial_matvec(i, cb->aux, kappa_pq_dqi, stab);

    double  *hi = hmat->val + i*n_ent;

    hi[i] += dbetac2 * stab[i];

    /* Compute the extra-diagonal terms */

    const double  kai = kappa[i];
    for (int j = i+1; j < n_ent; j++) {
      double  contrib = ovc * stab[j];
      contrib -= kai*dq_pq[j*n_ent+i] + kappa[j]*dqi_pq[j];
      contrib *= dbetac;
      hi[j] += contrib;
    }

  } /* Loop on rows (entities i) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the discrete EpFd Hodge operator (the upper right part).
 *          Co+St algo. in case of anisotropic material property.
 *
 * \param[in]      n_ent    number of local entities
 * \param[in]      dbeta2   space dim * squared value of the stabilization coef.
 * \param[in]      vol_c    cell volume
 * \param[in]      pty      values of the tensor related to the material pty
 * \param[in]      pq       pointer to the first set of quantities
 * \param[in]      dq       pointer to the second set of quantities
 * \param[in, out] cb       temporary buffers
 * \param[in, out] hmat     pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_aniso_hodge_ur(const int               n_ent,
                        const double            dbeta2,
                        const double            ovc,
                        const cs_real_t         pty[3][3],
                        const cs_real_3_t      *pq,
                        const cs_real_3_t      *dq,
                        cs_cell_builder_t      *cb,
                        cs_sdm_t               *hmat)
{
  double  *kappa = cb->values;                /* size = n_ent */
  double  *kappa_pq_dqi = cb->values + n_ent; /* size = n_ent */
  double  *stab = cb->values + 2*n_ent;       /* size = n_ent */
  double  *dq_pq = cb->aux->val;              /* size = n_ent*n_ent */

  cs_sdm_square_init(n_ent, cb->aux);

  const double  dbetac = dbeta2*ovc;
  const double  dbetac2 = dbetac*ovc;
  const double  beta_coef = ovc * (1 - 2*dbeta2);

  /* Initialize the upper right part of the discrete Hodge op and store useful
     quantities */

  for (int i = 0; i < n_ent; i++) {

    const double  dqi[3] = {dq[i][0], dq[i][1], dq[i][2]};

    double  *dqi_pq = dq_pq + i*n_ent;
    for (int j = 0; j < n_ent; j++)
      dqi_pq[j] = _dp3(dqi, pq[j]);

    const double  mdqi[3]
      = { pty[0][0] * dqi[0] + pty[0][1] * dqi[1] + pty[0][2] * dqi[2],
          pty[1][0] * dqi[0] + pty[1][1] * dqi[1] + pty[1][2] * dqi[2],
          pty[2][0] * dqi[0] + pty[2][1] * dqi[1] + pty[2][2] * dqi[2] };

    const double  dqi_m_dqi = _dp3(dqi, mdqi);

    kappa[i] = dqi_m_dqi/dqi_pq[i];

    double  *hi = hmat->val + i*n_ent;

    hi[i] = dqi_m_dqi*beta_coef + dbeta2*kappa[i];
    for (int j = i+1; j < n_ent; j++)
      hi[j] = ovc * _dp3(dq[j], mdqi);

  }

  /* Build the upper right part of the discrete Hodge operator */

  for (int i = 0; i < n_ent; i++) {

    const double  *dqi_pq = dq_pq + n_ent*i;

    for (int k = 0; k < n_ent; k++)
      kappa_pq_dqi[k] = kappa[k] * dqi_pq[k];

    _partial_matvec(i, cb->aux, kappa_pq_dqi, stab);

    double  *hi = hmat->val + i*n_ent;

    hi[i] += dbetac2 * stab[i];

    /* Compute the extra-diagonal terms */

    const double  kai = kappa[i];
    for (int j = i+1; j < n_ent; j++) {
      double  contrib = ovc*stab[j];
      contrib -= kai*dq_pq[j*n_ent+i] + kappa[j]*dqi_pq[j];
      contrib *= dbetac;
      hi[j] += contrib;
    }

  } /* Loop on rows (entities i) */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the discrete Hodge operator (the upper right part).
 *          Co+St algo. with bubble stabilization in case of isotropic
 *          material property.
 *
 * \param[in]      n_ent    number of local entities
 * \param[in]      beta     the stabilization coef.
 * \param[in]      ovc      reciprocal of the cell volume
 * \param[in]      pty_val  value of the property inside this cell
 * \param[in]      pq       pointer to the first set of quantities
 * \param[in]      dq       pointer to the second set of quantities
 * \param[in, out] cb       temporary buffers
 * \param[in, out] hmat     pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_iso_bubble_hodge_ur(const int               n_ent,
                             const double            beta,
                             const double            ovc,
                             const cs_real_t         pty_val,
                             const cs_real_3_t      *pq,
                             const cs_real_3_t      *dq,
                             cs_cell_builder_t      *cb,
                             cs_sdm_t               *hmat)
{
  double  *kappa = cb->values;              /* size = n_ent */
  double  *dq_pq = cb->aux->val;            /* size = n_ent*n_ent */

  /* Initialize the upper right part of the discrete Hodge op and store useful
     quantities */

  for (int i = 0; i < n_ent; i++) {

    const double  dqi[3] = {dq[i][0], dq[i][1], dq[i][2]};
    const double  dqi_m_dqi = pty_val * _dp3(dqi, dqi);
    const double  dqi_pqi = _dp3(dqi, pq[i]);

    kappa[i] = dqi_m_dqi/dqi_pqi;

    double  *dqi_pq = dq_pq + i*n_ent;
    for (int j = 0; j < i; j++)
      dqi_pq[j] = _dp3(dqi, pq[j])*ovc;
    dqi_pq[i] = dqi_pqi*ovc;
    for (int j = i+1; j < n_ent; j++)
      dqi_pq[j] = _dp3(dqi, pq[j])*ovc;

    /* Consistent part */

    double  *hi = hmat->val + i*n_ent;
    hi[i] = dqi_m_dqi*ovc;
    const double  coef = ovc * pty_val;
    for (int j = i+1; j < n_ent; j++)
      hi[j] = coef * _dp3(dq[j], dqi);

  }

  /* Add the stabilization part */

  /* \int_{p_{ec}} \theta_e*\theta_e = 0.1*|p_{ec}| and d=3 */

  const double  stab_coef = 0.3*beta*beta;

  /* Build the upper right part of the discrete Hodge operator */

  for (int i = 0; i < n_ent; i++) {

    const double  *dqi_pq = dq_pq + n_ent*i;
    double  *hi = hmat->val + i*n_ent;

    /* Diagonal term */

    double stab = 0;
    for (int k = 0; k < n_ent; k++) {
      const double  a_ik = (i == k) ? 1 - dqi_pq[k] : -dqi_pq[k];
      stab += kappa[k] * a_ik * a_ik;
    }
    hi[i] += stab_coef * stab;

    /* Extra-diag term */

    for (int j = i+1; j < n_ent; j++) {

      const double  *dqj_pq = dq_pq + n_ent*j;

      stab = 0;
      for (int k = 0; k < n_ent; k++) {
        const double  a_ik = (i == k) ? 1 - dqi_pq[k] : -dqi_pq[k];
        const double  a_jk = (j == k) ? 1 - dqj_pq[k] : -dqj_pq[k];
        stab += kappa[k] * a_ik * a_jk;
      }
      hi[j] += stab_coef * stab;

    } /* Loop on row elements */

  } /* Loop on partition elements */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the discrete Hodge operator (the upper right part).
 *          Co+St algo. with bubble stabilization in case of anisotropic
 *          material property.
 *
 * \param[in]      n_ent    number of local entities
 * \param[in]      beta     the stabilization coef.
 * \param[in]      ovc      reciprocal of the cell volume
 * \param[in]      pty      value of the property inside this cell
 * \param[in]      pq       pointer to the first set of quantities
 * \param[in]      dq       pointer to the second set of quantities
 * \param[in, out] cb       temporary buffers
 * \param[in, out] hmat     pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_aniso_bubble_hodge_ur(const int               n_ent,
                               const double            beta,
                               const double            ovc,
                               const cs_real_t         pty[3][3],
                               const cs_real_3_t      *pq,
                               const cs_real_3_t      *dq,
                               cs_cell_builder_t      *cb,
                               cs_sdm_t               *hmat)
{
  double  *kappa = cb->values;              /* size = n_ent */
  double  *dq_pq = cb->aux->val;            /* size = n_ent*n_ent */

  /* Initialize the upper right part of the discrete Hodge op and store useful
     quantities */

  for (int i = 0; i < n_ent; i++) {

    const double  dqi[3] = {dq[i][0], dq[i][1], dq[i][2]};
    const double  mdqi[3]
      = { pty[0][0] * dqi[0] + pty[0][1] * dqi[1] + pty[0][2] * dqi[2],
          pty[1][0] * dqi[0] + pty[1][1] * dqi[1] + pty[1][2] * dqi[2],
          pty[2][0] * dqi[0] + pty[2][1] * dqi[1] + pty[2][2] * dqi[2] };
    const double  dqi_m_dqi = _dp3(dqi, mdqi);
    const double  dqi_pqi = _dp3(dqi, pq[i]);

    kappa[i] = dqi_m_dqi/dqi_pqi;

    double  *dqi_pq = dq_pq + i*n_ent;
    for (int j = 0; j < i; j++)
      dqi_pq[j] = _dp3(dqi, pq[j])*ovc;
    dqi_pq[i] = dqi_pqi*ovc;
    for (int j = i+1; j < n_ent; j++)
      dqi_pq[j] = _dp3(dqi, pq[j])*ovc;

    /* Consistent part */

    double  *hi = hmat->val + i*n_ent;
    hi[i] = dqi_m_dqi*ovc;
    for (int j = i+1; j < n_ent; j++)
      hi[j] = ovc * _dp3(dq[j], mdqi);

  }

  /* Add the stabilization part */

  /* \int_{p_{ec}} \theta_e*\theta_e = 0.1*|p_{ec}| and d=3 */

  const double  stab_coef = 0.3*beta*beta;

  /* Build the upper right part of the discrete Hodge operator */

  for (int i = 0; i < n_ent; i++) {

    const double  *dqi_pq = dq_pq + n_ent*i;
    double  *hi = hmat->val + i*n_ent;

    /* Diagonal term */
    double stab = 0;
    for (int k = 0; k < n_ent; k++) {
      const double  a_ik = (i == k) ? 1 - dqi_pq[k] : -dqi_pq[k];
      stab += kappa[k] * a_ik * a_ik;
    }
    hi[i] += stab_coef * stab;

    /* Extra-diag term */

    for (int j = i+1; j < n_ent; j++) {

      const double  *dqj_pq = dq_pq + n_ent*j;

      stab = 0;
      for (int k = 0; k < n_ent; k++) {
        const double  a_ik = (i == k) ? 1 - dqi_pq[k] : -dqi_pq[k];
        const double  a_jk = (j == k) ? 1 - dqj_pq[k] : -dqj_pq[k];
        stab += kappa[k] * a_ik * a_jk;
      }
      hi[j] += stab_coef * stab;

    } /* Loop on row elements */

  } /* Loop on partition elements */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the discrete EpFd Hodge operator (the upper right part)
 *          from primal edges to dual faces with the algorithm called:
 *          Orthogonal Consistent/Sub-Stabilization decomposition (OCS2) with a
 *          subdivision of pvol_{e,c}
 *          Case of anisotropic material property.
 *
 * \param[in]      dbeta2   space dim * squared value of the stabilization coef.
 * \param[in]      pty      values of the tensor related to the material pty
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb       temporary buffers
 * \param[in, out] hmat     pointer to a cs_sdm_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_compute_aniso_hepfd_ocs2_ur(const double            dbeta2,
                             const cs_real_t         pty[3][3],
                             const cs_cell_mesh_t   *cm,
                             cs_cell_builder_t      *cb,
                             cs_sdm_t               *hmat)
{
  const double  ovc = 1./cm->vol_c;
  const int  n_ent = cm->n_ec;

  /* Store the consistent part of the reconstruction of the basis element */

  cs_real_3_t  *consist = cb->vectors;
  for (int i = 0; i < n_ent; i++) {
    const  double  fd_coef = ovc * cm->dface[i].meas;
    for (int k = 0; k < 3; k++)
      consist[i][k] = fd_coef * cm->dface[i].unitv[k];
  }

  /* Initialize the upper right part of the discrete Hodge op and store useful
     quantities */

  /* Consistency part */

  for (int i = 0; i < n_ent; i++) {

    double  pty_fd_i[3] = {0, 0, 0};
    for (int k = 0; k < 3; k++) {
      pty_fd_i[0] += pty[0][k] * cm->dface[i].unitv[k];
      pty_fd_i[1] += pty[1][k] * cm->dface[i].unitv[k];
      pty_fd_i[2] += pty[2][k] * cm->dface[i].unitv[k];
    }
    for (int k = 0; k < 3; k++)
      pty_fd_i[k] *= cm->dface[i].meas;

    double  *h_i = hmat->val + i*n_ent;
    for (int j = i; j < n_ent; j++)
      h_i[j] = _dp3(consist[j], pty_fd_i);

  }

  /* Compute the contribution part of each edge for the stabilization term */

  cs_real_t  *contrib_pe = cb->values;
  memset(contrib_pe, 0, n_ent*sizeof(cs_real_t));

  for (int f = 0; f < cm->n_fc; f++) {

    for (int i = cm->f2e_idx[f]; f < cm->f2e_idx[f+1]; f++) {

      const short int  k = cm->f2e_ids[i]; /* edge k */
      const cs_real_t  ep_k[3] = { cm->edge[k].meas * cm->edge[k].unitv[0],
                                   cm->edge[k].meas * cm->edge[k].unitv[1],
                                   cm->edge[k].meas * cm->edge[k].unitv[2] };

      cs_real_3_t  pty_fdk = {0, 0, 0};
      for (int kk = 0; kk < 3; kk++) {
        pty_fdk[0] += pty[0][kk] * cm->sefc[i].unitv[kk];
        pty_fdk[1] += pty[1][kk] * cm->sefc[i].unitv[kk];
        pty_fdk[2] += pty[2][kk] * cm->sefc[i].unitv[kk];
      }

      const cs_real_t  coef = cm->sefc[i].meas * dbeta2;
      contrib_pe[k] +=
        coef * _dp3(cm->sefc[i].unitv, pty_fdk)/_dp3(cm->sefc[i].unitv, ep_k) ;

    } /* Loop on face edges */

  } /* Loop on cell faces */

  /* Stabilization part */

  for (int k = 0; k < n_ent; k++) {

    const cs_real_t  contrib_pek = contrib_pe[k];
    const cs_real_t  ep_k[3] = { cm->edge[k].meas * cm->edge[k].unitv[0],
                                 cm->edge[k].meas * cm->edge[k].unitv[1],
                                 cm->edge[k].meas * cm->edge[k].unitv[2] };


    for (int i = 0; i < n_ent; i++) {

      double  contrib_i = -_dp3(consist[i], ep_k);
      if (i == k) contrib_i += 1;
      const double  contrib_ik = contrib_pek * contrib_i;

      double  *h_i = hmat->val + i*n_ent;
      h_i[i] += contrib_ik * contrib_i;

      for (int j = i+1; j < n_ent; j++) {
        if (j != k)
          h_i[j] += - _dp3(consist[j], ep_k) * contrib_ik;
        else
          h_i[j] += (1 -_dp3(consist[j], ep_k)) * contrib_ik;
      }

    } /* i */

  } /* Loop on sub-volume k */
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

    mi[i] += beta2 * stab_part; /* Consistency part has already been computed */

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
      mj[i] = mi[j]; /* Symmetric by construction */

    } /* End of loop on J entities */

  } /* End of loop on I entities */
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize a pointer to a cs_hodge_t structure
 *
 * \param[in] connect        pointer to cs_cdo_connect_t structure
 * \param[in] property       pointer to a property structure
 * \param[in] hp             pointer to a cs_hodge_param_t structure
 * \param[in] need_tensor    true if one needs a tensor otherwise false
 * \param[in] need_eigen     true if one needs to compute eigen valuese
 *
 * \return a pointer to the new allocated cs_hodge_t structure
 */
/*----------------------------------------------------------------------------*/

cs_hodge_t *
cs_hodge_create(const cs_cdo_connect_t   *connect,
                const cs_property_t      *property,
                const cs_hodge_param_t   *hp,
                bool                      need_tensor,
                bool                      need_eigen)
{
  cs_hodge_t  *hdg = NULL;

  BFT_MALLOC(hdg, 1, cs_hodge_t);

  hdg->param = hp;

  switch (hp->type) {

  case CS_HODGE_TYPE_VPCD:
    hdg->matrix = cs_sdm_square_create(connect->n_max_vbyc);
    break;
  case CS_HODGE_TYPE_EPFD:
    hdg->matrix = cs_sdm_square_create(connect->n_max_ebyc);
    break;
  case CS_HODGE_TYPE_FPED:
  case CS_HODGE_TYPE_EDFP:
    hdg->matrix = cs_sdm_square_create(connect->n_max_fbyc);
    break;
  case CS_HODGE_TYPE_CPVD:
    hdg->matrix = cs_sdm_square_create(1);
    break;
  case CS_HODGE_TYPE_FB:
    hdg->matrix = cs_sdm_square_create(connect->n_max_fbyc + 1);
    break;
  case CS_HODGE_TYPE_VC:
    hdg->matrix = cs_sdm_square_create(connect->n_max_vbyc + 1);
    break;

  default:
    hdg->matrix = NULL;
    break;
  }

  BFT_MALLOC(hdg->pty_data, 1, cs_property_data_t);
  cs_property_data_init(need_tensor, need_eigen, property, hdg->pty_data);
  if (hdg->pty_data->is_unity == false && connect->n_cells > 0)
    cs_hodge_set_property_value(0, 0, 0, hdg);

  return hdg;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create and initialize an array of pointers to a cs_hodge_t
 *         structures. This array is of size the number of OpenMP threads.
 *         Only the one associated to the current thread is set.
 *
 * \param[in] connect        pointer to cs_cdo_connect_t structure
 * \param[in] property       pointer to a property structure
 * \param[in] hp             pointer to a cs_hodge_param_t structure
 * \param[in] need_tensor    true if one needs a tensor otherwise false
 * \param[in] need_eigen     true if one needs to compute eigen valuese
 *
 * \return an array of pointers of cs_hodge_t structures
 */
/*----------------------------------------------------------------------------*/

cs_hodge_t **
cs_hodge_init_context(const cs_cdo_connect_t   *connect,
                      const cs_property_t      *property,
                      const cs_hodge_param_t   *hp,
                      bool                      need_tensor,
                      bool                      need_eigen)
{
  cs_hodge_t  **hodge_array = NULL;

  BFT_MALLOC(hodge_array, cs_glob_n_threads, cs_hodge_t *);
  for (int i = 0; i < cs_glob_n_threads; i++)
    hodge_array[i] = NULL;

#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
#pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    assert(t_id < cs_glob_n_threads);

    hodge_array[t_id] = cs_hodge_create(connect, property, hp,
                                        need_tensor, need_eigen);
  }
#else
  assert(cs_glob_n_threads == 1);
  hodge_array[0] = cs_hodge_create(connect, property, hp,
                                   need_tensor, need_eigen);
#endif /* openMP */

  return hodge_array;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_hodge_t structure
 *
 * \param[in, out] p_hodge    double pointer to a cs_hodge_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_free(cs_hodge_t    **p_hodge)
{
  cs_hodge_t  *hdg = *p_hodge;

  if (hdg == NULL)
    return;

  /* Other pointers are shared */

  hdg->matrix = cs_sdm_free(hdg->matrix);

  BFT_FREE(hdg->pty_data);

  BFT_FREE(hdg);
  *p_hodge = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a set of cs_hodge_t structures
 *
 * \param[in, out] p_hodges    triple pointer to cs_hodge_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_free_context(cs_hodge_t    ***p_hodges)
{
  cs_hodge_t  **hodge_array = *p_hodges;

  if (hodge_array == NULL)
    return;

#if defined(HAVE_OPENMP)
# pragma omp parallel
  {
    int t_id = omp_get_thread_num();
    cs_hodge_free(&(hodge_array[t_id]));
  }
#else  /* No OpenMP */
  cs_hodge_free(&(hodge_array[0]));
#endif /* openMP */

  BFT_FREE(hodge_array);
  *p_hodges = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a function pointer to compute a discrete Hodge operator
 *
 * \param[in] calling_func    name of the calling function
 * \param[in] hp              a cs_hodge_param_t structure
 *
 * \return a pointer to the corresponding function
 */
/*----------------------------------------------------------------------------*/

cs_hodge_compute_t *
cs_hodge_get_func(const char               *calling_func,
                  const cs_hodge_param_t    hp)
{
  cs_hodge_compute_t  *hdg_func = NULL;

  switch (hp.type) {

  case CS_HODGE_TYPE_VPCD:
    {
      switch (hp.algo) {

      case CS_HODGE_ALGO_COST:
      case CS_HODGE_ALGO_OCS2:
      case CS_HODGE_ALGO_BUBBLE:
      case CS_HODGE_ALGO_VORONOI:
        return cs_hodge_vpcd_voro_get;
      case CS_HODGE_ALGO_WBS:
        return cs_hodge_vpcd_wbs_get;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid algorithm to compute a Vp-Cd Hodge operator.\n"
                  " The calling function is %s\n", __func__, calling_func);
        break;
      }
    }
    break;

  case CS_HODGE_TYPE_EPFD:
    {
      switch (hp.algo) {

      case CS_HODGE_ALGO_COST:
        return cs_hodge_epfd_cost_get;
      case CS_HODGE_ALGO_OCS2:
        return cs_hodge_epfd_ocs2_get;
      case CS_HODGE_ALGO_BUBBLE:
      case CS_HODGE_ALGO_WBS:   /* By default, one should define a specific
                                   algorithm for WBS */
        return cs_hodge_epfd_bubble_get;
      case CS_HODGE_ALGO_VORONOI:
        return cs_hodge_epfd_voro_get;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid algorithm to compute a Ep-Fd Hodge operator.\n"
                  " The calling function is %s\n", __func__, calling_func);
        break;
      }
    }
    break;

  case CS_HODGE_TYPE_EDFP:
    {
      switch (hp.algo) {

      case CS_HODGE_ALGO_COST:
        return cs_hodge_edfp_cost_get_opt;
      case CS_HODGE_ALGO_BUBBLE:
      case CS_HODGE_ALGO_WBS:   /* By default, one should define a specific
                                   algorithm for WBS */
        return cs_hodge_edfp_bubble_get;
      case CS_HODGE_ALGO_VORONOI:
        return cs_hodge_edfp_voro_get;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid algorithm to compute a Ed-Fp Hodge operator.\n"
                  " The calling function is %s\n", __func__, calling_func);
        break;
      }
    }
    break;

  case CS_HODGE_TYPE_FPED:
    {
      switch (hp.algo) {

      case CS_HODGE_ALGO_COST:
        return cs_hodge_fped_cost_get;
      case CS_HODGE_ALGO_BUBBLE:
      case CS_HODGE_ALGO_WBS:   /* By default, one should define a specific
                                   algorithm for WBS */
        return cs_hodge_fped_bubble_get;
      case CS_HODGE_ALGO_VORONOI:
        return cs_hodge_fped_voro_get;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  "%s: Invalid algorithm to compute a Fp-Ed Hodge operator.\n"
                  " The calling function is %s\n", __func__, calling_func);
        break;
      }
    }
    break;

  case CS_HODGE_TYPE_FB:
    return cs_hodge_fb_get;

  case CS_HODGE_TYPE_VC:
    switch (hp.algo) {

    case CS_HODGE_ALGO_VORONOI:
      return cs_hodge_vcb_voro_get;
    case CS_HODGE_ALGO_WBS:
      return cs_hodge_vcb_wbs_get;
    default:
      bft_error(__FILE__, __LINE__, 0,
                "%s: Invalid algorithm to compute a Fp-Ed Hodge operator.\n"
                " The calling function is %s\n", __func__, calling_func);
      break;
    }
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              " %s: Invalid type of Hodge operator called by %s\n",
              __func__, calling_func);
    break;

  } /* Switch on the type of Hodge operator */

  return hdg_func;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check the consistency of the settings between terms related to a
 *         mass matrix and define the common algorithm to use.
 *         If a term should not be considered, set the algorithm to
 *         CS_HODGE_N_ALGOS
 *
 * \param[in] eqname     name of the equation to check
 * \param[in] reac_algo  optional algo. used for the reaction term
 * \param[in] time_algo  optional algo. used for the unsteady term
 * \param[in] srct_algo  optional algo. used for the source term
 *
 * \return the common algorithm to use
 */
/*----------------------------------------------------------------------------*/

cs_hodge_algo_t
cs_hodge_set_mass_algo(const char         *eqname,
                       cs_hodge_algo_t     reac_algo,
                       cs_hodge_algo_t     time_algo,
                       cs_hodge_algo_t     srct_algo)
{
  cs_hodge_algo_t  return_algo = CS_HODGE_ALGO_VORONOI;

  if (reac_algo != CS_HODGE_N_ALGOS) { /* Hodge algo. is set for reaction */

    return_algo = reac_algo;

    if (time_algo != CS_HODGE_N_ALGOS) {

      if (reac_algo != time_algo)
        bft_error(__FILE__, __LINE__, 0,
                  " %s: The configuration of the Hodge algorithm between the"
                  " reaction and unsteady term is not consistent.\n"
                  " Please check your settings for equation \"%s\"\n",
                  __func__, eqname);

      if (srct_algo != CS_HODGE_N_ALGOS)
        if (time_algo != srct_algo)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: The configuration of the Hodge algorithm between the"
                    " source term and unsteady term is not consistent.\n"
                    " Please check your settings for equation \"%s\"\n",
                    __func__, eqname);

    }
    else { /* Hodge algo not set for the unsteady term */

      if (srct_algo != CS_HODGE_N_ALGOS)
        if (reac_algo != srct_algo)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: The configuration of the Hodge algorithm between the"
                    " reaction and source term is not consistent.\n"
                    " Please check your settings for equation \"%s\"\n",
                    __func__, eqname);

    }

  }
  else { /* Hodge algo not set for the reaction term */

    if (time_algo != CS_HODGE_N_ALGOS) {

      return_algo = time_algo;

      if (srct_algo != CS_HODGE_N_ALGOS)
        if (time_algo != srct_algo)
          bft_error(__FILE__, __LINE__, 0,
                    " %s: The configuration of the Hodge algorithm between the"
                    " source term and unsteady term is not consistent.\n"
                    " Please check your settings for equation \"%s\"\n",
                    __func__, eqname);

    }
    else { /* Neither time_algo nor reac_algo is set */

      if (srct_algo != CS_HODGE_N_ALGOS)
        return_algo = srct_algo;

    }

  }

  return return_algo;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Output the settings related to a cs_hodge_param_t structure
 *
 * \param[in] prefix    optional string
 * \param[in] property  optional pointer to a property structure
 * \param[in] hp        a cs_hodge_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_param_log(const char               *prefix,
                   const cs_property_t      *property,
                   const cs_hodge_param_t    hp)
{
  const char  *_p;
  const char _empty_prefix[2] = "";
  if (prefix == NULL)
    _p = _empty_prefix;
  else
    _p = prefix;

  cs_log_printf(CS_LOG_SETUP, "%s | Type: %s\n",
                _p, cs_hodge_type_desc[hp.type]);
  cs_log_printf(CS_LOG_SETUP, "%s | Algo: %s\n",
                _p, cs_hodge_algo_desc[hp.algo]);
  if (hp.algo == CS_HODGE_ALGO_COST ||
      hp.algo == CS_HODGE_ALGO_OCS2 ||
      hp.algo == CS_HODGE_ALGO_BUBBLE)
    cs_log_printf(CS_LOG_SETUP, "%s | Algo.Coef: %.3e\n",
                  _p, hp.coef);

  if (property != NULL)
    cs_log_printf(CS_LOG_SETUP, "%s | Associated property: %s\n",
                  _p, cs_property_get_name(property));
  cs_log_printf(CS_LOG_SETUP, "%s | Property inversion: %s\n",
                _p, cs_base_strtf(hp.inv_pty));
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy the set of parameters associated to a discrete Hodge operator
 *         to another one
 *
 * \param[in]       h_ref   reference set of parameters
 * \param[in, out]  h_cpy   set of parameters to update
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_copy_parameters(const cs_hodge_param_t   *h_ref,
                         cs_hodge_param_t         *h_cpy)
{
  if (h_ref == NULL || h_cpy == NULL)
    return;

  h_cpy->inv_pty = h_ref->inv_pty;
  h_cpy->type = h_ref->type;
  h_cpy->algo = h_ref->algo;
  h_cpy->coef = h_ref->coef;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the property value (scalar- or tensor-valued) related to a
 *         discrete Hodge operator inside a cell and if needed other related
 *         quantities
 *
 * \param[in]      c_id    id of the cell to deal with
 * \param[in]      t_eval  time at which one performs the evaluation
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_set_property_value(const cs_lnum_t       c_id,
                            const cs_real_t       t_eval,
                            const cs_flag_t       c_flag,
                            cs_hodge_t           *hodge)
{
  assert(hodge != NULL);

  cs_property_data_t  *ptyd = hodge->pty_data;

  if (ptyd->property == NULL)
    return;  /* The default initialization corresponds to what is needed */

  if (ptyd->need_tensor) {

    cs_property_get_cell_tensor(c_id,
                                t_eval,
                                ptyd->property,
                                hodge->param->inv_pty,
                                ptyd->tensor);

    if (ptyd->is_iso)
      ptyd->value = ptyd->tensor[0][0];

  }
  else {

    if (ptyd->is_iso) {

      ptyd->value = cs_property_get_cell_value(c_id,
                                               t_eval,
                                               ptyd->property);

      if (hodge->param->inv_pty) {
        assert(fabs(ptyd->value) > FLT_MIN);
        ptyd->value = 1/ptyd->value;
      }

    }
    else { /* Anisotropic so a tensor is needed */

      ptyd->need_tensor = true;
      cs_property_get_cell_tensor(c_id,
                                  t_eval,
                                  ptyd->property,
                                  hodge->param->inv_pty,
                                  ptyd->tensor);

    }

  } /* Tensor-valued evaluation of the property is not needed */

  /* Set additional quantities in case of more advanced way of enforcing the
     essential BCs */

  if (c_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {
    if (ptyd->need_eigen) {

      if (ptyd->need_tensor)
        cs_math_33_eigen(ptyd->tensor,
                         &(ptyd->eigen_ratio), &(ptyd->eigen_max));
      else
        ptyd->eigen_ratio = 1.0, ptyd->eigen_max = ptyd->value;

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the property value (scalar- or tensor-valued) related to a
 *         discrete Hodge operator inside a cell and if needed ohter related
 *         quantities.
 *         Cell-wise variant (usage of cs_cell_mesh_t structure)
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in]      t_eval  time at which one performs the evaluation
 * \param[in]      c_flag  flag related to this cell
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_set_property_value_cw(const cs_cell_mesh_t   *cm,
                               const cs_real_t         t_eval,
                               const cs_flag_t         c_flag,
                               cs_hodge_t             *hodge)
{
  assert(hodge != NULL);

  cs_property_data_t  *ptyd = hodge->pty_data;

  if (ptyd->property == NULL)
    return;  /* The default initialization corresponds to what is needed */

  if (ptyd->need_tensor) {

    cs_property_tensor_in_cell(cm,
                               ptyd->property,
                               t_eval,
                               hodge->param->inv_pty,
                               ptyd->tensor);

    if (ptyd->is_iso)
      ptyd->value = ptyd->tensor[0][0];

  }
  else {

    if (ptyd->is_iso) {

      ptyd->value = cs_property_value_in_cell(cm, ptyd->property, t_eval);

      if (hodge->param->inv_pty) {
        assert(fabs(ptyd->value) > FLT_MIN);
        ptyd->value = 1/ptyd->value;
      }

    }
    else { /* Anisotropic so a tensor is needed */

      ptyd->need_tensor = true;
      cs_property_tensor_in_cell(cm,
                                 ptyd->property,
                                 t_eval,
                                 hodge->param->inv_pty,
                                 ptyd->tensor);

    }

  } /* Tensor-valued evaluation of the property is not needed */

  /* Set additional quantities in case of more advanced way of enforcing the
     essential BCs */

  if (c_flag & CS_FLAG_BOUNDARY_CELL_BY_FACE) {
    if (ptyd->need_eigen) {

      if (ptyd->need_tensor)
        cs_math_33_eigen(ptyd->tensor,
                         &(ptyd->eigen_ratio), &(ptyd->eigen_max));
      else
        ptyd->eigen_ratio = 1.0, ptyd->eigen_max = ptyd->value;

    }
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fb_cost_get_stiffness(const cs_cell_mesh_t     *cm,
                               cs_hodge_t               *hodge,
                               cs_cell_builder_t        *cb)
{
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  /* Initialize the local stiffness matrix */

  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_fc + 1, sloc);

  /* Compute the local discrete Hodge operator */

  bool  continue_computation = cs_hodge_edfp_cost_get_opt(cm, hodge, cb);

  if (!continue_computation)
    return continue_computation;

  const cs_sdm_t  *hmat = hodge->matrix;

  double  *sval_crow = sloc->val + cm->n_fc*sloc->n_rows;
  double  full_sum = 0.;

  for (int i = 0; i < hmat->n_rows; i++) {

    const short int  fi_sgn = cm->f_sgn[i];
    const double  *hval_i = hmat->val + i*hmat->n_rows;

    double  *sval_i = sloc->val + i*sloc->n_rows;
    double  row_sum = 0.;
    for (int j = 0; j < hmat->n_rows; j++) {
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

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, sloc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          with the usage of bubble stabilization.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fb_bubble_get_stiffness(const cs_cell_mesh_t    *cm,
                                 cs_hodge_t              *hodge,
                                 cs_cell_builder_t       *cb)
{
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  /* Initialize the local stiffness matrix */

  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_fc + 1, sloc);

  /* Compute the local discrete Hodge operator */

  bool  continue_computation = cs_hodge_edfp_bubble_get(cm, hodge, cb);

  if (!continue_computation)
    return continue_computation;

  const cs_sdm_t  *hmat = hodge->matrix;

  double  *sval_crow = sloc->val + cm->n_fc*sloc->n_rows;
  double  full_sum = 0.;

  for (int i = 0; i < hmat->n_rows; i++) {

    const short int  fi_sgn = cm->f_sgn[i];
    const double  *hval_i = hmat->val + i*hmat->n_rows;

    double  *sval_i = sloc->val + i*sloc->n_rows;
    double  row_sum = 0.;
    for (int j = 0; j < hmat->n_rows; j++) {
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

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, sloc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Voronoi algorithm
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fb_voro_get_stiffness(const cs_cell_mesh_t     *cm,
                               cs_hodge_t               *hodge,
                               cs_cell_builder_t        *cb)
{
  assert(hodge->param->type == CS_HODGE_TYPE_EDFP);
  assert(hodge->param->algo == CS_HODGE_ALGO_VORONOI);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PF | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  /* Initialize the local stiffness matrix */

  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_fc + 1, sloc);

  /* Compute the local discrete Hodge operator */

  bool  continue_computation = cs_hodge_edfp_voro_get(cm, hodge, cb);

  if (!continue_computation)
    return  continue_computation;

  const cs_sdm_t  *hmat = hodge->matrix;

  double  full_sum = 0.;
  double  *sval_crow = sloc->val + cm->n_fc*sloc->n_rows;

  for (int i = 0; i < hmat->n_rows; i++) {

    /* Hodge operator is diagonal */

    const double  *hval_i = hmat->val + i*hmat->n_rows;
    const double  row_sum = hval_i[i];

    double  *sval_i = sloc->val + i*sloc->n_rows;

    sval_i[i] = hval_i[i];
    sval_i[cm->n_fc] = -row_sum;
    sval_crow[i] = -row_sum;
    full_sum += row_sum;

  } /* Loop on rows (faces) */

  /* (c, c) diagonal entry */

  sval_crow[cm->n_fc] = full_sum;

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, sloc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes and an isotropic property
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_cost_get_iso_stiffness(const cs_cell_mesh_t   *cm,
                                   cs_hodge_t             *hodge,
                                   cs_cell_builder_t      *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(hodgep->type == CS_HODGE_TYPE_EPFD);
  assert(hodgep->algo == CS_HODGE_ALGO_COST);
  assert(ptyd->is_iso);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ |
                       CS_FLAG_COMP_EV));

  if (!(fabs(ptyd->value) > 0))
    return false; /* One avoids to compute the Hodge op. for nothing */

  /* Set numbering and geometrical quantities Hodge builder */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_ec;

  _init_vb_geom_quant(cm, pq, dq);

  /* Compute the upper right part of the local Hodge matrix
   *  Rk: Switch arguments between discrete Hodge operator from PRIMAL->DUAL
   *  or DUAL->PRIMAL space */

  cs_sdm_square_init(cm->n_ec, hodge->matrix);

  _compute_iso_hodge_ur(cm->n_ec,
                        3*hodgep->coef*hodgep->coef,
                        1./cm->vol_c,
                        ptyd->value,
                        (const cs_real_t (*)[3])pq,
                        (const cs_real_t (*)[3])dq,
                        cb, hodge->matrix);

  _define_vb_stiffness(cm, hodge->matrix, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, cb->loc);
#endif

  return true; /* Something has been done */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes and an anistropic property
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge    pointer to a cs_hodge_t structure
 * \param[in, out] cb       pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_cost_get_aniso_stiffness(const cs_cell_mesh_t    *cm,
                                     cs_hodge_t              *hodge,
                                     cs_cell_builder_t       *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(hodgep->type == CS_HODGE_TYPE_EPFD);
  assert(hodgep->algo == CS_HODGE_ALGO_COST);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ |
                       CS_FLAG_COMP_EV));

  /* Set numbering and geometrical quantities Hodge builder */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_ec;

  _init_vb_geom_quant(cm, pq, dq);

  /* Compute the upper right part of the local Hodge matrix
   * Remark: Switch arguments between discrete Hodge operator from PRIMAL->DUAL
   * or DUAL->PRIMAL space */

  cs_sdm_square_init(cm->n_ec, hodge->matrix);

  if (_tensor_norm_l1(ptyd->tensor) > 0)
    _compute_aniso_hodge_ur(cm->n_ec,
                            3*hodgep->coef*hodgep->coef,
                            1./cm->vol_c,
                            ptyd->tensor,
                            (const cs_real_t (*)[3])pq,
                            (const cs_real_t (*)[3])dq,
                            cb, hodge->matrix);
  else
    return false; /* No need to perform a computation */

  _define_vb_stiffness(cm, hodge->matrix, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, cb->loc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic Bubble algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes and isotropic material property
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_bubble_get_iso_stiffness(const cs_cell_mesh_t    *cm,
                                     cs_hodge_t              *hodge,
                                     cs_cell_builder_t       *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(hodgep->type == CS_HODGE_TYPE_EPFD);
  assert(hodgep->algo == CS_HODGE_ALGO_BUBBLE);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ |
                       CS_FLAG_COMP_EV));

  if (!(fabs(ptyd->value) > 0))
    return false; /* One avoids to compute the Hodge op. for nothing */

  /* Set numbering and geometrical quantities Hodge builder */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_ec;

  _init_vb_geom_quant(cm, pq, dq);

  /* Compute the upper right part of the local Hodge matrix
   *  Rk: Switch arguments between discrete Hodge operator from PRIMAL->DUAL
   *  or DUAL->PRIMAL space */

  cs_sdm_square_init(cm->n_ec, hodge->matrix);

  _compute_iso_bubble_hodge_ur(cm->n_ec,
                               hodgep->coef,
                               1./cm->vol_c,
                               ptyd->value,
                               (const cs_real_t (*)[3])pq,
                               (const cs_real_t (*)[3])dq,
                               cb, hodge->matrix);

  _define_vb_stiffness(cm, hodge->matrix, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, cb->loc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic Bubble algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes and anisotropic material property
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_bubble_get_aniso_stiffness(const cs_cell_mesh_t    *cm,
                                       cs_hodge_t              *hodge,
                                       cs_cell_builder_t       *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(hodgep->type == CS_HODGE_TYPE_EPFD);
  assert(hodgep->algo == CS_HODGE_ALGO_BUBBLE);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ |
                       CS_FLAG_COMP_EV));

  if (!(_tensor_norm_l1(ptyd->tensor) > 0))
    return false;

  /* Set numbering and geometrical quantities Hodge builder */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_ec;

  _init_vb_geom_quant(cm, pq, dq);

  /* Compute the upper right part of the local Hodge matrix
   *  Rk: Switch arguments between discrete Hodge operator from PRIMAL->DUAL
   *  or DUAL->PRIMAL space */

  cs_sdm_square_init(cm->n_ec, hodge->matrix);

  _compute_aniso_bubble_hodge_ur(cm->n_ec,
                                 hodgep->coef,
                                 1./cm->vol_c,
                                 ptyd->tensor,
                                 (const cs_real_t (*)[3])pq,
                                 (const cs_real_t (*)[3])dq,
                                 cb, hodge->matrix);

  _define_vb_stiffness(cm, hodge->matrix, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, cb->loc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Orthogonal
 *          Consistent/Sub-Stabilization decomposition (OCS2) with a
 *          subdivision of pvol_{e,c}.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case Vb schemes and an anisotropic material property
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_ocs2_get_aniso_stiffness(const cs_cell_mesh_t     *cm,
                                     cs_hodge_t               *hodge,
                                     cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(hodgep->type == CS_HODGE_TYPE_EPFD);
  assert(hodgep->algo == CS_HODGE_ALGO_OCS2);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ |
                       CS_FLAG_COMP_EV | CS_FLAG_COMP_SEF));

  /* Initialize the hodge matrix */

  cs_sdm_square_init(cm->n_ec, hodge->matrix);

  /* Compute the upper right part of the local Hodge matrix
   *  Rk: Switch arguments between discrete Hodge operator from PRIMAL->DUAL
   *  or DUAL->PRIMAL space */

  _compute_aniso_hepfd_ocs2_ur(3*hodgep->coef*hodgep->coef, ptyd->tensor, cm,
                               cb, hodge->matrix);

  _define_vb_stiffness(cm, hodge->matrix, cb->loc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, cb->loc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic COST algo.
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_cost_get_stiffness(const cs_cell_mesh_t     *cm,
                               cs_hodge_t               *hodge,
                               cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(hodgep->type == CS_HODGE_TYPE_EPFD);
  assert(hodgep->algo == CS_HODGE_ALGO_COST);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ |
                       CS_FLAG_COMP_EV));

  /* Set numbering and geometrical quantities Hodge builder */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_ec;

  _init_vb_geom_quant(cm, pq, dq);

  /* Compute additional geometrical quantities.
     Initial the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  double  *kappa = cb->values;
  double  *alpha = cb->values + cm->n_ec;

  /* Initialize the hodge matrix */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_ec, hmat);

  const double  invcvol = 1/cm->vol_c;
  const double  beta2 = hodgep->coef*hodgep->coef;

  if (ptyd->is_iso || ptyd->is_unity) {

    if (fabs(ptyd->value) > 0)
      _compute_cost_quant_iso(cm->n_ec, invcvol, ptyd->value,
                              (const cs_real_t (*)[3])pq,
                              (const cs_real_t (*)[3])dq,
                              alpha, kappa, hmat);
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0)
      _compute_cost_quant(cm->n_ec, invcvol, ptyd->tensor,
                          (const cs_real_t (*)[3])pq,
                          (const cs_real_t (*)[3])dq,
                          alpha, kappa, hmat);
    else
      return false; /* No need to perform a computation */

  }

  /* Initialize the local stiffness matrix */

  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_vc, sloc);

  for (int ei = 0; ei < cm->n_ec; ei++) { /* Loop on cell edges I */

    const int  shift_i = ei*cm->n_ec;
    const double  *alpha_i = alpha + shift_i;
    const short int  i1ei = cm->e2v_sgn[ei];
    const short int  i1 = cm->e2v_ids[2*ei];
    const short int  i2 = cm->e2v_ids[2*ei+1];
    const double  *hi = hmat->val + shift_i;

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

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, sloc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the Voronoi algorithm
 *          The computed matrix is stored in cb->loc and the related discrete
 *          hodge operator in hodge->matrix
 *          Case of CDO vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_voro_get_stiffness(const cs_cell_mesh_t     *cm,
                               cs_hodge_t               *hodge,
                               cs_cell_builder_t        *cb)
{
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(hodge->param->type == CS_HODGE_TYPE_EPFD);
  assert(hodge->param->algo == CS_HODGE_ALGO_VORONOI);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ |
                       CS_FLAG_COMP_EV));

  /* Initialize the local stiffness matrix */

  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_vc, sloc);

  if (ptyd->is_iso || ptyd->is_unity) {

    double  dpty_val = 1.0;  /* is_unity */
    if (ptyd->is_iso)
      dpty_val = ptyd->value;

    if (!(fabs(dpty_val) > 0))
      return false; /* One avoids to compute the Hodge op. for nothing */

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
      si[vj] = sj[vi] = -dval; /* sgn_i * sgn_j = -1 */

    } /* End of loop on cell edges */

  }
  else { /* Diffusion property is anisotropic */

    if (_tensor_norm_l1(ptyd->tensor) > 0) {

      cs_real_3_t  mv;

      /* Loop on cell edges */

      for (int ii = 0; ii < cm->n_ec; ii++) {

        cs_nvec3_t  dfq = cm->dface[ii];
        cs_quant_t  peq = cm->edge[ii];

        cs_math_33_3_product((const cs_real_3_t *)ptyd->tensor, dfq.unitv, mv);

        /* Only a diagonal term */

        const double  dval = _dp3(mv, dfq.unitv) * dfq.meas/peq.meas;
        const short int  vi = cm->e2v_ids[2*ii];
        const short int  vj = cm->e2v_ids[2*ii+1];

        double  *si = sloc->val + vi*sloc->n_rows;
        double  *sj = sloc->val + vj*sloc->n_rows;

        si[vi] += dval;
        sj[vj] += dval;
        si[vj] = sj[vi] = -dval; /* sgn_j * sgn_i = -1 */

      } /* End of loop on cell edges */

    }
    else
      return false;

  } /* Tensor-valued property */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, sloc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic WBS algo.
 *          WBS standing for Whitney Barycentric Subdivision (WBS)
 *          algo.
 *          The computed matrix is stored in cb->loc
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vb_wbs_get_stiffness(const cs_cell_mesh_t     *cm,
                              cs_hodge_t               *hodge,
                              cs_cell_builder_t        *cb)
{
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(ptyd->need_tensor);
  assert(hodge->param->type == CS_HODGE_TYPE_EPFD);
  assert(hodge->param->algo == CS_HODGE_ALGO_WBS);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PVQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ |
                       CS_FLAG_COMP_EV  | CS_FLAG_COMP_HFQ | CS_FLAG_COMP_FEQ));

  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg;

  cs_real_3_t  *uvc = cb->vectors;
  cs_real_3_t  *glv = cb->vectors + cm->n_vc;
  cs_real_t  *lvc = cb->values;
  cs_real_t  *wvf = cb->values + cm->n_vc;
  cs_real_t  *wef = cb->values + 2*cm->n_vc;

  /* Initialize the local stiffness matrix */

  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(cm->n_vc, sloc);

  if (!(_tensor_norm_l1(ptyd->tensor) > 0))
    return false; /* One avoids to compute the Hodge op. for nothing */

  /* Define the length and unit vector of the segment x_c --> x_v */

  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, lvc + v, uvc[v]);

  /* Loop on cell faces */

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute the gradient of the lagrange function related to a cell
       in each p_{f,c} and the weights for each vertex related to this face */

    cs_compute_grdfc_cw(f, cm, grd_c);
    cs_compute_wef_wvf(f, cm, wvf, wef);

    /* Loop on face edges to scan p_{e,f,c} subvolumes */

    const short int  *f2e_idx = cm->f2e_idx + f;
    const short int  *f2e_ids = cm->f2e_ids + f2e_idx[0];
    for (int i = 0; i < f2e_idx[1] - f2e_idx[0]; i++) {

      const short int  ee = 2*f2e_ids[i];
      const double  subvol = wef[i]*cm->pvol_f[f];
      const short int  v1 = cm->e2v_ids[ee];
      const short int  v2 = cm->e2v_ids[ee+1];

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

        if (wvf[si] > 0) /* Face contrib. */
          for (int k = 0; k < 3; k++)
            glv[si][k] += wvf[si]*grd_f[k];

        if (si == v1) /* Vertex 1 contrib */
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v1[k];

        if (si == v2) /* Vertex 2 contrib */
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v2[k];

      } /* Loop on cell vertices */

      /* Build the upper right part */

      for (int si = 0; si < sloc->n_rows; si++) {

        cs_math_33_3_product(ptyd->tensor, glv[si], matg);

        /* Diagonal contribution */

        double  *mi = sloc->val + si*sloc->n_rows;
        mi[si] += subvol * _dp3(matg, glv[si]);

        /* Loop on vertices v_j (j > i) */
        for (int sj = si+1; sj < sloc->n_rows; sj++)
          mi[sj] += subvol * _dp3(matg, glv[sj]);

      } /* Loop on vertices v_i */

    }

  } /* Loop on cell faces */

  /* Matrix is symmetric by construction */

  cs_sdm_symm_ur(sloc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, sloc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local stiffness matrix using the generic WBS algo.
 *          WBS standing for Whitney Barycentric Subdivision (WBS) algo.
 *          The computed matrix is stored in cb->loc
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vcb_get_stiffness(const cs_cell_mesh_t     *cm,
                           cs_hodge_t               *hodge,
                           cs_cell_builder_t        *cb)
{
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(ptyd->need_tensor);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ |
                       CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_PFC));

  cs_real_3_t  grd_c, grd_f, grd_v1, grd_v2, matg, matg_c;

  cs_real_3_t  *uvc = cb->vectors;
  cs_real_3_t  *glv = cb->vectors + cm->n_vc;
  cs_real_t  *lvc = cb->values;
  cs_real_t  *wvf = cb->values + cm->n_vc;
  cs_real_t  *wef = cb->values + 2*cm->n_vc;

  /* Initialize the local stiffness matrix */

  const int  nc_dofs = cm->n_vc + 1;
  /* index to the (cell,cell) entry */
  const int  cc = nc_dofs*cm->n_vc + cm->n_vc;

  cs_sdm_t  *sloc = cb->loc;
  cs_sdm_square_init(nc_dofs, sloc);

  if (!(_tensor_norm_l1(ptyd->tensor) > 0))
    return false; /* One avoids to compute the Hodge op. for nothing */

  /* Define the length and unit vector of the segment x_c --> x_v */

  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, lvc + v, uvc[v]);

  /* Loop on cell faces */

  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute for the current face:
       - the gradient of the Lagrange function related xc in p_{f,c}
       - weights related to vertices
       - subvolume p_{ef,c} related to edges */

    cs_compute_grdfc_cw(f, cm, grd_c);
    cs_compute_wef_wvf(f, cm, wvf, wef);

    /* Compute the contribution to the entry A(c,c) */

    cs_math_33_3_product(ptyd->tensor, grd_c, matg_c);
    sloc->val[cc] += cm->pvol_f[f] * _dp3(grd_c, matg_c);

    /* Loop on face edges to scan p_{e,f,c} subvolumes */

    const short int  *f2e_idx = cm->f2e_idx + f;
    const short int  *f2e_ids = cm->f2e_ids + f2e_idx[0];
    for (short int i = 0; i < f2e_idx[1] - f2e_idx[0]; i++) {

      const short int  ee = 2*f2e_ids[i];
      const double  subvol = wef[i]*cm->pvol_f[f];
      const short int  v1 = cm->e2v_ids[ee];
      const short int  v2 = cm->e2v_ids[ee+1];

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

        if (wvf[si] > 0) /* Face contrib. */
          for (int k = 0; k < 3; k++)
            glv[si][k] += wvf[si]*grd_f[k];

        if (si == v1) /* Vertex 1 contrib */
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v1[k];

        if (si == v2) /* Vertex 2 contrib */
          for (int k = 0; k < 3; k++)
            glv[si][k] += grd_v2[k];

      } /* Loop on cell vertices */

      /* Build the upper right part (v-v and v-c)
         Be careful: sloc->n_rows = cm->n_vc + 1 */

      for (int si = 0; si < cm->n_vc; si++) {

        double  *mi = sloc->val + si*sloc->n_rows;

        /* Add v-c contribution */

        mi[cm->n_vc] += subvol * _dp3(matg_c, glv[si]);

        /* Add v-v contribution */

        cs_math_33_3_product(ptyd->tensor, glv[si], matg);

        /* Loop on vertices v_j (j >= i) */

        for (int sj = si; sj < cm->n_vc; sj++)
          mi[sj] += subvol * _dp3(matg, glv[sj]);

      } /* Loop on vertices v_i */

    }

  } /* Loop on cell faces */

  /* Matrix is symmetric by construction */

  cs_sdm_symm_ur(sloc);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_stiffness(cm->c_id, sloc);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a local Hodge operator on a given cell which is equivalent of
 *         a mass matrix. It relies on a CO+ST algo. and is specific to CDO-Fb
 *         schemes.
 *         The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fb_get(const cs_cell_mesh_t     *cm,
                cs_hodge_t               *hodge,
                cs_cell_builder_t        *cb)
{
  CS_UNUSED(cb);

  assert(hodge != NULL);
  assert(hodge->param->type == CS_HODGE_TYPE_FB);
  assert(hodge->param->algo == CS_HODGE_ALGO_COST);
  assert(hodge->pty_data->is_unity);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ));

  const int n_cols = cm->n_fc + 1;
  const cs_real_t  over_cell = 1./(cm->vol_c*cm->vol_c);

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(n_cols, hmat);

  /* cell-cell entry (cell-face and face-cell block are null) */

  hmat->val[cm->n_fc*(n_cols + 1)] = cm->vol_c;

  /* Compute the inertia (useful for the face-face block) */

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

    hmat->val[fi*(n_cols+1)] = over_cell * pfq_i.meas*pfq_i.meas * dval;

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

      hmat->val[fi*n_cols+fj] = eval;
      hmat->val[fj*n_cols+fi] = eval;

    } /* Loop on cell faces (fj) */

  } /* Loop on cell faces (fi) */

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Voronoi
 *          algo. This leads to a diagonal operator.
 *          This function is specific for vertex+cell-based schemes
 *          The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vcb_voro_get(const cs_cell_mesh_t     *cm,
                      cs_hodge_t               *hodge,
                      cs_cell_builder_t        *cb)
{
  CS_UNUSED(cb);

  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cm != NULL && hodge != NULL);
  assert(hodge->param->type == CS_HODGE_TYPE_VC);
  assert(hodge->param->algo == CS_HODGE_ALGO_VORONOI);
  assert(ptyd->is_iso);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  cs_sdm_t  *hmat = hodge->matrix;

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_square_init(cm->n_vc + 1, hmat);

  const int  msize = cm->n_vc + 1;

  if (ptyd->is_unity) {

    /* H(c,c) = 0.25*|c| */

    hmat->val[msize*cm->n_vc] = 0.25*cm->vol_c;

    /* H(c,c) = 0.75*|dcell(v) \cap c| */

    const double  vol_coef = 0.75*cm->vol_c;
    for (short int vi = 0; vi < cm->n_vc; vi++)
      hmat->val[msize*vi] = vol_coef*cm->wvc[vi];

  }
  else {

    if (fabs(ptyd->value) > 0) {

      /* H(c,c) = 0.25*|c| */

      hmat->val[msize*cm->n_vc] = ptyd->value*0.25*cm->vol_c;

      /* H(c,c) = 0.75*|dcell(v) \cap c| */

      const double  vol_coef = 0.75*cm->vol_c*ptyd->value;
      for (short int vi = 0; vi < cm->n_vc; vi++)
        hmat->val[msize*vi] = vol_coef*cm->wvc[vi];

    }
    else
      return false; /* No need to perform a computation */

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hmat);
  }
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the WBS algo.
 *          This function is specific for vertex+cell-based schemes
 *          The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vcb_wbs_get(const cs_cell_mesh_t     *cm,
                     cs_hodge_t               *hodge,
                     cs_cell_builder_t        *cb)
{
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL && hodge != NULL);
  assert(hodge->param->type == CS_HODGE_TYPE_VC);
  assert(hodge->param->algo == CS_HODGE_ALGO_WBS);
  assert(ptyd->is_iso);
  assert(cs_eflag_test(cm->flag,
                      CS_FLAG_COMP_PVQ | CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ |
                      CS_FLAG_COMP_EV  | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ));

  cs_sdm_t  *hmat = hodge->matrix;

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_square_init(cm->n_vc + 1, hmat);

  if (!(fabs(ptyd->value) > 0))
    return false; /* One avoids to compute the Hodge op. for nothing */

  double  *wvf = cb->values;
  double  *wef = cb->values + cm->n_vc;

  const int  msize = cm->n_vc + 1;
  const double  c_coef1 = 0.2*cm->vol_c;
  const double  c_coef2 = cs_hodge_vc_coef * cm->vol_c;

  /* H(c,c) = 0.1*|c| */

  hmat->val[msize*cm->n_vc + cm->n_vc] = 0.1*cm->vol_c;

  /* Initialize the upper part of the local Hodge matrix
     diagonal and cell column entries */

  for (short int vi = 0; vi < cm->n_vc; vi++) {

    double  *mi = hmat->val + vi*msize;

    mi[vi] = c_coef1 * cm->wvc[vi];       /* Diagonal entry */
    for (short int vj = vi+1; vj < cm->n_vc; vj++)
      mi[vj] = 0.;
    mi[cm->n_vc] = c_coef2 * cm->wvc[vi]; /* Cell column */

  } /* Loop on cell vertices */

  /* Loop on each pef and add the contribution */

  for (short int f = 0; f < cm->n_fc; f++) {

    /* Define useful quantities for WBS algo. */

    cs_compute_wef_wvf(f, cm, wvf, wef);

    const double f_coef = 0.3 * cm->pvol_f[f];
    const double ef_coef = 0.05 * cm->pvol_f[f];

    /* Add face contribution:
       Diagonal entry    H(i,i) += 0.3*wif*wif*pfc_vol
       Extra-diag. entry H(i,j) += 0.3*wjf*wif*pfc_vol */

    for (short int vi = 0; vi < cm->n_vc; vi++) {

      const double  coef_if = f_coef * wvf[vi];
      double  *mi = hmat->val + vi*msize;

      /* Diagonal and Extra-diagonal entries: Add face contribution */

      for (short int vj = vi; vj < cm->n_vc; vj++)
        mi[vj] += coef_if * wvf[vj];

    } /* Extra-diag entries */

    /* Add edge-face contribution (only extra-diag) = 0.05 * |p_{ef,c}| */

    const short int  *f2e_idx = cm->f2e_idx + f;
    const short int  *f2e_ids = cm->f2e_ids + f2e_idx[0];
    for (short int i = 0; i < f2e_idx[1] - f2e_idx[0]; i++) {

      const short int  ee = 2*f2e_ids[i];
      const short int  v1 = cm->e2v_ids[ee];
      const short int  v2 = cm->e2v_ids[ee+1];

      assert(v1 > -1 && v2 > -1);
      if (v1 < v2)
        hmat->val[v1*msize+v2] += ef_coef * wef[i];
      else
        hmat->val[v2*msize+v1] += ef_coef * wef[i];

    } /* Loop on face edges */

  } /* Loop on cell faces */

  /* Take into account the value of the associated property */

  if (!ptyd->is_unity) {
    for (short int vi = 0; vi < msize; vi++) {
      double  *mi = hmat->val + vi*msize;
      for (short int vj = vi; vj < msize; vj++)
        mi[vj] *= ptyd->value;
    }
  }

  /* Local matrix is symmetric by construction. Set the lower part. */

  for (short int vj = 0; vj < msize; vj++) {
    double  *mj = hmat->val + vj*msize;
    for (short int vi = vj+1; vi < msize; vi++)
      hmat->val[vi*msize + vj] = mj[vi];
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hmat);
  }
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a local Hodge operator for a given cell using WBS algo.
 *         Hodge op. from primal vertices to dual cells.
 *         This function is specific for vertex-based schemes
 *         The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vpcd_wbs_get(const cs_cell_mesh_t    *cm,
                      cs_hodge_t              *hodge,
                      cs_cell_builder_t       *cb)
{
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL);
  assert(hodge->param->type == CS_HODGE_TYPE_VPCD);
  assert(hodge->param->algo == CS_HODGE_ALGO_WBS);
  assert(ptyd->is_iso);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PVQ |CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ |
                       CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_PFC));

  double  *wvf = cb->values;
  double  *wef = cb->values + cm->n_vc;

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_vc, hmat);

  if (!(fabs(ptyd->value) > 0))
    return false; /* One avoids to compute the Hodge op. for nothing */

  const double  c_coef = 0.1*cm->vol_c;

  /* Initialize the upper part of the local Hodge matrix */

  for (short int vi = 0; vi < cm->n_vc; vi++) {

    double  *mi = hmat->val + vi*cm->n_vc;
    const double  vi_coef = 4 * c_coef * cm->wvc[vi];

    /* Diag. entry has an additional contrib */
    mi[vi] = vi_coef * (0.5 + cm->wvc[vi]);
    for (short int vj = vi + 1; vj < cm->n_vc; vj++)
      mi[vj] = vi_coef * cm->wvc[vj]; /* Extra-diagonal entries */

  } /* Loop on cell vertices */

  /* Loop on each pef and add the contribution */

  for (short int f = 0; f < cm->n_fc; f++) {

    /* Define useful quantities for WBS algo. */

    cs_compute_wef_wvf(f, cm, wvf, wef);

    const double  f_coef = 0.3 * cm->pvol_f[f];
    const double  ef_coef = 0.05 * cm->pvol_f[f];

    /* Add face contribution:
       Diagonal entry    H(i,i) += 0.3*wif*wif*pfc_vol
       Extra-diag. entry H(i,j) += 0.3*wjf*wif*pfc_vol */

    for (short int vi = 0; vi < cm->n_vc; vi++) {

      double  *mi = hmat->val + vi*cm->n_vc;

      /* Diagonal and Extra-diagonal entries: Add face contribution */
      const double  coef_if = f_coef * wvf[vi];
      for (short int vj = vi; vj < cm->n_vc; vj++)
        mi[vj] += coef_if * wvf[vj];

    } /* Face contribution */

    /* Add edge-face contribution (only extra-diag) = 0.05 * |p_{ef,c}| */

    const short int  *f2e_idx = cm->f2e_idx + f;
    const short int  *f2e_ids = cm->f2e_ids + f2e_idx[0];
    for (short int i = 0; i < f2e_idx[1] - f2e_idx[0]; i++) {

      const short int  ee = 2*f2e_ids[i];
      const short int  v1 = cm->e2v_ids[ee];
      const short int  v2 = cm->e2v_ids[ee+1];

      assert(v1 > -1 && v2 > -1);
      if (v1 < v2)
        hmat->val[v1*cm->n_vc+v2] += ef_coef * wef[i];
      else
        hmat->val[v2*cm->n_vc+v1] += ef_coef * wef[i];

    } /* Loop on face edges */

  } /* Loop on cell faces */

  /* Take into account the value of the associated property */

  if (!ptyd->is_unity) {
    for (short int vi = 0; vi < cm->n_vc; vi++) {
      double  *mi = hmat->val + vi*cm->n_vc;
      for (short int vj = vi; vj < cm->n_vc; vj++)
        mi[vj] *= ptyd->value;
    }
  }

  /* Local matrix is symmetric by construction. Set the lower part. */

  for (short int vj = 0; vj < cm->n_vc; vj++) {
    double  *mj = hmat->val + vj*cm->n_vc;
    for (short int vi = vj+1; vi < cm->n_vc; vi++)
      hmat->val[vi*cm->n_vc + vj] = mj[vi];
  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hmat);
  }
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a local Hodge operator for a given cell using VORONOI algo.
 *         Hodge op. from primal vertices to dual cells.
 *         This function is specific for vertex-based schemes
 *         The discrete Hodge operator is stored in hodge->matrix
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_vpcd_voro_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb)
{
  CS_UNUSED(cb);

  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(hodge->param->type == CS_HODGE_TYPE_VPCD);
  assert(hodge->param->algo == CS_HODGE_ALGO_VORONOI);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PVQ));

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_vc, hmat);

  const int stride = cm->n_vc + 1;
  if (ptyd->is_unity) {

    for (int v = 0; v < cm->n_vc; v++)
      hmat->val[v*stride] = cm->wvc[v] * cm->vol_c;

  }
  else {

    if (fabs(ptyd->value) > 0) {

      const cs_real_t  coef = ptyd->value * cm->vol_c;
      for (int v = 0; v < cm->n_vc; v++)
        hmat->val[v*stride] = coef * cm->wvc[v];

    }
    else
      return false; /* No need to perform a computation */

  }

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hmat);
  }
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal edges to dual faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_epfd_voro_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb)
{
  CS_UNUSED(cb);

  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(hodge->param->type == CS_HODGE_TYPE_EPFD);
  assert(hodge->param->algo == CS_HODGE_ALGO_VORONOI);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ | CS_FLAG_COMP_SEF));

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_ec, hmat);

  /* Set the diagonal entries */

  if (ptyd->is_iso) {

    if (fabs(ptyd->value) > 0) {

      for (short int e = 0; e < cm->n_ec; e++)
        hmat->val[e*cm->n_ec+e] =
          ptyd->value * cm->dface[e].meas/cm->edge[e].meas;

    }
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0) {

      cs_real_3_t  mv;
      for (short int f = 0; f < cm->n_fc; f++) {
        for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {
          const short int  e = cm->f2e_ids[i];
          const cs_nvec3_t  *sefc = cm->sefc + i;
          cs_math_33_3_product(ptyd->tensor, sefc->unitv, mv);
          hmat->val[e*cm->n_ec+e] += sefc->meas * _dp3(mv, sefc->unitv);
        }
      }

      for (short int e = 0; e < cm->n_ec; e++)
        hmat->val[e*cm->n_ec+e] /= cm->edge[e].meas;

    }
    else
      return false; /* No need to perform a computation */

  } /* Anisotropic */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hmat);
  }
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from primal edges to dual faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_epfd_cost_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL && hodgep != NULL);
  assert(hodgep->type == CS_HODGE_TYPE_EPFD);
  assert(hodgep->algo == CS_HODGE_ALGO_COST);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ));

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_ec, hmat);

  /* Set numbering and geometrical quantities Hodge builder */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_ec;

  _init_vb_geom_quant(cm, pq, dq);

  /* Compute additional geometrical quantities: qmq and T
     Initial the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  double  *kappa = cb->values;
  double  *alpha = cb->values + cm->n_ec;

  if (ptyd->is_unity)
    _compute_cost_quant_iso(cm->n_ec, 1/cm->vol_c, 1.0,
                            (const cs_real_t (*)[3])pq,
                            (const cs_real_t (*)[3])dq,
                            alpha, kappa, hmat);

  else if (ptyd->is_iso) {

    if (fabs(ptyd->value) > 0)
      _compute_cost_quant_iso(cm->n_ec, 1/cm->vol_c, ptyd->value,
                              (const cs_real_t (*)[3])pq,
                              (const cs_real_t (*)[3])dq,
                              alpha, kappa, hmat);
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0)
      _compute_cost_quant(cm->n_ec, 1/cm->vol_c,
                          (const cs_real_3_t *)ptyd->tensor,
                          (const cs_real_t (*)[3])pq,
                          (const cs_real_t (*)[3])dq,
                          alpha, kappa, hmat);
    else
      return false; /* No need to perform a computation */

  }

  double  beta2 = hodgep->coef*hodgep->coef;
  _compute_hodge_cost(cm->n_ec, beta2, alpha, kappa, hmat->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_vector_hodge(cm->c_id,
                      (const cs_real_t (*)[3])pq, (const cs_real_t (*)[3])dq,
                      hodge, cb);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          with a bubble stabilization.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          Hodge op. from primal edges to dual faces. This function is
 *          specific for vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_epfd_bubble_get(const cs_cell_mesh_t     *cm,
                         cs_hodge_t               *hodge,
                         cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL && hodgep != NULL && ptyd != NULL);
  assert(hodgep->type == CS_HODGE_TYPE_EPFD);
  assert(hodgep->algo == CS_HODGE_ALGO_BUBBLE);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ));

  /* Set numbering and geometrical quantities Hodge builder */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_ec;

  _init_vb_geom_quant(cm, pq, dq);

  /* Compute the upper right part of the local Hodge matrix
   *  Rk: Switch arguments between discrete Hodge operator from PRIMAL->DUAL
   *  or DUAL->PRIMAL space */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_ec, hmat);

  if (ptyd->is_iso || ptyd->is_unity) {

    if (fabs(ptyd->value) > 0)
      _compute_iso_bubble_hodge_ur(cm->n_ec,
                                   hodgep->coef,
                                   1./cm->vol_c,
                                   ptyd->value,
                                   (const cs_real_t (*)[3])pq,
                                   (const cs_real_t (*)[3])dq,
                                   cb, hmat);
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0)
      _compute_aniso_bubble_hodge_ur(cm->n_ec,
                                     hodgep->coef,
                                     1./cm->vol_c,
                                     ptyd->tensor,
                                     (const cs_real_t (*)[3])pq,
                                     (const cs_real_t (*)[3])dq,
                                     cb, hmat);
    else
      return false; /* No need to perform a computation */

  }

  /* Hodge matrix is symmetric by construction */

  cs_sdm_symm_ur(hmat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_vector_hodge(cm->c_id,
                      (const cs_real_t (*)[3])pq, (const cs_real_t (*)[3])dq,
                      hodge, cb);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Orthogonal
 *          Consistent/Sub-Stabilization decomposition (OCS2) with a
 *          subdivision of pvol_{e,c}.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          Hodge op. from primal edges to dual faces.
 *          This function is specific for vertex-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_epfd_ocs2_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL && hodgep != NULL);
  assert(hodgep->type == CS_HODGE_TYPE_EPFD);
  assert(hodgep->algo == CS_HODGE_ALGO_OCS2);
  assert(cs_eflag_test(cm->flag,
                       CS_FLAG_COMP_PV | CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ |
                       CS_FLAG_COMP_EV | CS_FLAG_COMP_SEF));

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_ec, hmat);

  /* Compute the upper right part of the local Hodge matrix */

  if (_tensor_norm_l1(ptyd->tensor) > 0)
    _compute_aniso_hepfd_ocs2_ur(3*hodgep->coef*hodgep->coef, ptyd->tensor,
                                 cm, cb, hmat);
  else
    return false; /* No need to perform a computation */

  /* Hodge operator leads to a symmetric matrix by construction */

  cs_sdm_symm_ur(hmat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hmat);
  }
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from primal faces to dual edges.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to cell-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fped_voro_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb)
{
  CS_UNUSED(cb);

  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(ptyd != NULL);
  assert(hodge->param->type == CS_HODGE_TYPE_FPED);
  assert(hodge->param->algo == CS_HODGE_ALGO_VORONOI);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ));

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_fc, hmat);

  if (ptyd->is_iso) {

    if (fabs(ptyd->value) > 0) {

      for (short int f = 0; f < cm->n_fc; f++)
        hmat->val[f*cm->n_fc+f] =
          ptyd->value*cm->face[f].meas/cm->dedge[f].meas;

    }
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0) {

      cs_real_3_t  mv;
      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_nvec3_t  deq = cm->dedge[f];
        cs_math_33_3_product(ptyd->tensor, deq.unitv, mv);
        hmat->val[f*cm->n_fc+f] = deq.meas*_dp3(mv, deq.unitv)/cm->face[f].meas;

      } /* Loop on cell faces */

    }
    else
      return false; /* No need to perform a computation */

  } /* Isotropic or anisotropic */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hodge->matrix);
  }
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from primal faces to dual edges.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to cell-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fped_cost_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL && hodgep != NULL && ptyd != NULL);
  assert(hodgep->type == CS_HODGE_TYPE_FPED);
  assert(hodgep->algo == CS_HODGE_ALGO_COST);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_fc;

  /* Initialize the geometrical quantities related to this Hodge operator */

  _init_fb_geom_quant(cm, pq, dq);

  /* Compute additional geometrical quantities:
     Initialize the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  double  *kappa = cb->values;
  double  *alpha = cb->values + cm->n_fc;

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_fc, hmat);

  if (ptyd->is_unity)
    _compute_cost_quant_iso(cm->n_fc, 1/cm->vol_c, 1.0,
                            (const cs_real_t (*)[3])pq,
                            (const cs_real_t (*)[3])dq,
                            alpha, kappa, hmat);

  else if (ptyd->is_iso) {

    if (fabs(ptyd->value) > 0)
      _compute_cost_quant_iso(cm->n_fc, 1/cm->vol_c, ptyd->value,
                              (const cs_real_t (*)[3])pq,
                              (const cs_real_t (*)[3])dq,
                              alpha, kappa, hmat);
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0)
      _compute_cost_quant(cm->n_fc, 1/cm->vol_c,
                          (const cs_real_3_t *)ptyd->tensor,
                          (const cs_real_t (*)[3])pq,
                          (const cs_real_t (*)[3])dq,
                          alpha, kappa, hmat);
    else
      return false; /* No need to perform a computation */

  }

  _compute_hodge_cost(cm->n_fc, hodgep->coef*hodgep->coef, alpha, kappa,
                      hmat->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_vector_hodge(cm->c_id,
                      (const cs_real_t (*)[3])pq, (const cs_real_t (*)[3])dq,
                      hodge, cb);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Bubble algo.
 *          Hodge op. from primal faces to dual edges.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to cell-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_fped_bubble_get(const cs_cell_mesh_t     *cm,
                         cs_hodge_t               *hodge,
                         cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL && hodgep != NULL && ptyd != NULL);
  assert(hodgep->type == CS_HODGE_TYPE_EDFP);
  assert(hodgep->algo == CS_HODGE_ALGO_BUBBLE);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  /* Initialize the geometrical quantities related to this Hodge operator */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_fc;

  _init_fb_geom_quant(cm, pq, dq);

  /* Compute additional geometrical quantities:
     Initialize the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_fc, hmat);

  if (ptyd->is_iso || ptyd->is_unity) {

    if (fabs(ptyd->value) > 0)
      _compute_iso_bubble_hodge_ur(cm->n_fc,
                                   hodgep->coef,
                                   1./cm->vol_c,
                                   ptyd->value,
                                   (const cs_real_t (*)[3])dq,
                                   (const cs_real_t (*)[3])pq,
                                   cb, hmat);
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0)
      _compute_aniso_bubble_hodge_ur(cm->n_fc,
                                     hodgep->coef,
                                     1./cm->vol_c,
                                     ptyd->tensor,
                                     (const cs_real_t (*)[3])dq,
                                     (const cs_real_t (*)[3])pq,
                                     cb, hmat);
    else
      return false; /* No need to perform a computation */

  }

  /* Hodge operator is symmetric */

  cs_sdm_symm_ur(hmat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_vector_hodge(cm->c_id,
                      (const cs_real_t (*)[3])dq, (const cs_real_t (*)[3])pq,
                      hodge, cb);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using VORONOI algo.
 *          Hodge op. from dual edges to primal faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_edfp_voro_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb)
{
  CS_UNUSED(cb);

  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(ptyd != NULL);
  assert(hodge->param->type == CS_HODGE_TYPE_EDFP);
  assert(hodge->param->algo == CS_HODGE_ALGO_VORONOI);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ));

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_fc, hmat);

  if (ptyd->is_iso) {

    if (fabs(ptyd->value) > 0) {

      for (short int f = 0; f < cm->n_fc; f++)
        hmat->val[f*cm->n_fc+f] =
          ptyd->value*cm->face[f].meas/cm->dedge[f].meas;

    }
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0) {

      cs_real_3_t  mv;
      for (short int f = 0; f < cm->n_fc; f++) {

        const cs_quant_t  pfq = cm->face[f];
        cs_math_33_3_product(ptyd->tensor, pfq.unitv, mv);
        hmat->val[f*cm->n_fc+f] = pfq.meas * _dp3(mv, pfq.unitv)/cm->edge[f].meas;

      } /* Loop on cell faces */

    }
    else
      return false; /* No need to perform a computation */

  } /* Isotropic or anistropic */

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (cm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Hodge op.   ");
    cs_sdm_dump(cm->c_id, NULL, NULL, hdg);
  }
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from dual edges to primal faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_edfp_cost_get(const cs_cell_mesh_t     *cm,
                       cs_hodge_t               *hodge,
                       cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL && hodgep != NULL && ptyd != NULL);
  assert(hodgep->type == CS_HODGE_TYPE_EDFP);
  assert(hodgep->algo == CS_HODGE_ALGO_COST);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  /* Initialize the geometrical quantities related to this Hodge operator */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_fc;

  _init_fb_geom_quant(cm, pq, dq);

  /* Compute additional geometrical quantities:
     Initialize the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  double  *kappa = cb->values;
  double  *alpha = cb->values + cm->n_fc;

  /* Initialize the local matrix related to this discrete Hodge operator */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_fc, hmat);

  if (ptyd->is_unity)
    _compute_cost_quant_iso(cm->n_fc, 1/cm->vol_c, 1.0,
                            (const cs_real_t (*)[3])dq,
                            (const cs_real_t (*)[3])pq,
                            alpha, kappa, hmat);
  else if (ptyd->is_iso) {

    if (fabs(ptyd->value) > 0)
      _compute_cost_quant_iso(cm->n_fc, 1/cm->vol_c, ptyd->value,
                              (const cs_real_t (*)[3])dq,
                              (const cs_real_t (*)[3])pq,
                              alpha, kappa, hmat);
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0)
      _compute_cost_quant(cm->n_fc, 1/cm->vol_c,
                          (const cs_real_3_t *)ptyd->tensor,
                          (const cs_real_t (*)[3])dq,
                          (const cs_real_t (*)[3])pq,
                          alpha, kappa, hmat);
    else
      return false; /* No need to perform a computation */

  }

  _compute_hodge_cost(cm->n_fc, hodgep->coef*hodgep->coef, alpha, kappa,
                      hmat->val);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_vector_hodge(cm->c_id,
                      (const cs_real_t (*)[3])dq, (const cs_real_t (*)[3])pq,
                      hodge, cb);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the COST algo.
 *          Hodge op. from dual edges to primal faces.
 *          This function is related to face-based schemes
 *
 * \param[in]      cm        pointer to a cs_cell_mesh_t struct.
 * \param[in, out] hodge     pointer to a cs_hodge_t structure
 * \param[in, out] cb        pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_edfp_cost_get_opt(const cs_cell_mesh_t     *cm,
                           cs_hodge_t               *hodge,
                           cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL && hodge != NULL && ptyd != NULL);
  assert(hodgep->type == CS_HODGE_TYPE_EDFP);
  assert(hodgep->algo == CS_HODGE_ALGO_COST);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  /* Initialize the geometrical quantities related to this Hodge operator */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_fc;

  _init_fb_geom_quant(cm, pq, dq);

  /* Compute additional geometrical quantities:
     Initialize the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_fc, hmat);

  if (ptyd->is_iso || ptyd->is_unity) {

    if (fabs(ptyd->value) > 0)
      _compute_iso_hodge_ur(cm->n_fc,
                            3*hodgep->coef*hodgep->coef,
                            1./cm->vol_c,
                            ptyd->value,
                            (const cs_real_t (*)[3])dq,
                            (const cs_real_t (*)[3])pq,
                            cb, hmat);
    else
      return false; /* No need to perform a computation */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0)
      _compute_aniso_hodge_ur(cm->n_fc,
                              3*hodgep->coef*hodgep->coef,
                              1./cm->vol_c,
                              ptyd->tensor,
                              (const cs_real_t (*)[3])dq,
                              (const cs_real_t (*)[3])pq,
                              cb, hmat);
    else
      return false; /* No need to perform a computation */

  }

  /* Hodge operator is symmetric */

  cs_sdm_symm_ur(hmat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_vector_hodge(cm->c_id,
                      (const cs_real_t (*)[3])dq, (const cs_real_t (*)[3])pq,
                      hodge, cb);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build a local Hodge operator for a given cell using the Bubble algo.
 *          Hodge op. from dual edges to primal faces.
 *          The discrete Hodge operator is stored in hodge->matrix
 *          This function is related to face-based schemes
 *
 * \param[in]      cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out] hodge   pointer to a cs_hodge_t structure
 * \param[in, out] cb      pointer to a cs_cell_builder_t structure
 *
 * \return true if something has been computed or false otherwise
 */
/*----------------------------------------------------------------------------*/

bool
cs_hodge_edfp_bubble_get(const cs_cell_mesh_t     *cm,
                         cs_hodge_t               *hodge,
                         cs_cell_builder_t        *cb)
{
  const cs_hodge_param_t  *hodgep = hodge->param;
  const cs_property_data_t  *ptyd = hodge->pty_data;

  assert(cb != NULL && hodgep != NULL && ptyd != NULL);
  assert(hodgep->type == CS_HODGE_TYPE_EDFP);
  assert(hodgep->algo == CS_HODGE_ALGO_BUBBLE);
  assert(cs_eflag_test(cm->flag, CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ));

  /* Initialize the geometrical quantities related to this Hodge operator */

  cs_real_3_t  *pq = cb->vectors;
  cs_real_3_t  *dq = cb->vectors + cm->n_fc;

  _init_fb_geom_quant(cm, pq, dq);

  /* Compute additional geometrical quantities:
     Initialize the local Hodge matrix with the consistency part which is
     constant over a cell.
     Switch arguments between discrete Hodge operator from PRIMAL->DUAL space
     and discrete Hodge operator from DUAL->PRIMAL space */

  cs_sdm_t  *hmat = hodge->matrix;
  cs_sdm_square_init(cm->n_fc, hmat);

  if (ptyd->is_iso || ptyd->is_unity) {

    if (fabs(ptyd->value) > 0)
      _compute_iso_bubble_hodge_ur(cm->n_fc,
                                   hodgep->coef,
                                   1./cm->vol_c,
                                   ptyd->value,
                                   (const cs_real_t (*)[3])dq,
                                   (const cs_real_t (*)[3])pq,
                                   cb, hmat);
    else
      return false; /* Nothing else to compute */

  }
  else {

    if (_tensor_norm_l1(ptyd->tensor) > 0)
      _compute_aniso_bubble_hodge_ur(cm->n_fc,
                                     hodgep->coef,
                                     1./cm->vol_c,
                                     ptyd->tensor,
                                     (const cs_real_t (*)[3])dq,
                                     (const cs_real_t (*)[3])pq,
                                     cb, hmat);
    else
      return false;

  }

  /* Hodge operator is symmetric */

  cs_sdm_symm_ur(hmat);

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 0
  _check_vector_hodge(cm->c_id,
                      (const cs_real_t (*)[3])dq, (const cs_real_t (*)[3])pq,
                      hodge, cb);
#endif

  return true;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute cellwise a discrete hodge operator and multiple it with
 *          a vector
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      hodgep    cs_hodge_param_t structure
 * \param[in]      pty       pointer to a cs_property_t structure or NULL
 * \param[in]      in_vals   vector to multiply with the discrete Hodge op.
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in, out] result    array storing the resulting matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_matvec(const cs_cdo_connect_t       *connect,
                const cs_cdo_quantities_t    *quant,
                const cs_hodge_param_t        hodgep,
                const cs_property_t          *pty,
                const cs_real_t               in_vals[],
                cs_real_t                     t_eval,
                cs_real_t                     result[])
{
  if (in_vals == NULL)
    return;

  const char *func_name = __func__;

  if (result == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              "%s: Resulting vector must be allocated", __func__);
    return; /* Avoid a warning */
  }
  assert(connect != NULL && quant != NULL); /* Sanity checks */

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)             \
  shared(quant, connect, in_vals, t_eval, result, pty, func_name)  \
  firstprivate(hodgep)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    cs_eflag_t  msh_flag = 0;
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = NULL;
    double  *_in = NULL;

    /* Each thread has its own pointer to a cs_hodge_t structure */

    cs_hodge_t  *hodge = cs_hodge_create(connect, pty, &hodgep, true, false);
    cs_hodge_compute_t  *compute = cs_hodge_get_func(func_name, hodgep);
    bool  pty_uniform = cs_property_is_uniform(pty);

    switch (hodgep.type) {

    case CS_HODGE_TYPE_VPCD:

      msh_flag |= CS_FLAG_COMP_PVQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOVB, connect);
      BFT_MALLOC(_in, connect->n_max_vbyc, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_vertices; i++) result[i] = 0;

      switch (hodgep.algo) {

      case CS_HODGE_ALGO_WBS:
        msh_flag |= CS_FLAG_COMP_PVQ |CS_FLAG_COMP_DEQ | CS_FLAG_COMP_PFQ |
          CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ | CS_FLAG_COMP_HFQ;
        break;

      default: /* Minimal requirement for other Hodge algorithms */
        msh_flag |= CS_FLAG_COMP_PVQ;
        break;

      }
      break;

    case CS_HODGE_TYPE_EPFD:

      msh_flag |= CS_FLAG_COMP_PEQ | CS_FLAG_COMP_DFQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOVB, connect);
      BFT_MALLOC(_in, connect->n_max_ebyc, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_edges; i++) result[i] = 0;

      switch (hodgep.algo) {

      case CS_HODGE_ALGO_OCS2:
        msh_flag |= CS_FLAG_COMP_EV | CS_FLAG_COMP_SEF;
        break;
      case CS_HODGE_ALGO_VORONOI:
        msh_flag |= CS_FLAG_COMP_SEF;
        break;

      default: /* Minimal requirement for other Hodge algorithms */
        break;

      }
      break;

    case CS_HODGE_TYPE_EDFP:

      msh_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOFB, connect);
      BFT_MALLOC(_in, connect->n_max_fbyc, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_faces; i++) result[i] = 0;

      break;

    case CS_HODGE_TYPE_FPED:

      msh_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOFB, connect);
      BFT_MALLOC(_in, connect->n_max_fbyc, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_faces; i++) result[i] = 0;

      break;

    case CS_HODGE_TYPE_VC:

      msh_flag |= CS_FLAG_COMP_PVQ | CS_FLAG_COMP_PFQ | CS_FLAG_COMP_HFQ |
         CS_FLAG_COMP_DEQ | CS_FLAG_COMP_EV | CS_FLAG_COMP_FEQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOVCB, connect);
      BFT_MALLOC(_in, connect->n_max_vbyc + 1, double);

#     pragma omp for
      for (cs_lnum_t i = 0; i < quant->n_vertices + quant->n_cells; i++)
        result[i] = 0;

      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of discrete Hodge operator", func_name);
    }

    const cs_flag_t  cell_flag = 0;
    if (pty_uniform) /* Get the value from the first cell */
      cs_hodge_set_property_value(0, t_eval, cell_flag, hodge);

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Retrieve the value of the property inside the current cell */

      if (!pty_uniform)
        cs_hodge_set_property_value_cw(cm, t_eval, cell_flag, hodge);

      /* Build the local discrete Hodge operator */

      compute(cm, hodge, cb);

      /* Define a local vector to multiply for the current cell */

      switch (hodgep.type) {

      case CS_HODGE_TYPE_VPCD:
        for (short int v = 0; v < cm->n_vc; v++)
          _in[v] = in_vals[cm->v_ids[v]];

        /* Local matrix-vector operation */

        cs_sdm_square_matvec(hodge->matrix, _in, cb->values);

        /* Assemble the resulting vector */

        for (short int v = 0; v < cm->n_vc; v++)
#         pragma omp atomic
          result[cm->v_ids[v]] += cb->values[v];
        break;

      case CS_HODGE_TYPE_EPFD:
        for (short int e = 0; e < cm->n_ec; e++)
          _in[e] = in_vals[cm->e_ids[e]];

        /* Local matrix-vector operation */

        cs_sdm_square_matvec(hodge->matrix, _in, cb->values);

        /* Assemble the resulting vector */

        for (short int e = 0; e < cm->n_ec; e++)
#         pragma omp atomic
          result[cm->e_ids[e]] += cb->values[e];
        break;

      case CS_HODGE_TYPE_FPED:
        for (short int f = 0; f < cm->n_fc; f++)
          _in[f] = in_vals[cm->f_ids[f]];

        /* Local matrix-vector operation */

        cs_sdm_square_matvec(hodge->matrix, _in, cb->values);

        /* Assemble the resulting vector */

        for (short int f = 0; f < cm->n_fc; f++)
#         pragma omp atomic
          result[cm->f_ids[f]] += cb->values[f];
        break;

      default:
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid type of discrete Hodge operator", func_name);

      } /* Hodge type */

    } /* Main loop on cells */

    BFT_FREE(_in);
    cs_cell_builder_free(&cb);
    cs_hodge_free(&hodge);

  } /* OpenMP Block */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute cellwise a discrete hodge operator in order to define
 *          a circulation array from a flux array
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 * \param[in]      quant     pointer to a cs_cdo_quantities_t structure
 * \param[in]      t_eval    time at which one performs the evaluation
 * \param[in]      hodgep    cs_hodge_param_t structure
 * \param[in]      pty       pointer to a cs_property_t structure or NULL
 * \param[in]      flux      vector to multiply with the discrete Hodge op.
 * \param[in, out] circul    array storing the resulting matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_hodge_circulation_from_flux(const cs_cdo_connect_t       *connect,
                               const cs_cdo_quantities_t    *quant,
                               cs_real_t                     t_eval,
                               const cs_hodge_param_t        hodgep,
                               const cs_property_t          *pty,
                               const cs_real_t               flux[],
                               cs_real_t                     circul[])
{
  if (flux == NULL)
    return;

  const char *func_name = __func__;

  if (circul == NULL) {
    bft_error(__FILE__, __LINE__, 0,
              "%s: Resulting vector must be allocated", __func__);
    return; /* Avoid a warning */
  }
  assert(connect != NULL && quant != NULL); /* Sanity checks */

# pragma omp parallel if (quant->n_cells > CS_THR_MIN)                  \
  shared(quant, connect, flux, t_eval, circul, pty, func_name)          \
  firstprivate(hodgep)
  {
#if defined(HAVE_OPENMP) /* Determine default number of OpenMP threads */
    int  t_id = omp_get_thread_num();
#else
    int  t_id = 0;
#endif

    const cs_adjacency_t  *c2f = connect->c2f;
    const cs_flag_t  cell_flag = 0;

    cs_eflag_t  msh_flag = 0;
    cs_cell_mesh_t  *cm = cs_cdo_local_get_cell_mesh(t_id);
    cs_cell_builder_t  *cb = NULL;
    double  *_fluxes = NULL;

    /* Each thread has its own pointer to a cs_hodge_t structure */
    cs_hodge_t  *hodge = cs_hodge_create(connect, pty, &hodgep, true, false);
    cs_hodge_compute_t  *compute = cs_hodge_get_func(func_name, hodgep);
    bool  pty_uniform = cs_property_is_uniform(pty);

    switch (hodgep.type) {

    /* Only this type of discrete Hodge operator makes sense for this
       operation */

    case CS_HODGE_TYPE_FPED:

      msh_flag |= CS_FLAG_COMP_PFQ | CS_FLAG_COMP_DEQ;
      cb = _cell_builder_create(CS_SPACE_SCHEME_CDOFB, connect);
      BFT_MALLOC(_fluxes, connect->n_max_fbyc, double);

      break;

    default:
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid type of discrete Hodge operator", func_name);
    }

    if (pty_uniform) /* Get the value from the first cell */
      cs_hodge_set_property_value(0, t_eval, cell_flag, hodge);

#   pragma omp for CS_CDO_OMP_SCHEDULE
    for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

      /* Set the local mesh structure for the current cell */

      cs_cell_mesh_build(c_id, msh_flag, connect, quant, cm);

      /* Retrieve the value of the property inside the current cell */

      if (!pty_uniform)
        cs_hodge_set_property_value_cw(cm, t_eval, cell_flag, hodge);

      /* Build the local discrete Hodge operator */

      compute(cm, hodge, cb);

      /* Define a local vector to multiply for the current cell */

      for (short int f = 0; f < cm->n_fc; f++)
        _fluxes[f] = flux[cm->f_ids[f]];
      cs_real_t  *_circ = circul + c2f->idx[c_id];

      /* Local matrix-vector operation */

      cs_sdm_square_matvec(hodge->matrix, _fluxes, _circ);

    } /* Main loop on cells */

    BFT_FREE(_fluxes);
    cs_cell_builder_free(&cb);
    cs_hodge_free(&hodge);

  } /* OpenMP Block */
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

#if defined(DEBUG) && !defined(NDEBUG) && CS_HODGE_DBG > 1
  if (fm->c_id % CS_HODGE_MODULO == 0) {
    cs_log_printf(CS_LOG_DEFAULT, " Surfacic Hodge op.   ");
    cs_sdm_dump(fm->f_id, NULL, NULL, hf);
  }
#endif
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
