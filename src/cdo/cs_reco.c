/*============================================================================
 * Routines to handle the reconstruction of fields
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_math.h"
#include "cs_cdo_scheme_geometry.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_reco.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro and structure definitions
 *============================================================================*/

/* Redefined the name of functions from cs_math to get shorter names */
#define _dp3  cs_math_3_dot_product

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct at cell centers and face centers a vertex-based field
 *         Linear interpolation. If p_crec and/or p_frec are not allocated, this
 *         done in this subroutine.
 *
 *  \param[in]      connect  pointer to the connectivity struct.
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      dof      pointer to the field of vtx-based DoFs
 *  \param[in, out] p_crec   reconstructed values at cell centers
 *  \param[in, out] p_frec   reconstructed values at face centers
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_conf_vtx_dofs(const cs_cdo_connect_t     *connect,
                      const cs_cdo_quantities_t  *quant,
                      const double               *dof,
                      double                     *p_crec[],
                      double                     *p_frec[])
{
  int  i, j, l, eid, vid;
  double  lef, lve;
  cs_real_3_t  uef, uve, cp;
  cs_quant_t  fq, eq;

  double  *crec = *p_crec, *frec = *p_frec;

  const cs_connect_index_t  *c2v = connect->c2v;
  const double  *dcv = quant->dcell_vol;
  const cs_sla_matrix_t  *f2e = connect->f2e;
  const cs_sla_matrix_t  *e2v = connect->e2v;

  if (dof == NULL)
    return;

  /* Allocate reconstruction arrays if necessary */
  if (crec == NULL)
    BFT_MALLOC(crec, quant->n_cells, double);
  if (frec == NULL)
    BFT_MALLOC(frec, quant->n_faces, double);

  /* Reconstruction at cell centers */
  for (i = 0; i < quant->n_cells; i++) {

    crec[i] = 0;
    for (j = c2v->idx[i]; j < c2v->idx[i+1]; j++)
      crec[i] += dcv[j]*dof[c2v->ids[j]];
    crec[i] /= quant->cell_vol[i];

  } /* End of loop on cells */

  /* Reconstruction at face centers */
  for (i = 0; i < quant->n_faces; i++) {

    frec[i] = 0;
    fq = quant->face[i];

    for (j = f2e->idx[i]; j < f2e->idx[i+1]; j++) {

      eid = f2e->col_id[j];
      eq = quant->edge[eid];
      cs_math_3_length_unitv(eq.center, fq.center, &lef, uef);

      for (l = e2v->idx[eid]; l < e2v->idx[eid+1]; l++) {

        vid = e2v->col_id[l];
        cs_math_3_length_unitv(quant->vtx_coord + 3*vid, eq.center, &lve, uve);
        cs_math_3_cross_product(uve, uef, cp);
        frec[i] += 0.5 * lve * lef * cs_math_3_norm(cp) * dof[vid];

      }

    } /* End of loop on face edges */
    frec[i] /= fq.meas;

  } /* End of loop on faces */

  /* Return pointers */
  *p_crec = crec;
  *p_frec = frec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct by a constant vector a field of edge-based DoFs
 *         in a volume surrounding an edge
 *
 *  \param[in]      cid     cell id
 *  \param[in]      e1_id   sub-volume related to this edge id
 *  \param[in]      c2e     cell -> edges connectivity
 *  \param[in]      quant   pointer to the additional quantities struct.
 *  \param[in]      dof     pointer to the field of edge-based DoFs
 *  \param[in, out] reco    value of the reconstructed field in this sub-volume
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cost_edge_dof(cs_lnum_t                    cid,
                      cs_lnum_t                    e1_id,
                      const cs_connect_index_t    *c2e,
                      const cs_cdo_quantities_t   *quant,
                      const double                *dof,
                      double                       reco[])
{
  int  i, k;
  double  inv_e1df1, e1df2;

  double  sum_vale1df2 = 0.0;
  double  t1[3] = {0., 0., 0.}, t2[3] = {0., 0., 0.}, t3[3] = {0., 0., 0.};

  const double  invvol = 1/quant->cell_vol[cid];
  const cs_quant_t e1q = quant->edge[e1_id]; /* Edge quantities */

  if (dof == NULL)
    return;

  for (i = c2e->idx[cid]; i < c2e->idx[cid+1]; i++) {

    const cs_dface_t  df2q = quant->dface[i];   /* Dual face quantities */
    const cs_lnum_t  e2_id = c2e->ids[i];
    const double  val = dof[e2_id];             /* Edge value */

    for (k = 0; k < 3; k++)
      t2[k] += val * df2q.vect[k];

    /* Better accuracy for the dot product with normalized vectors */
    e1df2  =  df2q.sface[0].meas * _dp3(e1q.unitv, df2q.sface[0].unitv);
    e1df2 +=  df2q.sface[1].meas * _dp3(e1q.unitv, df2q.sface[1].unitv);
    e1df2 *=  e1q.meas;
    sum_vale1df2 += e1df2 * val;

    if (e1_id == e2_id) {
      inv_e1df1 = 1./e1df2;
      for (k = 0; k < 3; k++) {
        t3[k] = inv_e1df1 * df2q.vect[k];
        t1[k] = val*t3[k];
      }
    }

  } /* End of loop on cell edges */

  /* Divide by cell volume */
  for (k = 0; k < 3; k++) {
    t2[k] *= invvol;
    t3[k] *= -invvol*sum_vale1df2;
    reco[k] = t1[k] + t2[k] + t3[k];
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at all cell centers from an array of values
 *         defined on primal vertices.
 *
 *  \param[in]      c2v      cell -> vertices connectivity
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      array    pointer to the array of values
 *  \param[in, out] val_xc   values of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_pv_at_cell_centers(const cs_connect_index_t    *c2v,
                           const cs_cdo_quantities_t   *quant,
                           const double                *array,
                           cs_real_t                   *val_xc)
{
  if (array == NULL)
    return;

  /* Sanity checks */
  assert(c2v != NULL && quant != NULL);

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    const double  invvol = 1/quant->cell_vol[c_id];
    const cs_real_t  *dcvol = quant->dcell_vol;

    cs_real_t  reco_val = 0;
    for (cs_lnum_t jv = c2v->idx[c_id]; jv < c2v->idx[c_id+1]; jv++) {

      const cs_lnum_t  v_id = c2v->ids[jv];

      reco_val += dcvol[jv] * array[v_id];

    } // Loop on cell vertices;

    val_xc[c_id] = invvol * reco_val;

  } // Loop on cells

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the cell center from an array of values
 *         defined on primal vertices.
 *
 *  \param[in]      c_id     cell id
 *  \param[in]      c2v      cell -> vertices connectivity
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      array    pointer to the array of values
 *  \param[in, out] val_xc   value of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_pv_at_cell_center(cs_lnum_t                    c_id,
                          const cs_connect_index_t    *c2v,
                          const cs_cdo_quantities_t   *quant,
                          const double                *array,
                          cs_real_t                   *val_xc)
{
  cs_real_t  reco_val = 0;

  if (array == NULL) {
    *val_xc = reco_val;
    return;
  }

  /* Sanity checks */
  assert(c2v != NULL && quant != NULL && c_id > -1);

  const double  invvol = 1/quant->cell_vol[c_id];
  const cs_real_t  *dcvol = quant->dcell_vol;

  for (cs_lnum_t jv = c2v->idx[c_id]; jv < c2v->idx[c_id+1]; jv++) {

    const cs_lnum_t  v_id = c2v->ids[jv];

    reco_val += dcvol[jv] * array[v_id];

  } // Loop on cell vertices;

  *val_xc = invvol * reco_val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the face center from an array of values
 *         defined on primal vertices.
 *
 *  \param[in]      f_id     face id (interior and border faces)
 *  \param[in]      connect  pointer to a cs_cdo_connect_t structure
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      pdi      pointer to the array of values
 *  \param[in, out] pdi_f    value of the reconstruction at the face center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_pf_from_pv(cs_lnum_t                     f_id,
                   const cs_cdo_connect_t       *connect,
                   const cs_cdo_quantities_t    *quant,
                   const double                 *pdi,
                   cs_real_t                    *pdi_f)
{
  *pdi_f = 0.;

  if (pdi == NULL)
    return;

  const cs_quant_t  qf = quant->face[f_id];
  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_sla_matrix_t  *e2v = connect->e2v;
  const cs_sla_matrix_t  *f2e = connect->f2e;

  for (cs_lnum_t i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

    const cs_lnum_t  e_id = f2e->col_id[i];
    const cs_lnum_t  shift_e = 2*e_id;
    const cs_lnum_t  v1_id = e2v->col_id[shift_e];
    const cs_lnum_t  v2_id = e2v->col_id[shift_e+1];
    const double  pdi_e = 0.5*(pdi[v1_id] + pdi[v2_id]);
    const cs_real_t  *xv1 = xyz + 3*v1_id;
    const cs_real_t  *xv2 = xyz + 3*v2_id;

    *pdi_f += pdi_e * cs_math_surftri(xv1, xv2, qf.center);

  } // Loop on face edges

  *pdi_f /= qf.meas;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a constant vector at the cell center from an array of
 *         values defined on dual faces lying inside each cell.
 *         This array is scanned thanks to the c2e connectivity.
 *
 *  \param[in]      c_id     cell id
 *  \param[in]      c2e      cell -> edges connectivity
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      array    pointer to the array of values
 *  \param[in, out] val_xc   value of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dfbyc_at_cell_center(cs_lnum_t                    c_id,
                             const cs_connect_index_t    *c2e,
                             const cs_cdo_quantities_t   *quant,
                             const double                *array,
                             cs_real_3_t                  val_xc)
{
  /* Initialization */
  val_xc[0] = val_xc[1] = val_xc[2] = 0.;

  if (array == NULL)
    return;

  const double  invvol = 1/quant->cell_vol[c_id];

  for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++) {

    const cs_lnum_t  e_id = c2e->ids[j];
    const cs_quant_t  peq = quant->edge[e_id];
    const cs_real_t  edge_contrib = array[j]*peq.meas;

    for (int k = 0; k < 3; k++)
      val_xc[k] += edge_contrib * peq.unitv[k];

  } // Loop on cell edges

  for (int k = 0; k < 3; k++)
    val_xc[k] *= invvol;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a constant vector inside the cell c.
 *         array is scanned thanks to the c2e connectivity. Pointer is already
 *         located at the beginning of the cell sequence.
 *         Reconstruction used is based on DGA (stabilization = 1/d where d is
 *         the space dimension)
 *
 *  \param[in]      cm        pointer to a cs_cell_mesh_t structure
 *  \param[in]      array     local pointer to the array of values
 *  \param[in, out] val_c     value of the reconstructed vector in the cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dfbyc_in_cell(const cs_cell_mesh_t        *cm,
                      const cs_real_t             *array,
                      cs_real_3_t                  val_c)
{
  /* Initialization */
  val_c[0] = val_c[1] = val_c[2] = 0.;

  if (array == NULL)
    return;

  /* Sanity check */
  assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_PEQ));

  const double  invvol = 1/cm->vol_c;

  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_quant_t  peq = cm->edge[e];
    const cs_real_t  edge_contrib = array[e]*peq.meas;

    for (int k = 0; k < 3; k++)
      val_c[k] += edge_contrib * peq.unitv[k];

  } // Loop on cell edges

  for (int k = 0; k < 3; k++)
    val_c[k] *= invvol;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a constant vector inside pec which is a volume
 *         surrounding the edge e inside the cell c.
 *         array is scanned thanks to the c2e connectivity. Pointer is already
 *         located at the beginning of the cell sequence.
 *         Reconstruction used is based on DGA (stabilization = 1/d where d is
 *         the space dimension)
 *
 *  \param[in]      cm        pointer to a cs_cell_mesh_t structure
 *  \param[in]      e         local edge id
 *  \param[in]      array     local pointer to the array of values
 *  \param[in, out] val_pec   value of the reconstruction in pec
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_dfbyc_in_pec(const cs_cell_mesh_t        *cm,
                     short int                    e,
                     const cs_real_t             *array,
                     cs_real_3_t                  val_pec)
{
  /* Initialize values */
  val_pec[0] = val_pec[1] = val_pec[2] = 0.;

  if (array == NULL)
    return;

  /* Sanity check */
  assert(cs_test_flag(cm->flag, CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ));

  cs_real_3_t  val_c = {0., 0., 0.};
  /* Compute val_c */
  for (short int _e = 0; _e < cm->n_ec; _e++) {

    const cs_quant_t  _peq = cm->edge[_e];

    for (int k = 0; k < 3; k++)
      val_c[k] += array[_e] * _peq.meas * _peq.unitv[k];

  } // Loop on cell edges

  const double  invvol = 1/cm->vol_c;
  /* Compute the constency part related to this cell */
  for (int k = 0; k < 3; k++)
    val_c[k] *= invvol;

  /* Compute the reconstruction inside pec */
  const cs_quant_t  peq = cm->edge[e];
  const cs_nvec3_t  dfq = cm->dface[e];
  const double ecoef = (array[e] - dfq.meas * _dp3(dfq.unitv, val_c))
    / (dfq.meas * _dp3(dfq.unitv, peq.unitv));

  for (int k = 0; k < 3; k++)
    val_pec[k] = val_c[k] + ecoef * peq.unitv[k];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at the cell center of the gradient of a field
 *         defined on primal vertices.
 *
 * \param[in]      c_id     cell id
 * \param[in]      connect  pointer to a cs_cdo_connect_t structure
 * \param[in]      quant    pointer to the additional quantities struct.
 * \param[in]      pdi      pointer to the array of values
 * \param[in, out] val_xc   value of the reconstructed gradient at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_grd_cell_from_pv(cs_lnum_t                    c_id,
                         const cs_cdo_connect_t      *connect,
                         const cs_cdo_quantities_t   *quant,
                         const double                *pdi,
                         cs_real_t                    val_xc[])
{
  val_xc[0] = val_xc[1] = val_xc[2] = 0.;

  if (pdi == NULL)
    return;

  const cs_connect_index_t  *c2e = connect->c2e;
  const cs_sla_matrix_t  *e2v = connect->e2v;

  for (cs_lnum_t i = c2e->idx[c_id]; i < c2e->idx[c_id+1]; i++) {

    const cs_lnum_t  e_id = c2e->ids[i];
    const cs_lnum_t  shift_e = 2*e_id;
    const cs_lnum_t  v1_id = e2v->col_id[shift_e];
    const short int  sgn_v1 = e2v->sgn[shift_e];
    const cs_lnum_t  v2_id = e2v->col_id[shift_e+1];
    const short int  sgn_v2 = e2v->sgn[shift_e+1];
    const double  gdi_e = sgn_v1*pdi[v1_id] + sgn_v2*pdi[v2_id];
    const cs_dface_t  dfq = quant->dface[i];  /* Dual face quantities */

    for (int k = 0; k < 3; k++)
      val_xc[k] += gdi_e*dfq.vect[k];

  } // Loop on cell edges

  /* Divide by cell volume */
  const double  invvol = 1/quant->cell_vol[c_id];
  for (int k = 0; k < 3; k++)
    val_xc[k] *= invvol;

}
/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct at the cell center a field of edge-based DoFs
 *
 *  \param[in]      cid     cell id
 *  \param[in]      c2e     cell -> edges connectivity
 *  \param[in]      quant   pointer to the additional quantities struct.
 *  \param[in]      dof     pointer to the field of edge-based DoFs
 *  \param[in, out] reco    value of the reconstructed field at cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_ccen_edge_dof(cs_lnum_t                    cid,
                      const cs_connect_index_t    *c2e,
                      const cs_cdo_quantities_t   *quant,
                      const double                *dof,
                      double                       reco[])
{
  if (dof == NULL)
    return;

  /* Initialize value */
  reco[0] = reco[1] = reco[2] = 0.0;

  for (cs_lnum_t i = c2e->idx[cid]; i < c2e->idx[cid+1]; i++) {

    const cs_dface_t  dfq = quant->dface[i];  /* Dual face quantities */
    const double  val = dof[c2e->ids[i]];     /* Edge value */

    for (int k = 0; k < 3; k++)
      reco[k] += val*dfq.vect[k];

  } /* End of loop on cell edges */

  /* Divide by cell volume */
  const double  invvol = 1/quant->cell_vol[cid];
  for (int k = 0; k < 3; k++)
    reco[k] *= invvol;

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct at each cell center a field of edge-based DoFs
 *
 *  \param[in]      connect   pointer to the connectivity struct.
 *  \param[in]      quant     pointer to the additional quantities struct.
 *  \param[in]      dof       pointer to the field of edge-based DoFs
 *  \param[in, out] p_ccrec   pointer to the reconstructed values
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_ccen_edge_dofs(const cs_cdo_connect_t     *connect,
                       const cs_cdo_quantities_t  *quant,
                       const double               *dof,
                       double                     *p_ccrec[])
{
  assert(connect->c2e != NULL); /* Sanity check */

  double  *ccrec = *p_ccrec;

  if (dof == NULL)
    return;

  /* Allocate reconstructed vector field at each cell bary. */
  if (ccrec == NULL)
    BFT_MALLOC(ccrec, 3*quant->n_cells, double);

# pragma omp parallel for if (quant->n_cells > CS_THR_MIN)
  for (int c_id = 0; c_id < quant->n_cells; c_id++)
    cs_reco_ccen_edge_dof(c_id,
                          connect->c2e,
                          quant,
                          dof,
                          ccrec + 3*c_id);

  /* Return pointer */
  *p_ccrec = ccrec;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value of a scalar potential at a point inside a cell
 *         The scalar potential has DoFs located at primal vertices
 *
 * \param[in]      cm             pointer to a cs_cell_mesh_t structure
 * \param[in]      pdi            array of DoF values at vertices
 * \param[out]     cell_gradient  gradient inside the cell
 *
 * \return the value of the reconstructed potential at the evaluation point
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cell_grad_from_scalar_pv(const cs_cell_mesh_t    *cm,
                                    const cs_real_t          pdi[],
                                    cs_real_t               *cell_gradient)
{
  /* Sanity checks */
  assert(cm != NULL && pdi != NULL);
  assert(cs_test_flag(cm->flag,
                      CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_DFQ));

  /* Reconstruct a constant gradient inside the current cell */
  cell_gradient[0] = cell_gradient[1] = cell_gradient[2] = 0;
  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_lnum_t  v0 = cm->v_ids[cm->e2v_ids[2*e]];
    const cs_lnum_t  v1 = cm->v_ids[cm->e2v_ids[2*e+1]];
    const cs_real_t  coeff = cm->e2v_sgn[e]*(pdi[v0]-pdi[v1])*cm->dface[e].meas;
    for (int k = 0; k < 3; k++)
      cell_gradient[k] += coeff * cm->dface[e].unitv[k];

  } // Loop on cell edges

  const cs_real_t  invcell = 1./cm->vol_c;
  for (int k = 0; k < 3; k++) cell_gradient[k] *= invcell;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value of a scalar potential at a point inside a cell
 *         The scalar potential has DoFs located at primal vertices
 *
 * \param[in]      cm           pointer to a cs_cell_mesh_t structure
 * \param[in]      pdi          array of DoF values at vertices
 * \param[in]      length_xcxp  lenght of the segment [x_c, x_p]
 * \param[in]      unitv_xcxp   unitary vector pointed from x_c to x_p
 * \param[in, out] wbuf         pointer to a temporary buffer
 *
 * \return the value of the reconstructed potential at the evaluation point
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_reco_cw_scalar_pv_inside_cell(const cs_cell_mesh_t    *cm,
                                 const cs_real_t          pdi[],
                                 const cs_real_t          length_xcxp,
                                 const cs_real_t          unitv_xcxp[],
                                 cs_real_t                wbuf[])
{
  /* Sanity checks */
  assert(cm != NULL && pdi != NULL && wbuf != NULL);
  assert(cs_test_flag(cm->flag,
                      CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_DFQ));

  cs_real_t  *_pv = wbuf;  // Local value of the potential field

  /* Reconstruct the value at the cell center */
  cs_real_t  pc = 0.;
  for (short int v = 0; v < cm->n_vc; v++) {
    _pv[v] = pdi[cm->v_ids[v]];
    pc += cm->wvc[v] * _pv[v];
  }

  /* Reconstruct a constant gradient inside the current cell */
  cs_real_3_t  gcell = {0., 0., 0.};
  for (short int e = 0; e < cm->n_ec; e++) {
    const cs_real_t  ge =
      cm->e2v_sgn[e]*( _pv[cm->e2v_ids[2*e]] - _pv[cm->e2v_ids[2*e+1]]);
    const cs_real_t  coef_e = ge * cm->dface[e].meas;
    for (int k = 0; k < 3; k++)
      gcell[k] += coef_e * cm->dface[e].unitv[k];
  }

  const cs_real_t  invcell = 1./cm->vol_c;
  for (int k = 0; k < 3; k++) gcell[k] *= invcell;

  /* Evaluation at the given point */
  cs_real_t  p_rec = pc;
  p_rec += length_xcxp * cs_math_3_dot_product(gcell, unitv_xcxp);

  return  p_rec;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
