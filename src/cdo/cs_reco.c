/*============================================================================
 * Routines to handle the reconstruction of fields
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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
#include "cs_scheme_geometry.h"

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
  double  *crec = *p_crec, *frec = *p_frec;

  const cs_adjacency_t  *c2v = connect->c2v;
  const double  *dcv = quant->dcell_vol;
  const cs_adjacency_t  *f2e = connect->f2e;
  const cs_adjacency_t  *e2v = connect->e2v;

  if (dof == NULL)
    return;

  /* Allocate reconstruction arrays if necessary */
  if (crec == NULL)
    BFT_MALLOC(crec, quant->n_cells, double);
  if (frec == NULL)
    BFT_MALLOC(frec, quant->n_faces, double);

  /* Reconstruction at cell centers */
  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {

    crec[c_id] = 0;
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++)
      crec[c_id] += dcv[j]*dof[c2v->ids[j]];
    crec[c_id] /= quant->cell_vol[c_id];

  }

  /* Reconstruction at face centers */
  for (cs_lnum_t f_id = 0; f_id < quant->n_faces; f_id++) {

    const cs_real_t  *xf = cs_quant_get_face_center(f_id, quant);

    frec[f_id] = 0;
    double  f_surf = 0.;
    for (cs_lnum_t j = f2e->idx[f_id]; j < f2e->idx[f_id+1]; j++) {

      const cs_lnum_t  e_id = f2e->ids[j];
      const cs_lnum_t  v1_id = e2v->ids[2*e_id];
      const cs_lnum_t  v2_id = e2v->ids[2*e_id+1];
      const cs_real_t  *xv1 = quant->vtx_coord + 3*v1_id;
      const cs_real_t  *xv2 = quant->vtx_coord + 3*v2_id;

      cs_real_3_t  xe;
      for (int k = 0; k < 3; k++)
        xe[k] = 0.5 * (xv1[k] + xv2[k]);

      double  lef, lve;
      cs_real_3_t  uef, uve, cp;
      cs_math_3_length_unitv(xe, xf, &lef, uef);
      cs_math_3_length_unitv(xv1, xv2, &lve, uve);
      cs_math_3_cross_product(uve, uef, cp);

      const double  tef = 0.5 * lve * lef * cs_math_3_norm(cp);

      f_surf += tef;
      frec[f_id] += 0.5 * tef * (dof[v1_id] + dof[v2_id]);

    } /* End of loop on face edges */

    frec[f_id] /= f_surf;

  } /* End of loop on faces */

  /* Return pointers */
  *p_crec = crec;
  *p_frec = frec;
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
cs_reco_pv_at_cell_centers(const cs_adjacency_t        *c2v,
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

    } /* Loop on cell vertices; */

    val_xc[c_id] = invvol * reco_val;

  } /* Loop on cells */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct the value at all cell centers from an array of values
 *         defined on primal vertices.
 *         Case of vector-valued fields.
 *
 *  \param[in]      c2v      cell -> vertices connectivity
 *  \param[in]      quant    pointer to the additional quantities struct.
 *  \param[in]      array    pointer to the array of values
 *  \param[in, out] val_xc   values of the reconstruction at the cell center
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_vect_pv_at_cell_centers(const cs_adjacency_t        *c2v,
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

    cs_real_t  reco_val[3] = {0, 0, 0};
    for (cs_lnum_t jv = c2v->idx[c_id]; jv < c2v->idx[c_id+1]; jv++) {

      const cs_real_t  *_array = array + 3*c2v->ids[jv];
      const cs_real_t  vc_vol = dcvol[jv];

      reco_val[0] += vc_vol * _array[0];
      reco_val[1] += vc_vol * _array[1];
      reco_val[2] += vc_vol * _array[2];

    } /* Loop on cell vertices */

    for (int k = 0; k < 3; k++)
      val_xc[3*c_id+k] = invvol * reco_val[k];

  } /* Loop on cells */

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
                          const cs_adjacency_t        *c2v,
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

  } /* Loop on cell vertices; */

  *val_xc = invvol * reco_val;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reconstruct a vector-valued array at vertices from a vector-valued
 *         array at cells.
 *
 *  \param[in]      c2v       cell -> vertices connectivity
 *  \param[in]      quant     pointer to the additional quantities struct.
 *  \param[in]      val       pointer to the array of values
 *  \param[in, out] reco_val  values of the reconstruction at vertices
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_vect_pv_from_pc(const cs_adjacency_t        *c2v,
                        const cs_cdo_quantities_t   *quant,
                        const double                *val,
                        cs_real_t                   *reco_val)
{
  if (val == NULL || reco_val == NULL)
    return;

  memset(reco_val, 0, 3*quant->n_vertices*sizeof(cs_real_t));

  for (cs_lnum_t c_id = 0; c_id < quant->n_cells; c_id++) {
    for (cs_lnum_t j = c2v->idx[c_id]; j < c2v->idx[c_id+1]; j++) {

      const cs_real_t  vc_vol = quant->dcell_vol[j];
      cs_real_t  *rval = reco_val + 3*c2v->ids[j];

      rval[0] += vc_vol * val[3*c_id];
      rval[1] += vc_vol * val[3*c_id + 1];
      rval[2] += vc_vol * val[3*c_id + 2];

    }
  } /* Loop on cells */

  cs_real_t  *dual_vol = NULL;
  BFT_MALLOC(dual_vol, quant->n_vertices, cs_real_t);
  cs_cdo_quantities_compute_dual_volumes(quant, c2v, dual_vol);

# pragma omp parallel for if (quant->n_vertices > CS_THR_MIN)
  for (cs_lnum_t v_id = 0; v_id < quant->n_vertices; v_id++) {
    const cs_real_t  invvol = 1./dual_vol[v_id];
    for (int k = 0; k < 3; k++)
      reco_val[3*v_id+k] *= invvol;
  }

  BFT_FREE(dual_vol);
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

  const cs_real_t  *xf = cs_quant_get_face_center(f_id, quant);
  const cs_real_t  *xyz = quant->vtx_coord;
  const cs_adjacency_t  *e2v = connect->e2v;
  const cs_adjacency_t  *f2e = connect->f2e;

  double f_surf = 0.;
  for (cs_lnum_t i = f2e->idx[f_id]; i < f2e->idx[f_id+1]; i++) {

    const cs_lnum_t  e_id = f2e->ids[i];
    const cs_lnum_t  shift_e = 2*e_id;
    const cs_lnum_t  v1_id = e2v->ids[shift_e];
    const cs_lnum_t  v2_id = e2v->ids[shift_e+1];
    const double  pdi_e = 0.5*(pdi[v1_id] + pdi[v2_id]);
    const cs_real_t  *xv1 = xyz + 3*v1_id;
    const cs_real_t  *xv2 = xyz + 3*v2_id;
    const cs_real_t  tef = cs_math_surftri(xv1, xv2, xf);

    f_surf += tef;
    *pdi_f += pdi_e * tef;

  } /* Loop on face edges */

  *pdi_f /= f_surf;
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
                             const cs_adjacency_t        *c2e,
                             const cs_cdo_quantities_t   *quant,
                             const double                *array,
                             cs_real_3_t                  val_xc)
{
  /* Initialization */
  val_xc[0] = val_xc[1] = val_xc[2] = 0.;

  if (array == NULL)
    return;

  for (cs_lnum_t j = c2e->idx[c_id]; j < c2e->idx[c_id+1]; j++) {

    const cs_real_t  *e_vect = quant->edge_vector + 3*c2e->ids[j];
    for (int k = 0; k < 3; k++)
      val_xc[k] += array[j] * e_vect[k];

  } /* Loop on cell edges */

  const double  invvol = 1/quant->cell_vol[c_id];
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
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PEQ));

  const double  invvol = 1/cm->vol_c;

  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_quant_t  peq = cm->edge[e];
    const cs_real_t  edge_contrib = array[e]*peq.meas;

    for (int k = 0; k < 3; k++)
      val_c[k] += edge_contrib * peq.unitv[k];

  } /* Loop on cell edges */

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
  assert(cs_flag_test(cm->flag, CS_CDO_LOCAL_PEQ | CS_CDO_LOCAL_DFQ));

  cs_real_3_t  val_c = {0., 0., 0.};
  /* Compute val_c */
  for (short int _e = 0; _e < cm->n_ec; _e++) {

    const cs_quant_t  _peq = cm->edge[_e];

    for (int k = 0; k < 3; k++)
      val_c[k] += array[_e] * _peq.meas * _peq.unitv[k];

  } /* Loop on cell edges */

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
                      const cs_adjacency_t        *c2e,
                      const cs_cdo_quantities_t   *quant,
                      const double                *dof,
                      double                       reco[])
{
  if (dof == NULL)
    return;

  /* Initialize value */
  reco[0] = reco[1] = reco[2] = 0.0;

  for (cs_lnum_t i = c2e->idx[cid]; i < c2e->idx[cid+1]; i++) {

    /* Dual face quantities */
    const cs_real_t  *sf0 = quant->sface_normal + 6*i;
    const cs_real_t  *sf1 = sf0 + 3;
    const double  val = dof[c2e->ids[i]];     /* Edge value */

    for (int k = 0; k < 3; k++)
      reco[k] += val * (sf0[k] + sf1[k]);

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
cs_reco_grad_cell_from_pv(cs_lnum_t                    c_id,
                          const cs_cdo_connect_t      *connect,
                          const cs_cdo_quantities_t   *quant,
                          const cs_real_t             *pdi,
                          cs_real_t                    val_xc[])
{
  val_xc[0] = val_xc[1] = val_xc[2] = 0.;

  if (pdi == NULL)
    return;

  const cs_adjacency_t  *c2e = connect->c2e;
  const cs_adjacency_t  *e2v = connect->e2v;

  for (cs_lnum_t i = c2e->idx[c_id]; i < c2e->idx[c_id+1]; i++) {

    const cs_lnum_t  shift_e = 2*c2e->ids[i];
    const short int  sgn_v1 = e2v->sgn[shift_e];
    const cs_real_t  pv1 = pdi[e2v->ids[shift_e]];
    const cs_real_t  pv2 = pdi[e2v->ids[shift_e+1]];
    const cs_real_t  gdi_e = sgn_v1*(pv1 - pv2);

    /* Dual face quantities */
    const cs_real_t  *sf0 = quant->sface_normal + 6*i;
    const cs_real_t  *sf1 = sf0 + 3;

    for (int k = 0; k < 3; k++)
      val_xc[k] += gdi_e * (sf0[k] + sf1[k]);

  } /* Loop on cell edges */

  /* Divide by cell volume */
  const double  invvol = 1/quant->cell_vol[c_id];
  for (int k = 0; k < 3; k++)
    val_xc[k] *= invvol;

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
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_DFQ));

  /* Reconstruct a constant gradient inside the current cell */
  cell_gradient[0] = cell_gradient[1] = cell_gradient[2] = 0;
  for (short int e = 0; e < cm->n_ec; e++) {

    const cs_lnum_t  v0 = cm->v_ids[cm->e2v_ids[2*e]];
    const cs_lnum_t  v1 = cm->v_ids[cm->e2v_ids[2*e+1]];
    const cs_real_t  coeff = cm->e2v_sgn[e]*(pdi[v0]-pdi[v1])*cm->dface[e].meas;
    for (int k = 0; k < 3; k++)
      cell_gradient[k] += coeff * cm->dface[e].unitv[k];

  } /* Loop on cell edges */

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
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PVQ | CS_CDO_LOCAL_EV | CS_CDO_LOCAL_DFQ));

  cs_real_t  *_pv = wbuf;  /* Local value of the potential field */

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
/*!
 * \brief   Compute the weighted (by volume) gradient inside a given primal
 *          cell for the related vertices.
 *          Use the WBS algo. for approximating the gradient.
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices + cell
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] vgrd     gradient at vertices inside this cell
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_vgrd_wbs_from_pvc(const cs_cell_mesh_t   *cm,
                             const cs_real_t        *pot,
                             cs_cell_builder_t      *cb,
                             cs_real_t              *vgrd)
{
  /* Sanity checks */
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV  | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_HFQ));

  cs_real_3_t  grd_c, grd_v1, grd_v2;

  /* Temporary buffers */
  cs_real_3_t  *u_vc = cb->vectors;
  double  *l_vc = cb->values;

  const double  *p_v = pot;
  const double  p_c = pot[cm->n_vc];

  /* Reset local fluxes */
  for (int i = 0; i < 3*cm->n_vc; i++) vgrd[i] = 0.;

  /* Store segments xv --> xc for this cell */
  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];
    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute geometrical quantities for the current face:
       - the weighting of each triangle defined by a base e and an apex f
       - volume of each sub-tetrahedron pef_c
       - the gradient of the Lagrange function related xc in p_{f,c} */
    const cs_real_t  ohf = -cm->f_sgn[f]/cm->hfc[f];
    for (int k = 0; k < 3; k++) grd_c[k] = ohf * pfq.unitv[k];

    /* Compute the reconstructed value of the potential at p_f */
    double  p_f = 0.;
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  ee = 2*cm->f2e_ids[i];

      p_f += cm->tef[i]*(  p_v[cm->e2v_ids[ee]]      /* p_v1 */
                         + p_v[cm->e2v_ids[ee+1]] ); /* p_v2 */
    }
    p_f *= 0.5/pfq.meas;

    const double  dp_cf = p_c - p_f;

    /* Loop on face edges to scan p_{ef,c} subvolumes */
    const cs_real_t  hf_coef = cs_math_onethird * cm->hfc[f];
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  ee = 2*cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[ee];
      const short int  v2 = cm->e2v_ids[ee+1];

      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                        grd_v1, grd_v2);

      /* Gradient of the Lagrange function related to a face.
         grd_f = -(grd_c + grd_v1 + grd_v2)
         This formula is a consequence of the Partition of the Unity.
         This yields the following formula for grd(Lv^conf)|_p_{ef,c} */
      const cs_real_t  sefv_vol = 0.5 * hf_coef * cm->tef[i]; /* 0.5 pef_vol */
      cs_real_3_t  _grd;
      for (int k = 0; k < 3; k++) {
        _grd[k] = sefv_vol * ( dp_cf          * grd_c[k]  +
                              (p_v[v1] - p_f) * grd_v1[k] +
                              (p_v[v2] - p_f) * grd_v2[k]);
        vgrd[3*v1 + k] += _grd[k];
        vgrd[3*v2 + k] += _grd[k];
      }

    } /* Loop on face edges */

  } /* Loop on cell faces */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the mean value of a gradient inside a given primal cell.
 *          Use the WBS algo. for approximating the gradient.
 *          The computation takes into account a subdivision into tetrahedra of
 *          the current cell based on p_{ef,c}
 *
 * \param[in]      cm       pointer to a cs_cell_mesh_t structure
 * \param[in]      pot      values of the potential fields at vertices + cell
 * \param[in, out] cb       auxiliary structure for computing the flux
 * \param[in, out] cgrd     mean-value of the cell gradient
 */
/*----------------------------------------------------------------------------*/

void
cs_reco_cw_cgrd_wbs_from_pvc(const cs_cell_mesh_t   *cm,
                             const cs_real_t        *pot,
                             cs_cell_builder_t      *cb,
                             cs_real_t              *cgrd)
{
  /* Sanity checks */
  assert(cs_flag_test(cm->flag,
                      CS_CDO_LOCAL_PV  | CS_CDO_LOCAL_PFQ | CS_CDO_LOCAL_DEQ |
                      CS_CDO_LOCAL_FEQ | CS_CDO_LOCAL_EV  | CS_CDO_LOCAL_HFQ));

  cs_real_3_t  grd_c, grd_v1, grd_v2;

  /* Temporary buffers */
  cs_real_3_t  *u_vc = cb->vectors;
  double  *l_vc = cb->values;

  cgrd[0] = cgrd[1] = cgrd[2] = 0.;

  const double  *p_v = pot;
  const double  p_c = pot[cm->n_vc];

  /* Store segments xv --> xc for this cell */
  for (short int v = 0; v < cm->n_vc; v++)
    cs_math_3_length_unitv(cm->xc, cm->xv + 3*v, l_vc + v, u_vc[v]);

  /* Loop on cell faces */
  for (short int f = 0; f < cm->n_fc; f++) {

    const cs_quant_t  pfq = cm->face[f];
    const cs_nvec3_t  deq = cm->dedge[f];

    /* Compute geometrical quantities for the current face:
       - the weighting of each triangle defined by a base e and an apex f
       - volume of each sub-tetrahedron pef_c
       - the gradient of the Lagrange function related xc in p_{f,c} */
    const cs_real_t  ohf = -cm->f_sgn[f]/cm->hfc[f];
    for (int k = 0; k < 3; k++) grd_c[k] = ohf * pfq.unitv[k];

    /* Compute the reconstructed value of the potential at p_f */
    double  p_f = 0.;
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  ee = 2*cm->f2e_ids[i];

      p_f += cm->tef[i]*(  p_v[cm->e2v_ids[ee]]      /* p_v1 */
                         + p_v[cm->e2v_ids[ee+1]] ); /* p_v2 */
    }
    p_f *= 0.5/pfq.meas;

    const double  dp_cf = p_c - p_f;

    /* Loop on face edges to scan p_{ef,c} subvolumes */
    const cs_real_t  hf_coef = cs_math_onethird * cm->hfc[f];
    for (int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

      const short int  ee = 2*cm->f2e_ids[i];
      const short int  v1 = cm->e2v_ids[ee];
      const short int  v2 = cm->e2v_ids[ee+1];

      cs_compute_grd_ve(v1, v2, deq, (const cs_real_t (*)[3])u_vc, l_vc,
                        grd_v1, grd_v2);

      /* Gradient of the Lagrange function related to a face.
         grd_f = -(grd_c + grd_v1 + grd_v2)
         This formula is a consequence of the Partition of the Unity.
         This yields the following formula for grd(Lv^conf)|_p_{ef,c} */
      const cs_real_t  pefv_vol = hf_coef * cm->tef[i];
      for (int k = 0; k < 3; k++)
        cgrd[k] += pefv_vol * ( dp_cf          * grd_c[k]  +
                               (p_v[v1] - p_f) * grd_v1[k] +
                               (p_v[v2] - p_f) * grd_v2[k]);

    } /* Loop on face edges */

  } /* Loop on cell faces */

  const double  invvol = 1/cm->vol_c;
  for (int k = 0; k < 3; k++) cgrd[k] *= invvol;
}

/*----------------------------------------------------------------------------*/

#undef _dp3

END_C_DECLS
