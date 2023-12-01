/*============================================================================
 * Face viscosity
 *============================================================================*/

/* This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_math.h"
#include "cs_blas.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_internal_coupling.h"
#include "cs_log.h"
#include "cs_mesh.h"
#include "cs_field.h"
#include "cs_field_default.h"
#include "cs_field_pointer.h"
#include "cs_gradient.h"
#include "cs_ext_neighborhood.h"
#include "cs_mesh_quantities.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_porous_model.h"
#include "cs_prototypes.h"
#include "cs_timer.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_face_viscosity.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file  cs_face_viscosity.c
 *
 *  \brief Face viscosity.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_viscosity
 *----------------------------------------------------------------------------*/

void CS_PROCF (viscfa, VISCFA)
(
 const int  *const   visc_mean_type,
 cs_real_t           c_visc[],
 cs_real_t           i_visc[],
 cs_real_t           b_visc[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_face_viscosity(m,
                    fvq,
                    *visc_mean_type,
                    c_visc,
                    i_visc,
                    b_visc);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_anisotropic_viscosity_vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (vistnv, VISTNV)
(
 const int  *const   visc_mean_type,
 cs_real_6_t         c_visc[],
 cs_real_33_t        i_visc[],
 cs_real_t           b_visc[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_face_anisotropic_viscosity_vector(m,
                                       fvq,
                                       *visc_mean_type,
                                       c_visc,
                                       i_visc,
                                       b_visc);
}

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_anisotropic_viscosity_scalar
 *----------------------------------------------------------------------------*/

void CS_PROCF (vitens, VITENS)
(
 cs_real_6_t         c_visc[],
 const int    *const iwarnp,
 cs_real_2_t         weighf[],
 cs_real_t           weighb[],
 cs_real_t           i_visc[],
 cs_real_t           b_visc[]
)
{
  const cs_mesh_t  *m = cs_glob_mesh;
  cs_mesh_quantities_t  *fvq = cs_glob_mesh_quantities;

  cs_face_anisotropic_viscosity_scalar(m,
                                       fvq,
                                       c_visc,
                                       *iwarnp,
                                       weighf,
                                       weighb,
                                       i_visc,
                                       b_visc);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Computes the secondary viscosity contribution \f$\kappa
 * -\dfrac{2}{3} \mu\f$ in order to compute:
 * \f[
 * \grad\left( (\kappa -\dfrac{2}{3} \mu) \trace( \gradt(\vect{u})) \right)
 * \f]
 * with:
 *   - \f$ \mu = \mu_{laminar} + \mu_{turbulent} \f$
 *   - \f$ \kappa \f$ is the volume viscosity (generally zero)
 *
 * \remark
 * In LES, the tensor
 * \f$\overline{\left(\vect{u}-\overline{\vect{u}}\right)\otimes\left(\vect{u}
 *-\overline{\vect{u}}\right)}\f$
 * is modeled by \f$\mu_t \overline{\tens{S}}\f$
 * and not by
 * \f$\mu_t\overline{\tens{S}}-\dfrac{2}{3}\mu_t
 * \trace\left(\overline{\tens{S}}\right)\tens{1}+\dfrac{2}{3}k\tens{1}\f$
 * so that no term
 * \f$\mu_t \dive \left(\overline{\vect{u}}\right)\f$ is needed.
 *
 * Please refer to the
 * <a href="../../theory.pdf#visecv"><b>visecv</b></a> section
 * of the theory guide for more informations.
 *
 * \param[in,out] secvif        lambda*surface at interior faces
 * \param[in,out] secvib        lambda*surface at boundary faces
 */
/*----------------------------------------------------------------------------*/

void
cs_face_viscosity_secondary(cs_real_t  secvif[],
                            cs_real_t  secvib[])
{
  const cs_mesh_t *mesh = cs_glob_mesh;
  const cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  const cs_lnum_t n_cells_ext = mesh->n_cells_with_ghosts;
  const cs_lnum_t n_cells = mesh->n_cells;
  const cs_lnum_t n_b_faces = mesh->n_b_faces;
  const cs_lnum_t n_i_faces = mesh->n_i_faces;

  const cs_lnum_2_t *i_face_cells
    = (const cs_lnum_2_t *)mesh->i_face_cells;
  const cs_lnum_t *b_face_cells = mesh->b_face_cells;
  const cs_real_t *weight = fvq->weight;

  const int itytur = cs_glob_turb_model->itytur;

  /* Initialization
     -------------- */

  /* Allocate temporary arrays */
  cs_real_t *secvis;
  BFT_MALLOC(secvis, n_cells_ext, cs_real_t);

  cs_field_t *vel = CS_F_(vel);
  cs_equation_param_t *eqp_vel = cs_field_get_equation_param(vel);
  const cs_real_t *viscl = CS_F_(mu)->val;
  const cs_real_t *visct = CS_F_(mu_t)->val;

  cs_real_t *cpro_viscv = NULL;
  cs_field_t *f_viscv = cs_field_by_name_try("volume_viscosity");

  if (f_viscv != NULL)
    cpro_viscv = cs_field_by_name_try("volume_viscosity")->val;

  /* Time extrapolation ? */
  int key_t_ext_id = cs_field_key_id("time_extrapolated");

  /* Computation of the second viscosity: lambda = K -2/3 mu
     ------------------------------------------------------- */

  /* For order 2 in time, everything should be taken in n...*/

  const cs_real_t d2s3m = -2.0/3.0;

  const int iviext = cs_field_get_key_int(CS_F_(mu), key_t_ext_id);
  const int iviext_t = cs_field_get_key_int(CS_F_(mu_t), key_t_ext_id);

# pragma omp parallel
  {
    cs_lnum_t t_s_id, t_e_id;
    cs_parall_thread_range(mesh->n_cells, sizeof(cs_real_t), &t_s_id, &t_e_id);

    /* Laminar viscosity */

    if (cs_glob_time_scheme->isno2t > 0 && iviext > 0) {
      cs_real_t *cpro_viscl_pre = CS_F_(mu)->val_pre;
      for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++)
        secvis[c_id] = d2s3m * cpro_viscl_pre[c_id];
    }
    else {
      for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++)
        secvis[c_id] = d2s3m * viscl[c_id];
    }

    /* Volume viscosity if present */
    if (cs_glob_physical_model_flag[CS_COMPRESSIBLE] >= 0) {
      if (f_viscv != NULL) {
        for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++)
          secvis[c_id] += + cpro_viscv[c_id];
      }
      else {
        const cs_real_t viscv0 = cs_glob_fluid_properties->viscv0;
        for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++)
          secvis[c_id] += viscv0;
      }
    }

    /* Turbulent viscosity (if not in Rij or LES) */

    if (itytur != 3 && itytur != 4) {

      if (cs_glob_time_scheme->isno2t > 0 && iviext_t > 0) {
        cs_real_t *cpro_visct_pre = CS_F_(mu_t)->val_pre;
        for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++)
          secvis[c_id] += d2s3m * cpro_visct_pre[c_id];
      }
      else {
        for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++)
          secvis[c_id] += d2s3m * visct[c_id];
      }

    }

    /* With porosity */
    if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
      cs_real_t *porosity = CS_F_(poro)->val;
      for (cs_lnum_t c_id = t_s_id; c_id < t_e_id; c_id++)
        secvis[c_id] *= porosity[c_id];
    }

  }

  /* Parallelism and periodicity process */

  if (mesh->halo != NULL)
    cs_halo_sync_var(mesh->halo, CS_HALO_STANDARD, secvis);

  /* Interior faces
     TODO we should (re)test the weight walue for imvisf=0 */

  if (eqp_vel->imvisf == 0) {
#   pragma omp parallel for if (n_i_faces > CS_THR_MIN)
    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

      const cs_lnum_t c_id1 = i_face_cells[f_id][0];
      const cs_lnum_t c_id2 = i_face_cells[f_id][1];

      const cs_real_t secvsi = secvis[c_id1];
      const cs_real_t secvsj = secvis[c_id2];

      secvif[f_id] = 0.5 * (secvsi + secvsj);
    }
  }
  else {
#   pragma omp parallel for if (n_i_faces > CS_THR_MIN)
    for (cs_lnum_t f_id = 0; f_id < n_i_faces; f_id++) {

      const cs_lnum_t c_id1 = i_face_cells[f_id][0];
      const cs_lnum_t c_id2 = i_face_cells[f_id][1];

      const cs_real_t secvsi = secvis[c_id1];
      const cs_real_t secvsj = secvis[c_id2];
      const cs_real_t pnd = weight[f_id];

      secvif[f_id] = secvsi * secvsj / (pnd * secvsi + (1.0 - pnd) * secvsj);
    }
  }

  /* Boundary faces
     TODO shall we extrapolate this value? */

#  pragma omp parallel for if (n_b_faces > CS_THR_MIN)
  for (cs_lnum_t f_id = 0; f_id < n_b_faces; f_id++) {
    const cs_lnum_t c_id = b_face_cells[f_id];
    secvib[f_id] = secvis[c_id];
  }

  /* TODO stresses at the wall? */

  /* Free memory */
  BFT_FREE(secvis);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the diffusion velocity at faces.
 *
 * i_visc,b_visc = viscosity*surface/distance, homogeneous to a rate of flow
 * in kg/s.
 *
 * <a name="viscfa"></a>
 *
 * Please refer to the
 * <a href="../../theory.pdf#viscfa"><b>viscfa</b></a> section of the theory
 * guide for more informations.
 *
 * \remark: a priori, no need of reconstruction techniques
 * (to improve if necessary).
 *
 * \param[in]     m              pointer to mesh
 * \param[in]     fvq            pointer to finite volume quantities
 * \param[in]     visc_mean_type method to compute the viscosity at faces:
 *                                - 0 arithmetical
 *                                - 1 harmonic
 * \param[in]     c_visc         cell viscosity (scalar)
 * \param[out]    i_visc         inner face viscosity
 *                                (times surface divided by distance)
 * \param[out]    b_visc         boundary face viscosity
 *                                (surface, must be consistent with flux BCs)
 */
/*----------------------------------------------------------------------------*/

void
cs_face_viscosity(const cs_mesh_t               *m,
                  const cs_mesh_quantities_t    *fvq,
                  const int                      visc_mean_type,
                  cs_real_t            *restrict c_visc,
                  cs_real_t            *restrict i_visc,
                  cs_real_t            *restrict b_visc)
{
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_f_face_surf = fvq->i_f_face_surf;
  const cs_real_t *restrict b_f_face_surf = fvq->b_f_face_surf;

  /* Porosity field */
  cs_field_t *fporo = cs_field_by_name_try("porosity");

  cs_real_t *porosi = NULL;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
  }

  /* ---> Periodicity and parallelism treatment */
  if (halo != NULL) {
    cs_halo_type_t halo_type = CS_HALO_STANDARD;
    cs_halo_sync_var(halo, halo_type, c_visc);
    if (porosi != NULL) cs_halo_sync_var(halo, halo_type, porosi);
  }

  /* Without porosity */
  if (porosi == NULL) {

    if (visc_mean_type == 0) {

      for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        double visci = c_visc[ii];
        double viscj = c_visc[jj];

        i_visc[face_id] = 0.5*(visci+viscj)
                         *i_f_face_surf[face_id]/i_dist[face_id];

      }

    } else {

      for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        double visci = c_visc[ii];
        double viscj = c_visc[jj];
        double pnd   = weight[face_id];

        i_visc[face_id] = visci*viscj/CS_MAX(pnd*visci+(1.-pnd)*viscj,
                                             DBL_MIN)
                         *i_f_face_surf[face_id]/i_dist[face_id];

      }

    }

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      b_visc[face_id] = b_f_face_surf[face_id];

    }

  /* With porosity */
  } else {

    if (visc_mean_type == 0) {

      for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        double visci = c_visc[ii] * porosi[ii];
        double viscj = c_visc[jj] * porosi[jj];

        i_visc[face_id] = 0.5*(visci+viscj)
                         *i_f_face_surf[face_id]/i_dist[face_id];

      }

    } else {

      for (cs_lnum_t face_id = 0; face_id <m->n_i_faces; face_id++) {

        cs_lnum_t ii = i_face_cells[face_id][0];
        cs_lnum_t jj = i_face_cells[face_id][1];

        double visci = c_visc[ii] * porosi[ii];
        double viscj = c_visc[jj] * porosi[jj];
        double pnd   = weight[face_id];

        i_visc[face_id] = visci*viscj/CS_MAX(pnd*visci+(1.-pnd)*viscj,
                                             DBL_MIN)
                         *i_f_face_surf[face_id]/i_dist[face_id];

      }

    }

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];

      b_visc[face_id] = b_f_face_surf[face_id]*porosi[ii];

    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the equivalent tensor viscosity at faces for a 3x3 symetric
 * tensor.
 *
 * \param[in]     m              pointer to mesh
 * \param[in]     fvq            pointer to finite volume quantities
 * \param[in]     visc_mean_type method to compute the viscosity at faces:
 *                                - 0: arithmetic
 *                                - 1: harmonic
 * \param[in]     c_visc         cell viscosity symmetric tensor
 * \param[out]    i_visc         inner face tensor viscosity
 *                                (times surface divided by distance)
 * \param[out]    b_visc         boundary face viscosity
 *                                (surface, must be consistent with flux BCs)
 */
/*----------------------------------------------------------------------------*/

void
cs_face_anisotropic_viscosity_vector(const cs_mesh_t             *m,
                                     const cs_mesh_quantities_t  *fvq,
                                     const int                    visc_mean_type,
                                     cs_real_6_t        *restrict c_visc,
                                     cs_real_33_t       *restrict i_visc,
                                     cs_real_t          *restrict b_visc)
{
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict i_f_face_surf = fvq->i_f_face_surf;
  const cs_real_t *restrict b_f_face_surf = fvq->b_f_face_surf;

  double visci[3][3], viscj[3][3], s1[6], s2[6];

  cs_real_6_t *c_poro_visc = NULL;
  cs_real_6_t *w2 = NULL;

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = NULL;
  cs_real_6_t *porosf = NULL;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != NULL) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Without porosity */
  if (porosi == NULL) {

    c_poro_visc = c_visc;

  /* With porosity */
  } else if (porosi != NULL && porosf == NULL) {

    BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (int isou = 0; isou < 6; isou++) {
        w2[cell_id][isou] = porosi[cell_id]*c_visc[cell_id][isou];
      }
    }
    c_poro_visc = w2;

  /* With tensorial porosity */
  } else if (porosi != NULL && porosf != NULL) {

    BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      cs_math_sym_33_product(porosf[cell_id],
                             c_visc[cell_id],
                             w2[cell_id]);
    }
    c_poro_visc = w2;

  }

  /* ---> Periodicity and parallelism treatment */
  if (halo != NULL) {
    cs_halo_type_t halo_type = CS_HALO_STANDARD;
    cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)c_poro_visc, 6);
    if (m->n_init_perio > 0)
      cs_halo_perio_sync_var_sym_tens(halo,
                                      halo_type,
                                      (cs_real_t *)c_poro_visc);
  }

  /* Arithmetic mean */
  if (visc_mean_type == 0) {

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      visci[0][0] = c_poro_visc[ii][0];
      visci[1][1] = c_poro_visc[ii][1];
      visci[2][2] = c_poro_visc[ii][2];
      visci[1][0] = c_poro_visc[ii][3];
      visci[0][1] = c_poro_visc[ii][3];
      visci[2][1] = c_poro_visc[ii][4];
      visci[1][2] = c_poro_visc[ii][4];
      visci[2][0] = c_poro_visc[ii][5];
      visci[0][2] = c_poro_visc[ii][5];

      viscj[0][0] = c_poro_visc[jj][0];
      viscj[1][1] = c_poro_visc[jj][1];
      viscj[2][2] = c_poro_visc[jj][2];
      viscj[1][0] = c_poro_visc[jj][3];
      viscj[0][1] = c_poro_visc[jj][3];
      viscj[2][1] = c_poro_visc[jj][4];
      viscj[1][2] = c_poro_visc[jj][4];
      viscj[2][0] = c_poro_visc[jj][5];
      viscj[0][2] = c_poro_visc[jj][5];

      for (int isou = 0; isou < 3; isou++) {
        for (int jsou = 0; jsou < 3; jsou++) {
          i_visc[face_id][jsou][isou] =  0.5*(visci[jsou][isou]
                                             +viscj[jsou][isou])
                                       * i_f_face_surf[face_id]/i_dist[face_id];
        }
      }

    }

    /* Harmonic mean: Kf = Ki . (pnd Ki +(1-pnd) Kj)^-1 . Kj */
  } else {

    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

      cs_lnum_t ii = i_face_cells[face_id][0];
      cs_lnum_t jj = i_face_cells[face_id][1];

      double pnd = weight[face_id];

      for (int isou = 0; isou < 6; isou++) {
        s1[isou] = pnd*c_poro_visc[ii][isou] + (1.-pnd)*c_poro_visc[jj][isou];
      }

      cs_math_sym_33_inv_cramer(s1, s2);

      cs_math_sym_33_product(s2, c_poro_visc[jj], s1);

      cs_math_sym_33_product(c_poro_visc[ii], s1, s2);

      double srfddi = i_f_face_surf[face_id]/i_dist[face_id];

      i_visc[face_id][0][0] = s2[0]*srfddi;
      i_visc[face_id][1][1] = s2[1]*srfddi;
      i_visc[face_id][2][2] = s2[2]*srfddi;
      i_visc[face_id][1][0] = s2[3]*srfddi;
      i_visc[face_id][0][1] = s2[3]*srfddi;
      i_visc[face_id][2][1] = s2[4]*srfddi;
      i_visc[face_id][1][2] = s2[4]*srfddi;
      i_visc[face_id][2][0] = s2[5]*srfddi;
      i_visc[face_id][0][2] = s2[5]*srfddi;

    }

  }

  /* Without porosity */
  if (porosi == NULL) {

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      b_visc[face_id] = b_f_face_surf[face_id];

    }

  /* With porosity */
  } else if (porosi != NULL && porosf == NULL) {

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];

      b_visc[face_id] = b_f_face_surf[face_id]*porosi[ii];

    }

  /* With anisotropic porosity */
  } else if (porosi != NULL && porosf != NULL) {

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];

      b_visc[face_id] = b_f_face_surf[face_id]*porosi[ii];

    }

  }

  BFT_FREE(w2);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the equivalent viscosity at faces for a 3x3 symetric tensor,
 * always using a harmonic mean.
 *
 * \param[in]     m             pointer to mesh
 * \param[in]     fvq           pointer to finite volume quantities
 * \param[in]     c_visc        cell viscosity symmetric tensor
 * \param[in]     iwarnp        verbosity
 * \param[out]    weighf        inner face weight between cells i and j
 *                              \f$ \frac{\vect{IF} \cdot \tens{K}_\celli}
 *                               {\norm{\tens{K}_\celli \cdot \vect{S}}^2} \f$
 *                              and
 *                              \f$ \frac{\vect{FJ} \cdot \tens{K}_\cellj}
 *                               {\norm{\tens{K}_\cellj \cdot \vect{S}}^2} \f$
 * \param[out]    weighb        boundary face weight
 *                              \f$ \frac{\vect{IF} \cdot \tens{K}_\celli}
 *                               {\norm{\tens{K}_\celli \cdot \vect{S}}^2} \f$
 * \param[out]    i_visc        inner face viscosity
 *                               (times surface divided by distance)
 * \param[out]    b_visc        boundary face viscosity
 *                               (surface, must be consistent with flux BCs)
 */
/*----------------------------------------------------------------------------*/

void
cs_face_anisotropic_viscosity_scalar(const cs_mesh_t               *m,
                                     const cs_mesh_quantities_t    *fvq,
                                     cs_real_6_t          *restrict c_visc,
                                     const int                      iwarnp,
                                     cs_real_2_t          *restrict weighf,
                                     cs_real_t            *restrict weighb,
                                     cs_real_t            *restrict i_visc,
                                     cs_real_t            *restrict b_visc)
{
  const cs_halo_t  *halo = m->halo;

  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_lnum_2_t *restrict i_face_cells
    = (const cs_lnum_2_t *restrict)m->i_face_cells;
  const cs_lnum_t *restrict b_face_cells
    = (const cs_lnum_t *restrict)m->b_face_cells;
  const cs_real_t *restrict weight = fvq->weight;
  const cs_real_t *restrict i_dist = fvq->i_dist;
  const cs_real_t *restrict b_dist = fvq->b_dist;
  const cs_real_t *restrict b_f_face_surf = fvq->b_f_face_surf;
  const cs_real_3_t *restrict cell_cen
    = (const cs_real_3_t *restrict)fvq->cell_cen;
  const cs_real_3_t *restrict i_face_normal
    = (const cs_real_3_t *restrict)fvq->i_face_normal;
  const cs_real_t *restrict i_face_surf
    = (const cs_real_t *restrict)fvq->i_face_surf;
  const cs_real_t *restrict i_f_face_surf
    = (const cs_real_t *restrict)fvq->i_f_face_surf;
  const cs_real_3_t *restrict b_face_normal
    = (const cs_real_3_t *restrict)fvq->b_face_normal;
  const cs_real_3_t *restrict i_face_cog
    = (const cs_real_3_t *restrict)fvq->i_face_cog;
  const cs_real_3_t *restrict b_face_cog
    = (const cs_real_3_t *restrict)fvq->b_face_cog;

  cs_gnum_t nclipf = 0, nclipb = 0;
  const double eps = 0.1;

  cs_real_6_t *c_poro_visc = NULL;
  cs_real_6_t *w2 = NULL;

  /* Porosity fields */
  cs_field_t *fporo = cs_field_by_name_try("porosity");
  cs_field_t *ftporo = cs_field_by_name_try("tensorial_porosity");

  cs_real_t *porosi = NULL;
  cs_real_6_t *porosf = NULL;

  if (cs_glob_porous_model == 1 || cs_glob_porous_model == 2) {
    porosi = fporo->val;
    if (ftporo != NULL) {
      porosf = (cs_real_6_t *)ftporo->val;
    }
  }

  /* Without porosity */
  if (porosi == NULL) {

    c_poro_visc = c_visc;

  /* With porosity */
  } else if (porosi != NULL && porosf == NULL) {

    BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      for (int isou = 0; isou < 6; isou++) {
        w2[cell_id][isou] = porosi[cell_id]*c_visc[cell_id][isou];
      }
    }
    c_poro_visc = w2;

  /* With tensorial porosity */
  } else if (porosi != NULL && porosf != NULL) {

    BFT_MALLOC(w2, n_cells_ext, cs_real_6_t);
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
      cs_math_sym_33_product(porosf[cell_id],
                             c_visc[cell_id],
                             w2[cell_id]);
    }
    c_poro_visc = w2;

  }

  /* ---> Periodicity and parallelism treatment */
  if (halo != NULL) {
    cs_halo_type_t halo_type = CS_HALO_STANDARD;
    cs_halo_sync_var_strided(halo, halo_type, (cs_real_t *)c_poro_visc, 6);
    if (m->n_init_perio > 0)
      cs_halo_perio_sync_var_sym_tens(halo,
                                      halo_type,
                                      (cs_real_t *)c_poro_visc);
  }

  /* Always Harmonic mean */
  for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {

    cs_lnum_t ii = i_face_cells[face_id][0];
    cs_lnum_t jj = i_face_cells[face_id][1];

    /* ||Ki.S||^2 */
    cs_real_3_t viscisv;
    cs_math_sym_33_3_product(c_poro_visc[ii], i_face_normal[face_id], viscisv);
    cs_real_t viscis = cs_math_3_square_norm(viscisv);

    /* IF */
    cs_real_3_t fi;
    for (int kk = 0; kk < 3; kk++)
      fi[kk] = i_face_cog[face_id][kk]-cell_cen[ii][kk];

    /* IF.Ki.S */
    cs_real_3_t fiki;
    cs_math_sym_33_3_product(c_poro_visc[ii], fi, fiki);
    cs_real_t fikis = cs_math_3_dot_product(fiki, i_face_normal[face_id]);

    double distfi = (1. - weight[face_id])*i_dist[face_id];

    /* Take I" so that I"F= eps*||FI||*Ki.n when I" is in cell rji */
    double temp = eps*sqrt(viscis)*distfi;
    if (fikis < temp) {
      fikis = temp;
      nclipf++;
    }

    /* ||Kj.S||^2 */
    cs_real_3_t viscjsv;
    cs_math_sym_33_3_product(c_poro_visc[jj], i_face_normal[face_id], viscjsv);
    cs_real_t viscjs = cs_math_3_square_norm(viscjsv);

    /* FJ */
    cs_real_3_t fj;
    for (int kk = 0; kk < 3; kk++)
      fj[kk] = cell_cen[jj][kk]-i_face_cog[face_id][kk];

    /* FJ.Kj.S */
    cs_real_3_t fjkj;
    cs_math_sym_33_3_product(c_poro_visc[jj], fj, fjkj);
    cs_real_t fjkjs = cs_math_3_dot_product(fjkj, i_face_normal[face_id]);

    double distfj = weight[face_id]*i_dist[face_id];

    /* Take J" so that FJ"= eps*||FJ||*Kj.n when J" is in cell i */
    temp = eps*sqrt(viscjs)*distfj;
    if (fjkjs < temp) {
      fjkjs = temp;
      nclipf++;
    }

    weighf[face_id][0] = fikis/viscis;
    weighf[face_id][1] = fjkjs/viscjs;

    i_visc[face_id] = 1./(weighf[face_id][0] + weighf[face_id][1]);

  }

  /* For the porous modelling based on integral formulation Section and fluid
   * Section are different */
  if (cs_glob_porous_model == 3) {
    for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++)
      i_visc[face_id] *= i_f_face_surf[face_id] / i_face_surf[face_id];
  }

  for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

    cs_lnum_t ii = b_face_cells[face_id];

    /* ||Ki.S||^2 */
    cs_real_3_t viscisv;
    cs_math_sym_33_3_product(c_poro_visc[ii], b_face_normal[face_id], viscisv);
    cs_real_t viscis = cs_math_3_square_norm(viscisv);

    /* IF */
    cs_real_3_t fi;
    for (int kk = 0; kk < 3; kk++)
      fi[kk] = b_face_cog[face_id][kk]-cell_cen[ii][kk];

    /* IF.Ki.S */
    cs_real_3_t fiki;
    cs_math_sym_33_3_product(c_poro_visc[ii], fi, fiki);
    cs_real_t fikis = cs_math_3_dot_product(fiki, b_face_normal[face_id]);

    double distfi = b_dist[face_id];

    /* Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji */
    double temp = eps*sqrt(viscis)*distfi;
    if (fikis < temp) {
      fikis = temp;
      nclipb++;
    }

    weighb[face_id] = fikis/viscis;

  }

  /* Without porosity */
  if (porosi == NULL) {

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      /* Warning: hint must be ||Ki.n||/I"F */
      b_visc[face_id] = b_f_face_surf[face_id];
    }

  /* With porosity */
  } else if (porosi != NULL && porosf == NULL) {

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];

      /* Warning: hint must be ||Ki.n||/I"F */
      b_visc[face_id] = b_f_face_surf[face_id]*porosi[ii];

    }

  /* With tensorial porosity */
  } else if (porosi != NULL && porosf != NULL) {

    for (cs_lnum_t face_id = 0; face_id < m->n_b_faces; face_id++) {

      cs_lnum_t ii = b_face_cells[face_id];

      /* Warning: hint must be ||Ki.n||/I"F */
      b_visc[face_id] = b_f_face_surf[face_id]*porosi[ii];

    }

  }

  if (halo != NULL) {
    cs_parall_counter(&nclipf, 1);
    cs_parall_counter(&nclipb, 1);
  }

  if (iwarnp >= 3) {
    bft_printf("Computing the face viscosity from the tensorial viscosity:\n"
               "   Number of internal clippings: %lu\n"
               "   Number of boundary clippings: %lu\n",
               (unsigned long)nclipf, (unsigned long)nclipb);
  }

  BFT_FREE(w2);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
