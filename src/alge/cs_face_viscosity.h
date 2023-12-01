#ifndef __CS_FACE_VISCOSITY_H__
#define __CS_FACE_VISCOSITY_H__

/*============================================================================
 * Face viscosity
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

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
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_viscosity
 *----------------------------------------------------------------------------*/

void CS_PROCF (viscfa, VISCFA)
(
 const int   *const visc_mean_type,
 cs_real_t          c_visc[],
 cs_real_t          i_visc[],
 cs_real_t          b_visc[]
);

/*----------------------------------------------------------------------------
 * Wrapper to cs_face_anisotropic_viscosity_vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (vistnv, VISTNV)
(
 const int     *const visc_mean_type,
 cs_real_6_t          c_visc[],
 cs_real_33_t         i_visc[],
 cs_real_t            b_visc[]
);

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
);

/*============================================================================
 * Public function prototypes
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
                            cs_real_t  secvib[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute the diffusion velocity at faces.
 * i_visc,b_visc = viscosity*surface/distance, homogeneous to a rate of flow
 * in kg/s.
 *
 * Remark: a priori, no need of reconstruction techniques
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
                  cs_real_t            *restrict b_visc);

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
                                     cs_real_t          *restrict b_visc);

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
                                     cs_real_t            *restrict b_visc);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FACE_VISCOSITY_H__ */
