#ifndef __CS_CDO_BLAS_H__
#define __CS_CDO_BLAS_H__

/*============================================================================
 * Functions computing BLAS 1 operations (like square norms and dot products
 * from CDO quantities)
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for computing a dot product. Parallel
 *         synchronization is performed.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

typedef cs_real_t
(cs_cdo_blas_dotprod_t) (const cs_real_t     *a,
                         const cs_real_t     *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for computing a square norm. Parallel
 *         synchronization is performed.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

typedef cs_real_t
(cs_cdo_blas_square_norm_t) (const cs_real_t     *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic function pointer for computing a square norm of the
 *         difference between two arrays (defined at the same location and of
 *         the same dimension). The result may be normalized by the norm of the
 *         second array. Parallel synchronization is performed.
 *
 * \param[in]  a     first array
 * \param[in]  b     second array
 *
 * \return the square weighted L2-norm of the difference
 */
/*----------------------------------------------------------------------------*/

typedef cs_real_t
(cs_cdo_blas_square_norm_diff_t) (const cs_real_t     *a,
                                  const cs_real_t     *b);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a cs_cdo_connect_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_blas_init_sharing(const cs_cdo_quantities_t    *quant,
                         const cs_cdo_connect_t       *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a scalar-valued array defined as a potential at primal
 *         cells. Thus, the weigth is the cell volume. The computed quantities
 *         are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pcsp(const cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the norm ||b - a||**2
 *         Case of two scalar-valued arrays a and b defined as a potential at
 *         primal cells. Thus, the weigth is the cell volume. The computed
 *         quantities are synchronized in parallel.
 *
 * \param[in]  a   first array
 * \param[in]  b   second array
 *
 * \return the value of ||b - a||**2
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pcsp_diff(const cs_real_t        *a,
                                  const cs_real_t        *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the norm  ||a - ref||**2 / || ref||**2
 *         Case of two scalar-valued arrays a and ref defined as a potential at
 *         primal cells. Thus, the weigth is the cell volume. The computed
 *         quantities are synchronized in parallel. "ndiff" stands for
 *         "normalized difference"
 *
 * \param[in]  a     array to analyze
 * \param[in]  ref   array used for normalization and difference
 *
 * \return the normalized square weighted L2-norm of the difference between the
 *         two arrays
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pcsp_ndiff(const cs_real_t        *a,
                                   const cs_real_t        *ref);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using the classical Euclidean
 *         dot product (without weight).
 *         Case of a scalar-valued arrays defined at primal vertices.
 *         The computed quantity is synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_dotprod_vertex(const cs_real_t        *a,
                           const cs_real_t        *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array using an Euclidean 2-norm.
 *         Case of a scalar-valued array defined at primal vertices.
 *         The computed quantities are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_square_norm_vertex(const cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using a weighted Euclidean dot
 *         product relying on CDO quantities.
 *         Case of a scalar-valued arrays defined as a potential at primal
 *         vertices. Thus, the weigth is the portion of dual cell (associated
 *         to a primal vertex) inside a primal cell.  The computed quantity is
 *         synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_dotprod_pvsp(const cs_real_t        *a,
                         const cs_real_t        *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a scalar-valued array defined as a potential at primal
 *         vertices. Thus, the weigth is the portion of dual cell inside each
 *         (primal cell). The computed quantities are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pvsp(const cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the norm ||b - a||**2
 *         Case of two scalar-valued arrays a and b defined as a potential at
 *         primal vertices. Thus, the weigth is the portion of dual cell in a
 *         primal cell. The computed quantities are synchronized in parallel.
 *
 * \param[in]  a   first array
 * \param[in]  b   second array
 *
 * \return the value  of ||b - a||**2
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pvsp_diff(const cs_real_t        *a,
                                  const cs_real_t        *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a non-interlaced scalar-valued array of stride = 2 defined as
 *         a potential at primal vertices. Thus, the weigth is the portion of
 *         dual cell (associated to a primal vertex) inside a primal cell. The
 *         computed quantity is synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_square_norm_2pvsp(const cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using a weighted Euclidean dot
 *         product relying on CDO quantities.
 *         Case of non-interlaced scalar-valued arrays of stride = 2 defined as
 *         a potential at primal vertices. Thus, the weigth is the portion of
 *         dual cell (associated to a primal vertex) inside a primal cell. The
 *         computed quantity is synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

double
cs_cdo_blas_dotprod_2pvsp(const cs_real_t        *a,
                          const cs_real_t        *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using the classical Euclidean
 *         dot product (without weight).
 *         Case of a scalar-valued arrays defined at primal faces.
 *         The computed quantity is synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_dotprod_face(const cs_real_t        *a,
                         const cs_real_t        *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array using an Euclidean 2-norm.
 *         Case of a scalar-valued array defined at primal faces.
 *         The computed quantities are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_face(const cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a scalar-valued array defined as a potential at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. The computed quantities are synchronized in
 *         parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pfsp(const cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a vector-valued array defined as a potential at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. The computed quantities are synchronized in
 *         parallel.
 *
 * \param[in]  array   array to analyze (vector-valued)
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pfvp(const cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the dot product of two arrays using a weighted Euclidean
 *         dot product relying on CDO quantities.
 *         Case of a scalar-valued arrays defined as a flux at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. Each face quantity is normalized by the face
 *         surface. The computed quantity is synchronized in parallel.
 *
 * \param[in]  a   first array to analyze
 * \param[in]  b   second array to analyze
 *
 * \return the value of the dot product
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_dotprod_pfsf(const cs_real_t        *a,
                         const cs_real_t        *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the square norm of an array
 *         Case of a scalar-valued array defined as a flux at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. Each face quantity is normalized by the face
 *         surface. The computed quantities are synchronized in parallel.
 *
 * \param[in]  array   array to analyze
 *
 * \return the square weighted L2-norm
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pfsf(const cs_real_t        *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the norm ||b - a||**2
 *         Case of a scalar-valued array defined as a flux at primal
 *         faces. Thus, the weigth is the pyramid of apex the cell center and
 *         of basis the face. Each face quantity is normalized by the face
 *         surface. The computed quantities are synchronized in parallel.
 *
 * \param[in]  a   first array
 * \param[in]  b   second array
 *
 * \return the value of ||b - a||**2
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_cdo_blas_square_norm_pfsf_diff(const cs_real_t        *a,
                                  const cs_real_t        *b);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_BLAS_H__ */
