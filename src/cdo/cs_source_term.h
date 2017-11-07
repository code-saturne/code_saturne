#ifndef __CS_SOURCE_TERM_H__
#define __CS_SOURCE_TERM_H__

/*============================================================================
 * Functions and structures to deal with source term computation
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_cdo.h"
#include "cs_cdo_quantities.h"
#include "cs_cdo_local.h"
#include "cs_param.h"
#include "cs_quadrature.h"
#include "cs_time_step.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_N_MAX_SOURCE_TERMS   8 // Max number of source terms in an equation

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *
 * \param[in]      source     pointer to a cs_xdef_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_source_term_cellwise_t)(const cs_xdef_t         *source,
                            const cs_cell_mesh_t    *cm,
                            cs_cell_builder_t       *cb,
                            void                    *input,
                            double                  *values);


/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers to main domain members
 *
 * \param[in]      quant      additional mesh quantities struct.
 * \param[in]      connect    pointer to a cs_cdo_connect_t struct.
 * \param[in]      time_step  pointer to a time step structure
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_shared_pointers(const cs_cdo_quantities_t    *quant,
                                   const cs_cdo_connect_t       *connect,
                                   const cs_time_step_t         *time_step);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the default flag related to a source term according to the
 *         numerical scheme chosen for discretizing an equation
 *
 * \param[in]       scheme    numerical scheme used for the discretization
 *
 * \return a default flag
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_source_term_set_default_flag(cs_space_scheme_t  scheme);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set advanced parameters which are defined by default in a
 *         source term structure.
 *
 * \param[in, out]  st        pointer to a cs_xdef_t structure
 * \param[in]       flag      CS_FLAG_DUAL or CS_FLAG_PRIMAL
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_reduction(cs_xdef_t     *st,
                             cs_flag_t      flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get metadata related to the given source term structure
 *
 * \param[in, out]  st          pointer to a cs_xdef_t structure
 *
 * \return the value of the flag related to this source term
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_source_term_get_flag(const cs_xdef_t  *st);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize data to build the source terms
 *
 * \param[in]      space_scheme    scheme used to discretize in space
 * \param[in]      n_source_terms  number of source terms
 * \param[in]      source_terms    pointer to the definitions of source terms
 * \param[in, out] compute_source  array of function pointers
 * \param[in, out] sys_flag        metadata about the algebraic system
 * \param[in, out] source_mask     pointer to an array storing in a compact way
 *                                 which source term is defined in a given cell
 *
 * \return a flag which indicates what to build in a cell mesh structure
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_source_term_init(cs_space_scheme_t            space_scheme,
                    const int                    n_source_terms,
                    const cs_xdef_t            **source_terms,
                    cs_source_term_cellwise_t   *compute_source[],
                    cs_flag_t                   *sys_flag,
                    cs_mask_t                   *source_mask[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the local contributions of source terms in a cell
 *
 * \param[in]      n_source_terms  number of source terms
 * \param[in]      source_terms    pointer to the definitions of source terms
 * \param[in]      cm              pointer to a cs_cell_mesh_t structure
 * \param[in]      source_mask     array storing in a compact way which source
 *                                 term is defined in a given cell
 * \param[in]      compute_source  array of function pointers
 * \param[in, out] input           pointer to an element cast on-the-fly
 * \param[in, out] cb              pointer to a cs_cell_builder_t structure
 * \param[in, out] csys            cellwise algebraic system
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute_cellwise(const int                    n_source_terms,
                                const cs_xdef_t            **source_terms,
                                const cs_cell_mesh_t        *cm,
                                const cs_mask_t             *source_mask,
                                cs_source_term_cellwise_t   *compute_source[],
                                void                        *input,
                                cs_cell_builder_t           *cb,
                                cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a source term in the case of
 *         an input data which is a density
 *
 * \param[in]      loc        describe where is located the associated DoF
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in, out] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute_from_density(cs_flag_t                loc,
                                    const cs_xdef_t         *source,
                                    double                  *p_values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a source term in the case of
 *         an input data which is a potential
 *
 * \param[in]      loc        describe where is located the associated DoF
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in, out] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute_from_potential(cs_flag_t                loc,
                                      const cs_xdef_t         *source,
                                      double                  *p_values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar potential defined at primal vertices by a constant
 *         value.
 *         A discrete Hodge operator has to be computed before this call and
 *         stored inside a cs_cell_builder_t structure
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar potential defined at primal vertices by an
 *         analytical function.
 *         A discrete Hodge operator has to be computed before this call and
 *         stored inside a cs_cell_builder_t structure
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_analytic(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by a value.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_bary_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use a the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q1o1_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use a ten-point quadrature rule to evaluate the integral.
 *         Exact for quadratic function.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q10o2_by_analytic(const cs_xdef_t           *source,
                                      const cs_cell_mesh_t      *cm,
                                      cs_cell_builder_t         *cb,
                                      void                      *input,
                                      double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by an analytical
 *         function.
 *         Use a five-point quadrature rule to evaluate the integral.
 *         Exact for cubic function.
 *         This function may be expensive since many evaluations are needed.
 *         Please use it with care.
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q5o3_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar potential defined at primal vertices and cells
 *         by a constant value.
 *         A discrete Hodge operator has to be computed before this call and
 *         stored inside a cs_cell_builder_t structure
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_vcsp_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar potential defined at primal vertices and cells by
 *         an analytical function.
 *         A discrete Hodge operator has to be computed before this call and
 *         stored inside a cs_cell_builder_t structure
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_vcsp_by_analytic(const cs_xdef_t           *source,
                                const cs_cell_mesh_t      *cm,
                                cs_cell_builder_t         *cb,
                                void                      *input,
                                double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density (sd) defined on primal cells by a value.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fbsd_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a vector-valued density (vd) defined on primal cells
 *         by a value.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fbvd_by_value(const cs_xdef_t           *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
                             void                      *input,
                             double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at primal cells by an analytical
 *         function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fbsd_bary_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a vector-valued density defined at primal cells by an
 *         analytical function.
 *         Use the barycentric approximation as quadrature to evaluate the
 *         integral. Exact for linear function.
 *         Case of face-based schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_fbvd_bary_by_analytic(const cs_xdef_t           *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
                                     void                      *input,
                                     double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of a source term for a cell and add it to
 *         the given array of values.
 *         Case of a scalar density (sd) defined on primal cells by a value.
 *         Case of HHO schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_hhosd_by_value(const cs_xdef_t           *source,
                              const cs_cell_mesh_t      *cm,
                              cs_cell_builder_t         *cb,
                              void                      *input,
                              double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution of a source term for a cell and add it to
 *         the given array of values.
 *         Case of a scalar density (sd) defined on primal cells by an analytic
 *         function.
 *         Case of HHO schemes
 *
 * \param[in]      source     pointer to a cs_xdef_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] input      pointer to an element cast on-the-fly (or NULL)
 * \param[in, out] values     pointer to the computed value
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_hhosd_by_analytic(const cs_xdef_t           *source,
                                 const cs_cell_mesh_t      *cm,
                                 cs_cell_builder_t         *cb,
                                 void                      *input,
                                 double                    *values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOURCE_TERM_H__ */
