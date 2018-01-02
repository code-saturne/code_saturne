#ifndef __CS_SOURCE_TERM_H__
#define __CS_SOURCE_TERM_H__

/*============================================================================
 * Functions and structures to deal with source term computation
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_N_MAX_SOURCE_TERMS   8 // Max number of source terms in an equation

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  char *restrict name;      // short description of the source term
  int            ml_id;     // id of the related mesh location structure

  /* Specification related to the way of computing the source term */
  cs_flag_t                flag;       // Metadata
  cs_quadra_type_t         quad_type;  // barycentric, higher, highest
  cs_param_def_type_t      def_type;   // by value, by function...
  cs_def_t                 def;

  /* Useful buffers to deal with more complex definitions */
  cs_desc_t      array_desc;  // short description of the related array
  cs_real_t     *array;       // if the source term hinges on an array
  const void    *struc;       // if one needs a structure. Only shared

} cs_source_term_t;

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_source_term_cellwise_t)(const cs_source_term_t    *source,
                            const cs_cell_mesh_t      *cm,
                            cs_cell_builder_t         *cb,
                            double                    *values);


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
 * \brief  Destroy an array of cs_source_term_t structures
 *
 * \param[in]      n_source_terms   number of source terms
 * \param[in, out] source_terms     pointer to a cs_source_term_t structure
 *
 * \return NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_source_term_t *
cs_source_term_destroy(int                 n_source_terms,
                       cs_source_term_t   *source_terms);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic way to define a cs_source_term_t structure by value
 *
 * \param[in, out] st        pointer to the cs_source_term_t structure to set
 * \param[in]      st_id     id related to source term to define
 * \param[in]      name      name of the source term
 * \param[in]      var_type  type of variable (scalar, vector...)
 * \param[in]      ml_id     id related to the mesh location
 * \param[in]      flag      metadata related to this source term
 * \param[in]      val       accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_def_by_value(cs_source_term_t    *st,
                            int                  st_id,
                            const char          *name,
                            cs_param_var_type_t  var_type,
                            int                  ml_id,
                            cs_flag_t            flag,
                            const char          *val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_source_term_t structure thanks to an analytic function
 *
 * \param[in, out] st        pointer to the cs_source_term_t structure to set
 * \param[in]      st_id     id related to source term to define
 * \param[in]      name      name of the source term
 * \param[in]      var_type  type of variable (scalar, vector...)
 * \param[in]      ml_id     id related to the mesh location
 * \param[in]      flag      metadata related to this source term
 * \param[in]      func      pointer to a function
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_def_by_analytic(cs_source_term_t     *st,
                               int                   st_id,
                               const char           *name,
                               cs_param_var_type_t   var_type,
                               int                   ml_id,
                               cs_flag_t             flag,
                               cs_analytic_func_t   *func);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_source_term_t structure thanks to an array of values
 *
 * \param[in, out] st        pointer to the cs_source_term_t structure to set
 * \param[in]      st_id     id related to source term to define
 * \param[in]      name      name of the source term
 * \param[in]      var_type  type of variable (scalar, vector...)
 * \param[in]      ml_id     id related to the mesh location
 * \param[in]      flag      metadata related to this source term
 * \param[in]      desc      description of the main feature of this array
 * \param[in]      array     pointer to an array
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_def_by_array(cs_source_term_t     *st,
                            int                   st_id,
                            const char           *name,
                            cs_param_var_type_t   var_type,
                            int                   ml_id,
                            cs_flag_t             flag,
                            cs_desc_t             desc,
                            cs_real_t            *array);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the type of quadrature to use for computing the source term
 *
 * \param[in, out]  st          pointer to a cs_source_term_t structure
 * \param[in]       quad_type   type of quadrature to use
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_quadrature(cs_source_term_t  *st,
                              cs_quadra_type_t   quad_type);

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
 * \param[in, out]  st        pointer to a cs_source_term_t structure
 * \param[in]       flag      CS_FLAG_DUAL or CS_FLAG_PRIMAL
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_set_reduction(cs_source_term_t     *st,
                             cs_flag_t             flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get metadata related to the given source term structure
 *
 * \param[in, out]  st          pointer to a cs_source_term_t structure
 *
 * \return the value of the flag related to this source term
 */
/*----------------------------------------------------------------------------*/

cs_flag_t
cs_source_term_get_flag(const cs_source_term_t  *st);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the name related to a cs_source_term_t structure
 *
 * \param[in] st      pointer to a cs_source_term_t structure
 *
 * \return the name of the source term
 */
/*----------------------------------------------------------------------------*/

const char *
cs_source_term_get_name(const cs_source_term_t   *st);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summarize the content of a cs_source_term_t structure
 *
 * \param[in] eqname  name of the related equation
 * \param[in] st      pointer to a cs_source_term_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_summary(const char               *eqname,
                       const cs_source_term_t   *st);

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
                    const cs_source_term_t      *source_terms,
                    cs_source_term_cellwise_t   *compute_source[],
                    cs_flag_t                   *sys_flag,
                    cs_mask_t                   *source_mask[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution related to a source term
 *
 * \param[in]      dof_desc   description of the associated DoF
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in, out] p_values   pointer to the computed values (allocated if NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute(cs_desc_t                     dof_desc,
                       const cs_source_term_t       *source,
                       double                       *p_values[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute the local contributions of source terms in a cell
 *
 * \param[in]      n_source_terms  number of source terms
 * \param[in]      source_terms    pointer to the definitions of source terms
 * \param[in]      cm              pointer to a cs_cell_mesh_t structure
 * \param[in]      sys_flag        metadata about the algebraic system
 * \param[in]      source_mask     array storing in a compact way which source
 *                                 term is defined in a given cell
 * \param[in]      compute_source  array of function pointers
 * \param[in, out] cb              pointer to a cs_cell_builder_t structure
 * \param[in, out] csys            cellwise algebraic system
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_compute_cellwise(const int                    n_source_terms,
                                const cs_source_term_t      *source_terms,
                                const cs_cell_mesh_t        *cm,
                                const cs_flag_t              sys_flag,
                                const cs_mask_t             *source_mask,
                                cs_source_term_cellwise_t   *compute_source[],
                                cs_cell_builder_t           *cb,
                                cs_cell_sys_t               *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar potential defined at primal vertices by a constant
 *         value.
 *         A discrete Hodge operator has to be computed before this call and
 *         stored inside a cs_cell_builder_t structure
 *
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_value(const cs_source_term_t    *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
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
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_pvsp_by_analytic(const cs_source_term_t    *source,
                                const cs_cell_mesh_t      *cm,
                                cs_cell_builder_t         *cb,
                                double                    *values);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the contribution for a cell related to a source term and
 *         add it the given array of values.
 *         Case of a scalar density defined at dual cells by a value.
 *
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_by_value(const cs_source_term_t    *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
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
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_bary_by_analytic(const cs_source_term_t    *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
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
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q1o1_by_analytic(const cs_source_term_t    *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
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
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q10o2_by_analytic(const cs_source_term_t    *source,
                                      const cs_cell_mesh_t      *cm,
                                      cs_cell_builder_t         *cb,
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
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_dcsd_q5o3_by_analytic(const cs_source_term_t    *source,
                                     const cs_cell_mesh_t      *cm,
                                     cs_cell_builder_t         *cb,
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
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_vcsp_by_value(const cs_source_term_t    *source,
                             const cs_cell_mesh_t      *cm,
                             cs_cell_builder_t         *cb,
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
 * \param[in]      source     pointer to a cs_source_term_t structure
 * \param[in]      cm         pointer to a cs_cell_mesh_t structure
 * \param[in, out] cb         pointer to a cs_cell_builder_t structure
 * \param[in, out] values     pointer to the computed values
 */
/*----------------------------------------------------------------------------*/

void
cs_source_term_vcsp_by_analytic(const cs_source_term_t    *source,
                                const cs_cell_mesh_t      *cm,
                                cs_cell_builder_t         *cb,
                                double                    *values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SOURCE_TERM_H__ */
