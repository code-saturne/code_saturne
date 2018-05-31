#ifndef __CS_EQUATION_PRIV_H__
#define __CS_EQUATION_PRIV_H__

/*============================================================================
 * Routines to handle cs_equation_t structure and its related structures
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_equation_param.h"
#include "cs_equation_common.h"
#include "cs_field.h"
#include "cs_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a scheme data structure used during the building of the
 *         algebraic system
 *
 * \param[in]      eqp    pointer to a cs_equation_param_t structure
 * \param[in, out] eqb    pointer to a cs_equation_builder_t structure
 *
 * \return a pointer to a new allocated scheme builder structure
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_init_context_t)(const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a scheme data structure
 *
 * \param[in, out]  scheme_context   pointer to a builder structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_free_context_t)(void  *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create the matrix of the current algebraic system.
 *         Allocate and initialize the right-hand side associated to the given
 *         builder structure
 *
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to generic data structure
 * \param[in, out] system_matrix  pointer of pointer to a cs_matrix_t struct.
 * \param[in, out] system_rhs     pointer of pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_initialize_system_t)(const cs_equation_param_t   *eqp,
                                  cs_equation_builder_t       *eqb,
                                  void                        *data,
                                  cs_matrix_t                **system_matrix,
                                  cs_real_t                  **system_rhs);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the Dirichlet boundary stemming from the settings.
 *
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in, out] eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      t_eval      time at which one evaluates BCs
 * \param[in, out] field_val   pointer to the values of the variable field
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_set_dir_bc_t)(const cs_mesh_t              *mesh,
                           const cs_equation_param_t    *eqp,
                           cs_equation_builder_t        *eqb,
                           cs_real_t                     t_eval,
                           cs_real_t                     field_val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system within the CDO framework
 *
 * \param[in]      m          pointer to a cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      dt_cur     current value of the time step
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data      pointer to a scheme builder structure
 * \param[in, out] rhs        right-hand side to compute
 * \param[in, out] matrix     pointer to cs_matrix_t structure to compute
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_build_system_t)(const cs_mesh_t            *mesh,
                             const cs_real_t            *field_val,
                             double                      dt_cur,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data,
                             cs_real_t                  *rhs,
                             cs_matrix_t                *matrix);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Carry out operations for allocating and/or initializing the solution
 *         array and the right hand side of the linear system to solve.
 *         Handle parallelism thanks to cs_range_set_t structure.
 *
 * \param[in, out] eq_cast    pointer to generic builder structure
 * \param[in, out] p_x        pointer of pointer to the solution array
 * \param[in, out] p_rhs      pointer of pointer to the RHS array
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_prepare_solve_t)(void              *eq_to_cast,
                              cs_real_t         *p_x[],
                              cs_real_t         *p_rhs[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] data       pointer to a data structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_update_field_t)(const cs_real_t            *solu,
                             const cs_real_t            *rhs,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *data,
                             cs_real_t                  *field_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain between time t_cur and t_cur + dt_cur
 *
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in, out] eqb             pointer to a cs_equation_builder_t structure
 * \param[in, out] context         pointer to a scheme builder structure
 * \param[in]      var_field_id    id of the variable field
 * \param[in]      bflux_field_id  id of the variable field
 * \param[in]      dt_cur          current value of the time step
 *
 * \return a pointer to a cs_equation_balance_t structure
 */
/*----------------------------------------------------------------------------*/

typedef cs_equation_balance_t *
(cs_equation_get_balance_t)(const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context,
                            int                           var_field_id,
                            int                           bflux_field_id,
                            cs_real_t                     dt_cur);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the diffusive and convective flux across a list of faces
 *
 * \param[in]       normal     indicate in which direction flux is > 0
 * \param[in]       pdi        pointer to an array of field values
 * \param[in]       ml_id      id related to a cs_mesh_location_t struct.
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to data specific for this scheme
 * \param[in, out]  d_flux     pointer to the value of the diffusive flux
 * \param[in, out]  c_flux     pointer to the value of the convective flux
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_flux_plane_t)(const cs_real_t             normal[],
                           const cs_real_t            *pdi,
                           int                         ml_id,
                           const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *data,
                           double                     *d_flux,
                           double                     *c_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Cellwise computation of the diffusive flux across all faces.
 *         Primal or dual faces are considered according to the space scheme
 *
 * \param[in]       fvals       pointer to an array of field values
 * \param[in]       eqp         pointer to a cs_equation_param_t structure
 * \param[in]       t_eval      time at which one performs the evaluation
 * \param[in, out]  eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out]  data        pointer to a generic data structure
 * \param[in, out]  location    where the flux is defined
 * \param[in, out]  diff_flux   pointer to the value of the diffusive flux
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_cell_difflux_t)(const cs_real_t            *fvals,
                             const cs_equation_param_t  *eqp,
                             cs_real_t                   t_eval,
                             cs_equation_builder_t      *eqb,
                             void                       *data,
                             cs_flag_t                   location,
                             cs_real_t                  *d_flux);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Extra-operation related to this equation
 *
 * \param[in]       eqname     name of the equation
 * \param[in]       field      pointer to a field structure
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  data       pointer to a generic data structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_extra_op_t)(const char                 *eqname,
                         const cs_field_t           *field,
                         const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *data);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the computed values at a different location than that of the
 *         field associated to this equation
 *
 * \param[in]  scheme_context  pointer to a data structure cast on-the-fly
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

typedef double *
(cs_equation_get_extra_values_t)(const void      *scheme_context);

/*----------------------------------------------------------------------------
 * Structure type
 *----------------------------------------------------------------------------*/

/*! \struct cs_equation_t
 *  \brief Main structure to handle the discretization and the resolution of
 *  an equation.
 *
 */

struct _cs_equation_t {

  char *restrict         name;    /* Short description */
  int                    id;

  cs_equation_param_t   *param;   /* Set of parameters related to an equation */

  /* Variable attached to this equation is defined as a cs_field_t structure */
  char *restrict         varname;
  int                    field_id;
  int                    boundary_flux_id;

  /* Algebraic system */
  /* ---------------- */

  /* There are possibly two different sizes for the linear system to handle
     - One for "scatter"-type operations based on the number of geometrical
       entities owned by the local instance of the mesh
     - One for "gather"-type operations based on a balance of the number of
     DoFs from a algebraic point of view. In parallel runs, these two sizes
     can be different.
     n_sles_gather_elts <= n_sles_scatter_elts
  */

  cs_lnum_t                n_sles_scatter_elts;
  cs_lnum_t                n_sles_gather_elts;

  /* Right-hand side defined by a local cellwise building. This may be
     different from the rhs given to cs_sles_solve() in parallel mode. */
  cs_real_t               *rhs;

  /* Matrix to inverse with cs_sles_solve() The matrix size can be different
     from the rhs size in parallel mode since the decomposition is different */
  cs_matrix_t             *matrix;

  /* Range set to handle parallelism. Shared with cs_cdo_connect_t struct.*/
  const cs_range_set_t    *rset;

  /* \var builder
   * Common members for building the algebraic system between the numerical
   * schemes
   */
  cs_equation_builder_t   *builder;

  /* Data depending on the numerical scheme (cast on-the-fly) */
  void                    *scheme_context;

  /* Pointer to functions (see prototypes just above) */
  cs_equation_init_context_t       *init_context;
  cs_equation_free_context_t       *free_context;
  cs_equation_initialize_system_t  *initialize_system;
  cs_equation_set_dir_bc_t         *set_dir_bc;
  cs_equation_build_system_t       *build_system;
  cs_equation_prepare_solve_t      *prepare_solving;
  cs_equation_update_field_t       *update_field;

  cs_equation_get_balance_t        *compute_balance;
  cs_equation_flux_plane_t         *compute_flux_across_plane;
  cs_equation_cell_difflux_t       *compute_cellwise_diff_flux;
  cs_equation_extra_op_t           *postprocess;

  cs_equation_get_extra_values_t   *get_extra_values;

  /* Timer statistic for a "light" profiling */
  int     main_ts_id;   /* Id of the main timer states structure related
                           to this equation */
  int     solve_ts_id;  /* Id of the timer stats structure related to the
                           inversion of the linear system */

  bool    do_build;     /* false => keep the system as it is */

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_PRIV_H__ */
