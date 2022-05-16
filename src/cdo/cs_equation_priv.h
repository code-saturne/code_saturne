#ifndef __CS_EQUATION_PRIV_H__
#define __CS_EQUATION_PRIV_H__

/*============================================================================
 * Functions to handle cs_equation_t structure and its related structures
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_toolbox.h"
#include "cs_equation_param.h"
#include "cs_equation_builder.h"
#include "cs_field.h"
#include "cs_restart.h"

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
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id      id of the variable field
 * \param[in]      bflux_id    id of the boundary flux field
 * \param[in, out] eqb         pointer to a \ref cs_equation_builder_t struct.
 *
 * \return a pointer to a new allocated scheme context structure
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_init_context_t)(const cs_equation_param_t  *eqp,
                             int                         var_id,
                             int                         bflux_id,
                             cs_equation_builder_t      *eqb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a scheme data structure
 *
 * \param[in, out] scheme_context    pointer to a structure cast on-the-fly
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

typedef void *
(cs_equation_free_context_t)(void  *scheme_context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize the variable field values related to an equation
 *
 * \param[in]      t_eval     time at which one performs the evaluation
 * \param[in]      field_id   id related to the variable field of this equation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to the scheme context (cast on-the-fly)
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_init_values_t)(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build and solve a linear system within the CDO framework
 *
 * \param[in]      cur2prev   true="current to previous" operation is performed
 * \param[in]      mesh       pointer to a \ref cs_mesh_t structure
 * \param[in]      field_id   id related to the variable field of this equation
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] eqc        pointer to a scheme context structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_solve_t)(bool                        cur2prev,
                      const cs_mesh_t            *mesh,
                      const int                   field_id,
                      const cs_equation_param_t  *eqp,
                      cs_equation_builder_t      *eqb,
                      void                       *eqc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the Dirichlet boundary stemming from the settings.
 *
 * \param[in]      t_eval      time at which one evaluates BCs
 * \param[in]      mesh        pointer to a cs_mesh_t structure
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in, out] eqb         pointer to a cs_equation_builder_t structure
 * \param[in, out] context     pointer to the scheme context (cast on-the-fly)
 * \param[in, out] field_val   pointer to the values of the variable field
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_set_dir_bc_t)(cs_real_t                     t_eval,
                           const cs_mesh_t              *mesh,
                           const cs_equation_param_t    *eqp,
                           cs_equation_builder_t        *eqb,
                           void                         *context,
                           cs_real_t                     field_val[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system within the CDO framework
 *
 * \param[in]      m          pointer to a \ref cs_mesh_t structure
 * \param[in]      field_val  pointer to the current value of the field
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context    pointer to a scheme context structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_build_system_t)(const cs_mesh_t            *mesh,
                             const cs_real_t            *field_val,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Carry out operations for allocating and/or initializing the solution
 *         array and the right hand side of the linear system to solve.
 *         Handle parallelism thanks to cs_range_set_t structure.
 *
 * \param[in, out] eq_cast    pointer to generic builder structure
 * \param[in, out] p_x        pointer of pointer to the solution array
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_prepare_solve_t)(void              *eq_to_cast,
                              cs_real_t         *p_x[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Store solution(s) of the linear system into a field structure
 *         Update extra-field values if required (for hybrid discretization)
 *
 * \param[in]      solu       solution array
 * \param[in]      rhs        rhs associated to this solution array
 * \param[in]      eqp        pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context    pointer to a context structure
 * \param[in, out] field_val  pointer to the current value of the field
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_update_field_t)(const cs_real_t            *solu,
                             const cs_real_t            *rhs,
                             const cs_equation_param_t  *eqp,
                             cs_equation_builder_t      *eqb,
                             void                       *context,
                             cs_real_t                  *field_val);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain between time t_cur and t_cur + dt_cur
 *
 * \param[in]      eqp             pointer to a \ref cs_equation_param_t
 * \param[in, out] eqb             pointer to a \ref cs_equation_builder_t
 * \param[in, out] context         pointer to a scheme context structure
 *
 * \return a pointer to a cs_cdo_balance_t structure
 */
/*----------------------------------------------------------------------------*/

typedef cs_cdo_balance_t *
(cs_equation_get_balance_t)(const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the cellwise stiffness matrix associated to the property
 *         given as a parameter and apply it to the pot array to define
 *         the resulting array associated to entities defined at loc_res
 *
 * \param[in]      eqp      pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb      pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context  pointer to a scheme builder structure
 * \param[in]      property pointer to the property related to the stiffness op.
 * \param[in]      pot      array to multiply with the stiffness matrix
 * \param[in]      loc_res  location of entities in the resulting array
 * \param[in, out] res      resulting array
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_apply_sitffness_t)(const cs_equation_param_t     *eqp,
                                cs_equation_builder_t         *eqb,
                                void                          *context,
                                const cs_property_t           *property,
                                const cs_real_t               *pot,
                                cs_flag_t                      loc_res,
                                cs_real_t                     *res);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic prototype for extra-operations related to an equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to a generic data structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_extra_op_t)(const cs_equation_param_t  *eqp,
                         cs_equation_builder_t      *eqb,
                         void                       *context);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve cellwise structure including work buffers used to build
 *         a CDO system cellwise. Generic prototype for all CDO schemes.
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_get_builders_t)(cs_cell_sys_t       **csys,
                             cs_cell_builder_t   **cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute or retrieve an array of values at a given mesh location
 *         Currently, vertices, cells or faces are possible locations
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  scheme_context  pointer to a data structure cast on-the-fly
 * \param[in]       previous        retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of \ref cs_real_t
 */
/*----------------------------------------------------------------------------*/

typedef cs_real_t *
(cs_equation_get_values_t)(void      *scheme_context,
                           bool       previous);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Generic prototype dedicated to read or write additional arrays (not
 *         defined as fields) useful for the checkpoint/restart process
 *
 * \param[in, out]  restart         pointer to \ref cs_restart_t structure
 * \param[in]       eqname          name of the related equation
 * \param[in, out]  scheme_context  pointer to a data structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

typedef void
(cs_equation_restart_t)(cs_restart_t    *restart,
                        const char      *eqname,
                        void            *scheme_context);

/*----------------------------------------------------------------------------
 * Structure type
 *----------------------------------------------------------------------------*/

/*! \struct cs_equation_t
 *  \brief Main structure to handle the discretization and the resolution of
 *         an equation.
 */

struct _cs_equation_t {

  /*!
   * @name Metadata
   * @{
   *
   * \var varname
   * Name of the field of type variable associated to this equation
   *
   * \var field_id
   * Id of the variable field related to this equation (cf. \ref cs_field_t)
   *
   * \var boundary_flux_id
   * Id to the field storing the boundary flux asociated to variable field
   */

  int                    id;    /*!< Id of the current equation  */
  char *restrict         varname;
  int                    field_id;
  int                    boundary_flux_id;

  /* Timer statistic for a coarse profiling */

  int     main_ts_id;   /*!< Id of the main timer stats for this equation */

  /*! \var param
   *  Set of parameters related to an equation
   */

  cs_equation_param_t   *param;

  /*!
   * @}
   * @name Generic pointers to manage a cs_equation_t structure
   * @{
   *
   * Since each space discretization has its specificities, an operation is
   * associated to a function pointer which is set during the setup phase
   * according to the space discretization and then one calls this operation
   * without the need to test in which case we are.
   *
   * \var builder
   * Common members for building the algebraic system between the numerical
   * schemes
   */

  cs_equation_builder_t   *builder;

  /*! \var scheme_context
   * Data depending on the numerical scheme (cast on-the-fly)
   */

  void                    *scheme_context;

  /*!
   * \var init_context
   * Pointer of function given by the prototype cs_equation_init_context_t
   *
   * \var free_context
   * Pointer of function given by the prototype cs_equation_free_context_t
   *
   * \var init_field_values
   * Pointer of function given by the prototype cs_equation_init_values_t
   *
   * \var solve_steady_state
   * Case of steady-state solution. Pointer of function given by the prototype
   * cs_equation_solve_t
   *
   * \var solve
   * Pointer to a function given by the generic prototype cs_equation_solve_t
   * Case of an unsteady solution
   *
   * \var compute_balance
   * Pointer of function given by the prototype cs_equation_get_balance_t
   *
   * \var apply_stiffness
   * Pointer of function given by the prototype cs_equation_apply_sitffness_t
   *
   * \var postprocess
   * Additional predefined post-processing. Pointer of function given by the
   * prototype cs_equation_extra_op_t
   *
   * \var current_to_previous
   * Pointer of function given by the prototype cs_equation_extra_op_t
   *
   * \var read_restart
   * Pointer of function given by the prototype cs_equation_restart_t
   *
   * \var write_restart
   * Pointer of function given by the prototype cs_equation_restart_t
   *
   * \var get_cell_values
   * Pointer of function given by the prototype cs_equation_get_values_t
   *
   * \var get_face_values
   * Pointer of function given by the prototype cs_equation_get_values_t
   *
   * \var get_edge_values
   * Pointer of function given by the prototype cs_equation_get_values_t
   *
   * \var get_vertex_values
   * Pointer of function given by the prototype cs_equation_get_values_t
   *
   * \var get_cw_build_structures
   * Retrieve local structures used to build the algebraic system. Pointer of
   * function given by the prototype cs_equation_get_builders_t
   */

  cs_equation_init_context_t       *init_context;
  cs_equation_free_context_t       *free_context;

  cs_equation_init_values_t        *init_field_values;
  cs_equation_solve_t              *solve_steady_state;
  cs_equation_solve_t              *solve;

  cs_equation_get_balance_t        *compute_balance;
  cs_equation_apply_sitffness_t    *apply_stiffness;
  cs_equation_extra_op_t           *postprocess;
  cs_equation_extra_op_t           *current_to_previous;

  cs_equation_restart_t            *read_restart;
  cs_equation_restart_t            *write_restart;

  cs_equation_get_values_t         *get_cell_values;
  cs_equation_get_values_t         *get_face_values;
  cs_equation_get_values_t         *get_edge_values;
  cs_equation_get_values_t         *get_vertex_values;

  cs_equation_get_builders_t       *get_cw_build_structures;

  /* Deprecated functions --> use rather solve() and solve_steady_state()
   * Only HHO schemes rely on these functions
   */

  cs_equation_set_dir_bc_t         *set_dir_bc;
  cs_equation_build_system_t       *build_system;
  cs_equation_prepare_solve_t      *prepare_solving;
  cs_equation_update_field_t       *update_field;

  /*!
   * @}
   */

};

/*============================================================================
 * Public function prototypes
 *============================================================================*/


/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EQUATION_PRIV_H__ */
