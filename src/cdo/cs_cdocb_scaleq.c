/*============================================================================
 * Build an algebraic CDO cell-based system for the diffusion equation
 * and solved it as one block (monolithic approach of the flux-potential
 * coupling)
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
#include <assert.h>
#include <string.h>

#if defined(HAVE_OPENMP)
#include <omp.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>

#include "cs_blas.h"
#include "cs_domain.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_blas.h"
#include "cs_cdocb_monolithic_sles.h"
#if defined(DEBUG) && !defined(NDEBUG)
#include "cs_dbg.h"
#endif
#include "cs_equation_bc.h"
#include "cs_equation_builder.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_field_operator.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_sdm.h"
#include "cs_solid_selection.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_timer.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cdocb_priv.h"
#include "cs_cdocb_scaleq.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
 * \file cs_cdocb_scaleq.c
 *
 * \brief Build an algebraic CDO cell-based system for the diffusion
 *        equation and solved it with a monolithic approach
 */

/*=============================================================================
 * Local Macro definitions and structure definitions
 *============================================================================*/

#define CS_CDOCB_SCALEQ_DBG      0

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private variables
 *============================================================================*/

/* Pointer to shared structures */

static const cs_cdo_quantities_t    *cs_shared_quant;
static const cs_cdo_connect_t       *cs_shared_connect;
static const cs_time_step_t         *cs_shared_time_step;

static cs_cell_sys_t               **_scb_cell_system = NULL;
static cs_cell_builder_t           **_scb_cell_builder = NULL;

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local builder structure used for building the system
 *          cellwise
 *
 * \param[in]      connect     pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static cs_cell_builder_t *
_cell_builder_create(const cs_cdo_connect_t   *connect)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define several structures such as the cs_range_set_t,
 *        cs_interface_set_t, cs_matrix_assembler_t and cs_matrix_structure_t
 *        in case of a full assembly (i.e. the full saddle-point matrix is
 *        built).
 *        A variant, activated with add_pressure_diag, is available in order to
 *        enforce the pressure.
 *
 * \param[in, out]  block              pointer to a block structure
 * \param[in]       add_pressure_diag  true or false (pressure diagonal block)
 */
/*----------------------------------------------------------------------------*/

static void
_build_shared_full_structures(cs_cdo_system_block_t     *block,
                              bool                       add_pressure_diag)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Apply the boundary conditions to the cell system when this should
 *          be done before the static condensation.
 *          Case of CDO-Cb schemes with a monolithic flux-potential coupling
 *
 * \param[in]       eqp       pointer to a \ref cs_equation_param_t struct.
 * \param[in]       eqc              context structure for a scalar-valued Cb
 * \param[in]       cm        pointer to a \ref cs_cell_mesh_t structure
 * \param[in]       fm        pointer to a cs_face_mesh_t structure
 * \param[in]       hodge     pointer to a \ref cs_hodge_t structure
 * \param[in, out]  csys      structure storing the cell-wise system
 * \param[in, out]  cb        pointer to a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_scb_apply_bc(const cs_equation_param_t     *eqp,
              const cs_cdocb_scaleq_t       *eqc,
              const cs_cell_mesh_t          *cm,
              cs_face_mesh_t                *fm,
              cs_hodge_t                    *hodge,
              cs_cell_sys_t                 *csys,
              cs_cell_builder_t             *cb)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Perform the assembly stage for a scalar-valued system obtained
 *         with CDO-Cb schemes
 *
 * \param[in]      csys             pointer to a cs_cell_sys_t structure
 * \param[in]      cm               pointer to a cs_cell_mesh_t structure
 * \param[in]      eqb              pointer to a cs_equation_builder_t structure
 * \param[in, out] eqc              context structure for a scalar-valued Cb
 * \param[in, out] asb              pointer to cs_cdo_assembly_t
 */
/*----------------------------------------------------------------------------*/

static void
_full_assembly(const cs_cell_sys_t            *csys,
               const cs_cell_mesh_t           *cm,
               cs_equation_builder_t          *eqb,
               cs_cdocb_scaleq_t              *eqc,
               cs_cdo_assembly_t              *asb)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Initialize the local structure for the current cell
 *
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         pointer to a cs_cdocb_scaleq_t strucutre
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

static void
_cdocb_scaleq_init_cell_system(const cs_cell_mesh_t         *cm,
                               const cs_equation_param_t    *eqp,
                               const cs_equation_builder_t  *eqb,
                               const cs_cdocb_scaleq_t      *eqc,
                               cs_cell_sys_t                *csys,
                               cs_cell_builder_t            *cb)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Build a linear system in the steady-state case
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in]      eqb       pointer to a cs_equation_builder_t structure
 * \param[in,out]  eqc       pointer to a cs_cdocb_scaleq_t strucutre
 */
/*----------------------------------------------------------------------------*/

static void
_steady_build(const cs_equation_param_t      *eqp,
              const cs_equation_builder_t    *eqb,
              cs_cdocb_scaleq_t              *eqc)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current content of the related variables-fields
 *         to previous values
 *
 * \param[in, out] eqc       pointer to a context structure
 */
/*----------------------------------------------------------------------------*/

static void
_cdocb_scaleq_current_to_previous(cs_cdocb_scaleq_t       *eqc)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the divergence in a cell of a flux array defined at
 *         faces (values are defined both at interior and border faces).
 *         Variant based on the usage of \ref cs_cdo_quantities_t structure.
 *
 * \param[in]     cm           pointer to a cellwise view of the mesh
 * \param[in]     f_vals       values of the face DoFs
 *
 * \return the divergence for the corresponding cell
 */
/*----------------------------------------------------------------------------*/

static double
_cdocb_cell_divergence(const cs_cell_mesh_t         *cm,
                       const cs_real_t              *f_vals)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Update the variables related to CDO-Cb system after a resolution
 *
 * \param[in, out] tce       pointer to a timer counter
 * \param[in, out] fld       pointer to a cs_field_t structure
 * \param[in, out] eqc       pointer to a context structure
 * \param[in]      cur2prev  true if one performs "current to previous" op.
 */
/*----------------------------------------------------------------------------*/

static inline void
_update_cell_fields(cs_timer_counter_t      *tce,
                    cs_field_t              *fld,
                    cs_cdocb_scaleq_t       *eqc,
                    bool                     cur2prev)
{
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve work buffers used for building a CDO system cellwise
 *
 * \param[out]  csys   pointer to a pointer on a cs_cell_sys_t structure
 * \param[out]  cb     pointer to a pointer on a cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_get(cs_cell_sys_t       **csys,
                    cs_cell_builder_t   **cb)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Copy current content of the related variables-fields
 *         to previous values
 *         Case of the monolithic coupling algorithm.
 *
 * \param[in]      eqp       pointer to a cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a cs_equation_builder_t structure
 * \param[in, out] context   pointer to a scheme context structure
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_current_to_previous(const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t       *eqb,
                                   void                        *context)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Build the local matrices arising from the diffusion term in the
 *          scalar-valued CDO-Cb schemes.
 *
 * \param[in]      eqp         pointer to a cs_equation_param_t structure
 * \param[in]      eqb         pointer to a cs_equation_builder_t structure
 * \param[in]      eqc         context for this kind of discretization
 * \param[in]      cm          pointer to a cellwise view of the mesh
 * \param[in, out] mass_hodge  pointer to a cs_hodge_t structure for diffusion
 * \param[in, out] csys        pointer to a cellwise view of the system
 * \param[in, out] cb          pointer to a cellwise builder
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_diffusion(const cs_equation_param_t     *eqp,
                          const cs_equation_builder_t   *eqb,
                          const cs_cdocb_scaleq_t       *eqc,
                          const cs_cell_mesh_t          *cm,
                          cs_hodge_t                    *mass_hodge,
                          cs_cell_sys_t                 *csys,
                          cs_cell_builder_t             *cb)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predefined extra-operations related to this equation
 *
 * \param[in]       eqp        pointer to a cs_equation_param_t structure
 * \param[in, out]  eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out]  context    pointer to cs_cdocb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_extra_post(const cs_equation_param_t  *eqp,
                           cs_equation_builder_t      *eqb,
                           void                       *context)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh faces for the current context.
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size n_faces)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdocb_scaleq_get_face_values(void    *context,
                                bool     previous)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve an array of values at mesh vertices for the variable field
 *         associated to the given context
 *         The lifecycle of this array is managed by the code. So one does not
 *         have to free the return pointer.
 *
 * \param[in, out]  context    pointer to a data structure cast on-the-fly
 * \param[in]       previous   retrieve the previous state (true/false)
 *
 * \return  a pointer to an array of cs_real_t (size: n_cells)
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdocb_scaleq_get_cell_values(void        *context,
                                bool         previous)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set shared pointers from the main domain members
 *
 * \param[in]  quant       additional mesh quantities struct.
 * \param[in]  connect     pointer to a \ref cs_cdo_connect_t struct.
 * \param[in]  time_step   pointer to a \ref cs_time_step_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_init_sharing(const cs_cdo_quantities_t     *quant,
                             const cs_cdo_connect_t        *connect,
                             const cs_time_step_t          *time_step)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free shared pointers with lifecycle dedicated to this file
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_finalize_sharing(void)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize a cs_cdocb_scaleq_t structure storing data useful
 *         for building and managing such a scheme
 *
 * \param[in]      eqp         pointer to a \ref cs_equation_param_t structure
 * \param[in]      var_id      id of the variable field
 * \param[in]      bflux_id    id of the boundary flux field
 * \param[in, out] eqb         pointer to a \ref cs_equation_builder_t struct.
 *
 * \return a pointer to a new allocated cs_cdocb_scaleq_t structure
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdocb_scaleq_init_context(const cs_equation_param_t   *eqp,
                             int                          var_id,
                             int                          bflux_id,
                             cs_equation_builder_t       *eqb)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Destroy a \ref cs_cdocb_scaleq_t structure
 *
 * \param[in] scheme_context   pointer to a scheme context structure to free
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

void *
cs_cdocb_scaleq_free_context(void       *scheme_context)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the boundary conditions known from the settings
 *         Define an indirection array for the enforcement of internal DoFs
 *         only if needed.
 *         Case of scalar-valued CDO-Cb schemes
 *
 * \param[in]      t_eval          time at which one evaluates BCs
 * \param[in]      mesh            pointer to a cs_mesh_t structure
 * \param[in]      eqp             pointer to a cs_equation_param_t structure
 * \param[in, out] eqb             pointer to a cs_equation_builder_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_setup(cs_real_t                      t_eval,
                      const cs_mesh_t               *mesh,
                      const cs_equation_param_t     *eqp,
                      cs_equation_builder_t         *eqb)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the initial values of the variable field taking into account
 *         the boundary conditions.
 *         Case of scalar-valued CDO-Cb schemes.
 *
 * \param[in]      t_eval     time at which one evaluates BCs
 * \param[in]      field_id   id related to the variable field of this equation
 * \param[in]      mesh       pointer to a cs_mesh_t structure
 * \param[in]      eqp        pointer to a cs_equation_param_t structure
 * \param[in, out] eqb        pointer to a cs_equation_builder_t structure
 * \param[in, out] context    pointer to the scheme context (cast on-the-fly)
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_init_values(cs_real_t                     t_eval,
                            const int                     field_id,
                            const cs_mesh_t              *mesh,
                            const cs_equation_param_t    *eqp,
                            cs_equation_builder_t        *eqb,
                            void                         *context)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the balance for an equation over the full computational
 *         domain
 *         Case of scalar-valued CDO Cell-based scheme
 *
 * \param[in]      eqp       pointer to a \ref cs_equation_param_t structure
 * \param[in, out] eqb       pointer to a \ref cs_equation_builder_t structure
 * \param[in, out] context   pointer to a scheme builder structure
 *
 * \return a pointer to a \ref cs_cdo_balance_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_balance_t *
cs_cdocb_scaleq_balance(const cs_equation_param_t     *eqp,
                        cs_equation_builder_t         *eqb,
                        void                          *context)
{
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Solve the steady Stokes or Oseen system with a CDO cell-based scheme
 *         using a monolithic approach and GKB algorithm
 *
 * \param[in]      mesh            pointer to a \ref cs_mesh_t structure
 * \param[in, out] scheme_context  pointer to a structure cast on-the-fly
 */
/*----------------------------------------------------------------------------*/

void
cs_cdocb_scaleq_solve_steady_state(bool                        cur2prev,
                                   const cs_mesh_t            *mesh,
                                   const int                   field_id,
                                   const cs_equation_param_t  *eqp,
                                   cs_equation_builder_t      *eqb,
                                   void                       *context)
{
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
