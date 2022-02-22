#ifndef __CS_CDO_TOOLBOX_H__
#define __CS_CDO_TOOLBOX_H__

/*============================================================================
 * Set of toolbox functions: shared buffer, balance, synchronization
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

#include "cs_cdo_connect.h"
#include "cs_flag.h"
#include "cs_xdef.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* Strategy of synchronization of values shared across several cells
 * This applies to vertices, edges and faces
 *
 * CS_CDO_SYNC_ZERO_VALUE
 * If zero is a possible value then set this value, otherwise one takes
 * the mean-value
 *
 * CS_CDO_SYNC_MEAN_VALUE
 * Compute the mean-value across values to set
 *
 */

#define  CS_CDO_SYNC_ZERO_VALUE    1
#define  CS_CDO_SYNC_MEAN_VALUE    2

/*============================================================================
 * Type definitions
 *============================================================================*/

/*
 * Structure used to store information generated during the analysis
 * of the balance of each term of an equation
 */

typedef struct {

  /* Where the balance is computed: primal vertices or primal cells */

  cs_flag_t       location;
  cs_lnum_t       size;
  cs_real_t      *balance;

  /* Balance for each main term */

  cs_real_t      *unsteady_term;
  cs_real_t      *reaction_term;
  cs_real_t      *diffusion_term;
  cs_real_t      *advection_term;
  cs_real_t      *source_term;
  cs_real_t      *boundary_term;

} cs_cdo_balance_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a pointer to a buffer of size at least the n_cells for
 *         managing temporary usage of memory when dealing with equations
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *         Set also shared pointers from the main domain members
 *
 * \param[in]  connect      pointer to a cs_cdo_connect_t structure
 * \param[in]  eb_flag      metadata for Edge-based schemes
 * \param[in]  fb_flag      metadata for Face-based schemes
 * \param[in]  vb_flag      metadata for Vertex-based schemes
 * \param[in]  vcb_flag     metadata for Vertex+Cell-basde schemes
 * \param[in]  hho_flag     metadata for HHO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_toolbox_init(const cs_cdo_connect_t       *connect,
                    cs_flag_t                     eb_flag,
                    cs_flag_t                     fb_flag,
                    cs_flag_t                     vb_flag,
                    cs_flag_t                     vcb_flag,
                    cs_flag_t                     hho_flag);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free buffers shared among the equations solved with CDO schemes
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_toolbox_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve a pointer to a buffer of size at least the 2*n_cells
 *         The size of the temporary buffer can be bigger according to the
 *         numerical settings
 *
 * \return  a pointer to an array of double
 */
/*----------------------------------------------------------------------------*/

cs_real_t *
cs_cdo_toolbox_get_tmpbuf(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the allocation size of the temporary buffer
 *
 * \return  the size of the temporary buffer
 */
/*----------------------------------------------------------------------------*/

size_t
cs_cdo_toolbox_get_tmpbuf_size(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cdo_balance_t structure
 *
 * \param[in]  location   where the balance is performed
 * \param[in]  size       size of arrays in the structure
 *
 * \return  a pointer to the new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_balance_t *
cs_cdo_balance_create(cs_flag_t    location,
                      cs_lnum_t    size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset a cs_cdo_balance_t structure
 *
 * \param[in, out] b     pointer to a cs_cdo_balance_t to reset
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_balance_reset(cs_cdo_balance_t   *b);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize balance terms if this is a parallel computation
 *
 * \param[in]      connect    pointer to a cs_cdo_connect_t structure
 * \param[in, out] balance    pointer to a cs_cdo_balance_t to sync
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_balance_sync(const cs_cdo_connect_t    *connect,
                    cs_cdo_balance_t          *balance);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cdo_balance_t structure
 *
 * \param[in, out]  p_balance  pointer to the pointer to free
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_balance_destroy(cs_cdo_balance_t   **p_balance);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the volumetric definitions to consider at each vertex
 *
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2v_idx   index array  to define
 * \param[in, out]  def2v_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_sync_vol_def_at_vertices(int                      n_defs,
                                cs_xdef_t              **defs,
                                cs_lnum_t                def2v_idx[],
                                cs_lnum_t                def2v_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the volumetric definitions to consider at each edge
 *
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2e_idx   index array  to define
 * \param[in, out]  def2e_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_sync_vol_def_at_edges(int                      n_defs,
                             cs_xdef_t              **defs,
                             cs_lnum_t                def2e_idx[],
                             cs_lnum_t                def2e_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Synchronize the volumetric definitions to consider at each face
 *
 * \param[in]       n_defs      number of definitions
 * \param[in]       defs        number of times the values has been updated
 * \param[in, out]  def2f_idx   index array  to define
 * \param[in, out]  def2f_ids   array of ids to define
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_sync_vol_def_at_faces(int                        n_defs,
                             cs_xdef_t                **defs,
                             cs_lnum_t                  def2f_idx[],
                             cs_lnum_t                  def2f_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the mean-value across ranks at each vertex
 *
 * \param[in]       dim         number of entries for each vertex
 * \param[in]       counter     number of occurences on this rank
 * \param[in, out]  values      array to update
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_sync_vertex_mean_values(int                         dim,
                               int                        *counter,
                               cs_real_t                  *values);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_TOOLBOX_H__ */
