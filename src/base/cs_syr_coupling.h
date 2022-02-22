#ifndef __CS_SYR_COUPLING_H__
#define __CS_SYR_COUPLING_H__

/*============================================================================
 * SYRTHES coupling
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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_defs.h"

#include "cs_base.h"
#include "cs_zone.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Structure definition
 *============================================================================*/

/*============================================================================
 *  Global variables definition
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define new SYRTHES coupling.
 *
 * In the case of a single code_saturne and single SYRTHES instance, the
 * syrthes_name argument is ignored.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances based on the syrthes_name argument.
 *
 * arguments:
 *   syrthes_name      <-- name of SYRTHES instance
 *   boundary_criteria <-- boundary face selection criteria, or NULL
 *   volume_criteria   <-- volume cell selection criteria, or NULL
 *   projection_axis   <-- 'x', 'y', or 'y' for 2D projection axis (case
 *                         independent), or ' ' for standard 3D coupling
 *   allow_nonmatching <-- allow nearest-neighbor mapping where matching
 *                         within tolerance is not available
 *   tolerance         <-- addition to local extents of each element
 *                         extent = base_extent * (1 + tolerance)
 *   verbosity         <-- verbosity level
 *   visualization     <-- visualization output level (0 or 1)
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_define(const char  *syrthes_name,
                       const char  *boundary_criteria,
                       const char  *volume_criteria,
                       char         projection_axis,
                       bool         allow_nonmatching,
                       float        tolerance,
                       int          verbosity,
                       int          visualization);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Associated a zone to a defined SYRTHES coupling.
 *
 * \param[in] syrthes_name  matching SYRTHES application name
 * \param[in] z             pointer to matching zone
 *                          (boundary or volume)
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_add_zone(const char       *syrthes_name,
                         const cs_zone_t  *z);

/*----------------------------------------------------------------------------
 * Initialize SYRTHES couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_init(void);

/*----------------------------------------------------------------------------
 * Finalize all SYRTHES couplings.
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_finalize(void);

/*----------------------------------------------------------------------------
 * Return number of SYRTHES couplings.
 *
 * return:
 *   number of SYRTHES couplings defined
 *----------------------------------------------------------------------------*/

int
cs_syr_coupling_n_couplings(void);

/*----------------------------------------------------------------------------
 * Set conservativity forcing flag to True (1) or False (0) for all defined
 * SYRTHES couplings
 *
 * parameter:
 *   flag     <--  Conservativity forcing flag to set
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_set_conservativity(int  flag);

/*----------------------------------------------------------------------------
 * Set explicit treatment for the source terms in SYRTHES volume couplings
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_set_explicit_treatment(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log SYRTHES coupling setup information.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create coupled meshes and setup PLE locator for Syrthes couplings.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_init_meshes(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Check if the given SYRTHES coupling number is a surface couplings.
 *
 * \param[in] cpl_id   matching SYRTHES coupling id
 *
 * \return 1 if the coupling includes the surface, 0 otherwise.
 */
/*----------------------------------------------------------------------------*/

int
cs_syr_coupling_is_surf(int  cpl_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read boundary field/variable values relative to a SYRTHES coupling.
 *
 * \param[in]       nvar     number of variables
 * \param[in]       bc_type  boundary condition type
 * \param[in, out]  icodcl   boundary condition codes
 * \param[in, out]  rcodcl   boundary condition values
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_recv_boundary(int        nvar,
                              int        bc_type[],
                              int        icodcl[],
                              cs_real_t  rcodcl[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Send field/variable values relative to a SYRTHES coupling.
 *
 * \param[in]  h_wall   wall thermal exchange coefficient
 * \param[in]  v_fluid  near-wall fluid thermal variable
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_send_boundary(const cs_real_t  h_wall[],
                              cs_real_t        v_fluid[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Exchange volume values relative to a SYRTHES coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_exchange_volume(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the source term (implicit and/or explicit part) for a
 *         volume coupling with SYRTHES.
 *
 * \param[in]       field_id  field id
 * \param[in, out]  st_exp    explicit source term
 * \param[in, out]  st_imp    implicit source term
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_volume_source_terms(int        field_id,
                                    cs_real_t  st_exp[],
                                    cs_real_t  st_imp[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get number of coupled elements with SYRTHES.
 *
 * \param[in]   cpl_id  coupling id
 * \param[in]   mode    0 for boundary, 1 for volume
 *
 * \return  number of coupled elements for this coupling
 */
/*----------------------------------------------------------------------------*/

cs_lnum_t
cs_syr_coupling_n_elts(int  cpl_id,
                       int  mode);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get local ids of elements coupled with SYRTHES
 *
 * \param[in]    cpl_id   coupling id
 * \param[in]    mode     0 for boundary, 1 for volume
 * \param[out]   elt_ids  ids of coupled elements (preallocated)
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_elt_ids(int        cpl_id,
                        int        mode,
                        cs_lnum_t  elt_ids[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Receive coupling variables from SYRTHES.
 *
 * \param[in]    cpl_id   coupling id
 * \param[in]    mode     0 for boundary, 1 for volume
 * \param[out]   t_solid  solid temperature
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_recv_tsolid(int        cpl_id,
                            int        mode,
                            cs_real_t  t_solid[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Send coupling variables to SYRTHES.
 *
 * \param[in]    cpl_id   coupling id
 * \param[in]    mode     0 for boundary, 1 for volume
 * \param[in]    elt_ids  ids of coupled elements
 * \param[in]    t_fluid  fluid temperature
 * \param[in]    h_fluid  fluid exchange coefficient
 */
/*----------------------------------------------------------------------------*/

void
cs_syr_coupling_send_tf_hf(int              cpl_id,
                           int              mode,
                           const cs_lnum_t  elt_ids[],
                           cs_real_t        t_fluid[],
                           cs_real_t        h_fluid[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SYR_COUPLING_H__ */
