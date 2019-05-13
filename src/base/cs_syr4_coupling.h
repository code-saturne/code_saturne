#ifndef __CS_SYR4_COUPLING_H__
#define __CS_SYR4_COUPLING_H__

/*============================================================================
 * Syrthes 4 coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Structure definition
 *============================================================================*/

/* Structure associated to Syrthes coupling */

typedef struct _cs_syr4_coupling_t  cs_syr4_coupling_t;

/*============================================================================
 *  Global variables definition
 *============================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute the implicit/explicit contribution for source terms in a SYRTHES
 * volume coupling
 *
 * Fortran Interface:
 *
 * SUBROUTINE CTBVSY (NUMSYR, TFLUID, CTBIMP, CTBEXP)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * DOUBLE PRECISION TFLUID      : --> : Fluid temperature
 * DOUBLE PRECISION CTBIMP      : --> : Implicit contribution
 * DOUBLE PRECISION CTBEXP      : --> : Explicit contribution
 *----------------------------------------------------------------------------*/

void CS_PROCF (ctbvsy, CTBVSY)
(
 cs_int_t   *numsyr,
 cs_real_t  *tfluid,
 cs_real_t  *ctbimp,
 cs_real_t  *ctbexp
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get number of SYRTHES couplings.
 *
 * returns:
 *   number of SYRTHES couplings
 *----------------------------------------------------------------------------*/

int
cs_syr4_coupling_n_couplings(void);

/*----------------------------------------------------------------------------
 * Get pointer to SYRTHES coupling.
 *
 * parameters:
 *   coupling_id <-- Id (0 to n-1) of SYRTHES coupling
 *
 * returns:
 *   pointer to SYRTHES coupling structure
 *----------------------------------------------------------------------------*/

cs_syr4_coupling_t *
cs_syr4_coupling_by_id(cs_int_t coupling_id);

/*----------------------------------------------------------------------------
 * Create a syr4_coupling_t structure.
 *
 * parameters:
 *   dim                <-- spatial mesh dimension
 *   ref_axis           <-- reference axis
 *   face_sel_criterion <-- criterion for selection of boundary faces
 *   cell_sel_criterion <-- criterion for selection of cells
 *   app_name           <-- SYRTHES application name
 *   allow_nonmatching  <-- nearest-neighbor search for non-matching faces flag
 *   tolerance          <-- addition to local extents of each element
 *                          extent = base_extent * (1 + tolerance)
 *   verbosity          <-- verbosity level
 *   visualization      <-- visualization output flag
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_add(int          dim,
                     int          ref_axis,
                     const char  *face_sel_criterion,
                     const char  *cell_sel_criterion,
                     const char  *app_name,
                     bool         allow_nonmatching,
                     float        tolerance,
                     int          verbosity,
                     int          visualization);

/*----------------------------------------------------------------------------
 * Destroy cs_syr4_coupling_t structures
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_all_destroy(void);

/*----------------------------------------------------------------------------
 * Set conservativity forcing flag to True (1) or False (0) for all defined
 * SYRTHES couplings
 *
 * parameter:
 *   flag     <--  Conservativity forcing flag to set
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_set_conservativity(int  flag);

/*----------------------------------------------------------------------------
 * Set explicit treatment for the source terms in SYRTHES volume couplings
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_set_explicit_treatment(void);

/*----------------------------------------------------------------------------
 * Initialize communicator for SYRTHES coupling
 *
 * parameters:
 *   syr_coupling  <-> Syrthes coupling structure
 *   coupling_id   <-- id of this coupling (for log file message)
 *   syr_root_rank <-- SYRTHES root rank
 *   n_syr_ranks   <-- Number of ranks associated with SYRTHES
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_init_comm(cs_syr4_coupling_t *syr_coupling,
                           int                 coupling_id,
                           int                 syr_root_rank,
                           int                 n_syr_ranks);

/*----------------------------------------------------------------------------
 * Define coupled mesh and send it to SYRTHES
 *
 * Optional post-processing output is also built at this stage.
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_init_mesh(cs_syr4_coupling_t  *syr_coupling);

/*----------------------------------------------------------------------------
 * Return 1 if this coupling is a surface coupling else 0
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *
 * returns:
 *   1 or 0
 *----------------------------------------------------------------------------*/

int
cs_syr4_coupling_is_surf(const cs_syr4_coupling_t  *syr_coupling);

/*----------------------------------------------------------------------------
 * Return 1 if this coupling is a volume coupling else 0
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *
 * returns:
 *   1 or 0
 *----------------------------------------------------------------------------*/

int
cs_syr4_coupling_is_vol(const cs_syr4_coupling_t  *syr_coupling);

/*----------------------------------------------------------------------------
 * Get number of associated coupled elements in main mesh
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *   mode          <-- 0 (surface); 1 (volume)
 *
 * returns:
 *   number of vertices in coupled mesh
 *----------------------------------------------------------------------------*/

cs_lnum_t
cs_syr4_coupling_get_n_elts(const cs_syr4_coupling_t *syr_coupling,
                            int                       mode);

/*----------------------------------------------------------------------------
 * Get local numbering of coupled elements
 *
 * parameters:
 *   syr_coupling  <-- SYRTHES coupling structure
 *   cpl_elt_lst   --> List of coupled elements (1 to n)
 *   mode          <-- 0 (surface); 1 (volume)
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_get_elt_list(const cs_syr4_coupling_t  *syr_coupling,
                              cs_int_t                   cpl_elt_lst[],
                              int                        mode);

/*----------------------------------------------------------------------------
 * Receive coupling variables from SYRTHES
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *   tsolid       --> solid temperature
 *   mode         <-- 0: surface coupling; 1: volume coupling
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_recv_tsolid(cs_syr4_coupling_t  *syr_coupling,
                             cs_real_t            tsolid[],
                             int                  mode);

/*----------------------------------------------------------------------------
 * Send coupling variables to SYRTHES
 *
 * parameters:
 *   syr_coupling  <-- SYRTHES coupling structure
 *   cpl_elt_list  <-- list of coupled boundary faces
 *   tf            <-- fluid temperature
 *   hf            <-- fluid heat exchange coef. (numerical or user-defined)
 *   mode          <-- 0 (surface); 1 (volume)
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_send_tf_hf(cs_syr4_coupling_t  *syr_coupling,
                            const cs_lnum_t      cpl_elt_lst[],
                            cs_real_t            tf[],
                            cs_real_t            hf[],
                            int                  mode);

/*----------------------------------------------------------------------------
 * Compute the explicit/implicit contribution to source terms in case of
 * volume coupling with SYRTHES4
 *
 * parameters:
 *   syr_coupling  <-- SYRTHES coupling structure
 *   tf            <-- fluid temperature
 *   ctbimp        <-> implicit contribution
 *   ctbexp        <-> explicit contribution
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_ts_contrib(const cs_syr4_coupling_t  *syr_coupling,
                            const cs_real_t            tf[],
                            cs_real_t                  ctbimp[],
                            cs_real_t                  ctbexp[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SYR4_COUPLING_H__ */
