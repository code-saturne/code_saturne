/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_SYR4_COUPLING_H__
#define __CS_SYR4_COUPLING_H__

/*============================================================================
 * Syrthes 4 coupling
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
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

/*============================================================================
 * Public function definitions
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
 *   app_num            <-- SYRTHES application number, or -1
 *   app_name           <-- SYRTHES application name, or NULL
 *   verbosity          <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_add(fvm_lnum_t   dim,
                     fvm_lnum_t   ref_axis,
                     const char  *face_sel_criterion,
                     const char  *cell_sel_criterion,
                     int          app_num,
                     const char  *app_name,
                     int          verbosity);

/*----------------------------------------------------------------------------
 * Destroy cs_syr4_coupling_t structures
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_all_destroy(void);

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
 * Exchange new iteration or stop message with SYRTHES.
 *
 * If nt_cur_abs < nt_max_abs, a new iteration message is sent;
 * otherwise, a stop message is sent. A corresponding message
 * is received, and if it is a stop message, nt_max_abs is
 * set to nt_cur_abs.
 *
 * parameters:
 *   nt_cur_abs <-- current iteration number
 *   nt_max_abs <-> maximum iteration number
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_sync_iter(int   nt_cur_abs,
                           int  *nt_max_abs);

/*----------------------------------------------------------------------------
 * Synchronize new time step
 *
 * parameters:
 *   nt_cur_abs <-- current iteration number
 *   nt_max_abs --> maximum iteration number
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_new_time_step(int  nt_cur_abs,
                               int  nt_max_abs);

/*----------------------------------------------------------------------------
 * Define coupled mesh and send it to SYRTHES
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_init_mesh(cs_syr4_coupling_t  *syr_coupling);

/*----------------------------------------------------------------------------
 * Get number of associated coupled faces in main mesh
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *
 * returns:
 *   number of vertices in coupled mesh
 *----------------------------------------------------------------------------*/

fvm_lnum_t
cs_syr4_coupling_get_n_faces(const cs_syr4_coupling_t  *syr_coupling);

/*----------------------------------------------------------------------------
 * Get local numbering of coupled faces
 *
 * parameters:
 *   syr_coupling    <-- SYRTHES coupling structure
 *   coupl_face_list --> List of coupled faces (1 to n)
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_get_face_list(const cs_syr4_coupling_t  *syr_coupling,
                               cs_int_t                   coupl_face_list[]);

/*----------------------------------------------------------------------------
 * Receive coupling variables from SYRTHES
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *   twall        --> wall temperature
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_recv_twall(cs_syr4_coupling_t *syr_coupling,
                            cs_real_t           twall[]);

/*----------------------------------------------------------------------------
 * Send coupling variables to SYRTHES
 *
 * parameters:
 *   syr_coupling <-- SYRTHES coupling structure
 *   tf           <-- fluid temperature
 *   hwall        <-- wall heat exchange coefficient (numerical, not physical)
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_send_tf_hwall(cs_syr4_coupling_t *syr_coupling,
                               cs_real_t           tf[],
                               cs_real_t           hwall[]);

/*----------------------------------------------------------------------------
 * Initialize post-processing of a SYRTHES coupling
 *
 * parameters:
 *   coupling_id --> Id of SYRTHES coupling
 *   writer_id   --> Id of associated writer
 *----------------------------------------------------------------------------*/

void
cs_syr4_coupling_post_init(int       coupling_id,
                           cs_int_t  writer_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SYR4_COUPLING_H__ */
