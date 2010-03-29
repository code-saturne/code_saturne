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

#ifndef __CS_SYR_COUPLING_H__
#define __CS_SYR_COUPLING_H__

/*============================================================================
 * SYRTHES coupling
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

/*============================================================================
 *  Global variables definition
 *============================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define new SYRTHES coupling.
 *
 * In the case of a single Code_Saturne and single SYRTHES instance, the
 * syrthes_app_num and syrthes_name arguments are ignored.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances prioritarily based on the syrthes_name argument, then
 * on the syrthes_app_num argument. If syrthes_name is empty, matching will
 * be based on syrthes_app_num only.
 *
 * Fortran Interface:
 *
 * SUBROUTINE DEFSY1
 * *****************
 *
 * INTEGER        syrthes_app_num   : <-- : application number of coupled
 *                                  :     : SYRTHES instance, or -1
 * CHARACTER*     syrthes_name      : <-- : name of coupled SYRTHES instance
 * CHARACTER      projection_axis   : <-- : ' ' for 3D, 'x', 'y', or 'z'
 *                                  :     : for 2D projection
 * CHARACTER*     boundary_criteria : <-- : boundary face selection criteria,
 *                                  :     : empty if no boundary coupling
 * CHARACTER*     volume_criteria   : <-- : volume cell selection criteria,
 *                                  :     : empty if no volume coupling
 * INTEGER        verbosity         : <-- : verbosity level
 * INTEGER        syrthes_n_len     : <-- : length of syrthes_name
 * INTEGER        boundary_c_len    : <-- : length of boundary_criteria
 * INTEGER        volume_c_len      : <-- : length of volume_criteria
 *----------------------------------------------------------------------------*/

void CS_PROCF(defsy1, DEFSY1)
(
 cs_int_t    *syrthes_app_num,
 const char  *syrthes_name,
 char        *projection_axis,
 const char  *boundary_criteria,
 const char  *volume_criteria,
 cs_int_t    *verbosity,
 cs_int_t    *syrthes_n_len,
 cs_int_t    *boundary_c_len,
 cs_int_t    *volume_c_len
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Get number of SYRTHES couplings.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBCSYR
 * *****************
 *
 * INTEGER          n_couplings     : <-- : number of SYRTHES couplings
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbcsyr, NBCSYR)
(
 cs_int_t  *const n_couplings
);

/*----------------------------------------------------------------------------
 * Create nodal coupled mesh.
 * Send vertices's coordinates and connectivity of coupled mesh.
 *
 * Fortran Interface:
 *
 * SUBROUTINE GEOSYR
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(geosyr, GEOSYR)
(
 void
);

/*----------------------------------------------------------------------------
 * Check if SYRTHES coupling continues or if we must stop.
 *
 * For SYRTHES 4 couplings, if nt_cur_abs < nt_max_abs, a new iteration
 * message is sent; otherwise, a stop message is sent. A corresponding
 * message is received, and if it is a stop message, nt_max_abs is
 * set to nt_cur_abs.
 *
 * For SYRTHES 3 couplings, A message (stop or new iteration) is
 * received. No iteration start message is sent, as this is done
 * by ITDSYR.
 *
 * Fortran Interface:
 *
 * SUBROUTINE TSTSYR (IMSFIN)
 * *****************
 *
 * INTEGER          NTMABS      : <-> : Maximum iteration number
 * INTEGER          NTCABS      : --> : Current iteration numbern
 *----------------------------------------------------------------------------*/

void CS_PROCF(tstsyr, TSTSYR)
(
 cs_int_t *ntmabs,
 cs_int_t *ntcabs
);

/*----------------------------------------------------------------------------
 * Synchronize new time step message for SYRTHES 3 couplings.
 *
 * This function is not necessary for SYRTHES 4, as all synchronization
 * is done by TSTSYR. For SYRTHES 3, it is necessary to distinguish
 * the last iteration from other iterations (to allow for SYRTHES 3 to
 * determine in advance that it will need to output postprocessing/restart
 * data), so using this separate function allows it to be placed after
 * MODPAR in the main time loop, in case NTMABS is changed by that function.
 *
 * Fortran Interface:
 *
 * SUBROUTINE ITDSYR (NTCABS, NTMABS)
 * *****************
 *
 * INTEGER          NTMABS      : --> : Maximum iteration number
 * INTEGER          NTCABS      : --> : Current iteration number
 *----------------------------------------------------------------------------*/

void CS_PROCF(itdsyr, ITDSYR)
(
 cs_int_t   *ntcabs,
 cs_int_t   *ntmabs
);

/*----------------------------------------------------------------------------
 * Get number of boundary faces coupled with SYRTHES.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBFSYR
 * *****************
 *
 * INTEGER          coupl_num       : --> : coupling number
 * INTEGER          n_coupl_faces   : <-- : number of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbfsyr, NBFSYR)
(
 const cs_int_t  *coupl_num,
       cs_int_t  *n_coupl_faces
);

/*----------------------------------------------------------------------------
 * Get local numbering of coupled faces
 *
 * Fortran interface:
 *
 * SUBROUTINE LFASYR
 * *****************
 *
 * INTEGER      coupl_num       : --> : coupling number
 * INTEGER      coupl_face_list : <-- : list of coupled boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF(lfasyr, LFASYR)
(
 const cs_int_t    *coupl_num,
       fvm_lnum_t  *coupl_face_list
);

/*----------------------------------------------------------------------------
 * Receive coupling variables from SYRTHES
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYI (NUMSYR, TWALL)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * DOUBLE PRECISION TWALL       : <-- : Wall temerature
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyi, VARSYI)
(
 cs_int_t   *numsyr,
 cs_real_t  *twall
);

/*----------------------------------------------------------------------------
 * Send coupling variables to SYRTHES
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYO (NUMSYR, TFLUID, HWALL)
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * DOUBLE PRECISION TFLUID      : --> : Fluid temperature
 * DOUBLE PRECISION HWALL       : --> : Exchange coefficient
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyo, VARSYO)
(
 cs_int_t   *numsyr,
 cs_real_t  *tfluid,
 cs_real_t  *hwall
);

/*----------------------------------------------------------------------------
 * Initialize post processing of SYRTHES couplings.
 *
 * Fortran Interface:
 *
 * SUBROUTINE PSTISY
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(pstisy, PSTISY)
(
 void
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define new SYRTHES coupling.
 *
 * In the case of a single Code_Saturne and single SYRTHES instance, the
 * syrthes_app_num and syrthes_name arguments are ignored.
 *
 * In case of multiple couplings, a coupling will be matched with available
 * SYRTHES instances prioritarily based on the syrthes_name argument, then
 * on the syrthes_app_num argument. If syrthes_name is empty, matching will
 * be based on syrthes_app_num only.
 *
 * arguments:
 *   syrthes_app_num   <-- number of SYRTHES application, or -1
 *   syrthes_name      <-- name of SYRTHES instance, or NULL
 *   boundary_criteria <-- boundary face selection criteria, or NULL
 *   volume_criteria   <-- volume cell selection criteria, or NULL
 *   projection_axis   <-- 'x', 'y', or 'y' for 2D projection axis (case
 *                         independent), or ' ' for standard 3D coupling
 *   verbosity         <-- verbosity level
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_define(int          syrthes_app_num,
                       const char  *syrthes_name,
                       const char  *boundary_criteria,
                       const char  *volume_criteria,
                       char         projection_axis,
                       int          verbosity);

/*----------------------------------------------------------------------------
 * Initialize SYRTHES couplings.
 *
 * This function may be called once all couplings have been defined,
 * and it will match defined couplings with available applications.
 *
 * parameters:
 *   port_num <-- port number for rank 0 to enable sockets,
 *                < 0 to disable sockets
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_init(int  port_num);

/*----------------------------------------------------------------------------
 * Finalize all SYRTHES couplings.
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_all_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SYR_COUPLING_H__ */
