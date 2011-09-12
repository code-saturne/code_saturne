/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2011 EDF S.A., France
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

#include "fvm_defs.h"

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
 * Test if the given SYRTHES coupling number is a surface coupling
 * Return 1 if true else 0
 *
 * Fortran Interface:
 *
 * SUBROUTINE TSURSY
 * *****************
 *
 * INTEGER          cplnum     : <-- : number of the SYRTHES coupling
 * INTEGER          issurf     : --> : 1 if surface coupling else 0
 *----------------------------------------------------------------------------*/

void CS_PROCF(tsursy, TSURSY)
(
 cs_int_t  *const cplnum,
 cs_int_t  *issurf
);

/*----------------------------------------------------------------------------
 * Test if the given SYRTHES coupling number is a volume coupling
 * Return 1 if true else 0
 *
 * Fortran Interface:
 *
 * SUBROUTINE TVOLSY
 * *****************
 *
 * INTEGER          cplnum     : <-- : number of the SYRTHES coupling
 * INTEGER          issurf     : --> : 1 if volume coupling else 0
 *----------------------------------------------------------------------------*/

void CS_PROCF(tvolsy, TVOLSY)
(
 cs_int_t  *const cplnum,
 cs_int_t  *isvol
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
 * Check if SYRTHES 3 couplings continue or if we must stop.
 *
 * For each SYRTHES 3 coupling, A message (stop or new iteration) is
 * received. No iteration start message is sent, as this is done
 * by ITDSYR.
 *
 * Fortran Interface:
 *
 * SUBROUTINE TSTSY3 (IMSFIN)
 * *****************
 *
 * INTEGER          NTMABS      : <-> : Maximum iteration number
 * INTEGER          NTCABS      : --> : Current iteration numbern
 *----------------------------------------------------------------------------*/

void CS_PROCF(tstsy3, TSTSY3)
(
 cs_int_t *ntmabs,
 cs_int_t *ntcabs
);

/*----------------------------------------------------------------------------
 * Synchronize new time step message for SYRTHES 3 couplings.
 *
 * For SYRTHES 3, it is necessary to distinguish the last iteration from
 * other iterations (to allow for SYRTHES 3 to determine in advance that it
 * will need to output postprocessing/restart data), so using this separate
 * function allows it to be placed after MODPAR in the main time loop,
 * in case NTMABS is changed by that function.
 *
 * Fortran Interface:
 *
 * SUBROUTINE ITDSY3 (NTCABS, NTMABS)
 * *****************
 *
 * INTEGER          NTMABS      : --> : Maximum iteration number
 * INTEGER          NTCABS      : --> : Current iteration number
 *----------------------------------------------------------------------------*/

void CS_PROCF(itdsy3, ITDSY3)
(
 cs_int_t   *ntcabs,
 cs_int_t   *ntmabs
);

/*----------------------------------------------------------------------------
 * Get number of coupled elements with SYRTHES.
 *
 * Fortran Interface:
 *
 * SUBROUTINE NBESYR
 * *****************
 *
 * INTEGER          coupl_num       : --> : coupling number
 * INTEGER          mode            : --> : 0 (surface); 1 (volume)
 * INTEGER          n_coupl_elts    : <-- : number of coupled elements
 *----------------------------------------------------------------------------*/

void CS_PROCF(nbesyr, NBESYR)
(
 const cs_int_t  *coupl_num,
 const cs_int_t  *mode,
       cs_int_t  *n_coupl_elts
);

/*----------------------------------------------------------------------------
 * Get local numbering of coupled elements
 *
 * Fortran interface:
 *
 * SUBROUTINE LELTSY
 * *****************
 *
 * INTEGER      coupl_num       : --> : coupling number
 * INTEGER      mode            : --> : 0 (surface); 1 (volume)
 * INTEGER      coupl_elt_list  : <-- : list of coupled elements
 *----------------------------------------------------------------------------*/

void CS_PROCF(leltsy, LELTSY)
(
 const cs_int_t    *coupl_num,
 const cs_int_t    *mode,
       cs_lnum_t   *coupl_elt_list
);

/*----------------------------------------------------------------------------
 * User function wrapper for definition of SYRTHES couplings
 *
 * Fortran Interface:
 *
 * SUBROUTINE USSYRC
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (ussyrc, USSYRC)
(
 void
);

/*----------------------------------------------------------------------------
 * Receive coupling variables from SYRTHES
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYI
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * INTEGER          MODE        : --> : 0 (surface); 1 (volume)
 * DOUBLE PRECISION TSOLID      : <-- : Solid temperature
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyi, VARSYI)
(
 cs_int_t   *numsyr,
 cs_int_t   *mode,
 cs_real_t  *tsolid
);

/*----------------------------------------------------------------------------
 * Send coupling variables to SYRTHES
 *
 * Fortran Interface:
 *
 * SUBROUTINE VARSYO
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * INTEGER          MODE        : --> : 0 (surface); 1 (volume)
 * INTEGER          LSTELT      : --> : List of coupled elements
 * DOUBLE PRECISION TFLUID      : --> : Fluid temperature
 * DOUBLE PRECISION HFLUID      : --> : Exchange coefficient
 *----------------------------------------------------------------------------*/

void CS_PROCF (varsyo, VARSYO)
(
 cs_int_t   *numsyr,
 cs_int_t   *mode,
 cs_int_t   *lstelt,
 cs_real_t  *tfluid,
 cs_real_t  *hfluid
);

/*----------------------------------------------------------------------------
 * Compute the explicit/implicit contribution to source terms in case of
 * volume coupling with SYRTHES4
 *
 * Fortran Interface:
 *
 * SUBROUTINE CTBVSY
 * *****************
 *
 * INTEGER          NUMSYR      : --> : Number of SYRTHES coupling
 * DOUBLE PRECISION TFLUID      : --> : Fluid temperature
 * DOUBLE PRECISION CTBIMP      : <-> : Implicit contribution
 * DOUBLE PRECISION CTBEXP      : <-> : Explicit contribution
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
 * Define new SYRTHES coupling.
 *
 * In the case of a single Code_Saturne and single SYRTHES instance, the
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
 *   verbosity         <-- verbosity level
 *   visualization     <-- visualization output level (0 or 1)
 *----------------------------------------------------------------------------*/

void
cs_syr_coupling_define(const char  *syrthes_name,
                       const char  *boundary_criteria,
                       const char  *volume_criteria,
                       char         projection_axis,
                       int          verbosity,
                       int          visualization);

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

END_C_DECLS

#endif /* __CS_SYR_COUPLING_H__ */
