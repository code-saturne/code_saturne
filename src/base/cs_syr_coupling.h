#ifndef __CS_SYR_COUPLING_H__
#define __CS_SYR_COUPLING_H__

/*============================================================================
 * SYRTHES coupling
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

#include "fvm_defs.h"

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

END_C_DECLS

#endif /* __CS_SYR_COUPLING_H__ */
