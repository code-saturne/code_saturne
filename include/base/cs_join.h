/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 2008-2009 EDF S.A., France
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

#ifndef __CS_JOIN_H__
#define __CS_JOIN_H__

/*============================================================================
 * Structure and function headers handling with ghost cells
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *===========================================================================*/

/*============================================================================
 * Type definition
 *===========================================================================*/

/*=============================================================================
 * Global variables
 *===========================================================================*/

/*============================================================================
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define new boundary faces joining.
 *
 * Fortran Interface:
 *
 * SUBROUTINE DEFJO1
 * *****************
 *
 * CHARACTER*     joining_criteria : <-- : boundary face selection criteria,
 * REAL           fraction         : <-- : parameter for merging vertices
 * REAL           plane            : <-- : parameter for splitting faces
 * REAL           rtf              : <-- : reduction of tolerance factor
 * REAL           mtf              : <-- : merge tolerance coefficient
 * REAL           etf              : <-- : equivalence tolerance coefficient
 * INTEGER        verbosity        : <-- : verbosity level
 * INTEGER        joining_c_len    : <-- : length of joining_criteria
 *----------------------------------------------------------------------------*/

void CS_PROCF(defjo1, DEFJO1)
(
 const char  *joining_criteria,
 cs_real_t   *fraction,
 cs_real_t   *plane,
 cs_real_t   *rtf,
 cs_real_t   *mtf,
 cs_real_t   *etf,
 cs_int_t    *verbosity,
 cs_int_t    *joining_c_len
 CS_ARGF_SUPP_CHAINE
 );

/*=============================================================================
 * Public function prototypes
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Create and initialize a cs_join_t structure.
 *
 * parameters:
 *   sel_criteria  <-- boundary face selection criteria
 *   fraction      <-- value of the fraction parameter
 *   plane         <-- value of the plane parameter
 *   rtf           <-- value of the "reduction tolerance factor" parameter
 *   mtf           <-- value of the "merge tolerance factor" parameter
 *   etf           <-- value of the "edge equiv. tolerance factor" parameter
 *   max_sub_faces <-- value of the max_sub_faces
 *   tml           <-- value of the "tree max level" parameter
 *   tmb           <-- value of the "tree max boxes" parameter
 *   tmr           <-- value of the "tree max ratio" parameter
 *   verbosity     <-- level of verbosity required
 *
 * returns:
 *   a pointer to a cs_join_t structure.
 *---------------------------------------------------------------------------*/

void
cs_join_add(char   *sel_criteria,
            float   fraction,
            float   plane,
            float   rtf,
            float   mtf,
            float   etf,
            int     max_sub_faces,
            int     tml,
            int     tmb,
            float   tmr,
            int     verbosity);

/*----------------------------------------------------------------------------
 * Apply all the defined joining operations.
 *---------------------------------------------------------------------------*/

void
cs_join_all(void);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_H__ */
