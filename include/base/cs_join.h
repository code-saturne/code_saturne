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
 * Structure and function headers handling with joining operation
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_join_util.h"

/*---------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *===========================================================================*/

/*============================================================================
 * Type definition
 *===========================================================================*/

typedef struct {

  cs_join_param_t   param;      /* Set of parameters used to control
                                   the joining operations */

  char             *criteria;   /* Criteria used to select border faces
                                   implied in the joining operation */

} cs_join_t;

/*=============================================================================
 * Global variables
 *===========================================================================*/

extern cs_int_t  cs_glob_n_joinings;
extern cs_join_t  **cs_glob_join_array;

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
 * INTEGER        verbosity        : <-- : verbosity level
 * INTEGER        joining_c_len    : <-- : length of joining_criteria
 *----------------------------------------------------------------------------*/

void CS_PROCF(defjo1, DEFJO1)
(
 const char  *joining_criteria,
 cs_real_t   *fraction,
 cs_real_t   *plane,
 cs_int_t    *verbosity,
 cs_int_t    *joining_c_len
 CS_ARGF_SUPP_CHAINE
);

/*----------------------------------------------------------------------------
 * Set advanced parameters for the joining algorithm.
 *
 * Fortran Interface:
 *
 * SUBROUTINE SETAJP
 * *****************
 *
 * INTEGER      join_num          : <-- : join number
 * REAL         mtf               : <-- : merge tolerance coefficient
 * REAL         pmf               : <-- : pre-merge factor
 * INTEGER      tcm               : <-- : tolerance computation mode
 * INTEGER      icm               : <-- : intersection computation mode
 * INTEGER      maxbrk            : <-- : max number of tolerance reduction
 * INTEGER      max_sub_faces     : <-- : max. possible number of sub-faces
 *                                        by splitting a selected face
 * INTEGER      tml               : <-- : tree max level
 * INTEGER      tmb               : <-- : tree max boxes
 * REAL         tmr               : <-- : tree max ratio
 *---------------------------------------------------------------------------*/

void CS_PROCF(setajp, SETAJP)
(
 cs_int_t    *join_num,
 cs_real_t   *mtf,
 cs_real_t   *pmf,
 cs_int_t    *tcm,
 cs_int_t    *icm,
 cs_int_t    *maxbrk,
 cs_int_t    *max_sub_faces,
 cs_int_t    *tml,
 cs_int_t    *tmb,
 cs_real_t   *tmr
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
 *   verbosity     <-- level of verbosity required
 *---------------------------------------------------------------------------*/

void
cs_join_add(char   *sel_criteria,
            float   fraction,
            float   plane,
            int     verbosity);

/*----------------------------------------------------------------------------
 * Apply all the defined joining operations.
 *---------------------------------------------------------------------------*/

void
cs_join_all(void);

/*---------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_JOIN_H__ */
