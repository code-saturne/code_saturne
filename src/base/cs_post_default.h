#ifndef __CS_POST_DEFAULT_H__
#define __CS_POST_DEFAULT_H__

/*============================================================================
 * Post-processing management
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2012 EDF S.A.

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

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public Fortran function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Configure the post-processing output so that a mesh displacement field
 * may be output automatically for meshes based on the global volume mesh/
 *
 * Fortran interface:
 *
 * subroutine pstdfm
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstdfm, PSTDFM)
(
 void
);

/*----------------------------------------------------------------------------
 * Update the "active" or "inactive" flag for writers based on the current
 * time step and their default output frequency.
 *
 * Fortran interface:
 *
 * subroutine pstntc (ntmabs, ntcabs, ttcabs)
 * *****************
 *
 * integer          ntmabs      : <-- : maximum time step number
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : absolute time at the current time step
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstntc, PSTNTC)
(
 const cs_int_t  *ntmabs,
 const cs_int_t  *ntcabs,
 const cs_real_t *ttcabs
);

/*----------------------------------------------------------------------------
 * Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * Fortran interface:
 *
 * subroutine pstact (numwri, indact)
 * *****************
 *
 * integer          numwri      : <-- : writer number, or 0 for all writers
 * integer          indact      : <-- : 0 to deactivate, 1 to activate
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstact, PSTACT)
(
 const cs_int_t  *numwri,
 const cs_int_t  *indact
);

/*----------------------------------------------------------------------------
 * Output post-processing meshes using associated writers.
 *
 * Fortran interface:
 *
 * subroutine pstema (ntcabs, ttcabs)
 * *****************
 *
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : current physical time
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstema, PSTEMA)
(
 const cs_int_t   *ntcabs,
 const cs_real_t  *ttcabs
);

/*----------------------------------------------------------------------------
 * Loop on post-processing meshes to output variables
 *
 * Fortran interface:
 *
 * subroutine pstvar
 * *****************
 *                  ( ntcabs,
 *                    nvar,   nscal,  nvlsta, nvisbr,
 *                    ttcabs,
 *                    dt,     rtpa,   rtp,    propce, propfa, propfb,
 *                    coefa,  coefb,
 *                    statce, stativ, statfb)
 *
 * integer          ntcabs      : --> : current time step number
 * integer          nvar        : <-- : number of variables
 * integer          nscal       : <-- : number of scalars
 * integer          nvlsta      : <-- : number of statistical variables (lagr)
 * integer          nvisbr      : <-- : number of boundary stat. variables (lagr)
 * double precision ttcabs      : <-- : current physical time
 * double precision dt          : <-- : local time step
 * double precision rtpa        : <-- : cell variables at previous time step
 * double precision rtp         : <-- : cell variables
 * double precision propce      : <-- : cell physical properties
 * double precision propfa      : <-- : interior face physical properties
 * double precision propfb      : <-- : boundary face physical properties
 * double precision coefa       : <-- : boundary conditions array
 * double precision coefb       : <-- : boundary conditions array
 * double precision statce      : <-- : cell statistics (lagrangian)
 * double precision stativ      : <-- : cell variance statistics (lagrangian)
 * double precision statfb      : <-- : boundary face statistics (lagrangian)
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *ntcabs,
 const cs_int_t   *nvar,
 const cs_int_t   *nscal,
 const cs_int_t   *nvlsta,
 const cs_int_t   *nvisbr,
 const cs_real_t  *ttcabs,
 const cs_real_t   dt[],
 const cs_real_t   rtpa[],
 const cs_real_t   rtp[],
 const cs_real_t   propce[],
 const cs_real_t   propfa[],
 const cs_real_t   propfb[],
 const cs_real_t   coefa[],
 const cs_real_t   coefb[],
 const cs_real_t   statce[],
 const cs_real_t   stativ[],
 const cs_real_t   statfb[]
);

/*----------------------------------------------------------------------------
 * Post-processing output of a variable defined on cells or faces of a mesh
 * using associated writers.
 *
 * fortran interface; use psteva (see cs_post_f2c.f90)
 *
 * subroutine pstev1 (nummai, nomvar, lnmvar, idimt,  ientla, ivarpr,
 * *****************
 *                    ntcabs, ttcabs, varcel, varfac, varfbr)
 *
 * integer          nummai      : <-- : number of associated output mesh
 * character        nomvar      : <-- : name of associated variable
 * integer          lnmvar      : <-- : variable name length
 * integer          idimt       : <-- : 1 for scalar, 3 for vector
 * integer          ientla      : <-- : if a vector, 1 for interlaced values
 *                              :     : (x1, y1, z1, x2, y2, ..., yn, zn),
 *                              :     : 0 otherwise (x1, x2, ...xn, y1, y2, ...)
 * integer          ivarpr      : <-- : 1 if variable is defined on "parent"
 *                              :     : mesh, 2 if defined on output mesh
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : current physical time
 * double precision varcel(*)   : <-- : cell values
 * double precision varfac(*)   : <-- : interior face values
 * double precision varfbo(*)   : <-- : boundary face values
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstev1, PSTEV1)
(
 const cs_int_t   *nummai,
 const char       *nomvar,
 const cs_int_t   *lnmvar,
 const cs_int_t   *idimt,
 const cs_int_t   *ientla,
 const cs_int_t   *ivarpr,
 const cs_int_t   *ntcabs,
 const cs_real_t  *ttcabs,
 const cs_real_t   varcel[],
 const cs_real_t   varfac[],
 const cs_real_t   varfbr[]
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * User override of default frequency or calculation end based output.
 *
 * Fortran interface:
 *
 * subroutine pstusn (ntmabs, ntcabs, ttcabs)
 * *****************
 *
 * integer          ntmabs      : <-- : maximum time step number
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : absolute time at the current time step
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstusn, PSTUSN)
(
 const cs_int_t  *ntmabs,
 const cs_int_t  *ntcabs,
 const cs_real_t *ttcabs
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POST_DEFAULT_H__ */
