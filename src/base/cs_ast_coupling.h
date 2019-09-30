#ifndef __CS_AST_COUPLING_H__
#define __CS_AST_COUPLING_H__

/*============================================================================
 * Code_Aster coupling
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

typedef struct _cs_ast_coupling_t  cs_ast_coupling_t;

/*============================================================================
 * Global variables definition
 *============================================================================*/

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Send nodes coordinates and structure numbering of coupled mesh.
 *
 * Fortran Interface:
 *
 * SUBROUTINE ASTGEO
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(astgeo, ASTGEO)
(
 cs_int_t   *nbfast,
 cs_int_t   *lstfac,
 cs_int_t   *idfast,
 cs_int_t   *idnast,
 cs_real_t  *almax
);

/*----------------------------------------------------------------------------
 * Send stresses acting on the fluid/structure interface.
 *
 * Fortran Interface:
 *
 * SUBROUTINE ASTFOR
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(astfor, ASTFOR)
(
 cs_int_t    *ntcast,
 cs_int_t    *nbfast,
 cs_real_t   *forast
);

/*----------------------------------------------------------------------------
 * Receive displacement values of the fluid/structure interface
 *
 * Fortran Interface:
 *
 * SUBROUTINE ASTCIN
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(astcin, ASTCIN)
(
 cs_int_t    *ntcast,
 cs_real_3_t *disale
);

/*----------------------------------------------------------------------------
 * Receive coupling parameters
 *
 * Fortran Interface:
 *
 * SUBROUTINE ASTPAR
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(astpar, ASTPAR)
(
 cs_int_t  *nbpdt,
 cs_int_t  *nbsspdt,
 cs_real_t *delta,
 cs_real_t *tt,
 cs_real_t *dt
);

/*----------------------------------------------------------------------------
 * Exchange time step
 *
 * Fortran Interface:
 *
 * SUBROUTINE ASTPDT
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(astpdt, ASTPDT)
(
 cs_real_t *dttab,
 cs_int_t  *ncelet,
 cs_int_t  *nbpdt
);

/*----------------------------------------------------------------------------
 * Receive convergence value of Code_Saturne/Code_Aster coupling
 *
 * Fortran Interface:
 *
 * SUBROUTINE ASTCV1
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(astcv1, ASTCV1)
(
 cs_int_t  *ntcast,
 cs_int_t  *icv
);

/*-----------------------------------------------------------------------------
 * Send global convergence value of IFS calculations
 * (Internal and external structures)
 *
 * Fortran Interface:
 *
 * SUBROUTINE ASTCV2
 * *****************
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF(astcv2, ASTCV2)
(
 cs_int_t  *ntcast,
 cs_int_t  *icv
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Free coupling structure.
 *----------------------------------------------------------------------------*/

void
cs_ast_coupling_finalize(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_AST_COUPLING_H__ */
