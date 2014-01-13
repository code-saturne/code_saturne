#ifndef __CS_LAGR_PERIO_H__
#define __CS_LAGR_PERIO_H__

/*============================================================================
 * Management of the periodicity for particles
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Public functions definition for API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Build buffers to keep the link between a given halo cell and :
 *  - the related real cell
 *  - the transformation id
 *
 * Fortran Interface :
 *
 * SUBROUTINE PERLOC
 * *****************
 *
 * INTEGER ICELCR(NCELET-NCEL) : <-  : related real cell buffer
 * INTEGER IPERCR(NCELET-NCEL) : <-  : transformation id buffer
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (perloc, PERLOC)(cs_int_t   *icelcr,
                          cs_int_t   *ipercr);

/*----------------------------------------------------------------------------
 * Apply rotation to the location of a particle.
 *
 * Fortran Interface :
 *
 * SUBROUTINE LAGPER
 * *****************
 *
 * INTEGER          ITRANS        :  -> : transformation id buffer
 * DOUBLE PRECISION VTX_A         :  -> : location of vertex before transform.
 * DOUBLE PRECISION VTX_B         : <-  : location of the vertex after transform.
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lagper, LAGPER)(const cs_int_t   *itrans,
                          const cs_real_t   vtx_a[],
                                cs_real_t   vtx_b[]);

/*----------------------------------------------------------------------------
 * Apply rotation on the velocity vector of a particle.
 *
 * Fortran Interface :
 *
 * SUBROUTINE LAGVEC
 * *****************
 *
 * INTEGER          ITRANS        :  -> : transformation id
 * DOUBLE PRECISION VECTI         :  -> : vector before transformation
 * DOUBLE PRECISION VECTF         : <-  : vector after transformation
 *
 * Returns:
 *----------------------------------------------------------------------------*/

void
CS_PROCF (lagvec, LAGVEC)(const cs_int_t   *itrans,
                          const cs_real_t   vecti[],
                                cs_real_t   vectf[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_PERIO_H__ */
