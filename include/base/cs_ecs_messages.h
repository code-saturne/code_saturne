/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_ECS_MESSAGES_H__
#define __CS_ECS_MESSAGES_H__

/*============================================================================
 * Exchange of data between Code_Saturne Kernel and Preprocessor
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 *  Public functions definition for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Receive messages from the pre-processor about the dimensions of mesh
 * parameters
 *
 * FORTRAN Interface:
 *
 * SUBROUTINE LEDEVI(NOMRUB, TYPENT, NBRENT, TABENT)
 * *****************
 *
 * INTEGER          NDIM        : <-- : Spacial dimension (3)
 * INTEGER          NFML        : <-- : Number of families (group classes)
 * INTEGER          NPRFML      : <-- : Number of properties per family
 * INTEGER          IPERIO      : <-- : Periodicity indicator
 * INTEGER          IPEROT      : <-- : Number of rotation periodicities
 *----------------------------------------------------------------------------*/

void CS_PROCF(ledevi, LEDEVI)
(
 cs_int_t   *ndim,
 cs_int_t   *nfml,
 cs_int_t   *nprfml,
 cs_int_t   *iperio,
 cs_int_t   *iperot
);

/*============================================================================
 *  Public functions definition
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read data from the pre-processor and finalize pre-processor input.
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-- pointer to mesh builder structure
 *
 * returns:
 *----------------------------------------------------------------------------*/

void
cs_ecs_messages_read_data(cs_mesh_t          *mesh,
                          cs_mesh_builder_t  *mesh_builder);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_ECS_MESSAGES_H__ */

