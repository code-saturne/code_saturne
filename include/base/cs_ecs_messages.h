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

#include <stdarg.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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
 * INTEGER          NDIM        : <-- : Dimension de l'espace (3)
 * INTEGER          NFML        : <-- : Nombre de familles des faces de bord
 * INTEGER          NPRFML      : <-- : Nombre de propriétés max par famille
 * INTEGER          IPERIO      : <-- : Indicateur de périodicité
 * INTEGER          IPEROT      : <-- : Nombre de périodicités de rotation
 *----------------------------------------------------------------------------*/

void CS_PROCF(ledevi, LEDEVI)
(
 cs_int_t   *const ndim,    /* <-- dimension de l'espace                      */
 cs_int_t   *const nfml,    /* <-- nombre de familles des faces de bord       */
 cs_int_t   *const nprfml,  /* <-- nombre de propriétés max par famille       */
 cs_int_t   *const iperio,  /* <-- indicateur de périodicité                  */
 cs_int_t   *const iperot   /* <-- nombre de périodicités de rotation         */
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

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_ECS_MESSAGES_H__ */

