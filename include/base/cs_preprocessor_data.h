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

#ifndef __CS_PRE_PROCESSOR_DATA_H__
#define __CS_PRE_PROCESSOR_DATA_H__

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
 *  Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query or modification of the option for domain partitioning when no
 * partitioning file is present.
 *
 * This function returns 1 or 2 according to the selected algorithm.
 *
 * Fortran interface :
 *
 * subroutine algdom (iopt)
 * *****************
 *
 * integer          iopt        : <-> : Choice of the partitioning base
 *                                        0: query
 *                                        1: initial numbering
 *                                        2: space-filling curve (default)
 *----------------------------------------------------------------------------*/

void
CS_PROCF (algdom, ALGDOM)(cs_int_t  *iopt);

/*----------------------------------------------------------------------------
 * Receive messages from the pre-processor about the dimensions of mesh
 * parameters
 *
 * Fortran Interface:
 *
 * subroutine ledevi(ndim   , nfml  , nprfml, iperio, iperot)
 * *****************
 *
 * integer          ndim        : <-- : Spacial dimension (3)
 * integer          nfml        : <-- : Number of families
 * integer          nprfml      : <-- : Number of properties per family
 * integer          iperio      : <-- : Periodicity indicator
 * integer          iperot      : <-- : Number of rotation periodicities
 *----------------------------------------------------------------------------*/

void
CS_PROCF(ledevi, LEDEVI)(cs_int_t   *ndim,
                         cs_int_t   *nfml,
                         cs_int_t   *nprfml,
                         cs_int_t   *iperio,
                         cs_int_t   *iperot);

/*============================================================================
 *  Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Query or modification of the option for domain partitioning when no
 * partitioning file is present.
 *
 *  0 : query
 *  1 : partition based on initial numbering
 *  2 : partition based on space-filling curve (default)
 *
 * choice <-- of partitioning algorithm.
 *
 * returns:
 *   1 or 2 according to the selected algorithm.
 *----------------------------------------------------------------------------*/

int
cs_preprocessor_data_part_choice(int choice);

/*----------------------------------------------------------------------------
 * Read pre-processor mesh data and finalize input.
 *
 * parameters:
 *   mesh         <-- pointer to mesh structure
 *   mesh_builder <-- pointer to mesh builder structure
 *
 * returns:
 *----------------------------------------------------------------------------*/

void
cs_preprocessor_data_read_mesh(cs_mesh_t          *mesh,
                               cs_mesh_builder_t  *mesh_builder);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PRE_PROCESSOR_DATA_H__ */

