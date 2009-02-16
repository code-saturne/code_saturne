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

#ifndef __CS_MESH_SOLCOM_H__
#define __CS_MESH_SOLCOM_H__

/*============================================================================
 * Read a mesh in "SolCom" format
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Public functions prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Update mesh size information after reading a "SolCom" format file header.
 *
 * Fortran interface:
 *
 * SUBROUTINE DIMGEO (NDIM  , NCELET, NCEL  , NFAC  , NFABOR, NSOM  ,
 * *****************
 *                    LNDFAC, LNDFBR, NFML  , NPRFML,
 *                    NTETRA, NPYRAM, NPRISM, NHEXAE )
 *
 * INTEGER          NDIM        : <-- : spatial dimension (3)
 * INTEGER          NCELET      : <-- : number of extended cells
 * INTEGER          NCEL        : <-- : number of true cells
 * INTEGER          NFAC        : <-- : number of interior faces
 * INTEGER          NFABOR      : <-- : number of boundary faces
 * INTEGER          NSOM        : <-- : number of vertices (optional)
 * INTEGER          LNDFAC      : <-- : length of SOMFAC (optional)
 * INTEGER          LNDFBR      : <-- : length of SOMFBR (optional)
 * INTEGER          NFML        : <-- : number of families
 * INTEGER          NPRFML      : <-- : max. number of properties per family
 * INTEGER          NTETRA      : <-- : number of tetrahedra
 * INTEGER          NPYRAM      : <-- : number of pyramids
 * INTEGER          NPRISM      : <-- : number of prisms
 * INTEGER          NHEXAE      : <-- : number of hexahedra
 *----------------------------------------------------------------------------*/

void CS_PROCF (dimgeo, DIMGEO)
(
 const cs_int_t   *ndim,
 const cs_int_t   *ncelet,
 const cs_int_t   *ncel,
 const cs_int_t   *nfac,
 const cs_int_t   *nfabor,
 const cs_int_t   *nsom,
 const cs_int_t   *lndfac,
 const cs_int_t   *lndfbr,
 const cs_int_t   *nfml,
 const cs_int_t   *nprfml,
 const cs_int_t   *ntetra,
 const cs_int_t   *npyram,
 const cs_int_t   *nprism,
 const cs_int_t   *nhexae
);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Read a mesh in "SolCom" format (prior to Code_Saturne 1.0)
 *
 * parameters:
 *   mesh            <-- associated mesh
 *   mesh_quantities <-- associated quantities
 *----------------------------------------------------------------------------*/

void
cs_mesh_solcom_read(cs_mesh_t             *mesh,
                    cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MESH_SOLCOM_H__ */
