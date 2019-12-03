#ifndef __CS_EXT_NEIGHBOR_H__
#define __CS_EXT_NEIGHBOR_H__

/*============================================================================
 * Fortran interfaces of functions needing a synchronization of the extended
 * neighborhood.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Extended neighborhood type
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_EXT_NEIGHBORHOOD_NONE,                  /* No extended neighborhood */
  CS_EXT_NEIGHBORHOOD_COMPLETE,              /* Full extended neighborhood */
  CS_EXT_NEIGHBORHOOD_CELL_CENTER_OPPOSITE,  /* Cell centers best aligned
                                                opposite to adjacent
                                                cell centers */
  CS_EXT_NEIGHBORHOOD_NON_ORTHO_MAX          /* Cells adjacent to faces
                                                whose non-orthogonality exceeds
                                                a given threshold */

} cs_ext_neighborhood_type_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/* Short names for extended neighborhood types */

extern const char *cs_ext_neighborhood_type_name[];

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the extended neighborhood type.
 *
 * \return  extended neighborhood type
 */
/*----------------------------------------------------------------------------*/

cs_ext_neighborhood_type_t
cs_ext_neighborhood_get_type(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the extended neighborhood type.
 *
 * \param[in]  enh_type  extended neighborhood type
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_set_type(cs_ext_neighborhood_type_t  enh_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get the non_orthogonality threshold  (in degrees) associated with the
 *        CS_EXT_NEIGHBORHOOD_NON_ORTHO_MAX neighborhood type.
 *
 * \return  non-orthogonality threshold
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_ext_neighborhood_get_non_ortho_max(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set the non_orthogonality threshold (in degrees) associated with the
 *        CS_EXT_NEIGHBORHOOD_NON_ORTHO_MAX neighborhood type.
 *
 * \param[in]  non_ortho_max  non-orthogonality threshold
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_set_non_ortho_max(cs_real_t  non_ortho_max);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Reduce the "cell -> cells" connectivity for the
 *        extended neighborhood using a non-orthogonality criterion.
 *
 * Note: Only cells sharing only a vertex or vertices (not a face)
 *       belong to the "cell -> cells" connectivity.
 *
 * \param[in]  mesh             pointer to mesh structure
 * \param[in]  mesh_quantities  associated mesh quantities
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_reduce(cs_mesh_t             *mesh,
                           cs_mesh_quantities_t  *mesh_quantities);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create the  "cell -> cells" connectivity.
 *
 * \param[in, out]  mesh  pointer to a mesh structure
 */
/*----------------------------------------------------------------------------*/

void
cs_ext_neighborhood_define(cs_mesh_t  *mesh);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_EXT_NEIGHBOR_H__ */
