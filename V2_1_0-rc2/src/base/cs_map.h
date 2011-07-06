/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2010 EDF S.A., France
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

#ifndef __CS_MAP_H__
#define __CS_MAP_H__

/*============================================================================
 * Map helper structures
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_map_name_to_id_t  cs_map_name_to_id_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create empty name to id map.
 *
 * returns:
 *   pointer to newly initialized map structure.
 *----------------------------------------------------------------------------*/

cs_map_name_to_id_t *
cs_map_name_to_id_create(void);

/*----------------------------------------------------------------------------
 * Destroy name to id map structure.
 *
 * parameters:
 *   m <-> pointer to map structure.
 *----------------------------------------------------------------------------*/

void
cs_map_name_to_id_destroy(cs_map_name_to_id_t **m);

/*----------------------------------------------------------------------------
 * Find id matching a key, inserting key if not already present.
 *
 * parameters:
 *   m     <-> pointer to map structure
 *   key   <-- character string (key)
 *
 * returns:
 *   id matching key (already present or newly inserted)
 *----------------------------------------------------------------------------*/

int
cs_map_name_to_id(cs_map_name_to_id_t  *m,
                  const char           *key);

/*----------------------------------------------------------------------------
 * Return the size of a map.
 *
 * parameters:
 *   m <-- pointer to map structure.
 *
 * returns:
 *   number of entries in map.
 *----------------------------------------------------------------------------*/

size_t
cs_map_name_to_id_size(const cs_map_name_to_id_t *m);

/*----------------------------------------------------------------------------
 * Return a given key in a map.
 *
 * parameters:
 *   m  <-- pointer to map structure.
 *   id <-- key id
 *
 * returns:
 *   pointer to key.
 *----------------------------------------------------------------------------*/

const char *
cs_map_name_to_id_key(const cs_map_name_to_id_t  *m,
                      size_t                      id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MAP_H__ */
