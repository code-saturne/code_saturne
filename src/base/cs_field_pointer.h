#ifndef __CS_FIELD_POINTER_H__
#define __CS_FIELD_POINTER_H__

/*============================================================================
 * Field pointers and ids for standard and model fields
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2013 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/* Macro used for scoping of field pointer enums */

#define CS_ENUMF_(e) CS_FIELD_POINTER_ ## e

/* Macro used to return field a field pointer by its enumerated value */

#define CS_F_(e) cs_glob_field_pointers[CS_FIELD_POINTER_ ## e].f

#define CS_FI_(e, i) cs_glob_field_pointers[CS_FIELD_POINTER_ ## e].a->p[i]

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! Enumerated field pointer ids */

typedef enum {

  CS_ENUMF_(p),            /*!< pressure */
  CS_ENUMF_(u),            /*!< velocity */

  CS_ENUMF_(rho),          /*!< density (at cells) */
  CS_ENUMF_(rho_b),        /*!< density (at boundary faces) */

  /* End of attributes */

  CS_FIELD_N_POINTERS

} cs_field_pointer_id_t;

/*! Field pointer array type */

struct cs_field_pointer_array_t {
  int           n;    /*!< number of elements */
  cs_field_t  *p[];   /*!< array of field pointers */
};

/*! Field pointer value or array type */

union cs_field_pointer_val_t {
  cs_field_t                        *f;    /*!< pointer to single field */
  struct cs_field_pointer_array_t   *a;    /*!< pointer to array of fields */
};

/*============================================================================
 * Global variables
 *============================================================================*/

/* Pointers */

extern union cs_field_pointer_val_t  *cs_glob_field_pointers;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Free all field pointer data.
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_destroy_all(void);

/*----------------------------------------------------------------------------
 * Map a simple field to an enumerated pointer.
 *
 * The associated field pointer may then be retreived using \ref CS_F_(e).
 *
 * parameters:
 *   e <--  field enumerator value
 *   f <--  pointer to field structure
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map(cs_field_pointer_id_t   e,
                     cs_field_t             *f);

/*----------------------------------------------------------------------------
 * Map a field to an (enumerated pointer, index) couple.
 *
 * This sort of mapping may be used for sets of fields whose size
 * is not known in advance.
 *
 * The associated field pointer may then be retreived using \ref CS_F_(e, i).
 *
 * parameters:
 *   e     <-- field enumerator value
 *   index <-- field enumerator index
 *   f     <-- pointer to field structure
 *----------------------------------------------------------------------------*/

void
cs_field_pointer_map_indexed(cs_field_pointer_id_t   e,
                             int                     index,
                             cs_field_t             *f);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_base(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_FIELD_POINTER_H__ */
