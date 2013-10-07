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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "cs_field.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_field_pointer.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_field_pointer.c
        Field pointers and ids for standard and model fields.
*/

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Number of pointers (initially fixed, but extensible in case
   user fields should be added after the model fields) */

static cs_field_pointer_id_t _n_pointers = 0;
static union cs_field_pointer_val_t  *_field_pointer = NULL;

/* Handling of sublists */

static bool  *_is_sublist = NULL;

/*============================================================================
 * Global variables
 *============================================================================*/

/* Pointers */

union cs_field_pointer_val_t  *cs_glob_field_pointers = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize field pointers
 *----------------------------------------------------------------------------*/

static void
_init_pointers(void)
{
  assert(_field_pointer == NULL);

  _n_pointers = CS_FIELD_N_POINTERS;
  BFT_MALLOC(_field_pointer, _n_pointers, union cs_field_pointer_val_t);
  BFT_MALLOC(_is_sublist, _n_pointers, bool);

  for (cs_field_pointer_id_t i = 0; i < _n_pointers; i++) {
    _field_pointer[i].f = NULL;
    _is_sublist[i] = false;
  }

  cs_glob_field_pointers = _field_pointer;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all field pointer data.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_destroy_all(void)
{
  for (cs_field_pointer_id_t i = 0; i < _n_pointers; i++) {
    if (_is_sublist[i])
      BFT_FREE(_field_pointer[i].a);
  }
  BFT_FREE(_field_pointer);
  BFT_FREE(_is_sublist);

  cs_glob_field_pointers = NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map a simple field to an enumerated pointer.
 *
 * The associated field pointer may then be retreived using \ref CS_F_(e).
 *
 * \param[in]  e   field enumerator value
 * \param[in]  f   pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map(cs_field_pointer_id_t   e,
                     cs_field_t             *f)
{
  if (_field_pointer == NULL)
    _init_pointers();

  assert(e < _n_pointers);

 union cs_field_pointer_val_t v;
  v.f = f;
  _field_pointer[e] = v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map a field to an (enumerated pointer, index) couple.
 *
 * This sort of mapping may be used for sets of fields whose size
 * is not known in advance.
 *
 * The associated field pointer may then be retreived using \ref CS_F_(e, i).
 *
 * \param[in]  e      field enumerator value
 * \param[in]  index  field enumerator index
 * \param[in]  f      pointer to field structure
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_indexed(cs_field_pointer_id_t   e,
                             int                     index,
                             cs_field_t             *f)
{
  assert(index >= 0);

  if (_field_pointer == NULL)
    _init_pointers();

  struct cs_field_pointer_array_t *a;
  union cs_field_pointer_val_t v;

  int i;
  int _sub_size = index + 1;
  int _sub_size_prev = 0;

  /* Check for previous size */

  assert(e < _n_pointers);

  v = _field_pointer[e];

  if (v.f != NULL) {

    /* Check also that we did not already use in incompatible mapping */
    if (! _is_sublist[e]) {
      cs_field_t *_f = v.f;
      bft_error(__FILE__, __LINE__, 0,
                _("%s: field enum %d is already mapped as non-indexed\n"
                  "to field id %d (%s), so it cannot be mapped as indexed."),
                __func__, (int)e, _f->id, _f->name);
    }

    else {

      a = v.a;
      _sub_size_prev = a->n;

      if (_sub_size_prev < _sub_size) {
        /* BFT_MALLOC does not directly handle C flexible array members,
           (which we use here to minimize type width), so use bytes */
        void *p = a;
        size_t _s =   sizeof(struct cs_field_pointer_array_t)
                    + sizeof(cs_field_t *) * _sub_size;
        BFT_REALLOC(p, _s, unsigned char);
        a = p;
        v.a = a;
        for (i = _sub_size_prev; i < index; i++)
          a->p[i] = NULL;
      }

      _is_sublist[e] = true;

    }

  }

  v.a->p[index] = f;
  _field_pointer[e] = v;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Map base fields to enumerated pointers.
 */
/*----------------------------------------------------------------------------*/

void
cs_field_pointer_map_base(void)
{
  cs_field_pointer_map(CS_ENUMF_(p),
                       cs_field_by_name("pressure"));
  cs_field_pointer_map(CS_ENUMF_(u),
                       cs_field_by_name("velocity"));

  cs_field_pointer_map(CS_ENUMF_(rho),
                       cs_field_by_name("density"));
  cs_field_pointer_map(CS_ENUMF_(rho_b),
                       cs_field_by_name("boundary_density"));
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
