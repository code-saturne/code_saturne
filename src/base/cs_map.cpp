/*============================================================================
 * Map helper structure
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft/bft_error.h"
#include "bft/bft_printf.h"

#include "base/cs_base.h"
#include "base/cs_mem.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "base/cs_map.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Align to 64-bit size for better performance */

#define CS_ALIGN_SIZE(s) (((s-1)/8+1)*8)

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/* Basic Map structure */
/*---------------------*/

struct _cs_map_name_to_id_t {

  size_t      size;             /* Number of entries */
  size_t      max_size;         /* Maximum number of entries */

  unsigned    n_key_blocks;
  char      **key_blocks;       /* Key buffer (pointers to available space
                                   and end-of buffer occupy first bytes).
                                   Each block is at least twice as large
                                   as the previous one */

  char      **key;              /* Pointer to keys */
  int        *id;               /* Matching id */
  int        *reverse_id;       /* Reverse id mapping */
};

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Insert new key.
 *
 * parameters:
 *   m     <-> pointer to map structure
 *   key   <-- character string (key)
 *   id    <-- id associated with key
 *   index <-- index of key in map
 *----------------------------------------------------------------------------*/

static void
_name_to_id_insert_key(cs_map_name_to_id_t  *m,
                       const char           *key,
                       int                   id,
                       size_t                index)
{
  size_t i;
  size_t key_size = CS_ALIGN_SIZE(strlen(key) + 1);

  /* Resize map arrays if necessary */

  if (m->size >= m->max_size) {

    size_t prev_size = m->max_size;

    m->max_size*= 2;
    CS_REALLOC(m->key, m->max_size, char *);
    CS_REALLOC(m->id, m->max_size, int);
    CS_REALLOC(m->reverse_id, m->max_size, int);

    for (i = prev_size; i < m->max_size; i++) {
      m->key[i] = nullptr;
      m->id[i] = -1;
      m->reverse_id[i] = -1;
    }
  }

  /* Check for buffer space */

  assert(m->n_key_blocks > 0);

  char *key_block = m->key_blocks[m->n_key_blocks - 1];

  char *s = ((char **)key_block)[0];
  char *e = ((char **)key_block)[1];

  if (s + key_size > e) {

    /* The new block should be twice as large as the current one
       (to limit reallocations and fragmenting), and be able to contain at
       least several keys of similar size as the current one. */
    const size_t n_keys_mini = 8;
    size_t b_size = e - key_block;
    while (b_size < (key_size * n_keys_mini) + 2*sizeof(char*))
      b_size *= 2;

    CS_MALLOC(key_block, b_size, char);
    memset(key_block, 0, b_size);

    CS_REALLOC(m->key_blocks, m->n_key_blocks+1, char *);
    m->key_blocks[m->n_key_blocks] = key_block;
    m->n_key_blocks++;

    s = key_block + CS_ALIGN_SIZE(2*sizeof(char *));
    e = key_block + b_size;

    ((char **)key_block)[0] = s;
    ((char **)key_block)[1] = e;

  }

  /* Shift previous data */

  for (i = m->size; i > index; i--) {
    m->key[i] = m->key[i-1];
    m->id[i] = m->id[i-1];
    m->reverse_id[m->id[i]] = i;
  }

  /* Insert data */

  strcpy(s, key);

  m->key[index] = s;
  m->id[index] = id;
  m->reverse_id[m->size] = index;

  s += key_size;
  ((char **)key_block)[0] = s;

  m->size += 1;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create empty name to id map.
 *
 * returns:
 *   pointer to newly initialized map structure.
 *----------------------------------------------------------------------------*/

cs_map_name_to_id_t *
cs_map_name_to_id_create(void)
{
  cs_map_name_to_id_t *m = nullptr;

  CS_MALLOC(m, 1, cs_map_name_to_id_t);

  m->size = 0;
  m->max_size = 8;

  m->n_key_blocks = 1;
  CS_MALLOC(m->key_blocks, 1, char *);

  const size_t b_size = 128;
  char *key_block;

  CS_MALLOC(key_block, b_size, char);
  memset(key_block, 0, b_size);

  char *s = key_block + CS_ALIGN_SIZE(2*sizeof(char *));
  char *e = key_block + b_size;

  ((char **)key_block)[0] = s;
  ((char **)key_block)[1] = e;

  m->key_blocks[0] = key_block;
  m->n_key_blocks = 1;

  CS_MALLOC(m->key, m->max_size, char *);
  CS_MALLOC(m->id, m->max_size, int);
  CS_MALLOC(m->reverse_id, m->max_size, int);

  return m;
}

/*----------------------------------------------------------------------------
 * Destroy name to id map structure.
 *
 * parameters:
 *   m <-> pointer to map structure.
 *----------------------------------------------------------------------------*/

void
cs_map_name_to_id_destroy(cs_map_name_to_id_t **m)
{
  if (m != nullptr) {

    if (*m != nullptr) {

      cs_map_name_to_id_t *_m = *m;

      CS_FREE(_m->reverse_id);
      CS_FREE(_m->id);
      CS_FREE(_m->key);

      for (size_t i = 0; i < _m->n_key_blocks; i++)
        CS_FREE(_m->key_blocks[i]);
      CS_FREE(_m->key_blocks);

      CS_FREE(*m);

    }
  }
}

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
                  const char           *key)
{
  int start_id, end_id, mid_id;
  int cmp_ret = 1;

  /* Use binary search to find entry */

  start_id = 0;
  end_id = m->size - 1;
  mid_id = start_id + ((end_id -start_id) / 2);

  while (start_id <= end_id) {
    cmp_ret = strcmp(m->key[mid_id], key);
    if (cmp_ret < 0)
      start_id = mid_id + 1;
    else if (cmp_ret > 0)
      end_id = mid_id - 1;
    else
      break;
    mid_id = start_id + ((end_id -start_id) / 2);
  }

  /* If not found, insert key */

  if (cmp_ret != 0)
    _name_to_id_insert_key(m, key, m->size, mid_id);

  return m->id[mid_id];
}

/*----------------------------------------------------------------------------
 * Return id matching a key, or -1 if not present.
 *
 *
 * parameters:
 *   m      <-- pointer to map structure
 *   key    <-- character string (key)
 *
 * returns:
 *   id matching key, or -1.
 *----------------------------------------------------------------------------*/

int
cs_map_name_to_id_try(const cs_map_name_to_id_t  *m,
                      const char                 *key)
{
  int start_id, end_id, mid_id;
  int cmp_ret = 1;

  int retval = -1;

  /* Use binary search to find entry */

  if (m != nullptr) {

    start_id = 0;
    end_id = m->size - 1;
    mid_id = start_id + ((end_id -start_id) / 2);

    while (start_id <= end_id) {
      cmp_ret = strcmp(m->key[mid_id], key);
      if (cmp_ret < 0)
        start_id = mid_id + 1;
      else if (cmp_ret > 0)
        end_id = mid_id - 1;
      else
        break;
      mid_id = start_id + ((end_id -start_id) / 2);
    }

    if (cmp_ret == 0)
      retval = m->id[mid_id];

  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Return a key name in a map matching a given id.
 *
 * parameters:
 *   m  <-- pointer to map structure.
 *   id <-- key id
 *
 * returns:
 *   pointer to key.
 *----------------------------------------------------------------------------*/

const char *
cs_map_name_to_id_reverse(const cs_map_name_to_id_t  *m,
                          size_t                      id)
{
  const char *retval = nullptr;

  if (m == nullptr)
    return retval;

  if (id < m->size)
    retval = (const char *)(m->key[m->reverse_id[id]]);

  return retval;
}

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
cs_map_name_to_id_size(const cs_map_name_to_id_t *m)
{
  size_t retval = 0;

  if (m != nullptr)
    retval = m->size;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return key in a map for a given index position.
 *
 * parameters:
 *   m     <-- pointer to map structure.
 *   index <-- key index
 *
 * returns:
 *   pointer to key.
 *----------------------------------------------------------------------------*/

const char *
cs_map_name_to_id_key(const cs_map_name_to_id_t  *m,
                      size_t                      index)
{
  const char *retval = nullptr;

  if (m == nullptr)
    return retval;

  if (index < m->size)
    retval = (const char *)(m->key[index]);

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
