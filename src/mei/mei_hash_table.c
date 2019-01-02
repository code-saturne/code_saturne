/*!
 * \file mei_hash_table.c
 *
 * \brief Hash table, intended to provide a symbol table
 */

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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "mei_hash_table.h"

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private functions definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Hash function.
 *
 * \param [in] s key
 * \param [in] modulo internal parameter for hash table
 *
 * \return value associated to the key s
 */
/*----------------------------------------------------------------------------*/

static unsigned
_hash(const char *const s, const int modulo)
{
  unsigned h = 0;
  int i;

  for (i = 0; s[i] != '\0'; i++) {
    h = h*256 + (unsigned) s[i];
    if (h >= (unsigned) modulo)
      h %= modulo;
  }
  return h;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Allocate a structure \em item which contain a record in the hash
 *        table.
 *
 * \param [in] key key
 *
 * \return a structure item
 */
/*----------------------------------------------------------------------------*/

static struct item*
_new_item(const char *const key)
{
  struct item* new;

  /* Allocate item proper */
  BFT_MALLOC(new, 1, struct item);

  /* Allocate memory for string */
  BFT_MALLOC(new->key, strlen(key)+1, char);

  BFT_MALLOC(new->data, 1, data_t);

  return new;
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public functions definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the hash table to the size (modulo) asked for.
 *        Allocates space for the correct number of pointers and sets them
 *        to NULL.
 *
 * \param [in] htable hash table
 * \param [in] modulo size of the hash table
 */
/*----------------------------------------------------------------------------*/

void
mei_hash_table_create(hash_table_t *const htable, const int modulo)
{
  int i;

  /* htable length and record number */
  htable->n_inter = 0;
  htable->length  = modulo;
  htable->record  = 0;
  htable->table   = NULL;

  BFT_MALLOC(htable->table, modulo, struct item*);

  /*htable default initialization: at the beginning all lists are empty */
  for (i = 0; i < modulo; i++)
    htable->table[i] = NULL;
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Find a record in a hash table.
 *
 * \param [in] htable hash table
 * \param [in] key key
 *
 * \return a pointer containing the record
 */
/*----------------------------------------------------------------------------*/

struct item*
mei_hash_table_find(hash_table_t *const htable, const char *const key)
{
  unsigned v;
  struct item* l;
  struct item* p;

  /* List where the string is to be found */
  v = _hash(key, htable->length);
  l = htable->table[v];

  /* Parcours de la liste */
  for ( p = l; p != 0; p = p->next)
    if (strcmp(p->key, key) == 0)
      return p;

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Insert a record in a hash table.
 *
 * \param [in] htable hash table
 * \param [in] key key associated to the record
 * \param [in] type flag associated to the record
 * \param [in] value store a value if the record if a real
 * \param [in] f1 pointer on a one argument function
 * \param [in] f2 pointer on a two argument function
 */
/*----------------------------------------------------------------------------*/

void
mei_hash_table_insert(hash_table_t *const htable,
                      const char *const  key,
                      const mei_flag_t type,
                      const double value,
                      const func1_t f1,
                      const func2_t f2)
{
  unsigned v;
  struct item* item = mei_hash_table_find(htable, key);

  if (item == NULL) {

    /* The string was not found. Create a new item to insert it */
    item = _new_item(key);

    item->type = type;
    if (type == FUNC1) {
      item->data->func = f1;

    } else if (type == FUNC2) {
      item->data->f2 = f2;

    } else if (type == FUNC3) {
      bft_error(__FILE__, __LINE__, 0, "not implemented yet \n");

    } else if (type == FUNC4) {
      bft_error(__FILE__, __LINE__, 0, "not implemented yet \n");

    } else {
      item->data->value = value;
    }
    strcpy(item->key, key);

    htable->record++;

    /* Insert a cell at the head of a list */
    v = _hash(key, htable->length);
    item->next = htable->table[v];
    htable->table[v]= item;
  }
  else {
    item->data->value = value;
  }
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief Find a record in a hash table.
 *
 * \param [in] htable hash table
 * \param [in] key key
 *
 * \return a pointer containing the record
 */
/*----------------------------------------------------------------------------*/


struct item*
mei_hash_table_lookup(hash_table_t *const htable, const char *const key)
{
  unsigned v;
  struct item *item, *next;   /* Pointers to current and next record
                               * while traversing hash table bucket.  */

  v = _hash(key, htable->length);

  /* Lookup for name in hash table record corresponding to above hash value. */

  item = htable->table[v];

  while (item) {
    next = item->next;
    if (!strcmp(item->key, key))
      return item;
    item = next;
  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a hash table.
 *
 * \param [in] htable hash table
 */
/*----------------------------------------------------------------------------*/


void
mei_hash_table_free(hash_table_t *const htable)
{
  int i;
  struct item *item, *next;   /* Pointers to current and next record
                               * while traversing hash table bucket.  */

  if (htable == NULL) return;

  /* Delete hash table as well as data structure representing symbol table. */
  for (i = 0; i < htable->length; i++) {

    item = htable->table[i];

    while (item) {
      next = item->next;
      BFT_FREE(item->key);
      BFT_FREE(item->data);
      BFT_FREE(item);
      item = next;
    }
  }
  BFT_FREE(htable->table);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialize the hash table with default symbols
 *
 * \param [in] htable hash table
 */
/*----------------------------------------------------------------------------*/


void
mei_hash_table_init(hash_table_t *const htable)
{
  int i, j;

  /* predefined constants names*/
  static const char *constants_names[] = { "e", "pi" };

  /* predefined constants values */
  static const double constants[] = {2.7182818284590452354,
                                     3.14159265358979323846};

  /* predefined functions names with one arguments*/
  static const char *functions_names[] = { "exp",  "log",   "sqrt",
                                           "sin",  "cos",   "tan",
                                           "asin", "acos",  "atan",
                                           "sinh", "cosh",  "tanh",
                                           "abs" , "int" };

  /* predefined functions pointers to functions to calculate them */
  static double (*functions[]) (double) = { exp,  log,  sqrt,
                                            sin,  cos,  tan,
                                            asin, acos, atan,
                                            sinh, cosh, tanh,
                                            fabs, floor };

  /* predefined functions names with two arguments*/
  static const char *functions2_names[] = { "atan2", "min", "max", "mod" };

  /* predefined functions pointers to functions to calculate them */
  static double (*functions2[]) (double, double) = { atan2, fmin, fmax, fmod };

  j = sizeof(constants_names) / sizeof(constants_names[0]);
  for (i = 0; i < j; i++)
    mei_hash_table_insert(htable,
                          constants_names[i],
                          CONSTANT,
                          constants[i],
                          NULL,
                          NULL);

  j = sizeof(functions_names) / sizeof(functions_names[0]);
  for (i = 0; i < j; i++)
    mei_hash_table_insert(htable,
                          functions_names[i],
                          FUNC1,
                          0,
                          functions[i],
                          NULL);

  j = sizeof(functions2_names) / sizeof(functions2_names[0]);
  for (i = 0; i < j; i++)
    mei_hash_table_insert(htable,
                          functions2_names[i],
                          FUNC2,
                          0,
                          NULL,
                          functions2[i]);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump function of a single record.
 *
 * \param [in] item record
 */
/*----------------------------------------------------------------------------*/

void
mei_hash_table_item_print(struct item *item)
{
  /* Affichage d'une liste */
  while (item != NULL) {
    printf("%s -> %i \n", item->key, item->type);
    if (item->type != FUNC1 &&
        item->type != FUNC2 &&
        item->type != FUNC3 &&
        item->type != FUNC4)
      printf("valeur : %f\n", item->data->value);
    item = item->next;
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump of table contents for debuging purpose.
 *
 * \param [in] htable hash table
 */
/*----------------------------------------------------------------------------*/

void
mei_hash_table_dump(hash_table_t *const htable)
{
  int i;
  for (i = 0; i < htable->length; i++) {
    if (htable->table[i] != NULL) {
      printf("Entry %d \n", i);
      mei_hash_table_item_print(htable->table[i]);
    }
  }
}

/*============================================================================
 * Test
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

#undef _TEST_
#ifdef _TEST_

int
main_test(int argc, char *argv[])
{
#define HASHSIZE 701 // Modulo by default for hash tables (101, 211, 701,...)
  hash_table_t htable;
  struct item *item;
  char *toto = "pi";
  printf("Hash Table Testing\n");

  mei_hash_table_create(&htable, HASHSIZE);
  mei_hash_table_init(&htable);
  mei_hash_table_dump(&htable);
  item = mei_hash_table_lookup(&htable, toto);
  if (item) mei_hash_table_item_print(item);
  mei_hash_table_free(&htable);

  return 0;
}

#endif

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

