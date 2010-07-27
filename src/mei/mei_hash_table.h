#ifndef __MEI_HASH_TABLE_H__
#define __MEI_HASH_TABLE_H__

/*!
 * \file mei_hash_table.h
 *
 * \brief Hash table, intended to provide a symbol table
 *
 * A hash table consists of an array of container. Each container
 * holds a copy of the key, a pointer to the data associated with the
 * key, and a pointer to the next container that associated with this one,
 * if there was one.
 */

/*
  This file is part of the "Mathematical Expression Interpreter" library.

  Copyright (C) 2008-2009  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Enum definition
 *============================================================================*/

/*!
 * \brief List of the different type of symbol
 */

typedef enum { CONSTANT, ID, FUNC1, FUNC2, FUNC3, FUNC4, OPR } mei_flag_t;

/*============================================================================
 * Type definition
 *============================================================================*/

/*!
 * \brief Type definition for a pointer to a function of one argument
 */

typedef double (*func1_t) (double);

/*!
 * \brief Type definition for pointer to a function of two arguments
 */

typedef double (*func2_t) (double, double);

/*!
 * \brief Type definition for pointer to a function of three arguments
 */

typedef double (*func3_t) (double, double, double);

/*!
 * \brief Type definition for pointer to a function of for arguments
 */

typedef double (*func4_t) (double, double, double, double);

/*!
 * \brief Type definition for data of each element contained in the hash table
 */

typedef union {
  double value;    /*!< Constant or variable value.  */
  func1_t func;    /*!< Pointer to function with one argument */
  func2_t f2;      /*!< Pointer to function with two argument */
} data_t;

/*!
 * \brief Type definition for each bucket of the hash table
 */

struct item {
  char        *key;  /*!< Pointeur to string */
  mei_flag_t   type; /*!< Constant, variable, function,... */
  data_t      *data; /*!< Data of the current bucket */
  struct item *next; /*!< Pointer to next element */
};

/*!
 * \brief Structure defining a hash table
 */

struct HashTable {
  int           n_inter; /*!< number of interpreters associated with
                           the current table of symbols */
  int           record;  /*!< number of buckets of the hash table*/
  int           length;  /*!< length of the hash table */
  struct item **table;   /*!< 'table' is a list of pointers on 'item' */
};

/*!
 * \brief Type definition for a hash table
 */

typedef struct HashTable hash_table_t;


/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui cree une table de hachage
 *----------------------------------------------------------------------------*/


void mei_hash_table_create(hash_table_t *const htable, const int modulo);


/*----------------------------------------------------------------------------
 *  Fonction qui initialise la table de hachage
 *----------------------------------------------------------------------------*/


void mei_hash_table_init(hash_table_t *htable);


void mei_hash_table_dump(hash_table_t *htable);


void mei_hash_table_item_print(struct item *item);


void mei_hash_table_free(hash_table_t *htable);


struct item * mei_hash_table_lookup(hash_table_t *htable, const char *key);


void mei_hash_table_insert(hash_table_t *const htable,
                           const char *const  key,
                           const mei_flag_t type,
                           const double value,
                           const func1_t func,
                           const func2_t f2,
                           const func3_t f3,
                           const func4_t f4);


struct item* mei_hash_table_find(hash_table_t *htable, const char *key);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MEI_HASH_TABLE_H__ */

