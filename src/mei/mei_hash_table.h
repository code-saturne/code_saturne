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

/*============================================================================
 * Enum definition
 *============================================================================*/

/*!
 * \brief List of the different type of symbol
 */

typedef enum {

  CONSTANT,
  ID,
  FUNC1,
  FUNC2,
  FUNC3,
  FUNC4,
  OPR

} mei_flag_t;

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
  func2_t f2;      /*!< Pointer to function with two arguments */
} data_t;

/*!
 * \brief Type definition for each record of the hash table
 */

struct item {
  char        *key;  /*!< Pointeur to string */
  mei_flag_t   type; /*!< Constant, variable, function,... */
  data_t      *data; /*!< Data of the current record */
  struct item *next; /*!< Pointer to next element */
};

/*!
 * \brief Structure defining a hash table
 */

struct HashTable {
  int           n_inter; /*!< number of interpreters associated with
                           the current table of symbols */
  int           record;  /*!< number of records of the hash table*/
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
 * Initialize the hash table to the size (modulo) asked for.
 * Allocates space for the correct number of pointers and sets them to NULL.
 *
 * param [in] htable hash table
 * param [in] modulo size of the hash table
 *----------------------------------------------------------------------------*/

void mei_hash_table_create(hash_table_t *const htable,
                           const int modulo);

/*----------------------------------------------------------------------------*/
/*
 * Initialize the hash table with default symbols
 *
 * param [in] htable hash table
 */
/*----------------------------------------------------------------------------*/

void mei_hash_table_init(hash_table_t *htable);

/*----------------------------------------------------------------------------*/
/*
 * Destroy a hash table.
 *
 * param [in] htable hash table
 */
/*----------------------------------------------------------------------------*/

void mei_hash_table_free(hash_table_t *htable);

/*----------------------------------------------------------------------------*/
/*
 * Find a record in a hash table.
 *
 * param [in] htable hash table
 * param [in] key key
 *
 * return a pointer containing the record
 */
/*----------------------------------------------------------------------------*/

struct item * mei_hash_table_lookup(hash_table_t *htable,
                                    const char *key);

/*----------------------------------------------------------------------------*/
/*
 * Insert a record in a hash table.
 *
 * param [in] htable hash table
 * param [in] key key associated to the record
 * param [in] type flag associated to the record
 * param [in] value store a value if the record if a real
 * param [in] f1 pointer on a one argument function
 * param [in] f2 pointer on a two argument function
 */
/*----------------------------------------------------------------------------*/

void
mei_hash_table_insert(hash_table_t *const htable,
                      const char *const  key,
                      const mei_flag_t type,
                      const double value,
                      const func1_t f1,
                      const func2_t f2);

/*----------------------------------------------------------------------------*/
/*
 * Find a record in a hash table.
 *
 * param [in] htable hash table
 * param [in] key key
 *
 * return a pointer containing the record
 */
/*----------------------------------------------------------------------------*/

struct item* mei_hash_table_find(hash_table_t *htable,
                                 const char *key);

/*----------------------------------------------------------------------------*/
/*
 * Dump of table contents for debuging purpose.
 *
 * param [in] htable hash table
 */
/*----------------------------------------------------------------------------*/

void mei_hash_table_dump(hash_table_t *htable);

/*----------------------------------------------------------------------------*/
/*
 * Dump function of a single record.
 *
 * param [in] item record
 */
/*----------------------------------------------------------------------------*/

void mei_hash_table_item_print(struct item *item);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MEI_HASH_TABLE_H__ */

