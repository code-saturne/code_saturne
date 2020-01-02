#ifndef __CS_NUMBERING_H__
#define __CS_NUMBERING_H__

/*============================================================================
 * Numbering information for vectorization or multithreading
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
#include "cs_base.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* SIMD unit length */
/*------------------*/

#if defined(__NEC__) && defined(__ve__)  /* For NEC Aurora series */

#  define CS_NUMBERING_SIMD_SIZE 256

#elif defined(SX) && defined(_SX)        /* For NEC SX series (at least SX-9) */

#  define CS_NUMBERING_SIMD_SIZE 256

#elif defined(__AVX512F__)               /* For Intel with AVX 512 */

#  define CS_NUMBERING_SIMD_SIZE 64

#else

#  define CS_NUMBERING_SIMD_SIZE 4       /* Most current platforms */

#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Renumbering types */

typedef enum {

  CS_NUMBERING_DEFAULT,    /* Default numbering */
  CS_NUMBERING_VECTORIZE,  /* Numbered for vectorization */
  CS_NUMBERING_THREADS     /* Numbered for threads */

} cs_numbering_type_t;

/* Renumbering structure */

typedef struct {

  cs_numbering_type_t   type;     /* Numbering type */

  int   vector_size;              /* Vector size if vectorized, 1 otherwise */

  int   n_threads;                /* Number of threads */
  int   n_groups;                 /* Number of associated groups */

  int   n_no_adj_halo_groups;     /* number of groups for which only non-ghost
                                     entities are adjacent */

  cs_lnum_t  n_no_adj_halo_elts;  /* Number of elements not adjacent to
                                     halo elements */

  cs_lnum_t *group_index;         /* For thread t and group g, the start and
                                     past-the-end ids for entities in a given
                                     group and thread are respectively:
                                     group_index[t*n_groups*2 + g] and
                                     group_index[t*n_groups*2 + g + 1].
                                     (size: n_groups * n_threads * 2) */

} cs_numbering_t;

/*=============================================================================
 * Global variable definitions
 *============================================================================*/

/* Names for numbering types */

extern const char  *cs_numbering_type_name[];

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a default numbering information structure.
 *
 * \param[in]  n_elts  number of associated elements
 *
 * \return  pointer to created cs_numbering_t structure
 */
/*----------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_default(cs_lnum_t  n_elts);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a default numbering information structure
 *        in case of vectorization.
 *
 * \param[in]  n_elts       number of associated elements
 * \param[in]  vector_size  vector size used for this vectorization
 *
 * \return  pointer to created cs_numbering_t structure
 */
/*----------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_vectorized(cs_lnum_t  n_elts,
                               int        vector_size);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create a default numbering information structure
 *        in case of threading.
 *
 * \param[in]  n_groups     number of groups
 * \param[in]  group_index  group_index[thread_id*group_id*2 + group_id*2] and
 *                          group_index[thread_id*group_id*2 + group_id*2 +1]
 *                          define the start and end ids for entities in a
 *                          given group and thread;
 *                          (size: n_groups *2 * n_threads)
 *
 * \return  pointer to created cs_numbering_t structure
 */
/*----------------------------------------------------------------------------*/

cs_numbering_t *
cs_numbering_create_threaded(int        n_threads,
                             int        n_groups,
                             cs_lnum_t  group_index[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Destroy a numbering information structure.
 *
 * \param[in, out]  numbering  pointer to cs_numbering_t structure pointer
 *                             (or NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_numbering_destroy(cs_numbering_t  **numbering);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log information relative to a cs_numbering_t structure.
 *
 * \param[in]  log          log type
 * \param[in]  description  description of numbering type
 * \param[in]  numbering    pointer to cs_numbering_t structure (or NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_numbering_log_info(cs_log_t               log,
                      const char            *description,
                      const cs_numbering_t  *numbering);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump a cs_numbering_t structure.
 *
 * \param[in]  numbering    pointer to cs_numbering_t structure (or NULL)
 */
/*----------------------------------------------------------------------------*/

void
cs_numbering_dump(const cs_numbering_t  *numbering);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_NUMBERING_H__ */
