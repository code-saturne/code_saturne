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

#ifndef _SYR_COMM_H_
#define _SYR_COMM_H_

/*============================================================================*
 * Definitions of base communication functions
 *
 * Library: Code_Saturne                               Copyright EDF 2006-2009
 *============================================================================*/

/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "syr_defs.h"

/*----------------------------------------------------------------------------
 * Macro definitions
 *----------------------------------------------------------------------------*/

#define SYR_COMM_L_SEC_NAME  32

/*----------------------------------------------------------------------------
 * Structure definitions
 *----------------------------------------------------------------------------*/

typedef struct _syr_comm_t syr_comm_t;

/* Public structure used to save data from a message header, simplifying
   the transfer of this data to processing functions */

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a communicator
 *
 * arguments:
 *   coupling_num  <-- Coupling number
 *   cs_root_rank  <-- Root rank associated with Code_Saturne
 *   cs_n_ranks    <-- Number of MPI ranks associated with Code_Saturne
 *   echo          <-- Echo on main output
 *----------------------------------------------------------------------------*/

syr_comm_t *
syr_comm_initialize(int               coupling_num,
                    int               cs_root_rank,
                    int               cs_n_ranks,
                    int               echo);

/*----------------------------------------------------------------------------
 * Finalize a communicator
 *----------------------------------------------------------------------------*/

syr_comm_t *
syr_comm_finalize(syr_comm_t *comm);

/*----------------------------------------------------------------------------
 * Get a communicator's name
 *
 * This function returns a pointer to an existing name, so the string
 * returned should not be freed by the user.
 *----------------------------------------------------------------------------*/

const char *
syr_comm_get_name(const syr_comm_t *comm);

/*----------------------------------------------------------------------------
 * Get a communicator's number of distant processes
 *----------------------------------------------------------------------------*/

int
syr_comm_get_n_procs(const syr_comm_t *comm);

/*----------------------------------------------------------------------------
 * Get the number of values to be read for each sub-domain in a communicator.
 *
 * This function should only be called between syr_comm_read_header()
 * and syr_comm_read_data(). It returns a pointer to an array belonging to
 * the communicator, which should not be freed.
 *----------------------------------------------------------------------------*/

const int *
syr_comm_get_n_section_elements(const syr_comm_t *comm);

/*----------------------------------------------------------------------------
 * Write a section to the communication interface.
 *
 * A section contains:
 *
 * 1. section name
 * 2. number of elements
 * 3. element type
 * 4. element values
 *----------------------------------------------------------------------------*/

void
syr_comm_write_section(const char        *sec_name,
                       int                n_elts,
                       void              *elts,
                       syr_type_t         elt_type,
                       const syr_comm_t  *comm,
                       int                proc_id);

/*----------------------------------------------------------------------------
 * Read a section header from the communication interface.
 *
 * A section contains:
 *
 * 1. section name (size: SYR_COMM_L_SEC_NAME + 1)
 * 2. number of elements
 * 3. element type
 * 4. element values
 *
 * The first 3 of which constitute its header.
 *----------------------------------------------------------------------------*/

void
syr_comm_read_header(char              *sec_name,
                     int               *n_elts,
                     syr_type_t        *elt_type,
                     const syr_comm_t  *comm,
                     int                proc_id);

/*----------------------------------------------------------------------------
 * Read a section body from the communication interface.
 *
 * A section contains:
 *
 * 1. section name
 * 2. number of elements
 * 3. element type
 * 4. element values
 *
 * The last of which constitutes its body.
 *
 * If a buffer destined to receive data already exists, we give a pointer to
 * this buffer with the sec_elts argument, and this same pointer is returned.
 * Otherwise (setting the argument to NULL), memory is allocated, and a pointer
 * to this newly allocated memory is returned (and should be freed by the
 * caller once finished).
 *----------------------------------------------------------------------------*/

void *
syr_comm_read_body(int                n_sec_elts,
                   void              *sec_elts,
                   syr_type_t         elt_type,
                   const syr_comm_t  *comm,
                   int                proc_id);

/*----------------------------------------------------------------------------*/

#endif /* _SYR_COMM_H_ */
