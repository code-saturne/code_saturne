#ifndef __CS_CONTROL_H__
#define __CS_CONTROL_H__

/*============================================================================
 * Interactive control management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
#include "cs_time_step.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

typedef enum {

  CS_CONTROL_COMM_TYPE_SOCKET,    /* Communicate through sockets */
  CS_CONTROL_COMM_TYPE_NULL       /* Null communicator */

} cs_control_comm_type_t;

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize controller structures.
 */
/*----------------------------------------------------------------------------*/

void
cs_control_finalize(void);

/*----------------------------------------------------------------------------
 * Check the presence of a control file and deal with the interactive
 * control.
 *----------------------------------------------------------------------------*/

void
cs_control_check_file(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Establish a connection to a client.
 *
 * \param[in]  port_name  name of server port (host:port for IP sockets)
 * \param[in]  key        key for authentification
 * \param[in]  type       communication type
 */
/*----------------------------------------------------------------------------*/

void
cs_control_comm_initialize(const char              *port_name,
                           const char              *key,
                           cs_control_comm_type_t   type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize a connection to a client.
 */
/*----------------------------------------------------------------------------*/

void
cs_control_comm_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write a record to a client.
 *
 * \param[in]  rec    pointer to data to write
 * \param[in]  size   size of each data element, in bytes
 * \param[in]  count  number of data elements
 */
/*----------------------------------------------------------------------------*/

void
cs_control_comm_write(const void  *rec,
                      size_t       size,
                      size_t       count);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read data from a client into a command queue
 *
 * The function updates a pointer (view) to the data.
 *
 * \return number of useable elements read
 *         (i.e. elements before the last separator)
 */
/*----------------------------------------------------------------------------*/

size_t
cs_control_comm_read_to_queue(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CONTROL_H__ */
