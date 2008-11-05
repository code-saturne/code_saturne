/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_SYR3_COMM_H__
#define __CS_SYR3_COMM_H__

/*============================================================================
 * Communication with SYRTHES 3
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

#define CS_SYR3_COMM_FIN_FICHIER                           "EOF"

#define CS_SYR3_COMM_H_LEN       32   /* Length of a header name */

/* Socket communications: we suppose a maximum of 8 coupled SYRTHES instances;
   this value may be modified through the CS_SYR3_COMM_SOCKET_NBR_MAX
   environment variable */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Message type
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_SYR3_COMM_TYPE_NONE,     /* No communication (pre-initialization) */
  CS_SYR3_COMM_TYPE_MPI,      /* MPI messages */
  CS_SYR3_COMM_TYPE_SOCKET    /* IP sockets */

} cs_syr3_comm_type_t;

/*----------------------------------------------------------------------------
 * Send or receive a message
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_SYR3_COMM_MODE_RECEPTION,   /* Receive  */
  CS_SYR3_COMM_MODE_EMISSION     /* Send */

} cs_syr3_comm_mode_t;

/* Pointer associated with an opaque communicator structure. */

typedef struct _cs_syr3_comm_t cs_syr3_comm_t;

/* Structure used to save message header data, to simplify its use. */

typedef struct {

  char       nom_rub[CS_SYR3_COMM_H_LEN + 1];
  cs_int_t   nbr_elt;
  cs_type_t  typ_elt;

} cs_syr3_comm_msg_entete_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Function initializing a communication
 *
 * parameters:
 *   numero,       <-- coupling number
 *   rang_proc,    <-- communicating process rank (< 0 if using sockets)
 *   mode,         <-- send or receive
 *   type,         <-- communication type
 *   echo          <-- echo on main output (< 0 if none, header if 0,
 *                     n first and last elements if n)
 *
 * returns:
 *   pointer to communication structure
 *----------------------------------------------------------------------------*/

cs_syr3_comm_t *
cs_syr3_comm_initialise(const cs_int_t             numero,
#if defined(_CS_HAVE_MPI)
                        const cs_int_t             rang_proc,
#endif
                        const cs_syr3_comm_mode_t  mode,
                        const cs_syr3_comm_type_t  type,
                        const cs_int_t             echo);

/*----------------------------------------------------------------------------
 * Function finalizing a communication
 *----------------------------------------------------------------------------*/

cs_syr3_comm_t *
cs_syr3_comm_termine(cs_syr3_comm_t *comm);

/*----------------------------------------------------------------------------
 * Return a pointer to a communicator name
 *
 * parameters:
 *   comm <-- communicator
 *
 * returns:
 *   pointer to communicator name
 *----------------------------------------------------------------------------*/

const char *
cs_syr3_comm_ret_nom(const cs_syr3_comm_t  *comm);

/*----------------------------------------------------------------------------
 * Send message
 *
 * parameters:
 *   nom_rub <-- section name
 *   nbr_elt <-- number of elemeents
 *   typ_elt <-- element type if nbr_elt > 0
 *   elt     <-- elements if nbr_elt > 0
 *   comm    <-- communicator
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_envoie_message(const char             nom_rub[CS_SYR3_COMM_H_LEN],
                            cs_int_t               nbr_elt,
                            cs_type_t              typ_elt,
                            void                  *elt,
                            const cs_syr3_comm_t  *comm);

/*----------------------------------------------------------------------------
 * Receive message header
 *
 * parameters:
 *   entete --> message header
 *   comm   <-- communicator
 *
 * returns
 *   number of elements in message body
 *----------------------------------------------------------------------------*/

cs_int_t
cs_syr3_comm_recoit_entete(cs_syr3_comm_msg_entete_t  *entete,
                           const cs_syr3_comm_t       *comm);

/*----------------------------------------------------------------------------
 * Receive a message body
 *
 * parameters:
 *   entete <-- message header
 *   elt    --> received body values
 *   comm   <-- communicator
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_recoit_corps(const cs_syr3_comm_msg_entete_t  *entete,
                          void                             *elt,
                          const cs_syr3_comm_t             *comm);

#if defined(_CS_HAVE_SOCKET)

/*----------------------------------------------------------------------------
 * Open an IP socket to prepare for this communication mode
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_init_socket(void);

/*----------------------------------------------------------------------------
 * Close an IP socket associated with this communication mode
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_termine_socket(void);

#endif /* _CS_HAVE_SOCKET */

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_SYR3_COMM_H__ */
