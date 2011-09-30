/*============================================================================
 *  Communication with SYRTHES 3
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2011 EDF S.A.

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

/*
  Define _GNU_SOURCE if necessary before including any headers, to ensure
  the correct feature macros are defined first.
*/

#if defined(__linux__)
# define _GNU_SOURCE 1
#endif

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_syr3_comm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define CS_SYR3_COMM_ELT_TYPE_NAME_LEN       2    /* Length of type name */

/*----------------------------------------------------------------------------
 * Send or receive a message
 *----------------------------------------------------------------------------*/

typedef enum {

  CS_SYR3_COMM_MODE_RECEIVE,  /* Receive  */
  CS_SYR3_COMM_MODE_SEND      /* Send */

} _cs_syr3_comm_mode_t;

/*============================================================================
 * Local Structure Definitions
 *============================================================================*/

struct _cs_syr3_comm_t {

  char                 *nom;          /* Communicator name */

  int                   echo;         /* Data transfer verbosity level */

#if defined(HAVE_MPI)
  MPI_Comm              intracomm;    /* Associated MPI intracommunicator */
  int                   syr_rank;     /* Syrthes root rank */
#endif

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static char  cs_syr3_comm_elt_type_name_char[] = "c ";  /* String */
static char  cs_syr3_comm_elt_type_name_int[]  = "i ";  /* Integer */
static char  cs_syr3_comm_elt_type_name_real[] = "r8";  /* Real */

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Print an error message in case of MPI communication problem
 *----------------------------------------------------------------------------*/

static void
_comm_mpi_error_msg(const cs_syr3_comm_t  *comm,
                    int                    error)
{
  char buffer[MPI_MAX_ERROR_STRING];
  int  buffer_len;

  MPI_Error_string(error, buffer, &buffer_len);

  bft_error(__FILE__, __LINE__, 0,
            _("MPI error for communication:  %s\n"
              "Error type: %s"), comm->nom, buffer);
}

/*----------------------------------------------------------------------------
 * Initialize an MPI communication by receiving then sending a "magic string"
 * to check the correct data format.
 *----------------------------------------------------------------------------*/

static void
_comm_mpi_open(cs_syr3_comm_t  *comm,
               const char      *magic_string)
{
  int ierror;

  MPI_Status status;

  char * comm_magic_string;

  cs_int_t magic_string_len = strlen(magic_string);

  /* Create intracommunicator */
  /*--------------------------*/

  int local_range[2] = {-1, -1};
  int distant_range[2] = {-1, -1};

  bft_printf(_(" Initializing MPI communication: %s ... "), comm->nom);
  bft_printf_flush();

  ple_coupling_mpi_intracomm_create(MPI_COMM_WORLD,
                                    cs_glob_mpi_comm,
                                    comm->syr_rank,
                                    &(comm->intracomm),
                                    local_range,
                                    distant_range);

  bft_printf(_("[ok]\n"));
  bft_printf(_("  Local ranks = [%d..%d], distant ranks = [%d..%d].\n\n"),
             local_range[0], local_range[1] - 1,
             distant_range[0], distant_range[1] - 1);
  bft_printf_flush();

  comm->syr_rank = distant_range[0];

  /* Receive magic string */
  /*----------------------*/

  BFT_MALLOC(comm_magic_string, magic_string_len + 1, char);

  ierror = MPI_Recv(comm_magic_string, magic_string_len, MPI_CHAR,
                    comm->syr_rank,
                    MPI_ANY_TAG, comm->intracomm, &status);

  if (ierror != MPI_SUCCESS)
    _comm_mpi_error_msg(comm, ierror);

  comm_magic_string[magic_string_len] = '\0';

  /* If the magic string does not match, we have an error */

  if (strcmp(comm_magic_string, magic_string) != 0)
    bft_error
      (__FILE__, __LINE__, 0,
       _("Error for communication: \"%s\".\n"
         "The interface version is not correct.\n"
         "The magic string indicates an incorrect interface format version.\n"
         "magic string read:     \"%s\"\n"
         "magic string expected: \"%s\"\n"),
       comm->nom, comm_magic_string, magic_string);

  /* Send magic string */
  /*-------------------*/

  strncpy(comm_magic_string, magic_string, magic_string_len);

  ierror = MPI_Send(comm_magic_string, magic_string_len, MPI_CHAR,
                    comm->syr_rank,
                    0, comm->intracomm);

  if (ierror != MPI_SUCCESS)
    _comm_mpi_error_msg(comm, ierror);

  BFT_FREE(comm_magic_string);
}

/*----------------------------------------------------------------------------
 * Free MPI intracommunicator
 *----------------------------------------------------------------------------*/

static void
_comm_mpi_close(cs_syr3_comm_t  *comm)
{
  if (comm->intracomm != MPI_COMM_NULL) {
    MPI_Comm_free(&(comm->intracomm));
    comm->intracomm = MPI_COMM_NULL;
  }
}

/*----------------------------------------------------------------------------
 * Exchange a section header through MPI
 *----------------------------------------------------------------------------*/

static void
_comm_mpi_header(char                  *sec_name,
                 cs_int_t              *n_sec_elts,
                 char                  *elt_type_name,
                 _cs_syr3_comm_mode_t   mode,
                 const cs_syr3_comm_t  *comm)
{
#undef  CS_SYR3_COMM_MPI_PACK_SIZE
#define CS_SYR3_COMM_MPI_PACK_SIZE    CS_SYR3_COMM_H_LEN \
                                    + CS_SYR3_COMM_ELT_TYPE_NAME_LEN \
                                    + (sizeof(int) * 2)

  char buffer[CS_SYR3_COMM_MPI_PACK_SIZE];

  int position, ierror;

  MPI_Status  status;

  /* Instructions */

  assert(comm != NULL);
  assert(*n_sec_elts >= 0);
  assert(sizeof(int) == sizeof(cs_int_t));

  /* Receive mode */
  /*--------------*/

  if (mode == CS_SYR3_COMM_MODE_RECEIVE) {

    /* Receive message */

    ierror = MPI_Recv(buffer, CS_SYR3_COMM_MPI_PACK_SIZE, MPI_PACKED,
                      comm->syr_rank,
                      MPI_ANY_TAG, comm->intracomm, &status);

    if (ierror != MPI_SUCCESS)
      _comm_mpi_error_msg(comm, ierror);

    /* Extract buffer elements */

    position = 0;

    MPI_Unpack(buffer, CS_SYR3_COMM_MPI_PACK_SIZE, &position, sec_name,
               CS_SYR3_COMM_H_LEN, MPI_CHAR, comm->intracomm);

    MPI_Unpack(buffer, CS_SYR3_COMM_MPI_PACK_SIZE, &position, n_sec_elts,
               1, CS_MPI_INT, comm->intracomm);

    if (*n_sec_elts > 0)
      MPI_Unpack(buffer, CS_SYR3_COMM_MPI_PACK_SIZE, &position, elt_type_name,
                 CS_SYR3_COMM_ELT_TYPE_NAME_LEN, MPI_CHAR, comm->intracomm);

  }

  /* Send mode */
  /*-----------*/

  else { /* if (mode == CS_SYR3_COMM_MODE_SEND) */

    /* Pack buffer */

    position = 0;

    MPI_Pack(sec_name, CS_SYR3_COMM_H_LEN, MPI_CHAR, buffer,
             CS_SYR3_COMM_MPI_PACK_SIZE, &position, comm->intracomm);

    MPI_Pack(n_sec_elts, 1, CS_MPI_INT, buffer, CS_SYR3_COMM_MPI_PACK_SIZE,
             &position, comm->intracomm);

    if (*n_sec_elts > 0)
      MPI_Pack(elt_type_name, CS_SYR3_COMM_ELT_TYPE_NAME_LEN, MPI_CHAR, buffer,
               CS_SYR3_COMM_MPI_PACK_SIZE, &position, comm->intracomm);

    /* Send message */

    ierror = MPI_Send(buffer, position, MPI_PACKED, comm->syr_rank,
                      0, comm->intracomm);

    if (ierror != MPI_SUCCESS)
      _comm_mpi_error_msg(comm, ierror);
  }
}

/*----------------------------------------------------------------------------
 * Exchange a section body through MPI
 *----------------------------------------------------------------------------*/

static void
_comm_mpi_body(void                  *sec_elts,
               cs_int_t               n_sec_elts,
               cs_type_t              elt_type,
               _cs_syr3_comm_mode_t   mode,
               const cs_syr3_comm_t  *comm)
{
  int ierror = 0;
  int n_elts = n_sec_elts;

  MPI_Status  status;

  /* Instructions */

  assert(comm != NULL);
  assert(n_sec_elts >= 0);

  /* Receive mode */
  /*--------------*/

  if (mode == CS_SYR3_COMM_MODE_RECEIVE) {

    switch (elt_type) {

    case CS_TYPE_cs_int_t:

      ierror = MPI_Recv(sec_elts, n_elts, CS_MPI_INT,
                        comm->syr_rank,
                        MPI_ANY_TAG, comm->intracomm, &status);
      break;

    case CS_TYPE_cs_real_t:
      ierror = MPI_Recv(sec_elts, n_elts, CS_MPI_REAL,
                        comm->syr_rank,
                        MPI_ANY_TAG, comm->intracomm, &status);
      break;

    case CS_TYPE_char:
      ierror = MPI_Recv(sec_elts, n_elts, MPI_CHAR,
                        comm->syr_rank,
                        MPI_ANY_TAG, comm->intracomm, &status);
      break;

    default:
      assert (   elt_type == CS_TYPE_char
              || elt_type == CS_TYPE_cs_int_t
              || elt_type == CS_TYPE_cs_real_t);
    }

  }

  /* Send mode */
  /*-----------*/

  else { /* if (mode == CS_SYR3_COMM_MODE_SEND) */

    switch (elt_type) {

    case CS_TYPE_cs_int_t:
      ierror = MPI_Send(sec_elts, n_elts, CS_MPI_INT,
                        comm->syr_rank,
                        0, comm->intracomm);
      break;

    case CS_TYPE_cs_real_t:
      ierror = MPI_Send(sec_elts, n_elts, CS_MPI_REAL,
                        comm->syr_rank,
                        0, comm->intracomm);
      break;

    case CS_TYPE_char:
      ierror = MPI_Send(sec_elts, n_elts, MPI_CHAR,
                        comm->syr_rank,
                        0, comm->intracomm);
      break;

    default:
      assert(   elt_type == CS_TYPE_char
             || elt_type == CS_TYPE_cs_int_t
             || elt_type == CS_TYPE_cs_real_t);
    }

  }

  if (ierror != MPI_SUCCESS)
    _comm_mpi_error_msg(comm, ierror);
}

#endif /* (HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Print information on waiting for a message
 *----------------------------------------------------------------------------*/

static void
_comm_echo_pre(const cs_syr3_comm_t  *comm,
               _cs_syr3_comm_mode_t   mode)
{
  assert(comm != NULL);

  if (mode == CS_SYR3_COMM_MODE_RECEIVE)
    bft_printf(_("\nMessage received on \"%s\":\n"), comm->nom);

  else /* if (mode == CS_SYR3_COMM_MODE_SEND) */
    bft_printf(_("\nMessage sent on \"%s\":\n"), comm->nom);

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Print a message header
 *----------------------------------------------------------------------------*/

static void
_comm_echo_header(const char  *sec_name,
                  cs_int_t     n_elts,
                  cs_type_t    elt_type
)
{
  char sec_name_write[CS_SYR3_COMM_H_LEN + 1];

  /* instructions */

  strncpy(sec_name_write, sec_name,  CS_SYR3_COMM_H_LEN);
  sec_name_write[CS_SYR3_COMM_H_LEN] = '\0';

  bft_printf(_("    section name:          \"%s\"\n"
               "    number of elements:    %d\n"),
             sec_name_write, n_elts);

  if (n_elts > 0) {

    char *nom_typ = NULL;

    switch(elt_type) {
    case CS_TYPE_char:
      nom_typ = cs_syr3_comm_elt_type_name_char;
      break;
    case CS_TYPE_cs_int_t:
      nom_typ = cs_syr3_comm_elt_type_name_int;
      break;
    case CS_TYPE_cs_real_t:
      nom_typ = cs_syr3_comm_elt_type_name_real;
      break;
    default:
      assert(   elt_type == CS_TYPE_char
             || elt_type == CS_TYPE_cs_int_t
             || elt_type == CS_TYPE_cs_real_t);
    }

    bft_printf(_("    element type name:      \"%s\"\n"), nom_typ);

  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Print (partial) message content
 *----------------------------------------------------------------------------*/

static void
_comm_echo_body(cs_int_t     echo,
                cs_int_t     n_elts,
                cs_type_t    elt_type,
                const void  *sec_elts
)
{
  cs_int_t  echo_start = 0;
  cs_int_t  echo_end;
  cs_int_t  ii;

  /* Instructions */

  if (n_elts == 0) return;

  if (echo * 2 < n_elts) {
    echo_end = echo;
    bft_printf(_("    %d first and last elements:\n"), echo);
  }
  else {
    echo_end = n_elts;
    bft_printf(_("    elements:\n"));
  }

  do {

    switch (elt_type) {

    case CS_TYPE_cs_int_t:
      {
        const cs_int_t *sec_elts_int = (const cs_int_t *) sec_elts;

        for (ii = echo_start ; ii < echo_end ; ii++)
          bft_printf("    %10d : %12d\n", ii + 1, *(sec_elts_int + ii));
      }
      break;

    case CS_TYPE_cs_real_t:
      {
        const cs_real_t *sec_elts_real = (const cs_real_t *) sec_elts;

        for (ii = echo_start ; ii < echo_end ; ii++)
          bft_printf("    %10d : %12.5e\n", ii + 1, *(sec_elts_real + ii));
      }
      break;

    case CS_TYPE_char:
      {
        const char *sec_elts_char = (const char *) sec_elts;

        for (ii = echo_start ; ii < echo_end ; ii++) {
          if (*(sec_elts_char + ii) != '\0')
            bft_printf("    %10d : '%c'\n", ii + 1, *(sec_elts_char + ii));
          else
            bft_printf("    %10d : '\\0'\n", ii + 1);
        }
      }
      break;

    default:

      assert(   elt_type == CS_TYPE_cs_int_t
             || elt_type == CS_TYPE_cs_real_t
             || elt_type == CS_TYPE_char);

    }

    if (echo_end < n_elts) {
      bft_printf("    ..........   ............\n");
      echo_start = n_elts - echo;
      echo_end = n_elts;
    }
    else {
      assert(echo_end == n_elts);
      echo_end = n_elts + 1;
    }

  } while (echo_end <= n_elts);

  bft_printf_flush();

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a communication
 *
 * parameters:
 *   number,       <-- coupling number
 *   proc_rank,    <-- communicating process rank
 *   mode,         <-- send or receive
 *   echo          <-- echo on main output (< 0 if none, header if 0,
 *                     n first and last elements if n)
 *
 * returns:
 *   pointer to communication structure
 *----------------------------------------------------------------------------*/

cs_syr3_comm_t *
cs_syr3_comm_initialize(int                  number,
#if defined(HAVE_MPI)
                        int                  proc_rank,
#endif
                        cs_int_t             echo)
{
  const char base_name[] = "SYRTHES_";
  const char magic_string[] = "CFD_SYRTHES_COUPLING_2.2";
  cs_syr3_comm_t  *comm = NULL;

  BFT_MALLOC(comm, 1, cs_syr3_comm_t);

  /* Build communicator name */

  BFT_MALLOC(comm->nom, strlen(base_name) + 4 + 1, char);

  sprintf(comm->nom, "%s%04d", base_name, number);

  /* Initialize other fields */

  comm->echo = echo;

#if defined(HAVE_MPI)
  comm->syr_rank = proc_rank;
  comm->intracomm = MPI_COMM_NULL;
#endif

  /* Create interface file descriptor */
  /*----------------------------------*/

#if defined(HAVE_MPI)
  _comm_mpi_open(comm, magic_string);
#endif

  return comm;
}

/*----------------------------------------------------------------------------
 * Function finalizing a communication
 *----------------------------------------------------------------------------*/

cs_syr3_comm_t *
cs_syr3_comm_finalize(cs_syr3_comm_t *comm)
{
  /* Info on interface finalization */

  bft_printf(_("\n  Closing communication:  %s\n"), comm->nom);
  bft_printf_flush();

#if defined(HAVE_MPI)
  _comm_mpi_close(comm);
#endif

  BFT_FREE(comm->nom);
  BFT_FREE(comm);

  return NULL;
}

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
cs_syr3_comm_get_name(const cs_syr3_comm_t  *comm)
{
  assert(comm != NULL);

  return(comm->nom);
}

/*----------------------------------------------------------------------------
 * Send message
 *
 * parameters:
 *   sec_name  <-- section name
 *   n_elts    <-- number of elements
 *   elt_type  <-- element type if n_elts > 0
 *   elts      <-- elements if n_elts > 0
 *   comm      <-- communicator
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_send_message(const char             sec_name[CS_SYR3_COMM_H_LEN],
                          cs_int_t               n_elts,
                          cs_type_t              elt_type,
                          void                  *elts,
                          const cs_syr3_comm_t  *comm)
{
  cs_int_t  n_sec_elts_ecr;
  char   sec_name_write[CS_SYR3_COMM_H_LEN + 1];
  char   elt_type_name_write[CS_SYR3_COMM_ELT_TYPE_NAME_LEN + 1];

  char  *elt_type_name = NULL;


  assert(comm != NULL);
  assert(n_elts >= 0);

  /* section name */

  sprintf(sec_name_write,
          "%-*.*s",
          CS_SYR3_COMM_H_LEN,
          CS_SYR3_COMM_H_LEN,
          sec_name);

  /* element type name */

  if (n_elts != 0) {

    switch(elt_type) {

    case CS_TYPE_cs_int_t:
      elt_type_name = cs_syr3_comm_elt_type_name_int;
      break;

    case CS_TYPE_cs_real_t:
      elt_type_name = cs_syr3_comm_elt_type_name_real;
      break;

    case CS_TYPE_char:
      elt_type_name = cs_syr3_comm_elt_type_name_char;
      break;

    default:
      assert(   elt_type == CS_TYPE_cs_int_t
             || elt_type == CS_TYPE_cs_real_t
             || elt_type == CS_TYPE_char);

    }

    sprintf(elt_type_name_write,
            "%-*.*s",
            CS_SYR3_COMM_ELT_TYPE_NAME_LEN,
            CS_SYR3_COMM_ELT_TYPE_NAME_LEN,
            elt_type_name);

  }

  if (comm->echo  >= 0)
    _comm_echo_pre(comm, CS_SYR3_COMM_MODE_SEND);


#if defined(HAVE_MPI)

  /* MPI communication */
  /*-------------------*/

  n_sec_elts_ecr = n_elts;

  _comm_mpi_header(sec_name_write,
                   &n_sec_elts_ecr,
                   elt_type_name_write,
                   CS_SYR3_COMM_MODE_SEND,
                   comm);

  if (n_elts > 0)
    _comm_mpi_body((void *) elts,
                   n_elts,
                   elt_type,
                   CS_SYR3_COMM_MODE_SEND,
                   comm);

#endif /* (HAVE_MPI) */

  /* Possibly print to log file */

  if (comm->echo  >= 0)
    _comm_echo_header(sec_name,
                      n_elts,
                      elt_type);

  if (comm->echo > 0)
    _comm_echo_body(comm->echo,
                    n_elts,
                    elt_type,
                    elts);
}

/*----------------------------------------------------------------------------
 * Receive message header
 *
 * parameters:
 *   header --> message header
 *   comm   <-- communicator
 *
 * returns
 *   number of elements in message body
 *----------------------------------------------------------------------------*/

cs_int_t
cs_syr3_comm_receive_header(cs_syr3_comm_msg_header_t  *header,
                            const cs_syr3_comm_t       *comm)
{
  char   elt_type_name[CS_SYR3_COMM_ELT_TYPE_NAME_LEN + 1];

  assert(comm  != NULL);

  header->n_elts = 0;

  if (comm->echo >= 0)
    _comm_echo_pre(comm, CS_SYR3_COMM_MODE_RECEIVE);


#if defined(HAVE_MPI)

  /* MPI communication */
  /*-------------------*/

  _comm_mpi_header(header->sec_name,
                   &(header->n_elts),
                   elt_type_name,
                   CS_SYR3_COMM_MODE_RECEIVE,
                   comm);

#endif /* (HAVE_MPI) */

  header->sec_name[CS_SYR3_COMM_H_LEN] = '\0';

  if (header->n_elts != 0) {

    elt_type_name[CS_SYR3_COMM_ELT_TYPE_NAME_LEN] = '\0';

    if (strcmp(elt_type_name, cs_syr3_comm_elt_type_name_int) == 0)
      header->elt_type = CS_TYPE_cs_int_t;

    else if (strcmp(elt_type_name, cs_syr3_comm_elt_type_name_real) == 0)
      header->elt_type = CS_TYPE_cs_real_t;

    else if (strcmp(elt_type_name, cs_syr3_comm_elt_type_name_char) == 0)
      header->elt_type = CS_TYPE_char;

    else {
      assert(0);
    }
  }

  /* Possibly print to log file */

  if (comm->echo >= 0)
    _comm_echo_header(header->sec_name,
                      header->n_elts,
                      header->elt_type);

  /* Return number of elements to read */

  return header->n_elts;
}

/*----------------------------------------------------------------------------
 * Receive a message body
 *
 * parameters:
 *   header <-- message header
 *   elt    --> received body values
 *   comm   <-- communicator
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_receive_body(const cs_syr3_comm_msg_header_t  *header,
                          void                             *elt,
                          const cs_syr3_comm_t             *comm)
{
  cs_int_t    ii;
  void      *_sec_elts;

  assert(comm  != NULL);
  assert(header->n_elts >= 0);

  _sec_elts = elt;

  if (_sec_elts == NULL && header->n_elts != 0) {

    switch(header->elt_type) {

    case CS_TYPE_cs_int_t:
      {
        cs_int_t  *sec_elts_int;

        BFT_MALLOC(sec_elts_int, header->n_elts, cs_int_t);
        _sec_elts = (void *) sec_elts_int;
      }
      break;

    case CS_TYPE_cs_real_t:
      {
        cs_real_t  *sec_elts_rea;

        BFT_MALLOC(sec_elts_rea, header->n_elts, cs_real_t);
        _sec_elts = (void *)sec_elts_rea;
      }
      break;

    case CS_TYPE_char:
      {
        char  *sec_elts_cha;

        BFT_MALLOC(sec_elts_cha, header->n_elts + 1, char);
        _sec_elts = (void *)sec_elts_cha;
      }
      break;

    default:
      assert(   header->elt_type == CS_TYPE_cs_int_t
             || header->elt_type == CS_TYPE_cs_real_t
             || header->elt_type == CS_TYPE_char);
    }

  }

  /* element values */

  if (header->n_elts != 0) {

#if defined(HAVE_MPI)

    _comm_mpi_body((void *)_sec_elts,
                   header->n_elts,
                   header->elt_type,
                   CS_SYR3_COMM_MODE_RECEIVE,
                   comm);

#endif /* (HAVE_MPI) */

    /* Verification */

    if (header->elt_type == CS_TYPE_char) {
      for (ii = 0 ;
           ii < header->n_elts && ((char *)_sec_elts)[ii] != '\0' ;
           ii++);
      ((char *)_sec_elts)[ii] = '\0';
    }

    /* Possibly print to log file */

    if (comm->echo > 0)
      _comm_echo_body(comm->echo,
                      header->n_elts,
                      header->elt_type,
                      _sec_elts);

  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
