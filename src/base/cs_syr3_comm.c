/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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

/*============================================================================
 *  Communication with SYRTHES 3
 *============================================================================*/

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

#if defined(HAVE_SOCKET)
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
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

#define CS_SYR3_COMM_SOCKET_HEADER            "CS_comm_socket"

#define CS_SYR3_COMM_SOCKET_N_MAX            8
#define CS_LOC_SYR3_COMM_LNG_HOSTNAME      256

/*
  If SSIZE_MAX has not been defined through system headers, we take the
  minimal value required by POSIX (for low level read/write used with sockets).
*/

#if !defined(SSIZE_MAX)
#define SSIZE_MAX  32767
#endif

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

  cs_syr3_comm_type_t   type;         /* Type of data encoding */
  bool                  swap_endian;  /* Swap bytes ? */
  cs_int_t              echo;         /* Data transfer verbosity level */

#if defined(HAVE_MPI)
  MPI_Comm              intracomm;    /* Associated MPI intracommunicator */
  cs_int_t              syr_rank;     /* Syrthes root rank */
#endif

#if defined(HAVE_SOCKET)
  int                   sock;         /* Socket number */
#endif

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static char  cs_syr3_comm_elt_type_name_char[] = "c ";  /* String */
static char  cs_syr3_comm_elt_type_name_int[]  = "i ";  /* Integer */
static char  cs_syr3_comm_elt_type_name_real[] = "r8";  /* Real */

#if defined(HAVE_SOCKET)

static bool  cs_glob_comm_little_endian = false;

static char  cs_glob_comm_sock_hostname[CS_LOC_SYR3_COMM_LNG_HOSTNAME + 1];
static int   cs_glob_comm_sock_port_num = -1;

static int             cs_glob_comm_socket = 0;
struct sockaddr_in     cs_glob_comm_sock_addr;

static char  cs_glob_comm_socket_err[]
= N_("Error in socket communication:  %s (node %4d)\n");

#endif /* HAVE_SOCKET */

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

#if defined(HAVE_SOCKET)

/*----------------------------------------------------------------------------
 * Convert data from "little-endian" to "big-endian" or the reverse.
 *
 * The memory areas pointed to by src and dest should overlap either
 * exactly or not at all.
 *
 * parameters:
 *   dest <-- pointer to converted data location.
 *   src  --> pointer to source data location.
 *   size <-- size of each item of data in bytes.
 *   ni   <-- number of data items.
 *----------------------------------------------------------------------------*/

static void
_swap_endian(void        *dest,
             const void  *src,
             size_t       size,
             size_t       ni)
{
  size_t   i, ib, shift;
  unsigned char  tmpswap;

  unsigned char  *pdest = (unsigned char *)dest;
  const unsigned char  *psrc = (const unsigned char *)src;

  for (i = 0; i < ni; i++) {

    shift = i * size;

    for (ib = 0; ib < (size / 2); ib++) {

      tmpswap = *(psrc + shift + ib);
      *(pdest + shift + ib) = *(psrc + shift + (size - 1) - ib);
      *(pdest + shift + (size - 1) - ib) = tmpswap;

    }

  }

  if (dest != src && size == 1)
    memcpy(dest, src, ni);
}

/*----------------------------------------------------------------------------
 * Read a record from the interface socket
 *----------------------------------------------------------------------------*/

static void
_comm_read_sock(const cs_syr3_comm_t  *comm,
                cs_byte_t             *rec,
                const size_t           nbr,
                cs_type_t              type)
{
  size_t   start_id;
  size_t   end_id;
  size_t   n_loc;
  size_t   n_bytes;
  size_t   count = 0;
  ssize_t  ret;

  assert(rec  != NULL);
  assert(comm != NULL);

  /* Determine the number of bytes to receive */

  switch(type) {
  case CS_TYPE_cs_int_t:
    count = sizeof(cs_int_t);
    break;
  case CS_TYPE_cs_real_t:
    count = sizeof(cs_real_t);
    break;
  case CS_TYPE_char:
    count = sizeof(char);
    break;
  default:
    assert(type == CS_TYPE_cs_int_t  ||
           type == CS_TYPE_cs_real_t ||
           type == CS_TYPE_char);
  }

  n_bytes = count * nbr;

  /* Read data from socket */
  /*-----------------------*/

  start_id = 0;

  while (start_id < n_bytes) {

    end_id = CS_MIN(start_id + SSIZE_MAX, n_bytes);

    n_loc = end_id - start_id;

    ret = read(comm->sock, (void *)(rec + start_id), n_loc);

    if (ret < 1)
      bft_error(__FILE__, __LINE__, errno,
                _("Communication %s:\n"
                  "Error while receiving data by socket.\n"),
                comm->nom);

    start_id += ret;

  }

  if (comm->swap_endian == true)
    _swap_endian(rec, rec, count, nbr);

}

/*----------------------------------------------------------------------------
 * Write a record to the interface socket
 *----------------------------------------------------------------------------*/

static void
_comm_write_sock(const cs_syr3_comm_t  *comm,
                 const cs_byte_t       *rec,
                 size_t                 nbr,
                 cs_type_t              type)
{
  size_t   start_id;
  size_t   end_id;
  size_t   n_loc;
  size_t   n_bytes;
  size_t   count = 0;
  ssize_t  ret;

  cs_byte_t   * rec_tmp;

  assert(rec  != NULL);
  assert(comm != NULL);

  /* Determine the number of bytes to send */

  switch(type) {
  case CS_TYPE_cs_int_t:
    count = sizeof(cs_int_t);
    break;
  case CS_TYPE_cs_real_t:
    count = sizeof(cs_real_t);
    break;
  case CS_TYPE_char:
    count = sizeof(char);
    break;
  default:
    assert(type == CS_TYPE_cs_int_t  ||
           type == CS_TYPE_cs_real_t ||
           type == CS_TYPE_char);
  }

  n_bytes = count * nbr;

  /* Convert to "big-endian" */

  if (comm->swap_endian == true && count != 1) {
    BFT_MALLOC(rec_tmp, n_bytes, cs_byte_t);
    _swap_endian(rec_tmp, rec, count, nbr);
  }
  else
    rec_tmp = NULL;

  /* write data to socket */
  /*----------------------*/

  start_id = 0;

  while (start_id < n_bytes) {

    end_id = CS_MIN(start_id + SSIZE_MAX, n_bytes);

    n_loc = end_id - start_id;

    if (rec_tmp == NULL)
      ret = write(comm->sock, (const void *)(rec + start_id), n_loc);
    else
      ret = write(comm->sock, (const void *)(rec_tmp + start_id), n_loc);

    if (ret < 1)
      bft_error(__FILE__, __LINE__, errno,
                _("Communication %s:\n"
                  "Error sending data by socket.\n"),
                comm->nom);

    start_id += ret;

  }

  if (rec_tmp != NULL)
    BFT_FREE(rec_tmp);
}

/*----------------------------------------------------------------------------
 * Initialize socket communication
 *----------------------------------------------------------------------------*/

static void
_comm_sock_connect(cs_syr3_comm_t  *comm)
{
  int ii;

#if defined(_CS_ARCH_Linux)
  socklen_t sock_len;
#else
  int       sock_len;  /* size_t according to SUS-v2 standard, but according
                          to "man gethostbyname" under Linux, the standard
                          is bad, and we should have an int (or socklen_t) */
#endif

  char   size_str[6] = "     ";
  char  *host_names = NULL;
  int   *port_num_array = NULL;

#if defined(HAVE_MPI)
  int ierror = MPI_SUCCESS;
#endif
  int rank = (cs_glob_rank_id == -1 ? 0 : cs_glob_rank_id);

  const int lng_hostname = CS_LOC_SYR3_COMM_LNG_HOSTNAME + 1;

  /* Connect to server socket */

  sock_len = sizeof(cs_glob_comm_sock_addr);

  if (rank == 0) {

    comm->sock = accept(cs_glob_comm_socket,
                        (struct sockaddr *)&cs_glob_comm_sock_addr,
                        &sock_len);

    /* Send number of ranks */

    sprintf(size_str, "%5d", (int)cs_glob_n_ranks);

    if (write(comm->sock, size_str, 6) < 6)
      bft_error(__FILE__, __LINE__, errno,
                _("Error in socket communication\n"));
  }

  /* Obtains the name of the host machine and its port number on rank 0 */

  if (cs_glob_n_ranks > 1) {

    BFT_MALLOC(host_names,
               lng_hostname * cs_glob_n_ranks,
               char);

    BFT_MALLOC(port_num_array, cs_glob_n_ranks, int);

#if defined(HAVE_MPI)
    ierror = MPI_Gather(cs_glob_comm_sock_hostname, lng_hostname, MPI_CHAR,
                        host_names, lng_hostname, MPI_CHAR, 0,
                        cs_glob_mpi_comm);

    if (ierror < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error while sending the host name through MPI in sockets "
                  "initialization.\n"));

    /* Send the port number */

    ierror = MPI_Gather(&cs_glob_comm_sock_port_num, 1, MPI_INT,
                        port_num_array, 1, MPI_INT, 0, cs_glob_mpi_comm);

    if (ierror < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error while sending the port number through MPI in sockets "
                  "initialization.\n"));

    if (rank != 0)
      comm->sock = accept(cs_glob_comm_socket,
                          (struct sockaddr *)&cs_glob_comm_sock_addr,
                          &sock_len);

#else
    bft_error(__FILE__, __LINE__, 0,
              _("MPI is needed for socket initialization.\n"));
#endif

    /* rank 0 sends the number of ranks, hostnames, and port numbers */

    if (rank == 0) {

      /* Send max hostname size */

      sprintf(size_str, "%3d", lng_hostname);

      if (write(comm->sock, size_str, 4) < 4)
        bft_error(__FILE__, __LINE__, errno,
                  _("Error in socket communication\n"));

      for (ii = 1; ii < cs_glob_n_ranks; ii++) {

        /* Send host machine name */

        if (write(comm->sock, &(host_names[lng_hostname*ii]), lng_hostname)
            < lng_hostname)
          bft_error(__FILE__, __LINE__, errno,
                    _("Error in socket communication\n"));

        /* Send port number */

        sprintf(size_str, "%5d", port_num_array[ii]);

        if (write(comm->sock, size_str, 6) < 6)
          bft_error(__FILE__, __LINE__, errno,
                    _("Error in socket communication\n"));

      }

    } /* End of rank-0 specific operations */

    BFT_FREE(host_names);
    BFT_FREE(port_num_array);

  } /* End for cs_glob_n_ranks > 1 */

}

/*----------------------------------------------------------------------------
 * Ensure exchange of magic string through sockets
 *----------------------------------------------------------------------------*/

static void
_comm_sock_open(cs_syr3_comm_t   *comm,
                const char       *magic_string)
{
  char nom_tmp[32 + 1];

  char *magic_string_read = NULL;

  int rank = (cs_glob_rank_id == -1 ? 0 : cs_glob_rank_id);
  int len = strlen(CS_SYR3_COMM_SOCKET_HEADER);
  int magic_string_len = strlen(magic_string);

  /* Information on interface creation */

  bft_printf(_("\n  Opening communication:  %s ..."), comm->nom);
  bft_printf_flush();

  /* Initialize socket communication */

  _comm_sock_connect(comm);

  /* Check that the connection is from the correct application type */

  if (read(comm->sock, nom_tmp, len) < len)
    bft_error(__FILE__, __LINE__, errno,
              _(cs_glob_comm_socket_err), comm->nom,
              rank + 1);

  if (strncmp(nom_tmp, CS_SYR3_COMM_SOCKET_HEADER, len != 0))
    bft_error(__FILE__, __LINE__, 0,
              _("Attempt to connect to the communication port with\n"
                "an unknown message format\n"));

  /* Read magic string */
  /*-------------------*/

  BFT_MALLOC(magic_string_read, magic_string_len + 1, char);

  _comm_read_sock(comm,
                  (void *)(magic_string_read),
                  strlen(magic_string),
                  CS_TYPE_char);

  magic_string_read[magic_string_len] = '\0';

  /* If the magic string does not match, we have an error */

  if (strcmp(magic_string_read, magic_string) != 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error while initializating communication: \"%s\".\n"
                "The interface version is not correct.\n"
                "The magic string indicates the interface format version:\n"
                "magic string read:     \"%s\"\n"
                "magic string expected: \"%s\"\n"),
              comm->nom, magic_string_read, magic_string);

  BFT_FREE(magic_string_read);

  /* Write magic string */
  /*--------------------*/

  _comm_write_sock(comm,
                   (const void *)(magic_string),
                   strlen(magic_string),
                   CS_TYPE_char);

  /* Info on interface creation success */

  bft_printf(" [ok]\n");
  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Close the connection with the interface socket
 *----------------------------------------------------------------------------*/

static void
_comm_sock_close(cs_syr3_comm_t  *comm)
{
  if (close(comm->sock) != 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Communication %s):\n"
                "Error closing the socket.\n"),
              comm->nom);

  comm->sock = -1;
}

#endif /* (HAVE_SOCKET) */

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
 *   proc_rank,    <-- communicating process rank (< 0 if using sockets)
 *   mode,         <-- send or receive
 *   type,         <-- communication type
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
                        cs_syr3_comm_type_t  type,
                        cs_int_t             echo)
{
  unsigned    int_endian;

  const char base_name[] = "SYRTHES_";
  const char magic_string[] = "CFD_SYRTHES_COUPLING_2.2";
  cs_syr3_comm_t  *comm = NULL;

  BFT_MALLOC(comm, 1, cs_syr3_comm_t);

  /* Build communicator name */

  BFT_MALLOC(comm->nom, strlen(base_name) + 4 + 1, char);

  sprintf(comm->nom, "%s%04d", base_name, number);

  /* Initialize other fields */

  comm->type = type;
  comm->echo = echo;

#if defined(HAVE_MPI)
  comm->syr_rank = proc_rank;
  comm->intracomm = MPI_COMM_NULL;
#endif

#if defined(HAVE_SOCKET)
  comm->sock = -1;
#endif

  /* Test if system is "big-endian" or "little-endian" */

  comm->swap_endian = false;

  int_endian = 0;
  *((char *)(&int_endian)) = '\1';

  if (int_endian == 1)
    comm->swap_endian = true;

#if defined(DEBUG) && !defined(NDEBUG)

  else {
    int_endian = 0;
    *((char *)(&int_endian) + sizeof(unsigned) - 1) = '\1';
    assert(int_endian == 1);
  }

#endif

  /* Create interface file descriptor */
  /*----------------------------------*/

#if defined(HAVE_MPI)
  if (comm->type == CS_SYR3_COMM_TYPE_MPI)
    _comm_mpi_open(comm, magic_string);
#endif

#if defined(HAVE_SOCKET)
  if (comm->type == CS_SYR3_COMM_TYPE_SOCKET)
    _comm_sock_open(comm, magic_string);
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
  if (comm->type == CS_SYR3_COMM_TYPE_MPI)
    _comm_mpi_close(comm);
#endif

#if defined(HAVE_SOCKET)
  if (comm->type == CS_SYR3_COMM_TYPE_SOCKET)
    _comm_sock_close(comm);
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

  if (comm->type == CS_SYR3_COMM_TYPE_MPI) {

    cs_int_t  n_sec_elts_ecr = n_elts;

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

  }

#endif /* (HAVE_MPI) */

#if defined(HAVE_SOCKET)

  /* socket communication */
  /*----------------------*/

  if (comm->type == CS_SYR3_COMM_TYPE_SOCKET) {

    /* section name */

    _comm_write_sock(comm,
                     (const void *) sec_name_write,
                     CS_SYR3_COMM_H_LEN,
                     CS_TYPE_char);

    /* number of elements */

    _comm_write_sock(comm,
                     (const void *)(&n_elts),
                     1,
                     CS_TYPE_cs_int_t);

    if (n_elts != 0) {

      /* element type name */

      _comm_write_sock(comm,
                       (const void *) elt_type_name_write,
                       CS_SYR3_COMM_ELT_TYPE_NAME_LEN,
                       CS_TYPE_char);

      /* element values */

      _comm_write_sock(comm,
                       (const void *) elts,
                       (size_t) n_elts,
                       elt_type);

    }

  }

#endif /* (HAVE_SOCKET) */

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

  if (comm->type == CS_SYR3_COMM_TYPE_MPI) {

    _comm_mpi_header(header->sec_name,
                     &(header->n_elts),
                     elt_type_name,
                     CS_SYR3_COMM_MODE_RECEIVE,
                     comm);

  }

#endif /* (HAVE_MPI) */

#if defined(HAVE_SOCKET)

  /* socket communication */
  /*----------------------*/

  if (comm->type == CS_SYR3_COMM_TYPE_SOCKET) {

    /* section type name */

    _comm_read_sock(comm,
                    (void *) &(header->sec_name),
                    CS_SYR3_COMM_H_LEN,
                    CS_TYPE_char);

    /* number of elements */

    _comm_read_sock(comm,
                    (void *) &(header->n_elts),
                    1,
                    CS_TYPE_cs_int_t);


    if (header->n_elts != 0) {

      /* element type name */

      _comm_read_sock(comm,
                      (void *) elt_type_name,
                      CS_SYR3_COMM_ELT_TYPE_NAME_LEN,
                      CS_TYPE_char);

    }

  }

#endif /* (HAVE_SOCKET) */

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

    if (comm->type == CS_SYR3_COMM_TYPE_MPI)
      _comm_mpi_body((void *)_sec_elts,
                     header->n_elts,
                     header->elt_type,
                     CS_SYR3_COMM_MODE_RECEIVE,
                     comm);

#endif /* (HAVE_MPI) */

#if defined(HAVE_SOCKET)

    if (comm->type == CS_SYR3_COMM_TYPE_SOCKET)
      _comm_read_sock(comm,
                      (void *)_sec_elts,
                      (size_t) header->n_elts,
                      header->elt_type);

#endif /* (HAVE_SOCKET) */

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

#if defined(HAVE_SOCKET)

/*----------------------------------------------------------------------------
 * Open an IP socket to prepare for this communication mode
 *
 * parameters:
 *   port_num <-- port number (only used for rank 0; automatic on others)
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_init_socket(int port_num)
{
  char       s[CS_LOC_SYR3_COMM_LNG_HOSTNAME + 1];

  int        n_connect_max;

#if defined(_CS_ARCH_Linux)
  socklen_t sock_len;
#else
  int       sock_len;  /* size_t according to SUS-v2 standard, but acording
                          to "man gethostbyname" under Linux, the standard
                          is bad, and we should have an int (or socklen_t) */
#endif

  unsigned  int_endian;

  struct sockaddr_in   sock_addr;
  struct hostent      *host_ent;

  int _port_num = port_num;
  int rank = (cs_glob_rank_id == -1 ? 0 : cs_glob_rank_id);

  /* Initialization */

  n_connect_max = 0;

  if (getenv("CS_SYR3_COMM_SOCKET_N_MAX") != NULL)
    n_connect_max = atoi(getenv("CS_SYR3_COMM_SOCKET_N_MAX"));

  if (n_connect_max == 0)
    n_connect_max = CS_SYR3_COMM_SOCKET_N_MAX;

  /* Test if system is "big-endian" (network reference) or "little-endian" */

  cs_glob_comm_little_endian = false;

  int_endian = 0;
  *((char *) (&int_endian)) = '\1';

  if (int_endian == 1)
    cs_glob_comm_little_endian = true;

#if defined(DEBUG) && !defined(NDEBUG)
  else {
    int_endian = 0;
    *((char *) (&int_endian) + sizeof(unsigned) - 1) = '\1';
    assert (int_endian == 1);
  }
#endif

  /* Create server socket */

  cs_glob_comm_socket = socket(AF_INET, SOCK_STREAM, 0);

  if (cs_glob_comm_socket == -1)
    bft_error(__FILE__, __LINE__, errno,
              _("Initialization error for socket communication support.\n"));

  /* Prepare for use */

  sock_len = sizeof(sock_addr);

  memset((char *) &sock_addr, 0, sock_len);

  sock_addr.sin_family = AF_INET;
  if (rank > 0)    /* port number automatic on ranks > 0) */
    sock_addr.sin_addr.s_addr = INADDR_ANY;
  else
    sock_addr.sin_addr.s_addr = _port_num;
  sock_addr.sin_port = 0;

  if (cs_glob_comm_little_endian == true) {
    _swap_endian(&(sock_addr.sin_addr.s_addr),
                 &(sock_addr.sin_addr.s_addr),
                 sizeof(sock_addr.sin_addr.s_addr),
                 1);
    _swap_endian(&(sock_addr.sin_port),
                 &(sock_addr.sin_port),
                 sizeof(sock_addr.sin_port),
                 1);
  }

  if (gethostname(s, CS_LOC_SYR3_COMM_LNG_HOSTNAME) < 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Error obtaining computer's name"));
  s[CS_LOC_SYR3_COMM_LNG_HOSTNAME] = '\0';

  host_ent = gethostbyname(s);
  memcpy(host_ent->h_addr_list[0], &sock_addr.sin_addr, host_ent->h_length);

  if (bind(cs_glob_comm_socket,
           (struct sockaddr *)&sock_addr,
           sock_len) != 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Initialization error for socket communication support.\n"));

  if (listen(cs_glob_comm_socket, n_connect_max) < 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Initialization error for socket communication support.\n"));

  /* Obtain assigned service number */

  if (getsockname(cs_glob_comm_socket,
                  (struct sockaddr *)&sock_addr,
                  &sock_len) != 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Initialization error for socket communication support.\n"));

  _port_num = sock_addr.sin_port;
  if (cs_glob_comm_little_endian == true) {
    _swap_endian(&(sock_addr.sin_port),
                 &(sock_addr.sin_port),
                 sizeof(sock_addr.sin_port), 1);
    _port_num = sock_addr.sin_port;
    _swap_endian(&(sock_addr.sin_port),
                 &(sock_addr.sin_port),
                 sizeof(sock_addr.sin_port), 1);
  }

  /* Save the structure in the associated global variable */

  cs_glob_comm_sock_addr = sock_addr;

  /* Write host and port names in process order */

  if (rank == 0) {

    /* Print available socket information to log for rank  0
       (do not internationalize this string so that scripts
       my use it more easily). */

    bft_printf("\n  SYRTHES server port initialized\n\n");
    bft_printf_flush();

  }

  memcpy(cs_glob_comm_sock_hostname, s, CS_LOC_SYR3_COMM_LNG_HOSTNAME);
  cs_glob_comm_sock_hostname[CS_LOC_SYR3_COMM_LNG_HOSTNAME] = '\0';
  cs_glob_comm_sock_port_num = _port_num;
}

/*----------------------------------------------------------------------------
 * Close an IP socket associated with this communication mode
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_finalize_socket(void)
{
  if (cs_glob_comm_socket == 0)
    return;

  close(cs_glob_comm_socket);

  bft_printf(_("\nClosing socket...\t [ok]\n"));
  bft_printf_flush();
}

#endif /* HAVE_SOCKET */

/*----------------------------------------------------------------------------*/

END_C_DECLS
