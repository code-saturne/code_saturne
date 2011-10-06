/*===========================================================================
 * Definitions of base communication functions
 *===========================================================================*/

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

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

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

#include <ple_defs.h>
#if defined(HAVE_MPI)
#include <ple_coupling.h>
#endif

/*---------------------------------------------------------------------------
 *  Local headers
 *---------------------------------------------------------------------------*/

#include "syr_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "syr_comm.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*===========================================================================
 *  Constants and Macros
 *===========================================================================*/

#define CS_COMM_SOCKET_HEADER       "CS_comm_socket"

#define CS_COMM_SOCKET_N_MAX         65
#define CS_COMM_L_HOSTNAME          256
#define CS_COMM_L_NOM_MAX           256

#define SYR_COMM_L_TYPE_NAME         2

/* If SSIZE_MAX is not defined by the sytem headers, we take the minimum value
   required by POSIX (for low level reads/writes with sockets). */

#if !defined(SSIZE_MAX)
# define SSIZE_MAX  32767
#endif

typedef enum {

  SYR_COMM_MODE_RECEIVE,   /* Receive */
  SYR_COMM_MODE_SEND       /* Send */

} syr_comm_mode_t;

/*===========================================================================
 * Structure definition
 *===========================================================================*/

struct _syr_comm_t {

  char            *name;         /* Communicator name */

  int              swap_endian;  /* Force big-endian communications */

  syr_comm_type_t  type;         /* Communicator type */
  int              n_procs;      /* Number of communicating processes */
  int              echo;         /* Log (printout) level of communications */

  int             *n_sec_elts;   /* Number of elements in a section for each
                                    proc (when reading) */

#if defined(HAVE_SOCKET)
  int             *socket;       /* Array of socket numbers */
#endif

#if defined(HAVE_MPI)
  int              proc_rank;    /* Rank of first distant process */
  MPI_Comm         intracomm;    /* Intracommunicator */
#endif

};

/*===========================================================================
 * Global variables
 *===========================================================================*/

#if defined(HAVE_SOCKET)

static char  syr_glob_comm_err_socket[]
  = "Error in socket communication:\n"
    "%s (proc %4d)\n";

#endif

/*===========================================================================
 * Private function definitions
 *===========================================================================*/

/*--------------------------------------------------------------------------
 * Convert data from "little-endian" to "big-endian" or the reverse.
 *
 * The memory areas pointed to by src and dest should overlap either
 * exactly or not at all.
 *--------------------------------------------------------------------------*/

static void
_swap_endian(void          *dest,
             const void    *src,
             const size_t   size,
             const size_t   ni)
{
  size_t   i, ib, shift;
  unsigned char  tmpswap;

  unsigned char  *pdest = (unsigned char *)dest;
  const unsigned char  *psrc = (const unsigned char *)src;

  for (i = 0 ; i < ni ; i++) {

    shift = i * size;

    for (ib = 0 ; ib < (size / 2) ; ib++) {

      tmpswap = *(psrc + shift + ib);
      *(pdest + shift + ib) = *(psrc + shift + (size - 1) - ib);
      *(pdest + shift + (size - 1) - ib) = tmpswap;

    }

  }

  if (dest != src && size == 1)
    memcpy(dest, src, ni);
}

/*--------------------------------------------------------------------------
 * Print "wait on message"
 *--------------------------------------------------------------------------*/

static void
_comm_echo_pre(const syr_comm_t  *comm,
               int                proc_id,
               syr_comm_mode_t    mode)
{
  assert(comm != NULL);

  switch(mode) {

  case SYR_COMM_MODE_RECEIVE:
    if (comm->n_procs == 1)
      printf("\nMessage recu sur \"%s\":\n", comm->name);
    else
      printf("\nMessage recu sur \"%s\" (proc %d):\n",
             comm->name, proc_id + 1);

    break;

  case SYR_COMM_MODE_SEND:
    if (comm->n_procs == 1)
      printf("\nMessage envoye sur \"%s\":\n", comm->name);
    else
      printf("\nMessage envoye sur \"%s\" (proc %d):\n",
             comm->name, proc_id + 1);

    break;

  default:
    assert (mode == SYR_COMM_MODE_RECEIVE ||
            mode == SYR_COMM_MODE_SEND    );

  }

  fflush(stdout);
}

/*--------------------------------------------------------------------------
 * Print a section header
 *
 *  1. section name
 *  2. number of elements
 *  3. element type name
 *--------------------------------------------------------------------------*/

static void
_comm_echo_header(const syr_comm_t  *comm,
                  const char        *sec_name,
                  int                n_elts,
                  const char        *type_name)
{
  char sec_name_out[SYR_COMM_L_SEC_NAME + 1];

  assert(comm != NULL);

  strncpy(sec_name_out, sec_name, SYR_COMM_L_SEC_NAME);
  sec_name_out[SYR_COMM_L_SEC_NAME] = '\0';

  printf("    nom de rubrique :    \"%s\"\n"
         "    nombre de valeurs :  %d\n",
         sec_name_out, n_elts);

  if (n_elts > 0) {

    char type_name_wri[SYR_COMM_L_TYPE_NAME + 1];

    strncpy(type_name_wri, type_name, SYR_COMM_L_TYPE_NAME);
    type_name_wri[SYR_COMM_L_TYPE_NAME] = '\0';

    printf("    type de valeurs :    \"%s\"\n", type_name_wri);

  }

  fflush(stdout);

}

/*--------------------------------------------------------------------------
 * Print (part of) a section's body
 *--------------------------------------------------------------------------*/

static void
_comm_echo_body(const syr_comm_t  *comm,
                int                n_elts,
                syr_type_t         type,
                void              *sec_elts)
{
  int        start_id;
  int        end_id;
  int        id;

  if (n_elts == 0) return;

  assert(comm != NULL);

  start_id = 0;

  if (comm->echo * 2 < n_elts) {

    end_id = comm->echo;
    printf("    %d premieres et dernieres valeurs :\n", comm->echo);

  }
  else {

    end_id = n_elts;
    printf("    valeurs :\n");

  }

  do {

    switch(type) {

    case SYR_TYPE_int: {

      int *sec_elts_int = (int *)sec_elts;

      for (id = start_id; id < end_id; id++)
        printf("    %10d : %12d\n", id + 1, *(sec_elts_int + id));

    }
    break;

    case SYR_TYPE_float: {

      float *sec_elts_float = (float *)sec_elts;

      for (id = start_id; id < end_id; id++)
        printf("    %10d : %12.5e\n", id + 1,
               (double)(*(sec_elts_float + id)));

    }
      break;

    case SYR_TYPE_double: {

      double *sec_elts_double = (double *)sec_elts;

      for (id = start_id; id < end_id; id++)
        printf("    %10d : %14.7e\n", id + 1, *(sec_elts_double + id));

    }
    break;

    case SYR_TYPE_char: {

      char *sec_elts_char = (char *) sec_elts;

      for (id = start_id; id < end_id; id++) {
        if (*(sec_elts_char + id) != '\0')
          printf("    %10d : '%c'\n", id + 1, *(sec_elts_char + id));
        else
          printf("    %10d : '\\0'\n", id + 1);
      }

    }
    break;

    default:

      assert(type == SYR_TYPE_int    ||
             type == SYR_TYPE_double ||
             type == SYR_TYPE_float  ||
             type == SYR_TYPE_char    );

    } /* End of switch on element type */

    if (end_id < n_elts) {

      printf("    ..........   ............\n");

      start_id = n_elts - comm->echo;
      end_id = n_elts;

    }
    else {

      assert(end_id == n_elts);
      end_id = n_elts + 1;

    }

  } while (end_id <= n_elts);

  fflush(stdout);
}

#if defined(HAVE_MPI)

/*--------------------------------------------------------------------------
 *  Print message in case of MPI communication error
 *--------------------------------------------------------------------------*/

static void
_comm_mpi_msg_err(const syr_comm_t  *comm,
                  int                proc_id,
                  int                error)
{
  char buffer[MPI_MAX_ERROR_STRING];
  int  buffer_len;

  MPI_Error_string(error, buffer, &buffer_len);

  ple_error(__FILE__, __LINE__, 0,
            "Erreur MPI pour la communication : %s (proc %4d)\n"
            "Type d'erreur : %s", comm->name, proc_id + 1, buffer);
}

/*--------------------------------------------------------------------------
 * Exchange a section header using MPI
 *
 * A section header contains:
 *
 * 1. section name
 * 2. number of elements
 * 3. element type
 *--------------------------------------------------------------------------*/

static void
_comm_mpi_header(char              *sec_name,
                 int               *n_sec_elts,
                 char              *type_name,
                 const syr_comm_t  *comm,
                 syr_comm_mode_t    mode,
                 int                proc_id)
{

#undef  SYR_COMM_MPI_PACK_SIZE
#define SYR_COMM_MPI_PACK_SIZE    SYR_COMM_L_SEC_NAME \
                                + SYR_COMM_L_TYPE_NAME \
                                + (sizeof(int) * 2)

  char buffer[SYR_COMM_MPI_PACK_SIZE];

  int position, ierror;

  MPI_Status  status;

  assert(comm != NULL);

  /* Receive mode */
  /*------------- */

  if (mode == SYR_COMM_MODE_RECEIVE) {

    /* Receive message */

    ierror = MPI_Recv(buffer,
                      SYR_COMM_MPI_PACK_SIZE,
                      MPI_PACKED,
                      comm->proc_rank + proc_id,
                      MPI_ANY_TAG,
                      comm->intracomm,
                      &status);

    if (ierror != MPI_SUCCESS)
      _comm_mpi_msg_err(comm,
                        proc_id,
                        ierror);

    /* Extract elements from buffer */

    position = 0;

    MPI_Unpack(buffer, SYR_COMM_MPI_PACK_SIZE, &position, sec_name,
               SYR_COMM_L_SEC_NAME, MPI_CHAR, comm->intracomm);

    MPI_Unpack(buffer, SYR_COMM_MPI_PACK_SIZE, &position, n_sec_elts,
               1, MPI_INT, comm->intracomm);

    if (*n_sec_elts > 0)
      MPI_Unpack(buffer, SYR_COMM_MPI_PACK_SIZE, &position, type_name,
                 SYR_COMM_L_TYPE_NAME, MPI_CHAR, comm->intracomm);

  }

  /* Send mode */
  /*---------- */

  else { /* if (mode == SYR_COMM_MODE_SEND) */

    /* Assemble buffer */

    position = 0;

    MPI_Pack(sec_name, SYR_COMM_L_SEC_NAME, MPI_CHAR, buffer,
             SYR_COMM_MPI_PACK_SIZE, &position, comm->intracomm);

    MPI_Pack(n_sec_elts, 1, MPI_INT, buffer, SYR_COMM_MPI_PACK_SIZE,
             &position, comm->intracomm);

    if (*n_sec_elts > 0)
      MPI_Pack(type_name, SYR_COMM_L_TYPE_NAME, MPI_CHAR, buffer,
               SYR_COMM_MPI_PACK_SIZE, &position, comm->intracomm);

    /* Send message */

    ierror = MPI_Send(buffer,
                      position,
                      MPI_PACKED,
                      comm->proc_rank + proc_id,
                      0,
                      comm->intracomm);

    if (ierror != MPI_SUCCESS)
      _comm_mpi_msg_err(comm,
                        proc_id,
                        ierror);
  }
}

/*--------------------------------------------------------------------------
 * Exchange section data using MPI
 *--------------------------------------------------------------------------*/

static void
_comm_mpi_body(void              *sec_elts,
               int                n_sec_elts,
               syr_type_t         type,
               const syr_comm_t  *comm,
               syr_comm_mode_t    mode,
               int                proc_id)
{
  int ierror = MPI_SUCCESS;
  int n_elts = n_sec_elts;

  MPI_Status    status;
  MPI_Datatype  datatype;

  assert(comm != NULL);
  assert(n_sec_elts >= 0);

  /* Set datatype */

  switch (type) {
  case SYR_TYPE_int:
    datatype = MPI_INT;
    break;
  case SYR_TYPE_float:
    datatype = MPI_FLOAT;
    break;
  case SYR_TYPE_double:
    datatype = MPI_DOUBLE;
    break;
  case SYR_TYPE_char:
    datatype = MPI_CHAR;
    break;
  default:
    assert(0);
  }

  /* Receive or send */
  /*---------------- */

  if (mode == SYR_COMM_MODE_RECEIVE)
    ierror = MPI_Recv(sec_elts,
                      n_elts,
                      datatype,
                      comm->proc_rank + proc_id,
                      MPI_ANY_TAG,
                      comm->intracomm,
                      &status);

  else /* if (mode == SYR_COMM_MODE_SEND) */
    ierror = MPI_Send(sec_elts,
                      n_elts,
                      datatype,
                      comm->proc_rank + proc_id,
                      0,
                      comm->intracomm);

  if (ierror != MPI_SUCCESS)
    _comm_mpi_msg_err(comm, proc_id, ierror);
}

#endif /* defined(HAVE_MPI) */

#if defined(HAVE_SOCKET)

/*--------------------------------------------------------------------------
 * Read a record from an interface socket
 *--------------------------------------------------------------------------*/

static void
_comm_read_sock(const syr_comm_t  *comm,
                int                proc_id,
                char              *rec,
                size_t             n,
                syr_type_t         type)
{
  size_t   start_id;
  size_t   end_id;
  size_t   n_loc;
  size_t   n_bytes;
  size_t   size;
  ssize_t  ret;

  assert(rec  != NULL);
  assert(comm != NULL);
  assert(comm->socket != NULL);

  /* Determine associated size */

  switch (type) {
  case SYR_TYPE_int:
    size = sizeof(int);
    break;
  case SYR_TYPE_float:
    size = sizeof(float);
    break;
  case SYR_TYPE_double:
    size = sizeof(double);
    break;
  case SYR_TYPE_char:
    size = sizeof(char);
    break;
  default:
    assert(0);
  }

  n_bytes = size * n;

  /* Read record from socket */

  start_id = 0;

  while (start_id < n_bytes) {

    end_id = SYR_MIN(start_id + SSIZE_MAX, n_bytes);
    n_loc = end_id - start_id;

    ret = read(comm->socket[proc_id], (void *)(rec + start_id), n_loc);

    if (ret < 1)
      ple_error(__FILE__, __LINE__, errno,
                "Communication %s (proc %d) :\n"
                "Erreur a la reception via un socket",
                comm->name, proc_id + 1);

    start_id += ret;

  }

  if (comm->swap_endian == 1 && size > 1)
    _swap_endian(rec, rec, size, n);
}

/*--------------------------------------------------------------------------
 * Write a record to an interface socket
 *--------------------------------------------------------------------------*/

static void
_comm_write_sock(const syr_comm_t  *comm,
                 int                proc_id,
                 const char        *rec,
                 size_t             n,
                 syr_type_t         type)
{
  size_t   start_id;
  size_t   end_id;
  size_t   n_loc;
  size_t   n_bytes;
  size_t   size;
  ssize_t  ret;

  char    *rec_tmp;

  assert(rec  != NULL);
  assert(comm != NULL);
  assert(comm->socket != NULL);

  /* Determine associated size */

  switch(type) {
  case SYR_TYPE_int:
    size = sizeof(int);
    break;
  case SYR_TYPE_double:
    size = sizeof(double);
    break;
  case SYR_TYPE_float:
    size = sizeof(float);
    break;
  case SYR_TYPE_char:
    size = sizeof(char);
    break;
  default:
    assert(0);
  }

  n_bytes = size * n;

  /* Convert if "little-endian" */

  if (comm->swap_endian == 1 && size != 1) {
    PLE_MALLOC(rec_tmp, n_bytes, char);
    _swap_endian(rec_tmp, rec, size, n);
  }
  else
    rec_tmp = NULL;

  /* Write record to socket */

  start_id = 0;

  while (start_id < n_bytes) {

    end_id = SYR_MIN(start_id + SSIZE_MAX, n_bytes);
    n_loc = end_id - start_id;

    if (rec_tmp == NULL)
      ret = write(comm->socket[proc_id], (const void *)(rec + start_id),
                  n_loc);
    else
      ret = write(comm->socket[proc_id], (const void *)(rec_tmp + start_id),
                  n_loc);

    if (ret < 1)
      ple_error(__FILE__, __LINE__, errno,
                "Communication %s (proc %d) :\n"
                "Erreur a l'envoi via un socket",
                comm->name, proc_id + 1);

    start_id += ret;
  }

  if (rec_tmp != NULL)
    PLE_FREE(rec_tmp);
}

#endif /* defined(HAVE_SOCKET) */

/*--------------------------------------------------------------------------
 * Write then read a "magic string" to check file format
 *--------------------------------------------------------------------------*/

static void
_comm_magic_string(syr_comm_t  *comm,
                   int          proc_id,
                   const char  *magic_string)
{
  char *comm_magic_string;

  int lng_magic_string = strlen(magic_string);

  PLE_MALLOC(comm_magic_string, lng_magic_string + 1, char);

  /* Write magic string */
  /*--------------------*/

  strncpy(comm_magic_string, magic_string, lng_magic_string);

  switch (comm->type) {

#if defined(HAVE_MPI)
  case SYR_COMM_TYPE_MPI:
    {
      int ierror = MPI_Send(comm_magic_string,
                            lng_magic_string, MPI_CHAR,
                            comm->proc_rank + proc_id,
                            0,
                            comm->intracomm);

      if (ierror != MPI_SUCCESS)
        _comm_mpi_msg_err(comm,
                          proc_id,
                          ierror);

    }
    break;
#endif /* HAVE_MPI */

#if defined(HAVE_SOCKET)
  case SYR_COMM_TYPE_SOCKET:
    {
      _comm_write_sock(comm,
                       proc_id,
                       (const char *)(comm_magic_string),
                       strlen(magic_string),
                       SYR_TYPE_char);
    }
    break;
#endif

  default:
    break;

  }

  /* Read magic string */
  /*--------------------*/

  switch (comm->type) {

#if defined(HAVE_MPI)
  case SYR_COMM_TYPE_MPI:
    {
      int ierror;
      MPI_Status status;

      ierror = MPI_Recv(comm_magic_string,
                        lng_magic_string,
                        MPI_CHAR,
                        comm->proc_rank + proc_id,
                        MPI_ANY_TAG,
                        comm->intracomm,
                        &status);

      if (ierror != MPI_SUCCESS)
        _comm_mpi_msg_err(comm,
                          proc_id,
                          ierror);

    }
    break;
#endif /* HAVE_MPI */

#if defined(HAVE_SOCKET)
  case SYR_COMM_TYPE_SOCKET:
    {
      _comm_read_sock(comm,
                      proc_id,
                      (char *)(comm_magic_string),
                      strlen(magic_string),
                      SYR_TYPE_char);
    }
    break;
#endif

  default:
    break;

  }
  comm_magic_string[lng_magic_string] = '\0';

  /* If the magic string does not correspond, we have an error */

  if (strcmp(comm_magic_string, magic_string) != 0)
    ple_error(__FILE__, __LINE__, 0,
              "Error a la lecture de : \"%s\".\n"
              "La version du format d'interface est incompatible.\n"
              "La chaine magique indique la version du format d'interface :\n"
              "chaine magique lue :      \"%s\"\n"
              "chaine magique attendue : \"%s\"",
              comm->name, comm_magic_string, magic_string);

  PLE_FREE(comm_magic_string);
}

#if defined(HAVE_MPI)

/*--------------------------------------------------------------------------
 * Create an MPI intracommunicator
 *--------------------------------------------------------------------------*/

static void
_comm_mpi_init(syr_comm_t  *comm)
{
  int local_range[2] = {-1, -1};
  int distant_range[2] = {-1, -1};

  printf(" Initialisation de la communication MPI: %s ... ", comm->name);
  fflush(stdout);

  ple_coupling_mpi_intracomm_create(MPI_COMM_WORLD,
                                    syr_glob_mpi_comm,
                                    comm->proc_rank,
                                    &(comm->intracomm),
                                    local_range,
                                    distant_range);

  printf("[ok]\n");
  printf("  Rangs locaux = [%d..%d], rangs distants = [%d..%d].\n\n",
         local_range[0], local_range[1] - 1,
         distant_range[0], distant_range[1] - 1);
  fflush(stdout);

  comm->proc_rank = distant_range[0];
}

/*--------------------------------------------------------------------------
 * Initialize an MPI communication by sending or reading a "magic string"
 * allowing format verification
 *--------------------------------------------------------------------------*/

static void
_comm_mpi_open(syr_comm_t  *comm,
               int          proc_id,
               const char  *magic_string)
{
  /* Write and read "magic string" */
  /*------------------------------ */

  _comm_magic_string(comm,
                     proc_id,
                     magic_string);
}

#endif /* (HAVE_MPI) */

#if defined(HAVE_SOCKET)

/*--------------------------------------------------------------------------
 * Connection for socket connection initialization
 *--------------------------------------------------------------------------*/

static void
_comm_sock_connect(syr_comm_t  *comm,
                   int          proc_id,
                   char        *host_name,
                   int          port_num)
{

#if defined(__linux__) || defined(__linux) || defined(linux)
  socklen_t  sock_len;
#else
  size_t     sock_len;
#endif

  struct sockaddr_in   sock_addr;
  struct hostent      *host_ent;

  /* Create socket interface descriptor */

  comm->socket[proc_id] = socket(AF_INET, SOCK_STREAM, 0);

  if (comm->socket[proc_id] == -1)
    ple_error(__FILE__, __LINE__, errno,
              "Erreur a l'initialisation de la communication par socket "
              "(proc %d).\n",
              proc_id + 1);

  /* Prepare connection */

  sock_len = sizeof(sock_addr);

  memset((char *) &sock_addr, 0, sock_len);

  sock_addr.sin_family = AF_INET;
  sock_addr.sin_addr.s_addr = inet_addr(host_name);

  if (sock_addr.sin_addr.s_addr == INADDR_NONE) {
    host_ent = gethostbyname(host_name);

    if (!host_ent)
      ple_error(__FILE__, __LINE__, 0,
                "Communication par socket : hote \"%s\" inconnu.",
                host_name);

    memcpy(&sock_addr.sin_addr, host_ent->h_addr_list[0], host_ent->h_length);
  }

  sock_addr.sin_port = port_num;

  if (comm->swap_endian == 1)
    _swap_endian((char *)&(sock_addr.sin_port),
                 (char *)&(sock_addr.sin_port),
                 sizeof(sock_addr.sin_port),
                 1);

  if ( connect(comm->socket[proc_id],
               (struct sockaddr *)&sock_addr, sock_len) < 0)
    ple_error(__FILE__, __LINE__, errno,
              "Communication par socket : erreur de connexion a\n"
              "%s (port %d)\n", host_name, port_num);

}

/*--------------------------------------------------------------------------
 *  Initialize socket communication
 *--------------------------------------------------------------------------*/

static void
_comm_sock_init(syr_comm_t  *comm,
                const char  *sock_serv)
{
  int  i, port_num, id, proc_id, name_len;

  char   str_len[7] = "      ";
  char  *host_name = NULL;

  /* Decode sock_serv string */

  for (id = strlen(sock_serv) - 1;
       id > 0 && sock_serv[id] != ':'; id--);

  port_num = atoi(sock_serv + id + 1);

  PLE_MALLOC(host_name, id + 1, char);
  strncpy(host_name, sock_serv, id);
  host_name[id] = '\0';

  /* Establish communication with rank 0 on Code_Saturne side */

  _comm_sock_connect(comm, 0, host_name, port_num);

  /* Receive number of ranks */

  if (read(comm->socket[0], str_len, 6) < 6)
    ple_error(__FILE__, __LINE__, errno,
              syr_glob_comm_err_socket, comm->name,  1);

  str_len[6] = '\0';
  comm->n_procs = atoi(str_len);

  if (comm->n_procs > 1) {
    PLE_REALLOC(comm->socket, comm->n_procs, int);
    for (i = 1; i < comm->n_procs; i++)
      comm->socket[i] = 0;
  }

  /* Receive max hostname size */

  if (comm->n_procs > 1) {

    /* Receive max hostname size */

    if (read(comm->socket[0], str_len, 4) < 4)
      ple_error(__FILE__, __LINE__, errno,
                syr_glob_comm_err_socket, comm->name,  1);

    str_len[4] = '\0';
    name_len = atoi(str_len);

    PLE_REALLOC(host_name, name_len + 1, char);

    for (proc_id = 1; proc_id < comm->n_procs; proc_id++) {

      /* Receive hostname */

      if (read(comm->socket[0], host_name, name_len) < name_len)
        ple_error(__FILE__, __LINE__, errno,
                  syr_glob_comm_err_socket, comm->name, proc_id + 1);

      host_name[name_len] = '\0';

      /* Receive port number */

      if (read(comm->socket[0], str_len, 6) < 6)
        ple_error(__FILE__, __LINE__, errno,
                  syr_glob_comm_err_socket, comm->name, proc_id + 1);

      str_len[6] = '\0';
      port_num = atoi(str_len);

      /* Establish connection */

      _comm_sock_connect(comm, proc_id, host_name, port_num);

    }

  } /* End of n_procs > 1 case */

  PLE_FREE(host_name);
}

/*--------------------------------------------------------------------------
 * Initialize a socket for communication
 *--------------------------------------------------------------------------*/

static void
_comm_sock_open(syr_comm_t  *comm,
                int          proc_id,
                const char  *magic_string)
{
  char  comm_header[] = CS_COMM_SOCKET_HEADER;

  /* Send initialization message to each proc */

  if (write(comm->socket[proc_id], comm_header, strlen(comm_header))
      < (ssize_t)strlen(comm_header))
    ple_error(__FILE__, __LINE__, errno,
              "Erreur dans la communication par socket.");

  /* Write or read "magic string" */

  _comm_magic_string(comm,
                     proc_id,
                     magic_string);
}

/*--------------------------------------------------------------------------
 * Close an interface socket
 *--------------------------------------------------------------------------*/

static void
_comm_sock_close(syr_comm_t  *comm,
                 int          proc_id)
{
  if (close(comm->socket[proc_id]) != 0)
    ple_error(__FILE__, __LINE__, errno,
              "Communication %s (proc %d) :\n"
              "Erreur a la fermeture d'un socket.",
              comm->name, proc_id + 1);

  comm->socket[proc_id] = -1;
}

#endif /* defined(HAVE_SOCKET) */

/*===========================================================================
 * Public function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Initialize a communicator
 *
 * arguments:
 *   coupling_num  <-- Coupling number
 *   cs_root_rank  <-- Root rank associated with Code_Saturne
 *   cs_n_ranks    <-- Number of MPI ranks associated with Code_Saturne
 *   sock_str      <-- Communicating host and port (hostname:port)
 *   type          <-- Communication type
 *   echo          <-- Echo on main output
 *----------------------------------------------------------------------------*/

syr_comm_t *
syr_comm_initialize(int               coupling_num,
                    int               cs_root_rank,
                    int               cs_n_ranks,
                    const char       *sock_str,
                    syr_comm_type_t   type,
                    int               echo)
{
  unsigned  int_endian;

  syr_comm_t *comm;
  int proc_id;

  const char magic_string[] = "CFD_SYRTHES_COUPLING_2.2";
  const char base_name[] = "CFD_";

  PLE_MALLOC(comm, 1, syr_comm_t);

  /* Build communicator name */
  /*------------------------ */

  PLE_MALLOC(comm->name, strlen(base_name) + 4 + 1, char);

  sprintf(comm->name, "%s%04d", base_name, coupling_num);

  /* Initialize other fields */

  comm->type = type;
  comm->echo = echo;

  comm->n_sec_elts = NULL;

  /* Test if system is big-endian */

  comm->swap_endian = 0; /* Use "big-endian" mode to communicate */

  int_endian = 0;
  *((char *) (&int_endian)) = '\1';

  if (int_endian == 1)
    comm->swap_endian = 1;

#if defined(DEBUG) && !defined(NDEBUG)
  else {
    int_endian = 0;
    *((char *) (&int_endian) + sizeof(unsigned) - 1) = '\1';
    assert (int_endian == 1);
  }
#endif

#if defined(HAVE_MPI)
  comm->n_procs = -1;
  comm->proc_rank = -1;
#endif

#if defined(HAVE_SOCKET)
  comm->socket = NULL;
#endif

  if (type == SYR_COMM_TYPE_MPI) {
#if defined(HAVE_MPI)
    comm->proc_rank = cs_root_rank;
    comm->n_procs = cs_n_ranks;
#else
    ple_error
      (__FILE__, __LINE__, 0,
       "Librarie compilee sans support MPI, donc le type communicateur\n"
       "doit etre different de SYR_COMM_TYPE_MPI (%d).",
       (int)SYR_COMM_TYPE_MPI);
#endif
  }

  if (type == SYR_COMM_TYPE_SOCKET) {
#if defined(HAVE_SOCKET)
    PLE_MALLOC(comm->socket, 1, int);
    comm->socket[0] = 0;

    /* Using sockets communication, the number of associated ranks
       may be updated during communication initialization */

    _comm_sock_init(comm, sock_str);
#else
    ple_error
      (__FILE__, __LINE__, 0,
       "Librarie compilee sans support \"socket\", donc le type communicateur\n"
       "doit etre different de SYR_COMM_TYPE_SOCKET (%d).",
       (int)SYR_COMM_TYPE_SOCKET);
#endif
  }

  /* Jump one line in log file */

  printf("\n");

#if defined(HAVE_MPI)
  if (comm->type == SYR_COMM_TYPE_MPI)
    _comm_mpi_init(comm);
#endif

  /* Build interfaces */

  for (proc_id = 0; proc_id < comm->n_procs; proc_id++) {

    /* Build interface file name */
    /*-------------------------- */

    /* Info on interface creation */

    if (comm->n_procs == 1)
      printf("  Ouverture de la communication :  %s ...",
             comm->name);
    else
      printf("  Ouverture de la communication :  %s (proc %4d) ...",
             comm->name, proc_id + 1);

    fflush(stdout);

    /* Create a descriptor for the interface file */
    /*------------------------------------------- */

#if defined(HAVE_MPI)
    if (comm->type == SYR_COMM_TYPE_MPI)
      _comm_mpi_open(comm, proc_id, magic_string);
#endif

#if defined(HAVE_SOCKET)
  if (comm->type == SYR_COMM_TYPE_SOCKET)
    _comm_sock_open(comm, proc_id, magic_string);
#endif

    /* Info on interface creation success */

    printf(" [ok]\n");
    fflush(stdout);

  }

  PLE_MALLOC(comm->n_sec_elts, comm->n_procs, int);
  for (proc_id = 0; proc_id < comm->n_procs; proc_id++)
    comm->n_sec_elts[proc_id] = 0;

  /* Fin */

  return comm;
}


/*--------------------------------------------------------------------------
 * Finalize a communicator
 *--------------------------------------------------------------------------*/

syr_comm_t *
syr_comm_finalize(syr_comm_t *comm)
{
  int proc_id;

  printf("\n");

  if (comm->n_procs == 1) {

    /* Info on closing interface files */

    printf("  Fermeture de la communication: %s\n",
           comm->name);

#if defined(HAVE_SOCKET)
    if (comm->socket != NULL)
      _comm_sock_close(comm, 0);
#endif

  }
  else {

    for (proc_id = 0; proc_id < comm->n_procs; proc_id++) {

      /* Info on closing interface files */

      printf("  Fermeture de la communication: %s (proc %4d)\n",
             comm->name, proc_id + 1);

#if defined(HAVE_SOCKET)
      if (comm->socket != NULL)
        _comm_sock_close(comm, proc_id);
#endif

    }

  }

#if defined(HAVE_SOCKET)
  if (comm->socket != NULL)
    PLE_FREE(comm->socket);
#endif

  PLE_FREE(comm->name);
  PLE_FREE(comm->n_sec_elts);

  PLE_FREE(comm);

  return NULL;
}

/*--------------------------------------------------------------------------
 * Get a communicator's name
 *
 * This function returns a pointer to an existing name, so the string
 * returned should not be freed by the user.
 *--------------------------------------------------------------------------*/

const char *
syr_comm_get_name(const syr_comm_t *comm)
{
  return comm->name;
}

/*--------------------------------------------------------------------------
 * Get a communicator's number of distant processes
 *--------------------------------------------------------------------------*/

int
syr_comm_get_n_procs(const syr_comm_t *comm)
{
  return comm->n_procs;
}

/*--------------------------------------------------------------------------
 * Get the number of values to be read for each sub-domain in a communicator.
 *
 * This function should only be called between syr_comm_read_header()
 * and syr_comm_read_data(). It returns a pointer to an array belonging to
 * the communicator, which should not be freed.
 *--------------------------------------------------------------------------*/

const int *
syr_comm_get_n_section_elements(const syr_comm_t *comm)
{
 return comm->n_sec_elts;
}

/*--------------------------------------------------------------------------
 * Write a section to the communication interface.
 *
 * A section contains:
 *
 * 1. section name
 * 2. number of elements
 * 3. element type
 * 4. element values
 *--------------------------------------------------------------------------*/

void
syr_comm_write_section(const char        *sec_name,
                       int                n_elts,
                       void              *elts,
                       syr_type_t         elt_type,
                       const syr_comm_t  *comm,
                       int                proc_id)
{
  char sec_name_out [SYR_COMM_L_SEC_NAME + 1];
  char type_name    [SYR_COMM_L_TYPE_NAME + 1];
  char type_name_out[SYR_COMM_L_TYPE_NAME + 1];

  assert(comm != NULL);
  assert(n_elts >= 0);

  /* Section name */

  sprintf(sec_name_out,
          "%-*.*s",
          SYR_COMM_L_SEC_NAME,
          SYR_COMM_L_SEC_NAME,
          sec_name);

  /* Element type name */

  if (n_elts != 0) {

    switch(elt_type) {

    case SYR_TYPE_int:
      strcpy(type_name, "i ");
      break;

    case SYR_TYPE_double:
      strcpy(type_name, "r8");
      break;

    case SYR_TYPE_float:
      strcpy(type_name, "r4");
      break;

    case SYR_TYPE_char:
      strcpy(type_name, "c ");
      break;

    default:
      assert(0);
    }

    sprintf(type_name_out,
            "%-*.*s",
            SYR_COMM_L_TYPE_NAME,
            SYR_COMM_L_TYPE_NAME,
            type_name);

  }

  if (comm->echo >= 0)
    _comm_echo_pre(comm, proc_id, SYR_COMM_MODE_SEND);

#if defined(HAVE_MPI)

  /* MPI communication */
  /*------------------ */

  if (comm->type == SYR_COMM_TYPE_MPI) {

    int  n_elts_wri = n_elts;

    _comm_mpi_header(sec_name_out,
                     &n_elts_wri,
                     type_name_out,
                     comm,
                     SYR_COMM_MODE_SEND,
                     proc_id);

    if (n_elts > 0)
      _comm_mpi_body((void *)elts,
                     n_elts,
                     elt_type,
                     comm,
                     SYR_COMM_MODE_SEND,
                     proc_id);

  }

#endif /* (HAVE_MPI) */

#if defined(HAVE_SOCKET)

  /* Socket communication */
  /*--------------------- */

  if (comm->type == SYR_COMM_TYPE_SOCKET) {

    /* 1. section name */

    _comm_write_sock(comm,
                     proc_id,
                     (const char *)sec_name_out,
                     SYR_COMM_L_SEC_NAME,
                     SYR_TYPE_char);

    /* 2. number of elements */

    _comm_write_sock(comm,
                     proc_id,
                     (const char *)(&n_elts),
                     1,
                     SYR_TYPE_int);

    if (n_elts != 0) {

      /* 3. element type name */

      _comm_write_sock(comm,
                       proc_id,
                       (const char *)type_name_out,
                       SYR_COMM_L_TYPE_NAME,
                       SYR_TYPE_char);

      /* 4. element values */

      _comm_write_sock(comm,
                       proc_id,
                       (const char *)elts,
                       (size_t)n_elts,
                       elt_type);

    }

  }

#endif /* (HAVE_SOCKET) */

  /* Optional logging */
  /*----------------- */

  if (comm->echo  >= 0)
    _comm_echo_header(comm,
                      sec_name_out,
                      n_elts,
                      type_name_out);

  if (comm->echo > 0)
    _comm_echo_body(comm,
                    n_elts,
                    elt_type,
                    elts);
}


/*--------------------------------------------------------------------------
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
 *--------------------------------------------------------------------------*/

void
syr_comm_read_header(char              *sec_name,
                     int               *n_elts,
                     syr_type_t        *elt_type,
                     const syr_comm_t  *comm,
                     int                proc_id)
{
  char   type_name[SYR_COMM_L_TYPE_NAME + 1];

  assert(comm  != NULL);

  *n_elts = 0;

  if (comm->echo >= 0)
    _comm_echo_pre(comm,
                   proc_id,
                   SYR_COMM_MODE_RECEIVE);

#if defined(HAVE_MPI)

  /* MPI communication */
  /*------------------ */

  if (comm->type == SYR_COMM_TYPE_MPI) {

    _comm_mpi_header(sec_name,
                     &(comm->n_sec_elts[proc_id]),
                     type_name,
                     comm,
                     SYR_COMM_MODE_RECEIVE,
                     proc_id);

    *n_elts = comm->n_sec_elts[proc_id];

  }

#endif /* (HAVE_MPI) */

#if defined(HAVE_SOCKET)

  /* Socket communication */
  /*--------------------- */

  if (comm->type == SYR_COMM_TYPE_SOCKET) {

    /* 1. section name */

    _comm_read_sock(comm,
                    proc_id,
                    (char *)sec_name,
                    SYR_COMM_L_SEC_NAME,
                    SYR_TYPE_char);

    sec_name[SYR_COMM_L_SEC_NAME] = '\0';

    /* 2. number of elements */

    _comm_read_sock(comm,
                    proc_id,
                    (char *) &(comm->n_sec_elts[proc_id]),
                    1,
                    SYR_TYPE_int);

    *n_elts = comm->n_sec_elts[proc_id];

    if (comm->n_sec_elts[proc_id] != 0) {

      /* 3. element type name */

      _comm_read_sock(comm,
                      proc_id,
                      (char *)type_name,
                      SYR_COMM_L_TYPE_NAME,
                      SYR_TYPE_char);

    }

  }

#endif /* (HAVE_SOCKET) */

  sec_name[SYR_COMM_L_SEC_NAME] = '\0';
  type_name[SYR_COMM_L_TYPE_NAME] = '\0';

  /* Optional logging */
  /*----------------- */

  if (comm->echo >= 0)
    _comm_echo_header(comm,
                      sec_name,
                      comm->n_sec_elts[proc_id],
                      type_name);

  if (comm->n_sec_elts[proc_id] != 0) {

    type_name[SYR_COMM_L_TYPE_NAME] = '\0';

    if (strcmp(type_name, "i ") == 0)
      *elt_type = SYR_TYPE_int;

    else if (strcmp(type_name, "r4") == 0)
      *elt_type = SYR_TYPE_float;

    else if (strcmp(type_name, "r8") == 0)
      *elt_type = SYR_TYPE_double;

    else if (strcmp(type_name, "c ") == 0)
      *elt_type = SYR_TYPE_char;

    else
      assert(0);
  }

}

/*--------------------------------------------------------------------------
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
 *--------------------------------------------------------------------------*/

void *
syr_comm_read_body(int                n_sec_elts,
                   void              *sec_elts,
                   syr_type_t         elt_type,
                   const syr_comm_t  *comm,
                   int                proc_id)
{
  void   *_sec_elts;

  assert(comm  != NULL);
  assert(n_sec_elts >= 0);

  _sec_elts = sec_elts;

  if (_sec_elts == NULL && n_sec_elts != 0) {

    switch(elt_type) {

    case SYR_TYPE_int:
      {
        int *sec_elts_int;

        PLE_MALLOC(sec_elts_int, n_sec_elts, int);
        _sec_elts = (void *)sec_elts_int;
      }
      break;

    case SYR_TYPE_float:
      {
        float *sec_elts_flo;

        PLE_MALLOC(sec_elts_flo, n_sec_elts, float);
        _sec_elts = (void *)sec_elts_flo;
      }
      break;

    case SYR_TYPE_double:
      {
        double *sec_elts_dou;

        PLE_MALLOC(sec_elts_dou, n_sec_elts, double);
        _sec_elts = (void *)sec_elts_dou;
      }
      break;

    case SYR_TYPE_char:
      {
        char *sec_elts_cha;

        PLE_MALLOC(sec_elts_cha, n_sec_elts + 1, char);
        _sec_elts = (void *)sec_elts_cha;
      }
      break;

    default:
      assert(elt_type == SYR_TYPE_int    ||
             elt_type == SYR_TYPE_float  ||
             elt_type == SYR_TYPE_double ||
             elt_type == SYR_TYPE_char   );

    }

  }

  /* 5. Element values */

  if (n_sec_elts != 0) {

    /* Effective read */

#if defined(HAVE_MPI)
    if (comm->type == SYR_COMM_TYPE_MPI)
      _comm_mpi_body((void *)_sec_elts,
                     n_sec_elts,
                     elt_type,
                     comm,
                     SYR_COMM_MODE_RECEIVE,
                     proc_id);
#endif /* (HAVE_MPI) */

#if defined(HAVE_SOCKET)
    if (comm->type == SYR_COMM_TYPE_SOCKET)
      _comm_read_sock(comm,
                      proc_id,
                      (char *)_sec_elts,
                      (size_t)n_sec_elts,
                      elt_type);
#endif /* (HAVE_SOCKET) */

    /* Verification */
    /*------------- */

    if (elt_type == SYR_TYPE_char) {

      char  *sec_elts_cha = (char *)_sec_elts;

      sec_elts_cha[comm->n_sec_elts[proc_id]] = '\0';

    }

    /* Optional logging */
    /*----------------- */

    if (comm->echo > 0)
      _comm_echo_body(comm,
                      n_sec_elts,
                      elt_type,
                      _sec_elts);

  }

  /* Return pointer to values read */

  return _sec_elts;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
