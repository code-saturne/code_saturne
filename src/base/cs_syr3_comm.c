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

/*============================================================================
 *  Communication with SYRTHES 3
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_CS_HAVE_MPI)
#include <mpi.h>
#endif

#if defined(_CS_HAVE_SOCKET)
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_file.h>
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

#define CS_SYR3_COMM_LNG_NOM_TYPE_ELT        2    /* Length of type name */

#define CS_SYR3_COMM_SOCKET_ENTETE            "CS_comm_socket"

#define CS_SYR3_COMM_SOCKET_NBR_MAX          8
#define CS_LOC_SYR3_COMM_LNG_HOSTNAME      256
#define CS_LOC_SYR3_COMM_LNG_NOM_MAX       256

/*
  If SSIZE_MAX has not been defined through system headers, we take the
  minimal value required by POSIX (for low level read/write used with sockets).
*/

#if !defined(SSIZE_MAX)
#define SSIZE_MAX  32767
#endif

/*============================================================================
 * Local Structure Definitions
 *============================================================================*/

struct _cs_syr3_comm_t {

  char                 *nom;          /* Communicator name */

  cs_int_t              rang_proc;    /* Communication process name */
  int                   sock;         /* Socket number */

  cs_syr3_comm_mode_t   mode;         /* Communication mode */
  cs_syr3_comm_type_t   type;         /* Type of data encoding */
  cs_bool_t             swap_endian;  /* Swap bytes ? */
  cs_int_t              echo;         /* Data transfer verbosity level */

};

/*============================================================================
 *  Global variables
 *============================================================================*/

static char  cs_syr3_comm_nom_typ_elt_char[] = "c ";  /* String */
static char  cs_syr3_comm_nom_typ_elt_int[]  = "i ";  /* Integer */
static char  cs_syr3_comm_nom_typ_elt_real[] = "r8";  /* Real */

#if defined(_CS_HAVE_SOCKET)

static cs_bool_t       cs_glob_comm_little_endian = false;

static char  cs_glob_comm_sock_nom_hote[CS_LOC_SYR3_COMM_LNG_HOSTNAME + 1];
static int   cs_glob_comm_sock_num_port = -1;

static int             cs_glob_comm_socket = 0;
struct sockaddr_in     cs_glob_comm_addr_sock;

static char  cs_glob_comm_err_socket[]
= N_("Error in socket communication:  %s (node %4d)\n");

#endif /* _CS_HAVE_SOCKET */

/*============================================================================
 * Private function definitions
 *============================================================================*/

#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Print an error message in case of MPI communication problem
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_mpi_msg_err(const cs_syr3_comm_t  *comm,
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
 * Initialize an MPI communication by sending or receiving an eventual
 * "magic string" used to check the correct data format.
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_mpi_ouvre(cs_syr3_comm_t  *const  comm,
                           const char      *const  chaine_magique)
{
  int ierror, comm_size;

  MPI_Status status;

  char * chaine_magique_comm;

  cs_int_t lng_chaine_magique = strlen(chaine_magique);

  /*--------------------------*/
  /* Initialize communication */
  /*--------------------------*/

  assert(   comm->mode == CS_SYR3_COMM_MODE_RECEPTION
         || comm->mode == CS_SYR3_COMM_MODE_EMISSION);

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  if (comm->rang_proc >= comm_size)

    bft_error(__FILE__, __LINE__, 0,
              _("Impossible to establish the communication: %s\n"
                "because the requested process rank (%d)\n"
                "is greater than or equal to the number of MPI processes (%d)"),
              comm->nom, comm->rang_proc, comm_size);

  BFT_MALLOC(chaine_magique_comm, lng_chaine_magique + 1, char);

  /*------------------------------*/
  /* Receive or send magic string */
  /*------------------------------*/

  if (comm->mode == CS_SYR3_COMM_MODE_RECEPTION) {

    ierror = MPI_Recv(chaine_magique_comm, lng_chaine_magique, MPI_CHAR,
                      comm->rang_proc,
                      MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (ierror != MPI_SUCCESS)
      cs_loc_syr3_comm_mpi_msg_err(comm, ierror);

    chaine_magique_comm[lng_chaine_magique] = '\0';

    /* If the magic string does not match, we have an error */

    if (strcmp(chaine_magique_comm, chaine_magique) != 0)
      bft_error
        (__FILE__, __LINE__, 0,
         _("Error for communication: \"%s\".\n"
           "The interface version is not correct.\n"
           "The magic string indicates an incorrect interface format version.\n"
           "magic string read:     \"%s\"\n"
           "magic string expected: \"%s\"\n"),
         comm->nom, chaine_magique_comm, chaine_magique);

  }
  else if (comm->mode == CS_SYR3_COMM_MODE_EMISSION) {

    strncpy(chaine_magique_comm, chaine_magique, lng_chaine_magique);

    ierror = MPI_Send(chaine_magique_comm, lng_chaine_magique, MPI_CHAR,
                      comm->rang_proc,
                      0, MPI_COMM_WORLD);

    if (ierror != MPI_SUCCESS)
      cs_loc_syr3_comm_mpi_msg_err(comm, ierror);

  }

  BFT_FREE(chaine_magique_comm);

}

/*----------------------------------------------------------------------------
 * Exchange a section header through MPI
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_mpi_entete(char                  *nom_rub,
                            cs_int_t              *nbr_elt_rub,
                            char                  *nom_typ_elt,
                            const cs_syr3_comm_t  *comm)
{
#undef  CS_SYR3_COMM_MPI_PACK_SIZE
#define CS_SYR3_COMM_MPI_PACK_SIZE    CS_SYR3_COMM_H_LEN \
                                    + CS_SYR3_COMM_LNG_NOM_TYPE_ELT \
                                    + (sizeof(int) * 2)

  char buffer[CS_SYR3_COMM_MPI_PACK_SIZE];

  int position, ierror;

  MPI_Status  status;

  /* Instructions */

  assert(comm != NULL);
  assert(*nbr_elt_rub >= 0);
  assert(sizeof(int) == sizeof(cs_int_t));

  /* Receive mode */
  /*--------------*/

  if (comm->mode == CS_SYR3_COMM_MODE_RECEPTION) {

    /* Receive message */

    ierror = MPI_Recv(buffer, CS_SYR3_COMM_MPI_PACK_SIZE, MPI_PACKED,
                      comm->rang_proc,
                      MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if (ierror != MPI_SUCCESS)
      cs_loc_syr3_comm_mpi_msg_err(comm, ierror);

    /* Extract buffer elements */

    position = 0;

    MPI_Unpack(buffer, CS_SYR3_COMM_MPI_PACK_SIZE, &position, nom_rub,
               CS_SYR3_COMM_H_LEN, MPI_CHAR, MPI_COMM_WORLD);

    MPI_Unpack(buffer, CS_SYR3_COMM_MPI_PACK_SIZE, &position, nbr_elt_rub,
               1, CS_MPI_INT, MPI_COMM_WORLD);

    if (*nbr_elt_rub > 0)
      MPI_Unpack(buffer, CS_SYR3_COMM_MPI_PACK_SIZE, &position, nom_typ_elt,
                 CS_SYR3_COMM_LNG_NOM_TYPE_ELT, MPI_CHAR, MPI_COMM_WORLD);

  }

  /* Send mode */
  /*-----------*/

  else if (comm->mode == CS_SYR3_COMM_MODE_EMISSION) {

    /* Pack buffer */

    position = 0;

    MPI_Pack(nom_rub, CS_SYR3_COMM_H_LEN, MPI_CHAR, buffer,
             CS_SYR3_COMM_MPI_PACK_SIZE, &position, MPI_COMM_WORLD);

    MPI_Pack(nbr_elt_rub, 1, CS_MPI_INT, buffer, CS_SYR3_COMM_MPI_PACK_SIZE,
             &position, MPI_COMM_WORLD);

    if (*nbr_elt_rub > 0)
      MPI_Pack(nom_typ_elt, CS_SYR3_COMM_LNG_NOM_TYPE_ELT, MPI_CHAR, buffer,
               CS_SYR3_COMM_MPI_PACK_SIZE, &position, MPI_COMM_WORLD);

    /* Send message */

    ierror = MPI_Send(buffer, position, MPI_PACKED, comm->rang_proc,
                      0, MPI_COMM_WORLD);

    if (ierror != MPI_SUCCESS)
      cs_loc_syr3_comm_mpi_msg_err(comm, ierror);
  }

  else
    assert(   comm->mode == CS_SYR3_COMM_MODE_RECEPTION
           || comm->mode == CS_SYR3_COMM_MODE_EMISSION);
}

/*----------------------------------------------------------------------------
 * Exchange a section body through MPI
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_mpi_corps(void                  *elt_rub,
                           cs_int_t               nbr_elt_rub,
                           cs_type_t              typ_elt,
                           const cs_syr3_comm_t  *comm)
{
  int ierror;
  int nbr_elt = nbr_elt_rub;

  MPI_Status  status;

  /* Instructions */

  assert(comm != NULL);
  assert(nbr_elt_rub >= 0);

  /* Receive mode */
  /*--------------*/

  if (comm->mode == CS_SYR3_COMM_MODE_RECEPTION) {

    switch (typ_elt) {

    case CS_TYPE_cs_int_t:

      ierror = MPI_Recv(elt_rub, nbr_elt, CS_MPI_INT,
                        comm->rang_proc,
                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      break;

    case CS_TYPE_cs_real_t:
      ierror = MPI_Recv(elt_rub, nbr_elt, CS_MPI_REAL,
                        comm->rang_proc,
                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      break;

    case CS_TYPE_char:
      ierror = MPI_Recv(elt_rub, nbr_elt, MPI_CHAR,
                        comm->rang_proc,
                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      break;

    default:
      assert (   typ_elt == CS_TYPE_char
              || typ_elt == CS_TYPE_cs_int_t
              || typ_elt == CS_TYPE_cs_real_t);
    }

  }

  /* Send mode */
  /*-----------*/

  else if (comm->mode == CS_SYR3_COMM_MODE_EMISSION) {

    switch (typ_elt) {

    case CS_TYPE_cs_int_t:
      ierror = MPI_Send(elt_rub, nbr_elt, CS_MPI_INT,
                        comm->rang_proc,
                        0, MPI_COMM_WORLD);
      break;

    case CS_TYPE_cs_real_t:
      ierror = MPI_Send(elt_rub, nbr_elt, CS_MPI_REAL,
                        comm->rang_proc,
                        0, MPI_COMM_WORLD);
      break;

    case CS_TYPE_char:
      ierror = MPI_Send(elt_rub, nbr_elt, MPI_CHAR,
                        comm->rang_proc,
                        0, MPI_COMM_WORLD);
      break;

    default:
      assert(   typ_elt == CS_TYPE_char
             || typ_elt == CS_TYPE_cs_int_t
             || typ_elt == CS_TYPE_cs_real_t);
    }

  }

  else
    assert(   comm->mode == CS_SYR3_COMM_MODE_RECEPTION
           || comm->mode == CS_SYR3_COMM_MODE_EMISSION);

  if (ierror != MPI_SUCCESS)
    cs_loc_syr3_comm_mpi_msg_err(comm, ierror);
}

#endif /* (_CS_HAVE_MPI) */

#if defined(_CS_HAVE_SOCKET)

/*----------------------------------------------------------------------------
 * Read a record from the interface socket
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_lit_sock(const cs_syr3_comm_t  *comm,
                          cs_byte_t             *rec,
                          const size_t           nbr,
                          cs_type_t              type)
{
  size_t   ind_deb;
  size_t   ind_fin;
  size_t   nbr_loc;
  size_t   nbr_octet;
  size_t   taille;
  ssize_t  ret;

  assert(rec  != NULL);
  assert(comm != NULL);

  /* Determine the number of bytes to receive */

  switch(type) {
  case CS_TYPE_cs_int_t:
    taille = sizeof(cs_int_t);
    break;
  case CS_TYPE_cs_real_t:
    taille = sizeof(cs_real_t);
    break;
  case CS_TYPE_char:
    taille = sizeof(char);
    break;
  default:
    assert(type == CS_TYPE_cs_int_t  ||
           type == CS_TYPE_cs_real_t ||
           type == CS_TYPE_char);
  }

  nbr_octet = taille * nbr;

  /* Read data from socket */
  /*-----------------------*/

  ind_deb = 0;

  while (ind_deb < nbr_octet) {

    ind_fin = CS_MIN(ind_deb + SSIZE_MAX, nbr_octet);

    nbr_loc = ind_fin - ind_deb;

    ret = read(comm->sock, (void *)(rec + ind_deb), nbr_loc);

    if (ret < 1)
      bft_error(__FILE__, __LINE__, errno,
                _("Communication %s:\n"
                  "Error while receiving data by socket.\n"),
                comm->nom);

    ind_deb += ret;

  }

  if (comm->swap_endian == true)
    bft_file_swap_endian(rec, rec, taille, nbr);

}

/*----------------------------------------------------------------------------
 * Write a record to the interface socket
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_ecrit_sock(const cs_syr3_comm_t  *comm,
                            const cs_byte_t       *rec,
                            const size_t           nbr,
                            cs_type_t              type)
{
  size_t   ind_deb;
  size_t   ind_fin;
  size_t   nbr_loc;
  size_t   nbr_octet;
  size_t   taille;
  ssize_t  ret;

  cs_byte_t   * rec_tmp;

  assert(rec  != NULL);
  assert(comm != NULL);

  /* Determine the number of bytes to send */

  switch(type) {
  case CS_TYPE_cs_int_t:
    taille = sizeof(cs_int_t);
    break;
  case CS_TYPE_cs_real_t:
    taille = sizeof(cs_real_t);
    break;
  case CS_TYPE_char:
    taille = sizeof(char);
    break;
  default:
    assert(type == CS_TYPE_cs_int_t  ||
           type == CS_TYPE_cs_real_t ||
           type == CS_TYPE_char);
  }

  nbr_octet = taille * nbr;

  /* Convert to "big-endian" */

  if (comm->swap_endian == true && taille != 1) {
    BFT_MALLOC(rec_tmp, nbr_octet, cs_byte_t);
    bft_file_swap_endian(rec_tmp, rec, taille, nbr);
  }
  else
    rec_tmp = NULL;

  /* write data to socket */
  /*----------------------*/

  ind_deb = 0;

  while (ind_deb < nbr_octet) {

    ind_fin = CS_MIN(ind_deb + SSIZE_MAX, nbr_octet);

    nbr_loc = ind_fin - ind_deb;

    if (rec_tmp == NULL)
      ret = write(comm->sock, (const void *)(rec + ind_deb), nbr_loc);
    else
      ret = write(comm->sock, (const void *)(rec_tmp + ind_deb), nbr_loc);

    if (ret < 1)
      bft_error(__FILE__, __LINE__, errno,
                _("Communication %s:\n"
                  "Error sending data by socket.\n"),
                comm->nom);

    ind_deb += ret;

  }

  if (rec_tmp != NULL)
    BFT_FREE(rec_tmp);
}

/*----------------------------------------------------------------------------
 * Initialize socket communication
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_sock_connect(cs_syr3_comm_t  *comm)
{
  int ind;

#if defined(_CS_ARCH_Linux)
  socklen_t long_sock;
#else
  int       long_sock;  /* size_t according to SUS-v2 standard, but acording
                           to "man gethostbyname" under Linux, the standard
                           is bad, and we should have an int (or socklen_t) */
#endif

  char   str_taille[6] = "     ";
  char  *host_names = NULL;
  int   *tab_num_port = NULL;

#if defined (_CS_HAVE_MPI)
  int ierror = MPI_SUCCESS;
#endif
  int rang = (cs_glob_base_rang == -1 ? 0 : cs_glob_base_rang);

  const int lng_hostname = CS_LOC_SYR3_COMM_LNG_HOSTNAME + 1;

  /* Connect to server socket */

  long_sock = sizeof(cs_glob_comm_addr_sock);

  if (rang == 0) {

    comm->sock = accept(cs_glob_comm_socket,
                        (struct sockaddr *)&cs_glob_comm_addr_sock,
                        &long_sock);

    /* Send number of ranks */

    sprintf(str_taille, "%5d", (int)cs_glob_base_nbr);

    if (write(comm->sock, str_taille, 6) < 6)
      bft_error(__FILE__, __LINE__, errno,
                _("Error in socket communication\n"));
  }

  /* Obtains the name of the host machine and its port number on rank 0 */

  if (cs_glob_base_nbr > 1) {

    BFT_MALLOC(host_names,
               lng_hostname * cs_glob_base_nbr,
               char);

    BFT_MALLOC(tab_num_port, cs_glob_base_nbr, int);

#if defined(_CS_HAVE_MPI)
    ierror = MPI_Gather(cs_glob_comm_sock_nom_hote, lng_hostname, MPI_CHAR,
                        host_names, lng_hostname, MPI_CHAR, 0,
                        cs_glob_base_mpi_comm);

    if (ierror < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error while sending the host name through MPI in sockets "
                  "initialization.\n"));

    /* Send the port number */

    ierror = MPI_Gather(&cs_glob_comm_sock_num_port, 1, MPI_INT,
                        tab_num_port, 1, MPI_INT, 0, cs_glob_base_mpi_comm);

    if (ierror < 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error while sending the port number through MPI in sockets "
                  "initialization.\n"));

    if (rang != 0)
      comm->sock = accept(cs_glob_comm_socket,
                          (struct sockaddr *)&cs_glob_comm_addr_sock,
                          &long_sock);

#else
    bft_error(__FILE__, __LINE__, 0,
              _("MPI is needed for socket initialization.\n"));
#endif

    /* rank 0 sends the number of ranks, hostnames, and port numbers */

    if (rang == 0) {

      /* Send max hostname size */

      sprintf(str_taille, "%3d", lng_hostname);

      if (write(comm->sock, str_taille, 4) < 4)
        bft_error(__FILE__, __LINE__, errno,
                  _("Error in socket communication\n"));

      for (ind = 1; ind < cs_glob_base_nbr; ind++) {

        /* Send host machine name */

        if (write(comm->sock, &(host_names[lng_hostname*ind]), lng_hostname)
            < lng_hostname)
          bft_error(__FILE__, __LINE__, errno,
                    _("Error in socket communication\n"));

        /* Send port number */

        sprintf(str_taille, "%5d", tab_num_port[ind]);

        if (write(comm->sock, str_taille, 6) < 6)
          bft_error(__FILE__, __LINE__, errno,
                    _("Error in socket communication\n"));

      }

    } /* End of rank-0 specific operations */

    BFT_FREE(host_names);
    BFT_FREE(tab_num_port);

  } /* End for cs_glob_base_nbr > 1 */

}

/*----------------------------------------------------------------------------
 * Ensure exchange of magic string through sockets
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_sock_ouvre(cs_syr3_comm_t   *comm,
                            const char       *chaine_magique)
{
  char nom_tmp[32 + 1];
  int taille;

  int rang = (cs_glob_base_rang == -1 ? 0 : cs_glob_base_rang);

  taille = strlen(CS_SYR3_COMM_SOCKET_ENTETE);

  if (read(comm->sock, nom_tmp, taille) < taille)
    bft_error(__FILE__, __LINE__, errno,
              _(cs_glob_comm_err_socket), comm->nom,
              rang + 1);

  /* Check that the connection is from the correct application type */

  if (strncmp(nom_tmp, CS_SYR3_COMM_SOCKET_ENTETE, taille != 0))
    bft_error(__FILE__, __LINE__, 0,
              _("Attempt to connect to the communication port with\n"
                "an unknown message format\n"));

  /*----------------------------*/
  /* Write or read magic string */
  /*----------------------------*/

  if (comm->mode == CS_SYR3_COMM_MODE_RECEPTION) {

    char      *chaine_magique_lue;
    cs_int_t   lng_chaine_magique = strlen(chaine_magique);

    BFT_MALLOC(chaine_magique_lue, lng_chaine_magique + 1, char);

    cs_loc_syr3_comm_lit_sock(comm,
                         (void *)(chaine_magique_lue),
                         strlen(chaine_magique),
                         CS_TYPE_char);

    chaine_magique_lue[lng_chaine_magique] = '\0';

    /* If the magic string does not match, we have an error */

    if (strcmp(chaine_magique_lue, chaine_magique) != 0)
      bft_error(__FILE__, __LINE__, 0,
                _("Error while initializating communication: \"%s\".\n"
                  "The interface version is not correct.\n"
                  "The magic string indicates the interface format version:\n"
                  "magic string read:     \"%s\"\n"
                  "magic string expected: \"%s\"\n"),
                comm->nom, chaine_magique_lue, chaine_magique);

    BFT_FREE(chaine_magique_lue);

  }
  else if (comm->mode == CS_SYR3_COMM_MODE_EMISSION)

    cs_loc_syr3_comm_ecrit_sock(comm,
                           (const void *)(chaine_magique),
                           strlen(chaine_magique),
                           CS_TYPE_char);
}

/*----------------------------------------------------------------------------
 * Close the connection with the interface socket
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_sock_ferme(cs_syr3_comm_t  *comm)
{
  if (close(comm->sock) != 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Communication %s):\n"
                "Error closing the socket.\n"),
              comm->nom);

  comm->sock = -1;
}

#endif /* (_CS_HAVE_SOCKET) */

/*----------------------------------------------------------------------------
 * Print information on waiting for a message
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_echo_pre(const cs_syr3_comm_t  *comm)
{
  assert(comm != NULL);

  switch(comm->mode) {

  case CS_SYR3_COMM_MODE_RECEPTION:
    bft_printf(_("\nMessage received on \"%s\":\n"), comm->nom);
    break;

  case CS_SYR3_COMM_MODE_EMISSION:
    bft_printf(_("\nMessage sent on \"%s\":\n"), comm->nom);
    break;

  default:
    assert(   comm->mode == CS_SYR3_COMM_MODE_RECEPTION
           || comm->mode == CS_SYR3_COMM_MODE_EMISSION);
  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Print a message header
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_echo_entete(const char  *nom_rub,
                             cs_int_t     nbr_elt,
                             cs_type_t    typ_elt
)
{
  char nom_rub_ecr[CS_SYR3_COMM_H_LEN + 1];

  /* instructions */

  strncpy(nom_rub_ecr, nom_rub,  CS_SYR3_COMM_H_LEN);
  nom_rub_ecr[CS_SYR3_COMM_H_LEN] = '\0';

  bft_printf(_("    section name:          \"%s\"\n"
               "    number of elements:    %d\n"),
             nom_rub_ecr, nbr_elt);

  if (nbr_elt > 0) {

    char *nom_typ;

    switch(typ_elt) {
    case CS_TYPE_char:
      nom_typ = cs_syr3_comm_nom_typ_elt_char;
      break;
    case CS_TYPE_cs_int_t:
      nom_typ = cs_syr3_comm_nom_typ_elt_int;
      break;
    case CS_TYPE_cs_real_t:
      nom_typ = cs_syr3_comm_nom_typ_elt_real;
      break;
    default:
      assert(   typ_elt == CS_TYPE_char
             || typ_elt == CS_TYPE_cs_int_t
             || typ_elt == CS_TYPE_cs_real_t);
    }

    bft_printf(_("    element type name:      \"%s\"\n"), nom_typ);

  }

  bft_printf_flush();
}

/*----------------------------------------------------------------------------
 * Print (partial) message content
 *----------------------------------------------------------------------------*/

static void
cs_loc_syr3_comm_echo_donnees(cs_int_t     echo,
                              cs_int_t     nbr_elt,
                              cs_type_t    typ_elt,
                              const void  *elt_rub
)
{
  cs_int_t  echo_deb = 0;
  cs_int_t  echo_fin;
  cs_int_t  ind;

  /* Instructions */

  if (nbr_elt == 0) return;

  if (echo * 2 < nbr_elt) {
    echo_fin = echo;
    bft_printf(_("    %d first and last elements:\n"), echo);
  }
  else {
    echo_fin = nbr_elt;
    bft_printf(_("    elements:\n"));
  }

  do {

    switch (typ_elt) {

    case CS_TYPE_cs_int_t:
      {
        const cs_int_t *elt_rub_int = (const cs_int_t *) elt_rub;

        for (ind = echo_deb ; ind < echo_fin ; ind++)
          bft_printf("    %10d : %12d\n", ind + 1, *(elt_rub_int + ind));
      }
      break;

    case CS_TYPE_cs_real_t:
      {
        const cs_real_t *elt_rub_real = (const cs_real_t *) elt_rub;

        for (ind = echo_deb ; ind < echo_fin ; ind++)
          bft_printf("    %10d : %12.5e\n", ind + 1, *(elt_rub_real + ind));
      }
      break;

    case CS_TYPE_char:
      {
        const char *elt_rub_char = (const char *) elt_rub;

        for (ind = echo_deb ; ind < echo_fin ; ind++) {
          if (*(elt_rub_char + ind) != '\0')
            bft_printf("    %10d : '%c'\n", ind + 1, *(elt_rub_char + ind));
          else
            bft_printf("    %10d : '\\0'\n", ind + 1);
        }
      }
      break;

    default:

      assert(   typ_elt == CS_TYPE_cs_int_t
             || typ_elt == CS_TYPE_cs_real_t
             || typ_elt == CS_TYPE_char);

    }

    if (echo_fin < nbr_elt) {
      bft_printf("    ..........   ............\n");
      echo_deb = nbr_elt - echo;
      echo_fin = nbr_elt;
    }
    else {
      assert(echo_fin == nbr_elt);
      echo_fin = nbr_elt + 1;
    }

  } while (echo_fin <= nbr_elt);

  bft_printf_flush();

}

/*============================================================================
 * Public function definitions
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
                        const cs_int_t             echo)
{
  unsigned    int_endian;

  const char *base_name[] = {"receive_from_SYRTHES_",
                             "send_to_SYRTHES_"};
  const char *chaine_magique[] = {"SYRTHES_TO_SATURNE_2.2",
                                  "SATURNE_TO_SYRTHES_2.2"};
  cs_syr3_comm_t  *comm = NULL;

  BFT_MALLOC(comm, 1, cs_syr3_comm_t);

  /* Build communicator name */

  BFT_MALLOC(comm->nom, strlen(base_name[mode]) + 4 + 1, char);

  sprintf(comm->nom, "%s%04d", base_name[mode], numero);

  /* Initialize other fields */

  comm->mode = mode;
  comm->type = type;
  comm->echo = echo;

#if defined(_CS_HAVE_MPI)
  comm->rang_proc = rang_proc;
#else
  comm->rang_proc = -1;
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

  /* Information on interface creation */

  bft_printf(_("\n  Opening communication:  %s..."), comm->nom);
  bft_printf_flush();

#if defined(_CS_HAVE_SOCKET)
  if (comm->type == CS_SYR3_COMM_TYPE_SOCKET)
    cs_loc_syr3_comm_sock_connect(comm);
#endif /* (_CS_HAVE_SOCKET) */

  /* Create interface file descriptor */
  /*----------------------------------*/

  if (comm->type == CS_SYR3_COMM_TYPE_MPI) {

#if defined(_CS_HAVE_MPI)
    cs_loc_syr3_comm_mpi_ouvre(comm, chaine_magique[mode]);
#else
    assert(comm->rang_proc < 0);
#endif

  }
  else {

#if defined(_CS_HAVE_SOCKET)
    if (comm->type == CS_SYR3_COMM_TYPE_SOCKET)
      cs_loc_syr3_comm_sock_ouvre(comm, chaine_magique[mode]);
#endif

  }

  /* Info on interface creation success */

  bft_printf(" [ok]\n");
  bft_printf_flush();

  return comm;
}

/*----------------------------------------------------------------------------
 * Function finalizing a communication
 *----------------------------------------------------------------------------*/

cs_syr3_comm_t *
cs_syr3_comm_termine(cs_syr3_comm_t *comm)
{
  /* Info on interface finalization */

  bft_printf(_("\n  Closing communication:  %s\n"), comm->nom);
  bft_printf_flush();

#if defined(_CS_HAVE_SOCKET)

  if (comm->type == CS_SYR3_COMM_TYPE_SOCKET)
    cs_loc_syr3_comm_sock_ferme(comm);

#endif /* (_CS_HAVE_SOCKET) */

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
cs_syr3_comm_ret_nom(const cs_syr3_comm_t  *comm)
{
  assert(comm != NULL);

  return(comm->nom);
}

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
                            const cs_syr3_comm_t  *comm)
{
  char   nom_rub_ecr[CS_SYR3_COMM_H_LEN + 1];

  char  *nom_typ_elt;
  char   nom_typ_elt_ecr[CS_SYR3_COMM_LNG_NOM_TYPE_ELT + 1];


  assert(comm != NULL);
  assert(nbr_elt >= 0);

  /* section name */

  sprintf(nom_rub_ecr,
          "%-*.*s",
          CS_SYR3_COMM_H_LEN,
          CS_SYR3_COMM_H_LEN,
          nom_rub);

  /* element type name */

  if (nbr_elt != 0) {

    switch(typ_elt) {

    case CS_TYPE_cs_int_t:
      nom_typ_elt = cs_syr3_comm_nom_typ_elt_int;
      break;

    case CS_TYPE_cs_real_t:
      nom_typ_elt = cs_syr3_comm_nom_typ_elt_real;
      break;

    case CS_TYPE_char:
      nom_typ_elt = cs_syr3_comm_nom_typ_elt_char;
      break;

    default:
      assert(   typ_elt == CS_TYPE_cs_int_t
             || typ_elt == CS_TYPE_cs_real_t
             || typ_elt == CS_TYPE_char);

    }

    sprintf(nom_typ_elt_ecr,
            "%-*.*s",
            CS_SYR3_COMM_LNG_NOM_TYPE_ELT,
            CS_SYR3_COMM_LNG_NOM_TYPE_ELT,
            nom_typ_elt);

  }

  if (comm->echo  >= 0)
    cs_loc_syr3_comm_echo_pre(comm);


#if defined(_CS_HAVE_MPI)

  /* MPI communication */
  /*-------------------*/

  if (comm->type == CS_SYR3_COMM_TYPE_MPI) {

    cs_int_t  nbr_elt_rub_ecr = nbr_elt;

    cs_loc_syr3_comm_mpi_entete(nom_rub_ecr,
                                &nbr_elt_rub_ecr,
                                nom_typ_elt_ecr,
                                comm);

    if (nbr_elt > 0)
      cs_loc_syr3_comm_mpi_corps((void *) elt,
                                 nbr_elt,
                                 typ_elt,
                                 comm);

  }

#endif /* (_CS_HAVE_MPI) */

#if defined(_CS_HAVE_SOCKET)

  /* socket communication */
  /*----------------------*/

  if (comm->type == CS_SYR3_COMM_TYPE_SOCKET) {

    /* section name */

    cs_loc_syr3_comm_ecrit_sock(comm,
                                (const void *) nom_rub_ecr,
                                CS_SYR3_COMM_H_LEN,
                                CS_TYPE_char);

    /* number of elements */

    cs_loc_syr3_comm_ecrit_sock(comm,
                                (const void *)(&nbr_elt),
                                1,
                                CS_TYPE_cs_int_t);

    if (nbr_elt != 0) {

      /* element type name */

      cs_loc_syr3_comm_ecrit_sock(comm,
                                  (const void *) nom_typ_elt_ecr,
                                  CS_SYR3_COMM_LNG_NOM_TYPE_ELT,
                                  CS_TYPE_char);

      /* element values */

      cs_loc_syr3_comm_ecrit_sock(comm,
                                  (const void *) elt,
                                  (size_t) nbr_elt,
                                  typ_elt);

    }

  }

#endif /* (_CS_HAVE_SOCKET) */

  /* Possibly print to log file */

  if (comm->echo  >= 0)
    cs_loc_syr3_comm_echo_entete(nom_rub,
                                 nbr_elt,
                                 typ_elt);

  if (comm->echo > 0)
    cs_loc_syr3_comm_echo_donnees(comm->echo,
                                  nbr_elt,
                                  typ_elt,
                                  elt);
}

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
                           const cs_syr3_comm_t       *comm)
{
  char   nom_typ_elt[CS_SYR3_COMM_LNG_NOM_TYPE_ELT + 1];

  assert(comm  != NULL);

  entete->nbr_elt = 0;

  if (comm->echo >= 0)
    cs_loc_syr3_comm_echo_pre(comm);


#if defined(_CS_HAVE_MPI)

  /* MPI communication */
  /*-------------------*/

  if (comm->type == CS_SYR3_COMM_TYPE_MPI) {

    cs_loc_syr3_comm_mpi_entete(entete->nom_rub,
                                &(entete->nbr_elt),
                                nom_typ_elt,
                                comm);

  }

#endif /* (_CS_HAVE_MPI) */

#if defined(_CS_HAVE_SOCKET)

  /* socket communication */
  /*----------------------*/

  if (comm->type == CS_SYR3_COMM_TYPE_SOCKET) {

    /* section type name */

    cs_loc_syr3_comm_lit_sock(comm,
                              (void *) &(entete->nom_rub),
                              CS_SYR3_COMM_H_LEN,
                              CS_TYPE_char);

    /* number of elements */

    cs_loc_syr3_comm_lit_sock(comm,
                              (void *) &(entete->nbr_elt),
                              1,
                              CS_TYPE_cs_int_t);


    if (entete->nbr_elt != 0) {

      /* element type name */

      cs_loc_syr3_comm_lit_sock(comm,
                                (void *) nom_typ_elt,
                                CS_SYR3_COMM_LNG_NOM_TYPE_ELT,
                                CS_TYPE_char);

    }

  }

#endif /* (_CS_HAVE_SOCKET) */

  entete->nom_rub[CS_SYR3_COMM_H_LEN] = '\0';

  if (entete->nbr_elt != 0) {

    nom_typ_elt[CS_SYR3_COMM_LNG_NOM_TYPE_ELT] = '\0';

    if (strcmp(nom_typ_elt, cs_syr3_comm_nom_typ_elt_int) == 0)
      entete->typ_elt = CS_TYPE_cs_int_t;

    else if (strcmp(nom_typ_elt, cs_syr3_comm_nom_typ_elt_real) == 0)
      entete->typ_elt = CS_TYPE_cs_real_t;

    else if (strcmp(nom_typ_elt, cs_syr3_comm_nom_typ_elt_char) == 0)
      entete->typ_elt = CS_TYPE_char;

    else {
      assert(0);
    }
  }

  /* Possibly print to log file */

  if (comm->echo >= 0)
    cs_loc_syr3_comm_echo_entete(entete->nom_rub,
                                 entete->nbr_elt,
                                 entete->typ_elt);

  /* Return number of elements to read */

  return entete->nbr_elt;
}

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
                          const cs_syr3_comm_t             *comm)
{
  cs_int_t    ind;
  void      *_elt_rub;

  assert(comm  != NULL);
  assert(entete->nbr_elt >= 0);

  _elt_rub = elt;

  if (_elt_rub == NULL && entete->nbr_elt != 0) {

    switch(entete->typ_elt) {

    case CS_TYPE_cs_int_t:
      {
        cs_int_t  *elt_rub_int;

        BFT_MALLOC(elt_rub_int, entete->nbr_elt, cs_int_t);
        _elt_rub = (void *) elt_rub_int;
      }
      break;

    case CS_TYPE_cs_real_t:
      {
        cs_real_t  *elt_rub_rea;

        BFT_MALLOC(elt_rub_rea, entete->nbr_elt, cs_real_t);
        _elt_rub = (void *)elt_rub_rea;
      }
      break;

    case CS_TYPE_char:
      {
        char  *elt_rub_cha;

        BFT_MALLOC(elt_rub_cha, entete->nbr_elt + 1, char);
        _elt_rub = (void *)elt_rub_cha;
      }
      break;

    default:
      assert(   entete->typ_elt == CS_TYPE_cs_int_t
             || entete->typ_elt == CS_TYPE_cs_real_t
             || entete->typ_elt == CS_TYPE_char);
    }

  }

  /* element values */

  if (entete->nbr_elt != 0) {

#if defined(_CS_HAVE_MPI)

    if (comm->type == CS_SYR3_COMM_TYPE_MPI)
      cs_loc_syr3_comm_mpi_corps((void *)_elt_rub,
                                 entete->nbr_elt,
                                 entete->typ_elt,
                                 comm);

#endif /* (_CS_HAVE_MPI) */

#if defined(_CS_HAVE_SOCKET)

    if (comm->type == CS_SYR3_COMM_TYPE_SOCKET)
      cs_loc_syr3_comm_lit_sock(comm,
                                (void *)_elt_rub,
                                (size_t) entete->nbr_elt,
                                entete->typ_elt);

#endif /* (_CS_HAVE_SOCKET) */

    /* Verification */

    if (entete->typ_elt == CS_TYPE_char) {
      for (ind = 0 ;
           ind < entete->nbr_elt && ((char *)_elt_rub)[ind] != '\0' ;
           ind++);
      ((char *)_elt_rub)[ind] = '\0';
    }

    /* Possibly print to log file */

    if (comm->echo > 0)
      cs_loc_syr3_comm_echo_donnees(comm->echo,
                                    entete->nbr_elt,
                                    entete->typ_elt,
                                    _elt_rub);

  }

}

#if defined(_CS_HAVE_SOCKET)

/*----------------------------------------------------------------------------
 * Open an IP socket to prepare for this communication mode
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_init_socket(void)
{
  char       chaine[CS_LOC_SYR3_COMM_LNG_HOSTNAME + 1];

  int        nbr_connect_max;
  int        num_port;

#if defined(_CS_ARCH_Linux)
  socklen_t long_sock;
#else
  int       long_sock;  /* size_t according to SUS-v2 standard, but acording
                           to "man gethostbyname" under Linux, the standard
                           is bad, and we should have an int (or socklen_t) */
#endif

  unsigned  int_endian;

  struct sockaddr_in   addr_sock;
  struct hostent      *ent_hote;


  int rang = (cs_glob_base_rang == -1 ? 0 : cs_glob_base_rang);

  /* Initialization */

  nbr_connect_max = 0;

  if (getenv("CS_SYR3_COMM_SOCKET_NBR_MAX") != NULL)
    nbr_connect_max = atoi(getenv("CS_SYR3_COMM_SOCKET_NBR_MAX"));

  if (nbr_connect_max == 0)
    nbr_connect_max = CS_SYR3_COMM_SOCKET_NBR_MAX;

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

  long_sock = sizeof(addr_sock);

  memset((char *) &addr_sock, 0, long_sock);

  addr_sock.sin_family = AF_INET;
  addr_sock.sin_addr.s_addr = INADDR_ANY;
  addr_sock.sin_port = 0;

  if (cs_glob_comm_little_endian == true) {
    bft_file_swap_endian(&(addr_sock.sin_addr.s_addr),
                         &(addr_sock.sin_addr.s_addr),
                         sizeof(addr_sock.sin_addr.s_addr),
                         1);
    bft_file_swap_endian(&(addr_sock.sin_port),
                         &(addr_sock.sin_port),
                         sizeof(addr_sock.sin_port),
                         1);
  }

  if (gethostname(chaine, CS_LOC_SYR3_COMM_LNG_HOSTNAME) < 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Error obtaining computer's name"));
  chaine[CS_LOC_SYR3_COMM_LNG_HOSTNAME] = '\0';

  ent_hote = gethostbyname(chaine);
  memcpy(ent_hote->h_addr_list[0], &addr_sock.sin_addr, ent_hote->h_length);

  if (bind(cs_glob_comm_socket,
           (struct sockaddr *)&addr_sock,
           long_sock) != 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Initialization error for socket communication support.\n"));

  if (listen(cs_glob_comm_socket, nbr_connect_max) < 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Initialization error for socket communication support.\n"));

  /* Obtain assigned service number */

  if (getsockname(cs_glob_comm_socket,
                  (struct sockaddr *)&addr_sock,
                  &long_sock) != 0)
    bft_error(__FILE__, __LINE__, errno,
              _("Initialization error for socket communication support.\n"));

  num_port = addr_sock.sin_port;
  if (cs_glob_comm_little_endian == true) {
    bft_file_swap_endian(&(addr_sock.sin_port),
                         &(addr_sock.sin_port),
                         sizeof(addr_sock.sin_port), 1);
    num_port = addr_sock.sin_port;
    bft_file_swap_endian(&(addr_sock.sin_port),
                         &(addr_sock.sin_port),
                         sizeof(addr_sock.sin_port), 1);
  }

  /* Save the structure in the associated global variable */

  cs_glob_comm_addr_sock = addr_sock;

  /* Write host and port names in process order */

  if (rang == 0) {

    /* Print available socket information to log for rank  0
       (do not internationalize this string so that scripts
       my use it more easily). */

    bft_printf("\n  SYRTHES: port <%s:%d>\n\n", chaine, num_port);
    bft_printf_flush();

  }

  memcpy(cs_glob_comm_sock_nom_hote, chaine, CS_LOC_SYR3_COMM_LNG_HOSTNAME);
  cs_glob_comm_sock_nom_hote[CS_LOC_SYR3_COMM_LNG_HOSTNAME] = '\0';
  cs_glob_comm_sock_num_port = num_port;
}

/*----------------------------------------------------------------------------
 * Close an IP socket associated with this communication mode
 *----------------------------------------------------------------------------*/

void
cs_syr3_comm_termine_socket(void)
{
  if (cs_glob_comm_socket == 0)
    return;

  close(cs_glob_comm_socket);

  bft_printf(_("\nClosing socket...\t [ok]\n"));
  bft_printf_flush();
}

#endif /* _CS_HAVE_SOCKET */

/*----------------------------------------------------------------------------*/

END_C_DECLS
