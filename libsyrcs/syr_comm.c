/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2011 EDF S.A., France
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

/*===========================================================================
 * Definitions of base communication functions
 *===========================================================================*/

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

#include <mpi.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>
#include <ple_coupling.h>

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

#define SYR_COMM_L_TYPE_NAME         2

typedef enum {

  SYR_COMM_MODE_RECEIVE,   /* Receive */
  SYR_COMM_MODE_SEND       /* Send */

} syr_comm_mode_t;

/*===========================================================================
 * Structure definition
 *===========================================================================*/

struct _syr_comm_t {

  char            *name;         /* Communicator name */

  int              n_procs;      /* Number of communicating processes */
  int              echo;         /* Log (printout) level of communications */

  int             *n_sec_elts;   /* Number of elements in a section for each
                                    proc (when reading) */

  int              proc_rank;    /* Rank of first distant process */
  MPI_Comm         intracomm;    /* Intracommunicator */

};

/*===========================================================================
 * Global variables
 *===========================================================================*/

/*===========================================================================
 * Private function definitions
 *===========================================================================*/

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

/*--------------------------------------------------------------------------
 * Write then read a "magic string" to check file format
 *--------------------------------------------------------------------------*/

static void
_comm_magic_string(syr_comm_t  *comm,
                   int          proc_id,
                   const char  *magic_string)
{
  int ierror;
  MPI_Status status;
  char *comm_magic_string;

  int lng_magic_string = strlen(magic_string);

  PLE_MALLOC(comm_magic_string, lng_magic_string + 1, char);

  /* Write magic string */
  /*--------------------*/

  strncpy(comm_magic_string, magic_string, lng_magic_string);

  ierror = MPI_Send(comm_magic_string,
                    lng_magic_string, MPI_CHAR,
                    comm->proc_rank + proc_id,
                    0,
                    comm->intracomm);

  if (ierror != MPI_SUCCESS)
    _comm_mpi_msg_err(comm,
                      proc_id,
                      ierror);

  /* Read magic string */
  /*--------------------*/

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
 *   type          <-- Communication type
 *   echo          <-- Echo on main output
 *----------------------------------------------------------------------------*/

syr_comm_t *
syr_comm_initialize(int               coupling_num,
                    int               cs_root_rank,
                    int               cs_n_ranks,
                    int               echo)
{
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

  comm->echo = echo;

  comm->n_sec_elts = NULL;

  comm->n_procs = -1;
  comm->proc_rank = -1;

  comm->proc_rank = cs_root_rank;
  comm->n_procs = cs_n_ranks;

  /* Jump one line in log file */

  printf("\n");

  _comm_mpi_init(comm);

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

    _comm_mpi_open(comm, proc_id, magic_string);

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

  }
  else {

    for (proc_id = 0; proc_id < comm->n_procs; proc_id++) {

      /* Info on closing interface files */

      printf("  Fermeture de la communication: %s (proc %4d)\n",
             comm->name, proc_id + 1);

    }

  }

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
  int  n_elts_wri;
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

  /* MPI communication */
  /*------------------ */

  n_elts_wri = n_elts;

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

  /* MPI communication */
  /*------------------ */

  _comm_mpi_header(sec_name,
                   &(comm->n_sec_elts[proc_id]),
                   type_name,
                   comm,
                   SYR_COMM_MODE_RECEIVE,
                   proc_id);

  *n_elts = comm->n_sec_elts[proc_id];

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

    _comm_mpi_body((void *)_sec_elts,
                   n_sec_elts,
                   elt_type,
                   comm,
                   SYR_COMM_MODE_RECEIVE,
                   proc_id);

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
