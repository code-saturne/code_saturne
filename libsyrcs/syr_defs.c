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

/*============================================================================
 * Definitions, Global variables, and basic functions
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_GETCWD)
#include <unistd.h>
#include <errno.h>
#endif

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>

#if defined(HAVE_MPI)
#include <ple_coupling.h>
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "syr_defs.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Global variables
 *============================================================================*/

char syr_glob_build_date[] = __DATE__;  /* Build date */

#if defined(HAVE_MPI)
MPI_Comm  syr_glob_mpi_comm = MPI_COMM_NULL;
ple_coupling_mpi_set_t *syr_glob_coupling_world = NULL;
#endif

/*===========================================================================
 * Private function definitions
 *===========================================================================*/

/*----------------------------------------------------------------------------
 * Basic error handler
 *----------------------------------------------------------------------------*/

static void _syr_error_handler
(
 const char     *const filename,
 const int             linenum,
 const int             sys_err_code,
 const char     *const format,
       va_list         arg_ptr
)
{
  fflush(stdout);

  fprintf(stderr, "\n"
                  "Erreur a l'execution de Syrthes\n"
                  "===============================\n");

  if (sys_err_code != 0)
    fprintf(stderr, "\nErreur Systeme : %s\n", strerror(sys_err_code));

  fprintf(stderr, "\n%s:%d: Erreur fatale.\n\n", filename, linenum);

  vfprintf(stderr, format, arg_ptr);

  fprintf(stderr, "\n\n");

  assert(0);   /* Use of assert() here allows interception by a debugger. */

#if defined(HAVE_MPI)
  {
    int mpi_flag;

    MPI_Initialized(&mpi_flag);

    if (mpi_flag != 0) {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }
#endif /* HAVE_MPI */

  exit(EXIT_FAILURE);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * First analysis of the command line to determine if we require MPI,
 * and initialization if necessary
 *
 * parameters:
 *   argc  <-> number of command line arguments
 *   argv  <-> array of command line arguments
 *----------------------------------------------------------------------------*/

static void
_mpi_test_and_initialize(int    *argc,
                         char  **argv[])
{
  int arg_id = 0, flag = 0;
  int use_mpi = 0;

#if   defined(__blrts__) || defined(__bgp__) \
   || defined(__CRAYXT_COMPUTE_LINUX_TARGET)
  use_mpi = 1;

#elif defined(MPICH2)
  if (getenv("PMI_RANK") != NULL)
    use_mpi = 1;

#elif defined(MPICH_NAME)

  /*
    Using standard MPICH1 1.2.x with the p4 (default) mechanism,
    the information required by MPI_Init() are transferred through
    the command line, which is then modified by MPI_Init();
    in this case, only rank 0 knows the "user" command line arguments
    at program startup, the other processes obtaining them only upon
    calling  MPI_Init(). In this case, it is thus necessary to initialize
    MPI before parsing the command line.
  */

  for (arg_id = 0; arg_id < *argc; arg_id++) {
    if (   !strcmp((*argv)[arg_id], "-p4pg")         /* For process 0 */
        || !strcmp((*argv)[arg_id], "-p4rmrank")) {  /* For other processes */
      use_mpi = 1;
      break;
    }
  }

  if (getenv("GMPI_ID") != NULL) /* In case we are using MPICH-GM */
    use_mpi = 1;

#elif defined(LAM_MPI)
  if (getenv("LAMRANK") != NULL)
    use_mpi = 1;

#elif defined(OPEN_MPI)
  if (getenv("OMPI_MCA_ns_nds_vpid") != NULL)
    use_mpi = 1;
  else if (getenv("OMPI_COMM_WORLD_RANK") != NULL)
    use_mpi = 1;

#endif /* Tests for known MPI variants */

  /* If we have determined from known MPI environment variables
     of command line arguments that we are running under MPI,
     initialize MPI */

  if (use_mpi == 1) {
    MPI_Initialized(&flag);
    if (!flag)
      MPI_Init(argc, argv);
  }

  /* Loop on command line arguments */

  arg_id = 0;

  while (++arg_id < *argc) {

    if (strcmp((*argv)[arg_id], "-comm-mpi") == 0)
      use_mpi = 1;

  } /* End of loop on command line arguments */

  if (use_mpi == 1) {
    MPI_Initialized(&flag);
    if (!flag)
      MPI_Init(argc, argv);
  }
}

/*----------------------------------------------------------------------------
 * First analysis of the command line to determine an application name.
 *
 * If no name is defined by the command line, a name is determined based
 * on the working directory.
 *
 * The caller is responsible for freeing the returned string.
 *
 * parameters:
 *   argc  <-- number of command line arguments
 *   argv  <-- array of command line arguments
 *
 * returns:
 *   pointer to character string with application name
 *----------------------------------------------------------------------------*/

static char *
_app_name(int    argc,
          char  *argv[])
{
  char *app_name = NULL;
  int arg_id = 0;

  /* Loop on command line arguments */

  arg_id = 0;

  while (++arg_id < argc) {
    char *s = argv[arg_id];
    if (strcmp(s, "--app-name") == 0) {
      if (arg_id + 1 < argc) {
        PLE_MALLOC(app_name, strlen(argv[arg_id + 1]) + 1, char);
        strcpy(app_name, argv[arg_id + 1]);
      }
    }
  }

  /* Use execution directory if name is unavailable */

  if (app_name == NULL) {

#if defined(HAVE_GETCWD)

    int i;
    int buf_size = 128;
    char *wd = NULL, *buf = NULL;

    while (wd == NULL) {
      buf_size *= 2;
      PLE_REALLOC(buf, buf_size, char);
      wd = getcwd(buf, buf_size);
      if (wd == NULL && errno != ERANGE)
        ple_error(__FILE__, __LINE__, errno,
                  "Erreur d'interrogation du repertoire d'execution.\n");
    }
    for (i = strlen(buf) - 1; i > 0 && buf[i-1] != '/'; i--);
    PLE_MALLOC(app_name, strlen(buf + i) + 1, char);
    strcpy(app_name, buf + i);
    PLE_FREE(buf);

#else
    ple_error(__FILE__, __LINE__, 0,
              "Nom d'instance SYRTHES (option --app-name) non fourni.\n");
#endif
  }

  return app_name;
}

#endif /* HAVE_MPI */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Error initialization and handling
 *----------------------------------------------------------------------------*/

void
syr_errhandler_initialize(void)
{
  ple_error_handler_set(_syr_error_handler);
}

#if defined (HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize MPI communication if necessary
 *----------------------------------------------------------------------------*/

void
syr_mpi_initialize(int    *argc,
                   char  **argv[])
{
  int mpi_init_flag, rank, syr_comm_rank, syr_comm_size;
  int app_num = -1, inactive = 0;
  int sync_flag = PLE_COUPLING_NO_SYNC;
  char *app_name = NULL;
  const char *app_type[2] = {"SYRTHES 3.4", "Unused (SYRTHES 3.4)"};

  /* Initialize MPI */

  _mpi_test_and_initialize(argc, argv);

  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    return;

  app_name = _app_name(*argc, *argv);

  app_num = ple_coupling_mpi_name_to_id(MPI_COMM_WORLD, app_name);

  /*
    Split MPI_COMM_WORLD to separate different coupled applications
    (collective operation, like all MPI communicator creation operations).

    app_num is equal to -1 if all applications have the same instance
    name, in which case no communicator split is necessary.
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (app_num > -1) {
    int ierror = MPI_Comm_split(MPI_COMM_WORLD, app_num, rank, &syr_glob_mpi_comm);
    if (ierror != 0)
      ple_error(__FILE__, __LINE__, 0,
                "Erreur a la creation d'un communicateur local a SYRTHES.");
  }
  else
    ple_error(__FILE__, __LINE__, 0,
              "Une seule application SYRTHES presente dans MPI_COMM_WORLD.");

  MPI_Comm_rank(syr_glob_mpi_comm, &syr_comm_rank);
  MPI_Comm_size(syr_glob_mpi_comm, &syr_comm_size);

  /* If extra ranks are assigned to SYRTHES,
     separate rank 0 (useful) from others (ignored) */

  if (syr_comm_size > 1) {
    MPI_Comm  mpi_comm_syr_ini = syr_glob_mpi_comm;
    syr_glob_mpi_comm = MPI_COMM_NULL;
    if (syr_comm_rank > 0)
      inactive = 1;
    if (MPI_Comm_split(mpi_comm_syr_ini, inactive, syr_comm_rank,
                       &syr_glob_mpi_comm) != 0)
      ple_error(__FILE__, __LINE__, 0,
                "Erreur a la subdivision d'un communicateur local a SYRTHES.\n");
    MPI_Comm_free(&mpi_comm_syr_ini);
  }

  /* Discover other applications in the same MPI root communicator
     (and participate in corresponding communication). */

  if (inactive > 0) {
    char *new_name = NULL;
    sync_flag = PLE_COUPLING_NO_SYNC;
    PLE_MALLOC(new_name, strlen(app_name) + strlen("Unused ()") + 1, char);
    sprintf(new_name, "Unused (%s)", app_name);
    PLE_FREE(app_name);
    app_name = new_name;
  }

  syr_glob_coupling_world = ple_coupling_mpi_set_create(sync_flag,
                                                        app_type[inactive],
                                                        app_name,
                                                        MPI_COMM_WORLD,
                                                        syr_glob_mpi_comm);

  PLE_FREE(app_name);


  /* If more than 1 rank was assigned to SYRTHES, only 1 is active
     (this may be the case in IBM Blue Gene/P, where a whole pset
     must be assigned to a given executable) */

  if (syr_comm_rank > 0) {
    syr_mpi_finalize();
    syr_exit(EXIT_SUCCESS);
  }
}

/*----------------------------------------------------------------------------
 * Finalize MPI communication
 *----------------------------------------------------------------------------*/

void
syr_mpi_finalize(void)
{
  int ierror = 0;

  ple_coupling_mpi_set_destroy(&syr_glob_coupling_world);

  assert(syr_glob_coupling_world == NULL);

  MPI_Comm_free(&syr_glob_mpi_comm);

  ierror = MPI_Barrier(MPI_COMM_WORLD);

  if (ierror != 0)
    ple_error(__FILE__, __LINE__, 0,
              "Erreur dans MPI_Barrier lors de la finalisation du\n"
              "communicateur global cote SYRTHES.");

  ierror = MPI_Finalize();

  if (ierror != 0)
    ple_error(__FILE__, __LINE__, 0,
              "Erreur lors de la finalisation du\n"
              "communicateur global cote SYRTHES.");
}

/*----------------------------------------------------------------------------
 * Force abort of MPI communication (for atexit)
 *
 * This function only forces finalization if an MPI coupling world is defined,
 * so it should do nothing in case of a normal exit, in which the MPI coupling
 * world info structure should have been destroyed through a regular call
 * to syr_mpi_finalize.
 *----------------------------------------------------------------------------*/

void
syr_mpi_exit_force(void)
{
  if (syr_glob_coupling_world != NULL)
    MPI_Abort(MPI_COMM_WORLD, 1);
}

/*----------------------------------------------------------------------------
 * Recover rank information on a given application number
 *
 * parameters:
 *   app_name  <-- application name
 *   root_rank --> associated root rank
 *   n_ranks   --> number of associated ranks
 *----------------------------------------------------------------------------*/

void
syr_mpi_appinfo(const char  *app_name,
                int         *root_rank,
                int         *n_ranks)
{
  int n_apps = 0;

  *root_rank = -1;
  *n_ranks = -1;

  if (syr_glob_coupling_world != NULL) {

    int i;

    n_apps = ple_coupling_mpi_set_n_apps(syr_glob_coupling_world);

    for (i = 0; i < n_apps; i++) {

      int match = 0;

      const ple_coupling_mpi_set_info_t
        ai = ple_coupling_mpi_set_get_info(syr_glob_coupling_world, i);

      /* If app name not given, assume any non-SYRTHES application
         is a CFD code (so as to allow packages based on Code_Saturne but
         using their own identities) */

      if (app_name != NULL) {
        if (!strcmp(ai.app_name, app_name))
          match = 1;
      }
      else if (strncmp(ai.app_type, "SYRTHES", 7))
        match = 1;

      if (match) {

        *root_rank = ai.root_rank;
        *n_ranks = ai.n_ranks;

        printf("  Couplage CFD:\n"
               "   type d'application :        \"%s\"\n"
               "   nom de l'instance :         \"%s\"\n"
               "   rang racine MPI :           %d\n"
               "   nombre de rangs MPI :       %d\n\n",
               ai.app_type, ai.app_name,
               ai.root_rank, ai.n_ranks);

        break;
      }
    }
  }

  if (*root_rank < 0)
    ple_error(__FILE__, __LINE__, 0,
              "Application MPI \"%s\" non trouvee.", app_name);
}

#endif /* defined (HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Exit / Stop
 *----------------------------------------------------------------------------*/

void
syr_exit(int status)
{
  if (status == EXIT_FAILURE) {

    fprintf(stdout, "\n\n %s \n\n", "SYRTHES : erreur(s) rencontree(s).");

#if defined(DEBUG) || !defined(NDEBUG)
    assert(0);
#endif

#if defined(HAVE_MPI)
    {
      int mpi_flag;

      MPI_Initialized(&mpi_flag);

      if (mpi_flag != 0)
        MPI_Abort(MPI_COMM_WORLD, status);
    }
#endif

  }

  exit(status);
}

/*----------------------------------------------------------------------------
 * Print warning
 *----------------------------------------------------------------------------*/

void
syr_warn(void)
{
  printf("\n"
         "Attention (SYRTHES)\n"
         "===================\n");
  fflush(stdout);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
