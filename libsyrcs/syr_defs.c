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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_mem.h>
#include <bft_mem_usage.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
#include <fvm_coupling.h>
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

int syr_glob_base_rank = -1;            /* Parallel rank; -1 if serial */

char syr_glob_build_date[] = __DATE__;  /* Build date */

#if defined(HAVE_MPI)
fvm_coupling_mpi_world_t *syr_glob_coupling_world = NULL;
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
 * First analysis of the command line to determine if we require MPI,
 * and initialization if necessary
 *
 * parameters:
 *   argc  <-- number of command line arguments
 *   argv  <-- array of command line arguments
 *
 * returns:
 *   -1 if MPI is not needed, or application number in MPI_COMM_WORLD of
 *   processes associated with this instance of Code_Saturne
 *----------------------------------------------------------------------------*/

static int
_mpi_app_num(int    argc,
             char  *argv[])
{
  char *s;

  int arg_id = 0, flag = 0;
  int appnum = -1;

  /* Loop on command line arguments */

  arg_id = 0;

  while (++arg_id < argc) {

    s = argv[arg_id];

    if (strcmp(s, "-app-num") == 0) {
      int _appnum = 0;
      if (arg_id + 1 < argc) {
        if (sscanf(argv[arg_id + 1], "%d", &_appnum) == 1)
          appnum = _appnum;
      }
    }

  }

  /*
    If appnum was not given through the command line but we
    are running under MPI-2, appnum may be available through
    the MPI_APPNUM attribute.
  */

#if defined(MPI_VERSION) && (MPI_VERSION >= 2)
  if (appnum < 0) {
    void *attp = NULL;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_APPNUM, &attp, &flag);
    if (flag != 0)
      appnum = *(int *)attp;
  }
#endif

  if (appnum < 0)
    appnum = 0;

  return appnum;
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
  bft_error_handler_set(_syr_error_handler);
}

/*----------------------------------------------------------------------------
 * Initialize memory management
 *----------------------------------------------------------------------------*/

void
syr_mem_initialize(void)
{
  char  *base_name;
  char  *full_name = NULL;

  /* Initialize memory counting */

  bft_mem_usage_init();

  /* Initialize memory management */

  if ((base_name = getenv("SYR_FIC_MEM")) != NULL) {

    full_name = (char *)(malloc((strlen(base_name) + 6) * sizeof (char)));

    if (full_name != NULL)
        strcpy(full_name, base_name);

  }

  bft_mem_init(full_name);

  if (full_name != NULL)
    free(full_name);
}

/*----------------------------------------------------------------------------
 * Finalize memory management
 *----------------------------------------------------------------------------*/

void
syr_mem_finalize(void)
{
  double  rval;
  int ii;

  const char  unit[] = {'k', 'm', 'g', 't', 'p'};

  /* Memory summary */

  printf("\nBilan de l'occupation memoire :\n\n");

  rval = (double) bft_mem_usage_max_pr_size();

  for (ii = 0; rval > 1024. && unit[ii] != 'p'; ii++)
    rval /= 1024.;

  /* Printout */

  printf ("  Consommation memoire totale mesuree :  %12.3f %co\n",
          rval, unit[ii]);

  /* Finalize memory managment */

  bft_mem_end();

  /* Finalize memory count */

  bft_mem_usage_end();
}

#if defined (HAVE_MPI)

/*----------------------------------------------------------------------------
 * Initialize MPI communication if necessary
 *----------------------------------------------------------------------------*/

void
syr_mpi_initialize(int    *argc,
                   char  **argv[])
{
  int mpi_init_flag, rank;
  MPI_Comm  mpi_comm_syr = MPI_COMM_NULL;
  int app_num = -1, app_num_max = -1;
  int ierror = 0;

  /* Initialize MPI */

  _mpi_test_and_initialize(argc, argv);

  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    return;

  app_num = _mpi_app_num(*argc, *argv);

  /*
    Split MPI_COMM_WORLD to separate different coupled applications
    (collective operation, like all MPI communicator creation operations).

    We suppose the color argument to MPI_Comm_split is equal to the
    application number, given through the command line or through
    mpiexec and passed here as an argument.
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Allreduce(&app_num, &app_num_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (app_num_max > 0) {
    ierror = MPI_Comm_split(MPI_COMM_WORLD, app_num, rank, &mpi_comm_syr);
    MPI_Comm_rank(mpi_comm_syr, &syr_glob_base_rank);
  }
  else
    ierror = 1;

  if (ierror != 0)
    bft_error(__FILE__, __LINE__, 0,
              "Erreur a la creation d'un communicateur local a SYRTHES.\n");

  /* Discover other applications in the same MPI root communicator
     (and participate in corresponding communication). */

  syr_glob_coupling_world = fvm_coupling_mpi_world_create(app_num,
                                                          "SYRTHES 3.4",
                                                          NULL,
                                                          mpi_comm_syr);

  /* mpi_comm_syr is not used anymore */

  if (mpi_comm_syr != MPI_COMM_NULL)
    MPI_Comm_free(&mpi_comm_syr);

  /* If more than 1 rank was assigned to SYRTHES, only 1 is active
     (this may be the case in IBM BlueGene/P, where a whole pset
     must be assigned to a given executable) */

  if (syr_glob_base_rank > 0) {
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

  fvm_coupling_mpi_world_destroy(&syr_glob_coupling_world);

  assert(syr_glob_coupling_world == NULL);

  ierror = MPI_Barrier(MPI_COMM_WORLD);

  if (ierror != 0)
    bft_error(__FILE__, __LINE__, 0,
              "Erreur dans MPI_Barrier lors de la finalisation du\n"
              "communicateur global cote SYRTHES.");

  ierror = MPI_Finalize();

  if (ierror != 0)
    bft_error(__FILE__, __LINE__, 0,
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
 *   app_num   <-- application number
 *   root_rank --> associated root rank
 *   n_ranks   --> number of associated ranks
 *----------------------------------------------------------------------------*/

void
syr_mpi_appinfo(int    app_num,
                int   *root_rank,
                int   *n_ranks)
{
  int n_apps = 0;

  *root_rank = -1;
  *n_ranks = -1;

  if (syr_glob_coupling_world != NULL) {

    int i;

    n_apps = fvm_coupling_mpi_world_n_apps(syr_glob_coupling_world);

    for (i = 0; i < n_apps; i++) {

      const fvm_coupling_mpi_world_info_t
        ai = fvm_coupling_mpi_world_get_info(syr_glob_coupling_world, i);

      if (ai.app_num == app_num) {

        *root_rank = ai.root_rank;
        *n_ranks = ai.n_ranks;

        printf("  Couplage CFD:\n"
               "    Numero d'application MPI :  %d\n"
               "   type d'application :        \"%s\"\n"
               "   nom de l'instance :         \"%s\"\n"
               "   rang racine MPI :           %d\n"
               "   nombre de rangs MPI :       %d\n\n",
               ai.app_num, ai.app_type, ai.app_name,
               ai.root_rank, ai.n_ranks);

        break;
      }
    }
  }

  if (*root_rank < 0)
    bft_error(__FILE__, __LINE__, 0,
              "Application MPI numero %d non trouvee.", app_num);
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
