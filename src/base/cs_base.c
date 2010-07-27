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

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

#if defined(HAVE_GETCWD) || defined(HAVE_SLEEP)
#include <unistd.h>
#endif

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)
#include <pwd.h>
#endif

#if defined(HAVE_UNAME)
#include <sys/utsname.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_backtrace.h>
#include <bft_mem_usage.h>
#include <bft_mem.h>
#include <bft_printf.h>
#include <bft_sys_info.h>
#include <bft_timer.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>
#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Fortran API */
/*-------------*/

/*
 * 'usual' maximum name length; a longer name is possible, but will
 * provoque a dynamic memory allocation.
 */

#define CS_BASE_N_STRINGS                               5
#define CS_BASE_STRING_LEN                             64

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

#if defined(HAVE_MPI)

typedef struct
{
  long val;
  int  rank;
} _cs_base_mpi_long_int_t;

typedef struct
{
  double val;
  int    rank;
} _cs_base_mpi_double_int_t;

#endif

/* Type to backup signal handlers */

typedef void (*_cs_base_sighandler_t) (int);

/*============================================================================
 *  Global variables
 *============================================================================*/

int  cs_glob_n_threads = 1;    /* Number of threads */

int  cs_glob_rank_id = -1;     /* Rank of process in communicator */
int  cs_glob_n_ranks =  1;     /* Number of processes in communicator */

#if defined(HAVE_MPI)
MPI_Comm  cs_glob_mpi_comm = MPI_COMM_NULL;   /* Intra-communicator */
#endif

static bft_error_handler_t  *cs_glob_base_err_handler_save = NULL;

/* Static (private) global variables */
/*-----------------------------------*/

static cs_bool_t  cs_glob_base_bft_mem_init = false;

static cs_bool_t  cs_glob_base_str_init = false;
static cs_bool_t  cs_glob_base_str_is_free[CS_BASE_N_STRINGS];
static char       cs_glob_base_str[CS_BASE_N_STRINGS]
                                     [CS_BASE_STRING_LEN + 1];

/* Global variables associated with signal handling */

#if defined(SIGHUP)
static _cs_base_sighandler_t cs_glob_base_sighup_save = SIG_DFL;
#endif

static _cs_base_sighandler_t cs_glob_base_sigint_save = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigterm_save = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigfpe_save = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigsegv_save = SIG_DFL;

#if defined(SIGXCPU)
static _cs_base_sighandler_t cs_glob_base_sigcpu_save = SIG_DFL;
#endif

/* Global variables for handling of Fortran work arrays */

static char      _cs_glob_srt_ia_peak[7] = "------";
static char      _cs_glob_srt_ra_peak[7] = "------";
static cs_int_t  _cs_glob_mem_ia_peak = 0;
static cs_int_t  _cs_glob_mem_ra_peak = 0;
static cs_int_t  _cs_glob_mem_ia_size = 0;
static cs_int_t  _cs_glob_mem_ra_size = 0;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * False print of a message to standard output for discarded logs
 *----------------------------------------------------------------------------*/

static int
_cs_base_bft_printf_null(const char  *format,
                         va_list      arg_ptr)
{
  return 0;
}

/*----------------------------------------------------------------------------
 * Flush of log output buffer
 *----------------------------------------------------------------------------*/

static int
_cs_base_bft_printf_flush(void)
{
  return fflush(stdout);
}

/*----------------------------------------------------------------------------
 * False flush of log output buffer for discarded logs
 *----------------------------------------------------------------------------*/

static int
_cs_base_bft_printf_flush_null(void)
{
  return 0;
}

/*----------------------------------------------------------------------------
 * Print a message to the error output
 *
 * The message is repeated on the standard output and an error file.
 *----------------------------------------------------------------------------*/

static void
_cs_base_err_vprintf(const char  *format,
                     va_list      arg_ptr)
{
  static cs_bool_t  initialized = false;

  /* message to the standard output */

#if defined(va_copy) || defined(__va_copy)
  {
    va_list arg_ptr_2;
    bft_printf_proxy_t  *_bft_printf_proxy = bft_printf_proxy_get();

#if defined(va_copy)
    va_copy(arg_ptr_2, arg_ptr);
#else
    __va_copy(arg_ptr_2, arg_ptr);
#endif
    _bft_printf_proxy(format, arg_ptr_2);
    va_end(arg_ptr_2);
  }
#endif

  /* message on a specific error output, initialized only if the
     error output is really necessary */

  if (initialized == false) {

    char err_file_name[81];

    if (cs_glob_rank_id < 1)
      strcpy(err_file_name, "error");

    else {
#if defined(HAVE_SLEEP)
      /* Wait a few seconds, so that if rank 0 also has encountered an error,
         it may kill other ranks through MPI_Abort, so that only rank 0 will
         generate an error file. If rank 0 has not encountered the error,
         proceed normally after the wait.
         As sleep() may be interrupted by a signal, repeat as long as the wait
         time is not elapsed; */
      int wait_time = (cs_glob_n_ranks < 64) ? 1: 10;
      double stime = bft_timer_wtime();
      double etime = 0.0;
      do {
        sleep(wait_time);
        etime = bft_timer_wtime();
      }
      while (etime > -0.5 && etime - stime < wait_time); /* etime = -1 only if
                                                            bft_timer_wtime()
                                                            is unusable. */
#endif
      if (cs_glob_n_ranks > 9999)
        sprintf(err_file_name, "error_n%07d", cs_glob_rank_id + 1);
      else
        sprintf(err_file_name, "error_n%04d", cs_glob_rank_id + 1);
    }

    freopen(err_file_name, "w", stderr);

    initialized = true;
  }

  vfprintf(stderr, format, arg_ptr);
}

/*----------------------------------------------------------------------------
 * Print a message to error output
 *
 * The message is repeated on the standard output and an error file.
 *----------------------------------------------------------------------------*/

static void
_cs_base_err_printf(const char  *format,
                    ...)
{
  /* Initialisation de la liste des arguments */

  va_list  arg_ptr;

  va_start(arg_ptr, format);

  /* message sur les sorties */

  _cs_base_err_vprintf(format, arg_ptr);

  /* Finalisation de la liste des arguments */

  va_end(arg_ptr);
}

/*----------------------------------------------------------------------------
 * Exit function
 *----------------------------------------------------------------------------*/

static void
_cs_base_exit(int status)
{
#if defined(HAVE_MPI)
  {
    int mpi_flag;

    MPI_Initialized(&mpi_flag);

    if (mpi_flag != 0) {

      if (status != EXIT_SUCCESS)
        MPI_Abort(cs_glob_mpi_comm, EXIT_FAILURE);

      else { /*  if (status == EXIT_SUCCESS) */

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();

      }
    }
  }
#endif /* HAVE_MPI */

  exit(status);
}

/*----------------------------------------------------------------------------
 * Stop the code in case of error
 *----------------------------------------------------------------------------*/

static void
_cs_base_error_handler(const char  *nom_fic,
                       int          num_ligne,
                       int          code_err_sys,
                       const char  *format,
                       va_list      arg_ptr)
{
  bft_printf_flush();

  _cs_base_err_printf("\n");

  if (code_err_sys != 0)
    _cs_base_err_printf(_("\nSystem error: %s\n"), strerror(code_err_sys));

  _cs_base_err_printf(_("\n%s:%d: Fatal error.\n\n"), nom_fic, num_ligne);

  _cs_base_err_vprintf(format, arg_ptr);

  _cs_base_err_printf("\n\n");

  bft_backtrace_print(3);

  _cs_base_exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------
 * Print memory usage summary in case of error
 *----------------------------------------------------------------------------*/

static void
_error_mem_summary(void)
{
  size_t mem_usage;

  _cs_base_err_printf(_("\n\n"
                        "Memory allocation summary\n"
                        "-------------------------\n\n"));

  /* Available memory usage information */

  _cs_base_err_printf
    (_("Theoretical current allocated memory:   %llu kB\n"),
     (unsigned long long)(bft_mem_size_current()));

  _cs_base_err_printf
    (_("Theoretical maximum allocated memory:   %llu kB\n"),
     (unsigned long long)(bft_mem_size_max()));

  if (bft_mem_usage_initialized() == 1) {

    /* Maximum measured memory */

    mem_usage = bft_mem_usage_max_pr_size();
    if (mem_usage > 0)
      _cs_base_err_printf
        (_("Maximum program memory measure:         %llu kB\n"),
         (unsigned long long)mem_usage);

    /* Current measured memory */

    mem_usage = bft_mem_usage_pr_size();
    if (mem_usage > 0)
      _cs_base_err_printf
        (_("Current program memory measure:         %llu kB\n"),
         (unsigned long long)mem_usage);
  }
}

/*----------------------------------------------------------------------------
 * Memory allocation error handler.
 *
 * Memory status is written to the error output, and the general error
 * handler used by bft_error() is called (which results in the termination
 * of the current process).
 *
 * parameters:
 *   file_name      <-- name of source file from which error handler called.
 *   line_num       <-- line of source file from which error handler called.
 *   sys_error_code <-- error code if error in system or libc call, 0 otherwise.
 *   format         <-- format string, as printf() and family.
 *   arg_ptr        <-> variable argument list based on format string.
 *----------------------------------------------------------------------------*/

static void
_cs_mem_error_handler(const char  *file_name,
                      int          line_num,
                      int          sys_error_code,
                      const char  *format,
                      va_list      arg_ptr)
{
  bft_error_handler_t * general_err_handler;

  _error_mem_summary();

  general_err_handler = bft_error_handler_get();
  general_err_handler(file_name, line_num, sys_error_code, format, arg_ptr);
}

/*----------------------------------------------------------------------------
 * Print a stack trace
 *----------------------------------------------------------------------------*/

static void
_cs_base_backtrace_print(int  niv_debut)
{
  size_t  ind;
  bft_backtrace_t  *tr = NULL;

  tr = bft_backtrace_create();

  if (tr != NULL) {

    char s_func_buf[67];

    const char *s_file;
    const char *s_func;
    const char *s_addr;

    const char s_inconnu[] = "?";
    const char s_vide[] = "";
    const char *s_prefix = s_vide;

    size_t nbr = bft_backtrace_size(tr);

    if (nbr > 0)
      _cs_base_err_printf(_("\nCall stack:\n"));

    for (ind = niv_debut ; ind < nbr ; ind++) {

      s_file = bft_backtrace_file(tr, ind);
      s_func = bft_backtrace_function(tr, ind);
      s_addr = bft_backtrace_address(tr, ind);

      if (s_file == NULL)
        s_file = s_inconnu;
      if (s_func == NULL)
        strcpy(s_func_buf, "?");
      else {
        s_func_buf[0] = '<';
        strncpy(s_func_buf + 1, s_func, 64);
        strcat(s_func_buf, ">");
      }
      if (s_addr == NULL)
        s_addr = s_inconnu;

      _cs_base_err_printf("%s%4d: %-12s %-32s (%s)\n", s_prefix,
                          ind-niv_debut+1, s_addr, s_func_buf, s_file);

    }

    bft_backtrace_destroy(tr);

    if (nbr > 0)
      _cs_base_err_printf(_("End of stack\n\n"));
  }

}

/*----------------------------------------------------------------------------
 * Handle a fatal signal (such as SIGFPE or SIGSEGV)
 *----------------------------------------------------------------------------*/

static void
_cs_base_sig_fatal(int  signum)
{
  bft_printf_flush();

  switch (signum) {

#if defined(SIGHUP)
  case SIGHUP:
    _cs_base_err_printf(_("SIGHUP signal (hang-up) intercepted.\n"
                          "--> computation interrupted.\n"));
    break;
#endif

  case SIGINT:
    _cs_base_err_printf(_("SIGINT signal (Control+C or equivalent) received.\n"
                          "--> computation interrupted by user.\n"));
    break;

  case SIGTERM:
    _cs_base_err_printf(_("SIGTERM signal (termination) received.\n"
                          "--> computation interrupted by environment.\n"));
    break;

  case SIGFPE:
    _cs_base_err_printf(_("SIGFPE signal (floating point exception) "
                          "intercepted!\n"));
    break;

  case SIGSEGV:
    _cs_base_err_printf(_("SIGSEGV signal (forbidden memory area access) "
                          "intercepted!\n"));
    break;

#if defined(SIGXCPU)
  case SIGXCPU:
    _cs_base_err_printf(_("SIGXCPU signal (CPU time limit reached) "
                          "intercepted.\n"));
    break;
#endif

  default:
    _cs_base_err_printf(_("Signal %d intercepted!\n"), signum);
  }

  bft_backtrace_print(3);

  _cs_base_exit(EXIT_FAILURE);
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 *  Finalisation MPI
 *----------------------------------------------------------------------------*/

static void
_cs_base_mpi_fin(void)
{
#if defined(HAVE_MPI)
  fvm_parall_set_mpi_comm(MPI_COMM_NULL);
#endif

  bft_error_handler_set(cs_glob_base_err_handler_save);
  ple_error_handler_set(cs_glob_base_err_handler_save);

  if (   cs_glob_mpi_comm != MPI_COMM_NULL
      && cs_glob_mpi_comm != MPI_COMM_WORLD)
    MPI_Comm_free(&cs_glob_mpi_comm);
}


#if defined(DEBUG) || !defined(NDEBUG)

/*----------------------------------------------------------------------------
 * MPI error handler
 *----------------------------------------------------------------------------*/

static void
_cs_base_erreur_mpi(MPI_Comm  *comm,
                    int       *errcode,
                    ...)
{
  int err_len;
  char err_string[MPI_MAX_ERROR_STRING + 1];

#if defined MPI_MAX_OBJECT_NAME
  int name_len = 0;
  char comm_name[MPI_MAX_OBJECT_NAME + 1];
#endif

  bft_printf_flush();

  _cs_base_err_printf("\n");

  MPI_Error_string(*errcode, err_string, &err_len);
  err_string[err_len] = '\0';

#if defined MPI_MAX_OBJECT_NAME
  MPI_Comm_get_name(*comm, comm_name, &name_len);
  comm_name[name_len] = '\0';
  _cs_base_err_printf(_("\nMPI error (communicator %s):\n"
                        "%s\n"), comm_name, err_string);
#else
  _cs_base_err_printf(_("\nMPI error:\n"
                        "%s\n"), err_string);
#endif

  _cs_base_err_printf("\n\n");

  bft_backtrace_print(3);

  _cs_base_exit(EXIT_FAILURE);
}

#endif

/*----------------------------------------------------------------------------
 * Complete MPI setup.
 *
 * MPI should have been initialized by cs_base_mpi_init().
 *
 * The application name is used to build subgroups of processes with
 * identical name from the MPI_COMM_WORLD communicator, thus separating
 * this instance of Code_Saturne from other coupled codes. It may be
 * defined using the --app-num argument, and is based on the working
 * directory's base name otherwise.
 *
 * parameters:
 *   app_name <-- pointer to application instance name.
 *
 * returns:
 *   application number (-1 without coupled codes, >= 0 if coupled codes)
 *----------------------------------------------------------------------------*/

static int
_cs_base_mpi_setup(const char *app_name)
{
  int nbr, rank;

  int app_num = -1;

#if defined(DEBUG) || !defined(NDEBUG)
  MPI_Errhandler errhandler;
#endif

  app_num = ple_coupling_mpi_name_to_id(MPI_COMM_WORLD, app_name);

  /*
    Split MPI_COMM_WORLD to separate different coupled applications
    (collective operation, like all MPI communicator creation operations).

    app_num is equal to -1 if all applications have the same instance
    name, in which case no communicator split is necessary.
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (app_num > -1)
    MPI_Comm_split(MPI_COMM_WORLD, app_num, rank, &cs_glob_mpi_comm);
  else
    cs_glob_mpi_comm = MPI_COMM_WORLD;

  MPI_Comm_size(cs_glob_mpi_comm, &nbr);
  MPI_Comm_rank(cs_glob_mpi_comm, &rank);

  cs_glob_n_ranks = nbr;

  if (cs_glob_n_ranks > 1)
    cs_glob_rank_id = rank;

  /* cs_glob_mpi_comm may not be freed at this stage, as it
     it may be needed to build intercommunicators for couplings,
     but we may set cs_glob_rank_id to its serial value if
     we are only using MPI for coupling. */

  if (cs_glob_n_ranks == 1 && app_num > -1)
    cs_glob_rank_id = -1;

  /* Initialize associated libraries */

#if defined(HAVE_MPI)
  if (cs_glob_rank_id > -1)
    fvm_parall_set_mpi_comm(cs_glob_mpi_comm);
  else
    fvm_parall_set_mpi_comm(MPI_COMM_NULL);
#endif

#if defined(DEBUG) || !defined(NDEBUG)
  if (nbr > 1 || cs_glob_mpi_comm != MPI_COMM_NULL) {
    MPI_Errhandler_create(&_cs_base_erreur_mpi, &errhandler);
    MPI_Errhandler_set(MPI_COMM_WORLD, errhandler);
    if (   cs_glob_mpi_comm != MPI_COMM_WORLD
        && cs_glob_mpi_comm != MPI_COMM_NULL)
      MPI_Errhandler_set(cs_glob_mpi_comm, errhandler);
    MPI_Errhandler_free(&errhandler);
  }
#endif

  return app_num;
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Maximum value of a counter and associated 6 character string
 * (used for Fortan maximum memory count in IA/RA).
 *
 * parameters:
 *   mem_peak  <-> maximum value reached
 *   srt_peak  <-> associated subroutine name
 *----------------------------------------------------------------------------*/

static void
_cs_base_work_mem_max(cs_int_t  *mem_peak,
                      char       srt_peak[6])
{
#if defined(HAVE_MPI)

  _cs_base_mpi_long_int_t  val_in, val_max;

  assert(sizeof(double) == sizeof(cs_real_t));

  val_in.val  = *mem_peak;
  val_in.rank = cs_glob_rank_id;

  MPI_Allreduce(&val_in, &val_max, 1, MPI_LONG_INT, MPI_MAXLOC,
                cs_glob_mpi_comm);

  *mem_peak = val_max.val;

  MPI_Bcast(srt_peak, 6, MPI_CHAR, val_max.rank, cs_glob_mpi_comm);
#endif
}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Call exit routine from Fortran code
 *
 * Fortran interface:
 *
 * SUBROUTINE CSEXIT (STATUS)
 * *****************
 *
 * INTEGER          STATUS      : --> : 0 for success, 1+ for error
 *----------------------------------------------------------------------------*/

void CS_PROCF (csexit, CSEXIT)
(
  const cs_int_t  *status
)
{
  cs_exit (*status);
}

/*----------------------------------------------------------------------------
 * CPU time used since execution start
 *
 * Fortran interface:
 *
 * SUBROUTINE DMTMPS (TCPU)
 * *****************
 *
 * DOUBLE PRECISION TCPU        : --> : CPU time (user + system)
 *----------------------------------------------------------------------------*/

void CS_PROCF (dmtmps, DMTMPS)
(
  cs_real_t  *tcpu
)
{
  *tcpu = bft_timer_cpu_time();
}

/*----------------------------------------------------------------------------
 * Check that main integer working array memory reservation fits within
 * the allocated size of IA.
 *
 * Fortran interface:
 *
 * SUBROUTINE IASIZE (CALLER, MEMINT)
 * *****************
 *
 * CHARACTER*6      CALLER      : --> : Name of calling subroutine
 * INTEGER          MEMINT      : --> : Last required element in IA
 *----------------------------------------------------------------------------*/

void CS_PROCF (iasize, IASIZE)
(
 const char   caller[6],
 cs_int_t    *memint
)
{
  /* Test if enough memory is available */

  if (*memint > _cs_glob_mem_ia_size) {
    char _caller[7];
    strncpy(_caller, caller, 6);
    _caller[6] = '\0';
    bft_error
      (__FILE__, __LINE__, 0,
       _(" Sub-routine calling iasize:                %s\n"
         " Memory needed in ia (number of integers):  %d\n"
         "         available:                         %d\n\n"
         " ----> Define iasize to a value at least equal to %d integers)."),
       _caller, *memint, _cs_glob_mem_ia_size, *memint);
  }

  /* Update _cs_glob_mem_ia_peak and _cs_glob_srt_ia_peak */

  else if (*memint > _cs_glob_mem_ia_peak) {

    _cs_glob_mem_ia_peak = *memint;
    strncpy(_cs_glob_srt_ia_peak, caller, 6);
    _cs_glob_srt_ia_peak[6] = '\0';
  }
}

/*----------------------------------------------------------------------------
 * Check that main floating-point working array memory reservation fits
 * within the allocated size of RA.
 *
 * Fortran interface:
 *
 * SUBROUTINE RASIZE (CALLER, MEMINT)
 * *****************
 *
 * CHARACTER*6      CALLER      : --> : Name of calling subroutine
 * INTEGER          MEMRDP      : --> : Last required element in RA
 *----------------------------------------------------------------------------*/

void CS_PROCF (rasize, RASIZE)
(
 const char   caller[6],
 cs_int_t    *memrdp
)
{
  /* Test if enough memory is available */

  if (*memrdp > _cs_glob_mem_ra_size) {
    char _caller[7];
    strncpy(_caller, caller, 6);
    _caller[6] = '\0';
    bft_error
      (__FILE__, __LINE__, 0,
       _(" Sub-routine calling rasize:             %s\n"
         " Memory needed in ra (number of reals):  %d\n"
         "         available:                      %d\n\n"
         " ----> Define rasize to a value at least equal to %d reals)."),
       _caller, *memrdp, _cs_glob_mem_ra_size, *memrdp);
  }

  /* Update _cs_glob_mem_ra_peak and _cs_glob_srt_ra_peak */

  else if (*memrdp > _cs_glob_mem_ra_peak) {

    _cs_glob_mem_ra_peak = *memrdp;
    strncpy(_cs_glob_srt_ra_peak, caller, 6);
    _cs_glob_srt_ra_peak[6] = '\0';
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

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

char *
cs_base_get_app_name(int          argc,
                     const char  *argv[])
{
  char *app_name = NULL;
  int arg_id = 0;

  /* Loop on command line arguments */

  arg_id = 0;

  while (++arg_id < argc) {
    const char *s = argv[arg_id];
    if (strcmp(s, "--app-name") == 0) {
      if (arg_id + 1 < argc) {
        BFT_MALLOC(app_name, strlen(argv[arg_id + 1]) + 1, char);
        strcpy(app_name, argv[arg_id + 1]);
      }
    }
  }

  /* Use execution directory if name is unavailable */

#if defined(HAVE_GETCWD)

  if (app_name == NULL) {

    int i;
    int buf_size = 128;
    char *wd = NULL, *buf = NULL;

    while (wd == NULL) {
      buf_size *= 2;
      BFT_REALLOC(buf, buf_size, char);
      wd = getcwd(buf, buf_size);
      if (wd == NULL && errno != ERANGE)
        bft_error(__FILE__, __LINE__, errno,
                  _("Error querying working directory.\n"));
    }

    for (i = strlen(buf) - 1; i > 0 && buf[i-1] != '/'; i--);
    BFT_MALLOC(app_name, strlen(buf + i) + 1, char);
    strcpy(app_name, buf + i);
    BFT_FREE(buf);
  }

#endif /* defined(HAVE_GETCWD) */

  return app_name;
}

/*----------------------------------------------------------------------------
 * Print logfile header
 *
 * parameters:
 *   argc  <-- number of command line arguments
 *   argv  <-- array of command line arguments
 *----------------------------------------------------------------------------*/

void
cs_base_logfile_head(int    argc,
                     char  *argv[])
{
  char str[81];
  int ii;
  char date_str[] = __DATE__;
  char time_str[] = __TIME__;
  const char mon_name[12][4]
    = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
  struct tm time_cnv;

  /* Define MPI Information */

#if defined(MPI_VERSION) && defined(MPI_SUBVERSION)
#if defined(OPEN_MPI)
  const char mpi_lib[] = "Open MPI";
#elif defined(MPICH2)
  const char mpi_lib[] = "MPICH2";
#elif defined(LAM_MPI)
  const char mpi_lib[] = "LAM/MPI";
#elif defined(MPICH_NAME)
  const char mpi_lib[] = "MPICH";
#elif defined(HP_MPI)
  const char mpi_lib[] = "HP-MPI";
#elif defined(MPI_VERSION) && defined(MPI_SUBVERSION)
  const char *mpi_lib = NULL;
#endif
#endif /* defined(MPI_VERSION) && defined(MPI_SUBVERSION) */

  /* Determine compilation date */

  for (ii = 0; ii < 12; ii++) {
    if (strncmp(date_str, mon_name[ii], 3) == 0) {
      time_cnv.tm_mon = ii ;
      break;
    }
  }

  sscanf(date_str + 3, "%d", &(time_cnv.tm_mday)) ;
  sscanf(date_str + 6, "%d", &(time_cnv.tm_year)) ;

  time_cnv.tm_year -= 1900 ;

  sscanf(time_str    , "%d", &(time_cnv.tm_hour)) ;
  sscanf(time_str + 3, "%d", &(time_cnv.tm_min)) ;
  sscanf(time_str + 6, "%d", &(time_cnv.tm_sec)) ;

  time_cnv.tm_isdst = -1 ;

  /* Re-compute and internationalize build date */

  mktime(&time_cnv) ;
  strftime(str, 80, "%c", &time_cnv) ;

  /* Now print info */

  bft_printf(_("command: \n"));

  for (ii = 0 ; ii < argc ; ii++)
    bft_printf(" %s", argv[ii]);

  bft_printf("\n");
  bft_printf("\n************************************"
             "***************************\n\n");
  bft_printf("                                  (R)\n"
             "                      Code_Saturne\n\n"
             "                      Version %s\n\n",
             CS_APP_VERSION);

  bft_printf("\n  Copyright (C) 1998-2010 EDF S.A., France\n\n");

  bft_printf(_("  build %s\n"), str);

#if defined(MPI_VERSION) && defined(MPI_SUBVERSION)
  if (mpi_lib != NULL)
    bft_printf(_("  MPI version %d.%d (%s)\n\n"),
               MPI_VERSION, MPI_SUBVERSION, mpi_lib);
  else
    bft_printf(_("  MPI version %d.%d\n\n"),
               MPI_VERSION, MPI_SUBVERSION);
#endif

  bft_printf("\n");
  bft_printf("  The Code_Saturne CFD tool  is free software;\n"
             "  you can redistribute it and/or modify it under the terms\n"
             "  of the GNU General Public License as published by the\n"
             "  Free Software Foundation; either version 2 of the License,\n"
             "  or (at your option) any later version.\n\n");

  bft_printf("  The Code_Saturne CFD tool is distributed in the hope that\n"
             "  it will be useful, but WITHOUT ANY WARRANTY; without even\n"
             "  the implied warranty of MERCHANTABILITY or FITNESS FOR A\n"
             "  PARTICULAR PURPOSE.  See the GNU General Public License\n"
             "  for more details.\n");

  bft_printf("\n************************************"
             "***************************\n\n");
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * First analysis of the command line and environment variables to determine
 * if we require MPI, and initialization if necessary.
 *
 * parameters:
 *   argc  <-> number of command line arguments
 *   argv  <-> array of command line arguments
 *
 * Global variables `cs_glob_n_ranks' (number of Code_Saturne processes)
 * and `cs_glob_rank_id' (rank of local process) are set by this function.
 *
 * returns:
 *   -1 if MPI is not or no other code is coupled, or application number
 *   in MPI_COMM_WORLD of processes associated with this instance of
 *   Code_Saturne
 *----------------------------------------------------------------------------*/

int
cs_base_mpi_init(int    *argc,
                 char  **argv[])
{
#if defined(HAVE_MPI)

  char *s;

  int arg_id = 0, flag = 0;
  int use_mpi = false;

  int app_num = -1;

#if   defined(__bg__) || defined(__CRAYXT_COMPUTE_LINUX_TARGET)

  /* Notes: Blue Gene/L also defines the BGLMPI_SIZE environment variable.
   *        Blue Gene/P defines BG_SIZE (plus BG_MAPPING, and BG_RELEASE). */

  use_mpi = true;

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

  for (arg_id = 0 ; arg_id < *argc ; arg_id++) {
    if (   !strcmp((*argv)[arg_id], "-p4pg")         /* For process 0 */
        || !strcmp((*argv)[arg_id], "-p4rmrank")) {  /* For other processes */
      use_mpi = true;
      break;
    }
  }

  if (getenv("GMPI_ID") != NULL) /* In case we are using MPICH-GM */
    use_mpi = true;

#elif defined(LAM_MPI)
  if (getenv("LAMRANK") != NULL)
    use_mpi = true;

#elif defined(OPEN_MPI)
  if (getenv("OMPI_MCA_ns_nds_vpid") != NULL)         /* OpenMPI 1.2 */
    use_mpi = true;
  else if (getenv("OMPI_COMM_WORLD_RANK") != NULL)    /* OpenMPI 1.3 */
    use_mpi = true;

#elif defined(MPICH2)
  if (getenv("PMI_RANK") != NULL)
    use_mpi = true;

#endif /* Tests for known MPI variants */

  /* If we have determined from known MPI environment variables
     of command line arguments that we are running under MPI,
     initialize MPI */

  if (use_mpi == true) {
    MPI_Initialized(&flag);
    if (!flag) {
#if defined(MPI_VERSION) && (MPI_VERSION >= 2) && defined(HAVE_OPENMP)
      int mpi_threads;
      MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &mpi_threads);
#else
      MPI_Init(argc, argv);
#endif
    }
  }

  /* Loop on command line arguments */

  arg_id = 0;

  while (++arg_id < *argc) {

    s = (*argv)[arg_id];

    /* Parallel run */

    if (strcmp(s, "--mpi") == 0)
      use_mpi = true;

  } /* End of loop on command line arguments */

  if (use_mpi == true) {

    MPI_Initialized(&flag);
    if (!flag) {
#if defined(MPI_VERSION) && (MPI_VERSION >= 2) && defined(HAVE_OPENMP)
      int mpi_threads;
      MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &mpi_threads);
#else
      MPI_Init(argc, argv);
#endif
    }

  }

  /* Now setup global variables and communicators */

  if (use_mpi == true) {

    char *app_name = cs_base_get_app_name(*argc, (const char **)(*argv));

    app_num = _cs_base_mpi_setup(app_name);

    BFT_FREE(app_name);
  }

  return app_num;

#else /* if defined(HAVE_MPI) */

  return -1;

#endif
}

#endif /* HAVE_MPI */

/*----------------------------------------------------------------------------
 * Exit, with handling for both normal and error cases.
 *
 * Finalize MPI if necessary.
 *
 * parameters:
 *   status <-- value to be returned to the parent:
 *              EXIT_SUCCESS / 0 for the normal case,
 *              EXIT_FAILURE or other nonzero code for error cases.
 *----------------------------------------------------------------------------*/

void
cs_exit(int  status)
{
  if (status == EXIT_FAILURE) {

    bft_printf_flush();
    bft_backtrace_print(2);

  }

#if defined(HAVE_MPI)

  {
    int mpi_flag;

    MPI_Initialized(&mpi_flag);

    if (mpi_flag != 0) {

      if (status != EXIT_FAILURE) {
        _cs_base_mpi_fin();
      }
    }
  }

#endif /* HAVE_MPI */

  _cs_base_exit(status);
}

/*----------------------------------------------------------------------------
 * Initialize error and signal handlers.
 *----------------------------------------------------------------------------*/

void
cs_base_error_init(void)
{
  /* Error handler */

  cs_glob_base_err_handler_save = bft_error_handler_get();
  bft_error_handler_set(_cs_base_error_handler);
  ple_error_handler_set(_cs_base_error_handler);

  /* Signal handlers */

  bft_backtrace_print_set(_cs_base_backtrace_print);

#if defined(SIGHUP)
  if (cs_glob_rank_id <= 0)
    cs_glob_base_sighup_save  = signal(SIGHUP, _cs_base_sig_fatal);
#endif

  if (cs_glob_rank_id <= 0) {
    cs_glob_base_sigint_save  = signal(SIGINT, _cs_base_sig_fatal);
    cs_glob_base_sigterm_save = signal(SIGTERM, _cs_base_sig_fatal);
  }

  cs_glob_base_sigfpe_save  = signal(SIGFPE, _cs_base_sig_fatal);
  cs_glob_base_sigsegv_save = signal(SIGSEGV, _cs_base_sig_fatal);

#if defined(SIGXCPU)
  if (cs_glob_rank_id <= 0)
    cs_glob_base_sigcpu_save = signal(SIGXCPU, _cs_base_sig_fatal);
#endif
}

/*----------------------------------------------------------------------------
 * Initialize management of memory allocated through BFT.
 *----------------------------------------------------------------------------*/

void
cs_base_mem_init(void)
{
  char  *nom_base;
  char  *nom_complet = NULL;

  /* Set error handler */

  bft_mem_error_handler_set(_cs_mem_error_handler);

  /* Set PLE library memory handler */

  ple_mem_functions_set(bft_mem_malloc,
                        bft_mem_realloc,
                        bft_mem_free);

  /* Memory usage measure initialization */

  bft_mem_usage_init();

  /* Memory management initialization */

  if ((nom_base = getenv("CS_MEM_LOG")) != NULL) {

    /* We may not use BFT_MALLOC here as memory management has
       not yet been initialized using bft_mem_init() */

    nom_complet = malloc((strlen(nom_base) + 6) * sizeof (char));

    if (nom_complet != NULL) {

      /* In parallel, we will have one trace file per MPI process */
      if (cs_glob_rank_id >= 0)
        sprintf(nom_complet, "%s.%04d", nom_base, cs_glob_rank_id + 1);
      else
        strcpy(nom_complet, nom_base);

    }

  }

  if (bft_mem_initialized())
    cs_glob_base_bft_mem_init = false;

  else {
    cs_glob_base_bft_mem_init = true;
    bft_mem_init(nom_complet);
  }

  if (nom_complet != NULL)
    free (nom_complet);
}

/*----------------------------------------------------------------------------
 * Allocate Fortran work arrays and prepare for their use.
 *
 * parameters:
 *   iasize <-- integer working array size (maximum number of values)
 *   rasize <-- floating-point working array size (maximum number of values)
 *   ia     --> pointer to integer working array
 *   ra     --> pointer to floating-point working array
 *----------------------------------------------------------------------------*/

void
cs_base_mem_init_work(size_t       iasize,
                      size_t       rasize,
                      cs_int_t   **ia,
                      cs_real_t  **ra)
{
  /* Allocate work arrays */

  BFT_MALLOC(*ia, iasize, cs_int_t);
  BFT_MALLOC(*ra, rasize, cs_real_t);

  _cs_glob_mem_ia_size = iasize;
  _cs_glob_mem_ra_size = rasize;
}

/*----------------------------------------------------------------------------
 * Finalize management of memory allocated through BFT.
 *
 * A summary of the consumed memory is given.
 *----------------------------------------------------------------------------*/

void
cs_base_mem_finalize(void)
{
  int    ind_bil, itot;
  double valreal[2];

#if defined(HAVE_MPI)
  int  imax = 0, imin = 0;
  double val_somme[2];
  int  ind_min[2];
  _cs_base_mpi_double_int_t  val_in[2], val_min[2], val_max[2];
#endif

  int   ind_val[2] = {1, 1};
  char  unite[]    = {'k', 'm', 'g', 't', 'p'};

  const char  * type_bil[] = {N_("Total memory used:                       "),
                              N_("Theoretical instrumented dynamic memory: ")};

  /* Memory summary */

  bft_printf(_("\nMemory use summary:\n\n"));

  valreal[0] = (double)bft_mem_usage_max_pr_size();
  valreal[1] = (double)bft_mem_size_max();

  /* Ignore inconsistent measurements */

  for (ind_bil = 0; ind_bil < 2; ind_bil++) {
    if (valreal[ind_bil] < 1.0)
      ind_val[ind_bil] = 0;
  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    MPI_Reduce(ind_val, ind_min, 2, MPI_INT, MPI_MIN,
               0, cs_glob_mpi_comm);
    MPI_Reduce(valreal, val_somme, 2, MPI_DOUBLE, MPI_SUM,
               0, cs_glob_mpi_comm);
    for (ind_bil = 0; ind_bil < 2; ind_bil++) {
      val_in[ind_bil].val = valreal[ind_bil];
      val_in[ind_bil].rank = cs_glob_rank_id;
    }
    MPI_Reduce(&val_in, &val_min, 2, MPI_DOUBLE_INT, MPI_MINLOC,
               0, cs_glob_mpi_comm);
    MPI_Reduce(&val_in, &val_max, 2, MPI_DOUBLE_INT, MPI_MAXLOC,
               0, cs_glob_mpi_comm);
    if (cs_glob_rank_id == 0) {
      for (ind_bil = 0; ind_bil < 2; ind_bil++) {
        ind_val[ind_bil]  = ind_min[ind_bil];
        valreal[ind_bil] = val_somme[ind_bil];
      }
    }
  }
#endif

  /* Similar handling of several instrumentation methods */

  for (ind_bil = 0 ; ind_bil < 2 ; ind_bil++) {

    /* If an instrumentation method returns an apparently consistent
       result, print it. */

    if (ind_val[ind_bil] == 1) {

      for (itot = 0;
           valreal[ind_bil] > 1024. && unite[itot] != 'p';
           itot++)
        valreal[ind_bil] /= 1024.;
#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1 && cs_glob_rank_id == 0) {
        for (imin = 0;
             val_min[ind_bil].val > 1024. && unite[imin] != 'p';
             imin++)
          val_min[ind_bil].val /= 1024.;
        for (imax = 0;
             val_max[ind_bil].val > 1024. && unite[imax] != 'p';
             imax++)
          val_max[ind_bil].val /= 1024.;
      }
#endif

      /* Print to log file */

      bft_printf(_("  %s %12.3f %cb\n"),
                 _(type_bil[ind_bil]), valreal[ind_bil], unite[itot]);

#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1 && cs_glob_rank_id == 0) {
        bft_printf(_("                             "
                     "local minimum: %12.3f %cb  (rank %d)\n"),
                   val_min[ind_bil].val, unite[imin], val_min[ind_bil].rank);
        bft_printf(_("                             "
                     "local maximum: %12.3f %cb  (rank %d)\n"),
                   val_max[ind_bil].val, unite[imax], val_max[ind_bil].rank);
      }
#endif
    }

  }

  /* Information on Fortran working arrays */

  if (cs_glob_n_ranks > 1) {
    _cs_base_work_mem_max(&_cs_glob_mem_ia_peak, _cs_glob_srt_ia_peak);
    _cs_base_work_mem_max(&_cs_glob_mem_ra_peak, _cs_glob_srt_ra_peak);
  }

  if (_cs_glob_mem_ia_size > 0 || _cs_glob_mem_ra_size > 0) {

    size_t wk_unit[2] = {0, 0};
    double wk_size[2] = {0., 0.};

    wk_size[0] = (  sizeof(cs_int_t)*_cs_glob_mem_ia_size
                  + sizeof(cs_real_t)*_cs_glob_mem_ra_size) / 1000;
    wk_size[1] = (  sizeof(cs_int_t)*_cs_glob_mem_ia_peak
                  + sizeof(cs_real_t)*_cs_glob_mem_ra_peak) / 1000;

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {
      double _wk_size_loc = wk_size[0];
      MPI_Allreduce(&_wk_size_loc, &(wk_size[0]), 1, MPI_DOUBLE, MPI_MAX,
                    cs_glob_mpi_comm);
    }
#endif

    for (ind_bil = 0; ind_bil < 2; ind_bil++) {
      for (itot = 0; wk_size[ind_bil] > 1024. && unite[itot] != 'p'; itot++)
        wk_size[ind_bil] /= 1024.;
      wk_unit[ind_bil] = itot;
    }

    bft_printf(_("\n"
                 "  Fortran work arrays memory use:\n"
                 "   %-12llu integers needed (maximum reached in %s)\n"
                 "   %-12llu reals    needed (maximum reached in %s)\n\n"
                 "   Local maximum work memory requested %12.3f %cb\n"
                 "                                  used %12.3f %cb\n"),
               (unsigned long long)_cs_glob_mem_ia_peak, _cs_glob_srt_ia_peak,
               (unsigned long long)_cs_glob_mem_ra_peak, _cs_glob_srt_ra_peak,
               wk_size[0], unite[wk_unit[0]],
               wk_size[1], unite[wk_unit[1]]);

  }

  /* Finalize memory handling */

  if (cs_glob_base_bft_mem_init == true)
    bft_mem_end();

  /* Finalize memory usage count */

  bft_mem_usage_end();
}

/*----------------------------------------------------------------------------
 * Print summary of running time, including CPU and elapsed times.
 *----------------------------------------------------------------------------*/

void
cs_base_time_summary(void)
{
  double  utime;
  double  stime;
  double  time_cpu;
  double  time_tot;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  bft_printf(_("\nCalculation time summary:\n"));

  bft_timer_cpu_times(&utime, &stime);

  if (utime > 0. || stime > 0.)
    time_cpu = utime + stime;

  else
    time_cpu = bft_timer_cpu_time();


  /* CPU time */

  if (utime > 0. || stime > 0.) {
    bft_printf (_("\n  User CPU time:       %12.3f s\n"),
                (float)utime);
    bft_printf (_("  System CPU time:     %12.3f s\n"),
                (float)stime);
  }

  else if (time_cpu > 0.)
    bft_printf (_("\n  CPU time:            %12.3f s\n"),
                (float)time_cpu);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    double time_cumul;
    MPI_Reduce (&time_cpu, &time_cumul, 1, MPI_DOUBLE, MPI_SUM,
                0, cs_glob_mpi_comm);
    if (cs_glob_rank_id == 0)
      bft_printf (_("  Total CPU time:      %12.3f s\n"),
                  time_cumul);
  }
#endif

  /* Elapsed (wall-clock) time */

  time_tot = bft_timer_wtime();

  if (time_tot > 0.) {

    bft_printf (_("\n  Elapsed time:        %12.3f s\n"),
                time_tot);

    bft_printf (_("  CPU / elapsed time   %12.3f\n"),
                (float)(time_cpu/time_tot));

  }

}

/*----------------------------------------------------------------------------
 * Print available system information.
 *----------------------------------------------------------------------------*/

void
cs_base_system_info(void)
{
  time_t          date;
  size_t          ram;

#if defined(HAVE_UNAME)
  struct utsname  sys_config;
#endif

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)
  struct passwd   *pwd_user = NULL;
#endif

#if !defined(PATH_MAX)
#define PATH_MAX 1024
#endif

  char  str_date[81];
  char  str_directory[PATH_MAX] = "";
  char  *str_user = NULL;

  /* Date */

  if (   time(&date) == -1
      || strftime(str_date, 80, "%c", localtime(&date)) == 0)
    strcpy(str_date, "");

  /* Available memory */

  ram = bft_sys_info_mem_ram();

  /* User */

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)

  /* Functions not available on IBM Blue Gene or Cray XT,
     but a stub may exist, so we make sure we ignore it */
#if   defined(__blrts__) || defined(__bgp__) \
   || defined(__CRAYXT_COMPUTE_LINUX_TARGET)
  pwd_user = NULL;
#else
  pwd_user = getpwuid(geteuid());
#endif

  if (pwd_user != NULL) {

    size_t l_user = strlen(pwd_user->pw_name);
    size_t l_info = 0;

    if (pwd_user->pw_gecos != NULL) {
      for (l_info = 0;
           (   pwd_user->pw_gecos[l_info] != '\0'
            && pwd_user->pw_gecos[l_info] != ',');
           l_info++);
    }

    BFT_MALLOC(str_user, l_info + l_user + 3 + 1, char);
    strcpy(str_user, pwd_user->pw_name);

    if (pwd_user->pw_gecos != NULL) {
      strcat(str_user, " (");
      strncpy(str_user + l_user + 2, pwd_user->pw_gecos, l_info);
      str_user[l_user + 2 + l_info]     = ')';
      str_user[l_user + 2 + l_info + 1] = '\0';
    }

  }

#endif /* defined(HAVE_GETPWUID) && defined(HAVE_GETEUID) */

  /* Working directory */

#if defined(HAVE_GETCWD)
  if (getcwd(str_directory, 1024) == NULL)
    strcpy(str_directory, "");
#endif

  /* Print local configuration */
  /*---------------------------*/

  bft_printf("\n%s\n", _("Local case configuration:\n"));

  bft_printf("  %s%s\n", _("Date:              "), str_date);


  /* System and machine */

#if defined(HAVE_UNAME)

  if (uname(&sys_config) != -1) {
    bft_printf("  %s%s %s\n", _("System:            "),
               sys_config.sysname, sys_config.release);
    bft_printf("  %s%s\n", _("Machine:           "),  sys_config.nodename);
  }
#endif

  bft_printf("  %s%s\n", _("Processor:         "), bft_sys_info_cpu());

  if (ram > 0)
    bft_printf("  %s%llu\n", _("Memory:            "), (unsigned long long)ram);

  if (str_user != NULL) {
    bft_printf("  %s%s\n", _("User:              "), str_user);
    BFT_FREE(str_user);
  }

  bft_printf("  %s%s\n", _("Directory:         "), str_directory);

#if defined(HAVE_MPI)
#  if defined(MPI_VERSION) && (MPI_VERSION >= 2)
  {
    int flag = 0;
    void *attp = NULL;
    MPI_Initialized(&flag);
    if (flag != 0) {
      flag = 0;
      MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_APPNUM, &attp, &flag);
    }
    if (flag != 0)
      bft_printf("  %s%d\n", _("MPI appnum:        "), *(int *)attp);
  }
#  endif
  bft_printf("  %s%d\n", _("MPI ranks:         "), cs_glob_n_ranks);
#endif

#if defined(HAVE_OPENMP)
  bft_printf("  %s%d\n", _("OpenMP threads:    "), cs_glob_n_threads);
#endif
}

/*----------------------------------------------------------------------------
 * Replace default bft_printf() mechanism with internal mechanism.
 *
 * This allows redirecting or suppressing logging for different ranks.
 *
 * parameters:
 *   log_name    <-- base file name for log, or NULL for stdout
 *   r0_log_flag <-- redirection for rank 0 log;
 *                   0: not redirected; 1: redirected to <log_name> file
 *   rn_log_flag <-- redirection for ranks > 0 log:
 *                   0: not redirected; 1: redirected to <log_name>_n*" file;
 *                   2: redirected to "/dev/null" (suppressed)
 *----------------------------------------------------------------------------*/

void
cs_base_bft_printf_set(const char  *log_name,
                       int          r0_log_flag,
                       int          rn_log_flag)
{
  /* Non-suppressed logs */

  if (log_name != NULL && (cs_glob_rank_id < 1 || rn_log_flag != 2)) {

    char *filename = NULL;
    BFT_MALLOC(filename, strlen(log_name) + 10, char);

    bft_printf_proxy_set(vprintf);
    bft_printf_flush_proxy_set(_cs_base_bft_printf_flush);
    ple_printf_function_set(vprintf);

    filename[0] = '\0';

    if (cs_glob_rank_id < 1) {
      if (r0_log_flag != 0)
        strcpy(filename, log_name);
    }
    else {
      if (rn_log_flag != 0) {
        if (cs_glob_n_ranks > 9999)
          sprintf(filename, "%s_n%07d", log_name, cs_glob_rank_id+1);
        else
          sprintf(filename, "%s_n%04d", log_name, cs_glob_rank_id+1);
      }
    }

    /* Redirect log */

    if (filename[0] != '\0') {

      FILE *fp = freopen(filename, "w", stdout);

      if (fp == NULL)
        bft_error(__FILE__, __LINE__, errno,
                  _("It is impossible to redirect the standard output "
                    "to file:\n%s"), filename);

#if defined(HAVE_DUP2)
      if (dup2(fileno(fp), fileno(stderr)) == -1)
        bft_error(__FILE__, __LINE__, errno,
                  _("It is impossible to redirect the standard error "
                    "to file:\n%s"), filename);
#endif
    }

    BFT_FREE(filename);
  }

  /* Suppressed logs */

  else if (cs_glob_rank_id > 0) {
    bft_printf_proxy_set(_cs_base_bft_printf_null);
    bft_printf_flush_proxy_set(_cs_base_bft_printf_flush_null);
    ple_printf_function_set(_cs_base_bft_printf_null);
  }
}

/*----------------------------------------------------------------------------
 * Print a warning message header.
 *
 * parameters:
 *   file_name <-- name of source file
 *   line_nume <-- line number in source file
 *----------------------------------------------------------------------------*/

void
cs_base_warn(const char  *file_name,
             int          line_num)
{
  bft_printf(_("\n\nCode_Saturne: %s:%d: Warning\n"),
             file_name, line_num);
}

/*----------------------------------------------------------------------------
 * Convert a character string from the Fortran API to the C API.
 *
 * Eventual leading and trailing blanks are removed.
 *
 * parameters:
 *   f_str <-- Fortran string
 *   f_len <-- Fortran string length
 *
 * returns:
 *   pointer to C string
 *----------------------------------------------------------------------------*/

char *
cs_base_string_f_to_c_create(const char  *f_str,
                             int          f_len)
{
  char * c_str = NULL;
  int    i, i1, i2, l;

  /* Initialization if necessary */

  if (cs_glob_base_str_init == false) {
    for (i = 0 ; i < CS_BASE_N_STRINGS ; i++)
      cs_glob_base_str_is_free[i] = true;
    cs_glob_base_str_init = true;
  }

  /* Handle name for C API */

  for (i1 = 0 ;
       i1 < f_len && (f_str[i1] == ' ' || f_str[i1] == '\t') ;
       i1++);

  for (i2 = f_len - 1 ;
       i2 > i1 && (f_str[i2] == ' ' || f_str[i2] == '\t') ;
       i2--);

  l = i2 - i1 + 1;

  /* Allocation if necessary */

  if (l < CS_BASE_STRING_LEN) {
    for (i = 0 ; i < CS_BASE_N_STRINGS ; i++) {
      if (cs_glob_base_str_is_free[i] == true) {
        c_str = cs_glob_base_str[i];
        cs_glob_base_str_is_free[i] = false;
        break;
      }
    }
  }

  if (c_str == NULL)
    BFT_MALLOC(c_str, l + 1, char);

  for (i = 0 ; i < l ; i++, i1++)
    c_str[i] = f_str[i1];

  c_str[l] = '\0';

  return c_str;
}

/*----------------------------------------------------------------------------
 * Free a string converted from the Fortran API to the C API.
 *
 * parameters:
 *   str <-> pointer to C string
 *----------------------------------------------------------------------------*/

void
cs_base_string_f_to_c_free(char  **c_str)
{
  cs_int_t ind;

  for (ind = 0 ; ind < CS_BASE_N_STRINGS ; ind++) {
    if (*c_str == cs_glob_base_str[ind]) {
      cs_glob_base_str_is_free[ind] = true;
      *c_str = NULL;
      break;
    }
  }

  if (ind == CS_BASE_N_STRINGS && *c_str != NULL)
    BFT_FREE(*c_str);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
