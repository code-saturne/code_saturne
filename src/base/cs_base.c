/*============================================================================
 * Low-level functions and global variables definition.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(HAVE_DLOPEN)
#include <dlfcn.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>
#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_backtrace.h"
#include "bft_mem_usage.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "cs_file.h"
#include "cs_fp_exception.h"
#include "cs_log.h"
#include "cs_timer.h"
#include "cs_version.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Macro Definitions
 *============================================================================*/

#define DIR_SEPARATOR '/'

/* Fortran API */
/*-------------*/

/*
 * 'usual' maximum name length; a longer name is possible, but will
 * provoque a dynamic memory allocation.
 */

#define CS_BASE_N_STRINGS                               5

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

static bft_error_handler_t  *cs_glob_base_err_handler_save = NULL;

static bool  cs_glob_base_bft_mem_init = false;

static bool  cs_glob_base_str_init = false;
static bool  cs_glob_base_str_is_free[CS_BASE_N_STRINGS];
static char  cs_glob_base_str[CS_BASE_N_STRINGS][CS_BASE_STRING_LEN + 1];


/* Global variables associated with signal handling */

#if defined(SIGHUP)
static _cs_base_sighandler_t cs_glob_base_sighup_save = SIG_DFL;
#endif

static _cs_base_sighandler_t cs_glob_base_sigint_save = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigterm_save = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigfpe_save = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigsegv_save = SIG_DFL;

#if defined(__bgq__)
static _cs_base_sighandler_t cs_glob_base_sigabrt_save = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigtrap_save = SIG_DFL;
#endif

#if defined(SIGXCPU)
static _cs_base_sighandler_t cs_glob_base_sigcpu_save = SIG_DFL;
#endif

#if defined(HAVE_DLOPEN)
#if defined(CS_DLOPEN_USE_RTLD_GLOBAL)
static int _cs_dlopen_flags = RTLD_LAZY | RTLD_GLOBAL;
#else
static int _cs_dlopen_flags = RTLD_LAZY;
#endif
#endif

/* Installation paths */

static const char _cs_base_build_localedir[] = LOCALEDIR;
static const char _cs_base_build_pkgdatadir[] = PKGDATADIR;
static const char _cs_base_build_pkglibdir[] = PKGLIBDIR;
static char *_cs_base_env_localedir = NULL;
static char *_cs_base_env_pkgdatadir = NULL;
static char *_cs_base_env_pkglibdir = NULL;

/* Log file */

static FILE  *_bft_printf_file = NULL;
static char  *_bft_printf_file_name = NULL;
static bool   _bft_printf_suppress = false;
static bool   _cs_trace = false;

/* Additional cleanup steps */

static cs_base_atexit_t  * _cs_base_atexit = NULL;

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
  CS_UNUSED(format);
  CS_UNUSED(arg_ptr);

  return 0;
}

/*----------------------------------------------------------------------------
 * False print of a message to standard output for discarded logs
 *----------------------------------------------------------------------------*/

static int
_cs_base_bft_printf_file(const char  *format,
                         va_list      arg_ptr)
{
  return  vfprintf(_bft_printf_file, format, arg_ptr);
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
 * False flush of log output buffer for discarded logs
 *----------------------------------------------------------------------------*/

static int
_cs_base_bft_printf_flush_file(void)
{
  return fflush(_bft_printf_file);
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
  static bool  initialized = false;

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
    int i;
    int n_dec = 1;

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
      double stime = cs_timer_wtime();
      double etime = 0.0;
      do {
        sleep(wait_time);
        etime = cs_timer_wtime();
      }
      while (etime > -0.5 && etime - stime < wait_time); /* etime = -1 only if
                                                            cs_timer_wtime()
                                                            is unusable. */
#endif
      for (i = cs_glob_n_ranks; i >= 10; i /= 10, n_dec += 1);
      sprintf(err_file_name, "error_r%0*d", n_dec, cs_glob_rank_id);
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
  /* Initialize arguments list */

  va_list  arg_ptr;
  va_start(arg_ptr, format);

  /* message sur les sorties */

  _cs_base_err_vprintf(format, arg_ptr);

  /* Finalize arguments list */

  va_end(arg_ptr);
}

/*----------------------------------------------------------------------------
 * Exit function
 *----------------------------------------------------------------------------*/

static void
_cs_base_exit(int status)
{
  if (status == EXIT_SUCCESS)
    cs_base_update_status(NULL);

#if defined(HAVE_MPI)
  {
    int mpi_flag;

    MPI_Initialized(&mpi_flag);

#if (MPI_VERSION >= 2)
    if (mpi_flag != 0) {
      int finalized_flag;
      MPI_Finalized(&finalized_flag);
      if (finalized_flag != 0)
        mpi_flag = 0;
    }
#endif

    if (mpi_flag != 0) {

      /* For safety, flush all streams before calling MPI_Abort
       * (should be done by exit, but in case we call MPI_Abort
       * due to a SIGTERM received from another rank's MPI_Abort,
       * make sure we avoid ill-defined behavior) */

      fflush(NULL);

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
  if (_cs_base_atexit != NULL) {
    _cs_base_atexit();
    _cs_base_atexit = NULL;
  }

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

  mem_usage = bft_mem_size_max();

  if (mem_usage > 0)
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
  if (_cs_base_atexit != NULL) {
    _cs_base_atexit();
    _cs_base_atexit = NULL;
  }

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

#if defined(__bgq__)
  case SIGABRT:
    _cs_base_err_printf(_("SIGABRT signal (Abort) "
                          "intercepted.\n"));
    cs_glob_base_sigabrt_save  = signal(SIGABRT, cs_glob_base_sigabrt_save);
    break;
#endif

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

/*----------------------------------------------------------------------------
 * Return a string providing path information.
 *
 * This is normally the path determined upon configuration, but may be
 * adapted for movable installs using the CS_ROOT_DIR environment variable
 * or by a guess on the assumed relative path.
 *----------------------------------------------------------------------------*/

static const char *
_get_path(const char   *dir_path,
          const char   *build_path,
          char        **env_path)
{
#if defined(HAVE_RELOCATABLE)
  {
    const char *cs_root_dir = NULL;
    const char *rel_path = NULL;

    /* Allow for displacable install */

    if (*env_path != NULL)
      return *env_path;

    /* First try with an environment variable CS_ROOT_DIR */

    if (getenv("CS_ROOT_DIR") != NULL) {
      cs_root_dir = getenv("CS_ROOT_DIR");
      rel_path = "/";
    }

    /* Second try with an environment variable CFDSTUDY_ROOT_DIR */

    else if (getenv("CFDSTUDY_ROOT_DIR") != NULL) {
      cs_root_dir = getenv("CFDSTUDY_ROOT_DIR");
      rel_path = "/";
    }

#if defined(HAVE_GETCWD)

    /*
      Then, try to guess a relative path, knowing that executables are
      located in libexecdir/code_saturne
    */

    else {

      int buf_size = 128;
      char *buf = NULL;

      while (cs_root_dir == NULL) {
        buf_size *= 2;
        BFT_REALLOC(buf, buf_size, char);
        cs_root_dir = getcwd(buf, buf_size);
        if (cs_root_dir == NULL && errno != ERANGE)
          bft_error(__FILE__, __LINE__, errno,
                    _("Error querying working directory.\n"));
      }

      rel_path = "/../../";

    }
#endif /* defined(HAVE_GETCWD) */

    BFT_MALLOC(*env_path,
               strlen(cs_root_dir) + strlen(rel_path) + strlen(dir_path) + 1,
               char);
    strcpy(*env_path, cs_root_dir);
    strcat(*env_path, rel_path);
    strcat(*env_path, dir_path);

    return *env_path;
  }
#else

  CS_UNUSED(dir_path);
  CS_UNUSED(env_path);

#endif /* defined(HAVE_RELOCATABLE) */

  /* Standard install */

  return build_path;
}

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 *  Finalisation MPI
 *----------------------------------------------------------------------------*/

static void
_cs_base_mpi_fin(void)
{
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
_cs_base_mpi_error(MPI_Comm  *comm,
                   int       *errcode,
                   ...)
{
  int err_len;
  char err_string[MPI_MAX_ERROR_STRING + 1];

#if defined MPI_MAX_OBJECT_NAME
  int name_len = 0;
  char comm_name[MPI_MAX_OBJECT_NAME + 1];
#endif

  if (_cs_base_atexit != NULL) {
    _cs_base_atexit();
    _cs_base_atexit = NULL;
  }

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
 * Ensure Code_Saturne to MPI datatype conversion has correct values.
 *----------------------------------------------------------------------------*/

static void
_cs_datatype_to_mpi_init(void)
{
  int size_short, size_int, size_long, size_long_long;

  MPI_Type_size(MPI_SHORT, &size_short);
  MPI_Type_size(MPI_INT,   &size_int);
  MPI_Type_size(MPI_LONG,  &size_long);

#if defined(MPI_LONG_LONG)
  MPI_Type_size(MPI_LONG_LONG, &size_long_long);
#else
  size_long_long = 0;
#endif

  if (size_int == 4) {
    cs_datatype_to_mpi[CS_INT32] = MPI_INT;
    cs_datatype_to_mpi[CS_UINT32] = MPI_UNSIGNED;
  }
  else if (size_short == 4) {
    cs_datatype_to_mpi[CS_INT32] = MPI_SHORT;
    cs_datatype_to_mpi[CS_UINT32] = MPI_UNSIGNED_SHORT;
  }
  else if (size_long == 4) {
    cs_datatype_to_mpi[CS_INT32] = MPI_LONG;
    cs_datatype_to_mpi[CS_UINT32] = MPI_UNSIGNED_LONG;
  }

  if (size_int == 8) {
    cs_datatype_to_mpi[CS_INT64] = MPI_INT;
    cs_datatype_to_mpi[CS_UINT64] = MPI_UNSIGNED;
  }
  else if (size_long == 8) {
    cs_datatype_to_mpi[CS_INT64] = MPI_LONG;
    cs_datatype_to_mpi[CS_UINT64] = MPI_UNSIGNED_LONG;
  }
#if defined(MPI_LONG_LONG)
  else if (size_long_long == 8) {
    cs_datatype_to_mpi[CS_INT64] = MPI_LONG_LONG;
#if defined(MPI_UNSIGNED_LONG_LONG)
    cs_datatype_to_mpi[CS_UINT64] = MPI_UNSIGNED_LONG_LONG;
#else
    cs_datatype_to_mpi[CS_UINT64] = MPI_LONG_LONG;
#endif
  }
#endif
}

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
 *----------------------------------------------------------------------------*/

static void
_cs_base_mpi_setup(const char *app_name)
{
  int nbr, rank;

  int app_num = -1;

#if (defined(DEBUG) || !defined(NDEBUG)) && (MPI_VERSION >= 2)
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

  /* Initialize datatype conversion */

  _cs_datatype_to_mpi_init();

  /* Initialize error handlers */

#if (defined(DEBUG) || !defined(NDEBUG)) && (MPI_VERSION >= 2)
  if (nbr > 1 || cs_glob_mpi_comm != MPI_COMM_NULL) {
    MPI_Comm_create_errhandler(&_cs_base_mpi_error, &errhandler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, errhandler);
    if (   cs_glob_mpi_comm != MPI_COMM_WORLD
        && cs_glob_mpi_comm != MPI_COMM_NULL)
      MPI_Comm_set_errhandler(cs_glob_mpi_comm, errhandler);
    MPI_Errhandler_free(&errhandler);
  }
#endif
}

#endif /* HAVE_MPI */

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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

#if defined(MPI_SUBVERSION)

  char mpi_vendor_lib[32] = "";
  char mpi_lib[32] = "";

  /* Base MPI library information */

#if defined(MPI_VENDOR_NAME)

#if defined(OMPI_MAJOR_VERSION)
  snprintf(mpi_lib, 31, "%s %d.%d.%d",
           MPI_VENDOR_NAME,
           OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION, OMPI_RELEASE_VERSION);
#elif defined(MPICH2_VERSION)
  snprintf(mpi_lib, 31, "%s %s", MPI_VENDOR_NAME, MPICH2_VERSION);
#elif defined(MPICH_VERSION)
  snprintf(mpi_lib, 31, "%s %s", MPI_VENDOR_NAME, MPICH_VERSION);
#else
  snprintf(mpi_lib, 31, "%s", MPI_VENDOR_NAME);
#endif

#elif defined(OPEN_MPI)
#if defined(OMPI_MAJOR_VERSION)
  snprintf(mpi_lib, 31, "Open MPI %d.%d.%d",
           OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION, OMPI_RELEASE_VERSION);
#else
  snprintf(mpi_lib, 31, "Open MPI");
#endif

#elif defined(MPICH2)
#if defined(MPICH2_VERSION)
  snprintf(mpi_lib, 31, "MPICH2 %s", MPICH2_VERSION);
#else
  snprintf(mpi_lib, 31, "MPICH2");
#endif
#elif defined(MPICH_NAME)
#if defined(MPICH_VERSION)
  snprintf(mpi_lib, 31, "MPICH %s", MPICH_VERSION);
#else
  snprintf(mpi_lib, 31, "MPICH");
#endif

#elif defined(LAM_MPI)
#if defined(LAM_VERSION)
  snprintf(mpi_lib, 31, "LAM/MPI %s", LAM_VERSION);
#else
  snprintf(mpi_lib, 31, "LAM/MPI");
#endif
#endif

  mpi_lib[31] = '\0';

  /* Possible additional MPI vendor information */

#if defined(MVAPICH2_VERSION)
  snprintf(mpi_vendor_lib, 31, "MVAPICH2 %s", MVAPICH2_VERSION);
#elif defined(MSMPI_VER)
  snprintf(mpi_vendor_lib, 31, "MS-MPI");
#elif defined(PLATFORM_MPI)
  {
    int v, v0, v1, v2, v3;
    v0 =  PLATFORM_MPI>>24;
    v =  (PLATFORM_MPI - (v0<<24));
    v1 =  v>>16;
    v =  (v - (v1<<16));
    v2 =  v>>8;
    v3 =  (v - (v2<<8));
    snprintf(mpi_vendor_lib, 31, "Platform MPI %x.%x.%x.%x\n",
             v0, v1, v2, v3);
  }
#endif

  mpi_vendor_lib[31] = '\0';

#endif /* defined(MPI_SUBVERSION) */

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

  bft_printf("\n  Copyright (C) 1998-2019 EDF S.A., France\n\n");

#if defined(CS_REVISION)
  if (strlen(CS_REVISION) > 0)
    bft_printf(_("  revision %s\n"), CS_REVISION);
#endif

  bft_printf(_("  build %s\n"), str);

#if defined(MPI_SUBVERSION)
  if (mpi_vendor_lib[0] != '\0') {
    if (mpi_lib[0] != '\0')
      bft_printf(_("  MPI version %d.%d (%s, based on %s)\n\n"),
                 MPI_VERSION, MPI_SUBVERSION, mpi_vendor_lib, mpi_lib);
    else
      bft_printf(_("  MPI version %d.%d (%s)\n\n"),
                 MPI_VERSION, MPI_SUBVERSION, mpi_vendor_lib);
  }
  else {
    if (mpi_lib[0] != '\0')
      bft_printf(_("  MPI version %d.%d (%s)\n\n"),
                 MPI_VERSION, MPI_SUBVERSION, mpi_lib);
    else
      bft_printf(_("  MPI version %d.%d\n\n"),
                 MPI_VERSION, MPI_SUBVERSION);
  }
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
 *----------------------------------------------------------------------------*/

void
cs_base_mpi_init(int    *argc,
                 char  **argv[])
{
#if defined(HAVE_MPI)

  char *s;

  int arg_id = 0, flag = 0;
  int use_mpi = false;

#if   defined(__bg__) || defined(__CRAYXT_COMPUTE_LINUX_TARGET)

  /* Blue Gene/Q or Cray: assume MPI is always used. */

  use_mpi = true;

#elif defined(MPICH) || defined(MPICH2) || defined(MSMPI_VER)

  /* Notes: Microsoft MPI is based on MPICH */

  if (getenv("PMI_RANK") != NULL)
    use_mpi = true;

  else if (getenv("PCMPI") != NULL) /* Platform MPI */
    use_mpi = true;

#elif defined(OPEN_MPI)
  if (getenv("OMPI_MCA_ns_nds_vpid") != NULL)         /* OpenMPI 1.2 */
    use_mpi = true;
  else if (getenv("OMPI_COMM_WORLD_RANK") != NULL)    /* OpenMPI 1.3 + */
    use_mpi = true;

#endif /* Tests for known MPI variants */

  /* Test for run through SLURM's srun */

  if (getenv("SLURM_SRUN_COMM_HOST") != NULL)
    use_mpi = true;

  /* If we have determined from known MPI environment variables
     of command line arguments that we are running under MPI,
     initialize MPI */

  if (use_mpi == true) {
    MPI_Initialized(&flag);
    if (!flag) {
#if (MPI_VERSION >= 2) && defined(HAVE_OPENMP)
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
#if (MPI_VERSION >= 2) && defined(HAVE_OPENMP)
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

    _cs_base_mpi_setup(app_name);

    BFT_FREE(app_name);
  }

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
  if (_cs_base_atexit != NULL) {
    _cs_base_atexit();
    _cs_base_atexit = NULL;
  }

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
 *
 * parameters:
 *   signal_defaults <-- leave default signal handlers in place if true
 *----------------------------------------------------------------------------*/

void
cs_base_error_init(bool  signal_defaults)
{
  /* Error handler */

  cs_glob_base_err_handler_save = bft_error_handler_get();
  bft_error_handler_set(_cs_base_error_handler);
  ple_error_handler_set(_cs_base_error_handler);

  /* Signal handlers */

  if (signal_defaults == false) {

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

#if defined(__bgq__)
    cs_glob_base_sigabrt_save  = signal(SIGABRT, _cs_base_sig_fatal);
    cs_glob_base_sigtrap_save  = signal(SIGTRAP, _cs_base_sig_fatal);
#endif

#if defined(SIGXCPU)
    if (cs_glob_rank_id <= 0)
      cs_glob_base_sigcpu_save = signal(SIGXCPU, _cs_base_sig_fatal);
#endif

  }
}

/*----------------------------------------------------------------------------
 * Initialize management of memory allocated through BFT.
 *----------------------------------------------------------------------------*/

void
cs_base_mem_init(void)
{
  char  *base_name;
  char  *file_name = NULL;

  /* Set error handler */

  bft_mem_error_handler_set(_cs_mem_error_handler);

  /* Set PLE library memory handler */

  ple_mem_functions_set(bft_mem_malloc,
                        bft_mem_realloc,
                        bft_mem_free);

  /* Memory usage measure initialization */

  bft_mem_usage_init();

  /* Memory management initialization */

  if (bft_mem_initialized())
    cs_glob_base_bft_mem_init = false;

  else {

    if ((base_name = getenv("CS_MEM_LOG")) != NULL) {

      /* We may not use BFT_MALLOC here as memory management has
         not yet been initialized using bft_mem_init() */

      if (base_name != NULL) {

        /* In parallel, we will have one trace file per MPI process */
        if (cs_glob_rank_id >= 0) {
          int i;
          int n_dec = 1;
          for (i = cs_glob_n_ranks; i >= 10; i /= 10, n_dec += 1);
          file_name = malloc((strlen(base_name) + n_dec + 2) * sizeof (char));
          sprintf(file_name, "%s.%0*d", base_name, n_dec, cs_glob_rank_id);
        }
        else {
          file_name = malloc((strlen(base_name) + 1) * sizeof (char));
          strcpy(file_name, base_name);
        }

        /* Actually initialize bft_mem instrumentation only when
           CS_MEM_LOG is defined (for better performance) */

        bft_mem_init(file_name);

        free (file_name);

      }

    }

    cs_glob_base_bft_mem_init = true;

  }
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
  double val_sum[2];
  int  ind_min[2];
  _cs_base_mpi_double_int_t  val_in[2], val_min[2], val_max[2];
#endif

  int   ind_val[2] = {1, 1};
  const char  unit[8] = {'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'};

  const char  * type_bil[] = {N_("Total memory used:                       "),
                              N_("Theoretical instrumented dynamic memory: ")};

  /* Memory summary */

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\nMemory use summary:\n\n"));

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
    MPI_Reduce(valreal, val_sum, 2, MPI_DOUBLE, MPI_SUM,
               0, cs_glob_mpi_comm);
    for (ind_bil = 0; ind_bil < 2; ind_bil++) {
      val_in[ind_bil].val = valreal[ind_bil];
      val_in[ind_bil].rank = cs_glob_rank_id;
    }
    MPI_Reduce(val_in, val_min, 2, MPI_DOUBLE_INT, MPI_MINLOC,
               0, cs_glob_mpi_comm);
    MPI_Reduce(val_in, val_max, 2, MPI_DOUBLE_INT, MPI_MAXLOC,
               0, cs_glob_mpi_comm);
    if (cs_glob_rank_id == 0) {
      for (ind_bil = 0; ind_bil < 2; ind_bil++) {
        ind_val[ind_bil]  = ind_min[ind_bil];
        valreal[ind_bil] = val_sum[ind_bil];
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
           valreal[ind_bil] > 1024. && itot < 8;
           itot++)
        valreal[ind_bil] /= 1024.;
#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1 && cs_glob_rank_id == 0) {
        for (imin = 0;
             val_min[ind_bil].val > 1024. && imin < 8;
             imin++)
          val_min[ind_bil].val /= 1024.;
        for (imax = 0;
             val_max[ind_bil].val > 1024. && imax < 8;
             imax++)
          val_max[ind_bil].val /= 1024.;
      }
#endif

      /* Print to log file */

      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("  %s %12.3f %ciB\n"),
                    _(type_bil[ind_bil]), valreal[ind_bil], unit[itot]);

#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1 && cs_glob_rank_id == 0) {
        cs_log_printf(CS_LOG_PERFORMANCE,
                      _("                             "
                        "local minimum: %12.3f %ciB  (rank %d)\n"),
                      val_min[ind_bil].val, unit[imin], val_min[ind_bil].rank);
        cs_log_printf(CS_LOG_PERFORMANCE,
                      _("                             "
                        "local maximum: %12.3f %ciB  (rank %d)\n"),
                      val_max[ind_bil].val, unit[imax], val_max[ind_bil].rank);
      }
#endif
    }

  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);

  /* Finalize memory handling */

  if (cs_glob_base_bft_mem_init == true) {

    BFT_FREE(_cs_base_env_localedir);
    BFT_FREE(_cs_base_env_pkgdatadir);
    BFT_FREE(_cs_base_env_pkglibdir);
    BFT_FREE(_bft_printf_file_name);

    bft_mem_end();

  }

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

  cs_log_printf(CS_LOG_PERFORMANCE,
                _("\nCalculation time summary:\n"));

  cs_timer_cpu_times(&utime, &stime);

  if (utime > 0. || stime > 0.)
    time_cpu = utime + stime;

  else
    time_cpu = cs_timer_cpu_time();

  /* CPU time */

  if (utime > 0. || stime > 0.) {
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n  User CPU time:       %12.3f s\n"),
                  (float)utime);
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  System CPU time:     %12.3f s\n"),
                  (float)stime);
  }

  else if (time_cpu > 0.)
    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n  CPU time:            %12.3f s\n"),
                  (float)time_cpu);

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {
    double time_cumul;
    MPI_Reduce (&time_cpu, &time_cumul, 1, MPI_DOUBLE, MPI_SUM,
                0, cs_glob_mpi_comm);
    if (cs_glob_rank_id == 0)
      cs_log_printf(CS_LOG_PERFORMANCE,
                    _("  Total CPU time:      %12.3f s\n"),
                    time_cumul);
  }
#endif

  /* Elapsed (wall-clock) time */

  time_tot = cs_timer_wtime();

  if (time_tot > 0.) {

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("\n  Elapsed time:        %12.3f s\n"),
                  time_tot);

    cs_log_printf(CS_LOG_PERFORMANCE,
                  _("  CPU / elapsed time   %12.3f\n"),
                  (float)(time_cpu/time_tot));

  }

  cs_log_printf(CS_LOG_PERFORMANCE, "\n");
  cs_log_separator(CS_LOG_PERFORMANCE);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update status file.
 *
 * If the format string is NULL, the file is removed.

 * \param[in]  format  format string, or NULL
 * \param[in]  ...     format arguments
 */
/*----------------------------------------------------------------------------*/

void
cs_base_update_status(const char  *format,
                      ...)
{
  static const char _status_file_name[] = "run_status.running";
  static FILE *_status_file = NULL;

  if (cs_glob_rank_id < 1) {

    if (format == NULL) {
      if (_status_file != NULL) {
        if (fclose(_status_file) == 0) {
          _status_file = NULL;
          remove(_status_file_name);
        }
      }
    }

    else {

      va_list  arg_ptr;
      va_start(arg_ptr, format);

      /* Output to trace */

#if defined(va_copy) || defined(__va_copy)
      if (_cs_trace && format != NULL) {
        va_list arg_ptr_2;
#if defined(va_copy)
        va_copy(arg_ptr_2, arg_ptr);
#else
        __va_copy(arg_ptr_2, arg_ptr);
#endif
        vprintf(format, arg_ptr_2);
        va_end(arg_ptr_2);
      }
#endif

      /* Status file */

      if (_status_file == NULL)
        _status_file = fopen(_status_file_name, "w");

      if (_status_file != NULL) {
        long p_size = ftell(_status_file);
        fseek(_status_file, 0, SEEK_SET);
        vfprintf(_status_file, format, arg_ptr);
        long c_size = ftell(_status_file);

        while (p_size > c_size) {
          size_t l = 0;
          char buf[64];
          while (l < 64 && p_size > c_size) {
            buf[l++] = ' ';
            p_size--;
          }
          fwrite(buf, 1, l, _status_file);
        }
      }

      va_end(arg_ptr);

    }

  }
}

/*----------------------------------------------------------------------------
 * Set tracing of progress on or off.
 *
 * This function should be called before cs_base_bft_printf_set() if tracing
 * is activated.
 *
 * parameters:
 *   trace  <-- trace progress to stdout
 *----------------------------------------------------------------------------*/

void
cs_base_trace_set(bool  trace)
{
  if (_bft_printf_file_name == NULL)
    _cs_trace = trace;
}

/*----------------------------------------------------------------------------
 * Set output file name and suppression flag for bft_printf().
 *
 * This allows redirecting or suppressing logging for different ranks.
 *
 * parameters:
 *   log_name    <-- base file name for log
 *   rn_log_flag <-- redirection for ranks > 0 log:
 *                   false: to "/dev/null" (suppressed)
 *                   true: to <log_name>_r*.log" file;
 *----------------------------------------------------------------------------*/

void
cs_base_bft_printf_init(const char  *log_name,
                        bool         rn_log_flag)
{
  BFT_FREE(_bft_printf_file_name);
  _bft_printf_suppress = false;

  const char ext[] = ".log";

  /* Allow bypassing this with environment variable to accomodate
     some debug habits */

  bool log_to_stdout = false;
  const char *p = getenv("CS_LOG_TO_STDOUT");
  if (p != NULL) {
    if (atoi(p) > 0)
      log_to_stdout = true;
  }

  /* Rank 0 */

  if (   cs_glob_rank_id < 1
      && log_name != NULL
      && log_to_stdout == false) {

    BFT_MALLOC(_bft_printf_file_name,
               strlen(log_name) + strlen(ext) + 1,
               char);
    strcpy(_bft_printf_file_name, log_name);
    strcat(_bft_printf_file_name, ext);

  }

  /* Other ranks */

  else if (cs_glob_rank_id > 0) {

    if (log_name != NULL && rn_log_flag > 0) { /* Non-suppressed logs */

      if (log_to_stdout == false) {
        int n_dec = 1;
        for (int i = cs_glob_n_ranks; i >= 10; i /= 10, n_dec += 1);
        BFT_MALLOC(_bft_printf_file_name,
                   strlen(log_name) + n_dec + 3 + strlen(ext), char);
        sprintf(_bft_printf_file_name,
                "%s_r%0*d%s",
                log_name,
                n_dec,
                cs_glob_rank_id,
                ext);
      }

    }

    else { /* Suppressed logs */

      _bft_printf_suppress = true;
      bft_printf_proxy_set(_cs_base_bft_printf_null);
      bft_printf_flush_proxy_set(_cs_base_bft_printf_flush_null);
      ple_printf_function_set(_cs_base_bft_printf_null);

    }

  }
}

/*----------------------------------------------------------------------------
 * Replace default bft_printf() mechanism with internal mechanism.
 *
 * This allows redirecting or suppressing logging for different ranks.
 *
 * parameters:
 *   log_name    <-- base file name for log
 *   rn_log_flag <-- redirection for ranks > 0 log:
 *                   false: to "/dev/null" (suppressed)
 *                   true: to <log_name>_r*.log" file;
 *----------------------------------------------------------------------------*/

void
cs_base_bft_printf_set(const char  *log_name,
                       bool         rn_log_flag)
{
  cs_base_bft_printf_init(log_name, rn_log_flag);

  if (_bft_printf_file_name != NULL && _bft_printf_suppress == false) {

    /* Redirect log */

    if (_bft_printf_file_name != NULL) {

      bft_printf_proxy_set(vprintf);
      bft_printf_flush_proxy_set(_cs_base_bft_printf_flush);
      ple_printf_function_set(vprintf);

      if (cs_glob_rank_id > 0 || _cs_trace == false) {

        FILE *fp = freopen(_bft_printf_file_name, "w", stdout);

        if (fp == NULL)
          bft_error(__FILE__, __LINE__, errno,
                    _("It is impossible to redirect the standard output "
                      "to file:\n%s"), _bft_printf_file_name);

#if defined(HAVE_DUP2)
        if (dup2(fileno(fp), fileno(stderr)) == -1)
          bft_error(__FILE__, __LINE__, errno,
                    _("It is impossible to redirect the standard error "
                      "to file:\n%s"), _bft_printf_file_name);
#endif

      }
      else {

        _bft_printf_file = fopen(_bft_printf_file_name, "w");
        if (_bft_printf_file == NULL)
          bft_error(__FILE__, __LINE__, errno,
                    _("Error opening log file:\n%s"),
                    _bft_printf_file_name);

        bft_printf_proxy_set(_cs_base_bft_printf_file);
        bft_printf_flush_proxy_set(_cs_base_bft_printf_flush_file);
        ple_printf_function_set(_cs_base_bft_printf_file);

      }

    }

  }

}

/*----------------------------------------------------------------------------
 * Return name of default log file.
 *
 * cs_base_bft_printf_set or cs_base_c_bft_printf_set() must have
 * been called before this.
 *
 * returns:
 *   name of default log file
 *----------------------------------------------------------------------------*/

const char *
cs_base_bft_printf_name(void)
{
  return _bft_printf_file_name;
}

/*----------------------------------------------------------------------------
 * Return flag indicating if the default log file output is suppressed.
 *
 * cs_base_bft_printf_set or cs_base_c_bft_printf_set() must have
 * been called before this.
 *
 * returns:
 *   name of default log file
 *----------------------------------------------------------------------------*/

bool
cs_base_bft_printf_suppressed(void)
{
  return _bft_printf_suppress;
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
 * Define a function to be called when entering cs_exit() or bft_error().
 *
 * Compared to the C atexit(), only one function may be called (latest
 * setting wins), but the function is called slighty before exit,
 * so it is well adapted to cleanup such as flushing of non-C API logging.
 *
 * parameters:
 *   fct <-- pointer tu function to be called
 *----------------------------------------------------------------------------*/

void
cs_base_atexit_set(cs_base_atexit_t  *const fct)
{
  _cs_base_atexit = fct;
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

/*----------------------------------------------------------------------------
 * Clean a string representing options.
 *
 * Characters are converted to lowercase, leading and trailing whitespace
 * is removed, and multiple whitespaces or tabs are replaced by single
 * spaces.
 *
 * parameters:
 *   s <-> string to be cleaned
 *----------------------------------------------------------------------------*/

void
cs_base_option_string_clean(char  *s)
{
  if (s != NULL) {

    int i, j;

    int l = strlen(s);

    for (i = 0, j = 0 ; i < l ; i++) {
      s[j] = tolower(s[i]);
      if (s[j] == ',' || s[j] == ';' || s[j] == '\t')
        s[j] = ' ';
      if (s[j] != ' ' || (j > 0 && s[j-1] != ' '))
        j++;
    }
    if (j > 0 && s[j-1] == ' ')
      j--;

    s[j] = '\0';
  }
}

/*----------------------------------------------------------------------------
 * Return a string providing locale path information.
 *
 * returns:
 *   locale path
 *----------------------------------------------------------------------------*/

const char *
cs_base_get_localedir(void)
{
  return _get_path("share/locale",
                   _cs_base_build_localedir,
                   &_cs_base_env_localedir);
}

/*----------------------------------------------------------------------------
 * Return a string providing package data path information.
 *
 * returns:
 *   package data path
 *----------------------------------------------------------------------------*/

const char *
cs_base_get_pkgdatadir(void)
{
  return _get_path("share/" PACKAGE_NAME,
                   _cs_base_build_pkgdatadir,
                   &_cs_base_env_pkgdatadir);
}

/*----------------------------------------------------------------------------
 * Return a string providing loadable library path information.
 *
 * This is normally the path determined upon configuration, but may be
 * adapted for movable installs using the CS_ROOT_DIR environment variable.
 *
 * returns:
 *   package loadable library (plugin) path
 *----------------------------------------------------------------------------*/

const char *
cs_base_get_pkglibdir(void)
{
  return _get_path("lib/" PACKAGE_NAME,
                   _cs_base_build_pkglibdir,
                   &_cs_base_env_pkglibdir);
}

/*----------------------------------------------------------------------------
 * Ensure bool argument has value 0 or 1.
 *
 * This allows working around issues with Intel compiler C bindings,
 * which seem to pass incorrect values in some cases.
 *
 * parameters:
 *   b <-> pointer to bool
 *----------------------------------------------------------------------------*/

void
cs_base_check_bool(bool *b)
{
  if (sizeof(bool) == 1) {
    char *pb = (char *)b;
    int i = *pb;
    if (i != 0 && i != 1)
      *b = true;
  }
  else if (sizeof(bool) == sizeof(int)) {
    int *pb = (int *)b;
    if (*pb != 0 && *pb != 1)
      *b = true;
  }
}

/*----------------------------------------------------------------------------
 * Open a data file in read mode.
 *
 * If a file of the given name in the working directory is found, it
 * will be opened. Otherwise, it will be searched for in the "data/thch"
 * subdirectory of pkgdatadir.
 *
 * parameters:
 *   base_name      <-- base file name
 *
 * returns:
 *   pointer to opened file
 *----------------------------------------------------------------------------*/

FILE *
cs_base_open_properties_data_file(const char  *base_name)
{
  FILE *f = NULL;

  char *_f_name = NULL;
  const char *file_name = base_name;

  /* choose local file if present, default otherwise */

  if (! cs_file_isreg(file_name)) {
    const char *datadir = cs_base_get_pkgdatadir();
    const char subdir[] = "/data/thch/";
    BFT_MALLOC(_f_name,
               strlen(datadir) + strlen(subdir) + strlen(base_name) + 1,
               char);
    sprintf(_f_name, "%s%s%s", datadir, subdir, base_name);
    file_name = _f_name;
  }

  f = fopen(file_name, "r");

  if (f == NULL)
    bft_error(__FILE__, __LINE__, errno,
              _("Error opening data file \"%s\""), file_name);

  BFT_FREE(_f_name);

  return f;
}

#if defined(HAVE_DLOPEN)

/*----------------------------------------------------------------------------*/
/*!
 * \brief Load a dynamic library.
 *
 * \param[in]  filename  path to shared library file
 *
 * \return  handle to shared library
 */
/*----------------------------------------------------------------------------*/

void*
cs_base_dlopen(const char *filename)
{
  void *retval = NULL;

  /* Disable floating-point traps as the initialization of some libraries
     may interfere with this (for example, embree, and optional Paraview
     depedency) */

  cs_fp_exception_disable_trap();

  /* Load symbols from shared library */

  retval = dlopen(filename, _cs_dlopen_flags);

  if (retval == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Error loading %s: %s."), filename, dlerror());

  /* Restore floating-point trap behavior */

  cs_fp_exception_restore_trap();

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Load a plugin's dynamic library
 *
 * This function is similar to \ref cs_base_dlopen, execpt that only
 * the base plugin file name (with no extension) needs to be given.
 * It is assumed the file is available in the code's "pkglibdir" directory,
 *
 * \param[in]  name  path to shared library file
 *
 * \return  handle to shared library
 */
/*----------------------------------------------------------------------------*/

void *
cs_base_dlopen_plugin(const char *name)
{
  void *retval = NULL;

  char  *lib_path = NULL;
  const char *pkglibdir = cs_base_get_pkglibdir();

  /* Open shared library */

  BFT_MALLOC(lib_path,
             strlen(pkglibdir) + 1 + 3 + strlen(name) + 3 + 1,
             char);

  sprintf(lib_path, "%s%c%s.so", pkglibdir, DIR_SEPARATOR, name);

  retval = cs_base_dlopen(lib_path);

  BFT_FREE(lib_path);

  return retval;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get flags for dlopen.
 *
 * \return  flags used for dlopen.
 */
/*----------------------------------------------------------------------------*/

int
cs_base_dlopen_get_flags(void)
{
  return _cs_dlopen_flags;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Set flags for dlopen.
 *
 * \param[in]  flags  flags to set
 */
/*----------------------------------------------------------------------------*/

void
cs_base_dlopen_set_flags(int flags)
{
  _cs_dlopen_flags = flags;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Unload a dynamic library.
 *
 * Note that the dlopen underlying mechanism uses a reference count, so
 * a library is really unloaded only one \ref cs_base_dlclose has been called
 * the same number of times as \ref cs_base_dlopen.
 *
 * \param[in]  filename  optional path to shared library file name for error
 *                       logging, or NULL
 * \param[in]  handle    handle to shared library
 */
/*----------------------------------------------------------------------------*/

void
cs_base_dlclose(const char  *filename,
                void        *handle)
{
  int retval = 0;

  if (handle != NULL)
    retval = dlclose(handle);

  if (retval != 0) {
    if (filename != NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("Error decrementing count or unloading %s: %s."),
                filename, dlerror());
    else
      bft_error(__FILE__, __LINE__, 0,
                _("Error decrementing count or unloading %s."),
                dlerror());
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get a shared library function pointer
 *
 * \param[in]  handle            handle to shared library
 * \param[in]  name              name of function symbol in library
 * \param[in]  errors_are_fatal  abort if true, silently ignore if false
 *
 * \return  pointer to function in shared library
 */
/*----------------------------------------------------------------------------*/

void *
cs_base_get_dl_function_pointer(void        *handle,
                                const char  *name,
                                bool         errors_are_fatal)
{
  void  *retval = NULL;
  char  *error = NULL;

  dlerror();    /* Clear any existing error */

  retval = dlsym(handle, name);
  error = dlerror();

  if (error != NULL && errors_are_fatal)
    bft_error(__FILE__, __LINE__, 0,
              _("Error calling dlsym for %s: %s\n"), name, error);

  return retval;
}

#endif /* defined(HAVE_DLOPEN) */

/*----------------------------------------------------------------------------*/

END_C_DECLS
