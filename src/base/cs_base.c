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

#if defined(__bgp__)
#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#endif

#if defined(__bgq__)
#include <spi/include/kernel/location.h>
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

#include <fvm_config.h>
#if !defined(HAVE_MPI) && defined(FVM_HAVE_MPI)
#error "Either both or neither Code_Saturne and FVM must be configured with MPI"
#endif

#include <fvm_parall.h>

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

static bft_error_handler_t  *cs_glob_base_gest_erreur_save = NULL;

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

#if defined(__bgq__)
static _cs_base_sighandler_t cs_glob_base_sigtrap_save = SIG_DFL;
#endif

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
 * Print a message to standard output
 *----------------------------------------------------------------------------*/

static int
_cs_base_bft_printf(const char     *const format,
                    va_list         arg_ptr)
{
 cs_int_t  line;
 cs_int_t  msgsize;

 /* Buffer for printing from C code: print to a character string, which will
    be printed to file by Fortran code.
    If Fortran output is completely replaced by C output in the future,
    we will be able to remove this step, but it is currently necessary
    so as to ensure that the same output files may be used and output
    remains ordered. */

#undef CS_BUF_PRINT_F_SIZE
#define CS_BUF_PRINT_F_SIZE 16384

 static char cs_buf_print_f[CS_BUF_PRINT_F_SIZE];

 /* Write to buffer */

#if (_CS_STDC_VERSION < 199901L)
  msgsize = vsprintf (cs_buf_print_f, format, arg_ptr);
#else
  msgsize = vsnprintf (cs_buf_print_f, CS_BUF_PRINT_F_SIZE, format, arg_ptr);
#endif

  line = __LINE__ - 1;

  if (msgsize == -1 || msgsize > CS_BUF_PRINT_F_SIZE - 1) {
    fprintf(stderr,
            _("Fatal error: bft_printf() called on a message of size %d\n"
              "whereas the print buffer is of size %d."),
            msgsize, CS_BUF_PRINT_F_SIZE);

    /* Try to force segmentation fault (to call signal handlers);
       as stack has most likely been corrupted, this is the most
       "similar" error that allows for portable handling. */
    {
      int *_force_err = NULL;
      *_force_err = 0;
    }
    cs_exit(EXIT_FAILURE);
  }

  /* Effective output by Fortran code */

  CS_PROCF (csprnt, CSPRNT) (cs_buf_print_f, &msgsize);

  return msgsize;
}

/*----------------------------------------------------------------------------
 * Flush log output buffer
 *----------------------------------------------------------------------------*/

static int
_cs_base_bft_printf_flush(void)
{
  CS_PROCF (csflsh, CSFLSH) ();

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

#if defined(va_copy)
    va_copy(arg_ptr_2, arg_ptr);
#else
    __va_copy(arg_ptr_2, arg_ptr);
#endif
    _cs_base_bft_printf(format, arg_ptr_2);
    va_end(arg_ptr_2);
  }
#endif

  /* message on a specific error output, initialized only if the
     error output is really necessary */

  if (initialized == false) {

    char nom_fic_err[81];

    if (cs_glob_rank_id < 1)
      strcpy(nom_fic_err, "error");

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
        sprintf(nom_fic_err, "error_n%07d", cs_glob_rank_id + 1);
      else
        sprintf(nom_fic_err, "error_n%04d", cs_glob_rank_id + 1);
    }

    freopen(nom_fic_err, "w", stderr);

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
_cs_base_gestion_erreur(const char  *nom_fic,
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
#if defined(FVM_HAVE_MPI)
  fvm_parall_set_mpi_comm(MPI_COMM_NULL);
#endif

  bft_error_handler_set(cs_glob_base_gest_erreur_save);

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

/*----------------------------------------------------------------------------
 * Copy a Fortan string buffer to a C string buffer
 *
 * The aim of this function is to aviod issues with Fortran array bounds
 * checking when compilers such as icc 11 consider a character array from C
 * as an array of 1-character length strings.
 *
 * Fortran interface
 *
 * SUBROUTINE CSSF2C (LEN, CSTR, FSTR)
 * *****************
 *
 * INTEGER          LEN         : --> : String length
 * CHARACTER*       FSTR        : --> : Fortran string
 * CHARACTER*       CSTR        : <-- : C string
 *----------------------------------------------------------------------------*/

void CS_PROCF (cssf2c, CSSF2C)
(
 const cs_int_t   *len,
 const char       *fstr,
 char             *cstr
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
)
{
  memcpy(cstr, fstr, *len);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(HAVE_MPI)

/*----------------------------------------------------------------------------
 * Complete MPI setup.
 *
 * MPI should have been initialized by cs_opts_mpi_init().
 *
 * Global variables `cs_glob_n_ranks' (number of Code_Saturne processes)
 * and `cs_glob_rank_id' (rank of local process) are set by this function.
 *
 * parameters:
 *   app_num <-- -1 if MPI is not needed, or application number in
 *               MPI_COMM_WORLD of this instance of Code_Saturne.
 *----------------------------------------------------------------------------*/

void
cs_base_mpi_init(int  app_num)
{
  int nbr, rank;

  int app_num_l = app_num, app_num_max = -1;

#if defined(DEBUG) || !defined(NDEBUG)
  MPI_Errhandler errhandler;
#endif

  /*
    If necessary, split MPI_COMM_WORLD to separate different coupled
    applications (collective operation, like all MPI communicator
    creation operations).

    We suppose the color argument to MPI_Comm_split is equal to the
    application number, given through the command line or through
    mpiexec and passed here as an argument.
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Allreduce(&app_num_l, &app_num_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (app_num_max > 0)
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

  if (cs_glob_n_ranks == 1 && app_num_max > 0)
    cs_glob_rank_id = -1;

  /* Initialize associated libraries */

#if defined(FVM_HAVE_MPI)
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

  CS_PROCF(csclli, CSCLLI)(); /* Close log files */

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

  cs_glob_base_gest_erreur_save = bft_error_handler_get();
  bft_error_handler_set(_cs_base_gestion_erreur);

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

#if defined(__bgq__)
  cs_glob_base_sigtrap_save  = signal(SIGTRAP, _cs_base_sig_fatal);
#endif

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

  /* Memory usage measure initialization */

  bft_mem_usage_init();

  /* Memory management initialization */

  if ((nom_base = getenv("CS_FIC_MEM")) != NULL) {

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
cs_base_mem_fin(void)
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
  char  unite[]    = {'K', 'M', 'G', 'T', 'P'};

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
    MPI_Reduce(val_in, val_min, 2, MPI_DOUBLE_INT, MPI_MINLOC,
               0, cs_glob_mpi_comm);
    MPI_Reduce(val_in, val_max, 2, MPI_DOUBLE_INT, MPI_MAXLOC,
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
           valreal[ind_bil] > 1024. && unite[itot] != 'P';
           itot++)
        valreal[ind_bil] /= 1024.;
#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1 && cs_glob_rank_id == 0) {
        for (imin = 0;
             val_min[ind_bil].val > 1024. && unite[imin] != 'P';
             imin++)
          val_min[ind_bil].val /= 1024.;
        for (imax = 0;
             val_max[ind_bil].val > 1024. && unite[imax] != 'P';
             imax++)
          val_max[ind_bil].val /= 1024.;
      }
#endif

      /* Print to log file */

      bft_printf(_("  %s %12.3f %cB\n"),
                 _(type_bil[ind_bil]), valreal[ind_bil], unite[itot]);

#if defined(HAVE_MPI)
      if (cs_glob_n_ranks > 1 && cs_glob_rank_id == 0) {
        bft_printf(_("                             "
                     "local minimum: %12.3f %cB  (rank %d)\n"),
                   val_min[ind_bil].val, unite[imin], val_min[ind_bil].rank);
        bft_printf(_("                             "
                     "local maximum: %12.3f %cB  (rank %d)\n"),
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
                 "   Local maximum work memory requested %12.3f %cB\n"
                 "                                  used %12.3f %cB\n"),
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
cs_base_bilan_temps(void)
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

#if defined(__bgp__)
  _BGP_Personality_t personality;
  Kernel_GetPersonality(&personality, sizeof(personality));
#elif defined(__bgq__)
  Personality_t personality;
  Kernel_GetPersonality(&personality, sizeof(personality));
#endif

  /* Date */

  if (   time(&date) == -1
      || strftime(str_date, 80, "%c", localtime(&date)) == 0)
    strcpy(str_date, "");

  /* Available memory */

#if defined(__bgp__)
  ram = personality.DDR_Config.DDRSizeMB;
#elif defined(__bgq__)
  ram = personality.DDR_Config.DDRSizeMB;
#else
  ram = bft_sys_info_mem_ram() / (size_t)1024;
#endif

  /* User */

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)

  /* Functions not available on IBM Blue Gene or Cray XT,
     but a stub may exist, so we make sure we ignore it */
#if   defined(__blrts__) || defined(__bg__) \
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
    bft_printf("  %s%llu %s\n", _("Memory:            "),
               (unsigned long long)ram, _("MB"));

  if (str_user != NULL) {
    bft_printf("  %s%s\n", _("User:              "), str_user);
    BFT_FREE(str_user);
  }

  bft_printf("  %s%s\n", _("Directory:         "), str_directory);

  /* MPI Info */

#if defined(__bgp__)
  {
    int node_config = personality.Kernel_Config.ProcessConfig;
    const char *_mode_names[] = {N_("SMP mode"),
                                 N_("Virtual-node mode"),
                                 N_("Dual mode"),
                                 N_("Unknown mode")};
    const char *mode_name = _mode_names[3];

    if (node_config == _BGP_PERS_PROCESSCONFIG_SMP)
      mode_name = _mode_names[0];
    else if (node_config == _BGP_PERS_PROCESSCONFIG_VNM)
      mode_name = _mode_names[1];
    else if (node_config == _BGP_PERS_PROCESSCONFIG_2x2)
      mode_name = _mode_names[2];

    bft_printf("  %s%d (%s)\n", _("MPI ranks:         "),
               cs_glob_n_ranks, _(mode_name));
    bft_printf("  %s%d\n", _("pset size:         "),
               personality.Network_Config.PSetSize);
    bft_printf("  %s<%d,%d,%d>\n", _("Torus dimensions:  "),
               personality.Network_Config.Xnodes,
               personality.Network_Config.Ynodes,
               personality.Network_Config.Znodes);
  }
#elif defined(__bgq__)
  {
    int a_torus, b_torus, c_torus, d_torus, e_torus;
    int n_flags = personality.Network_Config.NetFlags;
    int n_world_ranks = cs_glob_n_ranks;

    MPI_Comm_size(MPI_COMM_WORLD, &n_world_ranks);

    bft_printf("  %s%d\n", _("MPI ranks:           "), cs_glob_n_ranks);
    if (n_world_ranks > cs_glob_n_ranks)
      bft_printf("  %s%d\n", _("MPI_COMM_WORLD size: "),
                 n_world_ranks);
    if (n_flags & ND_ENABLE_TORUS_DIM_A) a_torus = 1; else a_torus = 0;
    if (n_flags & ND_ENABLE_TORUS_DIM_B) b_torus = 1; else b_torus = 0;
    if (n_flags & ND_ENABLE_TORUS_DIM_C) c_torus = 1; else c_torus = 0;
    if (n_flags & ND_ENABLE_TORUS_DIM_D) d_torus = 1; else d_torus = 0;
    if (n_flags & ND_ENABLE_TORUS_DIM_E) e_torus = 1; else e_torus = 0;

    bft_printf("  %s<%d,%d,%d,%d,%d>\n", _("Block shape:         "),
               personality.Network_Config.Anodes,
               personality.Network_Config.Bnodes,
               personality.Network_Config.Cnodes,
               personality.Network_Config.Dnodes,
               personality.Network_Config.Enodes);

    bft_printf("  %s<%d,%d,%d,%d,%d>\n", _("Torus links enabled: "),
               a_torus, b_torus, c_torus, d_torus, e_torus);
    }
#elif defined(HAVE_MPI)
  bft_printf("  %s%d\n", _("MPI ranks:         "), cs_glob_n_ranks);
#endif
}

/*----------------------------------------------------------------------------
 * Replace default bft_printf() mechanism with internal mechanism.
 *
 * This is necessary for good consistency of messages output from C or
 * from Fortran, and to handle parallel and serial logging options.
 *----------------------------------------------------------------------------*/

void
cs_base_bft_printf_set(void)
{
  bft_printf_proxy_set(_cs_base_bft_printf);
  bft_printf_flush_proxy_set(_cs_base_bft_printf_flush);
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
