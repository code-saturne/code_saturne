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

/*============================================================================
 * Parsing of program arguments and associated initializations
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_gui_util.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_opts.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print command line help
 *
 * parameters:
 *   name  <-- name of executable program
 *----------------------------------------------------------------------------*/

static void
_arg_env_help(const char  *name)
{
  FILE *e = stderr;

  if (cs_glob_rank_id >= 1)
    return;

  fprintf (e, _("Usage: %s [options]\n"), name);

  fprintf (e, _("\nCommand line options:\n\n"));
  fprintf
    (e, _(" --solcom          stand-alone kernel with \"geomet\" mesh in\n"
          "                   SolCom format (obsolete)\n"));
#if defined(HAVE_MPI)
  fprintf
    (e, _(" --mpi             use MPI for parallelism or coupling\n"
          "                   [appnum]: number of this application in\n"
          "                             case of code coupling (default: 0)\n"));
  fprintf
    (e, _(" --mpi-io          <mode> set parallel I/O behavior\n"
          "                     off: do not use MPI-IO\n"
          "                     eo:  MPI-IO with explicit offsets\n"
          "                          (default if available)\n"
          "                     ip:  MPI-IO with individual file pointers\n"));
#endif
  fprintf
    (e, _(" -q, --quality     mesh quality verification mode\n"));
  fprintf
    (e, _(" --benchmark       elementary operations performance\n"
          "                   [--mpitrace] operations done only once\n"
          "                                for light MPI traces\n"));
  fprintf
    (e, _(" --log             output redirection for rank -1 or 0:\n"
          "                     0: standard output\n"
          "                     1: output in \"listing\" (default)\n"));
  fprintf
    (e, _(" --logp            output redirection for rank > 0:\n"
          "                    -1: remove output (default)\n"
          "                     0: no redirection (if independant\n"
          "                        terminals, debugger type)\n"
          "                     1: output in \"listing_n<rang>\"\n"));
#if defined(HAVE_LIBXML2)
  fprintf
    (e, _(" -p, --param       <file_name> parameter file\n"));
#endif

#if defined(HAVE_SOCKET)
  fprintf
    (e, _(" --syr-socket      enable sockets for SYRTHES 3 coupling\n"
          "                   <port_num> port number on rank 0\n"));
#endif

  fprintf
    (e, _(" --version         print version number\n"));
  fprintf
    (e, _(" -h, --help        this help message\n\n"));
}

/*----------------------------------------------------------------------------
 * Print version number
 *----------------------------------------------------------------------------*/

static void
_print_version(void)
{
  if (cs_glob_rank_id >= 1)
    return;

  printf(_("%s version %s\n"), CS_APP_NAME, CS_APP_VERSION);
}

/*----------------------------------------------------------------------------
 * Convert an argument to an integer and check its validity
 *
 * parameters:
 *   arg_id  <-- index of argument in argv
 *   argc    <-- number of command line arguments
 *   argv    <-- array of command line arguments
 *   argerr  --> error indicator
 *
 * returns:
 *   integer value
 *----------------------------------------------------------------------------*/

static int
_arg_to_int(int    arg_id,
            int    argc,
            char  *argv[],
            int   *argerr)
{
  char  *start = NULL;
  char  *end = NULL;
  int  retval = 0;

  *argerr = 0;

  if (arg_id < argc) {
    start = argv[arg_id];
    end = start + strlen(start);
    retval = strtol(start, &end, 0);
    if (end != start + strlen(start)) *argerr = 1;
  }
  else {
    *argerr = 1;
  }

  return retval;
}

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define options and call some associated initializations
 * based on command line arguments
 *
 * parameters:
 *   argc  <-- number of command line arguments
 *   argv  <-- array of command line arguments
 *   opts  --> options structure
 *----------------------------------------------------------------------------*/

void
cs_opts_define(int         argc,
               char       *argv[],
               cs_opts_t  *opts)
{
  /* Local variables */

  const char *s;
  int arg_id = 0, argerr = 0;

  const char moduleoptbase[] = "--yacs-module=";
  const char socketoptbase[] = "--proxy-socket=";
  const char keyoptbase[] = "--proxy-key=";

  /* Default initialization */

  opts->ifoenv = 1;

  opts->ilisr0 = 1;
  opts->ilisrp = 2;

  opts->mpi_io_mode = -1;

  opts->verif = false;
  opts->benchmark = 0;

  opts->syr_socket = -1;

  opts->yacs_module = NULL;

  opts->proxy_socket = NULL;
  opts->proxy_key = -1;

  /* Parse command line arguments */

  while (++arg_id < argc && argerr == 0) {

    s = argv[arg_id];

    if (strcmp(s, "--solcom") == 0)
      opts->ifoenv = 0;

#if defined(HAVE_MPI)

    else if (strcmp(s, "--mpi") == 0) {
      int tmperr = 0;
      (void)_arg_to_int(arg_id + 1, argc, argv, &tmperr);
      if (tmperr == 0) {
        arg_id++;
      }
    }

    else if (strcmp(s, "--mpi-io") == 0) {
      if (arg_id + 1 < argc) {
        const char *s_n = argv[arg_id + 1];
        if (strcmp(s_n, "off") == 0)
          opts->mpi_io_mode = 0;
        else if (strcmp(s_n, "eo") == 0)
          opts->mpi_io_mode = 1;
        else if (strcmp(s_n, "ip") == 0)
          opts->mpi_io_mode = 2;
        else
          argerr = 1;
        if (argerr == 0)
          arg_id++;
      }
      else
        argerr = 1;
    }

#endif /* defined(HAVE_MPI) */

    else if (strcmp(s, "-q") == 0 || strcmp(s, "--quality") == 0)
      opts->verif = true;

    else if (strcmp(s, "--log") == 0) {
      int n1 = 0;
      n1 = _arg_to_int(++arg_id, argc, argv, &argerr);
      if (n1 == 0)
        opts->ilisr0 = 0;
      else if (n1 == 1)
        opts->ilisr0 = 1;
      else
        argerr = 1;
    }

    else if (strcmp(s, "--logp") == 0) {
      int n1 = 0;
      n1 = _arg_to_int(++arg_id, argc, argv, &argerr);
      if (n1 == -1)
        opts->ilisrp = 2;
      else if (n1 == 0)
        opts->ilisrp = 0;
      else if (n1 == 1)
        opts->ilisrp = 1;
      else
        argerr = 1;
    }

    else if (strcmp(s, "--benchmark") == 0) {
      opts->benchmark = 1;
      if (arg_id + 1 < argc) {
        if (strcmp(argv[arg_id + 1], "--mpitrace") == 0) {
          opts->benchmark = 2;
          arg_id++;
        }
      }
    }

#if defined(HAVE_LIBXML2)
    else if (strcmp(s, "-p") == 0 || strcmp(s, "--param") == 0) {
      s = argv[++arg_id];
      argerr = cs_gui_load_file(s);
    }
#endif

#if defined(HAVE_SOCKET)

    /* SYRTHES 3 coupling socket */

    else if (strcmp(s, "--syr-socket") == 0) {
      opts->syr_socket = _arg_to_int(arg_id + 1, argc, argv, &argerr);
      if (argerr == 0)
        arg_id++;
    }

    /* Proxy connection options (do not appear in help as they
       are not destined to be used directly by a user) */

    else if (strncmp(s, socketoptbase, strlen(socketoptbase)) == 0) {
      const char *_s = s + strlen(socketoptbase);
      BFT_MALLOC(opts->proxy_socket, strlen(_s) + 1, char);
      strcpy(opts->proxy_socket, _s);
    }

    else if (strncmp(s, keyoptbase, strlen(keyoptbase)) == 0) {
      const char *_start = s + strlen(keyoptbase);
      char *_end = NULL;
      opts->proxy_key = strtol(_start, &_end, 0);
      if (_end != _start + strlen(_start))
        argerr = 1;
    }

#endif /* defined(HAVE_SOCKET) */

#if defined(HAVE_DLOPEN)

    /* Library loader options (do not appear in help as they
       are not destined to be used directly by a user) */

    else if (strncmp(s, moduleoptbase, strlen(moduleoptbase)) == 0) {
      const char *_s = s + strlen(moduleoptbase);
      BFT_MALLOC(opts->yacs_module, strlen(_s) + 1, char);
      strcpy(opts->yacs_module, _s);
    }

#endif /* defined(HAVE_DLOPEN) */

    /* Version number */

    else if (strcmp(s, "--version") == 0)
      argerr = 3;

    /* Usage */

    else if (strcmp(s, "-h") == 0 || strcmp(s, "--help") == 0)
      argerr = 2;
    else
      argerr = 1;

  } /* End parsing command line */

  /* Print version/help and exit if required or in case of command line error */
  if (argerr != 0) {
    if (cs_glob_rank_id <= 0) {
      switch (argerr) {
      case 1:
      case 2:
        cs_base_logfile_head(argc, argv);
        _arg_env_help(argv[0]);
        break;
      case 3:
        _print_version();
        break;
      default:
        break;
      }
    }
    if (argerr == 1)
      cs_exit(EXIT_FAILURE);
    else
      cs_exit(EXIT_SUCCESS);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
