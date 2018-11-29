/*============================================================================
 * Parsing of program arguments and associated initializations
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_mem_usage.h"
#include "bft_printf.h"

#include "cs_base.h"
#include "cs_gui_util.h"
#include "cs_notebook.h"
#include "cs_parameters.h"
#include "cs_partition.h"
#include "cs_system_info.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_opts.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

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

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

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
    (e, _(" --app-name        <app_name> to name this code instance\n"
          "                   (default: working directory base name)\n"));
  fprintf
    (e, _(" --benchmark       elementary operations performance\n"
          "                   [--mpitrace] operations done only once\n"
          "                                for light MPI traces\n"));
  fprintf
    (e, _(" -h, --help        this help message\n\n"));

  fprintf
    (e, _(" --mpi             force use of MPI for parallelism or coupling\n"
          "                   (usually automatic, only required for\n"
          "                   undetermined MPI libraries)\n"));
  fprintf
    (e, _(" --log             output redirection for rank -1 or 0:\n"
          "                     0: standard output\n"
          "                     1: output in \"listing\" (default)\n"));
  fprintf
    (e, _(" --logp            output redirection for rank > 0:\n"
          "                    -1: remove output (default)\n"
          "                     0: no redirection (if independant\n"
          "                        terminals, debugger type)\n"
          "                     1: output in \"listing_n<rank>\"\n"));
  fprintf
    (e, _(" -p, --param       <file_name> parameter file\n"));

  fprintf
    (e, _(" --preprocess      mesh preprocessing mode\n"));

  fprintf
    (e, _(" -q, --quality     mesh quality verification mode\n"));

  fprintf
    (e, _(" --sig-defaults    use default runtime behavior when signals\n"
          "                   are received\n"));

  fprintf
    (e, _(" --system-info     print system information and exit\n"));

  fprintf
    (e, _(" --version         print version number\n"));

#if defined(HAVE_UNISTD_H)
  fprintf
    (e, _(" -wdir, --wdir     <path> working directory\n"));
#endif
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

  /* Default initialization */

  opts->app_name = NULL;

  opts->ilisr0 = 1;
  opts->ilisrp = 2;

  opts->sig_defaults = false;

  opts->preprocess = false;
  opts->verif = false;
  opts->benchmark = 0;

  opts->yacs_module = NULL;

  /* Parse command line arguments */

  while (++arg_id < argc && argerr == 0) {

    s = argv[arg_id];

    if (strcmp(s, "--app-name") == 0) {
      if (arg_id + 1 < argc) {
        BFT_MALLOC(opts->app_name, strlen(argv[arg_id + 1]) + 1, char);
        strcpy(opts->app_name, argv[arg_id + 1]);
        arg_id++;
      }
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

#if defined(HAVE_UNISTD_H)

    else if (strcmp(s, "-wdir") == 0 || strcmp(s, "--wdir") == 0) {
      if (arg_id + 1 < argc) {
        s = argv[++arg_id];
        if (chdir(s) != 0) {
          fprintf(stderr, _("Error switching to directory \"%s\":\n\n"
                            "%s\n"),
                  s, strerror(errno));
          cs_exit(EXIT_FAILURE);
        }
      }
      else
        argerr = 1;
    }

#endif

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

#if defined(HAVE_MPI)

    else if (strcmp(s, "--mpi") == 0) {
      /* Handled in pre-reading stage */
    }

#else /* !defined(HAVE_MPI) */

    else if (strcmp(s, "--mpi") == 0) {
      fprintf(stderr, _("%s was built without MPI support,\n"
                        "so option \"%s\" may not be used.\n"),
              argv[0], s);
      cs_exit(EXIT_FAILURE);
    }

#endif /* defined(HAVE_MPI) */

    else if (strcmp(s, "--preprocess") == 0)
      opts->preprocess = true;

    else if (strcmp(s, "-q") == 0 || strcmp(s, "--quality") == 0)
      opts->verif = true;

    /* Library loader options (do not appear in help as they
       are not destined to be used directly by a user) */

    else if (strncmp(s, moduleoptbase, strlen(moduleoptbase)) == 0) {
      if (cs_glob_rank_id <= 0) {
        const char *_s = s + strlen(moduleoptbase);
#if defined(HAVE_DLOPEN)
        BFT_MALLOC(opts->yacs_module, strlen(_s) + 1, char);
        strcpy(opts->yacs_module, _s);
#else
        fprintf(stderr, _("%s was built without dynamic loader support,\n"
                          "so module file \"%s\" may not be loaded.\n"),
                argv[0], _s);
        cs_exit(EXIT_FAILURE);
#endif /* defined(HAVE_DLOPEN) */
      }
    }

    /* signal handling */

    else if (strcmp(s, "--sig-defaults") == 0)
      opts->sig_defaults = true;

    /* system information */

    else if (strcmp(s, "--system-info") == 0) {
#if defined(HAVE_MPI)
      cs_system_info_no_log(cs_glob_mpi_comm);
#else
      cs_system_info_no_log();
#endif
      cs_partition_external_library_info();
      cs_exit(EXIT_SUCCESS);
    }

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

  /* If application name has not been defined, use working directory
     base name as default. */

  if (opts->app_name == NULL)
    opts->app_name = cs_base_get_app_name(0, NULL);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
