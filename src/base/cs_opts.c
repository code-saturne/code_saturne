/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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

#include <bft_config.h>
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
 *   name  --> name of executable program
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
    (e, _(" -solcom           stand-alone kernel with \"geomet\" mesh in\n"
          "                   SolCom format (obsolete)\n"));
#if defined(HAVE_MPI)
  fprintf
    (e, _(" -mpi, --mpi       use MPI for parallelism or coupling\n"
          "                   [appnum]: number of this application in\n"
          "                             case of code coupling (default: 0)\n"));

#endif
  fprintf
    (e, _(" -q, --quality     mesh quality verification mode\n"));
  fprintf
    (e, _(" -cwf              <criterion> cut warped faces\n"
          "                    -post: activate the post-processing related\n"
          "                           to the cutting of warped faces\n"));
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
 *   arg_id  --> index of argument in argv
 *   argc    --> number of command line arguments
 *   argv    --> array of command line arguments
 *   argerr  <-- error indicator
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

/*----------------------------------------------------------------------------
 * Convert an argument to a double and check its validity
 *
 * parameters:
 *   arg_id  --> index of argument in argv
 *   argc    --> number of command line arguments
 *   argv    --> array of command line arguments
 *   argerr  <-- error indicator
 *
 * returns:
 *   integer value
 *----------------------------------------------------------------------------*/

static double
_arg_to_double(int    arg_id,
               int    argc,
               char  *argv[],
               int   *argerr)
{
  char  *start = NULL;
  char  *end =  NULL;
  double  retval = 0.;

  *argerr = 0;

  if (arg_id < argc) {
    start = argv[arg_id];
    end = start + strlen(start);
    retval = strtod(start, &end);
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
 * Print logfile header
 *
 * parameters:
 *   argc  --> number of command line arguments
 *   argv  --> array of command line arguments
 *----------------------------------------------------------------------------*/

void
cs_opts_logfile_head(int    argc,
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

  bft_printf("\n  Copyright (C) 1998-2008 EDF S.A., France\n\n");

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

/*----------------------------------------------------------------------------
 * First analysis of the command line to determine if we require MPI,
 * and initialization if necessary
 *
 * parameters:
 *   argc  <-> number of command line arguments
 *   argv  <-> array of command line arguments
 *
 * returns:
 *   -1 if MPI is not needed, or application number in MPI_COMM_WORLD of
 *   processes associated with this instance of Code_Saturne
 *----------------------------------------------------------------------------*/

int
cs_opts_mpi_app_num(int    *argc,
                    char  **argv[])
{
#if defined(HAVE_MPI)

  char *s;

  int arg_id = 0, flag = 0;
  int appnum = -1;
  cs_int_t use_mpi = false;

#if defined(MPICH_NAME)

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

#elif   defined(__blrts__) || defined(__bgp__) \
   || defined(__CRAYXT_COMPUTE_LINUX_TARGET)
  use_mpi = true;

#elif defined(LAM_MPI)
  if (getenv("LAMRANK") != NULL)
    use_mpi = true;

#elif defined(OPEN_MPI)
  if (getenv("OMPI_MCA_ns_nds_vpid") != NULL)
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
    if (!flag)
      MPI_Init(argc, argv);
  }

  /* Loop on command line arguments */

  arg_id = 0;

  while (++arg_id < *argc) {

    s = (*argv)[arg_id];

    /* Parallel run */

    if (strcmp(s, "-mpi") == 0 || strcmp(s, "--mpi") == 0) {
      cs_int_t tmperr = 0;
      int _appnum = 0;
      use_mpi = true;
      _appnum = _arg_to_int(arg_id + 1, *argc, *argv, &tmperr);
      if (tmperr == 0) {
        arg_id++;
        appnum = _appnum;
      }
    }

  } /* End of loop on command line arguments */

  if (use_mpi == true) {

    MPI_Initialized(&flag);
    if (!flag)
      MPI_Init(argc, argv);

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
  }

  return appnum;

#else /* if defined(HAVE_MPI) */

  return -1;

#endif
}

/*----------------------------------------------------------------------------
 * Define options and call some associated initializations
 * based on command line arguments
 *
 * parameters:
 *   argc  --> number of command line arguments
 *   argv  --> array of command line arguments
 *   opts  <-- options structure
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

  opts->verif = false;
  opts->benchmark = 0;

  opts->cwf = false;
  opts->cwf_post = false;
  opts->cwf_criterion = 0.01;

  opts->yacs_module = NULL;

  opts->proxy_socket = NULL;
  opts->proxy_key = -1;

  /* Parse command line arguments */

  while (++arg_id < argc && argerr == 0) {

    s = argv[arg_id];

    if (strcmp(s, "-solcom") == 0)
      opts->ifoenv = 0;

#if defined(HAVE_MPI)

    else if (strcmp(s, "-mpi") == 0 || strcmp(s, "--mpi") == 0) {
      cs_int_t tmperr = 0;
      (void)_arg_to_int(arg_id + 1, argc, argv, &tmperr);
      if (tmperr == 0) {
        arg_id++;
      }
    }

#endif /* defined(HAVE_MPI) */

    else if (strcmp(s, "-q") == 0 || strcmp(s, "--quality") == 0)
      opts->verif = true;

    else if (strcmp(s, "--log") == 0) {
      cs_int_t n1 = 0;
      n1 = (cs_int_t) _arg_to_int(++arg_id, argc, argv, &argerr);
      if (n1 == 0)
        opts->ilisr0 = 0;
      else if (n1 == 1)
        opts->ilisr0 = 1;
      else
        argerr = 1;
    }

    else if (strcmp(s, "--logp") == 0) {
      cs_int_t n1 = 0;
      n1 = (cs_int_t) _arg_to_int(++arg_id, argc, argv, &argerr);
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
      argerr = cs_gui_file_loading(s);
    }
#endif

    else if (strcmp(s, "-cwf") == 0) {
      opts->cwf = true;
      if (arg_id + 1 < argc) {
        if (*(argv[arg_id+1]) != '-') {
          opts->cwf_criterion = _arg_to_double(arg_id + 1, argc, argv, &argerr);
          if (argerr == 0)
            arg_id++;
        }
      }
      if (arg_id + 1 < argc) {
        if (strcmp((argv[arg_id+1]), "-post") == 0) {
          opts->cwf_post = true;
          arg_id++;
        }
      }
    }

#if defined(HAVE_SOCKET)

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
        cs_opts_logfile_head(argc, argv);
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
