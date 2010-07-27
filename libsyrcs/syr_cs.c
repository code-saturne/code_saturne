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

/*============================================================================*/
/* SYRTHES wrapper for coupling with Code_Saturne                             */
/*============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_defs.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "syr_defs.h"
#include "syr_coupling.h"
#include "syr_comm.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*=============================================================================
 * Local Structure Definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Prototypes for SYRTHES coupling functions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * SYRTHES initialization for Code_Saturne coupling
 *----------------------------------------------------------------------------*/

void
proc(syrtc1, SYRTC1)(int *ndim_,
                     int *npoinf,
                     int *nodebf,
                     int *nelebf,
                     double *xyzf,
                     double *tf,
                     double *hht);

/*----------------------------------------------------------------------------
 * SYRTHES solver for Code_Saturne coupling
 *----------------------------------------------------------------------------*/

void
proc(syrtc2, SYRTC2)(int *fin,
                     int *npoinf,
                     double *dtfluid,
                     double *tf,
                     double *hht);

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Print help and information on command-line options
 *
 * parameters:
 *   nom    <-- name of executable program
 *   niveau <-- 1: options list; 2: description
 *----------------------------------------------------------------------------*/

static void
syr_cs_loc_aidelc(char  *nom,
                  int    niveau)
{
  FILE *e = stderr;

#ifdef HAVE_MPI
  char opt_mpi[]  = "              [--comm-mpi <a1> <a2> ...]\n";
#else
  char opt_mpi[]  = "";
#endif

#ifdef HAVE_SOCKET
  char opt_sock[]
    = "              [--comm-socket <machine1:port> <machine2:port> ...]\n";
#else
  char opt_sock[] = "";
#endif

  fprintf
    (e,
     "\nUtilisation : %s\n%s%s"
     "              [--app_name nom] [--echo-comm <n>] [-h]\n",
     nom, opt_sock, opt_mpi);

  if (niveau == 2) {
    fprintf
      (e,
       "\nOptions de la ligne de commandes :\n\n");
    fprintf
      (e,
       " --app-name <nom> : nom d'application SYRTHES (defaut:\n"
       "                    nom du repertoire d'execution)\n");
#ifdef HAVE_MPI
    fprintf
      (e,
       " --comm-mpi :       communication par MPI\n"
       "                    <a1 a2 ...> : noms d'applications couplees\n");
#endif
#ifdef HAVE_SOCKET
    fprintf
      (e,
       " --comm-socket <machine:port> : communication par sockets IP\n");
#endif
    fprintf
      (e,
       " --echo-comm <n> :  echo de la communication ;\n"
       "                    -1 : erreur seulement (defaut)\n"
       "                     0 : impression des entetes des messages\n"
       "                     n : impression des entetes des messages ainsi\n"
       "                         que des n premiers et derniers elements\n");
    fprintf
      (e,
       " -h, --help :       appel de l'aide (cet affichage)\n");
    fprintf(e, "\n");
  }
}

/*----------------------------------------------------------------------------*
 * Check presence of string argument
 *
 * parameters:
 *   numarg <-- number of argument to convert
 *   argc   <-- number of arguments in command line
 *   argv   <-- array of command-line arguments
 *   argerr <-- error code
 *----------------------------------------------------------------------------*/

static char *
syr_cs_loc_argstr(int   numarg,
                  int   argc,
                  char *argv[],
                  int  *argerr)
{
  char *retval = NULL;

  if (numarg < argc)
    retval = argv[numarg];
  else
    *argerr = 2;

  return retval;
}

/*----------------------------------------------------------------------------*
 * Convert integer argument to integer with validity check
 *
 * parameters:
 *   numarg <-- number of argument to convert
 *   argc   <-- number of arguments in command line
 *   argv   <-- array of command-line arguments
 *   argerr <-- error code
 *----------------------------------------------------------------------------*/

static int
syr_cs_loc_argint(int   numarg,
                  int   argc,
                  char *argv[],
                  int  *argerr)
{
  char *argdeb;
  char *argfin;
  int   valint = 0;

  if (numarg < argc) {
    argdeb = argv[numarg];
    argfin = argdeb + strlen(argdeb);
    valint = strtol(argdeb, &argfin, 0);
    if (argfin != argdeb + strlen(argdeb)) *argerr = 1;
  }
  else {
    *argerr = 2;
  }

  return valint;
}

/*============================================================================
 * Main function
 *============================================================================*/

int
main(int argc,
     char *argv[])
{
  int  i, i_cas_sat, j_cas_sat;
  int  icoo, isom, ielt;

  char  *s = NULL;

  int  *dernier = NULL;
  int  *fin = NULL;

  int  *idx_som = NULL;
  int  *idx_elt = NULL;

  int  fin_syr = 0;

  int  *_ndim_   = NULL;
  int  *_npoinf  = NULL;
  int  *_nelebf  = NULL;
  int  **_nodebf  = NULL;
  double  **_xyzf = NULL;

  int  ndim_   = 0;
  int  npoinf  = 0;
  int  nelebf  = 0;
  int  *nodebf  = NULL;
  double  *xyzf = NULL;

  double  *tf  = NULL;
  double  *hht = NULL;

  double  *_dtfluid = NULL;
  double dtfluid = -1;

  int fin_total = 0;
  int numarg = 0, argerr = 0;

  /* ------------------------------------------ */
  /* Default initialization of coupling options */
  /* ------------------------------------------ */

  syr_coupling_t **syrcoupl = NULL;

  int nbr_cas_sat = 0;           /* number of Code_Saturne - SYRTHES couplings */
  int echo_comm = -1;            /* Communication verbosity */
  syr_comm_type_t type_comm = SYR_COMM_TYPE_NULL; /* Communication type */
  char **app_sat = NULL;         /* application names of coupled Code_Saturne
                                    instances in the global communicator */
  char **sock_str = NULL;        /* strings for server sockets description */

  const char *app_name = NULL;   /* name of this SYRTHES application */

  /* Initialize error handler */

  syr_errhandler_initialize();

  /* Initialize MPI if necessary (pre-analyze command line) */

#if defined (HAVE_MPI)
  syr_mpi_initialize(&argc, &argv);
  atexit(syr_mpi_exit_force);
#endif

  /* ---------------------------- */
  /* Parse command-line arguments */
  /* ---------------------------- */

  printf("\n*** Interpretation de la ligne de commande ***\n");
  for (numarg = 0; numarg < argc; numarg++)
    printf("%s ",argv[numarg]);
  printf("\n");
  fflush(stdout);

  numarg = 0;

  while (++numarg < argc) {

    s = argv[numarg];

    if (strcmp(s, "--app-name") == 0)
      app_name = syr_cs_loc_argstr(++numarg, argc, argv, &argerr);

#ifdef HAVE_MPI
    else if (strcmp(s, "--comm-mpi") == 0) {
      type_comm = SYR_COMM_TYPE_MPI;

      while (numarg + 1 < argc && *(argv[numarg + 1]) != '-') {
        PLE_REALLOC(app_sat, nbr_cas_sat + 1, char *);
        PLE_REALLOC(sock_str, nbr_cas_sat + 1, char *);
        app_sat[nbr_cas_sat] = syr_cs_loc_argstr(++numarg, argc, argv, &argerr);
        sock_str[nbr_cas_sat] = NULL;
        nbr_cas_sat++;
      }

    }
#endif
#ifdef HAVE_SOCKET
    else if (strcmp(s, "--comm-socket") == 0) {
      type_comm = SYR_COMM_TYPE_SOCKET;

      while (numarg + 1 < argc && *(argv[numarg + 1]) != '-') {
        PLE_REALLOC(app_sat, nbr_cas_sat + 1, char *);
        PLE_REALLOC(sock_str, nbr_cas_sat + 1, char *);
        app_sat[nbr_cas_sat] = NULL;
        PLE_MALLOC(sock_str[nbr_cas_sat], strlen(argv[numarg + 1]) + 1, char);
        strcpy(sock_str[nbr_cas_sat], argv[++numarg]);
        nbr_cas_sat++;
      }

    }
#endif
    else if (strcmp(s, "--echo-comm") == 0 || strcmp(s, "-ec") == 0) {
      if (numarg + 1 < argc && *(argv[numarg + 1]) != '-')
        echo_comm = (int)syr_cs_loc_argint(++numarg, argc, argv, &argerr);
    }
    else if (strcmp(s, "-h") == 0 || strcmp(s, "--help") == 0) {
      syr_cs_loc_aidelc(argv[0], 2);
      syr_exit(EXIT_SUCCESS);
    }
    else
      argerr = 1;
  }

  /* Final default initializations */

  if (echo_comm < -1) argerr = 2;

  /* Print help if necessary */

  if (argerr != 0) {
    syr_cs_loc_aidelc(argv[0], argerr);
    ple_error(__FILE__, __LINE__, 0,
              "Erreur lors de la lecture de la ligne de commande.\n");
  }

  /* ----------------------------------------*/
  /* Initialize of syr_coupling_t structures */
  /* ----------------------------------------*/

  if (echo_comm >= 0) {
    printf
      ("\n*** Initialisation des structures SYRTHES pour le couplage\n");
    fflush(stdout);
  }

  PLE_MALLOC(syrcoupl, nbr_cas_sat, syr_coupling_t *);

  for (i_cas_sat = 0; i_cas_sat < nbr_cas_sat; i_cas_sat++) {

    syrcoupl[i_cas_sat] = syr_coupling_initialize(i_cas_sat,
                                                  app_sat[i_cas_sat],
                                                  sock_str[i_cas_sat],
                                                  type_comm,
                                                  echo_comm);

#if defined (HAVE_SOCKET)
    PLE_FREE(sock_str[i_cas_sat]);
#endif

  }

#if defined (HAVE_SOCKET)
  PLE_FREE(sock_str);
#endif

  PLE_FREE(app_sat);

  /* Allocate working arrays */
  /* ----------------------- */

  PLE_MALLOC(_ndim_,  nbr_cas_sat, int);

  PLE_MALLOC(_npoinf, nbr_cas_sat, int);
  PLE_MALLOC(_nelebf, nbr_cas_sat, int);

  PLE_MALLOC(idx_som, nbr_cas_sat + 1, int);
  PLE_MALLOC(idx_elt, nbr_cas_sat + 1, int);

  for (i = 0; i < nbr_cas_sat; i++) {
    _npoinf[i] = 0;
    _nelebf[i] = 0;
    idx_som[i] = 0;
    idx_elt[i] = 0;
  }

  /* Element connectivity and vertex coordinates */

  PLE_MALLOC(_nodebf, nbr_cas_sat, int *);
  PLE_MALLOC(_xyzf,   nbr_cas_sat, double *);

  for (i = 0; i < nbr_cas_sat; i++) {
    _xyzf[i] = NULL;
    _nodebf[i] = NULL;
  }

  if (echo_comm >= 0) {
    printf("\n*** Reception du maillage couple depuis le(s) noyau(x)\n");
    fflush(stdout);
  }

  /* ----------------------------------------- */
  /* Receive data necessary required by syrtc1 */
  /* ----------------------------------------- */

  idx_som[0] = 0;
  idx_elt[0] = 0;

  for (i_cas_sat = 0; i_cas_sat < nbr_cas_sat; i_cas_sat++) {

    syr_coupling_receive_bc_mesh(syrcoupl[i_cas_sat],
                                 &_ndim_[i_cas_sat],
                                 &_npoinf[i_cas_sat],
                                 &_nelebf[i_cas_sat],
                                 &_xyzf[i_cas_sat],
                                 &_nodebf[i_cas_sat]);

    if (echo_comm >= 0) {
      printf("\n------------------------------------------------\n");
      printf("\tCouplage numero: %d\n", i_cas_sat);
      printf("\tNombre de sommets couples : %9d\n",
             _npoinf[i_cas_sat]);
      printf("\tNombre d'elements couples : %9d\n",
             _nelebf[i_cas_sat]);
      printf("------------------------------------------------\n\n");
      fflush(stdout);
    }

    idx_som[i_cas_sat + 1] = idx_som[i_cas_sat] + _npoinf[i_cas_sat];
    idx_elt[i_cas_sat + 1] = idx_elt[i_cas_sat] + _nelebf[i_cas_sat];

  }

  npoinf = idx_som[nbr_cas_sat];
  nelebf = idx_elt[nbr_cas_sat];

  /* Check that all couplings returned the same value for ndim */

  ndim_ = _ndim_[0];

  /* Append arrays in case of multiple couplings */

  if (nbr_cas_sat > 1) {

    PLE_MALLOC(xyzf, npoinf * ndim_, double);
    PLE_MALLOC(nodebf, nelebf * ndim_, int);

    for (i_cas_sat = 0; i_cas_sat < nbr_cas_sat; i_cas_sat++) {

      for (isom = 0; isom < _npoinf[i_cas_sat]; isom++) {

        for (icoo = 0; icoo < ndim_; icoo++)
          xyzf[isom + idx_som[i_cas_sat] + (npoinf*icoo)]
            = _xyzf[i_cas_sat][isom + (_npoinf[i_cas_sat]*icoo)];

      }

      for (ielt = 0 ; ielt < _nelebf[i_cas_sat] ; ielt++) {

        for (icoo = 0; icoo < ndim_; icoo++) {

          int this_value = _nodebf[i_cas_sat][ielt + (_nelebf[i_cas_sat]*icoo)];

          nodebf[ielt + idx_elt[i_cas_sat] + (nelebf*icoo)]
            = this_value + idx_som[i_cas_sat];

        }

      }

    } /* End of loop on Code_Saturne instance */

    /* Free arrays not needed anymore */

    for (i_cas_sat = 0; i_cas_sat < nbr_cas_sat; i_cas_sat++) {
      PLE_FREE(_xyzf[i_cas_sat]);
      PLE_FREE(_nodebf[i_cas_sat]);
    }

  }
  else {

    xyzf = _xyzf[0];
    nodebf = _nodebf[0];

  }

  PLE_FREE(idx_elt);

  PLE_FREE(_ndim_);
  PLE_FREE(_xyzf);
  PLE_FREE(_nodebf);
  PLE_FREE(_npoinf);
  PLE_FREE(_nelebf);

  /* Allocate arrays for exchanged variables:   */
  /* fluid temperature and exchange coefficient */

  PLE_MALLOC(tf,  npoinf, double);
  PLE_MALLOC(hht, npoinf, double);

  PLE_MALLOC(_dtfluid, nbr_cas_sat, double);

  /* ----------- */
  /* Call syrtc1 */
  /* ----------- */

  proc(syrtc1, SYRTC1)(&ndim_,
                       &npoinf,
                       nodebf,
                       &nelebf,
                       xyzf,
                       tf,
                       hht);

  /* Free nodebf and xyzf, which are not needed anymore */

  PLE_FREE(xyzf);
  PLE_FREE(nodebf);

  /* Prepare for time loop */
  /* --------------------- */

  PLE_MALLOC(dernier, nbr_cas_sat, int);
  PLE_MALLOC(fin,     nbr_cas_sat, int);

  for (i_cas_sat = 0; i_cas_sat < nbr_cas_sat; i_cas_sat++) {
    dernier[i_cas_sat] = 0;
    fin[i_cas_sat] = 0;
  }

  /* ------------------- */
  /* Start time stepping */
  /* ------------------- */

  for (i_cas_sat = 0; i_cas_sat < nbr_cas_sat; i_cas_sat++) {

    syr_coupling_supervise(syrcoupl[i_cas_sat],
                           &dernier[i_cas_sat],
                           &fin[i_cas_sat]);

    /* Send wall temperature and receive
       fluid temperature and exchange coefficient */

    if (fin[i_cas_sat] == 0)
      syr_coupling_exchange_var(syrcoupl[i_cas_sat],
                                tf  + idx_som[i_cas_sat],
                                hht + idx_som[i_cas_sat],
                                _dtfluid + i_cas_sat);

    if (dernier[i_cas_sat] == 1 || fin[i_cas_sat] == 1)
      fin_syr = 1;

  }

  /* ------------- */
  /* Time stepping */
  /* ------------- */

  while (fin_total == 0) {

    /* Handle time step if given by CFD code */

    dtfluid = -1;

    for (i_cas_sat = 0; i_cas_sat < nbr_cas_sat; i_cas_sat++) {
      if (   _dtfluid[i_cas_sat] > 0
          && (dtfluid < 0 || dtfluid > _dtfluid[i_cas_sat]))
        dtfluid = _dtfluid[i_cas_sat];
    }

    /* Call to syrtc2 => SYRTHES calculation in solid mesh */

    proc(syrtc2, SYRTC2)(&fin_syr,
                         &npoinf,
                         &dtfluid,
                         tf,
                         hht);

    for (i_cas_sat = 0; i_cas_sat < nbr_cas_sat; i_cas_sat++) {

      if (fin_syr == 1)
        fin[i_cas_sat] = 1;

      syr_coupling_supervise(syrcoupl[i_cas_sat],
                             &dernier[i_cas_sat],
                             &fin[i_cas_sat]);

      /* Send wall temperature and receive
         fluid temperature and exchange coefficient */

      if (fin[i_cas_sat] == 0)
        syr_coupling_exchange_var(syrcoupl[i_cas_sat],
                                  tf  + idx_som[i_cas_sat],
                                  hht + idx_som[i_cas_sat],
                                  _dtfluid + i_cas_sat);

      if (dernier[i_cas_sat] == 1 || fin[i_cas_sat] == 1)
        fin_syr = 1;

      for (j_cas_sat = 0; j_cas_sat < nbr_cas_sat; j_cas_sat++) {
        if (fin[j_cas_sat] == 1)
          fin_total = 1;
        else {
          fin_total = 0;
          break;
        }
      }

    } /* Loop on couplings with Code_Saturne */

  } /* End of loop on time steps */

  /* ------- */
  /* Cleanup */
  /* ------- */

  for (i_cas_sat = 0; i_cas_sat < nbr_cas_sat; i_cas_sat++) {

    /* Destroy structures used for coupling */

    if (echo_comm >= 0) {
      printf
        ("\n*** Destruction des structures SYRTHES lie au couplage %d\n",
         i_cas_sat);
      fflush(stdout);
    }

    syrcoupl[i_cas_sat] = syr_coupling_finalize(syrcoupl[i_cas_sat]);

  }

  PLE_FREE(syrcoupl);

  PLE_FREE(idx_som);

  PLE_FREE(fin);
  PLE_FREE(dernier);

  PLE_FREE(tf);
  PLE_FREE(hht);

  PLE_FREE(_dtfluid);

  /* Close MPI communications if necessary */

#if defined (HAVE_MPI)
  if (type_comm == SYR_COMM_TYPE_MPI)
    syr_mpi_finalize();
#endif

  /* Normal program exit */

  syr_exit(EXIT_SUCCESS);

  /* The next instruction is never called, but is used to avoid
     compiler warnings about a function returning no value */

  return 0;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
