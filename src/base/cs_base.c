/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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
 * Définitions, variables globales, et fonctions de base
 *============================================================================*/

/* includes système */

#include <assert.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#if defined(_POSIX_SOURCE)
#include <time.h>
#include <unistd.h>
#include <unistd.h>      /* getcwd(), getpid() */
#include <pwd.h>         /* getpwuid()         */
#include <sys/utsname.h>
#endif

/* Includes librairie BFT et FVM */

#include <bft_backtrace.h>
#include <bft_mem_usage.h>
#include <bft_mem.h>
#include <bft_printf.h>
#include <bft_sys_info.h>
#include <bft_timer.h>

#include <fvm_config.h>
#if !defined(_CS_HAVE_MPI) && defined(FVM_HAVE_MPI)
#error "Either both or neither Code_Saturne and FVM must be configured with MPI"
#endif

#include <fvm_parall.h>

/* Includes librairie */

#include "cs_base.h"


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fausse accolade pour corriger l'indentation automatique d'Emacs */
#endif
#endif /* __cplusplus */

/*============================================================================
 *  Constantes et Macros
 *============================================================================*/

/* API Fortran */
/*-------------*/

/*
 * (longueur max 'usuelle' de nom ; un nom plus long est possible,
 * mais provoquera une allocation de mémoire dynamique).
 */

#define CS_BASE_NBR_CHAINE                               5
#define CS_BASE_LNG_CHAINE                               64

/*============================================================================
 * Types privés
 *============================================================================*/

#if defined(_CS_HAVE_MPI)

typedef struct
{
  long val;
  int  rang;
} _cs_base_mpi_long_int_t;

#endif

/* Type pour la sauvegarde des signaux */

typedef void (*_cs_base_sighandler_t) (int);

/*----------------------------------------------------------------------------
 * Prototypes de fonctions privées nécessaires aux variables globales.
 *----------------------------------------------------------------------------*/

static void
_cs_base_exit(int status);

/*============================================================================
 * Variables globales associées a la librairie
 *============================================================================*/

cs_int_t  cs_glob_base_rang = -1;     /* Rang du processus dans le groupe     */
cs_int_t  cs_glob_base_nbr  =  1;     /* Nombre de processus dans le groupe   */

#if defined(_CS_HAVE_MPI)
MPI_Comm  cs_glob_base_mpi_comm = MPI_COMM_NULL;      /* Intra-communicateur  */
#endif

bft_error_handler_t  *cs_glob_base_gest_erreur_sauve = NULL;

/* Variables globales statiques (variables privées de cs_base.c) */

static cs_bool_t  cs_glob_base_bft_mem_init = false;

static cs_bool_t  cs_glob_base_chaine_init = false;
static cs_bool_t  cs_glob_base_chaine_libre[CS_BASE_NBR_CHAINE];
static char       cs_glob_base_chaine[CS_BASE_NBR_CHAINE]
                                     [CS_BASE_LNG_CHAINE + 1];

/* Variables globales associée à la gestion des signaux */

#if defined(SIGHUP)
static _cs_base_sighandler_t cs_glob_base_sighup_sauve = SIG_DFL;
#endif

static _cs_base_sighandler_t cs_glob_base_sigint_sauve = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigterm_sauve = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigfpe_sauve = SIG_DFL;
static _cs_base_sighandler_t cs_glob_base_sigsegv_sauve = SIG_DFL;

#if defined(SIGXCPU)
static _cs_base_sighandler_t cs_glob_base_sigcpu_sauve = SIG_DFL;
#endif

/* Variables globales associées à des pointeurs de fonction */

static cs_exit_t  *_cs_glob_exit = (_cs_base_exit);

/* Variables globales associées à la gestion des tableaux de travail Fortran */

static char      _cs_glob_srt_ia_peak[7] = "------";
static char      _cs_glob_srt_ra_peak[7] = "------";
static cs_int_t  _cs_glob_mem_ia_peak = 0;
static cs_int_t  _cs_glob_mem_ra_peak = 0;
static cs_int_t  _cs_glob_mem_ia_size = 0;
static cs_int_t  _cs_glob_mem_ra_size = 0;

/* Variables globales associées à l'instrumentation */

#if defined(_CS_HAVE_MPI) && defined(_CS_HAVE_MPE)
int  cs_glob_mpe_broadcast_a = 0;
int  cs_glob_mpe_broadcast_b = 0;
int  cs_glob_mpe_synchro_a = 0;
int  cs_glob_mpe_synchro_b = 0;
int  cs_glob_mpe_send_a = 0;
int  cs_glob_mpe_send_b = 0;
int  cs_glob_mpe_rcv_a = 0;
int  cs_glob_mpe_rcv_b = 0;
int  cs_glob_mpe_reduce_a = 0;
int  cs_glob_mpe_reduce_b = 0;
int  cs_glob_mpe_compute_a = 0;
int  cs_glob_mpe_compute_b = 0;
#endif

/*============================================================================
 * Prototypes de fonctions Fortran associées aux impressions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction Fortran d'impression d'un message sur la sortie standard
 *----------------------------------------------------------------------------*/

void CS_PROCF (csprnt, CSPRNT)
(
  char       *cs_buf_print,
  cs_int_t   *msgsize
);


/*----------------------------------------------------------------------------
 * Fonction Fortran de vidage du tampon du fichier d'impression
 *----------------------------------------------------------------------------*/

void CS_PROCF (csflsh, CSFLSH)
(
 void
);


/*----------------------------------------------------------------------------
 * SUBROUTINE CSCLLI : sous-programme de CLose LIsting Fortran
 *----------------------------------------------------------------------------*/

extern void CS_PROCF(csclli, CSCLLI)
(
 void
);


/*============================================================================
 * Définitions de fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction d'impression d'un message sur la sortie standard
 *----------------------------------------------------------------------------*/

static int
_cs_base_bft_printf(const char     *const format,
                    va_list         arg_ptr)
{
 cs_int_t  line;
 cs_int_t  msgsize;

 /* Tampon pour impressions depuis du code C : on imprime dans un chaîne
    de caractères, qui sera imprimée vers un fichier par du code Fortran.
    Une fois les impressions Fortran totalement remplacées par des impressions
    C, on pourra supprimer cette étape, mais elle est nécessaire pour l'instant
    afin de pouvoir utiliser les mêmes fichiers de sortie */

#undef CS_BUF_PRINT_F_SIZE
#define CS_BUF_PRINT_F_SIZE 16384

 static char cs_buf_print_f[CS_BUF_PRINT_F_SIZE];

 /* Impression dans le tampon */

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

  /* Impression effective par le code Fortran */

  CS_PROCF (csprnt, CSPRNT) (cs_buf_print_f, &msgsize);

  return msgsize;
}

/*----------------------------------------------------------------------------
 * Fonction de vidage du tampon d'impression sur la sortie standard
 *----------------------------------------------------------------------------*/

static int
_cs_base_bft_printf_flush(void)
{
  CS_PROCF (csflsh, CSFLSH) ();

  return 0;
}

/*----------------------------------------------------------------------------
 * Fonction d'impression d'un message sur les sorties erreur
 *
 * On répète le message sur la sortie standard et sur un fichier erreur.
 *----------------------------------------------------------------------------*/

static void
_cs_base_err_vprintf(const char  *format,
                     va_list      arg_ptr)
{
  static cs_bool_t  initialise = false;

  /* message sur la sortie standard */

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

  /* message sur une sortie erreur spécifique, à n'initialiser
     que si la sortie erreur est effectivement nécessaire */

  if (initialise == false) {

    char nom_fic_err[81];

    if (cs_glob_base_rang < 1)
      strcpy(nom_fic_err, "error");
    else
      sprintf(nom_fic_err, "error_n%04d", cs_glob_base_rang + 1);

    freopen(nom_fic_err, "w", stderr);

    initialise = true;

  }

  vfprintf(stderr, format, arg_ptr);
}

/*----------------------------------------------------------------------------
 * Fonction d'impression d'un message sur les sorties erreur
 *
 * On répète le message sur la sortie standard et sur un fichier erreur.
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
 * Fonction de sortie
 *----------------------------------------------------------------------------*/

static void
_cs_base_exit(int status)
{
#if defined(_CS_HAVE_MPI)
  {
    int mpi_flag;

    MPI_Initialized(&mpi_flag);

    if (mpi_flag != 0) {

      if (status != EXIT_SUCCESS)
        MPI_Abort(cs_glob_base_mpi_comm, EXIT_FAILURE);

      else { /*  if (status == EXIT_SUCCESS) */

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();

      }
    }
  }
#endif /* _CS_HAVE_MPI */

  exit(status);
}

/*----------------------------------------------------------------------------
 * Fonction d'arret du code en cas d'erreur
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

  _cs_glob_exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------
 * Fonction d'impression d'un "backtrace"
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
 * Fonction de gestion d'un signal fatal (de type SIGFPE ou SIGSEGV)
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

  _cs_glob_exit(EXIT_FAILURE);
}

#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 *  Finalisation MPI
 *----------------------------------------------------------------------------*/

static void
_cs_base_mpi_fin(void)
{
#if defined(_CS_HAVE_MPE)
  if (cs_glob_base_nbr > 1)
    _cs_base_prof_mpe_fin();
#endif

#if defined(FVM_HAVE_MPI)
  fvm_parall_set_mpi_comm(MPI_COMM_NULL);
#endif

  bft_error_handler_set(cs_glob_base_gest_erreur_sauve);

  if (   cs_glob_base_mpi_comm != MPI_COMM_NULL
      && cs_glob_base_mpi_comm != MPI_COMM_WORLD)
    MPI_Comm_free(&cs_glob_base_mpi_comm);
}


#if defined(DEBUG) || !defined(NDEBUG)

/*----------------------------------------------------------------------------
 * Gestionnaire d'erreur MPI
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

  _cs_glob_exit(EXIT_FAILURE);
}

#endif
#endif /* _CS_HAVE_MPI */


#if defined(_CS_HAVE_MPI) && defined(_CS_HAVE_MPE)

/*----------------------------------------------------------------------------
 * Initialisation de l'instrumentation MPE.
 *----------------------------------------------------------------------------*/

static void
_cs_base_prof_mpe_init(int  rang)
{
  int flag;

  MPI_Initialized(&flag);

  /* MPE_Init_log() & MPE_finish_log() ne sont PAS nécessaires si la
     libraire liblmpe.a est "linkée" avec saturne. Dans ce cas, c'est
     MPI_Init() qui se charge de l'appel de MPE_Init_log() */

  if (flag) {

    MPE_Init_log();

    MPE_Log_get_state_eventIDs(&cs_glob_mpe_broadcast_a,
                               &cs_glob_mpe_broadcast_b);
    MPE_Log_get_state_eventIDs(&cs_glob_mpe_synchro_a,
                               &cs_glob_mpe_synchro_b);
    MPE_Log_get_state_eventIDs(&cs_glob_mpe_send_a,
                               &cs_glob_mpe_send_b);
    MPE_Log_get_state_eventIDs(&cs_glob_mpe_rcv_a,
                               &cs_glob_mpe_rcv_b);
    MPE_Log_get_state_eventIDs(&cs_glob_mpe_reduce_a,
                               &cs_glob_mpe_reduce_b);
    MPE_Log_get_state_eventIDs(&cs_glob_mpe_compute_a,
                               &cs_glob_mpe_compute_b);

    if (rang == 0) {

      MPE_Describe_state(cs_glob_mpe_broadcast_a,
                         cs_glob_mpe_broadcast_b,
                         "Broadcast", "orange");
      MPE_Describe_state(cs_glob_mpe_synchro_a,
                         cs_glob_mpe_synchro_b,
                         "MPI Barrier", "blue");
      MPE_Describe_state(cs_glob_mpe_send_a,
                         cs_glob_mpe_send_b,
                         "Send", "yellow");
      MPE_Describe_state(cs_glob_mpe_rcv_a,
                         cs_glob_mpe_rcv_b,
                         "Receive", "red");
      MPE_Describe_state(cs_glob_mpe_reduce_a,
                         cs_glob_mpe_reduce_b,
                         "Reduce", "white");
      MPE_Describe_state(cs_glob_mpe_compute_a,
                         cs_glob_mpe_compute_b,
                         "Compute", "green");

      MPE_Start_log();
      MPE_Log_event(cs_glob_mpe_compute_a, 0, NULL);

    }

  }

}

/*----------------------------------------------------------------------------
 *  Finalisation de l'instrumentation MPE
 *----------------------------------------------------------------------------*/

static void
_cs_base_prof_mpe_fin(void)
{
  int flag, rang;

  MPI_Initialized(&flag);

  if (flag)
    MPI_Comm_rank(MPI_COMM_WORLD, &rang);

  MPE_Log_event(cs_glob_mpe_compute_b, 0, NULL);
  MPE_Log_sync_clocks();

  if (rang == 0)
    MPE_Finish_log("Code_Saturne");
}

#endif /* defined(_CS_HAVE_MPI) && defined(_CS_HAVE_MPE) */

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
#if defined(_CS_HAVE_MPI)

  _cs_base_mpi_long_int_t  val_in, val_max;

  assert(sizeof(double) == sizeof(cs_real_t));

  val_in.val  = *mem_peak;
  val_in.rang = cs_glob_base_rang;

  MPI_Allreduce(&val_in, &val_max, 1, MPI_LONG_INT, MPI_MAXLOC,
                cs_glob_base_mpi_comm);

  *mem_peak = val_max.val;

  MPI_Bcast(srt_peak, 6, MPI_CHAR, val_max.rang, cs_glob_base_mpi_comm);
#endif
}

/*============================================================================
 *  Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction d'arret depuis du code Fortran
 *
 * Interface Fortran :
 *
 * SUBROUTINE CSEXIT (STATUT)
 * *****************
 *
 * INTEGER          STATUT      : --> : 0 pour succès, 1 ou + pour erreur
 *----------------------------------------------------------------------------*/

void CS_PROCF (csexit, CSEXIT)
(
  const cs_int_t  *const statut
)
{
  cs_exit (*statut);
}


/*----------------------------------------------------------------------------
 * Temps CPU écoulé depuis le début de l'exécution
 *
 * Interface Fortran :
 *
 * SUBROUTINE DMTMPS (TCPU)
 * *****************
 *
 * DOUBLE PRECISION TCPU        : --> : temps CPU (utilisateur + système)
 *----------------------------------------------------------------------------*/

void CS_PROCF (dmtmps, DMTMPS)
(
  cs_real_t  *const tcpu
)
{
  *tcpu = bft_timer_cpu_time();
}


/*----------------------------------------------------------------------------
 * Verifier que la réservation en mémoire effectuée par souspg ne dépasse
 * pas  LONGIA.
 *
 * Interface Fortran :
 *
 * SUBROUTINE IASIZE (SOUSPG, MEMINT)
 * *****************
 *
 * CHARACTER*6      SOUSPG      : --> : Nom du sous-programme appelant
 * INTEGER          MEMINT      : --> : Indice de la dernière case utilisée
 *                              :     : dans IA
 *----------------------------------------------------------------------------*/

void CS_PROCF (iasize, IASIZE)
(
 const char   souspg[6],
 cs_int_t    *memint
)
{
  /* Test if enough memory is available */

  if (*memint > _cs_glob_mem_ia_size) {
    char _souspg[7];
    strncpy(_souspg, souspg, 6);
    _souspg[6] = '\0';
    bft_error(__FILE__, __LINE__, 0,
              _(" Sub-routine calling IASIZE:                %s\n"
                " Memory needed in IA (number of integers):  %d\n"
                "         available:                         %d\n\n"
                " ----> Define IASIZE to a value at least equal to %d integers)."),
              _souspg, *memint, _cs_glob_mem_ia_size, *memint);
  }

  /* Update _cs_glob_mem_ia_peak and _cs_glob_srt_ia_peak */

  else if (*memint > _cs_glob_mem_ia_peak) {

    _cs_glob_mem_ia_peak = *memint;
    strncpy(_cs_glob_srt_ia_peak, souspg, 6);
    _cs_glob_srt_ia_peak[6] = '\0';
  }
}

/*----------------------------------------------------------------------------
 * Verifier que la réservation en mémoire effectuée par souspg ne dépasse
 * pas  LONGRA.
 *
 * Interface Fortran :
 *
 * SUBROUTINE RASIZE (SOUSPG, MEMINT)
 * *****************
 *
 * CHARACTER*6      SOUSPG      : --> : Nom du sous-programme appelant
 * INTEGER          MEMRDP      : --> : Indice de la dernière case utilisée
 *                              :     : dans RA
 *----------------------------------------------------------------------------*/

void CS_PROCF (rasize, RASIZE)
(
 const char   souspg[6],
 cs_int_t    *memrdp
)
{
  /* Test if enough memory is available */

  if (*memrdp > _cs_glob_mem_ra_size) {
    char _souspg[7];
    strncpy(_souspg, souspg, 6);
    _souspg[6] = '\0';
    bft_error(__FILE__, __LINE__, 0,
              _(" Sub-routine calling RASIZE:             %s\n"
                " Memory needed in RA (number of reals):  %d\n"
                "         available:                      %d\n\n"
                " ----> Define RASIZE to a value at least equal to %d reals)."),
              _souspg, *memrdp, _cs_glob_mem_ra_size, *memrdp);
  }

  /* Update _cs_glob_mem_ia_peak and _cs_glob_srt_ia_peak */

  else if (*memrdp > _cs_glob_mem_ra_peak) {

    _cs_glob_mem_ra_peak = *memrdp;
    strncpy(_cs_glob_srt_ra_peak, souspg, 6);
    _cs_glob_srt_ra_peak[6] = '\0';
  }
}

/*============================================================================
 * Fonctions publiques
 *============================================================================*/

#if defined(_CS_HAVE_MPI)

/*----------------------------------------------------------------------------
 *  Initialisation MPI ; les variables globales `cs_glob_base_nbr' indiquant
 *  le nombre de processus Code_Saturne et `cs_glob_base_rang' indiquant le
 *  rang du processus courant parmi les processus Code_Saturne sont
 *  (re)positionnées par cette fonction.
 *----------------------------------------------------------------------------*/

void cs_base_mpi_init
(
 int      *argc,          /* --> Nombre d'arguments ligne de commandes        */
 char   ***argv,          /* --> Tableau des arguments ligne de commandes     */
 int       rang_deb       /* --> Rang du premier processus du groupe
                           *     dans MPI_COMM_WORLD                          */
)
{
  int flag, nbr, rang;

#if defined(DEBUG) || !defined(NDEBUG)
  MPI_Errhandler errhandler;
#endif

  MPI_Initialized(&flag);
  if (!flag)
    MPI_Init(argc, argv);

  /*
    Création si nécessaire d'un groupe interne au noyau et éventuellement
    de groupes internes à d'autres processus reliés (opération collective,
    comme toute opération de création de communicateur MPI).

    On suppose que le deuxième argument de MPI_Comm_split (la couleur) est
    égal au rang du premier processus de ce communicateur dans MPI_COMM_WORLD
    (et que ses membres sont contigus dans MPI_COMM_WORLD).
  */

  MPI_Comm_rank(MPI_COMM_WORLD, &rang);

  MPI_Comm_split(MPI_COMM_WORLD, rang_deb, rang - rang_deb + 1,
                 &cs_glob_base_mpi_comm);

  MPI_Comm_size(cs_glob_base_mpi_comm, &nbr);
  MPI_Comm_rank(cs_glob_base_mpi_comm, &rang);

  cs_glob_base_nbr = nbr;

  if (cs_glob_base_nbr > 1)
    cs_glob_base_rang = rang;

  /* Si l'on n'a besoin que de MPI_COMM_WORLD, il est préférable de
     n'utiliser que ce communicateur (potentillement mieux optimisé
     sous certaines architectures) */

  MPI_Comm_size(MPI_COMM_WORLD, &nbr);

  if (cs_glob_base_nbr == 1) {
    MPI_Comm_free(&cs_glob_base_mpi_comm);
    cs_glob_base_mpi_comm = MPI_COMM_NULL;
  }
  else if (nbr == cs_glob_base_nbr) {
    MPI_Comm_free(&cs_glob_base_mpi_comm);
    cs_glob_base_mpi_comm = MPI_COMM_WORLD;
  }

  /* Initialisation des librairies associées */

#if defined(FVM_HAVE_MPI)
  fvm_parall_set_mpi_comm(cs_glob_base_mpi_comm);
  /* IBM Blue Gene ou Cray XT */
#if   defined(__blrts__) || defined(__bgp__) \
   || defined(__CRAYXT_COMPUTE_LINUX_TARGET)
  fvm_parall_set_safe_gather_mode(1);
#endif
#endif

#if defined(DEBUG) || !defined(NDEBUG)
  if (nbr > 1 || cs_glob_base_mpi_comm != MPI_COMM_NULL) {
    MPI_Errhandler_create(&_cs_base_erreur_mpi, &errhandler);
    MPI_Errhandler_set(MPI_COMM_WORLD, errhandler);
    if (   cs_glob_base_mpi_comm != MPI_COMM_WORLD
        && cs_glob_base_mpi_comm != MPI_COMM_NULL)
      MPI_Errhandler_set(cs_glob_base_mpi_comm, errhandler);
    MPI_Errhandler_free(&errhandler);
  }
#endif

#if defined(_CS_HAVE_MPE)
  if (cs_glob_base_nbr > 1)
    _cs_base_prof_mpe_init(rang);
#endif
}

#endif /* _CS_HAVE_MPI */


/*----------------------------------------------------------------------------
 * Fonction d'arrêt
 *----------------------------------------------------------------------------*/

void cs_exit
(
  const cs_int_t  statut
)
{
  if (statut == EXIT_FAILURE) {

    bft_printf_flush();
    bft_backtrace_print(2);

  }

  CS_PROCF(csclli, CSCLLI)(); /* Fermeture des listings */

#if defined(_CS_HAVE_MPI)

  {
    int mpi_flag;

    MPI_Initialized(&mpi_flag);

    if (mpi_flag != 0) {

      if (statut != EXIT_FAILURE) {
        _cs_base_mpi_fin();
      }
    }
  }

#endif /* _CS_HAVE_MPI */

  _cs_glob_exit(statut);
}


/*----------------------------------------------------------------------------
 * Fonction initialisant la gestion des erreurs et des signaux
 *----------------------------------------------------------------------------*/

void cs_base_erreur_init
(
 void
)
{
  /* Gestionnaire d'erreurs */

  cs_glob_base_gest_erreur_sauve = bft_error_handler_get();
  bft_error_handler_set(_cs_base_gestion_erreur);

  /* Gestion des signaux */

  bft_backtrace_print_set(_cs_base_backtrace_print);

#if defined(SIGHUP)
  if (cs_glob_base_rang <= 0)
    cs_glob_base_sighup_sauve  = signal(SIGHUP, _cs_base_sig_fatal);
#endif

  if (cs_glob_base_rang <= 0) {
    cs_glob_base_sigint_sauve  = signal(SIGINT, _cs_base_sig_fatal);
    cs_glob_base_sigterm_sauve = signal(SIGTERM, _cs_base_sig_fatal);
  }

  cs_glob_base_sigfpe_sauve  = signal(SIGFPE, _cs_base_sig_fatal);
  cs_glob_base_sigsegv_sauve = signal(SIGSEGV, _cs_base_sig_fatal);

#if defined(SIGXCPU)
  if (cs_glob_base_rang <= 0)
    cs_glob_base_sigcpu_sauve  = signal(SIGXCPU, _cs_base_sig_fatal);
#endif
}


/*----------------------------------------------------------------------------
 * Fonction initialisant la gestion de contrôle de la mémoire allouée
 *----------------------------------------------------------------------------*/

void cs_base_mem_init
(
 void
)
{
  char  *nom_base;
  char  *nom_complet = NULL;

  /* Initialisation du comptage mémoire */

#if defined(_CS_ARCH_Linux)
  bft_mem_usage_set_options(BFT_MEM_USAGE_TRACK_PR_SIZE |
                            BFT_MEM_USAGE_TRACK_ALLOC_SIZE);
#endif
  bft_mem_usage_init();

  /* Initialisation de la gestion mémoire */

  if ((nom_base = getenv("CS_FIC_MEM")) != NULL) {

    nom_complet = malloc((strlen(nom_base) + 6) * sizeof (char));

    if (nom_complet != NULL) {

      /* En cas de parallélisme, on aura un fichier de trace par processus */
      if (cs_glob_base_rang >= 0)
        sprintf(nom_complet, "%s.%04d", nom_base, cs_glob_base_rang + 1);
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
 * Allocate Fortran work arrays and prepare for their use
 *----------------------------------------------------------------------------*/

void cs_base_mem_init_work
(
 size_t       iasize,
 size_t       rasize,
 cs_int_t   **ia,
 cs_real_t  **ra
)
{
  /* Allocation des tableaux de travail */

  BFT_MALLOC(*ia, iasize, cs_int_t);
  BFT_MALLOC(*ra, rasize, cs_real_t);

  _cs_glob_mem_ia_size = iasize;
  _cs_glob_mem_ra_size = rasize;
}


/*----------------------------------------------------------------------------
 * Fonction terminant la gestion de contrôle de la mémoire allouée
 * et affichant le bilan de la mémoire consommée.
 *----------------------------------------------------------------------------*/

void cs_base_mem_fin
(
 void
)
{
  int        ind_bil, itot;
  cs_real_t  valreal[3];

#if defined(_CS_HAVE_MPI)
  int                imax, imin;
  cs_mpi_real_int_t  val_in[3], val_min[3], val_max[3];
  cs_real_t          val_somme[3];
  cs_int_t           ind_min[3];
#endif

  cs_int_t   ind_val[3] = {1, 1, 1};
  char       unite[]    = {'k', 'm', 'g', 't', 'p'};

  const char  * type_bil[] = {N_("Total memory used:                       "),
                              N_("Memory reported by C library:            "),
                              N_("Theoretical instrumented dynamic memory: ")};

  /* Bilan mémoire */

  bft_printf(_("\nMemory use summary:\n\n"));

  valreal[0] = (cs_real_t) bft_mem_usage_max_pr_size();
  valreal[1] = (cs_real_t) bft_mem_usage_max_alloc_size();
  valreal[2] = (cs_real_t) bft_mem_size_max();

  /* On ignorera les mesures non cohérentes */

  if (valreal[1] < valreal[2])
    ind_val[1] = 0;

  for (ind_bil = 0; ind_bil < 3; ind_bil++) {
    if (valreal[ind_bil] < 1.0)
      ind_val[ind_bil] = 0;
  }

#if defined(_CS_HAVE_MPI)
  if (cs_glob_base_nbr > 1) {
    MPI_Reduce(ind_val, ind_min, 3, CS_MPI_INT, MPI_MIN,
               0, cs_glob_base_mpi_comm);
    MPI_Reduce(valreal, val_somme, 3, CS_MPI_REAL, MPI_SUM,
               0, cs_glob_base_mpi_comm);
    for (ind_bil = 0; ind_bil < 3; ind_bil++) {
      val_in[ind_bil].val = valreal[ind_bil];
      val_in[ind_bil].rang = cs_glob_base_rang;
    }
    MPI_Reduce(&val_in, &val_min, 3, CS_MPI_REAL_INT, MPI_MINLOC,
               0, cs_glob_base_mpi_comm);
    MPI_Reduce(&val_in, &val_max, 3, CS_MPI_REAL_INT, MPI_MAXLOC,
               0, cs_glob_base_mpi_comm);
    if (cs_glob_base_rang == 0) {
      for (ind_bil = 0; ind_bil < 3; ind_bil++) {
        ind_val[ind_bil]  = ind_min[ind_bil];
        valreal[ind_bil] = val_somme[ind_bil];
      }
    }
  }
#endif


  /* Traitement semblable pour les diverses méthodes d'instrumentation */

  for (ind_bil = 0 ; ind_bil < 3 ; ind_bil++) {

    /* Si une méthode d'instrumentation fournit un résultat
       qui semble cohérent, on l'affiche */

    if (ind_val[ind_bil] == 1) {

      for (itot = 0 ;
           valreal[ind_bil] > 1024. && unite[itot] != 'p' ;
           itot++)
        valreal[ind_bil] /= 1024.;
#if defined(_CS_HAVE_MPI)
      if (cs_glob_base_nbr > 1 && cs_glob_base_rang == 0) {
        for (imin = 0 ;
             val_min[ind_bil].val > 1024. && unite[imin] != 'p' ;
             imin++)
          val_min[ind_bil].val /= 1024.;
        for (imax = 0 ;
             val_max[ind_bil].val > 1024. && unite[imax] != 'p' ;
             imax++)
          val_max[ind_bil].val /= 1024.;
      }
#endif

      /* Impressions */

      bft_printf (_("  %s %12.3f %co\n"),
                  _(type_bil[ind_bil]), valreal[ind_bil], unite[itot]);

#if defined(_CS_HAVE_MPI)
      if (cs_glob_base_nbr > 1 && cs_glob_base_rang == 0) {
        bft_printf (_("                             "
                      "local minimum: %12.3f %co  (rank %d)\n"),
                    val_min[ind_bil].val, unite[imin], val_min[ind_bil].rang);
        bft_printf (_("                             "
                      "local maximum: %12.3f %co  (rank %d)\n"),
                    val_max[ind_bil].val, unite[imax], val_max[ind_bil].rang);
      }
#endif
    }

  }

  /* Information sur la gestion des tableaux de travail Fortran */

  if (cs_glob_base_nbr > 1) {
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

#if defined(_CS_HAVE_MPI)
    if (cs_glob_base_nbr > 1) {
      double _wk_size_loc = wk_size[0];
      MPI_Allreduce(&_wk_size_loc, &(wk_size[0]), 1, MPI_DOUBLE, MPI_MAX,
                    cs_glob_base_mpi_comm);
    }
#endif

    for (ind_bil = 0; ind_bil < 2; ind_bil++) {
      for (itot = 0; wk_size[ind_bil] > 1024. && unite[itot] != 'p'; itot++)
        wk_size[ind_bil] /= 1024.;
      wk_unit[ind_bil] = itot;
    }

    bft_printf(_("\n"
                 "  Fortran work arrays memory use:\n"
                 "   %-12lu integers needed (maximum reached in %s)\n"
                 "   %-12lu reals    needed (maximum reached in %s)\n\n"
                 "   Local maximum work memory requested %12.3f %co\n"
                 "                                  used %12.3f %co\n"),
               (unsigned long)_cs_glob_mem_ia_peak, _cs_glob_srt_ia_peak,
               (unsigned long)_cs_glob_mem_ra_peak, _cs_glob_srt_ra_peak,
               wk_size[0], unite[wk_unit[0]],
               wk_size[1], unite[wk_unit[1]]);

  }

  /* Arrêt de la gestion mémoire */

  if (cs_glob_base_bft_mem_init == true)
    bft_mem_end();

  /* Arrêt du comptage mémoire */

  bft_mem_usage_end();
}


/*----------------------------------------------------------------------------
 * Fonction affichant le bilan du temps de calcul et temps écoulé.
 *----------------------------------------------------------------------------*/

void cs_base_bilan_temps
(
 void
)
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


  /* Temps CPU */

  if (utime > 0. || stime > 0.) {
    bft_printf (_("\n  User CPU time:       %12.3f s\n"),
                (float)utime);
    bft_printf (_("  System CPU time:     %12.3f s\n"),
                (float)stime);
  }

  else if (time_cpu > 0.)
    bft_printf (_("\n  CPU time:            %12.3f s\n"),
                (float)time_cpu);

#if defined(_CS_HAVE_MPI)
  if (cs_glob_base_nbr > 1) {
    double time_cumul;
    MPI_Reduce (&time_cpu, &time_cumul, 1, MPI_DOUBLE, MPI_SUM,
                0, cs_glob_base_mpi_comm);
    if (cs_glob_base_rang == 0)
      bft_printf (_("  Total CPU time:      %12.3f s\n"),
                  time_cumul);
  }
#endif


  /* Durée d'exécution  */

  time_tot = bft_timer_wtime();

  if (time_tot > 0.) {

    bft_printf (_("\n  Elapsed time:        %12.3f s\n"),
                time_tot);

    bft_printf (_("  CPU / elapsed time   %12.3f\n"),
                (float)(time_cpu/time_tot));

  }

}


/*----------------------------------------------------------------------------
 * Fonction affichant le bilan du temps de calcul et temps écoulé.
 *----------------------------------------------------------------------------*/

void cs_base_info_systeme
(
 void
)
{
#if defined(_POSIX_SOURCE)

  int             l_user;
  int             l_info;
  time_t          date;
  size_t          ram;
  struct utsname  sys_config;
  struct passwd   *pwd_user;

#undef  _CS_INFO_SYS_STR_SIZE
#define _CS_INFO_SYS_STR_SIZE 80

  char  str_date     [81];
  char  str_system   [_CS_INFO_SYS_STR_SIZE + 1];
  char  str_machine  [_CS_INFO_SYS_STR_SIZE + 1];
  char  str_ram      [_CS_INFO_SYS_STR_SIZE + 1];
  char  str_user     [_CS_INFO_SYS_STR_SIZE + 1];
  char  str_directory[1024];


  /* Date */

  if (time(&date) == -1 ||
      strftime(str_date, _CS_INFO_SYS_STR_SIZE,
               "%c", localtime(&date)) == 0)
    strcpy(str_date, "");


  /* Système et machine */

  if (uname(&sys_config) != -1) {
    strcpy(str_system  , sys_config.sysname );
    strcat(str_system  , " "                );
    strcat(str_system  , sys_config.release );
    strcpy(str_machine , sys_config.nodename);
  }
  else {
    strcpy(str_system  , "");
    strcpy(str_machine , "");
  }

  /* Mémoire vive disponible */

  ram = bft_sys_info_mem_ram();
  if (ram > 1)
    sprintf(str_ram, "%lu", (unsigned long)ram);

  /* Utilisateur */

  /* Fonctions non disponibles sur IBM Blue Gene ou Cray XT */
#if   defined(__blrts__) || defined(__bgp__) \
   || defined(__CRAYXT_COMPUTE_LINUX_TARGET)
  pwd_user = NULL;
#else
  pwd_user = getpwuid(geteuid());
#endif

  if (pwd_user != NULL) {

    str_user[_CS_INFO_SYS_STR_SIZE] = '\0';
    strncpy(str_user, pwd_user->pw_name, _CS_INFO_SYS_STR_SIZE);

    if (pwd_user->pw_gecos != NULL) {

      l_user = strlen(str_user);
      for (l_info = 0;
           (   pwd_user->pw_gecos[l_info] != '\0'
            && pwd_user->pw_gecos[l_info] != ',');
           l_info++);

      if (l_user + l_info + 3 < _CS_INFO_SYS_STR_SIZE) {
        strcat(str_user, " (");
        strncpy(str_user + l_user + 2, pwd_user->pw_gecos, l_info);
        str_user[l_user + 2 + l_info]     = ')';
        str_user[l_user + 2 + l_info + 1] = '\0';
      }

    }

  }
  else
    strcpy(str_user, "");

  /* Répertoire courant */

  if (getcwd(str_directory, 1024) == NULL)
    strcpy(str_directory, "");


  /* Affichage de la configuration locale */
  /*--------------------------------------*/

  bft_printf("\n%s\n", _("Local case configuration:\n"));

  bft_printf("  %s%s\n", _("Date:              "), str_date);

  bft_printf("  %s%s\n", _("System:            "),  str_system);
  bft_printf("  %s%s\n", _("Machine:           "),  str_machine);
  bft_printf("  %s%s\n", _("Processor:         "),  bft_sys_info_cpu());
  if (ram > 0)
    bft_printf("  %s%s\n", _("Memory:            "),   str_ram);
  bft_printf("  %s%s\n", _("User:                "), str_user);
  bft_printf("  %s%s\n", _("Directory:           "),  str_directory);

  bft_printf("\n");

#undef  _CS_INFO_SYS_STR_SIZE

#endif /* _POSIX_SOURCE */
}


/*----------------------------------------------------------------------------
 * Modification du comportement des fonctions bft_printf() par défaut
 *----------------------------------------------------------------------------*/

void cs_base_bft_printf_set
(
 void
)
{
  bft_printf_proxy_set(_cs_base_bft_printf);
  bft_printf_flush_proxy_set(_cs_base_bft_printf_flush);
}


/*----------------------------------------------------------------------------
 * Fonction d'impression d'un message "avertissement"
 *----------------------------------------------------------------------------*/

void cs_base_warn
(
 const char  *file_name,
 const int    line_num
)
{
  bft_printf(_("\n\nCode_Saturne: %s:%d: Warning\n"),
             file_name, line_num);
}


/*----------------------------------------------------------------------------
 * Conversion d'une chaîne de l'API Fortran vers l'API C,
 * (avec suppression des blancs en début ou fin de chaîne).
 *----------------------------------------------------------------------------*/

char  * cs_base_chaine_f_vers_c_cree
(
 const char      * f_str,                   /* --> Chaîne Fortran             */
 int               f_len                    /* --> Longueur de la chaîne      */
)
{
  char * c_str = NULL;
  int    i, i1, i2, l;

  /* Initialisation si nécessaire */

  if (cs_glob_base_chaine_init == false) {
    for (i = 0 ; i < CS_BASE_NBR_CHAINE ; i++)
      cs_glob_base_chaine_libre[i] = true;
    cs_glob_base_chaine_init = true;
  }

  /* Traitement du nom pour l'API C */

  for (i1 = 0 ;
       i1 < f_len && (f_str[i1] == ' ' || f_str[i1] == '\t') ;
       i1++);

  for (i2 = f_len - 1 ;
       i2 > i1 && (f_str[i2] == ' ' || f_str[i2] == '\t') ;
       i2--);

  l = i2 - i1 + 1;

  /* Allocation si nécessaire */

  if (l < CS_BASE_LNG_CHAINE) {
    for (i = 0 ; i < CS_BASE_NBR_CHAINE ; i++) {
      if (cs_glob_base_chaine_libre[i] == true) {
        c_str = cs_glob_base_chaine[i];
        cs_glob_base_chaine_libre[i] = false;
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
 *  Libération d'une chaîne convertie de l'API Fortran vers l'API C
 *----------------------------------------------------------------------------*/

char  * cs_base_chaine_f_vers_c_detruit
(
 char  * c_str                              /* --> Chaîne C                   */
)
{
  cs_int_t ind;

  for (ind = 0 ; ind < CS_BASE_NBR_CHAINE ; ind++) {
    if (c_str == cs_glob_base_chaine[ind]) {
      cs_glob_base_chaine_libre[ind] = true;
      c_str = NULL;
      break;
    }
  }

  if (ind == CS_BASE_NBR_CHAINE)
    BFT_FREE(c_str);

  return c_str;
}


/*----------------------------------------------------------------------------
 * Modification du comportement de sortie du programme par défaut
 *----------------------------------------------------------------------------*/

void
cs_base_exit_set(cs_exit_t *exit_func)
{
  _cs_glob_exit = exit_func;
}


/*============================================================================
 * Fonctions privées
 *============================================================================*/
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
