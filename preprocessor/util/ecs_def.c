/*============================================================================
 * Définitions, variables globales, et fonctions de base
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

#include "cs_config.h"


/* includes système */

#include <assert.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Includes librairie */

#include "ecs_backtrace.h"
#include "ecs_def.h"


/*============================================================================
 *  Constantes et Macros
 *============================================================================*/

/* Type pour la sauvegarde des signaux */

typedef void (*ecs_loc_def__sighandler_t) (int);

/*=============================================================================
 * Définitions de variables globales
 *============================================================================*/

/* Date de compilation */

char  ecs_glob_build_date[]   = __DATE__;

/* Variables globales associée à la gestion des signaux */

#if defined(SIGHUP)
static ecs_loc_def__sighandler_t ecs_loc_def__sighup_sauve  = SIG_DFL;
#endif

static ecs_loc_def__sighandler_t ecs_loc_def__sigabrt_sauve  = SIG_DFL;
static ecs_loc_def__sighandler_t ecs_loc_def__sigint_sauve  = SIG_DFL;
static ecs_loc_def__sighandler_t ecs_loc_def__sigterm_sauve = SIG_DFL;
static ecs_loc_def__sighandler_t ecs_loc_def__sigfpe_sauve  = SIG_DFL;
static ecs_loc_def__sighandler_t ecs_loc_def__sigsegv_sauve = SIG_DFL;

#if defined(SIGXCPU)
static ecs_loc_def__sighandler_t ecs_loc_def__sigcpu_sauve  = SIG_DFL;
#endif

/* Détermination du type d'élément en connectivité nodale en fonction du
   nombre de ses sommets, pour des éléments de dimension 2 (faces)
   ou 3 (cellules) ;
   Au dessus de 8 sommets, on a toujours le type ECS_ELT_TYP_FAC_POLY pour
   les faces et ECS_ELT_TYP_FAC_POLY pour les cellules */

const ecs_elt_typ_t  ecs_glob_typ_elt[2][9] = {

  /* Faces ; au dessus de 5 sommets, on a toujours des polygones */

  {ECS_ELT_TYP_NUL,
   ECS_ELT_TYP_NUL,
   ECS_ELT_TYP_NUL,
   ECS_ELT_TYP_FAC_TRIA,
   ECS_ELT_TYP_FAC_QUAD,
   ECS_ELT_TYP_FAC_POLY,
   ECS_ELT_TYP_FAC_POLY,
   ECS_ELT_TYP_FAC_POLY,
   ECS_ELT_TYP_FAC_POLY},

  /* Cellules ; les polyèdres sont décrits en connectivité descendante,
     avec un sommet supplémentaire par face comme marqueur de fin de face,
     donc un tétraèdre nécessiterait 16 sommets sous cette forme ; tous
     les polyèdres ont donc plus de 8 sommets. */

  {ECS_ELT_TYP_NUL,
   ECS_ELT_TYP_NUL,
   ECS_ELT_TYP_NUL,
   ECS_ELT_TYP_NUL,
   ECS_ELT_TYP_CEL_TETRA,
   ECS_ELT_TYP_CEL_PYRAM,
   ECS_ELT_TYP_CEL_PRISM,
   ECS_ELT_TYP_NUL,
   ECS_ELT_TYP_CEL_HEXA},
};

/*----------------------------------------------------------------------------
 * Fonction de gestion d'un signal fatal (de type SIGFPE ou SIGSEGV)
 *----------------------------------------------------------------------------*/

static void
ecs_loc_def__sig_fatal(int  signum)
{
  fflush(stdout);

  switch (signum) {

#if defined(SIGHUP)
  case SIGHUP:
    fprintf(stderr, _("SIGHUP signal (hang-up) intercepted.\n"
                      "--> computation interrupted.\n"));
    break;
#endif

  case SIGABRT:
    fprintf(stderr, _("SIGABRT signal (abort) intercepted !\n"));
    break;

  case SIGINT:
    fprintf(stderr, _("SIGINT signal (Control+C or equivalent) received.\n"
                      "--> computation interrupted by user.\n"));
    break;

  case SIGTERM:
    fprintf(stderr, _("SIGTERM signal (termination) received.\n"
                      "--> computation interrupted by environment.\n"));
    break;

  case SIGFPE:
    fprintf(stderr, _("SIGFPE signal (floating-point exception) "
                      "intercepted !\n"));
    break;

  case SIGSEGV:
    fprintf(stderr, _("SIGSEGV signal (access to forbidden memory area) "
                      " intercepted !\n"));
    break;

#if defined(SIGXCPU)
  case SIGXCPU:
    fprintf(stderr, _("SIGXCPU signal (CPU time limit exceeded) "
                      "intercepted.\n"));
    break;
#endif

  default:
    fprintf(stderr, _("Signal %d intercepted !\n"), signum);
  }

  ecs_backtrace_print(3);

  assert(0);   /* Use assert to avoit exiting under debugger */

  fflush(stderr);

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------
 *  Fonction qui compte la largeur d'une chaîne de caractères, en tenant
 *  compte de la possibilité que cette chaîne soit de type UTF-8.
 *----------------------------------------------------------------------------*/

static int
ecs_loc_def__strlen(const char  *chaine)
{
  static int mode_utf8 = -1;

  int lng = 0;
  int retval = 0;

  if (mode_utf8 == -1) {
    char *lang = getenv("LANG");
    mode_utf8 = 0;
    if (lang != NULL) {
      if (   strcmp(lang + strlen(lang) - 5, "UTF-8") == 0
          || strcmp(lang + strlen(lang) - 4, "utf8") == 0)
        mode_utf8 = 1;
    }
  }

  if (chaine != NULL) {

    lng = strlen(chaine);

    if (mode_utf8 == 0)
      retval = lng;

    else if (mode_utf8 == 1) {

      int ind;

      for (ind = 0; ind < lng; ind++) {

        char c = chaine[ind];

        if (c < 0x80 || c > 0xBF) { /* Single byte or first byte in UTF-8 */
          retval++;
        }

      }
    }
  }

  return retval;
}

/*============================================================================
 * Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction d'initialisation de la gestion des erreurs
 *----------------------------------------------------------------------------*/

void
ecs_init_gestion_erreur(void)
{
  /* Gestion des signaux */

#if defined(SIGHUP)
  ecs_loc_def__sighup_sauve  = signal(SIGHUP, ecs_loc_def__sig_fatal);
#endif

  ecs_loc_def__sigabrt_sauve = signal(SIGABRT, ecs_loc_def__sig_fatal);
  ecs_loc_def__sigint_sauve  = signal(SIGINT, ecs_loc_def__sig_fatal);
  ecs_loc_def__sigterm_sauve = signal(SIGTERM, ecs_loc_def__sig_fatal);
  ecs_loc_def__sigfpe_sauve  = signal(SIGFPE, ecs_loc_def__sig_fatal);
  ecs_loc_def__sigsegv_sauve = signal(SIGSEGV, ecs_loc_def__sig_fatal);

#if defined(SIGXCPU)
  ecs_loc_def__sigcpu_sauve  = signal(SIGXCPU, ecs_loc_def__sig_fatal);
#endif
}

/*----------------------------------------------------------------------------
 * Fonction d'arrêt
 *----------------------------------------------------------------------------*/

void
ecs_exit(int  statut)
{
  if (statut == EXIT_FAILURE) {

    fprintf (stdout, "\n\n %s\n\n", _("Abnormal end"));
    ecs_backtrace_print(2);

#if defined(DEBUG) || !defined(NDEBUG)
    assert(0);
#endif

  }
  exit (statut);
}

/*----------------------------------------------------------------------------
 * Fonction d'impression d'un avertissement
 *----------------------------------------------------------------------------*/

void
ecs_warn(void)
{
  printf(_("\n"
           "Warning\n"
           "=======\n"));
}

/*----------------------------------------------------------------------------
 * Fonction d'arrêt sur erreur
 *----------------------------------------------------------------------------*/

void
ecs_error(const char  *file_name,
          const int    line_num,
          const int    sys_error_code,
          const char  *format,
          ...)
{
  va_list  arg_ptr;

  fflush(stdout);

  fprintf(stderr, _("\n"
                    "Error in preprocessor execution\n"
                    "===============================\n"));

  if (sys_error_code != 0)
    fprintf(stderr, _("\nSystem error: %s\n"), strerror(sys_error_code));

  fprintf(stderr, _("\n%s:%d: Fatal error.\n\n"), file_name, line_num);

  va_start(arg_ptr, format);
  vfprintf(stderr, format, arg_ptr);
  va_end(arg_ptr);

  fprintf(stderr, "\n\n");

  ecs_backtrace_print(2);

  assert(0);   /* Utilisation de assert pour interception sous un debugger */

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------
 *  Fonction qui imprime une chaîne de caractères avec une largeur
 *  de colonne donnée.
 *----------------------------------------------------------------------------*/

void
ecs_print_padded_str(const char  *str,
                     int          width)
{
  int lng = ecs_loc_def__strlen(str);

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (str != NULL)
    printf("%s", str);

  if (width > lng)
    printf("%-*s", width-lng, "" );
}

/*----------------------------------------------------------------------------*/
