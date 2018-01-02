/*============================================================================
 * Methods for particle restart reading
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

/*============================================================================
 * Functions dealing with restart reading
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "cs_log.h"

#include "cs_field.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_parameters_check.h"
#include "cs_physical_model.h"
#include "cs_restart.h"
#include "cs_restart_default.h"

#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_post.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_restart.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_lec.h"

BEGIN_C_DECLS

/*============================================================================
 * Fortran wrapper function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *\brief Fortran wrapper for restart files reading.
 */
/*----------------------------------------------------------------------------*/

void
CS_PROCF(laglec, LAGLEC)(void)
{
  cs_restart_lagrangian_checkpoint_read();
}

/*--------------------------------------------------------------------*/
/*! \brief Fortran wrapper for restart files output.
 */
/*--------------------------------------------------------------------*/

void
CS_PROCF(lagout, LAGOUT)(void)
{
  cs_restart_lagrangian_checkpoint_write();
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *\brief  Read Lagrangian restart files.
 */
/*----------------------------------------------------------------------------*/

void
cs_restart_lagrangian_checkpoint_read(void)
{
  cs_lnum_t mstits, musbor, itysup, nbval, nberro;
  cs_lnum_t mstist, jdstnt, nclsto, mstbor, jsttio;
  cs_lnum_t jturb , jtytur;

  cs_int_t  ierror = CS_RESTART_SUCCESS;

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  int nvplmx = 50+4*cs_glob_lagr_const_dim->nlayer;

  typedef char cs_char_64_t[64];

  cs_lnum_t nfabor      = cs_glob_mesh->n_b_faces;
  cs_lnum_t ncel        = cs_glob_mesh->n_cells;
  cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

  /* Default initializations */

  if (cs_glob_lagr_time_scheme->iilagr == 2) {

    for (int ivar = 0; ivar < cs_glob_lagr_dim->ntersl; ivar++) {

      cs_real_t *st_val = cs_glob_lagr_source_terms->st_val + ivar*n_cells_ext;

      for (cs_lnum_t iel = 0; iel < ncel; iel++)
        st_val[iel] = 0.;

    }

  }

  if (cs_glob_lagr_dim->nvisbr > 0) {

    for (cs_lnum_t ivar = 0; ivar < cs_glob_lagr_dim->nvisbr; ivar++) {

      for (cs_lnum_t ifac = 0; ifac < nfabor; ifac++)
        bound_stat[ifac + nfabor * ivar]   = 0.0;

    }

  }

  if (cs_glob_lagr_time_scheme->isuila == 0)
    return;

  cs_char_64_t *nomtsl = NULL;

  BFT_MALLOC(nomtsl, nvplmx, cs_char_64_t);

  /* ====================================================================   */
  /* 2. LECTURE DU FICHIER SUITE STATISTIQUES ET TERMES SOURCES   */
  /*    DE COUPLAGE RETOUR     */
  /* ====================================================================   */

  static cs_restart_t  *cs_lag_stat_restart = NULL;

  if (cs_glob_lagr_stat_options->isuist == 1) {

    cs_log_printf
      (CS_LOG_DEFAULT,
       "   ** INFORMATION ON THE LAGRANGIAN COMPUTATION\n"
       "-------------------------------------\n"
       "   Read restart file for statistics and return coupling source terms\n");

    char const *ficsui = "lagrangian_stats";
    cs_lag_stat_restart = cs_restart_create(ficsui, NULL, CS_RESTART_MODE_READ);
    if (cs_lag_stat_restart == NULL)
      bft_error(__FILE__, __LINE__, 0,
                _("Abort while opening the lagrangian module restart file in read mode.\n"
                  "Verify the existence and the name of the restart file: %s\n"),
                ficsui);

    cs_log_printf(CS_LOG_DEFAULT,
                  "      Debut de la lecture\n");

    /*  ---> Type de fichier suite    */
    /*        Pourrait porter le numero de version si besoin.  */
    /*        On ne se sert pas de IVERS pour le moment   */
    {
      itysup = 0;
      nbval  = 1;
      char rubriq[] = "version_fichier_suite_Lagrangien_statistiques";
      cs_int_t *tabvar;

      BFT_MALLOC(tabvar, nbval, cs_int_t);

      ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq,
                                       itysup, nbval, CS_TYPE_cs_int_t, tabvar);
      BFT_FREE(tabvar);
      if (ierror != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("Abort while opening the lagrangian module restart file: %s\n"
                    "This file does not seem to be a Lagrangian checkpoint file."),
                  ficsui);
    }

    {
      itysup = 0;
      nbval  = 1;
      char rubriq[] = "indicateur_ecoulement_stationnaire";
      cs_int_t *tabvar;

      BFT_MALLOC(tabvar, nbval, cs_int_t);

      ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq,
                                       itysup, nbval, CS_TYPE_cs_int_t, tabvar);
      jsttio = tabvar[0];
      BFT_FREE(tabvar);

      if (ierror != 0)
        bft_error(__FILE__, __LINE__, 0,
                  _("Abort while opening the lagrangian module restart file: %s\n"
                    "Error reading section: %s"),
                  ficsui, rubriq);
    }

    /* Mesh dimensions*/
    bool ncelok, nfaiok, nfabok, nsomok;
    cs_restart_check_base_location (cs_lag_stat_restart, &ncelok, &nfaiok,
                                    &nfabok, &nsomok);

    if (! ncelok)
      cs_parameters_error
        (CS_ABORT_DELAYED,
         _("in Lagrangian module"),
         _("The number of cells in restart file: %s\n"
           "is different from that of the current mesh.\n"),
         cs_restart_get_name(cs_lag_stat_restart));

    if (! nfaiok)
      cs_parameters_error
        (CS_WARNING,
         _("in Lagrangian module"),
         _("The number of interior faces in restart file: %s\n"
           "is different from that of the current mesh.\n\n"
           "interior face data will be reinitialized.\n"),
         cs_restart_get_name(cs_lag_stat_restart));

    if (! nfabok)
      cs_parameters_error
        (CS_WARNING,
         _("in Lagrangian module"),
         _("The number of boundary faces in restart file: %s\n"
           "is different from that of the current mesh.\n\n"
           "boundary face data will be reinitialized.\n"),
         cs_restart_get_name(cs_lag_stat_restart));

    /* --> Est-on cense lire une suite de stats volumiques ?   */
    if (cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->idstnt) {

      nberro  = 0;

      {
        itysup  = 0;
        nbval   = 1;
        cs_int_t *tabvar;

        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[] = "iteration_debut_statistiques";
        ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup,
                                         nbval, CS_TYPE_cs_int_t, tabvar);
        jdstnt  = tabvar[0];
        nberro += ierror;
        BFT_FREE(tabvar);
      }

      {
        itysup  = 0;
        nbval   = 1;
        cs_int_t *tabvar;

        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[]  = "iteration_debut_statistiques_stationnaires";
        ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup,
                                         nbval, CS_TYPE_cs_int_t, tabvar);
        mstist  = tabvar[0];
        nberro += ierror;
        BFT_FREE(tabvar);
      }

      /*  ---> S'il y a des erreurs, on suppose que c'est parce que le fichier  */
      /*         suite ne contient pas d'infos sur les stats volumiques.   */
      /*         Dans ce cas, si on est en instationnaire on se dit que c'est   */
      /*         pas grave, on saute l'etape et on continue. Par contre si on   */
      /*         est dans une configuration de calcul de stats volumiques  */
      /*         en stationnaire on stoppe.  */

      if (nberro != 0) {

        if (   cs_glob_lagr_time_scheme->isttio == 0
            || (   cs_glob_lagr_time_scheme->isttio == 1
                && cs_glob_time_step->nt_cur < cs_glob_lagr_stat_options->nstist))
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      LES STATISTIQUES VOLUMIQUES DU CALCUL AMONT\n"
             "@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES\n"
             "@        DU FICHIER SUITE\n"
             "@\n"
             "@    Les indicateurs concernant le calcul\n"
             "@      instationnaire/stationnaire des\n"
             "@      statistiques volumiques sont positionnes en mode\n"
             "@      instationnaire ou en debut de calcul stationnaire :\n"
             "@\n"
             "@          ISTTIO    IDSTNT    NSTIST    Iter de redemarrage\n"
             "@      %10d%10d%10d%10d\n"
             "@\n"
             "@    Elles seront initialisees par des valeurs par defaut.\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui, cs_glob_lagr_time_scheme->isttio,
             cs_glob_lagr_stat_options->idstnt,
             cs_glob_lagr_stat_options->nstist,
             cs_glob_time_step->nt_cur+1);

        else
          bft_error
            (__FILE__, __LINE__, 0,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      LES STATISTIQUES VOLUMIQUES DU CALCUL AMONT\n"
             "@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES\n"
             "@        DU FICHIER SUITE\n"
             "@\n"
             "@    Les indicateurs concernant le calcul\n"
             "@      instationnaire/stationnaire des\n"
             "@      statistiques volumiques sont positionnes\n"
             "@      en mode stationnaire :\n"
             "@\n"
             "@          ISTTIO    IDSTNT    NSTIST    Iter de redemarrage\n"
             "@      %10d%10d%10d%10d\n"
             "@\n"
             "@    Le calcul ne peut pas etre execute.\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui, cs_glob_lagr_time_scheme->isttio,
             cs_glob_lagr_stat_options->idstnt,
             cs_glob_lagr_stat_options->nstist,
             cs_glob_time_step->nt_cur+1);

      }
      /* --> A partir d'ici on considere que le fichier suite contient     */
      /*       des stats volumiques     */
      else {

        if (   jsttio != cs_glob_lagr_time_scheme->isttio
            || jdstnt != cs_glob_lagr_stat_options->idstnt
            || mstist != cs_glob_lagr_stat_options->nstist)
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
             "@\n"
             "@    Les indicateurs concernant le calcul\n"
             "@      instationnaire/stationnaire des\n"
             "@      statistiques volumiques sont modifies :\n"
             "@\n"
             "@              ISTTIO    IDSTNT    NSTIST\n"
             "@  AMONT : %10d%10d%10d'\n"
             "@  ACTUEL: %10d%10d%10d\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             jsttio,
             jdstnt,
             mstist,
             cs_glob_lagr_time_scheme->isttio,
             cs_glob_lagr_stat_options->idstnt,
             cs_glob_lagr_stat_options->nstist);

        {
          itysup  = 0;
          nbval   = 1;
          cs_int_t *tabvar;

          BFT_MALLOC(tabvar, nbval, cs_int_t);

          char rubriq[] = "classe_statistique_particules";
          ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup,
                                           nbval, CS_TYPE_cs_int_t, tabvar);
          nclsto  = tabvar[0];
          BFT_FREE(tabvar);
        }
        if (ierror != 0)
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
             "@ %s\n"
             "@\n"
             "@    Le mot cle est initialise avec sa valeur par defaut\n"
             "@      ou celle donnee dans le sous-programme cs_user_lagr_model.c :\n"
             "@        %s  = %d\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             "classes_statistiques",
             "NBCLST",
             nclsto);

        /*  --> Verif de coherence de l'avancement du calcul avec les   */
        /*       indicateurs de calcul de la suite actuelle : */

        if (cs_glob_lagr_model->n_stat_classes != nclsto)
          bft_error
            (__FILE__, __LINE__, 0,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
             "@\n"
             "@    LE CALCUL SE POURSUIT AVEC UN CALCUL DE\n"
             "@      STATISTIQUES VOLUMIQUES EN MODE STATIONNAIRE\n"
             "@      MAIS LES INDICATEURS DE CONTROLES DES STATISTIQUES\n"
             "@      ON ETE MODIFIEES.\n"
             "@\n"
             "@    Pour eviter les incoherences dans le calcul\n"
             "@      NBCLST ne devrait pas etre modifies entre\n"
             "@      deux calculs de statistiques volumiques.\n"
             "@\n"
             "@    Le calcul ne peut pas etre execute.\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui);

      }

    }

    if (   cs_glob_lagr_dim->nvisbr > 0
        && nfabok) {

      {
        itysup  = 0;
        nbval   = 1;

        cs_int_t tabvar;

        char rubriq[]  = "iteration_debut_stats_frontieres_stationnaires";
        ierror = cs_restart_read_section(cs_lag_stat_restart,
                                         rubriq, itysup, nbval, CS_TYPE_cs_int_t,
                                         &tabvar);
        mstbor = tabvar;
      }

      /*  ---> S'il y a une erreur, on suppose que c'est parce que le fichier   */
      /*         suite ne contient pas d'infos sur les stats aux frontieres.    */
      /*         Dans ce cas, si on est en instationnaire on se dit que c'est   */
      /*         pas grave, on saute l'etape et on continue. Par contre si on   */
      /*         est dans une configuration de calcul de stats aux frontieres   */
      /*         en stationnaire on stoppe.  */

      if (ierror != 0) {

        if (   cs_glob_lagr_time_scheme->isttio == 0
            || (   cs_glob_lagr_time_scheme->isttio == 1
                &&   cs_glob_time_step->nt_cur
                   < cs_glob_lagr_stat_options->nstist)) {

          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      LES STATISTIQUES AUX FRONTIERES DU CALCUL AMONT\n"
             "@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES\n"
             "@        DU FICHIER SUITE\n"
             "@\n"
             "@    Les indicateurs concernant le calcul\n"
             "@      instationnaire/stationnaire des\n"
             "@      statistiques aux frontieres sont positionnes en mode\n"
             "@      instationnaire ou en debut de calcul stationnaire :\n"
             "@\n"
             "@          ISTTIO    NSTBOR    Iter de redemarrage\n"
             "@      %10d%10d%10d\n"
             "@\n"
             "@    Elles seront initialisees par des valeurs par defaut.\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             cs_glob_lagr_time_scheme->isttio,
             cs_glob_lagr_stat_options->nstist,
             cs_glob_time_step->nt_cur+1);

        }
        else
          bft_error
            (__FILE__, __LINE__, 0,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      LES STATISTIQUES AUX FRONTIERES DU CALCUL AMONT\n"
             "@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTES\n"
             "@        DU FICHIER SUITE\n"
             "@\n"
             "@    Les indicateurs concernant le calcul\n"
             "@      instationnaire/stationnaire des\n"
             "@      statistiques aux frontieres sont positionnes\n"
             "@      en mode stationnaire :\n"
             "@\n"
             "@          ISTTIO    NSTBOR    Iter de redemarrage\n"
             "@      %10d%10d%10d\n"
             "@\n"
             "@    Le calcul ne peut pas etre execute.\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             cs_glob_lagr_time_scheme->isttio,
             cs_glob_lagr_stat_options->nstist,
             cs_glob_time_step->nt_cur+1);

      }
      /* --> A partir d'ici on considere que le fichier suite contient     */
      /*       des stats volumiques     */
      else {

        if (   jsttio != cs_glob_lagr_time_scheme->isttio
            || mstbor != cs_glob_lagr_stat_options->nstist)
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
             "@\n"
             "@    Les indicateurs concernant le calcul\n"
             "@      instationnaire/stationnaire des\n"
             "@      statistiques aux frontieres sont modifies :\n"
             "@\n"
             "@                        ISTTIO    NSTBOR\n"
             "@            AMONT : %10d%10d\n"
             "@            ACTUEL: %10d%10d\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             jsttio,
             mstbor,
             cs_glob_lagr_time_scheme->isttio,
             cs_glob_lagr_stat_options->nstist);

        /*  --> Lecture de l'avancement du calcul stats aux frontieres  */
        {
          char rubriq[] = "nombre_iterations_stats_frontieres";
          cs_int_t *tabvar;

          BFT_MALLOC(tabvar, nbval, cs_int_t);
          ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                           CS_TYPE_cs_int_t, tabvar);
          cs_glob_lagr_boundary_interactions->npstft  = tabvar[0];
          BFT_FREE(tabvar);
        }
        if (ierror != 0)
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
             "@ %s\n"
             "@\n"
             "@    Le mot cle est initialise avec sa valeur par defaut\n"
             "@      ou celle donnee dans le sous-programme cs_user_lagr_model :\n"
             "@        %s  = %d\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             "nombre_iterations_stats_frontieres",
             "NPSTFT",
             cs_glob_lagr_boundary_interactions->npstft);

        {
          char rubriq[]  = "nombre_iterations_stats_frontieres_stationnaires";
          cs_int_t *tabvar;

          BFT_MALLOC(tabvar, nbval, cs_int_t);
          ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                           CS_TYPE_cs_int_t, tabvar);
          cs_glob_lagr_boundary_interactions->npstf   = tabvar[0];
          BFT_FREE(tabvar);
        }
        if (ierror != 0)
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
             "@ %s\n"
             "@\n"
             "@    Le mot cle est initialise avec sa valeur par defaut\n"
             "@      ou celle donnee dans le sous-programme cs_user_lagr_model :\n"
             "@        %s  = %d\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             "nombre_iterations_stats_frontieres_stationnaires",
             "NPSTF",
             cs_glob_lagr_boundary_interactions->npstf);

        {
          char rubriq[] = "temps_stats_frontieres_stationnaires";
          cs_real_t *tabvar;

          BFT_MALLOC(tabvar, nbval, cs_real_t);
          ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                           CS_TYPE_cs_real_t, tabvar);
          cs_glob_lagr_boundary_interactions->tstatp  = tabvar[0];
          BFT_FREE(tabvar);
        }
        if (ierror != 0)
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
             "@ %s\n"
             "@\n"
             "@    Le mot cle est initialise avec sa valeur par defaut\n"
             "@      ou celle donnee dans le sous-programme cs_user_lagr_model :\n"
             "@        %s  = %14.5E\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             "temps_stats_frontieres_stationnaires",
             "TSTATP",
             cs_glob_lagr_boundary_interactions->tstatp);

        /*  --> Verif de coherence de l'avancement du calcul avec les   */
        /*       indicateurs de calcul de la suite actuelle : */

        if (   cs_glob_lagr_boundary_interactions->npstf == 0
            && (   cs_glob_lagr_time_scheme->isttio == 1
                &&    cs_glob_lagr_stat_options->nstist
                   <= cs_glob_time_step->nt_cur))
          bft_error
            (__FILE__, __LINE__, 0,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
             "@\n"
             "@    LE CALCUL ACTUEL DES STATISTIQUES AUX FRONTIERES\n"
             "@      EST EN MODE STATIONNAIRE, ALORS QUE LE FICHIER\n"
             "@      SUITE CONTIENT DES STATISTIQUES INSTATIONNAIRES.\n"
             "@\n"
             "@    NSTBOR devrait etre un entier superieur ou egal\n"
             "@      a l'iteration Lagrangienne absolue de redemarrage\n"
             "@      du calcul (iteration : %d)\n"
             "@\n"
             "@      Il vaut ici NSTBOR = %d\n"
             "@\n"
             "@    Le calcul ne peut pas etre execute.\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             cs_glob_time_step->nt_cur+1,
             cs_glob_lagr_stat_options->nstist);

        if (   cs_glob_lagr_boundary_interactions->npstf > 0
            && (   (   cs_glob_lagr_time_scheme->isttio == 1
                    &&    cs_glob_time_step->nt_cur
                       <= cs_glob_lagr_stat_options->nstist)
                || cs_glob_lagr_time_scheme->isttio == 0))
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
             "@\n"
             "@    LE CALCUL ACTUEL DES STATISTIQUES AUX FRONTIERES\n"
             "@      EST EN MODE INSTATIONNAIRE, ALORS QUE LE FICHIER\n"
             "@      SUITE CONTIENT DES STATISTIQUES STATIONNAIRES.\n"
             "@\n"
             "@    Les statistiques aux frontieres stationnaires amont\n"
             "@      seront remises a zero.\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui);

        if (   cs_glob_lagr_boundary_interactions->npstf > 0
            && (   cs_glob_lagr_time_scheme->isttio == 1
                &&    cs_glob_time_step->nt_cur
                   >= cs_glob_lagr_stat_options->nstist)) {

          if (mstbor != cs_glob_lagr_stat_options->nstist)
            bft_error
              (__FILE__, __LINE__, 0,
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
               "@    =========     LAGRANGIEN %s\n"
               "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
               "@\n"
               "@    LE CALCUL SE POURSUIT AVEC UN CALCUL DE\n"
               "@      STATISTIQUES AUX FRONTIERES EN MODE STATIONNAIRE\n"
               "@      MAIS LES INDICATEURS DE CONTROLES DES STATISTIQUES\n"
               "@      ON ETE MODIFIEES.\n"
               "@\n"
               "@    Pour eviter les incoherences dans le calcul\n"
               "@      NSTBOR ne devrait pas etre modifie entre deux calculs\n"
               "@      de statistiques aux frontieres stationnaires.\n"
               "@\n"
               "@    Le calcul ne peut pas etre execute.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               ficsui);

        }

        /* --> Stats supplementaires utilisateurs   */
        {
          cs_int_t *tabvar;
          BFT_MALLOC(tabvar, nbval, cs_int_t);
          char rubriq[] = "nombre_stats_frontieres_utilisateur";
          ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval, CS_TYPE_cs_int_t, tabvar);
          musbor  = tabvar[0];
          BFT_FREE(tabvar);
        }
        if (cs_glob_lagr_boundary_interactions->nusbor < musbor)
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
             "@\n"
             "@    L'indicateur du  nombre de statistiques aux frontieres\n"
             "@      supplementaires utilisateur est modifie,\n"
             "@      ou n'a pas pu etre relu.\n"
             "@\n"
             "@              NUSBOR\n"
             "@    AMONT : %d      ACTUEL : %d\n"
             "@\n"
             "@    Si ACTUEL > AMONT, on initialise les %d 1eres\n"
             "@      statistiques supplementaires actuelles avec celles\n"
             "@      du fichier suite, les autres sont initialisees a zero.\n"
             "@\n"
             "@    Si ACTUEL < AMONT, on initialise les %d 1eres\n"
             "@      statistiques supplementaires actuelles avec les 1eres\n"
             "@      du fichier suite, le reste des statistiques du fichier\n"
             "@      suite sont perdues.\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             musbor,
             cs_glob_lagr_boundary_interactions->nusbor,
             cs_glob_lagr_boundary_interactions->nusbor,
             cs_glob_lagr_boundary_interactions->nusbor);

        /*  --> Lecture des stats aux frontieres. Pas de traitement d'erreurs,    */
        /*        on suppose qu'elles sont dues a un changement de physique. */
        itysup  = 3;
        nbval   = 1;

        for (cs_lnum_t ivar = 0; ivar < cs_glob_lagr_dim->nvisbr; ivar++) {

          char rubriq[32];
          sprintf(rubriq, "stat_bord_%s",
                  cs_glob_lagr_boundary_interactions->nombrd[ivar]);
          ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup,
                                           nbval, CS_TYPE_cs_real_t,
                                           &bound_stat[nfabor * ivar]);

        }

      }

    }

    if (cs_glob_lagr_time_scheme->iilagr == 2) {

      itysup  = 0;
      nbval   = 1;
      {
        cs_int_t *tabvar;

        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[] = "iteration_debut_termes_sources_stationnaires";
        ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                         CS_TYPE_cs_int_t, tabvar);
        mstits  = tabvar[0];
        BFT_FREE(tabvar);
      }
      if (ierror != 0)
        cs_log_printf
          (CS_LOG_DEFAULT,
           "@\n"
           "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
           "@\n"
           "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
           "@    =========     LAGRANGIEN %s\n"
           "@      TYPE DE FICHIER INCORRECT\n"
           "@\n"
           "@    Ce fichier ne semble pas etre un fichier\n"
           "@      suite Lagrangien.\n"
           "@\n"
           "@    Le calcul ne peut etre execute.\n"
           "@\n"
           "@    Verifier que le fichier suite utilise correspond bien\n"
           "@        a un fichier suite Lagrangien.\n"
           "@\n"
           "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
           "@\n",
           ficsui);
      // "iteration_debut_termes_sources_stationnaires", "NSTITS", mstits);


      /*  ---> S'il y a une erreur, on suppose que c'est parce que le fichier   */
      /*         suite ne contient pas d'infos sur les TS de couplage retour.   */
      /*         Dans ce cas, si on est en instationnaire on se dit que c'est   */
      /*         pas grave, on saute l'etape et on continue. Par contre si on   */
      /*         est dans une configuration de calcul de stats aux frontieres   */
      /*         en stationnaire on stoppe.  */

      if (ierror != 0) {

        if (   cs_glob_lagr_time_scheme->isttio == 0
            || (   cs_glob_lagr_time_scheme->isttio == 1
                && cs_glob_time_step->nt_cur < cs_glob_lagr_source_terms->nstits))
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      LES TERMES SOURCES DE COUPLAGE RETOUR DU CALCUL AMONT\n"
             "@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTS\n"
             "@        DU FICHIER SUITE\n"
             "@\n"
             "@    Les indicateurs concernant le calcul\n"
             "@      instationnaire/stationnaire des termes sources\n"
             "@      de couplage retour sont positionnes en mode\n"
             "@      instationnaire ou en debut de calcul stationnaire :\n"
             "@\n"
             "@          ISTTIO    NSTITS    Iter de redemarrage\n"
             "@      %10d%10d%10d\n"
             "@\n"
             "@    Ils seront initialisees par des valeurs par defaut.\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             cs_glob_lagr_time_scheme->isttio,
             cs_glob_lagr_source_terms->nstits,
             cs_glob_time_step->nt_cur+1);

        else
          bft_error
            (__FILE__, __LINE__, 0,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      LES TERMES SOURCES DE COUPLAGE RETOUR DU CALCUL AMONT\n"
             "@        NE PEUVENT PAS ETRE RELUES OU SONT ABSENTS\n"
             "@        DU FICHIER SUITE\n"
             "@\n"
             "@    Les indicateurs concernant le calcul\n"
             "@      instationnaire/stationnaire des termes sources\n"
             "@      de couplage retour sont positionnes\n"
             "@      en mode stationnaire :\n"
             "@\n"
             "@          ISTTIO    NSTITS    Iter de redemarrage\n"
             "@      %10d%10d%10d\n"
             "@\n"
             "@    Le calcul ne peut pas etre execute.\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             cs_glob_lagr_time_scheme->isttio,
             cs_glob_lagr_source_terms->nstits,
             cs_glob_time_step->nt_cur+1);

      }
      /* --> A partir d'ici on considere que le fichier suite contient     */
      /*       des stats volumiques     */
      else {

        {
          cs_int_t *tabvar;

          BFT_MALLOC(tabvar, nbval, cs_int_t);

          char rubriq[] = "modele_turbulence_termes_sources";
          ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                           CS_TYPE_cs_int_t, tabvar);
          jturb   = tabvar[0];
          jtytur  = jturb / 10;
          BFT_FREE(tabvar);
        }

        if (   jsttio != cs_glob_lagr_time_scheme->isttio
            || mstits != cs_glob_lagr_source_terms->nstits) {

          char car8[9], kar8[9];
          if (jtytur == 2)
            sprintf(car8, "k-eps");
          if (jtytur == 3)
            sprintf(car8, "Rij-eps");
          if (jturb == 50)
            sprintf(car8, "v2f");
          if (jturb == 60)
            sprintf(car8, "k-omega");
          if (extra->itytur == 2)
            sprintf(kar8, "k-eps");
          if (extra->itytur == 3)
            sprintf(kar8, "Rij-eps");
          if (extra->iturb == 50)
            sprintf(kar8, "v2f");
          if (extra->iturb == 60)
            sprintf(kar8, "k-omega");

          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
             "@\n"
             "@    Les indicateurs concernant le calcul\n"
             "@      instationnaire/stationnaire des termes sources\n"
             "@      de couplage retour sont modifies :\n"
             "@\n"
             "@                   ISTTIO    NSTITS    Turbulence\n"
             "@       AMONT : %10d%10d%13s\n"
             "@       ACTUEL: %10d%10d%13s\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             jsttio,
             mstits,
             car8,
             cs_glob_lagr_time_scheme->isttio,
             cs_glob_lagr_source_terms->nstits,
             kar8);

        }

        /*  --> Lecture de l'avancement du couplage retour    */
        {
          cs_int_t *tabvar;

          BFT_MALLOC(tabvar, nbval, cs_int_t);

          char rubriq[] = "nombre_iterations_termes_sources_stationnaires";
          ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                           CS_TYPE_cs_int_t, tabvar);
          cs_glob_lagr_source_terms->npts = tabvar[0];
          BFT_FREE(tabvar);
        }
        if (ierror != 0)
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@\n"
             "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
             "@ %s\n"
             "@\n"
             "@    Le mot cle est initialise avec sa valeur par defaut\n"
             "@      ou celle donnee dans le sous-programme cs_user_lagr_model :\n"
             "@        %s  = %d\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             "nombre_iterations_termes_sources_stationnaires",
             "NPTS",
             cs_glob_lagr_source_terms->npts);

        /*  --> Verif de coherence de l'avancement du calcul avec les   */
        /*       indicateurs de calcul de la suite actuelle : */

        if (   cs_glob_lagr_source_terms->npts == 0
               && (   cs_glob_lagr_time_scheme->isttio == 1
                   && cs_glob_lagr_source_terms->nstits <= cs_glob_time_step->nt_cur))
          bft_error
            (__FILE__, __LINE__, 0,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
             "@\n"
             "@    LE CALCUL ACTUEL DES TERMES SOURCES DE COUPLAGE RETOUR\n"
             "@      EST EN MODE STATIONNAIRE, ALORS QUE LE FICHIER\n"
             "@      SUITE CONTIENT DES TERMES SOURCES INSTATIONNAIRES.\n"
             "@\n"
             "@    NSTITS devrait etre un entier superieur ou egal\n"
             "@      a l'iteration Lagrangienne absolue de redemarrage\n"
             "@      du calcul (iteration : %d)\n"
             "@\n"
             "@      Il vaut ici NSTITS = %d\n"
             "@\n"
             "@    Le calcul ne peut pas etre execute.\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui,
             cs_glob_time_step->nt_cur+1,
             cs_glob_lagr_source_terms->nstits);

        if (cs_glob_lagr_source_terms->npts > 0
            && (   (   cs_glob_lagr_time_scheme->isttio == 1
                    && cs_glob_time_step->nt_cur <= cs_glob_lagr_source_terms->nstits)
                || cs_glob_lagr_time_scheme->isttio == 0))
          cs_log_printf
            (CS_LOG_DEFAULT,
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n"
             "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
             "@    =========     LAGRANGIEN %s\n"
             "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
             "@\n"
             "@    LE CALCUL ACTUEL DES TERMES SOURCES DE COUPLAGE RETOUR\n"
             "@      EST EN MODE INSTATIONNAIRE, ALORS QUE LE FICHIER\n"
             "@      SUITE CONTIENT DES TERMES SOURCES STATIONNAIRES.\n"
             "@\n"
             "@    Les termes sources de couplage retour stationnaires\n"
             "@      amont seront remises a zero.\n"
             "@\n"
             "@    Le calcul se poursuit...\n"
             "@\n"
             "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
             "@\n",
             ficsui);

        if (cs_glob_lagr_source_terms->npts > 0
            && (   cs_glob_lagr_time_scheme->isttio == 1
                && cs_glob_time_step->nt_cur >= cs_glob_lagr_source_terms->nstits)) {

          if (mstits != cs_glob_lagr_source_terms->nstits)
            bft_error
              (__FILE__, __LINE__, 0,
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n"
               "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
               "@    =========     LAGRANGIEN %s\n"
               "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
               "@\n"
               "@    LE CALCUL SE POURSUIT AVEC UN CALCUL DES TERMES\n"
               "@      SOURCES DE COUPLAGE RETOUR EN MODE STATIONNAIRE\n"
               "@      MAIS LES INDICATEURS DE CONTROLES DES TERMES SOURCES\n"
               "@      ON ETE MODIFIEES.\n"
               "@\n"
               "@    Pour eviter les incoherences dans le calcul\n"
               "@      NSTITS ne devrait pas etre modifie entre deux calculs\n"
               "@      de termes sources de couplage retour stationnaires.\n"
               "@\n"
               "@    Le calcul ne peut pas etre execute.\n"
               "@\n"
               "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
               "@\n",
               ficsui);

        }

        /* On donne des labels au different TS pour les noms de rubriques
         * On donne le meme label au keps, au v2f et au k-omega (meme variable k)*/
        if (cs_glob_lagr_source_terms->ltsdyn == 1) {

          sprintf(nomtsl[cs_glob_lagr_source_terms->itsli],
                  "terme_source_vitesse_implicite");

          if (extra->itytur == 2 || extra->iturb == 50 || extra->iturb == 60)
            sprintf(nomtsl[cs_glob_lagr_source_terms->itske],
                    "terme_source_turbulence_keps");

        }

        if (cs_glob_lagr_source_terms->ltsmas == 1)
          sprintf(nomtsl[cs_glob_lagr_source_terms->itsmas], "terme_source_masse");

        if (cs_glob_lagr_source_terms->ltsthe == 1) {

          if (   cs_glob_lagr_model->physical_model == 1
              && cs_glob_lagr_specific_physics->itpvar == 1) {

            sprintf(nomtsl[cs_glob_lagr_source_terms->itste],
                    "terme_source_thermique_explicite");
            sprintf(nomtsl[cs_glob_lagr_source_terms->itsti],
                    "terme_source_thermique_implicite");

          }
          else if (cs_glob_lagr_model->physical_model == 2) {

            sprintf(nomtsl[cs_glob_lagr_source_terms->itste],
                    "terme_source_thermique_explicite");
            sprintf(nomtsl[cs_glob_lagr_source_terms->itsti],
                    "terme_source_thermique_implicite");

            for (int icha = 0; icha < extra->ncharb; icha++) {

              sprintf(nomtsl[cs_glob_lagr_source_terms->itsmv1[icha]],
                      "terme_source_legeres_F1_%04d", icha);
              sprintf(nomtsl[cs_glob_lagr_source_terms->itsmv2[icha]],
                      "terme_source_lourdes_F2_%04d", icha);

            }

            sprintf(nomtsl[cs_glob_lagr_source_terms->itsco],
                    "terme_source_F3");
            sprintf(nomtsl[cs_glob_lagr_source_terms->itsfp4],
                    "terme_source_variance_traceur_air");

          }

        }

        /* Old style return coupling terms */

        itysup  = 1;
        nbval   = 1;
        for (cs_lnum_t ivar = 0; ivar < cs_glob_lagr_dim->ntersl; ivar++) {

          cs_real_t *st_val = cs_glob_lagr_source_terms->st_val + ivar*n_cells_ext;

          cs_restart_read_section(cs_lag_stat_restart,
                                  nomtsl[ivar+1], itysup, nbval,
                                  CS_TYPE_cs_real_t, st_val);
        }

        /* New style return coupling terms */

        {
          cs_field_t *f = cs_field_by_name_try("velocity_st_lagr");
          if (f != NULL) {
            cs_restart_read_real_3_t_compat(cs_lag_stat_restart,
                                            f->name,
                                            "terme_source_vitesseX",
                                            "terme_source_vitesseY",
                                            "terme_source_vitesseZ",
                                            f->location_id,
                                            (cs_real_3_t *)(f->vals));
          }
        }

        {
          cs_field_t *f = cs_field_by_name_try("rij_st_lagr");
          if (f != NULL) {
            cs_restart_read_real_6_t_compat(cs_lag_stat_restart,
                                            f->name,
                                            "terme_source_turbulence_R11",
                                            "terme_source_turbulence_R22",
                                            "terme_source_turbulence_R33",
                                            "terme_source_turbulence_R12",
                                            "terme_source_turbulence_R23",
                                            "terme_source_turbulence_R13",
                                            f->location_id,
                                            (cs_real_6_t *)(f->vals));
          }
        }

      }

    }

    cs_restart_read_fields(cs_lag_stat_restart, CS_RESTART_LAGR_STAT);

    cs_log_printf(CS_LOG_DEFAULT,
                  _("      Fin   de la lecture\n"));

    /*  ---> Fermeture du fichier suite    */
    cs_restart_destroy(&cs_lag_stat_restart);

    /* ---> En cas d'erreur, on continue quand meme  */
    cs_log_printf(CS_LOG_DEFAULT,
                  _("    Fin de la lecture du fichier suite\n"
                    "      sur les statistiques et TS couplage retour\n"));

  }

  cs_log_printf(CS_LOG_DEFAULT,
                "\n-------------------------------------------------------------\n");

  BFT_FREE(nomtsl);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Read Lagrangian particle and statistics restart files.
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_restart_read_p(void)
{
  /* counters */
  cs_lnum_t ierror, itysup, nbval, iok;
  iok = 0;
  /* read variables */
  cs_lnum_t mvls, jphyla, jtpvar, jdpvar, jmpvar;

  cs_lagr_particle_counter_t *pc = cs_lagr_get_particle_counter();

  if (cs_glob_lagr_time_scheme->isuila == 0)
    return;

  /* ====================================================================   */
  /* 1. LECTURE DU FICHIER SUITE : VARIABLES LIEES AUX PARTICULES */
  /* ====================================================================   */

  /*  ---> Ouverture */
  cs_log_printf(CS_LOG_DEFAULT,
                _("   ** INFORMATIONS SUR LE CALCUL LAGRANGIEN\n"
                  "      -------------------------------------\n"
                  "    Lecture d'un fichier suite\n"
                  "    sur les variables liees aux particules\n"));

  static cs_restart_t  *cs_lag_stat_restart = NULL;
  char ficsui[] = "lagrangian";

  cs_lag_stat_restart = cs_restart_create(ficsui, NULL, CS_RESTART_MODE_READ);

  cs_log_printf(CS_LOG_DEFAULT,
                "      Debut de la lecture\n");

  /*  ---> Type de fichier suite    */
  /*        Pourrait porter le numero de version si besoin.  */

  itysup = 0;
  nbval  = 1;
  {
    cs_int_t *tabvar;

    BFT_MALLOC(tabvar, nbval, cs_int_t);

    char rubriq[] = "version_fichier_suite_Lagrangien_variables";
    ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                     CS_TYPE_cs_int_t, tabvar);
    BFT_FREE(tabvar);
  }

  if (ierror != 0)
    bft_error
      (__FILE__, __LINE__, 0,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      TYPE DE FICHIER INCORRECT\n"
       "@\n"
       "@    Ce fichier ne semble pas etre un fichier\n"
       "@      suite Lagrangien.\n"
       "@\n"
       "@    Le calcul ne peut etre execute.\n"
       "@\n"
       "@    Verifier que le fichier suite utilise correspond bien\n"
       "@        a un fichier suite Lagrangien.\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui);

  /*     Dimensions des supports    */
  bool ncelok, nfaiok, nfabok, nsomok;
  cs_restart_check_base_location(cs_lag_stat_restart, &ncelok, &nfaiok,
                                 &nfabok, &nsomok);

  if (! ncelok) {

    cs_log_printf
      (CS_LOG_DEFAULT,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
       "@\n"
       "@    Le nombre de cellules a ete modifie\n"
       "@\n"
       "@    Le calcul ne peut etre execute.\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui);
    iok++;

  }

  if (!nfaiok)
    cs_log_printf
      (CS_LOG_DEFAULT,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
       "@\n"
       "@    Le nombre de faces %s a ete modifie.\n"
       "@\n"
       "@    Le calcul peut etre execute mais les donnees\n"
       "@      sur les faces %s ne seront pas relues\n"
       "@      dans le fichier suite.\n"
       "@    Elles seront initialisees par des valeurs par defaut.\n"
       "@\n"
       "@    Le calcul se poursuit...\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui,
       "internes",
       "internes");

  if (!nfabok)
    cs_log_printf
      (CS_LOG_DEFAULT,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
       "@\n"
       "@    Le nombre de faces %s a ete modifie.\n"
       "@\n"
       "@    Le calcul peut etre execute mais les donnees\n"
       "@      sur les faces %s ne seront pas relues\n"
       "@      dans le fichier suite.\n"
       "@    Elles seront initialisees par des valeurs par defaut.\n"
       "@\n"
       "@    Le calcul se poursuit...\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui,
       "de bord",
       "de bord");

  itysup = 0;
  nbval  = 1;

  /*     Physique associee aux particules     */
  {
    cs_int_t *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    char rubriq[] = "indicateur_physique_particules";
    ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                     CS_TYPE_cs_int_t, tabvar);
    jphyla = tabvar[0];
    BFT_FREE(tabvar);
    if (ierror != 0) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n"
         "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
         "@    =========     LAGRANGIEN %s\n"
         "@\n"
         "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
         "@ %s\n"
         "@\n"
         "@    Le calcul ne peut pas etre execute.\n"
         "@\n"
         "@    Verifier que ce fichier suite\n"
         "@        correspond bien a un fichier suite Lagrangien,\n"
         "@        et qu'il n'a pas ete endommage.\n"
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n",
         ficsui,
         rubriq);
      iok++;

    }

  }

  {
    cs_int_t *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    char rubriq[] = "indicateur_temperature_particules";
    ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                     CS_TYPE_cs_int_t, tabvar);
    jtpvar = tabvar[0];
    BFT_FREE(tabvar);

    if (ierror != 0) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n"
         "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
         "@    =========     LAGRANGIEN %s\n"
         "@\n"
         "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
         "@ %s\n"
         "@\n"
         "@    Le calcul ne peut pas etre execute.\n"
         "@\n"
         "@    Verifier que ce fichier suite\n"
         "@        correspond bien a un fichier suite Lagrangien,\n"
         "@        et qu'il n'a pas ete endommage.\n"
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n",
         ficsui,
         rubriq);
      iok++;

    }

  }

  /*     Arret  */
  if (iok != 0)
    cs_exit(1);

  {
    cs_int_t *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    char rubriq[] = "indicateur_diametre_particules";
    ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                     CS_TYPE_cs_int_t, tabvar);
    jdpvar      = tabvar[0];
    BFT_FREE(tabvar);

    if (ierror != 0) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n"
         "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
         "@    =========     LAGRANGIEN %s\n"
         "@\n"
         "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
         "@ %s\n"
         "@\n"
         "@    Le calcul se poursuit...\n"
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n",
         ficsui,
         rubriq);
      jdpvar    = cs_glob_lagr_specific_physics->idpvar;

    }

  }

  {
    cs_int_t *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    char rubriq[] = "indicateur_masse_particules";
    ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                     CS_TYPE_cs_int_t, tabvar);
    jmpvar = tabvar[0];
    BFT_FREE(tabvar);

    if (ierror != 0) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n"
         "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
         "@    =========     LAGRANGIEN %s\n"
         "@\n"
         "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
         "@ %s\n"
         "@\n"
         "@    Le calcul se poursuit...\n"
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n",
         ficsui,
         rubriq);
      jmpvar = cs_glob_lagr_specific_physics->impvar;

    }

  }

  /* ---> On previent si des parametres sont differents */

  if (   jphyla != cs_glob_lagr_model->physical_model
      || jtpvar != cs_glob_lagr_specific_physics->itpvar
      || jdpvar != cs_glob_lagr_specific_physics->idpvar
      || jmpvar != cs_glob_lagr_specific_physics->impvar)

    cs_log_printf
      (CS_LOG_DEFAULT,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
       "@\n"
       "@    Les indicateurs concernant la physique associee\n"
       "@      aux particules sont modifies :\n"
       "@\n"
       "@              IPHYLA    ITPVAR    IDPVAR    IMPVAR\n"
       "@  AMONT : %10d%10d%10d%10d\n"
       "@  ACTUEL: %10d%10d%10d%10d\n"
       "@\n"
       "@    Le calcul se poursuit...\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui,
       jphyla,
       jtpvar,
       jdpvar,
       jmpvar,
       cs_glob_lagr_model->physical_model,
       cs_glob_lagr_specific_physics->itpvar,
       cs_glob_lagr_specific_physics->idpvar,
       cs_glob_lagr_specific_physics->impvar);

  /* ---> Verification de la compatibilite si changement de thermique  */

  if (jphyla != 0 && cs_glob_lagr_model->physical_model == 0)
    cs_log_printf
      (CS_LOG_DEFAULT,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
       "@\n"
       "@    Aucune selection de physique associee aux particules\n"
       "@      n'est active. Les donnees amont sont perdues.\n"
       "@\n"
       "@    Le calcul se poursuit...\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui);

  if (cs_glob_lagr_specific_physics->itpvar == 1 && jtpvar == 0)
    cs_log_printf
      (CS_LOG_DEFAULT,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
       "@\n"
       "@    Une equation sur la temperature des particules est\n"
       "@      enclenchee en cours de calcul.\n"
       "@    Initialisation par defaut :\n"
       "@       Temperature TPART = %14.5E\n"
       "@       Chaleur massique CPPART = %14.5E\n"
       "@\n"
       "@    Le calcul se poursuit...\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui,
       cs_glob_lagr_specific_physics->tpart,
       cs_glob_lagr_specific_physics->cppart);

  if (cs_glob_lagr_model->physical_model == 2 && jphyla != 2)
    bft_error
      (__FILE__, __LINE__, 0,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
       "@\n"
       "@    L'indicateur d'un calcul Lagrangien de grains\n"
       "@      de charbon est enclenche (IPHYLA = 2).\n"
       "@    Ce fichier suite ne correspond pas\n"
       "@      a un calcul Lagrangien de grains de charbon.\n"
       "@\n"
       "@    Le calcul ne peut etre execute.\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui);

  if (   (   jphyla == 2
          && cs_glob_lagr_model->physical_model == 1)
      || (   jphyla == 1
          && cs_glob_lagr_model->physical_model == 2))
    bft_error
      (__FILE__, __LINE__, 0,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
       "@\n"
       "@    Ce fichier suite correspond\n"
       "@      a un calcul Lagrangien de grains de charbon.\n"
       "@    L'indicateur de physique actuel associee aux particules\n"
       "@      a une valeur non permise dans le cadre d'une suite\n"
       "@      d'un calcul Lagrangien de grains de charbon.\n"
       "@\n"
       "@    Le calcul ne peut etre execute.\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui);

  if (   cs_glob_lagr_stat_options->isuist == 0
      && cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->idstnt)
    bft_error
      (__FILE__, __LINE__, 0,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
       "@\n"
       "@    L'INDICATEUR DE CALCUL DES STATISTIQUES VOLUMIQUES\n"
       "@       A UNE VALEUR NON PERMISE (LAGLEC_P).\n"
       "@\n"
       "@    LORSQU'IL N'Y A PAS DE SUITE DE CALCUL SUR LES\n"
       "@    STATISTIQUES VOLUMIQUES (ISUIST = %d)\n"
       "@    IDSTNT DOIT ETRE UN ENTIER SUPERIEUR AU NUMERO\n"
       "@       DE L'ITERATION LAGRANGIENNE DE REDEMARRAGE %d\n"
       "@       IL VAUT ICI IDSTNT = %d\n"
       "@\n"
       "@  Le calcul ne sera pas execute.\n"
       "@\n"
       "@  Verifier la valeur de IDSTNT.\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui,
       cs_glob_lagr_stat_options->isuist,
       cs_glob_time_step->nt_cur +1,
       cs_glob_lagr_stat_options->idstnt);

  if (   cs_glob_lagr_stat_options->isuist == 0
      && cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->nstist)
    bft_error
      (__FILE__, __LINE__, 0,
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n"
       "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
       "@    =========     LAGRANGIEN %s\n"
       "@      DONNEES AMONT ET ACTUELLES INCOHERENTES\n"
       "@\n"
       "@    L'INDICATEUR DE CALCUL STATIONNAIRES DES STATISTIQUES\n"
       "@       AUX FRONTIERES A UNE VALEUR NON PERMISE (LAGLEC_P).\n"
       "@\n"
       "@    LORSQU'IL N'Y A PAS DE SUITE DE CALCUL SUR LES\n"
       "@    STATISTIQUES AUX FRONTIERES (ISUIST = %d),\n"
       "@    NSTBOR DOIT ETRE UN ENTIER SUPERIEUR AU NUMERO\n"
       "@       DE L'ITERATION LAGRANGIENNE DE REDEMARRAGE %d\n"
       "@       IL VAUT ICI NSTBOR = %d\n"
       "@\n"
       "@  Le calcul ne sera pas execute.\n"
       "@\n"
       "@  Verifier la valeur de NSTBOR.\n"
       "@\n"
       "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       "@\n",
       ficsui,
       cs_glob_lagr_stat_options->isuist,
       cs_glob_time_step->nt_cur +1,
       cs_glob_lagr_stat_options->nstist);

  {
    cs_real_t *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_real_t);

    char rubriq[] = "temps_physique_Lagrangien";
    ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                     CS_TYPE_cs_real_t, tabvar);
    cs_glob_lagr_time_step->ttclag = tabvar[0];
    BFT_FREE(tabvar);

    if (ierror != 0)
      cs_log_printf
        (CS_LOG_DEFAULT,
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n"
         "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
         "@    =========     LAGRANGIEN %s\n"
         "@\n"
         "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
         "@ %s\n"
         "@\n"
         "@    Le mot cle est initialise avec sa valeur par defaut\n"
         "@      ou celle donnee dans le sous-programme cs_user_lagr_model :\n"
         "@        %s  = %14.5E\n"
         "@\n"
         "@    Le calcul se poursuit...\n"
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n",
         ficsui,
         "temps_physique_Lagrangien",
         "TTCLAG",
         cs_glob_lagr_time_step->ttclag);
  }

  {
    cs_int_t *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    char rubriq[] = "nombre_total_particules";
    ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                     CS_TYPE_cs_int_t, tabvar);
    pc->n_g_cumulative_total = tabvar[0];
    BFT_FREE(tabvar);

    if (ierror != 0)
      cs_log_printf
        (CS_LOG_DEFAULT,
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n"
         "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
         "@    =========     LAGRANGIEN %s\n"
         "@\n"
         "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
         "@ %s\n"
         "@\n"
         "@    Le mot cle est initialise avec sa valeur par defaut\n"
         "@      ou celle donnee dans le sous-programme cs_user_lagr_model :\n"
         "@        %s  = %llu\n"
         "@\n"
         "@    Le calcul se poursuit...\n"
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n",
         ficsui,
         "nombre_total_particules",
         "NBPTOT",
         (unsigned long long)(pc->n_g_cumulative_total));
  }

  {
    cs_int_t *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    char rubriq[] = "nombre_particules_perdues";
    ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                     CS_TYPE_cs_int_t, tabvar);
    pc->n_g_cumulative_failed = tabvar[0];
    BFT_FREE(tabvar);

    if (ierror != 0)
      cs_log_printf
        (CS_LOG_DEFAULT,
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n"
         "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
         "@    =========     LAGRANGIEN %s\n"
         "@\n"
         "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
         "@ %s\n"
         "@\n"
         "@    Le mot cle est initialise avec sa valeur par defaut\n"
         "@      ou celle donnee dans le sous-programme cs_user_lagr_model :\n"
         "@        %s  = %llu\n"
         "@\n"
         "@    Le calcul se poursuit...\n"
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n",
         ficsui,
         "nombre_particules_perdues",
         "NBPERT",
         (unsigned long long)(pc->n_g_cumulative_failed));
  }

  {
    cs_int_t *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);
    char rubriq[] = "nombre_variables_utilisateur";
    ierror = cs_restart_read_section(cs_lag_stat_restart, rubriq, itysup, nbval,
                                     CS_TYPE_cs_int_t, tabvar);
    mvls   = tabvar[0];
    BFT_FREE(tabvar);

    if (ierror != 0) {

      mvls = 0;

      if (cs_glob_lagr_model->n_user_variables > 0)
        cs_log_printf
          (CS_LOG_DEFAULT,
           "@\n"
           "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
           "@\n"
           "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
           "@    =========     LAGRANGIEN %s\n"
           "@\n"
           "@      ERREUR A LA LECTURE DE LA RUBRIQUE\n"
           "@ %s\n"
           "@\n"
           "@    Le calcul se poursuit...\n"
           "@\n"
           "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
           "@\n",
           ficsui, "nombre_variables_utilisateur");
    }

    if (cs_glob_lagr_model->n_user_variables < mvls) {

      cs_log_printf
        (CS_LOG_DEFAULT,
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n"
         "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
         "@    =========     LAGRANGIEN %s'\n"
         "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
         "@\n"
         "@    L'indicateur du  nombre de variables supplementaires\n"
         "@      utilisateur est modifie, ou n'a pas pu etre relu.\n"
         "@\n"
         "@              NVLS\n"
         "@    AMONT : %d      ACTUEL : %d\n"
         "@\n"
         "@    Si ACTUEL > AMONT, on initialise les %d 1eres\n"
         "@      variables supplementaires actuelles avec celles\n"
         "@      du fichier suite, les autres sont initialisees a zero.\n"
         "@\n"
         "@    Si ACTUEL < AMONT, on initialise les %d 1eres\n"
         "@      variables supplementaires actuelles avec les 1eres\n"
         "@      du fichier suite, le reste des variables du fichier\n"
         "@      suite sont perdues.\n"
         "@\n"
         "@    Le calcul se poursuit...\n"
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n",
         ficsui,
         mvls,
         cs_glob_lagr_model->n_user_variables,
         cs_glob_lagr_model->n_user_variables,
         cs_glob_lagr_model->n_user_variables);
      mvls = cs_glob_lagr_model->n_user_variables;

    }
    else if (cs_glob_lagr_model->n_user_variables > mvls)
      cs_log_printf
        (CS_LOG_DEFAULT,
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n"
         "@ @@ ATTENTION : A LA LECTURE DU FICHIER SUITE\n"
         "@    =========     LAGRANGIEN %s'\n"
         "@      DONNEES AMONT ET ACTUELLES DIFFERENTES\n"
         "@\n"
         "@    L'indicateur du  nombre de variables supplementaires\n"
         "@      utilisateur est modifie, ou n'a pas pu etre relu.\n"
         "@\n"
         "@              NVLS\n"
         "@    AMONT : %d      ACTUEL : %d\n"
         "@\n"
         "@    Si ACTUEL > AMONT, on initialise les %d 1eres\n"
         "@      variables supplementaires actuelles avec celles\n"
         "@      du fichier suite, les autres sont initialisees a zero.\n"
         "@\n"
         "@    Si ACTUEL < AMONT, on initialise les %d 1eres\n"
         "@      variables supplementaires actuelles avec les 1eres\n"
         "@      du fichier suite, le reste des variables du fichier\n"
         "@      suite sont perdues.\n"
         "@\n"
         "@    Le calcul se poursuit...\n"
         "@\n"
         "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
         "@\n",
         ficsui,
         mvls,
         cs_glob_lagr_model->n_user_variables,
         cs_glob_lagr_model->n_user_variables,
         cs_glob_lagr_model->n_user_variables);
  }

  /* Particle data */
  cs_lagr_restart_read_particle_data(cs_lag_stat_restart);

  cs_restart_read_fields(cs_lag_stat_restart, CS_RESTART_LAGR);

  cs_log_printf(CS_LOG_DEFAULT,
                _("    End reading particle data restart file\n"));

  /* Close restart file */
  cs_restart_destroy(&cs_lag_stat_restart);

  cs_log_printf(CS_LOG_DEFAULT,
                _("    End reading particle statistics restart file\n"));
}

/*--------------------------------------------------------------------*/
/*! \brief Output Lagrangian restart files
 */
/*--------------------------------------------------------------------*/

void
cs_restart_lagrangian_checkpoint_write(void)
{
  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;

  int nvplmx = 50+4*cs_glob_lagr_const_dim->nlayer;

  typedef char cs_char_64_t[64];

  cs_char_64_t *nomtsl = NULL;

  BFT_MALLOC(nomtsl, nvplmx, cs_char_64_t);

  /* Output restart file: variables related to particles */
  /*-----------------------------------------------------*/

  static cs_restart_t  *cs_lag_stat_restart = NULL;
  /* Open restart file    */
  char const *ficsui = "lagrangian";
  cs_lag_stat_restart = cs_restart_create(ficsui, NULL, CS_RESTART_MODE_WRITE);

  cs_log_printf(CS_LOG_DEFAULT,
                _("   ** Writing the Lagrangian restart file\n"
                  "-----------------------------------\n"));

  /* Entete et Infos sur le calcul ou on saute si erreur
   *   On inclut une rubrique destinee a distinguer ce fichier
   *     d'un autre fichier suite
   *   Pour le moment, IVERS n'est pas utilise */

  cs_lnum_t itysup  = 0;
  cs_lnum_t nbval  = 1;

  {
    int ivers = 32000;

    cs_restart_write_section(cs_lag_stat_restart,
                             "version_fichier_suite_Lagrangien_variables",
                             itysup, 1, CS_TYPE_cs_int_t, &ivers);
  }

  {
    cs_real_t val = cs_glob_lagr_time_step->ttclag;

    cs_restart_write_section(cs_lag_stat_restart,
                             "temps_physique_Lagrangien",
                             itysup, 1, CS_TYPE_cs_real_t, &val);
  }

  /* Infos sur le suivi du calcul   */
  {
    cs_int_t   *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    tabvar[0] = cs_glob_lagr_particle_counter->n_g_cumulative_total;
    char rubriq[] = "nombre_total_particules";

    cs_restart_write_section(cs_lag_stat_restart, rubriq,
                             itysup, nbval, CS_TYPE_cs_int_t, tabvar);
    BFT_FREE(tabvar);
  }

  {
    cs_int_t   *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    tabvar[0] = cs_glob_lagr_particle_counter->n_g_cumulative_failed;
    char rubriq[] = "nombre_particules_perdues";

    cs_restart_write_section(cs_lag_stat_restart, rubriq,
                             itysup, nbval, CS_TYPE_cs_int_t, tabvar);
    BFT_FREE(tabvar);
  }

  {
    cs_int_t   *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    tabvar[0] = cs_glob_lagr_model->physical_model;
    char rubriq[] = "indicateur_physique_particules";

    cs_restart_write_section(cs_lag_stat_restart, rubriq,
                             itysup, nbval, CS_TYPE_cs_int_t, tabvar);
    BFT_FREE(tabvar);
  }

  {
    cs_int_t   *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    tabvar[0] = cs_glob_lagr_specific_physics->itpvar;
    char rubriq[] = "indicateur_temperature_particules";

    cs_restart_write_section(cs_lag_stat_restart, rubriq,
                             itysup, nbval, CS_TYPE_cs_int_t, tabvar);
    BFT_FREE(tabvar);
  }

  {
    cs_int_t   *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    tabvar[0] = cs_glob_lagr_specific_physics->idpvar;
    char rubriq[]= "indicateur_diametre_particules";

    cs_restart_write_section(cs_lag_stat_restart, rubriq,
                             itysup, nbval, CS_TYPE_cs_int_t, tabvar);
    BFT_FREE(tabvar);
  }

  {
    cs_int_t   *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    tabvar[0] = cs_glob_lagr_specific_physics->impvar;
    char rubriq[] = "indicateur_masse_particules";

    cs_restart_write_section(cs_lag_stat_restart, rubriq,
                             itysup, nbval, CS_TYPE_cs_int_t, tabvar);
    BFT_FREE(tabvar);
  }

  {
    cs_int_t   *tabvar;
    BFT_MALLOC(tabvar, nbval, cs_int_t);

    tabvar[0] = cs_glob_lagr_model->n_user_variables;
    char rubriq[] = "nombre_variables_utilisateur";

    cs_restart_write_section(cs_lag_stat_restart, rubriq,
                             itysup, nbval, CS_TYPE_cs_int_t, tabvar);
    BFT_FREE(tabvar);
  }

  cs_restart_write_fields(cs_lag_stat_restart, CS_RESTART_LAGR);

  cs_log_printf(CS_LOG_DEFAULT,
                _("      End writing info on calculation\n"));

  cs_lagr_restart_write_particle_data(cs_lag_stat_restart);

  cs_log_printf(CS_LOG_DEFAULT,
                _("      End writing of specific info\n"));

  /* ---> Fermeture du fichier suite     */
  cs_restart_destroy(&cs_lag_stat_restart);

  cs_log_printf(CS_LOG_DEFAULT,
                _("    End writing of restart file\n"
                  "      on particle-based variables\n"));

  /*====================================================================
   * 2. ECRITURE DU FICHIER SUITE STATISTIQUES ET TERMES SOURCES
   *    DE COUPLAGE RETOUR
   *====================================================================*/

  if (   cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->idstnt
      || cs_glob_lagr_time_scheme->iilagr == 2
      || cs_glob_lagr_dim->nvisbr > 0) {

    /* ---> Ouverture (et on saute si erreur)   */
    cs_log_printf(CS_LOG_DEFAULT,
                  _("   ** INFORMATION ON LAGRANGIAN CALCULATION\n"
                    "-------------------------------------\n"
                    "    Writing a restart file for volume\n"
                    "    and boundary statistics and for\n"
                    "    return coupling source terms\n"));

    /* Open restart file    */
    char const *ficsuist = "lagrangian_stats";
    cs_lag_stat_restart = cs_restart_create(ficsuist, NULL, CS_RESTART_MODE_WRITE);

    cs_log_printf(CS_LOG_DEFAULT,
                  _("      Start writing statistics and ST\n"));

    /* Entete et Infos sur le calcul ou on saute si erreur
     *   On inclut une rubrique destinee a distinguer ce fichier
     *     d'un autre fichier suite
     *   Pour le moment, IVERS n'est pas utilise*/

    itysup = 0;
    nbval = 1;
    {
      //cs_lnum_t ivers = 111;
      cs_lnum_t ivers = 112;
      cs_int_t  *tabvar;
      BFT_MALLOC(tabvar, nbval, cs_int_t);

      char rubriq[] = "version_fichier_suite_Lagrangien_statistiques";
      tabvar[0] = ivers;

      cs_restart_write_section(cs_lag_stat_restart, rubriq,
                               itysup, nbval, CS_TYPE_cs_int_t, tabvar);
      BFT_FREE(tabvar);
    }


    /* ---> On ecrit ISTTIO c'est utile dans tous les cas */
    {
      cs_int_t *tabvar;
      BFT_MALLOC(tabvar, nbval, cs_int_t);
      char rubriq[] = "indicateur_ecoulement_stationnaire";
      tabvar[0] = cs_glob_lagr_time_scheme->isttio;

      cs_restart_write_section(cs_lag_stat_restart, rubriq,
                               itysup, nbval, CS_TYPE_cs_int_t, tabvar);
      BFT_FREE(tabvar);
    }

    /* --> En premier, on ecrit les statistiques volumiques    */
    if (cs_glob_time_step->nt_cur >= cs_glob_lagr_stat_options->idstnt) {

      {
        cs_int_t *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);
        char rubriq[] = "iteration_debut_statistiques";
        tabvar[0] = cs_glob_lagr_stat_options->idstnt;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      {
        cs_int_t *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);
        char rubriq[] = "iteration_debut_statistiques_stationnaires";
        tabvar[0] = cs_glob_lagr_stat_options->nstist;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      {
        cs_int_t *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);
        char rubriq[] = "classe_statistique_particules";
        tabvar[0] = cs_glob_lagr_model->n_stat_classes;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      /*  Statistiques volumiques  */
      cs_lagr_stat_restart_write(cs_lag_stat_restart);

    }

    /* --> En second, c'est le tour des statistiques aux frontieres */
    if (cs_glob_lagr_dim->nvisbr > 0) {

      itysup = 0;
      nbval = 1;

      {
        cs_int_t   *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[] = "iteration_debut_stats_frontieres_stationnaires";
        tabvar[0] = cs_glob_lagr_stat_options->nstist;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      {
        cs_int_t *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[] = "nombre_iterations_stats_frontieres";
        tabvar[0] = cs_glob_lagr_boundary_interactions->npstft;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      {
        cs_int_t *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[] = "nombre_iterations_stats_frontieres_stationnaires";
        tabvar[0] = cs_glob_lagr_boundary_interactions->npstf;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      {
        cs_real_t *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_real_t);

        char rubriq[] = "temps_stats_frontieres_stationnaires";
        tabvar[0] = cs_glob_lagr_boundary_interactions->tstatp;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_real_t, tabvar);
        BFT_FREE(tabvar);
      }

      {
        cs_int_t *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[] = "nombre_stats_frontieres_utilisateur";
        tabvar[0] = cs_glob_lagr_boundary_interactions->nusbor;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      /*  Statistiques aux frontieres   */
      itysup = 3;
      nbval = 1;

      for (cs_lnum_t ii = 0; ii < cs_glob_lagr_dim->nvisbr; ii++) {

        char rubriq[32];
        sprintf(rubriq, "stat_bord_%s", cs_glob_lagr_boundary_interactions->nombrd[ii]);

        cs_lnum_t nfabor = cs_glob_mesh->n_b_faces;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_real_t,
                                 (void *)(&bound_stat[ii * nfabor]));

      }

    }


    /* --> Enfin, en cas de couplage retour, on ecrit les termes sources */

    if (cs_glob_lagr_time_scheme->iilagr == 2) {

      itysup = 0;
      nbval = 1;

      {
        cs_int_t   *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[] = "iteration_debut_termes_sources_stationnaires";
        tabvar[0] = cs_glob_lagr_source_terms->nstits;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      {
        cs_int_t *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[] = "nombre_iterations_termes_sources_stationnaires";
        tabvar[0] = cs_glob_lagr_source_terms->npts;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      {
        cs_int_t *tabvar;
        BFT_MALLOC(tabvar, nbval, cs_int_t);

        char rubriq[] = "modele_turbulence_termes_sources";
        tabvar[0] = extra->iturb;

        cs_restart_write_section(cs_lag_stat_restart, rubriq,
                                 itysup, nbval, CS_TYPE_cs_int_t, tabvar);
        BFT_FREE(tabvar);
      }

      /* On donne des labels au different TS pour les noms de rubriques    */
      /* On donne le meme label au keps, au v2f et au k-omega (meme variable k) */
      if (cs_glob_lagr_source_terms->ltsdyn == 1) {
        sprintf(nomtsl[cs_glob_lagr_source_terms->itsli], "terme_source_vitesse_implicite");
        if (extra->itytur == 2 || extra->iturb == 50 || extra->iturb == 60) {
          sprintf(nomtsl[cs_glob_lagr_source_terms->itske], "terme_source_turbulence_keps");
        }
      }
      if (cs_glob_lagr_source_terms->ltsmas == 1) {
        sprintf(nomtsl[cs_glob_lagr_source_terms->itsmas], "terme_source_masse");
      }
      if (cs_glob_lagr_source_terms->ltsthe == 1) {
        if (   cs_glob_lagr_model->physical_model == 1
            && cs_glob_lagr_specific_physics->itpvar == 1) {
          sprintf(nomtsl[cs_glob_lagr_source_terms->itste],
                  "terme_source_thermique_explicite");
          sprintf(nomtsl[cs_glob_lagr_source_terms->itsti],
                  "terme_source_thermique_implicite");
        }
        else if (cs_glob_lagr_model->physical_model == 2) {
          sprintf(nomtsl[cs_glob_lagr_source_terms->itste],
                  "terme_source_thermique_explicite");
          sprintf(nomtsl[cs_glob_lagr_source_terms->itsti],
                  "terme_source_thermique_implicite");
          for (int icha = 0; icha < extra->ncharb; icha++) {
            sprintf(nomtsl[cs_glob_lagr_source_terms->itsmv1[icha]],
                    "terme_source_legeres_F1_%04d", icha);
            sprintf(nomtsl[cs_glob_lagr_source_terms->itsmv2[icha]],
                    "terme_source_lourdes_F2_%04d", icha);
          }
          sprintf(nomtsl[cs_glob_lagr_source_terms->itsco],
                  "terme_source_F3");
          sprintf(nomtsl[cs_glob_lagr_source_terms->itsfp4],
                  "terme_source_variance_traceur_air");
        }
      }

      /* Old style return source terms */

      itysup = 1;
      nbval = 1;

      const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;

      for (int ii = 0; ii < cs_glob_lagr_dim->ntersl; ii++) {

        cs_real_t *st_val = cs_glob_lagr_source_terms->st_val + ii*n_cells_ext;

        cs_restart_write_section(cs_lag_stat_restart, nomtsl[ii+1],
                                 itysup, nbval, CS_TYPE_cs_real_t,
                                 st_val);

      }

      /* New style return source terms */

      const char *st_names[] = {"velocity_st_lagr",
                                "rij_st_lagr"};

      for (int st_id = 0; st_id < 2; st_id++) {
        cs_field_t *f = cs_field_by_name_try(st_names[st_id]);
        if (f != NULL)
          cs_restart_write_field_vals(cs_lag_stat_restart, f->id, 0);
      }

    }

    cs_restart_write_fields(cs_lag_stat_restart, CS_RESTART_LAGR_STAT);

    cs_log_printf(CS_LOG_DEFAULT,
                  _("      End writing statistics and ST\n"));

    cs_restart_destroy(&cs_lag_stat_restart);

    cs_log_printf(CS_LOG_DEFAULT,
                  _("    End writing of restart file\n"
                    "      on statistics and return coupling ST\n"));

  }

  BFT_FREE(nomtsl);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
