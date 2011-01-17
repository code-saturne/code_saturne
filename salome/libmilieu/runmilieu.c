/*============================================================================
 *
 *     This file is part of the Code_Saturne CFD tool.
 *
 *     Copyright (C) 2008-2011 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne CFD tool is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne CFD tool is distributed in the hope that it will be
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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "utilitaires.h"
#include "donnees.h"
#include "communication.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "runmilieu.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 *  Global variables
 *============================================================================*/

int     nb_dyn = 0;
int     nb_for = 0;
int     ntcast = 0;
double  lref = 0.0;

double  *xast = NULL;
double  *xvast = NULL;
double  *xvasa = NULL;
double  *xastp = NULL;

double  *foras = NULL;
double  *foaas = NULL;
double  *fopas = NULL;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction runmilieu
 *----------------------------------------------------------------------------*/

void
runmilieu(void *icompo)
{
  /* variables locales */

  int i,j;
  int ierr = 0;
  int icv = 0;
  double c1, c2, c3;
  double alpha, beta;
  double dt_ast, dt_sat;
  double dtold = 0.;
  double dt = dtref;

  /* Les données d'entrée du calcul couplé sont définies dans SALOME dans le
     XML du cas d'étude */

  /* Initialisation de la communication */
  if ((ierr = inicom(icompo)) >= 0) {
    printf(" Initialisation de la communication\n");
  }

  /* Envoi des paramètres aux codes */
  if ((ierr = send_param(icompo)) >= 0) {
    printf(" Envoi des parametres de calcul aux codes\n");
  }

  /* Dimensionnement des tableaux et initialisation */

  /* Réception des données géométriques (nb_for et nb_dyn) */
  if (ierr >= 0) {
    ierr = recv_geom(icompo);

    printf("----------------------------------\n");
    printf(" Parametres geometriques\n");
    printf(" nombre de faces couplees : %i\n", nb_for);
    printf(" nombre de noeuds couples : %i\n", nb_dyn);
    printf(" longueur de reference (m): %4.2le\n", lref  );
    printf("----------------------------------\n");

    /* dynamique */
    alldyn();

    /* efforts  */
    allfor();

    /* Coefficients de prediction */
    c1    = 0.;
    c2    = 0.;
    c3    = 0.;
    beta  = 0.;
    alpha = 0.;
  }

  /* Envoi des donnees geometriques a Code_Aster (a supprimer a l avenir)*/
  if (ierr >= 0) {
    ierr = send_geom(icompo);
  }

  /* initialisation pas de temps */
  dt = 0.;
  dt_ast = 0.;
  dt_sat = 0.;

  /* initialisation iteration de couplage */
  ntcast = 0;

 /* Boucle principale */
  i = 1;
  while (ierr >= 0) {
    printf("\n");
    printf("----------------------------------------------------\n");
    printf("\n");
    printf("*********************************\n");
    printf("*         iteration %i          *\n", i);
    printf("*********************************\n");
    printf("\n");
    printf("\n");

    /* Info sur les schémas en temps */
    if (nbssit <= 1) {
      printf("Schema explicite de marche en temps\n");
    }
    else {
      printf("Schema implicite de marche en temps\n");
      printf("Nombre de sous-iterations: %i\n", nbssit);
    }

    /* Gestion pas de temps */

    /* Réception pdt Code_Aster et pdt Code_Saturne */
    ierr = recv_pdt(icompo,&(dt_ast), &(dt_sat), i);

    printf("----------------------------------\n");
    printf("pas de temps de reference: %4.21e \n", dtref );
    printf("pas de temps saturne     : %4.2le \n", dt_sat);
    printf("pas de temps aster       : %4.2le \n", dt_ast);

    /* choix du pdt minimum : dt = dt_ast; */
    dt = dtref;
    if (dt > dt_ast) {
      dt = dt_ast;
    }
    if (dt > dt_sat) {
      dt = dt_sat;
    }

    /* envoi du pas de temps retenu */
    if (ierr >= 0) ierr = send_pdt(icompo, dt, i);

    printf("pas de temps retenu      : %4.2le \n", dt);
    printf("----------------------------------\n");
    printf("\n\n");

    icv = 0;

    j = 1;
    while (ierr >= 0) {

      printf("*********************************\n");
      printf("*    sous - iteration %i        *\n", j);
      printf("*********************************\n");
      printf("\n\n");

      /* incrémentation itération de couplage */
      printf("milieu\n");
      ntcast = ntcast + 1;
      printf("ntcast = %i\n", ntcast);

      /* printf("***************************************\n"); */
      /* printf("*     prediction des deplacements     *\n"); */
      /* printf("***************************************\n"); */

      /* Prédiction des déplacements */

      c1 = 0.;
      c2 = 0.;
      c3 = 0.;

      /* distinction cas explicite / cas implicite  pour la prédiction */
      if (j == 1) {
        alpha = 0.5;
        beta  = 0.;
        c1    = 1.;
        c2    = (alpha + beta) * dt ;
        c3    = -beta * dtold ;
        pred(xastp, xast, xvast, xvasa, c1, c2, c3, nb_dyn);
      }

      if (j > 1) {
        alpha = 0.5;
        c1    = alpha;
        c2    = 1. - alpha ;
        c3    = 0.;
        pred(xastp, xast, xastp, xast, c1, c2, c3, nb_dyn);
      }

      printf("--------------------------------------------\n");
      printf("Coefficients de prediction des deplacements \n");
      printf(" C1: %4.2le\n", c1);
      printf(" C2: %4.2le\n", c2);
      printf(" C3: %4.2le\n", c3);
      printf("--------------------------------------------\n");
      printf("\n\n");

      /* envoi des déplacements prédits */
      if (ierr >= 0) ierr = send_dyn(icompo);

      /* cas explicite: pas besoin de faire un test de convergence */

      /* cas implicite: nécessite un test de convergence */

      /* printf("***************************************\n"); */
      /* printf("*  fin prediction des deplacements    *\n"); */
      /* printf("***************************************\n"); */

      /* printf("*********************************\n"); */
      /* printf("*     prediction des forces     *\n"); */
      /* printf("*********************************\n"); */

      /* réception des  forces */
      ierr = recv_for(icompo);

      /* pas de distinction cas explicite et cas implicite pour les forces */
      alpha = 2.0;
      c1    = alpha;
      c2    = 1-alpha;
      c3    = 0.;
      pred(fopas, foras, foaas, foaas, c1, c2, c3, nb_for);
      printf("--------------------------------------\n");
      printf("Coefficients de prediction des forces \n");
      printf(" C1: %4.2le\n",c1);
      printf(" C2: %4.2le\n",c2);
      printf(" C3: %4.2le\n",c3);
      printf("--------------------------------------\n");
      printf("\n\n");

      /* envoi des forces */
      if (ierr >= 0) ierr = send_for(icompo);

      /* printf("*********************************\n"); */
      /* printf("*   fin prediction des forces   *\n"); */
      /* printf("*********************************\n"); */

      printf("\n");

      /* cas explicite: pas besoin de faire un test de convergence */
      if (nbssit <= 1) {
        /* gestion de la convergence même si pas de test */
        icv =1;
        if (ierr >= 0) ierr = send_icv1(icompo,icv);
        if (ierr >= 0) ierr = recv_icv(icompo,&(icv));
        icv = 1;
        if (ierr >= 0) ierr = send_icv2(icompo,icv);

        /* réception des déplacements calculés effectivement par Code_Aster */
        if (ierr >= 0) ierr = recv_dyn(icompo);

        /* enregistrement des valeurs antérieures */
        val_ant();

        break;
      }

      /* cas implicite: necessite un test de convergence */
      else {
        /* calcul de icv */
        if (ierr >= 0) ierr = conv(&(icv));
        if (ierr >= 0) ierr = send_icv1(icompo,icv);
        if (ierr >= 0) ierr = recv_icv(icompo,&(icv));
        if (ierr >= 0) ierr = send_icv2(icompo,icv);

        if((j>=nbssit) || (icv == 1)) {
          /* réception des déplacements calculés effectivement par Code_Aster */
          /* Réception des déplacements */
          if (ierr >= 0) ierr = recv_dyn(icompo);

          /* et envoi vers Code_Saturne ? la question demeure ? */
          /* si nécessaire routine d'envoi de ces deps à créer dans milieu */
          /* et homologue de réception dans Code_Saturne */
          /* if (ierr >= 0) ierr = send2_dyn(); */

          /* réception des déplacements calculés effectivement par Code_Aster */
          if (ierr >= 0) ierr = recv_dyn(icompo);
          break;
        }
        else {
          j = j+1;
        }
      }
    } /* fin boucle sous itération */

    /* test iterations */
    if (i >= nbpdtm) {
      ierr = -1;
    }
    /* fin test iterations */

    i = i+1;

    /* sauvegarde pas de temps */
    dtold = dt;

  } /* fin boucle itération */

  ierr = calfin(icompo);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif
