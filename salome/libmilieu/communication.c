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
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * SALOME headers
 *----------------------------------------------------------------------------*/

#include <calcium.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "runmilieu.h"
#include "donnees.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "communication.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction recv_geom
 *
 * Réceptionne les données geométriques et renseigne les variables
 * nb_for et nb_dyn
 *----------------------------------------------------------------------------*/

int
recv_geom(void *component)
{
  /* variables locales */
  int    i = 0;
  int    ii = 0;
  int    iret = 0;
  char   nomvar[144];
  int    geom[2] = {0, 0};
  float  tps = 0.;
  double tps2=0.;
  double almloc = 0.;

  printf("dans recv_geom \n");

  /* Initialisations */
  nb_for = 0;
  nb_dyn = 0;
  lref   = 0.0;

  strcpy(nomvar, "DONGEO");

  /* commande de reception des variables geometriques
   * 1 tableau de taille 1 * 2 est recu
   * geom[0] = nb_for; geom[1] = nb_dyn */

  i = 0;
  iret = cp_len(component,
                CP_ITERATION,
                &(tps),
                &(tps),
                &(i),
                nomvar,
                2,
                &(ii),
                &(geom[0]));

  if (iret < 1) {
    strcpy(nomvar, "ALMAXI");
    /* commande de reception de la variable lref */
    i = 0;
    iret = cp_ldb(component,
                  CP_ITERATION,
                  &(tps2),
                  &(tps2),
                  &(i),
                  nomvar,
                  1,
                  &(ii),
                  &(almloc));
    lref = almloc;
  }

  nb_for = geom[0];
  nb_dyn = geom[1];

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction send_geom
 *
 * Envoi les données géometriques à Code_Aster et renseigne
 * À supprimer lors de la phase d'initialisation
 *----------------------------------------------------------------------------*/

int
send_geom(void* component)
{
  /* Variables locales */
  int iret = 0;
  char nomvar[] = "NB_DYN";

  /* Test impression passage dans la routine */
  printf("dans recv_geom\n");

  /* envoi du nombre de points Code_Saturne à Code_Aster */
  iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(nb_dyn));
  if (iret < 1) {
    strcpy(nomvar,"NB_FOR");
    /* commande d'envoi de la variable nbssit */
    iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(nb_for));
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction send_pdt
 *
 * Envoie le pas de temps calculé par le composant Milieu aux codes
 *----------------------------------------------------------------------------*/

int
send_pdt(void *component,
         double dt,
         int numpdt)
{
  /* Variables locales */
  int iret = 0;
  char nomvar[] = "DTCALC";

  /* Test impression passage dans la routine */
  printf("dans send_pdt\n");

  /* envoi du pas de temps dtref aux codes */
  iret = cp_edb(component, CP_ITERATION, 0.0, numpdt, nomvar, 1, &(dt));

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction send_pdt
 *
 * Reçoit les pas de temps venant de Code_Aster et de Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_pdt(void *component,
         double *dt_ast,
         double *dt_sat,
         int numpdt)
{
  /* Variables locales */
  int    iret = 0;
  char   nomvar[] = "DTSAT";
  int    i, ii;
  double tps = 0.;
  double dtloc = 0.;

  /* Test impression passage dans la routine */
  printf("dans recv_pdt\n");

  i = numpdt;

  if (iret < 1) {
    strcpy(nomvar, "DTSAT");
    /* commande de réception du pas de temps Code_Saturne */
    i = numpdt;
    iret = cp_ldb(component,
                  CP_ITERATION,
                  &(tps),
                  &(tps),
                  &(i),
                  nomvar,
                  1,
                  &(ii),
                  &(dtloc));
    *(dt_sat) = dtloc;
  }

  if (iret < 1) {
    strcpy(nomvar, "DTSAT");
    /* commande de reception du pas de temps aster */
    i = numpdt;
    iret = cp_ldb(component,
                  CP_ITERATION,
                  &(tps),
                  &(tps),
                  &(i),
                  nomvar,
                  1,
                  &(ii),
                  &(dtloc));
    *(dt_ast)=dtloc;
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction send_param
 *
 * Envoie les donnees suivantes:
 *                               nbpdtm
 *                               nbssit
 *                               epsilo
 *                               isyncp
 *                               ntchr
 *                               ttpabs
 *----------------------------------------------------------------------------*/

int
send_param(void* component)
{
  /* Variables locales */
  char nomvar[] = "NBPDTM";
  int iret = 0;

  printf("dans send_param \n");

  /* commande d'envoi de la variable nbpdtm */
  iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(nbpdtm));

  if (iret < 1) {
    strcpy(nomvar, "NBSSIT");

    /* commande d'envoi de la variable nbssit */
    iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(nbssit));
  }

  if (iret < 1) {
    strcpy(nomvar, "EPSILO");
    /* commande d'envoi de la variable epsilo */
    iret = cp_edb(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(epsilo));
  }

  if (iret < 1) {
    strcpy(nomvar, "ISYNCP");
    /* commande d'envoi de la variable isyncp */
    iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(isyncp));
  }

  if (iret < 1) {
    strcpy(nomvar, "NTCHRO");
    /* commande d'envoi de la variable ntchr  */
    iret = cp_een(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(ntchr));
  }

  if (iret < 1) {
    strcpy(nomvar,"TTINIT");
    /* commande d'envoi de la variable ttinit */
    iret = cp_edb(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(ttinit));
  }

  if (iret < 1) {
    strcpy(nomvar, "PDTREF");
    /* commande d'envoi de la variable dtref */
    iret = cp_edb(component, CP_ITERATION, 0.0, 0, nomvar, 1, &(dtref));
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction recv_dyn
 *
 * Reçoit les déplacements et les vitesses venant de Code_Aster
 * au pas de temps courant
 *----------------------------------------------------------------------------*/

int
recv_dyn(void *component)
{
  /* Variables locales */
  char nomvar[] = "DEPAST";
  int iret = 0;
  int i, ii;
  double tps = 0.;

  printf("dans recv_dyn \n");

  /* commande de reception de la variable depast */
  i = ntcast;
  iret = cp_ldb(component,
                CP_ITERATION,
                &(tps),
                &(tps),
                &(i),
                nomvar,
                3*nb_dyn,
                &(ii),
                xast);

  if (iret < 1) {
    strcpy(nomvar, "VITAST");
    /* commande de reception de la variable vitast */
    i = ntcast;
    iret = cp_ldb(component,
                  CP_ITERATION,
                  &(tps),
                  &(tps),
                  &(i),
                  nomvar,
                  3*nb_dyn,
                  &(ii),
                  xvast);
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction send_dyn
 *
 * Envoie les déplacements prédits à Code_Saturne
 *----------------------------------------------------------------------------*/

int
send_dyn(void *component)
{
  /* Variables locales */
  char nomvar[] = "DEPSAT";
  int iret = 0;

  printf("dans send_dyn \n");

  iret = cp_edb(component,
                CP_ITERATION,
                0.0,
                ntcast,
                nomvar,
                3*nb_dyn,
                xastp);

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction recv_for
 *
 * Reçoit les efforts venant de Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_for(void *component)
{
  /* Variables locales */
  char nomvar[] = "FORSAT";
  int iret = 0;
  double tps;
  int i, ii;

  printf("dans recv_for \n");

  i = ntcast;
  iret = cp_ldb(component,
                CP_ITERATION,
                &(tps),
                &(tps),
                &(i),
                nomvar,
                3*nb_for,
                &(ii),
                foras);

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction recv_for
 *
 * Envoie les efforts prédits vers Code_Aster
 *----------------------------------------------------------------------------*/

int
send_for(void *component)
{
  /* Variables locales */
  char nomvar[] = "FORAST";
  int iret = 0;
  int i;

  printf("dans send_for \n");

  i = ntcast;
  iret = cp_edb(component, CP_ITERATION, 0.0, i, nomvar, 3*nb_for, fopas);

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction send_icv1
 *
 * Envoie l'indice de convergence à Code_Saturne
 *----------------------------------------------------------------------------*/

int
send_icv1(void *component,
          int icv)
{
  /* Variables locales */
  int iret;
  char nomvar[] = "ICVEXT";

  printf("dans send_icv1 \n");

  /* commande d'envoi de la variable icvext */
  iret =  cp_een(component, CP_ITERATION, 0.0, ntcast, nomvar, 1, &(icv));

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction recv_icv
 *
 * Reçoit l'indice de convergence de Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_icv(void *component,
         int *icv)
{
  /* Variables locales */
  int iret;
  char nomvar[] = "ICV";
  float tps=0.;
  int i, ii;

  /* commande d'envoi de la variable icv */
  i = ntcast;
  iret = cp_len(component,
                CP_ITERATION,
                &(tps),
                &(tps),
                &(i),
                nomvar,
                1,
                &(ii),
                icv);
  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction send_icv2
 *
 * Envoie l'indice de convergence à Code_Aster
 *----------------------------------------------------------------------------*/

int
send_icv2(void *component,
          int icv)
{
  /* Variables locales */
  int iret;
  char nomvar[] = "ICVAST";

  printf("dans send_icv2 %d \n",icv);

  /* commande d'envoi de la variable icvext */
  iret = cp_een(component, CP_ITERATION, 0.0, ntcast, nomvar, 1, &(icv));

  printf("iret = %d\n", iret);

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction inicom
 *
 * Initialisation de la communication avec Calcium
 *----------------------------------------------------------------------------*/

int
inicom(void *component)
{
  char instance[200];

  int iret = cp_cd(component, instance);

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction calfin
 *
 * Fin de la communication avec Calcium et arrêt du calcul
 *----------------------------------------------------------------------------*/

int
calfin(void *component)
{
  int iret = cp_fin(component, CP_ARRET);

  exit(0);

  return iret;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif
