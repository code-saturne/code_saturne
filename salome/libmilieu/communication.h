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

#ifndef __COMMUNICATION_H__
#define __COMMUNICATION_H__

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction recv_geom
 *
 * Réceptionne les données geométriques et renseigne les variables
 * nb_for, nb_dyn et lref
 *----------------------------------------------------------------------------*/

int
recv_geom(void *component);

/*----------------------------------------------------------------------------
 * Fonction send_geom
 *
 * Envoi les données géometriques à Code_Aster et renseigne
 * À supprimer lors de la phase d'initialisation
 *----------------------------------------------------------------------------*/

int
send_geom(void* component);

/*----------------------------------------------------------------------------
 * Fonction send_pdt
 *
 * Envoie le pas de temps calculé par le composant Milieu aux codes
 *----------------------------------------------------------------------------*/

int
send_pdt(void *component,
         double dt,
         int numpdt);

/*----------------------------------------------------------------------------
 * Fonction send_pdt
 *
 * Reçoit les pas de temps venant de Code_Aster et de Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_pdt(void *component,
         double *dt_ast,
         double *dt_sat,
         int numpdt);

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
send_param(void* component);

/*----------------------------------------------------------------------------
 * Fonction recv_dyn
 *
 * Reçoit les déplacements et les vitesses venant de Code_Aster
 * au pas de temps courant
 *----------------------------------------------------------------------------*/

int
recv_dyn(void *component);

/*----------------------------------------------------------------------------
 * Fonction send_dyn
 *
 * Envoie les déplacements prédits à Code_Saturne
 *----------------------------------------------------------------------------*/

int
send_dyn(void *component);

/*----------------------------------------------------------------------------
 * Fonction recv_for
 *
 * Reçoit les efforts venant de Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_for(void *component);

/*----------------------------------------------------------------------------
 * Fonction recv_for
 *
 * Envoie les efforts prédits vers Code_Aster
 *----------------------------------------------------------------------------*/

int
send_for(void *component);

/*----------------------------------------------------------------------------
 * Fonction send_icv1
 *
 * Envoie l'indice de convergence à Code_Saturne
 *----------------------------------------------------------------------------*/

int
send_icv1(void *component,
          int icv);

/*----------------------------------------------------------------------------
 * Fonction recv_icv
 *
 * Reçoit l'indice de convergence de Code_Saturne
 *----------------------------------------------------------------------------*/

int
recv_icv(void *component,
         int *icv);

/*----------------------------------------------------------------------------
 * Fonction send_icv2
 *
 * Envoie l'indice de convergence à Code_Aster
 *----------------------------------------------------------------------------*/

int
send_icv2(void *component,
          int icv);

/*----------------------------------------------------------------------------
 * Fonction inicom
 *
 * Initialisation de la communication avec Calcium
 *----------------------------------------------------------------------------*/

int
inicom(void *component);

/*----------------------------------------------------------------------------
 * Fonction calfin
 *
 * Fin de la communication avec Calcium et arrêt du calcul
 *----------------------------------------------------------------------------*/

int
calfin(void *component);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif /* __COMMUNICATION_H__ */
