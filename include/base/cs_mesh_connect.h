/*============================================================================
*
*                    Code_Saturne version 1.3
*                    ------------------------
*
*
*     This file is part of the Code_Saturne Kernel, element of the
*     Code_Saturne CFD tool.
*
*     Copyright (C) 1998-2007 EDF S.A., France
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

#ifndef __CS_MESH_CONNECT_H__
#define __CS_MESH_CONNECT_H__

/*============================================================================
 * Passage d'une connectivité noyau à une connecitvité nodale de la
 * structure principale associée à un maillage
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fichiers `include' locaux
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_mesh.h"

#include "fvm_nodal.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Définitions de macros
 *============================================================================*/


/*============================================================================
 * Définitions de types
 *============================================================================*/


/*=============================================================================
 * Variables globales_statiques
 *============================================================================*/


/*============================================================================
 *  Fonctions publiques pour API Fortran
 *============================================================================*/


/*=============================================================================
 * Prototypes de fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Extraction de la connectivité "cellules -> faces" d'un maillage.
 *
 * On considère une numérotation commune des faces internes et des
 * faces de bord, dans laquelle les faces de bord sont définies en
 * premier. L'indice commun de la i-ème face de bord est donc égal à i,
 * et celui de la j-ième face interne à nbr_fbr + j.
 *
 * Si ind_cel_extr != NULL, alors :
 * --- ind_cel_extr[icel] = indice dans la liste à extraire (0 à n-1)
 *     si icel correspond à une cellule à extraire
 * --- ind_cel_extr[icel] = -1 si la cellule icel est à ignorer
 *----------------------------------------------------------------------------*/

void cs_maillage_ret_cel_fac
(
 const cs_mesh_t       *const maillage,       /* --> Maillage */
 const cs_int_t               nbr_cel_extr,   /* --> Taille de ind_cel_extr[] */
 const cs_int_t               ind_cel_extr[], /* --> ind_cel_extr[cellule]
                                               *     = indice cellule extraite
                                               *       ou -1 */
 cs_int_t            * *const p_pos_cel_fac,  /* <-- idx cellule -> face */
 cs_int_t            * *const p_val_cel_fac   /* <-- val cellule -> face */
);


/*----------------------------------------------------------------------------
 * Extraction et conversion en connectivité nodale externe d'un sous-ensemble
 * des cellules d'un maillage.
 *
 * La liste des cellules à traiter est optionnelle ; elle peut ne pas
 * être ordonnée en entrée, elle le sera toujours en sortie (les cellules
 * étant extraites au cours d'un parcours en ordre croissant, la liste
 * est réordonnée pour assurer la cohérence des liens des cellules extraites
 * vers leurs cellules parentes, construits à partir de cette liste).
 *----------------------------------------------------------------------------*/

fvm_nodal_t  * cs_maillage_extrait_cel_nodal
(
 const cs_mesh_t      *const maillage,      /* --> maillage                   */
 const char           *const nom,           /* --> nom à affecter             */
 const cs_int_t              nbr_liste_cel, /* --> taille de liste_cel[]      */
       cs_int_t              liste_cel[]    /* <-> liste optionnelle des
                                             *     cellules à traiter (1 à n) */
);


/*----------------------------------------------------------------------------
 * Extraction et conversion en connectivité nodale externe d'un sous-ensemble
 * des faces d'un maillage.
 *
 * Les listes des faces à traiter sont optionnelles (si aucune des deux
 * n'est fournie, on extrait les faces de bord par défaut); elle peuvent
 * ne pas être ordonnées en entrée, elle le seront toujours en sortie
 * (les faces étant extraites au cours d'un parcours en ordre croissant,
 * la liste est réordonnée pour assurer la cohérence des liens des faces
 * extraites vers leurs faces parentes, construites à partir de cette liste).
 *----------------------------------------------------------------------------*/

fvm_nodal_t  * cs_maillage_extrait_fac_nodal
(
 const cs_mesh_t      *const maillage,      /* --> maillage                   */
 const char           *const nom,           /* --> nom à affecter             */
 const cs_int_t              nbr_liste_fac, /* --> taille de liste_fac[]      */
 const cs_int_t              nbr_liste_fbr, /* --> taille de liste_fbr[]      */
       cs_int_t              liste_fac[],   /* <-> liste optionnelle des faces
                                             *     internes à traiter (1 à n) */
       cs_int_t              liste_fbr[]    /* <-> liste optionnelle des faces
                                             *     de bord à traiter (1 à n)  */
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_MESH_CONNECT_H__ */
