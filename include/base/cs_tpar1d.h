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

#ifndef __CS_TPAR1D_H__
#define __CS_TPAR1D_H__

/*============================================================================
 *  Gestion des fichiers suite
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*----------------------------------------------------------------------------
 *  Fichiers `include' librairie standard C
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Fichiers `include' locaux
 *----------------------------------------------------------------------------*/

#include "cs_base.h"


/*============================================================================
 *  Prototype de fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation des maillages de chaque face et initialisation de la temperature
 *----------------------------------------------------------------------------*/

void CS_PROCF (mait1d,MAIT1D)(cs_int_t * , cs_int_t *,
                              cs_real_t *, cs_real_t *, cs_real_t* );


/*----------------------------------------------------------------------------
 * Resolution de l'equation 1D pour une face donnee
 *----------------------------------------------------------------------------*/

void CS_PROCF (tpar1d,TPAR1D)(cs_int_t * , cs_int_t * , cs_real_t *,
                              cs_real_t *, cs_real_t *, cs_real_t *,
                              cs_real_t *, cs_real_t *, cs_real_t *,
                              cs_real_t *, cs_real_t *);


/*----------------------------------------------------------------------------
 * Lecture du fichier suite du module thermique 1D en paroi
 *----------------------------------------------------------------------------*/

void CS_PROCF (lect1d,LECT1D)
(
 const char       *const nomsui,  /* <- Nom du fichier suite               */
 const cs_int_t   *const lngnom,  /* <- Longueur du nom                       */
 const cs_int_t   *const ifovt1,  /* <- Indicateur binaire (0) / ascii (1)    */
 const cs_int_t   *const nfpt1d,  /* <- Nbr de  faces avec couplage           */
 const cs_int_t   *const nfpt1t,  /* <- Nbr de  faces avec couplage cumule sur
                                        tous les processeurs                  */
 const cs_int_t   *const nmxt1d,  /* <- Nbr max de pts sur les maillages 1D   */
 const cs_int_t   *const nfabor,  /* <- Nbr de faces de bord                  */
 const cs_int_t   *const nppt1d,  /* <- Nbr de points de discretisation des
                                                faces avec module 1D                  */
 const cs_int_t   *const ifpt1d,  /* <- Tableau d'indirection des faces avec
                                                module 1D                             */
 const cs_real_t  *const eppt1d,  /* <- Epaisseur de paroi des faces          */
 const cs_real_t  *const rgpt1d,  /* <- Raison geometrique associee aux faces */
       cs_real_t  *const tppt1d   /* -> Temperature de paroi avec module 1D   */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' eventuels F77, */
                                  /*     inutilises lors de l'appel mais      */
                                  /*     places par de nombreux compilateurs) */
 );


/*----------------------------------------------------------------------------
 * Ecriture du fichier suite du module thermique 1D en paroi
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecrt1d,ECRT1D)
(
 const char       *const nomsui,  /* <- Nom du fichier suite                  */
 const cs_int_t   *const lngnom,  /* <- Longueur du nom                       */
 const cs_int_t   *const ifovt1,  /* <- Indicateur binaire (0) / ascii (1)    */
 const cs_int_t   *const nfpt1d,  /* <- Nbr de  faces avec couplage           */
 const cs_int_t   *const nmxt1d,  /* <- Nbr max de pts sur les maillages 1D   */
 const cs_int_t   *const nfabor,  /* <- Nbr de faces de bord                  */
 const cs_real_t  *const tppt1d,  /* <- Temperature de paroi avec module 1D   */
 const cs_int_t   *const ifpt1d   /* <- Tableau d'indirection des faces avec
                                     module 1D                                */
 CS_ARGF_SUPP_CHAINE              /*     (arguments 'longueur' eventuels F77, */
                                  /*     inutilises lors de l'appel mais      */
                                  /*     places par de nombreux compilateurs) */
 );


/*----------------------------------------------------------------------------
 * Liberation de la memoire
 *----------------------------------------------------------------------------*/

void CS_PROCF (lbrt1d,LBRT1D)(void);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_TPAR1D_H__ */

