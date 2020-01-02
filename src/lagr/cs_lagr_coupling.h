#ifndef __CS_LAGR_COUPLING_H__
#define __CS_LAGR_COUPLING_H__

/*============================================================================
 * Functions and types for the Lagrangian module
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 *     CALCUL DES TERMES SOURCES DU COUPLAGE RETOUR
 *     Remarque : les termes sources sont calcules pour
 *                la cellule de depart de la particule
 *                lors de l'iteration courante. Attention, meme
 *                si la particule est sortante du domaine de
 *                calcul (peu importe la maniere) on doit calculer
 *                un terme source qui correspond a ce qu'echange le
 *                fluide porteur et la particule au debut du pas de
 *                temps. Si NORDRE = 2 et que la particule est en
 *                interaction avec la frontiere, alors les termes
 *                source sont calcules comme si NORDRE=1
 *                (on oublie le pre-remplissage de TSFEXT dans
 * ONFC                 LAGES2).
 * --------------------------------------------------------------------
 * Arguments
 *
 *   ntersl            <--  nbr termes sources de couplage retour
 *
 *   taup(nbpart)      <--  temps caracteristique dynamique
 *   tsfext(nbpart)    <--  forces externes
 *   tempct            <--  temps caracteristique thermique
 *    (nbpart,2)
 *   cpgd1,cpgd2,      <--  termes de devolatilisation 1 et 2 et
 *    cpght(nbpart)           de combusion heterogene (charbon
 *                            avec couplage retour thermique)
 *   volp(ncelet)      ---  fraction volumique des particules
 *   volm(ncelet)      ---  fraction massique des particules
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_coupling(cs_real_t taup[],
                 cs_real_t tempct[],
                 cs_real_t tsfext[],
                 cs_real_t cpgd1[],
                 cs_real_t cpgd2[],
                 cs_real_t cpght[],
                 cs_real_t volp[],
                 cs_real_t volm[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LAGR_COUPLING_H__ */
