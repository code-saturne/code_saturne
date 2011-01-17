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

#ifndef __UTILITAIRES_H__
#define __UTILITAIRES_H__

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fonction pred
 *
 * Réalise la prédiction du déplacement ou des forces à partir des valeurs
 * aux pas de temps courant et precedent(s)
 *
 * valpre = c1 * val1 + c2 * val2 + c3 * val3
 *----------------------------------------------------------------------------*/

void
pred(double *valpre,
     double *val1,
     double *val2,
     double *val3,
     double c1,
     double c2,
     double c3,
     int nbpts);

/*----------------------------------------------------------------------------
 * Fonction dinorm
 *
 * Calcule la norme de la différence entre les vecteurs vect1 et vect2
 *
 * dinorm = sqrt(somme sur nbpts i
 *                 (somme sur composante j
 *                    ((vect1[i,j]-vect2[i,j])^2)))
 *----------------------------------------------------------------------------*/

double
dinorm(double *vect1,
       double *vect2,
       double nbpts);

/*----------------------------------------------------------------------------
 * Fonction alldyn
 *
 * Réalise l'allocation et l'initialisation des vecteurs dynamiques (double)
 * sur la base du nombre de points 'nb_dyn'.
 *----------------------------------------------------------------------------*/

void
alldyn(void);

/*----------------------------------------------------------------------------
 * Fonction allfor
 *
 * Réalise l'allocation et l'initialisation des vecteurs dynamiques (double)
 * sur la base du nombre de points 'nb_for'.
 *----------------------------------------------------------------------------*/

void
allfor(void);

/*----------------------------------------------------------------------------
 * Fonction conv
 *
 * Réalise le test de convergence en cas de calcul implicite
 *
 * renvoie: 0 si non-convergence
 *          1 si     convergence
 *----------------------------------------------------------------------------*/

int
conv(int *icv);

/*----------------------------------------------------------------------------
 * Fonction val_ant
 *
 * Écrase les données de la sous iter k-1 avec les données de la sous iter k
 * données dynamiques : vitesses
 * efforts            : forces
 *----------------------------------------------------------------------------*/

void
val_ant(void);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

#endif /* __UTILITAIRES_H__ */
