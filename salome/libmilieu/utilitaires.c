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
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "runmilieu.h"
#include "donnees.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "utilitaires.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * Public function definitions
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
     int nbpts)
{
  /* Variables locales */
  int i;

  /* Mise à jour du tableau de prediction */
  for (i = 0; i < nbpts; i++) {
    valpre[3*i]     = c1*val1[3*i]     + c2*val2[3*i]     + c3*val3[3*i];
    valpre[(3*i)+1] = c1*val1[(3*i)+1] + c2*val2[(3*i)+1] + c3*val3[(3*i)+1];
    valpre[(3*i)+2] = c1*val1[(3*i)+2] + c2*val2[(3*i)+2] + c3*val3[(3*i)+2];
  }
}

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
       double nbpts)
{
  /* Variables locales */
  int i;
  double norme;

  /* Calcul de la norme de la difference */
  norme = 0.;
  for (i = 0; i < nbpts; i++) {
    norme = norme + (vect1[3*i]-vect2[3*i])*(vect1[3*i]-vect2[3*i]);
    norme = norme + (vect1[3*i+1]-vect2[3*i+1])*(vect1[3*i+1]-vect2[3*i+1]);
    norme = norme + (vect1[3*i+2]-vect2[3*i+2])*(vect1[3*i+2]-vect2[3*i+2]);
  }
  norme = sqrt(norme/nbpts);
  return norme;
}

/*----------------------------------------------------------------------------
 * Fonction alldyn
 *
 * Réalise l'allocation et l'initialisation des vecteurs dynamiques (double)
 * sur la base du nombre de points 'nb_dyn'.
 *----------------------------------------------------------------------------*/

void
alldyn(void)
{
  /* Variables locales */
  int k;

  xast =(double *)calloc(3*nb_dyn,sizeof(double));
  xvast=(double *)calloc(3*nb_dyn,sizeof(double));
  xvasa=(double *)calloc(3*nb_dyn,sizeof(double));
  xastp=(double *)calloc(3*nb_dyn,sizeof(double));

  for (k = 0; k < nb_dyn; k++) {

    xast[3*k]   = 0.;
    xast[3*k+1] = 0.;
    xast[3*k+2] = 0.;

    xvast[3*k]   = 0.;
    xvast[3*k+1] = 0.;
    xvast[3*k+2] = 0.;

    xvasa[3*k]   = 0.;
    xvasa[3*k+1] = 0.;
    xvasa[3*k+2] = 0.;

    xastp[3*k]   = 0.;
    xastp[3*k+1] = 0.;
    xastp[3*k+2] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * Fonction allfor
 *
 * Réalise l'allocation et l'initialisation des vecteurs dynamiques (double)
 * sur la base du nombre de points 'nb_for'.
 *----------------------------------------------------------------------------*/

void
allfor(void)
{
  /* Variables locales */
  int k;

  foras =(double *)calloc(3*nb_for,sizeof(double));
  foaas =(double *)calloc(3*nb_for,sizeof(double));
  fopas =(double *)calloc(3*nb_for,sizeof(double));

  for (k = 0; k < nb_for; k++) {

    foras[3*k]   = 0.;
    foras[3*k+1] = 0.;
    foras[3*k+2] = 0.;

    foaas[3*k]   = 0.;
    foaas[3*k+1] = 0.;
    foaas[3*k+2] = 0.;

    fopas[3*k]   = 0.;
    fopas[3*k+1] = 0.;
    fopas[3*k+2] = 0.;
  }
}

/*----------------------------------------------------------------------------
 * Fonction conv
 *
 * Réalise le test de convergence en cas de calcul implicite
 *
 * renvoie: 0 si non-convergence
 *          1 si     convergence
 *----------------------------------------------------------------------------*/

int
conv(int *icv)
{
  /* Variables locales */
  int iret;
  double delast = 0.;

  if (lref > 0.) {

    delast = (dinorm(xast, xastp, nb_dyn))/lref;

    printf("-------------------------------------\n");
    printf("test de convergence:\n");
    printf("delast = %4.2le\n", delast);

    if (delast <= epsilo) {
      *(icv) = 1;
      printf("icv = %i\n", *(icv));
      printf("convergence de la sous iteration\n");
      printf("-------------------------------------\n");
    }
    else {
      printf("icv = %i\n", *(icv));
      printf("non convergence de la sous iteration\n");
      printf("-------------------------------------\n");
    }

    iret = 0;
  }
  else {
    printf("la valeur de lref negative ou nulle\n");
    printf("le calcul est interrompu\n");
    printf("------------------------------------------\n");
    iret = -1;
  }

  return iret;
}

/*----------------------------------------------------------------------------
 * Fonction val_ant
 *
 * Écrase les données de la sous iter k-1 avec les données de la sous iter k
 * données dynamiques : vitesses
 * efforts            : forces
 *----------------------------------------------------------------------------*/

void
val_ant(void)
{
  /* Variables locales */
  int i;

  /* enregistrement des efforts */
  for (i = 0; i< nb_for; i++) {
    foaas[3*i]   = foras[3*i];
    foaas[3*i+1] = foras[3*i+1];
    foaas[3*i+2] = foras[3*i+2];
  }

  /* enregistrement des donnees dynamiques */
  for (i = 0; i< nb_dyn; i++) {
    xvasa[3*i]   = xvast[3*i];
    xvasa[3*i+1] = xvast[3*i+1];
    xvasa[3*i+2] = xvast[3*i+2];
  }
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif
