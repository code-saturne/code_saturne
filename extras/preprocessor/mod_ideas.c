/*============================================================================
 *  Modification des coordonnées d'un fichier de maillage
 *  IDEAS-MS au format "universel"
 *============================================================================*/

/*
  This file is part of the code_saturne Preprocessor, element of the
  code_saturne CFD tool.

  Copyright (C) 1999-2007 EDF S.A., France

  contact: saturne-support@edf.fr

  The code_saturne Preprocessor is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The code_saturne Preprocessor is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the code_saturne Preprocessor; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
*/


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Pour une lecture de 80 caracteres par ligne  */
/* auxquels il faut ajouter le `\n' et le `\0'  */
/* pour l'affectation dans la chaine receptrice */
#define LNG_MAX_CHAINE_IDEAS  82                /* Dimension des chaînes */

typedef enum {
  HORS_DATASET,
  MARQUE_DEB,
  NUM_DATASET,
  CONTENU,
  MARQUE_FIN
} ligne_data_t;


/*----------------------------------------------------------------------------
 *  Fonction utilitaire pour la lecture/écriture des réels d'un fichier
 *  I-DEAS.
 *
 *  La chaîne contenant un réel avec (ou non) un exposant `d' ou `D'
 *  est convertie en réel avec un exposant `e'
 *----------------------------------------------------------------------------*/

static void transf_geo
(
 double *x,
 double *y,
 double *z
)
{
  *x = (*x) * 1.0;
  *y = (*y) * 1.0;
  *z = (*z) * 1.0;
}

/*============================================================================
 *                             Fonction principale
 *============================================================================*/

int main
(
 int    argc   ,
 char * argv[]
)
{

  char  chaine[LNG_MAX_CHAINE_IDEAS + 2];
  char *s;
  FILE *fic_lec, *fic_ecr;
  int   ret;

  char         *res = NULL;
  ligne_data_t ind_rub = HORS_DATASET;
  int          ds_nod = 0;
  int          ind_pair = 0;


  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if  (argc < 3) {
    printf("Utilisation :\n%s nom_fic nom_fic_conv\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  fic_lec = fopen(argv[1], "rb");
  if (fic_lec == NULL) {
    fprintf(stderr, "Erreur à l'ouverture en lecture du fichier %s\n",
	    argv[1]);
    return EXIT_FAILURE;
  }
  fic_ecr = fopen(argv[2], "wb");
  if (fic_ecr == NULL) {
    fprintf(stderr, "Erreur à l'ouverture en écriture du fichier %s\n",
	    argv[2]);
    return EXIT_FAILURE;
  }

  while (ds_nod == 0 && feof(fic_lec) == 0) {

    res = fgets(chaine, LNG_MAX_CHAINE_IDEAS, fic_lec) ;

    if (res != NULL) {

      /* Début ou fin de rubrique */
      if ((ind_rub == HORS_DATASET || ind_rub == CONTENU)
	  && strlen(chaine) < 8) {
	if (strncmp(chaine, "    -1", 6) == 0) {
	  if (ind_rub == HORS_DATASET)
	    ind_rub = MARQUE_DEB;
	  else
	    ind_rub = MARQUE_FIN;
	}
      }

      /* Numéro de rubrique */
      else if (ind_rub == NUM_DATASET) {
	if (strncmp(chaine, "  2411", 6) == 0) {
	  fprintf(fic_ecr, "%s", chaine);
	  break;
	}
      }

      fprintf(fic_ecr, "%s", chaine);

      /* Préparation ligne suivante */
      switch(ind_rub) {
      case MARQUE_DEB:
	ind_rub = NUM_DATASET;
	break;
      case NUM_DATASET:
	ind_rub = CONTENU;
	break;
      case MARQUE_FIN:
	ind_rub = HORS_DATASET;
	break;
      default:
	break;
      }

    }

  }

  while (feof(fic_lec) == 0) {

    res = fgets(chaine, LNG_MAX_CHAINE_IDEAS, fic_lec) ;

    if (res != NULL) {

      if (strlen(chaine) < 8 && strncmp(chaine, "    -1", 6) == 0) {
	fprintf(fic_ecr, "%s", chaine);
	break;
      }

      /* Ligne paire ou impaire ? */

      if (ind_pair == 0)
	ind_pair = 1;

      else {

	char   ch_coord[3][LNG_MAX_CHAINE_IDEAS];
	double coord[3];
	ind_pair = 0;

	s = chaine ;

	while(*s != '\0' ) {
	  if (*s == 'D' || *s == 'd')
	    *s = 'e' ;
	  s++ ;
	}

	ret = sscanf(chaine, " %lf %lf %lf", coord, coord+1, coord+2);

	transf_geo(coord, coord+1, coord+2);

	sprintf(chaine, " %24.16E %24.16E %24.16E\n",
		coord[0], coord[1], coord[2]);

	s = chaine ;

	while(*s != '\0' ) {
	  if (*s == 'E' || *s == 'e')
	    *s = 'D' ;
	  s++ ;
	}

      }

      fprintf(fic_ecr, "%s", chaine);

    }

  }

  /* Recopie simple de la suite du fichier */
  while (feof(fic_lec) == 0) {
    res = fgets(chaine, LNG_MAX_CHAINE_IDEAS, fic_lec) ;
    if (res != NULL)
      fprintf(fic_ecr, "%s", chaine);
  }

  /* Fin */
  fclose(fic_lec);
  fclose(fic_ecr);

  return EXIT_SUCCESS ;

}
