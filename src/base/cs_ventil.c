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

/*============================================================================
 * Définitions, variables globales, et fonctions associées aux ventilateurs
 *============================================================================*/

/* includes système */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Includes librairie BFT */

#include <bft_mem.h>

/* Includes librairie */

#include "cs_ventil.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#if 0 /* Fausse "}" pour corriger l'auto-indentation d'Emacs */
}
#endif


/*============================================================================
 * Définitions de types
 *============================================================================*/

/* Structure associée à un ventilateur */

struct _cs_ventil_t {

  int                 num;              /* Numéro du ventilateur */
  int                 dim_modele;       /* Modèle 1D, 2D, ou 3D */
  int                 dim_ventil;       /* Géométrie 2D ou 3D */

  cs_real_t           coo_axe_amont[3]; /* Coordonnées du point de l'axe
                                           de la face amont */
  cs_real_t           coo_axe_aval[3];  /* Coordonnées du point de l'axe
                                           de la face amont */
  cs_real_t           dir_axe[3];       /* Vecteur directeur unitaire de
                                           l'axe (amont vers aval) */
  cs_real_t           epaisseur;        /* Épaisseur du ventilateur */
  cs_real_t           surface;          /* Surface totale du ventilateur */

  cs_real_t           ray_ventil;       /* Rayon du ventilateur */
  cs_real_t           ray_pales;        /* Rayon des pales */
  cs_real_t           ray_moyeu;        /* Rayon du moyeu */
  cs_real_t           coeff_carac[3];   /* Coefficients des termes de
                                           degré 0, 1, et 2 de la courbe
                                           caractéristique */
  cs_real_t           couple_axial;     /* Couple axial du ventilateur*/

  cs_int_t            nbr_cel;

  cs_int_t           *lst_cel;          /* Liste des cellules appartenant au
                                           ventilateur */

  cs_real_t           debit_entrant;    /* Débit entrant courant */
  cs_real_t           debit_sortant;    /* Débit sortant courant */

} ;


/*============================================================================
 *  Variables globales statiques
 *============================================================================*/

/* Tableau des ventilateurs */

cs_int_t         cs_glob_ventil_nbr_max = 0;

cs_int_t         cs_glob_ventil_nbr = 0;
cs_ventil_t  * * cs_glob_ventil_tab = NULL;


/*============================================================================
 *  Définitions de macros
 *============================================================================*/

/* Définition de macros locales */

enum {X, Y, Z} ;

#define CS_LOC_PRODUIT_VECTORIEL(prod_vect, vect1, vect2)  \
  (prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z], \
   prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z], \
   prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y])

#define CS_LOC_PRODUIT_SCALAIRE(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z])

#define CS_LOC_MODULE(vect) \
  sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])


/*============================================================================
 * Prototypes de fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Marquage des cellules appartenant aux différents ventilateurs
 * (par le numéro de ventilateur, 0 sinon)
 *----------------------------------------------------------------------------*/

static void cs_loc_ventil_marque_cellules
(
 const cs_mesh_t  *const mesh,     /* <-- structure maillage associée */
       cs_int_t          num_vtl_cel[] /* --> indicateur par cellule      */
);


/*============================================================================
 *  Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Récupération du nombre de ventilateurs
 *
 * Interface Fortran :
 *
 * SUBROUTINE TSTVTL
 * *****************
 *
 * INTEGER          NBRVTL         : --> : nombre de ventilateurs
 *----------------------------------------------------------------------------*/

void CS_PROCF (tstvtl, TSTVTL)
(
 cs_int_t  *const nbrvtl              /* <-- nombre de ventilateurs           */
)
{
  *nbrvtl = cs_glob_ventil_nbr ;
}


/*----------------------------------------------------------------------------
 * Ajout d'un ventilateur
 *
 * Interface Fortran :
 *
 * SUBROUTINE DEFVTL (XYZVT1, XYZVT2, RVVT  , RPVT  , RMVT  ,
 * *****************
 *                    CCARAC, TAUVT)
 *
 * INTEGER          DIMMOD     : <-- : Dimension du modèle de ventilateur :
 *                             :     : f_constante ; 2 : profil_force ;
 *                             :     : 3 : profil_force + couple tangentiel
 *                  DIMVTL     : <-- : Dimension du ventilateur :
 *                             :     : 2 : pseudo-2D (maillage extrudé)
 *                             :     : 3 : 3D (standard)
 * DOUBLE PRECISION XYZVT1(3)  : <-- : Coord. point de l'axe en face amont
 * DOUBLE PRECISION XYZVT2(3)  : <-- : Coord. point de l'axe en face aval
 * DOUBLE PRECISION RVVT       : <-- : Rayon du ventilateur
 * DOUBLE PRECISION RPVT       : <-- : Rayon des pales
 * DOUBLE PRECISION RMVT       : <-- : Rayon du moyeu
 * DOUBLE PRECISION CCARAC(3)  : <-- : Coefficients de degré 0, 1, et 2
 *                             :     : de la courbe caractéristique
 * DOUBLE PRECISION TAUVT      : <-- : Couple axial du ventilateur
 *----------------------------------------------------------------------------*/

void CS_PROCF (defvtl, DEFVTL)
(
 const cs_int_t   *const dimmod,     /* Dimension du modèle de ventilateur :
                                        1 : f_constante ; 2 : profil_force ;
                                        3 : profil_force + couple tangentiel */
 const cs_int_t   *const dimvtl,     /* Dimension du ventilateur :
                                        2 : pseudo-2D (maillage extrudé)
                                        3 : 3D (standard) */
 const cs_real_t         xyzvt1[3],  /* Coord. point de l'axe en face amont */
 const cs_real_t         xyzvt2[3],  /* Coord. point de l'axe en face aval */
 const cs_real_t  *const rvvt,       /* Rayon du ventilateur */
 const cs_real_t  *const rpvt,       /* Rayon des pales */
 const cs_real_t  *const rmvt,       /* Rayon du moyeu */
 const cs_real_t         ccarac[3],  /* Coefficients courbe caractéristique */
 const cs_real_t  *const tauvt       /* Couple axial du ventilateur*/
)
{
  cs_ventil_definit(*dimmod,
                    *dimvtl,
                    xyzvt1,
                    xyzvt2,
                    *rvvt,
                    *rpvt,
                    *rmvt,
                    ccarac,
                    *tauvt);
}


/*----------------------------------------------------------------------------
 * Construction des listes de cellules associées aux ventilateurs
 *
 * Interface Fortran :
 *
 * SUBROUTINE INIVTL
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (inivtl, INIVTL)
(
 void
)
{
  cs_ventil_cree_listes(cs_glob_mesh,
                        cs_glob_mesh_quantities);
}


/*----------------------------------------------------------------------------
 * Marquage des ventilateurs ; affecte le numéro de ventilateur aux cellules
 * appartenant à un ventilateur, 0 sinon
 *
 * Interface Fortran :
 *
 * SUBROUTINE NUMVTL (INDIC)
 * *****************
 *
 * INTEGER INDIC(NCELET)       : --> : Numéro du ventilateur d'appartenance
 *                             :     : de la cellule (0 si hors ventilateur)
 *----------------------------------------------------------------------------*/

void CS_PROCF (numvtl, NUMVTL)
(
 cs_int_t  indic[]              /* Numéro de ventilateur d'appartenance, ou 0 */
)
{
  cs_loc_ventil_marque_cellules(cs_glob_mesh,
                                indic);
}


/*----------------------------------------------------------------------------
 * Calcul des débits à travers les ventilateurs
 *
 * Interface Fortran :
 *
 * SUBROUTINE DEBVTL (NBRVTL, FLUMAS, FLUMAB, RHO, RHOFAB, DEBENT, DEBSOR)
 * *****************
 *
 * DOUBLE PRECISION FLUMAS(*)      : <-- : Flux de masse aux faces intérieures
 * DOUBLE PRECISION FLUMAB(*)      : <-- : Flux de masse aux faces de bord
 * DOUBLE PRECISION RHO   (*)      : <-- : Densité aux cellules
 * DOUBLE PRECISION RHOFAB(*)      : <-- : Densité aux faces de bord
 * DOUBLE PRECISION DEBENT(NBRVTL) : --> : Débit entrant par ventilateur
 * DOUBLE PRECISION DEBSOR(NBRVTL) : --> : Débit sortant par ventilateur
 *----------------------------------------------------------------------------*/

void CS_PROCF (debvtl, DEBVTL)
(
 cs_real_t  flumas[],           /* <-- Flux de masse aux faces intérieures    */
 cs_real_t  flumab[],           /* <-- Flux de masse aux faces de bord        */
 cs_real_t  rho[],              /* <-- Densité aux cellules                   */
 cs_real_t  rhofab[],           /* <-- Densité aux faces de bord              */
 cs_real_t  debent[],           /* <-- Débit entrant par ventilateur          */
 cs_real_t  debsor[]            /* <-- Débit sortant par ventilateur          */
)
{
  int i;

  cs_ventil_calcul_debits(cs_glob_mesh,
                          cs_glob_mesh_quantities,
                          flumas,
                          flumab,
                          rho,
                          rhofab);

  for (i = 0 ; i < cs_glob_ventil_nbr ; i++) {
    debent[i] = cs_glob_ventil_tab[i]->debit_entrant;
    debsor[i] = cs_glob_ventil_tab[i]->debit_sortant;
  }

}


/*----------------------------------------------------------------------------
 * Calcul de la force induite par les ventilateurs (nécessite le
 * calcul préalable des débits à travers chaque ventilateur) ;
 * La force induite est ajoutée au tableau CRVXEP (qui peut contenir
 * d'autres contributions)
 *
 * Interface Fortran :
 *
 * SUBROUTINE TSVVTL (DEBENT, DEBSOR)
 * *****************
 *
 * INTEGER          IDIMTS         : <-- : Dimension associée au terme source
 *                                 :     : de vitesse (1 : X ; 2 : Y ; 3 : Z)
 * DOUBLE PRECISION CRVEXP(NCELET) : <-> : Terme source explicite de vitesse
 *----------------------------------------------------------------------------*/

void CS_PROCF (tsvvtl, TSVVTL)
(
 cs_int_t  *idimts,             /* <-- dimension associée au
                                 *     terme source de vitesse :
                                 *     0 (X), 1 (Y), ou 2 (Z)                 */
 cs_real_t  crvexp[]            /* <-- Terme source explicite de vitesse      */
)
{
  cs_ventil_calcul_force(cs_glob_mesh_quantities,
                         (*idimts) - 1,
                         crvexp);
}


/*============================================================================
 * Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Définition d'un ventilateur (qui est ajouté à ceux déjà définis)
 *----------------------------------------------------------------------------*/

void cs_ventil_definit
(
 const cs_int_t   dim_modele,       /* Dimension du modèle de ventilateur :
                                       1 : f_constante ; 2 : profil_force ;
                                       3 : profil_force + couple tangentiel */
 const cs_int_t   dim_ventil,       /* Dimension du ventilateur :
                                       2 : pseudo-2D (maillage extrudé)
                                       3 : 3D (standard) */
 const cs_real_t  coo_axe_amont[3], /* Coord. point de l'axe en face amont */
 const cs_real_t  coo_axe_aval[3],  /* Coord. point de l'axe en face aval */
 const cs_real_t  ray_ventil,       /* Rayon du ventilateur */
 const cs_real_t  ray_pales,        /* Rayon des pales */
 const cs_real_t  ray_moyeu,        /* Rayon du moyeu */
 const cs_real_t  coeff_carac[3],   /* Coefficients des termes de degré 0,
                                       1, et 2 de la caractéristique */
 const cs_real_t  couple_axial      /* Couple axial du ventilateur*/
)
{
  int  i;

  cs_ventil_t  *ventil = NULL;

  /* Définition d'un nouveau ventilateur */

  BFT_MALLOC(ventil, 1, cs_ventil_t);

  ventil->num = cs_glob_ventil_nbr + 1;

  ventil->dim_modele = dim_modele;
  ventil->dim_ventil = dim_ventil;

  for (i = 0 ; i < 3 ; i++) {
    ventil->coo_axe_amont[i] = coo_axe_amont[i];
    ventil->coo_axe_aval[i] = coo_axe_aval[i];
  }

  ventil->ray_ventil = ray_ventil;
  ventil->ray_pales  = ray_pales;
  ventil->ray_moyeu  = ray_moyeu;

  for (i = 0 ; i < 3 ; i++)
    ventil->coeff_carac[i] = coeff_carac[i];
  ventil->couple_axial = couple_axial;

  ventil->nbr_cel = 0;
  ventil->lst_cel = NULL;

  /* Calcul de la normale directrice de l'axe */

  ventil->epaisseur = 0.0;

  for (i = 0 ; i < 3 ; i++) {
    ventil->dir_axe[i] = coo_axe_aval[i] - coo_axe_amont[i];
    ventil->epaisseur += (ventil->dir_axe[i] * ventil->dir_axe[i]);
  }
  ventil->epaisseur = sqrt(ventil->epaisseur);

  for (i = 0 ; i < 3 ; i++)
    ventil->dir_axe[i] /= ventil->epaisseur;

  /* Surface initialisée à 0, sera initialisée par cs_ventil_cree_listes */

  ventil->surface = 0.0;

  /* Débits initialisés à 0 */

  ventil->debit_entrant = 0.0;
  ventil->debit_sortant = 0.0;

  /* Redimensionnement du tableau des ventilateurs si nécessaire */

  if (cs_glob_ventil_nbr == cs_glob_ventil_nbr_max) {
    cs_glob_ventil_nbr_max = (cs_glob_ventil_nbr_max + 1) * 2;
    BFT_REALLOC(cs_glob_ventil_tab, cs_glob_ventil_nbr_max, cs_ventil_t *);
  }

  /* Ajout dans le tableau des ventilateurs */

  cs_glob_ventil_tab[cs_glob_ventil_nbr] = ventil;
  cs_glob_ventil_nbr += 1;

}


/*----------------------------------------------------------------------------
 * Destruction des structures associées aux ventilateurs
 *----------------------------------------------------------------------------*/

void cs_ventil_detruit_tous
(
 void
)
{
  int i;

  cs_ventil_t  *ventil = NULL;

  for (i = 0 ; i < cs_glob_ventil_nbr ; i++) {

    ventil = cs_glob_ventil_tab[i];

    BFT_FREE(ventil->lst_cel);

    BFT_FREE(ventil);

  }

  cs_glob_ventil_nbr_max = 0;
  cs_glob_ventil_nbr = 0;
  BFT_FREE(cs_glob_ventil_tab);

}


/*----------------------------------------------------------------------------
 * Recherche des cellules appartenant aux différents ventilateurs
 *----------------------------------------------------------------------------*/

void cs_ventil_cree_listes
(
 const cs_mesh_t              *mesh,             /* <-- structure maillage associée  */
 const cs_mesh_quantities_t   *mesh_quantities   /* <-- grandeurs du maillage        */
)
{
  cs_int_t  icel, icel_1, icel_2;
  cs_int_t  ifac;
  cs_int_t  ivtl;
  cs_int_t  idim;

  cs_real_t  coo_axe;
  cs_real_t  d_2_axe;
  cs_real_t  d_cel_axe[3];
  cs_real_t  surf_loc;

  cs_ventil_t  *ventil = NULL;
  cs_int_t  *cpt_cel_vtl = NULL;
  cs_int_t  *num_vtl_cel = NULL;

  const cs_int_t  nbr_cel_et = mesh->n_cells_with_ghosts;
  const cs_int_t  *fac_cel = mesh->i_face_cells;
  const cs_int_t  *fbr_cel = mesh->b_face_cells;
  const cs_real_t  *coo_cen  = mesh_quantities->cell_cen;
  const cs_real_t  *surf_fac = mesh_quantities->i_face_normal;
  const cs_real_t  *surf_fbr = mesh_quantities->b_face_normal;

  /* Création d'un tableau de marquage des cellules */
  /*------------------------------------------------*/

  BFT_MALLOC(num_vtl_cel, nbr_cel_et, cs_int_t);

  for (icel = 0 ; icel < nbr_cel_et ; icel++)
    num_vtl_cel[icel] = 0;

  /* Boucle principale sur les cellules */

  for (icel = 0 ; icel < nbr_cel_et ; icel++) {

    /* Boucle sur les ventilateurs */

    for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {

      ventil = cs_glob_ventil_tab[ivtl];

      /* Vecteur allant du point de l'axe face amont au centre cellule */

      for (idim = 0 ; idim < 3 ; idim++) {
        d_cel_axe[idim] =   (coo_cen[icel*3 + idim])
                          - ventil->coo_axe_amont[idim];
      }

      /* Produit scalaire avec le vecteur directeur de l'axe */

      coo_axe = (  d_cel_axe[0] * ventil->dir_axe[0]
                 + d_cel_axe[1] * ventil->dir_axe[1]
                 + d_cel_axe[2] * ventil->dir_axe[2]);

      /* Cellule potentiellement dans le ventilateur si la
         projection de son centre sur l'axe est bien dans l'épaisseur */

      if (coo_axe >= 0.0 && coo_axe <= ventil->epaisseur) {

        /* Projection du vecteur allant du point de l'axe face amont au
           centre cellule dans le plan du ventilateur */

        for (idim = 0 ; idim < 3 ; idim++)
          d_cel_axe[idim] -= coo_axe * ventil->dir_axe[idim];

        /* Distance au carré à l'axe (carrés moins chers à calculer
           que racines carrés, donc on passe tout au carré) */

        d_2_axe = (  d_cel_axe[0] * d_cel_axe[0]
                   + d_cel_axe[1] * d_cel_axe[1]
                   + d_cel_axe[2] * d_cel_axe[2]);

        /* Si la cellule est dans le ventilateur */

        if (d_2_axe <= ventil->ray_ventil * ventil->ray_ventil) {

          num_vtl_cel[icel] = ivtl + 1;
          ventil->nbr_cel += 1;
          break;

        }

      }

    } /* Fin de la boucle sur les ventilateurs */

  } /* Fin de la boucle principale sur les cellules */

  /* Création des listes de cellules appartenant à chaque ventilateur */
  /*------------------------------------------------------------------*/

  BFT_MALLOC(cpt_cel_vtl, cs_glob_ventil_nbr, cs_int_t);

  for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {

    ventil = cs_glob_ventil_tab[ivtl];
    BFT_MALLOC(ventil->lst_cel, ventil->nbr_cel, cs_int_t);

    cpt_cel_vtl[ivtl] = 0;
  }

  for (icel = 0 ; icel < nbr_cel_et ; icel++) {

    if (num_vtl_cel[icel] > 0) {
      ivtl = num_vtl_cel[icel] - 1;
      ventil = cs_glob_ventil_tab[ivtl];
      ventil->lst_cel[cpt_cel_vtl[ivtl]] = icel + 1;
      cpt_cel_vtl[ivtl] += 1;
    }

  }

#if defined(DEBUG) && !defined(NDEBUG)
  for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {
    ventil = cs_glob_ventil_tab[ivtl];
    assert(cpt_cel_vtl[ivtl] == ventil->nbr_cel);
  }
#endif

  /* Calcul de la surface de chaque ventilateur */
  /*--------------------------------------------*/

  /* Contribution à l'intérieur du domaine */

  for (ifac = 0 ; ifac < mesh->n_i_faces ; ifac++) {

    icel_1 = fac_cel[ifac * 2]     - 1;
    icel_2 = fac_cel[ifac * 2 + 1] - 1;

    if (   icel_1 < mesh->n_cells /* Assure contrib. par un seul domaine */
        && num_vtl_cel[icel_1] != num_vtl_cel[icel_2]) {

      surf_loc = CS_LOC_MODULE((surf_fac + 3*ifac));
      if (num_vtl_cel[icel_1] > 0) {
        ivtl = num_vtl_cel[icel_1] - 1;
        ventil = cs_glob_ventil_tab[ivtl];
        ventil->surface += surf_loc;
      }
      if (num_vtl_cel[icel_2] > 0) {
        ivtl = num_vtl_cel[icel_2] - 1;
        ventil = cs_glob_ventil_tab[ivtl];
        ventil->surface += surf_loc;
      }
    }

  }

  /* Contribution au bord du domaine */

  for (ifac = 0 ; ifac < mesh->n_b_faces ; ifac++) {

    if (num_vtl_cel[fbr_cel[ifac] - 1] > 0) {
      surf_loc = CS_LOC_MODULE((surf_fbr + 3*ifac));
      ivtl = num_vtl_cel[fbr_cel[ifac] - 1] - 1;
      ventil = cs_glob_ventil_tab[ivtl];
      ventil->surface += surf_loc;
    }

  }

#if defined(_CS_HAVE_MPI)
  if (cs_glob_base_nbr > 1) {

    for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {
      cs_real_t surf_glob;

      surf_loc = (cs_glob_ventil_tab[ivtl])->surface;
      MPI_Allreduce (&surf_loc, &surf_glob, 1, CS_MPI_REAL, MPI_SUM,
                     cs_glob_base_mpi_comm);
      (cs_glob_ventil_tab[ivtl])->surface = surf_glob;
    }

  }
#endif

  /* Libération mémoire */

  BFT_FREE(cpt_cel_vtl);
  BFT_FREE(num_vtl_cel);
}


#if 0
/*----------------------------------------------------------------------------
 * Création d'une coupe correspondant aux faces de bord des ventilateurs
 *----------------------------------------------------------------------------*/

void cs_ventil_cree_coupe
(
 const cs_mesh_t      *mesh          /* <-- structure maillage        */
)
{
  cs_int_t   icel_1, icel_2;
  cs_int_t   ifac;

  cs_int_t  cpt_fac = 0;
  cs_int_t  cpt_fbr = 0;

  cs_int_t  *num_vtl_cel = NULL;
  cs_int_t  *liste_fac = NULL;
  cs_int_t  *liste_fbr = NULL;

  const cs_int_t  nbr_cel_et = mesh->n_cells_with_ghosts;
  const cs_int_t  nbr_fac = mesh->n_i_faces;
  const cs_int_t  nbr_fbr = mesh->n_b_faces;
  const cs_int_t  *fac_cel = mesh->i_face_cells;
  const cs_int_t  *fbr_cel = mesh->b_face_cells;

  /* Marquage des cellules et faces */

  BFT_MALLOC(num_vtl_cel, nbr_cel_et, cs_int_t);
  BFT_MALLOC(liste_fac, nbr_fac, cs_int_t);
  BFT_MALLOC(liste_fbr, nbr_fbr, cs_int_t);

  cs_loc_ventil_marque_cellules(mesh,
                                num_vtl_cel);

  /* Contribution à l'intérieur du domaine */

  for (ifac = 0 ; ifac < nbr_fac ; ifac++) {

    icel_1 = fac_cel[ifac * 2]     - 1;
    icel_2 = fac_cel[ifac * 2 + 1] - 1;

    if (num_vtl_cel[icel_1] != num_vtl_cel[icel_2])
      liste_fac[cpt_fac++] = ifac;

  }

  /* Contribution au bord du domaine */

  for (ifac = 0 ; ifac < nbr_fbr ; ifac++) {

    if (num_vtl_cel[fbr_cel[ifac] - 1] > 0)
      liste_fbr[cpt_fbr++] = ifac;

  }

  /* Libération mémoire */

  BFT_FREE(num_vtl_cel);
  BFT_FREE(liste_fac);
  BFT_FREE(liste_fbr);
}
#endif


/*----------------------------------------------------------------------------
 * Calcul des debits à travers les ventilateurs
 *----------------------------------------------------------------------------*/

void cs_ventil_calcul_debits
(
 const cs_mesh_t      *mesh,         /* <-- structure maillage        */
 const cs_mesh_quantities_t  *mesh_quantities,     /* <-- grandeurs du maillage     */
 const cs_real_t           flux_masse_fac[], /* <-- flux masse faces internes */
 const cs_real_t           flux_masse_fbr[], /* <-- flux masse faces de bord  */
 const cs_real_t           densite_cel[],    /* <-- densité aux cellules      */
 const cs_real_t           densite_fbr[]     /* <-- densité aux faces de bord */
)
{
  cs_int_t   icel, icel_1, icel_2;
  cs_int_t   ifac;
  cs_int_t   ivtl;
  cs_int_t   idim;
  cs_int_t   i, sens;

  cs_real_t  debit;
  cs_real_t  orient[3];

  cs_ventil_t  *ventil = NULL;
  cs_int_t  *num_vtl_cel = NULL;

  const cs_int_t  nbr_cel_et = mesh->n_cells_with_ghosts;
  const cs_int_t  nbr_fac = mesh->n_i_faces;
  const cs_int_t  nbr_fbr = mesh->n_b_faces;
  const cs_real_t  *coo_cen = mesh_quantities->cell_cen;
  const cs_int_t   *fac_cel = mesh->i_face_cells;
  const cs_int_t   *fbr_cel = mesh->b_face_cells;

  /* Marquage des cellules */

  BFT_MALLOC(num_vtl_cel, nbr_cel_et, cs_int_t);

  cs_loc_ventil_marque_cellules(mesh,
                                num_vtl_cel);

  /* Mise à zéro des débits par ventilateur */

  for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {
    ventil = cs_glob_ventil_tab[ivtl];
    ventil->debit_entrant = 0.0;
    ventil->debit_sortant = 0.0;
  }

  /* Contribution à l'intérieur du domaine */

  for (ifac = 0 ; ifac < nbr_fac ; ifac++) {

    icel_1 = fac_cel[ifac * 2]     - 1;
    icel_2 = fac_cel[ifac * 2 + 1] - 1;

    if (   icel_1 < mesh->n_cells /* Assure contrib. par un seul domaine */
        && num_vtl_cel[icel_1] != num_vtl_cel[icel_2]) {

      for (idim = 0 ; idim < 3 ; idim++)
        orient[idim] = coo_cen[icel_2*3 + idim] - coo_cen[icel_1*3 + idim];

      for (i = 0 ; i < 2 ; i++) {

        icel = fac_cel[ifac * 2 + i] - 1;
        ivtl = num_vtl_cel[icel] - 1;

        if (ivtl > -1) {
          ventil = cs_glob_ventil_tab[ivtl];
          debit = CS_ABS(flux_masse_fac[ifac]/densite_cel[icel]);
          sens = (i == 0 ? 1 : - 1);
          if (CS_LOC_PRODUIT_SCALAIRE(ventil->dir_axe, orient) * sens > 0.0)
            ventil->debit_sortant += debit;
          else
            ventil->debit_entrant += debit;
        }

      }

    }

  }

  /* Contribution au bord du domaine */

  for (ifac = 0 ; ifac < nbr_fbr ; ifac++) {

    ivtl = num_vtl_cel[fbr_cel[ifac] - 1] - 1;

    if (ivtl > -1) {

      ventil = cs_glob_ventil_tab[ivtl];

      for (idim = 0 ; idim < 3 ; idim++)
        orient[idim] = mesh_quantities->b_face_normal[ifac * 3 + idim];

      debit = CS_ABS(flux_masse_fbr[ifac]/densite_fbr[ifac]);
      if (CS_LOC_PRODUIT_SCALAIRE(ventil->dir_axe, orient) > 0.0)
        ventil->debit_sortant += debit;
      else
        ventil->debit_entrant += debit;

    }

  }

#if defined(_CS_HAVE_MPI)
  if (cs_glob_base_nbr > 1) {

    for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {

      cs_real_t debit_glob[2];
      cs_real_t debit_loc[2];

      ventil = cs_glob_ventil_tab[ivtl];

      debit_loc[0] = ventil->debit_sortant;
      debit_loc[1] = ventil->debit_entrant;

      MPI_Allreduce (debit_loc, debit_glob, 2, CS_MPI_REAL, MPI_SUM,
                     cs_glob_base_mpi_comm);

      ventil->debit_sortant = debit_glob[0];
      ventil->debit_entrant = debit_glob[1];

    }

  }
#endif

  /* En 2D, on ramène le débit a l'unité de longueur */

  if (ventil->dim_ventil == 2) {
    cs_real_t  surf_2d;
    surf_2d =   (0.5*ventil->surface - 2*ventil->ray_ventil*ventil->epaisseur)
              /                       (2*ventil->ray_ventil+ventil->epaisseur);
    ventil->debit_sortant = ventil->debit_sortant / surf_2d;
    ventil->debit_entrant = ventil->debit_entrant / surf_2d;
  }

  /* Libération mémoire */

  BFT_FREE(num_vtl_cel);
}


/*----------------------------------------------------------------------------
 * Calcul de la force induite par les ventilateurs (nécessite le
 * calcul préalable des débits à travers chaque ventilateur).
 * La force induite est ajoutée au tableau t_source (qui peut contenir
 * d'autres contributions)
 *----------------------------------------------------------------------------*/

void cs_ventil_calcul_force
(
 const cs_mesh_quantities_t  *mesh_quantities,  /* <-- grandeurs du maillage     */
 const cs_int_t            idim_source,         /* --> dimension associée au
                                                 *     terme source de vitesse :
                                                 *     0 (X), 1 (Y), ou 2 (Z)    */
       cs_real_t           t_source[]           /* --> terme source de vitesse   */
)
{
  cs_int_t  icel, iloc;
  cs_int_t  ivtl;
  cs_int_t  idim;

  cs_real_t  f_z, f_theta;
  cs_real_t  f_rot[3];

  const cs_real_t  *coo_cen = mesh_quantities->cell_cen;
  const cs_real_t  pi = 3.14159265358979323846;

  /* Calcul de la force induite par les ventilateurs */

  /* Boucle sur les ventilateurs */
  /*-----------------------------*/

  for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {

    const cs_ventil_t  *ventil = cs_glob_ventil_tab[ivtl];

    const cs_real_t  ray_moyeu  = ventil->ray_moyeu;
    const cs_real_t  ray_pales  = ventil->ray_pales;
    const cs_real_t  ray_ventil = ventil->ray_ventil;

    const cs_real_t  debit_moy = 0.5 * (  ventil->debit_entrant
                                        + ventil->debit_sortant);

    const cs_real_t  delta_p = - (ventil->coeff_carac[2] * debit_moy*debit_moy)
                               + (ventil->coeff_carac[1] * debit_moy)
                               + (ventil->coeff_carac[0]);

    /* Boucle sur les cellules du ventilateur */
    /*----------------------------------------*/

    for (iloc = 0 ; iloc < ventil->nbr_cel ; iloc++) {

      icel = ventil->lst_cel[iloc] - 1;

      f_z = 0.0;
      f_theta = 0.0;
      f_rot[0] = 0.0, f_rot[1] = 0.0, f_rot[2] = 0.0;

      if (ray_pales < 1.0e-12 && ray_moyeu < 1.0e-12) {

        f_z = delta_p / ventil->epaisseur;
        f_theta = 0.0;

      }
      else if (ray_moyeu < ray_pales) {

        cs_real_t  r_1, r_2, aux, aux_1, aux_2, coo_axe, d_axe, d_cel_axe[3];

        r_1 = 0.7  * ventil->ray_pales;
        r_2 = 0.85 * ventil->ray_pales;

        if (ventil->dim_ventil == 2) {
          aux_1 =   (delta_p * 2.0 * ray_ventil)
                  / (ventil->epaisseur * (1.15*ray_pales - ray_moyeu));
          aux_2 = 0.0;
        }
        else {
          cs_real_t f_base;
          const cs_real_t ray_moyeu3 = ray_moyeu * ray_moyeu * ray_moyeu;
          const cs_real_t ray_pales3 = ray_pales * ray_pales * ray_pales;
          const cs_real_t ray_pales2 = ray_pales * ray_pales;
          const cs_real_t ray_ventil2 = ray_ventil * ray_ventil;
          f_base =   (0.7*ray_pales - ray_moyeu)
                   / (1.0470*ventil->epaisseur * (  ray_moyeu3
                                                  + 1.4560*ray_pales3
                                                  - 2.570*ray_pales2*ray_moyeu));
          aux_1 = f_base * delta_p * pi * ray_ventil2;
          aux_2 = f_base * ventil->couple_axial;
        }

        /* Vecteur allant du point de l'axe face amont au centre cellule */

        for (idim = 0 ; idim < 3 ; idim++) {
          d_cel_axe[idim] =   (coo_cen[icel*3 + idim])
                            - ventil->coo_axe_amont[idim];
        }

        /* Projection du centre de la cellule sur l'axe du ventilateur */

        coo_axe = (  d_cel_axe[0] * ventil->dir_axe[0]
                   + d_cel_axe[1] * ventil->dir_axe[1]
                   + d_cel_axe[2] * ventil->dir_axe[2]);

        /* Projection du vecteur allant du point de l'axe face amont au
           centre cellule dans le plan du ventilateur */

        for (idim = 0 ; idim < 3 ; idim++)
          d_cel_axe[idim] -= coo_axe * ventil->dir_axe[idim];

        d_axe = CS_LOC_MODULE(d_cel_axe); /* Distance à l'axe */

        CS_LOC_PRODUIT_VECTORIEL(f_rot, ventil->dir_axe, d_cel_axe);

        aux = CS_LOC_MODULE(f_rot);
        for (idim = 0 ; idim < 3 ; idim++)
          f_rot[idim] /= aux;

        if (d_axe < ray_moyeu) {
          f_z     = 0.0;
          f_theta = 0.0;
        }
        else if (d_axe < r_1) {
          f_z     = aux_1 * (d_axe - ray_moyeu) / (r_1 - ray_moyeu);
          f_theta = aux_2 * (d_axe - ray_moyeu) / (r_1 - ray_moyeu);
        }
        else if (d_axe < r_2) {
          f_z     = aux_1;
          f_theta = aux_2;
        }
        else if (d_axe < ray_pales) {
          f_z     = aux_1 * (ray_pales - d_axe) / (ray_pales - r_2);
          f_theta = aux_2 * (ray_pales - d_axe) / (ray_pales - r_2);
        }
        else {
          f_z     = 0.0;
          f_theta = 0.0;
        }

      }

      t_source[icel] +=   (f_z * ventil->dir_axe[idim_source])
                        + (f_theta * f_rot[idim_source]);

    }  /* Fin de la boucle sur les cellules du ventilateur */

  } /* Fin de la boucle sur les ventilateurs */

}


/*============================================================================
 * Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Marquage des cellules appartenant aux différents ventilateurs
 * (par le numéro de ventilateur, 0 sinon)
 *----------------------------------------------------------------------------*/

static void cs_loc_ventil_marque_cellules
(
 const cs_mesh_t  *const mesh,     /* <-- structure maillage associée */
       cs_int_t          num_vtl_cel[] /* --> indicateur par cellule      */
)
{
  cs_int_t   icel;
  cs_int_t   ivtl;
  cs_int_t   iloc;

  cs_ventil_t  *ventil;

  const cs_int_t  nbr_cel_et = mesh->n_cells_with_ghosts;

  /* Marquage des cellules */
  /*-----------------------*/

  for (icel = 0 ; icel < nbr_cel_et ; icel++)
    num_vtl_cel[icel] = 0;

  for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {

    ventil = cs_glob_ventil_tab[ivtl];

    for (iloc = 0 ; iloc < ventil->nbr_cel ; iloc++) {
      icel = ventil->lst_cel[iloc] - 1;
      num_vtl_cel[icel] = ivtl + 1;
    }

  }

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
