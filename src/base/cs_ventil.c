/*============================================================================
 * Management of fans
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_base.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ventil.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/* Structure associated to a fan */

struct _cs_ventil_t {

  int                 num;              /* Fan number */
  int                 dim_modele;       /* 1D, 2D, or 3D modelling */
  int                 dim_ventil;       /* 2D or 3D geometry */

  cs_real_t           coo_axe_amont[3]; /* Axis point coordinates of the
                                           upstream face */
  cs_real_t           coo_axe_aval[3];  /* Axis point coordinates of the
                                           downstrem face */
  cs_real_t           dir_axe[3];       /* Unit vector of the axis
                                           (upstream to downstream) */
  cs_real_t           epaisseur;        /* Fan thickness */
  cs_real_t           surface;          /* Fan total surface */

  cs_real_t           ray_ventil;       /* Fan radius */
  cs_real_t           ray_pales;        /* Blades radius */
  cs_real_t           ray_moyeu;        /* Hub radius */
  cs_real_t           coeff_carac[3];   /* Coefficients of the terms of
                                           degree 0, 1 and 2 of the
                                           characteristic curve */
  cs_real_t           couple_axial;     /* Fan axial couple */

  cs_int_t            nbr_cel;          /* Number of cells */

  cs_int_t           *lst_cel;          /* List of the cells belonging
                                           to the fan */

  cs_real_t           debit_entrant;    /* Current inlet flow */
  cs_real_t           debit_sortant;    /* Current outlet flow */

} ;


/*============================================================================
 * Global variables
 *============================================================================*/

/* Fans array */

cs_int_t         cs_glob_ventil_nbr_max = 0;

cs_int_t         cs_glob_ventil_nbr = 0;
cs_ventil_t  * * cs_glob_ventil_tab = NULL;


/*============================================================================
 * Macro definitions
 *============================================================================*/

enum {X, Y, Z} ;

#define CS_LOC_VECTORIAL_PRODUCT(prod_vect, vect1, vect2)  \
  (prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z], \
   prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z], \
   prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y])

#define CS_LOC_DOT_PRODUCT(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z])

#define CS_LOC_MODULE(vect) \
  sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])


/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Mark the cells belonging to the different fans
 * (by the fan number, 0 otherwise)
 *
 * parameters:
 *   mesh        <-- associated mesh structure
 *   num_vtl_cel --> indicator by cell
 *----------------------------------------------------------------------------*/

static void
cs_loc_ventil_marque_cellules(const cs_mesh_t  *const mesh,
                              cs_int_t          num_vtl_cel[])
{
  cs_int_t   icel;
  cs_int_t   ivtl;
  cs_int_t   iloc;

  cs_ventil_t  *ventil;

  const cs_int_t  nbr_cel_et = mesh->n_cells_with_ghosts;

  /* Mark the cells */

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

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Get the number of fans.
 *
 * Fortran interface:
 *
 * SUBROUTINE TSTVTL
 * *****************
 *
 * INTEGER          NBRVTL         : --> : number of fans
 *----------------------------------------------------------------------------*/

void CS_PROCF (tstvtl, TSTVTL)
(
 cs_int_t  *const nbrvtl
)
{
  *nbrvtl = cs_glob_ventil_nbr;
}

/*----------------------------------------------------------------------------
 * Adds a fan.
 *
 * Fortran interface:
 *
 * SUBROUTINE DEFVTL
 * *****************
 *
 * INTEGER          DIMMOD     : <-- : Fan model dimension:
 *                             :     : 1: constant_f; 2: force_profile;
 *                             :     : 3: force_profile + tangential couple
 *                  DIMVTL     : <-- : Fan dimension:
 *                             :     : 2: pseudo-2D (extruded mesh)
 *                             :     : 3: 3D (standard)
 * DOUBLE PRECISION XYZVT1(3)  : <-- : Coo. of the axis point in upstream face
 * DOUBLE PRECISION XYZVT2(3)  : <-- : Coo. of the axis point in downstream face
 * DOUBLE PRECISION RVVT       : <-- : Fan radius
 * DOUBLE PRECISION RPVT       : <-- : Blades radius
 * DOUBLE PRECISION RMVT       : <-- : Hub radius
 * DOUBLE PRECISION CCARAC(3)  : <-- : Coefficients of degre 0, 1 and 2
 *                             :     : of the characteristic curve
 * DOUBLE PRECISION TAUVT      : <-- : Fan axial couple
 *----------------------------------------------------------------------------*/

void CS_PROCF (defvtl, DEFVTL)
(
 const cs_int_t   *const dimmod,
 const cs_int_t   *const dimvtl,
 const cs_real_t         xyzvt1[3],
 const cs_real_t         xyzvt2[3],
 const cs_real_t  *const rvvt,
 const cs_real_t  *const rpvt,
 const cs_real_t  *const rmvt,
 const cs_real_t         ccarac[3],
 const cs_real_t  *const tauvt
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
 * Build the list of cells associated to the fans
 *
 * Fotrtran interface:
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
 * Mark the fans and associate the fan number to the cells belonging to
 * thus fan, 0 otherwise.
 *
 * Fortran interface:
 *
 * SUBROUTINE NUMVTL (INDIC)
 * *****************
 *
 * INTEGER INDIC(NCELET)       : --> : Fan number (0 if outside the fan)
 *----------------------------------------------------------------------------*/

void CS_PROCF (numvtl, NUMVTL)
(
 cs_int_t  indic[]
)
{
  cs_loc_ventil_marque_cellules(cs_glob_mesh,
                                indic);
}

/*----------------------------------------------------------------------------
 * Calculate the flows through the fans
 *
 * Fortran interface:
 *
 * subroutine debvtl
 * *****************
 *
 * parameters:
 *  flumas         <-- Interior faces mass flux
 *  flumab         <-- Boundary faces mass flux
 *  rho            <-- Density at cells
 *  rhofab         <-- Density at boundary faces
 *  debent         --> Inlet flow through the fan (size is nbrvtl)
 *  debsor         --> Outlet flow through the fan (size is nbrvtl)
 *----------------------------------------------------------------------------*/

void CS_PROCF (debvtl, DEBVTL)
(
 cs_real_t  flumas[],
 cs_real_t  flumab[],
 cs_real_t  rho[],
 cs_real_t  rhofab[],
 cs_real_t  debent[],
 cs_real_t  debsor[]
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
 * Calculate the force induced by the fans (needs a previous calculation
 * of the flows through each fan).
 *
 * The induced force is added to the array crvxep (which can have other
 * contributions).
 *
 * Fortran interface:
 *
 * subroutine tsvvtl
 * *****************
 *
 * parameters:
 *  idimts         <-- Dimension associated to the source
 *                     term of velocity (1: X; 2: Y; 3: Z)
 *  crvexp         <-> Explicit source term (velocity)
 *----------------------------------------------------------------------------*/

void CS_PROCF (tsvvtl, TSVVTL)
(
 cs_int_t  *idimts,
 cs_real_t  crvexp[]
)
{
  cs_ventil_calcul_force(cs_glob_mesh_quantities,
                         (*idimts) - 1,
                         crvexp);
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Fan definition (added to the ones previously defined)
 *
 * parameters:
 *   dim_modele    <-- Fan model dimension:
 *                     1: constant_f
 *                     2: force_profile
 *                     3: force_profile + tangential couple
 *   dim_ventil    <-- Fan dimension:
 *                     2: pseudo-2D (extruded mesh)
 *                     3: 3D (standard)
 *   coo_axe_amont <-- Coo. of the axis point in upstream face
 *   coo_axe_aval  <-- Coo. of the axis point in downstream face
 *   ray_ventil    <-- Fan radius
 *   ray_pales     <-- Blades radius
 *   ray_moyeu     <-- Hub radius
 *   coeff_carac   <-- Coefficients of degre 0, 1 and 2 of
                       the characteristic curve
 *   couple_axial  <-- Fan axial couple
 *----------------------------------------------------------------------------*/

void
cs_ventil_definit(const cs_int_t   dim_modele,
                  const cs_int_t   dim_ventil,
                  const cs_real_t  coo_axe_amont[3],
                  const cs_real_t  coo_axe_aval[3],
                  const cs_real_t  ray_ventil,
                  const cs_real_t  ray_pales,
                  const cs_real_t  ray_moyeu,
                  const cs_real_t  coeff_carac[3],
                  const cs_real_t  couple_axial)
{
  int  i;

  cs_ventil_t  *ventil = NULL;

  /* Define a new fan */

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

  /* Compute the axis vector */

  ventil->epaisseur = 0.0;

  for (i = 0 ; i < 3 ; i++) {
    ventil->dir_axe[i] = coo_axe_aval[i] - coo_axe_amont[i];
    ventil->epaisseur += (ventil->dir_axe[i] * ventil->dir_axe[i]);
  }
  ventil->epaisseur = sqrt(ventil->epaisseur);

  for (i = 0 ; i < 3 ; i++)
    ventil->dir_axe[i] /= ventil->epaisseur;

  /* Surface initialized to 0, will be set by cs_ventil_cree_listes */

  ventil->surface = 0.0;

  /* Flows initialized to 0 */

  ventil->debit_entrant = 0.0;
  ventil->debit_sortant = 0.0;

  /* Increase the fans array if necessary */

  if (cs_glob_ventil_nbr == cs_glob_ventil_nbr_max) {
    cs_glob_ventil_nbr_max = (cs_glob_ventil_nbr_max + 1) * 2;
    BFT_REALLOC(cs_glob_ventil_tab, cs_glob_ventil_nbr_max, cs_ventil_t *);
  }

  /* Adds in the fans array */

  cs_glob_ventil_tab[cs_glob_ventil_nbr] = ventil;
  cs_glob_ventil_nbr += 1;

}

/*----------------------------------------------------------------------------
 * Destroy the structures associated to fans
 *----------------------------------------------------------------------------*/

void
cs_ventil_detruit_tous(void)
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
 * Looks for the cells belonging to the different fans.
 *
 * parameters:
 *   mesh            <-- associated mesh structure
 *   mesh_quantities <-- mesh quantities
 *----------------------------------------------------------------------------*/

void
cs_ventil_cree_listes(const cs_mesh_t              *mesh,
                      const cs_mesh_quantities_t   *mesh_quantities)
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

  /* Create an array for cells marking */
  /*-----------------------------------*/

  BFT_MALLOC(num_vtl_cel, nbr_cel_et, cs_int_t);

  for (icel = 0 ; icel < nbr_cel_et ; icel++)
    num_vtl_cel[icel] = 0;

  /* Main loop on cells */

  for (icel = 0 ; icel < nbr_cel_et ; icel++) {

    /* Loop on fans */

    for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {

      ventil = cs_glob_ventil_tab[ivtl];

      /* Vector from the downstream face axis point to the cell centre */

      for (idim = 0 ; idim < 3 ; idim++) {
        d_cel_axe[idim] =   (coo_cen[icel*3 + idim])
                          - ventil->coo_axe_amont[idim];
      }

      /* Dot product with the axis vector */

      coo_axe = (  d_cel_axe[0] * ventil->dir_axe[0]
                 + d_cel_axe[1] * ventil->dir_axe[1]
                 + d_cel_axe[2] * ventil->dir_axe[2]);

      /* Cell potentially in the fan if its centre projection on the axis
         is within the thickness */

      if (coo_axe >= 0.0 && coo_axe <= ventil->epaisseur) {

        /* Projection of the vector from the downstream face axis point
           to the cell centre in the fan plane */

        for (idim = 0 ; idim < 3 ; idim++)
          d_cel_axe[idim] -= coo_axe * ventil->dir_axe[idim];

        /* Square distance to the axis */

        d_2_axe = (  d_cel_axe[0] * d_cel_axe[0]
                   + d_cel_axe[1] * d_cel_axe[1]
                   + d_cel_axe[2] * d_cel_axe[2]);

        /* If the cell is in the fan */

        if (d_2_axe <= ventil->ray_ventil * ventil->ray_ventil) {

          num_vtl_cel[icel] = ivtl + 1;
          ventil->nbr_cel += 1;
          break;

        }

      }

    } /* End of loop on fans */

  } /* End of main loop on cells */

  /* Create the lists of cells belonging to each fan */
  /*-------------------------------------------------*/

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

  /* Compute each fan surface */
  /*--------------------------*/

  /* Contribution to the domain interior */

  for (ifac = 0 ; ifac < mesh->n_i_faces ; ifac++) {

    icel_1 = fac_cel[ifac * 2]     - 1;
    icel_2 = fac_cel[ifac * 2 + 1] - 1;

    if (   icel_1 < mesh->n_cells /* Make sure the contrib is from one domain */
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

  /* Contribution to the domain boundary */

  for (ifac = 0 ; ifac < mesh->n_b_faces ; ifac++) {

    if (num_vtl_cel[fbr_cel[ifac] - 1] > 0) {
      surf_loc = CS_LOC_MODULE((surf_fbr + 3*ifac));
      ivtl = num_vtl_cel[fbr_cel[ifac] - 1] - 1;
      ventil = cs_glob_ventil_tab[ivtl];
      ventil->surface += surf_loc;
    }

  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {
      cs_real_t surf_glob;

      surf_loc = (cs_glob_ventil_tab[ivtl])->surface;
      MPI_Allreduce (&surf_loc, &surf_glob, 1, CS_MPI_REAL, MPI_SUM,
                     cs_glob_mpi_comm);
      (cs_glob_ventil_tab[ivtl])->surface = surf_glob;
    }

  }
#endif

  /* Free memory */

  BFT_FREE(cpt_cel_vtl);
  BFT_FREE(num_vtl_cel);
}


#if 0
/*----------------------------------------------------------------------------
 * Creates a post-processing mesh corresponding to the fans boundary faces
 *
 * parameters:
 *   mesh <-- mesh structure
 *----------------------------------------------------------------------------*/

void
cs_ventil_cree_coupe(const cs_mesh_t  *mesh)
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

  /* Mark cells and faces */

  BFT_MALLOC(num_vtl_cel, nbr_cel_et, cs_int_t);
  BFT_MALLOC(liste_fac, nbr_fac, cs_int_t);
  BFT_MALLOC(liste_fbr, nbr_fbr, cs_int_t);

  cs_loc_ventil_marque_cellules(mesh,
                                num_vtl_cel);

  /* Contribution to the domain interior */

  for (ifac = 0 ; ifac < nbr_fac ; ifac++) {

    icel_1 = fac_cel[ifac * 2]     - 1;
    icel_2 = fac_cel[ifac * 2 + 1] - 1;

    if (num_vtl_cel[icel_1] != num_vtl_cel[icel_2])
      liste_fac[cpt_fac++] = ifac;

  }

  /* Contribution to the domain boundary */

  for (ifac = 0 ; ifac < nbr_fbr ; ifac++) {

    if (num_vtl_cel[fbr_cel[ifac] - 1] > 0)
      liste_fbr[cpt_fbr++] = ifac;

  }

  /* Free memory */

  BFT_FREE(num_vtl_cel);
  BFT_FREE(liste_fac);
  BFT_FREE(liste_fbr);
}
#endif

/*----------------------------------------------------------------------------
 * Calculate the flows through the fans
 *
 * parameters:
 *   mesh           <-- mesh structure
 *   mesh_qantities <-- mesh quantities
 *   flux_masse_fac <-- interior faces mass flux
 *   flux_masse_fbr <-- boundary faces mass flux
 *   densite_cel    <-- density at cells
 *   densite_fbr    <-- density at boundary faces
 *----------------------------------------------------------------------------*/

void
cs_ventil_calcul_debits(const cs_mesh_t             *mesh,
                        const cs_mesh_quantities_t  *mesh_quantities,
                        const cs_real_t              flux_masse_fac[],
                        const cs_real_t              flux_masse_fbr[],
                        const cs_real_t              densite_cel[],
                        const cs_real_t              densite_fbr[])
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

  /* Mark the cells */

  BFT_MALLOC(num_vtl_cel, nbr_cel_et, cs_int_t);

  cs_loc_ventil_marque_cellules(mesh,
                                num_vtl_cel);

  /* Set the fans flows to zero */

  for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {
    ventil = cs_glob_ventil_tab[ivtl];
    ventil->debit_entrant = 0.0;
    ventil->debit_sortant = 0.0;
  }

  /* Contribution to the domain interior */

  for (ifac = 0 ; ifac < nbr_fac ; ifac++) {

    icel_1 = fac_cel[ifac * 2]     - 1;
    icel_2 = fac_cel[ifac * 2 + 1] - 1;

    if (   icel_1 < mesh->n_cells /* Make sure the contrib is from one domain */
        && num_vtl_cel[icel_1] != num_vtl_cel[icel_2]) {

      for (idim = 0 ; idim < 3 ; idim++)
        orient[idim] = coo_cen[icel_2*3 + idim] - coo_cen[icel_1*3 + idim];

      for (i = 0 ; i < 2 ; i++) {

        icel = fac_cel[ifac * 2 + i] - 1;
        ivtl = num_vtl_cel[icel] - 1;

        if (ivtl > -1) {
          ventil = cs_glob_ventil_tab[ivtl];
          debit = flux_masse_fac[ifac]/densite_cel[icel];
          sens = (i == 0 ? 1 : - 1);
          if (CS_LOC_DOT_PRODUCT(ventil->dir_axe, orient) * sens > 0.0)
            ventil->debit_sortant += debit;
          else
            ventil->debit_entrant += debit;
        }

      }

    }

  }

  /* Contribution to the domain boundary */

  for (ifac = 0 ; ifac < nbr_fbr ; ifac++) {

    ivtl = num_vtl_cel[fbr_cel[ifac] - 1] - 1;

    if (ivtl > -1) {

      ventil = cs_glob_ventil_tab[ivtl];

      for (idim = 0 ; idim < 3 ; idim++)
        orient[idim] = mesh_quantities->b_face_normal[ifac * 3 + idim];

      debit = flux_masse_fbr[ifac]/densite_fbr[ifac];
      if (CS_LOC_DOT_PRODUCT(ventil->dir_axe, orient) > 0.0)
        ventil->debit_sortant += debit;
      else
        ventil->debit_entrant += debit;

    }

  }

#if defined(HAVE_MPI)
  if (cs_glob_n_ranks > 1) {

    for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {

      cs_real_t debit_glob[2];
      cs_real_t debit_loc[2];

      ventil = cs_glob_ventil_tab[ivtl];

      debit_loc[0] = ventil->debit_sortant;
      debit_loc[1] = ventil->debit_entrant;

      MPI_Allreduce (debit_loc, debit_glob, 2, CS_MPI_REAL, MPI_SUM,
                     cs_glob_mpi_comm);

      ventil->debit_sortant = debit_glob[0];
      ventil->debit_entrant = debit_glob[1];

    }

  }
#endif

  /* In 2D, the flow is normalized by the surface */

  if (ventil->dim_ventil == 2) {
    cs_real_t  surf_2d;
    surf_2d =   (0.5*ventil->surface - 2*ventil->ray_ventil*ventil->epaisseur)
              /                       (2*ventil->ray_ventil+ventil->epaisseur);
    ventil->debit_sortant = ventil->debit_sortant / surf_2d;
    ventil->debit_entrant = ventil->debit_entrant / surf_2d;
  }

  /* Free memory */

  BFT_FREE(num_vtl_cel);
}

/*----------------------------------------------------------------------------
 * Calculate the force induced by the fans (needs a previous calculation
 * of the flows through each fan).
 *
 * The induced force is added to the array CRVXEP (which can have other
 * other contributions).
 *
 * parameters:
 *   mesh_quantities <-- mesh quantities
 *   idim_source     <-- Dimension associated to the source term of velocity
 *                       (1: X; 2: Y; 3: Z)
 *   t_source        <-> Explicit source term for the velocity
 *----------------------------------------------------------------------------*/

void
cs_ventil_calcul_force(const cs_mesh_quantities_t  *mesh_quantities,
                       const cs_int_t               idim_source,
                       cs_real_t                    t_source[])
{
  cs_int_t  icel, iloc;
  cs_int_t  ivtl;
  cs_int_t  idim;

  cs_real_t  f_z, f_theta;
  cs_real_t  f_rot[3];

  const cs_real_t  *coo_cen = mesh_quantities->cell_cen;
  const cs_real_t  pi = 3.14159265358979323846;

  /* Compute the force induced by fans */

  /* Loop on fans */
  /*--------------*/

  for (ivtl = 0 ; ivtl < cs_glob_ventil_nbr ; ivtl++) {

    const cs_ventil_t  *ventil = cs_glob_ventil_tab[ivtl];

    const cs_real_t  ray_moyeu  = ventil->ray_moyeu;
    const cs_real_t  ray_pales  = ventil->ray_pales;
    const cs_real_t  ray_ventil = ventil->ray_ventil;

    const cs_real_t  debit_moy = 0.5 * (  ventil->debit_sortant
                                        - ventil->debit_entrant);

    const cs_real_t  delta_p = - (ventil->coeff_carac[2] * debit_moy*debit_moy)
                               + (ventil->coeff_carac[1] * debit_moy)
                               + (ventil->coeff_carac[0]);

    /* Loop on fan cells */
    /*-------------------*/

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

        /* Vector from the downstream face axis point to the cell centre */

        for (idim = 0 ; idim < 3 ; idim++) {
          d_cel_axe[idim] =   (coo_cen[icel*3 + idim])
                            - ventil->coo_axe_amont[idim];
        }

        /* Projection of the cell centre on the fan axis */

        coo_axe = (  d_cel_axe[0] * ventil->dir_axe[0]
                   + d_cel_axe[1] * ventil->dir_axe[1]
                   + d_cel_axe[2] * ventil->dir_axe[2]);

        /* Projection of the vector from the downstream face axis point
           to the cell centre in the fan plane */

        for (idim = 0 ; idim < 3 ; idim++)
          d_cel_axe[idim] -= coo_axe * ventil->dir_axe[idim];

        d_axe = CS_LOC_MODULE(d_cel_axe); /* Distance to the axis */

        CS_LOC_VECTORIAL_PRODUCT(f_rot, ventil->dir_axe, d_cel_axe);

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

    }  /* End of loop on fan cells */

  } /* End of loop on fans */

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
