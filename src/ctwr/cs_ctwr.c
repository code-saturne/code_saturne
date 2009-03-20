/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
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
 * Definitions, Global variables variables, and functions associated with the
 * exchange zones
 *============================================================================*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif


/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_error.h>
#include <bft_printf.h>
#include <bft_mem.h>
#include <bft_file.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_interface.h>
#include <fvm_nodal_extract.h>
#include <fvm_nodal_extrude.h>
#include <fvm_nodal_project.h>
#include <fvm_writer.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"
#include "cs_ctwr_air_props.h"
#include "cs_ctwr_halo.h"
#include "cs_halo.h"
#include "cs_mesh_connect.h"
#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

enum {X, Y, Z} ;

#define CS_LOC_PRODUIT_SCALAIRE(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z])

#define CS_LOC_MODULE(vect) \
  sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])

#define CS_PRODUIT_VECTORIEL(prod_vect, vect1, vect2)  \
  (prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z], \
   prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z], \
   prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y])

#define CS_CT_MPI_TAG    (int)('C'+'S'+'Z'+'E') /* MPI tag for FVM operations */


/*============================================================================
 * Static global variables
 *============================================================================*/

static double _epsilon_denom = 1.e-14; /* Minimum denominator */

/* array of exchanges area */

cs_int_t            cs_glob_ct_nbr_max = 0;
static cs_int_t     cs_ctwr_nmaxvoi       = 50;
cs_int_t            cs_glob_ct_nbr = 0;
cs_ctwr_zone_t     ** cs_glob_ct_tab = NULL;

/* Start and end (negative) numbers associated with
   dedicated post processing meshes */

static int  cs_glob_ct_post_mesh_ext[2] = {0, 1};

/* array containing the stacking of the exchange area*/
cs_int_t  *  cs_stack_ct    = NULL;

/*array containing the treatment order of the exchanges areas */
cs_int_t  *  cs_chain_ct = NULL;



#if defined(HAVE_MPI)

MPI_Status status;

#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 *  Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Récupération du nombre de zones d'echanges
 *
 * Interface Fortran :
 *
 * SUBROUTINE NBZECT
 * *****************
 *
 * INTEGER          NBRct         : --> : nombre de zones d'echange
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbzect, NBZECT)
(
  cs_int_t  *const nbrct             /* <-- nombre de zones d'echanges        */
)
{
  *nbrct = cs_glob_ct_nbr ;
}

/*----------------------------------------------------------------------------
 * Récupération du modele de Poppe ou de Merkel
 *
 * Interface Fortran :
 *
 * SUBROUTINE AEMODEL
 * ******************
 *
 * INTEGER          IMctCH         : --> : type de ct (Poppe ou Merkel)
 *----------------------------------------------------------------------------*/

void CS_PROCF (aemode, AEMODE)
(
 cs_int_t  *const  imctch              /* <-- type de ct (Poppe ou Merkel)  */
)
{
  cs_ctwr_zone_t   *ct = cs_glob_ct_tab[0];

  *imctch = ct->imctch;
}

/*----------------------------------------------------------------------------
 * Addition d'une constante au vecteur de temperature pour toutes les zones
 * d'echange
 *
 * Interface Fortran :
 *
 * SUBROUTINE AEPROTP
 * ******************
 *
 * INTEGER          IMctCH         : --> : type de ct (Poppe ou Merkel)
 *----------------------------------------------------------------------------*/

void CS_PROCF (aeprot, AEPROT)
(
 cs_real_t  *const   cons             /* <-- constante */
)
{
  cs_int_t ct_id, i, j, ii;
  cs_ctwr_zone_t  *ct;

  for (ct_id = 0; ct_id < cs_glob_ct_nbr; ct_id++) {

    ct = cs_glob_ct_tab[ct_id];

    for (i = 0; i < ct->nnpsct; i++)
      for (j = 0; j < ct->nelect; j++) {
          ii = i*ct->nelect + j;
          ct->teau[ii] += (*cons);
    }

  }
}

/*----------------------------------------------------------------------------
 * Ajout d'une zone d'echange
 *
 * Interface Fortran :
 *
 * SUBROUTINE DEFct
 *
 * INTEGER          ICOUL       : <-- : Couleur des elements de la zone d'echange
 * INTEGER          IMctCH      : <-- :
 * INTEGER          NTYPct      : <-- :
 * INTEGER          NELEct      : <-- :
 * REAL             XAP         : <-- :
 * REAL             XNP         : <-- :
 * REAL             SURFACE     : <-- : Surface de la face superieure de la ct
 *----------------------------------------------------------------------------*/

void CS_PROCF (defct, DEFCT)
(
  const cs_int_t   *const idimct,   /* Dimemsion du probleme 2:2D  3:3D       */
  const cs_int_t   *const icoul,    /* Couleur des elements de la ct          */
  const cs_int_t   *const imctch,   /* 1: Modele de Poppe
                                       2: Merkel 0: Rien                      */
  const cs_int_t   *const ntypct,   /* 1: Contre courant  2: Courant croises
                                       3: Zone de pluie                      */
  const cs_int_t   *const nelect,   /* Nombre d'elements sur chaque ligne du
                                       maillage eau pour la zone de noeuds par
                                       segment eau */
  const cs_real_t  *const deltat,   /* Ecart de temperature imposé en entree
                                       de la zone d'echange */
  const cs_real_t  *const teau,     /* Teau en entree de la zone d'echange    */
  const cs_real_t  *const fem,      /* fem en entree de la zone d'echange     */
  const cs_real_t  *const xap,      /* coefficient lambda de la loi d'echange */
  const cs_real_t  *const xnp,      /* exposant n de la loi d'echange         */

  const cs_real_t  *const surface,  /* Surface totale arrivee d eau de la ct  */

  const cs_real_t   *const   cpa,   /* Capacite calorifique de l air          */
  const cs_real_t   *const   cpv,   /* Capacite calorifique de la vapeur      */
  const cs_real_t   *const   cpe,   /* Capacite calorifique de l eau          */
  const cs_real_t   *const   hv0,   /* Chaleur latente                        */
  const cs_real_t   *const   rhoe,  /* Masse volumique de l eau               */
  const cs_real_t   *const   dgout, /* Diametre de goutte pour
                                       les zones de pluie                     */
  const cs_real_t   *const   visc,  /* Viscosite Dynamique                    */
  const cs_real_t   *const   conduc /* Conductivite                           */
)
{
  cs_ctwr_definit(*idimct,*icoul,*imctch,*ntypct,*nelect,*deltat,*teau,*fem,*xap,*xnp,*surface,
                *cpa,*cpv,*cpe,*hv0,*rhoe,*dgout,*visc,*conduc);
}


/*----------------------------------------------------------------------------
 * Resolution des variables eau
 *
 * Interface Fortran :
 *
 * SUBROUTINE AETEAU ( )
 *----------------------------------------------------------------------------*/

void CS_PROCF (aeteau, AETEAU)
(
  cs_real_t   temp[],         /* Temperature air */
  cs_real_t   xa[],           /* humidite air */
  cs_real_t   rho[],          /* masse volumique air */
  cs_real_t   vitx[],         /* vitesse air suivant x */
  cs_real_t   vity[],         /* vitesse air suivant y */
  cs_real_t   vitz[],         /* vitesse air suivant z */
  const cs_real_t   *const gx,      /* composante x de la gravite */
  const cs_real_t   *const gy,      /* composante y de la gravite */
  const cs_real_t   *const gz       /* composante z de la gravite */

)
{
   cs_ctwr_aeteau(temp,xa,rho,vitx,vity,vitz,*gx,*gy,*gz);
}


/*----------------------------------------------------------------------------
 * Calcul des termes sources pour l'air
 *
 * Interface Fortran :
 *
 * SUBROUTINE AETSSC ( )
 *----------------------------------------------------------------------------*/

void CS_PROCF (aetssc, AETSSC)
(
  const cs_int_t    *const iscal,       /*   */

  cs_real_t   temp[],             /* Temperature air */
  cs_real_t   xa[],               /* humidite air */
  cs_real_t   rho[],              /* masse volumique air */
  cs_real_t   utsim[],                  /* vitesse verticale air */
  cs_real_t   utsex[],                  /* vitesse horizontale air */
  cs_real_t   vitx[],             /* vitesse air suivant x */
  cs_real_t   vity[],             /* vitesse air suivant y */
  cs_real_t   vitz[],             /* vitesse air suivant z */
  const cs_real_t   *const gx,          /* composante x de la gravite */
  const cs_real_t   *const gy,          /* composante y de la gravite */
  const cs_real_t   *const gz           /* composante z de la gravite */
)
{
  cs_ctwr_aetssc(*iscal, temp,xa,rho,utsim,utsex,vitx,vity,vitz,*gx,*gy,*gz);
}




/*----------------------------------------------------------------------------
 * Calcul des PdC induites dans les zones de pluie
 *
 * Interface Fortran :
 *
 * SUBROUTINE AETSVI ( )
 *----------------------------------------------------------------------------*/

void CS_PROCF (aetsvi, AETSVI)
(
  const cs_int_t    *const idim,
  const cs_real_t   rho[],              /* masse volumique air    */
  const cs_real_t   vitx[],             /* vitesse air suivant x  */
  const cs_real_t   vity[],             /* vitesse air suivant y  */
  const cs_real_t   vitz[],             /* vitesse air suivant z  */
  const cs_real_t   xair[],             /* humidite de l'air      */
  const cs_real_t   *const gx,          /*                        */
  const cs_real_t   *const gy,          /*                        */
  const cs_real_t   *const gz,          /*                        */
  cs_real_t   utsex[]                   /* terme source explicite */

)
{
  cs_ctwr_aetsvi(*idim,rho,vitx,vity,vitz,xair,*gx,*gy,*gz,utsex);

}



/*----------------------------------------------------------------------------
 * Bilan dans les ct
 *
 * Interface Fortran :
 *
 * SUBROUTINE BILANct ( )
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilanct, BILANCT)
(
  const cs_real_t   *const time,
  cs_real_t   fem_entree[],           /* debit eau entree           */
  cs_real_t   fem_sortie[],           /* debit eau sortie           */
  cs_real_t   teau_entree[],          /* temperature eau entree     */
  cs_real_t   teau_sortie[],          /* temperature eau sortie     */
  cs_real_t   heau_entree[],          /* enthalpie eau entree       */
  cs_real_t   heau_sortie[],          /* enthalpie eau sortie       */
  cs_real_t   tair_entree[],          /* temperature air entree     */
  cs_real_t   tair_sortie[],          /* temperature air sortie     */
  cs_real_t   xair_entree[],          /*                            */
  cs_real_t   xair_sortie[],          /*                            */
  cs_real_t   hair_entree[],          /*                            */
  cs_real_t   hair_sortie[],          /*                            */
  cs_real_t   debit_entree[],         /*                            */
  cs_real_t   debit_sortie[],         /*                            */

  const cs_real_t   temp[],           /* Temperature air            */
  const cs_real_t   xa[],             /* humidite air               */
  const cs_real_t   flux_masse_fac[], /* vitesse verticale air      */
  const cs_real_t   flux_masse_fbr[], /* vitesse horizontale air    */
  const cs_real_t   vitx[],           /* vitesse air suivant x      */
  const cs_real_t   vity[],           /* vitesse air suivant y      */
  const cs_real_t   vitz[]            /* vitesse air suivant z      */
)
{
  cs_ctwr_bilanct(*time,fem_entree,fem_sortie,teau_entree,teau_sortie,
                  heau_entree,heau_sortie,tair_entree,tair_sortie,
                  xair_entree,xair_sortie,hair_entree,hair_sortie,
                  debit_entree,debit_sortie,
                  temp,xa,flux_masse_fac,flux_masse_fbr,vitx,vity,vitz,
                  cs_glob_mesh, cs_glob_mesh_quantities);

}

/*----------------------------------------------------------------------------
 * Initialict post processing.
 *
 * Fortran Interface:
 *
 * SUBROUTINE PSTIct
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(pstict, PSTICT)
(
 void
)
{
  cs_int_t ct_id;

  for (ct_id = 0; ct_id < cs_glob_ct_nbr; ct_id++)
    cs_ctwr_post_init(ct_id, -1);
}

/*----------------------------------------------------------------------------
 * Get the local (negative) numbers associated with the first and last
 * post processing meshes dedicated to exchange area
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTEct
 * *****************
 *
 * INTEGER          first_id        : <-- : id of first post processing mesh
 * INTEGER          last_id         : <-- : id of last post processing mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstect, PSTECT)
(
 cs_int_t  *const first_id,
 cs_int_t  *const last_id
)
{
  cs_ctwr_post_id_extents(first_id, last_id);
}

/*----------------------------------------------------------------------------
 * Create nodal coupled mesh.
 * Send vertices's coordinates and connectivity of coupled mesh.
 *
 * Fortran Interface:
 *
 * SUBROUTINE GEOct
 * *****************
 *
 * INTEGER          n_ct     : <-- : number of exchange area
 *----------------------------------------------------------------------------*/

void CS_PROCF(geoct, GEOCT)
(
  const cs_real_t  *const gx,           /* composante x de la gravite */
  const cs_real_t  *const gy,           /* composante y de la gravite */
  const cs_real_t  *const gz            /* composante z de la gravite */
)
{
  /* construction du maillage eau*/
  cs_ctwr_maille(*gx, *gy, *gz, cs_glob_mesh, cs_glob_mesh_quantities );

  /* chainage des ct*/
  cs_ctwr_stacking(*gx, *gy, *gz);
  /* construction de l'interpolation  AIR -> EAU    */
  cs_ctwr_adeau(cs_glob_mesh, cs_glob_mesh_quantities);
  /* construction de l'interpolation  EAU -> AIR   */
  cs_ctwr_adair(*gx, *gy, *gz);

}


/*============================================================================
 * Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Définition d'une zone d'echange (qui est ajoutée à celles déjà définies)
 *----------------------------------------------------------------------------*/

void cs_ctwr_definit
(
  const cs_int_t   idimct,    /* Dimemsion du probleme 2:2D  3:3D */
  const cs_int_t   icoul,     /* Couleur des elements de la ct */
  const cs_int_t   imctch,    /* 1: Modele de Poppe
                                 2: Merkel
                                 0: Rien */
  const cs_int_t   ntypct,    /* 1: Contre courant
                                 2: Courant croises
                                 3: Zone de pluie */
  const cs_int_t   nelect,    /* Nombre d'elements sur chaque ligne du maillage
                                 eau pour la zone de noeuds par segment eau */
  const cs_real_t  deltat,    /* Ecart de temperature imposé en entree de la
                                 zone d'echange */
  const cs_real_t  teau_cl,   /* Teau en entree de la zone d'echange */
  const cs_real_t  fem_cl,    /* debit en entree de la zone d'echange */
  const cs_real_t  xap,       /* coefficient lambda de la loi d'echange */
  const cs_real_t  xnp,       /* exposant n de la loi d'echange */
  const cs_real_t  surface,   /* Surface totale arrive d eau de la ct */
  const cs_real_t   cpa,      /* Capacite calorifique de l air */
  const cs_real_t   cpv,      /* Capacite calorifique de la vapeur */
  const cs_real_t   cpe,      /* Capacite calorifique de l eau */
  const cs_real_t   hv0,      /* Chaleur latente */
  const cs_real_t   rhoe,     /* Masse volumique de l eau*/
  const cs_real_t   dgout,    /* Diametre de goutte pour les zones de pluie */
  const cs_real_t   visc,     /* Viscosite Dynamique */
  const cs_real_t   conduc    /* Conductivite */
)
{
  cs_ctwr_zone_t  *ct;
  cs_int_t length;
  bft_file_t *f;
  char  *file_name = NULL;
  bft_file_type_t file_type;

  file_type = BFT_FILE_TYPE_TEXT;
  /* Définition d'une nouvelle zone d'echange */

  BFT_MALLOC(ct, 1, cs_ctwr_zone_t);

  ct->num = cs_glob_ct_nbr + 1;

  ct->idimct = idimct;
  ct->icoul = icoul;
  ct->imctch = imctch;
  ct->ntypct = ntypct;
  ct->nelect = nelect;

  ct->hmin  =  10000.;
  ct->hmax  = -10000.;

  ct->voiseau  = NULL;
  ct->pvoiseau = NULL;
  ct->voisair  = NULL;
  ct->pvoisair = NULL;

  ct->fac_sup_connect_idx = NULL;
  ct->fac_sup_connect_lst = NULL;

  ct->surf_fac_sup = NULL;
  ct->coefeau  = NULL;
  ct->coefair  = NULL;

  ct->deltat   = deltat;
  ct->cl_teau  = teau_cl;
  ct->cl_fem = fem_cl;
  ct->xap = xap;
  ct->xnp = xnp;

  ct->surface_in  = 0.;
  ct->surface_out = 0.;
  ct->surface     = surface;
  ct->nnpsct = 0;

  ct->nbfac_sct = 0;
  ct->nbfac_ict = 0;
  ct->nbfac_lct = 0;
  ct->nbfac_ct  = 0;
  ct->nbfbr_sct = 0;
  ct->nbfbr_ict = 0;
  ct->nbfbr_lct = 0;

  ct->nbevct = 0;

  ct->id_amont = 999;

  ct->face_sup_mesh = NULL;
  ct->face_inf_mesh = NULL;
  ct->face_lat_mesh = NULL;

  ct->cell_mesh    = NULL;
  ct->water_mesh   = NULL;

  ct->teau    = NULL;
  ct->fem     = NULL;
  ct->vgoutte = NULL;

  ct->fem_e   = 0.0 ;
  ct->fem_s   = 0.0 ;
  ct->teau_e  = 0.0 ;
  ct->teau_s  = 0.0 ;
  ct->heau_e  = 0.0 ;
  ct->heau_s  = 0.0 ;
  ct->tair_e  = 0.0 ;
  ct->tair_s  = 0.0 ;
  ct->xair_e  = 0.0 ;
  ct->xair_s  = 0.0 ;
  ct->hair_e  = 0.0 ;
  ct->hair_s  = 0.0 ;
  ct->debit_e = 0.0 ;
  ct->debit_s = 0.0 ;

  ct->cpa   = cpa ;
  ct->cpv   = cpv ;
  ct->cpe   = cpe ;
  ct->hv0   = hv0 ;
  ct->rhoe  = rhoe ;
  ct->dgout = dgout ;
  ct->visc  = visc ;
  ct->conduc  = conduc ;
  /* Redimensionnement du tableau des zones d'echange si necessaire */

  if (cs_glob_ct_nbr == cs_glob_ct_nbr_max) {
    cs_glob_ct_nbr_max = (cs_glob_ct_nbr_max + 1) ;
    BFT_REALLOC(cs_glob_ct_tab, cs_glob_ct_nbr_max, cs_ctwr_zone_t *);
  }

  /* Ajout dans le tableau des zones d'echanges */

  cs_glob_ct_tab[cs_glob_ct_nbr] = ct;
  cs_glob_ct_nbr += 1;

#if defined(HAVE_MPI)
  ct->cs_array_rank = NULL;
#endif
  ct->locat_air_water     = NULL;
  ct->locat_water_air     = NULL;
  ct->locat_cell_ct_upwind= NULL;

  ct->post_mesh_id = 0;

  ct->nnpsct_with_ghosts = 0;
  ct->n_ghost_npsct = 0;
  ct->water_halo    = NULL;


  if (cs_glob_rank_id <= 0) {
    length = strlen("bltctc.") + 3 ;
    BFT_MALLOC(file_name, length, char);
    sprintf(file_name, "bltctc.%02d", ct->num);

    f = bft_file_open(file_name,
                      BFT_FILE_MODE_APPEND,
                      file_type);

    bft_file_printf(f, "# BILANS POUR LA ZONE D\'ECHANGES \n");
    bft_file_printf(f, "# =================================\n");
    bft_file_printf(f,"\tTEMP\tFLUX A/E\tTA MOY SOR\t TE MOY SOR");
    bft_file_printf(f,"\tXA MOY SOR\tDEBI A ENT\tDEBI A SOR \n");
    f = bft_file_free(f);
    BFT_FREE(file_name);

  }



}

/*---------------------------------------------------------------------------
 * Solve the equation "matrix.x = b" with Cramer's rule.
 *
 * parameters:
 *   m[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   1 if matrix is singular, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_inverse_3x3(double  m[3][3],
             double  b[3],
             double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det =   m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - m[1][0]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + m[2][0]*(m[0][1]*m[1][2] - m[1][1]*m[0][2]);

  if (FVM_ABS(det) < _epsilon_denom)
    return 1;
  else
    det_inv = 1./det;

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - b[1]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + b[2]*(m[0][1]*m[1][2] - m[1][1]*m[0][2])) * det_inv;

  x1 = (  m[0][0]*(b[1]*m[2][2] - b[2]*m[1][2])
        - m[1][0]*(b[0]*m[2][2] - b[2]*m[0][2])
        + m[2][0]*(b[0]*m[1][2] - b[1]*m[0][2])) * det_inv;

  x2 = (  m[0][0]*(m[1][1]*b[2] - m[2][1]*b[1])
        - m[1][0]*(m[0][1]*b[2] - m[2][1]*b[0])
        + m[2][0]*(m[0][1]*b[1] - m[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */

  x[0] = x0; x[1] = x1; x[2] = x2;

  return 0;
}

/*----------------------------------------------------------------------------
 * Test of coplanarity
 *----------------------------------------------------------------------------*/

static int
_is_coplanar(const cs_real_t *coord,
             const cs_int_t   nvoi[cs_ctwr_nmaxvoi],
             const cs_int_t   nbvoi,
             cs_real_t        mat[4][4],
             cs_real_t        vectBase[3][3],
             const cs_int_t   decF,
             const cs_real_t  dh)
{
  cs_real_t det,norme1, norme2, min;
  cs_real_t tmpRes[3];
  cs_int_t  i,ii,iii,ind,i0,i1,i2;
  cs_int_t numi, numii, numiii;

  det = 0.0;
  i0 = 0;
  i1 = 0;
  i2 = 0;
  min = 1000.0;

  for(i = 0; i < 4 ; i++){
    ind = 0;
    for(ii = 0; ii < 4 ; ii++){
      if(ii != i){
        vectBase[0][ind] = mat[ii][1];
        vectBase[1][ind] = mat[ii][2];
        vectBase[2][ind] = mat[ii][3];
        ind++;
      }
    }
    CS_PRODUIT_VECTORIEL(tmpRes, vectBase[0],vectBase[1]);
    det += pow(-1,i) * mat[i][0] * CS_LOC_PRODUIT_SCALAIRE(vectBase[2], tmpRes);
  }

  if (CS_ABS(det) <= (0.000001*dh) ) {
    /*2D*/
    for(i=0; i< nbvoi ; i++){
      for(ii=i+1; ii< nbvoi ; ii++){
        for(iii=ii+1; iii< nbvoi ; iii++){
          numi   = nvoi[i]   + decF;
          numii  = nvoi[ii]  + decF;
          numiii = nvoi[iii] + decF;
          vectBase[0][0] = coord[3*numii   ] - coord[3*numi  ];
          vectBase[0][1] = coord[3*numii +1] - coord[3*numi+1];
          vectBase[0][2] = coord[3*numii +2] - coord[3*numi+2];
          vectBase[1][0] = coord[3*numiii  ] - coord[3*numi  ];
          vectBase[1][1] = coord[3*numiii+1] - coord[3*numi+1];
          vectBase[1][2] = coord[3*numiii+2] - coord[3*numi+2];
          CS_PRODUIT_VECTORIEL(tmpRes, vectBase[0],vectBase[1]);
          if(   (CS_LOC_MODULE(tmpRes) > (0.000001*dh))
           && (CS_ABS(CS_LOC_PRODUIT_SCALAIRE(vectBase[0],vectBase[1]))< min)
          ){
            i0 = i;
            i1 = ii;
            i2 = iii;
          }
        }
      }
    }
  }
  else{
    /*3D*/
    vectBase[0][0] = 1.0; vectBase[1][0] = 0.0; vectBase[2][0] = 0.0;
    vectBase[0][1] = 0.0; vectBase[1][1] = 1.0; vectBase[2][1] = 0.0;
    vectBase[0][2] = 0.0; vectBase[1][2] = 0.0; vectBase[2][2] = 1.0;
    return 3;
  }

  if (i0 == 0 && i1 == 0 && i2 == 0){
    /*1D*/
    vectBase[0][0] = coord[3*(nvoi[1]+ decF)  ] - coord[3*(nvoi[0]+ decF)  ];
    vectBase[0][1] = coord[3*(nvoi[1]+ decF)+1] - coord[3*(nvoi[0]+ decF)+1];
    vectBase[0][2] = coord[3*(nvoi[1]+ decF)+2] - coord[3*(nvoi[0]+ decF)+2];
    return 1;
  }
  else{
    vectBase[0][0] = coord[3*(nvoi[i1]+ decF)  ] - coord[3*(nvoi[i0]+ decF)  ];
    vectBase[0][1] = coord[3*(nvoi[i1]+ decF)+1] - coord[3*(nvoi[i0]+ decF)+1];
    vectBase[0][2] = coord[3*(nvoi[i1]+ decF)+2] - coord[3*(nvoi[i0]+ decF)+2];
    vectBase[1][0] = coord[3*(nvoi[i2]+ decF)  ] - coord[3*(nvoi[i0]+ decF)  ];
    vectBase[1][1] = coord[3*(nvoi[i2]+ decF)+1] - coord[3*(nvoi[i0]+ decF)+1];
    vectBase[1][2] = coord[3*(nvoi[i2]+ decF)+2] - coord[3*(nvoi[i0]+ decF)+2];
    CS_PRODUIT_VECTORIEL(tmpRes, vectBase[0],vectBase[1]);
    CS_PRODUIT_VECTORIEL(vectBase[1] ,tmpRes,vectBase[0]);
    norme1= CS_LOC_MODULE(vectBase[0]);
    norme2= CS_LOC_MODULE(vectBase[1]);
    vectBase[0][0] /= norme1 ; vectBase[1][0] /= norme2;
    vectBase[0][1] /= norme1 ; vectBase[1][1] /= norme2;
    vectBase[0][2] /= norme1 ; vectBase[1][2] /= norme2;
    vectBase[2][0] = 0.0; vectBase[2][1] = 0.0; vectBase[2][2] = 0.0;
    return 2;
  }

}

/*----------------------------------------------------------------------------
 * Inversion matrice (methode de Jordan)
 *----------------------------------------------------------------------------*/

static int
_invmat(cs_real_t mat[4][4],
        cs_real_t matInv[4][4],
        cs_int_t  idim)
{
  cs_real_t aux;
  cs_int_t i,j,k,err;

  err=1;

  for(i = 0; i < idim + 1; i++)
    for(j = 0; j < idim + 1; j++)
      matInv[i][j] = mat[i][j];


  i = 0;
  while (err == 1 && (i < (idim+1))){

    if (CS_ABS(matInv[i][i]) > 1.e-15){

      aux = 1.0/matInv[i][i];

      for(j = 0; j < idim + 1; j++)
         matInv[i][j] *=  aux;

      matInv[i][i]= aux;

      for(k = 0; k < i; k++){
        aux = matInv[k][i];
        for(j = 0; j < idim + 1; j++){
           matInv[k][j] -= aux * matInv[i][j];
        }
        matInv[k][i] = -aux *matInv[i][i];
      }

      for(k = i + 1; k < idim + 1; k++){
        aux = matInv[k][i];
        for(j = 0; j < idim + 1; j++){
           matInv[k][j] -= aux * matInv[i][j];
        }
        matInv[k][i] = -aux *matInv[i][i];
      }
      i++;
    }
    else{
      err =0;
    }
  }

  return err;
}

/*----------------------------------------------------------------------------
 * Calculation of the ponderation function
 *----------------------------------------------------------------------------*/

static cs_real_t
_weighting(const cs_real_t  dx,
           const cs_real_t  dy,
           const cs_real_t  dz,
           const cs_real_t  ouv,
           const cs_int_t   lf,
           const cs_real_t  epgauss,
           const cs_real_t  cx,
           const cs_real_t  cy,
           const cs_real_t  cz)
{
  cs_real_t pi, lambda;
  cs_real_t poids = 0.0;
  cs_real_t xy = sqrt(pow(dx/cx,2.) + pow(dy/cy,2.) + pow(dz/cz,2.));

  if (xy < ouv)
  {
    switch(lf) {
    case 1:
      poids=1.-xy/ouv;
      break;
    case 2:
      pi = acos(-1.);
      poids = 0.5*(1.+cos(pi*xy/ouv));
      break;
    case 3:
      lambda = ouv/sqrt(epgauss*log(10.));
      poids = exp(-pow((xy/lambda),2.));
      break;
    default:
      assert(lf == 1 || lf == 2 || lf == 3);
    }
  }

  return poids ;
}

/*----------------------------------------------------------------------------
 * Scalar product between the face normal vector and the gravity vector
 *----------------------------------------------------------------------------*/

static cs_real_t
_dot_product_ng(const cs_int_t   ifac,
                const cs_int_t   dim,
                const cs_real_t *surf_f,
                const cs_real_t  gravite[3],
                const cs_real_t  direction)
{
  cs_real_t n_sortant[3], g_uni[3];
  cs_int_t idim;
  cs_real_t aux1,aux2;

  n_sortant[2] = 0.0 ;
  g_uni[2] = 0.0;

  aux1 = CS_LOC_MODULE(gravite);

  for (idim = 0 ; idim < 3 ; idim++) {
    n_sortant[idim] = direction * surf_f[ifac*3+idim];
    g_uni[idim]= gravite[idim];
  }

  aux2 = CS_LOC_MODULE(n_sortant);
  for (idim = 0 ; idim < dim ; idim++) {
    n_sortant[idim] /= aux2;
    g_uni[idim] /= aux1 ;

  }

  return CS_LOC_PRODUIT_SCALAIRE(n_sortant , g_uni );

}


/*---------------------------------------------------------------------------
 *
 *---------------------------------------------------------------------------*/

static void
_search_height(cs_ctwr_zone_t   *ct,
               const cs_real_t  *gravite,
               cs_real_t        *hmin,
               cs_real_t        *hmax)
{
  cs_int_t    i, ifac, nb,nb_dist, axe, idx, ii, jj;
  cs_real_t   *lst_xyz_sup ;
  cs_real_t   *lst_xyz_inf ;
  cs_real_t   *lst_xyz_fi ;
  cs_real_t   *lst_xyz_fs ;
  cs_real_t   *hmin_dist;
  const fvm_coord_t *lst_xyz_dist = NULL;

  cs_real_t   aux;
  const fvm_lnum_t  *location_fac = NULL;

  cs_int_t  *faces_vtx_idx   = NULL;
  cs_int_t  *faces_vtx_lst   = NULL;

  const double tolerance = 0.1  ;
  fvm_nodal_t *fs_tmp_mesh = NULL;
  fvm_nodal_t *fi_tmp_mesh = NULL;



  double coeff[3], v_aux[3], vertex_coords[2];
  double  v_x, v_y ;
  double  v_f_x = 0., v_f_y = 0., v_f_z = 0.;
  double a[3][3] = {{0., 0., 0.} ,
                    {0., 0., 0.} ,
                    {0., 0., 0.} };

  double b_x[3] = {0., 0., 0. };
  double b_y[3] = {0., 0., 0. };
  double b_z[3] = {0., 0., 0. };

  double matrice[6] = {0., 0., 0. ,0., 0., 0. };



  nb = (cs_int_t) fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);
  BFT_MALLOC(lst_xyz_sup, nb*3, fvm_coord_t);
  fvm_nodal_get_vertex_coords(ct->face_sup_mesh,
                              FVM_INTERLACE,
                              lst_xyz_sup);

  nb = (cs_int_t) fvm_nodal_get_n_entities(ct->face_inf_mesh, 0);
  BFT_MALLOC(lst_xyz_inf, nb*3, fvm_coord_t);
  fvm_nodal_get_vertex_coords(ct->face_inf_mesh,
                              FVM_INTERLACE,
                              lst_xyz_inf);

  fs_tmp_mesh = fvm_nodal_copy(ct->face_sup_mesh);
  fi_tmp_mesh = fvm_nodal_copy(ct->face_inf_mesh);

  aux = 0.;

  for (i = 0 ; i < 3 ; i++)
    if(CS_ABS(gravite[i]) > aux) {
      axe = i;
      aux = CS_ABS (gravite [ i ]);
    }

  if (axe == 0) {
    matrice[1] = 1;
    matrice[5] = 1;
  }
  else {
    matrice[0] = 1;
    if (axe == 1)
      matrice[5] = 1;
    else
      matrice[4] = 1;
  }


  fvm_nodal_project_coords(fs_tmp_mesh, matrice);
  fvm_nodal_project_coords(fi_tmp_mesh, matrice);

  nb = (cs_int_t) fvm_nodal_get_n_entities(fs_tmp_mesh, 0);

  BFT_MALLOC(lst_xyz_fs, nb*2, fvm_coord_t);

  fvm_nodal_get_vertex_coords(fs_tmp_mesh, FVM_INTERLACE, lst_xyz_fs);


  nb = (cs_int_t) fvm_nodal_get_n_entities(fi_tmp_mesh, 0);

  BFT_MALLOC(lst_xyz_fi , nb*2 , fvm_coord_t );

  fvm_nodal_get_vertex_coords(fi_tmp_mesh, FVM_INTERLACE, lst_xyz_fi);

  /* Create locator for interpolate */

  fvm_locator_t   *locator = NULL;

#if defined(FVM_HAVE_MPI)
  locator = fvm_locator_create(tolerance,
                               cs_glob_mpi_comm,
                               cs_glob_n_ranks,
                               0);
#else
  locator = fvm_locator_create(tolerance);
#endif

  nb = (cs_int_t) fvm_nodal_get_n_entities(fs_tmp_mesh, 0);


  fvm_locator_set_nodal(locator,
                        fi_tmp_mesh,
                        0,
                        2,
                        nb,
                        NULL,
                        lst_xyz_fs);

  nb_dist = fvm_locator_get_n_dist_points(locator );

  BFT_MALLOC(hmin_dist ,  nb_dist, fvm_coord_t );

  cs_reverse_vtx_faces_connect(fi_tmp_mesh        ,
                                &(faces_vtx_idx ) ,
                                &(faces_vtx_lst ) );

  location_fac = fvm_locator_get_dist_locations(locator );
  lst_xyz_dist = fvm_locator_get_dist_coords(   locator );

  for (i = 0 ; i < nb_dist ; i++) {

    ifac = location_fac[ i ] - 1;

    vertex_coords [0] = lst_xyz_dist [ i*2     ];
    vertex_coords [1] = lst_xyz_dist [ i*2 + 1 ];

    for(ii = 0; ii < 3 ; ii++){
      b_x[ ii ] = 0. ;
      b_y[ ii ] = 0. ;
      b_z[ ii ] = 0. ;
      for(jj = 0; jj < 3 ; jj++)
        a[ ii ][ jj ] = 0.;

    }


    for (idx = faces_vtx_idx[ ifac    ];
         idx < faces_vtx_idx[ ifac +1 ]; idx++) {

      v_x = lst_xyz_fi[faces_vtx_lst[ idx ]* 2     ];
      v_y = lst_xyz_fi[faces_vtx_lst[ idx ]* 2 + 1 ];

      v_f_x = lst_xyz_inf[faces_vtx_lst[ idx ]* 3     ];
      v_f_y = lst_xyz_inf[faces_vtx_lst[ idx ]* 3 + 1 ];
      v_f_z = lst_xyz_inf[faces_vtx_lst[ idx ]* 3 + 2 ];

      a[0][0] += v_x * v_x;
      a[0][1] += v_x * v_y;
      a[0][2] += v_x;

      a[1][1] += v_y * v_y;
      a[1][2] += v_y;


      a[2][2] += 1.;

      b_x[0] += v_x * v_f_x;
      b_x[1] += v_y * v_f_x;
      b_x[2] += v_f_x;

      b_y[0] += v_x * v_f_y;
      b_y[1] += v_y * v_f_y;
      b_y[2] += v_f_y;

      b_z[0] += v_x * v_f_z;
      b_z[1] += v_y * v_f_z;
      b_z[2] += v_f_z;


    }

    /* Matrix is symmetric */

    a[1][0] = a[0][1];
    a[2][0] = a[0][2];
    a[2][1] = a[1][2];

    if (_inverse_3x3(a, b_x, coeff) == 0 ) {

      v_aux[0] = -(  coeff[0]*vertex_coords[0]
                   + coeff[1]*vertex_coords[1]
                   + coeff[2]);
    }
    else
       v_aux[0] = -v_f_x;

    if (_inverse_3x3(a, b_y, coeff) == 0) {

      v_aux[1] = -(  coeff[0]*vertex_coords[0]
                   + coeff[1]*vertex_coords[1]
                   + coeff[2]);
    }
    else
      v_aux[1] = -v_f_y;

    if (_inverse_3x3(a, b_z, coeff) == 0) {

      v_aux[2] = -(  coeff[0]*vertex_coords[0]
                   + coeff[1]*vertex_coords[1]
                   + coeff[2]);
    }
    else
      v_aux[2] = -v_f_z;

    hmin_dist[i] =  CS_LOC_PRODUIT_SCALAIRE(v_aux , gravite)
                  / CS_LOC_MODULE(gravite);
  }


  fvm_locator_exchange_point_var(locator,
                                 hmin_dist, hmin, NULL, sizeof(cs_real_t),1,0);

  for (i = 0 ; i < nb ; i++) {

      v_aux[0] = -lst_xyz_sup[i*3    ];/* Opposite Vector to g */
      v_aux[1] = -lst_xyz_sup[i*3 + 1];
      v_aux[2] = -lst_xyz_sup[i*3 + 2];

      aux = CS_LOC_PRODUIT_SCALAIRE(v_aux , gravite);
      hmax[i] = aux / CS_LOC_MODULE(gravite); /* project on "g" axis */

  }

   BFT_FREE(lst_xyz_inf);
   BFT_FREE(lst_xyz_sup);
   BFT_FREE(lst_xyz_fi);
   BFT_FREE(lst_xyz_fs);
   BFT_FREE(hmin_dist);

   locator = fvm_locator_destroy(locator);
   fs_tmp_mesh = fvm_nodal_destroy(fs_tmp_mesh );
   fi_tmp_mesh = fvm_nodal_destroy(fi_tmp_mesh );
}


/*----------------------------------------------------------------------------
 * Function cs_ctwr_maille
 * Construction du maillage eau
 *----------------------------------------------------------------------------*/

void cs_ctwr_maille
(
  const cs_real_t          gx,             /* composante x de la gravite      */
  const cs_real_t          gy,             /* composante y de la gravite      */
  const cs_real_t          gz,             /* composante z de la gravite      */
  const cs_mesh_t       *mesh,             /* <-- structure maillage associée */
  const cs_mesh_quantities_t *mesh_quantities   /* <-- grandeurs du maillage  */
)
{


  cs_int_t   icel_1, icel_2, ii, icel, length, nb, rank,
             dist_rank, res_loc, res_dist;
  cs_int_t   ifac, ict, icpt, icpti, icptl, icptla, icptfac,
             iaux, i, dim, j;
  cs_real_t  aux       , gravite[3], v_aux[3] , alpha ;
  cs_int_t    *lst_aux;
  fvm_coord_t *extrusion_vectors, *lst_xyz_cel, *lst_xyz;
  fvm_lnum_t  *lst_par_fac_sup;
  fvm_gnum_t  *fsup_gb_vt_num = NULL;
  cs_real_t   *hmin_vect ;
  cs_real_t   *hmax_vect ;

  char  *mesh_name       = NULL ;
  char  *export_name     = NULL ;
  const double tolerance = 0.1  ;


  fvm_gnum_t   n_vertices;

  cs_int_t   *face_sup;      /* liste des faces internes superieures de la ct */
                             /* de taille  (nbfac_sct )                       */
  cs_int_t   *fbr_sup;       /* liste des faces de bord superieures de la ct  */
                             /* de taille  (nbfbr_sct)                       */
  cs_int_t   *face_inf;      /* liste des faces internes inferieures de la ct */
                             /* de taille  (nbfac_ict)                        */
  cs_int_t   *fbr_inf;       /* liste des faces de bord  inferieures de la ct */
                             /* de taille  (nbfac_ict)                        */
  cs_int_t   *face_lat;      /* liste des faces internes laterales de la ct   */
                             /* de taille  (nbfac_lct )                       */
  cs_int_t   *fbr_lat;       /* liste des faces de bord laterales de la ct    */
                             /* de taille  (nbfbr_lct)                       */
  cs_int_t   *face_ct;       /* liste des faces interne de la ct              */
                             /* de taille  (nbfac_ct)                        */


  const cs_int_t  *i_face_cells  = mesh->i_face_cells;
  const cs_int_t  *b_face_cells  = mesh->b_face_cells;
  const cs_real_t *i_face_normal = mesh_quantities->i_face_normal;
  const cs_real_t *b_face_normal = mesh_quantities->b_face_normal;
  const cs_int_t  *family_item   = mesh->family_item ;
  const cs_int_t  *cell_family   = mesh->cell_family;
  const cs_int_t  nbr_cel        = mesh->n_cells;

  fvm_interface_set_t  *interface_set = NULL;
  cs_ctwr_zone_t  *ct;


  iaux = 0;
  alpha = 0.875;

  /* Vecteur gravite */
  gravite[0] = gx;
  gravite[1] = gy;
  gravite[2] = gz;

  /*--------------------------------------------*/
  /* Numbers of air noode in the EA: nbevct     */
  /*--------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {
   ct = cs_glob_ct_tab[ict];
   dim = ct->idimct;
    for (icel = 0 ; icel <  nbr_cel ; icel++) {
      if (ct->icoul == family_item[cell_family[icel]-1]){
        ct->nbevct++;
      }
    }
  }
  /*--------------------------------------------*/
  /* End Numbers of air noode in the EA: nbevct */
  /*--------------------------------------------*/

  /*--------------------------------------------*/
  /* List of air nodes for each Exchange Area   */
  /*--------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {

    icpt = 0;
    ct = cs_glob_ct_tab[ict];
    BFT_MALLOC(lst_aux,ct->nbevct,cs_int_t);

    for (icel = 0 ; icel <  nbr_cel; icel++) {
      if (ct->icoul == family_item[cell_family[icel]-1]) {
        assert(icpt < ct->nbevct);
        lst_aux[++icpt-1] = icel + 1;
      }
    }

    length = strlen("cell_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "cell_mesh_ct_%d", ict);


    ct->cell_mesh = cs_mesh_connect_cells_to_nodal(mesh,
                                                   mesh_name,
                                                   ct->nbevct,
                                                   lst_aux);
  }
  /*---------------------------------------------*
   * End list of air nodes for each Exchange Area*
   *---------------------------------------------*/

  /*--------------------------------------------------------*
   * Calcul du nombre de noeuds eau des faces superieures   *
   * des zones d'echanges et du nombre de faces superieures *
   * et inferieures                                         *
   *--------------------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {
    ct = cs_glob_ct_tab[ict];
    /* Contribution faces internes */
    for (ifac = 0 ; ifac < mesh->n_i_faces ; ifac++) {
      assert((ifac * 2 + 1) < (2*mesh->n_i_faces));
      icel_1 = i_face_cells[ifac * 2]     - 1;/* indice de la cellule 1 */
      icel_2 = i_face_cells[ifac * 2 + 1] - 1;/* indice  de la cellule 2 */
      /* Comparaison  des couleurs des cellules 1 et 2 */
      if((family_item[cell_family[icel_1]-1] == ct->icoul ) ||
         (family_item[cell_family[icel_2]-1] == ct->icoul)) {
        if  (family_item[cell_family[icel_1]-1] != family_item[cell_family[icel_2]-1]) {
          if (family_item[cell_family[icel_1]-1] == ct->icoul ) {
            aux = _dot_product_ng(ifac ,ct->idimct, i_face_normal, gravite, 1);
          }
          if (family_item[cell_family[icel_2]-1] == ct->icoul ) {
            aux = _dot_product_ng(ifac ,ct->idimct, i_face_normal, gravite, -1);
          }

          if (aux < (-alpha) ){
            ct->nnpsct++;
            ct->nbfac_sct++;
          }else{
            if (aux > alpha ){
              ct->nbfac_ict++;
            }else{
              ct->nbfac_lct++;
            }
          }
        }else{
          ct->nbfac_ct++;
        }
      }

    }  /* fin contribution faces internes */

    /* Contribution faces externes */
    for (ifac = 0 ; ifac < mesh->n_b_faces ; ifac++) {
      icel_1 = b_face_cells[ifac] - 1; /* indice de la cellule  */
      if (family_item[cell_family[icel_1]-1] == ct->icoul ) {

        aux = _dot_product_ng(ifac,ct->idimct, b_face_normal, gravite, 1);

        if (aux < (-alpha) ){
          ct->nnpsct++;
          ct->nbfbr_sct++;
        }else{
          if (aux > alpha  ){
            ct->nbfbr_ict++;
          }else{
            ct->nbfbr_lct++;
          }
        }
      }
    }/* fin contribution faces externes */


    /* allocation memoire pour la liste des faces superieures et inferieures
    * des ct */
    BFT_MALLOC(face_sup,ct->nbfac_sct ,cs_int_t );
    BFT_MALLOC(face_inf,ct->nbfac_ict ,cs_int_t );
    BFT_MALLOC(face_lat,ct->nbfac_lct ,cs_int_t );
    BFT_MALLOC(fbr_sup ,ct->nbfbr_sct ,cs_int_t );
    BFT_MALLOC(fbr_inf ,ct->nbfbr_ict ,cs_int_t );
    BFT_MALLOC(fbr_lat ,ct->nbfbr_lct ,cs_int_t );
    BFT_MALLOC(face_ct ,ct->nbfac_ct  ,cs_int_t );


  /* --------------------------------------------------------*
   * Fin Calcul du nombre de noeuds eau des faces superieures*
   * des zones d'echanges et du nombre de faces superieures  *
   * et inferieures                                          *
   *---------------------------------------------------------*/


  /*-----------------------------------------------------------------*
   * Liste des faces superieures et inferieures des zones d echanges *
   * et liste des noeuds eau des faces sup de la ct sans ct amont    *
   *-----------------------------------------------------------------*/

    /* Contribution faces internes */
    icpt   = 0 ; /*indice tableau des faces  sup */
    icpti  = 0 ; /*indice tableau des faces  inf */
    icptla = 0 ; /*indice tableau des faces  laterales */
    icptl  = 0 ; /*indice tableau des noeuds sup ct */
    icptfac  = 0 ; /*indice tableau des noeuds sup ct */
    /* Boucle sur les faces internes du domaine */
    for (ifac = 0 ; ifac < mesh->n_i_faces ; ifac++) {
      icel_1 = i_face_cells[ifac * 2]     - 1; /* indice de la cellule 1 */
      icel_2 = i_face_cells[ifac * 2 + 1] - 1; /* indice  de la cellule 2 */
      /* Comparaison  couleur de la ct et couleur des cellules 1 et 2 */
      if((family_item[cell_family[icel_1]-1] == ct->icoul )||
         (family_item[cell_family[icel_2]-1] == ct->icoul)){

        if (family_item[cell_family[icel_1]-1] != family_item[cell_family[icel_2]-1] ){

          if (family_item[cell_family[icel_1]-1] == ct->icoul ) {
            aux = _dot_product_ng(ifac ,ct->idimct, i_face_normal, gravite, 1);
          }
          if (family_item[cell_family[icel_2]-1] == ct->icoul ) {
            aux = _dot_product_ng(ifac ,ct->idimct, i_face_normal, gravite, -1);
          }

          if (aux < (-alpha) ){
            /*ajout d'une face sup de la ct*/
            face_sup[icpt] = ifac + 1;
            icpt ++;
          }else{
            if (aux > alpha ){
            /*ajout d'un face inf de la ct*/
            face_inf[icpti]  = ifac + 1;
            icpti ++;

            }else{
              assert(icptla < ct->nbfac_lct+ct->nbfbr_lct);
              face_lat[icptla] = ifac  + 1;
              icptla ++;
            }
          }
        }else{
          face_ct[icptfac] = ifac  + 1;
          icptfac ++;
        }
      }
    }/* fin contribution faces internes */

    /* Contribution faces de bords */
    /* initialisation des indices */
    icpt   = 0 ; /*indice tableau des faces  sup */
    icpti  = 0 ; /*indice tableau des faces  inf */
    icptla = 0 ; /*indice tableau des faces  laterales */



    for (ifac = 0 ; ifac < mesh->n_b_faces ; ifac++) {

      icel_1 = b_face_cells[ifac] - 1;/* indice de la cellule  */
      if (family_item[cell_family[icel_1]-1] == ct->icoul ) {

        aux = _dot_product_ng(ifac, ct->idimct, b_face_normal, gravite, 1);

        if (aux < (-alpha) ){
          /* ajout d'une face sup de la ct */
          fbr_sup[icpt]= ifac + 1;
          icpt ++;
        }else{
          if (aux > alpha ){
            /*ajout d'un face inf de la ct*/
            fbr_inf[icpti]= ifac + 1;
            icpti ++;
          }else{
            fbr_lat[icptla]= ifac + 1;
            icptla ++;
          }
        }
      }
    } /* fin contribution faces externes */


    /*---------------------------------------------------------*
    * Creation des maillages surfacique en connectivité nodale*
    *---------------------------------------------------------*/

    /* mesh for superiors faces */

    BFT_FREE(mesh_name);
    length = strlen("face_sup_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "face_sup_mesh_ct_%d", ict);


    ct->face_sup_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                       mesh_name,
                                                       ct->nbfac_sct,
                                                       ct->nbfbr_sct,
                                                       face_sup,
                                                       fbr_sup);

    /* mesh for inferiors faces*/
    BFT_FREE(mesh_name);
    length = strlen("face_inf_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "face_inf_mesh_ct_%d", ict);


    ct->face_inf_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                       mesh_name,
                                                       ct->nbfac_ict,
                                                       ct->nbfbr_ict,
                                                       face_inf,
                                                       fbr_inf);
    /* mesh for laterals faces*/
    BFT_FREE(mesh_name);
    length = strlen("face_lat_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "face_lat_mesh_ct_%d", ict);


    ct->face_lat_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                       mesh_name,
                                                       ct->nbfac_lct,
                                                       ct->nbfbr_lct,
                                                       face_lat,
                                                       fbr_lat);
    /* mesh for laterals faces*/
    BFT_FREE(mesh_name);
    length = strlen("face_mesh_ct_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "face_mesh_ct_%d", ict);


    ct->fac_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                  mesh_name,
                                                  ct->nbfac_ct,
                                                  0,
                                                  face_ct,
                                                  NULL);
    /* water mesh*/
    BFT_FREE(mesh_name);
    length = strlen("water_mesh_") + 1 + 1;
    BFT_MALLOC(mesh_name, length, char);
    sprintf(mesh_name, "water_mesh_%d", ict);


    ct->water_mesh = cs_mesh_connect_faces_to_nodal(mesh,
                                                    mesh_name,
                                                    ct->nbfac_sct,
                                                    ct->nbfbr_sct,
                                                    face_sup,
                                                    fbr_sup);


    /*--------------------------------------------------------------*
    *  Fin creation des maillages surfacique en connectivité nodale*
    *--------------------------------------------------------------*/

    /*--------------------------------------------------------------*
    * Construct cs_array_rank                                  *
    *--------------------------------------------------------------*/

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {

      nb   = cs_glob_n_ranks;
      rank = cs_glob_rank_id;
      BFT_MALLOC(ct->cs_array_rank, nb, cs_int_t );


      ct->cs_array_rank[ rank ] = res_loc = ct->nbevct;


      for(dist_rank = 0; dist_rank <  nb; dist_rank++ )
        if(dist_rank != rank ){
          MPI_Sendrecv(&res_loc,  1, CS_MPI_INT, dist_rank, CS_CT_MPI_TAG,
                        &res_dist, 1, CS_MPI_INT, dist_rank, CS_CT_MPI_TAG,
                        cs_glob_mpi_comm, &status);

          ct->cs_array_rank[ dist_rank ] = res_dist;

        }
    }
#endif

    /*--------------------------------------------------------------*
    * End of Construction cs_array_rank                         *
    *--------------------------------------------------------------*/

    /*--------------------------------------------------------------*
    * Reseach of  hmin and hmax                                    *
    *--------------------------------------------------------------*/

    /* loop on the superior faces for hmax */
    nb = (cs_int_t) fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);


    BFT_MALLOC(hmax_vect , nb , fvm_coord_t );
    BFT_MALLOC(hmin_vect , nb , fvm_coord_t );

    _search_height(ct          ,
                     gravite    ,
                     hmin_vect  ,
                     hmax_vect  );


    for (i = 0 ; i < nb ; i++) {

        aux = hmax_vect[ i ];
        if (aux >= ct->hmax )
          ct->hmax = aux;

        aux = hmin_vect[ i ];;
        if (aux <= ct->hmin )
          ct->hmin = aux;
    }

    /* loop on the sup faces for surface_in and surface_out */
    BFT_MALLOC(lst_par_fac_sup , ct->nnpsct , fvm_lnum_t );

    fvm_nodal_get_parent_num(ct->face_sup_mesh, 2, lst_par_fac_sup);

    BFT_MALLOC(ct->surf_fac_sup , ct->nnpsct , cs_real_t );

    for (ifac = 0 ; ifac < ct->nnpsct ; ifac++) {
      if(ifac< ct->nbfbr_sct ){
        for (ii = 0 ; ii < 3 ; ii++){
          v_aux[ii] = b_face_normal[3 * (cs_int_t) (lst_par_fac_sup[ ifac ] -1)
                               + ii ];
        }
      }
      else{
        for (ii = 0 ; ii < 3 ; ii++){
          v_aux[ii] = i_face_normal[ 3 * (cs_int_t)
                                (lst_par_fac_sup[ifac] - mesh->n_b_faces - 1)
                                + ii ];
          }
      }
      aux = CS_LOC_MODULE(v_aux);
      ct->surface_in += aux;
      ct->surface_out += aux;
      ct->surf_fac_sup[ifac] = aux;

    }



#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {

      nb   = cs_glob_n_ranks;
      rank = cs_glob_rank_id;

      // TODO : changer ce bordel !!!!!!!!!!!!!!!!!
      // sans doute equivalent a MPI_Allreduce(ct-hmax, ..., MPI_MAX)

      if(ct->cs_array_rank[ rank ] != 0){
        for(dist_rank = 0; dist_rank < nb; dist_rank++){
          if(dist_rank != rank ){
            if(ct->cs_array_rank [ dist_rank ] != 0 ){

              MPI_Sendrecv(&ct->hmax, 1, CS_MPI_REAL, dist_rank, CS_CT_MPI_TAG,
                           &aux, 1, CS_MPI_REAL, dist_rank, CS_CT_MPI_TAG,
                           cs_glob_mpi_comm, &status);

              if ( aux > ct->hmax ) ct->hmax = aux;

              MPI_Sendrecv(&ct->hmin, 1, CS_MPI_REAL, dist_rank, CS_CT_MPI_TAG,
                           &aux  , 1, CS_MPI_REAL, dist_rank, CS_CT_MPI_TAG,
                         cs_glob_mpi_comm, &status);

              if ( aux < ct->hmin ) ct->hmin = aux;

            }
          }
        }
      }

      MPI_Allreduce (&ct->surface_in, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->surface_in = aux;

      MPI_Allreduce (&ct->surface_out, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->surface_out = aux;
    }
#endif

   /* -------------------------------------------------------------*
    * End of Reseach  hmin et hmax                                  *
    *---------------------------------------------------------------*/


    nb = fvm_nodal_get_n_entities(ct->water_mesh, 0);



    BFT_MALLOC(extrusion_vectors, (nb*3), fvm_coord_t);


    for (i=0 ; i < nb ; i++){

      aux =CS_ABS(hmax_vect[ i ] -  hmin_vect[ i ])/CS_LOC_MODULE(gravite);

      extrusion_vectors[ i*3 ]    =  gravite[0] * aux;
      extrusion_vectors[ i*3 + 1] =  gravite[1] * aux;
      extrusion_vectors[ i*3 + 2] =  gravite[2] * aux;
    }



    fvm_nodal_extrude(ct->water_mesh,
                      ct->nelect,
                      extrusion_vectors,
                      NULL);

    BFT_FREE(extrusion_vectors);


    /* Set halo structure for the water mesh */

    n_vertices = fvm_nodal_get_n_entities(ct->face_sup_mesh, 0);

    BFT_MALLOC(fsup_gb_vt_num, n_vertices, fvm_gnum_t);

    fvm_nodal_get_global_vertex_num(ct->face_sup_mesh, fsup_gb_vt_num);

    interface_set = fvm_interface_set_create(n_vertices, NULL, fsup_gb_vt_num,
                                             NULL, 0, NULL, NULL, NULL);

    /*fvm_interface_set_dump(interface_set);*/

    /* Creation of the cs_mesh_halo_t structure. */

    bft_printf(_(" Création des halos\n"));
    bft_printf_flush();

    ct->water_halo = cs_halo_create(interface_set);

    bft_printf(_(" Définition des halos\n"));
    bft_printf_flush();

    cs_ctwr_halo_define(ct, interface_set);

    fvm_interface_set_destroy(interface_set);

    /* Create locator for interpolate */

#if defined(FVM_HAVE_MPI)
    ct->locat_water_air = fvm_locator_create(tolerance,
                                             cs_glob_mpi_comm,
                                             cs_glob_n_ranks,
                                             0);
#else
    ct->locat_water_air = fvm_locator_create(tolerance);
#endif

    BFT_MALLOC(lst_xyz_cel , ct->nbevct*3, fvm_coord_t);

    fvm_nodal_get_element_centers(ct->cell_mesh, FVM_INTERLACE, 3, lst_xyz_cel);

    fvm_locator_set_nodal(ct->locat_water_air,
                          ct->water_mesh,
                          0,
                          3,
                          ct->nbevct,
                          NULL,
                          lst_xyz_cel);


#if defined(FVM_HAVE_MPI)
    ct->locat_air_water = fvm_locator_create(tolerance,
                                             cs_glob_mpi_comm,
                                             cs_glob_n_ranks,
                                             0);
#else
    ct->locat_air_water = fvm_locator_create(tolerance);
#endif

    BFT_MALLOC(lst_xyz, ct->nnpsct*ct->nelect*3, fvm_coord_t );

    fvm_nodal_get_element_centers(ct->water_mesh, FVM_INTERLACE, 3, lst_xyz);


    fvm_locator_set_nodal(ct->locat_air_water,
                          ct->cell_mesh,
                          0,
                          3,
                          ct->nnpsct*ct->nelect,
                          NULL,
                          lst_xyz);


    BFT_FREE(mesh_name       );
    BFT_FREE(export_name     );
    BFT_FREE(face_sup        );
    BFT_FREE(face_inf        );
    BFT_FREE(face_lat        );
    BFT_FREE(fbr_sup         );
    BFT_FREE(lst_par_fac_sup );
    BFT_FREE(fbr_inf         );
    BFT_FREE(fbr_lat         );
    BFT_FREE(face_ct         );
    BFT_FREE(lst_aux         );
    BFT_FREE(lst_xyz         );
    BFT_FREE(lst_xyz_cel     );
    BFT_FREE(fsup_gb_vt_num  );

  }
  /*--------------------------------------------*
   * Fin Liste des faces superieures et inferieures des zones d echanges
   * et liste des noeuds eau des faces sup de la ct*
   *--------------------------------------------*/

  /*--------------------------------------------*
   * Initialization of the water variables      *
   *--------------------------------------------*/

  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {
    ct = cs_glob_ct_tab[ict];
    /* Te */
    BFT_MALLOC(ct->teau   ,(ct->nnpsct_with_ghosts*ct->nelect), cs_real_t);
    /* Fe */
    BFT_MALLOC(ct->fem    ,(ct->nnpsct_with_ghosts*ct->nelect), cs_real_t);
    /* vg */
    BFT_MALLOC(ct->vgoutte,(ct->nnpsct_with_ghosts*ct->nelect), cs_real_t);

    /* initialisation*/
    for (iaux = 0 ; iaux < (ct->nnpsct_with_ghosts*ct->nelect) ; iaux++) {
      /* temperature de l eau*/
      ct->teau[iaux]    = ct->cl_teau;
      /* debit massique par unite de surface */
      ct->fem[iaux]     = ct->cl_fem/ct->surface;
      /* vitesse des gouttes */
      ct->vgoutte[iaux] = 0.0;
    }
    /* Initialisation en tenant compte de l'écart de température imposé*/
    aux = ct->deltat / ( (cs_real_t) (ct->nelect - 1)  ) ;
    for (i = 0 ; i < ct->nnpsct_with_ghosts ; i++) {
      for (j = 1 ; j < ct->nelect ; j++) {
          ii = i*ct->nelect + j ;
          ct->teau[ ii ] =  ct->teau[ ii - 1 ] - aux;
        }
    }


  }/*  fin de la boucle sur les zones d'echanges */


  /*----------------------------------------------*
   * End of fnitialization of the water variables *
   *----------------------------------------------*/


  /*--------------------------------------------*
   * Initialisation des tableaux d interpolation*
   *--------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {

    ct = cs_glob_ct_tab[ ict ];
    /* Liste des voisins eau des cellules air */
    nb = (int) fvm_locator_get_n_dist_points(ct->locat_air_water );
    BFT_MALLOC(ct->voiseau,(nb * cs_ctwr_nmaxvoi ), cs_int_t );
    /* Coefficients d interpolation eau pour l air*/
    BFT_MALLOC(ct->coefeau , (nb * cs_ctwr_nmaxvoi ), cs_real_t );
    /* Positions dans la liste des voisins eau */
    BFT_MALLOC(ct->pvoiseau, (nb + 1) , cs_int_t );
    /* Liste des voisins air des noeuds eau */

    ct->pvoiseau[ 0 ] = 0 ;

    for (iaux = 0 ; iaux < (nb *cs_ctwr_nmaxvoi) ; iaux++){
      ct->voiseau[iaux] = -1 ;
    }

    nb = (int) fvm_locator_get_n_dist_points(ct->locat_water_air );

    BFT_MALLOC(ct->voisair  ,(nb * cs_ctwr_nmaxvoi ), cs_int_t  );
    /* Positions dans la liste voisins air */
    BFT_MALLOC(ct->pvoisair ,(    nb + 1   ), cs_int_t  );
    /* Coefficients d interpolation air pour l eau */
    BFT_MALLOC(ct->coefair  ,( nb * cs_ctwr_nmaxvoi ), cs_real_t );

    ct->pvoisair[0] = 0 ;

    for (iaux = 0 ; iaux < (nb* cs_ctwr_nmaxvoi) ; iaux++){
      ct->voisair[iaux] = -1 ;
    }
  }
  /*------------------------------------------------------*
   * Fin de l initialisation des tableaux d interpolation *
   *------------------------------------------------------*/


}

/*----------------------------------------------------------------------------
 * Function cs_ctwr_adeau
 * Interpolation AIR -> EAU
 *----------------------------------------------------------------------------*/

void cs_ctwr_adeau
(
  const cs_mesh_t             *mesh,
  const cs_mesh_quantities_t  *mesh_quantities
)
{
  /* Coordonnées des centres des cellules  */
  const cs_real_t *coo_cel        = mesh_quantities->cell_cen;
  const cs_int_t  *i_face_cells   = mesh->i_face_cells       ;
  const cs_int_t  *cell_cells_idx = mesh->cell_cells_idx     ;
  const cs_int_t  *cell_cells_lst = mesh->cell_cells_lst     ;
  const cs_int_t  *family_item    = mesh->family_item        ;
  const cs_int_t  *cell_family    = mesh->cell_family        ;

  cs_int_t   ict, iwat,nb_node_water, ii, jj, iair, nbvois,
             nbn, nvois[ cs_ctwr_nmaxvoi ], ifac, icel_1, icel_2,icel, lf, indice, dim;
  cs_real_t  dhi, dmin;
  cs_real_t  xwat, ywat, zwat,
             dx, dy, dz, dxx, dyy, dzz, coeff[ cs_ctwr_nmaxvoi ], ouv, aux;
  cs_real_t  vectBase[3][3];
  cs_real_t  cx, cy, cz, epgauss, w,
             pp[4][4],ppInv[4][4];

  const fvm_coord_t *lst_xyz_water = NULL;
  fvm_coord_t *lst_xyz_cel  ;
  fvm_lnum_t  *lst_par_fac  ;
  fvm_lnum_t  *lst_par_cel  ;

  const fvm_lnum_t  *location_cel  = NULL;

  /*--------------------------------------------*
   * parametres et initialisation               *
   *--------------------------------------------*/
  cs_ctwr_zone_t  *ct;

  /*--------------------------------------------*
   * fin parametres et initialisation           *
   *--------------------------------------------*/

  /* Make sure with have extended neighborhood */

  assert(cell_cells_idx != NULL);

  /*---------------------------------------------*
   * Construction des coefficient d'interpolation*
   * sur chaque zone d echange ict               *
   *---------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {

    ct = cs_glob_ct_tab[ict];

    nbn = 3 ;
    if (ct->idimct==3) {
      nbn = 4 ;
    }

    /* Calcul de dh */
    dhi = (ct->hmax-ct->hmin)/(ct->nelect-1);

    /* Copy element centers of the water mesh to an array.*/

    nb_node_water = (int) fvm_locator_get_n_dist_points(ct->locat_air_water);

    /* */
    location_cel = fvm_locator_get_dist_locations(ct->locat_air_water);

    /* */
    lst_xyz_water   = fvm_locator_get_dist_coords(ct->locat_air_water);

    BFT_MALLOC(lst_xyz_cel , ct->nbevct*3, fvm_coord_t );
    fvm_nodal_get_element_centers(ct->cell_mesh, FVM_INTERLACE , 3 , lst_xyz_cel);

    BFT_MALLOC(lst_par_cel , ct->nbevct, fvm_lnum_t );
    fvm_nodal_get_parent_num(ct->cell_mesh, 3, lst_par_cel);

    BFT_MALLOC(lst_par_fac , ct->nbfac_ct, fvm_lnum_t );
    fvm_nodal_get_parent_num(ct->fac_mesh, 2, lst_par_fac);


    /* boucle sur les noeuds eau */
    for (iwat = 0 ; iwat < nb_node_water ; iwat++) {

      /*--------------------------------------------*
       * Calcul des coord. du noeud a partir de     *
       *  celles du noeud de la face sup            *
       *   Noeud  = NoeudSup + iloc * dh *g / ||g|| *
       *--------------------------------------------*/
      xwat = (cs_real_t) lst_xyz_water[ iwat*3     ] ;
      ywat = (cs_real_t) lst_xyz_water[ iwat*3 + 1 ] ;
      zwat = (cs_real_t) lst_xyz_water[ iwat*3 + 2 ] ;

      /*--------------------------------------------*
       * boucle sur les cellules appartenant a la ct*
       * recherche du noeud air le plus proche      *
       *--------------------------------------------*/
      dmin = 1000. ;
      iair = location_cel[iwat] -1;

      /*--------------------------------------------*
       * initialiation matrice d interpolation et   *
       *  tableau nvois[] et coeff[]          *
       *--------------------------------------------*/

      for (ii=0;ii<4;ii++){
        for (jj=0;jj<4;jj++) {
          pp[jj][ii]   = 0.0;
          ppInv[jj][ii]= 0.0;
        }
      }

      for (jj=0;jj< cs_ctwr_nmaxvoi;jj++) {
        coeff[jj]= -1.0 ;
        nvois[jj]= -1  ;
      }
      /* fin initialisation */

      /*-------------------------------------------------*
       * Recherche des voisins du noeuds air le + proche *
       *  boucle sur les faces internes du maillage        *
       *-------------------------------------------------*/
      nbvois = 1 ;
      nvois[0] = iair ;

      for (ifac = 0 ; ifac < ct->nbfac_ct ; ifac++) {
        icel_1 = i_face_cells[ (lst_par_fac[ifac]- mesh->n_b_faces - 1) * 2    ]-1 ;
        icel_2 = i_face_cells[ (lst_par_fac[ifac]- mesh->n_b_faces - 1) * 2 + 1]-1 ;
        if (icel_1==iair) {
          nvois[nbvois]=icel_2 ;
          nbvois += 1 ;
        }
        else if (icel_2==iair) {
          nvois[nbvois]=icel_1 ;
          nbvois += 1 ;
        }
      }


      for (icel = cell_cells_idx[ iair     ];
            icel < cell_cells_idx[ iair + 1 ]; icel++) {

         indice = cell_cells_lst[ icel ] - 1;
         if (ct->icoul == family_item[ cell_family[ indice ]]){
            nvois[nbvois]= indice ;
            nbvois += 1 ;
         }
      }

      /* fin Recherche */

      /*--------------------------------------------*
       *nombre de voisins insuffisant               *
       *--------------------------------------------*/

      if (nbvois<nbn){
        nbvois = 1 ;
        nvois[0] = iair ;
        coeff[0] = 1.0 ;
        goto enregistre ;
      }

      dim = ct->idimct;
      vectBase[0][0] = 1.0; vectBase[1][0] = 0.0; vectBase[2][0] = 0.0;
      vectBase[0][1] = 0.0; vectBase[1][1] = 1.0; vectBase[2][1] = 0.0;
      vectBase[0][2] = 0.0; vectBase[1][2] = 0.0; vectBase[2][2] = 1.0;
      passage2D : ;


      /*--------------------------------------------*/
      /* Calcul de l'ouverture de la fonction  de pondération*/
      /*  egale au max de la distance entre le noeud eau et les voisins air */
      /*--------------------------------------------*/
      ouv = 0. ;
      for (ii = 0 ; ii < nbvois ; ii++) {
        iair = nvois[ii];
        dxx  =  (cs_real_t) (coo_cel[iair*3+0] - xwat );
        dyy  =  (cs_real_t) (coo_cel[iair*3+1] - ywat );
        dzz  =  (cs_real_t) (coo_cel[iair*3+2] - zwat );

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        aux =  pow(dx,2.)+pow(dy,2.)+pow(dz,2.) ;
        if (ouv<aux) {
          ouv = aux;
        }
      }
      ouv = sqrt(ouv)*1.1 ;
      /* fin calcul de l'ouverture */

      /*--------------------------------------------*/
      /*Construction de la matrice A               */
      /*--------------------------------------------*/
      for (ii = 0 ; ii < nbvois ; ii++) {

        indice = nvois[ii];
        dxx = (cs_real_t) (coo_cel[indice*3+0] - xwat) ;
        dyy = (cs_real_t) (coo_cel[indice*3+1] - ywat) ;
        dzz = (cs_real_t) (coo_cel[indice*3+2] - zwat) ;

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        /* parametre de la fonction de ponderation*/
        cx = 1.0 ;
        cy = 1.0 ;
        cz = 1.e10 ;
        lf = 3;
        epgauss = 5.0 ;
        if (dim==3) {
          cz = 1.0 ;
        }
        if (dim==1) {
          cy = 1.e10 ;
          cz = 1.e10 ;
        }
        /*fonction de ponderation*/
        w  = _weighting(dx,dy,dz,ouv,lf,epgauss,cx,cy,cz) ;
        if (dim == 1)/* 1D */{
          pp[0][0] = w      + pp[0][0] ;
          pp[0][1] = w*dx    + pp[0][1] ;
          pp[1][0] = w*dx    + pp[1][0] ;
          pp[1][1] = w*dx*dx + pp[1][1] ;
        }else
          if (dim == 2) /* 2D */{
            pp[0][0] = w       + pp[0][0] ;
            pp[0][1] = w*dx    + pp[0][1] ;
            pp[0][2] = w*dy    + pp[0][2] ;
            pp[1][0] = w*dx    + pp[1][0] ;
            pp[1][1] = w*dx*dx + pp[1][1] ;
            pp[1][2] = w*dx*dy + pp[1][2] ;
            pp[2][0] = w*dy    + pp[2][0] ;
            pp[2][1] = w*dx*dy + pp[2][1] ;
            pp[2][2] = w*dy*dy + pp[2][2] ;
          }
          else if (dim == 3)/* 3D */{
            pp[0][0] = w       + pp[0][0] ;
            pp[0][1] = w*dx    + pp[0][1] ;
            pp[0][2] = w*dy    + pp[0][2] ;
            pp[0][3] = w*dz    + pp[0][3] ;
            pp[1][0] = w*dx    + pp[1][0] ;
            pp[1][1] = w*dx*dx + pp[1][1] ;
            pp[1][2] = w*dy*dx + pp[1][2] ;
            pp[1][3] = w*dz*dx + pp[1][3] ;
            pp[2][0] = w*dy    + pp[2][0] ;
            pp[2][1] = w*dx*dy + pp[2][1] ;
            pp[2][2] = w*dy*dy + pp[2][2] ;
            pp[2][3] = w*dz*dy + pp[2][3] ;
            pp[3][0] = w*dz    + pp[3][0] ;
            pp[3][1] = w*dx*dz + pp[3][1] ;
            pp[3][2] = w*dy*dz + pp[3][2] ;
            pp[3][3] = w*dz*dz + pp[3][3] ;
        }
      }
      /*Fin Construction de la matrice A */
      if(ct->idimct == 3 && dim == 3){
          dim  = _is_coplanar(coo_cel, nvois, nbvois, pp, vectBase, 0, dhi);

            if(dim!= 3){
              for (ii=0;ii<3;ii++)
                for (jj=0;jj<3;jj++)
                  pp[jj][ii]=0.0;
              goto passage2D ;
            }

      }

      /*--------------------------------------------*/
      /* inversion de la matrice par la methode de  */
      /* jordan                                     */
      /*--------------------------------------------*/

      if (_invmat(pp, ppInv, dim) == 0) cs_exit(EXIT_FAILURE);

      /*--------------------------------------------*/
      /* Calcul des coefficients                    */
      /*--------------------------------------------*/
      for (ii = 0 ; ii < nbvois ; ii++){

        indice = nvois[ii];
        dxx = (cs_real_t) (coo_cel[indice*3+0] - xwat );
        dyy = (cs_real_t) (coo_cel[indice*3+1] - ywat );
        dzz = (cs_real_t) (coo_cel[indice*3+2] - zwat );

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        /* parametre de la fonction de ponderation*/
        cx = 1.0 ;
        cy = 1.0 ;
        cz = 1.e10 ;
        lf = 3;
        epgauss = 5.0 ;
        if (dim==3) {
          cz = 1.0 ;
        }
        if (dim==1) {
          cy = 1.e10 ;
          cz = 1.e10 ;
        }

        w  = _weighting(dx,dy,dz,ouv,lf,epgauss,cx,cy,cz) ;

        if (dim == 1){
          coeff[ii] = w*(ppInv[0][0]+ppInv[0][1]*dx);
        }
        else if (dim ==2){
          coeff[ii] = w*(ppInv[0][0]+ppInv[0][1]*dx+ppInv[0][2]*dy);
        }
        else if (dim ==3){
          coeff[ii] = w*(  ppInv[0][0]
                          +ppInv[0][1]*dx
                          +ppInv[0][2]*dy
                          +ppInv[0][3]*dz );
        }

      }
      /* Fin Calcul des coefficients */

      enregistre : ;

      /*--------------------------------------------*/
      /* boucle while sur pvoiseau pour trouver le  */
      /* dernier indice                             */
      /*--------------------------------------------*/
      indice = 0 ;
      while (ct->voiseau[indice]!=-1) {
        indice += 1 ;
      }
      /*--------------------------------------------*
       * Ajout des voisins et des coefficients      *
       *--------------------------------------------*/
      for (ii = 0 ; ii < nbvois ; ii++) {
        ct->voiseau[indice+ii]  = nvois[ii] ;
        ct->coefeau[indice+ii]  = coeff[ii];
      }

      ct->pvoiseau[iwat+1] = ct->pvoiseau[iwat] + nbvois ;

    }/* fin boucle sur iseg */
    BFT_FREE(lst_par_fac);
    BFT_FREE(lst_xyz_cel);
    BFT_FREE(lst_par_cel);

  } /* fin boucle sur les zones d echanges ict */

}

/*-----------------------------------------------------------------------------*
 * Function cs_ctwr_adair                                                        *
 * Interpolation EAU -> AIR                                                    *
 *-----------------------------------------------------------------------------*/
void cs_ctwr_adair
(
  const cs_real_t          gx,            /* composante x de la gravite */
  const cs_real_t          gy,            /* composante y de la gravite */
  const cs_real_t          gz             /* composante z de la gravite */
)
{
  /* Coordonnées des centres des cellules  */

  const fvm_coord_t  *lst_xyz_cel   = NULL;
  fvm_coord_t  *lst_xyz_water = NULL;
  const fvm_lnum_t   *location_cel  = NULL;
  cs_int_t   ict,icol, ilig,ieau,ii,jj,iair,nbvois,
             nbn,nvois[ cs_ctwr_nmaxvoi ],lf,indice;
  cs_int_t   dim, nb_air_node;
  cs_real_t  norme_g,dhi,dmin,dist,coeff[ cs_ctwr_nmaxvoi ],ouv,aux,gravite[3];
  cs_real_t  dx,dy,dz,dxx,dyy,dzz;
  cs_real_t  cx,cy,cz,epgauss,w;
  cs_real_t  pp[4][4], ppInv[4][4] ;
  cs_real_t  vectBase[3][3];
  fvm_lnum_t loca_cel;
  cs_ctwr_zone_t  *ct;
  /*--------------------------------------------*
   * parametres et initialisation               *
   *--------------------------------------------*/
  /* Vecteur gravite */
  gravite[0] = gx;
  gravite[1] = gy;
  gravite[2] = gz;

  norme_g = sqrt( pow(gravite[0],2.)
                  +pow(gravite[1],2.)
                  +pow(gravite[2],2.) );



  /* fin parametres et initialisation */

  /*---------------------------------------------*
   * Construction des coefficient d'interpolation*
   * sur chaque zone d echange ict               *
   *---------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {

    ct = cs_glob_ct_tab[ict];

    nbn = 3 ;
    if (ct->idimct==3) {
      nbn = 4 ;
    }


    /* Calcul de dh */
    dhi = (ct->hmax-ct->hmin)/(ct->nelect-1);

    /* Copy element centers of the water mesh to an array.*/


    nb_air_node = (int)fvm_locator_get_n_dist_points(ct->locat_water_air);

    lst_xyz_cel = fvm_locator_get_dist_coords(ct->locat_water_air);

    location_cel = fvm_locator_get_dist_locations(ct->locat_water_air);


    BFT_MALLOC(lst_xyz_water, (3*ct->nelect*ct->nnpsct) ,fvm_coord_t);
    fvm_nodal_get_element_centers(ct->water_mesh,FVM_INTERLACE,3, lst_xyz_water );

    if(ct->n_ghost_npsct >0  ) {
      BFT_REALLOC(lst_xyz_water, (3*ct->nelect*ct->nnpsct_with_ghosts) ,fvm_coord_t);
      cs_halo_sync_var_strided(ct->water_halo, ct->halo_type, lst_xyz_water, 3);

    }
    /*--------------------------------------------*
     * Loops on the air nodes of teh exchange area*
     *--------------------------------------------*/

    for (iair = 0 ; iair < nb_air_node ; iair++)  {

       loca_cel = location_cel[iair] -1;
      /*--------------------------------------------*
       * initialiation matrice d interpolation et   *
       * tableau nvois[] et coeff[]                 *
       *--------------------------------------------*/

      for (ii=0;ii<4;ii++) {
        for (jj=0;jj<4;jj++) {
          pp[jj][ii]=0.0;
          ppInv[jj][ii]=0.0;
        }
      }

      for (jj=0;jj< cs_ctwr_nmaxvoi;jj++) {
        coeff[jj]=0.0 ;
        nvois[jj]= -1  ;
      }/* fin initialisation */

      /*--------------------------------------------*
       * Traitement particulier pour les noeuds air *
       * en bordure inferieure ou superieure des ct *
       *--------------------------------------------*/

      /* indice du noeud air dans le maillage */
      dmin   = 1000. ;

        if ((loca_cel%(ct->nelect) == 0 ) ||
                  (loca_cel% (ct->nelect) == (ct->nelect-1) ))  {
          for (jj = 0 ; jj < (ct->nelect*ct->nnpsct_with_ghosts); jj++) {
            dx = (cs_real_t) (lst_xyz_water[3*jj  ] - lst_xyz_cel[iair*3  ]) ;
            dy = (cs_real_t) (lst_xyz_water[3*jj+1] - lst_xyz_cel[iair*3+1]) ;
            dz = (cs_real_t) (lst_xyz_water[3*jj+2] - lst_xyz_cel[iair*3+2]) ;
            dist = (pow(dx,2.)+pow(dy,2.)+pow(dz,2.)) ;
            if (dmin>dist ) {
              dmin = dist ;
              ieau = jj ;
            }
          }
        }

      if (dmin<1000.) {
        /* Cellule air est en bordure inf ou sup, on saute l'etape suivante */
        nbvois   = 1 ;
        nvois[0] = ieau ;
        coeff[0] = 1.0 ;
        goto enregistre ;
      }
      /*------------------------------------------------*
       * Fin Traitement particulier pour les noeuds air *
       * en bordure inferieure ou superieure des ct     *
       *------------------------------------------------*/

      /*--------------------------------------------*
       * On continue avec les noeuds air qui ne sont*
       * pas en bordure                             *
       * recherche des cellules air voisines pour   *
       * les noeuds air qui ne sont                 *
       *  pas en bordure inferieure ou superieure   *
       *--------------------------------------------*/
      nbvois = 1 ;
      nvois[0] = loca_cel;
      /*---------------------------------------------*
       * Recherche du nombre de voisins du noeuds air*
       * boucle sur les faces internes du maillage   *
       *---------------------------------------------*/


      indice=1;
      nvois[indice++] = loca_cel + 1;
      nvois[indice++] = loca_cel - 1;
      nbvois+=2;
      icol = loca_cel/(ct->nelect);
      ilig = loca_cel%(ct->nelect);
      for (ii = ct->fac_sup_connect_idx[ icol ] ;
           ii < ct->fac_sup_connect_idx[ icol + 1 ] ; ii++) {

           nvois[indice++] = ct->nelect*ct->fac_sup_connect_lst[ ii ] + ilig-1;
           nvois[indice++] = ct->nelect*ct->fac_sup_connect_lst[ ii ] + ilig  ;
           nvois[indice++] = ct->nelect*ct->fac_sup_connect_lst[ ii ] + ilig+1;
           nbvois += 3;

      }

      /*---------------------------------------------*
       * nombre de voisin eau insuffisants           *
       * meme que noeuds en bordure                  *
       *-------------------------------------------- */

      if (nbvois<nbn) {
        nbvois   = 1 ;
        coeff[0] = 1.0 ;
        goto enregistre ;
      }

       dim = ct->idimct;
       vectBase[0][0] = 1.0; vectBase[1][0] = 0.0; vectBase[2][0] = 0.0;
       vectBase[0][1] = 0.0; vectBase[1][1] = 1.0; vectBase[2][1] = 0.0;
       vectBase[0][2] = 0.0; vectBase[1][2] = 0.0; vectBase[2][2] = 1.0;

       passage2D : ;

      /*--------------------------------------------*
       * Calcul de l'ouverture de la fonction  de   *
       * pondération egale au max de la distance    *
       *  entre le noeud air et les voisins eau     *
       *--------------------------------------------*/
      ouv = 0. ;
      for (ii = 0 ; ii < nbvois ; ii++) {
        ieau = nvois[ii] ;
        dxx = (cs_real_t) (lst_xyz_water[3*ieau]   - lst_xyz_cel[iair*3+0]);
        dyy = (cs_real_t) (lst_xyz_water[3*ieau+1] - lst_xyz_cel[iair*3+1]);
        dzz = (cs_real_t) (lst_xyz_water[3*ieau+2] - lst_xyz_cel[iair*3+2]);

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        aux =  pow(dx,2.)+pow(dy,2.)+pow(dz,2.) ;
        if (ouv<aux) {
          ouv = aux;
        }
      }
      ouv = sqrt(ouv)*1.1 ;
      /* Fin de calcul de l'ouverture */

      /*--------------------------------------------*
       *Construction de la matrice A                *
       *--------------------------------------------*/

      for (ii = 0 ; ii < nbvois ; ii++) {

        ieau = nvois[ii] ;
        dxx = (cs_real_t) (lst_xyz_water[3*ieau +0] - lst_xyz_cel[iair*3+0]) ;
        dyy = (cs_real_t) (lst_xyz_water[3*ieau +1] - lst_xyz_cel[iair*3+1]) ;
        dzz = (cs_real_t) (lst_xyz_water[3*ieau +2] - lst_xyz_cel[iair*3+2]) ;

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];
        /* parametre de la fonction de ponderation*/
        cx = 1.0 ;
        cy = 1.0 ;
        cz = 1.e10 ;
        lf = 3;
        epgauss = 5.0 ;
        if (dim==3) {
          cz = 1.0 ;
        }
        if (dim==1) {
          cy = 1.e10 ;
          cz = 1.e10 ;
        }
        /*fonction de ponderation*/
        w  = _weighting(dx,dy,dz,ouv,lf,epgauss,cx,cy,cz) ;

        if (dim == 1)/* 1D */{
            pp[0][0] = w       + pp[0][0] ;
            pp[0][1] = w*dx    + pp[0][1] ;
            pp[1][0] = w*dx    + pp[1][0] ;
            pp[1][1] = w*dx*dx + pp[1][1] ;
        }else if (dim == 2)/* 2D */{
                  pp[0][0] = w     + pp[0][0] ;
                  pp[0][1] = w*dx    + pp[0][1] ;
                  pp[0][2] = w*dy    + pp[0][2] ;
                  pp[1][0] = w*dx    + pp[1][0] ;
                  pp[1][1] = w*dx*dx + pp[1][1] ;
                  pp[1][2] = w*dx*dy + pp[1][2] ;
                  pp[2][0] = w*dy    + pp[2][0] ;
                  pp[2][1] = w*dx*dy + pp[2][1] ;
                  pp[2][2] = w*dy*dy + pp[2][2] ;
             }else if (dim == 3)/* 3D */{
                        pp[0][0] = w       + pp[0][0] ;
                        pp[0][1] = w*dx    + pp[0][1] ;
                        pp[0][2] = w*dy    + pp[0][2] ;
                        pp[0][3] = w*dz    + pp[0][3] ;
                        pp[1][0] = w*dx    + pp[1][0] ;
                        pp[1][1] = w*dx*dx + pp[1][1] ;
                        pp[1][2] = w*dy*dx + pp[1][2] ;
                        pp[1][3] = w*dz*dx + pp[1][3] ;
                        pp[2][0] = w*dy    + pp[2][0] ;
                        pp[2][1] = w*dx*dy + pp[2][1] ;
                        pp[2][2] = w*dy*dy + pp[2][2] ;
                        pp[2][3] = w*dz*dy + pp[2][3] ;
                        pp[3][0] = w*dz    + pp[3][0] ;
                        pp[3][1] = w*dx*dz + pp[3][1] ;
                        pp[3][2] = w*dy*dz + pp[3][2] ;
                        pp[3][3] = w*dz*dz + pp[3][3] ;

             }

      }/* Fin de construction de la matrice A*/
      if(ct->idimct == 3 && dim == 3){
        dim  = _is_coplanar(lst_xyz_water, nvois, nbvois, pp, vectBase,0,dhi);
        if(dim!= 3) {
          for (ii=0;ii<3;ii++)
            for (jj=0;jj<3;jj++)
              pp[jj][ii]=0.;
          goto passage2D ;
        }
      }

      /*--------------------------------------------*
       * inversion de la matrice par la methode de  *
       * jordan                                     *
       *--------------------------------------------*/
      if (_invmat(pp, ppInv, dim) == 0) cs_exit(EXIT_FAILURE);

      /*--------------------------------------------*
       * Calcul des coefficients                    *
       *--------------------------------------------*/
      for (ii = 0 ; ii < nbvois ; ii++) {
        ieau = nvois[ii] ;
        dxx = (cs_real_t) (lst_xyz_water[3*ieau   ] - lst_xyz_cel[iair*3+0]);
        dyy = (cs_real_t) (lst_xyz_water[3*ieau +1] - lst_xyz_cel[iair*3+1]);
        dzz = (cs_real_t) (lst_xyz_water[3*ieau +2] - lst_xyz_cel[iair*3+2]);

        dx = dxx * vectBase[0][0] + dyy * vectBase[0][1] + dzz * vectBase[0][2];
        dy = dxx * vectBase[1][0] + dyy * vectBase[1][1] + dzz * vectBase[1][2];
        dz = dxx * vectBase[2][0] + dyy * vectBase[2][1] + dzz * vectBase[2][2];

        /*parametre de la fonction de ponderation*/
        cx = 1.0 ;
        cy = 1.0 ;
        cz = 1.e10 ;
        lf = 3;
        epgauss = 5.0 ;
        if (dim==3) {
          cz = 1.0 ;
        }
        if (dim==1) {
          cy = 1.e10 ;
          cz = 1.e10 ;
        }

        w = _weighting(dx,dy,dz,ouv,lf,epgauss,cx,cy,cz) ;



        if (dim == 1){
          coeff[ii] = w*(ppInv[0][0]+ppInv[0][1]*dx);
        }
        else if (dim == 2){
          coeff[ii] = w*(ppInv[0][0]+ppInv[0][1]*dx+ppInv[0][2]*dy);
        }
        else if (dim == 3){
          coeff[ii] = w*(ppInv[0][0]
                         +ppInv[0][1]*dx
                         +ppInv[0][2]*dy
                         +ppInv[0][3]*dz );
        }

      }
      /* Fin de calcul des coefficients */

      /*---------------------------------------------*
       * note :Reprise pour les noeuds air en bordure*
       * ou avec un nbre de voisin insuffisant       *
       *---------------------------------------------*/
      enregistre : ;

      /*--------------------------------------------*
       * trouver le dernier indice sur pvoisair     *
       *--------------------------------------------*/
      indice = 0 ;
      while (ct->voisair[indice]!=-1) {
        indice += 1 ;
      }

      /*--------------------------------------------*
       * Ajout des voisins et des coefficients      *
       *--------------------------------------------*/
      for (icol = 0 ; icol < nbvois ; icol++)
        {
          ct->voisair[indice + icol]  = nvois[ icol ];
          ct->coefair[indice + icol]  = coeff[ icol ];
         }


      ct->pvoisair[iair+1] = ct->pvoisair[iair] + nbvois ;
    }
    /*---------------------------------------------*
     * fin de la boucle sur les noeuds air de la ct*
     *---------------------------------------------*/

    BFT_FREE(lst_xyz_water);


  }/* fin de la boucle sur les ct */

}


/*----------------------------------------------------------------------------*
 * Chaining of the exchange area                                              *
 *----------------------------------------------------------------------------*/

void
cs_ctwr_stacking(const cs_real_t  gx,
                 const cs_real_t  gy,
                 const cs_real_t  gz)
{
  cs_int_t i, j, rank, dist_rank, nb, nb_ct, itmp, ict, ict_uw;
  cs_int_t * aux;
  cs_ctwr_zone_t  *ct, *ct_upw;
  cs_real_t tmp;
  cs_real_t gravite[3];
  fvm_coord_t * lst_xyz;
  const double tolerance = 0.1;

  nb = cs_glob_ct_nbr  * cs_glob_ct_nbr;

  BFT_MALLOC(cs_stack_ct, nb, cs_int_t);
  BFT_MALLOC(cs_chain_ct, cs_glob_ct_nbr, cs_int_t);

  gravite[0]= gx;
  gravite[1]= gy;
  gravite[2]= gz;

  for (i=0 ; i < cs_glob_ct_nbr ; i++)
    for (j=0 ; j < cs_glob_ct_nbr ; j++)
      cs_stack_ct[i*cs_glob_ct_nbr + j]=0;

  for (i=0 ; i < cs_glob_ct_nbr; i++)
    for (j=0 ; j < cs_glob_ct_nbr ; j++)
      if (CS_ABS(cs_glob_ct_tab[i]->hmax - cs_glob_ct_tab[j]->hmin)< 1.e-6)
        cs_stack_ct[i*cs_glob_ct_nbr + j] =1;

#if defined(HAVE_MPI)

  if (cs_glob_n_ranks > 1) {

    BFT_MALLOC(aux, nb, cs_int_t);
    rank = cs_glob_rank_id;

    for (dist_rank = 0; dist_rank < cs_glob_n_ranks; dist_rank++)
      if (dist_rank != rank){

        MPI_Sendrecv(cs_stack_ct, nb, CS_MPI_INT , dist_rank, CS_CT_MPI_TAG,
                     aux, nb, CS_MPI_INT, dist_rank, CS_CT_MPI_TAG,
                     cs_glob_mpi_comm, &status);
        for (i=0 ; i < cs_glob_ct_nbr ; i++)
          for (j=0 ; j < cs_glob_ct_nbr ; j++){
            if (aux[i*cs_glob_ct_nbr + j] > cs_stack_ct[i*cs_glob_ct_nbr + j])
              cs_stack_ct[i*cs_glob_ct_nbr + j] = aux[i*cs_glob_ct_nbr + j];
          }
      }

    BFT_FREE(aux);

  }
#endif

  /* to order the exchange area */
    /*Init the chaining array */
  for(i = 0; i < cs_glob_ct_nbr ; i++)
    cs_chain_ct[i] = i;

  for (i = 0; i < cs_glob_ct_nbr ; i++)
    for (j = i+1; j < cs_glob_ct_nbr ; j++)
      if (cs_stack_ct[cs_chain_ct[i]*cs_glob_ct_nbr + cs_chain_ct[j]] == 1 ){
        itmp = cs_chain_ct [i];
        cs_chain_ct [i] = cs_chain_ct [j];
        cs_chain_ct [j] = itmp;
      }

  for(ict = 0; ict< cs_glob_ct_nbr ; ict++){

    ct = cs_glob_ct_tab[cs_chain_ct[ict]];
    nb_ct = 0;

    for (ict_uw = 0 ; ict_uw < cs_glob_ct_nbr ; ict_uw++)
      if (cs_stack_ct[cs_chain_ct[ict]*cs_glob_ct_nbr + cs_chain_ct[ict_uw]]
          == 1){

        nb_ct++;
        ct_upw = cs_glob_ct_tab[cs_chain_ct[ict_uw]];

        BFT_MALLOC( lst_xyz ,
                    3*(ct_upw->nbfac_ict+ct_upw->nbfbr_ict) ,fvm_coord_t);

        fvm_nodal_get_element_centers
                    ( ct_upw->face_inf_mesh,FVM_INTERLACE,2,lst_xyz );

        tmp  = CS_ABS(ct_upw->hmax - ct_upw->hmin)/(ct_upw->nelect-1);
        tmp /= CS_LOC_MODULE(gravite);

        for (i=0 ; i < (ct_upw->nbfac_ict+ct_upw->nbfbr_ict) ; i++){
          lst_xyz[3*i + 0 ] -= tmp * gravite[0];
          lst_xyz[3*i + 1 ] -= tmp * gravite[1];
          lst_xyz[3*i + 2 ] -= tmp * gravite[2];
        }

        BFT_REALLOC(ct->locat_cell_ct_upwind, nb_ct, fvm_locator_t *);

#if defined(FVM_HAVE_MPI)
        ct->locat_cell_ct_upwind[nb_ct-1] =
          fvm_locator_create(tolerance,
                             cs_glob_mpi_comm,
                             cs_glob_n_ranks,
                             0);
#else
        ct->locat_cell_ct_upwind[nb_ct-1] = fvm_locator_create(tolerance);
#endif

        fvm_locator_set_nodal(ct->locat_cell_ct_upwind[nb_ct-1],
                              ct_upw->water_mesh,
                              0,
                              3,
                              ct_upw->nbfac_ict+ct_upw->nbfbr_ict,
                              NULL,
                              lst_xyz);
        BFT_FREE(lst_xyz);

      }
  }
}

/*----------------------------------------------------------------------------
 * Destroy cs_ctwr_t structures
 *----------------------------------------------------------------------------*/

void
cs_ctwr_all_destroy(void)
{
  int i;
  cs_ctwr_zone_t  *ct;

  for (i = 0 ; i < cs_glob_ct_nbr; i++) {

    ct = cs_glob_ct_tab[i];
    BFT_FREE(ct);

  }

  cs_glob_ct_nbr_max = 0;
  cs_glob_ct_nbr = 0;

  BFT_FREE(cs_stack_ct);
  BFT_FREE(cs_chain_ct);


  BFT_FREE(cs_glob_ct_tab);
}


/*----------------------------------------------------------------------------
 * Water variables resolution
 *----------------------------------------------------------------------------*/

void
cs_ctwr_aeteau(cs_real_t   temp[],      /* Temperature air */
               cs_real_t   xa[],        /* humidite air */
               cs_real_t   rho[],       /* masse volumique air */
                 cs_real_t   vitx[],      /* vitesse air suivant x */
               cs_real_t   vity[],      /* vitesse air suivant y */
               cs_real_t   vitz[],      /* vitesse air suivant z */
               const cs_real_t   gx,          /* composante x de la gravite */
               const cs_real_t   gy,          /* composante y de la gravite */
               const cs_real_t   gz           /* composante z de la gravite */

)
{
  cs_int_t  ict,iseg,iloc,i, ii,j ,ieau, nb_dist_water, nb_dist_upw, ind;
  cs_real_t dhi,vvai,vhai,norme_g;
  cs_real_t gravite[3];
  cs_real_t faxn,bxan,xsata,xsate,cfen,ff1,ff2,xlew,eta,aux;
  cs_real_t vgin,dvg,cpx,rre,rpr,anu;
  cs_real_t   cpe, cpv, cpa, hv0, dgout, visc, conduc, rhoe;


  fvm_lnum_t *lst_par_fac_sup_ct, *lst_par_fac_inf_ct_upw;
  const fvm_lnum_t *locat_cel_upw = NULL;

  cs_ctwr_zone_t  *ct;
  cs_ctwr_zone_t  *ct_upw;
  cs_real_t *tai_inter, *xai_inter, *rhoai_inter,*vx_inter, *vy_inter,*vz_inter;
  cs_real_t *tai, *xai, *rhoai,*vx, *vy, *vz, *teau_upw_rec, *teau_upw_send ;
  cs_real_t *fem_upw_rec, *fem_upw_send;

  gravite[0] = -gx;
  gravite[1] = -gy;
  gravite[2] = -gz;

  norme_g = sqrt( pow(gravite[0],2.)
                  +pow(gravite[1],2.)
                  +pow(gravite[2],2.) );

  gravite[0] /= norme_g;
  gravite[1] /= norme_g;
  gravite[2] /= norme_g;


  /*--------------------------------------------*
   * Résolution des variable eau                *
   * sur chaque  ct                             *
   *--------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {

    ct = cs_glob_ct_tab[cs_chain_ct[ict]];

    if ((ct->ntypct>=2) && ( ct->idimct==2) )
     gravite[2] = 1.0;

    cpa    = ct->cpa ;
    cpv    = ct->cpv ;
    cpe    = ct->cpe ;
    hv0    = ct->hv0 ;
    rhoe   = ct->rhoe ;
    dgout  = ct->dgout ;
    visc   = ct->visc ;
    conduc = ct->conduc  ;

    /*--------------------------------------------*
     * synchronisation   Halo                             *
     *--------------------------------------------*/

    if( ct->n_ghost_npsct >0 ){

      cs_halo_t *halo = ct->water_halo;

      cs_halo_sync_var(halo, ct->halo_type, temp);
      cs_halo_sync_var(halo, ct->halo_type, xa);
      cs_halo_sync_var(halo, ct->halo_type, rho);
      cs_halo_sync_var(halo, ct->halo_type, vitx);
      cs_halo_sync_var(halo, ct->halo_type, vity);
      cs_halo_sync_var(halo, ct->halo_type, vitz);

    }

    /*--------------------------------------------*
     * interpolation  air->eau                    *
     *--------------------------------------------*/
    nb_dist_water = (int) fvm_locator_get_n_dist_points(ct->locat_air_water);

    BFT_MALLOC( tai_inter  , nb_dist_water, cs_real_t);
    BFT_MALLOC( xai_inter  , nb_dist_water, cs_real_t);
    BFT_MALLOC( rhoai_inter, nb_dist_water, cs_real_t);
    BFT_MALLOC( vx_inter   , nb_dist_water, cs_real_t);
    BFT_MALLOC( vy_inter   , nb_dist_water, cs_real_t);
    BFT_MALLOC( vz_inter   , nb_dist_water, cs_real_t);

    for (ieau= 0 ; ieau < nb_dist_water ; ieau++) {
       tai_inter[ieau]   = 0.;
       xai_inter[ieau]   = 0.;
       rhoai_inter[ieau] = 0.;
       vx_inter[ieau]    = 0.;
       vy_inter[ieau]    = 0.;
       vz_inter[ieau]    = 0.;
      for (i = (ct->pvoiseau[ieau]) ; i < (ct->pvoiseau[ieau+1]) ; i++) {
        tai_inter[ieau]  += ct->coefeau[i] * temp[ct->voiseau[i]];
        xai_inter[ieau]  += ct->coefeau[i] * xa[ct->voiseau[i]];
        rhoai_inter[ieau]+= ct->coefeau[i] * rho[ct->voiseau[i]];
        vx_inter[ieau]   += ct->coefeau[i] * vitx[ct->voiseau[i]];
        vy_inter[ieau]   += ct->coefeau[i] * vity[ct->voiseau[i]];
        vz_inter[ieau]   += ct->coefeau[i] * vitz[ct->voiseau[i]];
      }
    }
    BFT_MALLOC( tai  , ct->nnpsct*ct->nelect, cs_real_t );
    BFT_MALLOC( xai  , ct->nnpsct*ct->nelect, cs_real_t );
    BFT_MALLOC( rhoai, ct->nnpsct*ct->nelect, cs_real_t );
    BFT_MALLOC( vx   , ct->nnpsct*ct->nelect, cs_real_t );
    BFT_MALLOC( vy   , ct->nnpsct*ct->nelect, cs_real_t );
    BFT_MALLOC( vz   , ct->nnpsct*ct->nelect, cs_real_t );

    fvm_locator_exchange_point_var(ct->locat_air_water,
                                   tai_inter, tai, NULL, sizeof(cs_real_t),1,0);
    fvm_locator_exchange_point_var(ct->locat_air_water,
                                   xai_inter, xai, NULL, sizeof(cs_real_t),1,0);
    fvm_locator_exchange_point_var(ct->locat_air_water,
                                   rhoai_inter,rhoai, NULL, sizeof(cs_real_t),1,0);
    fvm_locator_exchange_point_var(ct->locat_air_water,
                                   vx_inter,vx, NULL, sizeof(cs_real_t),1,0);
    fvm_locator_exchange_point_var(ct->locat_air_water,
                                   vy_inter,vy, NULL, sizeof(cs_real_t),1,0);
    fvm_locator_exchange_point_var(ct->locat_air_water,
                                   vz_inter,vz, NULL, sizeof(cs_real_t),1,0);

    BFT_FREE( tai_inter  );
    BFT_FREE( xai_inter  );
    BFT_FREE( rhoai_inter);
    BFT_FREE( vx_inter  );
    BFT_FREE( vy_inter  );
    BFT_FREE( vz_inter  );
    /*--------------------------------------------*
     *  end interpolation  air->eau               *
     *--------------------------------------------*/



   /*--------------------------------------------*
    * Calcul pour la face superieure ,           *
    * Introduction des conditions aux limites ct *
    *--------------------------------------------*/
    BFT_MALLOC( lst_par_fac_sup_ct , ct->nnpsct, fvm_lnum_t );

    fvm_nodal_get_parent_num(ct->face_sup_mesh,
                                      2,lst_par_fac_sup_ct);
    ind = 0;
    for (j=0 ; j < cs_glob_ct_nbr ; j++)
      if(cs_stack_ct[cs_chain_ct[ict]*cs_glob_ct_nbr + cs_chain_ct[j]] == 1){
        ct_upw = cs_glob_ct_tab[ cs_chain_ct[j]];

        nb_dist_upw =
              (int)fvm_locator_get_n_dist_points(ct->locat_cell_ct_upwind[ind]);

        BFT_MALLOC( teau_upw_send  , nb_dist_upw, cs_real_t);
        BFT_MALLOC( fem_upw_send  , nb_dist_upw, cs_real_t);
        BFT_MALLOC( lst_par_fac_inf_ct_upw , (ct_upw->nbfac_ict+ct_upw->nbfbr_ict), fvm_lnum_t );

        fvm_nodal_get_parent_num(ct_upw->face_inf_mesh,
                                      2,lst_par_fac_inf_ct_upw);
        locat_cel_upw =
                  fvm_locator_get_dist_locations(ct->locat_cell_ct_upwind[ind]);

        for (i=0 ; i < nb_dist_upw ; i++){
          teau_upw_send[i] =  ct_upw->teau[(cs_int_t) locat_cel_upw[i]-1];
          fem_upw_send[i]  =  ct_upw->fem[(cs_int_t) locat_cel_upw[i]-1];
        }

        BFT_MALLOC( teau_upw_rec,
                   (ct_upw->nbfac_ict+ct_upw->nbfbr_ict), cs_real_t );
        BFT_MALLOC( fem_upw_rec,
                   (ct_upw->nbfac_ict+ct_upw->nbfbr_ict), cs_real_t );

        fvm_locator_exchange_point_var(ct->locat_cell_ct_upwind[ind],
                                       teau_upw_send,
                                       teau_upw_rec,
                                       NULL,
                                       sizeof(cs_real_t),
                                       1,0);
        fvm_locator_exchange_point_var(ct->locat_cell_ct_upwind[ind],
                                       fem_upw_send,
                                       fem_upw_rec,
                                       NULL,
                                       sizeof(cs_real_t),
                                       1,0);

        for (i=0 ; i < ct->nnpsct ; i++){
          ii = 0;
          while (ii < (ct_upw->nbfac_ict+ct_upw->nbfbr_ict) ){
            if( lst_par_fac_sup_ct[i] == lst_par_fac_inf_ct_upw[ii]){
              ct->teau[i*ct->nelect] = teau_upw_rec[ii];
              ct->fem[i*ct->nelect]  = fem_upw_rec[ii];
              ii = ct_upw->nbfac_ict+ct_upw->nbfbr_ict;
            }
            ii++;
          }
        }
        BFT_FREE( teau_upw_rec );
        BFT_FREE( teau_upw_send );
        BFT_FREE( fem_upw_rec );
        BFT_FREE( fem_upw_send );
        BFT_FREE( lst_par_fac_inf_ct_upw );
        ind++;
      }

    BFT_FREE( lst_par_fac_sup_ct );

    /* Pas d'espace */
    dhi = -(ct->hmax-ct->hmin)/(ct->nelect-1);

    /*--------------------------------------------*/
    /* Modele de Poppe                            */
    /*--------------------------------------------*/
    if (ct->imctch==1) {
      /*--------------------------------------------*/
      /* courant-croise ou contre-courant           */
      /*--------------------------------------------*/
      if (ct->ntypct<=2) {

        for (iseg = 0 ; iseg < ct->nnpsct ; iseg++) {
          /*--------------------------------------------*/
          /* Resolution Fe                              */
          /*--------------------------------------------*/
          for (iloc = 1 ; iloc < ct->nelect ; iloc++) {

            ieau = iseg*ct->nelect + iloc ;


            vvai = sqrt(  pow((vx[ieau]*gravite[0]),2.)
                         +pow((vy[ieau]*gravite[1]),2.)
                         +pow((vz[ieau]*gravite[2]),2.));
            vhai = sqrt(pow((vx[ieau]*(1.-gravite[0])),2.)
                       +pow((vy[ieau]*(1.-gravite[1])),2.)
                       +pow((vz[ieau]*(1.-gravite[2])),2.));
            /* fin interpolation air->eau */

            xsate = cs_ctwr_xsath(ct->teau[ieau]);
            xsata = cs_ctwr_xsath(tai[ieau]);
            if (ct->ntypct==1) {
              faxn=rhoai[ieau]*vvai;
            }
            if (ct->ntypct==2) {
              faxn=rhoai[ieau]*vhai;
            }
            bxan=ct->xap*ct->fem[ieau]*pow((faxn/ct->fem[ieau]),ct->xnp);
            if ( xai[ieau]>xsata ) {
              aux=xsata;
            }else{
              aux=xai[ieau];
            }
            cfen=bxan*(xsate- aux )/(ct->fem[ieau]);
            ct->fem[ieau]=ct->fem[ieau-1]/(1.0-cfen*dhi);
          }
          /* Fin de résolution de Fe */

          /*--------------------------------------------*/
          /* Resolution Te                              */
          /*--------------------------------------------*/
          for (iloc = 1 ; iloc < ct->nelect ; iloc++) {

            ieau = iseg*ct->nelect + iloc ;

            vvai = sqrt( pow((vx[ieau]*gravite[0]),2.)
                        +pow((vy[ieau]*gravite[1]),2.)
                        +pow((vz[ieau]*gravite[2]),2.));
            vhai = sqrt( pow((vx[ieau]*(1.-gravite[0])),2.)
                        +pow((vy[ieau]*(1.-gravite[1])),2.)
                        +pow((vz[ieau]*(1.-gravite[2])),2.));
            /* fin interpolation air->eau */

            xsate = cs_ctwr_xsath(ct->teau[ieau]);
            xsata = cs_ctwr_xsath(tai[ieau]);
            if (ct->ntypct==1) {
              faxn=rhoai[ieau]*vvai;
            }
            if (ct->ntypct==2) {
              faxn=rhoai[ieau]*vhai;
            }
            bxan = ct->xap*ct->fem[ieau]*pow((faxn/ct->fem[ieau]),ct->xnp);
            cfen = bxan/ct->fem[ieau]/cpe;
            eta  = (0.622+xsate)/(0.622+xai[ieau]);
            xlew = pow(0.866,(2./3.))*(eta-1.)/log(eta);
            if (xai[ieau]<=xsata) {
              ff1 = (cpa+cpv*xsate)+(cpa+cpv*xai[ieau])*(xlew-1.);
              ff2 = xlew*(cpa+cpv*xai[ieau])*tai[ieau] -(xsate-xai[ieau])*hv0;
            }
            else {
              ff1 = xlew*(cpa+cpv*xsata)+(xsate-xsata)*(cpv-cpe);
              ff2 = xlew*(cpa+cpv*xsata+cpe*(xai[ieau]-xsata))*tai[ieau]
                    -(xsate-xsata)*(cpe*tai[ieau]+hv0);
            }
            ct->teau[ieau] = (ct->teau[ieau-1]-cfen*ff2*dhi)
                             /(1.0-cfen*ff1*dhi);

          }
          /* Fin de resolution Te  */


        }
      }
      /* Fin courant-croise ou contre-courant  */


      /*--------------------------------------------*/
      /* zone de pluie                              */
      /*--------------------------------------------*/
      else if (ct->ntypct==3){
        for (iseg = 0 ; iseg < ct->nnpsct ; iseg++) {
          /*--------------------------------------------*/
          /* Resolution Fe                              */
          /*--------------------------------------------*/
          for (iloc = 1 ; iloc < ct->nelect ; iloc++) {

            ieau = iseg*ct->nelect + iloc ;
            vgin=ct->vgoutte[ieau];

            if (CS_ABS(vgin)>=0.1){

              vvai = sqrt( pow((vx[ieau]*gravite[0]),2.)
                          +pow((vy[ieau]*gravite[1]),2.)
                          +pow((vz[ieau]*gravite[2]),2.));
              vhai = sqrt( pow((vx[ieau]*(1.-gravite[0])),2.)
                          +pow((vy[ieau]*(1.-gravite[1])),2.)
                          +pow((vz[ieau]*(1.-gravite[2])),2.));
              /* fin interpolation air->eau */

              xsate = cs_ctwr_xsath(ct->teau[ieau]);
              xsata = cs_ctwr_xsath(tai[ieau]);
              dvg=sqrt(pow((vgin+vvai),2.)+pow(vhai,2.));
              if (xai[ieau]<=xsata) {
                cpx = cpa+xai[ieau]*cpv;
              }
              else {
                cpx = cpa+xsata*cpv+(xai[ieau]-xsata)*cpe;
              }
              rre  = dvg*rhoai[ieau]*(1.+xai[ieau])*dgout/visc;
              rpr  = cpx*visc/conduc;
              anu  = 2.+0.6*sqrt(rre)*pow(rpr,(1./3.));
              bxan = 6.*conduc*anu*ct->fem[ieau]
                    /(0.92*rhoe*vgin*pow(dgout,2.)*cpx);
              if (xai[ieau]>xsata ) {
                aux = xsata;
              }else{
                aux = xai[ieau];
              }

              cfen=bxan*(xsate - aux)/(ct->fem[ieau]);
              ct->fem[ieau]=ct->fem[ieau-1]/(1.0-cfen*dhi);
            }
            else {
              ct->fem[ieau]=ct->fem[ieau-1];
            }
          }
          /* Fin de resolution Fe */

          /*--------------------------------------------*/
          /* Resolution Te                              */
          /*--------------------------------------------*/
          for (iloc = 1 ; iloc < ct->nelect ; iloc++) {

            ieau = iseg*ct->nelect + iloc ;

            vgin=ct->vgoutte[ieau];

            if (CS_ABS(vgin)>=0.1){

              vvai = sqrt(pow((vx[ieau]*gravite[0]),2.)
                         +pow((vy[ieau]*gravite[1]),2.)
                         +pow((vz[ieau]*gravite[2]),2.));
              vhai = sqrt(pow((vx[ieau]*(1.-gravite[0])),2.)
                         +pow((vy[ieau]*(1.-gravite[1])),2.)
                         +pow((vz[ieau]*(1.-gravite[2])),2.));
              /* fin interpolation air->eau */

              xsate = cs_ctwr_xsath(ct->teau[ieau]);
              xsata = cs_ctwr_xsath(tai[ieau]);
              dvg=sqrt(pow((vgin+vvai),2.)+pow(vhai,2.));
              if (xai[ieau]<=xsata) {
                cpx = cpa+xai[ieau]*cpv;
              }
              else  {
                cpx = cpa+xsata*cpv+(xai[ieau]-xsata)*cpe;
              }
              rre  = dvg*rhoai[ieau]*(1.+xai[ieau])*dgout/visc;
              rpr  = cpx*visc/conduc;
              anu  = 2.+0.6*sqrt(rre)*pow(rpr,(1./3.));
              bxan = 6.*conduc*anu*ct->fem[ieau]
                    /(0.92*rhoe*vgin*pow(dgout,2.)*cpx);
              cfen = bxan/ct->fem[ieau]/cpe;
              eta = (0.622+xsate)/(0.622+xai[ieau]);
              xlew = pow(0.866,(2./3.))*(eta-1.)/log(eta);
              if (xai[ieau]<=xsata) {
                ff1 = (cpa+cpv*xsate)+(cpa+cpv*xai[ieau])*(xlew-1.)-(xsate-xai[ieau])*cpe;
                ff2 = xlew*(cpa+cpv*xai[ieau])*tai[ieau] -(xsate-xai[ieau])*hv0;
              }
              else {
                ff1 = xlew*(cpa+cpv*xsata)+(xsate-xsata)*(cpv-cpe)
                      -(xsate-xsata)*cpe;
                ff2 = xlew*(cpa+cpv*xsata+cpe*(xai[ieau]-xsata))*tai[ieau]
                      -(xsate-xsata)*(cpe*tai[ieau]+hv0);
              }
              ct->teau[ieau] = (ct->teau[ieau-1]-cfen*ff2*dhi)
                             /(1.0-cfen*ff1*dhi);
            }
            else {
              ct->teau[ieau]=ct->teau[ieau-1];
            }

          }
          /* Fin de resolution Te */
        }
      }
      /*--------------------------------------------*/
      /* fin sur la zone de pluie ntype=3           */
      /*--------------------------------------------*/
    }
    /*--------------------------------------------*/
    /* Fin Modele de Poppe                        */
    /*--------------------------------------------*/

    /*--------------------------------------------*/
    /* Modele de Merkel                           */
    /*--------------------------------------------*/
    if (ct->imctch==2) {

      if (ct->ntypct<=2) {

      for (iseg = 0 ; iseg < ct->nnpsct ; iseg++) {

        for (iloc = 1 ; iloc < ct->nelect ; iloc++) {

          ieau = iseg*ct->nelect + iloc ;

          vvai = sqrt( pow((vx[ieau]*gravite[0]),2.)
                      +pow((vy[ieau]*gravite[1]),2.)
                      +pow((vz[ieau]*gravite[2]),2.));
          vhai = sqrt( pow((vx[ieau]*(1.-gravite[0])),2.)
                      +pow((vy[ieau]*(1.-gravite[1])),2.)
                      +pow((vz[ieau]*(1.-gravite[2])),2.));
          /* fin interpolation air->eau */

          ct->fem[ieau]=ct->fem[ieau-1];
          xsate = cs_ctwr_xsath(ct->teau[ieau]);
          xsata = cs_ctwr_xsath(tai[ieau]);
          if (ct->ntypct==1) {
           faxn = rhoai[ieau]*vvai;
          }
          if (ct->ntypct==2) {
            faxn = rhoai[ieau]*vhai;
          }

          bxan = ct->xap*ct->fem[ieau]*pow((faxn/ct->fem[ieau]),ct->xnp);
          cfen = bxan/ct->fem[ieau]/cpe;
          ff1 = cpa+cpv*xsate;
          ff2 = (cpa+xsata*cpv)*tai[ieau]-(xsate-xsata)*hv0;
          ct->teau[ieau] = (ct->teau[ieau-1]-cfen*ff2*dhi)
                      /(1.0-cfen*ff1*dhi);
        }
      } /* fin sur iseg */
    } /* fin if pour le ntypct<=2 */

    else if (ct->ntypct==3){ /* zone de pluie */

      for (iseg = 0 ; iseg < ct->nnpsct ; iseg++) {

        for (iloc = 1 ; iloc < ct->nelect ; iloc++) {
          ieau = iseg*ct->nelect + iloc ;
          ct->fem[ieau]=ct->fem[ieau-1];
          vgin=ct->vgoutte[ieau];

          if (CS_ABS(vgin)>=0.1){

            vvai = sqrt( pow((vx[ieau]*gravite[0]),2.)
                        +pow((vy[ieau]*gravite[1]),2.)
                        +pow((vz[ieau]*gravite[2]),2.));
            vhai = sqrt(pow((vx[ieau]*(1.-gravite[0])),2.)
                       +pow((vy[ieau]*(1.-gravite[1])),2.)
                       +pow((vz[ieau]*(1.-gravite[2])),2.));
            /* fin interpolation air->eau */

            xsate = cs_ctwr_xsath(ct->teau[ieau]);
            xsata = cs_ctwr_xsath(tai[ieau]);
            dvg=sqrt(pow((vgin+vvai),2.)+pow(vhai,2.));
            cpx = cpa+xai[ieau]*cpv;
            rre  = dvg*rhoai[ieau]*(1.+xai[ieau])*dgout/visc;
            rpr  = cpx*visc/conduc;
            anu  = 2.+0.6*sqrt(rre)*pow(rpr,(1./3.));
            bxan = 6.*conduc*anu*ct->fem[ieau]
                  /(0.92*rhoe*vgin*pow(dgout,2.)*cpx);
            cfen = bxan/ct->fem[ieau]/cpe;
            ff1 = cpa+cpv*xsate;
            ff2 = (cpa+xsata*cpv)*tai[ieau]-(xsate-xsata)*hv0;
            ct->teau[ieau]=(ct->teau[ieau-1]-cfen*ff2*dhi)
                            /(1.0-cfen*ff1*dhi);
            }
            else {
              ct->teau[ieau]=ct->teau[ieau-1];
            }
          }
        } /* fin sur iseg */
      }   /* fin sur la zone de pluie ntype=3 */
    }
    /*--------------------------------------------*/
    /* Fin du modele de Merkel                    */
    /*--------------------------------------------*/

    BFT_FREE( tai  );
    BFT_FREE( xai  );
    BFT_FREE( rhoai );
    BFT_FREE( vx  );
    BFT_FREE( vy  );
    BFT_FREE( vz  );

  }
  /*--------------------------------------------*/
  /* Fin de résolution des variable eau         */
  /* sur chaque  ct                             */
  /*--------------------------------------------*/
}

/*---------------------------------------------------------------------------*
* Function cs_ctwr_aeteau                                                      *
* Resolution des variables eau                                               *
*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------
* Function cs_ctwr_aetssc
* Calcul des termes source pour l'air
*----------------------------------------------------------------------------*/
void cs_ctwr_aetssc
(
  const cs_int_t    iscal,       /*   */

  cs_real_t   temp[],      /* Temperature air */
  cs_real_t   xa[],        /* humidite air */
  cs_real_t   rho[],       /* masse volumique air */
  cs_real_t   utsim[],           /* vitesse verticale air */
  cs_real_t   utsex[],           /* vitesse horizontale air */
  cs_real_t   vitx[],      /* vitesse air suivant x */
  cs_real_t   vity[],      /* vitesse air suivant y */
  cs_real_t   vitz[],      /* vitesse air suivant z */
  const cs_real_t   gx,          /* composante x de la gravite */
  const cs_real_t   gy,          /* composante y de la gravite */
  const cs_real_t   gz           /* composante z de la gravite */
)
{
  cs_int_t  ict,iseg,iloc,ieau,iair,i,nb,nb_dist_water,nb_dist_air;
  cs_real_t cd1,ain,bin,dhi,dvga,gravite[3],norme_g;
  cs_real_t fax,fx0,vvai,vhai,tim,tex,xim,xex;
  cs_real_t bxa,xsata,xsate,ff1,xlew,eta;
  cs_real_t dvg,cpx,rre,rpr,anu;
  cs_real_t   cpe, cpv, cpa, hv0, dgout, visc, conduc, rhoe;
  cs_ctwr_zone_t  *ct;
  cs_real_t *tai_inter, *xai_inter, *rhoai_inter,*vx_inter, *vy_inter,*vz_inter,
            *tei_inter, *femei_inter, *vgin_inter;
  cs_real_t *tai, *xai, *rhoai,*vx, *vy, *vz, *tei, *femei, *vgin;
  fvm_lnum_t  *lst_par_cel;

  gravite[0] = -gx;
  gravite[1] = -gy;
  gravite[2] = -gz;

  norme_g = sqrt( pow(gravite[0],2.)
                  +pow(gravite[1],2.)
                  +pow(gravite[2],2.) );

  gravite[0] /= norme_g;
  gravite[1] /= norme_g;
  gravite[2] /= norme_g;


  /*--------------------------------------------*/
  /* Calcul de la vitesse des gouttes pour les  */
  /* zones de pluie                             */
  /*--------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {

    ct = cs_glob_ct_tab[cs_chain_ct[ict]];
    cpa    = ct->cpa ;
    cpv    = ct->cpv ;
    cpe    = ct->cpe ;
    hv0    = ct->hv0 ;
    rhoe   = ct->rhoe ;
    dgout  = ct->dgout ;
    visc   = ct->visc ;
    conduc = ct->conduc  ;
    if (ct->ntypct==3){
      /*--------------------------------------------*
      * synchronisation   Halo                            *
      *--------------------------------------------*/
      if( ct->n_ghost_npsct >0 ){

        cs_halo_t *halo = ct->water_halo;

        cs_halo_sync_var(halo, ct->halo_type, temp);
        cs_halo_sync_var(halo, ct->halo_type, xa);
        cs_halo_sync_var(halo, ct->halo_type, rho);
        cs_halo_sync_var(halo, ct->halo_type, vitx);
        cs_halo_sync_var(halo, ct->halo_type, vity);
        cs_halo_sync_var(halo, ct->halo_type, vitz);

      }

      /*--------------------------------------------*
      * interpolation  air->eau                   *
      *--------------------------------------------*/
      nb_dist_water = (int) fvm_locator_get_n_dist_points(ct->locat_air_water);

      BFT_MALLOC( tai_inter  , nb_dist_water, cs_real_t);
      BFT_MALLOC( xai_inter  , nb_dist_water, cs_real_t);
      BFT_MALLOC( rhoai_inter, nb_dist_water, cs_real_t);
      BFT_MALLOC( vx_inter   , nb_dist_water, cs_real_t);
      BFT_MALLOC( vy_inter   , nb_dist_water, cs_real_t);
      BFT_MALLOC( vz_inter   , nb_dist_water, cs_real_t);

      for (ieau= 0 ; ieau < nb_dist_water ; ieau++) {
        tai_inter[ieau]   = 0.;
        xai_inter[ieau]   = 0.;
        rhoai_inter[ieau] = 0.;
        vx_inter[ieau]    = 0.;
        vy_inter[ieau]    = 0.;
        vz_inter[ieau]    = 0.;
        for (i = (ct->pvoiseau[ieau]) ; i < (ct->pvoiseau[ieau+1]) ; i++) {
          tai_inter[ieau]  += ct->coefeau[i] * temp[ct->voiseau[i]];
          xai_inter[ieau]  += ct->coefeau[i] * xa[ct->voiseau[i]];
          rhoai_inter[ieau]+= ct->coefeau[i] * rho[ct->voiseau[i]];
          vx_inter[ieau]   += ct->coefeau[i] * vitx[ct->voiseau[i]];
          vy_inter[ieau]   += ct->coefeau[i] * vity[ct->voiseau[i]];
          vz_inter[ieau]   += ct->coefeau[i] * vitz[ct->voiseau[i]];
        }
      }
      BFT_MALLOC( tai  , ct->nnpsct*ct->nelect, cs_real_t );
      BFT_MALLOC( xai  , ct->nnpsct*ct->nelect, cs_real_t );
      BFT_MALLOC( rhoai, ct->nnpsct*ct->nelect, cs_real_t );
      BFT_MALLOC( vx   , ct->nnpsct*ct->nelect, cs_real_t );
      BFT_MALLOC( vy   , ct->nnpsct*ct->nelect, cs_real_t );
      BFT_MALLOC( vz   , ct->nnpsct*ct->nelect, cs_real_t );

      fvm_locator_exchange_point_var(ct->locat_air_water,
                                   tai_inter, tai, NULL, sizeof(cs_real_t),1,0);
      fvm_locator_exchange_point_var(ct->locat_air_water,
                                   xai_inter, xai, NULL, sizeof(cs_real_t),1,0);
      fvm_locator_exchange_point_var(ct->locat_air_water,
                                   rhoai_inter,rhoai, NULL, sizeof(cs_real_t),1,0);
      fvm_locator_exchange_point_var(ct->locat_air_water,
                                   vx_inter,vx, NULL, sizeof(cs_real_t),1,0);
      fvm_locator_exchange_point_var(ct->locat_air_water,
                                   vy_inter,vy, NULL, sizeof(cs_real_t),1,0);
      fvm_locator_exchange_point_var(ct->locat_air_water,
                                   vz_inter,vz, NULL, sizeof(cs_real_t),1,0);
      /*--------------------------------------------*
      *  end interpolation  air->eau              *
      *--------------------------------------------*/


      dhi = -(ct->hmax-ct->hmin)/(ct->nelect-1);

      nb = (int) fvm_nodal_get_n_entities(ct->cell_mesh, 3);

      BFT_MALLOC( lst_par_cel , nb, fvm_lnum_t );
      fvm_nodal_get_parent_num( ct->cell_mesh, 3, lst_par_cel);

      for (iseg = 0 ; iseg < ct->nnpsct ; iseg++) {

        for (iloc = 1 ; iloc < ct->nelect ; iloc++) {

          ieau = iseg*ct->nelect + iloc ;

          vvai = sqrt(pow((vx[ieau]*gravite[0]),2.)
                   +pow((vy[ieau]*gravite[1]),2.)
                   +pow((vz[ieau]*gravite[2]),2.));
          /* fin interpolation air->eau */

          dvga = CS_ABS(ct->vgoutte[ieau]+vvai);

          rre  = dvga*rhoai[ieau]*dgout/visc;
          cd1 = (1.+0.15*pow(rre,0.687));
          ain = (18.*visc*cd1)/(rhoe*pow(dgout,2.));
          bin = -ain*dvga + 9.81 ;
          if (bin>0.) {
            ff1 = 2.*bin*dhi;
          }
          else {
            ff1 = 0.;
          }
          ct->vgoutte[ieau] = sqrt((pow(ct->vgoutte[ieau-1],2.)-ff1));
        }
      }
      BFT_FREE( lst_par_cel);
      BFT_FREE( tai  );
      BFT_FREE( xai  );
      BFT_FREE( rhoai );
      BFT_FREE( vx  );
      BFT_FREE( vy  );
      BFT_FREE( vz  );
      BFT_FREE( tai_inter  );
      BFT_FREE( xai_inter  );
      BFT_FREE( rhoai_inter );
      BFT_FREE( vx_inter  );
      BFT_FREE( vy_inter  );
      BFT_FREE( vz_inter  );
    }

  }  /* fin boucle ict sur les ct */
  /* Fin du calcul de la vitesse des gouttes pour les zones de pluie */

  /*--------------------------------------------*/
  /* Calcul des termes sources pour T et x      */
  /* pour chaque ct                             */
  /*--------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {
    ct = cs_glob_ct_tab[cs_chain_ct[ict]];

    if ((ct->ntypct >= 2) && ( ct->idimct==2) )
     gravite[2] = 1.0;

    cpa    = ct->cpa ;
    cpv    = ct->cpv ;
    cpe    = ct->cpe ;
    hv0    = ct->hv0 ;
    rhoe   = ct->rhoe ;
    dgout  = ct->dgout ;
    visc   = ct->visc ;
    conduc = ct->conduc  ;

    nb = (int) fvm_nodal_get_n_entities(ct->cell_mesh, 3);

    /*--------------------------------------------*
    * synchronisation Halo                        *
    *--------------------------------------------*/
    if (ct->n_ghost_npsct > 0) {
      cs_halo_t *halo = ct->water_halo;
      cs_halo_sync_var(halo, ct->halo_type, ct->teau);
      cs_halo_sync_var(halo, ct->halo_type, ct->fem);
      cs_halo_sync_var(halo, ct->halo_type, ct->vgoutte);
    }


    BFT_MALLOC( lst_par_cel , nb, fvm_lnum_t );
    fvm_nodal_get_parent_num( ct->cell_mesh, 3, lst_par_cel);
    /*--------------------------------------------*
     * interpolation  eau->air                    *
     *--------------------------------------------*/
    nb_dist_air = (int) fvm_locator_get_n_dist_points(ct->locat_water_air);

    BFT_MALLOC( tei_inter   , nb_dist_air, cs_real_t );
    BFT_MALLOC( femei_inter , nb_dist_air, cs_real_t );
    BFT_MALLOC( vgin_inter  , nb_dist_air, cs_real_t );

    for (iair= 0 ; iair < nb_dist_air ; iair++) {
       tei_inter  [ iair ] = 0.;
       femei_inter[ iair ] = 0.;
       vgin_inter [ iair ] = 0.;

      for (i = (ct->pvoisair[iair]) ; i < (ct->pvoisair[iair+1]) ; i++) {
        tei_inter[iair]    += ct->coefair[ i ]* ct->teau   [ ct->voisair[i] ];
        femei_inter[iair]  += ct->coefair[ i ]* ct->fem    [ ct->voisair[i] ];
        vgin_inter[iair]   += ct->coefair[ i ]* ct->vgoutte[ ct->voisair[i] ];
      }
    }
    BFT_MALLOC( tei   , ct->nbevct, cs_real_t );
    BFT_MALLOC( femei , ct->nbevct, cs_real_t );
    BFT_MALLOC( vgin  , ct->nbevct, cs_real_t );

    fvm_locator_exchange_point_var(ct->locat_water_air,
                                   tei_inter,     tei, NULL, sizeof(cs_real_t),1,0);
    fvm_locator_exchange_point_var(ct->locat_water_air,
                                   femei_inter, femei, NULL, sizeof(cs_real_t),1,0);
    fvm_locator_exchange_point_var(ct->locat_water_air,
                                   vgin_inter,   vgin, NULL, sizeof(cs_real_t),1,0);



    /*--------------------------------------------*
     *  end interpolation  air->eau               *
     *--------------------------------------------*/

    /*--------------------------------------------*
     * Modele de Poppe                            *
     *--------------------------------------------*/
    if (ct->imctch==1)  {
      /*--------------------------------------------*
       * courant-croise ou contre-courant           *
       *--------------------------------------------*/
      if (ct->ntypct<=2) {
        for (iloc = 0 ; iloc < ct->nbevct ; iloc++){
          iair = lst_par_cel[iloc]-1;
          /* fin interpolation eau->air */

          if (femei[iloc]>1.e-6) {
            vvai = sqrt(pow((vitx[iair]*gravite[0]),2.)
                       +pow((vity[iair]*gravite[1]),2.)
                       +pow((vitz[iair]*gravite[2]),2.));
            vhai = sqrt(pow((vitx[iair]*(1.-gravite[0])),2.)
                       +pow((vity[iair]*(1.-gravite[1])),2.)
                       +pow((vitz[iair]*(1.-gravite[2])),2.));
            if (ct->ntypct==1) {
              fax = rho[iair]*vvai;
            }
            else {
              fax = rho[iair]*vhai;
            }
            bxa = ct->xap*femei[iloc]*pow((fax/femei[iloc]),ct->xnp);
            xsata = cs_ctwr_xsath(temp[iair]);
            xsate = cs_ctwr_xsath(tei[iloc]);
            eta = (0.622+xsate)/(0.622+xa[iair]);
            xlew = pow(0.866,(2./3.))*(eta-1.)/log(eta);
            if (xa[iair]<=xsata) {
              tex = ((cpa+xsate*cpv)+(xlew-1.)*(cpa+xa[iair]*cpv))*tei[iloc];
              tim = ((cpa+xsate*cpv)+(xlew-1.)*(cpa+xa[iair]*cpv));
              xex = xsate;
              xim = 1.;
            }
            else {
              tex=xlew*(cpa+cpv*xsata+(xa[iair]-xsata)*cpe)*tei[iloc]
                  +(xsate-xsata)*cpv*tei[iloc]+(xsate-xsata)*hv0;
              tim = xlew*(cpa+cpv*xsata+(xa[iair]-xsata)*cpe)+ (xsate-xsata)*cpe;
              xex = xsate - xsata;
              xim = 0.;
            }
            /* termes sources pour T */
            if (iscal==1) {
              utsex[iair] = bxa*tex;
              utsim[iair] = bxa*tim;
            }
            /* termes sources pour x */
            if (iscal==2) {
              utsex[iair] = bxa*xex;
              utsim[iair] = bxa*xim;
            }
          }
        }
      }
      /* Fin courant-croise ou contre-courant  */

      /*--------------------------------------------*/
      /* zone de pluie                              */
      /*--------------------------------------------*/
      else if (ct->ntypct==3) {

        for (iloc = 0 ; iloc < ct->nbevct ; iloc++) {
          iair = lst_par_cel[iloc]-1;

          if (CS_ABS(vgin[iloc])>=0.1) {
            vvai = sqrt(pow(( vitx[iair]*gravite[0]),2.)
                       +pow(( vity[iair]*gravite[1]),2.)
                       +pow(( vitz[iair]*gravite[2]),2.));
            vhai = sqrt(pow(( vitx[iair]*(1.-gravite[0])),2.)
                       +pow(( vity[iair]*(1.-gravite[1])),2.)
                       +pow(( vitz[iair]*(1.-gravite[2])),2.));
            dvg = sqrt(pow((vvai+vgin[iloc]),2.)+pow(vhai,2.));

            bxa = ct->xap*femei[iloc]*pow((fax/femei[iloc]),ct->xnp);
            xsata = cs_ctwr_xsath(temp[iair]);
            xsate = cs_ctwr_xsath(tei[iloc]);
            if (xa[iair]<=xsata) {
              cpx = cpa+xa[iair]*cpv;
            }
            else {
            cpx = cpa+xsata*cpv+(xa[iair]-xsata)*cpe;
            }
            rre = dvg*rho[iair]*(1.+xsata)*dgout/visc;
            rpr = cpx*visc/conduc;
            anu = 2.+0.6*sqrt(rre)*pow(rpr,(1./3.));
            bxa = (6.*conduc*anu*femei[iloc])/(0.92*rhoe*vgin[iloc]*pow(dgout,2.)*cpx);
            eta = (0.622 + xsate)/(0.622 + xa[iair]);
            xlew = pow(0.866,(2./3.))*(eta-1.)/log(eta);
            if (xa[iair]<=xsata) {
              tex = ((cpa+xsate*cpv)+(xlew-1.)*(cpa+xa[iair]*cpv))*tei[iloc];
              tim = ((cpa+xsate*cpv)+(xlew-1.)*(cpa+xa[iair]*cpv));
              xex = xsate;
              xim = 1.;
            }
            else {
              tex = xlew*(cpa+cpv*xsata+(xa[iair]-xsata)*cpe)*tei[iloc]
                    +(xsate-xsata)*cpv*tei[iloc]+(xsate-xsata)*hv0;
              tim = xlew*(cpa+cpv*xsata+(xa[iair]-xsata)*cpe)+ (xsate-xsata)*cpe;
              xex = xsate - xsata;
              xim = 0.;
            }
            /* termes sources pour T */
            if (iscal==1){
              utsex[iair] = bxa*tex;
              utsim[iair] = bxa*tim;
            }
            /* termes sources pour x */
            if (iscal==2){
              utsex[iair] = bxa*xex;
              utsim[iair] = bxa*xim;
            }
          }
        }
      }
      /* Fin de la zone de pluie  */
    }
    /*--------------------------------------------*/
    /* Fin du modele de Poppe                     */
    /*--------------------------------------------*/

    /*--------------------------------------------*/
    /* Modele de Merkel                           */
    /*--------------------------------------------*/
    if (ct->imctch==2)  {
      if (ct->ntypct<=2){
        for (iloc = 0 ; iloc < ct->nbevct ; iloc++){
          iair = lst_par_cel[iloc]-1;
          if (femei[iloc]>1.e-6) {
            vvai = sqrt(pow(( vitx[iair]*gravite[0]),2.)
                       +pow(( vity[iair]*gravite[1]),2.)
                       +pow(( vitz[iair]*gravite[2]),2.));
            vhai = sqrt(pow(( vitx[iair]*(1.-gravite[0])),2.)
                       +pow(( vity[iair]*(1.-gravite[1])),2.)
                       +pow(( vitz[iair]*(1.-gravite[2])),2.));

            if (ct->ntypct==1) {
              fax=rho[iair]*vvai;
            }
            else {
              fax=rho[iair]*vhai;
            }

            xsata = cs_ctwr_xsath(temp[iair]);
            xsate = cs_ctwr_xsath(tei[iloc]);
            bxa = ct->xap*femei[iloc]*pow((fax/femei[iloc]),ct->xnp);
            fx0 = (xsate-xsata)*(cpv*tei[iloc]+hv0);
            if (iscal==1) {
              utsex[iair] = bxa*tei[iloc]*(cpa+cpv*xsata)+bxa*fx0;
              utsim[iair] = bxa*(cpa+cpv*xsata);
            }
            if (iscal==2) {
              utsex[iair] = 1.e20*xsata;
              utsim[iair] = 1.e20;
            }
          }
        }
      } /*Fin du ntypct<=2 */

      else if (ct->ntypct==3){ /* zone de pluie */

        for (iloc = 0 ; iloc < ct->nbevct ; iloc++){
          iair = lst_par_cel[iloc]-1;

          if (CS_ABS(vgin[iloc])>=0.1) {

            vvai = sqrt(pow(( vitx[iair]*gravite[0]),2.)
                       +pow(( vity[iair]*gravite[1]),2.)
                       +pow(( vitz[iair]*gravite[2]),2.));
            vhai = sqrt(pow(( vitx[iair]*(1.-gravite[0])),2.)
                       +pow(( vity[iair]*(1.-gravite[1])),2.)
                       +pow(( vitz[iair]*(1.-gravite[2])),2.));

            dvg = sqrt(pow((vvai+vgin[iloc]),2.)+pow(vhai,2.));
            xsata = cs_ctwr_xsath(temp[iair]);
            xsate = cs_ctwr_xsath(tei[iloc]);
            cpx = cpa+xsata*cpv;
            rre = dvg*rho[iair]*(1.+xsata)*dgout/visc;
            rpr = cpx*visc/conduc;
            anu = 2.+0.6*sqrt(rre)*pow(rpr,(1./3.));
            bxa = (6.*conduc*anu*femei[iloc])/(0.92*rhoe*vgin[iloc]*pow(dgout,2.)*cpx);
            fx0 = (xsate-xsata)*(cpv*tei[iloc]+hv0);
            if (iscal==1){
              utsex[iair] = bxa*tei[iloc]*(cpa+cpv*xsata)+bxa*fx0;
              utsim[iair] = bxa*(cpa+cpv*xsata);
            }
            if (iscal==2) {
              utsex[iair] = 1.e20*xsata;
              utsim[iair] = 1.e20;
            }
          }
        }
      }/* Fin de la zone de pluie ntypct=3 */
    }
    /*--------------------------------------------*/
    /* Fin pour le modele de Merkel */
    /*--------------------------------------------*/
    BFT_FREE( lst_par_cel);
    BFT_FREE( tei_inter );
    BFT_FREE( femei_inter );
    BFT_FREE( vgin_inter );
    BFT_FREE( tei );
    BFT_FREE( femei );
    BFT_FREE( vgin );
  }
  /*--------------------------------------------*/
  /* Fin  calcul des termes sources pour T et x*/
  /* pour chaque ct                             */
  /*--------------------------------------------*/
}



/*----------------------------------------------------------------------------
* Function cs_ctwr_aetsvi
* Calcul des PdC induites dans les zones de pluie
*----------------------------------------------------------------------------*/

void cs_ctwr_aetsvi
(
  const cs_int_t    idim,
  const cs_real_t   rho[],       /* masse volumique air */
  const cs_real_t   vitx[],      /* vitesse air suivant x */
  const cs_real_t   vity[],      /* vitesse air suivant y */
  const cs_real_t   vitz[],      /* vitesse air suivant z */
  const cs_real_t   xair[],      /* humidite de l'air */
  const cs_real_t   gx,          /*   */
  const cs_real_t   gy,          /*   */
  const cs_real_t   gz,          /*   */
  cs_real_t   utsex[]            /* terme source explicite */
)
{
  cs_int_t  ict, iloc, iair, i, *lst_par_cel, nb,nb_dist_air;
  cs_real_t dgout, visc, rhoe;
  cs_real_t absgrv, vginu,vginv,vginw,dvg,qer,rre,cdd1,cff0;

  cs_real_t  *femei_inter, *vgin_inter;
  cs_real_t  *femei, *vgin;
  cs_ctwr_zone_t  *ct;

  absgrv = sqrt(pow(gx,2.)+pow(gy,2.)+pow(gz,2.));

  /*--------------------------------------------*/
  /* Calcul de Kg pour chaque ct                */
  /*--------------------------------------------*/
  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {
    ct = cs_glob_ct_tab[cs_chain_ct[ict]];
    rhoe   = ct->rhoe ;
    dgout  = ct->dgout ;
    visc   = ct->visc ;

     /*--------------------------------------------*
    * synchronisation Halo                        *
    *--------------------------------------------*/
    if (ct->n_ghost_npsct > 0) {
      cs_halo_t *halo = ct->water_halo;
      cs_halo_sync_var(halo, ct->halo_type, ct->teau);
      cs_halo_sync_var(halo, ct->halo_type, ct->fem);
      cs_halo_sync_var(halo, ct->halo_type, ct->vgoutte);
    }

    nb = (int) fvm_nodal_get_n_entities(ct->cell_mesh, 3);

    BFT_MALLOC( lst_par_cel , (nb*3), fvm_lnum_t );
    fvm_nodal_get_parent_num( ct->cell_mesh, 3, lst_par_cel);
    /*--------------------------------------------*
     * interpolation  eau->air                    *
     *--------------------------------------------*/
    nb_dist_air = (int) fvm_locator_get_n_dist_points(ct->locat_water_air);

    BFT_MALLOC( femei_inter  , nb_dist_air, cs_real_t);
    BFT_MALLOC( vgin_inter  , nb_dist_air, cs_real_t);

    for (iair= 0 ; iair < nb_dist_air ; iair++) {

       femei_inter[iair] = 0.;
       vgin_inter[iair]  = 0.;

      for (i = (ct->pvoisair[iair]) ; i < (ct->pvoisair[iair+1]) ; i++) {

        femei_inter[ iair ]  += ct->coefair[ i ]* ct->fem    [ ct->voisair[i] ];
        vgin_inter [ iair ]  += ct->coefair[ i ]* ct->vgoutte[ ct->voisair[i] ];
      }
    }

    BFT_MALLOC( femei , ct->nbevct, cs_real_t );
    BFT_MALLOC( vgin , ct->nbevct, cs_real_t );

    fvm_locator_exchange_point_var(ct->locat_water_air,
                                   femei_inter, femei, NULL, sizeof(cs_real_t),1,0);
    fvm_locator_exchange_point_var(ct->locat_water_air,
                                   vgin_inter, vgin, NULL, sizeof(cs_real_t),1,0);

    /*--------------------------------------------*/
    /* zone de pluie                              */
    /*--------------------------------------------*/
    if (ct->ntypct==3)  {
      for (iloc = 0 ; iloc < ct->nbevct ; iloc++) {
        iair = lst_par_cel[iloc]-1;

        vginu  = - gx / absgrv * vgin[iloc];
        vginv  = - gy / absgrv * vgin[iloc];
        vginw  = - gz / absgrv * vgin[iloc];
        dvg = sqrt( pow((vitx[iair]+vginu ),2.)
              + pow((vity[iair]+vginv ),2.)
              + pow((vitz[iair]+vginw ),2.) );
        if (vgin[iloc] > 0.1) {
          qer = femei[iloc]/rhoe;
          rre = dvg*rho[iair]*(1 + xair[iair])*dgout/visc;
          cdd1 = (1.+0.15*pow(rre,0.687));
          cff0 = 18.*cdd1*visc*qer/(vgin[iloc]*pow(dgout,2.));
          if (idim==1){ utsex[iair] = -cff0 *( vitx[iair]+vginu ) ; }
          if (idim==2){ utsex[iair] = -cff0 *( vity[iair]+vginv ) ; }
          if (idim==3){ utsex[iair] = -cff0 *( vitz[iair]+vginw ) ; }
        }
      }
    }
    /* Fin  zones de pluie */
    BFT_FREE(lst_par_cel);
    BFT_FREE(femei_inter);
    BFT_FREE(vgin_inter);
    BFT_FREE(femei);
    BFT_FREE(vgin);
  }
  /* Fin de calcul de Kg pour chaque ct */
}

/*----------------------------------------------------------------------------
* Bilan dans les ct
*----------------------------------------------------------------------------*/

void cs_ctwr_bilanct
(
  const cs_real_t   time,                /*   */
  cs_real_t   fem_entree[],       /* debit eau entree */
  cs_real_t   fem_sortie[],       /* debit eau sortie */
  cs_real_t   teau_entree[],      /* temperature eau entree */
  cs_real_t   teau_sortie[],      /* temperature eau sortie */
  cs_real_t   heau_entree[],      /* enthalpie eau entree */
  cs_real_t   heau_sortie[],      /* enthalpie eau sortie */
  cs_real_t   tair_entree[],      /* temperature air entree */
  cs_real_t   tair_sortie[],      /* temperature air sortie */
  cs_real_t   xair_entree[],      /*   */
  cs_real_t   xair_sortie[],      /*   */
  cs_real_t   hair_entree[],      /*   */
  cs_real_t   hair_sortie[],      /*   */
  cs_real_t   debit_entree[],     /*   */
  cs_real_t   debit_sortie[],     /*   */

  const cs_real_t   temp[],             /* Temperature air */
  const cs_real_t   xa[],               /* humidite air */
  const cs_real_t   flux_masse_fac[],   /* vitesse verticale air */
  const cs_real_t   flux_masse_fbr[],   /* vitesse horizontale air */
  const cs_real_t   vitx[],             /* vitesse air suivant x */
  const cs_real_t   vity[],             /* vitesse air suivant y */
  const cs_real_t   vitz[],             /* vitesse air suivant z */

  const cs_mesh_t             *mesh,      /* <-- structure maillage associée  */
  const cs_mesh_quantities_t  *mesh_quantities   /* <-- grandeurs du maillage        */
)
{
  const cs_real_t  *i_face_normal = mesh_quantities->i_face_normal;
  const cs_real_t  *b_face_normal = mesh_quantities->b_face_normal;
  const cs_int_t   *family_item = mesh->family_item; /* Propriétés des familles */
  const cs_int_t   *cell_family  = mesh->cell_family;    /* Familles des cellules */
  cs_int_t         icel_1, icel_2, icel, ifac, ict, idim, i, j, ieau_Sup,
    ieau_inf, length;
  cs_real_t   cpe, cpv, cpa, hv0;
  const cs_int_t   *i_face_cells = mesh->i_face_cells;
  const cs_int_t   *b_face_cells = mesh->b_face_cells;
  const cs_real_t  *coo_cen  = mesh_quantities->cell_cen;
  cs_real_t        xsata,debit,hair,n_sortant[3],vitair[3],aux,
                   surf,surf_e,surf_s;

  fvm_lnum_t   *face_sup;      /* liste des faces  superieures de la ct */
  fvm_lnum_t   *face_inf;      /* liste des faces  inferior de la ct */
  fvm_lnum_t   *face_lat;      /* liste des faces  inferior de la ct */

  cs_ctwr_zone_t  *ct;
  bft_file_t *f;
  char  *file_name = NULL;
  bft_file_type_t file_type;

  file_type = BFT_FILE_TYPE_TEXT;

  for (ict=0 ; ict < cs_glob_ct_nbr ; ict++) {

    ct = cs_glob_ct_tab[cs_chain_ct[ict]];
    cpa    = ct->cpa ;
    cpv    = ct->cpv ;
    cpe    = ct->cpe ;
    hv0    = ct->hv0 ;
    cs_int_t nbr_fbr_air[3][2] = {{ct->nnpsct,ct->nbfbr_sct},
                                 {(ct->nbfbr_ict + ct->nbfac_ict),ct->nbfbr_ict},
                                 {(ct->nbfbr_lct + ct->nbfac_lct),ct->nbfbr_lct}};


    BFT_MALLOC( face_sup ,(ct->nbfac_sct + ct->nbfbr_sct) ,fvm_lnum_t );
    fvm_nodal_get_parent_num( ct->face_sup_mesh, 2, face_sup);
    BFT_MALLOC( face_inf ,(ct->nbfac_ict + ct->nbfbr_ict) ,fvm_lnum_t );
    fvm_nodal_get_parent_num( ct->face_inf_mesh, 2, face_inf);
    BFT_MALLOC( face_lat ,(ct->nbfbr_lct + ct->nbfac_lct) ,fvm_lnum_t );
    fvm_nodal_get_parent_num( ct->face_lat_mesh, 2, face_lat);

    ct->fem_e   = 0.0 ;
    ct->fem_s   = 0.0 ;
    ct->teau_e  = 0.0 ;
    ct->heau_s  = 0.0 ;
    ct->heau_e  = 0.0 ;
    ct->teau_s  = 0.0 ;
    ct->tair_e  = 0.0 ;
    ct->tair_s  = 0.0 ;
    ct->xair_e  = 0.0 ;
    ct->xair_s  = 0.0 ;
    ct->hair_e  = 0.0 ;
    ct->hair_s  = 0.0 ;
    ct->debit_e = 0.0 ;
    ct->debit_s = 0.0 ;

    /* calcul des valeurs eau */

    for (i = 0 ; i < ct->nnpsct ; i++) {
       ieau_Sup = i*ct->nelect;
       ieau_inf = (i+1)*ct->nelect - 1 ;

       surf = ct->surf_fac_sup[i];

       ct->teau_e += ct->teau[ieau_Sup]*ct->fem[ieau_Sup]*surf;
       ct->fem_e  += ct->fem[ieau_Sup]*surf ;
       ct->heau_e += ct->teau[ieau_Sup]*ct->fem[ieau_Sup]*surf;

       ct->teau_s += ct->teau[ieau_inf]*ct->fem[ieau_inf]*surf;
       ct->fem_s  += ct->fem[ieau_inf]*surf ;
       ct->heau_s += ct->teau[ieau_inf]*ct->fem[ieau_inf]*surf;

    }

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {

      MPI_Allreduce (&ct->teau_e, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->teau_e = aux;

      MPI_Allreduce (&ct->fem_e, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->fem_e = aux;

      MPI_Allreduce (&ct->heau_e, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->heau_e = aux;

      MPI_Allreduce (&ct->teau_s, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->teau_s = aux;

      MPI_Allreduce (&ct->fem_s, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->fem_s = aux;

      MPI_Allreduce (&ct->heau_s, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->heau_s = aux;

    }
#endif

    ct->teau_e /= ct->fem_e ;
    ct->fem_e  /= ct->surface_in;
    ct->heau_e *= ct->cpe ;

    ct->teau_s /= ct->fem_s ;
    ct->fem_s  /= ct->surface_out  ;
    ct->heau_s *= ct->cpe ;


    /* calcul des valeurs air */

    surf_e = 0. ;
    surf_s = 0. ;

    for (j = 0 ; j < 3 ; j++)
    for (i = 0 ; i < nbr_fbr_air[j][0] ; i++) {
      if( i< nbr_fbr_air[j][1] ){
        if( j==0) ifac = (cs_int_t) face_sup[i]-1;
        if( j==1) ifac = (cs_int_t) face_inf[i]-1;
        if( j==2) ifac = (cs_int_t) face_lat[i]-1;
        icel = b_face_cells[ifac] - 1 ;
        for (idim = 0 ; idim<3 ; idim++ )
          n_sortant[idim] =  mesh_quantities->b_face_normal[ifac*3+idim];
        debit = CS_ABS(flux_masse_fbr[ifac]);
        surf  = CS_LOC_MODULE((b_face_normal + 3*ifac));
      }else{
        if( j==0) ifac = (cs_int_t) face_sup[i] - mesh->n_b_faces - 1;
        if( j==1) ifac = (cs_int_t) face_inf[i] - mesh->n_b_faces - 1;
        if( j==2) ifac = (cs_int_t) face_lat[i] - mesh->n_b_faces - 1;
        icel_1 = i_face_cells[ifac * 2]     - 1;
        icel_2 = i_face_cells[ifac * 2 + 1] - 1;
        if ( family_item[cell_family[icel_1]-1] == ct->icoul ) {

          icel = icel_2 ;
          for (idim = 0 ; idim < 3 ; idim++) {
            n_sortant[idim] =  coo_cen[icel_2*3 + idim] - coo_cen[icel_1*3 + idim];
          }
        }
        if ( family_item[cell_family[icel_2]-1] == ct->icoul ) {

          icel = icel_1 ;
          for (idim = 0 ; idim < 3 ; idim++) {
            n_sortant[idim] =  coo_cen[icel_1*3 + idim] - coo_cen[icel_2*3 + idim];
          }
        }
        debit = CS_ABS(flux_masse_fac[ifac]);
        surf  = CS_LOC_MODULE((i_face_normal + 3*ifac));
      }
      xsata = cs_ctwr_xsath(temp[icel]);
      hair = (cpa+xa[icel]*cpv)*temp[icel]+xa[icel]*hv0;
      vitair[0] = vitx[icel] ;
      vitair[1] = vity[icel] ;
      vitair[2] = vitz[icel] ;
      if (CS_LOC_PRODUIT_SCALAIRE(n_sortant, vitair)>0.) {
        surf_s += surf;
        ct->hair_s  += hair*debit;
        ct->xair_s  += debit*xa[icel];
        ct->tair_s  += debit*temp[icel];
        ct->debit_s += debit;
      }
      else {
        surf_e += surf;
        ct->hair_e  += hair*debit;
        ct->xair_e  += debit*xa[icel];
        ct->tair_e  += debit*temp[icel];
        ct->debit_e += debit;
      }
    }

#if defined(HAVE_MPI)
    if (cs_glob_n_ranks > 1) {


      MPI_Allreduce (&ct->tair_e, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->tair_e = aux;

      MPI_Allreduce (&ct->xair_e, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->xair_e = aux;

      MPI_Allreduce (&ct->debit_e, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->debit_e = aux;

      MPI_Allreduce (&ct->hair_e, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->hair_e = aux;

      MPI_Allreduce (&ct->tair_s, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->tair_s = aux;

      MPI_Allreduce (&ct->xair_s, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->xair_s = aux;

      MPI_Allreduce (&ct->debit_s, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->debit_s = aux;

      MPI_Allreduce (&ct->hair_s, &aux, 1, CS_MPI_REAL, MPI_SUM,
                      cs_glob_mpi_comm);
      ct->hair_s = aux;

    }
#endif

    if (CS_ABS( ct->debit_e )> 1e-10 ){
      ct->tair_e /= ct->debit_e ;
      ct->xair_e /= ct->debit_e ;
    }

    if (CS_ABS( ct->debit_s )> 1e-10 ){
      ct->tair_s /= ct->debit_s ;
      ct->xair_s /= ct->debit_s ;
    }


    fem_entree[ict]   = ct->fem_e ;
    fem_sortie[ict]   = ct->fem_s ;
    teau_entree[ict]  = ct->teau_e ;
    teau_sortie[ict]  = ct->teau_s ;
    heau_entree[ict]  = ct->heau_e ;
    heau_sortie[ict]  = ct->heau_s ;
    tair_entree[ict]  = ct->tair_e ;
    tair_sortie[ict]  = ct->tair_s ;
    xair_entree[ict]  = ct->xair_e ;
    xair_sortie[ict]  = ct->xair_s ;
    hair_entree[ict]  = ct->hair_e ;
    hair_sortie[ict]  = ct->hair_s ;

    ct->debit_e *= (ct->surface/ct->surface_in) ;
    ct->debit_s *= (ct->surface/ct->surface_out);

    debit_entree[ict] = ct->debit_e ;
    debit_sortie[ict] = ct->debit_s ;



    if (cs_glob_rank_id <= 0) {
      length = strlen("bltctc.") + 3 ;
      BFT_MALLOC(file_name, length, char);
      sprintf(file_name, "bltctc.%02d", ct->num);

      if (CS_ABS(ct->heau_e-ct->heau_s)> 1.e-6){
        f = bft_file_open(file_name,
                        BFT_FILE_MODE_APPEND,
                        file_type);

        aux =CS_ABS( (ct->hair_s - ct->hair_e)/(ct->heau_e - ct->heau_s) );
        bft_file_printf(f,
            "%10f\t%10f\t%10f\t%10f\t%12.5e\t%10f\t%10f\n"
                                          ,time
                                          ,aux
                                          ,ct->tair_s
                                          ,ct->teau_s
                                          ,ct->xair_s
                                          ,ct->debit_e
                                          ,ct->debit_s);

        f = bft_file_free(f);
      }
    }

    BFT_FREE(file_name);
    BFT_FREE(face_sup);
    BFT_FREE(face_inf);
    BFT_FREE(face_lat);
  } /* fin de la boucle sur les zones d'echanges */

}

/*----------------------------------------------------------------------------
 * Initialict post-processing
 *
 * parameters:
 *   ct_id         -->  Id of exchange area
 *   writer_id           -->  Id of associated writer
 *----------------------------------------------------------------------------*/

void
cs_ctwr_post_init(cs_int_t  ct_id,
                  cs_int_t  writer_id)
{
  cs_int_t  mesh_id = cs_post_get_free_mesh_id();

  cs_ctwr_zone_t * ct = cs_ctwr_by_id(ct_id);

  assert(ct != NULL);

  /* Exit silently if associated writer is not available */

  if (cs_post_writer_exists(writer_id) != CS_TRUE)
    return;

  /* Initialict post processing flag, and free previous arrays in
     case this function is called more than once */

  ct->post_mesh_id = mesh_id;

  /* Associate external mesh description with post processing subsystem */

  cs_post_add_existing_mesh(mesh_id,
                            ct->water_mesh,
                            CS_FALSE);

  cs_post_associate(mesh_id, writer_id);

  /* Register post processing function */

  cs_post_add_time_dep_var(cs_ctwr_post_function, ct_id);

  /* Update start and end (negative) numbers associated with
     dedicated post processing meshes */

  if (cs_glob_ct_post_mesh_ext[0] == 0)
    cs_glob_ct_post_mesh_ext[0] = mesh_id;

  cs_glob_ct_post_mesh_ext[1] = mesh_id;
}

/*----------------------------------------------------------------------------
 * Post process variables associated with exchange area
 *
 * parameters:
 *   coupling_id         -->  Id of exchange area
 *   nt_cur_abs          -->  Current time step
 *   t_cur_abs           -->  Current time value
 *----------------------------------------------------------------------------*/

void
cs_ctwr_post_function(cs_int_t   ct_id,
                      cs_int_t   nt_cur_abs,
                      cs_real_t  t_cur_abs)
{
  cs_ctwr_zone_t * ct = cs_ctwr_by_id(ct_id);

  if (ct->post_mesh_id != 0) {

    cs_post_write_var(ct->post_mesh_id,
                      _("T water"),
                      1,
                      CS_FALSE,
                      CS_FALSE,
                      CS_POST_TYPE_cs_real_t,
                      nt_cur_abs,
                      t_cur_abs,
                      ct->teau,
                      NULL,
                      NULL);

    cs_post_write_var(ct->post_mesh_id,
                      _("Flux water"),
                      1,
                      CS_FALSE,
                       CS_FALSE,
                      CS_POST_TYPE_cs_real_t,
                      nt_cur_abs,
                      t_cur_abs,
                      ct->fem,
                       NULL,
                      NULL);

  }

}

/*----------------------------------------------------------------------------
 * Get the local (negative) numbers associated with the first and last
 * post processing meshes dedicated to exchange area
 *
 * parameters:
 *   first_mesh_id       <--  Id of first post processing mesh
 *   last_mesh_id        <--  Id of last post processing mesh
 *----------------------------------------------------------------------------*/

void
cs_ctwr_post_id_extents(cs_int_t  *const id_mesh_start,
                        cs_int_t  *const id_mesh_end)
{
  *id_mesh_start = cs_glob_ct_post_mesh_ext[0];
  *id_mesh_end   = cs_glob_ct_post_mesh_ext[1];
}


/*----------------------------------------------------------------------------
 * Get pointer to exchange area.
 *
 * parameters:
 *   ct_id  -->  Id (0 to n-1) of exchange area
 *
 * returns:
 *   pointer to exchange area structure
 *----------------------------------------------------------------------------*/

cs_ctwr_zone_t *
cs_ctwr_by_id(cs_int_t ct_id)
{
  cs_ctwr_zone_t  *retval = NULL;

  if (   ct_id > -1
      && ct_id <  cs_glob_ct_nbr)
    retval = cs_glob_ct_tab[ct_id];

  return retval;
}


#undef CS_LOC_PRODUIT_SCALAIRE
#undef CS_LOC_MODULE
#undef CS_PRODUIT_VECTORIEL

/*----------------------------------------------------------------------------*/

END_C_DECLS
