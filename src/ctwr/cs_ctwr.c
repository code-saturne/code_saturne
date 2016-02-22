/*============================================================================
 * Definitions, Global variables variables, and functions associated with the
 * exchange zones
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_printf.h"
#include "bft_mem.h"

#include "fvm_nodal_extract.h"

#include "cs_base.h"
#include "cs_ctwr_air_props.h"
#include "cs_ctwr_halo.h"
#include "cs_halo.h"
#include "cs_math.h"
#include "cs_mesh_location.h"
#include "cs_post.h"
#include "cs_restart.h"
#include "cs_selector.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_ctwr.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/* array of exchanges area */

cs_lnum_t            cs_glob_ct_nbr_max = 0;

cs_lnum_t            cs_glob_ct_nbr     = 0;
cs_ctwr_zone_t     ** cs_glob_ct_tab   = NULL;

/* Start and end (negative) numbers associated with
   dedicated post processing meshes */

static int  cs_glob_ct_post_mesh_ext[2] = {0, 1};

/* array containing the stacking of the exchange area*/
cs_lnum_t  *  cs_stack_ct    = NULL;

/* array containing the treatment order of the exchanges areas */
cs_lnum_t  *  cs_chain_ct = NULL;

/* Restart file */

static cs_restart_t *cs_glob_ctwr_suite = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Fonctions publiques pour API Fortran
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Ajout d'une zone d'echange
 *
 * Interface Fortran :
 *
 * SUBROUTINE DEFCT1
 *
 * INTEGER          ICOUL       : <-- : Couleur des elements de la zone d'echange
 * INTEGER          IMctCH      : <-- :
 * INTEGER          NTYPct      : <-- :
 * INTEGER          NELEct      : <-- :
 * REAL             XAP         : <-- :
 * REAL             XNP         : <-- :
 * REAL             SURFACE     : <-- : Surface de la face superieure de la ct
 *----------------------------------------------------------------------------*/

void CS_PROCF (defct1, DEFCT1)
(
  const cs_int_t   *const idimct,   /* Dimemsion du probleme 2:2D  3:3D       */
  const char       *ze_name,        /* Name of Ct area */
  cs_int_t         *ze_n_len,       /* lenght of Name of Ct area */
  const cs_int_t   *const imctch,   /* 1: Modele de Poppe
                                       2: Merkel 0: Rien                      */
  const cs_int_t   *const ntypct,   /* 1: Contre courant  2: Courant croises
                                       3: Zone de pluie                      */
  const cs_int_t   *const nelect,   /* Nombre d'elements sur chaque ligne du
                                       maillage eau pour la zone de noeuds par
                                       segment eau */
  const cs_real_t  *const deltat,   /* Ecart de temperature impose en entree
                                       de la zone d'echange */
  const cs_real_t  *const teau,     /* Teau en entree de la zone d'echange    */
  const cs_real_t  *const fem,      /* fem en entree de la zone d'echange     */
  const cs_real_t  *const xap,      /* coefficient lambda de la loi d'echange */
  const cs_real_t  *const xnp,      /* exposant n de la loi d'echange         */
  const cs_real_t  *const surface,  /* Surface totale arrivee d eau de la ct  */
  const cs_real_t  *const   dgout   /* Diametre de goutte pour
                                       les zones de pluie                     */
)
{
  char *_ze_name = NULL;

  if (ze_name != NULL && *ze_n_len > 0)
    _ze_name = cs_base_string_f_to_c_create(ze_name,
                                                      *ze_n_len);
  if (_ze_name != NULL && strlen(_ze_name) == 0)
    cs_base_string_f_to_c_free(&_ze_name);

  cs_ctwr_definit(*idimct,_ze_name, *imctch,*ntypct,*nelect,
                  *deltat,*teau,*fem,*xap,*xnp,*surface,*dgout);

  if (_ze_name != NULL)
    cs_base_string_f_to_c_free(&_ze_name);
}

/*----------------------------------------------------------------------------
 * Recuperation du nombre de zones d'echanges
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
  *nbrct = cs_glob_ct_nbr;
}

/*----------------------------------------------------------------------------
 * Recuperation du modele de Poppe ou de Merkel
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
 * Resolution des variables eau
 *
 * Interface Fortran :
 *
 * SUBROUTINE AETEAU ( )
 *----------------------------------------------------------------------------*/

void CS_PROCF (aeteau, AETEAU)
(
  cs_real_t          temp[],              /* Temperature air */
  cs_real_t          xa[]  ,              /* humidite air */
  cs_real_t          rho[] ,              /* masse volumique air */
  cs_real_t          vitx[],              /* vitesse air suivant x */
  cs_real_t          vity[],              /* vitesse air suivant y */
  cs_real_t          vitz[]               /* vitesse air suivant z */

)
{
  cs_ctwr_aeteau(temp,xa,rho,vitx,vity,vitz);
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

  cs_real_t             temp[] ,   /* Temperature air            */
  cs_real_t             xa[]   ,   /* humidite air               */
  cs_real_t             rho[]  ,   /* masse volumique air        */
  cs_real_t             utsim[],   /* vitesse verticale air      */
  cs_real_t             utsex[],   /* vitesse horizontale air    */
  cs_real_t             vitx[] ,   /* vitesse air suivant x      */
  cs_real_t             vity[] ,   /* vitesse air suivant y      */
  cs_real_t             vitz[]     /* vitesse air suivant z      */
)
{
  cs_ctwr_aetssc(*iscal, temp,xa,rho,utsim,utsex,vitx,vity,vitz);
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
  cs_real_t   utsex[]                   /* terme source explicite */

)
{
  cs_ctwr_aetsvi(*idim,rho,vitx,vity,vitz,xair,utsex);
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
 * Write the restart file of the cooling tower module
 *
 * Fortran interface:
 *
 * subroutine ecrctw
 * *****************
 *
 * character(kind=c_char)  nomsui : <-- : Name of the restart file
 *----------------------------------------------------------------------------*/

void CS_PROCF (ecrctw, ECRCTW)
(
 const char  *nomsui
)
{
  int  nbvent;
  int  ict;

  cs_restart_t            *suite;
  cs_mesh_location_type_t  location_id,support;
  cs_restart_val_type_t    typ_val;

  cs_ctwr_zone_t  *ct;
  char            *location_name = NULL;
  cs_int_t         length        = 0;
  cs_gnum_t       *g_elt_num     = NULL;

  cs_lnum_t n_g_elements, n_elements;

  /* Open the restart file */

  cs_glob_ctwr_suite
    = cs_restart_create(nomsui, NULL, CS_RESTART_MODE_WRITE);

  /* Pointer to the global restart structure */
  suite = cs_glob_ctwr_suite;

  if (cs_glob_ctwr_suite == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the cooling tower module restart "
                "file in write mode.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              nomsui);

  for (ict=0; ict < cs_glob_ct_nbr; ict++) {

    ct = cs_glob_ct_tab[ict];

    length = strlen("Cooling_Tower_restart_") + 3;
    BFT_MALLOC(location_name, length, char);
    sprintf(location_name, "Cooling_Tower_restart_%02d", ct->num);

    n_g_elements = fvm_nodal_get_n_g_elements(ct->water_mesh, FVM_CELL_HEXA);
    n_elements = fvm_nodal_get_n_elements(ct->water_mesh, FVM_CELL_HEXA);

    BFT_MALLOC(g_elt_num, n_g_elements, cs_gnum_t);

    fvm_nodal_get_global_element_num(ct->water_mesh,
                                     FVM_CELL_HEXA,
                                     g_elt_num);

    location_id = cs_restart_add_location(suite,
                                          location_name,
                                          n_g_elements,
                                          n_elements,
                                          g_elt_num);


    { /* Write the header */
      char * nomrub;
      cs_int_t   *tabvar;

      length = strlen("Parametres_int_ctwr_") + 3;
      BFT_MALLOC(nomrub, length, char);
      sprintf(nomrub, "Parametres_int_ctwr_%02d", ct->num);


      BFT_MALLOC(tabvar, 3, cs_int_t);

      tabvar[ 0 ] = ct->imctch; /* modele*/
      tabvar[ 1 ] = ct->ntypct; /* Type*/
      tabvar[ 2 ] = ct->nelect; /* nb of node per segment*/

      nbvent  = 3;
      support = CS_MESH_LOCATION_NONE;
      typ_val = CS_TYPE_cs_int_t;

      cs_restart_write_section(suite,
                               nomrub,
                               support,
                               nbvent,
                               typ_val,
                               tabvar);

      BFT_FREE(tabvar);
    }

    {/* Write the header */
      char * nomrub;
      cs_real_t   *tabvar;

      length = strlen("Parametres_real_ctwr_") + 3;
      BFT_MALLOC(nomrub, length, char);
      sprintf(nomrub, "Parametres_real_ctwr_%02d", ct->num);

      BFT_MALLOC(tabvar, 4, cs_real_t);

      tabvar[ 0 ] = ct->cl_teau; /*  Water entry temperature*/
      tabvar[ 1 ] = ct->cl_fem;  /*  Water flow */
      tabvar[ 2 ] = ct->xap;     /* xap         */
      tabvar[ 3 ] = ct->xnp;     /* xnp         */


      nbvent  = 4;
      support = CS_MESH_LOCATION_NONE;
      typ_val = CS_TYPE_cs_real_t;

      cs_restart_write_section(suite,
                               nomrub,
                               support,
                               nbvent,
                               typ_val,
                               tabvar);

      BFT_FREE(tabvar);
    }


    { /* Write the temperature */
      char       nomrub[] = "Temperature_eau";

      typ_val = CS_TYPE_cs_real_t;
      nbvent  = 1;
      cs_restart_write_section(suite,
                               nomrub,
                               location_id,
                               nbvent,
                               typ_val,
                               ct->teau);

    }

    { /* Write the  */
      char       nomrub[] = "Flux_eau";

      typ_val = CS_TYPE_cs_real_t;
      nbvent  = 1;
      cs_restart_write_section(suite,
                               nomrub,
                               location_id,
                               nbvent,
                               typ_val,
                               ct->fem);

    }

    {/* Write the  */
      char       nomrub[] = "vitesse_goutte";

      typ_val = CS_TYPE_cs_real_t;
      nbvent  = 1;
      cs_restart_write_section(suite,
                               nomrub,
                               location_id,
                               nbvent,
                               typ_val,
                               ct->vgoutte);

    }

  }

  /* Close the restart file and free structures */
  cs_restart_destroy(&cs_glob_ctwr_suite);
}

/*----------------------------------------------------------------------------
 * Read the restart file of the cooling tower module
 *
 * Fortran interface:
 *
 * SUBROUTINE LECTWR
 * *****************
 *
 * character(kind=c_char)  nomsui : <-- : Name of the restart file
 *----------------------------------------------------------------------------*/

void CS_PROCF (lecctw, LECCTW)
(
 const char  *nomsui
)
{
  bool                corresp_cel, corresp_fac, corresp_fbr, corresp_som;
  cs_int_t            nbvent;
  cs_int_t            i, ict,indfac, ierror;

  cs_restart_t             *suite;
  cs_mesh_location_type_t   location_id,support;
  cs_restart_val_type_t     typ_val;

  cs_lnum_t n_g_elements, n_elements;

  cs_ctwr_zone_t  *ct;

  char        *location_name = NULL;
  cs_int_t     length        = 0;
  cs_gnum_t   *g_elt_num     = NULL;

  ierror = CS_RESTART_SUCCESS;

  /* Open the restart file */

  cs_glob_ctwr_suite
    = cs_restart_create(nomsui, NULL, CS_RESTART_MODE_READ);

  if (cs_glob_ctwr_suite == NULL)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while opening the cooling tower AERO module restart file "
                "in read mode.\n"
                "Verify the existence and the name of the restart file: %s\n"),
              nomsui);


  /* Pointer to the global restart structure */
  suite = cs_glob_ctwr_suite;

  /* Verification of the associated "support" to the restart file */
  cs_restart_check_base_location(suite, &corresp_cel, &corresp_fac,
                                 &corresp_fbr, &corresp_som);

  /* Only boundary faces are of interest */
  indfac = (corresp_fbr == true ? 1 : 0);
  if (indfac == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Abort while reading the 1D-wall thermal module restart file.\n"
                "The number of boundary faces has been modified\n"
                "Verify that the restart file corresponds to "
                "the present study.\n"));

  for (ict=0; ict < cs_glob_ct_nbr; ict++) {

    ct = cs_glob_ct_tab[ict];
    length = strlen("Cooling_Tower_restart_") + 3;
    BFT_MALLOC(location_name, length, char);
    sprintf(location_name, "Cooling_Tower_restart_%02d", ct->num);

    n_g_elements = fvm_nodal_get_n_g_elements(ct->water_mesh, FVM_CELL_HEXA);
    n_elements   = fvm_nodal_get_n_elements  (ct->water_mesh, FVM_CELL_HEXA);

    BFT_MALLOC(g_elt_num , n_g_elements, cs_gnum_t);

    fvm_nodal_get_global_element_num(ct->water_mesh ,
                                     FVM_CELL_HEXA  ,
                                     g_elt_num);

    location_id = cs_restart_add_location(suite,
                                          location_name,
                                          n_g_elements,
                                          n_elements,
                                          g_elt_num);


    {
      char * nomrub;
      cs_int_t   *tabvar;
      length = strlen("Parametres_int_ctwr_") + 3;
      BFT_MALLOC(nomrub, length, char);
      sprintf(nomrub, "Parametres_int_ctwr_%02d", ct->num);

      BFT_MALLOC(tabvar, 3, cs_int_t);

      nbvent  = 3;
      support = CS_MESH_LOCATION_NONE;
      typ_val = CS_TYPE_cs_int_t;

      ierror = cs_restart_read_section(suite,
                                       nomrub,
                                       support,
                                       nbvent,
                                       typ_val,
                                       tabvar);

      if (ierror < CS_RESTART_SUCCESS)
        bft_error(__FILE__, __LINE__, 0,
                  _("Problem while reading section in the restart file\n"
                    "for the cooling tower module:\n"
                    "<%s>\n"
                    "The calculation will not be run.\n"), nomrub);

      /* Coherency checks between the read   and the one from usctdz */

      if (tabvar[ 0 ] != ct->imctch) /* modele*/
        bft_printf(_("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                     "********               cooling tower MODULE\n"
                     "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                     "\n"
                     "The model is different \n"
                     "PREVIOUS: %d \n"
                     "CURRENT:  %d \n"), tabvar[ 0 ], ct->imctch);

      if (tabvar[ 1 ] != ct->ntypct) /* Type*/
        bft_error(  __FILE__, __LINE__, 0,
                     _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                       "********               cooling tower MODULE\n"
                       "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                       "\n"
                       "The type is different \n"
                       "PREVIOUS: %d \n"
                       "CURRENT:  %d \n"), tabvar[ 1 ], ct->ntypct);

      if (tabvar[ 2 ] != ct->nelect) /* nb of node per segment*/
        bft_error(__FILE__, __LINE__, 0,
                  _("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                    "********               cooling tower MODULE\n"
                    "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                    "\n"
                    "The number of nodes on each vertical mesh for \n"
                    "the water mesh has been modified.\n"
                    "PREVIOUS: %d nodes\n"
                    "CURRENT:  %d nodes\n"
                    "\n"
                    "The calculation will not be run.\n"
                    "\n"
                    "Verify that the restart file corresponds to a\n"
                    "restart file for the cooling tower  module.\n"
                    "Verify usctdz.\n"), tabvar[ 2 ], ct->nelect);

      BFT_FREE(tabvar);

    }

    {
      char * nomrub;
      cs_real_t   *tabvar;
      length = strlen("Parametres_real_ctwr_") + 3;
      BFT_MALLOC(nomrub, length, char);
      sprintf(nomrub, "Parametres_real_ctwr_%02d", ct->num);


      BFT_MALLOC(tabvar, 4, cs_real_t);

      nbvent  = 4;
      support = CS_MESH_LOCATION_NONE;
      typ_val = CS_TYPE_cs_real_t;

      ierror = cs_restart_read_section(suite,
                                       nomrub,
                                       support,
                                       nbvent,
                                       typ_val,
                                       tabvar);

      if (ierror < CS_RESTART_SUCCESS)
        bft_error(__FILE__, __LINE__, 0,
                  _("Problem while reading section in the restart file\n"
                    "for the cooling tower module:\n"
                    "<%s>\n"
                    "The calculation will not be run.\n"), nomrub);

      /* Coherency checks between the read   and the one from usctdz */

      if (CS_ABS(tabvar[ 0 ] - ct->cl_teau) > 1e-10)
        bft_printf(_("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                     "********               cooling tower MODULE\n"
                     "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                     "\n"
                     "The Water entry temperature  is different \n"
                     "PREVIOUS: %f \n"
                     "CURRENT:  %f \n"), tabvar[ 0 ], ct->cl_teau);

      if (CS_ABS(tabvar[ 1 ] - ct->cl_fem) > 1e-10)
        bft_printf(_("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                     "********               cooling tower MODULE\n"
                     "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                     "\n"
                     "The Water entry flow is different \n"
                     "PREVIOUS: %f \n"
                     "CURRENT:  %f \n"), tabvar[ 1 ], ct->cl_fem);

      if (CS_ABS(tabvar[ 2 ] - ct->xap) > 1e-10)
        bft_printf(_("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                     "********               cooling tower MODULE\n"
                     "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                     "\n"
                     "The value of Exchange law lambda coefficient is different \n"
                     "PREVIOUS: %f \n"
                     "CURRENT:  %f \n"), tabvar[ 2 ], ct->xap);

      if (CS_ABS(tabvar[ 3 ] - ct->xnp) > 1e-10)
        bft_printf(_("WARNING: ABORT WHILE READING THE RESTART FILE\n"
                     "********               cooling tower MODULE\n"
                     "       CURRENT AND PREVIOUS DATA ARE DIFFERENT\n"
                     "\n"
                     "The value of Exchange law lambda coefficient is different \n"
                     "PREVIOUS: %f \n"
                     "CURRENT:  %f \n"), tabvar[ 3 ], ct->xnp);

      BFT_FREE(tabvar);
    }


    { /* Read the wall thickness and check the coherency with USPT1D*/
      char        nomrub[] = "Temperature_eau";
      cs_real_t   *tabvar;

      BFT_MALLOC(tabvar, n_elements, cs_real_t);

      nbvent  = 1;
      typ_val = CS_TYPE_cs_real_t;

      ierror = cs_restart_read_section(suite,
                                       nomrub,
                                       location_id,
                                       nbvent,
                                       typ_val,
                                       tabvar);

      if (ierror < CS_RESTART_SUCCESS)
        bft_error(__FILE__, __LINE__, 0,
                  _("Problem while reading section in the restart file\n"
                    "for the cooling tower module:\n"
                    "<%s>\n"
                    "The calculation will not be run.\n"), nomrub);

      for (i = 0; i < n_elements; i++)
        ct->teau[i]= tabvar[i];

      BFT_FREE(tabvar);

    }

    { /* Read the wall thickness and check the coherency with USPT1D*/
      char        nomrub[] = "Flux_eau";
      cs_real_t   *tabvar;

      BFT_MALLOC(tabvar, n_elements, cs_real_t);

      nbvent  = 1;
      typ_val = CS_TYPE_cs_real_t;

      ierror = cs_restart_read_section(suite,
                                       nomrub,
                                       location_id,
                                       nbvent,
                                       typ_val,
                                       tabvar);

      if (ierror < CS_RESTART_SUCCESS)
        bft_error(__FILE__, __LINE__, 0,
                  _("Problem while reading section in the restart file\n"
                    "for the cooling tower module:\n"
                    "<%s>\n"
                    "The calculation will not be run.\n"), nomrub);

      for (i = 0; i < n_elements; i++)
        ct->fem[i]= tabvar[i];

      BFT_FREE(tabvar);
    }

    { /* Read the wall thickness and check the coherency with USPT1D*/
      char        nomrub[] = "vitesse_goutte";
      cs_real_t   *tabvar;


      BFT_MALLOC(tabvar, n_elements, cs_real_t);

      nbvent  = 1;
      typ_val = CS_TYPE_cs_real_t;

      ierror = cs_restart_read_section(suite,
                                       nomrub,
                                       location_id,
                                       nbvent,
                                       typ_val,
                                       tabvar);

      if (ierror < CS_RESTART_SUCCESS)
        bft_error(__FILE__, __LINE__, 0,
                  _("Problem while reading section in the restart file\n"
                    "for the cooling tower module:\n"
                    "<%s>\n"
                    "The calculation will not be run.\n"), nomrub);

      for (i = 0; i < n_elements; i++)
        ct->vgoutte[i]= tabvar[i];

      BFT_FREE(tabvar);
    }

  }

  /* Close the restart file and free structures */
  cs_restart_destroy(&cs_glob_ctwr_suite);
}

/*----------------------------------------------------------------------------
 * Post process variables associated with exchange area
 *
 * parameters:
 *   ct  <--  Void pointer to cooling tower function
 *   ts  <--  time step status structure, or NULL
 *----------------------------------------------------------------------------*/

static void
_cs_ctwr_post_function(void                  *ct,
                       const cs_time_step_t  *ts)
{
  const cs_ctwr_zone_t  *_ct = ct;

  if (_ct->post_mesh_id != 0) {

    cs_post_write_var(_ct->post_mesh_id,
                      _("T water"),
                      1,
                      false,
                      false,
                      CS_POST_TYPE_cs_real_t,
                      _ct->teau,
                      NULL,
                      NULL,
                      ts);

    cs_post_write_var(_ct->post_mesh_id,
                      _("Flux water"),
                      1,
                      false,
                      false,
                      CS_POST_TYPE_cs_real_t,
                      _ct->fem,
                      NULL,
                      NULL,
                      ts);

  }

}

/*============================================================================
 * Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Definition d'une zone d'echange (qui est ajoutee a celles deja definies)
 *----------------------------------------------------------------------------*/

void cs_ctwr_definit
(
  const int        idimct,    /* Dimemsion du probleme 2:2D  3:3D */
  const char      *ze_name,   /* Nom de la zone aero */
  const int        imctch,    /* 1: Modele de Poppe
                                 2: Merkel
                                 0: Rien */
  const int        ntypct,    /* 1: Contre courant
                                 2: Courant croises
                                 3: Zone de pluie */
  const cs_lnum_t  nelect,    /* Nombre d'elements sur chaque ligne du maillage
                                 eau pour la zone de noeuds par segment eau */
  const cs_real_t  deltat,    /* Ecart de temperature impose en entree de la
                                 zone d'echange */
  const cs_real_t  teau_cl,   /* Teau en entree de la zone d'echange */
  const cs_real_t  fem_cl,    /* debit en entree de la zone d'echange */
  const cs_real_t  xap,       /* coefficient lambda de la loi d'echange */
  const cs_real_t  xnp,       /* exposant n de la loi d'echange */
  const cs_real_t  surface,   /* Surface totale arrive d eau de la ct */
  const cs_real_t   dgout     /* Diametre de goutte pour les zones de pluie */
)
{
  cs_ctwr_zone_t  *ct;
  int length;
  FILE *f;
  char  *file_name = NULL;

  /* Definition d'une nouvelle zone d'echange */

  BFT_MALLOC(ct, 1, cs_ctwr_zone_t);

  ct->num = cs_glob_ct_nbr + 1;

  ct->idimct = idimct;
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

  ct->fem_e   = 0.0;
  ct->fem_s   = 0.0;
  ct->teau_e  = 0.0;
  ct->teau_s  = 0.0;
  ct->heau_e  = 0.0;
  ct->heau_s  = 0.0;
  ct->tair_e  = 0.0;
  ct->tair_s  = 0.0;
  ct->xair_e  = 0.0;
  ct->xair_s  = 0.0;
  ct->hair_e  = 0.0;
  ct->hair_s  = 0.0;
  ct->debit_e = 0.0;
  ct->debit_s = 0.0;

  ct->dgout = dgout;

  /* Selection des cellules */

  BFT_MALLOC(ct->ze_cell_list, cs_glob_mesh->n_b_faces, cs_lnum_t);

  cs_selector_get_cell_list(ze_name, &(ct->nbevct), ct->ze_cell_list);

  BFT_REALLOC(ct->ze_cell_list, ct->nbevct, cs_lnum_t);

  /* Redimensionnement du tableau des zones d'echange si necessaire */

  if (cs_glob_ct_nbr == cs_glob_ct_nbr_max) {
    cs_glob_ct_nbr_max = (cs_glob_ct_nbr_max + 1);
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
  ct->water_halo    = NULL;


  if (cs_glob_rank_id <= 0) {
    length = strlen("bltctc.") + 3;
    BFT_MALLOC(file_name, length, char);
    sprintf(file_name, "bltctc.%02d", ct->num);

    f = fopen(file_name, "a");

    fprintf(f, "# BILANS POUR LA ZONE D'ECHANGES \n");
    fprintf(f, "# ===============================\n");
    fprintf(f, "\tTEMP\tFLUX A/E\tTA MOY SOR\t TE MOY SOR");
    fprintf(f, "\tXA MOY SOR\tDEBI A ENT\tDEBI A SOR \n");
    fclose(f);
    BFT_FREE(file_name);
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

  for (i = 0; i < cs_glob_ct_nbr; i++) {

    ct = cs_glob_ct_tab[i];
    BFT_FREE(ct);

  }

  cs_glob_ct_nbr_max = 0;
  cs_glob_ct_nbr = 0;

  BFT_FREE(cs_stack_ct);
  BFT_FREE(cs_chain_ct);
  BFT_FREE(cs_glob_ctwr_props);

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
               cs_real_t   vitz[]       /* vitesse air suivant z */

)
{
  cs_lnum_t  ict,iseg,iloc,i, ii,j ,ieau, nb_dist_water, nb_dist_upw, ind;
  cs_real_t dhi,vvai,vhai,norme_g;
  cs_real_t gravite[3];
  cs_real_t faxn,bxan,xsata,xsate,cfen,ff1,ff2,xlew,eta,aux;
  cs_real_t vgin,dvg,cpx,rre,rpr,anu;
  cs_real_t   cpe, cpv, cpa, hv0, dgout, visc, conduc, rhoe;


  cs_lnum_t *lst_par_fac_sup_ct, *lst_par_fac_inf_ct_upw;
  const cs_lnum_t *locat_cel_upw = NULL;

  cs_ctwr_zone_t  *ct;
  cs_ctwr_zone_t  *ct_upw;
  cs_real_t *tai_inter, *xai_inter, *rhoai_inter,*vx_inter, *vy_inter,*vz_inter;
  cs_real_t *tai, *xai, *rhoai,*vx, *vy, *vz, *teau_upw_rec, *teau_upw_send;
  cs_real_t *fem_upw_rec, *fem_upw_send;
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;

  faxn = 0.;

  gravite[0] = -ct_prop->gravx;
  gravite[1] = -ct_prop->gravy;
  gravite[2] = -ct_prop->gravz;

  norme_g = sqrt( pow(gravite[0],2.)
                  +pow(gravite[1],2.)
                  +pow(gravite[2],2.));

  gravite[0] /= norme_g;
  gravite[1] /= norme_g;
  gravite[2] /= norme_g;


  /*--------------------------------------------*
   * Resolution des variable eau                *
   * sur chaque  ct                             *
   *--------------------------------------------*/
  for (ict=0; ict < cs_glob_ct_nbr; ict++) {

    ct = cs_glob_ct_tab[cs_chain_ct[ict]];

    if ((ct->ntypct>=2) && (ct->idimct==2))
     gravite[2] = 1.0;

    cpa    = ct_prop->cpa;
    cpv    = ct_prop->cpv;
    cpe    = ct_prop->cpe;
    hv0    = ct_prop->hv0;
    rhoe   = ct_prop->rhoe;
    visc   = ct_prop->visc;
    conduc = ct_prop->cond ;

    dgout  = ct->dgout;

    /*--------------------------------------------*
     * synchronisation   Halo                             *
     *--------------------------------------------*/

    if (ct->water_halo != NULL) {

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
    nb_dist_water = (int) ple_locator_get_n_dist_points(ct->locat_air_water);

    BFT_MALLOC(tai_inter  , nb_dist_water, cs_real_t);
    BFT_MALLOC(xai_inter  , nb_dist_water, cs_real_t);
    BFT_MALLOC(rhoai_inter, nb_dist_water, cs_real_t);
    BFT_MALLOC(vx_inter   , nb_dist_water, cs_real_t);
    BFT_MALLOC(vy_inter   , nb_dist_water, cs_real_t);
    BFT_MALLOC(vz_inter   , nb_dist_water, cs_real_t);

    for (ieau= 0; ieau < nb_dist_water; ieau++) {
       tai_inter[ieau]   = 0.;
       xai_inter[ieau]   = 0.;
       rhoai_inter[ieau] = 0.;
       vx_inter[ieau]    = 0.;
       vy_inter[ieau]    = 0.;
       vz_inter[ieau]    = 0.;
      for (i = (ct->pvoiseau[ieau]); i < (ct->pvoiseau[ieau+1]); i++) {
        tai_inter[ieau]  += ct->coefeau[i] * temp[ct->voiseau[i]];
        xai_inter[ieau]  += ct->coefeau[i] * xa[ct->voiseau[i]];
        rhoai_inter[ieau]+= ct->coefeau[i] * rho[ct->voiseau[i]];
        vx_inter[ieau]   += ct->coefeau[i] * vitx[ct->voiseau[i]];
        vy_inter[ieau]   += ct->coefeau[i] * vity[ct->voiseau[i]];
        vz_inter[ieau]   += ct->coefeau[i] * vitz[ct->voiseau[i]];
      }
    }
    BFT_MALLOC(tai  , ct->nnpsct*ct->nelect, cs_real_t);
    BFT_MALLOC(xai  , ct->nnpsct*ct->nelect, cs_real_t);
    BFT_MALLOC(rhoai, ct->nnpsct*ct->nelect, cs_real_t);
    BFT_MALLOC(vx   , ct->nnpsct*ct->nelect, cs_real_t);
    BFT_MALLOC(vy   , ct->nnpsct*ct->nelect, cs_real_t);
    BFT_MALLOC(vz   , ct->nnpsct*ct->nelect, cs_real_t);

    ple_locator_exchange_point_var(ct->locat_air_water,
                                   tai_inter, tai, NULL, sizeof(cs_real_t),1,0);
    ple_locator_exchange_point_var(ct->locat_air_water,
                                   xai_inter, xai, NULL, sizeof(cs_real_t),1,0);
    ple_locator_exchange_point_var(ct->locat_air_water,
                                   rhoai_inter,rhoai, NULL, sizeof(cs_real_t),1,0);
    ple_locator_exchange_point_var(ct->locat_air_water,
                                   vx_inter,vx, NULL, sizeof(cs_real_t),1,0);
    ple_locator_exchange_point_var(ct->locat_air_water,
                                   vy_inter,vy, NULL, sizeof(cs_real_t),1,0);
    ple_locator_exchange_point_var(ct->locat_air_water,
                                   vz_inter,vz, NULL, sizeof(cs_real_t),1,0);

    BFT_FREE(tai_inter);
    BFT_FREE(xai_inter);
    BFT_FREE(rhoai_inter);
    BFT_FREE(vx_inter);
    BFT_FREE(vy_inter);
    BFT_FREE(vz_inter);
    /*--------------------------------------------*
     *  end interpolation  air->eau               *
     *--------------------------------------------*/



   /*--------------------------------------------*
    * Calcul pour la face superieure ,           *
    * Introduction des conditions aux limites ct *
    *--------------------------------------------*/
    BFT_MALLOC(lst_par_fac_sup_ct , ct->nnpsct, cs_lnum_t);

    fvm_nodal_get_parent_num(ct->face_sup_mesh,
                                      2,lst_par_fac_sup_ct);
    ind = 0;
    for (j=0; j < cs_glob_ct_nbr; j++)
      if (cs_stack_ct[cs_chain_ct[ict]*cs_glob_ct_nbr + cs_chain_ct[j]] == 1) {
        ct_upw = cs_glob_ct_tab[ cs_chain_ct[j]];

        nb_dist_upw =
              (int)ple_locator_get_n_dist_points(ct->locat_cell_ct_upwind[ind]);

        BFT_MALLOC(teau_upw_send, nb_dist_upw, cs_real_t);
        BFT_MALLOC(fem_upw_send, nb_dist_upw, cs_real_t);
        BFT_MALLOC(lst_par_fac_inf_ct_upw,
                   (ct_upw->nbfac_ict+ct_upw->nbfbr_ict),
                   cs_lnum_t);

        fvm_nodal_get_parent_num(ct_upw->face_inf_mesh,
                                 2, lst_par_fac_inf_ct_upw);
        locat_cel_upw
          = ple_locator_get_dist_locations(ct->locat_cell_ct_upwind[ind]);

        for (i=0; i < nb_dist_upw; i++) {
          teau_upw_send[i] =  ct_upw->teau[(cs_lnum_t) locat_cel_upw[i]-1];
          fem_upw_send[i]  =  ct_upw->fem[(cs_lnum_t) locat_cel_upw[i]-1];
        }

        BFT_MALLOC(teau_upw_rec,
                   (ct_upw->nbfac_ict+ct_upw->nbfbr_ict), cs_real_t);
        BFT_MALLOC(fem_upw_rec,
                   (ct_upw->nbfac_ict+ct_upw->nbfbr_ict), cs_real_t);

        ple_locator_exchange_point_var(ct->locat_cell_ct_upwind[ind],
                                       teau_upw_send,
                                       teau_upw_rec,
                                       NULL,
                                       sizeof(cs_real_t),
                                       1,0);
        ple_locator_exchange_point_var(ct->locat_cell_ct_upwind[ind],
                                       fem_upw_send,
                                       fem_upw_rec,
                                       NULL,
                                       sizeof(cs_real_t),
                                       1,0);

        for (i=0; i < ct->nnpsct; i++) {
          ii = 0;
          while (ii < (ct_upw->nbfac_ict+ct_upw->nbfbr_ict)) {
            if (lst_par_fac_sup_ct[i] == lst_par_fac_inf_ct_upw[ii]) {
              ct->teau[i*ct->nelect] = teau_upw_rec[ii];
              ct->fem[i*ct->nelect]  = fem_upw_rec[ii];
              ii = ct_upw->nbfac_ict+ct_upw->nbfbr_ict;
            }
            ii++;
          }
        }
        BFT_FREE(teau_upw_rec);
        BFT_FREE(teau_upw_send);
        BFT_FREE(fem_upw_rec);
        BFT_FREE(fem_upw_send);
        BFT_FREE(lst_par_fac_inf_ct_upw);
        ind++;
      }

    BFT_FREE(lst_par_fac_sup_ct);

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

        for (iseg = 0; iseg < ct->nnpsct; iseg++) {
          /*--------------------------------------------*/
          /* Resolution Fe                              */
          /*--------------------------------------------*/
          for (iloc = 1; iloc < ct->nelect; iloc++) {

            ieau = iseg*ct->nelect + iloc;


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
            if (xai[ieau]>xsata) {
              aux=xsata;
            }else{
              aux=xai[ieau];
            }
            cfen=bxan*(xsate- aux)/(ct->fem[ieau]);
            ct->fem[ieau]=ct->fem[ieau-1]/(1.0-cfen*dhi);
          }
          /* Fin de resolution de Fe */

          /*--------------------------------------------*/
          /* Resolution Te                              */
          /*--------------------------------------------*/
          for (iloc = 1; iloc < ct->nelect; iloc++) {

            ieau = iseg*ct->nelect + iloc;

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
      else if (ct->ntypct==3) {
        for (iseg = 0; iseg < ct->nnpsct; iseg++) {
          /*--------------------------------------------*/
          /* Resolution Fe                              */
          /*--------------------------------------------*/
          for (iloc = 1; iloc < ct->nelect; iloc++) {

            ieau = iseg*ct->nelect + iloc;
            vgin=ct->vgoutte[ieau];

            if (CS_ABS(vgin)>=0.1) {

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
              if (xai[ieau]>xsata) {
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
          for (iloc = 1; iloc < ct->nelect; iloc++) {

            ieau = iseg*ct->nelect + iloc;

            vgin=ct->vgoutte[ieau];

            if (CS_ABS(vgin)>=0.1) {

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

      for (iseg = 0; iseg < ct->nnpsct; iseg++) {

        for (iloc = 1; iloc < ct->nelect; iloc++) {

          ieau = iseg*ct->nelect + iloc;

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

    else if (ct->ntypct==3) { /* zone de pluie */

      for (iseg = 0; iseg < ct->nnpsct; iseg++) {

        for (iloc = 1; iloc < ct->nelect; iloc++) {
          ieau = iseg*ct->nelect + iloc;
          ct->fem[ieau]=ct->fem[ieau-1];
          vgin=ct->vgoutte[ieau];

          if (CS_ABS(vgin)>=0.1) {

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

    BFT_FREE(tai);
    BFT_FREE(xai);
    BFT_FREE(rhoai);
    BFT_FREE(vx);
    BFT_FREE(vy);
    BFT_FREE(vz);

  }
  /*--------------------------------------------*/
  /* Fin de resolution des variable eau         */
  /* sur chaque  ct                             */
  /*--------------------------------------------*/
}

/*----------------------------------------------------------------------------
* Function cs_ctwr_aetssc
* Calcul des termes source pour l'air
*----------------------------------------------------------------------------*/
void cs_ctwr_aetssc
(
  int         iscal,       /*   */

  cs_real_t   temp[],      /* Temperature air */
  cs_real_t   xa[],        /* humidite air */
  cs_real_t   rho[],       /* masse volumique air */
  cs_real_t   utsim[],     /* vitesse verticale air */
  cs_real_t   utsex[],     /* vitesse horizontale air */
  cs_real_t   vitx[],      /* vitesse air suivant x */
  cs_real_t   vity[],      /* vitesse air suivant y */
  cs_real_t   vitz[]       /* vitesse air suivant z */
)
{
  cs_lnum_t  ict,iseg,iloc,ieau,iair,i,nb,nb_dist_water,nb_dist_air;
  cs_real_t cd1,ain,bin,dhi,dvga,gravite[3],norme_g;
  cs_real_t fax,fx0,vvai,vhai,tim,tex,xim,xex;
  cs_real_t bxa,xsata,xsate,ff1,xlew,eta;
  cs_real_t dvg,cpx,rre,rpr,anu;
  cs_real_t   cpe, cpv, cpa, hv0, dgout, visc, conduc, rhoe;
  cs_ctwr_zone_t  *ct;
  cs_real_t *tai_inter, *xai_inter, *rhoai_inter,*vx_inter, *vy_inter,*vz_inter,
            *tei_inter, *femei_inter, *vgin_inter;
  cs_real_t *tai, *xai, *rhoai,*vx, *vy, *vz, *tei, *femei, *vgin;
  cs_lnum_t  *lst_par_cel;
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;

  fax = 0.;

  gravite[0] = -ct_prop->gravx;
  gravite[1] = -ct_prop->gravy;
  gravite[2] = -ct_prop->gravz;

  norme_g = sqrt(  pow(gravite[0],2.)
                  +pow(gravite[1],2.)
                  +pow(gravite[2],2.));

  gravite[0] /= norme_g;
  gravite[1] /= norme_g;
  gravite[2] /= norme_g;


  /*--------------------------------------------*/
  /* Calcul de la vitesse des gouttes pour les  */
  /* zones de pluie                             */
  /*--------------------------------------------*/
  for (ict=0; ict < cs_glob_ct_nbr; ict++) {

    ct = cs_glob_ct_tab[cs_chain_ct[ict]];
    cpa    = ct_prop->cpa;
    cpv    = ct_prop->cpv;
    cpe    = ct_prop->cpe;
    hv0    = ct_prop->hv0;
    rhoe   = ct_prop->rhoe;
    visc   = ct_prop->visc;
    conduc = ct_prop->cond ;

    dgout  = ct->dgout;

    if (ct->ntypct==3) {
      /*--------------------------------------------*
      * synchronisation   Halo                            *
      *--------------------------------------------*/

      if (ct->water_halo != NULL) {

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
      nb_dist_water = (int) ple_locator_get_n_dist_points(ct->locat_air_water);

      BFT_MALLOC(tai_inter  , nb_dist_water, cs_real_t);
      BFT_MALLOC(xai_inter  , nb_dist_water, cs_real_t);
      BFT_MALLOC(rhoai_inter, nb_dist_water, cs_real_t);
      BFT_MALLOC(vx_inter   , nb_dist_water, cs_real_t);
      BFT_MALLOC(vy_inter   , nb_dist_water, cs_real_t);
      BFT_MALLOC(vz_inter   , nb_dist_water, cs_real_t);

      for (ieau= 0; ieau < nb_dist_water; ieau++) {
        tai_inter[ieau]   = 0.;
        xai_inter[ieau]   = 0.;
        rhoai_inter[ieau] = 0.;
        vx_inter[ieau]    = 0.;
        vy_inter[ieau]    = 0.;
        vz_inter[ieau]    = 0.;
        for (i = (ct->pvoiseau[ieau]); i < (ct->pvoiseau[ieau+1]); i++) {
          tai_inter[ieau]  += ct->coefeau[i] * temp[ct->voiseau[i]];
          xai_inter[ieau]  += ct->coefeau[i] * xa[ct->voiseau[i]];
          rhoai_inter[ieau]+= ct->coefeau[i] * rho[ct->voiseau[i]];
          vx_inter[ieau]   += ct->coefeau[i] * vitx[ct->voiseau[i]];
          vy_inter[ieau]   += ct->coefeau[i] * vity[ct->voiseau[i]];
          vz_inter[ieau]   += ct->coefeau[i] * vitz[ct->voiseau[i]];
        }
      }
      BFT_MALLOC(tai  , ct->nnpsct*ct->nelect, cs_real_t);
      BFT_MALLOC(xai  , ct->nnpsct*ct->nelect, cs_real_t);
      BFT_MALLOC(rhoai, ct->nnpsct*ct->nelect, cs_real_t);
      BFT_MALLOC(vx   , ct->nnpsct*ct->nelect, cs_real_t);
      BFT_MALLOC(vy   , ct->nnpsct*ct->nelect, cs_real_t);
      BFT_MALLOC(vz   , ct->nnpsct*ct->nelect, cs_real_t);

      ple_locator_exchange_point_var(ct->locat_air_water,
                                   tai_inter, tai, NULL, sizeof(cs_real_t),1,0);
      ple_locator_exchange_point_var(ct->locat_air_water,
                                   xai_inter, xai, NULL, sizeof(cs_real_t),1,0);
      ple_locator_exchange_point_var(ct->locat_air_water,
                                   rhoai_inter,rhoai, NULL, sizeof(cs_real_t),1,0);
      ple_locator_exchange_point_var(ct->locat_air_water,
                                   vx_inter,vx, NULL, sizeof(cs_real_t),1,0);
      ple_locator_exchange_point_var(ct->locat_air_water,
                                   vy_inter,vy, NULL, sizeof(cs_real_t),1,0);
      ple_locator_exchange_point_var(ct->locat_air_water,
                                   vz_inter,vz, NULL, sizeof(cs_real_t),1,0);
      /*--------------------------------------------*
      *  end interpolation  air->eau              *
      *--------------------------------------------*/


      dhi = -(ct->hmax-ct->hmin)/(ct->nelect-1);

      nb = (int) fvm_nodal_get_n_entities(ct->cell_mesh, 3);

      BFT_MALLOC(lst_par_cel , nb, cs_lnum_t);
      fvm_nodal_get_parent_num(ct->cell_mesh, 3, lst_par_cel);

      for (iseg = 0; iseg < ct->nnpsct; iseg++) {

        for (iloc = 1; iloc < ct->nelect; iloc++) {

          ieau = iseg*ct->nelect + iloc;

          vvai = sqrt(pow((vx[ieau]*gravite[0]),2.)
                   +pow((vy[ieau]*gravite[1]),2.)
                   +pow((vz[ieau]*gravite[2]),2.));
          /* fin interpolation air->eau */

          dvga = CS_ABS(ct->vgoutte[ieau]+vvai);

          rre  = dvga*rhoai[ieau]*dgout/visc;
          cd1 = (1.+0.15*pow(rre,0.687));
          ain = (18.*visc*cd1)/(rhoe*pow(dgout,2.));
          bin = -ain*dvga + 9.81;
          if (bin>0.) {
            ff1 = 2.*bin*dhi;
          }
          else {
            ff1 = 0.;
          }
          ct->vgoutte[ieau] = sqrt((pow(ct->vgoutte[ieau-1],2.)-ff1));
        }
      }
      BFT_FREE(lst_par_cel);
      BFT_FREE(tai);
      BFT_FREE(xai);
      BFT_FREE(rhoai);
      BFT_FREE(vx);
      BFT_FREE(vy);
      BFT_FREE(vz);
      BFT_FREE(tai_inter);
      BFT_FREE(xai_inter);
      BFT_FREE(rhoai_inter);
      BFT_FREE(vx_inter);
      BFT_FREE(vy_inter);
      BFT_FREE(vz_inter);
    }

  }  /* fin boucle ict sur les ct */
  /* Fin du calcul de la vitesse des gouttes pour les zones de pluie */

  /*--------------------------------------------*/
  /* Calcul des termes sources pour T et x      */
  /* pour chaque ct                             */
  /*--------------------------------------------*/
  for (ict=0; ict < cs_glob_ct_nbr; ict++) {
    ct = cs_glob_ct_tab[cs_chain_ct[ict]];

    if ((ct->ntypct >= 2) && (ct->idimct==2))
     gravite[2] = 1.0;

    cpa    = ct_prop->cpa;
    cpv    = ct_prop->cpv;
    cpe    = ct_prop->cpe;
    hv0    = ct_prop->hv0;
    rhoe   = ct_prop->rhoe;
    visc   = ct_prop->visc;
    conduc = ct_prop->cond ;

    dgout  = ct->dgout;

    nb = (int) fvm_nodal_get_n_entities(ct->cell_mesh, 3);

    /*--------------------------------------------*
    * synchronisation Halo                        *
    *--------------------------------------------*/

    if (ct->water_halo != NULL) {
      cs_halo_t *halo = ct->water_halo;
      cs_halo_sync_var(halo, ct->halo_type, ct->teau);
      cs_halo_sync_var(halo, ct->halo_type, ct->fem);
      cs_halo_sync_var(halo, ct->halo_type, ct->vgoutte);
    }


    BFT_MALLOC(lst_par_cel , nb, cs_lnum_t);
    fvm_nodal_get_parent_num(ct->cell_mesh, 3, lst_par_cel);
    /*--------------------------------------------*
     * interpolation  eau->air                    *
     *--------------------------------------------*/
    nb_dist_air = (int) ple_locator_get_n_dist_points(ct->locat_water_air);

    BFT_MALLOC(tei_inter   , nb_dist_air, cs_real_t);
    BFT_MALLOC(femei_inter , nb_dist_air, cs_real_t);
    BFT_MALLOC(vgin_inter  , nb_dist_air, cs_real_t);

    for (iair= 0; iair < nb_dist_air; iair++) {
       tei_inter  [ iair ] = 0.;
       femei_inter[ iair ] = 0.;
       vgin_inter [ iair ] = 0.;

      for (i = (ct->pvoisair[iair]); i < (ct->pvoisair[iair+1]); i++) {
        tei_inter[iair]    += ct->coefair[ i ]* ct->teau   [ ct->voisair[i] ];
        femei_inter[iair]  += ct->coefair[ i ]* ct->fem    [ ct->voisair[i] ];
        vgin_inter[iair]   += ct->coefair[ i ]* ct->vgoutte[ ct->voisair[i] ];
      }
    }
    BFT_MALLOC(tei   , ct->nbevct, cs_real_t);
    BFT_MALLOC(femei , ct->nbevct, cs_real_t);
    BFT_MALLOC(vgin  , ct->nbevct, cs_real_t);

    ple_locator_exchange_point_var(ct->locat_water_air,
                                   tei_inter,     tei, NULL, sizeof(cs_real_t),1,0);
    ple_locator_exchange_point_var(ct->locat_water_air,
                                   femei_inter, femei, NULL, sizeof(cs_real_t),1,0);
    ple_locator_exchange_point_var(ct->locat_water_air,
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
        for (iloc = 0; iloc < ct->nbevct; iloc++) {
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

        for (iloc = 0; iloc < ct->nbevct; iloc++) {
          iair = lst_par_cel[iloc]-1;

          if (CS_ABS(vgin[iloc])>=0.1) {
            vvai = sqrt(pow((vitx[iair]*gravite[0]),2.)
                       +pow((vity[iair]*gravite[1]),2.)
                       +pow((vitz[iair]*gravite[2]),2.));
            vhai = sqrt(pow((vitx[iair]*(1.-gravite[0])),2.)
                       +pow((vity[iair]*(1.-gravite[1])),2.)
                       +pow((vitz[iair]*(1.-gravite[2])),2.));
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
      /* Fin de la zone de pluie  */
    }
    /*--------------------------------------------*/
    /* Fin du modele de Poppe                     */
    /*--------------------------------------------*/

    /*--------------------------------------------*/
    /* Modele de Merkel                           */
    /*--------------------------------------------*/
    if (ct->imctch==2)  {
      if (ct->ntypct<=2) {
        for (iloc = 0; iloc < ct->nbevct; iloc++) {
          iair = lst_par_cel[iloc]-1;
          if (femei[iloc]>1.e-6) {
            vvai = sqrt(pow((vitx[iair]*gravite[0]),2.)
                       +pow((vity[iair]*gravite[1]),2.)
                       +pow((vitz[iair]*gravite[2]),2.));
            vhai = sqrt(pow((vitx[iair]*(1.-gravite[0])),2.)
                       +pow((vity[iair]*(1.-gravite[1])),2.)
                       +pow((vitz[iair]*(1.-gravite[2])),2.));

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

      else if (ct->ntypct==3) { /* zone de pluie */

        for (iloc = 0; iloc < ct->nbevct; iloc++) {
          iair = lst_par_cel[iloc]-1;

          if (CS_ABS(vgin[iloc])>=0.1) {

            vvai = sqrt(pow((vitx[iair]*gravite[0]),2.)
                       +pow((vity[iair]*gravite[1]),2.)
                       +pow((vitz[iair]*gravite[2]),2.));
            vhai = sqrt(pow((vitx[iair]*(1.-gravite[0])),2.)
                       +pow((vity[iair]*(1.-gravite[1])),2.)
                       +pow((vitz[iair]*(1.-gravite[2])),2.));

            dvg = sqrt(pow((vvai+vgin[iloc]),2.)+pow(vhai,2.));
            xsata = cs_ctwr_xsath(temp[iair]);
            xsate = cs_ctwr_xsath(tei[iloc]);
            cpx = cpa+xsata*cpv;
            rre = dvg*rho[iair]*(1.+xsata)*dgout/visc;
            rpr = cpx*visc/conduc;
            anu = 2.+0.6*sqrt(rre)*pow(rpr,(1./3.));
            bxa = (6.*conduc*anu*femei[iloc])/(0.92*rhoe*vgin[iloc]*pow(dgout,2.)*cpx);
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
      }/* Fin de la zone de pluie ntypct=3 */
    }
    /*--------------------------------------------*/
    /* Fin pour le modele de Merkel */
    /*--------------------------------------------*/
    BFT_FREE(lst_par_cel);
    BFT_FREE(tei_inter);
    BFT_FREE(femei_inter);
    BFT_FREE(vgin_inter);
    BFT_FREE(tei);
    BFT_FREE(femei);
    BFT_FREE(vgin);
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
  const int         idim,
  const cs_real_t   rho[],       /* masse volumique air */
  const cs_real_t   vitx[],      /* vitesse air suivant x */
  const cs_real_t   vity[],      /* vitesse air suivant y */
  const cs_real_t   vitz[],      /* vitesse air suivant z */
  const cs_real_t   xair[],      /* humidite de l'air */
  cs_real_t         utsex[]      /* terme source explicite */
)
{
  cs_lnum_t  ict, iloc, iair, i, *lst_par_cel, nb,nb_dist_air;
  cs_real_t dgout, visc, rhoe;
  cs_real_t absgrv, vginu,vginv,vginw,dvg,qer,rre,cdd1,cff0;

  cs_real_t  *femei_inter, *vgin_inter;
  cs_real_t  *femei, *vgin;
  cs_ctwr_zone_t  *ct;
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;

  absgrv = sqrt(pow(ct_prop->gravx,2.)+pow(ct_prop->gravy,2.)+pow(ct_prop->gravz,2.));

  /*--------------------------------------------*/
  /* Calcul de Kg pour chaque ct                */
  /*--------------------------------------------*/
  for (ict=0; ict < cs_glob_ct_nbr; ict++) {
    ct = cs_glob_ct_tab[cs_chain_ct[ict]];
    rhoe   = ct_prop->rhoe;
    dgout  = ct->dgout;
    visc   = ct_prop->visc;

     /*--------------------------------------------*
    * synchronisation Halo                        *
    *--------------------------------------------*/

    if (ct->water_halo != NULL) {
      cs_halo_t *halo = ct->water_halo;
      cs_halo_sync_var(halo, ct->halo_type, ct->teau);
      cs_halo_sync_var(halo, ct->halo_type, ct->fem);
      cs_halo_sync_var(halo, ct->halo_type, ct->vgoutte);
    }

    nb = (int) fvm_nodal_get_n_entities(ct->cell_mesh, 3);

    BFT_MALLOC(lst_par_cel , (nb*3), cs_lnum_t);
    fvm_nodal_get_parent_num(ct->cell_mesh, 3, lst_par_cel);
    /*--------------------------------------------*
     * interpolation  eau->air                    *
     *--------------------------------------------*/
    nb_dist_air = (int) ple_locator_get_n_dist_points(ct->locat_water_air);

    BFT_MALLOC(femei_inter  , nb_dist_air, cs_real_t);
    BFT_MALLOC(vgin_inter  , nb_dist_air, cs_real_t);

    for (iair= 0; iair < nb_dist_air; iair++) {

       femei_inter[iair] = 0.;
       vgin_inter[iair]  = 0.;

      for (i = (ct->pvoisair[iair]); i < (ct->pvoisair[iair+1]); i++) {

        femei_inter[ iair ]  += ct->coefair[ i ]* ct->fem    [ ct->voisair[i] ];
        vgin_inter [ iair ]  += ct->coefair[ i ]* ct->vgoutte[ ct->voisair[i] ];
      }
    }

    BFT_MALLOC(femei , ct->nbevct, cs_real_t);
    BFT_MALLOC(vgin , ct->nbevct, cs_real_t);

    ple_locator_exchange_point_var(ct->locat_water_air,
                                   femei_inter, femei, NULL, sizeof(cs_real_t),1,0);
    ple_locator_exchange_point_var(ct->locat_water_air,
                                   vgin_inter, vgin, NULL, sizeof(cs_real_t),1,0);

    /*--------------------------------------------*/
    /* zone de pluie                              */
    /*--------------------------------------------*/
    if (ct->ntypct==3)  {
      for (iloc = 0; iloc < ct->nbevct; iloc++) {
        iair = lst_par_cel[iloc]-1;

        vginu  = -ct_prop->gravx/ absgrv * vgin[iloc];
        vginv  = -ct_prop->gravy/ absgrv * vgin[iloc];
        vginw  = -ct_prop->gravz/ absgrv * vgin[iloc];
        dvg = sqrt(pow((vitx[iair]+vginu),2.)
              + pow((vity[iair]+vginv),2.)
              + pow((vitz[iair]+vginw),2.));
        if (vgin[iloc] > 0.1) {
          qer = femei[iloc]/rhoe;
          rre = dvg*rho[iair]*(1 + xair[iair])*dgout/visc;
          cdd1 = (1.+0.15*pow(rre,0.687));
          cff0 = 18.*cdd1*visc*qer/(vgin[iloc]*pow(dgout,2.));
          if (idim==1) {utsex[iair] = -cff0 *(vitx[iair]+vginu);}
          if (idim==2) {utsex[iair] = -cff0 *(vity[iair]+vginv);}
          if (idim==3) {utsex[iair] = -cff0 *(vitz[iair]+vginw);}
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

  const cs_mesh_t             *mesh,      /* <-- structure maillage associee  */
  const cs_mesh_quantities_t  *mesh_quantities   /* <-- grandeurs du maillage */
)
{
  const cs_real_t  *i_face_normal = mesh_quantities->i_face_normal;
  const cs_real_t  *b_face_normal = mesh_quantities->b_face_normal;
  cs_lnum_t         icel_1, icel_2, ict, idim, i, j, ieau_Sup, ieau_inf, length;
  cs_lnum_t         icel = -1, ifac = -1;
  cs_real_t        cpv, cpa, hv0;
  const cs_real_t  *coo_cen  = mesh_quantities->cell_cen;
  cs_real_t        debit,hair,n_sortant[3],vitair[3],aux,
                   surf,surf_e,surf_s;

  cs_lnum_t  *face_sup;      /* liste des faces  superieures de la ct */
  cs_lnum_t  *face_inf;      /* liste des faces  inferior de la ct */
  cs_lnum_t  *face_lat;      /* liste des faces  inferior de la ct */

  cs_ctwr_zone_t  *ct;
  FILE *f;
  char  *file_name = NULL;
  cs_ctwr_fluid_props_t  *ct_prop = cs_glob_ctwr_props;

  for (ict=0; ict < cs_glob_ct_nbr; ict++) {

    cs_lnum_t nbr_fbr_air[3][2];

    ct = cs_glob_ct_tab[cs_chain_ct[ict]];
    cpa    = ct_prop->cpa;
    cpv    = ct_prop->cpv;
    hv0    = ct_prop->hv0;

    nbr_fbr_air[0][0] = ct->nnpsct;
    nbr_fbr_air[0][1] = ct->nbfbr_sct;
    nbr_fbr_air[1][0] = ct->nbfbr_ict + ct->nbfac_ict;
    nbr_fbr_air[1][1] = ct->nbfbr_ict;
    nbr_fbr_air[2][0] = ct->nbfbr_lct + ct->nbfac_lct;
    nbr_fbr_air[2][1] = ct->nbfbr_lct;

    BFT_MALLOC(face_sup ,(ct->nbfac_sct + ct->nbfbr_sct) ,cs_lnum_t);
    fvm_nodal_get_parent_num(ct->face_sup_mesh, 2, face_sup);
    BFT_MALLOC(face_inf ,(ct->nbfac_ict + ct->nbfbr_ict) ,cs_lnum_t);
    fvm_nodal_get_parent_num(ct->face_inf_mesh, 2, face_inf);
    BFT_MALLOC(face_lat ,(ct->nbfbr_lct + ct->nbfac_lct) ,cs_lnum_t);
    fvm_nodal_get_parent_num(ct->face_lat_mesh, 2, face_lat);

    ct->fem_e   = 0.0;
    ct->fem_s   = 0.0;
    ct->teau_e  = 0.0;
    ct->heau_s  = 0.0;
    ct->heau_e  = 0.0;
    ct->teau_s  = 0.0;
    ct->tair_e  = 0.0;
    ct->tair_s  = 0.0;
    ct->xair_e  = 0.0;
    ct->xair_s  = 0.0;
    ct->hair_e  = 0.0;
    ct->hair_s  = 0.0;
    ct->debit_e = 0.0;
    ct->debit_s = 0.0;

    /* calcul des valeurs eau */

    for (i = 0; i < ct->nnpsct; i++) {
       ieau_Sup = i*ct->nelect;
       ieau_inf = (i+1)*ct->nelect - 1;

       surf = ct->surf_fac_sup[i];

       ct->teau_e += ct->teau[ieau_Sup]*ct->fem[ieau_Sup]*surf;
       ct->fem_e  += ct->fem[ieau_Sup]*surf;
       ct->heau_e += ct->teau[ieau_Sup]*ct->fem[ieau_Sup]*surf;

       ct->teau_s += ct->teau[ieau_inf]*ct->fem[ieau_inf]*surf;
       ct->fem_s  += ct->fem[ieau_inf]*surf;
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

    ct->teau_e /= ct->fem_e;
    ct->fem_e  /= ct->surface_in;
    ct->heau_e *= ct_prop->cpe;

    ct->teau_s /= ct->fem_s;
    ct->fem_s  /= ct->surface_out ;
    ct->heau_s *= ct_prop->cpe;

    /* calcul des valeurs air */

    surf_e = 0.;
    surf_s = 0.;

    for (j = 0; j < 3; j++)
    for (i = 0; i < nbr_fbr_air[j][0]; i++) {
      if (i< nbr_fbr_air[j][1]) {
        if (j==0) ifac = (cs_lnum_t) face_sup[i]-1;
        if (j==1) ifac = (cs_lnum_t) face_inf[i]-1;
        if (j==2) ifac = (cs_lnum_t) face_lat[i]-1;
        icel = mesh->b_face_cells[ifac];
        for (idim = 0; idim<3; idim++)
          n_sortant[idim] =  mesh_quantities->b_face_normal[ifac*3+idim];
        debit = CS_ABS(flux_masse_fbr[ifac]);
        surf  = cs_math_3_norm((b_face_normal + 3*ifac));
      } else {
        if (j==0) ifac = (cs_lnum_t) face_sup[i] - mesh->n_b_faces - 1;
        if (j==1) ifac = (cs_lnum_t) face_inf[i] - mesh->n_b_faces - 1;
        if (j==2) ifac = (cs_lnum_t) face_lat[i] - mesh->n_b_faces - 1;
        icel_1 = mesh->i_face_cells[ifac][0];
        icel_2 = mesh->i_face_cells[ifac][1];
        if (ct->mark_ze[icel_1] == 1) {

          icel = icel_2;
          for (idim = 0; idim < 3; idim++) {
            n_sortant[idim] =  coo_cen[icel_2*3 + idim] - coo_cen[icel_1*3 + idim];
          }
        }
        if (ct->mark_ze[icel_2] == 1) {

          icel = icel_1;
          for (idim = 0; idim < 3; idim++) {
            n_sortant[idim] =  coo_cen[icel_1*3 + idim] - coo_cen[icel_2*3 + idim];
          }
        }
        debit = CS_ABS(flux_masse_fac[ifac]);
        surf  = cs_math_3_norm((i_face_normal + 3*ifac));
      }
      hair = (cpa+xa[icel]*cpv)*temp[icel]+xa[icel]*hv0;
      vitair[0] = vitx[icel];
      vitair[1] = vity[icel];
      vitair[2] = vitz[icel];
      if (cs_math_3_dot_product(n_sortant, vitair)>0.) {
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

    if (CS_ABS(ct->debit_e) > 1e-10) {
      ct->tair_e /= ct->debit_e;
      ct->xair_e /= ct->debit_e;
    }

    if (CS_ABS(ct->debit_s) > 1e-10) {
      ct->tair_s /= ct->debit_s;
      ct->xair_s /= ct->debit_s;
    }

    fem_entree[ict]   = ct->fem_e;
    fem_sortie[ict]   = ct->fem_s;
    teau_entree[ict]  = ct->teau_e;
    teau_sortie[ict]  = ct->teau_s;
    heau_entree[ict]  = ct->heau_e;
    heau_sortie[ict]  = ct->heau_s;
    tair_entree[ict]  = ct->tair_e;
    tair_sortie[ict]  = ct->tair_s;
    xair_entree[ict]  = ct->xair_e;
    xair_sortie[ict]  = ct->xair_s;
    hair_entree[ict]  = ct->hair_e;
    hair_sortie[ict]  = ct->hair_s;

    ct->debit_e *= (ct->surface/ct->surface_in);
    ct->debit_s *= (ct->surface/ct->surface_out);

    debit_entree[ict] = ct->debit_e;
    debit_sortie[ict] = ct->debit_s;

    if (cs_glob_rank_id <= 0) {
      length = strlen("bltctc.") + 3;
      BFT_MALLOC(file_name, length, char);
      sprintf(file_name, "bltctc.%02d", ct->num);

      if (CS_ABS(ct->heau_e-ct->heau_s)> 1.e-6) {
        f = fopen(file_name, "a");

        aux = CS_ABS((ct->hair_s - ct->hair_e)/(ct->heau_e - ct->heau_s));
        fprintf(f, "%10f\t%10f\t%10f\t%10f\t%12.5e\t%10f\t%10f\n",
                time,
                aux,
                ct->tair_s,
                ct->teau_s,
                ct->xair_s,
                ct->debit_e,
                ct->debit_s);

        fclose(f);
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
cs_ctwr_post_init(int  ct_id,
                  int  writer_id)
{
  int  mesh_id = cs_post_get_free_mesh_id();
  int  writer_ids[] = {writer_id};

  cs_ctwr_zone_t * ct = cs_ctwr_by_id(ct_id);

  assert(ct != NULL);

  /* Exit silently if associated writer is not available */

  if (cs_post_writer_exists(writer_id) != true)
    return;

  /* Initialict post processing flag, and free previous arrays in
     case this function is called more than once */

  ct->post_mesh_id = mesh_id;

  /* Associate external mesh description with post processing subsystem */

  cs_post_define_existing_mesh(mesh_id,
                               ct->water_mesh,
                               0,
                               false,
                               false,
                               1,
                               writer_ids);

  /* Register post processing function */

  cs_post_add_time_dep_output(_cs_ctwr_post_function, (void *)ct);

  /* Update start and end (negative) numbers associated with
     dedicated post processing meshes */

  if (cs_glob_ct_post_mesh_ext[0] == 0)
    cs_glob_ct_post_mesh_ext[0] = mesh_id;

  cs_glob_ct_post_mesh_ext[1] = mesh_id;
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
cs_ctwr_by_id(int ct_id)
{
  cs_ctwr_zone_t  *retval = NULL;

  if (   ct_id > -1
      && ct_id <  cs_glob_ct_nbr)
    retval = cs_glob_ct_tab[ct_id];

  return retval;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
