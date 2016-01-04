#ifndef __CS_CTWR_H__
#define __CS_CTWR_H__

/*============================================================================
 * Main for cooling towers related functions
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

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_locator.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvm_nodal.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _cs_ctwr_zone_t cs_ctwr_zone_t;

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/* Cooling tower exchange zone structure definition */
/*--------------------------------------------------*/

struct _cs_ctwr_zone_t {

  int        idimct;          /* Problem dimension (2 or 3) */
  int        num;             /* Exchange zone number */
  char      *ze_name;         /* Exchange zone elements name */
  int        imctch;          /* 0: None; 1: Poppe's model; 2: Merkel's model */
  int        ntypct;          /* 1: Counter currents; 2: Crossed-currents;
                                 3: Rain zone */
  int        nelect;          /* Number of nodes on each vertical mesh for
                                 the water mesh */

  cs_real_t  hmin;            /* Minimum vertical height of exchange zone */
  cs_real_t  hmax;            /* Maximum height of exchange zone */
  cs_real_t  deltat;          /* Temperature delta required for exchange zone */

  cs_real_t  cl_teau;         /* Water entry temperature */
  cs_real_t  cl_fem;          /* Water flow */

  cs_real_t  xap;             /* Exchange law lambda coefficient */
  cs_real_t  xnp;             /* Exchange law n exponent */

  cs_real_t  surface_in;      /* Water inlet surface */
  cs_real_t  surface_out;     /* Water outlet surface */
  cs_real_t  surface;         /* Total surface */

  cs_int_t   nnpsct;          /* Number of points on top face */

  cs_int_t   nbfac_sct;       /* Number of top "interior faces" */
  cs_int_t   nbfac_ict;       /* Number of bottom "interior faces" */
  cs_int_t   nbfac_lct;       /* Number of lateral "interior faces" */
  cs_int_t   nbfac_ct;        /* Number of inside "interior faces" */
  cs_int_t   nbfbr_sct;       /* Number of top "boundary faces" */
  cs_int_t   nbfbr_ict;       /* Number of bottom "boundary faces" */
  cs_int_t   nbfbr_lct;       /* Number of lateral "boundary faces" */

  cs_int_t   nbevct;          /* Number of air cells belonging to the zone */

  cs_int_t   id_amont;        /* Number of upstream exchange zone (if any) */

  fvm_nodal_t *face_sup_mesh; /* Nodal Mesh of the zone's top faces */
  fvm_nodal_t *face_inf_mesh; /* Nodal mesh of the zone's bottom faces */
  fvm_nodal_t *face_lat_mesh; /* Nodal mesh of the zone's lateral faces */
  fvm_nodal_t *cell_mesh;     /* Nodal mesh of cells in the zone */
  fvm_nodal_t *fac_mesh;      /* Nodal mesh of internal faces in the zone */

  fvm_nodal_t *water_mesh;    /* Nodal mesh of water in the exchange area */

  cs_int_t   *ze_cell_list;        /* List of cells of ct criteria */
  cs_int_t   *voiseau;        /* List of water neighbors of air cells */
  cs_int_t   *pvoiseau;       /* Positions in the list of water neighbors */
  cs_int_t   *voisair;        /* List of air neighbors of water points */
  cs_int_t   *pvoisair;       /* Positions in the list of air neighbors */
  cs_int_t   *mark_ze;       /* Cell marker for ct */

  cs_int_t   *fac_sup_connect_idx; /* Top faces point connectivity index */
  cs_int_t   *fac_sup_connect_lst; /* Top faces point connectivity */

  cs_real_t  *surf_fac_sup;   /* Top faces surfaces */

  cs_real_t  *coefeau;        /* Water -> air interpolation coefficients */
  cs_real_t  *coefair;        /* Air -> water interpolation coefficients */

  cs_real_t  *teau;           /* Water temperature field */
  cs_real_t  *fem;            /* Water flow field */
  cs_real_t  *vgoutte;        /* Water drop velocity field (rain zones) */

  cs_real_t  fem_e;           /* Water entry flow */
  cs_real_t  fem_s;           /* Water exit flow */
  cs_real_t  teau_e;          /* Mean water entry temperature */
  cs_real_t  teau_s;          /* Mean water exit temperature */
  cs_real_t  heau_e;          /* Mean water entry enthalpy */
  cs_real_t  heau_s;          /* Mean water exit enthalpy */
  cs_real_t  tair_e;          /* Mean air entry temperature */
  cs_real_t  tair_s;          /* Mean air exit temperature */
  cs_real_t  xair_e;          /* Mean air entry humidity */
  cs_real_t  xair_s;          /* Mean air exit humidity */
  cs_real_t  hair_e;          /* Mean air entry enthalpy */
  cs_real_t  hair_s;          /* Mean air exit enthalpy */
  cs_real_t  debit_e;         /* Air entry flow */
  cs_real_t  debit_s;         /* Air exit flow */

  cs_real_t  dgout;           /* Drop diameter for rain zones */

  ple_locator_t   *locat_air_water; /* Locator water -> air interpolation */
  ple_locator_t   *locat_water_air; /* Locator for air -> water interpolation */
  ple_locator_t * *locat_cell_ct_upwind;

  cs_int_t   post_mesh_id;    /* 0 if post-processing is not active,
                                 mesh_id if post-processing is active */

  /* Parallelism and/or periodic features */

  cs_halo_type_t  halo_type;       /* Halo type */

  cs_int_t   *cs_array_rank;       /* Array of process ranks which
                                      have an exchange area */
  cs_int_t    nnpsct_with_ghosts;  /* Total number of water nodes
                                     (nnpsct + n_ghost_wcells) */

  cs_halo_t  *water_halo;          /* Structure used to manage ghost cells */

};

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* array of exchanges area */


extern cs_int_t            cs_glob_ct_nbr_max;
extern cs_int_t            cs_glob_ct_nbr;
extern cs_ctwr_zone_t     ** cs_glob_ct_tab;

/* array containing the stacking of the exchange area*/
extern cs_int_t  *  cs_stack_ct;

/* array containing the treatment order of the exchanges areas */
extern cs_int_t  *  cs_chain_ct;

/*============================================================================
 * Public function prototypes for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Define an exchange zone
 *
 * Fortran interface:
 *
 * SUBROUTINE DEFCT1
 * *****************
 *
 * INTEGER          IDIMCT        : --> : problem dimension (2 or 3)
 * INTEGER          IMCTCH        : --> : 1: Poppe's Model; 2: Merkel's model
 * INTEGER          NTYPCT        : --> : 1: Counter-currents
 *                                :     : 2: Crossed-currents
 *                                :     : 3: Rain zones
 * INTEGER          NELECT        : --> : number of nodes on each water line
 *                                :     : (i.e. vertical resolution)
 * DOUBLE PRECISION DELTAT        : --> : imposed water infow/outflow
 *                                : --> : temperature delta
 * DOUBLE PRECISION TEAU          : --> : water inlet mean temperature
 * DOUBLE PRECISION FEM           : --> : water inlet flow
 * DOUBLE PRECISION XAP           : --> : exchange law lambda coeffcient
 * DOUBLE PRECISION XNP           : --> : exchange law n exponent
 * DOUBLE PRECISION DGOUT         : --> : drop diameter for rain zones
 *----------------------------------------------------------------------------*/

void CS_PROCF (defct1, DEFCT1)
(
  const cs_int_t   *idimct,
  const char       *zecrit,
  cs_int_t         *ze_n_len,
  const cs_int_t   *imctch,
  const cs_int_t   *ntypct,
  const cs_int_t   *nelect,
  const cs_real_t  *deltat,
  const cs_real_t  *teau,
  const cs_real_t  *fem,
  const cs_real_t  *xap,
  const cs_real_t  *xnp,
  const cs_real_t  *surface,
  const cs_real_t  *dgout
);

/*----------------------------------------------------------------------------
 * Get number of cooling tower exchange zones.
 *
 * Fortran interface:
 *
 * SUBROUTINE NBZECT
 * *****************
 *
 * INTEGER          NBRCTZ        : --> : number of exchange zones
 *----------------------------------------------------------------------------*/

void CS_PROCF (nbzect, NBZECT)
(
 cs_int_t  *nbrctz
);

/*----------------------------------------------------------------------------
 * Indicate if the cooling tower model used is that of Poppe or Merkel.
 *
 * Fortran interface:
 *
 * SUBROUTINE AEMODE
 * *****************
 *
 * INTEGER          IMCTCH        : --> : model type (1: Poppe; 2: Merkel)
 *----------------------------------------------------------------------------*/

void CS_PROCF (aemode, AEMODE)
(
 cs_int_t  *imctch
);

/*----------------------------------------------------------------------------
 * Add a constant to the temperature vector for all exchange zones.
 *
 * Fortran interface:
 *
 * SUBROUTINE AEPROT
 * *****************
 *
 * DOUBLE PRECISION DELTA         : --> : type de ct (Poppe ou Merkel)
 *----------------------------------------------------------------------------*/

void CS_PROCF (aeprot, AEPROT)
(
 cs_real_t  *delta
);

/*----------------------------------------------------------------------------
 * Resolution of water variables
 *
 * Fortran interface:
 *
 * SUBROUTINE AETEAU ( )
 *
 * DOUBLE PRECISION TEMP(*)       : --> : air temperature
 * DOUBLE PRECISION XA(*)         : --> : air humidity
 * DOUBLE PRECISION RHO(*)        : --> : air density
 * DOUBLE PRECISION VITX(*)       : --> : air velocity component (x)
 * DOUBLE PRECISION VITY(*)       : --> : air velocity component (y)
 * DOUBLE PRECISION VITZ(*)       : --> : air velocity component (z)
 * DOUBLE PRECISION GX            : --> : gravity component x
 * DOUBLE PRECISION GY            : --> : gravity component y
 * DOUBLE PRECISION GZ            : --> : gravity component z
 *----------------------------------------------------------------------------*/

void CS_PROCF (aeteau, AETEAU)
(
  cs_real_t   temp[],
  cs_real_t   xa[],
  cs_real_t   rho[],
  cs_real_t   vitx[],
  cs_real_t   vity[],
  cs_real_t   vitz[]
);

/*----------------------------------------------------------------------------
 * Calculation of source terms for air equations
 *
 * Fortran interface:
 *
 * SUBROUTINE AETSSC
 * *****************
 *
 * INTEGER          ISCAL  : : scalar number
 * DOUBLE PRECISION TEMP   : : air temperature
 * DOUBLE PRECISION XA     : : air humidity
 * DOUBLE PRECISION RHO    : : air density
 * DOUBLE PRECISION UTSIM  : : implicite source term
 * DOUBLE PRECISION UTSEX  : : explicite source term
 * DOUBLE PRECISION VITX   : : air velocity along x
 * DOUBLE PRECISION VITY   : : air velocity along y
 * DOUBLE PRECISION VITZ   : : air velocity along z
 * DOUBLE PRECISION GX     : : x component of the gravity vector
 * DOUBLE PRECISION GY     : : y component of the gravity vector
 * DOUBLE PRECISION GZ     : : z component of the gravity vector
 *----------------------------------------------------------------------------*/

void CS_PROCF (aetssc, AETSSC)
(
  const cs_int_t   *iscal,
  cs_real_t         temp[],
  cs_real_t         xa[],
  cs_real_t         rho[],
  cs_real_t         utsim[],
  cs_real_t         utsex[],
  cs_real_t         vitx[],
  cs_real_t         vity[],
  cs_real_t         vitz[]
);

/*----------------------------------------------------------------------------
 * Calculation of induced head loss in rain zones
 *
 * Fortran interface:
 *
 * SUBROUTINE AETSVI
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF (aetsvi, AETSVI)
(
  const cs_int_t    *const idim,
  const cs_real_t   rho[],              /* masse volumique air */
  const cs_real_t   vitx[],             /* vitesse air suivant x */
  const cs_real_t   vity[],             /* vitesse air suivant y */
  const cs_real_t   vitz[],             /* vitesse air suivant z */
  const cs_real_t   xair[],             /* humidite de l'air */
  cs_real_t   utsex[]                   /* terme source explicite */
);

/*----------------------------------------------------------------------------
 * Bilan dans les ct
 *
 * Fortran interface:
 *
 * Interface Fortran :
 *
 * SUBROUTINE BILANct ( )
 *----------------------------------------------------------------------------*/

void CS_PROCF (bilanct, BILANCT)
(
  const cs_real_t   *const time,
  cs_real_t   fem_entree[],       /* debit eau entree */
  cs_real_t   fem_sortie[],       /* debit eau sortie */
  cs_real_t   teau_entree[],      /* temperature eau entree */
  cs_real_t   teau_sortie[],      /* temperature eau sortie */
  cs_real_t   heau_entree[],      /* enthalpie eau entree */
  cs_real_t   heau_sortie[],      /* enthalpie eau sortie */
  cs_real_t   tair_entree[],      /* temperature air entree */
  cs_real_t   tair_sortie[],      /* temperature air sortie */
  cs_real_t   xair_entree[],      /*  */
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
  const cs_real_t   vitz[]              /* vitesse air suivant z */
);

/*----------------------------------------------------------------------------
 * Initialize post processing.
 *
 * Fortran interface:
 *
 * Fortran Interface:
 *
 * SUBROUTINE PSTICT
 * *****************
 *----------------------------------------------------------------------------*/

void CS_PROCF(pstict, PSTICT)
(
 void
);

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
);

/*----------------------------------------------------------------------------
 * Read the restart file of the cooling tower module
 *
 * Fortran interface:
 *
 * subroutine lecctw
 * *****************
 *
 * character(kind=c_char)  nomsui : <-- : Name of the restart file
 *----------------------------------------------------------------------------*/

void CS_PROCF (lecctw, LECCTW)
(
 const char  *nomsui
);

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Definition d'une zone d'echange (qui est ajoute a celles deja definies)
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
  const cs_real_t  dgout      /* Diametre de goutte pour les zones de pluie */
);

/*----------------------------------------------------------------------------
 * Destruction des structures associees aux ventilateurs
 *----------------------------------------------------------------------------*/

void
cs_ctwr_all_destroy(void);

/*----------------------------------------------------------------------------
 * Resolution des variables eau
 *----------------------------------------------------------------------------*/

void cs_ctwr_aeteau
(
  cs_real_t   temp[],             /* Temperature air */
  cs_real_t   xa[],               /* humidite air */
  cs_real_t   rho[],              /* masse volumique air */
  cs_real_t   vitx[],             /* vitesse air suivant x */
  cs_real_t   vity[],             /* vitesse air suivant y */
  cs_real_t   vitz[]              /* vitesse air suivant z */
);

/*----------------------------------------------------------------------------
 * Calcul des termes source pour l'air
 *----------------------------------------------------------------------------*/

void cs_ctwr_aetssc
(
  int         iscal,               /*   */

  cs_real_t   temp[],             /* Temperature air */
  cs_real_t   xa[],               /* humidite air */
  cs_real_t   rho[],              /* masse volumique air */
  cs_real_t   utsim[],                  /* vitesse verticale air */
  cs_real_t   utsex[],                  /* vitesse horizontale air */
  cs_real_t   vitx[],             /* vitesse air suivant x */
  cs_real_t   vity[],             /* vitesse air suivant y */
  cs_real_t   vitz[]              /* vitesse air suivant z */
);

/*----------------------------------------------------------------------------
 * Calcul des PdC induites dans les zones de pluie
 *----------------------------------------------------------------------------*/

void cs_ctwr_aetsvi
(
  const int         idim,
  const cs_real_t   rho[],       /* masse volumique air */
  const cs_real_t   vitx[],      /* vitesse air suivant x */
  const cs_real_t   vity[],      /* vitesse air suivant y */
  const cs_real_t   vitz[],      /* vitesse air suivant z */
  const cs_real_t   xair[],             /* humidite de l'air */
  cs_real_t   utsex[]            /* terme source explicite */
);

/*----------------------------------------------------------------------------
 * Bilan dans les ct
 *----------------------------------------------------------------------------*/

void cs_ctwr_bilanct
(
  const cs_real_t   time,                /*   */
  cs_real_t   fem_entree[],             /* debit eau entree */
  cs_real_t   fem_sortie[],             /* debit eau sortie */
  cs_real_t   teau_entree[],            /* temperature eau entree */
  cs_real_t   teau_sortie[],            /* temperature eau sortie */
  cs_real_t   heau_entree[],            /* enthalpie eau entree */
  cs_real_t   heau_sortie[],            /* enthalpie eau sortie */
  cs_real_t   tair_entree[],            /* temperature air entree */
  cs_real_t   tair_sortie[],            /* temperature air sortie */
  cs_real_t   xair_entree[],            /*  */
  cs_real_t   xair_sortie[],            /*   */
  cs_real_t   hair_entree[],            /*   */
  cs_real_t   hair_sortie[],            /*   */
  cs_real_t   debit_entree[],           /*   */
  cs_real_t   debit_sortie[],           /*   */

  const cs_real_t   temp[],             /* Temperature air */
  const cs_real_t   xa[],               /* humidite air */
  const cs_real_t   flux_masse_fac[],   /* vitesse verticale air */
  const cs_real_t   flux_masse_fbr[],   /* vitesse horizontale air */
  const cs_real_t   vitx[],             /* vitesse air suivant x */
  const cs_real_t   vity[],             /* vitesse air suivant y */
  const cs_real_t   vitz[],             /* vitesse air suivant z */
  const cs_mesh_t      *mesh,      /* <-- structure maillage associee  */
  const cs_mesh_quantities_t  *mesh_quantities   /* <-- grandeurs du maillage */
);

/*----------------------------------------------------------------------------
 * Initialict post-processing
 *
 * parameters:
 *   ct_id         -->  Id of exchange area
 *   writer_id           -->  Id of associated writer
 *----------------------------------------------------------------------------*/

void
cs_ctwr_post_init(int  ct_id,
                  int  writer_id);

/*----------------------------------------------------------------------------
 * Get pointer to exchange area.
 *
 * parameters:
 *   ct_id  <--  Id (0 to n-1) of exchange area
 *
 * returns:
 *   pointer to exchange area structure
 *----------------------------------------------------------------------------*/

cs_ctwr_zone_t *
cs_ctwr_by_id(int ct_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CTWR_H__ */
