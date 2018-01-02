/*============================================================================
 * Handling of new particles.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*============================================================================
 * Functions dealing with particle tracking
 *============================================================================*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_periodicity.h"

#include "cs_base.h"
#include "cs_halo.h"
#include "cs_log.h"
#include "cs_interface.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_adjacencies.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_physical_model.h"
#include "cs_prototypes.h"
#include "cs_random.h"
#include "cs_search.h"
#include "cs_timer_stats.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_lagr_clogging.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_dlvo.h"
#include "cs_lagr_stat.h"
#include "cs_lagr_geom.h"
#include "cs_lagr.h"
#include "cs_lagr_tracking.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_new.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */
/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Inject a series of particles at random positions on a given zone.
 *
 * \warning Currently works only for tri and quadrangular faces.
 *
 * The fluid velocity and other variables and attributes are computed here.
 *
 * \param[inout]   npt      current number of particles
 * \param[in]      nznew    number of added particles
 * \param[in]      zone_id  zone number (zone index + 1 in arrays)
 * \param[in]      ifrlag   boundary zone number for lagrangian
 * \param[in,out]  iworkp   array containing injection face number
 *----------------------------------------------------------------------------*/

void
cs_lagr_new(cs_lnum_t  *npt,
            cs_lnum_t   nznew,
            int         zone_id,
            const int   ifrlag[],
            cs_lnum_t   iworkp[])
{
  cs_real_t eps = 0.001;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq  = cs_glob_mesh_quantities;

  cs_real_t *xyzcen = fvq->cell_cen;
  cs_real_t *xyznod = mesh->vtx_coord;

  cs_lagr_particle_set_t  *particles = NULL;
  particles = cs_glob_lagr_particle_set;

  const cs_lnum_t n_b_faces = mesh->n_b_faces;

  /* CALCUL DE LA SURFACE MAX DE L'ENTREE :   */
  cs_real_t surfm  = -10.0;
  cs_lnum_t minfac = n_b_faces;
  cs_lnum_t maxfac = -1 ;

  for (cs_lnum_t ifac = 0; ifac < n_b_faces; ifac++) {

    if (ifrlag[ifac] == zone_id) {

      surfm  = CS_MAX(surfm, fvq->b_face_surf[ifac]);
      minfac = CS_MIN(ifac, minfac);
      maxfac = CS_MAX(ifac, maxfac);

    }

  }

  /* STOP PAR SECURITE    */
  if (maxfac == 0 || minfac == n_b_faces + 1) {
    /* write(nfecra,9000) zone_id   */
    cs_exit (1);
  }

  cs_lnum_t n_vertices;
  cs_lnum_t ifac;
  cs_real_t random;

  cs_real_t ctr[18], vec[3];
  int       iconfo[4];

  /* BOUCLE SUR LES NOUVELLES PARTICULES :    */
  for (cs_lnum_t newpart = 0; newpart < nznew; newpart++) {

    /* tirage aleatoire d'une face :  */
    bool goon = true;
    while (goon) {

      cs_random_uniform(1, &random);
      random = random * (cs_real_t) ((maxfac - minfac + 1) - eps);

      ifac    = minfac + (int) (random);

      if ( ifac >= minfac && ifac <= maxfac) {

        if (ifrlag[ifac] == zone_id) {

          /* tirage aleatoire pour determiner si cette face convient */
          /* plus la fa7 est grande plus elle a des chance d'etre choisie */
          cs_random_uniform(1, &random);

          if (random <= (fvq->b_face_surf[ifac] / surfm)) {

            /* ATTENTION :     */
            /* type de face : 3 ou 4 points supports    */
            /* pour l'instant je ne sais pas traiter les autres   */
            /* avec plus de points supports...     */
            n_vertices = mesh->b_face_vtx_idx[ifac + 1] - mesh->b_face_vtx_idx[ifac];
            if (n_vertices <= 4)
              goon = false;

          }

        }

      }

    }

    /* si face a 4 points, on choisit l'un des deux triangles  */
    if (n_vertices == 4) {

      n_vertices = 0;

      for (cs_lnum_t i = mesh->b_face_vtx_idx[ifac]; i < mesh->b_face_vtx_idx[ifac + 1] ; i++) {

        iconfo[n_vertices] = mesh->b_face_vtx_lst[i] ;
        n_vertices         = n_vertices + 1;

      }

      cs_real_t surftr[2];
      cs_real_t are[6];
      /* longueur des arretes 1 et 2 du premier triangle :  */
      for (cs_lnum_t i = 0; i < 2; i++){
        are[i*3 + 0] = xyznod[3 * iconfo[1 + i]    ] - xyznod[3 * iconfo[0]    ] ;
        are[i*3 + 1] = xyznod[3 * iconfo[1 + i] + 1] - xyznod[3 * iconfo[0] + 1] ;
        are[i*3 + 2] = xyznod[3 * iconfo[1 + i] + 2] - xyznod[3 * iconfo[0] + 2] ;
      }

      /* surface du premier triangle    */
      vec[0]    = are[1] * are[5] - are[2] * are[4];
      vec[1]    = are[2] * are[3] - are[0] * are[5];
      vec[2]    = are[0] * are[4] - are[1] * are[3];
      surftr[0] = sqrt (pow (vec[0], 2) + pow (vec[1], 2) + pow (vec[2], 2));

      /* longueur des arretes 1 et 2 du deuxieme triangle : */
      for (cs_lnum_t i = 0; i < 2; i++){
        are[i*3 + 0] = xyznod[3 * iconfo[2+i]    ] - xyznod[3 * iconfo[0]    ] ;
        are[i*3 + 1] = xyznod[3 * iconfo[2+i] + 1] - xyznod[3 * iconfo[0] + 1] ;
        are[i*3 + 2] = xyznod[3 * iconfo[2+i] + 2] - xyznod[3 * iconfo[0] + 2] ;
      }
      for (cs_lnum_t i = 0; i < 4; i++)

      /* surface du deuxieme triangle   */
      vec[0]    = are[1] * are[5] - are[2] * are[4];
      vec[1]    = are[2] * are[3] - are[0] * are[5];
      vec[2]    = are[0] * are[4] - are[1] * are[3];
      surftr[1] = sqrt (vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

      /* tirage d'un nombre aleatoire entre 0 et 1 */
      /* pour determiner quel triangle choisir    */
      cs_random_uniform(1, &random);

      /* si le deuxieme triangle est choisit, on reorganise */
      /* les points : 4 <--> 2     */
      if (random <= (surftr[1] / (surftr[0] + surftr[1]))) {

        cs_lnum_t i  = iconfo[3];
        iconfo[3]  = iconfo[1];
        iconfo[1]  = i;

      }

    }

    /* dans le cas ou la face est un triangle...     */
    else if (n_vertices == 3) {

      n_vertices = 0;

      for (cs_lnum_t i = mesh->b_face_vtx_idx[ifac];
           i < mesh->b_face_vtx_idx[ifac + 1];
           i++) {

        iconfo[n_vertices] = mesh->b_face_vtx_lst[i];
        n_vertices = n_vertices + 1;

      }

    }

    /* constitution des coordonnees du triangle */
    for (cs_lnum_t i = 0; i < 3; i++) {

      ctr[i * 3 + 0] = xyznod[3 * iconfo[i]    ];
      ctr[i * 3 + 1] = xyznod[3 * iconfo[i] + 1];
      ctr[i * 3 + 2] = xyznod[3 * iconfo[i] + 2];

    }

    /* tirage aleatoire d'un point dans le triangle  */
    /* constitue des points 1,2,3     */
    /* 6 est dans le triangle si PM1*PM6>=0     */
    cs_real_t pm1  =  -1.0;
    cs_real_t pm6  = 1.0;
    while (pm1 * pm6 < 0.0) {

      /* 1) tirage du point 4 sur l'arete 12 */
      random = 0.0;
      while (random <= 0.0 || random >= 1.0)
        cs_random_uniform(1, &random);

      for (cs_lnum_t i = 0; i < 3; i++)
        ctr[3 * 3 + i]     = random * ctr[i] + (1.0 - random) * ctr[3 + i];

      /* 2) tirage du point 5 sur l'arete 13 */
      random = 0.0;
      while (random <= 0.0 || random >= 1.0)
        cs_random_uniform(1, &random);

      for (cs_lnum_t i = 0; i < 3; i++)
        ctr[12 + i]     = random * ctr[i] + (1.0 - random) * ctr[6 + i];


      /* 3) le point 6 est le sommet du parallelogramme 1465     */
      for (cs_lnum_t i = 0; i < 3; i++)
        ctr[15 + i]    = ctr[9 + i] + ctr[12 + i] - ctr[i];

      /* 4) reste a verifier que le point 6 appartient au triangle 123     */
      /* 4.1) vecteur normal au triangle : 12^13  */
      vec[0]  = (ctr[4] - ctr[1]) * (ctr[8] - ctr[2])
              - (ctr[5] - ctr[2]) * (ctr[7] - ctr[1]);
      vec[1]  = (ctr[5] - ctr[2]) * (ctr[6] - ctr[0])
              - (ctr[3] - ctr[0]) * (ctr[8] - ctr[2]);
      vec[2]  = (ctr[3] - ctr[0]) * (ctr[7] - ctr[1])
              - (ctr[4] - ctr[1]) * (ctr[6] - ctr[0]);

      /* 4.2) produit mixte pour le point 1 :     */
      pm1     = 0.0;
      pm1     = pm1 + vec[0] * ((ctr[4] - ctr[1]) * (ctr[8] - ctr[5])
                              - (ctr[5] - ctr[2]) * (ctr[7] - ctr[4]));
      pm1     = pm1 + vec[1] * ((ctr[5] - ctr[2]) * (ctr[6] - ctr[3])
                              - (ctr[3] - ctr[0]) * (ctr[8] - ctr[5]));
      pm1     = pm1 + vec[2] * ((ctr[3] - ctr[0]) * (ctr[7] - ctr[4])
                              - (ctr[4] - ctr[1]) * (ctr[6] - ctr[3]));

      /* 4.3) produit mixte pour le point 6 :     */
      pm6     = 0.0;
      pm6     = pm6 + vec[0] * ((ctr[4] - ctr[16]) * (ctr[8] - ctr[5])
                              - (ctr[5] - ctr[17]) * (ctr[7] - ctr[4]));
      pm6     = pm6 + vec[1] * ((ctr[5] - ctr[17]) * (ctr[6] - ctr[3])
                              - (ctr[3] - ctr[15]) * (ctr[8] - ctr[5]));
      pm6     = pm6 + vec[2] * ((ctr[3] - ctr[15]) * (ctr[7] - ctr[4])
                              - (ctr[4] - ctr[16]) * (ctr[6] - ctr[3]));
    }

    /* 5) POUR PLUS DE SECURITE, ON DEPLACE LE POINT */
    /* D'UN EPSILON EN DIRECTION DU CENTRE CELLULE   */
    /* ATTENTION : CE DECALAGE PEUT ETRE DANGEREUX DANS LE CAS */
    /* DE CELULES CONCAVES  */
    cs_lnum_t cell = mesh->b_face_cells[ifac];
    ctr[15] = ctr[15] + (xyzcen[3 * cell    ] - ctr[15]) * eps;
    ctr[16] = ctr[16] + (xyzcen[3 * cell + 1] - ctr[16]) * eps;
    ctr[17] = ctr[17] + (xyzcen[3 * cell + 2] - ctr[17]) * eps;

    /* LE TRAITEMENT EST TERMINE POUR LE POINT NPT,  */
    /* ON REMPLIT LES TABLEAUX POUR LE LAGRANGIEN :  */
    cs_real_t *part_coord = cs_lagr_particles_attr(particles, *npt, CS_LAGR_COORDS);
    part_coord[0] = ctr[15];
    part_coord[1] = ctr[16];
    part_coord[2] = ctr[17];

    cs_lagr_particles_set_lnum(particles, *npt, CS_LAGR_CELL_NUM, cell + 1);

    iworkp[*npt] = ifac ;

    /* incrementation du pointeur sur les particules */
    *npt = *npt + 1;

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Initialization for new particles.
 *
 * The fluid velocity seen is computed here.
 *
 * \param[in]   p_id_l     lower particle id bound (included)
 * \param[in]   p_id_u     uppder particle id bound (excluded)
 * \param[in]   time_id    associated time id (0: current, 1: previous)
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_new_particle_init(cs_lnum_t   p_id_l,
                          cs_lnum_t   p_id_u,
                          cs_lnum_t   time_id,
                          cs_real_t   vislen[])
{
  cs_lagr_particle_set_t  *pset = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = pset->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  cs_lagr_bdy_condition_t  *_bdy_conditions = cs_lagr_get_bdy_conditions();
  int *ifrlag = _bdy_conditions->b_face_zone_id;

  const cs_lnum_t  n_cells = cs_glob_mesh->n_cells;
  const cs_lnum_t  ncelet = cs_glob_mesh->n_cells_with_ghosts;

  /* Map field arrays */

  cs_real_3_t *vel;
  cs_real_t *cvar_k, *cvar_r11, *cvar_r22, *cvar_r33;
  cs_real_t *cvar_rij;

  vel = (cs_real_3_t *)extra->vel->vals[time_id];

  if (   extra->itytur == 2
      || extra->iturb == 50
      || extra->iturb == 60) {

    cvar_k = extra->cvar_k->vals[time_id];

  }

  else if (extra->itytur == 3) {

    if (extra->cvar_rij == NULL) {
      cvar_r11 = extra->cvar_r11->vals[time_id];
      cvar_r22 = extra->cvar_r22->vals[time_id];
      cvar_r33 = extra->cvar_r33->vals[time_id];
    } else {
      cvar_rij = extra->cvar_rij->vals[time_id];
    }

  }

  /* Initialization */

  cs_real_t d2s3 = 2.0 / 3.0;

  /* Allocate a work array     */
  cs_real_t *w1 = NULL;
  BFT_MALLOC(w1, ncelet, cs_real_t);

  /* ==============================================================================
   * 2. SIMULATION DES VITESSES TURBULENTES FLUIDES INSTANTANNEES VUES
   *    PAR LES PARTICULES SOLIDES LE LONG DE LEUR TRAJECTOIRE.
   * ============================================================================== */

  if (cs_glob_lagr_time_scheme->idistu == 1) {

    if (   extra->itytur == 2
        || extra->iturb  == 50
        || extra->iturb  == 60) {

      for (cs_lnum_t iel = 0; iel < n_cells; iel++)
        w1[iel] = cvar_k[iel];

    }

    else if (extra->itytur == 3) {

      if (extra->cvar_rij == NULL) {
        for (cs_lnum_t iel = 0; iel < n_cells; iel++)
          w1[iel] = 0.5 * (cvar_r11[iel] + cvar_r22[iel] + cvar_r33[iel]);
      } else {
        for (cs_lnum_t iel = 0; iel < n_cells; iel++)
          w1[iel] = 0.5*(cvar_rij[6*iel] + cvar_rij[6*iel+1] + cvar_rij[6*iel+2]);
      }

    }

    else {
      cs_log_printf(CS_LOG_DEFAULT, "   %d %d %d\n",
                    cs_glob_lagr_time_scheme->iilagr,
                    cs_glob_lagr_time_scheme->idistu,
                    extra->iturb);
      bft_error
        (__FILE__, __LINE__, 0,
         _("The Lagrangian module is incompatible with the selected\n"
           " turbulence model.\n\n"
           "Turbulent dispersion is used with:\n"
           "  cs_glob_lagr_time_scheme->idistu = %d\n"
           "And the turbulence model is iturb = %d\n\n"
           "The only turbulence models compatible with the Lagrangian model's\n"
           "turbulent dispersion are k-epsilon, Rij-epsilon, v2f, and k-omega."),
         cs_glob_lagr_time_scheme->idistu,
         extra->iturb);
    }

  }

  else {

    for (cs_lnum_t iel = 0; iel < n_cells; iel++)
      w1[iel] = 0.0;

  }

  /* --> CALCUL DES TIRAGES ALEATOIRES   */
  /*     CALCUL DU TEMPS CARACTERISTIQUE DES PARTICULES */

  cs_lnum_t nomb = p_id_u - p_id_l;
  cs_real_t *vagaus[3];

  for (int i = 0; i < 3; i++)
    BFT_MALLOC(vagaus[i], nomb, cs_real_t);

  if (cs_glob_lagr_time_scheme->idistu == 1 && nomb > 0) {

    cs_random_normal(nomb, vagaus[0]);
    cs_random_normal(nomb, vagaus[1]);
    cs_random_normal(nomb, vagaus[2]);

  }

  else {

    for (cs_lnum_t npt = 0; npt < nomb; npt++) {
      vagaus[0][npt] = 0.0;
      vagaus[1][npt] = 0.0;
      vagaus[2][npt] = 0.0;
    }

  }

  for (cs_lnum_t npt = p_id_l; npt < p_id_u; npt++) {

    unsigned char *particle = pset->p_buffer + p_am->extents * npt;

    cs_lnum_t iel  = cs_lagr_particle_get_cell_id(particle, p_am);

    cs_real_t  *particle_velocity_seen
      = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

    cs_real_t tu = sqrt(d2s3 * w1[iel]);

    for (cs_lnum_t i = 0; i < 3; i++)
      particle_velocity_seen[i] = vel[iel][i] + vagaus[i][npt - p_id_l] * tu;

    cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_REBOUND_ID, -1);

    cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TR_TRUNCATE, 0);

  }

  for (int i = 0; i < 3; i++)
    BFT_FREE(vagaus[i]);

  /* Calcul de la fluctuation de vitesse si le modèle de dépôt est activé */

  if (cs_glob_lagr_model->deposition == 1) {

    const cs_mesh_adjacencies_t  *ma = cs_glob_mesh_adjacencies;

    for (cs_lnum_t npt = p_id_l; npt < p_id_u; npt++) {

      unsigned char *particle = pset->p_buffer + p_am->extents * npt;

      cs_lnum_t iel = cs_lagr_particle_get_cell_id(particle, p_am);

      /* Calculation of the normalized wall-normal particle distance (y^+)  */

      cs_real_t yplus  = 1000.0;
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_YPLUS, yplus);

      for (cs_lnum_t il = ma->cell_b_faces_idx[iel];
           il < ma->cell_b_faces_idx[iel+1];
           il++) {

        cs_lnum_t ifac = ma->cell_b_faces[il];

        int zone_id = ifrlag[ifac];

        /* Test if the particle is located in a boundary cell */

        if (   _bdy_conditions->b_zone_natures[zone_id] == CS_LAGR_DEPO1
            || _bdy_conditions->b_zone_natures[zone_id] == CS_LAGR_DEPO2
            || _bdy_conditions->b_zone_natures[zone_id] == CS_LAGR_DEPO_DLVO
            || _bdy_conditions->b_zone_natures[zone_id] == CS_LAGR_REBOUND) {

          /* Calculation of the wall units  */

          cs_lnum_t  *neighbor_face_id;
          cs_real_t  *particle_yplus;

          if (cs_glob_lagr_model->deposition > 0) {
            neighbor_face_id
              = cs_lagr_particle_attr(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID);
            particle_yplus
              = cs_lagr_particle_attr(particle, p_am, CS_LAGR_YPLUS);
          }
          else {
            neighbor_face_id = NULL;
            particle_yplus = 0;  /* allow tests even without particle y+ */
          }

          _test_wall_cell(particle, p_am, vislen,
                          particle_yplus, neighbor_face_id);

        }

      }

      if (yplus < cs_lagr_particle_get_real(particle, p_am, CS_LAGR_INTERF)) {

        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_MARKO_VALUE, 10);

      }

      else if (yplus > 100.0) {

        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_MARKO_VALUE, -1);

      }

      else {

        cs_real_t random;
        cs_random_uniform(1, &random);

        if (random < 0.25) {

          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_MARKO_VALUE, 12);

        }

        else if (random > 0.625) {

          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_MARKO_VALUE, 1);

        }

        else { //if ((random > 0.25) && (random < 0.625)) { JB: problem if random=0.25 or 0.625 !

          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_MARKO_VALUE, 3);

        }

      }

      if (yplus <= cs_lagr_particle_get_real(particle, p_am, CS_LAGR_INTERF)) {

        cs_real_t *vel_seen
          = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);

        for (cs_lnum_t i = 0; i < 3; i++)
          vel_seen[i] = vel[iel][i];

      }

      /* No deposited particles at the injection */
      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                CS_LAGR_PART_IN_FLOW);

      /* Initialization of additional "pointers"
       * for the resuspension model              */

      if (cs_glob_lagr_model->resuspension > 0) {

        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_N_LARGE_ASPERITIES, 0);
        cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_N_SMALL_ASPERITIES, 0);
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DISPLACEMENT_NORM, 0.0);

      }

    }

  }

  /* Free memory     */
  BFT_FREE(w1);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
