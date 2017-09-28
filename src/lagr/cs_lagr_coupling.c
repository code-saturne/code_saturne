/*============================================================================
 * Methods for particle coupling
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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
 * Functions dealing with lagrangian coupling
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

#include "cs_base.h"
#include "cs_math.h"
#include "cs_prototypes.h"

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_physical_constants.h"

#include "cs_lagr.h"
#include "cs_lagr_particle.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_coupling.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Constants */

static const cs_real_t  _c_stephan = 5.6703e-8;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/* --------------------------------------------------------------------
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
 * --------------------------------------------------------------------   */

void
cs_lagr_coupling(cs_real_t taup[],
                 cs_real_t tempct[],
                 cs_real_t tsfext[],
                 cs_real_t cpgd1[],
                 cs_real_t cpgd2[],
                 cs_real_t cpght[],
                 cs_real_t volp[],
                 cs_real_t volm[])
{
  cs_real_3_t *st_vel = NULL, *t_st_vel = NULL;
  cs_real_6_t *st_rij = NULL, *t_st_rij = NULL;

  /* ====================================================================   */
  /* 1. INITIALISATION    */
  /* ====================================================================   */

  {
    cs_field_t *f = cs_field_by_name_try("velocity_st_lagr");
    if (f != NULL)
      st_vel = (cs_real_3_t *)(f->val);
  }

  {
    cs_field_t *f = cs_field_by_name_try("rij_st_lagr");
    if (f != NULL)
      st_rij = (cs_real_6_t *)(f->val);
  }

  cs_lagr_extra_module_t *extra = cs_glob_lagr_extra_module;
  cs_lagr_source_terms_t *lag_st = cs_glob_lagr_source_terms;

  cs_real_3_t grav    = {cs_glob_physical_constants->gravity[0],
                         cs_glob_physical_constants->gravity[1],
                         cs_glob_physical_constants->gravity[2]};

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lnum_t ncelet = cs_glob_mesh->n_cells_with_ghosts;
  cs_lnum_t ncel = cs_glob_mesh->n_cells;
  cs_lnum_t nbpart = p_set->n_particles;

  cs_real_t dtp = cs_glob_lagr_time_step->dtp;

  cs_lnum_t ntersl = cs_glob_lagr_dim->ntersl;

  cs_real_t *tslag, *auxl1, *auxl2, *auxl3;
  BFT_MALLOC (tslag, ncelet * ntersl, cs_real_t);
  BFT_MALLOC (auxl1, nbpart, cs_real_t);
  BFT_MALLOC (auxl2, nbpart, cs_real_t);
  BFT_MALLOC (auxl3, nbpart, cs_real_t);

  /*   Nombre de passage pour les termes sources en stationnaire  */
  if (   cs_glob_lagr_time_scheme->isttio == 1
      && cs_glob_time_step->nt_cur >= cs_glob_lagr_source_terms->nstits)
    lag_st->npts += 1;

  cs_glob_lagr_source_terms->ntxerr = 0;
  cs_glob_lagr_source_terms->vmax = 0.0;
  cs_glob_lagr_source_terms->tmamax = 0.0;

  for (cs_lnum_t iel = 0; iel < ncel; iel++) {
    volp[iel]      = 0.0;
    volm[iel]      = 0.0;
  }

  for (cs_lnum_t ivar = 0; ivar < ntersl; ivar++) {

    for (cs_lnum_t iel = 0; iel < ncel; iel++)
      tslag[ncelet * ivar + iel]  = 0.0;

  }

  /* ====================================================================   */
  /* 2. CALCULS PRELIMINAIRES  */
  /* ====================================================================   */
  /* Finalisation des forces externes (Si la particule a interagit avec     */
  /* une frontiere du domaine de calcul, on degenere a l'ordre 1).     */

  for (cs_lnum_t npt = 0; npt < nbpart; npt++) {

    cs_real_t aux1 = dtp / taup[npt];
    cs_real_t p_mass = cs_lagr_particles_get_real(p_set, npt, CS_LAGR_MASS);

    if (   cs_glob_lagr_time_scheme->t_order == 1
        || cs_lagr_particles_get_lnum(p_set, npt, CS_LAGR_REBOUND_ID) == 0)
      tsfext[npt] = (1.0 - exp(-aux1)) * p_mass * taup[npt];

    else
      tsfext[npt] +=  (1.0 - (1.0 - exp (-aux1)) / aux1) * taup[npt]
                    * p_mass;

  }

  for (cs_lnum_t npt = 0; npt < nbpart; npt++) {

    cs_real_t  p_stat_w = cs_lagr_particles_get_real(p_set, npt, CS_LAGR_STAT_WEIGHT);
    cs_real_t  p_mass   = cs_lagr_particles_get_real(p_set, npt, CS_LAGR_MASS);
    cs_real_t *p_vel    = cs_lagr_particles_attr(p_set, npt, CS_LAGR_VELOCITY);

    cs_real_t  prev_p_mass = cs_lagr_particles_get_real_n(p_set, npt, 1, CS_LAGR_MASS);
    cs_real_t *prev_p_vel  = cs_lagr_particles_attr_n(p_set, npt, 1, CS_LAGR_VELOCITY);

    auxl1[npt] = p_stat_w * (p_mass * p_vel[0] - prev_p_mass * prev_p_vel[0]
                                               - grav[0] * tsfext[npt]) / dtp;
    auxl2[npt] = p_stat_w * (p_mass * p_vel[1] - prev_p_mass * prev_p_vel[1]
                                               - grav[1] * tsfext[npt]) / dtp;
    auxl3[npt] = p_stat_w * (p_mass * p_vel[2] - prev_p_mass * prev_p_vel[2]
                                               - grav[2] * tsfext[npt]) / dtp;

  }

  /* ====================================================================   */
  /* 3. TERMES SOURCES DE QUANTITE DE MOUVEMENT    */
  /* ====================================================================   */

  if (cs_glob_lagr_source_terms->ltsdyn == 1) {

    if (cs_glob_lagr_time_scheme->isttio == 1 && lag_st->npts > 0)
      BFT_MALLOC(t_st_vel, ncel, cs_real_3_t);
    else
      t_st_vel = st_vel;

    for (cs_lnum_t i = 0; i < ncel; i++) {
      for (cs_lnum_t j = 0; j < 3; j++)
        t_st_vel[i][j] = 0;
    }

    for (cs_lnum_t npt = 0; npt < nbpart; npt++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

      cs_real_t  p_stat_w = cs_lagr_particle_get_real(particle, p_am,
                                                      CS_LAGR_STAT_WEIGHT);

      cs_real_t  prev_p_diam = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                           CS_LAGR_DIAMETER);
      cs_real_t  prev_p_mass = cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                           CS_LAGR_MASS);
      cs_real_t  p_mass = cs_lagr_particle_get_real(particle, p_am,
                                                    CS_LAGR_MASS);

      cs_lnum_t iel = cs_lagr_particle_get_cell_id(particle, p_am);

      /* Volume et masse des particules dans la maille */
      volp[iel] += p_stat_w * cs_math_pi * pow(prev_p_diam, 3) / 6.0;
      volm[iel] += p_stat_w * prev_p_mass;

      /* TS de QM   */
      t_st_vel[iel][0] += - auxl1[npt];
      t_st_vel[iel][1] += - auxl2[npt];
      t_st_vel[iel][2] += - auxl3[npt];
      tslag[iel + (lag_st->itsli-1) * ncelet] += - 2.0 * p_stat_w * p_mass / taup[npt];

    }

    /* ====================================================================   */
    /* 4. TERMES SOURCES SUR LA TURBULENCE */
    /* ====================================================================   */

    if (extra->itytur == 2 || extra->iturb == 50 || extra->iturb == 60) {

      /* En v2f (ITURB=50) les TS lagrangiens influent uniquement sur k et eps  */
      /* (difficile d'ecrire quoi que ce soit sur v2, qui perd son sens de */
      /*  "composante de Rij")     */

      for (cs_lnum_t npt = 0; npt < nbpart; npt++) {

        unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

        cs_lnum_t  iel         = cs_lagr_particle_get_cell_id(particle, p_am);
        cs_real_t *prev_f_vel  = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                         CS_LAGR_VELOCITY_SEEN);
        cs_real_t *f_vel       = cs_lagr_particle_attr(particle, p_am,
                                                       CS_LAGR_VELOCITY_SEEN);

        cs_real_t uuf = 0.5 * (prev_f_vel[0] + f_vel[0]);
        cs_real_t vvf = 0.5 * (prev_f_vel[1] + f_vel[1]);
        cs_real_t wwf = 0.5 * (prev_f_vel[2] + f_vel[2]);

        tslag[iel + (lag_st->itske-1) * ncelet] += - uuf * auxl1[npt]
                                                   - vvf * auxl2[npt]
                                                   - wwf * auxl3[npt];

      }

      for (cs_lnum_t iel = 0; iel < ncel; iel++)
        tslag[iel + (lag_st->itske-1) * ncelet]
          += - extra->vel->val[iel * 3    ] * t_st_vel[iel][0]
             - extra->vel->val[iel * 3 + 1] * t_st_vel[iel][1]
             - extra->vel->val[iel * 3 + 2] * t_st_vel[iel][2];

    }
    else if (extra->itytur == 3) {

      if (cs_glob_lagr_time_scheme->isttio == 1 && lag_st->npts > 0)
        BFT_MALLOC(t_st_rij, ncel, cs_real_6_t);
      else
        t_st_rij = st_rij;

      for (cs_lnum_t i = 0; i < ncel; i++) {
        for (cs_lnum_t j = 0; j < 6; j++)
          t_st_rij[i][j] = 0;
      }

      for (cs_lnum_t npt = 0; npt < nbpart; npt++) {

        unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

        cs_lnum_t  iel         = cs_lagr_particle_get_cell_id(particle, p_am);

        cs_real_t *prev_f_vel  = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                         CS_LAGR_VELOCITY_SEEN);
        cs_real_t *f_vel       = cs_lagr_particle_attr(particle, p_am,
                                                       CS_LAGR_VELOCITY_SEEN);

        cs_real_t uuf = 0.5 * (prev_f_vel[0] + f_vel[0]);
        cs_real_t vvf = 0.5 * (prev_f_vel[1] + f_vel[1]);
        cs_real_t wwf = 0.5 * (prev_f_vel[2] + f_vel[2]);

        t_st_rij[iel][0] += - 2.0 * uuf * auxl1[npt];
        t_st_rij[iel][1] += - 2.0 * vvf * auxl2[npt];
        t_st_rij[iel][2] += - 2.0 * wwf * auxl3[npt];
        t_st_rij[iel][3] += - uuf * auxl2[npt] - vvf * auxl1[npt];
        t_st_rij[iel][4] += - vvf * auxl3[npt] - wwf * auxl2[npt];
        t_st_rij[iel][5] += - uuf * auxl3[npt] - wwf * auxl1[npt];

      }
      for (cs_lnum_t iel = 0; iel < ncel; iel++) {

        t_st_rij[iel][0] += - 2.0 * extra->vel->val[iel * 3    ]
                                  * t_st_vel[iel][0];

        t_st_rij[iel][1] += - 2.0 * extra->vel->val[iel * 3 + 1]
                                  * t_st_vel[iel][1];

        t_st_rij[iel][2] += - 2.0 * extra->vel->val[iel * 3 + 2]
                                  * t_st_vel[iel][2];

        t_st_rij[iel][3] += - extra->vel->val[iel * 3    ] * t_st_vel[iel][1]
                            - extra->vel->val[iel * 3 + 1] * t_st_vel[iel][0];

        t_st_rij[iel][4] += - extra->vel->val[iel * 3 + 1] * t_st_vel[iel][2]
                            - extra->vel->val[iel * 3 + 2] * t_st_vel[iel][1];

        t_st_rij[iel][5] += - extra->vel->val[iel * 3    ] * t_st_vel[iel][2]
                            - extra->vel->val[iel * 3 + 2] * t_st_vel[iel][0];

      }

    }

  }

  /* ====================================================================   */
  /* 5. TERME SOURCE MASSIQUES */
  /* ====================================================================   */

  if (    cs_glob_lagr_source_terms->ltsmas == 1
      && (   cs_glob_lagr_specific_physics->impvar == 1
          || cs_glob_lagr_specific_physics->idpvar == 1)) {

    for (cs_lnum_t npt = 0; npt < nbpart; npt++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

      cs_real_t  p_stat_w = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);
      cs_real_t  prev_p_mass = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_MASS);
      cs_real_t  p_mass = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_MASS);

      /* Dans saturne TSmasse > 0 ===> Apport de masse sur le fluide  */
      cs_lnum_t iel = cs_lagr_particle_get_cell_id(particle, p_am);

      tslag[iel + (lag_st->itsmas-1) * ncelet] += - p_stat_w * (p_mass - prev_p_mass) / dtp;

    }

  }

  /* ====================================================================   */
  /* 6. TERMES SOURCES THERMIQUE    */
  /* ====================================================================   */

  if (cs_glob_lagr_source_terms->ltsthe == 1) {

    if (   cs_glob_lagr_model->physical_model == 1
        && cs_glob_lagr_specific_physics->itpvar == 1) {

      for (cs_lnum_t npt = 0; npt < nbpart; npt++) {

        unsigned char *particle = p_set->p_buffer + p_am->extents * npt;
        cs_lnum_t iel = cs_lagr_particle_get_cell_id(particle, p_am);
        cs_real_t  p_mass = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_MASS);
        cs_real_t  prev_p_mass = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_MASS);
        cs_real_t  p_cp = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_CP);
        cs_real_t  prev_p_cp = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_CP);
        cs_real_t  p_tmp = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_TEMPERATURE);
        cs_real_t  prev_p_tmp = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_TEMPERATURE);
        cs_real_t  p_stat_w = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);

        tslag[iel + (lag_st->itste-1) * ncelet] += - (p_mass * p_tmp * p_cp
                                - prev_p_mass * prev_p_tmp * prev_p_cp) / dtp * p_stat_w;
        tslag[iel + (lag_st->itsti-1) * ncelet] += tempct[nbpart + npt] * p_stat_w;

      }
      if (extra->radiative_model > 0) {

        for (cs_lnum_t npt = 0; npt < nbpart; npt++) {

          unsigned char *particle = p_set->p_buffer + p_am->extents * npt;
          cs_lnum_t iel = cs_lagr_particle_get_cell_id(particle, p_am);
          cs_real_t  p_diam = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_DIAMETER);
          cs_real_t  p_eps = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_EMISSIVITY);
          cs_real_t  p_tmp = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_TEMPERATURE);
          cs_real_t  p_stat_w = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);

          cs_real_t aux1 = cs_math_pi * p_diam * p_diam * p_eps
                          * (extra->luminance->val[iel] - 4.0 * _c_stephan * pow (p_tmp, 4));

          tslag[iel + (lag_st->itste-1) * ncelet] += aux1 * p_stat_w;

        }

      }

    }
    else if (cs_glob_lagr_model->physical_model == 2) {
      if (cs_glob_lagr_const_dim->nlayer > 1)
        bft_error(__FILE__, __LINE__, 0,
                  _("Couplage thermique non-fonctionnel en multi-layer"));

      else {

        for (cs_lnum_t npt = 0; npt < nbpart; npt++) {

          unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

          cs_lnum_t iel  = cs_lagr_particle_get_cell_id(particle, p_am);
          cs_lnum_t icha = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_COAL_NUM);

          cs_real_t  p_mass = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_MASS);
          cs_real_t  p_tmp = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_TEMPERATURE);
          cs_real_t  p_cp = cs_lagr_particle_get_real_n(particle, p_am, 0, CS_LAGR_CP);

          cs_real_t  prev_p_mass = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_MASS);
          cs_real_t  prev_p_tmp  = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_TEMPERATURE);
          cs_real_t  prev_p_cp   = cs_lagr_particle_get_real_n(particle, p_am, 1, CS_LAGR_CP);

          cs_real_t  p_stat_w = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);

          tslag[iel + (lag_st->itste-1) * ncelet]  += - (  p_mass * p_tmp * p_cp
                                             - prev_p_mass * prev_p_tmp * prev_p_cp) / dtp * p_stat_w;
          tslag[iel + (lag_st->itsti-1) * ncelet] += tempct[nbpart + npt] * p_stat_w;
          tslag[iel + (lag_st->itsmv1[icha]-1) * ncelet] += p_stat_w * cpgd1[npt];
          tslag[iel + (lag_st->itsmv2[icha]-1) * ncelet] += p_stat_w * cpgd2[npt];
          tslag[iel + (lag_st->itsco-1) * ncelet] += p_stat_w * cpght[npt];
          tslag[iel + (lag_st->itsfp4-1) * ncelet] = 0.0;

        }

      }

    }

  }

  /* ====================================================================   */
  /* 7. Verif que le taux volumique maximal TVMAX admissible de particules  */
  /*    ne soit pas depasse dans quelques cellules.     */
  /* ====================================================================   */

  cs_real_t *st_val = cs_glob_lagr_source_terms->st_val; /* short alias */

  const cs_real_t tvmax = 0.8;
  const cs_real_t *cell_vol = cs_glob_mesh_quantities->cell_vol;

  for (cs_lnum_t iel = 0; iel < ncel; iel++) {

    cs_real_t mf   = cell_vol[iel] * extra->cromf->val[iel];
    cs_real_t tauv = volp[iel] / cell_vol[iel];
    cs_real_t taum = volm[iel] / mf;

    if (tauv > tvmax) {

      cs_glob_lagr_source_terms->ntxerr++;;

      for (int ivar = 0; ivar < ntersl; ivar++)
        st_val[iel + ivar * ncelet] = 0.0;

      if (t_st_vel != NULL) {
        for (cs_lnum_t j = 0; j < 3; j++)
          t_st_vel[iel][j] = 0.0;
      }

      if (t_st_rij != NULL) {
        for (cs_lnum_t j = 0; j < 6; j++)
          t_st_rij[iel][j] = 0.0;
      }

    }

    cs_glob_lagr_source_terms->vmax
      = CS_MAX(tauv, cs_glob_lagr_source_terms->vmax);
    cs_glob_lagr_source_terms->tmamax
      = CS_MAX(taum, cs_glob_lagr_source_terms->tmamax);

  }

  /* ====================================================================   */
  /* 8. MOYENNE TEMPORELLE DES TERMES SOURCES */
  /* ====================================================================   */

  if (cs_glob_lagr_time_scheme->isttio == 1 && lag_st->npts > 0) {

    for (int ivar = 0; ivar < ntersl; ivar++) {

      for (cs_lnum_t iel = 0; iel < ncel; iel++)
        st_val[iel + ivar * ncelet]
          =  (  tslag[iel + ncelet * ivar]
              + (lag_st->npts - 1.0) * st_val[iel + ncelet * ivar])
            / lag_st->npts;

    }

    if (st_vel != NULL) {
      for (cs_lnum_t iel = 0; iel < ncel; iel++) {
        for (cs_lnum_t j = 0; j < 3; j++) {
          st_vel[iel][j]
            =    (t_st_vel[iel][j] + (lag_st->npts - 1.0) * st_vel[iel][j])
               / lag_st->npts;
        }
      }
    }

    if (st_rij != NULL) {
      for (cs_lnum_t iel = 0; iel < ncel; iel++) {
        for (cs_lnum_t j = 0; j < 6; j++) {
          st_rij[iel][j]
            =    (t_st_rij[iel][j] + (lag_st->npts - 1.0) * st_rij[iel][j])
               / lag_st->npts;
        }
      }
    }

  }
  else {

    for (int ivar = 0; ivar < ntersl; ivar++) {
      for (cs_lnum_t iel = 0; iel < ncel; iel++)
        st_val[iel + ncelet * ivar] = tslag[iel + ncelet * ivar];
    }

  }

  if (t_st_vel != st_vel)
    BFT_FREE(t_st_vel);

  if (t_st_rij != st_rij)
    BFT_FREE(t_st_rij);

  BFT_FREE(auxl1);
  BFT_FREE(auxl2);
  BFT_FREE(auxl3);
  BFT_FREE(tslag);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
