/*============================================================================
 * Methods for particle precipitation
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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
 * Functions dealing with particle precipitation
 *============================================================================*/

#include "cs_defs.h"
#include "cs_math.h"

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
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_order.h"
#include "cs_parall.h"
#include "cs_random.h"
#include "cs_search.h"
#include "cs_timer_stats.h"

#include "cs_field.h"
#include "cs_field_pointer.h"

#include "cs_lagr_clogging.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_dlvo.h"
#include "cs_lagr_stat.h"
#include "cs_lagr.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_precipitation_model.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * File variables
 *============================================================================*/

/*============================================================================
 * Public function definitions for Fortran API
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Mass source term due to precipitation
 *
 * \param[in]   dtref      pointer to reference time step
 * \param[in]   crom       pointer to density value
 * \param[in]   cvar_scal  pointer to scal  value
 * \param[out]  crvexp     explicit part of the source term
 */
/*----------------------------------------------------------------------------*/
void
CS_PROCF (precst,PRECST) (cs_real_t *dtref,
                          cs_real_t *crom,
                          cs_real_t *cvar_scal,
                          cs_real_t  crvexp[])
{
  cs_real_t pis6 = cs_math_pi / 6.0;

  cs_lagr_precipitation_model_t *preci = cs_get_lagr_precipitation_model();
  cs_real_t *mp_diss = preci->mp_diss;
  cs_real_t *solub = preci->solub;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;

  assert(cs_glob_lagr_model->precipitation == 1);

  /* ====================================================================
   * 1. Initialization
   * ====================================================================   */

  // TODO : move to a global initialization routine for alloc and free

  if (mp_diss == NULL)
    BFT_MALLOC(mp_diss, mesh->n_cells_with_ghosts * preci->nbrclas, cs_real_t);
  if (solub == NULL)
    BFT_MALLOC(solub, mesh->n_cells_with_ghosts, cs_real_t);

  cs_real_t *mp_preci;
  cs_lnum_t *part_tot;
  BFT_MALLOC(mp_preci, mesh->n_cells_with_ghosts, cs_real_t);
  BFT_MALLOC(part_tot, mesh->n_cells_with_ghosts, cs_lnum_t);

  /* reference diameter taken from first injection (FIXME) */

  cs_real_t ref_diameter = 0;
  {
    const cs_lagr_zone_data_t *bcs = cs_glob_lagr_boundary_conditions;
    for (int z_id = 0; z_id < bcs->n_zones; z_id++) {
      if (bcs->n_injection_sets[z_id] > 0) {
        ref_diameter = bcs->injection_set[z_id][0].diameter;
        break;
      }
    }
  }

  /* ====================================================================
   * 2. Calculation of the mass source terms due to
   *    precipitation and dissolution phenomena
   * ==================================================================== */

  if (preci->nbrclas > 0) {

    if (p_set->n_particles > 0) {

      for (cs_lnum_t iel = 0; iel < mesh->n_cells; iel++) {

        for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

          unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

          cs_real_t part_mass
            =   preci->rho * pis6
              * pow(cs_lagr_particle_get_real(particle, p_am,
                                              CS_LAGR_DIAMETER),3.0);

          if (   cs_lagr_particle_get_cell_id(particle, p_am) == iel
              &&   cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS)
                 - part_mass < 1e-12)

            /* number of magnetite particles in the cell iel */
            part_tot[iel] += 1;

        }

      }

    }

    /* Source term applied to second scalar */

    for (cs_lnum_t iel = 0; iel < mesh->n_cells; iel++) {

      preci->nbprec[iel]  = 0;

      /* PRECIPITATION   */
      if (cvar_scal[iel] >= solub[iel]) {

        cs_real_t mass = pis6 * pow(preci->diameter, 3) * preci->rho;
        preci->nbprec[iel] =   (cvar_scal[iel] - solub[iel])
                             * fvq->cell_vol[iel] / mass;
        mp_preci[iel] =  preci->nbprec[iel] * mass;
        crvexp[iel]   = -crom[iel] * mp_preci[iel] / *dtref;

      }

      /* to do:  impose a limit on  nbprec   */
      /* DISSOLUTION     */
      if (cvar_scal[iel] < solub[iel] && part_tot[iel] >= 1) {

        if (p_set->n_particles > 0) {

          for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

            unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

            for (cs_lnum_t iclas = 0; iclas < preci->nbrclas; iclas++) {

              cs_real_t p_diam = cs_lagr_particle_get_real(particle, p_am,
                                                           CS_LAGR_DIAMETER);
              cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am,
                                                           CS_LAGR_MASS);
              cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);
              cs_real_t mass = preci->rho * pis6 * pow(p_diam,3.0);
              if (   cell_id == iel
                  && p_diam - ref_diameter < 1e-12
                  && p_mass - mass < 1e-12) {

                cs_real_t p_weight = cs_lagr_particle_get_real(particle, p_am,
                                                               CS_LAGR_STAT_WEIGHT);

                if (   ((solub[iel] - cvar_scal[iel]) * fvq->cell_vol[iel])
                    >= (mp_diss[iel * preci->nbrclas + iclas] + p_weight * mass) )
                  mp_diss[iel * preci->nbrclas + iclas] += p_weight * mass;

              }

            }

          }

        }

        for (cs_lnum_t iclas = 0; iclas < preci->nbrclas; iclas++)
          crvexp[iel] += crom[iel] * mp_diss[iel * preci->nbrclas + iclas] / *dtref;

      }

    }

  }

  BFT_FREE(mp_preci);
  BFT_FREE(part_tot);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Management of the injection of particles formed by precipitation.
 *
 * \param[in]   vela       pointer to fluid velocity array (per cell)
 * \param[out]  val        number of particles to inject (with weight)
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_precipitation_injection(cs_real_t   *vela,
                                cs_real_t   *val)
{
  cs_real_t pis6 = cs_math_pi / 6.0;

  cs_lagr_precipitation_model_t *preci = cs_get_lagr_precipitation_model();
  cs_real_t *mp_diss = preci->mp_diss;

  cs_mesh_t  *mesh = cs_glob_mesh;
  cs_mesh_quantities_t *fvq = cs_glob_mesh_quantities;

  cs_lagr_particle_set_t  *p_set = cs_lagr_get_particle_set();
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;

  /* ============================================================ */
  /* 1. INITIALIZATION    */
  /* ============================================================ */

  /* number of dissolved particles  */
  cs_lnum_t *nbdiss;
  BFT_MALLOC(nbdiss, preci->nbrclas, cs_lnum_t);

  cs_real_t *mp;
  BFT_MALLOC(mp, preci->nbrclas, cs_real_t);

  /* mass of dissolved particles    */
  cs_real_t *mp_diss_t;
  BFT_MALLOC(mp_diss_t, mesh->n_cells_with_ghosts, cs_real_t);

  /* number of precipated particles */
  cs_lnum_t nbprec_tot = 0;
  cs_lnum_t nbprec2    = 0;

  /* reference diameter taken from first injection (FIXME) */

  cs_real_t ref_diameter = 0;
  {
    const cs_lagr_zone_data_t *bcs = cs_glob_lagr_boundary_conditions;
    for (int z_id = 0; z_id < bcs->n_zones; z_id++) {
      if (bcs->n_injection_sets[z_id] > 0) {
        ref_diameter = bcs->injection_set[z_id][0].diameter;
        break;
      }
    }
  }

  /* ============================================================ */
  /* 2. GESTION DES PARTICULES */
  /* ============================================================ */

  for (cs_lnum_t iel = 0; iel < mesh->n_cells; iel++) {
    nbprec2   = nbprec2 + preci->nbprec[iel];
  }
  if (nbprec2 >= 1000000.0) {

    /* Write(nfecra,1000) nbprec2 */

    cs_exit(1);
  }

  cs_lnum_t *cell;
  BFT_MALLOC(cell, nbprec2, cs_lnum_t);

  for (cs_lnum_t iel = 0; iel < mesh->n_cells; iel++) {

    /* Precipitation (Add particles)   */
    if (preci->nbprec[iel] > 0) {

      for (cs_lnum_t i = 0; i < preci->nbprec[iel]; i++)
        cell[nbprec_tot + i]      = iel;

      nbprec_tot += preci->nbprec[iel];

    }

    for (cs_lnum_t k = 0; k < preci->nbrclas; k++)
      mp_diss_t[iel] = mp_diss_t[iel] + mp_diss[iel * preci->nbrclas + k];

    /* Dissolution (Remove particles)  */

    if (mp_diss_t[iel] > 0) {

      mp[iel] = 0.0;

      for (cs_lnum_t npt = 0; npt < p_set->n_particles; npt++) {

        unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

        for (cs_lnum_t iclas = 0; iclas < preci->nbrclas; iclas++) {

          if (   cs_lagr_particle_get_cell_id(particle, p_am) == iel
              && (  cs_lagr_particle_get_real(particle, p_am, CS_LAGR_DIAMETER)
                  - ref_diameter < 1e-12)
              && (mp[iclas] < mp_diss[iel * preci->nbrclas + iclas])) {

            /* Removing of particles due to dissolution */

            cs_lagr_particle_set_cell_id(particle, p_am, -1);
            cs_real_t d3 = pow (cs_lagr_particle_get_real(particle, p_am,
                                                          CS_LAGR_DIAMETER), 3);
            mp[iclas] += cs_lagr_particle_get_real(particle, p_am,
                                                   CS_LAGR_STAT_WEIGHT)
              * (pis6 * d3 * preci->rho);
            nbdiss[iclas] += 1;

          }

        }

      }

    }

  }

  cs_lnum_t npt = p_set->n_particles;
  p_set->n_part_new += nbprec_tot;

  cs_lnum_t vv = cs_lagr_particle_set_resize(p_set->n_particles + p_set->n_part_new);
  assert(vv == 0);

  if (nbprec_tot >= 1) {

    for (cs_lnum_t ip = npt; ip < npt + nbprec_tot; ip++) {

      /* TODO: place particle at random location in the cell iel
         (not always at the cog) */

      unsigned char *particle = p_set->p_buffer + p_am->extents * npt;

      /* Random value associated with each particle */

      cs_real_t part_random = -1;
      cs_random_uniform(1, &part_random);
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_RANDOM_VALUE,
                                part_random);

      cs_real_t *part_coord = cs_lagr_particle_attr(particle, p_am, CS_LAGR_COORDS);

      for (cs_lnum_t i = 0; i <  3; i++)
        part_coord[i] = fvq->cell_cen[cell[ip - npt] * 3 + i];

      cs_lagr_particle_set_cell_id(particle, p_am, cell[ip - npt]);

      cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_REBOUND_ID, -1);

      cs_real_t *part_vel_seen = cs_lagr_particle_attr(particle, p_am,
                                                       CS_LAGR_VELOCITY_SEEN);
      for (cs_lnum_t i = 0; i < 3; i++)
        part_vel_seen[i] = vela[cell[ip - npt] * 3 + i];

      cs_real_t *part_vel = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY);
      for (cs_lnum_t i = 0; i < 3; i++)
        part_vel[i] = vela[cell[ip - npt] * 3 + i];

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, preci->diameter);

      cs_real_t mass =   pow(preci->diameter, 3.0) * preci->rho * pis6;
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS, mass);

      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_STAT_WEIGHT, 1.0);

      /* Residence time (may be negative to ensure continuous injection) */

      cs_real_t res_time = - part_random *cs_glob_lagr_time_step->dtp;
      cs_lagr_particle_set_real(particle, p_am, CS_LAGR_RESIDENCE_TIME,
                                res_time);

      if (cs_glob_lagr_model->deposition == 1) {
        cs_real_t random;
        cs_random_uniform(1, &random);
        cs_lagr_particle_set_real(particle, p_am,
                                  CS_LAGR_INTERF, 5.0 + 15.0 * random);
        cs_lagr_particle_set_real(particle, p_am,
                                  CS_LAGR_YPLUS, 1000.0);
        cs_lagr_particle_set_lnum(particle, p_am,
                                  CS_LAGR_MARKO_VALUE, -1);
        cs_lagr_particle_set_lnum(particle, p_am,
                                  CS_LAGR_NEIGHBOR_FACE_ID, -1);
        cs_lagr_particle_set_lnum(particle, p_am,
                                  CS_LAGR_DEPOSITION_FLAG, CS_LAGR_PART_IN_FLOW);

      }

    }

  }

  *val = 0.;

  for (cs_lnum_t ip = npt; ip < npt + nbprec_tot; ip++) {

    unsigned char *particle = p_set->p_buffer + p_am->extents * ip;
    *val += cs_lagr_particle_get_real(particle, p_am, CS_LAGR_STAT_WEIGHT);

  }

  p_set->n_particles += nbprec_tot;

  BFT_FREE(cell);
  BFT_FREE(nbdiss);
  BFT_FREE(mp);
  BFT_FREE(mp_diss_t);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
