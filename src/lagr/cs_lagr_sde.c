/*============================================================================
 * Methods for particle deposition
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
 * Functions dealing with particle deposition
 *============================================================================*/

#include "cs_defs.h"
#include "cs_math.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"

#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_physical_constants.h"
#include "cs_physical_model.h"
#include "cs_random.h"
#include "cs_thermal_model.h"

#include "cs_lagr.h"
#include "cs_lagr_adh.h"
#include "cs_lagr_deposition_model.h"
#include "cs_lagr_roughness.h"
#include "cs_lagr_tracking.h"
#include "cs_lagr_prototypes.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_lagr_sde.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Boltzmann constant */
static const double _k_boltz = 1.38e-23;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Private function prototypes (definitions follow)
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*! \brief Integration of SDEs by 1st order time scheme
 *
 * \param[in]  dtp       time step
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      lagrangian fluid characteristic time
 * \param[in]  piil      term in integration of UP SDEs
 * \param[in]  bx        caracteristiques de la turbulence
 * \param[in]  vagaus    gaussian random variables
 * \param[in]  brgaus    gaussian random variables
 * \param[in]  gradpr    pressure gradient
 * \param[in]  romp      particles associated density
 * \param[in]  force_p   taup times forces on particles (m/s)
 * \param[out] terbru
 */
/*------------------------------------------------------------------------------*/

static void
_lages1(cs_real_t           dtp,
        const cs_real_t     taup[],
        const cs_real_3_t   tlag[],
        const cs_real_3_t   piil[],
        const cs_real_33_t  bx[],
        const cs_real_33_t  vagaus[],
        const cs_real_t     brgaus[],
        const cs_real_3_t   gradpr[],
        const cs_real_t     romp[],
        const cs_real_3_t   force_p[],
        cs_real_t          *terbru)
{
  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  /* Initialisations*/

  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t vitf = 0.0;

  cs_real_t aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9, aux10, aux11;
  cs_real_t ter1f, ter2f, ter3f;
  cs_real_t ter1p, ter2p, ter3p, ter4p, ter5p;
  cs_real_t ter1x, ter2x, ter3x, ter4x, ter5x;
  cs_real_t p11, p21, p22, p31, p32, p33;
  cs_real_t omega2, gama2, omegam;
  cs_real_t grga2, gagam, gaome;
  cs_real_t tbrix1, tbrix2, tbriu;

  cs_lnum_t nor = cs_glob_lagr_time_step->nor;

  /* Integrate SDE's over particles */

  for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

    unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

    cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);

    if (cell_id >= 0) {

      cs_real_t *old_part_vel      = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                             CS_LAGR_VELOCITY);
      cs_real_t *old_part_vel_seen = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                             CS_LAGR_VELOCITY_SEEN);
      cs_real_t *part_vel          = cs_lagr_particle_attr(particle, p_am,
                                                           CS_LAGR_VELOCITY);
      cs_real_t *part_vel_seen     = cs_lagr_particle_attr(particle, p_am,
                                                           CS_LAGR_VELOCITY_SEEN);
      cs_real_t *part_coords       = cs_lagr_particle_attr(particle, p_am,
                                                           CS_LAGR_COORDS);
      cs_real_t *old_part_coords   = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                             CS_LAGR_COORDS);

      cs_real_t rom = extra->cromf->val[cell_id];

      for (cs_lnum_t id = 0; id < 3; id++) {

        vitf = extra->vel->vals[1][cell_id * 3 + id];

        /* --> (2.1) Calcul preliminaires :    */
        /* ----------------------------   */
        /* calcul de II*TL+<u> et [(grad<P>/rhop+g)*tau_p+<Uf>] ?  */

        cs_real_t tci = piil[ip][id] * tlag[ip][id] + vitf;
        cs_real_t force = force_p[ip][id];

        /* --> (2.2) Calcul des coefficients/termes deterministes */
        /* ----------------------------------------------------    */

        aux1 = exp(-dtp / taup[ip]);
        aux2 = exp(-dtp / tlag[ip][id]);
        aux3 = tlag[ip][id] / (tlag[ip][id] - taup[ip]);
        aux4 = tlag[ip][id] / (tlag[ip][id] + taup[ip]);
        aux5 = tlag[ip][id] * (1.0 - aux2);
        aux6 = cs_math_pow2(bx[ip][id][nor-1]) * tlag[ip][id];
        aux7 = tlag[ip][id] - taup[ip];
        aux8 = cs_math_pow2(bx[ip][id][nor-1]) * cs_math_pow2(aux3);

        /* --> trajectory terms */
        cs_real_t aa = taup[ip] * (1.0 - aux1);
        cs_real_t bb = (aux5 - aa) * aux3;
        cs_real_t cc = dtp - aa - bb;

        ter1x = aa * old_part_vel[id];
        ter2x = bb * old_part_vel_seen[id];
        ter3x = cc * tci;
        ter4x = (dtp - aa) * force;

        /* --> flow-seen velocity terms   */
        ter1f = old_part_vel_seen[id] * aux2;
        ter2f = tci * (1.0 - aux2);

        /* --> termes pour la vitesse des particules     */
        cs_real_t dd = aux3 * (aux2 - aux1);
        cs_real_t ee = 1.0 - aux1;

        ter1p = old_part_vel[id] * aux1;
        ter2p = old_part_vel_seen[id] * dd;
        ter3p = tci * (ee - dd);
        ter4p = force * ee;

        /* --> integrale sur la vitesse du fluide vu     */
        gama2  = 0.5 * (1.0 - aux2 * aux2);
        p11   = sqrt(gama2 * aux6);
        ter3f = p11 * vagaus[ip][id][0];

        /* --> integral for the particles velocity  */
        aux9  = 0.5 * tlag[ip][id] * (1.0 - aux2 * aux2);
        aux10 = 0.5 * taup[ip] * (1.0 - aux1 * aux1);
        aux11 =   taup[ip] * tlag[ip][id]
                * (1.0 - aux1 * aux2)
                / (taup[ip] + tlag[ip][id]);

        grga2 = (aux9 - 2.0 * aux11 + aux10) * aux8;
        gagam = (aux9 - aux11) * (aux8 / aux3);

        if (CS_ABS(p11) > cs_math_epzero) {

          p21 = gagam / p11;
          p22 = grga2 - cs_math_pow2(p21);
          p22 = sqrt(CS_MAX(0.0, p22));

        }
        else {

          p21 = 0.0;
          p22 = 0.0;

        }

        ter5p = p21 * vagaus[ip][id][0] + p22 * vagaus[ip][id][1];

        /* --> (2.3) Calcul des coefficients pour les integrales stochastiques :  */
        /* --> integrale sur la position des particules  */
        gaome = ( (tlag[ip][id] - taup[ip]) * (aux5 - aa)
                  - tlag[ip][id] * aux9
                  - taup[ip] * aux10
                  + (tlag[ip][id] + taup[ip]) * aux11)
                * aux8;
        omegam = aux3 * ( (tlag[ip][id] - taup[ip]) * (1.0 - aux2)
                          - 0.5 * tlag[ip][id] * (1.0 - aux2 * aux2)
                          + cs_math_pow2(taup[ip]) / (tlag[ip][id] + taup[ip]) * (1.0 - aux1 * aux2)
                          ) * aux6;
        omega2 =   aux7 * (aux7 * dtp - 2.0 * (tlag[ip][id] * aux5 - taup[ip] * aa))
                 + 0.5 * tlag[ip][id] * tlag [ip][id] * aux5 * (1.0 + aux2)
                 + 0.5 * taup[ip] * taup[ip] * aa * (1.0 + aux1)
                 - 2.0 * aux4 * tlag[ip][id] * taup[ip] * taup[ip] * (1.0 - aux1* aux2);
        omega2 = aux8 * omega2;

        if (p11 > cs_math_epzero)
          p31 = omegam / p11;
        else
          p31 = 0.0;

        if (p22 > cs_math_epzero)
          p32 = (gaome - p31 * p21) / p22;
        else
          p32 = 0.0;

        p33 = omega2 - cs_math_pow2(p31) - cs_math_pow2(p32);
        p33 = sqrt(CS_MAX(0.0, p33));
        ter5x = p31 * vagaus[ip][id][0] + p32 * vagaus[ip][id][1] + p33 * vagaus[ip][id][2];

        /* --> (2.3) Calcul des Termes dans le cas du mouvement Brownien :   */
        if (cs_glob_lagr_brownian->lamvbr == 1) {

          /* Calcul de la temperature du fluide en fonction du type  */
          /* d'ecoulement    */
          cs_real_t tempf;
          if (   cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0
              || cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] >= 0)
            tempf = extra->t_gaz->val[cell_id];

          else if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] >= 0
                   || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] >= 0
                   || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 0
                   || cs_glob_physical_model_flag[CS_JOULE_EFFECT] >= 0)
            tempf = extra->temperature->val[cell_id];

          else if (   cs_glob_thermal_model->itherm ==
                              CS_THERMAL_MODEL_TEMPERATURE
                   && cs_glob_thermal_model->itpscl ==
                              CS_TEMPERATURE_SCALE_CELSIUS)
            tempf = extra->scal_t->val[cell_id] + tkelvi;

          else if (   cs_glob_thermal_model->itherm ==
                              CS_THERMAL_MODEL_TEMPERATURE
                   && cs_glob_thermal_model->itpscl ==
                              CS_TEMPERATURE_SCALE_KELVIN)
            tempf = extra->scal_t->val[cell_id];

          else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

            cs_lnum_t mode  = 1;
            CS_PROCF (usthht,USTHHT) (&mode, &(extra->scal_t->val[cell_id]), &tempf);

            tempf = tempf + tkelvi;

          }

          else
            tempf = cs_glob_fluid_properties->t0;

          cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS);

          cs_real_t ddbr = sqrt(2.0 * _k_boltz * tempf / (p_mass * taup[ip]));

          cs_real_t tix2 = cs_math_pow2((taup[ip] * ddbr)) * (dtp - taup[ip] * (1.0 - aux1) * (3.0 - aux1) / 2.0);
          cs_real_t tiu2 = ddbr * ddbr * taup[ip] * (1.0 - exp( -2.0 * dtp / taup[ip])) / 2.0;

          cs_real_t tixiu  = cs_math_pow2((ddbr * taup[ip] * (1.0 - aux1))) / 2.0;

          tbrix2 = tix2 - (tixiu * tixiu) / tiu2;

          if (tbrix2 > 0.0)
            tbrix2    = sqrt(tbrix2) * brgaus[ip * 6 + id];
          else
            tbrix2    = 0.0;

          if (tiu2 > 0.0)
            tbrix1    = tixiu / sqrt(tiu2) * brgaus[ip * 6 + id + 3];
          else
            tbrix1    = 0.0;

          if (tiu2 > 0.0) {

            tbriu      = sqrt(tiu2) * brgaus[ip * 6 + id + 3];
            terbru[ip] = sqrt(tiu2);

          }
          else {

            tbriu     = 0.0;
            terbru[ip]  = 0.0;

          }

        }
        else {

          tbrix1  = 0.0;
          tbrix2  = 0.0;
          tbriu   = 0.0;

        }

        /* Finalisation des ecritures */

        /* --> trajectoire */
        part_coords[id] = old_part_coords[id] + ter1x + ter2x + ter3x
                                              + ter4x + ter5x + tbrix1 + tbrix2;

        /* --> vitesse fluide vu     */
        part_vel_seen[id] = ter1f + ter2f + ter3f;

        /* --> vitesse particules    */
        part_vel[id] = ter1p + ter2p + ter3p + ter4p + ter5p + tbriu;

      }


    }

  }
}

/*----------------------------------------------------------------------------*/
/*! \brief Integration of SDEs by 2nd order scheme
 *
 * When there has beed interaction with a boundary face, the velocity and
 * velocity seen computations are forced to 1st order.
 *
 * \param[in]  taup      temps caracteristique dynamique
 * \param[in]  tlag      temps caracteristique fluide
 * \param[in]  piil      terme dans l'integration des eds up
 * \param[in]  bx        caracteristiques de la turbulence
 * \param[in]  tsfext    infos pour couplage retour dynamique
 * \param[in]  vagaus    variables aleatoires gaussiennes
 * \param[in]  brgaus    gaussian variable for brownian movement
 * \param[in]  gradpr    pressure gradient
 * \param[in]  romp      masse volumique des particules
 * \param[in]  force_p   taup times forces on particles (m/s)
 * \param[out] terbru
 */
/*----------------------------------------------------------------------------*/

static void
_lages2(cs_real_t           dtp,
        const cs_real_t     taup[],
        const cs_real_3_t   tlag[],
        const cs_real_3_t   piil[],
        const cs_real_33_t  bx[],
        cs_real_t           tsfext[],
        const cs_real_33_t  vagaus[],
        const cs_real_t     brgaus[],
        const cs_real_3_t   gradpr[],
        const cs_real_t     romp[],
        const cs_real_3_t   force_p[],
        cs_real_t          *terbru)
{
  cs_real_t  aux0, aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9;
  cs_real_t  aux10, aux11, aux12, aux17, aux18, aux19;
  cs_real_t  aux20;
  cs_real_t  ter1, ter2, ter3, ter4, ter5;
  cs_real_t  sige, tapn, gamma2;
  cs_real_t  grgam2, gagam;
  cs_real_t  p11, p21, p22;
  cs_real_t  tbriu;

  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  /* ==============================================================================*/
  /* 1. Initialisations                                                            */
  /* ==============================================================================*/

  cs_lnum_t nor = cs_glob_lagr_time_step->nor;

  cs_real_t *auxl;
  BFT_MALLOC(auxl, p_set->n_particles*6, cs_real_t);

  /* =============================================================================
   * 2. Integration of the SDE on partilces
   * =============================================================================*/

  /* =============================================================================
   * 2.1 CALCUL A CHAQUE SOUS PAS DE TEMPS
   * ==============================================================================*/
  /* --> Compute tau_p*A_p and II*TL+<u> :
   *     -------------------------------------*/

  for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

    unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

    cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);

    if (cell_id >= 0) {

      for (cs_lnum_t id = 0; id < 3; id++) {

        auxl[ip * 6 + id] = force_p[ip][id];

        if (nor == 1)
          auxl[ip * 6 + id + 3] =   piil[ip][id] * tlag[ip][id]
                                  + extra->vel->vals[1][cell_id * 3 + id];
        else
          auxl[ip * 6 + id + 3] =   piil[ip][id] * tlag[ip][id]
                                  + extra->vel->vals[0][cell_id * 3 + id];

      }

    }

  }

  /* ==============================================================================*/
  /* 2.2 ETAPE DE PREDICTION : */
  /* ==============================================================================*/

  if (nor == 1) {

    /* --> Sauvegarde de tau_p^n */
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      if (cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_NUM) >= 0)
        cs_lagr_particle_set_real(particle, p_am, CS_LAGR_TAUP_AUX, taup[ip]);

    }

    /* --> Sauvegarde couplage   */
    if (cs_glob_lagr_time_scheme->iilagr == CS_LAGR_TWOWAY_COUPLING) {

      for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

        unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

        if (cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_CELL_NUM) >= 0) {

          aux0     = -dtp / taup[ip];
          aux1     =  exp(aux0);
          tsfext[ip] =   taup[ip]
                            * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS)
                            * (-aux1 + (aux1 - 1.0) / aux0);

        }

      }

    }

    /* --> Chargement des termes a t = t_n :    */
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);

      if (cell_id >= 0) {

        cs_real_t *old_part_vel      = cs_lagr_particle_attr_n(particle, p_am,
                                                               1, CS_LAGR_VELOCITY);
        cs_real_t *old_part_vel_seen = cs_lagr_particle_attr_n(particle, p_am,
                                                               1, CS_LAGR_VELOCITY_SEEN);
        cs_real_t *pred_part_vel_seen = cs_lagr_particle_attr(particle, p_am,
                                                              CS_LAGR_PRED_VELOCITY_SEEN);
        cs_real_t *pred_part_vel = cs_lagr_particle_attr(particle, p_am,
                                                         CS_LAGR_PRED_VELOCITY);

        for (cs_lnum_t id = 0; id < 3; id++) {

          aux0    =  -dtp / taup[ip];
          aux1    =  -dtp / tlag[ip][id];
          aux2    = exp(aux0);
          aux3    = exp(aux1);
          aux4    = tlag[ip][id] / (tlag[ip][id] - taup[ip]);
          aux5    = aux3 - aux2;

          pred_part_vel_seen[id] =   0.5 * old_part_vel_seen[id]
                                   * aux3 + auxl[ip * 6 + id + 3]
                                   * (-aux3 + (aux3 - 1.0) / aux1);

          ter1    = 0.5 * old_part_vel[id] * aux2;
          ter2    = 0.5 * old_part_vel_seen[id] * aux4 * aux5;
          ter3    = auxl[ip * 6 + id + 3] * (-aux2 + ((tlag[ip][id] + taup[ip]) / dtp) * (1.0 - aux2)
                                             - (1.0 + tlag[ip][id] / dtp) * aux4 * aux5);
          ter4    = auxl[ip * 6 + id] * (-aux2 + (aux2 - 1.0) / aux0);
          pred_part_vel[id] = ter1 + ter2 + ter3 + ter4;

        }

      }

    }

    /* Euler scheme */

    _lages1(dtp,
            taup,
            tlag,
            piil,
            bx,
            vagaus,
            brgaus,
            gradpr,
            romp,
            force_p,
            terbru);
  }

  else {

    /* Correction stage
       ---------------- */

    /* Compute Us */

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);

      if (    cell_id >= 0
           && cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_REBOUND_ID) != 0) {

        cs_real_t *part_vel
          = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY);
        cs_real_t *part_vel_seen
          = cs_lagr_particle_attr(particle, p_am, CS_LAGR_VELOCITY_SEEN);
        cs_real_t *old_part_vel
          = cs_lagr_particle_attr_n(particle, p_am, 1, CS_LAGR_VELOCITY);
        cs_real_t *old_part_vel_seen
          = cs_lagr_particle_attr_n(particle, p_am, 1, CS_LAGR_VELOCITY_SEEN);
        cs_real_t *pred_part_vel_seen
          = cs_lagr_particle_attr(particle, p_am, CS_LAGR_PRED_VELOCITY_SEEN);
        cs_real_t *pred_part_vel
          = cs_lagr_particle_attr(particle, p_am, CS_LAGR_PRED_VELOCITY);

        for (cs_lnum_t id = 0; id < 3; id++) {

          aux0    =  -dtp / taup[ip];
          aux1    =  -dtp / tlag[ip][id];
          aux2    = exp(aux0);
          aux3    = exp(aux1);
          aux4    = tlag[ip][id] / (tlag[ip][id] - taup[ip]);
          aux5    = aux3 - aux2;
          aux6    = aux3 * aux3;

          ter1    = 0.5 * old_part_vel_seen[id] * aux3;
          ter2    = auxl[ip * 6 + id + 3] * (1.0 - (aux3 - 1.0) / aux1);
          ter3    =  -aux6 + (aux6 - 1.0) / (2.0 * aux1);
          ter4    = 1.0 - (aux6 - 1.0) / (2.0 * aux1);

          sige    =   (  ter3 * bx[ip][id][0]
                       + ter4 * bx[ip][id][1] )
                    * (1.0 / (1.0 - aux6));

          ter5    = 0.5 * tlag[ip][id] * (1.0 - aux6);

          part_vel_seen[id] = pred_part_vel_seen[id] + ter1 + ter2 + sige * sqrt(ter5) * vagaus[ip][id][0];

          /* --> Calcul de Up :   */
          ter1    = 0.5 * old_part_vel[id] * aux2;
          ter2    = 0.5 * old_part_vel_seen[id] * aux4 * aux5;
          ter3    = auxl[ip * 6 + id + 3]
            * (1.0 - ((tlag[ip][id] + taup[ip]) / dtp) * (1.0 - aux2) + (tlag[ip][id] / dtp) * aux4 * aux5)
            + auxl[ip * 6 + id] * (1.0 - (aux2 - 1.0) / aux0);

          tapn    = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_TAUP_AUX);

          aux7    = exp(-dtp / tapn);
          aux8    = 1.0 - aux3 * aux7;
          aux9    = 1.0 - aux6;
          aux10   = 1.0 - aux7 * aux7;
          aux11   = tapn / (tlag[ip][id] + tapn);
          aux12   = tlag[ip][id] / (tlag[ip][id] - tapn);
          aux17   = sige * sige * aux12 * aux12;
          aux18   = 0.5 * tlag[ip][id] * aux9;
          aux19   = 0.5 * tapn * aux10;
          aux20   = tlag[ip][id] * aux11 * aux8;

          /* ---> calcul de la matrice de correlation */
          gamma2  = sige * sige * aux18;
          grgam2  = aux17 * (aux18 - 2.0 * aux20 + aux19);
          gagam   = sige * sige * aux12 * (aux18 - aux20);

          /* ---> simulation du vecteur Gaussien */

          p11     = sqrt(CS_MAX(0.0, gamma2));
          if (p11 > cs_math_epzero) {

            p21  = gagam / p11;
            p22  = grgam2 - p21 * p21;
            p22  = sqrt(CS_MAX(0.0, p22));

          }
          else {

            p21  = 0.0;
            p22  = 0.0;

          }

          ter4    = p21 * vagaus[ip][id][0] + p22 * vagaus[ip][id][1];

          /* ---> Calcul des Termes dans le cas du mouvement Brownien     */
          if (cs_glob_lagr_brownian->lamvbr == 1)
            tbriu = terbru[ip] * brgaus[ip * 6 + id + 3];
          else
            tbriu = 0.0;

          /* ---> finalisation de l'ecriture     */

          part_vel[id] = pred_part_vel[id] + ter1 + ter2 + ter3 + ter4 + tbriu;

        }

      }

    }

  }

  BFT_FREE(auxl);
}

/*----------------------------------------------------------------------------*/
/*! \brief Deposition submodel
 *
 *  1/ Modification of the coordinate system (global ->local)
 *  2/ Call of subroutine lagcli
 *  3/ Integration of the stochastic differential equations
 *     in the 2 directions different from the normal to the boundary face
 *  4/ Modification of the coordinate system (local ->global)
 *  5/ Update of the particle position
 *
 * \param[in]  dtp       time step
 * \param[in]  ip        particle id
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      fluid characteristic time
 * \param[in]  piil      term in integration of UP SDEs
 * \param[in]  vagaus    gaussian random variables
 * \param[in]  gradpr    pressure gradient
 * \param[in]  romp      particles associated density
 * \param[in]  force_p   taup times forces on particles (m/s)
 * \param[in]  tempf     temperature of the fluid (K)
 * \param[in]  vislen    FIXME
 * \param[in]  depint    interface location near-wall/core-flow
 */
/*----------------------------------------------------------------------------*/

static void
_lagesd(cs_real_t           dtp,
        cs_lnum_t           ip,
        const cs_real_t     taup[],
        const cs_real_3_t   piil[],
        const cs_real_33_t  vagaus[],
        const cs_real_3_t   gradpr[],
        const cs_real_t     romp[],
        const cs_real_3_t   force_p[],
        cs_real_t           tempf,
        const cs_real_t     vislen[],
        cs_real_t          *depint,
        cs_lnum_t          *nresnew)
{
  /* mesh and mesh quantities */
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;

  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  /* Hydrodynamic drag and torque on a deposited particle     */
  cs_real_t drag_force[3];
  cs_real_t drag_torque[3];

  /* Lift force and torque on a deposited particle   */
  cs_real_t lift_force[1];
  cs_real_t lift_torque[3];

  /* Gravity force and torque on a deposited particle   */
  cs_real_t grav_force[3];
  cs_real_t grav_torque[3];

  /* Adhesion force and torque on a deposited particle   */
  cs_real_t adhes_torque[3];

  /* Map field arrays     */
  cs_real_t *vela = extra->vel->vals[1];

  /* ===========================================================================
   * 1. INITIALISATIONS
   * ======================================================================== */

  const cs_real_t  *grav  = cs_glob_physical_constants->gravity;

  /* particle data */
  unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

  cs_real_t p_mass = cs_lagr_particle_get_real(particle, p_am,
                                               CS_LAGR_MASS);
  cs_real_t p_diam = cs_lagr_particle_get_real(particle, p_am,
                                               CS_LAGR_DIAMETER);
  cs_real_t p_stat_w = cs_lagr_particle_get_real(particle, p_am,
                                                 CS_LAGR_STAT_WEIGHT);

  cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);

  cs_lnum_t face_id = cs_lagr_particle_get_lnum(particle, p_am,
                                                CS_LAGR_NEIGHBOR_FACE_ID);

  cs_real_t ustar = extra->ustar->val[face_id];
  cs_real_t lvisq = vislen[face_id];

  cs_real_t tvisq;
  if (ustar > 0.0)
    tvisq = lvisq / ustar;
  else
    tvisq = cs_math_big_r;

  /* Friction velocity    */

  face_id  = cs_lagr_particle_get_lnum(particle, p_am,
                                       CS_LAGR_NEIGHBOR_FACE_ID);

  /* Constants for the calculation of bxp and tlp  */
  cs_real_t c0 = 2.1;
  cs_real_t cl = 1.0 / (0.5 + (3.0 / 4.0) * c0);

  /* Pointer on the density w.r.t the flow    */
  cs_real_t romf = extra->cromf->val[cell_id];

  cs_real_t visccf = extra->viscl->val[cell_id] / romf;

  cs_real_t yplus = cs_lagr_particle_get_real(particle, p_am,
                                              CS_LAGR_YPLUS);

  /* Turbulent kinetic energy and dissipation w.r.t y+  */
  cs_real_t energi, dissip;
  if (yplus <= 5.0) {

    energi = 0.1 * cs_math_pow2(yplus) * cs_math_pow2(ustar);
    dissip = 0.2 * cs_math_pow4(ustar) / visccf;

  }

  else if (yplus <= 30.0) {

    energi = cs_math_pow2(ustar) / sqrt(0.09);
    dissip = 0.2 * cs_math_pow4(ustar) / visccf;

  }
  else {

    assert(yplus <= 100.0); /* should not arrive here otherwise */

    energi   = cs_math_pow2(ustar) / sqrt(0.09);
    dissip   = cs_math_pow4(ustar) / (0.41 * yplus * visccf);

  }

  /* ===========================================================================
   * 2. Reference frame change:
   * --------------------------
   * global reference frame --> local reference frame for the boundary face
   * ======================================================================== */

  const cs_real_3_t *rot_m
    = (const cs_real_3_t *)cs_glob_lagr_b_face_proj[face_id];

  /* 2.1 - particle velocity   */

  cs_real_t *old_part_vel = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                    CS_LAGR_VELOCITY);
  cs_real_t vpart[3];

  cs_math_33_3_product(rot_m, old_part_vel, vpart);

  cs_real_t vpart0[3] = {vpart[0], vpart[1], vpart[2]};

  /* 2.2 - flow-seen velocity  */

  cs_real_t *old_part_vel_seen = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                         CS_LAGR_VELOCITY_SEEN);
  cs_real_t vvue[3];

  cs_math_33_3_product(rot_m, old_part_vel_seen, vvue);

  /* 2.3 - Gravity vector */

  cs_real_t ggp[3];

  cs_math_33_3_product(rot_m, grav, ggp);

  /* 2.4 - flow velocity  */

  cs_real_t vflui[3];

  cs_math_33_3_product(rot_m, vela, vflui);

  cs_real_t norm = sqrt(cs_math_pow2(vflui[1]) + cs_math_pow2(vflui[2]));

  /* Velocity norm w.r.t y+ */
  cs_real_t norm_vit;

  if (yplus <= 5.0)
    norm_vit = yplus * ustar;

  else if (yplus <= 30.0)
    norm_vit = ( -3.05 + 5.0 * log(yplus)) * ustar;

  else {
    assert(yplus < 100.0);
    norm_vit = (2.5 * log (yplus) + 5.5) * ustar;
  }

  if (norm_vit > 0.) {
    vflui[1] = norm_vit * vflui[1] / norm;
    vflui[2] = norm_vit * vflui[2] / norm;
  }
  else {
    vflui[1] = 0.;
    vflui[2] = 0.;
  }

  /* 2.5 Particle force: - pressure gradient/romp + external force + g   */

  cs_real_t force_pn[3];

  cs_math_33_3_product(rot_m, force_p[ip], force_pn);

  /* 2.6 - "piil" term    */

  cs_real_t piilp[3];

  cs_math_33_3_product(rot_m, piil[ip], piilp);

  /* 2.7 - tlag */

  cs_real_t tlp = cs_math_epzero;
  if (energi > 0.0) {

    tlp = cl * energi / dissip;
    tlp = CS_MAX(tlp, cs_math_epzero);

  }

  /* 2.8 - bx   */
  cs_real_t bxp = sqrt(c0 * dissip);

  /* =========================================================================
   * 3. Integration of the EDS on the particles
   * =========================================================================*/

  /* Retrieve of the turbulent kinetic energy */
  cs_real_t  enertur;
  if (extra->itytur == 2 || extra->iturb == 50 || extra->iturb == 60)
    enertur  = extra->cvar_k->vals[1][cell_id];

  else if (extra->itytur == 3) {
    if (extra->cvar_rij == NULL) {
      enertur  = 0.5 * (  extra->cvar_r11->vals[1][cell_id]
                        + extra->cvar_r22->vals[1][cell_id]
                        + extra->cvar_r33->vals[1][cell_id]);
    } else {
      enertur  = 0.5 * (  extra->cvar_rij->vals[1][6*cell_id    ]
                        + extra->cvar_rij->vals[1][6*cell_id + 1]
                        + extra->cvar_rij->vals[1][6*cell_id + 2]);
    }
  }

  cs_lnum_t marko  = cs_lagr_particle_get_lnum(particle, p_am,
                                               CS_LAGR_MARKO_VALUE);
  cs_real_t interf = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_INTERF);
  cs_real_3_t depl;

  cs_lagr_deposition(dtp,
                     &marko,
                     tempf,
                     lvisq,
                     tvisq,
                     vpart,
                     vvue,
                     depl,
                     &p_diam,
                     romp[ip],
                     taup[ip],
                     &yplus,
                     &interf,
                     &enertur,
                     ggp,
                     vflui,
                     force_pn,
                     piilp,
                     depint);

  cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_MARKO_VALUE, marko);

  if (   cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
      != CS_LAGR_PART_IN_FLOW) {
    depl[0]  = 0.0;
    vpart[0] = 0.0;
  }

  /* Integration in the 2 other directions */
  if (   cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
      == CS_LAGR_PART_IN_FLOW) {

    for (cs_lnum_t id = 1; id < 3; id++) {

      cs_lnum_t i0 = id-1;//FIXME strange

      cs_real_t tci   = piilp[id] * tlp + vflui[id];
      cs_real_t force = force_pn[id];
      cs_real_t aux1  = exp(-dtp / taup[ip]);
      cs_real_t aux2  = exp(-dtp / tlp);
      cs_real_t aux3  = tlp / (tlp - taup[ip]);
      cs_real_t aux4  = tlp / (tlp + taup[ip]);
      cs_real_t aux5  = tlp * (1.0 - aux2);
      cs_real_t aux6  = bxp * bxp * tlp;
      cs_real_t aux7  = tlp - taup[ip];
      cs_real_t aux8  = bxp * bxp * cs_math_pow2(aux3);

      /* --> Terms for the trajectory   */
      cs_real_t aa    = taup[ip] * (1.0 - aux1);
      cs_real_t bb    = (aux5 - aa) * aux3;
      cs_real_t cc    = dtp - aa - bb;
      cs_real_t ter1x = aa * vpart[id];
      cs_real_t ter2x = bb * vvue[id];
      cs_real_t ter3x = cc * tci;
      cs_real_t ter4x = (dtp - aa) * force;

      /* --> Terms for the flow-seen velocity     */
      cs_real_t ter1f = vvue[id] * aux2;
      cs_real_t ter2f = tci * (1.0 - aux2);

      /* --> Terms for the particles velocity     */
      cs_real_t dd    = aux3 * (aux2 - aux1);
      cs_real_t ee    = 1.0 - aux1;
      cs_real_t ter1p = vpart[id] * aux1;
      cs_real_t ter2p = vvue[id] * dd;
      cs_real_t ter3p = tci * (ee - dd);
      cs_real_t ter4p = force * ee;

      /* --> (2.3) Coefficients computation for the stochastic integrals:  */
      cs_real_t gama2  = 0.5 * (1.0 - aux2 * aux2);
      cs_real_t omegam = aux3 * ( (tlp - taup[ip]) * (1.0 - aux2)
                                  - 0.5 * tlp * (1.0 - aux2 * aux2)
                                  + cs_math_pow2(taup[ip]) / (tlp + taup[ip])
                                  * (1.0 - aux1 * aux2)
                                  ) * aux6;
      cs_real_t omega2 =   aux7
                         * (aux7 * dtp - 2.0 * (tlp * aux5 - taup[ip] * aa))
                         + 0.5 * tlp * tlp * aux5 * (1.0 + aux2)
                         + 0.5 * cs_math_pow2(taup[ip]) * aa * (1.0 + aux1)
                         - 2.0 * aux4 * tlp * cs_math_pow2(taup[ip]) * (1.0 - aux1 * aux2);
      omega2 *= aux8 ;

      cs_real_t  p11, p21, p22, p31, p32, p33;

      if (CS_ABS(gama2) >cs_math_epzero) {

        p21    = omegam / sqrt (gama2);
        p22    = omega2 - cs_math_pow2(p21);
        p22    = sqrt(CS_MAX(0.0, p22));

      }
      else {

        p21    = 0.0;
        p22    = 0.0;

      }

      cs_real_t ter5x = p21 * vagaus[ip][i0][0] + p22 * vagaus[ip][i0][1];

      /* --> Integral for the flow-seen velocity  */

      p11   = sqrt(gama2 * aux6);
      cs_real_t ter3f = p11 * vagaus[ip][i0][0];

      /* --> Integral for particles velocity */

      cs_real_t aux9  = 0.5 * tlp * (1.0 - aux2 * aux2);
      cs_real_t aux10 = 0.5 * taup[ip] * (1.0 - aux1 * aux1);
      cs_real_t aux11 = taup[ip] * tlp * (1.0 - aux1 * aux2) / (taup[ip] + tlp);
      cs_real_t grga2 = (aux9 - 2.0 * aux11 + aux10) * aux8;
      cs_real_t gagam = (aux9 - aux11) * (aux8 / aux3);
      cs_real_t gaome = (  (tlp - taup[ip]) * (aux5 - aa)
                         - tlp * aux9
                         - taup[ip] * aux10
                         + (tlp + taup[ip]) * aux11) * aux8;

      if (p11 > cs_math_epzero)
        p31 = gagam / p11;

      else
        p31 = 0.0;

      if (p22 > cs_math_epzero)
        p32 = (gaome - p31 * p21) / p22;

      else
        p32 = 0.0;

      p33 = grga2 - cs_math_pow2(p31) - cs_math_pow2(p32);
      p33 = sqrt (CS_MAX(0.0, p33));

      cs_real_t ter5p =   p31 * vagaus[ip][i0][0]
                        + p32 * vagaus[ip][i0][1]
                        + p33 * vagaus[ip][i0][2];

      /* --> trajectory  */
      depl[id] = ter1x + ter2x + ter3x + ter4x + ter5x;

      /* --> flow-seen velocity    */
      vvue[id] = ter1f + ter2f + ter3f;

      /* --> particles velocity    */
      vpart[id] = ter1p + ter2p + ter3p + ter4p + ter5p;

    }

  }
  else {

    for (cs_lnum_t id = 1; id < 3; id++) {

      cs_lnum_t i0 = id-1;

      cs_real_t tci   = piilp[id] * tlp + vflui[id];
      cs_real_t aux2  = exp(-dtp / tlp);
      cs_real_t aux6  = bxp * bxp * tlp;

      /* --> Terms for the flow-seen velocity     */
      cs_real_t ter1f = vvue[id] * aux2;
      cs_real_t ter2f = tci * (1.0 - aux2);

      /* --> (2.3) Coefficients computation for the stochastic integrals:  */
      cs_real_t gama2 = 0.5 * (1.0 - aux2 * aux2);

      /* --> Integral for the flow-seen velocity  */
      cs_real_t p11   = sqrt (gama2 * aux6);
      cs_real_t ter3f = p11 * vagaus[ip][i0][0];

      /* --> flow-seen velocity    */
      vvue[id] = ter1f + ter2f + ter3f;

    }

  }

  if (cs_glob_lagr_model->resuspension == 1) {

    cs_real_t p_height = cs_lagr_particle_get_real(particle, p_am,
                                                   CS_LAGR_HEIGHT);

    cs_lnum_t iresusp = 0;

    if (cs_lagr_particle_get_lnum(particle, p_am,
                                  CS_LAGR_DEPOSITION_FLAG) != CS_LAGR_PART_IN_FLOW) {

      cs_lnum_t n_f_id
        = cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_NEIGHBOR_FACE_ID);

      cs_lnum_t nfabor = cs_glob_mesh->n_b_faces;

      /* Resuspension model
       * differentiation between monolayer and multilayer resuspension
       * based on the mean deposit height
       * (another possibility is to use the jamming limit...) */
      cs_real_t diam_mean = cs_glob_lagr_clogging_model->diam_mean;

      if (   (cs_glob_lagr_model->clogging == 0 && iresusp == 0)
          || (   cs_glob_lagr_model->clogging == 1
              &&   bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepm]
                 < diam_mean && iresusp == 0)) {

        /* Monolayer resuspension */

        /* Calculation of the hydrodynamic drag and torque
         * applied on the deposited particle   */

        drag_force[0] =  3.0 * cs_math_pi * p_diam
          * (vvue[0] - vpart[0]) * visccf * romf * 3.39;
        drag_torque[0] = 0.0;

        for (cs_lnum_t id = 1; id < 3; id++) {

          drag_force[id]   =  3.0 * cs_math_pi * p_diam
            * (vvue[id] - vpart[id]) * visccf * romf * 1.7;
          drag_torque[id] = 1.4 * drag_force[id] * p_diam * 0.5;

        }

        /* Calculation of lift force and torque */
        lift_force[0] = - 20.0 * cs_math_pow2(visccf) * romf *
          pow(ustar * p_diam * 0.5 / visccf, 2.31);

        for (cs_lnum_t id = 1; id < 3; id++) {

          lift_torque[id] = -lift_force[0] * p_diam * 0.5
            * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

        }

        /* Calculation of gravity force and torque */
        for (cs_lnum_t id = 0; id < 3; id++) {
          grav_force[id] = 4./3. * cs_math_pi * cs_math_pow3(p_diam/2) * romp[ip] * ggp[id];
        }

        for (cs_lnum_t id = 1; id < 3; id++) {

          cs_real_t sign11, sign22;
          if (grav_force[1]*vvue[1] < 0 )
            sign11 = -1.0;
          else
            sign11 = 1.0;
          if (grav_force[2]*vvue[2] < 0 )
            sign22 = -1.0;
          else
            sign22 = 1.0;

          grav_torque[id] = (-grav_force[0] +
                              grav_force[1]*sign11 +
                              grav_force[2]*sign22 ) * p_diam * 0.5
            * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]));

        }

        /* Case with consolidation */
        if (    cs_glob_lagr_consolidation_model->iconsol > 0
             &&   cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT)
                > 0.01 * diam_mean) {

          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE,
                                    cs_glob_lagr_consolidation_model->force_consol);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE,
                                    cs_glob_lagr_consolidation_model->force_consol
                                    * p_diam * 0.5);

        }

        cs_real_t adhes_force = cs_lagr_particle_get_real(particle, p_am,
                                                          CS_LAGR_ADHESION_FORCE);

        /* Is there direct wall-normal lift-off of the particle ? */
        if (   (adhes_force + grav_force[0] + lift_force[0] + drag_force[0]) < 0
            && iresusp == 0 ) {

          /* Update of the number and weight of resuspended particles     */
          p_set->n_part_resusp += 1;
          p_set->weight_resusp += p_stat_w;

          if (cs_glob_lagr_boundary_interactions->iflmbd == 1) {

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ires]
              += p_stat_w;

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->iflres]
              += p_stat_w + (p_stat_w * p_mass / mq->b_f_face_surf[n_f_id]);

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->iflm]
              += - (p_stat_w * p_mass / mq->b_f_face_surf[n_f_id]);

          }

          /* Update of surface covered and deposit height
           * (for non-rolling particles) */
          if (   cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
              == CS_LAGR_PART_DEPOSITED ) {

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->iscovc]
              -=   cs_math_pi * cs_math_pow2(p_diam) * p_stat_w
                 * 0.25 / mq->b_f_face_surf[n_f_id];

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepm]
              -=   cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                 * 0.25 / mq->b_f_face_surf[n_f_id];

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepv]
              -=   cs_math_pow2(cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                 * 0.25 / mq->b_f_face_surf[n_f_id]);

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->inclg]
              -= p_stat_w;

          }

          /* The particle is resuspended  */
          vpart[0] = CS_MIN(-1.0 / p_mass * dtp
                            * CS_ABS(drag_force[0] - adhes_force),
                            0.001);
          if (   cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
              == CS_LAGR_PART_DEPOSITED ) {
            vpart[1] = 0.0;
            vpart[2] = 0.0;
          }

          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG, 0);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);

          if (cs_glob_lagr_model->clogging == 1) {
            if (CS_ABS(p_height-p_diam)/p_diam > 1.0e-6) {
              cs_real_t d_resusp = pow(0.75 * cs_math_pow2(p_diam) * p_height, 1.0/3.0);
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, d_resusp);
              cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, d_resusp);
            }
          }

          if (p_am->count[0][CS_LAGR_N_LARGE_ASPERITIES] > 0)
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_N_LARGE_ASPERITIES, 0);

          if (p_am->count[0][CS_LAGR_N_SMALL_ASPERITIES] > 0)
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_N_SMALL_ASPERITIES, 0);

          if (p_am->count[0][CS_LAGR_DISPLACEMENT_NORM] > 0)
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_DISPLACEMENT_NORM, 0.0);

          iresusp = 1;
        }
        /* No direct normal lift-off */
        else if (iresusp == 0) {

          /* Calculation of the norm of the hydrodynamic
           * torque and drag (tangential) */

          for (cs_lnum_t id = 1; id < 3; id++) {

            adhes_torque[id]  = - cs_lagr_particle_get_real(particle, p_am,
                                                            CS_LAGR_ADHESION_TORQUE)
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) ) ;
          }

          cs_real_t iner_tor = (7.0 / 5.0) * p_mass * cs_math_pow2((p_diam * 0.5));
          cs_real_t cst_1, cst_4;
          cst_4 =   6 * cs_math_pi * visccf
            * romf * 1.7 * 1.4
            * cs_math_pow2(p_diam * 0.5);
          cst_1 = cst_4 * (p_diam * 0.5) / iner_tor;

          for (cs_lnum_t id = 1; id < 3; id++) {

            vpart0[id] = vpart[id];
            vpart[id]  =  vpart0[id] * exp(-cst_1 * dtp)
                        + (vvue[id] + (adhes_torque[id] + lift_torque[id]
                                       + grav_torque[id]) / cst_4)
                        * (1.0 - exp(-cst_1 * dtp) );

          }

          cs_real_t scalax  = vpart[1] * vvue[1] + vpart[2] * vvue[2];

          if (scalax > 0.0) {

            /* The calculated particle velocity is aligned   */
            /* with the flow seen   */
            /* --> The particle starts or keep on rolling    */

            /* If the particle stars rolling:
             * update of the surface covered and deposit height */

            if (cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
                == CS_LAGR_PART_DEPOSITED ) {

              bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->iscovc]
                -=   cs_math_pi * cs_math_pow2(p_diam) * p_stat_w
                   * 0.25 / mq->b_f_face_surf[n_f_id];

              bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepm]
                -=   cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                   * 0.25 / mq->b_f_face_surf[n_f_id];

              bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepv]
                -=   cs_math_pow2(cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                   * 0.25 / mq->b_f_face_surf[n_f_id]);

              bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->inclg]
                -= p_stat_w;

            }

            cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG,
                                      CS_LAGR_PART_ROLLING);

            vpart[0] = 0.0;

            for (cs_lnum_t id = 1; id < 3; id++) {

              if (CS_ABS (vpart[id]) > CS_ABS (vvue[id]))
                /* The velocity of the rolling particle cannot   */
                /* exceed the surrounding fluid velocity    */
                vpart[id] = vvue[id];

              cs_real_t kkk = vvue[id] + adhes_torque[id] / cst_4;
              cs_real_t kk  = vpart0[id] - kkk;

              depl[id] = (kkk * dtp + kk / cst_1 * (1. - exp(-cst_1 * dtp)));

            }

            iresusp = 1;

          }
          /* if (scalax..) */
          else {

            /* The particle is not set into motion or stops
             * the flag is set to 1 and velocity and displacement are null */

            /* Simplified treatment:
             * no update of iscovc, ihdepm for particles that stop */
            if (   cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
                == CS_LAGR_PART_ROLLING)
              iresusp = 1;

            for (cs_lnum_t id = 1; id < 3; id++) {

              depl[id]     = 0.0;
              vpart[id]    = 0.0;

            }

          } /* if (scalax..)   */

        }

      } /* End of monolayer resuspension*/

      else if (   cs_glob_lagr_model->clogging == 1
               &&    bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepm]
                  >= diam_mean
               && iresusp == 0) {

        /* Multilayer resuspension model */

        cs_real_t mean_depo_height
          = bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepm];
        /* Calculation of the hydrodynamic forces and torques
         * applied on the deposited cluster :
         * Principle: drag>0 for protruding clusters, 0 otherwise */

        if (   p_height > mean_depo_height
            &&    cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
               == CS_LAGR_PART_DEPOSITED) {

          /* Calculation of drag force and torque*/
          drag_force[0] =   3.0 * cs_math_pi * (p_height - mean_depo_height)
                          * (vvue[0] - vpart[0]) * visccf * romf * 3.39;
          drag_torque[0] = 0.0;

          for (cs_lnum_t id = 1; id < 3; id++) {

            drag_force[id]  =   3.0 * cs_math_pi * (p_height - mean_depo_height)
                              * (vvue[id] - vpart[id]) * visccf * romf * 1.7;
            drag_torque[id] = 1.4 * drag_force[id] * p_diam * 0.5;

          }

          /* Calculation of lift force and torque */
          lift_force[0] = - 20.0 * cs_math_pow2(visccf) * romf *
            pow(ustar * (p_height - mean_depo_height) * 0.5 / visccf, 2.31);

          for (cs_lnum_t id = 1; id < 3; id++) {

            lift_torque[id] = -lift_force[0] * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }

          /* Calculation of gravity force and torque */
          for (cs_lnum_t id = 0; id < 3; id++) {
            grav_force[id] = p_diam * ggp[id]
              *(p_height - mean_depo_height) / p_height;
          }

          for (cs_lnum_t id = 1; id < 3; id++) {

            cs_real_t sign11, sign22;
            if (grav_force[1]*vvue[1] < 0 )
              sign11 = -1.0;
            else
              sign11 = 1.0;
            if (grav_force[2]*vvue[2] < 0 )
              sign22 = -1.0;
            else
              sign22 = 1.0;

            grav_torque[id] = (-grav_force[0] +
                                grav_force[1]*sign11 +
                                grav_force[2]*sign22 ) * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }

        }
        else if (   p_height <= mean_depo_height
                 &&    cs_lagr_particle_get_lnum(particle, p_am,
                                                 CS_LAGR_DEPOSITION_FLAG)
                    == CS_LAGR_PART_DEPOSITED ) {

          /* Calculation of drag force and torque*/
          for (cs_lnum_t id = 0; id < 3; id++) {
            drag_force[id]  = 0.0;
            drag_torque[id] = 0.0;
          }

          /* Calculation of lift force and torque */
          lift_force[0] = 0.0;

          for (cs_lnum_t id = 1; id < 3; id++) {
            lift_torque[id] = 0.0;
          }

          /* Calculation of gravity force and torque */
          for (cs_lnum_t id = 0; id < 3; id++) {
            grav_force[id] = p_diam * ggp[id];
          }

          for (cs_lnum_t id = 1; id < 3; id++) {

            cs_real_t sign11, sign22;
            if (grav_force[1]*vvue[1] < 0 )
              sign11 = -1.0;
            else
              sign11 = 1.0;
            if (grav_force[2]*vvue[2] < 0 )
              sign22 = -1.0;
            else
              sign22 = 1.0;

            grav_torque[id] = (-grav_force[0] +
                                grav_force[1]*sign11 +
                                grav_force[2]*sign22 ) * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }

        }
        else if (   cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
                 == CS_LAGR_PART_ROLLING ) {

          /* Calculation of drag force and torque*/
          drag_force[0] =  3.0 * cs_math_pi * p_diam
            * (vvue[0] - vpart[0]) * visccf * romf * 3.39;
          drag_torque[0] = 0.0;

          for (cs_lnum_t id = 1; id < 3; id++) {

            drag_force[id]   =  3.0 * cs_math_pi * p_diam
              * (vvue[id] - vpart[id]) * visccf * romf * 1.7;
            drag_torque[id] = 1.4 * drag_force[id] * p_diam * 0.5;

          }

          /* Calculation of lift force and torque */
          lift_force[0] = - 20.0 * cs_math_pow2(visccf) * romf *
            pow(ustar * p_diam * 0.5 / visccf, 2.31);

          for (cs_lnum_t id = 1; id < 3; id++) {

            lift_torque[id] = -lift_force[0] * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }

          /* Calculation of gravity force and torque */
          for (cs_lnum_t id = 0; id < 3; id++) {
            grav_force[id] = p_diam * ggp[id];
          }

          for (cs_lnum_t id = 1; id < 3; id++) {

            cs_real_t sign11, sign22;
            if (grav_force[1]*vvue[1] < 0 )
              sign11 = -1.0;
            else
              sign11 = 1.0;
            if (grav_force[2]*vvue[2] < 0 )
              sign22 = -1.0;
            else
              sign22 = 1.0;

            grav_torque[id] = (-grav_force[0] +
                                grav_force[1]*sign11 +
                                grav_force[2]*sign22 ) * p_diam * 0.5
              * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );

          }
        }

        /* Adhesion force and torque for multilayered structures:
         * equal to adhesion force between single particles
         * times the number of particle-particle contacts  */

        cs_real_t adhes_energ, adhes_force, adhes_force_ps;
        cs_lagr_adh_pp(p_diam, tempf, &adhes_energ, &adhes_force);
        /* Average number of contact in a cluster */
        cs_real_t ncont_pp = cs_math_pow2(p_diam/diam_mean);

        /* Case with consolidation */
        if (cs_glob_lagr_consolidation_model->iconsol > 0 &&
            cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG) == 1) {

          const cs_real_t consol_height
            = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT);

          if (consol_height > 0.01 * diam_mean) {
            adhes_force = cs_glob_lagr_consolidation_model->force_consol +
              (adhes_force - cs_glob_lagr_consolidation_model->force_consol) * 0.5
            * (1.0 + tanh((mean_depo_height - consol_height)
                          /(0.1 * consol_height)));
            adhes_force_ps = cs_glob_lagr_consolidation_model->force_consol;
          }
          else {
            adhes_force *= ncont_pp;
            adhes_force_ps = adhes_force;
          }
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_N_SMALL_ASPERITIES,
                                    ncont_pp);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE,
                                    adhes_force);
          cs_real_t adhes_tor = adhes_force * p_diam * 0.5;
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE,
                                    adhes_tor);
        }
        else {
          /* Case without consolidation */

          cs_lnum_t ncont = 1;

          if (ncont_pp > 600.0) {

            cs_real_t rtmp;

            cs_random_normal(1, &rtmp);

            ncont = (int)ncont_pp + sqrt(ncont_pp) * rtmp;

          }
          else {

            cs_random_poisson(1, ncont_pp, &ncont);

          }
          ncont = CS_MAX(1, ncont);
          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_N_SMALL_ASPERITIES,
                                    ncont);

          adhes_energ *= ncont;
          adhes_force *= ncont ;
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE,
                                    adhes_force);
          adhes_force_ps = adhes_force;

          cs_real_t adhes_tor = adhes_force * p_diam * 0.5;
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE,
                                    adhes_tor);

        }

        for (cs_lnum_t id = 1; id < 3; id++) {
          adhes_torque[id] =
            - cs_lagr_particle_get_real(particle, p_am, CS_LAGR_ADHESION_TORQUE)
            * vvue[id] / sqrt(cs_math_pow2(vvue[1]) + cs_math_pow2(vvue[2]) );
        }


        /* Is there direct wall-normal lift-off of the cluster ?  */
        if ((adhes_force + grav_force[0] + lift_force[0] + drag_force[0]) < 0
             && iresusp == 0 ) {

          /* Update of the number and weight of resuspended particles     */
          p_set->n_part_resusp += 1;
          p_set->weight_resusp += p_stat_w;

          if (cs_glob_lagr_boundary_interactions->iflmbd == 1) {

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ires]
              += p_stat_w;

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->iflres]
              += p_stat_w + ( p_stat_w * p_mass / mq->b_f_face_surf[n_f_id]);

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->iflm]
              += - ( p_stat_w * p_mass / mq->b_f_face_surf[n_f_id]);

          }

          /* Update of surface covered and deposit height
           * (for non-rolling particles) */
          if (cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG) ==1 ) {

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->iscovc] -=
              cs_math_pi * cs_math_pow2(p_diam) * p_stat_w
              * 0.25 / mq->b_f_face_surf[n_f_id];

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepm] -=
              cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
              * 0.25 / mq->b_f_face_surf[n_f_id];

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepv] -=
              cs_math_pow2(cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                   * 0.25 / mq->b_f_face_surf[n_f_id]);

            bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->inclg] -=
              p_stat_w;

          }

          /* The particle is resuspended    */
          vpart[0] = CS_MIN(-1.0 / p_mass * dtp
                            * CS_ABS(-lift_force[0] - drag_force[0]
                                     - adhes_force -grav_force[0] ),
                            0.001);
          if (cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG) == 1 ) {
            vpart[1] = 0.0;
            vpart[2] = 0.0;
          }

          cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG, 0);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_FORCE, 0.0);
          cs_lagr_particle_set_real(particle, p_am, CS_LAGR_ADHESION_TORQUE, 0.0);

          if (CS_ABS(p_height-p_diam)/p_diam > 1.0e-6) {
            cs_real_t d_resusp = pow(0.75 * cs_math_pow2(p_diam) * p_height, 1.0/3.0);
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, d_resusp);
            cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, d_resusp);
          }

          if (p_am->count[0][CS_LAGR_N_LARGE_ASPERITIES] > 0)
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_N_LARGE_ASPERITIES, 0);

          if (p_am->count[0][CS_LAGR_N_SMALL_ASPERITIES] > 0)
            cs_lagr_particle_set_lnum(particle, p_am,
                                      CS_LAGR_N_SMALL_ASPERITIES, 0);

          if (p_am->count[0][CS_LAGR_DISPLACEMENT_NORM] > 0)
            cs_lagr_particle_set_real(particle, p_am,
                                      CS_LAGR_DISPLACEMENT_NORM, 0.0);

          iresusp = 1;
        }

        /* No direct normal lift-off */
        else if (iresusp == 0) {

          /* Calculation of the norm of the hydrodynamic
           * torque and drag (tangential) */

          cs_real_t drag_tor_norm  = sqrt (cs_math_pow2(drag_torque[1]) + cs_math_pow2(drag_torque[2]));
          cs_real_t lift_tor_norm  = sqrt (cs_math_pow2(lift_torque[1]) + cs_math_pow2(lift_torque[2]));
          cs_real_t grav_tor_norm  = sqrt (cs_math_pow2(grav_torque[1]) + cs_math_pow2(grav_torque[2]));

          /* Differentiation between two cases:
           *  a) Already rolling clusters
           *  b) Deposited clusters being broken */

          if (cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG) == 2 ) {

            cs_real_t iner_tor = (7.0 / 5.0) * p_mass * cs_math_pow2((p_diam * 0.5));
            cs_real_t cst_1, cst_4;
            cst_4 =   6 * cs_math_pi * visccf
              * romf * 1.7 * 1.4
              * cs_math_pow2(p_diam * 0.5);
            cst_1 = cst_4 * (p_diam * 0.5) / iner_tor;

            for (cs_lnum_t id = 1; id < 3; id++) {

              vpart0[id] = vpart[id];
              vpart[id]  =  vpart0[id] * exp(-cst_1 * dtp)
                + (vvue[id] + (adhes_torque[id] + lift_torque[id] + grav_torque[id]) / cst_4)
                * (1.0 - exp(-cst_1 * dtp) );

            }

            cs_real_t scalax  = vpart[1] * vvue[1] + vpart[2] * vvue[2];

            if (scalax > 0.0) {

              /* The cluster continues to roll (depo_flag = 2)  */

              vpart[0] = 0.0;

              for (cs_lnum_t id = 1; id < 3; id++) {

                if (CS_ABS (vpart[id]) > CS_ABS (vvue[id]))
                  /* The velocity of the rolling particle cannot   */
                  /* exceed the surrounding fluid velocity    */
                  vpart[id] = vvue[id];

                cs_real_t kkk = vvue[id] + adhes_torque[id] / cst_4;
                cs_real_t kk  = vpart0[id] - kkk;

                depl[id] =   (kkk * dtp + kk / cst_1 * (1. - exp(-cst_1 * dtp)));

              }

              iresusp = 1;

            }
            /* if (scalax..) */
            else {

              /* The cluster stops moving
               * the flag is set to 1 and velocity and displacement are null */

              /* Simplified treatment:
               * no update of iscovc, ihdepm for particles that stop */
              iresusp = 1;

              for (cs_lnum_t id = 1; id < 3; id++) {

                depl[id]     = 0.0;
                vpart[id]    = 0.0;

              }

            } /* if (scalax..)   */

          }
          else if (cs_lagr_particle_get_lnum(particle, p_am,
                                             CS_LAGR_DEPOSITION_FLAG) == 1 &&
                   iresusp == 0 ) {
            /* Cluster being broken:
             * we first check if there is resuspension */
            cs_real_t cond_resusp[2];
            for (cs_lnum_t id = 0; id < 2; id++) {
              cond_resusp[id] = (adhes_torque[id] + grav_torque[id] + lift_torque[id]
                                  + drag_torque[id] ) * vvue[id];
            }

            if (cond_resusp[0] > 0.0 || cond_resusp[1] > 0.0) {
              iresusp = 1;
              cs_real_t clust_resusp_height;
              cs_real_t clust_consol_height;
              cs_real_t height_reent;
              cs_real_t random;
              /* Sample of a possible break line */
              if (  cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT)
                  < diam_mean) {
                cs_random_uniform(1, &random);
                clust_resusp_height = random * p_height;
                clust_consol_height = 0.0;
              }
              else {
                cs_real_t param = 1.0 - ((drag_tor_norm + lift_tor_norm - grav_tor_norm)
                                          /p_diam - 2.0 * adhes_force ) /
                  (cs_glob_lagr_consolidation_model->force_consol - adhes_force);
                if (param >= 1.)
                  clust_consol_height = p_height; // Very high consolidation
                else if (param <= -1.)
                  clust_consol_height = 0.; // Very high hydrodynamic forces
                else
                  clust_consol_height =
                    cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CONSOL_HEIGHT)
                    * (1 + cs_glob_lagr_consolidation_model->slope_consol
                       * 0.5 * log((1.0+param)/(1.0-param) ) );
              }
              cs_random_uniform(1, &random);
              height_reent = random * (p_height - clust_consol_height);
              height_reent = CS_MIN(height_reent, p_height);
              height_reent = CS_MAX(height_reent, 0.0);

              /* Treatment of the new rolling particle */
              cs_lnum_t itreated = 0;
              cs_lnum_t nb_part_reent = height_reent / p_height *
                cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART);

              if (nb_part_reent < 1.0 && itreated == 0) {
                /* No resuspension (cluster too small)*/
                itreated = 1;

                for (cs_lnum_t id = 0; id < 3; id++) {
                  vpart[id] = 0.0;
                  depl[id]  = 0.0;
                }
              }
              else if ((cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART)
                         -nb_part_reent) < 1.0 && itreated == 0) {
                /* The whole cluster starts rolling*/
                cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG, 2);

                /* Update of surface covered and deposit height */
                bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->iscovc] -=
                  cs_math_pi * cs_math_pow2(p_diam) * p_stat_w
                  * 0.25 / mq->b_f_face_surf[n_f_id];

                bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepm] -=
                  cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                  * 0.25 / mq->b_f_face_surf[n_f_id];

                bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepv] -=
                  cs_math_pow2(cs_math_pi * p_height * cs_math_pow2(p_diam) * p_stat_w
                       * 0.25 / mq->b_f_face_surf[n_f_id]);

                bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->inclg] -=
                  p_stat_w;

                itreated = 1;
                cs_real_t d_resusp = pow(0.75 * cs_math_pow2(p_diam) * p_height, 1.0/3.0);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_DIAMETER, d_resusp);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, d_resusp);

                /* Treatment of cluster motion */
                cs_real_t iner_tor = (7.0 / 5.0) * p_mass * cs_math_pow2((p_diam * 0.5));
                cs_real_t cst_1, cst_4;
                cst_4 =   6 * cs_math_pi * visccf
                  * romf * 1.7 * 1.4
                  * cs_math_pow2(p_diam * 0.5);
                cst_1 = cst_4 * (p_diam * 0.5) / iner_tor;

                for (cs_lnum_t id = 1; id < 3; id++) {
                  vpart0[id] = vpart[id];
                  vpart[id]  =  vpart0[id] * exp(-cst_1 * dtp)
                    + (vvue[id] + (adhes_torque[id] + lift_torque[id] + grav_torque[id]) / cst_4)
                    * (1.0 - exp(-cst_1 * dtp) );
                }

                vpart[0] = 0.0;
                for (cs_lnum_t id = 1; id < 3; id++) {
                  if (CS_ABS (vpart[id]) > CS_ABS (vvue[id]))
                    /* The velocity of the rolling particle cannot   */
                    /* exceed the surrounding fluid velocity    */
                    vpart[id] = vvue[id];
                  cs_real_t kkk = vvue[id] + adhes_torque[id] / cst_4;
                  cs_real_t kk  = vpart0[id] - kkk;
                  depl[id] =   (kkk * dtp + kk / cst_1 * (1. - exp(-cst_1 * dtp)));
                }

              }
              else if (itreated == 0) {
                /* Duplication of the particle */
                *nresnew = *nresnew + 1;
                cs_lagr_particle_set_resize(p_set->n_particles + *nresnew);
                cs_lagr_part_copy(p_set->n_particles + *nresnew, ip);

                /* We split both particles:
                * Part ip stays while the new one starts rolling */
                unsigned char *new_part = p_set->p_buffer
                  + p_am->extents * p_set->n_particles+*nresnew;
                cs_real_t nb_resusp =  height_reent / p_height
                  * cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_CLUSTER_NB_PART, nb_resusp);
                cs_real_t m_resusp = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_MASS)
                  * cs_lagr_particle_get_real(new_part, p_am, CS_LAGR_CLUSTER_NB_PART)
                  / cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART) ;
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_MASS, m_resusp);
                cs_real_t d_resusp = pow(0.75 * cs_math_pow2(p_diam) * p_height, 1.0/3.0);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_DIAMETER, d_resusp);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_HEIGHT, d_resusp);

                cs_lagr_particle_set_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG, 1);
                vpart[0] = 0.0;
                for (cs_lnum_t id = 1; id < 3; id++) {
                  vpart[id] = 0.0;
                  depl[id]  = 0.0;
                }

                /* Update of deposit height */
                cs_real_t d_stay = cs_lagr_particle_get_real(particle, p_am, CS_LAGR_HEIGHT);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_HEIGHT, d_stay);

                bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepm] -=
                  cs_math_pi * height_reent * cs_math_pow2(p_diam) * p_stat_w
                  * 0.25 / mq->b_f_face_surf[n_f_id];

                bound_stat[n_f_id + nfabor * cs_glob_lagr_boundary_interactions->ihdepv] -=
                  cs_math_pow2(cs_math_pi * height_reent * cs_math_pow2(p_diam) * p_stat_w
                       * 0.25 / mq->b_f_face_surf[n_f_id]);

                cs_real_t nb_stay =
                  cs_lagr_particle_get_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART) -
                  cs_lagr_particle_get_real(new_part, p_am, CS_LAGR_CLUSTER_NB_PART);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_CLUSTER_NB_PART, nb_stay);

                cs_real_t mp_stay = p_mass - cs_lagr_particle_get_real(new_part, p_am, CS_LAGR_MASS);
                cs_lagr_particle_set_real(particle, p_am, CS_LAGR_MASS, mp_stay);

                /* The new particle starts rolling */
                cs_lagr_particle_set_lnum(new_part, p_am, CS_LAGR_DEPOSITION_FLAG, 2);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_ADHESION_FORCE, adhes_force);
                cs_lagr_particle_set_lnum(new_part, p_am, CS_LAGR_N_SMALL_ASPERITIES, CS_MAX(1,ncont_pp));
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_ADHESION_TORQUE, adhes_force * d_resusp * 0.5);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_DEPO_TIME, 0.0);
                cs_lagr_particle_set_real(new_part, p_am, CS_LAGR_CONSOL_HEIGHT, 0.0);

                itreated = 1;
              }

            } /* End of condition for resuspension */

            else {
              /* No rolling occurs */
              iresusp = 1;
              for (cs_lnum_t id = 1; id < 3; id++) {
                vpart[id] = 0.0;
                depl[id]  = 0.0;
              }
            }

          }

        }

      } /* End of multilayer resuspension */

    }

  }
  /* if ireent.eq.0 --> Motionless deposited particle   */
  else {

    if (cs_lagr_particle_get_lnum(particle, p_am,
                                  CS_LAGR_DEPOSITION_FLAG) != CS_LAGR_PART_IN_FLOW) {

      for (cs_lnum_t id = 1; id < 3; id++) {

        vpart[id] = 0.0;
        vvue[id]  = 0.0;
        depl[id]  = 0.0;

      }

    }

  }

  /* ===========================================================================
   * 3. Reference frame change:
   * --------------------------
   * local reference frame for the boundary face --> global reference frame
   * ======================================================================== */

  /* 3.1 - Displacement   */

  cs_real_t depg[3];

  cs_math_33t_3_product(rot_m, depl, depg);

  /* 3.2 - Particle velocity   */

  cs_real_t *part_vel = cs_lagr_particle_attr(particle, p_am,
                                              CS_LAGR_VELOCITY);

  cs_math_33t_3_product(rot_m, vpart, part_vel);

  /* 3.3 - flow-seen velocity  */

  cs_real_t *part_vel_seen = cs_lagr_particle_attr(particle, p_am,
                                                   CS_LAGR_VELOCITY_SEEN);

  cs_math_33t_3_product(rot_m, vvue, part_vel_seen);

  /* ===========================================================================
   * 5. Computation of the new particle position
   * ======================================================================== */

  cs_real_t *part_coords = cs_lagr_particle_attr(particle, p_am,
                                                 CS_LAGR_COORDS);
  for (cs_lnum_t id = 0 ; id < 3; id++)
    part_coords[id] += depg[id];
}

/*----------------------------------------------------------------------------*/
/*! \brief Deposition submodel
 *
 *  Main subroutine of the submodel
 *   1/ Calculation of the normalized wall-normal distance of
 *           the boundary-cell particles
 *   2/ Sorting of the particles with respect to their normalized
 *           wall-normal distance
 *         * If y^+ > depint : the standard Langevin model is applied
 *         * If y^+ < depint : the deposition submodel is applied
 *
 * \param[in] dtp       time step
 * \param[in] taup      dynamic characteristic time
 * \param[in] tlag      fluid characteristic time
 * \param[in] piil      term in integration of UP SDEs
 * \param[in] bx        turbulence characteristics
 * \param[in] vagaus    gaussian random variables
 * \param[in] gradpr    pressure gradient
 * \param[in] romp      particles associated density
 * \param[in] force_p   taup times forces on particles (m/s)
 * \param[in] vislen    FIXME
 */
/*----------------------------------------------------------------------------*/

static void
_lagdep(cs_real_t           dtp,
        const cs_real_t     taup[],
        const cs_real_3_t   tlag[],
        const cs_real_3_t   piil[],
        const cs_real_33_t  bx[],
        const cs_real_33_t  vagaus[],
        const cs_real_3_t   gradpr[],
        const cs_real_t     romp[],
        const cs_real_3_t   force_p[],
        const cs_real_t     vislen[],
        cs_lnum_t          *nresnew)
{
  /* Particles management */
  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am = p_set->p_am;

  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();

  /* Initialisations*/

  cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;

  cs_real_t vitf = 0.0;

  cs_real_t aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9, aux10, aux11;
  cs_real_t ter1f, ter2f, ter3f;
  cs_real_t ter1p, ter2p, ter3p, ter4p, ter5p;
  cs_real_t ter1x, ter2x, ter3x, ter4x, ter5x;
  cs_real_t p11, p21, p22, p31, p32, p33;
  cs_real_t omega2, gama2, omegam;
  cs_real_t grga2, gagam, gaome;

  cs_lnum_t nor = cs_glob_lagr_time_step->nor;

  /* Interface location between near-wall region   */
  /* and core of the flow (normalized units)  */

  cs_real_t depint      = 100.0;

  /* loop on the particles  */
  for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

    unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

    cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);

    if (cell_id >= 0 &&
        cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
        != CS_LAGR_PART_IMPOSED_MOTION) {

      cs_real_t *old_part_vel      = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                             CS_LAGR_VELOCITY);
      cs_real_t *old_part_vel_seen = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                             CS_LAGR_VELOCITY_SEEN);
      cs_real_t *part_vel          = cs_lagr_particle_attr(particle, p_am,
                                                           CS_LAGR_VELOCITY);
      cs_real_t *part_vel_seen     = cs_lagr_particle_attr(particle, p_am,
                                                           CS_LAGR_VELOCITY_SEEN);
      cs_real_t *part_coords       = cs_lagr_particle_attr(particle, p_am,
                                                           CS_LAGR_COORDS);
      cs_real_t *old_part_coords   = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                             CS_LAGR_COORDS);

      cs_real_t romf = extra->cromf->val[cell_id];

      /* Fluid temperature computation depending on the type of flow  */
      cs_real_t tempf;

      if (   cs_glob_physical_model_flag[CS_COMBUSTION_COAL] >= 0
          || cs_glob_physical_model_flag[CS_COMBUSTION_PCLC] >= 0
          || cs_glob_physical_model_flag[CS_COMBUSTION_FUEL] >= 0)
        tempf = extra->t_gaz->val[cell_id];

      else if (   cs_glob_physical_model_flag[CS_COMBUSTION_3PT] >= 0
               || cs_glob_physical_model_flag[CS_COMBUSTION_EBU] >= 0
               || cs_glob_physical_model_flag[CS_ELECTRIC_ARCS] >= 0
               || cs_glob_physical_model_flag[CS_JOULE_EFFECT] >= 0)
        tempf = extra->temperature->val[cell_id];

      else if (   cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE
               && cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_CELSIUS)
        tempf = extra->scal_t->val[cell_id] + tkelvi;

      else if (   cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_TEMPERATURE
               && cs_glob_thermal_model->itpscl == CS_TEMPERATURE_SCALE_KELVIN)
        tempf = extra->scal_t->val[cell_id];

      else if (cs_glob_thermal_model->itherm == CS_THERMAL_MODEL_ENTHALPY) {

        cs_lnum_t mode  = 1;
        CS_PROCF (usthht,USTHHT) (&mode, &(extra->scal_t->val[cell_id]), &tempf);

        tempf = tempf + tkelvi;

      }

      else
        tempf = cs_glob_fluid_properties->t0;

      /* If y^+ is greater than the interface location,
         the standard model is applied
         ============================================== */

      if (cs_lagr_particle_get_real(particle, p_am, CS_LAGR_YPLUS) > depint &&
          cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
          == CS_LAGR_PART_IN_FLOW) {

        cs_lagr_particle_set_lnum(particle,
                                  p_am,
                                  CS_LAGR_MARKO_VALUE,
                                  CS_LAGR_COHERENCE_STRUCT_BULK);

        for (cs_lnum_t id = 0; id < 3; id++) {

          vitf = extra->vel->vals[1][cell_id * 3 + id];

          /* --> (2.1) Calcul preliminaires :    */
          /* ----------------------------   */
          /* compute II*TL+<u> and [(grad<P>/rhop+g)*tau_p+<Uf>] ?  */

          cs_real_t tci = piil[ip][id] * tlag[ip][id] + vitf;
          cs_real_t force = force_p[ip][id];

          /* --> (2.2) Calcul des coefficients/termes deterministes */
          /* ----------------------------------------------------    */

          aux1 = exp(-dtp / taup[ip]);
          aux2 = exp(-dtp / tlag[ip][id]);
          aux3 = tlag[ip][id] / (tlag[ip][id] - taup[ip]);
          aux4 = tlag[ip][id] / (tlag[ip][id] + taup[ip]);
          aux5 = tlag[ip][id] * (1.0 - aux2);
          aux6 = cs_math_pow2(bx[ip][id][nor-1]) * tlag[ip][id];
          aux7 = tlag[ip][id] - taup[ip];
          aux8 = cs_math_pow2(bx[ip][id][nor-1]) * cs_math_pow2(aux3);

          /* --> trajectory terms */
          cs_real_t aa = taup[ip] * (1.0 - aux1);
          cs_real_t bb = (aux5 - aa) * aux3;
          cs_real_t cc = dtp - aa - bb;

          ter1x = aa * old_part_vel[id];
          ter2x = bb * old_part_vel_seen[id];
          ter3x = cc * tci;
          ter4x = (dtp - aa) * force;

          /* --> flow-seen velocity terms   */
          ter1f = old_part_vel_seen[id] * aux2;
          ter2f = tci * (1.0 - aux2);

          /* --> termes pour la vitesse des particules     */
          cs_real_t dd = aux3 * (aux2 - aux1);
          cs_real_t ee = 1.0 - aux1;

          ter1p = old_part_vel[id] * aux1;
          ter2p = old_part_vel_seen[id] * dd;
          ter3p = tci * (ee - dd);
          ter4p = force * ee;

          /* --> (2.3) Coefficients computation for the stochastic integral    */
          /* --> Integral for particles position */
          gama2  = 0.5 * (1.0 - aux2 * aux2);
          omegam = aux3 * ( (tlag[ip][id] - taup[ip]) * (1.0 - aux2)
                                  - 0.5 * tlag[ip][id] * (1.0 - aux2 * aux2)
                                  + cs_math_pow2(taup[ip]) / (tlag[ip][id] + taup[ip])
                                  * (1.0 - aux1 * aux2)
                                  ) * aux6;
          omega2 =  aux7 * (aux7 * dtp - 2.0 * (tlag[ip][id] * aux5 - taup[ip] * aa))
                   + 0.5 * tlag[ip][id] * tlag[ip][id] * aux5 * (1.0 + aux2)
                   + 0.5 * taup[ip] * taup[ip] * aa * (1.0 + aux1)
                   - 2.0 * aux4 * tlag[ip][id] * taup[ip] * taup[ip] * (1.0 - aux1 * aux2);
          omega2 = aux8 * omega2;

          if (CS_ABS(gama2) > cs_math_epzero) {

            p21 = omegam / sqrt(gama2);
            p22 = omega2 - cs_math_pow2(p21);
            p22 = sqrt(CS_MAX(0.0, p22));

          }
          else {

            p21 = 0.0;
            p22 = 0.0;

          }

          ter5x = p21 * vagaus[ip][id][0] + p22 * vagaus[ip][id][1];

          /* --> integral for the flow-seen velocity  */
          p11   = sqrt(gama2 * aux6);
          ter3f = p11 * vagaus[ip][id][0];

          /* --> integral for the particles velocity  */
          aux9  = 0.5 * tlag[ip][id] * (1.0 - aux2 * aux2);
          aux10 = 0.5 * taup[ip] * (1.0 - aux1 * aux1);
          aux11 =   taup[ip] * tlag[ip][id]
                  * (1.0 - aux1 * aux2)
                  / (taup[ip] + tlag[ip][id]);

          grga2 = (aux9 - 2.0 * aux11 + aux10) * aux8;
          gagam = (aux9 - aux11) * (aux8 / aux3);
          gaome = ( (tlag[ip][id] - taup[ip]) * (aux5 - aa)
                    - tlag[ip][id] * aux9
                    - taup[ip] * aux10
                    + (tlag[ip][id] + taup[ip]) * aux11)
                  * aux8;

          if (p11 > cs_math_epzero)
            p31 = gagam / p11;
          else
            p31 = 0.0;

          if (p22 > cs_math_epzero)
            p32 = (gaome - p31 * p21) / p22;
          else
            p32 = 0.0;

          p33 = grga2 - cs_math_pow2(p31) - cs_math_pow2(p32);
          p33 = sqrt(CS_MAX(0.0, p33));
          ter5p =   p31 * vagaus[ip][id][0]
                  + p32 * vagaus[ip][id][1]
                  + p33 * vagaus[ip][id][2];

          /*  Update of the particle state-vector     */

          part_coords[id] =   old_part_coords[id]
                            + ter1x + ter2x + ter3x + ter4x + ter5x;

          part_vel_seen[id] =  ter1f + ter2f + ter3f;

          part_vel[id]      = ter1p + ter2p + ter3p + ter4p + ter5p;

        }

      }

      /* Otherwise, the deposition submodel is applied
       * ============================================= */

      else if (cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_DEPOSITION_FLAG)
          != CS_LAGR_PART_TO_DELETE) {

        if (  cs_lagr_particle_get_real(particle, p_am,
                                        CS_LAGR_YPLUS)
            < cs_lagr_particle_get_real(particle, p_am,
                                        CS_LAGR_INTERF) ) {

          if (cs_lagr_particle_get_lnum(particle,
                                        p_am,
                                        CS_LAGR_MARKO_VALUE) < 0)
            cs_lagr_particle_set_lnum(particle,
                                      p_am,
                                      CS_LAGR_MARKO_VALUE,
                                      CS_LAGR_COHERENCE_STRUCT_DEGEN_INNER_ZONE_DIFF);
          else
            cs_lagr_particle_set_lnum(particle,
                                      p_am,
                                      CS_LAGR_MARKO_VALUE,
                                      CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF);

        }
        else {

          if (cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_MARKO_VALUE) < 0)
            cs_lagr_particle_set_lnum(particle,
                                      p_am,
                                      CS_LAGR_MARKO_VALUE,
                                      CS_LAGR_COHERENCE_STRUCT_DEGEN_SWEEP);

          else if (cs_lagr_particle_get_lnum(particle,
                                             p_am,
                                             CS_LAGR_MARKO_VALUE)
              == CS_LAGR_COHERENCE_STRUCT_INNER_ZONE_DIFF)
            cs_lagr_particle_set_lnum(particle,
                                      p_am,
                                      CS_LAGR_MARKO_VALUE,
                                      CS_LAGR_COHERENCE_STRUCT_DEGEN_EJECTION);

        }

        _lagesd(dtp,
                ip,
                taup,
                piil,
                vagaus,
                gradpr,
                romp,
                force_p,
                tempf,
                vislen,
                &depint,
                nresnew);

      }

    }

    /* Specific treatment for particles with
     * DEPOSITION_FLAG == CS_LAGR_PART_IMPOSED_MOTION */
    else if (cell_id >= 0) {
      cs_real_t disp[3] = {0., 0., 0.};

      cs_real_t *old_part_coords = cs_lagr_particle_attr_n(particle, p_am, 1,
                                                           CS_LAGR_COORDS);
      cs_real_t *part_coords = cs_lagr_particle_attr(particle, p_am,
                                                     CS_LAGR_COORDS);

      cs_real_t *part_vel_seen = cs_lagr_particle_attr(particle, p_am,
                                                       CS_LAGR_VELOCITY_SEEN);

      cs_real_t *part_vel = cs_lagr_particle_attr(particle, p_am,
                                                  CS_LAGR_VELOCITY);

      cs_user_lagr_imposed_motion(old_part_coords,
                                  dtp,
                                  disp);

      for (cs_lnum_t id = 0; id < 3; id++) {

        part_coords[id] = old_part_coords[id] + disp[id];

        part_vel_seen[id] =  0.0;

        part_vel[id] = disp[id]/dtp;

      }
    }
  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of particle equations of motion:
 *
 * - Standard Model : First or second order
 * - Deposition submodel (Guingo & Minier, 2008) if needed
 *
 * \param[in]  dt_p      lagrangian time step
 * \param[in]  taup      dynamic characteristic time
 * \param[in]  tlag      fluid characteristic time
 * \param[in]  piil      term in integration of U-P SDEs
 * \param[in]  bx        turbulence characteristics
 * \param[out] tsfext    info for return coupling source terms
 * \param[in]  gradpr    pressure gradient
 * \param[in]  gradvf    fluid velocity gradient
 * \param[out] terbru    FIXME
 * \param[in]  vislen    FIXME
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_sde(cs_real_t           dt_p,
            const cs_real_t     taup[],
            const cs_real_3_t   tlag[],
            const cs_real_3_t   piil[],
            const cs_real_33_t  bx[],
            cs_real_t           tsfext[],
            const cs_real_3_t   gradpr[],
            const cs_real_33_t  gradvf[],
            cs_real_t           terbru[],
            const cs_real_t     vislen[],
            cs_lnum_t          *nresnew)
{
  cs_real_t *romp;

  cs_lagr_particle_set_t  *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t *p_am = p_set->p_am;

  BFT_MALLOC(romp, p_set->n_particles, cs_real_t);

  /* Allocate temporay arrays  */
  cs_real_33_t *vagaus;
  BFT_MALLOC(vagaus, p_set->n_particles, cs_real_33_t);

  /* Random values */

  if (cs_glob_lagr_time_scheme->idistu == 1) {
    if (cs_glob_lagr_time_step->nor > 1) {
      for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
        unsigned char *particle = p_set->p_buffer + p_am->extents * ip;
        cs_real_t *_v_gauss = cs_lagr_particle_attr(particle, p_am,
                                                    CS_LAGR_V_GAUSS);
        for (cs_lnum_t id = 0; id < 3; id++) {
          for (cs_lnum_t ivf = 0; ivf < 3; ivf++)
            vagaus[ip][id][ivf] = _v_gauss[id*3 + ivf];
        }
      }
    }
    else {
      for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++)
        cs_random_normal(9, &(vagaus[ip][0][0]));
    }
  }

  else {
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
      for (cs_lnum_t id = 0; id < 3; id++) {
        for (cs_lnum_t ivf = 0; ivf < 3; ivf++)
          vagaus[ip][id][ivf] = 0.0;
      }
    }
  }

  /* Brownian movement */
  cs_real_t *brgaus = NULL;

  if (cs_glob_lagr_brownian->lamvbr == 1) {
    BFT_MALLOC(brgaus, p_set->n_particles*6, cs_real_t);
    if (cs_glob_lagr_time_step->nor > 1) {
      for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
        unsigned char *particle = p_set->p_buffer + p_am->extents * ip;
        cs_real_t *_br_gauss = cs_lagr_particle_attr(particle, p_am,
                                                     CS_LAGR_BR_GAUSS);
        for (cs_lnum_t id = 0; id < 6; id++)
          brgaus[ip*6 + id] = _br_gauss[id];
      }
    }
    else {
      for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++)
        cs_random_normal(6, &(brgaus[6 * ip]));
    }
  }

  /* Computation of particle density */
  cs_real_t aa = 6.0 / cs_math_pi;

  for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

    unsigned char *particle = p_set->p_buffer + p_am->extents * ip;
    if (cs_lagr_particle_get_cell_id(particle, p_am) >= 0) {

      cs_real_t d3 = cs_math_pow3(cs_lagr_particle_get_real(particle, p_am,
                                                   CS_LAGR_DIAMETER));
      romp[ip] = aa * cs_lagr_particle_get_real(particle, p_am,
                                                CS_LAGR_MASS) / d3;

    }

  }

  /* Management of user external force fields
     ---------------------------------------- */

  cs_real_3_t *force_p;
  BFT_MALLOC(force_p, p_set->n_particles, cs_real_3_t);

  for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
    force_p[ip][0] = 0.0;
    force_p[ip][1] = 0.0;
    force_p[ip][2] = 0.0;
  }

  cs_user_lagr_ef(dt_p,
                  (const cs_real_t *)taup,
                  (const cs_real_3_t *)tlag,
                  (const cs_real_3_t *)piil,
                  (const cs_real_33_t *)bx,
                  (const cs_real_t *)tsfext,
                  (const cs_real_33_t *)vagaus,
                  (const cs_real_3_t *)gradpr,
                  (const cs_real_33_t *)gradvf,
                  romp,
                  force_p);

  const cs_real_t  *grav  = cs_glob_physical_constants->gravity;
  cs_lagr_extra_module_t *extra = cs_get_lagr_extra_module();
  cs_real_t added_mass_const = cs_glob_lagr_time_scheme->added_mass_const;

  /* Finalize forces on particles:
   *
   *  (- pressure gradient /romp +ext forces + g) . taup
   *
   * */
  if (cs_glob_lagr_time_scheme->iadded_mass == 0) {
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;
      cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);
      for (int id = 0; id < 3; id++) {
        force_p[ip][id] = (- gradpr[cell_id][id] / romp[ip]
          + grav[id] + force_p[ip][id]) * taup[ip];

      }
    }
  }
  /* Added-mass term?     */
  else {
    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;
      cs_lnum_t cell_id = cs_lagr_particle_get_cell_id(particle, p_am);
      cs_real_t romf = extra->cromf->val[cell_id];
      for (int id = 0; id < 3; id++) {
        force_p[ip][id] = (- gradpr[cell_id][id] / romp[ip]
          * (1.0 + 0.5 * added_mass_const)
          / (1.0 + 0.5 * added_mass_const * romf / romp[ip])
          + grav[id] + force_p[ip][id]) * taup[ip];
      }
    }
  }

  /* First order
     ----------- */

  if (cs_glob_lagr_time_scheme->t_order == 1) {

    /* If no deposition sub-model is activated, call of subroutine lages1
       for every particle */

    if (cs_glob_lagr_model->deposition <= 0)
      _lages1(dt_p,
              taup,
              tlag,
              piil,
              bx,
              (const cs_real_33_t *)vagaus,
              brgaus,
              gradpr,
              romp,
              force_p,
              terbru);

    /* Management of the deposition submodel */

    else
      _lagdep(dt_p,
              taup,
              tlag,
              piil,
              bx,
              (const cs_real_33_t *)vagaus,
              gradpr,
              romp,
              force_p,
              vislen,
              nresnew);

  }

  /* Second order
     ------------ */

  else {

    _lages2(dt_p,
            taup,
            tlag,
            piil,
            bx,
            tsfext,
            (const cs_real_33_t *)vagaus,
            brgaus,
            gradpr,
            romp,
            force_p,
            terbru);

    /* Save Gaussian variable if needed */
    if (cs_glob_lagr_time_step->nor == 1) {

      for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {
        unsigned char *particle = p_set->p_buffer + p_am->extents * ip;
        if (cs_glob_lagr_time_scheme->idistu == 1) {
          cs_real_t *_v_gauss
            = cs_lagr_particle_attr(particle, p_am, CS_LAGR_V_GAUSS);
          for (cs_lnum_t id = 0; id < 3; id++) {
            for (cs_lnum_t ivf = 0; ivf < 3; ivf++)
              _v_gauss[id*3 + ivf] = vagaus[ip][id][ivf];
          }
        }
        if (cs_glob_lagr_brownian->lamvbr == 1) {
          cs_real_t *_br_gauss
            = cs_lagr_particle_attr(particle, p_am, CS_LAGR_BR_GAUSS);
          for (cs_lnum_t id = 0; id < 6; id++)
            _br_gauss[id] = brgaus[ip*6 + id];
        }
      }

    }

  }

  BFT_FREE(force_p);

  BFT_FREE(brgaus);
  BFT_FREE(vagaus);
  BFT_FREE(romp);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Integration of a stochastic differential equation (SDE) for
 *        a user particle variable (attribute).
 *
 * \f[
 *  \frac{dV}{dt} = \frac{V - PIP}{TCARAC}
 * \f]
 *
 * When there is interaction with a boundary face, the integration
 * degenerates to order 1 (even if the 2nd order scheme is active).
 *
 * \param[in]  attr    attribute/variable
 * \param[in]  tcarac  variable characteristic time
 * \param[in]  pip     right-hand side associated with SDE
 */
/*----------------------------------------------------------------------------*/

void
cs_lagr_sde_attr(cs_lagr_attribute_t   attr,
                 cs_real_t            *tcarac,
                 cs_real_t            *pip)
{
  /* Particles management */
  cs_lagr_particle_set_t         *p_set = cs_glob_lagr_particle_set;
  const cs_lagr_attribute_map_t  *p_am  = p_set->p_am;

  int ltsvar = 0;

  if (p_set->p_am->source_term_displ != NULL) {
    if (p_set->p_am->source_term_displ[attr] >= 0)
      ltsvar = 1;
  }

  int nor = cs_glob_lagr_time_step->nor;

  assert(nor == 1 || nor == 2);

  if (nor == 1) {

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      if (cs_lagr_particle_get_cell_id(particle, p_am) >= 0) {

        if (tcarac[ip] <= 0.0)
          bft_error
            (__FILE__, __LINE__, 0,
             _("The characteristic time for the stochastic differential equation\n"
               "of variable %d should be > 0.\n\n"
               "Here, for particle %d, its value is %e11.4."),
             attr, ip, tcarac[ip]);

        cs_real_t aux1 = cs_glob_lagr_time_step->dtp/tcarac[ip];
        cs_real_t aux2 = exp(-aux1);
        cs_real_t ter1 = cs_lagr_particle_get_real_n(particle, p_am, 1, attr)*aux2;
        cs_real_t ter2 = pip[ip] * (1.0 - aux2);

        /* Pour le cas NORDRE= 1 ou s'il y a rebond,     */
        /* le ETTP suivant est le resultat final    */
        cs_lagr_particle_set_real(particle, p_am, attr, ter1 + ter2);

        /* Pour le cas NORDRE= 2, on calcule en plus TSVAR pour NOR= 2  */
        if (ltsvar) {
          cs_real_t *part_ptsvar = cs_lagr_particles_source_terms(p_set, ip, attr);
          cs_real_t ter3 = (-aux2 + (1.0 - aux2) / aux1) * pip[ip];
          *part_ptsvar = 0.5 * ter1 + ter3;

        }

      }

    }

  }
  else if (nor == 2) {

    for (cs_lnum_t ip = 0; ip < p_set->n_particles; ip++) {

      unsigned char *particle = p_set->p_buffer + p_am->extents * ip;

      if (   cs_lagr_particle_get_cell_id(particle, p_am) >= 0
          && cs_lagr_particle_get_lnum(particle, p_am, CS_LAGR_REBOUND_ID) != 0) {

        if (tcarac [ip] <= 0.0)
          bft_error
            (__FILE__, __LINE__, 0,
             _("The characteristic time for the stochastic differential equation\n"
               "of variable %d should be > 0.\n\n"
               "Here, for particle %d, its value is %e11.4."),
             attr, ip, tcarac[ip]);

        cs_real_t aux1   = cs_glob_lagr_time_step->dtp / tcarac [ip];
        cs_real_t aux2   = exp(-aux1);
        cs_real_t ter1   = 0.5 * cs_lagr_particle_get_real_n(particle, p_am, 1,
                                                             attr) * aux2;
        cs_real_t ter2   = pip [ip] * (1.0 - (1.0 - aux2) / aux1);

        /* Pour le cas NORDRE= 2, le ETTP suivant est le resultat final */
        cs_real_t *part_ptsvar = cs_lagr_particles_source_terms(p_set, ip, attr);
        cs_lagr_particle_set_real(particle, p_am, attr,
                                  *part_ptsvar + ter1 + ter2 );

      }

    }

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
