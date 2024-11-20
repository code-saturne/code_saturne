#ifndef __CS_WALL_CONDENSATION_H__
#define __CS_WALL_CONDENSATION_H__

/*============================================================================
 * Base wall condensation model.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

typedef enum {
  CS_WALL_COND_MODEL_COPAIN    = 0,
  CS_WALL_COND_MODEL_COPAIN_BD = 1,
  CS_WALL_COND_MODEL_UCHIDA    = 2,
  CS_WALL_COND_MODEL_DEHBI     = 3,
} cs_wall_cond_natural_conv_model_t;

typedef enum {
  CS_WALL_COND_MODEL_WALL_LAW    = 0,
  CS_WALL_COND_MODEL_SCHLICHTING = 1
} cs_wall_cond_forced_conv_model_t;

typedef enum {
  CS_WALL_COND_MIXED_MAX       = 0,
  CS_WALL_COND_MIXED_INCROPERA = 1
} cs_wall_cond_mixed_conv_model_t;

typedef struct {

  int icondb; /*! Switch used to activate wall condensation (0 : activated) */
  int icondv; /*! Switch used to activate wall condensation
                  with metal structures(0 : activated) */
  int nztag1d; /*! Indicate if the thermal 1D model of severe accident is used.
                   1 a 1D thermal equation is used, 0 user defined dirichlet */

  // Model type
  cs_wall_cond_natural_conv_model_t natural_conv_model;
  cs_wall_cond_forced_conv_model_t  forced_conv_model;
  cs_wall_cond_mixed_conv_model_t   mixed_conv_model;

  /* Surface wall condensation */

  // Mesh related quantities
  /*! number of faces in which a condensation source terms */
  cs_lnum_t  nfbpcd;
  /*! faces in which a condensation source terms */
  cs_lnum_t *ifbpcd;
  /*! type of condensation source terms for each variable
       - 0 for an variable at ambient value,
       - 1 for an variable at imposed value.
       See \ref cs_user_wall_condensation*/
  cs_lnum_t *itypcd;
  /*! list of the zones associated to the faces whith a
    condensation source term */
  cs_lnum_t *izzftcd;
  /*! value of the condensation source terms for pressure.
      For the other variables, eventual imposed specific value.
      See \ref cs_user_wall_condensation*/
  cs_real_t *spcond;
  /*! value of the thermal exchange coefficient associated to
      the condensation model used.
       \ref cs_user_wall_condensation */
  cs_real_t *hpcond;
  /*! Temperature at condensing wall faces (for post-processing purposes) */
  cs_real_t *twall_cond;
  /*! value of the thermal flux for the condensation model.
      See \ref cs_user_wall_condensation */
  cs_real_t *thermal_condensation_flux;
  /*! */
  cs_real_t *convective_htc;
  /*! */
  cs_real_t *condensation_htc;
  /*! */
  cs_real_t *total_htc;
  /*! external heat flux used as flux conditions
      of the 1d thermal model (in unit \f$W.m^{-2}\f$). */
  cs_real_t *flthr;
  /*! external heat flux derivative used as flux conditions
      of the 1d thermal model (in unit \f$W.m^{-3}\f$) */
  cs_real_t *dflthr;

  // Zone related quantities
  /*! number of the zones with a specific condensation source terms
      depending on the wall temperature and material properties.
       by default (nzones = 1) if the user does not specified different zones
       in the user \ref cs_user_wall_condensation.c .*/
  cs_lnum_t  nzones;
  /*! compute method for the exchange coefficient of the
      condensation source term used by the copain model.
      - 1: the turbulent exchange coefficient of the flow
      - 2: the exchange coefficient of the copain correlation
      - 3: the maximal value between the two previous exchange coefficients */
  cs_lnum_t *izcophc;
  /*! compute method for the thermal exchange coefficient associated
      to the heat transfer to the wall due to the condensation phenomenon.
        - 2: the thermal exchange coefficient of the copain correlation
        - 3: the maximal value between the current and previous thermal
             exchange coefficient evaluated by the copain correlation */
  cs_lnum_t *izcophg;
  /*! compute method for the wall temperature at the solid/fluid interface
      coupled with condensation to the wall
      - 1: the wall temperature is computed with a 1-D thermal model
           with implicit numerical scheme
      - 0: the wall temperature is imposed as constant by the user (default)
           exchange coefficient evaluated by the copain correlation */
  cs_lnum_t *iztag1d;
  /*! Constant value of the wall temperature given by the user when
      the thermal 1D model is not activated for the condensation model with
      different zones specified in the user \ref cs_user_wall_condensation. */
  cs_real_t *ztpar;
  /*! Coordinates of the reference point for forced and mixed
    convection regimes */
  cs_real_t *zxrefcond;
  cs_real_t *zprojcond;

  /* Volume wall condensation */

  /*! number of the cells in which a condensation source terms is imposed.
      See  cs_user_wall_condensation */
  cs_lnum_t  ncmast;
  /*! number of the volume strutures with a specific condensation source terms
      depending on the wall temperature and material properties.
      by default (nvolumes = 1) if the user does not specified different volumes
      in the user function \ref cs_user_wall_condensation */
  cs_lnum_t  nvolumes;
  /*! list on the ncmast cells in which a condensation source terms is imposed.
       See  the user function \ref cs_user_wall_condensation. */
  cs_lnum_t *ltmast;
  /*! type of condensation source terms for each variable
          - 0 for a variable at ambient value,
          - 1 for a variable at imposed value.
      See the user function \ref  cs_user_wall_condensation. */
  cs_lnum_t *itypst;
  /*! zone type where a condensation source terms is imposed to model
      the metal structures condensation on a volumic zone. */
  cs_lnum_t *izmast;
  /*! value of the condensation source terms for pressure
      associated to the metal structures modelling.
      For the other variables, eventual imposed specific value.
      See the user function \ref cs_user_wall_condensation. */
  cs_real_t *svcond;
  /*! value of the thermal flux for the condensation model
      associated to the metal structures modelling.
      See the user function \ref cs_user_wall_condensation. */
  cs_real_t *flxmst;
  /*! compute method for the wall temperature at the solid/fluid
    interface coupled with condensation to the metal mass structures wall
    - 1: the wall temperature is computed with a 0-D thermal model
         with explicit numerical scheme
    - 0: the wall temperature is imposed as constant by the user (default)
         and passed to the copain correlation to evaluate the
         exchange coefficient */
  cs_lnum_t *itagms;

} cs_wall_condensation_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Pointer to wall condensation descriptor structure */
extern const cs_wall_condensation_t  *cs_glob_wall_condensation;

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*
 * \brief Provide writable access to _wall_cond structure.
 *
 * \return pointer to global wall_cond structure
 */
/*----------------------------------------------------------------------------*/

cs_wall_condensation_t *
cs_get_glob_wall_condensation(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Set the wall condensation model
 *
 * \param[in] model    integer corresponding to the desired model
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_set_model(cs_wall_cond_natural_conv_model_t  model);

/*----------------------------------------------------------------------------*/
/*
 * \brief Set the onoff state of wall condensation modeling
 *
 * \param[in] icondb  integer corresponding to the onoff state (-1: off, 0: on)
 * \param[in] icondv  integer corresponding to the onoff state with
 *                    metal structures (-1: off, 0: on)
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_set_onoff_state(int  icondb,
                                     int  icondv);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Create the context for wall condensation models.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_create(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Free all structures related to wall condensation models
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_free(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief  Initialize wall condensation models.
 *
 * This includes building the associated meshes.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_initialize(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Reset values
 *
 * \param[in]  wall_cond  wall_condensation strucutre
 * \param[in]  n_var number of variable
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_reset(cs_wall_condensation_t *wall_cond, const int n_var);

/*----------------------------------------------------------------------------*/
/*
 * \brief Return condensing volume structures surface at each cell.
 *
 * \param[out]  surf  array of volume structure surfaces at each cell
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_volume_exchange_surf_at_cells(cs_real_t  *surf);

/*----------------------------------------------------------------------------*/
/*
 * \brief Compute the wall condensation source terms.
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_compute(cs_real_t  total_htc[]);

/*----------------------------------------------------------------------------*/
/*
 * \brief Output statistics about wall condensation source terms (for user log)
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_log(void);

/*----------------------------------------------------------------------------*/
/*
 * \brief Explicit and implicit sources terms from sources
 *        condensation computation.
 *
 * \param[in]      f         pointer to field structure
 * \param[in]      xcpp      array of specific heat (Cp)
 * \param[in]      pvara     variable value at time step beginning
 * \param[in,out]  st_exp    explicit source term part linear in the variable
 * \param[in,out]  st_imp    associated value with \c tsexp
 *                           to be stored in the matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_wall_condensation_source_terms(const cs_field_t  *f,
                                  const cs_real_t    xcpp[],
                                  const cs_real_t    pvara[],
                                  cs_real_t          st_exp[],
                                  cs_real_t          st_imp[]);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_WALL_CONDENSATION_H__ */
