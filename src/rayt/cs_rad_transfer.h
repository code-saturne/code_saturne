#ifndef __CS_RAD_TRANSFER_H__
#define __CS_RAD_TRANSFER_H__

/*============================================================================
 * Radiation solver operations.
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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_time_control.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definition
 *============================================================================*/

/* Radiative transfer model type */

typedef enum {

  CS_RAD_TRANSFER_NONE = 0, /* Set a value to avoid compilation warnings */
  CS_RAD_TRANSFER_DOM  = 1,
  CS_RAD_TRANSFER_P1   = 2

} cs_rad_transfer_model_t;

/* Quadrature types */

typedef enum {

  CS_RAD_QUADRATURE_S4 = 1,
  CS_RAD_QUADRATURE_S6,
  CS_RAD_QUADRATURE_S8,
  CS_RAD_QUADRATURE_T2,
  CS_RAD_QUADRATURE_T4,
  CS_RAD_QUADRATURE_TN,
  CS_RAD_QUADRATURE_LC11,
  CS_RAD_QUADRATURE_DCT020_2468

} cs_rad_quadrature_type_t;

/* Radiative transfer model boundary condition type.
   We use a naming similar to that of cs_boundary.h here, to ease
   later migration. */

enum {

  /*! Grey or black wall with temperature based on fluid BCs */
  CS_BOUNDARY_RAD_WALL_GRAY = 1,

  /*! Grey or black wall with 1D wall model */
  CS_BOUNDARY_RAD_WALL_GRAY_1D_T = 4,

  /*! Grey or black wall with imposed exterior temperature;
    interior wall temperature computed through a flux balance */
  CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T = 21,

  /*! Reflecting wall with imposed exterior temperature;
    interior wall temperature computed through a flux balance
    (same as CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T with
    zero emissivity) */
  CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T = 22,

  /*! Grey or black wall with imposed conduction flux;
    interior wall temperature computed through a flux balance */
  CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX = 31,

  /*! Reflecting wall face to which a conduction flux is imposed,
    which is equivalent to impose this flux directly to the fluid. */
  CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX = 32,

};

/* Radiative transfer model type for atmospheric module */

enum {

  CS_RAD_ATMO_3D_NONE = 0,
  CS_RAD_ATMO_3D_DIRECT_SOLAR = 1 << 0,/* Solar IR band (SIR) absobed by H2O */
  CS_RAD_ATMO_3D_DIRECT_SOLAR_O3BAND = 1 << 1, /* UV-visible band (SUV)
                                                  absobed by H2O */
  CS_RAD_ATMO_3D_DIFFUSE_SOLAR = 1 << 2,/* Solar IR band (SIR) absobed by H2O */
  CS_RAD_ATMO_3D_DIFFUSE_SOLAR_O3BAND = 1 << 3, /* UV-visible band (SUV)
                                                   absobed by H2O */
  CS_RAD_ATMO_3D_INFRARED = 1 << 4

};

/*============================================================================
 *  Global variables
 *============================================================================*/

/*! Model name */
extern const char *cs_rad_transfer_model_name[];

/*! Quadrature name */
extern const char *cs_rad_transfer_quadrature_name[];

typedef struct {

  cs_rad_transfer_model_t   type;  /*!< model activation and type */

  int           nrphas;
  int           iimpar;
  int           verbosity; /*!< Radiance resolution verbosity */
  int           imodak; /*!< Absorption coefficent computation
                             0: without modak
                             1: with modak */
  int           imoadf; /*!< ADF model:
                             0: not used
                             1: with wavelength interval of a 8
                             2: with wavelength interval of a 50 */
  int           iwrp1t;
  int           imfsck; /*!< FSCK model (0: off, 1: on) */
  double        xnp1mx;
  int           idiver; /*!< Explicit radiative source term computation mode
                             -1: no renormalization
                             0: Semi-analytic (mandatory if transparent)
                             1: Conservative
                             2: Corrected semi-analytic (to be conservative)
                             REMARK: if transparent, idiver = -1
                             automatically in DOM */
  int           i_quadrature;
  int           ndirec;
  int           ndirs;
  cs_real_3_t  *vect_s;
  cs_real_t    *angsol;
  int           restart;
  int           nwsgg;
  cs_real_t    *wq;
  int           itpimp;
  int           ipgrno;
  int           iprefl;
  int           ifgrno;
  int           ifrefl;
  int           itpt1d;
  int           ifinfe;

  int           atmo_model;          /*!< Atmospheric radiation model:
                                          - Direct Solar the first bit
                                          - diFfuse Solar for the second bit
                                          - InfraRed for the third bit */
  int           atmo_dr_id;          /*!< Atmospheric radiation model:
                                          id of the Direct Solar band
                                          or -1 if not activated
                                          (automatically computed) */
  int           atmo_dr_o3_id;       /*!< Atmospheric radiation model:
                                          id of the Direct Solar O3 band
                                          or -1 if not activated
                                          (automatically computed) */
  int           atmo_df_id;          /*!< Atmospheric radiation model:
                                          id of the Diffuse Solar band
                                          or -1 if not activated
                                          (automatically computed) */
  int           atmo_df_o3_id;       /*!< Atmospheric radiation model:
                                          id of the Diffuse Solar O3 band (SUV)
                                          or -1 if not activated
                                          (automatically computed) */
  int           atmo_ir_id;          /*!< Atmospheric radiation model:
                                          id of the InfraRed band
                                          or -1 if not activated
                                          (automatically computed) */
  bool          dispersion;          /*!< add dispersion (through diffusion) */
  cs_real_t     dispersion_coeff;    /*!< dispersion coefficient.
                                       The dispersion coefficient leading to the
                                       best precision may depend on the chosen
                                       quadrature, and has been observed to be
                                       3 for 128 directions (T4) and
                                       5 for 32  directions (T2) on a (cube
                                       with point source) test case; the default
                                       value of 1 already improves precision in
                                       both cases. */

  cs_time_control_t  time_control;   /* Time control for radiation updates */


} cs_rad_transfer_params_t;

extern cs_rad_transfer_params_t *cs_glob_rad_transfer_params;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Finalize radiative transfer module.
 */
/*----------------------------------------------------------------------------*/

void
cs_rad_transfer_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Is time step for radiative transfer active?
 *
 * \return  true if active, false otherwise.
 */
/*----------------------------------------------------------------------------*/

bool
cs_rad_time_is_active(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_RAD_TRANSFER_H__ */
